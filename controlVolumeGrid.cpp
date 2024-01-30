#include <iostream>
#include "controlVolumeGrid.hpp"
#include "nanoflann.hpp"

namespace vSPHPost {
    void controlVolumeGrid::initializeStuff() {
        // build a cartesian mesh with bounds and use a resolution lower than inter particle distance
        // find the particles contained in the control volume
        // Write 3 functions called interpolate scalar field, vector field and tensor field
        // Interpolate all existing quantities from simulation onto the nodes of cartesian mesh
        buildCartesianMesh();
        // build a kdtree on input points (with proper ids)
        // for every point in cartesian mesh, look for nearest neighbours within smoothing length radius
        neighbourSearch();
        // average their values (component wise and assign them)
        interpolateAllScalarFields();
        interpolateAllVectorFields();
        estimateMissingFields();
    }

    void controlVolumeGrid::buildCartesianMesh() {
        // Calculate the number of points in each direction
        auto minCorner = control_volume_2d[0];
        auto maxCorner = control_volume_2d[2];
        double delta = smoothing_length / 5.0;
        size_t numPointsX = static_cast<size_t>((maxCorner[0] - minCorner[0]) / delta) + 1;
        size_t numPointsY = static_cast<size_t>((maxCorner[1] - minCorner[1]) / delta) + 1;
        meshDimensions[0] = numPointsX;
        meshDimensions[1] = numPointsY;

        // Create the mesh points
        for (size_t i = 0; i < numPointsX; ++i) {
            for (size_t j = 0; j < numPointsY; ++j) {
                meshPoints.emplace_back(Eigen::Vector3d(minCorner[0] + i * delta, minCorner[1] + j * delta, 0.0));
            }
        }

        // Estimate areas too (it can be used later) -- The implementation below is a bit wonky. The right approach would be
        // to calculate a dual Cartesian mesh and calculate areas of quads. But the approximation below is a decent starting
        // point
        areas.resize(numPointsX * numPointsY);
        for (size_t i = 0; i < numPointsX; ++i) {
            for (size_t j = 0; j < numPointsY; ++j) {
                double area = delta * delta; // Default area for interior points

                bool isOnVerticalBoundary = (i == 0 || i == numPointsX - 1);
                bool isOnHorizontalBoundary = (j == 0 || j == numPointsY - 1);

                if (isOnVerticalBoundary && isOnHorizontalBoundary) {
                    // Corner point
                    area *= 0.25;
                } else if (isOnVerticalBoundary || isOnHorizontalBoundary) {
                    // Edge point, but not a corner
                    area *= 0.5;
                }

                areas(i * numPointsY + j) = area;
            }
        }
    }

    double controlVolumeGrid::euclideanDistance(const Eigen::Vector3d &p1, const Eigen::Vector3d &p2) {
        return (p1 - p2).norm();
    }

    void controlVolumeGrid::neighbourSearch() {
        // Build a kdtree and use smoothing length as the search radius and
        // find all neighbours for every point within smoothing length
        // For some reason, radius search is weird, so i do index search and filter manually (still fast :D)
        flannAdaptor adaptor(frame->particles);
        typedef nanoflann::KDTreeSingleIndexAdaptor<
                nanoflann::L2_Simple_Adaptor<double, flannAdaptor>,
                flannAdaptor, 2> KDTree;
        KDTree kdTree(2, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(10));
        kdTree.buildIndex();
        unsigned int num_neighbors = 25;
        for (unsigned int i = 0; i < meshPoints.size(); ++i) {
            auto query_point = meshPoints[i];
            nanoflann::KNNResultSet<double> resultSet(num_neighbors);
            std::vector<size_t> indices(num_neighbors);
            std::vector<double> distances(num_neighbors);
            resultSet.init(&indices[0], &distances[0]);
            kdTree.findNeighbors(resultSet, &query_point[0], nanoflann::SearchParameters(10));
            std::vector<std::pair<unsigned int, double>> current_matches;
            for (unsigned int neigh = 0; neigh < indices.size(); ++neigh) {
                auto current_index = indices[neigh];
                auto nX = frame->particles(current_index, 0);
                auto nY = frame->particles(current_index, 1);
                auto nZ = frame->particles(current_index, 2);
                Eigen::Vector3d neigh_pt(nX, nY, nZ);
                double distance = euclideanDistance(query_point, neigh_pt);
                if (distance <= smoothing_length) {
                    current_matches.push_back(std::make_pair(current_index, distance));
                }
            }
            if (current_matches.size() > 0) {
                meshPointNeighbours[i] = current_matches;
            }
        }
    }

    double controlVolumeGrid::cubicSplineKernel(unsigned int dimension, double r, double h) {
        double sigma;
        // Dimension is capped at 3
        if (dimension > 3) {
            dimension = 3;
        }
        if (dimension == 1) {
            sigma = 2.0 / (3.0 * M_PI * h);
        } else if (dimension == 2) {
            sigma = 10.0 / (7.0 * M_PI * h * h);
        } else {
            sigma = 1.0 / (M_PI * h * h * h);
        }
        // Dimensionless distance
        double q = r / h;
        // Cubic spline kernel calculation
        if (q >= 0 && q <= 1) {
            return sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q);
        } else if (q > 1 && q <= 2) {
            double term = 2.0 - q;
            return sigma * 0.25 * term * term * term;
        }
        return 0.0;
    }

    Eigen::VectorXd controlVolumeGrid::interpolateScalarQuantity(const std::string _scalar_field) {
        auto scalar_field_vector = frame->scalar_fields[_scalar_field];
        Eigen::VectorXd converted_one(scalar_field_vector.size());
        for (unsigned int i = 0; i < scalar_field_vector.size(); ++i) {
            converted_one(i) = scalar_field_vector[i];
        }
        Eigen::VectorXd interpolated_thing(meshPoints.size());
        for (size_t i = 0; i < meshPoints.size(); ++i) {
            double weightedSum = 0.0;
            double weightSum = 0.0;
            for (auto iter: meshPointNeighbours[i]) {
                double distance = iter.second;
                double weight = cubicSplineKernel(2, distance, smoothing_length);
                weightedSum += converted_one(iter.first) * weight;
                weightSum += weight;
            }
            interpolated_thing(i) = weightSum > 0 ? weightedSum / weightSum : 0;
        }
        return interpolated_thing;
    }

    std::vector<Eigen::Vector3d> controlVolumeGrid::interpolateVectorQuantity(const std::string field) {
        auto vector_field_vector = frame->vector_fields[field];
        std::vector<Eigen::Vector3d> converted_one(vector_field_vector.size());
        for (unsigned int i = 0; i < vector_field_vector.size(); ++i) {
            auto temp_vec = vector_field_vector.at(i);
            converted_one[i] = Eigen::Vector3d(temp_vec[0], temp_vec[1], temp_vec[2]);
        }
        std::vector<Eigen::Vector3d> interpolated_thing(meshPoints.size());
        for (size_t i = 0; i < meshPoints.size(); ++i) {
            Eigen::Vector3d weightedSum = Eigen::Vector3d::Zero();
            double weightSum = 0.0;
            for (auto iter: meshPointNeighbours[i]) {
                double distance = iter.second;
                double weight = cubicSplineKernel(2, distance, smoothing_length);
                weightedSum += converted_one[iter.first] * weight;
                weightSum += weight;
            }
            if (weightSum > 0) {
                auto interpolated_value = weightedSum / weightSum;
                interpolated_thing[i] = interpolated_value;
//                std::cout << interpolated_value << std::endl;
            } else {
                interpolated_thing[i] = Eigen::Vector3d::Zero();
            }
        }
        return interpolated_thing;
    }


    void controlVolumeGrid::interpolateAllScalarFields() {
        // mass, density, pressure
        density = interpolateScalarQuantity("density");
        mass = interpolateScalarQuantity("mass");
        pressure = interpolateScalarQuantity("pressure");
    }

    void controlVolumeGrid::interpolateAllVectorFields() {
        // velocity, normals, acceleration
        velocity = interpolateVectorQuantity("velocity");
        normals = interpolateVectorQuantity("normal");
        acceleration = interpolateVectorQuantity("acceleration");
    }

    void controlVolumeGrid::estimateMissingFields() {
        estimateMassFlux();
        estimateShearRate();
        estimateMomentumFlux();
    }

    void controlVolumeGrid::estimateMassFlux() {
        mass_flux.resize(meshPoints.size());
        for (unsigned int i = 0; i < meshPoints.size(); ++i) {
            mass_flux[i] = density[i] * velocity[i];
        }
    }

    void controlVolumeGrid::estimateShearRate() {
        auto numX = meshDimensions[0];
        auto numY = meshDimensions[1];
        shearRate.resize(numX * numY, Eigen::Vector3d::Zero());
        auto getVelocity = [=](int i, int j) -> Eigen::Vector3d {
            if (i >= 0 && i < numX && j >= 0 && j < numX * numY) {
                return velocity[j * numX + i];
            }
            return Eigen::Vector3d::Zero();  // Handle out-of-bounds by returning zero velocity
        };
        auto dx = smoothing_length / 5.0;
        auto dy = smoothing_length / 5.0;

        for (int j = 0; j < numY; ++j) {
            for (int i = 0; i < numX; ++i) {
                double dvdx, dudy;

                // Calculate dvdx
                if (i == 0) {
                    dvdx = (getVelocity(i + 1, j)[1] - getVelocity(i, j)[1]) / dx;
                } else if (i == numX - 1) {
                    dvdx = (getVelocity(i, j)[1] - getVelocity(i - 1, j)[1]) / dx;
                } else {
                    dvdx = (getVelocity(i + 1, j)[1] - getVelocity(i - 1, j)[1]) / (2 * dx);
                }

                // Calculate dudy
                if (j == 0) {
                    dudy = (getVelocity(i, j + 1)[0] - getVelocity(i, j)[0]) / dy;
                } else if (j == numY - 1) {
                    dudy = (getVelocity(i, j)[0] - getVelocity(i, j - 1)[0]) / dy;
                } else {
                    dudy = (getVelocity(i, j + 1)[0] - getVelocity(i, j - 1)[0]) / (2 * dy);
                }

                // The shear rate is the sum of dvdx and dudy
                // We use the x component of the Vector3d to store the shear rate
                shearRate[j * numX + i][0] = dvdx + dudy;
                // y and z components are kept zero
            }
        }
    }

    void controlVolumeGrid::estimateMomentumFlux() {
        momentumFlux.resize(meshPoints.size());

        for (int i = 0; i < meshPoints.size(); ++i) {
            // Compute convective momentum flux for the i-th point
            Eigen::Matrix3d convectiveFlux = density[i] * velocity[i] * velocity[i].transpose();

            // Compute pressure term for the i-th point
            Eigen::Matrix3d pressureTerm = pressure[i] * Eigen::Matrix3d::Identity();

            // Compute viscous term for the i-th point
            Eigen::Matrix3d viscousTerm;
            viscousTerm << 0, shearRate[i].x(), 0,
                    shearRate[i].y(), 0, 0,
                    0, 0, 0; // Assuming 2D problem, so no shear in the z-direction

            viscousTerm *= 0.001;

            // Combine the terms to get the momentum flux for the i-th point
            Eigen::Matrix3d totalMomentumFlux = convectiveFlux + pressureTerm - viscousTerm;

            momentumFlux[i] = totalMomentumFlux;
        }
    }

    double controlVolumeGrid::estimateTotalMassChange() {
        double total_mass_change = 0.0;
        for(unsigned int i = 0; i < meshPoints.size(); ++i) {
            // ignore stuff with negative density
            if(density(i) > 0) {
                total_mass_change += (mass_flux[i].dot(areas[i] * normals[i]));
            }
        }
        return total_mass_change;
    }

    Eigen::Vector3d controlVolumeGrid::estimateTotalMomentumChange() {
        Eigen::Vector3d total_momentum_change = Eigen::Vector3d::Zero();
        for(unsigned int i = 0; i < meshPoints.size(); ++i) {
            // ignore stuff with negative density
            if(density(i) > 0) {
                total_momentum_change +=  (momentumFlux[i] * (areas[i] * normals[i]));
            }
        }
        return total_momentum_change;
    }
}