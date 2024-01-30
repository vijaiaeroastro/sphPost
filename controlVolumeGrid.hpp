#ifndef VIJAI_SPH_POST_CONTROL_VOLUME_HPP
#define VIJAI_SPH_POST_CONTROL_VOLUME_HPP
#pragma once

#include "sphFrame.hpp"

namespace vSPHPost {
    // An euclidean distance based adaptor for flann
    class flannAdaptor {
    public:
        const Eigen::MatrixXd &data;

        flannAdaptor(const Eigen::MatrixXd &data_) : data(data_) {}

        inline std::size_t kdtree_get_point_count() const { return data.rows(); }

        inline double kdtree_get_pt(const std::size_t idx, int dim) const {
            return data(idx, dim);
        }

        template<class BBOX>
        bool kdtree_get_bbox(BBOX &) const { return false; }
    };

    class controlVolumeGrid {
    public:
        controlVolumeGrid(std::vector<std::array<double, 3>> &_control_volume_2d,
                          sphFrame *_frame_ptr,
                          double _smoothing_length) :
                control_volume_2d(_control_volume_2d),
                frame(_frame_ptr),
                smoothing_length(_smoothing_length) {
            initializeStuff();
        }

        ~controlVolumeGrid() = default;

    public:
        double estimateTotalMassChange();

        Eigen::Vector3d estimateTotalMomentumChange();

    private:
        void estimateMassFlux();

        void estimateMomentumFlux();

        void estimateShearRate();

    private:
        void initializeStuff();

    private:
        void buildCartesianMesh();

        void neighbourSearch();

        void interpolateAllScalarFields();

        void interpolateAllVectorFields();

        void estimateMissingFields();

    private:
        double euclideanDistance(const Eigen::Vector3d &p1, const Eigen::Vector3d &p2);

    private:
        // Methods here are implemented for a regular grid explicitly

        Eigen::VectorXd interpolateScalarQuantity(const std::string field);

        std::vector<Eigen::Vector3d> interpolateVectorQuantity(const std::string field);

    private:
        // Diffusion methods listed here will actually work on cartesian mesh / triangle / quad / poly mesh since they operate on
        // raw point cloud and always constructs an intermediate representation itself

        // Direct diffusion
        Eigen::VectorXd
        extendScalarFieldUsingDiffusionAlgorithm(
                const std::string field);

        // Direct diffusion
        std::vector<Eigen::Vector3d>
        extendVectorFieldUsingDiffusionAlgorithm(
                const std::string field);

    private:
        double cubicSplineKernel(unsigned int dimension, double r, double h);

    public:
        std::vector<Eigen::Vector3d> meshPoints;
        std::array<unsigned int, 3> meshDimensions;
        std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, double>>> meshPointNeighbours;
        std::vector<Eigen::Vector3d> velocity;
        std::vector<Eigen::Vector3d> normals;
        std::vector<Eigen::Vector3d> acceleration;
        std::vector<Eigen::Vector3d> mass_flux;
        std::vector<Eigen::Vector3d> shearRate;
        std::vector<Eigen::Matrix3d> momentumFlux;
        Eigen::VectorXd areas;
        Eigen::VectorXd mass;
        Eigen::VectorXd pressure;
        Eigen::VectorXd density;

    private:
        std::vector<std::array<double, 3>> control_volume_2d;
        sphFrame *frame;
        unsigned int id;
        double smoothing_length;
    };
}

#endif