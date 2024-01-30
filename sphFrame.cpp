#include "sphFrame.hpp"
#include "vtpReaderWrapper.hpp"
#include <pugixml.hpp>
#include <iostream>
#include <filesystem>
#include <stdexcept>

namespace vSPHPost {

    sphFrame::sphFrame(unsigned int _frame_id,
                       std::string _vtp_file,
                       std::unordered_map<unsigned int, std::vector<std::array<double, 3>>> *_control_volumes) {
        frame_id = _frame_id;
        vtp_file = std::move(_vtp_file);
        control_volumes = _control_volumes;
        populateData();
    }

    std::vector<std::array<double, 3>> sphFrame::getVectorField(const std::string &_field_name) {
        return vector_fields[_field_name];
    }

    std::vector<double> sphFrame::getScalarField(const std::string &_field_name) {
        return scalar_fields[_field_name];
    }

    std::vector<std::array<double, 3>> sphFrame::getBoundingBox() {
        auto tempPositions = vector_fields["positions"];
        double minX, minY, maxX, maxY;
        std::vector<double> allX, allY;
        for (auto iter: tempPositions) {
            allX.push_back(iter[0]);
            allY.push_back(iter[1]);
        }
        std::sort(allX.begin(), allX.end());
        std::sort(allY.begin(), allY.end());
        minX = allX[0];
        maxX = allX[allX.size() - 1];
        minY = allY[0];
        maxY = allY[allY.size() - 1];
        std::vector<std::array<double, 3>> bounding_box;
        bounding_box.push_back({minX, minY, 0.0});
        bounding_box.push_back({maxX, minY, 0.0});
        bounding_box.push_back({maxX, maxY, 0.0});
        bounding_box.push_back({minX, maxY, 0.0});
        std::cout << minX << "," << minY << " --> " << maxX << "," << maxY << std::endl;
        return bounding_box;
    }

    unsigned int sphFrame::getFrameId() const {
        return frame_id;
    }

    std::unordered_map<unsigned int, std::vector<std::array<double, 3>>> *sphFrame::getCV() {
        return control_volumes;
    }

    void sphFrame::populateData() {
        pugi::xml_document vtkDocument;
        if (vtkDocument.load_file(vtp_file.c_str())) {
            std::string rawVTKPath = vtkDocument
                    .child("VTKFile")
                    .child("PPolyData")
                    .child("Piece").attribute("Source")
                    .as_string();
            std::filesystem::path original_path(vtp_file);
            vtp_file = original_path.parent_path().generic_string() + "/" + rawVTKPath;
            extractVTPData();
        } else {
            throw std::invalid_argument("Error: Unable to Parse VTP file");
        }
    }

    void sphFrame::extractVTPData() {
        auto vtpDataPair = readVTPDataset(vtp_file);
        scalar_fields = vtpDataPair.first;
        vector_fields = vtpDataPair.second;

        for (auto iter: scalar_fields) {
            auto arrayName = iter.first;
            auto arrayEntries = iter.second;
            if (arrayName == "mass") {
                mass.resize(arrayEntries.size());
                for (unsigned int i = 0; i < arrayEntries.size(); ++i) {
                    mass[i] = arrayEntries[i];
                }
            }
            if (arrayName == "pressure") {
                pressure.resize(arrayEntries.size());
                for (unsigned int i = 0; i < arrayEntries.size(); ++i) {
                    pressure[i] = arrayEntries[i];
                }
            }
            if (arrayName == "density") {
                density.resize(arrayEntries.size());
                for (unsigned int i = 0; i < arrayEntries.size(); ++i) {
                    density[i] = arrayEntries[i];
                }
            }
            if (arrayName == "particle_id") {
                particle_ids.resize(arrayEntries.size());
                for (unsigned int i = 0; i < arrayEntries.size(); ++i) {
                    particle_ids[i] = arrayEntries[i];
                }
            }
        }

        for (auto iter: vector_fields) {
            auto arrayName = iter.first;
            auto arrayEntries = iter.second;
            if (arrayName == "positions") {
                particles.resize(arrayEntries.size(), 3);
                for (unsigned int i = 0; i < arrayEntries.size(); ++i) {
                    auto currentEntry = arrayEntries.at(i);
                    for (unsigned int dim = 0; dim < currentEntry.size(); ++dim) {
                        particles(i, dim) = currentEntry[dim];
                    }
                }
//                std::cout << "PARTICLE SIZE : " << particles.rows() << std::endl;
            }
            if (arrayName == "velocity") {
                velocity.resize(arrayEntries.size(), 3);
                for (unsigned int i = 0; i < arrayEntries.size(); ++i) {
                    auto currentEntry = arrayEntries.at(i);
                    for (unsigned int dim = 0; dim < currentEntry.size(); ++dim) {
                        velocity(i, dim) = currentEntry[dim];
                    }
                }
            }
            if (arrayName == "acceleration") {
                acceleration.resize(arrayEntries.size(), 3);
                for (unsigned int i = 0; i < arrayEntries.size(); ++i) {
                    auto currentEntry = arrayEntries.at(i);
                    for (unsigned int dim = 0; dim < currentEntry.size(); ++dim) {
                        acceleration(i, dim) = currentEntry[dim];
                    }
                }
            }
            if (arrayName == "normal") {
                normals.resize(arrayEntries.size(), 3);
                for (unsigned int i = 0; i < arrayEntries.size(); ++i) {
                    auto currentEntry = arrayEntries.at(i);
                    for (unsigned int dim = 0; dim < currentEntry.size(); ++dim) {
                        normals(i, dim) = currentEntry[dim];
                    }
                }
            }
        }
    }
}