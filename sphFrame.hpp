#ifndef VIJAI_SPH_POST_SPH_FRAME_HPP
#define VIJAI_SPH_POST_SPH_FRAME_HPP
#pragma once

#include <utility>
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>

namespace vSPHPost {
    class sphFrame {
    public:
        sphFrame(unsigned int _frame_id,
                 std::string _vtp_file,
                 std::unordered_map<unsigned int, std::vector<std::array<double, 3>>> *_control_volumes);

        ~sphFrame() = default;

    public:
        std::vector<std::array<double, 3>> getVectorField(const std::string &_field_name);

        std::vector<double> getScalarField(const std::string &_field_name);

        std::vector<std::array<double, 3>> getBoundingBox();

        unsigned int getFrameId() const;

        std::unordered_map<unsigned int, std::vector<std::array<double, 3>>> *getCV();

    private:
        void populateData();

        void extractVTPData();

    private:
        unsigned int frame_id;
        std::string vtp_file;

    public:
        Eigen::MatrixXd particles;
        Eigen::MatrixXd velocity;
        Eigen::MatrixXd normals;
        Eigen::MatrixXd acceleration;
        Eigen::VectorXd mass;
        Eigen::VectorXd pressure;
        Eigen::VectorXd density;
        Eigen::VectorXi particle_ids;
        std::unordered_map<std::string, std::vector<std::array<double, 3>>> vector_fields;
        std::unordered_map<std::string, std::vector<double>> scalar_fields;
        std::unordered_map<unsigned int, std::vector<std::array<double, 3>>> *control_volumes;
    };
}

#endif