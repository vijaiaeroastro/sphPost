#pragma once

#include <string>
#include <vector>
#include <array>
#include <Eigen/Dense>

class vRender {
public:
    vRender();

    ~vRender();

public:
    void show();

public:
    std::string addPointCloud(const std::string _name, std::vector<std::array<double,3>> &_point_cloud);

    std::string addPointCloud(const std::string _name, std::vector<Eigen::Vector3d> &_point_cloud);

    std::string addScalarQuantityToPointCloud(const std::string _cloud_name, const std::string _scalar_field_name,
                                              std::vector<double> &_scalar_field);

    std::string addVectorQuantityToPointCloud(const std::string _cloud_name, const std::string _vector_field_name,
                                              std::vector<std::array<double,3>> &_vector_field);

    std::string addVectorQuantityToPointCloud(const std::string _cloud_name, const std::string _vector_field_name,
                                              std::vector<Eigen::Vector3d> &_vector_field);

    std::string drawCloseLoopCurve(const std::string _name, std::vector<std::array<double,3>> &_curve_loop);

    void addResidualPlot(const std::string _name, const std::string _x_name, const std::string _y_name,
                    std::vector<double> &quantity
    );
};