#include "sphInputDataParser.hpp"
#include <vRender/vRender.hpp>
#include "controlVolumeGrid.hpp"

std::vector<vSPHPost::sphFrame> simulation_data;
std::unordered_map<unsigned int, std::vector<std::array<double, 3>>> control_volumes;
std::vector<std::array<double, 3>> bounding_boxes;
std::unordered_map<unsigned int, std::vector<Eigen::Vector3d>> meshPoints;
std::vector<double> total_mass_change;

int main() {
    vSPHPost::sphInputDataParser vtp_parser("../input_data.json");
    vRender renderer;
    simulation_data = vtp_parser.getSPHDataFrame();
    control_volumes = vtp_parser.getControlVolumes();
    auto positions = simulation_data[simulation_data.size() - 1].getVectorField("positions");
    auto velocity = simulation_data[simulation_data.size() - 1].getVectorField("velocity");
    auto density = simulation_data[simulation_data.size() - 1].getScalarField("density");
    auto pressure = simulation_data[simulation_data.size() - 1].getScalarField("pressure");
    auto cloud = renderer.addPointCloud("SPH_POSITIONS", positions);
    renderer.addScalarQuantityToPointCloud(cloud, "density", density);
    renderer.addScalarQuantityToPointCloud(cloud, "pressure", pressure);
    renderer.addVectorQuantityToPointCloud(cloud, "velocity", velocity);
//    bounding_boxes = simulation_data[0].getBoundingBox();
    for(auto cv: control_volumes) {
        std::string name =  "CV_" + std::to_string(cv.first);
        renderer.drawCloseLoopCurve(name, cv.second);
    }

    for (unsigned int i = 0; i < simulation_data.size(); ++i) {
        auto frame = simulation_data[i];
        vSPHPost::controlVolumeGrid cv_grid(control_volumes[0],
                                            &frame,
                                            vtp_parser.getSmoothingLength());
        auto totalMassChange = cv_grid.estimateTotalMassChange();
        total_mass_change.push_back(totalMassChange);
        std::cout << frame.getFrameId() << " --> " << totalMassChange << std::endl;
        if(i == simulation_data.size() - 1) {
            auto cv_temp = renderer.addPointCloud("CV_0_CLOUD", cv_grid.meshPoints);
            renderer.addVectorQuantityToPointCloud(cv_temp, "INTERP_VELOCITY", cv_grid.velocity);
        }
    }
    renderer.show();
    return 0;
}