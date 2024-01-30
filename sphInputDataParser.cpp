#include "sphInputDataParser.hpp"
#include <algorithm>
#include <fstream>

#include "json.hpp"

namespace vSPHPost {

    sphInputDataParser::sphInputDataParser(const std::string &_configuration_file) : configuration_file(_configuration_file) {
        std::ifstream configFile(configuration_file);
        nlohmann::json configData;
        if (configFile.is_open()) {
            try {
                configFile >> configData;
                // Parse simulation results directory
                directoryPath = configData["SIMULATION_RESULTS"];
                // Parse control volumes for post-processing
                for (auto iter: configData["CONTROL_VOLUMES"]) {
                    // Order list of points forming a polygon (here just a box)
                    auto vertexId = iter["ID"].get<unsigned int>();
                    auto vertexData = iter["VERTICES"].get<std::vector<std::array<double, 3>>>();
                    control_volumes[vertexId] = vertexData;
                }
                smoothing_length = configData["SMOOTHING_LENGTH"];
            } catch (const std::exception &e) {
                throw std::runtime_error("Error parsing JSON");
            }
            configFile.close();
        } else {
            throw std::runtime_error("Error parsing JSON");
        }
    }

    std::vector<sphFrame> sphInputDataParser::getSPHDataFrame() {
        auto vtp_files = processPVTPFiles();
        std::vector<sphFrame> sph_data;
        for (unsigned int i = 0; i < vtp_files.size(); ++i) {
            sphFrame current_frame(i, vtp_files.at(i), &control_volumes);
            sph_data.push_back(current_frame);
        }
        return sph_data;
    }

    std::unordered_map<unsigned int, std::vector<std::array<double, 3>>> sphInputDataParser::getControlVolumes() {
        return control_volumes;
    }

    double sphInputDataParser::getSmoothingLength() {
        return smoothing_length;
    }

    int sphInputDataParser::extractNumber(const std::filesystem::path &path) {
        std::string filename = path.stem().string(); // Get the filename without the extension
        // Remove the non-digit characters and convert the remainder to an integer
        filename.erase(std::remove_if(filename.begin(), filename.end(),
                                      [](char c) { return !std::isdigit(c); }),
                       filename.end());
        return std::stoi(filename);
    }

    std::vector<std::filesystem::path> sphInputDataParser::sortPVTPFiles(const std::string &directoryPath) {
        std::vector<std::filesystem::path> paths;

        // Check if provided path is a directory
        if (!std::filesystem::is_directory(directoryPath)) {
            std::cerr << "Error: The provided path is not a directory." << std::endl;
            return paths;
        }

        // Populate the vector with paths to .pvtp files
        for (const auto &entry: std::filesystem::directory_iterator(directoryPath)) {
            if (entry.is_regular_file() && entry.path().extension() == ".pvtp") {
                paths.push_back(entry.path());
            }
        }

        // Sort the vector based on the numerical value in the file names
        std::sort(paths.begin(), paths.end(), [](const std::filesystem::path &a, const std::filesystem::path &b) {
            return extractNumber(a) < extractNumber(b);
        });

        return paths;
    }

    std::vector<std::string> sphInputDataParser::processPVTPFiles() {
        // Check if provided path is a directory
        if (!std::filesystem::is_directory(directoryPath)) {
            throw std::invalid_argument("Error: The provided path is not a directory.");
        } else {
            std::vector<std::string> absolute_path_vtp_files;
            for (const auto &iter: sortPVTPFiles(directoryPath)) {
                absolute_path_vtp_files.push_back(std::filesystem::absolute(iter).generic_string());
            }
            return absolute_path_vtp_files;
        }
    }
}
