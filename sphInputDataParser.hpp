#ifndef VIJAI_SPH_POST_INPUT_PARSER_HPP
#define VIJAI_SPH_POST_INPUT_PARSER_HPP
#pragma once


#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include "sphFrame.hpp"

namespace vSPHPost {
    class sphInputDataParser {
    public:
        sphInputDataParser(const std::string &_configuration_file);

        ~sphInputDataParser() = default;

    public:
        std::vector<sphFrame> getSPHDataFrame();

        std::unordered_map<unsigned int, std::vector<std::array<double, 3>>> getControlVolumes();

        double getSmoothingLength();

    private:
        static int extractNumber(const std::filesystem::path &path);

        static std::vector<std::filesystem::path> sortPVTPFiles(const std::string &directoryPath);

        std::vector<std::string> processPVTPFiles();

    private:
        std::string directoryPath;
        std::string configuration_file;
        std::unordered_map<unsigned int, std::vector<std::array<double, 3>>> control_volumes;
        double smoothing_length;
    };
}

#endif