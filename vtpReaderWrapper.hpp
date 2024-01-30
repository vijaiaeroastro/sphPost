#ifndef VTPREADERWRAPPER_LIBRARY_H
#define VTPREADERWRAPPER_LIBRARY_H

#include <unordered_map>
#include <vector>
#include <array>
#include <string>
#include <utility>

typedef std::unordered_map<std::string, std::vector<std::array<double,3>>> vFieldMapType;
typedef std::unordered_map<std::string, std::vector<double>> sFieldMapType;

std::pair<sFieldMapType, vFieldMapType> readVTPDataset(const std::string _vtp_file);


#endif //VTPREADERWRAPPER_LIBRARY_H
