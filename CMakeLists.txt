cmake_minimum_required(VERSION 3.27)
project(SPHPost)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED)
add_subdirectory(deps/pugixml-1.14)

# VTP reader (only capabilities for vtp file reading, since python package does not work)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/deps)
set(VTK_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../DEPS/VTK_930_INSTALL/lib/cmake/vtk-9.3)
find_package(VTK COMPONENTS
        CommonCore
        CommonDataModel
        IOXML)

# Post processor headers and sources
set(HDRS
        sphFrame.hpp
        sphInputDataParser.hpp
        controlVolumeGrid.hpp
        vtpReaderWrapper.hpp
)

set(SRCS
        sphFrame.cpp
        sphInputDataParser.cpp
        vtpReaderWrapper.cpp
        controlVolumeGrid.cpp)

# Library first
add_library(SPHPost SHARED
        ${HDRS}
        ${SRCS})

target_link_libraries(SPHPost
        pugixml
        Eigen3::Eigen
        ${VTK_LIBRARIES})

# Executable is next
add_executable(SPHMain main.cpp)

target_link_libraries(SPHMain
        SPHPost
        ${CMAKE_CURRENT_SOURCE_DIR}/deps/vRender/libvRenderLib.so)
