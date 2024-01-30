#include "vtpReaderWrapper.hpp"
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>

std::pair<sFieldMapType, vFieldMapType> readVTPDataset(const std::string vtp_file) {
    // Declare the fields
    sFieldMapType scalar_fields;
    vFieldMapType vector_fields;

    // Create a reader for the VTP file
    vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName(vtp_file.c_str());
    reader->Update(); // Read the file

    // The output of the reader is the polydata object
    vtkSmartPointer<vtkPolyData> polyData = reader->GetOutput();

    // Access point data
    auto pointData = polyData->GetPointData();
    if (pointData) {
        for (unsigned int arrayIdx = 0; arrayIdx < pointData->GetNumberOfArrays(); arrayIdx++) {
            vtkDataArray *array = pointData->GetArray(arrayIdx);
            if (array) {
                std::string arrayName = array->GetName();
                int numComponents = array->GetNumberOfComponents();

                // Inspect each array based on its name and properties
//                std::cout << "Array Name: " << arrayName << std::endl;
//                std::cout << "Number of Components: " << numComponents << std::endl;

                // Add more inspections or processing as needed
                // mass, pressure, density, particle_id (scalar fields)
                if ((arrayName == "mass") || (arrayName == "pressure") || (arrayName == "density") ||
                    (arrayName == "particle_id")) {
                    std::vector<double> current_scalar_field;
                    for (vtkIdType i = 0; i < array->GetNumberOfTuples(); ++i) {
                        double value;
                        array->GetTuple(i, &value);
                        current_scalar_field.push_back(value);
                    }
                    scalar_fields[arrayName] = current_scalar_field;
                }

                // velocity, acceleration, normal
                if ((arrayName == "velocity") || (arrayName == "acceleration") || (arrayName == "normal")) {
                    std::vector<std::array<double,3>> current_vector_field;
                    for (vtkIdType i = 0; i < array->GetNumberOfTuples(); ++i) {
                        double values[3];
                        array->GetTuple(i, values);
                        current_vector_field.push_back(std::array<double,3>({values[0], values[1], values[2]}));
                    }
                    vector_fields[arrayName] = current_vector_field;
                }
            } else {
                std::cerr << "This array seems invalid : "  << arrayIdx << std::endl;
            }
        }
    }

    // Retrieve the points from the polyData
    vtkSmartPointer<vtkPoints> points = polyData->GetPoints();
    if (points) {
        std::vector<std::array<double,3>> position_data;
        for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
            double p[3];
            points->GetPoint(i, p);
            position_data.push_back(std::array<double,3>({p[0], p[1], p[2]}));
        }
        vector_fields["positions"] = position_data;
    } else {
        std::cerr << "Couldn't find position data" << std::endl;
    }

    return std::make_pair(scalar_fields, vector_fields);
}