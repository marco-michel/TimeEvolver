#pragma once

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#ifdef USE_HDF
#include <H5Cpp.h>
#endif

namespace test_utils {

class cleanup_guard {
public:
    cleanup_guard() = default;

    explicit cleanup_guard(std::vector<std::string> files) : files_(std::move(files)) {}

    ~cleanup_guard()
    {
        for (const std::string& file : files_) {
            std::error_code error;
            std::filesystem::remove(file, error);
        }
    }

private:
    std::vector<std::string> files_;
};

inline void run_command(const std::string& command)
{
    const int exit_code = std::system(command.c_str());
    ASSERT_EQ(exit_code, 0) << "Command failed: " << command;
}

inline std::vector<double> read_csv_values(const std::string& path)
{
    std::ifstream input(path);
    if (!input.is_open()) {
        throw std::runtime_error("Failed to open CSV file `" + path + "`.");
    }

    std::vector<double> values;
    std::string line;
    while (std::getline(input, line)) {
        std::stringstream row(line);
        std::string cell;
        while (std::getline(row, cell, ',')) {
            values.push_back(std::stod(cell));
        }
    }

    return values;
}

inline void expect_csv_file_near(const std::string& reference_path, const std::string& test_path, double tolerance)
{
    const std::vector<double> reference_values = read_csv_values(reference_path);
    const std::vector<double> test_values = read_csv_values(test_path);

    ASSERT_EQ(reference_values.size(), test_values.size()) << "File size mismatch for " << test_path;
    for (std::size_t i = 0; i < reference_values.size(); ++i) {
        EXPECT_NEAR(test_values[i], reference_values[i], tolerance) << "Mismatch at index " << i << " for " << test_path;
    }
}

#ifdef USE_HDF
inline std::vector<double> read_hdf_dataset(const std::string& path, const std::string& dataset_name)
{
    H5::H5File file(path.c_str(), H5F_ACC_RDONLY);
    H5::DataSet dataset = file.openDataSet(dataset_name);
    H5::DataSpace dataspace = dataset.getSpace();

    hsize_t dims[1] = { 0 };
    const int rank = dataspace.getSimpleExtentDims(dims);
    if (rank != 1) {
        throw std::runtime_error("Expected one-dimensional dataset `" + dataset_name + "`.");
    }

    std::vector<double> values(static_cast<std::size_t>(dims[0]));
    dataset.read(values.data(), H5::PredType::NATIVE_DOUBLE);
    return values;
}

inline void expect_hdf_dataset_near(
    const std::string& reference_path,
    const std::string& test_path,
    const std::string& dataset_name,
    double tolerance)
{
    const std::vector<double> reference_values = read_hdf_dataset(reference_path, dataset_name);
    const std::vector<double> test_values = read_hdf_dataset(test_path, dataset_name);

    ASSERT_EQ(reference_values.size(), test_values.size()) << "Dataset size mismatch for " << dataset_name;
    for (std::size_t i = 0; i < reference_values.size(); ++i) {
        EXPECT_NEAR(test_values[i], reference_values[i], tolerance)
            << "Mismatch in dataset `" << dataset_name << "` at index " << i;
    }
}
#endif

} // namespace test_utils
