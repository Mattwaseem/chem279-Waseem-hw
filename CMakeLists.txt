cmake_minimum_required(VERSION 3.10)

# Project Name
project(LJ_Project)

# Set C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include directories
include_directories(include)

# Set the output directory for executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Add the executable with source files
add_executable(ljpotenial src/main.cpp src/file_io.cpp src/LJ_PE.cpp)

# Custom target to clean the bin directory
add_custom_target(
    clean-bin
    COMMAND ${CMAKE_COMMAND} -E rm -rf ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/*
    COMMENT "Cleaning up bin directory"
)

# Custom target to clean all build files (including bin)
add_custom_target(
    clean-all
    COMMAND ${CMAKE_COMMAND} -E rm -rf ${CMAKE_BINARY_DIR}/* ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/*
    COMMENT "Cleaning all build and bin directories"
)

# Step 1: Generate the output file by running ljpotenial executable
add_custom_command(
    OUTPUT truncation_error.txt
    COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ljpotenial
    COMMENT "Running ljpotenial to generate truncation_error.txt"
)

# Step 2: Run the Python script to generate the plot
add_custom_command(
    OUTPUT plot.png
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/plot_truncation_error.py
    DEPENDS truncation_error.txt
    COMMENT "Running Python script to generate plot.png"
)

# Step 3: Add a custom target to make sure the Python script runs after the C++ executable
add_custom_target(run_plot ALL DEPENDS plot.png)

# Custom target to build everything, including the plot
add_custom_target(
    all_plot
    DEPENDS run_plot ljpotenial
    COMMENT "Building ljpotenial and generating the plot"
)
