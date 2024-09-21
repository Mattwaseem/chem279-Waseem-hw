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
