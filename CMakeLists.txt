# This file is part of GRAMPC-S - (https://github.com/grampc/grampc-s)
#
# GRAMPC-S -- A software framework for stochastic model predictive control (SMPC)
#
# Copyright 2024 Daniel Landgraf, Andreas Voelz, Knut Graichen.
# All rights reserved.
#
# GRAMPC-S is distributed under the BSD-3-Clause license, see LICENSE.txt


cmake_minimum_required(VERSION 3.13)
project(grampc-s VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "-march=native")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# link_directories should treat paths as relative to the source directory
cmake_policy(SET CMP0015 NEW)

# Eigen 3.4.90 is required
find_package(Eigen3 3.4.90 REQUIRED NO_MODULE)

include_directories(
	include
	libs/grampc/include
	${EIGEN3_INCLUDE_DIR}
)

# Add library and link Eigen
add_library(grampc-s)
target_link_libraries(grampc-s Eigen3::Eigen)

add_subdirectory(examples)
add_subdirectory(libs)
add_subdirectory(src)