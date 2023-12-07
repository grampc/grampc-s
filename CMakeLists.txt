cmake_minimum_required(VERSION 3.13)
set(CMAKE_CXX_STANDARD 11)

project(grampc-s)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# link_directories should treat paths as relative to the source directory
cmake_policy(SET CMP0015 NEW)

# Use timestamps of the time of extraction
cmake_policy(SET CMP0135 NEW)

include_directories(
	include
	libs/grampc/include
	libs/eigen
)

add_library(grampc-s)

add_subdirectory(examples)
add_subdirectory(libs)
add_subdirectory(src)
