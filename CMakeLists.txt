#cmake_minimum_required(VERSION 3.16)
cmake_minimum_required(VERSION 3.10)

project(cgshop2022_simple LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED true)

set(SRC ${CMAKE_SOURCE_DIR}/src)
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
include_directories("${SRC}")

set(CMAKE_CXX_FLAGS "-Wall -pthread -O3 -mavx -fopenmp -g")

find_package(Threads)
find_package(X11)

file(GLOB SOURCES_dimacs ${SRC}/dimacs.cpp)
file(GLOB SOURCES_simple ${SRC}/simple.cpp)

set(JSON_BuildTests OFF CACHE INTERNAL "")

add_subdirectory(json)

add_executable(simple ${SOURCES_simple})
add_executable(dimacs ${SOURCES_dimacs})
target_link_libraries(simple PUBLIC nlohmann_json::nlohmann_json)
target_link_libraries(dimacs PUBLIC nlohmann_json::nlohmann_json)
