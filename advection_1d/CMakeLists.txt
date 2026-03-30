cmake_minimum_required(VERSION 3.20)
project(dg_solver LANGUAGES CXX)

# set standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# find TNL and Eigen3
find_package(TNL REQUIRED)
find_package(Eigen3 REQUIRED)

# create executable
add_executable(dg_solver src/main.cpp)

# DGM headers live in headers/
target_include_directories(dg1d PRIVATE headers)

# link TNL and Eigen3
target_link_libraries(dg_solver
  TNL::TNL
  Eigen3::Eigen
)

# compile options: Release, Debug, Diagnostics
target_compile_options(dg_solver PRIVATE
    $<$<CONFIG:Release>:-O3 -march=native>
    $<$<CONFIG:Debug>:  -O0 -g>
    $<$<CONFIG:Diagnostics>:  -O3 -march=native -DDIAGNOSTICS >
)
