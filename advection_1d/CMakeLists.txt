cmake_minimum_required(VERSION 3.20)
project(dg1d LANGUAGES CXX)

# ── C++ standard ─────────────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ── Find TNL ──────────────────────────────────────────────────────────────────
# TNL installs a TNLConfig.cmake alongside its headers.
# REQUIRED means CMake aborts with a clear error if TNL is not found.
find_package(TNL REQUIRED)
find_package(Eigen3 REQUIRED)
# ── Executable ────────────────────────────────────────────────────────────────
add_executable(dg1d src/main.cpp)

# Our own headers live in src/
target_include_directories(dg1d PRIVATE src)

# Linking against TNL::TNL automatically adds TNL's include paths and any
# compile definitions it needs (e.g. -DHAVE_ZLIB if zlib was found at
# TNL install time). Nothing extra is needed because TNL is header-only.
target_link_libraries(dg1d
  TNL::TNL
  Eigen3::Eigen
)

# ── Per-configuration compile options ────────────────────────────────────────
target_compile_options(dg1d PRIVATE
    $<$<CONFIG:Release>:-O3 -march=native>
    $<$<CONFIG:Debug>:  -O0 -g>
    $<$<CONFIG:Diagnostics>:  -O3 -march=native -DDIAGNOSTICS >
)
