cmake_minimum_required(VERSION 3.10)
project(advection)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# define sources
set(SOURCES main.cpp src/ )
# Create executable
add_executable(advection ${SOURCES})
