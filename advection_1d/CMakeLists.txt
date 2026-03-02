cmake_minimum_required(VERSION 3.10)
project(advection CXX)

# set standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# locate TNL
find_package(TNL REQUIRED)

# define sources
set(SOURCES main.cpp src/ )

# create executable
add_executable(advection ${SOURCES})

# link the TNL library
target_link_libraries(advection PUBLIC TNL::TNL)
