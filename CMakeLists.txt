cmake_minimum_required(VERSION 3.22)
project(Optimization_Assignment_2)

set(CMAKE_CXX_STANDARD 14)

include_directories(headers)

add_executable(Optimization_Assignment_2 main.cpp
        headers/Matrix.hpp
        headers/Vector.hpp
        headers/InteriorPoint.hpp
        headers/SimplexMethod.hpp
        source/Matrix.cpp
        source/Vector.cpp
        source/InteriorPoint.cpp
        source/SimplexMethod.cpp)
