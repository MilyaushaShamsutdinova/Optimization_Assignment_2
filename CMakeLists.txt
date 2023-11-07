cmake_minimum_required(VERSION 3.22)
project(Optimization_Assignment_2)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")

include_directories(headers)

file(GLOB SOURCES "source/*.cpp")
file(GLOB HEADERS "headers/*.hpp")

add_executable(Optimization_Assignment_2
        ${SOURCES}
        ${HEADERS}
        main.cpp)
