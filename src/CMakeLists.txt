cmake_minimum_required(VERSION 3.14)

project(Matrix VERSION 0.0.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(Matrix SolverMatrix.cpp Matrix.cpp Array.cpp)
