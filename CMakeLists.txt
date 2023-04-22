cmake_minimum_required(VERSION 3.14)

project(Matrix VERSION 0.0.1) #name of your project

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(ExternalProject)
ExternalProject_Add(
    gtest
    GIT_REPOSITORY https://github.com/google/googletest
    GIT_TAG release-1.12.1
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gtest
    INSTALL_COMMAND ""
)
ExternalProject_Get_Property(gtest source_dir binary_dir)

add_executable(runUnitTests ArrayTest.cpp BasicMatrixTest.cpp MatrixTest.cpp Matrix.cpp BasicMatrix.cpp Array.cpp)
add_library(Matrix Matrix.cpp BasicMatrix.cpp Array.cpp)
add_dependencies(runUnitTests gtest)
message(source_dir="${source_dir}")
include_directories(${source_dir}/googletest/include)
message(binary_dir="${binary_dir}")
target_link_libraries(runUnitTests ${binary_dir}/lib/libgtest.a)
target_link_libraries(runUnitTests ${binary_dir}/lib/libgtest_main.a)

enable_testing()
add_test(NAME runUnitTests COMMAND runUnitTests)
