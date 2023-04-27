cmake_minimum_required(VERSION 3.14)

set(libName Matrix)
set(testLibName ${libName}Test)

if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(TOP_LEVEL ON)
else()
  set(TOP_LEVEL OFF)
endif()

project(${libName} VERSION 0.0.4)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/lib;${PROJECT_BINARY_DIR}/lib/cmake")

if (TOP_LEVEL) 
  message(STATUS "Top level project - enabling code coverage")
  add_subdirectory(lib)
  set(CODE_COVERAGE ON)
  include(code-coverage)
  add_code_coverage()
endif()

add_subdirectory(src)

if(TOP_LEVEL)
  message(STATUS "Top level project - enabling tests")
  enable_testing()
  add_subdirectory(test)
  target_code_coverage(${testLibName} EXCLUDE build/_deps/googletest-src/* test/*) 
endif()
