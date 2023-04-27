cmake_minimum_required(VERSION 3.14)

set(libName Matrix)
set(testLibName ${libName}Test)

project(${libName} VERSION 0.0.4)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CODE_COVERAGE ON)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/lib;${PROJECT_BINARY_DIR}/lib/cmake")

add_subdirectory(lib)

include(code-coverage)
add_code_coverage()

add_subdirectory(src)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  enable_testing()
  add_subdirectory(test)
  target_code_coverage(${testLibName} EXCLUDE build/_deps/googletest-src/* test/*) 
endif()
