cmake_minimum_required(VERSION 3.14)

set(matrixName Matrix)
set(matrixTestName ${matrixName}Test)

project(${matrixName} VERSION 0.2.7 LANGUAGES CXX)

# We only want to enable testing and code coverage for the top level project
# On Windows, these are also only used when building a vcpkg package.
if ((CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR) AND NOT DEFINED VCPKG_TARGET_TRIPLET)
  set(TOP_LEVEL ON)
#  if(WIN32)
#    set(CMAKE_GENERATOR "MinGW Makefiles")
#  endif()
else()
  set(TOP_LEVEL OFF)
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/lib;${PROJECT_BINARY_DIR}/lib/cmake")

# haven't made code coverage work for Windows yet
if (TOP_LEVEL AND NOT DEFINED WIN32) 
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
  if (CODE_COVERAGE)
    setup_target_code_coverage(${matrixTestName} EXCLUDE build/_deps/googletest-src/* test/*)
  endif()
endif()

# On Windows, build as follows:
# cmake -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX="../install" ..
# cmake --build . --target install
