include(tools)

include(FindGit)
find_package(Git)

assertVariableSet(matrixName matrixTestName Git_FOUND)

message(STATUS "CMAKE_CURRENT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "PROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}")
message(STATUS "CMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}")
include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest
    GIT_TAG release-1.12.1
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gtest
    INSTALL_COMMAND ""
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable_With_Check(googletest)

add_executable(${matrixTestName} "")

#file(GLOB_RECURSE myHeaders LIST_DIRECTORIES true *.h)
#file(GLOB_RECURSE mySources LIST_DIRECTORIES true *.cpp)

set(myHeaders ArrayTest.h MatrixTest.h)
set(mySources ArrayTest.cpp MatrixTest.cpp SolverMatrixTest.cpp)
target_sources (${matrixTestName} PRIVATE ${myHeaders} PRIVATE ${mySources})

target_link_libraries(${matrixTestName} ${matrixName} gtest_main)

add_test(NAME ${matrixTestName} COMMAND ${matrixTestName})
