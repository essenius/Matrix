include(tools)

include(FindGit)
find_package(Git)

assertVariableSet(matrixName matrixTestName Git_FOUND)

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

message(CMAKE_BINARY_DIR=${CMAKE_BINARY_DIR})

add_executable(${matrixTestName} "")

file(GLOB_RECURSE myHeaders LIST_DIRECTORIES true *.h)
file(GLOB_RECURSE mySources LIST_DIRECTORIES true *.cpp)
target_sources (${matrixTestName} PUBLIC ${myHeaders} PRIVATE ${mySources})

#message(STATUS "CMAKE_PROJECT_NAME=${CMAKE_PROJECT_NAME} = ${matrixName}")
target_link_libraries(${matrixTestName} ${matrixName} gtest_main)

add_test(NAME ${matrixTestName} COMMAND ${matrixTestName})
