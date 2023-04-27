include(tools)

include(FindGit)
find_package(Git)

assertVariableSet(testLibName Git_FOUND)

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

add_executable(${testLibName} "")

file(GLOB_RECURSE myHeaders LIST_DIRECTORIES true *.h)
file(GLOB_RECURSE mySources LIST_DIRECTORIES true *.cpp)
target_sources (${testLibName} PUBLIC ${myHeaders} PRIVATE ${mySources})

target_link_libraries(${testLibName} ${CMAKE_PROJECT_NAME} gtest_main)

add_test(NAME ${testLibName} COMMAND ${testLibName})
