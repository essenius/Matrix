set (testName ${CMAKE_PROJECT_NAME}Test)

include(FindGit)
find_package(Git)

if (NOT Git_FOUND)
	message(FATAL_ERROR "Git not found!")
endif ()

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

FetchContent_MakeAvailable(googletest)

message(CMAKE_BINARY_DIR=${CMAKE_BINARY_DIR})

file(GLOB_RECURSE myHeaders LIST_DIRECTORIES true *.h)
file(GLOB_RECURSE mySources LIST_DIRECTORIES true *.cpp)
target_sources (${testName} PUBLIC ${myHeaders} PRIVATE ${mySources})

target_link_libraries(${testName} ${CMAKE_PROJECT_NAME} gtest_main)

message(source_dir="${source_dir}")
message(binary_dir="${binary_dir}")

add_test(NAME ${testName} COMMAND ${testName})
