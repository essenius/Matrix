include(tools)

assertVariableSet(libName)

message(STATUS "${libName} src: CMAKE_PROJECT_NAME=${CMAKE_PROJECT_NAME}")

add_library(${libName} "")

file(GLOB_RECURSE myHeaders LIST_DIRECTORIES true *.h)
file(GLOB_RECURSE mySources LIST_DIRECTORIES true *.cpp)
target_sources (${libName} PUBLIC ${myHeaders} PRIVATE ${mySources})
target_include_directories(${libName} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ${libName} DESTINATION ${CMAKE_BINARY_DIR})
