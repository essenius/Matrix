include(tools)

assertVariableSet(libName)

message(STATUS "In Matrix ${libName} src: CMAKE_PROJECT_NAME=${CMAKE_PROJECT_NAME}")

#get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

add_library(${libName} "")

file(GLOB_RECURSE myHeaders LIST_DIRECTORIES true *.h)
file(GLOB_RECURSE mySources LIST_DIRECTORIES true *.cpp)
target_sources (${libName} PUBLIC ${myHeaders} PRIVATE ${mySources})
target_include_directories(${libName} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ${libName} DESTINATION ${CMAKE_BINARY_DIR})
