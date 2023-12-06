include(tools)

assertVariableSet(matrixName)

file(APPEND ${PROJECT_BINARY_DIR}/output.txt "In Matrix ${matrixName} src: CMAKE_PROJECT_NAME=${CMAKE_PROJECT_NAME}")

file(APPEND ${PROJECT_BINARY_DIR}/output.txt "CMAKE_CURRENT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}")
file(APPEND ${PROJECT_BINARY_DIR}/output.txt "PROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}")
file(APPEND ${PROJECT_BINARY_DIR}/output.txt "CMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}")

get_cmake_property(_variableNames VARIABLES)
list (SORT _variableNames)
foreach (_variableName ${_variableNames})
    file(APPEND ${PROJECT_BINARY_DIR}/output.txt "${_variableName}=${${_variableName}}\n")
endforeach()

add_library(${matrixName} "")

#file(GLOB_RECURSE myHeaders LIST_DIRECTORIES true *.h)
#file(GLOB_RECURSE mySources LIST_DIRECTORIES true *.cpp)

set(myHeaders Array.h Matrix.h SolverMatrix.h)
set(mySources Array.cpp Matrix.cpp SolverMatrix.cpp)

target_sources (${matrixName} PUBLIC ${myHeaders} PRIVATE ${mySources})
target_include_directories(${matrixName} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

file(APPEND ${PROJECT_BINARY_DIR}/output.txt "Installing to ${CMAKE_CURRENT_BINARY_DIR}")
install(TARGETS ${matrixName} DESTINATION lib) # ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${myHeaders} DESTINATION include) 
