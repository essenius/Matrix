set(myHeaders Array.h Matrix.h SolverMatrix.h)
set(mySources Array.cpp Matrix.cpp SolverMatrix.cpp)

if (ESP_PLATFORM AND DEFINED ENV{IDF_PATH})
    idf_component_register(SRCS ${mySources}
                    INCLUDE_DIRS ".")
else()
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


    target_sources (${matrixName} PUBLIC ${myHeaders} PRIVATE ${mySources})
    target_include_directories(${matrixName} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

    message(STATUS "CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}")
    message(STATUS "CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}")
    message(STATUS "CMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}")
    message(STATUS "CMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}")

    install(TARGETS ${matrixName} DESTINATION lib)
    install(FILES ${myHeaders} DESTINATION include) 
endif()