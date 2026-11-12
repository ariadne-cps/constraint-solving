cmake_minimum_required(VERSION 3.15)

project(AriadneConstraintSolving VERSION 3.0)

option(WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)

if(NOT TARGET ariadne-build)
    find_package(Git REQUIRED)
    # Keep initial Build checkout output in English.
    execute_process(
        COMMAND "${CMAKE_COMMAND}" -E env "LC_ALL=C"
                "${GIT_EXECUTABLE}" -C "${CMAKE_CURRENT_SOURCE_DIR}"
                submodule update --init -- ss/bui
        RESULT_VARIABLE _result)
    if(_result)
        message(FATAL_ERROR "Cannot initialize build submodule")
    endif()
    include("${CMAKE_CURRENT_SOURCE_DIR}/ss/bui/cmake/ProjectOptions.cmake")
endif()

require_same_dependency_commit(Build CHAINS
    "build" "nonlinear-programming|build")

require_same_dependency_commit(Logging CHAINS
    "nonlinear-programming|logging" "benchmarking|logging")

ariadne_setup_module(constraint-solving)
