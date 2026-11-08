cmake_minimum_required(VERSION 3.15)

project(ConstraintSolving VERSION 1.0)

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

if(EXISTS "${PROJECT_SOURCE_DIR}/.git")
    # Compare direct occurrences; child commits guarantee internal consistency.
    require_same_dependency_commit(Build CHAINS
        "build" "nonlinear-programming|build" "testing|build")
endif()

setup_project()

include_directories(${PROJECT_SOURCE_DIR}/include)

add_subdirectory(ss)

add_library(constraint_solving INTERFACE)
target_include_directories(constraint_solving INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/ariadne-constraint-solving>
)
target_link_libraries(constraint_solving INTERFACE nonlinear_programming)

add_subdirectory(src)
target_link_libraries(CONSTRAINT_SOLVING_SRC PRIVATE constraint_solving)

ariadne_register_public_headers(
    TARGET constraint_solving
    DIRECTORIES "${PROJECT_SOURCE_DIR}/include/constraint-solving"
    FILES "${PROJECT_SOURCE_DIR}/include/ariadne-constraint-solving.hpp"
)

setup_project_library(
    TARGET ariadne-constraint-solving
    OBJECTS CONSTRAINT_SOLVING_SRC
    LINK_LIBRARIES constraint_solving
)

setup_standalone_project_tests(
    TARGET ariadne-constraint-solving
    TEST_DIRECTORY tests
    TEST_TARGET tests
    SOURCES
        "${PROJECT_SOURCE_DIR}/include"
        "${PROJECT_SOURCE_DIR}/src"
)

# Reuse the benchmarks aggregator provided by dependencies, as in paving.
# The diagnostic executable is opt-in and excluded from the library/test build.
if(NOT TARGET benchmarks)
    add_custom_target(benchmarks)
endif()

add_executable(benchmark_constraint_reduce
    benchmarks/benchmark_constraint_reduce.cpp)
target_link_libraries(benchmark_constraint_reduce
    PRIVATE ariadne-constraint-solving constraint_solving logging)
set_target_properties(benchmark_constraint_reduce PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/benchmarks")
add_dependencies(benchmarks benchmark_constraint_reduce)

if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(python)
endif()

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    ariadne_install_dependency_bundle(
        TARGET constraint_solving
        DESTINATION include/ariadne-constraint-solving
    )
endif()
