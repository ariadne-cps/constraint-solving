# Loaded only for standalone builds by setup_standalone_project_tests().
if(NOT TARGET benchmarks)
    add_custom_target(benchmarks)
endif()

add_executable(benchmark_constraint_reduce
    benchmark_constraint_reduce.cpp)
target_link_libraries(benchmark_constraint_reduce
    PRIVATE ariadne-constraint-solving constraint_solving logging)
set_target_properties(benchmark_constraint_reduce PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
add_dependencies(benchmarks benchmark_constraint_reduce)
