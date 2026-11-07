# Diagnostic only: build explicitly, not as part of ordinary tests/CI.
if(NOT TARGET benchmarks)
    add_custom_target(benchmarks)
endif()

add_executable(benchmark_constraint_reduce benchmark_constraint_reduce.cpp)
target_link_libraries(benchmark_constraint_reduce
    PRIVATE ariadne-constraint-solving constraint_solving)
add_dependencies(benchmarks benchmark_constraint_reduce)
