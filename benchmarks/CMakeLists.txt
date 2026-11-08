# The framework registers this directory and the aggregate benchmarks target.
# Link production objects plus any Logging objects not already bundled in
# ariadne-constraint-solving; do not duplicate transitive object libraries.
get_target_property(_benchmark_objects logging ARIADNE_PROJECT_OBJECT_TARGETS)
get_target_property(_constraint_objects constraint_solving ARIADNE_PROJECT_OBJECT_TARGETS)
if(NOT _benchmark_objects OR NOT _constraint_objects)
    message(FATAL_ERROR "Missing registered ConstraintSolving/Logging objects")
endif()
foreach(_object IN LISTS _constraint_objects)
    list(REMOVE_ITEM _benchmark_objects "${_object}")
endforeach()

add_executable(benchmark_constraint_reduce benchmark_constraint_reduce.cpp)
target_link_libraries(benchmark_constraint_reduce PRIVATE
    ariadne-constraint-solving constraint_solving ariadne-benchmarking)
foreach(_object IN LISTS _benchmark_objects)
    target_sources(benchmark_constraint_reduce PRIVATE $<TARGET_OBJECTS:${_object}>)
endforeach()
add_dependencies(benchmarks benchmark_constraint_reduce)
