include(CTest)
add_executable(test_constraint_solver test_constraint_solver.cpp)
target_link_libraries(test_constraint_solver PRIVATE ariadne-constraint-solving constraint_solving)
add_test(NAME test_constraint_solver COMMAND $<TARGET_FILE:test_constraint_solver>)
set_property(TEST test_constraint_solver PROPERTY LABELS "constraint-solving")
if(WIN32)
    set_tests_properties(test_constraint_solver PROPERTIES
        ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:$<TARGET_FILE_DIR:ariadne-constraint-solving>"
    )
endif()
add_custom_target(tests)
add_dependencies(tests test_constraint_solver)
