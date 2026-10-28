execute_process(
    COMMAND "${Python3_EXECUTABLE}" -m pytest --version
    RESULT_VARIABLE PYTEST_result
    OUTPUT_QUIET
    ERROR_QUIET
)

if(PYTEST_result EQUAL 0)
    set(PYTHON_TESTS ON)
    message(STATUS "pytest package found for ${Python3_EXECUTABLE}; Python tests will be performed.")
else()
    if(PYTHON_ENFORCE_TESTS)
        message(FATAL_ERROR "pytest is not available for ${Python3_EXECUTABLE}; can not perform Python tests.")
    else()
        message(STATUS "pytest is not available for ${Python3_EXECUTABLE}; Python tests will not be performed.")
    endif()
endif()

if(PYTHON_TESTS)
    add_custom_target(
        test_python_constraint_solving
        COMMAND "${Python3_EXECUTABLE}" -m pytest "${CMAKE_CURRENT_SOURCE_DIR}/test_constraint_solving.py"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    )
    add_test(
        NAME test_python_constraint_solving
        COMMAND "${Python3_EXECUTABLE}" -m pytest "${CMAKE_CURRENT_SOURCE_DIR}/test_constraint_solving.py"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    )
    set_property(TEST test_python_constraint_solving PROPERTY
        ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:pyariadne>"
    )
    set_property(TEST test_python_constraint_solving PROPERTY LABELS "python")

    add_dependencies(test_python_constraint_solving pyariadne)
    add_dependencies(tests pyariadne)
endif()
