if(NOT TARGET nonlinear_programming)
    ensure_submodule_initialized(nonlinear-programming)
    add_subdirectory(nonlinear-programming)
endif()
