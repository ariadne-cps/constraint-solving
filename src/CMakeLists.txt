add_library(CONSTRAINT_SOLVING_SRC OBJECT
    constraint_solver.cpp
)
target_include_directories(CONSTRAINT_SOLVING_SRC PRIVATE
    "${PROJECT_SOURCE_DIR}/include/constraint-solving"
)
