include_directories(${Python3_INCLUDE_DIRS})

get_target_property(PYTHON_COMMON_SOURCE_DIR pyariadne-python-common SOURCE_DIR)

add_library(pyariadne-constraint_solving-submodule-obj OBJECT
    constraint_solving_submodule.cpp
)
target_include_directories(pyariadne-constraint_solving-submodule-obj PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${PYTHON_COMMON_SOURCE_DIR}/include"
    $<TARGET_PROPERTY:pyariadne-nonlinear_programming,INTERFACE_INCLUDE_DIRECTORIES>
)
target_link_libraries(pyariadne-constraint_solving-submodule-obj PRIVATE
    pyariadne-python-common
    pybind11::module
    constraint_solving
)
set_target_properties(pyariadne-constraint_solving-submodule-obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(MSVC)
    target_compile_options(pyariadne-constraint_solving-submodule-obj PRIVATE
        /bigobj
        /wd4702
    )
endif()

add_library(pyariadne-constraint_solving INTERFACE)
target_sources(pyariadne-constraint_solving INTERFACE
    $<TARGET_OBJECTS:pyariadne-constraint_solving-submodule-obj>
)
target_include_directories(pyariadne-constraint_solving INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<TARGET_PROPERTY:pyariadne-nonlinear_programming,INTERFACE_INCLUDE_DIRECTORIES>
)
target_link_libraries(pyariadne-constraint_solving INTERFACE
    pyariadne-nonlinear_programming
)

ariadne_add_python_module(
    NAME constraint_solving
    MODULE_SOURCE constraint_solving_module.cpp
)
