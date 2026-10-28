#include "pybind11.hpp"

void constraint_solving_submodule(pybind11::module& module);

void constraint_solving_module(pybind11::module& module) {
    constraint_solving_submodule(module);
}
