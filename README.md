# Ariadne Constraint Solving

[![macOS Status](https://github.com/ariadne-cps/constraint-solving/actions/workflows/macos.yml/badge.svg)](https://github.com/ariadne-cps/constraint-solving/actions/workflows/macos.yml)
[![Ubuntu Status](https://github.com/ariadne-cps/constraint-solving/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/ariadne-cps/constraint-solving/actions/workflows/ubuntu.yml)
[![Windows Status](https://github.com/ariadne-cps/constraint-solving/actions/workflows/win.yml/badge.svg)](https://github.com/ariadne-cps/constraint-solving/actions/workflows/win.yml)
[![Sanitizers Status](https://github.com/ariadne-cps/constraint-solving/actions/workflows/sanitizers.yml/badge.svg)](https://github.com/ariadne-cps/constraint-solving/actions/workflows/sanitizers.yml)
[![Coverage Status](https://github.com/ariadne-cps/constraint-solving/actions/workflows/coverage.yml/badge.svg)](https://github.com/ariadne-cps/constraint-solving/actions/workflows/coverage.yml)
[![codecov](https://codecov.io/gh/ariadne-cps/constraint-solving/branch/main/graph/badge.svg)](https://codecov.io/gh/ariadne-cps/constraint-solving)

Ariadne Constraint Solving is a standalone C++20 library for validated numerical constraint solving. `ConstraintSolver` provides hull propagation, coordinate shaving, derivative-assisted monotone contraction, incremental propagation with cached procedures, and validated feasibility and infeasibility checks using nonlinear interior-point methods.

## Repository layout

- `include/constraint-solving/`: public `ConstraintSolver`, `ConstraintSolverInterface`, `ConstraintPropagationConstraint` and `ConstraintPropagationStatistics` APIs.
- `src/`: constraint-solving implementation.
- `tests/`: C++ regression tests for contraction, propagation, feasibility and validation.
- `python/bindings/`, `python/tests/`, `python/tutorials/`: Python interface, tests and tutorial.
- `tutorials/tutorial_constraint_solving/`: standalone installed C++ tutorial.

The umbrella header is `ariadne-constraint-solving.hpp`.

## Dependencies

- [ariadne-cps/nonlinear-programming](https://github.com/ariadne-cps/nonlinear-programming): nonlinear feasibility and interior-point methods.
- [ariadne-cps/configuration](https://github.com/ariadne-cps/configuration): common CMake build infrastructure.

## Build and test

```sh
git clone --recurse-submodules https://github.com/ariadne-cps/constraint-solving.git
cd constraint-solving
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel --target tests
ctest --test-dir build --output-on-failure
```

To install the C++ library and Python bindings (with a supported Python 3.11+ environment):

```sh
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release --parallel --target ariadne-constraint-solving bindings
cmake --install build-release --prefix /path/to/install
```

See `tutorials/tutorial_constraint_solving/tutorial_constraint_solving.cpp` and `python/tutorials/tutorial_constraint_solving.py` for examples of both propagation and feasibility.

**Validation semantics:** propagation returns true when inconsistency is proved, even if the closed interval box remains nonempty for strict inequalities. Feasibility results may be indeterminate, and a numerical candidate is not automatically a certified existence witness.

## License

GNU General Public License v3.0 (or later where noted in source files).
