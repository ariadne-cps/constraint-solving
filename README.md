# Ariadne Constraint Solving

[![Unix Status](https://github.com/ariadne-cps/constraint-solving/actions/workflows/unix.yml/badge.svg)](https://github.com/ariadne-cps/constraint-solving/actions/workflows/unix.yml)
[![Windows Status](https://github.com/ariadne-cps/constraint-solving/actions/workflows/win.yml/badge.svg)](https://github.com/ariadne-cps/constraint-solving/actions/workflows/win.yml)
[![Coverage Status](https://github.com/ariadne-cps/constraint-solving/actions/workflows/coverage.yml/badge.svg)](https://github.com/ariadne-cps/constraint-solving/actions/workflows/coverage.yml)

Ariadne Constraint Solving is a standalone C++20 library for rigorous numerical
constraint solving. It restores the unified `ConstraintSolver` and
`ConstraintSolverInterface` from `ariadne-cps/solving`, based on the
pre-split commit [`4502fea11dd7b50306d2113af40fcc60db14ed0c`](https://github.com/ariadne-cps/solving/commit/4502fea11dd7b50306d2113af40fcc60db14ed0c).

The solver provides validated interval contraction (hull, shaving, monotone
reduction, cached/incremental propagation), numerical feasibility checks and
point validation, and interior-point-based candidate search. Inconsistent
strict inequalities can be rejected even if the closed interval box is not empty.

## Dependencies

The only direct **computational** dependency is
[nonlinear-programming](https://github.com/ariadne-cps/nonlinear-programming).
Its Function dependency is deliberately inherited transitively rather than
declared or checked out independently. The repository also pins
[configuration](https://github.com/ariadne-cps/configuration) for CMake.

## Public API

`#include "ariadne-constraint-solving.hpp"`, or
`#include "constraint-solving/constraint_solver.hpp"`.

The single public `ConstraintSolver` offers `reduce`, `propagate`,
`propagate_hull_once`, `hull_reduce`, `box_reduce`,
`monotone_reduce`, `lyapunov_reduce`, `feasible`,
`check_feasibility` and `split`. The numerical algorithms are maintained as
separate operations with distinct correctness contracts. A numerical candidate
is not automatically a certified feasibility witness.

## Build and test

```sh
git clone --recurse-submodules https://github.com/ariadne-cps/constraint-solving.git
cd constraint-solving
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel --target tests
ctest --test-dir build --output-on-failure
```

To build and install the standalone C++ library and Python bindings:

```sh
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release --parallel --target ariadne-constraint-solving bindings
cmake --install build-release --prefix /path/to/install
```

The C++ and Python examples each demonstrate **both propagation and feasibility**:

- `tutorials/tutorial_constraint_solving/tutorial_constraint_solving.cpp`
- `python/tutorials/tutorial_constraint_solving.py`

The original unified regression suite is in `tests/test_constraint_solver.cpp`.
Python tests exercise hull/interval propagation, strict-bound rejection,
feasibility, and validated point checking.

## Scope

This is a standalone extraction, not yet a migration of consumers. `set`
and `smt` still require separate dependency and include-path changes.
The existing epsilon-SMT result contract is not changed by this extraction.
