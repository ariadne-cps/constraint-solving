#!/usr/bin/env python3
from pyariadne import *

def exact_box(bounds):
    return FloatDPExactBox([
        {exact(lo): exact(hi)} for lo, hi in bounds
    ])

def tutorial_constraint_solving():
    solver = ConstraintSolver()
    x0 = ValidatedScalarMultivariateFunction.coordinate(2, 0)
    x1 = ValidatedScalarMultivariateFunction.coordinate(2, 1)
    function = ValidatedVectorMultivariateFunction([x0 + x1])
    domain = exact_box([(0, 2), (0, 2)])
    bounds = exact_box([(0.875, 1.125)])

    rejected, contracted = solver.reduce_box(domain, function, bounds)
    print("rejected:", rejected, "reduced domain:", contracted)

    constraints = [ConstraintPropagationConstraint(
        x0 + x1, bounds[0]
    )]
    sweep_rejected, swept, sweep_stats = solver.propagate_hull_once(
        domain, constraints
    )
    print("one hull pass: rejected:", sweep_rejected,
          "box:", swept, "contractor calls:", sweep_stats.hull_contractor_calls)

    propagation_rejected, propagated, stats = solver.propagate_precompiled(
        domain, constraints
    )
    print("propagation: rejected:", propagation_rejected,
          "box:", propagated, "hull rounds:", stats.hull_rounds,
          "shaving rounds:", stats.shaving_rounds)

    y = ValidatedScalarMultivariateFunction.coordinate(1, 0)
    boundary = exact_box([(0, 0)])
    strict = [ConstraintPropagationConstraint(
        y, exact_box([(0, 1)])[0], True, False
    )]
    strict_rejected, strict_domain, _ = solver.propagate_precompiled(
        boundary, strict, monotone=False, shaving=False, hull=False
    )
    print("strict x>0 at x=0: rejected:", strict_rejected,
          "domain:", strict_domain)

    # Feasibility search, validated point checks and rejection.
    solver = ConstraintSolver()
    x0 = ValidatedScalarMultivariateFunction.coordinate(2, 0)
    x1 = ValidatedScalarMultivariateFunction.coordinate(2, 1)
    function = ValidatedVectorMultivariateFunction([x0 + x1])
    domain = exact_box([(0, 1), (0, 1)])
    codomain = exact_box([(0.19, 0.21)])

    status, candidate = solver.feasible(domain, function, codomain)
    print("feasibility:", status, "candidate:", candidate)

    dp = DoublePrecision()
    witness = FloatDPPoint([FloatDP(exact(0.1), dp), FloatDP(exact(0.1), dp)])
    print("point feasibility:",
          solver.check_feasibility(domain, function, codomain, witness))

    impossible = exact_box([(3, 4)])
    inconsistent, _ = solver.feasible(domain, function, impossible)
    print("inconsistent feasibility:", inconsistent)


if __name__ == "__main__":
    tutorial_constraint_solving()
