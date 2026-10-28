from pyariadne import (
    ConstraintSolver, ConstraintPropagationConstraint,
    ConstraintPropagationStatistics, ValidatedConstraint,
    ValidatedScalarMultivariateFunction, ValidatedVectorMultivariateFunction,
    FloatDPExactBox, exact
)


def box(bounds):
    return FloatDPExactBox([{exact(lo): exact(hi)} for lo, hi in bounds])


def test_contract_and_propagate():
    solver = ConstraintSolver()
    x0 = ValidatedScalarMultivariateFunction.coordinate(2, 0)
    x1 = ValidatedScalarMultivariateFunction.coordinate(2, 1)
    function = ValidatedVectorMultivariateFunction([x0 + x1])
    domain = box([(0, 2), (0, 2)])
    bounds = box([(0.875, 1.125)])

    rejected, reduced = solver.reduce_box(domain, function, bounds)
    assert not rejected
    assert reduced is not None

    constraints = [ConstraintPropagationConstraint(x0 + x1, bounds[0])]
    rejected, reduced, stats = solver.propagate_precompiled(
        domain, constraints, monotone=False, shaving=True, hull=True
    )
    assert not rejected
    assert isinstance(stats, ConstraintPropagationStatistics)
    assert stats.hull_rounds >= 1
    assert reduced is not None
    rejected_once, swept, sweep_stats = solver.propagate_hull_once(
        domain, constraints
    )
    assert not rejected_once
    assert sweep_stats.hull_contractor_calls == 1
    assert swept is not None


def test_strict_boundary_rejection():
    solver = ConstraintSolver()
    x = ValidatedScalarMultivariateFunction.coordinate(1, 0)
    strict = [ConstraintPropagationConstraint(x, box([(0, 1)])[0], True, False)]
    rejected, unchanged, stats = solver.propagate_precompiled(
        box([(0, 0)]), strict, monotone=False, shaving=False, hull=False
    )
    assert rejected
    assert unchanged is not None
    assert stats.hull_contractor_calls == 0


from pyariadne import *

def exact_box(bounds):
    return FloatDPExactBox([{exact(a): exact(b)} for a,b in bounds])

def test_point_validation_and_infeasibility():
    solver=ConstraintSolver()
    x=ValidatedScalarMultivariateFunction.coordinate(1,0)
    function=ValidatedVectorMultivariateFunction([x])
    domain=exact_box([(0,2)])
    codomain=exact_box([(0,1)])
    point=FloatDPPoint([FloatDP(exact(0.5),DoublePrecision())])
    assert definitely(solver.check_feasibility(domain,function,codomain,point))
    impossible=exact_box([(3,4)])
    status,_=solver.feasible(domain,function,impossible)
    assert not possibly(status)
