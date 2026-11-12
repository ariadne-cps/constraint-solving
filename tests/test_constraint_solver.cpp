/***************************************************************************
 *            tests/test_constraint_solver.cpp
 *
 *  Copyright  2010-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne Constraint Solving.
 *
 *  Ariadne Constraint Solving is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne Constraint Solving is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Constraint Solving. If not, see <https://www.gnu.org/licenses/>.
 */

/** @file tests/test_constraint_solver.cpp
 *  @brief Test validated constraint propagation and reduction.
 */

#include <fstream>
#include <iostream>
#include <vector>

#include "algebra/vector.hpp"
#include "foundation/logical.hpp"
#include "function/constraint.hpp"
#include "function/function.hpp"
#include "function/procedure.hpp"
#include "function/taylor_function.hpp"
#include "geometry/box.hpp"
#include "logging/logging.hpp"
#include "numeric/numeric.hpp"

#include "ariadne-testing.hpp"
#include "constraint-solving/constraint_solver.hpp"

using namespace std;
using namespace Ariadne;

class TestConstraintSolver : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    Void test_empty_reduce_inequality() {
        List<EffectiveScalarMultivariateFunction> x = EffectiveScalarMultivariateFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{
            {0.0_x, 1.0_x},
            {0.0_x, 1.0_x}
        };
        List<EffectiveConstraint> c = {4 <= 2 * x[0] + x[1]};

        ConstraintSolver propagator;

        ARIADNE_TESTING_NOTE("Executing: propagator.reduce(D,c)");
        propagator.reduce(D, c);
        ARIADNE_TESTING_NOTE("D = " << (D));
        ARIADNE_TESTING_PROVE(D.is_empty());
    }

    Void test_empty_reduce_equality() {
        List<EffectiveScalarMultivariateFunction> x = EffectiveScalarMultivariateFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{
            {0.0_x, 1.0_x},
            {0.0_x, 1.0_x}
        };
        List<EffectiveConstraint> c = {2 * x[0] + x[1] == 4};

        ConstraintSolver propagator;

        ARIADNE_TESTING_NOTE("Executing: propagator.reduce(D,c)");
        propagator.reduce(D, c);
        ARIADNE_TESTING_NOTE("D = " << (D));
        ARIADNE_TESTING_PROVE(D.is_empty());
    }

    Void test_empty_reduce_mixed() {
        List<EffectiveScalarMultivariateFunction> x = EffectiveScalarMultivariateFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{
            {0.0_x, 0.25_x},
            {0.0_x,  2.0_x}
        };
        List<EffectiveConstraint> c = {x[1] <= 1, x[0] + x[1] == 2};

        ConstraintSolver propagator;

        ARIADNE_TESTING_NOTE("Executing: propagator.reduce(D,c)");
        propagator.reduce(D, c);
        ARIADNE_TESTING_NOTE("D = " << (D));
        ARIADNE_TESTING_PROVE(D.is_empty());
    }

    Void test_empty_hull_reduce() {
        List<EffectiveScalarMultivariateFunction> x = EffectiveScalarMultivariateFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{
            {0.0_x, 0.25_x},
            {0.0_x,  2.0_x}
        };
        List<EffectiveConstraint> c = {x[1] <= 1, x[0] + x[1] == 2};

        ConstraintSolver propagator;

        ARIADNE_TESTING_NOTE("Executing: propagator.hull_reduce(D,c[0])");
        propagator.hull_reduce(D, c[0]);
        ARIADNE_TESTING_NOTE("Executing: propagator.hull_reduce(D,c[1])");
        propagator.hull_reduce(D, c[1]);
        ARIADNE_TESTING_NOTE("D = " << (D));
        ARIADNE_TESTING_PROVE(D.is_empty());
    }

    Void test_empty_box_reduce() {
        List<EffectiveScalarMultivariateFunction> x = EffectiveScalarMultivariateFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{
            {0.0_x, 0.25_x},
            {0.0_x,  2.0_x}
        };
        List<EffectiveConstraint> c = {x[1] <= 1, x[0] + x[1] == 2};

        ConstraintSolver propagator;

        ARIADNE_TESTING_NOTE("Executing: propagator.box_reduce(D,c[0],0)");
        propagator.box_reduce(D, c[0], 0);
        ARIADNE_TESTING_NOTE("Executing: propagator.box_reduce(D,c[1],0)");
        propagator.box_reduce(D, c[1], 0);
        ARIADNE_TESTING_NOTE("Executing: propagator.box_reduce(D,c[0],1)");
        propagator.box_reduce(D, c[0], 1);
        ARIADNE_TESTING_NOTE("Executing: propagator.box_reduce(D,c[1],1)");
        propagator.box_reduce(D, c[1], 1);
        ARIADNE_TESTING_NOTE("Executing: propagator.hull_reduce(D,c[0])");
        propagator.hull_reduce(D, c[0]);
        ARIADNE_TESTING_NOTE("Executing: propagator.hull_reduce(D,c[1])");
        propagator.hull_reduce(D, c[1]);
        ARIADNE_TESTING_NOTE("D = " << (D));
        ARIADNE_TESTING_PROVE(D.is_empty());
    }

    Void test_hull_reduce() {
        List<EffectiveScalarMultivariateFunction> x = EffectiveScalarMultivariateFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{
            {0.0_x, 2.0_x},
            {0.0_x, 2.0_x}
        };
        List<EffectiveConstraint> c = {-2 <= 2 * x[0] + x[1] <= 1};

        ConstraintSolver propagator;

        ARIADNE_TESTING_NOTE("Executing: propagator.hull_reduce(D,c[0])");
        propagator.hull_reduce(D, c[0]);
        ARIADNE_TESTING_EXPECT_SAME(D, UpperBoxType({
                                           {0.0_x, 0.5_x},
                                           {0.0_x, 1.0_x}
        }));
    }

    Void test_box_reduce() {
        List<EffectiveScalarMultivariateFunction> x = EffectiveScalarMultivariateFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{
            {0.0_x, 2.0_x},
            {0.0_x, 2.0_x}
        };
        EffectiveConstraint c = (-2 <= 2 * x[0] + x[1] <= 1);

        ConstraintSolver propagator;

        ARIADNE_TESTING_NOTE("Executing: propagator.box_reduce(D,c,0)");
        propagator.box_reduce(D, c, 0);
        ARIADNE_TESTING_EXPECT_SAME(D, UpperBoxType({
                                           {0.0_x, 0.75_x},
                                           {0.0_x,  2.0_x}
        }));
        ARIADNE_TESTING_NOTE("Executing: propagator.box_reduce(D,c,1)");
        propagator.box_reduce(D, c, 1);
        ARIADNE_TESTING_EXPECT_SAME(D, UpperBoxType({
                                           {0.0_x, 0.75_x},
                                           {0.0_x, 1.25_x}
        }));
    }

    Void test_box_reduce_boundary_slices() {
        ConstraintSolver solver;
        auto x = ValidatedScalarMultivariateFunction::coordinates(1);

        // A zero-width variable has no slice to shave; its value is preserved.
        ARIADNE_TESTING_NOTE("[constraint-shaving] singleton domain remains unchanged");
        UpperBoxType singleton = ExactBoxType({
            {0.5_x, 0.5_x}
        });
        ARIADNE_TESTING_PROVE(
            not solver.box_reduce(singleton, x[0], ExactIntervalType(0.5_x, 0.5_x), 0u));
        ARIADNE_TESTING_EXPECT_SAME(singleton, UpperBoxType({
                                                   {0.5_x, 0.5_x}
        }));

        // The only feasible eighth is the last one. The upper scanning loop
        // must terminate immediately, retaining the feasible upper endpoint.
        ARIADNE_TESTING_NOTE("[constraint-shaving] final slice is the first feasible slice");
        UpperBoxType last_slice = ExactBoxType({
            {0.0_x, 1.0_x}
        });
        ARIADNE_TESTING_PROVE(
            not solver.box_reduce(last_slice, x[0], ExactIntervalType(0.9375_x, 1.0_x), 0u));
        ARIADNE_TESTING_EXPECT_SAME(last_slice, UpperBoxType({
                                                    {0.875_x, 1.0_x}
        }));
    }

    Void test_pruning_well_definedness() {
        ConstraintSolver contractor;
        auto x = ValidatedScalarMultivariateFunction::coordinates(2);

        {
            ARIADNE_TESTING_NOTE("[constraint-prune] hull reduction is reductive and preserves a known solution");
            UpperBoxType before = ExactBoxType({
                {-2.0_x, 2.0_x},
                {-2.0_x, 2.0_x}
            });
            UpperBoxType after = before;
            ValidatedScalarMultivariateFunction function = x[0] + x[1];
            ExactIntervalType codomain(1.0_x, 1.0_x);
            UpperBoxType solution = ExactBoxType({
                {0.25_x, 0.25_x},
                {0.75_x, 0.75_x}
            });

            contractor.hull_reduce(after, function, codomain);
            ARIADNE_TESTING_PROVE(refines(after, before));
            ARIADNE_TESTING_PROVE(not definitely(disjoint(solution, after)));
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-prune] nonlinear hull reduction preserves a known solution");
            UpperBoxType before = ExactBoxType({
                {-2.0_x, 2.0_x},
                {-2.0_x, 2.0_x}
            });
            UpperBoxType after = before;
            ValidatedScalarMultivariateFunction function = sqr(x[0]) + x[1];
            ExactIntervalType codomain(1.0_x, 1.0_x);
            UpperBoxType solution = ExactBoxType({
                {0.0_x, 0.0_x},
                {1.0_x, 1.0_x}
            });

            contractor.hull_reduce(after, function, codomain);
            ARIADNE_TESTING_PROVE(refines(after, before));
            ARIADNE_TESTING_PROVE(not definitely(disjoint(solution, after)));
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-prune] box shaving is reductive and preserves a known solution");
            UpperBoxType before = ExactBoxType({
                {-2.0_x, 2.0_x},
                {-2.0_x, 2.0_x}
            });
            UpperBoxType after = before;
            ValidatedScalarMultivariateFunction function = sqr(x[0]) + x[1];
            ExactIntervalType codomain(1.0_x, 1.0_x);
            UpperBoxType solution = ExactBoxType({
                {0.0_x, 0.0_x},
                {1.0_x, 1.0_x}
            });

            contractor.box_reduce(after, function, codomain, 0u);
            contractor.box_reduce(after, function, codomain, 1u);
            ARIADNE_TESTING_PROVE(refines(after, before));
            ARIADNE_TESTING_PROVE(not definitely(disjoint(solution, after)));
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-prune] definitely infeasible hull reduction closes the box");
            auto y = ValidatedScalarMultivariateFunction::coordinates(1);
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 1.0_x}
            });
            contractor.hull_reduce(
                domain, y[0], ExactIntervalType(2.0_x, 3.0_x));
            ARIADNE_TESTING_PROVE(definitely(domain.is_empty()));
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-prune] definitely infeasible box shaving closes the box");
            auto y = ValidatedScalarMultivariateFunction::coordinates(1);
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 1.0_x}
            });
            contractor.box_reduce(
                domain, y[0], ExactIntervalType(2.0_x, 3.0_x), 0u);
            ARIADNE_TESTING_PROVE(definitely(domain.is_empty()));
        }
    }

    Void test_composite_pruning_witness_preservation() {
        ConstraintSolver contractor;
        auto x = ValidatedScalarMultivariateFunction::coordinates(2);

        auto check = [&](const String& label,
                         const UpperBoxType& domain,
                         const ValidatedScalarMultivariateFunction& function,
                         const ExactIntervalType& codomain,
                         const UpperBoxType& witness) {
            ARIADNE_TESTING_NOTE("[constraint-prune-composite] " << label);
            UpperBoxType reduced = domain;
            contractor.hull_reduce(reduced, function, codomain);
            ARIADNE_TESTING_PROVE(refines(reduced, domain));
            ARIADNE_TESTING_PROVE(not definitely(disjoint(witness, reduced)));
        };

        check(
            "negative square branch survives",
            ExactBoxType({
                {-2.0_x, 2.0_x},
                { 0.0_x, 0.0_x}
        }),
            sqr(x[0]) + x[1],
            ExactIntervalType(1.0_x, 1.0_x),
            ExactBoxType({{-1.0_x, -1.0_x}, {0.0_x, 0.0_x}}));

        check(
            "abs of repeated nonlinear expression survives",
            ExactBoxType({
                {-2.0_x, 2.0_x},
                {-1.0_x, 1.0_x}
        }),
            abs(sqr(x[0]) - 1) + sqr(x[1]),
            ExactIntervalType(0.0_x, 0.0_x),
            ExactBoxType({{-1.0_x, -1.0_x}, {0.0_x, 0.0_x}}));

        check(
            "product and division composition survives",
            ExactBoxType({
                {-2.0_x, 2.0_x},
                { 1.0_x, 3.0_x}
        }),
            (x[0] * x[1]) / (x[0] + 3),
            ExactIntervalType(-1.0_x, -1.0_x),
            ExactBoxType({{-1.0_x, -1.0_x}, {2.0_x, 2.0_x}}));

        check(
            "nested elementary composition survives",
            ExactBoxType({
                {-1.0_x, 1.0_x},
                {-1.0_x, 1.0_x}
        }),
            exp(sqr(x[0])) + abs(x[1]),
            ExactIntervalType(1.0_x, 1.0_x),
            ExactBoxType({{0.0_x, 0.0_x}, {0.0_x, 0.0_x}}));

        check(
            "repeated coordinate cancellation does not lose witness",
            ExactBoxType({
                {-2.0_x, 2.0_x},
                {-2.0_x, 2.0_x}
        }),
            (x[0] - x[0]) + x[1],
            ExactIntervalType(0.5_x, 0.5_x),
            ExactBoxType({{-1.5_x, -1.5_x}, {0.5_x, 0.5_x}}));
    }

    Void test_partial_domain_pruning_witness_preservation() {
        ConstraintSolver contractor;
        auto x = ValidatedScalarMultivariateFunction::coordinates(2);

        auto check = [&](const String& label,
                         const UpperBoxType& domain,
                         const ValidatedScalarMultivariateFunction& function,
                         const ExactIntervalType& codomain,
                         const UpperBoxType& witness) {
            ARIADNE_TESTING_NOTE("[constraint-prune-domain] " << label);
            UpperBoxType reduced = domain;
            contractor.hull_reduce(reduced, function, codomain);
            ARIADNE_TESTING_PROVE(refines(reduced, domain));
            ARIADNE_TESTING_PROVE(not definitely(disjoint(witness, reduced)));
        };

        check(
            "reciprocal domain straddles zero",
            ExactBoxType({
                {-2.0_x, 2.0_x},
                { 0.0_x, 0.0_x}
        }),
            rec(x[0]) + x[1],
            ExactIntervalType(1.0_x, 1.0_x),
            ExactBoxType({{1.0_x, 1.0_x}, {0.0_x, 0.0_x}}));

        check(
            "division denominator straddles zero",
            ExactBoxType({
                { 1.0_x, 2.0_x},
                {-2.0_x, 2.0_x}
        }),
            x[0] / x[1],
            ExactIntervalType(1.0_x, 1.0_x),
            ExactBoxType({{1.0_x, 1.0_x}, {1.0_x, 1.0_x}}));

        check(
            "sqrt argument crosses negative side",
            ExactBoxType({
                {-1.0_x, 4.0_x},
                { 0.0_x, 0.0_x}
        }),
            sqrt(x[0]) + x[1],
            ExactIntervalType(1.0_x, 1.0_x),
            ExactBoxType({{1.0_x, 1.0_x}, {0.0_x, 0.0_x}}));

        check(
            "log argument crosses zero",
            ExactBoxType({
                {-1.0_x, 3.0_x},
                { 0.0_x, 0.0_x}
        }),
            log(x[0]) + x[1],
            ExactIntervalType(0.0_x, 0.0_x),
            ExactBoxType({{1.0_x, 1.0_x}, {0.0_x, 0.0_x}}));

        check(
            "asin argument crosses both boundaries",
            ExactBoxType({
                {-2.0_x, 2.0_x},
                { 0.0_x, 0.0_x}
        }),
            asin(x[0]) + x[1],
            ExactIntervalType(0.0_x, 0.0_x),
            ExactBoxType({{0.0_x, 0.0_x}, {0.0_x, 0.0_x}}));

        check(
            "acos argument crosses both boundaries",
            ExactBoxType({
                {-2.0_x, 2.0_x},
                {-2.0_x, 2.0_x}
        }),
            acos(x[0]) + x[1],
            ExactIntervalType(0.0_x, 0.0_x),
            ExactBoxType({{1.0_x, 1.0_x}, {-1.5707963267948966_pr, -1.5707963267948966_pr}}));
    }

    Void test_monotone_reduce() {
        ConstraintSolver propagator;
        auto x = ValidatedScalarMultivariateFunction::coordinates(1);

        UpperBoxType linear_domain = ExactBoxType{
            {0.0_x, 2.0_x}
        };
        UpperBoxType linear_original = linear_domain;
        Bool linear_empty = propagator.monotone_reduce(
            linear_domain,
            2 * x[0] + 1,
            ExactIntervalType(3.0_x, 3.0_x),
            0u);
        ARIADNE_TESTING_PROVE(not linear_empty);
        ARIADNE_TESTING_PROVE(refines(linear_domain, linear_original));
        ARIADNE_TESTING_PROVE(
            possibly(contains(linear_domain[0], ExactDouble(1.0_x))));

        UpperBoxType smooth_domain = ExactBoxType{
            {0.0_x, 2.0_x}
        };
        UpperBoxType smooth_original = smooth_domain;
        auto smooth_function = exp(x[0]) + x[0];
        Bool smooth_empty = propagator.monotone_reduce(
            smooth_domain,
            smooth_function,
            smooth_function.derivative(0u),
            ExactIntervalType(3.0_x, 3.0_x),
            0u);
        ARIADNE_TESTING_PROVE(not smooth_empty);
        ARIADNE_TESTING_PROVE(refines(smooth_domain, smooth_original));
        ARIADNE_TESTING_PROVE(
            smooth_domain[0].width().raw() <= smooth_original[0].width().raw());

        UpperBoxType singleton_domain = ExactBoxType{
            {1.0_x, 1.0_x}
        };
        Bool singleton_empty = propagator.monotone_reduce(
            singleton_domain,
            x[0],
            ExactIntervalType(1.0_x, 1.0_x),
            0u);
        ARIADNE_TESTING_PROVE(not singleton_empty);
        ARIADNE_TESTING_EXPECT_SAME(
            singleton_domain, UpperBoxType(ExactBoxType{
                                  {1.0_x, 1.0_x}
        }));

        UpperBoxType lower_boundary_domain = ExactBoxType{
            {0.0_x, 2.0_x}
        };
        Bool lower_boundary_empty = propagator.monotone_reduce(
            lower_boundary_domain, x[0], ExactIntervalType(0.0_x, 0.0_x), 0u);
        ARIADNE_TESTING_PROVE(not lower_boundary_empty);
        ARIADNE_TESTING_PROVE(
            possibly(contains(lower_boundary_domain[0], ExactDouble(0.0_x))));

        UpperBoxType upper_boundary_domain = ExactBoxType{
            {0.0_x, 2.0_x}
        };
        Bool upper_boundary_empty = propagator.monotone_reduce(
            upper_boundary_domain, x[0], ExactIntervalType(2.0_x, 2.0_x), 0u);
        ARIADNE_TESTING_PROVE(not upper_boundary_empty);
        ARIADNE_TESTING_PROVE(
            possibly(contains(upper_boundary_domain[0], ExactDouble(2.0_x))));

        UpperBoxType constraint_domain = ExactBoxType{
            {0.0_x, 2.0_x}
        };
        ValidatedConstraint constraint(
            ValidatedNumber(1.0_x),
            x[0],
            ValidatedNumber(1.0_x));
        Bool constraint_empty = propagator.monotone_reduce(
            constraint_domain, constraint, 0u);
        ARIADNE_TESTING_PROVE(not constraint_empty);
        ARIADNE_TESTING_PROVE(
            possibly(contains(constraint_domain[0], ExactDouble(1.0_x))));
    }

    Void test_lyapunov_reduce() {
        ConstraintSolver solver;
        Sweeper<FloatDP> sweeper;

        {
            ARIADNE_TESTING_NOTE("[constraint-lyapunov] exact centre contracts identity to equality");
            auto x = ValidatedScalarMultivariateFunction::coordinates(1);
            ExactBoxType exact_domain({
                {0.0_x, 2.0_x}
            });
            ValidatedVectorMultivariateFunction function({x[0]});
            ValidatedVectorMultivariateTaylorFunctionModelDP model(
                exact_domain, function, sweeper);
            UpperBoxType domain = exact_domain;
            ExactBoxType codomain({
                {1.0_x, 1.0_x}
            });
            Vector<FloatDP> centre({FloatDP(1.0_x, dp)});
            Vector<FloatDP> multipliers({FloatDP(1.0_x, dp)});

            Bool empty = solver.lyapunov_reduce(
                domain, model, codomain, centre, multipliers);
            ARIADNE_TESTING_PROVE(not empty);
            ARIADNE_TESTING_PROVE(
                possibly(contains(domain[0], ExactDouble(1.0_x))));
            ARIADNE_TESTING_PROVE(
                domain[0].width().raw() <= UpperIntervalType(ExactIntervalType(0.0_x, 2.0_x)).width().raw());
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-lyapunov] approximate centre overload preserves witness");
            auto x = ValidatedScalarMultivariateFunction::coordinates(2);
            ExactBoxType exact_domain({
                {0.0_x, 2.0_x},
                {1.0_x, 3.0_x}
            });
            ValidatedVectorMultivariateFunction function({x[0] + x[1]});
            ValidatedVectorMultivariateTaylorFunctionModelDP model(
                exact_domain, function, sweeper);
            UpperBoxType domain = exact_domain;
            ExactBoxType codomain({
                {3.0_x, 3.0_x}
            });
            Vector<FloatDPApproximation> centre({FloatDPApproximation(1.0_x, dp),
                                                 FloatDPApproximation(2.0_x, dp)});
            Vector<FloatDPApproximation> multipliers({FloatDPApproximation(1.0_x, dp)});

            Bool empty = solver.lyapunov_reduce(
                domain, model, codomain, centre, multipliers);
            ARIADNE_TESTING_PROVE(not empty);
            ARIADNE_TESTING_PROVE(refines(domain, UpperBoxType(exact_domain)));
            UpperBoxType witness = ExactBoxType({
                {1.0_x, 1.0_x},
                {2.0_x, 2.0_x}
            });
            ARIADNE_TESTING_PROVE(not definitely(disjoint(domain, witness)));
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-lyapunov] impossible equality closes identity domain");
            auto x = ValidatedScalarMultivariateFunction::coordinates(1);
            ExactBoxType exact_domain({
                {0.0_x, 2.0_x}
            });
            ValidatedVectorMultivariateFunction function({x[0]});
            ValidatedVectorMultivariateTaylorFunctionModelDP model(
                exact_domain, function, sweeper);
            UpperBoxType domain = exact_domain;
            ExactBoxType codomain({
                {3.0_x, 3.0_x}
            });
            Vector<FloatDP> centre({FloatDP(1.0_x, dp)});
            Vector<FloatDP> multipliers({FloatDP(1.0_x, dp)});

            ARIADNE_TESTING_PROVE(
                solver.lyapunov_reduce(
                    domain, model, codomain, centre, multipliers));
            ARIADNE_TESTING_PROVE(definitely(domain.is_empty()));
        }
    }

    Void test_propagate() {
        ConstraintSolver solver;
        auto x = ValidatedScalarMultivariateFunction::coordinates(1);

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate] hull then shaving reaches fixed point");
            UpperBoxType domain = ExactBoxType({
                {-2.0_x, 2.0_x}
            });
            List<ValidatedConstraint> constraints({ValidatedConstraint(
                ValidatedNumber(0.75_x),
                sqr(x[0]) + x[0],
                ValidatedNumber(0.75_x))});
            ConstraintPropagationStatistics statistics;
            Bool empty = solver.propagate(domain, constraints, statistics);
            ARIADNE_TESTING_PROVE(not empty);
            ARIADNE_TESTING_PROVE(statistics.hull_rounds >= 1u);
            ARIADNE_TESTING_PROVE(statistics.hull_procedure_builds >= 1u);
            ARIADNE_TESTING_PROVE(statistics.hull_procedure_build_seconds >= 0.0);
            ARIADNE_TESTING_PROVE(statistics.hull_contraction_seconds >= 0.0);
            ARIADNE_TESTING_PROVE(statistics.hull_direct_rejection_seconds >= 0.0);
            ARIADNE_TESTING_PROVE(statistics.shaving_rounds >= 1u);
            ARIADNE_TESTING_PROVE(statistics.hull_effective <= statistics.hull_rounds);
            ARIADNE_TESTING_PROVE(statistics.shaving_effective <= statistics.shaving_rounds);
            ARIADNE_TESTING_PROVE(not definitely(domain.is_empty()));
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate] disabling shaving stops after hull stall");
            UpperBoxType domain = ExactBoxType({
                {-2.0_x, 2.0_x}
            });
            List<ValidatedConstraint> constraints({ValidatedConstraint(
                ValidatedNumber(0.75_x),
                sqr(x[0]) + x[0],
                ValidatedNumber(0.75_x))});
            ConstraintPropagationStatistics statistics;
            Bool empty = solver.propagate(domain, constraints, statistics, false);
            ARIADNE_TESTING_PROVE(not empty);
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.shaving_rounds, 0u);
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.shaving_function_evaluations, 0u);
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate] disabling hull keeps direct rejection");
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 1.0_x}
            });
            List<ValidatedConstraint> constraints({ValidatedConstraint(
                ValidatedNumber(2.0_x),
                sqr(x[0]) - x[0],
                ValidatedNumber(2.0_x))});
            ConstraintPropagationStatistics statistics;
            Bool empty = solver.propagate(domain, constraints, statistics, false, false);
            ARIADNE_TESTING_PROVE(empty);
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.hull_procedure_builds, 0u);
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.shaving_rounds, 0u);
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate] direct validated rejection proves empty");
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 1.0_x}
            });
            List<ValidatedConstraint> constraints({ValidatedConstraint(
                ValidatedNumber(2.0_x),
                sqr(x[0]) - x[0],
                ValidatedNumber(2.0_x))});
            ConstraintPropagationStatistics statistics;
            Bool empty = solver.propagate(domain, constraints, statistics);
            ARIADNE_TESTING_PROVE(empty);
            ARIADNE_TESTING_PROVE(statistics.hull_rounds >= 1u);
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate] shaving proves dependency-hidden infeasibility");
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 1.0_x}
            });
            auto function = x[0] * (1 - x[0]);
            List<ValidatedConstraint> constraints({ValidatedConstraint(
                ValidatedNumber(0.32_pr),
                function,
                ValidatedNumber(0.32_pr))});
            ConstraintPropagationStatistics statistics;
            Bool empty = solver.propagate(
                domain, constraints, statistics, true, false);
            ARIADNE_TESTING_PROVE(empty);
            ARIADNE_TESTING_PROVE(definitely(domain.is_empty()));
            ARIADNE_TESTING_PROVE(statistics.shaving_function_evaluations >= 1u);
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate] shaving reduction records a nonempty effective contraction");
            auto y = ValidatedScalarMultivariateFunction::coordinates(2);
            UpperBoxType domain = ExactBoxType({
                {-1.0_x, 1.0_x},
                {-1.0_x, 1.0_x}
            });
            List<ValidatedConstraint> constraints({ValidatedConstraint(
                ValidatedNumber(0.0_x),
                sqr(y[0]),
                ValidatedNumber(0.25_x))});
            ConstraintPropagationStatistics statistics;
            Bool empty = solver.propagate(
                domain, constraints, statistics, true, false);
            ARIADNE_TESTING_PROVE(not empty);
            ARIADNE_TESTING_PROVE(statistics.shaving_coordinate_effective >= 1u);
            ARIADNE_TESTING_PROVE(statistics.shaving_effective >= 1u);
            ARIADNE_TESTING_PROVE(not definitely(domain.is_empty()));
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate] later hull contraction exposes earlier infeasibility");
            auto y = ValidatedScalarMultivariateFunction::coordinates(1);
            auto dependency_hidden = y[0] * (1 - y[0]);
            List<ValidatedConstraint> constraints({ValidatedConstraint(
                                                       ValidatedNumber(0.32_pr),
                                                       dependency_hidden,
                                                       ValidatedNumber(0.32_pr)),
                                                   ValidatedConstraint(
                                                       ValidatedNumber(0.5_x),
                                                       y[0],
                                                       ValidatedNumber(0.5_x))});
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 1.0_x}
            });
            ConstraintPropagationStatistics statistics;
            Bool empty = solver.propagate(
                domain, constraints, statistics, false, true);
            ARIADNE_TESTING_PROVE(empty);
            ARIADNE_TESTING_PROVE(
                domain[0].lower_bound().raw() == domain[0].upper_bound().raw());
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate] shaving can contract only the upper endpoint");
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 1.0_x}
            });
            List<ValidatedConstraint> constraints({ValidatedConstraint(
                ValidatedNumber(0.0_x),
                x[0],
                ValidatedNumber(0.5_x))});
            ConstraintPropagationStatistics statistics;
            Bool empty = solver.propagate(
                domain, constraints, statistics, true, false);
            ARIADNE_TESTING_PROVE(not empty);
            ARIADNE_TESTING_PROVE(
                domain[0].lower_bound().raw() == UpperIntervalType(ExactIntervalType(0.0_x, 1.0_x)).lower_bound().raw());
            ARIADNE_TESTING_PROVE(
                domain[0].upper_bound().raw() < UpperIntervalType(ExactIntervalType(0.0_x, 1.0_x)).upper_bound().raw());
            ARIADNE_TESTING_PROVE(statistics.shaving_coordinate_effective >= 1u);
        }
    }

    Void test_precompiled_propagate() {
        ConstraintSolver solver;
        auto x = ValidatedScalarMultivariateFunction::coordinates(1);

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] disabling hull keeps strict direct rejection");
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 0.0_x}
            });
            std::vector<ConstraintPropagationConstraint> constraints({
                {x[0], ExactIntervalType(0.0_x, +infty), {}, true, false}
            });
            ConstraintPropagationStatistics statistics;
            ARIADNE_TESTING_PROVE(
                solver.propagate(domain, constraints, false, statistics, false, false));
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.hull_procedure_builds, 0u);
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.shaving_rounds, 0u);
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] strict lower endpoint rejects equality");
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 0.0_x}
            });
            std::vector<ConstraintPropagationConstraint> constraints({
                {x[0], ExactIntervalType(0.0_x, +infty), {}, true, false}
            });
            ConstraintPropagationStatistics statistics;
            ARIADNE_TESTING_PROVE(solver.propagate(domain, constraints, false, statistics));
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.hull_procedure_builds, 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.monotone_rounds, 0u);
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] strict upper endpoint rejects equality");
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 0.0_x}
            });
            std::vector<ConstraintPropagationConstraint> constraints({
                {x[0], ExactIntervalType(-infty, 0.0_x), {}, false, true}
            });
            ConstraintPropagationStatistics statistics;
            ARIADNE_TESTING_PROVE(solver.propagate(domain, constraints, false, statistics));
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] strict upper feasible side survives");
            UpperBoxType domain = ExactBoxType({
                {-1.0_x, -0.5_x}
            });
            std::vector<ConstraintPropagationConstraint> constraints({
                {x[0], ExactIntervalType(-infty, 0.0_x), {}, false, true}
            });
            ConstraintPropagationStatistics statistics;
            ARIADNE_TESTING_PROVE(not solver.propagate(domain, constraints, false, statistics));
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] disabling shaving reaches monotone phase");
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 2.0_x}
            });
            auto function = exp(x[0]) + x[0];
            std::vector<ConstraintPropagationConstraint> constraints({
                {function,
                 ExactIntervalType(3.0_x, 3.0_x),
                 {function.derivative(0u)},
                 false,
                 false}
            });
            ConstraintPropagationStatistics statistics;
            Bool empty = solver.propagate(domain, constraints, true, statistics, false);
            ARIADNE_TESTING_PROVE(not empty);
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.shaving_rounds, 0u);
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.shaving_function_evaluations, 0u);
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.monotone_rounds, 1u);
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] monotone derivative contracts after stall");
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 2.0_x}
            });
            auto function = exp(x[0]) + x[0];
            std::vector<ConstraintPropagationConstraint> constraints({
                {function,
                 ExactIntervalType(3.0_x, 3.0_x),
                 {function.derivative(0u)},
                 false,
                 false}
            });
            ConstraintPropagationStatistics statistics;
            Bool empty = solver.propagate(domain, constraints, true, statistics);
            ARIADNE_TESTING_PROVE(not empty);
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.monotone_rounds, 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.monotone_effective, 1u);
            ARIADNE_TESTING_PROVE(refines(
                domain, UpperBoxType(ExactBoxType{
                            {0.0_x, 2.0_x}
            })));
            UpperIntervalType image = Ariadne::apply(function, domain);
            ARIADNE_TESTING_PROVE(
                possibly(intersect(image, ExactIntervalType(3.0_x, 3.0_x))));
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] absent derivative vector is skipped");
            UpperBoxType domain = ExactBoxType({
                {-1.0_x, 1.0_x}
            });
            std::vector<ConstraintPropagationConstraint> constraints({
                {sqr(x[0]), ExactIntervalType(0.0_x, 1.0_x), {}, false, false}
            });
            ConstraintPropagationStatistics statistics;
            ARIADNE_TESTING_PROVE(not solver.propagate(domain, constraints, true, statistics));
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.monotone_rounds, 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.monotone_effective, 0u);
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] negative monotone derivative contracts");
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 2.0_x}
            });
            auto function = -exp(x[0]) - x[0];
            std::vector<ConstraintPropagationConstraint> constraints({
                {function,
                 ExactIntervalType(-3.0_x, -3.0_x),
                 {function.derivative(0u)},
                 false,
                 false}
            });
            ConstraintPropagationStatistics statistics;
            ARIADNE_TESTING_PROVE(not solver.propagate(domain, constraints, true, statistics));
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.monotone_rounds, 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.monotone_effective, 1u);
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] missing derivative is skipped");
            UpperBoxType domain = ExactBoxType({
                {-1.0_x, 1.0_x}
            });
            std::vector<ConstraintPropagationConstraint> constraints({
                {sqr(x[0]),
                 ExactIntervalType(0.0_x, 1.0_x),
                 {std::nullopt},
                 false,
                 false}
            });
            ConstraintPropagationStatistics statistics;
            ARIADNE_TESTING_PROVE(not solver.propagate(domain, constraints, true, statistics));
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.monotone_rounds, 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.monotone_effective, 0u);
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] nonmonotone derivative is skipped");
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 4.0_x}
            });
            auto function = sin(x[0]);
            std::vector<ConstraintPropagationConstraint> constraints({
                {function,
                 ExactIntervalType(0.0_x, 0.0_x),
                 {function.derivative(0u)},
                 false,
                 false}
            });
            ConstraintPropagationStatistics statistics;
            ARIADNE_TESTING_PROVE(not solver.propagate(domain, constraints, true, statistics));
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.monotone_rounds, 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.monotone_effective, 0u);
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] cached agenda reaches sparse-chain fixpoint");
            auto y = ValidatedScalarMultivariateFunction::coordinates(2);
            auto link = y[1] - y[0];
            auto root = y[0];
            std::vector<ConstraintPropagationConstraint> constraints({
                {link,
                 ExactIntervalType(0.0_x, 0.0_x),
                 {},
                 false,
                 false,
                 std::make_shared<ValidatedProcedure>(link)},
                {root,
                 ExactIntervalType(0.0_x, 0.0_x),
                 {},
                 false,
                 false,
                 std::make_shared<ValidatedProcedure>(root)}
            });
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 1.0_x},
                {0.0_x, 1.0_x}
            });
            ConstraintPropagationStatistics statistics;
            Bool empty = solver.propagate(
                domain, constraints, false, statistics, false, true);
            ARIADNE_TESTING_PROVE(not empty);
            ARIADNE_TESTING_EXPECT_SAME(
                domain,
                UpperBoxType(ExactBoxType({
                    {0.0_x, 0.0_x},
                    {0.0_x, 0.0_x}
            })));
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.hull_procedure_builds, 0u);
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.hull_rounds, 1u);
            ARIADNE_TESTING_PROVE(
                statistics.hull_agenda_pushes > constraints.size());
            ARIADNE_TESTING_EXPECT_EQUAL(
                statistics.hull_agenda_pops, statistics.hull_agenda_pushes);
            ARIADNE_TESTING_PROVE(statistics.hull_agenda_effective_calls >= 2u);
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.hull_effective, 1u);
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] one hull sweep uses cached procedures");
            auto y = ValidatedScalarMultivariateFunction::coordinates(2);
            auto link = y[1] - y[0];
            auto root = y[0];
            std::vector<ConstraintPropagationConstraint> constraints({
                {link,
                 ExactIntervalType(0.0_x, 0.0_x),
                 {},
                 false,
                 false,
                 std::make_shared<ValidatedProcedure>(link)},
                {root,
                 ExactIntervalType(0.0_x, 0.0_x),
                 {},
                 false,
                 false,
                 std::make_shared<ValidatedProcedure>(root)}
            });
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 1.0_x},
                {0.0_x, 1.0_x}
            });
            ConstraintPropagationStatistics statistics;
            Bool empty = solver.propagate_hull_once(
                domain, constraints, statistics);
            ARIADNE_TESTING_PROVE(not empty);
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.hull_rounds, 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.hull_procedure_builds, 0u);
            ARIADNE_TESTING_PROVE(
                domain[0].lower_bound().raw() == domain[0].upper_bound().raw());
            ARIADNE_TESTING_PROVE(
                domain[1].width().raw() == UpperIntervalType(
                                               ExactIntervalType(0.0_x, 1.0_x))
                                               .width()
                                               .raw());
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] uncached hull-only propagation falls back to full scan");
            auto y = ValidatedScalarMultivariateFunction::coordinates(1);
            std::vector<ConstraintPropagationConstraint> constraints({
                {y[0],
                 ExactIntervalType(0.0_x, 0.0_x),
                 {},
                 false,
                 false}
            });
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 1.0_x}
            });
            ConstraintPropagationStatistics statistics;
            Bool empty = solver.propagate(
                domain, constraints, false, statistics, false, true);
            ARIADNE_TESTING_PROVE(not empty);
            ARIADNE_TESTING_PROVE(statistics.hull_procedure_builds >= 1u);
            ARIADNE_TESTING_PROVE(
                domain[0].lower_bound().raw() == domain[0].upper_bound().raw());
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] cached agenda stable constraint has no effective call");
            auto y = ValidatedScalarMultivariateFunction::coordinates(1);
            auto function = y[0];
            std::vector<ConstraintPropagationConstraint> constraints({
                {function,
                 ExactIntervalType(0.0_x, 1.0_x),
                 {},
                 false,
                 false,
                 std::make_shared<ValidatedProcedure>(function)}
            });
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 1.0_x}
            });
            ConstraintPropagationStatistics statistics;
            ARIADNE_TESTING_PROVE(
                not solver.propagate(
                    domain, constraints, false, statistics, false, true));
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.hull_agenda_effective_calls, 0u);
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.hull_effective, 0u);
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] cached agenda hull contraction proves empty");
            auto y = ValidatedScalarMultivariateFunction::coordinates(1);
            auto function = y[0];
            std::vector<ConstraintPropagationConstraint> constraints({
                {function,
                 ExactIntervalType(2.0_x, 3.0_x),
                 {},
                 false,
                 false,
                 std::make_shared<ValidatedProcedure>(function)}
            });
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 1.0_x}
            });
            ConstraintPropagationStatistics statistics;
            ARIADNE_TESTING_PROVE(
                solver.propagate(
                    domain, constraints, false, statistics, false, true));
            ARIADNE_TESTING_PROVE(definitely(domain.is_empty()));
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] cached agenda applies strict endpoint rejection");
            auto y = ValidatedScalarMultivariateFunction::coordinates(1);
            auto function = y[0];
            std::vector<ConstraintPropagationConstraint> constraints({
                {function,
                 ExactIntervalType(0.0_x, +infty),
                 {},
                 true,
                 false,
                 std::make_shared<ValidatedProcedure>(function)}
            });
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 0.0_x}
            });
            ConstraintPropagationStatistics statistics;
            ARIADNE_TESTING_PROVE(
                solver.propagate(
                    domain, constraints, false, statistics, false, true));
            ARIADNE_TESTING_PROVE(not definitely(domain.is_empty()));
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] strict lower feasible side survives");
            auto y = ValidatedScalarMultivariateFunction::coordinates(1);
            UpperBoxType domain = ExactBoxType({
                {0.5_x, 1.0_x}
            });
            std::vector<ConstraintPropagationConstraint> constraints({
                {y[0], ExactIntervalType(0.0_x, +infty), {}, true, false}
            });
            ConstraintPropagationStatistics statistics;
            ARIADNE_TESTING_PROVE(
                not solver.propagate(
                    domain, constraints, false, statistics, false, false));
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] one uncached hull sweep stays stable");
            auto y = ValidatedScalarMultivariateFunction::coordinates(1);
            std::vector<ConstraintPropagationConstraint> constraints({
                {y[0], ExactIntervalType(0.0_x, 1.0_x), {}, false, false}
            });
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 1.0_x}
            });
            ConstraintPropagationStatistics statistics;
            ARIADNE_TESTING_PROVE(
                not solver.propagate_hull_once(domain, constraints, statistics));
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.hull_procedure_builds, 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.hull_effective, 0u);
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] one hull sweep can prove empty");
            auto y = ValidatedScalarMultivariateFunction::coordinates(1);
            std::vector<ConstraintPropagationConstraint> constraints({
                {y[0], ExactIntervalType(2.0_x, 3.0_x), {}, false, false}
            });
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 1.0_x}
            });
            ConstraintPropagationStatistics statistics;
            ARIADNE_TESTING_PROVE(
                solver.propagate_hull_once(domain, constraints, statistics));
            ARIADNE_TESTING_PROVE(definitely(domain.is_empty()));
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] one cached hull sweep applies strict rejection");
            auto y = ValidatedScalarMultivariateFunction::coordinates(1);
            auto function = y[0];
            std::vector<ConstraintPropagationConstraint> constraints({
                {function,
                 ExactIntervalType(0.0_x, +infty),
                 {},
                 true,
                 false,
                 std::make_shared<ValidatedProcedure>(function)}
            });
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 0.0_x}
            });
            ConstraintPropagationStatistics statistics;
            ARIADNE_TESTING_PROVE(
                solver.propagate_hull_once(domain, constraints, statistics));
            ARIADNE_TESTING_PROVE(not definitely(domain.is_empty()));
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] cached full scan reuses procedure");
            auto y = ValidatedScalarMultivariateFunction::coordinates(1);
            auto function = y[0];
            std::vector<ConstraintPropagationConstraint> constraints({
                {function,
                 ExactIntervalType(0.0_x, 1.0_x),
                 {},
                 false,
                 false,
                 std::make_shared<ValidatedProcedure>(function)}
            });
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 1.0_x}
            });
            ConstraintPropagationStatistics statistics;
            ARIADNE_TESTING_PROVE(
                not solver.propagate(
                    domain, constraints, true, statistics, false, true));
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.hull_procedure_builds, 0u);
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.monotone_rounds, 1u);
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] full scan hull contraction proves empty");
            auto y = ValidatedScalarMultivariateFunction::coordinates(1);
            std::vector<ConstraintPropagationConstraint> constraints({
                {y[0], ExactIntervalType(2.0_x, 3.0_x), {}, false, false}
            });
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 1.0_x}
            });
            ConstraintPropagationStatistics statistics;
            ARIADNE_TESTING_PROVE(
                solver.propagate(
                    domain, constraints, true, statistics, false, true));
            ARIADNE_TESTING_PROVE(definitely(domain.is_empty()));
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] shaving can prove dependency-hidden infeasibility");
            auto y = ValidatedScalarMultivariateFunction::coordinates(1);
            auto function = y[0] * (1 - y[0]);
            std::vector<ConstraintPropagationConstraint> constraints({
                {function,
                 ExactIntervalType(0.32_pr, 0.32_pr),
                 {},
                 false,
                 false}
            });
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 1.0_x}
            });
            ConstraintPropagationStatistics statistics;
            ARIADNE_TESTING_PROVE(
                solver.propagate(
                    domain, constraints, false, statistics, true, false));
            ARIADNE_TESTING_PROVE(definitely(domain.is_empty()));
            ARIADNE_TESTING_PROVE(statistics.shaving_function_evaluations >= 1u);
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] one hull sweep rechecks strict constraints after contraction");
            auto y = ValidatedScalarMultivariateFunction::coordinates(1);
            std::vector<ConstraintPropagationConstraint> constraints({
                {y[0], ExactIntervalType(0.0_x, +infty), {},  true, false},
                {y[0], ExactIntervalType(0.0_x,  0.0_x), {}, false, false}
            });
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 1.0_x}
            });
            ConstraintPropagationStatistics statistics;
            ARIADNE_TESTING_PROVE(
                solver.propagate_hull_once(domain, constraints, statistics));
            ARIADNE_TESTING_PROVE(not definitely(domain.is_empty()));
            ARIADNE_TESTING_PROVE(
                domain[0].lower_bound().raw() == domain[0].upper_bound().raw());
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] full scan rechecks strict constraints after contraction");
            auto y = ValidatedScalarMultivariateFunction::coordinates(1);
            std::vector<ConstraintPropagationConstraint> constraints({
                {y[0], ExactIntervalType(0.0_x, +infty), {},  true, false},
                {y[0], ExactIntervalType(0.0_x,  0.0_x), {}, false, false}
            });
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 1.0_x}
            });
            ConstraintPropagationStatistics statistics;
            ARIADNE_TESTING_PROVE(
                solver.propagate(
                    domain, constraints, true, statistics, false, true));
            ARIADNE_TESTING_PROVE(not definitely(domain.is_empty()));
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] adaptive shaving skips inactive coordinates");
            auto y = ValidatedScalarMultivariateFunction::coordinates(2);
            std::vector<ConstraintPropagationConstraint> constraints({
                {sqr(y[0]),
                 ExactIntervalType(0.0_x, 0.25_x),
                 {},
                 false,
                 false}
            });
            UpperBoxType domain = ExactBoxType({
                {-1.0_x, 1.0_x},
                {-1.0_x, 1.0_x}
            });
            ConstraintPropagationStatistics statistics;
            ARIADNE_TESTING_PROVE(
                not solver.propagate(
                    domain, constraints, false, statistics, true, false));
            ARIADNE_TESTING_PROVE(statistics.shaving_effective >= 1u);
            ARIADNE_TESTING_PROVE(statistics.shaving_active_rounds >= 1u);
            ARIADNE_TESTING_PROVE(statistics.shaving_adaptive_skipped >= 1u);
            ARIADNE_TESTING_PROVE(not definitely(domain.is_empty()));
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] cached agenda records lower-only coordinate change");
            auto y = ValidatedScalarMultivariateFunction::coordinates(1);
            auto function = y[0];
            std::vector<ConstraintPropagationConstraint> constraints({
                {function,
                 ExactIntervalType(0.5_x, 1.0_x),
                 {},
                 false,
                 false,
                 std::make_shared<ValidatedProcedure>(function)}
            });
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 1.0_x}
            });
            ConstraintPropagationStatistics statistics;
            ARIADNE_TESTING_PROVE(
                not solver.propagate(
                    domain, constraints, false, statistics, false, true));
            ARIADNE_TESTING_PROVE(
                domain[0].lower_bound().raw() > UpperIntervalType(ExactIntervalType(0.0_x, 1.0_x)).lower_bound().raw());
            ARIADNE_TESTING_PROVE(
                domain[0].upper_bound().raw() == UpperIntervalType(ExactIntervalType(0.0_x, 1.0_x)).upper_bound().raw());
            ARIADNE_TESTING_EXPECT_EQUAL(statistics.hull_effective, 1u);
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-propagate-precompiled] cached contraction coalesces duplicate wakeups");
            auto y = ValidatedScalarMultivariateFunction::coordinates(2);
            auto function = y[0] + y[1];
            std::vector<ConstraintPropagationConstraint> constraints({
                {function,
                 ExactIntervalType(0.0_x, 0.0_x),
                 {},
                 false,
                 false,
                 std::make_shared<ValidatedProcedure>(function)}
            });
            UpperBoxType domain = ExactBoxType({
                {0.0_x, 1.0_x},
                {0.0_x, 1.0_x}
            });
            ConstraintPropagationStatistics statistics;
            ARIADNE_TESTING_PROVE(
                not solver.propagate(
                    domain, constraints, false, statistics, false, true));
            ARIADNE_TESTING_EXPECT_SAME(
                domain,
                UpperBoxType(ExactBoxType({
                    {0.0_x, 0.0_x},
                    {0.0_x, 0.0_x}
            })));
            ARIADNE_TESTING_PROVE(statistics.hull_agenda_effective_calls >= 1u);
            ARIADNE_TESTING_PROVE(statistics.hull_agenda_pushes >= 2u);
        }
    }

    Void test_reduce_certified_containment() {
        ConstraintSolver solver;
        auto x = ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction identity({x[0]});

        // Satisfied outputs cannot remove domain points.
        UpperBoxType contained = ExactBoxType({
            {0.0_x, 1.0_x}
        });
        UpperBoxType original = contained;
        ARIADNE_TESTING_PROVE(
            not solver.reduce(contained, identity, ExactBoxType({
                                                       {-1.0_x, 2.0_x}
        })));
        ARIADNE_TESTING_EXPECT_SAME(contained, original);

        // Equality at the codomain endpoints is still containment.
        UpperBoxType touching = original;
        ARIADNE_TESTING_PROVE(
            not solver.reduce(touching, identity, ExactBoxType({
                                                      {0.0_x, 1.0_x}
        })));
        ARIADNE_TESTING_EXPECT_SAME(touching, original);

        // The original 0.75 stopping rule and a subsequent reduction round
        // must still produce a sound contraction.
        UpperBoxType narrowed = original;
        ARIADNE_TESTING_PROVE(
            not solver.reduce(narrowed, identity, ExactBoxType({
                                                      {0.5_x, 0.5625_x}
        })));
        ARIADNE_TESTING_PROVE(refines(narrowed, original));
        ARIADNE_TESTING_PROVE(
            refines(narrowed, UpperBoxType(ExactBoxType({
                                  {0.5_x, 0.5625_x}
        }))));
        ARIADNE_TESTING_PROVE(not definitely(disjoint(
            narrowed, UpperBoxType(ExactBoxType({
                          {0.5_x, 0.5_x}
        })))));

        // Neither early-certified infeasibility nor empty/singleton domains
        // should be affected by the optimization.
        UpperBoxType disjoint_domain = original;
        ARIADNE_TESTING_PROVE(solver.reduce(
            disjoint_domain, identity, ExactBoxType({
                                           {2.0_x, 3.0_x}
        })));
        ARIADNE_TESTING_PROVE(definitely(disjoint_domain.is_empty()));

        UpperBoxType singleton = ExactBoxType({
            {0.5_x, 0.5_x}
        });
        UpperBoxType singleton_before = singleton;
        ARIADNE_TESTING_PROVE(not solver.reduce(
            singleton, identity, ExactBoxType({
                                     {0.5_x, 0.5_x}
        })));
        ARIADNE_TESTING_EXPECT_SAME(singleton, singleton_before);

        UpperBoxType empty = ExactBoxType({ExactIntervalType::empty_interval()});
        ARIADNE_TESTING_PROVE(solver.reduce(
            empty, identity, ExactBoxType({
                                 {0.0_x, 1.0_x}
        })));
        ARIADNE_TESTING_PROVE(definitely(empty.is_empty()));

        // Natural interval over-approximation cannot dismiss a nonlinear
        // constraint that actually requires shaving.
        ValidatedVectorMultivariateFunction nonlinear({x[0] * (1 - x[0])});
        UpperBoxType impossible = original;
        ARIADNE_TESTING_PROVE(solver.reduce(
            impossible, nonlinear, ExactBoxType({
                                       {0.32_pr, 0.32_pr}
        })));
        ARIADNE_TESTING_PROVE(definitely(impossible.is_empty()));

        // A certified-satisfied output coexists with a coupled nonlinear
        // constraint. Preserve a mathematically valid witness.
        auto xy = ValidatedScalarMultivariateFunction::coordinates(2);
        ValidatedVectorMultivariateFunction mixed({xy[0], sqr(xy[0]) + xy[1]});
        UpperBoxType mixed_original = ExactBoxType({
            {-1.0_x, 1.0_x},
            { 0.0_x, 1.0_x}
        });
        UpperBoxType mixed_domain = mixed_original;
        ARIADNE_TESTING_PROVE(not solver.reduce(
            mixed_domain, mixed,
            ExactBoxType({
                {-1.0_x,  1.0_x},
                { 0.5_x, 0.75_x}
        })));
        ARIADNE_TESTING_PROVE(refines(mixed_domain, mixed_original));
        ARIADNE_TESTING_PROVE(not definitely(disjoint(
            mixed_domain,
            UpperBoxType(ExactBoxType({
                {  0.0_x,   0.0_x},
                {0.625_x, 0.625_x}
        })))));

        // Independent deterministic 2D inputs (not the paving traces).
        ExactIntervalType slices[] = {
            ExactIntervalType(0.0_x, 0.25_x),
            ExactIntervalType(0.25_x, 0.5_x),
            ExactIntervalType(0.5_x, 0.75_x),
            ExactIntervalType(0.75_x, 1.0_x)};
        ValidatedVectorMultivariateFunction coupled({xy[0] + 2 * xy[1], xy[0] * xy[1]});
        for (SizeType i = 0u; i != 4u; ++i) {
            for (SizeType j = 0u; j != 4u; ++j) {
                ARIADNE_TESTING_NOTE("[constraint-contained] 2D grid " << i << "," << j);
                UpperBoxType before = ExactBoxType({slices[i], slices[j]});
                UpperBoxType after = before;
                const Bool proven_empty = solver.reduce(
                    after, coupled,
                    ExactBoxType({
                        {0.5_x,  1.5_x},
                        {0.0_x, 0.25_x}
                }));
                ARIADNE_TESTING_PROVE(refines(after, before));
                ARIADNE_TESTING_PROVE(proven_empty == definitely(after.is_empty()));
                if (i == 3u and j == 3u) {
                    ARIADNE_TESTING_PROVE(proven_empty);
                }
                if (i == 0u and j == 0u) {
                    ARIADNE_TESTING_PROVE(not proven_empty);
                    ARIADNE_TESTING_PROVE(not definitely(disjoint(
                        after, UpperBoxType(ExactBoxType({
                                   {0.25_x, 0.25_x},
                                   {0.25_x, 0.25_x}
                    })))));
                }
            }
        }

        // Three-dimensional mixed polynomial constraints preserve a
        // known feasible point despite intermediate contractions.
        auto xyz = ValidatedScalarMultivariateFunction::coordinates(3);
        ValidatedVectorMultivariateFunction three_dimensional({xyz[0], xyz[0] + xyz[1] + xyz[2], sqr(xyz[1]) + xyz[2]});
        UpperBoxType initial_3d = ExactBoxType({
            {0.0_x, 1.0_x},
            {0.0_x, 1.0_x},
            {0.0_x, 1.0_x}
        });
        UpperBoxType contracted_3d = initial_3d;
        ARIADNE_TESTING_PROVE(not solver.reduce(
            contracted_3d, three_dimensional,
            ExactBoxType({
                { 0.0_x,  1.0_x},
                { 1.0_x, 1.25_x},
                {0.25_x,  1.0_x}
        })));
        ARIADNE_TESTING_PROVE(refines(contracted_3d, initial_3d));
        ARIADNE_TESTING_PROVE(not definitely(disjoint(
            contracted_3d, UpperBoxType(ExactBoxType({
                               { 0.5_x,  0.5_x},
                               {0.25_x, 0.25_x},
                               { 0.5_x,  0.5_x}
        })))));
    }

    Void test_prepared_reduce_matches_production() {
        ConstraintSolver solver;
        auto x = ValidatedScalarMultivariateFunction::coordinates(2);
        ValidatedVectorMultivariateFunction function({x[0] + x[1], x[0] * x[1], x[0] - x[1]});
        const PreparedConstraintSolverFunction prepared(function);
        const PreparedConstraintSolverFunction shared(prepared);
        const ExactBoxType initial({
            {0.0_x, 1.0_x},
            {0.0_x, 1.0_x}
        });
        std::vector<ExactBoxType> codomains({
            ExactBoxType({ {-1.0_x, 3.0_x}, {-1.0_x, 2.0_x},   {-2.0_x, 2.0_x}}
            ),
            ExactBoxType({{0.75_x, 1.25_x}, {0.0_x, 0.25_x}, {-0.25_x, 0.25_x}}
            ),
            ExactBoxType({  {3.0_x, 4.0_x},  {0.0_x, 1.0_x},   {-1.0_x, 1.0_x}}
            )
        });
        // Alternate codomains and repeat, covering reused, independently
        // copied prepared contexts and certified-empty outcomes.
        for (SizeType repetition = 0u; repetition != 3u; ++repetition) {
            for (SizeType index = 0u; index < codomains.size(); ++index) {
                UpperBoxType reference = initial;
                UpperBoxType candidate = initial;
                const Bool baseline = solver.reduce(reference, function, codomains[index]);
                const Bool prepped = solver.reduce(
                    candidate, (repetition % 2u == 0u ? prepared : shared),
                    codomains[index]);
                ARIADNE_TESTING_EXPECT_EQUAL(baseline, prepped);
                ARIADNE_TESTING_EXPECT_SAME(reference, candidate);
            }
        }

        // The prepared object owns a function value beyond the source's
        // lexical lifetime; a second function cannot contaminate it.
        auto make_prepared = []() {
            auto y = ValidatedScalarMultivariateFunction::coordinates(1);
            ValidatedVectorMultivariateFunction temporary({y[0] * y[0]});
            return PreparedConstraintSolverFunction(temporary);
        };
        const PreparedConstraintSolverFunction owned = make_prepared();
        auto y = ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction square({y[0] * y[0]});
        const ExactBoxType target({
            {0.25_x, 0.5_x}
        });
        UpperBoxType reference = ExactBoxType({
            {-1.0_x, 1.0_x}
        });
        UpperBoxType candidate = reference;
        solver.reduce(reference, square, target);
        solver.reduce(candidate, owned, target);
        ARIADNE_TESTING_EXPECT_SAME(reference, candidate);
    }

    Void test_reduce_per_query_procedure_lifetime() {
        ConstraintSolver solver;
        auto x = ValidatedScalarMultivariateFunction::coordinates(2);
        ValidatedVectorMultivariateFunction function({x[0] + x[1], x[0] * x[1]});
        const ExactBoxType initial({
            {0.0_x, 1.0_x},
            {0.0_x, 1.0_x}
        });
        const ExactBoxType broad_codomain({
            {-1.0_x, 3.0_x},
            {-1.0_x, 2.0_x}
        });
        const ExactBoxType narrow_codomain({
            {0.75_x, 1.25_x},
            { 0.0_x, 0.25_x}
        });
        const ExactBoxType impossible_codomain({
            {3.0_x, 4.0_x},
            {0.0_x, 1.0_x}
        });

        // Alternate codomains while reusing the same function and solver.
        // Neither the certified-containment state nor the procedure may
        // accidentally depend on the previous reduce() query.
        UpperBoxType broad_first = initial;
        ARIADNE_TESTING_PROVE(not solver.reduce(
            broad_first, function, broad_codomain));
        ARIADNE_TESTING_EXPECT_SAME(broad_first, UpperBoxType(initial));

        UpperBoxType narrowed = initial;
        ARIADNE_TESTING_PROVE(not solver.reduce(
            narrowed, function, narrow_codomain));
        ARIADNE_TESTING_PROVE(refines(narrowed, UpperBoxType(initial)));
        ARIADNE_TESTING_PROVE(not definitely(disjoint(
            narrowed, UpperBoxType(ExactBoxType({
                          {0.5_x, 0.5_x},
                          {0.5_x, 0.5_x}
        })))));

        UpperBoxType rejected = initial;
        ARIADNE_TESTING_PROVE(solver.reduce(
            rejected, function, impossible_codomain));
        ARIADNE_TESTING_PROVE(definitely(rejected.is_empty()));

        UpperBoxType broad_again = initial;
        ARIADNE_TESTING_PROVE(not solver.reduce(
            broad_again, function, broad_codomain));
        ARIADNE_TESTING_EXPECT_SAME(broad_again, broad_first);

        UpperBoxType narrow_again = initial;
        ARIADNE_TESTING_PROVE(not solver.reduce(
            narrow_again, function, narrow_codomain));
        ARIADNE_TESTING_EXPECT_SAME(narrow_again, narrowed);

        // A separate function of the same dimension and result count
        // must compile its own vector procedure, never the previous one.
        ValidatedVectorMultivariateFunction other({x[0], x[1]});
        UpperBoxType other_domain = initial;
        ARIADNE_TESTING_PROVE(not solver.reduce(
            other_domain, other, ExactBoxType({
                                     {0.25_x, 0.75_x},
                                     {0.25_x, 0.75_x}
        })));
        ARIADNE_TESTING_PROVE(refines(other_domain, UpperBoxType(initial)));
        ARIADNE_TESTING_PROVE(not definitely(disjoint(
            other_domain, UpperBoxType(ExactBoxType({
                              {0.5_x, 0.5_x},
                              {0.5_x, 0.5_x}
        })))));
    }

    Void test_reduce_edge_cases() {
        ConstraintSolver solver;
        auto x = ValidatedScalarMultivariateFunction::coordinates(1);

        ARIADNE_TESTING_NOTE("[constraint-reduce] already empty vector domain");
        UpperBoxType empty_domain = ExactBoxType({ExactIntervalType::empty_interval()});
        ValidatedVectorMultivariateFunction vector_function({x[0]});
        ExactBoxType vector_codomain({ExactIntervalType(0.0_x, 1.0_x)});
        ARIADNE_TESTING_PROVE(
            solver.reduce(empty_domain, vector_function, vector_codomain));

        ARIADNE_TESTING_NOTE("[constraint-reduce] already empty deprecated list domain");
        UpperBoxType empty_list_domain = ExactBoxType({ExactIntervalType::empty_interval()});
        List<ValidatedConstraint> one_constraint({ValidatedConstraint(
            ValidatedNumber(0.0_x),
            x[0],
            ValidatedNumber(1.0_x))});
        ARIADNE_TESTING_PROVE(
            solver.reduce(empty_list_domain, one_constraint));

        ARIADNE_TESTING_NOTE("[constraint-reduce] nonempty vector domain reaches fixed point");
        UpperBoxType stable_domain = ExactBoxType({
            {0.0_x, 1.0_x}
        });
        ARIADNE_TESTING_PROVE(
            not solver.reduce(stable_domain, vector_function, vector_codomain));
        ARIADNE_TESTING_PROVE(not definitely(stable_domain.is_empty()));

        // The vector reduction must stop as soon as hull propagation proves
        // that the coordinate range cannot intersect the requested codomain.
        ARIADNE_TESTING_NOTE("[constraint-reduce] vector hull reduction proves infeasibility");
        UpperBoxType disjoint_domain = ExactBoxType({
            {0.0_x, 1.0_x}
        });
        ExactBoxType disjoint_codomain({
            {2.0_x, 3.0_x}
        });
        ARIADNE_TESTING_PROVE(
            solver.reduce(disjoint_domain, vector_function, disjoint_codomain));
        ARIADNE_TESTING_PROVE(definitely(disjoint_domain.is_empty()));

        // Hull reduction contracts the coordinate to the narrow allowed
        // range. The vector reducer must repeat after this large contraction
        // rather than treating the first sweep as a fixed point.
        ARIADNE_TESTING_NOTE("[constraint-reduce] vector reduction repeats after significant contraction");
        UpperBoxType contracting_vector_domain = ExactBoxType({
            {0.0_x, 1.0_x}
        });
        ExactBoxType narrow_codomain({
            {0.5_x, 0.5625_x}
        });
        ARIADNE_TESTING_PROVE(
            not solver.reduce(contracting_vector_domain, vector_function, narrow_codomain));
        ARIADNE_TESTING_PROVE(not definitely(contracting_vector_domain.is_empty()));
        ARIADNE_TESTING_PROVE(
            contracting_vector_domain[0].width().raw() < FloatDP(0.75_x, dp));
        ARIADNE_TESTING_PROVE(
            not definitely(disjoint(
                contracting_vector_domain, UpperBoxType({
                                               {0.5_x, 0.5_x}
        }))));

        // Interval evaluation of x*(1-x) over [0,1] contains 0.32,
        // although the true maximum is 0.25. Hull propagation alone does
        // not establish infeasibility; coordinate shaving must reject it.
        ARIADNE_TESTING_NOTE("[constraint-reduce] shaving rejects dependency-hidden infeasibility");
        auto nonlinear_x = ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction nonlinear_function({nonlinear_x[0] * (1 - nonlinear_x[0])});
        UpperBoxType nonlinear_domain = ExactBoxType({
            {0.0_x, 1.0_x}
        });
        ExactBoxType impossible_codomain({
            {0.32_pr, 0.32_pr}
        });
        ARIADNE_TESTING_PROVE(
            solver.reduce(nonlinear_domain, nonlinear_function, impossible_codomain));
        ARIADNE_TESTING_PROVE(definitely(nonlinear_domain.is_empty()));

        ARIADNE_TESTING_NOTE("[constraint-reduce] nonempty deprecated list reaches fixed point");
        UpperBoxType stable_list_domain = ExactBoxType({
            {0.0_x, 1.0_x}
        });
        List<ValidatedConstraint> stable_constraints({ValidatedConstraint(
            ValidatedNumber(0.0_x),
            x[0],
            ValidatedNumber(1.0_x))});
        ARIADNE_TESTING_PROVE(
            not solver.reduce(stable_list_domain, stable_constraints));
        ARIADNE_TESTING_PROVE(not definitely(stable_list_domain.is_empty()));

        ARIADNE_TESTING_NOTE("[constraint-reduce] deprecated list repeats after strong contraction");
        UpperBoxType contracting_list_domain = ExactBoxType({
            {0.0_x, 1.0_x}
        });
        List<ValidatedConstraint> contracting_constraints({ValidatedConstraint(
            ValidatedNumber(0.0_x),
            x[0],
            ValidatedNumber(0.0_x))});
        ARIADNE_TESTING_PROVE(
            not solver.reduce(contracting_list_domain, contracting_constraints));
        ARIADNE_TESTING_PROVE(
            contracting_list_domain[0].lower_bound().raw() == contracting_list_domain[0].upper_bound().raw());

        ARIADNE_TESTING_NOTE("[constraint-reduce] reject argument/domain dimension mismatch");
        UpperBoxType wrong_domain = ExactBoxType({
            {0.0_x, 1.0_x},
            {0.0_x, 1.0_x}
        });
        ARIADNE_TESTING_EXPECT_EXCEPTION(
            solver.reduce(wrong_domain, vector_function, vector_codomain));

        ARIADNE_TESTING_NOTE("[constraint-reduce] reject result/codomain dimension mismatch");
        UpperBoxType one_dimensional_domain = ExactBoxType({
            {0.0_x, 1.0_x}
        });
        ExactBoxType empty_codomain(0u);
        ARIADNE_TESTING_EXPECT_EXCEPTION(
            solver.reduce(
                one_dimensional_domain, vector_function, empty_codomain));
    }

    Void test_check_feasibility() {
        ConstraintSolver solver;
        auto x = ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction function({x[0]});
        ExactBoxType domain({
            {0.0_x, 2.0_x}
        });
        ExactBoxType codomain({
            {0.0_x, 1.0_x}
        });

        ARIADNE_TESTING_NOTE("[constraint-feasible] reject point outside domain");
        ConstraintSolver::ExactPointType outside({FloatDP(-1.0_x, dp)});
        ARIADNE_TESTING_PROVE(
            not possibly(solver.check_feasibility(
                domain, function, codomain, outside)));

        ARIADNE_TESTING_NOTE("[constraint-feasible] reject image outside codomain");
        ConstraintSolver::ExactPointType image_outside({FloatDP(2.0_x, dp)});
        ARIADNE_TESTING_PROVE(
            not possibly(solver.check_feasibility(
                domain, function, codomain, image_outside)));

        ARIADNE_TESTING_NOTE("[constraint-feasible] boundary point is indeterminate");
        ConstraintSolver::ExactPointType boundary({FloatDP(0.0_x, dp)});
        ARIADNE_TESTING_PROVE(
            is_indeterminate(solver.check_feasibility(
                domain, function, codomain, boundary)));

        ARIADNE_TESTING_NOTE("[constraint-feasible] strict interior point is validated");
        ConstraintSolver::ExactPointType interior({FloatDP(0.5_x, dp)});
        ARIADNE_TESTING_PROVE(
            definitely(solver.check_feasibility(
                domain, function, codomain, interior)));
    }

    Void test_logging_branches() {
        ARIADNE_TESTING_NOTE("[constraint-coverage] exercise logging-enabled solver paths");
        auto& configuration = Logger::instance().configuration();
        const unsigned int old_verbosity = configuration.verbosity();
        configuration.set_verbosity(10u);

        ConstraintSolver solver;
        auto x = ValidatedScalarMultivariateFunction::coordinates(1);
        ExactIntervalType bounds(0.0_x, 1.0_x);

        UpperBoxType scalar_domain = ExactBoxType({
            {0.0_x, 1.0_x}
        });
        solver.hull_reduce(scalar_domain, x[0], bounds);

        UpperBoxType vector_domain = ExactBoxType({
            {0.0_x, 1.0_x}
        });
        ValidatedVectorMultivariateFunction vector_function({x[0]});
        ExactBoxType vector_bounds({bounds});
        solver.hull_reduce(vector_domain, vector_function, vector_bounds);

        UpperBoxType monotone_domain = ExactBoxType({
            {0.0_x, 2.0_x}
        });
        solver.monotone_reduce(
            monotone_domain, x[0], ExactIntervalType(1.0_x, 1.0_x), 0u);

        UpperBoxType shaving_domain = ExactBoxType({
            {0.0_x, 1.0_x}
        });
        solver.box_reduce(
            shaving_domain, x[0], ExactIntervalType(0.0_x, 0.5_x), 0u);

        ConstraintSolver::ExactPointType point({FloatDP(0.5_x, dp)});
        ARIADNE_TESTING_PROVE(definitely(solver.check_feasibility(
            ExactBoxType({
                {0.0_x, 1.0_x}
        }),
            vector_function,
            ExactBoxType({{0.0_x, 1.0_x}}),
            point)));

        {
            auto y = ValidatedScalarMultivariateFunction::coordinates(1);
            ValidatedVectorMultivariateFunction function({y[0]});
            ExactBoxType domain({
                {0.0_x, 1.0_x}
            });
            ARIADNE_TESTING_PROVE(not possibly(
                solver.feasible(
                          domain, function, ExactBoxType({
                                                {2.0_x, 3.0_x}
            }))
                    .first));
        }

        {
            auto y = ValidatedScalarMultivariateFunction::coordinates(1);
            ValidatedVectorMultivariateFunction function({y[0] * (1 - y[0])});
            ExactBoxType domain({
                {0.0_x, 1.0_x}
            });
            auto feasibility_result = solver.feasible(
                domain, function, ExactBoxType({
                                      {0.32_pr, 0.32_pr}
            }));
            ARIADNE_TESTING_PROVE(not definitely(feasibility_result.first));
        }

        configuration.set_verbosity(old_verbosity);
    }

    Void test_virtual_destruction() {
        ARIADNE_TESTING_NOTE("[constraint-lifetime] destroy ConstraintSolver through interface");
        ConstraintSolverInterface* solver = new ConstraintSolver();
        delete solver;
    }

    Void test_split() {
        ConstraintSolver solver;
        auto x = ValidatedScalarMultivariateFunction::coordinates(2);
        UpperBoxType domain = ExactBoxType({
            {0.0_x, 4.0_x},
            {0.0_x, 1.0_x}
        });
        ValidatedVectorMultivariateFunction function({x[0] + x[1]});
        ExactBoxType codomain({
            {1.0_x, 2.0_x}
        });

        ARIADNE_TESTING_NOTE("[constraint-split] delegates to box bisection");
        auto children = solver.split(domain, function, codomain);
        ARIADNE_TESTING_PROVE(refines(children.first, domain));
        ARIADNE_TESTING_PROVE(refines(children.second, domain));
        ARIADNE_TESTING_PROVE(not same(children.first, children.second));
        ARIADNE_TESTING_PROVE(
            children.first[0].upper_bound().raw() == children.second[0].lower_bound().raw());
    }

    Void test_feasible() {
        ConstraintSolver contractor;
        {
            ARIADNE_TESTING_NOTE("[constraint-feasible] reject empty codomain");
            auto coordinate = ValidatedScalarMultivariateFunction::coordinates(1);
            ValidatedVectorMultivariateFunction empty_function(
                0u, coordinate[0].domain());
            ARIADNE_TESTING_EXPECT_EXCEPTION(
                contractor.feasible(
                    ExactBoxType({
                        {0.0_x, 1.0_x}
            }),
                    empty_function,
                    ExactBoxType(0u)));
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-feasible] empty constraint system");
            List<ValidatedConstraint> no_constraints;
            ExactBoxType nonempty_domain({
                {0.0_x, 1.0_x}
            });
            auto nonempty = contractor.feasible(nonempty_domain, no_constraints);
            ARIADNE_TESTING_PROVE(definitely(nonempty.first));
            ARIADNE_TESTING_EXPECT_EQUAL(
                nonempty.second.dimension(), nonempty_domain.dimension());

            ExactBoxType empty_domain({ExactIntervalType::empty_interval()});
            auto empty = contractor.feasible(empty_domain, no_constraints);
            ARIADNE_TESTING_PROVE(not possibly(empty.first));
        }

        List<EffectiveScalarMultivariateFunction> x = EffectiveScalarMultivariateFunction::coordinates(1);
        EffectiveConstraint c = (x[0] - 2 <= 0);

        List<ValidatedConstraint> constraints;
        constraints.append(c);

        ARIADNE_TESTING_NOTE("constraints = " << (constraints));

        ExactBoxType domain1({
            {1.9375_x, 2.0_x}
        });
        auto feasible1 = contractor.feasible(domain1, constraints);
        ARIADNE_TESTING_PROVE(definitely(feasible1.first));
        ARIADNE_TESTING_EXPECT_EQUAL(feasible1.second.dimension(), domain1.dimension());
        ARIADNE_TESTING_PROVE(definitely(contractor.check_feasibility(
            domain1,
            ValidatedVectorMultivariateFunction(
                {constraints[0].function()}),
            ExactBoxType({constraints[0].bounds()}),
            feasible1.second)));

        ExactBoxType domain2({
            {2.015625_x, 2.5_x}
        });
        ARIADNE_TESTING_PROVE(!possibly(contractor.feasible(domain2, constraints).first));

        ExactBoxType domain3({
            {2.0_x, 2.015625_x}
        });
        ARIADNE_TESTING_PROVE(is_indeterminate(contractor.feasible(domain3, constraints).first));

        {
            auto xy = ValidatedScalarMultivariateFunction::coordinates(2);
            ValidatedVectorMultivariateFunction function({xy[0] + xy[1]});
            ExactBoxType domain({
                {0.0_x, 1.0_x},
                {0.0_x, 1.0_x}
            });
            ExactBoxType codomain({
                {0.19_pr, 0.21_pr}
            });
            ARIADNE_TESTING_PROVE(not definitely(contractor.check_feasibility(
                domain, function, codomain, domain.midpoint())));
            auto feasibility_result = contractor.feasible(domain, function, codomain);
            ARIADNE_TESTING_PROVE(possibly(feasibility_result.first));
            if (definitely(feasibility_result.first)) {
                ARIADNE_TESTING_EXPECT_EQUAL(feasibility_result.second.dimension(), domain.dimension());
                ARIADNE_TESTING_PROVE(definitely(contractor.check_feasibility(
                    domain, function, codomain, feasibility_result.second)));
            }
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-feasible] coupled constraints expose optimiser infeasibility");
            auto xy = ValidatedScalarMultivariateFunction::coordinates(2);
            ValidatedVectorMultivariateFunction function({xy[0],
                                                          xy[1],
                                                          xy[0] + xy[1]});
            ExactBoxType domain({
                {0.0_x, 1.0_x},
                {0.0_x, 1.0_x}
            });
            ExactBoxType codomain({
                { 0.0_x, 0.1_pr},
                { 0.0_x, 0.1_pr},
                {1.8_pr,  2.0_x}
            });
            UpperBoxType direct_image = Ariadne::apply(function, domain);
            for (SizeType i = 0u; i != codomain.size(); ++i) {
                ARIADNE_TESTING_PROVE(
                    possibly(intersect(direct_image[i], codomain[i])));
            }
            auto feasibility_result = contractor.feasible(domain, function, codomain);
            ARIADNE_TESTING_PROVE(not possibly(feasibility_result.first));
        }

        {
            ARIADNE_TESTING_NOTE("[constraint-feasible] dependency-hidden nonlinear infeasibility");
            auto validated_x = ValidatedScalarMultivariateFunction::coordinates(1);
            ValidatedVectorMultivariateFunction function({validated_x[0] * (1 - validated_x[0])});
            ExactBoxType domain({
                {0.0_x, 1.0_x}
            });
            ExactBoxType codomain({
                {0.32_pr, 0.32_pr}
            });
            UpperBoxType direct_image = Ariadne::apply(function, domain);
            ARIADNE_TESTING_PROVE(
                possibly(intersect(direct_image[0], codomain[0])));
            auto feasibility_result = contractor.feasible(domain, function, codomain);
            ARIADNE_TESTING_PROVE(not definitely(feasibility_result.first));
        }

        {
            auto coordinates = ValidatedScalarMultivariateFunction::coordinates(7);
            ValidatedScalarMultivariateFunction sum = coordinates[0];
            for (SizeType i = 1u; i != 7u; ++i) {
                sum = sum + coordinates[i];
            }
            ValidatedVectorMultivariateFunction function({sum});
            ExactBoxType domain({
                {0.0_x, 1.0_x},
                {0.0_x, 1.0_x},
                {0.0_x, 1.0_x},
                {0.0_x, 1.0_x},
                {0.0_x, 1.0_x},
                {0.0_x, 1.0_x},
                {0.0_x, 1.0_x}
            });
            ExactBoxType codomain({
                {1.175_pr, 1.425_pr}
            });
            auto feasibility_result = contractor.feasible(domain, function, codomain);
            ARIADNE_TESTING_PROVE(possibly(feasibility_result.first));
            if (definitely(feasibility_result.first)) {
                ARIADNE_TESTING_EXPECT_EQUAL(feasibility_result.second.dimension(), domain.dimension());
                ARIADNE_TESTING_PROVE(definitely(contractor.check_feasibility(
                    domain, function, codomain, feasibility_result.second)));
            }
        }
    }
};

int main(int argc, char** argv) {
    // Keep solver logging code exercised while routing test diagnostics
    // exclusively through Ariadne Testing.
#ifdef _WIN32
    Logger::instance().redirect_to_file("NUL");
#else
    Logger::instance().redirect_to_file("/dev/null");
#endif
    ARIADNE_TESTING_SUITE("test_constraint_solver");
    ARIADNE_TESTING_ADD_TEST(TestConstraintSolver, test_empty_reduce_inequality);
    ARIADNE_TESTING_ADD_TEST(TestConstraintSolver, test_empty_reduce_equality);
    ARIADNE_TESTING_ADD_TEST(TestConstraintSolver, test_empty_reduce_mixed);
    ARIADNE_TESTING_ADD_TEST(TestConstraintSolver, test_empty_hull_reduce);
    ARIADNE_TESTING_ADD_TEST(TestConstraintSolver, test_empty_box_reduce);
    ARIADNE_TESTING_ADD_TEST(TestConstraintSolver, test_hull_reduce);
    ARIADNE_TESTING_ADD_TEST(TestConstraintSolver, test_box_reduce);
    ARIADNE_TESTING_ADD_TEST(TestConstraintSolver, test_box_reduce_boundary_slices);
    ARIADNE_TESTING_ADD_TEST(TestConstraintSolver, test_pruning_well_definedness);
    ARIADNE_TESTING_ADD_TEST(TestConstraintSolver, test_composite_pruning_witness_preservation);
    ARIADNE_TESTING_ADD_TEST(TestConstraintSolver, test_partial_domain_pruning_witness_preservation);
    ARIADNE_TESTING_ADD_TEST(TestConstraintSolver, test_monotone_reduce);
    ARIADNE_TESTING_ADD_TEST(TestConstraintSolver, test_lyapunov_reduce);
    ARIADNE_TESTING_ADD_TEST(TestConstraintSolver, test_propagate);
    ARIADNE_TESTING_ADD_TEST(TestConstraintSolver, test_precompiled_propagate);
    ARIADNE_TESTING_ADD_TEST(TestConstraintSolver, test_reduce_certified_containment);
    ARIADNE_TESTING_ADD_TEST(TestConstraintSolver, test_reduce_per_query_procedure_lifetime);
    ARIADNE_TESTING_ADD_TEST(TestConstraintSolver, test_prepared_reduce_matches_production);
    ARIADNE_TESTING_ADD_TEST(TestConstraintSolver, test_reduce_edge_cases);
    ARIADNE_TESTING_ADD_TEST(TestConstraintSolver, test_check_feasibility);
    ARIADNE_TESTING_ADD_TEST(TestConstraintSolver, test_logging_branches);
    ARIADNE_TESTING_ADD_TEST(TestConstraintSolver, test_feasible);
    ARIADNE_TESTING_ADD_TEST(TestConstraintSolver, test_split);
    ARIADNE_TESTING_ADD_TEST(TestConstraintSolver, test_virtual_destruction);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
