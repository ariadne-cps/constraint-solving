#include <iostream>
#include "ariadne-constraint-solving.hpp"

using namespace Ariadne;

void propagation_demo() {
    ConstraintSolver propagator;
    auto x=ValidatedScalarMultivariateFunction::coordinates(2);
    ValidatedVectorMultivariateFunction function({x[0]+x[1]});
    ExactBoxType domain({{0.0_x,2.0_x},{0.0_x,2.0_x}});
    ExactBoxType bounds({{0.875_x,1.125_x}});

    // Basic interval contraction, retaining all feasible solutions.
    UpperBoxType contracted(domain);
    Bool rejected=propagator.reduce(contracted,function,bounds);
    std::cout<<"rejected="<<rejected<<", reduced domain="<<contracted<<"\n";

    // A single cached/full-scan hull pass and its statistics.
    std::vector<ConstraintPropagationConstraint> constraints({
        {x[0]+x[1],ExactIntervalType(0.875_x,1.125_x)}
    });
    UpperBoxType swept(domain);
    ConstraintPropagationStatistics sweep_stats;
    Bool sweep_rejected=propagator.propagate_hull_once(
        swept,constraints,sweep_stats);
    std::cout<<"one hull pass: rejected="<<sweep_rejected
             <<", box="<<swept
             <<", contractor calls="<<sweep_stats.hull_contractor_calls<<"\n";

    // Full fixed-point propagation with the same numerical constraint.
    UpperBoxType propagated(domain);
    ConstraintPropagationStatistics stats;
    Bool propagation_rejected=propagator.propagate(
        propagated,constraints,false,stats);
    std::cout<<"propagation: rejected="<<propagation_rejected
             <<", box="<<propagated
             <<", hull rounds="<<stats.hull_rounds
             <<", shaving rounds="<<stats.shaving_rounds<<"\n";

    // Strict inequalities can prove rejection without emptying the box.
    auto y=ValidatedScalarMultivariateFunction::coordinates(1);
    UpperBoxType boundary=ExactBoxType({{0.0_x,0.0_x}});
    std::vector<ConstraintPropagationConstraint> strict({
        {y[0],ExactIntervalType(0.0_x,1.0_x),{},true,false}
    });
    ConstraintPropagationStatistics strict_stats;
    Bool strict_rejected=propagator.propagate(
        boundary,strict,false,strict_stats,false,false);
    std::cout<<"strict x>0 at x=0: rejected="<<strict_rejected
             <<", domain="<<boundary<<"\n";
}

// The same solver also handles numerical feasibility and point checks.
void feasibility_demo() {
    ConstraintSolver solver;
    auto x=ValidatedScalarMultivariateFunction::coordinates(2);
    ValidatedVectorMultivariateFunction function({x[0]+x[1]});
    ExactBoxType domain({{0.0_x,1.0_x},{0.0_x,1.0_x}});
    ExactBoxType codomain({{0.19_pr,0.21_pr}});

    // A solver may validate existence near a numerical candidate.
    auto result=solver.feasible(domain,function,codomain);
    std::cout<<"feasibility="<<result.first
             <<", candidate="<<result.second<<"\n";

    // Checking a specific interior point is a separate validation task.
    ConstraintSolver::ExactPointType witness({
        // Decimal 0.1 is not exactly representable as a binary64 value.
        FloatDP(0.1_pr,dp),FloatDP(0.1_pr,dp)
    });
    std::cout<<"point feasibility="
             <<solver.check_feasibility(domain,function,codomain,witness)<<"\n";

    // Direct interval evaluation certifies this inconsistent codomain.
    ExactBoxType impossible({{3.0_x,4.0_x}});
    auto inconsistent=solver.feasible(domain,function,impossible);
    std::cout<<"inconsistent feasibility="<<inconsistent.first<<"\n";
}

int main() {
    propagation_demo();
    feasibility_demo();
    return 0;
}
