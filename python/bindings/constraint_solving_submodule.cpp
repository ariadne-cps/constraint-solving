#include "pybind11.hpp"
#include "function-utilities.hpp"

#include <optional>
#include <utility>
#include <vector>

#include "algebra/algebra.hpp"
#include "numeric/numeric.hpp"
#include "function/function.hpp"
#include "function/constraint.hpp"
#include "function/procedure.hpp"
#include "function/taylor_function.hpp"
#include "geometry/box.hpp"
#include "constraint-solving/constraint_solver.hpp"

using namespace Ariadne;

Void export_constraint(pybind11::module& module)
{
    pybind11::class_<EffectiveConstraint> effective_nonlinear_constraint_class(module,"EffectiveConstraint");
    effective_nonlinear_constraint_class.def(pybind11::init<EffectiveConstraint>());
    effective_nonlinear_constraint_class.def(pybind11::init<EffectiveNumber,EffectiveScalarMultivariateFunction,EffectiveNumber>());
    effective_nonlinear_constraint_class.def(pybind11::init([](EffectiveScalarMultivariateFunction const& f, EffectiveNumber const& u){return EffectiveConstraint(EffectiveNumber(-infty),f,u);}));
    effective_nonlinear_constraint_class.def(pybind11::init([](EffectiveNumber const& l, EffectiveScalarMultivariateFunction const& f){return EffectiveConstraint(l,f,EffectiveNumber(infty));}));
    effective_nonlinear_constraint_class.def(pybind11::init<EffectiveScalarMultivariateFunction,EffectiveNumber>());

    effective_nonlinear_constraint_class.def("__str__",&__cstr__<EffectiveConstraint>);

    pybind11::class_<ValidatedConstraint> validated_nonlinear_constraint_class(module,"ValidatedConstraint");
    validated_nonlinear_constraint_class.def(pybind11::init<ValidatedNumber,ValidatedScalarMultivariateFunction,ValidatedNumber
    >());
    validated_nonlinear_constraint_class.def(pybind11::init<ValidatedConstraint>());
    validated_nonlinear_constraint_class.def(pybind11::init<EffectiveConstraint>());
    validated_nonlinear_constraint_class.def("lower_bound", &ValidatedConstraint::lower_bound);
    validated_nonlinear_constraint_class.def("upper_bound", &ValidatedConstraint::upper_bound);
    validated_nonlinear_constraint_class.def("function", (const ValidatedScalarMultivariateFunction&(ValidatedConstraint::*)()const) &ValidatedConstraint::function);
    validated_nonlinear_constraint_class.def("__str__",&__cstr__<ValidatedConstraint>);
}

Void export_constraint_solving(pybind11::module& module)
{
    // Value objects make incremental propagation usable without exposing
    // mutable C++ reference parameters to Python.
    pybind11::class_<ConstraintPropagationStatistics>(module,"ConstraintPropagationStatistics")
        .def(pybind11::init<>())
        .def_readonly("hull_rounds", &ConstraintPropagationStatistics::hull_rounds)
        .def_readonly("hull_effective", &ConstraintPropagationStatistics::hull_effective)
        .def_readonly("hull_contractor_calls", &ConstraintPropagationStatistics::hull_contractor_calls)
        .def_readonly("hull_agenda_pushes", &ConstraintPropagationStatistics::hull_agenda_pushes)
        .def_readonly("hull_agenda_pops", &ConstraintPropagationStatistics::hull_agenda_pops)
        .def_readonly("hull_agenda_effective_calls", &ConstraintPropagationStatistics::hull_agenda_effective_calls)
        .def_readonly("hull_procedure_builds", &ConstraintPropagationStatistics::hull_procedure_builds)
        .def_readonly("hull_procedure_build_seconds", &ConstraintPropagationStatistics::hull_procedure_build_seconds)
        .def_readonly("hull_contraction_seconds", &ConstraintPropagationStatistics::hull_contraction_seconds)
        .def_readonly("hull_direct_rejection_seconds", &ConstraintPropagationStatistics::hull_direct_rejection_seconds)
        .def_readonly("hull_temporary_allocation_seconds", &ConstraintPropagationStatistics::hull_temporary_allocation_seconds)
        .def_readonly("hull_forward_execution_seconds", &ConstraintPropagationStatistics::hull_forward_execution_seconds)
        .def_readonly("hull_backward_propagation_seconds", &ConstraintPropagationStatistics::hull_backward_propagation_seconds)
        .def_readonly("shaving_rounds", &ConstraintPropagationStatistics::shaving_rounds)
        .def_readonly("shaving_effective", &ConstraintPropagationStatistics::shaving_effective)
        .def_readonly("shaving_coordinate_attempts", &ConstraintPropagationStatistics::shaving_coordinate_attempts)
        .def_readonly("shaving_coordinate_effective", &ConstraintPropagationStatistics::shaving_coordinate_effective)
        .def_readonly("shaving_dependency_skipped", &ConstraintPropagationStatistics::shaving_dependency_skipped)
        .def_readonly("shaving_adaptive_skipped", &ConstraintPropagationStatistics::shaving_adaptive_skipped)
        .def_readonly("shaving_refresh_rounds", &ConstraintPropagationStatistics::shaving_refresh_rounds)
        .def_readonly("shaving_active_rounds", &ConstraintPropagationStatistics::shaving_active_rounds)
        .def_readonly("shaving_function_evaluations", &ConstraintPropagationStatistics::shaving_function_evaluations)
        .def_readonly("shaving_seconds", &ConstraintPropagationStatistics::shaving_seconds)
        .def_readonly("monotone_rounds", &ConstraintPropagationStatistics::monotone_rounds)
        .def_readonly("monotone_effective", &ConstraintPropagationStatistics::monotone_effective);

    pybind11::class_<ConstraintPropagationConstraint>(module,"ConstraintPropagationConstraint")
        .def(pybind11::init([](ValidatedScalarMultivariateFunction function,
                                ExactIntervalType bounds,
                                Bool strict_lower, Bool strict_upper) {
            return ConstraintPropagationConstraint{
                std::move(function),std::move(bounds),{},strict_lower,strict_upper};
        }),pybind11::arg("function"),pybind11::arg("bounds"),
           pybind11::arg("strict_lower")=false,
           pybind11::arg("strict_upper")=false)
        .def_readwrite("function",&ConstraintPropagationConstraint::function)
        .def_readwrite("bounds",&ConstraintPropagationConstraint::bounds)
        .def_readwrite("derivatives",&ConstraintPropagationConstraint::derivatives)
        .def_readwrite("strict_lower",&ConstraintPropagationConstraint::strict_lower)
        .def_readwrite("strict_upper",&ConstraintPropagationConstraint::strict_upper);

    pybind11::class_<ConstraintSolver> solver_class(module,"ConstraintSolver");
    solver_class.def(pybind11::init<>());
    solver_class.def("feasible",
        (Pair<ValidatedKleenean,ConstraintSolver::ExactPointType>
        (ConstraintSolver::*)(const ExactBoxType&,const ValidatedVectorMultivariateFunction&,const ExactBoxType&)const)
        &ConstraintSolver::feasible);
    solver_class.def("feasible",
        (Pair<ValidatedKleenean,ConstraintSolver::ExactPointType>
        (ConstraintSolver::*)(const ExactBoxType&,const List<ValidatedConstraint>&)const)
        &ConstraintSolver::feasible);
    solver_class.def("check_feasibility",&ConstraintSolver::check_feasibility);


    // The C++ reducer modifies an UpperBoxType by reference; the Python
    // convenience overload returns both the proof flag and reduced box.
    solver_class.def("reduce_box",
        [](ConstraintSolver const& solver, ExactBoxType const& domain,
           ValidatedVectorMultivariateFunction const& function,
           ExactBoxType const& codomain) {
            UpperBoxType reduced(domain);
            Bool rejected=solver.reduce(reduced,function,codomain);
            return pybind11::make_tuple(rejected,reduced);
        });

    solver_class.def("propagate_precompiled",
        [](ConstraintSolver const& solver, ExactBoxType const& domain,
           std::vector<ConstraintPropagationConstraint> const& constraints,
           Bool monotone, Bool shaving, Bool hull) {
            UpperBoxType reduced(domain);
            ConstraintPropagationStatistics statistics;
            Bool rejected=solver.propagate(
                reduced,constraints,monotone,statistics,shaving,hull);
            return pybind11::make_tuple(rejected,reduced,statistics);
        },pybind11::arg("domain"),pybind11::arg("constraints"),
          pybind11::arg("monotone")=false,pybind11::arg("shaving")=true,
          pybind11::arg("hull")=true);

    solver_class.def("propagate_constraints",
        [](ConstraintSolver const& solver, ExactBoxType const& domain,
           List<ValidatedConstraint> const& constraints,
           Bool shaving, Bool hull) {
            UpperBoxType reduced(domain);
            ConstraintPropagationStatistics statistics;
            Bool rejected=solver.propagate(
                reduced,constraints,statistics,shaving,hull);
            return pybind11::make_tuple(rejected,reduced,statistics);
        },pybind11::arg("domain"),pybind11::arg("constraints"),
          pybind11::arg("shaving")=true,pybind11::arg("hull")=true);

    solver_class.def("propagate_hull_once",
        [](ConstraintSolver const& solver, ExactBoxType const& domain,
           std::vector<ConstraintPropagationConstraint> const& constraints) {
            UpperBoxType reduced(domain);
            ConstraintPropagationStatistics statistics;
            Bool rejected=solver.propagate_hull_once(
                reduced,constraints,statistics);
            return pybind11::make_tuple(rejected,reduced,statistics);
        });

    solver_class.def("hull_reduce_box",
        [](ConstraintSolver const& solver, ExactBoxType const& domain,
           ValidatedVectorMultivariateFunction const& function,
           ExactBoxType const& codomain) {
            UpperBoxType reduced(domain);
            Bool rejected=solver.hull_reduce(reduced,function,codomain);
            return pybind11::make_tuple(rejected,reduced);
        });
    solver_class.def("box_reduce_box",
        [](ConstraintSolver const& solver, ExactBoxType const& domain,
           ValidatedScalarMultivariateFunction const& function,
           ExactIntervalType const& bounds, SizeType coordinate) {
            UpperBoxType reduced(domain);
            Bool rejected=solver.box_reduce(reduced,function,bounds,coordinate);
            return pybind11::make_tuple(rejected,reduced);
        });
    solver_class.def("monotone_reduce_box",
        [](ConstraintSolver const& solver, ExactBoxType const& domain,
           ValidatedScalarMultivariateFunction const& function,
           ExactIntervalType const& bounds, SizeType coordinate) {
            UpperBoxType reduced(domain);
            Bool rejected=solver.monotone_reduce(reduced,function,bounds,coordinate);
            return pybind11::make_tuple(rejected,reduced);
        });

    // Keep C++-style in-place calls available for existing clients.
    solver_class.def("hull_reduce",
        (Bool(ConstraintSolver::*)(UpperBoxType&,const ValidatedScalarMultivariateFunction&,const ExactIntervalType&)const)
        &ConstraintSolver::hull_reduce);
    solver_class.def("box_reduce",
        (Bool(ConstraintSolver::*)(UpperBoxType&,const ValidatedScalarMultivariateFunction&,const ExactIntervalType&,SizeType)const)
        &ConstraintSolver::box_reduce);
    solver_class.def("monotone_reduce",
        (Bool(ConstraintSolver::*)(UpperBoxType&,const ValidatedScalarMultivariateFunction&,const ExactIntervalType&,SizeType)const)
        &ConstraintSolver::monotone_reduce);
    solver_class.def("reduce",
        (Bool(ConstraintSolver::*)(UpperBoxType&,const List<ValidatedConstraint>&)const)
        &ConstraintSolver::reduce);
    solver_class.def("reduce",
        (Bool(ConstraintSolver::*)(UpperBoxType&,const ValidatedVectorMultivariateFunction&,const ExactBoxType&)const)
        &ConstraintSolver::reduce);
}


void constraint_solving_submodule(pybind11::module& module) {
    export_constraint(module);
    export_constraint_solving(module);
}
