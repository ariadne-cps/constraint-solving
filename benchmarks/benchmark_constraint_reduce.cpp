/***************************************************************************
 *            benchmark_constraint_reduce.cpp
 *
 *  Standalone solver profiling; no dependency on set or paving.
 *  This file is part of Ariadne Constraint Solving, GPL-3.0-or-later.
 ***************************************************************************/

#include "constraint-solving/constraint_solver.hpp"
#include "function/function.hpp"
#include "function/constraint.hpp"
#include "function/procedure.hpp"
#include "geometry/box.hpp"
#include "foundation/logical.hpp"
#include "logging/logging.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace Ariadne;

namespace {

constexpr char const* known_shapes[]={
    "aligned-box","rotated-box","aligned-ellipsoid",
    "rotated-ellipsoid","ball"
};

struct Selection {
    std::vector<std::string> shapes{
        std::begin(known_shapes),std::end(known_shapes)};
    unsigned samples=7u;
    unsigned min_sample_ms=12u;
    unsigned max_iterations=256u;
    Bool verify_only=false;
    Bool compare_contained_skip=false;
    Bool compare_dependency_skip=false;
    Bool profile_residual_shaving=false;
    Bool compare_adaptive_shaving=false;
    Bool compare_adaptive_cache=false;
    unsigned adaptive_depth=8u;
    unsigned adaptive_budget=8u;
    unsigned adaptive_resolution_budget=64u;
    std::string trace_prefix;
};

unsigned unsigned_argument(std::string const& text) {
    if(text.empty() or text.front()=='-')
        throw std::invalid_argument("expected nonnegative integer");
    std::size_t used=0u;
    unsigned long long number=std::stoull(text,&used);
    if(used!=text.size() or number>std::numeric_limits<unsigned>::max())
        throw std::invalid_argument("integer out of range");
    return static_cast<unsigned>(number);
}

std::vector<std::string> shapes_argument(std::string const& input) {
    if(input.empty() or input.front()==',' or input.back()==',')
        throw std::invalid_argument("invalid shape list");
    std::vector<std::string> result;
    std::size_t start=0u;
    do {
        std::size_t end=input.find(',',start);
        std::string name=input.substr(
            start,end==std::string::npos?end:end-start);
        if(std::find(std::begin(known_shapes),std::end(known_shapes),name)
                ==std::end(known_shapes) or
           std::find(result.begin(),result.end(),name)!=result.end())
            throw std::invalid_argument("unknown or duplicate shape: "+name);
        result.push_back(name);
        if(end==std::string::npos) break;
        start=end+1u;
    } while(true);
    return result;
}

// This component has no threading dependency. Keep its tiny CLI local,
// while accepting --concurrency 0 for parity with the paving benchmarks.
void print_usage(char const* program) {
    std::cout<<"Usage: "<<program<<" [options]\n"
             <<"  -g, --shapes NAMES         Comma-separated shape list\n"
             <<"  -n, --samples N            Odd number of samples (>=3)\n"
             <<"  --min-sample-ms N          Calibration threshold in ms\n"
             <<"  --max-iterations N         Maximum calibrated repetitions\n"
             <<"  --concurrency 0            Sequential mode (only supported value)\n"
             <<"  --trace-prefix PREFIX      Replay exported paving d4/f3 boxes\n"
             <<"  --compare-contained-skip   Compare legacy full scan vs production\n"
             <<"  --compare-dependency-skip  Try static procedure-dependency shaving\n"
             <<"  --profile-residual-shaving Profile P1 surviving box calls (0..8 evals)\n"
             <<"  --compare-adaptive-shaving Compare 8-slice vs adaptive shaving\n"
             <<"  --compare-adaptive-cache   Benchmark reuse of tested dyadic cells\n"
             <<"  --adaptive-depth N         Maximum dyadic depth (1..12; default 8)\n"
             <<"  --adaptive-budget N        Tight evaluation budget (>=2; default 8)\n"
             <<"  --adaptive-resolution-budget N Evaluation cap for depth-limited mode (default 64)\n"
             <<"  --verify-only              Validate without timing\n"
             <<"  -h, --help                 Print this help\n";
}

Bool parse_arguments(int argc,char const* argv[],Selection& selection) {
    for(int i=1;i<argc;++i) {
        std::string argument(argv[i]);
        if(argument=="-h" or argument=="--help") {
            print_usage(argv[0]);
            return false;
        }
        if(argument=="--verify-only") {
            selection.verify_only=true;
            continue;
        }
        if(argument=="--compare-contained-skip") {
            selection.compare_contained_skip=true;
            continue;
        }
        if(argument=="--compare-dependency-skip") {
            selection.compare_dependency_skip=true;
            continue;
        }
        if(argument=="--profile-residual-shaving") {
            selection.profile_residual_shaving=true;
            continue;
        }
        if(argument=="--compare-adaptive-shaving") {
            selection.compare_adaptive_shaving=true;
            continue;
        }
        if(argument=="--compare-adaptive-cache") {
            selection.compare_adaptive_cache=true;
            continue;
        }
        if(i+1>=argc) {
            throw std::invalid_argument("missing value for "+argument);
        }
        std::string value(argv[++i]);
        if(argument=="-g" or argument=="--shapes") {
            selection.shapes=shapes_argument(value);
        } else if(argument=="-n" or argument=="--samples") {
            unsigned number=unsigned_argument(value);
            if(number<3u or number%2u==0u)
                throw std::invalid_argument("samples must be odd and >=3");
            selection.samples=number;
        } else if(argument=="--min-sample-ms") {
            unsigned number=unsigned_argument(value);
            if(number==0u or number>10000u)
                throw std::invalid_argument("min-sample-ms out of range");
            selection.min_sample_ms=number;
        } else if(argument=="--max-iterations") {
            unsigned number=unsigned_argument(value);
            if(number==0u)
                throw std::invalid_argument("max-iterations must be positive");
            selection.max_iterations=number;
        } else if(argument=="--adaptive-depth") {
            unsigned depth=unsigned_argument(value);
            if(depth<1u or depth>12u)
                throw std::invalid_argument("adaptive-depth must be in 1..12");
            selection.adaptive_depth=depth;
        } else if(argument=="--adaptive-budget") {
            unsigned budget=unsigned_argument(value);
            if(budget<2u or budget>4096u)
                throw std::invalid_argument("adaptive-budget must be in 2..4096");
            selection.adaptive_budget=budget;
        } else if(argument=="--adaptive-resolution-budget") {
            unsigned budget=unsigned_argument(value);
            if(budget<2u or budget>4096u)
                throw std::invalid_argument(
                    "adaptive-resolution-budget must be in 2..4096");
            selection.adaptive_resolution_budget=budget;
        } else if(argument=="--trace-prefix") {
            if(value.empty())
                throw std::invalid_argument("trace prefix must not be empty");
            selection.trace_prefix=value;
        } else if(argument=="--concurrency") {
            if(value!="0") {
                throw std::invalid_argument(
                    "only --concurrency 0 is supported by this sequential benchmark");
            }
        } else {
            throw std::invalid_argument("unknown argument: "+argument);
        }
    }
    return true;
}

List<EffectiveConstraint> shape_constraints(std::string const& shape) {
    auto x=EffectiveScalarMultivariateFunction::coordinates(4u);
    auto first=x[0]-0.5_x;
    auto second=x[1]-0.5_x;
    auto rotated_first=3*first+4*second;
    auto rotated_second=-4*first+3*second;
    List<EffectiveConstraint> result;

    if(shape=="aligned-box") {
        result.append(first>=-0.25_x);
        result.append(first<=0.25_x);
        result.append(second>=-0.125_x);
        result.append(second<=0.125_x);
    } else if(shape=="rotated-box") {
        result.append(rotated_first>=-1.25_x);
        result.append(rotated_first<=1.25_x);
        result.append(rotated_second>=-0.625_x);
        result.append(rotated_second<=0.625_x);
    } else {
        EffectiveScalarMultivariateFunction quadratic=
            first*first*16+second*second*16;
        if(shape=="aligned-ellipsoid") {
            quadratic=first*first*16+second*second*64;
        } else if(shape=="rotated-ellipsoid") {
            quadratic=rotated_first*rotated_first*16
                +rotated_second*rotated_second*64;
        } else if(shape!="ball") {
            throw std::invalid_argument("unknown geometry");
        }
        for(unsigned i=2u;i<4u;++i) {
            auto centered=x[i]-0.5_x;
            quadratic=quadratic+centered*centered*
                (shape=="rotated-ellipsoid"?400:16);
        }
        result.append(
            quadratic<=(shape=="rotated-ellipsoid"?25.0_x:1.0_x));
    }
    if(shape=="aligned-box" or shape=="rotated-box") {
        for(unsigned i=2u;i<4u;++i) {
            result.append(x[i]>=0.25_x);
            result.append(x[i]<=0.75_x);
        }
    }
    return result;
}

// The same 4D identity map and parameter constraints as paving's five
// geometric workloads, but expressed directly as a solver problem.
struct SolverProblem {
    ValidatedVectorMultivariateFunction function;
    UpperBoxType initial_domain;
    std::vector<ExactBoxType> codomains;
    SizeType constraints_count;
};

// Fixed, reproducible dyadic boxes. These are NOT the adaptive candidate
// cells from paving; any conclusion applies first to this input sample.
SolverProblem make_problem(std::string const& shape) {
    List<EffectiveConstraint> constraints=shape_constraints(shape);
    auto x=EffectiveScalarMultivariateFunction::coordinates(4u);
    EffectiveVectorMultivariateFunction effective(
        4u+constraints.size(),EuclideanDomain(4u));
    for(SizeType i=0u;i<4u;++i) effective[i]=x[i];
    for(SizeType i=0u;i<constraints.size();++i)
        effective[i+4u]=constraints[i].function();
    // The public ConstrainedImageSet::separated() passes this effective
    // vector to reduce(), whose parameter is a validated vector.
    ValidatedVectorMultivariateFunction function=effective;

    ExactIntervalType intervals[]={
        ExactIntervalType(0.0_x,0.25_x),
        ExactIntervalType(0.25_x,0.5_x),
        ExactIntervalType(0.5_x,0.75_x),
        ExactIntervalType(0.75_x,1.0_x)
    };
    std::vector<ExactBoxType> codomains;
    codomains.reserve(32u);
    auto append=[&](ExactBoxType const& box) {
        ExactBoxType codomain(4u+constraints.size());
        for(SizeType i=0u;i<4u;++i) codomain[i]=box[i];
        for(SizeType i=0u;i<constraints.size();++i)
            codomain[4u+i]=constraints[i].bounds();
        codomains.push_back(std::move(codomain));
    };
    for(unsigned i=0u;i<4u;++i) {
        for(unsigned j=0u;j<4u;++j) {
            append(ExactBoxType({
                intervals[i],intervals[j],
                intervals[(i+j)%4u],intervals[(i+2u*j)%4u]}));
            append(ExactBoxType({
                intervals[i],intervals[j],
                ExactIntervalType(0.25_x,0.75_x),
                ExactIntervalType(0.25_x,0.75_x)}));
        }
    }
    UpperBoxType initial=ExactBoxType({
        {0.0_x,1.0_x},{0.0_x,1.0_x},
        {0.0_x,1.0_x},{0.0_x,1.0_x}});
    return {std::move(function),std::move(initial),
            std::move(codomains),constraints.size()};
}


struct TraceRecord {
    ExactBoxType box;
    Bool separated;
};

std::vector<TraceRecord> read_adaptive_trace(
        std::string const& prefix,std::string const& shape) {
    std::string path=prefix+"-"+shape+"-d4-f3.trace";
    std::ifstream in(path);
    if(not in) throw std::runtime_error("cannot open paving trace: "+path);
    std::string header,recorded_shape;
    unsigned dimension=0u,fineness=0u;
    SizeType count=0u;
    if(not (in>>header>>recorded_shape>>dimension>>fineness>>count)
       or header!="ARIADNE_PAVING_SOLVER_TRACE_V1"
       or recorded_shape!=shape or dimension!=4u or fineness!=3u
       or count==0u or count>1000000u) {
        throw std::runtime_error("invalid paving trace header: "+path);
    }
    std::vector<TraceRecord> result;
    result.reserve(count);
    auto parse_number=[](std::string const& token)->double {
        if(token.empty()) throw std::runtime_error("empty trace endpoint");
        char* end=nullptr;
        double value=std::strtod(token.c_str(),&end);
        if(end==token.c_str() or *end!='\0' or not std::isfinite(value))
            throw std::runtime_error("invalid trace endpoint: "+token);
        return value;
    };
    // Parse strings with strtod; standard iostream floating input does not
    // reliably accept C99 hexadecimal-float notation across platforms.
    for(SizeType i=0u;i<count;++i) {
        unsigned separated=0u;
        if(not (in>>separated) or separated>1u)
            throw std::runtime_error("invalid trace separation bit at query "
                                    +std::to_string(i));
        ExactBoxType box(4u);
        for(SizeType j=0u;j<4u;++j) {
            std::string lower,upper;
            if(not (in>>lower>>upper))
                throw std::runtime_error("truncated paving trace");
            double l=parse_number(lower);
            double u=parse_number(upper);
            if(l>u) throw std::runtime_error("reversed trace interval");
            box[j]=ExactIntervalType(ExactDouble(l),ExactDouble(u));
        }
        result.push_back({std::move(box),separated!=0u});
    }
    std::string trailing;
    if(in>>trailing)
        throw std::runtime_error("unexpected trailing trace payload");
    return result;
}

// A call is "changed" when its targeted interval differs from the input,
// including certified emptiness. Unchanged calls are not necessarily
// redundant: they may still be required to establish a validated result.
struct BoxCallStatistics {
    SizeType calls=0u;
    SizeType evaluations=0u;
    SizeType changed=0u;
    SizeType unchanged=0u;
    SizeType emptied=0u;
    SizeType zero_evaluations=0u;
    double seconds=0.0;

    void add(SizeType n,Bool modified,Bool empty,double duration) {
        ++calls;
        evaluations+=n;
        changed+=static_cast<SizeType>(modified);
        unchanged+=static_cast<SizeType>(not modified);
        emptied+=static_cast<SizeType>(empty);
        zero_evaluations+=static_cast<SizeType>(n==0u);
        seconds+=duration;
    }
};

struct ReduceTraceStatistics {
    SizeType queries=0u;
    SizeType proven_empty=0u;
    SizeType rounds=0u;
    SizeType hull_calls=0u;
    SizeType box_calls=0u;
    SizeType stopped_empty_hull=0u;
    SizeType stopped_empty_box=0u;
    SizeType stopped_threshold=0u;
    double hull_seconds=0.0;
    double box_seconds=0.0;
    double total_seconds=0.0;

    BoxCallStatistics box_activity;
    std::vector<BoxCallStatistics> by_output;
    std::vector<BoxCallStatistics> by_variable;
    std::vector<BoxCallStatistics> by_pair;

    void prepare(SizeType outputs,SizeType variables) {
        by_output.resize(outputs);
        by_variable.resize(variables);
        by_pair.resize(outputs*variables);
    }
    void record_box(SizeType output,SizeType variable,SizeType count,
                    SizeType evaluations,Bool changed,Bool empty,
                    double elapsed_seconds) {
        box_activity.add(evaluations,changed,empty,elapsed_seconds);
        by_output[output].add(
            evaluations,changed,empty,elapsed_seconds);
        by_variable[variable].add(
            evaluations,changed,empty,elapsed_seconds);
        by_pair[output*count+variable].add(
            evaluations,changed,empty,elapsed_seconds);
    }
};

// Frozen pre-containment full-scan reducer. Benchmark-only reference:
// the production reducer must not be used for this baseline.
Bool reduce_full_scan_reference(ConstraintSolver const& solver,
                                UpperBoxType& domain,
                                ValidatedVectorMultivariateFunction const& function,
                                ExactBoxType const& codomain) {
    const ExactDouble MINIMUM_REDUCTION=0.75_x;
    if(definitely(domain.is_empty())) return true;
    FloatDPUpperBound magnitude={0u,dp};
    for(SizeType j=0u;j<domain.size();++j)
        magnitude+=domain[j].width();
    FloatDPUpperBound previous=magnitude;
    do {
        solver.hull_reduce(domain,function,codomain);
        if(definitely(domain.is_empty())) return true;
        for(SizeType i=0u;i<codomain.size();++i) {
            for(SizeType j=0u;j<domain.size();++j) {
                solver.box_reduce(domain,function[i],codomain[i],j);
                if(definitely(domain.is_empty())) return true;
            }
        }
        previous=magnitude;
        magnitude=0u;
        for(SizeType j=0u;j<domain.size();++j)
            magnitude+=domain[j].width();
    } while(magnitude.raw()
            < mul(near,previous.raw(),FloatDP(MINIMUM_REDUCTION,dp)));
    return false;
}

void traced_reduce(ConstraintSolver const& solver,UpperBoxType& domain,
                   ValidatedVectorMultivariateFunction const& function,
                   ExactBoxType const& codomain,ReduceTraceStatistics& stats) {
    using Clock=std::chrono::steady_clock;
    auto elapsed=[](Clock::time_point started)->double {
        return std::chrono::duration<double>(Clock::now()-started).count();
    };
    ++stats.queries;
    auto start=Clock::now();
    const ExactDouble MINIMUM_REDUCTION=0.75_x;
    if(definitely(domain.is_empty())) {
        ++stats.proven_empty;
        stats.total_seconds+=elapsed(start);
        return;
    }
    FloatDPUpperBound magnitude={0u,dp};
    for(SizeType j=0u;j<domain.size();++j)
        magnitude+=domain[j].width();
    FloatDPUpperBound previous=magnitude;
    do {
        ++stats.rounds;
        auto hstart=Clock::now();
        ++stats.hull_calls;
        solver.hull_reduce(domain,function,codomain);
        stats.hull_seconds+=elapsed(hstart);
        if(definitely(domain.is_empty())) {
            ++stats.stopped_empty_hull;
            ++stats.proven_empty;
            stats.total_seconds+=elapsed(start);
            return;
        }
        for(SizeType i=0u;i<codomain.size();++i) {
            for(SizeType j=0u;j<domain.size();++j) {
                // Use the solver's existing evaluation-count overload.
                // Its numerical algorithm is identical to the four-argument
                // overload used by production reduce().
                UpperIntervalType before=domain[j];
                SizeType function_evaluations=0u;
                auto bstart=Clock::now();
                ++stats.box_calls;
                solver.box_reduce(
                    domain,function[i],codomain[i],j,function_evaluations);
                double seconds=elapsed(bstart);
                stats.box_seconds+=seconds;
                Bool empty=definitely(domain.is_empty());
                stats.record_box(
                    i,j,domain.size(),function_evaluations,
                    not same(before,domain[j]),empty,seconds);
                if(empty) {
                    ++stats.stopped_empty_box;
                    ++stats.proven_empty;
                    stats.total_seconds+=elapsed(start);
                    return;
                }
            }
        }
        previous=magnitude;
        magnitude=0u;
        for(SizeType j=0u;j<domain.size();++j)
            magnitude+=domain[j].width();
    } while(magnitude.raw()
            < mul(near,previous.raw(),FloatDP(MINIMUM_REDUCTION,dp)));
    ++stats.stopped_threshold;
    stats.total_seconds+=elapsed(start);
}

void profile_adaptive_trace(std::string const& prefix,
                            std::string const& shape,Bool verify_only) {
    SolverProblem problem=make_problem(shape);
    std::vector<TraceRecord> records=read_adaptive_trace(prefix,shape);
    ConstraintSolver solver;
    ReduceTraceStatistics stats;
    stats.prepare(problem.function.result_size(),problem.initial_domain.size());
    SizeType mismatch=0u;
    auto sweep_start=std::chrono::steady_clock::now();
    for(TraceRecord const& entry:records) {
        // Reuse exactly the geometry and constraints from this shape; only
        // the query-box portion varies with the exported adaptive trace.
        ExactBoxType codomain=problem.codomains.front();
        for(SizeType i=0u;i<4u;++i)
            codomain[i]=entry.box[i];
        UpperBoxType original=problem.initial_domain;
        reduce_full_scan_reference(solver,original,problem.function,codomain);
        UpperBoxType instrumented=problem.initial_domain;
        traced_reduce(solver,instrumented,problem.function,codomain,stats);
        if(not same(original,instrumented))
            throw std::runtime_error(
                "instrumented legacy reducer differs from full-scan reference for "+shape);
        UpperBoxType production=problem.initial_domain;
        solver.reduce(production,problem.function,codomain);
        if(not same(original,production))
            throw std::runtime_error(
                "production reducer differs from full-scan reference for "+shape);
        Bool empty=definitely(original.is_empty());
        if(empty!=entry.separated) ++mismatch;
    }
    double sweep_ms=1000.0*std::chrono::duration<double>(
        std::chrono::steady_clock::now()-sweep_start).count();
    if(mismatch!=0u) throw std::runtime_error(
        "solver replay disagrees with paving separation trace: "+shape+
        " mismatches="+std::to_string(mismatch));
    if(stats.queries!=records.size()
       or stats.stopped_empty_hull+stats.stopped_empty_box+
          stats.stopped_threshold!=stats.queries)
        throw std::runtime_error("inconsistent solver trace counters");
    auto const& a=stats.box_activity;
    if(a.calls!=stats.box_calls
       or a.changed+a.unchanged!=a.calls
       or a.emptied>a.changed or a.zero_evaluations>a.calls)
        throw std::runtime_error("inconsistent box contractor counters");
    auto check_group=[&](std::vector<BoxCallStatistics> const& group) {
        SizeType calls=0u,evals=0u,changed=0u,unchanged=0u;
        SizeType emptied=0u,without_evaluation=0u;
        for(BoxCallStatistics const& cell:group) {
            calls+=cell.calls;
            evals+=cell.evaluations;
            changed+=cell.changed;
            unchanged+=cell.unchanged;
            emptied+=cell.emptied;
            without_evaluation+=cell.zero_evaluations;
        }
        if(calls!=a.calls or evals!=a.evaluations
           or changed!=a.changed or unchanged!=a.unchanged
           or emptied!=a.emptied
           or without_evaluation!=a.zero_evaluations)
            throw std::runtime_error("inconsistent box group totals");
    };
    check_group(stats.by_output);
    check_group(stats.by_variable);
    check_group(stats.by_pair);
    std::cout<<"TRACE_VERIFIED shape="<<shape
             <<" d=4 f=3 queries="<<stats.queries
             <<" definitely_separated="<<stats.proven_empty
             <<" rounds="<<stats.rounds
             <<" hull_calls="<<stats.hull_calls
             <<" box_calls="<<stats.box_calls
             <<" stop_hull="<<stats.stopped_empty_hull
             <<" stop_box="<<stats.stopped_empty_box
             <<" stop_threshold="<<stats.stopped_threshold<<"\n";
    if(verify_only) return;
    std::cout<<"SOLVER_TRACE shape="<<shape
             <<" d=4 f=3"
             <<" trace_queries="<<stats.queries
             <<" instrumented_total_ms="<<std::fixed<<std::setprecision(4)
                 <<1000.0*stats.total_seconds
             <<" hull_ms="<<1000.0*stats.hull_seconds
             <<" box_ms="<<1000.0*stats.box_seconds
             <<" other_ms="<<1000.0*(
                 stats.total_seconds-stats.hull_seconds-stats.box_seconds)
             <<" paired_replay_wall_ms="<<sweep_ms
             <<"\n";
    std::cout<<"SOLVER_TRACE measures an instrumented legacy full-scan replay and"
                " independent old-reference and production checks per input;"
                " it is NOT a direct timing of paving::separated().\n";

    // The benchmark already verifies final-domain equivalence before
    // reporting the activity counters. Slice evaluation counts are supplied
    // by the production contractor, not estimated from a copy of its loop.
    std::cout<<"BOX_ACTIVITY shape="<<shape
             <<" calls="<<a.calls
             <<" evaluations="<<a.evaluations
             <<" changed="<<a.changed
             <<" unchanged="<<a.unchanged
             <<" proven_empty="<<a.emptied
             <<" zero_evaluations="<<a.zero_evaluations
             <<" evaluations_per_call="<<std::fixed<<std::setprecision(3)
                 <<(a.calls?static_cast<double>(a.evaluations)/a.calls:0.0)
             <<" changed_percent="
                 <<(a.calls?100.0*static_cast<double>(a.changed)/a.calls:0.0)
             <<" unchanged_percent="
                 <<(a.calls?100.0*static_cast<double>(a.unchanged)/a.calls:0.0)
             <<"\n";
    std::cout<<"BOX_ACTIVITY changed means an interval was contracted"
                " or proved empty; unchanged does NOT imply that a"
                " call was safely skippable.\n";

    std::cout<<std::string(118,'=')<<"\n"
             <<std::left<<std::setw(19)<<"Shape"
             <<std::setw(14)<<"Group"
             <<std::right<<std::setw(6)<<"Index"
             <<std::setw(13)<<"Calls"
             <<std::setw(15)<<"Evaluations"
             <<std::setw(13)<<"Changed"
             <<std::setw(13)<<"Unchanged"
             <<std::setw(12)<<"Empty"
             <<std::setw(13)<<"Time (ms)"<<"\n"
             <<std::string(118,'-')<<"\n";
    auto row=[&](char const* group,SizeType index,BoxCallStatistics const& b) {
        std::cout<<std::left<<std::setw(19)<<shape
                 <<std::setw(14)<<group
                 <<std::right<<std::setw(6)<<index
                 <<std::setw(13)<<b.calls
                 <<std::setw(15)<<b.evaluations
                 <<std::setw(13)<<b.changed
                 <<std::setw(13)<<b.unchanged
                 <<std::setw(12)<<b.emptied
                 <<std::fixed<<std::setprecision(3)
                 <<std::setw(13)<<(1000.0*b.seconds)<<"\n";
    };
    for(SizeType i=0u;i<stats.by_output.size();++i) {
        row(i<4u?"image-axis":"constraint",i,stats.by_output[i]);
    }
    std::cout<<std::string(118,'-')<<"\n";
    for(SizeType j=0u;j<stats.by_variable.size();++j) {
        row("variable",j,stats.by_variable[j]);
    }
    std::cout<<std::string(118,'=')<<"\n";

    // Show the hottest output-variable pairs without producing 48 lines
    // per shape. Group totals above remain exhaustive.
    std::vector<SizeType> ranks;
    ranks.reserve(stats.by_pair.size());
    for(SizeType index=0u;index<stats.by_pair.size();++index)
        ranks.push_back(index);
    std::sort(ranks.begin(),ranks.end(),
        [&](SizeType a,SizeType b) {
            return stats.by_pair[a].seconds>stats.by_pair[b].seconds;
        });
    SizeType dimensions=stats.by_variable.size();
    for(SizeType rank=0u;rank<std::min(SizeType(8u),ranks.size());++rank) {
        SizeType index=ranks[rank];
        BoxCallStatistics const& b=stats.by_pair[index];
        std::cout<<"BOX_HOT_PAIR shape="<<shape
                 <<" rank="<<(rank+1u)
                 <<" output="<<(index/dimensions)
                 <<" variable="<<(index%dimensions)
                 <<" calls="<<b.calls
                 <<" evaluations="<<b.evaluations
                 <<" changed="<<b.changed
                 <<" unchanged="<<b.unchanged
                 <<" ms="<<std::fixed<<std::setprecision(4)
                     <<1000.0*b.seconds<<"\n";
    }
}

// Experimental contractor scheduling. A certified image enclosure
// included in the constraint bounds proves that constraint everywhere
// on the current domain. Since reduce only contracts its domain, the
// certificate stays valid for the rest of this query. Neither a static
// shape assumption nor a nonvalidated floating tolerance is used.
struct ContainedSkipStatistics {
    SizeType tests=0u;
    SizeType certificates=0u;
    SizeType reused_certificates=0u;
    SizeType saved_box_calls=0u;
    SizeType executed_box_calls=0u;
};

void reduce_contained_skip(ConstraintSolver const& solver,
                           UpperBoxType& domain,
                           ValidatedVectorMultivariateFunction const& function,
                           ExactBoxType const& codomain,
                           ContainedSkipStatistics* statistics=nullptr) {
    const ExactDouble MINIMUM_REDUCTION=0.75_x;
    if(definitely(domain.is_empty())) return;

    std::vector<Bool> proved(codomain.size(),false);
    FloatDPUpperBound magnitude={0u,dp};
    for(SizeType j=0u;j<domain.size();++j)
        magnitude+=domain[j].width();
    FloatDPUpperBound previous=magnitude;

    do {
        solver.hull_reduce(domain,function,codomain);
        if(definitely(domain.is_empty())) return;

        for(SizeType i=0u;i<codomain.size();++i) {
            if(proved[i]) {
                if(statistics) {
                    ++statistics->reused_certificates;
                    statistics->saved_box_calls+=domain.size();
                }
                continue;
            }
            if(statistics) ++statistics->tests;
            if(definitely(subset(apply(function[i],domain),codomain[i]))) {
                proved[i]=true;
                if(statistics) {
                    ++statistics->certificates;
                    statistics->saved_box_calls+=domain.size();
                }
                continue;
            }
            for(SizeType j=0u;j<domain.size();++j) {
                if(statistics) ++statistics->executed_box_calls;
                solver.box_reduce(domain,function[i],codomain[i],j);
                if(definitely(domain.is_empty())) return;
            }
        }
        previous=magnitude;
        magnitude=0u;
        for(SizeType j=0u;j<domain.size();++j)
            magnitude+=domain[j].width();
    } while(magnitude.raw()
            < mul(near,previous.raw(),FloatDP(MINIMUM_REDUCTION,dp)));
}

template<class T>
void keep_alive(T const& value) {
#if defined(__GNUC__) || defined(__clang__)
    asm volatile("" : : "g"(&value) : "memory");
#else
    static void const* volatile sink=nullptr;
    sink=&value;
#endif
}

struct Timing {
    double milliseconds;
    unsigned iterations;
};

template<class Work>
Timing median_ms(Work&& work,Selection const& selection) {
    using Clock=std::chrono::steady_clock;
    unsigned iterations=1u;
    auto run=[&](unsigned count) {
        auto start=Clock::now();
        for(unsigned i=0u;i<count;++i) {
            auto value=work();
            keep_alive(value);
        }
        return std::chrono::duration<double,std::milli>(
            Clock::now()-start).count();
    };
    double elapsed=run(iterations);
    while(elapsed<static_cast<double>(selection.min_sample_ms)
          and iterations<selection.max_iterations) {
        unsigned next=iterations>selection.max_iterations/2u
            ?selection.max_iterations:iterations*2u;
        if(next==iterations) break;
        iterations=next;
        elapsed=run(iterations);
    }
    std::vector<double> samples;
    samples.reserve(selection.samples);
    for(unsigned i=0u;i<selection.samples;++i) {
        double ms=run(iterations)/iterations;
        if(not std::isfinite(ms) or ms<=0.0)
            throw std::runtime_error("nonfinite/nonpositive timing");
        samples.push_back(ms);
    }
    std::sort(samples.begin(),samples.end());
    return {samples[selection.samples/2u],iterations};
}

// Same immutable adaptive query sequence and contract for both variants.
// Preparation and correctness checks are deliberately outside timed loops.
void compare_contained_skip(Selection const& selection,
                            std::string const& shape) {
    SolverProblem problem=make_problem(shape);
    std::vector<TraceRecord> const trace=
        read_adaptive_trace(selection.trace_prefix,shape);
    std::vector<ExactBoxType> targets;
    targets.reserve(trace.size());
    for(TraceRecord const& entry:trace) {
        ExactBoxType codomain=problem.codomains.front();
        for(SizeType j=0u;j<4u;++j) codomain[j]=entry.box[j];
        targets.push_back(std::move(codomain));
    }
    ConstraintSolver solver;
    ContainedSkipStatistics counters;
    SizeType oracle_empty=0u;
    for(SizeType q=0u;q<targets.size();++q) {
        UpperBoxType baseline=problem.initial_domain;
        reduce_full_scan_reference(solver,baseline,problem.function,targets[q]);
        UpperBoxType candidate=problem.initial_domain;
        reduce_contained_skip(
            solver,candidate,problem.function,targets[q],&counters);
        UpperBoxType production=problem.initial_domain;
        solver.reduce(production,problem.function,targets[q]);
        if(not same(candidate,production))
            throw std::runtime_error(
                "production differs from experimental containment reduction"
                " shape="+shape+" query="+std::to_string(q));
        if(not same(baseline,candidate)) {
            throw std::runtime_error(
                "certified containment pruning changed the contracted domain"
                " shape="+shape+" query="+std::to_string(q));
        }
        Bool empty=definitely(baseline.is_empty());
        if(empty!=trace[q].separated) {
            throw std::runtime_error(
                "baseline separation changed relative to paving trace"
                " shape="+shape+" query="+std::to_string(q));
        }
        oracle_empty+=static_cast<SizeType>(empty);
    }
    std::cout<<"CONTAINED_SKIP_VERIFIED shape="<<shape
             <<" d=4 f=3 queries="<<targets.size()
             <<" separated="<<oracle_empty
             <<" guard_tests="<<counters.tests
             <<" guard_proofs="<<counters.certificates
             <<" reused_proofs="<<counters.reused_certificates
             <<" saved_box_calls="<<counters.saved_box_calls
             <<" executed_box_calls="<<counters.executed_box_calls
             <<"\n";

    if(selection.verify_only) return;
    auto original=[&] {
        SizeType checksum=0u;
        for(ExactBoxType const& target:targets) {
            UpperBoxType domain=problem.initial_domain;
            reduce_full_scan_reference(
                solver,domain,problem.function,target);
            checksum+=static_cast<SizeType>(definitely(domain.is_empty()));
            keep_alive(domain);
        }
        return checksum;
    };
    auto optimized=[&] {
        SizeType checksum=0u;
        for(ExactBoxType const& target:targets) {
            UpperBoxType domain=problem.initial_domain;
            solver.reduce(domain,problem.function,target);
            checksum+=static_cast<SizeType>(definitely(domain.is_empty()));
            keep_alive(domain);
        }
        return checksum;
    };
    if(original()!=oracle_empty or optimized()!=oracle_empty)
        throw std::runtime_error("contained skip sweep checksum mismatch");
    Timing original_time=median_ms(original,selection);
    Timing optimized_time=median_ms(optimized,selection);
    double ratio=original_time.milliseconds/optimized_time.milliseconds;
    if(not std::isfinite(ratio) or ratio<=0.0)
        throw std::runtime_error("invalid contained skip speedup");
    std::cout<<"CONTAINED_SKIP_COMPARE shape="<<shape
             <<" d=4 f=3 queries="<<targets.size()
             <<" baseline_median_ms="<<std::fixed<<std::setprecision(5)
                <<original_time.milliseconds
             <<" guarded_median_ms="<<optimized_time.milliseconds
             <<" speedup="<<std::setprecision(3)<<ratio
             <<" baseline_iters="<<original_time.iterations
             <<" guarded_iters="<<optimized_time.iterations
             <<" samples="<<selection.samples
             <<" min_sample_ms="<<selection.min_sample_ms
             <<"\n";
    std::cout<<"CONTAINED_SKIP_COMPARE uses independent calibrated medians"
                " and requires exact final-domain and separation equality"
                " for every adaptive query. Production compared to frozen full scan.\n";
}

// Experimental P2: a Procedure makes the syntactic input-coordinate
// dependency set explicit via IndexProcedureInstruction. This is a
// conservative mask: unused coordinates are absent even when the
// non-linear interval extension is inconclusive about containment.
using DependencyMasks=std::vector<std::vector<unsigned char>>;

DependencyMasks compile_dependency_masks(
        ValidatedVectorMultivariateFunction const& function) {
    DependencyMasks masks(
        function.result_size(),
        std::vector<unsigned char>(function.argument_size(),0u));
    for(SizeType i=0u;i<function.result_size();++i) {
        ValidatedProcedure procedure(function[i]);
        for(auto const& instruction:procedure._instructions) {
            auto const* index=std::get_if<IndexProcedureInstruction>(
                &instruction.base());
            if(index==nullptr) continue;
            if(index->_ind>=function.argument_size())
                throw std::runtime_error("invalid procedure dependency index");
            masks[i][index->_ind]=1u;
        }
    }
    return masks;
}

struct DependencySkipStatistics {
    SizeType guard_tests=0u;
    SizeType certified_outputs=0u;
    SizeType reused_certificates=0u;
    SizeType skipped_pairs=0u;
    SizeType executed_box_calls=0u;
    SizeType disjoint_fallback_outputs=0u;
};

void reduce_dependency_skip(
        ConstraintSolver const& solver,UpperBoxType& domain,
        ValidatedVectorMultivariateFunction const& function,
        ExactBoxType const& codomain,DependencyMasks const& masks,
        DependencySkipStatistics* stats=nullptr) {
    const ExactDouble MINIMUM_REDUCTION=0.75_x;
    if(masks.size()!=codomain.size())
        throw std::invalid_argument("dependency output dimension mismatch");
    for(auto const& mask:masks) {
        if(mask.size()!=domain.size())
            throw std::invalid_argument("dependency input dimension mismatch");
    }
    if(definitely(domain.is_empty())) return;

    std::vector<Bool> certified_satisfied(codomain.size(),false);
    FloatDPUpperBound magnitude={0u,dp};
    for(SizeType j=0u;j<domain.size();++j)
        magnitude+=domain[j].width();
    FloatDPUpperBound previous=magnitude;
    do {
        solver.hull_reduce(domain,function,codomain);
        if(definitely(domain.is_empty())) return;
        for(SizeType i=0u;i<codomain.size();++i) {
            if(certified_satisfied[i]) {
                if(stats) ++stats->reused_certificates;
                continue;
            }

            if(stats) ++stats->guard_tests;
            UpperIntervalType image=Ariadne::apply(function[i],domain);
            if(definitely(subset(image,codomain[i]))) {
                certified_satisfied[i]=true;
                if(stats) ++stats->certified_outputs;
                continue;
            }

            // If direct interval evaluation already *proves* exclusion,
            // retain the original complete sweep. Otherwise a contractor
            // on an unused coordinate has no useful shaving effect.
            // This conservative fallback avoids losing an emptiness
            // detection when the hull and direct interval evaluations
            // happen to have different pruning strength.
            Bool const any_independent=
                std::find(masks[i].begin(),masks[i].end(),0u)
                !=masks[i].end();
            Bool const disjoint_fallback=any_independent
                and definitely(disjoint(image,codomain[i]));
            if(disjoint_fallback and stats)
                ++stats->disjoint_fallback_outputs;
            for(SizeType j=0u;j<domain.size();++j) {
                if(not disjoint_fallback and masks[i][j]==0u) {
                    if(stats) ++stats->skipped_pairs;
                    continue;
                }
                if(stats) ++stats->executed_box_calls;
                solver.box_reduce(domain,function[i],codomain[i],j);
                if(definitely(domain.is_empty())) return;
            }
        }
        previous=magnitude;
        magnitude=0u;
        for(SizeType j=0u;j<domain.size();++j)
            magnitude+=domain[j].width();
    } while(magnitude.raw()
            < mul(near,previous.raw(),FloatDP(MINIMUM_REDUCTION,dp)));
}

// Compare the production contained-skip reducer against an additional
// dependency-aware candidate. Compiling the masks is measured both as
// a one-time per-function setup and *inside* the per-query hot path;
// the latter is the realistic cost of using this optimization directly
// in the current reduce() API, which receives uncompiled functions.
void compare_dependency_skip(Selection const& selection,
                             std::string const& shape) {
    SolverProblem problem=make_problem(shape);
    std::vector<TraceRecord> const trace=
        read_adaptive_trace(selection.trace_prefix,shape);
    std::vector<ExactBoxType> targets;
    targets.reserve(trace.size());
    for(TraceRecord const& entry:trace) {
        ExactBoxType codomain=problem.codomains.front();
        for(SizeType j=0u;j<4u;++j) codomain[j]=entry.box[j];
        targets.push_back(std::move(codomain));
    }

    ConstraintSolver solver;
    DependencyMasks const prepared=compile_dependency_masks(problem.function);
    SizeType total_edges=0u;
    for(auto const& mask:prepared) {
        for(unsigned char bit:mask)
            total_edges+=static_cast<SizeType>(bit!=0u);
    }
    DependencySkipStatistics stats;
    SizeType separated=0u;
    for(SizeType q=0u;q<targets.size();++q) {
        UpperBoxType baseline=problem.initial_domain;
        solver.reduce(baseline,problem.function,targets[q]);
        UpperBoxType candidate=problem.initial_domain;
        reduce_dependency_skip(
            solver,candidate,problem.function,targets[q],prepared,&stats);
        if(not same(candidate,baseline)) {
            throw std::runtime_error(
                "dependency pruning changed production contracted domain"
                " shape="+shape+" query="+std::to_string(q));
        }
        UpperBoxType legacy=problem.initial_domain;
        reduce_full_scan_reference(solver,legacy,problem.function,targets[q]);
        if(not same(legacy,baseline)) {
            throw std::runtime_error(
                "production differs from legacy full-scan oracle"
                " shape="+shape+" query="+std::to_string(q));
        }
        Bool const empty=definitely(baseline.is_empty());
        if(empty!=trace[q].separated) {
            throw std::runtime_error(
                "replayed separation does not match exported paving trace"
                " shape="+shape+" query="+std::to_string(q));
        }
        separated+=static_cast<SizeType>(empty);
    }
    std::cout<<"DEPENDENCY_SKIP_VERIFIED shape="<<shape
             <<" d=4 f=3 queries="<<targets.size()
             <<" separated="<<separated
             <<" edges="<<total_edges
             <<" possible_edges="<<(prepared.size()*problem.initial_domain.size())
             <<" guard_tests="<<stats.guard_tests
             <<" guard_proofs="<<stats.certified_outputs
             <<" reused_proofs="<<stats.reused_certificates
             <<" dependency_skipped="<<stats.skipped_pairs
             <<" executed_box_calls="<<stats.executed_box_calls
             <<" disjoint_fallbacks="<<stats.disjoint_fallback_outputs
             <<"\n";
    if(selection.verify_only) return;

    auto original=[&]() {
        SizeType checksum=0u;
        for(ExactBoxType const& codomain:targets) {
            UpperBoxType domain=problem.initial_domain;
            solver.reduce(domain,problem.function,codomain);
            checksum+=static_cast<SizeType>(definitely(domain.is_empty()));
            keep_alive(domain);
        }
        return checksum;
    };
    auto prepared_work=[&]() {
        SizeType checksum=0u;
        for(ExactBoxType const& codomain:targets) {
            UpperBoxType domain=problem.initial_domain;
            reduce_dependency_skip(
                solver,domain,problem.function,codomain,prepared);
            checksum+=static_cast<SizeType>(definitely(domain.is_empty()));
            keep_alive(domain);
        }
        return checksum;
    };
    auto rebuilt_work=[&]() {
        SizeType checksum=0u;
        for(ExactBoxType const& codomain:targets) {
            DependencyMasks masks=compile_dependency_masks(problem.function);
            UpperBoxType domain=problem.initial_domain;
            reduce_dependency_skip(
                solver,domain,problem.function,codomain,masks);
            checksum+=static_cast<SizeType>(definitely(domain.is_empty()));
            keep_alive(domain);
        }
        return checksum;
    };

    if(original()!=separated or prepared_work()!=separated
       or rebuilt_work()!=separated) {
        throw std::runtime_error(
            "dependency pruning sweep checksum mismatch for "+shape);
    }
    Timing production=median_ms(original,selection);
    Timing prepared_time=median_ms(prepared_work,selection);
    Timing rebuilt_time=median_ms(rebuilt_work,selection);
    double prepared_ratio=production.milliseconds/prepared_time.milliseconds;
    double rebuilt_ratio=production.milliseconds/rebuilt_time.milliseconds;
    if(not std::isfinite(prepared_ratio) or not std::isfinite(rebuilt_ratio)
       or prepared_ratio<=0.0 or rebuilt_ratio<=0.0)
        throw std::runtime_error("invalid dependency pruning speed ratio");

    std::cout<<"DEPENDENCY_SKIP_COMPARE shape="<<shape
             <<" d=4 f=3 queries="<<targets.size()
             <<" production_median_ms="<<std::fixed<<std::setprecision(5)
                <<production.milliseconds
             <<" prepared_median_ms="<<prepared_time.milliseconds
             <<" rebuilt_median_ms="<<rebuilt_time.milliseconds
             <<" prepared_speedup="<<std::setprecision(3)<<prepared_ratio
             <<" rebuilt_speedup="<<rebuilt_ratio
             <<" production_iters="<<production.iterations
             <<" prepared_iters="<<prepared_time.iterations
             <<" rebuilt_iters="<<rebuilt_time.iterations
             <<" samples="<<selection.samples
             <<"\n";
    std::cout<<"DEPENDENCY_SKIP_COMPARE: production is the certified-"
                "containment reducer; prepared excludes dependency compilation"
                " from the timed sweep, rebuilt includes it for every query."
                " All three produce identical contracted domains and"
                " separation results; no production code changed.\n";
}

// Instrument the *current production P1 schedule*, not the historical
// full-scan loop. The 5-argument box_reduce overload runs exactly the
// same contractor, but reports the number of validated interval-function
// evaluations. Nothing in this profiler changes the production solver.
struct ResidualShavingStatistics {
    SizeType queries=0u;
    SizeType separated=0u;
    SizeType rounds=0u;
    SizeType hull_calls=0u;
    SizeType stopped_empty_hull=0u;
    SizeType stopped_empty_box=0u;
    SizeType stopped_threshold=0u;
    SizeType guard_tests=0u;
    SizeType guard_proofs=0u;
    SizeType reused_proofs=0u;
    SizeType unproved_outputs=0u;
    SizeType unchanged_evaluations=0u;
    SizeType changed_evaluations=0u;
    double total_seconds=0.0;
    double hull_seconds=0.0;
    double guard_seconds=0.0;
    double box_seconds=0.0;
    BoxCallStatistics all;
    std::vector<BoxCallStatistics> by_output;
    std::vector<BoxCallStatistics> by_variable;
    std::vector<BoxCallStatistics> by_pair;
    std::array<SizeType,9u> histogram{};
    std::array<SizeType,9u> unchanged_histogram{};

    void prepare(SizeType outputs,SizeType variables) {
        by_output.resize(outputs);
        by_variable.resize(variables);
        by_pair.resize(outputs*variables);
    }

    void record_box(SizeType i,SizeType j,SizeType variables,
                    SizeType evaluations,Bool changed,Bool empty,
                    double seconds) {
        if(evaluations>=histogram.size())
            throw std::runtime_error("box_reduce used more than eight slices");
        all.add(evaluations,changed,empty,seconds);
        by_output[i].add(evaluations,changed,empty,seconds);
        by_variable[j].add(evaluations,changed,empty,seconds);
        by_pair[i*variables+j].add(evaluations,changed,empty,seconds);
        ++histogram[evaluations];
        if(changed) {
            changed_evaluations+=evaluations;
        } else {
            unchanged_evaluations+=evaluations;
            ++unchanged_histogram[evaluations];
        }
    }

    void validate() const {
        if(queries!=stopped_empty_hull+stopped_empty_box+stopped_threshold
           or queries!=separated+stopped_threshold
           or all.changed+all.unchanged!=all.calls
           or all.zero_evaluations!=histogram[0]
           or all.evaluations!=changed_evaluations+unchanged_evaluations)
            throw std::runtime_error("inconsistent residual summary");
        SizeType histogram_calls=0u, histogram_evaluations=0u;
        SizeType unchanged_histogram_calls=0u;
        for(SizeType n=0u;n<histogram.size();++n) {
            histogram_calls+=histogram[n];
            histogram_evaluations+=n*histogram[n];
            unchanged_histogram_calls+=unchanged_histogram[n];
            if(unchanged_histogram[n]>histogram[n])
                throw std::runtime_error("invalid unchanged residual histogram");
        }
        if(histogram_calls!=all.calls
           or histogram_evaluations!=all.evaluations
           or unchanged_histogram_calls!=all.unchanged)
            throw std::runtime_error("invalid residual evaluation histogram");
        for(auto const* group:{&by_output,&by_variable,&by_pair}) {
            SizeType calls=0u,evals=0u,changed=0u,unchanged=0u;
            SizeType emptied=0u,zero=0u;
            for(auto const& r:*group) {
                calls+=r.calls;
                evals+=r.evaluations;
                changed+=r.changed;
                unchanged+=r.unchanged;
                emptied+=r.emptied;
                zero+=r.zero_evaluations;
            }
            if(calls!=all.calls or evals!=all.evaluations
               or changed!=all.changed or unchanged!=all.unchanged
               or emptied!=all.emptied or zero!=all.zero_evaluations)
                throw std::runtime_error("invalid residual group statistics");
        }
    }
};

void traced_reduce_residual(
        ConstraintSolver const& solver,UpperBoxType& domain,
        ValidatedVectorMultivariateFunction const& function,
        ExactBoxType const& codomain,ResidualShavingStatistics& stats) {
    using Clock=std::chrono::steady_clock;
    auto elapsed=[](Clock::time_point begin) {
        return std::chrono::duration<double>(
            Clock::now()-begin).count();
    };
    auto const begin=Clock::now();
    ++stats.queries;
    if(definitely(domain.is_empty())) {
        ++stats.separated;
        ++stats.stopped_empty_hull;
        stats.total_seconds+=elapsed(begin);
        return;
    }

    const ExactDouble MINIMUM_REDUCTION=0.75_x;
    std::vector<Bool> certified_satisfied(codomain.size(),false);
    FloatDPUpperBound magnitude={0u,dp};
    for(SizeType j=0u;j<domain.size();++j)
        magnitude+=domain[j].width();
    FloatDPUpperBound previous=magnitude;

    do {
        ++stats.rounds;
        ++stats.hull_calls;
        auto const hstart=Clock::now();
        solver.hull_reduce(domain,function,codomain);
        stats.hull_seconds+=elapsed(hstart);
        if(definitely(domain.is_empty())) {
            ++stats.separated;
            ++stats.stopped_empty_hull;
            stats.total_seconds+=elapsed(begin);
            return;
        }

        for(SizeType i=0u;i<codomain.size();++i) {
            if(certified_satisfied[i]) {
                ++stats.reused_proofs;
                continue;
            }
            ++stats.guard_tests;
            auto const gstart=Clock::now();
            Bool const proved=definitely(subset(
                Ariadne::apply(function[i],domain),codomain[i]));
            stats.guard_seconds+=elapsed(gstart);
            if(proved) {
                certified_satisfied[i]=true;
                ++stats.guard_proofs;
                continue;
            }

            ++stats.unproved_outputs;
            for(SizeType j=0u;j<domain.size();++j) {
                UpperIntervalType const before=domain[j];
                SizeType evaluations=0u;
                auto const bstart=Clock::now();
                solver.box_reduce(
                    domain,function[i],codomain[i],j,evaluations);
                double const seconds=elapsed(bstart);
                stats.box_seconds+=seconds;
                Bool const empty=definitely(domain.is_empty());
                stats.record_box(i,j,domain.size(),evaluations,
                                 not same(before,domain[j]),empty,seconds);
                if(empty) {
                    ++stats.separated;
                    ++stats.stopped_empty_box;
                    stats.total_seconds+=elapsed(begin);
                    return;
                }
            }
        }
        previous=magnitude;
        magnitude=0u;
        for(SizeType j=0u;j<domain.size();++j)
            magnitude+=domain[j].width();
    } while(magnitude.raw()
            < mul(near,previous.raw(),FloatDP(MINIMUM_REDUCTION,dp)));
    ++stats.stopped_threshold;
    stats.total_seconds+=elapsed(begin);
}

void profile_residual_shaving(Selection const& selection,
                              std::string const& shape) {
    SolverProblem problem=make_problem(shape);
    auto records=read_adaptive_trace(selection.trace_prefix,shape);
    ConstraintSolver solver;
    ResidualShavingStatistics stats;
    stats.prepare(problem.function.result_size(),problem.initial_domain.size());
    for(SizeType q=0u;q<records.size();++q) {
        auto const& record=records[q];
        ExactBoxType codomain=problem.codomains.front();
        for(SizeType j=0u;j<4u;++j) codomain[j]=record.box[j];

        UpperBoxType production=problem.initial_domain;
        solver.reduce(production,problem.function,codomain);
        UpperBoxType traced=problem.initial_domain;
        traced_reduce_residual(solver,traced,problem.function,codomain,stats);
        if(not same(production,traced))
            throw std::runtime_error(
                "residual shaving instrumentation differs from production"
                " shape="+shape+" query="+std::to_string(q));

        UpperBoxType legacy=problem.initial_domain;
        reduce_full_scan_reference(solver,legacy,problem.function,codomain);
        if(not same(legacy,traced))
            throw std::runtime_error(
                "residual shaving differs from legacy full scan"
                " shape="+shape+" query="+std::to_string(q));
        if(definitely(traced.is_empty())!=record.separated)
            throw std::runtime_error(
                "residual shaving disagrees with paving separation trace"
                " shape="+shape+" query="+std::to_string(q));
    }
    stats.validate();

    std::cout<<"RESIDUAL_SHAVING_VERIFIED shape="<<shape
             <<" d=4 f=3 queries="<<stats.queries
             <<" separated="<<stats.separated
             <<" rounds="<<stats.rounds
             <<" hull_calls="<<stats.hull_calls
             <<" guard_tests="<<stats.guard_tests
             <<" guard_proofs="<<stats.guard_proofs
             <<" reused_proofs="<<stats.reused_proofs
             <<" unproved_outputs="<<stats.unproved_outputs
             <<" box_calls="<<stats.all.calls
             <<" box_evaluations="<<stats.all.evaluations
             <<" box_changed="<<stats.all.changed
             <<" box_unchanged="<<stats.all.unchanged
             <<" box_empty="<<stats.all.emptied
             <<"\n";
    if(selection.verify_only) return;

    std::cout<<"RESIDUAL_SHAVING_TIMES shape="<<shape
             <<" instrumented_total_ms="<<std::fixed<<std::setprecision(4)
                <<1000.0*stats.total_seconds
             <<" hull_ms="<<1000.0*stats.hull_seconds
             <<" guard_ms="<<1000.0*stats.guard_seconds
             <<" box_ms="<<1000.0*stats.box_seconds
             <<" other_ms="<<1000.0*(
                stats.total_seconds-stats.hull_seconds
                -stats.guard_seconds-stats.box_seconds)
             <<" unchanged_box_evaluations="<<stats.unchanged_evaluations
             <<" changed_box_evaluations="<<stats.changed_evaluations
             <<"\n";

    std::cout<<"RESIDUAL_SHAVING_EVAL_HIST shape="<<shape;
    for(SizeType n=0u;n<stats.histogram.size();++n)
        std::cout<<" eval_"<<n<<"="<<stats.histogram[n];
    std::cout<<"\n";

    std::cout<<"RESIDUAL_SHAVING_UNCHANGED_HIST shape="<<shape;
    for(SizeType n=0u;n<stats.unchanged_histogram.size();++n)
        std::cout<<" eval_"<<n<<"="<<stats.unchanged_histogram[n];
    std::cout<<"\n";

    auto print_group=[&](char const* name,SizeType index,
                         BoxCallStatistics const& s) {
        std::cout<<"RESIDUAL_SHAVING_GROUP shape="<<shape
                 <<" group="<<name<<" index="<<index
                 <<" calls="<<s.calls
                 <<" evaluations="<<s.evaluations
                 <<" changed="<<s.changed
                 <<" unchanged="<<s.unchanged
                 <<" empty="<<s.emptied
                 <<" zero_evaluations="<<s.zero_evaluations
                 <<" elapsed_ms="<<std::fixed<<std::setprecision(4)
                    <<1000.0*s.seconds<<"\n";
    };
    for(SizeType i=0u;i<stats.by_output.size();++i) {
        if(stats.by_output[i].calls!=0u)
            print_group("output",i,stats.by_output[i]);
    }
    for(SizeType j=0u;j<stats.by_variable.size();++j)
        print_group("variable",j,stats.by_variable[j]);

    // Rank output/variable pairs by measured box-contractor time.
    // This is an instrumented single replay, *not* a calibrated speedup.
    std::vector<SizeType> ranking;
    ranking.reserve(stats.by_pair.size());
    for(SizeType pair=0u;pair<stats.by_pair.size();++pair) {
        if(stats.by_pair[pair].calls!=0u) ranking.push_back(pair);
    }
    std::sort(ranking.begin(),ranking.end(),[&](SizeType a,SizeType b) {
        if(stats.by_pair[a].seconds!=stats.by_pair[b].seconds)
            return stats.by_pair[a].seconds>stats.by_pair[b].seconds;
        return a<b;
    });
    SizeType const variables=problem.initial_domain.size();
    SizeType const limit=std::min(ranking.size(),SizeType(12u));
    for(SizeType rank=0u;rank<limit;++rank) {
        SizeType const pair=ranking[rank];
        auto const& s=stats.by_pair[pair];
        std::cout<<"RESIDUAL_SHAVING_HOT_PAIR shape="<<shape
                 <<" rank="<<(rank+1u)
                 <<" output="<<(pair/variables)
                 <<" variable="<<(pair%variables)
                 <<" calls="<<s.calls
                 <<" evaluations="<<s.evaluations
                 <<" changed="<<s.changed
                 <<" unchanged="<<s.unchanged
                 <<" empty="<<s.emptied
                 <<" zero_evaluations="<<s.zero_evaluations
                 <<" elapsed_ms="<<std::fixed<<std::setprecision(4)
                    <<1000.0*s.seconds
                 <<"\n";
    }
    std::cout<<"RESIDUAL_SHAVING_NOTE: instrumented P1 replay only;"
                " timings include profiling overhead. Histograms count"
                " actual interval-function evaluations from box_reduce;"
                " unchanged does not prove a call was dispensable."
                " Full domain equality with production and the frozen"
                " legacy reference was checked for every query.\n";
}

// P5, benchmark-only: a sound boundary contractor using a dyadic
// subdivision tree rather than the production eight equal slices.
// Every removed prefix/suffix consists entirely of cells for which a
// validated interval image is certified disjoint from the codomain.
// An inconclusive leaf is NEVER removed; bisection only refines it.
struct AdaptiveShavingSettings {
    unsigned max_depth;
    unsigned evaluation_budget;
};

struct AdaptiveShavingStatistics {
    SizeType calls=0u;
    SizeType evaluations=0u;
    SizeType changed_calls=0u;
    SizeType emptied_calls=0u;
    SizeType rejected_cells=0u;
    SizeType unresolved_leaves=0u;
    SizeType budget_stops=0u;
    SizeType depth_stops=0u;
    SizeType certified_outputs=0u;
    SizeType rounds=0u;
    SizeType cell_visits=0u;
    SizeType cache_hits=0u;
};

// The production solver's affine() interpolation helper is local to
// src/constraint_solver.cpp and is not visible to this benchmark.
// Keep the same arithmetic here; do not add a shared benchmark header.
template<class X>
inline Approximation<X> adaptive_affine(
        Approximation<X> l,Approximation<X> u,Nat i,Nat n) {
    return (l*(n-i)+u*i)/n;
}

Bool adaptive_box_reduce(UpperBoxType& domain,
                         ValidatedScalarMultivariateFunction const& function,
                         ExactIntervalType const& bounds,SizeType variable,
                         AdaptiveShavingSettings const& settings,
                         AdaptiveShavingStatistics* stats=nullptr,
                         bool reuse_evaluations=false) {
    if(stats) ++stats->calls;
    if(definitely(domain.is_empty())) return true;
    if(definitely(domain[variable].lower_bound()
                  >=domain[variable].upper_bound())) return false;

    UpperIntervalType const original=domain[variable];
    FloatDPApproximation const lower=original.lower_bound();
    FloatDPApproximation const upper=original.upper_bound();
    UpperIntervalType narrowed=original;
    struct Cell { unsigned depth; unsigned index; };
    struct CachedProof { unsigned depth; unsigned index; bool disjoint; };
    // Per-box_reduce cache only: the other input coordinates are fixed
    // throughout these two boundary scans, but change in later calls.
    std::vector<CachedProof> cache;
    if(reuse_evaluations)
        cache.reserve(std::min(settings.evaluation_budget,64u));
    bool entire_domain_excluded=false;

    auto scan=[&](bool from_lower,unsigned budget)->bool {
        std::vector<Cell> todo({Cell{0u,0u}});
        bool unresolved=false;
        unsigned used=0u;
        while(not todo.empty() and used<budget) {
            Cell const cell=todo.back();
            todo.pop_back();
            Nat const parts=static_cast<Nat>(1u<<cell.depth);
            Nat const index=static_cast<Nat>(cell.index);
            ExactIntervalType const interval(
                cast_exact(adaptive_affine(lower,upper,index,parts)),
                cast_exact(adaptive_affine(lower,upper,index+1u,parts)));
            ++used;
            if(stats) ++stats->cell_visits;
            bool disjoint_proof=false;
            bool cache_hit=false;
            if(reuse_evaluations) {
                for(CachedProof const& proof:cache) {
                    if(proof.depth==cell.depth and proof.index==cell.index) {
                        disjoint_proof=proof.disjoint;
                        cache_hit=true;
                        break;
                    }
                }
            }
            if(cache_hit) {
                if(stats) ++stats->cache_hits;
            } else {
                UpperBoxType slice=domain;
                slice[variable]=interval;
                if(stats) ++stats->evaluations;
                disjoint_proof=definitely(disjoint(
                    Ariadne::apply(function,slice),bounds));
                if(reuse_evaluations)
                    cache.push_back({cell.depth,cell.index,disjoint_proof});
            }
            if(disjoint_proof) {
                if(stats) ++stats->rejected_cells;
                if(from_lower) {
                    narrowed.set_lower_bound(interval.upper_bound());
                } else {
                    narrowed.set_upper_bound(interval.lower_bound());
                }
                continue;
            }
            if(cell.depth<settings.max_depth) {
                unsigned const child_depth=cell.depth+1u;
                unsigned const first=cell.index*2u;
                // The last appended cell is processed first, preserving
                // contiguous left-to-right or right-to-left exploration.
                if(from_lower) {
                    todo.push_back({child_depth,first+1u});
                    todo.push_back({child_depth,first});
                } else {
                    todo.push_back({child_depth,first});
                    todo.push_back({child_depth,first+1u});
                }
            } else {
                unresolved=true;
                if(stats) {
                    ++stats->unresolved_leaves;
                    ++stats->depth_stops;
                }
                break;
            }
        }
        if(not todo.empty() and used==budget and stats)
            ++stats->budget_stops;
        // No unresolved leaf and nothing left on the stack means
        // the entire original interval has been certified excluded.
        return not unresolved and todo.empty();
    };

    // A symmetric budget ensures neither boundary starves the other.
    unsigned const left_budget=(settings.evaluation_budget+1u)/2u;
    unsigned const right_budget=settings.evaluation_budget-left_budget;
    entire_domain_excluded=scan(true,left_budget);
    if(not entire_domain_excluded)
        entire_domain_excluded=scan(false,right_budget);

    if(entire_domain_excluded or definitely(narrowed.is_empty())) {
        domain[variable]=ExactIntervalType(+inf,-inf);
        if(stats) {
            ++stats->changed_calls;
            ++stats->emptied_calls;
        }
        return true;
    }
    domain[variable]=narrowed;
    if(stats and not same(original,narrowed))
        ++stats->changed_calls;
    return false;
}

Bool reduce_adaptive_shaving(
        ConstraintSolver const& solver,UpperBoxType& domain,
        ValidatedVectorMultivariateFunction const& function,
        ExactBoxType const& codomain,
        AdaptiveShavingSettings const& settings,
        AdaptiveShavingStatistics* stats=nullptr,
        bool reuse_evaluations=false) {
    const ExactDouble MINIMUM_REDUCTION=0.75_x;
    if(definitely(domain.is_empty())) return true;
    std::vector<Bool> certified(codomain.size(),false);
    FloatDPUpperBound magnitude={0u,dp};
    for(SizeType j=0u;j<domain.size();++j)
        magnitude+=domain[j].width();
    FloatDPUpperBound previous=magnitude;

    do {
        if(stats) ++stats->rounds;
        solver.hull_reduce(domain,function,codomain);
        if(definitely(domain.is_empty())) return true;
        for(SizeType i=0u;i<codomain.size();++i) {
            if(certified[i]) continue;
            if(definitely(subset(
                    Ariadne::apply(function[i],domain),codomain[i]))) {
                certified[i]=true;
                if(stats) ++stats->certified_outputs;
                continue;
            }
            for(SizeType j=0u;j<domain.size();++j) {
                if(adaptive_box_reduce(
                        domain,function[i],codomain[i],j,settings,stats,
                        reuse_evaluations))
                    return true;
            }
        }
        previous=magnitude;
        magnitude=0u;
        for(SizeType j=0u;j<domain.size();++j)
            magnitude+=domain[j].width();
    } while(magnitude.raw()
            < mul(near,previous.raw(),FloatDP(MINIMUM_REDUCTION,dp)));
    return false;
}

struct AdaptiveComparison {
    SizeType exactly_same=0u;
    SizeType stronger=0u;
    SizeType weaker=0u;
    SizeType incomparable=0u;
    SizeType candidate_only_empty=0u;
    SizeType production_only_empty=0u;
    SizeType both_empty=0u;
    SizeType feasible_witnesses=0u;
    AdaptiveShavingStatistics activity;
};

void compare_adaptive_shaving(Selection const& selection,
                              std::string const& shape) {
    SolverProblem problem=make_problem(shape);
    auto const records=read_adaptive_trace(selection.trace_prefix,shape);
    std::vector<ExactBoxType> codomains;
    codomains.reserve(records.size());
    for(auto const& record:records) {
        ExactBoxType codomain=problem.codomains.front();
        for(SizeType i=0u;i<4u;++i) codomain[i]=record.box[i];
        codomains.push_back(std::move(codomain));
    }

    ConstraintSolver solver;
    AdaptiveShavingSettings const settings[]={
        {selection.adaptive_depth,selection.adaptive_budget},
        {selection.adaptive_depth,selection.adaptive_resolution_budget}
    };
    char const* labels[]={"budget","resolution"};
    AdaptiveComparison comparisons[2];
    SizeType production_separated=0u;
    for(SizeType q=0u;q<records.size();++q) {
        UpperBoxType production=problem.initial_domain;
        solver.reduce(production,problem.function,codomains[q]);
        UpperBoxType legacy=problem.initial_domain;
        reduce_full_scan_reference(solver,legacy,problem.function,codomains[q]);
        if(not same(production,legacy))
            throw std::runtime_error(
                "production differs from frozen full-scan"
                " shape="+shape+" query="+std::to_string(q));
        Bool const empty=definitely(production.is_empty());
        if(empty!=records[q].separated)
            throw std::runtime_error(
                "baseline differs from exported paving separation"
                " shape="+shape+" query="+std::to_string(q));
        production_separated+=static_cast<SizeType>(empty);

        // The centre is a known feasible point of all five geometries,
        // but is a witness for this *query* only if validated evaluation
        // proves its entire image lies within the query codomain.
        UpperBoxType const witness=ExactBoxType({
            {0.5_x,0.5_x},{0.5_x,0.5_x},
            {0.5_x,0.5_x},{0.5_x,0.5_x}
        });
        bool verified_witness=true;
        for(SizeType i=0u;i<codomains[q].size();++i) {
            if(not definitely(subset(
                    Ariadne::apply(problem.function[i],witness),
                    codomains[q][i]))) {
                verified_witness=false;
                break;
            }
        }

        for(unsigned method=0u;method!=2u;++method) {
            AdaptiveComparison& result=comparisons[method];
            UpperBoxType candidate=problem.initial_domain;
            reduce_adaptive_shaving(
                solver,candidate,problem.function,codomains[q],
                settings[method],&result.activity);
            if(not definitely(refines(candidate,problem.initial_domain)))
                throw std::runtime_error(
                    "adaptive contractor enlarged its initial domain"
                    " shape="+shape+" method="+labels[method]
                    +" query="+std::to_string(q));
            if(verified_witness) {
                ++result.feasible_witnesses;
                if(definitely(disjoint(candidate,witness)))
                    throw std::runtime_error(
                        "adaptive contractor lost a certified feasible witness"
                        " shape="+shape+" method="+labels[method]
                        +" query="+std::to_string(q));
            }
            Bool const candidate_empty=definitely(candidate.is_empty());
            result.candidate_only_empty+=
                static_cast<SizeType>(candidate_empty and not empty);
            result.production_only_empty+=
                static_cast<SizeType>(not candidate_empty and empty);
            result.both_empty+=
                static_cast<SizeType>(candidate_empty and empty);
            if(same(candidate,production)) {
                ++result.exactly_same;
            } else if(definitely(refines(candidate,production))) {
                ++result.stronger;
            } else if(definitely(refines(production,candidate))) {
                ++result.weaker;
            } else {
                ++result.incomparable;
            }
        }
    }

    for(unsigned method=0u;method!=2u;++method) {
        auto const& result=comparisons[method];
        auto const& stats=result.activity;
        if(result.exactly_same+result.stronger+result.weaker
           +result.incomparable!=records.size())
            throw std::runtime_error("inconsistent adaptive comparison counters");
        if(stats.evaluations>stats.calls*
                static_cast<SizeType>(settings[method].evaluation_budget))
            throw std::runtime_error("adaptive shaving exceeded per-call budget");
        std::cout<<"ADAPTIVE_SHAVING_VERIFIED shape="<<shape
                 <<" method="<<labels[method]
                 <<" queries="<<records.size()
                 <<" baseline_separated="<<production_separated
                 <<" exactly_same="<<result.exactly_same
                 <<" stronger="<<result.stronger
                 <<" weaker="<<result.weaker
                 <<" incomparable="<<result.incomparable
                 <<" candidate_only_empty="<<result.candidate_only_empty
                 <<" production_only_empty="<<result.production_only_empty
                 <<" both_empty="<<result.both_empty
                 <<" certified_witnesses="<<result.feasible_witnesses
                 <<" box_calls="<<stats.calls
                 <<" box_evaluations="<<stats.evaluations
                 <<" changed_box_calls="<<stats.changed_calls
                 <<" excluded_cells="<<stats.rejected_cells
                 <<" depth_stops="<<stats.depth_stops
                 <<" budget_stops="<<stats.budget_stops
                 <<" max_depth="<<settings[method].max_depth
                 <<" evaluation_budget="<<settings[method].evaluation_budget
                 <<"\n";
    }
    if(selection.verify_only) return;

    auto production_work=[&] {
        SizeType result=0u;
        for(auto const& codomain:codomains) {
            UpperBoxType domain=problem.initial_domain;
            solver.reduce(domain,problem.function,codomain);
            result+=static_cast<SizeType>(definitely(domain.is_empty()));
            keep_alive(domain);
        }
        return result;
    };
    Timing const baseline=median_ms(production_work,selection);
    for(unsigned method=0u;method!=2u;++method) {
        auto adaptive_work=[&] {
            SizeType result=0u;
            for(auto const& codomain:codomains) {
                UpperBoxType domain=problem.initial_domain;
                reduce_adaptive_shaving(
                    solver,domain,problem.function,codomain,settings[method]);
                result+=static_cast<SizeType>(definitely(domain.is_empty()));
                keep_alive(domain);
            }
            return result;
        };
        Timing const adaptive=median_ms(adaptive_work,selection);
        double const ratio=baseline.milliseconds/adaptive.milliseconds;
        if(not std::isfinite(ratio) or ratio<=0.0)
            throw std::runtime_error("invalid adaptive shaving speed ratio");
        std::cout<<"ADAPTIVE_SHAVING_COMPARE shape="<<shape
                 <<" method="<<labels[method]
                 <<" production_median_ms="<<std::fixed<<std::setprecision(5)
                    <<baseline.milliseconds
                 <<" adaptive_median_ms="<<adaptive.milliseconds
                 <<" speedup="<<std::setprecision(3)<<ratio
                 <<" samples="<<selection.samples
                 <<" production_iters="<<baseline.iterations
                 <<" adaptive_iters="<<adaptive.iterations
                 <<"\n";
    }
    std::cout<<"ADAPTIVE_SHAVING_NOTE: candidate reductions are allowed to"
                " differ from production; no equality or paving separation"
                " assumption is imposed on the adaptive variants."
                " Every excluded interval has a validated disjointness proof."
                " Witness checks are restricted to certified feasible points."
                " Only benchmark code changed.\n";
}

// The cached and uncached dyadic contractors traverse exactly the same
// cells, spending the same *visit* budget. A memoized image-exclusion
// proof only avoids another apply() to an identical cell in one box call.
// The exact final domain must therefore be identical for each query.
void compare_adaptive_cache(Selection const& selection,
                            std::string const& shape) {
    SolverProblem problem=make_problem(shape);
    std::vector<TraceRecord> const records=
        read_adaptive_trace(selection.trace_prefix,shape);
    std::vector<ExactBoxType> codomains;
    codomains.reserve(records.size());
    for(auto const& record:records) {
        ExactBoxType codomain=problem.codomains.front();
        for(SizeType i=0u;i<4u;++i) codomain[i]=record.box[i];
        codomains.push_back(std::move(codomain));
    }

    ConstraintSolver solver;
    AdaptiveShavingSettings const settings[]={
        {selection.adaptive_depth,selection.adaptive_budget},
        {selection.adaptive_depth,selection.adaptive_resolution_budget}
    };
    char const* labels[]={"budget","resolution"};
    for(unsigned method=0u;method!=2u;++method) {
        AdaptiveShavingStatistics uncached_stats;
        AdaptiveShavingStatistics cached_stats;
        SizeType changed_from_production=0u;
        SizeType extra_empty=0u;
        SizeType missed_empty=0u;
        SizeType separated=0u;
        for(SizeType q=0u;q<codomains.size();++q) {
            UpperBoxType production=problem.initial_domain;
            solver.reduce(production,problem.function,codomains[q]);
            UpperBoxType frozen=problem.initial_domain;
            reduce_full_scan_reference(
                solver,frozen,problem.function,codomains[q]);
            if(not same(production,frozen))
                throw std::runtime_error(
                    "adaptive cache baseline differs from frozen full scan"
                    " shape="+shape+" query="+std::to_string(q));
            if(definitely(production.is_empty())!=records[q].separated)
                throw std::runtime_error(
                    "adaptive cache baseline differs from paving separation"
                    " shape="+shape+" query="+std::to_string(q));

            UpperBoxType uncached=problem.initial_domain;
            reduce_adaptive_shaving(
                solver,uncached,problem.function,codomains[q],
                settings[method],&uncached_stats,false);
            UpperBoxType cached=problem.initial_domain;
            reduce_adaptive_shaving(
                solver,cached,problem.function,codomains[q],
                settings[method],&cached_stats,true);
            if(not same(cached,uncached))
                throw std::runtime_error(
                    "memoized adaptive shaving changed contracted domain"
                    " shape="+shape+" method="+labels[method]
                    +" query="+std::to_string(q));
            if(not same(uncached,production)) ++changed_from_production;
            Bool const before_empty=definitely(production.is_empty());
            Bool const candidate_empty=definitely(cached.is_empty());
            extra_empty+=static_cast<SizeType>(
                candidate_empty and not before_empty);
            missed_empty+=static_cast<SizeType>(
                not candidate_empty and before_empty);
            separated+=static_cast<SizeType>(candidate_empty);
        }
        if(uncached_stats.cache_hits!=0u
           or uncached_stats.evaluations!=uncached_stats.cell_visits
           or cached_stats.cell_visits!=uncached_stats.cell_visits
           or cached_stats.evaluations+cached_stats.cache_hits
              !=cached_stats.cell_visits
           or cached_stats.calls!=uncached_stats.calls
           or cached_stats.rounds!=uncached_stats.rounds
           or cached_stats.changed_calls!=uncached_stats.changed_calls
           or cached_stats.emptied_calls!=uncached_stats.emptied_calls
           or cached_stats.rejected_cells!=uncached_stats.rejected_cells
           or cached_stats.depth_stops!=uncached_stats.depth_stops
           or cached_stats.budget_stops!=uncached_stats.budget_stops)
            throw std::runtime_error(
                "memoized adaptive shaving changed visit sequence or counters");
        if(cached_stats.evaluations>cached_stats.calls*
              static_cast<SizeType>(settings[method].evaluation_budget))
            throw std::runtime_error("cached adaptive evaluation budget exceeded");

        std::cout<<"ADAPTIVE_CACHE_VERIFIED shape="<<shape
                 <<" method="<<labels[method]
                 <<" queries="<<records.size()
                 <<" contracted_boxes_changed_vs_production="
                    <<changed_from_production
                 <<" candidate_separated="<<separated
                 <<" candidate_only_empty="<<extra_empty
                 <<" production_only_empty="<<missed_empty
                 <<" box_calls="<<cached_stats.calls
                 <<" cell_visits="<<cached_stats.cell_visits
                 <<" uncached_evaluations="<<uncached_stats.evaluations
                 <<" cached_evaluations="<<cached_stats.evaluations
                 <<" cache_hits="<<cached_stats.cache_hits
                 <<" reused_pct="<<std::fixed<<std::setprecision(2)
                    <<(cached_stats.cell_visits
                        ?100.0*static_cast<double>(cached_stats.cache_hits)
                             /cached_stats.cell_visits:0.0)
                 <<" budget_stops="<<cached_stats.budget_stops
                 <<" depth_stops="<<cached_stats.depth_stops
                 <<"\n";
        if(selection.verify_only) continue;

        auto production_work=[&] {
            SizeType checksum=0u;
            for(auto const& codomain:codomains) {
                UpperBoxType domain=problem.initial_domain;
                solver.reduce(domain,problem.function,codomain);
                checksum+=static_cast<SizeType>(definitely(domain.is_empty()));
                keep_alive(domain);
            }
            return checksum;
        };
        auto adaptive_work=[&](bool cache_enabled) {
            SizeType checksum=0u;
            for(auto const& codomain:codomains) {
                UpperBoxType domain=problem.initial_domain;
                reduce_adaptive_shaving(
                    solver,domain,problem.function,codomain,
                    settings[method],nullptr,cache_enabled);
                checksum+=static_cast<SizeType>(definitely(domain.is_empty()));
                keep_alive(domain);
            }
            return checksum;
        };
        if(adaptive_work(false)!=separated
           or adaptive_work(true)!=separated)
            throw std::runtime_error("adaptive cache timed checksum mismatch");

        Timing const production_time=median_ms(production_work,selection);
        Timing const uncached_time=median_ms(
            [&]{return adaptive_work(false);},selection);
        Timing const cached_time=median_ms(
            [&]{return adaptive_work(true);},selection);
        double const cache_speedup=
            uncached_time.milliseconds/cached_time.milliseconds;
        double const production_speedup=
            production_time.milliseconds/cached_time.milliseconds;
        if(not std::isfinite(cache_speedup)
           or not std::isfinite(production_speedup)
           or cache_speedup<=0.0 or production_speedup<=0.0)
            throw std::runtime_error("invalid adaptive cache speedup");

        std::cout<<"ADAPTIVE_CACHE_COMPARE shape="<<shape
                 <<" method="<<labels[method]
                 <<" production_median_ms="<<std::fixed<<std::setprecision(5)
                     <<production_time.milliseconds
                 <<" uncached_median_ms="<<uncached_time.milliseconds
                 <<" cached_median_ms="<<cached_time.milliseconds
                 <<" cache_speedup="<<std::setprecision(3)<<cache_speedup
                 <<" production_to_cached_speedup="<<production_speedup
                 <<" samples="<<selection.samples
                 <<"\n";
    }
    std::cout<<"ADAPTIVE_CACHE_NOTE: each mode keeps the identical"
                " cell-visit budget and order for cached and uncached runs."
                " Memoization is local to one box_reduce; final boxes must"
                " match exactly. Different domains from P1 production are"
                " still diagnostic and are not endorsed for promotion."
                " Independent medians are fixed-order measurements.\n";
}

struct Datum {
    std::string shape;
    std::string operation;
    double milliseconds;
    SizeType queries;
    SizeType constraints;
    unsigned iterations;
};

void append_case(std::vector<Datum>& data,Selection const& selection,
                 std::string const& shape) {
    SolverProblem problem=make_problem(shape);
    ConstraintSolver solver;
    // The library Vector<ValidatedProcedure> uses UniformArray, which
    // requires element characteristics unavailable for Procedure. Cache
    // procedures in std::vector and compare identical scalar hull passes.
    std::vector<ValidatedProcedure> compiled;
    compiled.reserve(problem.function.result_size());
    for(SizeType i=0u;i<problem.function.result_size();++i) {
        compiled.emplace_back(problem.function[i]);
    }
    std::vector<UpperBoxType> after_hull;
    after_hull.reserve(problem.codomains.size());
    SizeType full_empty=0u;
    SizeType hull_empty=0u;
    SizeType box_calls=0u;

    for(ExactBoxType const& target:problem.codomains) {
        UpperBoxType direct=problem.initial_domain;
        solver.reduce(direct,problem.function,target);
        UpperBoxType repeat=problem.initial_domain;
        solver.reduce(repeat,problem.function,target);
        if(not same(direct,repeat))
            throw std::runtime_error("nondeterministic full reduction");
        full_empty+=static_cast<SizeType>(definitely(direct.is_empty()));

        UpperBoxType ordinary=problem.initial_domain;
        solver.hull_reduce(ordinary,problem.function,target);
        hull_empty+=static_cast<SizeType>(definitely(ordinary.is_empty()));
        after_hull.push_back(std::move(ordinary));

        // Both scalar paths use the same sequence and stop condition.
        UpperBoxType scalar_fresh=problem.initial_domain;
        UpperBoxType scalar_cached=problem.initial_domain;
        for(SizeType i=0u;i<target.size();++i) {
            solver.hull_reduce(scalar_fresh,problem.function[i],target[i]);
            solver.hull_reduce(scalar_cached,compiled[i],target[i]);
            if(not same(scalar_fresh,scalar_cached))
                throw std::runtime_error(
                    "fresh and cached scalar hull steps disagree");
            if(definitely(scalar_fresh.is_empty())) break;
        }
    }
    // Precomputed post-hull domains isolate a single subsequent box sweep.
    // The full reduce() may execute more than one hull/box round.
    auto box_sweep=[&] {
        SizeType empty=0u;
        for(SizeType q=0u;q<problem.codomains.size();++q) {
            UpperBoxType domain=after_hull[q];
            if(not definitely(domain.is_empty())) {
                ExactBoxType const& target=problem.codomains[q];
                Bool done=false;
                for(SizeType i=0u;i<target.size() and not done;++i) {
                    for(SizeType j=0u;j<domain.size();++j) {
                        ++box_calls;
                        solver.box_reduce(domain,problem.function[i],target[i],j);
                        if(definitely(domain.is_empty())) {
                            done=true;
                            break;
                        }
                    }
                }
            }
            empty+=static_cast<SizeType>(definitely(domain.is_empty()));
            keep_alive(domain);
        }
        return empty;
    };
    // Verify determinism and preconditioned behavior before timing.
    SizeType box_empty=box_sweep();
    if(box_sweep()!=box_empty)
        throw std::runtime_error("nondeterministic box sweep");
    SizeType initial_box_calls=box_calls;
    if(selection.verify_only) {
        std::cout<<"VERIFIED shape="<<shape
                 <<" queries="<<problem.codomains.size()
                 <<" constraints="<<problem.constraints_count
                 <<" full_empty="<<full_empty
                 <<" hull_empty="<<hull_empty
                 <<" box_sweep_empty="<<box_empty
                 <<" box_calls_per_sweep="<<initial_box_calls/2u
                 <<"\n";
        return;
    }

    auto full=[&] {
        SizeType empty=0u;
        for(ExactBoxType const& target:problem.codomains) {
            UpperBoxType domain=problem.initial_domain;
            solver.reduce(domain,problem.function,target);
            empty+=static_cast<SizeType>(definitely(domain.is_empty()));
            keep_alive(domain);
        }
        return empty;
    };
    auto hull_fresh=[&] {
        SizeType empty=0u;
        for(ExactBoxType const& target:problem.codomains) {
            UpperBoxType domain=problem.initial_domain;
            solver.hull_reduce(domain,problem.function,target);
            empty+=static_cast<SizeType>(definitely(domain.is_empty()));
            keep_alive(domain);
        }
        return empty;
    };
    auto hull_scalar_fresh=[&] {
        SizeType empty=0u;
        for(ExactBoxType const& target:problem.codomains) {
            UpperBoxType domain=problem.initial_domain;
            for(SizeType i=0u;i<target.size();++i) {
                solver.hull_reduce(domain,problem.function[i],target[i]);
                if(definitely(domain.is_empty())) break;
            }
            empty+=static_cast<SizeType>(definitely(domain.is_empty()));
            keep_alive(domain);
        }
        return empty;
    };
    auto hull_scalar_cached=[&] {
        SizeType empty=0u;
        for(ExactBoxType const& target:problem.codomains) {
            UpperBoxType domain=problem.initial_domain;
            for(SizeType i=0u;i<target.size();++i) {
                solver.hull_reduce(domain,compiled[i],target[i]);
                if(definitely(domain.is_empty())) break;
            }
            empty+=static_cast<SizeType>(definitely(domain.is_empty()));
            keep_alive(domain);
        }
        return empty;
    };
    SizeType scalar_empty=hull_scalar_fresh();
    if(full()!=full_empty or hull_fresh()!=hull_empty
       or hull_scalar_cached()!=scalar_empty or box_sweep()!=box_empty)
        throw std::runtime_error("reduction checksum mismatch");

    auto record=[&](char const* operation,auto&& workload) {
        Timing timing=median_ms(workload,selection);
        data.push_back({shape,operation,timing.milliseconds,
                        problem.codomains.size(),
                        problem.constraints_count,timing.iterations});
    };
    record("reduce_full",full);
    record("hull_vector_fresh",hull_fresh);
    record("hull_scalar_fresh",hull_scalar_fresh);
    record("hull_scalar_cached",hull_scalar_cached);
    record("box_sweep_post_hull",box_sweep);
}

void report(std::vector<Datum> const& data,Selection const& selection) {
    constexpr int width=121;
    std::cout<<"ConstraintSolver::reduce() phase benchmark\n"
             <<"Workload: 4D identity map plus five paving-equivalent"
                " constraint formulas; 32 fixed dyadic query boxes/shape.\n"
             <<"Prepared input functions and codomains are outside"
                " all timed loops. The box sample is NOT paving's adaptive"
                " collection.\n"
             <<"Measured workloads overlap: full reduce, vector hull,"
                " scalar hull (fresh/cached) and one post-hull box sweep.\n"
             <<"Fresh/cached scalar hulls have matching steps; they need"
                " not match the full-vector hull contractor.\n"
             <<"Independent within-case medians. DO NOT add/subtract"
                " phase medians to claim exact attribution.\n"
             <<"samples="<<selection.samples
             <<" min_sample_ms="<<selection.min_sample_ms
             <<" max_iterations="<<selection.max_iterations
             <<" manager=none  concurrency=0 (solver sequential)\n";
    std::cout<<std::string(width,'=')<<"\n"
             <<std::left<<std::setw(20)<<"Shape"<<"  "
             <<std::setw(22)<<"Operation"
             <<std::right<<std::setw(16)<<"Median (ms)"
             <<std::setw(10)<<"Queries"
             <<std::setw(14)<<"Constraints"
             <<std::setw(15)<<"us/query"
             <<std::setw(10)<<"Iter"<<"\n"
             <<std::string(width,'-')<<"\n";
    std::string previous;
    for(Datum const& row:data) {
        if(not previous.empty() and row.shape!=previous)
            std::cout<<std::string(width,'-')<<"\n";
        previous=row.shape;
        std::cout<<std::left<<std::setw(20)<<row.shape<<"  "
                 <<std::setw(22)<<row.operation
                 <<std::right<<std::fixed<<std::setprecision(5)
                 <<std::setw(16)<<row.milliseconds
                 <<std::setw(10)<<row.queries
                 <<std::setw(14)<<row.constraints
                 <<std::setprecision(3)
                 <<std::setw(15)<<(1000.0*row.milliseconds/row.queries)
                 <<std::setw(10)<<row.iterations<<"\n";
    }
    std::cout<<std::string(width,'=')<<"\n"
             <<"Speed ratios: phase median / full_reduce median"
                " (<1 means phase faster; not an additive cost breakdown).\n"
             <<std::string(67,'=')<<"\n"
             <<std::left<<std::setw(23)<<"Group"
             <<std::right<<std::setw(12)<<"Cases"
             <<std::setw(17)<<"Geomean"<<"\n"
             <<std::string(67,'-')<<"\n";
    for(char const* operation:{
            "hull_vector_fresh","hull_scalar_fresh",
            "hull_scalar_cached","box_sweep_post_hull"}) {
        SizeType cases=0u;
        double total_log=0.0;
        for(Datum const& row:data) {
            if(row.operation!=operation) continue;
            auto full=std::find_if(data.begin(),data.end(),
                [&](Datum const& candidate) {
                    return candidate.shape==row.shape
                        and candidate.operation=="reduce_full";
                });
            if(full==data.end())
                throw std::runtime_error("missing full reduce datum");
            double ratio=row.milliseconds/full->milliseconds;
            if(not std::isfinite(ratio) or ratio<=0.0)
                throw std::runtime_error("invalid phase ratio");
            total_log+=std::log(ratio);
            ++cases;
        }
        std::cout<<std::left<<std::setw(23)<<operation
                 <<std::right<<std::setw(12)<<cases
                 <<std::fixed<<std::setprecision(3)
                 <<std::setw(16)<<std::exp(total_log/cases)<<"x\n";
    }
    std::cout<<std::string(67,'=')<<"\n";
}

} // namespace

int main(int argc,char const* argv[]) {
    try {
        Logger::instance().configuration().set_verbosity(1u);
        Selection selection;
        if(not parse_arguments(argc,argv,selection)) return 0;
        unsigned const modes=
            static_cast<unsigned>(selection.compare_contained_skip)
            +static_cast<unsigned>(selection.compare_dependency_skip)
            +static_cast<unsigned>(selection.profile_residual_shaving)
            +static_cast<unsigned>(selection.compare_adaptive_shaving)
            +static_cast<unsigned>(selection.compare_adaptive_cache);
        if(modes>1u)
            throw std::invalid_argument(
                "comparison and profiling modes are mutually exclusive");
        if(not selection.trace_prefix.empty()) {
            for(std::string const& shape:selection.shapes) {
                if(selection.compare_contained_skip) {
                    compare_contained_skip(selection,shape);
                } else if(selection.compare_dependency_skip) {
                    compare_dependency_skip(selection,shape);
                } else if(selection.profile_residual_shaving) {
                    profile_residual_shaving(selection,shape);
                } else if(selection.compare_adaptive_shaving) {
                    compare_adaptive_shaving(selection,shape);
                } else if(selection.compare_adaptive_cache) {
                    compare_adaptive_cache(selection,shape);
                } else {
                    profile_adaptive_trace(
                        selection.trace_prefix,shape,selection.verify_only);
                }
            }
            return 0;
        }
        if(modes!=0u)
            throw std::invalid_argument(
                "comparison or profiling mode requires --trace-prefix");
        std::vector<Datum> results;
        for(std::string const& shape:selection.shapes) {
            std::cout<<"Measuring "<<shape
                     <<" | d=4 | queries=32\n";
            append_case(results,selection,shape);
        }
        if(not selection.verify_only) report(results,selection);
        return 0;
    } catch(std::exception const& error) {
        std::cerr<<"benchmark_constraint_reduce: "<<error.what()<<"\n";
        return 1;
    }
}
