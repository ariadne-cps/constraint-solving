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
        if(selection.compare_contained_skip and selection.compare_dependency_skip)
            throw std::invalid_argument(
                "--compare-contained-skip and --compare-dependency-skip are exclusive");
        if(not selection.trace_prefix.empty()) {
            for(std::string const& shape:selection.shapes) {
                if(selection.compare_contained_skip) {
                    compare_contained_skip(selection,shape);
                } else if(selection.compare_dependency_skip) {
                    compare_dependency_skip(selection,shape);
                } else {
                    profile_adaptive_trace(
                        selection.trace_prefix,shape,selection.verify_only);
                }
            }
            return 0;
        }
        if(selection.compare_contained_skip or selection.compare_dependency_skip)
            throw std::invalid_argument(
                "comparison mode requires --trace-prefix");
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
