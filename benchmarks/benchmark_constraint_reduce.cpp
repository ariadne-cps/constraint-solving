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
};

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
                auto bstart=Clock::now();
                ++stats.box_calls;
                solver.box_reduce(domain,function[i],codomain[i],j);
                stats.box_seconds+=elapsed(bstart);
                if(definitely(domain.is_empty())) {
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
    SizeType mismatch=0u;
    auto sweep_start=std::chrono::steady_clock::now();
    for(TraceRecord const& entry:records) {
        // Reuse exactly the geometry and constraints from this shape; only
        // the query-box portion varies with the exported adaptive trace.
        ExactBoxType codomain=problem.codomains.front();
        for(SizeType i=0u;i<4u;++i)
            codomain[i]=entry.box[i];
        UpperBoxType original=problem.initial_domain;
        solver.reduce(original,problem.function,codomain);
        UpperBoxType instrumented=problem.initial_domain;
        traced_reduce(solver,instrumented,problem.function,codomain,stats);
        if(not same(original,instrumented))
            throw std::runtime_error(
                "instrumented reduce changed contracted domain for "+shape);
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
    std::cout<<"SOLVER_TRACE measures an instrumented diagnostic replay and"
                " a separate uninstrumented oracle call per input;"
                " it is NOT a direct timing of paving::separated().\n";
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
        if(not selection.trace_prefix.empty()) {
            for(std::string const& shape:selection.shapes) {
                profile_adaptive_trace(
                    selection.trace_prefix,shape,selection.verify_only);
            }
            return 0;
        }
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
