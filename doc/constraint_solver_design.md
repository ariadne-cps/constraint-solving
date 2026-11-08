# Constraint solver design

**Status:** living design document, initialized 2026-11-08.  
**Owner:** `ariadne-cps/constraint-solving`.  
**Primary implementation:** [`src/constraint_solver.cpp`](../src/constraint_solver.cpp).  
**Public interfaces:** [`include/constraint-solving/constraint_solver.hpp`](../include/constraint-solving/constraint_solver.hpp).  
**Diagnostic benchmark:** [`benchmarks/benchmark_constraint_reduce.cpp`](../benchmarks/benchmark_constraint_reduce.cpp).

This is the **source of truth for solver design, optimization hypotheses,
experiments, and decisions**. Update it in the same commit as any future
nontrivial solver or diagnostic change. Performance numbers from other
repositories belong here when they inform the solver design; links to their
journals preserve workload provenance. This document describes existing
behavior separately from proposed behavior: an experimental benchmark
variant is **not** a production feature.

## 1. Scope and correctness contract

For a box `D`, a validated vector-valued function `f`, and an exact
codomain box `C`, the core problem is to find a tighter box `D'` while
preserving all possible solutions of

`x in D` and `f(x) in C`.

The active entry point is

```cpp
Bool ConstraintSolver::reduce(
    UpperBoxType& domain,
    ValidatedVectorMultivariateFunction const& function,
    ExactBoxType const& codomain) const;
```

`domain` is modified in place. A result of `true` means the
resulting domain has been certified empty. A result of `false`
**does not establish feasibility**, uniqueness, or even existence
of a solution. In particular, callers must not interpret a failure
to prove emptiness as an intersection proof.

The essential invariant for any proposed optimization is
`{x in original D | f(x) in C} subseteq final D`.
Validated enclosure and directed-rounding semantics must be preserved:
never drop a candidate using an approximate numerical comparison.
If a new method changes the contracted box but still preserves this
invariant, that is a change to contraction strength and may affect
downstream behavior. Our first experiments therefore require the
**stronger**, directly testable oracle of identical final upper boxes
for every replayed query, in addition to the same certified emptiness
result.

`ConstraintSolverInterface` also exposes `feasible()` and
`check_feasibility()`, which have distinct contracts from `reduce()`.
The solver additionally provides `hull_reduce`, `box_reduce`,
`monotone_reduce`, `lyapunov_reduce`, and `split`, as well as
the separate `propagate` family described below.

## 2. Current production reduction design

### 2.1 Full `reduce()`

The current production method:

1. Validates matching function/domain and result/codomain sizes,
   and returns immediately if the incoming domain is proved empty.
2. Computes a magnitude as the **sum of upper bounds** on coordinate
   interval widths.
3. Runs vector `hull_reduce(domain,function,codomain)`.
   If the domain is proved empty, returns `true`.
4. Runs `box_reduce(domain,function[i],codomain[i],j)` for every
   output `i` and input coordinate `j` in their original order,
   stopping immediately if the domain becomes empty.
5. Recomputes the magnitude and repeats the hull+box round while
   the magnitude is *strictly below* 0.75 times the previous
   magnitude (using the existing validated/rounded comparison).
   Otherwise returns `false`.

The 0.75 stopping factor and the original contractor order are part
of the current **operational** behavior. Do not silently change either
during a performance experiment. More than one full round can be needed
to prove emptiness or reach the present stopping threshold.

The vector `hull_reduce` path constructs a
`Vector<Procedure<ValidatedNumber>>` for the current function and
runs `simple_hull_reduce`. The scalar hull overload also has a
procedure-based counterpart. Procedure compilation/caching can change
cost, but scalar and vector hull paths must not be assumed
interchangeable without testing the actual contracted boxes.

### 2.2 One `box_reduce()` call

The contractor performs coordinate-wise interval *shaving* for
one scalar output:

- Return without work for an already degenerate target coordinate.
- Split that coordinate into eight slices.
- Evaluate the function on slices from the lower end until the first
  slice not proved disjoint from the constraint bounds. Certified
  disjoint slices may be shaved off.
- If **all** lower slices are proved disjoint, certify emptiness.
- Otherwise scan remaining upper slices backwards, shaving any
  provably disjoint suffix.
- Replace the target interval by the resulting enclosure.

The scan is **not a generic binary search**: membership/disjointness
of a nonlinear function's interval image is not generally monotonic
across coordinate slices. The sequence and early stopping protect
soundness, but may perform many evaluations.

The existing five-argument overload accumulates actual function
evaluations:

```cpp
Bool box_reduce(UpperBoxType& domain,
                ValidatedScalarMultivariateFunction const& function,
                ExactIntervalType const& bounds,
                SizeType variable,
                SizeType& function_evaluations) const;
```

The four-argument production call delegates to this implementation,
discarding the counter. The benchmark uses the counting overload
without duplicating the contractor algorithm.

### 2.3 Other existing propagation mechanisms

The `propagate()` APIs are **separate** from `reduce()`; their
results, stopping policies, and optional contractors must not be
conflated.

- A legacy `List<ValidatedConstraint>` overload performs hull,
  direct rejection checks, and shaving to a fixed point.
- The newer `std::vector<ConstraintPropagationConstraint>` path
  supports optional precompiled `ValidatedProcedure` objects,
  derivative-assisted monotone contraction, validated
  strict-bound rejection, and `ConstraintPropagationStatistics`.
- Its **cached hull-only** configuration has a dependency-driven
  agenda. Changed coordinates wake dependent constraints.
- Its general shaving path tracks function-coordinate dependencies,
  learns an active set of productive pairs, and refreshes all
  genuine dependencies before declaring a stalled fixed point.
  Missing cached procedures conservatively imply dependence on
  every coordinate.
- `propagate_hull_once()` provides one standalone hull sweep.
- `split()` is available separately. It is **not** currently
  invoked inside the active `reduce()` loop.

These are candidate sources of reusable design ideas, not evidence
that a new schedule is automatically equivalent to full `reduce()`.

## 3. Baseline measurements and provenance

### 3.1 Workloads and reproducibility

The profiling input comes from
[`ariadne-cps/paving`](https://github.com/ariadne-cps/paving).
The paving benchmark exports the **actual adaptive candidate boxes**
at dimension 4, grid-relative fineness 3, with the outcome of each
separation query. The constraint-solving benchmark reconstructs
the same five effective constraint formulas independently, then
replays the exact exported boxes. This avoids a reverse dependency
from constraint-solving to paving or set.

A *stdout log is not a trace file*. `--outer-trace-prefix PREFIX`
writes five `PREFIX-SHAPE-d4-f3.trace` files. Consumption requires
`--trace-prefix PREFIX` in the constraint-solving executable.

The earlier standalone test with **32 fixed dyadic query boxes per
shape** is useful for local contractor comparisons, but it does
not reproduce the adaptive box distribution. The measurements below
are from actual adaptive traces supplied by the operator on
2026-11-08; they are not a portable performance promise.

### 3.2 Full adaptive replay: hull versus box

On the verified 4D/f3 traces, the diagnostic replay matched
the final contracted upper box from the original production
`reduce()` for **every query** and matched all paving separation
bits. The replay times are instrumentation totals from a single
local run, **not** repeated calibrated medians or direct
end-to-end paving times.

| Shape | Queries | Solver rounds | Box calls | Hull ms | Box ms | Total ms |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| aligned-box | 2927 | 5254 | 223392 | 44.8297 | 457.2457 | 514.6919 |
| rotated-box | 2591 | 4733 | 203174 | 91.4492 | 494.3138 | 596.9049 |
| aligned-ellipsoid | 1127 | 4679 | 92145 | 63.8011 | 424.5983 | 494.9640 |
| rotated-ellipsoid | 1055 | 5841 | 115414 | 105.3823 | 661.0145 | 774.5302 |
| ball | 1467 | 5795 | 114344 | 79.8555 | 522.7505 | 611.3022 |

Across five replay runs: **9,167 queries**, **26,302 reduction
rounds**, **748,469 box calls**. Contractor calls consume roughly
**86%** of the instrumented solver time. Full replay also
uses an uninstrumented oracle run per query; the paired wall-clock
column is consequently *not* a standalone solver measurement.
Do not add/subtract overlapping microbenchmark phase medians to
estimate CPU attribution.

### 3.3 Box-contractor activity

The five-argument `box_reduce` overload records actual interval
function evaluations; the benchmark compares the targeted interval
before and after each call.

| Shape | Calls | Function evaluations | Changed | Unchanged | No evaluations |
| --- | ---: | ---: | ---: | ---: | ---: |
| aligned-box | 223392 | 278592 | 0 | 223392 | 84096 |
| rotated-box | 203174 | 315092 | 3540 | 199634 | 52284 |
| aligned-ellipsoid | 92145 | 214460 | 11351 | 80794 | 16370 |
| rotated-ellipsoid | 115414 | 291702 | 19082 | 96332 | 938 |
| ball | 114344 | 263612 | 13528 | 100816 | 19834 |
| **Total** | **748469** | **1363458** | **47501** | **700968** | **173522** |

**Measured:** 93.65% of calls leave their target interval unchanged;
23.18% perform zero function evaluations. **Unchanged is not
synonymous with useless**: those calls may still be needed for
validated exclusion or to confirm that no interval tail can be
removed.

For these specific workloads the four identity-image outputs
(indices 0–3) produce **no** box contractions. They account for
403,248 calls (53.88% of all box calls) and approximately 779 ms
(30.43% of box-contractor time across these replays). The
nonlinear constraint is output index 4 for the ellipsoids and
ball; the rotated-box productive constraints are outputs 4–7.
These observations motivate scheduling changes but do not
justify hard-coding output indices or omitting checks globally.

An earlier fixed-32-query study found about **2.44x** geometric-mean
speedup for cached versus freshly constructed **scalar hull**
procedures. That is not a measured full-`reduce()` speedup.

### 3.4 Benchmark contracts

The current benchmark can:

- Time `reduce_full`, vector/scalar hull paths and one post-hull
  box sweep on fixed input boxes.
- Replay adaptive traces and count full rounds, hull and box calls,
  evaluations, interval changes, certified empty results, and
  per-output/per-coordinate activity.
- Verify the replayed final upper box against unmodified
  `ConstraintSolver::reduce()` and the certified separation bit
  from paving.
- Compare an experimental certified-containment scheduling variant
  against the original on the *same adaptive* inputs.

No official README or CI workflow is to be edited merely for
experimental profiling. Keep benchmark construction opt-in and
benchmark measurements on the operator's local build unless
separately requested.

## 4. Optimization options, decisions, and status

Priority is based on the *measured bottleneck* and proof burden,
not on expected headline speedups. None of the following is yet
an accepted production-`reduce()` optimization.

### P1 — Certified containment pruning (benchmark prototype)

**Status:** implemented in
`benchmark_constraint_reduce.cpp` at commit
[`a62046b`](https://github.com/ariadne-cps/constraint-solving/commit/a62046bb3f62f53cbebe851403ba844f80855766);
**operator-verified on macOS Release, 2026-11-08:** all five adaptive
traces passed exact final-domain and separation-bit checks; calibrated
seven-sample medians show a speedup in every shape. The prototype
remains diagnostic only, with no production solver changes.

For output `i`, first evaluate a validated enclosure
`apply(function[i], domain)`. If
`definitely(subset(image, codomain[i]))`, the entire current
domain satisfies that output constraint, hence no slice of the
current domain can be certified inconsistent *by that constraint*
alone. Skip its coordinate-wise box contractors and cache the
certificate for the **remainder of this query only**.
Contractors only shrink the domain, so certified containment
remains valid on subsequent subdomains for the same function
and codomain. Do not cache certificates across queries or
different function/codomain pairs.

The precheck has a cost. Compare it against box work saved,
particularly on nonlinear constraints, and avoid claiming a
win before observing `CONTAINED_SKIP_VERIFIED` and
`CONTAINED_SKIP_COMPARE` on all five traces.

The experimental algorithm is intentionally **duplicated locally**
inside `tests/test_constraint_solver.cpp`; the benchmark keeps its
own implementation. Do not introduce a shared header between the
benchmark and the tests. This trades some maintenance duplication
for independent build boundaries; when the experimental algorithm
changes, synchronize the two copies deliberately and compare
against the unmodified production oracle. The new trace-free
test exercises both variants on independent scalar, coupled 2D,
and mixed 3D cases, including certified containment, boundary
contact, reuse across rounds, degenerate and empty domains,
and nonlinear shaving-dependent infeasibility. The initial local macOS compilation failed because the unqualified `apply` selected `std::apply` (the test file imports `std`). The test-local copy now explicitly calls `Ariadne::apply` with unchanged numerical behavior. Compilation and test results after this correction are still pending.

The prototype retains the original hull call, loop stopping
threshold, and iteration structure. It compares full final
contracted boxes and separation bits against the production
oracle for all traced queries; then times independent calibrated
full-sweep medians (rather than deriving a speedup from call
counts). This is an experiment, not a production code change.

### P2 — Static variable dependencies

**Status:** proposed; related machinery already implemented
in `propagate()`, not yet transferred to `reduce()`.

A contractor for `f_i` with respect to coordinate `j`
cannot shave based on changing `j` if `f_i` provably does
not depend on that coordinate. A dependency set derived from
the validated expression/procedure can skip such
output-coordinate pairs conservatively. If dependency analysis
is uncertain, fall back to all coordinates. Inspect
`propagation_procedure_dependencies()` before reuse;
it reads cached procedure instructions. Never use the
absence of observed contractions in a trace as a proof of
independence.

First test the dependency filter alone against the original
contractor order; only then compose it with P1.

### P3 — Incremental worklist / wake-up scheduling

**Status:** proposed for `reduce()`; an agenda exists already
for the cached hull-only `propagate()` configuration and an
adaptive shaving working set exists in its general path.

Maintain a dependency graph from changed variables to
constraints, enqueue only potentially affected contractors,
and retain an initial full activation. For any heuristic
working-set or scheduled round, define and test **refresh
semantics** before claiming a fixed point or stopping by the
current 0.75 reduction rule.

**Risk:** `reduce()` stops on shrinkage between complete
rounds, not solely on a mathematical fixed point. Reordering,
omitting, or postponing contractor calls may change the final
box even if the resulting enclosure remains sound. Require
explicit baseline comparison and re-establish the intended
stopping and strength contract before production adoption.

### P4 — Adaptive contractor ordering and cheap rejection

**Status:** hypothesis only.

Try ordering productive constraints first using online
cost/benefit (`new contraction / time`) or cheaply
certifiable direct rejection. Separating trivial image
coordinate constraints from nonlinear constraints may improve
latency, but must be implemented based on certified properties,
never by shape-specific indices. A new order may produce
different final domains under the existing stopping rule;
benchmark and validate before adopting.

### P5 — Better shaving search and subdivision

**Status:** exploratory, highest proof and maintenance cost.

`box_reduce` currently scans up to eight lower slices and
the remaining upper slices. Adaptive slice counts, reuse of
evaluations, sound bounding methods, and interval-specific
contractors may reduce evaluations for hard nonlinear cases.

Naive binary search of a boundary is **not** sound without
an established monotonicity predicate: separate slices can
alternate between provably excluded and undecided. Any
hierarchical branch-and-bound variant must conservatively
retain uncertain slices, preserve all feasible points, and
have a measured complexity budget. Subdivision is a distinct
algorithm from the present box-contractor sweep and may
increase work or change contraction strength.

### Cross-cutting — Precompiled procedures

**Status:** existing functionality in `propagate()`;
microbenchmark evidence for scalar hull paths only.

Evaluate caching across a whole `reduce()` call or across
compatible query contexts separately, verifying ownership,
function identity, precision, and thread safety.
Avoid introducing the invalid `Vector<ValidatedProcedure>`
generator construction; `UniformArray` requires
characteristics not available on `ValidatedProcedure`.
The scalar benchmark stores precompiled procedures in
`std::vector<ValidatedProcedure>`. Do not assume this is
directly suitable for the vector hull API.

### P1 measurement: certified-containment replay, 2026-11-08

The operator compiled and executed the opt-in
`--compare-contained-skip` benchmark on macOS Release with all five
d4/f3 adaptive traces. The correctness-only run emitted
`CONTAINED_SKIP_VERIFIED` for every shape. The subsequent timed run
reverified exact final-domain equality against production
`ConstraintSolver::reduce` **for every one of 9,167 queries**, plus
agreement with all exported paving separation bits. Seven calibrated
samples were taken for each independent whole-trace sweep, with
`--min-sample-ms 12` and one iteration per measured sweep.

| Shape | Baseline median (ms) | Guarded median (ms) | Speedup | Guard tests | New certificates | Reused certificates | Box calls saved | Box calls executed |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| aligned-box | 489.98808 | 51.74458 | 9.469x | 27924 | 27924 | 27924 | 223392 | 0 |
| rotated-box | 586.77229 | 206.16996 | 2.846x | 28266 | 22732 | 22644 | 181504 | 21670 |
| aligned-ellipsoid | 490.85296 | 354.58704 | 1.384x | 9082 | 4424 | 14208 | 74528 | 17617 |
| rotated-ellipsoid | 769.07454 | 602.28217 | 1.277x | 9913 | 4108 | 19112 | 92880 | 22534 |
| ball | 599.47400 | 431.18262 | 1.390x | 11521 | 5812 | 17344 | 92624 | 21720 |
| **Total counters** | | | | **86706** | **65000** | **101232** | **664928** | **83541** |

The new proof guard skips **88.84%** (664928/748469) of
the original box calls in this workload, and still gives identical
contracted boxes. The **geometric mean of the five per-shape median
speed ratios is 2.313x**. The separate sweep timings sum to 2936.16
and 1645.97 ms respectively; their ratio (1.784x) is a different,
workload-weighted summary and is **not** the geometric mean.

The aligned-box improvement is extreme because the containment
guard proves all outputs for all relevant rounds, eliminating
all box calls (the original 600 empty proofs came from hull
contraction). Other geometries still perform many nonlinear
box contractions. All timings are end-to-end **solver sweep**
times on the trace workload, not timings of `paving::separated()`.
Independent baseline and guarded measurements were taken in a fixed
order, so their medians do not eliminate drift/thermal/order effects;
a future benchmarking refinement can alternate order and report
multiple repeat batches. The proof check itself is not free.

**Design decision:** P1 is a promising optimization candidate,
especially for the axis-aligned workloads. Do not claim a generic
9.47x speedup or silently promote this prototype. Before promotion,
add direct unit/regression coverage (strict interval containment,
non-contained and nonlinear constraints, degenerate dimensions,
repeated rounds, empty domains), compile on supported platforms,
and measure against independent workloads. Keep current stopping
policy and validated semantics.

**Compiler hygiene:** the diagnostic comparison function contained
an unused local type alias `using Clock=std::chrono::steady_clock;`
(previously at benchmark line 709). Remove that alias only;
the timing utility and traced replay still use their own
`Clock` aliases. This fixes the source of the reported Clang
unused-local-typedef warning but has **no runtime or algorithmic
effect**; verify the warning disappears in the next local build.

## 5. Acceptance criteria for every experiment

1. **Baseline fidelity:** The unmodified `reduce()` result
   is the reference, with unchanged inputs and identical
   adaptive query order; compare each final domain, not only
   its empty/nonempty flag. Also compare all exported
   separation bits. Distinguish an algorithmic regression
   from an incorrect test expectation.
2. **Soundness:** Preserve validated outward enclosures and
   early emptiness proofs. When the purpose of a change is
   to alter contraction strength, document that separately,
   retain coverage of previously tested behavior, and
   prove the solution-enclosure invariant.
3. **Fair timing:** Build Release with matching compiler,
   dependencies, flags, hardware and concurrency. Exclude
   trace parsing, setup, correctness verification, and
   diagnostic output from paired timed kernels. Calibrate
   odd sample counts and report per-shape medians and
   geometric means of paired ratios, as in paving.
   Alternate execution order if drift matters; never
   infer an additive CPU budget from overlapping phases.
4. **Diagnostics:** Record number of queries, contractor
   calls, function evaluations, useful contractions,
   the fraction skipped, and reasons for termination.
   Keep trace details and all unverified results out of
   production performance claims.
5. **Regression coverage:** Use unit tests for containment,
   independent coordinates, degenerate domains,
   disjoint inputs, coupled nonlinear constraints,
   different variable dimensions, and multiple rounds.
   Check macOS, Ubuntu, Windows and sanitizer workflows
   when relevant. Do not modify test expectations merely
   to make a failed optimization pass.
6. **Promotion gate:** A proposal enters production
   `reduce()` only after demonstrated meaningful
   end-to-end gains on representative adaptive queries,
   preserved correctness, and maintained or explicitly
   justified contraction strength. Update this document
   in the same commit.

## 6. Reproducing the adaptive-trace experiments

These are commands for two **independent sibling repositories**,
assuming both have Release builds and starting from
`ariadne-cps/constraint-solving`:

```sh
mkdir -p "$HOME/ariadne-benchmark-traces"
(
  cd ../paving
  cmake --build build-release --target benchmark_paving_backends
  ./build-release/benchmarks/benchmark_paving_backends \
    --dimensions 4 --finenesses 3 \
    --shapes aligned-box,rotated-box,aligned-ellipsoid,rotated-ellipsoid,ball \
    --outer-trace-prefix "$HOME/ariadne-benchmark-traces/outer" \
    --concurrency 0
)
ls -lh "$HOME/ariadne-benchmark-traces/"*.trace
cmake --build build-release --target benchmark_constraint_reduce
./build-release/benchmarks/benchmark_constraint_reduce \
  --trace-prefix "$HOME/ariadne-benchmark-traces/outer" \
  --concurrency 0
./build-release/benchmarks/benchmark_constraint_reduce \
  --trace-prefix "$HOME/ariadne-benchmark-traces/outer" \
  --compare-contained-skip --concurrency 0
```

The **five actual `.trace` files** must exist before the
consumer starts. Redirecting `paving`'s stdout into a file
does not create a valid replay trace. `--verify-only`
can be added for correctness-only replay.

## 7. Decision log (append; never silently rewrite history)

| Date | Decision / evidence | Status |
| --- | --- | --- |
| 2026-11-07 | Move solver-level diagnosis to `constraint-solving`, retain paving workload provenance; avoid `set`/paving reverse dependencies. | Adopted architecture |
| 2026-11-07 | Fixed 32-query benchmark separated full reduce, hull variants and box sweep; fresh/cached *scalar* hull speed ratio ~2.44x geomean. | Measured microbenchmark; not an end-to-end speedup |
| 2026-11-07 | Actual d4/f3 adaptive replay verified original final domains and separated bits on 9,167 queries. Box contractor is dominant (~86% of instrumented time). | Measured |
| 2026-11-08 | 748,469 box calls, 1,363,458 evaluations; 700,968 unchanged and 173,522 evaluation-free calls. Identity outputs account for ~53.9% of calls but no observed contractions. | Measured |
| 2026-11-08 | Try per-query validated containment certificates before changing `reduce()`. Diagnostic-only variant at `a62046b`; no operator timing result yet. | Experiment implemented; evaluation pending |
| 2026-11-08 | Plan dependency-filter, incremental agenda, adaptive ordering, and sound shaving search only after isolated correctness and cost studies. | Proposed; not authorized as production changes |
| 2026-11-08 | macOS Release exact-domain replay passed 9,167 queries for certified-containment prototype; 664,928/748,469 box calls skipped, geometric-mean median speedup 2.313x (per shape 1.277-9.469x). Removed benchmark-only unused Clock alias warning. | Measured on adaptive paving traces; candidate for tested production promotion, not yet promoted |

| 2026-11-08 | Added trace-independent regression cases for certified containment, with a deliberate local copy of the benchmark algorithm in the unit-test file; no shared benchmark/test header or production change. | Tests submitted; compilation and execution pending |

| 2026-11-08 | Fixed test-only `apply` lookup ambiguity with libc++ `std::apply` by qualifying `Ariadne::apply`; benchmark and production algorithm unchanged. | Correction committed; test rerun pending |

### Update protocol

For every follow-up investigation: add a dated row to the
decision log; revise **Current production reduction design**
only when production code actually changes; record new
measurements with environment/trace provenance; and update
the **status and risks** of the affected proposal.
Commit these documentation changes alongside the relevant
benchmark/solver changes, rather than leaving design history
in a chat or solely in `paving/doc`.
