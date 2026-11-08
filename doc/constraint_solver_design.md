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

### 2.1 Production `reduce()`: certified containment scheduling

The production method, revised on 2026-11-08:

1. Validates matching function/domain and result/codomain sizes,
   returning immediately if the incoming domain is proved empty.
2. Initializes a Boolean cache of outputs that have been certified
   satisfied throughout a domain. The cache is **local to one call**.
3. Computes a magnitude as the sum of upper bounds on interval widths.
4. Runs vector `hull_reduce(domain,function,codomain)`; returns on
   certified emptiness.
5. Visits outputs in their original order. Previously certified
   outputs need no shaving. Otherwise a validated enclosure
   `Ariadne::apply(function[i],domain)` is checked against
   `codomain[i]`. Only if `definitely(subset(...))` succeeds is
   that output marked satisfied and its coordinate-wise
   `box_reduce` calls skipped. An uncertain or failed check
   executes the **complete original coordinate sweep**.
6. Recomputes the magnitude and repeats the hull+box round while
   strictly below the original 0.75 reduction threshold.

This leaves the hull contractor, the output and coordinate order,
the early-empty behavior, and the existing stop criterion unchanged.
Soundness follows because validated containment `f_i(D) subseteq C_i`
continues to hold after future contractions `D' subseteq D`.
No certificate is reused for a different `reduce()` query.

The **old full-scan algorithm is retained exclusively in**
`benchmarks/benchmark_constraint_reduce.cpp` as
`reduce_full_scan_reference()`, to preserve an independent baseline.
Unit tests exercise **only the production `ConstraintSolver::reduce()`
API**, with expected outcomes, feasible witnesses and reductivity
properties; they contain no copy of either full reduction algorithm
and do not include benchmark headers.

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
- Compare production certified containment against an independent
  frozen full-scan reference on the *same adaptive* inputs.
- Compare the additional static dependency mask as a diagnostic
  against production and the frozen reference, including both
  prepared-mask and per-query rebuilt-mask timings.
- Profile the *residual* calls after production P1 containment
  (not the legacy full-scan counts), including a 0–8 evaluation
  histogram, contraction status, time by output and coordinate,
  and full equality checks against production and legacy.

No official README or CI workflow is to be edited merely for
experimental profiling. Keep benchmark construction opt-in and
benchmark measurements on the operator's local build unless
separately requested.

## 4. Optimization options, decisions, and status

Priority is based on the *measured bottleneck* and proof burden,
not on expected headline speedups. P1 has been promoted to
production; the remaining items are experimental or proposed.

### P1 — Certified containment pruning (promoted to production)

**Status:** implemented in
`benchmark_constraint_reduce.cpp` at commit
[`a62046b`](https://github.com/ariadne-cps/constraint-solving/commit/a62046bb3f62f53cbebe851403ba844f80855766);
**operator-verified on macOS Release, 2026-11-08:** all five adaptive
traces passed exact final-domain and separation-bit checks; calibrated
seven-sample medians show a speedup in every shape. The prototype
was promoted to production only after the trace-independent tests
were confirmed green locally. Fresh production CI and adaptive
trace replay results are pending on the promotion commit.

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

The trace-independent unit tests now call only production
`ConstraintSolver::reduce()`. They verify contained and boundary
identity cases, contraction and feasible witness preservation,
empty and singleton domains, nonlinear shaving-dependent rejection,
coupled constraints over a deterministic 2D grid, and a 3D mixed
polynomial system. No legacy or experimental full-loop copy remains
in the test file, and there is no shared benchmark/test header.
The operator reported the previous regression suite green before
promotion; the newly revised production tests require a fresh run.

Production preserves the original hull call, loop stopping
threshold, and contractor order. The benchmark keeps a separate
frozen full-scan reference and the independently measured
experimental schedule. It compares the old full scan, experimental
contained-skip implementation, and production on every adaptive
query, then times old full scan versus production on calibrated
full-sweep medians. This ensures the original speedup experiment
remains meaningful after the promotion.

### P2 — Static variable dependencies

**Status: rejected as a general optimization for the uncompiled
production `reduce()` API (2026-11-08).** The complete adaptive
trace replay passed, but reconstructing static dependencies for
each call slowed all five shapes. The opt-in diagnostic remains
available in the benchmark; a future prepared-function interface
could enable a different tradeoff.

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

P1 is already in production, so the new isolated experiment
compares **the current production contained-skip reducer** against
a candidate that additionally avoids contractor pairs where the
validated procedure has no `IndexProcedureInstruction` referring to
that input coordinate. A variable is never ruled out on the basis
of observed ineffectiveness or approximate gradients.

The diagnostic constructs a conservative dependency mask from each
`ValidatedProcedure(function[i])`, then checks exact final-box
identity on every query against both current production and the
frozen old full-scan reference. It retains the original contractor
order and 0.75 round threshold. If the direct interval image is
certifiably disjoint from the codomain, it conservatively retains
the complete contractor sweep for that output, rather than losing
any immediate empty-domain proof in a coordinate with no dependency.

**Two timing scenarios are mandatory:** `prepared` builds masks
once for the collection of repeated queries with the same functions
(outside the timed loop); `rebuilt` includes new procedure
compilation for **every query**, modeling adoption through the
current uncompiled `reduce()` interface. The latter can lose all
nominal contractor savings to setup cost; do not extrapolate the
prepared result to production without a supported prepared-function
API and a lifecycle/thread-safety analysis.

Run from `constraint-solving` with existing d4/f3 traces:

```sh
cmake --build build-release --target benchmark_constraint_reduce
./build-release/benchmarks/benchmark_constraint_reduce \
  --trace-prefix "$HOME/Downloads/traces/outer" \
  --compare-dependency-skip --verify-only --concurrency 0
./build-release/benchmarks/benchmark_constraint_reduce \
  --trace-prefix "$HOME/Downloads/traces/outer" \
  --compare-dependency-skip --concurrency 0
```

Do not interpret a zero or small `dependency_skipped` count as a
bug: the remaining nonlinear outputs may use all input coordinates.
P2 is **not** to be promoted through the current API.

#### P2 actual adaptive replay results — macOS Release, 2026-11-08

The operator completed the correctness-only and timed
`--compare-dependency-skip` runs against five d4/f3 traces
at `$HOME/Downloads/traces/outer`. All five runs reported
`DEPENDENCY_SKIP_VERIFIED`: the candidate produced exactly
the same final upper boxes as the current P1 production solver
and frozen legacy full scan, and matched the exported paving
separation bits, on **9,167 adaptive queries**.

| Shape | Queries | Static dependency edges | Additional box calls skipped | Box calls executed | P1 median ms | Prepared median ms | Rebuilt median ms | Prepared speedup | Rebuilt speedup |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| aligned-box | 2927 | 12/48 | 0 | 0 | 52.57417 | 52.02363 | 84.74917 | 1.011x | 0.620x |
| rotated-box | 2591 | 16/48 | 10756 | 10914 | 208.11783 | 179.28250 | 219.88029 | 1.161x | 0.947x |
| aligned-ellipsoid | 1127 | 8/20 | 0 | 17617 | 354.84821 | 354.85913 | 365.53246 | 1.000x | 0.971x |
| rotated-ellipsoid | 1055 | 8/20 | 0 | 22534 | 583.24437 | 583.42346 | 598.61592 | 1.000x | 0.974x |
| ball | 1467 | 8/20 | 0 | 21720 | 427.80558 | 428.24942 | 442.52763 | 0.999x | 0.967x |

The filter eliminated **10,756** more box calls, exclusively
for rotated-box, reducing that shape from 21,670 to 10,914
surviving calls. That workload also triggered 156 conservative
disjointness fallback events. There were no extra eliminations
for the other four shapes, because P1 already skipped the
certified-satisfied constraints and the surviving nonlinear
ones use all four input coordinates. Static sparse edge counts
are not interchangeable with actual runtime savings.

The **geometric mean of the five per-shape median speedups**
is **1.032x prepared**, versus **0.883x rebuilt per query**.
Sums of medians are 1626.59 ms for P1, 1597.84 ms prepared,
and 1711.31 ms rebuilt; the summed ratio is a different metric.
Each complete sweep used seven samples, one iteration per sample.
These are independent fixed-order measurements and may include
order or thermal effects; they are not a universal runtime claim.

**Decision:** do **not** compile dependency masks inside each
production `reduce()` call. Retain P2 as an isolated benchmark
for future prepared-function interfaces or other sparse workloads.
Prioritize profiling the **surviving nonlinear shaving work**
before attempting adaptive search, contractor ordering, or new
subdivision strategies.

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

**Status:** residual shaving profiler added to the benchmark on
2026-11-08; compilation and operator results pending.
No production contractor change. The goal is to measure actual
post-P1 evaluations and contractor outcomes before selecting a
search algorithm.

The opt-in `--profile-residual-shaving` benchmark follows the
**exact production P1 hull/containment/box schedule**. It uses
the existing counting overload of `box_reduce`, instruments
individual calls, and checks all final boxes against both
production `reduce()` and the frozen pre-P1 full-scan reference,
as well as exported paving separation bits. It reports:

- `RESIDUAL_SHAVING_VERIFIED`: queries, rounds, P1 guard
  proofs/reuse, surviving box calls/evaluations and changes;
- `RESIDUAL_SHAVING_TIMES`: instrumented totals for hull,
  guard, and box calls, and evaluations associated with changed
  versus unchanged targeted intervals;
- `RESIDUAL_SHAVING_EVAL_HIST` and
  `RESIDUAL_SHAVING_UNCHANGED_HIST`: exact call counts with
  **0–8** interval function evaluations, both for all calls
  and for calls with no contracted coordinate;
- `RESIDUAL_SHAVING_GROUP`: counts and timings by output
  index or input coordinate; and
- `RESIDUAL_SHAVING_HOT_PAIR`: twelve most expensive
  output/coordinate pairs by measured box-contractor time.

The current eight-slice algorithm searches from both boundaries
and performs **at most eight** interval-function evaluations
per nondegenerate call. A count of zero reflects an immediate
return, usually for a degenerate target coordinate. Counting
an unchanged box does not prove the call was redundant.
The histogram distinguishes cheap early-stopping calls from
expensive inconclusive scans. These are **instrumented single
replay times**, not calibrated comparative speedups, and no
mathematical feasibility claim follows from a nonempty box.

From the `constraint-solving` directory, using the already
exported traces:

```sh
cmake --build build-release --target benchmark_constraint_reduce
./build-release/benchmarks/benchmark_constraint_reduce \
  --trace-prefix "$HOME/Downloads/traces/outer" \
  --profile-residual-shaving --verify-only --concurrency 0
./build-release/benchmarks/benchmark_constraint_reduce \
  --trace-prefix "$HOME/Downloads/traces/outer" \
  --profile-residual-shaving --concurrency 0
```

The first operator-provided `RESIDUAL_SHAVING_VERIFIED`
summary shows the following **post-P1** activity:

| Shape | Residual box calls | Interval evaluations | Evaluations per call | Contracting calls | Unchanged calls |
| --- | ---: | ---: | ---: | ---: | ---: |
| aligned-box | 0 | 0 | n/a | 0 | 0 |
| rotated-box | 21670 | 48000 | 2.215 | 3540 | 18130 |
| aligned-ellipsoid | 17617 | 91644 | 5.202 | 11351 | 6266 |
| rotated-ellipsoid | 22534 | 107526 | 4.772 | 19082 | 3452 |
| ball | 21720 | 110412 | 5.083 | 13528 | 8192 |

All 9,167 query results passed the replay equivalence checks.
These summary figures alone do **not** give the 0–8 evaluation
histograms. The fractions of costly unchanged calls and the
relative cost of each output-coordinate pair remain unknown
until the full profiler output is supplied.

#### P5 adaptive subdivision experiment (benchmark-only)

The eight-slice ceiling in production is heuristic: it is an
**interval-grid resolution**, not eight mandatory evaluations.
The existing contractor already stops each boundary scan at
the first undecidable slice. The new
`--compare-adaptive-shaving` benchmark explores a dyadic
boundary subdivision tree instead of a fixed eight-cell grid.
A cell is eliminated **only after a validated disjointness
proof**; inconclusive cells are subdivided until a leaf at
the configured resolution or an evaluation budget is reached.
An unproved leaf is never removed. The left and right boundary
passes each receive approximately half the evaluation budget,
ensuring a bounded amount of work on both sides.

Two independently calibrated experiments are compared with
production P1, on the same d4/f3 trace inputs:

- `method=budget`: total eight evaluations per
  `box_reduce` call and maximum dyadic depth eight;
- `method=resolution`: a provisional 64-evaluation ceiling
  with maximum dyadic depth eight (a nominal relative grid
  resolution of 1/256 of the original coordinate width).

These values are **configurable experimental parameters**,
not recommended defaults: `--adaptive-depth`,
`--adaptive-budget`, and
`--adaptive-resolution-budget`. Different adaptive
contractions can be stronger, weaker, or incomparable with
the current eight-slice contractor. Therefore correctness
cannot be defined by *equality* with production. The
benchmark checks that each candidate remains inside its
initial domain, preserves independently certified feasible
centre witnesses, and every discarded cell has a validated
proof of disjointness. It compares candidate/baseline domain
relations and separate certified separation bits, but does
not require them to match.

The benchmark reports `ADAPTIVE_SHAVING_VERIFIED` for each
shape/method (domain-strength relation, rejection differences,
feasible witness checks, evaluations, and budget/depth stops)
and `ADAPTIVE_SHAVING_COMPARE` for seven-sample calibrated
median timings. Production and the frozen pre-P1 full-scan
reference must still agree **exactly** on every trace query.

```sh
cmake --build build-release --target benchmark_constraint_reduce
./build-release/benchmarks/benchmark_constraint_reduce \
  --trace-prefix "$HOME/Downloads/traces/outer" \
  --compare-adaptive-shaving --verify-only --concurrency 0
./build-release/benchmarks/benchmark_constraint_reduce \
  --trace-prefix "$HOME/Downloads/traces/outer" \
  --compare-adaptive-shaving --samples 7 --concurrency 0
```

#### P5 correctness-only replay, 2026-11-08

The operator compiled the corrected benchmark and provided
ADAPTIVE_SHAVING_VERIFIED for both methods on all five adaptive
d4/f3 traces. Production P1 and the frozen legacy full scan
were checked for identical final boxes and paving separation bits
on all **9,167 queries**. The diagnostic prototype was checked
for domain reductivity and for preservation of the **certified
centre-point witnesses** in the recorded trace boxes. This is
a useful regression check, **not** a proof that all feasible
points are preserved or that each additional empty box is valid.
Separate correctness tests, including adversarial witnesses and
a careful audit of boundary-interval coverage, remain necessary.

| Shape | P1 evaluations | Budget-8 evaluations | Depth-8 / cap-64 evaluations | Budget-8 same / stronger / weaker / incomparable | Resolution same / stronger / weaker / incomparable |
| --- | ---: | ---: | ---: | --- | --- |
| aligned-box | 0 | 0 | 0 | 2927 / 0 / 0 / 0 | 2927 / 0 / 0 / 0 |
| rotated-box | 48000 | 141500 | 326468 | 2439 / 120 / 32 / 0 | 2223 / 304 / 64 / 0 |
| aligned-ellipsoid | 91644 | 147855 | 191479 | 1003 / 0 / 124 / 0 | 788 / 331 / 0 / 8 |
| rotated-ellipsoid | 107526 | 85208 | 448994 | 671 / 24 / 324 / 36 | 425 / 618 / 0 / 12 |
| ball | 110412 | 181704 | 243708 | 1387 / 0 / 80 / 0 | 1043 / 424 / 0 / 0 |
| **Total** | **357582** | **556267** | **1210649** | **8427 / 144 / 560 / 36** | **7406 / 1677 / 64 / 20** |

The budget-8 version performs **1.556x** as many interval
evaluations as P1, despite having the same *maximum per-call*
evaluation count: it induces different contractions, rounds
and calls, and frequently evaluates both halves of the
subdivision tree. It also produces weaker boxes in 560
queries versus stronger boxes in only 144. The resolution
variant performs **3.386x** as many evaluations, mostly for
stronger domain reductions (1677 queries), but 64 weaker
and 20 incomparable cases remain. The depth limit is
relative to the initial coordinate width, not an absolute
or scale-invariant stopping tolerance.

The resolution variant reports **28 additional empty boxes**
for rotated-ellipsoid, while P1 certifies none of those 28
empty. These are **candidate-only rejection claims**, not
independently established infeasible queries: the benchmark
checks a single feasible centre candidate per input only
when its image is independently certified inside the
codomain. The two methods agree on all recorded P1 empty
queries. Do not use these 28 new emptiness decisions in
production paving without stronger validation.

**Decision:** neither prototype is ready for production.
The next diagnostic is to distinguish repeated interval
evaluations of the same dyadic cell, redundant testing from
opposite ends, and the cost of stopping at uncertain leaves.
A reuse-aware/shallow-first variant should be compared
against P1 only after checking its exclusion proofs and
separation claims on independent cases. The **timed**
ADAPTIVE_SHAVING_COMPARE output was supplied later on
2026-11-08; see the seven-sample medians below. Evaluation
counts remain useful to explain the measured behavior, but
are not interchangeable with elapsed time.

#### P5 timed adaptive replay — macOS Release, 2026-11-08

The operator measured independently calibrated medians of
whole-trace sweeps with seven samples and one iteration per
sample. The following are **speedups relative to current P1
production `reduce()`**, not to the pre-P1 legacy contractor.

| Shape | Production median (ms) | Budget-8 median (ms) | Budget-8 speedup | Resolution median (ms) | Resolution speedup |
| --- | ---: | ---: | ---: | ---: | ---: |
| aligned-box | 51.92567 | 51.80229 | 1.002x | 51.61212 | 1.006x |
| rotated-box | 206.25858 | 282.87367 | 0.729x | 457.71250 | 0.451x |
| aligned-ellipsoid | 359.88067 | 493.93383 | 0.729x | 520.59354 | 0.691x |
| rotated-ellipsoid | 589.35817 | 363.47425 | **1.621x** | 1629.36462 | 0.362x |
| ball | 425.65896 | 592.16896 | 0.719x | 643.68737 | 0.661x |
| **Five-shape geometric mean** | | | **0.909x** | | **0.596x** |

The sum of the five medians is 1633.08 ms (production P1),
1784.25 ms (budget-8), and 3302.97 ms (resolution).
Their ratios differ from the geometric-mean speedups; neither
should be presented as a direct measurement of end-to-end
`paving::outer_approximation`.

The **only material performance win** is
`rotated-ellipsoid` under budget-8: **1.621x** faster.
However, on this shape its contracted box is weaker in
**324/1055** queries, stronger in only **24/1055**, and
incomparable in another 36. The reduced elapsed time is
therefore not an equivalent-strength speedup. This shape
also shows a drop from 22,534 to 10,994 box calls and from
107,526 to 85,208 interval evaluations, partly due to
the different sequence of contractions and solver rounds.
In contrast, the other nonlinear shapes become slower.

The depth-8, cap-64 resolution approach is slower on every
nontrivial shape; its extra contractions do not justify
calling it a performance optimization. The **28 additional
empty results** on rotated-ellipsoid should remain
flagged as candidate-only until independently validated.
The user-provided logs prove the benchmark's programmed
invariants passed, not universal soundness or parity in
contraction strength.

**Decision:** reject both adaptive variants for production
as-is. Do not tune the global slice count or switch to
dyadic search solely because of the rotated-ellipsoid
speedup. If investigating this family further, instrument
repeated cell evaluations, overlapping passes, and
evaluation-cost-to-contraction benefit on the affected
output-coordinate pairs, and strengthen adversarial
correctness tests before changing the contractor.

#### P5 evaluation reuse experiment (benchmark only, 2026-11-08)

The original adaptive prototype traverses dyadic cells separately
from the lower and upper boundaries, revisiting the root and
sometimes other dyadic nodes. Its new opt-in mode,
`--compare-adaptive-cache`, records the identity of a node
(`depth,index`) and its validated exclusion decision within
a **single coordinate contractor invocation**. On a repeated
visit, the contractor reuses that decision without repeating
`apply(function,slice)`. It must never reuse a proof after
the domain or another coordinate has changed.

The evaluation **budget is counted in visited cells, not
cached misses**. This deliberately preserves the old candidate's
traversal order and stopping points, yielding *exactly the
same contracted boxes* between cached and uncached adaptive
variants. The benchmark checks this equality for every adaptive
query and confirms matching call, round, contraction and
budget/depth-stop counters. It also checks the untouched
production solver against the frozen legacy reference and
the original paving separation bits.

For each of `budget` (8 visits) and `resolution` (64 visits)
the benchmark prints `ADAPTIVE_CACHE_VERIFIED` with visits,
actual function evaluations, cache hits, and separation counts;
and `ADAPTIVE_CACHE_COMPARE` with seven-sample median
times for production, uncached dyadic traversal and cached
dyadic traversal. The instrumentation uses a per-call vector
of memoized decisions, so lookup and allocation overhead are
included in the real timing; a smaller evaluation count does
not imply a wall-time improvement.

```sh
cmake --build build-release --target benchmark_constraint_reduce
./build-release/benchmarks/benchmark_constraint_reduce \
  --trace-prefix "$HOME/Downloads/traces/outer" \
  --compare-adaptive-cache --verify-only --concurrency 0
./build-release/benchmarks/benchmark_constraint_reduce \
  --trace-prefix "$HOME/Downloads/traces/outer" \
  --compare-adaptive-cache --samples 7 --concurrency 0
```

**Scope and decision:** the cache is an implementation optimization
of the *diagnostic dyadic algorithm*, not of production's
eight-slice `box_reduce()`. The operator has now verified and
timed the cached and uncached versions; their final boxes match
exactly on all 9,167 adaptive queries. The extra cost of dyadic
subdivision, weaker contractions and 28 candidate-only empty
results remain; **do not promote** this cache or the dyadic
algorithm to production. The current production eight-slice
scan already avoids rescanning the examined lower slices from
its upper pass, so these dyadic hit rates cannot simply be
transferred to production.

#### P5 cached dyadic replay results — macOS Release, 2026-11-08

The operator ran `--compare-adaptive-cache` with both
`--verify-only` and seven-sample timed sweeps, using the five
d4/f3 adaptive traces at `$HOME/Downloads/traces/outer`.
Every mode passed **exact final-box equality between cached
and uncached dyadic contractors** on all 9,167 queries, with
identical cell visits, calls, rounds, contractor outcomes,
budget stops, and separation decisions. The production P1
baseline was also checked against the frozen old full scan
and the paving separation trace.

| Shape | Budget-8 hits / visits | Budget-8 cache speedup | Budget-8 P1-to-cache speedup | Resolution hits / visits | Resolution cache speedup | Resolution P1-to-cache speedup |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| aligned-box | 0 / 0 | 1.027x | 1.002x | 0 / 0 | 1.026x | 1.081x |
| rotated-box | 19268 / 141500 | 1.074x | 0.797x | 25428 / 326468 | 1.045x | 0.466x |
| aligned-ellipsoid | 32864 / 147855 | 1.180x | 0.857x | 39040 / 191479 | 1.208x | 0.825x |
| rotated-ellipsoid | 15988 / 85208 | 1.168x | **1.877x** | 34190 / 448994 | 1.073x | 0.387x |
| ball | 39758 / 181704 | 1.179x | 0.842x | 47486 / 243708 | 1.196x | 0.790x |
| **Totals / geometric mean speedups** | **107878 / 556267** | **1.124x** | **1.016x** | **146144 / 1210649** | **1.107x** | **0.662x** |

Budget-8 memoization avoids **107,878 (19.39%)** repeated
validated function evaluations, reducing 556,267 to
448,389, and accelerates the original adaptive prototype
by 1.124x geometric mean. Resolution memoization avoids
**146,144 (12.07%)**, reducing 1,210,649 to 1,064,505
evaluations, for a 1.107x cached-vs-uncached speedup.
The medians include cache allocation and linear lookup costs.
The aligned-box measurements have **zero box work** and are
therefore just measurement variability, not cache benefits.

The budget-8 cached variant has a five-shape geomean
**1.016x** relative to P1, but this number is misleading
as a promotion criterion: four shapes are slower or
approximately equal, while the sole substantial win,
`rotated-ellipsoid` (**1.877x**), is still weaker on
324/1055 final boxes and differs in another 60.
The resolution variant remains substantially slower
than P1 (geomean 0.662x), and its 28 additional
`rotated-ellipsoid` empty results remain unverified
by an independent complete feasibility argument.

The production benchmark and diagnostic algorithms have
different contraction behavior. The equality achieved by
cached versus uncached dyadic traversal is a meaningful
**implementation-level isolation test**; it does not confer
equivalence between dyadic traversal and production.
Do not infer a statistically robust few-percent gain from
the fixed-order unpaired medians, especially where there
is no contractor work.

#### Integrated exact audit of rotated-ellipsoid extra empties

The exact geometry auditor is integrated into
`benchmarks/benchmark_constraint_reduce.cpp`, and runs automatically
when `--compare-adaptive-shaving` replays `rotated-ellipsoid` queries.
The separate Python script has been removed. Both the budget-8 and
resolution (default budget-64) adaptive candidates are audited. The
`ADAPTIVE_SHAVING_EXTRA_EMPTY` messages retain zero-based query IDs.

The auditor independently minimizes the *actual convex quadratic*
`Q=16(3u+4v)^2+64(-4u+3v)^2+400(w^2+z^2)`, with
`(u,v,w,z)=(x0-1/2,x1-1/2,x2-1/2,x3-1/2)` and `Q<=25`.
It intersects each trace box with `[0,1]^4`, converts the original
`strtod` binary-double endpoints into arbitrary-precision
`boost::multiprecision::cpp_rational` values *without rounding*,
and checks the stationary interior point, all four edges and their
corners. On each fixed-u edge the minimizing v is `9u/13`;
on each fixed-v edge the minimizing u is `36v/73`.
The w,z minima follow by clamping to `1/2`. All
comparisons with 25 are rational and exact. For a nonempty
intersection whose true minimum is at most 25, the oracle returns
a feasible rational witness.

The audit checks **every baseline paving separation bit** and **every
candidate-only empty result** (for both methods) against this
independent oracle, emits `AUDIT_ELLIPSOID_QUERY` and
`AUDIT_ELLIPSOID_SUMMARY` diagnostics, and exits nonzero
if any claim is contradicted. It takes place before timed sweeps;
it does not modify the production solver.

```sh
cmake --build build-release --target benchmark_constraint_reduce
./build-release/benchmarks/benchmark_constraint_reduce --audit-ellipsoid-self-test
./build-release/benchmarks/benchmark_constraint_reduce \
  --trace-prefix "$HOME/Downloads/traces/outer" \
  --shapes rotated-ellipsoid \
  --compare-adaptive-shaving --verify-only --concurrency 0
```

To also run the benchmark timing after the exact audit, remove
`--verify-only`. The seven built-in self-test cases cover
the interior minimum, an exterior point, the exact threshold,
an edge minimum, initial-domain disjointness, binary-rational
reconstruction and a just-outside-threshold value.

**Scope:** this proves or refutes separation claims only for the
specific rotated-ellipsoid geometry and the trace boxes actually
replayed. It does **not** prove adaptive dyadic contractor soundness
for arbitrary functions or coordinates. Previously reported
28 additional empties remain unverified until this C++ replay is
run on the operator's original traces.

**Decision:** retain the memoization diagnostic for reuse
studies but do not insert it into production. Before any
further shaving redesign, independently audit the dyadic
interval partition/exclusion proof (especially the 28
additional empty results), and seek a strategy that
minimizes evaluation cost **without sacrificing contraction
strength** on nonlinear inputs.

Do **not** promote an adaptive variant merely because its
loop is cheaper. Its contraction strength, certified
separation count, and feasible-witness preservation are
coequal acceptance criteria, and its evaluation budget is
part of the algorithm's practical semantics. Changing the
production contractor requires additional unit tests with
adversarial interval cases and independent workloads.

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
    --outer-trace-prefix "$HOME/Downloads/traces/outer" \
    --concurrency 0
)
ls -lh "$HOME/ariadne-benchmark-traces/"*.trace
cmake --build build-release --target benchmark_constraint_reduce
./build-release/benchmarks/benchmark_constraint_reduce \
  --trace-prefix "$HOME/Downloads/traces/outer" \
  --concurrency 0
./build-release/benchmarks/benchmark_constraint_reduce \
  --trace-prefix "$HOME/Downloads/traces/outer" \
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

| 2026-11-08 | After green local test run, promote validated containment caching into vector `ConstraintSolver::reduce()`. Keep the frozen old full-scan only in benchmark as an independent speed/equivalence baseline; unit tests exercise production behavior and invariants directly, with no experimental loop copy. | Production change committed; cross-platform CI and post-promotion trace comparison pending |

| 2026-11-08 | Added opt-in P2 static-dependency benchmark with certified procedure-index masks; verifies against production and frozen legacy full scan on adaptive traces, and measures both prepared-once and rebuild-per-query costs. | Diagnostic implementation committed; correctness and speed pending |

| 2026-11-08 | P2 replay validated identical final domains and certified separations on 9,167 queries; only rotated-box skipped 10,756 more box calls. Geometric-mean speedups were 1.032x prepared and 0.883x rebuilt per query. | Reject general per-call dependency compilation; retain P2 diagnostic and profile residual nonlinear shaving next |

| 2026-11-08 | Added opt-in P1-residual shaving profiler: 0–8 box evaluation histograms, unchanged-cost histogram, output/variable hot pairs and guard/hull/box timing. Verify every final box against production and frozen full scan. | Diagnostic implementation committed; correctness and measured results pending |

| 2026-11-08 | Add benchmark-only adaptive dyadic boundary shaving prototypes: 8-evaluation budget versus depth-limited resolution with a higher cap, validated disjointness proofs, and measurement of contraction-strength changes relative to current P1. | Experiment implemented; compile, independent correctness checks, and timings pending |

| 2026-11-08 | Fixed adaptive benchmark build failure: the solver's file-local `affine()` was not visible in the benchmark. Added a benchmark-local `adaptive_affine()` with identical interpolation arithmetic and removed an unused evaluation variable. | Source corrected; local rebuild/verification pending |

| 2026-11-08 | P5 adaptive correctness-only replay passed baseline equivalence on 9,167 queries; budget-8 used 556,267 versus 357,582 production evaluations, resolution-cap-64 used 1,210,649. The resolution candidate returned 28 additional empty results on rotated-ellipsoid, not independently proved. | Keep both prototypes diagnostic; inspect duplicated dyadic evaluations and strengthen correctness validation before any promotion |

| 2026-11-08 | Completed P5 seven-sample adaptive replay timings. Geomean P1-relative speedup: budget-8 0.909x, depth-8/cap-64 0.596x. Rotated-ellipsoid budget-8 alone achieved 1.621x but weakened 324/1055 boxes. | Reject both adaptive variants for production; retain diagnostic and investigate duplicated evaluations and contraction strength separately |

| 2026-11-08 | Added per-box-call dyadic exclusion-decision cache with identical visit order and visit budget; opt-in benchmark reports verified domain equivalence to uncached adaptive shaving, cache hits, evaluation reduction and timed comparison including cache overhead. | Prototype committed; build and operator results pending |

| 2026-11-08 | Adaptive dyadic per-box cache validated identical domains versus uncached on all 9,167 queries. Cache saves 107,878 evaluations (19.39%) for budget-8 and 146,144 (12.07%) for resolution; cache-only geomeans 1.124x and 1.107x. Relative to production: 1.016x and 0.662x, with the only major budget win on rotated-ellipsoid still weakening 324/1055 boxes. | Retain benchmark diagnostics, no production change; prioritize independent empty-proof audit and strength-preserving shaving |

| 2026-11-08 | Added a shape-specific independent exact-arithmetic verifier for the 28 candidate-only rotated-ellipsoid empty results. The benchmark now emits candidate-only zero-based query indices; the Python auditor computes the exact convex quadratic minimum on every query box and checks baseline separation bits. | Auditor self-test passed in development; operator replay on actual traces and C++ rebuild pending |
| 2026-11-08 | Moved the independent rotated-ellipsoid oracle into `benchmark_constraint_reduce.cpp`, using arbitrary-precision rational arithmetic on original binary-double trace endpoints; both adaptive modes and all baseline separation bits are audited before timing. Removed the Python script and added a C++ seven-case self-test mode. | Standalone oracle compiled and passed seven self-test cases; complete Ariadne build and replay on operator traces remain pending |

### Update protocol

For every follow-up investigation: add a dated row to the
decision log; revise **Current production reduction design**
only when production code actually changes; record new
measurements with environment/trace provenance; and update
the **status and risks** of the affected proposal.
Commit these documentation changes alongside the relevant
benchmark/solver changes, rather than leaving design history
in a chat or solely in `paving/doc`.
