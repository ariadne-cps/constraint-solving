#!/usr/bin/env python3
"""Independently audit candidate-only empty boxes for rotated-ellipsoid.

The benchmark's geometry is Q <= 25, where
  Q=16*(3*u+4*v)^2+64*(-4*u+3*v)^2+400*(w^2+z^2)
and (u,v,w,z)=(x0-1/2,x1-1/2,x2-1/2,x3-1/2).
Q is a strictly convex quadratic. Its exact minimum over a closed
4D box is obtained by minimizing the separable w,z terms and
checking the stationary interior point and the four boundary edges
of the u,v rectangle. Arithmetic uses stdlib Fraction throughout.
The result is a geometry-specific oracle, not a general theorem
about the experimental adaptive contractor. GPL-3.0-or-later.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
import math
from pathlib import Path
import re
import subprocess
import sys

H = Fraction(1, 2)
MAGIC = "ARIADNE_PAVING_SOLVER_TRACE_V1"
EXTRA = re.compile(
    r"^ADAPTIVE_SHAVING_EXTRA_EMPTY shape=rotated-ellipsoid "
    r"method=resolution query=(\d+)$", re.MULTILINE)
SUMMARY = re.compile(
    r"^ADAPTIVE_SHAVING_VERIFIED shape=rotated-ellipsoid "
    r"method=resolution\b[^\n]*\bcandidate_only_empty=(\d+)\b",
    re.MULTILINE)


def binary_fraction(token: str) -> Fraction:
    """Reproduce benchmark strtod -> binary double -> ExactDouble."""
    lower = token.lower()
    value = float.fromhex(token) if ("0x" in lower or "p" in lower) else float(token)
    if not math.isfinite(value):
        raise ValueError(f"nonfinite trace endpoint: {token}")
    return Fraction.from_float(value)


def clamp(x: Fraction, lo: Fraction, hi: Fraction) -> Fraction:
    return min(hi, max(lo, x))


def Q(u: Fraction, v: Fraction, w: Fraction, z: Fraction) -> Fraction:
    return 16*(3*u+4*v)**2+64*(-4*u+3*v)**2+400*(w*w+z*z)


def exact_minimum(box: tuple[tuple[Fraction, Fraction], ...]):
    """Return (minimum Q, feasible minimizer) or None for box disjoint from [0,1]^4.

    The u,v block has Hessian [[2336,-1152],[-1152,1664]], positive
    definite. On fixed-u edges the optimum is v=9*u/13; on fixed-v
    edges it is u=36*v/73. Include (0,0) when inside the rectangle.
    All comparisons, evaluations and witnesses are exact rationals.
    """
    if len(box) != 4:
        raise ValueError("expected 4-dimensional query box")
    domain = [(max(Fraction(0), lo), min(Fraction(1), hi)) for lo, hi in box]
    if any(lo > hi for lo, hi in domain):
        return None
    (ul, uu), (vl, vu) = [(lo-H, hi-H) for lo, hi in domain[:2]]
    candidates = []
    if ul <= 0 <= uu and vl <= 0 <= vu:
        candidates.append((Fraction(0), Fraction(0)))
    for u in (ul, uu):
        candidates.append((u, clamp(Fraction(9, 13)*u, vl, vu)))
    for v in (vl, vu):
        candidates.append((clamp(Fraction(36, 73)*v, ul, uu), v))
    for u in (ul, uu):
        for v in (vl, vu):
            candidates.append((u, v))
    w = clamp(H, *domain[2])-H
    z = clamp(H, *domain[3])-H
    u, v = min(candidates, key=lambda uv: Q(uv[0], uv[1], w, z))
    witness = (u+H, v+H, w+H, z+H)
    assert all(lo <= x <= hi for x, (lo, hi) in zip(witness, domain))
    return Q(u, v, w, z), witness


def read_trace(prefix: Path):
    path = Path(f"{prefix}-rotated-ellipsoid-d4-f3.trace")
    data = path.read_text(encoding="utf-8").split()
    if len(data) < 5:
        raise ValueError("truncated trace header")
    if tuple(data[:4]) != (MAGIC, "rotated-ellipsoid", "4", "3"):
        raise ValueError("unexpected trace geometry or header")
    n = int(data[4])
    if n <= 0 or len(data) != 5 + 9*n:
        raise ValueError("truncated or overlong trace payload")
    records = []
    for i in range(n):
        row = data[5+9*i:5+9*(i+1)]
        if row[0] not in ("0", "1"):
            raise ValueError(f"invalid separation bit on query {i}")
        box = tuple((binary_fraction(row[1+2*j]),
                     binary_fraction(row[2+2*j])) for j in range(4))
        if any(lo > hi for lo, hi in box):
            raise ValueError(f"reversed interval on query {i}")
        records.append((row[0] == "1", box))
    return records


def audit(content: str, records):
    indices = [int(q) for q in EXTRA.findall(content)]
    summary = SUMMARY.findall(content)
    if len(summary) != 1:
        raise ValueError("expected one resolution summary; rebuild the benchmark")
    expected = int(summary[0])
    if len(indices) != expected or len(set(indices)) != len(indices):
        raise ValueError(f"reported {expected} extra empties, but received {len(indices)} distinct indices")
    if any(q < 0 or q >= len(records) for q in indices):
        raise ValueError("out-of-range candidate-only query ID")

    baseline_invalid = []
    true_empty = 0
    for i, (baseline_empty, box) in enumerate(records):
        minimum = exact_minimum(box)
        actual_empty = minimum is None or minimum[0] > 25
        true_empty += actual_empty
        if baseline_empty and not actual_empty:
            baseline_invalid.append(i)

    candidate_invalid = []
    for q in sorted(indices):
        baseline_empty, box = records[q]
        if baseline_empty:
            raise ValueError(f"candidate-only query {q} is marked empty in baseline")
        minimum = exact_minimum(box)
        if minimum is None:
            print(f"AUDIT_ELLIPSOID_QUERY query={q} status=EXACTLY_EMPTY box_disjoint_from_initial_domain=true")
        elif minimum[0] > 25:
            value = minimum[0]
            print(f"AUDIT_ELLIPSOID_QUERY query={q} status=EXACTLY_EMPTY min_num={value.numerator} min_den={value.denominator} threshold=25")
        else:
            candidate_invalid.append(q)
            print(f"AUDIT_ELLIPSOID_QUERY query={q} status=COUNTEREXAMPLE min_Q={minimum[0]} feasible_point={minimum[1]}")
    print("AUDIT_ELLIPSOID_SUMMARY"
          f" queries={len(records)}"
          f" true_empty={true_empty}"
          f" baseline_separated={sum(int(b) for b, _ in records)}"
          f" baseline_false_positive={len(baseline_invalid)}"
          f" candidate_only_reported={len(indices)}"
          f" candidate_only_proved_empty={len(indices)-len(candidate_invalid)}"
          f" candidate_only_false_positive={len(candidate_invalid)}")
    if baseline_invalid:
        print(f"AUDIT_ELLIPSOID_BASELINE_INVALID queries={baseline_invalid}")
    return not candidate_invalid and not baseline_invalid


def self_test():
    point = lambda x: (Fraction(x), Fraction(x))
    full = ((Fraction(0), Fraction(1)),)*4
    assert exact_minimum(full) == (0, (H,)*4)
    assert exact_minimum((point(1), point(1), point(H), point(H)))[0] == 212
    assert exact_minimum((point(H), point(H), point(Fraction(3,4)), point(H)))[0] == 25
    edge = exact_minimum((point(1), (Fraction(0), Fraction(1)), point(H), point(H)))
    assert edge[1][1] == Fraction(11,13)
    assert edge[0] == Q(H, Fraction(9,26), 0, 0)
    assert exact_minimum((point(2), point(H), point(H), point(H))) is None
    print("AUDIT_ELLIPSOID_SELF_TEST_OK cases=5")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trace-prefix", type=Path)
    parser.add_argument("--benchmark", type=Path,
                        help="run corrected benchmark to retrieve candidate-only query IDs")
    parser.add_argument("--benchmark-output", type=Path,
                        help="read output from a corrected --verify-only benchmark run")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    if args.trace_prefix is None:
        parser.error("--trace-prefix is required")
    if bool(args.benchmark) == bool(args.benchmark_output):
        parser.error("specify exactly one of --benchmark and --benchmark-output")
    records = read_trace(args.trace_prefix)
    if args.benchmark:
        command = [str(args.benchmark.resolve()),
                   "--trace-prefix", str(args.trace_prefix),
                   "--shapes", "rotated-ellipsoid",
                   "--compare-adaptive-shaving", "--verify-only",
                   "--concurrency", "0"]
        completed = subprocess.run(command, capture_output=True,
                                   text=True, check=False)
        if completed.returncode:
            print(completed.stdout, file=sys.stderr)
            print(completed.stderr, file=sys.stderr)
            raise RuntimeError(f"benchmark failed with code {completed.returncode}")
        output = completed.stdout
    else:
        output = args.benchmark_output.read_text(encoding="utf-8")
    return 0 if audit(output, records) else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, OSError, RuntimeError) as error:
        print(f"AUDIT_ELLIPSOID_ERROR {error}", file=sys.stderr)
        sys.exit(2)
