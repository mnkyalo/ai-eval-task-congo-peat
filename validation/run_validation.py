#!/usr/bin/env python3
"""Run the oracle, honest variants and plausible wrong answers through the
task's REAL verifier code (task/tests/test_outputs.py), without Docker.

A task's grader is only as good as the mistakes it demonstrably catches, so
this script is the evidence behind the README's claims:

  * the oracle and every honest variant must score 1 (all checks pass)
  * every wrong answer must score 0, and fail on the section its mistake
    belongs to
  * a shortcut the grader is known NOT to catch is listed as such and
    asserted to pass, so a weakness is recorded rather than hidden

Usage:  python3 validation/run_validation.py        (needs only numpy)
Exit status is non-zero if any expectation is violated, so CI can gate on it.
"""
import json
import os
import sys
import tempfile
import types

HERE = os.path.dirname(os.path.abspath(__file__))
TASK = os.path.join(os.path.dirname(HERE), "task")
DATA = os.path.join(TASK, "environment", "data")
sys.path[:0] = [os.path.join(TASK, "tests"), os.path.join(TASK, "solution"), HERE]

# The verifier is written for pytest. If pytest is not installed, a two-line
# stand-in lets us import the very same test functions and call them directly.
try:
    import pytest  # noqa: F401
except ImportError:
    stub = types.ModuleType("pytest")
    stub.fixture = lambda *a, **k: (lambda f: f)
    stub.mark = types.SimpleNamespace(parametrize=lambda *a, **k: (lambda f: f))
    sys.modules["pytest"] = stub

import grading_lib as G          # noqa: E402
import solve as oracle           # noqa: E402
import test_outputs as T         # noqa: E402
import variants as V             # noqa: E402

# (name, builder, expected reward, what it models)
CASES = [
    ("oracle (solution/solve.py)", lambda: oracle.solve(DATA), 1,
     "reference derivation"),
    ("honest: 35% southern mixture", lambda: V.solve(DATA, p_south=0.35), 1,
     "defensible alternative mixing proportion"),
    ("honest: 65% southern mixture", lambda: V.solve(DATA, p_south=0.65), 1,
     "defensible alternative mixing proportion"),
    ("honest: 5-yr calibration grid", lambda: V.solve(DATA, grid_step=5.0), 1,
     "coarser numerical integration"),
    ("wrong: IntCal20 only", lambda: V.solve(DATA, curve="intcal"), 0,
     "ignores the ITCZ / hemispheric mixing"),
    ("wrong: Marine20", lambda: V.solve(DATA, curve="marine"), 0,
     "picks the wrong curve from the library"),
    ("wrong: keep rootlet-contaminated dates", lambda: V.solve(DATA, screen_dates=False), 0,
     "no contamination screen"),
    ("wrong: use superseded count at 55 cm", lambda: V.solve(DATA, use_superseded=True), 0,
     "ignores the counting record"),
    ("wrong: assume 2 tablets in 1 cm3", lambda: V.solve(DATA, constant_spike=True), 0,
     "ignores the per-level preparation record"),
    ("wrong: expansions from percentages", lambda: V.solve(DATA, expansions_on="percent"), 0,
     "closed-sum artefact"),
    # Documented weakness, kept here on purpose: this shortcut was meant to be
    # a trap, but the two breaks are strong enough that a naive largest-jump
    # rule finds the same boundaries. The grader cannot catch it, and this
    # row asserts that fact so the README's claim stays honest.
    ("NOT CAUGHT: boundaries from largest jumps", lambda: V.solve(DATA, zoning="jumps"), 1,
     "zonation by eye -- a trap that does not bite"),
    ("wrong: every shortcut at once",
     lambda: V.solve(DATA, curve="intcal", screen_dates=False, use_superseded=True,
                     constant_spike=True, expansions_on="percent", zoning="jumps"), 0,
     "careless end-to-end run"),
]


def checks():
    """Every check the verifier runs, as (section, callable(submission))."""
    key = T.KEY
    out = [("integrity", lambda s: T.test_frozen_key_matches_fresh_derivation()),
           ("schema", T.test_submission_present_and_well_formed)]
    out += [("chronology", lambda s, d=d, e=e: T.test_calibrated_age(s, d, e))
            for d, e in G.chronology_cases(key)]
    out += [("excluded", T.test_excluded_from_age_model),
            ("boundaries", T.test_zone_boundaries)]
    out += [("influx", lambda s, d=d, t=t, e=e: T.test_influx(s, d, t, e))
            for d, t, e in G.influx_cases(key)]
    out += [("expansions", lambda s, b=b, e=e: T.test_expansions(s, b, e))
            for b, e in G.expansion_cases(key)]
    out += [("consistency", T.test_influx_follows_from_the_reported_chronology)]
    return out


def grade(result):
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
        json.dump(result, fh)
    sub = G.load_submission(fh.name)      # same parser the verifier uses
    os.unlink(fh.name)
    failed = {}
    all_checks = checks()
    for section, fn in all_checks:
        try:
            fn(sub)
        except AssertionError:
            failed[section] = failed.get(section, 0) + 1
    return len(all_checks), failed


def main():
    ok = True
    print(f"{'submission':42s} {'reward':>6s} {'expect':>6s}  failed checks")
    print("-" * 100)
    for name, build, expected, _ in CASES:
        n, failed = grade(build())
        reward = int(not failed)
        detail = ", ".join(f"{k} {v}" for k, v in failed.items()) or f"all {n} passed"
        flag = "" if reward == expected else "   <-- UNEXPECTED"
        ok &= reward == expected
        print(f"{name:42s} {reward:>6d} {expected:>6d}  {detail}{flag}")
    print("-" * 100)
    print("ALL EXPECTATIONS MET" if ok else "EXPECTATION VIOLATED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
