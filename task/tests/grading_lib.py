"""Grading rules for the Core MBK-2 pollen influx task.

Same contract as the pipeline's standalone grader: five sections, per-field
rules, no human judgment. Tolerances were measured, not chosen -- each comment
records the honest spread it must absorb and the error it must catch.
"""
import json
import os

ANSWER_KEY = os.path.join(os.path.dirname(os.path.abspath(__file__)), "answer_key.json")
SUBMISSION = "/app/results.json"

# Calibrated ages, absolute, cal yr. Honest implementations -- mixing
# proportions 0.35 to 1.0 southern, grid resolution, sigma-mixing rule, CDF
# integration, and a third-party calibration engine -- span 31.1 yr.
# Calibrating on IntCal20 alone deviates 70.1 yr.
AGE_TOL_YR = 45.0

# Pollen influx, relative. Honest variants span 13.2%. Assuming a constant
# Lycopodium spike is 50% out, using the superseded count 35.3-70.5%, and
# keeping the contaminated dates in the age model 225%.
INFLUX_TOL_REL = 0.25

SECTIONS = ("chronology_cal_bp", "excluded_from_age_model", "zone_boundaries_cm",
            "influx_grains_per_cm2_per_yr", "expansions")


def norm_depth(k):
    return f"{float(k):.0f}"


def norm_taxon(k):
    return str(k).strip().lower()


def load_key():
    with open(ANSWER_KEY) as fh:
        return json.load(fh)


def load_submission(path=SUBMISSION):
    assert os.path.exists(path), f"no submission written to {path}"
    with open(path) as fh:
        raw = fh.read()
    try:
        sub = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise AssertionError(f"{path} is not valid JSON: {exc}") from None
    assert isinstance(sub, dict), f"top level of {path} must be a JSON object"
    for section in SECTIONS:
        assert section in sub, f"{path} is missing the required section '{section}'"
    return sub


def as_number(value, where):
    try:
        return float(value)
    except (TypeError, ValueError):
        raise AssertionError(f"{where}: expected a number, got {value!r}") from None


def chronology_cases(key):
    return sorted(key["chronology_cal_bp"].items(), key=lambda kv: float(kv[0]))


def influx_cases(key):
    out = []
    for depth, row in sorted(key["influx_grains_per_cm2_per_yr"].items(),
                             key=lambda kv: float(kv[0])):
        for taxon, expected in sorted(row.items()):
            out.append((depth, taxon, expected))
    return out


def expansion_cases(key):
    return sorted(key["expansions"].items(), key=lambda kv: float(kv[0]))


def check_age(sub, depth, expected):
    got = {norm_depth(k): v for k, v in sub["chronology_cal_bp"].items()}
    d = norm_depth(depth)
    assert d in got, f"chronology_cal_bp is missing depth {d} cm"
    value = as_number(got[d], f"chronology_cal_bp[{d}]")
    dev = abs(value - float(expected))
    assert dev <= AGE_TOL_YR, (
        f"{d} cm: {value:.1f} cal BP vs {float(expected):.1f}, "
        f"off by {dev:.1f} yr (tolerance {AGE_TOL_YR:.0f} yr)")


def check_depth_set(sub, key, field):
    value = sub[field]
    assert isinstance(value, list), f"{field} must be a list of depths in cm"
    have = {norm_depth(x) for x in value}
    want = {norm_depth(x) for x in key[field]}
    spurious, missing = sorted(have - want), sorted(want - have)
    assert have == want, (
        f"{field}: spurious={spurious or '-'} missing={missing or '-'}; "
        f"expected {sorted(want, key=float)}")


def check_influx(sub, depth, taxon, expected):
    got = {norm_depth(k): {norm_taxon(t): v for t, v in row.items()}
           for k, row in sub["influx_grains_per_cm2_per_yr"].items()}
    d, t = norm_depth(depth), norm_taxon(taxon)
    assert d in got, f"influx_grains_per_cm2_per_yr is missing depth {d} cm"
    assert t in got[d], f"influx at {d} cm is missing taxon {taxon}"
    value = as_number(got[d][t], f"influx[{d}][{taxon}]")
    exp = float(expected)
    rel = abs(value - exp) / exp if exp else (0.0 if value == 0 else float("inf"))
    assert rel <= INFLUX_TOL_REL, (
        f"{d} cm {taxon}: {value:.3f} vs {exp:.3f}, off by {100 * rel:.1f}% "
        f"(tolerance {100 * INFLUX_TOL_REL:.0f}%)")


def check_expansions(sub, boundary, expected):
    got = {norm_depth(k): v for k, v in sub["expansions"].items()}
    b = norm_depth(boundary)
    assert b in got, f"expansions is missing boundary {b} cm -- was it located?"
    value = got[b]
    assert isinstance(value, list), f"expansions[{b}] must be a list of taxon names"
    have = {norm_taxon(x) for x in value}
    want = {norm_taxon(x) for x in expected}
    spurious, missing = sorted(have - want), sorted(want - have)
    assert have == want, (
        f"{b} cm: spurious={spurious or '-'} missing={missing or '-'}; "
        f"expected {sorted(expected)}")
