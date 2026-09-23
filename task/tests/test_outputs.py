"""Verifier for the Core MBK-2 pollen influx task.

The answer is DERIVED here, not looked up. The raw data is sealed into the
verifier image at /tests/data and `reference.py` recomputes the whole analysis
at grading time, so nothing the agent writes can influence what it is compared
against. `answer_key.json` is kept only as an integrity check on the verifier
itself: if the derivation ever stops reproducing it, the run fails loudly
rather than silently grading against something new.

Three layers, 81 checks:
  * the frozen key still matches a fresh derivation (verifier integrity)
  * the submission matches that derivation, field by field
  * the submission is internally consistent -- the influx values must follow
    from the chronology the same submission reports, which a fabricated or
    hedged answer cannot satisfy
"""
import os

import numpy as np
import pytest

import grading_lib as G
import reference

SEALED_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

# Derived fresh, inside the verifier, from the sealed copy of the raw data.
DERIVED = reference.solve(SEALED_DATA)
KEY = G.load_key()


@pytest.fixture(scope="session")
def submission():
    return G.load_submission()


# ---------------------------------------------------- verifier integrity
def test_frozen_key_matches_fresh_derivation():
    """The stored key and a fresh derivation must agree, or the verifier is
    grading against something other than what the data implies."""
    assert DERIVED == KEY, (
        "the sealed derivation no longer reproduces answer_key.json; "
        "the verifier image is inconsistent and the run cannot be trusted")


def test_submission_present_and_well_formed(submission):
    assert set(G.SECTIONS).issubset(submission)


# ------------------------------------------------------- field-by-field
@pytest.mark.parametrize("depth,expected", G.chronology_cases(KEY),
                         ids=[f"age_{d}cm" for d, _ in G.chronology_cases(KEY)])
def test_calibrated_age(submission, depth, expected):
    G.check_age(submission, depth, DERIVED["chronology_cal_bp"][depth])


def test_excluded_from_age_model(submission):
    G.check_depth_set(submission, DERIVED, "excluded_from_age_model")


def test_zone_boundaries(submission):
    G.check_depth_set(submission, DERIVED, "zone_boundaries_cm")


@pytest.mark.parametrize("depth,taxon,expected", G.influx_cases(KEY),
                         ids=[f"influx_{d}cm_{t}" for d, t, _ in G.influx_cases(KEY)])
def test_influx(submission, depth, taxon, expected):
    G.check_influx(submission, depth, taxon,
                   DERIVED["influx_grains_per_cm2_per_yr"][depth][taxon])


@pytest.mark.parametrize("boundary,expected", G.expansion_cases(KEY),
                         ids=[f"expansions_{b}cm" for b, _ in G.expansion_cases(KEY)])
def test_expansions(submission, boundary, expected):
    G.check_expansions(submission, boundary, DERIVED["expansions"][boundary])


# --------------------------------------------------- internal consistency
def test_influx_follows_from_the_reported_chronology(submission):
    """Real execution evidence, not a lookup.

    Recompute each reported influx from the raw counts, the spike record, and
    the age-depth model implied by THIS submission's own chronology and
    excluded-date list. An answer that was derived by doing the work satisfies
    this automatically; one assembled field by field, or hedged, does not.
    """
    ages = {G.norm_depth(k): float(v)
            for k, v in submission["chronology_cal_bp"].items()}
    excluded = {G.norm_depth(x) for x in submission["excluded_from_age_model"]}
    kept = sorted((float(d) for d in ages if d not in excluded), key=float)
    assert len(kept) >= 2, "fewer than two dated levels left in the age model"

    dep = np.array(kept)
    age = np.array([ages[G.norm_depth(d)] for d in kept])
    assert np.all(np.diff(age) > 0), (
        "the reported chronology is not monotonic once the reported exclusions "
        "are removed; no age-depth model can be built from it")

    def acc(depth):
        i = int(np.clip(np.searchsorted(dep, depth, side="right") - 1,
                        0, len(dep) - 2))
        return (dep[i + 1] - dep[i]) / (age[i + 1] - age[i])

    counts = [r for r in reference.read_csv(
        os.path.join(SEALED_DATA, "pollen_counts.csv"))
        if "superseded by" not in r["notes"]]
    by_depth = {float(r["depth_cm"]): r for r in counts}
    spike = {float(r["depth_cm"]): r for r in reference.read_csv(
        os.path.join(SEALED_DATA, "spike_and_volume.csv"))}

    got = {G.norm_depth(k): {G.norm_taxon(t): float(v) for t, v in row.items()}
           for k, row in submission["influx_grains_per_cm2_per_yr"].items()}

    worst = 0.0
    for depth in reference.REPORT_DEPTHS:
        r, s = by_depth[depth], spike[depth]
        taxon_counts = {t: float(r[t]) for t in reference.TAXA}
        pollen_sum = sum(taxon_counts.values())
        conc = (pollen_sum / float(r["lycopodium_count"])) * (
            float(s["lycopodium_tablets"]) * float(s["spores_per_tablet"])
            / float(s["sample_volume_cm3"]))
        rate = acc(depth)
        for taxon in reference.TAXA:
            implied = conc * taxon_counts[taxon] / pollen_sum * rate
            reported = got[G.norm_depth(depth)][G.norm_taxon(taxon)]
            rel = abs(reported - implied) / implied if implied else 0.0
            worst = max(worst, rel)
            assert rel <= 0.05, (
                f"{depth:.0f} cm {taxon}: reported influx {reported:.3f} does "
                f"not follow from this submission's own chronology, which "
                f"implies {implied:.3f} ({100 * rel:.1f}% apart)")
