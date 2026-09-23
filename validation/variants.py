"""Parameterised solver used to build honest variants and plausible wrong answers.

Each wrong answer switches exactly ONE analytical decision away from the
reference, so a failure can be attributed to the mistake it models. Honest
variants change decisions a competent analyst could reasonably make
differently (mixing proportion, grid resolution) and must still pass.

Reads only the solver-visible data, exactly like the oracle.
"""
import csv
import os

import numpy as np

TAXA = ["Raphia", "Uapaca", "Alchornea", "Elaeis", "Poaceae",
        "Cyperaceae", "Syzygium", "Macaranga", "Pandanus", "Musanga"]
REPORT_DEPTHS = [27.0, 35.0, 55.0, 65.0, 69.0, 101.0]


def _read_csv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def _load_curve(path, grid):
    cal, age, sig = [], [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.strip().split(",")
            if len(p) < 3:
                continue
            cal.append(float(p[0])); age.append(float(p[1])); sig.append(float(p[2]))
    o = np.argsort(cal)
    cal = np.asarray(cal)[o]
    return (np.interp(grid, cal, np.asarray(age)[o]),
            np.interp(grid, cal, np.asarray(sig)[o]))


def _curve(env, kind, p_south, grid):
    c = os.path.join(env, "curves")
    if kind == "intcal":
        return _load_curve(os.path.join(c, "intcal20.14c"), grid)
    if kind == "marine":
        return _load_curve(os.path.join(c, "marine20.14c"), grid)
    mu_n, sd_n = _load_curve(os.path.join(c, "intcal20.14c"), grid)
    mu_s, sd_s = _load_curve(os.path.join(c, "shcal20.14c"), grid)
    pn = 1.0 - p_south
    return pn * mu_n + p_south * mu_s, np.sqrt((pn * sd_n) ** 2 + (p_south * sd_s) ** 2)


def _median(age14, err, mu, sd, grid):
    tau2 = err ** 2 + sd ** 2
    dens = np.exp(-((age14 - mu) ** 2) / (2 * tau2)) / np.sqrt(tau2)
    cdf = np.cumsum(dens / dens.sum())
    return float(np.interp(0.5, cdf, grid))


def _constrained_clustering(depths, pct, n=2):
    groups = [[i] for i in range(len(depths))]
    merges = []
    while len(groups) > 1:
        best = None
        for i in range(len(groups) - 1):
            a, b = pct[groups[i]].mean(0), pct[groups[i + 1]].mean(0)
            n1, n2 = len(groups[i]), len(groups[i + 1])
            cost = float(np.sqrt(((a - b) ** 2).sum())) * n1 * n2 / (n1 + n2)
            if best is None or cost < best[0]:
                best = (cost, i)
        i = best[1]
        merges.append((depths[groups[i][-1]] + depths[groups[i + 1][0]]) / 2)
        groups[i:i + 2] = [groups[i] + groups[i + 1]]
    return sorted(merges[::-1][:n])


def _largest_jumps(depths, pct, n=2):
    """'By eye': the two biggest level-to-level changes, no clustering."""
    jumps = np.sqrt((np.diff(pct, axis=0) ** 2).sum(1))
    idx = np.argsort(jumps)[::-1][:n]
    return sorted((depths[i] + depths[i + 1]) / 2 for i in idx)


def solve(env, curve="mixed", p_south=0.5, grid_step=1.0,
          screen_dates=True, use_superseded=False, constant_spike=False,
          expansions_on="influx", zoning="clustering"):
    grid = np.arange(0.0, 8001.0, grid_step)
    mu, sd = _curve(env, curve, p_south, grid)

    dates = _read_csv(os.path.join(env, "radiocarbon_dates.csv"))
    cal = {float(r["depth_cm"]): _median(float(r["c14_age_bp"]), float(r["c14_error"]),
                                         mu, sd, grid) for r in dates}
    cn = np.array([float(r["c_to_n"]) for r in dates])
    excluded = []
    if screen_dates:
        excluded = sorted(float(r["depth_cm"]) for r, v in zip(dates, cn)
                          if "insufficient material to sieve" in r["pretreatment"]
                          and v > 1.5 * np.median(cn))
    kept = sorted(d for d in cal if d not in excluded)
    dep, age = np.array(kept), np.array([cal[d] for d in kept])
    if np.any(np.diff(age) <= 0):
        # An analyst who keeps every date still needs a usable model; the
        # usual patch is to force monotonicity (1 yr minimum step).
        age = np.maximum.accumulate(age) + np.arange(len(age)) * 1.0

    def acc(depth):
        i = int(np.clip(np.searchsorted(dep, depth, side="right") - 1, 0, len(dep) - 2))
        return (dep[i + 1] - dep[i]) / (age[i + 1] - age[i])

    rows = _read_csv(os.path.join(env, "pollen_counts.csv"))
    if use_superseded:
        rows = [r for r in rows if "supersedes" not in r["notes"]]
    else:
        rows = [r for r in rows if "superseded by" not in r["notes"]]
    spike = {float(r["depth_cm"]): r for r in _read_csv(os.path.join(env, "spike_and_volume.csv"))}

    influx, pct = {}, {}
    for r in rows:
        d = float(r["depth_cm"])
        c = {t: float(r[t]) for t in TAXA}
        s = sum(c.values())
        sp = spike[d]
        tablets = 2.0 if constant_spike else float(sp["lycopodium_tablets"])
        vol = 1.0 if constant_spike else float(sp["sample_volume_cm3"])
        conc = (s / float(r["lycopodium_count"])) * tablets * float(sp["spores_per_tablet"]) / vol
        influx[d] = {t: conc * c[t] / s * acc(d) for t in TAXA}
        pct[d] = {t: 100 * c[t] / s for t in TAXA}

    levels = np.array(sorted(pct))
    P = np.array([[pct[d][t] for t in TAXA] for d in levels])
    bounds = (_constrained_clustering if zoning == "clustering" else _largest_jumps)(levels, P)

    measure = influx if expansions_on == "influx" else pct
    exp = {}
    for b in bounds:
        above = [measure[d] for d in measure if b - 10 < d < b]
        below = [measure[d] for d in measure if b < d < b + 10]
        exp[f"{b:.0f}"] = sorted(t for t in TAXA
                                 if np.mean([x[t] for x in above]) >= 1.5 * np.mean([x[t] for x in below]))

    return {
        "chronology_cal_bp": {f"{d:.0f}": round(cal[d], 1) for d in sorted(cal)},
        "excluded_from_age_model": [f"{d:.0f}" for d in excluded],
        "zone_boundaries_cm": [f"{b:.0f}" for b in bounds],
        "influx_grains_per_cm2_per_yr": {f"{d:.0f}": {t: round(influx[d][t], 3) for t in TAXA}
                                         for d in REPORT_DEPTHS},
        "expansions": exp,
    }
