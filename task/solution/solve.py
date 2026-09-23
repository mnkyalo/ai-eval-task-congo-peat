#!/usr/bin/env python3
"""Oracle solution -- Core MBK-2 pollen influx task.

Derives every number from the solver-visible files only. It has no access to
the generator or to any record of the true parameters.

Method. The parts marked [STATED] are given to the solver in the task
instructions; everything else is the domain reasoning the task tests.

  1. [reasoning] The core sits at 0.68 degrees N and the ITCZ crosses it twice a
     year, so the peat fixed carbon from both hemispheres' air. Calibrating
     against IntCal20 alone is wrong for such a site; SHCal20's guidance for the
     equatorial zone is an equal mixture of the two atmospheric curves. Marine20
     does not apply -- the site is rain-fed, with no marine or reservoir carbon.
  2. [STATED] Each date is reported as the median of its calibrated density.
  3. [reasoning] Three samples are unfit for the age model. Their pre-treatment
     record shows they could not be sieved, their C:N sits far above the rest of
     the core (toward fresh root tissue rather than peat), and each dates
     younger than the level above it -- which peat cannot do. They carry modern
     rootlet carbon and are excluded from the age-depth model. They are still
     calibrated and reported: calibration is well defined for a sample that is
     simply unfit to date its depth.
  4. [STATED] The age-depth model is linear interpolation between the retained
     dated levels, so accumulation rate is constant within each interval.
  5. [reasoning] Level 55 cm was counted twice. The counting record marks
     PC-0271 as abandoned when the slide cracked at 412 grains and superseded by
     PC-0316. Only the recount is used.
  6. [reasoning] Pollen percentages are a closed sum, so a percentage rise can
     be produced entirely by another taxon's collapse. Absolute delivery comes
     from the exotic marker:
         concentration = (pollen sum / Lycopodium counted)
                         * (tablets * spores per tablet / sample volume)
         influx        = concentration * accumulation rate
     Tablet count and volume are read per level, not assumed.
  7. [STATED] Two zone boundaries, located by stratigraphically constrained
     cluster analysis on the pollen percentage data.
  8. [STATED] A taxon expanded across a boundary when its mean influx over the
     10 cm above is at least 1.5x its mean influx over the 10 cm below.

Reads only /app/data, the same files the agent gets, and writes
/app/results.json. Run via solution/solve.sh, or directly: python3 solve.py
"""
import csv
import json
import os
import numpy as np

TAXA = ["Raphia", "Uapaca", "Alchornea", "Elaeis", "Poaceae",
        "Cyperaceae", "Syzygium", "Macaranga", "Pandanus", "Musanga"]
WINDOW_CM = 10.0
EXPANSION_FACTOR = 1.5
N_ZONE_BOUNDARIES = 2
REPORT_DEPTHS = [27.0, 35.0, 55.0, 65.0, 69.0, 101.0]
GRID = np.arange(0.0, 8001.0, 1.0)


# ------------------------------------------------------------------- curves
def load_curve(path):
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
    return (np.interp(GRID, cal, np.asarray(age)[o]),
            np.interp(GRID, cal, np.asarray(sig)[o]))


def atmospheric_curve(env, p_south=0.5):
    """Equal mixture of the two hemispheric curves, for an ITCZ-crossing site."""
    mu_n, sd_n = load_curve(os.path.join(env, "curves", "intcal20.14c"))
    mu_s, sd_s = load_curve(os.path.join(env, "curves", "shcal20.14c"))
    pn = 1.0 - p_south
    return (pn * mu_n + p_south * mu_s,
            np.sqrt((pn * sd_n) ** 2 + (p_south * sd_s) ** 2))


def calibrated_median(c14_age, c14_err, mu, sd):
    tau2 = c14_err ** 2 + sd ** 2
    dens = np.exp(-((c14_age - mu) ** 2) / (2.0 * tau2)) / np.sqrt(tau2)
    dens = dens / dens.sum()
    return float(np.interp(0.5, np.cumsum(dens), GRID))


def read_csv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


# ---------------------------------------------------------------- zonation
def constrained_clustering(depths, pct):
    """Stratigraphically constrained agglomerative clustering.

    Only vertically adjacent groups may merge, so every cluster is a
    contiguous run of the core. The last merges are the strongest breaks; the
    two deepest-splitting boundaries are the zone boundaries.
    """
    groups = [[i] for i in range(len(depths))]
    merge_order = []
    while len(groups) > 1:
        best = None
        for i in range(len(groups) - 1):
            a = pct[groups[i]].mean(axis=0)
            b = pct[groups[i + 1]].mean(axis=0)
            n1, n2 = len(groups[i]), len(groups[i + 1])
            cost = float(np.sqrt(((a - b) ** 2).sum())) * n1 * n2 / (n1 + n2)
            if best is None or cost < best[0]:
                best = (cost, i)
        _, i = best
        merge_order.append((depths[groups[i][-1]] + depths[groups[i + 1][0]]) / 2.0)
        groups[i:i + 2] = [groups[i] + groups[i + 1]]
    return sorted(merge_order[::-1][:N_ZONE_BOUNDARIES])


# -------------------------------------------------------------------- main
def solve(env):
    mu, sd = atmospheric_curve(env)

    dates = read_csv(os.path.join(env, "radiocarbon_dates.csv"))
    calibrated = {float(r["depth_cm"]):
                  calibrated_median(float(r["c14_age_bp"]),
                                    float(r["c14_error"]), mu, sd)
                  for r in dates}

    # Unfit for the age model: could not be sieved AND C:N far above the core
    # background. Both signals point at the same samples, and each of them also
    # dates younger than the level above it.
    cn = np.array([float(r["c_to_n"]) for r in dates])
    background = float(np.median(cn))
    excluded = sorted(
        float(r["depth_cm"]) for r, v in zip(dates, cn)
        if "insufficient material to sieve" in r["pretreatment"]
        and v > 1.5 * background)

    kept = sorted(d for d in calibrated if d not in excluded)
    dep = np.array(kept)
    age = np.array([calibrated[d] for d in kept])

    def accumulation_rate(depth):
        i = int(np.clip(np.searchsorted(dep, depth, side="right") - 1,
                        0, len(dep) - 2))
        return (dep[i + 1] - dep[i]) / (age[i + 1] - age[i])

    # Superseded counts are dropped; the record names which row replaces which.
    counts = [r for r in read_csv(os.path.join(env, "pollen_counts.csv"))
              if "superseded by" not in r["notes"]]
    spike = {float(r["depth_cm"]): r
             for r in read_csv(os.path.join(env, "spike_and_volume.csv"))}

    influx, percent = {}, {}
    for r in counts:
        d = float(r["depth_cm"])
        taxon_counts = {t: float(r[t]) for t in TAXA}
        pollen_sum = sum(taxon_counts.values())
        s = spike[d]
        added = float(s["lycopodium_tablets"]) * float(s["spores_per_tablet"])
        volume = float(s["sample_volume_cm3"])
        concentration = (pollen_sum / float(r["lycopodium_count"])) * (added / volume)
        acc = accumulation_rate(d)
        influx[d] = {t: concentration * taxon_counts[t] / pollen_sum * acc
                     for t in TAXA}
        percent[d] = {t: 100.0 * taxon_counts[t] / pollen_sum for t in TAXA}

    levels = np.array(sorted(percent))
    pct_matrix = np.array([[percent[d][t] for t in TAXA] for d in levels])
    boundaries = constrained_clustering(levels, pct_matrix)

    expansions = {}
    for b in boundaries:
        above = [influx[d] for d in influx if b - WINDOW_CM < d < b]
        below = [influx[d] for d in influx if b < d < b + WINDOW_CM]
        expansions[f"{b:.0f}"] = sorted(
            t for t in TAXA
            if np.mean([x[t] for x in above]) >=
               EXPANSION_FACTOR * np.mean([x[t] for x in below]))

    return {
        "chronology_cal_bp": {f"{d:.0f}": round(calibrated[d], 1)
                              for d in sorted(calibrated)},
        "excluded_from_age_model": [f"{d:.0f}" for d in excluded],
        "zone_boundaries_cm": [f"{b:.0f}" for b in boundaries],
        "influx_grains_per_cm2_per_yr": {
            f"{d:.0f}": {t: round(influx[d][t], 3) for t in TAXA}
            for d in REPORT_DEPTHS},
        "expansions": expansions,
    }


if __name__ == "__main__":
    result = solve("/app/data")
    with open("/app/results.json", "w") as fh:
        json.dump(result, fh, indent=1, sort_keys=True)
    print("wrote /app/results.json")
