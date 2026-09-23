# Core MBK-2, chronology, pollen accumulation and swamp forest change

A peat core from the Cuvette Centrale has been dated and counted. Everything that exists for it is in `/app/data`: fifteen AMS radiocarbon determinations with their laboratory record, seventy counted pollen levels, the sample preparation record, a site description, and a calibration curve library. See `/app/data/README.md` for the columns in each file.

Put the record on a calendar timescale, quantify how much pollen of each taxon was being delivered to the peat surface through time, locate where the assemblage changes, and say which taxa expanded at those changes.

## What to produce

Write your answer to `/app/results.json`, with exactly these five sections:

```
{
  "chronology_cal_bp": {
    "<dated depth in cm>": <calibrated age in cal BP, a number>
  },
  "excluded_from_age_model": ["<depth in cm>"],
  "zone_boundaries_cm": ["<depth in cm>", "<depth in cm>"],
  "influx_grains_per_cm2_per_yr": {
    "<depth in cm>": { "<taxon name>": <influx, a number> }
  },
  "expansions": {
    "<boundary depth in cm>": ["<taxon name>", "<taxon name>"]
  }
}
```

The angle brackets mark placeholders and show the shape only. The required contents are:

1. `chronology_cal_bp` -- the calibrated age of each of the fifteen dated depths, in calendar years before present, keyed by depth in cm.
2. `excluded_from_age_model` -- any dated depths you left out of the age-depth model. Give an empty list if you left none out.
3. `zone_boundaries_cm` -- the two zone boundaries you locate.
4. `influx_grains_per_cm2_per_yr` -- pollen influx for all ten taxa at six levels: 27, 35, 55, 65, 69 and 101 cm. Key by depth, then by taxon, using the taxon names as they appear in `pollen_counts.csv`.
5. `expansions` -- for each of the two zone boundaries, keyed by its depth in cm, the list of taxa that expanded across it.

## Conventions to follow

These are fixed, so that your numbers and the reference are comparable. Everything else about the analysis is yours to decide.

* Report each radiocarbon date as the median of its calibrated probability density, not an intercept or a range midpoint. Report all fifteen, whether or not you use them in the age-depth model.
* Build the age-depth model by linear interpolation between dated levels, so that the sediment accumulation rate is constant within each interval.
* The assemblage defines two zone boundaries. Locate them by stratigraphically constrained cluster analysis on the pollen percentage data, taking the two strongest breaks.
* A taxon counts as having expanded across a boundary when its mean influx over the 10 cm above the boundary is at least 1.5 times its mean influx over the 10 cm below it.

Depth keys may be written as `30` or `30.0`, taxon capitalisation and list order do not matter.

You have 10800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
