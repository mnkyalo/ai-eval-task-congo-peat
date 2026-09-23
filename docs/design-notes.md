# Design notes and post-mortem

This file records how the task was designed and why I retired it before submitting it for review. I kept both on purpose. Being able to explain why a task fails is part of the same skill as being able to build one that passes.

## 1. Design goals

A good task for evaluating AI agents on scientific computing has to meet four conditions at once:

1. **It is real work.** The workflow is one that a palaeoecologist actually runs on a core: calibration, screening, an age-depth model, the exotic-marker method, influx and zonation. Nothing in it is a puzzle made up for the benchmark.
2. **Its difficulty comes from the science.** No information is withheld. Every correction has evidence in the files, and an expert who reads them carefully has what they need to act on it.
3. **It is graded objectively.** Each field has a tolerance measured from honest implementations, with no human judgement anywhere in the grading.
4. **A plausible wrong answer fails for the right reason.** Each shortcut should fail only in the section of the output its mistake belongs to.

The task meets conditions 1, 3 and 4, and the validation harness shows it. It does not meet condition 2 the way that matters in practice: the difficulty has to hold up against a strong solver, not only be present in principle.

## 2. How the tolerances were set

The procedure for every numeric field was the same. First I measured the spread of honest implementations. Then I measured the smallest error the field has to catch, and put the band in the gap between the two.

**Calibrated ages** (the table shows the largest deviation from the reference across all 15 ages):

| Variant | Largest deviation from reference | Classed as |
|---|---|---|
| 5-yr calibration grid | 3.1 yr | honest |
| 65% southern mixture | 7.6 yr | honest |
| 35% southern mixture | 10.1 yr | honest |
| SHCal20 only (100% southern) | 31.1 yr | honest (defensible) |
| IntCal20 only | 70.1 yr | wrong |

A band of ±45 yr leaves 1.45× headroom above the largest honest deviation and 1.56× headroom below the error it has to catch.

**Influx:** the largest honest deviation is 13.2%. The wrong answers deviate by 35–70% (superseded count), 50% (constant spike) and far more (contaminated dates, where a forced-monotonic model gives near-infinite rates). A band of ±25% leaves about 2× margin on both sides.

## 3. Why I retired it

### 3.1 The evidence gives the correction away

Take each piece of evidence and ask: *would this sentence or column exist, worded the same way, in a study where the problem did not exist?*

- `pretreatment = "AAA only; insufficient material to sieve"` exists only to flag the three bad dates.
- `notes = "slide cracked at 412 grains; superseded by PC-0316"` spells out the correction in words.
- `spike_and_volume.csv` has its own line in the README, and its only purpose is to be read level by level.
- The site description states the latitude and "twice-yearly passage of the ITCZ", which cues any reader who has seen the SHCal20 guidance.

None of these has a second, genuine reason to be in the dataset. A strong model reads the documentation, writes down every correction and then carries them out. What is left is a checklist, and checklists are easy for strong models. I later made this into a routine screen, the **plan-from-documentation test**. A fresh model gets only the instruction and the methods text, with no data, and is asked for its analysis plan. If that plan already contains the correction, the trap is not hidden.

### 3.2 Two traps are thinner than they look

- **Calibration curve.** Calibrating on IntCal20 alone fails on only **2 of 15** ages: 70.1 yr at 10 cm and 52.4 yr at 70 cm. The next-largest miss is 27.0 yr, well inside the band. The influx values move by only 11.9%, and the expansion sets do not change. The whole "hemispheric mixing" layer therefore rests on two numbers near the edge of the band.
- **Zonation.** A naive rule, "take the two largest level-to-level changes in the percentages", gives the same boundaries (30 and 60 cm) as constrained clustering. Both breaks are strong enough to see by eye. The trap does not discriminate, and the validation harness records this as a **NOT CAUGHT** row instead of hiding it. The task's own `verification_explanation` says an eyeballed-boundary answer fails. That holds for the particular eyeball placement tested when the task was built, but this simple jump rule passes, so that claim is weaker than it reads.

### 3.3 One trap announces itself

Each contaminated date is younger than the date above it. With all 15 dates kept, the linear age-depth model has negative accumulation rates, which is physically impossible. Any solver is forced to deal with the three dates, and the C:N values and pretreatment notes then tell them which three to drop. A trap the data points at this loudly works as a signpost.

## 4. What I carried forward

1. **Difficulty has to come from a mechanism that leaves no visible pattern** in a standard first look at the data (the usual plots, residuals and correlations). A signal that is visible cannot be hidden by tuning its size. The task just gets noisier without getting harder.
2. **Evidence must serve two purposes.** Every file that makes a correction possible should also exist for an honest, independent reason, the way a qPCR dilution series is there to bring the analyte into range.
3. **Measure the margin of every trap**, not only whether it fails. A trap that fails on 2 of 79 fields at 1.2–1.6× the band is fragile.
4. **Run the plan-from-documentation test before building the data.** It costs minutes, and it predicts a probe failure that would otherwise cost a full review round.

These rules shaped my later tasks, two of which have since passed both an easiness screen and a multi-model difficulty screen in independent review.

## 5. If I rebuilt it

The workflow is worth keeping. The traps are what need replacing. One direction: drop the self-announcing evidence and make the chronology problem one that only appears when two independent records are compared. For example, the core's pollen record could be compared with a charcoal or δ¹³C record whose own chronology constrains the contaminated interval. Then no single file flags the problem, and the solver has to notice that the two records disagree. That idea would go through the plan-from-documentation test before any data is generated.
