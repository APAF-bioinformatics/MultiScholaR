# KIT-129 R0, part B: amendment of the warning's rule before the blind check (2026-11-18)

Authority: the user's decision of 2026-11-18 ("C and yes"), decision M23 in
`docs/qa/workflow-execution/kit/kit-124/design-review-20261118.md` section 6 (origin/performance-v2 `0eb89342`,
committed 2026-11-18 13:57 +11:00). This file is committed before the per-design thresholds are computed and before
the check set is launched; the commit times of this branch and of `sim-k129-r0b` give the order.

## State of the check set when this was written

The check set (runs 2001 to 2500 of C0, C3, C4, C5, S3, S11, simulators S1 to S4) had not been run, with one
exception: the lifecycle of the workflow's `first_run` offset (hosted-runner run 44225068063) computed runs 2001 to
2004 of C0 in each of the four simulators (16 runs). Of those tables only the cell, run and seed columns were read
(`offset-lifecycle-seeds.tsv`; sha256 of the tables in `offset-lifecycle-assoc.sha256`); no statistic of them was
read, listed or summarised. Those 16 runs are part of the check set as declared and stay in it: the check run
recomputes them from the same seeds.

## The original rule (ticket section 3; unchanged, still scored)

One threshold per omic, pooled over the four designs (arms 3 or 4, 5, 7, 8) and the three correlation-0 cells: the
value on the setting set that fires in at most 10% of the correlation-0 runs. Accepted, per omic, when on the check
set it fires in at most 10% of correlation-0 runs and at least 80% of the runs at 0.3 or more. Its thresholds are in
`threshold-S<n>.tsv` (set once, before this amendment): 0.0939, 0.3347, 0.3501, 0.3341.

## What the setting set showed, and why the rule is amended

Pooled, the thresholds fire on the setting set in 75.4% (S1), 40.6% (S2), 38.1% (S3) and 54.7% (S4) of the runs at
0.3 or more. The statistic behaves differently by design, so one pooled value serves no design well. Share fired
under the pooled threshold, correlation 0 / 0.3 or more (`setting-S<n>.tsv`):

| Omic | arm 3 or 4 (today's sweep) | arm 5 | arm 7 | arm 8 |
|---|---|---|---|---|
| S1 classic proteomics | 13.8% / 1.9% | 21.2% / 100.0% | 2.1% / 100.0% | 2.9% / 99.9% |
| S2 DIA_LIMPA | 0.1% / 0.1% | 7.9% / 44.0% | 15.8% / 59.6% | 16.2% / 58.9% |
| S3 transcriptomics | 12.6% / 15.0% | 17.1% / 82.0% | 0.0% / 25.0% | 10.3% / 30.5% |
| S4 lipidomics and metabolomics | 38.5% / 60.1% | 0.2% / 26.2% | 0.1% / 60.8% | 1.3% / 71.8% |

The sweep's factors (arm 4) stay near zero in S1 and S2 whatever the truth (median association 0.016 and 0.008 when
the true factor follows the groups), and are high without any tie in S4 (38.5% fire at correlation 0); either way
they move the pooled threshold away from where the other designs separate.

## The amended rule

- One threshold per design (arm) and omic: the value on the setting set that fires in at most 10% of that design's
  correlation-0 runs (C0, C3, C4).
- Conventions as in the original: fires = the largest factor association is strictly above the threshold; the
  threshold is the smallest observed value with at most 10% of the design's correlation-0 rows above it (0 when fewer
  rows than that have a W); a row without a W cannot fire and stays in the denominator; a crashed row is left out and
  reported; cells are classed by their declared correlation (C0, C3, C4: 0; C5, S3, S11: 0.3 or more); in S4 the two
  assays are taken together; a set must be complete as declared.
- Accepted, per design and omic, when on the check set the threshold fires in at most 10% of that design's
  correlation-0 runs and at least 80% of its runs at 0.3 or more.
- Both rules are scored on the check set, each once, and both are reported.
- One attempt: no further tuning after the check set is read, whatever it shows. A design that does not meet the
  marks gets "no verdict".
- The sets are as declared (unchanged). Ceiling for part B: about 220 core-hours on hosted runners (135.3 used when
  this was written; the check set is estimated at about 82).

Scoring: `threshold-by-design.R` beside `threshold.R` on `sim-k129-r0b` (the original script is not edited), with
the same separation: `set` refuses check-set rows and a second set; `check` refuses without thresholds and a second
run. Its per-design thresholds and setting-set shares are committed here before the check set is launched.
