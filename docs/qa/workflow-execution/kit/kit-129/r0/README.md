# KIT-129 R0, part B: the threshold of the factor-tied-to-groups warning (record of 2026-11-18)

Stage reached: **thresholds set on the setting set; check set NOT launched** (it does not fit under the ceiling;
the user decides). No verdict exists yet for any omic.

## The declaration (ticket section 3, frozen)

Setting set: runs 1 to 1000 of C0, C3, C4 (correlation 0) and C5 (0.3), runs 1 to 500 of stress S3 (0.6) and S11
(3 per group, 0.3), for arms 3 or 4, 5, 7 and 8, per simulator. Check set: runs 2001 to 2500 of the same cells, never
looked at while setting. Threshold: the value on the setting set that fires in at most 10% of the correlation-0 runs.
Accepted, per omic, when on the check set it fires in at most 10% of correlation-0 runs and at least 80% of the runs
at 0.3 or more. One attempt. If it does not separate, the GUI shows the value with no verdict.

Readings (R0 status note; in the header of `threshold.R` on `sim-k129-r0b`): one threshold per omic over the four arms
and the three correlation-0 cells (S4: both assays); a row without a W cannot fire and stays in the denominator; cells
classed by their declared correlation; fires = strictly above. Ceiling for part B: about 150 core-hours in all (user).

## What was run, and where

GitHub-hosted runners (4 cores per job), repo APAF-bioinformatics/MultiScholaR, branch `sim-k129-r0b`, workflow
`k124-sim.yml`, association mode, package tree `0b6fcb6d`. Core-hours = job duration x 4 cores (`jobs.tsv`).

| Run | What | Jobs | Attempts | Core-hours |
|---|---|---|---|---|
| 44201467446 | lifecycle of the association mode (C0, runs 1 to 8) | 5 + plan | all first | 3.16 |
| 44203265622 | the setting set (commit `33d5d7af`) | 28 + plan | all first | 129.62 |
| 44225068063 | lifecycle of `first_run` (C0, runs 2001 to 2004; commit `90a6685a`) | 4 + plan | all first | 2.49 |
| | **used** | | | **135.27** |

Setting set per simulator (job hours x 4): S1 43.0, S2 58.8, S3 11.9, S4 15.9. Every declared cell and run range is
present once in each simulator (S2's 1,000-run cells as runs 1 to 500 and 501 to 1000); no duplicate cell, run, arm
and assay; no crashed row. Rows per cell, arm, assay and status: `setting-status.tsv` (ok 90,787; no_factor 7,748;
refused 1,465). Raw tables are not committed (9 MB): `setting-assoc.sha256`; they are the run's artifacts
(`.github/k124-sim/assemble.sh 44203265622 r0b-setting <dir>`, kept 30 days by GitHub).

## Thresholds (`threshold.R set`, once per simulator; `threshold-S<n>.tsv`, shares in `setting-S<n>.tsv`)

| Omic | Threshold | Setting set: fires at correlation 0 | fires at 0.3 or more (needed on the check set: 80%) |
|---|---|---|---|
| S1 classic proteomics | 0.0939 | 10.0% | 75.4% |
| S2 DIA_LIMPA | 0.3347 | 10.0% | 40.6% |
| S3 transcriptomics | 0.3501 | 10.0% | 38.1% |
| S4 lipidomics and metabolomics | 0.3341 | 10.0% | 54.7% |

The setting-set shares are not the acceptance (that is the check set alone), but on them no omic reaches 80%.

Per arm, share fired (correlation 0 / 0.3 or more), from `setting-S<n>.tsv`:

| Omic | arm 3 or 4 (today's sweep) | arm 5 | arm 7 | arm 8 |
|---|---|---|---|---|
| S1 | 13.8% / 1.9% | 21.2% / 100.0% | 2.1% / 100.0% | 2.9% / 99.9% |
| S2 | 0.1% / 0.1% | 7.9% / 44.0% | 15.8% / 59.6% | 16.2% / 58.9% |
| S3 | 12.6% / 15.0% | 17.1% / 82.0% | 0.0% / 25.0% | 10.3% / 30.5% |
| S4 | 38.5% / 60.1% | 0.2% / 26.2% | 0.1% / 60.8% | 1.3% / 71.8% |

The statistic differs strongly between arms, and the one pooled threshold per omic carries that:
- S1 and S2, arm 4 (the sweep's W): the association stays near zero when the true factor follows the groups (median
  0.016 in S1 and 0.008 in S2, against 0.23 to 0.42 in the other arms). The warning cannot fire for that arm, and its
  quarter of the tied rows caps the pooled share at about 75% in S1.
- S1: the pooled threshold (0.094) is set by arms 4 and 5 at correlation 0; arms 5, 7 and 8 alone separate almost
  fully (at most 21% / at least 99.9%).
- S4, arm 4: fires in 38% of the correlation-0 rows (median association 0.26 with no tie), which pushes the pooled
  threshold up to 0.33, above the association most arms show at a true 0.3 (medians 0.24 to 0.37 in C5).
- S3: arm 5 reaches 0.34 median at correlation 0 in C0 (49% fire); arms 7 and 8 sit at the true value (0.29 to
  0.30 in C5 and S11), below the threshold 0.35, and fire in 0 to 14% there; they fire fully only at 0.6 (S3).
- S2: the true factor's own association is 0.21 (median) in the correlation-0 cells C0 and C4, since the correlation
  holds in expectation only; arms 7 and 8 follow it (22 to 25% fire there).
What the declared rule gives is what `threshold.R check` will say on the check set, per omic, pooled as above.

## The check set: not launched

Estimate from the setting run's own durations (`check-estimate.tsv`: seconds per run x 500 plus the job's setup, x 4
cores): **about 82 core-hours** (S1 27.9, S2 35.2, S3 8.4, S4 10.9). Used 135.3 + 82 = about 218, against about 150:
about 15 are left. Stopped as instructed.

The offset is ready: `first_run` in `.github/k124-sim/params.json` (optional, default 1: the matrix is identical to
the one before for the setting params, compared with jq). Lifecycle run 44225068063: the four simulators wrote runs
2001 to 2004 of C0 with seed 102001 to 102004 = 100000 + 10000 x 0 + run (`offset-lifecycle-seeds.tsv`). Those 4 runs
per simulator lie in the check set's range; only their cell, run and seed columns were read, after the thresholds were
written (`offset-lifecycle-assoc.sha256`). The check set recomputes them from the same seeds.

To launch (one push to `sim-k129-r0b`, params.json): label `r0b-check`, `runs_core` 500, `runs_stress` 500,
`first_run` 2001, cores 4, `timeout_minutes` 200, `assoc` true, no `chunk`, `arms` {}, the six cells per simulator.
Then, once: `assemble.sh <run> r0b-check <dir>` and, per simulator n,
`Rscript --vanilla threshold.R check S<n> <this directory> <dir>/s<n>-r0b-check/s<n>-assoc.tsv`
(it reads `threshold-S<n>.tsv` here and refuses a second time).
