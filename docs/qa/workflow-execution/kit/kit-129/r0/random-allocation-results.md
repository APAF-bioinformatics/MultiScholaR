# KIT-129 R0 part B, supplement: the unrelated cells with the factor allocated at random, read once (2026-11-19)

Declaration: `random-allocation-declaration.md` (decision M37; pushed alone at `29093154` before the launch at
`62c052bc`; order in `random-allocation-pushes.tsv`, GitHub's activity log). **Descriptive: no threshold, verdict or
earlier record is changed by this result.**

## The run

| Run | What | Jobs | Core-hours |
|---|---|---|---|
| 44419142081 | lifecycle (C0, C0R, runs 1 to 4) | 6, all first attempt | 3.74 |
| 44421002832 | scoring lifecycle | 1 | 0.04 |
| 44421250021 | declared run: C0R, C4R, S1R of S1, S3, S4, runs 3001 to 3500 | 9, all first attempt, no rerun | 21.44 |
| 44431538031 | declared scoring, once (`random-allocation.R`, sha256 `281b70d1...`) | 1 | 0.04 |
| | **total** (budget 40) | | **25.26** |

Every cell and run is present once, with no crashed or refused rows (`random-allocation-status.tsv`). Two S1R runs had
E = 0 (arms 7 and 8 `no_factor`: no W, cannot fire, kept in the denominator). The raw tables are not committed:
`random-allocation-inputs.sha256` (the run's artifacts, kept 30 days).

## Result 1: how often the check fires with the factor allocated at random

The true association is |correlation| for two groups, one value per run. The share fired is the share of rows above
the validated limit (S4: both assays). C0R+C4R are the calibration's own unrelated designs; S1R has 3 per group.

| Omic | Design | Limit (set at) | True assoc. median / P90 | Estimated median / P90 | Fires: C0R+C4R | Fires: S1R |
|---|---|---|---|---|---|---|
| classic proteomics | 5 | 0.164 (5%) | 0.231 / 0.493 | 0.175 / 0.367 | 52.9% | 64.0% |
| classic proteomics | 7 | 0.058 (5%) | 0.231 / 0.493 | 0.215 / 0.466 | 86.6% | 90.6% |
| classic proteomics | 8 | 0.070 (5%) | 0.231 / 0.493 | 0.216 / 0.469 | 81.9% | 89.4% |
| transcriptomics | 7 | 0.039 (10%) | 0.221 / 0.500 | 0.213 / 0.478 | 89.1% | 93.2% |
| lipid/metab | 7 | 0.089 (10%) | 0.226 / 0.503 | 0.180 / 0.431 | 78.4% | 75.9% |
| lipid/metab | 8 | 0.195 (5%) | 0.226 / 0.503 | 0.209 / 0.469 | 53.9% | 54.7% |

S1R's true association: median 0.36 (proteomics), 0.39 (transcriptomics), 0.34 (lipid/metab); P90 0.72 to 0.76.
Calibration, for comparison: the true association was exactly 0 and the same limits fired in 4.3% to 6.0% (5% point)
or 9.1% to 9.5% (10% point) of unrelated runs (M24, M27).

## Result 2: false discoveries (C0R, the centre with true effects)

Mean false discovery proportion at the 5% target; type I at 0.05 in brackets. "Fired" means the runs where the design's
check exceeded the limit.

| Omic | Design | Design: all runs | Design: fired | Arm 1 (no correction): all / fired | Arm 2 (ceiling): all / fired |
|---|---|---|---|---|---|
| classic proteomics | 5 | 0.090 (0.069) | 0.119 (0.084) | 0.077 / 0.135 | 0.051 / 0.050 |
| classic proteomics | 7 | 0.056 (0.053) | 0.056 (0.054) | 0.077 / 0.088 | 0.051 / 0.050 |
| classic proteomics | 8 | 0.052 (0.051) | 0.052 (0.052) | 0.077 / 0.091 | 0.051 / 0.050 |
| transcriptomics | 7 | 0.049 (0.056) | 0.049 (0.056) | 0.054 / 0.059 | 0.048 / 0.048 |
| lipid/metab | 7 | 0.069 (0.052) | 0.072 (0.053) | 0.033 / 0.037 | 0.054 / 0.055 |
| lipid/metab | 8 | 0.061 (0.046) | 0.069 (0.049) | 0.033 / 0.047 | 0.054 / 0.056 |

By thirds of the true association (C0R+C4R together; `random-allocation-errors.tsv`), the designs other than
proteomics arm 5 hold type I at 0.043 to 0.057 in every third. Without correction, type I rises from about 0.012 to
0.016 in the lowest third to 0.12 (proteomics), 0.18 (transcriptomics) and 0.05 (lipid/metab) in the highest. Proteomics
arm 5 does not hold: type I is 0.097 and mean FDP 0.176 in the highest third (C0R+C4R), and 0.142 and 0.265 at S1R.
At 3 per group (S1R) lipid/metab arms 7 and 8 also exceed the target in the highest third (type I 0.13 and 0.12, mean
FDP 0.35 for both). That third's ceiling arm has type I 0.056 but mean FDP 0.21.

In C4R (no true effects) any call is false, so the "FDP" is the share of runs with any call. It reads high for every
arm in lipid/metab (arm 2 0.134, arms 7 and 8 0.176) and is not a 5% target there. That is why Result 2 uses C0R.

## What this means for "fires in about 5% of unrelated data sets"

That statement holds only for data sets where the factor is exactly balanced across the groups. With the factor
allocated at random and 6 per group, the true association has a median of about 0.23 and a 90th percentile of about
0.50. The validated limits sit well below that (0.04 to 0.19), so the check fires in 53% to 89% of such unrelated data
sets (S1R, 3 per group: 55% to 93%). The check therefore detects chance imbalance; it does not signal a design that is
"related" by construction. For the recommended designs (proteomics 7 and 8, transcriptomics 7) false discoveries stay at
target whether the check fires or not, so a firing check there is not a sign of excess false discoveries. For
proteomics arm 5 and lipid/metab 7 and 8, firing runs carry somewhat higher false discoveries (C0R mean FDP 0.07 to 0.12
against 0.05 to 0.06 for the ceiling). This supports M36's detection-limit wording with measured consequences. The
wording itself is the user's.

## Not as declared

- The run used the study mode with K124_ARMS plus a side table, not the association mode, because the association mode
  runs no test. This is stated in the declaration.
- Lifecycle identity: S1's C0 rows matched the earlier run in every count and rate but not in the last digit of
  bias and rmse (read as CPU-dependent OpenBLAS kernels; the code path is textually unchanged). S3 and S4 were exact.
- No other deviation. Nothing was rerun or adjusted after the set was read.
