# KIT-184 L1: the simulation's lifecycle run

Run 45333818847 (evidence lane, branch `run-mofa-evidence-k184` at 96787082; package commit 21053957; one job of 21.5
minutes and the scorer, 0.9 minutes). Scripts (`l1/script-sha256.txt`): `simulate-k184.R` 4ea61d24..., `score-k184.R`
5487ab06..., KIT-183's `simulate-k183.R` 4cdf7105... and `score-k183.R` 7efd7c68... (unchanged). The lane's first push
(45333717225) was refused by GitHub as an invalid workflow file (a colon in the unquoted `name:`); no job ran.

**L1 ends at a stop rule of ticket section 11: a replicate takes 111 to 120 seconds (median per cell) on one core,
above the 70 seconds at which the ticket stops and returns the numbers to the owner with a re-split. Nothing has been
declared or counted.** Every other check of L1 passed. The rates below come from 5 replicates per cell: they show that
every arm, class and factor runs to completion, and they decide nothing.

## Observed

| Check (ticket 11, L1 row) | Result |
|---|---|
| 5 replicates of each of cells A, B and D to completion with every arm, class and factor | 15 of 15 completed (`l1/completion.tsv`); every arm (1, 2, 3, 3o, D0, D1), every test (`camera`, `camera_pr`, `camera_floor`, `gsea`, `permutation` where declared) and every factor (F_main, F_ind, F_ret; F_main and F_unw of the default model) produced rows (`l1/error-rates.tsv`) |
| Seconds per replicate on one core (stop above 70) | **cell A 111 (101 to 117), cell B 113 (82 to 118), cell D 120 (72 to 133): above the stop.** Of this, the two trainings take a median of 50 s (A, B) and 58 s (D) |
| Parity: KIT-183's 12-sample cell, replicates 1 to 3, through `simulate-k184.R` (blocks off) and through KIT-183's own file | identical: 10 tables of 10, the status table but for its seconds; 3 of 3 completed on each side (`l1/parity.txt`) |
| A declared invalid argument (7 biological samples) | observed and recorded as `did_not_complete` with its message (`l1/ruvfail-status.tsv`) |
| The lane's ruv and RUVIIIC (gate G3) | ruv 0.9.7.1, RUVIIIC 1.0.19, limma 3.68.4, MOFA2 1.22.0: the locked versions (`l1/versions.tsv`); no extra install step was needed |
| The written factors (`exact` expected) | `exact` in 45 of 45 data sets; largest s_(k+1)/s_1 3.2e-15; largest relative residual 2.9e-12 (`l1/exactness.tsv`) |
| The identity check (arm 3's `correlation` against the same model fitted to the corrected values; below 1e-8 expected) | largest difference 1.8e-14 over every data set, factor and replicate (135 comparisons, none above 1e-8); for proteomics and metabolomics, whose levels are restored, the p-values agree within 1.4e-13; for transcriptomics the p-values differ (its level is not restored, so the trend covariate moves), as ticket 5.2 says (`l1/identity.tsv`) |
| What the removed factors caught: the smallest canonical correlation between the written and the planted unwanted factors, median per data set (stop to the owner below 0.9) | transcriptomics 0.994, 0.987, 0.994 (cells A, B, D); proteomics 0.950, 0.965, 0.976; metabolomics 0.983, 0.974, 0.916: **no median below 0.9**. Single replicates were lower: proteomics 0.48 (cell A), transcriptomics 0.70 (cell B), metabolomics 0.69 (cell A) (`l1/caught.tsv`) |
| p_floor at least p_camera, set by set, under every arm that runs both | largest p_camera minus p_floor: 0 (`l1/floor.tsv`) |
| Pools: classes forming fewer than 90 % of their declared sets | proteomics C5 at size 100: 67 % (its pool is about 100 features); transcriptomics C4 (a reported class): 56 to 65 % at sizes 10, 30, 100. Every other class and size formed all its sets (`l1/pools.tsv`). Both are to be redefined in a dated section of the declaration before the counted run (ticket 6.2) |
| Feature-level rates | computed for every arm on F_ind, F_ret and F_main, except metabolomics on F_main: the model takes every metabolomics feature (250, fewer than its 1,000 per data set), so "the clean complete features the model did not select" is empty there. To be stated in the declaration |

C5's rate (`camera`, share of C5 sets with p at most 0.05) on these 5 replicates, as the L1 row asks (65 to 75 sets
per cell for transcriptomics and proteomics, 25 for metabolomics; **5 replicates: not a result**):

| Cell, factor | Data set | Arm 1 | Arm 2 | Arm 3 | Arm 3o |
|---|---|---|---|---|---|
| A, F_main | transcriptomics / proteomics / metabolomics | 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / 0.000 | 0.013 / 0.000 / 0.000 |
| A, F_ind | transcriptomics / proteomics / metabolomics | 0.013 / 0.040 / 0.040 | 0.173 / 0.000 / 0.000 | 0.040 / 0.093 / 0.040 | 0.067 / 0.173 / 0.080 |
| B, F_main | transcriptomics / proteomics / metabolomics | 0.240 / 0.167 / 0.120 | 0.000 / 0.127 / 0.280 | 0.200 / 0.107 / 0.080 | 0.013 / 0.013 / 0.000 |
| B, F_ind | transcriptomics / proteomics / metabolomics | 0.000 / 0.000 / 0.040 | 0.000 / 0.080 / 0.000 | 0.027 / 0.120 / 0.040 | 0.067 / 0.260 / 0.040 |
| D, F_main | transcriptomics / proteomics / metabolomics | 0.053 / 0.020 / 0.040 | 0.000 / 0.000 / 0.000 | 0.040 / 0.020 / 0.000 | 0.000 / 0.020 / 0.000 |
| D, F_ind | transcriptomics / proteomics / metabolomics | 0.027 / 0.113 / 0.080 | 0.000 / 0.000 / 0.000 | 0.027 / 0.147 / 0.200 | 0.040 / 0.207 / 0.080 |

Every figure is in `l1/error-rates.tsv`; this table is a reading aid.

## Inferred, not shown

- Why a replicate takes more than twice the 45 s the ticket estimated: two trainings (50 to 58 s against the 2 x 18 s
  inferred from KIT-183's replicates), and about 60 s of tests (fifteen analyses of a data set and factor per
  replicate, each with up to thirteen test calls, among them the pre-ranked test and the permutation). Not profiled.
- On cell B's main factor the adjusted test (arm 3) called C5 sets in 8 to 20 % on these 5 replicates and the oracle arm
  (the planted factors in the design) in 0 to 1 %. With 5 replicates this is an indication only; it is the declared
  check K-E4, which the counted run measures and which does not hold the merge (outcome R7).

## What the stop needs: a re-split, for the owner

The replicate count is not reduced (ticket 11). At 120 s per replicate (cell D's median; 133 s was the largest), four
single-thread workers per job:

| | Ticket (45 s) | Re-split at 120 s |
|---|---|---|
| Jobs for cells A (800), B (800), D (400) | 8 of 250 replicates | 16 of 125 replicates |
| Minutes per job (replicates x seconds / 4 / 60 + set-up) | 46.9 + 10 = 56.9 | 62.5 + 11 = 73.5 (set-up measured: 10.6 min); timeout 85 as now, which leaves 133 s per replicate |
| Determinism job (10 replicates) | 11.9 | 5 + 11 = 16; timeout 25 |
| S expected | 466.9 job-minutes = 31.13 core-hours | 16 x 73.5 + 16 = 1,192 job-minutes = 79.47 core-hours |
| S ceiling | 705 job-minutes = 47.0 core-hours | 16 x 85 + 25 = 1,385 job-minutes = 92.33 core-hours |
| S' (one instrumentation correction) ceiling | 65.67 core-hours | 92.33 core-hours |
| Ticket expected | 140.35 core-hours | 140.35 - 31.13 + 79.47 = 188.69 core-hours |
| Ticket ceiling | 303.33 core-hours (accepted, VAL-58 (2)) | 303.33 - 47.0 - 65.67 + 92.33 + 92.33 = 375.32 core-hours |

So the expected cost stays inside the accepted ceiling of 303.33 core-hours, and the ceiling itself would rise to
375.32 core-hours if an instrumentation correction of the whole simulation were ever needed. The owner decides: (a)
the re-split as above with the ceiling raised to 375.32; (b) the re-split with the ceiling kept at 303.33, which
leaves room for S and F but for a correction S' of at most 20.3 core-hours (not a whole repeat); (c) another split.
Nothing in the generator, the correction, the arms or the classes was changed to make it faster, and nothing is
proposed to be cut.

Also due before the declaration (ticket section 11 and 20.1): the coordinator sees these seconds with KIT-191's added
arms switched on (they were: `camera_floor` under arms 3, 2, 3o and D1, the permutation under arms 2 and D0), and the
auditor re-checks the budget's re-derivation.
