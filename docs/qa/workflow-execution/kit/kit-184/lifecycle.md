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
| Pools: classes forming fewer than 90 % of their declared sets | Counted per replicate (a size forms in a replicate when its pool holds at least that many features; `l1/pools.tsv`, and the workers' `formed.tsv`): proteomics C5 at size 100 formed in 10 of 15 replicates (4, 3 and 3 of 5 in cells A, B, D; the same under both models); transcriptomics C4 (a reported class) under the default model formed in 2 of 15 replicates at size 10 and in none at sizes 30 and 100, and under the model trained on corrected values in 15 of 15 at sizes 10 and 30 and in 14 of 15 at size 100 (4 of 5 in cell B); metabolomics C4 at size 20 formed in 14 of 15 replicates (4 of 5 in cell B; the same under both models). Every other class and size formed all its sets. To be redefined in a dated section of the declaration before the counted run (ticket 6.2): proposals in `declaration-draft.md`. *Corrected 2026-11-24: this row first gave transcriptomics C4 as "56 to 65 %" (a share pooled over models, factors and cells, which hid that the failures are almost all under the default model) and said every other class formed all its sets, which missed metabolomics C4 at size 20 in cell B.* |
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

## After the stop (2026-11-24)

The coordinator ruled (DISPATCH-STATE, engineering): the re-split to 16 jobs of 125 replicates is accepted (S expected
79.47 core-hours, ceiling 92.33); the ticket's ceiling stays 303.33 core-hours, so the one correction S' is limited to
re-running affected jobs within 20.33 core-hours (305 job-minutes), and a whole repeat would be a stop to the owner.
Ticket section 11, revision 1.5, has the re-derivation row by row. Two figures above are corrected there:

- The ticket's expected total after the re-split is 192.02 core-hours, not 188.69: the figure above started from
  140.35 and left out the run of the tests-only commit (R1a, 50 job-minutes), which revision 1.4 had already added
  (143.68).
- The medians above mix replicates run under four busy workers with replicates run under fewer. The lane deals five
  replicates into the ranges 1, 2, 3 and 4 to 5, so the last replicates of a cell ran while fewer workers were busy.
  Under four busy workers: cell A 110 to 117 s (replicates 1 to 4), cell B 110 to 118 s (1 to 4), cell D 120 to 133 s
  (1 to 3). Under fewer: A 101 s, B 82 s, D 72 and 77 s. The budget's 120 s per replicate stays a fair expectation for
  the whole (800 x 112.5 + 800 x 113.9 + 400 x 127.8 s = 64.5 core-hours of replicates against the 66.7 budgeted), but a
  job holding only cell D replicates would run its last replicate close to the job's limit (31 x 127.8 s = 66.0 minutes
  against 68.4 available before the last-start time); the job list in `declaration-draft.md` therefore gives every job
  100 replicates of cell A or B and 25 of cell D.

## L1b, the blind lifecycle extension: written before it is sent (2026-11-24)

Authority: declaration audit K184-D05 and its re-check 1 (K184-R02, "order before the freeze", steps 1 to 5); the
coordinator's word for L1b, given on the condition that this section is written first. Not a VAL-54 re-run.

**Purpose.** To run, before anything is frozen, the code L1 did not run: the simulation writing its error table per set
size and its package versions in every replicate; the scorer summing the sizes, reading proteomics C5 in two
compositions, and comparing a determinism pass; and to measure the scorer at the counted run's size (about 2,000
replicates). It shows mechanics and counts. It decides nothing about the simulation's results, and no rate leaves
the job.

**Commits.** Package `kit-184-ruv-factors-in-integration` at 9e525fd94ec585cc116919884cea83cc1d014091
(`simulate-k184.R` SHA-256 58333b1ef1b479a65f30089dd0307aad5ebab2495c1fb060621a8738a873a75b, `score-k184.R`
282b844e1c3427e8c99f296aa486c313957776bda8067d782a9a4511b43fa8b2; no file under `R/` differs from L1's commit
21053957). Lane `run-mofa-evidence-k184` at 41013a465eb477bf5a27f65c0685b41e996a9c68 (`k184-evidence.yml` SHA-256
2b680757025310996e0d35e70a11c3d0b94034f8173f9a03717bbc754826b6fe; `k184.R` the byte copy of `simulate-k184.R`,
58333b1e...). Sent by a commit to the lane's `k184-params.json` with `run = blind` and that package commit.

**Ceiling.** One job with a timeout of 50 minutes, and the parameters job (seconds): inside what is left of L1's row,
52.3 job-minutes (75 less the 22.7 L1 used). Expected about 25 job-minutes: set-up 11, the replicates 5, the
determinism pass 2, the scorer under 1, the scale test a few minutes (its own limit is 25). One attempt. A failure in
the mechanics is fixed and sent again under VAL-54 (at most three sends in all) only if what is left of the 52.3
job-minutes covers another job; if it does not, I stop and report.

**What it decides.** Whether the two scripts at 9e525fd9 may be the frozen ones (every item of the list below), and
the `score` job's timeout and budget line (three times the scale test's measured time, at least 30 and at most 90
minutes). A failed item is a defect of the instrument, not a result.

**The check that nothing but the listed files can leave the job** (the workflow at 41013a46, read again line by line
after the commit):

| Way out of the job | What the workflow does for `run = blind` |
|---|---|
| Artifacts | Two upload steps. The one for the whole result folder runs only when `run != 'blind'` (line 369). The other runs only when `run == 'blind'` and uploads `RUNNER_TEMP/staged` (line 377), a folder the L1b step creates and fills by name: `plan.txt`, `runner.txt`, `script-sha256.txt`, `status.tsv` and `formed.tsv` (the workers' rows), `completion.tsv`, `not-completed.tsv`, `pools.tsv`, `versions.tsv`, `determinism.tsv`, `determinism-selftest.tsv`, `l1-reference-check.tsv` (copied one by one from the scorer's folder if they exist), `scale-test.txt`, the workers' logs (`w*.log`, `det-w*.log`) and `tables.tsv`. The workers' tables, the scorer's folder and the scale test's folders are in other folders of `RUNNER_TEMP` that no step uploads; the scale test's output folder is deleted |
| The `score` job | Skipped: its condition is `run != 'blind'` (line 386). It is the only step that writes to the step summary (line 419) and the only other upload (line 420) |
| The step summary | Nothing in the evidence job writes to it |
| Job outputs | The evidence job has none; the parameters job's are the commit, the run kind and the job list |
| The log | The simulate step prints the set-up time, the workers' ranges, the seconds and a count of replicates per status. The L1b step prints the scorer's exit status, the scale test's line and the names of the staged files. The scorer runs with `K184_BLIND`, so it prints no summary, and its standard output and error go to a file that is not staged |

Not on the re-check's list, and in the job on purpose (declaration draft, section 6): `tables.tsv` has a column with the
rows per cell; if the scorer exits with an error, the lines of its standard error from `Error` on (at most twelve)
are printed, which is R's error message and call. Also possible: `withheld.txt`, the lane's standing token rule, which
names a staged file it removed for naming an internal package. The workers' logs hold R's warnings only: the
simulation prints nothing itself (`simulate-k184.R` has no `cat`, `message` or `print`).

**How it is judged**: from the staged files alone, against the re-check's list. Nothing else of the run is opened,
downloaded or printed.

| Item of the list | Read from |
|---|---|
| every cell A replicate completed | `status.tsv`, `completion.tsv` |
| at least one cell D replicate `not_started`, and counted | `status.tsv`, `completion.tsv`, `not-completed.tsv` |
| `ke1-checks` and `ke5-checks` 34 rows, `ke4-checks` 4 rows, per deciding cell (cell B with missing values) | `tables.tsv` (rows and rows per cell) |
| `error-rates-by-size` exists | `tables.tsv` |
| `versions.tsv` has four rows | `versions.tsv` |
| `determinism.tsv` 12 of 12 TRUE; the altered copy FALSE | `determinism.tsv`, `determinism-selftest.tsv` |
| `l1-reference-check` TRUE | `l1-reference-check.tsv` (with its counts if FALSE) |
| the scale test ends within 25 minutes and below 12 GiB | `scale-test.txt` |
