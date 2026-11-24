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
| The `score` job | Skipped: its condition is `run != 'blind'` (line 388). It is the only step that writes to the step summary (line 419) and the only other upload (line 420) |
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

## L1b: what came back (mechanics and counts only)

Run 45348572248 (lane commit 5f75281c, which changed only `k184-params.json`; the workflow as at 41013a46; package
9e525fd9). The parameters job took 5 seconds and the one evidence job 18.6 minutes: 18.7 job-minutes, 1.25 core-hours;
L1's row has now used 41.4 of its 75 job-minutes. The `score` job was skipped. The run has one artifact,
`k184-lifeb-staged`; I downloaded that and nothing else, and did not open the job's log. Its files, without the
workers' logs, are in `l1b/`. First send; no failure; nothing re-sent.

| Item of the re-check's list | Value | From |
|---|---|---|
| every cell A replicate completed | 5 of 5 | `l1b/status.tsv`, `l1b/completion.tsv` |
| at least one cell D replicate `not_started`, and counted | replicate 5 (`job deadline`); cell D: 4 completed, 1 not started, 0 failed | `l1b/status.tsv`, `l1b/completion.tsv`, `l1b/not-completed.tsv` |
| `ke1-checks`, `ke5-checks`: 34 rows per deciding cell | 68 rows each: A 34, B 34 | `l1b/tables.tsv` |
| `ke4-checks`: 4 rows per deciding cell | 8 rows: A 4, B 4 | `l1b/tables.tsv` |
| `error-rates-by-size` exists | 10,995 rows (A 5,556; D 5,439), with a `nominal` and a `replicates` column | `l1b/tables.tsv` |
| `versions.tsv` has four rows | ruv 0.9.7.1, RUVIIIC 1.0.19, limma 3.68.4, MOFA2 1.22.0; every worker wrote four rows per replicate | `l1b/versions.tsv`, `l1b/tables.tsv` |
| `determinism.tsv` 12 of 12 TRUE | 12 of 12 (errors 10,652 rows on each side, power 1,272, formed 760, counts 30, feature 72, identity 18, floor 30, exact 6, caught 6, factors 10, orientation 2, status 2) | `l1b/determinism.tsv` |
| the altered copy FALSE | FALSE | `l1b/determinism-selftest.tsv` |
| `l1-reference-check` TRUE | TRUE: 2,237 rows here, 2,237 in L1, 2,237 matched; sets equal, replicates equal, 0 rows with a different T | `l1b/l1-reference-check.tsv` |
| the scale test within 25 minutes and below 12 GiB | 250 copies, 2,250 completed replicates, 11,933,250 error-table lines: 119 seconds, peak 5,347,588 KiB (5.1 GiB), exit status 0 | `l1b/scale-test.txt` |

Also in the staged files: the scripts' SHA-256 as the job read them (`simulate-k184.R` 58333b1e..., `score-k184.R`
282b844e..., the workflow 2b680757...: those of the pre-send record); set-up 616 seconds; the replicates took 84 to 99
seconds in cell A and 100 to 110 in cell D under four workers (73 for the one that ran with fewer); the runner was
another processor than L1's (Intel Xeon Platinum 8573C; L1 ran on an AMD EPYC 7763), and cell A's rates still equal
L1's in every row. Because this runner was faster, cell D replicate 4 started before the deadline and completed; the
acceptance asks for at least one replicate not started, which replicate 5 is.

One word for the auditor: `operating-characteristic.tsv` has one row per deciding cell (2 rows), whose `checks`
column counts the rows that entered; the declaration's sentence says the file "counts the proteomics C5 check twice
(34 rows per cell)", which is true of what it counts and not of how many rows the file has.

**The scorer's limit and budget line (re-check K184-R03), set from the scale test.** Three times 119 seconds is 6
minutes, below the floor, so the `score` job's timeout is 30 minutes (lane commit ef6ce13d; `k184-evidence.yml` SHA-256
ae0dc5a52c2735d5b6092fca60455c4f1e9093a7de08b3a7b86c98f04e03eea7). Expected about 3 job-minutes (2 of scoring; L1's
score job took 0.9 minutes in all with its set-up), ceiling 30, taken from S', which becomes 305 - 30 = 275
job-minutes (18.33 core-hours): at most three of S's jobs at their 85-minute timeout, or two and the determinism job.
The scale test had two cells where S has three, and the same nine replicates copied, so the counted run's scorer may
take somewhat longer than 119 seconds; the limit leaves fifteen times that.

**Not checked by L1b:** any result (by design); the counted run's own job list (it has not run); the `score` job as a
separate job on downloaded artifacts (L1b ran the scorer inside its job).

## S, the counted run: the conditions for sending it, written before the send (2026-11-24)

Authority: re-check 2 of the declaration audit ("the first thing that must be true before S is sent") and the
coordinator's word, given for when all of these hold.

| Condition | Evidence (read just before the send) |
|---|---|
| The freeze commit is on the remote package branch | `git ls-remote origin kit-184-ruv-factors-in-integration` gave 81ee30527cbfab973d5cd77a9ed31207234643be, the freeze commit, whose `git show --stat` lists one path, `docs/qa/workflow-execution/kit/kit-184/declaration.md` (432 lines added). (The commit that adds this section comes after it; the freeze commit stays its ancestor and is the `package_commit` of the send.) |
| The freeze commit's scripts hash to the declared values | `git show 81ee3052:.../simulate-k184.R | sha256sum` = 58333b1ef1b479a65f30089dd0307aad5ebab2495c1fb060621a8738a873a75b; `score-k184.R` = 282b844e1c3427e8c99f296aa486c313957776bda8067d782a9a4511b43fa8b2; `git diff --stat 9e525fd9 81ee3052` over `R/`, `inst/`, `DESCRIPTION`, `NAMESPACE`, `tests/`, the two scripts and `kit-183/` is empty |
| The lane branch's head is exactly the lane commit the declaration names | `git ls-remote origin run-mofa-evidence-k184` gave 27c8766ca187b314b1951f31b30d43dadf0287fd; at that commit `k184-evidence.yml` hashes to 15378427e28bdbcd2066584167c74ae6e7b0b392d07041f6bf166ed0e3ce68d8 and `k184.R` to 58333b1e...a75b |
| The send commit changes only `k184-params.json`, with `run = counted` and `package_commit` = the freeze commit's full id | prepared and not yet pushed: lane commit 008852f82c41103a43a3b0e333febda81f7feee5, parent 27c8766c; `git show --stat`: `.github/mofa-smoke/k184-params.json`, 1 file changed, 3 insertions, 3 deletions; `run` is `counted`, `package_commit` is 81ee30527cbfab973d5cd77a9ed31207234643be |
| `gh run list` shows room for the 17 jobs | 19 runs of other tickets were not completed; their jobs: 77 in progress, 0 queued. So the account was running far more than 20 jobs at once and none was waiting for a runner; no full suite was queued. The 17 jobs are sent on that reading; if they queue, they wait (a job's limits count from its own start) |

Purpose: the counted simulation, once. Ceiling: 16 jobs of 85 minutes, the determinism job of 25, the score job of 30:
1,415 job-minutes. It decides K-E1 to K-E5 as `declaration.md` states them. While it runs nothing of it is read but the
status of its jobs, with one exception the coordinator asked for: the hashes in the run's `script-sha256.txt`, which
must be the declared ones (if they are not, the run is cancelled).

## S: what happened (2026-11-24)

Run 45353248796, sent by lane commit 008852f8 (only `k184-params.json`: `run = counted`, `package_commit` =
81ee30527cbfab973d5cd77a9ed31207234643be). All 17 evidence jobs started at once (none queued).

- **The hashes the run read.** A running job's log is not served (the request returned "BlobNotFound"), and its
  artifact exists only when the job has finished, so the run's hashes could not be read before any job finished. They
  were read from the log of the first job to finish, the determinism job (16 minutes), by extracting the five hash
  lines and nothing else: `simulate-k184.R` 58333b1e...a75b, `score-k184.R` 282b844e...a8b2, `k184-evidence.yml`
  15378427...68d8, KIT-183's two files 4cdf7105... and 7efd7c68...: the declared values. (Before the send they had
  been confirmed from git at the two commits; every job also stops if the lane's `k184.R` is not the package's
  `simulate-k184.R`.) No cancellation.
- **While it ran** nothing else of it was read but the status of its jobs.
- **Jobs.** All 17 evidence jobs and the score job succeeded at the first attempt: 42.9 to 77.5 minutes for the
  sixteen replicate jobs (limit 85; the last replicate may start until minute 79), 16.1 for the determinism job, 2.5
  for the score job, 0.1 for the parameters: 1,043.0 job-minutes, 69.53 core-hours (expected 1,197, ceiling 1,415). No
  replicate was left unstarted and none failed, so S' was not used and nothing was re-scored. The jobs ran on five
  processor models (from each job's log, the `Model name` line only): AMD EPYC 7763 (6 jobs, the determinism job among
  them), AMD EPYC 9V74 (4), Intel Xeon Platinum 8573C (4, job s01 among them), AMD EPYC 9V45 (2), Intel Xeon Platinum
  8370C (1).
- **After it finished** I downloaded the `score` artifact (the results; `s/`), and, to size the determinism
  difference the scorer reported, the artifacts of the determinism job and of job s01, of which I compared the twelve
  worker tables of cell A replicates 6 to 15 column by column and read `runner.txt` and the first lines of `plan.txt`.
- **Results**: `results.md`. **The texts**: `wording.md`. The full suite (A17) has not been run; the package's texts
  still hold the placeholders.
- **Where the per-replicate tables are** (results audit, for the coordinator): only in run 45353248796's artifacts
  (`k184-s01` to `k184-s16`, `k184-det`, `score`; about 86 MB compressed; kept 90 days from 2026-11-24). They are needed
  to recompute the criteria; the branch holds the scorer's tables (`s/`) and not the workers'. Whether a durable copy
  is kept is the coordinator's decision.
- **Determinism**: recorded in `results.md` as not met as declared, with its size and cause (the owner's decision,
  VAL-160 (3)); no tolerance adopted after the fact, nothing re-run, R3 and S' not invoked.

## After the owner's decision VAL-160: the regression runs, and F declared before it is sent (2026-11-24)

The owner decided that the adjustment ships off by default with a stated limit (VAL-160). The change and its texts
are commits a613b9cb (tests first), 74ce8305 (the change) and 5ba85c1f (one test call). Runs so far (`runs.tsv`):
R2 (45388084875: the restated expectations fail on the tests-only commit, 79 in 10 blocks, and pass on the change;
check identical; the four mutations caught; `r2/`), M2 (45398052678: K13 on real MOFA2), B2 (45398091058: the changed
file on the installed bundle, ubuntu-22.04 and windows-2025), and the same file on ubuntu-24.04 (45401782250,
running when this was written). The first R2 request (45387864719) was refused at `plan` for a requester line over 300
characters; no job ran.

**F, the one full suite base against head (acceptance criterion A17), declared before sending.**

- Purpose: the merge checkpoint. Every test file on the head against the same on the base, with R CMD check on both
  sides. It decides A17: no file may be worse on the head than on the base.
- Bindings: head = 0a51179569e2cf1d9c0774af038497a963b05f0d, the merge of this branch (5ba85c1f and the records after
  it) with `v0.9-staging` at 8cb536cc3bc09591ec1bfb6fef0cb15437a5ad4a (staging moved from a1b397b1 while the lanes
  ran; merged twice; the one conflict both times was the methods ledger, resolved as staging's rows followed by this
  ticket's three `VAL-<merge>` rows). Base = that staging commit. Relative to the base the head changes seven files
  under `R/`, the starter workbook, one manual page, NEWS, the ledger, ten test files, one reference table, and this
  ticket's records folder.
- Request: `tests = "all"`, 20 shards per side, `check: true`, the five long files at 120 minutes, as the recent full
  suites of other tickets were requested.
- Ceiling: row F of ticket section 11: expected 985.3 job-minutes (65.69 core-hours), ceiling 1,080 job-minutes (72
  core-hours). **The row's ceiling may be exceeded**: the coordinator's figure for recent full suites is 67 to 80
  core-hours (1,005 to 1,200 job-minutes), so up to 120 job-minutes (8 core-hours) above the row. The ticket's rule for
  that case stands (section 11): F beyond 1,080 job-minutes is not repeated or split; its times go to the owner. The
  ticket's ceiling is not at risk: 118 of 303.33 core-hours used before F.
- Attempts: one. A failure in the mechanics of a shard (runner lost, time) is a re-run of the failed shards only
  (VAL-54); a file worse on the head is a result and is reported, not re-run.
- The account before sending (`gh run list`): 13 runs not completed, 51 jobs in progress, none queued; two full
  suites of other tickets were finishing (7 jobs left between them) and none was queued.

## F and the last runs: what came back (2026-11-24)

- **F was sent once** (run 45402062487, request `k184-f-full-0a511795`, head 0a511795 against staging 8cb536cc) and ran
  in one attempt; every job succeeded. 1,268.4 job-minutes (84.56 core-hours): **188.4 above row F's ceiling of 1,080**
  (the declaration above allowed for up to 120). By ticket section 11 it is not repeated or split, and its times go to
  the owner. The ticket's ceiling is not reached (205.2 of 303.33 core-hours used).
- **Result** (`f/summary.md`, `f/compare.tsv`): 466 files; 455 `same`; the new file passes (493, 0 failed); this
  ticket's other test files pass with more expectations; check identical on both sides. **One file is worse on the
  head**: `test-gui-entry-proteomics.R`, lane `prot_dia_limpa`, four expectations (the steps `summary_save`,
  `summary_copy`, `summary_report`, `summary_export`). The stored reference
  `tests/testdata/gui-entry/prot_dia_limpa.json` lists the files of the proteomics integration input, and with RUV
  applied the head writes one more, `integration/proteomics/removed_proteomics.parquet` (`f/moved-references.txt`).
  That is this ticket's intended change (ticket 7.2). When the ticket was written no GUI-entry reference listed a
  Parquet file (ticket 8.5); KIT-199 re-froze this one on staging since (bf2b04d7). Ticket 8.5 foresees the case ("any
  reference its full suite finds moved (a capture that lists an integration input's files ...), with its size ... on
  merge"). The implementer's brief forbids re-recording a stored reference unless the ticket names the file, and the
  ticket names the kind, not the file: **not re-recorded; for the coordinator.** Its size: one artifact entry in each
  of four steps.
- The same file name appears in the failure messages of four captures that fail on the base already (`trx_main`,
  `trx_managed` in `test-omic-summary-module.R` and `test-omics-summary-core.R`): their references are stale on
  staging for other reasons; when they are re-frozen they will list `removed_transcriptomics.parquet`.
- **The same file on ubuntu-24.04** (45401782250): 493 passed, 0 failed. **Gate off at the final head** (45463933843):
  pages identical to staging's.
