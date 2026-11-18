# KIT-129 R0, part B: the follow-up for the four near-miss designs, read once (2026-11-18)

Declaration: `followup-20261118.md` (decision M24; copied from `sim-k129-r0b`, where it was pushed first). The
blind check's result (`results-20261118.md`) stands as read. Nothing was re-set, rerun or adjusted after the
validation set was read.

## Order (pushes of `sim-k129-r0b`, GitHub's activity log, UTC: `followup-pushes.tsv`)

| Push | Commit | Pushed |
|---|---|---|
| Declaration and `threshold-followup.R` | `c31adab2` | 04:33:35 |
| The four thresholds and their setting-set shares | `ac90d782` | 04:33:52 |
| Launch params (starts run 44238222543, created 04:34:08) | `5a1ee379` | 04:34:06 |

## The run

Hosted-runner run 44238222543 (`r0b-followup`): S1 and S4, runs 2501 to 3000 of C0, C3, C4, C5, S3, S11,
association mode; 12 jobs, every job on its first attempt, no rerun. Every cell and run present once, no duplicate;
the four designs' rows: ok 13,240, no_factor 1,760, none crashed or refused (`followup-status.tsv`; the other arms'
rows of this range were not tabulated or scored). Raw tables not committed (`followup-assoc.sha256`). Core-hours
(job duration x 4, `jobs.tsv`): 38.2 (S1 26.8, S4 11.3); part B in all **253.7** (215.5 before).

## Result: accepted for all four designs

Share fired on the validation set, percent. Marks: correlation 0 at most 10.0; 0.3 or more at least 80.0.

| Omic | Arm | Threshold | Corr. 0 | C0 | C3 | C4 | 0.3 or more | C5 (0.3) | S11 (0.3, 3 per group) | S3 (0.6) | Result |
|---|---|---|---|---|---|---|---|---|---|---|---|
| S1 | 5 | 0.1638 | 4.7 | 3.2 | 10.8 | 0.0 | 98.9 | 100.0 | 96.8 | 100.0 | accepted |
| S1 | 7 | 0.05808 | 5.1 | 7.2 | 7.8 | 0.4 | 100.0 | 100.0 | 100.0 | 100.0 | accepted |
| S1 | 8 | 0.07007 | 4.3 | 2.0 | 7.6 | 3.2 | 100.0 | 100.0 | 100.0 | 100.0 | accepted |
| S4 | 8 | 0.1949 | 6.0 | 5.4 | 6.9 | 5.7 | 91.5 | 97.3 | 78.6 | 98.6 | accepted |

Setting-set shares of the same thresholds (not results; `setting-followup-S<n>.tsv`): 5.0% at correlation 0 for each;
at 0.3 or more 98.7% (S1 arm 5), 100.0% (S1 arm 7), 100.0% (S1 arm 8), 92.7% (S4 arm 8; S11 alone 79.8%).

Single correlation-0 cells are not uniform under these thresholds (see the C0, C3, C4 columns): the 10% mark is
declared and met on the three cells together.

## The second-run refusal

After the results were written, `threshold-followup.R check` was invoked a second time for S1 and for S4: both
stopped ("result-followup-S<n>.tsv exists: the validation set is read once (one attempt)", exit 1;
`followup-second-check.txt`), and the sha256 of all eight follow-up files was the same before and after.
