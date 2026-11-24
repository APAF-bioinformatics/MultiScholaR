# KIT-197 lifecycle L1 (blind)

Run 45316713443 (request `k197-l1-5a702b7b`), package commit 5a702b7b6484e97e03ac515b2dd2c8b0b130c80e (the scripts as
declared), R 4.6.1, reference BLAS, runner AMD EPYC 7763 (4 cores). Scripts job 18.7 min (script 383 s); 1.25
core-hours against 1.33 expected and 2.67 ceiling. Run once; every check passed at the first attempt.

Blind: the lifecycle rows (cells A to D, r = 9001, 9002, six masks) and every score computed from them stayed in the
runner's temporary directory. Uploaded only: `l1/checks.tsv`, `l1/timing.tsv`, `l1/hashes.tsv`, the console
(`l1/console.txt`). No estimate was seen.

| Check | Pass | Observed |
|---|---|---|
| E0 | yes | local text check (`e0-k197.py`, `l1/e0.txt`): the six verbatim blocks equal `filter-effect.R` at 6416dbd6, one differing line (215, as stated) |
| E1 | yes | 48 of 48 (cell, mask, r): the in-script old rule equals `filterKept` on a worktree of 604ec4dd in a separate Rscript |
| E2 first half | yes | 48 of 48 DA calls: `logFC` equals the observed mean difference, max deviation 3.6e-14 |
| E2 second half | yes | 48 of 48: `id`, `p`, `q` of the modified `daTable` identical to S2's `daTable` (read from the committed file, hash checked) |
| E3 (F1) | yes | O ⊆ H in all 48 |
| Determinism mechanics | yes | `K197_DETERMINISM=1` on cell A r = 9001, 9002: 2 of 2 identical |
| E4 | yes | constructed table (`e4hand.py`): pooled and per-setting estimates of ρ_O, m6_NA, m6_ℓ, FDR_Δ, FDR_Δ^cal, E_cal, E_comp, both splits, D within relative 1e-12 of the hand values; jackknife SEs (pool and setting) likewise; labels bias, attribution "not specific to value dependence", outcome O1 |
| E5 | yes | the deciding scoring twice on the L1 rows: 11 files each, all byte-identical |
| G1 mechanics | yes | the descriptive mode on the L1 rows stops at G1 (exit 1, class `k197_g1`), 96 comparisons written, no estimand written |
| E6 | yes | killed (SIGKILL to the process tree) with 1 of 4 tasks complete; continued with the 3 missing; assembled rows and metrics (seconds excluded) identical to an uninterrupted run; each replicate once with six masks |

Departure in mechanics (no estimand or rule affected): E6 kills the run with processx's `kill_tree()` as soon as one
task is complete, instead of a `timeout` with a fixed delay, so that the kill falls mid-run whatever the load time.
The console's stray line "R 20131" is the printed return value of `kill_tree()` (a process id).

**G2.** L1: 5.71 s per replicate-mask (4 workers; 5.31 s in E6's uninterrupted run). Planned script wall per counted
job: 2,400 × 5.71 / 4 = 3,426 s = 57 min (+ load and scoring), below 150 min: no split; `timeout_minutes` 90 stands.
