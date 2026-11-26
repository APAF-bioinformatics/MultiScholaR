# KIT-262 L2 (the M2 lifecycle): what it showed

Blind (VAL-11): the lifecycle's row files stayed under `tempdir()`; only checks, timing, hashes and console counts left it. `l2-plan.md` states the purpose, ceiling and pass rule before the first send.

| Send | Run | Commit | Result |
|---|---|---|---|
| 1 | 45714964374 | d6f9c081 | scripts 1 (E5b), 2 (scorer self-test, E6 on constructed inputs) and 5 (failure check) passed; scripts 3 and 4 hit their timeouts after 80 of 198 cells and 3 of 10 cells. Instrument defect (no result produced, no expected value touched): the loop ran one cell at a time, so a cell with one replicate used one of four workers. Fix: cells run in batches of about 40 data sets. Recorded in `runs.tsv` (instrument rerun 1 of at most 3). Files: `l2-run1/`. |
| 2 | 45721854918 | f92521d9 | scripts 3 (M2 at r = 91 in all 198 cells; E8 at r = 91), 4 (E8 at r = 92 to 95) and 5 (failure check) all passed. Files: `l2/`. |

Checks of send 2 (`l2/s1/checks.tsv`, `l2/s2/checks.tsv`, `l2/s3/checks.tsv`): E1, E1', E2, E3, E4, E7 pass in all 198 data sets of the main script (units 198, failed 0); E8: 10 of 10 (r = 91) and 40 of 40 (r = 92 to 95), metabolomics trace identical to lipidomics; COMPLETION 198 of 198 (and 40 of 40); E6 (the scorer twice on the lifecycle's own rows, byte-identical): 6 of 6 files. Failure check: with `K262_FAULT=1` the first task signalled after its data were generated; exactly one error row was written (`FAULT_ROW` 1 of 1), the other data set completed, and the scorer's self-test shows an errored data set flags its cell and is never a denominator. Send 1: E5b 12 of 12 (A to E, the cause split, the reason parser), scorer self-test 12 of 12 (Wilson against `e5hand-k262.py`; gain shown at 7 of 40 and not at 6; skip applied and not applied; C2-only gain; error row; report-only cells; stress-cell tie break; E6). E0 and `e5hand-k262.py` ran locally (7 of 7 and exit 0).

Instrument observation: the instrument was shown to observe completion (every planned data set has a row) and failure (an injected error becomes a recorded error row and a flagged cell), not merely startup.

Timing (`l2/s1/timing.tsv`, seconds per data set with four workers running side by side; `seconds_sweep` is the hooked sweep, `seconds_total` includes the unhooked E1' sweep): F = 300: 6.0, 6.6, 6.3, 6.5, 7.1, 9.1 s at n = 4, 6, 12, 24, 48, 96; F = 1,000: 15.7, 16.3, 16.7, 17.9, 19.4, 22.2 s. The ticket's basis was 6.0 and 20.0 s at n = 12; the n-factors are smaller than assumed (1.5 and 1.3 at n = 96 against 2.0).

Cost: send 1 44.8 job-minutes (2.99 core-hours); send 2 38.55 job-minutes (2.57 core-hours). Together 5.56 core-hours against the L2 ceiling of 3.00 and the instrument reserve of 10.00: send 1 is the ticket's L2 row, send 2 is charged to the instrument reserve.

## G2 (the counted jobs' plan, recomputed from the lifecycle's seconds)

Per counted data set: its hooked sweep plus preparation (seconds_total minus twice seconds_sweep); plus two E1' replicates per cell (r = 7, 27) and the E8 replicates in the ten E8 cells (r <= 5). Total 120,479 CPU-seconds (the ticket's expected 131,976, ceiling 196,560). Four jobs of about 126 minutes each at four workers (`m2-plan.tsv`: cells per job), each below the 150-minute limit; `timeout_minutes` 190 (1.5 x plan). Expected 4 x (126 + about 10 set-up) = 544 job-minutes = 36.3 core-hours (ticket 39.33 expected, 58.60 ceiling): the total stays inside the ceiling; no split is needed and R is not reduced.
