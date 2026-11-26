# KIT-275 proving send: declaration, committed alone before the request

**Status.** Sent once after the owner's decision of the ceiling (VAL-382: 406.90 core-hours; the 1.5 x timeout rule with the continuation rule and the controls-only interim read). The scripts are the ones audited at d6b95a05: `git diff d6b95a05 HEAD` over `sim-k275.R`, `score-k275.R`, `fixture-k275.R`, `e0-diff.py` and `selftest/` is empty (verified at the declaration commit). SHA-256 (first 16 hex): `sim-k275.R` 2cc87ea569eef7fd, `score-k275.R` a43a2c10c9a41ec2, `fixture-k275.R` 407949a9a6c0f90f. `<commit>` in the request is this declaration's own commit; `K275_G1` = `docs/qa/workflow-execution/kit/kit-275/g1/g1.tsv`.
**Binding.** The package commit named in the request (the declaration's own commit, on staging ec27f69e0 as L1).

**One request, scripts in order (each its own explicit `timeout_minutes`):**
1. `selftest/run-selftest.R` (10): all cases (46), 0 differ; E0 0 unmarked.
2. **A** `sim-k275.R S,T D3 9201 9202 <p>/a` with `K275_JOB=k275-prove K275_WORKERS=4 OMP_NUM_THREADS=1`, counted mode, **`timeout_minutes` 2**.
3. **B** the same command into `<p>/b` with `K275_ONLY_UNITS=T-D3-9201,T-D3-9202` (8).
4. **Scorer 1** `score-k275.R <p> <s1>` with `K275_G1=<g1.tsv> K275_SCORE_R=200,9202,20 K275_BOOT_B=200` (5), uploading only `decision.tsv` and `settings.txt` (`out_pattern`).
5. **C** `<p>/c` with `K275_ONLY_UNITS=T-D3-9201` (8).
6. **Scorer 2** over `<p>` (5), with the same environment as scorer 1 (`K275_G1=<g1.tsv> K275_SCORE_R=200,9202,20 K275_BOOT_B=200`): without it RANGE-1 would fire before the double-source check.
The scorers' output directories lie outside `<p>`. Scorer 1 runs before C writes (the runner runs scripts in order).

**Predictions (declared before the send).**
- **A:** verdict `timeout` (exit 124, or 137 after the kill-after); files: `checks-k275-prove-start.tsv` with E3 for lipid TRUE; provenance `final: false`; four part files, none written as gz; no `checks-k275-prove.tsv`; no timing file; each T part file holds 1 to 5 data lines and each S part file 1 to 4 (4 = complete). Basis: L1's T unit at D3 reaches its rows at 20.9, 48.5, 73.2, 102.9, 151.0 and 181.6 s, so T cannot finish inside 120 s unless the runner is about 1.5 times faster than L1's load; S (132.9 s) may finish on a fast runner.
- **B:** exit 0; console `k275-prove-cont: 2 units`; `rows/T-D3-9201-9202.tsv.gz` with 12 rows; no part files; `timing-units-k275-prove-cont.tsv` with exactly the two tags and no error; provenance `final: true`.
- **Scorer 1:** exit 0; `interrupted data sets` = `T|D3|9201|a; T|D3|9202|a`, plus any S group of `a` with 1 to 3 rows; outcome `Stop (F1)` (no tier D rows: PC1 has no data: predicted, not a finding); `settings.txt`: rows = 12 + 4 x (the S data sets complete in `a`), from N + 1 files, N = the number of part files in `a` (the units with at least one row: 3 or 4; on a slow runner an S unit may write no row before the 120 s kill and then has no part file).
- **Declared alternative branch** (a runner more than about 1.5 times faster, so that any T unit finishes in A): scorer 1 stops "complete in two sources"; the interrupted path on real rows is then not shown; the send is reported as such and not re-run.
- **C:** exit 0; 6 rows. **Scorer 2:** exit non-zero; stderr contains "two sources".

**Ceiling.** Set-up 10 + self-test 10 + A (2 + 1 kill-after) 3 + B 8 + scorer 1 5 + C 8 + scorer 2 5 = **49 job-minutes = 3.27 core-hours** (36 = 2.40 without C and scorer 2; expected about 25 = 1.67). Running ceiling: used 12.43 + 3.27 + planned 391.20 = 406.90 (the owner's figure under the 1.5 x ruling); written in `runs.tsv` before sending: **12.43 used + 3.27 this request + 391.20 remaining = 406.90, equal to the owner's ceiling (VAL-382).**
**What it cannot inform.** R, seeds, margins, levels, the outcome table, the jobs and their timeouts are fixed before it; the ceiling is decided on L1's figures; it exposes S and T rows of two non-counted data sets, used for line-count checks only.

**Also predicted (J03).** The self-test step reports 46 cases and 3,128 predictions, 0 differences; E0 reports 41 hunks, 0 unmarked. The scripts table reads: script 1 ok, script 2 (A) timeout, scripts 3 to 5 ok, script 6 (scorer 2) failed; so **the run's overall outcome is "failure" by prediction** (a failed script). A, B and C are uploaded in full; scorer 2 writes no file.
