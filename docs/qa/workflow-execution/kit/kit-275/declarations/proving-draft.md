# KIT-275 proving send: DRAFT declaration (the auditor's design, RH5). NOT SENT: it waits for the owner's figure (running-ceiling rule) and the auditor's closure of the changed scripts.

**Purpose.** Show the continuation paths on real counted-mode directories before any counted job: counted mode, a job killed by its timeout, the part files, `K275_ONLY_UNITS`, source tagging, the interrupted-data-set path and the double-source stop of the scorer. It measures nothing and its timings are never used for any timeout or ceiling (L1 is the basis). Seeds 9201 and 9202 lie outside every counted range (1 to 200, 1 to 100, 1 to 20), L1's (9101 to 9103) and the smokes'; under the default `K275_SCORE_R` the scorer's RANGE-1 stop means these rows can never enter counted scoring. `cells.tsv` is not uploaded.
**Binding.** The kit-275 tip after the auditor's closure `<commit>`; scripts' SHA-256 `<hashes>`; `K275_G1` = the path of `g1/g1.tsv` from the package root.

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

**Ceiling.** Set-up 10 + self-test 10 + A (2 + 1 kill-after) 3 + B 8 + scorer 1 5 + C 8 + scorer 2 5 = **49 job-minutes = 3.27 core-hours** (36 = 2.40 without C and scorer 2; expected about 25 = 1.67). Running ceiling: used 12.43 + 3.27 + planned 391.20 = 406.90 (the owner's figure under the 1.5 x ruling); written in `runs.tsv` before sending.
**What it cannot inform.** R, seeds, margins, levels, the outcome table, the jobs and their timeouts are fixed before it; the ceiling is decided on L1's figures; it exposes S and T rows of two non-counted data sets, used for line-count checks only.

**Also predicted (J03).** The self-test step reports 46 cases and 3,128 predictions, 0 differences; E0 reports 41 hunks, 0 unmarked. The scripts table reads: script 1 ok, script 2 (A) timeout, scripts 3 to 5 ok, script 6 (scorer 2) failed; so **the run's overall outcome is "failure" by prediction** (a failed script). A, B and C are uploaded in full; scorer 2 writes no file.
