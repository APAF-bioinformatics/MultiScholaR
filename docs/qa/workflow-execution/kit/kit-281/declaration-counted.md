# KIT-281 counted run: declaration (committed ALONE, before any request; one run, scored once)

Authority: ticket KIT-281 r1.1 sections 5 to 10; owner's acceptance of the ceiling 74.67 core-hours and the replicate rule (VAL-371); script audit and Re-checks 1 and 2 (AUDIT-KIT-281-scripts-20261126, RC2-4, RC2-5); coordinator's rulings. L1 (run 45761773898) and the parse send (45759150472) are recorded in `l1/`, `smoke/`, `runs.tsv`.

## Binding
- `package_commit`: the commit of branch `kit-281` that carries this declaration and `g3-rule.tsv` (named in each request). Its `R/`, `src/`, `inst/`, `tests/`, `DESCRIPTION` and `NAMESPACE` tree equals a5de1c2e's (`git diff --stat a5de1c2eb <package_commit> -- R/ src/ inst/ tests/ DESCRIPTION NAMESPACE` is empty, checked). G0: `git diff --stat 894b9cc8 <package_commit> -- R/func_prot_qc_filtering_helpers.R R/func_prot_s4_qc_filter_methods.R R/func_general_parameter_helpers.R` is empty (checked; L1's G0 1/1). The gate must not include KIT-280: `threshold_values_in_filter` is absent from `R/` (grep: no hit) and every simulation process asserts it absent at load.
- Unchanged from L1 (SHA-256): `sim-k281.R` 10fa5bfda5cb05c3ee84def2fd528b64a2bb84b4c51ce242731f5eb0b8b4a236; `score-k281.R` 0f838f47e3eb21bf4fe389e6e776c0639fba751f73880d0e589e4b3ddf6d8628; `l1-k281.R` 7d5d4cb295c4ea3968d33e93ea737cd662ad1efead1accfc7beb7e31fdb88c32; `check-request-k281.py` efc28b36da40746e3adcc00134bd5f8cac95935e5e12b95c73759e5a5232ee9d.
- `g3-rule.tsv` (the path the scripts read) is a byte copy of `l1/g3-rule.tsv`: SHA-256 9568f599fccf85442a49da738fc4f2ebe25bbccc4c31a607aff2b84ff19d8497 (every request cites its first 16 hex digits, 9568f599fccf8544). Rule G3: N-hat = 175,600 -> R = 315 for every cell, r = 401 to 715.
- Requests: `base_commit` "none", `tests` [], `blas` "reference", `extra_packages` [], one `scripts` entry each, `out` `k281-out`, `out_pattern` `"\\.(tsv|gz|json|txt|md)$"`.

## Layout (coordinator's ruling): eight jobs, one cell per request, split evenly 158 + 157
Per cell: job 1 `sim-k281.R k281-out 401 558` (158 replicates), job 2 `sim-k281.R k281-out 559 715` (157), `env` `{"K211_CELLS": "<cell>"}`, plus `"K211_DETERMINISM": "1"` on cell A's first job. Every scripts entry carries `timeout_minutes` 70 (G2's plan from L1's seconds: 46.0 and 45.8 job-minutes, 1.5x). `python3 -I check-request-k281.py --g3 g3-rule.tsv --count l1/g3-count.tsv <the eight request files>` must exit 0 before the push (the checker needs every cell covered, so the eight together); outputs recorded in `runs.tsv`.

## Ceilings
Expected about 24.48 core-hours (4 x (55.6 + 36.2) min is the 200 + 115 layout's figure; this layout's expectation is the same, 24.48). Counted ceiling 8 x 70 min = 560 min = 37.33 core-hours (4.67 each). **Running-ceiling line before EACH request (all eight identical, the earlier jobs at their ceilings):** used 2.28 + counted 37.33 + scoring 2.67 + continuation reserve 5.00 (75 min) + re-scoring reserve 2.67 + Part C 16.00 = **65.95 <= 74.67**.

## The decision as coded (fixed; nothing changes after any data exist)
VAL-365 Q5 (a): rho_g = alpha_B1(H_B1 without H_B) / alpha_B(H_B) pooled over the twelve settings (A to D at vd5, vd10, vd25); Calibrated needs its bootstrap 0.975 quantile (type 7, B = 10,000, seed 280000001, unit = a cell's replicate) at most 1.10 and the guards G1a (FDP_BH(B1) <= 0.05), G1b (paired difference <= 0.0025), G2 (paired BH power cost <= 0.005); F5' (the pooled two-sided 95 % interval of the fold-change bias of flagged one-sided proteins in H_B1 excludes 0) gives Withheld; outcome order Withheld, Calibrated, Not calibrated. Boundary rule: a limit passes iff its computed double is <= the double of the written limit; no tolerance. VAL-371 rulings: a non-finite ESTIMATE of rho_g or a guard is a Stop (instrument); the q-method fallback count is reported, not a stop. The jackknife-normal limits, Phi, the positive control, the per-design, per-setting and by-flag tables are reported and decide nothing.

## Stops
F1 (an identity of F1 violated), F2 (E3, E10 or determinism), and every instrument status: incomplete or duplicated rows, a row outside r 401 to 715, a smoke flag, `instrument-stop.txt` or `errors.txt` in a job folder, a missing determinism table, a non-finite estimate, a G3 table not as declared, the continuation's second failure. A Stop claims nothing; the default stays 2.

## Continuation (once per job; VAL-54), declared verbatim before the send
- A job that ends without `provenance-sim.json` (timeout, lost runner) is "did not complete" for its range; no row of it is scored.
- If its uploaded `out` holds `parts/<cell>-<r>.metrics.tsv` files, those parts are committed at `docs/qa/workflow-execution/kit/kit-281/cont/<cell>-<first>/parts/`. One continuation request runs `sim-k281.R` with `out` = that folder, the same `<first r> <last r>` and the replicate list of the missing tasks only (the simulation keeps existing parts and assembles all).
- If nothing was uploaded, the continuation runs the whole range once. No row of it existed, so this is the range's first and only result.
- In every case each replicate is scored exactly once (the scorer stops on a duplicate). The continuation's ceiling comes from the reserve (75 minutes = 5.00 core-hours), and it is recorded in `runs.tsv`.
- A second failure of the same range stops the gate: Stop, nothing claimed, the default stays 2. No partial cell is ever re-run and no completed replicate is ever replaced.

## Handling of the data
Only completeness is read before scoring: units present per cell (counted from the uploaded tables' headers/lengths by file listing and `range.tsv`), statuses (exit, `provenance-sim.json` status, smoke), hashes, job-minutes. No row or metric is opened, printed or computed on. The scorer is NOT run until the coordinator says so.

## Predictions of the CHECKS only (never the outcome)
Every job: exit 0, `range.tsv` complete = expected (158 or 157 per cell), `provenance-sim.json` status complete, smoke false, counted true, no `instrument-stop.txt`, no `errors.txt`; the determinism table (cell A's first job) 5 of 5 identical. At scoring: `complete` 1 (15,120 metrics rows: 4 cells x 315 replicates x 3 masks x 4 arms; erratum of the first commit's 3,780, a slip in the count only), `ranges_incomplete` 0, F1 0, F2 0, `determinism` 5/5, `bootstrap_B` 10000, `range_declared` g3-table, q-method count reported. No prediction of the outcome beyond the ticket's section 10 expectation (guards pass; rho_g near the margin, resting mostly on cell A).

## Scoring plan
Once, after every counted job is complete (the scorer's completeness check refuses less), on the eight job folders, `Rscript docs/qa/workflow-execution/kit/kit-281/score-k281.R <rows root> <out>` with no environment variable (B = 10,000, range from the table), `outcome.tsv` must read `counted` TRUE, `timeout_minutes` 40 (expected a few minutes; the measured seconds are recorded). No re-scoring except an instrument defect (VAL-54, reserve 2.67).
