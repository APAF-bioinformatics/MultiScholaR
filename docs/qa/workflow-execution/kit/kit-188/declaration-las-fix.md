# KIT-188: the `las` job's exact-match repair (instrument defect, VAL-54) and its non-counted rehearsal. Declared alone; nothing is sent by this commit; the auditor reads the diff first.

**Authority.** Coordinator ruling R4 (2026-11-27). Defect: counted `las` job `k188-s-las-l64-01` (run 45875346789) exited 2 after 6.6 minutes; all 12 units ended `instrument`, `partial-match guard: partial match of 'LA' to 'LAs'`. No unit has a result; nothing is judged. A `$` on `journeys` or `choices` for `LA` matches `LAs` by prefix when `LA` is absent, which is the counted `las` job's shape (`ARMS_OF$L$las` is `LAs`, `OR`; the lifecycle `las` job adds `LA` for E11, so L1r never met it).

## 1. The fix: `sim-k188-las.R`, a copy of `sim-k188.R` with eight lines changed

`diff sim-k188.R sim-k188-las.R` is eight changed lines (1634, 1688, 1694, 1695, 1701, 1702, 1732, 1733); each change is `journeys$LA` or `choices$LA` becoming `journeys[["LA"]]` or `choices[["LA"]]` and nothing else (`$LAu` and `$LAs` stay: no arm name extends them). `sim-k188.R` is untouched (SHA-256 b569d1c69b6b20e4d168ee7884d35e490e613396b8a2e4b847bd2dcbb1aa7a17); **`sim-k188-las.R` SHA-256 64e8e3149cdbbd938f3ec1ef1fb2c161be8ab60f0654db27c28385c1acf4aef3.** The script hashes itself from its own path and re-launches children with that path, so the provenance of a `las` job records the new hash and every other job keeps b569d1c6…. A list's `[["LA"]]` returns the value when `LA` is present and `NULL` when absent; it never matches by prefix.

## 2. Every site read (line numbers in `sim-k188.R`; the same lines in the copy)

Cases per site: LA present / only LAs present (the counted `las` job) / neither.

| Line | Expression | Before | After | Changed |
|---|---|---|---|---|
| 1634 | `la <- choices$LA` (arm LCm, `reported` only) | present: value; only LAs: n/a (no LAs in `reported`); neither: NULL | same value; NULL; NULL | yes (listed, never reached with LAs) |
| 1688 | `choices$LA` inside `if (r$scoreset == "T_noLAc")` | present: value; LAs only: not evaluated (the scoreset is set only for `reported` LA and LCm rows); neither: not evaluated | the same | yes |
| 1694 | `!is.null(journeys$LA) ... && !is.null(choices$LA)` (E8 and M-amp) | present: value; **only LAs: the LAs value, with the guard's warning: a stop**; neither: NULL | present: value; **only LAs: NULL**; neither: NULL | yes |
| 1695 | `checkE8(sim, journeys$LA, choices$LA, ...)` | evaluated only inside the 1694 branch | same | yes |
| 1701, 1702 | the closed-form branch; `identical(job, "deciding")` is the first operand | `deciding` always has LA; in `las` the first operand is FALSE and `$LA` is not evaluated | the same | yes |
| 1732 | E10u: `e10u && ...` (lifecycle `deciding` only) | e10u FALSE in counted runs: not evaluated | the same | yes (`journeys$LAu` unchanged) |
| 1733 | E11: `lifecycle && identical(job, "las") && ...` | lifecycle FALSE in counted runs: not evaluated; lifecycle `las` has LA: exact | the same | yes (`journeys$LAs` unchanged) |

Read and not changed (no partial match is possible): `journeys$OR <-` (1622, 1847) is an assignment; `journeys$NO$da` (1735, 1741, 1901, 1908) is lifecycle `deciding` only, where `NO` is present, and `NOnone` is a local `j2`, never stored in `journeys`; `requests$LA`, `requests$NO`, `requests$NOnone`, `requests$LM` (1625, 1736, 1842, 1902, 1945 to 1965): the lists `REQUESTS_L` and `REQUESTS_P` hold each of these names exactly, and an exact name always wins over a prefix; `journey$e1`, `journey$da` and the other `journey$` fields are the journey record's own fields (the earlier prefix collisions, `da` against `da_reached`, `e1` against `e12`, were renamed and the guard has run clean since).

**The E8 / M-amp blocks in a `las` job.** Line 1694 also requires `job %in% c("deciding")`, and line 1701 `identical(job, "deciding")`: Part M's E8 and M-amp blocks run in `deciding` jobs only, in the lifecycle as in the counted run. In a `las` job they never run, with or without LA; the repair changes only that the missing LA no longer stops the unit before the job test is reached. Not a design gap: the declared design runs them in `deciding` jobs. (`reported` runs M-df, not M-amp; unchanged.)

## 3. Provenance by job type: what the scorer and the completeness check need (nothing changed now)

`score-k188.R` reads the rows files only; it has no provenance or script-hash check (grep for `sha`, `provenance`: none). The script-hash rule lives in `declaration-counted.md` section 6 (the completeness check, read from `provenance-*.json`). It needs one amendment, no scorer change: `script_sha256` equals b569d1c6… for every counted job except the 23 `las` jobs, which equal 64e8e314…; any other pairing stops. The scorer's E5 compares the OR rows of the `las` jobs with `reported` and `lp`; the OR path is the same code in both scripts.

The 23 `las` jobs' `invocation` in `counted-jobs.tsv` (rows with job type `las`) change only the script path to `docs/qa/workflow-execution/kit/kit-188/sim-k188-las.R` (arguments, environment and timeouts unchanged); `counted-jobs.tsv` is not edited here (the amendment is made when the coordinator releases the `las` block). The failed `k188-s-las-l64-01` is re-sent under its name as the same job.

## 4. Rehearsal (NOT counted; one request, declared alone; NOT SENT)

Purpose: (i) the counted `las` configuration (no LA arm, lifecycle off) completes on the fixed script; (ii) a non-`las` unit gives the same rows under the fixed and the original script. Rehearsal seed: data set 401 (the L1r data set; never a counted data set 1 to 400).

One request `k188-rh-las-fix`, package commit = the commit carrying this declaration and the script (documents and scripts only; package tree equal to e4a6d4ae), `blas` reference, `extra_packages` `["sva"]`, `OMP_NUM_THREADS=1`, **no `K188_LIFECYCLE`, no `K188_SMOKE`, no `K188_PREDICTIONS`**; four script entries in order (the timeout of each in minutes):

| Script | Run | Env `K188_JOBS` / `K188_JOB` | Timeout |
|---|---|---|---|
| 1 | `sim-k188-las.R L M12 401 401 k188-rh-las-m12` | las / k188-rh-las-m12 | 20 |
| 2 | `sim-k188-las.R L L12 401 401 k188-rh-las-l12` | las / k188-rh-las-l12 | 30 |
| 3 | `sim-k188-las.R L M12 401 401 k188-rh-lp-fixed` | lp / k188-rh-lp-fixed | 12 |
| 4 | `sim-k188.R L M12 401 401 k188-rh-lp-orig` | lp / k188-rh-lp-orig | 12 |

Total timeout 74 minutes; **ceiling (74 + 10) x 4 / 60 = 5.60 core-hours, expected about 2.9; inside the 6.0 allowed, charged to the lifecycle line** (64.43 + 5.60 = 70.03 at ceilings).

**Predictions (counts, statuses, hashes only; no estimand is opened or printed; the written rows are unblinded because lifecycle is off, and are not read):**
1. Each of the four scripts exits 0; `units_with_error` 0 in each provenance; every unit's timing `exit` 0, `error` empty; no `instrument` row; no `partial-match guard` anywhere.
2. Rows per unit: `k188-rh-las-m12` 72, `k188-rh-las-l12` 48 (L1r's las rows; the arms `LAs` and `OR` only), `k188-rh-lp-fixed` 72, `k188-rh-lp-orig` 72.
3. Provenance `script_sha256`: 64e8e314… for scripts 1 to 3, b569d1c6… for script 4; `smoke` false, `lifecycle` false.
4. Scripts 3 and 4 (the same lp unit, M12 data set 401, under the fixed and the original script): the rows, with `seconds` and `warnings` removed, equal: `python3 -I e5d-compare.py <rh-lp-fixed rows gz> <rh-lp-orig rows gz> 401 401` prints `identical`, with the counts and the two hashes only.
5. The comparison against L1r's recorded `rows_sha256` for the same unit is **not possible and not required**: that hash is over the full row including `seconds` and `warnings` (and over the lifecycle's unblinded rows), so it changes on every run; item 4 is the equivalent check, with the original script as the reference. Said now so no one reads a missing match as a finding.

Decision: all of 1 to 4 met: the `las` block may go to the coordinator for release on `sim-k188-las.R`. Any miss is an instrument finding (VAL-54: no result is judged; a script or test mechanics failure may be fixed and re-run under the three-correction rule), reported, not worked around.
