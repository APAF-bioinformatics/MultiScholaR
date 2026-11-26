# KIT-275 smoke send: declaration (written before the send)

**Instrument check. It judges nothing and measures nothing.** One request, sent once. Its purpose is to find parse and run-time defects in the
unrun scripts (`sim-k275.R`, `score-k275.R`, `fixture-k275.R`, `selftest/`, `e0-diff.py`) before the script audit. No estimand is written anywhere
in what is uploaded: the journeys run in the blind lifecycle (`K275_LIFECYCLE=1`), their row files stay under `tempdir()`, and the fixture writes no
reading (`K275_SMOKE=1`). What leaves the runner is: checks (id, units compared, passed, exactness deviations), seconds per arm class and design, hashes,
the scorer's parse and coverage lines, the E4 self-test's text, progress files and consoles.

## Binding

- Package under test: the commit named in the request (this branch, `kit-275`), which carries the scripts. The branch is based on `origin/v0.9-staging`
  `ec27f69e0` (the bound commit of every tier but V: it holds KIT-255 with its I01 repair, merge aa119308f). `base_commit` none.
- `fetch_refs`: `kit-188-ruv-default-calibration` (E0 reads `git show 0866982f:...`).
- Not the bound commit's evidence: a smoke result is never reused as a lifecycle or counted result (ticket 8.4: L1 is run once, blind, after the audit).

## What it must show (each item is met or not; none is a result)

1. All five R scripts parse (`parse: all scripts parse`, first lines of script 1's console).
2. E0: `python3 -I e0-diff.py` reports `0 unmarked` hunks against 0866982f on the runner (the commit is fetched).
3. E4: the scorer self-test runs: 19 cases, 1,859 predictions derived in Python (`make-selftest.py`), `selftest: 1860 predictions, 0 differ`. A difference is a defect
   of the scorer or of a prediction, recorded and not repaired here.
4. G1 fixture on `X3`, data sets 9101 and 9102 (smoke flag): runs to the end, writes `g1.tsv` with `rho` NA and `smoke` TRUE and `checks-g1.tsv`; every E6 and E10
   row passes (E6: the moderated SE reading, the unmoderated reading against the kernel's sigma and stdev.unscaled, LX's reading on the hand-made matrix; E10:
   kappa_c, PR and the dot product of pi on a constructed matrix).
5. Lifecycle job A: tiers D, K, S, E, T, design D3, data set 9101, `K275_FAULT=1`, `K275_REPEAT=S`: every unit completes (no `unit_completed` check row); every
   journey of every arm class (NO, LM, LG, LA, GA, GK, GSC, PEL2, PEG2) and both LX rows end `ok` (check `L1_status`), except the deliberate fault (`E7_crash`
   passes) and a `refused_after_exclusion` row of the emulation; `E7_refusal` (GK3 at D3 on VD is refused with class `ruv_k_above_limit`) passes; `E1`, `E2`, `E3`, `E5`,
   `E9` (tier S), `E11`, `E12` pass; the blind scorer exits 0 and writes coverage for every arm class.
6. Lifecycle job B: tiers D and ID, design D4, data set 9101: the metabolomics identity journeys (NO, LM2, GK2, LA, GA) complete; the blind scorer's `E8` line
   compares 5 rows (a difference would be a finding, not a stop, at counted scale; here it is recorded).
7. `SMOKE_counted_refusal`: script C (`K275_SMOKE=1`, no `K275_LIFECYCLE`) exits 1 within seconds with that message before any unit.
8. No `partial-match guard` stop in any console, progress file or table (lesson 2).
9. Provenance and progress exist: `provenance-<job>.json` (written before the first unit, again at the end; `"smoke": true`; script, scorer and declaration SHA-256),
   `checks-<job>-start.tsv`, `progress/checks-*.tsv` for every unit.
10. Timing by arm class and design (`timing-<job>.tsv`) is recorded for the audit's information; it sets nothing (G2 rests on L1).

## The request (one; scripts run in order from the package root)

| # | Script | Arguments and environment | Limit (min) |
|---|---|---|---|
| 1 | `docs/qa/workflow-execution/kit/kit-275/selftest/run-selftest.R` | `K275_SELFTEST_OUT=k275-smoke-selftest` | 10 |
| 2 | `.../fixture-k275.R` | `D3 9101 9102 k275-smoke-g1`; `K275_SMOKE=1` | 4 |
| 3 | `.../sim-k275.R` | `D,K,S,E,T D3 9101 9101 k275-smoke-a`; `K275_LIFECYCLE=1 K275_SMOKE=1 K275_FAULT=1 K275_REPEAT=S K275_JOB=smoke-a OMP_NUM_THREADS=1` | 22 |
| 4 | `.../sim-k275.R` | `D,ID D4 9101 9101 k275-smoke-b`; `K275_LIFECYCLE=1 K275_SMOKE=1 K275_JOB=smoke-b OMP_NUM_THREADS=1` | 10 |
| 5 | `.../sim-k275.R` | `D D3 9101 9101 k275-smoke-c`; `K275_SMOKE=1 K275_JOB=smoke-c` (refused) | 2 |

Four workers per sim job (the runner's four cores). Planning basis (ticket section 10): Automatic 22.8 s, Manual 17.8 s per journey; job A holds 47 journeys
and one refused journey (largest unit: tier K, 17 journeys, about 6 minutes of wall), job B 19.

## Ceiling, charge, decision

- Ceiling: (22 + 10 + 4 + 10 + 2 + 10 set-up) = 58 job-minutes x 4 / 60 = **3.87 core-hours**; expected about 31 job-minutes = 2.07 core-hours. Charged to L1's line (`runs.tsv`).
- Running-ceiling rule (ticket 10 (2)): used 0.00 + this ceiling 3.87 + everything still planned 321.21 = 325.08 core-hours, at most 350.73: met (`runs.tsv`).
- Attempts: one send. If it shows a defect, the defect is recorded in `lifecycle.md` with the file and line and **nothing is repaired and nothing is re-sent**; the coordinator rules
  on a second send. If it shows none, the scripts go to the script audit unchanged.
- No result of the measurement is seen by anyone: none is written. A number in the uploaded tables other than a check deviation, a second count or a hash is itself a defect.
