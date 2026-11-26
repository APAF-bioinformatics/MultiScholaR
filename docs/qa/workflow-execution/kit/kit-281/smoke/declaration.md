# KIT-281 parse-and-self-test send: declaration (DRAFT; not sent; no run of any kind before the coordinator says so)

Script audit AUDIT-KIT-281-scripts-20261126, S09. The three R scripts have never met R. L1 is blind and declared once, so the first contact is a small send that settles
what cannot be settled without R: three parses, the scorer self-test (S01 and S07 cases included), the partial-match guard's first real events (`jsonlite` and `digest`
writes under the scorer's global option, `globalCallingHandlers` in a driver, bare `Rscript` through the runner's PATH), and one smoke-flagged simulation task. It decides
nothing about calibration and produces no result: the simulated rows stay in `tempdir()` and the scorer refuses them (smoke flag).

## Binding
- `package_commit`: a commit of branch `kit-281` whose `R/`, `src/`, `inst/`, `tests/`, `DESCRIPTION` and `NAMESPACE` tree equals the branch point's a5de1c2e (the shipped
  filter, arm S); recorded in `runs.tsv`. If KIT-280 has merged when the commit is chosen, bind to the last staging commit before its merge.
- `base_commit` "none", `tests` [], `blas` "reference", `extra_packages` []; one `scripts` entry `smoke/smoke-k281.R smoke-out`, `out` `smoke-out`, `out_pattern`
  `"\\.(tsv|json|txt)$"`, `timeout_minutes` 30 (draft: `request-draft.json`).
- Seed: r = 9301 (cell A, three masks), declared here: outside KIT-197 (1-200, 9001-9004), KIT-211 (201-400, 9101-9104), KIT-281's counted block (401-971) and L1 (9201-9204,
  201-205).

## Ceiling
About 20 job-minutes = 1.33 core-hours expected; `timeout_minutes` 30 = 2.00 core-hours at the ceiling. The ticket's ceiling becomes 72.67 + 2.00 = **74.67 core-hours**
(not yet accepted by the owner; no request before it is).

## What it must show (all fixed now)
1. Three parses: `sim-k281.R`, `score-k281.R`, `l1-k281.R` all `parses`.
2. `selftest`: `SELFTEST: <n> checks, 0 differ` (n = 2,007 for the cases as committed: 64 cases, 1,936 predicted values, 7 request fixtures); the runner exits 0.
3. `probe stopped as instrument`: the probe case's `outcome.tsv` reads Stop / instrument (the deliberate `$` prefix match was caught by the guard).
4. `cal105 completed`: Calibrated (proves the scorer's first writes of JSON and the digests raise no partial-match warning under the global option).
5. The simulation task: exit 0; `checks-sim.tsv` E2, E3, E7, E10 each 3 of 3; the marker assertion passed (exit 0 after `load_all`); `provenance-sim.json` smoke true and complete.
6. `scorer refuses the smoke directory`: exit 1 with "smoke flag".
`smoke-checks.tsv` lists every line with its pass flag; the run is green only if all pass.

## Output (public)
`smoke-checks.tsv`, `selftest-summary.txt`, `checks-sim.tsv`, `provenance-sim-start.json`. Nothing derived from a p-value, a call or a fold change.

## Failure handling
A failure with no result produced (a parse error, a mechanics defect in a script or in the self-test's constructed data) is an instrument defect (VAL-54): fix the mechanics, re-hash,
re-send, with no expected value, margin, rule or seed changed; each send recorded in `runs.tsv` and `outcome.md`; at most three further sends, then stop and report. A self-test
difference is first a question about the PREDICTION (made in Python from the planted counts) or the scorer, decided by the audit, not by changing a number to pass.
