# KIT-159 outcome (running record; build phase)

Ticket: `tickets/kit/KIT-159.md` revision 1.7 (1.6 plus the third re-check's one note, applied as step 0:
the label of 2.1's decision rule; no contract change). Branch `kit-159-validation-harness` from
`v0.9-staging` at a60c4e9a (version 0.9.1). Worktree `/home/doktersmol/Projects/APAF/MultiScholaR-kit-159`.

## Scope of the build phase

Steps that need only the code (section 11; the auditor's T0, T1, T2, R1): the harness `msrval` with its
unit and classification tests (tests first, shown failing on a runner: T0; then passing: T1, with
`covr`, `lintr` and the f_sim timing), the coverage regimes T-COV (T2, sized on the measured f_sim),
and the first record request (R1). Not in this phase: G1D (another agent, branch `run-k159-g1d`),
G1(a) at the release commit, G1(r), S0, L1a, L1b and later, P1 (they wait for KIT-172, KIT-176,
KIT-177, KIT-180 and KIT-182).

## Work items, in order (each: tests first, then implementation)

The harness is built in three module groups so that each group's tests are shown failing on a runner
before its implementation (no R runs on the laptop, so a single tests-only commit for every test of
section 10 would be one unreadable mass failure). T0 and T1 are each split across these groups; their
jobs are counted against T0's and T1's ceilings (6 core-hours each), not against new budget lines.

| Group | Files | Tests |
|---|---|---|
| a. scoring | `score.R`, `metrics.R`, `intervals.R`, `baseline.R` | 10.1 (all hand cases), T-VIEW, T-KEY, T-BLIND-1 |
| b. records and guard | `inputs.R`, `licence.R`, `declaration.R`, `provenance.R`, `determinism.R`, `evidence.R`, `units.R` | T-IN, T-LIC, T-LOG, T-DECL, T-PRE, T-ENV, T-SCORE, T-ANN, T-PROV, T-DET, T-CHUNK, T-RUN, T-RNG |
| c. driver | `driver.R`, `arms.R`, `parity.R`, `integration.R` (Part B later) | T-DRV, T-HELP, T-CHAIN, T-KER, T-TAB |

## Log

### Step 0 (2026-11-23)

The third re-check's note applied in the main checkout's `tickets/kit/KIT-159.md`: 2.1's heading
"**Decision rule** (the register's, unchanged)" replaced by "**Decision rule** (the round-2 re-check's,
with re-check2 C13's conditions)"; revision 1.7; closure section 22; history line. The rule's text is
unchanged.

### G1D recorded (2026-11-23; not run by this agent)

Run 45121359461, script commit c9c4a412 (branch `run-k159-g1d`), about 1.1 core-hours. Tables and
script copied to `g1d/` (aggregates, column names and digests only, checked before commit); entry in
`gates.md` under G1(r). Observed: the same DA table on 28d5a4be and 12d68ca7 on the runner; head
reproduction holds, base reproduction fails (8a4c2dd5 on the runner against the oracle's 80cfedaa),
so `limpa_da_last_digits` is not admitted on G1D alone. `prot_dia_limpa`: awaiting D19 (owner).

### T0a (group a tests only): run 45122558749 at 525f9733

Observed: `msrval` (no R code) installs; all 53 group (a) tests error with "could not find function"
(`t0/t0a-harness-tests.tsv`, `t0/t0a-harness-failures.txt`). Scripts job 13.2 min, 13.7 job-minutes
with the results job: 0.91 core-hours of T0's 2 expected / 6 ceiling. The script's `out` comes back in
the run's `scripts` artifact (not `results`): `gh run download <id> -n scripts`.

Test corrected after T0a (recorded): T-M3 now passes `n_min = 1`, as 10.1's preamble states for every
hand case (the T0a version omitted it, so the 10-feature floor would have made the value NA).
`test-coverage.R` (pass bands, the regime summary, small regime runs) was written with the coverage
code and has no failing-first run.

### T1a (group a implementation): run 45124309425 at 5ab037ae

Observed: 57 tests, 49 pass (`t1/t1a-harness-tests.tsv`, `t1/t1a-harness-failures.txt`). Seven hand
cases (T-PT, T-JK-1, T-JK-2, T-JK-crit, T-CP, T-ONE, T-REF-FILL) fail with "actual != expected but
don't know how to show the difference": the tests used a relative tolerance (1e-6 or 1e-9) tighter
than the ticket's printed rounding (7 to 10 decimals), so values equal to print precision failed. One
test (small C1 run) recycled a filtered vector against the full level column. Inferred: test defects,
not implementation defects (the waldo message means the difference is below print precision; the
correction run shows the values). `covr` is not in the scripts job's library (it needs
`extra_packages`, as `lintr` does), so the script stopped before covr, lintr and f_sim.

Correction (T1a's one correction pass): the hand-case tests compare within one unit of the ticket's
last printed decimal (`expectClose`, which prints the actual value); the C1 test fixed. covr and lintr
move to the final T1 part over the complete package (A2 is over the whole `msrval`).

T0b (group b tests only, 51ed29f1) was sent before this; its result is recorded when it returns.

### T0b (group b tests only): run 45124916341 at 51ed29f1

Observed: 93 tests; all 36 group (b) tests fail: 23 with "could not find function", the others with
the missing internal constants `HARNESS_PATH` and `EVIDENCE_COLUMNS`; T-DET-1 also had a test defect
(`data.table::fwrite` has no `digits` argument; corrected after T0b, the test now writes with
`fwrite`'s round-trip form). The git-repository fixtures ran (git is in the job). T0 total: 1.46 of 6
core-hours.

## Next

T1a correction run (k159-t1a2-5e1c3202); then group (b) implementation. Former plan for T1a: the group (a) implementation (`utils.R`, `score.R`, `metrics.R`, `intervals.R`, `baseline.R`,
`coverage.R`), tests passing, `covr`, `lintr`, f_sim on regime C3 and a cost probe of every T-COV
family (T2 is sized on these).

## Not checked

(filled as work proceeds)
