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

## Next

Group a: tests-only commit, then the T0 request for it.

## Not checked

(filled as work proceeds)
