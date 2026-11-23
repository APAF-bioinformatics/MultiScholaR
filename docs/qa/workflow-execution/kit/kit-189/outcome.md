# KIT-189 outcome (in progress: stopped at gate G2)

Ticket revision 1.2. Branch `kit-189-qvalue-vs-bh` from `v0.9-staging` at 3ff60190448150fa13e06d624378cad1481d3820.

## State

L1 is complete and blind; G1 (the mechanism) passed; G2 (time) stops the work: re-split by its rule, the counted run needs about 239 core-hours expected and 447 at the ceiling against the accepted 151 and 281 (D3). Nothing is declared and no counted job has run. Details and the re-split table: `lifecycle.md`.

## Observed

- Runs 45275836889 (L1, commit 426cebb0) and 45281701515 (E7 instrument run 1 of 3, commit 9f78d703): `runs.tsv`.
- E0 to E4, E6, E7 pass (`checks.tsv`). The mechanism of 2.1 holds on the runner: the installed enrichit, qvalue and clusterProfiler functions are the fetched sources (E1); on every L1 call the q-value is one constant times the BH-adjusted p-value over the kept terms and the repeated chain reproduces it exactly, and the BH calls lie within the q-value calls (E2); KIT-173's 120-term list reproduces exactly, with pi0_hat 0.0716 from enrichit's first stage (E0).

## Inferred

- Measured time per replicate on four workers is about 1.6 times the ticket's basis (B1) and B3 is 1.7 times B1 rather than 1.5 (`lifecycle.md`).

## Acceptance criteria

| # | State | Evidence |
|---|---|---|
| A1 | met so far | `git diff --stat 3ff60190..kit-189-qvalue-vs-bh` lists only `docs/qa/workflow-execution/kit/kit-189/` |
| A2 | L1 complete and blind; E0 to E4, E6, E7 pass; G1 applied (pass); G2 applied (stop) | `lifecycle.md`, `checks.tsv`, runs 45275836889, 45281701515 |
| A3 to A8 | not reached | — |

## Not checked

- Everything after G2: the declaration, S, D (E5), scoring, the outcome and the readings, `results.md`, the ledger row.
- The scoring script's outputs (run blind on L1's rows; only their hashes are recorded).
