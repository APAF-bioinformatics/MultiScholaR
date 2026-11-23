# KIT-177 outcome: stopped at gate F3 (outcome 3), before R1

Ticket revision 1.6 (step 0 applied 2026-11-24 in the main checkout). Branch `kit-177-prot-plain-log2` from
`v0.9-staging` fc48b6ec (KIT-180 merged). Runs on `run-tests` (R 4.6.1, reference BLAS, ubuntu-22.04); ids, commits and
job minutes in `runs.tsv`.

| Commit | Content |
|---|---|
| 025c9542 | R0 scripts only (`route.R`, `effect.R`, `replay-check.R`, `gui-entry-control.R`); `R/` and `tests/` equal to fc48b6ec |
| (records) | this file, `runs.tsv`, `feasibility.md`, `r0/` |
| (not run) | the tests commit and the fix commit, prepared and pushed but never sent to a runner (see "State of the branch") |

## Observed

**R0 (run 45222712083, 122.9 job-minutes = 8.19 core-hours; expected 8.9, ceiling 18).** Three scripts, all `ok`.

- **F4** (S2 on `prot_lfq`, head 025c9542 against base fc48b6ec): the mechanism works from a script; every comparison
  empty, as expected with equal code (`feasibility.md`, `r0/f4-*.tsv`). R4f is not needed.
- **F3** (`r0/replay.tsv`): on the MaxQuant (LFQ) and Proteome Discoverer (TMT) e2e fixtures the records replay in
  order (history equal, the DA re-run), 5 of 6 later states `identical()`; `correlation_filtered` differs in
  `@args$daAnalysisUI` only; the DA long table is not `identical()` (its column names differ between the two reads);
  the kit fingerprint finds no mismatch. **Outcome (3), "replays in order, not identical"**: the ticket's stop rule
  (revision 1.6, from AUDIT-KIT-180-impl I06) stops here for the owner.
- **9.1** (`r0/effect.tsv`, per run `r0/effect-runs.tsv`, pairs `r0/effect-pairs.tsv`, fixtures
  `r0/effect-fixtures.tsv`): 400 of 400 route runs ok; first seed 4.92 minutes, projection 103.4 of 187; done in 103
  minutes. Means over 20 seeds (range in `effect.tsv`); bias = mean of sign(truth) x estimate over changed proteins
  tested (1 = unbiased); FDP and true-positive share at BH q < 0.05 on the reported p; "low" = lowest third of true
  abundance. `old` = log2(x + 1), `new` = plain log2.

| m | RUV | bias old / new | bias, low third old / new | FDP old / new | FDP, low third old / new | TP share old / new | TP share, low third old / new | calls differing |
|---|---|---|---|---|---|---|---|---|
| 10 | automatic | 0.815 / 0.941 | 0.551 / 0.745 | 0.190 / 0.174 | 0.646 / 0.711 | 0.528 / 0.456 | 0.145 / 0.084 | 42.5 |
| 10 | skip | 0.813 / 0.941 | 0.548 / 0.750 | 0.051 / 0.041 | 0.495 / 0.418 | 0.392 / 0.357 | 0.052 / 0.031 | 9.2 |
| 40 | automatic | 0.960 / 1.002 | 0.928 / 1.024 | 0.109 / 0.106 | 0.217 / 0.221 | 0.917 / 0.902 | 0.744 / 0.705 | 18.6 |
| 40 | skip | 0.956 / 1.002 | 0.919 / 1.023 | 0.043 / 0.039 | 0.115 / 0.107 | 0.895 / 0.881 | 0.684 / 0.647 | 3.9 |
| 150 | automatic | 0.984 / 0.993 | 0.980 / 1.000 | 0.113 / 0.112 | 0.148 / 0.151 | 0.969 / 0.968 | 0.901 / 0.898 | 4.3 |
| 150 | skip | 0.983 / 0.994 | 0.978 / 1.001 | 0.048 / 0.048 | 0.078 / 0.083 | 0.983 / 0.983 | 0.947 / 0.947 | 0.9 |
| 1,000 | automatic | 0.991 / 0.992 | 0.992 / 0.995 | 0.106 / 0.106 | 0.114 / 0.114 | 0.961 / 0.960 | 0.876 / 0.875 | 1.95 |
| 1,000 | skip | 0.991 / 0.992 | 0.992 / 0.995 | 0.049 / 0.050 | 0.050 / 0.051 | 0.979 / 0.979 | 0.934 / 0.934 | 0.05 |
| 10,000 | automatic | 0.992 / 0.992 | 0.994 / 0.994 | 0.110 / 0.111 | 0.117 / 0.117 | 0.961 / 0.961 | 0.876 / 0.876 | 0.75 |
| 10,000 | skip | 0.992 / 0.992 | 0.994 / 0.995 | 0.051 / 0.051 | 0.052 / 0.052 | 0.979 / 0.979 | 0.933 / 0.933 | 0 |

  The five e2e fixtures (both configurations): no call changes; largest change in a log fold change 1.5e-5 to 1.0e-4;
  on the two five-protein fixtures RUV automatic stops ("No valid percentage found") with either transform.

Inferred, not tested: at m = 10 and 40 the offset shrinks fold changes toward zero (bias 0.81 and 0.96) and plain log2
removes it (0.94 and 1.00); at low abundance plain log2 calls fewer true positives and, with RUV automatic, more false
ones in the lowest third, consistent with the unmodelled variance at low signal (ticket section 5, limitation). From
m = 150 on the two transforms agree to within the seed range. With RUV automatic the FDP is about 0.11 at every scale
for both transforms (the legacy RUV route's own behaviour, not this ticket's).

## Acceptance criteria (status at the stop)

| # | Status | Evidence |
|---|---|---|
| A1 | F1, F2 answered by reading; F3 answered (outcome 3, stop); F4 answered (yes) | `feasibility.md`, `r0/` |
| A2 | not reached (stopped before R1) | |
| A3 | not reached | |
| A5 | not reached (NEWS, Methods and ledger rows not written) | |
| A6 | not reached | |
| A7 | 9.1's table delivered (above, `r0/effect.tsv`); to be shown to the owner; T1b not run | `r0/effect.tsv` |
| A8 | not reached (R4 not run) | |
| A9 | not reached | |

## State of the branch

Prepared on the branch and not run (no runner has seen them): the tests commit (`tests/testthat/test-prot-design-log2.R`
with T1a, T2, T3, T3b, T4, T6, T7, T8, T9; the two contract rewrites in `test-workflow-native-session.R`; T5 not
written, because it depends on F3), then the fix commit (`R/`: the design step's rule, the counts and the refusal, the
notifications, the reopen notice, the `scale` warning, the texts) and the T1b script `t1b.R`. They wait for the owner's
answer on F3.

## Not checked

- Whether F3's non-identity is the moment of capture (instrument) or a replay difference (result): not tested.
- Why the DA long table's column names differ between the first run's read and the replay's read.
- The prepared tests and fix: never run; no claim about them.
- Proteome Discoverer's default "Reporter Abundance Based On" setting (F2: not found in Thermo's pages read).
