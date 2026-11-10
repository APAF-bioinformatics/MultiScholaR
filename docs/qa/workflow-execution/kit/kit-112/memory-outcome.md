# KIT-112 D8: the memory and time measurement (2026-11-10)

**Decision: "better and no worse". The claim is made.**
- On the 300-sample DIA corpus, the application's (controller's) `da` peak fell from 2.208 to 1.575 GiB, against a limit of ≤ 1.958 for "better".
- Every "no worse" rule holds on the large pair and on the screen members.
- Every result is identical: 16 states, 105 store comparisons including the `da_results` manifest, 0 different. The 18 verifier mismatches are the same.
- No contingent row is raised (P-36, P-38).

All members are recorded and uncounted, with one attempt each under the host rule (`power-state.txt`). Every member is valid: largest sample interval ≤ 0.26 s, swap 0, slice stopped. GiB = 2^30 bytes. Figures are the controller's anonymous memory per step unless named. Per-step figures for every pair are in `memory/steps.tsv`.

## The large pair (111q in process on `kit_111`; 112q on the worker on `kit_112`)

- **Null check first:** `compare-stores.R` on `111p` against `111q` (both in process): 16 states, 105 comparisons, **0 different** (`memory/compare-111p-111q.txt`). That count includes the seven new `da_results` manifest rows.
- **Exactness:**
  - Both verifier decisions are `rejected` on `oracle_mismatch`, recorded verbatim. The mismatch set is identical: the same 18 entries, KIT-065's set. The `oracle`, `populations` and canonical-form blocks are equal too.
  - Store against store (`memory/compare-111q-112q.{tsv,txt}`): 16 states, 105 comparisons, **0 different**. The manifest rows (present, contrasts, parameters, model, source, contrast names, the result `grouptreated-groupcontrol`) are all identical.
- **Dispatches:** 112q = 111q's 16 + `DA-contrast-output` 1/1. `da_publication` stays 1/1.

| Rule | 111q | 112q | Limit | Result |
|---|---|---|---|---|
| **Better:** `da` controller peak | 2.208 | **1.575** | ≤ 1.958 | **met** |
| No worse, application: every step's peak | — | every step lower (−0.050 to −0.633) | +0.10 | met |
| No worse, slice: whole-member anonymous peak | 5.813 | 5.492 (−0.321) | +0.25 | met |
| No worse, time: `da` | 100.9 s | 82.9 s (0.822×) | ≤ 1.10× | met |
| Reported: `da` slice peak | 3.461 | 2.832 | — | — |
| Reported: `report` peak | 2.319 | 1.889 | — | — |
| Reported: `resident_end_bytes` | 2.643 | 2.333 | — | — |
| Reported: whole journey | 2,132 s | 2,270 s (1.065×) | not judged | — |

**What the fall is made of (the traces, `candidate-1/memory-trace.tsv`, 50 ms samples):**
- **111q (in process):** the controller sits at 1.70 GiB at the click. It rises to 2.03 within 5 s and to 2.20 at about 30 s, where it holds to the step's end. That is +0.50 GiB from the kernel, as W0 found.
- **112q (on the worker):** the controller is flat at 1.555 GiB through the whole step (1.555 to 1.570). The slice rises from 1.92 to 2.67 while the worker runs the kernel, falls back to 1.92 when the worker exits, then briefly reaches 2.81 at about 79 s (the publication's figures). The controller's own rise in the step is **+0.015 GiB**, against **+0.500** in process.
- **A run-level offset:** 112q's controller is lower than 111q's from the first step, `import` −0.05 and `design` −0.36, and every step before `da` is 0.15 to 0.37 lower. KIT-112 changes none of those steps. This is the same kind of release floor, varying between runs, that KIT-111 saw on 111m (`kit-111/memory-outcome.md`), here in the favourable direction.
- **So the claim is stated on the step's own rise, not on the whole −0.633:** the move removes about 0.49 GiB from the controller's `da` step (+0.500 → +0.015 above `session_export`). "Better" also holds with the offset added back: 112q's `da` peak plus the `session_export` gap (1.575 + 0.148 = 1.723) is still ≤ 1.958.

**Time:**
- `da` is faster (100.9 → 82.9 s).
- The whole journey is 1.065× slower. That is in early steps KIT-112 does not change: `qvalue` 178 → 272 s, `protein_peptide` 121 → 207 s, `rollup` 44 → 73 s. Later steps (`imputation`, `replicate`) were faster.
- It is not judged, and it is not explained by the ticket's change, which runs only at `da`.
- The host during the window: AC power, the low-power profile (unchanged), load 1.94 at the start and 1.32 at the end. No agents ran and no R runs were made by the lead; the lead only wrote the chain C script.

## The screen members (against `111l/m/t`, reused)

The lead's decision before any run: the screen's "no worse" rule judges the `da` step only. As it turned out, every step is within +0.10.

| Member | Exact | `DA-contrast-output` | `da` peak | Largest step rise | Slice | `resident_end` | `da` seconds | Elapsed |
|---|---|---|---|---|---|---|---|---|
| 112l | yes | 1/1 | 0.657 → 0.657 | +0.003 | 1.749 → 1.751 | 1.101 → 1.079 | 14.3 → 19.8 | 139.7 → 165.6 s |
| 112m | yes | 1/1 | 1.018 → 0.983 | +0.017 | 2.096 → 2.062 | 1.366 → 1.329 | 16.0 → 22.6 | 250.5 → 300.5 s |
| 112t | yes | 1/1 | 0.660 → 0.658 | +0.005 | 1.750 → 1.749 | 1.099 → 1.078 | 13.3 → 18.4 | 136.3 → 167.0 s |

Screen times are reported, not judged. The `da` step's +5 to +7 s is a worker start-up on a small fixture, as `correlation` showed in KIT-111.

**The other omics (`112a/b/c` against `111a/b/c`):**
- exact, with equal dispatches;
- largest step change +0.014, +0.016, +0.000;
- slice 1.834 → 1.836, 1.853 → 1.859, 1.928 → 1.919.

## Lifecycle (M2) and the failure half

`112f` failed at `da` as declared:
- `DA-contrast-output` 1/0;
- no `da_publication`;
- head `correlation_filtered`;
- `da_results` NULL.

See `proof/chain-b.md`.

## Information

- **The verifier's file count** is +1 on 112l, 112t and 112q (255 → 256, 514 → 515) and equal on 112m (326).
- **The extra file** is one more native export: the worker's copy of the DA action's input (`.workflow-export/native-export-*/files`, 59 → 60 on 112t). It comes with one more managed action's records.
- **112m has no extra file:** that input was already exported earlier in its journey.
- **Spread between builds:** the count moved by ±1 between builds before (108t 256 → 111t 255).
- It is a workflow record, not a result. The stores are identical.

## Contingent rows

- **P-36 (not raised):** "better" is met.
- **P-38 (not raised):** no "no worse" rule fails. The correction pass was not used.
