# KIT-121: a deterministic memory release after import and design (2026-11-13)

**Outcome:** done, after one correction pass. The app's session owners carry out the existing memory reclaim after a completed import, design, QC or normalisation. No result changes. Authority: the user, 2026-11-13 (autonomy for the GUI journey; this changes no scientific output).

## Why

KIT-120's diagnosis: the application freed the imported table (about 0.7 GiB on the DIA full corpus) only when R's next full collection happened to fall after it became garbage, so its memory after import depended on unrelated code. The reclaim mechanism existed (the capacity trial built and tested it) but had been off since KIT-011/012, by scope, with no recorded hazard.

## Changes

| | Change | Commit |
|---|---|---|
| W1 | `reclaim_enabled` on the session's first owner (`createWorkflowSessionOwner`) | `57661c2b` |
| Correction pass | The switch moved into `attachDIAStateSession`, through which every app owner passes. The import replaces the first owner (`bindDIASessionOwnerType`), so W1's switch never reached the owner the import and design use. | `30514baf` |

Tests: `test-workflow-omics-owner.R` 79/0 (new; the attached-owner case fails on W1's code); `test-workflow-gui-memory-lifetimes.R` 81/0 (one assertion pinned the old default, off, and now checks that a switched-off owner refuses a request); the other reclaim and session files and the four entry lanes at KIT-120's exit counts.

## Measurement (criteria declared before any run, `tickets/kit/KIT-121.md`)

**First, on kit_931 (W1):** the release missed (the after-runs equal the before-runs, 2.65 against 2.62 GiB); stores identical; screen journeys exact. A diagnosis build (`kit_932`, side branch `kit121-diag` `b60d18ef`, messages only) showed `request-call import enabled=FALSE`: the cause above.

**Repeated with unchanged criteria on kit_933 (the correction):**

| Criterion | Result |
|---|---|
| Short runs (stopped after design), release ≥ 0.5 GiB below the before-build (923h first, 923i last: 2.74, 2.61 GiB) | 933u 1.21, 933v 1.23, 933w 1.96 GiB released: **met** |
| Short runs' design step ≤ 1.10× | 0.95×: met |
| Short-run stores, 923h against 933u | 2 states, 15 comparisons, 0 different: met |
| Full journey 933q, stores against 923q | 16 states, 105 comparisons, 0 different: met |
| Full journey, application level over `qvalue` to `duplicates` ≤ 922q + 0.10 (1.71 GiB) | 1.53 GiB (922q 1.61, 923q 1.92): **met** |
| Full journey time ≤ 1.10× 923q | 0.994×: met |
| Screen journeys exact (933m DIA_LIMPA, 933l LFQ, 933t TMT, 933a lipid, 933b metab) | all accepted: met |

**Reported, not judged:** from `normalisation_ruv` on, 933q's step peaks are about 0.1 GiB above 923q's (normalisation 2.64 against 2.55, correlation 1.86 against 1.75, report 2.24 against 2.16), while the middle steps are about 0.4 GiB lower. The release happens; this is later allocation, not retained import data.

## Notes

- The harness's `reclaim_memory` launch option is dead (nothing in `R/` reads it): for KIT-117.
- The diagnosis branch `kit121-diag` (`b60d18ef`) is deleted; its build stays in the bindings as a recorded member.
- Evidence: `exit/journeys/` (decisions, store comparisons, member results).
