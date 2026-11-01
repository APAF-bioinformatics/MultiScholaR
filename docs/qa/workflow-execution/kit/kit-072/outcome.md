# KIT-072 outcome: one path in the lipid and metab modules

Status 2026-11-01: **completed**. Stage 3b, ticket 1 of 10.
- **Code:** a93d4c14.
- **Bindings and records:** 99e2ddb3.
- **Evidence:** this commit, on `performance-v2`.
- **Base:** 97aa1861 (the stage-3b audit's repairs).
- **Install:** `kit_072` = `core/browser-smoke.kJ0lsx` (a93d4c14), over the kit_032 library.
- **Ticket:** `tickets/kit/KIT-072.md`, revision 1.1.

## What changed

**1. One owner path per module site.**
- **Lipid sites.** The double-only branches are gone from:
  - intensity, duplicates and QC finalise;
  - the normalisation action (ITSD, log2, between-sample, RUV);
  - correlation apply and skip;
  - DA run and DA session restore.

  The summary's record and bundle add-ons now run whenever the bus exists.
- **Metab sites.** The same applies to:
  - intensity and duplicates (the apply step only publishes);
  - QC finalise;
  - the ITSD, log2, between-sample and RUV steps (one `metab_ruv`);
  - correlation dispatch, outcome and skip;
  - DA run and restore;
  - design import.
- **Removed helpers.** Five helpers only the double path used are removed: `prepareMetabDuplicateResolutionState`, `saveMetabIntensityFilterState`, `runMetabNormRuvOptimizationStep`, `runMetabNormRuvCorrectionStep` and `runMetabNormApplyCorrelationWorkflow`. Grep before and after: [branch-grep-before.txt](branch-grep-before.txt), [branch-grep-after.txt](branch-grep-after.txt). The after-list holds only refusals and the norm exports, which were already one path.

**2. The shared helpers refuse a non-owner by name** (`saved_entity_scope`): `markWorkflowOmicsStoreStep`, `recordLipidDesignProvenance`, `recordMetabDesignProvenance`, `bindWorkflowOwnerTypeAtDesign`, `ingestWorkflowOmicsLegacySession`, `resolveWorkflowOmicsOwnedSession` and `workflowSummaryFinalStateName`.

**3. Tests.**
- **Lanes.** The twelve lipid/metab module-CI lanes run on real owners (`omicsTestOwner`, in `helper-workflow-omics-owner.R`). The lanes' injected kernels are bound in the namespace the declarations call, so their data is unchanged.
- **Module seams.** These run on real owners, with the declared wrapper recorded (`omicsMockWrapper`).
- **Contract changes.** The pass-through tests became refusal tests.
- **Quarantine.** Five tests of the removed helpers are quarantined. [affected-tests.txt](affected-tests.txt) lists them, together with the two dropped double-only assertions and the migration patterns.

**4. Bindings.** The generator binds the stage-3b exactness keys `0NNl`/`0NNm` and `KIT_0NN_LIFECYCLE` for KIT-072 to KIT-081, as one loop.

## Results

| Gate | Evidence | Result |
|---|---|---|
| Affected tests | [affected-tests.txt](affected-tests.txt) | 139 files, one process per file: **10,036 pass, 0 fail, 0 error**, 25 skips (5 quarantined here). The first run on the code change alone had 104 failing tests in 37 files, all of them tests of the removed path. |
| Module-CI | [module-ci-comparison.txt](module-ci-comparison.txt) | Every data product equal (bytes or content). The import, design, QC and DA lanes are byte-identical. The norm lanes differ only in store- and time-bearing files and the ITSD timestamp. `study_parameters.txt` now shows the owner's rendering. |
| R CMD check | [check-item1.log](check-item1.log) | Baseline: 3 WARNINGs, 3 NOTEs; the code note is identical to KIT-032's. |
| Exactness set (1): the build's records | [exactness-gates.txt](exactness-gates.txt) | Lipid and metab recorders on kit_072: 10 states each, records valid, no legacy step. Bound as `lipid72`/`metab72`. |
| Exactness set (2): the members | [kit072-concurrent-kit-journey-072l-stage.json](kit072-concurrent-kit-journey-072l-stage.json), [-072m](kit072-concurrent-kit-journey-072m-stage.json) | **`kit-journey-072l` accepted, exact** (273 s, 2.55 GiB peak). **`kit-journey-072m` accepted, exact** (273 s, 2.61 GiB). No sampler gap above 250 ms. |
| Exactness set (3): no change across builds | [exactness-gates.txt](exactness-gates.txt) | `lipid72` against `lipid31`: **identical**, 10 steps, 0 mismatches. `metab72` against `metab32`: **identical**, 10 steps. |
| Exactness set (4): records | [member-record-replay-072l.json](member-record-replay-072l.json), [-072m](member-record-replay-072m.json) | Both records valid, 10 steps, no legacy step. Both replays from the design state are **exact** (10 states), and a tampered input is refused. |
| Exactness set (5): lifecycle | [lipid](kit072-concurrent-kit-lifecycle-lipid-stage.json), [metab](kit072-concurrent-kit-lifecycle-metab-stage.json) | Both observed. Peaks: lipid 2.53 GiB, metab 2.61 GiB, equal to the kit_032 members'. |
| Exactness set (6): proteomics | [exactness-gates.txt](exactness-gates.txt), [replay-check-kit072.json](replay-check-kit072.json) | KIT-028 replay **exact**, 16 states, 279 s. LFQ (11), TMT (11) and DIA_LIMPA (14) records **identical** to the kit_065 oracles. |

## Deviations and discoveries (recorded)

- **The four members ran concurrently.** The runner returns at launch, and its "no active unit" precondition did not block the next three launches. It was not enforced by me.
  - Each member ran in its own slice with no sampler gap above 250 ms.
  - Exactness does not depend on concurrency.
  - The lifecycle peaks equal kit_032's to 0.01 GiB.
  - Elapsed time rose by about 12% (273 s against 244 s), which is reported only.
  - No rerun was spent. The attempt is kept as the evidence, with stages labelled `kit072-concurrent-*`.
- **Acceptance item 2 is met with a scoped divergence.** Tests whose subject code has no owner branch keep a reading stand-in that offers only `getState` or `setWorkflowType`: import type recording, ITSD candidate reads, QC tab renders, norm seams 03c and the audit carry-overs. A stand-in cannot reach a removed branch. Two summary-context tests keep a `getState` stand-in for an undeclared characterization class, with `workflowSummaryFinalStateName` recorded. Listed in [affected-tests.txt](affected-tests.txt).
- **Discovery, left untouched.** On the owner, `study_parameters.txt` prints `@args` read back from the store with 17-digit doubles (`q_cutoff = 0.10000000000000001`) and lists as JSON. This is production's path since KIT-070. It became visible now that the lanes run on the owner.
- **Discovery.** The recorders write their runtime project under `runtime/oracles/lipid31_lipid_screen` and `metab32_metab_screen` whatever the output file. A later recorder run therefore overwrites the store an earlier binding's populations point at. The comparison is by fingerprint, so no gate depends on it; KIT-070 recorded the same.

## Backward compatibility and rollback

Production behaviour is unchanged. Production always took the owner's path, and every data product is equal. No saved format changed, and proteomics is unchanged (gates above). Rollback: `git revert` the ticket's commits.

## Handoff

KIT-073 (the shared declaration core) now has one path per module site. From here the exactness set's references are `lipid72`/`metab72`, which equal `lipid31`/`metab32`.
