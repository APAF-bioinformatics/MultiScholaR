# KIT-010 pre-execution remediation (2026-10-13)

Register: [KIT-010-PREEXEC-AUDIT-20261013.md](KIT-010-PREEXEC-AUDIT-20261013.md). User: "repair everything first".

| ID | Repair | Closure |
|---|---|---|
| PE-01 | None (the user's decision at "go"; the ledger row records the one-build/two-arm pair at launch). | Ledger at launch. |
| PE-02, PE-04, PE-05, PE-06 | Accepted without change (recorded). | — |
| PE-03 | The worker's protein compute takes the experiment paths from the request's publication context when the controller wrote one beside the request (`workflowRequestExperimentPaths`), so the state configuration records the project's protein QC output path as the GUI does; without a context (a bare script) the working-directory fallback stays. Unit test added (`test-workflow-operation-registry.R`). Code changed in `R/`, so the tree was reinstalled and rebound. | Test green; install and binding below. |
| PE-07 | The KIT-009 memory-attribution probe relaunched (both runs); its record lands beside the KIT-009 script when it ends. | `kit-009/memory-attribution-probe.json` when written. |
| PE-08 | Test drift from accepted changes repaired: the two accession-cleanup tests now put the FASTA table on the workflow bus (KIT-013 moved the step off the global environment); the dia-state adapter check compares the sample filter's arguments order-insensitively (since KIT-008 the adapter restores the support entries into the arguments after the finishing step's bookkeeping; same names, same values). `test-prot-02c` 651, `test-workflow-checkpoint-dia-state.R` 236, `test-module-ci-prot-protein-qc.R` 65 green. Left as a recorded discovery: MCI-007.6 fails inside `plotDensity` because a ggplot2 object now carries the S7 class `ggplot2::ggplot`, which the method's `ggplot` signature no longer matches; a package compatibility change outside any KIT ticket. | Tests; discovery recorded. |

Install after the remediation: `kit_010` = `core/browser-smoke.CZRMll` from `8cf56195` (empty runtime patch), bound; the earlier install `browser-smoke.IaLhTL` (acc8687e) stays recorded, never run as a member; stage `kit-010` validates.
