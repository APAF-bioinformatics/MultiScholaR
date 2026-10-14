# KIT-019: DA publication through the managed worker

Status 2026-10-14: **complete** (no counted member, as the ledger says; one install recorded: kit_019 = `core/browser-smoke.cDoGhF`, `aa2c3031`, bound as KIT-017's candidate; the control stays `kit_018`). Authority: `tickets/kit/00-KIT-PLAN.md` (decision row 2026-10-14, KIT-016 post-hoc audit PH-01; run ledger KIT-019: none, one install recorded); protocol v1.3. Control: `kit_018` (`core/browser-smoke.yXPrcl`). Hand-offs: KIT-012 (Glimma's namespace tree, 0.44 GiB in the application process) and KIT-016 (the pages are saved assets and the render path loads no Glimma, but the DA publication built them in process).

## What was built (one mechanism)

1. **`da_publication`, a registered worker operation** (`R/utils_workflow_operations.R`; input `checkpoint`, class `ProteinQuantitativeData`, compute role `protein`, completion and publication roles `da`, one parameter `annotations`: the annotation table as a saved reference, or absent). The worker handler `computeWorkflowDaPublicationAction` (`R/utils_workflow_execution.R`) reads the accepted manifest from the head snapshot (`da_results`), refuses an input that is not the head or a manifest of another source, materialises the annotation table (an empty table when absent, as the in-process path), and runs the one publication body `produceDIADaPublication` (`R/utils_workflow_dia_da.R`: the annotated exports, the static volcano and NumSig figures through the unchanged `outputDaResultsAllContrasts` kernel, the interactive Glimma pages) into `request directory/da-publication/<role>`; the publication kernels resolve their tables through the session's active manager, so the worker runs them inside a mock session scope that names its owner. The worker saves no asset and publishes no head (the controller holds the head): its value is the head checkpoint and reference plus the inventory. Completion (`completeDIADaPublication`) commits nothing new and republishes the head unchanged so the controller's acceptance names the head state; publication (`publishWorkflowDaResult`) moves the staged files into the project through the owner's action-output journal (`publishDIAActionOutputs`, the same identity as the in-process path) and saves each contrast's page as the KIT-016 rendered asset. The DA kernel (limma, the KIT-004 contrast declaration) stays in process.
2. **The DA module requests it** (`R/mod_prot_da_handlers_run_analysis.R`): under a managed application the handler runs the DA step with the publication deferred (`runDIAReferenceDaAnalysis(..., publish = FALSE)`), saves the annotation table as a display table and applies the managed step (`installWorkflowManagedStep`, KIT-008/014: request, poll, publish, browser acknowledgement on the status text and the static volcano plot's identity, Revert cancels, the result read by `readManagedDaResult`); `deliver` hands the manifest to the module's holders (`completeDIADaDisplay`) and updates the contrast selectors. Without a managed application the in-process path publishes as before. The status text is one helper (`diaDaStatusText`) shared by the render and the acknowledgement.
3. **Probe** ([da-publication-probe.R](da-publication-probe.R) → [da-publication-probe.json](da-publication-probe.json)): on a reflinked copy of the KIT-008 screen member's project the application-side process runs the DA kernel with the publication deferred; a separate R process runs the registered handler on the same store; the application-side process completes and publishes the worker's files (the controller's two phases) and records whether Glimma loaded in it, its anonymous memory before and after, the published files and the page asset shown without a widget rebuild.

## Tests (working tree, proof library on `R_LIBS`)

| Suite | Result |
|---|---|
| test-workflow-dia-products.R | 92 pass (new: the worker handler stages the exports, figures and pages under the request directory with the annotations joined, saves no asset and leaves the head unchanged, refuses a non-head input; the completion republishes the head and names the head state; the controller's publication moves the staged files into the project and saves the page assets; `readManagedDaResult`; the worker entry point from a frozen request) |
| test-workflow-operation-registry.R | 93 pass (new: the `da_publication` entry and its parameter validator) |
| test-module-ci-prot-da.R | 245 pass (the DA render mock widened for the plot identity) |
| test-prot-07b-da-handlers-compat.R | 55 pass |
| test-prot-07-da-analysis.R / test-prot-08-volcano.R | 3 (1 skip) / 54 (1 skip) pass |
| test-prot-da-compat-shared.R / test-prot-da-render-direct-shared.R | 16 / 37 pass |
| test-workflow-partitions.R | 294 pass |
| test-perf107-audit.R | 193 pass |
| test-workflow-execution.R / test-workflow-prot-declarations.R | 64 (10 skips) / 112 pass |
| test-workflow-kit-harness.R | 326 pass |
| test-perf-003.R | 137 pass, 1 pre-existing error ("registered dispatch" needs the installed runtime; baseline) |
| test-workflow-native-session.R | 500 pass, 1 pre-existing failure+error ("native publication recovers after process loss": needs an installed package under load_all; baseline) |

## Probe

Record at `aa2c3031` on the KIT-008 screen project (8 samples): [da-publication-probe.json](da-publication-probe.json).

| Measure | Application-side process | Worker process |
|---|---|---|
| Glimma namespace loaded after the DA step | **no** (seven namespaces added by the limma kernel: BiocParallel, corpcor, ellipse, MASS, matrixStats, mixOmics, rARPACK) | yes |
| Anonymous memory | 0.613 GiB before the DA step, 0.677 GiB after kernel, completion and publication (+0.063 GiB) | 0.955 GiB at the end of the publication |
| Time | kernel 9.5 s; completion and publication 1.8 s | publication 9.7 s (10 files, 4.6 MB) |
| Project outputs | the two annotated exports, the three volcano files, the four NumSig files, the interactive page: `completed_outputs$da` TRUE | — |
| Interactive volcano | shown from the saved asset in 0.1 s, no widget rebuild | — |

Against the KIT-012 record (Glimma's tree: 0.44 GiB anonymous, 0.30 GiB live heap in the application process), the DA step now leaves 0.063 GiB in the application-side process on the screen project.

## Claim limits

- The probe exercises the handler, the completion and the publication in process and in a plain child process, not the systemd worker and monitor units of a managed application; the production dispatch of `da_publication` through those units is first observed by KIT-017's journey pair (the README pre-launch chain items 1 to 6 are covered: registry entry, controller compute helper, worker entry point from a frozen request, the adapters vehicle's input check (component checkpoint), completion and publication on the fixture, the probe on the real screen project).
- No counted member; the resident-memory judgement stays with KIT-017's `resident_end` rule.
- The in-memory owner path (no managed application) publishes in process as before, Glimma included.
