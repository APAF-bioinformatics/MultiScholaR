# KIT-001: full-journey driver, per-step populations and step map

Status 2026-10-09: harness built and verified by its tests and validate-only passes; three user decisions applied the same day (launch hook, RUV replicate structure, quantile ceiling); the control build `kit_001` and the two counted lifecycle members are declared below and wait for the user's go. Authority: `tickets/kit/00-KIT-PLAN.md` run ledger (KIT-001 row). Protocol v1.1 (`../protocol.md`, sha256 in `../bindings.json`).

## What was built (new harness files; the four CORE-pinned files byte-identical)

| File | Role |
|---|---|
| `tools/profiling/workflow_kit.R` | bindings, fixed cases, selectors, journey step map, phase measures (`phase_peak_anon_bytes$<step>`, `resident_end_bytes`), eligibility and lifecycle decisions, store and file fingerprints, independent verifier |
| `tools/profiling/workflow_kit_driver.R` | the journey: import, design (any number of groups), the 7 peptide and 5 protein QC actions, normalisation with RUV (`ruv_pearson` profile), correlation, session export, DA; `start`/`complete`/`failed` marks per step and a `teardown` mark in `journey-phases.tsv`; per-step store snapshot and fingerprints; the forced-failure injection |
| `tools/profiling/workflow_kit_coordinator.R` | stage coordinator: frozen source snapshot, bindings and hash checks, host precondition, managed `gui` launch with the KIT journey record from a fresh launcher on the case library, 50 ms slice sampler, browser driven from the coordinator (outside the slice), operation receipt, verifier launch, decision |
| `tools/profiling/run_workflow_kit.sh` | guarded `validate` / `run` entry (clean tree, 12 GiB, no other managed unit, no live stage) and the `--sample` loop with the 15.5 GiB / 2 GiB / deadline guards |
| `docs/qa/workflow-execution/kit/bindings-generator.R`, `bindings.json` | builds, copied CORE populations and oracles, the two corpora, frozen journey settings, lifecycle fault, harness, protocol and step-map hashes |
| `tests/testthat/test-workflow-kit-harness.R` | 51 assertions: step map, case expansion and digests, selectors, phase measures and resident end, lifecycle decisions, setup refusal, path masking |
| `kit-001/step-map.json` | 19 entries (18 journey steps plus enrichment): units determined from the kernel code with file:line evidence, global statistics, merge rules, holders, mode; 0 open flags |
| `kit-001/holder-inventory.json` | KIT-012 worklist: 121 holder sites in 40 files, 49 holding a full S4 or a DA list with `theObject`, 19 in-memory `saveState` calls, 26 `.GlobalEnv` assigns (KIT-013) |

## Three user decisions applied (2026-10-09)

1. **Launch hook (R/ change, bounded).** `launchWorkflowApplication` reads `gui-member.rds` through `workflowGuiLaunchRecord`; programme `kit-journey-v1` sets `config$journey` (native DIA sessions switched on in the controller before `runApp`, deadline 60 s to 90 min) and never a qualification, so `workflowQualificationVariant` stays `legacy`. The journey therefore runs the production native-session path from import onward (`initializeDIANativeImport`: on-disk reference states, in-process reference kernels, the GUI rollup on its reference path). Test: `tests/testthat/test-workflow-execution-launch-record.R` (16 assertions).
2. **RUV-III replicate structure.** `resolveProtRuvReplicateBlocks` chooses group × replicate blocks when the design has technical replicates (any duplicated block), else the RUV grouping variable; `ruvIII_C_Varying` builds M from it and records `ruv_replicate_structure` in its args. Control selection (ANOVA by group) and the k search keep the biological grouping. Test: `tests/testthat/test-prot-ruv-replicate-structure.R` (15 assertions).
3. **Quantile ceiling.** Interrogation: the ceiling has been present since the first S4 pass (fa80eb90, 2025-04) with no stated rationale; on the protein stage the values are log2 (after IQ), so `ceiling()` lifted the threshold by up to one log2 unit (up to twice as strict as the percentile); on the linear peptide scale the effect is negligible. Removed at the five proteomics sites (`func_prot_qc_peptide_methods.R` ×2, `func_prot_qc_filtering_helpers.R`, `func_prot_s4_qc_filter_methods.R`, `utils_workflow_dia_qc_rollup.R`). Three expectations in `test-prot-qc-peptide-methods-direct-shared.R` encoded the ceiling (2 for 1.03, 3 for 2.5) and were corrected. The lipidomics and metabolomics intensity filters keep theirs until their stage (KIT-031/032).

Consequence: the KIT control is no longer the `core_018` install but `kit_001`, built from the core_018 source plus these changes; protocol v1.1; KIT-005 gains one production-default baseline journey member (no record, in-memory state) as the memory baseline (ledger 35).

## Corpora (bound with sha256)

- Screen: 8 samples of the PERF-107 calibration recipe (same universe, proteins and cutoff; counts 500, 2999, 3000, 3001, 5000, 12000, 40000, 60000), 133,892 report rows, under `<proof>/kit/corpus/screen`. The sample-filter cutoff 3000 retains 6 samples (2 per group).
- Full: PERF-107 `calibration` (300 samples, 5,846,578 rows); `heldout` stays reserved.

## Journey settings frozen in the bindings (audit A13)

GUI defaults with `diagnostic_settings = FALSE`, except the sample-filter cutoff (3000, the corpus's designed cutoff); normalisation `cyclicloess`, RUV automatic (1–20 %, max_difference, penalty 0.5, adaptive, grouping `group`), correlation threshold 0.5; DA `~ 0 + group`, q 0.05, treat lfc 0.

## Claim limits (audit A14)

- The journey's state-manager mode is the native DIA session (reference states from import), entered through the launch record; every step is snapshotted after its own step.
- A production-default journey (no record) has no per-step references; it is measured for memory only (KIT-005's extra member).

## Step map: determinations (no open flag)

| Unit | Steps |
|---|---|
| store | import, design, session export |
| complete sample | P14 rollup (after the global annotation/evidence lookups and validation), sample filter (after design validation and the inclusion-absence set) |
| feature group | q-value (protein; global row offsets), min peptides per protein (protein), intensity filter (protein+peptide after the global quantile pre-pass), replicate filter (protein+peptide), imputation (protein+peptide within technical-replicate groups; deterministic), IQ (protein; surviving-sample union), duplicate removal (protein) |
| row | accession cleanup (DIA branch annotates `protein_id_table` only), protein replicate filter |
| whole matrix | protein intensity filter (dataset-wide quantile), normalisation + RUV, correlation, DA, enrichment |

Two defects in the existing partitioned reference path are recorded for KIT-009/KIT-014 (intensity helper aborts on absent design runs; non-additive summary fields copied from the first group).

## Declared work waiting for the user's go

1. **Install `kit_001`** (recorded, not counted): `run_workflow_checkpoint_proof.sh --service install-corrected` on the clean tree at the commit carrying the three changes; then `bindings-generator.R kit_001=<install run dir>` and a commit of the bindings.
2. **Lifecycle members** (ledger row KIT-001):

| Member | Stage | Build / corpus | Ceiling | Decision |
|---|---|---|---|---|
| lifecycle complete | `kit-lifecycle-complete` | `kit_001` on the screen corpus, journey record | 300 s | `lifecycle_complete_observed`: every step `complete`, a finite anonymous peak per step, `resident_end_bytes` observed, every state fingerprinted, verifier recomputes and matches |
| lifecycle fail | `kit-lifecycle-fail` | same | 300 s | `lifecycle_failure_observed`: guard injected after the intensity filter (the first post-rollup step), the monitor stops the slice (`injected_guard_stop`), the sample-filter step `failed`, head equal to the intensity-filter snapshot, slice stopped |

```
tools/profiling/run_workflow_kit.sh run kit-lifecycle-complete <kit_001 library>
tools/profiling/run_workflow_kit.sh run kit-lifecycle-fail     <kit_001 library>
```

Risk stated before launch: the 300 s ceiling is the ledger's; the screen journey has not been timed through this driver. A time guard would be a valid recorded result, not a retry; a demonstrated harness fault allows one correction.

## Evidence

Stage receipts, member receipts, traces, journey records and verifications land under `<proof>/kit/runtime/<stage>/`; the decision in `decision.rds`. This file records the outcome after the runs.
