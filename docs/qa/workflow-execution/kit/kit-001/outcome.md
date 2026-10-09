# KIT-001: full-journey driver, per-step populations and step map

Status 2026-10-09: **complete.** Lifecycle pair observed on the `kit_001` control build with screen corpus v2: `lifecycle_complete_observed` and `lifecycle_failure_observed`. Authority: `tickets/kit/00-KIT-PLAN.md` (run ledger row KIT-001, decision rows of 2026-10-09). Protocol v1.1 (`../protocol.md`, sha256 in `../bindings.json`).

## Results

| Member | Decision | Elapsed | Whole-slice peak | Anonymous peak | Resident end | Sampling gap | Evidence |
|---|---|---|---|---|---|---|---|
| `kit-lifecycle-complete` | `lifecycle_complete_observed` | 269 s (ceiling 300 s) | 2.76 GiB | 2.15 GiB | 1.86 GiB | 70 ms | [lifecycle-complete.json](lifecycle-complete.json) |
| `kit-lifecycle-fail` | `lifecycle_failure_observed` | 87 s | — | — | — | 70 ms | [lifecycle-fail.json](lifecycle-fail.json) |

Per step of the complete member (anonymous peak of the whole slice during the step):

| Step | s | GiB | Step | s | GiB |
|---|---|---|---|---|---|
| import | 9.0 | 1.08 | iq | 8.9 | 1.41 |
| design | 13.9 | 1.24 | accession | 6.1 | 1.40 |
| qvalue | 13.7 | 1.28 | protein_intensity | 8.9 | 1.43 |
| protein_peptide | 9.4 | 1.30 | duplicates | 6.1 | 1.45 |
| rollup | 9.4 | 1.35 | protein_replicate | 6.9 | 1.45 |
| intensity | 16.5 | 1.35 | normalisation_ruv | 70.2 | 2.15 |
| sample | 8.1 | 1.33 | correlation | 8.9 | 1.86 |
| replicate | 8.4 | 1.35 | session_export | 11.2 | 1.86 |
| imputation | 14.7 | 1.39 | da | 11.4 | 1.99 |

The failure member completed import to the intensity filter, injected the application's own guard (`application-guard-inject.rds`), the monitor stopped the slice (`injected_guard_stop`), the sample-filter step failed on disconnection, and the store head equalled the intensity-filter snapshot (`intensity_filtered`, head `d2c71061…`). Zero swap, slice stopped, in both members.

The complete member's store: 16 states (`initial` through `correlation_filtered`), 31 products (QC plot products and rendered assets, the DA fit, the DA result references), 527 publication files; every step has an in-run snapshot (head sha256, current state, entity id, rows) and the verifier recomputed and matched every fingerprint.

## Launch record (the ledger allowed one correction; the user authorised attempts 3 and 4)

| Attempt | Cause | Classification | Fix |
|---|---|---|---|
| 1 | the coordinator did not export the browser executable for chromote | instrumentation | `run_workflow_kit.sh` exports `CHROMOTE_CHROME` (`a3b40b45`) |
| 2 | native-path defect: the per-group reference loop handed the intensity helper the full design and it refused groups lacking runs (155 of 2,500 proteins present in every run after the rollup) | application defect, observed correctly | `require_complete_design_runs` flag, per-group result shown identical to the whole-table result on the member's data (`9075b361`); control reinstalled as `kit_001` (`browser-smoke.y8rDOi`) |
| 3 | screen corpus v1 too sparse: the group-count intensity filter left 365 of 60,000 features, the 3,000 cutoff then removed every sample and the GUI carried an empty table to IQ | fixture | screen corpus v2 (`ea731185`); GUI gap recorded for KIT-013/014 |
| 4 | journey complete; the driver's fingerprinting overran the C stack on the RUV diagnostic plot product | instrumentation, after a valid measurement | masking depth bounded (`99487853`); fingerprints recomputed from the immutable store by [attempt4-post-hoc-fingerprints.R](attempt4-post-hoc-fingerprints.R) with the in-run snapshots kept; the state-attach slip found by the verifier fixed (`dfd6aacd`); verifier rerun from the corrected harness; driver originals kept beside the member (`*.driver.rds`, `verification.driver-failed.rds`) |

The failure member ran once and verified on its first launch (harness at `3512f84d`: the masking fix was in, the state-attach fix `dfd6aacd` came after; the failure decision does not depend on it).

Provenance of the attempt-4 correction (audit PA-02): `journey$correction$harness_commit` records `3512f84d`, but the script and driver as run were the versions committed afterwards at `dfd6aacd` (the attach-loop fix was applied before the correction ran and committed after). The successor script [attempt4-post-hoc-correction-v2.R](attempt4-post-hoc-correction-v2.R) records the full correction (fingerprints, verifier rerun, decision) as it ran.

Install sequence (audit PA-04): `kit/browser-smoke.sLJ4aC` failed at configure (the exact-sum artifacts live under `core/`); `core/browser-smoke.uvS55B` (source `d2c73bfc`) was bound as `kit_001` and superseded by the native-path fix; `core/browser-smoke.y8rDOi` (source `9075b361`) is the bound `kit_001`.

Claim scope of the forced failure (audit PA-03): the guard was injected and the sample step marked `start` at the same uptime; the step failed 0.27 s later on disconnection before any sample-filter work, so the member proves the application's guard path and head stability across a stop between steps, not during a step's write. KIT-005's operation-member lifecycle uses a mid-step fault.

## What was built (new harness files; the four CORE-pinned files byte-identical)

| File | Role |
|---|---|
| `tools/profiling/workflow_kit.R` | bindings, fixed cases, selectors, journey step map, phase measures (`phase_peak_anon_bytes$<step>`, `resident_end_bytes`), eligibility and lifecycle decisions, store and file fingerprints, independent verifier |
| `tools/profiling/workflow_kit_driver.R` | the journey: import, design (any number of groups), the 7 peptide and 5 protein QC actions, normalisation with RUV (`ruv_pearson` profile), correlation, session export, DA; `start`/`complete`/`failed` marks per step and a `teardown` mark; per-step store snapshot and fingerprints; the forced-failure injection |
| `tools/profiling/workflow_kit_coordinator.R` | stage coordinator: frozen source snapshot, bindings and hash checks, host precondition, managed `gui` launch with the KIT journey record from a fresh launcher on the case library, 50 ms slice sampler, browser driven from the coordinator (outside the slice), operation receipt, verifier launch, decision |
| `tools/profiling/run_workflow_kit.sh` | guarded `validate` / `run` entry (clean tree, 12 GiB, no other managed unit, no live stage, browser present) and the `--sample` loop with the 15.5 GiB / 2 GiB / deadline guards |
| `docs/qa/workflow-execution/kit/bindings-generator.R`, `bindings.json` | builds (`kit_001` control, `core_018` superseded), copied CORE populations and oracles, the two corpora, frozen journey settings, lifecycle fault, harness, protocol and step-map hashes |
| `tests/testthat/test-workflow-kit-harness.R` | 56 assertions |
| `kit-001/step-map.json` | 19 entries: units from the kernel code with file:line evidence, global statistics, merge rules, holders, mode; 0 open flags |
| `kit-001/holder-inventory.json` | KIT-012 worklist: 121 holder sites in 40 files, 49 holding a full S4 or a DA list with `theObject`, 19 in-memory `saveState` calls, 26 `.GlobalEnv` assigns (KIT-013) |

## Three user decisions applied (2026-10-09)

1. **Launch hook (R/).** `workflowGuiLaunchRecord`: programme `kit-journey-v1` sets `config$journey` (native DIA sessions switched on in the controller, deadline 60 s to 90 min), never a qualification. Test: `test-workflow-execution-launch-record.R`.
2. **RUV-III replicate structure.** `resolveProtRuvReplicateBlocks`: group × replicate blocks when the design has technical replicates, else the grouping variable; recorded in `args$ruvIII_C_Varying$ruv_replicate_structure`. Test: `test-prot-ruv-replicate-structure.R`.
3. **Quantile ceiling.** Removed at the five proteomics sites (it lifted the log2 protein threshold by up to one unit); three ceiling-encoding test expectations corrected; lipidomics and metabolomics keep theirs until their stage.

Plus the native-path fix found by attempt 2 (item 2 of the launch record). The KIT control is `kit_001` (core_018 source plus these); protocol v1.1; KIT-005 gains one production-default baseline journey member (ledger 35).

## Corpora (bound with sha256)

- Screen v2: 8 samples of the PERF-107 calibration recipe (same universe, proteins and cutoff) with counts 3,000 to 60,000; 320,872 report rows; the 3,000 cutoff drops one sample. v1 (counts 500 to 60,000; 133,892 rows) stays on disk as the attempt-3 record.
- Full: PERF-107 `calibration` (300 samples, 5,846,578 rows); `heldout` stays reserved.

## Journey settings frozen in the bindings

GUI defaults with `diagnostic_settings = FALSE`, except the sample-filter cutoff (3,000); normalisation `cyclicloess`, RUV automatic (1–20 %, max_difference, penalty 0.5, adaptive, grouping `group`); correlation threshold 0.5; DA `~ 0 + group`, q 0.05, treat lfc 0. State-manager mode: native DIA session from import (reference states), entered through the journey record.

## Step map: determinations (no open flag)

| Unit | Steps |
|---|---|
| store | import, design, session export |
| complete sample | P14 rollup (after the global annotation/evidence lookups and validation), sample filter (after design validation and the inclusion-absence set) |
| feature group | q-value (protein; global row offsets), min peptides per protein, intensity filter (protein+peptide after the global quantile pre-pass), replicate filter, imputation (deterministic), IQ (surviving-sample union), duplicate removal |
| row | accession cleanup (DIA branch), protein replicate filter |
| whole matrix | protein intensity filter (dataset-wide quantile), normalisation + RUV, correlation, DA, enrichment |

Recorded for later tickets: the sample-filter summary merge (non-additive fields from the first group; KIT-009/014) and the empty-assay GUI gap (KIT-013/014).

## Hand-off

- **KIT-005.** The screen journey took 269 s at 8 samples with 70 s in normalisation/RUV; the 300-sample ceiling is set from this with the scaling recorded before launch (never above 90 min; the journey record allows up to 5,400 s). Populations and oracles come from the counted control journey on `kit_001`; the production-default baseline member measures memory only.
- **KIT-012.** `holder-inventory.json`. **KIT-013/014.** The empty-assay gap and the summary merge rule.
- Claim limit: the per-step fingerprints of the complete member were recomputed post hoc from the immutable store after the driver's fingerprint failure; the in-run snapshots are the driver's. The next journey member (KIT-005) exercises the corrected in-run path.
