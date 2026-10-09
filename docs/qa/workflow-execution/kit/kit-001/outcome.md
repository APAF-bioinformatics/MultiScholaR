# KIT-001: full-journey driver, per-step populations and step map

Status 2026-10-09: harness built and verified by its tests and validate-only passes; the two counted lifecycle members are declared below and wait for the user's go. Authority: `tickets/kit/00-KIT-PLAN.md` run ledger (KIT-001 row). Protocol v1 frozen (`../protocol.md`, sha256 in `../bindings.json`).

## What was built (new files only; CORE-pinned files byte-identical)

| File | Role |
|---|---|
| `tools/profiling/workflow_kit.R` | bindings, fixed cases, selectors, journey step map, phase measures (`phase_peak_anon_bytes$<step>`, `resident_end_bytes`), eligibility and lifecycle decisions, store and file fingerprints, independent verifier |
| `tools/profiling/workflow_kit_driver.R` | the journey: import, design (any number of groups), the 7 peptide and 5 protein QC actions, normalisation with RUV (`ruv_pearson` profile), correlation, session export, DA; per-step `start`/`complete`/`failed` marks and a `teardown` mark in `journey-phases.tsv`; per-step store snapshot and fingerprints; the forced-failure injection |
| `tools/profiling/workflow_kit_coordinator.R` | stage coordinator: frozen source snapshot, bindings and hash checks, host precondition, managed `gui` launch from a fresh launcher on the case library, 50 ms slice sampler, browser driven from the coordinator (outside the slice), operation receipt, verifier launch, decision |
| `tools/profiling/run_workflow_kit.sh` | guarded `validate` / `run` entry (clean tree, 12 GiB, no other managed unit, no live stage) and the `--sample` loop with the 15.5 GiB / 2 GiB / deadline guards |
| `docs/qa/workflow-execution/kit/bindings-generator.R`, `bindings.json` | core_018 build, copied CORE populations and oracles, the two corpora, frozen journey settings, lifecycle fault, harness and protocol/step-map hashes |
| `tests/testthat/test-workflow-kit-harness.R` | 51 assertions: step map, case expansion and digests, selectors, phase measures and resident end, lifecycle decisions, setup refusal, path masking |
| `kit-001/step-map.json` | 19 entries (18 journey steps plus enrichment), units determined from the kernel code with file:line evidence, global statistics, merge rules, holders, mode; 0 open flags |
| `kit-001/holder-inventory.json` | KIT-012 worklist: 121 holder sites in 40 files, 49 holding a full S4 or a DA list with `theObject`, 19 in-memory `saveState` calls, 26 `.GlobalEnv` assigns (KIT-013) |

## Corpora (bound with sha256)

- Screen: 8 samples of the PERF-107 calibration recipe (same universe, proteins and cutoff; counts 500, 2999, 3000, 3001, 5000, 12000, 40000, 60000), 133,892 report rows, generated once under `<proof>/kit/corpus/screen`. The sample-filter cutoff 3000 retains 6 samples (2 per group).
- Full: PERF-107 `calibration` (300 samples, 5,846,578 rows); `heldout` stays reserved.

## Journey settings frozen in the bindings (audit A13)

GUI defaults with `diagnostic_settings = FALSE`, except the sample-filter cutoff (3000, the corpus's designed cutoff); normalisation `cyclicloess`, RUV automatic (1–20 %, max_difference, penalty 0.5, adaptive, grouping `group`), correlation threshold 0.5; DA `~ 0 + group`, q 0.05, treat lfc 0.

## Claim limits (audit A14) and one constraint for KIT-005

1. **How the on-disk owner is entered.** The production `gui` launch has no switch for the on-disk adapter (`multischolar.dia_native_sessions` is an R option; the controller runs `Rscript --vanilla`). The only existing launch-side path is the PERF-107 phase-A member record (`programme perf107-phase-a`, variant `physical`, entry `resident`, no fixture test): the in-memory `WorkflowState` holds every state until the sample filter, where the history is converted into a `.dia-gui-adoption` store and the sample filter runs through the managed worker. The import to intensity states are fingerprinted after adoption (converted component checkpoints), the later ones after their own step; the record marks each with `snapshot_time`.
2. **Controller deadline.** That member record caps the controller at 1800 s. The 300-sample journey ceiling (KIT-005, scaled from the screen duration) must stay at or below 1800 s, or the user decides a different entry before KIT-005 launches. Recorded in `bindings.json` `journey.controller_deadline_seconds`.
3. **Production default journey.** Without the record, nothing is on disk per step; KIT-011 removes this limit.

## Step map: determinations (no open flag)

| Unit | Steps |
|---|---|
| store | import, design, session export |
| complete sample | P14 rollup (after the global annotation/evidence lookups and validation), sample filter (after design validation and the inclusion-absence set) |
| feature group | q-value (protein; global row offsets), min peptides per protein (protein), intensity filter (protein+peptide after the global quantile pre-pass), replicate filter (protein+peptide), imputation (protein+peptide within technical-replicate groups; deterministic, no random draws), IQ (protein; surviving-sample union), duplicate removal (protein) |
| row | accession cleanup (DIA branch annotates `protein_id_table` only), protein replicate filter |
| whole matrix | protein intensity filter (dataset-wide quantile), normalisation + RUV, correlation, DA, enrichment |

Two defects in the existing partitioned reference path are recorded for KIT-009/KIT-014 (intensity helper aborts on absent design runs; non-additive summary fields copied from the first group). Five science observations are listed for the user; none changes a determination.

## Declared runs (ledger row KIT-001; waiting for the user's go)

| Member | Stage | Build / corpus | Ceiling | Decision |
|---|---|---|---|---|
| lifecycle complete | `kit-lifecycle-complete` | core_018 (`browser-smoke.gQ40C7`) on the screen corpus | 300 s | `lifecycle_complete_observed`: every step `complete`, a finite anonymous peak per step, `resident_end_bytes` observed, every state fingerprinted, verifier recomputes and matches |
| lifecycle fail | `kit-lifecycle-fail` | same | 300 s | `lifecycle_failure_observed`: guard injected after the sample filter (`application-guard-inject.rds`), the monitor stops the slice (`injected_guard_stop`), the replicate step `failed`, head equal to the sample-filter snapshot, slice stopped |

Commands (from the repository root, tree clean, bindings committed):

```
tools/profiling/run_workflow_kit.sh run kit-lifecycle-complete /home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/core/browser-smoke.gQ40C7/library
tools/profiling/run_workflow_kit.sh run kit-lifecycle-fail     /home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/core/browser-smoke.gQ40C7/library
```

Risk stated before launch: the 300 s ceiling is the ledger's; the screen journey (import, 12 QC actions, RUV search, DA with plots) has not been timed through this driver. A time guard would be a valid recorded result, not a retry; a demonstrated harness fault allows one correction.

## Evidence

Stage receipts, member receipts, traces, journey records and verifications land under `<proof>/kit/runtime/<stage>/`; the decision in `decision.rds`. This file records the outcome after the runs.
