# KIT-009: feature-group partitioner

Status 2026-10-12: **accepted and complete.** Stage `kit-009` decision `accepted` (user decision: one build, two arms; user amendment: both members rerun after the attempt 1 harness failure); `kit_009` (`core/browser-smoke.s292Da`, source `a7a46909`) is the control from KIT-010 onward. Authority: `tickets/kit/00-KIT-PLAN.md` (run ledger KIT-009: control and candidate on the imputation population, 10 min each); protocol v1.3. Control: `kit_011` (`core/browser-smoke.AFaqC1`). Predecessor cases: PERF-113, PERF-114, PERF-116; audit findings A1, B22, B23.

## Unit gate (simulation setup 1)

`kit-001/step-map.json` (bound by hash): `protein_peptide`, `intensity`, `replicate` and `imputation` carry unit `feature_group` with their grouping column and code evidence; `open_flags` is 0 and none of the four carries a flag. The intensity filter's dataset-wide quantile is its declared `global_statistics` entry; the imputation's `source_row_order` is the writer's global ordering.

## What was built

1. **`partitionByFeatureGroup(scope, columns, group, kernel, consume, worker_cap, working_bytes_per_row, overhead_bytes, policy)`** in `R/utils_workflow_partitions.R`: batches of whole groups from the group cursor (`openWorkflowTableGroups`, which now exposes the directory's `group_sizes`), the kernel once per batch, `consume` per kernel value in batch order; the admission estimate is the working set of the largest batch the cursor can deliver (its delivery bound is 10,000 rows of whole groups, or one larger group alone, never above 64 MiB) from the declared working bytes per hydrated row plus the step's fixed overhead, refused above the cap and recorded. No scientific table, column or operation is named in the file (test). Divergence from the ticket's sketch: the batches come from the existing cursor's delivery bound, not from a directory-driven count chooser, so the batches a step takes are the ones the per-group loop already took (the admission can only refuse, never enlarge a batch beyond the codec's bound); the ticket's `declaration$merge` rules are declared on the operation and applied by the family's existing writer, which already merges by those rules.
2. **The peptide-group family on the partitioner.** `runDIAReferencePeptideGroups` (GUI entry, contract unchanged) now calls `computeDIAReferencePeptideGroups(manager, reference, section, parameters, worker_cap, variant, state_name, description)`, a compute-only function with the finishing step's contract: the kernel body and the writer body of the loop are its `kernel` and `consume` closures. Variant `partitioned` (production) runs them through `partitionByFeatureGroup` with the engine bounded by the admitted cap (an eighth, at most the loop's 1 GB) instead of the fixed `SET memory_limit = '1GB'`; variant `grouped` is the loop as it stood before KIT-009 (fixed 1 GB engine limit, no admission), kept as the control arm. The intensity filter's quantile stays the one global pre-pass before the batches (audit A1). The q-value step and the sample filter are not in this family (ticket item 2). A section computed outside a session and without a global configuration keeps a private configuration copy seeded from the saved state (the worker process; before, the loop required the GUI's `config_list`).
3. **Imputation as a worker operation** (audit B23): `diaPeptideImputationOperation()` (unit `feature_group`, parameters `proportion_missing_values`, `exclusion_column`, merge rules design/concatenate/sum, kernel = the compute-only entry point) registered as `peptide_imputation` (handler `computeWorkflowPeptideImputationAction`, validator `validateWorkflowPeptideImputationParameters`, step label `8_imputed`, variants `grouped`/`partitioned`, default `partitioned`); the peptide completion and publication handlers are reused unchanged. The GUI's imputation module keeps its in-process path (dispatching it is KIT-014's work).
4. **Harness.** An operation pair may bind each arm's registry `variant`, its `parameters` and a journey `oracle`; the coordinator hands the variant and the population's head state to the adapters vehicle, whose variant list now comes from the registry entry and whose input state comes from the member record (the CORE populations keep `protein_peptide_filtered`); the verifier judges a journey-oracle member by `kitOracleStateMismatches` (class, rows, every component and support entry, the args summary; products reported). Generator: `tickets$009` once `kit_009` is bound.

## Per-section exactness (acceptance 1)

`test-workflow-partitions.R`: on the fixture (the protein-evidence filter on the precursor table as in production; intensity, replicate and imputation on the rolled-up, sample-filtered table), for every section the partitioned and grouped variants produce the same table, summary and support entries, and the table equals the unchanged kernel run on each protein group alone (the per-group independence check of the declared unit). No section is flagged; none stays on the old loop by exception.

## In-process check on the 300-sample population (before the members)

[in-process-imputation-check.R](in-process-imputation-check.R) → [in-process-imputation-check.json](in-process-imputation-check.json): the KIT-005 `journey_replicate` population (5,515,154 rows, 300 samples) materialised as a snapshot store per variant; proportion 0.5, no exclusion column (the journey's recorded parameters); the state committed through the shared finishing step and fingerprinted under the current canonical form against the `journey_imputation` oracle. Both arms take the same 625 batches (the largest protein group, 55,508 rows, is delivered alone; the admission estimate at the normal cap is 1.61 GiB); the products (`qc_progress`, `qvalue_progress`) are reported, not judged (the GUI's progress context is absent in process).

| Variant | Elapsed | Output rows | Batches | Groups | Exact against the oracle |
|---|---|---|---|---|---|
| grouped (the loop before KIT-009) | 361 s | 11,826,409 | 625 | 2,500 | yes: every component, support entry and the args summary |
| partitioned | 369 s | 11,826,409 | 625 | 2,500 | yes: the same record |

## Tests (working tree, proof library on `R_LIBS`)

partitions 132 (new: partitioner, per-section exactness, declaration, private configuration), operation registry 46, operations 66, execution 64 (10 installed-only skips), adapters (6 installed-only skips), kit harness 252, perf107 consumers 29 / audit 78 / perf-107 21, peptide module contracts 772, impute 51, q-value 57, peptide QC module CI 43, GUI memory lifetimes 252 (3 skips), DIA baseline compatibility 42, saved QC rollup 59 (1 skip); all green.

## Attempt 1 (harness failure before any work; user "go" 2026-10-12)

Stage `kit-009` launched as bound (user decision: one build, two arms). Both members were refused by the adapters vehicle within seconds, before any kernel ran (`controller_failed`; no worker request, no measurement): the vehicle checks the member's input by its `peptide_data` component id, while a journey population is bound by its state descriptor's id (the CORE populations of KIT-006 are bound by the component itself), so `input_entity_id` never matched. The stage directory is kept as `kit-009.attempt1-harness`. Correction: the coordinator records the bound state's assay component (`kitPopulationInputComponent`) as `input_entity_id` and the state descriptor as `input_state_entity_id`; the verifier accepts either; harness test added (255 green). A rerun of both members needs the user's ledger amendment (as KIT-006 attempt 1 did).

## Results (stage `kit-009`, 2026-10-12; tracked summary [kit-009-stage.json](kit-009-stage.json))

| Member | Arm | Status | Elapsed | Slice peak | Anonymous peak | Max gap | Batches | Verified |
|---|---|---|---|---|---|---|---|---|
| control-1 | kit_009, `grouped` (the loop before KIT-009, fixed 1 GB engine) | complete | 438 s | 3.498 GiB | 2.678 GiB | 70 ms | 625 | exact against `journey_imputation`: 2 components, 8 support entries, args summary; products reported |
| candidate-1 | kit_009, `partitioned` (admission 1.43 GiB at the 14.3 GiB cap, engine 1 GB) | complete | 380 s | 3.365 GiB | 2.549 GiB | 70 ms | 625 | exact, same record |

Pair ratios (candidate / control): slice peak 0.962, anonymous 0.952, elapsed 0.866. `not_worse` holds; no reduction claimed (0.962 is above the 0.95 rule and the batches did not change; the difference is the engine's share and the loop's bookkeeping). Zero swap, every gap under 250 ms, slices stopped. Predeclared batch rule held (625 = 625). Decision `accepted` (`kitOperationDecision`); placement: `kit_009` becomes the control (protocol section 2).

### Claim limits

- One member per arm; the arms share the build, so the comparison isolates the partitioner, its admission and the engine bound, not a build-versus-build difference (the control build has no imputation operation; recorded deviation, user decision).
- Exactness covers the state's components, support entries and args summary; the two progress products are reported (the worker has no GUI progress context).
- The protein-evidence, intensity and replicate sections run on the partitioner in the GUI from this build but were measured only on the fixture (per-section exactness); their counted members belong to KIT-014.
- Times are reported, not judged.

## Declared runs (as authorised; recorded above)

The ledger names the pair "control and candidate on the imputation population". The imputation step is a worker operation only from this ticket's build onward: the control build `kit_011` has no `peptide_imputation` operation, so a build-versus-build pair cannot run it. The smallest pair that isolates the mechanism is **one build, two arms**: `control-1` = `kit_009` with variant `grouped` (the per-group loop exactly as it stood, fixed 1 GB engine limit) and `candidate-1` = `kit_009` with variant `partitioned` (the partitioner, admission, cap-bounded engine), both on `journey_replicate` with the journey's parameters, verified against `journey_imputation`, 600 s each. Both arms take the same batches (the cursor's delivery bound), so the expected outcome is `not_worse` within noise with no reduction (as KIT-006's pairs were); the decision stays `kitOperationDecision` (kernel peak and anonymous peak ≤ 1.02 × the control's, exact). The alternative is a rebuilt control (`kit_011` source plus the operation registration), which measures the same thing at the cost of an extra install.

1. **Install** (recorded, not counted): the KIT-009 tree at its commit into `<proof>/core`, bound as `kit_009`; bindings regenerated (`tickets$009`); `run_workflow_kit.sh validate kit-009` must pass.
2. **Stage `kit-009`, two counted operation members, 600 s ceiling each**: `control-1` (kit_009, `peptide_imputation`, variant `grouped`), then `candidate-1` (kit_009, variant `partitioned`), each in its own slice (16 GiB cap, 15.5 GiB early stop, 12 GiB start headroom, 2 GiB host floor, 250 ms gap rule); the verifier runs outside the slice after each member.
3. **Decision (predeclared, `kitOperationDecision`)**: both members eligible and exact against `journey_imputation`; the candidate's kernel peak and anonymous peak ≤ 1.02 × the control's → `accepted` (placement: the kit_009 build becomes the control per protocol section 2). Predeclared batch rule: the candidate's batch count equals the control's (625 at the normal cap). Any rule failing → `rejected` or `invalid`, recorded; one correction rerun of a single member only for a demonstrated harness failure.
