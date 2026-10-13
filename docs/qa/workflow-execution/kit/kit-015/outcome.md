# KIT-015: protein QC declarations

Status 2026-10-14: **accepted and complete.** Stage `kit-015` decision `accepted` (four single members, one launch, no correction); `kit_015` (`core/browser-smoke.hXymGJ`, source `3ad41274`, after the pre-execution audit repairs; the earlier install `Zi0y8A` from `31603afc` is recorded, never run) is the control from KIT-016 onward (protocol section 2). Pre-execution audit: [KIT-015-PREEXEC-AUDIT-20261014.md](KIT-015-PREEXEC-AUDIT-20261014.md) (6 findings, none blocking) and its [remediation](KIT-015-PREEXEC-REMEDIATION-20261014.md); the IQ step rerun in process after the repairs: exact, 198 s ([in-process-protein-qc-check-iq.json](in-process-protein-qc-check-iq.json)). Authority: `tickets/kit/00-KIT-PLAN.md` (run ledger KIT-015: up to five single members, exactness and completion, 10 min each); protocol v1.3 (single members: oracle hash, peak recorded, no ratio; the memory decision is KIT-017's). Control: `kit_014` (`core/browser-smoke.Dl22nP`). Predecessor cases: PERF-114, PERF-115, PERF-117, PERF-118; audit findings A2, A4, A5, A6, A11.

## Unit gate (simulation setup 1)

`kit-001/step-map.json` (bound by hash): `iq` `feature_group`, `accession` `row`, `protein_intensity` `whole_matrix`, `duplicates` `feature_group`, `protein_replicate` `row`, each with its code evidence; `open_flags` is 0 and none carries a flag. The protein table at 300 samples is 2,500 proteins wide, so every step is placed whole (one hydrated object, the kernel once); the declarations carry the step map's units and the merge rules a partitioned placement would need, none of which is exercised here (claim limit).

## What was built

1. **`R/utils_workflow_prot_declarations.R`** (new): compute-only entry points with the protein contract (`checkpoint`, `reference`, `state_name`, `description`, `config`, `section`, `protein_count`, `execution`, display `text`, the step's QC record `qc_name`/`qc`, the product `field`):
   - `computeDIAReferenceIQRollup(manager, reference, parameters(peptide_qc_dir, protein_qc_dir), worker_cap, variant = "whole", seams)`: the module's reference path as a function (IQ input streamed by `writeDIAIQInput`, `iq::process_long_format` with `normalization = "none"`, the protein object built as the module builds it, checkpointed through the class declaration). Declaration `diaProteinIQRollupOperation()` (unit `feature_group`, global statistic `surviving_sample_set`).
   - `computeDIAReferenceAccessionCleanup(..., parameters(delimiter, aggregation_method, fasta_table), seqinr, fasta_metadata)`: the unchanged `chooseBestProteinAccession` through the module's alias with the FASTA table from the caller (the GUI's bus), the bound `fasta_table` path or the project's source directory, else the module's skip. Declaration `diaProteinAccessionCleanupOperation()` (unit `row`).
   - `computeDIAReferenceProteinIntensityFilter(..., parameters(strict_mode, min_reps_per_group, min_groups, proteins_intensity_cutoff_percentile))`: the cutoffs through `updateMissingValueParameters`/`updateConfigParameter` on a private configuration (KIT-009), the unchanged `removeRowsWithMissingValuesPercent` through its alias. Declaration `diaProteinIntensityFilterOperation()` (`whole_matrix`, audit A2).
   - `resolveProteinDuplicateRows(table, protein_id_column, aggregation_fn)`: the module's inline dplyr extracted verbatim (audit A5), and `computeDIAReferenceProteinDuplicateRemoval(..., parameters(aggregation_method))` on it. Declaration `diaProteinDuplicateRemovalOperation()` (`feature_group`).
   - `finishDIAReferenceProteinStep(workflow, value, publish)`: the shared protein finishing step (commit, QC record merged into the bus, the output product through `finishDIAProteinAction`); `completeDIAProteinAction` (KIT-010) delegates to it for any value that carries a product field and keeps the replicate completion otherwise, now with a `publish` argument for the GUI's in-process replicate path. `readManagedProteinResult` recovers the product for a module's managed step.
2. **Registry**: `protein_iq_rollup` (input `peptide`, compute/completion/publication `protein`, label `9_protein_s4_created`), `protein_accession_cleanup` (`10_protein_accession_cleaned`), `protein_intensity_filter` (`11_protein_intensity_filtered`), `protein_duplicate_removal` (`12_duplicates_removed`), each with its validator and handler, variant `whole`; `protein_replicate_filter` (KIT-010) keeps `whole`/`rows` with `rows` the default and `whole` the adapters-only control arm (hand-off decision).
3. **Modules**: on the on-disk default each protein QC module runs its declaration in process (compute, then the shared finishing step; the IQ module keeps its science request directories and dropped-sample notification; the replicate module runs KIT-010's declaration and records its QC entry); each server installs the managed step for its operation (`read_result = readManagedProteinResult`, the plots carry the acknowledgement identity); every protein module takes the experiment paths third (`mod_prot_qc_protein.R`: `needsExperimentPaths` for all five). Limpa's path is untouched.
4. **Harness**: `tickets$015` (kind `single`, four members, the FASTA table bound as the journey project's `aa_seq_tbl_final.RDS`, read-only).

## Grep evidence (acceptance 1, claim limit)

In each protein QC module the kernel call sits after the on-disk early return, so the on-disk default never calls a kernel: IQ `R/mod_prot_qc_protein_rollup.R:243` (reference branch returns at the declaration) vs `:336` (in-memory path); accession `cleanup.R:127` vs `:162`; intensity `intensity.R:130` vs `:181`; duplicates `dedup.R:78` vs `:131`; replicate `replicate.R:92` vs `:131`. The in-memory owner path (unreachable on the DIA default since KIT-011) is the fallback the module tests cover.

## Tests (working tree, proof library on `R_LIBS`)

`test-workflow-prot-declarations.R` 100 (new: the four declarations validate with the step map's units; the dedup kernel equals the module's inline code on a table with duplicates and the module's in-memory step; the duplicate removal and intensity declarations reproduce the module's kernel path on the saved fixture; the accession declaration skips without FASTA data as the module does; the IQ declaration reproduces the module's reference path on the imputed fixture and IQ per protein group equals IQ on the whole table for the first five groups (partition check); the worker entry point computes the IQ rollup from a frozen request and the controller's completion finishes it through the protein finishing step; the protein result reader merges QC entries), operation registry 79, perf107 audit 193 (new: the five protein modules under a managed application and in process; the replicate module's two observers recorded), protein module contracts 651, cleanup 53, rollup 54 (1 skip), replicate 42, orchestrator 19, module CI 68, partitions 293, execution 64 (10 installed-only skips), harness 292, Limpa e2e (1 skip, unchanged); all green. Test drift of this ticket: the protein suites' `renderPlot`/plot-binder mocks widened to the acknowledgement identity; the orchestrator mocks take the experiment paths for every module.

Pre-existing failures, identical at the committed baseline (`619fa1ad`, checked in a clean worktree) and outside this ticket: `test-prot-02d-qc-protein-intensity-shared.R` (5 failures, 2 errors) and `test-prot-02l-qc-protein-dedup-shared.R` (1 error): the modules' in-memory path reads `@protein_id_column` on a fake state class without that slot (commit `6662b785`).

## In-process check on the 300-sample populations (before the members)

[in-process-protein-qc-check.R](in-process-protein-qc-check.R) → [in-process-protein-qc-check.json](in-process-protein-qc-check.json): each step on the KIT-005 control journey's state that fed it, materialised as a snapshot store by the harness, through `runProjected` on the registry's declaration with its production variant, committed through the protein finishing step and fingerprinted under the current canonical form against the journey's next state (components, support entries, args summary; products reported). The accession cleanup read the journey project's `aa_seq_tbl_final.RDS` (cleanup applied, as the oracle records). Not counted.

| Step | Operation | Input | Oracle | Elapsed | Proteins | Exact |
|---|---|---|---|---|---|---|
| IQ rollup | `protein_iq_rollup` | `journey_imputation` (11,826,409 rows) | `journey_iq` | 210 s | 2,500 | yes |
| accession cleanup | `protein_accession_cleanup` | `journey_iq` (2,500 × 297) | `journey_accession` | 19 s | 2,500 | yes |
| protein intensity filter | `protein_intensity_filter` | `journey_accession` | `journey_protein_intensity` | 21 s | 2,500 | yes |
| duplicate removal | `protein_duplicate_removal` | `journey_protein_intensity` | `journey_duplicates` | 26 s | 2,500 | yes |

Every step exact on the first run (the parameters were bound from the oracle states' configurations: delimiter `;`, aggregation `mean`; flexible intensity filter with 2 replicates, 2 groups, percentile 1; KIT-014 lesson).

## Results (stage `kit-015`, 2026-10-14, user "go"; tracked summary [kit-015-stage.json](kit-015-stage.json))

| Member | Operation (variant `whole`) | Population → oracle | Status | Elapsed | Slice peak | Anonymous peak | Max gap | Verified |
|---|---|---|---|---|---|---|---|---|
| single-1 | `protein_iq_rollup` | `journey_imputation` (11,826,409 rows) → `journey_iq` | complete | 241 s | 4.759 GiB | 2.569 GiB | 130 ms | exact |
| single-2 | `protein_accession_cleanup` (FASTA table bound, hash checked) | `journey_iq` → `journey_accession` | complete | 31 s | 1.205 GiB | 1.186 GiB | 60 ms | exact |
| single-3 | `protein_intensity_filter` | `journey_accession` → `journey_protein_intensity` | complete | 32 s | 1.193 GiB | 1.177 GiB | 60 ms | exact |
| single-4 | `protein_duplicate_removal` | `journey_protein_intensity` → `journey_duplicates` | complete | 39 s | 1.366 GiB | 1.340 GiB | 60 ms | exact |

Every member eligible (zero swap, every gap under 250 ms, slices stopped, under the 600 s ceiling) and exact against its journey oracle (class, rows, every component and support entry, args summary; products reported). Decision `accepted` (`kitSingleDecision`, rule `exactness_and_completion`); no ratio and no memory claim (protocol: single members; the memory decision is KIT-017's). Placement: `kit_015` becomes the control (latest accepted build). The IQ member's slice peak (4.76 GiB against 2.57 GiB anonymous) is file cache: IQ's 11.8-million-row input file written and read back inside the slice.

### Claim limits

- Exactness and completion only: peaks recorded, not judged, on the operation-member floor; no comparison with any control.
- Every step placed whole; the step map's units and the declarations' merge rules are carried, not exercised (2,500 proteins).
- The accession applied path is proven on the real population and the member (cleanup applied with the journey's FASTA table), not on the fixture.
- The GUI's managed dispatch for the five protein modules is covered by mocked module tests; the closure journey pair (KIT-017) is the first real-browser observation.
- The modules' in-memory owner path still calls the kernels (grep evidence above), unreachable on the DIA default.

## Declared runs (as authorised; ledger row KIT-015: up to five single members, 10 min each)

Four single members on the one `kit_015` build (the replicate filter was measured by KIT-010), each on the KIT-005 population the journey fed it with the parameters its oracle state records, verified against the journey's next state; exactness and completion only (protocol: peak recorded, no ratio, no control member; the memory decision is KIT-017's).

1. **Install** (recorded, not counted): the KIT-015 tree at its commit into `<proof>/core`, bound as `kit_015` (`hXymGJ` from `3ad41274`; `Zi0y8A` superseded by the audit repairs); bindings regenerated (`tickets$015`, kind `single`; the accession member binds `fasta_table` to the control journey project's `aa_seq_tbl_final.RDS`, read-only); `run_workflow_kit.sh validate kit-015` must pass.
2. **Stage `kit-015`, four counted single members, 600 s ceiling each**, in order: `single-1` (`protein_iq_rollup` on `journey_imputation`, oracle `journey_iq`; IQ's working files under the worker's request directory, published by the controller), `single-2` (`protein_accession_cleanup` on `journey_iq`, oracle `journey_accession`), `single-3` (`protein_intensity_filter` on `journey_accession`, oracle `journey_protein_intensity`), `single-4` (`protein_duplicate_removal` on `journey_protein_intensity`, oracle `journey_duplicates`); each in its own slice (16 GiB cap, 15.5 GiB early stop, 12 GiB start headroom, 2 GiB host floor, 250 ms gap rule); the verifier runs outside the slice after each member.
3. **Decision (predeclared, `kitSingleDecision`)**: every member eligible and exact against its oracle → `accepted` (placement: `kit_015` becomes the control, protocol section 2); an oracle mismatch → `rejected`; any other rule → `invalid`; recorded either way; one correction rerun of a single member only for a demonstrated harness failure. Expected: all four exact (the in-process check above); the IQ member is the closest to the ceiling (210 s in process) and the first member to run a peptide-input operation with the protein compute role on the worker (the launch chain was exercised in process: registry, handler, worker entry point from a frozen request, completion, the result reader; the publication phase for `protein_s4_created` copies IQ's files from the request directory).
