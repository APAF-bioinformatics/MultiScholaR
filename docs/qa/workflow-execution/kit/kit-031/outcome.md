# KIT-031 outcome: lipidomics declarations (the kit completeness test)

Status 2026-10-31: **completed**. Code 6506af27 (item 1: the declarations, the modules on them, record, replay, overlay, the recorder) + 4d6ac0a4 (bindings, the record on the install) + 5d2b1a39 (the member's record and replay) on `performance-v2`; base 843ff460 (KIT-071 closed). Install `kit_031` = `core/browser-smoke.6oW1oZ` (6506af27) over the kit_071 library. No counted member. Plan: [plan.md](plan.md) (the planning gate, revision 2.1 of the ticket; the user chose the full report overlay on 2026-10-31).

## Kit completeness verdict

The lipidomics route went onto the kit with **declarations, registry entries, record and replay branches and module code only**. No new compute role, enum, partitioner or harness function was needed; the ticket's stop-and-addendum triggers did not fire. Three shared changes were needed and are recorded as divergences (below): the scientific identity reads the assay through the class declaration (a `lists` assay had no single entity id), a declared parameter whose value is NULL is now recorded as null rather than dropped, and the record reads the DA operation and long table from the manifest's declaration rather than the proteomics names. Each is a generalisation of an existing function with the proteomics behaviour unchanged (the KIT-028 replay exact and the three proteomics routes' records identical to their oracles, before and after).

## What changed

| Where | Change |
|---|---|
| `R/utils_workflow_lipid_declarations.R` (new) | twelve `newWorkflowOperation` declarations on `LipidomicsAssayData` (`fields = lipid_data, design_matrix`): `lipid_intensity_filter` (row; `global_statistics` the per-assay quantile), `lipid_duplicate_resolution` (feature group), `lipid_qc_finalize`, `lipid_norm_complete`, `lipid_reset` (store), `lipid_itsd_normalisation` (complete sample; the mean standard factor; `fixed` method ITSD), `lipid_log2` (row), `lipid_between_sample_normalisation`, `lipid_ruv`, `lipid_pearson_pairs`, `lipid_pearson_selection`, `lipid_da_output` (whole matrix; DA product-only, `fixed` eBayes trend/robust). Kernels unchanged, called through their legacy `…Fn` aliases so the recorded `@args` keep their names. Wrappers `runLipidReference*` (operation → products on the checkpoint → `commitWorkflowReferenceState` under the module's state name/config/description), `commitLipidReferenceStoreStep`, `recordLipidDesignProvenance`, `runLipidReferenceDaAnalysis` (the manifest `lipid.da_results.v1`), `replayLipidomicsStep` |
| registry (`utils_workflow_operations.R`) | a second in-process loop with `omic_type = "lipidomics"` and the QC tracking labels; defaults and labels for every parameter; `workflowSectionOperationsByType$LIPIDOMICS` (the `@args` sections `lipidIntensityFilteringFn`, `ITSDNormalization`, `ruvIII_C_Varying`) |
| `workflowOutputDeclarations` | `lipid.da_results.v1` (`isLipidResultManifest`); both manifests name their `operation` and `long_table` |
| the lipid modules | intensity, duplicates, QC finalise, the normalisation action's four kernels, correlation apply/skip, reset, DA run: on the owner each runs its declaration; a test double keeps the direct kernel path (the module tests' seams); the design builder and design import record `design_builder` provenance and the inputs (`lipid_assay_<n>` per import copy, `design_file`) |
| `lipidIntensityFiltering` | the threshold is `quantile(...)[1]`; the `ceiling()` removed (user decision 2026-10-28) |
| record (`utils_workflow_analysis_record.R`) | `analyses$differential_abundance` from the manifest's declaration; the lipid RUV outcomes (k and control counts per assay and in total, the mean percentage) and the ITSD outcome; the input allow-list takes `lipid_assay_<n>` / `metab_assay_<n>`; the replay's lipid branch (`entry$omic_type == "lipidomics"` → `replayLipidomicsStep`) and the lipid DA; the overlay takes a report variable's first record source the record carries, and gains the lipid template's variables |
| the lipid summary | `analysis_record.json` written beside `study_parameters.txt` at Save Workflow Args; the record bundle (JSON, Word, Excel) exported to `Study_report/` after the report renders (both the KIT-027/029 proteomics patterns) |
| `inst/reports/lipidomics/lipidomics_report.rmd` | the parameter chunk overlays the record (`workflowRecordReportParameters`), unchanged without one |
| harness (`tools/profiling/workflow_kit.R`) | `kitOperationIdentityProducts` also matches `lipid_*_science` (the operation products on each state are reported across builds, the states judged) |
| `kit-031/oracle-record.R`, the generator | the in-process recorder for the lipid route; ticket `031l` and the `lipid31_<step>` oracles; `KIT_031_LIFECYCLE=1` |

### Divergences from the ticket's wording, recorded

- **Shared code beyond "declarations only"**, each a generalisation with proteomics unchanged: `workflowComponentScientificIdentity` takes the assay from the class declaration (a `lists` assay identified by the digest of its sorted element ids, the KIT-070 harness rule); `boundDIAProvenanceValues` keeps a declared NULL as null (the ITSD selection when none is made); the record's DA block reads its operation and long table from the manifest's declaration.
- **Store steps** (QC finalise, correlation skip, reset) commit the head's record under the new name with the declaration's provenance; on the owner they no longer re-serialise the module's object (the validated object is the head's).
- **The DA declaration** runs `runLipidsDA` once over every contrast and assay (the kernel's own loop) rather than once per contrast as the proteomics declaration does; the manifest's `results` are keyed by assay. Its parameters carry the shared DA names (`formula`, `qvalue`, `lfc`) so the overlay's `DA` sources serve both routes.
- **The intermediate operation products** (`lipid_*_science` on each state) carry the ITSD kernel's `@args` timestamp in their checkpoint bytes; the harness reports them, as it reports the proteomics operation-identity products.
- **The `ceiling()` removal** changed no retained row on the e2e lane, the module-CI lanes, the declarations fixture or the screen corpus (every DA table and workbook of the kit_031 lipid member identical to kit_071's); recorded in [ceiling-effect.txt](ceiling-effect.txt).
- **Not changed**: the summary's `session_state_<date>.RDS` and `_final_s4.RDS` writers (KIT-071's note); no worker role for a non-proteomics compute; metabolomics (KIT-032).

## Tests

- New `tests/testthat/test-workflow-lipid-declarations.R` (223 expectations; fixture `helper-workflow-lipid-declarations.R`, 24 samples × 2 assays with standards and duplicates): registration with the step map's units, merges and pre-pass statistics; the refusals for a non-owner; the whole journey on a real owner: every commit a component checkpoint with its declaration's provenance, the record valid with nothing incomplete but the external-resources limit, the inputs hashed, the DA manifest the head's and its display list equal to the kernel's direct output (a second run a cache hit), the overlay's lipid variables, the replay on a copy cut to the design state exact (10 states, the DA re-run), a tampered input refused.
- Updated, none deleted: `test-workflow-operation-registry.R` (the two id lists), `test-lipid-15-summary-module-seams.R` (the report observer passes the bus).
- [affected-tests.txt](affected-tests.txt): the lipid module, module-CI and e2e files, the dependents of every shared symbol changed, and the ticket's regression files; all passing after the two updates.

## Results

| Gate | Evidence | Result |
|---|---|---|
| Item 1: KIT-028 proteomics replay | [replay-check-kit031-item1.json](replay-check-kit031-item1.json) | **exact, 16 states, 0 mismatches, 290 s**; the edited, inputs and legacy refusals as before |
| Item 1: LFQ / TMT / DIA_LIMPA records against the kit_065 oracles | [compare-oracles-item1.txt](compare-oracles-item1.txt) | **identical** (11, 11, 14 steps; 0 reproduce and 0 per-state mismatches) |
| Item 1: lipid module-CI | [module-ci-comparison.txt](module-ci-comparison.txt) | design, QC, DA lanes byte-identical; import differs only in test-process scratch; norm only in store entities and export journals; summary only in date lines of `study_parameters.txt` |
| R CMD check | [check-6506af27.log](check-6506af27.log) | baseline: 3 WARNINGs, 3 NOTEs; the code note identical to KIT-071's |
| The recorder against the GUI (before the install) | [recorder-vs-kit071-member.txt](recorder-vs-kit071-member.txt) | the tree's recorder equal to the kit_071 member's store state by state (9 journey states and the 3 intermediate normalisation states), after mirroring two GUI facts: whole-number tab inputs arrive as integers; the internal-standard table preselects its candidates |
| The record on the install | `oracles-lipid31_lipid_screen.rds` (bound `lipid31_<step>`) | 10 states; the analysis record valid, 0 legacy steps; 48 s |
| **Acceptance: the screen member exact through the report** | [kit031-journey-031l-stage.json](kit031-journey-031l-stage.json); `<proof>/kit/kit-journey-031l.log` | **accepted, exact**: 0 oracle mismatches, 0 canonical-form differences, 10/10 steps; 239 s, 2.54 GiB peak |
| **Acceptance: the member's record valid and its replay exact** | [member-record-replay-031l.json](member-record-replay-031l.json) ([member-record-replay.R](member-record-replay.R)) | the record valid, 10 steps, 0 legacy, nothing incomplete but the external-resources limit; the replay from the design state **exact** (10 states, 0 mismatches, the DA re-run); a tampered input refused |
| The two KIT-068 lifecycle members on kit_031 | [kit031-lifecycle-lipid-stage.json](kit031-lifecycle-lipid-stage.json), [kit031-lifecycle-metab-stage.json](kit031-lifecycle-metab-stage.json); [verify-member-kit031-lipid.txt](verify-member-kit031-lipid.txt), [verify-member-kit031-metab.txt](verify-member-kit031-metab.txt) | both `lifecycle_complete_observed` at the first attempt (lipid 252 s / 2.52 GiB, metab 231 s / 2.43 GiB); the lipid store's ten entries all declared (`import`, `design_builder`, the eight lipid operations); the record bundle in `Study_report/`; metabolomics still `legacy` (KIT-032) |
| Member projects against the kit_071 members' | [member-project-comparison-kit071-vs-kit031.txt](member-project-comparison-kit071-vs-kit031.txt) | both omics: every DA long table, publication workbook, NumSigDeMolecules table and `study_parameters.txt` identical; the session RDS by the ITSD timestamp only; the lipid store gains the operation products and the record bundle |
| Closure: KIT-028 replay | [replay-check-kit031-final.json](replay-check-kit031-final.json) | **exact, 16 states, 0 mismatches, 266 s**; the inputs and legacy refusals hold; within the ceiling |

## Backward compatibility and rollback

- Public API, S4 classes, checkpoint and saved formats unchanged; file names unchanged; the proteomics routes' records identical to their oracles; the fastduckdb pin untouched. A lipid store written before this ticket reads as before (its steps `legacy`); an older RDS session adopts as in KIT-071.
- Rollback: `git revert` of 5d2b1a39, 4d6ac0a4 and 6506af27; the stores already written keep their states (the provenance entries name operations an older build does not register and the record then lists them unregistered).

## Handoff

KIT-032 mirrors this for metabolomics: the declarations on `MetaboliteAssayData` (`metabolite_data`), `metab.da_results.v1`, the `METABOLOMICS` section map, the metab `ceiling()` (`R/func_metab_s4_qc_methods.R:102-105`), the recorder and member `032m`, and the stage closure (the closure full suite). The wrappers, the replay branch and the overlay's first-source rule are written to take a second omic by a second set of declarations.
