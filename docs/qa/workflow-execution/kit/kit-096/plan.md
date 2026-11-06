# KIT-096 plan (rev 1.1), written 2026-11-06 before any code

Bounded item: the transcriptomics module (import, design, QC) on the managed runner, built by giving the existing lipid
import, design-builder and QC modules one `profile` argument. Exit: the acceptance criteria with test evidence, the
regression checks (test first, targeted, dependents, `R CMD check` at baseline, install) and `outcome.md`. No app launch, no
kit journeys (KIT-098).

## How the profile is threaded (F01)

- `profile = NULL` on each lipid entry point is the lipid profile: every lipid call is made exactly as today (same
  functions, same arguments). The lipid helpers already take their collaborators as arguments (`buildModulePanel`,
  `handleProcessRequest`, `runBuilderObserverShell`, `validateImportState`, ...). A non-NULL profile only supplies the
  transcriptomics collaborators at those seams, plus a few labels with `%||%` lipid fallbacks.
- The profile is a plain list, `trxOmicProfile()`, in `R/mod_transcriptomics.R`: labels and test-id prefix, the id and
  name column labels, class, workflow type `TRANSCRIPTOMICS`, format `trx_counts`, no sanitisation, run order (the counts
  file), the QC sub-tabs, the count-policy check, and the binding's readers and wrappers.
- Metab and proteomics modules are not touched.

## Files and line budgets

| File | Change | Budget |
|---|---|---|
| `R/mod_transcriptomics.R` (new) | `mod_transcriptomics_ui`/`_server` (bus, owner, runner install as `mod_lipidomics.R:180-196`; tabs import, design, QC; the shared stepper); `trxOmicProfile()`; the import preview, count-policy check and processing collaborators; the design's object, validation (KIT-093's `validateTrxDesign` at save, plus at least two runs per group), level order and provenance collaborators | 300 |
| `R/mod_trx_qc_filter.R` (new) | the expression-filter sub-tab (UI, server: managed `omicsManagedCommit` of `trxExpressionFilterStep` + `trx_qc_figures` stage `raw`, or `runTrxReferenceExpressionFilter` in process), its metrics table, the filtering figure, raw-logCPM PCA/RLE; the finalise collaborators (metadata summary, store step `trx_qc_finalize`, text) | 220 |
| `R/mod_lipid_import_ui.R`, `_ui_helpers.R` | `profile` on the UI: header, format choices, id/name labels, test-id prefix, help text; sanitise, IS pattern, vendor-normalised and second-assay controls only for lipid | +40 |
| `R/mod_lipid_import_server.R` | `profile`: the preview loader and the processing request at the existing seams | +15 |
| `R/mod_lipid_import_processing_helpers.R` | the format label of the count table in the detection status; the summary file name as an argument (lipid default) | +3 |
| `R/mod_lipid_design_ui.R`, `_server.R`, `_builder_helpers.R`, `_builder_server.R` | `profile`: help text, no design import for transcriptomics (its import format is lipid's), run order, the level order and excluded runs on the result, the commit's class constructor, validation, owner type, provenance, state name | +40 |
| `R/mod_lipid_qc.R`, `R/mod_lipid_qc_s4.R` | `profile`: title, sub-tabs and their servers, class check, the finalise's runner omic, metrics operation, in-process read, summary | +30 |
| `R/app_ui.R` | the stale "(finding m-2)" comment removed | -1 |
| `DESCRIPTION` | Collate gains the two new files (outside the file list, needed to load them; as KIT-095) | +2 |
| tests: `helper-module-ci-trx-tabs.R`, `test-module-ci-trx-{import,design,qc,app}.R` | MCI-026 import, MCI-027 design, MCI-028 QC, MCI-029 module and stable ids; lane `transcriptomics-unit` | new |
| `tests/testdata/mci/manifest.json`, `generate-fixtures.R`, `fx/transcriptomics/` | four scenarios; the lane; the KIT-040 generator's `trx_screen` corpus (3 per group, 300 genes) and its one-value fixtures | new |
| `helper-module-ci-manifest.R`, `test-module-ci-manifest.R` | `transcriptomics` a valid omic; the scenario sets extended by the new ids only | +10 |
| `helper-workflow-omics-residency.R`, `test-workflow-omics-residency.R`, `test-workflow-omics-application-reads.R` | the walk's classes include `TranscriptomicsAssayData`; the trx journey (import, design, filter, finalise); the KIT-082 probe on the QC tab, managed and in process; positive controls | +150 |

`fixture_packs.json` and `helper-module-ci-fixtures.R` are not changed: the MCI-002 catalog is the vendor-format class
matrix (ten classes per format), and the transcriptomics fixtures are KIT-040's, pinned by their own `fixture.json` and
manifest hashes.

## Data flow

import (file read for one action) -> `data_tbl` reference -> builder (`data_cln` reference) -> commit
`trx_raw_data_s4` (`createTranscriptomicsAssayData(list(counts = ...), design, group_levels, roles =
list(gene_id_column, annotation_id_column, technical_replicate_id = "tech_rep_group"))`, owner bound to `TRANSCRIPTOMICS`,
`recordTrxDesignProvenance`) -> QC filter `trx_expression_filtered` (worker or wrapper) -> finalise `trx_qc_complete`
(store step). The tabs read state metadata or worker products; in process the filter reads the state for one action.

## Interpretations recorded before coding

- Run order follows the counts file (ticket); for the KIT-040 corpus that is `truth$column_order`. KIT-093's class orders
  the kernels' input by the design, so a GUI design runs in counts-file order. Hand-off note for KIT-098: the KIT-040
  reference is in design-file order.
- Group level order: the order of the factor levels the user added (`factors()`), restricted to the assigned groups; a
  group of several factors is ordered by its factors' positions. Recorded in the design state (`args$group_levels`).
- Unassigned runs: the builder drops them from the result; under the profile the commit validates the design against all
  imported count columns except the runs the user removed, so an unassigned run is refused (`trx_design_runs_extra`).
- A group with fewer than two runs is refused at save under `trx_design_group_singleton` (child of `trx_design_invalid`).
- The finalise's per-assay statistics panel is lipid's (one count table here); the data summary is read from the state
  metadata (row counts and design), never the counts.

## Tests first

Every new test file and the extended residency, probe and manifest tests are written and run before any `R/` change; the
failures are saved in `test-first.txt`.
