# KIT-096 outcome: the transcriptomics import, design and QC tabs (2026-11-06)

Ticket revision 1.1. Plan: [plan.md](plan.md), written before any code. In-process and module CI only: no app launch, no
kit journeys (KIT-098).

## What was done

- **Reuse by parameterisation (F01).** The lipid import, design-builder and QC modules take one `profile` argument.
  `NULL`, the default, is lipidomics: every lipid call is made as before (same functions, same arguments; the profile is
  passed to a collaborator only when it is not `NULL`). A profile supplies the omic's collaborators at the seams these
  modules already take (`buildModulePanel`, `handleProcessRequest`, `validateImportState`, `writeImportArtifacts`,
  `applyResultToWorkflow`, `runBuilderObserverShell`, the finalise's `get/validate/save/update/report` functions), plus
  labels with lipid fallbacks. Nothing is copied. Metab and proteomics modules are not edited.
- **The profile** (`trxOmicProfile()`, `R/mod_transcriptomics.R`): labels and the `trx` test-id prefix, class
  `TranscriptomicsAssayData`, workflow type `TRANSCRIPTOMICS`, format `trx_counts`, the id and name column labels, the
  preview loader and processing, run order (`identity`: the counts file), the level order, the design check, object,
  provenance and state name, the QC title, sub-tabs, sub-module servers, metadata and state readers, summary and the
  finalise's collaborators.
- **Wiring (F13).** `mod_transcriptomics_ui`/`_server`: the bus, the session owner and the runner install of
  `mod_lipidomics.R`, tabs import, design and QC, the shared stepper (Import, Design, QC), QC enabled when the design
  completes. `app_server.R` dispatches by name and is not edited; the stale "(finding m-2)" comment above the
  transcriptomics tile in `app_ui.R` is removed (the integration tile's comment stays).
- **Import (F05, F06, F14).** One format, the notebook's count table read with vroom (as the notebook and KIT-095's
  journey read it). The user picks the gene id and name columns (`GeneID`, `GeneName` proposed); every other column is a
  sample, in file order. KIT-093's count policy runs at import (`checkTrxImportCounts`): a refusal keeps KIT-093's class
  and names it in the message the tab shows, e.g. `NA or NaN count (...) (trx_counts_na)`; there is no rounding. No
  sanitisation (the control is not offered). No vendor-normalised, internal-standard or second-assay controls. The help
  text states the whole-number policy. The owner is bound to `TRANSCRIPTOMICS` at import; the summary file is
  `transcriptomics_import_summary.tsv`.
- **Design (F04, F06, F07, F10, F15).** The lipid builder and its commit under the profile:
  - a fresh design keeps the counts file's run order (no `mixedsort`) and names;
  - at save the builder's result carries the level order the user entered (`trxDesignGroupLevels`: the factor levels in
    the order added, restricted to assigned groups) and the runs the user removed;
  - the commit runs KIT-093's `validateTrxDesign` against every imported count column except the removed runs (so an
    unassigned run is refused as `trx_design_runs_extra`), by name only, plus at least two runs per group
    (`trx_design_group_singleton`, a child of `trx_design_invalid`); each refusal is shown with its class;
  - the state `trx_raw_data_s4` is `createTranscriptomicsAssayData(list(counts = ...), design, group_levels, roles =
    list(gene_id_column, annotation_id_column, technical_replicate_id = "tech_rep_group"))`, the owner bound to
    `TRANSCRIPTOMICS`, provenance through `recordTrxDesignProvenance`; the completion text is unchanged
    ("Design matrix and contrasts saved successfully!");
  - the "Import Existing Design" button is not offered under the profile (its directory format is the lipid import's);
    the help text says the design is matched by sample name only.
- **QC (F08, F09).** The QC tab under the profile has two sub-tabs:
  - **Expression Filter** (`R/mod_trx_qc_filter.R`, new): KIT-095's `trx_expression_filter` with no library-size
    argument; managed it is `omicsManagedCommit(trxExpressionFilterStep())` then `omicsManagedProduct("trx_qc_figures",
    list(stage = "raw"))`, in process `runTrxReferenceExpressionFilter` and the state read for one action. It shows
    genes before, after filterByExpr and after the zero-variance step, library size and genes detected per run (the
    filter's `filtering_metrics` product), the filtering figure and the raw-logCPM PCA and RLE. Revert reuses the lipid
    intensity sub-tab's revert observer.
  - **Finalize QC**: the lipid finalise under the profile: the store step `trx_qc_finalize` -> `trx_qc_complete`, the
    worker's `trx_filtering_metrics` on a managed session (the profile's runner and operation), a data summary from the
    state's metadata (row count and design; no count read), no per-assay table.
- **Stable ids (F11)**: module `transcriptomics_workflow`; tab input `transcriptomics_tabs` (`import`, `design`, `qc`);
  sub-modules `import`, `design`, `qc` (QC sub-tabs `filter`, `s4_finalize`; buttons `qc-filter-apply_filter`,
  `qc-s4_finalize-finalize_qc`); the builder's shared ids and completion text; test ids `trx-workflow-tabs`,
  `trx-tab-import`, `trx-tab-design`, `trx-tab-qc`, and `trx-import-format`, `trx-import-assay1-file`,
  `trx-import-id-col`, `trx-import-annotation-col`, `trx-import-process`, `trx-qc-apply-filter`.

EVIDENCE_PLACEHOLDER
