# KIT-003: list-of-tables slots and result classes in the registry

Status 2026-10-09: **complete** at `e723f30d`. No counted member. Authority: `tickets/kit/00-KIT-PLAN.md`; protocol v1.1.

## What changed

- `R/utils_workflow_class_declarations.R`: the declaration vocabulary gains `kind` (`assay` or `result`) and the slot kinds `lists` (named lists whose data.frame elements become one table component each, other elements one value each), `nested` (an S4 slot stored as a reference to its own component checkpoint) and `opaque` (a slot stored as one serialized value). Entries added: `LipidomicsAssayData` and `MetaboliteAssayData` (assay is the list slot, wide), `LipidomicsDifferentialAbundanceResults` and `MetabolomicsDifferentialAbundanceResults` (nested `theObject`; opaque `fit.eb` and the four plot lists; list `contrasts_results_table`; three tables), `da_results_for_enrichment` (tables `contrasts`, `design_matrix`; list `da_data`) and `EnrichmentResults` (table `contrasts`; list `enrichment_data`; opaque plots, plotly and summaries). `workflowOutputDeclarations()` declares the existing `dia.da_results.v1` manifest as the proteomics DA output (audit B11): no new class. `workflowAssayDimensions()` gives rows, columns and samples for table or list assays and NA for result classes.
- `R/utils_workflow_dia_state.R`: `checkpointDIAExecutionResult` stores list elements under `components$<slot>$<element>`, nested objects as a `dia.nested_checkpoint.v1` reference to their own checkpoint (never a second copy: the nested tables deduplicate to the same entity ids as the assay saved on its own), opaque slots as one value, and keeps the unchanged-previous rule for all kinds. `readComponentMetadata` accepts a list of references as a required component. A new `rebuildDIACheckpointObject(context, spec, fields, columns)` is the single rebuild path: `readDIAExecutionComponents` (projected; `fields` may name `"<slot>/<element>"`) and `materializeDIAStateRecord` (full, then support entries) both use it. An unrequested nested object is rebuilt with zero-row assay prototypes but with its design table, because the nested class's own validity requires sample ids.
- `readDIAStateMetadata` stays the DIA-only wrapper (KIT-002); non-DIA classes read through `readComponentMetadata`.

## Tests (`tests/testthat/test-workflow-class-declarations.R`, 161 assertions)

- Declarations and slot kinds for all six classes; the manifest declaration; an invalid declaration (a slot in two kinds) refuses.
- Lipid and metabolite objects (package fixtures `module_ci_lipid_norm_object`, `module_ci_metab_norm_object`): one component per element, exact full round trip, projected read of one element with zero-row prototypes for the others.
- DA result objects for both omics (nested fixture object, a limma `MArrayLM` fit, contrast tables, ggplot lists): tables exact; nested object equal and stored once (its table entity ids equal those of the separately saved assay state); opaque slots equal by serialize hash, plots by their built layer data (their environments are re-created on read, as the ticket allows); an empty opaque slot stays the class prototype; projected read hydrates one table and leaves the fit at its prototype.
- `da_results_for_enrichment` (a NULL list element survives) and `EnrichmentResults` (nested list elements as values, plots as opaque) round-trip.
- Regressions: checkpoint state 238, saved QC rollup 76, bulk relational 45 (the exact P14 cases, in-process), native session 425 (the one failure is the known source-load process-loss artefact, identical at baseline).

## Scientific contract

None; storage shape only. Existing DIA checkpoints are unchanged (the two DIA classes have no list, nested or opaque slots, so their records are written exactly as in KIT-002).

## Limit recorded for the consumers

The entity store refuses payloads whose environments hold active or lazy bindings (`checkSavedWorkflowPayload`). A plot object built inside a function whose arguments are still promises is such a payload, so an opaque slot holding it is refused with `saved_entity_unsupported`. The tests build plots in a helper that forces its arguments. KIT-016, KIT-030 and KIT-031 should keep producing plots as bounded products (plot data, as the DIA QC products do) rather than relying on opaque plot objects; plotly and htmlwidget objects are likely to carry such environments.

## Hand-off

KIT-004 consumes the DA output declaration and the registry; KIT-030 consumes the lipidomics and metabolomics declarations.
