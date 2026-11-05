# KIT-095 plan (rev 1.1 + consistency pass), written 2026-11-06 before any declaration

Bounded item: the transcriptomics binding on the KIT-094 core, the shared DA saving extended to an edgeR QL fit, the
publication operations on KIT-093's renderers, registry tables, record and replay. Exit: the acceptance criteria with
test evidence, the regression checks (test first, targeted, dependents, `R CMD check` at baseline, install, the full-suite
checkpoint) and `outcome.md`. No journeys or kit members (KIT-098).

Order (ticket, planning gate): `kit-001/step-map-trx.json` first and committed; it carries no open flag (the decisions it
applies are the ticket's); then the failing tests; then the declarations.

## State, product and dependency table (KIT095-F04)

| Step | State or product | Reads | Dependencies (parameters) | Replay reproduces |
|---|---|---|---|---|
| design | state `trx_raw_data_s4` | counts, design | – | precondition (never re-run) |
| `trx_expression_filter` | state `trx_expression_filtered` | head | – | the state |
| `trx_qc_finalize` | state `trx_qc_complete` | head | – | the state |
| `trx_tmm_logcpm` | product `trx.tmm_logcpm` | head | – | re-run from the RUVs step's `preceded_by` |
| `trx_ruv_grid` | product `trx.ruv_grid` | head | `logcpm` | re-run from `preceded_by`, `logcpm` the replay's |
| `trx_ruvs` | state `trx_ruvs`; product `trx.ruvs_corrected_logcpm` in its value | head | `logcpm`, `grid` | the state |
| `trx_da_contrast` | manifest `trx.da_results.v2` | head (RUVs) | – | the DA, formula checked |
| `trx_ruviii_export` | product `trx.ruviii` | head | `logcpm`, `grid` (from the RUVs product), `ruviii_m_source` | the export analysis |
| `trx_reset` | store | – | – | the state |

The RUVIII export is product only after DA (decision 4). It is accepted on the head as the DA manifest is
(`snapshot$exports`), and the record lists it under `analyses$trx_ruviii_export` with its parameters and outcomes; replay
re-runs it after the DA through the binding's replay entry.

## Files

- `R/utils_workflow_trx_declarations.R` (new): the binding (fields as lipid's, plus `da_fixed`, `da_design`,
  `da_model`, `da_settings`, `integration_export`); the steps (filter with its partitioned variant and edgeR-free
  template, store steps, TMM/logCPM with its partitioned variant, grid, RUVs, the shared contrast declaration, RUVIII
  export); the DA halves (fit, test) taken from `trxDaContrast`; the record-outcomes hook; the replay entries; the
  in-process wrappers (refuse a non-owner manager); the TRX publication handlers (`trx_ruv_publication`,
  `trx_qc_figures`, `trx_composite_figure`, `trx_da_publication`) and adaptors for the shared ones.
- `R/utils_workflow_omics_declarations.R`: `omicsDeclarationBindings()` gains `transcriptomics`; the contrast
  declaration reads `fixed` and the design from the binding; `acceptWorkflowOmicsDaManifest` reads the model from the
  binding; the integration export calls the binding's adaptor. Each defaults to today's lipid/metab value.
- `R/utils_workflow_saved_entities.R`: `newSavedWorkflowDaFit` accepts a `DGEGLM` (same identity checks, plus
  `engine = "edgeR_QL"`); the `MArrayLM` envelope unchanged; the consumer payload check accepts both shapes.
- `R/utils_workflow_operations.R`: defaults (`fixed`, `input`, `no_shipped_default`) and labels for every `trx_*`
  parameter. No `TRANSCRIPTOMICS` section map (KIT-093's kernels write no `@args`), so every worker step is `section = NULL`.
- `R/utils_workflow_class_declarations.R`: the `trx.da_results.v2` output declaration only.
- `R/utils_workflow_analysis_record.R`: report sources (`ruv_mode`, `ruviii_m_source`; the TRX rows of `ruv_applied`,
  `best_k`, `best_percentage`, `num_neg_ctrl`); the export analysis in the record; outcomes read from an analysis; the
  replay's DA formula check and export replay.
- Tests: `test-workflow-trx-declarations.R` and its helper (journey, manual mode, refusals, technical-replicate RUVIII,
  partitioned equals whole, worker harness for every step and publication, replay, registry shape, non-owner
  refusals); `test-workflow-saved-consumers.R` (edgeR envelope; `MArrayLM` envelope byte-identical);
  `test-workflow-operation-registry.R` (the trx names), `test-workflow-omics-declarations.R` (the added binding name),
  `test-workflow-omics-worker.R` (the trx worker steps).

## Compatibility

Lipid and metab: the frozen dump file and fixture unmodified. `test-workflow-omics-binding-steps.R` is not in this
ticket's files; its characterisation snapshot includes the binding list, the registry by omic and the
`TRANSCRIPTOMICS` type, which grow by the trx entries by construction. If it fails only by those entries, that is shown
by a script that removes them, and the test-file change is reported for a decision (it is outside the write scope).

## Checks

Test first (fails before, passes after) for the declarations test and the envelope test; targeted files one process
each with the kit_093 library; the dependents grep; `R CMD check` (3 WARNINGs, 3 NOTEs); install (kit_095); the full
suite, one process per file, in its own worktree. One bounded correction pass.
