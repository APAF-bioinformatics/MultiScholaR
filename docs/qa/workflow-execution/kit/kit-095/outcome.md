# KIT-095 outcome: the transcriptomics binding (2026-11-06)

Ticket revision 1.1 plus the consistency pass. In-process only (no journeys, no kit members: KIT-098). Plan:
[plan.md](plan.md). Step map: `kit-001/step-map-trx.json`, committed before any declaration (`41bf4904`), `open_flags: 0`.

## What was done

- **Step map first.** Every step's unit, its global statistics and line evidence from KIT-093's kernels and the installed
  edgeR 4.10.1, RUVSeq 1.46.0, EDASeq 2.46.0 and ruv 0.9.7.1 ([environment.txt](environment.txt)). Its
  `kernel_observations` record the determinism line (F13), the per-cell exactness of `filterByExpr` and `cpm` given the
  library sizes, TMM as a whole read, the edgeR-free template (F05), and two places where the declarations follow the
  kernels: the grid's k cap is the kernel's rule (a fixed setting, not a parameter), and the DE kernel is taken as its fit
  and test halves. No row carries a flag: every decision it applies is the ticket's. Its bindings hash is KIT-098's.
- **The binding** (`R/utils_workflow_trx_declarations.R`, new): prefix `trx`, class `TranscriptomicsAssayData`,
  `count_data`, input prefix `trx_counts_`, steps `trx_expression_filter` (row; `whole` and an exact `partitioned`),
  `trx_qc_finalize`, `trx_tmm_logcpm` (product; `whole` and an exact `partitioned`), `trx_ruv_grid` (product),
  `trx_ruvs` (state `4_RUVs`, W in the design; the display matrix a saved table of its product), `trx_da_contrast` (the
  shared declaration), `trx_ruviii_export` (product after DA), `trx_reset`. Worker steps are every compute step;
  default variant `whole` for all (KIT-099 decides). Dependency parameters `logcpm`, `grid`; native step the contrast;
  filtering step the filter. Ten publication operations (below). Record outcomes and replay entries are the binding's.
- **States, products and data flow (F04).** As the plan's table. The RUVs product keeps the logCPM and grid it used as
  references (an entity id is run-specific; a reference is compared by content), so the export after DA reads them from
  the head. The RUVIII export is accepted on the head as the DA manifest is (`snapshot$exports`), the record lists it as
  the analysis `trx_ruviii_export`, and replay re-runs it after the DA.
- **RUV selection (F11).** `ruv_mode` `automatic` (the notebook) or `manual`; `skip` is refused
  (`trx_ruv_mode_unsupported`); manual runs the kernel at `manual_percentage` and uses `manual_k`, refused above the cap
  (`trx_ruv_k_above_cap`). Outcomes: mode, grid settings, the realised cap, k, percentage, control count,
  `controls = "inferred"` with its rule, separation and composite scores, failed members, `ruv_applied`.
- **DA (F01, decisions 2 and 3).** `newSavedWorkflowDaFit` accepts an edgeR `DGEGLM` with the same identity checks and
  records `engine = "edgeR_QL"`; the `MArrayLM` envelope is the five-field list byte for byte (test). The consumer
  payload check accepts both shapes and refuses an engine that is not the fit's. The shared contrast declaration reads
  `da_fixed` and `da_design` from the binding, `acceptWorkflowOmicsDaManifest` reads `da_model`; lipid and metab keep the
  limma values (the frozen dump is identical). The formula is realised from the head's design
  (`~ 0 + group + W_1 + ... + W_k`) and recorded; replay refuses a different one (`workflow_replay_formula_mismatch`). One
  declaration per contrast on one saved fit; `lfc` a display threshold, no treat. Manifest `trx.da_results.v2`.
- **Publications (F02, F06-F08, KIT096-F09, KIT097-F06/F08).** Shared handlers: `trx_filtering_metrics`,
  `trx_da_tab_results`, `trx_final_pca`, `trx_volcano_page` (Glimma on the contrast's DGELRT; a failure is the page's
  recorded message), `trx_session_export`, `trx_integration_export` (the shared handler with a binding adaptor: the RUVIII
  matrix, samples by genes with `Run`, as the notebook writes it; lipid and metab keep the S4 RDS). TRX handlers on
  KIT-093's renderers: `trx_ruv_publication` (CCA panels, the workbook, `ruv_normalised_log_values.tsv`), `trx_qc_figures`
  (per stage `raw`, `tmm`, `ruvs`, `ruviii`: PCA and RLE statistics saved, figures staged), `trx_composite_figure`
  (`normalization_comparison_grid`, `ruviii_before_and_after`), `trx_da_publication` (per contrast the DE CSV, the volcano
  and bar plots; `top_de_genes.tsv`; `ruv_optimization_summary.txt` with the project's dataset line). Roles
  `gene_qc_dir`, `da_output_dir` (leaf `da_genes`), `publication_graphs_dir`. The integration export writes one path per
  action; both notebook paths are two actions (tested, same hash).
- **Registry.** Defaults `no_shipped_default` for the notebook's settings and `input` for the dependencies, the realised
  formula and the contrast; labels for every parameter; no `TRANSCRIPTOMICS` section map (KIT-093's kernels write no
  `@args`), so every worker step and publication declares `section = NULL`. Output declaration `trx.da_results.v2`
  (`operation = "trx_da_contrast"`, `long_table = "da_genes_long"`, `identifies = "isTrxResultManifest"`).
- **Report sources.** `ruv_applied`, `best_k`, `best_percentage`, `num_neg_ctrl` gain the `trx_ruvs` source; `ruv_mode`
  and `ruviii_m_source` are added; `da_model`, `s4_de_q_val_thresh` and `s4_formula_string` already read the DA.

## Evidence

- **Test first** ([test-first.txt](test-first.txt)): before any `R/` change the six declarations tests errored, the edgeR
  envelope test errored, and the registry and binding-name tests failed.
- **Journey** (`test-workflow-trx-declarations.R`, trx_edge_12): design, filter, QC finalise, TMM/logCPM, grid, RUVs, DA,
  RUVIII export. Every product named in the KIT-040 product map that the declarations produce has the reference's
  sha256: the filter's library sizes and statistics, the filtered counts, the samples table, logCPM, the ranking, the
  cap, every grid member and its curves and controls, the summary, the winner, W, the RUVs counts, the corrected logCPM,
  the design, the formula, the contrast, the contrast matrix, the fit, the DGELRT and its table, the top table, newY, M,
  the control mask and k (KIT-093's kernels equal the same references, so the declarations equal the kernels). The
  filter's counts equal the kernel's masks. The record is valid; its outcomes, its DA model and the export analysis are
  asserted, and the report overlay reads `ruv_mode`, `ruviii_m_source`, `da_model`, `s4_formula_string`, `best_k`,
  `ruv_applied`. **Replay** on a copy cut to the design state: every state identical (kit harness), the replayed DA's
  parameters and top table and the replayed export's newY identical; a tampered formula is refused by class.
- **Manual mode, refusals:** skip refused, manual k 4 above the cap 3 refused, the head unchanged; manual k 2 at 10%
  recorded (`ruv_mode = "manual"`, k 2, percentage 10, controls inferred, applied).
- **Technical-replicate RUVIII:** a distinct product (entity id) whose M is `replicate.matrix(tech_rep_group)` and whose
  newY is the direct `ruv::RUVIII` call; the record's export carries `ruviii_m_source = "technical_replicate"`, its
  parameter and its product id; without technical replicates it is refused (`trx_ruviii_replicates_invalid`).
- **Partitioned equals whole** (trx_edge_16; TMM factors not all 1, effective library sizes not the column sums): the
  filter's value and state, the logCPM table and the samples table identical; the logCPM equals the KIT-040 hash; the
  one `DGEList` call of the partitioned filter is the pre-pass on every gene (the template path makes none).
- **Worker harness** (`helper-workflow-omics-worker.R`): the six compute steps and the ten publication operations
  complete and are accepted; every worker product has the KIT-040 hash; the record is valid; the store steps refuse
  dispatch. Files: `de_results_Sera_vs_RPMI.csv`, `top_de_genes.tsv`, `ruv_normalised_log_values.tsv` and
  `transcriptomics_normalized_data.tsv` (both paths) equal KIT-040's stable hashes; the summary text equals KIT-040's
  `summary_text` with the timestamp excluded and the dataset line the declared difference; the workbook's grid sheet
  equals the grid summary; the PCA scores and RLE whiskers of the four stages equal the reference's; the Glimma page is
  saved; the final PCA's matrix equals the RUVs product's.
- **Registry shape and owners:** the 18 names, worker/in-process roles, `section = NULL`, defaults and labels for every
  parameter, step labels, dependency fields, the output declaration, the route; nine wrappers refuse a manager that is
  not the owner (`saved_entity_scope`).
- **Targeted and dependents** ([tests.txt](tests.txt), [dependents.txt](dependents.txt)): every file passes with 0 fail
  and 0 error except `test-workflow-omics-binding-steps.R` (below). The trx files show SKIP 0.
- **R CMD check** ([check.log](check.log), `--no-tests --no-manual --no-build-vignettes`, HEAD `815db998`, kit_093 library):
  3 WARNINGs and 3 NOTEs, the same sections as KIT-094's baseline; the dependencies WARNING and the code-problems NOTE
  bodies are identical (no new global or undeclared import).
- **Install** ([installed-test.txt](installed-test.txt)): `kit_095` = `$P/core/browser-smoke.Q1VVzp` (source `815db998`),
  installed through the kit runner from kit_093's library (RUVSeq, EDASeq and the rest carried). One run of
  `test-workflow-trx-declarations.R` against the installed build: 308 pass, 0 fail, 0 error, 0 skip.
- **Full suite checkpoint** ([suite.txt](suite.txt), `$P/kit/scripts-095/suite.sh`, the scratch worktree at `815db998`,
  one process per file, kit_095): **368 files, 32,701 pass, 4 fail, 2 error, 313 skip**. The two errors are the known
  diagnostics (`test-workflow-perf101-diagnostic.R`, `test-workflow-ruv-pearson-diagnostic.R`); the four failures are
  `test-workflow-omics-binding-steps.R`'s characterisation (Compatibility, below). The skips equal KIT-090's checkpoint
  (313); every trx file has SKIP 0. The ticket's full-suite box stays unchecked until the characterisation decision.

## Compatibility

- **Frozen dump:** `test-workflow-omics-declarations.R` passes (68/0/0/0) with its fixture unmodified; the file differs only
  by the added binding name.
- **Registry test:** the two name lists gain the `trx_*` names. Its RNG scan now reaches KIT-093's `validateTrxCounts`,
  whose text matches only in the message `"sample(s) with a zero library"`; the test checks with the parse data that no
  call in it draws a random number and excludes that one name.
- **`test-workflow-omics-binding-steps.R`** (not in this ticket's files): its characterisation fails in four parts
  (`report_sources`, `types`, `steps`, `registry`). [characterisation-diff.txt](characterisation-diff.txt) shows each
  differs only by the transcriptomics entries (the new report sources, the binding's steps and registry rows,
  `TRANSCRIPTOMICS` accepted by the reopen, the session and the route, and that type named in the allowlist message); the
  `input_roles` and the lipid and metab `records` are identical. Regenerating the snapshot or filtering it to lipid and
  metab is a test-file change outside the write scope: **a decision is needed**.
- **Report overlay:** a lipid or metab record's overlay now also returns `ruv_mode` and `ruviii_m_source` as "not
  recorded" (the overlay's rule for an outcome source the record lacks).

## Divergences and notes

- **Outside the file list:** `DESCRIPTION` Collate gains `utils_workflow_trx_declarations.R` (one line; required by the
  ticket's new file).
- **The template (F05):** the class validity refuses an object with empty counts (a zero library), so
  `readWorkflowExecutionComponents(..., fields = "design_matrix")` cannot build a transcriptomics shell. The core's
  partitioned runner takes the binding's `read_template` when it has one; the trx one builds the shell without validity.
  `readWorkflowOmicsShell`/`Rows`/`Ids` have the same limit for this class: a hand-off to KIT-096/097.
- **The DE halves:** the fit (`trxDaFit`) and the test (`trxDaApply`) repeat `trxDaContrast`'s calls in their order and
  with their arguments, so one saved fit serves every contrast; identical to the reference fit, DGELRT and table.
- **The final PCA** reads the state alone (shared handler), so its display matrix is recomputed with `trxRuvs`'s
  `removeBatchEffect` call on the same inputs; identical to the product (test).
- **The k cap** is the kernel's rule, a fixed setting of the grid (the ticket's `max_k` follows the kernel).
- **The store handle** exposes `remember`, so a table a partitioned product writes itself is rooted with the product.
- **Test harness observation:** in a whole-file run, the session export's staged file is gone by the time the sync
  runner's `done` callback runs, though it exists right after the action; the test reads it right after the action. Not
  investigated further (the shared session-export machinery, not this ticket's code).
- **Open user decision kept:** more than two groups stay refused (KIT-093's `trx_design_groups_unsupported`).

## Hand-offs

- **KIT-096/097:** the design state `trx_raw_data_s4` with `recordTrxDesignProvenance`; the wrappers
  `runTrxReferenceExpressionFilter`, `runTrxReferenceRuvs(manager, trxRuvGridParameters(...), config)`,
  `runTrxReferenceDaAnalysis`, `runTrxReferenceRuviiiExport`, `acceptTrxRuviiiExport` (worker path); the publication
  operations and their parameters (`stage`, `figure`, `dataset`, `contrast`/`assay`/`qvalue`, `path`/`state_name`); the
  report sources; the shell limit above.
- **KIT-098:** the step map's bindings hash; journeys.
- **KIT-099:** the default variants.

Rollback: `git revert` of this ticket's commits; the evidence is kept.

## Closure (coordinator, 2026-11-06)

**The full suite's one real failure is resolved.** The four failing checks were all in the KIT-094 characterisation test (`test-workflow-omics-binding-steps.R`), whose snapshot was taken before transcriptomics existed. The test now compares the characterisation with exactly the transcriptomics entries removed:
- the trx report sources, steps, registry and `TRANSCRIPTOMICS` type rows;
- the supported-types message, read without the added name.

It also asserts that transcriptomics is now served by the reopen, the route and the steps. The snapshot fixture is unchanged. The file passes 51/0/0/0 on `kit_095`.

**Suite result:**
- At 815db998: 368 files, 32,701 pass, 4 fail (this test), 2 errors (the two known diagnostics).
- With this test fix: no unexpected failures.

**Recorded and accepted:**
- **DESCRIPTION Collate:** one line was added, outside the file list, because the new R file needs it to load.
- **Report overlay:** lipid and metab records now also carry `ruv_mode` and `ruviii_m_source`, reading "not recorded". These are report variables the lipid and metab templates don't use. The lipid and metab declaration tests (259/0 each) and the frozen dump pass. The next lipid and metab journeys (KIT-098's regression members) confirm it end to end.

**Hand-off to KIT-096:** `readWorkflowOmicsShell`, `readWorkflowOmicsRows` and `readWorkflowOmicsIds` can't build a `TranscriptomicsAssayData` shell, because the class's validity check rejects empty counts as a zero library. The binding reads its template with `trxReadTemplate`, and the tabs must do the same.
