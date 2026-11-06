# KIT-098 plan (2026-11-06)

Ticket revision 1.1, with the user's decisions of 2026-11-06:
- **v1.8:** "Approve v1.8 as proposed".
- **Run order:** "Second notebook reference".

This plan is committed with protocol v1.8 before the first bindings commit. Nothing below changes `R/`.

## 1. Run order (settled)

**(a) is false.** Aligning runs by name gives 13 of 21 kernel products `identical` ([run-order.txt](run-order.txt)):
- edgeR's sums over samples change in the last digits;
- the RUVs SVD flips W_1's sign.

**(b) needs `R/`.**
- The lipid and metab builders order runs by `mixedsort`.
- Under the user's decision, (b) and (c) are not taken.

**(d), the user's decision: a second notebook reference.** [reference-gui-order.txt](reference-gui-order.txt) records it:
- **Render.** The unchanged notebook (`44d5f643…`) through KIT-040's driver (`8d1a2181…`, 900 s ceiling, one attempt):
  - exit 0, 13 chunks, 95 products, 48 s;
  - `git status` unchanged.
- **Input.** The bound `trx_screen` corpus, with `design_matrix.tab`'s rows permuted to the counts file's run order:
  - the permuted design has sha256 `f8763aaf…`;
  - every other file is byte-identical to the bound corpus.
- **Stores.**
  - The variant corpus is at `$P/kit/corpus/trx-screen-gui-order`.
  - The full products are at `$P/kit/trx-reference/reference-trx_screen-gui-order.rds`.
  - The fingerprint store is [reference-trx_screen_gui_order.rds](reference-trx_screen_gui_order.rds).
- **Confirmed before use.** KIT-093's kernels, run in counts-file order with no alignment, are `identical` to this reference on 21 of 21 products.

The journey is compared with this reference bit for bit:
- no alignment;
- no tolerance;
- KIT-040's reference unchanged.

## 2. Step list (`kitJourneySteps("transcriptomics")`)

The steps follow KIT-097's GUI: tab order, the managed sequences in `R/mod_trx_*.R`, and the e2e journey `tests/testthat/test-e2e-trx-journey.R`.

On DA and the session export: the DA tab reads the RUVs head directly; there is no session load (`R/mod_trx_da.R`, `runTrxDaAnalysis`). The session export is the summary tab's "Export Session State", after the report (e2e).

| Step | State it leaves | Action (`transcriptomics_workflow-` prefix) | Worker operations (managed) | Completion observed |
|---|---|---|---|---|
| `import` | none | `import-process_import` (shinyFiles `import-assay1_file`, format `trx_counts`) | none | "Data imported successfully!" |
| `design` | `trx_raw_data_s4` | `design-builder-*` (factors RPMI, Sera; runs by name; technical replicates per biological sample; contrast Sera vs RPMI; save) | none | "Design matrix and contrasts saved successfully!" |
| `expression_filter` | `trx_expression_filtered` | `qc-filter-apply_filter` | **`trx_expression_filter`**; `trx_qc_figures` (raw) | "Expression filter applied" |
| `qc_finalize` | `trx_qc_complete` | `qc-s4_finalize-finalize_qc` | `trx_filtering_metrics` | "QC complete!" |
| `normalisation` | `trx_ruvs` | `norm-run_normalization` | **`trx_tmm_logcpm`**, **`trx_ruv_grid`**, **`trx_ruvs`**; figures (raw, tmm, ruvs); `trx_ruv_publication`; `trx_composite_figure` | "Normalization complete!" |
| `da` | `trx_ruvs` (head unchanged) | `de-run_da_analysis` | **`trx_da_contrast`**; DA publications | `de-da_status` "Analysis Complete" |
| `ruviii_export` | `trx_ruvs` (head unchanged; a product) | `summary-ruviii_m_source`, `summary-run_ruviii_export` | **`trx_ruviii_export`**; figures (ruviii); composite | `summary-ruviii_status` "M = <source>" |
| `report` | none | `summary-save_workflow_args`, `copy_to_publication`, `generate_report` | integration export | "Study parameters saved", "Files copied to publication directory", "Report generated successfully!" |
| `session_export` | none | `summary-export_session_state` | none | "Session state exported to:" |

**One GUI action runs three worker steps.** The normalisation button runs `trx_tmm_logcpm`, `trx_ruv_grid` and `trx_ruvs` (`runTrxNormManaged`). The GUI cannot run them as separate journey steps. Per the ticket's alternative, two checks cover them:
- **Every history state by name.** The verifier compares the head's state history with the oracle's recorded history (bound as `<label>_history`). Every history state is also some step's state, so each is judged.
- **Products.** TMM/logCPM and the grid are products. They are judged as products attached to the `normalisation` step (`kitStepProducts`: `trx.tmm_logcpm:<sha12>`, `trx.ruv_grid:<sha12>`).

The RUVIII export product is attached to `ruviii_export`. The binding refuses an oracle for that step that carries no `trx.ruviii:` product.

**Worker steps (v1.8):** `trx_expression_filter`, `trx_tmm_logcpm`, `trx_ruv_grid`, `trx_ruvs`, `trx_da_contrast`, `trx_ruviii_export`.

A journey's decision needs one completed worker journal per listed step:
- `request.rds` and `journal.json` under `.workflow-export/wfx-*`;
- the scan of `kitGuiMemberDecision`, factored out and reused.

A listed step with no completed journal makes the member `invalid` (`worker_dispatch_missing`).

## 3. Settings

Settings are bound on the ticket keys and on `lifecycle_trx`, never the shared proteomics `journey$normalisation` or `journey$da`.

| Setting | Value | Source |
|---|---|---|
| Group levels, in order | RPMI, Sera | `design.tsv` order = notebook `c("RPMI", "Sera")` |
| Contrast | Sera vs RPMI (`Sera_vs_RPMI=groupSera-groupRPMI`) | `contrast_strings.tab` |
| RUV mode | `automatic` | notebook |
| Percentage range | 5 to 30 by 5 | notebook (rmd:195-258) |
| k penalty | 0.05 | notebook (rmd:243) |
| k cap | the kernel's rule (3 below 15 samples, else 4) | notebook (rmd:242) |
| Best k | automatic | notebook |
| DA FDR threshold | 0.05 | notebook |
| Fold-change line | 1 (display only) | notebook (rmd:685) |
| `ruviii_m_source` | `group` (098a, lifecycle), `technical_replicate` (098b) | user plan decision 2 |
| Technical replicates in the design | each biological sample's two runs (`design.tsv` `replicates`) | corpus; both journeys, so their upstream states are equal |

**Experiment label:** `kit_omics_lifecycle`, the driver's.

**The recorder** (`kit-098/oracle-record.R`) reads every setting from `bindings.json`, the same object the driver reads.

## 4. Lifecycle pair and fault

**Bindings:** `lifecycle_trx` with:
- build `kit_098`;
- corpus `trx_screen`;
- ceiling 1,800 s;
- the 098a settings;
- oracles `trx98a`;
- `in_memory` not declared (a trx in-memory declaration is refused).

**Stages and selectors:** stages `kit-lifecycle-trx` and `kit-lifecycle-trx-fail`; selectors `kit:lifecycle:trx:(complete|fail)`.

**The complete member:**
- runs with no launch record;
- is verified against the oracles;
- runs the decision's populations rule.

**The fail member:**
- carries the launch record (the guard needs it);
- uses fault `guard_after_qc_finalize`: `after_step = qc_finalize`, `failed_step = normalisation`.

The mechanism, after the QC finalise completes:
1. `application-guard-inject.rds` is written;
2. the application monitor stops the slice (`injected_guard_stop`);
3. the normalisation step, the first worker step after it, fails;
4. `head_unchanged` must be TRUE: the head stays `trx_qc_complete`.

## 5. Oracles: in-process recordings (protocol §7.1)

`kit-098/oracle-record.R`, modelled on `kit-032/oracle-record.R`, runs on the installed `kit_098`.

**Steps:**
1. **Import:** `runTrxImportProcessing` on the preview.
2. **Design:** the lipid builder server under `trxOmicProfile()`, with the driver's actions, then `runLipidDesignBuilderObserverShell`.
3. **Expression filter:** `runTrxReferenceExpressionFilter`.
4. **QC finalise:** `commitWorkflowOmicsStoreStep(trx_qc_finalize)`.
5. **Normalisation:** `runTrxReferenceRuvs(trxRuvGridParameters(...))`.
6. **DA:** `runTrxReferenceDaAnalysis`, with the managed tab's settings: the realised formula, FDR, lfc, and `trxDaContrastEntries`.
7. **RUVIII export:** `runTrxReferenceRuviiiExport(m_source)`.

The report and the session export commit nothing, so their populations are the head as it stands.

**Fingerprints:** `kitFingerprintStore` and `kitStepProducts(rendered = FALSE)`.

**The record:** schema `kit.oracle_record.v1`, plus:
- the `workflow_kit.R` sha256;
- R, BLAS and LAPACK;
- the variant per step (the registry default, `whole`);
- protocol v1.8;
- the history.

**Recordings:** two, `trx98a` (M = group) and `trx98b` (M = technical_replicate). Each takes at most 1,800 s.

**Publication products.** Products the managed GUI adds (figures, metrics, DA tables, the session and integration exports) are reported as additional, as on the lipid and metab routes. Every oracle product is judged.

**`kitOperationIdentityProducts`** stays a harness constant. It is extended for `trx_[a-z0-9_]+_science` only, with KIT-031's cause:
- the declaration's own scientific product is the run's entity, carrying the run's identity;
- the state's components and the step's attached `trx.*:<sha12>` products carry its content, and those are judged.

The two recordings, whose upstream states have equal content, show whether these differ. That is recorded in the outcome.

## 6. Notebook comparator (`kit-098/notebook-compare.R`)

**Input.** The comparator reads the **journey member's store**. It reflink-copies the member's project to a scratch directory first and opens the copy, so the member's runtime directory is never written. It compares each journey value with the GUI-order notebook reference (section 1). The comparison is `identical()` throughout, with no tolerance. Every product is reported as compared or excluded with its reason.

**Reads.** As KIT-095's `trxDeclJourneyProducts`:
- states by name through `readTrxWorkflowState` and `trxAssayInputs`;
- products by `product_type` from the head's products;
- the DA from the head's `da_results` manifest.

| Reference product | Journey source | Projection |
|---|---|---|
| `count_matrix_filtered` | state `trx_expression_filtered`, `trxAssayInputs()$counts` | none |
| `filter_effective_lib_size`, `filter_min_sample_size`, `filter_median_lib_size` | `trx.expression_filter_result` value | none |
| filter stage counts: `length(keep_expr)`, `sum(keep_expr)`, `sum(gene_vars_positive)` | that value's `counts` (before, after `filterByExpr`, after zero variance) | integers compared |
| `dge_filtered_samples` | `trx.tmm_logcpm` value `$samples` | none |
| `log_cpm_tmm` | `trx.tmm_logcpm` `$log_cpm` (`trxSavedMatrix`) | none |
| `ordered_features_idx`, `max_acceptable_k`, `df_summary`, `best_control_genes` | `trx.ruv_grid` value (`prelim_order`, `max_acceptable_k`, `grid_summary`, `control_genes`) | none |
| `grid_members`, `grid_cancor_data`, `grid_control_genes`, `best_res` | the grid value | KIT-095's shapes (`trxDeclGridMembers`, `pct_<n>` names, `best_res` fields) |
| `best_k` | `trx.ruviii` value `$k` | none |
| `w_matrix`, `log_cpm_corrected` | `trx.ruvs_result` value (`$W`, `$corrected` via `trxSavedMatrix`) | none |
| `ruvs_counts` | state `trx_ruvs` counts | none |
| `design_final`, `fit_final` | DA fit (`hydrateSavedWorkflowEntity`, `$fit`) | none |
| `design_formula_str` | the manifest's `parameters$formula` | none |
| `contrast_entry` | the first result's name | none |
| `contrasts_matrix`, `qlf`, `qlf_table`, `top_genes` | the first contrast's applied product | none |
| `ruviii_newY`, `ruviii_M`, `ruviii_ctl_mask` | `trx.ruviii` value | none; for 098b see below |
| PCA and RLE inputs (tmm, ruvs, post_ruviii) | the stage matrices above (`log_cpm_tmm`, `log_cpm_corrected`, `t(newY)`) | covered by those rows |

**Excluded, with reasons:**
- **`keep_expr`, `gene_vars`, `gene_vars_positive` as vectors:** not stored by the journey. The kept gene set is compared exactly through `count_matrix_filtered`'s row names, and the stage counts are compared.
- **`log_cpm_raw`:** not stored (KIT-095: `trx_qc_figures` computes it for stage raw).
- **`pca_*`/`rle_*` statistics, layers, boxes and limits:** a display value saved outside the head's products. Their inputs are compared above; KIT-093/095's tests cover the statistics.
- **Renderer products** (CCA, volcano, bar, grids, Glimma, summary text) and **`feature_meta`, `top_genes_annotated`, `top_up`/`top_down`/`top_de_combined`:** files and display values, not head products.
- **`tt_prelim`, `dge_raw_samples`, `ruvs_counts_identical`, `glimma_html_exists`:** notebook intermediates not produced by the journey.

**098b (M from technical replicates):**
- Every row above except the RUVIII rows is compared with the same GUI-order reference. The upstream states equal 098a's.
- **`ruviii_newY`** is compared with the direct call `ruv::RUVIII(Y = t(ref$log_cpm_tmm), M = ruv::replicate.matrix(tech_rep_group), ctl = ref$ruviii_ctl_mask, k = ref$best_k, return.info = TRUE)$newY` (KIT-093's oracle). `tech_rep_group` is the journey design's column in run order, and its partition must equal `design.tsv`'s.
- **`ruviii_M`** is compared with that `replicate.matrix`.
- **The post-RUVIII stage** is the newY row.

**Outcome:** `accepted` only when every compared product is `identical`.

## 7. Harness changes (test first, `tests/testthat/test-workflow-kit-harness.R`)

| Site | Change |
|---|---|
| `kitJourneySteps` | `transcriptomics` in the allowlist; the step list above |
| `kitMemberSeconds` | `kit-lifecycle-trx[-fail]`: `lifecycle_trx$ceiling_seconds` (at most 1,800 s) |
| `kitBindings` | carries `lifecycle_trx` |
| `kitCases` regex, `journey_case` | `step_map_trx` guard; the step-map switch refuses an unknown route; a trx case needs `oracles`, never `in_memory`, and every judged step and the RUVIII product have an oracle (`oracle_missing`), lifecycle included |
| `kitCases` | trx lifecycle block |
| `kitSelectorStage` | `kit:lifecycle:trx:(complete\|fail)` |
| `kitLifecycleDecision` | the stage by route (trx: `kit-lifecycle-trx`); an unknown route refused |
| `kitOperationIdentityProducts` | `trx_[a-z0-9_]+_science` (section 5) |
| `kitSingleJourneyDecision` | with `worker_steps` bound, every listed step needs a completed journal (`worker_dispatch_missing`); the scan is shared with `kitGuiMemberDecision` (`kitWorkerDispatches`) |
| `verifyKitMember` | trx: the history check; an oracle mismatch is raised as "oracle mismatch: …", the text the decisions classify as `rejected` |

**The `verifyKitMember` message change applies to every route.** Before it, a journey verifier's mismatch read "`!length(...)` is not TRUE", and the decisions classified it as `invalid`. Recorded as a divergence.

**Driver:**
- `runKitTrxImportActions` and `runKitTrxStep`, with the trx prefixes;
- `runKitJourney` dispatches the trx route;
- `runKitOmicsImportActions` and `runKitOmicsStep` refuse a route other than lipidomics or metabolomics.

**Coordinator:**
- the stage regex;
- trx takes `case$normalisation` and `case$driver_settings` (DA, M source) in place of the shared journey settings;
- `member_dir` is passed to `kitSingleJourneyDecision`.

**Runner:** usage text and stage regex.

**Generator:**
- **Corpus.** The `trx_screen` corpus with five files: counts, `design_matrix.tab`, `contrast_strings.tab`, `truth.json`, `design.tsv`. It asserts:
  - the manifest's sha256 equals `kit-040/reference-trx_screen.rds$corpus_manifest_sha256`;
  - the GUI-order store's manifest names `trx_screen`;
  - its counts and contrast hashes equal the bound corpus's.
- **Harness.** `step_map_trx`.
- **Lifecycle.** `lifecycle_trx`.
- **Tickets.** `098a` and `098b`, each with:
  - `kind = single_journey`, `whole_journey = TRUE`;
  - build `kit_098`;
  - ceiling 1,800 s;
  - `oracles` (`trx98a`, `trx98b`);
  - the normalisation and driver settings;
  - `worker_steps`;
  - `notebook_reference` (path, sha256 and store).
- **Oracle ingestion.** As at `:393-412`, labelled "in-process (protocol v1.8 §7.1)", plus `<label>_history`.

**Harness tests (each failing before its change):**
- an unknown route is refused;
- the trx step map is in the digest;
- the trx lifecycle ceiling is 1,800 s;
- the trx driver uses its own prefixes;
- a trx case with a missing oracle step, or without the RUVIII product, is refused (`oracle_missing`);
- a trx lifecycle case is not in memory, and its failure needs `head_unchanged`;
- a trx product outside the identity set is judged;
- a missing worker journal is detected.

## 8. Order, budget, decisions

1. Write this plan and protocol v1.8 (commit).
2. Make the test-first harness changes (commit).
3. Install `kit_098` from a pinned worktree at that commit.
4. Run the README chain check (items 1-6), reusing KIT-095's tests where they cover an item, and record which.
5. Make the two recordings (`trx98a`, `trx98b`), at most 1,800 s each.
6. Commit the bindings.
7. Run the members, in order, each with one attempt and a 1,800 s ceiling, recorded and uncounted:
   1. `kit:lifecycle:trx:complete`;
   2. `kit:lifecycle:trx:fail`;
   3. `kit:journey-098a:candidate:1`;
   4. `kit:journey-098b:candidate:1`.
8. Run the notebook comparison on 098a's and 098b's stores. It reads stores only and is not a run.

**Terminal outcomes:**
- `accepted`;
- `rejected` (an oracle or reference mismatch), recorded verbatim, never retried;
- `invalid`.

One bounded correction pass is allowed, only for a demonstrated instrumentation or harness failure.

**Claims:**
- The journey against the oracle is path equivalence only.
- The journey store against the notebook is the independent check.
- No memory, cross-build or real-data claim.
