# KIT-201 R2b predictions (written before the run)

Run: base = `origin/v0.9-staging` 63143aa3, head = the merge 0a16eb23 (KIT-201 merged with 63143aa3), 108 files (R2's
101, `test-omics-enrich-core.R`, and the DA tests merged since dd769671: `test-kit174-ruv-grouping.R`,
`test-kit212-da-current.R`, `test-kit215-not-tested-lists.R`, `test-kit218-da-wrapper-tail.R`,
`test-kit220-da-threshold-consumers.R`, `test-kit222-replicate-keys.R`), 3 shards per side, `check: true`.
Source of the predictions: AUDIT-KIT-201-impl-20261125 section 4 (4) and I13.

**Predicted to move (head's value differs from stored; never edited, D5):**

| File | Leaves |
|---|---|
| `test-gui-entry-proteomics.R` | `prot_dia.json`: the DA long table `da_proteins_groupKO-groupWT_long_annot.tsv` at da, enrichment, summary_refusal, summary_copy, export, report, save; `prot_lfq.json`, `prot_tmt.json`: the DA long table at da, summary_copy, summary_export, summary_refusal, summary_report, summary_save; `prot_lfq_fragpipe.json`: the DA long table (6 steps) and the 12 NumSig leaves (`da_proteins_num_sig_da_molecules.tab` (6), `Publication_figures/NumSigDaMolecules/...tab` (3), `Publication_tables/da_proteomics_num_sig_da_molecules.tab` (3)), re-frozen by KIT-220 and moving again |
| `test-workflow-saved-da.R` | `savedDaW1Frozen()` W1 (a) and W2: `dia_non_head record`, `dia_reopened_rerun record` newly failing; changed text in the already-failing `dependencies`, `display`, `files`, `head_manifest`, `manifest`, `products`, `record` of dia_changed_arguments, dia_non_head, dia_one_contrast, dia_qvalue_only, dia_refused, dia_reopened_rerun, dia_repeat, dia_two_contrasts, dia_retry; `fits`, `states` and every limpa case unchanged |
| `test-omic-da-module.R` | `omicDaW1Frozen`: `prot_main` (22 leaves: run files/rows; select, table_down, table_limit, volcano_labels, volcano_ticked, download, save_heatmap, tree_dynamic, tree_height, tree_k, tree_none, dendrogram_off: files/values); `prot_load` (run files, rows); `prot_no_contrasts` (run_no_contrasts files, rows); KIT-112 W3 (rows, files) |
| `test-omic-enrich-module.R` | `omicEnrichW1Frozen`: the 40 leaves of I01 (gprofiler 5, clusterprofiler 6, kernel_gprofiler 5, kernel_clusterprofiler 4, organism_filter 9, refusals 6, restore 5), holding only the DA long table digest, `da_results_data` and the kernel's `da_data` |
| `test-omics-enrich-core.R` | the 11 leaves of I01 (the `da_data` digest of the run step) |

**Predicted verdicts:** `test-kit201-treat-moderation.R` `only-on-head` (all pass); `test-prot-ruv-factors-in-model.R`
`fixed` (staging's legacy digest test fails on runners; head's rewritten test passes); the five files above
`new-failure` or `message-changed`; every other file `same`, including the six newly added DA test files (they compare
consumers with the same run's table, not stored literals; INFERRED), every lipid/metab file (cutoff 0), `prot_dia_limpa`
and transcriptomics. `R CMD check`: no finding only on head. No `compare-error`.

**Flaky leaves (named as such; neither a prediction nor a stop):** I13's three: `test-omic-enrich-module.R` "refusals
table_removed/{da, open, run_table_removed} values" (failed on R3's base only); plan section 11's `gprofiler_plot`.

**Stop rule:** any other moved key stops the train: any NumSig leaf outside `prot_lfq_fragpipe.json`, any lipid/metab
statistic at cutoff 0, any limpa or transcriptomics value.

**Departures recorded by the coordinator (2026-11-25):** R3b is not run; the next full-suite checkpoint on staging
stands in for it (as for KIT-204, KIT-218, KIT-222). M1 is not repeated: since dd769671 staging changed none of
`func_lipid_da_model.R`, `func_metab_da_core_engine.R`, `func_prot_da_model_stats.R` or `fitDIAExecutionBaseModel`; the
one executed function that changed, `resolveProtTechnicalReplicateBlock` (KIT-222), returns nothing new unless the
design has a `tech_reps` column, which M1's designs lack (AUDIT-KIT-201-impl section 4 (6)).
