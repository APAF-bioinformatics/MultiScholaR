# KIT-101 findings

Defects found while writing `design.md` (at `bb87a463`). Append only.

The user, 2026-11-07: "fix all findings". So F-01 to F-05 are fixed in the KIT-101 addendum, not left to their family tickets. Each fix is recorded under its entry, with its commit and its test.

| Id | Where | What | Repro | Was owned by | Status |
|---|---|---|---|---|---|
| F-01 | `R/mod_metab_design.R:228-252` | When the metab design-save preflight fails, `stop("Invalid metabolomics design: ...")` runs outside the `tryCatch` that would show a notification. The result is an observer error, and the "Processing Design Matrix" modal (:206-217) stays open. Lipid shows a notification and returns (`mod_lipid_design_builder_helpers.R:82-88`). Metab also writes the bus (:243-247) before it checks `source_dir` (:252). | `testServer` on `mod_metab_design_server` with a builder result whose design has a run that is not a sample column. Expected: a notification and the modal closed. Today: an observer error. | KIT-103 | fixing (addendum) |
| F-02 | `R/mod_prot_da_ui.R:346`; `R/mod_prot_da_render_volcano.R:41-42` | The proteomics DA UI defines `heatmap_show_dendro`, which no server code reads. The server reads `volcano_show_labels` and `volcano_label_top_n`, which the UI never defines, so their defaults (TRUE, 10) always apply. | grep for the three ids in `R/mod_prot_da*.R`. | KIT-106 | fixing (addendum) |
| F-03 | `R/mod_metab_da.R` (no tree-cut inputs; :358 refers to one); `R/mod_metab_da_registration_helpers.R:90-93,100,165-168` | The metab DA tree-cut defect. The server reads `heatmap_tree_cut_method`, `heatmap_n_clusters`, `heatmap_cut_height` and `heatmap_min_cluster_size`, but the UI defines none of them, so they are always NULL. Lipid defines them at `mod_lipid_da_ui.R:310-331`, with defaults none, 4, 10 and 5. | `grep -c 'ns("heatmap_tree_cut_method")' R/mod_metab_da.R` gives 0. | KIT-106 (parity) | fixing (addendum) |
| F-04 | `R/mod_lipid_summary_observer_helpers.R:35,64,94,251` | No R code calls the four `runLipidSummary*ObserverShell` functions. MCI-023 (`test-module-ci-lipid-summary.R:17,74`) tests them, not the live handlers, so it passes on dead code. | `grep -lw` for each name in `R/*.R` lists only the file that defines it. | KIT-108 | fixing (addendum) |
| F-05 | `R/mod_prot_norm_workflow_helpers.R:568,572-582,584-585`; `R/mod_prot_qc_peptide_log2.R:61-64` | Dead seams. `runProtNormNormalizationWorkflow` takes 14 parameters it never uses. `runPeptideLog2ApplyStep` takes four unused `Fn` parameters and calls only `runDIAReferencePeptideLog2Matrix`. Tests that inject these parameters never reach production code. | Read the function bodies; grep the tests that pass these arguments. | KIT-105 (norm), KIT-104 (log2) | fixing (addendum) |

F-02, F-03 and F-04 were checked by grep on 2026-11-07. F-01 and F-05 come from reading the code; the addendum's failing-first tests reproduce them.
