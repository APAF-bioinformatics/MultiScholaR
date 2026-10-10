.libPaths(c(Sys.getenv("R_LIBS"), .libPaths())); suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
P <- "/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5"
heads <- list.files(file.path(P, "kit/runtime/kit-journey.attempt4-ruv-k-degenerate/control-1/project"), "^dia-state-head[.]rds$", recursive = TRUE, full.names = TRUE, all.files = TRUE)
store_dir <- dirname(heads[grepl("[.]dia-native-", heads)][1])
identity <- readSavedStoreIdentity(file.path(store_dir, "store.json"))
manager <- newDIAStateAdapter(newWorkflowEntityStore(dirname(store_dir), identity$project_id, identity$workflow_id, max_cache_bytes = 0), reopen = TRUE)
ref <- manager$.context$snapshot$gui_references[["normalized-science"]]
obj <- withDIAStateScope(manager, function() readDIAExecutionComponents(manager, ref))
closeDIAStateAdapter(manager)
t0 <- proc.time()[["elapsed"]]
sweep <- suppressMessages(suppressWarnings(findBestNegCtrlPercentage(obj, percentage_range = 1:20, num_components_to_impute = 5, ruv_grouping_variable = "group", ruv_fdr_method = "BH", separation_metric = "max_difference", k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, verbose = FALSE)))
cat("RESULT sweep seconds:", round(proc.time()[["elapsed"]] - t0), "\n")
r <- sweep$optimization_results
for (i in seq_len(nrow(r))) cat(sprintf("RESULT pct %2d controls %3d k %2s sep %.4f comp %.4f %s\n", r$percentage_requested[i], r$realized_num_controls[i], as.character(r$best_k[i]), r$separation_score[i], r$composite_score[i], r$status[i]))
cat("RESULT WINNER percentage", sweep$best_percentage, "k", sweep$best_k, "controls", sum(sweep$best_control_genes_index), "\n")
corrected <- suppressMessages(suppressWarnings(ruvIII_C_Varying(obj, ruv_grouping_variable = "group", ruv_number_k = sweep$best_k, ctrl = sweep$best_control_genes_index)))
cat("RESULT RUV-III corrected table:", nrow(corrected@protein_quant_table), "x", ncol(corrected@protein_quant_table), "\n")
cleaned <- suppressMessages(removeRowsWithMissingValuesPercent(theObject = corrected))
cat("RESULT after missing-value cleanup:", nrow(cleaned@protein_quant_table), "x", ncol(cleaned@protein_quant_table), "\n")
