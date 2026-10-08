# CORE-008 premise probe for item 3 (2026-10-08, in-process; not a counted run, no performance claim).
# Profiles one protein-grouped peptide QC step (imputation) on the scaled small fixture.
# Usage: Rscript groups-probe.R <n_proteins> <output dir>   (from the repository root)
args <- commandArgs(TRUE); proteins <- as.integer(args[[1L]]); out <- args[[2L]]
.libPaths(c("/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/perf-108/ordered-groups/browser-smoke.jEmE2F/library", .libPaths()))
suppressMessages(pkgload::load_all("/home/doktersmol/Projects/APAF/MultiScholaR", quiet = TRUE, export_all = TRUE, helpers = FALSE))
source("/home/doktersmol/Projects/APAF/MultiScholaR/tests/testthat/helper-workflow-checkpoint-fixtures.R")
recipe <- readWorkflowCheckpointRecipe("small", "/home/doktersmol/Projects/APAF/MultiScholaR"); recipe$n_proteins <- proteins
object <- makeWorkflowCheckpointPeptide(recipe)
object <- rollUpPrecursorToPeptide(object, core_utilisation = NA)
root <- file.path(out, paste0("store-", proteins)); unlink(root, recursive = TRUE); dir.create(root)
manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "probe", "DIA"), WorkflowState$new(audit_enabled = FALSE))
manager$.context$reference_states <- TRUE
manager$saveState("input", object, list(), "input")
workflow <- list2env(list(state_manager = manager, qc_params = list(peptide_qc = list())))
prof <- file.path(out, paste0("prof-", proteins, ".out"))
assign("config_list", list(peptideMissingValueImputation = list(proportion_missing_values = 0.5, exclusion_column = NULL)), envir = .GlobalEnv)
Rprof(prof, interval = 0.01)
t <- system.time(runDIAReferencePeptideGroups(workflow, "peptideMissingValueImputation",
    list(proportion_missing_values = 0.5, exclusion_column = NULL), "imputed", "probe"))[["elapsed"]]
Rprof(NULL)
s <- summaryRprof(prof)$by.total
pick <- function(f) if (f %in% rownames(s)) s[f, "total.pct"] else 0
cat(sprintf("RESULT proteins=%d peptide_rows=%d seconds=%.2f append_pct=%.1f writer_pct=%.1f next_group_pct=%.1f kernel_pct=%.1f\n", proteins, nrow(object@peptide_data), t,
  pick("\"appendWorkflowTableBatch\""), pick("\"appendWorkflowTableWriter\""), pick("\"cursor$next_group\""),
  max(pick("\"peptideMissingValueImputationHelper\""), pick("\"peptideMissingValueImputation\""))))
print(head(s[, c("total.pct", "self.pct")], 25))
