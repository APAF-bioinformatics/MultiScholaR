# KIT-020 audit F-01: the control build's in-memory path on the same fixtures. Run with the control library
# (kit_036, `core/browser-smoke.h6QeDJ`) on R_LIBS and the installed package, not load_all:
#   R_LIBS=<proof>/core/browser-smoke.h6QeDJ/library Rscript docs/qa/workflow-execution/kit/kit-020/control-check.R
# It builds ProteinQuantitativeData in memory through the GUI's design builder on a plain WorkflowState and runs
# the same kernels the owner path's declarations wrap (protein intensity filter, duplicate removal, replicate
# filter, between-sample normalisation with RUV skipped, Pearson pairs and selection, DA), hashing the protein
# table at every state and the DA tables with the same canonical digest as owner-check.R. compare-check.R joins
# the two JSON files.
suppressMessages(library(MultiScholaR))
started <- Sys.time()
repo <- normalizePath(".")
ns <- asNamespace("MultiScholaR")
api <- function(name) get(name, ns)
fixtures <- file.path(repo, "tests/testdata/e2e")
cases <- list(
    maxquant = list(file = file.path(fixtures, "prot_lfq", "proteinGroups.txt"), fasta = file.path(fixtures, "prot_lfq", "proteins.fasta"),
        design = file.path(fixtures, "prot_lfq", "design.tsv"), read = function(p) importMaxQuantData(p), type = "LFQ"),
    fragpipe = list(file = file.path(fixtures, "prot_lfq", "seed_combined_protein.tsv"), fasta = file.path(fixtures, "prot_lfq", "proteins.fasta"),
        design = file.path(fixtures, "prot_lfq", "design.tsv"), read = function(p) importFragPipeData(p), type = "LFQ"),
    pd_tmt = list(file = file.path(fixtures, "prot_tmt", "pd_tmt_peptides.tsv"), fasta = file.path(fixtures, "prot_tmt", "proteins.fasta"),
        design = file.path(fixtures, "prot_tmt", "design.tsv"), read = function(p) importProteomeDiscovererTMTData(p), type = "TMT"))
only <- commandArgs(TRUE)
if (length(only)) cases <- cases[only]
table_hash <- function(table) {
    table <- as.data.frame(table)
    table <- table[do.call(order, unname(as.list(table[, seq_len(min(2L, ncol(table))), drop = FALSE]))), order(names(table)), drop = FALSE]
    rownames(table) <- NULL
    digest::digest(lapply(table, function(column) if (is.numeric(column)) signif(column, 12) else as.character(column)), algo = "sha256")
}
work <- file.path(Sys.getenv("HOME"), ".cache", "multischolar-kit020-control-check")
unlink(work, recursive = TRUE); dir.create(work, recursive = TRUE)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
config <- readConfigFile(system.file("config", "config.ini", package = "MultiScholaR"), create_workers = FALSE)
results <- list()
run_case <- function(name, case) {
    out <- list(format = name, type = case$type)
    base <- file.path(work, name, "project"); dir.create(base, recursive = TRUE)
    project_dirs <- setupDirectories(base_dir = base, omic_types = "proteomics", label = name, force = TRUE, reuse_existing = TRUE)
    assign("project_dirs", project_dirs, envir = globalenv())
    paths <- project_dirs[[1L]]
    paths$time_dir <- paths$time_dir %||% file.path(paths$publication_graphs_dir, "filtering_qc")
    dir.create(paths$time_dir, recursive = TRUE, showWarnings = FALSE)
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- WorkflowState$new(audit_enabled = FALSE)
    workflow$config_list <- config
    workflow$config_list$globalParameters$workflow_type <- case$type
    workflow$qc_params <- list(protein_qc = list())
    stopifnot(!inherits(workflow$state_manager, "DIAStateAdapter"))
    imported <- case$read(case$file)
    type <- api("applyProtImportResultToWorkflow")(workflow, imported, name, case$fasta, logInfo = function(...) NULL,
        showNotification = function(...) NULL)
    out$workflow_type <- type
    design_in <- utils::read.delim(case$design, stringsAsFactors = FALSE)
    runs <- as.character(unique(workflow$data_cln[[imported$column_mapping$run_col]]))
    group <- vapply(runs, function(run) design_in$group[[which(endsWith(run, design_in$sample))[[1L]]]], "")
    design <- data.frame(Run = runs, group = unname(group), stringsAsFactors = FALSE)
    design$replicates <- stats::ave(seq_len(nrow(design)), design$group, FUN = seq_along)
    workflow$design_matrix <- design
    out$design_state <- api("buildProtDesignStateCheckpoint")(workflow, case$type, "Import")
    object <- workflow$state_manager$getState("protein_s4_initial")
    out$state_hashes <- list(protein_s4_initial = table_hash(object@protein_quant_table))
    # The kernels, as the owner path's compute functions call them (utils_workflow_prot_declarations.R,
    # utils_workflow_dia_qc_rollup.R, utils_workflow_dia_products.R), with the owner check's parameters.
    env <- new.env(parent = emptyenv()); env$config_list <- workflow$config_list
    section <- "removeRowsWithMissingValuesPercent"
    object <- api("updateMissingValueParameters")(object, min_reps_per_group = 2L, min_groups = 2L, env = env)
    object <- api("updateConfigParameter")(object, function_name = section, parameter_name = "proteins_intensity_cutoff_percentile",
        new_value = 1L, env = env)
    object <- removeRowsWithMissingValuesPercent(object)
    out$state_hashes$protein_intensity_filtered <- table_hash(object@protein_quant_table)
    resolved <- api("resolveProteinDuplicateRows")(object@protein_quant_table, object@protein_id_column, base::get("mean"))
    object@protein_quant_table <- resolved$table
    out$state_hashes$duplicates_removed <- table_hash(object@protein_quant_table)
    object <- removeProteinsWithOnlyOneReplicate(object, NA, grouping_variable = "group")
    out$state_hashes$protein_replicate_filtered <- table_hash(object@protein_quant_table)
    object <- normaliseBetweenSamples(object, normalisation_method = "cyclicloess")
    out$state_hashes$normalized <- table_hash(object@protein_quant_table)
    pairs <- pearsonCorForSamplePairs(object, tech_rep_remove_regex = "pool", correlation_group = "group")
    object <- filterSamplesByProteinCorrelationThreshold(object, pearson_correlation_per_pair = pairs, min_pearson_correlation_threshold = 0.5)
    out$state_hashes$correlation_filtered <- table_hash(object@protein_quant_table)
    single <- data.frame(contrasts = "KO_vs_WT=groupKO-groupWT", stringsAsFactors = FALSE)
    arguments <- list(contrasts_tbl = single, formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0,
        eBayes_trend = TRUE, eBayes_robust = TRUE, args_group_pattern = "(\\d+)", args_row_id = object@protein_id_column,
        qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue")
    for (s in c("differentialAbundanceAnalysis", "differentialAbundanceAnalysisHelper")) if (is.null(object@args[[s]])) object@args[[s]] <- arguments
    result <- differentialAbundanceAnalysis(object, contrasts_tbl = single, formula_string = "~ 0 + group", da_q_val_thresh = 0.05,
        treat_lfc_cutoff = 0, qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue")
    out$da <- list(rows = nrow(result$da_proteins_long), table_hashes = list(`groupKO-groupWT` = lapply(
        result[c("da_proteins_long", "da_proteins_wide", "contrasts_results_table")], table_hash)))
    out
}
for (name in names(cases)) {
    results[[name]] <- tryCatch(run_case(name, cases[[name]]), error = function(e) list(format = name, error = conditionMessage(e)))
    r <- results[[name]]
    if (!is.null(r$error)) { cat(name, "ERROR:", r$error, "\n"); next }
    cat(name, "| type", r$workflow_type, "| design", r$design_state, "| hashes", length(r$state_hashes), "| DA rows", r$da$rows, "\n")
}
results$seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
results$library <- list(path = find.package("MultiScholaR"), source_commit = utils::packageDescription("MultiScholaR")$`Config/MultiScholaR/SourceCommit`)
jsonlite::write_json(results, file.path(repo, "docs/qa/workflow-execution/kit/kit-020", "control-check.json"), auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = I(17))
cat("done", round(results$seconds), "s\n")
