# KIT-021 acceptance criterion 2: the control build's in-memory path on the bound LFQ corpus (or any bound protein
# corpus). Run with the control library (kit_036, `core/browser-smoke.h6QeDJ`) on R_LIBS and the installed package:
#   R_LIBS=<proof>/core/browser-smoke.h6QeDJ/library Rscript docs/qa/workflow-execution/kit/kit-021/control-check.R lfq_300
# It imports through the same import helpers on a plain WorkflowState (no owner), builds the design as the builder's
# save does, runs the kernels the owner path's declarations wrap (accession cleanup, protein intensity filter,
# duplicate removal, replicate filter, between-sample normalisation with the journey's automatic RUV, Pearson pairs
# and selection) and the DA kernel per contrast as the DA declaration calls it, hashing the protein table at every
# state and the three DA tables per contrast with the canonical digest of kit-020/control-check.R.
# compare-check.R joins this with da-tables.R's hashes of the oracle run's saved tables (the kit_021 owner path).
suppressMessages(library(MultiScholaR))
started <- Sys.time()
args <- commandArgs(TRUE)
corpus_key <- if (length(args)) args[[1L]] else "lfq_300"
repo <- normalizePath(".")
ns <- asNamespace("MultiScholaR")
api <- function(name) get(name, ns)
`%||%` <- function(a, b) if (is.null(a)) b else a
for (f in ls(ns, all.names = FALSE)) if (!exists(f, envir = globalenv(), inherits = TRUE) && is.function(get(f, ns))) assign(f, get(f, ns), envir = globalenv())
bindings <- jsonlite::read_json(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
corpus <- bindings$corpora[[corpus_key]]
if (is.null(corpus)) stop("corpus ", corpus_key, " is not bound")
for (name in names(corpus$files)) stopifnot(identical(digest::digest(file = corpus$files[[name]]$path, algo = "sha256"), corpus$files[[name]]$sha256))
report <- corpus$files[[corpus$report_file]]$path; fasta <- corpus$files[["proteins.fasta"]]$path; design_path <- corpus$files[["design.tsv"]]$path
format <- corpus$format
journey <- bindings$journey
label <- journey$experiment_label
norm <- journey$normalisation
if (!is.null(corpus$normalisation)) norm <- utils::modifyList(norm, corpus$normalisation)
table_hash <- function(table) {
    table <- as.data.frame(table)
    table <- table[do.call(order, unname(as.list(table[, seq_len(min(2L, ncol(table))), drop = FALSE]))), order(names(table)), drop = FALSE]
    rownames(table) <- NULL
    digest::digest(lapply(table, function(column) if (is.numeric(column)) signif(column, 12) else as.character(column)), algo = "sha256")
}
work <- file.path(bindings$harness$root, "control", corpus_key); unlink(work, recursive = TRUE); dir.create(work, recursive = TRUE)
Sys.setenv(MULTISCHOLAR_TEST_MODE = "true")
project_dirs <- setupDirectories(base_dir = file.path(work, label), omic_types = "proteomics", label = label, force = TRUE, reuse_existing = TRUE)
assign("project_dirs", project_dirs, envir = globalenv())
paths <- project_dirs[[1L]]
paths$time_dir <- paths$time_dir %||% file.path(paths$publication_graphs_dir, "filtering_qc")
dir.create(paths$time_dir, recursive = TRUE, showWarnings = FALSE)
workflow <- new.env(parent = emptyenv())
workflow$state_manager <- WorkflowState$new(audit_enabled = FALSE)
workflow$tab_status <- list(); workflow$processing_log <- list()
stopifnot(!inherits(workflow$state_manager, "DIAStateAdapter"))
out <- list(corpus = corpus_key, corpus_sha256 = corpus$files[[1L]]$sha256, format = format, state_hashes = list(), rows = list(), steps = list())
step_started <- Sys.time()
tick <- function(id) { out$steps[[id]] <<- as.numeric(difftime(Sys.time(), step_started, units = "secs")); step_started <<- Sys.time(); cat("  ", id, "\n") }
# 1. Import on a plain WorkflowState through the import helpers the owner path also calls.
imported <- switch(format, maxquant = importMaxQuantData(report, use_lfq = TRUE, filter_contaminants = TRUE),
    pd_tmt = importProteomeDiscovererTMTData(report), stop("format ", format))
# The control build's MaxQuant reader keeps the flagged rows its sampled type guess misses (the reader repair of
# KIT-021, commit 88df9314, is on the candidate build only): the rows the repaired reader drops are dropped here by
# their MaxQuant prefixes so the kernels under comparison see the same table; the count is recorded.
if (identical(format, "maxquant")) {
    id_col <- names(imported$data)[[1L]]
    flagged <- grepl("^(CON__|REV__)", imported$data[[id_col]])
    out$control_reader_flagged_rows_dropped <- length(unique(imported$data[[id_col]][flagged]))
    imported$data <- imported$data[!flagged, , drop = FALSE]
}
out$workflow_type <- applyProtImportResultToWorkflow(workflow, imported, format, fasta, sanitizeNames = TRUE,
    logInfo = function(...) NULL, showNotification = function(...) NULL)
processProtImportFastaData(workflowData = workflow, fastaPath = fasta, organismName = "Homo sapiens", experimentPaths = paths,
    logInfo = function(...) NULL, logWarn = function(...) NULL)
load_args <- list(workflowData = workflow, experimentPaths = paths, useShinyFiles = FALSE,
    updateStatus = function(...) NULL, messageFn = function(...) NULL, showNotification = function(...) NULL,
    logInfo = function(...) NULL, logError = function(...) NULL)
if ("readConfig" %in% names(formals(loadProtImportConfigurationResources))) load_args$readConfig <- readProtGuiConfig
do.call(loadProtImportConfigurationResources, load_args)
if (is.null(workflow$config_list$globalParameters)) workflow$config_list$globalParameters <- list()
finalizeProtImportSetupState(workflow, imported, format, basename(report), basename(fasta), 9606, "Homo sapiens")
tick("import")
# 2. Design as the builder's save does it (every run assigned to its design group, replicates from 1), the checkpoint in memory.
design_in <- utils::read.delim(design_path, stringsAsFactors = FALSE, check.names = FALSE)
run_col <- imported$column_mapping$run_col %||% "Run"
runs <- gtools::mixedsort(as.character(unique(workflow$data_cln[[run_col]])))
raw <- tibble::tibble(Run = runs, group = NA_character_, factor1 = NA_character_, factor2 = NA_character_, factor3 = NA_character_,
    batch = NA_character_, replicates = NA_integer_, tech_reps = NA_integer_)
for (group in unique(design_in$group)) {
    selected <- design_in$sample[design_in$group == group]
    idx <- which(raw$Run %in% selected)
    if (length(idx) != length(selected)) stop("design samples not among the imported runs")
    raw$factor1[idx] <- group; raw$factor2[idx] <- ""; raw$factor3[idx] <- NA_character_
    raw$replicates[idx] <- seq.int(1L, length.out = length(selected))
    raw$group[idx] <- group
}
formula <- workflow$config_list$deAnalysisParameters$formula_string %||% "~ 0 + group"
# The one contrast the journey's design step adds (runKitDesignActions): the design file's second group minus its
# first, named as the builder's add-contrast action names it (`<group1>.vs.<group2>`).
design_groups <- unique(design_in$group)
contrast_data <- data.frame(contrast_name = paste0(design_groups[[2L]], ".vs.", design_groups[[1L]]),
    numerator = design_groups[[2L]], denominator = design_groups[[1L]], stringsAsFactors = FALSE)
results <- buildProtDesignSaveResultsPayload(designMatrix = raw, currentRemovedSamples = character(0),
    dataCln = workflow$data_cln, contrastData = contrast_data, configList = workflow$config_list, formulaString = formula)
hydrateProtDesignBuilderResults(results, workflow)
persistProtDesignBuilderArtifacts(results = results, workflowData = workflow, sourceDir = paths$source_dir)
out$design_state <- buildProtDesignStateCheckpoint(workflowData = workflow, workflowType = out$workflow_type, actionLabel = "Save Design")
object <- workflow$state_manager$getState("protein_s4_initial")
record_state <- function(state) { out$state_hashes[[state]] <<- table_hash(object@protein_quant_table); out$rows[[state]] <<- nrow(object@protein_quant_table) }
record_state("protein_s4_initial")
tick("design")
# 3. The four protein QC kernels as the declarations call them (computeDIAReferenceAccessionCleanup,
#    computeDIAReferenceProteinIntensityFilter, computeDIAReferenceProteinDuplicateRemoval, the replicate filter).
seqinr <- readRDS(file.path(paths$source_dir, "aa_seq_tbl_final.RDS"))
seqinr <- dplyr::rename(seqinr, uniprot_acc = "database_id")
object <- chooseBestProteinAccession(theObject = object, delim = ";", seqinr_obj = seqinr, seqinr_accession_column = "uniprot_acc",
    replace_zero_with_na = TRUE, aggregation_method = "mean")
record_state("protein_accession_cleaned")
tick("accession")
env <- new.env(parent = emptyenv()); env$config_list <- workflow$config_list
object <- api("updateMissingValueParameters")(object, min_reps_per_group = 2L, min_groups = 2L, env = env)
object <- api("updateConfigParameter")(object, function_name = "removeRowsWithMissingValuesPercent",
    parameter_name = "proteins_intensity_cutoff_percentile", new_value = 1L, env = env)
object <- removeRowsWithMissingValuesPercent(object)
record_state("protein_intensity_filtered")
tick("protein_intensity")
resolved <- api("resolveProteinDuplicateRows")(object@protein_quant_table, object@protein_id_column, base::get("mean"))
object@protein_quant_table <- resolved$table
record_state("duplicates_removed")
tick("duplicates")
object <- removeProteinsWithOnlyOneReplicate(object, NA, grouping_variable = "group")
record_state("protein_replicate_filtered")
tick("protein_replicate")
# 4. Normalisation and RUV as the normalization-ruv declaration's kernel runs them, with the journey's bound settings.
grouping <- getProtNormRuvGroupingVariable(norm$ruv$ruv_grouping_variable %||% "group")
settings <- c(list(norm_method = norm$norm_method, ruv_mode = norm$ruv_mode), norm$ruv[setdiff(names(norm$ruv), "ruv_grouping_variable")],
    list(ruv_percentage = 5, ruv_k = NA))
object <- normaliseBetweenSamples(object, normalisation_method = settings$norm_method)
record_state("normalized")
if (!identical(settings$ruv_mode, "skip")) {
    ruv <- api("resolveProtRuvParameters")(object, settings, grouping)
    out$ruv <- list(k = ruv$k, percentage = ruv$percentage, controls = length(ruv$controls))
    object <- api("ruvIII_C_Varying")(object, ruv_grouping_variable = grouping, ruv_number_k = ruv$k, ctrl = ruv$controls)
    object <- removeRowsWithMissingValuesPercent(theObject = object)
    record_state("ruv_corrected")
}
tick("normalisation_ruv")
pairs <- pearsonCorForSamplePairs(object, tech_rep_remove_regex = "pool", correlation_group = grouping)
object <- filterSamplesByProteinCorrelationThreshold(object, pearson_correlation_per_pair = pairs,
    min_pearson_correlation_threshold = norm$correlation_threshold)
record_state("correlation_filtered")
tick("correlation")
# 5. DA per contrast as the DA-contrast-output declaration's kernel calls the unchanged method (the DA tab's contrasts).
# The contrasts the design builder's save generated on the workflow (what the DA tab and the recorder run).
stopifnot(is.data.frame(workflow$contrasts_tbl), nrow(workflow$contrasts_tbl) > 0L)
contrasts <- as.character(workflow$contrasts_tbl$contrasts)
out$da <- list(contrasts = contrasts, table_hashes = list(), rows = list())
for (contrast in contrasts) {
    # The contrast as normaliseProtDaContrastsTable's full format names it (the kernel receives `full_format`).
    name <- gsub("-", "_vs_", gsub("-group", "-", sub("^group", "", contrast)))
    single <- data.frame(contrasts = paste0(name, "=", contrast), stringsAsFactors = FALSE)
    arguments <- list(contrasts_tbl = single, formula_string = journey$da$formula, da_q_val_thresh = journey$da$q_value,
        treat_lfc_cutoff = journey$da$treat_lfc, eBayes_trend = TRUE, eBayes_robust = TRUE, args_group_pattern = "(\\d+)",
        args_row_id = object@protein_id_column, qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue")
    this <- object
    for (section in c("differentialAbundanceAnalysis", "differentialAbundanceAnalysisHelper")) if (is.null(this@args[[section]])) this@args[[section]] <- arguments
    result <- differentialAbundanceAnalysis(this, contrasts_tbl = single, formula_string = journey$da$formula,
        da_q_val_thresh = journey$da$q_value, treat_lfc_cutoff = journey$da$treat_lfc, qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue")
    out$da$table_hashes[[contrast]] <- lapply(result[c("da_proteins_long", "da_proteins_wide", "contrasts_results_table")], table_hash)
    out$da$rows[[contrast]] <- nrow(result$da_proteins_long)
    cat("  DA", contrast, "rows", nrow(result$da_proteins_long), "\n")
}
tick("da")
out$seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
out$library <- list(path = find.package("MultiScholaR"), source_commit = utils::packageDescription("MultiScholaR")$`Config/MultiScholaR/SourceCommit` %||% NA_character_,
    version = as.character(utils::packageVersion("MultiScholaR")))
jsonlite::write_json(out, file.path(repo, "docs/qa/workflow-execution/kit/kit-021", paste0("control-check-", corpus_key, ".json")),
    auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = I(17))
cat("done", corpus_key, "| states", length(out$state_hashes), "| contrasts", length(contrasts), "|", round(out$seconds), "s\n")
