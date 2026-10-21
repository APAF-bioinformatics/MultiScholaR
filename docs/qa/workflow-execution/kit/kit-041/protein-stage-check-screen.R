# KIT-041 protein-stage check (in process, screen DIA-NN fixture). The DIA_LIMPA screen project is built through the entry
# points the GUI calls up to the DPC quantification exactly as KIT-047's check does (that part of this script is copied from
# kit-047/dpc-check-screen.R); then the route's protein stage runs on the owner through its in-process entry points (the
# between-sample normalisation with the shipped default `none`, the peptide count filter, the shared correlation filter at
# 0.5 and the RUV-III factors as covariates); the same kernels run in memory on the hydrated DPC state and the states are
# compared by content hash; the vendored script optimiser runs on the same correlation-filtered object; the record is read.
# Usage: [KIT041_INSTALLED=1] R_LIBS=<library> Rscript docs/qa/workflow-execution/kit/kit-041/protein-stage-check-screen.R <label>
# (The first part below is the KIT-047 build; its own checks are not repeated here.)
# The DIA_LIMPA screen project is built through the entry points
# the GUI calls (as the KIT-046 route check does: owner attached for the route, import, FASTA, configuration, design,
# the reference Q-value filter and precursor rollup); then the three KIT-047 steps run on the owner through their
# in-process entry points (the registered declarations; the DPC step publishes the limpa QC figures); the same three
# kernels run in memory on the hydrated rollup state and the states are compared by content hash (exactness); the
# DPC engine runs a second time in memory (determinism within a process; across thread settings and engines by
# running this script under RAYON_NUM_THREADS / KIT047_ENGINE and comparing the JSON files); limpa runs on the same
# input for the alignment numbers; the args sizes of every protein state are measured; the analysis record is read.
# Usage: [KIT047_INSTALLED=1] [KIT047_ENGINE=fastlimpadpc|limpa] [RAYON_NUM_THREADS=n] R_LIBS=<library> \
#        Rscript docs/qa/workflow-execution/kit/kit-047/dpc-check-screen.R <label>
if (identical(Sys.getenv("KIT041_INSTALLED"), "1")) suppressMessages(library(MultiScholaR)) else suppressMessages(devtools::load_all(quiet = TRUE))
started <- Sys.time()
label <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[[1L]] else "tree"
engine <- "limpa"
repo <- normalizePath(".")
`%||%` <- function(a, b) if (is.null(a)) b else a
ns <- asNamespace("MultiScholaR"); api <- function(name) get(name, ns)
for (f in ls(ns, all.names = FALSE)) if (!exists(f, envir = globalenv(), inherits = TRUE) && is.function(get(f, ns))) assign(f, get(f, ns), envir = globalenv())
harness <- new.env(parent = globalenv()); sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), harness)
bindings <- jsonlite::read_json(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
corpus <- bindings$corpora$screen
report <- corpus$files[["report.tsv"]]$path; fasta <- corpus$files[["proteins.fasta"]]$path; design_path <- corpus$files[["design.tsv"]]$path
work <- file.path(bindings$harness$root, "kit041", label); unlink(work, recursive = TRUE); dir.create(work, recursive = TRUE)
Sys.setenv(MULTISCHOLAR_TEST_MODE = "true")
sha <- function(x) digest::digest(x, algo = "sha256")
table_sha <- function(df) { df <- as.data.frame(df); sha(list(names(df), lapply(df, function(col) { attributes(col) <- NULL; col }))) }
matrix_sha <- function(m) sha(list(dimnames(m), as.vector(m)))
cpu_flags <- function() {
    lines <- readLines("/proc/cpuinfo", warn = FALSE)
    flags <- strsplit(sub("^flags\\s*:\\s*", "", grep("^flags", lines, value = TRUE)[[1L]]), " ")[[1L]]
    list(model = sub("^model name\\s*:\\s*", "", grep("^model name", lines, value = TRUE)[[1L]]), n_flags = length(flags),
        simd = intersect(c("sse4_2", "avx", "avx2", "fma", "avx512f", "avx512bw", "avx512vl", "avx_vnni", "avx512_bf16"), flags),
        sha256 = sha(sort(flags)))
}
fake_session <- function() {
    session <- new.env(parent = emptyenv()); session$userData <- new.env(parent = emptyenv())
    session$onFlush <- session$onFlushed <- session$onSessionEnded <- function(...) function() invisible(NULL)
    session
}
out <- list(schema = "kit.protein_stage_check.v1", ticket = "KIT-041", label = label, date = format(Sys.Date()),
    head_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    library = list(path = find.package("MultiScholaR"), source_commit = utils::packageDescription("MultiScholaR")$`Config/MultiScholaR/SourceCommit` %||% NA_character_),
    engine = engine, rayon_num_threads = Sys.getenv("RAYON_NUM_THREADS", unset = "<host default>"), host_threads = parallel::detectCores(),
    installed = list(limpa = requireNamespace("limpa", quietly = TRUE), fastlimpadpc = requireNamespace("fastlimpadpc", quietly = TRUE)),
    versions = list(limpa = tryCatch(as.character(utils::packageVersion("limpa")), error = function(e) NA_character_),
        fastlimpadpc = tryCatch(as.character(utils::packageVersion("fastlimpadpc")), error = function(e) NA_character_),
        fastlimpadpc_pinned_sha = sub(".*fast_limpa_dpc@([0-9a-f]+).*", "\\1", grep("fast_limpa_dpc@", readLines(file.path(repo, "DESCRIPTION")), value = TRUE)[[1L]]),
        rustc = tryCatch(system2("rustc", "--version", stdout = TRUE), error = function(e) NA_character_),
        cargo = tryCatch(system2("cargo", "--version", stdout = TRUE), error = function(e) NA_character_)),
    cpu = cpu_flags())
if (!out$installed$limpa) stop("limpa is not installed in this library")
if (identical(engine, "fastlimpadpc") && !out$installed$fastlimpadpc) stop("fastlimpadpc is not installed in this library")

# 1. The DIA_LIMPA screen project through the entry points (the KIT-046 route check's path).
label_dir <- "dia_baseline_smoke"
project_dirs <- setupDirectories(base_dir = file.path(work, "project", label_dir), omic_types = "proteomics", label = label_dir, force = TRUE, reuse_existing = TRUE)
assign("project_dirs", project_dirs, envir = globalenv())
paths <- project_dirs[[1L]]
paths$time_dir <- paths$time_dir %||% file.path(paths$publication_graphs_dir, "filtering_qc")
dir.create(paths$time_dir, recursive = TRUE, showWarnings = FALSE)
workflow <- new.env(parent = emptyenv())
workflow$state_manager <- WorkflowState$new(audit_enabled = FALSE)
workflow$tab_status <- list(); workflow$processing_log <- list(); workflow$qc_params <- list()
session <- fake_session()
initializeDIANativeImport(workflow, paths, session, "diann", "DIA_LIMPA")
manager <- workflow$state_manager
stopifnot(isDIANativeSession(manager), identical(manager$workflow_type, "DIA_LIMPA"))
withDIAStateScope(manager, function() {
    inputs <- c(report, fasta); names(inputs) <- c("diann_report", "fasta")
    registerDIAInputFiles(manager, inputs)
    imported <- importDIANNData(report)
    applyProtImportResultToWorkflow(workflow, imported, "diann", fasta, sanitizeNames = TRUE, logInfo = function(...) NULL, showNotification = function(...) NULL)
    releaseDIAWorkingReference(manager, workflow$data_tbl)
    processProtImportFastaData(workflowData = workflow, fastaPath = fasta, organismName = "Homo sapiens", experimentPaths = paths,
        logInfo = function(...) NULL, logWarn = function(...) NULL)
    loadProtImportConfigurationResources(workflowData = workflow, experimentPaths = paths, useShinyFiles = FALSE, readConfig = readProtGuiConfig,
        updateStatus = function(...) NULL, messageFn = function(...) NULL, showNotification = function(...) NULL,
        logInfo = function(...) NULL, logError = function(...) NULL)
    finalizeProtImportSetupState(workflow, getDIAImportMetadata(manager), "diann", basename(report), basename(fasta), 9606, "Homo sapiens")
    completeDIAInputFiles(manager)
})
publish <- function(result, step_name) publishDIAQvalueOutput(manager, result, list(
    data = result$filteringSummary %||% result$filteredS4@peptide_data, step_name = step_name, omic_type = "proteomics",
    experiment_label = label_dir, return_grid = TRUE, overwrite = TRUE, project_dirs = project_dirs), updateProteinFiltering)
withDIAStateScope(manager, function() {
    design_in <- utils::read.delim(design_path, stringsAsFactors = FALSE, check.names = FALSE)
    runs <- gtools::mixedsort(as.character(getDIAImportMetadata(manager)$samples$Run))
    raw <- tibble::tibble(Run = runs, group = NA_character_, factor1 = NA_character_, factor2 = NA_character_, factor3 = NA_character_,
        batch = NA_character_, replicates = NA_integer_, tech_reps = NA_integer_)
    for (group in unique(design_in$group)) {
        selected <- design_in$sample[design_in$group == group]
        idx <- which(raw$Run %in% selected)
        if (length(idx) != length(selected)) stop("design samples not among the imported runs")
        raw$factor1[idx] <- group; raw$factor2[idx] <- ""; raw$factor3[idx] <- NA_character_
        raw$replicates[idx] <- seq.int(1L, length.out = length(selected)); raw$group[idx] <- group
    }
    design_groups <- unique(design_in$group)
    contrast_data <- data.frame(contrast_name = paste0(design_groups[[2L]], ".vs.", design_groups[[1L]]),
        numerator = design_groups[[2L]], denominator = design_groups[[1L]], stringsAsFactors = FALSE)
    formula <- workflow$config_list$deAnalysisParameters$formula_string %||% "~ 0 + group"
    results <- buildProtDesignSaveResultsPayload(designMatrix = raw, currentRemovedSamples = character(0),
        dataCln = workflow$data_cln, contrastData = contrast_data, configList = workflow$config_list, formulaString = formula)
    results$data_reference <- saveDIAOwnedValue(results$data_cln, "design-assigned")
    results["data_cln"] <- list(NULL)
    hydrateProtDesignBuilderResults(results, workflow)
    persistProtDesignBuilderArtifacts(results = results, workflowData = workflow, sourceDir = paths$source_dir)
    buildProtDesignStateCheckpoint(workflowData = workflow, workflowType = manager$workflow_type, actionLabel = "Save Design")
    q <- workflow$config_list$srlQvalueProteotypicPeptideClean
    qv <- runDIAReferenceQvalue(workflow, as.numeric(q$qvalue_threshold %||% 0.01), as.numeric(q$global_qvalue_threshold %||% 0.01),
        as.logical(as.integer(q$choose_only_proteotypic_peptide %||% 1L)), as.numeric(q$global_pg_qvalue_threshold %||% 0.01),
        q$confidence_provenance_mode %||% "diann_main_report", Sys.time)
    invisible(publish(qv, "2_qval_filtered"))
    ru <- runDIAReferencePrecursorRollup(workflow)
    invisible(publish(ru, "4_precursor_rollup"))
})
cat("PROJECT:", manager$workflow_type, "| history", paste(manager$getHistory(), collapse = ","), "|", round(as.numeric(difftime(Sys.time(), started, units = "secs"))), "s\n")

# 2. The three steps on the owner through their in-process entry points.
dpc_parameters <- list(engine = engine, dpc_slope = 0.8, chunk = 1000)
accession_parameters <- list(delimiter = ";", aggregation_method = "mean")
timings <- list()
withDIAStateScope(manager, function() {
    t0 <- Sys.time(); lg <- runDIAReferencePeptideLog2Matrix(workflow); timings$log2_owner <<- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    invisible(publish(lg, "4_log2_matrix"))
    t0 <- Sys.time(); ac <- runDIAReferencePeptideAccessionCleanup(workflow, accession_parameters); timings$accession_owner <<- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    invisible(publish(ac, "7_annotation_cleanup"))
    t0 <- Sys.time(); dq <- suppressMessages(runDIAReferenceProteinDpcQuant(workflow, dpc_parameters)); timings$dpc_owner <<- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    out$owner_text <<- list(log2 = lg$resultText, accession = ac$resultText, dpc = dq$resultText)
})
out$history <- manager$getHistory()
cat("OWNER:", paste(out$history, collapse = ","), "| log2", round(timings$log2_owner, 1), "s | accession", round(timings$accession_owner, 1), "s | dpc", round(timings$dpc_owner, 1), "s\n")
limpa_figures <- list.files(paths$protein_qc_dir, pattern = "^limpa_.*[.]png$")
out$limpa_figures <- limpa_figures
cat("LIMPA FIGURES:", length(limpa_figures), paste(limpa_figures, collapse = ","), "\n")


# 3. The protein stage on the owner through the in-process entry points.
norm_data <- new.env(parent = emptyenv())
norm_data$qc_plot_paths <- NULL
aesthetics <- function() list(color_var = "group", shape_var = "group")
grouping <- function() "group"
ruv_qc <- function(reference) {
    norm_data$qc_plot_paths <- generateProtNormRuvCorrectedQcArtifacts(ruvCorrectedS4 = reference,
        qcDir = norm_data$publication_paths$protein_qc_dir %||% paths$protein_qc_dir, aesthetics = aesthetics(), qcPlotPaths = norm_data$qc_plot_paths)
}
stage <- list()
withDIAStateScope(manager, function() {
    t0 <- Sys.time(); nm <- runDIAReferenceProteinNormalisation(workflow, list(norm_method = "none")); timings$normalisation_owner <<- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    t0 <- Sys.time(); pc <- runDIAReferenceProteinPeptideCountFilter(workflow, list(num_peptides_per_protein_thresh = 1, num_peptidoforms_per_protein_thresh = 2)); timings$peptide_count_owner <<- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    norm_data$ruv_normalized_obj <- getDIAStateReference(manager)
    t0 <- Sys.time(); cr <- runDIAReferenceCorrelation(workflow, norm_data, paths, 0.5, FALSE, "proteomics", label_dir, grouping, aesthetics); timings$correlation_owner <<- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    t0 <- Sys.time(); rc <- runDIAReferenceRuvCovariates(list(ruv_mode = "covariates", color_variable = "group", shape_variable = "group"), workflow, norm_data, paths, "proteomics", label_dir, ruv_qc, grouping); timings$ruv_covariates_owner <<- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    stage$texts <<- list(normalisation = nm$resultText, peptide_count = pc$resultText)
    stage$ruv <<- norm_data$ruv_optimization_result[c("best_percentage", "best_k", "best_separation_score", "best_composite_score", "max_acceptable_k", "sample_size", "replicate_structure", "block_reason")]
    stage$control_count <<- sum(norm_data$control_genes_index, na.rm = TRUE)
})
out$history <- manager$getHistory()
cat("OWNER:", paste(out$history, collapse = ","), "| normalisation", round(timings$normalisation_owner, 1), "s | peptide count", round(timings$peptide_count_owner, 1), "s | correlation", round(timings$correlation_owner, 1), "s | ruv covariates", round(timings$ruv_covariates_owner, 1), "s\n")
cat("RUV:", paste(names(stage$ruv), unlist(lapply(stage$ruv, function(x) if (is.numeric(x)) signif(x, 6) else x)), collapse = " "), "| controls", stage$control_count, "\n")
out$ruv <- stage$ruv; out$control_count <- stage$control_count
published <- list.files(paths$protein_qc_dir, pattern = "^ruv_(covariates|adjusted|corrected)_.*[.](png|tsv)$")
out$published_files <- published
cat("PUBLISHED:", length(published), paste(published, collapse = ","), "\n")

# 4. Exactness: the owner's four states against the same kernels in memory on the hydrated DPC state.
protein_hash <- function(object) list(protein_quant_table = table_sha(object@protein_quant_table), design = table_sha(object@design_matrix),
    design_columns = names(object@design_matrix), proteins = nrow(object@protein_quant_table), samples = nrow(object@design_matrix),
    elist_E = if (!is.null(object@args$limpa_dpc_quant_results$quantified_elist$E)) matrix_sha(object@args$limpa_dpc_quant_results$quantified_elist$E) else NA_character_,
    elist_columns = colnames(object@args$limpa_dpc_quant_results$quantified_elist$E))
withDIAStateScope(manager, function() {
    owner <- lapply(stats::setNames(c("protein_normalised", "peptide_count_filtered", "correlation_filtered", "ruv_covariates"),
        c("protein_normalised", "peptide_count_filtered", "correlation_filtered", "ruv_covariates")), function(state) protein_hash(manager$getState(state)))
    dpc <- manager$getState("protein_dpc_quant")
    t0 <- Sys.time()
    normaliseBetweenSamplesFn <- normaliseBetweenSamples
    normalised <- normaliseBetweenSamplesFn(dpc, normalisation_method = "none")
    filterMinNumPeptidesPerProteinFn <- filterMinNumPeptidesPerProtein
    filtered <- suppressMessages(filterMinNumPeptidesPerProteinFn(normalised, num_peptides_per_protein_thresh = 1, num_peptidoforms_per_protein_thresh = 2, verbose = FALSE))
    pairs <- suppressMessages(pearsonCorForSamplePairs(filtered, tech_rep_remove_regex = "pool", correlation_group = "group"))
    filterSamplesFn <- filterSamplesByProteinCorrelationThreshold
    correlated <- suppressMessages(filterSamplesFn(filtered, pearson_correlation_per_pair = pairs, min_pearson_correlation_threshold = 0.5))
    correlated <- subsetLimpaResultsToSamples(correlated)
    sweep <- suppressMessages(suppressWarnings(findBestNegCtrlPercentageCapped(correlated, verbose = FALSE)))
    factors <- computeRuvCovariates(correlated, sweep$best_k, sweep$best_control_genes_index, "group")
    with_w <- correlated
    runs <- as.character(with_w@design_matrix[[with_w@sample_id]])
    for (i in seq_len(factors$k)) with_w@design_matrix[[paste0("W", i)]] <- unname(factors$W[runs, i])
    timings$memory_chain <<- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    memory <- list(protein_normalised = protein_hash(normalised), peptide_count_filtered = protein_hash(filtered),
        correlation_filtered = protein_hash(correlated), ruv_covariates = protein_hash(with_w))
    out$owner <<- owner; out$memory <<- memory
    same <- function(state, field) identical(owner[[state]][[field]], memory[[state]][[field]])
    out$exact <<- list(normalised_table = same("protein_normalised", "protein_quant_table"), normalised_design = same("protein_normalised", "design"),
        peptide_count_table = same("peptide_count_filtered", "protein_quant_table"), peptide_count_elist = same("peptide_count_filtered", "elist_E"),
        correlation_table = same("correlation_filtered", "protein_quant_table"), correlation_design = same("correlation_filtered", "design"),
        correlation_elist = same("correlation_filtered", "elist_E"),
        covariates_table = same("ruv_covariates", "protein_quant_table"), covariates_design = same("ruv_covariates", "design"),
        covariates_design_columns = identical(owner$ruv_covariates$design_columns, memory$ruv_covariates$design_columns))
    out$w_columns <<- setdiff(owner$ruv_covariates$design_columns, owner$correlation_filtered$design_columns)
    out$matrix_equals_correlation <<- identical(owner$ruv_covariates$protein_quant_table, owner$correlation_filtered$protein_quant_table)
    out$elist_in_step_with_design <<- setequal(owner$correlation_filtered$elist_columns, as.character(manager$getState("correlation_filtered")@design_matrix$Run))
    out$sweep_memory <<- list(best_percentage = sweep$best_percentage, best_k = sweep$best_k, controls = sum(sweep$best_control_genes_index), composite = sweep$best_composite_score)
    # The vendored script optimiser on the same object.
    vendored <- new.env(parent = globalenv())
    sys.source(file.path(repo, "docs/qa/workflow-execution/kit/limpa-workflow/findBestNegCtrlPercentageCustom.R"), vendored)
    t0 <- Sys.time()
    script <- suppressMessages(suppressWarnings(vendored$findBestNegCtrlPercentageCustom(correlated, percentage_range = seq(1, 30, by = 1), num_components_to_impute = 5,
        ruv_grouping_variable = "group", ruv_fdr_method = "BH", ruv_qval_cutoff = 0.05, separation_metric = "max_difference", k_penalty_weight = 0.5,
        adaptive_k_penalty = TRUE, verbose = FALSE)))
    timings$vendored_script <<- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    columns <- c("percentage", "separation_score", "best_k", "composite_score", "num_controls", "valid_plot")
    out$script_vs_capped <<- list(percentage = identical(script$best_percentage, sweep$best_percentage), k = identical(as.integer(script$best_k), sweep$best_k),
        controls = identical(script$best_control_genes_index, sweep$best_control_genes_index),
        trace = isTRUE(all.equal(as.data.frame(script$optimization_results[columns]), as.data.frame(sweep$optimization_results[columns]), check.attributes = FALSE)),
        script = list(percentage = script$best_percentage, k = script$best_k), script_commit = sub("^.*commit ([0-9a-f]+).*$", "\\1", readLines(file.path(repo, "docs/qa/workflow-execution/kit/limpa-workflow/findBestNegCtrlPercentageCustom.R"), n = 1L)))
    # Sizes and the record.
    sizes <- list()
    for (state in out$history) {
        spec <- tryCatch(readDIAStateMetadata(manager, getDIAStateReference(manager, state))$data, saved_entity_unsupported = function(e) NULL)
        if (is.null(spec) || !isDIAComponentCheckpoint(spec)) next
        owned <- spec$slots$args$limpa_dpc_quant_results
        sizes[[state]] <- list(inline_args_bytes = as.numeric(object.size(spec$slots$args)), owned_entry = isSavedWorkflowReference(owned),
            owned_entity = if (isSavedWorkflowReference(owned)) substr(owned$entity_id, 1, 12) else NA_character_,
            owned_file_bytes = if (isSavedWorkflowReference(owned)) as.numeric(owned$payload$bytes %||% NA_real_) else NA_real_)
    }
    out$sizes <<- sizes
    record <- buildWorkflowAnalysisRecord(manager)
    steps <- lapply(stats::setNames(c("protein_normalised", "peptide_count_filtered", "correlation_filtered", "ruv_covariates"), c("protein_normalised", "peptide_count_filtered", "correlation_filtered", "ruv_covariates")), function(st) {
        s <- Filter(function(x) identical(x$state, st), record$steps)[[1L]]
        list(operation = s$operation, variant = s$variant, outcomes = s$outcomes, parameters = lapply(s$parameters, function(p) list(name = p$name, provenance = p$provenance, value = decodeWorkflowRecordValue(p$value))))
    })
    out$record <<- list(steps = steps, valid = validateWorkflowAnalysisRecord(record)$valid)
})
cat("EXACT:", paste(names(out$exact), unlist(out$exact), collapse = " "), "\n")
cat("W COLUMNS:", paste(out$w_columns, collapse = ","), "| matrix equals correlation state:", out$matrix_equals_correlation, "| EList in step with design:", out$elist_in_step_with_design, "\n")
cat("SCRIPT vs CAPPED:", paste(names(out$script_vs_capped)[1:4], unlist(out$script_vs_capped[1:4]), collapse = " "), "| script", out$script_vs_capped$script$percentage, out$script_vs_capped$script$k, "| commit", out$script_vs_capped$script_commit, "\n")
cat("SIZES:", paste(vapply(names(out$sizes), function(s) sprintf("%s inline %s B owned %s", s, format(out$sizes[[s]]$inline_args_bytes, big.mark = ","), out$sizes[[s]]$owned_entity), ""), collapse = " | "), "\n")
cat("RECORD:", paste(vapply(names(out$record$steps), function(s) sprintf("%s=%s[%s]", s, out$record$steps[[s]]$operation, paste(vapply(out$record$steps[[s]]$parameters, function(p) paste0(p$name, ":", p$provenance), ""), collapse = ",")), ""), collapse = " | "), "| valid", out$record$valid, "\n")
out$timings <- timings
closeDIAStateAdapter(manager)
out$seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
jsonlite::write_json(out, file.path(repo, "docs/qa/workflow-execution/kit/kit-041", paste0("protein-stage-check-", label, ".json")), auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA)
cat("TOTAL", round(out$seconds), "s\n")
