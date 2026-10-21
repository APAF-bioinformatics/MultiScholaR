# KIT-047 DPC check (in process, screen DIA-NN fixture). The DIA_LIMPA screen project is built through the entry points
# the GUI calls (as the KIT-046 route check does: owner attached for the route, import, FASTA, configuration, design,
# the reference Q-value filter and precursor rollup); then the three KIT-047 steps run on the owner through their
# in-process entry points (the registered declarations; the DPC step publishes the limpa QC figures); the same three
# kernels run in memory on the hydrated rollup state and the states are compared by content hash (exactness); the
# DPC engine runs a second time in memory (determinism within a process; across thread settings and engines by
# running this script under RAYON_NUM_THREADS / KIT047_ENGINE and comparing the JSON files); limpa runs on the same
# input for the alignment numbers; the args sizes of every protein state are measured; the analysis record is read.
# Usage: [KIT047_INSTALLED=1] [KIT047_ENGINE=fastlimpadpc|limpa] [RAYON_NUM_THREADS=n] R_LIBS=<library> \
#        Rscript docs/qa/workflow-execution/kit/kit-047/dpc-check-screen.R <label>
if (identical(Sys.getenv("KIT047_INSTALLED"), "1")) suppressMessages(library(MultiScholaR)) else suppressMessages(devtools::load_all(quiet = TRUE))
started <- Sys.time()
label <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[[1L]] else "tree"
engine <- Sys.getenv("KIT047_ENGINE", unset = "limpa")
repo <- normalizePath(".")
`%||%` <- function(a, b) if (is.null(a)) b else a
ns <- asNamespace("MultiScholaR"); api <- function(name) get(name, ns)
for (f in ls(ns, all.names = FALSE)) if (!exists(f, envir = globalenv(), inherits = TRUE) && is.function(get(f, ns))) assign(f, get(f, ns), envir = globalenv())
harness <- new.env(parent = globalenv()); sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), harness)
bindings <- jsonlite::read_json(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
corpus <- bindings$corpora$screen
report <- corpus$files[["report.tsv"]]$path; fasta <- corpus$files[["proteins.fasta"]]$path; design_path <- corpus$files[["design.tsv"]]$path
work <- file.path(bindings$harness$root, "kit047", label); unlink(work, recursive = TRUE); dir.create(work, recursive = TRUE)
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
out <- list(schema = "kit.dpc_check.v1", ticket = "KIT-047", label = label, date = format(Sys.Date()),
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

# 3. Exactness: the owner's states against the same kernels run in memory on the hydrated rollup state.
hashes <- function(states) {
    stats::setNames(lapply(states, function(state) {
        object <- manager$getState(state)
        if (methods::is(object, "PeptideQuantitativeData")) list(
            peptide_data = table_sha(object@peptide_data), peptide_matrix = matrix_sha(object@peptide_matrix),
            quantity = sha(object@peptide_data[[object@norm_quantity_column]]), is_logged = object@is_logged_data,
            rows = nrow(object@peptide_data), matrix_dim = dim(object@peptide_matrix), protein_groups = length(unique(object@peptide_data[[object@protein_id_column]])))
        else {
            r <- object@args$limpa_dpc_quant_results
            list(protein_quant_table = table_sha(object@protein_quant_table), E = matrix_sha(r$quantified_elist$E),
                standard_error = matrix_sha(r$standard_errors), design = table_sha(object@design_matrix),
                proteins = nrow(object@protein_quant_table), samples = ncol(object@protein_quant_table) - 1L,
                dpc = unname(as.numeric(r$dpc_parameters_used)), interpretation = r$slope_interpretation,
                engine = r$engine, engine_version = r$engine_version, peptides_used = r$total_peptides_used)
        }
    }), states)
}
withDIAStateScope(manager, function() {
    out$owner <<- hashes(c("log2_matrix", "peptide_accession_cleaned", "protein_dpc_quant"))
    rolled <- manager$getState("precursor_rollup")
    t0 <- Sys.time(); log2_obj <- log2TransformPeptideMatrix(calcPeptideMatrix(rolled)); timings$log2_memory <<- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    seqinr <- dplyr::rename(shiny::isolate(workflow$aa_seq_tbl_final), uniprot_acc = "database_id")
    t0 <- Sys.time(); acc_obj <- chooseBestProteinAccession(theObject = log2_obj, delim = ";", seqinr_obj = seqinr, seqinr_accession_column = "uniprot_acc", aggregation_method = "mean"); timings$accession_memory <<- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    t0 <- Sys.time(); dpc_obj <- suppressMessages(proteinMissingValueImputationLimpa(acc_obj, dpc_slope = 0.8, verbose = FALSE, chunk = 1000, engine = engine)); timings$dpc_memory <<- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    t0 <- Sys.time(); dpc_again <- suppressMessages(proteinMissingValueImputationLimpa(acc_obj, dpc_slope = 0.8, verbose = FALSE, chunk = 1000, engine = engine)); timings$dpc_memory_again <<- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    mem <- function(object) {
        if (methods::is(object, "PeptideQuantitativeData")) list(peptide_data = table_sha(object@peptide_data), peptide_matrix = matrix_sha(object@peptide_matrix),
            quantity = sha(object@peptide_data[[object@norm_quantity_column]]), rows = nrow(object@peptide_data))
        else { r <- object@args$limpa_dpc_quant_results; list(protein_quant_table = table_sha(object@protein_quant_table), E = matrix_sha(r$quantified_elist$E),
            standard_error = matrix_sha(r$standard_errors), proteins = nrow(object@protein_quant_table), dpc = unname(as.numeric(r$dpc_parameters_used))) }
    }
    out$memory <<- list(log2_matrix = mem(log2_obj), peptide_accession_cleaned = mem(acc_obj), protein_dpc_quant = mem(dpc_obj), protein_dpc_quant_again = mem(dpc_again))
    same <- function(state, field) identical(out$owner[[state]][[field]], out$memory[[state]][[field]])
    out$exact <<- list(log2_peptide_data = same("log2_matrix", "peptide_data"), log2_peptide_matrix = same("log2_matrix", "peptide_matrix"),
        accession_peptide_data = same("peptide_accession_cleaned", "peptide_data"), accession_peptide_matrix = same("peptide_accession_cleaned", "peptide_matrix"),
        dpc_protein_quant_table = same("protein_dpc_quant", "protein_quant_table"), dpc_E = same("protein_dpc_quant", "E"), dpc_standard_error = same("protein_dpc_quant", "standard_error"))
    out$determinism_in_process <<- list(E = identical(out$memory$protein_dpc_quant$E, out$memory$protein_dpc_quant_again$E),
        standard_error = identical(out$memory$protein_dpc_quant$standard_error, out$memory$protein_dpc_quant_again$standard_error),
        dpc_fit = identical(dpc_obj@args$limpa_dpc_quant_results$dpc_parameters_used, dpc_again@args$limpa_dpc_quant_results$dpc_parameters_used))
    # Where a difference lives: the fit alone (two fits on the same matrix) and the quantification alone (a fixed fit).
    y <- dpc_obj@args$limpa_dpc_quant_results$y_peptide_for_dpc
    engine_fn <- function(name) get(name, envir = asNamespace(engine), mode = "function")
    fit1 <- engine_fn("dpc")(y); fit2 <- engine_fn("dpc")(y)
    ids <- api(".ensurePeptideFeatureKeyMap")(acc_obj)@args$peptide_feature_key_map
    ids <- ids$active_protein_id[match(rownames(y), ids$feature_key)]
    elist <- list(E = y, genes = data.frame(protein.id = ids, peptide.id = rownames(y), stringsAsFactors = FALSE)); class(elist) <- "EList"
    q1 <- engine_fn("dpcQuant")(y = elist, protein.id = "protein.id", dpc = fit1$dpc, verbose = FALSE, chunk = 1000)
    q2 <- engine_fn("dpcQuant")(y = elist, protein.id = "protein.id", dpc = fit1$dpc, verbose = FALSE, chunk = 1000)
    out$determinism_parts <<- list(dpc_fit_identical = identical(fit1$dpc, fit2$dpc), dpc_fit = unname(as.numeric(fit1$dpc)),
        dpc_fit_abs_diff = abs(unname(as.numeric(fit1$dpc)) - unname(as.numeric(fit2$dpc))),
        quant_given_fit_E_identical = identical(q1$E, q2$E), quant_given_fit_SE_identical = identical(q1$other$standard.error, q2$other$standard.error),
        quant_given_fit_E_sha256 = matrix_sha(q1$E), quant_given_fit_SE_sha256 = matrix_sha(q1$other$standard.error))
    # Alignment against the reference engine on the same input (recorded, not a gate).
    if (out$installed$fastlimpadpc && out$installed$limpa) {
        other <- if (identical(engine, "fastlimpadpc")) "limpa" else "fastlimpadpc"
        t0 <- Sys.time(); other_obj <- suppressMessages(proteinMissingValueImputationLimpa(acc_obj, dpc_slope = 0.8, verbose = FALSE, chunk = 1000, engine = other)); timings[[paste0("dpc_memory_", other)]] <<- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        a <- dpc_obj@args$limpa_dpc_quant_results; b <- other_obj@args$limpa_dpc_quant_results
        common <- intersect(rownames(a$quantified_elist$E), rownames(b$quantified_elist$E))
        dE <- abs(a$quantified_elist$E[common, , drop = FALSE] - b$quantified_elist$E[common, , drop = FALSE])
        dSE <- abs(a$standard_errors[common, , drop = FALSE] - b$standard_errors[common, , drop = FALSE])
        out$alignment <<- list(engine = engine, reference = other, proteins_engine = nrow(a$quantified_elist$E), proteins_reference = nrow(b$quantified_elist$E),
            proteins_common = length(common), max_abs_dE = max(dE), median_abs_dE = stats::median(dE), max_abs_dSE = max(dSE), median_abs_dSE = stats::median(dSE),
            dpc_engine = unname(as.numeric(a$dpc_parameters_used)), dpc_reference = unname(as.numeric(b$dpc_parameters_used)),
            dpc_abs_diff = abs(unname(as.numeric(a$dpc_parameters_used)) - unname(as.numeric(b$dpc_parameters_used))))
    }
    # The readers of the owned entry on the hydrated owner state.
    hydrated <- readDIAExecutionComponents(manager)
    kept <- suppressMessages(filterMinNumPeptidesPerProtein(hydrated, peptides_per_protein_cutoff = 1L, peptidoforms_per_protein_cutoff = 2L))
    plots <- suppressMessages(generateLimpaQCPlots(hydrated, save_plots = FALSE, verbose = FALSE))
    out$readers <<- list(filter_min_peptides_proteins = nrow(kept@protein_quant_table), limpa_qc_plots = names(plots))
})
cat("EXACT:", paste(names(out$exact), unlist(out$exact), collapse = " "), "\n")
cat("DETERMINISM (two runs in process): E", out$determinism_in_process$E, "| SE", out$determinism_in_process$standard_error, "| fit", out$determinism_in_process$dpc_fit,
    "| parts: fit identical", out$determinism_parts$dpc_fit_identical, "fit diff", signif(out$determinism_parts$dpc_fit_abs_diff, 3),
    "| quant given fit E/SE identical", out$determinism_parts$quant_given_fit_E_identical, out$determinism_parts$quant_given_fit_SE_identical, "\n")
if (!is.null(out$alignment)) cat("ALIGNMENT", out$alignment$engine, "vs", out$alignment$reference, ": proteins", out$alignment$proteins_engine, "/", out$alignment$proteins_reference,
    "| max|dE|", signif(out$alignment$max_abs_dE, 3), "| max|dSE|", signif(out$alignment$max_abs_dSE, 3), "| dpc diff", signif(out$alignment$dpc_abs_diff, 3), "\n")
cat("READERS: filter kept", out$readers$filter_min_peptides_proteins, "| limpa QC plots", length(out$readers$limpa_qc_plots), "\n")

# 4. The args sizes of every state of the route and the owned value's file; the record.
store_dir <- harness$kitFindStore(file.path(work, "project"))
store_root <- dirname(store_dir)
sizes <- list()
withDIAStateScope(manager, function() {
    for (state in out$history) {
        reference <- getDIAStateReference(manager, state)
        spec <- tryCatch(readDIAStateMetadata(manager, reference)$data, saved_entity_unsupported = function(e) NULL)
        if (is.null(spec) || !isDIAComponentCheckpoint(spec)) next
        owned <- spec$slots$args$limpa_dpc_quant_results
        entry <- list(class = spec$class, inline_args_bytes = as.numeric(object.size(spec$slots$args)), record_bytes = as.numeric(object.size(spec)),
            support_sections = names(spec$support), owned_entry = isSavedWorkflowReference(owned))
        if (isSavedWorkflowReference(owned)) entry$owned_file_bytes <- as.numeric(owned$payload$bytes %||% NA_real_)
        fasta_support <- spec$support$chooseBestProteinAccession$seqinr_obj
        if (!is.null(fasta_support)) entry$fasta_table_support_bytes <- as.numeric(fasta_support$payload$bytes %||% NA_real_)
        sizes[[state]] <<- entry
    }
    record <- buildWorkflowAnalysisRecord(manager)
    step <- Filter(function(s) identical(s$state, "protein_dpc_quant"), record$steps)[[1L]]
    acc <- Filter(function(s) identical(s$state, "peptide_accession_cleaned"), record$steps)[[1L]]
    lg <- Filter(function(s) identical(s$state, "log2_matrix"), record$steps)[[1L]]
    out$record <<- list(dpc = list(operation = step$operation, variant = step$variant, outcomes = step$outcomes,
            parameters = lapply(step$parameters, function(p) list(name = p$name, provenance = p$provenance, value = decodeWorkflowRecordValue(p$value)))),
        accession = list(operation = acc$operation, parameters = lapply(acc$parameters, function(p) list(name = p$name, provenance = p$provenance))),
        log2 = list(operation = lg$operation, fixed = lg$fixed),
        software_packages = local({ listed <- stats::setNames(lapply(record$software$packages, `[[`, "version"), vapply(record$software$packages, `[[`, character(1), "name"))
            listed[intersect(c("limpa", "fastlimpadpc", "limma"), names(listed))] }),
        valid = validateWorkflowAnalysisRecord(record)$valid)
    fp <- harness$kitFingerprintStore(store_dir, api, states = c("log2_matrix", "peptide_accession_cleaned", "protein_dpc_quant"))
    out$fingerprints <<- lapply(fp$states, function(s) list(components = lapply(s$components, `[[`, "sha256"), summary_sha256 = s$summary_sha256, rows = s$rows))
})
out$sizes <- sizes
cat("SIZES:", paste(vapply(names(sizes), function(s) sprintf("%s inline %s B%s%s", s, format(sizes[[s]]$inline_args_bytes, big.mark = ","),
    if (!is.null(sizes[[s]]$owned_file_bytes)) sprintf(" owned file %s B", format(sizes[[s]]$owned_file_bytes, big.mark = ",")) else "",
    if (!is.null(sizes[[s]]$fasta_table_support_bytes)) sprintf(" fasta support %s B", format(sizes[[s]]$fasta_table_support_bytes, big.mark = ",")) else ""), ""), collapse = " | "), "\n")
cat("RECORD: dpc", out$record$dpc$operation, "| outcomes engine", out$record$dpc$outcomes$engine, out$record$dpc$outcomes$engine_version, "| params",
    paste(vapply(out$record$dpc$parameters, function(p) paste0(p$name, "=", p$value, ":", p$provenance), ""), collapse = " "),
    "| accession", out$record$accession$operation, "| log2", out$record$log2$operation, "| packages", paste(names(out$record$software_packages), unlist(out$record$software_packages), collapse = " "), "| valid", out$record$valid, "\n")
out$timings <- timings
closeDIAStateAdapter(manager)
out$seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
jsonlite::write_json(out, file.path(repo, "docs/qa/workflow-execution/kit/kit-047", paste0("dpc-check-", label, ".json")), auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA)
cat("TOTAL", round(out$seconds), "s\n")
