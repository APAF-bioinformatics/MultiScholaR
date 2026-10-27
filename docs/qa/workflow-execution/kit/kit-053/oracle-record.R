# KIT-053 oracle recorder (protocol v1.5 item 1) for the DIA_LIMPA route: an in-process run of the bound build on a
# bound DIA-NN corpus through the entry points the GUI's limpa route calls (the owner attached for the route, import,
# FASTA, configuration, the design builder's save, the reference Q-value filter and precursor rollup, the route's
# peptide stage (log2 matrix, peptide accession cleanup, DPC quantification on the shipped engine), its protein stage
# (between-sample normalisation none, peptide count filter), the shared correlation filter, the RUV covariates, the
# native session export and the reference DA with dpcDE), each state fingerprinted right after its step with the
# harness's canonical form (kitFingerprintStore). Read-only with respect to the member: it writes
# docs/qa/workflow-execution/kit/kit-053/oracles-limpa_<corpus>.rds, bound as limpa_<corpus>_<step>.
# Usage (installed bound build on R_LIBS):
#   KIT053_INSTALLED=1 R_LIBS=<library> Rscript docs/qa/workflow-execution/kit/kit-053/oracle-record.R <corpus key> [output rds]
if (identical(Sys.getenv("KIT053_INSTALLED"), "1")) suppressMessages(library(MultiScholaR)) else suppressMessages(devtools::load_all(quiet = TRUE))
args <- commandArgs(TRUE)
corpus_key <- args[[1L]]
repo <- normalizePath(".")
output <- if (length(args) > 1L) args[[2L]] else file.path(repo, "docs/qa/workflow-execution/kit/kit-053", paste0("oracles-limpa_", corpus_key, ".rds"))
started <- Sys.time()
`%||%` <- function(a, b) if (is.null(a)) b else a
ns <- asNamespace("MultiScholaR"); api <- function(name) get(name, ns)
for (f in ls(ns, all.names = FALSE)) if (!exists(f, envir = globalenv(), inherits = TRUE) && is.function(get(f, ns))) assign(f, get(f, ns), envir = globalenv())
harness <- new.env(parent = globalenv()); sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), harness)
bindings <- jsonlite::read_json(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
corpus <- bindings$corpora[[corpus_key]]
if (is.null(corpus)) stop("corpus ", corpus_key, " is not bound")
for (name in names(corpus$files)) stopifnot(identical(digest::digest(file = corpus$files[[name]]$path, algo = "sha256"), corpus$files[[name]]$sha256))
report <- corpus$files[[corpus$report_file]]$path; fasta <- corpus$files[["proteins.fasta"]]$path; design_path <- corpus$files[["design.tsv"]]$path
journey <- bindings$journey
label_dir <- journey$experiment_label
# KIT-056: `KIT_LIMPA_WORK_SUFFIX` keeps a second record's store apart (the correction mode beside the covariates mode).
work <- file.path(bindings$harness$root, "oracles", paste0("limpa_", corpus_key, Sys.getenv("KIT_LIMPA_WORK_SUFFIX"))); unlink(work, recursive = TRUE); dir.create(work, recursive = TRUE)
Sys.setenv(MULTISCHOLAR_TEST_MODE = "true")
fake_session <- function() {
    session <- new.env(parent = emptyenv()); session$userData <- new.env(parent = emptyenv())
    session$onFlush <- session$onFlushed <- session$onSessionEnded <- function(...) function() invisible(NULL)
    session
}
steps <- harness$kitJourneySteps("DIA_LIMPA")
populations <- list()
step_started <- Sys.time()
snapshot_step <- function(step) {
    store_dir <- harness$kitFindStore(file.path(work, "project"))
    stopifnot(!is.null(store_dir))
    head <- harness$kitReadHead(store_dir, api)
    state <- head$current_state
    populations[[step$id]] <<- list(step = step$id, expected_state = step$state, state = state, state_as_expected = identical(state, step$state),
        store = store_dir, head_sha256 = head$head_sha256, entity_id = head$states[[state]], rows = head$rows[[state]],
        product_ids = as.character(unlist(head$products)), seconds = as.numeric(difftime(Sys.time(), step_started, units = "secs")))
    cat(sprintf("  %-28s %-26s rows %s  %5.0f s\n", step$id, state, format(head$rows[[state]] %||% NA), populations[[step$id]]$seconds))
}
run <- function(step_id, expr) {
    step <- steps[[match(step_id, vapply(steps, `[[`, "", "id"))]]
    step_started <<- Sys.time()
    force(expr)
    snapshot_step(step)
}
project_dirs <- setupDirectories(base_dir = file.path(work, "project", label_dir), omic_types = "proteomics", label = label_dir, force = TRUE, reuse_existing = TRUE)
assign("project_dirs", project_dirs, envir = globalenv())
paths <- project_dirs[[1L]]
paths$time_dir <- paths$time_dir %||% file.path(paths$publication_graphs_dir, "filtering_qc")
dir.create(paths$time_dir, recursive = TRUE, showWarnings = FALSE)
workflow <- new.env(parent = emptyenv())
workflow$state_manager <- WorkflowState$new(audit_enabled = FALSE)
workflow$tab_status <- list(); workflow$processing_log <- list(); workflow$qc_params <- list()
session <- fake_session()
# 1. Import on the limpa route (the entry point's radio), as the import observer runs it.
run("import", {
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
})
manager <- workflow$state_manager
# The GUI's filtering-progress owner is the workflow bus (KIT-013); outside Shiny the renderer writes the global
# environment, so the bus mirrors it after every step's display (the next step's product records it as prior progress).
mirror_progress <- function() workflow$filtering_progress <- get0("filtering_progress", envir = globalenv(), inherits = FALSE)
publish <- function(result, step_name) {
    out <- publishDIAQvalueOutput(manager, result, list(
        data = result$filteringSummary %||% result$filteredS4@peptide_data, step_name = step_name, omic_type = "proteomics",
        experiment_label = label_dir, return_grid = TRUE, overwrite = TRUE, project_dirs = project_dirs), updateProteinFiltering)
    mirror_progress()
    out
}
publish_protein <- function(result, step_name) {
    out <- publishDIAProteinProgress(result, step_name, "proteomics", label_dir, updateProteinFiltering)
    mirror_progress()
    out
}
# 2. Design, as the builder's save does it (one contrast: the design's second group minus its first).
run("design", withDIAStateScope(manager, function() {
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
    # The driver's technical-replicate actions (runKitDesignActions, manual mode), as the builder applies them
    # (mod_prot_design_builder_action_helpers.R, assign_tech_reps): for each group and each replicate value of the
    # corpus design with more than one sample, those samples take that biological replicate number and technical
    # replicates 1..n in the order selected; the group's later replicate numbers shift down by the number consolidated.
    # Without this the recorder's design differed from the GUI's (KIT-053 lifecycle attempt 4: every design_matrix
    # mismatched, the RUV and DA blocks resolved differently).
    if ("replicates" %in% names(design_in)) {
        for (group in unique(design_in$group)) for (replicate in unique(design_in$replicates[design_in$group == group])) {
            samples <- design_in$sample[design_in$group == group & design_in$replicates == replicate]
            if (length(samples) < 2L) next
            selected <- which(raw$Run %in% samples)
            original <- unique(raw$replicates[selected])
            consolidated <- length(original) - 1 # a double, as the builder computes it: the shift below makes `replicates` double whenever a later sample is renumbered (the full corpus; not the screen corpus)
            highest <- max(original, na.rm = TRUE)
            adjust <- which(raw$group == group & !is.na(raw$group) & raw$replicates > highest & !raw$Run %in% samples)
            if (length(adjust)) raw$replicates[adjust] <- raw$replicates[adjust] - consolidated
            raw$replicates[selected] <- as.integer(replicate)
            raw$tech_reps[selected] <- stats::setNames(seq_along(samples), samples)[raw$Run[selected]]
        }
        # The GUI's saved design carries names on `Run` (equal to its values) whenever technical replicates were
        # assigned (observed on kit-journey-028's and kit-journey-053 attempt 4's raw_data_s4; absent on
        # kit-journey-021's design without technical replicates); the attribute holds no data, but the table
        # fingerprint serialises it. Mirrored here as observed; the GUI line that attaches it is not yet located.
        if (any(!is.na(raw$tech_reps))) names(raw$Run) <- raw$Run
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
}))
# 3. The shared peptide steps and the route's peptide stage through their in-process entry points (each publishing its QC outputs).
q <- workflow$config_list$srlQvalueProteotypicPeptideClean
run("qvalue", withDIAStateScope(manager, function() invisible(publish(runDIAReferenceQvalue(workflow, as.numeric(q$qvalue_threshold %||% 0.01),
    as.numeric(q$global_qvalue_threshold %||% 0.01), as.logical(as.integer(q$choose_only_proteotypic_peptide %||% 1L)),
    as.numeric(q$global_pg_qvalue_threshold %||% 0.01), q$confidence_provenance_mode %||% "diann_main_report", Sys.time), "2_qval_filtered"))))
run("rollup", withDIAStateScope(manager, function() invisible(publish(runDIAReferencePrecursorRollup(workflow), "4_precursor_rollup"))))
run("peptide_log2", withDIAStateScope(manager, function() invisible(publish(runDIAReferencePeptideLog2Matrix(workflow), "4_log2_matrix"))))
run("peptide_accession", withDIAStateScope(manager, function() invisible(publish(runDIAReferencePeptideAccessionCleanup(workflow,
    list(delimiter = ";", aggregation_method = "mean")), "7_annotation_cleanup"))))
dpc_engine <- resolveLimpaDpcEngine()
run("dpc_quant", withDIAStateScope(manager, function() invisible(publish_protein(suppressMessages(runDIAReferenceProteinDpcQuant(workflow,
    list(engine = dpc_engine, dpc_slope = 0.8, chunk = 1000))), "9_protein_dpc_quant"))))
# 4. The protein stage, the correlation filter and the RUV covariates as the tab runs them (the bound settings).
# The route's tab settings: the ticket's binding when bound, else the route's shipped defaults (no between-sample
# normalisation, the RUV covariates); the DIA journey's own defaults (cyclic loess, RUV-III correction) never apply here.
norm <- journey$normalisation
norm$norm_method <- "none"; norm$ruv_mode <- "covariates"
# The ticket's binding applies to the corpus it binds (none on the screen corpus, cyclic loess on the full corpus).
# KIT-056: `KIT_LIMPA_TICKET` names the ticket whose binding applies (default 053; 056 binds ruv_mode = correction).
limpa_ticket <- Sys.getenv("KIT_LIMPA_TICKET", "053")
ticket_norm <- if (identical(bindings$tickets[[limpa_ticket]]$corpus, corpus_key)) bindings$tickets[[limpa_ticket]]$normalisation
if (!is.null(ticket_norm)) norm <- utils::modifyList(norm, ticket_norm)
# KIT-056: `KIT_LIMPA_RUV_MODE` overrides the mode (the tree's correction-mode record before the ticket's build is bound).
if (nzchar(Sys.getenv("KIT_LIMPA_RUV_MODE"))) norm$ruv_mode <- Sys.getenv("KIT_LIMPA_RUV_MODE")
stopifnot(norm$ruv_mode %in% c("covariates", "correction"))
steps <- harness$kitJourneySteps("DIA_LIMPA", norm$ruv_mode)
run("protein_normalisation_limpa", withDIAStateScope(manager, function() invisible(publish_protein(runDIAReferenceProteinNormalisation(workflow,
    list(norm_method = norm$norm_method)), "10_protein_normalised"))))
run("peptide_count", withDIAStateScope(manager, function() invisible(publish_protein(runDIAReferenceProteinPeptideCountFilter(workflow,
    list(num_peptides_per_protein_thresh = 1, num_peptidoforms_per_protein_thresh = 2)), "11_peptide_count_filtered"))))
norm_data <- new.env(parent = emptyenv())
norm_data$qc_plot_paths <- NULL
grouping <- function() getProtNormRuvGroupingVariable(norm$ruv$ruv_grouping_variable %||% "group")
aesthetics <- function() getProtNormPlotAesthetics(NULL, NULL)
ruv_qc <- function(reference) {
    norm_data$qc_plot_paths <- generateProtNormRuvCorrectedQcArtifacts(ruvCorrectedS4 = reference,
        qcDir = norm_data$publication_paths$protein_qc_dir %||% paths$protein_qc_dir, aesthetics = aesthetics(), qcPlotPaths = norm_data$qc_plot_paths)
}
invisible(runDIAReferenceLimpaNormalisationPassThrough(list(ruv_mode = norm$ruv_mode), workflow, norm_data, function(reference) NULL))
run("correlation", invisible(runDIAReferenceCorrelation(workflow, norm_data, paths, norm$correlation_threshold, FALSE, "proteomics", label_dir, grouping, aesthetics)))
norm_data$correlation_filtering_complete <- TRUE
run("ruv_covariates", invisible(runDIAReferenceRuvCovariates(list(ruv_mode = norm$ruv_mode, color_variable = NULL, shape_variable = NULL),
    workflow, norm_data, paths, "proteomics", label_dir, ruv_qc, grouping)))
# 5. The native session export (no state) and the DA on the head with the design builder's contrasts (dpcDE on the route).
run("session_export", invisible(runProtNormExportSessionWorkflow(workflow, norm_data, list(norm_method = norm$norm_method %||% "none", ruv_mode = norm$ruv_mode),
    paths$source_dir, withProgressFn = function(message, value, expr) force(expr), incProgressFn = function(...) NULL, messageFn = function(...) NULL)))
run("da", withDIAStateScope(manager, function() {
    stopifnot(is.data.frame(workflow$contrasts_tbl), nrow(workflow$contrasts_tbl) > 0L)
    da_data <- new.env(parent = emptyenv()); da_data$current_s4_object <- getDIAStateReference(manager)
    da_data$contrasts_available <- workflow$contrasts_tbl$contrasts
    invisible(runDIAReferenceDaAnalysis(list(formula_string = journey$da$formula, da_q_val_thresh = journey$da$q_value, treat_lfc_cutoff = journey$da$treat_lfc),
        da_data, workflow, paths, contrasts = workflow$contrasts_tbl, annotations = workflow$uniprot_dat_cln, publish = TRUE))
}))
# KIT-065: the journey's last step, the summary tab's publication copy and report render, in process: the study
# parameters the tab bootstraps, the route's shipped template placed under scripts/, the copy, the render. The report is
# a rendered file (path and size recorded); the head is unchanged.
report_file <- NULL
run("report", {
    report_template <- resolveProtSummaryExpectedTemplate("dia_limpa")
    summary_dirs <- stats::setNames(list(paths), "proteomics") # the summary tab's projectDirs are keyed by omic type
    bootstrapProtSummaryCopyFallbackStudyParams(values = new.env(parent = emptyenv()), projectDirs = summary_dirs,
        omicType = "proteomics", experimentLabel = label_dir, description = "kit record", catFn = function(...) NULL)
    ensureProtSummaryReportTemplate(projectDirs = summary_dirs, templateFilename = report_template, omicType = "proteomics",
        showNotificationFn = function(...) NULL)
    copyToResultsSummary(omic_type = "proteomics", experiment_label = label_dir, contrasts_tbl = workflow$contrasts_tbl,
        design_matrix = workflow$design_matrix, force = TRUE)
    rendered <- RenderReport(omic_type = "proteomics", experiment_label = label_dir, rmd_filename = report_template)
    stopifnot(is.character(rendered), length(rendered) == 1L, file.exists(rendered))
    report_file <<- list(template = report_template, path = sub(paste0("^", normalizePath(paths$base_dir), "/"), "", normalizePath(rendered)),
        bytes = file.size(rendered))
    cat(sprintf("  %-20s %-28s %s bytes\n", "report", report_template, format(report_file$bytes)))
})
store_dir <- harness$kitFindStore(file.path(work, "project"))
head <- harness$kitReadHead(store_dir, api)
fingerprints <- harness$kitFingerprintStore(store_dir, api)
stopifnot(identical(fingerprints$head_sha256, head$head_sha256))
seen <- character()
for (id in names(populations)) {
    populations[[id]]$fingerprint <- fingerprints$states[[populations[[id]]$state]]
    new_ids <- setdiff(populations[[id]]$product_ids, seen)
    seen <- union(seen, populations[[id]]$product_ids)
    products <- fingerprints$products[intersect(new_ids, names(fingerprints$products))]
    products <- products[!vapply(products, function(p) identical(p$kind, "rendered_asset"), logical(1))]
    products <- unname(products[order(vapply(products, `[[`, character(1), "sha256"))])
    names(products) <- vapply(products, function(p) paste0(p$product_type %||% p$kind, ":", substr(p$sha256, 1L, 12L)), character(1))
    populations[[id]]$fingerprint$products <- c(populations[[id]]$fingerprint$products, products)
}
manifest <- withDIAStateScope(manager, function() resolveDIAStateValue(manager$.context, manager$.context$snapshot$da_results))
recorded <- list(schema = "kit.oracle_record.v1", corpus = corpus_key, corpus_sha256 = corpus$files[[1L]]$sha256, route = "DIA_LIMPA",
    build = utils::packageDescription("MultiScholaR")$`Config/MultiScholaR/SourceCommit` %||% NA_character_, library = find.package("MultiScholaR"),
    store = store_dir, head_sha256 = head$head_sha256, history = head$history, populations = populations,
    settings = list(normalisation = norm, da = journey$da, dpc_engine = dpc_engine), da_model = manifest$model,
    report = report_file, recorded_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
closeDIAStateAdapter(manager)
saveRDS(recorded, output, compress = FALSE, version = 3)
cat("recorded limpa", corpus_key, "->", output, "| states", length(populations), "| history", paste(head$history, collapse = ","),
    "| DA model", manifest$model$model, "k", manifest$model$k %||% NA, "|", round(recorded$seconds), "s\n")
