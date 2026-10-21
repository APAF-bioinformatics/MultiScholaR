# KIT-021 oracle recorder (protocol v1.5 item 1): an in-process run of the bound build on a bound protein corpus,
# through the same entry points the GUI's protein route calls (import, the design builder's save, the four protein
# QC declarations through their registry handlers and the shared completion, the reference normalisation,
# correlation, the native session export and the reference DA), each state fingerprinted right after its step with
# the harness's canonical form (kitFingerprintStore). Read-only with respect to the member: it writes
# docs/qa/workflow-execution/kit/kit-021/oracles-<corpus>.rds, which the bindings generator binds as <corpus>_<step>.
# Usage (installed bound build on R_LIBS):
#   KIT021_INSTALLED=1 R_LIBS=<library> Rscript docs/qa/workflow-execution/kit/kit-021/oracle-record.R <corpus key> [output rds]
if (identical(Sys.getenv("KIT021_INSTALLED"), "1")) suppressMessages(library(MultiScholaR)) else suppressMessages(devtools::load_all(quiet = TRUE))
args <- commandArgs(TRUE)
corpus_key <- args[[1L]]
repo <- normalizePath(".")
output <- if (length(args) > 1L) args[[2L]] else file.path(repo, "docs/qa/workflow-execution/kit/kit-021", paste0("oracles-", corpus_key, ".rds"))
started <- Sys.time()
ns <- asNamespace("MultiScholaR")
api <- function(name) get(name, ns)
for (f in ls(ns, all.names = FALSE)) if (!exists(f, envir = globalenv(), inherits = TRUE) && is.function(get(f, ns))) assign(f, get(f, ns), envir = globalenv())
harness <- new.env(parent = globalenv()); sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), harness)
bindings <- jsonlite::read_json(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
corpus <- bindings$corpora[[corpus_key]]
if (is.null(corpus)) stop("corpus ", corpus_key, " is not bound")
for (name in names(corpus$files)) stopifnot(identical(digest::digest(file = corpus$files[[name]]$path, algo = "sha256"), corpus$files[[name]]$sha256))
report <- corpus$files[[corpus$report_file]]$path
fasta <- corpus$files[[corpus$fasta_file %||% "proteins.fasta"]]$path; design_path <- corpus$files[[corpus$design_file %||% "design.tsv"]]$path
stopifnot(file.exists(report), file.exists(fasta), file.exists(design_path))
format <- corpus$format
journey <- bindings$journey
label <- journey$experiment_label
work <- file.path(bindings$harness$root, "oracles", corpus_key); unlink(work, recursive = TRUE); dir.create(work, recursive = TRUE)
Sys.setenv(MULTISCHOLAR_TEST_MODE = "true")
fake_session <- function() {
    session <- new.env(parent = emptyenv()); session$userData <- new.env(parent = emptyenv())
    session$onFlush <- session$onFlushed <- session$onSessionEnded <- function(...) function() invisible(NULL)
    session
}
project_dirs <- setupDirectories(base_dir = file.path(work, label), omic_types = "proteomics", label = label, force = TRUE, reuse_existing = TRUE)
assign("project_dirs", project_dirs, envir = globalenv())
paths <- project_dirs[[1L]]
paths$time_dir <- paths$time_dir %||% file.path(paths$publication_graphs_dir, "filtering_qc")
dir.create(paths$time_dir, recursive = TRUE, showWarnings = FALSE)
workflow <- new.env(parent = emptyenv())
workflow$state_manager <- WorkflowState$new(audit_enabled = FALSE)
workflow$tab_status <- list(); workflow$processing_log <- list()
session <- fake_session()
steps <- harness$kitJourneySteps("protein")
populations <- list()
snapshot_step <- function(step) {
    # The driver's after-step snapshot (kitSnapshotAfterStep): every step records the head it left current
    # (session_export leaves the correlation state; a skipped RUV commits `normalized` where the route names
    # `ruv_corrected`) and the product ids the head holds; the fingerprints are taken once at the end as the driver does.
    store_dir <- harness$kitFindStore(file.path(work, label))
    stopifnot(!is.null(store_dir))
    head <- harness$kitReadHead(store_dir, api)
    state <- head$current_state
    populations[[step$id]] <<- list(step = step$id, expected_state = step$state, state = state, state_as_expected = identical(state, step$state),
        store = store_dir, head_sha256 = head$head_sha256, entity_id = head$states[[state]], rows = head$rows[[state]],
        product_ids = as.character(unlist(head$products)), seconds = as.numeric(difftime(Sys.time(), step_started, units = "secs")))
    cat(sprintf("  %-20s %-28s rows %s\n", step$id, state, format(head$rows[[state]] %||% NA)))
}
run <- function(step_id, expr) {
    step <- steps[[match(step_id, vapply(steps, `[[`, "", "id"))]]
    step_started <<- Sys.time()
    force(expr)
    snapshot_step(step)
}
# 1. Import, as the import observer runs it: the owner, the inputs by role, the reader, the owned table, the FASTA,
#    the configuration, the input completion.
run("import", {
    initializeDIANativeImport(workflow, paths, session, format, "iq")
    manager <- workflow$state_manager
    stopifnot(isDIANativeSession(manager))
    withDIAStateScope(manager, function() {
        inputs <- c(report, fasta); names(inputs) <- c(paste0(format, "_report"), "fasta")
        registerDIAInputFiles(manager, inputs)
        imported <- switch(format, maxquant = importMaxQuantData(report, use_lfq = TRUE, filter_contaminants = TRUE),
            pd_tmt = importProteomeDiscovererTMTData(report), stop("format ", format))
        applyProtImportResultToWorkflow(workflow, imported, format, fasta, sanitizeNames = TRUE,
            logInfo = function(...) NULL, showNotification = function(...) NULL)
        releaseDIAWorkingReference(manager, workflow$data_tbl)
        processProtImportFastaData(workflowData = workflow, fastaPath = fasta, organismName = "Homo sapiens", experimentPaths = paths,
            logInfo = function(...) NULL, logWarn = function(...) NULL)
        loadProtImportConfigurationResources(workflowData = workflow, experimentPaths = paths, useShinyFiles = FALSE, readConfig = readProtGuiConfig,
            updateStatus = function(...) NULL, messageFn = function(...) NULL, showNotification = function(...) NULL,
            logInfo = function(...) NULL, logError = function(...) NULL)
        finalizeProtImportSetupState(workflow, getDIAImportMetadata(manager), format, basename(report), basename(fasta), 9606, "Homo sapiens")
        completeDIAInputFiles(manager)
    })
})
manager <- workflow$state_manager
# 2. Design, as the builder's save does it: every run assigned to its design group in one action per group
#    (factor 1 = group, replicates numbered from 1 within the assignment), then the save payload and the checkpoint.
run("design", withDIAStateScope(manager, function() {
    design_in <- utils::read.delim(design_path, stringsAsFactors = FALSE, check.names = FALSE)
    runs <- gtools::mixedsort(as.character(getDIAImportMetadata(manager)$samples$Run))
    raw <- tibble::tibble(Run = runs, group = NA_character_, factor1 = NA_character_, factor2 = NA_character_, factor3 = NA_character_,
        batch = NA_character_, replicates = NA_integer_, tech_reps = NA_integer_)
    for (group in unique(design_in$group)) {
        selected <- design_in$sample[design_in$group == group]
        idx <- which(raw$Run %in% selected)
        if (length(idx) != length(selected)) stop("design samples not among the imported runs; runs: ", paste(head(runs, 8), collapse = ", "),
            " | design: ", paste(head(selected, 8), collapse = ", "))
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
    results$data_reference <- saveDIAOwnedValue(results$data_cln, "design-assigned")
    results["data_cln"] <- list(NULL)
    hydrateProtDesignBuilderResults(results, workflow)
    persistProtDesignBuilderArtifacts(results = results, workflowData = workflow, sourceDir = paths$source_dir)
    buildProtDesignStateCheckpoint(workflowData = workflow, workflowType = manager$workflow_type, actionLabel = "Save Design")
}))
# 3. The four protein QC steps through their registry handlers and the shared completion (the managed step's path,
#    in process), with the modules' default settings.
# Each module then publishes the step's QC progress outputs (the "qc" label the normalisation entry waits for).
protein_step <- function(operation, parameters, step_name) withDIAStateScope(manager, function() {
    entry <- workflowOperationEntry(operation)
    directory <- file.path(work, "requests", operation); dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    request <- list(operation = operation, input_references = list(state = getDIAStateReference(manager)),
        parameters = parameters, worker_cap = 12 * 1024^3, directory = directory)
    value <- get(entry$handler, mode = "function", envir = ns)(manager, request, entry$default_variant)
    value$operation <- operation
    setDIAStepProvenance(manager, operation, entry$default_variant, entry$parameters, entry$declaration$fixed, values = parameters)
    result <- completeDIAProteinAction(workflow, value, publish = TRUE)
    invisible(publishDIAProteinProgress(result, step_name, "proteomics", label, updateProteinFiltering))
})
run("accession", protein_step("protein_accession_cleanup", list(delimiter = ";", aggregation_method = "mean"), "10_protein_accession_cleaned"))
run("protein_intensity", protein_step("protein_intensity_filter", list(strict_mode = FALSE, min_reps_per_group = 2L, min_groups = 2L, proteins_intensity_cutoff_percentile = 1L), "11_protein_intensity_filtered"))
run("duplicates", protein_step("protein_duplicate_removal", list(aggregation_method = "mean"), "12_duplicates_removed"))
run("protein_replicate", protein_step("protein_replicate_filter", list(grouping_variable = "group", parallel_cores = 4L), "13_protein_replicate_filtered"))
# 4. Normalisation and RUV, correlation, as the reference entry points run them with the journey's bound settings.
norm_data <- new.env(parent = emptyenv())
norm <- journey$normalisation
if (!is.null(corpus$normalisation)) norm <- utils::modifyList(norm, corpus$normalisation)
settings <- c(list(norm_method = norm$norm_method, ruv_mode = norm$ruv_mode), norm$ruv[setdiff(names(norm$ruv), "ruv_grouping_variable")],
    list(ruv_percentage = 5L, ruv_k = NA, color_variable = NULL, shape_variable = NULL))
grouping <- function() getProtNormRuvGroupingVariable(norm$ruv$ruv_grouping_variable %||% "group")
aesthetics <- function() getProtNormPlotAesthetics(NULL, NULL)
run("normalisation_ruv", invisible(runDIAReferenceNormalization(settings, workflow, norm_data, paths, "proteomics", label,
    function(reference) NULL, function(reference) NULL, grouping)))
run("correlation", invisible(runDIAReferenceCorrelation(workflow, norm_data, paths, norm$correlation_threshold, FALSE, "proteomics", label, grouping, aesthetics)))
# 5. The native session export (no state) and the DA on the head with the contrasts the DA tab generates.
run("session_export", invisible(runProtNormExportSessionWorkflow(workflow, norm_data, list(norm_method = norm$norm_method, ruv_mode = norm$ruv_mode),
    paths$source_dir, withProgressFn = function(message, value, expr) force(expr), incProgressFn = function(...) NULL, messageFn = function(...) NULL)))
run("da", withDIAStateScope(manager, function() {
    # The contrasts the DA tab finds on the workflow: the design builder's save generated them (the TMT fixture's
    # design lists WT first, and the builder still names the contrast KO-WT; a pairwise combination of the design's
    # groups in order would not).
    stopifnot(is.data.frame(workflow$contrasts_tbl), nrow(workflow$contrasts_tbl) > 0L)
    da_data <- new.env(parent = emptyenv()); da_data$current_s4_object <- getDIAStateReference(manager)
    da_data$contrasts_available <- workflow$contrasts_tbl$contrasts
    invisible(runDIAReferenceDaAnalysis(list(formula_string = journey$da$formula, da_q_val_thresh = journey$da$q_value, treat_lfc_cutoff = journey$da$treat_lfc),
        da_data, workflow, paths, contrasts = workflow$contrasts_tbl, annotations = workflow$uniprot_dat_cln, publish = TRUE))
}))
store_dir <- harness$kitFindStore(file.path(work, label))
head <- harness$kitReadHead(store_dir, api)
# The fingerprint pass as the driver runs it after the journey: each step's state fingerprint from the immutable
# store, and the products first seen after the step appended under `<type>:<sha12>` in hash order.
fingerprints <- harness$kitFingerprintStore(store_dir, api)
stopifnot(identical(fingerprints$head_sha256, head$head_sha256))
seen <- character()
for (id in names(populations)) {
    populations[[id]]$fingerprint <- fingerprints$states[[populations[[id]]$state]]
    new_ids <- setdiff(populations[[id]]$product_ids, seen)
    seen <- union(seen, populations[[id]]$product_ids)
    products <- fingerprints$products[intersect(new_ids, names(fingerprints$products))]
    # Rendered assets (the interactive volcano page, the QC plot images) are display renderings whose bytes carry
    # widget ids and device state; the KIT-005 oracles carry none and the verifier reports a member's as additional
    # products, so the protein route's oracles leave them out too.
    products <- products[!vapply(products, function(p) identical(p$kind, "rendered_asset"), logical(1))]
    products <- unname(products[order(vapply(products, `[[`, character(1), "sha256"))])
    names(products) <- vapply(products, function(p) paste0(p$product_type %||% p$kind, ":", substr(p$sha256, 1L, 12L)), character(1))
    populations[[id]]$fingerprint$products <- c(populations[[id]]$fingerprint$products, products)
}
recorded <- list(schema = "kit.oracle_record.v1", corpus = corpus_key, corpus_sha256 = corpus$files[[1L]]$sha256, route = "protein",
    build = utils::packageDescription("MultiScholaR")$`Config/MultiScholaR/SourceCommit` %||% NA_character_, library = find.package("MultiScholaR"),
    store = store_dir, head_sha256 = head$head_sha256, history = head$history, populations = populations,
    settings = list(normalisation = norm, da = journey$da), recorded_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
closeDIAStateAdapter(manager)
saveRDS(recorded, output, compress = FALSE, version = 3)
cat("recorded", corpus_key, "->", output, "| states", length(populations), "| history", paste(head$history, collapse = ","), "|", round(recorded$seconds), "s\n")
