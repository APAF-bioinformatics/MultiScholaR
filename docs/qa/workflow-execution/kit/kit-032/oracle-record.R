# KIT-032 oracle recorder (protocol v1.5 item 1) for the metabolomics route (the KIT-031 recorder's twin): an in-process run of the bound build on the
# bound metab screen corpus through the entry points the GUI's metab module calls (the session owner from the start,
# the import tab's processing with the vendor reader's detected columns, the design builder's own assignment helpers
# and save shell, then the declared operations through the wrappers the tabs call: the intensity filter, the duplicate
# resolution, the QC finalise, the normalisation action's four kernels (ITSD, log2, between-sample, RUV-III in the
# bound mode), the correlation filter, the DA, and the summary tab's publication copy and report render), each state
# fingerprinted right after its step with the harness's canonical form (kitFingerprintStore). The tab inputs the
# driver leaves at their defaults are the UI defaults here (intensity 10 % / 0.5, ITSD median, log2 offset 1, the RUV
# sweep 1-20 %, max_difference, 0.5, adaptive, manual k 2 / 5 %, the grouping `group`); the bound normalisation
# settings are the lifecycle members' (bindings `lifecycle_omics$normalisation`, overridden by ticket 031l's binding).
# The owned session export commits no state and is not run; its population is the head unchanged.
# Read-only with respect to the member: it writes docs/qa/workflow-execution/kit/kit-032/oracles-metab32_<corpus>.rds,
# bound as metab32_<step>.
# Usage (installed bound build on R_LIBS):
#   KIT032_INSTALLED=1 R_LIBS=<library> Rscript docs/qa/workflow-execution/kit/kit-032/oracle-record.R metab_screen [output rds]
if (identical(Sys.getenv("KIT032_INSTALLED"), "1")) suppressMessages(library(MultiScholaR)) else suppressMessages(devtools::load_all(quiet = TRUE))
args <- commandArgs(TRUE)
corpus_key <- if (length(args)) args[[1L]] else "metab_screen"
repo <- normalizePath(".")
output <- if (length(args) > 1L) args[[2L]] else file.path(repo, "docs/qa/workflow-execution/kit/kit-032", paste0("oracles-metab32_", corpus_key, ".rds"))
started <- Sys.time()
`%||%` <- function(a, b) if (is.null(a)) b else a
ns <- asNamespace("MultiScholaR"); api <- function(name) get(name, ns)
for (f in ls(ns, all.names = FALSE)) if (!exists(f, envir = globalenv(), inherits = TRUE) && is.function(get(f, ns))) assign(f, get(f, ns), envir = globalenv())
harness <- new.env(parent = globalenv()); sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), harness)
bindings <- jsonlite::read_json(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
corpus <- bindings$corpora[[corpus_key]]
if (is.null(corpus)) stop("corpus ", corpus_key, " is not bound")
stopifnot(identical(corpus$route, "metabolomics"))
for (name in names(corpus$files)) stopifnot(identical(digest::digest(file = corpus$files[[name]]$path, algo = "sha256"), corpus$files[[name]]$sha256))
assay_files <- corpus$assay_files
assay_paths <- vapply(assay_files, function(file) corpus$files[[file]]$path, character(1))
design_path <- corpus$files[[corpus$design_file %||% "design.tsv"]]$path
journey <- bindings$journey
settings <- bindings$lifecycle_omics$normalisation
ticket_norm <- if (identical(bindings$tickets[["032m"]]$corpus, corpus_key)) bindings$tickets[["032m"]]$normalisation
if (!is.null(ticket_norm)) settings <- utils::modifyList(settings, ticket_norm)
stopifnot(settings$ruv_mode %in% c("automatic", "manual", "skip"))
project_label <- "kit_omics_lifecycle" # the driver's experiment label for the omics members (runKitOmicsImportActions)
summary_label <- journey$experiment_label # the summary tab's label the driver types (the report's file name)
work <- file.path(bindings$harness$root, "oracles", paste0("metab32_", corpus_key)); unlink(work, recursive = TRUE); dir.create(work, recursive = TRUE)
Sys.setenv(MULTISCHOLAR_TEST_MODE = "true")
session <- shiny::MockShinySession$new()
steps <- harness$kitJourneySteps("metabolomics", settings$ruv_mode)
populations <- list()
step_started <- Sys.time()
snapshot_step <- function(step) {
    store_dir <- harness$kitFindStore(file.path(work, "project"))
    stopifnot(!is.null(store_dir))
    head <- harness$kitReadHead(store_dir, api)
    state <- head$current_state
    populations[[step$id]] <<- list(step = step$id, expected_state = step$state, state = state, state_as_expected = is.null(step$state) || identical(state, step$state),
        store = store_dir, head_sha256 = head$head_sha256, entity_id = head$states[[state]], rows = head$rows[[state]],
        product_ids = as.character(unlist(head$products)), seconds = as.numeric(difftime(Sys.time(), step_started, units = "secs")))
    cat(sprintf("  %-16s %-28s rows %s  %5.0f s\n", step$id, state, format(head$rows[[state]] %||% NA), populations[[step$id]]$seconds))
}
run <- function(step_id, expr) {
    step <- steps[[match(step_id, vapply(steps, `[[`, "", "id"))]]
    step_started <<- Sys.time()
    force(expr)
    snapshot_step(step)
}
# The app sets the project up with no label (app_server.R: the experiment directory is the label); the omic keys the list.
project_dirs <- setupDirectories(base_dir = file.path(work, "project", project_label), omic_types = "metabolomics", label = NULL, force = TRUE, reuse_existing = TRUE)
assign("project_dirs", project_dirs, envir = globalenv())
paths <- project_dirs[["metabolomics"]]
workflow <- new.env(parent = emptyenv())
workflow$config_list <- list(); workflow$tab_status <- list(setup_import = "pending", design_matrix = "disabled", quality_control = "disabled", normalization = "disabled")
workflow$processing_log <- list(); workflow$state_update_trigger <- NULL
# The session owner from module start (mod_metabolomics.R), pending until the design step binds it.
workflow$state_manager <- createDIASessionOwner(workflow, paths, session, omic_type = "metabolomics")
quiet <- function(...) NULL
# 1. Import, as the tab's process_import runs it: the first assay's selection (the vendor reader's detected id and
# annotation columns and standard pattern, the tab's inputs updated from it), the sample columns resolved as the tab's
# accessor does (no pattern), names not sanitised (the driver's setting), the second assay read by the processing.
run("import", shiny::withReactiveDomain(session, {
    selection <- prepareMetabImportAssaySelectionState(assay1File = assay_paths[[1L]])
    assay1_data <- selection$importResult$data
    is_pattern <- selection$isPattern %||% ""
    result <- runMetabImportProcessing(assay1Data = assay1_data, assay1Name = names(assay_files)[[1L]], assay1File = assay_paths[[1L]],
        assay2File = if (length(assay_paths) > 1L) assay_paths[[2L]], assay2Name = if (length(assay_files) > 1L) names(assay_files)[[2L]],
        vendorFormat = corpus$vendor, detectedFormat = selection$formatInfo$format, sanitizeNames = FALSE, isPattern = is_pattern,
        getMetaboliteIdColFn = function() selection$selectedMetaboliteId, getAnnotationColFn = function() selection$selectedAnnotation,
        getSampleColumnsFn = function() resolveMetabImportSampleColumns(assay1_data, corpus$vendor, "", selection$importResult),
        workflowData = workflow, experimentPaths = paths, showNotificationFn = quiet, logInfoFn = quiet)
    stopifnot(!is.null(workflow$data_tbl), identical(workflow$state_manager$workflow_type, "METABOLOMICS"))
}))
# 2. Design, as the builder does it: the metab builder server itself (testServer) with the driver's actions (the factors,
# each group's samples assigned with factor 1 = the group and the replicate start the tab renders (1), the technical
# replicates in manual mode per corpus replicate value, one contrast (the second group minus the first), the builder's
# default formula, save), then the module's results flow (files, the S4, the owner bound, the design provenance and
# inputs, the `metab_raw_data_s4` commit).
run("design", {
    kit_design <- utils::read.delim(design_path, stringsAsFactors = FALSE, check.names = FALSE)
    kit_groups <- unique(kit_design$group)
    stopifnot(length(kit_groups) >= 2L)
    kit_holder <- new.env(parent = emptyenv())
    shiny::testServer(mod_metab_design_builder_server, args = list(data_tbl = shiny::reactive(workflow$data_tbl),
        config_list = shiny::reactive(workflow$config_list), column_mapping = shiny::reactive(workflow$column_mapping)), {
        for (kit_group in kit_groups) { session$setInputs(new_factor = kit_group); session$setInputs(add_factor = match(kit_group, kit_groups)) }
        kit_clicks <- 0L
        for (kit_group in kit_groups) {
            kit_clicks <- kit_clicks + 1L
            session$setInputs(selected_runs = kit_design$sample[kit_design$group == kit_group], factor1_select = kit_group,
                factor2_select = "", factor3_select = "", replicate_start = 1L)
            session$setInputs(assign_metadata = kit_clicks)
        }
        kit_tech <- 0L
        if ("replicates" %in% names(kit_design)) for (kit_group in kit_groups) for (kit_replicate in unique(kit_design$replicates[kit_design$group == kit_group])) {
            kit_samples <- kit_design$sample[kit_design$group == kit_group & kit_design$replicates == kit_replicate]
            if (length(kit_samples) < 2L) next
            kit_tech <- kit_tech + 1L
            session$setInputs(tech_rep_samples = kit_samples, tech_rep_assignment_mode = "manual", manual_replicate_number = as.integer(kit_replicate))
            session$setInputs(assign_tech_reps = kit_tech)
        }
        session$setInputs(contrast_group1 = kit_groups[[2L]], contrast_group2 = kit_groups[[1L]], formula_string = "~ 0 + group")
        session$setInputs(add_contrast = 1L)
        session$setInputs(save_results = 1L)
        kit_holder$results <- session$returned()
    })
    results <- kit_holder$results
    stopifnot(is.list(results), nrow(results$design_matrix) == nrow(kit_design))
    shiny::withReactiveDomain(session, registerMetabDesignBuilderResultsObserver(builderResultsRv = function() results,
        workflowData = workflow, experimentPaths = paths, qcTrigger = NULL,
        observeEventFn = function(eventExpr, handlerExpr, ...) { force(eventExpr); force(handlerExpr) }))
    stopifnot(identical(workflow$state_manager$current_state, "metab_raw_data_s4"))
})
owner <- workflow$state_manager
config <- function() workflow$config_list
# 3. The QC tab's actions through the declared operations (the intensity filter at its defaults, the duplicates, the finalise).
# Shiny delivers a whole-number numeric input as an integer; the kernels record it as such (the state summary hashes the type).
run("intensity", invisible(runMetabReferenceIntensityFilter(owner, 10L, 0.5, config(), hydrate = FALSE)))
run("duplicates", invisible(runMetabReferenceDuplicateResolution(owner, config(), hydrate = FALSE)))
run("qc_finalize", invisible(commitWorkflowOmicsStoreStep(owner, "metab_qc_finalize", "metab_qc_complete", config(), "QC processing complete - ready for normalization")))
# 4. The normalisation tab's one action (the module's order and state names): ITSD when bound, log2, between-sample, RUV-III.
ruv <- list(ruv_mode = settings$ruv_mode, ruv_grouping_variable = "group", percentage_min = 1L, percentage_max = 20L, separation_metric = "max_difference",
    k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, manual_k = 2L, manual_percentage = 5L)
for (name in names(journey$normalisation$ruv %||% list())) if (name %in% names(ruv)) ruv[[name]] <- journey$normalisation$ruv[[name]]
# The tab's internal-standard tables preselect the candidate rows (the metab normalisation tab's selection tracking records the
# selected rows; the run maps them to feature ids through the same selection table), so the kernel runs in its
# manual mode with the candidates' ids, as the GUI member does.
itsd_feature_ids <- function() {
    current <- owner$getState()
    ids <- lapply(stats::setNames(names(current@metabolite_data), names(current@metabolite_data)), function(assay) {
        table <- buildItsdSelectionTable(assay_data = current@metabolite_data[[assay]], metabolite_id_col = current@metabolite_id_column,
            annotation_cols = current@annotation_id_column)
        selected <- which(table$is_candidate)
        if (length(selected)) table$feature_id[selected] else NULL
    })
    ids <- Filter(Negate(is.null), ids)
    if (length(ids)) ids else NULL
}
run("normalisation", {
    if (isTRUE(settings$apply_itsd)) invisible(runMetabReferenceItsdNormalisation(owner, "median", itsd_feature_ids(), config(), hydrate = FALSE))
    invisible(runMetabReferenceLog2(owner, 1L, config(), hydrate = FALSE))
    invisible(runMetabReferenceBetweenSampleNormalisation(owner, settings$norm_method %||% "cyclicloess", config(), hydrate = FALSE))
    if (!identical(settings$ruv_mode, "skip")) invisible(runMetabReferenceRuv(owner, ruv, config(), hydrate = FALSE))
})
# 5. The correlation filter on the head, the session export (no state; the head unchanged), the DA on the head.
run("correlation", invisible(runMetabReferenceCorrelation(owner, ruv$ruv_grouping_variable, settings$correlation_threshold %||% 0.5, config(), hydrate = FALSE)))
run("session_export", invisible(NULL))
da_settings <- list(formula_string = journey$da$formula, da_q_val_thresh = journey$da$q_value, treat_lfc_cutoff = journey$da$treat_lfc)
run("da", {
    stopifnot(is.data.frame(workflow$contrasts_tbl), nrow(workflow$contrasts_tbl) > 0L)
    invisible(runMetabReferenceDaAnalysis(owner, getWorkflowStateReference(owner), da_settings, workflow$contrasts_tbl))
})
# 6. The summary tab: the publication copy (the fallback study parameters the copy writes when none were saved, as
# the member's driver runs it), the shipped template under scripts/, the render. The head is unchanged.
report_file <- NULL
run("report", {
    params_file <- file.path(paths$source_dir, "study_parameters.txt")
    if (!file.exists(params_file)) writeLines(c("Study Parameters", "================", "", paste("Workflow Name:", summary_label),
        "Description: KIT-068 lifecycle member", paste("Timestamp:", Sys.time()),
        "Note: Some parameters could not be saved due to serialization issues"), params_file)
    template <- file.path(paths$base_dir, "scripts", "metabolomics", "metabolomics_report.rmd")
    dir.create(dirname(template), recursive = TRUE, showWarnings = FALSE)
    if (!file.exists(template)) stopifnot(file.copy(system.file("reports", "metabolomics", "metabolomics_report.rmd", package = "MultiScholaR"), template))
    copyToResultsSummary(omic_type = "metabolomics", experiment_label = summary_label, contrasts_tbl = workflow$contrasts_tbl,
        design_matrix = workflow$design_matrix, force = TRUE)
    rendered <- RenderReport(omic_type = "metabolomics", experiment_label = summary_label, rmd_filename = "metabolomics_report.rmd")
    stopifnot(is.character(rendered), length(rendered) == 1L, file.exists(rendered))
    report_file <<- list(template = "metabolomics_report.rmd", path = sub(paste0("^", normalizePath(paths$base_dir), "/"), "", normalizePath(rendered)), bytes = file.size(rendered))
    cat(sprintf("  %-16s %-28s %s bytes\n", "report", "metabolomics_report.rmd", format(report_file$bytes)))
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
manifest <- withWorkflowStateScope(owner, function() resolveWorkflowStateValue(owner$.context, owner$.context$snapshot$da_results))
record <- buildWorkflowAnalysisRecord(owner)
check <- validateWorkflowAnalysisRecord(record)
recorded <- list(schema = "kit.oracle_record.v1", corpus = corpus_key, corpus_sha256 = corpus$files[[1L]]$sha256, route = "metabolomics",
    build = utils::packageDescription("MultiScholaR")$`Config/MultiScholaR/SourceCommit` %||% NA_character_, library = find.package("MultiScholaR"),
    store = store_dir, head_sha256 = head$head_sha256, history = head$history, populations = populations,
    settings = list(normalisation = settings, ruv = ruv, intensity = list(percentile = 10L, proportion = 0.5), itsd_aggregation = "median", itsd_selection = "the tab's preselected candidates per assay",
        log_offset = 1L, da = journey$da, session_export = "not run: the owned export commits no state"),
    da_model = manifest$model, analysis_record = list(valid = check$valid, incomplete = check$incomplete,
        operations = vapply(record$steps, function(s) s$operation, "")),
    report = report_file, recorded_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
closeDIAStateAdapter(owner)
saveRDS(recorded, output, compress = FALSE, version = 3)
cat("recorded metab32", corpus_key, "->", output, "| states", length(populations), "| history", paste(head$history, collapse = ","),
    "| record valid", check$valid, "| legacy steps", sum(recorded$analysis_record$operations == "legacy"), "|", round(recorded$seconds), "s\n")
