# KIT-098 oracle recorder (protocol v1.8 section 10 item 1, under section 7.1): an in-process run of the bound build on the
# bound transcriptomics corpus through the entry points KIT-096/097's tabs call (modelled on kit-032/oracle-record.R):
# - the session owner from the start;
# - the import tab's processing on the count table's preview;
# - the design builder (the lipid builder under trxOmicProfile()) with the driver's actions: factors in design.tsv order,
#   runs assigned by name, each biological sample's two runs as technical replicates, contrast second group vs first, save;
#   then the design tab's commit;
# - then the declared operations through the wrappers the tabs run in process: the expression filter, the QC finalise, the
#   normalisation (TMM/logCPM, grid, RUVs), the DA on the RUVs head (with the managed tab's settings: the realised formula,
#   FDR, fold-change line, trxDaContrastEntries) and the RUVIII export after DA.
# Each state is fingerprinted with the harness's canonical form (kitFingerprintStore); products attach to the step after which
# they first appear (kitStepProducts, rendered pages left out). The report and the session export commit nothing; their
# populations are the head as it stands. Every setting is read from bindings.json (the ticket key the driver reads).
# Read-only with respect to any member: writes docs/qa/workflow-execution/kit/kit-098/oracles-<label>_<corpus>.rds.
# Usage (installed bound build on R_LIBS):
#   KIT098_INSTALLED=1 R_LIBS=<library> Rscript docs/qa/workflow-execution/kit/kit-098/oracle-record.R trx98a|trx98b [output rds]
if (identical(Sys.getenv("KIT098_INSTALLED"), "1")) suppressMessages(library(MultiScholaR)) else suppressMessages(devtools::load_all(quiet = TRUE))
args <- commandArgs(TRUE)
label <- args[[1L]]
tickets <- c(trx98a = "098a", trx98b = "098b")
stopifnot(label %in% names(tickets))
repo <- normalizePath(".")
started <- Sys.time()
`%||%` <- function(a, b) if (is.null(a)) b else a
ns <- asNamespace("MultiScholaR"); api <- function(name) get(name, ns)
for (f in ls(ns, all.names = FALSE)) if (!exists(f, envir = globalenv(), inherits = TRUE) && is.function(get(f, ns))) assign(f, get(f, ns), envir = globalenv())
harness_path <- file.path(repo, "tools/profiling/workflow_kit.R")
harness <- new.env(parent = globalenv()); sys.source(harness_path, harness)
bindings_path <- file.path(repo, "docs/qa/workflow-execution/kit/bindings.json")
bindings <- jsonlite::read_json(bindings_path)
spec <- bindings$tickets[[tickets[[label]]]]
stopifnot(identical(spec$route, "transcriptomics"), identical(spec$oracles, label))
corpus_key <- spec$corpus
corpus <- bindings$corpora[[corpus_key]]
if (is.null(corpus)) stop("corpus ", corpus_key, " is not bound")
for (name in names(corpus$files)) stopifnot(identical(digest::digest(file = corpus$files[[name]]$path, algo = "sha256"), corpus$files[[name]]$sha256))
output <- if (length(args) > 1L) args[[2L]] else file.path(repo, "docs/qa/workflow-execution/kit/kit-098", sprintf("oracles-%s_%s.rds", label, corpus_key))
counts_path <- corpus$files[[corpus$assay_files$counts]]$path
design_path <- corpus$files[[corpus$design_file]]$path
settings <- spec$normalisation
driver_settings <- spec$driver_settings
project_label <- "kit_omics_lifecycle" # the driver's experiment label (runKitTrxImportActions)
work <- file.path(bindings$harness$root, "oracles", label); unlink(work, recursive = TRUE); dir.create(work, recursive = TRUE)
Sys.setenv(MULTISCHOLAR_TEST_MODE = "true")
session <- shiny::MockShinySession$new()
steps <- harness$kitJourneySteps("transcriptomics")
populations <- list()
variants <- list()
step_started <- Sys.time()
snapshot_step <- function(step) {
    store_dir <- harness$kitFindStore(file.path(work, "project"))
    stopifnot(!is.null(store_dir))
    head <- harness$kitReadHead(store_dir, api)
    state <- head$current_state
    populations[[step$id]] <<- list(step = step$id, expected_state = step$state, state = state,
        state_as_expected = is.null(step$state) || identical(state, step$state), store = store_dir, head_sha256 = head$head_sha256,
        entity_id = head$states[[state]], rows = head$rows[[state]], product_ids = as.character(unlist(head$products)),
        seconds = as.numeric(difftime(Sys.time(), step_started, units = "secs")))
    cat(sprintf("  %-18s %-26s rows %s  %5.0f s\n", step$id, state, format(head$rows[[state]] %||% NA), populations[[step$id]]$seconds))
}
run <- function(step_id, expr, operations = character()) {
    step <- steps[[match(step_id, vapply(steps, `[[`, "", "id"))]]
    step_started <<- Sys.time()
    force(expr)
    for (operation in operations) variants[[operation]] <<- tryCatch(workflowOperationEntry(operation)$default_variant %||% "whole",
        error = function(error) NA_character_)
    snapshot_step(step)
}
quiet <- function(...) invisible(NULL)
# The app sets the project up with no label (app_server.R: the experiment directory is the label).
project_dirs <- setupDirectories(base_dir = file.path(work, "project", project_label), omic_types = "transcriptomics", label = NULL,
    force = TRUE, reuse_existing = TRUE)
assign("project_dirs", project_dirs, envir = globalenv())
paths <- project_dirs[["transcriptomics"]]
workflow <- new.env(parent = emptyenv())
workflow$config_list <- list(); workflow$processing_log <- list(); workflow$state_update_trigger <- NULL
workflow$data_type <- "transcriptomics"; workflow$contrasts <- list()
workflow$tab_status <- list(setup_import = "pending", design_matrix = "disabled", quality_control = "disabled", normalization = "disabled",
    differential_analysis = "disabled", session_summary = "disabled")
# The session owner from module start (mod_transcriptomics_server), bound to the route at the import.
workflow$state_manager <- createWorkflowSessionOwner(workflow, paths, session, omic_type = "transcriptomics")
profile <- trxOmicProfile()
# 1. Import: the tab's preview (the id and name columns proposed) and its processing.
run("import", shiny::withReactiveDomain(session, {
    preview <- trxImportPreview(counts_path)
    runTrxImportProcessing(workflowData = workflow, assay1Name = profile$assay_name, assay1Data = preview$assayData, assay1File = counts_path,
        lipidIdCol = preview$updates$lipidId$selected, annotationCol = preview$updates$annotation$selected,
        sampleColumns = preview$importResult$sample_columns, vendorFormat = "trx_counts", detectedFormat = preview$detectedFormat,
        isPattern = NULL, sanitizeNames = NULL, assay2File = NULL, assay2Name = NULL, experimentPaths = paths,
        notify = quiet, removeNotify = quiet, logError = function(message) stop(message, call. = FALSE))
    stopifnot(identical(workflow$state_manager$workflow_type, "TRANSCRIPTOMICS"))
}))
# 2. Design: the builder server with the driver's actions (runKitDesignActions), then the design tab's commit.
run("design", {
    kit_design <- utils::read.delim(design_path, stringsAsFactors = FALSE, check.names = FALSE)
    kit_groups <- unique(kit_design$group)
    kit_holder <- new.env(parent = emptyenv())
    shiny::testServer(mod_lipid_design_builder_server, args = list(data_tbl = function() getWorkflowOmicsBusTables(workflow, "data_tbl"),
        config_list = shiny::reactive(workflow$config_list), column_mapping = shiny::reactive(workflow$column_mapping),
        state_owner = function() workflow$state_manager, profile = profile), {
        for (i in seq_along(kit_groups)) { session$setInputs(new_factor = kit_groups[[i]]); session$setInputs(add_factor = i) }
        for (i in seq_along(kit_groups)) {
            session$setInputs(selected_runs = kit_design$sample[kit_design$group == kit_groups[[i]]], factor1_select = kit_groups[[i]],
                factor2_select = "", factor3_select = "", replicate_start = 1L)
            session$setInputs(assign_metadata = i)
        }
        kit_tech <- 0L
        for (kit_group in kit_groups) for (kit_replicate in unique(kit_design$replicates[kit_design$group == kit_group])) {
            kit_samples <- kit_design$sample[kit_design$group == kit_group & kit_design$replicates == kit_replicate]
            if (length(kit_samples) < 2L) next
            kit_tech <- kit_tech + 1L
            session$setInputs(tech_rep_samples = kit_samples, tech_rep_assignment_mode = "manual", manual_replicate_number = as.integer(kit_replicate))
            session$setInputs(assign_tech_reps = kit_tech)
        }
        session$setInputs(contrast_group1 = kit_groups[[2L]], contrast_group2 = kit_groups[[1L]], formula_string = "~ 0 + group")
        session$setInputs(add_contrast = 1L)
        session$setInputs(save_results = 1L)
        kit_holder$results <- shiny::isolate(session$returned())
    })
    results <- kit_holder$results
    stopifnot(is.list(results), nrow(results$design_matrix) == nrow(kit_design))
    results$data_cln <- resolveWorkflowOwnedValue(results$data_cln, workflow$state_manager)
    shiny::withReactiveDomain(session, runLipidDesignBuilderObserverShell(results = results, workflowData = workflow, experimentPaths = paths,
        qcTrigger = function(value) invisible(NULL), profile = profile))
    stopifnot(identical(workflow$state_manager$current_state, "trx_raw_data_s4"))
})
owner <- workflow$state_manager
config <- function() workflow$config_list
# 3. QC: the expression filter and the finalise (the profile's collaborators' store step).
run("expression_filter", invisible(runTrxReferenceExpressionFilter(owner, config())), "trx_expression_filter")
run("qc_finalize", invisible(commitWorkflowOmicsStoreStep(owner, "trx_qc_finalize", "trx_qc_complete", config(),
    "QC processing complete - ready for normalization")))
# 4. The normalisation tab's one action, with the tab's parameters as trxNormGridParameters reads them.
grid_parameters <- trxRuvGridParameters(settings$ruv_mode, as.numeric(settings$auto_percentage_min), as.numeric(settings$auto_percentage_max),
    as.numeric(settings$auto_percentage_step), as.numeric(settings$k_penalty_weight))
stopifnot(identical(settings$ruv_mode, "automatic"))
run("normalisation", invisible(runTrxReferenceRuvs(owner, grid_parameters, config())), c("trx_tmm_logcpm", "trx_ruv_grid", "trx_ruvs"))
# 5. The DA on the RUVs head with the managed tab's settings (runWorkflowOmicsManagedDa's: the realised formula, the thresholds,
# the binding's contrast entries).
run("da", {
    source <- getWorkflowStateReference(owner)
    formula <- withWorkflowStateScope(owner, function() trxRealisedFormula(owner, source))
    da_settings <- list(formula_string = formula, da_q_val_thresh = driver_settings$da$q_value, treat_lfc_cutoff = driver_settings$da$treat_lfc)
    stopifnot(is.data.frame(workflow$contrasts_tbl), nrow(workflow$contrasts_tbl) > 0L)
    invisible(runTrxReferenceDaAnalysis(owner, source, da_settings, trxDaContrastEntries(workflow$contrasts_tbl)))
}, "trx_da_contrast")
# 6. The RUVIII export after DA (a product; the head stays the RUVs state).
run("ruviii_export", invisible(runTrxReferenceRuviiiExport(owner, driver_settings$ruviii_m_source)), "trx_ruviii_export")
run("report", invisible(NULL))
run("session_export", invisible(NULL))
store_dir <- harness$kitFindStore(file.path(work, "project"))
head <- harness$kitReadHead(store_dir, api)
fingerprints <- harness$kitFingerprintStore(store_dir, api)
stopifnot(identical(fingerprints$head_sha256, head$head_sha256))
attached <- harness$kitStepProducts(populations, fingerprints$products, rendered = FALSE)
for (id in names(populations)) {
    populations[[id]]$fingerprint <- fingerprints$states[[populations[[id]]$state]]
    populations[[id]]$fingerprint$products <- c(populations[[id]]$fingerprint$products, attached[[id]])
}
record <- buildWorkflowAnalysisRecord(owner)
check <- validateWorkflowAnalysisRecord(record)
software <- extSoftVersion()
recorded <- list(schema = "kit.oracle_record.v1", label = label, ticket = tickets[[label]], corpus = corpus_key,
    corpus_sha256 = corpus$files[[1L]]$sha256, route = "transcriptomics", protocol = "v1.8 (section 10, under section 7.1)",
    protocol_sha256 = digest::digest(file = file.path(repo, "docs/qa/workflow-execution/kit/protocol.md"), algo = "sha256"),
    bindings_sha256 = digest::digest(file = bindings_path, algo = "sha256"),
    workflow_kit_sha256 = digest::digest(file = harness_path, algo = "sha256"),
    build = utils::packageDescription("MultiScholaR")$`Config/MultiScholaR/SourceCommit` %||% NA_character_, library = find.package("MultiScholaR"),
    environment = list(r = R.version.string, blas = software[["BLAS"]], lapack = La_library(), lapack_version = La_version()),
    store = store_dir, head_sha256 = head$head_sha256, history = head$history, populations = populations, variants = variants,
    settings = list(normalisation = settings, grid_parameters = grid_parameters, driver = driver_settings, design = "design.tsv: factors in file order, runs by name, technical replicates per biological sample, contrast second vs first, formula ~ 0 + group"),
    analysis_record = list(valid = check$valid, incomplete = check$incomplete, operations = vapply(record$steps, function(s) s$operation, "")),
    recorded_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
closeDIAStateAdapter(owner)
saveRDS(recorded, output, compress = FALSE, version = 3)
cat("recorded", label, corpus_key, "->", output, "| states", length(populations), "| history", paste(head$history, collapse = ","),
    "| record valid", check$valid, "|", round(recorded$seconds), "s\n")
