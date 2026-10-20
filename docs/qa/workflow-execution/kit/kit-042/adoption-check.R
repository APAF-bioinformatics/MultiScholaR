# KIT-042 adoption check. Two phases, two libraries:
#   control  : R_LIBS=<proof>/core/browser-smoke.h6QeDJ/library Rscript docs/qa/workflow-execution/kit/kit-042/adoption-check.R control
#              The control build (kit_036, installed package) runs the in-memory LFQ (MaxQuant) and TMT (PD) paths on the
#              e2e fixtures on a plain WorkflowState, saving every state, and writes `filtered_session_data` RDS files with
#              the in-memory collector; it runs the in-memory DA on the final state and hashes the DA tables. A third RDS
#              is the MaxQuant session flagged as the retired in-memory limpa route (`use_limpa`, `limpa_applied`).
#   candidate: R_LIBS=<proof>/core/browser-smoke.3cm61C/library Rscript docs/qa/workflow-execution/kit/kit-042/adoption-check.R candidate
#              The working tree (load_all) opens each RDS through the DA load handler, once without an owner and once with an
#              owner attached (the import of the same fixture), checks the adoption (component checkpoints, head type, legacy
#              provenance at ingest, hydrated objects identical to the RDS objects), runs DA on the adopted owner and hashes
#              the tables, exports the owner natively and reopens the archive; the limpa-flagged RDS must be refused by name.
phase <- commandArgs(TRUE)[1]
if (!phase %in% c("control", "candidate")) stop("phase: control | candidate")
if (phase == "control") suppressMessages(library(MultiScholaR)) else suppressMessages(devtools::load_all(quiet = TRUE))
started <- Sys.time()
repo <- normalizePath(".")
ns <- asNamespace("MultiScholaR")
api <- function(name) get(name, ns)
fixtures <- file.path(repo, "tests/testdata/e2e")
evidence <- file.path(repo, "docs/qa/workflow-execution/kit/kit-042")
work <- file.path(Sys.getenv("HOME"), ".cache", "multischolar-kit042-adoption-check")
cases <- list(
    maxquant = list(file = file.path(fixtures, "prot_lfq", "proteinGroups.txt"), fasta = file.path(fixtures, "prot_lfq", "proteins.fasta"),
        design = file.path(fixtures, "prot_lfq", "design.tsv"), read = function(p) importMaxQuantData(p), type = "LFQ"),
    pd_tmt = list(file = file.path(fixtures, "prot_tmt", "pd_tmt_peptides.tsv"), fasta = file.path(fixtures, "prot_tmt", "proteins.fasta"),
        design = file.path(fixtures, "prot_tmt", "design.tsv"), read = function(p) importProteomeDiscovererTMTData(p), type = "TMT"))
table_hash <- function(table) {
    table <- as.data.frame(table)
    table <- table[do.call(order, unname(as.list(table[, seq_len(min(2L, ncol(table))), drop = FALSE]))), order(names(table)), drop = FALSE]
    rownames(table) <- NULL
    digest::digest(lapply(table, function(column) if (is.numeric(column)) signif(column, 12) else as.character(column)), algo = "sha256")
}
contrasts <- data.frame(contrasts = "groupKO-groupWT", friendly_names = "KO_vs_WT", full_format = "KO_vs_WT=groupKO-groupWT", stringsAsFactors = FALSE)
da_arguments <- function(object) list(contrasts_tbl = data.frame(contrasts = contrasts$full_format, stringsAsFactors = FALSE),
    formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, eBayes_trend = TRUE, eBayes_robust = TRUE,
    args_group_pattern = "(\\d+)", args_row_id = object@protein_id_column, qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue")
project_paths <- function(base, name) {
    project_dirs <- setupDirectories(base_dir = base, omic_types = "proteomics", label = name, force = TRUE, reuse_existing = TRUE)
    assign("project_dirs", project_dirs, envir = globalenv())
    paths <- project_dirs[[1L]]
    paths$time_dir <- paths$time_dir %||% file.path(paths$publication_graphs_dir, "filtering_qc")
    dir.create(paths$time_dir, recursive = TRUE, showWarnings = FALSE)
    paths
}
config <- readConfigFile(system.file("config", "config.ini", package = "MultiScholaR"), create_workers = FALSE)

if (phase == "control") {
    unlink(work, recursive = TRUE); dir.create(work, recursive = TRUE)
    results <- list()
    for (name in names(cases)) {
        case <- cases[[name]]
        paths <- project_paths(file.path(work, name, "control-project"), name)
        workflow <- new.env(parent = emptyenv())
        workflow$state_manager <- WorkflowState$new(audit_enabled = FALSE)
        workflow$config_list <- config
        workflow$config_list$globalParameters$workflow_type <- case$type
        workflow$qc_params <- list(protein_qc = list())
        imported <- case$read(case$file)
        api("applyProtImportResultToWorkflow")(workflow, imported, name, case$fasta, logInfo = function(...) NULL, showNotification = function(...) NULL)
        design_in <- utils::read.delim(case$design, stringsAsFactors = FALSE)
        runs <- as.character(unique(workflow$data_cln[[imported$column_mapping$run_col]]))
        group <- vapply(runs, function(run) design_in$group[[which(endsWith(run, design_in$sample))[[1L]]]], "")
        design <- data.frame(Run = runs, group = unname(group), stringsAsFactors = FALSE)
        design$replicates <- stats::ave(seq_len(nrow(design)), design$group, FUN = seq_along)
        workflow$design_matrix <- design
        workflow$contrasts_tbl <- contrasts
        api("buildProtDesignStateCheckpoint")(workflow, case$type, "Import")
        manager <- workflow$state_manager
        object <- manager$getState("protein_s4_initial")
        env <- new.env(parent = emptyenv()); env$config_list <- workflow$config_list
        section <- "removeRowsWithMissingValuesPercent"
        object <- api("updateMissingValueParameters")(object, min_reps_per_group = 2L, min_groups = 2L, env = env)
        object <- api("updateConfigParameter")(object, function_name = section, parameter_name = "proteins_intensity_cutoff_percentile", new_value = 1L, env = env)
        object <- removeRowsWithMissingValuesPercent(object)
        manager$saveState("protein_intensity_filtered", object, list(strict_mode = FALSE), "Protein intensity filter")
        resolved <- api("resolveProteinDuplicateRows")(object@protein_quant_table, object@protein_id_column, base::get("mean"))
        object@protein_quant_table <- resolved$table
        manager$saveState("duplicates_removed", object, list(aggregation_method = "mean"), "Removed duplicate proteins by aggregation")
        object <- removeProteinsWithOnlyOneReplicate(object, NA, grouping_variable = "group")
        manager$saveState("protein_replicate_filtered", object, list(grouping_variable = "group"), "Applied protein replicate filter")
        object <- normaliseBetweenSamples(object, normalisation_method = "cyclicloess")
        manager$saveState("normalized", object, list(norm_method = "cyclicloess", ruv_mode = "skip", ruv_applied = FALSE), "Post-normalization complete: RUV-III skipped by user")
        pairs <- pearsonCorForSamplePairs(object, tech_rep_remove_regex = "pool", correlation_group = "group")
        object <- filterSamplesByProteinCorrelationThreshold(object, pearson_correlation_per_pair = pairs, min_pearson_correlation_threshold = 0.5)
        manager$saveState("correlation_filtered", object, list(min_pearson_correlation_threshold = 0.5), "Applied final sample correlation filter")
        norm <- new.env(parent = emptyenv()); norm$correlation_filtered_obj <- object; norm$best_k <- NA; norm$correlation_threshold <- 0.5
        session_data <- api("collectProtNormExportSessionData")(workflow, norm, list(norm_method = "cyclicloess", ruv_mode = "skip"),
            messageFn = function(...) NULL)
        out <- file.path(work, name, "rds"); dir.create(out, recursive = TRUE)
        saveRDS(session_data, file.path(out, "filtered_session_data_latest.rds"))
        # In-memory DA on the final state (the control build's contrast tables).
        arguments <- da_arguments(object)
        for (s in c("differentialAbundanceAnalysis", "differentialAbundanceAnalysisHelper")) if (is.null(object@args[[s]])) object@args[[s]] <- arguments
        result <- differentialAbundanceAnalysis(object, contrasts_tbl = arguments$contrasts_tbl, formula_string = "~ 0 + group",
            da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue")
        results[[name]] <- list(type = case$type, rds = file.path(out, "filtered_session_data_latest.rds"),
            states = manager$getHistory(), current = manager$current_state,
            state_classes = lapply(stats::setNames(nm = manager$getHistory()), function(state) class(manager$states[[state]]$data)[1]),
            state_hashes = lapply(stats::setNames(nm = manager$getHistory()), function(state) {
                value <- manager$getState(state); if (isS4(value)) table_hash(value@protein_quant_table) else NA_character_ }),
            da = list(rows = nrow(result$da_proteins_long), table_hashes = lapply(result[c("da_proteins_long", "da_proteins_wide", "contrasts_results_table")], table_hash)),
            session_fields = names(session_data), workflow_type = session_data$workflow_type, use_limpa = session_data$use_limpa)
        cat(name, "| control states", paste(manager$getHistory(), collapse = ","), "| DA rows", nrow(result$da_proteins_long), "\n")
    }
    # The retired in-memory limpa route: the MaxQuant session flagged as that route (no limpa workbook is run here;
    # the refusal keys on the flags the in-memory collector wrote for that route).
    # A session saved after enrichment (audit F-02): the MaxQuant session with an `EnrichmentResults` state appended,
    # the class the in-memory enrichment step saves.
    enriched <- readRDS(results$maxquant$rds)
    enriched$r6_complete_states$enrichment_completed <- list(timestamp = Sys.time(), data = methods::new("EnrichmentResults"),
        config = list(), description = "enrichment", audit_metadata = list(provenance_status = "legacy_or_not_applicable"))
    enriched$r6_state_history <- c(enriched$r6_state_history, "enrichment_completed")
    dir.create(file.path(work, "enriched", "rds"), recursive = TRUE)
    saveRDS(enriched, file.path(work, "enriched", "rds", "filtered_session_data_latest.rds"))
    results$enriched <- list(rds = file.path(work, "enriched", "rds", "filtered_session_data_latest.rds"), appended = "enrichment_completed")
    limpa <- readRDS(results$maxquant$rds)
    limpa$use_limpa <- TRUE; limpa$limpa_applied <- TRUE
    dir.create(file.path(work, "limpa", "rds"), recursive = TRUE)
    saveRDS(limpa, file.path(work, "limpa", "rds", "filtered_session_data_latest.rds"))
    results$limpa <- list(rds = file.path(work, "limpa", "rds", "filtered_session_data_latest.rds"), flagged = c("use_limpa", "limpa_applied"))
    # The control install predates the DESCRIPTION's SourceCommit field; the bindings record its commit.
    bindings <- jsonlite::read_json(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
    library_path <- find.package("MultiScholaR")
    bound <- Filter(function(build) identical(build$library, dirname(library_path)), bindings$builds)
    results$library <- list(path = library_path,
        source_commit = utils::packageDescription("MultiScholaR")$`Config/MultiScholaR/SourceCommit` %||%
            (if (length(bound)) bound[[1L]]$source_commit else NA_character_),
        binding = if (length(bound)) names(bound)[[1L]] else NA_character_)
    results$seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    jsonlite::write_json(results, file.path(evidence, "adoption-check-control.json"), auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = I(17))
    cat("done", round(results$seconds), "s\n")
    quit(status = 0)
}

# ---- candidate phase ----
control <- jsonlite::read_json(file.path(evidence, "adoption-check-control.json"))
Sys.setenv(MULTISCHOLAR_TEST_MODE = "true")
fake_session <- function() {
    session <- new.env(parent = emptyenv()); session$userData <- new.env(parent = emptyenv())
    session$onFlush <- session$onFlushed <- session$onSessionEnded <- function(...) function() invisible(NULL)
    session
}
load_through_handler <- function(workflow_data, da_data, paths, then) {
    # The DA tab's load handler, driven as the module contract test drives it (test-prot-07b).
    testthat::local_mocked_bindings(withProgress = function(message, value, expr) force(expr),
        incProgress = function(amount, detail = NULL) invisible(NULL),
        showNotification = function(message, type = NULL, duration = NULL, ...) invisible(NULL),
        updateTextAreaInput = function(session, inputId, value = NULL, ...) invisible(NULL),
        updateSelectInput = function(...) invisible(NULL), updateNumericInput = function(...) invisible(NULL), .package = "shiny")
    testthat::local_mocked_bindings(log_info = function(...) invisible(NULL), log_warn = function(...) invisible(NULL),
        log_error = function(...) invisible(NULL), .package = "logger")
    # The checks run inside the server block: the adopted owner is attached to the module session and closes
    # with it (attachDIAStateSession's end-of-session hook).
    result <- NULL
    shiny::testServer(function(input, output, session) {
        api("da_server_load_session_handler")(input, output, session, da_data, workflow_data, paths)
    }, {
        withCallingHandlers(session$setInputs(load_filtered_session = 1), message = function(m) invokeRestart("muffleMessage"))
        session$flushReact()
        result <<- then()
    })
    result
}
run_adoption <- function(name, case, with_owner) {
    label <- paste0(name, if (with_owner) "-owner" else "-plain")
    out <- list(format = name, type = case$type, with_owner = with_owner)
    paths <- project_paths(file.path(work, name, label, "project"), name)
    file.copy(control[[name]]$rds, file.path(paths$source_dir, "filtered_session_data_latest.rds"))
    workflow <- shiny::reactiveValues(state_manager = WorkflowState$new(audit_enabled = FALSE), config_list = config,
        tab_status = list(normalization = "pending", differential_expression = "locked"), state_update_trigger = NULL)
    da_data <- shiny::reactiveValues(current_s4_object = NULL, contrasts_available = NULL, formula_from_s4 = NULL)
    if (with_owner) {
        # An owner attached by the import of the same fixture (KIT-020) before the session is opened.
        shiny::isolate({
            workflow$config_list$globalParameters$workflow_type <- case$type
            initializeDIANativeImport(workflow, paths, fake_session(), name, "iq")
            withDIAStateScope(workflow$state_manager, function() {
                imported <- case$read(case$file)
                applyProtImportResultToWorkflow(workflow, imported, name, case$fasta, logInfo = function(...) NULL, showNotification = function(...) NULL)
            })
        })
        out$owner_before <- shiny::isolate(workflow$state_manager$.context$store$identity)
    }
    expected <- readRDS(control[[name]]$rds)
    load_through_handler(workflow, da_data, paths, function() {
    manager <- shiny::isolate(workflow$state_manager)
    out$adopted <- isDIANativeSession(manager) && usesDIAReferenceStates(manager)
    if (!out$adopted) { out$error <- "the session was not adopted"; return(out) }
    out$store_identity <- manager$.context$store$identity[c("project_id", "workflow_id")]
    out$head_type <- manager$workflow_type
    out$history <- manager$getHistory(); out$current <- manager$current_state
    snapshot <- manager$.context$snapshot
    out$states <- lapply(stats::setNames(nm = out$history), function(state) {
        # The `initial` placeholder has no data and no component record (KIT-025); the S4 states are component checkpoints.
        record <- if (identical(state, "initial")) list(schema = NA_character_, class = NA_character_) else
            readDIAStateMetadata(manager, getDIAStateReference(manager, state))$data
        restored <- manager$getState(state)
        original <- expected$r6_complete_states[[state]]$data
        provenance <- snapshot$step_provenance[[state]]
        list(schema = record$schema, class = record$class, provenance = provenance$operation, backfilled = isTRUE(provenance$backfilled),
            adopted = isTRUE(provenance$adopted),
            identical = identical(restored, original),
            hash = if (isS4(restored)) table_hash(restored@protein_quant_table) else NA_character_, control_hash = control[[name]]$state_hashes[[state]])
    })
    out$all_identical <- all(vapply(out$states, function(s) isTRUE(s$identical), logical(1)))
    # The `initial` placeholder keeps the owner's own `import` entry; every adopted scientific state is `legacy`.
    out$initial_provenance <- out$states$initial$provenance
    out$all_legacy <- all(vapply(out$states[setdiff(names(out$states), "initial")], function(s) identical(s$provenance, "legacy") && !s$backfilled && s$adopted, logical(1)))
    out$all_component <- all(vapply(out$states[setdiff(names(out$states), "initial")], function(s) identical(s$schema, "dia.component_checkpoint.v2"), logical(1)))
    out$da_current <- shiny::isolate(isSavedWorkflowReference(da_data$current_s4_object))
    # DA on the adopted owner through the reference entry point; the tables hashed as the control's were.
    out$da <- tryCatch({
        da_input <- new.env(parent = emptyenv()); da_input$current_s4_object <- getDIAStateReference(manager); da_input$contrasts_available <- contrasts$contrasts
        da_workflow <- new.env(parent = emptyenv()); da_workflow$state_manager <- manager; da_workflow$contrasts_tbl <- contrasts
        runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), da_input, da_workflow, paths,
            contrasts = contrasts, publish = FALSE)
        manifest <- resolveDIAStateValue(manager$.context, manager$.context$snapshot$da_results)
        hashes <- withDIAStateScope(manager, function() lapply(readDIAContrastOutput(manifest$results[[1L]],
            fields = c("da_proteins_long", "da_proteins_wide", "contrasts_results_table")), table_hash))
        list(ok = TRUE, rows = manifest$results[[1L]]$rows, table_hashes = hashes,
            equal_to_control = mapply(identical, hashes, control[[name]]$da$table_hashes[names(hashes)]))
    }, error = function(e) list(ok = FALSE, error = conditionMessage(e)))
    # Native export of the adopted owner and its reopen (KIT-042 item 2).
    out$export <- tryCatch({
        norm <- new.env(parent = emptyenv()); norm$best_k <- NA; norm$correlation_threshold <- 0.5
        collected <- shiny::isolate(withDIAStateScope(manager, function() collectDIANativeSessionMetadata(workflow, norm, list(norm_method = "cyclicloess", ruv_mode = "skip"))))
        output <- file.path(work, name, label, "export"); dir.create(output, recursive = TRUE)
        exported <- withDIAStateScope(manager, function() runDIANativeSessionExport(manager, collected$metadata, output, extra_references = collected$references))
        archive <- exported$publication$archive_path
        extracted <- file.path(work, name, label, "manifest"); dir.create(extracted)
        utils::untar(archive, files = "manifest.json", exdir = extracted)
        manifest <- jsonlite::read_json(file.path(extracted, "manifest.json"))
        latest <- jsonlite::read_json(exported$publication$latest_path)
        reopen_root <- file.path(work, name, label, "reopen"); dir.create(reopen_root)
        reopened <- openDIANativeSession(archive, reopen_root)
        on.exit(closeDIAStateAdapter(reopened$manager), add = TRUE)
        list(ok = TRUE, archive = basename(archive), manifest_workflow = manifest$workflow, manifest_route = manifest$route,
            latest_fields = names(latest), reopened_type = reopened$manager$workflow_type, reopened_state = reopened$manager$current_state,
            metadata_workflow_type = collected$metadata$workflow_type)
    }, error = function(e) list(ok = FALSE, error = conditionMessage(e)))
    out
    })
}
unlink(list.files(work, pattern = "-(owner|plain)$", full.names = TRUE, include.dirs = TRUE, recursive = FALSE), recursive = TRUE)
results <- list()
for (name in names(cases)) for (with_owner in c(FALSE, TRUE)) {
    key <- paste0(name, if (with_owner) "_owner" else "_plain")
    results[[key]] <- tryCatch(run_adoption(name, cases[[name]], with_owner), error = function(e) list(format = name, with_owner = with_owner, error = conditionMessage(e)))
    r <- results[[key]]
    if (!is.null(r$error)) { cat(key, "ERROR:", r$error, "\n"); next }
    cat(key, "| store", r$store_identity$project_id, r$store_identity$workflow_id, "| head", r$head_type, "| states", length(r$states),
        "component", r$all_component, "identical", r$all_identical, "legacy-at-ingest", r$all_legacy,
        "| DA", if (isTRUE(r$da$ok)) paste(r$da$rows, "rows; equal to control", paste(r$da$equal_to_control, collapse = ",")) else paste("ERR", r$da$error),
        "| export", if (isTRUE(r$export$ok)) paste(r$export$manifest_workflow, r$export$manifest_route, "reopened", r$export$reopened_type, r$export$reopened_state) else paste("ERR", r$export$error), "\n")
}
# The limpa-flagged session is refused by name, with and without an owner.
limpa_value <- readRDS(control$limpa$rds)
paths <- project_paths(file.path(work, "limpa", "project"), "limpa")
workflow <- new.env(parent = emptyenv()); workflow$state_manager <- WorkflowState$new(audit_enabled = FALSE)
refusal <- tryCatch({ ingestDIANativeLegacySession(limpa_value, workflow, fake_session(), paths$source_dir, paths); "adopted" },
    error = function(e) paste(class(e)[1], conditionMessage(e)))
results$limpa <- list(refusal = refusal, refused = grepl("^workflow_session_unsupported", refusal), owner_unchanged = !inherits(workflow$state_manager, "DIAStateAdapter"))
cat("limpa |", refusal, "\n")
# The session saved after enrichment is refused by name (audit F-02).
enriched_value <- readRDS(control$enriched$rds)
refusal <- tryCatch({ ingestDIANativeLegacySession(enriched_value, workflow, fake_session(), paths$source_dir, paths); "adopted" },
    error = function(e) paste(class(e)[1], conditionMessage(e)))
results$enriched <- list(refusal = refusal, refused = grepl("^workflow_session_unsupported", refusal) && grepl("enrichment_completed", refusal),
    owner_unchanged = !inherits(workflow$state_manager, "DIAStateAdapter"))
cat("enriched |", refusal, "\n")
results$head_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
results$tree_clean <- length(system2("git", c("status", "--porcelain", "--", "R", "tests", "DESCRIPTION", "NAMESPACE"), stdout = TRUE)) == 0L
results$control_commit <- control$library$source_commit
results$seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
jsonlite::write_json(results, file.path(evidence, "adoption-check.json"), auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = I(17))
cat("done", round(results$seconds), "s\n")
