# KIT-046 route check (in process). The screen DIA-NN fixture imported on both DIA-NN routes through the entry points
# the GUI calls (the owner attached by `initializeDIANativeImport` with the chosen route, the reader, the import
# application, FASTA, configuration, input completion, the design builder's save, the reference Q-value filter and
# precursor rollup); the two routes' states compared by component hash; the manifest route by type; an old in-memory
# limpa session refused by name; the template selection by type (a DIA project with a limpa configuration section and
# limpa arguments selects the DIA-NN report).
# Usage: [KIT046_INSTALLED=1] R_LIBS=<library> Rscript docs/qa/workflow-execution/kit/kit-046/route-check.R <label>
if (identical(Sys.getenv("KIT046_INSTALLED"), "1")) suppressMessages(library(MultiScholaR)) else suppressMessages(devtools::load_all(quiet = TRUE))
started <- Sys.time()
label <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[[1L]] else "tree"
repo <- normalizePath(".")
`%||%` <- function(a, b) if (is.null(a)) b else a
ns <- asNamespace("MultiScholaR"); api <- function(name) get(name, ns)
for (f in ls(ns, all.names = FALSE)) if (!exists(f, envir = globalenv(), inherits = TRUE) && is.function(get(f, ns))) assign(f, get(f, ns), envir = globalenv())
harness <- new.env(parent = globalenv()); sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), harness)
bindings <- jsonlite::read_json(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
corpus <- bindings$corpora$screen
report <- corpus$files[["report.tsv"]]$path; fasta <- corpus$files[["proteins.fasta"]]$path; design_path <- corpus$files[["design.tsv"]]$path
work <- file.path(bindings$harness$root, "kit046", label); unlink(work, recursive = TRUE); dir.create(work, recursive = TRUE)
Sys.setenv(MULTISCHOLAR_TEST_MODE = "true")
fake_session <- function() {
    session <- new.env(parent = emptyenv()); session$userData <- new.env(parent = emptyenv())
    session$onFlush <- session$onFlushed <- session$onSessionEnded <- function(...) function() invisible(NULL)
    session
}
out <- list(schema = "kit.route_check.v1", ticket = "KIT-046", label = label, date = format(Sys.Date()),
    head_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    library = list(path = find.package("MultiScholaR"), source_commit = utils::packageDescription("MultiScholaR")$`Config/MultiScholaR/SourceCommit` %||% NA_character_),
    limpa_installed = requireNamespace("limpa", quietly = TRUE))
run_route <- function(route) {
    res <- list(route = route)
    label_dir <- "dia_baseline_smoke"
    project_dirs <- setupDirectories(base_dir = file.path(work, route, label_dir), omic_types = "proteomics", label = label_dir, force = TRUE, reuse_existing = TRUE)
    assign("project_dirs", project_dirs, envir = globalenv())
    paths <- project_dirs[[1L]]
    paths$time_dir <- paths$time_dir %||% file.path(paths$publication_graphs_dir, "filtering_qc")
    dir.create(paths$time_dir, recursive = TRUE, showWarnings = FALSE)
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- WorkflowState$new(audit_enabled = FALSE)
    workflow$tab_status <- list(); workflow$processing_log <- list()
    session <- fake_session()
    initializeDIANativeImport(workflow, paths, session, "diann", route)
    manager <- workflow$state_manager
    res$owner_attached <- isDIANativeSession(manager)
    res$store_identity <- if (res$owner_attached) manager$.context$store$identity[c("project_id", "workflow_id")] else NULL
    res$manager_type_after_init <- manager$workflow_type
    withDIAStateScope(manager, function() {
        inputs <- c(report, fasta); names(inputs) <- c("diann_report", "fasta")
        registerDIAInputFiles(manager, inputs)
        imported <- importDIANNData(report)
        res$applied_type <<- applyProtImportResultToWorkflow(workflow, imported, "diann", fasta, sanitizeNames = TRUE,
            logInfo = function(...) NULL, showNotification = function(...) NULL)
        releaseDIAWorkingReference(manager, workflow$data_tbl)
        processProtImportFastaData(workflowData = workflow, fastaPath = fasta, organismName = "Homo sapiens", experimentPaths = paths,
            logInfo = function(...) NULL, logWarn = function(...) NULL)
        loadProtImportConfigurationResources(workflowData = workflow, experimentPaths = paths, useShinyFiles = FALSE, readConfig = readProtGuiConfig,
            updateStatus = function(...) NULL, messageFn = function(...) NULL, showNotification = function(...) NULL,
            logInfo = function(...) NULL, logError = function(...) NULL)
        finalizeProtImportSetupState(workflow, getDIAImportMetadata(manager), "diann", basename(report), basename(fasta), 9606, "Homo sapiens")
        completeDIAInputFiles(manager)
    })
    res$manager_type_after_import <- manager$workflow_type
    res$config_workflow_type <- workflow$config_list$globalParameters$workflow_type
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
        res$design_state <<- buildProtDesignStateCheckpoint(workflowData = workflow, workflowType = manager$workflow_type, actionLabel = "Save Design")
        q <- workflow$config_list$srlQvalueProteotypicPeptideClean
        # Each peptide step publishes its QC outputs as the module does (publishDIAQvalueOutput); the next step's
        # export-ready assertion requires it.
        publish <- function(result, step_name) publishDIAQvalueOutput(manager, result, list(
            data = result$filteringSummary %||% result$filteredS4@peptide_data, step_name = step_name, omic_type = "proteomics",
            experiment_label = label_dir, return_grid = TRUE, overwrite = TRUE, project_dirs = project_dirs), updateProteinFiltering)
        qv <- runDIAReferenceQvalue(workflow, as.numeric(q$qvalue_threshold %||% 0.01), as.numeric(q$global_qvalue_threshold %||% 0.01),
            as.logical(as.integer(q$choose_only_proteotypic_peptide %||% 1L)), as.numeric(q$global_pg_qvalue_threshold %||% 0.01),
            q$confidence_provenance_mode %||% "diann_main_report", Sys.time)
        invisible(publish(qv, "2_qval_filtered"))
        ru <- runDIAReferencePrecursorRollup(workflow)
        invisible(publish(ru, "4_precursor_rollup"))
    })
    store_dir <- harness$kitFindStore(file.path(work, route))
    head <- harness$kitReadHead(store_dir, api)
    fp <- harness$kitFingerprintStore(store_dir, api, states = c("raw_data_s4", "qvalue_filtered", "precursor_rollup"))
    res$history <- head$history; res$rows <- head$rows[c("raw_data_s4", "qvalue_filtered", "precursor_rollup")]
    res$components <- lapply(fp$states, function(s) lapply(s$components, `[[`, "sha256"))
    res$summary_sha256 <- lapply(fp$states, `[[`, "summary_sha256")
    res$manifest_route <- workflowRouteForType(manager$workflow_type)
    res$export_use_limpa <- tryCatch(describeDIAProteinState, error = function(e) NULL)
    closeDIAStateAdapter(manager)
    res
}
for (route in c("DIA", "DIA_LIMPA")) {
    out[[route]] <- tryCatch(run_route(route), error = function(e) list(route = route, error = conditionMessage(e)))
    r <- out[[route]]
    cat(route, "|", if (!is.null(r$error)) paste("ERROR:", r$error) else paste("owner", r$owner_attached, "| store", paste(unlist(r$store_identity), collapse = "/"),
        "| type", r$manager_type_after_init, "->", r$applied_type, "->", r$manager_type_after_import, "| config", r$config_workflow_type,
        "| design", r$design_state, "| history", paste(r$history, collapse = ","), "| manifest route", r$manifest_route), "\n")
}
if (is.null(out$DIA$error) && is.null(out$DIA_LIMPA$error)) {
    same <- function(state) identical(out$DIA$components[[state]], out$DIA_LIMPA$components[[state]])
    out$routes_equal <- list(raw_data_s4 = same("raw_data_s4"), qvalue_filtered = same("qvalue_filtered"), precursor_rollup = same("precursor_rollup"),
        rows_equal = identical(out$DIA$rows, out$DIA_LIMPA$rows))
    cat("ROUTES EQUAL (components by state):", paste(names(out$routes_equal), unlist(out$routes_equal), collapse = " "), "\n")
}
# An old in-memory limpa session is refused by name; the template follows the type.
workflow <- new.env(parent = emptyenv()); workflow$state_manager <- WorkflowState$new(audit_enabled = FALSE)
legacy <- list(workflow_type = "DIA", use_limpa = TRUE, r6_complete_states = list(), r6_state_history = character(), r6_current_state_name = "x")
refusal <- tryCatch(ingestDIANativeLegacySession(legacy, workflow, fake_session(), work), workflow_session_unsupported = function(e) e)
out$legacy_limpa_refused <- inherits(refusal, "workflow_session_unsupported")
out$legacy_limpa_message <- if (inherits(refusal, "condition")) conditionMessage(refusal) else NA_character_
dia_with_limpa_config <- new.env(parent = emptyenv())
dia_with_limpa_config$config_list <- list(globalParameters = list(workflow_type = "DIA", use_limpa = TRUE), proteinMissingValueImputationLimpa = list(engine = "fastlimpadpc"))
out$template_dia_with_limpa_section <- resolveProtSummaryReportTemplate(dia_with_limpa_config, catFn = function(...) NULL)$templateFilename
limpa_type <- new.env(parent = emptyenv()); limpa_type$config_list <- list(globalParameters = list(workflow_type = "DIA_LIMPA"))
out$template_dia_limpa_type <- resolveProtSummaryReportTemplate(limpa_type, catFn = function(...) NULL)$templateFilename
out$route_map <- as.list(workflowRoutesByType)
out$checkpoint_route_refuses_type <- inherits(tryCatch(validateDIAPeptideCheckpointRoute(list(workflow_type = "DIA_LIMPA")), error = function(e) e), "saved_entity_unsupported")
out$import_ui <- list(dia_route_choices = regmatches(readLines(file.path(repo, "R/mod_prot_import_ui.R")), regexpr('"[^"]*" = "DIA(_LIMPA)?"', readLines(file.path(repo, "R/mod_prot_import_ui.R")))),
    quantification_radio_sites = length(grep("dia_quantification_method", list.files(file.path(repo, "R"), full.names = TRUE) |> lapply(readLines) |> unlist(), fixed = TRUE)))
cat("LEGACY limpa refused:", out$legacy_limpa_refused, "|", out$legacy_limpa_message, "\n")
cat("TEMPLATE: DIA with limpa section ->", out$template_dia_with_limpa_section, "| DIA_LIMPA type ->", out$template_dia_limpa_type, "\n")
cat("ROUTE MAP:", paste(names(out$route_map), unlist(out$route_map), collapse = " "), "| checkpoint route refuses DIA_LIMPA:", out$checkpoint_route_refuses_type, "\n")
cat("IMPORT UI choices:", paste(out$import_ui$dia_route_choices, collapse = "; "), "| dia_quantification_method sites in R/:", out$import_ui$quantification_radio_sites, "\n")
out$seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
out$DIA$export_use_limpa <- NULL; out$DIA_LIMPA$export_use_limpa <- NULL
jsonlite::write_json(out, file.path(repo, "docs/qa/workflow-execution/kit/kit-046", paste0("route-check-", label, ".json")), auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
cat("TOTAL", round(out$seconds), "s\n")
