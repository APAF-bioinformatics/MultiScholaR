# KIT-078: the omics worker on the declaration fixtures, in this process: a frozen request through `workflowWorkerMain`,
# as the worker's service runs it, against a project the controller's owner keeps open.

#' Run one omics action's frozen request through the worker's entry point.
#' @param fixture `lipidDeclarationsOwner` / `metabDeclarationsOwner`.
#' @param operation The registered operation.
#' @param parameters The request's parameters.
#' @param directory The request's directory; a new one under the fixture's source directory when NULL.
#' @return The worker's handoff (`worker-result.rds`), with the request's directory.
omicsWorkerRun <- function(fixture, operation, parameters = list(), directory = NULL, token = list(request_id = "wfx-test"), config = list(mode = "adapters")) {
  owner <- fixture$owner
  store <- owner$.context$store
  if (is.null(directory)) directory <- tempfile("request-", tmpdir = fixture$source_dir)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  request <- list(schema = 1L, build = workflowExecutionBuildIdentity(), operation = operation,
    input_references = list(state = getWorkflowStateReference(owner)), parameters = parameters,
    token = token, owner = store$identity, project_root = dirname(store$directory),
    directory = directory, config = config, worker_cap = 12 * 1024^3)
  path <- file.path(directory, "request.rds")
  writeWorkflowExecutionMessage(path, request)
  # A failed worker writes its handoff and stops; the test reads the handoff.
  withr::with_envvar(c(MULTISCHOLAR_WORKFLOW_WORKER = "1"), suppressWarnings(suppressMessages(utils::capture.output(try(
    workflowWorkerMain(c(path, digest::digest(file = path, algo = "sha256"))), silent = TRUE)))))
  result <- readRDS(file.path(directory, "worker-result.rds"))
  result$directory <- directory
  result
}

#' The binding, fixture owner and the module's parameters per tab step, one set for both omics.
omicsWorkerFixtures <- function(root) {
  list(lipidomics = list(binding = lipidomicsDeclarationBinding, fixture = lipidDeclarationsOwner(file.path(root, "lipid")),
      ruv = lipidDeclarationsRuvParameters()),
    metabolomics = list(binding = metabolomicsDeclarationBinding, fixture = metabDeclarationsOwner(file.path(root, "metab")),
      ruv = metabDeclarationsRuvParameters()))
}

#' A controller on the fixture's owner with no application (the service layer is the journeys' to exercise).
omicsWorkerController <- function(fixture) {
  owner <- fixture$owner
  store <- owner$.context$store
  controller <- WorkflowExecutionController$new()
  controller$manager <- owner
  controller$pid <- Sys.getpid()
  controller$lease <- filelock::lock(tempfile(fileext = ".lock"))
  controller$control <- list(schema = 1L, incarnation = store$identity$incarnation, generation = 1L, controller_id = "test-controller")
  store$execution_owner <- controller
  controller
}

#' One managed action on a controller, in this process: the frozen request through the worker's entry point, the
#' acceptance (the commit), the publication, and what the tab reads. Stops if the worker fails.
#' `finishWorkflowExecutionJob` (the service teardown) must be mocked by the test.
#' @param context The module's commit (`validateWorkflowOmicsCompletionContext`).
#' @return The tab's handoff (`readManagedOmicsResult`), with the acceptance's status.
omicsWorkerAction <- function(fixture, controller, operation, parameters = list(), context = NULL) {
  owner <- fixture$owner
  store <- owner$.context$store
  id <- paste0("wfx-", substr(digest::digest(list(tempfile(), operation), algo = "sha256"), 1L, 32L))
  paths <- workflowExportPaths(store$directory, id)
  token <- list(request_id = id, incarnation = store$identity$incarnation, expected_head = workflowExecutionHead(store)$entity_id,
    generation = controller$control$generation, controller_id = controller$control$controller_id)
  writeWorkflowExecutionMessage(file.path(paths$request_dir, "publication-context.rds"),
    list(experiment_label = "test", experiment_paths = fixture$paths, omic_type = workflowActionSpec(operation)$omic_type))
  writeWorkflowExportJournal(paths, list(schema = "multischolar.workflow_action.v1", request_id = id, status = "admitting",
    expected_head = token$expected_head, token_hash = digest::digest(token, algo = "sha256")))
  message <- omicsWorkerRun(fixture, operation, parameters, directory = paths$request_dir, token = token)
  if (!identical(message$status, "complete")) stop("worker failed: ", message$error)
  controller$job <- list(request = list(token = token, operation = operation, config = list(mode = "adapters"), parameters = parameters),
    paths = paths, completion_context = validateWorkflowOmicsCompletionContext(operation, context))
  receipt <- acceptWorkflowActionResult(controller, message[c("schema", "status", "token", "result")])
  accepted <- resumeWorkflowAction(controller, id)
  c(readManagedOmicsResult(owner, accepted, NULL), list(acceptance = receipt$status, paths = paths, accepted = accepted))
}

#' The whole journey of `lipidDeclarationsJourney` / `metabDeclarationsJourney` on the worker: the same values, in the shape
#' `omicsDeclarationsJourneyDump` reads.
omicsWorkerJourney <- function(fixture, binding, ruv) {
  owner <- fixture$owner
  config <- fixture$workflow$config_list
  controller <- omicsWorkerController(fixture)
  commit <- function(step) omicsWorkerAction(fixture, controller, step$operation, step$parameters,
    list(state_name = step$state_name, config = config, description = step$description))
  out <- list()
  out$intensity <- commit(omicsIntensityFilterStep(binding, 1, 0.5))
  out$duplicates <- commit(omicsDuplicateResolutionStep(binding))
  commitWorkflowOmicsStoreStep(owner, paste0(binding$prefix, "_qc_finalize"), paste0(binding$prefix, "_qc_complete"), config, "QC processing complete - ready for normalization")
  out$itsd <- commit(omicsItsdNormalisationStep(binding, "median", NULL))
  out$log2 <- commit(omicsLog2Step(binding, 1))
  out$between <- commit(omicsBetweenSampleNormalisationStep(binding, "cyclicloess"))
  out$ruv <- commit(omicsRuvStep(binding, ruv))
  pairs <- omicsWorkerAction(fixture, controller, paste0(binding$prefix, "_pearson_pairs"), list(correlation_group = "group"))
  out$correlation <- commit(omicsPearsonSelectionStep(binding, pairs$reference$entity_id, 0.5))
  out$correlation$correlations <- pairs$value
  source <- getWorkflowStateReference(owner)
  settings <- list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)
  contrasts <- fixture$workflow$contrasts_tbl
  values <- list()
  for (step in omicsDaContrastSteps(binding, settings, contrasts)) values[[step$contrast]] <- omicsWorkerAction(fixture, controller, step$operation, step$parameters)$value
  out$da <- withWorkflowStateScope(owner, function() acceptWorkflowOmicsDaManifest(binding, owner, source, settings, contrasts, values))
  out
}

#' A managed session for the tabs, in this process: `session`'s runner dispatches each action to the real worker entry point
#' and acceptance (`omicsWorkerAction`), and the probe records every assay read the controller makes between the action's
#' reservation and its delivery. The service layer (the reservation's systemd units, the poll) is mocked.
#' @param probe An environment; `tables` collects the entity ids the controller materialises during an action.
omicsManagedChainMocks <- function(fixture, probe, env = parent.frame()) {
  controller <- omicsWorkerController(fixture)
  application <- WorkflowExecutionApplication$new()
  table <- materializeWorkflowTable
  components <- readWorkflowExecutionComponents
  resume <- resumeWorkflowAction
  composite <- writeWorkflowOmicsCompositeFigure
  working <- function() isTRUE(probe$active) && !identical(Sys.getenv("MULTISCHOLAR_WORKFLOW_WORKER"), "1")
  testthat::local_mocked_bindings(
    workflowSessionApplication = function(session) application,
    workflowExecutionBuildIdentity = function() strrep("f", 64),
    finishWorkflowExecutionJob = function(owner, status) { owner$job <- NULL; list(status = status) },
    beginWorkflowGUIAction = function(application, manager, operation, publication) {
      probe$publication <- publication
      probe$active <- TRUE
      probe$operations <- c(probe$operations, operation)
      list2env(list(id = paste0("gui-", operation), closed = FALSE))
    },
    finishWorkflowGUIAction = function(reservation, status, ...) {
      probe$active <- FALSE
      probe$statuses <- c(probe$statuses, status)
      reservation$closed <- TRUE
      list()
    },
    # The composite's raster read, counted in the controller only.
    writeWorkflowOmicsCompositeFigure = function(...) {
      if (!identical(Sys.getenv("MULTISCHOLAR_WORKFLOW_WORKER"), "1")) probe$controller_png_reads <- (probe$controller_png_reads %||% 0L) + 1L
      composite(...)
    },
    runWorkflowAction = function(operation, input_references, parameters, owner, reservation, completion_context = NULL) {
      probe$handoff <- omicsWorkerAction(fixture, controller, operation, parameters, completion_context)
      paste0("worker-", operation)
    },
    pollWorkflowAction = function(owner) list(status = "complete"),
    # The step's resume returns what the in-test worker's acceptance already published.
    resumeWorkflowAction = function(owner, request_id) if (startsWith(request_id, "worker-")) probe$handoff$accepted else resume(owner, request_id),
    materializeWorkflowTable = function(context, reference, ...) {
      if (working()) probe$tables <- c(probe$tables, reference$entity_id)
      table(context, reference, ...)
    },
    readWorkflowExecutionComponents = function(manager, reference, fields = NULL, ...) {
      if (working()) probe$reads <- c(probe$reads, paste(fields %||% "<every field>", collapse = ","))
      components(manager, reference, fields = fields, ...)
    },
    .package = "MultiScholaR", .env = env)
  controller
}

#' The assay tables of every committed state in the owner, by entity id.
omicsAssayTableIds <- function(owner, binding) {
  ids <- character()
  for (name in setdiff(names(owner$.context$snapshot$states), "initial")) {
    spec <- readWorkflowStateMetadata(owner, owner$.context$snapshot$states[[name]])$data
    if (isDIAComponentCheckpoint(spec)) ids <- c(ids, vapply(spec$components[[binding$data_slot]], function(table) table$entity_id, ""))
  }
  unique(ids)
}

#' The inputs of the tabs' journey.
omicsTabInputs <- function() {
  list(intensity_cutoff_percentile = 1, proportion_below_cutoff = 0.5, apply_itsd = TRUE, itsd_aggregation = "median", log_offset = 1,
    norm_method = "cyclicloess", ruv_mode = "manual", auto_percentage_min = 1, auto_percentage_max = 20, ruv_grouping_variable = "group",
    separation_metric = "max_difference", k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, ruv_k = 1L, ruv_percentage = 20,
    min_pearson_correlation_threshold = 0.5)
}

#' One server with one omic's real tab observers: intensity, duplicates, normalisation, correlation and DA.
omicsTabServer <- function(fixture, omic, probe, managed) {
  force(fixture); force(omic); force(probe); force(managed)
  function(input, output, session) {
    workflow <- fixture$workflow
    paths <- fixture$paths
    if (managed) {
      runner <- installWorkflowOmicsActions(session, workflow, omic, "kit078", paths)
      runner$progress <- function(message) list(inc = function(...) invisible(NULL), close = function() invisible(NULL))
    }
    addLog <- function(entry) probe$log <- c(probe$log, entry)
    aesthetics <- function(...) list(color_var = "group", shape_var = NULL)
    stats <- shiny::reactiveVal(NULL); plot <- shiny::reactiveVal(NULL); info <- shiny::reactiveVal(NULL)
    daData <- list2env(list(current_state_name = NULL, analysis_complete = FALSE))
    session$userData$probe_da <- daData
    if (identical(omic, "lipidomics")) {
      normData <- createLipidNormReactiveState()
      session$userData$probe_norm <- normData
      registerLipidIntensityApplyFilterObserver(input, workflow, output, stats, plot, omic)
      registerLipidDuplicateResolveObserver(input, workflow, omic, stats, info, plot, output)
      registerLipidNormRunNormalizationObserver(input, workflow, paths, omic, normData, addLog, aesthetics, function(...) NULL)
      registerLipidNormApplyCorrelationFilterObserver(input, workflow, normData, addLog)
      shiny::observeEvent(input$run_da, {
        daData$current_state_name <- workflow$state_manager$current_state
        executeLipidDaRunAnalysis(NULL, workflow$contrasts_tbl, "~ 0 + group", 0.05, 0, daData, workflow, session, paths)
      })
    } else {
      normData <- createMetabNormReactiveState()
      session$userData$probe_norm <- normData
      mod_metab_qc_intensity_server("intensity", workflow, omic, "kit078")
      shiny::observeEvent(input$resolve_duplicates, runMetabDuplicateResolutionObserver(workflow, omic, output, info, stats, plot))
      shiny::observeEvent(input$run_normalization, runMetabNormNormalizationObserverWrapper(workflow, input, paths, omic, normData, addLog,
        getPlotAestheticsFn = aesthetics, generateCompositeFromFilesFn = function(...) NULL))
      shiny::observeEvent(input$apply_correlation_filter, runMetabNormApplyCorrelationObserverEntry(workflow, normData, input$min_pearson_correlation_threshold,
        input$ruv_grouping_variable, addLog))
      shiny::observeEvent(input$run_da, {
        daData$current_state_name <- workflow$state_manager$current_state
        runMetabDaAnalysisObserverShell(NULL, workflow$contrasts_tbl, "~ 0 + group", 0.05, 0, workflow, daData, session, paths)
      })
    }
  }
}

#' Click through the tabs' journey in the server.
omicsTabJourney <- function(session, omic, selections = NULL) {
  do.call(session$setInputs, omicsTabInputs())
  apply_intensity <- if (identical(omic, "metabolomics")) list(`intensity-intensity_cutoff_percentile` = 1, `intensity-proportion_below_cutoff` = 0.5, `intensity-apply_filter` = 1) else list(apply_filter = 1)
  do.call(session$setInputs, apply_intensity)
  session$flushReact()
  session$setInputs(resolve_duplicates = 1)
  session$flushReact()
  commitWorkflowOmicsStoreStep(session$userData$probe_owner, paste0(if (identical(omic, "metabolomics")) "metab" else "lipid", "_qc_finalize"),
    paste0(if (identical(omic, "metabolomics")) "metab" else "lipid", "_qc_complete"), list(), "QC processing complete - ready for normalization")
  if (!is.null(selections)) {
    norm <- session$userData$probe_norm
    norm$itsd_selections <- selections$rows
    norm$itsd_candidate_ids <- selections$candidates
  }
  session$userData$probe_norm$assay_names <- selections$assays
  session$setInputs(run_normalization = 1)
  session$flushReact()
  session$setInputs(apply_correlation_filter = 1)
  session$flushReact()
  session$setInputs(run_da = 1)
  session$flushReact()
  invisible(NULL)
}
