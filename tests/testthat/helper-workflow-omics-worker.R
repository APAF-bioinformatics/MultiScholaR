# KIT-078: the omics worker on the declaration fixtures, in this process: a frozen request through `workflowWorkerMain`,
# as the worker's service runs it, against a project the controller's owner keeps open.

#' Run one omics action's frozen request through the worker's entry point.
#' @param fixture `lipidDeclarationsOwner` / `metabDeclarationsOwner`.
#' @param operation The registered operation.
#' @param parameters The request's parameters.
#' @param directory The request's directory; a new one under the fixture's source directory when NULL.
#' @return The worker's handoff (`worker-result.rds`), with the request's directory.
omicsWorkerRun <- function(fixture, operation, parameters = list(), directory = NULL, token = list(request_id = "wfx-test")) {
  owner <- fixture$owner
  store <- owner$.context$store
  if (is.null(directory)) directory <- tempfile("request-", tmpdir = fixture$source_dir)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  request <- list(schema = 1L, build = workflowExecutionBuildIdentity(), operation = operation,
    input_references = list(state = getWorkflowStateReference(owner)), parameters = parameters,
    token = token, owner = store$identity, project_root = dirname(store$directory),
    directory = directory, config = list(mode = "adapters"), worker_cap = 12 * 1024^3)
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
  c(readManagedOmicsResult(owner, accepted, NULL), list(acceptance = receipt$status, paths = paths))
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
