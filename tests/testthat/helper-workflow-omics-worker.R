# KIT-078: the omics worker on the declaration fixtures, in this process: a frozen request through `workflowWorkerMain`,
# as the worker's service runs it, against a project the controller's owner keeps open.

#' Run one omics action's frozen request through the worker's entry point.
#' @param fixture `lipidDeclarationsOwner` / `metabDeclarationsOwner`.
#' @param operation The registered operation.
#' @param parameters The request's parameters.
#' @param directory The request's directory; a new one under the fixture's source directory when NULL.
#' @return The worker's handoff (`worker-result.rds`), with the request's directory.
omicsWorkerRun <- function(fixture, operation, parameters = list(), directory = NULL) {
  owner <- fixture$owner
  store <- owner$.context$store
  if (is.null(directory)) directory <- tempfile("request-", tmpdir = fixture$source_dir)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  request <- list(schema = 1L, build = workflowExecutionBuildIdentity(), operation = operation,
    input_references = list(state = getWorkflowStateReference(owner)), parameters = parameters,
    token = list(request_id = "wfx-test"), owner = store$identity, project_root = dirname(store$directory),
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
