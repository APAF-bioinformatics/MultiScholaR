# KIT-071: adopt the kit_065 legacy filtered-session fixtures (written by the in-memory modules on the InuLC5 build
# from the e2e lanes; kept under the proof root with their hashes in the outcome) through the DE tab's Open on a
# fresh owner, and report the adoption: the owner replaced, provenance `adopted`, the state and class, the head's
# component ids equal to a fresh checkpoint of the fixture's object (content-addressed).
# usage: Rscript adopt-legacy.R tree|installed <label>   (tree: load_all of the working tree; installed: library())
args <- commandArgs(TRUE)
mode <- args[[1]]; label <- args[[2]]
legacy_dir <- "/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/kit/kit-071-legacy"
if (identical(mode, "tree")) {
  suppressMessages(pkgload::load_all("/home/doktersmol/Projects/APAF/MultiScholaR", quiet = TRUE, export_all = TRUE))
} else {
  library(MultiScholaR)
}
f <- function(name) get(name, envir = asNamespace("MultiScholaR"))
session <- new.env(parent = emptyenv())
session$userData <- new.env(parent = emptyenv())
session$onFlush <- session$onFlushed <- session$onSessionEnded <- function(...) function() invisible(NULL)
quiet <- function(...) invisible(NULL)
ids <- function(components) sort(unname(vapply(components, function(element) element$entity_id, character(1))))
cases <- list(
  lipidomics = list(type = "LIPIDOMICS", slot = "lipid_data", file = "lipid_filtered_session_data_latest.rds"),
  metabolomics = list(type = "METABOLOMICS", slot = "metabolite_data", file = "metab_filtered_session_data_latest.rds"))
for (omic in names(cases)) {
  case <- cases[[omic]]
  fixture <- file.path(legacy_dir, case$file)
  root <- tempfile("kit-071-adopt-"); dir.create(root)
  source_dir <- file.path(root, "source"); dir.create(source_dir)
  file.copy(fixture, file.path(source_dir, case$file))
  workflow <- list2env(list(state_update_trigger = NULL))
  workflow$state_manager <- f("createDIASessionOwner")(workflow, list(source_dir = source_dir), session, omic_type = omic)
  shiny::withReactiveDomain(session, f("bindWorkflowOwnerTypeAtDesign")(workflow, case$type))
  previous <- workflow$state_manager
  da <- new.env(parent = emptyenv())
  path <- file.path(source_dir, case$file)
  result <- shiny::withReactiveDomain(session, if (identical(omic, "lipidomics")) {
    f("loadLipidDaSessionFromFile")(path, workflow, da, session, showLoadingFn = function() NULL,
      restorePostReadFn = function(...) NULL, finalizeSuccessFn = function(...) NULL,
      finalizeErrorFn = function(errorMessage, ...) list(errorMessage = errorMessage))
  } else {
    f("runMetabDaLoadSessionObserverShell")(path, workflow, da, session, showNotification = quiet,
      removeNotification = quiet, logInfo = quiet, logError = quiet,
      restoreState = function(sessionData, sessionFile, workflowData, daData, session, debugLog) {
        f("restoreMetabDaLoadedSessionState")(sessionData, sessionFile, workflowData, daData, session,
          updateSelectInput = quiet, updateTextAreaInput = quiet, logInfo = quiet, logWarn = quiet)
      })
  })
  adopted <- workflow$state_manager
  value <- readRDS(fixture)
  record <- f("readWorkflowStateMetadata")(adopted)
  checkpoint <- f("checkpointWorkflowExecutionResult")(adopted, value$current_s4_object)
  hydrated <- adopted$getState()
  provenance <- adopted$.context$snapshot$step_provenance[[adopted$current_state]]
  rows_fixture <- vapply(methods::slot(value$current_s4_object, case$slot), nrow, integer(1))
  rows_hydrated <- vapply(methods::slot(hydrated, case$slot), nrow, integer(1))
  cat(sprintf(paste0("%s: fixture=%s sha256=%s status=%s owner_replaced=%s previous_closed=%s identity=%s/%s omic=%s ",
      "state=%s da_state=%s class=%s provenance=%s adopted=%s head_ids_equal_fixture_checkpoint=%s rows_fixture=%s rows_hydrated=%s build=%s\n"),
    omic, case$file, digest::digest(file = fixture, algo = "sha256"), result$status, !identical(adopted, previous),
    isTRUE(previous$.context$store$closed), adopted$.context$store$identity$project_id, adopted$.context$store$identity$workflow_id,
    adopted$.context$omic_type, adopted$current_state, da$current_state_name, record$data$class, provenance$operation,
    isTRUE(provenance$adopted), identical(ids(record$data$components[[case$slot]]), ids(checkpoint$components[[case$slot]])),
    paste(names(rows_fixture), rows_fixture, sep = "=", collapse = ","), paste(names(rows_hydrated), rows_hydrated, sep = "=", collapse = ","),
    f("workflowSoftwareIdentity")()$source_commit))
  f("closeDIAStateAdapter")(adopted)
  unlink(root, recursive = TRUE)
}
cat(sprintf("ADOPT-DONE %s %s\n", mode, label))
