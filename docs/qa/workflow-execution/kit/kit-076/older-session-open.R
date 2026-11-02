# KIT-076 acceptance: an older session on the new build. A copy of a kit_075 journey member's project (its head holds a
# v1 manifest): the record, the DE tab's Open on the member's export, the state the tab shows, and a DA re-run on the
# same state, whose long table is compared with the table the old build saved. Prints aggregates only; read-only on the
# member. Usage: R_LIBS=<library> Rscript docs/qa/workflow-execution/kit/kit-076/older-session-open.R <member dir> <lipid|metab>
suppressMessages(library(MultiScholaR))
args <- commandArgs(TRUE); member <- args[[1L]]; omic <- args[[2L]]
ns <- asNamespace("MultiScholaR"); api <- function(name) get(name, ns)
h <- new.env(parent = globalenv()); sys.source("tools/profiling/workflow_kit.R", h)
work <- tempfile("kit076-older-"); dir.create(work)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
stopifnot(system2("cp", c("-a", "--reflink=auto", file.path(member, "project"), work)) == 0L)
store_dir <- h$kitFindStore(file.path(work, "project"))
store_root <- dirname(store_dir)
source_dir <- dirname(store_root)
identity <- api("readSavedStoreIdentity")(file.path(store_dir, "store.json"))
manager <- api("newDIAStateAdapter")(api("newWorkflowEntityStore")(store_root, identity$project_id, identity$workflow_id), reopen = TRUE)
quiet <- function(...) invisible(NULL)
session <- shiny::MockShinySession$new()
workflow <- new.env(parent = emptyenv()); workflow$state_manager <- manager; workflow$tab_status <- list()
da <- new.env(parent = emptyenv())
context <- manager$.context
old <- api("withWorkflowStateScope")(manager, function() api("resolveWorkflowStateValue")(context, context$snapshot$da_results))
old_results <- api("withWorkflowStateScope")(manager, function() api("hydrateSavedWorkflowEntity")(context$store, old$products$da_results)$value)
cat(sprintf("%s: head manifest %s, %d result entries\n", omic, old$schema, length(old$results)))
stopifnot(grepl("\\.v1$", old$schema))
record <- api("buildWorkflowAnalysisRecord")(manager)
check <- api("validateWorkflowAnalysisRecord")(record)
cat(sprintf("record valid %s; DA operation %s; entries labelled %s\n", check$valid, record$analyses$differential_abundance$operation,
    paste(vapply(record$analyses$differential_abundance$results, function(r) r$contrast, ""), collapse = ", ")))
stopifnot(isTRUE(check$valid), grepl("_da_output$", record$analyses$differential_abundance$operation))
files <- list.files(source_dir, pattern = "filtered_session_data_latest\\.rds$", full.names = TRUE)
stopifnot(length(files) == 1L)
result <- shiny::withReactiveDomain(session, if (identical(omic, "lipid")) {
    api("loadLipidDaSessionFromFile")(files[[1L]], workflow, da, session, showLoadingFn = function() NULL,
        restorePostReadFn = function(...) NULL, finalizeSuccessFn = function(...) NULL,
        finalizeErrorFn = function(errorMessage, ...) list(errorMessage = errorMessage))
} else {
    api("runMetabDaLoadSessionObserverShell")(files[[1L]], workflow, da, session, showNotification = quiet, removeNotification = quiet,
        logInfo = quiet, logError = quiet, restoreState = function(sessionData, sessionFile, workflowData, daData, session, debugLog) {
            api("restoreMetabDaLoadedSessionState")(sessionData, sessionFile, workflowData, daData, session,
                updateSelectInput = quiet, updateTextAreaInput = quiet, logInfo = quiet, logWarn = quiet)
        })
})
stopifnot(identical(result$status, "success"))
# The copied project's export is not this owner's own, so the Open takes the path an older session file takes: it adopts
# into a fresh owner on the bus and closes the reopened one.
adopted <- !identical(workflow$state_manager, manager)
manager <- workflow$state_manager
context <- manager$.context
shown <- api(if (identical(omic, "lipid")) "lipidDaCurrentS4" else "metabDaCurrentS4")(da, workflow)
cat(sprintf("Open: status %s (%s); the tab shows state %s (class %s, assays %s)\n", result$status, if (adopted) "adopted into a fresh owner" else "the owner reopened",
    da$current_state_name, shown$class, paste(shown$assay_names, collapse = ", ")))
stopifnot(identical(da$current_state_name, manager$current_state))
# The re-run on the same state: the list the display reads, against the table the old build saved.
wrapper <- api(if (identical(omic, "lipid")) "runLipidReferenceDaAnalysis" else "runMetabReferenceDaAnalysis")
parameters <- old$parameters
run <- suppressWarnings(suppressMessages(wrapper(manager, api("getWorkflowStateReference")(manager, da$current_state_name), list(formula_string = parameters$formula,
    da_q_val_thresh = parameters$qvalue, treat_lfc_cutoff = parameters$lfc), old$contrasts)))
long <- if (identical(omic, "lipid")) "da_lipids_long" else "da_metabolites_long"
same <- identical(run$results[[long]], old_results[[long]]) && identical(run$results$per_assay_results, old_results$per_assay_results) &&
    identical(run$results$significant_counts, old_results$significant_counts)
cat(sprintf("re-run: manifest %s, %d long rows; identical to the table the old build saved: %s\n", run$manifest$schema, nrow(run$results[[long]]), same))
stopifnot(same, grepl("\\.v2$", run$manifest$schema))
api("closeDIAStateAdapter")(manager)
cat("RESULT older session opens, records and re-runs\n")
