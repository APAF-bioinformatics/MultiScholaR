# KIT-071: the lipidomics and metabolomics sessions on the owner. The normalisation export is the owner's
# publication (the files named as every build before, with the journal and completion record); the DE tab's Open
# reopens the owner on its own export and adopts a legacy RDS from an older build into a fresh owner with `adopted`
# provenance; a session of another shape is refused by name; the summary's final state is the head's. The owner
# fixtures are in helper-workflow-omics-owner.R.
library(testthat)

omicsSessionCase <- function(omic) {
  if (identical(omic, "lipidomics")) {
    list(omic = "lipidomics", type = "LIPIDOMICS", prefix = "lipid", slot = "lipid_data", object = omicsOwnerLipidObject(),
      other = omicsOwnerMetabObject(), counts_field = "lipid_counts",
      states = c("lipid_raw_data_s4", "lipid_correlation_filtered"),
      data_states = c("lipid_correlation_filtered", "lipid_norm_complete", "lipid_ruv_corrected", "lipid_normalized",
        "loaded_for_de", "lipid_qc_complete"))
  } else {
    list(omic = "metabolomics", type = "METABOLOMICS", prefix = "metab", slot = "metabolite_data", object = omicsOwnerMetabObject(),
      other = omicsOwnerLipidObject(), counts_field = "metabolite_counts",
      states = c("metab_raw_data_s4", "metab_correlation_filtered"),
      data_states = c("metab_correlation_filtered", "metab_norm_complete", "metab_ruv_corrected", "metab_normalized",
        "loaded_for_de", "metab_qc_complete"))
  }
}

# A module bus with a bound owner in `source_dir`, the two states committed (the design's and the correlation filter's).
omicsSessionWorkflow <- function(case, source_dir, session, commit = TRUE) {
  workflow <- list2env(list(state_update_trigger = NULL, config_list = list(threshold = 0.5),
    contrasts_tbl = data.frame(friendly_names = "B_vs_A", contrasts = "groupB-groupA", stringsAsFactors = FALSE),
    design_matrix = case$object@design_matrix, qc_params = list(min_pearson_correlation_threshold = 0.75)))
  workflow[[case$counts_field]] <- list(n = 1L)
  workflow$state_manager <- createDIASessionOwner(workflow, list(source_dir = source_dir), session, omic_type = case$omic)
  shiny::withReactiveDomain(session, bindWorkflowOwnerTypeAtDesign(workflow, case$type))
  if (commit) {
    workflow$state_manager$saveState(case$states[[1]], case$object, workflow$config_list, "the design's object")
    workflow$state_manager$saveState(case$states[[2]], case$object, workflow$config_list, "the correlation filter's object")
  }
  workflow
}

omicsSessionNormData <- function(case) {
  list2env(list(normalization_complete = TRUE, ruv_complete = FALSE, correlation_filtering_complete = TRUE,
    itsd_selections = list(), ruv_optimization_results = list(), correlation_results = NULL,
    assay_names = names(methods::slot(case$object, case$slot))), parent = emptyenv())
}

omicsSessionInput <- function() {
  list(norm_method = "quantile", ruv_mode = "skip", apply_itsd = FALSE, itsd_aggregation = NA, log_offset = 1,
    min_pearson_correlation_threshold = 0.75, ruv_grouping_variable = "group")
}

quietFn <- function(...) invisible(NULL)
progressFn <- function(message = NULL, value = NULL, expr) expr

omicsSessionExport <- function(case, workflow, source_dir) {
  norm <- omicsSessionNormData(case)
  if (identical(case$omic, "lipidomics")) {
    handleLipidNormExportSession(omicsSessionInput(), workflow, list(source_dir = source_dir), "Owner test", norm,
      addLog = quietFn, reqFn = identity, showNotificationFn = quietFn, withProgressFn = progressFn,
      incProgressFn = quietFn, logInfoFn = quietFn, logWarnFn = quietFn, logErrorFn = function(message) stop(message))
  } else {
    runMetabNormExportSessionWorkflow(workflow, norm, omicsSessionInput(), "Owner test", source_dir,
      withProgressFn = progressFn, incProgressFn = quietFn, logInfoFn = quietFn)
  }
}

omicsSessionOpen <- function(case, workflow, da, session_file, session) {
  shiny::withReactiveDomain(session, if (identical(case$omic, "lipidomics")) {
    loadLipidDaSessionFromFile(session_file, workflow, da, session, showLoadingFn = function() NULL,
      restorePostReadFn = function(...) NULL, finalizeSuccessFn = function(...) NULL,
      finalizeErrorFn = function(errorMessage, ...) list(errorMessage = errorMessage))
  } else {
    runMetabDaLoadSessionObserverShell(session_file, workflow, da, session, showNotification = quietFn,
      removeNotification = quietFn, logInfo = quietFn, logError = quietFn,
      restoreState = function(sessionData, sessionFile, workflowData, daData, session, debugLog) {
        restoreMetabDaLoadedSessionState(sessionData, sessionFile, workflowData, daData, session,
          updateSelectInput = quietFn, updateTextAreaInput = quietFn, logInfo = quietFn, logWarn = quietFn)
      })
  })
}

# The component entity ids of an object's checkpoint on `owner`, by assay element (content-addressed).
omicsSessionComponentIds <- function(owner, object, slot) {
  checkpoint <- checkpointWorkflowExecutionResult(owner, object)
  sort(unname(vapply(checkpoint$components[[slot]], function(element) element$entity_id, character(1))))
}

omicsSessionHeadIds <- function(owner, slot) {
  record <- readWorkflowStateMetadata(owner)
  sort(unname(vapply(record$data$components[[slot]], function(element) element$entity_id, character(1))))
}

omicsSessionLegacyList <- function(case, object = case$object, omic = case$omic, state = case$states[[2]]) {
  list(r6_current_state_name = state, current_s4_object = object,
    contrasts_tbl = data.frame(friendly_names = "B_vs_A", contrasts = "groupB-groupA", stringsAsFactors = FALSE),
    design_matrix = case$object@design_matrix, config_list = list(threshold = 0.5),
    itsd_selections = list(), ruv_optimization_results = list(), correlation_results = NULL,
    assay_names = names(methods::slot(case$object, case$slot)),
    export_timestamp = as.POSIXct("2026-10-28 12:00:00", tz = "UTC"), omic_type = omic, experiment_label = "legacy",
    normalization_method = "quantile", ruv_mode = "skip", itsd_applied = FALSE, itsd_aggregation = NA, log_offset = 1,
    correlation_threshold = 0.75, ruv_grouping_variable = "group", feature_counts = NULL, qc_params = list(),
    normalization_complete = TRUE, ruv_complete = FALSE, correlation_filtering_complete = TRUE)
}

for (omic in c("lipidomics", "metabolomics")) {
  case <- omicsSessionCase(omic)

  test_that(sprintf("KIT-071: the %s export is the owner's publication and the DE tab's Open reopens the owner", omic), {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    source_dir <- file.path(root, "source"); dir.create(source_dir)
    session <- omicsOwnerSession()
    workflow <- omicsSessionWorkflow(case, source_dir, session)
    withr::defer(closeDIAStateAdapter(workflow$state_manager))
    owner <- workflow$state_manager
    head_before <- owner$.context$head$entity_id
    head_ids <- omicsSessionHeadIds(owner, case$slot)

    artifacts <- omicsSessionExport(case, workflow, source_dir)
    expect_false(isFALSE(artifacts))
    latest <- file.path(source_dir, paste0(case$prefix, "_filtered_session_data_latest.rds"))
    expect_identical(normalizePath(artifacts$latestFilepath), normalizePath(latest))
    expect_true(file.exists(artifacts$sessionFilepath))
    expect_match(basename(artifacts$sessionFilepath), paste0("^", case$prefix, "_filtered_session_data_[0-9]{8}_[0-9]{6}[.]rds$"))
    expect_true(file.exists(latest))
    expect_true(file.exists(file.path(source_dir, paste0(case$prefix, "_filtered_session_summary.txt"))))
    expect_true(file.exists(file.path(source_dir, paste0(case$prefix, "_qc_params.RDS"))))
    expect_identical(artifacts$publication$status, "complete")
    expect_identical(artifacts$publication$latest_name, basename(latest))
    expect_true(file.exists(file.path(source_dir, paste0("export_complete_", artifacts$publication$request_id, ".json"))))
    request <- owner$.context$owned_export
    expect_identical(request$generation, head_before)
    expect_identical(normalizePath(request$path), normalizePath(latest))
    expect_null(owner$.context$export_pending)
    # The file is the legacy shape: the session list with the head's state as its object, and it is the head's by content.
    exported <- readRDS(latest)
    expect_identical(exported$omic_type, omic)
    expect_identical(exported$r6_current_state_name, case$states[[2]])
    expect_s4_class(exported$current_s4_object, class(case$object)[[1]])
    expect_identical(omicsSessionComponentIds(owner, exported$current_s4_object, case$slot), head_ids)
    expect_identical(owner$.context$head$entity_id, head_before)

    # Open on the owner's own export: the owner is reopened at its head; no new commit, the same owner on the bus.
    da <- new.env(parent = emptyenv())
    result <- omicsSessionOpen(case, workflow, da, latest, session)
    expect_identical(result$status, "success")
    expect_identical(workflow$state_manager, owner)
    expect_identical(owner$.context$head$entity_id, head_before)
    expect_identical(as.character(unlist(owner$getHistory())), c("initial", case$states))
    expect_identical(da$current_state_name, case$states[[2]])
    expect_identical(omicsSessionHeadIds(owner, case$slot), head_ids)
    expect_s4_class(owner$getState(da$current_state_name), class(case$object)[[1]])
    # Another file is not the owner's export.
    expect_null(resolveWorkflowOmicsOwnedSession(owner, file.path(source_dir, "other.rds")))
  })

  test_that(sprintf("KIT-071: a legacy %s session adopts into a fresh owner with adopted provenance; other shapes are refused by name", omic), {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    source_dir <- file.path(root, "source"); dir.create(source_dir)
    session <- omicsOwnerSession()
    workflow <- omicsSessionWorkflow(case, source_dir, session, commit = FALSE)
    previous <- workflow$state_manager
    withr::defer(closeDIAStateAdapter(workflow$state_manager))
    legacy <- file.path(source_dir, paste0(case$prefix, "_filtered_session_data_latest.rds"))
    saveRDS(omicsSessionLegacyList(case), legacy)

    da <- new.env(parent = emptyenv())
    result <- omicsSessionOpen(case, workflow, da, legacy, session)
    expect_identical(result$status, "success")
    adopted <- workflow$state_manager
    expect_s3_class(adopted, "DIAStateAdapter")
    expect_false(identical(adopted, previous))
    expect_true(previous$.context$store$closed)
    expect_identical(adopted$.context$store$identity$workflow_id, case$type)
    expect_identical(adopted$.context$store$identity$project_id, paste0(case$type, "-native"))
    expect_identical(adopted$.context$omic_type, omic)
    expect_identical(adopted$.context$experiment_paths$source_dir, source_dir)
    expect_identical(adopted$current_state, case$states[[2]])
    expect_identical(da$current_state_name, case$states[[2]])
    expect_null(session$userData$dia_capacity_manager)
    provenance <- adopted$.context$snapshot$step_provenance[[case$states[[2]]]]
    expect_identical(provenance$operation, "legacy")
    expect_true(isTRUE(provenance$adopted))
    record <- readWorkflowStateMetadata(adopted)
    expect_identical(record$data$class, class(case$object)[[1]])
    expect_identical(record$config$loaded_from, legacy)
    expect_true(isTRUE(record$config$adopted))
    expect_identical(omicsSessionHeadIds(adopted, case$slot), omicsSessionComponentIds(adopted, case$object, case$slot))
    metadata <- resolveDIAStateValue(adopted$.context, adopted$.context$native_metadata_reference)
    expect_identical(metadata$omic_type, omic)
    expect_identical(metadata$experiment_label, "legacy")
    expect_null(metadata$current_s4_object)
    expect_s4_class(adopted$getState(), class(case$object)[[1]])

    # The other omic's object, the other omic's session, a descriptor and no state name are refused by name.
    expect_error(shiny::withReactiveDomain(session, ingestWorkflowOmicsLegacySession(
      omicsSessionLegacyList(case, object = case$other), workflow, source_dir, omic)), class = "workflow_session_unsupported")
    expect_error(shiny::withReactiveDomain(session, ingestWorkflowOmicsLegacySession(
      omicsSessionLegacyList(case, omic = setdiff(c("lipidomics", "metabolomics"), omic)), workflow, source_dir, omic)),
      class = "workflow_session_unsupported")
    expect_error(shiny::withReactiveDomain(session, ingestWorkflowOmicsLegacySession(
      omicsSessionLegacyList(case, object = list(schema = 1L, entity_id = strrep("a", 64), payload = list())), workflow, source_dir, omic)),
      class = "saved_entity_corrupt")
    expect_error(shiny::withReactiveDomain(session, ingestWorkflowOmicsLegacySession(
      omicsSessionLegacyList(case, state = NULL), workflow, source_dir, omic, state_name = NULL)), class = "workflow_session_unsupported")
    expect_identical(workflow$state_manager, adopted)
    expect_false(adopted$.context$store$closed)
    # Through the tab: a refused file leaves the owner in place and reports the error.
    wrong <- file.path(source_dir, "wrong.rds")
    saveRDS(omicsSessionLegacyList(case, object = case$other), wrong)
    refused <- omicsSessionOpen(case, workflow, new.env(parent = emptyenv()), wrong, session)
    expect_identical(refused$status, "error")
    expect_match(refused$errorMessage, "session holds a", fixed = TRUE)
    expect_identical(workflow$state_manager, adopted)
    # A manager that is not an owner is a test double: nothing is adopted.
    double <- list2env(list(state_manager = list(saveState = function(...) NULL)))
    expect_false(ingestWorkflowOmicsLegacySession(omicsSessionLegacyList(case), double, source_dir, omic))
  })

  test_that(sprintf("KIT-071: the %s summary's final state is the head's", omic), {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    source_dir <- file.path(root, "source"); dir.create(source_dir)
    session <- omicsOwnerSession()
    workflow <- omicsSessionWorkflow(case, source_dir, session)
    withr::defer(closeDIAStateAdapter(workflow$state_manager))
    owner <- workflow$state_manager
    expect_identical(workflowSummaryFinalStateName(owner, case$data_states), case$states[[2]])
    picked <- if (identical(omic, "lipidomics")) getLipidSummaryStateObject(workflow) else getMetabSummaryStateObject(workflow)
    expect_s4_class(picked, class(case$object)[[1]])
    expect_identical(omicsSessionComponentIds(owner, picked, case$slot), omicsSessionHeadIds(owner, case$slot))
    # The head reverted to a state that is not a data state: the scan over the head's history finds none, and the
    # site falls back to the current state as before.
    owner$revertToState(case$states[[1]])
    expect_null(workflowSummaryFinalStateName(owner, case$data_states))
    picked <- if (identical(omic, "lipidomics")) getLipidSummaryStateObject(workflow) else getMetabSummaryStateObject(workflow)
    expect_s4_class(picked, class(case$object)[[1]])
    # A double keeps the ordered scan of its history, through the site's detect seam.
    double <- list(getHistory = function() c(case$states[[1]], case$data_states[[4]], case$data_states[[2]]))
    expect_identical(workflowSummaryFinalStateName(double, case$data_states), case$data_states[[2]])
    expect_identical(workflowSummaryFinalStateName(double, case$data_states, detectFn = purrr::detect), case$data_states[[2]])
  })
}
