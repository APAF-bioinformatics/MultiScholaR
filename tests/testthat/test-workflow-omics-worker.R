# KIT-078: the omics steps on the worker.
library(testthat)

test_that("KIT-078: every omics compute declaration runs from a frozen request through the worker's entry point", {
  skip_if_not(workflowOwnerCapable())
  testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64), .package = "MultiScholaR")
  fixtures <- omicsWorkerFixtures(withr::local_tempdir())
  for (omic in names(fixtures)) {
    fixture <- fixtures[[omic]]$fixture
    binding <- fixtures[[omic]]$binding
    withr::defer(closeDIAStateAdapter(fixture$owner))
    steps <- list(omicsIntensityFilterStep(binding, 1, 0.5), omicsDuplicateResolutionStep(binding),
      omicsItsdNormalisationStep(binding, "median", NULL), omicsLog2Step(binding, 1),
      omicsBetweenSampleNormalisationStep(binding, "cyclicloess"), omicsRuvStep(binding, fixtures[[omic]]$ruv))
    for (step in steps) {
      result <- omicsWorkerRun(fixture, step$operation, step$parameters)
      expect_identical(result$status, "complete", info = paste(omic, step$operation, result$error))
      expect_identical(result$result$operation, step$operation, info = omic)
      expect_true(isDIAComponentCheckpoint(result$result$checkpoint), info = step$operation)
      expect_identical(result$result$checkpoint$class, binding$class, info = step$operation)
      expect_true(isSavedWorkflowReference(result$result$reference), info = step$operation)
      expect_true(length(result$result$roots) >= 1L, info = step$operation)
      expect_true(length(result$result$owned) >= 1L, info = step$operation)
      # The worker moved no head: the source is still the session's state.
      expect_false(step$state_name %in% unlist(fixture$owner$getHistory()), info = step$operation)
    }
    # The filtering steps hand over their metrics; a later step does not.
    filtering <- omicsWorkerRun(fixture, steps[[1L]]$operation, steps[[1L]]$parameters)
    expect_true(isSavedWorkflowReference(filtering$result$display$filtering_metrics), info = omic)
    expect_null(omicsWorkerRun(fixture, steps[[4L]]$operation, steps[[4L]]$parameters)$result$display, info = omic)
    # The product-only declarations return the source's checkpoint.
    source <- readWorkflowStateMetadata(fixture$owner)$data
    pairs <- omicsWorkerRun(fixture, paste0(binding$prefix, "_pearson_pairs"), list(correlation_group = "group"))
    expect_identical(pairs$status, "complete", info = pairs$error)
    expect_identical(pairs$result$checkpoint, source, info = omic)
    contrast <- omicsWorkerRun(fixture, paste0(binding$prefix, "_da_contrast"), list(formula = "~ 0 + group", qvalue = 0.05, lfc = 0,
      contrast = as.character(fixture$workflow$contrasts_tbl$contrasts[[1L]]),
      friendly_name = as.character(fixture$workflow$contrasts_tbl$friendly_names[[1L]])))
    expect_identical(contrast$status, "complete", info = contrast$error)
    expect_identical(contrast$result$checkpoint, source, info = omic)
    expect_true(length(contrast$result$native) >= 1L, info = omic)
    # Every declaration the registry puts on the role is covered above; the store steps refuse dispatch.
    expect_setequal(Filter(function(name) identical(workflowActionSpec(name)$omic_type, omic) && identical(workflowActionSpec(name)$compute, "omics") &&
      !is.null(workflowActionSpec(name)$declaration), workflowOperationNames()), c(vapply(steps, `[[`, "", "operation"),
      paste0(binding$prefix, c("_pearson_pairs", "_pearson_selection", "_da_contrast"))))
    for (name in paste0(binding$prefix, c("_qc_finalize", "_norm_complete", "_reset", "_da_output")))
      expect_error(assertWorkflowOperationDispatchable(name), class = "workflow_execution_unsupported")
  }
})

test_that("KIT-078: a worker's QC figures are staged and saved, and it leaves the head to the controller", {
  skip_if_not(workflowOwnerCapable())
  testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64), .package = "MultiScholaR")
  fixtures <- omicsWorkerFixtures(withr::local_tempdir())
  for (omic in names(fixtures)) {
    fixture <- fixtures[[omic]]$fixture
    binding <- fixtures[[omic]]$binding
    withr::defer(closeDIAStateAdapter(fixture$owner))
    # The fixture's duplicate ids would stop every plot; the module's QC figures follow the QC tabs.
    suppressWarnings(suppressMessages(runWorkflowOmicsReferenceDuplicateResolution(binding, fixture$owner, list())))
    head <- fixture$owner$.context$head
    result <- omicsWorkerRun(fixture, paste0(binding$prefix, "_qc_figures"), list(stage = "post_filter", grouping_variable = "group", shape_variable = NULL))
    expect_identical(result$status, "complete", info = result$error)
    staged <- list.files(file.path(result$directory, "omics-publication", binding$qc_dir_role))
    expect_true(any(grepl("_pre_norm_pca[.]png$", staged)), info = omic)
    expect_identical(fixture$owner$.context$head, head)
    # Each figure file's product and asset wait for the controller's root.
    expect_gte(length(result$result$roots), 2L * length(staged))
  }
})

test_that("KIT-078: updating the filtering progress from a worker's metrics equals updating it from the object", {
  object <- lipidDeclarationsObject()
  metrics <- suppressMessages(computeLipidFilteringMetrics(object))
  expect_identical(names(metrics), c("assay_names", "metrics_list_this_step", "total_lipids"))
  expect_identical(metrics$assay_names, names(object@lipid_data))
  object <- metabDeclarationsObject()
  metrics <- suppressMessages(computeMetabFilteringMetrics(object))
  expect_identical(names(metrics), c("assay_names", "metrics_list_this_step", "total_metabolites"))
  expect_identical(metrics$assay_names, names(object@metabolite_data))
})

test_that("KIT-078: every omics declaration gives the same products on the worker as in process (scientific contract)", {
  skip_if_not(workflowOwnerCapable())
  testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
    finishWorkflowExecutionJob = function(owner, status) { owner$job <- NULL; list(status = status) }, .package = "MultiScholaR")
  root <- withr::local_tempdir()
  for (omic in c("lipidomics", "metabolomics")) {
    lipid <- identical(omic, "lipidomics")
    binding <- omicsDeclarationBindings()[[omic]]
    owner_of <- function(path) if (lipid) lipidDeclarationsOwner(path) else metabDeclarationsOwner(path)
    in_process <- owner_of(file.path(root, omic, "in-process"))
    worker <- owner_of(file.path(root, omic, "worker"))
    withr::defer({ closeDIAStateAdapter(in_process$owner); closeDIAStateAdapter(worker$owner) })
    ruv <- if (lipid) lipidDeclarationsRuvParameters() else metabDeclarationsRuvParameters()
    journey_in_process <- if (lipid) lipidDeclarationsJourney else metabDeclarationsJourney
    expected <- suppressWarnings(suppressMessages(omicsDeclarationsJourneyDump(in_process, journey_in_process)))
    actual <- suppressWarnings(suppressMessages(omicsDeclarationsJourneyDump(worker, function(fixture) omicsWorkerJourney(fixture, binding, ruv))))
    expect_identical(actual$history, expected$history, info = omic)
    for (state in expected$history) expect_identical(actual$states[[state]], expected$states[[state]], info = paste(omic, state))
    expect_identical(actual$values, expected$values, info = omic)
    expect_identical(actual$correlations, expected$correlations, info = omic)
    expect_identical(actual$da, expected$da, info = omic)
    # The head's provenance is what the wrappers record: operation, variant, values and the pairs before the selection.
    provenance <- function(fixture) omicsDeclarationsScrub(lapply(fixture$owner$.context$snapshot$step_provenance, function(entry) {
      if (!is.null(entry$values$pairs)) entry$values$pairs <- "<entity>"
      entry$preceded_by <- lapply(entry$preceded_by, function(previous) previous[c("operation", "variant", "parameter_names", "values")])
      entry[setdiff(names(entry), "committed_at")]
    }))
    expect_identical(provenance(worker), provenance(in_process), info = omic)
    # Every product and table the session names is rooted in the head.
    snapshot <- worker$owner$.context$snapshot
    science <- paste0(binding$prefix, c("_intensity_filter", "_duplicate_resolution", "_itsd_normalisation", "_log2", "_between_sample_normalisation", "_ruv", "_pearson_selection"), "_science")
    for (state in expected$history[-(1:2)]) {
      products <- readWorkflowStateMetadata(worker$owner, snapshot$states[[state]])$data$products
      for (name in intersect(science, names(products))) expect_false(is.null(snapshot$products[[products[[name]]$entity_id]]), info = paste(state, name))
    }
    expect_gt(length(snapshot$owned_components), 0L)
    expect_true(all(vapply(snapshot$da_dependencies, length, 1L) > 0L), info = omic)
    # The record validates and replays as the in-process one does.
    expect_true(validateWorkflowAnalysisRecord(buildWorkflowAnalysisRecord(worker$owner))$valid, info = omic)
  }
})

test_that("KIT-078: the worker's QC figures reach the project through the journal, rooted, and a publication that fails is retried, not recomputed", {
  skip_if_not(workflowOwnerCapable())
  calls <- new.env()
  testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
    finishWorkflowExecutionJob = function(owner, status) { calls$status <- status; owner$job <- NULL; list(status = status) }, .package = "MultiScholaR")
  fixture <- lipidDeclarationsOwner(withr::local_tempdir())
  withr::defer(closeDIAStateAdapter(fixture$owner))
  suppressWarnings(suppressMessages(runWorkflowOmicsReferenceDuplicateResolution(lipidomicsDeclarationBinding, fixture$owner, list())))
  controller <- omicsWorkerController(fixture)
  state <- fixture$owner$current_state
  history <- fixture$owner$getHistory()
  provenance <- fixture$owner$.context$pending_provenance
  expect_length(list.files(fixture$paths$lipid_qc_dir, "[.]png$"), 0L)
  done <- omicsWorkerAction(fixture, controller, "lipid_qc_figures", list(stage = "post_filter", grouping_variable = "group", shape_variable = NULL))
  expect_identical(done$acceptance, "complete")
  published <- list.files(fixture$paths$lipid_qc_dir, "[.]png$", full.names = TRUE)
  expect_true(length(published) > 0L)
  # It commits nothing and leaves the next step's provenance as it was.
  expect_identical(fixture$owner$current_state, state)
  expect_identical(fixture$owner$getHistory(), history)
  expect_identical(fixture$owner$.context$pending_provenance, provenance)
  # Every published file is byte-identical to a saved, rooted asset of the controller's store.
  snapshot <- fixture$owner$.context$snapshot
  assets <- Filter(function(reference) identical(reference$kind, "rendered_asset"), snapshot$products)
  hashes <- vapply(assets, function(reference) digest::digest(hydrateSavedWorkflowEntity(fixture$owner$.context$store, reference)$bytes, algo = "sha256"), "")
  for (file in published) expect_true(digest::digest(readBin(file, "raw", n = file.size(file)), algo = "sha256") %in% hashes, info = basename(file))

  # A publication that fails is pending: the next action is refused until the retry publishes it, from the same request.
  unlink(published)
  failing <- TRUE
  testthat::local_mocked_bindings(workflowExportBoundary = function(phase) {
    if (identical(phase, "before_action_outputs") && failing) stop("disk full")
    invisible(NULL)
  }, .package = "MultiScholaR")
  id <- paste0("wfx-", strrep("b", 32))
  store <- fixture$owner$.context$store
  paths <- workflowExportPaths(store$directory, id)
  token <- list(request_id = id, incarnation = store$identity$incarnation, expected_head = workflowExecutionHead(store)$entity_id,
    generation = controller$control$generation, controller_id = controller$control$controller_id)
  writeWorkflowExecutionMessage(file.path(paths$request_dir, "publication-context.rds"), list(experiment_label = "test", experiment_paths = fixture$paths, omic_type = "lipidomics"))
  writeWorkflowExportJournal(paths, list(schema = "multischolar.workflow_action.v1", request_id = id, status = "admitting",
    expected_head = token$expected_head, token_hash = digest::digest(token, algo = "sha256")))
  message <- omicsWorkerRun(fixture, "lipid_qc_figures", list(stage = "post_filter", grouping_variable = "group", shape_variable = NULL), directory = paths$request_dir, token = token)
  controller$job <- list(request = list(token = token, operation = "lipid_qc_figures", config = list(mode = "adapters"), parameters = list()),
    paths = paths, completion_context = NULL)
  receipt <- acceptWorkflowActionResult(controller, message[c("schema", "status", "token", "result")])
  expect_identical(receipt$status, "publication_pending")
  expect_length(list.files(fixture$paths$lipid_qc_dir, "[.]png$"), 0L)
  # `assertWorkflowExportReady`'s action-output guard refuses dependent work; no state list is needed (recorded divergence).
  expect_error(assertWorkflowExportReady(fixture$owner), class = "saved_entity_outputs_pending")
  failing <- FALSE
  retried <- resumeWorkflowAction(controller, id)
  expect_identical(retried$operation, "lipid_qc_figures")
  expect_identical(sort(list.files(fixture$paths$lipid_qc_dir, "[.]png$", full.names = TRUE)), sort(published))
  expect_null(assertWorkflowExportReady(fixture$owner))
})

test_that("KIT-078: an emptied omics assay is saved as the kernel's result, not taken for an untouched one", {
  skip_if_not(workflowOwnerCapable())
  fixture <- lipidDeclarationsOwner(withr::local_tempdir())
  withr::defer(closeDIAStateAdapter(fixture$owner))
  owner <- fixture$owner
  source <- getWorkflowStateReference(owner)
  object <- readLipidWorkflowState(owner)
  emptied <- object
  emptied@lipid_data <- lapply(object@lipid_data, function(assay) assay[0L, , drop = FALSE])
  emptied <- inheritWorkflowExecutionSource(emptied, owner, source)
  expect_setequal(attr(emptied, "dia_execution_fields"), c("design_matrix", "lipid_data"))
  checkpoint <- checkpointWorkflowExecutionResult(owner, emptied)
  previous <- readWorkflowStateMetadata(owner, source)$data$components$lipid_data
  expect_false(identical(checkpoint$components$lipid_data[[1L]]$entity_id, previous[[1L]]$entity_id))
  # The protein class's fields are the same three as before.
  expect_setequal(with(workflowClassDeclaration("ProteinQuantitativeData"), c(tables, lists)), c("protein_quant_table", "design_matrix", "protein_id_table"))
})

test_that("KIT-078: the controller materialises no assay table of the omic during any managed step (the objective's probe); the same tabs in process do", {
  skip_if_not(workflowOwnerCapable())
  root <- withr::local_tempdir()
  for (omic in c("lipidomics", "metabolomics")) {
    lipid <- identical(omic, "lipidomics")
    binding <- omicsDeclarationBindings()[[omic]]
    owner_of <- function(path) if (lipid) lipidDeclarationsOwner(path) else metabDeclarationsOwner(path)
    for (managed in c(TRUE, FALSE)) {
      fixture <- owner_of(file.path(root, omic, if (managed) "managed" else "in-process"))
      withr::defer(closeDIAStateAdapter(fixture$owner))
      probe <- new.env()
      local({
        if (managed) omicsManagedChainMocks(fixture, probe)
        # The tabs in process: every action counts, from the click to the tab's completion.
        else {
          table <- materializeWorkflowTable
          components <- readWorkflowExecutionComponents
          testthat::local_mocked_bindings(materializeWorkflowTable = function(context, reference, ...) {
            probe$tables <- c(probe$tables, reference$entity_id)
            table(context, reference, ...)
          }, readWorkflowExecutionComponents = function(manager, reference, fields = NULL, ...) {
            probe$reads <- c(probe$reads, paste(fields %||% "<every field>", collapse = ","))
            components(manager, reference, fields = fields, ...)
          }, .package = "MultiScholaR")
        }
        assays <- names(methods::slot(fixture$object, binding$data_slot))
        # A manual internal-standard selection: the ids the shown table holds, resolved without the matrix on a managed run.
        id_column <- methods::slot(fixture$object, binding$id_slot)
        ids <- stats::setNames(lapply(assays, function(assay) methods::slot(fixture$object, binding$data_slot)[[assay]][[id_column]]), assays)
        selections <- list(assays = assays, candidates = ids, rows = stats::setNames(lapply(assays, function(assay) 1:2), assays))
        if (!managed) selections <- list(assays = assays)
        suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(omicsTabServer(fixture, omic, probe, managed), {
          session$userData$probe_owner <- fixture$owner
          omicsTabJourney(session, omic, selections)
        }))))
      })
      history <- unlist(fixture$owner$getHistory())
      prefix <- binding$prefix
      expect_true(all(paste0(prefix, c("_intensity_filtered", "_duplicates_resolved", "_qc_complete", "_itsd_norm", "_log2", "_normalized", "_ruv_corrected", "_correlation_filtered")) %in% history),
        info = paste(omic, managed, paste(history, collapse = " ")))
      expect_identical(fixture$owner$.context$snapshot$da_results$kind %||% "absent" != "absent", TRUE, info = paste(omic, managed))
      materialised <- intersect(probe$tables, omicsAssayTableIds(fixture$owner, binding))
      reads <- Filter(function(read) grepl(binding$data_slot, read, fixed = TRUE) || identical(read, "<every field>"), probe$reads)
      if (managed) {
        expect_length(materialised, 0L)
        expect_length(reads, 0L)
        # Every compute declaration was dispatched, with the figures and the export, each closed on delivery.
        expect_true(all(paste0(prefix, c("_intensity_filter", "_duplicate_resolution", "_qc_figures", "_itsd_normalisation", "_log2",
          "_between_sample_normalisation", "_ruv", "_pearson_pairs", "_pearson_selection", "_da_contrast", "_da_publication")) %in% probe$operations),
          info = paste(omic, paste(probe$operations, collapse = " ")))
        expect_true(all(probe$statuses == "delivered"), info = omic)
        expect_identical(sum(probe$operations == paste0(prefix, "_qc_figures")), 3L, info = omic)
        # The tab's completion ran: the normalisation's references, and the DA's files through the journal.
        qc_dir <- fixture$paths[[binding$qc_dir_role]]
        expect_true(any(grepl("_ruv_corrected_pca[.]png$", list.files(qc_dir))), info = omic)
        expect_true(length(list.files(fixture$paths$da_output_dir, recursive = TRUE)) > 0L, info = omic)
        expect_identical(fixture$owner$.context$snapshot$step_provenance[[paste0(prefix, "_itsd_norm")]]$values$itsd_feature_ids,
          stats::setNames(lapply(names(methods::slot(fixture$object, binding$data_slot)), function(assay) methods::slot(fixture$object, binding$data_slot)[[assay]][[methods::slot(fixture$object, binding$id_slot)]][1:2]),
            names(methods::slot(fixture$object, binding$data_slot))), info = omic)
      } else {
        # The positive control: the probe sees the in-process tabs read the matrix.
        expect_gt(length(reads), 0L)
      }
    }
  }
})
