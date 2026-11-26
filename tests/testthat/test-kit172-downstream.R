# KIT-172: plain log2 downstream: T7 (the steps after the log step on concentration-scale data, with the counts that
# must not move), T8 (the consumers F1b fixed: the filtering figures' detected and missing counts on log2 values, the
# interactive volcano's sentence) and T10 (the `scale` warning through the normalisation tab's run, in process and
# managed). Fixtures from helper-kit172.R.
library(testthat)

kit172Runs <- function() kit172Design()$Run

kit172NonMissing <- function(matrix, features) sum(!is.na(matrix[intersect(features, rownames(matrix)), kit172Runs(), drop = FALSE]))

kit172IntegrationSave <- function(omic, fixture, root, label = "kit172") {
  paths <- list(source_dir = file.path(root, "scripts", omic), integration_dir = file.path(root, "integration"), base_dir = root)
  dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
  failure <- NULL
  kit172Quiet(saveWorkflowSummaryParameters(omicGuiConfig(omic)$summary, fixture$workflow,
    projectDirs = stats::setNames(list(paths), omic), omicType = omic, label = label, description = "d",
    notify = function(message, ...) NULL, done = function(value) NULL, failed = function(error) failure <<- error,
    createWorkflowArgsFn = function(...) "study_parameters.txt", writeAnalysisRecordFn = function(...) NULL,
    catFn = function(...) NULL))
  if (!is.null(failure)) stop(failure)
  file.path(paths$integration_dir, omic)
}

for (class in c("lipid", "metab")) {
  omic <- if (identical(class, "lipid")) "lipidomics" else "metabolomics"

  test_that(paste("KIT-172 T7:", class, "every step runs on concentration data, the counts hold and nothing is imputed"), {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    journey <- kit172Journey(file.path(root, "j"), kit172Object(kit172ConcentrationValues(), class), class,
      ruv = kit172LegacyRuv("automatic"))
    withr::defer(closeDIAStateAdapter(journey$owner))
    prefix <- journey$prefix
    before <- kit172StateMatrix(journey, paste0(prefix, "_itsd_norm"))
    logged <- kit172StateMatrix(journey, paste0(prefix, "_log2"))
    normalised <- kit172StateMatrix(journey, paste0(prefix, "_normalized"))
    corrected <- kit172StateMatrix(journey, journey$owner$current_state)
    features <- rownames(corrected)
    # The count of non-missing cells entering the test equals the count leaving the log step (the features kept), and
    # no later step fills a missing cell.
    expect_identical(kit172NonMissing(corrected, features), kit172NonMissing(logged, features))
    expect_identical(is.na(normalised), is.na(logged[rownames(normalised), , drop = FALSE]))
    expect_identical(is.na(corrected), is.na(logged[features, , drop = FALSE]))
    # The QC tabs' detected and missing counts on the log state equal those on the state before it.
    assay_of <- function(state) methods::slot(journey$owner$getState(state), journey$binding$data_slot)$LCMS_Pos
    id <- methods::slot(journey$owner$getState(paste0(prefix, "_log2")), journey$binding$id_slot)
    # The filtering figures' helpers read the state's log flag (`args$log_transformed`), as their callers pass it.
    logged_flag <- function(state) isTRUE(journey$owner$getState(state)@args$log_transformed)
    detected <- function(state) {
      assay <- assay_of(state)
      if (identical(class, "lipid")) countLipidsPerSample(assay, "Run", id, sample_columns = kit172Runs(), logged = logged_flag(state)) else
        countMetabolitesPerSample(assay, "Run", id, sample_columns = kit172Runs(), logged = logged_flag(state))
    }
    missing <- function(state) {
      assay <- assay_of(state)
      if (identical(class, "lipid")) calculateLipidMissingness(assay, "Run", sample_columns = kit172Runs(), logged = logged_flag(state)) else
        calculateMissingness(assay, "Run", sample_columns = kit172Runs(), logged = logged_flag(state))
    }
    expect_true(logged_flag(paste0(prefix, "_log2")))
    expect_identical(as.data.frame(detected(paste0(prefix, "_log2"))), as.data.frame(detected(paste0(prefix, "_itsd_norm"))))
    expect_identical(missing(paste0(prefix, "_log2")), missing(paste0(prefix, "_itsd_norm")))
    # The DA summary counts the features not tested for a group with no observed value.
    groups <- kit172Design()$group
    unobserved <- sum(rowSums(!is.na(corrected[, groups == "A", drop = FALSE])) == 0 |
      rowSums(!is.na(corrected[, groups == "B", drop = FALSE])) == 0)
    expect_gte(unobserved, 1L)
    da <- omicGuiConfig(omic)$da
    text <- omicDaSummaryText(kit172DaTable(journey), "groupB-groupA", "LCMS_Pos", 0.05, da$summary, da$table_filter %||% "all_only")
    line <- regmatches(text, regexpr("Not tested \\(no value in a compared group\\): [0-9]+", text))
    expect_length(line, 1L)
    expect_identical(as.integer(sub(".*: ", "", line)), as.integer(unobserved))
    # The integration input writer runs, and its matrix holds the missing cells as missing.
    dir <- kit172IntegrationSave(omic, journey$fixture, file.path(root, "project"))
    input <- readIntegrationInputWithMatrices(dir)
    view <- input$views[[paste0(omic, "_LCMS_Pos")]]
    written <- as.matrix(view$matrix[-1L])
    expect_identical(sum(is.na(written)), sum(is.na(normalised[as.character(view$matrix$feature_id), colnames(written), drop = FALSE])))
  })

  test_that(paste("KIT-172 T8:", class, "QC counts read a log2 value of 0 or below as measured on log2 values, a zero as missing on intensities"), {
    logged <- tibble::tibble(ID = c("a", "b", "c"), A_1 = c(-1.5, 0, NA), A_2 = c(2, -0.2, 0.3))
    raw <- tibble::tibble(ID = c("a", "b", "c"), A_1 = c(5, 0, NA), A_2 = c(2, 7, 0.3))
    runs <- c("A_1", "A_2")
    if (identical(class, "lipid")) {
      count <- function(x, ...) as.data.frame(countLipidsPerSample(x, "Run", "ID", sample_columns = runs, ...))$n_detected
      missing <- function(x, ...) calculateLipidMissingness(x, "Run", sample_columns = runs, ...)
    } else {
      count <- function(x, ...) as.data.frame(countMetabolitesPerSample(x, "Run", "ID", sample_columns = runs, ...))$n_detected
      missing <- function(x, ...) calculateMissingness(x, "Run", sample_columns = runs, ...)
    }
    expect_identical(count(logged, logged = TRUE), c(2L, 3L))
    expect_equal(missing(logged, logged = TRUE), 100 / 6)
    # On intensities a zero is not detected, as before.
    expect_identical(count(raw), c(1L, 3L))
    expect_equal(missing(raw), 200 / 6)
  })

  test_that(paste("KIT-172 T8:", class, "the volcano says where it draws missing cells"), {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    journey <- kit172Journey(file.path(root, "j"), kit172Object(kit172ConcentrationValues(), class), class, itsd = FALSE)
    withr::defer(closeDIAStateAdapter(journey$owner))
    prefix <- journey$prefix
    skip_if_not_installed("Glimma")
    results <- journey$da$results
    results$theObject <- readWorkflowOmicsState(journey$binding, journey$owner, getWorkflowStateReference(journey$owner))
    build <- function() if (identical(class, "lipid")) generateLipidDAVolcanoPlotGlimma(results, "groupB-groupA", "LCMS_Pos", 0.05) else
      generateMetabDAVolcanoPlotGlimma(results, "groupB-groupA", "LCMS_Pos", 0.05)
    # The real widget on this fixture (missing cells and untested features): recorded, not asserted (R1b, run
    # 45121681643: glimmaXY stopped with "arguments imply differing number of rows: 1, 0" before this ticket's change).
    real <- kit172Quiet(build())
    cat(sprintf("KIT172-T8-VOLCANO\t%s\treal widget built: %s\n", class, !is.null(real)))
    # The sentence, on a widget built by a stand-in for glimmaXY.
    testthat::local_mocked_bindings(glimmaXY = function(...) htmlwidgets::createWidget("glimmaXY", list()), .package = "Glimma")
    widget <- kit172Quiet(build())
    expect_false(is.null(widget))
    prepended <- paste(vapply(widget$prepend, function(tag) as.character(tag), ""), collapse = "\n")
    expect_true(grepl(omicGuiLog2Texts()[["volcano_missing"]], prepended, fixed = TRUE))
  })
}

# ------------------------------------------------------------------------------------------------
# T10: the `scale` warning through the tab's run.
# ------------------------------------------------------------------------------------------------

#' An owner at the QC-complete state (intensity filter, QC finalise), as the normalisation tab finds it.
kit172PreNorm <- function(root, values, class) {
  fixture <- kit172Owner(root, kit172Object(values, class), class)
  binding <- kit172Binding(class)
  kit172Quiet({
    runWorkflowOmicsReferenceIntensityFilter(binding, fixture$owner, 1, 0.5, fixture$workflow$config_list)
    commitWorkflowOmicsStoreStep(fixture$owner, paste0(binding$prefix, "_qc_finalize"), paste0(binding$prefix, "_qc_complete"),
      fixture$workflow$config_list, "QC complete")
  })
  fixture
}

#' A runner that runs the managed actions in this process: each commit through its declaration's wrapper, its result
#' delivered as the worker's (`stateReference`, `value`); a product action delivers an empty display.
kit172InProcessRunner <- function(manager) {
  list(progress = function(...) list(inc = function(...) invisible(NULL), close = function() invisible(NULL)),
    run = function(actions, done = function() NULL, failed = function(error) stop(error)) {
      outcome <- tryCatch({
        for (action in actions) {
          if (is.function(action$before)) action$before()
          if (is.null(action$context)) {
            action$deliver(list(display = list()))
          } else {
            step <- list(operation = action$operation, parameters = action$parameters, state_name = action$context$state_name,
              description = action$context$description)
            result <- runWorkflowOmicsReferenceStep(manager, step, action$context$config)
            action$deliver(list(stateReference = result$reference, value = result$value))
          }
        }
        NULL
      }, error = function(e) e)
      if (inherits(outcome, "error")) failed(outcome) else done()
    })
}

#' One normalisation run of the tab (route `in_process` or `managed`) with the notifications recorded.
kit172TabRun <- function(fixture, class, route, method, offset = kit172DefaultOffset(class)) {
  binding <- kit172Binding(class)
  binding$qc_plots <- function(...) invisible(NULL)
  normData <- new.env(parent = emptyenv())
  log <- character()
  notes <- list()
  addLog <- function(entry) log <<- c(log, entry)
  testthat::local_mocked_bindings(showNotification = function(ui, action = NULL, duration = 5, closeButton = TRUE, id = NULL,
    type = "default", session = NULL) {
    notes[[length(notes) + 1L]] <<- list(type = type, message = paste(as.character(ui), collapse = " "))
    invisible(NULL)
  }, .package = "shiny")
  settings <- list(apply_itsd = TRUE, itsd_aggregation = "median", itsd_feature_ids = NULL, log_offset = offset,
    norm_method = method, ruv_mode = "skip", ruv_parameters = list(), grouping_variable = "group", shape_variable = "group",
    assay_names = NULL, qc_dir = NULL, experiment_paths = list())
  kit172Quiet(if (identical(route, "in_process")) {
    runWorkflowOmicsInProcessNormalisation(binding, fixture$owner, settings, normData, fixture$workflow$config_list, addLog,
      advance = function(detail) NULL, notify = function(...) shiny::showNotification(...), finish = function() NULL,
      itsd_table = NULL)
  } else {
    runWorkflowOmicsManagedNormalisation(binding, kit172InProcessRunner(fixture$owner), fixture$owner, settings, normData,
      fixture$workflow$config_list, addLog, finish = function() NULL, failed = function(error) addLog(paste("ERROR:", conditionMessage(error))))
  })
  list(log = log, notes = notes, binding = binding)
}

for (class in c("lipid", "metab")) {
  test_that(paste("KIT-172 T10:", class, "scale on log2 values below 0 is warned once in the log and once as a notification"), {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    owners <- list()
    withr::defer(for (owner in owners) closeDIAStateAdapter(owner))
    values <- kit172Values(172001L)$values
    for (route in c("in_process", "managed")) {
      fixture <- kit172PreNorm(file.path(root, route), values, class)
      owners[[length(owners) + 1L]] <- fixture$owner
      run <- kit172TabRun(fixture, class, route, "scale")
      prefix <- run$binding$prefix
      expect_false(any(startsWith(run$log, "ERROR")), info = paste(route, paste(run$log, collapse = " | ")))
      logged <- fixture$owner$getState(paste0(prefix, "_log2"))
      table <- methods::slot(logged, run$binding$data_slot)$LCMS_Pos
      below <- sum(as.matrix(as.data.frame(table)[, kit172Runs()]) < 0, na.rm = TRUE)
      expect_gt(below, 0L)
      text <- sprintf(omicGuiLog2Texts()[["scale_warning"]], omicsLog2CountWords(list(LCMS_Pos = below)))
      expect_identical(sum(grepl(text, run$log, fixed = TRUE)), 1L, info = route)
      # KIT-242 (VAL-242 OD-2 (b); owner 2026-11-27, VAL-373 (Q-B)): the fixture matches no internal standard, so the tab
      # also warns that ITSD was not applied; T10 counts the scale warning by its own text.
      warnings <- Filter(function(note) identical(note$type, "warning") && identical(note$message, text), run$notes)
      expect_length(warnings, 1L)
      expect_identical(warnings[[1L]]$message, text, info = route)
      # `scale` itself is unchanged: the state is normaliseBetweenSamples(scale) of the log2 state.
      direct <- kit172Quiet(normaliseBetweenSamples(theObject = logged, normalisation_method = "scale"))
      direct_table <- as.data.frame(methods::slot(direct, run$binding$data_slot)$LCMS_Pos)
      stored_table <- as.data.frame(methods::slot(fixture$owner$getState(paste0(prefix, "_normalized")), run$binding$data_slot)$LCMS_Pos)
      id <- methods::slot(logged, run$binding$id_slot)
      at <- match(direct_table[[id]], stored_table[[id]])
      expect_identical(unname(as.matrix(stored_table[at, kit172Runs()])), unname(as.matrix(direct_table[, kit172Runs()])), info = route)
    }
  })

  test_that(paste("KIT-172 T10:", class, "no scale warning with cyclic loess, on intensity-scale data or at a supplied offset of 1"), {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    owners <- list()
    withr::defer(for (owner in owners) closeDIAStateAdapter(owner))
    values <- kit172Values(172001L)$values
    intensity <- values
    intensity[!is.na(intensity)] <- 1e4 * intensity[!is.na(intensity)] / min(intensity, na.rm = TRUE)
    cases <- list(list("cyclicloess", values, kit172DefaultOffset(class)), list("scale", intensity, kit172DefaultOffset(class)),
      list("scale", values, 1))
    for (i in seq_along(cases)) {
      fixture <- kit172PreNorm(file.path(root, i), cases[[i]][[2L]], class)
      owners[[length(owners) + 1L]] <- fixture$owner
      run <- kit172TabRun(fixture, class, "in_process", cases[[i]][[1L]], offset = cases[[i]][[3L]])
      expect_false(any(startsWith(run$log, "ERROR")), info = paste(i, paste(run$log, collapse = " | ")))
      expect_false(any(grepl("Scale normalisation rescales", run$log, fixed = TRUE)), info = i)
      # KIT-242 (owner 2026-11-27, VAL-373 (Q-B)): no scale warning; the ITSD notice for this fixture is not counted here.
      expect_length(Filter(function(note) identical(note$type, "warning") &&
        grepl("Scale normalisation rescales", note$message, fixed = TRUE), run$notes), 0L)
    }
  })
}

for (class in c("lipid", "metab")) {
  test_that(paste("KIT-172 T9:", class, "a stored import holding zeros, normalised again at offset 0, gives missing cells and the counts in the log"), {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    values <- kit172Values(172001L)$values
    # Before KIT-118 M3 the import kept zeros; the stored import state is written with three of them.
    zeros <- rbind(c("F020", "A_1"), c("F021", "B_2"), c("F022", "B_3"))
    values[zeros] <- 0
    fixture <- kit172PreNorm(file.path(root, "p"), values, class)
    withr::defer(closeDIAStateAdapter(fixture$owner))
    run <- kit172TabRun(fixture, class, "in_process", "cyclicloess", offset = 0)
    expect_false(any(startsWith(run$log, "ERROR")), info = paste(run$log, collapse = " | "))
    prefix <- run$binding$prefix
    logged <- fixture$owner$getState(paste0(prefix, "_log2"))
    table <- as.data.frame(methods::slot(logged, run$binding$data_slot)$LCMS_Pos)
    id <- methods::slot(logged, run$binding$id_slot)
    kept <- zeros[zeros[, 1L] %in% table[[id]], , drop = FALSE]
    expect_gt(nrow(kept), 0L)
    for (i in seq_len(nrow(kept))) expect_true(is.na(table[[kept[i, 2L]]][match(kept[i, 1L], table[[id]])]), info = kept[i, 1L])
    line <- sprintf(omicGuiLog2Texts()[["nonpositive_log"]], sprintf("LCMS_Pos %d non-positive, 0 non-finite", nrow(kept)))
    expect_true(line %in% run$log, info = paste(run$log, collapse = " | "))
  })
}
