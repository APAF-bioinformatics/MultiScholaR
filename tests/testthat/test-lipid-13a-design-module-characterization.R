library(testthat)
library(shiny)

makeLipidDesignCharacterizationHarness <- function() {
  data_tbl <- shiny::reactive(
    list(
      LCMS_Pos = data.frame(
        database_identifier = c("L1", "L2"),
        annotation = c("Annot1", "Annot2"),
        Sample_1 = c(10, 20),
        Sample_2 = c(30, 40),
        stringsAsFactors = FALSE
      )
    )
  )
  config_list <- shiny::reactive(
    list(deAnalysisParameters = list(formula_string = "~ 0 + group"))
  )
  column_mapping <- shiny::reactive(
    list(
      lipid_id_col = "database_identifier",
      annotation_col = "annotation",
      sample_columns = c("Sample_1", "Sample_2")
    )
  )

  session <- shiny::MockShinySession$new()

  returned <- shiny::withReactiveDomain(session, {
    mod_lipid_design_builder_server(
      "builder",
      data_tbl = data_tbl,
      config_list = config_list,
      column_mapping = column_mapping
    )
  })

  list(
    session = session,
    returned = returned
  )
}

setLipidDesignInputs <- function(session, ...) {
  values <- list(...)
  names(values) <- paste0("builder-", names(values))
  do.call(session$setInputs, values)
  shiny:::flushReact()
}

test_that("mod_lipid_design_builder_server preserves the empty save-results public behavior", {
  harness <- makeLipidDesignCharacterizationHarness()

  expect_null(shiny::isolate(harness$returned()))

  setLipidDesignInputs(
    harness$session,
    formula_string = "~ 0 + group",
    save_results = 1
  )

  expect_null(shiny::isolate(harness$returned()))
})

test_that("buildLipidDesignInitialState preserves fresh initialization", {
  skip_if_not(
    exists("buildLipidDesignInitialState", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  build_initial_state <- get("buildLipidDesignInitialState", envir = asNamespace("MultiScholaR"))

  assay_list <- list(
    LCMS_Pos = data.frame(
      database_identifier = c("L1", "L2"),
      annotation = c("A1", "A2"),
      Sample_10 = c(1, 2),
      Sample_2 = c(3, 4),
      stringsAsFactors = FALSE
    )
  )

  state <- build_initial_state(
    assayList = assay_list,
    configList = list(deAnalysisParameters = list(formula_string = "~ 0 + group")),
    colMap = list(
      lipid_id_col = "database_identifier",
      annotation_col = "annotation",
      sample_columns = c("Sample_10", "Sample_2")
    ),
    messageFn = function(...) invisible(NULL),
    logErrorFn = function(...) stop("unexpected log error")
  )

  expect_identical(state$formula, "~ 0 + group")
  expect_identical(state$sample_names, c("Sample_10", "Sample_2"))
  expect_identical(state$design_matrix$Run, c("Sample_2", "Sample_10"))
  expect_identical(state$groups, character(0))
  expect_identical(state$factors, character(0))
  expect_s3_class(state$contrasts, "data.frame")
  expect_identical(nrow(state$contrasts), 0L)
})

test_that("buildLipidDesignInitialState preserves existing-design reuse", {
  skip_if_not(
    exists("buildLipidDesignInitialState", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  build_initial_state <- get("buildLipidDesignInitialState", envir = asNamespace("MultiScholaR"))

  assay_list <- list(
    LCMS_Pos = data.frame(
      database_identifier = c("L1", "L2"),
      annotation = c("A1", "A2"),
      Sample_1 = c(1, 2),
      Sample_2 = c(3, 4),
      stringsAsFactors = FALSE
    )
  )
  existing_design <- data.frame(
    Run = c("Sample_1", "Sample_2"),
    group = c("Control", "Treatment"),
    factor1 = c("Condition", "Condition"),
    factor2 = c("", ""),
    factor3 = c(NA_character_, NA_character_),
    replicates = c(1L, 2L),
    tech_reps = c(NA_integer_, NA_integer_),
    stringsAsFactors = FALSE
  )
  existing_contrasts <- data.frame(
    contrasts = "groupTreatment-groupControl",
    stringsAsFactors = FALSE
  )

  state <- build_initial_state(
    assayList = assay_list,
    configList = list(deAnalysisParameters = list(formula_string = "~ 0 + group")),
    colMap = list(
      lipid_id_col = "database_identifier",
      annotation_col = "annotation",
      sample_columns = c("Sample_1", "Sample_2")
    ),
    existingDesignMatrix = existing_design,
    existingContrasts = existing_contrasts,
    messageFn = function(...) invisible(NULL),
    logErrorFn = function(...) stop("unexpected log error")
  )

  expect_identical(state$design_matrix, existing_design)
  expect_identical(state$data_cln, assay_list)
  expect_identical(state$groups, c("Control", "Treatment"))
  expect_identical(state$factors, "Condition")
  expect_identical(state$contrasts$contrast_name, "Treatment.vs.Control")
  expect_identical(state$contrasts$numerator, "Treatment")
  expect_identical(state$contrasts$denominator, "Control")
})

test_that("runLipidDesignInitialStateShell preserves state hydration", {
  skip_if_not(
    exists("runLipidDesignInitialStateShell", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  run_initial_state_shell <- get("runLipidDesignInitialStateShell", envir = asNamespace("MultiScholaR"))
  state <- list(
    design_matrix = data.frame(
      Run = c("Sample_1", "Sample_2"),
      group = c("Control", "Treatment"),
      stringsAsFactors = FALSE
    ),
    data_cln = list(
      LCMS_Pos = data.frame(
        database_identifier = c("L1", "L2"),
        Sample_1 = c(1, 2),
        Sample_2 = c(3, 4),
        stringsAsFactors = FALSE
      )
    ),
    sample_names = c("Sample_1", "Sample_2"),
    groups = c("Control", "Treatment"),
    factors = c("Condition"),
    contrasts = data.frame(
      contrast_name = "Treatment.vs.Control",
      numerator = "Treatment",
      denominator = "Control",
      stringsAsFactors = FALSE
    ),
    formula = "~ 0 + group"
  )

  captured <- list(
    selectize = list(),
    text = list()
  )
  reactive_val <- function(initial = NULL) {
    value <- initial
    function(new_value) {
      if (missing(new_value)) {
        return(value)
      }
      value <<- new_value
      invisible(NULL)
    }
  }

  design_matrix <- reactive_val()
  data_cln_reactive <- reactive_val()
  sample_names_reactive <- reactive_val()
  groups <- reactive_val()
  factors <- reactive_val()
  contrasts <- reactive_val()
  removed_samples <- reactive_val("stale")

  result <- run_initial_state_shell(
    state = state,
    session = "builder-session",
    designMatrix = design_matrix,
    dataClnReactive = data_cln_reactive,
    sampleNamesReactive = sample_names_reactive,
    groups = groups,
    factors = factors,
    contrasts = contrasts,
    removedSamples = removed_samples,
    updateSelectizeInputFn = function(session, inputId, choices = NULL, selected = NULL) {
      captured$selectize[[length(captured$selectize) + 1L]] <<- list(
        session = session,
        inputId = inputId,
        choices = choices,
        selected = selected
      )
      invisible(NULL)
    },
    updateTextInputFn = function(session, inputId, value = NULL) {
      captured$text[[length(captured$text) + 1L]] <<- list(
        session = session,
        inputId = inputId,
        value = value
      )
      invisible(NULL)
    }
  )

  expect_identical(result, state)
  expect_identical(design_matrix(), state$design_matrix)
  expect_identical(data_cln_reactive(), state$data_cln)
  expect_identical(sample_names_reactive(), state$sample_names)
  expect_identical(groups(), state$groups)
  expect_identical(factors(), state$factors)
  expect_identical(contrasts(), state$contrasts)
  expect_identical(removed_samples(), character(0))
  expect_length(captured$selectize, 5L)
  expect_identical(captured$text, list(list(
    session = "builder-session",
    inputId = "formula_string",
    value = "~ 0 + group"
  )))
})

test_that("registerLipidDesignSummaryOutputShells preserves render handoff", {
  skip_if_not(
    exists("registerLipidDesignSummaryOutputShells", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  register_summary_outputs <- get("registerLipidDesignSummaryOutputShells", envir = asNamespace("MultiScholaR"))
  output <- new.env(parent = emptyenv())
  render_calls <- character()

  result <- register_summary_outputs(
    output = output,
    designMatrix = function() {
      data.frame(
        Run = c("Sample_1", "Sample_2"),
        group = c("Control", "Control"),
        replicates = c(1L, 1L),
        tech_reps = c(1L, 2L),
        stringsAsFactors = FALSE
      )
    },
    removedSamples = function() c("Sample_10", "Sample_2"),
    formulaString = function() "~ 0 + group",
    renderTextFn = function(expr) {
      render_calls <<- c(render_calls, force(expr))
      structure(list(kind = "mock_render", index = length(render_calls)), class = "mock_render")
    },
    reqFn = function(...) invisible(NULL)
  )

  expect_identical(result, output)
  expect_identical(
    render_calls,
    c(
      "Group: Control, Biological Replicate: 1\n  Samples: Sample_1, Sample_2\n  Technical Replicates: 1, 2",
      "Removed 2 sample(s):\nSample_2\nSample_10",
      "Note: Contrasts will use 'group' prefix based on formula: ~ 0 + group"
    )
  )
  expect_s3_class(output$tech_rep_summary, "mock_render")
  expect_s3_class(output$removed_samples_display, "mock_render")
  expect_s3_class(output$contrast_factors_info, "mock_render")
})

test_that("registerLipidDesignAdjacentOutputShells preserves render handoff", {
  skip_if_not(
    exists("registerLipidDesignAdjacentOutputShells", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  register_adjacent_outputs <- get("registerLipidDesignAdjacentOutputShells", envir = asNamespace("MultiScholaR"))
  output <- new.env(parent = emptyenv())
  ui_calls <- list()
  text_calls <- character()

  result <- register_adjacent_outputs(
    output = output,
    session = list(ns = function(id) paste0("ns-", id)),
    factors = function() c("Treatment", "Timepoint"),
    contrasts = function() {
      data.frame(
        contrast_name = "Treatment.vs.Control",
        numerator = "Treatment",
        denominator = "Control",
        stringsAsFactors = FALSE
      )
    },
    formulaString = function() "~ 0 + group",
    samplesToTransform = function() c("Sample_A", "Sample_B"),
    rangeStart = function() 1L,
    rangeEnd = function() 2L,
    selectedRuns = function() c("Sample_A", "Sample_B", "Sample_C"),
    renderUiFn = function(expr) {
      ui_calls[[length(ui_calls) + 1L]] <<- force(expr)
      structure(list(kind = "mock_ui", index = length(ui_calls)), class = "mock_ui")
    },
    renderTextFn = function(expr) {
      text_calls <<- c(text_calls, force(expr))
      structure(list(kind = "mock_text", index = length(text_calls)), class = "mock_text")
    },
    reqFn = function(...) invisible(NULL),
    paragraphFn = function(value) list(kind = "p", value = value),
    tagListFn = function(value) list(kind = "tagList", value = value),
    codeTagFn = function(value) list(kind = "code", value = value),
    extractExperimentFn = function(sampleName, mode, start, end) {
      paste(sampleName, mode, start, end, sep = "|")
    },
    numericInputFn = function(inputId, label, value, min) {
      list(
        kind = "numericInput",
        inputId = inputId,
        label = label,
        value = value,
        min = min
      )
    }
  )

  expect_identical(result, output)
  expect_identical(ui_calls[[1L]], list(kind = "p", value = "Treatment, Timepoint"))
  expect_identical(
    ui_calls[[2L]]$value[[1L]]$value,
    list(kind = "code", value = "Treatment_vs_Control=groupTreatment-groupControl")
  )
  expect_identical(
    ui_calls[[3L]],
    list(kind = "numericInput", inputId = "ns-replicate_start", label = "Starting replicate number for 3 selected samples:", value = 1, min = 1)
  )
  expect_identical(text_calls, "\"Sample_A\" -> \"Sample_A|range|1|2\"")
  expect_s3_class(output$available_factors_display, "mock_ui")
  expect_s3_class(output$defined_contrasts_display, "mock_ui")
  expect_s3_class(output$range_preview, "mock_text")
  expect_s3_class(output$replicate_inputs, "mock_ui")
})

test_that("runLipidDesignSaveResultsShell preserves the empty-design warning path", {
  skip_if_not(
    exists("runLipidDesignSaveResultsShell", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  run_save_results_shell <- get("runLipidDesignSaveResultsShell", envir = asNamespace("MultiScholaR"))
  notifications <- list()
  assigned_result <- "unset"

  result <- run_save_results_shell(
    designMatrix = data.frame(Run = c("Sample_1", "Sample_2"), stringsAsFactors = FALSE),
    currentRemovedSamples = character(0),
    dataCln = list(LCMS_Pos = data.frame(Run = c("Sample_1", "Sample_2"), stringsAsFactors = FALSE)),
    colMap = list(lipid_id_col = "database_identifier", annotation_col = "annotation"),
    contrastData = NULL,
    configList = list(deAnalysisParameters = list(formula_string = "~ old")),
    formulaString = "~ group",
    resultSetter = function(value) {
      assigned_result <<- value
    },
    buildSaveResultsPayloadFn = function(...) NULL,
    showNotificationFn = function(message, type = NULL, duration = NULL) {
      notifications[[length(notifications) + 1L]] <<- list(
        message = message,
        type = type,
        duration = duration
      )
      invisible(NULL)
    }
  )

  expect_null(result)
  expect_identical(assigned_result, "unset")
  expect_identical(
    notifications,
    list(list(
      message = "No samples have been assigned to groups.",
      type = "warning",
      duration = NULL
    ))
  )
})

test_that("runLipidDesignSaveResultsShell preserves the save-results success path", {
  skip_if_not(
    exists("runLipidDesignSaveResultsShell", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  run_save_results_shell <- get("runLipidDesignSaveResultsShell", envir = asNamespace("MultiScholaR"))
  notifications <- list()
  assigned_result <- NULL
  expected_result <- list(
    design_matrix = data.frame(Run = "Sample_1", group = "Control", stringsAsFactors = FALSE),
    data_cln = list(LCMS_Pos = data.frame(database_identifier = "L1", Sample_1 = 10, stringsAsFactors = FALSE)),
    contrasts_tbl = NULL,
    config_list = list(deAnalysisParameters = list(formula_string = "~ 0 + group"))
  )

  result <- run_save_results_shell(
    designMatrix = data.frame(Run = "Sample_1", group = "Control", stringsAsFactors = FALSE),
    currentRemovedSamples = "Sample_2",
    dataCln = list(LCMS_Pos = data.frame(Run = c("Sample_1", "Sample_2"), stringsAsFactors = FALSE)),
    colMap = list(lipid_id_col = "database_identifier", annotation_col = "annotation"),
    contrastData = NULL,
    configList = list(deAnalysisParameters = list(formula_string = "~ old")),
    formulaString = "~ 0 + group",
    resultSetter = function(value) {
      assigned_result <<- value
    },
    buildSaveResultsPayloadFn = function(
        designMatrix,
        currentRemovedSamples,
        dataCln,
        colMap,
        contrastData,
        configList,
        formulaString
    ) {
      expect_identical(designMatrix$Run, "Sample_1")
      expect_identical(currentRemovedSamples, "Sample_2")
      expect_identical(names(dataCln), "LCMS_Pos")
      expect_identical(colMap$lipid_id_col, "database_identifier")
      expect_null(contrastData)
      expect_identical(configList$deAnalysisParameters$formula_string, "~ old")
      expect_identical(formulaString, "~ 0 + group")
      expected_result
    },
    showNotificationFn = function(message, type = NULL, duration = NULL) {
      notifications[[length(notifications) + 1L]] <<- list(
        message = message,
        type = type,
        duration = duration
      )
      invisible(NULL)
    }
  )

  expect_identical(result, expected_result)
  expect_identical(assigned_result, expected_result)
  expect_identical(
    notifications,
    list(list(
      message = "Design saved successfully!",
      type = "message",
      duration = 5
    ))
  )
})

test_that("runLipidDesignResetConfirmationShell preserves the reset flow", {
  skip_if_not(
    exists("runLipidDesignResetConfirmationShell", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  run_reset_shell <- get("runLipidDesignResetConfirmationShell", envir = asNamespace("MultiScholaR"))
  call_log <- list()

  result <- run_reset_shell(
    scope = "formula",
    initialState = list(formula = "~ condition"),
    designMatrix = function(...) invisible(NULL),
    dataClnReactive = function(...) invisible(NULL),
    sampleNamesReactive = function(...) invisible(NULL),
    removedSamples = function(...) invisible(NULL),
    factors = function(...) invisible(NULL),
    groups = function(...) invisible(NULL),
    contrasts = function(...) invisible(NULL),
    session = "builder-session",
    applyResetStateFn = function(
        scope,
        initialState,
        designMatrix,
        dataClnReactive,
        sampleNamesReactive,
        removedSamples,
        factors,
        groups,
        contrasts,
        updateFormulaFn
    ) {
      call_log[[length(call_log) + 1L]] <<- list(
        kind = "applyResetState",
        scope = scope,
        initialState = initialState
      )
      updateFormulaFn("~ condition")
      invisible(NULL)
    },
    removeModalFn = function() {
      call_log[[length(call_log) + 1L]] <<- list(kind = "removeModal")
      invisible(NULL)
    },
    showNotificationFn = function(message, type) {
      call_log[[length(call_log) + 1L]] <<- list(kind = "showNotification", message = message, type = type)
      invisible(NULL)
    },
    updateTextInputFn = function(session, inputId, value) {
      call_log[[length(call_log) + 1L]] <<- list(
        kind = "updateTextInput",
        session = session,
        inputId = inputId,
        value = value
      )
      invisible(NULL)
    }
  )

  expect_null(result)
  expect_identical(call_log[[1L]], list(
    kind = "applyResetState",
    scope = "formula",
    initialState = list(formula = "~ condition")
  ))
  expect_identical(call_log[[2L]], list(
    kind = "updateTextInput",
    session = "builder-session",
    inputId = "formula_string",
    value = "~ condition"
  ))
  expect_identical(call_log[[3L]], list(kind = "removeModal"))
  expect_identical(call_log[[4L]], list(
    kind = "showNotification",
    message = "Reset of formula completed.",
    type = "message"
  ))
})

test_that("applyLipidDesignSampleRenameMap preserves rename propagation", {
  skip_if_not(
    exists("applyLipidDesignSampleRenameMap", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  apply_sample_rename_map <- get("applyLipidDesignSampleRenameMap", envir = asNamespace("MultiScholaR"))

  renamed <- apply_sample_rename_map(
    designMatrix = data.frame(
      Run = c("Sample_A", "Sample_B"),
      group = c("Control", "Treatment"),
      stringsAsFactors = FALSE
    ),
    assayList = list(
      LCMS_Pos = data.frame(
        database_identifier = "L1",
        Sample_A = 10,
        Sample_B = 20,
        stringsAsFactors = FALSE
      )
    ),
    sampleNames = c("Sample_A", "Sample_B"),
    renameMap = c(Sample_A = "Renamed_A")
  )

  expect_identical(renamed$designMatrix$Run, c("Renamed_A", "Sample_B"))
  expect_identical(names(renamed$assayList$LCMS_Pos), c("database_identifier", "Renamed_A", "Sample_B"))
  expect_identical(renamed$sampleNames, c("Renamed_A", "Sample_B"))
})

test_that("buildLipidDesignBulkRenameMap preserves supported transforms", {
  skip_if_not(
    exists("buildLipidDesignBulkRenameMap", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  build_bulk_rename_map <- get("buildLipidDesignBulkRenameMap", envir = asNamespace("MultiScholaR"))

  expect_identical(
    build_bulk_rename_map(
      samplesToTransform = c("Alpha_Beta_Gamma", "Delta_Epsilon"),
      transformMode = "before_underscore",
      extractExperimentFn = function(sampleName, mode, start = NULL, end = NULL) {
        expect_identical(mode, "start")
        paste(sampleName, mode, sep = "|")
      },
      reqFn = function(...) invisible(NULL)
    ),
    c(
      Alpha_Beta_Gamma = "Alpha_Beta_Gamma|start",
      Delta_Epsilon = "Delta_Epsilon|start"
    )
  )

  expect_identical(
    build_bulk_rename_map(
      samplesToTransform = c("Alpha_Beta_Gamma"),
      transformMode = "after_underscore",
      extractExperimentFn = function(sampleName, mode, start = NULL, end = NULL) {
        expect_identical(mode, "end")
        paste(sampleName, mode, sep = "|")
      },
      reqFn = function(...) invisible(NULL)
    ),
    c(Alpha_Beta_Gamma = "Alpha_Beta_Gamma|end")
  )

  expect_identical(
    build_bulk_rename_map(
      samplesToTransform = c("Alpha_Beta_Gamma"),
      transformMode = "range",
      rangeStart = 1L,
      rangeEnd = 2L,
      extractExperimentFn = function(sampleName, mode, start = NULL, end = NULL) {
        expect_identical(mode, "range")
        paste(sampleName, mode, start, end, sep = "|")
      },
      reqFn = function(...) invisible(NULL)
    ),
    c(Alpha_Beta_Gamma = "Alpha_Beta_Gamma|range|1|2")
  )
})

test_that("applyLipidDesignResetState preserves the all-scope reset behavior", {
  skip_if_not(
    exists("applyLipidDesignResetState", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  apply_reset_state <- get("applyLipidDesignResetState", envir = asNamespace("MultiScholaR"))

  reactive_val <- function(initial = NULL) {
    value <- initial
    function(new_value) {
      if (missing(new_value)) {
        return(value)
      }
      value <<- new_value
      invisible(NULL)
    }
  }

  design_matrix <- reactive_val(data.frame(
    Run = c("Sample_1", "Sample_2"),
    factor1 = c("Old", "Old"),
    factor2 = c("Old2", "Old2"),
    factor3 = c("Old3", "Old3"),
    group = c("Control", "Treatment"),
    tech_reps = c(1L, 2L),
    stringsAsFactors = FALSE
  ))
  data_cln_reactive <- reactive_val(list(LCMS_Pos = data.frame(x = 1)))
  sample_names_reactive <- reactive_val(c("Old_1", "Old_2"))
  removed_samples <- reactive_val(c("Drop_1"))
  factors <- reactive_val(c("OldFactor"))
  groups <- reactive_val(c("Control", "Treatment"))
  contrasts <- reactive_val(data.frame(contrast_name = "old", stringsAsFactors = FALSE))
  formula_value <- "unset"

  initial_state <- list(
    design_matrix = data.frame(
      Run = c("Sample_A", "Sample_B"),
      factor1 = c(NA_character_, NA_character_),
      factor2 = c(NA_character_, NA_character_),
      factor3 = c(NA_character_, NA_character_),
      group = c(NA_character_, NA_character_),
      tech_reps = c(NA_integer_, NA_integer_),
      stringsAsFactors = FALSE
    ),
    data_cln = list(LCMS_Pos = data.frame(y = 1)),
    sample_names = c("Sample_A", "Sample_B"),
    factors = c("Condition"),
    groups = character(0),
    contrasts = data.frame(
      contrast_name = "Treatment.vs.Control",
      numerator = "Treatment",
      denominator = "Control",
      stringsAsFactors = FALSE
    ),
    formula = "~ 0 + group"
  )

  result <- apply_reset_state(
    scope = "all",
    initialState = initial_state,
    designMatrix = design_matrix,
    dataClnReactive = data_cln_reactive,
    sampleNamesReactive = sample_names_reactive,
    removedSamples = removed_samples,
    factors = factors,
    groups = groups,
    contrasts = contrasts,
    updateFormulaFn = function(value) {
      formula_value <<- value
      invisible(NULL)
    }
  )

  expect_null(result)
  expect_identical(design_matrix()$Run, c("Sample_A", "Sample_B"))
  expect_true(all(is.na(design_matrix()$factor1)))
  expect_true(all(is.na(design_matrix()$factor2)))
  expect_true(all(is.na(design_matrix()$factor3)))
  expect_true(all(is.na(design_matrix()$group)))
  expect_true(all(is.na(design_matrix()$tech_reps)))
  expect_identical(data_cln_reactive(), initial_state$data_cln)
  expect_identical(sample_names_reactive(), initial_state$sample_names)
  expect_identical(removed_samples(), character(0))
  expect_identical(factors(), initial_state$factors)
  expect_identical(groups(), initial_state$groups)
  expect_identical(contrasts(), initial_state$contrasts)
  expect_identical(formula_value, "~ 0 + group")
})

test_that("getLipidDesignSampleColumns preserves empty, numeric, and explicit selection behavior", {
  skip_if_not(
    exists("getLipidDesignSampleColumns", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  get_sample_columns <- get("getLipidDesignSampleColumns", envir = asNamespace("MultiScholaR"))

  expect_identical(
    get_sample_columns(
      assayList = list(),
      colMap = NULL,
      messageFn = function(...) invisible(NULL)
    ),
    character(0)
  )

  expect_identical(
    get_sample_columns(
      assayList = list(data.frame(
        database_identifier = c("L1", "L2"),
        annotation = c("A1", "A2"),
        Sample_1 = c(10, 20),
        Sample_2 = c(30, 40),
        stringsAsFactors = FALSE
      )),
      colMap = NULL,
      messageFn = function(...) invisible(NULL)
    ),
    c("Sample_1", "Sample_2")
  )

  expect_identical(
    get_sample_columns(
      assayList = list(data.frame(
        database_identifier = c("L1", "L2"),
        annotation = c("A1", "A2"),
        Sample_1 = c(10, 20),
        Sample_2 = c(30, 40),
        stringsAsFactors = FALSE
      )),
      colMap = list(
        lipid_id_col = "database_identifier",
        annotation_col = "annotation",
        sample_columns = c("Sample_2", "Sample_1")
      ),
      messageFn = function(...) invisible(NULL)
    ),
    c("Sample_2", "Sample_1")
  )
})

test_that("buildLipidDesignSaveResults helpers preserve non-group-prefix contrast formatting", {
  skip_if_not(
    exists("buildLipidDesignSaveResultsDataList", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("buildLipidDesignSaveResultsContrastsTable", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  build_data_list <- get("buildLipidDesignSaveResultsDataList", envir = asNamespace("MultiScholaR"))
  build_contrasts_table <- get("buildLipidDesignSaveResultsContrastsTable", envir = asNamespace("MultiScholaR"))

  filtered_data <- build_data_list(
    dataCln = list(
      LCMS_Pos = data.frame(
        database_identifier = c("L1", "L2"),
        annotation = c("A1", "A2"),
        Sample_1 = c(10, 20),
        Sample_2 = c(30, 40),
        stringsAsFactors = FALSE
      )
    ),
    assignedSamples = "Sample_2",
    colMap = list(lipid_id_col = "database_identifier", annotation_col = "annotation")
  )

  expect_identical(
    names(filtered_data$LCMS_Pos),
    c("database_identifier", "annotation", "Sample_2")
  )

  contrast_table <- build_contrasts_table(
    contrastData = data.frame(
      contrast_name = "Treatment.vs.Control",
      numerator = "Treatment",
      denominator = "Control",
      stringsAsFactors = FALSE
    ),
    formulaString = "~ batch + group"
  )

  expect_identical(contrast_table$contrasts, "Treatment-Control")
  expect_identical(contrast_table$friendly_names, "Treatment_vs_Control")
  expect_identical(contrast_table$full_format, "Treatment_vs_Control=Treatment-Control")
})

test_that("registerLipidDesignContrastManagementShells and sample removal preserve observer handoff", {
  skip_if_not(
    exists("registerLipidDesignContrastManagementShells", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("registerLipidDesignSampleRemovalShells", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  register_contrast_shells <- get("registerLipidDesignContrastManagementShells", envir = asNamespace("MultiScholaR"))
  register_sample_removal_shells <- get("registerLipidDesignSampleRemovalShells", envir = asNamespace("MultiScholaR"))

  reactive_val <- function(initial = NULL) {
    value <- initial
    function(new_value) {
      if (missing(new_value)) {
        return(value)
      }
      value <<- new_value
      invisible(NULL)
    }
  }

  observed_events <- character()
  observe_event_fn <- function(eventExpr, handlerExpr, ignoreNULL = FALSE) {
    observed_events <<- c(observed_events, deparse(substitute(eventExpr), nlines = 1L))
    eval.parent(substitute(handlerExpr))
    invisible(NULL)
  }

  contrasts <- reactive_val(data.frame(
    contrast_name = character(),
    numerator = character(),
    denominator = character(),
    stringsAsFactors = FALSE
  ))
  register_contrast_shells(
    input = list(add_contrast = 1, contrast_group1 = "Treatment", contrast_group2 = "Control"),
    contrasts = contrasts,
    observeEventFn = observe_event_fn,
    reqFn = function(...) invisible(NULL)
  )

  expect_identical(contrasts()$contrast_name, "Treatment.vs.Control")

  removal_updates <- list()
  removal_notifications <- list()
  removed_samples <- reactive_val("Existing")
  register_sample_removal_shells(
    input = list(remove_samples = 1, samples_to_remove = c("Sample_1", "Sample_2")),
    removedSamples = removed_samples,
    session = "builder-session",
    observeEventFn = observe_event_fn,
    reqFn = function(...) invisible(NULL),
    updateSelectizeInputFn = function(session, inputId, selected = NULL, choices = NULL) {
      removal_updates[[length(removal_updates) + 1L]] <<- list(
        session = session,
        inputId = inputId,
        selected = selected,
        choices = choices
      )
      invisible(NULL)
    },
    showNotificationFn = function(message, type = NULL) {
      removal_notifications[[length(removal_notifications) + 1L]] <<- list(
        message = message,
        type = type
      )
      invisible(NULL)
    }
  )

  expect_identical(observed_events, c("input$add_contrast", "input$remove_samples"))
  expect_identical(removed_samples(), c("Existing", "Sample_1", "Sample_2"))
  expect_identical(removal_updates, list(list(
    session = "builder-session",
    inputId = "samples_to_remove",
    selected = "",
    choices = NULL
  )))
  expect_identical(removal_notifications, list(list(
    message = "Removed 2 sample(s) from analysis.",
    type = "message"
  )))
})
