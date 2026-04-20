library(testthat)
library(shiny)

makeProtDesignCharacterizationHarness <- function() {
  data_tbl <- shiny::reactive({
    data.frame(
      Run = c("S1", "S2"),
      Intensity = c(10, 20),
      stringsAsFactors = FALSE
    )
  })
  config_list <- shiny::reactive(
    list(deAnalysisParameters = list(formula_string = "~ 0 + group"))
  )
  column_mapping <- shiny::reactive(
    list(quantity_col = "Intensity")
  )

  session <- shiny::MockShinySession$new()

  returned <- shiny::withReactiveDomain(session, {
    mod_prot_design_builder_server(
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

setProtDesignInputs <- function(session, ...) {
  values <- list(...)
  names(values) <- paste0("builder-", names(values))
  do.call(session$setInputs, values)
  shiny:::flushReact()
}

test_that("mod_prot_design_builder_server preserves the empty save-results public behavior", {
  harness <- makeProtDesignCharacterizationHarness()

  expect_null(shiny::isolate(harness$returned()))

  setProtDesignInputs(
    harness$session,
    formula_string = "~ 0 + group",
    save_results = 1
  )

  expect_null(shiny::isolate(harness$returned()))
})

test_that("runProtDesignSaveResultsObserverShell preserves the empty-design warning path", {
  skip_if_not(
    exists("runProtDesignSaveResultsObserverShell", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  save_results_shell <- get("runProtDesignSaveResultsObserverShell", envir = asNamespace("MultiScholaR"))
  notifications <- list()
  assigned_result <- "unset"

  result <- save_results_shell(
    designMatrix = data.frame(Run = c("S1", "S2"), stringsAsFactors = FALSE),
    currentRemovedSamples = character(0),
    dataCln = data.frame(Run = c("S1", "S2"), stringsAsFactors = FALSE),
    contrastData = NULL,
    configList = list(deAnalysisParameters = list(formula_string = "~ old")),
    formulaString = "~ group",
    resultSetter = function(value) {
      assigned_result <<- value
    },
    buildSaveResultsPayload = function(...) NULL,
    showNotification = function(message, type = NULL, duration = NULL) {
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
      message = "No samples have been assigned to groups. Please assign metadata before saving.",
      type = "warning",
      duration = NULL
    ))
  )
})

test_that("runProtDesignSaveResultsObserverShell preserves the save-results success path", {
  skip_if_not(
    exists("runProtDesignSaveResultsObserverShell", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  save_results_shell <- get("runProtDesignSaveResultsObserverShell", envir = asNamespace("MultiScholaR"))
  notifications <- list()
  assigned_result <- NULL
  expected_result <- list(
    design_matrix = data.frame(Run = "S1", group = "GA", stringsAsFactors = FALSE),
    data_cln = data.frame(Run = "S1", value = 10, stringsAsFactors = FALSE)
  )

  result <- save_results_shell(
    designMatrix = data.frame(Run = "S1", stringsAsFactors = FALSE),
    currentRemovedSamples = "S2",
    dataCln = data.frame(Run = c("S1", "S2"), stringsAsFactors = FALSE),
    contrastData = NULL,
    configList = list(deAnalysisParameters = list(formula_string = "~ old")),
    formulaString = "~ 0 + group",
    resultSetter = function(value) {
      assigned_result <<- value
    },
    buildSaveResultsPayload = function(
        designMatrix,
        currentRemovedSamples,
        dataCln,
        contrastData,
        configList,
        formulaString
    ) {
      expect_identical(designMatrix$Run, "S1")
      expect_identical(currentRemovedSamples, "S2")
      expect_identical(dataCln$Run, c("S1", "S2"))
      expect_null(contrastData)
      expect_identical(configList$deAnalysisParameters$formula_string, "~ old")
      expect_identical(formulaString, "~ 0 + group")
      expected_result
    },
    showNotification = function(message, type = NULL, duration = NULL) {
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
      message = "Design saved successfully. You can close this builder.",
      type = "message",
      duration = 5
    ))
  )
})

test_that("runProtDesignResetConfirmationObserverShell preserves the reset flow", {
  skip_if_not(
    exists("runProtDesignResetConfirmationObserverShell", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  reset_confirmation_shell <- get("runProtDesignResetConfirmationObserverShell", envir = asNamespace("MultiScholaR"))
  call_log <- list()
  design_matrix_value <- data.frame(Run = "ResetSample1", stringsAsFactors = FALSE)

  result <- reset_confirmation_shell(
    scope = "formula",
    initialState = list(formula = "~ condition"),
    designMatrix = function(value) {
      if (missing(value)) {
        call_log[[length(call_log) + 1L]] <<- list(kind = "designMatrixGetter")
        return(design_matrix_value)
      }

      call_log[[length(call_log) + 1L]] <<- list(kind = "designMatrixSetter", value = value)
      design_matrix_value <<- value
      invisible(NULL)
    },
    dataClnReactive = function(value) {
      call_log[[length(call_log) + 1L]] <<- list(kind = "dataClnSetter", value = value)
      invisible(NULL)
    },
    removedSamples = function(value) {
      if (missing(value)) {
        call_log[[length(call_log) + 1L]] <<- list(kind = "removedSamplesGetter")
        return(c("Drop1", "Drop2"))
      }

      call_log[[length(call_log) + 1L]] <<- list(kind = "removedSamplesSetter", value = value)
      invisible(NULL)
    },
    factors = function(value) {
      call_log[[length(call_log) + 1L]] <<- list(kind = "factorsSetter", value = value)
      invisible(NULL)
    },
    groups = function(value) {
      call_log[[length(call_log) + 1L]] <<- list(kind = "groupsSetter", value = value)
      invisible(NULL)
    },
    contrasts = function(value) {
      call_log[[length(call_log) + 1L]] <<- list(kind = "contrastsSetter", value = value)
      invisible(NULL)
    },
    session = "builder-session",
    applyResetStateFn = function(
        scope,
        initialState,
        designMatrixGetter,
        designMatrixSetter,
        dataClnSetter,
        currentRemovedSamples,
        removedSamplesSetter,
        factorsSetter,
        groupsSetter,
        contrastsSetter,
        updateFormulaFn
    ) {
      call_log[[length(call_log) + 1L]] <<- list(
        kind = "applyResetState",
        scope = scope,
        initialState = initialState,
        currentRemovedSamples = currentRemovedSamples
      )
      expect_identical(scope, "formula")
      expect_identical(initialState$formula, "~ condition")
      expect_identical(currentRemovedSamples, c("Drop1", "Drop2"))
      expect_identical(designMatrixGetter(), design_matrix_value)
      designMatrixSetter(design_matrix_value)
      dataClnSetter(data.frame(Run = "ResetSample1", stringsAsFactors = FALSE))
      removedSamplesSetter(character(0))
      factorsSetter("condition")
      groupsSetter("GroupA")
      contrastsSetter(data.frame(contrast_name = "A.vs.B", stringsAsFactors = FALSE))
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
  expect_identical(call_log[[1L]]$kind, "removedSamplesGetter")
  expect_identical(call_log[[2L]]$kind, "applyResetState")
  expect_identical(call_log[[3L]]$kind, "designMatrixGetter")
  expect_identical(call_log[[4L]]$kind, "designMatrixSetter")
  expect_identical(call_log[[5L]]$kind, "dataClnSetter")
  expect_identical(call_log[[6L]]$kind, "removedSamplesSetter")
  expect_identical(call_log[[7L]]$kind, "factorsSetter")
  expect_identical(call_log[[8L]]$kind, "groupsSetter")
  expect_identical(call_log[[9L]]$kind, "contrastsSetter")
  expect_identical(call_log[[10L]], list(
    kind = "updateTextInput",
    session = "builder-session",
    inputId = "formula_string",
    value = "~ condition"
  ))
  expect_identical(call_log[[11L]], list(kind = "removeModal"))
  expect_identical(call_log[[12L]], list(
    kind = "showNotification",
    message = "Reset of formula completed.",
    type = "message"
  ))
})
