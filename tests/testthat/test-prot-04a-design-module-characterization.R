# fidelity-coverage-compare: shared
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

skipIfMissingProtDesignBinding <- function(name) {
  skip_if_not(
    exists(name, envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )
}

skipIfMissingProtDesignBindings <- function(...) {
  names <- unlist(list(...), use.names = FALSE)
  missing <- names[!vapply(names, function(name) {
    exists(name, envir = asNamespace("MultiScholaR"), inherits = FALSE)
  }, logical(1))]
  if (length(missing) > 0) {
    skip(paste("requires extracted design bindings:", paste(missing, collapse = ", ")))
  }
}

getProtDesignBinding <- function(name) {
  get(name, envir = asNamespace("MultiScholaR"), inherits = FALSE)
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

test_that("formatProtDesignTechRepSummary preserves empty and grouped summaries", {
  skipIfMissingProtDesignBinding("formatProtDesignTechRepSummary")

  empty_design <- data.frame(
    Run = c("Run1", "Run2"),
    group = c("GA", "GA"),
    replicates = c(1L, 1L),
    tech_reps = c(NA_integer_, NA_integer_),
    stringsAsFactors = FALSE
  )
  grouped_design <- data.frame(
    Run = c("Run1", "Run2", "Run3"),
    group = c("GA", "GA", "GB"),
    replicates = c(1L, 1L, 2L),
    tech_reps = c(1L, 2L, 1L),
    stringsAsFactors = FALSE
  )

  expect_identical(
    formatProtDesignTechRepSummary(empty_design),
    "No technical replicates assigned yet."
  )
  expect_identical(
    formatProtDesignTechRepSummary(grouped_design),
    paste(
      "Group: GA, Biological Replicate: 1",
      "  Samples: Run1, Run2",
      "  Technical Replicates: 1, 2",
      "",
      "Group: GB, Biological Replicate: 2",
      "  Samples: Run3",
      "  Technical Replicates: 1",
      sep = "\n"
    )
  )
})

test_that("createProtDesignBuilderState preserves mutable builder shells", {
  skipIfMissingProtDesignBinding("createProtDesignBuilderState")

  seeds <- list()
  reactive_val_fn <- function(value = NULL) {
    current_value <- value
    seeds[[length(seeds) + 1L]] <<- list(value)

    function(new_value) {
      if (!missing(new_value)) {
        current_value <<- new_value
      }

      current_value
    }
  }

  builder_state <- createProtDesignBuilderState(
    reactiveValFn = reactive_val_fn
  )

  expect_identical(
    names(builder_state),
    c(
      "resultRv",
      "designMatrix",
      "dataClnReactive",
      "groups",
      "factors",
      "contrasts",
      "removedSamples"
    )
  )
  expect_identical(
    seeds,
    list(list(NULL), list(NULL), list(NULL), list(NULL), list(NULL), list(NULL), list(character(0)))
  )
  expect_null(builder_state$resultRv())
  expect_null(builder_state$designMatrix())
  expect_null(builder_state$dataClnReactive())
  expect_null(builder_state$groups())
  expect_null(builder_state$factors())
  expect_null(builder_state$contrasts())
  expect_identical(builder_state$removedSamples(), character(0))
})

test_that("createProtDesignMutableStateShells preserves mutable aliases", {
  skipIfMissingProtDesignBinding("createProtDesignMutableStateShells")

  builder_state <- list(
    resultRv = function(value) value,
    designMatrix = function(value) value,
    dataClnReactive = function(value) value,
    groups = function(value) value,
    factors = function(value) value,
    contrasts = function(value) value,
    removedSamples = function(value) value
  )

  mutable_state <- createProtDesignMutableStateShells(
    createBuilderStateFn = function() builder_state
  )

  expect_identical(names(mutable_state), names(builder_state))
  expect_identical(mutable_state$resultRv, builder_state$resultRv)
  expect_identical(mutable_state$designMatrix, builder_state$designMatrix)
  expect_identical(mutable_state$dataClnReactive, builder_state$dataClnReactive)
  expect_identical(mutable_state$groups, builder_state$groups)
  expect_identical(mutable_state$factors, builder_state$factors)
  expect_identical(mutable_state$contrasts, builder_state$contrasts)
  expect_identical(mutable_state$removedSamples, builder_state$removedSamples)
})

test_that("createProtDesignDataTableProxy preserves the default proxy id", {
  skipIfMissingProtDesignBinding("createProtDesignDataTableProxy")

  proxy_calls <- list()
  proxy <- createProtDesignDataTableProxy(
    dataTableProxyFn = function(outputId) {
      proxy_calls[[length(proxy_calls) + 1L]] <<- outputId
      list(proxyId = outputId)
    }
  )

  expect_identical(proxy_calls, list("data_table"))
  expect_identical(proxy, list(proxyId = "data_table"))
})

test_that("createProtDesignInitialStateReactive preserves deferred bootstrap semantics", {
  skipIfMissingProtDesignBinding("createProtDesignInitialStateReactive")

  data_tbl_calls <- 0L
  config_list_calls <- 0L
  column_mapping_calls <- 0L
  build_initial_state_calls <- 0L
  req_calls <- 0L
  reactive_expr <- NULL
  reactive_env <- NULL

  data_tbl <- function(...) {
    data_tbl_calls <<- data_tbl_calls + 1L
    data.frame(Run = "S1", stringsAsFactors = FALSE)
  }
  config_list <- function(...) {
    config_list_calls <<- config_list_calls + 1L
    list(globalParameters = list(workflow = "design"))
  }
  column_mapping <- function(...) {
    column_mapping_calls <<- column_mapping_calls + 1L
    list(sample_id = "Run")
  }

  initial_state <- createProtDesignInitialStateReactive(
    dataTbl = data_tbl,
    configList = config_list,
    columnMapping = column_mapping,
    buildInitialState = function(dataTbl, configList, columnMapping) {
      build_initial_state_calls <<- build_initial_state_calls + 1L
      list(
        dataTbl = dataTbl,
        configList = configList,
        columnMapping = columnMapping
      )
    },
    reactiveFn = function(expr) {
      reactive_expr <<- substitute(expr)
      reactive_env <<- parent.frame()
      function() eval(reactive_expr, envir = reactive_env)
    },
    reqFn = function(value) {
      req_calls <<- req_calls + 1L
      value
    }
  )

  expect_true(is.function(initial_state))
  expect_identical(data_tbl_calls, 0L)
  expect_identical(config_list_calls, 0L)
  expect_identical(column_mapping_calls, 0L)
  expect_identical(build_initial_state_calls, 0L)
  expect_identical(req_calls, 0L)

  result <- initial_state()

  expect_identical(
    result,
    list(
      dataTbl = data.frame(Run = "S1", stringsAsFactors = FALSE),
      configList = list(globalParameters = list(workflow = "design")),
      columnMapping = list(sample_id = "Run")
    )
  )
  expect_identical(data_tbl_calls, 2L)
  expect_identical(config_list_calls, 2L)
  expect_identical(column_mapping_calls, 1L)
  expect_identical(build_initial_state_calls, 1L)
  expect_identical(req_calls, 2L)
  expect_true(is.language(reactive_expr))
  expect_true(is.environment(reactive_env))
})

test_that("completeProtDesignBuilderServerBootstrap preserves proxy bootstrap fan-out", {
  skipIfMissingProtDesignBinding("completeProtDesignBuilderServerBootstrap")

  call_log <- character()
  fake_input <- list(reset_builder = 1)
  fake_output <- structure(list(), class = "shinyoutput")
  fake_session <- list(ns = function(id) paste0("builder-shell-", id))
  proxy_seed <- structure(list(id = "data_table"), class = "prot_design_proxy")
  config_list <- function(...) list(globalParameters = list(workflow = "design"))

  builder_state <- list(
    initialState = function() "initial-state",
    resultRv = function() "builder-results",
    designMatrix = function() "design-matrix",
    dataClnReactive = function() "clean-data",
    groups = function() c("Control", "Treatment"),
    factors = function() c("Condition", "Batch"),
    contrasts = function() "contrasts",
    removedSamples = function() "Run2"
  )

  result <- completeProtDesignBuilderServerBootstrap(
    input = fake_input,
    output = fake_output,
    session = fake_session,
    builderState = builder_state,
    configList = config_list,
    createDataTableProxy = function() {
      call_log <<- c(call_log, "proxy")
      proxy_seed
    },
    registerBuilderServerShells = function(input, output, session, initialState, proxyDataTable, designMatrix, dataClnReactive, removedSamples, factors, groups, contrasts, configList, resultRv) {
      expect_identical(input, fake_input)
      expect_identical(output, fake_output)
      expect_identical(session, fake_session)
      expect_identical(initialState, builder_state$initialState)
      expect_identical(proxyDataTable, proxy_seed)
      expect_identical(designMatrix, builder_state$designMatrix)
      expect_identical(dataClnReactive, builder_state$dataClnReactive)
      expect_identical(removedSamples, builder_state$removedSamples)
      expect_identical(factors, builder_state$factors)
      expect_identical(groups, builder_state$groups)
      expect_identical(contrasts, builder_state$contrasts)
      expect_identical(configList, config_list)
      expect_identical(resultRv, builder_state$resultRv)
      call_log <<- c(call_log, "builderServerShells")
      invisible(NULL)
    }
  )

  expect_identical(result, builder_state$resultRv)
  expect_identical(call_log, c("proxy", "builderServerShells"))
})

test_that("runProtDesignBuilderServerEntryShell preserves module entry wiring", {
  skipIfMissingProtDesignBinding("runProtDesignBuilderServerEntryShell")

  call_log <- character()
  data_tbl <- function() data.frame(Run = "S1", stringsAsFactors = FALSE)
  config_list <- function() list(globalParameters = list(workflow = "design"))
  column_mapping <- function() list(sample_id = "Run")
  fake_input <- list(reset_builder = 1)
  fake_output <- structure(list(), class = "shinyoutput")
  fake_session <- list(ns = function(id) paste0("builder-entry-", id))
  entry_id <- "builder-entry"

  result <- runProtDesignBuilderServerEntryShell(
    id = entry_id,
    dataTbl = data_tbl,
    configList = config_list,
    columnMapping = column_mapping,
    moduleServer = function(id, module, ...) {
      expect_identical(id, entry_id)
      call_log <<- c(call_log, "moduleServer")
      module(fake_input, fake_output, fake_session)
    },
    runBuilderModuleServerShell = function(input, output, session, dataTbl, configList, columnMapping) {
      expect_identical(input, fake_input)
      expect_identical(output, fake_output)
      expect_identical(session, fake_session)
      expect_identical(dataTbl, data_tbl)
      expect_identical(configList, config_list)
      expect_identical(columnMapping, column_mapping)
      call_log <<- c(call_log, "builderModuleServerShell")
      "builder-results"
    }
  )

  expect_identical(result, "builder-results")
  expect_identical(call_log, c("moduleServer", "builderModuleServerShell"))
})

test_that("buildProtDesignAvailableFactorsDisplay preserves empty and populated displays", {
  skipIfMissingProtDesignBinding("buildProtDesignAvailableFactorsDisplay")

  expect_equal(
    as.character(htmltools::renderTags(buildProtDesignAvailableFactorsDisplay(character(0)))$html),
    "<p>No factors defined yet (use the 'Factors' tab).</p>"
  )
  expect_equal(
    as.character(htmltools::renderTags(buildProtDesignAvailableFactorsDisplay(c("Condition", "Batch", "Timepoint")))$html),
    "<p>Condition, Batch, Timepoint</p>"
  )
})

test_that("applyProtDesignFactorAppendReset preserves unique append semantics", {
  skipIfMissingProtDesignBinding("applyProtDesignFactorAppendReset")

  updated <- applyProtDesignFactorAppendReset(
    currentFactors = c("Condition", "Batch"),
    newFactorInput = "  Timepoint  "
  )
  blank_update <- applyProtDesignFactorAppendReset(
    currentFactors = c("Condition", "Batch"),
    newFactorInput = "   "
  )
  duplicate_update <- applyProtDesignFactorAppendReset(
    currentFactors = c("Condition", "Batch"),
    newFactorInput = "Batch"
  )

  expect_identical(updated$factors, c("Condition", "Batch", "Timepoint"))
  expect_identical(updated$newFactorValue, "")
  expect_identical(blank_update$factors, c("Condition", "Batch"))
  expect_identical(blank_update$newFactorValue, "")
  expect_identical(duplicate_update$factors, c("Condition", "Batch"))
  expect_identical(duplicate_update$newFactorValue, "")
})

test_that("buildProtDesignSaveResults helpers preserve formatted payload construction", {
  skipIfMissingProtDesignBinding("buildProtDesignSaveResultsContrastsTable")
  skipIfMissingProtDesignBinding("buildProtDesignSaveResultsPayload")

  contrast_data <- data.frame(
    contrast_name = c("GA.Control.vs.GA.Elevated", "GB.Control.vs.GB.Elevated"),
    numerator = c("GA_Control", "GB_Control"),
    denominator = c("GA_Elevated", "GB_Elevated"),
    stringsAsFactors = FALSE
  )
  design_matrix <- data.frame(
    Run = c("S1", "S2", "S3"),
    group = c("GA", "GB", "GC"),
    replicates = c(1L, 2L, 3L),
    stringsAsFactors = FALSE
  )
  data_cln <- data.frame(
    Run = c("S1", "S1", "S2", "S3"),
    value = c(10, 11, 20, 30),
    stringsAsFactors = FALSE
  )
  config_list <- list(
    deAnalysisParameters = list(formula_string = "~ old")
  )

  expect_null(
    buildProtDesignSaveResultsContrastsTable(
      contrastData = NULL,
      formulaString = "~ group"
    )
  )
  expect_identical(
    buildProtDesignSaveResultsContrastsTable(
      contrastData = contrast_data,
      formulaString = "~ 0 + group"
    ),
    data.frame(
      contrasts = c(
        "groupGA_Control-groupGA_Elevated",
        "groupGB_Control-groupGB_Elevated"
      ),
      friendly_names = c(
        "GA_Control_vs_GA_Elevated",
        "GB_Control_vs_GB_Elevated"
      ),
      full_format = c(
        "GA_Control_vs_GA_Elevated=groupGA_Control-groupGA_Elevated",
        "GB_Control_vs_GB_Elevated=groupGB_Control-groupGB_Elevated"
      ),
      stringsAsFactors = FALSE
    )
  )
  expect_null(
    buildProtDesignSaveResultsPayload(
      designMatrix = data.frame(
        Run = c("S1", "S2"),
        group = c(NA_character_, ""),
        replicates = c(1L, 2L),
        stringsAsFactors = FALSE
      ),
      currentRemovedSamples = character(0),
      dataCln = data.frame(Run = c("S1", "S2"), value = c(10, 20), stringsAsFactors = FALSE),
      contrastData = NULL,
      configList = list(deAnalysisParameters = list(formula_string = "~ old")),
      formulaString = "~ group"
    )
  )

  result <- buildProtDesignSaveResultsPayload(
    designMatrix = design_matrix,
    currentRemovedSamples = "S2",
    dataCln = data_cln,
    contrastData = data.frame(
      contrast_name = "GA.vs.GC",
      numerator = "GA",
      denominator = "GC",
      stringsAsFactors = FALSE
    ),
    configList = config_list,
    formulaString = "~ 0 + group"
  )

  expect_identical(result$design_matrix$Run, c("S1", "S3"))
  expect_identical(result$design_matrix$tech_rep_group, c("GA_1", "GC_3"))
  expect_identical(result$data_cln$Run, c("S1", "S1", "S3"))
  expect_identical(result$config_list$deAnalysisParameters$formula_string, "~ 0 + group")
  expect_identical(
    result$contrasts_tbl$full_format,
    "GA_vs_GC=groupGA-groupGC"
  )
})

test_that("showProtDesignResetConfirmationModal preserves the modal shell", {
  skipIfMissingProtDesignBinding("showProtDesignResetConfirmationModal")

  modal <- showProtDesignResetConfirmationModal(
    resetScope = "removed_samples",
    nsFn = function(id) paste0("builder-", id),
    showModalFn = function(modal) modal,
    modalDialogFn = function(title, body, footer, easyClose) {
      list(
        title = title,
        body = body,
        footer = footer,
        easyClose = easyClose
      )
    },
    htmlFn = function(text) text,
    tagListFn = function(...) list(...),
    modalButtonFn = function(label) list(type = "modalButton", label = label),
    actionButtonFn = function(inputId, label, class = NULL) {
      list(
        type = "actionButton",
        inputId = inputId,
        label = label,
        class = class
      )
    }
  )

  expect_identical(modal$title, "Confirm Reset")
  expect_match(modal$body, "removed_samples", fixed = TRUE)
  expect_true(modal$easyClose)
  expect_identical(modal$footer[[1]], list(type = "modalButton", label = "Cancel"))
  expect_identical(
    modal$footer[[2]],
    list(
      type = "actionButton",
      inputId = "builder-confirm_reset",
      label = "Reset",
      class = "btn-danger"
    )
  )
})

test_that("removed-sample and contrast formatting helpers preserve display behavior", {
  skipIfMissingProtDesignBinding("applyProtDesignRemovedSamplesUpdate")
  skipIfMissingProtDesignBinding("formatProtDesignRemovedSamplesDisplay")
  skipIfMissingProtDesignBinding("formatProtDesignContrastFactorsInfo")

  expect_identical(
    applyProtDesignRemovedSamplesUpdate(
      currentRemovedSamples = c("Sample1", "Sample3"),
      selectedSamples = c("Sample2", "Sample3", "Sample4")
    ),
    c("Sample1", "Sample3", "Sample2", "Sample4")
  )
  expect_identical(
    applyProtDesignRemovedSamplesUpdate(
      currentRemovedSamples = c("Sample1", "Sample2"),
      reset = TRUE
    ),
    character(0)
  )
  expect_identical(
    formatProtDesignRemovedSamplesDisplay(character(0)),
    "No samples have been removed."
  )
  expect_identical(
    formatProtDesignRemovedSamplesDisplay(c("Sample10", "Sample2", "Sample1")),
    paste(
      "Removed 3 sample(s):",
      paste(c("Sample1", "Sample2", "Sample10"), collapse = "\n"),
      sep = "\n"
    )
  )
  expect_identical(
    formatProtDesignContrastFactorsInfo("~ 0 + group"),
    paste(
      "Note: Contrasts will use 'group' prefix (e.g., groupGA_Control-groupGA_Elevated)",
      "based on current formula: ~ 0 + group",
      sep = "\n"
    )
  )
  expect_identical(
    formatProtDesignContrastFactorsInfo("~ condition + batch"),
    "Note: Contrasts will use group names as-is (e.g., GA_Control-GA_Elevated)"
  )
})

test_that("rename helper updates preserve design and data table rewrites", {
  skipIfMissingProtDesignBinding("applyProtDesignBulkRenameUpdates")
  skipIfMissingProtDesignBinding("applyProtDesignSingleRenameUpdate")

  design_matrix <- data.frame(
    Run = c("SampleA_T1_Rep1", "SampleB_T2_Rep1", "SampleC_T3_Rep1"),
    group = c("G1", "G2", "G3"),
    stringsAsFactors = FALSE
  )
  data_cln <- data.frame(
    Run = c("SampleA_T1_Rep1", "SampleB_T2_Rep1", "SampleA_T1_Rep1"),
    intensity = c(10, 20, 30),
    stringsAsFactors = FALSE
  )

  updated_bulk <- applyProtDesignBulkRenameUpdates(
    designMatrix = design_matrix,
    dataCln = data_cln,
    selectedSamples = c("SampleA_T1_Rep1", "SampleB_T2_Rep1"),
    newNames = c("SampleA", "SampleB")
  )
  updated_single <- applyProtDesignSingleRenameUpdate(
    designMatrix = design_matrix,
    dataCln = data_cln,
    originalName = "SampleA_T1_Rep1",
    newName = "SampleA"
  )

  expect_identical(
    updated_bulk$designMatrix$Run,
    c("SampleA", "SampleB", "SampleC_T3_Rep1")
  )
  expect_identical(updated_bulk$dataCln$Run, c("SampleA", "SampleB", "SampleA"))
  expect_identical(
    updated_single$designMatrix$Run,
    c("SampleA", "SampleB_T2_Rep1", "SampleC_T3_Rep1")
  )
  expect_identical(updated_single$dataCln$Run, c("SampleA", "SampleB_T2_Rep1", "SampleA"))
})

test_that("buildProtDesignInitialState preserves builder bootstrap defaults", {
  skip_if_not(
    exists("buildProtDesignInitialState", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  initial_state <- buildProtDesignInitialState(
    dataTbl = data.frame(
      Run = c("Run10", "Run2", "Run1", "Run2"),
      norm_qty = c("10.5", "11.5", "12.5", "13.5"),
      raw_qty = c("100", "101", "102", "103"),
      quantity = c("1000", "1001", "1002", "1003"),
      stringsAsFactors = FALSE
    ),
    configList = list(
      deAnalysisParameters = list(formula_string = "~ 0 + group")
    ),
    columnMapping = list(
      norm_quantity_col = "norm_qty",
      raw_quantity_col = "raw_qty",
      quantity_col = "quantity"
    )
  )

  expect_identical(initial_state$design_matrix$Run, c("Run1", "Run2", "Run10"))
  expect_true(is.double(initial_state$data_cln$norm_qty))
  expect_true(is.double(initial_state$data_cln$raw_qty))
  expect_true(is.double(initial_state$data_cln$quantity))
  expect_identical(initial_state$formula, "~ 0 + group")
})

test_that("initializeProtDesignBuilderServerState preserves builder state wiring order", {
  skip_if_not(
    exists("initializeProtDesignBuilderServerState", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  call_log <- character()
  fake_input <- list(reset_builder = 1)
  fake_session <- list(ns = function(id) paste0("builder-shell-", id))
  initial_state_seed <- NULL

  data_tbl <- function(...) data.frame(Run = "S1", stringsAsFactors = FALSE)
  config_list <- function(...) list(globalParameters = list(workflow = "design"))
  column_mapping <- function(...) list(sample_id = "Run")

  mutable_state <- list(
    resultRv = function() "builder-results",
    designMatrix = function() "design-matrix",
    dataClnReactive = function() "clean-data",
    groups = function() c("Control", "Treatment"),
    factors = function() c("Condition", "Batch"),
    contrasts = function() "contrasts",
    removedSamples = function() "Run2"
  )

  builder_state <- initializeProtDesignBuilderServerState(
    input = fake_input,
    session = fake_session,
    dataTbl = data_tbl,
    configList = config_list,
    columnMapping = column_mapping,
    buildInitialState = function(...) list(...),
    createInitialStateReactive = function(dataTbl, configList, columnMapping, buildInitialState) {
      expect_identical(dataTbl, data_tbl)
      expect_identical(configList, config_list)
      expect_identical(columnMapping, column_mapping)
      expect_true(is.function(buildInitialState))
      initial_state_seed <<- function() "initial-state"
      call_log <<- c(call_log, "initialStateReactive")
      initial_state_seed
    },
    createMutableStateShells = function() {
      call_log <<- c(call_log, "mutableState")
      mutable_state
    },
    registerInitialStateSyncObserver = function(...) {
      call_log <<- c(call_log, "initialStateSync")
      invisible("observer")
    }
  )

  expect_identical(call_log, c("initialStateReactive", "mutableState", "initialStateSync"))
  expect_identical(builder_state$initialState, initial_state_seed)
  expect_identical(builder_state$resultRv, mutable_state$resultRv)
})

test_that("registerProtDesignBuilderServerShells preserves top-level registration fan-out", {
  skip_if_not(
    exists("registerProtDesignBuilderServerShells", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  call_log <- list()
  input_values <- list(reset_builder = 1)
  output_env <- new.env(parent = emptyenv())
  session_obj <- list(ns = function(id) paste0("builder-", id))
  initial_state <- function() "initial-state"
  proxy_data_table <- list(proxy = "data-table")
  design_matrix <- function() "design-matrix"
  data_cln_reactive <- function() "clean-data"
  removed_samples <- function() "Run2"
  factors_rv <- function() c("Condition", "Batch")
  groups_rv <- function() c("Control", "Treatment")
  contrasts_rv <- function() "contrasts"
  config_list <- function() list(formula = "~ 0 + group")
  result_rv <- function(value) value

  result <- registerProtDesignBuilderServerShells(
    input = input_values,
    output = output_env,
    session = session_obj,
    initialState = initial_state,
    proxyDataTable = proxy_data_table,
    designMatrix = design_matrix,
    dataClnReactive = data_cln_reactive,
    removedSamples = removed_samples,
    factors = factors_rv,
    groups = groups_rv,
    contrasts = contrasts_rv,
    configList = config_list,
    resultRv = result_rv,
    registerInputSyncObserverShells = function(...) call_log[[length(call_log) + 1L]] <<- "inputSync",
    registerRenderOutputShells = function(...) call_log[[length(call_log) + 1L]] <<- "renderOutputs",
    registerEventObserverShells = function(...) call_log[[length(call_log) + 1L]] <<- "eventShells"
  )

  expect_null(result)
  expect_identical(call_log, list("inputSync", "renderOutputs", "eventShells"))
})

test_that("runProtDesignBuilderModuleServerShell preserves bootstrap ordering", {
  skip_if_not(
    exists("runProtDesignBuilderModuleServerShell", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  call_log <- character()
  fake_input <- list(reset_builder = 1)
  fake_output <- structure(list(), class = "shinyoutput")
  fake_session <- list(ns = function(id) paste0("builder-shell-", id))
  data_tbl <- function(...) data.frame(Run = "S1", stringsAsFactors = FALSE)
  config_list <- function(...) list(globalParameters = list(workflow = "design"))
  column_mapping <- function(...) list(sample_id = "Run")

  builder_state <- list(
    initialState = function() "initial-state",
    resultRv = function() "builder-results",
    designMatrix = function() "design-matrix",
    dataClnReactive = function() "clean-data",
    groups = function() c("Control", "Treatment"),
    factors = function() c("Condition", "Batch"),
    contrasts = function() "contrasts",
    removedSamples = function() "Run2"
  )

  result <- runProtDesignBuilderModuleServerShell(
    input = fake_input,
    output = fake_output,
    session = fake_session,
    dataTbl = data_tbl,
    configList = config_list,
    columnMapping = column_mapping,
    initializeBuilderServerState = function(...) {
      call_log <<- c(call_log, "initialize")
      builder_state
    },
    completeBuilderServerBootstrap = function(...) {
      call_log <<- c(call_log, "complete")
      "builder-results"
    }
  )

  expect_identical(call_log, c("initialize", "complete"))
  expect_identical(result, "builder-results")
})

test_that("buildProtDesignDefinedContrastsDisplay preserves formatted grouped contrasts", {
  skip_if_not(
    exists("buildProtDesignDefinedContrastsDisplay", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  contrast_data <- data.frame(
    contrast_name = c("GA.Control.vs.GA.Elevated", "GB.Control.vs.GB.Elevated"),
    numerator = c("GA_Control", "GB_Control"),
    denominator = c("GA_Elevated", "GB_Elevated"),
    stringsAsFactors = FALSE
  )

  display <- buildProtDesignDefinedContrastsDisplay(
    contrastData = contrast_data,
    formulaString = "~ 0 + group"
  )
  html <- htmltools::renderTags(display)$html

  expect_match(html, "GA_Control_vs_GA_Elevated=groupGA_Control-groupGA_Elevated", fixed = TRUE)
  expect_match(html, "GB_Control_vs_GB_Elevated=groupGB_Control-groupGB_Elevated", fixed = TRUE)
})

test_that("transformProtDesignSampleNames preserves supported rename modes", {
  skip_if_not(
    exists("transformProtDesignSampleNames", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  extracted_calls <- list()
  req_calls <- list()
  format_optional <- function(value) if (is.null(value)) "" else as.character(value)

  extract_fn <- function(sampleName, mode, start = NULL, end = NULL) {
    extracted_calls[[length(extracted_calls) + 1L]] <<- list(
      sampleName = sampleName,
      mode = mode,
      start = start,
      end = end
    )
    paste(sampleName, mode, format_optional(start), format_optional(end))
  }

  req_fn <- function(...) {
    req_calls[[length(req_calls) + 1L]] <<- list(...)
    invisible(NULL)
  }

  expect_identical(
    transformProtDesignSampleNames(
      selectedSamples = c("SampleA_T1_Rep1", "SampleB_T2_Rep1"),
      transformMode = "range",
      rangeStart = 2,
      rangeEnd = 2,
      extractExperimentFn = extract_fn,
      reqFn = req_fn
    ),
    c("SampleA_T1_Rep1 range 2 2", "SampleB_T2_Rep1 range 2 2")
  )
  expect_length(req_calls, 2)
  expect_identical(
    transformProtDesignSampleNames(
      selectedSamples = "SampleA_T1_Rep1",
      transformMode = "before_underscore",
      extractExperimentFn = extract_fn
    ),
    "SampleA_T1_Rep1 start  "
  )
  expect_identical(
    transformProtDesignSampleNames(
      selectedSamples = "SampleA_T1_Rep1",
      transformMode = "after_underscore",
      extractExperimentFn = extract_fn
    ),
    "SampleA_T1_Rep1 end  "
  )

  expect_identical(extracted_calls[[1]], list(
    sampleName = "SampleA_T1_Rep1",
    mode = "range",
    start = 2,
    end = 2
  ))
  expect_identical(extracted_calls[[4]], list(
    sampleName = "SampleA_T1_Rep1",
    mode = "end",
    start = NULL,
    end = NULL
  ))
})

test_that("transformProtDesignSampleNames rejects unsupported transform modes", {
  skip_if_not(
    exists("transformProtDesignSampleNames", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  expect_error(
    transformProtDesignSampleNames(
      selectedSamples = "SampleA_T1_Rep1",
      transformMode = "invalid",
      extractExperimentFn = function(...) "ignored"
    ),
    "Unsupported transform mode: invalid"
  )
})

test_that("prot design import artifact helpers preserve sidecar and table loading behavior", {
  skipIfMissingProtDesignBindings(
    "readProtDesignImportedContrasts",
    "resolveProtDesignImportArtifacts",
    "loadProtDesignImportedConfigAndTables",
    "hydrateProtDesignImportedUniprotSidecar",
    "hydrateProtDesignImportedFastaSidecar"
  )

  read_contrasts <- getProtDesignBinding("readProtDesignImportedContrasts")
  resolve_artifacts <- getProtDesignBinding("resolveProtDesignImportArtifacts")
  load_artifacts <- getProtDesignBinding("loadProtDesignImportedConfigAndTables")
  hydrate_uniprot <- getProtDesignBinding("hydrateProtDesignImportedUniprotSidecar")
  hydrate_fasta <- getProtDesignBinding("hydrateProtDesignImportedFastaSidecar")

  cleanup_globals <- intersect(c("contrasts_tbl", "config_list", "uniprot_dat_cln", "aa_seq_tbl_final"), ls(envir = .GlobalEnv))
  withr::defer(rm(list = cleanup_globals, envir = .GlobalEnv), teardown_env())

  import_dir <- tempfile("prot-design-import-")
  source_dir <- tempfile("prot-design-source-")
  dir.create(import_dir, recursive = TRUE)
  dir.create(source_dir, recursive = TRUE)

  contrast_file <- file.path(import_dir, "contrast_strings.tab")
  writeLines(c("groupA-groupB", "groupC-groupD"), contrast_file)
  contrasts <- read_contrasts(contrast_file)
  expect_identical(contrasts$friendly_names, c("A_vs_B", "C_vs_D"))
  expect_null(read_contrasts(file.path(import_dir, "missing.tab")))

  design_file <- file.path(import_dir, "design_matrix.tab")
  data_file <- file.path(import_dir, "data_cln.tab")
  writeLines("Run\tgroup\nS1\tA", design_file)
  writeLines("Run\tProtein.Ids\tAbundance\nS1\tP1\t10", data_file)
  fasta_file <- file.path(import_dir, "proteome.fasta")
  writeLines(">P1\nAAAA", fasta_file)

  resolved <- resolve_artifacts(import_dir)
  expect_true(resolved$ok)
  expect_identical(basename(resolved$fastaPath), "proteome.fasta")
  missing_resolved <- resolve_artifacts(tempfile("missing-design-"), listFiles = function(...) character(), fileExists = function(path) FALSE)
  expect_false(missing_resolved$ok)

  writeLines("[globalParameters]\nworkflow_type = DIA", file.path(source_dir, "config.ini"))
  workflow_data <- new.env(parent = emptyenv())
  loaded <- load_artifacts(
    workflowData = workflow_data,
    experimentPaths = list(source_dir = source_dir),
    designFile = design_file,
    dataClnFile = data_file,
    contrastFile = contrast_file,
    readConfig = function(file) list(globalParameters = list(workflow_type = "DIA")),
    readTabular = function(path) data.frame(source = basename(path), stringsAsFactors = FALSE),
    showNotification = function(...) invisible(NULL),
    assignFn = function(name, value, envir) assign(name, value, envir = envir)
  )
  expect_identical(workflow_data$config_list$globalParameters$workflow_type, "DIA")
  expect_identical(loaded$importedDesign$source, "design_matrix.tab")
  expect_identical(loaded$importedDataCln$source, "data_cln.tab")
  expect_identical(loaded$importedContrasts$friendly_names, c("A_vs_B", "C_vs_D"))

  uniprot_rows <- data.frame(Protein.Ids = "P1", Gene = "GENE1", stringsAsFactors = FALSE)
  saveRDS(uniprot_rows, file.path(import_dir, "uniprot_dat_cln.RDS"))
  fasta_rows <- data.frame(Protein.Ids = "P1", Sequence = "AAAA", stringsAsFactors = FALSE)
  fasta_metadata <- list(fasta_format = "mock", num_sequences = 1)
  saveRDS(fasta_rows, file.path(import_dir, "aa_seq_tbl_final.RDS"))
  saveRDS(fasta_metadata, file.path(import_dir, "fasta_metadata.RDS"))

  session <- shiny::MockShinySession$new()
  shiny::withReactiveDomain(session, {
    expect_identical(
      hydrate_uniprot(workflow_data, importPath = import_dir, sourceDir = source_dir),
      uniprot_rows
    )
    expect_identical(
      hydrate_fasta(workflow_data, importPath = import_dir, sourceDir = source_dir),
      fasta_rows
    )
  })
  expect_true(file.exists(file.path(source_dir, "uniprot_dat_cln.RDS")))
  expect_true(file.exists(file.path(source_dir, "aa_seq_tbl_final.RDS")))
  expect_identical(workflow_data$fasta_metadata, fasta_metadata)
})

test_that("prot design import helpers preserve fallback and bootstrap edge branches", {
  skipIfMissingProtDesignBindings(
    "readProtDesignImportedContrasts",
    "resolveProtDesignImportArtifacts",
    "loadProtDesignImportedConfigAndTables",
    "hydrateProtDesignImportedUniprotSidecar",
    "hydrateProtDesignImportedFastaSidecar"
  )

  read_contrasts <- getProtDesignBinding("readProtDesignImportedContrasts")
  resolve_artifacts <- getProtDesignBinding("resolveProtDesignImportArtifacts")
  load_artifacts <- getProtDesignBinding("loadProtDesignImportedConfigAndTables")
  hydrate_uniprot <- getProtDesignBinding("hydrateProtDesignImportedUniprotSidecar")
  hydrate_fasta <- getProtDesignBinding("hydrateProtDesignImportedFastaSidecar")

  cleanup_globals <- intersect(c("config_list", "uniprot_dat_cln", "aa_seq_tbl_final"), ls(envir = .GlobalEnv))
  withr::defer(rm(list = cleanup_globals, envir = .GlobalEnv), teardown_env())

  import_dir <- tempfile("prot-design-import-fallback-")
  source_dir <- tempfile("prot-design-source-fallback-")
  results_dir <- tempfile("prot-design-results-fallback-")
  dir.create(import_dir, recursive = TRUE)
  dir.create(source_dir, recursive = TRUE)
  dir.create(results_dir, recursive = TRUE)

  empty_contrasts <- file.path(import_dir, "empty-contrast_strings.tab")
  file.create(empty_contrasts)
  expect_null(read_contrasts(empty_contrasts))

  design_file <- file.path(import_dir, "design_matrix.tab")
  data_file <- file.path(import_dir, "data_cln.tab")
  contrast_file <- file.path(import_dir, "contrast_strings.tab")
  writeLines("Run\tgroup\nS1\tA", design_file)
  writeLines("Run\tProtein.Ids\tAbundance\nS1\tP1\t10", data_file)
  writeLines("groupA-groupB", contrast_file)

  selected_fasta <- file.path(import_dir, "selected.fa")
  writeLines(">P1\nAAAA", selected_fasta)
  selected_resolved <- resolve_artifacts(import_dir, selectedFastaPath = selected_fasta)
  expect_true(selected_resolved$ok)
  expect_identical(selected_resolved$fastaPath, selected_fasta)

  uniprot_rows <- data.frame(Protein.Ids = "P1", Gene = "GENE1", stringsAsFactors = FALSE)
  fasta_rows <- data.frame(Protein.Ids = "P1", Sequence = "AAAA", stringsAsFactors = FALSE)
  fasta_metadata <- list(fasta_format = "mock", num_sequences = 1)
  saveRDS(uniprot_rows, file.path(source_dir, "uniprot_dat_cln.RDS"))
  saveRDS(fasta_rows, file.path(source_dir, "aa_seq_tbl_final.RDS"))
  saveRDS(fasta_metadata, file.path(source_dir, "fasta_metadata.RDS"))

  session <- shiny::MockShinySession$new()
  source_only_workflow <- new.env(parent = emptyenv())
  shiny::withReactiveDomain(session, {
    expect_identical(
      hydrate_uniprot(source_only_workflow, importPath = import_dir, sourceDir = source_dir),
      uniprot_rows
    )
    expect_identical(
      hydrate_fasta(source_only_workflow, importPath = import_dir, sourceDir = source_dir),
      fasta_rows
    )
  })
  expect_identical(source_only_workflow$fasta_metadata, fasta_metadata)

  no_sidecar_source <- tempfile("prot-design-no-sidecars-source-")
  no_sidecar_import <- tempfile("prot-design-no-sidecars-import-")
  dir.create(no_sidecar_source, recursive = TRUE)
  dir.create(no_sidecar_import, recursive = TRUE)
  no_sidecar_workflow <- new.env(parent = emptyenv())
  shiny::withReactiveDomain(session, {
    expect_null(hydrate_uniprot(no_sidecar_workflow, importPath = no_sidecar_import, sourceDir = no_sidecar_source))
    expect_null(hydrate_fasta(no_sidecar_workflow, importPath = no_sidecar_import, sourceDir = no_sidecar_source))
  })
  expect_null(no_sidecar_workflow$uniprot_dat_cln)
  expect_null(no_sidecar_workflow$aa_seq_tbl_final)

  fasta_process_workflow <- new.env(parent = emptyenv())
  processed_fasta <- data.frame(Protein.Ids = "P2", Sequence = "BBBB", stringsAsFactors = FALSE)
  processed_metadata <- list(fasta_format = "processed", num_sequences = 1)
  local_mocked_bindings(
    processFastaFile = function(fasta_file_path, ...) {
      if (basename(fasta_file_path) == "bad.fasta") {
        stop("bad fasta")
      }
      list(aa_seq_tbl_final = processed_fasta, fasta_metadata = processed_metadata)
    },
    .env = asNamespace("MultiScholaR")
  )
  shiny::withReactiveDomain(session, {
    expect_identical(
      hydrate_fasta(
        fasta_process_workflow,
        importPath = no_sidecar_import,
        sourceDir = no_sidecar_source,
        resultsDir = results_dir,
        fastaPath = file.path(no_sidecar_import, "good.fasta"),
        organismName = "Homo sapiens"
      ),
      processed_fasta
    )
  })
  expect_true(file.exists(file.path(no_sidecar_source, "aa_seq_tbl_final.RDS")))
  expect_true(file.exists(file.path(no_sidecar_source, "fasta_metadata.RDS")))

  fasta_error_workflow <- new.env(parent = emptyenv())
  shiny::withReactiveDomain(session, {
    expect_null(hydrate_fasta(
      fasta_error_workflow,
      importPath = no_sidecar_import,
      sourceDir = tempfile("prot-design-fasta-error-source-"),
      fastaPath = file.path(no_sidecar_import, "bad.fasta")
    ))
  })
  expect_null(fasta_error_workflow$aa_seq_tbl_final)

  config_design <- file.path(import_dir, "design_matrix.tab")
  config_data <- file.path(import_dir, "data_cln.tab")
  read_tabular <- function(path) data.frame(source = basename(path), stringsAsFactors = FALSE)
  read_config <- function(file) list(globalParameters = list(workflow_type = "DIA"), file = file)

  package_config <- tempfile("package-config-", fileext = ".ini")
  writeLines("[globalParameters]\nworkflow_type = DIA", package_config)
  package_source <- tempfile("prot-design-package-config-source-")
  dir.create(package_source, recursive = TRUE)
  package_notifications <- character()
  package_loaded <- load_artifacts(
    workflowData = new.env(parent = emptyenv()),
    experimentPaths = list(source_dir = package_source),
    designFile = config_design,
    dataClnFile = config_data,
    contrastFile = contrast_file,
    readConfig = read_config,
    readTabular = read_tabular,
    systemFileFn = function(...) package_config,
    fileCopy = function(from, to) file.copy(from, to),
    showNotification = function(message, ...) package_notifications <<- c(package_notifications, message),
    assignFn = function(...) invisible(NULL)
  )
  expect_identical(package_loaded$importedDesign$source, "design_matrix.tab")
  expect_true(file.exists(file.path(package_source, "config.ini")))
  expect_true(any(grepl("default config.ini", package_notifications)))

  download_source <- tempfile("prot-design-download-config-source-")
  dir.create(download_source, recursive = TRUE)
  download_dest <- NULL
  download_loaded <- load_artifacts(
    workflowData = new.env(parent = emptyenv()),
    experimentPaths = list(source_dir = download_source),
    designFile = config_design,
    dataClnFile = config_data,
    contrastFile = contrast_file,
    readConfig = read_config,
    readTabular = read_tabular,
    systemFileFn = function(...) "",
    downloadFile = function(url, destfile, quiet = TRUE) {
      download_dest <<- destfile
      writeLines("[globalParameters]\nworkflow_type = DIA", destfile)
    },
    showNotification = function(...) invisible(NULL),
    assignFn = function(...) invisible(NULL)
  )
  expect_identical(download_loaded$importedDataCln$source, "data_cln.tab")
  expect_identical(download_dest, file.path(download_source, "config.ini"))

  error_source <- tempfile("prot-design-error-config-source-")
  dir.create(error_source, recursive = TRUE)
  expect_error(
    load_artifacts(
      workflowData = new.env(parent = emptyenv()),
      experimentPaths = list(source_dir = error_source),
      designFile = config_design,
      dataClnFile = config_data,
      contrastFile = contrast_file,
      readConfig = read_config,
      readTabular = read_tabular,
      systemFileFn = function(...) "",
      downloadFile = function(...) stop("offline"),
      showNotification = function(...) invisible(NULL),
      assignFn = function(...) invisible(NULL)
    ),
    "Could not obtain a configuration file"
  )
})

test_that("prot design state checkpoint helpers preserve DIA/TMT state and post-checkpoint behavior", {
  skipIfMissingProtDesignBindings(
    "buildProtDesignStateCheckpoint",
    "completeProtDesignPostCheckpoint"
  )

  build_checkpoint <- getProtDesignBinding("buildProtDesignStateCheckpoint")
  complete_checkpoint <- getProtDesignBinding("completeProtDesignPostCheckpoint")

  saved_states <- list()
  state_manager <- new.env(parent = emptyenv())
  state_manager$saveState <- function(...) {
    saved_states[[length(saved_states) + 1L]] <<- list(...)
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- state_manager
  workflow_data$config_list <- list(globalParameters = list(workflow_type = "DIA"))
  workflow_data$design_matrix <- data.frame(Run = c("S1", "S2"), group = c("A", "B"), replicates = c(1, 1))
  workflow_data$data_cln <- data.frame(
    Precursor.Id = c("p1", "p2"),
    Precursor.Charge = c(2, 2),
    Precursor.Quantity = c(10, 20),
    Stripped.Sequence = c("AAA", "BBB"),
    Run = c("S1", "S2")
  )

  local_mocked_bindings(
    PeptideQuantitativeDataDiann = function(...) structure(list(...), class = "mock_peptide_s4"),
    ProteinQuantitativeData = function(...) structure(list(...), class = "mock_protein_s4"),
    .capture_checkpoint = function(...) invisible(NULL),
    .env = asNamespace("MultiScholaR")
  )

  dia_state <- build_checkpoint(
    workflowData = workflow_data,
    workflowType = "DIA",
    actionLabel = "Import"
  )
  expect_identical(dia_state, "raw_data_s4")
  expect_identical(saved_states[[1L]]$state_name, "raw_data_s4")
  expect_identical(workflow_data$design_matrix$tech_rep_group, c("A_1", "B_1"))

  workflow_data$design_matrix <- data.frame(Run = c("S1", "S2"), group = c("A", "B"), replicates = c(1, 1))
  workflow_data$data_cln <- data.frame(
    Protein.Ids = c("P1", "P1"),
    Run = c("S1", "S2"),
    Abundance = c(10, 30)
  )
  workflow_data$column_mapping <- list(protein_col = "Protein.Ids", run_col = "Run", quantity_col = "Abundance")
  tmt_state <- build_checkpoint(
    workflowData = workflow_data,
    workflowType = "TMT",
    actionLabel = "Import",
    validateColumnMapping = TRUE
  )
  expect_identical(tmt_state, "protein_s4_initial")
  expect_identical(saved_states[[2L]]$state_name, "protein_s4_initial")

  qc_value <- FALSE
  qc_trigger <- function(value) {
    if (missing(value)) {
      qc_value
    } else {
      qc_value <<- value
    }
  }
  workflow_data$tab_status <- list(design_matrix = "pending")
  workflow_data$taxon_id <- 9606
  workflow_data$column_mapping <- list(protein_col = "Protein.Ids")
  workflow_data$data_cln <- data.frame(Protein.Ids = "P1", Run = "S1", Abundance = 10)
  source_dir <- tempfile("prot-design-post-source-")
  results_dir <- tempfile("prot-design-post-results-")
  dir.create(source_dir, recursive = TRUE)
  dir.create(results_dir, recursive = TRUE)

  local_mocked_bindings(
    getUniprotAnnotationsFull = function(...) data.frame(Protein.Ids = "P1", Gene = "GENE1"),
    .env = asNamespace("MultiScholaR")
  )
  session <- shiny::MockShinySession$new()
  shiny::withReactiveDomain(session, {
    complete_checkpoint(
      workflowData = workflow_data,
      experimentPaths = list(source_dir = source_dir, results_dir = results_dir),
      session = session,
      qcTrigger = qc_trigger,
      successMessage = "Design complete",
      successNotificationId = "working"
    )
  })
  expect_true(qc_value)
  expect_identical(workflow_data$tab_status$design_matrix, "complete")
  expect_true(file.exists(file.path(source_dir, "uniprot_dat_cln.RDS")))
})

test_that("prot design builder host helpers preserve persistence, preview, and observer wiring", {
  skipIfMissingProtDesignBindings(
    "persistProtDesignBuilderArtifacts",
    "hydrateProtDesignBuilderResults",
    "runProtDesignBuilderSaveFlow",
    "runProtDesignBuilderObserverShell",
    "registerProtDesignPreviewOutputs",
    "registerProtDesignBuilderModule",
    "registerProtDesignBuilderResultsObserver"
  )

  persist_artifacts <- getProtDesignBinding("persistProtDesignBuilderArtifacts")
  hydrate_results <- getProtDesignBinding("hydrateProtDesignBuilderResults")
  run_save_flow <- getProtDesignBinding("runProtDesignBuilderSaveFlow")
  run_observer_shell <- getProtDesignBinding("runProtDesignBuilderObserverShell")
  register_preview <- getProtDesignBinding("registerProtDesignPreviewOutputs")
  register_builder_module <- getProtDesignBinding("registerProtDesignBuilderModule")
  register_results_observer <- getProtDesignBinding("registerProtDesignBuilderResultsObserver")

  source_dir <- tempfile("prot-design-builder-source-")
  dir.create(source_dir, recursive = TRUE)
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$config_list <- list(globalParameters = list(workflow_type = "DIA"))
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$workflow_type <- "DIA"

  results <- list(
    design_matrix = data.frame(Run = "S1", group = "A", stringsAsFactors = FALSE),
    data_cln = data.frame(Run = "S1", Protein.Ids = "P1", stringsAsFactors = FALSE),
    contrasts_tbl = data.frame(contrasts = "groupA-groupB", stringsAsFactors = FALSE),
    config_list = list(globalParameters = list(workflow_type = "DIA"))
  )

  persist_artifacts(results, workflow_data, source_dir)
  expect_true(file.exists(file.path(source_dir, "design_matrix.tab")))
  expect_true(file.exists(file.path(source_dir, "data_cln.tab")))
  expect_true(file.exists(file.path(source_dir, "contrast_strings.tab")))
  expect_true(file.exists(file.path(source_dir, "manifest.json")))

  hydrate_results(results, workflow_data)
  expect_identical(workflow_data$design_matrix$Run, "S1")
  expect_identical(workflow_data$config_list$globalParameters$workflow_type, "DIA")

  flow_calls <- character()
  local_mocked_bindings(
    persistProtDesignBuilderArtifacts = function(...) flow_calls <<- c(flow_calls, "persist"),
    buildProtDesignStateCheckpoint = function(...) {
      flow_calls <<- c(flow_calls, "checkpoint")
      "raw_data_s4"
    },
    completeProtDesignPostCheckpoint = function(...) flow_calls <<- c(flow_calls, "complete"),
    .env = asNamespace("MultiScholaR")
  )
  session <- shiny::MockShinySession$new()
  shiny::withReactiveDomain(session, {
    expect_true(run_save_flow(
      results = results,
      workflowData = workflow_data,
      experimentPaths = list(source_dir = source_dir),
      session = session,
      qcTrigger = function(...) TRUE
    ))
    expect_false(run_save_flow(
      results = results,
      workflowData = workflow_data,
      experimentPaths = list(source_dir = tempfile("missing-source-")),
      session = session
    ))
  })
  expect_identical(flow_calls, c("persist", "checkpoint", "complete"))

  observer_calls <- character()
  run_observer_shell(
    results = results,
    workflowData = workflow_data,
    experimentPaths = list(source_dir = source_dir),
    session = session,
    qcTrigger = function(...) TRUE,
    showProcessingModal = function() observer_calls <<- c(observer_calls, "modal"),
    hydrateBuilderResults = function(...) observer_calls <<- c(observer_calls, "hydrate"),
    runBuilderSaveFlow = function(...) observer_calls <<- c(observer_calls, "save"),
    logInfo = function(...) invisible(NULL)
  )
  expect_identical(observer_calls, c("modal", "hydrate", "save"))

  output <- new.env(parent = emptyenv())
  workflow_data$data_tbl <- data.frame(Run = "S1")
  workflow_data$config_list <- list()
  workflow_data$design_matrix <- data.frame(Run = "S1")
  workflow_data$contrasts_tbl <- data.frame(contrast = "A-B")
  register_preview(
    output = output,
    workflowData = workflow_data,
    reactiveFn = function(expr) eval.parent(substitute(expr)),
    outputOptionsFn = function(...) invisible(NULL),
    renderDT = function(expr, ...) eval.parent(substitute(expr)),
    reqFn = force
  )
  expect_true(output$data_available)
  expect_true(output$design_matrix_exists)
  expect_identical(output$design_matrix_preview$Run, "S1")
  expect_identical(output$contrasts_preview$contrast, "A-B")

  builder_args <- NULL
  builder_result <- register_builder_module(
    workflowData = workflow_data,
    moduleId = "builder",
    builderServerExists = TRUE,
    builderServerFn = function(id, data_tbl, config_list, column_mapping) {
      builder_args <<- list(id = id, data_tbl = data_tbl(), config_list = config_list(), column_mapping = column_mapping())
      "builder-rv"
    },
    reactiveFn = function(value) function() value
  )
  expect_identical(builder_result, "builder-rv")
  expect_identical(builder_args$id, "builder")

  fallback <- register_builder_module(
    workflowData = workflow_data,
    builderServerExists = FALSE,
    reactiveValFn = function(value) structure(list(value = value), class = "mock_reactive_val")
  )
  expect_s3_class(fallback, "mock_reactive_val")

  observed <- NULL
  register_results_observer(
    builderResultsRv = function() results,
    workflowData = workflow_data,
    experimentPaths = list(source_dir = source_dir),
    session = session,
    observeEventFn = function(eventExpr, handlerExpr, ..., ignoreNULL = TRUE) {
      force(eventExpr)
      force(handlerExpr)
    },
    reqFn = force,
    runBuilderObserverShell = function(results, ...) observed <<- results
  )
  expect_identical(observed$design_matrix$Run, "S1")
})

test_that("prot design import workflow helpers preserve flow, observers, and bootstrap wiring", {
  skipIfMissingProtDesignBindings(
    "initializeProtDesignImportedWorkflowState",
    "runProtDesignImportConfirmationFlow",
    "runProtDesignImportObserverShell",
    "registerProtDesignImportConfirmationObserver",
    "initializeProtDesignImportBootstrap",
    "registerProtDesignImportModalShell"
  )

  initialize_state <- getProtDesignBinding("initializeProtDesignImportedWorkflowState")
  run_flow <- getProtDesignBinding("runProtDesignImportConfirmationFlow")
  run_observer_shell <- getProtDesignBinding("runProtDesignImportObserverShell")
  register_confirmation <- getProtDesignBinding("registerProtDesignImportConfirmationObserver")
  initialize_bootstrap <- getProtDesignBinding("initializeProtDesignImportBootstrap")
  register_modal <- getProtDesignBinding("registerProtDesignImportModalShell")

  workflow_data <- new.env(parent = emptyenv())
  workflow_type <- NULL
  workflow_data$state_manager <- list(setWorkflowType = function(value) workflow_type <<- value)
  workflow_data$config_list <- list()
  workflow_data$column_mapping <- NULL
  imported_design <- data.frame(Run = "S1", group = "A", replicates = 1)
  imported_data <- data.frame(
    Precursor.Id = "p1",
    Precursor.Charge = 2,
    Precursor.Quantity = 10,
    Stripped.Sequence = "AAA",
    Protein.Ids = "P1",
    Run = "S1"
  )
  imported_contrasts <- data.frame(contrasts = "groupA-groupB", stringsAsFactors = FALSE)

  detected_type <- initialize_state(
    workflowData = workflow_data,
    importedDesign = imported_design,
    importedDataCln = imported_data,
    importedContrasts = imported_contrasts,
    taxonId = 9606,
    organismName = "Homo sapiens"
  )
  expect_identical(detected_type, "DIA")
  expect_identical(workflow_type, "DIA")
  expect_identical(workflow_data$column_mapping$quantity_col, "Precursor.Quantity")

  lfq_type <- NULL
  lfq_workflow_data <- new.env(parent = emptyenv())
  lfq_workflow_data$state_manager <- list(setWorkflowType = function(value) lfq_type <<- value)
  lfq_workflow_data$config_list <- list(globalParameters = list(workflow_type = "LFQ"))
  lfq_workflow_data$column_mapping <- NULL
  expect_identical(
    initialize_state(
      workflowData = lfq_workflow_data,
      importedDesign = imported_design,
      importedDataCln = data.frame(Protein.Ids = "P1", Run = "S1", Intensity = 42),
      importedContrasts = NULL,
      taxonId = 9606,
      organismName = "Homo sapiens"
    ),
    "LFQ"
  )
  expect_identical(lfq_type, "LFQ")
  expect_identical(lfq_workflow_data$column_mapping$quantity_col, "Intensity")

  tmt_type <- NULL
  tmt_workflow_data <- new.env(parent = emptyenv())
  tmt_workflow_data$state_manager <- list(setWorkflowType = function(value) tmt_type <<- value)
  tmt_workflow_data$config_list <- list(globalParameters = list(workflow_type = "TMT"))
  tmt_workflow_data$column_mapping <- NULL
  expect_identical(
    initialize_state(
      workflowData = tmt_workflow_data,
      importedDesign = imported_design,
      importedDataCln = data.frame(Protein.Ids = "P1", Run = "S1", Abundance = 42),
      importedContrasts = NULL,
      taxonId = 9606,
      organismName = "Homo sapiens"
    ),
    "TMT"
  )
  expect_identical(tmt_type, "TMT")
  expect_identical(tmt_workflow_data$column_mapping$quantity_col, "Abundance")

  flow_calls <- character()
  flow_result <- run_flow(
    workflowData = workflow_data,
    experimentPaths = list(source_dir = tempdir(), results_dir = tempdir()),
    importArtifacts = list(importPath = tempdir(), fastaPath = "selected.fasta"),
    importedArtifacts = list(
      importedDesign = imported_design,
      importedDataCln = imported_data,
      importedContrasts = imported_contrasts
    ),
    taxonId = 9606,
    organismName = "Homo sapiens",
    session = shiny::MockShinySession$new(),
    qcTrigger = function(...) TRUE,
    hydrateFastaSidecar = function(...) flow_calls <<- c(flow_calls, "fasta"),
    hydrateUniprotSidecar = function(...) flow_calls <<- c(flow_calls, "uniprot"),
    initializeWorkflowState = function(...) {
      flow_calls <<- c(flow_calls, "state")
      "DIA"
    },
    buildStateCheckpointFn = function(...) {
      flow_calls <<- c(flow_calls, "checkpoint")
      "raw_data_s4"
    },
    completePostCheckpointFn = function(...) flow_calls <<- c(flow_calls, "complete")
  )
  expect_identical(flow_result$workflowType, "DIA")
  expect_identical(flow_result$stateName, "raw_data_s4")
  expect_identical(flow_calls, c("fasta", "uniprot", "state", "checkpoint", "complete"))

  notifications <- list()
  fail_result <- run_observer_shell(
    workflowData = workflow_data,
    experimentPaths = list(source_dir = tempdir()),
    importPath = tempdir(),
    taxonId = 9606,
    organismName = "Homo sapiens",
    session = shiny::MockShinySession$new(),
    resolveImportArtifacts = function(...) list(ok = FALSE, errorMessage = "missing import"),
    showNotification = function(message, type = NULL, duration = NULL, id = NULL) {
      notifications[[length(notifications) + 1L]] <<- list(message = message, type = type, duration = duration, id = id)
    },
    removeNotification = function(...) invisible(NULL)
  )
  expect_false(fail_result$ok)
  expect_identical(notifications[[1L]]$message, "missing import")

  success_result <- run_observer_shell(
    workflowData = workflow_data,
    experimentPaths = list(source_dir = tempdir()),
    importPath = tempdir(),
    taxonId = 9606,
    organismName = "Homo sapiens",
    session = shiny::MockShinySession$new(),
    resolveImportArtifacts = function(...) list(ok = TRUE, importPath = "import", designFile = "design", dataClnFile = "data", contrastFile = "contrast"),
    loadImportedArtifacts = function(...) list(importedDesign = imported_design, importedDataCln = imported_data, importedContrasts = imported_contrasts),
    runImportConfirmationFlow = function(...) list(workflowType = "DIA", stateName = "raw_data_s4"),
    showNotification = function(...) invisible(NULL),
    removeNotification = function(...) invisible(NULL)
  )
  expect_true(success_result$ok)
  expect_identical(success_result$flowResult$stateName, "raw_data_s4")

  input <- list(
    confirm_import = 1,
    import_dir = "encoded-dir",
    import_taxon_id = 9606,
    import_organism_name = "Homo sapiens"
  )
  observed_import <- NULL
  register_confirmation(
    input = input,
    resolvedVolumes = c(Home = tempdir()),
    importFastaPath = function() "selected.fasta",
    workflowData = workflow_data,
    experimentPaths = list(source_dir = tempdir()),
    session = shiny::MockShinySession$new(),
    observeEventFn = function(eventExpr, handlerExpr, ...) {
      force(eventExpr)
      force(handlerExpr)
    },
    reqFn = force,
    parseDirPathFn = function(roots, selection) "/tmp/import",
    removeModalFn = function() invisible(NULL),
    runImportObserverShell = function(...) observed_import <<- list(...)
  )
  expect_identical(observed_import$importPath, "/tmp/import")
  expect_identical(observed_import$selectedFastaPath, "selected.fasta")

  bootstrap_calls <- list()
  bootstrap <- initialize_bootstrap(
    input = list(),
    session = shiny::MockShinySession$new(),
    experimentPaths = list(base_dir = tempdir()),
    volumes = function() c(Home = tempdir()),
    dirChooseFn = function(input, name, roots, session) bootstrap_calls$dir <<- list(name = name, roots = roots),
    fileChooseFn = function(input, name, roots, session, filetypes) bootstrap_calls$file <<- list(name = name, roots = roots, filetypes = filetypes),
    reactiveValFn = function(value) {
      stored <- value
      function(newValue) {
        if (!missing(newValue)) stored <<- newValue
        stored
      }
    },
    isolateFn = force,
    dirExistsFn = function(path) TRUE,
    logInfo = function(...) invisible(NULL)
  )
  expect_identical(bootstrap_calls$dir$name, "import_dir")
  expect_identical(bootstrap_calls$file$filetypes, c("fasta", "fa", "faa"))
  expect_true("Project Base Dir" %in% names(bootstrap$resolvedVolumes))
  expect_null(bootstrap$importFastaPath())

  session <- shiny::MockShinySession$new()
  selected_fasta <- "preexisting"
  import_fasta_path <- function(value) {
    if (!missing(value)) selected_fasta <<- value
    selected_fasta
  }
  shiny::withReactiveDomain(session, {
    register_modal(
      input = session$input,
      output = session$output,
      session = session,
      resolvedVolumes = c(Home = tempdir()),
      importFastaPath = import_fasta_path
    )
    session$setInputs(show_import_modal = 1)
    shiny:::flushReact()
  })
  expect_null(import_fasta_path())
})
