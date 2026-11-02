# Extracted from test-module-ci-lipid-da.R:192

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
object <- module_ci_lipid_da_object(layout = "pos_gcms")
valid_cases <- list(
    two_group = list(
      contrasts = module_ci_lipid_da_contrasts("two_group"),
      formula = "~ 0 + group"
    ),
    multi_group = list(
      contrasts = module_ci_lipid_da_contrasts("multi_group"),
      formula = "~ 0 + group"
    ),
    batch_aware = list(
      contrasts = module_ci_lipid_da_contrasts("batch_aware"),
      formula = "~ 0 + group + batch"
    ),
    reversed = list(
      contrasts = module_ci_lipid_da_contrasts("reversed"),
      formula = "~ 0 + group"
    )
  )
for (case in valid_cases) {
    preflight <- validateLipidDaAnalysisPreflight(
      currentS4 = object,
      contrastsTbl = case$contrasts,
      formulaString = case$formula
    )
    expect_true(preflight$valid)
    expect_true(all(case$contrasts$contrasts %in% preflight$contrast_validation$contrasts_tbl$contrasts))
  }
duplicate <- validateLipidDaAnalysisPreflight(
    currentS4 = object,
    contrastsTbl = module_ci_lipid_da_contrasts("duplicate"),
    formulaString = "~ 0 + group"
  )
expect_false(duplicate$valid)
expect_true(any(grepl("duplicate friendly contrast labels", duplicate$errors, fixed = TRUE)))
invalid_term <- validateLipidDaAnalysisPreflight(
    currentS4 = object,
    contrastsTbl = module_ci_lipid_da_contrasts("invalid_term"),
    formulaString = "~ 0 + group"
  )
expect_false(invalid_term$valid)
expect_true(any(grepl("references terms absent", invalid_term$errors, fixed = TRUE)))
invalid_formula <- validateLipidDaAnalysisPreflight(
    currentS4 = object,
    contrastsTbl = module_ci_lipid_da_contrasts("two_group"),
    formulaString = "~ 0 + missing_factor"
  )
expect_false(invalid_formula$valid)
expect_true(any(grepl("cannot produce a model matrix", invalid_formula$errors, fixed = TRUE)))
no_contrasts_context <- prepareLipidDaRunAnalysisContext(
    daData = list2env(list(current_state_name = "lipid_correlation_filtered"), parent = emptyenv()),
    workflowData = module_ci_lipid_da_workflow(object, contrasts = module_ci_lipid_da_contrasts("no_contrasts")),
    notify = function(...) invisible(NULL),
    logInfo = function(...) invisible(NULL)
  )
expect_null(no_contrasts_context)
invalid_workflow <- module_ci_lipid_da_workflow(
    object,
    contrasts = module_ci_lipid_da_contrasts("invalid_term")
  )
capture <- module_ci_lipid_da_capture_ui()
shell_called <- FALSE
analysis_state <- handleLipidDaRunAnalysisPreflight(
    analysisContext = list(
      currentS4 = object,
      contrastsTbl = invalid_workflow$contrasts_tbl,
      runningNotification = list(id = "da_running")
    ),
    formulaString = "~ 0 + group",
    daQValThresh = 0.05,
    treatLfcCutoff = 0,
    daData = module_ci_lipid_da_data_env(),
    workflowData = invalid_workflow,
    session = "module-ci-lipid-da",
    experimentPaths = list(),
    removeNotificationFn = capture$removeNotification,
    notify = capture$showNotification,
    executeRunAnalysisFn = function(...) {
      shell_called <<- TRUE
      stop("analysis shell must not run for invalid preflight")
    }
  )
