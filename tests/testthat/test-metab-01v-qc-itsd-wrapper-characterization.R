# fidelity-coverage-compare: shared
library(testthat)

test_that("metabolomics ITSD QC server preserves module registration handoff", {
  namespace <- asNamespace("MultiScholaR")
  skip_if_not(
    exists("mod_metab_qc_itsd_server", envir = namespace, inherits = FALSE),
    "mod_metab_qc_itsd_server is unavailable in this ref."
  )

  server <- get("mod_metab_qc_itsd_server", envir = namespace, inherits = FALSE)
  server_env <- environment(server)
  has_body_seam <- exists("runMetabQcItsdServerBody", envir = server_env, inherits = FALSE)
  workflow_data <- list(state_manager = list(getState = function() "current-s4"))
  captured <- new.env(parent = emptyenv())

  if (has_body_seam) {
    testthat::local_mocked_bindings(
      runMetabQcItsdServerBody = function(input, output, session, workflowData, omicType, experimentLabel, ...) {
        captured$body <- list(
          input = input,
          output = output,
          session = session,
          workflowData = workflowData,
          omicType = omicType,
          experimentLabel = experimentLabel
        )
        invisible(NULL)
      },
      .env = server_env
    )
  }

  testthat::local_mocked_bindings(
    moduleServer = function(id, module, ...) {
      captured$id <- id
      captured$module <- module

      if (isTRUE(has_body_seam)) {
        output <- new.env(parent = emptyenv())
        module(
          input = list(analyze_is = FALSE, is_pattern = "^IS_"),
          output = output,
          session = list(ns = function(value) paste(id, value, sep = "-"))
        )
      }

      invisible("registered")
    },
    .package = "shiny"
  )

  result <- server(
    id = "itsd",
    workflow_data = workflow_data,
    omic_type = "metabolomics",
    experiment_label = "Experiment A"
  )

  expect_identical(result, "registered")
  expect_identical(captured$id, "itsd")
  expect_true(is.function(captured$module))

  if (has_body_seam) {
    expect_identical(captured$body$workflowData, workflow_data)
    expect_identical(captured$body$omicType, "metabolomics")
    expect_identical(captured$body$experimentLabel, "Experiment A")
    expect_identical(captured$body$session$ns("cv_plot"), "itsd-cv_plot")
  }
})
