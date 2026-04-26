# fidelity-coverage-compare: shared
library(testthat)

msr <- function(name) {
  get(name, envir = asNamespace("MultiScholaR"), inherits = FALSE)
}

mod_prot_qc_ui <- msr("mod_prot_qc_ui")
mod_prot_qc_server <- msr("mod_prot_qc_server")

makeQcStateManager <- function(workflow_type, s4_obj = NULL, final_state = TRUE) {
  structure(
    list(
      workflow_type = workflow_type,
      states = list(protein_replicate_filtered = final_state),
      getState = function(name) {
        if (identical(name, "protein_s4_initial")) {
          return(s4_obj)
        }
        NULL
      }
    ),
    class = "fake_qc_state_manager"
  )
}

newQcWorkflowData <- function(workflow_type, s4_obj = NULL, final_state = TRUE) {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- makeQcStateManager(
    workflow_type = workflow_type,
    s4_obj = s4_obj,
    final_state = final_state
  )
  workflow_data
}

localQcServerMocks <- function(captured, .local_envir = parent.frame()) {
  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      output <- new.env(parent = emptyenv())
      session <- list(ns = function(value) paste(id, value, sep = "-"))
      module(list(), output, session)
      captured$output <- output
      invisible(NULL)
    },
    renderUI = function(expr) {
      eval(substitute(expr), parent.frame())
    },
    observeEvent = function(eventExpr, handlerExpr, ..., ignoreNULL = FALSE, once = FALSE) {
      value <- eval(substitute(eventExpr), parent.frame())
      trigger <- TRUE

      if (is.null(value) || identical(value, FALSE)) {
        trigger <- FALSE
      }

      if (trigger) {
        eval(substitute(handlerExpr), parent.frame())
      }
      invisible(NULL)
    },
    req = function(...) {
      values <- list(...)
      for (value in values) {
        if (is.null(value) || identical(value, FALSE)) {
          stop("required value missing", call. = FALSE)
        }
      }
      invisible(values[[1]])
    },
    isolate = function(expr) eval(substitute(expr), parent.frame()),
    tabsetPanel = function(..., id = NULL) {
      list(kind = "tabsetPanel", id = id, tabs = list(...))
    },
    tabPanel = function(title, ...) {
      list(title = title, contents = list(...))
    },
    NS = function(id) function(value) paste(id, value, sep = "-"),
    tagList = function(...) list(...),
    fluidRow = function(...) list(...),
    column = function(width, ...) list(width = width, children = list(...)),
    wellPanel = function(...) list(...),
    h3 = function(...) paste(...),
    uiOutput = function(id) list(kind = "uiOutput", id = id),
    .package = "shiny",
    .env = .local_envir
  )
}

test_that("proteomics QC coordinator preserves the TMT existing-S4 routing path", {
  captured <- new.env(parent = emptyenv())
  captured$logs <- character()
  captured$calls <- list()

  localQcServerMocks(captured)
  testthat::local_mocked_bindings(
    mod_prot_qc_protein_ui = function(id, workflow_type) paste("protein-ui", id, workflow_type, sep = "::"),
    mod_prot_qc_protein_s4_ui = function(id) paste("protein-s4-ui", id, sep = "::"),
    mod_prot_qc_peptide_ui = function(id) paste("peptide-ui", id, sep = "::"),
    mod_prot_qc_protein_server = function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
      captured$calls[[length(captured$calls) + 1L]] <<- list(
        fn = "protein_server",
        id = id,
        workflow_data = workflow_data,
        experiment_paths = experiment_paths,
        omic_type = omic_type,
        experiment_label = experiment_label
      )
      invisible(NULL)
    },
    mod_prot_qc_protein_s4_server = function(...) {
      captured$calls[[length(captured$calls) + 1L]] <<- list(fn = "protein_s4_server")
      invisible(NULL)
    },
    mod_prot_qc_peptide_server = function(...) {
      captured$calls[[length(captured$calls) + 1L]] <<- list(fn = "peptide_server")
      invisible(NULL)
    },
    .package = "MultiScholaR"
  )
  testthat::local_mocked_bindings(
    log_info = function(message, ...) {
      captured$logs <- c(captured$logs, message)
      invisible(NULL)
    },
    log_warn = function(message, ...) {
      captured$logs <- c(captured$logs, message)
      invisible(NULL)
    },
    .package = "logger"
  )

  workflow_data <- newQcWorkflowData(workflow_type = "TMT", s4_obj = structure(list(), class = "fake_s4"))
  experiment_paths <- list(protein_qc_dir = tempdir())

  mod_prot_qc_server(
    id = "qc",
    workflow_data = workflow_data,
    experiment_paths = experiment_paths,
    omic_type = "proteomics",
    experiment_label = "TMT Experiment",
    qc_trigger = function() TRUE
  )

  tabs <- captured$output$dynamic_qc_tabs
  expect_identical(tabs$id, "qc-qc_tabs_tmt")
  expect_length(tabs$tabs, 1L)
  expect_identical(tabs$tabs[[1L]]$title, "Protein QC")
  expect_identical(tabs$tabs[[1L]]$contents[[1L]], "protein-ui::qc-protein_qc::TMT")
  expect_identical(vapply(captured$calls, `[[`, character(1), "fn"), c("protein_server"))
  expect_true(any(grepl("Skipping S4 creation UI", captured$logs, fixed = TRUE)))
  expect_true(any(grepl("Final protein state reached", captured$logs, fixed = TRUE)))
})

test_that("proteomics QC coordinator preserves the TMT missing-S4 fallback path", {
  captured <- new.env(parent = emptyenv())
  captured$logs <- character()
  captured$calls <- list()

  localQcServerMocks(captured)
  testthat::local_mocked_bindings(
    mod_prot_qc_protein_ui = function(id, workflow_type) paste("protein-ui", id, workflow_type, sep = "::"),
    mod_prot_qc_protein_s4_ui = function(id) paste("protein-s4-ui", id, sep = "::"),
    mod_prot_qc_peptide_ui = function(id) paste("peptide-ui", id, sep = "::"),
    mod_prot_qc_protein_server = function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
      captured$calls[[length(captured$calls) + 1L]] <<- list(fn = "protein_server", id = id)
      invisible(NULL)
    },
    mod_prot_qc_protein_s4_server = function(id, workflow_data, omic_type, experiment_label) {
      captured$calls[[length(captured$calls) + 1L]] <<- list(fn = "protein_s4_server", id = id)
      invisible(NULL)
    },
    mod_prot_qc_peptide_server = function(...) {
      captured$calls[[length(captured$calls) + 1L]] <<- list(fn = "peptide_server")
      invisible(NULL)
    },
    .package = "MultiScholaR"
  )
  testthat::local_mocked_bindings(
    log_info = function(message, ...) {
      captured$logs <- c(captured$logs, message)
      invisible(NULL)
    },
    log_warn = function(message, ...) {
      captured$logs <- c(captured$logs, message)
      invisible(NULL)
    },
    .package = "logger"
  )

  workflow_data <- newQcWorkflowData(workflow_type = "LFQ", s4_obj = NULL)

  mod_prot_qc_server(
    id = "qc",
    workflow_data = workflow_data,
    experiment_paths = list(protein_qc_dir = tempdir()),
    omic_type = "proteomics",
    experiment_label = "LFQ Experiment",
    qc_trigger = function() TRUE
  )

  tabs <- captured$output$dynamic_qc_tabs
  expect_identical(tabs$id, "qc-qc_tabs_tmt")
  expect_identical(vapply(tabs$tabs, `[[`, character(1), "title"), c("Initial Protein Processing", "Protein QC"))
  expect_identical(
    vapply(captured$calls, `[[`, character(1), "fn"),
    c("protein_s4_server", "protein_server")
  )
  expect_true(any(grepl("S4 object not found", captured$logs, fixed = TRUE)))
})

test_that("proteomics QC coordinator preserves the DIA peptide-plus-protein routing path", {
  captured <- new.env(parent = emptyenv())
  captured$logs <- character()
  captured$calls <- list()

  localQcServerMocks(captured)
  testthat::local_mocked_bindings(
    mod_prot_qc_protein_ui = function(id, workflow_type) paste("protein-ui", id, workflow_type, sep = "::"),
    mod_prot_qc_protein_s4_ui = function(id) paste("protein-s4-ui", id, sep = "::"),
    mod_prot_qc_peptide_ui = function(id) paste("peptide-ui", id, sep = "::"),
    mod_prot_qc_protein_server = function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
      captured$calls[[length(captured$calls) + 1L]] <<- list(fn = "protein_server", id = id)
      invisible(NULL)
    },
    mod_prot_qc_protein_s4_server = function(...) {
      captured$calls[[length(captured$calls) + 1L]] <<- list(fn = "protein_s4_server")
      invisible(NULL)
    },
    mod_prot_qc_peptide_server = function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
      captured$calls[[length(captured$calls) + 1L]] <<- list(fn = "peptide_server", id = id)
      invisible(NULL)
    },
    .package = "MultiScholaR"
  )
  testthat::local_mocked_bindings(
    log_info = function(message, ...) {
      captured$logs <- c(captured$logs, message)
      invisible(NULL)
    },
    log_warn = function(message, ...) {
      captured$logs <- c(captured$logs, message)
      invisible(NULL)
    },
    .package = "logger"
  )

  workflow_data <- newQcWorkflowData(workflow_type = "DIA", s4_obj = NULL)

  mod_prot_qc_server(
    id = "qc",
    workflow_data = workflow_data,
    experiment_paths = list(protein_qc_dir = tempdir()),
    omic_type = "proteomics",
    experiment_label = "DIA Experiment",
    qc_trigger = function() TRUE
  )

  tabs <- captured$output$dynamic_qc_tabs
  expect_identical(tabs$id, "qc-qc_tabs_lfq")
  expect_identical(vapply(tabs$tabs, `[[`, character(1), "title"), c("Peptide QC", "Protein QC"))
  expect_identical(
    vapply(captured$calls, `[[`, character(1), "fn"),
    c("peptide_server", "protein_server")
  )
  expect_true(any(grepl("Routing for workflow type: DIA", captured$logs, fixed = TRUE)))
})
