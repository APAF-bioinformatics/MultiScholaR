# Extracted from test-lipid-01b-qc-module-contracts.R:108

# prequel ----------------------------------------------------------------------
library(testthat)
suppressPackageStartupMessages(library(shiny))
loadLipidQcModuleHarness <- function(
    stateFactory = function() NULL,
    qcTriggerValue = FALSE
) {
    source_lines <- readLines(
        test_path("..", "..", "R", "mod_lipid_qc.R"),
        warn = FALSE
    )
    source_lines <- sub(
        "shiny::moduleServer",
        "moduleServer",
        source_lines,
        fixed = TRUE
    )

    module_env <- new.env(parent = globalenv())
    module_env$moduleServer <- function(id, module, session = shiny::getDefaultReactiveDomain()) {
        assign("capturedModule", module, envir = module_env)
        invisible(NULL)
    }

    eval(parse(text = source_lines), envir = module_env)

    module_env$mod_lipid_qc_intensity_ui <- function(id) shiny::tabPanel("Intensity")
    module_env$mod_lipid_qc_duplicates_ui <- function(id) shiny::tabPanel("Duplicates")
    module_env$mod_lipid_qc_itsd_ui <- function(id) shiny::tabPanel("ITSD")
    module_env$mod_lipid_qc_s4_ui <- function(id) shiny::tabPanel("Finalize")

    module_env$serverCalls <- character()
    module_env$mod_lipid_qc_intensity_server <- function(...) {
        module_env$serverCalls <- c(module_env$serverCalls, "intensity")
    }
    module_env$mod_lipid_qc_duplicates_server <- function(...) {
        module_env$serverCalls <- c(module_env$serverCalls, "duplicates")
    }
    module_env$mod_lipid_qc_itsd_server <- function(...) {
        module_env$serverCalls <- c(module_env$serverCalls, "itsd")
    }
    module_env$mod_lipid_qc_s4_server <- function(...) {
        module_env$serverCalls <- c(module_env$serverCalls, "s4")
    }

    qc_trigger <- shiny::reactiveVal(qcTriggerValue)
    workflow_data <- list(
        state_manager = list(getState = stateFactory)
    )

    module_env$mod_lipid_qc_server(
        "qc",
        workflow_data = workflow_data,
        experiment_paths = NULL,
        omic_type = "lipidomics",
        experiment_label = "Lipidomics",
        qc_trigger = qc_trigger
    )

    session <- shiny:::MockShinySession$new()
    shiny::callModule(module_env$capturedModule, "qc", session = session)
    session$flushReact()

    list(
        moduleEnv = module_env,
        session = session
    )
}
renderLipidQcTabs <- function(harness) {
    htmltools::renderTags(
        harness$session$getOutput("qc-dynamic_qc_tabs")
    )$html
}

# test -------------------------------------------------------------------------
harness <- loadLipidQcModuleHarness(
        stateFactory = function() structure(list(), class = "LipidomicsAssayData"),
        qcTriggerValue = FALSE
    )
rendered <- renderLipidQcTabs(harness)
expect_match(rendered, "qc-lipid_qc_tabs", fixed = TRUE)
