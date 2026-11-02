# Extracted from test-lipid-02b-da-module-contracts.R:102

# prequel ----------------------------------------------------------------------
library(testthat)
source(test_path("..", "..", "R", "mod_lipid_da_display_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_da_results_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_da_session_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_da_analysis_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_da_export_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_da_ui.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_da_server.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_da.R"), local = environment())
if (!methods::isClass("MockLipidDaArgsCarrier")) {
    methods::setClass("MockLipidDaArgsCarrier", slots = c(args = "list"))
}
newMockLipidDaArgsCarrier <- function(args = list()) {
    methods::new("MockLipidDaArgsCarrier", args = args)
}

# test -------------------------------------------------------------------------
omicsStandInReaders()
output <- new.env(parent = emptyenv())
workflow_data <- new.env(parent = emptyenv())
da_data <- new.env(parent = emptyenv())
render_calls <- character()
captured <- list()
workflow_data$contrasts_tbl <- data.frame(
        friendly_names = c("B vs A", "C vs A"),
        stringsAsFactors = FALSE
    )
da_data$analysis_complete <- TRUE
da_data$da_results_list <- list(marker = "results")
mock_s4 <- structure(list(marker = "s4"), class = "mock_s4")
da_data$current_state_name <- "loaded_for_de"
workflow_data$state_manager <- list(getState = function(state_name = NULL) {
        if (identical(state_name, "loaded_for_de")) mock_s4 else NULL
    })
result <- registerLipidDaPrimaryTextOutputs(
        output = output,
        workflowData = workflow_data,
        daData = da_data,
        renderPrintFn = function(expr) {
            render_calls <<- c(render_calls, force(expr))
            structure(list(kind = "mock_render", index = length(render_calls)), class = "mock_render")
        },
        buildContrastsTextFn = function(contrastsTbl) {
            captured$contrastsTbl <<- contrastsTbl
            "contrasts-text"
        },
        buildStatusTextFn = function(analysisComplete, daResultsList, currentS4Object) {
            captured$statusArgs <<- list(
                analysisComplete = analysisComplete,
                daResultsList = daResultsList,
                currentS4Object = currentS4Object
            )
            "status-text"
        },
        catFn = function(...) paste0(..., collapse = "")
    )
expect_identical(result, output)
expect_identical(render_calls, c("contrasts-text", "status-text"))
expect_identical(captured$contrastsTbl, workflow_data$contrasts_tbl)
expect_true(captured$statusArgs$analysisComplete)
expect_identical(captured$statusArgs$daResultsList, da_data$da_results_list)
expect_identical(captured$statusArgs$currentS4Object, mock_s4)
