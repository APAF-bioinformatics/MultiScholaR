# Extracted from test-lipid-02b-da-module-contracts.R:2032

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
da_data <- new.env(parent = emptyenv())
da_data$current_state_name <- "loaded_for_de"
workflow_data <- new.env(parent = emptyenv())
workflow_data$state_manager <- list(getState = function(state_name = NULL) {
        if (identical(state_name, "loaded_for_de")) structure(list(marker = TRUE), class = "LipidomicsAssayData") else NULL
    })
workflow_data$contrasts_tbl <- data.frame()
notifications <- list()
resolve_calls <- 0L
result <- prepareLipidDaRunAnalysisContext(
        daData = da_data,
        workflowData = workflow_data,
        notify = function(message, type = NULL, duration = NULL, id = NULL) {
            notifications[[length(notifications) + 1]] <<- list(
                message = message,
                type = type,
                duration = duration,
                id = id
            )
        },
        logInfo = function(...) NULL,
        resolveCurrentStateFn = function(...) {
            resolve_calls <<- resolve_calls + 1L
            stop("resolver should not run when da_data already has an S4 object")
        }
    )
