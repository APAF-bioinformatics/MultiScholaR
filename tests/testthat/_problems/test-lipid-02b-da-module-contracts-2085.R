# Extracted from test-lipid-02b-da-module-contracts.R:2085

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
workflow_data <- new.env(parent = emptyenv())
workflow_data$contrasts_tbl <- data.frame(
        contrasts = c("B_vs_A", "C_vs_A"),
        stringsAsFactors = FALSE
    )
resolved_s4 <- structure(list(marker = TRUE), class = "LipidomicsAssayData")
notifications <- list()
info_messages <- character()
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
        logInfo = function(message) {
            info_messages <<- c(info_messages, message)
        },
        resolveCurrentStateFn = function(workflowData) {
            resolve_calls <<- resolve_calls + 1L
            expect_identical(
                workflowData$contrasts_tbl$contrasts,
                c("B_vs_A", "C_vs_A")
            )
            resolved_s4
        }
    )
expect_identical(resolve_calls, 1L)
expect_identical(result$currentS4, resolved_s4)
