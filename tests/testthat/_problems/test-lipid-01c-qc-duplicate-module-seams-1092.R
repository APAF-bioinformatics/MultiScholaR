# Extracted from test-lipid-01c-qc-duplicate-module-seams.R:1092

# prequel ----------------------------------------------------------------------
library(testthat)
source(test_path("..", "..", "R", "mod_lipid_qc_duplicates_server_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_qc_duplicates.R"), local = environment())
ensureLipidDuplicateWorkflowTestClass <- function() {
    if (!methods::isClass("LipidomicsAssayData")) {
        methods::setClass(
            "LipidomicsAssayData",
            slots = c(
                lipid_data = "list",
                lipid_id_column = "character"
            )
        )
    }
}

# test -------------------------------------------------------------------------
reverted_state <- NULL
info_messages <- character()
result <- handleLipidDuplicateRevert(
        workflowData = list(
            state_manager = list(
                getHistory = function() c("initial", "lipid_duplicates_resolved", "lipid_filtered"),
                revertToState = function(state_name) {
                    reverted_state <<- state_name
                    list(reverted = state_name)
                }
            )
        ),
        reqFn = function(...) invisible(NULL),
        logInfoFn = function(msg) {
            info_messages <<- c(info_messages, msg)
        }
    )
