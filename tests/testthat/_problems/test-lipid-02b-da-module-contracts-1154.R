# Extracted from test-lipid-02b-da-module-contracts.R:1154

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
update_calls <- list()
debug_messages <- character()
warn_message <- NULL
result <- restoreLipidDaFormulaFromSession(
        sessionData = list(current_s4_object = NULL),
        daData = da_data,
        session = structure(list(kind = "mock-session"), class = "mock_session"),
        updateTextAreaInputFn = function(...) {
            update_calls[[length(update_calls) + 1]] <<- list(...)
        },
        debugLog = function(...) {
            debug_messages <<- c(debug_messages, paste0(..., collapse = ""))
        },
        hasArgsMethodFn = function() FALSE,
        logWarn = function(message) {
            warn_message <<- message
        }
    )
expect_null(result)
expect_length(update_calls, 0)
expect_false(exists("formula_from_s4", envir = da_data, inherits = FALSE))
expect_identical(debug_messages[1:2], c(
        "  STEP 9a: Checking S4 @args slot...",
        "    hasSlot 'args': FALSE"
    ))
expect_match(debug_messages[3], "^  STEP 9 ERROR: ")
expect_match(warn_message, "^Could not extract formula from S4: ")
