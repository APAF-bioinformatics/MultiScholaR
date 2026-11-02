# Extracted from test-lipid-01c-qc-duplicate-module-seams.R:397

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
ensureLipidDuplicateWorkflowTestClass()
current_s4 <- methods::new(
        "LipidomicsAssayData",
        lipid_data = list(),
        lipid_id_column = "lipid_id",
        design_matrix = data.frame(Sample_ID = "sample_a", group = "A"),
        sample_id = "Sample_ID",
        group_id = "group"
    )
result <- handleLipidDuplicateDetection(
        workflowData = list(
            state_manager = list(
                getState = function() current_s4
            )
        ),
        reqFn = function(...) invisible(NULL),
        findDuplicatesFn = function(...) {
            list(
                AssayA = NULL,
                AssayB = data.frame(lipid_id = character(), count = integer())
            )
        },
        logInfoFn = function(...) invisible(NULL)
    )
