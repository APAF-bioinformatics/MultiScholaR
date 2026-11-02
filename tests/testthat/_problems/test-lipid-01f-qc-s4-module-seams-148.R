# Extracted from test-lipid-01f-qc-s4-module-seams.R:148

# prequel ----------------------------------------------------------------------
library(testthat)
suppressPackageStartupMessages(library(shiny))
source(test_path("..", "..", "R", "mod_lipid_qc_s4.R"), local = environment())
ensureLipidQcS4WorkflowTestClass <- function() {
    if (!methods::isClass("LipidomicsAssayData")) {
        methods::setClass(
            "LipidomicsAssayData",
            slots = c(
                lipid_data = "list",
                design_matrix = "data.frame",
                group_id = "character",
                lipid_id_column = "character",
                sample_id = "character"
            )
        )
    }
}

# test -------------------------------------------------------------------------
ensureLipidQcS4WorkflowTestClass()
current_s4 <- methods::new(
        "LipidomicsAssayData",
        lipid_data = list(),
        design_matrix = data.frame(sample = "SampleA", condition = "Control"),
        group_id = "condition",
        lipid_id_column = "lipid_id",
        sample_id = "sample"
    )
output_env <- new.env(parent = emptyenv())
captured_state <- NULL
registerLipidQcS4DataSummaryOutput(
        output = output_env,
        workflowData = list(
            state_manager = list(
                getState = function() current_s4
            )
        ),
        renderUiFn = function(expr) htmltools::renderTags(expr)$html,
        reqFn = identity,
        buildDataSummaryUiFn = function(x) {
            captured_state <<- x
            shiny::div("summary shell")
        }
    )
expect_identical(captured_state, current_s4)
