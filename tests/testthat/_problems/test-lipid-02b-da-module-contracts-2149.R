# Extracted from test-lipid-02b-da-module-contracts.R:2149

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
execute_call <- NULL
returned <- handleLipidDaRunAnalysisPreflight(
        analysisContext = list(
            currentS4 = "mock-s4",
            contrastsTbl = data.frame(
                contrasts = "B_vs_A",
                stringsAsFactors = FALSE
            )
        ),
        formulaString = "~ condition",
        daQValThresh = 0.05,
        treatLfcCutoff = 1,
        daData = "mock-da-data",
        workflowData = "mock-workflow-data",
        session = "mock-session",
        experimentPaths = list(
            da_output_dir = "/tmp/da-output",
            publication_graphs_dir = "/tmp/publication-graphs"
        ),
        executeRunAnalysisFn = function(...) {
            execute_call <<- list(...)
            list(marker = "executed")
        }
    )
