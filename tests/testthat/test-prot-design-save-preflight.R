# KIT-118 (P-10): the proteomics design save refuses a contrast the design's model cannot hold,
# before any file is written, and lets a valid design through to the persistence step.

protDesignSavePreflightCase <- function(contrasts, design = data.frame(Run = c("S1", "S2", "S3", "S4"),
        group = c("A", "A", "B", "B"), stringsAsFactors = FALSE)) {
    source_dir <- withr::local_tempdir(.local_envir = parent.frame())
    notes <- character()
    steps <- character()
    testthat::local_mocked_bindings(
        persistProtDesignBuilderArtifacts = function(...) {
            steps <<- c(steps, "persist")
            stop("stop after the preflight")
        }, .package = "MultiScholaR", .env = parent.frame())
    testthat::local_mocked_bindings(showNotification = function(ui, ...) {
        notes <<- c(notes, as.character(ui))
        invisible(NULL)
    }, .package = "shiny", .env = parent.frame())
    saved <- suppressMessages(runProtDesignBuilderSaveFlow(
        results = list(design_matrix = design, contrasts_tbl = contrasts, config_list = list()),
        workflowData = new.env(), experimentPaths = list(source_dir = source_dir), session = NULL))
    list(saved = saved, notes = notes, steps = steps, files = list.files(source_dir))
}

test_that("a contrast naming a group absent from the design is refused before anything is written", {
    case <- protDesignSavePreflightCase(data.frame(contrasts = "groupA-groupC", stringsAsFactors = FALSE))
    expect_false(case$saved)
    expect_identical(case$steps, character())
    expect_identical(case$files, character())
    expect_identical(case$notes, paste("Invalid proteomics design: Proteomics DA contrast 'groupA-groupC'",
        "references terms absent from the model matrix: groupC"))
})

test_that("a valid design, with or without contrasts, reaches the persistence step", {
    for (contrasts in list(data.frame(contrasts = "groupA-groupB", stringsAsFactors = FALSE), NULL)) {
        case <- protDesignSavePreflightCase(contrasts)
        expect_identical(case$steps, "persist")
        expect_match(case$notes, "^Error saving design matrix results: stop after the preflight$")
    }
})
