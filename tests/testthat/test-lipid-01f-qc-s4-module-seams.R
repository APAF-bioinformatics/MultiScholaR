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

test_that("buildLipidQcS4StateHistoryUi keeps the empty-state prompt stable", {
    rendered <- htmltools::renderTags(
        buildLipidQcS4StateHistoryUi(character(0))
    )$html

    expect_match(rendered, "No processing history available.", fixed = TRUE)
    expect_match(rendered, "circle-info", fixed = TRUE)
})

test_that("buildLipidQcS4StateHistoryUi keeps ordered current-state markup stable", {
    rendered <- htmltools::renderTags(
        buildLipidQcS4StateHistoryUi(c(
            "lipid_imported",
            "lipid_qc_duplicates_resolved",
            "lipid_qc_complete"
        ))
    )$html

    expect_match(rendered, "<ol", fixed = TRUE)
    expect_match(rendered, "1. lipid_imported", fixed = TRUE)
    expect_match(rendered, "2. lipid_qc_duplicates_resolved", fixed = TRUE)
    expect_match(rendered, "3. lipid_qc_complete", fixed = TRUE)
    expect_match(rendered, "(current)", fixed = TRUE)
    expect_match(rendered, "arrow-right", fixed = TRUE)
})

test_that("registerLipidQcS4StateHistoryOutput keeps the history fallback shell stable", {
    output_env <- new.env(parent = emptyenv())

    registerLipidQcS4StateHistoryOutput(
        output = output_env,
        workflowData = list(
            state_manager = list(
                getHistory = function() stop("history unavailable")
            )
        ),
        renderUiFn = function(expr) htmltools::renderTags(expr)$html,
        reqFn = identity
    )

    expect_match(output_env$state_history, "No processing history available.", fixed = TRUE)
})

test_that("buildLipidQcS4DataSummaryUi keeps the missing-state prompt stable", {
    rendered <- htmltools::renderTags(
        buildLipidQcS4DataSummaryUi(NULL)
    )$html

    expect_match(rendered, "No LipidomicsAssayData object available.", fixed = TRUE)
    expect_match(rendered, "triangle-exclamation", fixed = TRUE)
})

test_that("buildLipidQcS4DataSummaryUi keeps the summary table shell stable", {
    ensureLipidQcS4WorkflowTestClass()

    current_s4 <- methods::new(
        "LipidomicsAssayData",
        lipid_data = list(
            AssayA = data.frame(
                lipid_id = c("L1", "L2"),
                SampleA = c(10, 20),
                SampleB = c(30, 40),
                annotation = c("a", "b"),
                stringsAsFactors = FALSE
            ),
            AssayB = data.frame(
                lipid_id = c("L3", "L3"),
                SampleA = c(50, 60),
                SampleB = c(70, 80),
                stringsAsFactors = FALSE
            )
        ),
        design_matrix = data.frame(
            sample = c("SampleA", "SampleB"),
            condition = c("Control", "Treatment"),
            stringsAsFactors = FALSE
        ),
        group_id = "condition",
        lipid_id_column = "lipid_id",
        sample_id = "sample"
    )

    rendered <- htmltools::renderTags(
        buildLipidQcS4DataSummaryUi(current_s4)
    )$html

    expect_match(rendered, "Number of Assays:", fixed = TRUE)
    expect_match(rendered, "Total Lipids:", fixed = TRUE)
    expect_match(rendered, "Number of Samples:", fixed = TRUE)
    expect_match(rendered, "Experimental Groups:", fixed = TRUE)
    expect_match(rendered, "Lipid ID Column:", fixed = TRUE)
    expect_match(rendered, "Sample ID Column:", fixed = TRUE)
    expect_match(rendered, ">2<")
    expect_match(rendered, ">3<")
    expect_match(rendered, "lipid_id", fixed = TRUE)
    expect_match(rendered, "sample", fixed = TRUE)
})

test_that("registerLipidQcS4DataSummaryOutput keeps the builder handoff stable", {
    ensureLipidQcS4WorkflowTestClass()

    current_s4 <- methods::new(
        "LipidomicsAssayData",
        lipid_data = list(),
        design_matrix = data.frame(),
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
    expect_match(output_env$data_summary, "summary shell", fixed = TRUE)
})

test_that("buildLipidQcS4AssayStatsTable keeps the invalid-state fallback stable", {
    expect_null(buildLipidQcS4AssayStatsTable(NULL))
})

test_that("buildLipidQcS4AssayStatsTable keeps the assay-stats shell stable", {
    ensureLipidQcS4WorkflowTestClass()

    current_s4 <- methods::new(
        "LipidomicsAssayData",
        lipid_data = list(
            AssayA = data.frame(
                lipid_id = c("L1", "L2"),
                SampleA = c(10, 0),
                SampleB = c(NA, 40),
                annotation = c("a", "b"),
                stringsAsFactors = FALSE
            ),
            AssayB = data.frame(
                lipid_id = c("L3", "L3"),
                SampleA = c(50, 60),
                SampleB = c(70, 80),
                stringsAsFactors = FALSE
            )
        ),
        design_matrix = data.frame(),
        group_id = "condition",
        lipid_id_column = "lipid_id",
        sample_id = "sample"
    )

    captured <- buildLipidQcS4AssayStatsTable(
        current_s4,
        datatableFn = function(data, options, rownames, class) {
            list(
                data = data,
                options = options,
                rownames = rownames,
                class = class
            )
        }
    )

    expect_named(
        captured$data,
        c("Assay", "Lipids", "Samples", "Missingness")
    )
    expect_identical(captured$data$Assay, c("AssayA", "AssayB"))
    expect_identical(captured$data$Lipids, c(2L, 1L))
    expect_identical(captured$data$Samples, c(2L, 2L))
    expect_equal(captured$data$Missingness, c(50, 0))
    expect_identical(captured$options$dom, "t")
    expect_false(captured$options$paging)
    expect_false(captured$options$ordering)
    expect_false(captured$rownames)
    expect_identical(captured$class, "compact stripe")
})

test_that("registerLipidQcS4AssayStatsOutput keeps the builder handoff stable", {
    ensureLipidQcS4WorkflowTestClass()

    current_s4 <- methods::new(
        "LipidomicsAssayData",
        lipid_data = list(),
        design_matrix = data.frame(),
        group_id = "condition",
        lipid_id_column = "lipid_id",
        sample_id = "sample"
    )
    output_env <- new.env(parent = emptyenv())
    captured_state <- NULL

    registerLipidQcS4AssayStatsOutput(
        output = output_env,
        workflowData = list(
            state_manager = list(
                getState = function() current_s4
            )
        ),
        renderDtFn = identity,
        reqFn = identity,
        buildAssayStatsTableFn = function(x) {
            captured_state <<- x
            "assay stats shell"
        }
    )

    expect_identical(captured_state, current_s4)
    expect_identical(output_env$assay_stats_table, "assay stats shell")
})

test_that("mod_lipid_qc_s4_server routes state-history registration through the seam", {
    source_lines <- readLines(
        test_path("..", "..", "R", "mod_lipid_qc_s4.R"),
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

    seam_calls <- list()
    module_env$registerLipidQcS4StateHistoryOutput <- function(
        output,
        workflowData,
        renderUiFn = shiny::renderUI,
        reqFn = shiny::req,
        buildStateHistoryUiFn = module_env$buildLipidQcS4StateHistoryUi
    ) {
        seam_calls <<- c(seam_calls, list(list(
            output = output,
            workflowData = workflowData,
            renderUiFn = renderUiFn,
            reqFn = reqFn,
            buildStateHistoryUiFn = buildStateHistoryUiFn
        )))
        output
    }

    module_env$mod_lipid_qc_s4_server(
        "s4_finalize",
        workflow_data = list(
            state_manager = list(
                getState = function() NULL,
                getHistory = function() character(0)
            ),
            tab_status = list(),
            config_list = list()
        ),
        omic_type = "lipidomics",
        experiment_label = "Lipidomics"
    )

    session <- shiny:::MockShinySession$new()
    shiny::callModule(module_env$capturedModule, "s4_finalize", session = session)
    session$flushReact()

    expect_length(seam_calls, 1L)
    expect_identical(
        seam_calls[[1]]$workflowData$state_manager$getHistory(),
        character(0)
    )
    expect_identical(
        seam_calls[[1]]$buildStateHistoryUiFn,
        module_env$buildLipidQcS4StateHistoryUi
    )
})

test_that("mod_lipid_qc_s4_server routes data-summary registration through the seam", {
    source_lines <- readLines(
        test_path("..", "..", "R", "mod_lipid_qc_s4.R"),
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

    seam_calls <- list()
    module_env$registerLipidQcS4DataSummaryOutput <- function(
        output,
        workflowData,
        renderUiFn = shiny::renderUI,
        reqFn = shiny::req,
        buildDataSummaryUiFn = module_env$buildLipidQcS4DataSummaryUi
    ) {
        seam_calls <<- c(seam_calls, list(list(
            output = output,
            workflowData = workflowData,
            renderUiFn = renderUiFn,
            reqFn = reqFn,
            buildDataSummaryUiFn = buildDataSummaryUiFn
        )))
        output
    }

    module_env$mod_lipid_qc_s4_server(
        "s4_finalize",
        workflow_data = list(
            state_manager = list(
                getState = function() NULL,
                getHistory = function() character(0)
            ),
            tab_status = list(),
            config_list = list()
        ),
        omic_type = "lipidomics",
        experiment_label = "Lipidomics"
    )

    session <- shiny:::MockShinySession$new()
    shiny::callModule(module_env$capturedModule, "s4_finalize", session = session)
    session$flushReact()

    expect_length(seam_calls, 1L)
    expect_identical(
        seam_calls[[1]]$workflowData$state_manager$getState(),
        NULL
    )
    expect_identical(
        seam_calls[[1]]$buildDataSummaryUiFn,
        module_env$buildLipidQcS4DataSummaryUi
    )
})

test_that("mod_lipid_qc_s4_server routes assay-stats registration through the seam", {
    source_lines <- readLines(
        test_path("..", "..", "R", "mod_lipid_qc_s4.R"),
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

    seam_calls <- list()
    module_env$registerLipidQcS4AssayStatsOutput <- function(
        output,
        workflowData,
        renderDtFn = DT::renderDT,
        reqFn = shiny::req,
        buildAssayStatsTableFn = module_env$buildLipidQcS4AssayStatsTable
    ) {
        seam_calls <<- c(seam_calls, list(list(
            output = output,
            workflowData = workflowData,
            renderDtFn = renderDtFn,
            reqFn = reqFn,
            buildAssayStatsTableFn = buildAssayStatsTableFn
        )))
        output
    }

    module_env$mod_lipid_qc_s4_server(
        "s4_finalize",
        workflow_data = list(
            state_manager = list(
                getState = function() NULL,
                getHistory = function() character(0)
            ),
            tab_status = list(),
            config_list = list()
        ),
        omic_type = "lipidomics",
        experiment_label = "Lipidomics"
    )

    session <- shiny:::MockShinySession$new()
    shiny::callModule(module_env$capturedModule, "s4_finalize", session = session)
    session$flushReact()

    expect_length(seam_calls, 1L)
    expect_identical(
        seam_calls[[1]]$workflowData$state_manager$getState(),
        NULL
    )
    expect_identical(
        seam_calls[[1]]$buildAssayStatsTableFn,
        module_env$buildLipidQcS4AssayStatsTable
    )
})
