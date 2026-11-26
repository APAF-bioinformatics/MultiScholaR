# KIT-269: formatConfigList lists every vector setting comma-joined and in full, and the study-parameters
# writer's handler writes its failure line into the file. The no-record fixture below is KIT-263's first tests-only
# fixture (8c83cf05): a final object with the shipped configuration as @args and `workflow_data` without a state
# manager, so the writer takes the path without an analysis record.

k269Config <- function() {
    readConfigFile(system.file("config", "config.ini", package = "MultiScholaR"), create_workers = FALSE)
}

k269Protein <- function(config) {
    ProteinQuantitativeData(
        protein_quant_table = data.frame(Protein.Ids = c("P1", "P2", "P3"), S1 = c(10, 11, 12), S2 = c(11, 12, 13),
            S3 = c(20, 21, 22), S4 = c(21, 22, 23), stringsAsFactors = FALSE),
        design_matrix = data.frame(Run = c("S1", "S2", "S3", "S4"), group = c("G1", "G1", "G2", "G2"),
            stringsAsFactors = FALSE),
        sample_id = "Run", group_id = "group", protein_id_column = "Protein.Ids", args = config
    )
}

#' The writer's file lines on the path without an analysis record, with the console lines captured.
k269Lines <- function(config = k269Config()) {
    dir <- withr::local_tempdir(.local_envir = parent.frame())
    utils::capture.output(
        file <- createWorkflowArgsFromConfig(workflow_name = "kit269", description = "d", source_dir_path = dir,
            final_s4_object = k269Protein(config), workflow_data = list(config_list = config))
    )
    readLines(file.path(dir, "study_parameters.txt"), warn = FALSE)
}

#' The lines after the dashes of the "Additional Configuration Parameters" header.
k269After <- function(lines) {
    header <- which(lines == "Additional Configuration Parameters:")
    expect_length(header, 1L)
    lines[-seq_len(header + 1L)]
}

test_that("(KIT-269 T1) a numeric pair is listed comma-joined", {
    expect_identical(formatConfigList(list(plotRle = list(yaxis_limit = c(-3, 3)))),
        c("plotRle:", "  Yaxis Limit: -3, 3"))
})

test_that("(KIT-269 T2) integer, logical and double vectors, and the RUV control shape, are listed in full", {
    expect_identical(formatConfigList(list(alpha = c(1L, 2L, 3L))), "Alpha: 1, 2, 3")
    expect_identical(formatConfigList(list(alpha = c(TRUE, FALSE))), "Alpha: TRUE, FALSE")
    expect_identical(formatConfigList(list(alpha = c(0.5, 1.5))), "Alpha: 0.5, 1.5")
    expect_identical(
        formatConfigList(list(ctrl = list(assay_pos = c(TRUE, FALSE), assay_neg = c(TRUE, TRUE, FALSE)))),
        c("Ctrl:", "  Assay Pos: TRUE, FALSE", "  Assay Neg: TRUE, TRUE, FALSE")
    )
    # The real shape of the RUV control vector: one logical per feature, named by the feature id (as the analysis
    # record prints it, `id=value`), so the controls can be identified from the file.
    ctrl <- c(Feat_A = TRUE, Feat_B = FALSE, Feat_C = TRUE)
    expect_identical(formatConfigList(list(ctrl = list(assay_pos = ctrl))),
        c("Ctrl:", "  Assay Pos: Feat_A=TRUE, Feat_B=FALSE, Feat_C=TRUE"))
})

test_that("(KIT-269 T3) scalars, pairs, empty lists and skipped names are unchanged; text vectors print in full", {
    expect_identical(formatConfigList(list(alpha = 0.7)), "Alpha: 0.7")
    expect_identical(formatConfigList(list(alpha = c("pdf", "png"))), "Alpha: pdf, png")
    expect_identical(formatConfigList(list(alpha = list())), "Alpha: [Empty List]")
    expect_identical(formatConfigList(list(internal_workflow_source_dir = "/x", alpha = 1)), "Alpha: 1")
    seven <- paste0("v", 1:7)
    expect_identical(formatConfigList(list(alpha = seven)), paste0("Alpha: ", paste(seven, collapse = ", ")))
})

test_that("(KIT-269 T4) the shipped configuration is listed, the vector and the nine column ids in full", {
    lines <- formatConfigList(k269Config())
    expect_true("  Yaxis Limit: -3, 3" %in% lines)
    ids <- grep("^ +Input Matrix Column Ids: ", lines, value = TRUE)
    expect_length(ids, 1L)
    expect_identical(sub("^ +Input Matrix Column Ids: ", "", ids), paste(
        c("Run", "Precursor.Id", "Protein.Group", "Protein.Ids", "Stripped.Sequence", "Modified.Sequence",
            "Precursor.Charge", "Precursor.Quantity", "Precursor.Normalised"), collapse = ", "))
})

test_that("(KIT-269 T5) the writer without an analysis record lists the settings under the Additional header", {
    lines <- k269Lines()
    expect_false(any(lines == "Analysis Record:"))
    after <- k269After(lines)
    cat(sprintf("\nKIT-269 T5: %d lines, the first six after the Additional dashes:\n%s\n", length(lines),
        paste(utils::head(after, 6L), collapse = "\n")))
    expect_true("  Yaxis Limit: -3, 3" %in% after)
    expect_false(any(grepl("Error formatting config", lines, fixed = TRUE)))
})

test_that("(KIT-269 T6) a listing that fails says so in the file, with the condition's message on one line", {
    run <- function(fault) {
        local_mocked_bindings(formatConfigList = function(...) fault(), .package = "MultiScholaR")
        dir <- withr::local_tempdir(.local_envir = parent.frame())
        config <- k269Config()
        console <- utils::capture.output(
            file <- createWorkflowArgsFromConfig(workflow_name = "kit269", description = "d", source_dir_path = dir,
                final_s4_object = k269Protein(config), workflow_data = list(config_list = config))
        )
        list(file = file, console = console, lines = readLines(file.path(dir, "study_parameters.txt"), warn = FALSE))
    }
    line <- function(lines) {
        expect_false(any(lines == "Analysis Record:"))
        k269After(lines)[1L]
    }
    boom <- run(function() stop("boom"))
    expect_true(any(grepl("WORKFLOW ARGS: formatConfigList failed: boom", boom$console, fixed = TRUE)))
    expect_true(file.exists(boom$file))
    expect_identical(line(boom$lines), "The configuration settings could not be listed here: boom")
    parent <- run(function() rlang::abort("outer", parent = simpleError("inner cause")))
    expect_match(line(parent$lines), "^The configuration settings could not be listed here: .*inner cause")
    expect_identical(sum(lengths(regmatches(parent$lines, gregexpr("inner cause", parent$lines, fixed = TRUE)))), 1L)
    multi <- run(function() stop("line one\nline   two"))
    expect_identical(line(multi$lines), "The configuration settings could not be listed here: line one line two")
})

test_that("(KIT-269 T7) createStudyParametersFile writes the file with the shipped configuration", {
    config_env <- new.env(parent = emptyenv())
    config_env$config_list <- k269Config()
    dir <- withr::local_tempdir()
    study_fun <- local({
        fun <- createStudyParametersFile
        environment(fun) <- list2env(list(requireNamespace = function(package, quietly = TRUE) FALSE),
            parent = environment(createStudyParametersFile))
        fun
    })
    file <- suppressMessages(study_fun(workflow_name = "kit269", description = "d", source_dir_path = dir,
        config_list_name = "config_list", env = config_env))
    expect_identical(basename(file), "study_parameters.txt")
    expect_true("  Yaxis Limit: -3, 3" %in% readLines(file, warn = FALSE))
})
