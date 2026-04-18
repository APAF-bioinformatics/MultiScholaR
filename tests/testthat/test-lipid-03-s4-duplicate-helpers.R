library(testthat)

test_that("resolveDuplicateFeatures delegates through the active duplicate-helper source", {
    wrapper_source_lines <- readLines(file.path("..", "..", "R", "func_lipid_s4_objects.R"))
    helper_source_lines <- readLines(file.path("..", "..", "R", "func_lipid_s4_duplicate_helpers.R"))

    expect_true(any(grepl('setMethod\\("resolveDuplicateFeatures"', wrapper_source_lines)))
    expect_true(any(grepl("resolveDuplicateFeaturesForLipidObject <- function", helper_source_lines, fixed = TRUE)))
    expect_true(any(grepl("resolveDuplicateFeaturesForLipidObject", wrapper_source_lines, fixed = TRUE)))
})

test_that("findLipidDuplicateFeatureIDs reports duplicate IDs per assay", {
    assay_one <- data.frame(
        lipid_id = c("L1", "L1", "L2"),
        S1 = c(10, 40, 5),
        S2 = c(11, 41, 6),
        check.names = FALSE
    )
    assay_two <- data.frame(
        lipid_id = c("L3", "L4"),
        S1 = c(20, 30),
        S2 = c(21, 31),
        check.names = FALSE
    )
    design_matrix <- data.frame(
        Run = c("S1", "S2"),
        group = c("A", "B"),
        stringsAsFactors = FALSE
    )

    lipid_object <- createLipidomicsAssayData(
        lipid_data = list(AssayA = assay_one, AssayB = assay_two),
        design_matrix = design_matrix,
        lipid_id_column = "lipid_id",
        sample_id = "Run",
        group_id = "group"
    )

    duplicates <- findLipidDuplicateFeatureIDs(lipid_object)

    expect_named(duplicates, c("AssayA", "AssayB"))
    expect_s3_class(duplicates$AssayA, "data.frame")
    expect_equal(nrow(duplicates$AssayA), 1)
    expect_identical(duplicates$AssayA$lipid_id, "L1")
    expect_equal(duplicates$AssayA$count, 2)
    expect_null(duplicates$AssayB)
})

test_that("findLipidDuplicateFeatureIDs preserves the empty-assay warning contract", {
    lipid_object <- createLipidomicsAssayData(
        lipid_data = list(),
        design_matrix = data.frame(
            Run = "S1",
            group = "A",
            stringsAsFactors = FALSE
        ),
        sample_id = "Run",
        group_id = "group"
    )

    duplicates <- NULL
    expect_warning(
        duplicates <- findLipidDuplicateFeatureIDs(lipid_object),
        "No assays found"
    )

    expect_equal(duplicates, list())
})

test_that("resolveDuplicateFeaturesByIntensity keeps the highest-average duplicate row", {
    assay_tibble <- data.frame(
        lipid_id = c("L1", "L1", "L2"),
        annotation = c("lower", "higher", "unique"),
        S1 = c(10, 50, 20),
        S2 = c(12, 60, 22),
        S3 = c(14, 70, 24),
        check.names = FALSE
    )

    resolved <- resolveDuplicateFeaturesByIntensity(
        assay_tibble = assay_tibble,
        id_col = "lipid_id",
        sample_cols = c("S1", "S2", "S3")
    )

    expect_equal(nrow(resolved), 2)
    expect_identical(sort(resolved$lipid_id), c("L1", "L2"))
    expect_identical(
        resolved$annotation[resolved$lipid_id == "L1"],
        "higher"
    )
    expect_equal(
        resolved$S1[resolved$lipid_id == "L1"],
        50
    )
})

test_that("resolveDuplicateFeaturesByIntensity preserves the no-sample warning contract", {
    assay_tibble <- data.frame(
        lipid_id = c("L1", "L1"),
        S1 = c(1, 2),
        check.names = FALSE
    )

    resolved <- NULL
    expect_warning(
        resolved <- resolveDuplicateFeaturesByIntensity(
            assay_tibble = assay_tibble,
            id_col = "lipid_id",
            sample_cols = character()
        ),
        "No sample columns provided"
    )

    expect_identical(resolved, assay_tibble)
})

test_that("resolveDuplicateFeatures routes the S4 method through the duplicate helper seam", {
    lipid_object <- createLipidomicsAssayData(
        lipid_data = list(
            AssayA = data.frame(
                lipid_id = c("L1", "L1"),
                lipid = c("feature_one", "feature_two"),
                annotation = c("ann1", "ann2"),
                S1 = c(10, 12),
                S2 = c(11, 15),
                check.names = FALSE
            )
        ),
        design_matrix = data.frame(
            Run = c("S1", "S2"),
            group = c("A", "B"),
            stringsAsFactors = FALSE
        ),
        lipid_id_column = "lipid_id",
        annotation_id_column = "annotation",
        sample_id = "Run",
        group_id = "group"
    )

    had_helper <- exists("resolveDuplicateFeaturesForLipidObject", envir = .GlobalEnv, inherits = FALSE)
    old_helper <- if (had_helper) {
        get("resolveDuplicateFeaturesForLipidObject", envir = .GlobalEnv, inherits = FALSE)
    } else {
        NULL
    }
    on.exit({
        if (had_helper) {
            assign("resolveDuplicateFeaturesForLipidObject", old_helper, envir = .GlobalEnv)
        } else if (exists("resolveDuplicateFeaturesForLipidObject", envir = .GlobalEnv, inherits = FALSE)) {
            rm("resolveDuplicateFeaturesForLipidObject", envir = .GlobalEnv)
        }
    }, add = TRUE)

    captured_call <- NULL
    assign(
        "resolveDuplicateFeaturesForLipidObject",
        function(theObject, itsd_pattern_columns = NULL) {
            captured_call <<- list(
                lipid_id_column = theObject@lipid_id_column,
                itsd_pattern_columns = itsd_pattern_columns
            )
            theObject@args$delegated_duplicate_helper <- TRUE
            theObject
        },
        envir = .GlobalEnv
    )

    resolved_object <- resolveDuplicateFeatures(
        lipid_object,
        itsd_pattern_columns = c("annotation", "lipid")
    )

    expect_s4_class(resolved_object, "LipidomicsAssayData")
    expect_identical(captured_call$lipid_id_column, "lipid_id")
    expect_identical(captured_call$itsd_pattern_columns, c("annotation", "lipid"))
    expect_true(isTRUE(resolved_object@args$delegated_duplicate_helper))
})

test_that("resolveDuplicateFeatures keeps the highest-intensity non-ITSD duplicate rows", {
    lipid_object <- createLipidomicsAssayData(
        lipid_data = list(
            AssayA = data.frame(
                lipid_id = c("L1", "L1", "L2"),
                lipid = c("feature_low", "feature_high", "feature_unique"),
                annotation = c("ann1", "ann2", "ann3"),
                S1 = c(10, 50, 20),
                S2 = c(12, 60, 22),
                S3 = c(14, 70, 24),
                check.names = FALSE
            )
        ),
        design_matrix = data.frame(
            Run = c("S1", "S2", "S3"),
            group = c("A", "A", "B"),
            stringsAsFactors = FALSE
        ),
        lipid_id_column = "lipid_id",
        annotation_id_column = "annotation",
        sample_id = "Run",
        group_id = "group"
    )

    expect_message(
        resolved_object <- resolveDuplicateFeatures(lipid_object),
        "Defaulting to check column specified in slot `annotation_id_column`"
    )

    resolved_assay <- resolved_object@lipid_data$AssayA

    expect_equal(nrow(resolved_assay), 2)
    expect_identical(sort(resolved_assay$lipid_id), c("L1", "L2"))
    expect_identical(
        resolved_assay$lipid[resolved_assay$lipid_id == "L1"],
        "feature_high"
    )
    expect_equal(
        resolved_assay$S1[resolved_assay$lipid_id == "L1"],
        50
    )
})
