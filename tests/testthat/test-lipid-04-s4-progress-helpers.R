# fidelity-coverage-compare: shared
library(testthat)

test_that("getFilteringProgressLipidomics initializes the global progress object once", {
    had_progress <- exists("filtering_progress_lipidomics", envir = .GlobalEnv)
    if (had_progress) {
        original_progress <- get("filtering_progress_lipidomics", envir = .GlobalEnv)
    }
    if (exists("filtering_progress_lipidomics", envir = .GlobalEnv)) {
        rm("filtering_progress_lipidomics", envir = .GlobalEnv)
    }
    on.exit({
        if (had_progress) {
            assign("filtering_progress_lipidomics", original_progress, envir = .GlobalEnv)
        } else if (exists("filtering_progress_lipidomics", envir = .GlobalEnv)) {
            rm("filtering_progress_lipidomics", envir = .GlobalEnv)
        }
    }, add = TRUE)

    progress_object <- getFilteringProgressLipidomics()

    expect_s4_class(progress_object, "FilteringProgressLipidomics")
    expect_true(exists("filtering_progress_lipidomics", envir = .GlobalEnv))
    expect_identical(
        progress_object,
        get("filtering_progress_lipidomics", envir = .GlobalEnv)
    )
    expect_identical(getFilteringProgressLipidomics(), progress_object)
})

test_that("FilteringProgressLipidomics exposes the expected slot contract", {
    progress_object <- methods::new("FilteringProgressLipidomics")

    expect_identical(methods::slotNames(progress_object), c(
        "steps",
        "assay_names",
        "n_lipids_per_assay",
        "n_lipids_total",
        "detected_per_sample",
        "missingness_per_assay",
        "sum_intensity_per_sample",
        "cv_distribution_per_assay",
        "is_metrics_per_assay"
    ))
    expect_identical(progress_object@steps, character())
    expect_identical(progress_object@assay_names, list())
    expect_identical(progress_object@n_lipids_per_assay, list())
    expect_identical(progress_object@n_lipids_total, numeric())
})
