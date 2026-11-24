# Helper extraction (KIT-159 7.3, G1(c); section 10.3: T-HELP-1, T-HELP-2).

helperFiles <- function() {
    dir <- withr::local_tempdir(.local_envir = parent.frame())
    writeLines(c(
        "assign('MSRVAL_PLANTED', TRUE, envir = globalenv())",
        "named <- function(x) unnamedHelper(x) + 1",
        "unnamedHelper <- function(x) x * 2",
        "broken <- function() definedNowhere()",
        "notWanted <- function() stop('never extracted')"
    ), file.path(dir, "helper-a.R"))
    writeLines("other <- function() named(1)", file.path(dir, "helper-b.R"))
    file.path(dir, c("helper-a.R", "helper-b.R"))
}

test_that("T-HELP-1: extraction runs no top-level code", {
    if (exists("MSRVAL_PLANTED", envir = globalenv())) rm("MSRVAL_PLANTED", envir = globalenv())
    env <- valExtractHelpers(helperFiles(), "named", parent = baseenv())
    expect_false(exists("MSRVAL_PLANTED", envir = globalenv()))
    expect_identical(env$named(2), 5)
})

test_that("T-HELP-2: the closure is extracted and hashed; a call to nothing defined is refused", {
    files <- helperFiles()
    env <- valExtractHelpers(files, c("named", "other"), parent = baseenv())
    h <- attr(env, "helpers")
    expect_setequal(h$`function`, c("named", "unnamedHelper", "other"))
    expect_false("notWanted" %in% h$`function`)
    expect_true(all(grepl("^[0-9a-f]{64}$", h$sha256_of_source)))
    expect_identical(h$file[h$`function` == "other"], "helper-b.R")
    expect_error(valExtractHelpers(files, "broken", parent = baseenv()),
                 class = "val_driver_helpers_invalid")
    expect_error(valExtractHelpers(files, "absent", parent = baseenv()),
                 class = "val_driver_helpers_invalid")
})
