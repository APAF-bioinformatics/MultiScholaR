# run-win-k147 lifecycle probe (instrumentation, not a package test): one passing, one failing, one erroring and one
# skipped expectation, so the lane is seen to report each kind through to its SUITE and SUITEFAIL lines. It also
# leaves a folder in the temporary root on purpose, so the SUITETMP listing is seen to report a leftover.
test_that("lifecycle probe: passes", {
    expect_true(exists("writeIntegrationInput", envir = asNamespace("MultiScholaR"), inherits = FALSE))
})
test_that("lifecycle probe: fails on purpose", {
    expect_identical("lifecycle", "probe")
})
test_that("lifecycle probe: errors on purpose", {
    stop("lifecycle probe error on purpose")
})
test_that("lifecycle probe: skips on purpose", {
    skip("lifecycle probe skip on purpose")
})
dir.create(file.path(Sys.getenv("TMPDIR"), "lifecycle-probe-leftover"), showWarnings = FALSE)
