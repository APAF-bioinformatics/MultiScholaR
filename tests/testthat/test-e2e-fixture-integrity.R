library(testthat)

# Runs alphabetically before test-e2e-*-browser.R files ('f' < 'p').
# Catches broken fixtures early so browser tests don't silently use bad data.

test_that("E2E fixtures pass structural integrity checks", {
    skip_if(
        !dir.exists(testthat::test_path("..", "testdata", "e2e"))
        , "E2E fixture directory not found — skipping integrity check"
    )

    expect_e2e_fixtures_valid()
})
