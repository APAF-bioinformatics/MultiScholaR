# The coverage table's ends and failing condition (KIT-159 5.7, r1.12; build-points2 point 3):
# coverage_min, coverage_max over the generator conditions and, for a row not offered,
# failed_condition and failed_side, computed from the recorded detail columns.

covRow <- function(condition, coverage, kind = "approx_within", level = 0.95) {
    data.frame(cell = paste0("C6-", condition), regime = "C6", family = "newcombe_paired", n = 21L,
               condition = condition, statistic = "newcombe_paired", class = NA_character_,
               level = level, nsim = 2000L, coverage = coverage, mc_se = 0.003,
               achieved = NA_real_, undefined = 0, target = 0,
               status = valCoverageStatus(kind, coverage, level), kind = kind,
               stringsAsFactors = FALSE)
}

test_that("T-COV-ENDS: both ends and the failing condition and side are in the regime row", {
    tab <- valSummariseCoverage(rbind(covRow("p0.9_0.8", 0.9765), covRow("p0.9_0.9", 0.9925)))
    expect_identical(tab$status, "fail")
    expect_equal(tab$coverage_min, 0.9765)
    expect_equal(tab$coverage_max, 0.9925)
    expect_identical(tab$failed_condition, "p0.9_0.9")
    expect_identical(tab$failed_side, "high")
    low <- valSummariseCoverage(rbind(covRow("a", 0.90, "approx_run"), covRow("b", 0.95, "approx_run")))
    expect_identical(low$failed_condition, "a")
    expect_identical(low$failed_side, "low")
    ok <- valSummariseCoverage(rbind(covRow("a", 0.95), covRow("b", 0.96)))
    expect_true(is.na(ok$failed_condition))
    expect_true(is.na(ok$failed_side))
})

test_that("T-COV-ENDS-2: the committed T2 table is what the summary gives on the committed detail", {
    harness <- Sys.getenv("VAL_HARNESS_DIR")
    skip_if(!nzchar(harness), "VAL_HARNESS_DIR not set")
    t2 <- file.path(harness, "..", "t2")
    detail <- do.call(rbind, lapply(c("a", "b", "c"), function(r) {
        utils::read.delim(file.path(t2, sprintf("t2%s-coverage-detail.tsv", r)),
                          stringsAsFactors = FALSE)
    }))
    got <- valSummariseCoverage(detail)
    want <- utils::read.delim(file.path(t2, "coverage-table.tsv"), stringsAsFactors = FALSE)
    key <- function(x) paste(x$regime, x$class, x$design, x$transform, x$n, sprintf("%.6f", x$level))
    ends <- c("coverage_min", "coverage_max", "failed_condition", "failed_side")
    expect_true(all(ends %in% names(want)))
    expect_true(all(ends %in% names(got)))
    got <- got[match(key(want), key(got)), ]
    expect_identical(got$status, want$status)
    expect_equal(got$coverage_min, want$coverage_min)
    expect_equal(got$coverage_max, want$coverage_max)
    expect_identical(ifelse(is.na(got$failed_condition), "", got$failed_condition),
                     ifelse(is.na(want$failed_condition), "", want$failed_condition))
    expect_identical(ifelse(is.na(got$failed_side), "", got$failed_side),
                     ifelse(is.na(want$failed_side), "", want$failed_side))
})
