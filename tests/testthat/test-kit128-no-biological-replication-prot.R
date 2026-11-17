# KIT-128 K128-1 (the user's rule, 2026-11-17): a design whose every group is one biological unit measured more than once
# (the rule's `spanned_by_design`) is refused in every omic, proteomics included: no biological replication.
library(testthat)

kit128ProtData <- function() {
  data <- matrix(c(10, 10.2, 10.5, 10.4, 12, 12.1, 11.8, 11.9, 15, 14.8, 15.1, 15.3, 13, 13.2, 12.9, 13.1), nrow = 2L,
    byrow = TRUE, dimnames = list(c("P1", "P2"), paste0("S", 1:8)))
  data
}

kit128ProtDesign <- function(replicates) {
  design <- data.frame(Run = paste0("S", 1:8), group = rep(c("A", "B"), each = 4L), replicates = replicates,
    stringsAsFactors = FALSE)
  rownames(design) <- design$Run
  design
}

test_that("KIT-128 K128-1: proteomics refuses a design whose every group is one unit measured four times", {
  expect_error(suppressMessages(runTestsContrasts(data = kit128ProtData(), contrast_strings = "groupB-groupA",
    design_matrix = kit128ProtDesign(c(1, 1, 1, 1, 2, 2, 2, 2)), formula_string = "~ 0 + group")),
    "groups A, B", class = "omics_da_no_biological_replication")
  # Two biological units per group, each measured twice: fitted with the block, as before.
  results <- suppressMessages(runTestsContrasts(data = kit128ProtData(), contrast_strings = "groupB-groupA",
    design_matrix = kit128ProtDesign(c(1, 1, 2, 2, 3, 3, 4, 4)), formula_string = "~ 0 + group"))
  expect_identical(nrow(results$results[[1L]]), 2L)
})
