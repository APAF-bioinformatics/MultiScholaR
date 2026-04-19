library(testthat)

test_that("log2Transformation keeps the pseudo-count log2 contract", {
  input_matrix <- matrix(
    c(1, 0, 2, NA),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("pep1", "pep2"), c("S1", "S2"))
  )

  pseudo_count <- min(input_matrix[input_matrix > 0], na.rm = TRUE) / 100
  positive_values <- input_matrix > 0 & !is.na(input_matrix)
  expected <- input_matrix
  expected[positive_values] <- expected[positive_values] + pseudo_count
  expected <- log2(expected)

  result <- log2Transformation(input_matrix)

  expect_identical(dimnames(result), dimnames(input_matrix))
  expect_equal(result, expected)
  expect_equal(result["pep1", "S1"], log2(1 + pseudo_count))
  expect_true(is.infinite(result["pep1", "S2"]) && result["pep1", "S2"] < 0)
  expect_true(is.na(result["pep2", "S2"]))
})
