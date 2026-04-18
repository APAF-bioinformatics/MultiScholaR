library(testthat)

repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)

loadSelectedFunctions <- function(paths, symbols, env) {
  for (path in paths[file.exists(paths)]) {
    exprs <- parse(file = path, keep.source = TRUE)

    for (expr in exprs) {
      is_assignment <- is.call(expr) &&
        length(expr) >= 3 &&
        as.character(expr[[1]]) %in% c("<-", "=")
      if (!is_assignment || !is.symbol(expr[[2]])) {
        next
      }

      symbol_name <- as.character(expr[[2]])
      if (symbol_name %in% symbols) {
        eval(expr, envir = env)
      }
    }
  }
}

findSelectedExpression <- function(paths, matcher) {
  for (path in paths[file.exists(paths)]) {
    exprs <- parse(file = path, keep.source = TRUE)

    for (expr in exprs) {
      if (matcher(expr)) {
        return(expr)
      }
    }
  }

  NULL
}

isTargetSymbolAssignment <- function(expr, symbol_name) {
  is.call(expr) &&
    length(expr) >= 3 &&
    as.character(expr[[1]]) %in% c("<-", "=") &&
    is.symbol(expr[[2]]) &&
    identical(as.character(expr[[2]]), symbol_name)
}

target_paths <- c(
  file.path(repo_root, "R", "func_metab_s4_qc_methods.R"),
  file.path(repo_root, "R", "func_metab_s4_objects.R")
)

loadSelectedFunctions(
  paths = target_paths,
  symbols = c("calculateMetabolitePairCorrelation"),
  env = environment()
)

test_that("metabolomics S4 pair-correlation helper preserves pairwise Pearson behavior", {
  input_pair_table <- data.frame(
    metabolite_id = c("M1", "M1", "M2", "M2", "M3", "M3"),
    Run = c("51581", "51582", "51581", "51582", "51581", "51582"),
    abundance = c(1, 2, 2, 4, 3, 6),
    check.names = FALSE
  )

  correlation <- calculateMetabolitePairCorrelation(
    input_pair_table = input_pair_table,
    feature_id_column = "metabolite_id",
    sample_id_column = "Run",
    value_column = "abundance"
  )

  expect_equal(correlation, 1)
})

test_that("metabolomics S4 pair-correlation helper keeps NA fallbacks for invalid and non-finite results", {
  single_sample_table <- data.frame(
    metabolite_id = c("M1", "M2"),
    Run = c("51581", "51581"),
    abundance = c(10, 20),
    check.names = FALSE
  )

  expect_warning(
    invalid_pair_result <- calculateMetabolitePairCorrelation(
      input_pair_table = single_sample_table,
      feature_id_column = "metabolite_id",
      sample_id_column = "Run",
      value_column = "abundance"
    ),
    "does not contain exactly two samples"
  )
  expect_true(is.na(invalid_pair_result))

  constant_pair_table <- data.frame(
    metabolite_id = c("M1", "M1", "M2", "M2", "M3", "M3"),
    Run = c("51581", "51582", "51581", "51582", "51581", "51582"),
    abundance = c(10, 30, 10, 30, 10, 30),
    check.names = FALSE
  )

  non_finite_result <- suppressWarnings(
    calculateMetabolitePairCorrelation(
      input_pair_table = constant_pair_table,
      feature_id_column = "metabolite_id",
      sample_id_column = "Run",
      value_column = "abundance"
    )
  )
  expect_true(is.na(non_finite_result))
})

test_that("metabolomics S4 pair-correlation helper source retains pivot and guard logic", {
  target_expr <- findSelectedExpression(
    paths = target_paths,
    matcher = function(expr) {
      isTargetSymbolAssignment(expr, "calculateMetabolitePairCorrelation")
    }
  )

  expect_false(is.null(target_expr))

  target_text <- paste(deparse(target_expr), collapse = "\n")

  expect_match(target_text, "tidyr::pivot_wider\\(")
  expect_match(target_text, "expected_colnames <- as.character\\(sample_ids\\)")
  expect_match(
    target_text,
    "stats::cor\\(values_x, values_y, use = \"pairwise.complete.obs\"\\)"
  )
  expect_match(target_text, "!is.finite\\(cor_result\\)")
})
