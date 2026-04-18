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
  symbols = c("resolveDuplicateFeaturesByIntensity"),
  env = environment()
)

test_that("metabolomics S4 duplicate-intensity helper keeps the highest average-intensity feature", {
  assay_tibble <- data.frame(
    database_identifier = c("DB_DUP", "DB_DUP", "DB_KEEP"),
    metabolite = c("low", "high", "unique"),
    Sample_1 = c(5, 50, 2),
    Sample_2 = c(15, 60, 4),
    check.names = FALSE
  )

  resolved <- resolveDuplicateFeaturesByIntensity(
    assay_tibble = assay_tibble,
    id_col = "database_identifier",
    sample_cols = c("Sample_1", "Sample_2")
  )
  resolved <- resolved[order(resolved$database_identifier), , drop = FALSE]

  expect_equal(resolved$database_identifier, c("DB_DUP", "DB_KEEP"))
  expect_equal(resolved$metabolite, c("high", "unique"))
  expect_equal(resolved$Sample_1, c(50, 2))
  expect_equal(resolved$Sample_2, c(60, 4))
})

test_that("metabolomics S4 duplicate-intensity helper preserves passthrough fallbacks", {
  duplicate_free <- data.frame(
    database_identifier = c("DB_1", "DB_2"),
    Sample_1 = c(10, 20),
    Sample_2 = c(30, 40),
    check.names = FALSE
  )

  expect_identical(
    resolveDuplicateFeaturesByIntensity(
      assay_tibble = duplicate_free,
      id_col = "database_identifier",
      sample_cols = c("Sample_1", "Sample_2")
    ),
    duplicate_free
  )

  expect_warning(
    missing_id_result <- resolveDuplicateFeaturesByIntensity(
      assay_tibble = duplicate_free,
      id_col = "missing_id",
      sample_cols = c("Sample_1", "Sample_2")
    ),
    "ID column 'missing_id' not found"
  )
  expect_identical(missing_id_result, duplicate_free)

  expect_warning(
    missing_sample_result <- resolveDuplicateFeaturesByIntensity(
      assay_tibble = duplicate_free,
      id_col = "database_identifier",
      sample_cols = character()
    ),
    "No sample columns provided"
  )
  expect_identical(missing_sample_result, duplicate_free)
})

test_that("metabolomics S4 duplicate-intensity helper source retains ranking-based deduplication", {
  target_expr <- findSelectedExpression(
    paths = target_paths,
    matcher = function(expr) {
      isTargetSymbolAssignment(expr, "resolveDuplicateFeaturesByIntensity")
    }
  )

  expect_false(is.null(target_expr))

  target_text <- paste(deparse(target_expr), collapse = "\n")

  expect_match(target_text, "dplyr::mutate\\(dplyr::across\\(dplyr::any_of\\(sample_cols\\),")
  expect_match(target_text, "as.numeric\\)\\)")
  expect_match(target_text, "dplyr::rowwise\\(\\)")
  expect_match(target_text, "avg_intensity = mean\\(dplyr::c_across\\(dplyr::any_of\\(sample_cols\\)\\),")
  expect_match(target_text, "na.rm = TRUE\\)\\)")
  expect_match(target_text, "dplyr::slice_max\\(order_by = avg_intensity, n = 1, with_ties = FALSE\\)")
  expect_match(target_text, "dplyr::select\\(-avg_intensity\\)")
})
