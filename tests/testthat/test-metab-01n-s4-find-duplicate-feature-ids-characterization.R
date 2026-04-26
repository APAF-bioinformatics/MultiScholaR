library(methods)
library(dplyr)
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

if (!methods::isClass("MetaboliteAssayData")) {
  methods::setClass(
    "MetaboliteAssayData",
    slots = c(
      metabolite_data = "list",
      metabolite_id_column = "character",
      design_matrix = "data.frame",
      sample_id = "character"
    ),
    prototype = list(
      metabolite_data = list(),
      metabolite_id_column = "database_identifier",
      design_matrix = data.frame(),
      sample_id = "Sample_ID"
    )
  )
}

target_paths <- c(
  file.path(repo_root, "R", "func_metab_s4_qc_methods.R"),
  file.path(repo_root, "R", "func_metab_s4_objects.R")
)

loadSelectedFunctions(
  paths = target_paths,
  symbols = c("findMetabDuplicateFeatureIDs"),
  env = environment()
)

newMetabDuplicateIdsObject <- function(assay_list, feature_id_col = "database_identifier") {
  methods::new(
    "MetaboliteAssayData",
    metabolite_data = assay_list,
    metabolite_id_column = feature_id_col,
    design_matrix = data.frame(
      Sample_ID = c("Sample_1", "Sample_2"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Sample_ID"
  )
}

test_that("metabolomics S4 duplicate-id helper preserves assay-specific counts and NULL passthroughs", {
  duplicate_ids <- suppressMessages(
    findMetabDuplicateFeatureIDs(
      newMetabDuplicateIdsObject(
        assay_list = list(
          LCMS_Pos = data.frame(
            database_identifier = c("DB_1", "DB_1", "DB_2"),
            intensity = c(10, 20, 30),
            check.names = FALSE
          ),
          LCMS_Neg = data.frame(
            database_identifier = c("DB_3", "DB_4"),
            intensity = c(40, 50),
            check.names = FALSE
          )
        )
      )
    )
  )

  expect_named(duplicate_ids, c("LCMS_Pos", "LCMS_Neg"))
  expect_s3_class(duplicate_ids$LCMS_Pos, "data.frame")
  expect_identical(as.character(duplicate_ids$LCMS_Pos$database_identifier), "DB_1")
  expect_identical(as.integer(duplicate_ids$LCMS_Pos$count), 2L)
  expect_null(duplicate_ids$LCMS_Neg)
})

test_that("metabolomics S4 duplicate-id helper keeps unnamed-assay fallback names", {
  expect_warning(
    duplicate_ids <- suppressMessages(
      findMetabDuplicateFeatureIDs(
        newMetabDuplicateIdsObject(
          assay_list = unname(list(
            data.frame(
              database_identifier = c("DB_1", "DB_1"),
              intensity = c(10, 20),
              check.names = FALSE
            ),
            data.frame(
              database_identifier = c("DB_2", "DB_3"),
              intensity = c(30, 40),
              check.names = FALSE
            )
          ))
        )
      )
    ),
    "Assay list was unnamed"
  )

  expect_named(duplicate_ids, c("Assay_1", "Assay_2"))
  expect_identical(as.character(duplicate_ids$Assay_1$database_identifier), "DB_1")
  expect_null(duplicate_ids$Assay_2)
})

test_that("metabolomics S4 duplicate-id helper preserves missing-column warnings and class guard", {
  expect_error(
    findMetabDuplicateFeatureIDs("bad-object"),
    "Input must be a MetaboliteAssayData object."
  )

  expect_warning(
    duplicate_ids <- suppressMessages(
      findMetabDuplicateFeatureIDs(
        newMetabDuplicateIdsObject(
          assay_list = list(
            Broken = data.frame(
              other_identifier = c("DB_1", "DB_1"),
              intensity = c(10, 20),
              check.names = FALSE
            ),
            Valid = data.frame(
              database_identifier = c("DB_2", "DB_2"),
              intensity = c(30, 40),
              check.names = FALSE
            )
          )
        )
      )
    ),
    "Feature ID column 'database_identifier' not found"
  )

  expect_null(duplicate_ids$Broken)
  expect_identical(as.character(duplicate_ids$Valid$database_identifier), "DB_2")
  expect_identical(as.integer(duplicate_ids$Valid$count), 2L)
})

test_that("metabolomics S4 duplicate-id helper source retains count-filter-naming flow", {
  target_expr <- findSelectedExpression(
    paths = target_paths,
    matcher = function(expr) {
      isTargetSymbolAssignment(expr, "findMetabDuplicateFeatureIDs")
    }
  )

  expect_false(is.null(target_expr))

  target_text <- paste(deparse(target_expr), collapse = "\n")

  expect_match(target_text, "inherits\\(theObject, \"MetaboliteAssayData\"\\)")
  expect_match(target_text, "warning\\(\"Assay list was unnamed")
  expect_match(
    target_text,
    "dplyr::count\\(!!rlang::sym\\(feature_id_col\\),\\s+name = \"count\"\\)"
  )
  expect_match(target_text, "dplyr::filter\\(.data\\$count >\\s+1\\)")
  expect_match(target_text, "purrr::set_names\\(assay_names\\)")
})
