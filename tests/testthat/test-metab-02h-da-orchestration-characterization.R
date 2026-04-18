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

if (!methods::isClass("MetaboliteAssayData")) {
  methods::setClass(
    "MetaboliteAssayData",
    slots = c(
      metabolite_data = "list",
      design_matrix = "data.frame",
      sample_id = "character",
      group_id = "character",
      metabolite_id_column = "character",
      annotation_id_column = "character"
    )
  )
}

loadSelectedFunctions(
  paths = c(
    file.path(repo_root, "R", "func_metab_da_orchestration.R"),
    file.path(repo_root, "R", "func_metab_da.R")
  ),
  symbols = "runMetabolitesDA",
  env = environment()
)

test_that("metabolomics DA orchestration preserves the input object guard", {
  expect_error(
    runMetabolitesDA(
      theObject = list(),
      contrasts_tbl = data.frame(contrasts = "groupB-groupA", stringsAsFactors = FALSE)
    ),
    "theObject must be a MetaboliteAssayData S4 object"
  )
})

test_that("metabolomics DA orchestration preserves assay aggregation and warning capture", {
  observed <- new.env(parent = emptyenv())

  environment(runMetabolitesDA)$runTestsContrastsMetabDA <- function(
    data,
    contrast_strings,
    design_matrix,
    formula_string,
    sample_id_col = "Run",
    treat_lfc_cutoff = NA,
    eBayes_trend = TRUE,
    eBayes_robust = TRUE
  ) {
    observed$sample_cols <- colnames(data)
    observed$contrast_strings <- contrast_strings
    observed$sample_id_col <- sample_id_col
    observed$formula_string <- formula_string
    observed$treat_lfc_cutoff <- treat_lfc_cutoff

    list(
      results = list(
        "groupB-groupA" = data.frame(
          logFC = 1.5,
          P.Value = 0.001,
          fdr_qvalue = 0.01,
          fdr_value_bh = 0.01,
          raw_pvalue = 0.001,
          row.names = "M1",
          check.names = FALSE
        )
      ),
      fit.eb = structure(list(marker = "fit"), class = "mock_fit"),
      qvalue_warnings = list("groupB-groupA" = TRUE)
    )
  }

  environment(runMetabolitesDA)$createMetabDaResultsLongFormat <- function(
    lfc_qval_tbl,
    expr_matrix,
    design_matrix,
    sample_id_col,
    group_id_col,
    metabolite_id_col
  ) {
    observed$group_id_col <- group_id_col
    observed$metabolite_id_col <- metabolite_id_col

    lfc_qval_tbl$intensity.Sample_1.groupA <- expr_matrix[lfc_qval_tbl[[metabolite_id_col]], "Sample_1"]
    lfc_qval_tbl$intensity.Sample_2.groupB <- expr_matrix[lfc_qval_tbl[[metabolite_id_col]], "Sample_2"]
    lfc_qval_tbl
  }

  assay_with_matches <- data.frame(
    metabolite_id = c("M1", "M2"),
    metabolite_name = c("Met One", "Met Two"),
    Sample_1 = c(10, 11),
    Sample_2 = c(20, 21),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  assay_without_matches <- data.frame(
    metabolite_id = c("N1"),
    metabolite_name = c("No Match"),
    Other_1 = 5,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  design_matrix <- data.frame(
    sample_id = c("Sample_1", "Sample_2"),
    group = c("groupA", "groupB"),
    stringsAsFactors = FALSE
  )

  theObject <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      LCMS_Pos = assay_with_matches,
      LCMS_Neg = assay_without_matches
    ),
    design_matrix = design_matrix,
    sample_id = "sample_id",
    group_id = "group",
    metabolite_id_column = "metabolite_id",
    annotation_id_column = "metabolite_name"
  )

  result <- runMetabolitesDA(
    theObject = theObject,
    contrasts_tbl = data.frame(
      contrasts = "groupB-groupA",
      friendly_names = "B vs A",
      stringsAsFactors = FALSE
    ),
    da_q_val_thresh = 0.05,
    treat_lfc_cutoff = 0
  )

  expect_identical(observed$sample_cols, c("Sample_1", "Sample_2"))
  expect_identical(observed$contrast_strings, "groupB-groupA")
  expect_identical(observed$sample_id_col, "sample_id")
  expect_identical(observed$formula_string, "~ 0 + group")
  expect_true(is.na(observed$treat_lfc_cutoff))
  expect_identical(observed$group_id_col, "group")
  expect_identical(observed$metabolite_id_col, "metabolite_id")

  expect_identical(names(result$contrasts_results), "LCMS_Pos")
  expect_identical(names(result$per_assay_results), "LCMS_Pos")
  expect_identical(names(result$significant_counts), "LCMS_Pos")
  expect_identical(names(result$qvalue_warnings), "LCMS_Pos")
  expect_true(result$qvalue_warnings$LCMS_Pos[["groupB-groupA"]])

  expect_equal(nrow(result$da_metabolites_long), 1)
  expect_identical(result$da_metabolites_long$assay, "LCMS_Pos")
  expect_identical(result$da_metabolites_long$comparison, "groupB-groupA")
  expect_identical(result$da_metabolites_long$friendly_name, "B vs A")
  expect_identical(result$da_metabolites_long$metabolite_name, "Met One")
  expect_identical(result$da_metabolites_long$significant, "Up")
  expect_identical(result$da_metabolites_long$intensity.Sample_1.groupA, 10)
  expect_identical(result$da_metabolites_long$intensity.Sample_2.groupB, 20)

  expect_identical(result$significant_counts$LCMS_Pos$up, 1L)
  expect_identical(result$significant_counts$LCMS_Pos$down, 0L)
  expect_identical(result$significant_counts$LCMS_Pos$ns, 0L)
  expect_identical(result$da_q_val_thresh, 0.05)
  expect_identical(result$treat_lfc_cutoff, 0)
})
