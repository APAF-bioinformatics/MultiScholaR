# fidelity-coverage-compare: shared
library(testthat)

localBinding <- function(env, name, value, .local_envir = parent.frame()) {
  had_binding <- exists(name, envir = env, inherits = FALSE)
  old_value <- if (had_binding) get(name, envir = env, inherits = FALSE) else NULL
  was_locked <- had_binding && bindingIsLocked(name, env)

  if (was_locked) {
    unlockBinding(name, env)
  }
  assign(name, value, envir = env)
  if (was_locked) {
    lockBinding(name, env)
  }

  withr::defer({
    if (exists(name, envir = env, inherits = FALSE) && bindingIsLocked(name, env)) {
      unlockBinding(name, env)
    }
    if (had_binding) {
      assign(name, old_value, envir = env)
    } else if (exists(name, envir = env, inherits = FALSE)) {
      rm(list = name, envir = env)
    }
    if (was_locked && exists(name, envir = env, inherits = FALSE)) {
      lockBinding(name, env)
    }
  }, envir = .local_envir)
}

makeMetabDaObject <- function(assays) {
  methods::new(
    "MetaboliteAssayData",
    metabolite_data = assays,
    metabolite_id_column = "database_identifier",
    annotation_id_column = "Name",
    database_identifier_type = "database_identifier",
    internal_standard_regex = "^IS_",
    design_matrix = data.frame(
      Run = c("S1", "S2"),
      group = c("groupA", "groupB"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    technical_replicate_id = NA_character_,
    args = list()
  )
}

test_that("metabolomics DA orchestration preserves object and contrast validation", {
  expect_error(
    runMetabolitesDA(
      theObject = list(),
      contrasts_tbl = data.frame(contrasts = "groupB-groupA", stringsAsFactors = FALSE)
    ),
    "theObject must be a MetaboliteAssayData S4 object",
    fixed = TRUE
  )

  expect_error(
    runMetabolitesDA(
      theObject = makeMetabDaObject(list()),
      contrasts_tbl = data.frame(other = "groupB-groupA", stringsAsFactors = FALSE)
    ),
    "contrasts_tbl must have a 'contrasts' or 'contrast_string' column",
    fixed = TRUE
  )
})

test_that("metabolomics DA orchestration preserves assay aggregation and skip branches", {
  helper_env <- environment(runMetabolitesDA)
  observed <- new.env(parent = emptyenv())

  localBinding(helper_env, "runTestsContrastsMetabDA", function(
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
  })

  localBinding(helper_env, "createMetabDaResultsLongFormat", function(
    lfc_qval_tbl,
    expr_matrix,
    design_matrix,
    sample_id_col,
    group_id_col,
    metabolite_id_col
  ) {
    observed$group_id_col <- group_id_col
    observed$metabolite_id_col <- metabolite_id_col

    lfc_qval_tbl$intensity.S1.groupA <- expr_matrix[lfc_qval_tbl[[metabolite_id_col]], "S1"]
    lfc_qval_tbl$intensity.S2.groupB <- expr_matrix[lfc_qval_tbl[[metabolite_id_col]], "S2"]
    lfc_qval_tbl
  })

  assay_with_matches <- data.frame(
    database_identifier = c("M1", "M2"),
    Name = c("Met One", "Met Two"),
    S1 = c(10, 11),
    S2 = c(20, 21),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  result <- runMetabolitesDA(
    theObject = makeMetabDaObject(list(LCMS_Pos = assay_with_matches)),
    contrasts_tbl = data.frame(
      contrast_string = "groupB-groupA",
      stringsAsFactors = FALSE
    ),
    da_q_val_thresh = 0.05,
    treat_lfc_cutoff = 0
  )

  expect_identical(observed$sample_cols, c("S1", "S2"))
  expect_identical(observed$contrast_strings, "groupB-groupA")
  expect_identical(observed$sample_id_col, "Run")
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
  expect_identical(result$da_metabolites_long$friendly_name, "groupB-groupA")
  expect_identical(result$da_metabolites_long$metabolite_name, "Met One")
  expect_identical(result$da_metabolites_long$significant, "Up")
  expect_identical(result$da_metabolites_long$intensity.S1.groupA, 10)
  expect_identical(result$da_metabolites_long$intensity.S2.groupB, 20)
  expect_identical(result$significant_counts$LCMS_Pos$up, 1L)
  expect_identical(result$significant_counts$LCMS_Pos$down, 0L)
  expect_identical(result$significant_counts$LCMS_Pos$ns, 0L)
  expect_identical(result$da_q_val_thresh, 0.05)
  expect_identical(result$treat_lfc_cutoff, 0)
})

test_that("metabolomics DA orchestration preserves assay error handling", {
  helper_env <- environment(runMetabolitesDA)

  localBinding(helper_env, "runTestsContrastsMetabDA", function(...) {
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
      qvalue_warnings = list()
    )
  })

  localBinding(helper_env, "createMetabDaResultsLongFormat", function(...) {
    stop("mock long-format failure")
  })

  result <- runMetabolitesDA(
    theObject = makeMetabDaObject(list(
      LCMS_Pos = data.frame(
        database_identifier = c("M1", "M2"),
        Name = c("Met One", "Met Two"),
        S1 = c(10, 11),
        S2 = c(20, 21),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    )),
    contrasts_tbl = data.frame(
      contrasts = "groupB-groupA",
      friendly_names = "B vs A",
      stringsAsFactors = FALSE
    )
  )

  expect_identical(result$per_assay_results$LCMS_Pos, NULL)
  expect_identical(length(result$contrasts_results), 1L)
  expect_equal(nrow(result$da_metabolites_long), 0)
  expect_identical(length(result$significant_counts), 0L)
})
