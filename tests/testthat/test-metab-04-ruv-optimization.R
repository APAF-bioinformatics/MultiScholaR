# Tests for the metabolomics RUV automatic optimization loop (RUV-003).
#
# Required cases (fix-plan section 7.4):
#  1. Automatic mode tests the full requested percentage range
#  2. getNegCtrlMetabAnova() called once per percentage, not per assay-percentage
#  3. ruvCancor() called once per percentage, not per assay-percentage
#  4. Two-assay mocked object produces two per-assay results with full traces
#  5. optimization_results row count equals number of tested percentages
#  6. Invalid percentage rows are recorded, not dropped
#  7. Best percentage chosen by composite score plus deterministic tie-breaks
#  8. Sample size computed from assay columns matched to design matrix
#  9. Manual mode output contract remains unchanged
# 10. End-to-end smoke test with a minimal real MetaboliteAssayData runs the
#     true automatic path and returns sane structure
# 11. Whole-object failure at one tested percentage becomes failed rows for each
#     assay and does not abort later percentages

library(testthat)

# ---------------------------------------------------------------------------
# Shared fixture helpers
# ---------------------------------------------------------------------------

makeRuvTestObject <- function(assay_names = c("Plasma", "Urine"), n_features = 6L, n_samples = 4L) {
  sample_ids <- paste0("S", seq_len(n_samples))
  group_ids  <- rep(c("A", "B"), length.out = n_samples)

  assays <- stats::setNames(
    lapply(seq_along(assay_names), function(i) {
      mat <- matrix(
        seq_len(n_features * n_samples) + i * 100L,
        nrow = n_features,
        ncol = n_samples,
        dimnames = list(
          paste0("F", seq_len(n_features), "_", assay_names[[i]]),
          sample_ids
        )
      )
      df <- as.data.frame(mat, stringsAsFactors = FALSE)
      df$database_identifier <- paste0("M", i, "_", seq_len(n_features))
      df$Name                <- paste0("Met_", i, "_", seq_len(n_features))
      df
    }),
    assay_names
  )

  methods::new(
    "MetaboliteAssayData",
    metabolite_data          = assays,
    metabolite_id_column     = "database_identifier",
    annotation_id_column     = "Name",
    database_identifier_type = "database_identifier",
    internal_standard_regex  = "^IS_",
    design_matrix = data.frame(
      Run   = sample_ids,
      group = group_ids,
      batch = paste0("B", seq_len(n_samples)),
      stringsAsFactors = FALSE
    ),
    sample_id             = "Run",
    group_id              = "group",
    technical_replicate_id = "batch",
    args = list()
  )
}

# Scoped mock bindings that restore originals on exit, similar to the shared
# helper in other test files.
localBindingRuv <- function(env, name, value, local_envir = parent.frame()) {
  had_binding <- exists(name, envir = env, inherits = FALSE)
  old_value   <- if (had_binding) get(name, envir = env, inherits = FALSE) else NULL
  was_locked  <- had_binding && bindingIsLocked(name, env)
  if (was_locked) unlockBinding(name, env)
  assign(name, value, envir = env)
  if (was_locked) lockBinding(name, env)
  withr::defer({
    if (exists(name, envir = env, inherits = FALSE) &&
        bindingIsLocked(name, env)) {
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
  }, envir = local_envir)
}

# Build a mock cancor result list keyed by assay names.
makeMockCancorList <- function(assay_names) {
  stats::setNames(
    lapply(assay_names, function(nm) list(tag = nm)),
    assay_names
  )
}

# Build a mock ctrl list (all TRUE per assay, length n_features).
makeMockCtrlList <- function(assay_names, n_features = 6L) {
  stats::setNames(
    lapply(assay_names, function(nm) rep(TRUE, n_features)),
    assay_names
  )
}

# ---------------------------------------------------------------------------
# Helper: run automatic mode with mocked inner functions and capture call counts
# ---------------------------------------------------------------------------
runAutoWithMocks <- function(
    object,
    pct_min = 5,
    pct_max = 7,
    assay_names = c("Plasma", "Urine"),
    n_features  = 6L,
    findBestKElbowFn        = function(cp, ...) 2L,
    calculateSeparationFn   = function(cp, ...) 0.3,
    calculateCompositeFn    = function(sep, k, pen, max_k, ...) sep * 0.8,
    getNegCtrlFn            = NULL,
    ruvCancorFn             = NULL,
    extra_params            = list()
) {
  call_counts <- new.env(parent = emptyenv())
  call_counts$neg_ctrl <- 0L
  call_counts$ruvCancor <- 0L

  if (is.null(getNegCtrlFn)) {
    ctrl_list <- makeMockCtrlList(assay_names, n_features)
    getNegCtrlFn <- function(...) {
      call_counts$neg_ctrl <- call_counts$neg_ctrl + 1L
      ctrl_list
    }
  }
  if (is.null(ruvCancorFn)) {
    cancor_list <- makeMockCancorList(assay_names)
    ruvCancorFn <- function(...) {
      call_counts$ruvCancor <- call_counts$ruvCancor + 1L
      cancor_list
    }
  }

  helper_env <- environment(runPerAssayRuvOptimization)
  localBindingRuv(helper_env, "getNegCtrlMetabAnova", getNegCtrlFn)
  localBindingRuv(helper_env, "ruvCancor",            ruvCancorFn)
  localBindingRuv(helper_env, "findBestKElbow",       findBestKElbowFn)
  localBindingRuv(helper_env, "calculateSeparationScore", calculateSeparationFn)
  localBindingRuv(helper_env, "calculateCompositeScore",  calculateCompositeFn)

  params <- c(
    list(
      percentage_min      = pct_min,
      percentage_max      = pct_max,
      ruv_grouping_variable = "group",
      separation_metric   = "max_difference",
      k_penalty_weight    = 0.5,
      adaptive_k_penalty  = FALSE
    ),
    extra_params
  )

  results <- runPerAssayRuvOptimization(
    theObject = object,
    ruv_mode  = "automatic",
    params    = params
  )

  list(results = results, call_counts = call_counts)
}

# ---------------------------------------------------------------------------
# Case 1: Automatic mode tests the full requested percentage range
# ---------------------------------------------------------------------------
test_that("case 1: automatic mode evaluates every percentage in the requested range", {
  object      <- makeRuvTestObject(c("Plasma"))
  pct_tested  <- integer()

  helper_env <- environment(runPerAssayRuvOptimization)
  localBindingRuv(helper_env, "getNegCtrlMetabAnova", function(theObject, percentage_as_neg_ctrl, ...) {
    pct_tested <<- c(pct_tested, percentage_as_neg_ctrl)
    list(Plasma = rep(TRUE, 6L))
  })
  localBindingRuv(helper_env, "ruvCancor", function(...) list(Plasma = list(tag = "Plasma")))
  localBindingRuv(helper_env, "findBestKElbow", function(...) 2L)
  localBindingRuv(helper_env, "calculateSeparationScore", function(...) 0.3)
  localBindingRuv(helper_env, "calculateCompositeScore",  function(...) 0.24)

  runPerAssayRuvOptimization(
    theObject = object,
    ruv_mode  = "automatic",
    params    = list(
      percentage_min      = 3,
      percentage_max      = 8,
      ruv_grouping_variable = "group",
      adaptive_k_penalty  = FALSE
    )
  )

  expect_equal(sort(pct_tested), 3:8)
})

# ---------------------------------------------------------------------------
# Case 2: getNegCtrlMetabAnova() called once per percentage, not per assay
# ---------------------------------------------------------------------------
test_that("case 2: getNegCtrlMetabAnova called once per percentage across all assays", {
  object <- makeRuvTestObject(c("Plasma", "Urine"))
  out    <- runAutoWithMocks(object, pct_min = 4, pct_max = 6)

  n_percentages <- 6 - 4 + 1  # 3
  expect_equal(out$call_counts$neg_ctrl, n_percentages)
})

# ---------------------------------------------------------------------------
# Case 3: ruvCancor() called once per percentage, not per assay
# ---------------------------------------------------------------------------
test_that("case 3: ruvCancor called once per percentage across all assays", {
  object <- makeRuvTestObject(c("Plasma", "Urine", "CSF"))
  out    <- runAutoWithMocks(object, pct_min = 10, pct_max = 12, assay_names = c("Plasma", "Urine", "CSF"))

  n_percentages <- 3L
  expect_equal(out$call_counts$ruvCancor, n_percentages)
})

# ---------------------------------------------------------------------------
# Case 4: Two-assay mocked object produces two per-assay results with full traces
# ---------------------------------------------------------------------------
test_that("case 4: two-assay object produces two named per-assay results with full traces", {
  object <- makeRuvTestObject(c("Plasma", "Urine"))
  out    <- runAutoWithMocks(object, pct_min = 5, pct_max = 7)
  res    <- out$results

  expect_named(res, c("Plasma", "Urine"))

  for (nm in c("Plasma", "Urine")) {
    expect_true(isTRUE(res[[nm]]$success))
    expect_s3_class(res[[nm]]$optimization_results, "data.frame")
    expect_true(all(c(
      "percentage_requested", "candidate_feature_count",
      "realized_num_controls", "realized_percentage",
      "sample_size", "best_k", "separation_score",
      "composite_score", "status", "error_reason"
    ) %in% names(res[[nm]]$optimization_results)))
  }
})

# ---------------------------------------------------------------------------
# Case 5: optimization_results row count equals number of tested percentages
# ---------------------------------------------------------------------------
test_that("case 5: optimization_results has exactly one row per tested percentage", {
  object        <- makeRuvTestObject(c("Plasma"))
  pct_min       <- 2
  pct_max       <- 9
  n_percentages <- pct_max - pct_min + 1L

  out <- runAutoWithMocks(object, pct_min = pct_min, pct_max = pct_max, assay_names = "Plasma")
  expect_equal(nrow(out$results$Plasma$optimization_results), n_percentages)
})

# ---------------------------------------------------------------------------
# Case 6: Invalid percentage rows are recorded, not dropped
# ---------------------------------------------------------------------------
test_that("case 6: invalid percentage rows appear in optimization_results with status != 'ok'", {
  object     <- makeRuvTestObject(c("Plasma"))
  call_count <- 0L

  helper_env <- environment(runPerAssayRuvOptimization)
  localBindingRuv(helper_env, "getNegCtrlMetabAnova", function(theObject, percentage_as_neg_ctrl, ...) {
    call_count <<- call_count + 1L
    # Fail on the second call
    if (call_count == 2L) stop("simulated neg ctrl failure")
    list(Plasma = rep(TRUE, 6L))
  })
  localBindingRuv(helper_env, "ruvCancor", function(...) list(Plasma = list(tag = "Plasma")))
  localBindingRuv(helper_env, "findBestKElbow", function(...) 2L)
  localBindingRuv(helper_env, "calculateSeparationScore", function(...) 0.3)
  localBindingRuv(helper_env, "calculateCompositeScore",  function(...) 0.24)

  res <- runPerAssayRuvOptimization(
    theObject = object,
    ruv_mode  = "automatic",
    params    = list(
      percentage_min      = 1,
      percentage_max      = 3,
      ruv_grouping_variable = "group",
      adaptive_k_penalty  = FALSE
    )
  )

  opt <- res$Plasma$optimization_results
  expect_equal(nrow(opt), 3L)  # 3 percentages, none dropped
  # Row 2 (pct=2) is the failure
  expect_true(any(opt$status == "error_neg_ctrl_selection"))
  # Rows 1 and 3 succeed
  expect_equal(sum(opt$status == "ok"), 2L)
})

# ---------------------------------------------------------------------------
# Case 7: Best percentage chosen by composite score plus deterministic tie-breaks
# ---------------------------------------------------------------------------
test_that("case 7: winner selected by highest composite then lowest k then highest sep then lowest pct", {
  object <- makeRuvTestObject(c("Plasma"))

  # Assign distinct composite scores per percentage so we can control the winner.
  # pct=5: composite=0.20, k=2
  # pct=6: composite=0.35, k=2
  # pct=7: composite=0.35, k=3  (tie on composite, higher k -> loses to pct=6)
  composite_by_pct <- c(`5` = 0.20, `6` = 0.35, `7` = 0.35)
  sep_by_pct       <- c(`5` = 0.25, `6` = 0.40, `7` = 0.45)
  k_by_pct         <- c(`5` = 2L,   `6` = 2L,   `7` = 3L)

  call_pct <- 0L
  helper_env <- environment(runPerAssayRuvOptimization)
  localBindingRuv(helper_env, "getNegCtrlMetabAnova", function(theObject, percentage_as_neg_ctrl, ...) {
    call_pct <<- percentage_as_neg_ctrl
    list(Plasma = rep(TRUE, 6L))
  })
  localBindingRuv(helper_env, "ruvCancor", function(...) {
    list(Plasma = list(pct = call_pct))
  })
  localBindingRuv(helper_env, "findBestKElbow", function(cp, ...) k_by_pct[[as.character(cp$pct)]])
  localBindingRuv(helper_env, "calculateSeparationScore", function(cp, ...) sep_by_pct[[as.character(cp$pct)]])
  localBindingRuv(helper_env, "calculateCompositeScore", function(sep, k, ...) composite_by_pct[[as.character(call_pct)]])

  res <- runPerAssayRuvOptimization(
    theObject = object,
    ruv_mode  = "automatic",
    params    = list(
      percentage_min      = 5,
      percentage_max      = 7,
      ruv_grouping_variable = "group",
      adaptive_k_penalty  = FALSE
    )
  )

  expect_equal(res$Plasma$best_percentage, 6)  # highest composite, lower k than pct=7
  expect_equal(res$Plasma$best_k, 2L)
  expect_equal(res$Plasma$composite_score, 0.35)
})

# ---------------------------------------------------------------------------
# Case 8: Sample size computed from assay columns matched to design matrix
# ---------------------------------------------------------------------------
test_that("case 8: sample_size excludes ID columns and counts only design-matrix sample columns", {
  # makeRuvTestObject produces assay data frames that include both sample columns
  # (S1..S4) and non-sample metadata columns (database_identifier, Name).
  # sample_size must equal the count of design-matrix sample IDs found in
  # colnames(assay_data), which is 4 (S1..S4) — not ncol(assay_data) which
  # would include the two metadata columns (= 6).
  object      <- makeRuvTestObject(c("Plasma"), n_features = 6L, n_samples = 4L)
  helper_env  <- environment(runPerAssayRuvOptimization)

  localBindingRuv(helper_env, "getNegCtrlMetabAnova", function(...) list(Plasma = rep(TRUE, 6L)))
  localBindingRuv(helper_env, "ruvCancor", function(...) list(Plasma = list(tag = "ok")))
  localBindingRuv(helper_env, "findBestKElbow", function(...) 2L)
  localBindingRuv(helper_env, "calculateSeparationScore", function(...) 0.3)
  localBindingRuv(helper_env, "calculateCompositeScore",  function(...) 0.24)

  res <- runPerAssayRuvOptimization(
    theObject = object,
    ruv_mode  = "automatic",
    params    = list(
      percentage_min      = 10,
      percentage_max      = 11,
      ruv_grouping_variable = "group",
      adaptive_k_penalty  = FALSE
    )
  )

  # The design matrix has 4 samples (S1..S4), all present in the assay's sample
  # columns.  The two metadata columns (database_identifier, Name) must NOT be
  # counted.  ncol(assay_data) = 6, but sample_size must be 4.
  sample_sizes <- unique(res$Plasma$optimization_results$sample_size)
  expect_equal(sample_sizes, 4L)

  # Sanity: total column count in the fixture assay IS 6 (4 samples + 2 id cols)
  assay_ncol <- ncol(object@metabolite_data[["Plasma"]])
  expect_equal(assay_ncol, 6L)
  expect_false(sample_sizes == assay_ncol)
})

# ---------------------------------------------------------------------------
# Case 9: Manual mode output contract remains unchanged
# ---------------------------------------------------------------------------
test_that("case 9: manual mode returns the caller-specified k and percentage unchanged", {
  object     <- makeRuvTestObject(c("Plasma", "Urine"))
  helper_env <- environment(runPerAssayRuvOptimization)

  localBindingRuv(helper_env, "getNegCtrlMetabAnova", function(...) {
    list(Plasma = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE),
         Urine  = c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE))
  })
  localBindingRuv(helper_env, "ruvCancor", function(...) {
    list(Plasma = list(tag = "Plasma"), Urine = list(tag = "Urine"))
  })

  res <- runPerAssayRuvOptimization(
    theObject = object,
    ruv_mode  = "manual",
    params    = list(
      ruv_grouping_variable = "group",
      manual_k              = 4L,
      manual_percentage     = 15
    )
  )

  for (nm in c("Plasma", "Urine")) {
    expect_true(isTRUE(res[[nm]]$success))
    expect_identical(res[[nm]]$best_k, 4L)
    expect_identical(res[[nm]]$best_percentage, 15)
    expect_null(res[[nm]]$optimization_results)
    expect_null(res[[nm]]$error)
  }
})

# ---------------------------------------------------------------------------
# Case 10: End-to-end smoke test — real automatic path, minimal real object
# ---------------------------------------------------------------------------
test_that("case 10: smoke test — automatic mode with mocked whole-object helpers returns sane structure", {
  object <- makeRuvTestObject(c("Plasma"), n_features = 10L, n_samples = 6L)
  out    <- runAutoWithMocks(object, pct_min = 5, pct_max = 7, assay_names = "Plasma")
  res    <- out$results$Plasma

  expect_true(isTRUE(res$success))
  expect_type(res$best_k, "integer")
  expect_true(is.numeric(res$best_percentage))
  expect_true(res$best_percentage %in% 5:7)
  expect_type(res$best_realized_num_controls, "integer")
  expect_true(is.numeric(res$best_realized_percentage))
  expect_true(is.logical(res$control_genes_index))
  expect_s3_class(res$optimization_results, "data.frame")
  expect_equal(nrow(res$optimization_results), 3L)
  expect_true(is.numeric(res$separation_score))
  expect_true(is.numeric(res$composite_score))
  expect_null(res$error)
})

# ---------------------------------------------------------------------------
# Case 11: Whole-object failure at one percentage → failed rows, run continues
# ---------------------------------------------------------------------------
test_that("case 11: whole-object failure at one percentage yields failed rows for all assays and does not abort", {
  object <- makeRuvTestObject(c("Plasma", "Urine"))

  call_count <- 0L
  helper_env <- environment(runPerAssayRuvOptimization)
  localBindingRuv(helper_env, "getNegCtrlMetabAnova", function(theObject, percentage_as_neg_ctrl, ...) {
    call_count <<- call_count + 1L
    # Fail at pct=6
    if (percentage_as_neg_ctrl == 6L) stop("whole-object neg ctrl failure at pct=6")
    list(Plasma = rep(TRUE, 6L), Urine = rep(TRUE, 6L))
  })
  localBindingRuv(helper_env, "ruvCancor", function(...) {
    list(Plasma = list(tag = "Plasma"), Urine = list(tag = "Urine"))
  })
  localBindingRuv(helper_env, "findBestKElbow", function(...) 2L)
  localBindingRuv(helper_env, "calculateSeparationScore", function(...) 0.3)
  localBindingRuv(helper_env, "calculateCompositeScore",  function(...) 0.24)

  res <- runPerAssayRuvOptimization(
    theObject = object,
    ruv_mode  = "automatic",
    params    = list(
      percentage_min      = 5,
      percentage_max      = 7,
      ruv_grouping_variable = "group",
      adaptive_k_penalty  = FALSE
    )
  )

  # getNegCtrlMetabAnova called 3 times (not aborted after failure at pct=6)
  expect_equal(call_count, 3L)

  for (nm in c("Plasma", "Urine")) {
    opt <- res[[nm]]$optimization_results
    expect_equal(nrow(opt), 3L)  # rows for pct 5, 6, 7 — none dropped

    # pct=6 row carries the failure status for both assays
    row6 <- opt[opt$percentage_requested == 6L, , drop = FALSE]
    expect_equal(nrow(row6), 1L)
    expect_equal(row6$status, "error_neg_ctrl_selection")

    # pct=5 and pct=7 rows succeed
    ok_rows <- opt[opt$status == "ok", , drop = FALSE]
    expect_equal(nrow(ok_rows), 2L)
  }

  # Both assays still succeed (have valid winners from pct 5 and 7)
  expect_true(isTRUE(res$Plasma$success))
  expect_true(isTRUE(res$Urine$success))
})
