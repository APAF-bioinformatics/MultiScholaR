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

makeMetabNormObject <- function(assay_names = c("Plasma", "Urine")) {
  assays <- stats::setNames(lapply(seq_along(assay_names), function(i) {
    data.frame(
      database_identifier = c(paste0("M", i, "_1"), paste0("M", i, "_2")),
      Name = c(paste0("Met_", i, "_1"), paste0("Met_", i, "_2")),
      S1 = c(10 + i, 20 + i),
      S2 = c(30 + i, 40 + i),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }), assay_names)

  methods::new(
    "MetaboliteAssayData",
    metabolite_data = assays,
    metabolite_id_column = "database_identifier",
    annotation_id_column = "Name",
    database_identifier_type = "database_identifier",
    internal_standard_regex = "^IS_",
    design_matrix = data.frame(
      Run = c("S1", "S2"),
      group = c("A", "B"),
      batch = c("B1", "B2"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    technical_replicate_id = "batch",
    args = list()
  )
}

test_that("metabolomics standalone QC plot helper preserves per-assay output routing", {
  object <- makeMetabNormObject(c("LCMS Pos", "GCMS"))
  saved_paths <- character()

  localBinding(
    asNamespace("ggplot2"),
    "ggsave",
    function(filename, plot, ...) {
      saved_paths <<- c(saved_paths, basename(filename))
      invisible(filename)
    }
  )

  helper_env <- environment(generateMetabQcPlots)
  localBinding(helper_env, "plotPca", function(...) {
    list(
      `LCMS Pos` = ggplot2::ggplot(),
      GCMS = ggplot2::ggplot()
    )
  })
  localBinding(helper_env, "plotRle", function(...) ggplot2::ggplot())
  localBinding(helper_env, "plotDensity", function(...) list(ggplot2::ggplot()))
  localBinding(helper_env, "plotPearson", function(...) {
    list(
      `LCMS Pos` = ggplot2::ggplot(),
      GCMS = ggplot2::ggplot()
    )
  })

  output <- generateMetabQcPlots(
    theObject = object,
    experiment_paths = list(metabolite_qc_dir = tempdir()),
    stage = "post_filter",
    shape_variable = "batch"
  )

  expect_length(output, 8L)
  expect_equal(
    sort(saved_paths),
    sort(c(
      "gcms_pre_norm_correlation.png",
      "gcms_pre_norm_density.png",
      "gcms_pre_norm_pca.png",
      "gcms_pre_norm_rle.png",
      "lcms_pos_pre_norm_correlation.png",
      "lcms_pos_pre_norm_density.png",
      "lcms_pos_pre_norm_pca.png",
      "lcms_pos_pre_norm_rle.png"
    ))
  )

  empty_object <- makeMetabNormObject(character())
  expect_identical(
    generateMetabQcPlots(
      theObject = empty_object,
      experiment_paths = list(metabolite_qc_dir = tempdir())
    ),
    list()
  )
})

test_that("metabolomics standalone RUV optimization helper preserves automatic, manual, and error outputs", {
  object <- makeMetabNormObject(c("Plasma", "Urine"))
  helper_env <- environment(runPerAssayRuvOptimization)

  localBinding(helper_env, "getNegCtrlMetabAnova", function(...) {
    list(
      Plasma = c(TRUE, FALSE),
      Urine = c(FALSE, TRUE)
    )
  })
  localBinding(helper_env, "ruvCancor", function(...) {
    list(
      Plasma = list(tag = "Plasma"),
      Urine = list(tag = "Urine")
    )
  })

  localBinding(helper_env, "findBestK", function(cancor_plot) {
    if (identical(cancor_plot$tag, "Urine")) {
      stop("mock cancor failure")
    }
    3L
  })
  localBinding(helper_env, "calculateSeparationScore", function(cancor_plot, metric = "max_difference") {
    if (identical(cancor_plot$tag, "Plasma")) {
      0.42
    } else {
      0.17
    }
  })

  automatic <- runPerAssayRuvOptimization(
    theObject = object,
    ruv_mode = "automatic",
    params = list(
      percentage_min = 2,
      percentage_max = 12,
      ruv_grouping_variable = "missing_column",
      separation_metric = "mean_difference",
      k_penalty_weight = 0.75,
      adaptive_k_penalty = FALSE
    )
  )

  expect_true(isTRUE(automatic$Plasma$success))
  expect_identical(automatic$Plasma$best_k, 3L)
  expect_identical(automatic$Plasma$best_percentage, 12)
  expect_equal(automatic$Plasma$separation_score, 0.42)
  expect_false(automatic$Urine$success)
  expect_match(automatic$Urine$error, "mock cancor failure", fixed = TRUE)

  localBinding(helper_env, "findBestK", function(cancor_plot) 2L)

  manual <- runPerAssayRuvOptimization(
    theObject = object,
    ruv_mode = "manual",
    params = list(
      ruv_grouping_variable = "batch",
      manual_k = 5,
      manual_percentage = 18
    )
  )

  expect_true(isTRUE(manual$Plasma$success))
  expect_true(isTRUE(manual$Urine$success))
  expect_identical(manual$Plasma$best_k, 5)
  expect_identical(manual$Plasma$best_percentage, 18)
  expect_true(is.null(manual$Plasma$optimization_results))
})
