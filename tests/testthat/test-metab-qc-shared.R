# fidelity-coverage-compare: shared
library(testthat)

localNamespaceBinding <- function(env, name, value, .local_envir = parent.frame()) {
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

localNamespaceBindings <- function(env, bindings, .local_envir = parent.frame()) {
  for (name in names(bindings)) {
    localNamespaceBinding(
      env = env,
      name = name,
      value = bindings[[name]],
      .local_envir = .local_envir
    )
  }
}

makeFunctionWithOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

test_that("updateMetaboliteFiltering preserves metrics collection, save-path derivation, and grid return", {
  recorded <- new.env(parent = emptyenv())
  recorded$update_args <- NULL
  recorded$saved <- character()

  save_dir <- tempfile("metab-qc-save-")
  dir.create(save_dir, recursive = TRUE)
  withr::defer(unlink(save_dir, recursive = TRUE, force = TRUE))

  had_project_dirs <- exists("project_dirs", envir = .GlobalEnv, inherits = FALSE)
  old_project_dirs <- if (had_project_dirs) get("project_dirs", envir = .GlobalEnv, inherits = FALSE) else NULL
  had_omic_type <- exists("omic_type", envir = .GlobalEnv, inherits = FALSE)
  old_omic_type <- if (had_omic_type) get("omic_type", envir = .GlobalEnv, inherits = FALSE) else NULL
  had_experiment_label <- exists("experiment_label", envir = .GlobalEnv, inherits = FALSE)
  old_experiment_label <- if (had_experiment_label) get("experiment_label", envir = .GlobalEnv, inherits = FALSE) else NULL
  withr::defer({
    if (had_project_dirs) {
      assign("project_dirs", old_project_dirs, envir = .GlobalEnv)
    } else if (exists("project_dirs", envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = "project_dirs", envir = .GlobalEnv)
    }
    if (had_omic_type) {
      assign("omic_type", old_omic_type, envir = .GlobalEnv)
    } else if (exists("omic_type", envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = "omic_type", envir = .GlobalEnv)
    }
    if (had_experiment_label) {
      assign("experiment_label", old_experiment_label, envir = .GlobalEnv)
    } else if (exists("experiment_label", envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = "experiment_label", envir = .GlobalEnv)
    }
  })

  assign(
    "project_dirs",
    list(metabolomics_demo = list(time_dir = save_dir)),
    envir = .GlobalEnv
  )
  assign("omic_type", "metabolomics", envir = .GlobalEnv)
  assign("experiment_label", "demo", envir = .GlobalEnv)

  update_fun <- makeFunctionWithOverrides(
    updateMetaboliteFiltering,
    list(
      getFilteringProgressMetabolomics = function() "progress-object",
      countUniqueMetabolites = function(current_assay_data, metabolite_id_col) 2L,
      countMetabolitesPerSample = function(current_assay_data, sample_id_col, metabolite_id_col, sample_columns) {
        data.frame(Run = sample_columns, n_detected = c(2L, 1L), stringsAsFactors = FALSE)
      },
      calculateMissingness = function(current_assay_data, sample_id_col, sample_columns) 0.25,
      calculateSumIntensityPerSample = function(current_assay_data, sample_id_col, sample_columns) {
        data.frame(Run = sample_columns, sum_intensity = c(10, 20), stringsAsFactors = FALSE)
      },
      calculateMetaboliteCVs = function(current_assay_data, design_matrix, group_id_col, unused, sample_id_col, metabolite_id_col, sample_columns) {
        data.frame(metabolite_id = "m1", group = "A", cv = 0.1, stringsAsFactors = FALSE)
      },
      getInternalStandardMetrics = function(current_assay_data, is_pattern, metabolite_id_col, sample_id_col, sample_columns) {
        data.frame(is_id = "IS_1", mean_intensity = 5, cv = 0.2, stringsAsFactors = FALSE)
      },
      calculateTotalUniqueMetabolitesAcrossAssays = function(assay_list, metabolite_id_col) 3L,
      updateFilteringProgressMetabolomics = function(prog_met, step_name, assay_names, metrics_list_this_step, total_metabolites, overwrite) {
        recorded$update_args <- list(
          prog_met = prog_met,
          step_name = step_name,
          assay_names = assay_names,
          metrics_list_this_step = metrics_list_this_step,
          total_metabolites = total_metabolites,
          overwrite = overwrite
        )
        invisible(NULL)
      },
      generateMetaboliteFilteringPlots = function(progress) {
        list(
          detected = ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) + ggplot2::geom_point(),
          missingness = ggplot2::ggplot(data.frame(x = 1, y = 2), ggplot2::aes(x, y)) + ggplot2::geom_col()
        )
      },
      ggsave = function(filename, plot, ...) {
        recorded$saved <- c(recorded$saved, filename)
        invisible(filename)
      }
    )
  )

  the_object <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      LC = data.frame(
        metabolite = c("m1", "m2"),
        S1 = c(1, 0),
        S2 = c(2, 3),
        stringsAsFactors = FALSE
      )
    ),
    design_matrix = data.frame(
      Run = c("S1", "S2"),
      group = c("A", "A"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    metabolite_id_column = "metabolite",
    internal_standard_regex = "^IS"
  )

  grid_plot <- update_fun(
    theObject = the_object,
    step_name = "filter_step",
    return_grid = TRUE
  )

  expect_s3_class(grid_plot, "gtable")
  expect_identical(recorded$update_args$prog_met, "progress-object")
  expect_identical(recorded$update_args$step_name, "filter_step")
  expect_identical(recorded$update_args$assay_names, "LC")
  expect_identical(recorded$update_args$total_metabolites, 3L)
  expect_true(any(grepl("filter_step_detected.png", recorded$saved, fixed = TRUE)))
  expect_true(any(grepl("filter_step_combined_plots.png", recorded$saved, fixed = TRUE)))

  expect_error(update_fun(list(), "bad"), "`theObject` must be an S4 object.")
})

test_that("metaboliteIntensityFiltering preserves numeric parsing, per-assay skipping, and helper-driven filtering", {
  package_ns <- asNamespace("MultiScholaR")
  recorded <- new.env(parent = emptyenv())
  recorded$params <- character()
  recorded$thresholds <- list()

  localNamespaceBindings(
    package_ns,
    list(
      checkParamsObjectFunctionSimplify = function(theObject, key, default = NULL) {
        recorded$params <- c(recorded$params, key)
        if (identical(key, "metabolites_intensity_cutoff_percentile")) {
          return("50 # comment")
        }
        if (identical(key, "metabolites_proportion_of_samples_below_cutoff")) {
          return("0.5")
        }
        default
      },
      updateParamInObject = function(theObject, key) theObject,
      metaboliteIntensityFilteringHelper = function(assay_table,
                                                    min_metabolite_intensity_threshold,
                                                    metabolites_proportion_of_samples_below_cutoff,
                                                    metabolite_id_column) {
        recorded$thresholds[[length(recorded$thresholds) + 1L]] <- list(
          threshold = min_metabolite_intensity_threshold,
          cutoff = metabolites_proportion_of_samples_below_cutoff,
          id_col = metabolite_id_column
        )
        assay_table[1, , drop = FALSE]
      }
    )
  )

  the_object <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      valid = data.frame(
        metabolite = c("m1", "m2"),
        S1 = c(10, 1),
        S2 = c(20, 2),
        stringsAsFactors = FALSE
      ),
      no_id = data.frame(
        metabolite = "m3",
        S1 = 5,
        S2 = 7,
        stringsAsFactors = FALSE
      ),
      no_numeric = data.frame(
        metabolite = "m4",
        S1 = "x",
        S2 = "y",
        stringsAsFactors = FALSE
      ),
      all_na = data.frame(
        metabolite = "m5",
        S1 = NA_real_,
        S2 = NA_real_,
        stringsAsFactors = FALSE
      )
    ),
    design_matrix = data.frame(
      Run = c("S1", "S2"),
      group = c("A", "B"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    metabolite_id_column = "metabolite"
  )
  the_object@metabolite_data$no_id <- data.frame(
    other = "m3",
    S1 = 5,
    S2 = 7,
    stringsAsFactors = FALSE
  )

  filtered <- suppressWarnings(metaboliteIntensityFiltering(the_object))

  expect_equal(recorded$params, c(
    "metabolites_intensity_cutoff_percentile",
    "metabolites_proportion_of_samples_below_cutoff"
  ))
  expect_equal(length(recorded$thresholds), 1L)
  expect_identical(unname(recorded$thresholds[[1]]$threshold), 6)
  expect_identical(recorded$thresholds[[1]]$cutoff, 0.5)
  expect_equal(nrow(filtered@metabolite_data$valid), 1L)
  expect_identical(filtered@metabolite_data$valid$metabolite[[1]], "m1")
  expect_identical(filtered@metabolite_data$no_id$other[[1]], "m3")
  expect_identical(filtered@metabolite_data$no_numeric$S1[[1]], "x")
  expect_true(all(is.na(filtered@metabolite_data$all_na$S1)))

  empty_object <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(),
    design_matrix = data.frame(
      Run = c("S1", "S2"),
      group = c("A", "B"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    metabolite_id_column = "metabolite"
  )
  expect_warning(metaboliteIntensityFiltering(empty_object), "no assays")

  localNamespaceBinding(
    package_ns,
    "checkParamsObjectFunctionSimplify",
    function(theObject, key, default = NULL) {
      if (identical(key, "metabolites_intensity_cutoff_percentile")) {
        return("not-a-number")
      }
      "0.5"
    }
  )
  expect_error(
    suppressWarnings(metaboliteIntensityFiltering(the_object)),
    "Failed to convert cleaned metabolites_intensity_cutoff_percentile"
  )
})
