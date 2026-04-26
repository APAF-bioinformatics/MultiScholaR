# fidelity-coverage-compare: shared
library(testthat)

buildMetabDaOutputResultsShared <- function() {
  data.frame(
    metabolite_id = c("M1", "M2", "M3", "M4"),
    metabolite_name = c("Met One", "Met Two", "Met Three", "Met Four"),
    logFC = c(1.8, -0.7, 2.1, -1.6),
    raw_pvalue = c(0.001, 0.20, 0.003, 0.02),
    fdr_qvalue = c(0.01, 0.20, 0.02, 0.04),
    significant = c("Up", "NS", "Up", "Down"),
    comparison = c("groupB-groupA", "groupB-groupA", "groupC-groupA", "groupC-groupA"),
    friendly_name = c(
      "groupB-groupA = B vs A",
      "groupB-groupA = B vs A",
      "groupC-groupA = C vs A",
      "groupC-groupA = C vs A"
    ),
    numerator = c("groupB", "groupB", "groupC", "groupC"),
    denominator = c("groupA", "groupA", "groupA", "groupA"),
    assay = c("LCMS_Pos", "LCMS_Pos", "LCMS_Neg", "LCMS_Neg"),
    intensity.S1.groupA = c(10, 12, 9, 15),
    intensity.S2.groupB = c(20, 25, 18, 30),
    stringsAsFactors = FALSE
  )
}

localMetabDaOutputBinding <- function(env, name, value, .local_envir = parent.frame()) {
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

test_that("metabolomics DA results output writes tables, volcanoes, and summary artifacts", {
  tmp_root <- tempfile("metab-da-output-shared-")
  da_output_dir <- file.path(tmp_root, "de-results")
  publication_graphs_dir <- file.path(tmp_root, "publication")
  dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)
  withr::defer(unlink(tmp_root, recursive = TRUE, force = TRUE))

  package_ns <- asNamespace("MultiScholaR")
  localMetabDaOutputBinding(
    package_ns,
    "generateMetabDAVolcanoStatic",
    function(...) {
      ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
        ggplot2::geom_point()
    }
  )
  localMetabDaOutputBinding(
    package_ns,
    "generateMetabDAHeatmap",
    function(...) NULL
  )

  expect_true(outputMetabDaResultsAllContrasts(
    da_results_list = list(da_metabolites_long = buildMetabDaOutputResultsShared()),
    da_output_dir = da_output_dir,
    publication_graphs_dir = publication_graphs_dir,
    da_q_val_thresh = 0.05,
    lfc_threshold = 1
  ))

  expect_true(file.exists(file.path(
    da_output_dir,
    "de_posmode_metabolites_groupB-groupA_long_annot.tsv"
  )))
  expect_true(file.exists(file.path(
    da_output_dir,
    "de_negmode_metabolites_groupC-groupA_long_annot.xlsx"
  )))
  expect_true(file.exists(file.path(
    publication_graphs_dir,
    "Volcano_Plots",
    "posmode_groupB-groupA.png"
  )))
  expect_true(file.exists(file.path(
    publication_graphs_dir,
    "Volcano_Plots",
    "negmode_groupC-groupA.pdf"
  )))
  expect_true(file.exists(file.path(
    publication_graphs_dir,
    "Volcano_Plots",
    "all_volcano_plots_combined.pdf"
  )))
  expect_false(file.exists(file.path(
    publication_graphs_dir,
    "Heatmaps",
    "all_heatmaps_combined.pdf"
  )))

  numsig_tsv <- file.path(
    publication_graphs_dir,
    "NumSigDeMolecules",
    "metabolites_num_sig_de_molecules.tab"
  )
  expect_true(file.exists(numsig_tsv))
  numsig <- read.delim(numsig_tsv, check.names = FALSE)
  expect_setequal(numsig$mode, c("posmode", "negmode"))
  expect_equal(numsig$significant, c(1, 2))
  expect_true(file.exists(file.path(
    publication_graphs_dir,
    "NumSigDeMolecules",
    "metabolites_num_sig_barplot.png"
  )))
})

test_that("metabolomics DA results output covers heatmap rendering and empty-result exit", {
  tmp_root <- tempfile("metab-da-output-heatmap-")
  da_output_dir <- file.path(tmp_root, "de-results")
  publication_graphs_dir <- file.path(tmp_root, "publication")
  dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)
  withr::defer(unlink(tmp_root, recursive = TRUE, force = TRUE))

  package_ns <- asNamespace("MultiScholaR")
  localMetabDaOutputBinding(
    package_ns,
    "generateMetabDAVolcanoStatic",
    function(...) {
      ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
        ggplot2::geom_point()
    }
  )
  localMetabDaOutputBinding(
    package_ns,
    "generateMetabDAHeatmap",
    function(...) {
      ComplexHeatmap::Heatmap(matrix(c(1, 2, 3, 4), nrow = 2))
    }
  )

  expect_true(outputMetabDaResultsAllContrasts(
    da_results_list = list(da_metabolites_long = buildMetabDaOutputResultsShared()),
    da_output_dir = da_output_dir,
    publication_graphs_dir = publication_graphs_dir
  ))

  expect_true(file.exists(file.path(
    publication_graphs_dir,
    "Heatmaps",
    "posmode_groupB-groupA_heatmap.png"
  )))
  expect_true(file.exists(file.path(
    publication_graphs_dir,
    "Heatmaps",
    "negmode_groupC-groupA_heatmap.pdf"
  )))
  expect_true(file.exists(file.path(
    publication_graphs_dir,
    "Heatmaps",
    "all_heatmaps_combined.pdf"
  )))

  expect_false(outputMetabDaResultsAllContrasts(
    da_results_list = NULL,
    da_output_dir = da_output_dir,
    publication_graphs_dir = publication_graphs_dir
  ))
})
