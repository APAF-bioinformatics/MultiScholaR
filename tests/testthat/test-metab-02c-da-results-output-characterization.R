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

loadSelectedFunctions(
  paths = c(
    file.path(repo_root, "R", "func_metab_da_results_output.R"),
    file.path(repo_root, "R", "func_metab_da_volcano_static.R"),
    file.path(repo_root, "R", "func_metab_da.R")
  ),
  symbols = c(
    "outputMetabDaResultsAllContrasts",
    "generateMetabDAVolcanoStatic"
  ),
  env = environment()
)

generateMetabDAHeatmap <- function(...) {
  NULL
}

buildMetabDaOutputResults <- function() {
  data.frame(
    metabolite_id = c("M1", "M2"),
    metabolite_name = c("Met One", "Met Two"),
    logFC = c(1.8, -0.7),
    raw_pvalue = c(0.001, 0.20),
    fdr_qvalue = c(0.01, 0.20),
    significant = c("Up", "NS"),
    comparison = c("groupB-groupA", "groupB-groupA"),
    friendly_name = c("groupB-groupA = B vs A", "groupB-groupA = B vs A"),
    numerator = c("groupB", "groupB"),
    denominator = c("groupA", "groupA"),
    assay = c("LCMS_Pos", "LCMS_Pos"),
    intensity.S1.groupA = c(10, 12),
    intensity.S2.groupB = c(20, 25),
    stringsAsFactors = FALSE
  )
}

test_that("metabolomics DA results output helper preserves tabular and volcano exports", {
  tmp_root <- file.path(tempdir(), paste0("metab-da-output-", Sys.getpid()))
  unlink(tmp_root, recursive = TRUE, force = TRUE)
  dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

  da_output_dir <- file.path(tmp_root, "de-results")
  publication_graphs_dir <- file.path(tmp_root, "publication")

  expect_true(outputMetabDaResultsAllContrasts(
    da_results_list = list(da_metabolites_long = buildMetabDaOutputResults()),
    da_output_dir = da_output_dir,
    publication_graphs_dir = publication_graphs_dir,
    da_q_val_thresh = 0.05,
    lfc_threshold = 1
  ))

  tsv_path <- file.path(
    da_output_dir,
    "de_posmode_metabolites_groupB-groupA_long_annot.tsv"
  )
  expect_true(file.exists(tsv_path))

  exported <- read.delim(tsv_path, check.names = FALSE)
  expect_equal(
    colnames(exported)[1:6],
    c(
      "metabolite_id",
      "metabolite_name",
      "logFC",
      "raw_pvalue",
      "fdr_qvalue",
      "significant"
    )
  )
  expect_equal(exported$metabolite_id, c("M1", "M2"))
  expect_equal(exported$intensity.S2.groupB, c(20, 25))

  expect_true(file.exists(file.path(
    da_output_dir,
    "de_posmode_metabolites_groupB-groupA_long_annot.xlsx"
  )))
  expect_true(file.exists(file.path(
    publication_graphs_dir,
    "Volcano_Plots",
    "posmode_groupB-groupA.png"
  )))
  expect_true(file.exists(file.path(
    publication_graphs_dir,
    "Volcano_Plots",
    "posmode_groupB-groupA.pdf"
  )))
  expect_true(file.exists(file.path(
    publication_graphs_dir,
    "Volcano_Plots",
    "all_volcano_plots_combined.pdf"
  )))
  expect_true(dir.exists(file.path(publication_graphs_dir, "Heatmaps")))
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
  expect_equal(numsig$significant, 1)
  expect_equal(numsig$mode, "posmode")
  expect_true(file.exists(file.path(
    publication_graphs_dir,
    "NumSigDeMolecules",
    "metabolites_num_sig_de_molecules.xlsx"
  )))
  expect_true(file.exists(file.path(
    publication_graphs_dir,
    "NumSigDeMolecules",
    "metabolites_num_sig_barplot.png"
  )))
})

test_that("metabolomics DA results output helper preserves null-result exit", {
  expect_false(outputMetabDaResultsAllContrasts(
    da_results_list = NULL,
    da_output_dir = tempdir(),
    publication_graphs_dir = tempdir()
  ))
})
