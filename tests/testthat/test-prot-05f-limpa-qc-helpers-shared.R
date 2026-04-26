# fidelity-coverage-compare: shared
library(testthat)

localFakeLimpaQc <- function(env = parent.frame()) {
  old_lib_paths <- .libPaths()
  fake_lib <- file.path(tempdir(), "multischolar-fake-limpa-qc-lib")
  fake_pkg <- file.path(tempdir(), "multischolar-fake-limpa-qc-src")
  installed_pkg <- file.path(fake_lib, "limpa")

  if (!dir.exists(installed_pkg)) {
    unlink(fake_pkg, recursive = TRUE)
    dir.create(file.path(fake_pkg, "R"), recursive = TRUE, showWarnings = FALSE)
    dir.create(fake_lib, recursive = TRUE, showWarnings = FALSE)

    writeLines(
      c(
        "Package: limpa",
        "Version: 0.0.1",
        "Title: Fake limpa QC test double",
        "Description: Minimal namespace used by MultiScholaR QC coverage tests.",
        "License: MIT",
        "Encoding: UTF-8"
      ),
      file.path(fake_pkg, "DESCRIPTION")
    )
    writeLines("export(dpc,plotDPC)", file.path(fake_pkg, "NAMESPACE"))
    writeLines(
      c(
        "dpc <- function(y) {",
        "  list(dpc = c(0.25, 0.75), y = y)",
        "}",
        "",
        "plotDPC <- function(dpc_obj) {",
        "  ggplot2::ggplot(data.frame(x = c(0, 1), y = c(0, 1)), ggplot2::aes(x, y)) +",
        "    ggplot2::geom_line() +",
        "    ggplot2::ggtitle('Fake limpa DPC')",
        "}"
      ),
      file.path(fake_pkg, "R", "limpa.R")
    )
    utils::install.packages(fake_pkg, lib = fake_lib, repos = NULL, type = "source", quiet = TRUE)
  }

  .libPaths(c(fake_lib, old_lib_paths))
  withr::defer(.libPaths(old_lib_paths), envir = env)
  invisible(fake_lib)
}

localSharedProjectDirs <- function(value, .local_envir = parent.frame()) {
  had_project_dirs <- exists("project_dirs", envir = .GlobalEnv, inherits = FALSE)
  old_project_dirs <- if (had_project_dirs) get("project_dirs", envir = .GlobalEnv, inherits = FALSE) else NULL

  assign("project_dirs", value, envir = .GlobalEnv)

  withr::defer({
    if (had_project_dirs) {
      assign("project_dirs", old_project_dirs, envir = .GlobalEnv)
    } else if (exists("project_dirs", envir = .GlobalEnv, inherits = FALSE)) {
      rm("project_dirs", envir = .GlobalEnv)
    }
  }, envir = .local_envir)
}

makeSharedPeptideLimpaObject <- function(values, args = list(), is_logged_data = TRUE) {
  stopifnot(is.matrix(values))

  if (is.null(rownames(values))) {
    rownames(values) <- paste0("P", seq_len(nrow(values)), "%pep", seq_len(nrow(values)))
  }
  if (is.null(colnames(values))) {
    colnames(values) <- paste0("S", seq_len(ncol(values)))
  }

  peptide_rows <- do.call(rbind, lapply(seq_len(nrow(values)), function(i) {
    ids <- strsplit(rownames(values)[i], "%", fixed = TRUE)[[1]]
    data.frame(
      Protein.Ids = ids[[1]],
      Stripped.Sequence = ids[[2]],
      Run = colnames(values),
      Q.Value = rep(0.01, ncol(values)),
      Global.Q.Value = rep(0.01, ncol(values)),
      Proteotypic = rep(TRUE, ncol(values)),
      Precursor.Quantity = as.numeric(values[i, ]),
      Precursor.Normalised = as.numeric(values[i, ]),
      stringsAsFactors = FALSE
    )
  }))

  design_matrix <- data.frame(
    Run = colnames(values),
    group = rep(c("A", "B"), length.out = ncol(values)),
    replicate = paste0("r", seq_len(ncol(values))),
    stringsAsFactors = FALSE
  )

  methods::new(
    "PeptideQuantitativeData",
    peptide_data = peptide_rows,
    peptide_matrix = values,
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    q_value_column = "Q.Value",
    global_q_value_column = "Global.Q.Value",
    proteotypic_peptide_sequence_column = "Proteotypic",
    raw_quantity_column = "Precursor.Quantity",
    norm_quantity_column = "Precursor.Normalised",
    is_logged_data = is_logged_data,
    design_matrix = design_matrix,
    sample_id = "Run",
    group_id = "group",
    technical_replicate_id = "replicate",
    args = args
  )
}

makeSharedProteinLimpaQcObject <- function(values, args = list()) {
  stopifnot(is.matrix(values))

  if (is.null(rownames(values))) {
    rownames(values) <- paste0("P", seq_len(nrow(values)))
  }
  if (is.null(colnames(values))) {
    colnames(values) <- paste0("S", seq_len(ncol(values)))
  }

  protein_table <- data.frame(
    Protein.Ids = rownames(values),
    as.data.frame(values, check.names = FALSE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  design_matrix <- data.frame(
    Run = colnames(values),
    group = rep(c("A", "B"), length.out = ncol(values)),
    replicate = paste0("r", seq_len(ncol(values))),
    stringsAsFactors = FALSE
  )

  methods::new(
    "ProteinQuantitativeData",
    protein_quant_table = protein_table,
    protein_id_column = "Protein.Ids",
    design_matrix = design_matrix,
    protein_id_table = data.frame(
      Protein.Ids = rownames(values),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    technical_replicate_id = "replicate",
    args = args
  )
}

test_that("generateLimpaQCPlots requires stored limpa results", {
  after_object <- makeSharedProteinLimpaQcObject(matrix(
    c(1, 2,
      3, 4),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("P1", "P2"), c("S1", "S2"))
  ))

  expect_error(
    generateLimpaQCPlots(after_object, save_plots = FALSE, verbose = FALSE),
    "No limpa DPC results found"
  )
})

test_that("generateLimpaQCPlots covers peptide-level fallback and save-dir resolution", {
  before_object <- makeSharedPeptideLimpaObject(matrix(
    c(NA_real_, 4,
      2, NA_real_),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("P1%pep1", "P2%pep2"), c("S1", "S2"))
  ))
  after_object <- makeSharedPeptideLimpaObject(
    matrix(
      c(1, 4,
        2, 3),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("P1%pep1", "P2%pep2"), c("S1", "S2"))
    ),
    args = list(
      limpa_dpc_results = list(
        dpc_parameters = c(0.25, 0.75),
        slope_interpretation = "MNAR-like",
        missing_percentage_before = 50,
        dpc_method = "limpa_dpc"
      )
    )
  )

  plot_dir <- tempfile("shared-limpa-qc-")
  dir.create(plot_dir)
  localSharedProjectDirs(list(proteomics_exp1 = list(peptide_qc_dir = plot_dir)))

  plot_list <- generateLimpaQCPlots(
    after_object = after_object,
    before_object = before_object,
    save_plots = TRUE,
    save_dir = "peptide_qc",
    plot_prefix = "peptide_limpa",
    verbose = FALSE
  )

  expect_identical(
    names(plot_list),
    c("dpc_curve", "missing_comparison", "intensity_distribution", "summary")
  )
  expect_true(all(vapply(plot_list, inherits, logical(1), what = "ggplot")))
  expect_true(all(file.exists(file.path(plot_dir, c(
    "peptide_limpa_dpc_curve.png",
    "peptide_limpa_missing_comparison.png",
    "peptide_limpa_intensity_distribution.png",
    "peptide_limpa_summary.png",
    "peptide_limpa_composite.png"
  )))))
})

test_that("generateLimpaQCPlots covers explicit DPC plotting without saving", {
  localFakeLimpaQc()

  before_object <- makeSharedPeptideLimpaObject(matrix(
    c(NA_real_, 4,
      2, NA_real_),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("P1%pep1", "P2%pep2"), c("S1", "S2"))
  ))
  after_object <- makeSharedPeptideLimpaObject(
    matrix(
      c(1, 4,
        2, 3),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("P1%pep1", "P2%pep2"), c("S1", "S2"))
    ),
    args = list(
      limpa_dpc_results = list(
        dpc_parameters = c(0.25, 0.75),
        dpc_object = list(tag = "peptide"),
        slope_interpretation = "MNAR-like",
        missing_percentage_before = 50,
        dpc_method = "limpa_dpc"
      )
    )
  )

  plot_list <- generateLimpaQCPlots(
    after_object = after_object,
    before_object = before_object,
    save_plots = FALSE,
    save_dir = "peptide_qc",
    plot_prefix = "peptide_limpa",
    verbose = FALSE
  )

  expect_identical(
    names(plot_list),
    c("dpc_curve", "missing_comparison", "intensity_distribution", "summary")
  )
  expect_true(all(vapply(plot_list, inherits, logical(1), what = "ggplot")))
})

test_that("generateLimpaQCPlots covers protein-level comparison and no-save branch", {
  localFakeLimpaQc()

  before_object <- makeSharedPeptideLimpaObject(matrix(
    c(NA_real_, 4,
      2, NA_real_,
      3, 5),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("P1%pep1", "P2%pep2", "P3%pep3"), c("S1", "S2"))
  ))
  after_object <- makeSharedProteinLimpaQcObject(
    matrix(
      c(1, 4,
        2, 3),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("P1", "P2"), c("S1", "S2"))
    ),
    args = list(
      limpa_dpc_quant_results = list(
        dpc_parameters_used = c(0.1, 0.9),
        dpc_object_used = list(tag = "protein"),
        slope_interpretation = "MNAR",
        missing_percentage_before = 33.3,
        dpc_method = "limpa_dpc_quant"
      )
    )
  )

  plot_list <- generateLimpaQCPlots(
    after_object = after_object,
    before_object = before_object,
    save_plots = FALSE,
    save_dir = "not_used",
    plot_prefix = "protein_limpa",
    verbose = FALSE
  )

  expect_identical(
    names(plot_list),
    c("dpc_curve", "missing_comparison", "intensity_distribution", "summary")
  )
  expect_true(all(vapply(plot_list, inherits, logical(1), what = "ggplot")))
})
