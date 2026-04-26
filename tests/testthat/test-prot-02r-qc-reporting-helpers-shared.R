# fidelity-coverage-compare: shared
library(testthat)

captureSharedPrintedResult <- function(fn, ...) {
  output <- capture.output(value <- fn(...))
  list(value = value, output = output)
}

localSharedProteinFilteringProgress <- function(.local_envir = parent.frame()) {
  had_progress <- exists("filtering_progress", envir = .GlobalEnv, inherits = FALSE)
  old_progress <- if (had_progress) get("filtering_progress", envir = .GlobalEnv, inherits = FALSE) else NULL

  if (exists("filtering_progress", envir = .GlobalEnv, inherits = FALSE)) {
    rm("filtering_progress", envir = .GlobalEnv)
  }

  withr::defer({
    if (had_progress) {
      assign("filtering_progress", old_progress, envir = .GlobalEnv)
    } else if (exists("filtering_progress", envir = .GlobalEnv, inherits = FALSE)) {
      rm("filtering_progress", envir = .GlobalEnv)
    }
  }, envir = .local_envir)
}

localSharedGraphicsState <- function(.local_envir = parent.frame()) {
  old_device <- getOption("device")

  while (grDevices::dev.cur() > 1) {
    grDevices::dev.off()
  }

  options(device = grDevices::pdf)

  withr::defer({
    while (grDevices::dev.cur() > 1) {
      grDevices::dev.off()
    }
    options(device = old_device)
  }, envir = .local_envir)
}

withSharedPdfDevice <- function(code) {
  localSharedGraphicsState()
  pdf_path <- tempfile(fileext = ".pdf")
  dev_before <- grDevices::dev.list()
  grDevices::pdf(pdf_path)
  opened_dev <- grDevices::dev.cur()
  on.exit({
    active_devices <- grDevices::dev.list()
    if (is.null(active_devices)) {
      invisible(NULL)
    } else {
      active_ids <- as.integer(active_devices)
      if (opened_dev %in% active_ids && opened_dev != 1L) {
        grDevices::dev.off(which = opened_dev)
      } else {
        new_devices <- setdiff(active_ids, c(as.integer(dev_before), 1L))
        if (length(new_devices) > 0) {
          grDevices::dev.off(which = max(new_devices))
        }
      }
    }
  }, add = TRUE)
  force(code)
}

makeSharedProteinReportingObject <- function(values, groups = NULL) {
  stopifnot(is.matrix(values))

  if (is.null(rownames(values))) {
    rownames(values) <- paste0("P", seq_len(nrow(values)))
  }

  if (is.null(colnames(values))) {
    colnames(values) <- paste0("S", seq_len(ncol(values)))
  }

  if (is.null(groups)) {
    groups <- rep(c("A", "B"), length.out = ncol(values))
  }

  protein_table <- data.frame(
    Protein.Ids = rownames(values),
    as.data.frame(values, check.names = FALSE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  design_matrix <- data.frame(
    Run = colnames(values),
    group = groups,
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
    technical_replicate_id = "replicate"
  )
}

makeSharedPeptideFilteringData <- function() {
  data.frame(
    Protein.Ids = c("P1", "P1", "P1", "P2", "P2", "P3"),
    Stripped.Sequence = c("pep1", "pep1", "pep2", "pep3", "pep3", "pep4"),
    Run = c("S1", "S2", "S1", "S1", "S2", "S2"),
    Intensity = c(10, 11, 12, 13, 14, 15),
    stringsAsFactors = FALSE
  )
}

makeSharedProteinQuantFilteringData <- function() {
  data.frame(
    Protein.Ids = c("P1", "P2", "P3", "P4"),
    S1 = c(1, 4, 7, 10),
    S2 = c(2, 5, NA_real_, 11),
    S3 = c(3, NA_real_, 9, 12),
    check.names = FALSE
  )
}

test_that("protein NA helper rejects non-S4 inputs", {
  expect_error(
    checkProteinNAPercentages(list()),
    "ProteinQuantitativeData"
  )
})

test_that("protein NA recommendations cover all recommendation tiers", {
  cases <- list(
    list(
      object = makeSharedProteinReportingObject(matrix(
        c(1, 2, 3,
          4, 5, 6,
          7, 8, 9,
          10, 11, NA_real_),
        nrow = 4,
        byrow = TRUE,
        dimnames = list(NULL, c("S1", "S2", "S3"))
      )),
      recommendation = "complete_case",
      message = "Complete Case Analysis"
    ),
    list(
      object = makeSharedProteinReportingObject(matrix(
        c(1, 2, NA_real_,
          4, NA_real_, 6,
          7, 8, 9,
          10, 11, NA_real_),
        nrow = 4,
        byrow = TRUE,
        dimnames = list(NULL, c("S1", "S2", "S3"))
      )),
      recommendation = "imputation_or_filtering",
      message = "Consider Protein-Level Imputation"
    ),
    list(
      object = makeSharedProteinReportingObject(matrix(
        c(1, NA_real_, NA_real_,
          4, NA_real_, NA_real_,
          7, 8, NA_real_,
          10, 11, NA_real_),
        nrow = 4,
        byrow = TRUE,
        dimnames = list(NULL, c("S1", "S2", "S3"))
      )),
      recommendation = "strict_filtering",
      message = "Strict Filtering + Targeted Imputation"
    ),
    list(
      object = makeSharedProteinReportingObject(matrix(
        c(NA_real_, NA_real_, NA_real_,
          NA_real_, NA_real_, 6,
          NA_real_, NA_real_, NA_real_,
          10, NA_real_, NA_real_),
        nrow = 4,
        byrow = TRUE,
        dimnames = list(NULL, c("S1", "S2", "S3"))
      )),
      recommendation = "data_quality_review",
      message = "Review Data Quality"
    )
  )

  for (case in cases) {
    captured <- captureSharedPrintedResult(
      getProteinNARecommendations,
      case$object,
      include_code = FALSE
    )

    expect_identical(
      captured$value$primary_recommendation,
      case$recommendation
    )
    expect_true(any(grepl(case$message, captured$output, fixed = TRUE)))
  }
})

test_that("post-imputation protein validation reports expected warning tiers", {
  typical_object <- makeSharedProteinReportingObject(matrix(
    c(1, 2, NA_real_,
      4, NA_real_, 6,
      7, 8, 9,
      10, 11, NA_real_),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(NULL, c("S1", "S2", "S3"))
  ))
  typical_capture <- captureSharedPrintedResult(
    validatePostImputationProteinData,
    typical_object,
    expected_na_percent = NULL,
    tolerance = 10
  )

  expect_true(isTRUE(typical_capture$value$is_valid))
  expect_identical(typical_capture$value$expected_na_percent, 35)
  expect_true(any(grepl("NA percentage is typical for protein-level data", typical_capture$output, fixed = TRUE)))

  high_missing_object <- makeSharedProteinReportingObject(matrix(
    c(NA_real_, NA_real_, 3,
      NA_real_, NA_real_, 6,
      NA_real_, NA_real_, 9,
      NA_real_, 11, 12),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(NULL, c("S1", "S2", "S3"))
  ))
  high_missing_capture <- captureSharedPrintedResult(
    validatePostImputationProteinData,
    high_missing_object,
    expected_na_percent = NULL,
    tolerance = 10
  )

  expect_false(high_missing_capture$value$is_valid)
  expect_true(any(grepl("Very high NA percentage (>50%)", high_missing_capture$output, fixed = TRUE)))
  expect_true(any(grepl("Large variation in NA% between samples detected", high_missing_capture$output, fixed = TRUE)))
  expect_true(any(grepl("Samples with >80% missing proteins detected", high_missing_capture$output, fixed = TRUE)))

  low_missing_object <- makeSharedProteinReportingObject(matrix(
    c(1, 2, 3,
      4, 5, 6,
      7, 8, 9,
      10, 11, NA_real_),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(NULL, c("S1", "S2", "S3"))
  ))
  low_missing_capture <- captureSharedPrintedResult(
    validatePostImputationProteinData,
    low_missing_object,
    expected_na_percent = NULL,
    tolerance = 10
  )

  expect_false(low_missing_capture$value$is_valid)
  expect_true(any(grepl("Very low NA percentage (<10%) - excellent protein coverage!", low_missing_capture$output, fixed = TRUE)))
  expect_true(any(grepl("Possible over-imputation", low_missing_capture$output, fixed = TRUE)))
})

test_that("sample correlation helper excludes HEK runs and zeroes undefined correlations", {
  input_table <- data.frame(
    S1 = c(1, 1, 1),
    S2 = c(1, 2, 3),
    S3 = c(3, 4, 5),
    check.names = FALSE
  )
  metadata_tbl <- data.frame(
    Run = c("S1", "S2", "S3"),
    is_HEK = c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  correlation_matrix <- suppressWarnings(
    getSamplesCorrelationMatrix(input_table, metadata_tbl)
  )

  expect_identical(sort(colnames(correlation_matrix)), c("S1", "S2"))
  expect_identical(sort(rownames(correlation_matrix)), c("S1", "S2"))
  expect_identical(unname(correlation_matrix["S1", "S1"]), 0)
  expect_identical(unname(correlation_matrix["S1", "S2"]), 0)
  expect_identical(unname(correlation_matrix["S2", "S2"]), 1)
})

test_that("updateProteinFiltering records peptide metrics and saves shared plot artifacts", {
  localSharedProteinFilteringProgress()

  peptide_data <- makeSharedPeptideFilteringData()
  plot_dir <- tempfile("shared-prot-qc-reporting-")
  dir.create(plot_dir)

  saved_grid <- withSharedPdfDevice(
    updateProteinFiltering(
      peptide_data,
      step_name = "peptide_step",
      omic_type = "proteomics",
      experiment_label = "exp1",
      return_grid = TRUE,
      overwrite = FALSE,
      formats = "png",
      project_dirs = list(
        proteomics_exp1 = list(
          publication_graphs_dir = plot_dir,
          time_dir = plot_dir
        )
      )
    )
  )

  progress <- get("filtering_progress", envir = .GlobalEnv)

  expect_false(is.null(saved_grid))
  expect_s4_class(progress, "FilteringProgress")
  expect_identical(progress@steps, "peptide_step")
  expect_identical(unname(progress@proteins), 3)
  expect_identical(unname(progress@total_peptides), 4)
  expect_identical(nrow(progress@proteins_per_run[[1]]), 2L)
  expect_identical(nrow(progress@peptides_per_protein[[1]]), 3L)
  expect_true(all(file.exists(file.path(plot_dir, c(
    "peptide_step_proteins_total.png",
    "peptide_step_proteins_per_run.png",
    "peptide_step_peptides_total.png",
    "peptide_step_peptides_per_protein.png",
    "peptide_step_peptides_per_run.png",
    "peptide_step_combined_plots.png"
  )))))
})

test_that("updateProteinFiltering carries peptide metrics into protein steps and enforces overwrite", {
  localSharedProteinFilteringProgress()

  peptide_data <- makeSharedPeptideFilteringData()
  protein_data <- makeSharedProteinQuantFilteringData()

  withSharedPdfDevice(
    updateProteinFiltering(
      peptide_data,
      step_name = "peptide_step",
      return_grid = TRUE
    )
  )
  withSharedPdfDevice(
    suppressWarnings(
      updateProteinFiltering(
        protein_data,
        step_name = "protein_step",
        return_grid = TRUE
      )
    )
  )

  progress <- get("filtering_progress", envir = .GlobalEnv)

  expect_identical(progress@steps, c("peptide_step", "protein_step"))
  expect_identical(length(progress@total_peptides), 2L)
  expect_identical(progress@total_peptides[[2]], progress@total_peptides[[1]])
  expect_identical(unname(progress@proteins), c(3, 4))

  expect_error(
    suppressWarnings(
      updateProteinFiltering(
        protein_data,
        step_name = "protein_step",
        return_grid = FALSE
      )
    ),
    "already exists"
  )

  withSharedPdfDevice(
    suppressWarnings(
      updateProteinFiltering(
        protein_data,
        step_name = "protein_step",
        overwrite = TRUE,
        return_grid = TRUE
      )
    )
  )

  overwritten_progress <- get("filtering_progress", envir = .GlobalEnv)
  expect_identical(overwritten_progress@steps, c("peptide_step", "protein_step"))
  expect_identical(unname(overwritten_progress@proteins), c(3, 4))
})

test_that("updateProteinFiltering initializes protein-only progress with placeholder peptide metrics", {
  localSharedProteinFilteringProgress()

  protein_object <- makeSharedProteinReportingObject(matrix(
    c(1, 2, 3,
      4, 5, 6,
      7, NA_real_, 9,
      10, 11, 12),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(NULL, c("S1", "S2", "S3"))
  ))

  saved_grid <- withSharedPdfDevice(
    suppressWarnings(
      updateProteinFiltering(
        protein_object,
        step_name = "protein_only",
        return_grid = TRUE
      )
    )
  )

  progress <- get("filtering_progress", envir = .GlobalEnv)

  expect_false(is.null(saved_grid))
  expect_identical(progress@steps, "protein_only")
  expect_true(is.na(progress@total_peptides))
  expect_identical(nrow(progress@peptides_per_protein[[1]]), 0L)
  expect_identical(nrow(progress@peptides_per_run[[1]]), 0L)
})
