# fidelity-coverage-compare: shared
# testthat for Proteomics DA Analysis
# Phase 4 of Proteomics GUI Test Strategy

repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)

skipIfMissingProtDaSplitTargetFiles <- function() {
  required_paths <- c(
    "R/func_prot_da_model_methods.R",
    "R/func_prot_da_results_methods.R",
    "R/func_prot_da_heatmap.R"
  )
  missing <- required_paths[!file.exists(file.path(repo_root, required_paths))]
  if (length(missing) > 0) {
    testthat::skip(
      sprintf(
        "Target-only prot DA split file(s) not present: %s",
        paste(basename(missing), collapse = ", ")
      )
    )
  }
}

skipIfMissingProtDaSplitTargetFiles()

localSharedDaGraphicsState <- function(.local_envir = parent.frame()) {
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

localSharedDaMethodMocks <- function(.local_envir = parent.frame()) {
  testthat::local_mocked_bindings(
    plotPca = function(...) {
      ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
        ggplot2::geom_point()
    },
    plotNumOfValuesNoLog = function(...) {
      ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
        ggplot2::geom_col()
    },
    plotVolcano = function(...) {
      ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
        ggplot2::geom_point()
    },
    printPValuesDistribution = function(...) {
      ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
        ggplot2::geom_point()
    },
    printCountDaGenesTable = function(...) {
      data.frame(
        comparison = "G2_vs_G1",
        status = "Significant and Up",
        counts = 1L,
        stringsAsFactors = FALSE
      )
    },
    plotOneVolcanoNoVerticalLines = function(...) {
      ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
        ggplot2::geom_point()
    },
    .package = "MultiScholaR",
    .env = .local_envir
  )
}

test_that("DA analysis snapshot is valid", {
  cp_file <- test_path("..", "testdata", "sepsis", "proteomics", "cp07_da_results.rds")
  
  if (file.exists(cp_file)) {
    first_line <- tryCatch(readLines(cp_file, n = 1, warn = FALSE), error = function(e) "")
    if (length(first_line) > 0 && identical(first_line[[1]], "version https://git-lfs.github.com/spec/v1")) {
      skip("Snapshot cp07 is a Git LFS pointer and the binary artifact is not present")
    }

    results <- readRDS(cp_file)
    expect_type(results, "list")
    expect_true("da_proteins_long" %in% names(results))
    expect_true(nrow(results$da_proteins_long) > 0)
  } else {
    skip("Snapshot cp07 not found")
  }
})

test_that("differentialAbundanceAnalysis works with mock data", {
  localSharedDaGraphicsState()
  localSharedDaMethodMocks()

  # Mock ProteinQuantitativeData
  n_prot <- 20
  pqd <- new("ProteinQuantitativeData",
    protein_quant_table = data.frame(
      Protein.Ids = paste0("P", 1:n_prot),
      # Use deterministic group-shifted values so limma stays stable under the
      # larger shared-suite load order.
      S1 = seq(10, 10 + n_prot - 1),
      S2 = seq(10.5, 10.5 + n_prot - 1),
      S3 = seq(15, 15 + n_prot - 1),
      S4 = seq(15.5, 15.5 + n_prot - 1),
      stringsAsFactors = FALSE
    ),
    design_matrix = data.frame(
      Run = c("S1", "S2", "S3", "S4"),
      Group = c("G1", "G1", "G2", "G2"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "Group",
    protein_id_column = "Protein.Ids",
    protein_id_table = data.frame(
      Protein.Ids = paste0("P", 1:n_prot),
      Gene.Names = paste0("GENE", 1:n_prot),
      stringsAsFactors = FALSE
    ),
    args = list()
  )
  
  # Mock contrasts table
  contrasts <- data.frame(
    full_format = "G2_vs_G1=GroupG2-GroupG1",
    contrasts = "GroupG2-GroupG1",
    friendly_names = "G2_vs_G1",
    stringsAsFactors = FALSE
  )
  
  # Run analysis
  # Note: internal functions might expect specific column names in design matrix
  result <- differentialAbundanceAnalysis(
    theObject = pqd,
    contrasts_tbl = contrasts,
    formula_string = "~ 0 + Group",
    args_row_id = "Protein.Ids",
    group_id = "Group"
  )
  
  expect_type(result, "list")
  expect_true("da_proteins_long" %in% names(result))
  expect_true("G2_vs_G1" %in% result$da_proteins_long$comparison)
})

# APAF Bioinformatics | test-prot-07-da-analysis.R | Approved | 2026-03-13
