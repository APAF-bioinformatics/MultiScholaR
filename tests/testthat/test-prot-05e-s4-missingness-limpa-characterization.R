# fidelity-coverage-compare: shared
library(testthat)

localFakeLimpa <- function(env = parent.frame()) {
  old_lib_paths <- .libPaths()
  old_options <- options(
    multischolar.fake_limpa.dpc = NULL,
    multischolar.fake_limpa.dpcImpute = NULL
  )

  fake_lib <- file.path(tempdir(), "multischolar-fake-limpa-lib")
  fake_pkg <- file.path(tempdir(), "multischolar-fake-limpa-src")
  installed_pkg <- file.path(fake_lib, "limpa")

  if (!dir.exists(installed_pkg)) {
    unlink(fake_pkg, recursive = TRUE)
    dir.create(file.path(fake_pkg, "R"), recursive = TRUE, showWarnings = FALSE)
    dir.create(fake_lib, recursive = TRUE, showWarnings = FALSE)

    writeLines(
      c(
        "Package: limpa",
        "Version: 0.0.1",
        "Title: Fake limpa test double",
        "Description: Minimal namespace used by MultiScholaR characterization tests.",
        "License: MIT",
        "Encoding: UTF-8"
      ),
      file.path(fake_pkg, "DESCRIPTION")
    )
    writeLines("export(dpc,dpcImpute)", file.path(fake_pkg, "NAMESPACE"))
    writeLines(
      c(
        ".call_handler <- function(handler, args) {",
        "  formal_names <- names(formals(handler))",
        "  if (is.null(formal_names) || '...' %in% formal_names) {",
        "    return(do.call(handler, args))",
        "  }",
        "  do.call(handler, args[intersect(names(args), formal_names)])",
        "}",
        "",
        "dpc <- function(y) {",
        "  handler <- getOption('multischolar.fake_limpa.dpc')",
        "  if (is.function(handler)) {",
        "    return(.call_handler(handler, list(y = y)))",
        "  }",
        "  list(dpc = c(0.25, 0.75), y = y)",
        "}",
        "",
        "dpcImpute <- function(y, dpc, dpc.slope, verbose = TRUE, chunk = 1000) {",
        "  handler <- getOption('multischolar.fake_limpa.dpcImpute')",
        "  if (is.function(handler)) {",
        "    dpc_arg <- if (missing(dpc)) NULL else dpc",
        "    slope_arg <- if (missing(dpc.slope)) NULL else dpc.slope",
        "    return(.call_handler(handler, list(y = y, dpc = dpc_arg, dpc.slope = slope_arg, verbose = verbose, chunk = chunk)))",
        "  }",
        "  E <- y",
        "  replacement <- if (missing(dpc.slope)) 1 else dpc.slope",
        "  E[is.na(E)] <- replacement",
        "  list(E = E)",
        "}"
      ),
      file.path(fake_pkg, "R", "limpa.R")
    )
    utils::install.packages(fake_pkg, lib = fake_lib, repos = NULL, type = "source", quiet = TRUE)
  }

  .libPaths(c(fake_lib, old_lib_paths))

  withr::defer({
    options(old_options)
    .libPaths(old_lib_paths)
  }, envir = env)

  invisible(fake_lib)
}

newProteinLimpaObject <- function(values,
                                  row_samples = colnames(values),
                                  args = list()) {
  stopifnot(length(row_samples) == nrow(values))
  protein_table <- data.frame(
    Protein.Ids = paste0("P", seq_len(nrow(values))),
    description = row_samples,
    values,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  new(
    "ProteinQuantitativeData",
    protein_quant_table = protein_table,
    protein_id_column = "Protein.Ids",
    design_matrix = data.frame(
      description = c("description", colnames(values)),
      group = "G1",
      replicates = "R1",
      stringsAsFactors = FALSE
    ),
    sample_id = "description",
    group_id = "group",
    technical_replicate_id = "replicates",
    args = args
  )
}

test_that("ProteinQuantitativeData limpa imputation preserves provided DPC behavior on raw-scale data", {
  localFakeLimpa()
  captured <- new.env(parent = emptyenv())
  options(
    multischolar.fake_limpa.dpcImpute = function(y, dpc, dpc.slope, verbose, chunk) {
      captured$y <- y
      captured$dpc <- dpc
      captured$dpc.slope <- dpc.slope
      captured$verbose <- verbose
      captured$chunk <- chunk
      E <- y
      E[is.na(E)] <- log2(9)
      list(E = E)
    }
  )

  protein_object <- newProteinLimpaObject(
    values = data.frame(
      S1 = c(NA_real_, 100),
      S2 = c(200, NA_real_),
      check.names = FALSE
    ),
    row_samples = c("S1", "S2")
  )

  imputed <- proteinMissingValueImputationLimpa(
    protein_object,
    dpc_results = list(dpc = c(0.1, 0.9), source = "provided"),
    quantified_protein_column = "Protein.Imputed.Test",
    verbose = TRUE,
    chunk = 7
  )

  expect_equal(captured$y["P1", "S2"], log2(201))
  expect_equal(captured$dpc$dpc, c(0.1, 0.9))
  expect_null(captured$dpc.slope)
  expect_true(captured$verbose)
  expect_identical(captured$chunk, 7)
  expect_equal(imputed@protein_quant_table$Protein.Imputed.Test, c(8, 8))
  expect_equal(
    imputed@args$limpa_protein_imputation_results$dpc_parameters_used,
    c(0.1, 0.9)
  )
  expect_identical(
    imputed@args$limpa_protein_imputation_results$imputation_method,
    "limpa_dpc_protein_imputation"
  )
  expect_identical(
    imputed@args$limpa_protein_imputation_results$total_proteins_imputed,
    2L
  )
})

test_that("ProteinQuantitativeData limpa imputation estimates DPC when no result is supplied", {
  localFakeLimpa()
  captured <- new.env(parent = emptyenv())
  options(
    multischolar.fake_limpa.dpc = function(y) {
      captured$dpc_input <- y
      list(dpc = c(0.2, 0.8), fitted = TRUE)
    },
    multischolar.fake_limpa.dpcImpute = function(y, dpc, dpc.slope, verbose, chunk) {
      captured$impute_dpc <- dpc
      captured$impute_slope <- dpc.slope
      E <- y
      E[is.na(E)] <- 6
      list(E = E)
    }
  )

  protein_object <- newProteinLimpaObject(
    values = data.frame(
      S1 = c(10, NA_real_),
      S2 = c(11, 12),
      check.names = FALSE
    ),
    row_samples = c("S1", "S2")
  )

  imputed <- proteinMissingValueImputationLimpa(protein_object, verbose = TRUE)

  expect_equal(captured$dpc_input["P1", "S1"], 10)
  expect_equal(captured$impute_dpc$dpc, c(0.2, 0.8))
  expect_null(captured$impute_slope)
  expect_identical(names(imputed@protein_quant_table)[5], "Protein.Imputed.Limpa")
  expect_equal(imputed@protein_quant_table$Protein.Imputed.Limpa, c(10, 12))
  expect_equal(
    imputed@args$limpa_protein_imputation_results$dpc_parameters_used,
    c(0.2, 0.8)
  )
  expect_true(is.list(imputed@args$limpa_protein_imputation_results$dpc_object_used))
})

test_that("ProteinQuantitativeData limpa imputation falls back to slope and cleans nonfinite values", {
  localFakeLimpa()
  captured <- new.env(parent = emptyenv())
  options(
    multischolar.fake_limpa.dpc = function(y) {
      captured$dpc_input <- y
      stop("mock dpc failure")
    },
    multischolar.fake_limpa.dpcImpute = function(y, dpc, dpc.slope, verbose, chunk) {
      captured$impute_input <- y
      captured$dpc <- dpc
      captured$dpc.slope <- dpc.slope
      E <- y
      E[is.na(E)] <- dpc.slope
      list(E = E)
    }
  )

  protein_object <- newProteinLimpaObject(
    values = data.frame(
      S1 = c(NaN, 15),
      S2 = c(14, NA_real_),
      check.names = FALSE
    ),
    row_samples = c("S1", "S2"),
    args = list(existing = "kept")
  )

  imputed <- proteinMissingValueImputationLimpa(
    protein_object,
    dpc_slope = 0.6,
    quantified_protein_column = "Fallback.Imputed",
    verbose = TRUE
  )

  expect_true(is.na(captured$dpc_input["P1", "S1"]))
  expect_true(is.na(captured$impute_input["P1", "S1"]))
  expect_null(captured$dpc)
  expect_identical(captured$dpc.slope, 0.6)
  expect_equal(imputed@protein_quant_table$Fallback.Imputed, c(0.6, 0.6))
  expect_identical(imputed@args$existing, "kept")
  expect_equal(
    imputed@args$limpa_protein_imputation_results$dpc_parameters_used,
    c(NA, 0.6)
  )
  expect_null(imputed@args$limpa_protein_imputation_results$dpc_object_used)
})

test_that("ProteinQuantitativeData limpa imputation records numeric DPC parameters", {
  localFakeLimpa()
  captured <- new.env(parent = emptyenv())
  options(
    multischolar.fake_limpa.dpcImpute = function(y, dpc, dpc.slope, verbose, chunk) {
      captured$dpc <- dpc
      E <- y
      E[is.na(E)] <- 3
      list(E = E)
    }
  )

  protein_object <- newProteinLimpaObject(
    values = data.frame(
      S1 = c(1, 2),
      S2 = c(NA_real_, 4),
      check.names = FALSE
    ),
    row_samples = c("S1", "S2")
  )

  imputed <- proteinMissingValueImputationLimpa(
    protein_object,
    dpc_results = c(0.3, 0.7),
    verbose = FALSE
  )

  expect_equal(captured$dpc, c(0.3, 0.7))
  expect_equal(imputed@protein_quant_table$Protein.Imputed.Limpa, c(1, 4))
  expect_equal(
    imputed@args$limpa_protein_imputation_results$dpc_parameters_used,
    c(0.3, 0.7)
  )
  expect_null(imputed@args$limpa_protein_imputation_results$dpc_object_used)
})

test_that("ProteinQuantitativeData limpa imputation wraps dpcImpute failures", {
  localFakeLimpa()
  options(
    multischolar.fake_limpa.dpcImpute = function(y, dpc, dpc.slope, verbose, chunk) {
      stop("mock dpcImpute failure")
    }
  )

  protein_object <- newProteinLimpaObject(
    values = data.frame(
      S1 = c(1, NA_real_),
      S2 = c(2, 3),
      check.names = FALSE
    ),
    row_samples = c("S1", "S2")
  )

  expect_error(
    proteinMissingValueImputationLimpa(
      protein_object,
      dpc_results = list(dpc = c(0.2, 0.8)),
      verbose = FALSE
    ),
    "limpa protein imputation failed: mock dpcImpute failure",
    fixed = TRUE
  )
})
