# fidelity-coverage-compare: shared
library(testthat)

PeptideQuantitativeData <- get(
  "PeptideQuantitativeData",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)

ProteinQuantitativeData <- get(
  "ProteinQuantitativeData",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)

preservePeptideNaValuesHelper <- get(
  "preservePeptideNaValuesHelper",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)

makeFunctionWithOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

localFakeLimpa <- function(env = parent.frame()) {
  old_lib_paths <- .libPaths()
  old_options <- options(
    multischolar.fake_limpa.dpc = NULL,
    multischolar.fake_limpa.dpcImpute = NULL,
    multischolar.fake_limpa.dpcQuant = NULL
  )

  fake_lib <- file.path(tempdir(), "multischolar-fake-limpa-prot-s4-lib")
  fake_pkg <- file.path(tempdir(), "multischolar-fake-limpa-prot-s4-src")
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
        "Description: Minimal namespace used by MultiScholaR missingness coverage tests.",
        "License: MIT",
        "Encoding: UTF-8"
      ),
      file.path(fake_pkg, "DESCRIPTION")
    )
    writeLines("export(dpc,dpcImpute,dpcQuant)", file.path(fake_pkg, "NAMESPACE"))
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
        "  E[is.na(E)] <- if (missing(dpc.slope)) 1 else dpc.slope",
        "  list(E = E)",
        "}",
        "",
        "dpcQuant <- function(y, protein.id, dpc, dpc.slope, verbose = TRUE, chunk = 1000) {",
        "  handler <- getOption('multischolar.fake_limpa.dpcQuant')",
        "  if (is.function(handler)) {",
        "    dpc_arg <- if (missing(dpc)) NULL else dpc",
        "    slope_arg <- if (missing(dpc.slope)) NULL else dpc.slope",
        "    return(.call_handler(handler, list(y = y, protein.id = protein.id, dpc = dpc_arg, dpc.slope = slope_arg, verbose = verbose, chunk = chunk)))",
        "  }",
        "  protein_ids <- unique(y$genes[[protein.id]])",
        "  E <- vapply(",
        "    protein_ids,",
        "    function(pid) {",
        "      keep <- y$genes[[protein.id]] == pid",
        "      values <- colMeans(y$E[keep, , drop = FALSE], na.rm = TRUE)",
        "      values[is.nan(values)] <- 0",
        "      values",
        "    },",
        "    numeric(ncol(y$E))",
        "  )",
        "  E <- t(E)",
        "  rownames(E) <- protein_ids",
        "  colnames(E) <- colnames(y$E)",
        "  other <- list(",
        "    standard.error = matrix(0.1, nrow(E), ncol(E), dimnames = dimnames(E)),",
        "    n.observations = matrix(1L, nrow(E), ncol(E), dimnames = dimnames(E))",
        "  )",
        "  list(E = E, genes = data.frame(protein.id = protein_ids, stringsAsFactors = FALSE), other = other)",
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

newProteinMissingnessObject <- function(values,
                                        row_samples = colnames(values),
                                        args = list()) {
  ProteinQuantitativeData(
    protein_quant_table = data.frame(
      Protein.Ids = paste0("P", seq_len(nrow(values))),
      description = row_samples,
      values,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
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

newPeptideMissingnessObject <- function(values = matrix(
                                          c(NA_real_, 120, 240, NA_real_),
                                          nrow = 2,
                                          byrow = TRUE,
                                          dimnames = list(c("P1%PEP1", "P2%PEP2"), c("S1", "S2"))
                                        ),
                                        calculate_matrix = FALSE,
                                        is_logged_data = FALSE,
                                        args = list()) {
  peptide_rows <- expand.grid(
    peptide_index = seq_len(nrow(values)),
    Run = colnames(values),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  peptide_rows$Protein.Ids <- sub("%.*$", "", rownames(values)[peptide_rows$peptide_index])
  peptide_rows$Stripped.Sequence <- sub("^.*%", "", rownames(values)[peptide_rows$peptide_index])
  peptide_rows$Precursor.Quantity <- 1000 + peptide_rows$peptide_index * 10
  peptide_rows$Precursor.Normalised <- as.numeric(values[cbind(
    peptide_rows$peptide_index,
    match(peptide_rows$Run, colnames(values))
  )])
  peptide_rows$Q.Value <- 0.001
  peptide_rows$Global.Q.Value <- 0.001
  peptide_rows$Proteotypic <- TRUE
  peptide_rows <- peptide_rows[
    c(
      "Protein.Ids",
      "Stripped.Sequence",
      "Run",
      "Precursor.Quantity",
      "Precursor.Normalised",
      "Q.Value",
      "Global.Q.Value",
      "Proteotypic"
    )
  ]

  peptide_object <- PeptideQuantitativeData(
    peptide_data = peptide_rows,
    peptide_matrix = values,
    design_matrix = data.frame(
      Run = colnames(values),
      group = rep(c("G1", "G2"), length.out = ncol(values)),
      replicates = rep(c("R1", "R2"), length.out = ncol(values)),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    technical_replicate_id = "replicates",
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    q_value_column = "Q.Value",
    global_q_value_column = "Global.Q.Value",
    proteotypic_peptide_sequence_column = "Proteotypic",
    raw_quantity_column = "Precursor.Quantity",
    norm_quantity_column = "Precursor.Normalised",
    is_logged_data = is_logged_data,
    args = args
  )

  if (calculate_matrix) {
    peptide_object <- calcPeptideMatrix(peptide_object)
  }

  peptide_object
}

test_that("protein missingness methods cover direct package routes and DPC fallbacks", {
  localFakeLimpa()

  quiet_replacements <- list(
    log_info = function(...) invisible(NULL),
    log_warn = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    checkParamsObjectFunctionSimplify = function(theObject, name, default) {
      value <- theObject@args[[name]]
      if (is.null(value)) default else value
    }
  )

  protein_method <- makeFunctionWithOverrides(
    methods::selectMethod("proteinMissingValueImputationLimpa", "ProteinQuantitativeData"),
    quiet_replacements
  )
  peptide_method <- makeFunctionWithOverrides(
    methods::selectMethod("proteinMissingValueImputationLimpa", "PeptideQuantitativeData"),
    quiet_replacements
  )

  protein_captured <- new.env(parent = emptyenv())
  peptide_captured <- new.env(parent = emptyenv())

  options(
    multischolar.fake_limpa.dpc = function(y) {
      protein_captured$dpc_input <- y
      list(dpc = c(0.2, 0.6), y = y)
    },
    multischolar.fake_limpa.dpcImpute = function(y, dpc = NULL, dpc.slope = NULL, verbose, chunk) {
      protein_captured$impute_input <- y
      protein_captured$dpc <- dpc
      protein_captured$dpc_slope <- dpc.slope
      E <- y
      E[is.na(E)] <- log2(9)
      list(E = E)
    },
    multischolar.fake_limpa.dpcQuant = function(y, protein.id, dpc = NULL, dpc.slope = NULL, verbose, chunk) {
      peptide_captured$y <- y
      peptide_captured$protein.id <- protein.id
      peptide_captured$dpc <- dpc
      peptide_captured$dpc_slope <- dpc.slope

      protein_ids <- unique(y$genes[[protein.id]])
      E <- matrix(
        c(4, 5, 6, 7),
        nrow = length(protein_ids),
        byrow = TRUE,
        dimnames = list(protein_ids, colnames(y$E))
      )
      other <- list(
        standard.error = matrix(0.1, nrow(E), ncol(E), dimnames = dimnames(E)),
        n.observations = matrix(1L, nrow(E), ncol(E), dimnames = dimnames(E))
      )

      list(
        E = E,
        genes = data.frame(protein.id = protein_ids, stringsAsFactors = FALSE),
        other = other
      )
    }
  )

  protein_object <- newProteinMissingnessObject(
    values = data.frame(
      S1 = c(NA_real_, 120),
      S2 = c(Inf, NaN),
      check.names = FALSE
    ),
    row_samples = c("S1", "S2")
  )

  protein_imputed <- protein_method(
    protein_object,
    dpc_results = NULL,
    quantified_protein_column = "Protein.Imputed.Direct",
    verbose = FALSE
  )

  expect_s4_class(protein_imputed, "ProteinQuantitativeData")
  expect_equal(protein_captured$dpc_input["P2", "S1"], log2(121))
  expect_true(all(is.na(protein_captured$impute_input[, "S2"])))
  expect_equal(protein_imputed@protein_quant_table$Protein.Imputed.Direct, c(8, 8))
  expect_equal(
    protein_imputed@args$limpa_protein_imputation_results$dpc_parameters_used,
    c(0.2, 0.6)
  )

  protein_supplied <- protein_method(
    newProteinMissingnessObject(
      values = data.frame(
        S1 = c(50, NA_real_),
        S2 = c(75, 25),
        check.names = FALSE
      ),
      row_samples = c("S1", "S2")
    ),
    dpc_results = c(0.3, 0.9),
    quantified_protein_column = "Protein.Imputed.Supplied",
    verbose = FALSE
  )
  expect_null(protein_captured$dpc_slope)
  expect_equal(
    protein_supplied@args$limpa_protein_imputation_results$dpc_parameters_used,
    c(0.3, 0.9)
  )

  options(multischolar.fake_limpa.dpc = function(y) stop("mock dpc failure"))

  peptide_object <- newPeptideMissingnessObject(
    values = matrix(
      c(
        NA_real_, 120,
        NA_real_, NA_real_,
        240, 360
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("P1%PEP1", "P2%PEP2", "P1%PEP3"), c("S1", "S2"))
    ),
    is_logged_data = FALSE,
    args = list(
      dpc_slope = 0.8,
      verbose = FALSE,
      chunk = 1000
    )
  )

  peptide_imputed <- peptide_method(
    peptide_object,
    dpc_results = NULL,
    quantified_protein_column = "Protein.Quantified.Direct",
    verbose = FALSE
  )

  expect_s4_class(peptide_imputed, "ProteinQuantitativeData")
  expect_null(peptide_captured$dpc)
  expect_identical(peptide_captured$dpc_slope, 0.8)
  expect_equal(nrow(peptide_captured$y$E), 2)
  expect_identical(peptide_captured$protein.id, "protein.id")
  expect_equal(peptide_imputed@args$limpa_dpc_quant_results$total_peptides_used, 2)
  expect_equal(peptide_imputed@args$limpa_dpc_quant_results$total_proteins_quantified, 2)
  expect_identical(
    peptide_imputed@args$limpa_dpc_quant_results$slope_interpretation,
    "unable to determine (DPC estimation failed)"
  )
})

test_that("preservePeptideNaValuesHelper preserves direct NA synchronization and validation", {
  peptide_object <- PeptideQuantitativeData(
    peptide_data = data.frame(
      Protein.Ids = c("P1", "P1", "P1", "P1", "P2", "P2"),
      Stripped.Sequence = c("pep1", "pep2", "pep1", "pep2", "pep3", "pep3"),
      Run = c("S1", "S1", "S2", "S2", "S1", "S2"),
      Peptide.Normalised = c(NA_real_, NA_real_, 4, 5, 6, NA_real_),
      Q.Value = rep(0.001, 6),
      Global.Q.Value = rep(0.001, 6),
      stringsAsFactors = FALSE
    ),
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    q_value_column = "Q.Value",
    global_q_value_column = "Global.Q.Value",
    raw_quantity_column = "Peptide.Normalised",
    norm_quantity_column = "Peptide.Normalised",
    design_matrix = data.frame(
      Run = c("S1", "S2"),
      group = c("A", "A"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group"
  )

  protein_object <- ProteinQuantitativeData(
    protein_quant_table = data.frame(
      S1 = c(10, 30),
      S2 = c(20, 40),
      row.names = c("1", "2"),
      check.names = FALSE
    ),
    design_matrix = data.frame(
      Run = c("S1", "S2"),
      group = c("A", "A"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    protein_id_column = "Protein.Ids"
  )

  preserved <- preservePeptideNaValuesHelper(peptide_object, protein_object)
  expect_true(is.na(preserved@protein_quant_table[1, "S1"]))
  expect_identical(preserved@protein_quant_table[1, "S2"], 20)
  expect_identical(preserved@protein_quant_table[2, "S1"], 30)
  expect_true(is.na(preserved@protein_quant_table[2, "S2"]))

  bad_rows <- protein_object
  rownames(bad_rows@protein_quant_table) <- c("1", "3")
  expect_error(
    preservePeptideNaValuesHelper(peptide_object, bad_rows),
    "do not match",
    fixed = TRUE
  )

  bad_cols <- protein_object
  colnames(bad_cols@protein_quant_table) <- c("S1", "S9")
  expect_error(
    preservePeptideNaValuesHelper(peptide_object, bad_cols),
    "do not match",
    fixed = TRUE
  )
})
