# fidelity-coverage-compare: shared
library(testthat)

localFakeLimpa <- function(env = parent.frame()) {
  old_lib_paths <- .libPaths()
  old_options <- options(
    multischolar.fake_limpa.dpc = NULL,
    multischolar.fake_limpa.dpcImpute = NULL,
    multischolar.fake_limpa.dpcQuant = NULL
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
        "  replacement <- if (missing(dpc.slope)) 1 else dpc.slope",
        "  E[is.na(E)] <- replacement",
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
        "  stop('No fake dpcQuant handler configured')",
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

newPeptideLimpaObject <- function(values = matrix(
                                    c(NA_real_, 120, 240, NA_real_),
                                    nrow = 2,
                                    byrow = TRUE,
                                    dimnames = list(c("PEP1", "PEP2"), c("S1", "S2"))
                                  ),
                                  calculate_matrix = TRUE,
                                  is_logged_data = FALSE,
                                  args = list()) {
  peptide_rows <- expand.grid(
    peptide_index = seq_len(nrow(values)),
    Run = colnames(values),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  peptide_rows$Protein.Ids <- paste0("P", peptide_rows$peptide_index)
  peptide_rows$Stripped.Sequence <- rownames(values)[peptide_rows$peptide_index]
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

  peptide_object <- new(
    "PeptideQuantitativeData",
    peptide_data = peptide_rows,
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

test_that("PeptideQuantitativeData limpa imputation transforms raw data and stores DPC metadata", {
  localFakeLimpa()
  captured <- new.env(parent = emptyenv())
  options(
    multischolar.fake_limpa.dpc = function(y) {
      captured$dpc_input <- y
      list(dpc = c(0.2, 0.8), y = y)
    },
    multischolar.fake_limpa.dpcImpute = function(y, dpc) {
      captured$impute_input <- y
      captured$dpc <- dpc
      E <- y
      E[is.na(E)] <- log2(9)
      list(E = E)
    }
  )

  peptide_object <- newPeptideLimpaObject()

  imputed <- peptideMissingValueImputationLimpa(
    peptide_object,
    imputed_value_column = "Peptide.Imputed.Test",
    use_log2_transform = TRUE,
    verbose = TRUE,
    ensure_matrix = TRUE
  )

  expect_equal(captured$dpc_input["P1%PEP1", "S2"], log2(121))
  expect_equal(captured$dpc$dpc, c(0.2, 0.8))
  expect_equal(
    imputed@peptide_data$Peptide.Imputed.Test[
      imputed@peptide_data$Protein.Ids == "P1" & imputed@peptide_data$Run == "S1"
    ],
    8
  )
  expect_identical(imputed@norm_quantity_column, "Peptide.Imputed.Test")
  expect_equal(imputed@args$peptideMissingValueImputationLimpa$imputed_value_column, "Peptide.Imputed.Test")
  expect_true(imputed@args$peptideMissingValueImputationLimpa$use_log2_transform)
  expect_true(imputed@args$peptideMissingValueImputationLimpa$verbose)
  expect_equal(imputed@args$limpa_dpc_results$dpc_parameters, c(0.2, 0.8))
  expect_identical(
    imputed@args$limpa_dpc_results$slope_interpretation,
    "strong intensity-dependent missing"
  )
  expect_identical(imputed@args$limpa_dpc_results$dpc_method, "limpa_dpc")
  expect_equal(imputed@args$limpa_dpc_results$missing_percentage_after, 0)
})

test_that("PeptideQuantitativeData limpa imputation calculates a missing matrix and respects logged data", {
  localFakeLimpa()
  captured <- new.env(parent = emptyenv())
  options(
    multischolar.fake_limpa.dpc = function(y) {
      captured$dpc_input <- y
      list(dpc = c(0.1, 0.2), y = y)
    },
    multischolar.fake_limpa.dpcImpute = function(y, dpc) {
      E <- y
      E[is.na(E)] <- 3
      list(E = E)
    }
  )

  peptide_object <- newPeptideLimpaObject(
    values = matrix(
      c(NA_real_, 4, 5, 6),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("PEP1", "PEP2"), c("S1", "S2"))
    ),
    calculate_matrix = FALSE,
    is_logged_data = TRUE
  )

  imputed <- peptideMissingValueImputationLimpa(
    peptide_object,
    use_log2_transform = TRUE,
    verbose = TRUE,
    ensure_matrix = TRUE
  )

  expect_equal(captured$dpc_input["P1%PEP1", "S2"], 4)
  expect_equal(imputed@peptide_matrix["P1%PEP1", "S1"], 3)
  expect_identical(
    imputed@args$limpa_dpc_results$slope_interpretation,
    "nearly random missing"
  )
})

test_that("PeptideQuantitativeData limpa imputation preserves no-extra-transform branches", {
  localFakeLimpa()
  captured <- new.env(parent = emptyenv())
  call_index <- 0L
  options(
    multischolar.fake_limpa.dpc = function(y) {
      call_index <<- call_index + 1L
      captured[[paste0("dpc_input_", call_index)]] <- y
      slope <- if (call_index == 1L) 0.5 else 1.5
      list(dpc = c(0.2, slope), y = y)
    },
    multischolar.fake_limpa.dpcImpute = function(y, dpc) {
      E <- y
      E[is.na(E)] <- 2
      list(E = E)
    }
  )

  raw_object <- newPeptideLimpaObject(
    values = matrix(
      c(NaN, 15, Inf, NA_real_),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("PEP1", "PEP2"), c("S1", "S2"))
    )
  )
  raw_imputed <- peptideMissingValueImputationLimpa(
    raw_object,
    use_log2_transform = FALSE,
    verbose = TRUE
  )

  logged_object <- newPeptideLimpaObject(
    values = matrix(
      c(1, NA_real_, 3, 4),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("PEP1", "PEP2"), c("S1", "S2"))
    ),
    is_logged_data = TRUE
  )
  logged_imputed <- peptideMissingValueImputationLimpa(
    logged_object,
    use_log2_transform = FALSE,
    verbose = TRUE
  )

  expect_true(is.na(captured$dpc_input_1["P1%PEP1", "S1"]))
  expect_true(is.na(captured$dpc_input_1["P2%PEP2", "S1"]))
  expect_equal(raw_imputed@peptide_matrix["P1%PEP1", "S1"], 2)
  expect_equal(captured$dpc_input_2["P2%PEP2", "S2"], 4)
  expect_equal(logged_imputed@peptide_matrix["P1%PEP1", "S2"], 2)
  expect_identical(
    raw_imputed@args$limpa_dpc_results$slope_interpretation,
    "moderate intensity-dependent missing"
  )
  expect_identical(
    logged_imputed@args$limpa_dpc_results$slope_interpretation,
    "very strong intensity-dependent missing (approaching left-censoring)"
  )
})

test_that("PeptideQuantitativeData limpa imputation wraps limpa failures", {
  localFakeLimpa()
  options(
    multischolar.fake_limpa.dpc = function(y) {
      stop("mock dpc failure")
    }
  )

  peptide_object <- newPeptideLimpaObject()

  expect_error(
    peptideMissingValueImputationLimpa(peptide_object, verbose = TRUE),
    "limpa imputation failed: mock dpc failure",
    fixed = TRUE
  )
})

test_that("limpa imputation round-trips DIA group identities through opaque feature keys", {
  localFakeLimpa()
  options(
    multischolar.fake_limpa.dpc = function(y) list(dpc = c(0.2, 0.8), y = y),
    multischolar.fake_limpa.dpcImpute = function(y, dpc) {
      E <- y
      E[is.na(E)] <- log2(33)
      list(E = E)
    }
  )

  peptide_object <- newPeptideLimpaObject()
  peptide_object@peptide_data$Protein.Group <- ifelse(
    peptide_object@peptide_data$Stripped.Sequence == "PEP1",
    "GROUP%ONE",
    "GROUP_TWO"
  )
  peptide_object@peptide_data$Protein.Ids <- "P_SHARED"
  peptide_object@peptide_data$Stripped.Sequence[
    peptide_object@peptide_data$Protein.Group == "GROUP%ONE"
  ] <- "PEP%ONE"
  peptide_object@protein_id_column <- "Protein.Group"
  peptide_object@args <- list()
  peptide_object <- calcPeptideMatrix(peptide_object)

  imputed <- peptideMissingValueImputationLimpa(
    peptide_object,
    imputed_value_column = "Peptide.Imputed.Group",
    use_log2_transform = TRUE,
    verbose = FALSE
  )

  expect_true(all(grepl("^.multischolar_peptide_", rownames(imputed@peptide_matrix))))
  expect_identical(
    imputed@args$peptide_feature_key_map,
    peptide_object@args$peptide_feature_key_map
  )
  expect_equal(
    imputed@peptide_data$Peptide.Imputed.Group[
      imputed@peptide_data$Protein.Group == "GROUP%ONE" &
        imputed@peptide_data$Run == "S1"
    ],
    32
  )
  expect_setequal(unique(imputed@peptide_data$Protein.Group), c("GROUP%ONE", "GROUP_TWO"))
  expect_identical(unique(imputed@peptide_data$Protein.Ids), "P_SHARED")
})

test_that("limpa DPC-Quant receives Protein.Group from the feature map without parsing row names", {
  localFakeLimpa()
  captured <- new.env(parent = emptyenv())
  options(
    multischolar.fake_limpa.dpc = function(y) list(dpc = c(0.2, 0.8), y = y),
    multischolar.fake_limpa.dpcQuant = function(y,
                                               protein.id,
                                               dpc = NULL,
                                               dpc.slope = NULL,
                                               verbose,
                                               chunk) {
      captured$genes <- y$genes
      protein_ids <- unique(y$genes[[protein.id]])
      E <- matrix(
        seq_len(length(protein_ids) * ncol(y$E)),
        nrow = length(protein_ids),
        dimnames = list(protein_ids, colnames(y$E))
      )
      list(
        E = E,
        genes = data.frame(protein.id = protein_ids, stringsAsFactors = FALSE),
        other = list(
          standard.error = matrix(0, nrow(E), ncol(E), dimnames = dimnames(E)),
          n.observations = matrix(1L, nrow(E), ncol(E), dimnames = dimnames(E))
        )
      )
    }
  )

  peptide_object <- newPeptideLimpaObject(
    values = matrix(
      c(10, 20, 30, 40),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("FEATURE_ONE", "FEATURE_TWO"), c("S1", "S2"))
    )
  )
  peptide_object@peptide_data$Protein.Group <- ifelse(
    peptide_object@peptide_data$Stripped.Sequence == "FEATURE_ONE",
    "GROUP%ONE",
    "GROUP_TWO"
  )
  peptide_object@peptide_data$Protein.Ids <- "P_SHARED"
  peptide_object@peptide_data$Stripped.Sequence <- "SHARED_SEQUENCE"
  peptide_object@protein_id_column <- "Protein.Group"
  peptide_object@args <- list()
  peptide_object <- calcPeptideMatrix(peptide_object)

  protein_method <- methods::selectMethod(
    "proteinMissingValueImputationLimpa",
    "PeptideQuantitativeData"
  )
  protein_method_env <- list2env(
    list(
      ProteinQuantitativeData = function(...) {
        methods::new("ProteinQuantitativeData", ...)
      }
    ),
    parent = environment(protein_method)
  )
  environment(protein_method) <- protein_method_env

  quantified <- protein_method(
    peptide_object,
    dpc_slope = 0.8,
    verbose = FALSE
  )

  expect_identical(captured$genes$protein.id, c("GROUP%ONE", "GROUP_TWO"))
  expect_false(any(grepl("SHARED_SEQUENCE", captured$genes$protein.id, fixed = TRUE)))
  expect_identical(quantified@protein_id_column, "Protein.Group")
  expect_identical(quantified@protein_quant_table$Protein.Group, c("GROUP%ONE", "GROUP_TWO"))
  expect_true(all(c("Protein.Group", "Protein.Ids") %in% names(quantified@protein_id_table)))
  expect_identical(unique(quantified@protein_id_table$Protein.Ids), "P_SHARED")
})
