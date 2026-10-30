# KIT-069: the lipidomics and metabolomics findings the 2026-10-27 static audit carried to stage 3
# (kit-064/AUDIT-ALIGNMENT-20261028.md, the rows marked "stage 3"). Each test failed on the tree
# before its repair and passes after; the repairs are outside the kernels.
library(testthat)

skip_if_not_installed("shiny")

omicsAuditNamespace <- function() asNamespace("MultiScholaR")

# ---------------------------------------------------------------------------
# BUG-10: the QC orchestrators' auto-initialisation observers no longer invalidate themselves.
# The observer depends on the state manager only; the flag it writes is read in isolation, so the
# sub-modules are initialised once and the state is read twice in all: once by the tab renderer and
# once by the observer (three before the repair: the observer re-ran after writing its own flag).
# ---------------------------------------------------------------------------
countQcAutoInitialisation <- function(server_name, init_name, state_class) {
  ns <- omicsAuditNamespace()
  server_fn <- get(server_name, envir = ns, inherits = FALSE)
  counts <- new.env(parent = emptyenv())
  counts$state_reads <- 0L
  counts$inits <- 0L

  state_manager <- list(getState = function(state_name = NULL) {
    counts$state_reads <- counts$state_reads + 1L
    structure(list(), class = state_class)
  })

  init_mock <- function(...) {
    counts$inits <- counts$inits + 1L
    invisible(NULL)
  }
  mock_args <- stats::setNames(list(init_mock), init_name)
  do.call(testthat::local_mocked_bindings, c(mock_args, list(.package = "MultiScholaR", .env = environment())))

  shiny::testServer(
    server_fn,
    args = list(
      workflow_data = shiny::reactiveValues(state_manager = state_manager),
      experiment_paths = list(),
      omic_type = "test",
      experiment_label = "Audit",
      qc_trigger = function() NULL
    ),
    {
      session$flushReact()
      session$flushReact()
    }
  )
  list(state_reads = counts$state_reads, inits = counts$inits)
}

test_that("BUG-10: the lipidomics QC auto-initialisation observer runs once per state-manager appearance", {
  counts <- countQcAutoInitialisation("mod_lipid_qc_server", "initializeLipidQcSubmodules", "LipidomicsAssayData")
  expect_identical(counts$inits, 1L)
  expect_identical(counts$state_reads, 2L)
})

test_that("BUG-10: the metabolomics QC auto-initialisation observer runs once per state-manager appearance", {
  counts <- countQcAutoInitialisation("mod_metab_qc_server", "initializeMetabQcSubmodules", "MetaboliteAssayData")
  expect_identical(counts$inits, 1L)
  expect_identical(counts$state_reads, 2L)
})

# ---------------------------------------------------------------------------
# BUG-16: the metabolomics composite generator with no plot files. `letters[1:0]` produced a phantom
# "a)" label and `rep(1:0, ...)` a six-element grouping, so `split` warned from base R and a row was
# labelled that had no file. With `seq_len` no label is built and the empty call ends, as in the
# proteomics and lipidomics copies, in the generator's own error handler (a zero-file call cannot
# come from the module, which lists the files it wrote).
# ---------------------------------------------------------------------------
test_that("BUG-16: the metabolomics composite generator builds no phantom label for an empty file list", {
  ns <- omicsAuditNamespace()
  generate <- get("generateMetabNormCompositeFromFiles", envir = ns, inherits = FALSE)
  messages <- character()
  result <- NULL
  expect_no_warning(
    result <- generate(
      plot_files = character(0),
      ncol = 3,
      logInfoFn = function(...) invisible(NULL),
      logErrorFn = function(m, ...) messages <<- c(messages, paste("error:", m)),
      warningFn = function(m, ...) messages <<- c(messages, paste("warning:", m)),
      requireNamespaceFn = function(...) TRUE,
      buildLabelPlotFn = function(...) stop("no label plot expected"),
      buildTitlePlotFn = function(...) stop("no title plot expected"),
      loadImageAsPlotFn = function(...) stop("no image expected"),
      wrapPlotsFn = function(...) stop("nothing to wrap"),
      plotLayoutFn = function(...) NULL,
      fileExistsFn = function(...) FALSE,
      gcFn = function() NULL
    )
  )
  expect_null(result)
  expect_length(messages, 1L)
  expect_false(any(grepl("label plot|title plot|image|wrap", messages)))
})

# ---------------------------------------------------------------------------
# BUG-17: the Glimma volcano builders keep every row when every FDR is zero.
# ---------------------------------------------------------------------------
buildAllZeroFdrVolcanoInputs <- function(omic) {
  skip_if_not_installed("limma")
  skip_if_not_installed("Glimma")
  set.seed(69)
  n <- 20L
  ids <- paste0("F", seq_len(n))
  id_col <- if (omic == "lipid") "lipid_id" else "metabolite_id"
  name_col <- if (omic == "lipid") "lipid_name" else "metabolite_name"
  assay <- data.frame(ids, paste("Feature", seq_len(n)),
    S1 = stats::rnorm(n), S2 = stats::rnorm(n), S3 = stats::rnorm(n), S4 = stats::rnorm(n),
    check.names = FALSE, stringsAsFactors = FALSE)
  names(assay)[1:2] <- c(id_col, name_col)
  design <- data.frame(Run = c("S1", "S2", "S3", "S4"), group = c("A", "A", "B", "B"), stringsAsFactors = FALSE)
  long <- data.frame(ids, paste("Feature", seq_len(n)), assay = "Assay1", comparison = "groupB - groupA",
    friendly_name = "B_vs_A", logFC = stats::rnorm(n), raw_pvalue = 0, fdr_qvalue = 0,
    significant = "Up", stringsAsFactors = FALSE)
  names(long)[1:2] <- c(id_col, name_col)
  mat <- matrix(stats::rnorm(4 * n), n, 4, dimnames = list(ids, design$Run))
  mm <- stats::model.matrix(~ 0 + group, data = design)
  fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(mat, mm), limma::makeContrasts(groupB - groupA, levels = mm)))
  if (omic == "lipid") {
    object <- methods::new("LipidomicsAssayData", lipid_data = list(Assay1 = assay), design_matrix = design,
      sample_id = "Run", group_id = "group", lipid_id_column = id_col, annotation_id_column = name_col)
    list(theObject = object, da_lipids_long = long, contrasts_results = list(Assay1 = list(fit.eb = fit)))
  } else {
    object <- methods::new("MetaboliteAssayData", metabolite_data = list(Assay1 = assay), design_matrix = design,
      sample_id = "Run", group_id = "group", metabolite_id_column = id_col, annotation_id_column = name_col)
    list(theObject = object, da_metabolites_long = long, contrasts_results = list(Assay1 = list(fit.eb = fit)))
  }
}

glimmaRowCount <- function(widget) {
  if (!inherits(widget, "htmlwidget")) return(NA_integer_)
  NROW(widget$x$data$table)
}

test_that("BUG-17: the lipidomics Glimma volcano keeps every feature when every FDR is zero", {
  inputs <- buildAllZeroFdrVolcanoInputs("lipid")
  widget <- suppressMessages(generateLipidDAVolcanoPlotGlimma(inputs, "B_vs_A", "Assay1"))
  expect_s3_class(widget, "glimmaXY")
  expect_identical(glimmaRowCount(widget), 20L)
  expect_true(all(is.finite(widget$x$data$table$negLog10FDR)))
})

test_that("BUG-17: the metabolomics Glimma volcano keeps every feature when every FDR is zero", {
  inputs <- buildAllZeroFdrVolcanoInputs("metab")
  widget <- suppressMessages(generateMetabDAVolcanoPlotGlimma(inputs, "B_vs_A", "Assay1"))
  expect_s3_class(widget, "glimmaXY")
  expect_identical(glimmaRowCount(widget), 20L)
  expect_true(all(is.finite(widget$x$data$table$negLog10FDR)))
})

# ---------------------------------------------------------------------------
# BUG-23: sapply in subscript and sum contexts returns a list on empty input and aborts; vapply does not.
# ---------------------------------------------------------------------------
test_that("BUG-23: the sample-column resolvers accept an assay table with no columns", {
  ns <- omicsAuditNamespace()
  empty <- data.frame()
  lipid <- get("resolveLipidImportSampleColumns", envir = ns, inherits = FALSE)
  metab <- get("resolveMetabImportSampleColumns", envir = ns, inherits = FALSE)
  expect_identical(lipid(empty, vendorFormat = "lipidsearch"), character(0))
  expect_identical(metab(empty, vendorFormat = "msdial"), character(0))
})

test_that("BUG-23: the lipidomics import finaliser logs an empty assay list", {
  ns <- omicsAuditNamespace()
  record <- get("finalizeLipidImportSetupState", envir = ns, inherits = FALSE)
  workflow <- new.env(parent = emptyenv())
  workflow$processing_log <- list()
  workflow$tab_status <- list(design_matrix = "disabled")
  logged <- character()
  record(workflow, assayList = list(), detectedFormat = "lipidsearch", lipidIdCol = "lipid_id",
    sampleColumns = character(0), now = function() 0, logInfo = function(m) logged <<- c(logged, m))
  expect_identical(logged, "Lipidomics import complete: 0 assays, 0 total lipids")
})

test_that("BUG-23: the metabolomics duplicate resolver reports an assay with no columns instead of aborting", {
  ns <- omicsAuditNamespace()
  resolve <- get("resolveMetabDuplicateAssayData", envir = ns, inherits = FALSE)
  warned <- character()
  result <- resolve(list(Empty = data.frame()), metaboliteIdCol = "metabolite_id",
    resolveDuplicateFeaturesByIntensityFn = function(...) stop("not reached"),
    logWarnFn = function(m) warned <<- c(warned, m))
  expect_identical(warned, "No numeric columns found in assay: Empty")
})
