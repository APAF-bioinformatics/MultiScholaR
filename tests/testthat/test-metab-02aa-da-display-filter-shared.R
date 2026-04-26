# fidelity-coverage-compare: shared
library(testthat)

repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)

skipIfMissingMetabDaDisplaySplitFiles <- function() {
  required_paths <- c(
    "R/mod_metab_da_display_helpers.R",
    "R/mod_metab_da_registration_helpers.R",
    "R/mod_metab_da_server_helpers.R"
  )
  missing <- required_paths[!file.exists(file.path(repo_root, required_paths))]
  if (length(missing) > 0) {
    testthat::skip(
      sprintf(
        "Target-only metab DA display file(s) not present: %s",
        paste(basename(missing), collapse = ", ")
      )
    )
  }
}

skipIfMissingMetabDaDisplaySplitFiles()

buildMetabDaDisplayResults <- function() {
  data.frame(
    metabolite_id = c("M1", "M2", "M3", "M4"),
    metabolite_name = c("Met One", "Met Two", "Met Three", "Met Four"),
    comparison = c("groupB-groupA", "groupB-groupA", "groupC-groupA", "groupB-groupA"),
    friendly_name = c("B vs A", "B vs A", "C vs A", "B vs A"),
    assay = c("LCMS_Pos", "LCMS_Neg", "LCMS_Pos", "LCMS_Pos"),
    logFC = c(1.2, -1.1, -0.9, 0.2),
    raw_pvalue = c(0.001, 0.004, 0.02, 0.3),
    fdr_qvalue = c(0.01, 0.02, 0.03, 0.20),
    significant = c("Up", "Down", "Down", "NS"),
    stringsAsFactors = FALSE
  )
}

test_that("metabolomics DA display filter preserves contrast and assay narrowing", {
  results <- buildMetabDaDisplayResults()

  filtered <- filterMetabDaDisplayResults(
    results = results,
    selectedContrast = "B vs A",
    selectedAssay = "LCMS_Pos"
  )

  expect_equal(filtered$metabolite_id, c("M1", "M4"))
  expect_true(all(filtered$friendly_name == "B vs A"))
  expect_true(all(filtered$assay == "LCMS_Pos"))
})

test_that("metabolomics DA display filter preserves significance branches and row limits", {
  results <- buildMetabDaDisplayResults()

  significant <- filterMetabDaDisplayResults(
    results = results,
    significanceFilter = "significant",
    daQValThresh = 0.05
  )
  expect_equal(significant$metabolite_id, c("M1", "M2", "M3"))

  up <- filterMetabDaDisplayResults(
    results = results,
    significanceFilter = "up",
    daQValThresh = 0.05,
    treatLfcCutoff = 0.5
  )
  expect_equal(up$metabolite_id, "M1")

  down <- filterMetabDaDisplayResults(
    results = results,
    significanceFilter = "down",
    daQValThresh = 0.05,
    treatLfcCutoff = 0.5
  )
  expect_equal(down$metabolite_id, c("M2", "M3"))

  limited <- filterMetabDaDisplayResults(
    results = results,
    maxRows = 2
  )
  expect_equal(limited$metabolite_id, c("M1", "M2"))
})

test_that("metabolomics DA display summary preserves counts and threshold", {
  summary_stats <- summarizeMetabDaDisplayResults(
    results = buildMetabDaDisplayResults(),
    daQValThresh = 0.05
  )

  expect_equal(summary_stats$total, 4)
  expect_equal(summary_stats$significant, 3)
  expect_equal(summary_stats$upRegulated, 1)
  expect_equal(summary_stats$downRegulated, 2)
  expect_equal(summary_stats$significantPct, 75)
  expect_equal(summary_stats$qValueThreshold, 0.05)
})

test_that("metabolomics DA summary text preserves empty and formatted branches", {
  expect_identical(
    buildMetabDaSummaryStatsText(NULL),
    "No results available."
  )

  empty <- buildMetabDaDisplayResults()[0, ]
  expect_identical(
    buildMetabDaSummaryStatsText(empty),
    "No results available."
  )

  summary_text <- buildMetabDaSummaryStatsText(
    results = buildMetabDaDisplayResults(),
    daQValThresh = 0.05
  )

  expect_identical(
    summary_text,
    paste(
      c(
        "Total metabolites: 4",
        "Significant (Q < 0.050): 3 (75.0%)",
        "  Up-regulated: 1",
        "  Down-regulated: 2"
      ),
      collapse = "\n"
    )
  )
})

test_that("metabolomics DA results table helpers preserve display columns and widget options", {
  display_results <- shapeMetabDaTableDisplayResults(buildMetabDaDisplayResults())

  expect_equal(
    colnames(display_results),
    c(
      "metabolite_id", "metabolite_name", "assay", "logFC",
      "raw_pvalue", "fdr_qvalue", "significant"
    )
  )

  table_widget <- buildMetabDaResultsDatatable(display_results)
  column_defs <- table_widget$x$options$columnDefs
  named_defs <- Filter(function(def) !is.null(def$name), column_defs)

  expect_s3_class(table_widget, "datatables")
  expect_equal(table_widget$x$filter, "top")
  expect_equal(unlist(table_widget$x$extensions), "Buttons")
  expect_equal(table_widget$x$options$pageLength, 25)
  expect_true(isTRUE(table_widget$x$options$scrollX))
  expect_equal(table_widget$x$options$dom, "Bfrtip")
  expect_equal(unlist(table_widget$x$options$buttons), c("copy", "csv", "excel"))
  expect_setequal(
    vapply(named_defs, `[[`, character(1), "name"),
    c(
      "metabolite_id", "metabolite_name", "assay", "logFC",
      "raw_pvalue", "fdr_qvalue", "significant"
    )
  )
})

test_that("metabolomics DA render helpers preserve req and handoff contracts", {
  captured_req <- NULL
  captured_filter <- NULL
  captured_text <- NULL
  da_results_list <- list(da_metabolites_long = "raw-results")

  summary_output <- buildMetabDaSummaryStatsRenderOutput(
    daResultsList = da_results_list,
    selectedContrast = "B vs A",
    selectedAssay = "LCMS_Pos",
    daQValThresh = 0.01,
    requireInputs = function(...) {
      captured_req <<- list(...)
      invisible(TRUE)
    },
    filterResults = function(...) {
      captured_filter <<- list(...)
      "filtered-results"
    },
    buildSummaryText = function(results, daQValThresh = 0.05) {
      captured_text <<- list(results = results, daQValThresh = daQValThresh)
      "summary-stats-output"
    },
    writeOutput = function(text, ...) {
      paste0("written:", text)
    }
  )

  expect_identical(captured_req, list(da_results_list))
  expect_identical(captured_filter$results, "raw-results")
  expect_identical(captured_filter$selectedContrast, "B vs A")
  expect_identical(captured_filter$selectedAssay, "LCMS_Pos")
  expect_identical(captured_text$results, "filtered-results")
  expect_identical(captured_text$daQValThresh, 0.01)
  expect_identical(summary_output, "written:summary-stats-output")

  captured_shape <- NULL
  captured_datatable <- NULL
  table_output <- buildMetabDaResultsTableRenderOutput(
    daResultsList = da_results_list,
    selectedContrast = "B vs A",
    selectedAssay = "LCMS_Pos",
    significanceFilter = "up",
    daQValThresh = 0.01,
    treatLfcCutoff = 1.5,
    maxRows = 7,
    requireInputs = function(...) invisible(TRUE),
    filterResults = function(...) data.frame(filtered = "yes", stringsAsFactors = FALSE),
    shapeResults = function(results) {
      captured_shape <<- results
      data.frame(shaped = "yes", stringsAsFactors = FALSE)
    },
    buildDatatable = function(results) {
      captured_datatable <<- results
      "datatable-output"
    }
  )

  expect_identical(captured_shape, data.frame(filtered = "yes", stringsAsFactors = FALSE))
  expect_identical(captured_datatable, data.frame(shaped = "yes", stringsAsFactors = FALSE))
  expect_identical(table_output, "datatable-output")
})

test_that("metabolomics DA cluster summary helpers preserve grouped text output", {
  clusters <- c(MetC = 2, MetA = 1, MetB = 2)

  summary_text <- buildMetabDaClusterSummaryText(clusters)
  expect_match(summary_text, "Total Clusters: 2", fixed = TRUE)
  expect_match(summary_text, "Cluster 1 \\(1 metabolites\\):")
  expect_match(summary_text, "MetA", fixed = TRUE)
  expect_match(summary_text, "Cluster 2 \\(2 metabolites\\):")

  captured_req <- NULL
  captured_text <- NULL
  output <- buildMetabDaClusterSummaryRenderOutput(
    treeCutMethod = "dynamic",
    clusters = c(MetA = 1, MetB = 2),
    requireInputs = function(...) {
      captured_req <<- list(...)
      invisible(TRUE)
    },
    buildSummaryText = function(clusters, maxMembers = 20) {
      captured_text <<- list(clusters = clusters, maxMembers = maxMembers)
      "cluster-summary-output"
    },
    writeOutput = function(text, ...) {
      paste0("written:", text)
    }
  )

  expect_identical(captured_req, list(TRUE))
  expect_identical(captured_text$clusters, c(MetA = 1, MetB = 2))
  expect_identical(captured_text$maxMembers, 20)
  expect_identical(output, "written:cluster-summary-output")
})
