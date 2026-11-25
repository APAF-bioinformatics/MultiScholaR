# KIT-226 (VAL-196): the g:Profiler results table says what its columns hold. g:Profiler returns only an adjusted
# p-value (`p_value`); the built table shows it once, as `p.adjust`, with `p_adjust_method` naming the correction as the
# answer states it (`meta$query_metadata$significance_threshold_method`, copied verbatim; "not recorded" when the
# answer does not state one), and `Count` is `intersection_size` (the list's genes in the term). The published
# per-direction TSV carries the same method column; the figure's axis and hover text say "adjusted p-value".
# T1 to T5 and G1 to G3 of the ticket, section 4. No service is called: the answers are the package's deterministic
# g:Profiler double (R/utils_workflow_dia_enrichment.R), its metadata given the method field in this file only.

#' The deterministic g:Profiler answer for one direction; its metadata states `method` when one is given.
k226Answer <- function(direction, method = NULL) {
    answer <- buildProtEnrichDeterministicGprofilerResult(direction)
    if (!is.null(method)) answer$meta$query_metadata$significance_threshold_method <- method
    answer
}

#' `EnrichmentResults` with one comparison, its up and down answers stating `up` and `down` as their method.
k226Results <- function(up = NULL, down = up) {
    results <- createEnrichmentResults(tibble::tibble(contrasts = "groupB-groupA"))
    results@enrichment_data <- list(`groupB-groupA` = list(up = k226Answer("positive", up),
        down = k226Answer("negative", down)))
    results
}

#' The table the tab builds for the comparison (the run and the restore build it this way).
k226Built <- function(up = NULL, down = up) {
    built <- buildProtEnrichAllContrastResults(enrichmentResults = k226Results(up, down),
        methodInfo = list(method = "gprofiler2"))
    built[["groupB-groupA"]]$gprofiler_results
}

k226Fixture <- function() buildProtEnrichDeterministicGprofilerResult()$result

test_that("(T1) Count is the number of the list's genes in the term (intersection_size), not the term's size", {
    built <- k226Built("g_SCS")
    expect_identical(built$directionality, rep(c("positive", "negative"), each = 5L))
    expect_identical(built$Count, rep(c(4L, 5L, 3L, 3L, 2L), 2L))
    expect_identical(built$Count, built$intersection_size)
})

test_that("(T2) the table has no pvalue or qvalue column (g:Profiler gives no unadjusted p-value or q-value)", {
    built <- k226Built("g_SCS")
    expect_false("pvalue" %in% names(built))
    expect_false("qvalue" %in% names(built))
})

test_that("(T3) p.adjust is the service's adjusted value; p_adjust_method is the correction the answer states", {
    built <- k226Built("g_SCS")
    expect_identical(built$p.adjust, built$p_value)
    expect_identical(built$p_adjust_method, rep("g_SCS", 10L))
    expect_identical(k226Built("fdr")$p_adjust_method, rep("fdr", 10L))
    expect_identical(k226Built()$p_adjust_method, rep("not recorded", 10L))
    # Each direction's rows carry its own answer's statement; an empty statement is not a method.
    expect_identical(k226Built("g_SCS", NULL)$p_adjust_method, rep(c("g_SCS", "not recorded"), each = 5L))
    expect_identical(k226Built("")$p_adjust_method, rep("not recorded", 10L))
})

test_that("(T4) an Open rebuilds the stored answer into the same columns", {
    results <- k226Results("g_SCS")
    manager <- list(getHistory = \() list("kit226_enrichment"), getState = \(state_name) results)
    saved <- list()
    testthat::local_mocked_bindings(saveDIADisplayTable = \(manager, table) {
        saved[[length(saved) + 1L]] <<- table
        table
    }, .package = "MultiScholaR")
    out <- restoreWorkflowEnrichment(manager, list(state_name = "kit226_enrichment"))
    restored <- out$all_enrichment_results[["groupB-groupA"]]$gprofiler_results
    expect_identical(out$method$method, "gprofiler2")
    expect_identical(restored$Count, rep(c(4L, 5L, 3L, 3L, 2L), 2L))
    expect_false(any(c("pvalue", "qvalue") %in% names(restored)))
    expect_identical(restored$p.adjust, restored$p_value)
    expect_identical(restored$p_adjust_method, rep("g_SCS", 10L))
    expect_identical(saved[[1L]], restored)
})

test_that("(T5) the tab's archive gprofiler2_results.tsv carries the corrected columns", {
    written <- list()
    out <- writeProtEnrichResultsDownloadArchive(file = tempfile(fileext = ".zip"), selectedContrast = "groupB-groupA",
        methodInfo = list(method = "gprofiler2", species_name = "Homo sapiens"), organismTaxid = "9606", upCutoff = 0,
        downCutoff = 0, qCutoff = 0.05, gprofilerResults = k226Built("g_SCS"), tempDir = withr::local_tempdir(),
        writeTsvFn = \(x, path) written[[basename(path)]] <<- x,
        writeLinesFn = \(text, path) written[[basename(path)]] <<- text,
        zipFn = \(...) invisible(NULL))
    expect_identical(out$status, "ok")
    tsv <- written[["gprofiler2_results.tsv"]]
    expect_false(any(c("pvalue", "qvalue") %in% names(tsv)))
    expect_identical(tsv$Count, tsv$intersection_size)
    expect_identical(unique(tsv$p_adjust_method), "g_SCS")
})

test_that("(G1) the clusterProfiler route's table is unchanged", {
    up <- buildProtEnrichDeterministicClusterProfilerResult("up")
    down <- buildProtEnrichDeterministicClusterProfilerResult("down")
    results <- createEnrichmentResults(tibble::tibble(contrasts = "groupB-groupA"))
    results@enrichment_data <- list(`groupB-groupA` = list(up = up, down = down))
    built <- buildProtEnrichAllContrastResults(enrichmentResults = results,
        methodInfo = list(method = "clusterprofiler"))[["groupB-groupA"]]
    expect_identical(as.list(built$clusterprofiler_results), as.list(rbind(up@result, down@result)))
    expect_null(built$gprofiler_results)
})

test_that("(G3) the built table keeps g:Profiler's own p_value, term_size and intersection_size", {
    built <- k226Built("g_SCS")
    fixture <- k226Fixture()
    expect_identical(built$p_value, rep(fixture$p_value, 2L))
    expect_identical(built$term_size, rep(fixture$term_size, 2L))
    expect_identical(built$intersection_size, rep(fixture$intersection_size, 2L))
    expect_identical(built$term_id, rep(fixture$term_id, 2L))
})

test_that("(G2) the published TSV gains only p_adjust_method; the figure's data and mappings are unchanged, its labels say adjusted", {
    answer <- k226Answer("positive", "g_SCS")
    fixture <- answer$result
    written <- NULL
    testthat::local_mocked_bindings(write_tsv = \(x, file, ...) written <<- x, .package = "readr")
    plots <- suppressWarnings(generate_enrichment_plots(answer, "groupB-groupA", "up", withr::local_tempdir()))
    expected <- fixture
    expected$parents <- sapply(expected$parents, paste, collapse = ", ")
    expect_identical(written[setdiff(names(written), "p_adjust_method")], expected)
    expect_identical(written$p_adjust_method, rep("g_SCS", 5L))
    written <- NULL
    suppressWarnings(generate_enrichment_plots(k226Answer("positive"), "groupB-groupA", "up", withr::local_tempdir()))
    expect_identical(written$p_adjust_method, rep("not recorded", 5L))

    plot <- plots$static
    expect_equal(plot$data$neg_log10_p, -log10(fixture$p_value))
    expect_identical(rlang::as_label(plot$mapping$x), "source")
    expect_identical(rlang::as_label(plot$mapping$y), "neg_log10_p")
    expect_identical(rlang::as_label(plot$layers[[2L]]$mapping$size), "term_size")
    expect_equal(ggplot2::layer_data(plot, 1)$yintercept, -log10(0.05))
    expect_identical(plot$labels$x, "Annotation Source")
    expect_identical(plot$labels$y, "-log10(adjusted p-value)")
    tooltips <- rlang::eval_tidy(plot$mapping$text, plot$data)
    expect_identical(tooltips, paste0("Term: ", fixture$term_name, "\n", "ID: ", fixture$term_id, "\n",
        "Adjusted p-value: ", signif(fixture$p_value, 3), "\n", "Genes: ", fixture$intersection_size))
})
