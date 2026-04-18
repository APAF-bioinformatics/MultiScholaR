buildLipidDaSummaryStatsText <- function(
    daResults,
    selectedContrast,
    selectedAssay,
    daQValThresh
) {
    if (is.null(daResults) || nrow(daResults) == 0) {
        return("No results available.")
    }

    results <- daResults

    if (!is.null(selectedContrast)) {
        results <- results[results$comparison == selectedContrast |
            results$friendly_name == selectedContrast, ]
    }

    if (!is.null(selectedAssay) && selectedAssay != "All") {
        results <- results[results$assay == selectedAssay, ]
    }

    total <- nrow(results)
    sig <- sum(results$significant != "NS", na.rm = TRUE)
    up <- sum(results$significant == "Up", na.rm = TRUE)
    down <- sum(results$significant == "Down", na.rm = TRUE)

    paste(
        c(
            sprintf("Total lipids: %d", total),
            sprintf(
                "Significant (Q < %.3f): %d (%.1f%%)",
                daQValThresh, sig, 100 * sig / max(total, 1)
            ),
            sprintf("  Up-regulated: %d", up),
            sprintf("  Down-regulated: %d", down)
        ),
        collapse = "\n"
    )
}

buildLipidDaSummaryStatsRenderText <- function(
    daResultsList,
    selectedContrast,
    selectedAssay,
    daQValThresh,
    reqFn = shiny::req,
    buildSummaryStatsTextFn = buildLipidDaSummaryStatsText
) {
    reqFn(daResultsList)

    buildSummaryStatsTextFn(
        daResults = daResultsList$da_lipids_long,
        selectedContrast = selectedContrast,
        selectedAssay = selectedAssay,
        daQValThresh = daQValThresh
    )
}

buildLipidDaResultsTableRenderWidget <- function(
    daResultsList,
    selectedContrast,
    selectedAssay,
    tableSignificance,
    daQValThresh,
    lfcThreshold,
    tableMaxRows,
    reqFn = shiny::req,
    buildResultsTableWidgetFn = buildLipidDaResultsTableWidget
) {
    reqFn(daResultsList)

    buildResultsTableWidgetFn(
        daResults = daResultsList$da_lipids_long,
        selectedContrast = selectedContrast,
        selectedAssay = selectedAssay,
        tableSignificance = tableSignificance,
        daQValThresh = daQValThresh,
        lfcThreshold = lfcThreshold,
        tableMaxRows = tableMaxRows
    )
}

buildLipidDaResultsTableWidget <- function(
    daResults,
    selectedContrast,
    selectedAssay,
    tableSignificance,
    daQValThresh,
    lfcThreshold,
    tableMaxRows,
    datatableFactory = DT::datatable,
    formatRoundFactory = DT::formatRound,
    formatStyleFactory = DT::formatStyle,
    styleEqualFactory = DT::styleEqual
) {
    if (is.null(daResults) || nrow(daResults) == 0) {
        return(NULL)
    }

    results <- daResults

    if (!is.null(selectedContrast)) {
        results <- results[results$comparison == selectedContrast |
            results$friendly_name == selectedContrast, ]
    }

    if (!is.null(selectedAssay) && selectedAssay != "All") {
        results <- results[results$assay == selectedAssay, ]
    }

    if (identical(tableSignificance, "significant")) {
        results <- results[results$fdr_qvalue < daQValThresh, ]
    } else if (identical(tableSignificance, "up")) {
        results <- results[results$fdr_qvalue < daQValThresh &
            results$logFC > lfcThreshold, ]
    } else if (identical(tableSignificance, "down")) {
        results <- results[results$fdr_qvalue < daQValThresh &
            results$logFC < -lfcThreshold, ]
    }

    if (nrow(results) > tableMaxRows) {
        results <- results[seq_len(tableMaxRows), ]
    }

    displayCols <- c(
        "lipid_id", "lipid_name", "assay", "logFC",
        "raw_pvalue", "fdr_qvalue", "significant"
    )
    displayCols <- intersect(displayCols, colnames(results))
    results <- results[, displayCols, drop = FALSE]

    datatableFactory(
        results,
        options = list(
            pageLength = 25,
            scrollX = TRUE,
            dom = "Bfrtip",
            buttons = c("copy", "csv", "excel")
        ),
        extensions = "Buttons",
        rownames = FALSE,
        filter = "top"
    ) |>
        formatRoundFactory(columns = c("logFC"), digits = 3) |>
        formatRoundFactory(columns = c("raw_pvalue", "fdr_qvalue"), digits = 6) |>
        formatStyleFactory(
            "significant",
            backgroundColor = styleEqualFactory(
                c("Up", "Down", "NS"),
                c("#ffcccc", "#cce5ff", "white")
            )
        )
}

buildLipidDaResultsDownloadHandler <- function(
    daResultsReactive,
    downloadHandlerFactory = shiny::downloadHandler,
    currentDate = Sys.Date,
    csvWriter = utils::write.csv,
    reqFn = shiny::req
) {
    downloadHandlerFactory(
        filename = function() {
            paste0("lipidomics_da_results_", currentDate(), ".csv")
        },
        content = function(file) {
            results <- daResultsReactive()
            reqFn(results)
            csvWriter(results, file, row.names = FALSE)
        }
    )
}

resolveLipidDaResultsDownloadData <- function(daResultsList) {
    daResultsList$da_lipids_long
}

buildLipidDaResultsDownloadOutputHandler <- function(
    daResultsList,
    buildResultsDownloadHandlerFn = buildLipidDaResultsDownloadHandler,
    resolveResultsDownloadDataFn = resolveLipidDaResultsDownloadData
) {
    buildResultsDownloadHandlerFn(
        daResultsReactive = function() {
            resolveResultsDownloadDataFn(daResultsList)
        }
    )
}

