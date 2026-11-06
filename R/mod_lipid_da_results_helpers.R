buildLipidDaSummaryStatsText <- function(
    daResults,
    selectedContrast,
    selectedAssay,
    daQValThresh,
    featureNoun = "lipids",
    significanceLabel = "Q"
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
            sprintf("Total %s: %d", featureNoun, total),
            sprintf(
                "Significant (%s < %.3f): %d (%.1f%%)",
                significanceLabel, daQValThresh, sig, 100 * sig / max(total, 1)
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
    buildSummaryStatsTextFn = buildLipidDaSummaryStatsText,
    longTable = "da_lipids_long"
) {
    reqFn(daResultsList)

    buildSummaryStatsTextFn(
        daResults = daResultsList[[longTable]],
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
    buildResultsTableWidgetFn = buildLipidDaResultsTableWidget,
    longTable = "da_lipids_long"
) {
    reqFn(daResultsList)

    buildResultsTableWidgetFn(
        daResults = daResultsList[[longTable]],
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
    styleEqualFactory = DT::styleEqual,
    columns = NULL
) {
    if (is.null(daResults) || nrow(daResults) == 0) {
        return(NULL)
    }
    # KIT-097: an omic profile's identifier, label, p-value and significance columns (lipidomics by default); the
    # significance column keeps its name in the data and is shown under its label.
    columns <- utils::modifyList(list(id = "lipid_id", label = "lipid_name", pvalue = "raw_pvalue", significance = "fdr_qvalue",
        significance_label = NULL), columns %||% list())
    significance <- columns$significance

    results <- daResults

    if (!is.null(selectedContrast)) {
        results <- results[results$comparison == selectedContrast |
            results$friendly_name == selectedContrast, ]
    }

    if (!is.null(selectedAssay) && selectedAssay != "All") {
        results <- results[results$assay == selectedAssay, ]
    }

    if (identical(tableSignificance, "significant")) {
        results <- results[results[[significance]] < daQValThresh, ]
    } else if (identical(tableSignificance, "up")) {
        results <- results[results[[significance]] < daQValThresh &
            results$logFC > lfcThreshold, ]
    } else if (identical(tableSignificance, "down")) {
        results <- results[results[[significance]] < daQValThresh &
            results$logFC < -lfcThreshold, ]
    }

    if (nrow(results) == 0) {
        return(NULL)
    }

    if (nrow(results) > tableMaxRows) {
        results <- results[seq_len(tableMaxRows), ]
    }

    displayCols <- c(
        columns$id, columns$label, "assay", "logFC",
        columns$pvalue, significance, "significant"
    )
    displayCols <- intersect(displayCols, colnames(results))
    results <- results[, displayCols, drop = FALSE]
    widgetArgs <- list(
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
    )
    if (!is.null(columns$significance_label)) {
        widgetArgs$colnames <- ifelse(displayCols == significance, columns$significance_label, displayCols)
    }

    do.call(datatableFactory, widgetArgs) |>
        formatRoundFactory(columns = c("logFC"), digits = 3) |>
        formatRoundFactory(columns = c(columns$pvalue, significance), digits = 6) |>
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
    reqFn = shiny::req,
    filePrefix = "lipidomics"
) {
    downloadHandlerFactory(
        filename = function() {
            paste0(filePrefix, "_da_results_", currentDate(), ".csv")
        },
        content = function(file) {
            results <- daResultsReactive()
            reqFn(results)
            # KIT-082: a managed tab holds no per-sample columns and writes the download from the store.
            if (is.function(results)) return(results(file))
            csvWriter(results, file, row.names = FALSE)
        }
    )
}

resolveLipidDaResultsDownloadData <- function(daResultsList, longTable = "da_lipids_long") {
    daResultsList$write_download %||% daResultsList[[longTable]]
}

buildLipidDaResultsDownloadOutputHandler <- function(
    daResultsList = NULL,
    buildResultsDownloadHandlerFn = buildLipidDaResultsDownloadHandler,
    resolveResultsDownloadDataFn = resolveLipidDaResultsDownloadData,
    readResultsFn = function() daResultsList
) {
    buildResultsDownloadHandlerFn(
        daResultsReactive = function() {
            resolveResultsDownloadDataFn(readResultsFn())
        }
    )
}
