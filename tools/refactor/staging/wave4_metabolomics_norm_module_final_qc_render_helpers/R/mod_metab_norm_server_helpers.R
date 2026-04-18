buildMetabNormCorrelationFilterSummary <- function(corrResults, filteredObject = NULL, originalObject = NULL) {
    if (is.null(corrResults) || length(corrResults) == 0) {
        return("No correlation results available.")
    }

    summaryLines <- c("=== Correlation Filtering Summary ===\n")

    for (assay_name in names(corrResults)) {
        assay_corr <- corrResults[[assay_name]]
        if (!is.null(assay_corr) && nrow(assay_corr) > 0) {
            n_pairs <- nrow(assay_corr)
            mean_corr <- round(mean(assay_corr$pearson_correlation, na.rm = TRUE), 3)
            min_corr <- round(min(assay_corr$pearson_correlation, na.rm = TRUE), 3)
            max_corr <- round(max(assay_corr$pearson_correlation, na.rm = TRUE), 3)

            summaryLines <- c(summaryLines, sprintf(
                "\n[%s]\n  Sample pairs: %d\n  Correlation: mean=%.3f, min=%.3f, max=%.3f",
                assay_name, n_pairs, mean_corr, min_corr, max_corr
            ))
        }
    }

    if (!is.null(originalObject) && !is.null(filteredObject)) {
        original_samples <- nrow(originalObject@design_matrix)
        filtered_samples <- nrow(filteredObject@design_matrix)
        removed <- original_samples - filtered_samples

        summaryLines <- c(summaryLines, sprintf(
            "\n\n[Sample Filtering]\n  Original: %d samples\n  After filtering: %d samples\n  Removed: %d samples",
            original_samples, filtered_samples, removed
        ))
    }

    paste(summaryLines, collapse = "")
}

resolveMetabNormFinalQcRenderState <- function(correlationFilteredObject = NULL, ruvCorrectedObject = NULL, postNormObject = NULL) {
    sourceObject <- if (!is.null(correlationFilteredObject)) {
        correlationFilteredObject
    } else if (!is.null(ruvCorrectedObject)) {
        ruvCorrectedObject
    } else {
        postNormObject
    }

    if (is.null(sourceObject)) {
        return(list(
            sourceObject = NULL
            , sourceStage = "empty"
            , isFallback = TRUE
            , plot = ggplot2::ggplot() +
                ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data", size = 6) +
                ggplot2::theme_void()
        ))
    }

    list(
        sourceObject = sourceObject
        , sourceStage = if (!is.null(correlationFilteredObject)) {
            "correlation_filter"
        } else if (!is.null(ruvCorrectedObject)) {
            "ruv_corrected"
        } else {
            "post_norm"
        }
        , isFallback = FALSE
        , plot = NULL
    )
}

buildMetabNormFinalQcPcaPlot <- function(sourceObject, colorVar = NULL, shapeVar = NULL, plotPcaFn = plotPca, wrapPlotsFn = patchwork::wrap_plots) {
    tryCatch({
        pcaPlots <- plotPcaFn(
            sourceObject
            , grouping_variable = colorVar
            , shape_variable = shapeVar
            , title = "Final QC - PCA"
        )

        if (is.list(pcaPlots) && length(pcaPlots) > 1) {
            wrapPlotsFn(pcaPlots, ncol = 1)
        } else if (is.list(pcaPlots) && length(pcaPlots) == 1) {
            pcaPlots[[1]]
        } else if (inherits(pcaPlots, "ggplot")) {
            pcaPlots
        } else {
            ggplot2::ggplot() + ggplot2::theme_void()
        }
    }, error = function(e) {
        ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message), size = 4) +
            ggplot2::theme_void()
    })
}

