# Keep plot-aesthetics resolution top-level so a later extraction wave can move
# it without reopening the module server wrapper.
getPlotAesthetics <- function(
    colorVariable = NULL,
    shapeVariable = NULL
) {
    list(
        color_var = if (is.null(colorVariable) || colorVariable == "") {
            "group"
        } else {
            colorVariable
        }
        , shape_var = if (is.null(shapeVariable) || shapeVariable == "") {
            "group"
        } else {
            shapeVariable
        }
    )
}

# Keep assay-label rendering top-level so a later extraction wave can move it
# without reopening the module server wrapper's static label bindings.
renderMetabNormAssayLabel <- function(
    assaySlot,
    getAssayNamesFn,
    renderTextFn = shiny::renderText
) {
    renderTextFn({
        assayNames <- getAssayNamesFn()

        if (!is.null(assayNames) && length(assayNames) >= assaySlot) {
            paste("Assay:", assayNames[[assaySlot]])
        } else {
            paste0("Assay ", assaySlot, ": (detecting...)")
        }
    })
}

# Keep QC image rendering top-level so a later extraction wave can move it
# without reopening the module server wrapper's static QC output bindings.
renderMetabNormQcImageForAssay <- function(
    assaySlot,
    plotType,
    stagePrefix,
    normData,
    qcDir,
    renderImageFn = shiny::renderImage,
    fileExistsFn = file.exists,
    filePathFn = file.path,
    sanitizeAssayNameFn = function(assayName) {
        gsub("[^A-Za-z0-9]", "_", tolower(assayName))
    }
) {
    renderImageFn({
        normData$plot_refresh_trigger
        assayNames <- normData$assay_names

        if (is.null(assayNames) || length(assayNames) < assaySlot) {
            return(list(src = "", alt = "Assay not detected yet"))
        }

        assayName <- assayNames[[assaySlot]]
        safeName <- sanitizeAssayNameFn(assayName)
        filename <- paste0(safeName, "_", stagePrefix, "_", plotType, ".png")

        if (is.null(qcDir)) {
            return(list(src = "", alt = "QC directory not configured"))
        }

        imgPath <- filePathFn(qcDir, filename)

        if (fileExistsFn(imgPath)) {
            list(
                src = imgPath
                , contentType = "image/png"
                , width = "100%"
                , height = "auto"
                , alt = paste(plotType, "-", assayName)
            )
        } else {
            list(src = "", alt = paste("Plot not generated yet:", filename))
        }
    }, deleteFile = FALSE)
}

# Keep normalization-log mutation top-level so a later extraction wave can move
# it without reopening the module server wrapper's log state shell.
appendMetabNormNormalizationLog <- function(
    normData,
    message,
    timestampFn = function() {
        format(Sys.time(), "%H:%M:%S")
    }
) {
    timestamp <- timestampFn()
    updatedLog <- c(
        normData$normalization_log
        , sprintf("[%s] %s", timestamp, message)
    )
    normData$normalization_log <- updatedLog

    invisible(updatedLog)
}

# Keep normalization-log rendering top-level so a later extraction wave can move
# it without reopening the module server wrapper's output binding shell.
renderMetabNormNormalizationLog <- function(
    normData,
    renderTextFn = shiny::renderText,
    emptyMessage = "Normalization log will appear here as you apply steps..."
) {
    renderTextFn({
        normalizationLog <- normData$normalization_log

        if (length(normalizationLog) == 0) {
            return(emptyMessage)
        }

        paste(normalizationLog, collapse = "\n")
    })
}

# Keep correlation-filter summary rendering top-level so a later extraction
# wave can move it without reopening the module server wrapper's output shell.
renderMetabNormCorrelationFilterSummary <- function(
    normData,
    renderTextFn = shiny::renderText,
    buildSummaryFn = buildMetabNormCorrelationFilterSummary,
    incompleteMessage = "Apply correlation filter to see results..."
) {
    renderTextFn({
        if (!normData$correlation_filtering_complete) {
            return(incompleteMessage)
        }

        corrResults <- normData$correlation_results
        filteredObject <- normData$correlation_filtered_obj
        originalObject <- if (!is.null(normData$ruv_corrected_obj)) {
            normData$ruv_corrected_obj
        } else {
            normData$post_norm_obj
        }

        buildSummaryFn(
            corrResults = corrResults
            , filteredObject = filteredObject
            , originalObject = originalObject
        )
    })
}

# Keep final-QC rendering top-level so a later extraction wave can move it
# without reopening the module server wrapper's output shell.
renderMetabNormFinalQcPlot <- function(
    normData,
    colorVariableFn = function() NULL,
    shapeVariableFn = function() NULL,
    renderPlotFn = shiny::renderPlot,
    reqFn = shiny::req,
    resolveRenderStateFn = resolveMetabNormFinalQcRenderState,
    getPlotAestheticsFn = getPlotAesthetics,
    buildPcaPlotFn = buildMetabNormFinalQcPcaPlot
) {
    renderPlotFn({
        reqFn(normData$correlation_filtering_complete || normData$ruv_complete)

        finalQcState <- resolveRenderStateFn(
            correlationFilteredObject = normData$correlation_filtered_obj
            , ruvCorrectedObject = normData$ruv_corrected_obj
            , postNormObject = normData$post_norm_obj
        )

        if (isTRUE(finalQcState$isFallback)) {
            return(finalQcState$plot)
        }

        aesthetics <- getPlotAestheticsFn(
            colorVariable = colorVariableFn()
            , shapeVariable = shapeVariableFn()
        )

        buildPcaPlotFn(
            sourceObject = finalQcState$sourceObject
            , colorVar = aesthetics$color_var
            , shapeVar = aesthetics$shape_var
        )
    })
}

