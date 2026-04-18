registerMetabDuplicateTableRenderers <- function(
    dupList
    , output
    , outputIdBuilderFn = function(assayName) {
        paste0("dup_table_", gsub("[^a-zA-Z0-9]", "_", assayName))
    }
    , renderDtFn = DT::renderDT
    , datatableFn = DT::datatable
    , datatableOptions = list(
        pageLength = 10
        , scrollX = TRUE
        , dom = 'frtip'
    )
) {
    if (is.null(dupList)) {
        return(invisible(character()))
    }

    registeredOutputIds <- unlist(lapply(names(dupList), function(assayName) {
        dupDf <- dupList[[assayName]]
        if (is.null(dupDf) || nrow(dupDf) == 0) {
            return(NULL)
        }

        outputId <- outputIdBuilderFn(assayName)
        output[[outputId]] <- renderDtFn({
            datatableFn(
                dupDf
                , options = datatableOptions
                , rownames = FALSE
                , class = "compact stripe"
            )
        })

        outputId
    }), use.names = FALSE)

    invisible(registeredOutputIds)
}

renderMetabDuplicateFilterPlot <- function(
    filterPlot
    , reqFn = shiny::req
    , inheritsFn = inherits
    , gridDrawFn = grid::grid.draw
    , printFn = print
) {
    plotObject <- filterPlot()
    reqFn(plotObject)

    if (inheritsFn(plotObject, "grob") || inheritsFn(plotObject, "gtable")) {
        gridDrawFn(plotObject)
    } else if (inheritsFn(plotObject, "ggplot")) {
        printFn(plotObject)
    }

    invisible(NULL)
}

