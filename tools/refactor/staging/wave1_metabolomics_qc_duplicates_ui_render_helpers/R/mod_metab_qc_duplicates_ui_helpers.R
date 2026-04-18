buildMetabDuplicateSummaryUi <- function(
    dupList
    , infoText = " Click 'Detect Duplicates' to scan for duplicate features."
    , paragraphFn = shiny::p
    , iconFn = shiny::icon
    , listTagFn = shiny::tags$ul
    , itemTagFn = shiny::tags$li
) {
    if (is.null(dupList)) {
        return(paragraphFn(
            iconFn("info-circle")
            , infoText
            , style = "color: #666;"
        ))
    }

    summaryItems <- lapply(names(dupList), function(assayName) {
        dupDf <- dupList[[assayName]]
        duplicateCount <- if (is.null(dupDf)) 0 else nrow(dupDf)

        iconClass <- if (duplicateCount == 0) "check-circle" else "exclamation-triangle"
        iconColor <- if (duplicateCount == 0) "green" else "orange"

        itemTagFn(
            iconFn(iconClass, style = paste("color:", iconColor))
            , sprintf(" %s: %d duplicate IDs", assayName, duplicateCount)
        )
    })

    listTagFn(summaryItems, style = "list-style: none; padding-left: 0;")
}

buildMetabDuplicateTablesUi <- function(
    dupList
    , nsFn = identity
    , wellPanelFn = shiny::wellPanel
    , iconFn = shiny::icon
    , headerFn = shiny::h5
    , paragraphFn = shiny::p
    , tabPanelFn = shiny::tabPanel
    , breakFn = shiny::br
    , tableOutputFn = DT::DTOutput
    , tabsetPanelFn = shiny::tabsetPanel
) {
    if (is.null(dupList)) {
        return(NULL)
    }

    assaysWithDups <- names(dupList)[vapply(dupList, function(dupDf) {
        !is.null(dupDf) && nrow(dupDf) > 0
    }, logical(1))]

    if (length(assaysWithDups) == 0) {
        return(wellPanelFn(
            iconFn("check-circle", style = "color: green; font-size: 24px;")
            , headerFn("No duplicates found in any assay!")
            , paragraphFn("All metabolite IDs are unique. No resolution needed.")
        ))
    }

    tabList <- lapply(assaysWithDups, function(assayName) {
        tabPanelFn(
            assayName
            , breakFn()
            , tableOutputFn(
                nsFn(paste0("dup_table_", gsub("[^a-zA-Z0-9]", "_", assayName)))
            )
        )
    })

    do.call(tabsetPanelFn, c(list(id = nsFn("dup_tables_tabs")), tabList))
}

