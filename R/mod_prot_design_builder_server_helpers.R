registerProtDesignInitialStateSyncObserver <- function(
    input,
    session,
    initialState,
    dataTbl,
    designMatrix,
    dataClnReactive,
    groups,
    factors,
    contrasts,
    removedSamples,
    observeFn = shiny::observe,
    bindEventFn = shiny::bindEvent,
    reqFn = shiny::req,
    hasBatchColumnFn = function(dataCln) "Batch" %in% names(dataCln),
    buildBatchAssignmentsFn = function(dataCln) {
        dataCln |>
            dplyr::distinct(Run, Batch)
    },
    mergeBatchAssignmentsFn = function(currentDesignMatrix, batchAssignments) {
        currentDesignMatrix |>
            dplyr::left_join(batchAssignments, by = "Run") |>
            dplyr::mutate(batch = dplyr::coalesce(batch, Batch)) |>
            dplyr::select(-Batch)
    },
    updateSelectizeInputFn = shiny::updateSelectizeInput,
    updateTextInputFn = shiny::updateTextInput,
    logInfoFn = logger::log_info,
    resultRv = NULL
) {
    observer <- observeFn({
        withDIAStateScope(getDIAActiveManager(), function() {
        restored <- NULL
        state <- initialState()
        reqFn(state)
        manager <- getDIAActiveManager()
        if (inherits(manager, "DIAStateAdapter")) {
            reference <- getDIAReactiveReference(dataTbl)
            pending <- manager$.context$design_adoption_pending
            restored <- manager$.context$peptide_design_restore
            if (!is.null(pending) || !is.null(restored)) {
                # A Design Save result has already been consumed. Drop its native
                # assay alias without re-saving Design or publishing another head.
                cache <- attr(resultRv, "dia_reference", exact = TRUE)
                if (is.function(cache)) cache(NULL)
            }
            if (!is.null(restored)) {
                if (!isSavedWorkflowReference(reference) ||
                    !identical(restored$source_reference, reference$entity_id)) {
                    stop("restored Design source differs")
                }
                for (field in c("design_matrix", "groups", "factors", "contrasts", "formula")) {
                    state[field] <- restored[field]
                }
                state$data_reference <- manager$.context$snapshot$gui_references$data_cln
                manager$.context$peptide_design_restore <- NULL
            }
            if (!is.null(pending) && isSavedWorkflowReference(reference) &&
                identical(pending, reference$entity_id)) {
                # Evaluate initialState above to replace its old resident cache,
                # but preserve the current design selections on owner adoption.
                dataClnReactive(manager$.context$snapshot$gui_references$data_cln)
                manager$.context$design_adoption_pending <- NULL
                return(invisible(NULL))
            }
        }

        batch_source <- state$sample_metadata %||% state$data_cln
        if (is.null(restored) && hasBatchColumnFn(batch_source)) {
            logInfoFn(
                "Design Matrix: 'Batch' column detected. Auto-assigning to batch column."
            )

            batchAssignments <- buildBatchAssignmentsFn(batch_source)
            state$design_matrix <- mergeBatchAssignmentsFn(
                state$design_matrix,
                batchAssignments
            )

            logInfoFn(
                "Design Matrix: Pre-populated 'batch' column with batch assignments."
            )
            logInfoFn(
                "Note: Batch is kept separate from experimental factors to enable cross-batch comparisons."
            )
        }

        designMatrix(state$design_matrix)
        dataClnReactive(state$data_reference %||% state$data_cln)
        groups(state$groups)
        factors(state$factors)
        contrasts(state$contrasts)
        removedSamples(if (is.null(restored)) character(0) else restored$removed_samples)

        sortedRuns <- state$design_matrix$Run
        updateSelectizeInputFn(
            session,
            "sample_to_rename",
            choices = sortedRuns,
            selected = ""
        )
        updateSelectizeInputFn(
            session,
            "selected_runs",
            choices = sortedRuns,
            selected = ""
        )
        updateSelectizeInputFn(
            session,
            "samples_to_transform",
            choices = sortedRuns,
            selected = ""
        )
        updateSelectizeInputFn(
            session,
            "tech_rep_samples",
            choices = sortedRuns,
            selected = ""
        )
        updateSelectizeInputFn(
            session,
            "samples_to_remove",
            choices = sortedRuns,
            selected = ""
        )
        updateTextInputFn(session, "formula_string", value = state$formula)
        })
    })

    bindEventFn(observer, getDIAReactiveReference(dataTbl), input$reset_changes)
}

createProtDesignBuilderState <- function(
    reactiveValFn = shiny::reactiveVal
) {
    list(
        resultRv = newDIASavedReactiveVal(NULL, "design-result", reactiveValFn),
        designMatrix = reactiveValFn(NULL),
        dataClnReactive = newDIASavedReactiveVal(NULL, "design-cleaned-data", reactiveValFn,
            snapshot_field = "design_current"),
        groups = reactiveValFn(NULL),
        factors = reactiveValFn(NULL),
        contrasts = reactiveValFn(NULL),
        removedSamples = reactiveValFn(character(0))
    )
}

createProtDesignInitialStateReactive <- function(
    dataTbl,
    configList,
    columnMapping,
    buildInitialState = buildProtDesignInitialState,
    reactiveFn = shiny::reactive,
    reqFn = shiny::req
) {
    reactiveFn({
        manager <- getDIAActiveManager()
        reference <- if (inherits(manager, "DIAStateAdapter")) getDIAReactiveReference(dataTbl) else NULL
        if (isSavedWorkflowReference(reference)) {
            reqFn(reference)
            reqFn(configList())
            return(withDIAStateScope(manager, function() {
                buildDIAReferenceDesignState(manager, reference, configList(), columnMapping())
            }))
        }
        reqFn(dataTbl())
        reqFn(configList())
        buildInitialState(
            dataTbl = dataTbl(),
            configList = configList(),
            columnMapping = columnMapping()
        )
    })
}

#' The peptide checkpoint's reader of the builder's current design (`peptide_checkpoint_design`,
#' read by `utils_workflow_dia_sessions.R`), cleared when the session ends. KIT-103 W6: the
#' `dia_long` adapter's start, after the initial-state observer, as the proteomics builder set it.
#' @param state The shared builder's reactive values.
#' @noRd
registerProtDesignCheckpointReader <- function(input, session, state) {
    if (!is.environment(session$userData)) return(invisible(NULL))
    metadata <- function() shiny::isolate(list(
        design_matrix = state$design_matrix(), groups = state$groups(),
        factors = state$factors(), contrasts = state$contrasts(),
        removed_samples = state$removed_samples(), formula = input$formula_string))
    session$userData$peptide_checkpoint_design <- metadata
    session$onSessionEnded(function() {
        if (identical(session$userData$peptide_checkpoint_design, metadata)) {
            session$userData$peptide_checkpoint_design <- NULL
        }
    })
    invisible(NULL)
}
