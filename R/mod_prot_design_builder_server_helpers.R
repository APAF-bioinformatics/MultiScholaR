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
    logInfoFn = logger::log_info
) {
    observer <- observeFn({
        withDIAStateScope(getDIAActiveManager(), function() {
        state <- initialState()
        reqFn(state)

        batch_source <- state$sample_metadata %||% state$data_cln
        if (hasBatchColumnFn(batch_source)) {
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
        removedSamples(character(0))

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

createProtDesignMutableStateShells <- function(
    createBuilderStateFn = createProtDesignBuilderState
) {
    builderState <- createBuilderStateFn()

    list(
        resultRv = builderState$resultRv,
        designMatrix = builderState$designMatrix,
        dataClnReactive = builderState$dataClnReactive,
        groups = builderState$groups,
        factors = builderState$factors,
        contrasts = builderState$contrasts,
        removedSamples = builderState$removedSamples
    )
}

createProtDesignDataTableProxy <- function(
    proxyId = "data_table",
    dataTableProxyFn = DT::dataTableProxy
) {
    dataTableProxyFn(proxyId)
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

#' Initialize the private DIA design from bounded metadata and immutable data.
#' @param manager Current session owner.
#' @param reference Imported-table reference.
#' @param config_list,column_mapping Existing configuration and mapped columns.
#' @return Initial design metadata with an explicit data reference, never an assay.
#' @noRd
buildDIAReferenceDesignState <- function(manager, reference, config_list, column_mapping) {
    context <- manager$.context
    validateSavedWorkflowReference(context$store, reference, payload = FALSE)
    metadata <- getDIAImportMetadata(manager)
    if (!identical(metadata$reference_id, reference$entity_id) ||
        !identical(metadata$column_mapping, column_mapping)) {
        if (usesDIAReferenceStates(manager) && isWorkflowTableReference(reference)) {
            cli::cli_abort("saved import metadata is missing or mismatched", class = "saved_entity_corrupt")
        }
        # Older validated native sessions have no import summary. Derive it once,
        # in this operation, rather than forcing future metadata events to read.
        table <- resolveDIAOwnedValue(reference, manager)
        metadata <- recordDIAImportMetadata(manager, table,
            list(data = table, data_type = "peptide", column_mapping = column_mapping), reference)
        table <- NULL
    }
    initial <- context$snapshot$design_initial
    if (!identical(initial$source_reference, reference) ||
        !identical(initial$column_mapping, column_mapping)) {
        cleaned_reference <- reference
        if (length(metadata$needs_numeric)) {
            if (usesDIAReferenceStates(manager)) {
                cleaned_reference <- transformDIADesignTable(manager, reference,
                    numeric_columns = metadata$needs_numeric, strip_parser = TRUE)
            } else {
                converted <- buildProtDesignInitialState(resolveDIAOwnedValue(reference, manager),
                    config_list, column_mapping)
                cleaned_reference <- saveDIAOwnedValue(converted$data_cln, "design-numeric", manager)
                converted <- NULL
            }
        }
        initial <- list(source_reference = reference, column_mapping = column_mapping,
            data_reference = cleaned_reference)
        snapshot <- context$snapshot
        snapshot$design_initial <- initial
        publishDIAStateHead(context, snapshot)
        releaseDIAReplacedOwner(manager, reference, cleaned_reference)
    }
    state <- buildProtDesignInitialState(metadata$samples, config_list, list())
    state["data_cln"] <- list(NULL)
    state$data_reference <- initial$data_reference
    state$sample_metadata <- metadata$samples
    state
}

initializeProtDesignBuilderServerState <- function(
    input,
    session,
    dataTbl,
    configList,
    columnMapping,
    buildInitialState = buildProtDesignInitialState,
    createInitialStateReactive = createProtDesignInitialStateReactive,
    createMutableStateShells = createProtDesignMutableStateShells,
    registerInitialStateSyncObserver = registerProtDesignInitialStateSyncObserver
) {
    initialState <- createInitialStateReactive(
        dataTbl = dataTbl,
        configList = configList,
        columnMapping = columnMapping,
        buildInitialState = buildInitialState
    )

    mutableState <- createMutableStateShells()

    registerInitialStateSyncObserver(
        input = input,
        session = session,
        initialState = initialState,
        dataTbl = dataTbl,
        designMatrix = mutableState$designMatrix,
        dataClnReactive = mutableState$dataClnReactive,
        groups = mutableState$groups,
        factors = mutableState$factors,
        contrasts = mutableState$contrasts,
        removedSamples = mutableState$removedSamples
    )

    c(
        list(initialState = initialState),
        mutableState
    )
}

registerProtDesignRenderOutputShells <- function(
    input,
    output,
    session,
    proxyDataTable,
    designMatrix,
    removedSamples,
    factors,
    contrasts,
    registerDataTableOutput = registerProtDesignDataTableOutput,
    registerDataTableProxyRefreshObserver = registerProtDesignDataTableProxyRefreshObserver,
    registerAvailableFactorsDisplayOutput = registerProtDesignAvailableFactorsDisplayOutput,
    registerDefinedContrastsDisplayOutput = registerProtDesignDefinedContrastsDisplayOutput,
    registerRangePreviewOutput = registerProtDesignRangePreviewOutput,
    registerTechRepSummaryOutput = registerProtDesignTechRepSummaryOutput,
    registerRemovedSamplesDisplayOutput = registerProtDesignRemovedSamplesDisplayOutput,
    registerReplicateInputsOutput = registerProtDesignReplicateInputsOutput,
    registerContrastFactorsInfoOutput = registerProtDesignContrastFactorsInfoOutput
) {
    registerDataTableOutput(
        output = output,
        designMatrix = designMatrix,
        removedSamples = removedSamples
    )

    registerDataTableProxyRefreshObserver(
        proxyDataTable = proxyDataTable,
        designMatrix = designMatrix,
        removedSamples = removedSamples
    )

    registerAvailableFactorsDisplayOutput(
        output = output,
        factors = factors
    )

    registerDefinedContrastsDisplayOutput(
        output = output,
        contrasts = contrasts,
        formulaString = function() input$formula_string
    )

    registerRangePreviewOutput(
        output = output,
        input = input
    )

    registerTechRepSummaryOutput(
        output = output,
        designMatrix = designMatrix
    )

    registerRemovedSamplesDisplayOutput(
        output = output,
        removedSamples = removedSamples
    )

    registerReplicateInputsOutput(
        output = output,
        input = input,
        nsFn = session$ns
    )

    registerContrastFactorsInfoOutput(
        output = output,
        formulaString = function() input$formula_string
    )

    invisible(NULL)
}

registerProtDesignRenameObserverShells <- function(
    input,
    designMatrix,
    dataCln,
    session,
    registerRenameSampleObserver = registerProtDesignRenameSampleObserver,
    registerBulkRenameObserver = registerProtDesignBulkRenameObserver
) {
    registerRenameSampleObserver(
        input = input,
        designMatrix = designMatrix,
        dataCln = dataCln,
        session = session
    )

    registerBulkRenameObserver(
        input = input,
        designMatrix = designMatrix,
        dataCln = dataCln
    )

    invisible(NULL)
}

registerProtDesignFactorMetadataObserverShells <- function(
    input,
    factors,
    designMatrix,
    groups,
    session,
    registerAddFactorObserver = registerProtDesignAddFactorObserver,
    registerAssignMetadataObserver = registerProtDesignAssignMetadataObserver
) {
    registerAddFactorObserver(
        input = input,
        factors = factors,
        session = session
    )

    registerAssignMetadataObserver(
        input = input,
        designMatrix = designMatrix,
        groups = groups
    )

    invisible(NULL)
}

registerProtDesignActionObserverShells <- function(
    input,
    designMatrix,
    contrasts,
    removedSamples,
    session,
    registerAssignTechRepsObserver = registerProtDesignAssignTechRepsObserver,
    registerAddContrastObserver = registerProtDesignAddContrastObserver,
    registerRemoveSamplesObserver = registerProtDesignRemoveSamplesObserver
) {
    registerAssignTechRepsObserver(
        input = input,
        designMatrix = designMatrix
    )

    registerAddContrastObserver(
        input = input,
        contrasts = contrasts
    )

    registerRemoveSamplesObserver(
        input = input,
        removedSamples = removedSamples,
        session = session
    )

    invisible(NULL)
}

registerProtDesignResetAndSaveObserverShells <- function(
    input,
    session,
    initialState,
    designMatrix,
    dataClnReactive,
    removedSamples,
    factors,
    groups,
    contrasts,
    configList,
    resultRv,
    registerResetRequestObserver = registerProtDesignResetRequestObserver,
    registerResetConfirmationObserver = registerProtDesignResetConfirmationObserver,
    registerSaveResultsObserver = registerProtDesignSaveResultsObserver
) {
    registerResetRequestObserver(
        input = input,
        session = session
    )

    registerResetConfirmationObserver(
        input = input,
        initialState = initialState,
        designMatrix = designMatrix,
        dataClnReactive = dataClnReactive,
        removedSamples = removedSamples,
        factors = factors,
        groups = groups,
        contrasts = contrasts,
        session = session
    )

    registerSaveResultsObserver(
        input = input,
        designMatrix = designMatrix,
        removedSamples = removedSamples,
        dataCln = dataClnReactive,
        contrastData = contrasts,
        configList = configList,
        resultRv = resultRv
    )

    invisible(NULL)
}

registerProtDesignEventObserverShells <- function(
    input,
    session,
    initialState,
    designMatrix,
    dataClnReactive,
    removedSamples,
    factors,
    groups,
    contrasts,
    configList,
    resultRv,
    registerRenameObserverShells = registerProtDesignRenameObserverShells,
    registerFactorMetadataObserverShells = registerProtDesignFactorMetadataObserverShells,
    registerActionObserverShells = registerProtDesignActionObserverShells,
    registerResetAndSaveObserverShells = registerProtDesignResetAndSaveObserverShells
) {
    registerRenameObserverShells(
        input = input,
        designMatrix = designMatrix,
        dataCln = dataClnReactive,
        session = session
    )

    registerFactorMetadataObserverShells(
        input = input,
        factors = factors,
        designMatrix = designMatrix,
        groups = groups,
        session = session
    )

    registerActionObserverShells(
        input = input,
        designMatrix = designMatrix,
        contrasts = contrasts,
        removedSamples = removedSamples,
        session = session
    )

    registerResetAndSaveObserverShells(
        input = input,
        session = session,
        initialState = initialState,
        designMatrix = designMatrix,
        dataClnReactive = dataClnReactive,
        removedSamples = removedSamples,
        factors = factors,
        groups = groups,
        contrasts = contrasts,
        configList = configList,
        resultRv = resultRv
    )

    invisible(NULL)
}

registerProtDesignBuilderServerShells <- function(
    input,
    output,
    session,
    initialState,
    proxyDataTable,
    designMatrix,
    dataClnReactive,
    removedSamples,
    factors,
    groups,
    contrasts,
    configList,
    resultRv,
    registerInputSyncObserverShells = registerProtDesignInputSyncObserverShells,
    registerRenderOutputShells = registerProtDesignRenderOutputShells,
    registerEventObserverShells = registerProtDesignEventObserverShells
) {
    registerInputSyncObserverShells(
        input = input,
        session = session,
        designMatrix = designMatrix,
        removedSamples = removedSamples,
        factors = factors,
        groups = groups
    )

    registerRenderOutputShells(
        input = input,
        output = output,
        session = session,
        proxyDataTable = proxyDataTable,
        designMatrix = designMatrix,
        removedSamples = removedSamples,
        factors = factors,
        contrasts = contrasts
    )

    registerEventObserverShells(
        input = input,
        session = session,
        initialState = initialState,
        designMatrix = designMatrix,
        dataClnReactive = dataClnReactive,
        removedSamples = removedSamples,
        factors = factors,
        groups = groups,
        contrasts = contrasts,
        configList = configList,
        resultRv = resultRv
    )

    invisible(NULL)
}

completeProtDesignBuilderServerBootstrap <- function(
    input,
    output,
    session,
    builderState,
    configList,
    createDataTableProxy = createProtDesignDataTableProxy,
    registerBuilderServerShells = registerProtDesignBuilderServerShells
) {
    proxyDataTable <- createDataTableProxy()

    registerBuilderServerShells(
        input = input,
        output = output,
        session = session,
        initialState = builderState$initialState,
        proxyDataTable = proxyDataTable,
        designMatrix = builderState$designMatrix,
        dataClnReactive = builderState$dataClnReactive,
        removedSamples = builderState$removedSamples,
        factors = builderState$factors,
        groups = builderState$groups,
        contrasts = builderState$contrasts,
        configList = configList,
        resultRv = builderState$resultRv
    )

    builderState$resultRv
}

runProtDesignBuilderModuleServerShell <- function(
    input,
    output,
    session,
    dataTbl,
    configList,
    columnMapping,
    buildInitialState = buildProtDesignInitialState,
    createInitialStateReactive = createProtDesignInitialStateReactive,
    createMutableStateShells = createProtDesignMutableStateShells,
    registerInitialStateSyncObserver = registerProtDesignInitialStateSyncObserver,
    initializeBuilderServerState = initializeProtDesignBuilderServerState,
    completeBuilderServerBootstrap = completeProtDesignBuilderServerBootstrap,
    createDataTableProxy = createProtDesignDataTableProxy,
    registerBuilderServerShells = registerProtDesignBuilderServerShells
) {
    builderState <- initializeBuilderServerState(
        input = input,
        session = session,
        dataTbl = dataTbl,
        configList = configList,
        columnMapping = columnMapping,
        buildInitialState = buildInitialState,
        createInitialStateReactive = createInitialStateReactive,
        createMutableStateShells = createMutableStateShells,
        registerInitialStateSyncObserver = registerInitialStateSyncObserver
    )

    completeBuilderServerBootstrap(
        input = input,
        output = output,
        session = session,
        builderState = builderState,
        configList = configList,
        createDataTableProxy = createDataTableProxy,
        registerBuilderServerShells = registerBuilderServerShells
    )
}

runProtDesignBuilderServerEntryShell <- function(
    id,
    dataTbl,
    configList,
    columnMapping,
    moduleServer = shiny::moduleServer,
    runBuilderModuleServerShell = runProtDesignBuilderModuleServerShell
) {
    moduleServer(id, function(input, output, session) {
        runBuilderModuleServerShell(
            input = input,
            output = output,
            session = session,
            dataTbl = dataTbl,
            configList = configList,
            columnMapping = columnMapping
        )
    })
}
