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
        state <- initialState()
        reqFn(state)

        if (hasBatchColumnFn(state$data_cln)) {
            logInfoFn(
                "Design Matrix: 'Batch' column detected. Auto-assigning to batch column."
            )

            batchAssignments <- buildBatchAssignmentsFn(state$data_cln)
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
        dataClnReactive(state$data_cln)
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

    bindEventFn(observer, dataTbl(), input$reset_changes)
}

createProtDesignBuilderState <- function(
    reactiveValFn = shiny::reactiveVal
) {
    list(
        resultRv = reactiveValFn(NULL),
        designMatrix = reactiveValFn(NULL),
        dataClnReactive = reactiveValFn(NULL),
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
        reqFn(dataTbl())
        reqFn(configList())
        buildInitialState(
            dataTbl = dataTbl(),
            configList = configList(),
            columnMapping = columnMapping()
        )
    })
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

