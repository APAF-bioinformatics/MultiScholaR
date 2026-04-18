runLipidNormModuleServerShell <- function(
    input,
    output,
    session,
    id,
    workflowData,
    experimentPaths,
    omicType,
    experimentLabel,
    selectedTab = NULL,
    logInfoFn = logger::log_info,
    createReactiveStateFn = createLipidNormReactiveState,
    createStartupRuntimeFn = createLipidNormStartupRuntime,
    registerServerRuntimeFn = registerLipidNormServerRuntime
) {
    ns <- session$ns

    logInfoFn("=== METABOLOMICS NORMALIZATION MODULE STARTED ===")
    logInfoFn(paste("Module ID:", id))

    normData <- createReactiveStateFn()

    startupRuntime <- createStartupRuntimeFn(
        input = input
        , experimentPaths = experimentPaths
        , normData = normData
        , ns = ns
    )

    addLog <- startupRuntime$addLog
    getPlotAesthetics <- startupRuntime$getPlotAesthetics
    generateCompositeFromFiles <- startupRuntime$generateCompositeFromFiles

    registerServerRuntimeFn(
        input = input
        , output = output
        , session = session
        , workflowData = workflowData
        , experimentPaths = experimentPaths
        , omicType = omicType
        , experimentLabel = experimentLabel
        , normData = normData
        , startupRuntime = startupRuntime
        , addLog = addLog
        , getPlotAestheticsFn = getPlotAesthetics
        , generateCompositeFromFilesFn = generateCompositeFromFiles
        , selectedTab = selectedTab
    )

    invisible(normData)
}

runLipidNormModuleServerEntryShell <- function(
    id,
    workflowData,
    experimentPaths,
    omicType,
    experimentLabel,
    selectedTab = NULL,
    moduleServerFn = shiny::moduleServer,
    runModuleServerShellFn = runLipidNormModuleServerShell
) {
    moduleServerFn(id, function(input, output, session) {
        runModuleServerShellFn(
            input = input
            , output = output
            , session = session
            , id = id
            , workflowData = workflowData
            , experimentPaths = experimentPaths
            , omicType = omicType
            , experimentLabel = experimentLabel
            , selectedTab = selectedTab
        )
    })
}

runLipidNormModuleServerPublicWrapper <- function(
    id,
    workflow_data,
    experiment_paths,
    omic_type,
    experiment_label,
    selected_tab = NULL,
    runModuleServerEntryShellFn = runLipidNormModuleServerEntryShell
) {
    runModuleServerEntryShellFn(
        id = id
        , workflowData = workflow_data
        , experimentPaths = experiment_paths
        , omicType = omic_type
        , experimentLabel = experiment_label
        , selectedTab = selected_tab
    )
}

