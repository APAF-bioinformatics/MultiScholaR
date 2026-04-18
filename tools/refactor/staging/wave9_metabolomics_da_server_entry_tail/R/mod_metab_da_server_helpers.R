registerMetabDaServerBodyOutputs <- function(
    input,
    output,
    session,
    workflowData,
    experimentPaths,
    daData,
    registerOverviewOutputs = registerMetabDaOverviewOutputs,
    registerMainOutputsAndObservers = registerMetabDaMainOutputsAndObservers
) {
    registerOverviewOutputs(
        output = output,
        workflowData = workflowData,
        daData = daData
    )

    mainRegistrations <- registerMainOutputsAndObservers(
        input = input,
        output = output,
        workflowData = workflowData,
        daData = daData,
        session = session,
        experimentPaths = experimentPaths
    )

    invisible(mainRegistrations)
}

startMetabDaServerBody <- function(
    input,
    output,
    session,
    workflowData,
    experimentPaths,
    createServerState = createMetabDaServerState,
    registerServerBodyOutputs = registerMetabDaServerBodyOutputs
) {
    daData <- createServerState()

    registerServerBodyOutputs(
        input = input,
        output = output,
        session = session,
        workflowData = workflowData,
        experimentPaths = experimentPaths,
        daData = daData
    )
}

initializeMetabDaServerBody <- function(
    input,
    output,
    session,
    workflowData,
    experimentPaths,
    createServerState = createMetabDaServerState,
    registerServerBodyOutputs = registerMetabDaServerBodyOutputs,
    startServerBody = startMetabDaServerBody
) {
    startServerBody(
        input = input,
        output = output,
        session = session,
        workflowData = workflowData,
        experimentPaths = experimentPaths,
        createServerState = createServerState,
        registerServerBodyOutputs = registerServerBodyOutputs
    )
}

createMetabDaServerState <- function(
    createReactiveValues = shiny::reactiveValues
) {
    createReactiveValues(
        da_results_list = NULL,
        contrasts_available = NULL,
        assays_available = NULL,
        analysis_complete = FALSE,
        current_s4_object = NULL,
        formula_from_s4 = NULL,
        current_row_clusters = NULL,
        current_col_clusters = NULL,
        current_heatmap_plot = NULL
    )
}

runMetabDaServerModule <- function(
    id,
    workflowData,
    experimentPaths,
    moduleServer = shiny::moduleServer,
    runServerCallback = runMetabDaServerModuleCallback,
    initializeServerBody = initializeMetabDaServerBody,
    runServerModuleShell = runMetabDaServerModuleShell
) {
    runServerModuleShell(
        id = id,
        workflowData = workflowData,
        experimentPaths = experimentPaths,
        moduleServer = moduleServer,
        runServerCallback = runServerCallback,
        initializeServerBody = initializeServerBody
    )
}

runMetabDaServerModuleShell <- function(
    id,
    workflowData,
    experimentPaths,
    moduleServer = shiny::moduleServer,
    runServerCallback = runMetabDaServerModuleCallback,
    initializeServerBody = initializeMetabDaServerBody,
    createServerModuleHandler = createMetabDaServerModuleHandler
) {
    serverModuleHandler <- createServerModuleHandler(
        workflowData = workflowData,
        experimentPaths = experimentPaths,
        runServerCallback = runServerCallback,
        initializeServerBody = initializeServerBody
    )

    moduleServer(id, serverModuleHandler)
}

createMetabDaServerModuleHandler <- function(
    workflowData,
    experimentPaths,
    runServerCallback = runMetabDaServerModuleCallback,
    initializeServerBody = initializeMetabDaServerBody
) {
    function(input, output, session) {
        runServerCallback(
            input = input,
            output = output,
            session = session,
            workflowData = workflowData,
            experimentPaths = experimentPaths,
            initializeServerBody = initializeServerBody
        )
    }
}

runMetabDaServerModuleCallback <- function(
    input,
    output,
    session,
    workflowData,
    experimentPaths,
    initializeServerBody = initializeMetabDaServerBody
) {
    initializeServerBody(
        input = input,
        output = output,
        session = session,
        workflowData = workflowData,
        experimentPaths = experimentPaths
    )
}

runMetabDaServerEntry <- function(
    id,
    workflowData,
    experimentPaths,
    omicType = NULL,
    experimentLabel = NULL,
    runServerModule = runMetabDaServerModule
) {
    runServerModule(
        id = id,
        workflowData = workflowData,
        experimentPaths = experimentPaths
    )
}

