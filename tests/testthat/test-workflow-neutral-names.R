# KIT-074: the kit functions the omics code calls (or a later stage-3b ticket names) carry omic-neutral names; each
# DIA name stays as the alias the proteomics modules, the harness and the tests keep. The rows are the inventory
# `docs/qa/workflow-execution/kit/kit-074/inventory.tsv`: its fourteen new names and the five KIT-045 names the
# omics call sites moved to.

test_that("every KIT-074 inventory row resolves under its neutral name and an identical DIA alias", {
    ns <- asNamespace("MultiScholaR")
    aliases <- c(saveDIAStateValue = "saveWorkflowStateValue", recordDIAWorkflowInputs = "recordWorkflowInputs",
        clearDIAStepProvenance = "clearWorkflowStepProvenance", collectDIANativeReferences = "collectWorkflowNativeReferences",
        runDIAOwnedSessionExport = "runWorkflowOwnedSessionExport", createDIASessionOwner = "createWorkflowSessionOwner",
        readDIAStateGuiMetadata = "readWorkflowStateGuiMetadata", getDIAOwnedPayload = "getWorkflowOwnedPayload",
        setDIAOwnedPayload = "setWorkflowOwnedPayload", getDIAOwnedInput = "getWorkflowOwnedInput",
        resolveDIAOwnedValue = "resolveWorkflowOwnedValue", getDIAReactiveReference = "getWorkflowReactiveReference",
        assertDIAExportReady = "assertWorkflowExportReady", publishDIAActionOutputs = "publishWorkflowActionOutputs",
        saveDIAOwnedValue = "saveWorkflowOwnedValue", publishDIAStateHead = "publishWorkflowStateHead",
        withDIAStateScope = "withWorkflowStateScope", selectDIAStateReference = "selectWorkflowStateReference",
        closeDIAStateAdapter = "closeWorkflowStateOwner")
    expect_length(aliases, 19L)
    for (old in names(aliases)) {
        expect_true(is.function(get0(aliases[[old]], ns, inherits = FALSE)), info = aliases[[old]])
        expect_identical(get(old, ns), get(aliases[[old]], ns), info = old)
    }
})
