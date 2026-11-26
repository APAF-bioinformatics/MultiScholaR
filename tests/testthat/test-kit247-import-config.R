# KIT-247 departure 3 (fix commit 0da49542): on proteomics Import Existing Design, the folder's config.ini is read into
# the bus before the design comparison; while a changed design waits for the user's answer, the session's own
# configuration is kept (Cancel then changes nothing), and the imported one applies only on Save. A first save or an
# unchanged design keeps the imported configuration, as before. Function level, the import frame's injected steps.

kit247ImportShell <- function(flowValue) {
    workflow <- new.env(parent = emptyenv())
    workflow$config_list <- list(source = "session")
    runProtDesignImportObserverShell(workflowData = workflow, experimentPaths = list(), importPath = "unused",
        taxonId = 9606L, organismName = "Homo sapiens", session = NULL,
        resolveImportArtifacts = function(importPath, selectedFastaPath) list(ok = TRUE, designFile = "d", dataClnFile = "c",
            contrastFile = "k"),
        loadImportedArtifacts = function(workflowData, experimentPaths, designFile, dataClnFile, contrastFile) {
            workflowData$config_list <- list(source = "imported")
            list()
        },
        runImportConfirmationFlow = function(...) flowValue,
        showNotification = function(...) invisible(NULL), removeNotification = function(...) invisible(NULL))
    workflow$config_list
}

test_that("KIT-247 D3: a changed design waiting for confirmation keeps the session's configuration", {
    expect_identical(kit247ImportShell("changed"), list(source = "session"))
})

test_that("KIT-247 D3 guard: a first or unchanged import keeps the imported configuration, as before", {
    expect_identical(kit247ImportShell(list(stateName = "protein_s4_initial")), list(source = "imported"))
})
