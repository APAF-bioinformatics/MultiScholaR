test_that("adoption clears consumed Design results without changing selections or the head", {
    skip("quarantined: Q-22 (KIT-048) the legacy-mode manager path these tests drove is retired; the owner is the one path")
    recipe <- readWorkflowCheckpointRecipe("small")
    recipe$n_proteins <- 8L; recipe$n_samples <- 8L
    recipe$peptides_per_protein <- 4L; recipe$forms_per_peptide <- 2L
    object <- makeWorkflowCheckpointPeptide(recipe)
    root <- withr::local_tempdir()
    native <- WorkflowState$new(audit_enabled = FALSE)
    native$saveState("intensity_filtered", object, list(), "input")
    workflow <- shiny::reactiveValues(state_manager = native, data_tbl = object@peptide_data,
        data_cln = object@peptide_data, config_list = list(globalParameters = list(workflow_type = "DIA")),
        data_type = "peptide", column_mapping = list(run_col = "Run", protein_col = "Protein.Group",
            peptide_col = "Stripped.Sequence", raw_quantity_col = "Precursor.Quantity",
            norm_quantity_col = "Precursor.Normalised"))
    shiny::testServer(function(input, output, session) {}, {
        builder <- createProtDesignBuilderState()
        runProtDesignSaveResultsObserverShell(object@design_matrix, dataCln = object@peptide_data,
            contrastData = NULL, configList = workflow$config_list, formulaString = "~ 0 + group",
            resultSetter = builder$resultRv, showNotification = function(...) NULL)
        expect_true(is.data.frame(getDIAReactiveReference(builder$resultRv)$data_cln))
        manager <- stageDIAGUIReferenceOwner(native, list(source_dir = root))
        adoptDIAGUIReferenceOwner(manager, workflow, session)
        head <- manager$.context$head
        preserve <- function(...) stop("unexpected Design reset")
        observer <- registerProtDesignInitialStateSyncObserver(input, session,
            initialState = function() buildDIAReferenceDesignState(manager, workflow$data_tbl,
                workflow$config_list, workflow$column_mapping), dataTbl = function() workflow$data_tbl,
            designMatrix = preserve, dataClnReactive = builder$dataClnReactive,
            groups = preserve, factors = preserve, contrasts = preserve, removedSamples = preserve,
            resultRv = builder$resultRv)
        session$flushReact()
        expect_null(getDIAReactiveReference(builder$resultRv))
        expect_true(isSavedWorkflowReference(getDIAReactiveReference(builder$dataClnReactive)))
        expect_identical(manager$.context$head, head)
        expect_identical(diaStateOwnership(manager)$full_s4_reads, 0L)
        observer$destroy()
        saved_design <- object@design_matrix
        saved_design$batch <- NA_character_
        manager$.context$peptide_design_restore <- list(
            source_reference = workflow$data_tbl$entity_id, design_matrix = saved_design,
            groups = c("control", "treated", "unused"), factors = character(),
            contrasts = data.frame(), formula = "~ 0 + group", removed_samples = "removed")
        restored_observer <- registerProtDesignInitialStateSyncObserver(input, session,
            initialState = function() list(sample_metadata = data.frame(Run = saved_design$Run,
                Batch = "import-default")), dataTbl = function() workflow$data_tbl,
            designMatrix = builder$designMatrix, dataClnReactive = builder$dataClnReactive,
            groups = builder$groups, factors = builder$factors, contrasts = builder$contrasts,
            removedSamples = builder$removedSamples, resultRv = builder$resultRv)
        session$flushReact()
        expect_identical(builder$designMatrix(), saved_design)
        expect_identical(builder$groups(), c("control", "treated", "unused"))
        expect_identical(builder$removedSamples(), "removed")
        clean <- function() diaStateOwnership(manager)$full_s4_reads == 0L &&
            !is.data.frame(getDIAReactiveReference(builder$resultRv)$data_cln)
        expect_true(clean())
        backing <- attr(builder$resultRv, "dia_reference")
        backing(list(data_cln = object@peptide_data))
        expect_false(clean())
        backing(NULL)
        manager$getState("intensity_filtered")
        expect_false(clean())
        restored_observer$destroy()
        closeDIAStateAdapter(manager)
    })
})

test_that("peptide checkpoints preserve resident history, relocation and failed-save recovery", {
    recipe <- readWorkflowCheckpointRecipe("small")
    recipe$n_proteins <- 8L; recipe$n_samples <- 8L
    recipe$peptides_per_protein <- 4L; recipe$forms_per_peptide <- 2L
    object <- makeWorkflowCheckpointPeptide(recipe)
    root <- withr::local_tempdir()
    paths <- list(source_dir = file.path(root, "source"), peptide_qc_dir = file.path(root, "peptide"))
    lapply(paths, dir.create)
    writeLines("accepted output", file.path(paths$peptide_qc_dir, "peptide.tsv"))
    native <- WorkflowState$new(audit_enabled = FALSE)
    native$setWorkflowType("DIA")
    native$saveState("intensity_filtered", object, list(), "before")
    changed <- object; changed@peptide_data$Precursor.Quantity[[1L]] <- 42
    native$saveState("sample_filtered", changed, list(), "after")
    workflow <- shiny::reactiveValues(state_manager = native, data_tbl = object@peptide_data,
        data_cln = object@peptide_data, data_type = "peptide", data_format = "diann",
        config_list = list(globalParameters = list(workflow_type = "DIA"),
            rollUpPrecursorToPeptide = list(core_utilisation = structure(
                rep(list(methods::new("externalptr")), 8L), class = "multidplyr_cluster"))),
        design_matrix = object@design_matrix,
        column_mapping = list(run_col = "Run", protein_col = "Protein.Group",
            peptide_col = "Stripped.Sequence", raw_quantity_col = "Precursor.Quantity",
            norm_quantity_col = "Precursor.Normalised"))
    original <- digest::digest(native$states)
    shiny::testServer(function(input, output, session) {}, {
        session$userData$peptide_checkpoint_design <- function() list(
            design_matrix = object@design_matrix, groups = unique(object@design_matrix$group),
            factors = character(), contrasts = data.frame(), removed_samples = character(),
            formula = "~ 0 + group")
        publication <- saveDIAPeptideCheckpoint(workflow, session, paths)
        expect_true(file.exists(publication$archive_path))
        expect_identical(workflow$state_manager, native)
        expect_identical(digest::digest(native$states), original)
        hash <- digest::digest(file = publication$archive_path)
        with_mocked_bindings({
            expect_error(saveDIAPeptideCheckpoint(workflow, session, paths), "injected save failure")
        }, workflowExportBoundary = function(stage) {
            if (stage == "before_native_archive") stop("injected save failure")
        }, .package = "MultiScholaR")
        expect_identical(digest::digest(file = publication$archive_path), hash)
        expect_identical(workflow$state_manager, native)
        invalid <- file.path(root, "bad.msdia.tar"); writeLines("invalid", invalid)
        expect_error(openDIAPeptideCheckpoint(invalid, workflow, session, paths))
        expect_identical(workflow$state_manager, native)
        relocated <- list(source_dir = file.path(root, "new-source"), peptide_qc_dir = file.path(root, "new-peptide"))
        lapply(relocated, dir.create)
        for (boundary in c("peptide_checkpoint_after_owner_swap", "peptide_checkpoint_after_metadata")) {
            old_owner <- session$userData$dia_capacity_manager
            with_mocked_bindings({
                expect_error(openDIAPeptideCheckpoint(publication$archive_path, workflow, session, relocated),
                    "injected restore failure")
            }, workflowExecutionBoundary = function(stage) {
                if (stage == boundary) stop("injected restore failure")
            }, .package = "MultiScholaR")
            expect_identical(workflow$state_manager, native)
            expect_identical(session$userData$dia_capacity_manager, old_owner)
            expect_identical(workflow$design_matrix, object@design_matrix)
            expect_identical(digest::digest(native$states), original)
        }
        expect_identical(openDIAPeptideCheckpoint(publication$archive_path, workflow, session, relocated), "sample_filtered")
        manager <- workflow$state_manager
        expect_true(usesDIAReferenceStates(manager))
        expect_identical(manager$getState("sample_filtered")@peptide_data, changed@peptide_data)
        expect_identical(manager$getState("intensity_filtered")@peptide_data, object@peptide_data)
        expect_identical(workflow$design_matrix, object@design_matrix)
        expect_s3_class(workflow$config_list$rollUpPrecursorToPeptide$core_utilisation, "DIAWorkerSpec")
        expect_identical(readLines(file.path(relocated$peptide_qc_dir, "peptide.tsv")), "accepted output")
        expect_identical(manager$.context$peptide_design_restore$design_matrix, object@design_matrix)
        closeDIAStateAdapter(manager)
    })
})

test_that("checkpoint route, metadata and pending-export guards reject unsupported work", {
    expect_error(validateDIAPeptideCheckpointRoute(list(workflow_type = "DIA_LIMPA")),
        class = "saved_entity_unsupported")
    expect_error(validateDIAPeptideCheckpointRoute(list(config_list = list(
        globalParameters = list(use_limpa = TRUE)))), class = "saved_entity_unsupported")
    expect_error(validateDIAPeptideCheckpointMetadata(list(checkpoint_kind = "peptide-qc")),
        class = "saved_entity_invalid")
    calls <- list()
    manager <- structure(list(), class = "DIAStateAdapter")
    with_mocked_bindings({
        checkDIAPeptideCheckpointIdle(list(state_manager = manager), list(userData = list()), saving = TRUE)
        checkDIAPeptideCheckpointIdle(list(state_manager = manager), list(userData = list()))
    }, assertDIAExportReady = function(manager, allow_action = NULL) {
        calls[[length(calls) + 1L]] <<- list(allow_action = allow_action)
    }, .package = "MultiScholaR")
    expect_identical(calls[[1L]]$allow_action, "session-export")
    expect_null(calls[[2L]]$allow_action)
})

test_that("checkpoint cutoff waits for its input binding and preserves later edits", {
    workflow <- shiny::reactiveValues(state_manager = list(current_state = "imputed"))
    messages <- list()
    mock_session <- shiny::MockShinySession$new()
    mock_session$sendInputMessage <- function(inputId, message) {
        messages[[length(messages) + 1L]] <<- list(id = inputId, message = message)
    }
    shiny::testServer(mod_prot_qc_peptide_sample_server, session = mock_session,
        args = list(workflow_data = workflow, omic_type = "proteomics", experiment_label = "restore"), {
        workflow$qc_params <- list(peptide_qc = list(sample_filter = list(min_peptides_per_sample = 100)))
        workflow$native_session_restored <- 1L
        session$flushReact()
        expect_length(messages, 0L)
        session$setInputs(min_peptides_per_sample = 500)
        expect_length(messages, 1L)
        expect_match(messages[[1L]]$id, "min_peptides_per_sample$")
        expect_identical(messages[[1L]]$message$value, "100")
        session$setInputs(min_peptides_per_sample = 100)
        session$setInputs(min_peptides_per_sample = 300)
        expect_length(messages, 1L)
        workflow$qc_params$peptide_qc$sample_filter$min_peptides_per_sample <- 200
        workflow$native_session_restored <- 2L
        session$flushReact()
        expect_length(messages, 2L)
        expect_identical(messages[[2L]]$message$value, "200")
    })
})
