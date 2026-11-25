# KIT-271 (VAL-317): on DIA with limpa the correlation filter's table is the protein table before any RUV step, in the
# none, covariates and correction modes alike; it is stored as normalised_results_cln_with_replicates.tsv (not
# ruv_normalised_...) and Copy delivers it as normalised_results.tsv (not RUV_normalised_results.tsv). A project written
# by an earlier build (only the ruv_ file on disk) keeps working: the export lists the file that exists and Copy
# delivers it under the new name; with both present the new one wins.
# The fixture and its helpers are test-k271-limpa-pass-through.R's (the limpa owner on protNormW1LimpaOwner's recipe,
# the correlation and the RUV covariates in process and through the worker), copied under the prefix k271 and never
# sourced; Copy is driven as test-kit257-norm-own-files.R drives it (k257Copy, k257CopyPaths, copied), here inside the
# owner's state scope because Copy reads the route from the active owner. The new function is called only where it
# exists (rows T6 b, T8b), so the file runs on the base and on the head and the failing rows fail at a value.
# Lines starting `KIT-271` are the evidence of docs/qa/workflow-execution/kit/kit-271/.

k271Quiet <- function(code) suppressWarnings(suppressMessages(code))
k271Noop <- function(...) invisible(NULL)
k271Aesthetics <- function() list(color_var = "group", shape_var = NULL)
k271Settings <- function(mode = "covariates") list(norm_method = "quantile", ruv_mode = mode, auto_percentage_min = NULL,
    auto_percentage_max = NULL, separation_metric = NULL, k_penalty_weight = NULL, adaptive_k_penalty = NULL,
    ruv_percentage = 50, ruv_k = 1, color_variable = "group", shape_variable = NULL)
k271Input <- function(mode = "covariates") list(norm_method = "quantile", ruv_mode = mode, ruv_percentage = 50,
    ruv_k = 1, color_variable = "group", shape_variable = NULL)
k271Paths <- function() c("in process", if (workflowOwnerCapable()) "worker")

k271Cache <- new.env(parent = emptyenv())
k271Objects <- function() {
    if (is.null(k271Cache$rolled)) {
        k271Cache$precursors <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
        k271Cache$rolled <- rollUpPrecursorToPeptide(k271Cache$precursors, core_utilisation = NA)
    }
    k271Cache
}
k271Quantified <- function() {
    cache <- k271Objects()
    if (is.null(cache$quantified)) cache$quantified <- k271Quiet(proteinMissingValueImputationLimpa(
        log2TransformPeptideMatrix(calcPeptideMatrix(cache$rolled)), dpc_slope = 0.8, verbose = FALSE, chunk = 1000,
        engine = "limpa"))
    cache$quantified
}

k271HeadRef <- function(fixture) withDIAStateScope(fixture$manager, function() getDIAStateReference(fixture$manager))
k271Name <- function(fixture, label, reference) {
    assign(label, reference, envir = fixture$refs)
    invisible(reference)
}
k271Label <- function(fixture, entity_id) {
    if (is.null(entity_id)) return("<null>")
    labels <- sort(ls(fixture$refs))
    hit <- labels[vapply(labels, function(label) identical(get(label, envir = fixture$refs)$entity_id, entity_id), NA)]
    if (length(hit)) paste(hit, collapse = "|") else paste0("other:", substr(entity_id, 1L, 8L))
}
k271Known <- function(fixture, entity_id) !startsWith(k271Label(fixture, entity_id), "other:")
k271HeadLabel <- function(fixture) k271Label(fixture, k271HeadRef(fixture)$entity_id)
k271Recorded <- function(fixture, norm) k271Label(fixture, norm$normalization_input_reference$entity_id)

k271Commit <- function(fixture, state, object) {
    k271Quiet(fixture$manager$saveState(state, object, list(), "k271"))
    k271Name(fixture, state, k271HeadRef(fixture))
}

#' A DIA with limpa owner. `peptide = FALSE`: protNormW1LimpaOwner's commits (protein_dpc_quant, protein_normalised,
#' peptide_count_filtered). `peptide = TRUE`: kit225LimpaOwner's peptide states first. Stops at `upto`.
k271LimpaOwner <- function(peptide = FALSE, upto = "peptide_count_filtered", envir = parent.frame()) {
    withr::local_options(multischolar.dia_native_sessions = TRUE, multischolar.capture_test_checkpoints = FALSE,
        .local_envir = envir)
    root <- withr::local_tempdir(.local_envir = envir)
    paths <- list(source_dir = file.path(root, "scripts", "proteomics"),
        protein_qc_dir = file.path(root, "results", "proteomics", "protein_qc"),
        peptide_qc_dir = file.path(root, "results", "proteomics", "peptide_qc"))
    paths$publication_graphs_dir <- file.path(root, "results", "proteomics", "publication_graphs")
    paths$time_dir <- file.path(paths$publication_graphs_dir, "filtering_qc")
    for (path in unlist(paths)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(paths$source_dir, "k271-limpa", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager), envir = envir)
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("DIA_LIMPA")
    manager$.context$experiment_paths <- paths
    workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list()), config_list = list()))
    fixture <- list(root = root, paths = paths, workflow = workflow, manager = manager,
        refs = new.env(parent = emptyenv()), seen = new.env(parent = emptyenv()))
    cache <- k271Objects()
    states <- if (peptide) c("raw_data_s4", "qvalue_filtered", "precursor_rollup", "log2_matrix",
        "peptide_accession_cleaned") else character()
    for (state in states) {
        object <- if (state %in% c("raw_data_s4", "qvalue_filtered")) cache$precursors else cache$rolled
        k271Commit(fixture, state, object)
        if (identical(state, upto)) return(fixture)
    }
    setDIAStepProvenance(manager, "design_builder")
    k271Commit(fixture, "protein_dpc_quant", k271Quantified())
    if (identical(upto, "protein_dpc_quant")) return(fixture)
    k271Quiet(withDIAStateScope(manager, function() {
        v1 <- computeDIAReferenceProteinNormalisation(manager, getDIAStateReference(manager), list(norm_method = "none"))
        setDIAStepProvenance(manager, "protein_between_sample_normalisation", "whole", "norm_method", NULL,
            values = list(norm_method = "none"))
        finishDIAReferenceProteinStep(workflow, v1, publish = FALSE)
        k271Name(fixture, manager$current_state, getDIAStateReference(manager))
        if (identical(upto, "protein_normalised")) return(invisible(NULL))
        v2 <- computeDIAReferenceProteinPeptideCountFilter(manager, getDIAStateReference(manager),
            list(num_peptides_per_protein_thresh = 1, num_peptidoforms_per_protein_thresh = 2))
        setDIAStepProvenance(manager, "protein_peptide_count_filter", "whole",
            c("num_peptides_per_protein_thresh", "num_peptidoforms_per_protein_thresh"), NULL, values = v2$config[1:2])
        finishDIAReferenceProteinStep(workflow, v2, publish = FALSE)
        k271Name(fixture, manager$current_state, getDIAStateReference(manager))
    }))
    fixture
}

#' protNormW2Mocks (test-prot-norm-core.R:1217-1221, file-local there), copied.
k271W2Mocks <- function(env = parent.frame()) {
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
        finishWorkflowExecutionJob = function(owner, status) { owner$job <- NULL; list(status = status) },
        .package = "MultiScholaR", .env = env)
}

#' The (A) owner on a path; on the worker path the W2 mocks in the test's frame and a controller.
k271Start <- function(path, envir = parent.frame()) {
    fixture <- k271LimpaOwner(envir = envir)
    if (identical(path, "worker")) {
        k271W2Mocks(env = envir)
        fixture$seen$controller <- omicsWorkerController(list(owner = fixture$manager))
    }
    fixture
}

k271Value <- function(fixture, reference) {
    manager <- fixture$manager
    state <- withWorkflowStateScope(manager, function() {
        object <- resolveDIAStateValue(manager$.context, reference)$data
        if (!isS4(object)) object <- rebuildDIACheckpointObject(manager$.context, object, TRUE, list())
        slots <- setdiff(methods::slotNames(object), "args")
        list(slots = lapply(stats::setNames(nm = slots), function(slot) methods::slot(object, slot)),
            args = object@args[intersect(c("limpa_dpc_quant_results", "ruvCovariates"), names(object@args))],
            samples = setdiff(colnames(object@protein_quant_table), object@protein_id_column))
    })
    list(samples = state$samples, sha = digest::digest(state[c("slots", "args")], algo = "sha256"))
}

#' A failure as "class: message", NA when the code completed.
k271Error <- function(code) tryCatch({
    k271Quiet(code)
    NA_character_
}, error = function(e) paste0(class(e)[[1L]], ": ", conditionMessage(e)))


#' The pass-through (in process on both paths, as protNormW2Clicks$pass_through); the tables its figures were drawn from.
k271Press <- function(fixture, norm) {
    drawn <- character()
    post <- function(reference) {
        drawn <<- c(drawn, k271Label(fixture, reference$entity_id))
        invisible(NULL)
    }
    error <- k271Error(runDIAReferenceLimpaNormalisationPassThrough(k271Settings(), fixture$workflow, norm, post))
    list(error = error, figures = drawn)
}

k271Dispatch <- function(fixture) function(operation, parameters, input, state_name = NULL) {
    context <- if (!is.null(state_name)) list(state_name = state_name, config = list(), description = "dispatch")
    fixture$seen$inputs <- c(fixture$seen$inputs, paste0(operation, ":", k271Label(fixture, input$entity_id)))
    protNormW2Dispatch(fixture, fixture$seen$controller, operation, parameters, input, context)
}

#' The managed correlation click (protNormW2Correlation, test-prot-norm-core.R:1135-1149).
k271WorkerCorrelation <- function(fixture, norm, threshold) {
    manager <- fixture$manager
    source <- withDIAStateScope(manager, function() {
        assertDIAExportReady(manager, allow_action = "correlation")
        source <- diaCorrelationSource(manager, norm)
        readDIAStateMetadata(manager, source)
        source
    })
    dispatch <- k271Dispatch(fixture)
    pairs <- dispatch("pearson-pairs", diaPearsonPairsParameters("group"), source)
    selection <- dispatch("pearson-selection", list(pairs = pairs$reference$entity_id, threshold = threshold), source,
        "correlation_filtered")
    withDIAStateScope(manager, function() publishDIAReferenceCorrelation(manager, list(pairs = pairs$reference,
        final = selection$value), fixture$workflow, norm, fixture$paths, threshold, FALSE, "proteomics", "core",
        k271Aesthetics))
}

#' Correlation Apply at `threshold`: its source, the worker's inputs, whether it committed, its result.
k271Correlate <- function(fixture, norm, path, threshold, label) {
    manager <- fixture$manager
    before <- k271HeadRef(fixture)$entity_id
    history <- diaStateHistory(manager)
    source <- tryCatch(withDIAStateScope(manager, function() diaCorrelationSource(manager, norm)), error = function(e) NULL)
    fixture$seen$inputs <- character()
    final <- NULL
    error <- k271Error(final <- if (identical(path, "worker")) k271WorkerCorrelation(fixture, norm, threshold) else
        runDIAReferenceCorrelation(fixture$workflow, norm, fixture$paths, threshold, FALSE, "proteomics", "core",
            function() "group", k271Aesthetics))
    head <- k271HeadRef(fixture)
    out <- list(error = error, source = k271Label(fixture, source$entity_id), inputs = fixture$seen$inputs,
        committed = !identical(before, head$entity_id), history_changed = !identical(history, diaStateHistory(manager)))
    if (is.null(final)) return(out)
    if (!k271Known(fixture, final$entity_id)) k271Name(fixture, label, final)
    if (!k271Known(fixture, head$entity_id)) k271Name(fixture, paste0(label, "-state"), head)
    c(out, list(head = k271Label(fixture, head$entity_id)), k271Value(fixture, final))
}

#' The RUV covariates (in process, or through the worker as protNormW2RuvCovariates, test-prot-norm-core.R:1153-1166).
k271RuvCovariates <- function(fixture, norm, path, mode) {
    manager <- fixture$manager
    settings <- list(ruv_mode = mode, color_variable = "group", shape_variable = NULL)
    if (!identical(path, "worker")) return(k271Error(runDIAReferenceRuvCovariates(settings, fixture$workflow, norm,
        fixture$paths, "proteomics", "core", k271Noop, function() "group")))
    k271Error({
        request <- withDIAStateScope(manager, function() {
            assertDIAExportReady(manager, allow_action = "ruv-covariates")
            list(source = diaRuvCovariatesSource(manager, norm),
                parameters = resolveDIARuvCovariatesParameters(fixture$workflow$config_list, "group", mode))
        })
        result <- k271Dispatch(fixture)("protein_ruv_covariates", request$parameters, request$source,
            if (identical(mode, "correction")) "ruv_corrected" else "ruv_covariates")
        withDIAStateScope(manager, function() publishDIAReferenceRuvCovariates(manager, settings,
            result[c("value", "reference")], fixture$workflow, norm, fixture$paths, k271Noop))
    })
}

#' RUV covariates set to Skip through the tab's observer, with the flag the tab sets after the correlation (R1-03).
k271RuvSkip <- function(fixture, norm) {
    notes <- character()
    note <- function(ui, ...) {
        notes <<- c(notes, paste(as.character(ui), collapse = " "))
        invisible(NULL)
    }
    norm$correlation_filtering_complete <- TRUE
    k271Quiet(runProtNormRuvCovariatesObserver(list(ruv_mode = "skip", color_variable = "group", shape_variable = NULL),
        new.env(), fixture$workflow, norm, fixture$paths, "proteomics", "core", k271Noop, function() "group",
        renderTextFn = function(x) x, showNotificationFn = note, messageFn = k271Noop))
    notes
}

# ---- KIT-271's own helpers ----

k271New <- "normalised_results_cln_with_replicates.tsv"
k271Old <- "ruv_normalised_results_cln_with_replicates.tsv"
k271Fn <- function(name) get0(name, envir = asNamespace("MultiScholaR"), inherits = FALSE)
k271HasRoute <- function() !is.null(k271Fn("resolveProtNormCorrelationTableNames"))
k271Log <- function(case, ...) cat(sprintf("KIT-271 %s %s\n", case, paste(..., collapse = " ")))
k271Need <- function() {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
}
k271Table <- function(fixture, name) file.path(fixture$paths$protein_qc_dir, name)
#' Which of the two names the protein QC folder holds.
k271Present <- function(fixture) c(k271New, k271Old)[file.exists(k271Table(fixture, c(k271New, k271Old)))]
k271Sha <- function(path) digest::digest(file = path, algo = "sha256")

#' The limpa owner after the pass-through and the correlation filter (Apply at 0.5) on `path`.
k271Project <- function(path, envir = parent.frame()) {
    fixture <- k271Start(path, envir = envir)
    norm <- new.env(parent = emptyenv())
    pressed <- k271Press(fixture, norm)
    correlated <- k271Correlate(fixture, norm, path, 0.5, "C1")
    k271Log("setup", path, "pass-through:", pressed$error, "correlation:", correlated$error, "state:",
        fixture$manager$current_state, "protein_qc:", paste(list.files(fixture$paths$protein_qc_dir, pattern = "tsv$"),
        collapse = ","))
    list(fixture = fixture, norm = norm, path = path, errors = c(pressed$error, correlated$error))
}

#' Copy (`copyToResultsSummary`) on `paths`: its `COPY:` lines, its failures and what reached the summary folders
#' (test-kit257-norm-own-files.R's k257Copy, copied).
k271CopyRun <- function(paths, omic) {
    for (path in unlist(paths)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    testthat::local_mocked_bindings(getProjectPaths = function(...) paths, .package = "MultiScholaR")
    failures <- NULL
    console <- utils::capture.output(failures <- suppressWarnings(suppressMessages(copyToResultsSummary(omic, "k271",
        contrasts_tbl = data.frame(contrasts = "groupB-groupA"), design_matrix = data.frame(Run = "S1", group = "A"),
        force = TRUE))))
    summary <- paths$results_summary_dir
    list(console = grep("^COPY:", console, value = TRUE),
        failed = vapply(failures, function(failure) basename(as.character(failure$source %||% "")), ""),
        tables = sort(list.files(file.path(summary, "Publication_tables"))),
        delivered = file.path(summary, "Publication_tables"))
}

#' Copy's paths beside the project's own (test-kit257-norm-own-files.R's k257CopyPaths, copied).
k271CopyPaths <- function(paths, root, feature_qc_dir) {
    defaults <- list(results_dir = file.path(root, "copy", "results"),
        results_summary_dir = file.path(root, "copy", "results_summary"), publication_graphs_dir = file.path(root, "copy", "pub"),
        time_dir = file.path(root, "copy", "time"), qc_dir = file.path(root, "copy", "qc"),
        da_output_dir = file.path(root, "copy", "da"), pathway_dir = file.path(root, "copy", "pathway"),
        source_dir = file.path(root, "copy", "source"))
    for (name in names(defaults)) if (is.null(paths[[name]])) paths[[name]] <- defaults[[name]]
    paths$feature_qc_dir <- feature_qc_dir
    paths
}

#' Copy as the tab's Summary runs it: inside the owner's state scope (the route is read from the active owner).
k271Copy <- function(fixture) {
    paths <- k271CopyPaths(fixture$paths, fixture$root, fixture$paths$protein_qc_dir)
    withDIAStateScope(fixture$manager, function() k271CopyRun(paths, "proteomics"))
}

#' The project as a native session (the flag the app's import and restore set), set after the steps ran.
k271Native <- function(fixture) {
    fixture$manager$.context$session_format <- "native_v1"
    invisible(fixture)
}

#' The export's metadata and file list for the project (the tab's `collectDIANativeSessionMetadata` and
#' `collectDIANativeGuiFiles`), and, with `export`, the native export itself (`runDIANativeSessionExport`, the session
#' summary written by its prepare step); failures come back as text.
k271Export <- function(fixture, norm, ruv_mode, export = FALSE, stored_results_file = NULL) {
    out <- list(error = NA_character_)
    tryCatch(k271Quiet(withDIAStateScope(fixture$manager, function() {
        collected <- collectDIANativeSessionMetadata(fixture$workflow, norm, list(norm_method = "quantile", ruv_mode = ruv_mode))
        out$results_file <<- collected$metadata$results_file
        metadata <- collected$metadata
        if (!is.null(stored_results_file)) metadata$results_file <- stored_results_file
        files <- collectDIANativeGuiFiles(fixture$manager, metadata)
        out$destinations <<- vapply(files, function(spec) spec$destination, "", USE.NAMES = FALSE)
        if (export) {
            source_dir <- fixture$paths$source_dir
            result <- runDIANativeSessionExport(fixture$manager, metadata, source_dir, extra_references = collected$references,
                dependencies = files, prepare_outputs = function(metadata, request) writeLines("kit271",
                    file.path(source_dir, "filtered_session_summary.txt")))
            out$archive <<- result$publication$archive_path
        }
    })), error = function(e) out$error <<- paste0(class(e)[[1L]], ": ", conditionMessage(e)))
    out
}

# The classic DIA route's clicks (test-kit257-norm-own-files.R's k257ProtSettings, k257ProtClick, k257ProtCorrelation).
k271ProtSettings <- function(mode = "skip") list(norm_method = "quantile", ruv_mode = mode, auto_percentage_min = NULL,
    auto_percentage_max = NULL, separation_metric = NULL, k_penalty_weight = NULL, adaptive_k_penalty = NULL,
    ruv_percentage = 50, ruv_k = 1, color_variable = "group", shape_variable = NULL)
k271ProtClick <- function(fixture, norm, mode = "skip") {
    protNormCoreQuiet(runDIAReferenceNormalization(k271ProtSettings(mode), fixture$workflow, norm, fixture$paths,
        "proteomics", "core", function(reference) invisible(NULL), function(reference) invisible(NULL),
        function() "group", pre_qc = NULL))
}
k271ProtCorrelation <- function(fixture, norm, threshold = 0.5) {
    final <- protNormCoreQuiet(runDIAReferenceCorrelation(fixture$workflow, norm, fixture$paths, threshold, FALSE,
        "proteomics", "core", function() "group", k271Aesthetics))
    norm$correlation_filtering_complete <- TRUE
    final
}

# ---- T1 to T4: the limpa project's table, its name and Copy's delivery ----

test_that("KIT-271 T1: DIA with limpa, after the correlation filter the protein QC folder holds normalised_ and no ruv_ table (fails before)", {
    k271Need()
    for (path in k271Paths()) {
        project <- k271Project(path)
        present <- k271Present(project$fixture)
        k271Log("T1", path, "present:", paste(present, collapse = ","))
        expect_identical(project$errors, c(NA_character_, NA_character_), info = path)
        expect_identical(present, k271New, info = path)
    }
})

test_that("KIT-271 T2: the RUV covariates action, in either mode, leaves the table's name and content unchanged (name row fails before)", {
    k271Need()
    for (mode in c("covariates", "correction")) {
        project <- k271Project("in process")
        fixture <- project$fixture
        before <- k271Present(fixture)
        sha_before <- k271Sha(k271Table(fixture, before[[1L]]))
        expect_identical(k271RuvCovariates(fixture, project$norm, "in process", mode), NA_character_, info = mode)
        expect_identical(fixture$manager$current_state, if (identical(mode, "correction")) "ruv_corrected" else "ruv_covariates",
            info = mode)
        after <- k271Present(fixture)
        k271Log("T2", mode, "present:", paste(after, collapse = ","), "table sha256 before/after:", sha_before,
            if (length(after)) k271Sha(k271Table(fixture, after[[1L]])) else "none",
            "adjusted matrix:", file.exists(k271Table(fixture, "ruv_adjusted_protein_matrix.tsv")))
        expect_identical(after, k271New, info = mode)
        expect_identical(after, before, info = mode)
        expect_identical(if (length(after)) k271Sha(k271Table(fixture, after[[1L]])) else NA_character_, sha_before, info = mode)
        expect_true(file.exists(k271Table(fixture, "ruv_adjusted_protein_matrix.tsv")), info = mode)
    }
})

test_that("KIT-271 T3: DIA with limpa, after Skip RUV the table is the same normalised_ table (fails before)", {
    k271Need()
    for (path in k271Paths()) {
        project <- k271Project(path)
        before <- k271Present(project$fixture)
        notes <- k271RuvSkip(project$fixture, project$norm)
        after <- k271Present(project$fixture)
        k271Log("T3", path, "notes:", paste(notes, collapse = " | "), "present:", paste(after, collapse = ","))
        expect_identical(after, before, info = path)
        expect_identical(after, k271New, info = path)
    }
})

test_that("KIT-271 T4: Copy after the correlation delivers normalised_results.tsv, not RUV_normalised_results.tsv, with the table's values (fails before)", {
    k271Need()
    for (path in k271Paths()) {
        project <- k271Project(path)
        fixture <- project$fixture
        table <- k271Present(fixture)
        sha_table <- k271Sha(k271Table(fixture, table[[1L]]))
        copy <- k271Copy(fixture)
        k271Log("T4", path, "copy:", paste(copy$console, collapse = " | "), "tables:", paste(copy$tables, collapse = ","))
        expect_true("normalised_results.tsv" %in% copy$tables, info = path)
        expect_false("RUV_normalised_results.tsv" %in% copy$tables, info = path)
        delivered <- if ("normalised_results.tsv" %in% copy$tables) file.path(copy$delivered, "normalised_results.tsv") else
            file.path(copy$delivered, "RUV_normalised_results.tsv")
        k271Log("A5", path, "table sha256:", sha_table, "delivered sha256:", k271Sha(delivered))
        expect_identical(k271Sha(delivered), sha_table, info = path)
    }
})

# ---- T5, T5b: the export's metadata and file list ----

test_that("KIT-271 T5: the export names the table normalised_ in its metadata and file list (fails before)", {
    k271Need()
    project <- k271Project("in process")
    k271Native(project$fixture)
    out <- k271Export(project$fixture, project$norm, "covariates")
    k271Log("T5", "error:", out$error, "results_file:", out$results_file, "destinations:",
        paste(grep("protein_qc/.*normalised", out$destinations, value = TRUE), collapse = ","))
    expect_identical(out$error, NA_character_)
    expect_identical(out$results_file, k271New)
    expect_true(paste0("protein_qc/", k271New) %in% out$destinations)
    expect_false(paste0("protein_qc/", k271Old) %in% out$destinations)
})

test_that("KIT-271 T5b: after the RUV covariates action, in either mode, the export names the same table and succeeds (name row fails before)", {
    k271Need()
    for (mode in c("covariates", "correction")) {
        project <- k271Project("in process")
        expect_identical(k271RuvCovariates(project$fixture, project$norm, "in process", mode), NA_character_, info = mode)
        k271Native(project$fixture)
        out <- k271Export(project$fixture, project$norm, mode, export = TRUE)
        k271Log("T5b", mode, "error:", out$error, "results_file:", out$results_file, "archive written:",
            !is.null(out$archive) && file.exists(out$archive))
        expect_identical(out$error, NA_character_, info = mode)
        expect_identical(out$results_file, k271New, info = mode)
        expect_true(paste0("protein_qc/", k271New) %in% out$destinations, info = mode)
        expect_true(!is.null(out$archive) && file.exists(out$archive), info = mode)
    }
})

# ---- T6: the names, by the resolver (unchanged) and by the new function ----

test_that("KIT-271 T6 a: the resolver is unchanged: no result and an applied result give ruv_, a skipped result gives normalised_ (guard)", {
    resolver <- k271Fn("resolveProtNormCorrelationResultFilenames")
    quiet <- function(...) invisible(NULL)
    expect_identical(resolver(NULL, NULL, messageFn = quiet)$tsvFilename, k271Old)
    expect_identical(resolver(list(ruv_skipped = FALSE), NULL, messageFn = quiet)$tsvFilename, k271Old)
    expect_identical(resolver(NULL, list(ruv_skipped = FALSE), messageFn = quiet)$tsvFilename, k271Old)
    expect_identical(resolver(list(ruv_skipped = TRUE), NULL, messageFn = quiet)$tsvFilename, k271New)
    expect_identical(resolver(NULL, list(ruv_skipped = TRUE), messageFn = quiet)$tsvFilename, k271New)
    expect_identical(resolver(NULL, NULL, messageFn = quiet)$rdsFilename, "ruv_normalised_results_cln_with_replicates.RDS")
})

test_that("KIT-271 T6 b: the new function gives the resolver's value on every route but DIA with limpa, and normalised_ names there", {
    skip_if_not(k271HasRoute(), "the new function does not exist on this commit")
    resolveTable <- k271Fn("resolveProtNormCorrelationTableNames")
    resolver <- k271Fn("resolveProtNormCorrelationResultFilenames")
    quiet <- function(...) invisible(NULL)
    results <- list(none = list(NULL, NULL), applied = list(list(ruv_skipped = FALSE), NULL),
        skipped = list(list(ruv_skipped = TRUE), NULL))
    for (route in c("DIA", "LFQ", "TMT")) {
        for (case in names(results)) {
            expect_identical(resolveTable(list(workflow_type = route), results[[case]][[1L]], results[[case]][[2L]],
                messageFn = quiet), resolver(results[[case]][[1L]], results[[case]][[2L]], messageFn = quiet),
                info = paste(route, case))
        }
    }
    expect_identical(resolveTable(NULL, NULL, NULL, messageFn = quiet), resolver(NULL, NULL, messageFn = quiet))
    for (case in names(results)) {
        expect_identical(resolveTable(list(workflow_type = "DIA_LIMPA"), results[[case]][[1L]], results[[case]][[2L]],
            messageFn = quiet)[c("tsvFilename", "rdsFilename")], list(tsvFilename = k271New,
            rdsFilename = "normalised_results_cln_with_replicates.RDS"), info = paste("DIA_LIMPA", case))
    }
    lines <- character()
    resolveTable(list(workflow_type = "DIA_LIMPA"), NULL, NULL, messagePrefix = "*** CORRELATION",
        messageFn = function(text) lines <<- c(lines, text))
    expect_length(lines, 1L)
    expect_match(lines, "DIA with limpa", fixed = TRUE)
    expect_false(grepl("RUV was applied", lines, fixed = TRUE))
})

# ---- T7: the classic DIA route keeps following the RUV choice ----

test_that("KIT-271 T7: classic DIA, the correlation's table is named by the RUV choice, skipped or applied (guard)", {
    for (mode in c("skip", "manual")) {
        fixture <- protNormCoreOwner()
        norm <- new.env(parent = emptyenv())
        k271ProtClick(fixture, norm, mode)
        k271ProtCorrelation(fixture, norm, 0.5)
        present <- k271Present(fixture)
        k271Log("T7", mode, "present:", paste(present, collapse = ","))
        expect_identical(present, if (identical(mode, "skip")) k271New else k271Old, info = mode)
    }
})

# ---- T8: a project written by an earlier build (only the ruv_ table on disk) ----

#' The project as an earlier build left it: the table under the ruv_ name only.
k271Earlier <- function(fixture) {
    if (file.exists(k271Table(fixture, k271New))) file.rename(k271Table(fixture, k271New), k271Table(fixture, k271Old))
    invisible(k271Present(fixture))
}

test_that("KIT-271 T8: an earlier build's project: the export lists the file that exists and succeeds (guard)", {
    k271Need()
    project <- k271Project("in process")
    on_disk <- k271Earlier(project$fixture)
    expect_identical(on_disk, k271Old)
    k271Native(project$fixture)
    out <- k271Export(project$fixture, project$norm, "covariates", export = TRUE)
    stored <- k271Export(project$fixture, project$norm, "covariates", stored_results_file = k271Old)
    k271Log("T8", "export error:", out$error, "results_file:", out$results_file, "destinations:",
        paste(grep("protein_qc/.*normalised", out$destinations, value = TRUE), collapse = ","),
        "stored ruv_ name destinations:", paste(grep("protein_qc/.*normalised", stored$destinations, value = TRUE), collapse = ","))
    expect_identical(out$error, NA_character_)
    expect_true(!is.null(out$archive) && file.exists(out$archive))
    expect_true(paste0("protein_qc/", k271Old) %in% out$destinations)
    expect_false(paste0("protein_qc/", k271New) %in% out$destinations)
    expect_true(paste0("protein_qc/", k271Old) %in% stored$destinations)
})

test_that("KIT-271 T8 (Copy): an earlier build's project: Copy delivers the table it holds as normalised_results.tsv (fails before)", {
    k271Need()
    project <- k271Project("in process")
    k271Earlier(project$fixture)
    copy <- k271Copy(project$fixture)
    k271Log("T8 copy:", paste(copy$console, collapse = " | "), "tables:", paste(copy$tables, collapse = ","))
    expect_true("normalised_results.tsv" %in% copy$tables)
    expect_false("RUV_normalised_results.tsv" %in% copy$tables)
})

test_that("KIT-271 T8b: both tables present: the export and Copy name the new one (head only)", {
    skip_if_not(k271HasRoute(), "the new function does not exist on this commit")
    k271Need()
    project <- k271Project("in process")
    fixture <- project$fixture
    sha_new <- k271Sha(k271Table(fixture, k271New))
    writeLines("an earlier build's table", k271Table(fixture, k271Old))
    copy <- k271Copy(fixture)
    delivered <- file.path(copy$delivered, "normalised_results.tsv")
    k271Native(fixture)
    out <- k271Export(fixture, project$norm, "covariates")
    k271Log("T8b", "copy:", paste(copy$console, collapse = " | "), "tables:", paste(copy$tables, collapse = ","),
        "destinations:", paste(grep("protein_qc/.*normalised", out$destinations, value = TRUE), collapse = ","))
    expect_true("normalised_results.tsv" %in% copy$tables)
    expect_false("RUV_normalised_results.tsv" %in% copy$tables)
    if (file.exists(delivered)) expect_identical(k271Sha(delivered), sha_new)
    expect_identical(out$error, NA_character_)
    expect_true(paste0("protein_qc/", k271New) %in% out$destinations)
    expect_false(paste0("protein_qc/", k271Old) %in% out$destinations)
})
