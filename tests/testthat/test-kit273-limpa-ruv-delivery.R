# KIT-273 (VAL-317 (2), VAL-332): on DIA with limpa, Copy to Results Summary delivers what the differential analysis used:
# in the RUV-III correction mode the corrected matrix (Publication_tables/ruv_corrected_results.tsv), in the covariates
# mode the whole design with the factors W1..Wk (Study_report/design_matrix_with_ruv_factors.tab), each written from the
# store by the current differential analysis's own record and source, and nothing when the RUV step was skipped or applied
# no correction, when the store holds no differential analysis, or when the analysis is not current.
# The limpa fixture is test-kit244-limpa-pass-through.R's (the limpa owner on protNormW1LimpaOwner's recipe; the
# correlation and the RUV covariates in process and through the worker), the DA click is test-workflow-saved-da.R's
# (savedDaW1Click), Copy's driver is test-kit271-limpa-table-name.R's (k271Copy, itself test-kit257-norm-own-files.R's);
# all copied under the prefix k273 and never sourced. The delivered names are two constants (VAL-332; the owner's wording
# pass changes the constants and the code together). Lines starting `KIT-273` are the evidence of
# docs/qa/workflow-execution/kit/kit-273/.
# Value comparisons are made on the text both sides are written as (the same vroom writer): a delivered file is text, so
# the DA's in-memory values are written by that writer and both texts compared with identical() (tolerance 0 on the
# representation the package delivers; no precision is lost to a second formatting).

k273Matrix <- "Publication_tables/ruv_corrected_results.tsv"
k273Design <- "Study_report/design_matrix_with_ruv_factors.tab"

k273Quiet <- function(code) suppressWarnings(suppressMessages(code))
k273Noop <- function(...) invisible(NULL)
k273Aesthetics <- function() list(color_var = "group", shape_var = NULL)
k273Paths <- function() c("in process", if (workflowOwnerCapable()) "worker")
k273Log <- function(case, ...) cat(sprintf("KIT-273 %s %s\n", case, paste(..., collapse = " ")))
k273Need <- function() {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    skip_if_not_installed("Glimma")
}
k273Sha <- function(path) digest::digest(file = path, algo = "sha256")

# ---- fixtures ----

k273Cache <- new.env(parent = emptyenv())
k273Objects <- function() {
    if (is.null(k273Cache$rolled)) {
        k273Cache$precursors <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
        k273Cache$rolled <- rollUpPrecursorToPeptide(k273Cache$precursors, core_utilisation = NA)
    }
    k273Cache
}
k273Quantified <- function() {
    cache <- k273Objects()
    if (is.null(cache$quantified)) cache$quantified <- k273Quiet(proteinMissingValueImputationLimpa(
        log2TransformPeptideMatrix(calcPeptideMatrix(cache$rolled)), dpc_slope = 0.8, verbose = FALSE, chunk = 1000,
        engine = "limpa"))
    cache$quantified
}

k273HeadRef <- function(fixture) withDIAStateScope(fixture$manager, function() getDIAStateReference(fixture$manager))
k273Name <- function(fixture, label, reference) {
    assign(label, reference, envir = fixture$refs)
    invisible(reference)
}
k273Label <- function(fixture, entity_id) {
    if (is.null(entity_id)) return("<null>")
    labels <- sort(ls(fixture$refs))
    hit <- labels[vapply(labels, function(label) identical(get(label, envir = fixture$refs)$entity_id, entity_id), NA)]
    if (length(hit)) paste(hit, collapse = "|") else paste0("other:", substr(entity_id, 1L, 8L))
}
k273Known <- function(fixture, entity_id) !startsWith(k273Label(fixture, entity_id), "other:")

k273Commit <- function(fixture, state, object) {
    k273Quiet(fixture$manager$saveState(state, object, list(), "k273"))
    k273Name(fixture, state, k273HeadRef(fixture))
}

#' A DIA with limpa owner (protNormW1LimpaOwner's commits: protein_dpc_quant, protein_normalised, peptide_count_filtered),
#' with the DA's publication folder and the tab's DA state (test-workflow-saved-da.R's savedDaW1LimpaOwner).
k273LimpaOwner <- function(envir = parent.frame()) {
    withr::local_options(multischolar.dia_native_sessions = TRUE, multischolar.capture_test_checkpoints = FALSE,
        .local_envir = envir)
    root <- withr::local_tempdir(.local_envir = envir)
    paths <- list(source_dir = file.path(root, "scripts", "proteomics"),
        protein_qc_dir = file.path(root, "results", "proteomics", "protein_qc"),
        peptide_qc_dir = file.path(root, "results", "proteomics", "peptide_qc"),
        da_output_dir = file.path(root, "results", "proteomics", "da_proteins"))
    paths$publication_graphs_dir <- file.path(root, "results", "proteomics", "publication_graphs")
    paths$time_dir <- file.path(paths$publication_graphs_dir, "filtering_qc")
    for (path in unlist(paths)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(paths$source_dir, "k273-limpa", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager), envir = envir)
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("DIA_LIMPA")
    manager$.context$experiment_paths <- paths
    workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list()), config_list = list(),
        tab_status = list()))
    fixture <- list(root = root, paths = paths, workflow = workflow, manager = manager,
        da = new.env(parent = emptyenv()), refs = new.env(parent = emptyenv()), seen = new.env(parent = emptyenv()))
    setDIAStepProvenance(manager, "design_builder")
    k273Commit(fixture, "protein_dpc_quant", k273Quantified())
    k273Quiet(withDIAStateScope(manager, function() {
        v1 <- computeDIAReferenceProteinNormalisation(manager, getDIAStateReference(manager), list(norm_method = "none"))
        setDIAStepProvenance(manager, "protein_between_sample_normalisation", "whole", "norm_method", NULL,
            values = list(norm_method = "none"))
        finishDIAReferenceProteinStep(workflow, v1, publish = FALSE)
        k273Name(fixture, manager$current_state, getDIAStateReference(manager))
        v2 <- computeDIAReferenceProteinPeptideCountFilter(manager, getDIAStateReference(manager),
            list(num_peptides_per_protein_thresh = 1, num_peptidoforms_per_protein_thresh = 2))
        setDIAStepProvenance(manager, "protein_peptide_count_filter", "whole",
            c("num_peptides_per_protein_thresh", "num_peptidoforms_per_protein_thresh"), NULL, values = v2$config[1:2])
        finishDIAReferenceProteinStep(workflow, v2, publish = FALSE)
        k273Name(fixture, manager$current_state, getDIAStateReference(manager))
    }))
    fixture
}

#' A native proteomics owner of another route (kit244Owner, copied): the classic DIA, label-free and TMT routes.
k273Owner <- function(type, envir = parent.frame()) {
    withr::local_options(multischolar.dia_native_sessions = TRUE, multischolar.capture_test_checkpoints = FALSE,
        .local_envir = envir)
    root <- withr::local_tempdir(.local_envir = envir)
    paths <- list(source_dir = file.path(root, "scripts", "proteomics"),
        protein_qc_dir = file.path(root, "results", "proteomics", "protein_qc"),
        peptide_qc_dir = file.path(root, "results", "proteomics", "peptide_qc"))
    paths$publication_graphs_dir <- file.path(root, "results", "proteomics", "publication_graphs")
    paths$time_dir <- file.path(paths$publication_graphs_dir, "filtering_qc")
    for (path in unlist(paths)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    session <- new.env(parent = emptyenv())
    session$userData <- new.env(parent = emptyenv())
    session$onFlush <- session$onFlushed <- session$onSessionEnded <- function(...) function() invisible(NULL)
    original <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- WorkflowState$new(audit_enabled = FALSE)
    workflow$design_matrix <- original@design_matrix
    workflow$config_list <- original@args
    initializeDIANativeImport(workflow, paths, session, "diann", "iq")
    manager <- workflow$state_manager
    withr::defer(closeDIAStateAdapter(manager), envir = envir)
    manager$setWorkflowType(type)
    list(root = root, paths = paths, workflow = workflow, manager = manager)
}

k273W2Mocks <- function(env = parent.frame()) {
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
        finishWorkflowExecutionJob = function(owner, status) { owner$job <- NULL; list(status = status) },
        .package = "MultiScholaR", .env = env)
}

#' The limpa owner on a path; on the worker path the W2 mocks in the test's frame and a controller.
k273Start <- function(path, envir = parent.frame()) {
    fixture <- k273LimpaOwner(envir = envir)
    if (identical(path, "worker")) {
        k273W2Mocks(env = envir)
        fixture$seen$controller <- omicsWorkerController(list(owner = fixture$manager))
    }
    fixture
}

#' A failure as "class: message", NA when the code completed.
k273Error <- function(code) tryCatch({
    k273Quiet(code)
    NA_character_
}, error = function(e) paste0(class(e)[[1L]], ": ", conditionMessage(e)))

#' The S4 object a state holds, every slot hydrated (kit244Value's reading).
k273Object <- function(fixture, reference) {
    manager <- fixture$manager
    withDIAStateScope(manager, function() {
        object <- resolveDIAStateValue(manager$.context, reference)$data
        if (!isS4(object)) object <- rebuildDIACheckpointObject(manager$.context, object, TRUE, list())
        object
    })
}
k273Samples <- function(fixture, reference) {
    object <- k273Object(fixture, reference)
    setdiff(colnames(object@protein_quant_table), object@protein_id_column)
}

# ---- the tab's actions (kit244's, copied) ----

k273Press <- function(fixture, norm) {
    settings <- list(norm_method = "quantile", ruv_mode = "covariates", auto_percentage_min = NULL,
        auto_percentage_max = NULL, separation_metric = NULL, k_penalty_weight = NULL, adaptive_k_penalty = NULL,
        ruv_percentage = 50, ruv_k = 1, color_variable = "group", shape_variable = NULL)
    k273Error(runDIAReferenceLimpaNormalisationPassThrough(settings, fixture$workflow, norm, k273Noop))
}

k273Dispatch <- function(fixture) function(operation, parameters, input, state_name = NULL) {
    context <- if (!is.null(state_name)) list(state_name = state_name, config = list(), description = "dispatch")
    protNormW2Dispatch(fixture, fixture$seen$controller, operation, parameters, input, context)
}

k273WorkerCorrelation <- function(fixture, norm, threshold) {
    manager <- fixture$manager
    source <- withDIAStateScope(manager, function() {
        assertDIAExportReady(manager, allow_action = "correlation")
        source <- diaCorrelationSource(manager, norm)
        readDIAStateMetadata(manager, source)
        source
    })
    dispatch <- k273Dispatch(fixture)
    pairs <- dispatch("pearson-pairs", diaPearsonPairsParameters("group"), source)
    selection <- dispatch("pearson-selection", list(pairs = pairs$reference$entity_id, threshold = threshold), source,
        "correlation_filtered")
    withDIAStateScope(manager, function() publishDIAReferenceCorrelation(manager, list(pairs = pairs$reference,
        final = selection$value), fixture$workflow, norm, fixture$paths, threshold, FALSE, "proteomics", "core",
        k273Aesthetics))
}

#' Correlation Apply at `threshold`: its error, and the samples it kept.
k273Correlate <- function(fixture, norm, path, threshold) {
    final <- NULL
    error <- k273Error(final <- if (identical(path, "worker")) k273WorkerCorrelation(fixture, norm, threshold) else
        runDIAReferenceCorrelation(fixture$workflow, norm, fixture$paths, threshold, FALSE, "proteomics", "core",
            function() "group", k273Aesthetics))
    list(error = error, samples = if (is.null(final)) NULL else k273Samples(fixture, final))
}

k273RuvCovariates <- function(fixture, norm, path, mode) {
    manager <- fixture$manager
    settings <- list(ruv_mode = mode, color_variable = "group", shape_variable = NULL)
    if (!identical(path, "worker")) return(k273Error(runDIAReferenceRuvCovariates(settings, fixture$workflow, norm,
        fixture$paths, "proteomics", "core", k273Noop, function() "group")))
    k273Error({
        request <- withDIAStateScope(manager, function() {
            assertDIAExportReady(manager, allow_action = "ruv-covariates")
            list(source = diaRuvCovariatesSource(manager, norm),
                parameters = resolveDIARuvCovariatesParameters(fixture$workflow$config_list, "group", mode))
        })
        result <- k273Dispatch(fixture)("protein_ruv_covariates", request$parameters, request$source,
            if (identical(mode, "correction")) "ruv_corrected" else "ruv_covariates")
        withDIAStateScope(manager, function() publishDIAReferenceRuvCovariates(manager, settings,
            result[c("value", "reference")], fixture$workflow, norm, fixture$paths, k273Noop))
    })
}

k273RuvSkip <- function(fixture, norm) {
    norm$correlation_filtering_complete <- TRUE
    k273Quiet(runProtNormRuvCovariatesObserver(list(ruv_mode = "skip", color_variable = "group", shape_variable = NULL),
        new.env(), fixture$workflow, norm, fixture$paths, "proteomics", "core", k273Noop, function() "group",
        renderTextFn = function(x) x, showNotificationFn = k273Noop, messageFn = k273Noop))
    invisible(NULL)
}

#' t_hi by KIT-244's G0 rule (kit244ThresholdHigh, copied): a threshold that removes exactly one sample the filter at 0.5
#' keeps. NULL when the fixture has none.
k273ThresholdHigh <- function() {
    if (!is.null(k273Cache$t_hi)) return(k273Cache$t_hi)
    fixture <- k273Start("in process", envir = environment())
    norm <- new.env(parent = emptyenv())
    k273Press(fixture, norm)
    first <- k273Correlate(fixture, norm, "in process", 0.5)
    hydrated <- hydrateSavedWorkflowEntity(fixture$manager$.context$store, norm$correlation_vector)
    pairs <- as.data.frame(if (is.data.frame(hydrated)) hydrated else hydrated$value)
    x <- grep("[.]x$", names(pairs), value = TRUE)[[1L]]
    y <- sub("[.]x$", ".y", x)
    in_pairs <- unique(c(as.character(pairs[[x]]), as.character(pairs[[y]])))
    best <- vapply(in_pairs, function(s) {
        r <- pairs$pearson_correlation[pairs[[x]] == s | pairs[[y]] == s]
        r <- r[!is.na(r)]
        if (length(r)) max(r) else NA_real_
    }, numeric(1))
    candidates <- best[!is.na(best) & best >= 0.5]
    if (!length(candidates)) return(NULL)
    m <- min(candidates)
    above <- best[!is.na(best) & best > m]
    if (!length(above)) return(NULL)
    k273Cache$samples_05 <- first$samples
    k273Cache$removed_hi <- names(best)[!is.na(best) & best == m]
    k273Cache$t_hi <- (m + min(above)) / 2
    k273Cache$t_hi
}

# ---- the DA, Copy and what they leave ----

k273Contrasts <- function() data.frame(contrasts = "treated.vs.control",
    full_format = "treated.vs.control=grouptreated-groupcontrol", friendly_names = "treated.vs.control",
    stringsAsFactors = FALSE)

#' The tab's DA click (savedDaW1Click, copied): one contrast, in process.
k273Da <- function(fixture) {
    contrasts <- k273Contrasts()
    fixture$workflow$contrasts_tbl <- contrasts
    fixture$da$contrasts_available <- contrasts$contrasts
    fixture$da$current_s4_object <- k273HeadRef(fixture)
    settings <- list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)
    error <- NA_character_
    utils::capture.output(error <- k273Error(runDIAReferenceDaAnalysis(settings, fixture$da, fixture$workflow,
        fixture$paths, contrasts, fixture$workflow$uniprot_dat_cln)))
    error
}

#' Copy's paths beside the project's own (k271CopyPaths, copied), the summary folders under `root/<tag>`.
k273CopyPaths <- function(fixture, tag) {
    paths <- fixture$paths
    base <- file.path(fixture$root, tag)
    defaults <- list(results_dir = file.path(base, "results"), results_summary_dir = file.path(base, "results_summary"),
        time_dir = file.path(base, "time"), qc_dir = file.path(base, "qc"), pathway_dir = file.path(base, "pathway"),
        da_output_dir = file.path(base, "da"), source_dir = file.path(base, "source"))
    for (name in names(defaults)) if (is.null(paths[[name]])) paths[[name]] <- defaults[[name]]
    paths$results_summary_dir <- defaults$results_summary_dir
    paths$feature_qc_dir <- fixture$paths$protein_qc_dir
    paths
}

#' Copy (`copyToResultsSummary`) as the tab's Summary runs it: inside the owner's state scope. Its `COPY:` lines, its
#' failures (type and source) and what reached the summary folders.
k273Copy <- function(fixture, tag, omic = "proteomics") {
    paths <- k273CopyPaths(fixture, tag)
    for (path in unlist(paths)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    withDIAStateScope(fixture$manager, function() {
        testthat::local_mocked_bindings(getProjectPaths = function(...) paths, .package = "MultiScholaR")
        failures <- NULL
        console <- utils::capture.output(failures <- suppressWarnings(suppressMessages(copyToResultsSummary(omic, "k273",
            contrasts_tbl = data.frame(contrasts = "groupB-groupA"), design_matrix = data.frame(Run = "S1", group = "A"),
            force = TRUE))))
        summary <- paths$results_summary_dir
        list(console = grep("^COPY:", console, value = TRUE), summary = summary,
            failed = vapply(failures, function(failure) basename(as.character(failure$source %||% "")), ""),
            types = vapply(failures, function(failure) as.character(failure$type %||% ""), ""),
            tables = sort(list.files(file.path(summary, "Publication_tables"))),
            study = sort(list.files(file.path(summary, "Study_report"))))
    })
}

k273Delivered <- function(copy, name) file.path(copy$summary, name)

#' The DA's own record and source: the manifest, its source state's object, and E, the EList the model fitted.
k273DaSource <- function(fixture) {
    manager <- fixture$manager
    manifest <- withDIAStateScope(manager, function()
        resolveWorkflowStateValue(manager$.context, manager$.context$snapshot$da_results))
    object <- k273Object(fixture, manifest$source)
    list(manifest = manifest, object = object, E = object@args$limpa_dpc_quant_results$quantified_elist$E)
}

#' A frame as the delivered text: written by the package's writer and read back as strings.
k273Text <- function(frame, writer = vroom::vroom_write) {
    path <- tempfile(fileext = ".tsv")
    on.exit(unlink(path), add = TRUE)
    writer(frame, path)
    utils::read.delim(path, check.names = FALSE, colClasses = "character", na.strings = character())
}
k273Read <- function(path) utils::read.delim(path, check.names = FALSE, colClasses = "character", na.strings = character())

#' E as the text a table of the package holds it (`writer` is the writer that table is written with: vroom for the RUV
#' files, readr for the correlation filter's table), ids in the first column named `id_column`.
k273EText <- function(E, id_column, writer = vroom::vroom_write) {
    frame <- data.frame(rownames(E), as.data.frame(E, check.names = FALSE), check.names = FALSE, stringsAsFactors = FALSE)
    names(frame)[[1L]] <- id_column
    k273Text(frame, writer)
}

#' Whether a delivered matrix (or the pre-RUV table) holds E: the protein ids and the sample columns are the same sets
#' (compared first), then every value matches by id and sample on the delivered text. Returns the findings.
k273CompareToE <- function(delivered, E, id_column, whole_table = TRUE, writer = vroom::vroom_write) {
    expected <- k273EText(E, id_column, writer)
    ids_delivered <- delivered[[id_column]]
    samples <- colnames(E)
    out <- list(ids = setequal(ids_delivered, rownames(E)) && !anyDuplicated(ids_delivered),
        samples = if (whole_table) setequal(setdiff(names(delivered), id_column), samples) else all(samples %in% names(delivered)))
    if (!out$ids || !out$samples) return(c(out, list(values = FALSE)))
    delivered <- delivered[match(rownames(E), delivered[[id_column]]), c(id_column, samples), drop = FALSE]
    out$values <- identical(unname(as.matrix(delivered[samples])), unname(as.matrix(expected[samples])))
    out
}

#' One project: the limpa owner, the pass-through, the correlation filter (Apply at `threshold`) and, with a `mode`, the RUV
#' covariates; setup failures are returned and asserted by the test.
k273Project <- function(path, mode, threshold = 0.5, envir = parent.frame()) {
    fixture <- k273Start(path, envir = envir)
    norm <- new.env(parent = emptyenv())
    errors <- c(press = k273Press(fixture, norm))
    errors[["correlation"]] <- k273Correlate(fixture, norm, path, threshold)$error
    if (!is.null(mode)) errors[["ruv"]] <- k273RuvCovariates(fixture, norm, path, mode)
    k273Log("setup", path, mode %||% "no RUV", "errors:", paste(errors, collapse = " / "), "state:",
        fixture$manager$current_state)
    list(fixture = fixture, norm = norm, path = path, mode = mode, errors = errors)
}

k273ExpectSetup <- function(project, info) {
    expect_identical(unname(project$errors), rep(NA_character_, length(project$errors)), info = info)
}

# ---- T1, T2b, T3, T7, T8, T9, T13: the covariates mode ----

test_that("KIT-273 T1, T2b, T3, T7, T8, T9, T13: covariates mode: the design with the factors is delivered, from the current DA's record", {
    k273Need()
    for (path in k273Paths()) {
        info <- paste(path, "covariates")
        project <- k273Project(path, "covariates")
        fixture <- project$fixture
        k273ExpectSetup(project, info)
        expect_identical(fixture$manager$current_state, "ruv_covariates", info = info)
        qc_design <- file.path(fixture$paths$protein_qc_dir, "ruv_covariates_design.tsv")
        expect_true(file.exists(qc_design), info = info)
        # T13: no differential analysis in the store: neither file.
        before <- k273Copy(fixture, "before_da")
        k273Log("T13 covariates", path, "tables:", paste(before$tables, collapse = ","), "study:", paste(before$study, collapse = ","))
        expect_false(file.exists(k273Delivered(before, k273Design)), info = paste("T13", info))
        expect_false(file.exists(k273Delivered(before, k273Matrix)), info = paste("T13", info))
        expect_identical(k273Da(fixture), NA_character_, info = info)
        expect_identical(workflowDaResultStatus(fixture$manager), "current", info = info)
        copy <- k273Copy(fixture, "copy1")
        k273Log("T1 covariates", path, "COPY:", paste(copy$console, collapse = " | "), "tables:", paste(copy$tables, collapse = ","),
            "study:", paste(copy$study, collapse = ","), "failed:", paste(copy$failed, collapse = ","))
        delivered <- k273Delivered(copy, k273Design)
        # T1: the design is delivered and holds the factors the model used.
        expect_true(file.exists(delivered), info = paste("T1", info))
        if (file.exists(delivered)) {
            source <- k273DaSource(fixture)
            table <- k273Read(delivered)
            w_columns <- grep("^W[0-9]+$", names(table), value = TRUE)
            expect_gte(length(w_columns), 1L)
            expect_identical(w_columns, source$manifest$model$w_terms, info = paste("T1", info))
            expect_identical(names(table), names(source$object@design_matrix), info = paste("T1 whole design", info))
            expected <- k273Text(source$object@design_matrix)
            expect_identical(table[w_columns], expected[w_columns], info = paste("T1 W values", info))
            k273Log("T1 covariates", path, "W terms:", paste(w_columns, collapse = ","), "digest:", k273Sha(delivered))
            # T8: the same bytes as the protein-QC file.
            expect_identical(k273Sha(delivered), k273Sha(qc_design), info = paste("T8", info))
        }
        # T3: the covariates mode delivers no matrix.
        expect_false(file.exists(k273Delivered(copy, k273Matrix)), info = paste("T3", info))
        # T2b: the pre-RUV table (KIT-271's name) is the DA's E in this mode.
        source <- k273DaSource(fixture)
        norm_table <- k273Delivered(copy, "Publication_tables/normalised_results.tsv")
        expect_true(file.exists(norm_table), info = paste("T2b", info))
        if (file.exists(norm_table)) {
            found <- k273CompareToE(k273Read(norm_table), source$E, source$object@protein_id_column, whole_table = FALSE,
                writer = function(x, path) readr::write_tsv(x, path))
            k273Log("T2b", path, "ids:", found$ids, "samples:", found$samples, "values:", found$values)
            expect_true(found$ids, info = paste("T2b", info))
            expect_true(found$samples, info = paste("T2b", info))
            expect_true(found$values, info = paste("T2b", info))
        }
        # T9: the existing deliveries are as before.
        expect_true("normalised_results.tsv" %in% copy$tables, info = paste("T9", info))
        expect_false(any(c("RUV_normalised_results.tsv", "ruv_normalised_results.RDS") %in% copy$tables), info = paste("T9", info))
        expect_true(all(c("contrasts_tbl.tab", "design_matrix.tab") %in% copy$study), info = paste("T9", info))
        expect_identical(setdiff(copy$study, "design_matrix_with_ruv_factors.tab"),
            setdiff(before$study, "design_matrix_with_ruv_factors.tab"), info = paste("T9", info))
        design_tab <- k273Read(file.path(copy$summary, "Study_report", "design_matrix.tab"))
        expect_identical(names(design_tab), c("Run", "group"), info = paste("T9", info))
        # T7: the protein-QC files gone (what an archive restore gives): the delivery is unchanged.
        file.remove(file.path(fixture$paths$protein_qc_dir, c("ruv_covariates_design.tsv", "ruv_adjusted_protein_matrix.tsv")))
        again <- k273Copy(fixture, "copy2")
        expect_true(file.exists(k273Delivered(again, k273Design)), info = paste("T7", info))
        if (file.exists(k273Delivered(again, k273Design)) && file.exists(delivered)) {
            expect_identical(k273Sha(k273Delivered(again, k273Design)), k273Sha(delivered), info = paste("T7", info))
        }
    }
})

# ---- T2, T3, T7, T8, T13: the correction mode ----

test_that("KIT-273 T2, T3, T7, T8, T13: correction mode: the corrected matrix is delivered, from the current DA's record", {
    k273Need()
    for (path in k273Paths()) {
        info <- paste(path, "correction")
        project <- k273Project(path, "correction")
        fixture <- project$fixture
        k273ExpectSetup(project, info)
        expect_identical(fixture$manager$current_state, "ruv_corrected", info = info)
        qc_matrix <- file.path(fixture$paths$protein_qc_dir, "ruv_adjusted_protein_matrix.tsv")
        expect_true(file.exists(qc_matrix), info = info)
        before <- k273Copy(fixture, "before_da")
        expect_false(file.exists(k273Delivered(before, k273Matrix)), info = paste("T13", info))
        expect_false(file.exists(k273Delivered(before, k273Design)), info = paste("T13", info))
        expect_identical(k273Da(fixture), NA_character_, info = info)
        expect_identical(workflowDaResultStatus(fixture$manager), "current", info = info)
        copy <- k273Copy(fixture, "copy1")
        k273Log("T2 correction", path, "COPY:", paste(copy$console, collapse = " | "), "tables:", paste(copy$tables, collapse = ","),
            "study:", paste(copy$study, collapse = ","), "failed:", paste(copy$failed, collapse = ","))
        delivered <- k273Delivered(copy, k273Matrix)
        expect_true(file.exists(delivered), info = paste("T2", info))
        if (file.exists(delivered)) {
            source <- k273DaSource(fixture)
            expect_identical(source$manifest$model$ruv_mode, "correction", info = info)
            found <- k273CompareToE(k273Read(delivered), source$E, source$object@protein_id_column)
            k273Log("T2 correction", path, "ids:", found$ids, "samples:", found$samples, "values:", found$values,
                "digest:", k273Sha(delivered))
            expect_true(found$ids, info = paste("T2 protein ids", info))
            expect_true(found$samples, info = paste("T2 samples", info))
            expect_true(found$values, info = paste("T2 values", info))
            # T8: the same bytes as the protein-QC file.
            expect_identical(k273Sha(delivered), k273Sha(qc_matrix), info = paste("T8", info))
        }
        # T3: the correction mode delivers no design with the factors.
        expect_false(file.exists(k273Delivered(copy, k273Design)), info = paste("T3", info))
        # T9: the existing deliveries are as before.
        expect_true("normalised_results.tsv" %in% copy$tables, info = paste("T9", info))
        expect_true(all(c("contrasts_tbl.tab", "design_matrix.tab") %in% copy$study), info = paste("T9", info))
        # T7: the protein-QC files gone: the delivery is unchanged.
        file.remove(file.path(fixture$paths$protein_qc_dir, c("ruv_covariates_design.tsv", "ruv_adjusted_protein_matrix.tsv")))
        again <- k273Copy(fixture, "copy2")
        expect_true(file.exists(k273Delivered(again, k273Matrix)), info = paste("T7", info))
        if (file.exists(k273Delivered(again, k273Matrix)) && file.exists(delivered)) {
            expect_identical(k273Sha(k273Delivered(again, k273Matrix)), k273Sha(delivered), info = paste("T7", info))
        }
    }
})

# ---- T4: RUV skipped ----

test_that("KIT-273 T4: RUV skipped: neither file is delivered (guard)", {
    k273Need()
    project <- k273Project("in process", NULL)
    fixture <- project$fixture
    k273ExpectSetup(project, "skip")
    k273RuvSkip(fixture, project$norm)
    expect_identical(fixture$manager$current_state, "correlation_filtered")
    expect_identical(k273Da(fixture), NA_character_)
    source <- k273DaSource(fixture)
    k273Log("T4", "DA record: model", source$manifest$model$model, "ruv_mode", source$manifest$model$ruv_mode, "W terms:",
        paste(source$manifest$model$w_terms, collapse = ","))
    copy <- k273Copy(fixture, "copy1")
    expect_false(file.exists(k273Delivered(copy, k273Matrix)))
    expect_false(file.exists(k273Delivered(copy, k273Design)))
    expect_true("normalised_results.tsv" %in% copy$tables)
})

# ---- T5: stale protein-QC files ----

test_that("KIT-273 T5: a covariates run, then the correlation filter again and a DA on its table: the earlier run's files are not delivered (guard)", {
    k273Need()
    project <- k273Project("in process", "covariates")
    fixture <- project$fixture
    k273ExpectSetup(project, "T5")
    expect_identical(k273Press(fixture, project$norm), NA_character_)
    expect_identical(k273Correlate(fixture, project$norm, "in process", 0.5)$error, NA_character_)
    expect_identical(fixture$manager$current_state, "correlation_filtered")
    qc <- file.path(fixture$paths$protein_qc_dir, c("ruv_covariates_design.tsv", "ruv_adjusted_protein_matrix.tsv"))
    expect_true(all(file.exists(qc)))
    expect_identical(k273Da(fixture), NA_character_)
    source <- k273DaSource(fixture)
    k273Log("T5", "DA record: ruv_mode", source$manifest$model$ruv_mode, "W terms:", paste(source$manifest$model$w_terms, collapse = ","))
    copy <- k273Copy(fixture, "copy1")
    expect_false(file.exists(k273Delivered(copy, k273Matrix)))
    expect_false(file.exists(k273Delivered(copy, k273Design)))
})

# ---- T6: a mode switch ----

test_that("KIT-273 T6: a correction run, then a covariates run, then the DA: the design is delivered, the matrix is not", {
    k273Need()
    project <- k273Project("in process", "correction")
    fixture <- project$fixture
    k273ExpectSetup(project, "T6")
    expect_identical(k273RuvCovariates(fixture, project$norm, "in process", "covariates"), NA_character_)
    expect_identical(fixture$manager$current_state, "ruv_covariates")
    # The protein-QC matrix file now holds the quality-control matrix of the covariates run, which the DA does not use.
    expect_true(file.exists(file.path(fixture$paths$protein_qc_dir, "ruv_adjusted_protein_matrix.tsv")))
    expect_identical(k273Da(fixture), NA_character_)
    copy <- k273Copy(fixture, "copy1")
    k273Log("T6", "tables:", paste(copy$tables, collapse = ","), "study:", paste(copy$study, collapse = ","))
    expect_true(file.exists(k273Delivered(copy, k273Design)))
    expect_false(file.exists(k273Delivered(copy, k273Matrix)))
})

# ---- T10: the other routes ----

test_that("KIT-273 T10: classic DIA, label-free and TMT owners: Copy lists no new file (guard)", {
    for (type in c("DIA", "LFQ", "TMT")) {
        fixture <- k273Owner(type)
        # Not a native session: Copy then writes no session archive (none was exported on this fixture).
        fixture$manager$.context$session_format <- "none"
        copy <- k273Copy(fixture, "copy1")
        k273Log("T10", type, "tables:", paste(copy$tables, collapse = ","), "study:", paste(copy$study, collapse = ","))
        expect_false(file.exists(k273Delivered(copy, k273Matrix)), info = type)
        expect_false(file.exists(k273Delivered(copy, k273Design)), info = type)
        expect_false("ruv_delivery" %in% copy$types, info = type)
    }
})

# ---- T11: a failure to write ----

test_that("KIT-273 T11: the writer fails: Copy records a ruv_delivery failure, still delivers the others, and no error escapes (the writer is mocked to throw)", {
    k273Need()
    project <- k273Project("in process", "covariates")
    fixture <- project$fixture
    k273ExpectSetup(project, "T11")
    expect_identical(k273Da(fixture), NA_character_)
    testthat::local_mocked_bindings(writeDIADesignTable = function(...) stop("kit273 injected write failure"),
        .package = "MultiScholaR")
    copy <- NULL
    expect_no_error(copy <- k273Copy(fixture, "copy1"))
    k273Log("T11", "types:", paste(copy$types, collapse = ","), "failed:", paste(copy$failed, collapse = ","), "tables:",
        paste(copy$tables, collapse = ","))
    expect_true("ruv_delivery" %in% copy$types)
    expect_false(file.exists(k273Delivered(copy, k273Design)))
    expect_true("normalised_results.tsv" %in% copy$tables)
    expect_true("contrasts_tbl.tab" %in% copy$study)
})

# ---- T12: a sample the correlation filter dropped ----

test_that("KIT-273 T12: the correction matrix has no column for a sample the correlation filter dropped", {
    k273Need()
    t_hi <- k273ThresholdHigh()
    skip_if(is.null(t_hi), "no t_hi on this fixture")
    expect_length(k273Cache$removed_hi, 1L)
    expect_true(all(k273Cache$removed_hi %in% k273Cache$samples_05))
    project <- k273Project("in process", "correction", threshold = t_hi)
    fixture <- project$fixture
    k273ExpectSetup(project, "T12")
    expect_identical(k273Da(fixture), NA_character_)
    copy <- k273Copy(fixture, "copy1")
    delivered <- k273Delivered(copy, k273Matrix)
    expect_true(file.exists(delivered))
    if (file.exists(delivered)) {
        columns <- names(k273Read(delivered))
        k273Log("T12", "removed:", paste(k273Cache$removed_hi, collapse = ","), "delivered columns:", length(columns))
        expect_false(any(k273Cache$removed_hi %in% columns))
        expect_true(all(setdiff(k273Cache$samples_05, k273Cache$removed_hi) %in% columns))
    }
})
