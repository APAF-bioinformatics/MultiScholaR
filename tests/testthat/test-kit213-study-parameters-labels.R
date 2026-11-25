# KIT-213 (VAL-169 (1)): in study_parameters.txt's analysis-record section, a parameter that holds an earlier result's
# identifier in the session's store (`pairs`, `logcpm`, `grid`, `results`, provenance `input`) is printed as a label naming
# the step that produced it, so two projects that run the same analysis write the same file apart from the run-bound
# lines (E1's masks). Written before the fix: on the tests-only commit E1 to E4 fail at their base values (the differing
# masked lines are exactly the L lines; the label lines are absent) and so does G5 (the L lines hold store ids); G3 and
# G4 pass on both commits. The step numbers in E3 and E4 are those G0 printed (run 45471271185). Lines starting
# `KIT-213` are the evidence the ticket's records copy (G2 compares each line's SHA-256 base against head).

library(testthat)

# ------------------------------------------------------------------------------------------------
# Copied unchanged (KIT-207 E-7 form) from tests/testthat/test-kit207-summary-head-guard.R at e305fab3: lines 16-29,
# 97-124, 126-127, 129-144, 148-175, 238-256, 279-284, 288-309 and 352-370; and from
# tests/testthat/test-workflow-analysis-record.R at e305fab3: lines 3-28 (`analysisRecordFixture`) and 955-968
# (`recordEnrichmentConfig`). Those files are not edited.
# ------------------------------------------------------------------------------------------------

k194Omic <- function(class) c(lipid = "lipidomics", metab = "metabolomics")[[class]]

#' Copy of `kit172PreNorm` (test-kit172-downstream.R:127-137; A11): an owner at the QC-complete state (intensity filter,
#' QC finalise), as the normalisation tab finds it.
k194PreNorm <- function(root, values, class) {
    fixture <- kit172Owner(root, kit172Object(values, class), class)
    binding <- kit172Binding(class)
    kit172Quiet({
        runWorkflowOmicsReferenceIntensityFilter(binding, fixture$owner, 1, 0.5, fixture$workflow$config_list)
        commitWorkflowOmicsStoreStep(fixture$owner, paste0(binding$prefix, "_qc_finalize"), paste0(binding$prefix, "_qc_complete"),
            fixture$workflow$config_list, "QC complete")
    })
    fixture
}

#' A lipid or metab project at QC complete on F172 (seed 172001), with the omic's GUI configuration (its QC figures
#' drawn by nothing) and the tab's log.
k194Fixture <- function(root, class) {
    fixture <- k194PreNorm(root, kit172Values(172001L)$values, class)
    fixture$class <- class
    fixture$omic <- k194Omic(class)
    fixture$binding <- kit172Binding(class)
    fixture$prefix <- fixture$binding$prefix
    config <- omicGuiConfig(fixture$omic)
    config$binding$binding$qc_plots <- function(...) invisible(NULL)
    fixture$config <- config
    log <- new.env(parent = emptyenv())
    log$lines <- character()
    fixture$log <- log
    fixture$addLog <- function(entry) log$lines <- c(log$lines, entry)
    fixture
}

k194Close <- function(fixture) try(closeDIAStateAdapter(fixture$workflow$state_manager), silent = TRUE)

#' Close the bus's owner of the moment (a load replaces it) when the calling test ends; `fixture` is forced here, so a
#' loop's fixtures are each closed.
k194Defer <- function(fixture, envir = parent.frame()) {
    force(fixture)
    withr::defer(k194Close(fixture), envir = envir)
}

k194Owner <- function(fixture) fixture$workflow$state_manager
k194Head <- function(fixture) k194Owner(fixture)$current_state
k194Ref <- function(fixture, state = NULL) getWorkflowStateReference(k194Owner(fixture), state)

#' Shiny's notifications recorded as "[type] text"; removal, progress and selector updates quiet.
k194Notes <- function(env = parent.frame()) {
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    testthat::local_mocked_bindings(
        showNotification = function(ui, ..., type = "default") {
            seen$notes <- c(seen$notes, paste0("[", type, "] ", paste(as.character(ui), collapse = " ")))
            invisible(NULL)
        },
        removeNotification = function(...) invisible(NULL),
        withProgress = function(expr, ...) expr,
        incProgress = function(...) invisible(NULL),
        updateSelectInput = function(...) invisible(NULL),
        .package = "shiny", .env = env)
    seen
}

#' The tab's inputs: ITSD off unless asked (by pattern), offset 0, `method`, RUV skip unless `ruv` (the legacy list).
k194Input <- function(method = "quantile", itsd = FALSE, ruv = NULL, threshold = 0.5) {
    list(apply_itsd = itsd, itsd_aggregation = "median", log_offset = 0, norm_method = method,
        ruv_mode = if (is.null(ruv)) "skip" else ruv$ruv_mode, auto_percentage_min = 1, auto_percentage_max = 20,
        ruv_grouping_variable = "group", separation_metric = "max_difference", k_penalty_weight = 0.5,
        adaptive_k_penalty = TRUE, ruv_k = ruv$manual_k %||% 1L, ruv_percentage = ruv$manual_percentage %||% 20,
        min_pearson_correlation_threshold = threshold)
}

#' "Run Normalization Pipeline": the frame's `runOmicNormPipeline`, in process (`runner` NULL) or managed.
k194Run <- function(fixture, nd, input = k194Input(), runner = NULL) {
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) runner, .package = "MultiScholaR")
    kit172Quiet(shiny::isolate(runOmicNormPipeline(input, fixture$workflow, list(), nd, fixture$addLog,
        function(...) list(color_var = "group", shape_var = NULL), fixture$config)))
}

#' The run a case builds on completed (an instrument check, VAL-54): head the normalised or RUV state, flag set.
k194ExpectRun <- function(fixture, nd, info = "") {
    expect_true(k194Head(fixture) %in% paste0(fixture$prefix, c("_normalized", "_ruv_corrected")),
        info = paste("run precondition:", info, k194Head(fixture), paste(fixture$log$lines, collapse = " | ")))
    expect_true(isTRUE(shiny::isolate(nd$normalization_complete)), info = paste("run precondition:", info))
}

k194Apply <- function(fixture, nd, runner = NULL, threshold = 0.5) {
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) runner, .package = "MultiScholaR")
    kit172Quiet(shiny::isolate(applyOmicNormCorrelation(k194Input(threshold = threshold), fixture$workflow, nd, fixture$addLog,
        fixture$config)))
}
#' trx's "Run Normalization Pipeline" (`runTrxNormalisation`, as test-trx-norm-core.R drives it) with the notifications
#' recorded in `seen$notes` and the tab's log in `nd`.
k194TrxRun <- function(project, nd, inputs, seen, runner = NULL) {
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) runner, .package = "MultiScholaR")
    notify <- function(ui, ..., type = "default", duration = 5) {
        seen$notes <- c(seen$notes, paste0("[", type, "] ", paste(as.character(ui), collapse = " ")))
        invisible(NULL)
    }
    set_status <- function(tab, value) {
        status <- project$workflow$tab_status
        status[[tab]] <- value
        project$workflow$tab_status <- status
    }
    withr::local_seed(20261124L)
    kit172Quiet(utils::capture.output(shiny::isolate(runTrxNormalisation(list(managed = runTrxNormManaged, in_process = runTrxNormInProcess),
        transcriptomicsDeclarationBinding, project$workflow$state_manager, trxNormGridParameters(inputs), nd, project$workflow$config_list,
        project$paths, omicNormAddLog(nd), notify, set_status))))
    invisible(NULL)
}

#' A fresh project's directories for the omic (the Summary's `projectDirs`).
k207Dirs <- function(root, omic) {
    paths <- list(source_dir = file.path(root, "scripts", omic), integration_dir = file.path(root, "integration"), base_dir = root)
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    stats::setNames(list(paths), omic)
}

#' "Save Workflow Arguments": the core's `saveWorkflowSummaryParameters` under the omic's group, in process (`runner`
#' NULL) or managed, with the real study-parameter and record writers; the notifications ("[type] text"), `done`,
#' `failed`, the in-process RDS writes and the console lines recorded.
k207Save <- function(workflow, omic, dirs, runner = NULL) {
    rec <- new.env(parent = emptyenv())
    rec$notes <- character()
    rec$done <- 0L
    rec$failed <- list()
    rec$rds <- list()
    rec$cat <- character()
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) runner, .package = "MultiScholaR")
    kit172Quiet(utils::capture.output(saveWorkflowSummaryParameters(omicGuiConfig(omic)$summary, workflow, dirs, omic, "k207", "",
        notify = function(message, type, duration) rec$notes <- c(rec$notes, paste0("[", type, "] ", message)),
        done = function(saved) rec$done <- rec$done + 1L,
        failed = function(error) rec$failed[[length(rec$failed) + 1L]] <- error,
        saveRDSFn = function(object, path) {
            rec$rds[[length(rec$rds) + 1L]] <- list(object = object, path = path)
            invisible(NULL)
        },
        catFn = function(...) rec$cat <- c(rec$cat, paste(..., collapse = " ")))))
    rec
}

#' The reference DA on the head (as `kit172Journey`); it commits no state.
k207Da <- function(fx) {
    kit172Quiet(runWorkflowOmicsReferenceDaAnalysis(fx$binding, k194Owner(fx), k194Ref(fx),
        list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), kit172Contrasts()))
    invisible(fx)
}

#' A lipid or metab fixture after "Run Normalization Pipeline" (`input`), Apply (unless `apply` is FALSE) and the DA
#' (when `da`); closed when the calling test ends.
k207Normalised <- function(root, class, input = k194Input(), apply = TRUE, da = FALSE, env = parent.frame()) {
    fx <- k194Fixture(root, class)
    k194Defer(fx, envir = env)
    nd <- omicNormReactiveState()
    k194Run(fx, nd, input)
    k194ExpectRun(fx, nd, paste("KIT-207", class, basename(root)))
    if (apply) expect_true(isTRUE(k194Apply(fx, nd)), info = paste("apply precondition", class, basename(root)))
    if (da) k207Da(fx)
    list(fx = fx, nd = nd)
}

analysisRecordFixture <- function(root) {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "record-test", "DIA"), WorkflowState$new(audit_enabled = FALSE))
    manager$setWorkflowType("DIA")
    manager$.context$reference_states <- TRUE
    peptide <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    peptide@args$filterMinNumPeptidesPerSample <- list(peptides_per_sample_cutoff = 500,
        inclusion_list = paste0("sample_", seq_len(40)))
    manager$saveState("raw_data_s4", peptide, list(), "imported")
    filtered <- peptide
    filtered@args$srlQvalueProteotypicPeptideClean <- list(qvalue_threshold = 0.05, global_qvalue_threshold = 0.01,
        global_pg_qvalue_threshold = 0.01, confidence_provenance_mode = "diann_main_report")
    entry <- workflowOperationEntry("peptide_qvalue_filter")
    setDIAStepProvenance(manager, "peptide_qvalue_filter", "whole", entry$parameters)
    manager$saveState("qvalue_filtered", filtered, list(proteotypic_only = TRUE), "q-value filter")
    fasta <- file.path(root, "proteins.fasta")
    writeLines(c(">sp|P1|ONE", "MPEPTIDE"), fasta)
    recordDIAWorkflowInputs(manager, c(fasta = fasta))
    protein <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    protein@args$chooseBestProteinAccession <- list(delim = ";", aggregation_method = "mean")
    entry <- workflowOperationEntry("protein_accession_cleanup")
    setDIAStepProvenance(manager, "protein_accession_cleanup", "whole", entry$parameters)
    manager$saveState("protein_accession_cleaned", protein, list(delimiter = ";", aggregation_method = "mean",
        fasta_sha256 = digest::digest(file = fasta, algo = "sha256")), "accession cleanup")
    recordWorkflowExternalResource(manager, "string-db", "12.0", "api_version", query_summary = "fixture")
    list(manager = manager, peptide = peptide, fasta = fasta)
}

#' The proteomics `enrichment` group with its function fields (KIT-107 D1): until the W2 merge wires them in
#' utils_omic_gui_config.R they are optional; this fills an absent one with the value the merge wires, and changes
#' nothing once they are wired.
recordEnrichmentConfig <- function(config) {
    config[["organism"]] <- config[["organism"]] %||% list(supported = buildProtEnrichSupportedOrganisms,
        method = resolveProtEnrichAnalysisMethod)
    config$annotation <- config$annotation %||% list(resolve = resolveProtEnrichUniprotAnnotations,
        match = resolveProtEnrichAnnotationMatching)
    config$kernel <- config$kernel %||% list(select = resolveProtEnrichProcessEnrichmentsFn,
        test_double = runProtEnrichDeterministicProcessEnrichments)
    config$organism_filter$mapping <- config$organism_filter$mapping %||% resolveProtEnrichOrganismMapping
    config$organism_filter$apply <- config$organism_filter$apply %||% applyProtEnrichOrganismFilter
    config
}

# ------------------------------------------------------------------------------------------------
# KIT-213's harness.
# ------------------------------------------------------------------------------------------------

# An L line: a record parameter that holds an earlier result of the analysis (ticket 2.2).
k213LPattern <- "^  ([^ ]+\\.)?(pairs|logcpm|grid|results) \\[input\\] = "

#' E1's masks: the lines holding `Timestamp:`, starting `  timestamp = ` (the proteomics DA settings' time) or naming the
#' version (repository, branch, commit, package version, source) dropped; the project roots masked.
k213Masked <- function(lines, roots) {
    lines <- lines[!grepl("Timestamp:", lines, fixed = TRUE) & !grepl("^  timestamp = ", lines) &
        !grepl("^(Repository|Branch|Commit|Package Version|Source):", lines)]
    roots <- unique(c(normalizePath(roots, mustWork = FALSE), roots))
    for (root in roots[order(-nchar(roots))]) lines <- gsub(root, "<root>", lines, fixed = TRUE)
    lines
}

k213Sha <- function(lines) digest::digest(paste(lines, collapse = "\n"), algo = "sha256", serialize = FALSE)

#' The store's incarnation and its entity ids: the store's index, and every `entity_id` under the head's states,
#' products, exports and DA results and under the DA manifest (G5).
k213StoreIds <- function(owner) {
    context <- owner$.context
    walk <- function(x) {
        if (!is.list(x)) return(character())
        c(if (is.character(x$entity_id)) x$entity_id, unlist(lapply(unname(x), walk), use.names = FALSE))
    }
    snapshot <- context$snapshot
    manifest <- tryCatch(resolveWorkflowStateValue(context, snapshot$da_results), error = function(e) NULL)
    ids <- c(ls(envir = context$store$records, all.names = TRUE),
        walk(snapshot[intersect(c("states", "products", "exports", "da_results"), names(snapshot))]), walk(manifest))
    list(incarnation = context$store$identity$incarnation, ids = unique(ids[grepl("^[0-9a-f]{64}$", ids)]))
}

#' One project's result: the study parameters' lines, the masks' roots, the store's ids, the save's record.
k213Result <- function(dirs, omic, owner, roots, rec) {
    file <- file.path(dirs[[omic]]$source_dir, "study_parameters.txt")
    list(lines = if (file.exists(file)) readLines(file, warn = FALSE) else character(), written = file.exists(file),
        roots = roots, ids = k213StoreIds(owner), done = rec$done, notes = paste(rec$notes, collapse = " | "))
}

#' Lipidomics or metabolomics in a fresh root: KIT-207 G1's journey (run, Apply, DA), then Save.
k213FeatureJourney <- function(class) {
    root <- withr::local_tempdir()
    seen <- k194Notes()
    omic <- k194Omic(class)
    q <- k207Normalised(file.path(root, "fx"), class, da = TRUE)
    dirs <- k207Dirs(file.path(root, "project"), omic)
    rec <- k207Save(q$fx$workflow, omic, dirs)
    k213Result(dirs, omic, k194Owner(q$fx), root, rec)
}

#' Transcriptomics in a fresh project (MCI-032's): the run with RUV automatic (G0: k 2 chosen, so no deviation), the DA, the RUVIII export (`group`),
#' then Save.
k213TrxJourney <- function() {
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    k194TrxRun(project, omicNormReactiveState(), list(ruv_mode = "automatic"), seen)
    expect_identical(owner$current_state, "trx_ruvs", info = paste("run precondition:", paste(seen$notes, collapse = " | ")))
    kit172Quiet(runTrxReferenceDaAnalysis(owner, getWorkflowStateReference(owner), list(da_q_val_thresh = 0.05,
        treat_lfc_cutoff = 0), trxDaContrastEntries(project$workflow$contrasts_tbl)))
    kit172Quiet(runTrxReferenceRuviiiExport(owner, "group"))
    expect_false(is.null(owner$.context$snapshot$exports$trx_ruviii_export), info = "export precondition")
    rec <- k207Save(project$workflow, "transcriptomics", project$project_dirs)
    base <- project$paths$base_dir
    k213Result(project$project_dirs, "transcriptomics", owner, c(dirname(dirname(base)), base), rec)
}

#' Proteomics in a fresh root: KIT-107 F-16's fixture (MCI-009's object, here as the RUV state), the correlation filter
#' at 0.5 in process (`computeDIAReferenceCorrelation`, committed as `runDIAReferenceCorrelation` commits it), the DA, its
#' publication and the enrichment (test mode), then Save.
k213ProtJourney <- function() {
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true")
    enrichment <- omicGuiProteomicsEnrichment
    testthat::local_mocked_bindings(omicGuiProteomicsEnrichment = function() recordEnrichmentConfig(enrichment()),
        .package = "MultiScholaR")
    root <- withr::local_tempdir()
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "kit213-prot", "DIA"), WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("DIA")
    design <- module_ci_prot_da_design()
    rownames(design) <- NULL
    quant <- module_ci_prot_da_table()
    rownames(quant) <- NULL
    setDIAStepProvenance(manager, "design_builder")
    manager$saveState("ruv_corrected", module_ci_prot_da_object(design = design, protein_quant_table = quant), list(), "fixture")
    paths <- list(source_dir = file.path(root, "source"), results_dir = file.path(root, "results"),
        da_output_dir = file.path(root, "results", "da_proteins"), pathway_dir = file.path(root, "results", "pathway_enrichment"),
        publication_graphs_dir = file.path(root, "publication_graphs"))
    for (p in paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
    manager$.context$experiment_paths <- paths
    contrasts <- module_ci_prot_da_contrasts("two_group")
    workflow <- list2env(list(state_manager = manager, contrasts_tbl = contrasts, tab_status = list(),
        qc_params = list(protein_qc = list()), config_list = list()))
    withDIAStateScope(manager, function() {
        run <- computeDIAReferenceCorrelation(manager, getDIAStateReference(manager), "group", 0.5)
        commit <- diaCorrelationCommit(manager, run$final, 0.5, FALSE)
        commitDIAReferenceState(manager, commit$state_name, commit$checkpoint, commit$config, commit$description)
        da_data <- new.env(parent = emptyenv())
        da_data$current_s4_object <- getDIAStateReference(manager)
        da_data$contrasts_available <- contrasts$contrasts
        manifest <- suppressMessages(suppressWarnings(runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group",
            da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), da_data, workflow, paths, contrasts = contrasts, publish = FALSE)))
        suppressMessages(suppressWarnings(produceDIADaPublication(manager, manifest, data.frame(), 0.05,
            paths[c("da_output_dir", "publication_graphs_dir")], save_assets = FALSE)))
        utils::capture.output(suppressMessages(suppressWarnings(runWorkflowEnrichment(workflow,
            omicGuiConfig("proteomics")$enrichment, manifest, list(selected_contrast = "B_vs_A", organism_taxid = "9606",
                up_cutoff = 0, down_cutoff = 0, q_cutoff = 0.05, correction_method = "gSCS", enable_organism_filter = FALSE),
            paths))))
    })
    expect_identical(manager$current_state, "enrichment_completed", info = "enrichment precondition")
    dirs <- k207Dirs(file.path(root, "project"), "proteomics")
    rec <- k207Save(workflow, "proteomics", dirs)
    k213Result(dirs, "proteomics", manager, root, rec)
}

#' E1 to E4, G2 and G5 for one omic's two projects: the evidence lines (the differing lines after the masks, each
#' line's SHA-256 for G2, the masked digests), then the expectations.
k213Expect <- function(row, omic, a, b, expected) {
    for (p in list(a, b)) expect_true(p$written, info = paste(row, "precondition: study_parameters.txt written;", p$notes))
    ma <- k213Masked(a$lines, a$roots)
    mb <- k213Masked(b$lines, b$roots)
    n <- max(length(ma), length(mb))
    la <- c(ma, rep("<no line>", n - length(ma)))
    lb <- c(mb, rep("<no line>", n - length(mb)))
    differing <- which(la != lb)
    for (i in differing) cat(sprintf("KIT-213 %s %s differs %d | %s | %s\n", row, omic, i, la[[i]], lb[[i]]))
    # G2: each masked line that holds no L parameter, by its SHA-256 (R2 compares them base against head).
    for (i in which(!grepl(k213LPattern, ma))) cat(sprintf("KIT-213 G2 %s %s %d %s\n", row, omic, i,
        substr(k213Sha(ma[[i]]), 1L, 16L)))
    cat(sprintf("KIT-213 %s %s masked lines=%d/%d sha256=%s/%s differing=%s\n", row, omic, length(ma), length(mb), k213Sha(ma),
        k213Sha(mb), paste(differing, collapse = ",")))
    info <- paste(row, omic, "differing masked lines:", paste(sprintf("%d: %s", differing, la[differing]), collapse = " || "))
    expect_identical(k213Sha(ma), k213Sha(mb), info = info)
    for (line in expected) {
        expect_true(line %in% a$lines, info = paste(row, omic, "project 1 lacks:", line, "| L lines:",
            paste(grep(k213LPattern, a$lines, value = TRUE), collapse = " || ")))
        expect_true(line %in% b$lines, info = paste(row, omic, "project 2 lacks:", line))
    }
    # G5: no line holds the store's incarnation or one of its entity ids.
    for (p in list(a, b)) {
        tokens <- regmatches(p$lines, gregexpr("[0-9a-f]{64}", p$lines))
        bound <- p$lines[vapply(seq_along(p$lines), function(i) any(tokens[[i]] %in% p$ids$ids) ||
            grepl(p$ids$incarnation, p$lines[[i]], fixed = TRUE), logical(1))]
        cat(sprintf("KIT-213 G5 %s %s store-bound lines=%d ids=%d\n", row, omic, length(bound), length(p$ids$ids)))
        expect_identical(bound, character(), info = paste("G5", row, omic))
    }
}

# ------------------------------------------------------------------------------------------------
# E1 to E4 (with G2 and G5): two projects of one analysis.
# ------------------------------------------------------------------------------------------------

for (class in c("lipid", "metab")) {
    local({
        class <- class
        row <- c(lipid = "E1", metab = "E2")[[class]]
        test_that(sprintf("KIT-213 %s: %s, two projects of one analysis write the same study parameters; the pairs line names step 7",
            row, k194Omic(class)), {
            skip_if_not(workflowOwnerCapable())
            a <- k213FeatureJourney(class)
            b <- k213FeatureJourney(class)
            k213Expect(row, k194Omic(class), a, b,
                "  pairs [input] = the sample-pair correlations computed in step 7 of this analysis")
        })
    })
}

test_that("KIT-213 E3: transcriptomics, two projects of one analysis write the same study parameters; logcpm and grid name their step", {
    skip_if_not_installed("edgeR")
    a <- k213TrxJourney()
    b <- k213TrxJourney()
    k213Expect("E3", "transcriptomics", a, b, c(
        "  trx_ruv_grid.logcpm [input] = the TMM log counts per million computed in step 5 of this analysis",
        "  logcpm [input] = the TMM log counts per million computed in step 5 of this analysis",
        "  grid [input] = the RUV choice made in step 5 of this analysis"))
})

test_that("KIT-213 E4: proteomics, two projects of one analysis write the same study parameters; pairs and results labelled", {
    skip_if_not(workflowOwnerCapable())
    a <- k213ProtJourney()
    b <- k213ProtJourney()
    k213Expect("E4", "proteomics", a, b, c(
        "  pairs [input] = the sample-pair correlations computed in step 3 of this analysis",
        "  results [input] = groupB-groupA=the differential abundance results of contrast groupB-groupA in this analysis"))
    # A contrast without a long table (`NA`, R/utils_workflow_dia_da.R:740-741) keeps `NA`, each element named as today.
    expect_true(exists("studyParameterRecordValueText", mode = "function"), info = "E4: the label helper")
    if (exists("studyParameterRecordValueText", mode = "function")) {
        parameter <- list(name = "results", provenance = "input",
            value = encodeWorkflowRecordValue(c(B_vs_A = strrep("a", 64), C_vs_A = NA_character_)))
        expect_identical(studyParameterRecordValueText(parameter, 9L),
            "B_vs_A=the differential abundance results of contrast B_vs_A in this analysis, C_vs_A=NA")
    }
})

# ------------------------------------------------------------------------------------------------
# G3, G4: what does not hold an L parameter is written as before.
# ------------------------------------------------------------------------------------------------

test_that("KIT-213 G3: an owner without a record block writes as before (masked text's digest, compared base against head)", {
    source_dir <- withr::local_tempdir()
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    object@args$filterMinNumPeptidesPerSample <- list(inclusion_list = paste0("s", 1:40))
    suppressMessages(utils::capture.output(createWorkflowArgsFromConfig("plain", source_dir_path = source_dir,
        final_s4_object = object, workflow_data = list(state_manager = WorkflowState$new(audit_enabled = FALSE),
            config_list = list())), type = "output"))
    text <- readLines(file.path(source_dir, "study_parameters.txt"))
    masked <- k213Masked(text, source_dir)
    cat(sprintf("KIT-213 G3 lines=%d sha256=%s\n", length(masked), k213Sha(masked)))
    expect_false("Analysis Record:" %in% text)
})

test_that("KIT-213 G4: a fasta_sha256 [input] line (protein_accession_cleanup) is written as before", {
    fixture <- analysisRecordFixture(withr::local_tempdir())
    manager <- fixture$manager
    withr::defer(closeDIAStateAdapter(manager))
    source_dir <- withr::local_tempdir()
    suppressMessages(utils::capture.output(createWorkflowArgsFromConfig("record-test", source_dir_path = source_dir,
        final_s4_object = getDIAStateReference(manager, "qvalue_filtered"),
        workflow_data = list(state_manager = manager, config_list = list())), type = "output"))
    text <- readLines(file.path(source_dir, "study_parameters.txt"))
    line <- sprintf("  fasta_sha256 [input] = %s", digest::digest(file = fixture$fasta, algo = "sha256"))
    cat(sprintf("KIT-213 G4 %s\n", paste(grep("fasta_sha256", text, value = TRUE), collapse = " || ")))
    expect_true(line %in% text, info = paste(grep("fasta_sha256", text, value = TRUE), collapse = " || "))
})
