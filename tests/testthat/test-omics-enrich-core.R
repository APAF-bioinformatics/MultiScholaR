# KIT-107 W2: the enrichment core (D3, D4), written before it.
#   - each moved function is in the package namespace and its source is its new file (utils_workflow_enrichment.R,
#     utils_workflow_dia_enrichment.R);
#   - the equivalence: `runWorkflowEnrichment` on W1 (b)'s fixture (test-omic-enrich-module.R, its helpers loaded below)
#     gives the start commit's run step for the two backends and the two kernel cases: the kernel's arguments (and the
#     fake services' query and background), the product, the state digests, the history, `tab_status`, the external
#     resources, the bus values, the project's files and store counts, the display tables and the progress steps.
#     The core side compares with W1 (b)'s frozen run step and stays live. The comparison side runs the start commit's
#     `runProtEnrichAnalysisBody` (mod_prot_enrich_analysis_helpers.R) on the same fixture and compares it with the core
#     side: one block, quarantined when W3 deletes that helper (Q-166);
#   - `restoreWorkflowEnrichment` gives W1 (b)'s `restore` capture;
#   - `refuse_missing_da` refuses a run whose DA table is missing before anything is written (F107-03's guard; W3 wires
#     it into the tab);
#   - F-16 (D6): test 4, a failed commit leaves no pending provenance; test 1, the module-level record after a Run (it
#     fails until W3 puts the tab on the core: the start commit's module still commits a `legacy` step).

# W1 (b)'s helpers and frozen captures: every top-level assignment of test-omic-enrich-module.R (its blocks are not run).
for (omicsEnrichCoreW1Expression in parse(testthat::test_path("test-omic-enrich-module.R"), keep.source = FALSE)) {
    if (is.call(omicsEnrichCoreW1Expression) && identical(omicsEnrichCoreW1Expression[[1L]], as.name("<-"))) {
        eval(omicsEnrichCoreW1Expression, environment())
    }
}
rm(omicsEnrichCoreW1Expression)

#' The proteomics `enrichment` group with its function fields (D1). Until the W2 merge wires them in
#' utils_omic_gui_config.R they are optional; this fills an absent one with the value the merge wires, and changes
#' nothing once they are wired.
omicsEnrichCoreConfig <- function(config = omicGuiConfig("proteomics")$enrichment) {
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

test_that("(W2) the moved functions are in the namespace, each in its new file", {
    core <- c("resolveProtEnrichAnalysisInputColumns", "buildProtEnrichProcessEnrichmentsArgs",
        "buildProtEnrichAllContrastResults", "resolveProtEnrichProcessEnrichmentsFn", "resolveProtEnrichAnalysisMethod",
        "persistProtEnrichOrganismFilterMetadata", "resolveProtEnrichRawContrastName", "resolveProtEnrichContrastsTable",
        "resolveProtEnrichSelectedContrastResults", "resolveProtEnrichSelectedDaResults", "resolveProtEnrichRunDependencies",
        "resolveProtEnrichOutputDirectories", "buildProtEnrichSupportedOrganisms", "normaliseProtEnrichTaxonId",
        "validateProtEnrichProcessParameters", "buildProtEnrichAnalysisResultsPayload", "propagateProtEnrichResultsArgs",
        "propagateProtEnrichUiParams", "updateProtEnrichStateManagerUiParams",
        "runWorkflowEnrichment", "restoreWorkflowEnrichment", "saveWorkflowEnrichmentState", "replayWorkflowEnrichmentStep")
    dia <- c("buildProtEnrichDeterministicUniprotAnnotations", "buildProtEnrichDeterministicGprofilerResult",
        "buildProtEnrichDeterministicClusterProfilerResult", "writeProtEnrichDeterministicResultTable",
        "buildProtEnrichDeterministicPlot", "runProtEnrichDeterministicProcessEnrichments", "resolveProtEnrichOrganismMapping",
        "applyProtEnrichOrganismFilter", "resolveProtEnrichUniprotAnnotations", "resolveProtEnrichAnnotationMatching")
    moved <- c(stats::setNames(rep("utils_workflow_enrichment.R", length(core)), core),
        stats::setNames(rep("utils_workflow_dia_enrichment.R", length(dia)), dia))
    expect_length(moved, 33L)
    ns <- asNamespace("MultiScholaR")
    for (name in names(moved)) {
        expect_true(exists(name, envir = ns, inherits = FALSE), info = name)
        if (!exists(name, envir = ns, inherits = FALSE)) next
        fn <- get(name, envir = ns, inherits = FALSE)
        expect_identical(environmentName(environment(fn)), "MultiScholaR", info = name)
        # Source files are known under load_all only; the installed package keeps no srcref (KIT-104 ca35c462's rule).
        if (!is.null(utils::getSrcref(fn))) expect_identical(utils::getSrcFilename(fn), moved[[name]], info = name)
    }
    expect_true(methods::isClass("ProtEnrichDeterministicResult"))
    if (!is.null(utils::getSrcref(get("runWorkflowEnrichment", envir = ns)))) {
        source <- readLines(testthat::test_path("..", "..", "R", "utils_workflow_dia_enrichment.R"))
        expect_true(any(grepl('methods::isClass("ProtEnrichDeterministicResult")', source, fixed = TRUE)))
        # The core calls no Shiny UI function (D3): the entries notify and report progress through callbacks.
        core_source <- readLines(testthat::test_path("..", "..", "R", "utils_workflow_enrichment.R"))
        expect_false(any(grepl("shiny::(showNotification|incProgress|withProgress|removeNotification|update[A-Z])",
            core_source)))
    }
})

# ------------------------------------------------------------------------------------------------
# The equivalence (D3): W1 (b)'s fixture, the DA through the module, then the run.
# ------------------------------------------------------------------------------------------------

#' The run's inputs as the tab holds them after W1 (b)'s select step (the UI's defaults for the rest).
omicsEnrichCoreInputs <- function(taxon) list(selected_contrast = "B_vs_A", organism_taxid = taxon, up_cutoff = 0,
    down_cutoff = 0, q_cutoff = 0.05, correction_method = "gSCS", enable_organism_filter = FALSE)

#' The four cases: W1 (b)'s case, its taxon and whether the kernel is the live one.
omicsEnrichCoreCases <- list(gprofiler = list(taxon = "9606", live = FALSE),
    clusterprofiler = list(taxon = "999999", live = FALSE), kernel_gprofiler = list(taxon = "9606", live = TRUE),
    kernel_clusterprofiler = list(taxon = "999999", live = TRUE))

#' The fields of a run step that the run itself produces (the tab's outputs, notices, selector updates and returned
#' values are the frame's, W3).
omicsEnrichCoreFields <- c("kernel", "services", "states", "history", "tab_status", "resources", "bus", "files",
    "errors")

omicsEnrichCoreRuns <- new.env(parent = emptyenv())

#' One session of W1 (b)'s case up to its select step, then `run(ctx, inputs)` as the run step.
omicsEnrichCoreDrive <- function(case, run, env = parent.frame()) {
    spec <- omicsEnrichCoreCases[[case]]
    w1 <- omicEnrichW1CaseSpec(case, env)[[1L]]
    inputs <- omicsEnrichCoreInputs(spec$taxon)
    # The steps before W1 (b)'s run step (a backend case starts with `start`), so each step's "<unchanged>" fields
    # are relative to the same previous step.
    before <- intersect(c("start", "da", "open", "select"), names(omicEnrichW1Frozen[[case]]))
    steps <- c(w1$steps[before], list(run = function(ctx) run(ctx, inputs)))
    omicEnrichW1Drive(w1$setup, steps, live = spec$live)
}

#' The core side: `runWorkflowEnrichment` on the bus and the manifest the DA published, with the frame's callbacks;
#' the product and display tables kept from its returned list.
omicsEnrichCoreRun <- function(case, env = parent.frame()) {
    if (is.null(omicsEnrichCoreRuns[[case]])) {
        config <- omicsEnrichCoreConfig()
        omicsEnrichCoreRuns[[case]] <- omicsEnrichCoreDrive(case, function(ctx, inputs) shiny::isolate({
            bus <- ctx$bus()
            manifest <- bus$da_analysis_results_list
            out <- runWorkflowEnrichment(bus, config, manifest, inputs, ctx$paths,
                notify = function(message, type, duration) shiny::showNotification(message, type = type, duration = duration),
                progress = function(value, detail) shiny::incProgress(value, detail = detail))
            ctx$keep("run_product", omicEnrichW1Product(out$reference, bus$state_manager, manifest))
            ctx$keep("run_display", omicEnrichW1Display(out$contrasts, bus$state_manager))
            ctx$keep("result", c(names = paste(names(out), collapse = ","), method = out$method$method,
                raw_contrast = out$raw_contrast, selected = out$selected$found, saved = out$saved$saved,
                saved_warning = out$saved$warning %||% "",
                manifest_unchanged = identical(bus$da_analysis_results_list, manifest)))
        }), env = env)
    }
    omicsEnrichCoreRuns[[case]]
}

for (case in names(omicsEnrichCoreCases)) {
    local({
        case <- case
        test_that(sprintf("(W2, core side) runWorkflowEnrichment gives W1 (b)'s %s run step", case), {
            seen <- omicsEnrichCoreRun(case)
            frozen <- omicEnrichW1Frozen[[case]]
            expect_identical(names(seen), c(intersect(c("start", "da", "open", "select"), names(frozen)), "run"))
            for (step in setdiff(names(seen), "run")) {
                expect_identical(seen[[step]][omicsEnrichCoreFields], frozen[[step]][omicsEnrichCoreFields],
                    info = paste(case, step))
            }
            run <- seen$run
            for (field in omicsEnrichCoreFields) {
                expect_identical(run[[field]], frozen$run[[field]], info = paste(case, "run", field))
            }
            expect_identical(run$run_product, frozen$run$product, info = case)
            expect_identical(run$run_display, frozen$run$display, info = case)
            # The progress steps (the frame's withProgress start line is the frame's).
            expect_identical(run$progress, frozen$run$progress[-1L], info = case)
            expect_identical(run$notes, character(), info = case)
            expect_identical(run$result[["names"]], paste(c("results", "method", "contrasts", "selected", "raw_contrast",
                "contrasts_tbl", "pathway_dir", "annotation_match", "organism_filter", "service", "reference", "parameters",
                "saved"), collapse = ","))
            expect_identical(run$result[["method"]], if (identical(omicsEnrichCoreCases[[case]]$taxon, "9606"))
                "gprofiler2" else "clusterprofiler")
            expect_identical(run$result[["raw_contrast"]], "groupB-groupA")
            expect_identical(run$result[["selected"]], "TRUE")
            # The state commits when it did at the start commit (W1 (b): the live kernel cases' commit fails today and is
            # reported as a warning; F-16 then leaves no provenance pending).
            expect_identical(run$result[["saved"]], if (identical(frozen$run$history, "<unchanged>")) "FALSE" else "TRUE",
                info = paste(case, run$result[["saved_warning"]]))
            # F107-01: the run never replaces the manifest the caller holds.
            expect_identical(run$result[["manifest_unchanged"]], "TRUE")
        })
    })
}

# Comparison side: quarantined when W3 deletes mod_prot_enrich_analysis_helpers.R (`runProtEnrichAnalysisBody`), Q-166.
test_that("(W2, comparison side) runWorkflowEnrichment gives the start commit's runProtEnrichAnalysisBody on the same fixture", {
    for (case in names(omicsEnrichCoreCases)) {
        core <- omicsEnrichCoreRun(case)
        body <- omicsEnrichCoreDrive(case, function(ctx, inputs) shiny::isolate({
            bus <- ctx$bus()
            tab <- ctx$returned()
            runProtEnrichAnalysisBody(input = inputs, enrichmentData = tab, workflowData = bus, experimentPaths = ctx$paths,
                currentAnalysisMethodFn = function() resolveProtEnrichAnalysisMethod(inputs$organism_taxid,
                    buildProtEnrichSupportedOrganisms())$methodInfo)
            ctx$keep("run_product", omicEnrichW1Product(tab$enrichment_product, bus$state_manager, bus$da_analysis_results_list))
            ctx$keep("run_display", omicEnrichW1Display(tab$all_enrichment_results, bus$state_manager))
        }))
        for (field in c(omicsEnrichCoreFields, "run_product", "run_display", "progress", "notes")) {
            expect_identical(core$run[[field]], body$run[[field]], info = paste(case, field))
        }
    }
})

test_that("(W2) restoreWorkflowEnrichment gives W1 (b)'s restore capture", {
    config <- omicsEnrichCoreConfig()
    w1 <- omicEnrichW1CaseSpec("restore", environment())[[1L]]
    steps <- w1$steps
    steps$restore <- function(ctx) {
        manager <- shiny::isolate(ctx$bus()$state_manager)
        out <- restoreWorkflowEnrichment(manager, config)
        ctx$keep("restore", c(restored = !is.null(out),
            enrichment_results_full = omicEnrichW1Digest(out$enrichment_results_full, ctx$roots),
            omicEnrichW1Display(out$all_enrichment_results, manager)))
        ctx$keep("method", out$method$method)
        ctx$keep("absent", is.null(restoreWorkflowEnrichment(manager, utils::modifyList(config,
            list(state_name = "not_a_state")))))
    }
    seen <- omicEnrichW1Drive(w1$setup, steps)
    frozen <- omicEnrichW1Frozen$restore$restore$restore
    expect_identical(seen$restore$restore, frozen[setdiff(names(frozen), c("again", "analysis_complete"))])
    expect_identical(seen$restore$method, "gprofiler2")
    expect_true(seen$restore$absent)
    expect_null(restoreWorkflowEnrichment(NULL, config))
})

test_that("(W2) refuse_missing_da refuses a run whose DA table is missing, before anything is written (F107-03's guard)", {
    config <- omicsEnrichCoreConfig()
    setup <- omicEnrichW1Setup(env = environment())
    steps <- list(da = omicEnrichW1Da, open = omicEnrichW1Open, refuse = function(ctx) shiny::isolate({
        unlink(file.path(ctx$paths$da_output_dir, "da_proteins_groupB-groupA_long_annot.tsv"))
        bus <- ctx$bus()
        manager <- bus$state_manager
        before <- list(history = diaStateHistory(manager), products = length(manager$.context$snapshot$products))
        error <- tryCatch(runWorkflowEnrichment(bus, config, bus$da_analysis_results_list, omicsEnrichCoreInputs("9606"),
            ctx$paths, refuse_missing_da = TRUE), error = function(e) e)
        ctx$keep("refused", c(class = paste(intersect(class(error), "workflow_enrichment_input"), collapse = ""),
            message = conditionMessage(error),
            history = identical(diaStateHistory(manager), before$history),
            products = identical(length(manager$.context$snapshot$products), before$products),
            filter = is.null(bus$enrichment_organism_filter), annotations = is.null(bus$uniprot_dat_cln),
            results = is.null(bus$enrichment_analysis_results)))
    }))
    seen <- suppressWarnings(omicEnrichW1Drive(setup, steps))
    expect_identical(seen$refuse$refused, c(class = "workflow_enrichment_input",
        message = "The DA results table is missing for contrast groupB-groupA.", history = "TRUE", products = "TRUE",
        filter = "TRUE", annotations = "TRUE", results = "TRUE"))
    expect_identical(seen$refuse$history, "<unchanged>")
})

# ------------------------------------------------------------------------------------------------
# F-16 (D4, D6).
# ------------------------------------------------------------------------------------------------

test_that("(F-16, D6 test 4) the enrichment state records its provenance, and a failed commit leaves none pending", {
    root <- withr::local_tempdir()
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "kit107-provenance", "DIA"), WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    manager$saveState("correlation_filtered", makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small")),
        list(), "fixture")
    config <- omicsEnrichCoreConfig()
    results <- runProtEnrichDeterministicProcessEnrichments(da_results = module_ci_prot_enrich_da_results(), taxon_id = 9606,
        up_cutoff = 0, down_cutoff = 0, q_cutoff = 0.05, pathway_dir = withr::local_tempdir(), go_annotations = NULL,
        exclude_iea = FALSE, protein_id_column = "uniprot_acc", contrast_names = "groupB-groupA", correction_method = "gSCS")
    declared <- list(taxon_id = 9606, up_cutoff = 0, down_cutoff = 0, q_cutoff = 0.05, exclude_iea = FALSE,
        correction_method = "gSCS", protein_id_column = "uniprot_acc", contrast_names = "groupB-groupA",
        results = c(`groupB-groupA` = strrep("a", 64L)))
    method <- list(method = "gprofiler2", supported = TRUE)
    parameters <- list(selected_contrast = "B_vs_A", organism_taxid = "9606", up_cutoff = 0, down_cutoff = 0, q_cutoff = 0.05)
    workflow <- list2env(list(state_manager = manager))
    # A failed commit: the provenance set before it is cleared, the next commit is not attributed to the enrichment.
    save <- function() {
        out <- NULL
        utils::capture.output(out <- saveWorkflowEnrichmentState(workflow, config, results, declared, method, parameters,
            "pathway"))
        out
    }
    failed <- (function() {
        testthat::local_mocked_bindings(saveDIAOwnerState = function(manager, ...) stop("disk full"))
        save()
    })()
    expect_identical(failed, list(attempted = TRUE, saved = FALSE, warning = "disk full"))
    expect_null(manager$.context$pending_provenance)
    manager$saveState("next_state", makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small")), list(), "next")
    expect_identical(manager$.context$snapshot$step_provenance$next_state$operation, "legacy")
    # The commit: the registry entry's operation, variant and names, and the declared values, on the head.
    saved <- save()
    expect_identical(saved, list(attempted = TRUE, saved = TRUE, warning = NULL))
    expect_null(manager$.context$pending_provenance)
    expect_identical(manager$current_state, "enrichment_completed")
    entry <- manager$.context$snapshot$step_provenance$enrichment_completed
    expect_identical(entry$operation, "enrichment")
    expect_identical(entry$variant, "whole")
    expect_identical(entry$parameter_names, diaEnrichmentOperation()$parameters)
    expect_identical(entry$values, declared)
    expect_null(entry$fixed)
    # The state's config and description as today's save writes them (a result-class state: KIT-101 F-15's read).
    metadata <- readComponentMetadata(manager, manager$.context$snapshot$states$enrichment_completed)
    expect_identical(metadata$config, list(selected_contrast = "B_vs_A", analysis_method = "gprofiler2",
        organism_supported = TRUE, up_cutoff = 0, down_cutoff = 0, q_cutoff = 0.05, organism_taxid = "9606",
        pathway_dir = "pathway"))
    expect_identical(metadata$description, "Enrichment analysis completed using gprofiler2 for contrast: B_vs_A")
})

test_that("(F-16, D6 test 1) after a Run the record's enrichment step carries its operation, variant and 9 declared parameters (passes from W3)", {
    w1 <- omicEnrichW1CaseSpec("gprofiler", environment())[[1L]]
    steps <- c(w1$steps[c("da", "open", "select", "run")], list(record = function(ctx) {
        manager <- shiny::isolate(ctx$bus()$state_manager)
        record <- jsonlite::fromJSON(jsonlite::toJSON(buildWorkflowAnalysisRecord(manager), auto_unbox = TRUE,
            digits = I(17), null = "null", na = "null"), simplifyVector = FALSE)
        ctx$keep("record", list(step = record$steps[[length(record$steps)]],
            incomplete = validateWorkflowAnalysisRecord(record)$incomplete))
    }))
    seen <- omicEnrichW1Drive(w1$setup, steps)
    step <- seen$record$record$step
    expect_identical(step$state, "enrichment_completed")
    expect_identical(step$operation, "enrichment")
    expect_identical(step$variant, "whole")
    parameters <- step$parameters
    expect_identical(vapply(parameters, `[[`, "", "name"), diaEnrichmentOperation()$parameters)
    expect_identical(vapply(parameters, `[[`, "", "source"), rep("provenance", 9L))
    expect_identical(vapply(parameters, `[[`, "", "provenance"), c("input", "no_shipped_default", "no_shipped_default",
        "no_shipped_default", "no_shipped_default", "no_shipped_default", "input", "input", "input"))
    expect_false(any(grepl("enrichment_completed", seen$record$record$incomplete, fixed = TRUE)))
})
