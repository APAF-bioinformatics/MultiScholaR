# GUI audit R4, P-19 / A105-01: a manual internal-standard (ITSD) selection is the rows the user ticked in the table
# the tab showed. Both paths resolve those rows from the shown table's feature ids (`normData$itsd_candidate_ids`).
# Before the repair the in-process path rebuilt the table from the head at Run time and applied the same row numbers,
# so a second Run without Reset (the head is then the normalised state, without the standards) used other features.

itsdSelectionSpec <- function(omic) {
    switch(omic,
        lipid = list(omic_type = "lipidomics", binding = lipidomicsDeclarationBinding, object = module_ci_lipid_norm_object,
            workflow = module_ci_lipid_norm_workflow, data = module_ci_lipid_norm_data, qc_key = "lipid_qc_dir",
            qc_kernel = "generateLipidQcPlots"),
        metab = list(omic_type = "metabolomics", binding = metabolomicsDeclarationBinding, object = module_ci_metab_norm_object,
            workflow = module_ci_metab_norm_workflow, data = module_ci_metab_norm_data, qc_key = "metabolite_qc_dir",
            qc_kernel = "generateMetabQcPlots"))
}

#' Open the tab's ITSD tables on the QC-complete state, tick the first row of each assay's table, then Run in process
#' `runs` times without a Reset. Returns the ticked rows' feature ids, the ids each Run handed to the ITSD step, each
#' Run's value, the log and the notifications.
itsdSelectionCase <- function(omic, runs = 2L, before_run = NULL) {
    spec <- itsdSelectionSpec(omic)
    object <- spec$object(positive_only = TRUE)
    root <- withr::local_tempdir(pattern = paste0("r4-itsd-", omic, "-"))
    paths <- stats::setNames(list(file.path(root, "qc"), file.path(root, "source"), file.path(root, "export")),
        c(spec$qc_key, "source_dir", "export_dir"))
    workflow <- spec$workflow(object)
    norm_data <- spec$data(object)
    norm_data$post_filter_obj <- NULL
    config <- omicGuiConfig(spec$omic_type)
    rec <- new.env(parent = emptyenv())
    rec$resolved <- list()
    rec$log <- character()
    rec$notes <- character()
    itsd_step <- runWorkflowOmicsReferenceItsdNormalisation
    stubs <- list(function(binding, manager, aggregation, feature_ids, config, ...) {
        rec$resolved[[length(rec$resolved) + 1L]] <- feature_ids
        itsd_step(binding, manager, aggregation, feature_ids, config, ...)
    }, function(...) invisible(NULL))
    names(stubs) <- c("runWorkflowOmicsReferenceItsdNormalisation", spec$qc_kernel)
    do.call(testthat::local_mocked_bindings, c(stubs, list(.package = "MultiScholaR")))
    # The tables the tab shows, built once from the head when the tab opens (omicNormItsdTables).
    omicNormItsdTables(new.env(parent = emptyenv()), workflow, norm_data, config)
    shown <- norm_data$itsd_candidate_ids
    norm_data$itsd_selections <- lapply(shown, function(ids) 1L)
    ticked <- lapply(shown, function(ids) ids[1L])
    settings <- list(apply_itsd = TRUE, itsd_aggregation = "median", itsd_feature_ids = NULL, log_offset = 1,
        norm_method = "none", ruv_mode = "skip", ruv_parameters = list(ruv_mode = "skip"), grouping_variable = "group",
        shape_variable = "batch", assay_names = norm_data$assay_names, qc_dir = paths[[1L]], experiment_paths = paths)
    values <- logical()
    for (run in seq_len(runs)) {
        if (!is.null(before_run)) before_run(run, workflow, spec)
        values[[run]] <- isTRUE(suppressWarnings(suppressMessages(runWorkflowOmicsInProcessNormalisation(spec$binding,
            workflow$state_manager, settings, norm_data, workflow$config_list, function(entry) rec$log <- c(rec$log, entry),
            advance = function(detail) invisible(NULL),
            notify = function(message, type = "default", ...) rec$notes <- c(rec$notes, paste0("[", type, "] ", message)),
            finish = function() invisible(NULL), itsd_table = config$norm$itsd_table))))
    }
    list(shown = shown, ticked = ticked, resolved = rec$resolved, values = values, log = rec$log, notes = rec$notes,
        history = as.character(workflow$state_manager$getHistory()))
}

for (omic in c("lipid", "metab")) local({
    omic <- omic
    test_that(sprintf("P-19: the first in-process Run uses the ticked rows' features (%s)", omic), {
        case <- itsdSelectionCase(omic, runs = 1L)
        expect_true(all(lengths(case$shown) > 1L))
        expect_true(case$values[[1L]])
        expect_identical(case$resolved[[1L]], case$ticked)
    })

    test_that(sprintf("P-19: a second in-process Run without Reset still resolves the ticked rows' features (%s)", omic), {
        case <- itsdSelectionCase(omic, runs = 2L)
        expect_length(case$resolved, 2L)
        expect_identical(case$resolved[[2L]], case$ticked)
    })
})
