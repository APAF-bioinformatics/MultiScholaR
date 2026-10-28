# KIT-030: the lipidomics and metabolomics modules and kernels hold nothing in the global environment.
# What they shared through it (config_list, contrasts_tbl, the filtering-progress trackers) travels on
# workflow_data; the figure writers take their directories as arguments; package helpers resolve by
# namespace; the DA holders keep a state name. The allowlist names the sites that remain by contract.
library(testthat)

omics_module_sources <- function() {
    root <- testthat::test_path("..", "..", "R")
    files <- list.files(root, pattern = "^(mod|func)_(lipid|metab)[A-Za-z0-9_]*\\.R$", full.names = TRUE)
    stats::setNames(lapply(files, readLines, warn = FALSE), basename(files))
}

omics_global_hits <- function(sources, pattern, allow = character()) {
    hits <- character()
    for (name in names(sources)) {
        lines <- grep(pattern, sources[[name]], perl = TRUE)
        lines <- lines[!vapply(sources[[name]][lines], function(line) any(vapply(allow, grepl, logical(1), x = line, fixed = TRUE)), logical(1))]
        if (length(lines)) hits <- c(hits, sprintf("%s:%d", name, lines))
    }
    hits
}

test_that("KIT-030: no lipid or metab module or kernel reads or writes the global environment", {
    sources <- omics_module_sources()
    expect_gt(length(sources), 100L)

    # The summaries' `globalEnv` parameter is the script/test seam the module-CI tests drive; its one
    # remaining use is the guarded `project_dirs` write the report templates read by name (the
    # lipidomics and metabolomics templates stop without a global `project_dirs`; KIT-067 owns the
    # templates, the proteomics summary does the same). Nothing else may name the global environment.
    allowed_global_env <- c("globalEnv = .GlobalEnv")
    hits <- omics_global_hits(sources, "\\.GlobalEnv|globalenv\\(\\)|as\\.environment\\(1\\)", allowed_global_env)
    expect_identical(hits, character(0))

    # No assign() or get0()/get()/exists() of the pairs that used to be shared globally.
    shared <- "assign(Fn)?\\(\"(config_list|contrasts_tbl|filtering_progress_[a-z]+|omic_type|experiment_label)\""
    expect_identical(omics_global_hits(sources, shared), character(0))
    lookups <- "(get0|get|exists)\\(\"(config_list|contrasts_tbl|filtering_progress_[a-z]+|project_dirs|omic_type|experiment_label)\""
    expect_identical(omics_global_hits(sources, lookups), character(0))
})

test_that("KIT-030: the super-assignments in the lipid and metab code are the listed local closures", {
    sources <- omics_module_sources()
    # Each is a `<<-` to a variable of the enclosing function (an accumulator inside a tryCatch handler,
    # purrr callback or local helper), never to a global: listed so a new one is a review item.
    allowed_files <- c(
        "mod_lipid_design_validation.R", "mod_metab_design_validation.R",
        "mod_lipid_qc_duplicates_server_helpers.R", "mod_metab_qc_duplicates_server_helpers.R",
        "func_lipid_da_model.R", "func_metab_da_core_engine.R",
        "func_lipid_da_results.R", "func_metab_da_orchestration.R",
        "func_lipid_s4_correlation_methods.R", "func_metab_s4_qc_methods.R"
    )
    hits <- omics_global_hits(sources, "<<-")
    files_hit <- unique(sub(":.*$", "", hits))
    expect_identical(setdiff(files_hit, allowed_files), character(0))
})

test_that("KIT-030: the filtering-progress trackers live on the omic's bus inside a session and in the global environment outside one", {
    withr::local_options(list(shiny.testmode = TRUE))
    tracked <- c("filtering_progress_lipidomics", "filtering_progress_metabolomics")
    had <- vapply(tracked, exists, logical(1), envir = globalenv(), inherits = FALSE)
    saved <- lapply(tracked[had], get, envir = globalenv(), inherits = FALSE)
    withr::defer({
        for (name in tracked) if (exists(name, envir = globalenv(), inherits = FALSE)) rm(list = name, envir = globalenv())
        for (name in names(saved)) assign(name, saved[[name]], envir = globalenv())
    })
    for (name in tracked) if (exists(name, envir = globalenv(), inherits = FALSE)) rm(list = name, envir = globalenv())

    # Outside a session: the documented script behaviour (a global object, initialised once).
    expect_identical(workflowOmicsProgressOwner("lipidomics"), globalenv())
    first <- suppressMessages(getFilteringProgressLipidomics())
    expect_s4_class(first, "FilteringProgressLipidomics")
    expect_true(exists("filtering_progress_lipidomics", envir = globalenv(), inherits = FALSE))

    # Inside a session with a registered bus: the bus owns it and the global environment is untouched.
    bus <- shiny::reactiveValues(filtering_progress_metabolomics = NULL)
    session <- shiny::MockShinySession$new()
    shiny::withReactiveDomain(session, {
        registerWorkflowOmicsBus(session, "metabolomics", bus)
        expect_identical(workflowOmicsProgressOwner("metabolomics"), bus)
        progress <- suppressMessages(getFilteringProgressMetabolomics())
        expect_s4_class(progress, "FilteringProgressMetabolomics")
        updated <- updateFilteringProgressMetabolomics(
            progress, "1_Raw_Data", "LCMS_Pos",
            list(LCMS_Pos = list(n_metabolites = 3L, detected_per_sample = NULL, missingness = NULL,
                                 sum_intensity_per_sample = NULL, cv_distribution = NULL, is_metrics = NULL)),
            total_metabolites = 3L
        )
        expect_identical(shiny::isolate(bus$filtering_progress_metabolomics)@steps, "1_Raw_Data")
        expect_identical(shiny::isolate(bus$filtering_progress_metabolomics)@n_metabolites_total, 3)
        # The figure writers never read a workbook's globals inside a session.
        expect_null(workflowScriptGlobal("project_dirs"))
    })
    expect_false(exists("filtering_progress_metabolomics", envir = globalenv(), inherits = FALSE))
})
