# KIT-116 W2, quarantine row Q-64: the frozen values of seven blocks of test-omics-qc-core.R, computed from the lipid
# QC functions KIT-104 W5 deleted (and the metab functions those blocks compared with), at their green commit 7b99e7e8
# (KIT-104 W2), where the comparison passed. Run from a checkout of 7b99e7e8:
#
#   Rscript q64-qc.R <checkout of 7b99e7e8> <output dir>
#
# The script loads the package and test helpers of that checkout and the test file's own helper definitions (qcCoreSpy,
# qcCoreNotes, qcCorePlotParts, qcCoreLipidItsd, ...), runs each block's old side on the block's inputs, and checks
# each old value against the core of the same commit as the block did. Values are literals, or sha256 of the value the
# block compared with expect_identical. Output: <output dir>/q64-frozen.R and q64-old-values.rds.

args <- commandArgs(trailingOnly = TRUE)
green <- normalizePath(args[[1]])
out_dir <- args[[2]]
stopifnot(identical(system2("git", c("-C", green, "rev-parse", "--short=8", "HEAD"), stdout = TRUE), "7b99e7e8"))
pkgload::load_all(green, quiet = TRUE)
withr::local_dir(green)
# The helpers defer their temporary owners to testthat's teardown environment, which test_file initialises.
testthat:::local_teardown_env()
testthat::source_test_helpers(file.path(green, "tests", "testthat"), env = globalenv())
library(testthat)

exprs <- parse(file.path(green, "tests", "testthat", "test-omics-qc-core.R"), keep.source = FALSE)
for (e in exprs) if (!(is.call(e) && identical(e[[1L]], as.name("test_that")))) eval(e, globalenv())

digestOf <- function(x) digest::digest(x, algo = "sha256")
frozen <- list()
old_values <- list()
unfaithful <- character()
record <- function(group, key, value, core_value) {
    if (is.null(frozen[[group]])) frozen[[group]] <<- list()
    frozen[[group]][key] <<- list(value)
    if (!identical(value, core_value)) unfaithful <<- c(unfaithful, paste(group, key))
}
frozenCode <- function(frozen, prefix) {
    vapply(names(frozen), function(group) {
        values <- frozen[[group]]
        keys <- names(values)
        keys <- ifelse(make.names(keys) == keys, keys, paste0("`", keys, "`"))
        lines <- vapply(seq_along(values), function(i) {
            paste0("    ", keys[[i]], " = ", paste(deparse(values[[i]], width.cutoff = 500L), collapse = " "))
        }, character(1))
        code <- paste0(prefix, group, " <- list(\n", paste(lines, collapse = ",\n"), ")")
        stopifnot(identical(eval(parse(text = code)), values))
        code
    }, character(1))
}

# The counts table: lipid's inline table (registerLipidIntensityApplyFilterObserver) and metab's
# buildMetabIntensityFilterStats, on the counts the observer's reference step returned.
(function() {
    workflow <- module_ci_lipid_qc_workflow()
    captured <- new.env(parent = emptyenv())
    qcCoreSpy("runLipidReferenceIntensityFilter", captured)
    stats <- shiny::reactiveVal(NULL)
    suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(function(input, output, session) {
        registerLipidIntensityApplyFilterObserver(input, workflow, output, stats, shiny::reactiveVal(NULL), "lipidomics",
            updateLipidFilteringFn = function(...) NULL)
    }, {
        session$setInputs(intensity_cutoff_percentile = 40L, proportion_below_cutoff = 0.5, apply_filter = 1L)
        captured$stats <- stats()
    }))))
    counts <- captured$runLipidReferenceIntensityFilter$value$counts
    stopifnot(is.data.frame(captured$stats))
    record("counts", "lipid", captured$stats, omicsQcCountsTable(counts))
    record("counts", "metab", buildMetabIntensityFilterStats(originalCounts = counts$before,
        filteredCounts = counts$after), omicsQcCountsTable(counts))
    odd <- list(before = c(A = 3L, B = 7L), after = c(A = 1L, B = 7L))
    record("counts", "odd_metab", buildMetabIntensityFilterStats(originalCounts = odd$before,
        filteredCounts = odd$after), omicsQcCountsTable(odd))
})()

# Duplicate detection: handleLipidDuplicateDetection and detectMetabDuplicateFeatures.
(function() {
    quiet <- function(...) invisible(NULL)
    lipid <- module_ci_lipid_qc_workflow(module_ci_lipid_qc_object(layout = "combined"))
    today <- handleLipidDuplicateDetection(lipid, logInfoFn = quiet)
    core <- detectWorkflowOmicsDuplicates(lipidomicsDeclarationBinding, lipid$state_manager)
    record("duplicates", "lipid_list", digestOf(today$duplicatesList), digestOf(core$duplicates))
    record("duplicates", "lipid_total", today$totalDuplicates, today$totalDuplicates)
    if (!isTRUE(all.equal(core$total, today$totalDuplicates))) unfaithful <<- c(unfaithful, "lipid_total")
    metab <- module_ci_metab_qc_workflow(module_ci_metab_qc_object(layout = "combined"))
    today <- detectMetabDuplicateFeatures(metab$state_manager)
    core <- detectWorkflowOmicsDuplicates(metabolomicsDeclarationBinding, metab$state_manager)
    record("duplicates", "metab_list", digestOf(today$duplicatesList), digestOf(core$duplicates))
    record("duplicates", "metab_total", today$totalDuplicates, today$totalDuplicates)
    if (!isTRUE(all.equal(core$total, today$totalDuplicates))) unfaithful <<- c(unfaithful, "metab_total")
    none <- module_ci_lipid_qc_workflow(module_ci_lipid_qc_object(include_duplicates = FALSE))
    none_total <- handleLipidDuplicateDetection(none, logInfoFn = quiet)$totalDuplicates
    record("duplicates", "none_total", none_total, none_total)
    if (!isTRUE(all.equal(detectWorkflowOmicsDuplicates(lipidomicsDeclarationBinding, none$state_manager)$total,
        none_total))) unfaithful <<- c(unfaithful, "none_total")
    record("duplicates", "other_omic", qcCoreError(handleLipidDuplicateDetection(metab, logInfoFn = quiet)),
        qcCoreError(detectWorkflowOmicsDuplicates(lipidomicsDeclarationBinding, metab$state_manager)))
})()

# Resolution statistics: handleLipidDuplicateResolution's statsList and runMetabDuplicateResolutionWorkflow's shown
# statistics.
(function() {
    lipid <- module_ci_lipid_qc_workflow(module_ci_lipid_qc_object(layout = "combined"))
    captured <- new.env(parent = emptyenv())
    qcCoreSpy("runLipidReferenceDuplicateResolution", captured)
    names_before <- readLipidWorkflowMetadata(lipid$state_manager)$assay_names
    today <- handleLipidDuplicateResolution(lipid, "lipidomics", updateFilteringFn = function(...) NULL,
        logWarnFn = function(...) invisible(NULL))
    value <- captured$runLipidReferenceDuplicateResolution$value
    old_values$resolution_lipid <<- list(today = today, value = value)
    record("resolution", "lipid", today$statsList, omicsDuplicateResolutionStats(value, names_before))
    metab <- module_ci_metab_qc_workflow(module_ci_metab_qc_object(layout = "combined"))
    qcCoreSpy("runMetabReferenceDuplicateResolution", captured)
    shown <- new.env(parent = emptyenv())
    runMetabDuplicateResolutionWorkflow(metab, "metabolomics", setResolutionStatsFn = function(stats) shown$stats <- stats,
        setFilterPlotFn = function(plot) invisible(NULL),
        applyResolutionStateFn = function(currentS4, statsList, setResolutionStatsFn, ...) {
            setResolutionStatsFn(statsList)
            list(stateName = "metab_duplicates_resolved")
        })
    value <- captured$runMetabReferenceDuplicateResolution$value
    record("resolution", "metab", shown$stats,
        omicsDuplicateResolutionStats(value, readMetabWorkflowMetadata(metab$state_manager)$assay_names))
})()

# Lipid's inline ITSD analysis (mod_lipid_qc_itsd_server before W5): values, text, the count notice, refusals.
(function() {
    itsd <- omicGuiConfig("lipidomics")$qc$itsd
    for (pattern in c("", "^IS_By_Name$", "^NO_MATCH$")) {
        workflow <- module_ci_lipid_qc_workflow(module_ci_lipid_qc_object(layout = "combined"))
        today <- qcCoreLipidItsd(workflow, pattern)
        core <- analyzeWorkflowOmicsItsd(lipidomicsDeclarationBinding, workflow$state_manager, pattern, itsd)
        old_values$itsd_lipid[[paste0("pattern", pattern)]] <<- today
        record("itsd_lipid", paste0("metrics", pattern), digestOf(today$metrics), digestOf(core$metrics))
        record("itsd_lipid", paste0("long", pattern), digestOf(today$long), digestOf(core$long_data))
        record("itsd_lipid", paste0("text", pattern), today$text, core$result_text)
        record("itsd_lipid", paste0("note", pattern), today$notes[[2L]],
            paste0("[message] ", sprintf("Found %d internal standards", core$total)))
    }
    cases <- list(no_standard = list(object = module_ci_lipid_qc_object(layout = "gc", include_itsd = FALSE),
        pattern = "^IS_"), no_pattern = list(object = module_ci_lipid_qc_object(layout = "gc",
        internal_standard_regex = NA_character_), pattern = ""))
    for (name in names(cases)) {
        case <- cases[[name]]
        workflow <- module_ci_lipid_qc_workflow(case$object)
        today <- qcCoreLipidItsd(workflow, case$pattern)
        message <- qcCoreError(analyzeWorkflowOmicsItsd(lipidomicsDeclarationBinding, workflow$state_manager,
            case$pattern, itsd))
        record("itsd_lipid", paste0("refusal_", name), today$notes[[2L]],
            paste("[error] Error analyzing internal standards:", message))
    }
    # F104: lipid's inline analysis on a NULL pattern (the test's literal, checked here at the green commit).
    workflow <- module_ci_lipid_qc_workflow(module_ci_lipid_qc_object(layout = "combined"))
    note <- qcCoreLipidItsd(workflow, NULL)$notes[[2L]]
    record("itsd_lipid", "null_pattern", note, "[error] Error analyzing internal standards: argument is of length zero")
})()

# The ITSD figures: metab's buildMetabQcItsdCvPlot and buildMetabQcItsdIntensityPlot on analyzeMetabQcItsdData's
# analysis, and the data lipid's registerLipidItsdCvPlotOutput drew.
(function() {
    workflow <- module_ci_metab_qc_workflow(module_ci_metab_qc_object(layout = "combined"))
    analysis <- analyzeMetabQcItsdData(readMetabQcItsdState(workflow$state_manager, ""), inputPattern = "",
        logInfoFn = function(...) invisible(NULL))
    core <- analyzeWorkflowOmicsItsd(metabolomicsDeclarationBinding, workflow$state_manager, "",
        omicGuiConfig("metabolomics")$qc$itsd)
    record("figures", "metab_metrics", digestOf(analysis$metrics), digestOf(core$metrics))
    record("figures", "metab_long", digestOf(analysis$longData), digestOf(core$long_data))
    record("figures", "cv_plot", digestOf(qcCorePlotParts(buildMetabQcItsdCvPlot(analysis$metrics))),
        digestOf(qcCorePlotParts(omicsItsdCvPlot(analysis$metrics))))
    record("figures", "intensity_plot", digestOf(qcCorePlotParts(buildMetabQcItsdIntensityPlot(analysis$longData))),
        digestOf(qcCorePlotParts(omicsItsdIntensityPlot(analysis$longData))))
    drawn <- new.env(parent = emptyenv())
    ggplot <- ggplot2::ggplot
    testthat::local_mocked_bindings(ggplot = function(data = NULL, ...) {
        drawn$data <- c(drawn$data, list(data))
        ggplot(data, ...)
    }, .package = "ggplot2")
    lipid <- module_ci_lipid_qc_workflow(module_ci_lipid_qc_object(layout = "combined"))
    metrics <- analyzeWorkflowOmicsItsd(lipidomicsDeclarationBinding, lipid$state_manager, "",
        omicGuiConfig("lipidomics")$qc$itsd)$metrics
    suppressWarnings(utils::capture.output(shiny::testServer(function(input, output, session) {
        registerLipidItsdCvPlotOutput(output, shiny::reactiveVal(metrics))
    }, {
        output$cv_plot
    })))
    stopifnot(length(drawn$data) == 1L)
    lipid_drawn <- drawn$data[[1L]]
    record("figures", "lipid_cv_data", digestOf(lipid_drawn), digestOf(omicsItsdCvPlot(metrics)$data))
})()

# Lipid's finalise validation (validateLipidQcS4FinalizeState): the pass and five refusals.
(function() {
    qc <- qcCoreConfig("lipidomics")
    aligned <- readLipidWorkflowShell(module_ci_lipid_qc_workflow()$state_manager, measures = "features")
    # The pass returned its input unchanged. The shell carries its owner's per-run entity reference
    # (dia_component_source), so the frozen value is that identity, not a digest of the shell.
    passed <- validateLipidQcS4FinalizeState(aligned)
    record("finalise_lipid", "aligned_returned_unchanged", identical(passed, aligned),
        identical(validateWorkflowOmicsQcFinalize(aligned, qc, lipidomicsDeclarationBinding), passed))
    misaligned <- readLipidWorkflowShell(module_ci_lipid_qc_workflow(qcCoreBlankGroup(module_ci_lipid_qc_object()))$state_manager)
    zero_rows <- aligned
    attr(zero_rows, "assay_rows")[[1L]] <- 0L
    no_design <- aligned
    no_design@design_matrix <- no_design@design_matrix[0L, , drop = FALSE]
    no_assays <- aligned
    no_assays@lipid_data <- list()
    attr(no_assays, "assay_rows") <- integer()
    other_class <- readMetabWorkflowShell(module_ci_metab_qc_workflow()$state_manager)
    cases <- list(misaligned = misaligned, zero_rows = zero_rows, no_design = no_design, no_assays = no_assays,
        other_class = other_class)
    for (name in names(cases)) {
        today <- qcCoreError(validateLipidQcS4FinalizeState(cases[[name]]))
        stopifnot(!is.null(today))
        record("finalise_lipid", name, today,
            qcCoreError(validateWorkflowOmicsQcFinalize(cases[[name]], qc, lipidomicsDeclarationBinding)))
    }
})()

if (length(unfaithful)) stop("the old value is not the green core's for: ", paste(unfaithful, collapse = ", "))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(old_values, file.path(out_dir, "q64-old-values.rds"))
writeLines(frozenCode(frozen, "qcCoreFrozen_"), file.path(out_dir, "q64-frozen.R"))
cat("frozen groups:", paste(names(frozen), lengths(frozen), sep = "=", collapse = " "), "\n")
