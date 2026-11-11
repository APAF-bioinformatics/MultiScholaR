# KIT-104 W2: the QC core for lipid, metab and trx, written before it exists. Each shared function is compared with
# today's lipid, metab and trx function on the MCI fixtures (helper-module-ci-{lipid,metab}-qc.R) and the trx_screen
# corpus: the counts table, the duplicate detection and resolution statistics, the ITSD analysis and its figures, the
# finalise validation and its alignment, and the QC action's two branches (D6). Back-reference #9
# (`resolveMetabDuplicateAssayData`) and trx's QC figure functions are checked where they moved.

#' The ggplot of a figure with its built layer data, for comparing two builders of the same figure.
qcCorePlotParts <- function(plot) {
    list(data = plot$data, labels = plot$labels, layers = ggplot2::ggplot_build(plot)$data)
}

#' The omic's qc configuration with the finalise's alignment (the lead wires it at the phase B merge).
qcCoreConfig <- function(omic_type) {
    qc <- omicGuiConfig(omic_type)$qc
    qc$finalize_alignment <- switch(omic_type,
        transcriptomics = trxDesignAlignment,
        function(shell) omicsQcFinalizeAlignment(shell, omic_type))
    qc
}

#' An object whose first design run has a blank group: valid for the class, refused by the alignment.
qcCoreBlankGroup <- function(object) {
    object@design_matrix[[object@group_id]][[1L]] <- NA_character_
    object
}

#' The message of an error, or NULL when `code` succeeds.
qcCoreError <- function(code) tryCatch({
    force(code)
    NULL
}, error = function(error) conditionMessage(error))

# ------------------------------------------------------------------------------------------------
# The intensity filter's counts table.
# ------------------------------------------------------------------------------------------------

# KIT-116 W2 (Q-64, Q-75): blocks that compared the core with functions KIT-104 W5 and W6 deleted now compare it with
# those functions' values, frozen where the comparison passed (7b99e7e8 for Q-64, 827f2328 for Q-75) by
# docs/qa/workflow-execution/kit/kit-116/freeze/q64-qc.R and q75-qc.R, which also checked each old value against the
# core of that commit. Small values are literals; others are the sha256 of the value the block compared with expect_identical. The
# core side stays live.
qcCoreDigest <- function(x) digest::digest(x, algo = "sha256")

# The table lipid's apply-filter observer (registerLipidIntensityApplyFilterObserver) showed for its reference step's
# counts (40th percentile, proportion 0.5), and metab's buildMetabIntensityFilterStats on the same counts and on `odd`,
# frozen from 7b99e7e8 (KIT-116).
qcCoreFrozenCounts <- list(
    lipid = structure(list(Assay = c("LCMS_Pos", "LCMS_Neg"), Original = c(11, 11), Filtered = c(8, 8),
        Removed = c(3, 3), Percent_Retained = c(72.7, 72.7)), row.names = c(NA, -2L), class = "data.frame"),
    metab = structure(list(Assay = c("LCMS_Pos", "LCMS_Neg"), Original = c(11, 11), Filtered = c(8, 8),
        Removed = c(3, 3), Percent_Retained = c(72.7, 72.7)), row.names = c(NA, -2L), class = "data.frame"),
    odd_metab = structure(list(Assay = c("A", "B"), Original = c(3, 7), Filtered = c(1, 7), Removed = c(2, 0),
        Percent_Retained = c(33.3, 100)), row.names = c(NA, -2L), class = "data.frame"))

test_that("omicsQcCountsTable is today's lipid (inline) and metab (buildMetabIntensityFilterStats) table (frozen from 7b99e7e8, KIT-116)", {
    workflow <- module_ci_lipid_qc_workflow()
    # The reference step the observer ran, with its inputs.
    counts <- runLipidReferenceIntensityFilter(workflow$state_manager, 40L, 0.5, workflow$config_list)$value$counts
    expect_identical(omicsQcCountsTable(counts), qcCoreFrozenCounts$lipid)
    expect_identical(omicsQcCountsTable(counts), qcCoreFrozenCounts$metab)
    # Rounding to one decimal, as both tabs round.
    odd <- list(before = c(A = 3L, B = 7L), after = c(A = 1L, B = 7L))
    expect_identical(omicsQcCountsTable(odd)$Percent_Retained, c(33.3, 100))
    expect_identical(omicsQcCountsTable(odd), qcCoreFrozenCounts$odd_metab)
})

# ------------------------------------------------------------------------------------------------
# Duplicates.
# ------------------------------------------------------------------------------------------------

# handleLipidDuplicateDetection and detectMetabDuplicateFeatures: the sha256 of the duplicates lists, the totals, and
# lipid's refusal of a metab state, frozen from 7b99e7e8 (KIT-116).
qcCoreFrozenDuplicates <- list(
    lipid_list = "97f418a7cf2cf55b4f66ee26c7b00f9b875d49d274f19462e8b4a2b174c4f2b1",
    lipid_total = 4,
    metab_list = "541bc986115f4812f522ee0de42666b1685abc4893e3874fd3509f9d558341a2",
    metab_total = 4L,
    none_total = 0,
    other_omic = "no slot of name \"lipid_id_column\" for this object of class \"MetaboliteAssayData\"")

test_that("detectWorkflowOmicsDuplicates finds today's duplicates (lipid and metab detection) (frozen from 7b99e7e8, KIT-116)", {
    lipid <- module_ci_lipid_qc_workflow(module_ci_lipid_qc_object(layout = "combined"))
    core <- detectWorkflowOmicsDuplicates(lipidomicsDeclarationBinding, lipid$state_manager)
    expect_identical(qcCoreDigest(core$duplicates), qcCoreFrozenDuplicates$lipid_list)
    expect_equal(core$total, qcCoreFrozenDuplicates$lipid_total)
    expect_gt(core$total, 0)

    metab <- module_ci_metab_qc_workflow(module_ci_metab_qc_object(layout = "combined"))
    core <- detectWorkflowOmicsDuplicates(metabolomicsDeclarationBinding, metab$state_manager)
    expect_identical(qcCoreDigest(core$duplicates), qcCoreFrozenDuplicates$metab_list)
    expect_equal(core$total, qcCoreFrozenDuplicates$metab_total)

    none <- module_ci_lipid_qc_workflow(module_ci_lipid_qc_object(include_duplicates = FALSE))
    expect_equal(detectWorkflowOmicsDuplicates(lipidomicsDeclarationBinding, none$state_manager)$total,
        qcCoreFrozenDuplicates$none_total)
    # Another omic's state is refused as today (the id reader fails on the binding's id slot).
    expect_identical(qcCoreError(detectWorkflowOmicsDuplicates(lipidomicsDeclarationBinding, metab$state_manager)),
        qcCoreFrozenDuplicates$other_omic)
})

# The statistics handleLipidDuplicateResolution returned and runMetabDuplicateResolutionWorkflow showed, frozen from
# 7b99e7e8 (KIT-116).
qcCoreFrozenResolution <- list(
    lipid = list(LCMS_Pos = list(original = 13L, resolved = 11L, removed = 2L),
        GCMS = list(original = 13L, resolved = 11L, removed = 2L)),
    metab = list(LCMS_Pos = list(original = 13L, resolved = 11L, removed = 2L),
        GCMS = list(original = 13L, resolved = 11L, removed = 2L)))

test_that("omicsDuplicateResolutionStats reads both products as today's lipid and metab resolutions do (frozen from 7b99e7e8, KIT-116)", {
    # The reference steps those callers ran, with their inputs.
    lipid <- module_ci_lipid_qc_workflow(module_ci_lipid_qc_object(layout = "combined"))
    names_before <- readLipidWorkflowMetadata(lipid$state_manager)$assay_names
    value <- runLipidReferenceDuplicateResolution(lipid$state_manager, lipid$config_list)$value
    expect_true("counts" %in% names(value) && !"stats" %in% names(value))
    expect_identical(omicsDuplicateResolutionStats(value, names_before), qcCoreFrozenResolution$lipid)

    metab <- module_ci_metab_qc_workflow(module_ci_metab_qc_object(layout = "combined"))
    value <- runMetabReferenceDuplicateResolution(metab$state_manager, metab$config_list)$value
    expect_true("stats" %in% names(value))
    expect_identical(omicsDuplicateResolutionStats(value, readMetabWorkflowMetadata(metab$state_manager)$assay_names),
        qcCoreFrozenResolution$metab)
    # An empty statistics list stays empty (metab's "No assay statistics available").
    expect_identical(omicsDuplicateResolutionStats(list(stats = list()), character()), list())
})

test_that("#9: metabolomicsDeclarationBinding$resolve_duplicates works through the moved resolveMetabDuplicateAssayData", {
    # Source references are kept under load_all only; the check layout's installed package has none.
    if (!is.null(utils::getSrcref(resolveMetabDuplicateAssayData))) {
        expect_identical(basename(utils::getSrcFilename(resolveMetabDuplicateAssayData, full.names = TRUE)),
            "func_metab_s4_qc_methods.R")
    }
    object <- module_ci_metab_qc_object(layout = "combined")
    resolved <- metabolomicsDeclarationBinding$resolve_duplicates(object)
    direct <- resolveMetabDuplicateAssayData(assayList = object@metabolite_data,
        metaboliteIdCol = object@metabolite_id_column, logWarnFn = function(...) NULL)
    expect_identical(resolved$value, list(stats = direct$statsList))
    expect_identical(resolved$object@metabolite_data, direct$resolvedAssayList)
    expect_true(all(vapply(direct$statsList, `[[`, numeric(1), "removed") >= 2))
    expect_false(any(vapply(resolved$object@metabolite_data, function(assay) anyDuplicated(assay$database_identifier) > 0L,
        logical(1))))
})

# ------------------------------------------------------------------------------------------------
# ITSD.
# ------------------------------------------------------------------------------------------------

# Lipid's inline ITSD analysis (mod_lipid_qc_itsd_server before W5) per pattern: the sha256 of its metrics and long
# data, its result text and its second notification; and its notifications for the two refusals. Frozen from 7b99e7e8
# (KIT-116).
qcCoreFrozenLipidItsdText <- function(pattern) {
    paste0("Internal Standard Analysis Complete\n====================================\nPattern used: ", pattern,
        "\nTotal IS detected: 4\n\nCV Statistics:\n  Median CV: 0.9%\n  Max CV: 0.9%\n  IS with CV > 30%: 0\n\n",
        "Per-Assay IS Counts:\n  LCMS_Pos: 2 internal standards\n  GCMS: 2 internal standards")
}
qcCoreFrozenLipidItsd <- list(
    metrics = "27cc90f8c1a941238c8c0c85e4cfa45369c8df0bb293e3616ddbb1a7d91916a8",
    long = "6b65c82ecf0f39a92294f505785ec0f9ee2b802f89c9426336830e03d1e593d0",
    # Lipid's inline analysis showed "^IS_" when the pattern was empty (its default).
    text = c(qcCoreFrozenLipidItsdText("^IS_"), qcCoreFrozenLipidItsdText("^IS_By_Name$"),
        qcCoreFrozenLipidItsdText("^NO_MATCH$")),
    note = "[message] Found 4 internal standards",
    refusal_no_standard = paste("[error] Error analyzing internal standards: No internal standards found. Pattern: ^IS_",
        "Columns searched: lipid_id, lipid", "Tip: Check if your data has a column with lipid names (e.g., 'Lipid name')",
        sep = "\n"),
    refusal_no_pattern = paste("[error] Error analyzing internal standards: No internal standard pattern provided.",
        "Please enter a regex pattern to identify IS."))

test_that("analyzeWorkflowOmicsItsd gives lipid's inline analysis: values, text and refusals (frozen from 7b99e7e8, KIT-116)", {
    itsd <- omicGuiConfig("lipidomics")$qc$itsd
    patterns <- c("", "^IS_By_Name$", "^NO_MATCH$")
    for (i in seq_along(patterns)) {
        pattern <- patterns[[i]]
        workflow <- module_ci_lipid_qc_workflow(module_ci_lipid_qc_object(layout = "combined"))
        core <- analyzeWorkflowOmicsItsd(lipidomicsDeclarationBinding, workflow$state_manager, pattern, itsd)
        expect_identical(qcCoreDigest(core$metrics), qcCoreFrozenLipidItsd$metrics, info = pattern)
        expect_identical(qcCoreDigest(core$long_data), qcCoreFrozenLipidItsd$long, info = pattern)
        expect_identical(core$result_text, qcCoreFrozenLipidItsd$text[[i]], info = pattern)
        expect_identical(sprintf("Found %d internal standards", core$total),
            sub("^\\[message\\] ", "", qcCoreFrozenLipidItsd$note))
    }
    # The refusals, with their texts: no standard found (the tip names lipid's column), no pattern at all.
    cases <- list(list(object = module_ci_lipid_qc_object(layout = "gc", include_itsd = FALSE), pattern = "^IS_",
        note = qcCoreFrozenLipidItsd$refusal_no_standard),
        list(object = module_ci_lipid_qc_object(layout = "gc", internal_standard_regex = NA_character_), pattern = "",
            note = qcCoreFrozenLipidItsd$refusal_no_pattern))
    for (case in cases) {
        workflow <- module_ci_lipid_qc_workflow(case$object)
        message <- qcCoreError(analyzeWorkflowOmicsItsd(lipidomicsDeclarationBinding, workflow$state_manager, case$pattern, itsd))
        expect_false(is.null(message))
        expect_identical(paste("[error] Error analyzing internal standards:", message), case$note)
    }
    expect_match(message, "No internal standard pattern provided", fixed = TRUE)
})

test_that("F104 (W2): a NULL pattern is read as an empty one, as metab reads it (lipid's inline analysis errors) (frozen from 7b99e7e8, KIT-116)", {
    itsd <- omicGuiConfig("lipidomics")$qc$itsd
    workflow <- module_ci_lipid_qc_workflow(module_ci_lipid_qc_object(layout = "combined"))
    # Lipid's inline analysis failed on a NULL input (nzchar(NULL)), notifying "[error] Error analyzing internal
    # standards: argument is of length zero" (checked by q64-qc.R at 7b99e7e8); only a test session leaves the input
    # NULL. The core reads it as an empty pattern.
    expect_identical(analyzeWorkflowOmicsItsd(lipidomicsDeclarationBinding, workflow$state_manager, NULL, itsd),
        analyzeWorkflowOmicsItsd(lipidomicsDeclarationBinding, workflow$state_manager, "", itsd))
})

# Metab's ITSD analysis (analyzeMetabQcItsdData on readMetabQcItsdState), frozen from 827f2328 (KIT-116): the same
# metrics, long data and per-assay metrics (sha256) and total for every pattern; the pattern it used and its text (an
# empty or NULL pattern used "^IS_"); and the three refusals' texts.
qcCoreFrozenMetabItsdText <- function(pattern) {
    paste0("Internal Standard Analysis Complete\n====================================\nPattern used: ", pattern,
        "\nTotal IS detected: 4\n\nCV Statistics:\n  Median CV: 0.9%\n  Max CV: 0.9%\n  IS with CV > 30%: 0\n\n",
        "Per-Assay IS Counts:\n  LCMS_Pos: 2 internal standards\n  GCMS: 2 internal standards")
}
qcCoreFrozenMetabItsd <- list(
    metrics = "4eb7294989b0e74ae594757e9f48804d1634b6e4378b351693afb742b4e7b790",
    long = "8a5f51eadfa26985695b67d26b5269707bd469c97d10118787e4b410f5a48f93",
    by_assay = "a7e4b7e7cebead1695fada1db483742ac50dd17b1c67be10b2df21a3dffad240",
    total = 4L,
    pattern = c("^IS_", "^IS_", "^NO_MATCH$", "^IS_"),
    refusal_no_standard = paste("No internal standards found. Pattern: ^IS_",
        "Columns searched: database_identifier, Name, metabolite",
        "Tip: Check if your data has a column with metabolite names (e.g., 'Metabolite name')", sep = "\n"),
    refusal_no_pattern = "No internal standard pattern provided. Please enter a regex pattern to identify IS.",
    refusal_other_omic = "no slot of name \"metabolite_id_column\" for this object of class \"LipidomicsAssayData\"")

test_that("analyzeWorkflowOmicsItsd gives metab's analysis (analyzeMetabQcItsdData on readMetabQcItsdState) (frozen from 827f2328, KIT-116)", {
    itsd <- omicGuiConfig("metabolomics")$qc$itsd
    patterns <- list("", "^IS_", "^NO_MATCH$", NULL)
    for (i in seq_along(patterns)) {
        pattern <- patterns[[i]]
        workflow <- module_ci_metab_qc_workflow(module_ci_metab_qc_object(layout = "combined"))
        manager <- workflow$state_manager
        core <- analyzeWorkflowOmicsItsd(metabolomicsDeclarationBinding, manager, pattern, itsd)
        info <- paste("pattern", format(pattern))
        expect_identical(qcCoreDigest(core$metrics), qcCoreFrozenMetabItsd$metrics, info = info)
        expect_identical(qcCoreDigest(core$long_data), qcCoreFrozenMetabItsd$long, info = info)
        expect_identical(core$result_text, qcCoreFrozenMetabItsdText(qcCoreFrozenMetabItsd$pattern[[i]]), info = info)
        expect_identical(core$total, qcCoreFrozenMetabItsd$total, info = info)
        expect_identical(core$pattern, qcCoreFrozenMetabItsd$pattern[[i]], info = info)
        expect_identical(qcCoreDigest(core$metrics_by_assay), qcCoreFrozenMetabItsd$by_assay, info = info)
    }
    for (case in list(list(object = module_ci_metab_qc_object(layout = "gc", include_itsd = FALSE), pattern = "^IS_",
        today = qcCoreFrozenMetabItsd$refusal_no_standard),
        list(object = module_ci_metab_qc_object(layout = "gc", internal_standard_regex = NA_character_), pattern = "",
            today = qcCoreFrozenMetabItsd$refusal_no_pattern))) {
        manager <- module_ci_metab_qc_workflow(case$object)$state_manager
        expect_identical(qcCoreError(analyzeWorkflowOmicsItsd(metabolomicsDeclarationBinding, manager, case$pattern, itsd)),
            case$today)
    }
    # Another omic's state is refused as today (the shell reader fails on the binding's id slot).
    manager <- module_ci_lipid_qc_workflow()$state_manager
    expect_identical(qcCoreError(analyzeWorkflowOmicsItsd(metabolomicsDeclarationBinding, manager, "", itsd)),
        qcCoreFrozenMetabItsd$refusal_other_omic)
})

# The ITSD figures, frozen from 7b99e7e8 (KIT-116): the sha256 of metab's analysis (analyzeMetabQcItsdData) metrics
# and long data on the combined fixture, of qcCorePlotParts of buildMetabQcItsdCvPlot and
# buildMetabQcItsdIntensityPlot on them, and of the data lipid's registerLipidItsdCvPlotOutput drew.
qcCoreFrozenFigures <- list(
    metab_metrics = "4eb7294989b0e74ae594757e9f48804d1634b6e4378b351693afb742b4e7b790",
    metab_long = "8a5f51eadfa26985695b67d26b5269707bd469c97d10118787e4b410f5a48f93",
    cv_plot = "e6185c5e3acd9ace3ddb152fea8dc5d03d9d63ae58d50aba2da0baa2df625c41",
    intensity_plot = "43cb1d8c54688a54fd4e43ea1738c26c344bf133d6d225c1f8f2ac3e9d454e34",
    lipid_cv_data = "1f47d6d60a4c04449a1da75ec3c5a349eb5e981077c8feb7efa1b192ecc3c92e")

test_that("omicsItsdCvPlot and omicsItsdIntensityPlot are today's ITSD figures (frozen from 7b99e7e8, KIT-116)", {
    workflow <- module_ci_metab_qc_workflow(module_ci_metab_qc_object(layout = "combined"))
    # The figures' input: the core's metab analysis, the same values as the deleted metab analysis gave.
    analysis <- analyzeWorkflowOmicsItsd(metabolomicsDeclarationBinding, workflow$state_manager, "",
        omicGuiConfig("metabolomics")$qc$itsd)
    expect_identical(qcCoreDigest(analysis$metrics), qcCoreFrozenFigures$metab_metrics)
    expect_identical(qcCoreDigest(analysis$long_data), qcCoreFrozenFigures$metab_long)
    expect_identical(qcCoreDigest(qcCorePlotParts(omicsItsdCvPlot(analysis$metrics))), qcCoreFrozenFigures$cv_plot)
    expect_identical(qcCoreDigest(qcCorePlotParts(omicsItsdIntensityPlot(analysis$long_data))),
        qcCoreFrozenFigures$intensity_plot)
    # Lipid's inline CV figure: the data its output drew (one ggplot call) for the same metrics.
    lipid <- module_ci_lipid_qc_workflow(module_ci_lipid_qc_object(layout = "combined"))
    metrics <- analyzeWorkflowOmicsItsd(lipidomicsDeclarationBinding, lipid$state_manager, "",
        omicGuiConfig("lipidomics")$qc$itsd)$metrics
    expect_identical(qcCoreDigest(omicsItsdCvPlot(metrics)$data), qcCoreFrozenFigures$lipid_cv_data)
})

# ------------------------------------------------------------------------------------------------
# The finalise.
# ------------------------------------------------------------------------------------------------

# validateLipidQcS4FinalizeState, frozen from 7b99e7e8 (KIT-116): the aligned shell passed and was returned unchanged
# (identical to its input; the shell carries its owner's per-run entity reference, so the frozen value is that
# identity, not a digest of the shell), and the five refusals' texts.
qcCoreFrozenLipidFinalise <- list(
    misaligned = paste("Current LipidomicsAssayData object has assay/design mismatch:",
        "Design matrix group column 'group' contains blank assignments."),
    zero_rows = "Current LipidomicsAssayData object has no assay data to finalize.",
    no_design = "Current LipidomicsAssayData object has no design matrix rows to finalize.",
    no_assays = "Current LipidomicsAssayData object has no assay data to finalize.",
    other_class = "Current state is not a LipidomicsAssayData object")

test_that("validateWorkflowOmicsQcFinalize is lipid's finalise validation: passes, refusals and their texts (frozen from 7b99e7e8, KIT-116)", {
    qc <- qcCoreConfig("lipidomics")
    aligned <- readLipidWorkflowShell(module_ci_lipid_qc_workflow()$state_manager, measures = "features")
    expect_identical(validateWorkflowOmicsQcFinalize(aligned, qc, lipidomicsDeclarationBinding), aligned)
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
        expect_identical(qcCoreError(validateWorkflowOmicsQcFinalize(cases[[name]], qc, lipidomicsDeclarationBinding)),
            qcCoreFrozenLipidFinalise[[name]], info = name)
    }
    expect_identical(qcCoreError(validateWorkflowOmicsQcFinalize(misaligned, qc, lipidomicsDeclarationBinding)),
        "Current LipidomicsAssayData object has assay/design mismatch: Design matrix group column 'group' contains blank assignments.")
})

# validateMetabQcS4FinalizeState, frozen from 827f2328 (KIT-116): the aligned shell and the zero-row shell passed and
# were returned unchanged (identical to their input; the shell carries its owner's per-run entity reference, so the
# frozen value is that identity, not a digest of the shell); the three refusals' texts; and a blank group passed
# (a MetaboliteAssayData was returned).
qcCoreFrozenMetabFinalise <- list(
    no_design = "Current MetaboliteAssayData object has no design matrix rows to finalize.",
    no_assays = "Current MetaboliteAssayData object has no assay data to finalize.",
    other_class = "Current state is not a MetaboliteAssayData object")

test_that("validateWorkflowOmicsQcFinalize for metab: today's passes and refusals, and the alignment (parity) (frozen from 827f2328, KIT-116)", {
    qc <- qcCoreConfig("metabolomics")
    aligned <- readMetabWorkflowShell(module_ci_metab_qc_workflow()$state_manager, measures = "features")
    expect_identical(validateWorkflowOmicsQcFinalize(aligned, qc, metabolomicsDeclarationBinding), aligned)
    # D2: metab refuses only an empty assay list; a zero-row assay passes, as today.
    zero_rows <- aligned
    attr(zero_rows, "assay_rows")[[1L]] <- 0L
    expect_identical(validateWorkflowOmicsQcFinalize(zero_rows, qc, metabolomicsDeclarationBinding), zero_rows)
    no_design <- aligned
    no_design@design_matrix <- no_design@design_matrix[0L, , drop = FALSE]
    no_assays <- aligned
    no_assays@metabolite_data <- list()
    attr(no_assays, "assay_rows") <- integer()
    other_class <- readLipidWorkflowShell(module_ci_lipid_qc_workflow()$state_manager)
    cases <- list(no_design = no_design, no_assays = no_assays, other_class = other_class)
    for (name in names(cases)) {
        expect_identical(qcCoreError(validateWorkflowOmicsQcFinalize(cases[[name]], qc, metabolomicsDeclarationBinding)),
            qcCoreFrozenMetabFinalise[[name]], info = name)
    }
    # Parity (D11): a blank group passed metab's finalise (frozen above) and is refused by the core with lipid's text.
    misaligned <- readMetabWorkflowShell(module_ci_metab_qc_workflow(qcCoreBlankGroup(module_ci_metab_qc_object()))$state_manager)
    expect_identical(qcCoreError(validateWorkflowOmicsQcFinalize(misaligned, qc, metabolomicsDeclarationBinding)),
        "Current MetaboliteAssayData object has assay/design mismatch: Design matrix group column 'group' contains blank assignments.")
    # Every omic's alignment is called: a configuration without one is refused, not skipped.
    qc$finalize_alignment <- NULL
    expect_error(validateWorkflowOmicsQcFinalize(aligned, qc, metabolomicsDeclarationBinding), class = "omic_gui_config_invalid")
})

test_that("omicsQcFinalizeAlignment is lipid's alignment call (validateLipidDesignAlignment with its arguments)", {
    for (object in list(module_ci_lipid_qc_object(), qcCoreBlankGroup(module_ci_lipid_qc_object()))) {
        shell <- readLipidWorkflowShell(module_ci_lipid_qc_workflow(object)$state_manager)
        runs <- as.character(shell@design_matrix[[shell@sample_id]])
        today <- validateLipidDesignAlignment(designMatrix = shell@design_matrix, assayList = shell@lipid_data,
            columnMapping = list(lipid_id_col = shell@lipid_id_column, annotation_col = shell@annotation_id_column,
                sample_columns = runs), sampleIdCol = shell@sample_id, groupCol = shell@group_id, requireReplicates = FALSE)
        expect_identical(omicsQcFinalizeAlignment(shell, "lipidomics"), today[c("valid", "errors")])
    }
    shell <- readMetabWorkflowShell(module_ci_metab_qc_workflow()$state_manager)
    expect_identical(omicsQcFinalizeAlignment(shell, "metabolomics"), list(valid = TRUE, errors = character()))
})

test_that("trxDesignAlignment passes on the trx_screen shell and refuses what validateTrxDesign refuses (frozen from 827f2328, KIT-116)", {
    skip_if_not(workflowOwnerCapable())
    fixture <- trxDeclarationsOwner(withr::local_tempdir(), module_ci_trx_corpus())
    withr::defer(closeDIAStateAdapter(fixture$owner))
    # The trx class's validity refuses a zero-row count table (W2 finding B: readWorkflowOmicsShell cannot read a trx
    # state), so the shell is the object with its zero-row prototype and stored row count, built without validity.
    shell <- fixture$object
    shell@count_data <- lapply(shell@count_data, function(table) table[0L, , drop = FALSE])
    attr(shell, "assay_rows") <- vapply(fixture$object@count_data, nrow, integer(1))
    expect_identical(trxDesignAlignment(shell), list(valid = TRUE, errors = character()))
    # The same check as every trx step (trxAssayInputs) on the whole object.
    expect_type(trxAssayInputs(fixture$object), "list")
    blank <- qcCoreBlankGroup(shell)
    expect_false(trxDesignAlignment(blank)$valid)
    expect_match(trxDesignAlignment(blank)$errors, "outside the declared levels", fixed = TRUE)
    expect_identical(trxDesignAlignment(blank)$errors, trxConditionText(tryCatch(trxAssayInputs(qcCoreBlankGroup(fixture$object)),
        trx_error = identity)))
    absent <- shell
    extra <- absent@design_matrix[1L, , drop = FALSE]
    extra[[absent@sample_id]] <- "S_absent"
    absent@design_matrix <- rbind(absent@design_matrix, extra)
    expect_match(trxDesignAlignment(absent)$errors, "S_absent", fixed = TRUE)

    # The finalise: today's trx checks (trxQcFinalizeCollaborators()$validateFinalizeStateFn, on the metadata) passed
    # the metadata unchanged and refused a zero-row assay with the text below (frozen from 827f2328, KIT-116); the
    # alignment is the parity.
    qc <- qcCoreConfig("transcriptomics")
    expect_identical(validateWorkflowOmicsQcFinalize(shell, qc, transcriptomicsDeclarationBinding), shell)
    zero_rows <- shell
    attr(zero_rows, "assay_rows")[[1L]] <- 0L
    expect_identical(qcCoreError(validateWorkflowOmicsQcFinalize(zero_rows, qc, transcriptomicsDeclarationBinding)),
        "Current TranscriptomicsAssayData object has no genes or no design rows to finalize.")
    expect_identical(qcCoreError(validateWorkflowOmicsQcFinalize(blank, qc, transcriptomicsDeclarationBinding)),
        paste("Current TranscriptomicsAssayData object has assay/design mismatch:", trxDesignAlignment(blank)$errors))
})

# ------------------------------------------------------------------------------------------------
# The QC action (D6).
# ------------------------------------------------------------------------------------------------

test_that("runWorkflowOmicsQcAction in process: today's reference step, then the figure from the state read for it", {
    workflow <- module_ci_lipid_qc_workflow()
    twin <- module_ci_lipid_qc_workflow()
    step <- omicsIntensityFilterStep(lipidomicsDeclarationBinding, 40L, 0.5)
    seen <- new.env(parent = emptyenv())
    runWorkflowOmicsQcAction(lipidomicsDeclarationBinding, workflow, step, deliver = function(result, figure) {
        seen$result <- result
        seen$history_at_delivery <- workflow$state_manager$getHistory()
        seen$figure <- figure()
    }, draw = function(arguments) arguments, failed = function(error) stop("not called in process"))
    today <- runLipidReferenceIntensityFilter(twin$state_manager, 40L, 0.5, twin$config_list)
    # The value's `final` is the owner's own product reference; the counts are the step's.
    expect_identical(names(seen$result$value), names(today$value))
    expect_identical(seen$result$value$counts, today$value$counts)
    expect_identical(workflow$state_manager$getHistory(), twin$state_manager$getHistory())
    expect_identical(seen$history_at_delivery, twin$state_manager$getHistory())
    expect_identical(seen$result$reference$entity_id, getWorkflowStateReference(workflow$state_manager)$entity_id)
    expect_named(seen$figure, "theObject")
    expect_identical(seen$figure$theObject@lipid_data, readLipidWorkflowState(twin$state_manager, today$reference)@lipid_data)
})

test_that("runWorkflowOmicsQcAction on a runner: the managed commit, then the figure from the worker's metrics", {
    seen <- new.env(parent = emptyenv())
    metrics <- new.env(parent = emptyenv())
    metrics$value <- list(metric = 1)
    runner <- list(run = function(actions, done = NULL, failed = NULL, dropped = NULL) {
        seen$actions <- actions
        seen$failed <- failed
        actions[[1L]]$deliver(list(value = list(counts = "counts"), stateReference = "reference"))
    })
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) {
        seen$omic_type <- omic_type
        runner
    }, readWorkflowOmicsFilteringMetrics = function(manager, result) metrics$value, .package = "MultiScholaR")
    workflow <- list2env(list(state_manager = "owner", config_list = list(config = TRUE)))
    step <- omicsDuplicateResolutionStep(metabolomicsDeclarationBinding)
    failed <- function(error) "failed"
    runWorkflowOmicsQcAction(metabolomicsDeclarationBinding, workflow, step, deliver = function(result, figure) {
        seen$result <- result
        seen$figure <- tryCatch(figure(), error = function(error) conditionMessage(error))
    }, draw = function(arguments) arguments, failed = failed)
    expect_identical(seen$omic_type, "metabolomics")
    expect_length(seen$actions, 1L)
    expected <- omicsManagedCommit(step, workflow$config_list, NULL)
    expect_identical(seen$actions[[1L]][c("operation", "parameters", "before", "context")], expected[c("operation", "parameters", "before", "context")])
    expect_identical(seen$failed, failed)
    expect_identical(seen$result, list(value = list(counts = "counts"), reference = "reference"))
    expect_identical(seen$figure, list(theObject = NULL, metrics = list(metric = 1)))
    # No metrics from the worker: the figure fails with today's message (the tab logs it and shows none).
    metrics$value <- NULL
    runWorkflowOmicsQcAction(metabolomicsDeclarationBinding, workflow, step, deliver = function(result, figure) {
        seen$figure <- tryCatch(figure(), error = function(error) conditionMessage(error))
    }, draw = function(arguments) arguments, failed = failed)
    expect_identical(seen$figure, "the worker computed no filtering metrics")
})

# ------------------------------------------------------------------------------------------------
# Trx's QC figure functions, moved with their names.
# ------------------------------------------------------------------------------------------------

test_that("trx's QC figure functions are in func_trx_qc_plots.R, text-identical to any module copy left", {
    definitions <- function(file) {
        exprs <- parse(testthat::test_path("..", "..", "R", file), keep.source = FALSE)
        found <- list()
        for (e in exprs) if (is.call(e) && identical(e[[1L]], as.name("<-")) && is.name(e[[2L]])) found[[as.character(e[[2L]])]] <- e[[3L]]
        found
    }
    moved <- c("trxQcRawStats", "trxQcFilterText", "trxQcMetricsTable", "trxQcFilteringFigure")
    core <- definitions("func_trx_qc_plots.R")
    expect_true(all(moved %in% names(core)))
    module_file <- testthat::test_path("..", "..", "R", "mod_trx_qc_filter.R")
    module <- if (file.exists(module_file)) definitions("mod_trx_qc_filter.R") else list()
    for (name in intersect(moved, names(module))) expect_identical(core[[name]], module[[name]], info = name)
    counts <- list(before = 10L, after_filter_by_expr = 8L, after_zero_variance = 7L)
    expect_match(trxQcFilterText(counts), "After the Zero-variance step: 7", fixed = TRUE)
    metrics <- list(library_size = c(R1 = 100, R2 = 200), detected_genes = c(R1 = 5L, R2 = 6L))
    expect_identical(names(trxQcMetricsTable(metrics)), c("Run", "Library size", "Genes detected"))
    expect_s3_class(trxQcFilteringFigure(counts, metrics), "ggplot")
})
