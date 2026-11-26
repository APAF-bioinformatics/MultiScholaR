# KIT-280 (VAL-354, VAL-365): the proteomics missing-value filter counts a value at or below the intensity threshold as
# present (rule B, KIT-211), so the filter and the test agree on what a value is. The step records the replay marker
# `threshold_values_in_filter`; the filter step writes `low_intensity_groups` into `protein_id_table` and the DA writes
# `low_intensity_value` per comparison; the three proteomics templates state the rule. Synthetic data only.
#
# Fixture F: 40 proteins x 6 samples (control c1 to c3, treated t1 to t3); 38 fillers, every value at least 20.00;
# P_low1: control (10.0, 10.1, 20.0), treated (20.0, 20.1, 20.2); P_low2: control (10.2, 20.0, 20.1), treated 20.3 to
# 20.5. n = 240, h = 3.39, the 1st percentile (type 7) is 10.2 + 0.39 x (20.0 - 10.2) = 14.022, and exactly three values
# are at or below it. Minimum replicates 2, minimum groups 2: the earlier rule removes P_low1 (two of three control
# values missing, one allowed); P_low2 stays. The fillers carry a level per protein (20 + (i - 1) / 10) so that the
# samples correlate for the route.
#
# The draft sentences are file-local constants: the owner's wording pass is a text edit here and in the templates
# (R/func_prot_qc_filtering_helpers.R, proteinThresholdValuesSentence).

k280Quiet <- function(expr) {
    value <- NULL
    utils::capture.output(value <- suppressWarnings(suppressMessages(expr)))
    value
}

k280Design <- function(sizes) {
    groups <- rep(names(sizes), sizes)
    runs <- paste0(groups, "_", sequence(sizes))
    data.frame(Run = runs, group = groups, replicates = runs, stringsAsFactors = FALSE)
}

# `special`: named list, protein to its full row of values (NA allowed).
k280Table <- function(design, special = list(), fillers = 38L) {
    ids <- c(sprintf("P%03d", seq_len(fillers)), names(special))
    values <- 20 + outer(seq_along(ids), seq_along(design$Run), \(i, j) (i - 1) / 10 + ((7 * i + 3 * j) %% 10) / 100)
    dimnames(values) <- list(ids, design$Run)
    for (id in names(special)) values[id, ] <- special[[id]]
    data.frame(Protein.Ids = ids, values, row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
}

k280Object <- function(table, design) {
    ProteinQuantitativeData(protein_quant_table = table, design_matrix = design, sample_id = "Run", group_id = "group",
        technical_replicate_id = "replicates", protein_id_column = "Protein.Ids",
        protein_id_table = data.frame(Protein.Ids = table$Protein.Ids, stringsAsFactors = FALSE),
        args = list(globalParameters = list(workflow_type = "LFQ")))
}

k280FSpecial <- list(P_low1 = c(10.0, 10.1, 20.0, 20.0, 20.1, 20.2), P_low2 = c(10.2, 20.0, 20.1, 20.3, 20.4, 20.5))
k280FDesign <- function() k280Design(c(control = 3L, treated = 3L))
k280FObject <- function() k280Object(k280Table(k280FDesign(), k280FSpecial), k280FDesign())
k280FIds <- c(sprintf("P%03d", 1:38), "P_low1", "P_low2")

k280Owner <- function(root, type = "DIA") {
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, paste0("k280-", basename(root)), type),
        WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType(type)
    paths <- list(source_dir = file.path(root, "src"), protein_qc_dir = file.path(root, "pqc"),
        da_output_dir = file.path(root, "da"), publication_graphs_dir = file.path(root, "pub"), time_dir = file.path(root, "time"),
        results_dir = file.path(root, "res"))
    for (p in paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
    manager$.context$experiment_paths <- paths
    manager
}

# The filter as the tab applies it, on a saved protein state: the declared kernel, then the commit.
k280Apply <- function(object, root = withr::local_tempdir(.local_envir = parent.frame()), percentile = 1, min_reps = 2L,
    min_groups = 2L, strict = FALSE, type = "DIA") {
    manager <- k280Owner(root, type)
    force(manager)
    withr::defer(closeDIAStateAdapter(manager), envir = parent.frame())
    manager$saveState("protein_accession_cleaned", object, list(), "fixture")
    workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list())))
    value <- k280Quiet(computeDIAReferenceProteinIntensityFilter(manager, getDIAStateReference(manager),
        list(strict_mode = strict, min_reps_per_group = min_reps, min_groups = min_groups,
            proteins_intensity_cutoff_percentile = percentile)))
    k280Quiet(finishDIAReferenceProteinStep(workflow, value, publish = FALSE))
    list(manager = manager, workflow = workflow, value = value, filtered = manager$getState("protein_intensity_filtered"))
}

k280Ids <- function(object) as.character(object@protein_quant_table$Protein.Ids)

k280MockPlots <- function(env = parent.frame()) {
    testthat::local_mocked_bindings(plotRle = function(...) ggplot2::ggplot(), plotPca = function(...) ggplot2::ggplot(),
        plotNumOfValuesNoLog = function(...) ggplot2::ggplot(), .package = "MultiScholaR", .env = env)
}

# differentialAbundanceAnalysis() on a saved state's object (limma route).
k280Da <- function(object, contrasts) {
    tryCatch(list(value = k280Quiet(differentialAbundanceAnalysis(object,
        contrasts_tbl = data.frame(contrasts = contrasts, stringsAsFactors = FALSE), formula_string = "~ 0 + group",
        da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, eBayes_trend = TRUE, eBayes_robust = TRUE,
        args_row_id = "Protein.Ids", qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue"))),
        error = function(e) list(error = conditionMessage(e)))
}

# The call completed (a failing call fails here and prints the caught message).
k280Completed <- function(run, label) {
    expect_null(run$error, label = label)
    is.null(run$error)
}

k280TvC <- "treated_vs_control=grouptreated-groupcontrol"

# The route after the filter: duplicate removal, replicate filter, normalisation (RUV skipped or manual), correlation
# filter and the declared DA (as the tabs run them). Returns the long table and the DA manifest.
k280Route <- function(manager, ruv = c("skip", "manual"), contrast = "treated_vs_control") {
    ruv <- match.arg(ruv)
    workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list()), config_list = list()))
    paths <- manager$.context$experiment_paths
    parameters <- list(norm_method = "quantile", ruv_mode = ruv, auto_percentage_min = NULL, auto_percentage_max = NULL,
        separation_metric = NULL, k_penalty_weight = NULL, adaptive_k_penalty = NULL,
        ruv_percentage = if (ruv == "manual") 50 else NULL, ruv_k = if (ruv == "manual") 1 else NULL, grouping = "group")
    k280Quiet(withDIAStateScope(manager, function() {
        finish <- function(value) finishDIAReferenceProteinStep(workflow, value, publish = FALSE)
        finish(computeDIAReferenceProteinDuplicateRemoval(manager, getDIAStateReference(manager),
            list(aggregation_method = "mean")))
        completeDIAProteinAction(workflow, computeDIAReferenceProteinReplicateFilter(manager,
            getDIAStateReference(manager), list(grouping_variable = "group", parallel_cores = 1L), variant = "rows"),
            publish = FALSE)
        run <- computeDIAReferenceNormalization(manager, getDIAStateReference(manager), parameters)
        commit <- diaNormalizationCommit(manager, parameters, run)
        if (!is.null(commit)) commitDIAReferenceState(manager, commit$state_name, commit$checkpoint, commit$config,
            commit$description)
        correlated <- computeDIAReferenceCorrelation(manager, run$value$final, "group", 0.5)
        commit <- diaCorrelationCommit(manager, correlated$final, 0.5, FALSE)
        if (!is.null(commit)) commitDIAReferenceState(manager, commit$state_name, commit$checkpoint, commit$config,
            commit$description)
        contrasts <- data.frame(contrasts = contrast,
            full_format = sprintf("%s=grouptreated-groupcontrol", contrast), friendly_names = contrast,
            stringsAsFactors = FALSE)
        da_data <- new.env(parent = emptyenv())
        da_data$current_s4_object <- getDIAStateReference(manager)
        da_data$contrasts_available <- contrasts$contrasts
        workflow$contrasts_tbl <- contrasts
        manifest <- runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05,
            treat_lfc_cutoff = 0), da_data, workflow, paths, contrasts = contrasts, publish = FALSE)
        tables <- readDIAContrastOutput(manifest$results[[1L]], c("da_proteins_long", "contrasts_results_table"))
        list(long = as.data.frame(tables$da_proteins_long), manifest = manifest, final = manager$getState(manager$current_state))
    }))
}

# ---- T1, T2, T3: fixture F --------------------------------------------------------------------------------------------

test_that("T1: rule B keeps every protein of F (the earlier rule removed P_low1)", {
    run <- k280Apply(k280FObject())
    expect_setequal(k280Ids(run$filtered), k280FIds)
    expect_identical(run$value$protein_count, 40L)
})

test_that("T2: P_low1 reaches the DA table with a finite p value and the logFC of its own values", {
    k280MockPlots()
    run <- k280Apply(k280FObject())
    da <- k280Da(run$filtered, k280TvC)
    if (!k280Completed(da, "T2 DA")) return(invisible())
    row <- da$value$da_proteins_long[da$value$da_proteins_long$Protein.Ids == "P_low1", , drop = FALSE]
    expect_identical(nrow(row), 1L)
    expect_equal(row$log2FC, 20.1 - 40.1 / 3, tolerance = 1e-9)
    expect_true(is.finite(row$raw_pvalue))
})

test_that("T3: low_intensity_value, treated against control, is TRUE for P_low1 and P_low2 and FALSE for the other 38", {
    k280MockPlots()
    run <- k280Apply(k280FObject())
    da <- k280Da(run$filtered, k280TvC)
    if (!k280Completed(da, "T3 DA")) return(invisible())
    long <- da$value$da_proteins_long
    expect_setequal(long$Protein.Ids[long$low_intensity_value %in% TRUE], c("P_low1", "P_low2"))
    expect_identical(sum(long$low_intensity_value %in% FALSE), 38L)
    expect_identical(nrow(long), 40L)
    pit <- run$filtered@protein_id_table
    expect_identical(pit$low_intensity_groups[pit$Protein.Ids == "P_low1"], "control")
    expect_identical(pit$low_intensity_groups[pit$Protein.Ids == "P_low2"], "control")
    expect_identical(sum(nzchar(pit$low_intensity_groups)), 2L)
})

test_that("T3b: the flag reaches the DA through duplicate removal, the replicate filter, normalisation and correlation (RUV skip and RUV manual, two chains)", {
    for (mode in c("skip", "manual")) {
        run <- k280Apply(k280FObject())
        result <- tryCatch(k280Route(run$manager, mode), error = function(e) list(error = conditionMessage(e)))
        expect_null(result$error, label = paste("T3b route", mode))
        if (!is.null(result$error)) next
        long <- result$long
        expect_true("low_intensity_value" %in% names(long), info = mode)
        expect_setequal(long$Protein.Ids[long$low_intensity_value %in% TRUE], c("P_low1", "P_low2"))
        expect_identical(sum(long$low_intensity_value %in% FALSE), 38L, info = mode)
    }
})

test_that("T3c: the threshold is taken from the pre-filter table (the filter removes the protein holding the minimum)", {
    # Pmin (control 5.0, then two missing; treated 20.0, then two missing) fails both groups and is removed. Before the
    # filter n = 242, h = 3.41 and the threshold is 11.5 + 0.41 x (20.0 - 11.5) = 14.985; after it n = 240 and the
    # threshold would be 20.0 or a little above, which would mark P_b's 20.0.
    design <- k280FDesign()
    special <- list(P_a = c(11.0, 11.5, 20.5, 20.6, 20.7, 20.8), P_b = c(20.9, 21.0, 21.1, 20.0, 21.2, 21.3),
        Pmin = c(5.0, NA, NA, 20.0, NA, NA))
    run <- k280Apply(k280Object(k280Table(design, special, fillers = 38L), design))
    expect_false("Pmin" %in% k280Ids(run$filtered))
    pit <- run$filtered@protein_id_table
    expect_identical(pit$low_intensity_groups[pit$Protein.Ids == "P_a"], "control")
    expect_identical(pit$low_intensity_groups[pit$Protein.Ids == "P_b"], "")
    # keepNA = TRUE: nzchar() is TRUE for NA under keepNA = NA, so the removed protein's NA would be counted.
    expect_identical(sum(nzchar(pit$low_intensity_groups, keepNA = TRUE), na.rm = TRUE), 1L)
    expect_true(is.na(pit$low_intensity_groups[pit$Protein.Ids == "Pmin"]))
})

test_that("T3d: digit-leading group labels (1h, 2h) are matched by the design's label in the comparison", {
    k280MockPlots()
    design <- k280Design(c(`1h` = 3L, `2h` = 3L))
    special <- list(P_low1 = c(10.0, 10.1, 20.0, 20.0, 20.1, 20.2), P_low2 = c(20.3, 20.4, 20.5, 10.2, 20.0, 20.1))
    run <- k280Apply(k280Object(k280Table(design, special), design))
    pit <- run$filtered@protein_id_table
    expect_identical(pit$low_intensity_groups[pit$Protein.Ids == "P_low1"], "1h")
    expect_identical(pit$low_intensity_groups[pit$Protein.Ids == "P_low2"], "2h")
    da <- k280Da(run$filtered, "h2_vs_h1=group2h-group1h")
    if (!k280Completed(da, "T3d DA")) return(invisible())
    long <- da$value$da_proteins_long
    expect_setequal(long$Protein.Ids[long$low_intensity_value %in% TRUE], c("P_low1", "P_low2"))
    expect_identical(sum(long$low_intensity_value %in% FALSE), 38L)
})

test_that("T4: three groups, a protein low only in g3: FALSE treated against control, TRUE g3 against control", {
    k280MockPlots()
    design <- k280Design(c(control = 3L, treated = 3L, g3 = 3L))
    special <- list(P_g3 = c(20.0, 20.1, 20.2, 20.3, 20.4, 20.5, 10.0, 10.1, 10.2))
    run <- k280Apply(k280Object(k280Table(design, special, fillers = 19L), design))
    da <- k280Da(run$filtered, c(k280TvC, "g3_vs_control=groupg3-groupcontrol"))
    if (!k280Completed(da, "T4 DA")) return(invisible())
    long <- da$value$da_proteins_long
    at <- function(comparison) long$low_intensity_value[long$comparison == comparison & long$Protein.Ids == "P_g3"]
    expect_identical(at("treated_vs_control"), FALSE)
    expect_identical(at("g3_vs_control"), TRUE)
})

# ---- T5, T5b: the replay marker --------------------------------------------------------------------------------------

k280Marker <- function(record, state, name = "threshold_values_in_filter") {
    step <- Filter(\(s) identical(s$state, state), record$steps)
    if (!length(step)) return(NULL)
    rows <- Filter(\(p) identical(p$name, name), step[[1L]]$parameters)
    if (!length(rows)) return(NULL)
    list(source = rows[[1L]]$source, value = decodeWorkflowRecordValue(rows[[1L]]$value))
}

k280Earlier <- list(MultiScholaR.missing_filter_rule = "rounded_percentage",
    MultiScholaR.missing_filter_threshold_rule = "missing")
k280CountsOnly <- list(MultiScholaR.missing_filter_threshold_rule = "missing")
k280SkipParameters <- list(norm_method = "quantile", ruv_mode = "skip", auto_percentage_min = NULL,
    auto_percentage_max = NULL, separation_metric = NULL, k_penalty_weight = NULL, adaptive_k_penalty = NULL,
    ruv_percentage = NULL, ruv_k = NULL, grouping = "group")

# An LFQ owner holding F: the filter made under `filter_options` (a step without fixed constants when `plain`), then
# optionally the normalisation (RUV skipped) under the current rules.
k280ReplayOwner <- function(root, filter_options = list(), plain = FALSE, normalise = FALSE) {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, basename(root), "LFQ"), WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("LFQ")
    setDIAStepProvenance(manager, "design_builder", fixed = list(pivot = "protein x run", transform = "log2(x + 1)"))
    manager$saveState("protein_s4_initial", k280FObject(), list(), "design")
    entry <- workflowOperationEntry("protein_intensity_filter")
    parameters <- list(strict_mode = FALSE, min_reps_per_group = 2L, min_groups = 2L, proteins_intensity_cutoff_percentile = 1)
    filter <- function() withDIAStateScope(manager, function() {
        value <- computeDIAReferenceProteinIntensityFilter(manager, getDIAStateReference(manager), parameters)
        setDIAStepProvenance(manager, "protein_intensity_filter", "whole", entry$parameters,
            if (!plain) entry$declaration$fixed, values = parameters)
        commitDIAReferenceState(manager, value$state_name, value$checkpoint, value$config, value$description)
    })
    normalisation <- function() withDIAStateScope(manager, function() {
        run <- computeDIAReferenceNormalization(manager, getDIAStateReference(manager), k280SkipParameters)
        commit <- diaNormalizationCommit(manager, k280SkipParameters, run)
        commitDIAReferenceState(manager, commit$state_name, commit$checkpoint, commit$config, commit$description)
    })
    utils::capture.output(suppressWarnings(suppressMessages({
        withr::with_options(filter_options, filter())
        if (normalise) normalisation()
    })))
    manager
}

# Write the owner's record, prepare a copy as the harness does and replay it (as test-prot-missing-value-counts.R).
k280Replay <- function(manager, root, copy, scratch) {
    json <- file.path(root, "analysis_record.json")
    writeWorkflowAnalysisRecord(manager, json)
    record <- jsonlite::read_json(json, simplifyVector = FALSE)
    states <- vapply(record$steps, \(step) step$state, "")
    later <- stats::setNames(states[-(1:2)], states[-(1:2)])
    origin <- withDIAStateScope(manager, function() lapply(later, \(state) manager$getState(state)))
    stopifnot(system2("cp", c("-a", file.path(root, "."), copy)) == 0L)
    h <- new.env(parent = globalenv())
    sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), h)
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    prepared <- h$kitPrepareReplayProject(copy, c("initial", "protein_s4_initial"), api)
    replayed <- list()
    replayed_record <- NULL
    utils::capture.output(replay <- suppressMessages(suppressWarnings(replayWorkflowAnalysisRecord(record, prepared,
        character(), list(), scratch, on_step = function(step, replay_manager) {
            replayed[[step$state]] <<- withDIAStateScope(replay_manager, function() replay_manager$getState(step$state))
            replayed_record <<- buildWorkflowAnalysisRecord(replay_manager)
        }))))
    list(record = record, replayed_record = replayed_record, states = states, replay_states = replay$states,
        origin = origin, replayed = replayed)
}

k280ExpectReplayIdentical <- function(result) {
    expect_identical(result$replay_states, result$states)
    for (state in names(result$origin)) expect_identical(result$replayed[[state]], result$origin[[state]], info = state)
}

k280HasFlag <- function(object) "low_intensity_groups" %in% names(object@protein_id_table)

test_that("T5: the marker row, a step made now, and records made before replay F with the rule each was made under", {
    marker <- workflowReplayRuleMarkers$threshold_values_in_filter
    expect_identical(marker$option, "MultiScholaR.missing_filter_threshold_rule")
    expect_identical(marker$current, "present")
    expect_identical(marker$earlier, "missing")
    expect_true("threshold_values_in_filter" %in% names(workflowOperationEntry("protein_intensity_filter")$declaration$fixed))
    expect_true("threshold_values_in_filter" %in% names(workflowOperationEntry("normalization-ruv")$declaration$fixed))
    filter_state <- "protein_intensity_filtered"
    kept <- function(result) sort(as.character(result$origin[[filter_state]]@protein_quant_table$Protein.Ids))
    # (a) a step made now lists the marker; it replays with the rule and writes the flag.
    root <- withr::local_tempdir()
    owner <- k280ReplayOwner(root)
    withr::defer(closeDIAStateAdapter(owner))
    a <- k280Replay(owner, root, withr::local_tempdir(), withr::local_tempdir())
    expect_identical(k280Marker(a$record, filter_state), list(source = "fixed", value = "present"))
    expect_identical(k280Marker(a$replayed_record, filter_state), list(source = "fixed", value = "present"))
    expect_identical(kept(a), sort(k280FIds))
    expect_true(k280HasFlag(a$origin[[filter_state]]))
    k280ExpectReplayIdentical(a)
    # (b) a record with only missing_value_comparison = counts replays with the earlier threshold rule.
    root <- withr::local_tempdir()
    owner <- k280ReplayOwner(root, filter_options = k280CountsOnly)
    withr::defer(closeDIAStateAdapter(owner))
    b <- k280Replay(owner, root, withr::local_tempdir(), withr::local_tempdir())
    expect_identical(k280Marker(b$record, filter_state, "missing_value_comparison"), list(source = "fixed", value = "counts"))
    expect_null(k280Marker(b$record, filter_state))
    expect_null(k280Marker(b$replayed_record, filter_state))
    expect_length(kept(b), 39L)
    expect_false("P_low1" %in% kept(b))
    expect_false(k280HasFlag(b$origin[[filter_state]]))
    expect_false(k280HasFlag(b$replayed[[filter_state]]))
    k280ExpectReplayIdentical(b)
    # (c) a record without any marker replays with both earlier rules.
    root <- withr::local_tempdir()
    owner <- k280ReplayOwner(root, filter_options = k280Earlier, plain = TRUE)
    withr::defer(closeDIAStateAdapter(owner))
    c_ <- k280Replay(owner, root, withr::local_tempdir(), withr::local_tempdir())
    expect_null(k280Marker(c_$record, filter_state, "missing_value_comparison"))
    expect_null(k280Marker(c_$record, filter_state))
    # Both earlier rules: the earlier comparison also removes P_low2 (one missing value of three is 33.3333 % against the
    # cutoff 33.333 %), so 38 proteins stay.
    expect_length(kept(c_), 38L)
    expect_false(any(c("P_low1", "P_low2") %in% kept(c_)))
    expect_false(k280HasFlag(c_$replayed[[filter_state]]))
    k280ExpectReplayIdentical(c_)
})

test_that("T5b: a session saved before this version and continued has no flag; each step keeps its own marker", {
    root <- withr::local_tempdir()
    owner <- k280ReplayOwner(root, filter_options = k280CountsOnly, normalise = TRUE)
    withr::defer(closeDIAStateAdapter(owner))
    result <- k280Replay(owner, root, withr::local_tempdir(), withr::local_tempdir())
    filter_state <- "protein_intensity_filtered"
    norm_state <- "normalized"
    expect_null(k280Marker(result$record, filter_state))
    expect_false(k280HasFlag(result$origin[[filter_state]]))
    expect_false(k280HasFlag(result$origin[[norm_state]]))
    expect_identical(k280Marker(result$record, norm_state), list(source = "fixed", value = "present"))
    k280ExpectReplayIdentical(result)
})

# ---- T6: strict mode and the clean-up after RUV ----------------------------------------------------------------------

test_that("T6: strict mode and the clean-up after RUV keep P_low1 (only an absent value counts as missing)", {
    run <- k280Apply(k280FObject(), strict = TRUE)
    expect_setequal(k280Ids(run$filtered), k280FIds)
    args <- list(removeRowsWithMissingValuesPercent = list(groupwise_percentage_cutoff = 33.333,
        max_groups_percentage_cutoff = 0, proteins_intensity_cutoff_percentile = 1, ruv_grouping_variable = "group"))
    corrected <- k280FObject()
    corrected@args <- c(corrected@args, args)
    filtered <- k280Quiet(removeRowsWithMissingValuesPercent(theObject = corrected))
    expect_setequal(k280Ids(filtered), k280FIds)
})

# ---- T7: routes and the tab's defaults -------------------------------------------------------------------------------

test_that("T7: the routes that hold the protein filter, and the tab's defaults (groups required stays 2)", {
    ids <- function(route) vapply(getProtQcProteinModuleSpecs(route), \(spec) spec$id, "")
    for (route in c("DIA", "LFQ", "TMT")) expect_true("intensity_filter" %in% ids(route), info = route)
    expect_false("intensity_filter" %in% ids("DIA_LIMPA"))
    html <- as.character(mod_prot_qc_protein_intensity_ui("k280"))
    expect_true(grepl("id=\"k280-min_groups\"[^>]*value=\"2\"", html))
    expect_true(grepl("id=\"k280-min_reps_per_group\"[^>]*value=\"2\"", html))
})

# ---- T8: the Apply text, tables without the flag ----------------------------------------------------------------------

test_that("T8: the Apply text states the marked count; a long table from an older object gains no column and is otherwise identical", {
    k280MockPlots()
    run <- k280Apply(k280FObject())
    expect_true(grepl("Proteins kept with a value at or below the threshold: 2 (marked in the results)", run$value$text,
        fixed = TRUE))
    plain <- k280Da(k280FObject(), k280TvC)
    flagged_object <- k280FObject()
    flagged_object@protein_id_table$low_intensity_groups <- ifelse(flagged_object@protein_id_table$Protein.Ids == "P_low1",
        "control", "")
    flagged <- k280Da(flagged_object, k280TvC)
    if (!k280Completed(plain, "T8 plain DA") || !k280Completed(flagged, "T8 flagged DA")) return(invisible())
    long_plain <- plain$value$da_proteins_long
    long_flagged <- flagged$value$da_proteins_long
    expect_false(any(c("low_intensity_groups", "low_intensity_value") %in% names(long_plain)))
    expect_identical(setdiff(names(long_flagged), names(long_plain)), c("low_intensity_groups", "low_intensity_value"))
    expect_identical(long_flagged[, names(long_plain)], long_plain)
    expect_identical(long_flagged$Protein.Ids[long_flagged$low_intensity_value %in% TRUE], "P_low1")
})

test_that("T8b: a saved DA table pages with the flag column when it holds one, and as before when it does not", {
    for (token in c("present", "missing")) {
        run <- withr::with_options(list(MultiScholaR.missing_filter_threshold_rule = token), k280Apply(k280FObject()))
        result <- tryCatch(k280Route(run$manager, "skip"), error = function(e) list(error = conditionMessage(e)))
        expect_null(result$error, label = paste("T8b route", token))
        if (!is.null(result$error)) next
        page <- withDIAStateScope(run$manager, function() k280Quiet(readDIAResultPage(result$manifest, "treated_vs_control",
            manager = run$manager)))
        expect_true(all(c("Protein.Ids", "log2FC", "raw_pvalue", "fdr_qvalue") %in% names(page$data)), info = token)
        expect_identical("low_intensity_value" %in% names(page$data), identical(token, "present"), info = token)
    }
})

# ---- T9: DIA with limpa never reaches the kernel ---------------------------------------------------------------------

test_that("T9: the DIA_LIMPA route's declared steps do not call the missing-value kernel", {
    for (name in c("computeDIAReferenceProteinDpcQuant", "computeDIAReferenceProteinNormalisation",
        "computeDIAReferenceProteinPeptideCountFilter")) {
        body_text <- paste(deparse(get(name, asNamespace("MultiScholaR"))), collapse = "\n")
        expect_false(grepl("removeRowsWithMissingValuesPercent", body_text, fixed = TRUE), info = name)
    }
})

# ---- T10: the report templates ---------------------------------------------------------------------------------------

k280Templates <- list(
    LFQ = list(file = "LFQ_report.rmd", extract = "ExtractStudyParametersLFQ", method = "FormatMethodDetailsLFQ"),
    TMT = list(file = "TMT_report.rmd", extract = "ExtractStudyParametersTMT", method = "FormatMethodDetailsTMT"),
    DIANN = list(file = "DIANN_report.rmd", extract = "ExtractStudyParametersNew", method = "FormatMethodDetails"))

k280Exprs <- function(path, label) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    start <- grep(paste0("^```\\{r ", label, "[,}]"), lines)
    stopifnot(length(start) == 1L)
    end <- start + which(lines[-seq_len(start)] == "```")[[1L]]
    chunk <- tempfile(fileext = ".Rmd")
    writeLines(lines[start:end], chunk, useBytes = TRUE)
    code <- tempfile(fileext = ".R")
    suppressMessages(knitr::purl(chunk, output = code, quiet = TRUE, documentation = 0))
    parse(code, keep.source = FALSE, encoding = "UTF-8")
}

k280Env <- function(name) {
    spec <- k280Templates[[name]]
    path <- system.file("reports", "proteomics", spec$file, package = "MultiScholaR")
    env <- new.env(parent = globalenv())
    env$params <- list(workflow_name = "maxquant", omic_type = "proteomics", experiment_label = "kit280")
    env$d66_log <- function(...) invisible(NULL)
    env$tibble <- tibble::tibble
    for (label in c("New Hierarchical Parser", "Method Section", "Result Section", "Sample Analysis Functions",
        "Main Execution Function")) {
        for (e in k280Exprs(path, label)) {
            if (is.call(e) && identical(e[[1L]], as.name("<-")) && is.call(e[[3L]]) &&
                identical(e[[3L]][[1L]], as.name("function"))) eval(e, env)
        }
    }
    env
}

k280Param <- function(name, value) list(name = name, value = list(type = typeof(value), value = as.list(value),
    length = length(value)))

k280Step <- function(operation, parameters = list()) {
    list(operation = operation, parameters = unname(Map(k280Param, names(parameters), parameters)), outcomes = list())
}

# The Methods text of one template for a file with no further sections and, optionally, a record.
k280Method <- function(name, steps = NULL) {
    spec <- k280Templates[[name]]
    env <- k280Env(name)
    dir <- withr::local_tempdir(pattern = "kit280-")
    env$report_source_dir <- dir
    params_path <- file.path(dir, "study_parameters.txt")
    writeLines(c("Study Parameters", "================", "", "Basic Information:", "-----------------",
        "* Workflow Name: kit280", "", "Protein Filtering Summary:", "-------------------------",
        "* Proteins after QC filtering: 1000", "* Proteins after RUV filtering: 950", ""), params_path, useBytes = TRUE)
    if (!is.null(steps)) jsonlite::write_json(list(steps = steps), file.path(dir, "analysis_record.json"),
        auto_unbox = TRUE, null = "null", digits = NA)
    study <- env$ReadStudyParamsNew(params_path)
    numbers <- env$AnalyzeSamples(NULL, study)$numbers
    stopifnot(length(numbers) > 0L)
    paste(as.character(env[[spec$method]](numbers)), collapse = "\n")
}

k280Rule <- paste0("A value at or below the %s percentile of all intensities in the table counted as present: ",
    "only an absent value counted as missing.")
k280NeutralSentence <- "Whether intensities at or below the %s percentile of all intensities counted as missing in this filter is not recorded for this analysis."
k280EarlierSentence <- "were treated as missing."
k280DiaSentence <- paste("On the DIA route the peptide intensity filter runs first, defaults to 2 groups and counts values",
    "below its own threshold as missing, so at the defaults this protein-level change removes almost no additional protein;",
    "the flag still marks proteins that hold a low value.")

k280Has <- function(text, pattern) grepl(pattern, text, fixed = TRUE)

test_that("T10: the Methods text states the rule the record says, the earlier rule without the marker, and neutral words without a record", {
    now <- list(k280Step("protein_intensity_filter", list(missing_value_comparison = "counts",
        threshold_values_in_filter = "present", proteins_intensity_cutoff_percentile = 5)))
    before <- list(k280Step("protein_intensity_filter", list(missing_value_comparison = "counts",
        proteins_intensity_cutoff_percentile = 5)))
    no_filter <- list(k280Step("normalization-ruv", list(ruv_mode = "automatic")))
    for (name in names(k280Templates)) {
        with_marker <- k280Method(name, now)
        expect_true(k280Has(with_marker, sprintf(k280Rule, "5th")), info = name)
        expect_true(k280Has(with_marker, "simulated data"), info = name)
        expect_true(k280Has(with_marker, "measured at the 1st percentile only"), info = name)
        expect_true(k280Has(with_marker, "low_intensity_value"), info = name)
        expect_false(k280Has(with_marker, k280EarlierSentence), info = name)
        expect_identical(k280Has(with_marker, k280DiaSentence), identical(name, "DIANN"), info = name)
        for (steps in list(before, no_filter)) {
            earlier <- k280Method(name, steps)
            expect_true(k280Has(earlier, k280EarlierSentence), info = name)
            expect_true(k280Has(earlier, "intensities in the table"), info = name)
            expect_false(k280Has(earlier, "counted as present"), info = name)
        }
        neutral <- k280Method(name)
        expect_true(k280Has(neutral, sprintf(k280NeutralSentence, "1st")), info = name)
        expect_false(k280Has(neutral, k280EarlierSentence), info = name)
        for (text in list(with_marker, neutral, k280Method(name, before))) {
            expect_true(k280Has(text, "exceeded the allowed share"), info = name)
        }
    }
    # DIA takes the percentile from the protein filter, not the peptide filter.
    expect_true(k280Has(k280Method("DIANN", before), "5th percentile of all"))
})

test_that("T10b: the percentile control says it marks values and does not remove proteins; the tab states the rule and the DIA sentence", {
    html <- as.character(mod_prot_qc_protein_intensity_ui("k280"))
    expect_true(grepl("Low-intensity marking percentile (%)", html, fixed = TRUE))
    expect_true(grepl("marking does not remove proteins", html, fixed = TRUE))
    expect_true(grepl("measured at 1 %", html, fixed = TRUE))
    expect_true(grepl("only an absent value counts as missing", html, fixed = TRUE))
    expect_true(grepl(k280DiaSentence, html, fixed = TRUE))
})
