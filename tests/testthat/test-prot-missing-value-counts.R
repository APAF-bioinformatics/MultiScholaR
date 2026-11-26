# KIT-187 (VAL-71, VAL-83): the protein missing-value filter decides by counts. Each stored cutoff is read as the
# three-decimal rounding of an exact share and compared in integers, so a group passes with at least min(n, m) values
# present and a protein with at least g passing groups. Tables are synthetic; every one holds `P_floor`, whose one value
# of 10 in the first run sets the threshold at percentile 0. Before KIT-280 present meant strictly above the threshold,
# so that value was missing; since KIT-280 a value counts as missing only when it is absent (the earlier rule is
# selected by `k187Legacy`).

k187Design <- function(sizes) {
    groups <- rep(names(sizes), sizes)
    runs <- paste0(groups, "_", sequence(sizes))
    data.frame(Run = runs, group = groups, replicates = runs, stringsAsFactors = FALSE)
}

# `missing`: named list, protein to the runs set to NA.
k187Table <- function(design, missing) {
    runs <- design$Run
    proteins <- c(names(missing), "P_floor")
    values <- 20 + outer(seq_along(proteins), seq_along(runs), \(i, j) ((7 * i + 3 * j) %% 10) / 100)
    dimnames(values) <- list(proteins, runs)
    for (protein in names(missing)) values[protein, missing[[protein]]] <- NA_real_
    values["P_floor", 1L] <- 10
    data.frame(Protein.Ids = proteins, values, check.names = FALSE, stringsAsFactors = FALSE)
}

k187Object <- function(table, design, args = list()) {
    ProteinQuantitativeData(protein_quant_table = table, design_matrix = design, sample_id = "Run", group_id = "group",
        technical_replicate_id = "replicates", protein_id_column = "Protein.Ids", args = args)
}

k187Quiet <- function(expr) {
    value <- NULL
    utils::capture.output(value <- suppressMessages(expr))
    value
}

# The tab's functions: updateMissingValueParameters, then the S4 method at percentile 0.
k187Filter <- function(object, min_reps, min_groups) {
    env <- new.env(parent = emptyenv())
    env$config_list <- list()
    updated <- k187Quiet(updateMissingValueParameters(object, min_reps_per_group = min_reps, min_groups = min_groups, env = env))
    filtered <- k187Quiet(removeRowsWithMissingValuesPercent(updated, proteins_intensity_cutoff_percentile = 0))
    list(object = updated, kept = as.character(filtered@protein_quant_table$Protein.Ids))
}

k187Helper <- function(table, design, groupwise, max_groups) {
    k187Quiet(removeRowsWithMissingValuesPercentHelper(table, cols = "Protein.Ids", design_matrix = design, sample_id = Run,
        row_id = Protein.Ids, grouping_variable = group, groupwise_percentage_cutoff = groupwise,
        max_groups_percentage_cutoff = max_groups, proteins_intensity_cutoff_percentile = 0))
}

k187T1Design <- function() k187Design(c(A = 3L, B = 3L))
k187T1Table <- function() k187Table(k187T1Design(), list(P_full = character(), P_one = "A_1", P_oneEach = c("A_1", "B_1"),
    P_two = c("A_1", "A_2")))
k187T1Kept <- c("P_full", "P_one", "P_oneEach", "P_floor")

# T3's grids. Per pair, TRUE when the count rule's outcome is observed: the protein with the allowed number of missing
# values (or failing groups) kept, the one with one more removed.
k187GridOutcome <- function() {
    outcome <- logical()
    for (n in 2:20) {
        design <- k187Design(stats::setNames(c(n, n), c("A", "B")))
        for (m in seq_len(n)) {
            allowed <- n - min(n, m)
            missing <- list(ok = if (allowed) paste0("A_", seq_len(allowed)) else character())
            if (n > m) missing$over <- paste0("A_", seq_len(allowed + 1L))
            kept <- k187Filter(k187Object(k187Table(design, missing), design), m, 2)$kept
            outcome[[sprintf("n%d_m%d", n, m)]] <- "ok" %in% kept && !"over" %in% kept
        }
    }
    for (G in 2:7) {
        groups <- paste0("G", seq_len(G))
        design <- k187Design(stats::setNames(rep(3L, G), groups))
        failing <- function(count) unlist(lapply(groups[seq_len(count)], \(g) paste0(g, "_", 1:3)))
        for (g in seq_len(G)) {
            missing <- list(ok = failing(G - g), over = failing(G - g + 1L))
            kept <- k187Filter(k187Object(k187Table(design, missing), design), 1, g)$kept
            outcome[[sprintf("G%d_g%d", G, g)]] <- "ok" %in% kept && !"over" %in% kept
        }
    }
    outcome
}

# The earlier rule restated in arithmetic (the order of operations of `updateMissingValueParameters` and of the
# helper's share): the pair fails when the exact allowed share exceeds its three-decimal rounding.
k187EarlierRuleFails <- function() {
    fails <- character()
    for (n in 2:20) for (m in seq_len(n)) {
        k <- n - min(n, m)
        if (k / n * 100 > round((k / n) * 100, 3)) fails <- c(fails, sprintf("n%d_m%d", n, m))
    }
    for (G in 2:7) for (g in seq_len(G)) {
        k <- G - g
        if (k / G * 100 > round((k / G) * 100, 3)) fails <- c(fails, sprintf("G%d_g%d", G, g))
    }
    fails
}

k187Manager <- function(object, root, state, type = "DIA") {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, basename(root), type), WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType(type)
    manager$saveState(state, object, list(), "fixture")
    manager
}

k187Marker <- function(record, state) {
    step <- Filter(\(s) identical(s$state, state), record$steps)
    if (!length(step)) return(NULL)
    rows <- Filter(\(p) identical(p$name, "missing_value_comparison"), step[[1L]]$parameters)
    if (!length(rows)) return(NULL)
    list(source = rows[[1L]]$source, value = decodeWorkflowRecordValue(rows[[1L]]$value))
}

# KIT-280: "the earlier rule" is the rule before both changes (the count comparison and the threshold-level values).
k187Legacy <- list(MultiScholaR.missing_filter_rule = "rounded_percentage",
    MultiScholaR.missing_filter_threshold_rule = "missing")

test_that("T1: groups of 3 at the tab's defaults keep a protein with one missing value per group", {
    result <- k187Filter(k187Object(k187T1Table(), k187T1Design()), 2, 2)
    expect_setequal(result$kept, k187T1Kept)
    expect_identical(result$object@args$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff, 33.333)
    expect_identical(result$object@args$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff, 0)
})

test_that("T2: three groups with two required keep a protein failing exactly one group", {
    design <- k187Design(c(A = 4L, B = 4L, C = 4L))
    table <- k187Table(design, list(P_full = character(), P_failOne = paste0("C_", 1:3),
        P_failTwo = c(paste0("B_", 1:3), paste0("C_", 1:3))))
    result <- k187Filter(k187Object(table, design), 2, 2)
    expect_setequal(result$kept, c("P_full", "P_failOne", "P_floor"))
    expect_identical(result$object@args$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff, 50)
    expect_identical(result$object@args$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff, 33.333)
})

test_that("T3: every balanced (n, m) up to 20 and every (G, g) up to 7 groups follows the count rule", {
    outcome <- k187GridOutcome()
    # 209 pairs (m from 1 to n; the 190 with m < n are 2.2's) and 27 (G, g) pairs.
    expect_length(outcome, 209L + 27L)
    # On the earlier rule the failure message is the failing set itself (acceptance criterion A2).
    expect_identical(names(outcome)[!outcome], character())
})

test_that("T4: unequal groups keep the largest group's share (3 and 6: unchanged), applied exactly (2 and 3)", {
    design <- k187Design(c(A = 3L, B = 6L))
    table <- k187Table(design, list(P_full = character(), P_twoA = c("A_1", "A_2")))
    expect_setequal(k187Filter(k187Object(table, design), 2, 2)$kept, c("P_full", "P_twoA", "P_floor"))
    design <- k187Design(c(A = 2L, B = 3L))
    table <- k187Table(design, list(P_full = character(), P_oneB = "B_1"))
    expect_setequal(k187Filter(k187Object(table, design), 2, 2)$kept, c("P_full", "P_oneB", "P_floor"))
})

test_that("T5: strict mode keeps only fully observed proteins", {
    manager <- k187Manager(k187Object(k187T1Table(), k187T1Design()), withr::local_tempdir(), "protein_accession_cleaned")
    withr::defer(closeDIAStateAdapter(manager))
    value <- k187Quiet(computeDIAReferenceProteinIntensityFilter(manager, getDIAStateReference(manager),
        list(strict_mode = TRUE, min_reps_per_group = 2L, min_groups = 2L, proteins_intensity_cutoff_percentile = 0)))
    expect_identical(value$protein_count, 2L)
    workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list())))
    k187Quiet(finishDIAReferenceProteinStep(workflow, value, publish = FALSE))
    expect_setequal(manager$getState("protein_intensity_filtered")@protein_quant_table$Protein.Ids, c("P_full", "P_floor"))
})

test_that("T6: explicit percentages give the same sets, and 55 % with 20 samples passes 11 missing values", {
    for (n in 4:5) {
        design <- k187Design(stats::setNames(c(n, n), c("A", "B")))
        missing <- lapply(0:n, \(k) if (k) c(paste0("A_", seq_len(k)), paste0("B_", seq_len(k))) else character())
        names(missing) <- paste0("k", 0:n)
        table <- k187Table(design, missing)
        expect_setequal(k187Helper(table, design, 50, 50)$Protein.Ids, c("k0", "k1", "k2", "P_floor"))
        expect_setequal(k187Helper(table, design, 60, 60)$Protein.Ids,
            c(if (n == 4L) c("k0", "k1", "k2") else c("k0", "k1", "k2", "k3"), "P_floor"))
    }
    design <- k187Design(c(A = 20L, B = 20L))
    table <- k187Table(design, list(P_full = character(), P_eleven = paste0("A_", 1:11)))
    expect_setequal(k187Helper(table, design, 55, 0)$Protein.Ids, c("P_full", "P_eleven", "P_floor"))
})

test_that("T7: the clean-up after RUV reads the stored cutoffs and applies the count rule", {
    args <- list(removeRowsWithMissingValuesPercent = list(groupwise_percentage_cutoff = 33.333,
        max_groups_percentage_cutoff = 0, proteins_intensity_cutoff_percentile = 0, ruv_grouping_variable = "group"))
    corrected <- k187Object(k187T1Table(), k187T1Design(), args)
    filtered <- k187Quiet(removeRowsWithMissingValuesPercent(theObject = corrected))
    expect_setequal(filtered@protein_quant_table$Protein.Ids, k187T1Kept)
})

test_that("T8 (i): a step replays with the count rule only when it lists the marker with source fixed", {
    expect_true(exists("workflowReplayRule", envir = asNamespace("MultiScholaR"), inherits = FALSE))
    if (exists("workflowReplayRule", envir = asNamespace("MultiScholaR"), inherits = FALSE)) {
        entry <- function(source) list(name = "missing_value_comparison", source = source,
            value = encodeWorkflowRecordValue("counts"), provenance = "fixed")
        for (operation in c("protein_intensity_filter", "normalization-ruv")) {
            step <- function(parameters) list(state = "s", operation = operation, parameters = parameters)
            expect_identical(workflowReplayRule(step(list(entry("fixed"))), "missing_value_comparison"), "counts",
                info = operation)
            expect_identical(workflowReplayRule(step(list()), "missing_value_comparison"), "rounded_percentage",
                info = operation)
            expect_identical(workflowReplayRule(step(list(entry("registry"))), "missing_value_comparison"),
                "rounded_percentage", info = operation)
        }
    }
})

test_that("T8 (ii): under the earlier rule the head gives the base's results", {
    withr::local_options(k187Legacy)
    expect_setequal(k187Filter(k187Object(k187T1Table(), k187T1Design()), 2, 2)$kept, "P_full")
    outcome <- k187GridOutcome()
    expect_identical(names(outcome)[!outcome], k187EarlierRuleFails())
})

test_that("T8 (iii): a step made in process lists the marker; a step without fixed constants lists none", {
    manager <- k187Manager(k187Object(k187T1Table(), k187T1Design()), withr::local_tempdir(), "protein_accession_cleaned")
    withr::defer(closeDIAStateAdapter(manager))
    workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list())))
    k187Quiet(runProteinIntensityFilterApplyStep(workflow, FALSE, 2, 2, 0))
    expect_identical(k187Marker(buildWorkflowAnalysisRecord(manager), "protein_intensity_filtered"),
        list(source = "fixed", value = "counts"))
    # The builder does not back-fill the marker into a step whose provenance has no fixed constants.
    other <- k187Manager(k187Object(k187T1Table(), k187T1Design()), withr::local_tempdir(), "protein_accession_cleaned")
    withr::defer(closeDIAStateAdapter(other))
    entry <- workflowOperationEntry("protein_intensity_filter")
    values <- list(strict_mode = FALSE, min_reps_per_group = 2L, min_groups = 2L, proteins_intensity_cutoff_percentile = 0)
    setDIAStepProvenance(other, "protein_intensity_filter", "whole", entry$parameters, values = values)
    other$saveState("protein_intensity_filtered", k187Object(k187T1Table(), k187T1Design()), values, "intensity filter")
    expect_null(k187Marker(buildWorkflowAnalysisRecord(other), "protein_intensity_filtered"))
})

# T8 (iv): an LFQ owner on the `small` recipe (2 groups of 6) with NA in 2 of 6 runs of group 1 for proteins 1 to 5 and
# in 3 of 6 for proteins 6 to 10; the filter at minimum 4 ((6, 4) is affected: 33.333 %), percentile 0 (the last run
# of the last protein is set to 10, below every other value, so it sets the threshold in group 2 only).
k187RuvParameters <- list(norm_method = "quantile", ruv_mode = "manual", auto_percentage_min = NULL,
    auto_percentage_max = NULL, separation_metric = NULL, k_penalty_weight = NULL, adaptive_k_penalty = NULL,
    ruv_percentage = 50, ruv_k = 1, grouping = "group")

k187ReplayOwner <- function(root, legacy_filter, normalisation = c("none", "shipped", "legacy")) {
    normalisation <- match.arg(normalisation)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, basename(root), "LFQ"), WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("LFQ")
    protein <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    table <- protein@protein_quant_table
    control <- as.character(protein@design_matrix$Run[protein@design_matrix$group == "control"])
    table[1:5, control[1:2]] <- NA_real_
    table[6:10, control[1:3]] <- NA_real_
    table[nrow(table), ncol(table)] <- 10
    protein@protein_quant_table <- table
    setDIAStepProvenance(manager, "design_builder", fixed = list(pivot = "protein x run", transform = "log2(x + 1)"))
    manager$saveState("protein_s4_initial", protein, list(), "design")
    entry <- workflowOperationEntry("protein_intensity_filter")
    parameters <- list(strict_mode = FALSE, min_reps_per_group = 4L, min_groups = 2L, proteins_intensity_cutoff_percentile = 0)
    filter <- function() withDIAStateScope(manager, function() {
        value <- computeDIAReferenceProteinIntensityFilter(manager, getDIAStateReference(manager), parameters)
        # A build before this change wrote no fixed constants for this step.
        setDIAStepProvenance(manager, "protein_intensity_filter", "whole", entry$parameters,
            if (!legacy_filter) entry$declaration$fixed, values = parameters)
        commitDIAReferenceState(manager, value$state_name, value$checkpoint, value$config, value$description)
    })
    normalise <- function() withDIAStateScope(manager, function() {
        run <- computeDIAReferenceNormalization(manager, getDIAStateReference(manager), k187RuvParameters)
        commit <- diaNormalizationCommit(manager, k187RuvParameters, run)
        commitDIAReferenceState(manager, commit$state_name, commit$checkpoint, commit$config, commit$description)
    })
    utils::capture.output(suppressWarnings(suppressMessages({
        if (legacy_filter) withr::with_options(k187Legacy, filter()) else filter()
        if (identical(normalisation, "shipped")) normalise()
        if (identical(normalisation, "legacy")) withr::with_options(k187Legacy, normalise())
    })))
    manager
}

# Write the owner's record, prepare a copy as the harness does and replay it, recording the rule in force at every
# call of the filter's helper and a separator after each replayed step.
k187Replay <- function(manager, root, copy, scratch) {
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
    calls <- character()
    replayed <- list()
    replayed_record <- NULL
    original <- get("removeRowsWithMissingValuesPercentHelper", envir = asNamespace("MultiScholaR"))
    testthat::local_mocked_bindings(removeRowsWithMissingValuesPercentHelper = function(...) {
        calls <<- c(calls, getOption("MultiScholaR.missing_filter_rule", "counts"))
        original(...)
    }, .package = "MultiScholaR")
    utils::capture.output(replay <- suppressMessages(suppressWarnings(replayWorkflowAnalysisRecord(record, prepared,
        character(), list(), scratch, on_step = function(step, replay_manager) {
            calls <<- c(calls, paste0("|", step$state))
            replayed[[step$state]] <<- withDIAStateScope(replay_manager, function() replay_manager$getState(step$state))
            replayed_record <<- buildWorkflowAnalysisRecord(replay_manager)
        }))))
    by_step <- list()
    current <- character()
    for (call in calls) {
        if (startsWith(call, "|")) {
            by_step[[substring(call, 2L)]] <- current
            current <- character()
        } else current <- c(current, call)
    }
    list(record = record, replayed_record = replayed_record, states = states, replay_states = replay$states,
        origin = origin, replayed = replayed, rules = by_step)
}

k187ExpectReplayIdentical <- function(result) {
    expect_identical(result$replay_states, result$states)
    for (state in names(result$origin)) expect_identical(result$replayed[[state]], result$origin[[state]], info = state)
}

test_that("T8 (iv): records made before and after the change replay, each step with its own rule", {
    proteins <- function(result) as.character(result$origin$protein_intensity_filtered@protein_quant_table$Protein.Group)
    first_five <- sprintf("PG_%06d", 1:5)
    filter_state <- "protein_intensity_filtered"
    # (a) before the change: no marker, the earlier rule.
    root <- withr::local_tempdir()
    owner <- k187ReplayOwner(root, legacy_filter = TRUE)
    withr::defer(closeDIAStateAdapter(owner))
    a <- k187Replay(owner, root, withr::local_tempdir(), withr::local_tempdir())
    expect_null(k187Marker(a$record, filter_state))
    expect_null(k187Marker(a$replayed_record, filter_state))
    expect_false(any(first_five %in% proteins(a)))
    k187ExpectReplayIdentical(a)
    # (b) as shipped: the marker, the count rule.
    root <- withr::local_tempdir()
    owner <- k187ReplayOwner(root, legacy_filter = FALSE)
    withr::defer(closeDIAStateAdapter(owner))
    b <- k187Replay(owner, root, withr::local_tempdir(), withr::local_tempdir())
    expect_identical(k187Marker(b$record, filter_state), list(source = "fixed", value = "counts"))
    expect_identical(k187Marker(b$replayed_record, filter_state), list(source = "fixed", value = "counts"))
    expect_true(all(first_five %in% proteins(b)))
    k187ExpectReplayIdentical(b)
    # (c) and (d): (a) continued through a normalisation with RUV applied, made as shipped (c) or before the change (d).
    root <- withr::local_tempdir()
    owner <- k187ReplayOwner(root, legacy_filter = TRUE, normalisation = "shipped")
    withr::defer(closeDIAStateAdapter(owner))
    c_ <- k187Replay(owner, root, withr::local_tempdir(), withr::local_tempdir())
    root <- withr::local_tempdir()
    owner <- k187ReplayOwner(root, legacy_filter = TRUE, normalisation = "legacy")
    withr::defer(closeDIAStateAdapter(owner))
    d <- k187Replay(owner, root, withr::local_tempdir(), withr::local_tempdir())
    norm_state <- "ruv_corrected"
    # (1)
    expect_identical(k187Marker(c_$record, norm_state), list(source = "fixed", value = "counts"))
    expect_identical(k187Marker(c_$replayed_record, norm_state), list(source = "fixed", value = "counts"))
    expect_null(k187Marker(c_$record, filter_state))
    expect_null(k187Marker(c_$replayed_record, filter_state))
    # (2)
    for (state in c(filter_state, norm_state)) {
        expect_null(k187Marker(d$record, state), label = paste("(d) record", state))
        expect_null(k187Marker(d$replayed_record, state), label = paste("(d) replayed record", state))
    }
    # (3)
    for (result in list(c_, d)) {
        expect_gte(length(result$rules[[filter_state]]), 1L)
        expect_true(all(result$rules[[filter_state]] == "rounded_percentage"))
    }
    # (4)
    expect_gte(length(c_$rules[[norm_state]]), 1L)
    expect_true(all(c_$rules[[norm_state]] == "counts"))
    expect_gte(length(d$rules[[norm_state]]), 1L)
    expect_true(all(d$rules[[norm_state]] == "rounded_percentage"))
    # (5)
    k187ExpectReplayIdentical(c_)
    k187ExpectReplayIdentical(d)
})

test_that("T9: the scripted entry keeps four of T1's five proteins", {
    manager <- k187Manager(k187Object(k187T1Table(), k187T1Design()), withr::local_tempdir(), "protein_accession_cleaned")
    withr::defer(closeDIAStateAdapter(manager))
    value <- k187Quiet(computeDIAReferenceProteinIntensityFilter(manager, getDIAStateReference(manager),
        list(strict_mode = FALSE, min_reps_per_group = 2L, min_groups = 2L, proteins_intensity_cutoff_percentile = 0)))
    expect_identical(value$protein_count, 4L)
})

test_that("T12: a cutoff that is not one finite number in [0, 100] is refused", {
    design <- k187T1Design()
    table <- k187T1Table()
    for (cutoff in list(NA_real_, -1, 101, c(10, 20), "50")) {
        outcome <- tryCatch(k187Helper(table, design, cutoff, 0), error = function(e) e)
        expect_s3_class(outcome, "multischolar_missing_filter_invalid")
    }
})
