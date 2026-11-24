# KIT-193 target 1's paired test: shared helpers of test-kit193-pairs-filter.R and test-kit193-pairs-filter-table.R
# (the shipped cutoffs and helper on synthetic designs; the Lean oracle table). See those files' headers.

kit193FilterTable <- function() {
    utils::read.delim(testthat::test_path("..", "testdata", "kit193", "filter-thresholds.tsv"), sep = "\t",
        stringsAsFactors = FALSE)
}

kit193FilterRow <- function(table, nstar, m, n) {
    row <- table[table$nstar == nstar & table$m == m & table$n == n, , drop = FALSE]
    if (nrow(row) != 1L) stop(sprintf("no unique oracle row for n* = %d, m = %d, n = %d", nstar, m, n))
    row
}

# I10 gap 1: the earlier comparison applies only when the marker's option holds its earlier token; the model describes
# the count rule, so the option is unset here.
kit193NoLegacy <- function(env = parent.frame()) {
    option <- workflowReplayRuleMarkers$missing_value_comparison$option
    withr::local_options(stats::setNames(list(NULL), option), .local_envir = env)
    invisible(option)
}

if (!methods::isClass("Kit193Design")) {
    methods::setClass("Kit193Design", slots = c(design_matrix = "data.frame", args = "list"))
}

kit193Design <- function(sizes) {
    groups <- rep(paste0("G", seq_along(sizes)), sizes)
    data.frame(Run = paste0(groups, "_", sequence(sizes)), group = groups, stringsAsFactors = FALSE)
}

#' The shipped cutoffs for a design: `updateMissingValueParameters()` on a minimal S4 object, grouping column `group`
kit193Cutoffs <- function(group_sizes, min_reps, min_groups) {
    object <- methods::new("Kit193Design", design_matrix = kit193Design(group_sizes), args = list())
    env <- new.env()
    assign("config_list", list(), envir = env)
    object <- suppressMessages(updateMissingValueParameters(object, min_reps_per_group = min_reps,
        min_groups = min_groups, grouping_variable = "group", env = env))
    object@args$removeRowsWithMissingValuesPercent
}

#' A table: `missing` maps a row id to the runs set to NA; every other value is 20. `P_floor`'s single value of 10 in
#' the first run sets the intensity threshold at percentile 0 (present means strictly above it).
kit193Table <- function(design, missing) {
    runs <- design$Run
    ids <- c(names(missing), "P_floor")
    values <- matrix(20, nrow = length(ids), ncol = length(runs), dimnames = list(ids, runs))
    for (id in names(missing)) values[id, missing[[id]]] <- NA_real_
    values["P_floor", 1L] <- 10
    data.frame(Protein.Ids = ids, values, check.names = FALSE, stringsAsFactors = FALSE)
}

#' The rows the shipped helper keeps (`P_floor` aside), grouping column `group` (I10 gap 4)
kit193Kept <- function(table, design, cutoffs) {
    kept <- NULL
    utils::capture.output(kept <- suppressMessages(removeRowsWithMissingValuesPercentHelper(table,
        cols = "Protein.Ids", design_matrix = design, sample_id = Run, row_id = Protein.Ids, grouping_variable = group,
        groupwise_percentage_cutoff = cutoffs$groupwise_percentage_cutoff,
        max_groups_percentage_cutoff = cutoffs$max_groups_percentage_cutoff, proteins_intensity_cutoff_percentile = 0)))
    setdiff(as.character(kept$Protein.Ids), "P_floor")
}

#' Rows `k0` to `kn`: row `kj` has the first j runs of group `group` missing
kit193Missing <- function(group, n) {
    stats::setNames(lapply(0:n, function(j) if (j) paste0(group, "_", seq_len(j)) else character()), paste0("k", 0:n))
}
