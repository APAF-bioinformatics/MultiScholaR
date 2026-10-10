# Complete-sample partitioner (KIT-006). A step whose independence unit is the complete sample runs
# as batches of whole samples: the admission estimate chooses the fewest batches that fit the admitted
# worker memory, the batches are visited in the declared sample order, and the step's own kernels
# stage and consume each batch. Nothing here names a scientific table, column or operation.

#' Upper bound of a per-batch quantity when `n` ordered groups are split into `batches` batches:
#' the sum of the largest `ceiling(n / batches)` values, so any one batch is bounded whatever the
#' assignment.
#' @param values Per-group quantities (one per group).
#' @param batches Number of batches.
#' @return The bound, or 0 for no groups.
#' @noRd
completeGroupBatchBound <- function(values, batches) {
    values <- as.numeric(values)
    if (!length(values)) return(0)
    n <- ceiling(length(values) / batches)
    sum(utils::head(sort(values, decreasing = TRUE), n))
}

#' Choose the fewest complete-group batches whose estimate fits the admitted worker memory.
#' @param counts Candidate batch counts in order of preference (the first that fits wins); the
#'   caller supplies its minimum (a profile's fixed count, or the fixed count alone when batching
#'   is not adaptive).
#' @param estimate_for `function(batches)` returning the estimated working set in bytes.
#' @param worker_cap Admitted worker memory in bytes.
#' @param policy Name recorded with the admission.
#' @return `batches`, `estimated_bytes`, `estimates` (a list indexed by batch count, NULL where a
#'   count was not a candidate), `worker_cap`, `policy` and the shared limitation note.
#' @noRd
chooseCompleteGroupBatches <- function(counts, estimate_for, worker_cap, policy) {
    counts <- as.integer(counts)
    stopifnot(length(counts) >= 1L, all(counts >= 1L), !anyDuplicated(counts), is.function(estimate_for))
    estimates <- vapply(counts, estimate_for, numeric(1))
    if (!is.numeric(worker_cap) || length(worker_cap) != 1L || !is.finite(worker_cap) ||
        worker_cap <= 0 || !any(estimates <= worker_cap)) cli::cli_abort(
        "working-memory estimate exceeds admitted worker capacity",
        class = "workflow_execution_unsupported")
    chosen <- which(estimates <= worker_cap)[[1L]]
    by_count <- vector("list", max(counts))
    by_count[counts] <- as.list(unname(estimates))
    list(batches = counts[[chosen]], estimated_bytes = estimates[[chosen]], estimates = by_count,
        worker_cap = worker_cap, policy = policy,
        limitation = "admission estimate; kernel cgroup guard remains authoritative")
}

#' Run a step over complete groups in batches, never splitting a group.
#' @param groups Ordered unique group identifiers (the step's output order).
#' @param batches Number of batches (from `chooseCompleteGroupBatches` or a fixed profile).
#' @param stage `function(selection, index)` run for each batch in order; its value is the batch's
#'   staged result (a file, a count table, ...).
#' @param between Optional `function()` run once after every batch is staged and before any is
#'   consumed (a step that must release its engine before reading its stages).
#' @param consume Optional `function(staged, index)` run for each staged result in batch order.
#' @return `batch_sizes` (groups per batch), `staged` (the stage values) and `consumed` (the consume
#'   values, or NULL).
#' @noRd
partitionBySamples <- function(groups, batches, stage, between = NULL, consume = NULL) {
    stopifnot(is.function(stage), is.null(between) || is.function(between), is.null(consume) || is.function(consume))
    selections <- partitionWorkflowCompleteGroups(groups, batches)
    staged <- vector("list", length(selections))
    for (index in seq_along(selections)) staged[[index]] <- stage(selections[[index]], index)
    if (!is.null(between)) between()
    consumed <- if (is.null(consume)) NULL else
        lapply(seq_along(staged), function(index) consume(staged[[index]], index))
    list(batch_sizes = as.integer(lengths(selections)), staged = staged, consumed = consumed)
}
