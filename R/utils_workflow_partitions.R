# Complete-sample partitioner (KIT-006) and feature-group partitioner (KIT-009). A step whose
# independence unit is the complete sample runs as batches of whole samples: the admission estimate
# chooses the fewest batches that fit the admitted worker memory, the batches are visited in the
# declared sample order, and the step's own kernels stage and consume each batch. A step whose unit is
# the feature group runs as batches of whole groups delivered by the group cursor, the kernel once per
# batch (KIT-009). A step whose unit is the row runs as pages of rows over a long or a wide table, the
# kernel once per page (KIT-010). Nothing here names a scientific table, column or operation.

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

#' Run a step over complete feature groups in batches, never splitting a group (KIT-009).
#'
#' The group cursor (`openWorkflowTableGroups`) delivers batches of whole groups in the declared
#' group order, up to the codec's delivery bound (10,000 rows of whole groups, or one larger group on
#' its own, never above 64 MiB). The admission estimate is the working set of the largest batch the
#' cursor can deliver, from the declared working bytes per hydrated row and the fixed overhead the
#' step names (engine limit, cohort read, transfer buffers); it is refused above the admitted cap and
#' recorded with the result. The kernel runs once per batch and `consume` receives each kernel value
#' in batch order; no batch is retained here.
#' @param scope Existing table read scope on the input.
#' @param columns Scientific columns the kernel reads, hydrated per batch.
#' @param group Grouping column of the independence unit.
#' @param kernel `function(part, index)`: `part` is one batch of whole groups with the attribute
#'   `workflow_source_order` (source ordinals of its rows).
#' @param consume `function(value, index)` run with each kernel value in batch order.
#' @param worker_cap Admitted worker memory in bytes.
#' @param working_bytes_per_row Declared working-set bytes per hydrated input row (the kernel's
#'   copies and intermediates), beside the hydrated row itself.
#' @param overhead_bytes Fixed bytes the step holds whatever the batch (engine limit, cohort read,
#'   transfer buffers, writers).
#' @param policy Name recorded with the admission.
#' @param group_order Passed to the cursor (`first` observation order or native `key` order).
#' @return `batches`, `batch_sizes` (rows per batch), `groups` (groups in the input), `admission`
#'   and the cursor's `dependencies`.
#' @noRd
partitionByFeatureGroup <- function(scope, columns, group, kernel, consume, worker_cap,
    working_bytes_per_row, overhead_bytes, policy, group_order = "first") {
    stopifnot(is.function(kernel), is.function(consume), is.numeric(working_bytes_per_row),
        length(working_bytes_per_row) == 1L, is.finite(working_bytes_per_row), working_bytes_per_row >= 0,
        is.numeric(overhead_bytes), length(overhead_bytes) == 1L, is.finite(overhead_bytes), overhead_bytes >= 0)
    cursor <- openWorkflowTableGroups(scope, columns, group, group_order)
    on.exit(cursor$close(), add = TRUE)
    row_bytes <- workflowProjectedRowBytes(scope$metadata, columns)
    # The cursor's delivery bound: 10,000 rows of whole groups, or one larger group alone, within 64 MiB.
    largest_group <- if (length(cursor$group_sizes)) max(cursor$group_sizes) else 0L
    batch_rows <- max(min(10000, floor(67108864 / max(1, row_bytes))), largest_group)
    estimate <- overhead_bytes + batch_rows * (row_bytes + working_bytes_per_row)
    if (!is.numeric(worker_cap) || length(worker_cap) != 1L || !is.finite(worker_cap) || worker_cap <= 0 ||
        estimate > worker_cap) cli::cli_abort(
        "feature-group working-memory estimate exceeds admitted worker capacity",
        class = "workflow_execution_unsupported")
    admission <- list(estimated_bytes = estimate, worker_cap = worker_cap, batch_rows = batch_rows,
        largest_group_rows = largest_group, row_bytes = row_bytes, working_bytes_per_row = working_bytes_per_row,
        overhead_bytes = overhead_bytes, policy = policy,
        limitation = "admission estimate; kernel cgroup guard remains authoritative")
    batch_sizes <- integer()
    index <- 0L
    repeat {
        part <- cursor$next_group()
        if (is.null(part)) break
        index <- index + 1L
        batch_sizes[[index]] <- nrow(part)
        value <- kernel(part, index)
        consume(value, index)
        part <- value <- NULL
    }
    cursor$close()
    list(batches = index, batch_sizes = batch_sizes, groups = length(cursor$group_sizes),
        admission = admission, dependencies = cursor$dependencies)
}

#' Run a step over pages of rows of a long or wide table, the kernel once per page (KIT-010).
#'
#' Pages come from the row streamer (`streamWorkflowTableRows`: 10,000 rows per page, every requested
#' column, a wide table read in 64-column blocks), in source order. The admission estimate is one page's
#' working set from the declared working bytes per row plus the step's fixed overhead, refused above the
#' cap and recorded. The kernel runs once per page and `consume` receives each kernel value in page
#' order; no page is retained here.
#' @param scope Existing table read scope on the input.
#' @param columns Columns the kernel reads (all of a wide assay's sample columns, when declared wide).
#' @param kernel `function(page, index, offset)`: `page` is one batch of rows, `offset` its first
#'   row's zero-based position in the source.
#' @param consume `function(value, index)` run with each kernel value in page order.
#' @param worker_cap Admitted worker memory in bytes.
#' @param working_bytes_per_row Declared working-set bytes per row (the kernel's copies and
#'   intermediates), beside the row itself.
#' @param overhead_bytes Fixed bytes the step holds whatever the page (transfer buffers, writers).
#' @param policy Name recorded with the admission.
#' @return `batches`, `batch_sizes` (rows per page), `rows` and `admission`.
#' @noRd
partitionByRows <- function(scope, columns, kernel, consume, worker_cap, working_bytes_per_row,
    overhead_bytes, policy) {
    stopifnot(is.function(kernel), is.function(consume), is.numeric(working_bytes_per_row),
        length(working_bytes_per_row) == 1L, is.finite(working_bytes_per_row), working_bytes_per_row >= 0,
        is.numeric(overhead_bytes), length(overhead_bytes) == 1L, is.finite(overhead_bytes), overhead_bytes >= 0)
    checkWorkflowTableReadScope(scope)
    if (!is.character(columns) || !length(columns) || anyDuplicated(columns) || !all(columns %in% scope$metadata$columns))
        cli::cli_abort("invalid row partition projection", class = "saved_entity_limit")
    row_bytes <- workflowProjectedRowBytes(scope$metadata, columns)
    page_rows <- min(10000L, max(1L, scope$metadata$rows))
    estimate <- overhead_bytes + page_rows * (row_bytes + working_bytes_per_row)
    if (!is.numeric(worker_cap) || length(worker_cap) != 1L || !is.finite(worker_cap) || worker_cap <= 0 ||
        estimate > worker_cap) cli::cli_abort(
        "row-page working-memory estimate exceeds admitted worker capacity",
        class = "workflow_execution_unsupported")
    admission <- list(estimated_bytes = estimate, worker_cap = worker_cap, page_rows = page_rows,
        row_bytes = row_bytes, working_bytes_per_row = working_bytes_per_row, overhead_bytes = overhead_bytes,
        columns = length(columns), policy = policy,
        limitation = "admission estimate; kernel cgroup guard remains authoritative")
    batch_sizes <- integer()
    index <- 0L
    streamWorkflowTableRows(scope, columns, function(page, offset, row_names) {
        if (!nrow(page) && index > 0L) return(invisible(NULL))
        index <<- index + 1L
        batch_sizes[[index]] <<- nrow(page)
        value <- kernel(page, index, offset)
        consume(value, index)
        page <- value <- NULL
        invisible(NULL)
    })
    list(batches = index, batch_sizes = batch_sizes, rows = sum(batch_sizes), admission = admission)
}
