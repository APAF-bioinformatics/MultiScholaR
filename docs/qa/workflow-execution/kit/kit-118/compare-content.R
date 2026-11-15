# KIT-118 M1, evidence item 2: the DIA pair by content. For each state, the object's main table in both stores is compared
# cell by cell: rows, numeric cells that differ, cells that became missing, cells that became observed, and the largest
# absolute and relative change among cells observed in both. Each table is read whole through the harness's pages and ordered by its non-numeric columns; a table whose
# non-numeric columns then differ between the two stores is marked not aligned and not compared.
# Usage (from the repository root): Rscript compare-content.R <library> <before member dir> <after member dir> <out tsv>
args <- commandArgs(TRUE)
stopifnot(length(args) == 4L)
.libPaths(c(args[[1L]], .libPaths()))
ns <- loadNamespace("MultiScholaR")
api <- function(name) get(name, ns)
env <- new.env(parent = globalenv())
sys.source("tools/profiling/workflow_kit.R", env)
open_member <- function(member) {
    dir <- readRDS(file.path(member, "journey-result.rds"))$store
    store <- env$kitOpenStore(dir, api)
    list(store = store, snapshot = api("hydrateSavedWorkflowCheckpoint")(store, readRDS(file.path(dir, "dia-state-head.rds"))))
}
component <- function(member, state) {
    record <- api("hydrateSavedWorkflowCheckpoint")(member$store, member$snapshot$states[[state]])
    for (name in c("peptide_data", "protein_quant_table")) if (!is.null(record$data$components[[name]])) return(record$data$components[[name]])
    NULL
}
a <- open_member(args[[2L]]); b <- open_member(args[[3L]])
rows <- list()
for (state in names(a$snapshot$states)) {
    ra <- component(a, state); rb <- component(b, state)
    if (is.null(ra) || is.null(rb)) next
    sa <- api("newWorkflowTableReadScope")(a$store, ra, reader = "fastduckdb")
    sb <- api("newWorkflowTableReadScope")(b$store, rb, reader = "fastduckdb")
    row <- data.frame(state = state, rows_before = sa$metadata$rows, rows_after = sb$metadata$rows, numeric_cells = NA_real_, differ = NA_real_,
        became_missing = NA_real_, became_observed = NA_real_, max_abs_change = NA_real_, max_rel_change = NA_real_, aligned = FALSE)
    if (row$rows_before == row$rows_after && identical(sa$metadata$columns, sb$metadata$columns) && row$rows_before > 0) {
        columns <- sa$metadata$columns
        whole <- function(store, reference, scope) {
            pages <- lapply(seq.int(0L, row$rows_before - 1L, by = 10000L), function(offset) as.data.frame(
                api("readSavedWorkflowTable")(store, reference, columns, offset = offset, limit = min(10000L, row$rows_before - offset),
                    max_columns = max(64L, length(columns)), scope = scope)$data))
            table <- as.data.frame(dplyr::bind_rows(pages))
            keys <- names(table)[!vapply(table, is.double, logical(1))]
            if (length(keys)) table <- table[do.call(order, c(unname(lapply(table[keys], as.character)), list(method = "radix"))), , drop = FALSE]
            table
        }
        x <- whole(a$store, ra, sa); y <- whole(b$store, rb, sb)
        numeric <- names(x)[vapply(x, is.double, logical(1))]
        keys <- setdiff(names(x), numeric)
        if (!length(keys) || isTRUE(all.equal(x[keys], y[keys], check.attributes = FALSE))) {
            u <- as.matrix(x[numeric]); v <- as.matrix(y[numeric])
            changed <- !is.na(u) & !is.na(v) & u != v
            row$aligned <- TRUE
            row$numeric_cells <- length(u)
            row$became_missing <- sum(!is.na(u) & is.na(v))
            row$became_observed <- sum(is.na(u) & !is.na(v))
            row$differ <- sum(changed) + row$became_missing + row$became_observed
            row$max_abs_change <- if (any(changed)) max(abs(v[changed] - u[changed])) else 0
            row$max_rel_change <- if (any(changed)) max(abs(v[changed] - u[changed]) / pmax(abs(u[changed]), .Machine$double.eps)) else 0
            row$median_abs_change <- if (any(changed)) stats::median(abs(v[changed] - u[changed])) else 0
        } else row$median_abs_change <- NA_real_
        rm(x, y); gc(verbose = FALSE)
    } else row$median_abs_change <- NA_real_
    api("releaseWorkflowTableReadScope")(sa); api("releaseWorkflowTableReadScope")(sb)
    rows[[length(rows) + 1L]] <- row
}
out <- do.call(rbind, rows)
utils::write.table(out, args[[4L]], sep = "\t", quote = FALSE, row.names = FALSE)
print(out, row.names = FALSE)
