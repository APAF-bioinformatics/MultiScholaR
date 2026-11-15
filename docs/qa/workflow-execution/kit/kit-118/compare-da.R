# KIT-118 M1: the DIA pair's differential-abundance tables compared by protein (the RUV-corrected intensities carry an
# arbitrary per-protein offset, so they are not comparable cell by cell; the contrasts are).
# Usage (from the repository root): Rscript compare-da.R <library> <before member dir> <after member dir>
args <- commandArgs(TRUE)
.libPaths(c(args[[1L]], .libPaths()))
ns <- loadNamespace("MultiScholaR"); api <- function(name) get(name, ns)
env <- new.env(parent = globalenv()); sys.source("tools/profiling/workflow_kit.R", env)
da <- function(member) {
    dir <- readRDS(file.path(member, "journey-result.rds"))$store
    store <- env$kitOpenStore(dir, api)
    snapshot <- api("hydrateSavedWorkflowCheckpoint")(store, readRDS(file.path(dir, "dia-state-head.rds")))
    value <- api("hydrateSavedWorkflowEntity")(store, snapshot$da_results)
    lapply(value$results, function(result) {
        tables <- result$tables["contrasts_results_table"]
        lapply(tables, function(reference) {
            scope <- api("newWorkflowTableReadScope")(store, reference, reader = "fastduckdb")
            rows <- scope$metadata$rows; columns <- scope$metadata$columns
            as.data.frame(dplyr::bind_rows(lapply(seq.int(0L, rows - 1L, by = 10000L), function(offset) as.data.frame(
                api("readSavedWorkflowTable")(store, reference, columns, offset = offset, limit = min(10000L, rows - offset),
                    max_columns = max(64L, length(columns)), scope = scope)$data))))
        })
    })
}
a <- da(args[[2L]]); b <- da(args[[3L]])
for (contrast in names(a)) for (name in names(a[[contrast]])) {
    x <- a[[contrast]][[name]]; y <- b[[contrast]][[name]]
    cat(contrast, name, ": rows", nrow(x), "->", nrow(y), "; columns", paste(names(x), collapse = ", "), "\n")
    id <- names(x)[[1L]]
    common <- intersect(x[[id]], y[[id]])
    cat("  proteins only before", length(setdiff(x[[id]], y[[id]])), ", only after", length(setdiff(y[[id]], x[[id]])), ", common", length(common), "\n")
    x <- x[match(common, x[[id]]), ]; y <- y[match(common, y[[id]]), ]
    for (column in intersect(c("log2FC", "raw_pvalue", "fdr_qvalue", "fdr_value_bh_adjustment"), names(x))) {
        d <- abs(y[[column]] - x[[column]])
        cat(sprintf("  %-22s changed %d of %d; |change| median %.3g, 99%% %.3g, max %.3g\n", column, sum(d > 0, na.rm = TRUE), length(d),
            stats::median(d, na.rm = TRUE), stats::quantile(d, 0.99, na.rm = TRUE), max(d, na.rm = TRUE)))
    }
    for (column in grep("qvalue|fdr|adj", names(x), value = TRUE, ignore.case = TRUE)) if (is.numeric(x[[column]])) {
        cat(sprintf("  significant at %s < 0.05: before %d, after %d, both %d\n", column, sum(x[[column]] < 0.05, na.rm = TRUE),
            sum(y[[column]] < 0.05, na.rm = TRUE), sum(x[[column]] < 0.05 & y[[column]] < 0.05, na.rm = TRUE)))
    }
}
