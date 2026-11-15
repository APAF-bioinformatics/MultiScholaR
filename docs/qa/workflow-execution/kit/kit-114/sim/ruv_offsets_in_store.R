# Row 14: per-protein offsets between the state before normalisation and ruv_corrected, read from a journey member store.
# Usage (repository root): Rscript ruv_offsets_in_store.R <library> <member dir> ...
args <- commandArgs(TRUE)
.libPaths(c(args[[1L]], .libPaths()))
ns <- loadNamespace("MultiScholaR"); api <- function(name) get(name, ns)
env <- new.env(parent = globalenv()); sys.source("tools/profiling/workflow_kit.R", env)
for (member in args[-1L]) {
    dir <- readRDS(file.path(member, "journey-result.rds"))$store
    store <- env$kitOpenStore(dir, api)
    snapshot <- api("hydrateSavedWorkflowCheckpoint")(store, readRDS(file.path(dir, "dia-state-head.rds")))
    read <- function(state) {
        record <- api("hydrateSavedWorkflowCheckpoint")(store, snapshot$states[[state]])
        ref <- record$data$components$protein_quant_table
        if (is.null(ref)) return(NULL)
        scope <- api("newWorkflowTableReadScope")(store, ref, reader = "fastduckdb"); on.exit(api("releaseWorkflowTableReadScope")(scope))
        n <- scope$metadata$rows; columns <- scope$metadata$columns
        table <- as.data.frame(dplyr::bind_rows(lapply(seq.int(0L, n - 1L, by = 10000L), function(offset) as.data.frame(
            api("readSavedWorkflowTable")(store, ref, columns, offset = offset, limit = min(10000L, n - offset),
                max_columns = max(64L, length(columns)), scope = scope)$data))))
        key <- names(table)[!vapply(table, is.double, logical(1))][[1L]]
        m <- as.matrix(table[vapply(table, is.double, logical(1))]); rownames(m) <- as.character(table[[key]]); m
    }
    cat(basename(dirname(member)), "states:", paste(names(snapshot$states), collapse = " "), "\n")
    if (!all(c("protein_replicate_filtered", "ruv_corrected") %in% names(snapshot$states))) next
    a <- read("protein_replicate_filtered"); b <- read("ruv_corrected")
    common <- intersect(rownames(a), rownames(b)); samples <- intersect(colnames(a), colnames(b))
    a <- a[common, samples]; b <- b[common, samples]
    offset <- rowMeans(b, na.rm = TRUE) - rowMeans(a, na.rm = TRUE)
    cat(sprintf("  proteins %d samples %d | input range %.1f to %.1f | corrected range %.1f to %.1f\n", length(common), length(samples),
        min(a, na.rm = TRUE), max(a, na.rm = TRUE), min(b, na.rm = TRUE), max(b, na.rm = TRUE)))
    cat(sprintf("  per-protein offset: median |%.2f|, 90%% |%.2f|, max |%.2f| log2; share of proteins shifted by more than 1 log2: %.0f%%\n",
        median(abs(offset)), quantile(abs(offset), 0.9), max(abs(offset)), 100 * mean(abs(offset) > 1)))
    cat(sprintf("  correlation of protein mean before and after: %.3f | within-protein SD ratio (median): %.2f\n",
        cor(rowMeans(a, na.rm = TRUE), rowMeans(b, na.rm = TRUE)), median(apply(b, 1, sd, na.rm = TRUE) / apply(a, 1, sd, na.rm = TRUE), na.rm = TRUE)))
}
