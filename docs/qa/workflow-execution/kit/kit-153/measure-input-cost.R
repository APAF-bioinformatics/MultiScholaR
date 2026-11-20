# KIT-153: the cost of carrying each data set's RUV-corrected matrix beside its normalised one. Two omics, one data set
# each, 10,000 features by 60 runs (random log2 values), written by the contract's writer without and with a corrected
# matrix, then imported into an integration store (the import keeps both matrices of a data set as tables). Reported:
# the input folders' bytes, the store's bytes after the import and the import's elapsed seconds, three repeats each, in
# alternating order.
#
# Usage (from the package root, MultiScholaR installed): Rscript measure-input-cost.R <out dir>
# Writes <out dir>/input-cost.tsv and <out dir>/input-cost-summary.tsv.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L)
out_dir <- arguments[[1L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ns <- asNamespace("MultiScholaR")
features <- 10000L
runs <- 60L

writeCostInput <- function(root, omic, view, corrected) {
    biological <- sprintf("%s_%d", rep(c("A", "B"), each = runs / 2L), rep(seq_len(runs / 2L), 2L))
    run_ids <- paste0(omic, "_", biological)
    values <- withr::with_seed(153L, matrix(stats::rnorm(features * runs, 20, 2), features, runs,
        dimnames = list(sprintf("F%05d", seq_len(features)), run_ids)))
    table <- data.frame(feature_id = rownames(values), values, check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL)
    samples <- data.frame(Run = run_ids, group = substr(biological, 1L, 1L), tech_rep_group = biological, stringsAsFactors = FALSE)
    annotation <- data.frame(feature_id = rownames(values), display_name = rownames(values), stringsAsFactors = FALSE)
    manifest <- list(omic_type = omic, workflow = "TEST", matrix_role = "normalised", scale = "log2",
        correction = if (corrected) list(id = "C1", method = "RUVIII_C_Varying", replicate_structure = "group"),
        source_state = list(name = "kit153-cost"), sample_id_column = "Run", group_column = "group",
        technical_replicate_column = "tech_rep_group", biological_sample_column = "tech_rep_group")
    shifted <- table
    shifted[-1L] <- shifted[-1L] - 0.1
    dir <- file.path(root, "integration", omic)
    ns$writeIntegrationInput(dir, manifest, samples, stats::setNames(list(table), view), stats::setNames(list(annotation), view),
        corrected = if (corrected) stats::setNames(list(shifted), view) else list())
    dir
}

folderBytes <- function(dir) sum(file.size(list.files(dir, recursive = TRUE, all.files = TRUE, full.names = TRUE)))

rows <- list()
for (repeat_index in 1:3) {
    for (corrected in if (repeat_index %% 2L) c(FALSE, TRUE) else c(TRUE, FALSE)) {
        root <- tempfile("kit153-cost-")
        dir.create(root)
        dirs <- c(writeCostInput(root, "proteomics", "proteomics", corrected),
            writeCostInput(root, "lipidomics", "lipidomics_LCMS_Pos", corrected))
        manager <- ns$newIntegrationStateOwner(file.path(root, "project"))
        seconds <- system.time(ns$importIntegrationInputs(manager, dirs))[["elapsed"]]
        ns$closeWorkflowStateOwner(manager)
        rows[[length(rows) + 1L]] <- data.frame(repeat_index = repeat_index, corrected = corrected,
            input_bytes = sum(vapply(dirs, folderBytes, numeric(1))), store_bytes = folderBytes(file.path(root, "project")),
            import_seconds = seconds)
        unlink(root, recursive = TRUE)
    }
}
table <- do.call(rbind, rows)
utils::write.table(table, file.path(out_dir, "input-cost.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
summary <- do.call(rbind, lapply(split(table, table$corrected), function(part) data.frame(corrected = part$corrected[[1L]],
    input_bytes = stats::median(part$input_bytes), store_bytes = stats::median(part$store_bytes),
    import_seconds_median = stats::median(part$import_seconds))))
utils::write.table(summary, file.path(out_dir, "input-cost-summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summary, row.names = FALSE)
