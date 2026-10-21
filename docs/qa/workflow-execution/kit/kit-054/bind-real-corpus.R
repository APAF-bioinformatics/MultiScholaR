# KIT-054: write the manifest of a real (non-synthetic) protein corpus beside its files, by hash and count only.
# Data handling rule (ticket): nothing of the data is printed, logged or committed; the manifest carries sha256s,
# byte sizes, the format and the counts the reader reports (rows, protein groups, runs, design rows and groups).
# Usage: Rscript docs/qa/workflow-execution/kit/kit-054/bind-real-corpus.R <directory> <format: maxquant|pd_tmt> <report relative path> [fasta] [design]
suppressMessages(library(MultiScholaR))
`%||%` <- function(a, b) if (is.null(a)) b else a
args <- commandArgs(TRUE)
stopifnot(length(args) >= 3L)
directory <- normalizePath(args[[1L]], mustWork = TRUE)
format <- args[[2L]]
report <- args[[3L]]
fasta <- if (length(args) > 3L) args[[4L]] else list.files(directory, pattern = "[.]fasta$", full.names = FALSE)[[1L]]
design <- if (length(args) > 4L) args[[5L]] else "design.tsv"
files <- c(report, fasta, design)
for (f in files) stopifnot(file.exists(file.path(directory, f)))
imported <- switch(format, maxquant = importMaxQuantData(file.path(directory, report), use_lfq = TRUE, filter_contaminants = TRUE),
    pd_tmt = importProteomeDiscovererTMTData(file.path(directory, report)), stop("format: maxquant | pd_tmt"))
id_col <- imported$column_mapping$protein_col %||% names(imported$data)[[1L]]
run_col <- imported$column_mapping$run_col %||% "Run"
design_table <- utils::read.delim(file.path(directory, design), stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("sample", "group") %in% names(design_table)))
runs <- unique(as.character(imported$data[[run_col]]))
matched <- sum(design_table$sample %in% janitor::make_clean_names(runs))
manifest <- list(schema = "kit.corpus_manifest.v1", kind = paste0(format, "_real"), synthetic = FALSE, format = format,
    sensitive = TRUE, note = "user-supplied dataset; bound by hash; never copied into the repository; counts only",
    samples = length(runs), rows = length(unique(imported$data[[id_col]])), long_rows = nrow(imported$data),
    design_rows = nrow(design_table), design_groups = length(unique(design_table$group)), design_samples_matched = matched,
    files = lapply(stats::setNames(files, files), function(f) list(path = file.path(directory, f),
        bytes = file.size(file.path(directory, f)), sha256 = digest::digest(file = file.path(directory, f), algo = "sha256"))),
    report_file = report, fasta_file = fasta, design_file = design, written = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
jsonlite::write_json(manifest, file.path(directory, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
cat(sprintf("manifest written: %s | format %s | protein groups %d | runs %d | design rows %d in %d groups, %d matched\n",
    file.path(directory, "manifest.json"), format, manifest$rows, manifest$samples, manifest$design_rows, manifest$design_groups, matched))
