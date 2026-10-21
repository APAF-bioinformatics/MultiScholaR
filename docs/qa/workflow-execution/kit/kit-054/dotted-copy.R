# KIT-054: a header-renamed copy of a real MaxQuant export for the legacy lane. The control build's reader
# (kit_036, before the KIT-054 reader repair) recognises only the dotted column names of the package's fixtures, so
# the legacy in-memory lane reads a copy whose header line carries the dotted names and whose rows are the original
# bytes. Data handling rule: nothing of the data is printed; the copy lives under the proof root (never the repository).
# Usage: Rscript docs/qa/workflow-execution/kit/kit-054/dotted-copy.R <source directory with manifest.json> <target directory>
`%||%` <- function(a, b) if (is.null(a)) b else a
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
source_dir <- normalizePath(args[[1L]], mustWork = TRUE); target <- args[[2L]]
manifest <- jsonlite::read_json(file.path(source_dir, "manifest.json"))
stopifnot(identical(manifest$format, "maxquant"))
dir.create(target, recursive = TRUE, showWarnings = FALSE)
# The same mapping the repaired reader applies (R/func_prot_import_readers.R, normaliseMaxQuantColumnNames).
normalise <- function(names) {
    fixed <- c("Protein IDs", "Majority protein IDs", "Protein names", "Gene names", "Potential contaminant",
        "Only identified by site", "Peptide counts (all)", "Peptide counts (razor+unique)", "Peptide counts (unique)",
        "Unique peptides", "MS/MS count")
    out <- names; at <- out %in% fixed; out[at] <- gsub(" ", ".", out[at], fixed = TRUE)
    for (prefix in c("LFQ intensity ", "Intensity ", "iBAQ ", "MS/MS count ", "Unique peptides ", "Peptides ")) {
        at <- startsWith(out, prefix)
        out[at] <- paste0(gsub(" ", ".", trimws(prefix), fixed = TRUE), ".", substring(out[at], nchar(prefix) + 1L))
    }
    out
}
report_in <- file.path(source_dir, manifest$report_file)
con <- file(report_in, "rb"); on.exit(close(con), add = TRUE)
header <- readLines(con, n = 1L, warn = FALSE)
columns <- strsplit(header, "\t", fixed = TRUE)[[1L]]
renamed <- normalise(columns)
report_out <- file.path(target, "proteinGroups.txt")
out <- file(report_out, "wb")
writeLines(paste(renamed, collapse = "\t"), out, sep = "\n")
# The rows are copied as bytes after the header line.
repeat { chunk <- readBin(con, "raw", 8L * 1024L^2); if (!length(chunk)) break; writeBin(chunk, out) }
close(out)
file.copy(file.path(source_dir, manifest$fasta_file), file.path(target, "proteins.fasta"), overwrite = TRUE)
file.copy(file.path(source_dir, manifest$design_file), file.path(target, "design.tsv"), overwrite = TRUE)
files <- c("proteinGroups.txt", "proteins.fasta", "design.tsv")
copy <- list(schema = "kit.corpus_manifest.v1", kind = "maxquant_real_dotted", synthetic = FALSE, format = "maxquant", sensitive = TRUE,
    note = "header-renamed copy of the user's MaxQuant export for the legacy lane (dotted column names; rows the original bytes); never copied into the repository",
    source = list(directory = source_dir, report_sha256 = manifest$files[[manifest$report_file]]$sha256),
    columns = length(columns), columns_renamed = sum(renamed != columns),
    samples = manifest$samples, rows = manifest$rows, design_rows = manifest$design_rows, design_groups = manifest$design_groups,
    files = lapply(stats::setNames(files, files), function(f) list(path = file.path(target, f), bytes = file.size(file.path(target, f)),
        sha256 = digest::digest(file = file.path(target, f), algo = "sha256"))),
    report_file = "proteinGroups.txt", fasta_file = "proteins.fasta", design_file = "design.tsv", written = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
jsonlite::write_json(copy, file.path(target, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
cat(sprintf("dotted copy written: %s | columns %d, renamed %d | bytes %d (source %d)\n", target, length(columns), copy$columns_renamed,
    file.size(report_out), file.size(report_in)))
