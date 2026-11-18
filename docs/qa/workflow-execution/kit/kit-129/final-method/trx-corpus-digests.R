# KIT-129 fifth pass item 2 (a): which files of the regenerated KIT-040 transcriptomics corpora differ off the laptop.
# Regenerates trx_edge_12, trx_edge_16 and trx_screen as the tests' helper does (tools/profiling/kit_omics_corpus.R,
# sourced, not changed) and writes one row per file: bytes and sha256 against kit-040/bindings.json; for truth.json also
# the sha256 of each line, to be located against the laptop's file by reading (digests only).
#
# Usage (from the package root): Rscript trx-corpus-digests.R <out dir>
# Writes <out dir>/trx-corpus-digests.tsv, truth-lines-<corpus>.tsv, trx-corpus-meta.txt.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L)
out_dir <- arguments[[1L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
bindings <- jsonlite::read_json("docs/qa/workflow-execution/kit/kit-040/bindings.json")$corpora
per_group <- c(trx_edge_12 = 3L, trx_edge_16 = 4L, trx_screen = 6L)
generator <- new.env(parent = globalenv())
sys.source("tools/profiling/kit_omics_corpus.R", generator)

rows <- list()
for (corpus in names(per_group)) {
    dir <- file.path(tempdir(), corpus)
    generator$writeKitTrxCorpus(dir, generator$kitOmicsRecipe("trx_screen", per_group[[corpus]]))
    bound <- bindings[[corpus]]
    for (file in c(names(bound$files), "manifest.json")) {
        path <- file.path(dir, file)
        sha <- digest::digest(file = path, algo = "sha256")
        bound_sha <- if (file == "manifest.json") bound$manifest_sha256 else bound$files[[file]]$sha256
        bound_bytes <- if (file == "manifest.json") NA_integer_ else as.integer(bound$files[[file]]$bytes)
        rows[[length(rows) + 1L]] <- data.frame(corpus = corpus, file = file, bytes = file.size(path), bound_bytes = bound_bytes,
            sha256 = sha, bound_sha256 = bound_sha, same = identical(sha, bound_sha), stringsAsFactors = FALSE)
    }
    # the line digests of truth.json, to locate a differing line against the laptop's file by reading
    lines <- readLines(file.path(dir, "truth.json"))
    utils::write.table(data.frame(line = seq_along(lines), sha256 = vapply(lines, digest::digest, "", algo = "sha256",
        serialize = FALSE, USE.NAMES = FALSE)), file.path(out_dir, sprintf("truth-lines-%s.tsv", corpus)), sep = "\t",
        quote = FALSE, row.names = FALSE)
}
table <- do.call(rbind, rows)
utils::write.table(table, file.path(out_dir, "trx-corpus-digests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(R.version.string, paste("jsonlite", utils::packageVersion("jsonlite")), paste("digest", utils::packageVersion("digest")),
    paste(names(extSoftVersion()), extSoftVersion(), collapse = "; "), paste("BLAS", extSoftVersion()[["BLAS"]]),
    paste("La", La_library())), file.path(out_dir, "trx-corpus-meta.txt"))
print(table[, c("corpus", "file", "bytes", "bound_bytes", "same")], row.names = FALSE)
