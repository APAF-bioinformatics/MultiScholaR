# KIT-209 (branch run-kit209-memory): generate one synthetic corpus on the runner with the committed seeded generators
# of the candidate's checkout (tools/profiling/kit_omics_corpus.R, kit_protein_corpus.R, workflow_perf107_synthetic.R),
# then compare each file's SHA-256 with the bindings' value (docs/qa/workflow-execution/kit/bindings.json). Real data
# (`kc_*`) is refused by name. Writes <out tsv>: corpus, file, bytes, sha256, bound sha256, match.
#
# Usage: Rscript --vanilla k209-corpus.R <corpus key> <directory> <checkout> <out tsv>
args <- commandArgs(trailingOnly = TRUE)
key <- args[[1L]]
directory <- args[[2L]]
repo <- normalizePath(args[[3L]], mustWork = TRUE)
out <- args[[4L]]
if (grepl("^kc", key)) stop("real data is refused by name: ", key)
tools <- file.path(repo, "tools", "profiling")
bindings <- jsonlite::read_json(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"), simplifyVector = TRUE)
bound <- bindings$corpora[[key]]
if (is.null(bound)) stop("corpus ", key, " is not bound")
if (!isTRUE(bound$synthetic)) stop("corpus ", key, " is not synthetic")
generators <- new.env(parent = globalenv())
`%||%` <- function(a, b) if (is.null(a)) b else a
omics <- function() sys.source(file.path(tools, "kit_omics_corpus.R"), generators)
protein <- function() sys.source(file.path(tools, "kit_protein_corpus.R"), generators)
recipe <- bound$recipe
manifest <- switch(key,
    lipid_screen = { omics(); generators$writeKitLipidCorpus(directory, generators$kitOmicsRecipe("lipid_screen")) },
    metab_screen = { omics(); generators$writeKitMetabCorpus(directory, generators$kitOmicsRecipe("metab_screen")) },
    trx_screen = { omics(); generators$writeKitTrxCorpus(directory, generators$kitOmicsRecipe("trx_screen")) },
    lipid_large_missing = { omics(); generators$writeKitLipidCorpus(directory, generators$kitOmicsRecipe("lipid_screen",
        recipe$biological_per_group, recipe$features_per_assay, recipe$standards_per_assay, recipe$missing_features)) },
    metab_large_missing = { omics(); generators$writeKitMetabCorpus(directory, generators$kitOmicsRecipe("metab_screen",
        recipe$biological_per_group, recipe$features_per_assay, recipe$standards_per_assay, recipe$missing_features)) },
    trx_60k_200 = { omics(); generators$writeKitTrxCorpus(directory, generators$kitOmicsRecipe("trx_screen",
        recipe$biological_per_group, recipe$features_per_assay, low_count_fraction = recipe$low_mean_fraction)) },
    lfq_fixture = { protein(); generators$writeKitFixtureCorpus(directory, "lfq", repo) },
    tmt_fixture = { protein(); generators$writeKitFixtureCorpus(directory, "tmt", repo) },
    lfq_300 = { protein(); generators$writeKitLfqCorpus(directory, generators$kitProteinRecipe("lfq", recipe$samples, recipe$proteins)) },
    tmt_300 = { protein(); generators$writeKitTmtCorpus(directory, generators$kitProteinRecipe("tmt", recipe$samples, recipe$proteins)) },
    screen = , full = {
        # The bound recipe itself (its counts are the bound ones); the encoder comes from the fixtures helper, sourced
        # (not edited) as runPerf107SyntheticCorpus sources it.
        sys.source(file.path(tools, "workflow_perf107_synthetic.R"), generators)
        sys.source(file.path(repo, "tests/testthat/helper-workflow-checkpoint-fixtures.R"), generators)
        recipe$seed <- as.integer(recipe$seed); recipe$samples <- as.integer(recipe$samples)
        recipe$universe <- as.integer(recipe$universe); recipe$counts <- as.integer(recipe$counts)
        recipe$cutoff <- as.integer(recipe$cutoff); recipe$proteins <- as.integer(recipe$proteins)
        recipe$feature_offset <- as.integer(recipe$feature_offset); recipe$annotation_tail_bytes <- as.integer(recipe$annotation_tail_bytes)
        generators$writePerf107SyntheticCohort(directory, recipe)
    },
    stop("no generator is declared for ", key))
files <- list.files(directory, full.names = TRUE)
files <- files[!basename(files) %in% c("manifest.json", "generation-status.json")]
rows <- do.call(rbind, lapply(files, function(path) {
    name <- basename(path)
    sha <- digest::digest(file = path, algo = "sha256")
    bound_sha <- bound$files[[name]]$sha256 %||% NA_character_
    data.frame(corpus = key, file = name, bytes = file.size(path), sha256 = sha, bound_sha256 = bound_sha,
        match = if (is.na(bound_sha)) NA else identical(sha, bound_sha), stringsAsFactors = FALSE)
}))
utils::write.table(rows, out, sep = "\t", row.names = FALSE, quote = FALSE, append = file.exists(out), col.names = !file.exists(out))
print(rows[, c("file", "bytes", "match")], row.names = FALSE)
