# KIT-173 declared control 2, capture process: the T7 project on the MultiScholaR that R_LIBS selects (the base),
# from the package root of the head tree (fixtures and helper). The kernel's DA tables and enricher's gene vectors are
# recorded through wrappers that call the real functions.
#
# Usage: Rscript t7-cross-capture.R <work dir>
# Writes into <work dir> (never uploaded): project/ (the store), record.json, capture.rds.

args <- commandArgs(trailingOnly = TRUE)
work <- args[[1L]]
dir.create(work, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(MultiScholaR))
Sys.setenv(MULTISCHOLAR_TEST_MODE = "false")
h <- new.env(parent = asNamespace("MultiScholaR"))
sys.source(file.path("tests", "testthat", "helper-go-propagation.R"), envir = h)
ns <- asNamespace("MultiScholaR")
rec <- new.env()
rec$da <- list()
rec$calls <- list()
real_process <- get("processEnrichments", envir = ns)
real_enricher <- get("enricher", envir = asNamespace("clusterProfiler"))
restore_process <- h$goPropReplaceBinding("MultiScholaR", "processEnrichments", function(da_results, ...) {
    column <- list(...)$protein_id_column %||% "Protein.IDs"
    rec$da[[length(rec$da) + 1L]] <- lapply(da_results@da_data, function(t) as.character(t[[column]]))
    real_process(da_results, ...)
})
restore_enricher <- h$goPropReplaceBinding("clusterProfiler", "enricher", function(gene, universe = NULL, TERM2GENE,
    TERM2NAME = NA, ...) {
    rec$calls[[length(rec$calls) + 1L]] <- list(gene = as.character(gene),
        term2gene_sha256 = h$goPropTableHash(TERM2GENE), term2name_sha256 = h$goPropTableHash(TERM2NAME))
    real_enricher(gene = gene, universe = universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, ...)
})
annotations <- h$goPropAnnotate(h$goPropT7Raw(), file.path(work, "cache"))
utils::capture.output(origin <- suppressMessages(suppressWarnings(h$goPropRunProject(file.path(work, "project"),
    annotations, process_fn = get("processEnrichments", envir = ns)))))
restore_process()
restore_enricher()
file.copy(origin$record, file.path(work, "record.json"), overwrite = TRUE)
saveRDS(list(da = rec$da, calls = rec$calls, evidence = h$goPropRunEvidence(origin$state, origin$paths$pathway_dir),
    provenance_names = names(origin$provenance$values), product = origin$product,
    collation = Sys.getlocale("LC_COLLATE"), package = find.package("MultiScholaR"),
    has_propagateGoAnnotations = exists("propagateGoAnnotations", envir = ns)), file.path(work, "capture.rds"))
cat("capture done\n")
