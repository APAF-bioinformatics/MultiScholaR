# KIT-173 declared control 2, replay process: the capture's record replayed on the installed head, in a second
# process, with the same wrappers; then T7(b)'s eleven comparisons and the order facts.
#
# Usage: Rscript t7-cross-replay.R <variant> <work dir> <out dir>
# Writes <out dir>/<variant>-expectations.tsv, <variant>-order.tsv, <variant>-geneid.tsv, <variant>-collation.tsv.
# No accession leaves the job: order facts are counts, geneID cells are given by term id, size and set equality.

args <- commandArgs(trailingOnly = TRUE)
variant <- args[[1L]]
work <- args[[2L]]
out <- args[[3L]]
suppressPackageStartupMessages(library(MultiScholaR))
Sys.setenv(MULTISCHOLAR_TEST_MODE = "false")
h <- new.env(parent = asNamespace("MultiScholaR"))
sys.source(file.path("tests", "testthat", "helper-go-propagation.R"), envir = h)
ns <- asNamespace("MultiScholaR")
capture <- readRDS(file.path(work, "capture.rds"))
record <- jsonlite::read_json(file.path(work, "record.json"), simplifyVector = FALSE)

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
utils::capture.output(replayed <- suppressMessages(suppressWarnings(h$goPropReplay(file.path(work, "project"), record,
    file.path(work, paste0("replay-", variant))))))
restore_process()
restore_enricher()
replay <- list(da = rec$da, calls = rec$calls, evidence = h$goPropRunEvidence(replayed$state, replayed$pathway_dir))

# --- T7(b)'s eleven comparisons -------------------------------------------------------------------------------------
step <- Filter(function(s) identical(s$operation, "enrichment"), record$steps)[[1L]]
declared <- vapply(step$parameters, function(p) p$name, "")
hit <- Filter(function(p) identical(p$name, "go_propagation"), step$parameters)
recorded_go <- if (length(hit)) MultiScholaR:::decodeWorkflowRecordValue(hit[[1L]]$value) else NULL
check <- function(name, actual, expected) {
    same <- identical(actual, expected)
    data.frame(variant = variant, expectation = name, result = if (same) "pass" else "fail")
}
genes_of <- function(calls, field) vapply(calls, function(call) call[[field]], "")
expectations <- rbind(
    check("states", replayed$states, c("initial", "correlation_filtered", "enrichment_completed")),
    # T7(b)'s nine names on a base record; a head record with parent terms also declares go_propagation (the
    # package's rule; added for t7-headhead.R after C3 and C4, which ran on base records)
    check("record declared names", declared, MultiScholaR:::diaEnrichmentDeclaredNames(
        MultiScholaR:::diaEnrichmentOperation()$parameters, list(go_propagation = recorded_go))),
    check("provenance names", names(replayed$provenance$values), capture$provenance_names),
    check("product parameters", replayed$product$parameters, capture$product$parameters),
    check("product service", replayed$product$service, capture$product$service),
    check("product contrasts", replayed$product$contrasts, capture$product$contrasts),
    check("TERM2GENE hashes", genes_of(replay$calls, "term2gene_sha256"), genes_of(capture$calls, "term2gene_sha256")),
    check("TERM2NAME hashes", genes_of(replay$calls, "term2name_sha256"), genes_of(capture$calls, "term2name_sha256")),
    check("result tables", replay$evidence$tables, capture$evidence$tables),
    check("plot data digests", replay$evidence$plot_inputs, capture$evidence$plot_inputs),
    check("pathway files", replay$evidence$files, capture$evidence$files))
utils::write.table(expectations, file.path(out, paste0(variant, "-expectations.tsv")), sep = "\t", quote = FALSE,
    row.names = FALSE)

# --- Order facts ----------------------------------------------------------------------------------------------------
moved <- function(a, b) if (length(a) == length(b)) sum(a != b) else NA_integer_
order_rows <- list()
for (i in seq_len(max(length(capture$da), length(replay$da)))) {
    a <- if (i <= length(capture$da)) capture$da[[i]] else list()
    b <- if (i <= length(replay$da)) replay$da[[i]] else list()
    for (contrast in union(names(a), names(b))) {
        x <- a[[contrast]] %||% character()
        y <- b[[contrast]] %||% character()
        order_rows[[length(order_rows) + 1L]] <- data.frame(variant = variant, what = "DA table read by the kernel",
            index = paste(i, contrast), rows_capture = length(x), rows_replay = length(y), same_set = setequal(x, y),
            same_order = identical(x, y), rows_at_another_position = moved(x, y))
    }
}
for (i in seq_len(max(length(capture$calls), length(replay$calls)))) {
    x <- if (i <= length(capture$calls)) capture$calls[[i]]$gene else character()
    y <- if (i <= length(replay$calls)) replay$calls[[i]]$gene else character()
    order_rows[[length(order_rows) + 1L]] <- data.frame(variant = variant, what = "gene vector handed to enricher",
        index = as.character(i), rows_capture = length(x), rows_replay = length(y), same_set = setequal(x, y),
        same_order = identical(x, y), rows_at_another_position = moved(x, y))
}
utils::write.table(do.call(rbind, order_rows), file.path(out, paste0(variant, "-order.tsv")), sep = "\t",
    quote = FALSE, row.names = FALSE)

# --- geneID cells that differ (term id, size, same set) -------------------------------------------------------------
cells <- list()
for (contrast in union(names(capture$evidence$tables), names(replay$evidence$tables))) {
    for (direction in c("up", "down")) {
        x <- capture$evidence$tables[[contrast]][[direction]]
        y <- replay$evidence$tables[[contrast]][[direction]]
        if (is.null(x) || is.null(y)) next
        for (id in intersect(as.character(x$ID), as.character(y$ID))) {
            gx <- as.character(x$geneID[match(id, x$ID)])
            gy <- as.character(y$geneID[match(id, y$ID)])
            if (identical(gx, gy)) next
            sx <- strsplit(gx, "/", fixed = TRUE)[[1L]]
            sy <- strsplit(gy, "/", fixed = TRUE)[[1L]]
            cells[[length(cells) + 1L]] <- data.frame(variant = variant, contrast = contrast, direction = direction,
                term = id, proteins = length(sx), same_set = setequal(sx, sy))
        }
    }
}
geneid <- if (length(cells)) do.call(rbind, cells) else data.frame(variant = variant, contrast = NA, direction = NA,
    term = "none differ", proteins = NA, same_set = NA)
utils::write.table(geneid, file.path(out, paste0(variant, "-geneid.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

utils::write.table(data.frame(variant = variant, process = c("capture (base)", "replay (head)"),
    collation = c(capture$collation, Sys.getlocale("LC_COLLATE")),
    package_has_propagateGoAnnotations = c(capture$has_propagateGoAnnotations,
        exists("propagateGoAnnotations", envir = ns))), file.path(out, paste0(variant, "-collation.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE)
cat("replay done\n")
