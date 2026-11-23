# KIT-173 declared control, one pair: capture a project and replay its record in one R process, on the MultiScholaR
# the library path gives (base or head), from the package root of the head tree (fixtures and helper).
#
# Usage: Rscript t7-control-pair.R <pair name> <out dir>   (R_LIBS selects the package under test)
#
# Writes <out dir>/<pair>-expectations.tsv (T7(b)'s eleven comparisons: pass or fail, the first 12 lines of a
# failure's difference) and <pair>-order.tsv (per DA table: same row order at capture and replay, rows at another
# position; per enricher call: same gene sequence, positions that differ). No feature-level row is written.

args <- commandArgs(trailingOnly = TRUE)
pair <- args[[1L]]
out <- args[[2L]]
suppressPackageStartupMessages(library(MultiScholaR))
Sys.setenv(MULTISCHOLAR_TEST_MODE = "false")
h <- new.env(parent = asNamespace("MultiScholaR"))
sys.source(file.path("tests", "testthat", "helper-go-propagation.R"), envir = h)
ns <- asNamespace("MultiScholaR")

# Recorders around the real functions: the DA tables' id order the kernel receives, the gene vectors and the hashes
# of TERM2GENE and TERM2NAME handed to enricher.
rec <- new.env()
rec$da <- list()
rec$calls <- list()
real_process <- get("processEnrichments", envir = ns)
real_enricher <- get("enricher", envir = asNamespace("clusterProfiler"))
restore_process <- h$goPropReplaceBinding("MultiScholaR", "processEnrichments", function(da_results, ...) {
    args <- list(...)
    column <- args$protein_id_column %||% "Protein.IDs"
    rec$da[[length(rec$da) + 1L]] <- lapply(da_results@da_data, function(t) as.character(t[[column]]))
    real_process(da_results, ...)
})
restore_enricher <- h$goPropReplaceBinding("clusterProfiler", "enricher", function(gene, universe = NULL, TERM2GENE,
    TERM2NAME = NA, ...) {
    rec$calls[[length(rec$calls) + 1L]] <- list(gene = as.character(gene),
        term2gene_sha256 = h$goPropTableHash(TERM2GENE), term2name_sha256 = h$goPropTableHash(TERM2NAME))
    real_enricher(gene = gene, universe = universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, ...)
})

work <- file.path(tempdir(), paste0("k173-control-", pair))
dir.create(work, recursive = TRUE)
annotations <- h$goPropAnnotate(h$goPropT7Raw(), file.path(work, "cache"))
utils::capture.output(origin <- suppressMessages(suppressWarnings(h$goPropRunProject(file.path(work, "project"),
    annotations, process_fn = get("processEnrichments", envir = ns)))))
capture <- list(da = rec$da, calls = rec$calls, evidence = h$goPropRunEvidence(origin$state, origin$paths$pathway_dir),
    provenance_names = names(origin$provenance$values), product = origin$product)
record <- jsonlite::read_json(origin$record, simplifyVector = FALSE)
rec$da <- list()
rec$calls <- list()
utils::capture.output(replayed <- suppressMessages(suppressWarnings(h$goPropReplay(file.path(work, "project"), record,
    file.path(work, "replay")))))
replay <- list(da = rec$da, calls = rec$calls, evidence = h$goPropRunEvidence(replayed$state, replayed$pathway_dir))
restore_process()
restore_enricher()

step <- Filter(function(s) identical(s$operation, "enrichment"), record$steps)[[1L]]
check <- function(name, actual, expected) {
    same <- identical(actual, expected)
    diff <- if (same) "" else gsub("\t", " ", paste(utils::head(waldo::compare(actual, expected,
        x_arg = "replay", y_arg = "capture"), 12L), collapse = " | "))
    data.frame(pair = pair, expectation = name, result = if (same) "pass" else "fail", difference = diff)
}
declared <- vapply(step$parameters, `[[`, "", "name")
hit <- Filter(function(p) identical(p$name, "go_propagation"), step$parameters)
recorded_go <- if (length(hit)) MultiScholaR:::decodeWorkflowRecordValue(hit[[1L]]$value) else NULL
expectations <- rbind(
    check("states", replayed$states, c("initial", "correlation_filtered", "enrichment_completed")),
    # T7(b)'s nine names; a head capture with parent terms also declares go_propagation (the package's own rule)
    check("record declared names", declared, MultiScholaR:::diaEnrichmentDeclaredNames(MultiScholaR:::diaEnrichmentOperation()$parameters,
        list(go_propagation = recorded_go))),
    check("provenance names", names(replayed$provenance$values), capture$provenance_names),
    check("product parameters", replayed$product$parameters, capture$product$parameters),
    check("product service", replayed$product$service, capture$product$service),
    check("product contrasts", replayed$product$contrasts, capture$product$contrasts),
    check("TERM2GENE hashes", vapply(replay$calls, `[[`, "", "term2gene_sha256"),
        vapply(capture$calls, `[[`, "", "term2gene_sha256")),
    check("TERM2NAME hashes", vapply(replay$calls, `[[`, "", "term2name_sha256"),
        vapply(capture$calls, `[[`, "", "term2name_sha256")),
    check("result tables", replay$evidence$tables, capture$evidence$tables),
    check("plot data digests", replay$evidence$plot_inputs, capture$evidence$plot_inputs),
    check("pathway files", replay$evidence$files, capture$evidence$files))
utils::write.table(expectations, file.path(out, paste0(pair, "-expectations.tsv")), sep = "\t", quote = FALSE,
    row.names = FALSE)

position_diff <- function(a, b) if (length(a) == length(b)) sum(a != b) else NA_integer_
order_rows <- list()
for (i in seq_len(max(length(capture$da), length(replay$da)))) {
    a <- capture$da[[i]] %||% list()
    b <- replay$da[[i]] %||% list()
    for (contrast in union(names(a), names(b))) {
        x <- a[[contrast]]
        y <- b[[contrast]]
        order_rows[[length(order_rows) + 1L]] <- data.frame(pair = pair, what = "DA table read by the kernel",
            index = paste(i, contrast), rows_capture = length(x), rows_replay = length(y),
            same_set = setequal(x, y), same_order = identical(x, y), rows_at_another_position = position_diff(x, y))
    }
}
for (i in seq_len(max(length(capture$calls), length(replay$calls)))) {
    x <- capture$calls[[i]]$gene %||% character()
    y <- replay$calls[[i]]$gene %||% character()
    order_rows[[length(order_rows) + 1L]] <- data.frame(pair = pair, what = "gene vector handed to enricher",
        index = as.character(i), rows_capture = length(x), rows_replay = length(y), same_set = setequal(x, y),
        same_order = identical(x, y), rows_at_another_position = position_diff(x, y))
}
utils::write.table(do.call(rbind, order_rows), file.path(out, paste0(pair, "-order.tsv")), sep = "\t", quote = FALSE,
    row.names = FALSE)
writeLines(c(paste("package", find.package("MultiScholaR")),
    paste("has_propagateGoAnnotations", exists("propagateGoAnnotations", envir = ns))),
    file.path(out, paste0(pair, "-package.txt")))
print(expectations[, 1:3])
