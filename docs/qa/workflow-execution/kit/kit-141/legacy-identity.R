# KIT-141, part 1: the exported runKeggEnrichment and runReactomeEnrichment now run their pathway test through the
# integration kernel (integrationEnrichmentGsea) and keep their own fallback. This script runs the base commit's
# versions of the two functions (from git) and the head's (the installed package) on the same small fixtures, with
# the real clusterProfiler, and compares the returned tables with identical(). Cases: Reactome test, Reactome with
# the test forced to fail (the fallback), KEGG with KEGG ids in the ranked list, KEGG through ChEBI ids, KEGG with the
# test forced to fail. KEGGREST::keggList is replaced by a fixed table (no network).
#
# Usage (from the package root, full git history, head installed): Rscript legacy-identity.R <base commit> <out dir>
# Writes <out dir>/legacy-identity.tsv (case, base rows, head rows, identical) and exits 1 on any difference.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
base <- arguments[[1L]]
out_dir <- arguments[[2L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(requireNamespace("MultiScholaR", quietly = TRUE), requireNamespace("clusterProfiler", quietly = TRUE),
    requireNamespace("KEGGREST", quietly = TRUE), requireNamespace("testthat", quietly = TRUE))
namespace <- asNamespace("MultiScholaR")

# --- Fixtures ---
data_dir <- file.path(tempdir(), "kit141-legacy-identity")
dir.create(data_dir, showWarnings = FALSE)
ids <- 1:60
values <- round(seq(3, -3, length.out = 60), 1) # rounded: several tied values, so the tie-breaking noise matters
sets <- list(top = 1:10, bottom = 51:60, spread = seq(1, 60, by = 6), middle = 21:35)
set_rows <- data.frame(set = rep(names(sets), lengths(sets)), id = unlist(sets), stringsAsFactors = FALSE)
utils::write.table(data.frame(chebi_id = sprintf("CHEBI:%d", set_rows$id), reactome_id = paste0("R-HSA-", set_rows$set),
    reactome_term = paste(set_rows$set, "pathway"), organism = "Homo sapiens"),
    file.path(data_dir, "chebi_to_reactome.tab"), sep = "\t", quote = FALSE, row.names = FALSE)
numbers <- c(top = "00010", bottom = "00020", spread = "00030", middle = "00040")
utils::write.table(data.frame(pathway = paste0("path:hsa", numbers[set_rows$set]),
    compound = sprintf("cpd:C%05d", set_rows$id)), file.path(data_dir, "species_specific_pathway_to_compound_tbl.tab"),
    sep = "\t", quote = FALSE, row.names = FALSE)
kegg_list <- function(database) stats::setNames(paste(names(numbers), "pathway"), paste0("path:map", numbers))
mapping_table <- data.frame(KEGG = sprintf("C%05d", ids), ChEBI = as.character(ids), stringsAsFactors = FALSE)
chebi_ranked <- stats::setNames(values, sprintf("CHEBI:%d", ids))
kegg_ranked <- stats::setNames(values, sprintf("cpd:C%05d", ids))

# --- The two sides: each function re-homed in an environment whose getProjectPaths points at the fixtures ---
side_env <- function() {
    env <- new.env(parent = namespace)
    env$getProjectPaths <- function(omic_type, experiment_label) list(data_dir = data_dir)
    env
}
base_env <- side_env()
for (file in c("R/func_multiomics_enrich_kegg_helpers.R", "R/func_multiomics_enrich_reactome_helpers.R")) {
    source_text <- system2("git", c("show", paste0(base, ":", file)), stdout = TRUE)
    eval(parse(text = source_text), envir = base_env)
}
head_env <- side_env()
for (name in c("runKeggEnrichment", "runReactomeEnrichment")) {
    f <- get(name, envir = namespace)
    environment(f) <- head_env
    assign(name, f, envir = head_env)
}

run_case <- function(env, case) {
    call <- switch(case$fn,
        reactome = function() env$runReactomeEnrichment(case$ranked, mapping_table, list(), "metabolomics", "demo",
            "assay", reactome_organism = "Homo sapiens"),
        kegg = function() env$runKeggEnrichment(case$ranked, mapping_table, list(), "metabolomics", "demo", "assay",
            kegg_species_code = "hsa"))
    run <- function() {
        set.seed(1) # the same random state before each call on each side
        suppressMessages(suppressWarnings(call()))
    }
    if (isTRUE(case$fail)) {
        testthat::with_mocked_bindings(testthat::with_mocked_bindings(run(), keggList = kegg_list, .package = "KEGGREST"),
            GSEA = function(...) stop("forced failure of the test"), .package = "clusterProfiler")
    } else {
        testthat::with_mocked_bindings(run(), keggList = kegg_list, .package = "KEGGREST")
    }
}

cases <- list(
    list(name = "reactome", fn = "reactome", ranked = chebi_ranked),
    list(name = "reactome_fallback", fn = "reactome", ranked = chebi_ranked, fail = TRUE),
    list(name = "kegg_direct", fn = "kegg", ranked = kegg_ranked),
    list(name = "kegg_via_chebi", fn = "kegg", ranked = chebi_ranked),
    list(name = "kegg_direct_fallback", fn = "kegg", ranked = kegg_ranked, fail = TRUE)
)
rows <- lapply(cases, function(case) {
    base_result <- run_case(base_env, case)
    head_result <- run_case(head_env, case)
    data.frame(case = case$name, base_rows = nrow(base_result), head_rows = nrow(head_result),
        identical = identical(base_result, head_result), stringsAsFactors = FALSE)
})
table <- do.call(rbind, rows)
utils::write.table(table, file.path(out_dir, "legacy-identity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(paste("base:", base), paste("head:", system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
    paste("clusterProfiler", as.character(utils::packageVersion("clusterProfiler"))), R.version.string),
    file.path(out_dir, "legacy-identity-meta.txt"))
print(table)
if (!all(table$identical)) quit(status = 1L)
