# KIT-193 target 3 (ticket section 6.3, 5.5): the shipped GO propagation against the Lean model's oracle table.
# tests/testdata/kit193/go-closure.tsv is written by proofs/lean/Kit193Oracles.lean from the Spec definitions the
# theorems are about (bfs = goAncestorPairs, propagateFin = propagateGoAnnotations, the root removal of
# buildUniprotGoGeneSets); the record job regenerates it and compares bytes. GO.db's reads are injected as KIT-173 T1
# injects them (helper-go-propagation.R). Kind: package (bindings.tsv rows T3). No service is called.

kit193GoTable <- function() {
    path <- testthat::test_path("..", "testdata", "kit193", "go-closure.tsv")
    utils::read.delim(path, sep = "\t", quote = "", stringsAsFactors = FALSE, colClasses = "character",
        na.strings = character())
}

#' One case's inputs and expected outputs from the long table
kit193GoCase <- function(rows) {
    pick <- function(record) rows[rows$record == record, , drop = FALSE]
    edges <- pick("edge")
    aspects <- pick("aspect")
    annots <- pick("annot")
    pairs_of <- function(x) if (nrow(x)) sort(paste(x$x, x$y, sep = "\t")) else character()
    list(
        labels = pick("labels")$x[[1L]],
        parents = data.frame(child = edges$x, parent = edges$y, relation = edges$z, stringsAsFactors = FALSE),
        aspects = stats::setNames(aspects$y, aspects$x),
        direct = data.frame(TERM = annots$x, GENE = annots$y, stringsAsFactors = FALSE),
        pair = pairs_of(pick("pair")),
        prop = pairs_of(pick("prop")),
        build = pairs_of(pick("build")))
}

#' The shipped functions' results for one case, GO.db's reads injected (BP gets the case's edges; MF and CC none)
kit193GoShipped <- function(case) {
    relations <- goPropagationRelations[[case$labels]]
    empty <- data.frame(child = character(), parent = character(), relation = character(), stringsAsFactors = FALSE)
    testthat::local_mocked_bindings(
        goParentTable = function(ontology) if (identical(ontology, "BP")) case$parents else empty,
        goDbTermAspects = function() case$aspects,
        goDbRelease = function() list(version = "kit193", date = "kit193"))
    restore <- goPropReplaceBinding("AnnotationDbi", "Term", function(x) rep("name", length(x)))
    on.exit(restore(), add = TRUE)
    pairs_of <- function(a, b) if (length(a)) sort(paste(a, b, sep = "\t")) else character()
    # the search's start terms are the case's distinct kept terms; every one is its own first pair
    starts <- unique(sub("\t.*", "", case$pair))
    pairs <- goAncestorPairs(starts, case$parents, relations)
    prop <- propagateGoAnnotations(case$direct, "BP", relations)
    by_entry <- split(case$direct$TERM, case$direct$GENE)
    annotations <- data.frame(Entry = names(by_entry),
        go_id_go_biological_process = vapply(by_entry, paste, character(1), collapse = "; "),
        go_id_go_molecular_function = NA_character_, go_id_go_cellular_compartment = NA_character_,
        stringsAsFactors = FALSE)
    sets <- buildUniprotGoGeneSets(annotations, case$labels)
    list(pair = pairs_of(pairs$term, pairs$node), prop = pairs_of(prop$TERM, prop$GENE),
        build = pairs_of(sets$term2gene$TERM, sets$term2gene$Entry))
}

kit193GoCompare <- function(table) {
    for (name in unique(table$case)) {
        case <- kit193GoCase(table[table$case == name, , drop = FALSE])
        got <- kit193GoShipped(case)
        expect_identical(got$pair, case$pair, info = paste(name, "goAncestorPairs"))
        expect_identical(got$prop, case$prop, info = paste(name, "propagateGoAnnotations"))
        expect_identical(got$build, case$build, info = paste(name, "buildUniprotGoGeneSets"))
    }
}

test_that("KIT-193 T3 pair: goAncestorPairs, propagateGoAnnotations and buildUniprotGoGeneSets equal the Lean model on every case", {
    expect_true(requireNamespace("GO.db", quietly = TRUE))
    expect_true(requireNamespace("AnnotationDbi", quietly = TRUE))
    table <- kit193GoTable()
    expect_identical(length(unique(table$case)), 211L)
    kit193GoCompare(table)
})

test_that("KIT-193 T3 pair, planted defect: a goAncestorPairs without the closure fails the comparison", {
    expect_true(requireNamespace("GO.db", quietly = TRUE))
    table <- kit193GoTable()
    case <- kit193GoCase(table[table$case == "chain", , drop = FALSE])
    testthat::local_mocked_bindings(
        goAncestorPairs = function(terms, parents, relations) data.frame(term = terms, node = terms,
            stringsAsFactors = FALSE))
    got <- kit193GoShipped(case)
    # one expectation, so expect_failure() sees exactly one failure
    testthat::expect_failure(expect_identical(got, case[c("pair", "prop", "build")]))
})
