# KIT-173: Gene Ontology parent terms on the proteomics clusterProfiler route (T1 to T11 of the ticket, section 10).
# Fixtures: tests/testdata/go-propagation/ (README.md: sources, releases, SHA-256, licences). No service is called.
# Helpers: helper-go-propagation.R.
# KIT-234: T10's loop over the five DIA workbooks removed (the workbooks are archived).

goPropSkipUnlessGoDb <- function() {
    testthat::skip_if_not_installed("GO.db")
    testthat::skip_if_not_installed("AnnotationDbi")
    testthat::skip_if_not_installed("clusterProfiler")
}

#' GO.db's own five-relation closure of one aspect (`GO*ANCESTOR`), the root `all` removed, with each term itself
goPropAncestorPairs <- function(aspect) {
    map <- switch(aspect, BP = GO.db::GOBPANCESTOR, MF = GO.db::GOMFANCESTOR, CC = GO.db::GOCCANCESTOR)
    pairs <- AnnotationDbi::toTable(map)
    pairs <- data.frame(term = pairs[[1L]], ancestor = pairs[[2L]], stringsAsFactors = FALSE)
    pairs[pairs$ancestor != "all", , drop = FALSE]
}

#' The recursive-SQL closure of one aspect over the given relation labels on GO.db's SQLite (root `all` removed)
goPropSqlPairs <- function(aspect, relations) {
    table <- sprintf("go_%s_parents", tolower(aspect))
    rel <- paste(sprintf("'%s'", relations), collapse = ", ")
    sql <- sprintf(paste("WITH RECURSIVE anc(c, a) AS (SELECT _id, _parent_id FROM %1$s WHERE relationship_type IN (%2$s)",
        "UNION SELECT anc.c, p._parent_id FROM anc JOIN %1$s p ON p._id = anc.a WHERE p.relationship_type IN (%2$s))",
        "SELECT t1.go_id AS term, t2.go_id AS ancestor FROM anc JOIN go_term t1 ON t1._id = anc.c",
        "JOIN go_term t2 ON t2._id = anc.a WHERE t2.go_id <> 'all'"), table, rel)
    DBI::dbGetQuery(GO.db::GO_dbconn(), sql)
}

goPropTwo <- c("isa", "part of")
goPropFive <- c("isa", "part of", "regulates", "positively regulates", "negatively regulates")

# --- Fixture hashes ---------------------------------------------------------------------------------------------------

test_that("KIT-173 fixtures: every committed fixture matches the SHA-256 its README states", {
    readme <- readLines(goPropFixtureDir("README.md"), warn = FALSE)
    rows <- regmatches(readme, regexec("^\\| `([^`]+)` \\| `([0-9a-f]{64})` \\|", readme))
    rows <- Filter(length, rows)
    expect_gte(length(rows), 5L)
    for (row in rows) {
        expect_true(file.exists(goPropFixtureDir(row[[2L]])), info = row[[2L]])
        expect_identical(digest::digest(file = goPropFixtureDir(row[[2L]]), algo = "sha256"), row[[3L]], info = row[[2L]])
    }
})

# --- T1: toy graph through the function's graph step -------------------------------------------------------------------

test_that("KIT-173 T1: a hand-built 8-term graph propagates to the hand-computed closure under each relation set", {
    goPropSkipUnlessGoDb()
    ids <- sprintf("GO:90000%02d", 1:8)
    parents <- data.frame(
        child = ids[c(1, 2, 3, 4, 5, 6, 7, 8, 8)],
        parent = c("all", ids[c(1, 2, 2, 3, 5, 4, 6, 4)]),
        relation = c("isa", "isa", "isa", "part of", "regulates", "isa", "positively regulates", "isa", "part of"),
        stringsAsFactors = FALSE)
    testthat::local_mocked_bindings(
        goParentTable = function(ontology) parents,
        goDbTermAspects = function() stats::setNames(rep("BP", 8L), ids),
        goDbRelease = function() list(version = "toy", date = "toy"))
    direct <- data.frame(TERM = ids[c(3, 5, 6, 7, 8)], GENE = c("g3", "g5", "g6", "g7", "g8"), stringsAsFactors = FALSE)
    expected <- function(sets) {
        out <- do.call(rbind, lapply(names(sets), function(t) data.frame(TERM = t, GENE = sets[[t]],
            stringsAsFactors = FALSE)))
        out <- out[order(out$TERM, out$GENE, method = "radix"), , drop = FALSE]
        rownames(out) <- NULL
        out
    }
    two <- expected(stats::setNames(list(c("g3", "g8"), c("g3", "g8"), "g3", "g8", c("g5", "g6", "g8"), c("g6", "g8"),
        "g7", "g8"), ids))
    five <- expected(stats::setNames(list(c("g3", "g5", "g6", "g7", "g8"), c("g3", "g5", "g6", "g7", "g8"),
        c("g3", "g5", "g6", "g8"), c("g7", "g8"), c("g5", "g6", "g8"), c("g6", "g8"), "g7", "g8"), ids))
    got_two <- propagateGoAnnotations(direct, "BP", goPropTwo)
    got_five <- propagateGoAnnotations(direct, "BP", goPropFive)
    expect_identical(as.data.frame(got_two)[, c("TERM", "GENE")], two)
    expect_identical(as.data.frame(got_five)[, c("TERM", "GENE")], five)
    # The regulates child's gene (g5, on GO:9000005, which regulates GO:9000003) reaches GO:9000003 only over five.
    expect_false("g5" %in% got_two$GENE[got_two$TERM == ids[[3L]]])
    expect_true("g5" %in% got_five$GENE[got_five$TERM == ids[[3L]]])
    expect_false("all" %in% c(got_two$TERM, got_five$TERM))
    expect_identical(attr(got_five, "relations"), goPropFive)
    expect_identical(attr(got_five, "pairs_direct"), 5L)
    expect_identical(attr(got_five, "pairs_propagated"), nrow(five))
    expect_identical(goPropagationRelations, list(is_a_part_of = goPropTwo, is_a_part_of_regulates = goPropFive))
    expect_identical(goPropagationShipped, "is_a_part_of_regulates")
    expect_error(propagateGoAnnotations(direct, "XX", goPropFive), class = "go_propagation_invalid")
    expect_error(propagateGoAnnotations(direct, "BP", "is_a"), class = "go_propagation_invalid")
    expect_error(propagateGoAnnotations(direct, "BP", character()), class = "go_propagation_invalid")
})

# --- T2: the real GO graph, exact ----------------------------------------------------------------------------------------

test_that("KIT-173 T2(a): over five relations the function's ancestors of every GO term equal GO.db's GO*ANCESTOR", {
    goPropSkipUnlessGoDb()
    total <- 0L
    for (aspect in c("BP", "MF", "CC")) {
        terms <- AnnotationDbi::Ontology(GO.db::GOTERM)
        terms <- sort(names(terms)[terms == aspect])
        total <- total + length(terms)
        got <- propagateGoAnnotations(data.frame(TERM = terms, GENE = terms, stringsAsFactors = FALSE), aspect,
            goPropFive)
        ours <- sort(paste(got$GENE, got$TERM)[got$GENE != got$TERM])
        reference <- goPropAncestorPairs(aspect)
        expect_identical(ours, sort(paste(reference$term, reference$ancestor)), info = aspect)
        expect_identical(attr(got, "not_in_go_db"), 0L, info = aspect)
    }
    expect_identical(total, 38739L)
})

test_that("KIT-173 T2(b): over is a and part of the term-ancestor pairs per aspect equal a recursive SQL count", {
    goPropSkipUnlessGoDb()
    expected <- c(BP = 352532L, MF = 64418L, CC = 42075L)
    for (aspect in names(expected)) {
        terms <- AnnotationDbi::Ontology(GO.db::GOTERM)
        terms <- sort(names(terms)[terms == aspect])
        got <- propagateGoAnnotations(data.frame(TERM = terms, GENE = terms, stringsAsFactors = FALSE), aspect,
            goPropTwo)
        sql <- goPropSqlPairs(aspect, goPropTwo)
        expect_identical(nrow(sql), expected[[aspect]], info = aspect)
        ours <- sort(paste(got$GENE, got$TERM)[got$GENE != got$TERM])
        expect_identical(ours, sort(paste(sql$term, sql$ancestor)), info = aspect)
    }
})

test_that("KIT-173 T2(c, d): three proteins on GO:0006110 reach exactly 19 terms over two relations and 76 over five", {
    goPropSkipUnlessGoDb()
    direct <- data.frame(TERM = "GO:0006110", GENE = c("Q1", "Q2", "Q3"), stringsAsFactors = FALSE)
    two_set <- c("GO:0006110", goPropSqlPairs("BP", goPropTwo)$ancestor[goPropSqlPairs("BP", goPropTwo)$term ==
        "GO:0006110"])
    five_set <- c("GO:0006110", setdiff(GO.db::GOBPANCESTOR[["GO:0006110"]], "all"))
    expect_length(two_set, 19L)
    expect_length(five_set, 76L)
    for (case in list(list(relations = goPropTwo, set = two_set), list(relations = goPropFive, set = five_set))) {
        got <- propagateGoAnnotations(direct, "BP", case$relations)
        for (gene in c("Q1", "Q2", "Q3")) expect_identical(sort(got$TERM[got$GENE == gene]), sort(case$set))
    }
    expect_true("GO:0008150" %in% two_set)
    expect_false(any(c("GO:0006096", "GO:0008152") %in% two_set))
    expect_true(all(c("GO:0006096", "GO:0008152") %in% five_set))
    # (d) the route: with "none" only GO:0006110 holds them; with parent terms, (c)'s sets without the BP root.
    # A fourth protein carries one MF and one CC term: the direct sets' code (unchanged from base) fails on an aspect
    # column with no annotation at all.
    annotations <- data.frame(Entry = c("Q1", "Q2", "Q3", "Q4"),
        go_id_go_biological_process = c(rep("GO:0006110", 3L), NA), go_term_go_biological_process = NA_character_,
        go_id_go_molecular_function = c(NA, NA, NA, "GO:0005515"), go_term_go_molecular_function = NA_character_,
        go_id_go_cellular_compartment = c(NA, NA, NA, "GO:0005739"), go_term_go_cellular_compartment = NA_character_,
        stringsAsFactors = FALSE)
    none <- buildUniprotGoGeneSets(annotations, "none")$term2gene
    expect_identical(sort(unique(none$TERM[none$Entry %in% c("Q1", "Q2", "Q3")])), "GO:0006110")
    for (case in list(list(value = "is_a_part_of", set = two_set), list(value = "is_a_part_of_regulates",
        set = five_set))) {
        sets <- buildUniprotGoGeneSets(annotations, case$value)
        for (gene in c("Q1", "Q2", "Q3")) expect_identical(sort(sets$term2gene$TERM[sets$term2gene$Entry == gene]),
            sort(setdiff(case$set, "GO:0008150")), info = case$value)
    }
})

# --- T3: lost identifiers ------------------------------------------------------------------------------------------------

test_that("KIT-173 T3(a): ids GO.db does not hold, or of another aspect, are dropped and counted by the function", {
    goPropSkipUnlessGoDb()
    direct <- data.frame(TERM = c("GO:0006110", "GO:9999999", "GO:0005739"), GENE = "Q1", stringsAsFactors = FALSE)
    got <- propagateGoAnnotations(direct, "BP", goPropFive)
    expect_identical(attr(got, "not_in_go_db"), 1L)
    expect_identical(attr(got, "not_in_go_db_examples"), "GO:9999999")
    expect_identical(attr(got, "not_in_aspect"), 1L)
    expect_false(any(c("GO:9999999", "GO:0005739") %in% got$TERM))
    expect_identical(attr(got, "go_db_version"), as.character(utils::packageVersion("GO.db")))
})

#' The fresh head project (T3(b), T7(c), T8, T9, T11(a)): the T7 fixture with three lost ids, made once per file.
goPropFresh <- local({
    cached <- NULL
    function() {
        if (is.null(cached)) {
            work <- tempfile("k173-fresh-")
            dir.create(work)
            annotations <- goPropAnnotate(goPropT7Raw(), file.path(work, "cache"))
            acc <- new.env()
            restore <- goPropRecordEnricher(acc)
            on.exit(restore(), add = TRUE)
            withr::with_envvar(c(MULTISCHOLAR_TEST_MODE = "false"), utils::capture.output(result <- suppressMessages(
                suppressWarnings(goPropRunProject(file.path(work, "project"), annotations)))))
            cached <<- c(result, list(annotations = annotations, enricher = acc$calls, work = work,
                record_list = jsonlite::read_json(result$record, simplifyVector = FALSE)))
        }
        cached
    }
})

goPropStep <- function(record) Filter(function(s) identical(s$operation, "enrichment"), record$steps)[[1L]]

test_that("KIT-173 T3(b): getUniprotAnnotations counts one obsolete, one secondary and one unknown id; the step records them", {
    goPropSkipUnlessGoDb()
    fresh <- goPropFresh()
    lost <- attr(fresh$annotations, "go_ids_not_in_go_db", exact = TRUE)
    expect_identical(lost[c("obsolete", "secondary", "unknown", "total")],
        list(obsolete = 1L, secondary = 1L, unknown = 1L, total = 3L))
    expect_identical(lost$go_db_version, as.character(utils::packageVersion("GO.db")))
    info <- GO.db::GO_dbInfo()
    expect_identical(lost$go_source_date, info$value[match("GOSOURCEDATE", info$name)])
    outcomes <- goPropStep(fresh$record_list)$outcomes
    expect_equal(unlist(outcomes$go_ids_not_in_go_db), c(obsolete = 1, secondary = 1, unknown = 1))
    expect_equal(outcomes$go_ids_not_in_go_db_total, 3)
    expect_identical(outcomes$go_annotation_go_db_version, as.character(utils::packageVersion("GO.db")))
    # A table made before the count records "not recorded".
    plain <- fresh$annotations
    attr(plain, "go_ids_not_in_go_db") <- NULL
    da <- methods::new("da_results_for_enrichment", contrasts = tibble::tibble(contrast = "x"),
        da_data = list(x = data.frame(uniprot_acc = plain$Entry, stringsAsFactors = FALSE)),
        design_matrix = data.frame(sample = "S1", stringsAsFactors = FALSE))
    old <- workflowEnrichmentGoOutcomes(plain, "is_a_part_of_regulates", da, "uniprot_acc")
    expect_identical(old$go_ids_not_in_go_db, "not recorded")
    expect_identical(old$go_ids_not_in_go_db_total, "not recorded")
    expect_identical(old$go_annotation_go_db_version, "not recorded")
})

# --- T4: row order -----------------------------------------------------------------------------------------------------

test_that("KIT-173 T4: shuffled annotation rows give identical TERM2GENE and identical enricher p-values", {
    goPropSkipUnlessGoDb()
    annotations <- goPropConvert(goPropReadTsv("t7-annotation-ecoli.tsv"))
    genes <- sort(annotations$Entry)[1:25]
    for (value in c("is_a_part_of", "is_a_part_of_regulates")) {
        reference <- buildUniprotGoGeneSets(annotations, value)
        p_reference <- clusterProfiler::enricher(genes, universe = annotations$Entry, TERM2GENE = reference$term2gene,
            TERM2NAME = reference$term2name, pvalueCutoff = 1, qvalueCutoff = 1)@result$pvalue
        for (seed in c(1L, 2L)) {
            set.seed(seed)
            shuffled <- annotations[sample.int(nrow(annotations)), , drop = FALSE]
            sets <- buildUniprotGoGeneSets(shuffled, value)
            expect_identical(sets$term2gene, reference$term2gene, info = paste(value, seed))
            p <- clusterProfiler::enricher(genes, universe = annotations$Entry, TERM2GENE = sets$term2gene,
                TERM2NAME = sets$term2name, pvalueCutoff = 1, qvalueCutoff = 1)@result$pvalue
            expect_identical(p, p_reference, info = paste(value, seed))
        }
    }
})

# --- T5: the model organism (section 5.2) ------------------------------------------------------------------------------

test_that("KIT-173 T5: yeast sets agree with g:Profiler's propagated sets (median C and P >= 0.90; C gains >= 0.20)", {
    goPropSkipUnlessGoDb()
    raw <- goPropReadTsv("uniprot-yeast-559292.tsv.gz")
    annotations <- goPropConvert(raw)
    gmt <- goPropReadGmt()
    pathway_dir <- withr::local_tempdir()
    # The shipped value through the kernel's default (the workbook-shaped call); "none" explicitly where the kernel
    # takes the argument (on the base commit the default is the direct sets).
    shipped <- goPropYeastTerm2Gene(annotations, NULL, pathway_dir)
    none <- if ("go_propagation" %in% names(formals(processEnrichments))) {
        goPropYeastTerm2Gene(annotations, "none", pathway_dir)
    } else shipped
    agreement <- goPropAgreement(list(none = none, is_a_part_of_regulates = shipped), raw, gmt, "is_a_part_of_regulates")
    by_value <- split(agreement, agreement$value)
    shipped_c <- stats::median(by_value$is_a_part_of_regulates$C)
    shipped_p <- stats::median(by_value$is_a_part_of_regulates$P, na.rm = TRUE)
    none_c <- stats::median(by_value$none$C)
    expect_gt(nrow(by_value$is_a_part_of_regulates), 1000L)
    expect_gte(shipped_c, 0.90)
    expect_gte(shipped_p, 0.90)
    expect_gte(shipped_c - none_c, 0.20)
})

# --- T6: memory and time of the route's gene-set build -------------------------------------------------------------------

test_that("KIT-173 T6: the yeast gene-set build stays under 200 MB of added peak RSS and 30 s", {
    goPropSkipUnlessGoDb()
    skip_if_not(file.exists("/proc/self/clear_refs") && file.exists("/proc/self/status"), "Linux /proc needed")
    annotations <- goPropConvert(goPropReadTsv("uniprot-yeast-559292.tsv.gz"))
    status_kb <- function(field) {
        line <- grep(paste0("^", field, ":"), readLines("/proc/self/status"), value = TRUE)
        as.numeric(sub("^[^0-9]*([0-9]+).*$", "\\1", line))
    }
    invisible(buildUniprotGoGeneSets(annotations[1:10, ], "is_a_part_of_regulates")) # GO.db reads made once
    gc()
    writeLines("5", "/proc/self/clear_refs")
    before <- status_kb("VmRSS")
    seconds <- system.time(sets <- buildUniprotGoGeneSets(annotations, "is_a_part_of_regulates"))[["elapsed"]]
    peak_mb <- (status_kb("VmHWM") - before) / 1024
    message(sprintf("KIT-173 T6: %.1f s, added peak RSS %.1f MB, %d pairs", seconds, peak_mb, nrow(sets$term2gene)))
    expect_lt(seconds, 30)
    expect_lt(peak_mb, 200)
})

# --- T7: replay ----------------------------------------------------------------------------------------------------------

test_that("KIT-173 T7(b): a record made on the base commit replays with the same names, sets, tables and files (regression)", {
    goPropSkipUnlessGoDb()
    # The base capture with its project store: script t7-replay.R makes it on the runner (base code installed beside
    # the head) and points KIT173_REPLAY_BASE at it. The committed replay-base/ holds the record and evidence without
    # the store (R1's upload withheld it), and a replay needs the record's own store (states are checked by entity id).
    base_dir <- Sys.getenv("KIT173_REPLAY_BASE", goPropFixtureDir("replay-base"))
    skip_if_not(file.exists(file.path(base_dir, "project.tar")), "the base project store is present only in t7-replay.R")
    scratch <- withr::local_tempdir()
    project <- file.path(scratch, "origin")
    dir.create(project)
    utils::untar(file.path(base_dir, "project.tar"), exdir = project, tar = "internal")
    record <- jsonlite::read_json(file.path(base_dir, "record.json"), simplifyVector = FALSE)
    base_evidence <- readRDS(file.path(base_dir, "evidence.rds"))
    base_product <- readRDS(file.path(base_dir, "product.rds"))
    base_calls <- utils::read.delim(file.path(base_dir, "enricher.tsv"), stringsAsFactors = FALSE)
    acc <- new.env()
    restore <- goPropRecordEnricher(acc)
    withr::defer(restore())
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "false")
    utils::capture.output(replayed <- suppressMessages(suppressWarnings(goPropReplay(project, record,
        file.path(scratch, "replay")))))
    expect_identical(replayed$states, c("initial", "correlation_filtered", "enrichment_completed"))
    step <- goPropStep(record)
    expect_identical(vapply(step$parameters, `[[`, "", "name"), utils::head(diaEnrichmentOperation()$parameters, 9L))
    expect_identical(names(replayed$provenance$values), base_product$provenance_names)
    expect_identical(replayed$product$parameters, base_product$parameters)
    expect_identical(replayed$product$service, base_product$service)
    expect_identical(replayed$product$contrasts, base_product$contrasts)
    expect_identical(vapply(acc$calls, `[[`, "", "term2gene_sha256"), base_calls$term2gene_sha256)
    expect_identical(vapply(acc$calls, `[[`, "", "term2name_sha256"), base_calls$term2name_sha256)
    evidence <- goPropRunEvidence(replayed$state, replayed$pathway_dir)
    # A4 as restated (VAL-61): a record made before the geneID fix is compared with its geneID cells as sets (their
    # order was set per R process); every other column exactly. Its results TSV files differ only in that text, so they
    # are compared by name; the plot files by existence as before.
    as_sets <- function(tables) lapply(tables, function(contrast) lapply(contrast, function(table) {
        if (is.null(table) || !"geneID" %in% names(table)) return(table)
        table$geneID <- vapply(strsplit(as.character(table$geneID), "/", fixed = TRUE),
            function(g) paste(sort(g, method = "radix"), collapse = "/"), character(1))
        table
    }))
    expect_identical(as_sets(evidence$tables), as_sets(base_evidence$tables))
    expect_identical(evidence$plot_inputs, base_evidence$plot_inputs)
    expect_identical(evidence$files$file, base_evidence$files$file)
    expect_identical(evidence$files$sha256[!grepl("\\.tsv$", evidence$files$file)],
        base_evidence$files$sha256[!grepl("\\.tsv$", base_evidence$files$file)])
})

test_that("KIT-173 T7(c): the same DA and annotation run fresh declare go_propagation and pass other sets", {
    goPropSkipUnlessGoDb()
    fresh <- goPropFresh()
    step <- goPropStep(fresh$record_list)
    names <- vapply(step$parameters, `[[`, "", "name")
    # KIT-203: the default significance rule is declared after it. KIT-205: then the fixed background rule marker.
    expect_identical(names, c(utils::head(diaEnrichmentOperation()$parameters, 9L), "go_propagation",
        "significance_rule", "enrichment_background"))
    value <- Filter(function(p) identical(p$name, "go_propagation"), step$parameters)[[1L]]$value
    expect_identical(decodeWorkflowRecordValue(value), "is_a_part_of_regulates")
    base_dir <- goPropFixtureDir("replay-base")
    skip_if_not(file.exists(file.path(base_dir, "enricher.tsv")), "replay-base not captured yet (R1 makes it)")
    base_calls <- utils::read.delim(file.path(base_dir, "enricher.tsv"), stringsAsFactors = FALSE)
    expect_false(any(vapply(fresh$enricher, `[[`, "", "term2gene_sha256") %in% base_calls$term2gene_sha256))
})

# --- T8: output labels and background ----------------------------------------------------------------------------------

test_that("KIT-173 T8: a propagated run labels parent terms from GO.db, tests no root, and drops the root-only protein from the background", {
    goPropSkipUnlessGoDb()
    fresh <- goPropFresh()
    roles <- goPropReadTsv("t7-fixture-entries.tsv")
    direct <- goPropConvert(goPropReadTsv("t7-annotation-ecoli.tsv"))
    direct_terms <- unique(unlist(strsplit(stats::na.omit(unlist(direct[c("go_id_go_biological_process",
        "go_id_go_molecular_function", "go_id_go_cellular_compartment")])), "; ", fixed = TRUE)))
    saved <- do.call(rbind, lapply(list.files(fresh$paths$pathway_dir, "_enrichment_results\\.tsv$", full.names = TRUE),
        function(f) readr::read_tsv(f, show_col_types = FALSE, col_types = readr::cols(.default = "c"))))
    expect_true(any(!saved$ID %in% direct_terms))
    expect_false(any(saved$source == "Other"))
    aspects <- AnnotationDbi::Ontology(GO.db::GOTERM)
    expect_identical(unname(saved$source), paste0("GO:", unname(aspects[saved$ID])))
    expect_identical(unname(saved$term), unname(AnnotationDbi::Term(saved$ID)))
    expect_false(any(saved$ID %in% c("GO:0008150", "GO:0003674", "GO:0005575")))
    # The background: the root-only protein is counted in "none"'s BgRatio denominator, not in the shipped value's.
    root_only <- roles$entry[roles$role == "root_only"]
    expect_length(root_only, 1L)
    table <- data.frame(Protein.IDs = sort(direct$Entry), fdr_qvalue = 0.01,
        log2FC = ifelse(sort(direct$Entry) %in% roles$entry[roles$role == "p_term"], 3, 0), stringsAsFactors = FALSE)
    da <- methods::new("da_results_for_enrichment", contrasts = tibble::tibble(contrast = "t8"),
        da_data = list(t8 = table), design_matrix = data.frame(sample = "S1", stringsAsFactors = FALSE))
    denominators <- vapply(c("none", "is_a_part_of_regulates"), function(value) {
        results <- NULL
        utils::capture.output(results <- suppressMessages(suppressWarnings(processEnrichments(da, "999999",
            q_cutoff = 0.05, pathway_dir = withr::local_tempdir(), go_annotations = direct, exclude_iea = FALSE,
            protein_id_column = "Protein.IDs", contrast_names = "t8", go_propagation = value))))
        bg <- results@enrichment_data$t8$up@result$BgRatio
        as.integer(unique(sub("^[0-9]+/", "", bg)))
    }, integer(1))
    expect_identical(unname(denominators[["none"]] - denominators[["is_a_part_of_regulates"]]), 1L)
    expect_equal(goPropStep(fresh$record_list)$outcomes$go_root_only_proteins, 1)
})

# --- T9: records and replay warnings -----------------------------------------------------------------------------------

test_that("KIT-173 T9: a 9606 run declares the nine names and passes no go_propagation; a clusterProfiler record declares it with its outcomes", {
    goPropSkipUnlessGoDb()
    # The g:Profiler route: the argument builder adds nothing unless given, and the declared names drop it.
    expect_identical(diaEnrichmentDeclaredNames(diaEnrichmentOperation()$parameters, list()),
        utils::head(diaEnrichmentOperation()$parameters, 9L))
    expect_identical(diaEnrichmentDeclaredNames(diaEnrichmentOperation()$parameters, list(go_propagation = "none")),
        utils::head(diaEnrichmentOperation()$parameters, 9L))
    method <- resolveProtEnrichAnalysisMethod("9606", buildProtEnrichSupportedOrganisms())
    expect_identical(method$analysisMethod, "gprofiler2")
    fresh <- goPropFresh()
    outcomes <- goPropStep(fresh$record_list)$outcomes
    expect_true(all(c("go_db_version", "go_source_date", "go_relations", "go_ids_not_in_go_db",
        "go_ids_not_in_go_db_total", "go_annotation_go_db_version", "go_annotation_go_source_date", "go_pairs_direct",
        "go_pairs_propagated", "go_root_only_proteins", "go_unreviewed_share") %in% names(outcomes)))
    expect_identical(unlist(outcomes$go_relations), c("is a", "part of", "regulates", "positively regulates",
        "negatively regulates"))
    expect_equal(outcomes$go_unreviewed_share, 0)
    expect_identical(outcomes$go_db_version, as.character(utils::packageVersion("GO.db")))
})

test_that("KIT-173 T9: a replay under another GO.db release warns with class workflow_replay_resource", {
    goPropSkipUnlessGoDb()
    fresh <- goPropFresh()
    testthat::local_mocked_bindings(goDbRelease = function() list(version = "0.0.0", date = "2000-01-01"))
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "false")
    classes <- character()
    utils::capture.output(suppressMessages(withCallingHandlers(goPropReplay(file.path(fresh$work, "project"),
        fresh$record_list, withr::local_tempdir()), warning = function(w) {
            classes <<- c(classes, class(w)[[1L]])
            invokeRestart("muffleWarning")
        })))
    expect_true("workflow_replay_resource" %in% classes)
})

# --- T10: script users and the electronic-annotation warning ---------------------------------------------------------

test_that("KIT-173 T10: exclude_iea = TRUE warns on the UniProt route; the workbook-shaped call propagates", {
    goPropSkipUnlessGoDb()
    annotations <- data.frame(Entry = sprintf("Q%02d", 1:12), go_id_go_biological_process = "GO:0006110",
        go_term_go_biological_process = "regulation of glycolytic process", go_id_go_molecular_function = NA_character_,
        go_term_go_molecular_function = NA_character_, go_id_go_cellular_compartment = NA_character_,
        go_term_go_cellular_compartment = NA_character_, stringsAsFactors = FALSE)
    table <- data.frame(Protein.IDs = annotations$Entry, fdr_qvalue = 0.01, log2FC = rep(c(2, 0), each = 6L),
        stringsAsFactors = FALSE)
    da <- methods::new("da_results_for_enrichment", contrasts = tibble::tibble(contrast = "t10"),
        da_data = list(t10 = table), design_matrix = data.frame(sample = "S1", stringsAsFactors = FALSE))
    acc <- new.env()
    restore <- goPropRecordEnricher(acc)
    withr::defer(restore())
    warned <- character()
    utils::capture.output(suppressMessages(withCallingHandlers(processEnrichments(da, "999999", q_cutoff = 0.05,
        pathway_dir = withr::local_tempdir(), go_annotations = annotations, exclude_iea = TRUE,
        protein_id_column = "Protein.IDs", contrast_names = "t10"), warning = function(w) {
            warned <<- c(warned, class(w))
            invokeRestart("muffleWarning")
        })))
    expect_true("multischolar_enrichment_iea_not_excludable" %in% warned)
    five_set <- c("GO:0006110", setdiff(GO.db::GOBPANCESTOR[["GO:0006110"]], c("all", "GO:0008150")))
    expect_identical(sort(unique(acc$calls[[1L]]$term2gene[[1L]])), sort(five_set))
})

# --- T11: report text ---------------------------------------------------------------------------------------------------

goPropTemplates <- c("LFQ_report.rmd", "TMT_report.rmd", "DIANN_report.rmd", "DIANN_limpa_report.rmd")

#' A template's `FormatGOEnrichment`, extracted with knitr::purl (an `.R` file, the base chunk, is read as it is) and
#' evaluated
goPropFormatter <- function(template_path) {
    code <- template_path
    if (grepl("\\.rmd$", template_path, ignore.case = TRUE)) {
        lines <- readLines(template_path, warn = FALSE, encoding = "UTF-8")
        start <- grep("^```\\{r GO Enrichment Section", lines)
        stopifnot(length(start) == 1L)
        end <- start + which(lines[-seq_len(start)] == "```")[[1L]]
        chunk <- tempfile(fileext = ".Rmd")
        writeLines(lines[start:end], chunk, useBytes = TRUE)
        code <- tempfile(fileext = ".R")
        suppressMessages(knitr::purl(chunk, output = code, quiet = TRUE, documentation = 0))
    }
    exprs <- parse(code, keep.source = FALSE, encoding = "UTF-8")
    at <- which(vapply(exprs, function(e) is.call(e) && identical(e[[1L]], as.name("<-")) &&
        identical(e[[2L]], as.name("FormatGOEnrichment")), logical(1)))
    stopifnot(length(at) == 1L)
    env <- new.env(parent = globalenv())
    eval(exprs[[at]], env)
    env$FormatGOEnrichment
}

goPropHeadFormatter <- function(template) {
    goPropFormatter(system.file("reports", "proteomics", template, package = "MultiScholaR"))
}

goPropNumbers <- function(source, up = 0, down = 0, q = 0.05, extra = list()) {
    numbers <- c(list(ui_enrichment_organism_selected = "999999", ui_enrichment_up_lfc = up, ui_enrichment_down_lfc = down,
        ui_enrichment_q_value_cutoff = q), if (!is.null(source)) list(ui_enrichment_database_source = source), extra)
    numbers
}

goPropD7 <- "Over-representation was tested with clusterProfiler (Yu et al., 2012) against Gene Ontology annotations from UniProt."

test_that("KIT-173 T11(a, b, c): the clusterProfiler paragraph states the record's sets, family, q-value and background", {
    goPropSkipUnlessGoDb()
    fresh <- goPropFresh()
    params <- workflowRecordReportParameters(fresh$record, list(ui_enrichment_database_source = "clusterprofiler"))
    expect_identical(params$enrichment_correction_method, "BH")
    info <- GO.db::GO_dbInfo()
    date <- info$value[match("GOSOURCEDATE", info$name)]
    version <- as.character(utils::packageVersion("GO.db"))
    d3 <- paste0("Gene Ontology protein sets were built from UniProt annotations. Each term's set also includes the ",
        "proteins annotated to its more specific terms (Gene Ontology release ", date, ", GO.db ", version,
        "; relations followed: is a, part of, regulates, positively regulates, negatively regulates).")
    for (template in goPropTemplates) {
        format <- goPropHeadFormatter(template)
        record_numbers <- goPropNumbers("clusterprofiler", extra = params[grepl("^enrichment_", names(params))])
        text <- as.character(format(record_numbers))
        expect_true(grepl(d3, text, fixed = TRUE), info = template)
        expect_true(grepl("3 annotated terms that this release does not contain were not used.", text, fixed = TRUE),
            info = template)
        expect_true(grepl(goPropD7, text, fixed = TRUE), info = template)
        # KIT-203: a fresh record carries the default rule, so the paragraph states the whole-family flag.
        expect_true(grepl(paste("A term was called significant when its Benjamini-Hochberg adjusted p-value over that",
            "family was below 0.05"), text, fixed = TRUE), info = template)
        expect_false(grepl("A term was called significant when its q-value (Storey, 2002)", text, fixed = TRUE),
            info = template)
        # KIT-205: the fresh record carries the background rule marker ("tested").
        expect_true(grepl(paste0("The background was the proteins tested in that comparison (those with a p-value) with ",
            "at least one Gene Ontology annotation (on a run with parent terms, at least one below the three top-level ",
            "terms)."), text, fixed = TRUE), info = template)
        expect_false(grepl("Benjamini-Hochberg method", text, fixed = TRUE), info = template)
        expect_false(grepl("all identified proteins", text, fixed = TRUE), info = template)
        # (b) a record without go_propagation: the D7 statements, no D3 sentence, no parent-term clause.
        plain <- as.character(format(goPropNumbers("clusterprofiler")))
        expect_true(grepl(goPropD7, plain, fixed = TRUE), info = template)
        expect_false(grepl("Gene Ontology protein sets were built", plain, fixed = TRUE), info = template)
        expect_false(grepl("on a run with parent terms", plain, fixed = TRUE), info = template)
        # KIT-203: a record without the rule keeps the q-value's paragraph and is told its measured failure.
        expect_true(grepl("A term was called significant when its q-value (Storey, 2002) was below 0.05", plain,
            fixed = TRUE), info = template)
        expect_true(grepl("Because terms without a list protein are left out before the adjustment", plain,
            fixed = TRUE), info = template)
        # (c) the count not recorded: no lost-term clause.
        unrecorded <- record_numbers
        unrecorded$enrichment_go_ids_not_in_go_db_total <- "not recorded"
        text_c <- as.character(format(unrecorded))
        expect_true(grepl(d3, text_c, fixed = TRUE), info = template)
        expect_false(grepl("annotated terms that this release does not contain", text_c, fixed = TRUE), info = template)
    }
})

test_that("KIT-173 T11(d): a g:Profiler run prints the corrected g:Profiler paragraph, the correction named from the record", {
    tool <- paste0("The enrichment analysis was performed with g:Profiler (Kolberg et al., 2023) through the gprofiler2 R ",
        "package (Kolberg et al., 2020). Gene Ontology biological process, molecular function and cellular component ",
        "terms, KEGG pathways and Reactome pathways were tested, each data source as its own family, ")
    cases <- list(
        gSCS = list(extra = list(enrichment_correction_method = "gSCS"), ending = "with g:SCS correction of the p-values."),
        fdr = list(extra = list(enrichment_correction_method = "fdr"),
            ending = "with Benjamini-Hochberg false discovery rate correction of the p-values."),
        none = list(extra = list(), ending = paste0("with the correction selected for the run (g:SCS unless the ",
            "Benjamini-Hochberg false discovery rate was chosen).")))
    for (template in goPropTemplates) {
        format <- goPropHeadFormatter(template)
        for (name in names(cases)) {
            text <- as.character(format(goPropNumbers("gprofiler2", extra = cases[[name]]$extra)))
            expect_true(grepl(paste0(tool, cases[[name]]$ending), text, fixed = TRUE), info = paste(template, name))
            expect_true(grepl("• A term was called significant when its adjusted p-value was below 0.05.", text,
                fixed = TRUE), info = paste(template, name))
            # KIT-205: without a record, the tested proteins; KIT-208's clause (E-7).
            expect_true(grepl(paste0("• The background was the identifiers of the proteins tested in each comparison ",
                "(those with a p-value in its differential-abundance table; the first gene name of each protein, from ",
                "the protein annotation or, where the annotation gives none, from the identifier table; proteins without ",
                "a gene name were left out of the lists and the background) that g:Profiler recognises (a custom ",
                "background)."), text, fixed = TRUE), info = paste(template, name))
            expect_false(grepl("clusterProfiler", text, fixed = TRUE), info = paste(template, name))
            expect_false(grepl("Benjamini-Hochberg method", text, fixed = TRUE), info = paste(template, name))
            expect_false(grepl("adjusted FDR cutoff", text, fixed = TRUE), info = paste(template, name))
            expect_false(grepl(goPropD7, text, fixed = TRUE), info = paste(template, name))
        }
    }
})

test_that("KIT-173 T11(e): N/A and an absent method print the base paragraph byte for byte (regression)", {
    for (template in goPropTemplates) {
        base <- goPropFormatter(goPropFixtureDir("report-base", sub("\\.rmd$", ".R", template)))
        head <- goPropHeadFormatter(template)
        for (numbers in list(goPropNumbers("N/A"), goPropNumbers(NULL), goPropNumbers("N/A", up = 1, down = 1))) {
            expect_identical(as.character(head(numbers)), as.character(base(numbers)), info = template)
        }
    }
})

test_that("KIT-173 T11(f): on both routes the fold-change bullet states the cutoffs the kernel applies", {
    for (template in goPropTemplates) {
        format <- goPropHeadFormatter(template)
        for (source in c("clusterprofiler", "gprofiler2")) {
            text <- as.character(format(goPropNumbers(source, up = 0, down = 1)))
            expect_true(grepl("above 0 and down-regulated proteins those below -1", text, fixed = TRUE),
                info = paste(template, source))
            expect_false(grepl("No fold change cutoff", text, fixed = TRUE), info = paste(template, source))
            both <- as.character(format(goPropNumbers(source, up = 1, down = 1)))
            expect_true(grepl("above 1 and down-regulated proteins those below -1", both, fixed = TRUE),
                info = paste(template, source))
            expect_false(grepl("< 1", both, fixed = TRUE), info = paste(template, source))
            none <- as.character(format(goPropNumbers(source)))
            expect_true(grepl("no fold-change cutoff was applied", none, fixed = TRUE), info = paste(template, source))
        }
    }
})

# --- T12 (VAL-61): geneID cells in a defined order ----------------------------------------------------------------

test_that("KIT-173 T12 (VAL-61): every geneID cell of the clusterProfiler route lists its proteins alphabetically", {
    goPropSkipUnlessGoDb()
    roles <- goPropReadTsv("t7-fixture-entries.tsv")
    direct <- goPropConvert(goPropReadTsv("t7-annotation-ecoli.tsv"))
    entries <- sort(direct$Entry)
    table <- data.frame(Protein.IDs = entries, fdr_qvalue = 0.01,
        log2FC = ifelse(entries %in% roles$entry[roles$role == "p_term"], 3,
            ifelse(entries %in% roles$entry[roles$role == "d_term"], -3, 0)), stringsAsFactors = FALSE)
    da <- methods::new("da_results_for_enrichment", contrasts = tibble::tibble(contrast = "t12"),
        da_data = list(t12 = table), design_matrix = data.frame(sample = "S1", stringsAsFactors = FALSE))
    results <- NULL
    utils::capture.output(results <- suppressMessages(suppressWarnings(processEnrichments(da, "999999",
        q_cutoff = 0.05, pathway_dir = withr::local_tempdir(), go_annotations = direct, exclude_iea = FALSE,
        protein_id_column = "Protein.IDs", contrast_names = "t12"))))
    cells <- unlist(lapply(c("up", "down"), function(direction) results@enrichment_data$t12[[direction]]@result$geneID))
    expect_gt(sum(grepl("/", cells, fixed = TRUE)), 10L)
    sorted <- vapply(strsplit(cells, "/", fixed = TRUE), function(g) identical(g, sort(g, method = "radix")), logical(1))
    expect_identical(sum(!sorted), 0L)
})
