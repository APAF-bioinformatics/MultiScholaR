# KIT-173 F1 (R0): feasibility values, fixture fetches, the T5 pilot and the term selections.
#
# Usage, from the package root on a run-tests runner, after the commit under test is installed:
#   Rscript docs/qa/workflow-execution/kit/kit-173/f1-feasibility.R <out dir>
#
# A throwaway script (ticket section 8): it never calls the propagation function under test. Its closures are
# GO.db's own ancestor tables (five relations) and a recursive SQL query on GO.db's SQLite (any relation set).
# Public data only: UniProt (CC BY 4.0), QuickGO, g:Profiler's GO GMT (CC BY 4.0), Europe PMC metadata.
# Stop rules (section 8): UniProt's go_id includes ancestors or NOT-qualified terms; a 2.2 value not reproduced;
# a fetch failing twice. Each is written to stop-rules.tsv; the decision is taken by reading it.

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[[1L]] else "k173-f1-out"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
fixture_dir <- file.path(out_dir, "fixtures")
dir.create(fixture_dir, showWarnings = FALSE, recursive = TRUE)
options(timeout = 1800, warn = 1)

SEED_TERMS <- 20261123L
ROOTS <- c(BP = "GO:0008150", MF = "GO:0003674", CC = "GO:0005575")
FIVE <- c("isa", "part of", "regulates", "positively regulates", "negatively regulates")
TWO <- c("isa", "part of")
UA <- "MultiScholaR-validation-runner (R httr)"

values <- list()
record <- function(measure, value) values[[measure]] <<- paste(as.character(value), collapse = ",")
stops <- data.frame(rule = character(), triggered = logical(), detail = character())
stopRule <- function(rule, triggered, detail) {
    stops[nrow(stops) + 1L, ] <<- list(rule, isTRUE(triggered), detail)
}
writeTsv <- function(x, name) {
    if (is.null(x)) x <- data.frame(none = character())
    utils::write.table(x, file.path(out_dir, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
writeAll <- function() {
    writeTsv(data.frame(measure = names(values), value = unlist(values, use.names = FALSE)), "f1-values.tsv")
    writeTsv(stops, "stop-rules.tsv")
}
sha256 <- function(path) sub("[[:space:]].*$", "", system2("sha256sum", shQuote(path), stdout = TRUE)[[1L]])
section <- function(text) message("\n=== ", text, " ===")

# One GET with one retry; the response, or NULL after two failures (recorded).
fetch <- function(name, url, query = list(), accept = NULL, path = NULL) {
    config <- list(httr::user_agent(UA), httr::timeout(900))
    if (!is.null(accept)) config <- c(config, list(httr::accept(accept)))
    if (!is.null(path)) config <- c(config, list(httr::write_disk(path, overwrite = TRUE)))
    for (attempt in 1:2) {
        response <- tryCatch(do.call(httr::GET, c(list(url = url, query = query), config)),
            error = function(e) e)
        ok <- !inherits(response, "error") && httr::status_code(response) == 200L
        record(paste0("fetch_", name, "_attempt", attempt),
            if (inherits(response, "error")) class(response)[[1L]] else httr::status_code(response))
        if (ok) return(response)
        Sys.sleep(10)
    }
    stopRule(paste("fetch failed twice:", name), TRUE, url)
    NULL
}

record("commit", tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)[[1L]], error = function(e) NA))
for (pkg in c("GO.db", "AnnotationDbi", "clusterProfiler", "enrichit", "gprofiler2", "MultiScholaR")) {
    record(paste0("version_", pkg), tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) "absent"))
}

# --- GO.db: info, labels, edges, closures (2.2) ---------------------------------------------------------------------
section("GO.db")
suppressPackageStartupMessages(library(GO.db))
info <- GO.db::GO_dbInfo()
writeTsv(info, "go-dbinfo.tsv")
con <- GO.db::GO_dbconn()

parent_maps <- list(BP = GO.db::GOBPPARENTS, MF = GO.db::GOMFPARENTS, CC = GO.db::GOCCPARENTS)
labels <- do.call(rbind, lapply(names(parent_maps), function(aspect) {
    edges <- AnnotationDbi::toTable(parent_maps[[aspect]])
    counts <- table(edges[[3L]])
    data.frame(aspect = aspect, label = names(counts), edges = as.integer(counts))
}))
writeTsv(labels, "go-parent-labels.tsv")
print(labels)
# as.list() names, the form a parent map's per-term vectors carry
list_names <- unique(unlist(lapply(parent_maps, function(m) unique(names(unlist(unname(as.list(m))))))))
record("parent_map_list_names", sort(list_names))
expected_edges <- data.frame(
    aspect = c("BP", "BP", "BP", "BP", "BP", "MF", "MF", "CC", "CC"),
    label = c("isa", "part of", "regulates", "positively regulates", "negatively regulates", "isa", "part of", "isa",
        "part of"),
    expected = c(41725L, 4493L, 2867L, 2506L, 2507L, 12413L, 8L, 4664L, 1824L))
edge_check <- merge(expected_edges, labels, all = TRUE)
stopRule("2.2 edge counts by label reproduced", !isTRUE(all(edge_check$edges == edge_check$expected)),
    paste(sprintf("%s %s %s/%s", edge_check$aspect, edge_check$label, edge_check$edges, edge_check$expected),
        collapse = "; "))

terms <- DBI::dbGetQuery(con, "SELECT ontology, count(*) AS n FROM go_term GROUP BY ontology")
writeTsv(terms, "go-term-counts.tsv")
record("terms_BP_MF_CC_all", terms$n[match(c("BP", "MF", "CC", "universal"), terms$ontology)])

closureSql <- function(aspect, relations) {
    table <- sprintf("go_%s_parents", tolower(aspect))
    rel <- paste(sprintf("'%s'", relations), collapse = ", ")
    sql <- sprintf(paste(
        "WITH RECURSIVE anc(c, a) AS (",
        "SELECT _id, _parent_id FROM %1$s WHERE relationship_type IN (%2$s)",
        "UNION SELECT anc.c, p._parent_id FROM anc JOIN %1$s p ON p._id = anc.a",
        "WHERE p.relationship_type IN (%2$s))",
        "SELECT t1.go_id AS term, t2.go_id AS ancestor FROM anc",
        "JOIN go_term t1 ON t1._id = anc.c JOIN go_term t2 ON t2._id = anc.a WHERE t2.go_id <> 'all'"),
        table, rel)
    DBI::dbGetQuery(con, sql)
}
ancestor_maps <- list(BP = GO.db::GOBPANCESTOR, MF = GO.db::GOMFANCESTOR, CC = GO.db::GOCCANCESTOR)
closure_five <- list()
closure_two <- list()
pairs <- data.frame(aspect = names(ROOTS), five_sql = NA_integer_, two_sql = NA_integer_, ancestor_map = NA_integer_,
    ancestor_map_equals_five = NA)
expected_pairs <- list(BP = c(464934L, 352532L), MF = c(64418L, 64418L), CC = c(42075L, 42075L))
for (i in seq_along(ROOTS)) {
    aspect <- names(ROOTS)[[i]]
    closure_five[[aspect]] <- closureSql(aspect, FIVE)
    closure_two[[aspect]] <- closureSql(aspect, TWO)
    anc <- AnnotationDbi::toTable(ancestor_maps[[aspect]])
    anc <- data.frame(term = anc[[1L]], ancestor = anc[[2L]])
    anc <- anc[anc$ancestor != "all", ]
    key <- function(d) sort(paste(d$term, d$ancestor))
    pairs$five_sql[[i]] <- nrow(closure_five[[aspect]])
    pairs$two_sql[[i]] <- nrow(closure_two[[aspect]])
    pairs$ancestor_map[[i]] <- nrow(anc)
    pairs$ancestor_map_equals_five[[i]] <- identical(key(anc), key(closure_five[[aspect]]))
}
writeTsv(pairs, "go-closure-pairs.tsv")
print(pairs)
pairs_ok <- all(pairs$five_sql == vapply(expected_pairs, `[[`, 1L, 1L)) &&
    all(pairs$two_sql == vapply(expected_pairs, `[[`, 1L, 2L)) && all(pairs$ancestor_map_equals_five)
stopRule("2.2 pair counts and GO*ANCESTOR = five-relation closure reproduced", !pairs_ok,
    paste(sprintf("%s five %s two %s map %s equal %s", pairs$aspect, pairs$five_sql, pairs$two_sql,
        pairs$ancestor_map, pairs$ancestor_map_equals_five), collapse = "; "))

cross <- do.call(rbind, lapply(names(ROOTS), function(aspect) {
    sql <- sprintf(paste("SELECT t.ontology AS parent_ontology, count(*) AS n FROM go_%s_parents p",
        "JOIN go_term t ON t._id = p._parent_id GROUP BY t.ontology"), tolower(aspect))
    cbind(aspect = aspect, DBI::dbGetQuery(con, sql))
}))
writeTsv(cross, "go-cross-aspect-edges.tsv")
cross_other <- cross[cross$parent_ontology != cross$aspect, ]
stopRule("2.2 no edge to another aspect except each root's edge to all",
    !(all(cross_other$parent_ontology == "universal") && all(cross_other$n == 1L)),
    paste(sprintf("%s->%s %s", cross$aspect, cross$parent_ontology, cross$n), collapse = "; "))

ex <- closure_two$BP$ancestor[closure_two$BP$term == "GO:0006110"]
ex5 <- closure_five$BP$ancestor[closure_five$BP$term == "GO:0006110"]
record("GO:0006110_ancestors_two_five", c(length(ex), length(ex5)))
record("GO:0006110_two_has_0006096_0008152", c("GO:0006096", "GO:0008152") %in% ex)
record("GO:0006110_five_has_0006096_0008152", c("GO:0006096", "GO:0008152") %in% ex5)
record("GO:0006096_ancestors_two_five", c(sum(closure_two$BP$term == "GO:0006096"),
    sum(closure_five$BP$term == "GO:0006096")))
stopRule("2.2 GO:0006110 18 vs 75 ancestors reproduced", !(length(ex) == 18L && length(ex5) == 75L),
    sprintf("%d / %d", length(ex), length(ex5)))
record("secondary_ids", DBI::dbGetQuery(con, "SELECT count(DISTINCT secondary) AS n FROM go_synonym
    WHERE secondary IS NOT NULL")$n)
record("obsolete_ids", DBI::dbGetQuery(con, "SELECT count(*) AS n FROM go_obsolete")$n)
record("GO:0000002_obsolete", DBI::dbGetQuery(con, "SELECT count(*) AS n FROM go_obsolete
    WHERE go_id = 'GO:0000002'")$n)
record("GO:0000021_secondary_of", DBI::dbGetQuery(con, "SELECT t.go_id FROM go_synonym s
    JOIN go_term t ON t._id = s._id WHERE s.secondary = 'GO:0000021'")$go_id)
writeAll()

# --- clusterProfiler::enricher on the lock ---------------------------------------------------------------------------
section("enricher")
enricher_text <- c(
    "## formals(clusterProfiler::enricher)", utils::capture.output(str(formals(clusterProfiler::enricher))),
    "## body(clusterProfiler::enricher)", deparse(body(clusterProfiler::enricher)))
cp_ns <- asNamespace("clusterProfiler")
for (fn in intersect(c("enricher_internal", ".enricher_internal"), ls(cp_ns, all.names = TRUE))) {
    enricher_text <- c(enricher_text, paste0("## clusterProfiler:::", fn), deparse(get(fn, cp_ns)))
}
if (requireNamespace("enrichit", quietly = TRUE)) {
    en_ns <- asNamespace("enrichit")
    enricher_text <- c(enricher_text, "## ls(asNamespace('enrichit'))", paste(ls(en_ns, all.names = TRUE),
        collapse = " "))
    for (fn in intersect(c("ora_gson", "ora", "calculate_qvalue", "enricher_internal"), ls(en_ns, all.names = TRUE))) {
        enricher_text <- c(enricher_text, paste0("## enrichit:::", fn), deparse(get(fn, en_ns)))
    }
}
writeLines(enricher_text, file.path(out_dir, "enricher-on-lock.txt"))
f <- formals(clusterProfiler::enricher)
record("enricher_minGSSize_maxGSSize", c(deparse(f$minGSSize), deparse(f$maxGSSize)))

# --- UniProt fixtures ------------------------------------------------------------------------------------------------
section("UniProt")
UNIPROT_FIELDS <- "accession,go_id,gene_oln,reviewed"
uniprot <- list()
for (organism in c(yeast = "559292", ecoli = "83333")) {
    name <- names(which(c(yeast = "559292", ecoli = "83333") == organism))
    query_text <- sprintf("(organism_id:%s) AND (reviewed:true)", organism)
    raw_path <- file.path(fixture_dir, sprintf("uniprot-%s-%s.tsv", name, organism))
    response <- fetch(paste0("uniprot_", name), "https://rest.uniprot.org/uniprotkb/stream",
        query = list(query = query_text, fields = UNIPROT_FIELDS, format = "tsv"), path = raw_path)
    if (is.null(response)) next
    headers <- httr::headers(response)
    record(paste0("uniprot_", name, "_query"), query_text)
    record(paste0("uniprot_", name, "_fields"), UNIPROT_FIELDS)
    record(paste0("uniprot_", name, "_release"), headers[["x-uniprot-release"]] %||% NA)
    record(paste0("uniprot_", name, "_release_date"), headers[["x-uniprot-release-date"]] %||% NA)
    record(paste0("uniprot_", name, "_fetched_utc"), format(Sys.time(), tz = "UTC", usetz = TRUE))
    gz_path <- paste0(raw_path, ".gz")
    zcon <- gzfile(gz_path, "w")
    writeLines(readLines(raw_path, warn = FALSE), zcon)
    close(zcon)
    tbl <- utils::read.delim(raw_path, sep = "\t", quote = "", stringsAsFactors = FALSE)
    record(paste0("uniprot_", name, "_columns"), names(tbl))
    record(paste0("uniprot_", name, "_rows"), nrow(tbl))
    record(paste0("uniprot_", name, "_sha256_tsv"), sha256(raw_path))
    record(paste0("uniprot_", name, "_sha256_gz"), sha256(gz_path))
    record(paste0("uniprot_", name, "_bytes_gz"), file.size(gz_path))
    unlink(raw_path)
    ids <- strsplit(ifelse(is.na(tbl$Gene.Ontology.IDs), "", tbl$Gene.Ontology.IDs), "; ", fixed = TRUE)
    root_only <- vapply(ids, function(x) length(x) > 0L && all(x %in% ROOTS), logical(1L))
    only_bp_root <- vapply(ids, function(x) identical(x, "GO:0008150"), logical(1L))
    record(paste0("uniprot_", name, "_entries_with_go"), sum(lengths(ids) > 0L))
    record(paste0("uniprot_", name, "_root_only_entries"), sum(root_only))
    record(paste0("uniprot_", name, "_only_GO:0008150_entries"), sum(only_bp_root))
    uniprot[[name]] <- tbl
}
writeAll()
gz_read <- function(name) uniprot[[name]]

# The route's conversion (unchanged package code, not the function under test)
convert <- function(tbl) {
    MultiScholaR::uniprotGoIdToTerm(tbl, uniprot_id_column = Entry, go_id_column = Gene.Ontology.IDs, sep = "; ")
}
directPairs <- function(conv) {
    cols <- c(BP = "go_id_go_biological_process", MF = "go_id_go_molecular_function",
        CC = "go_id_go_cellular_compartment")
    do.call(rbind, lapply(names(cols), function(aspect) {
        x <- conv[!is.na(conv[[cols[[aspect]]]]), c("Entry", cols[[aspect]])]
        if (nrow(x) == 0L) return(NULL)
        s <- strsplit(x[[2L]], "; ", fixed = TRUE)
        data.frame(Entry = rep(x$Entry, lengths(s)), TERM = unlist(s), aspect = aspect)
    }))
}
propagate <- function(direct, closures) {
    out <- lapply(names(ROOTS), function(aspect) {
        d <- direct[direct$aspect == aspect, c("Entry", "TERM")]
        cl <- closures[[aspect]]
        up <- merge(d, cl, by.x = "TERM", by.y = "term")
        p <- unique(rbind(d, data.frame(Entry = up$Entry, TERM = up$ancestor)))
        p <- p[p$TERM != ROOTS[[aspect]], ]
        if (nrow(p)) p$aspect <- aspect
        p
    })
    do.call(rbind, out)
}

# --- g:Profiler: version, GMT variants, GO GMT -------------------------------------------------------------------
section("g:Profiler")
version_info <- tryCatch(gprofiler2::get_version_info(organism = "scerevisiae"), error = function(e) conditionMessage(e))
jsonlite::write_json(version_info, file.path(out_dir, "gprofiler-version.json"), auto_unbox = TRUE, pretty = TRUE)
gmt_candidates <- c(
    "gprofiler_full_scerevisiae.ENSG.gmt", "gprofiler_full_scerevisiae.name.gmt",
    "gprofiler_full_scerevisiae_noiea.ENSG.gmt", "gprofiler_full_scerevisiae.ENSG.noiea.gmt",
    "gprofiler_full_scerevisiae_no_iea.ENSG.gmt", "gprofiler_scerevisiae.ENSG.gmt.zip")
variants <- data.frame(file = gmt_candidates, status = vapply(gmt_candidates, function(f) {
    r <- tryCatch(httr::HEAD(paste0("https://biit.cs.ut.ee/gprofiler/static/", f), httr::user_agent(UA),
        httr::timeout(60)), error = function(e) NULL)
    if (is.null(r)) NA_integer_ else httr::status_code(r)
}, integer(1L)))
writeTsv(variants, "gprofiler-gmt-variants.tsv")
print(variants)
for (page in c(docs = "https://biit.cs.ut.ee/gprofiler/page/docs", faq = "https://biit.cs.ut.ee/gprofiler/page/faq")) {
    name <- names(which(c(docs = "https://biit.cs.ut.ee/gprofiler/page/docs",
        faq = "https://biit.cs.ut.ee/gprofiler/page/faq") == page))
    r <- fetch(paste0("gprofiler_", name), page)
    if (!is.null(r)) writeLines(httr::content(r, "text", encoding = "UTF-8"),
        file.path(out_dir, sprintf("gprofiler-page-%s.html.txt", name)))
}
gmt_url <- "https://biit.cs.ut.ee/gprofiler/static/gprofiler_full_scerevisiae.ENSG.gmt"
gmt_full <- file.path(tempdir(), "gprofiler_full_scerevisiae.ENSG.gmt")
gmt_response <- fetch("gprofiler_gmt", gmt_url, path = gmt_full)
gmt <- NULL
if (!is.null(gmt_response)) {
    record("gmt_url", gmt_url)
    record("gmt_sha256_full", sha256(gmt_full))
    record("gmt_bytes_full", file.size(gmt_full))
    lines <- readLines(gmt_full, warn = FALSE)
    record("gmt_lines_full", length(lines))
    record("gmt_prefixes_full", names(table(sub(":.*$", "", sub("\t.*$", "", lines)))))
    go_lines <- lines[startsWith(lines, "GO:")]
    # R0b: the description field is replaced by the term id (the committed GMT keeps ids and genes only): a GO term
    # name such as "cargo receptor activity" matches the run's upload scan, which withheld the file in R0.
    go_lines <- sub("^([^\t]*)\t[^\t]*\t", "\\1\t\\1\t", go_lines)
    gmt_gz <- file.path(fixture_dir, "gprofiler-scerevisiae-go.ENSG.gmt.gz")
    zcon <- gzfile(gmt_gz, "w")
    writeLines(go_lines, zcon)
    close(zcon)
    record("gmt_go_lines", length(go_lines))
    record("gmt_sha256_go_gz", sha256(gmt_gz))
    record("gmt_bytes_go_gz", file.size(gmt_gz))
    parts <- strsplit(go_lines, "\t", fixed = TRUE)
    gmt <- data.frame(TERM = rep(vapply(parts, `[[`, "", 1L), lengths(parts) - 2L),
        GENE = unlist(lapply(parts, function(p) p[-(1:2)])))
}
writeAll()

# --- Europe PMC metadata of the two g:Profiler papers (F2) -------------------------------------------------------
papers <- c(service_2023 = "10.1093/nar/gkad347", package_2020 = "10.12688/f1000research.24956.2",
    service_2019 = "10.1093/nar/gkz369")
for (name in names(papers)) {
    r <- fetch(paste0("europepmc_", name), "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
        query = list(query = sprintf("DOI:\"%s\"", papers[[name]]), format = "json", resultType = "core"))
    if (is.null(r)) next
    text <- httr::content(r, "text", encoding = "UTF-8")
    writeLines(text, file.path(out_dir, sprintf("europepmc-%s.json", name)))
    # R0b: the open-access full text, for F2's check of the per-source correction statement
    pmcid <- tryCatch(jsonlite::fromJSON(text)$resultList$result$pmcid[[1L]], error = function(e) NA)
    if (is.character(pmcid) && !is.na(pmcid)) {
        full <- fetch(paste0("europepmc_full_", name), sprintf(
            "https://www.ebi.ac.uk/europepmc/webservices/rest/%s/fullTextXML", pmcid))
        if (!is.null(full)) writeLines(httr::content(full, "text", encoding = "UTF-8"),
            file.path(out_dir, sprintf("europepmc-%s-fulltext.xml.txt", name)))
    }
}

# --- QuickGO checks (stop rule: go_id holds ancestors or NOT-qualified terms) -------------------------------------
section("QuickGO")
quickgo <- function(query) {
    out <- list()
    page <- 1L
    repeat {
        r <- fetch(sprintf("quickgo_%s_p%d", gsub("[^A-Za-z0-9]", "", paste(query, collapse = "")), page),
            "https://www.ebi.ac.uk/QuickGO/services/annotation/search",
            query = c(query, list(limit = 200L, page = page)), accept = "application/json")
        if (is.null(r)) return(NULL)
        j <- jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"), simplifyVector = TRUE)
        if (length(j$results) == 0L) break
        out[[page]] <- data.frame(protein = sub("-.*$", "", sub("^UniProtKB:", "", j$results$geneProductId)),
            qualifier = j$results$qualifier, goId = j$results$goId)
        if (page >= j$pageInfo$total || page >= 25L) break
        page <- page + 1L
    }
    unique(do.call(rbind, out))
}
yeast <- gz_read("yeast")
qg_rows <- list()
if (!is.null(yeast)) {
    all_five <- do.call(rbind, closure_five)
    upGo <- function(acc) {
        x <- yeast$Gene.Ontology.IDs[yeast$Entry == acc]
        if (length(x) == 0L || is.na(x)) character() else strsplit(x, "; ", fixed = TRUE)[[1L]]
    }
    for (acc in c("P00330", "P00549", "P02829")) {
        q <- quickgo(list(geneProductId = acc))
        if (is.null(q)) next
        positive <- unique(q$goId[!startsWith(q$qualifier, "NOT")])
        u <- upGo(acc)
        anc_of_positive <- setdiff(unique(all_five$ancestor[all_five$term %in% positive]), positive)
        qg_rows[[acc]] <- data.frame(accession = acc, uniprot_go_ids = length(u), quickgo_positive = length(positive),
            both = length(intersect(u, positive)), uniprot_only = length(setdiff(u, positive)),
            quickgo_only = length(setdiff(positive, u)),
            uniprot_only_that_are_ancestors_of_quickgo = length(intersect(setdiff(u, positive), anc_of_positive)),
            ancestors_of_quickgo_terms = length(anc_of_positive),
            share_of_those_ancestors_in_uniprot = round(mean(anc_of_positive %in% u), 3))
    }
    not_rows <- list()
    for (qualifier in c("NOT|involved_in", "NOT|enables", "NOT|located_in")) {
        if (length(not_rows) == 3L) break
        nots <- quickgo(list(taxonId = "559292", geneProductType = "protein", qualifier = qualifier))
        if (is.null(nots)) next
        nots <- nots[startsWith(nots$qualifier, "NOT") & nots$protein %in% yeast$Entry, ]
        nots <- nots[order(nots$protein, nots$goId), ]
        for (k in seq_len(nrow(nots))) {
            if (length(not_rows) == 3L) break
            acc <- nots$protein[[k]]
            q <- quickgo(list(geneProductId = acc))
            if (is.null(q)) next
            positive <- unique(q$goId[!startsWith(q$qualifier, "NOT")])
            if (nots$goId[[k]] %in% positive) next
            not_rows[[length(not_rows) + 1L]] <- data.frame(accession = acc, not_term = nots$goId[[k]],
                qualifier = nots$qualifier[[k]], in_uniprot_go_id = nots$goId[[k]] %in% upGo(acc),
                ancestor_of_a_positive_term = nots$goId[[k]] %in% all_five$ancestor[all_five$term %in% positive])
        }
    }
    qg <- do.call(rbind, qg_rows)
    nt <- do.call(rbind, not_rows)
    writeTsv(qg, "quickgo-direct-vs-uniprot.tsv")
    writeTsv(nt, "quickgo-not-vs-uniprot.tsv")
    print(qg)
    print(nt)
    stopRule("go_id includes ancestors (UniProt holds most ancestors of QuickGO's direct terms)",
        is.null(qg) || any(qg$share_of_those_ancestors_in_uniprot > 0.5),
        if (is.null(qg)) "not measured" else paste(sprintf("%s %s", qg$accession, qg$share_of_those_ancestors_in_uniprot),
            collapse = "; "))
    stopRule("go_id includes NOT-qualified terms", is.null(nt) || nrow(nt) < 3L ||
        any(nt$in_uniprot_go_id & !nt$ancestor_of_a_positive_term),
        if (is.null(nt)) "not measured" else paste(sprintf("%s %s %s", nt$accession, nt$not_term, nt$in_uniprot_go_id),
            collapse = "; "))
}
writeAll()

# --- Pilot of T5 (50 terms) and mapping loss ------------------------------------------------------------------------
section("pilot")
if (!is.null(yeast) && !is.null(gmt)) {
    conv <- convert(yeast)
    direct <- directPairs(conv)
    oln <- strsplit(ifelse(is.na(yeast$Gene.Names..ordered.locus.), "", yeast$Gene.Names..ordered.locus.), "[ ;,/]+")
    entry_oln <- data.frame(Entry = rep(yeast$Entry, lengths(oln)), GENE = unlist(oln))
    entry_oln <- entry_oln[nzchar(entry_oln$GENE), ]
    record("yeast_entries_without_oln", sum(!yeast$Entry %in% entry_oln$Entry))
    gp_genes <- unique(gmt$GENE)
    record("gmt_go_genes", length(gp_genes))
    record("gmt_go_genes_unmapped_by_oln", sum(!gp_genes %in% entry_oln$GENE))
    record("gmt_go_genes_unmapped_share", round(mean(!gp_genes %in% entry_oln$GENE), 4))
    universe <- intersect(gp_genes, entry_oln$GENE)
    toGenes <- function(p) unique(merge(p, entry_oln, by = "Entry")[, c("TERM", "GENE")])
    sets <- list(
        none = toGenes(direct[direct$TERM != ROOTS[direct$aspect], c("Entry", "TERM")]),
        is_a_part_of = toGenes(propagate(direct, closure_two)),
        is_a_part_of_regulates = toGenes(propagate(direct, closure_five)))
    gp <- gmt[gmt$GENE %in% universe, ]
    gp_size <- table(gp$TERM)
    go_aspect <- DBI::dbGetQuery(con, "SELECT go_id, ontology FROM go_term")
    gmt_terms <- unique(gmt$TERM)
    record("gmt_go_terms", length(gmt_terms))
    record("gmt_go_terms_not_in_go_db", sum(!gmt_terms %in% go_aspect$go_id))
    in_scope <- names(gp_size)[gp_size >= 10L & gp_size <= 500L]
    in_scope <- intersect(in_scope, go_aspect$go_id)
    in_scope <- intersect(in_scope, unique(sets$is_a_part_of_regulates$TERM))
    record("pilot_terms_in_scope", length(in_scope))
    set.seed(SEED_TERMS)
    pilot_terms <- sort(sample(in_scope, min(50L, length(in_scope))))
    pilot <- do.call(rbind, lapply(names(sets), function(value) {
        s <- sets[[value]]
        s <- s[s$GENE %in% universe, ]
        do.call(rbind, lapply(pilot_terms, function(t) {
            ours <- s$GENE[s$TERM == t]
            ref <- gp$GENE[gp$TERM == t]
            both <- length(intersect(ours, ref))
            data.frame(value = value, term = t, ours = length(ours), reference = length(ref), both = both,
                C = both / length(ref), P = if (length(ours)) both / length(ours) else NA_real_)
        }))
    }))
    writeTsv(pilot, "pilot-t5-50-terms.tsv")
    summ <- do.call(rbind, lapply(split(pilot, pilot$value), function(d) data.frame(value = d$value[[1L]],
        terms = nrow(d), C_q1 = quantile(d$C, 0.25), C_median = median(d$C), C_q3 = quantile(d$C, 0.75),
        P_q1 = quantile(d$P, 0.25, na.rm = TRUE), P_median = median(d$P, na.rm = TRUE),
        P_q3 = quantile(d$P, 0.75, na.rm = TRUE), P_na = sum(is.na(d$P)))))
    writeTsv(summ, "pilot-t5-summary.tsv")
    print(summ)
}
writeAll()

# --- E. coli: T7 fixture selection and the 9.1 planted terms (rules fixed in the ticket) -------------------------
section("E. coli selections")
ecoli <- gz_read("ecoli")
if (!is.null(ecoli)) {
    conv <- convert(ecoli)
    direct <- directPairs(conv)
    prop <- propagate(direct, closure_five)
    bp <- prop[prop$aspect == "BP", ]
    size <- table(bp$TERM)
    # 9.1: one BP term per band, drawn with seed 20261123 from the sorted candidates (specific, mid, broad)
    set.seed(SEED_TERMS)
    bands <- list(specific = c(10L, 30L), mid = c(50L, 150L), broad = c(200L, 400L))
    planted <- do.call(rbind, lapply(names(bands), function(b) {
        cand <- sort(names(size)[size >= bands[[b]][[1L]] & size <= bands[[b]][[2L]]])
        pick <- cand[[sample.int(length(cand), 1L)]]
        data.frame(band = b, candidates = length(cand), term = pick, size = as.integer(size[[pick]]),
            name = AnnotationDbi::Term(GO.db::GOTERM[[pick]]))
    }))
    writeTsv(planted, "effect-planted-terms.tsv")
    print(planted)
    record("ecoli_annotated_entries", length(unique(direct$Entry)))
    record("ecoli_background_none", length(unique(direct$Entry)))
    record("ecoli_background_propagated", length(unique(prop$Entry)))

    # T7 fixture (rule of this script, stated in feasibility.md): P = a BP term with at most 3 direct proteins and
    # 15 to 40 propagated proteins; D = a term with 10 to 20 direct proteins; both drawn with seed 20261123 + 7
    # from sorted candidates. Fixture = P's propagated proteins, D's direct proteins, seeded fill from the proteins
    # with a non-root direct term to 69, plus one root-only protein. R0 found no reviewed entry whose GO ids are all
    # aspect roots (yeast 0, E. coli 0), so the root-only protein is synthetic: an entry without GO annotation, drawn
    # with seed 20261123 + 8, given GO:0008150 alone (as a curator's "no biological data available" annotation).
    direct_bp_size <- table(direct$TERM[direct$aspect == "BP"])
    direct_all_size <- table(direct$TERM)
    ids_raw <- strsplit(ifelse(is.na(ecoli$Gene.Ontology.IDs), "", ecoli$Gene.Ontology.IDs), "; ", fixed = TRUE)
    names(ids_raw) <- ecoli$Entry
    root_only <- names(ids_raw)[vapply(ids_raw, function(x) length(x) > 0L && all(x %in% ROOTS), logical(1L))]
    only_bp_root <- names(ids_raw)[vapply(ids_raw, function(x) identical(x, "GO:0008150"), logical(1L))]
    set.seed(SEED_TERMS + 7L)
    dsize <- direct_bp_size[names(size)]
    dsize[is.na(dsize)] <- 0L
    cand_p <- sort(names(size)[dsize <= 3L & size >= 15L & size <= 40L])
    p_term <- cand_p[[sample.int(length(cand_p), 1L)]]
    cand_d <- sort(setdiff(names(direct_all_size)[direct_all_size >= 10L & direct_all_size <= 20L], c(p_term, ROOTS)))
    d_term <- cand_d[[sample.int(length(cand_d), 1L)]]
    members <- union(unique(bp$Entry[bp$TERM == p_term]), unique(direct$Entry[direct$TERM == d_term]))
    members <- setdiff(members, root_only)
    has_nonroot <- unique(direct$Entry[!direct$TERM %in% ROOTS])
    members <- intersect(members, has_nonroot)
    pool <- sort(setdiff(has_nonroot, members))
    fill <- if (length(members) < 69L) sample(pool, 69L - length(members)) else character()
    no_go <- sort(ecoli$Entry[is.na(ecoli$Gene.Ontology.IDs) | !nzchar(ecoli$Gene.Ontology.IDs)])
    set.seed(SEED_TERMS + 8L)
    root_pick <- no_go[[sample.int(length(no_go), 1L)]]
    fixture_entries <- sort(c(members, fill, root_pick))
    t7 <- ecoli[ecoli$Entry %in% fixture_entries, ]
    t7$Gene.Ontology.IDs[t7$Entry == root_pick] <- "GO:0008150"
    ids_raw[[root_pick]] <- "GO:0008150"
    root_only <- c(root_only, root_pick)
    t7_path <- file.path(fixture_dir, "t7-annotation-ecoli.tsv")
    utils::write.table(t7, t7_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    # within-fixture checks
    fconv <- convert(t7)
    fdirect <- directPairs(fconv)
    fprop <- propagate(fdirect, closure_five)
    fprop_two <- propagate(fdirect, closure_two)
    record("t7_p_term", c(p_term, AnnotationDbi::Term(GO.db::GOTERM[[p_term]])))
    record("t7_d_term", c(d_term, AnnotationDbi::Term(GO.db::GOTERM[[d_term]])))
    record("t7_entries", length(fixture_entries))
    record("t7_root_only_entry", c(root_pick, paste(ids_raw[[root_pick]], collapse = " ")))
    record("t7_root_only_entries_in_fixture", sum(fixture_entries %in% root_only))
    record("t7_p_direct_propagated_five_two", c(sum(fdirect$TERM == p_term), sum(fprop$TERM == p_term),
        sum(fprop_two$TERM == p_term)))
    record("t7_d_direct", sum(fdirect$TERM == d_term))
    fsize <- table(fprop$TERM)
    fdsize <- table(fdirect$TERM)
    only_after <- names(fsize)[fsize >= 10L & !names(fsize) %in% names(fdsize)[fdsize >= 10L]]
    record("t7_terms_ge10_only_after_propagation", length(only_after))
    record("t7_terms_ge10_direct", sum(fdsize >= 10L))
    record("t7_terms_never_direct_ge10", sum(fsize >= 10L & !names(fsize) %in% fdirect$TERM))
    record("t7_sha256_tsv", sha256(t7_path))
    writeTsv(data.frame(entry = fixture_entries, role = ifelse(fixture_entries == root_pick, "root_only",
        ifelse(fixture_entries %in% bp$Entry[bp$TERM == p_term], "p_term", ifelse(fixture_entries %in%
            direct$Entry[direct$TERM == d_term], "d_term", "fill")))), "t7-fixture-entries.tsv")
}
writeAll()
print(stops)
message("F1 done")
