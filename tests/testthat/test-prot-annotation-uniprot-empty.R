# KIT-175: an empty or incomplete UniProt annotation becomes a visible, classed and recorded warning; enrichment refuses
# on no annotation; the request moves to httr2 with a time-out, three attempts and an outage stop; the query is
# field-qualified (D8, VAL-92). Tests T1 to T8 of the ticket (section 10). The HTTP client is mocked throughout
# (`httr2::local_mocked_responses`; T6 installs the httr double too) and the proxy points at a closed port, so a
# request that escapes a mock fails instead of reaching UniProt, on base and on head.

library(testthat)

k175Offline <- function(env = parent.frame()) {
    withr::local_envvar(http_proxy = "http://127.0.0.1:9", https_proxy = "http://127.0.0.1:9",
        HTTP_PROXY = "http://127.0.0.1:9", HTTPS_PROXY = "http://127.0.0.1:9", no_proxy = "", NO_PROXY = "",
        MULTISCHOLAR_TEST_MODE = "false", .local_envir = env)
    withr::local_options(multischolar.test_mode = FALSE, .local_envir = env)
}

k175Header <- paste("Entry", "Entry Name", "Protein names", "Gene Names", "Organism", "Length", "Gene Ontology IDs",
    "Reviewed", "Protein existence", "Annotation", sep = "\t")

#' Synthetic accessions of the UniProt pattern (P10001, P10002, ...)
k175Accessions <- function(n, from = 1L) sprintf("P%05d", 10000L + from - 1L + seq_len(n))

#' One TSV row per accession, as the search answers with the ticket's fields
k175Rows <- function(entries) {
    if (!length(entries)) return(character())
    paste(entries, paste0("NAME", seq_along(entries), "_HUMAN"), paste("Protein", entries), paste0("GENE", entries),
        "Homo sapiens (Human)", "100", "GO:0005575", "reviewed", "Evidence at protein level", "5.0", sep = "\t")
}

k175Body <- function(entries) paste(c(k175Header, k175Rows(entries)), collapse = "\n")

k175Response <- function(status = 200L, body = NULL, headers = list()) {
    httr2::response(status_code = status, url = "https://rest.uniprot.org/uniprotkb/search", method = "GET",
        headers = c(list(`Content-Type` = "text/plain; charset=utf-8", `x-uniprot-release` = "2026_04"), headers),
        body = if (is.null(body)) raw() else charToRaw(body))
}

#' The accessions a request asks about (either query form)
k175Requested <- function(req) {
    query <- httr2::url_parse(req$url)$query$query
    regmatches(query, gregexpr("[OPQ][0-9][A-Z0-9]{3}[0-9](-[0-9]+)?|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}(-[0-9]+)?",
        query))[[1L]]
}

#' Install an httr2 double: `handler(requested, i, req)` gives the i-th request's response (or raises). Returns the log.
k175Mock <- function(handler, env = parent.frame()) {
    log <- new.env(parent = emptyenv())
    log$calls <- 0L
    log$requests <- list()
    httr2::local_mocked_responses(function(req) {
        log$calls <- log$calls + 1L
        log$requests[[log$calls]] <- req
        handler(k175Requested(req), log$calls, req)
    }, env = env)
    log
}

#' The transport failure httr2 raises after a request's attempts (class `httr2_failure`, the curl error as parent)
k175TransportError <- function() {
    rlang::error_cnd(c("httr2_failure", "httr2_error"), message = "Failed to perform HTTP request.",
        parent = simpleError("Could not resolve host: rest.uniprot.org"))
}

#' Every answer holds the requested accessions it is given (`present`), with the headers of `headers(i)`
k175Answer <- function(absent = character(), headers = function(i, n) list(`X-Total-Results` = as.character(n))) {
    function(requested, i, req) {
        returned <- setdiff(unique(sub("-[0-9]+$", "", requested)), absent)
        k175Response(200L, k175Body(returned), headers(i, length(returned)))
    }
}

#' Run `expr` with console output, messages and the KIT-175 warnings captured; other warnings muffled
k175Capture <- function(expr) {
    conditions <- list()
    error <- NULL
    value <- NULL
    utils::capture.output(value <- tryCatch(withCallingHandlers(expr,
        warning = function(w) {
            if (inherits(w, c("prot_annotation_uniprot_empty", "prot_annotation_uniprot_incomplete"))) {
                conditions[[length(conditions) + 1L]] <<- w
            }
            invokeRestart("muffleWarning")
        },
        message = function(m) invokeRestart("muffleMessage")),
        error = function(e) {
            error <<- e
            NULL
        }))
    list(value = value, conditions = conditions, error = error,
        last = if (length(conditions)) conditions[[length(conditions)]] else NULL)
}

k175Manager <- function(label, env = parent.frame()) {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(.local_envir = env), label, "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    withr::defer(closeDIAStateAdapter(manager), envir = env)
    manager
}

k175UniprotRows <- function(manager) {
    Filter(function(row) identical(row$service, "uniprot"), manager$.context$snapshot$external_resources)
}

#' getUniprotAnnotations on `accessions` (25 per request, no delay) under the manager's scope
k175Annotate <- function(manager, accessions, cache_dir, taxon_id = 9606) {
    k175Capture(withDIAStateScope(manager, function() getUniprotAnnotations(
        data.frame(Protein.Ids = accessions, stringsAsFactors = FALSE), cache_dir = cache_dir, taxon_id = taxon_id,
        batch_size = 25, api_delay = 0)))
}

k175CacheFile <- function(cache_dir) file.path(cache_dir, "uniprot_annotations.RDS")

#' The design-step input: one protein group per row, saved as the owner's table
k175DesignTable <- function(groups) {
    text <- paste0("Run\tProtein.Group\tPrecursor.Normalised\n",
        paste0(rep(c("A", "B"), length.out = length(groups)), "\t", groups, "\t", seq_along(groups), "\n", collapse = ""))
    prepareDIAParserTable(readr::read_tsv(I(text), show_col_types = FALSE))
}

k175DesignPaths <- function(env = parent.frame()) {
    root <- withr::local_tempdir(.local_envir = env)
    paths <- list(results_dir = file.path(root, "results"), source_dir = file.path(root, "scripts"))
    dir.create(paths$results_dir)
    dir.create(paths$source_dir)
    paths
}

#' The design step with recording notification seams; returns the notices and removals
k175DesignStep <- function(manager, source, paths, env = parent.frame()) {
    log <- new.env(parent = emptyenv())
    log$shown <- list()
    log$removed <- character()
    testthat::local_mocked_bindings(
        showNotification = function(ui, action = NULL, duration = 5, closeButton = TRUE, id = NULL,
            type = "default", session = NULL) {
            log$shown[[length(log$shown) + 1L]] <- list(text = paste(as.character(ui), collapse = " "), type = type,
                duration = duration, id = id)
            invisible(id)
        },
        removeNotification = function(id, session = NULL) {
            log$removed <- c(log$removed, id)
            invisible(NULL)
        },
        .package = "shiny", .env = env)
    log$run <- function() {
        workflow_data <- new.env(parent = emptyenv())
        workflow_data$state_manager <- manager
        workflow_data$data_cln <- source
        workflow_data$column_mapping <- list(protein_col = "Protein.Group")
        workflow_data$taxon_id <- 9606
        workflow_data$tab_status <- list(design_matrix = "pending")
        session <- shiny::MockShinySession$new()
        captured <- k175Capture(withDIAStateScope(manager, function() shiny::withReactiveDomain(session,
            completeProtDesignPostCheckpoint(workflow_data, paths, session, successMessage = "saved"))))
        list(workflow_data = workflow_data, captured = captured)
    }
    log
}

#' The annotation notices of a log (the step's own "saved" message left out)
k175AnnotationNotices <- function(log) Filter(function(note) !identical(note$text, "saved"), log$shown)

#' No fixture accession appears in any of `texts`
k175ExpectNoAccession <- function(texts, accessions) {
    texts <- unlist(texts, use.names = FALSE)
    for (accession in accessions) expect_false(any(grepl(accession, texts, fixed = TRUE)), info = accession)
}

k175ConditionTexts <- function(condition) {
    if (is.null(condition)) return(character())
    fields <- unclass(condition)[setdiff(names(condition), c("trace", "parent", "call"))]
    c(conditionMessage(condition), unlist(lapply(fields, function(field) if (is.atomic(field)) as.character(field))))
}

# --- T1 request error ----------------------------------------------------------------------------------------------

test_that("KIT-175 T1: a request error is a classed empty warning with its reason and detail, recorded, not cached", {
    k175Offline()
    manager <- k175Manager("k175-t1")
    accessions <- k175Accessions(3L)
    source <- saveDIATableComponent(manager$.context, k175DesignTable(c("P10001;P10002", "P10003")))
    paths <- k175DesignPaths()
    log <- k175Mock(function(requested, i, req) stop(k175TransportError()))
    out <- k175Capture(withDIAStateScope(manager, function() retrieveProtDesignUniprotAnnotations(
        manager, source, "Protein.Group", 9606, resultsDir = paths$results_dir, sourceDir = paths$source_dir)))
    expect_null(out$error)
    expect_true(inherits(out$last, "prot_annotation_uniprot_empty"))
    expect_identical(out$last$reason, "request_failed")
    expect_identical(out$last$detail, "Could not resolve host: rest.uniprot.org")
    expect_identical(out$last$requested, 3L)
    expect_identical(nrow(out$value), 0L)
    expect_identical(attr(out$value, "annotation_outcome")$reason, "request_failed")
    expect_false(file.exists(file.path(paths$results_dir, "cache", "uniprot_annotations", "uniprot_annotations.RDS")))
    expect_false(file.exists(file.path(paths$source_dir, "uniprot_dat_cln.RDS")))
    rows <- k175UniprotRows(manager)
    expect_length(rows, 1L)
    expect_identical(rows[[1L]]$release_kind, "failed")
    expect_identical(rows[[1L]]$query_summary, "annotation failed: UniProt could not be reached; 0 of 3 accessions annotated")
    notice <- uniprotAnnotationNotice(attr(out$value, "annotation_outcome"))
    expect_identical(notice, paste("Protein annotation was not retrieved: UniProt could not be reached (Could not",
        "resolve host: rest.uniprot.org). The analysis can continue. Without annotation the results tables have no gene",
        "names, plots are labelled with accessions, and enrichment cannot be run. To try again, save the design again."))
    k175ExpectNoAccession(c(k175ConditionTexts(out$last), notice, rows[[1L]]$query_summary,
        unlist(attr(out$value, "annotation_outcome"))), accessions)
})

# --- T2 empty answers ------------------------------------------------------------------------------------------------

test_that("KIT-175 T2: every way of ending with no rows has its reason, and none is cached", {
    k175Offline()
    manager <- k175Manager("k175-t2")
    accessions <- k175Accessions(3L)
    cases <- list(
        list(name = "503", handler = function(requested, i, req) k175Response(503L), reason = "service_status",
            statuses = 503L),
        list(name = "429", handler = function(requested, i, req) k175Response(429L), reason = "service_status",
            statuses = 429L),
        list(name = "header only", handler = function(requested, i, req) k175Response(200L, k175Header,
            list(`X-Total-Results` = "0")), reason = "no_rows", statuses = integer()),
        list(name = "empty body", handler = function(requested, i, req) k175Response(200L, NULL),
            reason = "invalid_response", statuses = integer()),
        list(name = "html body", handler = function(requested, i, req) k175Response(200L,
            "<html><body>Sign in to the network</body></html>"), reason = "invalid_response", statuses = integer()))
    for (case in cases) {
        cache_dir <- withr::local_tempdir()
        k175Mock(case$handler)
        out <- k175Annotate(manager, accessions, cache_dir)
        expect_null(out$error, label = case$name)
        expect_true(inherits(out$last, "prot_annotation_uniprot_empty"), label = case$name)
        expect_identical(out$last$reason, case$reason, label = case$name)
        expect_identical(as.integer(out$last$statuses), case$statuses, label = case$name)
        expect_identical(nrow(out$value), 0L, label = case$name)
        expect_false(file.exists(k175CacheFile(cache_dir)), label = case$name)
    }

    # Identifiers that are not UniProt accessions: no request.
    cache_dir <- withr::local_tempdir()
    log <- k175Mock(k175Answer())
    out <- k175Annotate(manager, c("not-an-accession", "contaminant_1"), cache_dir)
    expect_identical(out$last$reason, "no_uniprot_accessions")
    expect_identical(out$last$requested, 2L)
    expect_identical(log$calls, 0L)
    expect_false(file.exists(k175CacheFile(cache_dir)))

    # Five requests, the third refused and the others header-only: service_status, one failed.
    cache_dir <- withr::local_tempdir()
    k175Mock(function(requested, i, req) if (i == 3L) k175Response(503L) else k175Response(200L, k175Header))
    out <- k175Annotate(manager, k175Accessions(125L), cache_dir)
    expect_identical(out$last$reason, "service_status")
    expect_identical(out$last$requests_failed, 1L)
    expect_identical(out$last$requests_total, 5L)
    expect_false(file.exists(k175CacheFile(cache_dir)))

    # Seven requests all refused: the outage stop after three; four not sent.
    cache_dir <- withr::local_tempdir()
    log <- k175Mock(function(requested, i, req) k175Response(503L))
    out <- k175Annotate(manager, k175Accessions(175L), cache_dir)
    expect_identical(out$last$reason, "service_status")
    expect_identical(log$calls, 3L)
    expect_identical(out$last$requests_failed, 7L)
    expect_identical(out$last$requests_not_sent, 4L)
    expect_false(file.exists(k175CacheFile(cache_dir)))

    # A second call with a success returns rows (nothing empty was cached).
    k175Mock(k175Answer())
    out <- k175Annotate(manager, k175Accessions(175L), cache_dir)
    expect_identical(nrow(out$value), 175L)
    expect_true(file.exists(k175CacheFile(cache_dir)))

    # A zero-row cache file from an earlier version: not returned; replaced on success, removed on failure.
    cache_dir <- withr::local_tempdir()
    saveRDS(createEmptyUniprotTable(), k175CacheFile(cache_dir))
    log <- k175Mock(k175Answer())
    out <- k175Annotate(manager, accessions, cache_dir)
    expect_identical(log$calls, 1L)
    expect_identical(nrow(out$value), 3L)
    expect_identical(nrow(readRDS(k175CacheFile(cache_dir))), 3L)
    cache_dir <- withr::local_tempdir()
    saveRDS(createEmptyUniprotTable(), k175CacheFile(cache_dir))
    k175Mock(function(requested, i, req) k175Response(503L))
    out <- k175Annotate(manager, accessions, cache_dir)
    expect_true(inherits(out$last, "prot_annotation_uniprot_empty"))
    expect_false(file.exists(k175CacheFile(cache_dir)))

    # An unreadable cache file: treated as no cache, removed, the query runs.
    cache_dir <- withr::local_tempdir()
    writeBin(charToRaw("not an RDS file"), k175CacheFile(cache_dir))
    log <- k175Mock(function(requested, i, req) k175Response(503L))
    out <- k175Annotate(manager, accessions, cache_dir)
    expect_null(out$error)
    expect_identical(log$calls, 1L)
    expect_false(file.exists(k175CacheFile(cache_dir)))
})

test_that("KIT-175 T2: Ensembl identifiers whose conversion fails give id_conversion_failed, recorded after not_queried", {
    k175Offline()
    manager <- k175Manager("k175-t2-ensembl")
    testthat::local_mocked_bindings(convertEnsemblToUniprot = function(...) stop("the conversion service is closed"),
        .package = "MultiScholaR")
    log <- k175Mock(k175Answer())
    data_tbl <- data.frame(Protein.Group = c("ENSP00000269305", "ENSP00000350283", "ENSP00000275493"),
        stringsAsFactors = FALSE)
    out <- k175Capture(withDIAStateScope(manager, function() getUniprotAnnotationsFull(data_tbl, "Protein.Group",
        withr::local_tempdir(), taxon_id = 9606)))
    expect_identical(log$calls, 0L)
    reasons <- vapply(out$conditions, function(condition) condition$reason, "")
    expect_identical(reasons, "id_conversion_failed")
    expect_identical(out$last$detail, "the conversion service is closed")
    expect_identical(nrow(out$value), 0L)
    rows <- k175UniprotRows(manager)
    expect_identical(vapply(rows, function(row) row$release_kind, ""), c("not_queried", "not_queried"))
    expect_identical(rows[[2L]]$query_summary,
        "not queried: the Ensembl identifiers could not be converted to UniProt accessions")
})

# --- T2b incomplete answers and the truncation signal ---------------------------------------------------------------

test_that("KIT-175 T2b: a failed request among answered ones is incomplete, flagged and not cached", {
    k175Offline()
    manager <- k175Manager("k175-t2b")
    accessions <- k175Accessions(125L)
    cache_dir <- withr::local_tempdir()
    answer <- k175Answer()
    k175Mock(function(requested, i, req) if (i == 3L) k175Response(503L) else answer(requested, i, req))
    out <- k175Annotate(manager, accessions, cache_dir)
    expect_true(inherits(out$last, "prot_annotation_uniprot_incomplete"))
    expect_identical(out$last$requests_failed, 1L)
    expect_identical(out$last$accessions_not_looked_up, 25L)
    expect_identical(out$last$accessions_requested, 125L)
    expect_identical(nrow(out$value), 100L)
    expect_false(file.exists(k175CacheFile(cache_dir)))
    rows <- k175UniprotRows(manager)
    expect_identical(rows[[length(rows)]]$release_kind, "header")
    expect_identical(rows[[length(rows)]]$query_summary,
        "1 of 5 requests failed (status 503); 25 of 125 accessions not looked up")
    k175Mock(k175Answer())
    out <- k175Annotate(manager, accessions, cache_dir)
    expect_identical(nrow(out$value), 125L)
    expect_length(out$conditions, 0L)
    expect_true(file.exists(k175CacheFile(cache_dir)))

    # Through the design step's core: the table carries its outcome as a plain list.
    source <- saveDIATableComponent(manager$.context, k175DesignTable(accessions))
    paths <- k175DesignPaths()
    k175Mock(function(requested, i, req) if (i == 3L) k175Response(503L) else answer(requested, i, req))
    out <- k175Capture(withDIAStateScope(manager, function() retrieveProtDesignUniprotAnnotations(
        manager, source, "Protein.Group", 9606, resultsDir = paths$results_dir, sourceDir = paths$source_dir)))
    outcome <- attr(out$value, "annotation_outcome")
    expect_type(outcome, "list")
    expect_identical(outcome$class, "prot_annotation_uniprot_incomplete")
    expect_identical(outcome$accessions_not_looked_up, 25L)
    expect_identical(attr(readRDS(file.path(paths$source_dir, "uniprot_dat_cln.RDS")), "annotation_outcome"), outcome)
})

test_that("KIT-175 T2b (i)-(v): an answer is cut short only when a requested accession is missing and a next page exists", {
    k175Offline()
    manager <- k175Manager("k175-t2b-signal")
    accessions <- k175Accessions(125L)
    next_link <- '<https://rest.uniprot.org/uniprotkb/search?cursor=abc&size=25>; rel="next"'

    # (i) F3's real shape: every requested accession, more hits than rows, a next link: complete.
    cache_dir <- withr::local_tempdir()
    log <- k175Mock(k175Answer(headers = function(i, n) list(`X-Total-Results` = "27", Link = next_link)))
    out <- k175Annotate(manager, accessions, cache_dir)
    expect_length(out$conditions, 0L)
    expect_identical(nrow(out$value), 125L)
    expect_true(file.exists(k175CacheFile(cache_dir)))
    rows <- k175UniprotRows(manager)
    expect_identical(rows[[length(rows)]]$release_kind, "header")
    expect_identical(rows[[length(rows)]]$query_summary, "uniprotkb search, organism_id 9606")
    # D8 (VAL-92): the query is field-qualified.
    query <- httr2::url_parse(log$requests[[1L]]$url)$query$query
    expect_identical(query, paste0("(", paste0("accession:", accessions[1:25], collapse = " OR "),
        ") AND organism_id:9606"))

    # (ii) One requested accession missing and a next link: incomplete, the other 24 rows kept.
    cache_dir <- withr::local_tempdir()
    missing <- accessions[[60L]]
    k175Mock(k175Answer(absent = missing, headers = function(i, n)
        if (i == 3L) list(`X-Total-Results` = "26", Link = next_link) else list(`X-Total-Results` = as.character(n))))
    out <- k175Annotate(manager, accessions, cache_dir)
    expect_true(inherits(out$last, "prot_annotation_uniprot_incomplete"))
    expect_identical(out$last$requests_truncated, 1L)
    expect_identical(out$last$requests_failed, 0L)
    expect_identical(out$last$accessions_not_looked_up, 1L)
    expect_identical(nrow(out$value), 124L)
    expect_false(file.exists(k175CacheFile(cache_dir)))

    # (iii) One requested accession missing, no next link: absent from UniProt, complete.
    cache_dir <- withr::local_tempdir()
    k175Mock(k175Answer(absent = missing))
    out <- k175Annotate(manager, accessions, cache_dir)
    expect_length(out$conditions, 0L)
    expect_identical(nrow(out$value), 124L)
    expect_true(file.exists(k175CacheFile(cache_dir)))

    # (iv) An isoform requested: the query asks for its canonical accession (no isoform term, which would make
    # UniProt return isoform entries for the whole request), and the canonical entry annotates it.
    cache_dir <- withr::local_tempdir()
    isoforms <- accessions
    isoforms[[60L]] <- paste0(accessions[[60L]], "-2")
    log <- k175Mock(k175Answer())
    out <- k175Annotate(manager, isoforms, cache_dir)
    queries <- vapply(log$requests, function(req) httr2::url_parse(req$url)$query$query, "")
    expect_false(any(grepl("accession:[A-Z0-9]+-[0-9]+", queries)))
    expect_true(any(grepl(paste0("accession:", accessions[[60L]], "\\b"), queries)))
    expect_length(out$conditions, 0L)
    expect_true(accessions[[60L]] %in% out$value$Entry)
    expect_true(file.exists(k175CacheFile(cache_dir)))

    # (v) A realistic absence under the field-qualified query: 2 of 25 absent, the total equal to the rows, no link.
    cache_dir <- withr::local_tempdir()
    k175Mock(k175Answer(absent = accessions[c(3L, 17L)]))
    out <- k175Annotate(manager, accessions[1:25], cache_dir)
    expect_length(out$conditions, 0L)
    expect_identical(nrow(out$value), 23L)
    expect_true(file.exists(k175CacheFile(cache_dir)))
    rows <- k175UniprotRows(manager)
    expect_identical(rows[[length(rows)]]$release_kind, "header")

    # Seven requests, two with rows then three refused: the outage stop; incomplete.
    cache_dir <- withr::local_tempdir()
    answer <- k175Answer()
    log <- k175Mock(function(requested, i, req) if (i <= 2L) answer(requested, i, req) else k175Response(503L))
    out <- k175Annotate(manager, k175Accessions(175L), cache_dir)
    expect_true(inherits(out$last, "prot_annotation_uniprot_incomplete"))
    expect_identical(log$calls, 5L)
    expect_identical(out$last$requests_failed, 5L)
    expect_identical(out$last$requests_not_sent, 2L)
    expect_identical(nrow(out$value), 50L)
    expect_false(file.exists(k175CacheFile(cache_dir)))
})

# --- T3 record -------------------------------------------------------------------------------------------------------

test_that("KIT-175 T3: the record holds one row per outcome, written once the outcome is known", {
    k175Offline()
    accessions <- k175Accessions(3L)
    last_row <- function(manager) {
        rows <- k175UniprotRows(manager)
        rows[[length(rows)]]
    }

    # No entry found: the release from the header.
    manager <- k175Manager("k175-t3-norows")
    k175Mock(function(requested, i, req) k175Response(200L, k175Header))
    k175Annotate(manager, accessions, withr::local_tempdir())
    expect_identical(last_row(manager)[c("release", "release_kind", "query_summary")],
        list(release = "2026_04", release_kind = "header",
            query_summary = "no entry found for 3 accessions under organism 9606"))

    # Refused: failed, release NA.
    manager <- k175Manager("k175-t3-refused")
    k175Mock(function(requested, i, req) k175Response(503L))
    k175Annotate(manager, accessions, withr::local_tempdir())
    expect_identical(last_row(manager)[c("release", "release_kind", "query_summary")],
        list(release = NA_character_, release_kind = "failed", query_summary = paste("annotation failed: UniProt",
            "refused 1 of 1 requests (status 503) and returned no entry for the others; 0 of 3 accessions annotated")))

    # Not queried.
    manager <- k175Manager("k175-t3-notqueried")
    k175Annotate(manager, c("not-an-accession"), withr::local_tempdir())
    expect_identical(last_row(manager)[c("release", "release_kind", "query_summary")],
        list(release = NA_character_, release_kind = "not_queried",
            query_summary = "not queried: none of 1 identifiers is a UniProt accession"))

    # First request answered, the second a transport error: empty, one failed row and no header row.
    manager <- k175Manager("k175-t3-midrun")
    answer <- k175Answer()
    k175Mock(function(requested, i, req) if (i == 1L) answer(requested, i, req) else stop(k175TransportError()))
    out <- k175Capture(withDIAStateScope(manager, function() getUniprotAnnotationsFull(
        data.frame(Protein.Group = k175Accessions(30L), stringsAsFactors = FALSE), "Protein.Group",
        withr::local_tempdir(), taxon_id = 9606)))
    expect_identical(nrow(out$value), 0L)
    expect_identical(out$last$reason, "request_failed")
    rows <- k175UniprotRows(manager)
    expect_identical(vapply(rows, function(row) row$release_kind, ""), "failed")

    # Not attempted (no protein column): release NA, not_queried, the reason in the summary (implementation audit I04).
    manager <- k175Manager("k175-t3-notattempted")
    source <- saveDIATableComponent(manager$.context, k175DesignTable("P10001"))
    out <- k175Capture(withDIAStateScope(manager, function() retrieveProtDesignUniprotAnnotations(
        manager, source, NULL, 9606)))
    expect_identical(out$error$reason, "not_attempted")
    expect_identical(last_row(manager)[c("release", "release_kind", "query_summary")],
        list(release = NA_character_, release_kind = "not_queried",
            query_summary = "not queried: annotation was not attempted (no protein identifier column is selected)"))

    # GO processing fails after a three-row answer: a header row, then a failed row.
    manager <- k175Manager("k175-t3-processing")
    k175Mock(k175Answer())
    testthat::local_mocked_bindings(uniprotGoIdToTerm = function(...) stop("GO terms could not be read"),
        .package = "MultiScholaR")
    out <- k175Capture(withDIAStateScope(manager, function() getUniprotAnnotationsFull(
        data.frame(Protein.Group = accessions, stringsAsFactors = FALSE), "Protein.Group", withr::local_tempdir(),
        taxon_id = 9606)))
    expect_identical(nrow(out$value), 0L)
    expect_identical(out$last$reason, "processing_failed")
    rows <- k175UniprotRows(manager)
    expect_identical(vapply(rows, function(row) row$release_kind, ""), c("header", "failed"))
    expect_identical(rows[[2L]]$query_summary, "annotation failed after retrieval: the answer could not be processed")
})

test_that("KIT-175 T3: the merge-back, the raw-results write and the scripts copy failing each record once", {
    k175Offline()
    accessions <- k175Accessions(3L)
    kinds <- function(manager) vapply(k175UniprotRows(manager), function(row) row$release_kind, "")

    # The Ensembl merge-back fails (the mapping lacks its id column): header, then failed.
    manager <- k175Manager("k175-t3-mergeback")
    k175Mock(k175Answer())
    testthat::local_mocked_bindings(convertEnsemblToUniprot = function(ids, ...) data.frame(converted_id = accessions,
        conversion_status = "success", stringsAsFactors = FALSE), .package = "MultiScholaR")
    out <- k175Capture(withDIAStateScope(manager, function() getUniprotAnnotationsFull(
        data.frame(Protein.Group = c("ENSP00000269305", "ENSP00000350283", "ENSP00000275493"),
            stringsAsFactors = FALSE), "Protein.Group", withr::local_tempdir(), taxon_id = 9606)))
    expect_identical(out$last$reason, "processing_failed")
    expect_identical(kinds(manager), c("header", "failed"))

    # The raw-results write fails: one failed row and no header row.
    manager <- k175Manager("k175-t3-rawwrite")
    cache_dir <- withr::local_tempdir()
    dir.create(file.path(cache_dir, "uniprot_results.tsv"))
    k175Mock(k175Answer())
    out <- k175Capture(withDIAStateScope(manager, function() getUniprotAnnotationsFull(
        data.frame(Protein.Group = accessions, stringsAsFactors = FALSE), "Protein.Group", cache_dir, taxon_id = 9606)))
    expect_identical(out$last$reason, "processing_failed")
    expect_identical(kinds(manager), "failed")

    # The scripts-directory copy fails: header then failed, the error re-raised with its reason; the notice says so.
    manager <- k175Manager("k175-t3-copy")
    source <- saveDIATableComponent(manager$.context, k175DesignTable(c("P10001;P10002", "P10003")))
    paths <- k175DesignPaths()
    paths$source_dir <- file.path(paths$source_dir, "missing", "folder")
    k175Mock(k175Answer())
    out <- k175Capture(withDIAStateScope(manager, function() retrieveProtDesignUniprotAnnotations(
        manager, source, "Protein.Group", 9606, resultsDir = paths$results_dir, sourceDir = paths$source_dir)))
    expect_true(inherits(out$error, "prot_annotation_uniprot_recorded"))
    expect_identical(out$error$reason, "processing_failed")
    expect_identical(kinds(manager), c("header", "failed"))
    expect_match(k175UniprotRows(manager)[[2L]]$query_summary, "^annotation retrieved but not saved to the project file")
    log <- k175DesignStep(manager, source, paths)
    unlink(file.path(paths$results_dir, "cache"), recursive = TRUE)
    log$run()
    notices <- k175AnnotationNotices(log)
    expect_match(notices[[length(notices)]]$text, "^Protein annotation was retrieved but could not be saved to the project")
})

# --- T4 app notice ---------------------------------------------------------------------------------------------------

test_that("KIT-175 T4: the design step's notice names the reason, stays until closed, and is replaced or removed", {
    k175Offline()
    manager <- k175Manager("k175-t4")
    source <- saveDIATableComponent(manager$.context, k175DesignTable(c("P10001;P10002", "P10003")))
    paths <- k175DesignPaths()
    log <- k175DesignStep(manager, source, paths)

    k175Mock(function(requested, i, req) k175Response(503L))
    log$run()
    log$run()
    notices <- k175AnnotationNotices(log)
    expect_length(notices, 2L)
    for (note in notices) {
        expect_identical(note$type, "warning")
        expect_null(note$duration)
        expect_identical(note$id, "prot_uniprot_annotation")
    }
    expect_identical(notices[[1L]]$text, paste("Protein annotation was not retrieved: UniProt answered with an error",
        "(status 503) for 1 of 1 requests. The analysis can continue. Without annotation the results tables have no gene",
        "names, plots are labelled with accessions, and enrichment cannot be run. To try again, save the design again."))

    k175Mock(k175Answer())
    log$run()
    expect_true("prot_uniprot_annotation" %in% log$removed)
    expect_identical(log$shown[[length(log$shown) - 1L]]$text, "UniProt annotations retrieved successfully.")

    # The error branch: no protein column mapped.
    run_without_column <- function() {
        workflow_data <- new.env(parent = emptyenv())
        workflow_data$state_manager <- manager
        workflow_data$data_cln <- source
        workflow_data$column_mapping <- list()
        workflow_data$taxon_id <- 9606
        workflow_data$tab_status <- list(design_matrix = "pending")
        session <- shiny::MockShinySession$new()
        k175Capture(shiny::withReactiveDomain(session,
            completeProtDesignPostCheckpoint(workflow_data, paths, session, successMessage = "saved")))
    }
    run_without_column()
    notices <- k175AnnotationNotices(log)
    last <- notices[[length(notices)]]
    expect_identical(last$text, paste("Protein annotation was not attempted: no protein identifier column is selected.",
        "The analysis can continue. Without annotation the results tables have no gene names, plots are labelled with",
        "accessions, and enrichment cannot be run."))
    expect_identical(last$id, "prot_uniprot_annotation")
    expect_null(last$duration)

    # An incomplete result shows its counts, not "retrieved successfully".
    unlink(file.path(paths$results_dir, "cache"), recursive = TRUE)
    source_many <- saveDIATableComponent(manager$.context, k175DesignTable(k175Accessions(50L)))
    log_many <- k175DesignStep(manager, source_many, paths)
    answer <- k175Answer()
    k175Mock(function(requested, i, req) if (i == 2L) k175Response(503L) else answer(requested, i, req))
    log_many$run()
    notices <- k175AnnotationNotices(log_many)
    expect_false(any(vapply(notices, function(note) grepl("retrieved successfully", note$text, fixed = TRUE), NA)))
    expect_identical(notices[[length(notices)]]$text, paste("Protein annotation is incomplete: 1 of 2 requests to",
        "UniProt failed (status 503), so 25 of 50 proteins were not looked up. The analysis can continue. The proteins",
        "not looked up have no gene names and are left out of enrichment. To try again, save the design again."))
})

# --- T5 downstream ---------------------------------------------------------------------------------------------------

#' The T5/T7 project's full annotation table (KIT-173's E. coli fixture through `getUniprotAnnotations`, no service)
k175ProjectAnnotations <- function(cache_dir) goPropAnnotate(goPropReadTsv("t7-annotation-ecoli.tsv"), cache_dir)

#' One proteomics project on the clusterProfiler route (goPropRunProject's, test-go-propagation.R's helper), with the
#' annotation's places set by the caller: `annotations` is written to `uniprot_dat_cln.RDS` (NULL: no file; a string:
#' those bytes), `bus` is put on `workflow$uniprot_dat_cln`, `da_annotations` publishes the DA, `cache` is written as
#' the project's UniProt cache. `before(manager)` runs inside the owner's scope before the DA. The run's error, the
#' notices and the step's outcomes are returned.
k175Project <- function(root, annotations = NULL, bus = NULL, da_annotations = NULL, cache = NULL, before = NULL,
    process_fn = processEnrichments, taxon = "999999", organism_filter = FALSE) {
    roles <- goPropReadTsv("t7-fixture-entries.tsv")
    paths <- goPropPaths(root)
    for (p in paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
    annotation_file <- file.path(paths$source_dir, "uniprot_dat_cln.RDS")
    if (is.character(annotations)) writeBin(charToRaw(annotations), annotation_file)
    else if (!is.null(annotations)) saveRDS(annotations, annotation_file)
    if (!is.null(cache)) {
        dir.create(file.path(paths$results_dir, "cache", "uniprot_annotations"), recursive = TRUE, showWarnings = FALSE)
        saveRDS(cache, file.path(paths$results_dir, "cache", "uniprot_annotations", "uniprot_annotations.RDS"))
    }
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "k175-project", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("DIA")
    entries <- sort(roles$entry)
    design <- data.frame(Run = paste0("S", 1:6), group = rep(c("A", "B"), each = 3L),
        batch = rep(c("b1", "b2", "b3"), times = 2L), stringsAsFactors = FALSE)
    quant <- goPropQuant(entries, roles$entry[roles$role == "p_term"], roles$entry[roles$role == "d_term"])
    object <- methods::new("ProteinQuantitativeData", protein_quant_table = quant, design_matrix = design,
        sample_id = "Run", group_id = "group", protein_id_column = "uniprot_acc",
        protein_id_table = data.frame(uniprot_acc = entries, gene_names = paste0("GENE", seq_along(entries)),
            stringsAsFactors = FALSE),
        args = list(globalParameters = list(workflow_type = "DIA", use_limpa = FALSE, report_template = "DIA_report"),
            deAnalysisParameters = list(formula_string = "~ 0 + group")))
    setDIAStepProvenance(manager, "design_builder")
    manager$saveState("correlation_filtered", object, list(), "fixture")
    manager$.context$experiment_paths <- paths
    contrasts <- data.frame(contrasts = "groupB-groupA", full_format = "B_vs_A=groupB-groupA",
        friendly_names = "B_vs_A", stringsAsFactors = FALSE)
    workflow <- list2env(list(state_manager = manager, contrasts_tbl = contrasts, tab_status = list(),
        qc_params = list(protein_qc = list()), config_list = list()))
    if (!is.null(bus)) workflow$uniprot_dat_cln <- bus
    parameters <- list(selected_contrast = "B_vs_A", organism_taxid = taxon, up_cutoff = 0, down_cutoff = 0,
        q_cutoff = 0.05, correction_method = "gSCS", enable_organism_filter = organism_filter)
    notes <- list()
    kernel_calls <- 0L
    kernel <- function(...) {
        kernel_calls <<- kernel_calls + 1L
        process_fn(...)
    }
    run <- NULL
    error <- NULL
    captured <- k175Capture(withDIAStateScope(manager, function() {
        if (!is.null(before)) before(manager)
        da_data <- new.env(parent = emptyenv())
        da_data$current_s4_object <- getDIAStateReference(manager)
        da_data$contrasts_available <- contrasts$contrasts
        manifest <- runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05,
            treat_lfc_cutoff = 0), da_data, workflow, paths, contrasts = contrasts, publish = FALSE)
        produceDIADaPublication(manager, manifest, da_annotations %||% annotations_or_empty(annotations), 0.05,
            paths[c("da_output_dir", "publication_graphs_dir")], save_assets = FALSE)
        run <<- runWorkflowEnrichment(workflow, omicGuiConfig("proteomics")$enrichment, manifest, parameters, paths,
            process_fn = kernel, notify = function(message, type, duration) {
                notes[[length(notes) + 1L]] <<- list(text = message, type = type, duration = duration)
            })
    }))
    record <- file.path(root, "analysis_record.json")
    writeWorkflowAnalysisRecord(manager, record)
    state_ref <- manager$.context$snapshot$states$enrichment_completed
    list(run = run, error = captured$error, conditions = captured$conditions, notes = notes, record = record,
        paths = paths, kernel_calls = kernel_calls, workflow = workflow,
        uniprot_rows = k175UniprotRows(manager),
        outcomes = if (!is.null(state_ref)) readComponentMetadata(manager, state_ref)$config$outcomes)
}

annotations_or_empty <- function(annotations) if (is.data.frame(annotations)) annotations else createEmptyUniprotTable()

#' A table through the peptide checkpoint's save and open (the session store, `saveDIAStateValue(codec = "rds")`)
k175SessionRoundTrip <- function(table) {
    recipe <- readWorkflowCheckpointRecipe("small")
    recipe$n_proteins <- 8L; recipe$n_samples <- 8L
    recipe$peptides_per_protein <- 4L; recipe$forms_per_peptide <- 2L
    object <- makeWorkflowCheckpointPeptide(recipe)
    root <- withr::local_tempdir()
    make_paths <- function(name) {
        paths <- list(source_dir = file.path(root, paste0(name, "-source")),
            peptide_qc_dir = file.path(root, paste0(name, "-peptide")))
        lapply(paths, dir.create)
        paths
    }
    make_workflow <- function(native) shiny::reactiveValues(state_manager = native, data_tbl = object@peptide_data,
        data_cln = object@peptide_data, data_type = "peptide", data_format = "diann",
        config_list = list(globalParameters = list(workflow_type = "DIA")), design_matrix = object@design_matrix,
        column_mapping = list(run_col = "Run", protein_col = "Protein.Group", peptide_col = "Stripped.Sequence",
            raw_quantity_col = "Precursor.Quantity", norm_quantity_col = "Precursor.Normalised"))
    save_paths <- make_paths("save")
    writeLines("accepted output", file.path(save_paths$peptide_qc_dir, "peptide.tsv"))
    holder <- new.env(parent = emptyenv())
    shiny::testServer(function(input, output, session) {}, {
        session$userData$peptide_checkpoint_design <- function() list(design_matrix = object@design_matrix,
            groups = unique(object@design_matrix$group), factors = character(), contrasts = data.frame(),
            removed_samples = character(), formula = "~ 0 + group")
        native <- WorkflowState$new(audit_enabled = FALSE)
        native$setWorkflowType("DIA")
        native$saveState("raw_data_s4", object, list(), "raw_data_s4")
        saving <- make_workflow(native)
        saving$uniprot_dat_cln <- table
        publication <- saveDIAPeptideCheckpoint(saving, session, save_paths)
        opening <- make_workflow(WorkflowState$new(audit_enabled = FALSE))
        opening$state_manager$setWorkflowType("DIA")
        openDIAPeptideCheckpoint(publication$archive_path, opening, session, make_paths("open"))
        holder$value <- opening$uniprot_dat_cln
        closeDIAStateAdapter(opening$state_manager)
    })
    holder$value
}

#' The annotation table with the first 25 of its entries removed and the outcome of a request that failed with 503
k175Incomplete <- function(full) {
    kept <- full[-seq_len(25L), , drop = FALSE]
    attr(kept, "annotation_outcome") <- list(class = "prot_annotation_uniprot_incomplete", reason = "incomplete",
        requested = nrow(full), taxon_id = 83333, requests_total = 3L, requests_failed = 1L, requests_not_sent = 0L,
        requests_truncated = 0L, accessions_requested = nrow(full), accessions_not_looked_up = 25L,
        hits_beyond_rows = 0L, statuses = 503L)
    kept
}

k175IncompleteNotice <- "Enrichment ran with incomplete protein annotation: 25 of 70 proteins have no annotation and were left out."
k175Refusal <- paste("Enrichment needs protein annotation and this project has none. Save the design again to",
    "retrieve it, then run the later steps again.")

test_that("KIT-175 T3: a record with a failed row validates and reads back with release NA", {
    k175Offline()
    root <- withr::local_tempdir()
    full <- k175ProjectAnnotations(withr::local_tempdir())
    k175Mock(function(requested, i, req) k175Response(503L))
    project <- k175Project(root, annotations = full, before = function(manager) {
        getUniprotAnnotations(data.frame(Protein.Ids = k175Accessions(3L), stringsAsFactors = FALSE),
            cache_dir = file.path(root, "uniprot-cache"), taxon_id = 9606, api_delay = 0)
    })
    expect_null(project$error)
    record <- jsonlite::read_json(project$record, simplifyVector = FALSE)
    expect_true(validateWorkflowAnalysisRecord(record)$valid)
    external <- workflowAnalysisRecordTables(record)$external_resources
    failed <- external[external$service == "uniprot" & external$release_kind == "failed", , drop = FALSE]
    expect_identical(nrow(failed), 1L)
    expect_true(is.na(failed$release))
})

test_that("KIT-175 T5 (a), (e): DA statistics do not read the annotation; integration features fall back to accessions", {
    k175Offline()
    full <- k175ProjectAnnotations(withr::local_tempdir())
    with_annotation <- k175Project(withr::local_tempdir(), annotations = full)
    without <- k175Project(withr::local_tempdir(), da_annotations = createEmptyUniprotTable(), cache = full)
    expect_null(with_annotation$error)
    expect_null(without$error)
    read_da <- function(project) {
        file <- list.files(project$paths$da_output_dir, "_long_annot\\.tsv$", full.names = TRUE)[[1L]]
        table <- utils::read.delim(file, stringsAsFactors = FALSE)
        table <- table[order(table[[1L]]), , drop = FALSE]
        rownames(table) <- NULL
        table
    }
    a <- read_da(with_annotation)
    b <- read_da(without)
    # The statistics: the numeric columns both tables hold (annotation columns such as a score exist on one side only).
    numeric <- intersect(names(a)[vapply(a, is.numeric, NA)], names(b)[vapply(b, is.numeric, NA)])
    expect_true(length(numeric) >= 2L)
    expect_identical(a[[1L]], b[[1L]])
    expect_identical(a[numeric], b[numeric])
    features <- protIntegrationInputFeatures(c("P10001", "P10002;P10003", "P10004-2"), createEmptyUniprotTable())
    expect_identical(features$gene_names, rep(NA_character_, 3L))
    expect_identical(features$display_name, c("P10001", "P10002;P10003", "P10004-2"))
})

test_that("KIT-175 T5 (b): enrichment refuses on an empty or missing annotation before any service or GO step", {
    k175Offline()
    # Empty table on the bus.
    empty <- k175Project(withr::local_tempdir(), bus = createEmptyUniprotTable())
    expect_true(inherits(empty$error, "workflow_enrichment_annotation"))
    expect_identical(conditionMessage(empty$error), k175Refusal)
    expect_identical(empty$kernel_calls, 0L)
    expect_null(empty$outcomes)
    # NULL (the step's own retrieval errors): refused before KIT-173's GO outcome step (I04).
    testthat::local_mocked_bindings(getUniprotAnnotations = function(...) stop("the cache directory cannot be read"),
        .package = "MultiScholaR")
    missing <- k175Project(withr::local_tempdir())
    expect_true(inherits(missing$error, "workflow_enrichment_annotation"))
    expect_identical(missing$kernel_calls, 0L)
    expect_identical(missing$notes[[1L]]$text, "Protein annotation was retrieved but could not be processed (the cache directory cannot be read).")
})

test_that("KIT-175 T5 (c): an incomplete table is reported on the bus, from the project file and after a session round trip", {
    k175Offline()
    incomplete <- k175Incomplete(k175ProjectAnnotations(withr::local_tempdir()))
    check <- function(project) {
        expect_null(project$error)
        expect_identical(project$kernel_calls, 1L)
        texts <- vapply(project$notes, function(note) note$text, "")
        expect_true(k175IncompleteNotice %in% texts)
        expect_identical(project$outcomes$annotation_error,
            "1 of 3 requests failed (status 503); 25 of 70 accessions not looked up")
    }
    check(k175Project(withr::local_tempdir(), bus = incomplete, da_annotations = incomplete))
    check(k175Project(withr::local_tempdir(), annotations = incomplete))
    restored <- k175SessionRoundTrip(incomplete)
    expect_identical(attr(restored, "annotation_outcome"), attr(incomplete, "annotation_outcome"))
    check(k175Project(withr::local_tempdir(), bus = restored, da_annotations = restored))
})

test_that("KIT-175 T5 (d): with the organism filter on, an incomplete annotation still runs and records the filter", {
    k175Offline()
    incomplete <- k175Incomplete(k175ProjectAnnotations(withr::local_tempdir()))
    project <- k175Project(withr::local_tempdir(), annotations = incomplete, organism_filter = TRUE)
    expect_null(project$error)
    expect_identical(project$outcomes$organism_filter_applied, FALSE)
    texts <- vapply(project$notes, function(note) note$text, "")
    expect_true(paste("Multi-species filtering enabled but no organism mapping available.",
        "Proceeding without filtering.") %in% texts)
})

test_that("KIT-175 T5 (f): the enrichment step's own retrieval on its real path (zero identifiers)", {
    k175Offline()
    # (f1) No cache file: the zero-identifier notice, then the refusal; nothing cached; one not_queried row.
    log <- k175Mock(k175Answer())
    project <- k175Project(withr::local_tempdir())
    expect_true(inherits(project$error, "workflow_enrichment_annotation"))
    expect_identical(project$notes[[1L]]$text, paste("Protein annotation could not be requested at this step: it has",
        "no protein identifiers to look up. Save the design again to retrieve the annotation."))
    expect_identical(log$calls, 0L)
    expect_null(project$outcomes)
    cache_file <- file.path(project$paths$results_dir, "cache", "uniprot_annotations", "uniprot_annotations.RDS")
    expect_false(file.exists(cache_file))
    expect_length(project$uniprot_rows, 1L)
    expect_identical(project$uniprot_rows[[1L]]$release_kind, "not_queried")
    expect_identical(project$uniprot_rows[[1L]]$query_summary,
        "not queried: this step had no protein identifiers (annotation is retrieved at the design step)")
    # (f2) A zero-row cache file from an earlier version: removed, then as (f1).
    project <- k175Project(withr::local_tempdir(), cache = createEmptyUniprotTable())
    expect_true(inherits(project$error, "workflow_enrichment_annotation"))
    expect_false(file.exists(file.path(project$paths$results_dir, "cache", "uniprot_annotations",
        "uniprot_annotations.RDS")))
    # (f3) A cache file with rows: the run proceeds on it, no record row (a pin).
    full <- k175ProjectAnnotations(withr::local_tempdir())
    project <- k175Project(withr::local_tempdir(), da_annotations = full, cache = full)
    expect_null(project$error)
    expect_length(project$uniprot_rows, 0L)
})

test_that("KIT-175 T5 (f): the resolver keeps the reason and the redacted detail of its own retrieval's error", {
    k175Offline()
    k175Mock(function(requested, i, req) stop(k175TransportError()))
    accessions <- k175Accessions(3L)
    object <- methods::new("ProteinQuantitativeData",
        protein_quant_table = data.frame(Protein.Ids = accessions, S1 = 1, S2 = 2, stringsAsFactors = FALSE),
        design_matrix = data.frame(Run = c("S1", "S2"), group = c("A", "B"), stringsAsFactors = FALSE),
        sample_id = "Run", group_id = "group", protein_id_column = "Protein.Ids",
        protein_id_table = data.frame(Protein.Ids = accessions, stringsAsFactors = FALSE))
    root <- withr::local_tempdir()
    resolved <- NULL
    utils::capture.output(suppressMessages(suppressWarnings(resolved <- resolveProtEnrichUniprotAnnotations(
        new.env(), list(source_dir = file.path(root, "source"), results_dir = file.path(root, "results")),
        currentS4Object = object, organismTaxid = "9606", getUniprotAnnotationsFn = function(...)
            getUniprotAnnotations(..., api_delay = 0)))))
    expect_identical(resolved$creationReason, "request_failed")
    expect_identical(resolved$creationDetail, "Could not resolve host: rest.uniprot.org")
})

test_that("KIT-175 T5 (g): a load error followed by a retrieval from the cache completes with the notice and outcome", {
    k175Offline()
    full <- k175ProjectAnnotations(withr::local_tempdir())
    project <- k175Project(withr::local_tempdir(), annotations = "this is not an RDS file", da_annotations = full,
        cache = full)
    expect_null(project$error)
    texts <- vapply(project$notes, function(note) note$text, "")
    expect_true(any(startsWith(texts, "The saved protein annotation could not be read (")))
    expect_match(project$outcomes$annotation_error, "^annotation file could not be read \\(")
})

test_that("KIT-175 T5 (h): a session saved with an empty annotation reopens with it, and enrichment refuses", {
    k175Offline()
    restored <- k175SessionRoundTrip(createEmptyUniprotTable())
    expect_identical(nrow(restored), 0L)
    project <- k175Project(withr::local_tempdir(), bus = restored)
    expect_true(inherits(project$error, "workflow_enrichment_annotation"))
    expect_identical(project$kernel_calls, 0L)
})

# --- T6 success unchanged (regression test; both doubles) -------------------------------------------------------------

k175T6Entries <- c("P10001", "P10002", "P10003")

#' Both doubles serve the same three-row body: httr's (base) and httr2's (head)
k175T6Doubles <- function(env = parent.frame()) {
    body <- k175Body(k175T6Entries)
    log <- new.env(parent = emptyenv())
    log$calls <- 0L
    testthat::local_mocked_bindings(
        GET = function(url, query, ...) {
            log$calls <- log$calls + 1L
            list(status = 200L, text = body)
        },
        status_code = function(x) x$status,
        content = function(x, as = "parsed", encoding = NULL) x$text,
        headers = function(x) list(`x-uniprot-release` = "2026_04", `x-total-results` = "3"),
        timeout = function(seconds) seconds,
        .package = "httr", .env = env)
    httr2::local_mocked_responses(function(req) {
        log$calls <- log$calls + 1L
        k175Response(200L, body, list(`X-Total-Results` = "3"))
    }, env = env)
    log
}

k175ExpectT6Table <- function(tbl) {
    expect_identical(names(tbl), c("Entry", "Entry.Name", "Protein_names", "gene_names", "Organism", "Length",
        "Reviewed", "Protein_existence", "Annotation", "From", "go_term_go_cellular_compartment",
        "go_id_go_cellular_compartment", "annotation_score"))
    expect_identical(tbl$Entry, k175T6Entries)
    expect_identical(tbl$Entry.Name, c("NAME1_HUMAN", "NAME2_HUMAN", "NAME3_HUMAN"))
    expect_identical(tbl$Protein_names, c("Protein P10001", "Protein P10002", "Protein P10003"))
    expect_identical(tbl$gene_names, c("GENEP10001", "GENEP10002", "GENEP10003"))
    expect_identical(tbl$Organism, rep("Homo sapiens (Human)", 3L))
    expect_identical(tbl$Length, rep(100L, 3L))
    expect_identical(tbl$Reviewed, rep("reviewed", 3L))
    expect_identical(tbl$Protein_existence, rep("Evidence at protein level", 3L))
    expect_identical(tbl$Annotation, rep(5, 3L))
    expect_identical(tbl$From, k175T6Entries)
    expect_identical(tbl$go_term_go_cellular_compartment, rep("cellular_component", 3L))
    expect_identical(tbl$go_id_go_cellular_compartment, rep("GO:0005575", 3L))
    expect_identical(tbl$annotation_score, rep(0, 3L))
    expect_null(attr(tbl, "annotation_outcome"))
}

test_that("KIT-175 T6 (regression): a successful retrieval, a cache hit and the step's own retrieval are unchanged", {
    k175Offline()
    manager <- k175Manager("k175-t6")
    log <- k175T6Doubles()

    # A fresh retrieval of three rows.
    cache_dir <- withr::local_tempdir()
    out <- k175Annotate(manager, k175T6Entries, cache_dir)
    expect_null(out$error)
    expect_length(out$conditions, 0L)
    k175ExpectT6Table(out$value)
    k175ExpectT6Table(readRDS(k175CacheFile(cache_dir)))
    rows <- k175UniprotRows(manager)
    expect_length(rows, 1L)
    expect_identical(rows[[1L]][c("service", "release", "release_kind", "query_summary")],
        list(service = "uniprot", release = "2026_04", release_kind = "header",
            query_summary = "uniprotkb search, organism_id 9606"))

    # A cache hit: no request, no record entry, the same table.
    calls <- log$calls
    out <- k175Annotate(manager, k175T6Entries, cache_dir)
    expect_identical(log$calls, calls)
    expect_length(k175UniprotRows(manager), 1L)
    k175ExpectT6Table(out$value)

    # The design step: the same table on the bus and in uniprot_dat_cln.RDS, and today's success notice.
    source <- saveDIATableComponent(manager$.context, k175DesignTable(c("P10001;P10002", "P10003")))
    paths <- k175DesignPaths()
    step <- k175DesignStep(manager, source, paths)
    result <- step$run()
    k175ExpectT6Table(result$workflow_data$uniprot_dat_cln)
    k175ExpectT6Table(readRDS(file.path(paths$source_dir, "uniprot_dat_cln.RDS")))
    notices <- k175AnnotationNotices(step)
    expect_length(notices, 1L)
    expect_identical(notices[[1L]][c("text", "type")],
        list(text = "UniProt annotations retrieved successfully.", type = "message"))

    # The enrichment step's own retrieval, called with an object that carries the three identifiers.
    object <- methods::new("ProteinQuantitativeData",
        protein_quant_table = data.frame(Protein.Ids = k175T6Entries, S1 = 1, S2 = 2, stringsAsFactors = FALSE),
        design_matrix = data.frame(Run = c("S1", "S2"), group = c("A", "B"), stringsAsFactors = FALSE),
        sample_id = "Run", group_id = "group", protein_id_column = "Protein.Ids",
        protein_id_table = data.frame(Protein.Ids = k175T6Entries, stringsAsFactors = FALSE))
    root <- withr::local_tempdir()
    resolved <- NULL
    utils::capture.output(suppressMessages(resolved <- resolveProtEnrichUniprotAnnotations(new.env(),
        list(source_dir = file.path(root, "source"), results_dir = file.path(root, "results")),
        currentS4Object = object, organismTaxid = "9606",
        getUniprotAnnotationsFn = function(...) getUniprotAnnotations(..., api_delay = 0))))
    expect_null(resolved$creationError)
    k175ExpectT6Table(resolved$uniprotDatCln)
})

# --- T7 replay -------------------------------------------------------------------------------------------------------

#' goPropReplay (test-go-propagation.R's helper) with the project's annotation file removed from the copy when asked
k175Replay <- function(project_root, record, scratch, drop_annotation = FALSE) {
    copy <- file.path(scratch, "project")
    dir.create(copy, recursive = TRUE)
    stopifnot(system2("cp", c("-a", file.path(project_root, "."), copy)) == 0L)
    if (drop_annotation) unlink(file.path(copy, "source", "uniprot_dat_cln.RDS"))
    h <- new.env(parent = globalenv())
    sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), h)
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    prepared <- h$kitPrepareReplayProject(copy, c("initial", "correlation_filtered"), api)
    out <- file.path(scratch, "replay")
    dir.create(out)
    seen <- new.env(parent = emptyenv())
    captured <- k175Capture(replayWorkflowAnalysisRecord(record, prepared, character(), goPropPaths(copy), out,
        on_step = function(step, manager) {
            if (identical(step$operation, "enrichment")) {
                reference <- manager$.context$snapshot$states$enrichment_completed
                seen$outcomes <- readComponentMetadata(manager, reference)$config$outcomes
            }
        }))
    c(captured, list(outcomes = seen$outcomes))
}

test_that("KIT-175 T7: replaying an enrichment with no annotation is refused; one with annotation replays unchanged", {
    k175Offline()
    log <- k175Mock(k175Answer())
    root <- withr::local_tempdir()
    full <- k175ProjectAnnotations(withr::local_tempdir())
    origin <- NULL
    utils::capture.output(origin <- suppressMessages(suppressWarnings(goPropRunProject(root, full))))
    record <- jsonlite::read_json(origin$record, simplifyVector = FALSE)
    without <- k175Replay(root, record, withr::local_tempdir(), drop_annotation = TRUE)
    expect_true(inherits(without$error, "workflow_enrichment_annotation"))
    expect_identical(log$calls, 0L)
    with <- k175Replay(root, record, withr::local_tempdir())
    expect_null(with$error)
    step <- Filter(function(s) identical(s$operation, "enrichment"), record$steps)[[1L]]
    expect_identical(as.character(unlist(with$outcomes$annotation_sha256)),
        as.character(unlist(step$outcomes$annotation_sha256)))
    expect_false(is.na(as.character(unlist(step$outcomes$annotation_sha256))))
})

# --- T8 the request's retry policy -----------------------------------------------------------------------------------

test_that("KIT-175 T8: each request has a 30 s time-out, three attempts, five transient statuses and a capped wait", {
    k175Offline()
    manager <- k175Manager("k175-t8")
    log <- k175Mock(k175Answer())
    k175Annotate(manager, k175Accessions(3L), withr::local_tempdir())
    expect_identical(log$calls, 1L)
    req <- log$requests[[1L]]
    expect_equal(req$policies$retry_max_tries, 3)
    expect_true(isTRUE(req$policies$retry_on_failure))
    for (status in c(429L, 500L, 502L, 503L, 504L)) {
        expect_true(isTRUE(req$policies$retry_is_transient(k175Response(status))), info = status)
    }
    for (status in c(400L, 404L)) {
        expect_false(isTRUE(req$policies$retry_is_transient(k175Response(status))), info = status)
    }
    expect_lte(req$policies$retry_after(k175Response(503L, headers = list(`Retry-After` = "3600"))), 10)
    expect_equal(req$options$timeout_ms, 30000)
    expect_identical(req$options$useragent, "MultiScholaR")
})
