# KIT-141: the integration enrichment kernels (R/utils_integration_enrichment.R). STRING is mocked with httr2's
# mocked responses, and httr2's transfer function is replaced by one that fails as a refused connection: no network
# and no real key. The key used here is an invented sentinel string.

SENTINEL_KEY <- "kit141-sentinel-not-a-real-key-Zq7w"
STRING_HOST <- "https://version-12-5.string-db.org"

stringRanked <- function(n = 12L) {
    data.frame(identifier = sprintf("P%05d", seq_len(n)), value = seq(2, -2, length.out = n),
        stringsAsFactors = FALSE)
}

stringEndpoint <- function(request) sub("[?].*$", "", sub("^.*/", "", httr2::req_get_url(request)))

stringFields <- function(request) {
    body <- httr2::req_get_body(request)
    if (is.null(body)) return(list())
    lapply(body, \(value) utils::URLdecode(as.character(value)))
}

# A STRING double: answers each endpoint, records every request, and (as a real call does) leaves the request in
# httr2's own record of the last request, which the client must clear.
newStringMock <- function(states = c("running", "success"), download_url = paste0(STRING_HOST,
    "/api/tsv/download?networkId=abc"), status_message = "done", results_tsv = paste0(
    "category\tterm\tdescription\tfalse_discovery_rate\n", "KEGG\thsa00010\tGlycolysis\t0.001\n"), overrides = list()) {
    log <- new.env()
    log$requests <- list()
    log$polls <- 0L
    json <- function(x, status = 200L) {
        httr2::response(status_code = status, headers = list(`Content-Type` = "application/json"),
            body = charToRaw(as.character(jsonlite::toJSON(x, auto_unbox = TRUE))))
    }
    # httr2 requires the mock's argument to be named `req`.
    handler <- function(req) {
        log$requests[[length(log$requests) + 1L]] <- req
        assign("last_request", req, envir = get("the", envir = asNamespace("httr2")))
        endpoint <- stringEndpoint(req)
        if (is.function(overrides[[endpoint]])) return(overrides[[endpoint]](req))
        switch(endpoint,
            version = json(list(list(string_version = "12.5", string_stable_address = STRING_HOST))),
            valuesranks_enrichment_submit = json(list(list(job_id = "job42", status = "submitted"))),
            valuesranks_enrichment_status = {
                log$polls <- log$polls + 1L
                json(list(list(job_id = "job42", status = states[[min(log$polls, length(states))]],
                    message = status_message, page_url = paste0(STRING_HOST, "/cgi/network?networkId=abc"),
                    download_url = download_url, graph_url = paste0(STRING_HOST, "/api/image/x"))))
            },
            download = httr2::response(status_code = 200L,
                headers = list(`Content-Type` = "text/tab-separated-values"), body = charToRaw(results_tsv)),
            stop("unexpected endpoint ", endpoint)
        )
    }
    list(handler = handler, log = log)
}

# A network failure one level below httr2's mocked responses (run-tests 44783719082: a double that only declined the
# request reached no failure inside the client). Two layers, both in place in every STRING test, so no request can leave
# the machine:
#   - httr2's send function `req_perform` is replaced: for the endpoints named failing it does what httr2 does on a
#     refused connection (it keeps the request as its last request, builds its own failure, class `httr2_failure`,
#     with the curl error as parent, attaches the request, and with it the key in the body, and signals it); for any
#     other endpoint it calls httr2's own `req_perform`, which answers from the STRING double;
#   - httr2's transfer function `curl_fetch` (the call into curl) fails as a refused connection does, should a request
#     ever get past the double.
# `network$sent` counts the failures the send function raised, `network$transfers` the requests that reached the
# transfer function (none expected).
network <- new.env()
network$failing <- character()
network$sent <- 0L
network$transfers <- 0L
curlConnectError <- function(url) {
    structure(class = c("curl_error_couldnt_connect", "curl_error", "error", "condition"), list(
        message = sprintf("Failed to connect to %s port 443 after 0 ms: Could not connect to server",
            httr2::url_parse(url)$hostname), call = NULL))
}
curlConnectFailure <- function(handle, url, path) {
    network$transfers <- network$transfers + 1L
    stop(curlConnectError(url))
}
httr2Perform <- httr2::req_perform
failingPerform <- function(req, path = NULL, verbosity = NULL, mock = getOption("httr2_mock", NULL),
    error_call = rlang::current_env()) {
    if (!stringEndpoint(req) %in% network$failing) {
        return(httr2Perform(req, path = path, verbosity = verbosity, mock = mock, error_call = error_call))
    }
    network$sent <- network$sent + 1L
    the <- get("the", envir = asNamespace("httr2"))
    assign("last_request", req, envir = the)
    assign("last_response", NULL, envir = the)
    condition <- get("curl_cnd", envir = asNamespace("httr2"))(curlConnectError(httr2::req_get_url(req)))
    condition$request <- req
    rlang::cnd_signal(condition)
}
localNetworkFailure <- function(failing = character(), env = parent.frame()) {
    network$failing <- failing
    network$sent <- 0L
    network$transfers <- 0L
    # A third guard: a proxy that refuses every connection, so curl cannot reach the network either.
    withr::local_envvar(c(https_proxy = "http://127.0.0.1:9", HTTPS_PROXY = "http://127.0.0.1:9",
        http_proxy = "http://127.0.0.1:9", ALL_PROXY = "http://127.0.0.1:9", no_proxy = "", NO_PROXY = ""), .local_envir = env)
    testthat::local_mocked_bindings(curl_fetch = curlConnectFailure, req_perform = failingPerform, .package = "httr2",
        .env = env)
}
networkFailure <- function(req) NULL

runString <- function(mock, ranked = stringRanked(), key = SENTINEL_KEY, options = integrationStringOptions(
    poll_interval_seconds = 0), failing = character()) {
    localNetworkFailure(failing)
    httr2::local_mocked_responses(mock$handler)
    integrationStringRankedEnrichment(ranked, key, 9606, options)
}

containsSentinel <- function(x) length(grepRaw(SENTINEL_KEY, serialize(x, NULL), fixed = TRUE)) > 0L

test_that("the network-failure double raises httr2's own failure, carrying the request and its key", {
    skip_if_not_installed("httr2")
    withr::defer(integrationForgetHttpRequest())
    request <- httr2::request(paste0(STRING_HOST, "/api/json/valuesranks_enrichment_submit")) |>
        httr2::req_body_form(api_key = SENTINEL_KEY)
    # The send function's failure, and httr2's own failure from a refused transfer, carry the request and its key.
    local({
        localNetworkFailure("valuesranks_enrichment_submit")
        condition <- tryCatch(httr2::req_perform(request), error = \(e) e)
        expect_identical(network$sent, 1L)
        expect_s3_class(condition, "httr2_failure")
        expect_s3_class(condition$parent, "curl_error_couldnt_connect")
        expect_true(containsSentinel(condition$request))
        expect_true(containsSentinel(httr2::last_request()))
    })
    local({
        localNetworkFailure()
        condition <- tryCatch(httr2::req_perform(request), error = \(e) e)
        expect_identical(network$transfers, 1L)
        expect_s3_class(condition, "httr2_failure")
        expect_s3_class(condition$parent, "curl_error_couldnt_connect")
        expect_true(containsSentinel(condition))
    })
})

test_that("the STRING client sends the documented requests and returns the release, a request summary and the response", {
    skip_if_not_installed("httr2")
    mock <- newStringMock()
    ranked <- stringRanked()
    result <- runString(mock, ranked)

    requests <- mock$log$requests
    expect_identical(vapply(requests, stringEndpoint, ""), c("version", "valuesranks_enrichment_submit",
        "valuesranks_enrichment_status", "valuesranks_enrichment_status", "download"))
    expect_true(all(startsWith(vapply(requests, httr2::req_get_url, ""), paste0(STRING_HOST, "/api/"))))
    expect_identical(vapply(requests, httr2::req_get_method, ""), c("GET", "POST", "POST", "POST", "GET"))

    submit <- stringFields(requests[[2L]])
    expect_setequal(names(submit), c("species", "caller_identity", "identifiers", "api_key", "ge_fdr",
        "ge_enrichment_rank_direction"))
    expect_identical(submit$species, "9606")
    expect_identical(submit$caller_identity, "MultiScholaR")
    expect_identical(submit$identifiers, paste(ranked$identifier, ranked$value, sep = "\t", collapse = "\n"))
    expect_identical(submit$ge_fdr, "0.05")
    expect_identical(submit$ge_enrichment_rank_direction, "-1")
    expect_setequal(names(stringFields(requests[[3L]])), c("api_key", "job_id"))
    expect_identical(stringFields(requests[[3L]])$job_id, "job42")

    expect_identical(result$service, "STRING")
    expect_identical(result$version, list(string_version = "12.5", string_stable_address = STRING_HOST))
    expect_identical(result$request$species, "9606")
    expect_identical(result$request$identifiers, 12L)
    expect_identical(result$request$polls, 2L)
    expect_identical(result$request$identifiers_sha256,
        digest::digest(submit$identifiers, algo = "sha256", serialize = FALSE))
    expect_false("api_key" %in% names(result$request))
    expect_identical(result$response$submission$job_id, "job42")
    expect_identical(result$response$status$status, "success")
    expect_s3_class(result$terms, "data.frame")
    expect_identical(result$terms$term, "hsa00010")
    expect_identical(result$terms$false_discovery_rate, 0.001)
})

test_that("the STRING key travels only in POST bodies, never in an address or a header", {
    skip_if_not_installed("httr2")
    # The service's download address carries a key: the client takes it out of the address and posts it instead.
    mock <- newStringMock(download_url = paste0(STRING_HOST, "/api/tsv/download?networkId=abc&api_key=",
        SENTINEL_KEY))
    runString(mock)
    for (request in mock$log$requests) {
        expect_false(grepl(SENTINEL_KEY, httr2::req_get_url(request), fixed = TRUE))
        expect_false(containsSentinel(unclass(request$headers)))
        expect_false(containsSentinel(request$options))
        fields <- stringFields(request)
        if (SENTINEL_KEY %in% unlist(fields)) {
            expect_identical(httr2::req_get_method(request), "POST")
            expect_identical(fields$api_key, SENTINEL_KEY)
        }
    }
    with_key <- vapply(mock$log$requests, \(request) identical(stringFields(request)$api_key, SENTINEL_KEY), TRUE)
    expect_identical(vapply(mock$log$requests, stringEndpoint, "")[with_key],
        c("valuesranks_enrichment_submit", "valuesranks_enrichment_status", "valuesranks_enrichment_status",
            "download"))
    download <- mock$log$requests[[5L]]
    expect_identical(httr2::req_get_url(download), paste0(STRING_HOST, "/api/tsv/download?networkId=abc"))
})

test_that("a planted sentinel key appears in no returned value, message, warning, condition, httr2 record or file", {
    skip_if_not_installed("httr2")
    work <- withr::local_tempdir()
    withr::local_dir(work)
    scanned_text <- character()
    scanned_conditions <- list()
    capture <- function(code) {
        withCallingHandlers(code,
            message = function(m) {
                scanned_text <<- c(scanned_text, conditionMessage(m))
                scanned_conditions <<- c(scanned_conditions, list(m))
                invokeRestart("muffleMessage")
            },
            warning = function(w) {
                scanned_text <<- c(scanned_text, conditionMessage(w))
                scanned_conditions <<- c(scanned_conditions, list(w))
                invokeRestart("muffleWarning")
            })
    }
    # The service echoes the key everywhere it could: in the download address, the status message, the table.
    echo <- newStringMock(download_url = paste0(STRING_HOST, "/api/tsv/download?networkId=abc&api_key=",
        SENTINEL_KEY), status_message = paste("job of key", SENTINEL_KEY, "finished"),
        results_tsv = paste0("term\tdescription\nhsa00010\tsubmitted with ", SENTINEL_KEY, "\n"))
    result <- capture(runString(echo))
    expect_false(containsSentinel(result))
    expect_true(grepl("[removed]", result$response$status$message, fixed = TRUE))

    # A refusal whose service message carries the key.
    refused <- newStringMock(overrides = list(valuesranks_enrichment_submit = function(request) {
        httr2::response(status_code = 200L, headers = list(`Content-Type` = "application/json"),
            body = charToRaw(sprintf('[{"status":"error","message":"key %s is not valid"}]', SENTINEL_KEY)))
    }))
    condition <- capture(tryCatch(runString(refused), integration_enrichment_service_failed = \(e) e))
    expect_s3_class(condition, "integration_enrichment_service_failed")
    expect_false(containsSentinel(condition))
    expect_false(grepl(SENTINEL_KEY, conditionMessage(condition), fixed = TRUE))
    expect_null(condition$parent)
    # A network failure at the submission, which carries the key: httr2's own failure holds the request.
    unreachable <- newStringMock(overrides = list(valuesranks_enrichment_submit = networkFailure))
    condition <- capture(tryCatch(runString(unreachable, failing = "valuesranks_enrichment_submit"),
        integration_enrichment_service_unreachable = \(e) e))
    expect_s3_class(condition, "integration_enrichment_service_unreachable")
    # The failure came from the send function, after the key was in the request, and nothing reached the transfer.
    expect_gte(network$sent, 1L)
    expect_identical(network$transfers, 0L)
    # The double answered the version; the submission never reached it (the send function failed first).
    expect_identical(vapply(unreachable$log$requests, stringEndpoint, ""), "version")
    expect_false(containsSentinel(condition))
    expect_null(condition$parent)
    expect_false(grepl(SENTINEL_KEY, conditionMessage(condition), fixed = TRUE))

    expect_false(containsSentinel(scanned_text))
    expect_false(containsSentinel(scanned_conditions))
    expect_null(httr2::last_request())
    expect_null(httr2::last_response())
    files <- c(list.files(work, recursive = TRUE, full.names = TRUE, all.files = TRUE),
        list.files(tempdir(), recursive = TRUE, full.names = TRUE, all.files = TRUE))
    sizes <- file.info(files)$size
    files <- files[!is.na(sizes) & sizes < 64 * 1024^2]
    on_disk <- vapply(files, \(path) length(grepRaw(SENTINEL_KEY, readBin(path, "raw", file.info(path)$size),
        fixed = TRUE)) > 0L, TRUE)
    expect_identical(sum(on_disk), 0L)
})

test_that("STRING failures are classed refusals with no partial result", {
    skip_if_not_installed("httr2")
    status <- function(code) function(request) httr2::response(status_code = code)
    expect_refusal <- function(mock, class, options = integrationStringOptions(poll_interval_seconds = 0),
        key = SENTINEL_KEY, failing = character()) {
        condition <- expect_error(runString(mock, key = key, options = options, failing = failing), class = class)
        expect_s3_class(condition, "integration_refused")
        expect_s3_class(condition, "integration_enrichment_refused")
        mock
    }
    mock <- expect_refusal(newStringMock(overrides = list(version = status(500L))),
        "integration_enrichment_service_failed")
    expect_length(mock$log$requests, 1L)
    mock <- expect_refusal(newStringMock(overrides = list(version = networkFailure)),
        "integration_enrichment_service_unreachable", failing = "version")
    expect_length(mock$log$requests, 0L)
    expect_identical(network$sent, 1L)
    expect_identical(network$transfers, 0L)
    expect_refusal(newStringMock(overrides = list(version = function(request) {
        httr2::response(status_code = 200L, body = charToRaw("not json"))
    })), "integration_enrichment_service_failed")
    mock <- expect_refusal(newStringMock(states = "failed"), "integration_enrichment_service_failed")
    expect_false("download" %in% vapply(mock$log$requests, stringEndpoint, ""))
    mock <- expect_refusal(newStringMock(states = "running"), "integration_enrichment_service_timeout",
        options = integrationStringOptions(poll_interval_seconds = 0, max_polls = 2L))
    expect_identical(sum(vapply(mock$log$requests, stringEndpoint, "") == "valuesranks_enrichment_status"), 2L)
    expect_refusal(newStringMock(overrides = list(download = status(404L))), "integration_enrichment_service_failed")
    expect_refusal(newStringMock(download_url = "http://version-12-5.string-db.org/api/tsv/download?x=1"),
        "integration_enrichment_service_failed")
    mock <- expect_refusal(newStringMock(), "integration_enrichment_key_missing", key = "")
    expect_length(mock$log$requests, 0L)
    mock <- newStringMock()
    localNetworkFailure()
    httr2::local_mocked_responses(mock$handler)
    expect_error(integrationStringRankedEnrichment(stringRanked(), SENTINEL_KEY, "human",
        integrationStringOptions(poll_interval_seconds = 0)), class = "integration_enrichment_parameters_invalid")
    expect_error(integrationStringRankedEnrichment(stringRanked(), SENTINEL_KEY, 9606,
        integrationStringOptions(host = "http://string-db.org")), class = "integration_enrichment_parameters_invalid")
    expect_length(mock$log$requests, 0L)
})

test_that("identifiers are mapped by the declared rule of each omic and counted", {
    accessions <- c("Q9Y261", "P04637", "P38398", "O15111", "Q13315", "P00533", "P01116", "P31749", "Q00987",
        "P42336")
    ranked <- data.frame(feature_id = c("P11111;Q22222", "P11111-2", "sp|P33333|ABC_HUMAN", "not_an_accession",
        "P44444", accessions), value = c(3, -5, 1, 2, NA, seq(1, 10)), stringsAsFactors = FALSE)
    mapped <- integrationEnrichmentMapIdentifiers(ranked, "proteomics")
    expect_identical(mapped$counts, list(features = 15L, mapped = 12L, unmapped = 1L, duplicated = 1L,
        missing_value = 1L))
    expect_identical(mapped$table$identifier, c("P11111", "P33333", accessions))
    expect_identical(mapped$table$feature_id[[1L]], "P11111-2")
    expect_identical(mapped$mapping$status, c("duplicated", "mapped", "mapped", "unmapped", "missing_value",
        rep("mapped", 10L)))
    expect_identical(mapped$rule$id_type, "uniprot")
    expect_identical(mapped$rule$minimum_mapped, 10L)

    # A tie of absolute values keeps the first in input order.
    tied <- data.frame(feature_id = c("P11111", "P11111-3", accessions), value = c(2, -2, seq(1, 10)))
    expect_identical(integrationEnrichmentMapIdentifiers(tied, "proteomics")$table$feature_id[[1L]], "P11111")

    genes <- data.frame(feature_id = c("ENSG00000141510.17", "TP53", sprintf("GENE%02d", 1:8), " "),
        value = c(1:10, 11))
    gene_map <- integrationEnrichmentMapIdentifiers(genes, "transcriptomics")
    expect_identical(gene_map$table$identifier[1:2], c("ENSG00000141510", "TP53"))
    expect_identical(gene_map$counts$unmapped, 1L)

    features <- sprintf("m%02d", 1:12)
    metabolites <- data.frame(feature_id = features, value = seq(-6, 5))
    annotation <- data.frame(feature_id = rev(features), chebi = rev(c("CHEBI:15377", "15903", "chebi:0017234", "",
        "CHEBI:15377; CHEBI:1", sprintf("CHEBI:%d", 100:106))), kegg = rev(c("cpd:C00031", sprintf("C%05d", 1:10),
        "x")), stringsAsFactors = FALSE)
    chebi <- integrationEnrichmentMapIdentifiers(metabolites, "metabolomics", "chebi", annotation, "chebi")
    expect_identical(chebi$mapping$identifier[1:5], c("CHEBI:15377", "CHEBI:15903", "CHEBI:17234", NA,
        "CHEBI:15377"))
    expect_identical(chebi$counts[c("mapped", "unmapped", "duplicated")], list(mapped = 10L, unmapped = 1L,
        duplicated = 1L))
    expect_identical(chebi$rule$id_column, "chebi")
    kegg <- integrationEnrichmentMapIdentifiers(metabolites, "lipidomics", "kegg", annotation, "kegg")
    expect_identical(kegg$table$identifier[[1L]], "C00031")
    expect_identical(kegg$counts$mapped, 11L)

    expect_error(integrationEnrichmentMapIdentifiers(metabolites, "metabolomics", "chebi"),
        class = "integration_enrichment_parameters_invalid")
    expect_error(integrationEnrichmentMapIdentifiers(ranked, "proteomics", "chebi"),
        class = "integration_enrichment_parameters_invalid")
    expect_error(integrationEnrichmentMapIdentifiers(ranked, "glycomics"),
        class = "integration_enrichment_parameters_invalid")
})

test_that("fewer than 10 mapped identifiers is refused with the counts", {
    ranked <- data.frame(feature_id = c(sprintf("P%05d", 10001:10009), "unknown"), value = 1:10)
    condition <- expect_error(integrationEnrichmentMapIdentifiers(ranked, "proteomics"),
        class = "integration_enrichment_too_few_mapped")
    expect_s3_class(condition, "integration_refused")
    expect_identical(condition$counts$mapped, 9L)
    expect_identical(condition$counts$unmapped, 1L)
    expect_match(conditionMessage(condition), "9 of 10 features matched", fixed = TRUE)
    ranked$feature_id[[10L]] <- "P10010"
    expect_identical(integrationEnrichmentMapIdentifiers(ranked, "proteomics")$counts$mapped, 10L)
})

pathwayFixture <- function() {
    ids <- sprintf("G%02d", 1:60)
    ranked <- data.frame(identifier = ids, value = seq(3, -3, length.out = 60), stringsAsFactors = FALSE)
    term2gene <- data.frame(term = rep(c("top", "bottom", "spread", "middle"), c(10, 10, 10, 15)),
        gene = c(ids[1:10], ids[51:60], ids[seq(1, 60, by = 6)], ids[21:35]), stringsAsFactors = FALSE)
    term2name <- data.frame(term = c("top", "bottom", "spread", "middle"), name = c("Top set", "Bottom set",
        "Spread set", "Middle set"), stringsAsFactors = FALSE)
    list(ranked = ranked, term2gene = term2gene, term2name = term2name)
}

test_that("the pathway test is deterministic with its seed and leaves the caller's random state alone", {
    skip_if_not_installed("clusterProfiler")
    fixture <- pathwayFixture()
    requireNamespace("clusterProfiler", quietly = TRUE)
    set.seed(7)
    before <- .Random.seed
    first <- integrationEnrichmentGsea(fixture$ranked, fixture$term2gene, fixture$term2name, seed = 11L)
    expect_identical(.Random.seed, before)
    second <- integrationEnrichmentGsea(fixture$ranked, fixture$term2gene, fixture$term2name, seed = 11L)
    expect_identical(first, second)
    expect_identical(first$method, "gsea")
    expect_identical(first$seed, 11L)
    expect_identical(first$settings, integrationGseaSettings())
    expect_true(all(c("ID", "pvalue", "p.adjust") %in% names(first$terms)))
    expect_setequal(first$terms$ID, c("top", "bottom", "spread", "middle"))
    expect_equal(first$terms$p.adjust, stats::p.adjust(first$terms$pvalue, "BH"))

    jittered <- integrationGseaSettings(p_adjust_method = "bonferroni", tie_jitter = 1e-4)
    third <- integrationEnrichmentGsea(fixture$ranked, fixture$term2gene, fixture$term2name, 5L, jittered)
    expect_identical(third, integrationEnrichmentGsea(fixture$ranked, fixture$term2gene, fixture$term2name, 5L,
        jittered))
    expect_identical(third$settings$p_adjust_method, "bonferroni")
    expect_equal(third$terms$p.adjust, stats::p.adjust(third$terms$pvalue, "bonferroni"))

    expect_error(integrationEnrichmentGsea(fixture$ranked, fixture$term2gene, fixture$term2name, seed = 1.5),
        class = "integration_enrichment_parameters_invalid")
    expect_error(integrationEnrichmentGsea(fixture$ranked, fixture$term2gene, fixture$term2name, 1L,
        integrationGseaSettings(p_adjust_method = "made-up")), class = "integration_enrichment_parameters_invalid")
})

test_that("an error of the pathway test is refused and no other test is run in its place", {
    skip_if_not_installed("clusterProfiler")
    fixture <- pathwayFixture()
    enricher_calls <- 0L
    local_mocked_bindings(
        GSEA = function(...) stop("the engine broke"),
        enricher = function(...) {
            enricher_calls <<- enricher_calls + 1L
            stop("the over-representation test must not run")
        },
        .package = "clusterProfiler"
    )
    condition <- expect_error(integrationEnrichmentGsea(fixture$ranked, fixture$term2gene, fixture$term2name, 1L),
        class = "integration_enrichment_method_failed")
    expect_s3_class(condition, "integration_refused")
    expect_match(condition$engine_message, "the engine broke", fixed = TRUE)
    expect_identical(enricher_calls, 0L)
})

test_that("the legacy Reactome helper keeps its fallback to over-representation on an error of the test", {
    skip_if_not_installed("clusterProfiler")
    fixture <- pathwayFixture()
    data_dir <- withr::local_tempdir()
    chebi <- sprintf("CHEBI:%d", 1:60)
    utils::write.table(data.frame(chebi_id = chebi[match(fixture$term2gene$gene, fixture$ranked$identifier)],
        reactome_id = paste0("R-HSA-", fixture$term2gene$term), reactome_term = fixture$term2name$name[
            match(fixture$term2gene$term, fixture$term2name$term)], organism = "Homo sapiens"),
        file.path(data_dir, "chebi_to_reactome.tab"), sep = "\t", quote = FALSE, row.names = FALSE)
    local_mocked_bindings(getProjectPaths = function(omic_type, experiment_label) list(data_dir = data_dir))
    local_mocked_bindings(GSEA = function(...) stop("the engine broke"), .package = "clusterProfiler")
    expect_message(
        result <- runReactomeEnrichment(stats::setNames(fixture$ranked$value, chebi), mapping_table = NULL,
            project_dirs = list(), omic_type = "metabolomics", experiment_label = "demo", assay_name = "assay",
            reactome_organism = "Homo sapiens"),
        "Trying enricher instead", fixed = TRUE
    )
    expect_identical(names(result), c("termDescription", "enrichmentScore", "falseDiscoveryRate", "genesMapped",
        "mappedIDs", "comparison"))
    expect_true(nrow(result) >= 1L)
    expect_true("Top set" %in% result$termDescription)
    expect_equal(result$enrichmentScore, -log10(result$falseDiscoveryRate))
})

test_that("the KEGG kernel declares organism and release and returns its mapping tables", {
    skip_if_not_installed("clusterProfiler")
    fixture <- pathwayFixture()
    compounds <- sprintf("C%05d", 1:60)
    ranked <- data.frame(identifier = compounds, value = fixture$ranked$value, stringsAsFactors = FALSE)
    numbers <- c(top = "00010", bottom = "00020", spread = "00030", middle = "00040")
    gene <- compounds[match(fixture$term2gene$gene, fixture$ranked$identifier)]
    pathway_compounds <- data.frame(
        pathway = c(paste0("path:hsa", numbers[fixture$term2gene$term]), "path:mmu00010", "path:hsa00050"),
        compound = c(paste0("cpd:", gene), "cpd:C00001", "cpd:C99999"), stringsAsFactors = FALSE)
    pathway_names <- data.frame(pathway = paste0("path:map", numbers), name = fixture$term2name$name,
        stringsAsFactors = FALSE)
    result <- integrationKeggEnrichment(ranked, pathway_compounds, pathway_names, "hsa", "KEGG 116.0", 3L)
    expect_identical(result$database, "KEGG")
    expect_identical(result$organism, "hsa")
    expect_identical(result$release, "KEGG 116.0")
    expect_true(all(grepl("^path:hsa", result$term2gene$term)))
    expect_true(all(result$term2gene$gene %in% compounds))
    expect_false("path:hsa00050" %in% result$term2gene$term)
    expect_identical(sort(result$term2name$name), sort(fixture$term2name$name))
    expect_identical(result, integrationKeggEnrichment(ranked, pathway_compounds, pathway_names, "hsa",
        "KEGG 116.0", 3L))
    expect_error(integrationKeggEnrichment(ranked, pathway_compounds, pathway_names, "HSA", "KEGG 116.0", 3L),
        class = "integration_enrichment_parameters_invalid")
    expect_error(integrationKeggEnrichment(ranked, pathway_compounds, pathway_names, "hsa", "", 3L),
        class = "integration_enrichment_parameters_invalid")
    expect_error(integrationKeggEnrichment(ranked, pathway_compounds, pathway_names, "eco", "KEGG 116.0", 3L),
        class = "integration_enrichment_no_overlap")
})

test_that("the Reactome kernel resolves the organism and returns its mapping tables", {
    skip_if_not_installed("clusterProfiler")
    fixture <- pathwayFixture()
    chebi <- sprintf("CHEBI:%d", 1:60)
    ranked <- data.frame(identifier = chebi, value = fixture$ranked$value, stringsAsFactors = FALSE)
    gene <- chebi[match(fixture$term2gene$gene, fixture$ranked$identifier)]
    mapping <- data.frame(chebi_id = c(gene, gene[[1L]], "CHEBI:1"),
        reactome_id = c(paste0("R-HSA-", fixture$term2gene$term), "R-HSA-top", "R-MMU-1"),
        reactome_term = c(fixture$term2name$name[match(fixture$term2gene$term, fixture$term2name$term)], "Top set",
            "Mouse set"), organism = c(rep("Homo sapiens", nrow(fixture$term2gene) + 1L), "Mus musculus"),
        stringsAsFactors = FALSE)
    result <- integrationReactomeEnrichment(ranked, mapping, "homo sapiens", "Reactome 94", 3L)
    expect_identical(result$database, "Reactome")
    expect_identical(result$organism, "Homo sapiens")
    expect_identical(result$release, "Reactome 94")
    expect_identical(nrow(result$term2gene), nrow(fixture$term2gene))
    expect_false(anyDuplicated(result$term2gene) > 0L)
    expect_false("R-MMU-1" %in% result$term2gene$term)
    expect_identical(result, integrationReactomeEnrichment(ranked, mapping, "Homo sapiens", "Reactome 94", 3L))
    expect_error(integrationReactomeEnrichment(ranked, mapping, "Klebsiella pneumoniae", "Reactome 94", 3L),
        class = "multischolar_reactome_organism")
    no_overlap <- data.frame(identifier = sprintf("CHEBI:%d", 1001:1012), value = 1:12)
    expect_error(integrationReactomeEnrichment(no_overlap, mapping, "Homo sapiens", "Reactome 94", 3L),
        class = "integration_enrichment_no_overlap")
})
