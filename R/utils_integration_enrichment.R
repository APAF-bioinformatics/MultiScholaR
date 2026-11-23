# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# ============================================================================
# utils_integration_enrichment.R
# ============================================================================
# KIT-141, part 1: the enrichment kernels of the integration module. Each function takes a plain ranked table
# (identifier, value) and declared settings, and returns plain tables and lists: identifier mapping counted per omic,
# the STRING values/ranks client (httr2; the key only in the POST body, never kept), and the KEGG and Reactome
# pathway tests (seeded, no global state, no swap of method on error). No S4 object, no global variable, no file
# written. The integration module's operations, tab, record and replay call these (part 2).
# ============================================================================

# --- Constants ---

#' Fewer mapped identifiers than this is refused (KIT-141).
#' @noRd
INTEGRATION_ENRICHMENT_MIN_MAPPED <- 10L

#' The STRING release queried by default (INT-18): the stable address of version 12.5, the current release named by
#' STRING's API help (https://string-db.org/help/api/, read 2026-11-21); the version STRING reports is recorded with
#' every result.
#' @noRd
INTEGRATION_STRING_DEFAULT_HOST <- "https://version-12-5.string-db.org"

#' The refusals' sentences, by condition class. Plain words for people who are not data-science specialists; the
#' final wording is the user's (KIT-141 report). A refusal's message is its sentence, then a short detail in brackets.
#' @noRd
INTEGRATION_ENRICHMENT_REFUSALS <- c(
    integration_enrichment_parameters_invalid = "An enrichment setting is not valid.",
    integration_enrichment_key_missing = "A STRING key is needed. Enter your own STRING key and try again.",
    integration_enrichment_too_few_mapped = paste(
        "Too few features could be matched to a database identifier for an enrichment test:",
        "at least 10 are needed."
    ),
    integration_enrichment_no_overlap = paste(
        "None of the matched features is in a pathway of the chosen database and organism,",
        "so the pathway test cannot be run."
    ),
    integration_enrichment_service_unreachable = paste(
        "The STRING service could not be reached. Check the internet connection and try again;",
        "no result was kept."
    ),
    integration_enrichment_service_failed = paste(
        "The STRING service did not return a usable answer; no result was kept."
    ),
    integration_enrichment_service_timeout = paste(
        "The STRING service did not finish the analysis in the time allowed; no result was kept.",
        "Try again later."
    ),
    integration_enrichment_method_failed = paste(
        "The pathway test could not be run on these data; no other test was run in its place."
    ),
    # KIT-141 part 2.
    integration_enrichment_no_model = "Train a model at the Model tab before running an enrichment.",
    integration_enrichment_no_testable_set = paste(
        "No pathway has enough measured features for the test: a pathway is tested only when it holds between the",
        "smallest and the largest number of measured features set below."
    ),
    integration_enrichment_download_failed = paste(
        "The pathway tables could not be downloaded. Check the internet connection and try again; nothing was kept."
    ),
    integration_enrichment_release_changed = paste(
        "The database published a new release while its tables were being downloaded, so the tables could not be",
        "labelled with one release. Download them again."
    ),
    integration_enrichment_response_mismatch = paste(
        "The saved answer of the service was for another list of features or other settings, so it was not used.",
        "Run the enrichment again."
    ),
    # KIT-183: the correlation-adjusted test (%d and %s are filled by the refusal's `values`).
    integration_enrichment_too_few_complete = paste(
        "Too few features have a value in every biological sample for the correlation-adjusted test: at least 10 are",
        "needed. The pre-ranked test can use features with missing values."
    ),
    integration_enrichment_too_few_samples_for_test = paste(
        "The correlation-adjusted test needs at least %d biological samples in this data set."
    ),
    integration_enrichment_test_not_offered = "This test is not offered here: %s.",
    # KIT-191 (D4 (c); %s: minutes; %d: features, biological samples, terms, reorderings).
    integration_enrichment_test_too_large = paste(
        "The permutation test would take about %s minutes for this data set (%d features, %d biological samples, %d",
        "terms, %d reorderings), longer than the 2 minutes allowed in the tab, where the application waits for it.",
        "Choose a correlation-adjusted test or fewer reorderings, or run it from a script."
    )
)

# --- Conditions and the key ---

#' Refuse with one of the enrichment's sentences: classed `<class>`, `integration_enrichment_refused` and
#' `integration_refused`. The message is built here from the sentence and a short detail, never from a request,
#' and any trace of the key is removed from it; no parent condition is attached (an httr2 failure carries its
#' request, and with it the key).
#' @param class The refusal's class (a name of `INTEGRATION_ENRICHMENT_REFUSALS`).
#' @param detail A short detail, or NULL.
#' @param key The key to remove from the message, or NULL.
#' @param ... Further fields of the condition (counts, the engine's message).
#' @param values KIT-183: the values a sentence with `%d` or `%s` takes, or NULL.
#' @noRd
integrationEnrichmentRefuse <- function(class, detail = NULL, key = NULL, ..., values = NULL) {
    sentence <- INTEGRATION_ENRICHMENT_REFUSALS[[class]]
    if (length(values)) sentence <- do.call(sprintf, c(list(sentence), values))
    message <- if (length(detail) == 1L && !is.na(detail) && nzchar(detail)) {
        paste0(sentence, " (", detail, ")")
    } else {
        sentence
    }
    fields <- integrationRedactKey(list(...), key)
    stop(structure(
        class = c(class, "integration_enrichment_refused", "integration_refused", "error", "condition"),
        c(list(message = integrationRedactKey(message, key), call = NULL), fields)
    ))
}

#' Remove every occurrence of the key from a value: strings, names, list and data frame elements.
#' @param x A value.
#' @param key The key, or NULL.
#' @return `x` with the key replaced by `[removed]`.
#' @noRd
integrationRedactKey <- function(x, key) {
    if (is.null(key) || !is.character(key) || length(key) != 1L || is.na(key) || !nzchar(key)) return(x)
    # The key as written and as encoded in an address.
    forms <- unique(c(key, utils::URLencode(key, reserved = TRUE)))
    redact <- function(text) {
        for (form in forms) text <- gsub(form, "[removed]", text, fixed = TRUE)
        text
    }
    if (is.character(x)) {
        x[] <- redact(x)
    } else if (is.list(x)) {
        x[] <- lapply(x, integrationRedactKey, key = key)
    }
    if (!is.null(names(x))) names(x) <- redact(names(x))
    x
}

#' Forget httr2's record of the last request and response. httr2 keeps both in its own environment after a call
#' (`httr2::last_request()`); the request holds the key in its body, so it must not outlive the call.
#' @noRd
integrationForgetHttpRequest <- function() {
    the <- get0("the", envir = asNamespace("httr2"), inherits = FALSE)
    if (!is.environment(the)) return(invisible(NULL))
    for (name in c("last_request", "last_response")) {
        if (exists(name, envir = the, inherits = FALSE)) assign(name, NULL, envir = the)
    }
    invisible(NULL)
}

# --- Identifier mapping ---

#' The identifier types per omic, and the rule each applies.
#' @noRd
INTEGRATION_ENRICHMENT_ID_TYPES <- list(
    proteomics = "uniprot",
    transcriptomics = "gene",
    metabolomics = c("chebi", "kegg"),
    lipidomics = c("chebi", "kegg")
)

#' One identifier from one raw value: the first token (split at `;`, `|` or `,`) that is a valid identifier of the
#' type, normalised. UniProt: an accession, isoform suffix removed. Gene: the first non-empty token, an Ensembl
#' version suffix removed. ChEBI: `CHEBI:<number>`. KEGG: a compound `C<5 digits>`, `cpd:` removed.
#' @param raw Character vector of raw values.
#' @param id_type One of `uniprot`, `gene`, `chebi`, `kegg`.
#' @return Character vector, NA where no valid identifier is found.
#' @noRd
integrationEnrichmentParseIdentifier <- function(raw, id_type) {
    uniprot <- "^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
    one <- function(value) {
        if (is.na(value)) return(NA_character_)
        tokens <- trimws(strsplit(value, "[;|,]")[[1L]])
        tokens <- tokens[nzchar(tokens)]
        candidates <- switch(id_type,
            uniprot = sub("-[0-9]+$", "", tokens),
            gene = ifelse(grepl("^ENS[A-Z]*[0-9]+\\.[0-9]+$", tokens), sub("\\.[0-9]+$", "", tokens), tokens),
            chebi = ifelse(grepl("^(chebi:)?[0-9]+$", tokens, ignore.case = TRUE),
                paste0("CHEBI:", sub("^(chebi:)?0*([0-9])", "\\2", tokens, ignore.case = TRUE)), NA_character_),
            kegg = ifelse(grepl("^(cpd:)?C[0-9]{5}$", tokens), sub("^cpd:", "", tokens), NA_character_)
        )
        valid <- switch(id_type,
            uniprot = grepl(uniprot, candidates),
            gene = !is.na(candidates) & nzchar(candidates) & !grepl("[[:space:]]", candidates),
            !is.na(candidates)
        )
        if (any(valid)) candidates[valid][[1L]] else NA_character_
    }
    vapply(as.character(raw), one, character(1), USE.NAMES = FALSE)
}

#' Map a ranked table's features to database identifiers, by the declared rule of the omic, and count the outcome.
#'
#' Proteomics: the UniProt accession in the feature id. Transcriptomics: the gene id in the feature id. Metabolomics
#' and lipidomics: the user-chosen ChEBI or KEGG column of the feature annotation. A feature whose value is missing
#' is dropped; one without a valid identifier is unmapped; of several features with one identifier the one with the
#' largest absolute value is kept (the first in input order on a tie), the others are counted as duplicated.
#' @param ranked Data frame with `feature_id` (character) and `value` (numeric).
#' @param omic One of `proteomics`, `transcriptomics`, `metabolomics`, `lipidomics`.
#' @param id_type `uniprot` (proteomics), `gene` (transcriptomics), `chebi` or `kegg` (metabolomics, lipidomics).
#' @param annotation For metabolomics and lipidomics: a data frame with `feature_id` and the identifier column.
#' @param id_column For metabolomics and lipidomics: the identifier column of `annotation`.
#' @return `list(table, mapping, counts, rule)`: `table` the mapped features (`identifier`, `value`, `feature_id`),
#'   in input order; `mapping` every feature with its identifier and status (`mapped`, `unmapped`, `duplicated`,
#'   `missing_value`); `counts` (features, mapped, unmapped, duplicated, missing_value); `rule` the declared rule.
#'   Fewer than 10 mapped identifiers is refused (`integration_enrichment_too_few_mapped`, with the counts).
#' @noRd
integrationEnrichmentMapIdentifiers <- function(ranked, omic, id_type = INTEGRATION_ENRICHMENT_ID_TYPES[[omic]][[1L]],
    annotation = NULL, id_column = NULL) {
    invalid <- function(detail) integrationEnrichmentRefuse("integration_enrichment_parameters_invalid", detail)
    if (!is.character(omic) || length(omic) != 1L || !omic %in% names(INTEGRATION_ENRICHMENT_ID_TYPES)) {
        invalid("the omic must be proteomics, transcriptomics, metabolomics or lipidomics")
    }
    if (!is.character(id_type) || length(id_type) != 1L || !id_type %in% INTEGRATION_ENRICHMENT_ID_TYPES[[omic]]) {
        invalid(sprintf("identifier type for %s must be one of: %s", omic,
            paste(INTEGRATION_ENRICHMENT_ID_TYPES[[omic]], collapse = ", ")))
    }
    if (!is.data.frame(ranked) || !all(c("feature_id", "value") %in% names(ranked)) || !is.numeric(ranked$value)) {
        invalid("the ranked table needs a feature_id column and a numeric value column")
    }
    feature_id <- as.character(ranked$feature_id)
    if (anyNA(feature_id) || anyDuplicated(feature_id)) invalid("feature ids must be present and unique")

    raw <- feature_id
    if (omic %in% c("metabolomics", "lipidomics")) {
        if (!is.data.frame(annotation) || !is.character(id_column) || length(id_column) != 1L ||
            !all(c("feature_id", id_column) %in% names(annotation))) {
            invalid("choose the annotation column that holds the ChEBI or KEGG identifiers")
        }
        raw <- as.character(annotation[[id_column]])[match(feature_id, as.character(annotation$feature_id))]
    }

    identifier <- integrationEnrichmentParseIdentifier(raw, id_type)
    value <- as.numeric(ranked$value)
    status <- ifelse(is.na(value), "missing_value", ifelse(is.na(identifier), "unmapped", "mapped"))
    candidate <- which(status == "mapped")
    keep_order <- candidate[order(-abs(value[candidate]), candidate)]
    status[keep_order[duplicated(identifier[keep_order])]] <- "duplicated"

    mapping <- data.frame(feature_id = feature_id, identifier = identifier, value = value, status = status,
        stringsAsFactors = FALSE)
    kept <- status == "mapped"
    table <- data.frame(identifier = identifier[kept], value = value[kept], feature_id = feature_id[kept],
        stringsAsFactors = FALSE)
    counts <- list(features = length(feature_id), mapped = sum(kept), unmapped = sum(status == "unmapped"),
        duplicated = sum(status == "duplicated"), missing_value = sum(status == "missing_value"))
    rule <- list(omic = omic, id_type = id_type, id_column = if (is.null(id_column)) NA_character_ else id_column,
        duplicates = "largest absolute value kept, first in input order on a tie",
        minimum_mapped = INTEGRATION_ENRICHMENT_MIN_MAPPED)
    if (counts$mapped < INTEGRATION_ENRICHMENT_MIN_MAPPED) {
        integrationEnrichmentRefuse("integration_enrichment_too_few_mapped",
            sprintf("%d of %d features matched", counts$mapped, counts$features), counts = counts, rule = rule)
    }
    list(table = table, mapping = mapping, counts = counts, rule = rule)
}

#' Check a mapped ranked table (`identifier`, `value`): unique present identifiers, numeric values without NA.
#' @noRd
integrationEnrichmentCheckRanked <- function(ranked) {
    ok <- is.data.frame(ranked) && nrow(ranked) > 0L && all(c("identifier", "value") %in% names(ranked)) &&
        is.character(ranked$identifier) && !anyNA(ranked$identifier) && !anyDuplicated(ranked$identifier) &&
        is.numeric(ranked$value) && !anyNA(ranked$value)
    if (!ok) {
        integrationEnrichmentRefuse("integration_enrichment_parameters_invalid",
            "the ranked table needs unique identifiers and a numeric value for each")
    }
    invisible(TRUE)
}

# --- STRING values/ranks enrichment ---

#' The STRING client's declared settings.
#' @param ge_fdr STRING's false discovery rate threshold for reported terms.
#' @param rank_direction STRING's `ge_enrichment_rank_direction`: -1, 0 or 1.
#' @param host The STRING release's stable address (its version is queried and recorded).
#' @param poll_interval_seconds Seconds between status queries.
#' @param max_polls Status queries before the job counts as not finished.
#' @param timeout_seconds Time allowed for one HTTP request.
#' @param caller_identity STRING's `caller_identity`.
#' @return A named list.
#' @noRd
integrationStringOptions <- function(ge_fdr = 0.05, rank_direction = -1L, host = INTEGRATION_STRING_DEFAULT_HOST,
    poll_interval_seconds = 10, max_polls = 30L, timeout_seconds = 60, caller_identity = "MultiScholaR") {
    list(ge_fdr = ge_fdr, rank_direction = rank_direction, host = host, poll_interval_seconds = poll_interval_seconds,
        max_polls = max_polls, timeout_seconds = timeout_seconds, caller_identity = caller_identity)
}

#' Check the STRING settings; refuse the first invalid one.
#' @noRd
integrationStringCheckOptions <- function(species, options) {
    invalid <- function(detail) integrationEnrichmentRefuse("integration_enrichment_parameters_invalid", detail)
    number <- function(x, lower, upper) is.numeric(x) && length(x) == 1L && !is.na(x) && x >= lower && x <= upper
    if (!(is.character(species) || is.numeric(species)) || length(species) != 1L || is.na(species) ||
        !grepl("^[0-9]+$", as.character(species))) {
        invalid("the species must be a STRING species number, such as 9606")
    }
    if (!number(options$ge_fdr, 0, 1)) invalid("the STRING false discovery rate must be between 0 and 1")
    if (!number(options$rank_direction, -1, 1) || options$rank_direction != round(options$rank_direction)) {
        invalid("the STRING rank direction must be -1, 0 or 1")
    }
    if (!is.character(options$host) || length(options$host) != 1L || !grepl("^https://[^/?#]+$", options$host)) {
        invalid("the STRING address must be an https address without a path")
    }
    if (!number(options$poll_interval_seconds, 0, 600)) invalid("the waiting time must be between 0 and 600 seconds")
    if (!number(options$max_polls, 1, 1000)) invalid("the number of status queries must be between 1 and 1000")
    if (!number(options$timeout_seconds, 1, 3600)) invalid("the request time limit must be between 1 and 3600 seconds")
    if (!is.character(options$caller_identity) || length(options$caller_identity) != 1L ||
        !nzchar(options$caller_identity)) {
        invalid("the caller identity must be a short text")
    }
    invisible(TRUE)
}

#' One request to STRING: POST with the fields in a form body (or GET with no body), a user agent, a time limit and
#' retries of a busy service. Every failure is a classed refusal with no parent condition; httr2's record of the
#' request is forgotten on exit.
#' @param url The address (never carrying the key).
#' @param fields Named list of form fields, or NULL for a GET.
#' @param options `integrationStringOptions()`.
#' @param key The key, to keep it out of any message.
#' @return `list(status, body)`: the HTTP status and the body as text.
#' @noRd
integrationStringRequest <- function(url, fields, options, key) {
    on.exit(integrationForgetHttpRequest(), add = TRUE)
    request <- httr2::request(url) |>
        httr2::req_user_agent("MultiScholaR") |>
        httr2::req_timeout(options$timeout_seconds) |>
        httr2::req_retry(max_tries = 3L, backoff = \(attempt) options$poll_interval_seconds) |>
        httr2::req_error(is_error = \(response) FALSE)
    if (!is.null(fields)) request <- do.call(httr2::req_body_form, c(list(request), fields))
    response <- tryCatch(httr2::req_perform(request), error = \(e) NULL)
    rm(request)
    if (is.null(response)) {
        integrationEnrichmentRefuse("integration_enrichment_service_unreachable", key = key)
    }
    status <- httr2::resp_status(response)
    body <- tryCatch(httr2::resp_body_string(response), error = \(e) "")
    rm(response)
    if (status >= 400L) {
        integrationEnrichmentRefuse("integration_enrichment_service_failed", sprintf("HTTP status %d", status),
            key = key)
    }
    list(status = status, body = integrationRedactKey(body, key))
}

#' Parse a STRING JSON reply to its first record (STRING answers with an array of one object).
#' @noRd
integrationStringRecord <- function(body, key, what) {
    # The reply is read before the parse's error handler is set: `body` is the request itself, passed unevaluated, and a
    # refusal raised by the request (the service unreachable) must not be taken for an unreadable reply (run-tests
    # 44783719082 and 44786840568: a refused connection was reported as "the reply could not be read").
    force(body)
    parsed <- tryCatch(jsonlite::fromJSON(body, simplifyVector = FALSE), error = \(e) NULL)
    if (is.list(parsed) && is.null(names(parsed)) && length(parsed) >= 1L) parsed <- parsed[[1L]]
    if (!is.list(parsed) || is.null(names(parsed))) {
        integrationEnrichmentRefuse("integration_enrichment_service_failed",
            sprintf("the %s reply could not be read", what), key = key)
    }
    integrationRedactKey(parsed, key)
}

#' A field of a STRING record as one string, or NA.
#' @noRd
integrationStringField <- function(record, name) {
    value <- record[[name]]
    if (is.null(value) || length(value) != 1L || is.na(value)) NA_character_ else as.character(value)
}

#' Fetch the results table at the address STRING gave. An `api_key` in that address is taken out of it and sent in
#' the POST body instead, so the key never travels in an address.
#' @noRd
integrationStringDownload <- function(download_url, api_key, options) {
    if (is.na(download_url) || !startsWith(download_url, "https://")) {
        integrationEnrichmentRefuse("integration_enrichment_service_failed",
            "STRING gave no https address for the results", key = api_key)
    }
    parsed <- httr2::url_parse(download_url)
    key_in_url <- "api_key" %in% names(parsed$query)
    parsed$query$api_key <- NULL
    fields <- if (key_in_url) list(api_key = api_key) else NULL
    reply <- integrationStringRequest(httr2::url_build(parsed), fields, options, api_key)
    terms <- integrationStringTermsTable(reply$body)
    if (!is.data.frame(terms)) {
        integrationEnrichmentRefuse("integration_enrichment_service_failed", "the results table could not be read",
            key = api_key)
    }
    list(text = reply$body, terms = integrationRedactKey(terms, api_key))
}

#' STRING values/ranks enrichment of a mapped ranked table (KIT-141).
#'
#' Queries the release's version, submits the identifiers and values, polls the job and fetches the results table.
#' The key is used only in the POST bodies of this call: never in an address, a setting, a message or the returned
#' value, and httr2's record of the last request is cleared. Any network or service failure is a classed refusal and
#' nothing is returned.
#' @param ranked Data frame with `identifier` and `value` (`integrationEnrichmentMapIdentifiers()$table`).
#' @param api_key The user's STRING key.
#' @param species STRING species number, such as 9606.
#' @param options `integrationStringOptions()`.
#' @return `list(service, method, version, request, response, terms)`: `version` the release STRING reports
#'   (`string_version`, `string_stable_address`); `request` a summary of what was sent (no key: the settings, the
#'   number of identifiers and the SHA-256 of the identifier/value text); `response` the submission reply, the last
#'   status reply and the results table as text; `terms` the results table.
#' @noRd
integrationStringRankedEnrichment <- function(ranked, api_key, species, options = integrationStringOptions()) {
    if (!is.character(api_key) || length(api_key) != 1L || is.na(api_key) || !nzchar(trimws(api_key))) {
        integrationEnrichmentRefuse("integration_enrichment_key_missing")
    }
    requireOptionalPackage("httr2", "STRING enrichment")
    integrationEnrichmentCheckRanked(ranked)
    integrationStringCheckOptions(species, options)
    species <- as.character(species)
    api <- paste0(options$host, "/api/json/")

    # --- The release ---
    version <- integrationStringRecord(integrationStringRequest(paste0(api, "version"), NULL, options,
        api_key)$body, api_key, "version")
    version <- list(string_version = integrationStringField(version, "string_version"),
        string_stable_address = integrationStringField(version, "string_stable_address"))
    if (is.na(version$string_version)) {
        integrationEnrichmentRefuse("integration_enrichment_service_failed", "STRING did not state its version",
            key = api_key)
    }

    # --- Submission ---
    identifiers <- integrationStringIdentifiersText(ranked)
    submission <- integrationStringRecord(integrationStringRequest(paste0(api, "valuesranks_enrichment_submit"),
        list(species = species, caller_identity = options$caller_identity, identifiers = identifiers,
            api_key = api_key, ge_fdr = options$ge_fdr,
            ge_enrichment_rank_direction = as.integer(options$rank_direction)),
        options, api_key)$body, api_key, "submission")
    job_id <- integrationStringField(submission, "job_id")
    if (is.na(job_id) || identical(integrationStringField(submission, "status"), "error")) {
        integrationEnrichmentRefuse("integration_enrichment_service_failed",
            integrationStringDetail(submission, "the job was not accepted"), key = api_key)
    }

    # --- Polling ---
    status <- NULL
    polls <- 0L
    repeat {
        if (polls >= options$max_polls) {
            integrationEnrichmentRefuse("integration_enrichment_service_timeout",
                sprintf("%d status queries", polls), key = api_key)
        }
        Sys.sleep(options$poll_interval_seconds)
        polls <- polls + 1L
        status <- integrationStringRecord(integrationStringRequest(paste0(api, "valuesranks_enrichment_status"),
            list(api_key = api_key, job_id = job_id), options, api_key)$body, api_key, "status")
        state <- integrationStringField(status, "status")
        if (identical(state, "success")) break
        if (!is.na(state) && state %in% c("failed", "error")) {
            integrationEnrichmentRefuse("integration_enrichment_service_failed",
                integrationStringDetail(status, "the job failed"), key = api_key)
        }
    }

    # --- Results ---
    results <- integrationStringDownload(integrationStringField(status, "download_url"), api_key, options)
    request <- list(host = options$host, method = "valuesranks_enrichment", species = species,
        ge_fdr = options$ge_fdr, rank_direction = as.integer(options$rank_direction),
        caller_identity = options$caller_identity, identifiers = nrow(ranked),
        identifiers_sha256 = digest::digest(identifiers, algo = "sha256", serialize = FALSE),
        poll_interval_seconds = options$poll_interval_seconds, max_polls = options$max_polls, polls = polls,
        queried_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3%z"))
    integrationRedactKey(list(service = "STRING", method = "valuesranks_enrichment", version = version,
        request = request, response = list(submission = submission, status = status, results_tsv = results$text),
        terms = results$terms), api_key)
}

#' The identifiers and values as STRING receives them: one `identifier<TAB>value` line per feature, the values written
#' with R's `as.character()` (15 significant digits). Its SHA-256 is the request's identity.
#' @param ranked Data frame with `identifier` and `value`.
#' @noRd
integrationStringIdentifiersText <- function(ranked) {
    paste(ranked$identifier, as.character(ranked$value), sep = "\t", collapse = "\n")
}

#' STRING's results table from its text (tab-separated, with a header), or NULL when it cannot be read.
#' @noRd
integrationStringTermsTable <- function(text) {
    tryCatch(utils::read.delim(text = text, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE,
        comment.char = ""), error = \(e) NULL)
}

#' STRING's own message of a record (at most 200 characters), or `fallback` when it gives none.
#' @noRd
integrationStringDetail <- function(record, fallback) {
    detail <- integrationStringField(record, "message")
    if (is.na(detail) || !nzchar(detail)) fallback else substr(detail, 1L, 200L)
}

# --- Pathway tests (KEGG, Reactome) ---

#' The pathway test's declared settings.
#' @param p_adjust_method Multiple-testing correction, one of `stats::p.adjust.methods`.
#' @param min_size,max_size Smallest and largest pathway size tested (features of the ranked table in it).
#' @param p_cutoff Pathways with a raw p-value above this are left out of the table (1 keeps every pathway).
#' @param tie_jitter Half-width of the uniform noise added to break tied values (0: values used as they are).
#' @return A named list.
#' @noRd
integrationGseaSettings <- function(p_adjust_method = "BH", min_size = 3L, max_size = 500L, p_cutoff = 1,
    tie_jitter = 0) {
    list(p_adjust_method = p_adjust_method, min_size = min_size, max_size = max_size, p_cutoff = p_cutoff,
        tie_jitter = tie_jitter)
}

#' Check the pathway test's seed and settings; refuse the first invalid one.
#' @noRd
integrationGseaCheckSettings <- function(seed, settings) {
    invalid <- function(detail) integrationEnrichmentRefuse("integration_enrichment_parameters_invalid", detail)
    whole <- function(x) is.numeric(x) && length(x) == 1L && !is.na(x) && x == round(x)
    if (!whole(seed)) invalid("the seed must be a whole number")
    if (!is.character(settings$p_adjust_method) || length(settings$p_adjust_method) != 1L ||
        !settings$p_adjust_method %in% stats::p.adjust.methods) {
        invalid("the multiple-testing correction is not known")
    }
    if (!whole(settings$min_size) || !whole(settings$max_size) || settings$min_size < 1 ||
        settings$max_size < settings$min_size) {
        invalid("the pathway sizes must be whole numbers, the smallest at least 1 and not above the largest")
    }
    if (!is.numeric(settings$p_cutoff) || length(settings$p_cutoff) != 1L || is.na(settings$p_cutoff) ||
        settings$p_cutoff <= 0 || settings$p_cutoff > 1) {
        invalid("the p-value cut-off must be above 0 and at most 1")
    }
    if (!is.numeric(settings$tie_jitter) || length(settings$tie_jitter) != 1L || is.na(settings$tie_jitter) ||
        settings$tie_jitter < 0) {
        invalid("the tie-breaking noise must be 0 or more")
    }
    invisible(TRUE)
}

#' The pathway table returned when no pathway passes: no rows, the engine's main columns.
#' @noRd
integrationEnrichmentEmptyTerms <- function() {
    data.frame(ID = character(), Description = character(), setSize = integer(), enrichmentScore = numeric(),
        NES = numeric(), pvalue = numeric(), p.adjust = numeric(), qvalue = numeric(), rank = numeric(),
        leading_edge = character(), core_enrichment = character(), stringsAsFactors = FALSE)
}

#' The pathway test: ranked gene-set enrichment (clusterProfiler `GSEA`) of a ranked table against pathway sets.
#'
#' Seeded: the values are sorted, tied values optionally broken with uniform noise, and the test run inside
#' `withr::with_seed(seed)` with R's default generator kinds (KIT-191 D9: Mersenne-Twister, Inversion, Rejection,
#' whatever the session's `RNGkind()`), so the caller's random state is untouched and the same input, seed and
#' settings give the same table. An error of the test is a classed refusal (`integration_enrichment_method_failed`, with the engine's
#' message in `engine_message`): no other test is run in its place.
#' @param ranked Data frame with `identifier` and `value`.
#' @param term2gene Data frame with `term` and `gene` (gene = an identifier of `ranked`).
#' @param term2name Data frame with `term` and `name` (may have no rows).
#' @param seed Whole number.
#' @param settings `integrationGseaSettings()`.
#' @return `list(method, engine, seed, settings, ranked, terms)`; `terms` the engine's table (no rows when no
#'   pathway passes the size limits and the cut-off), with raw and adjusted p-values.
#' @noRd
integrationEnrichmentGsea <- function(ranked, term2gene, term2name, seed, settings = integrationGseaSettings()) {
    requireOptionalPackage("clusterProfiler", "pathway enrichment")
    integrationEnrichmentCheckRanked(ranked)
    integrationGseaCheckSettings(seed, settings)
    if (!is.data.frame(term2gene) || !all(c("term", "gene") %in% names(term2gene)) ||
        !is.data.frame(term2name) || !all(c("term", "name") %in% names(term2name))) {
        integrationEnrichmentRefuse("integration_enrichment_parameters_invalid",
            "the pathway tables need term and gene, and term and name, columns")
    }
    values <- stats::setNames(ranked$value, ranked$identifier)
    outcome <- withr::with_seed(seed, {
        values <- sort(values, decreasing = TRUE)
        if (settings$tie_jitter > 0) {
            values <- values + stats::runif(length(values), min = -settings$tie_jitter, max = settings$tie_jitter)
            values <- sort(values, decreasing = TRUE)
        }
        tryCatch(
            clusterProfiler::GSEA(geneList = values, TERM2GENE = term2gene[, c("term", "gene")],
                TERM2NAME = term2name[, c("term", "name")],
                minGSSize = settings$min_size, maxGSSize = settings$max_size, pvalueCutoff = settings$p_cutoff,
                pAdjustMethod = settings$p_adjust_method, verbose = FALSE),
            error = \(e) e
        )
    }, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
    if (inherits(outcome, "error")) {
        integrationEnrichmentRefuse("integration_enrichment_method_failed", conditionMessage(outcome),
            engine_message = conditionMessage(outcome))
    }
    terms <- if (is.null(outcome)) integrationEnrichmentEmptyTerms() else methods::slot(outcome, "result")
    rownames(terms) <- NULL
    engine <- list(clusterProfiler = as.character(utils::packageVersion("clusterProfiler")),
        enrichit = tryCatch(as.character(utils::packageVersion("enrichit")), error = \(e) NA_character_))
    list(method = "gsea", engine = engine, seed = seed, settings = settings,
        ranked = data.frame(identifier = ranked$identifier, value = ranked$value, stringsAsFactors = FALSE),
        terms = terms)
}

#' Check a release identifier: one non-empty string.
#' @noRd
integrationEnrichmentCheckRelease <- function(release) {
    if (!is.character(release) || length(release) != 1L || is.na(release) || !nzchar(trimws(release))) {
        integrationEnrichmentRefuse("integration_enrichment_parameters_invalid",
            "the release of the pathway tables must be stated")
    }
    invisible(TRUE)
}

#' The pathway test over pathway sets restricted to the tested identifiers; refused when none overlaps. The entry
#' point a database's sets are passed to (KIT-183; KIT-181 calls it per Gene Ontology aspect).
#' @param ranked Data frame with `identifier` and `value` (the pre-ranked test's list).
#' @param term2gene,term2name The database's sets and names.
#' @param seed The seed of the pre-ranked and (KIT-191) permutation tests (the correlation-adjusted tests use none).
#' @param settings `integrationGseaSettings()`; the correlated tests read `min_size`, `max_size` and, for the
#'   simulation's sensitivity arm of `camera` only, `allow_neg_cor` (absent: TRUE, VAL-53 (2)); KIT-191: the permutation
#'   test reads `permutations` (absent: 1000) and `time_limit` (seconds; absent: no limit, a message instead).
#' @param test `gsea` (the pre-ranked test, unchanged), `camera` or `camera_pr` (KIT-183 5.1, 5.3); KIT-191:
#'   `camera_floor` or `permutation`.
#' @param values For the correlated tests: `integrationEnrichmentTestValues()`'s value; the sets are restricted
#'   to its identifiers (U_cam).
#' @return The test's list with `term2gene` and `term2name` as tested.
#' @noRd
integrationEnrichmentPathwayTest <- function(ranked, term2gene, term2name, seed, settings, test = "gsea", values = NULL) {
    if (identical(test, "gsea")) {
        term2gene <- unique(term2gene[term2gene$gene %in% ranked$identifier, c("term", "gene"), drop = FALSE])
        rownames(term2gene) <- NULL
        if (nrow(term2gene) == 0L) integrationEnrichmentRefuse("integration_enrichment_no_overlap")
        term2name <- unique(term2name[term2name$term %in% term2gene$term, c("term", "name"), drop = FALSE])
        rownames(term2name) <- NULL
        test <- integrationEnrichmentGsea(ranked, term2gene, term2name, seed, settings)
        return(c(test, list(term2gene = term2gene, term2name = term2name)))
    }
    if (!is.character(test) || length(test) != 1L || !test %in% INTEGRATION_ENRICHMENT_TESTS ||
        !is.list(values) || !is.matrix(values$y)) {
        integrationEnrichmentRefuse("integration_enrichment_parameters_invalid", paste("the test must be camera,",
            "camera_floor, permutation, camera_pr or gsea, and a correlated test needs the values of the data set"))
    }
    term2gene <- unique(term2gene[term2gene$gene %in% rownames(values$y), c("term", "gene"), drop = FALSE])
    rownames(term2gene) <- NULL
    if (nrow(term2gene) == 0L) integrationEnrichmentRefuse("integration_enrichment_no_overlap")
    term2name <- unique(term2name[term2name$term %in% term2gene$term, c("term", "name"), drop = FALSE])
    rownames(term2name) <- NULL
    result <- if (identical(test, "permutation")) {
        integrationEnrichmentPermutation(values, term2gene, term2name, settings, seed)
    } else {
        integrationEnrichmentCamera(values, term2gene, term2name, settings, test)
    }
    c(result, list(term2gene = term2gene, term2name = term2name))
}

#' KEGG pathway enrichment of a ranked table of KEGG compound identifiers (KIT-141).
#'
#' The caller supplies KEGG's tables and their release, so a replay can recompute offline from the returned ones.
#' @param ranked Data frame with `identifier` (KEGG compound, `C00031`) and `value`.
#' @param pathway_compounds Data frame with `pathway` (`path:<organism><5 digits>`) and `compound`
#'   (`cpd:C00031` or `C00031`).
#' @param pathway_names Data frame with `pathway` (`path:map00010`, `map00010` or the organism's id) and `name`.
#' @param organism KEGG organism code, such as `hsa`.
#' @param release KEGG release identifier of the tables.
#' @param seed Whole number.
#' @param settings `integrationGseaSettings()`.
#' @param test,values KIT-183: the test and, for a correlated test, the data set's values
#'   (`integrationEnrichmentPathwayTest()`).
#' @return The pathway test's list plus `database`, `organism`, `release`, `term2gene`, `term2name`.
#' @noRd
integrationKeggEnrichment <- function(ranked, pathway_compounds, pathway_names, organism, release, seed,
    settings = integrationGseaSettings(), test = "gsea", values = NULL) {
    integrationEnrichmentCheckRanked(ranked)
    integrationEnrichmentCheckRelease(release)
    if (!is.character(organism) || length(organism) != 1L || !grepl("^[a-z]{3,4}$", organism)) {
        integrationEnrichmentRefuse("integration_enrichment_parameters_invalid",
            "the KEGG organism must be a KEGG organism code, such as hsa")
    }
    if (!is.data.frame(pathway_compounds) || !all(c("pathway", "compound") %in% names(pathway_compounds)) ||
        !is.data.frame(pathway_names) || !all(c("pathway", "name") %in% names(pathway_names))) {
        integrationEnrichmentRefuse("integration_enrichment_parameters_invalid",
            "the KEGG tables need pathway and compound, and pathway and name, columns")
    }
    pathway <- sub("^(path:)?", "path:", as.character(pathway_compounds$pathway))
    in_organism <- grepl(paste0("^path:", organism, "[0-9]{5}$"), pathway)
    term2gene <- data.frame(term = pathway[in_organism],
        gene = sub("^cpd:", "", as.character(pathway_compounds$compound[in_organism])), stringsAsFactors = FALSE)
    number <- function(id) regmatches(id, regexpr("[0-9]{5}$", id))
    named <- grepl("[0-9]{5}$", as.character(pathway_names$pathway))
    names_by_number <- stats::setNames(as.character(pathway_names$name[named]),
        number(as.character(pathway_names$pathway[named])))
    terms <- unique(term2gene$term)
    term2name <- data.frame(term = terms, name = unname(names_by_number[number(terms)]), stringsAsFactors = FALSE)
    term2name$name[is.na(term2name$name)] <- term2name$term[is.na(term2name$name)]
    test <- integrationEnrichmentPathwayTest(ranked, term2gene, term2name, seed, settings, test, values)
    c(list(database = "KEGG", organism = organism, release = release), test)
}

#' Reactome pathway enrichment of a ranked table of ChEBI identifiers (KIT-141).
#'
#' The organism is matched in full against the mapping's organisms, ignoring case (`resolveReactomeOrganism`, which
#' refuses an unknown or ambiguous name with class `multischolar_reactome_organism`).
#' @param ranked Data frame with `identifier` (`CHEBI:<number>`) and `value`.
#' @param chebi_to_reactome Data frame with `chebi_id`, `reactome_id`, `reactome_term`, `organism` (Reactome's
#'   ChEBI-to-pathway file).
#' @param organism Organism name, such as `Homo sapiens`.
#' @param release Reactome release identifier of the mapping.
#' @param seed Whole number.
#' @param settings `integrationGseaSettings()`.
#' @param test,values KIT-183: the test and, for a correlated test, the data set's values
#'   (`integrationEnrichmentPathwayTest()`).
#' @return The pathway test's list plus `database`, `organism`, `release`, `term2gene`, `term2name`.
#' @noRd
integrationReactomeEnrichment <- function(ranked, chebi_to_reactome, organism, release, seed,
    settings = integrationGseaSettings(), test = "gsea", values = NULL) {
    integrationEnrichmentCheckRanked(ranked)
    integrationEnrichmentCheckRelease(release)
    columns <- c("chebi_id", "reactome_id", "reactome_term", "organism")
    if (!is.data.frame(chebi_to_reactome) || !all(columns %in% names(chebi_to_reactome))) {
        integrationEnrichmentRefuse("integration_enrichment_parameters_invalid",
            "the Reactome table needs chebi_id, reactome_id, reactome_term and organism columns")
    }
    organism <- resolveReactomeOrganism(organism, as.character(chebi_to_reactome$organism))
    rows <- chebi_to_reactome[as.character(chebi_to_reactome$organism) == organism, columns, drop = FALSE]
    chebi_number <- regmatches(as.character(rows$chebi_id), regexpr("[0-9]+", as.character(rows$chebi_id)))
    has_number <- grepl("[0-9]", as.character(rows$chebi_id))
    rows <- rows[has_number, , drop = FALSE]
    term2gene <- data.frame(term = as.character(rows$reactome_id), gene = paste0("CHEBI:", as.integer(chebi_number)),
        stringsAsFactors = FALSE)
    term2name <- data.frame(term = as.character(rows$reactome_id), name = as.character(rows$reactome_term),
        stringsAsFactors = FALSE)
    term2name <- term2name[!duplicated(term2name$term), , drop = FALSE]
    test <- integrationEnrichmentPathwayTest(ranked, term2gene, term2name, seed, settings, test, values)
    c(list(database = "Reactome", organism = organism, release = release), test)
}

# --- The legacy helpers' pathway step ---

#' The pathway step of the exported `runKeggEnrichment` and `runReactomeEnrichment`, kept as they behaved: the
#' kernel's test with their settings (seed 42, uniform tie-breaking noise of 1e-4, cut-off 0.05, `fdr`), and on an
#' error of the test their own over-representation test of the positively ranked features, as before. The kernel
#' itself never swaps the test; only these legacy functions keep their fallback.
#' @param ranked_sorted The ranked vector before the tie-breaking noise, sorted.
#' @param ranked_jittered The same vector after the legacy function's own noise (the fallback's input).
#' @param term2gene,term2name The legacy function's pathway tables.
#' @param label The legacy function's name, for its messages.
#' @return The kernel's pathway table (a data frame), the fallback's `enrichResult`, or NULL.
#' @noRd
integrationLegacyGseaStep <- function(ranked_sorted, ranked_jittered, term2gene, term2name, label) {
    test <- tryCatch(
        integrationEnrichmentGsea(
            data.frame(identifier = names(ranked_sorted), value = unname(ranked_sorted), stringsAsFactors = FALSE),
            term2gene, term2name, seed = 42L,
            settings = integrationGseaSettings(p_adjust_method = "fdr", p_cutoff = 0.05, tie_jitter = 0.0001)),
        # As before, a test that cannot run (also for want of clusterProfiler) goes to the fallback.
        integration_enrichment_method_failed = \(e) e,
        multischolar_package_missing = \(e) e
    )
    if (!inherits(test, "error")) return(test$terms)
    message(paste("GSEA error:", test$engine_message %||% conditionMessage(test)))
    message("Trying enricher instead...")
    positive_genes <- names(ranked_jittered)[ranked_jittered > 0]
    message(sprintf("   %s: Found %d positively ranked genes", label, length(positive_genes)))
    if (length(positive_genes) == 0L) return(NULL)
    tryCatch({
        enricher_result <- clusterProfiler::enricher(gene = positive_genes, universe = names(ranked_jittered),
            TERM2GENE = term2gene, TERM2NAME = term2name, pvalueCutoff = 0.05, pAdjustMethod = "fdr",
            minGSSize = 3, maxGSSize = 500)
        if (!is.null(enricher_result) && nrow(enricher_result@result) > 0) {
            message(sprintf("   %s: Enricher analysis successful with %d results", label,
                nrow(enricher_result@result)))
            enricher_result
        } else {
            message(sprintf("   %s: Enricher analysis returned no significant results", label))
            NULL
        }
    }, error = function(e2) {
        message(paste("Enricher error:", e2$message))
        NULL
    })
}

# ============================================================================
# KIT-141 part 2: enrichment of a trained model's factor (INT-18, INT-19). Two declared operations on the inputs state,
# run in the app process: the ranking of every measured feature of a data set by its correlation with the factor
# (`integration_enrichment_ranking`), and the enrichment of that ranked list (`integration_enrichment`). The STRING
# query runs between them, in the app process, with the key held only by the call; its answer and the KEGG or Reactome
# tables are saved as tables of the integration store, so the enrichment, its record and its replay need no network.
# The accepted results are kept under the model's entry on the head (references and plain values only).
# ============================================================================

INTEGRATION_ENRICHMENT_RANKING_OPERATION <- "integration_enrichment_ranking"
INTEGRATION_ENRICHMENT_OPERATION <- "integration_enrichment"
INTEGRATION_ENRICHMENT_RANKING_SCHEMA <- "multischolar.integration_enrichment_ranking.v1"
INTEGRATION_ENRICHMENT_SCHEMA <- "multischolar.integration_enrichment.v1"

#' The rule that fixes a factor's direction (INT-18): the factor is multiplied by +1 or -1 so that the weight of largest
#' absolute value over every data set of the model is positive (the first such feature in the model's data set and
#' feature order on a tie).
#' @noRd
INTEGRATION_ENRICHMENT_DIRECTION_RULE <- "largest_absolute_weight_positive"

#' The databases each omic is enriched against (INT-8).
#' @noRd
INTEGRATION_ENRICHMENT_DATABASES <- list(proteomics = "string", transcriptomics = "string",
    metabolomics = c("kegg", "reactome"), lipidomics = c("kegg", "reactome"))

#' KIT-183 (VAL-51, VAL-53): the tests of the shared kernel. `camera`: the test of record, limma's camera on each
#' feature's regression on the oriented factor, the correlation between a set's features estimated; `gsea`: the
#' pre-ranked test on each feature's correlation (KIT-141, unchanged); `camera_pr`: the labelled fallback with an assumed
#' correlation of 0.01, implemented so the simulation can measure it and offered only under outcome O2 of the ticket.
#' KIT-191 (VAL-85 (2)): `camera_floor`, camera with negative estimated correlations counted as zero
#' (`allow.neg.cor = FALSE`), and `permutation`, PCGSE's permutation of the factor's values as KIT-183's simulation
#' measured it (`integrationEnrichmentPermutation()`).
#' @noRd
INTEGRATION_ENRICHMENT_TESTS <- c("camera", "gsea", "camera_pr", "camera_floor", "permutation")

#' KIT-183: the tests offered (outcome O1 of the ticket's simulation, `docs/qa/workflow-execution/kit/kit-183/results.md`);
#' KIT-191 (VAL-85 (2)): the four tests the user chooses among, the default first.
#' @noRd
INTEGRATION_ENRICHMENT_TESTS_OFFERED <- c("camera", "camera_floor", "permutation", "gsea")

#' KIT-183: the tests that run on the data set's values (the correlated declaration) rather than on the ranked list.
#' @noRd
INTEGRATION_ENRICHMENT_CORRELATED_TESTS <- c("camera", "camera_pr", "camera_floor", "permutation")

#' KIT-183 (section 7): the canonical settings text of each correlated test, declared in its identity. KIT-191:
#' `camera_floor` is camera's text with `allow.neg.cor=FALSE`; the permutation's text holds its seed and number of
#' reorderings and is built by `integrationEnrichmentTestSettingsText()`.
#' @noRd
INTEGRATION_ENRICHMENT_TEST_SETTINGS <- c(
    camera = paste0("inter.gene.cor=NA;allow.neg.cor=TRUE;use.ranks=FALSE;trend.var=TRUE;directional=TRUE;",
        "statistic=zscoreT(moderated t, trend, approx=TRUE, method=hill);complete=all_samples;",
        "collapse=largest_mean_among_complete"),
    camera_pr = paste0("inter.gene.cor=0.01;use.ranks=FALSE;directional=TRUE;",
        "statistic=zscoreT(moderated t, trend, approx=TRUE, method=hill);complete=all_samples;",
        "collapse=largest_mean_among_complete"),
    camera_floor = paste0("inter.gene.cor=NA;allow.neg.cor=FALSE;use.ranks=FALSE;trend.var=TRUE;directional=TRUE;",
        "statistic=zscoreT(moderated t, trend, approx=TRUE, method=hill);complete=all_samples;",
        "collapse=largest_mean_among_complete"))

#' KIT-191 (ticket section 6): the permutation test's canonical settings text, with `%d` for the number of reorderings
#' and the seed. The factor's values are reordered with the model unchanged.
#' @noRd
INTEGRATION_ENRICHMENT_PERMUTATION_SETTINGS <- paste0("statistic=fisher_z(pearson);set=pcgse_t(vif=1);",
    "permute=factor_values_free;permutations=%d;p=(1+b)/(1+B);unchanged_draws=counted;two_sided=abs;seed=%d;",
    "rng=Mersenne-Twister/Inversion/Rejection;sums=by_index;complete=all_samples;collapse=largest_mean_among_complete")

#' KIT-191 (D1, VAL-97 (1)): the numbers of reorderings offered for the permutation test, the default first.
#' @noRd
INTEGRATION_ENRICHMENT_PERMUTATIONS <- c(1000L, 10000L)

#' KIT-191 (D4 (c), VAL-97 (4)): the longest predicted run of the permutation test the Enrichment tab accepts, in
#' seconds (the tab waits for the run).
#' @noRd
INTEGRATION_ENRICHMENT_TAB_SECONDS <- 120

#' KIT-191 (5.6): the permutation test's predicted time, a * n * G * B + c * G * B + s * sum(|set|) * B seconds (n
#' biological samples, G features, B reorderings). Each coefficient is twice the seconds per unit of its own term,
#' measured in the case its term dominates (L1, `docs/qa/workflow-execution/kit/kit-191/lifecycle.md`, run below).
#' @noRd
INTEGRATION_ENRICHMENT_PERMUTATION_COST <- list(a = 2e-9, c = 2e-8, s = 2e-8, run = "not measured yet")

#' KIT-191 (VAL-85 (2), ticket 2.3 and 5.4): how each offered test behaved in KIT-183's counted simulation (run
#' 45225854502, package fae77e5f, 800 replicates per cell), per kind of data set (`metabolites`: its metabolomics data
#' set, shown also for lipids, which were not simulated; `genes_proteins`: its proteomics and transcriptomics data
#' sets), number of biological samples and measure. Error rates are shares of null sets with an unadjusted p-value of
#' 0.05 or less (`sprintf("%.1f", 100 * x)`); power is the share of planted sets called at BH <= 0.05 within the call
#' in the right direction (`"%.2f"`, or `"%.3f"` when below 0.1). Measures: A unrelated sets against values unrelated
#' to the data; B sets of the model's most variable unrelated features, same values; C unrelated sets, main factor; D
#' correlated unrelated sets (0.3), main factor; E weaker factors, largest; F at least one false call in a result; C3
#' identifiers carrying 2 or 3 unrelated features; G, H, I power. Regenerated from the committed tables by
#' `docs/qa/workflow-execution/kit/kit-191/evidence-k191.R` (U14).
#' @noRd
INTEGRATION_ENRICHMENT_TEST_EVIDENCE <- local({
    rows <- matrix(c(
        "metabolites", "camera", "12", "A", "3.4 to 4.5 % (largest: 95 % interval 3.7 to 5.3 %)",
        "metabolites", "camera", "12", "B", "4.5 % (largest: 95 % interval 4.1 to 4.9 %)",
        "metabolites", "camera", "12", "C", "0.4 to 4.4 %",
        "metabolites", "camera", "12", "D", "0.6 %",
        "metabolites", "camera", "12", "E", "6.7 % (95 % interval 5.5 to 7.9 %)",
        "metabolites", "camera", "12", "F", "2.2 to 2.5 %",
        "metabolites", "camera", "12", "C3", "4.6 %",
        "metabolites", "camera", "12", "I", "0.020 to 0.022",
        "metabolites", "camera_floor", "12", "A", "1.7 to 3.6 % (largest: 95 % interval 2.7 to 4.5 %)",
        "metabolites", "camera_floor", "12", "B", "2.2 % (largest: 95 % interval 1.9 to 2.5 %)",
        "metabolites", "camera_floor", "12", "C", "0.1 to 2.3 %",
        "metabolites", "camera_floor", "12", "D", "0.5 %",
        "metabolites", "camera_floor", "12", "E", "6.7 % (95 % interval 5.5 to 7.9 %)",
        "metabolites", "camera_floor", "12", "F", "0.1 to 0.4 %",
        "metabolites", "camera_floor", "12", "C3", "2.1 %",
        "metabolites", "camera_floor", "12", "I", "0.001 to 0.002",
        "metabolites", "permutation", "12", "A", "4.2 to 5.2 % (largest: 95 % interval 4.2 to 6.3 %)",
        "metabolites", "permutation", "12", "B", "4.9 % (largest: 95 % interval 4.5 to 5.3 %)",
        "metabolites", "permutation", "12", "C", "0.6 to 4.7 %",
        "metabolites", "permutation", "12", "D", "0.6 %",
        "metabolites", "permutation", "12", "E", "8.6 % (95 % interval 7.2 to 9.9 %)",
        "metabolites", "permutation", "12", "F", "3.2 to 4.2 %",
        "metabolites", "permutation", "12", "C3", "5.0 %",
        "metabolites", "permutation", "12", "I", "0.010",
        "metabolites", "gsea", "12", "A", "3.2 to 3.9 % (largest: 95 % interval 3.7 to 4.2 %); terms whose features vary together: 8.6 to 28.7 %",
        "metabolites", "gsea", "12", "B", "3.2 % (largest: 95 % interval 2.9 to 3.6 %)",
        "metabolites", "gsea", "12", "C", "1.8 to 4.3 %; terms whose features vary together: 4.4 to 23.6 %",
        "metabolites", "gsea", "12", "D", "13.5 %",
        "metabolites", "gsea", "12", "E", "34.4 % (95 % interval 32.1 to 36.6 %)",
        "metabolites", "gsea", "12", "F", "4.2 to 4.4 %",
        "metabolites", "gsea", "12", "C3", "6.5 %",
        "metabolites", "gsea", "12", "I", "0.010 to 0.013",
        "genes_proteins", "camera", "12", "A", "2.2 to 6.6 % (largest: 95 % interval 5.8 to 7.4 %)",
        "genes_proteins", "camera", "12", "B", "9.4 to 10.1 % (largest: 95 % interval 9.5 to 10.7 %)",
        "genes_proteins", "camera", "12", "C", "0.4 to 4.2 %",
        "genes_proteins", "camera", "12", "D", "0.5 to 0.6 %",
        "genes_proteins", "camera", "12", "E", "12.8 % (95 % interval 11.8 to 13.8 %)",
        "genes_proteins", "camera", "12", "F", "1.9 to 2.8 %",
        "genes_proteins", "camera", "12", "G", "0.60 to 0.62",
        "genes_proteins", "camera", "12", "H", "0.99 to 1.00",
        "genes_proteins", "camera_floor", "12", "A", "1.1 to 5.1 % (largest: 95 % interval 4.5 to 5.7 %)",
        "genes_proteins", "camera_floor", "12", "B", "5.4 to 5.5 % (largest: 95 % interval 4.8 to 6.1 %)",
        "genes_proteins", "camera_floor", "12", "C", "0.1 to 1.8 %",
        "genes_proteins", "camera_floor", "12", "D", "0.5 to 0.6 %",
        "genes_proteins", "camera_floor", "12", "E", "12.8 % (95 % interval 11.8 to 13.8 %)",
        "genes_proteins", "camera_floor", "12", "F", "0.0 to 0.1 %",
        "genes_proteins", "camera_floor", "12", "G", "0.47 to 0.50",
        "genes_proteins", "camera_floor", "12", "H", "0.99",
        "genes_proteins", "permutation", "12", "A", "4.8 to 5.3 % (largest: 95 % interval 4.6 to 6.0 %)",
        "genes_proteins", "permutation", "12", "B", "4.9 to 5.3 % (largest: 95 % interval 4.8 to 5.7 %)",
        "genes_proteins", "permutation", "12", "C", "0.4 to 4.4 %",
        "genes_proteins", "permutation", "12", "D", "0.5 to 0.8 %",
        "genes_proteins", "permutation", "12", "E", "12.7 % (95 % interval 11.8 to 13.7 %)",
        "genes_proteins", "permutation", "12", "F", "2.8 to 4.9 %",
        "genes_proteins", "permutation", "12", "G", "0.39 to 0.52",
        "genes_proteins", "permutation", "12", "H", "0.80 to 0.96",
        "genes_proteins", "gsea", "12", "A", "3.9 to 5.5 % (largest: 95 % interval 4.9 to 6.0 %); terms whose features vary together: 28.3 to 62.9 %",
        "genes_proteins", "gsea", "12", "B", "4.1 to 6.1 % (largest: 95 % interval 5.6 to 6.6 %)",
        "genes_proteins", "gsea", "12", "C", "1.8 to 5.1 %; terms whose features vary together: 21.8 to 57.7 %",
        "genes_proteins", "gsea", "12", "D", "40.2 to 43.0 %",
        "genes_proteins", "gsea", "12", "E", "60.1 % (95 % interval 59.0 to 61.2 %)",
        "genes_proteins", "gsea", "12", "F", "2.6 to 7.1 %",
        "genes_proteins", "gsea", "12", "G", "0.63 to 0.72",
        "genes_proteins", "gsea", "12", "H", "1.00",
        "metabolites", "camera", "6", "A", "2.0 to 3.7 % (largest: 95 % interval 3.1 to 4.3 %)",
        "metabolites", "camera", "6", "B", "3.8 % (largest: 95 % interval 3.4 to 4.1 %)",
        "metabolites", "camera", "6", "I", "0.001 to 0.002",
        "metabolites", "camera_floor", "6", "A", "0.4 to 2.2 % (largest: 95 % interval 1.5 to 2.9 %)",
        "metabolites", "camera_floor", "6", "B", "0.3 % (largest: 95 % interval 0.2 to 0.4 %)",
        "metabolites", "camera_floor", "6", "I", "0.000",
        "metabolites", "permutation", "6", "A", "4.4 to 5.0 % (largest: 95 % interval 3.9 to 6.1 %)",
        "metabolites", "permutation", "6", "B", "5.0 % (largest: 95 % interval 4.6 to 5.4 %)",
        "metabolites", "permutation", "6", "I", "0.001 to 0.003",
        "metabolites", "gsea", "6", "A", "4.5 to 4.6 % (largest: 95 % interval 4.3 to 4.8 %); terms whose features vary together: 9.2 to 33.4 %",
        "metabolites", "gsea", "6", "B", "4.2 % (largest: 95 % interval 3.9 to 4.6 %)",
        "metabolites", "gsea", "6", "I", "0.007",
        "genes_proteins", "camera", "6", "A", "2.5 to 5.0 % (largest: 95 % interval 4.4 to 5.6 %)",
        "genes_proteins", "camera", "6", "B", "8.7 to 9.9 % (largest: 95 % interval 9.3 to 10.5 %)",
        "genes_proteins", "camera", "6", "G", "0.05 to 0.13",
        "genes_proteins", "camera", "6", "H", "0.23 to 0.44",
        "genes_proteins", "camera_floor", "6", "A", "0.3 to 4.1 % (largest: 95 % interval 3.6 to 4.7 %)",
        "genes_proteins", "camera_floor", "6", "B", "1.9 to 2.7 % (largest: 95 % interval 2.4 to 3.1 %)",
        "genes_proteins", "camera_floor", "6", "G", "0.000 to 0.001",
        "genes_proteins", "camera_floor", "6", "H", "0.00 to 0.12",
        "genes_proteins", "permutation", "6", "A", "4.8 to 5.3 % (largest: 95 % interval 4.8 to 5.8 %)",
        "genes_proteins", "permutation", "6", "B", "5.0 to 5.3 % (largest: 95 % interval 4.8 to 5.9 %)",
        "genes_proteins", "permutation", "6", "G", "0.003 to 0.008",
        "genes_proteins", "permutation", "6", "H", "0.006 to 0.019",
        "genes_proteins", "gsea", "6", "A", "4.5 to 5.2 % (largest: 95 % interval 4.7 to 5.7 %); terms whose features vary together: 30.3 to 64.7 %",
        "genes_proteins", "gsea", "6", "B", "5.9 to 6.2 % (largest: 95 % interval 5.7 to 6.7 %)",
        "genes_proteins", "gsea", "6", "G", "0.38 to 0.45",
        "genes_proteins", "gsea", "6", "H", "0.95 to 0.99",
        "metabolites", "camera", "20", "A", "4.1 to 5.4 % (largest: 95 % interval 4.3 to 6.5 %)",
        "metabolites", "camera", "20", "B", "5.2 % (largest: 95 % interval 4.8 to 5.7 %)",
        "metabolites", "camera", "20", "I", "0.075 to 0.078",
        "metabolites", "camera_floor", "20", "A", "3.0 to 5.1 % (largest: 95 % interval 4.0 to 6.2 %)",
        "metabolites", "camera_floor", "20", "B", "3.6 % (largest: 95 % interval 3.2 to 4.0 %)",
        "metabolites", "camera_floor", "20", "I", "0.025",
        "metabolites", "permutation", "20", "A", "4.7 to 5.8 % (largest: 95 % interval 4.6 to 6.9 %)",
        "metabolites", "permutation", "20", "B", "5.0 % (largest: 95 % interval 4.5 to 5.4 %)",
        "metabolites", "permutation", "20", "I", "0.019 to 0.023",
        "metabolites", "gsea", "20", "A", "3.6 to 4.0 % (largest: 95 % interval 3.8 to 4.3 %); terms whose features vary together: 8.1 to 27.7 %",
        "metabolites", "gsea", "20", "B", "3.3 % (largest: 95 % interval 2.9 to 3.6 %)",
        "metabolites", "gsea", "20", "I", "0.018 to 0.019",
        "genes_proteins", "camera", "20", "A", "1.9 to 6.4 % (largest: 95 % interval 5.7 to 7.1 %)",
        "genes_proteins", "camera", "20", "B", "8.3 to 10.3 % (largest: 95 % interval 9.3 to 11.4 %)",
        "genes_proteins", "camera", "20", "G", "0.88 to 0.90",
        "genes_proteins", "camera", "20", "H", "1.00",
        "genes_proteins", "camera_floor", "20", "A", "1.4 to 6.4 % (largest: 95 % interval 5.7 to 7.1 %)",
        "genes_proteins", "camera_floor", "20", "B", "5.5 to 7.2 % (largest: 95 % interval 6.3 to 8.1 %)",
        "genes_proteins", "camera_floor", "20", "G", "0.85 to 0.88",
        "genes_proteins", "camera_floor", "20", "H", "1.00",
        "genes_proteins", "permutation", "20", "A", "4.5 to 5.3 % (largest: 95 % interval 4.8 to 5.8 %)",
        "genes_proteins", "permutation", "20", "B", "4.9 to 5.1 % (largest: 95 % interval 4.4 to 5.9 %)",
        "genes_proteins", "permutation", "20", "G", "0.68 to 0.84",
        "genes_proteins", "permutation", "20", "H", "0.99 to 1.00",
        "genes_proteins", "gsea", "20", "A", "3.6 to 5.1 % (largest: 95 % interval 4.6 to 5.6 %); terms whose features vary together: 28.6 to 62.4 %",
        "genes_proteins", "gsea", "20", "B", "4.2 to 5.9 % (largest: 95 % interval 5.5 to 6.4 %)",
        "genes_proteins", "gsea", "20", "G", "0.84 to 0.90",
        "genes_proteins", "gsea", "20", "H", "1.00"
    ), ncol = 5L, byrow = TRUE)
    applies <- c(
        A = paste("C0, C1, C2 (gsea: C0, C2, with correlated sets C1 after the semicolon) against values unrelated to the",
            "data (F_ind), over the data sets of the kind; the 95 % interval of the largest cell"),
        B = "C4 (sets of the model's selected, zero-loading features) against values unrelated to the data (F_ind)",
        C = "C0, C1, C2 (gsea: C0, C2, with C1 after the semicolon) on the main factor (F_main)",
        D = "C1 at a correlation of 0.3 on the main factor (F_main)",
        E_genes_proteins = paste("C4, C4b, C1-out on F_weak and F_none, the largest cell; gsea's figure is C1-out (members",
            "the model did not use)"),
        E_metabolites = paste("C4, C4b on F_weak and F_none, the largest cell; metabolites: C4 and C4b only (no C1-out set",
            "was formed), every test's figure C4b"),
        F = "familywise: at least one false call among a result's C0 and C2 sets, on F_ind and F_main",
        C3 = "C3 (identifiers carrying 2 or 3 unrelated features) against values unrelated to the data (F_ind)",
        power = paste("power at BH <= 0.05 within result tables of about 87 terms (metabolites, 12 % truly linked) or 119",
            "to 152 terms (genes and proteins, 17 to 20 %); in larger tables with fewer linked terms every test finds",
            "fewer, the permutation test most (its floor, KIT-191 5.2 step 8)"))
    measure <- rows[, 4L]
    key <- ifelse(measure %in% c("G", "H", "I"), "power", ifelse(measure == "E", paste0("E_", rows[, 1L]), measure))
    source <- ifelse(measure %in% c("G", "H", "I"), "s/power.tsv", ifelse(measure == "F", "s/familywise.tsv",
        "s/error-rates.tsv"))
    data.frame(kind = rows[, 1L], test = rows[, 2L], samples = as.integer(rows[, 3L]), measure = measure,
        value_text = rows[, 5L], applies_to = unname(applies[key]),
        source = paste("KIT-183 simulation, run 45225854502,", source), stringsAsFactors = FALSE)
})

#' KIT-191: the canonical settings text of a correlated test (`test_settings`, declared in its identity).
#' @param test The test.
#' @param p The resolved settings (`seed`, `permutations`; the permutation test only).
#' @noRd
integrationEnrichmentTestSettingsText <- function(test, p) {
    if (identical(test, "permutation")) {
        return(sprintf(INTEGRATION_ENRICHMENT_PERMUTATION_SETTINGS, as.integer(p$permutations), as.integer(p$seed)))
    }
    INTEGRATION_ENRICHMENT_TEST_SETTINGS[[test]]
}

#' KIT-183: the correlated tests' declaration and product schema.
#' @noRd
INTEGRATION_ENRICHMENT_CORRELATED_OPERATION <- "integration_enrichment_correlated"
INTEGRATION_ENRICHMENT_CORRELATED_SCHEMA <- "multischolar.integration_enrichment_correlated.v1"

#' KIT-183: the smallest number of biological samples at which the test of record's error rate was checked (6.5);
#' below it the tab shows the few-samples sentence.
#' @noRd
INTEGRATION_ENRICHMENT_CHECKED_SAMPLES <- 11L

#' KIT-183: the test of record for a data set of `n_samples` biological samples (outcome O1 of the ticket's
#' simulation: `camera` whatever the number).
#' @param n_samples The data set's number of biological samples (S_d).
#' @return `"camera"`.
#' @noRd
integrationEnrichmentTestOfRecord <- function(n_samples) "camera"

#' The enrichment's settings with their defaults (KIT-141; INT-18, INT-19): the rules chosen in this ticket, for the
#' user to confirm (`settled_by`).
#' @return Named list by setting: `default`, `choices` (NULL for a number or a free value), `settled_by`.
#' @noRd
integrationEnrichmentParameterSpec <- function() {
    chosen <- "KIT-141 (to be confirmed by the user)"
    list(
        # The correlation of each feature with the factor over the biological samples (pairwise complete).
        correlation_method = list(default = "pearson", choices = c("pearson", "spearman"), settled_by = chosen),
        # A feature needs a value in at least this many of the model's biological samples to be ranked.
        min_samples = list(default = 5L, choices = NULL, type = "whole", settled_by = chosen),
        # STRING: the species (NCBI taxon), the false discovery rate of the terms it reports, the rank direction
        # (`ge_enrichment_rank_direction`, as STRING's own example: -1).
        string_species = list(default = 9606L, choices = NULL, type = "whole", settled_by = chosen),
        string_fdr = list(default = 0.05, choices = NULL, type = "number", settled_by = chosen),
        string_rank_direction = list(default = -1L, choices = NULL, type = "whole", settled_by = chosen),
        # KEGG and Reactome: the organism, the pathway sizes tested (measured, mapped members), the seed of the test.
        kegg_organism = list(default = "hsa", choices = NULL, type = "text", settled_by = chosen),
        reactome_organism = list(default = "Homo sapiens", choices = NULL, type = "text", settled_by = chosen),
        min_set_size = list(default = 5L, choices = NULL, type = "whole", settled_by = chosen),
        max_set_size = list(default = 500L, choices = NULL, type = "whole", settled_by = chosen),
        # The seed of the pre-ranked and (KIT-191) permutation tests; the correlation-adjusted tests use none.
        seed = list(default = 20261121L, choices = NULL, type = "whole", settled_by = chosen),
        # KIT-183: KEGG, Reactome (and Gene Ontology): the test; absent, the test of record for the data set
        # (`integrationEnrichmentTestOfRecord()`); STRING ignores it. KIT-191: four tests offered (VAL-85 (2)).
        test = list(default = NULL, choices = INTEGRATION_ENRICHMENT_TESTS_OFFERED, type = "text",
            default_text = paste("the test of record for the data set (camera); or camera_floor, permutation or gsea"),
            settled_by = "VAL-51; VAL-85; KIT-183 simulation, run 45225854502"),
        # KIT-191 (D1): the permutation test's number of random reorderings of the factor's values.
        permutations = list(default = INTEGRATION_ENRICHMENT_PERMUTATIONS[[1L]], choices = INTEGRATION_ENRICHMENT_PERMUTATIONS,
            type = "whole", default_text = "1000 (or 10000; the permutation test only)", settled_by = "VAL-85; VAL-97; KIT-191"))
}

#' Resolve and check an enrichment request against the model (KIT-141). Refused in plain words
#' (`integration_enrichment_parameters_invalid`): an unknown setting (the STRING key is never a setting), a factor or
#' data set the model lacks, a database the data set's omic is not enriched against, a metabolite or lipid data set
#' without its identifier column, a number out of range.
#' @param parameters `factor`, `view`, `database`, for metabolites and lipids `id_column`, and any setting of
#'   `integrationEnrichmentParameterSpec()` (absent ones take their default).
#' @param value The model product.
#' @param object The inputs state with `features` and `inputs`.
#' @param n_samples KIT-183: the data set's number of biological samples, for the test of record
#'   (`integrationEnrichmentTestOfRecord()`); NA when not known.
#' @return The resolved settings, with `omic` and `id_type`; KIT-183: `test` (NA for STRING, which ignores it) and
#'   `test_ignored` (TRUE when a test was given for STRING).
#' @noRd
integrationResolveEnrichmentParameters <- function(parameters, value, object, n_samples = NA_integer_) {
    invalid <- function(detail) integrationEnrichmentRefuse("integration_enrichment_parameters_invalid", detail)
    spec <- integrationEnrichmentParameterSpec()
    known <- c("factor", "view", "database", "id_column", names(spec))
    unknown <- setdiff(names(parameters), known)
    if (length(unknown)) invalid(sprintf("these are not enrichment settings: %s", paste(unknown, collapse = ", ")))
    one <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
    factors <- setdiff(names(value$factors), "biological_sample")
    views <- names(value$weights)
    p <- list(factor = parameters$factor, view = parameters$view, database = parameters$database)
    if (!one(p$factor) || !p$factor %in% factors) invalid(sprintf("the factor must be one of %s", paste(factors, collapse = ", ")))
    if (!one(p$view) || !p$view %in% views) invalid(sprintf("the data set must be one of %s", paste(views, collapse = ", ")))
    p$omic <- integrationViewOmics(object)[[p$view]]
    allowed <- INTEGRATION_ENRICHMENT_DATABASES[[p$omic]]
    if (!one(p$database) || !p$database %in% allowed) {
        invalid(sprintf("the data set %s can be enriched against %s", p$view, paste(allowed, collapse = " or ")))
    }
    for (name in names(spec)) {
        given <- parameters[[name, exact = TRUE]]
        p[name] <- list(if (is.null(given)) spec[[name]]$default else given)
    }
    whole <- function(x, lower, upper) is.numeric(x) && length(x) == 1L && !is.na(x) && x == round(x) && x >= lower &&
        x <= upper
    if (!one(p$correlation_method) || !p$correlation_method %in% spec$correlation_method$choices) {
        invalid("the correlation must be pearson or spearman")
    }
    if (!whole(p$min_samples, 3, 1e6)) invalid("the smallest number of samples per feature must be a whole number of 3 or more")
    if (!whole(p$string_species, 1, 1e9)) invalid("the species must be a STRING species number, such as 9606")
    if (!(is.numeric(p$string_fdr) && length(p$string_fdr) == 1L && !is.na(p$string_fdr) && p$string_fdr > 0 &&
        p$string_fdr <= 1)) invalid("the STRING false discovery rate must be above 0 and at most 1")
    if (!whole(p$string_rank_direction, -1, 1)) invalid("the STRING rank direction must be -1, 0 or 1")
    if (!one(p$kegg_organism) || !grepl("^[a-z]{3,4}$", p$kegg_organism)) {
        invalid("the KEGG organism must be a KEGG organism code, such as hsa")
    }
    if (!one(p$reactome_organism)) invalid("the Reactome organism must be named, such as Homo sapiens")
    if (!whole(p$min_set_size, 5, 1e6) || !whole(p$max_set_size, p$min_set_size, 1e6)) {
        invalid("the pathway sizes must be whole numbers, the smallest at least 5 and not above the largest")
    }
    if (!whole(p$seed, 0, .Machine$integer.max)) invalid("the seed must be a whole number")
    if (!whole(p$permutations, 1, 1e6) || !p$permutations %in% spec$permutations$choices) {
        invalid("the number of reorderings of the permutation test must be 1000 or 10000")
    }
    p$permutations <- as.integer(p$permutations)
    p$min_samples <- as.integer(p$min_samples)
    p$string_species <- as.integer(p$string_species)
    p$string_rank_direction <- as.integer(p$string_rank_direction)
    p$min_set_size <- as.integer(p$min_set_size)
    p$max_set_size <- as.integer(p$max_set_size)
    p$seed <- as.integer(p$seed)
    p$id_type <- switch(p$database, string = INTEGRATION_ENRICHMENT_ID_TYPES[[p$omic]][[1L]], kegg = "kegg",
        reactome = "chebi")
    p$id_column <- NA_character_
    if (p$database %in% c("kegg", "reactome")) {
        column <- parameters$id_column
        features <- object@features[[p$view]]
        if (!one(column) || !column %in% setdiff(names(features), "feature_id")) {
            invalid(sprintf("choose the column of the %s feature table that holds the %s identifiers", p$view,
                if (p$database == "kegg") "KEGG compound" else "ChEBI"))
        }
        p$id_column <- column
    }
    # KIT-183: STRING has one test, its own, and ignores the setting (never a refusal: the workbook and the replay
    # keep working whatever they send); KEGG and Reactome take the test given, or the test of record.
    if (identical(p$database, "string")) {
        p$test_ignored <- !is.null(p$test)
        p$test <- NA_character_
    } else {
        p$test_ignored <- FALSE
        if (is.null(p$test)) {
            p$test <- integrationEnrichmentTestOfRecord(n_samples)
        } else if (!one(p$test) || !p$test %in% INTEGRATION_ENRICHMENT_TESTS) {
            invalid(paste("the test must be camera (the correlation-adjusted test), camera_floor (the same, negative",
                "correlations counted as zero), permutation (the permutation test) or gsea (the pre-ranked test)"))
        } else if (!p$test %in% INTEGRATION_ENRICHMENT_TESTS_OFFERED) {
            integrationEnrichmentRefuse("integration_enrichment_test_not_offered",
                values = list("the correlation-adjusted test with an assumed correlation is not offered in this release"))
        }
    }
    p
}

#' The direction of a factor (INT-18; `INTEGRATION_ENRICHMENT_DIRECTION_RULE`): +1 or -1 so that the weight of
#' largest absolute value over every data set of the model is positive. The weights are read one data set at a time.
#' @param manager,value The store's owner and the model product.
#' @param factor The factor.
#' @return `list(rule, sign, view, feature_id, weight)`.
#' @noRd
integrationEnrichmentFactorDirection <- function(manager, value, factor) {
    best <- list(weight = 0, view = NA_character_, feature_id = NA_character_)
    for (view in names(value$weights)) {
        w <- integrationModelViewTable(manager, value, "weights", view)
        x <- w[[factor]]
        if (!length(x) || all(is.na(x))) next
        i <- which.max(abs(x))
        if (abs(x[[i]]) > abs(best$weight)) best <- list(weight = x[[i]], view = view, feature_id = as.character(w$feature_id[[i]]))
    }
    list(rule = INTEGRATION_ENRICHMENT_DIRECTION_RULE, sign = if (isTRUE(best$weight < 0)) -1L else 1L,
        view = best$view, feature_id = best$feature_id, weight = best$weight)
}

#' Each feature's correlation with the factor (INT-18), pairwise complete: over the samples where both are observed.
#' A feature with fewer than `min_samples` such samples, or without spread there, gets no correlation.
#' @param x Features x samples (the samples in the order of `f`).
#' @param f The factor's values, oriented.
#' @param method `pearson` or `spearman`.
#' @param min_samples The smallest number of samples.
#' @return `data.frame(n_samples, correlation, status)`; `status` `ranked`, `too_few_samples` or `no_variation`.
#' @noRd
integrationEnrichmentCorrelations <- function(x, f, method, min_samples) {
    out <- lapply(seq_len(nrow(x)), function(i) {
        ok <- !is.na(x[i, ]) & !is.na(f)
        n <- sum(ok)
        if (n < min_samples) return(list(n, NA_real_, "too_few_samples"))
        a <- x[i, ok]
        if (!(stats::var(a) > 0) || !(stats::var(f[ok]) > 0)) return(list(n, NA_real_, "no_variation"))
        list(n, stats::cor(a, f[ok], method = method), "ranked")
    })
    data.frame(n_samples = vapply(out, \(r) as.integer(r[[1L]]), integer(1)),
        correlation = vapply(out, \(r) as.numeric(r[[2L]]), numeric(1)),
        status = vapply(out, \(r) r[[3L]], character(1)), stringsAsFactors = FALSE)
}

#' Rank every measured feature of a data set by its correlation with a factor (INT-18): the data set is read from the
#' inputs state in row pages (`integrationStateViewPages`, at most 10,000 features at a time; never whole); the runs of
#' each biological sample averaged as the model does (`integrationAverageRuns`, the model's sample match applied); the
#' model's biological samples only; internal standards left out when the model left them out; the factor oriented by
#' `integrationEnrichmentFactorDirection()`. Every feature is kept in the table with its status.
#' @param object The inputs state (`features`, `design_matrix`, `inputs`; the data sets' column names).
#' @param manager,source The integration owner and the inputs state's reference.
#' @param value The model product; `model_parameters` its resolved parameters (sample match, internal standards).
#' @param factor,view,method,min_samples The resolved settings.
#' @param pages `function(view, columns, kernel)` (`integrationStateViewPages()`'s form; tests may pass another).
#' @return The ranking product: `schema`, `factor`, `view`, `omic`, `method`, `min_samples`, `direction`, `samples`,
#'   `ranking` (`feature_id`, `n_samples`, `correlation`, `status`), `counts`.
#' @noRd
integrationEnrichmentRank <- function(object, manager, source, value, model_parameters, factor, view, method, min_samples,
    pages = integrationStateViewPages(manager, source)) {
    design <- integrationApplySampleMatch(object@design_matrix, model_parameters$sample_match_table)
    omic <- integrationViewOmics(object)[[view]]
    columns <- names(object@views[[view]])
    rows <- design[design$omic == omic, , drop = FALSE]
    key <- rows$biological_sample[match(columns[-1L], rows$run)]
    if (anyNA(key)) integrationRefuse(sprintf("view %s has a run that its omic's design does not list", view), "integration_input_invalid")
    direction <- integrationEnrichmentFactorDirection(manager, value, factor)
    samples <- as.character(value$factors$biological_sample)
    f <- direction$sign * as.numeric(value$factors[[factor]])
    standards <- if (identical(model_parameters$internal_standards, "drop")) integrationInternalStandards(object@features[[view]])
    parts <- list()
    pages(view, columns, function(page, offset) {
        if (!nrow(page)) return(invisible(NULL))
        m <- as.matrix(page[-1L])
        storage.mode(m) <- "double"
        rownames(m) <- as.character(page[[1L]])
        averaged <- integrationAverageRuns(m, key)
        x <- matrix(NA_real_, nrow(m), length(samples))
        here <- intersect(samples, colnames(averaged))
        x[, match(here, samples)] <- averaged[, here, drop = FALSE]
        part <- data.frame(feature_id = rownames(m), integrationEnrichmentCorrelations(x, f, method, min_samples),
            stringsAsFactors = FALSE)
        if (!is.null(standards)) {
            flagged <- standards[offset + seq_len(nrow(m))]
            part$n_samples[flagged] <- NA_integer_
            part$correlation[flagged] <- NA_real_
            part$status[flagged] <- "internal_standard"
        }
        parts[[length(parts) + 1L]] <<- part
        invisible(NULL)
    })
    ranking <- do.call(rbind, parts)
    rownames(ranking) <- NULL
    counts <- list(features = nrow(ranking), ranked = sum(ranking$status == "ranked"),
        too_few_samples = sum(ranking$status == "too_few_samples"), no_variation = sum(ranking$status == "no_variation"),
        internal_standard = sum(ranking$status == "internal_standard"))
    list(schema = INTEGRATION_ENRICHMENT_RANKING_SCHEMA, factor = factor, view = view, omic = omic, method = method,
        min_samples = min_samples, direction = direction, samples = length(samples), ranking = ranking, counts = counts)
}

#' The ranking declaration (KIT-141): product only, on the inputs state; its identity is the state's, the model's entity,
#' the factor, data set, correlation, smallest number of samples and the direction rule. The model's value and
#' parameters are handed to the kernel beside the declared parameters (`model_value`, `model_parameters`).
#' @noRd
integrationEnrichmentRankingOperation <- function() {
    newWorkflowOperation(INTEGRATION_ENRICHMENT_RANKING_OPERATION, "IntegrationViewsData",
        fields = c("features", "design_matrix", "inputs"), unit = "whole_matrix", output_class = NULL,
        parameters = c("model", "factor", "view", "correlation_method", "min_samples", "direction_rule"),
        product_type = "integration.enrichment_ranking", factors = list(product_schema = INTEGRATION_ENRICHMENT_RANKING_SCHEMA),
        kernel = function(object, parameters, store) integrationEnrichmentRank(object, store$manager, store$source,
            parameters$model_value, parameters$model_parameters, parameters$factor, parameters$view,
            parameters$correlation_method, parameters$min_samples))
}

#' The ranked list of a ranking, mapped to the database's identifiers (INT-18, INT-19; `integrationEnrichmentMapIdentifiers`):
#' the ranked features only (every feature with a correlation), each identifier once (several features on one
#' identifier: the one with the largest absolute correlation kept, the first in the data set's order on a tie).
#' @param ranking The ranking product.
#' @param p The resolved settings.
#' @param features The data set's feature table (the identifier column for metabolites and lipids).
#' @noRd
integrationEnrichmentMappedList <- function(ranking, p, features) {
    r <- ranking$ranking[ranking$ranking$status == "ranked", , drop = FALSE]
    ranked <- data.frame(feature_id = r$feature_id, value = r$correlation, stringsAsFactors = FALSE)
    if (p$database %in% c("kegg", "reactome")) {
        integrationEnrichmentMapIdentifiers(ranked, p$omic, p$id_type, as.data.frame(features), p$id_column)
    } else {
        integrationEnrichmentMapIdentifiers(ranked, p$omic, p$id_type)
    }
}

# --- KIT-183: the correlation-adjusted test of record (VAL-51, VAL-53) ---

#' The biological samples of a data set (KIT-183 5.1 step 1, S_d): the model's biological samples that have at least
#' one run in the data set, after the model's sample match, in the model's order.
#' @param object The inputs state (`design_matrix`, `inputs`; the data sets' column names).
#' @param value The model product; `model_parameters` its resolved parameters (the sample match).
#' @param view The data set.
#' @return Character vector of biological samples.
#' @noRd
integrationEnrichmentViewSamples <- function(object, value, model_parameters, view) {
    design <- integrationApplySampleMatch(object@design_matrix, model_parameters$sample_match_table)
    rows <- design[design$omic == integrationViewOmics(object)[[view]], , drop = FALSE]
    key <- rows$biological_sample[match(names(object@views[[view]])[-1L], rows$run)]
    samples <- as.character(value$factors$biological_sample)
    samples[samples %in% key]
}

#' The design of the correlated tests (KIT-183 5.1 step 5, 8.8): the columns `intercept` and `factor` over S_d, then
#' one column per covariate, named by the caller (KIT-184 adds RUV's removed factors through this argument). With
#' `covariates = NULL` the design is exactly [1, f]. Every call takes the factor's column by its name, so added columns
#' never move the tested coefficient.
#' @param f The oriented factor scores over S_d.
#' @param covariates NULL, or a numeric matrix or data frame with one row per sample of S_d in S_d's order (aligned by
#'   `integrationEnrichmentTestValues()`), named columns other than `intercept` and `factor`, no missing value.
#' @return The design matrix. A design whose rank is below its number of columns is refused
#'   (`integration_enrichment_parameters_invalid`).
#' @noRd
integrationEnrichmentTestDesign <- function(f, covariates = NULL) {
    invalid <- function() integrationEnrichmentRefuse("integration_enrichment_parameters_invalid",
        "the covariates of the test's design are not valid")
    f <- as.numeric(f)
    design <- cbind(intercept = 1, factor = f)
    if (!is.null(covariates)) {
        w <- as.matrix(covariates)
        names <- colnames(w)
        if (!is.numeric(w) || nrow(w) != length(f) || anyNA(w) || is.null(names) || any(!nzchar(names)) ||
            anyDuplicated(names) || any(names %in% c("intercept", "factor"))) {
            invalid()
        }
        design <- cbind(design, w)
    }
    if (anyNA(design)) {
        integrationEnrichmentRefuse("integration_enrichment_parameters_invalid", "the factor has a missing value")
    }
    if (qr(design)$rank < ncol(design)) {
        if (is.null(covariates)) {
            integrationEnrichmentRefuse("integration_enrichment_parameters_invalid",
                "the factor does not vary over the biological samples")
        }
        invalid()
    }
    design
}

#' The values the correlated tests read (KIT-183 5.1 steps 1 to 5): the data set's values in the inputs state's slot
#' `views` (the normalised values before RUV, the slot the ranking reads), read in row pages as the ranking reads them
#' (only the rows of the ranked, mapped features kept), the runs of each biological sample averaged
#' (`integrationAverageRuns`), on S_d; the features with a value in every sample of S_d (no imputation); among the
#' complete features on one identifier, the one with the largest mean value over S_d (the first in the data set's order
#' on a tie; blind to the test's statistic, VAL-53 (3)); the design [1, f] (`integrationEnrichmentTestDesign`).
#' @param object The inputs state (`features`, `design_matrix`, `inputs`).
#' @param manager,source The integration owner and the inputs state's reference.
#' @param value,model_parameters The model product and its resolved parameters.
#' @param ranking The ranking product (factor, direction, smallest number of samples).
#' @param mapped `integrationEnrichmentMappedList()`'s value (the ranked features and their identifiers).
#' @param view The data set.
#' @param pages `function(view, columns, kernel)` (`integrationStateViewPages()`'s form; tests may pass another).
#' @param covariates NULL in this ticket; KIT-184: a numeric matrix or data frame whose row names are biological
#'   samples (a superset of S_d allowed), aligned here to S_d; a sample of S_d without a row is refused.
#' @return `list(y, f, design, samples, feature_id, counts, design_columns, design_covariates)`: `y` the identifiers of
#'   U_cam (rows, named by identifier, in the data set's order) by S_d; `feature_id` the feature kept per identifier;
#'   `counts` `samples`, `mapped`, `complete`, `incomplete`, `collapsed`, `universe`. Refused below the sample floor
#'   (`integration_enrichment_too_few_samples_for_test`) or with fewer than 10 identifiers
#'   (`integration_enrichment_too_few_complete`).
#' @noRd
integrationEnrichmentTestValues <- function(object, manager, source, value, model_parameters, ranking, mapped, view,
    pages = integrationStateViewPages(manager, source), covariates = NULL) {
    samples <- integrationEnrichmentViewSamples(object, value, model_parameters, view)
    n <- length(samples)
    needed <- max(4L, as.integer(ranking$min_samples), 2L + (if (is.null(covariates)) 0L else NCOL(covariates)) + 2L)
    if (n < needed) {
        integrationEnrichmentRefuse("integration_enrichment_too_few_samples_for_test",
            sprintf("%d biological samples in %s", n, view), values = list(needed), samples = n, needed = needed)
    }
    f <- ranking$direction$sign * as.numeric(value$factors[[ranking$factor]])[match(samples,
        as.character(value$factors$biological_sample))]
    if (!is.null(covariates)) {
        at <- match(samples, rownames(covariates))
        if (is.null(rownames(covariates)) || anyNA(at)) {
            integrationEnrichmentRefuse("integration_enrichment_parameters_invalid",
                "the covariates of the test's design are not valid")
        }
        covariates <- covariates[at, , drop = FALSE]
    }
    design <- integrationEnrichmentTestDesign(f, covariates)
    m <- mapped$mapping
    candidate <- m$status %in% c("mapped", "duplicated")
    ids <- stats::setNames(as.character(m$identifier[candidate]), as.character(m$feature_id[candidate]))
    runs <- integrationApplySampleMatch(object@design_matrix, model_parameters$sample_match_table)
    runs <- runs[runs$omic == integrationViewOmics(object)[[view]], , drop = FALSE]
    columns <- names(object@views[[view]])
    key <- runs$biological_sample[match(columns[-1L], runs$run)]
    parts <- list()
    pages(view, columns, function(page, offset) {
        keep <- as.character(page[[1L]]) %in% names(ids)
        if (!any(keep)) return(invisible(NULL))
        x <- as.matrix(page[keep, -1L, drop = FALSE])
        storage.mode(x) <- "double"
        rownames(x) <- as.character(page[[1L]])[keep]
        parts[[length(parts) + 1L]] <<- integrationAverageRuns(x, key)[, samples, drop = FALSE]
        invisible(NULL)
    })
    x <- if (length(parts)) do.call(rbind, parts) else matrix(numeric(), 0L, n, dimnames = list(NULL, samples))
    complete <- stats::complete.cases(x)
    x <- x[complete, , drop = FALSE]
    identifier <- unname(ids[rownames(x)])
    means <- rowMeans(x)
    by_identifier <- order(identifier, -means, seq_along(means))
    kept <- sort(by_identifier[!duplicated(identifier[by_identifier])])
    y <- x[kept, , drop = FALSE]
    feature_id <- rownames(y)
    rownames(y) <- identifier[kept]
    counts <- list(samples = n, mapped = length(ids), complete = nrow(x), incomplete = length(ids) - nrow(x),
        collapsed = nrow(x) - length(kept), universe = length(kept))
    if (counts$universe < INTEGRATION_ENRICHMENT_MIN_MAPPED) {
        integrationEnrichmentRefuse("integration_enrichment_too_few_complete",
            sprintf("%d of %d matched identifiers measured in all %d biological samples", counts$universe,
                length(unique(ids)), n), counts = counts)
    }
    # KIT-184 states the canonical text of a design with covariates.
    list(y = y, f = f, design = design, samples = samples, feature_id = feature_id, counts = counts,
        design_columns = colnames(design), design_covariates = if (is.null(covariates)) "none" else NA_character_)
}

#' The correlated tests on one database's sets (KIT-183 5.1 steps 6 and 7; 5.3). `camera`: limma's camera on the
#' values with the design and the factor's column by name, the correlation within each set estimated
#' (`inter.gene.cor = NA`), negative estimates kept (`allow.neg.cor = TRUE`, VAL-53 (2)), moderated variances with an
#' abundance trend, parametric, directional. `camera_pr`: limma's cameraPR on camera's own statistic (the moderated t
#' of the factor's coefficient, `eBayes(lmFit(y, design), trend = TRUE)`, as z-scores by `zscoreT(approx = TRUE,
#' method = "hill")`) with the correlation assumed to be 0.01. An error of limma is refused
#' (`integration_enrichment_method_failed`); no other test is run in its place (INT-19).
#' @param values `integrationEnrichmentTestValues()`'s value.
#' @param term2gene,term2name The sets (restricted to U_cam) and their names.
#' @param settings `min_size`, `max_size` (members of U_cam) and, for the simulation's sensitivity arm of `camera` only,
#'   `allow_neg_cor` (absent: TRUE).
#' @param test `camera`, `camera_pr` or (KIT-191 5.3) `camera_floor`: camera with `allow.neg.cor = FALSE`, the
#'   correlation column the estimate before the floor.
#' @return `list(method, engine, settings, terms, test_counts, design_columns, design_covariates)`; `terms` one row per
#'   tested set: `ID`, `Description`, `setSize`, `Direction`, `correlation`, `effect`, `pvalue`, `p.adjust` (BH over
#'   the sets of the call), `test`.
#' @noRd
integrationEnrichmentCamera <- function(values, term2gene, term2name, settings, test) {
    y <- values$y
    design <- values$design
    coefficient <- which(colnames(design) == "factor")
    # KIT-191 (5.3): the test sets the floor; the internal hook keeps its meaning for `camera` only.
    allow_neg_cor <- if (identical(test, "camera_floor")) FALSE else settings$allow_neg_cor %||% TRUE
    index <- split(as.character(term2gene$gene), as.character(term2gene$term))
    sizes <- lengths(index)
    index <- index[sizes >= settings$min_size & sizes <= settings$max_size]
    engine <- list(limma = as.character(utils::packageVersion("limma")))
    out <- list(method = test, engine = engine,
        settings = list(min_size = settings$min_size, max_size = settings$max_size, allow_neg_cor = allow_neg_cor,
            test_settings = INTEGRATION_ENRICHMENT_TEST_SETTINGS[[test]]),
        test_counts = values$counts, design_columns = colnames(design),
        design_covariates = values$design_covariates %||% "none")
    if (!length(index)) {
        out$terms <- data.frame(ID = character(), Description = character(), setSize = integer(), Direction = character(),
            correlation = numeric(), effect = numeric(), pvalue = numeric(), p.adjust = numeric(), test = character(),
            stringsAsFactors = FALSE)
        return(out)
    }
    outcome <- tryCatch({
        fit <- limma::eBayes(limma::lmFit(y, design), trend = TRUE)
        z <- limma::zscoreT(fit$t[, coefficient], df = fit$df.total, approx = TRUE, method = "hill")
        names(z) <- rownames(y)
        table <- if (test %in% c("camera", "camera_floor")) {
            limma::camera(y, index, design = design, contrast = coefficient, inter.gene.cor = NA,
                allow.neg.cor = allow_neg_cor, use.ranks = FALSE, trend.var = TRUE, sort = FALSE, directional = TRUE)
        } else {
            limma::cameraPR(z, index, use.ranks = FALSE, inter.gene.cor = 0.01, sort = FALSE, directional = TRUE)
        }
        list(z = z, table = table)
    }, error = \(e) e)
    if (inherits(outcome, "error")) {
        integrationEnrichmentRefuse("integration_enrichment_method_failed", conditionMessage(outcome),
            engine_message = conditionMessage(outcome))
    }
    z <- outcome$z
    table <- outcome$table
    # camera's own contrast: the mean of z over the set minus the mean over the other features of U_cam.
    effect <- vapply(index, function(members) {
        inside <- names(z) %in% members
        mean(z[inside]) - mean(z[!inside])
    }, numeric(1))
    names_by_term <- stats::setNames(as.character(term2name$name), as.character(term2name$term))
    description <- unname(names_by_term[names(index)])
    description[is.na(description)] <- names(index)[is.na(description)]
    terms <- data.frame(ID = names(index), Description = description, setSize = as.integer(table$NGenes),
        Direction = as.character(table$Direction),
        correlation = if (test %in% c("camera", "camera_floor")) as.numeric(table$Correlation) else 0.01,
        effect = unname(effect), pvalue = as.numeric(table$PValue), stringsAsFactors = FALSE)
    terms$p.adjust <- stats::p.adjust(terms$pvalue, "BH")
    terms$test <- test
    terms <- terms[order(terms$pvalue, terms$ID, method = "radix"), , drop = FALSE]
    rownames(terms) <- NULL
    out$terms <- terms
    out
}

#' KIT-191 (5.2 step 2): Fisher's z of correlations, clamped to +-(1 - 1e-12) as KIT-183's simulation did.
#' @noRd
integrationEnrichmentFisherZ <- function(r) atanh(pmax(pmin(r, 1 - 1e-12), -1 + 1e-12))

#' KIT-191 (5.2 steps 3 and 4): PCGSE's two-sample t of each set against the other features (variance inflation 1,
#' pooled over in and out, G - 2 degrees of freedom), from the sums and sums of squares in and out of each set, the
#' set's sums taken by indexing its members' rows (the arithmetic of KIT-183's `pcgseT()`, without its dense
#' set-by-feature matrix).
#' @param stat Features x columns (one column per reordering, or one).
#' @param members List of row indices, one element per set.
#' @return `list(t, effect)`: sets x columns; `effect` the mean in the set minus the mean outside.
#' @noRd
integrationEnrichmentSetT <- function(stat, members) {
    stat <- as.matrix(stat)
    g <- nrow(stat)
    squares <- stat^2
    sum_all <- colSums(stat)
    sq_all <- colSums(squares)
    t <- effect <- matrix(NA_real_, length(members), ncol(stat))
    for (i in seq_along(members)) {
        rows <- members[[i]]
        m1 <- length(rows)
        m2 <- g - m1
        sum_in <- colSums(stat[rows, , drop = FALSE])
        sq_in <- colSums(squares[rows, , drop = FALSE])
        mean_a <- sum_in / m1
        mean_b <- (sum_all - sum_in) / m2
        ss_a <- sq_in - m1 * mean_a^2
        ss_b <- (sq_all - sq_in) - m2 * mean_b^2
        pooled <- sqrt((ss_a + ss_b) / (m1 + m2 - 2))
        effect[i, ] <- mean_a - mean_b
        t[i, ] <- effect[i, ] / (pooled * sqrt(1 / m1 + 1 / m2))
    }
    list(t = t, effect = effect)
}

#' KIT-191 (5.6): the permutation test's predicted time in seconds (`INTEGRATION_ENRICHMENT_PERMUTATION_COST`).
#' @param n_samples,features,permutations n_d, G and B.
#' @param members The summed sizes of the sets tested.
#' @param cost The coefficients.
#' @noRd
integrationEnrichmentPermutationSeconds <- function(n_samples, features, permutations, members,
    cost = INTEGRATION_ENRICHMENT_PERMUTATION_COST) {
    cost$a * n_samples * features * permutations + cost$c * features * permutations + cost$s * members * permutations
}

#' KIT-191 (VAL-85 (2), VAL-97; ticket 5.2): the permutation test of PCGSE (Frost et al. 2015) as KIT-183's simulation
#' measured it (`docs/qa/workflow-execution/kit/kit-183/simulate-k183.R`, `comparators()`). Each feature's Pearson
#' correlation with the oriented factor over S_d, Fisher-transformed; each set's uncorrected two-sample t against the
#' other features of U_cam; B free reorderings of the factor's values across the biological samples, drawn once per
#' call with the declared seed and R's default generator kinds (the session's random state is restored), the model
#' unchanged; the same reorderings for every set. p = (1 + b) / (1 + B), b the reorderings whose |t*| is at least |t|,
#' a reordering that leaves the factor's values in their own order (the identity, or one among tied values) counted
#' without comparing (Phipson and Smyth 2010). Two-sided; `Direction` the sign of t; BH over the sets of the call.
#' Reorderings are processed in blocks of min(B, 1000, floor(5e6 / G)) columns. Refusals as `camera`; a run predicted
#' above `settings$time_limit` seconds is refused (`integration_enrichment_test_too_large`; the tab), without a limit the
#' predicted time is given in a message first (scripts). No other test is run in its place (INT-19).
#' @param values `integrationEnrichmentTestValues()`'s value (`y`, `f`, `counts`; covariates are never used).
#' @param term2gene,term2name The sets (restricted to U_cam) and their names.
#' @param settings `min_size`, `max_size`, `permutations` (absent: 1000), `time_limit` (seconds, or absent).
#' @param seed The declared seed.
#' @return `camera`'s form with `correlation` NA and `test = "permutation"`, plus `permutation`: `seed`,
#'   `permutations`, `identity_draws` (d, the reorderings that left the factor unchanged), `min_p` ((1 + d)/(1 + B)),
#'   `samples`.
#' @noRd
integrationEnrichmentPermutation <- function(values, term2gene, term2name, settings, seed) {
    invalid <- function(detail) integrationEnrichmentRefuse("integration_enrichment_parameters_invalid", detail)
    whole <- function(x, lower) is.numeric(x) && length(x) == 1L && !is.na(x) && x == round(x) && x >= lower
    permutations <- settings$permutations %||% INTEGRATION_ENRICHMENT_PERMUTATIONS[[1L]]
    if (!whole(permutations, 1)) invalid("the number of reorderings of the permutation test must be a whole number of 1 or more")
    if (!whole(seed, 0) || seed > .Machine$integer.max) invalid("the seed must be a whole number")
    permutations <- as.integer(permutations)
    seed <- as.integer(seed)
    y <- values$y
    f <- as.numeric(values$f)
    n <- ncol(y)
    g <- nrow(y)
    if (length(f) != n || anyNA(f) || !(stats::var(f) > 0)) invalid("the factor does not vary over the biological samples")
    index <- split(as.character(term2gene$gene), as.character(term2gene$term))
    sizes <- lengths(index)
    index <- index[sizes >= settings$min_size & sizes <= settings$max_size]
    test_settings <- integrationEnrichmentTestSettingsText("permutation", list(seed = seed, permutations = permutations))
    out <- list(method = "permutation", engine = list(R = as.character(getRversion())),
        settings = list(min_size = settings$min_size, max_size = settings$max_size, permutations = permutations, seed = seed,
            test_settings = test_settings),
        test_counts = values$counts, design_covariates = "none")
    if (!length(index)) {
        out$terms <- data.frame(ID = character(), Description = character(), setSize = integer(), Direction = character(),
            correlation = numeric(), effect = numeric(), pvalue = numeric(), p.adjust = numeric(), test = character(),
            stringsAsFactors = FALSE)
        out$permutation <- list(seed = seed, permutations = permutations, identity_draws = NA_integer_, min_p = NA_real_,
            samples = n)
        return(out)
    }
    # The guard (D4 (c)): refused in the tab above its limit; in scripts the predicted time is given first.
    predicted <- integrationEnrichmentPermutationSeconds(n, g, permutations, sum(lengths(index)))
    limit <- settings$time_limit
    if (!is.null(limit) && predicted > limit) {
        integrationEnrichmentRefuse("integration_enrichment_test_too_large",
            values = list(sprintf("%.1f", predicted / 60), g, n, length(index), permutations), predicted_seconds = predicted)
    }
    if (is.null(limit)) {
        rlang::inform(sprintf(integrationGuiTexts()[["enrich_test_permutation_time"]], sprintf("%.0f", predicted), g, n,
            length(index), permutations), class = "integration_enrichment_permutation_time")
    }
    members <- lapply(index, function(set) sort(match(set, rownames(y))))
    r <- as.numeric(stats::cor(t(y), f))
    observed <- integrationEnrichmentSetT(integrationEnrichmentFisherZ(r), members)
    t_obs <- observed$t[, 1L]
    perms <- withr::with_seed(seed, vapply(seq_len(permutations), function(i) sample(n), integer(n)),
        .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
    perms <- matrix(perms, n, permutations)
    # Reorderings that leave the factor's values in their own order are counted, never compared (5.2 step 5).
    unchanged <- vapply(seq_len(permutations), function(b) identical(f[perms[, b]], f), logical(1))
    ys <- t(scale(t(y)))
    width <- as.integer(min(permutations, 1000L, max(1L, floor(5e6 / g))))
    exceed <- numeric(length(members))
    for (start in seq.int(1L, permutations, by = width)) {
        columns <- start:min(permutations, start + width - 1L)
        fs <- scale(matrix(f[perms[, columns, drop = FALSE]], n))
        zp <- integrationEnrichmentFisherZ((ys %*% fs) / (n - 1))
        at_least <- abs(integrationEnrichmentSetT(zp, members)$t) >= abs(t_obs)
        at_least[, unchanged[columns]] <- TRUE
        exceed <- exceed + rowSums(at_least)
    }
    identity_draws <- sum(unchanged)
    names_by_term <- stats::setNames(as.character(term2name$name), as.character(term2name$term))
    description <- unname(names_by_term[names(index)])
    description[is.na(description)] <- names(index)[is.na(description)]
    terms <- data.frame(ID = names(index), Description = description, setSize = as.integer(lengths(index)),
        Direction = ifelse(t_obs > 0, "Up", "Down"), correlation = NA_real_, effect = observed$effect[, 1L],
        pvalue = (1 + exceed) / (1 + permutations), stringsAsFactors = FALSE)
    terms$p.adjust <- stats::p.adjust(terms$pvalue, "BH")
    terms$test <- "permutation"
    terms <- terms[order(terms$pvalue, terms$ID, method = "radix"), , drop = FALSE]
    rownames(terms) <- NULL
    out$terms <- terms
    out$permutation <- list(seed = seed, permutations = permutations, identity_draws = as.integer(identity_draws),
        min_p = (1 + identity_draws) / (1 + permutations), samples = n)
    out
}

#' The STRING options of the resolved settings.
#' @noRd
integrationEnrichmentStringOptions <- function(p, host = INTEGRATION_STRING_DEFAULT_HOST) {
    integrationStringOptions(ge_fdr = p$string_fdr, rank_direction = p$string_rank_direction, host = host)
}

#' Save a STRING answer in the integration store (KIT-141; audit ENG-07): one table of text, a row per part (the
#' service, the method, the release STRING reported, the request's summary, its submission and status replies as JSON,
#' the results table as STRING sent it). The key is in none of them (`integrationStringRankedEnrichment` removes it).
#' @param context The owner's context, in its state scope.
#' @param response `integrationStringRankedEnrichment()`'s value.
#' @return `list(reference, sha256)`.
#' @noRd
integrationStoreServiceResponse <- function(context, response) {
    json <- function(x) as.character(jsonlite::toJSON(x, auto_unbox = TRUE, digits = NA, null = "null", na = "null"))
    table <- data.frame(part = c("service", "method", "version", "request", "submission", "status", "results_tsv"),
        text = c(response$service, response$method, json(response$version), json(response$request),
            json(response$response$submission), json(response$response$status), response$response$results_tsv),
        stringsAsFactors = FALSE)
    list(reference = saveWorkflowOwnedTable(context, table), sha256 = integrationEnrichmentTableSha256(table))
}

#' The SHA-256 of a saved table's content (its columns, values as text), the identity the enrichment declares.
#' @noRd
integrationEnrichmentTableSha256 <- function(table) {
    digest::digest(lapply(as.list(table), as.character), algo = "sha256")
}

#' A STRING answer read back from its saved table, in `integrationStringRankedEnrichment()`'s form.
#' @param table The saved table (`part`, `text`).
#' @noRd
integrationServiceResponseFromTable <- function(table) {
    text <- stats::setNames(as.character(table$text), as.character(table$part))
    parse <- function(x) jsonlite::fromJSON(x, simplifyVector = FALSE)
    list(service = text[["service"]], method = text[["method"]], version = parse(text[["version"]]),
        request = parse(text[["request"]]), response = list(submission = parse(text[["submission"]]),
            status = parse(text[["status"]]), results_tsv = text[["results_tsv"]]),
        terms = integrationStringTermsTable(text[["results_tsv"]]))
}

#' The enrichment kernel (KIT-141): the ranked list mapped to the database's identifiers, then for STRING the terms of
#' the saved answer (checked to be the answer to this list and these settings), for KEGG and Reactome the pathway test
#' on the saved tables. Refused when no pathway holds between `min_set_size` and `max_set_size` mapped features
#' (`integration_enrichment_no_testable_set`). The multiple-testing family is one factor, one data set and one database:
#' STRING adjusts within each of its term categories, KEGG and Reactome over every pathway tested (Benjamini-Hochberg);
#' nothing is adjusted across factors, data sets or databases.
#' @param object The inputs state (`features`, `inputs`).
#' @param parameters The declared settings and, beside them, `ranking_value`, `p` (the resolved settings), `response`
#'   (STRING) or `mapping_table` and `mapping_meta` (KEGG, Reactome).
#' @return The enrichment product.
#' @noRd
integrationEnrichmentKernel <- function(object, parameters) {
    p <- parameters$p
    ranking <- parameters$ranking_value
    mapped <- integrationEnrichmentMappedList(ranking, p, object@features[[p$view]])
    ranked <- mapped$table[c("identifier", "value")]
    out <- list(schema = INTEGRATION_ENRICHMENT_SCHEMA, factor = p$factor, view = p$view, omic = p$omic, database = p$database,
        settings = p[c("correlation_method", "min_samples", "id_type", "id_column")], direction = ranking$direction,
        ranking_counts = ranking$counts, mapping_counts = mapped$counts, mapping_rule = mapped$rule, mapping = mapped$mapping)
    if (identical(p$database, "string")) {
        response <- parameters$response
        request <- response$request
        same <- identical(as.character(request$identifiers_sha256), digest::digest(integrationStringIdentifiersText(ranked),
            algo = "sha256", serialize = FALSE)) && identical(as.character(request$species), as.character(p$string_species)) &&
            isTRUE(all.equal(as.numeric(request$ge_fdr), as.numeric(p$string_fdr))) &&
            identical(as.integer(request$rank_direction), as.integer(p$string_rank_direction))
        if (!same) integrationEnrichmentRefuse("integration_enrichment_response_mismatch")
        terms <- response$terms
        if (!is.data.frame(terms)) integrationEnrichmentRefuse("integration_enrichment_service_failed", "the results table could not be read")
        rownames(terms) <- NULL
        out$service <- list(name = "STRING", method = response$method, version = response$version, host = request$host,
            species = as.character(request$species), ge_fdr = request$ge_fdr, rank_direction = request$rank_direction,
            identifiers = request$identifiers, identifiers_sha256 = request$identifiers_sha256, queried_at = request$queried_at)
        out$family <- "terms of one STRING category, for one factor and one data set (STRING's Benjamini-Hochberg adjustment)"
        out$tests <- nrow(terms)
        out$terms <- terms
        return(out)
    }
    table <- parameters$mapping_table
    meta <- parameters$mapping_meta
    settings <- integrationGseaSettings(p_adjust_method = "BH", min_size = p$min_set_size, max_size = p$max_set_size,
        p_cutoff = 1, tie_jitter = 0)
    # KIT-183: a correlated test counts and tests the sets on U_cam (the complete identifiers), never on the ranked list.
    test_name <- p$test %||% "gsea"
    correlated <- test_name %in% INTEGRATION_ENRICHMENT_CORRELATED_TESTS
    # KIT-191: the permutation test's number of reorderings and, from the tab, its time limit (never declared).
    if (identical(test_name, "permutation")) {
        settings$permutations <- p$permutations %||% INTEGRATION_ENRICHMENT_PERMUTATIONS[[1L]]
        settings$time_limit <- p$time_limit
    }
    values <- if (correlated) parameters$values
    members <- if (correlated) rownames(values$y) else ranked$identifier
    if (identical(p$database, "kegg")) {
        links <- table[c("pathway", "compound")]
        organism_links <- grepl(paste0("^path:", p$kegg_organism, "[0-9]{5}$"), links$pathway)
        sizes <- table(links$pathway[organism_links & sub("^cpd:", "", links$compound) %in% members])
    } else {
        rows <- table[as.character(table$organism) == resolveReactomeOrganism(p$reactome_organism, as.character(table$organism)), ,
            drop = FALSE]
        chebi <- paste0("CHEBI:", suppressWarnings(as.integer(regmatches(as.character(rows$chebi_id),
            regexpr("[0-9]+", as.character(rows$chebi_id))))))
        has <- grepl("[0-9]", as.character(rows$chebi_id))
        pairs <- unique(data.frame(term = rows$reactome_id[has], gene = chebi, stringsAsFactors = FALSE))
        sizes <- table(pairs$term[pairs$gene %in% members])
    }
    testable <- sum(sizes >= p$min_set_size & sizes <= p$max_set_size)
    if (!testable) {
        integrationEnrichmentRefuse("integration_enrichment_no_testable_set", sprintf("%d pathways hold a mapped feature",
            length(sizes)), pathways_with_mapped = length(sizes))
    }
    test <- if (identical(p$database, "kegg")) {
        integrationKeggEnrichment(ranked, links, unique(table[c("pathway", "name")]), p$kegg_organism, meta$release, p$seed,
            settings, test_name, values)
    } else {
        integrationReactomeEnrichment(ranked, table[c("chebi_id", "reactome_id", "reactome_term", "organism")],
            p$reactome_organism, meta$release, p$seed, settings, test_name, values)
    }
    out$service <- list(name = test$database, release = meta$release, sha256 = meta$sha256, url = meta$url,
        downloaded_at = meta$downloaded_at, organism = test$organism, method = test$method, engine = test$engine,
        seed = test$seed, min_set_size = p$min_set_size, max_set_size = p$max_set_size, p_adjust_method = "BH",
        pathways_testable = testable)
    out$family <- sprintf("pathways of %s tested for one factor and one data set (Benjamini-Hochberg)", test$database)
    out$tests <- nrow(test$terms)
    out$terms <- test$terms
    if (correlated) {
        # KIT-183 (7, 8.8): the test, its settings, the counts, the design; no seed (nothing in camera is random).
        out$schema <- INTEGRATION_ENRICHMENT_CORRELATED_SCHEMA
        out$service$seed <- NULL
        out$service$test_settings <- test$settings$test_settings
        out$test <- test_name
        out$test_counts <- test$test_counts
        out$design_columns <- test$design_columns
        out$design_covariates <- test$design_covariates
        # KIT-191: the permutation test keeps its seed, its number of reorderings, the reorderings that left the factor
        # unchanged and the smallest p-value the call could give.
        if (identical(test_name, "permutation")) {
            out$service$seed <- test$settings$seed
            out$permutation <- test$permutation
        }
    }
    out
}

#' The enrichment declaration (KIT-141): product only, on the inputs state. Its identity is the state's, the ranking's
#' entity, the database, the identifier column, the saved answer's or tables' SHA-256 and every setting that changes the
#' result; the ranking, the answer and the tables are the product's dependencies.
#' @noRd
integrationEnrichmentOperation <- function() {
    newWorkflowOperation(INTEGRATION_ENRICHMENT_OPERATION, "IntegrationViewsData",
        fields = c("features", "inputs"), unit = "whole_matrix", output_class = NULL,
        parameters = c("ranking", "database", "id_column", "service_response", "mapping", "organism", "string_species",
            "string_fdr", "string_rank_direction", "min_set_size", "max_set_size", "seed"),
        product_type = "integration.enrichment", factors = list(product_schema = INTEGRATION_ENRICHMENT_SCHEMA),
        kernel = function(object, parameters, store) integrationEnrichmentKernel(object, parameters))
}

#' The correlated tests' declaration (KIT-183 section 7): product only, on the inputs state, with the design and the
#' sample match's fields (5.1 steps 1 and 2). Its identity is the state's, the ranking's entity, the database, the
#' identifier column, the tables' SHA-256, the organism, the set sizes, the test, its settings text, the values' slot
#' and the design's covariates (`"none"` in this ticket; KIT-184 states its own). `service_response`, the STRING
#' settings and `seed` are not declared: the correlated tests never run on STRING and nothing in camera is random
#' (KIT-191: the permutation test's seed and number of reorderings enter its identity through `test_settings`). The
#' `integration_enrichment` and ranking declarations are unchanged, so every pre-ranked product keeps its identity.
#' The kernel reads the data set's values through the store (`integrationEnrichmentTestValues`).
#' @noRd
integrationEnrichmentCorrelatedOperation <- function() {
    newWorkflowOperation(INTEGRATION_ENRICHMENT_CORRELATED_OPERATION, "IntegrationViewsData",
        fields = c("features", "design_matrix", "inputs"), unit = "whole_matrix", output_class = NULL,
        parameters = c("ranking", "database", "id_column", "mapping", "organism", "min_set_size", "max_set_size", "test",
            "test_settings", "values_slot", "design_covariates"),
        product_type = "integration.enrichment",
        factors = list(product_schema = INTEGRATION_ENRICHMENT_CORRELATED_SCHEMA),
        kernel = function(object, parameters, store) {
            p <- parameters$p
            mapped <- integrationEnrichmentMappedList(parameters$ranking_value, p, object@features[[p$view]])
            parameters$values <- integrationEnrichmentTestValues(object, store$manager, store$source,
                parameters$model_value, parameters$model_parameters, parameters$ranking_value, mapped, p$view)
            integrationEnrichmentKernel(object, parameters)
        })
}

#' Run an enrichment of a trained model's factor and accept it on the head (KIT-141). In the app process, in order:
#' the ranking (a declared operation; a saved ranking of the same model and settings is reused); for STRING the query,
#' with the key held only by this call (`query`), unless a saved answer is given (`response`, the replay); the answer or
#' the KEGG/Reactome tables saved as tables of the store; the enrichment (a declared operation, no network); the result
#' kept under the model's entry on the head. The key is never a setting, never in the operations' parameters, their
#' identities, the head or a file.
#' @param manager The integration owner with an accepted model.
#' @param parameters The request (`integrationResolveEnrichmentParameters`).
#' @param api_key The user's STRING key (STRING only; not needed with `response`).
#' @param mapping KEGG and Reactome: the downloaded tables (`integrationDownloadKeggMapping()` or
#'   `integrationDownloadReactomeMapping()`'s value).
#' @param response STRING: a saved answer to use instead of a query (`integrationServiceResponseFromTable()`'s form).
#' @param query The STRING client (`integrationStringRankedEnrichment`).
#' @param host The STRING address queried.
#' @param time_limit KIT-191 (D4 (c)): the longest predicted run of the permutation test accepted, in seconds (the
#'   tab: `INTEGRATION_ENRICHMENT_TAB_SECONDS`); NULL (scripts) runs it after a message giving the predicted time.
#' @return `list(run, ranking, entry)`: the enrichment's run (`status`, `value`, `reference`), the ranking's run and the
#'   entry kept on the head.
#' @noRd
runIntegrationEnrichment <- function(manager, parameters, api_key = NULL, mapping = NULL, response = NULL,
    query = integrationStringRankedEnrichment, host = INTEGRATION_STRING_DEFAULT_HOST, time_limit = NULL) {
    if (!usesWorkflowReferenceStates(manager)) {
        cli::cli_abort("the integration enrichment needs the integration owner", class = "saved_entity_scope")
    }
    withWorkflowStateScope(manager, function() {
        context <- manager$.context
        export <- context$snapshot$exports[[INTEGRATION_MODEL_OPERATION]]
        if (is.null(export)) integrationEnrichmentRefuse("integration_enrichment_no_model")
        source <- export$source
        value <- hydrateSavedWorkflowEntity(context$store, export$reference)$value
        integrationAssertModelProduct(value)
        object <- readWorkflowExecutionComponents(manager, source, fields = c("features", "design_matrix", "inputs"))
        n_samples <- if (is.character(parameters$view) && length(parameters$view) == 1L &&
            parameters$view %in% names(value$weights)) {
            length(integrationEnrichmentViewSamples(object, value, export$parameters, parameters$view))
        } else {
            NA_integer_
        }
        p <- integrationResolveEnrichmentParameters(parameters, value, object, n_samples)
        p$time_limit <- time_limit
        if (p$database %in% c("kegg", "reactome") && !(is.list(mapping) && identical(tolower(mapping$database),
            p$database) && is.data.frame(mapping$table))) {
            integrationEnrichmentRefuse("integration_enrichment_parameters_invalid",
                sprintf("download the %s tables first", if (p$database == "kegg") "KEGG" else "Reactome"))
        }
        ranking <- withCallingHandlers(runSavedOperation(manager, integrationEnrichmentRankingOperation(), source,
            list(model = export$reference$entity_id, factor = p$factor, view = p$view, correlation_method = p$correlation_method,
                min_samples = p$min_samples, direction_rule = INTEGRATION_ENRICHMENT_DIRECTION_RULE, model_value = value,
                model_parameters = export$parameters), dependencies = list(export$reference)),
            error = function(error) clearWorkflowStepProvenance(manager))
        clearWorkflowStepProvenance(manager)
        stored <- NULL
        extra <- list()
        if (identical(p$database, "string")) {
            if (is.null(response)) {
                mapped <- integrationEnrichmentMappedList(ranking$value, p, object@features[[p$view]])
                response <- query(mapped$table[c("identifier", "value")], api_key, p$string_species,
                    integrationEnrichmentStringOptions(p, host))
                api_key <- NULL
                recordWorkflowExternalResource(manager, "string-db", response$version$string_version %||% NA_character_,
                    "database_version", query_summary = sprintf("values/ranks enrichment, %s, %s, species %s, %d identifiers",
                        p$factor, p$view, p$string_species, as.integer(response$request$identifiers)))
            }
            stored <- integrationStoreServiceResponse(context, response)
            extra <- list(response = response)
            organism <- as.character(p$string_species)
        } else {
            table <- as.data.frame(mapping$table, stringsAsFactors = FALSE)
            stored <- list(reference = saveWorkflowOwnedTable(context, table), sha256 = integrationEnrichmentTableSha256(table))
            meta <- mapping[c("database", "release", "sha256", "url", "downloaded_at")]
            extra <- list(mapping_table = table, mapping_meta = meta)
            organism <- if (p$database == "kegg") p$kegg_organism else p$reactome_organism
        }
        correlated <- p$database != "string" && p$test %in% INTEGRATION_ENRICHMENT_CORRELATED_TESTS
        if (correlated) {
            # KIT-183: the correlated tests' own declaration (section 7); the pre-ranked test's is unchanged.
            declared <- list(ranking = ranking$reference$entity_id, database = p$database, id_column = p$id_column,
                mapping = stored$sha256, organism = organism, min_set_size = p$min_set_size, max_set_size = p$max_set_size,
                test = p$test, test_settings = integrationEnrichmentTestSettingsText(p$test, p), values_slot = "views",
                design_covariates = "none")
            operation <- integrationEnrichmentCorrelatedOperation()
            extra <- c(extra, list(model_value = value, model_parameters = export$parameters))
        } else {
            declared <- list(ranking = ranking$reference$entity_id, database = p$database, id_column = p$id_column,
                service_response = if (p$database == "string") stored$sha256 else "none",
                mapping = if (p$database == "string") "none" else stored$sha256, organism = organism,
                string_species = p$string_species, string_fdr = p$string_fdr, string_rank_direction = p$string_rank_direction,
                min_set_size = p$min_set_size, max_set_size = p$max_set_size, seed = p$seed)
            operation <- integrationEnrichmentOperation()
        }
        run <- withCallingHandlers(runSavedOperation(manager, operation, source,
            c(declared, list(p = p, ranking_value = ranking$value), extra),
            dependencies = list(ranking$reference, stored$reference)), error = function(error) clearWorkflowStepProvenance(manager))
        clearWorkflowStepProvenance(manager)
        entry <- acceptIntegrationEnrichment(manager, export, run, ranking, p, stored, extra$mapping_meta,
            strength = if (p$database != "string") integrationEnrichmentFactorStrength(value, p$factor, p$view))
        list(run = run, ranking = ranking, entry = entry)
    })
}

#' The head's key of one enrichment: factor, data set and database; KIT-183: and the test for a correlated test
#' (`factor|view|database|camera`), so that both tests' results can be kept side by side. The pre-ranked test and
#' STRING keep the key without a suffix.
#' @noRd
integrationEnrichmentKey <- function(factor, view, database, test = NULL) {
    key <- paste(factor, view, database, sep = "|")
    if (!identical(database, "string") && isTRUE(test %in% INTEGRATION_ENRICHMENT_CORRELATED_TESTS)) {
        key <- paste(key, test, sep = "|")
    }
    key
}

#' A factor's variance explained in a data set, and the largest of any factor there (KIT-183 5.5,
#' `enrich_factor_strength`): the model product's `variance_explained` (`r2_percent`, the variance of the data set's
#' selected, scaled features, INT-12).
#' @param value The model product.
#' @param factor,view The factor and data set.
#' @return `list(factor_r2_percent, largest_r2_percent)`, or NULL when the table lacks them.
#' @noRd
integrationEnrichmentFactorStrength <- function(value, factor, view) {
    v <- value$variance_explained
    if (!is.data.frame(v) || !all(c("factor", "view", "r2_percent") %in% names(v))) return(NULL)
    rows <- v[as.character(v$view) == view, , drop = FALSE]
    own <- as.numeric(rows$r2_percent[as.character(rows$factor) == factor])
    if (length(own) != 1L || is.na(own)) return(NULL)
    list(factor_r2_percent = own, largest_r2_percent = max(as.numeric(rows$r2_percent), na.rm = TRUE))
}

#' Keep an enrichment under the model's entry of the head (KIT-141): its product's and ranking's references, the saved
#' answer's or tables' reference and SHA-256, and the settings, as plain values (a reopened store refuses data frames
#' among the head's roots). A later enrichment of the same factor, data set and database replaces it; a new model
#' starts without any (`acceptIntegrationModel` keeps them only for the same model).
#' @return The entry.
#' @noRd
acceptIntegrationEnrichment <- function(manager, export, run, ranking, p, stored, mapping_meta = NULL, strength = NULL) {
    context <- manager$.context
    snapshot <- context$snapshot
    current <- snapshot$exports[[INTEGRATION_MODEL_OPERATION]]
    if (!identical(current$reference$entity_id, export$reference$entity_id)) return(invisible(NULL))
    for (reference in list(run$reference, ranking$reference)) {
        if (is.null(snapshot$products[[reference$entity_id]])) rootWorkflowProduct(manager, reference, list(export$source))
    }
    snapshot <- context$snapshot
    current <- snapshot$exports[[INTEGRATION_MODEL_OPERATION]]
    settings <- p[c("factor", "view", "omic", "database", "id_type", "id_column", "correlation_method", "min_samples",
        "string_species", "string_fdr", "string_rank_direction", "kegg_organism", "reactome_organism", "min_set_size",
        "max_set_size", "seed", "test", "test_ignored")]
    # KIT-191: a permutation entry keeps its number of reorderings (with its seed, above), so it replays identically.
    if (identical(p$test, "permutation")) settings$permutations <- p$permutations
    entry <- list(factor = p$factor, view = p$view, database = p$database, reference = run$reference,
        ranking = ranking$reference, stored = list(reference = stored$reference, sha256 = stored$sha256,
            kind = if (p$database == "string") "service_response" else "mapping"),
        mapping = if (!is.null(mapping_meta)) lapply(mapping_meta, as.character), settings = settings,
        saved_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
    # KIT-183: the factor's strength in this data set (KEGG and Reactome), and the values the model was trained on.
    entry$factor_r2_percent <- strength$factor_r2_percent
    entry$largest_r2_percent <- strength$largest_r2_percent
    entry$input_matrix <- as.character(export$parameters$input_matrix %||% "uncorrected")
    current$enrichment[[integrationEnrichmentKey(p$factor, p$view, p$database, p$test)]] <- entry
    snapshot$exports[[INTEGRATION_MODEL_OPERATION]] <- current
    publishWorkflowStateHead(context, snapshot)
    invisible(entry)
}

#' The enrichment results kept for the latest model, read from the store (KIT-141): per entry its product, ranking and
#' saved answer or tables. Read before the store is closed (the replay).
#' @param manager The integration owner.
#' @return Named list by entry key of `list(entry, value, ranking, stored)`; empty when none.
#' @noRd
integrationRecordedEnrichment <- function(manager) {
    entries <- manager$.context$snapshot$exports[[INTEGRATION_MODEL_OPERATION]]$enrichment
    if (!length(entries)) return(list())
    withWorkflowStateScope(manager, function() {
        store <- manager$.context$store
        lapply(entries, function(entry) list(entry = entry, value = hydrateSavedWorkflowEntity(store, entry$reference)$value,
            ranking = hydrateSavedWorkflowEntity(store, entry$ranking)$value,
            stored = materializeWorkflowTable(manager$.context, entry$stored$reference)))
    })
}

# --- KEGG and Reactome tables, downloaded in the app (INT-19; audit ENG-26) ---

#' One download (KIT-141): GET with a user agent, a time limit and retries of a busy service; any failure is a classed
#' refusal (`integration_enrichment_download_failed`).
#' @param url The address.
#' @param timeout_seconds The time allowed.
#' @return The body as raw bytes.
#' @noRd
integrationEnrichmentDownload <- function(url, timeout_seconds = 300) {
    requireOptionalPackage("httr2", "the pathway tables' download")
    on.exit(integrationForgetHttpRequest(), add = TRUE)
    request <- httr2::request(url) |>
        httr2::req_user_agent("MultiScholaR") |>
        httr2::req_timeout(timeout_seconds) |>
        httr2::req_retry(max_tries = 3L) |>
        httr2::req_error(is_error = \(response) FALSE)
    response <- tryCatch(httr2::req_perform(request), error = \(e) NULL)
    if (is.null(response)) integrationEnrichmentRefuse("integration_enrichment_download_failed", url)
    status <- httr2::resp_status(response)
    if (status >= 400L) integrationEnrichmentRefuse("integration_enrichment_download_failed", sprintf("%s: HTTP status %d", url, status))
    httr2::resp_body_raw(response)
}

#' Reactome's ChEBI-to-pathway table (KIT-141; ENG-26): the release is read before and after the file, and a release
#' that changed in between is refused (`integration_enrichment_release_changed`), so the file is labelled with the
#' release it came from.
#' @param download `integrationEnrichmentDownload` (tests pass a stand-in).
#' @return `list(database, release, sha256, url, downloaded_at, table)`: `table` with `chebi_id`, `reactome_id`,
#'   `reactome_term`, `organism`; `sha256` of the file's bytes.
#' @noRd
integrationDownloadReactomeMapping <- function(download = integrationEnrichmentDownload) {
    version_url <- "https://reactome.org/ContentService/data/database/version"
    url <- "https://reactome.org/download/current/ChEBI2Reactome_All_Levels.txt"
    release <- function() trimws(rawToChar(download(version_url)))
    before <- release()
    bytes <- download(url)
    after <- release()
    if (!identical(before, after) || !grepl("^[0-9]+$", before)) {
        integrationEnrichmentRefuse("integration_enrichment_release_changed", sprintf("Reactome %s, then %s", before, after))
    }
    table <- utils::read.delim(text = rawToChar(bytes), header = FALSE, sep = "\t", quote = "", comment.char = "",
        stringsAsFactors = FALSE, colClasses = "character")
    if (ncol(table) < 6L) integrationEnrichmentRefuse("integration_enrichment_download_failed", "the Reactome file has another layout")
    table <- data.frame(chebi_id = table[[1L]], reactome_id = table[[2L]], reactome_term = table[[4L]], organism = table[[6L]],
        stringsAsFactors = FALSE)
    list(database = "Reactome", release = paste("Reactome", before), sha256 = digest::digest(bytes, algo = "sha256",
        serialize = FALSE), url = url, downloaded_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), table = table)
}

#' KEGG's compound-to-pathway table for an organism (KIT-141; ENG-26): the release (`info/kegg`) read before and after,
#' the reference pathways' compounds (`link/compound/pathway`) kept for the pathways the organism has
#' (`list/pathway/<organism>`), named as KEGG names them for the organism. KEGG's REST service is free for academic use
#' only (the tab says so before the download).
#' @param organism KEGG organism code.
#' @param download `integrationEnrichmentDownload` (tests pass a stand-in).
#' @return `list(database, release, sha256, url, downloaded_at, table)`: `table` with `pathway` (`path:<organism><5
#'   digits>`), `compound` (`cpd:C.....`), `name`; `sha256` of the two tables' bytes.
#' @noRd
integrationDownloadKeggMapping <- function(organism, download = integrationEnrichmentDownload) {
    if (!is.character(organism) || length(organism) != 1L || !grepl("^[a-z]{3,4}$", organism)) {
        integrationEnrichmentRefuse("integration_enrichment_parameters_invalid", "the KEGG organism must be a KEGG organism code, such as hsa")
    }
    release <- function() {
        text <- rawToChar(download("https://rest.kegg.jp/info/kegg"))
        hit <- regmatches(text, regexpr("Release [^,\r\n]+", text))
        if (length(hit)) trimws(hit) else NA_character_
    }
    lines <- function(bytes) {
        x <- strsplit(rawToChar(bytes), "\n", fixed = TRUE)[[1L]]
        x <- x[nzchar(trimws(x))]
        do.call(rbind, lapply(strsplit(x, "\t", fixed = TRUE), function(row) row[1:2]))
    }
    before <- release()
    link_url <- "https://rest.kegg.jp/link/compound/pathway"
    list_url <- paste0("https://rest.kegg.jp/list/pathway/", organism)
    link_bytes <- download(link_url)
    list_bytes <- download(list_url)
    after <- release()
    if (is.na(before) || !identical(before, after)) {
        integrationEnrichmentRefuse("integration_enrichment_release_changed", sprintf("KEGG %s, then %s", before, after))
    }
    links <- lines(link_bytes)
    pathway_col <- if (all(startsWith(links[, 1L], "path:"))) 1L else 2L
    pathway <- links[, pathway_col]
    compound <- links[, 3L - pathway_col]
    named <- lines(list_bytes)
    number <- function(id) sub("^.*([0-9]{5})$", "\\1", id)
    names_by_number <- stats::setNames(named[, 2L], number(named[, 1L]))
    keep <- number(pathway) %in% names(names_by_number) & startsWith(compound, "cpd:")
    table <- data.frame(pathway = paste0("path:", organism, number(pathway[keep])), compound = compound[keep],
        name = unname(names_by_number[number(pathway[keep])]), stringsAsFactors = FALSE)
    table <- unique(table)
    rownames(table) <- NULL
    list(database = "KEGG", release = paste("KEGG", before), sha256 = digest::digest(c(link_bytes, list_bytes), algo = "sha256",
        serialize = FALSE), url = paste(link_url, list_url), downloaded_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
        table = table)
}

# --- The record (audit ENG-15) ---

#' The enrichment entries of the record (KIT-141): per entry the settings, the direction, the counts, the service and
#' release, the family and number of tests, and the terms (as rows). Read from the store; plain values.
#' @param context The owner's context.
#' @param export The model's entry of the head.
#' @return A list of entries (empty when none).
#' @noRd
integrationRecordEnrichmentBlock <- function(context, export) {
    entries <- export$enrichment %||% list()
    unname(lapply(entries, function(entry) {
        value <- tryCatch(hydrateSavedWorkflowEntity(context$store, entry$reference)$value, error = function(error) NULL)
        if (!is.list(value)) return(list(factor = entry$factor, view = entry$view, database = entry$database, complete = FALSE))
        terms <- as.data.frame(value$terms)
        list(factor = value$factor, view = value$view, omic = value$omic, database = value$database, complete = TRUE,
            settings = entry$settings, direction = value$direction, ranking_counts = value$ranking_counts,
            mapping_counts = value$mapping_counts, mapping_rule = value$mapping_rule, service = value$service,
            stored_sha256 = entry$stored$sha256, family = value$family, tests = value$tests, entity_id = entry$reference$entity_id,
            # KIT-183: the test (absent: the pre-ranked test), its counts and design, the factor's strength, the values.
            test = value$test %||% (if (!identical(value$database, "string")) "gsea"), test_counts = value$test_counts,
            design_columns = value$design_columns, design_covariates = value$design_covariates,
            # KIT-191: the permutation test's seed, reorderings, unchanged draws and smallest p-value.
            permutation = value$permutation,
            factor_r2_percent = entry$factor_r2_percent, largest_r2_percent = entry$largest_r2_percent,
            input_matrix = entry$input_matrix,
            terms = unname(lapply(seq_len(nrow(terms)), function(i) lapply(as.list(terms[i, , drop = FALSE]), function(x)
                if (is.factor(x)) as.character(x) else x))))
    }))
}

#' The enrichment entries as the record's two sheets (KIT-141): one row per enrichment, and one row per term.
#' @param block The record's integration block.
#' @return `list(summary, terms)`, headers from `integrationRecordTexts()`.
#' @noRd
integrationRecordEnrichmentTables <- function(block, texts = integrationRecordTexts()) {
    s <- integrationRecordScalar
    entries <- Filter(function(entry) isTRUE(as.logical(s(entry$complete))), block$enrichment %||% list())
    text <- function(value) { value <- s(value); if (is.na(value)) texts[["not_recorded"]] else as.character(value) }
    summary <- lapply(entries, function(entry) {
        service <- entry$service
        release <- if (identical(s(entry$database), "string")) sprintf(texts[["enrich_string_release"]], text(service$version$string_version),
            text(service$host)) else text(service$release)
        data.frame(factor = text(entry$factor), view = text(entry$view), database = integrationEnrichmentDatabaseName(s(entry$database)),
            release = release, ranking = integrationRecordEnrichmentRanking(entry, texts),
            direction = sprintf(texts[["enrich_direction_value"]], text(entry$direction$feature_id), text(entry$direction$view),
                as.integer(s(entry$direction$sign))),
            features = sprintf(texts[["enrich_counts_value"]], as.integer(s(entry$ranking_counts$features)),
                as.integer(s(entry$ranking_counts$ranked)), as.integer(s(entry$mapping_counts$mapped))),
            tests = as.character(as.integer(s(entry$tests))), family = text(entry$family), checksum = text(entry$stored_sha256),
            strength = integrationEnrichmentStrengthText(entry, texts), stringsAsFactors = FALSE)
    })
    # Without an enrichment, one row saying so.
    none <- function(n) as.data.frame(as.list(c(texts[["none_recorded"]], rep("", n - 1L))), stringsAsFactors = FALSE)
    summary <- if (length(summary)) do.call(rbind, summary) else none(11L)
    names(summary) <- unname(texts[c("col_factor", "col_view", "col_enrich_database", "col_enrich_release", "col_enrich_ranking",
        "col_enrich_direction", "col_enrich_features", "col_enrich_tests", "col_enrich_family", "col_checksum",
        "col_enrich_strength")])
    terms <- lapply(entries, function(entry) {
        rows <- integrationRecordFrame(entry$terms)
        if (!nrow(rows)) return(NULL)
        string <- identical(s(entry$database), "string")
        pick <- function(name) if (name %in% names(rows)) as.character(rows[[name]]) else rep(NA_character_, nrow(rows))
        data.frame(factor = text(entry$factor), view = text(entry$view), database = integrationEnrichmentDatabaseName(s(entry$database)),
            term = if (string) pick("term") else pick("ID"), description = if (string) pick("description") else pick("Description"),
            category = if (string) pick("category") else "", score = if (string) pick("enrichment_score") else
                if ("effect" %in% names(rows)) pick("effect") else pick("NES"),
            fdr = if (string) pick("false_discovery_rate") else pick("p.adjust"), stringsAsFactors = FALSE)
    })
    terms <- do.call(rbind, Filter(Negate(is.null), terms))
    if (is.null(terms)) terms <- none(8L)
    names(terms) <- unname(texts[c("col_factor", "col_view", "col_enrich_database", "col_enrich_term", "col_enrich_description",
        "col_enrich_category", "col_enrich_score", "col_enrich_fdr")])
    rownames(summary) <- NULL
    rownames(terms) <- NULL
    list(summary = summary, terms = terms)
}

#' KIT-183: an entry's "How features were ranked" cell: STRING keeps today's text (with `enrich_string_test_note` when
#' a test was given and ignored); the pre-ranked test keeps today's; a correlated test names the regression, the
#' number of samples and of identifiers used, and the test.
#' @param entry An entry of the record's enrichment block (or of the head, with `settings`).
#' @noRd
integrationRecordEnrichmentRanking <- function(entry, texts = integrationRecordTexts()) {
    s <- integrationRecordScalar
    ranking <- sprintf(texts[["enrich_ranking_value"]], as.character(s(entry$settings$correlation_method)),
        as.integer(s(entry$settings$min_samples)))
    if (identical(s(entry$database), "string")) {
        if (isTRUE(as.logical(s(entry$settings$test_ignored)))) ranking <- paste0(ranking, "; ", texts[["enrich_string_test_note"]])
        return(ranking)
    }
    test <- s(entry$test)
    if (!isTRUE(test %in% INTEGRATION_ENRICHMENT_CORRELATED_TESTS)) return(ranking)
    sprintf(texts[[paste0("enrich_ranking_value_", test)]], as.integer(s(entry$test_counts$samples)),
        as.integer(s(entry$test_counts$universe)))
}

#' KIT-183: a percentage of variance explained as the texts print it: one decimal, below 0.05 % "less than 0.1 %".
#' @noRd
integrationEnrichmentPercent <- function(x) {
    x <- as.numeric(x)
    if (is.na(x)) return(NA_character_)
    if (x < 0.05) "less than 0.1 %" else sprintf("%.1f %%", x)
}

#' KIT-183: `enrich_factor_strength` with an entry's two numbers; empty for STRING and for an entry without them.
#' @param entry An entry with `database`, `factor_r2_percent` and `largest_r2_percent`.
#' @param texts Texts holding `enrich_factor_strength`.
#' @noRd
integrationEnrichmentStrengthText <- function(entry, texts) {
    s <- integrationRecordScalar
    own <- s(entry$factor_r2_percent)
    largest <- s(entry$largest_r2_percent)
    if (identical(s(entry$database), "string") || is.na(own) || is.na(largest)) return("")
    # KIT-191: the figure of the entry's own test and the data set's kind (row E of the evidence table).
    test <- integrationEnrichmentEntryTest(entry)
    omic <- s(entry$omic)
    if (is.na(omic)) omic <- s(entry$settings$omic)
    kind <- integrationEnrichmentEvidenceKind(omic)
    sprintf(texts[["enrich_factor_strength"]], integrationEnrichmentPercent(own), integrationEnrichmentPercent(largest),
        texts[[paste0("enrich_strength_kind_", if (identical(omic, "lipidomics")) "lipids" else kind)]],
        texts[[paste0("enrich_strength_terms_", kind)]], integrationEnrichmentEvidenceValue(kind, test, "E", short = TRUE),
        texts[[paste0("enrich_test_phrase_", test)]])
}

#' KIT-191: an entry's test among the offered ones, from the record's block (`test`) or the head (`settings$test`); an
#' entry without one (made before KIT-183) is the pre-ranked test.
#' @noRd
integrationEnrichmentEntryTest <- function(entry) {
    s <- integrationRecordScalar
    test <- s(entry$test)
    if (is.na(test)) test <- s(entry$settings$test)
    if (is.na(test)) return("gsea")
    if (!test %in% INTEGRATION_ENRICHMENT_TESTS_OFFERED) return("camera")
    test
}

#' KIT-191: the kind of simulated data whose figures apply to a data set of `omic` (lipids take the metabolites'
#' figures, labelled: no lipid data were simulated).
#' @noRd
integrationEnrichmentEvidenceKind <- function(omic) {
    if (isTRUE(omic %in% c("proteomics", "transcriptomics"))) "genes_proteins" else "metabolites"
}

#' KIT-191: one cell of `INTEGRATION_ENRICHMENT_TEST_EVIDENCE` (NA when absent); `short` drops the 95 % interval and
#' puts the pre-ranked test's figure on correlated terms in brackets.
#' @noRd
integrationEnrichmentEvidenceValue <- function(kind, test, measure, samples = 12L, short = FALSE) {
    e <- INTEGRATION_ENRICHMENT_TEST_EVIDENCE
    x <- e$value_text[e$kind == kind & e$test == test & e$samples == samples & e$measure == measure]
    if (length(x) != 1L) return(NA_character_)
    if (short) {
        x <- sub(" \\((largest: )?95 % interval [^)]*\\)", "", x)
        x <- sub("; (terms whose features vary together: .*)$", " (\\1)", x)
    }
    x
}

#' KIT-191 (5.4): the rows of the evidence table the tab shows for a kind of data set, at 12 biological samples.
#' @noRd
integrationEnrichmentEvidenceRows <- function(kind) {
    if (identical(kind, "genes_proteins")) c("A", "B", "D", "E", "G", "H") else c("A", "B", "D", "E", "C3", "I")
}

#' KIT-191 (5.4): the record's sheet of the evidence table, every row of both kinds, when the record holds a KEGG,
#' Reactome or Gene Ontology entry.
#' @param texts Texts holding `col_evidence_*` and the measures' and kinds' names.
#' @noRd
integrationRecordEvidenceTable <- function(texts = integrationRecordTexts()) {
    e <- INTEGRATION_ENRICHMENT_TEST_EVIDENCE
    measure <- ifelse(e$measure == "E", paste0("enrich_evidence_E_", e$kind), paste0("enrich_evidence_", e$measure))
    out <- data.frame(unname(texts[paste0("enrich_evidence_kind_", e$kind)]), unname(texts[paste0("enrich_test_phrase_",
        e$test)]), e$samples, unname(texts[measure]), e$value_text, e$applies_to, e$source, stringsAsFactors = FALSE)
    names(out) <- unname(texts[c("col_evidence_kind", "col_evidence_test", "col_evidence_samples", "col_evidence_measure",
        "col_evidence_value", "col_evidence_applies", "col_evidence_source")])
    out
}

#' A database's name as shown.
#' @noRd
integrationEnrichmentDatabaseName <- function(database) {
    switch(as.character(database), string = "STRING", kegg = "KEGG", reactome = "Reactome", as.character(database))
}

#' The methods paragraph's sentences on the enrichment (KIT-141), from the record: per database used, the ranking, the
#' direction, the mapping and the service or tables with their release, and the family of the multiple-testing
#' correction.
#' @param block The record's integration block.
#' @return A character vector (empty without enrichment).
#' @noRd
integrationMethodsEnrichment <- function(block, texts = integrationRecordTexts()) {
    s <- integrationRecordScalar
    entries <- Filter(function(entry) isTRUE(as.logical(s(entry$complete))), block$enrichment %||% list())
    if (!length(entries)) return(character())
    first <- entries[[1L]]
    sentences <- sprintf(texts[["methods_enrich_ranking"]], as.character(s(first$settings$correlation_method)),
        as.integer(s(first$settings$min_samples)))
    pathways <- FALSE
    for (entry in entries) {
        service <- entry$service
        database <- s(entry$database)
        test <- s(entry$test)
        pathways <- pathways || !identical(database, "string")
        sentences <- c(sentences, if (identical(database, "string")) sprintf(texts[["methods_enrich_string"]], s(entry$factor),
            s(entry$view), as.integer(s(entry$mapping_counts$mapped)), as.character(s(service$version$string_version)),
            as.character(s(service$species)), as.character(s(service$queried_at))) else
            if (isTRUE(test %in% INTEGRATION_ENRICHMENT_CORRELATED_TESTS)) {
                # KIT-183: the correlated test's sentence, composed from its name and settings fragments; KIT-191: the
                # permutation test's name holds its number of reorderings (twice) and its settings its seed first.
                counts <- entry$test_counts
                name <- texts[[paste0("methods_test_name_", test)]]
                numbers <- list(as.integer(s(counts$samples)), as.integer(s(counts$universe)),
                    as.integer(s(entry$mapping_counts$mapped)))
                if (identical(test, "permutation")) {
                    reorderings <- as.integer(s(entry$permutation$permutations))
                    name <- sprintf(name, reorderings, reorderings)
                    numbers <- c(list(as.integer(s(service$seed))), numbers)
                }
                sprintf(texts[[paste0("methods_enrich_", test)]], s(entry$factor), s(entry$view),
                    integrationEnrichmentDatabaseName(database), name,
                    do.call(sprintf, c(list(texts[[paste0("methods_test_settings_", test)]]), numbers)),
                    as.character(s(service$release)), as.character(s(service$organism)), as.integer(s(service$min_set_size)),
                    as.integer(s(service$max_set_size)), as.integer(s(service$pathways_testable)))
            } else
            sprintf(texts[["methods_enrich_pathway"]], s(entry$factor), s(entry$view), integrationEnrichmentDatabaseName(database),
                as.integer(s(entry$mapping_counts$mapped)), as.character(s(service$release)), as.character(s(service$organism)),
                as.integer(s(service$min_set_size)), as.integer(s(service$max_set_size)), as.integer(s(service$pathways_testable)),
                as.integer(s(service$seed))))
    }
    # KIT-191: a factor, data set and database with more than one test's entry says so, once.
    sentences <- c(sentences, integrationMethodsSeveralTests(entries, texts))
    # KIT-183: the limits, once, when the record holds a KEGG or Reactome entry; KIT-191: one clause per test used and
    # kind of data set.
    c(sentences, if (pathways) integrationMethodsEnrichmentLimits(entries, texts), texts[["methods_enrich_family"]])
}

#' KIT-191 (5.4): the limits paragraph of the Methods, composed from the tests and kinds of data set of the record's
#' KEGG, Reactome (and Gene Ontology) entries: the opening, one clause per test and kind with its rows A, B and E of
#' the evidence table, the closing.
#' @param entries The record's complete enrichment entries.
#' @noRd
integrationMethodsEnrichmentLimits <- function(entries, texts = integrationRecordTexts()) {
    s <- integrationRecordScalar
    pathway <- Filter(function(entry) !identical(s(entry$database), "string"), entries)
    pairs <- unique(do.call(rbind, lapply(pathway, function(entry) {
        omic <- s(entry$omic)
        data.frame(test = integrationEnrichmentEntryTest(entry), kind = integrationEnrichmentEvidenceKind(omic),
            lipids = identical(omic, "lipidomics"), stringsAsFactors = FALSE)
    })))
    order_tests <- match(pairs$test, INTEGRATION_ENRICHMENT_TESTS_OFFERED)
    pairs <- pairs[order(order_tests, pairs$kind, pairs$lipids), , drop = FALSE]
    clauses <- vapply(seq_len(nrow(pairs)), function(i) {
        test <- pairs$test[[i]]
        kind <- pairs$kind[[i]]
        figure <- function(measure) integrationEnrichmentEvidenceValue(kind, test, measure, short = TRUE)
        sprintf(texts[["methods_enrich_limits_clause"]], texts[[paste0("enrich_test_phrase_", test)]],
            texts[[paste0("enrich_strength_kind_", if (pairs$lipids[[i]]) "lipids" else kind)]], figure("A"), figure("B"),
            figure("E"))
    }, character(1))
    paste0(texts[["methods_enrich_limits_tests"]], " ", paste(clauses, collapse = "; "), ". ",
        texts[["methods_enrich_limits_end"]])
}

#' KIT-191 (5.4, A15): one sentence per factor, data set and database with more than one test's entry.
#' @noRd
integrationMethodsSeveralTests <- function(entries, texts = integrationRecordTexts()) {
    s <- integrationRecordScalar
    pathway <- Filter(function(entry) !identical(s(entry$database), "string"), entries)
    if (!length(pathway)) return(character())
    key <- vapply(pathway, function(entry) paste(s(entry$factor), s(entry$view), s(entry$database), sep = "|"), character(1))
    out <- character()
    for (k in unique(key)) {
        group <- pathway[key == k]
        tests <- unique(vapply(group, integrationEnrichmentEntryTest, character(1)))
        if (length(tests) < 2L) next
        tests <- tests[order(match(tests, INTEGRATION_ENRICHMENT_TESTS_OFFERED))]
        first <- group[[1L]]
        out <- c(out, sprintf(texts[["methods_enrich_several_tests"]], s(first$factor), s(first$view),
            integrationEnrichmentDatabaseName(s(first$database)), length(tests),
            paste(unname(texts[paste0("enrich_test_phrase_", tests)]), collapse = ", ")))
    }
    out
}

# --- The replay (audit ENG-15) ---

#' Repeat the recorded enrichments on a replayed model (KIT-141): each recorded enrichment is run again on the replay's
#' store with its recorded settings, by default from the recorded STRING answer or KEGG/Reactome tables (no network);
#' with `mode = "query_again"` STRING is queried again (`api_key`) and the tables downloaded again (`download`), and the
#' differences reported. The ranking and the terms are compared with the recorded ones (`integrationCompareTable`).
#' @param manager The replay's owner (a model accepted).
#' @param recorded `integrationRecordedEnrichment()`'s value.
#' @param mode `recorded` or `query_again`.
#' @param api_key The STRING key (`query_again` only).
#' @param query,download The STRING client and the downloader (tests pass stand-ins).
#' @param tolerance The relative tolerance.
#' @return `data.frame(enrichment, mode, ranking_equal, terms_equal, terms_identical, terms_only_recorded,
#'   terms_only_repeated, note)`.
#' @noRd
integrationReplayEnrichment <- function(manager, recorded, mode = c("recorded", "query_again"), api_key = NULL,
    query = integrationStringRankedEnrichment, download = integrationEnrichmentDownload, tolerance = 1e-5) {
    mode <- match.arg(mode)
    rows <- lapply(names(recorded), function(key) {
        item <- recorded[[key]]
        settings <- item$entry$settings
        parameters <- list(factor = settings$factor, view = settings$view, database = settings$database)
        for (name in names(integrationEnrichmentParameterSpec())) parameters[[name]] <- settings[[name]]
        if (!is.na(settings$id_column %||% NA)) parameters$id_column <- settings$id_column
        database <- settings$database
        # KIT-183: an entry replays with the test it was made with; a KEGG or Reactome entry without one (made before
        # KIT-183) is a pre-ranked entry, never the current default. STRING ignores the test (its stored NA removed).
        parameters$test <- if (identical(database, "string")) NULL else settings$test %||% "gsea"
        mapping <- NULL
        response <- NULL
        if (identical(database, "string")) {
            if (mode == "recorded") response <- integrationServiceResponseFromTable(item$stored)
        } else {
            mapping <- if (mode == "recorded") c(lapply(item$entry$mapping, identity), list(table = item$stored)) else
                if (database == "kegg") integrationDownloadKeggMapping(settings$kegg_organism, download) else
                    integrationDownloadReactomeMapping(download)
        }
        outcome <- tryCatch(runIntegrationEnrichment(manager, parameters, api_key = api_key, mapping = mapping,
            response = response, query = query, host = item$value$service$host %||% INTEGRATION_STRING_DEFAULT_HOST),
            integration_enrichment_refused = function(condition) condition)
        ranking_cmp <- NULL
        if (inherits(outcome, "condition")) {
            return(data.frame(enrichment = key, mode = mode, ranking_equal = NA, terms_equal = FALSE, terms_identical = FALSE,
                terms_only_recorded = NA_integer_, terms_only_repeated = NA_integer_, note = conditionMessage(outcome),
                stringsAsFactors = FALSE))
        }
        ranking_cmp <- integrationCompareTable("ranking", item$ranking$ranking, outcome$ranking$value$ranking, tolerance)
        a <- as.data.frame(item$value$terms)
        b <- as.data.frame(outcome$run$value$terms)
        id <- if (identical(database, "string")) "term" else "ID"
        terms_cmp <- integrationCompareTable("terms", a, b, tolerance)
        data.frame(enrichment = key, mode = mode, ranking_equal = ranking_cmp$within_tolerance,
            terms_equal = terms_cmp$within_tolerance, terms_identical = terms_cmp$identical,
            terms_only_recorded = length(setdiff(a[[id]], b[[id]])), terms_only_repeated = length(setdiff(b[[id]], a[[id]])),
            note = "", stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}
