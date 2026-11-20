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
#' @noRd
integrationEnrichmentRefuse <- function(class, detail = NULL, key = NULL, ...) {
    sentence <- INTEGRATION_ENRICHMENT_REFUSALS[[class]]
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
#' `withr::with_seed(seed)`, so the caller's random state is untouched and the same input, seed and settings give the
#' same table. An error of the test is a classed refusal (`integration_enrichment_method_failed`, with the engine's
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
    })
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

#' The pathway test over pathway sets restricted to the ranked identifiers; refused when none overlaps.
#' @noRd
integrationEnrichmentPathwayTest <- function(ranked, term2gene, term2name, seed, settings) {
    term2gene <- unique(term2gene[term2gene$gene %in% ranked$identifier, c("term", "gene"), drop = FALSE])
    rownames(term2gene) <- NULL
    if (nrow(term2gene) == 0L) integrationEnrichmentRefuse("integration_enrichment_no_overlap")
    term2name <- unique(term2name[term2name$term %in% term2gene$term, c("term", "name"), drop = FALSE])
    rownames(term2name) <- NULL
    test <- integrationEnrichmentGsea(ranked, term2gene, term2name, seed, settings)
    c(test, list(term2gene = term2gene, term2name = term2name))
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
#' @return The pathway test's list plus `database`, `organism`, `release`, `term2gene`, `term2name`.
#' @noRd
integrationKeggEnrichment <- function(ranked, pathway_compounds, pathway_names, organism, release, seed,
    settings = integrationGseaSettings()) {
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
    test <- integrationEnrichmentPathwayTest(ranked, term2gene, term2name, seed, settings)
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
#' @return The pathway test's list plus `database`, `organism`, `release`, `term2gene`, `term2name`.
#' @noRd
integrationReactomeEnrichment <- function(ranked, chebi_to_reactome, organism, release, seed,
    settings = integrationGseaSettings()) {
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
    test <- integrationEnrichmentPathwayTest(ranked, term2gene, term2name, seed, settings)
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

