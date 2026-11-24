# KIT-175 gate F3 (runner only, network open): does a UniProt search of 25 accessions, sent exactly
# as directUniprotDownload sends it, ever return more hits than rows (a cut first page), and does the
# TSV answer carry X-Total-Results? Writes counts, header presence and the release only; no accession
# is written or printed (UniProt content CC BY 4.0; licence guard ENG-07).
#
# Usage (from the package root): Rscript docs/qa/workflow-execution/kit/kit-175/f3-uniprot-pages.R <out dir>

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[[1L]] else "f3-out"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fields <- paste0(
  "accession,id,protein_name,gene_names,organism_name,length,go_id,reviewed,",
  "protein_existence,annotation_score"
)
accession_pattern <- "^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"

header_or_na <- function(response, name) {
  value <- httr::headers(response)[[name]]
  if (is.null(value)) NA_character_ else value
}

# One request for the list of a proteome, sorted by accession (the stream endpoint returns all rows).
fetch_list <- function(proteome, reviewed_only) {
  query <- paste0("proteome:", proteome, if (reviewed_only) " AND reviewed:true" else "")
  response <- httr::GET(
    url = "https://rest.uniprot.org/uniprotkb/stream",
    query = list(query = query, format = "tsv", fields = "accession,organism_id"),
    httr::timeout(120)
  )
  if (httr::status_code(response) != 200L) {
    stop(sprintf("list request for %s failed with status %d", proteome, httr::status_code(response)))
  }
  body <- httr::content(response, "text", encoding = "UTF-8")
  table <- utils::read.delim(text = body, sep = "\t", quote = "", stringsAsFactors = FALSE)
  table <- table[order(table$Entry), , drop = FALSE]
  table
}

query_chunks <- function(dataset, accessions, taxon_id) {
  chunks <- split(accessions, ceiling(seq_along(accessions) / 25))
  rows <- vector("list", length(chunks))
  for (i in seq_along(chunks)) {
    chunk <- chunks[[i]]
    query <- paste0("(", paste(chunk, collapse = " OR "), ") AND organism_id:", taxon_id)
    response <- tryCatch(
      httr::GET(
        url = "https://rest.uniprot.org/uniprotkb/search",
        query = list(query = query, format = "tsv", fields = fields),
        httr::timeout(30)
      ),
      error = function(e) e
    )
    Sys.sleep(1)
    if (inherits(response, "error")) {
      rows[[i]] <- data.frame(
        dataset = dataset, request = i, requested = length(chunk), status = NA_integer_,
        transport_error = TRUE, rows = NA_integer_, x_total_results = NA_integer_,
        x_total_results_present = NA, link_present = NA, link_next = NA, retry_after = NA_character_,
        entries_outside_request = NA_integer_, requested_not_returned = NA_integer_,
        release = NA_character_, stringsAsFactors = FALSE
      )
      next
    }
    status <- httr::status_code(response)
    total <- header_or_na(response, "x-total-results")
    link <- header_or_na(response, "link")
    n_rows <- NA_integer_
    outside <- NA_integer_
    missing_n <- NA_integer_
    if (status == 200L) {
      body <- httr::content(response, "text", encoding = "UTF-8")
      table <- tryCatch(
        suppressWarnings(utils::read.delim(text = body, sep = "\t", quote = "", stringsAsFactors = FALSE)),
        error = function(e) data.frame()
      )
      n_rows <- nrow(table)
      if ("Entry" %in% names(table)) {
        outside <- sum(!table$Entry %in% chunk)
        missing_n <- sum(!chunk %in% table$Entry)
      }
    }
    rows[[i]] <- data.frame(
      dataset = dataset, request = i, requested = length(chunk), status = status,
      transport_error = FALSE, rows = n_rows,
      x_total_results = suppressWarnings(as.integer(total)),
      x_total_results_present = !is.na(total), link_present = !is.na(link),
      link_next = !is.na(link) && grepl('rel="next"', link, fixed = TRUE),
      retry_after = header_or_na(response, "retry-after"),
      entries_outside_request = outside, requested_not_returned = missing_n,
      release = header_or_na(response, "x-uniprot-release"), stringsAsFactors = FALSE
    )
    message(sprintf("%s request %d of %d: status %d", dataset, i, length(chunks), status))
  }
  do.call(rbind, rows)
}

human_list <- fetch_list("UP000005640", reviewed_only = TRUE)
ecoli_list <- fetch_list("UP000000625", reviewed_only = FALSE)

human <- utils::head(human_list$Entry, 2000L)
ecoli <- utils::head(ecoli_list$Entry, 1000L)
ecoli_taxa <- table(utils::head(ecoli_list$Organism..ID., 1000L))

list_summary <- data.frame(
  dataset = c("human_UP000005640_reviewed", "ecoli_UP000000625"),
  list_rows = c(nrow(human_list), nrow(ecoli_list)),
  taken = c(length(human), length(ecoli)),
  taken_matching_accession_pattern = c(sum(grepl(accession_pattern, human)), sum(grepl(accession_pattern, ecoli))),
  organism_ids_of_taken = c(
    paste(names(table(utils::head(human_list$Organism..ID., 2000L))), collapse = ";"),
    paste(names(ecoli_taxa), collapse = ";")
  ),
  stringsAsFactors = FALSE
)
utils::write.table(list_summary, file.path(out_dir, "f3-lists.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

requests <- rbind(
  query_chunks("human_UP000005640_reviewed", human, 9606),
  query_chunks("ecoli_UP000000625", ecoli, 83333)
)
utils::write.table(requests, file.path(out_dir, "f3-requests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

ok <- requests[!requests$transport_error & requests$status == 200L, , drop = FALSE]
summary_lines <- c(
  sprintf("requests: %d (status 200: %d; other status: %d; transport errors: %d)",
          nrow(requests), nrow(ok), sum(!requests$transport_error & requests$status != 200L),
          sum(requests$transport_error)),
  sprintf("200 answers with X-Total-Results: %d of %d", sum(ok$x_total_results_present), nrow(ok)),
  sprintf("200 answers where X-Total-Results exceeds rows: %d",
          sum(ok$x_total_results_present & ok$x_total_results > ok$rows, na.rm = TRUE)),
  sprintf("200 answers with a Link header: %d (rel=next: %d)", sum(ok$link_present), sum(ok$link_next)),
  sprintf("answers with Retry-After: %d (values: %s)", sum(!is.na(requests$retry_after)),
          paste(unique(stats::na.omit(requests$retry_after)), collapse = ";")),
  sprintf("returned entries outside the requested accessions: %d (max per request %s)",
          sum(ok$entries_outside_request, na.rm = TRUE),
          max(c(0L, ok$entries_outside_request), na.rm = TRUE)),
  sprintf("requested accessions not returned: %d", sum(ok$requested_not_returned, na.rm = TRUE)),
  sprintf("max rows in one answer: %d", max(c(0L, ok$rows), na.rm = TRUE)),
  sprintf("releases seen: %s", paste(unique(stats::na.omit(requests$release)), collapse = ";"))
)
writeLines(summary_lines, file.path(out_dir, "f3-summary.txt"))
writeLines(summary_lines)
