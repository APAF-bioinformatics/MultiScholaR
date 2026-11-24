# KIT-175 measurement F4 (runner only, network open; run before the owner's decision D8, committing to neither
# option): the bare-term query that directUniprotDownload sends today against the field-qualified form
# `(accession:A OR accession:B ...) AND organism_id:<taxon>`.
#   (a) F3's 120 requests (the first 2,000 reviewed human accessions of UP000005640 and the first 1,000 of
#       UP000000625 by accession, 25 per request) under both forms: Entry sets, rows, X-Total-Results, next link.
#   (b) Absent kinds, five probes per request beside 20 current human primary accessions, under organism 9606, both
#       forms: up to 200 current secondary accessions whose primary is a reviewed human entry (UniProt's sec_ac.txt),
#       up to 200 deleted Swiss-Prot accessions (UniProt's delac_sp.txt), and 100 E. coli K-12 accessions (another
#       organism than the one queried).
# Writes counts, header facts and the release only; no accession is written or printed (licence guard ENG-07).
#
# Usage (from the package root): Rscript docs/qa/workflow-execution/kit/kit-175/f4-query-forms.R <out dir>

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[[1L]] else "f4-out"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
options(timeout = 600)

fields <- paste0(
  "accession,id,protein_name,gene_names,organism_name,length,go_id,reviewed,",
  "protein_existence,annotation_score"
)
accession_pattern <- "^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
docs_url <- "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/docs/"

header_or_na <- function(response, name) {
  value <- httr::headers(response)[[name]]
  if (is.null(value)) NA_character_ else value
}

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
  table <- utils::read.delim(text = httr::content(response, "text", encoding = "UTF-8"), sep = "\t",
    quote = "", stringsAsFactors = FALSE)
  sort(table$Entry)
}

build_query <- function(accessions, taxon_id, form) {
  terms <- if (identical(form, "qualified")) paste0("accession:", accessions) else accessions
  paste0("(", paste(terms, collapse = " OR "), ") AND organism_id:", taxon_id)
}

# One request as directUniprotDownload sends it (httr::GET, the same fields, no size, 30 s, one attempt), 1 s apart.
search <- function(accessions, taxon_id, form) {
  response <- tryCatch(
    httr::GET(
      url = "https://rest.uniprot.org/uniprotkb/search",
      query = list(query = build_query(accessions, taxon_id, form), format = "tsv", fields = fields),
      httr::timeout(30)
    ),
    error = function(e) e
  )
  Sys.sleep(1)
  if (inherits(response, "error")) {
    return(list(status = NA_integer_, entries = character(), total = NA_integer_, link_next = NA,
      release = NA_character_))
  }
  status <- httr::status_code(response)
  entries <- character()
  if (status == 200L) {
    table <- tryCatch(
      suppressWarnings(utils::read.delim(text = httr::content(response, "text", encoding = "UTF-8"), sep = "\t",
        quote = "", stringsAsFactors = FALSE)),
      error = function(e) data.frame()
    )
    if ("Entry" %in% names(table)) entries <- table$Entry
  }
  link <- header_or_na(response, "link")
  list(status = status, entries = entries,
    total = suppressWarnings(as.integer(header_or_na(response, "x-total-results"))),
    link_next = !is.na(link) && grepl('rel="next"', link, fixed = TRUE),
    release = header_or_na(response, "x-uniprot-release"))
}

human <- fetch_list("UP000005640", reviewed_only = TRUE)
ecoli <- fetch_list("UP000000625", reviewed_only = FALSE)
releases <- character()

# ---- (a) F3's 120 requests under both forms ----
compare_rows <- list()
form_rows <- list()
run_a <- function(dataset, accessions, taxon_id) {
  chunks <- split(accessions, ceiling(seq_along(accessions) / 25))
  for (i in seq_along(chunks)) {
    bare <- search(chunks[[i]], taxon_id, "bare")
    qualified <- search(chunks[[i]], taxon_id, "qualified")
    releases <<- c(releases, bare$release, qualified$release)
    for (form in c("bare", "qualified")) {
      r <- if (form == "bare") bare else qualified
      form_rows[[length(form_rows) + 1L]] <<- data.frame(dataset = dataset, request = i, form = form,
        status = r$status, rows = length(r$entries), x_total_results = r$total, link_next = r$link_next,
        requested_returned = sum(chunks[[i]] %in% r$entries), entries_outside_request = sum(!r$entries %in% chunks[[i]]),
        stringsAsFactors = FALSE)
    }
    compare_rows[[length(compare_rows) + 1L]] <<- data.frame(dataset = dataset, request = i,
      both_status_200 = identical(bare$status, 200L) && identical(qualified$status, 200L),
      entry_sets_identical = setequal(bare$entries, qualified$entries), stringsAsFactors = FALSE)
  }
}
run_a("human_UP000005640_reviewed", utils::head(human, 2000L), 9606)
run_a("ecoli_UP000000625", utils::head(ecoli, 1000L), 83333)
form_a <- do.call(rbind, form_rows)
compare_a <- do.call(rbind, compare_rows)
utils::write.table(form_a, file.path(out_dir, "f4a-forms.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(compare_a, file.path(out_dir, "f4a-compare.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

summary_lines <- c(
  sprintf("(a) requests: %d; both forms status 200: %d; Entry sets identical: %d", nrow(compare_a),
    sum(compare_a$both_status_200), sum(compare_a$entry_sets_identical)),
  vapply(c("bare", "qualified"), function(form) {
    x <- form_a[form_a$form == form, ]
    sprintf("(a) %s: next link %d of %d; X-Total-Results above rows %d; requested returned %d of %d; outside %d",
      form, sum(x$link_next, na.rm = TRUE), nrow(x), sum(x$x_total_results > x$rows, na.rm = TRUE),
      sum(x$requested_returned), 3000L, sum(x$entries_outside_request))
  }, character(1))
)
writeLines(summary_lines, file.path(out_dir, "f4-summary.txt"))

# ---- (b) absent kinds ----
read_doc <- function(name) {
  path <- tempfile(fileext = ".txt")
  utils::download.file(paste0(docs_url, name), path, quiet = TRUE, mode = "wb")
  readLines(path, warn = FALSE)
}
probe_rows <- list()
probe_status <- character()
tryCatch({
  sec_lines <- read_doc("sec_ac.txt")
  pairs <- regmatches(sec_lines, regexec("^([A-Z0-9]{6,10})\\s+([A-Z0-9]{6,10})\\s*$", sec_lines))
  pairs <- do.call(rbind, lapply(Filter(function(m) length(m) == 3L, pairs), function(m) m[2:3]))
  pairs <- pairs[grepl(accession_pattern, pairs[, 1]) & grepl(accession_pattern, pairs[, 2]), , drop = FALSE]
  probe_status <- c(probe_status, sprintf("sec_ac.txt: %d secondary-primary pairs read", nrow(pairs)))

  del_lines <- trimws(read_doc("delac_sp.txt"))
  deleted <- sort(unique(del_lines[grepl(accession_pattern, del_lines)]))
  probe_status <- c(probe_status, sprintf("delac_sp.txt: %d deleted accessions read", length(deleted)))

  # 20 current human primaries per probe request, from the reviewed list beyond F3's first 2,000.
  fillers_pool <- human[-seq_len(2000L)]
  secondary_human <- pairs[pairs[, 2] %in% human & !pairs[, 2] %in% fillers_pool[seq_len(20L * 40L)], , drop = FALSE]
  secondary_ids <- sort(unique(secondary_human[, 1]))
  # A secondary that is also a current primary of a reviewed human entry is left out (it is not absent).
  secondary_ids <- setdiff(secondary_ids, human)
  set.seed(175L)
  secondary_sample <- if (length(secondary_ids) > 200L) sort(sample(secondary_ids, 200L)) else secondary_ids
  deleted_sample <- if (length(deleted) > 200L) sort(sample(deleted, 200L)) else deleted
  other_sample <- utils::head(ecoli, 100L)
  probe_status <- c(probe_status, sprintf("samples: secondary %d, deleted %d, other organism %d",
    length(secondary_sample), length(deleted_sample), length(other_sample)))

  run_b <- function(kind, probes, filler_offset) {
    chunks <- split(probes, ceiling(seq_along(probes) / 5))
    for (i in seq_along(chunks)) {
      fillers <- fillers_pool[filler_offset + (i - 1L) * 20L + seq_len(20L)]
      request <- c(fillers, chunks[[i]])
      mapped <- if (identical(kind, "secondary")) secondary_human[secondary_human[, 1] %in% chunks[[i]], , drop = FALSE] else NULL
      for (form in c("bare", "qualified")) {
        r <- search(request, 9606, form)
        releases <<- c(releases, r$release)
        as_requested <- sum(chunks[[i]] %in% r$entries)
        under_other <- if (is.null(mapped)) 0L else sum(vapply(chunks[[i]], function(p)
          !(p %in% r$entries) && any(mapped[mapped[, 1] == p, 2] %in% r$entries), logical(1)))
        known <- c(request, if (!is.null(mapped)) mapped[, 2])
        probe_rows[[length(probe_rows) + 1L]] <<- data.frame(kind = kind, request = i, form = form,
          status = r$status, probes = length(chunks[[i]]), probes_found_as_requested = as_requested,
          probes_found_under_other_accession = under_other,
          probes_not_found = length(chunks[[i]]) - as_requested - under_other,
          fillers_returned = sum(fillers %in% r$entries), rows = length(r$entries),
          entries_outside_known = sum(!r$entries %in% known), x_total_results = r$total, link_next = r$link_next,
          stringsAsFactors = FALSE)
      }
    }
  }
  run_b("secondary", secondary_sample, 0L)
  run_b("deleted", deleted_sample, 0L)
  run_b("other_organism_ecoli_under_9606", other_sample, 0L)
}, error = function(e) {
  probe_status <<- c(probe_status, paste("(b) stopped:", conditionMessage(e)))
})

if (length(probe_rows)) {
  probes <- do.call(rbind, probe_rows)
  utils::write.table(probes, file.path(out_dir, "f4b-probes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  for (kind in unique(probes$kind)) for (form in c("bare", "qualified")) {
    x <- probes[probes$kind == kind & probes$form == form, ]
    probe_status <- c(probe_status, sprintf(paste("(b) %s, %s: requests %d; probes %d: as requested %d, under another",
      "accession %d, not found %d; fillers returned %d of %d; next link on %d requests (on %d of the %d requests with",
      "a probe not found as requested)"), kind, form, nrow(x), sum(x$probes), sum(x$probes_found_as_requested),
      sum(x$probes_found_under_other_accession), sum(x$probes_not_found), sum(x$fillers_returned), 20L * nrow(x),
      sum(x$link_next, na.rm = TRUE), sum(x$link_next & x$probes_found_as_requested < x$probes, na.rm = TRUE),
      sum(x$probes_found_as_requested < x$probes)))
  }
}
summary_lines <- c(summary_lines, probe_status,
  sprintf("releases seen: %s", paste(unique(stats::na.omit(releases)), collapse = ";")))
writeLines(summary_lines, file.path(out_dir, "f4-summary.txt"))
writeLines(summary_lines)
