# KIT-171 script S1: the MaxQuant reader on the real PXD000279 proteinGroups.txt (Cox et al. 2014, MaxQuant 1.3.10.9).
#
# Usage, from the package root on a run-tests runner, after the commit under test is installed:
#   Rscript docs/qa/workflow-execution/kit/kit-171/real-file-counts.R <out dir>
#
# Fetches proteomebenchmark.zip from PRIDE, checks its size, extracts proteinGroups.txt and refuses it unless its
# SHA-256 is the one in KIT-171 section 4. Runs MultiScholaR::importMaxQuantData with filter_contaminants TRUE and
# FALSE, and counts the same file read as text with base R. The input's licence class is unstated (VAL-7): only
# aggregate integers, the run names, hashes of the inputs and the commit are written to <out dir>/real-file-counts.tsv
# and printed; the reader's console and its log go to private files under RUNNER_TEMP, which are never uploaded; an
# error is reported by its class only, never by its message (a parse error can quote a field of the file). The data
# files are deleted before the script ends.

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[[1L]] else "k171-s1-out"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ZIP_URL <- "https://ftp.pride.ebi.ac.uk/pride/data/archive/2014/09/PXD000279/proteomebenchmark.zip"
ZIP_BYTES <- 38428469
MEMBER_BYTES <- 27782803
MEMBER_SHA256 <- "d245d4dc5b7ccc1a75eb21d36d89ed3fd1ace01ca858bd8eef9d65e3cc1411ea"
FLAG_COLUMNS <- c("Potential.contaminant", "Contaminant", "Reverse")

private_root <- Sys.getenv("RUNNER_TEMP")
if (!nzchar(private_root)) private_root <- tempdir()
work <- file.path(private_root, "kit171-s1")
dir.create(work, showWarnings = FALSE, recursive = TRUE)
zip_path <- file.path(work, "proteomebenchmark.zip")
member_path <- file.path(work, "proteinGroups.txt")

counts <- list()
record <- function(measure, value) counts[[measure]] <<- paste(as.character(value), collapse = ",")

writeCounts <- function() {
    table <- data.frame(measure = names(counts), value = unlist(counts, use.names = FALSE))
    utils::write.table(table, file.path(out_dir, "real-file-counts.tsv"), sep = "\t", quote = FALSE,
        row.names = FALSE)
    print(table, row.names = FALSE)
}

removeData <- function() unlink(c(zip_path, member_path))

fail <- function(stage, condition_class) {
    record("error_stage", stage)
    record("error_class", condition_class)
    removeData()
    writeCounts()
    quit(save = "no", status = 1L)
}

sha256 <- function(path) sub("[[:space:]].*$", "", system2("sha256sum", shQuote(path), stdout = TRUE)[[1L]])

record("commit", tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)[[1L]], error = function(e) NA))
record("installed_version", as.character(utils::packageVersion("MultiScholaR")))

# --- The input -----------------------------------------------------------------------------------------------------
options(timeout = 1800)
fetched <- tryCatch({
    utils::download.file(ZIP_URL, zip_path, mode = "wb", quiet = TRUE)
    TRUE
}, error = function(e) class(e)[[1L]], warning = function(w) class(w)[[1L]])
if (!isTRUE(fetched)) fail("download", fetched)
record("zip_bytes", file.size(zip_path))
if (!identical(as.numeric(file.size(zip_path)), ZIP_BYTES)) fail("zip_size", "size_mismatch")
record("zip_sha256", sha256(zip_path))

listing <- tryCatch(utils::unzip(zip_path, list = TRUE), error = function(e) class(e)[[1L]])
if (!is.data.frame(listing)) fail("zip_list", listing)
candidates <- listing$Name[basename(listing$Name) == "proteinGroups.txt"]
record("zip_members_named_proteinGroups", length(candidates))
if (length(candidates) != 1L) fail("zip_member", "member_count")
extracted <- tryCatch(utils::unzip(zip_path, files = candidates, exdir = work, junkpaths = TRUE),
    error = function(e) class(e)[[1L]], warning = function(w) class(w)[[1L]])
unlink(zip_path)
if (!file.exists(member_path)) fail("zip_extract", if (is.character(extracted)) extracted[[1L]] else "missing")
record("member_bytes", file.size(member_path))
member_hash <- sha256(member_path)
record("member_sha256", member_hash)
if (!identical(member_hash, MEMBER_SHA256)) fail("member_sha256", "hash_mismatch")

# --- The reader, its console and log kept private ------------------------------------------------------------------
callReader <- function(filter_contaminants, log_path) {
    console <- file(file.path(work, sprintf("reader-console-%s.txt", filter_contaminants)), open = "wt")
    sink(console)
    sink(console, type = "message")
    on.exit({
        sink(type = "message")
        sink()
        close(console)
    }, add = TRUE)
    logger::log_appender(function(lines) cat(lines, file = log_path, sep = "\n", append = TRUE))
    withCallingHandlers(
        tryCatch(
            MultiScholaR::importMaxQuantData(member_path, use_lfq = TRUE, filter_contaminants = filter_contaminants),
            error = function(e) structure(class(e)[[1L]], class = "k171_reader_error")
        ),
        warning = function(w) {
            message("warning: ", conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
}

for (filter in c(TRUE, FALSE)) {
    log_path <- file.path(work, sprintf("reader-log-%s.txt", filter))
    unlink(log_path)
    result <- callReader(filter, log_path)
    tag <- if (filter) "reader_filter_true" else "reader_filter_false"
    if (inherits(result, "k171_reader_error")) fail(tag, unclass(result))
    record(paste0(tag, "_distinct_protein_groups"), length(unique(result$data$Protein.Ids)))
    record(paste0(tag, "_missing_protein_ids"), sum(is.na(unique(result$data$Protein.Ids))))
    record(paste0(tag, "_runs"), sort(unique(result$data$Run)))
    log_lines <- if (file.exists(log_path)) readLines(log_path, warn = FALSE) else character()
    for (column in FLAG_COLUMNS) {
        pattern <- sprintf("MaxQuant import: removed ([0-9]+) rows marked '\\+' in %s$", gsub(".", "\\.", column,
            fixed = TRUE))
        hit <- regmatches(log_lines, regexec(pattern, log_lines))
        hit <- vapply(hit[lengths(hit) == 2L], `[[`, "", 2L)
        record(paste0(tag, "_logged_rows_removed_", column), if (length(hit) == 0L) "absent" else hit)
    }
    record(paste0(tag, "_logged_no_contaminant_column_lines"),
        sum(grepl("contaminant filtering requested but no contaminant column found", log_lines, fixed = TRUE)))
    rm(result)
}

# --- The same file read as text (no feature value is written) ------------------------------------------------------
cr_bytes <- sum(readBin(member_path, "raw", n = file.size(member_path)) == as.raw(13L))
record("text_cr_bytes", cr_bytes)
raw_lines <- readLines(member_path, warn = FALSE)
record("text_lines_including_header", length(raw_lines))
record("text_double_quote_characters",
    sum(lengths(regmatches(raw_lines, gregexpr("\"", raw_lines, fixed = TRUE)))))
field_counts <- lengths(regmatches(raw_lines, gregexpr("\t", raw_lines, fixed = TRUE))) + 1L
record("text_header_fields", field_counts[[1L]])
record("text_data_lines_fields_not_112", sum(field_counts[-1L] != 112L))

split <- strsplit(raw_lines, "\t", fixed = TRUE)
header <- split[[1L]]
rows <- split[-1L]
rm(raw_lines, split)
column <- function(name) {
    at <- match(name, header)
    if (is.na(at)) return(NULL)
    vapply(rows, function(x) if (length(x) >= at) x[[at]] else "", "")
}
for (name in c("Protein IDs", "Majority protein IDs", "Fasta headers", "Reverse", "Contaminant",
    "Potential contaminant", "Only identified by site")) {
    record(paste0("text_column_present_", gsub(" ", "_", name)), !is.na(match(name, header)))
}
majority <- column("Majority protein IDs")
protein_ids <- column("Protein IDs")
fasta <- column("Fasta headers")
reverse <- column("Reverse")
contaminant <- column("Contaminant")
by_site <- column("Only identified by site")
if (is.null(majority) || is.null(protein_ids) || is.null(fasta) || is.null(reverse) || is.null(contaminant) ||
    is.null(by_site)) {
    fail("text_columns", "column_missing")
}

marked_reverse <- reverse == "+"
marked_contaminant <- contaminant == "+"
record("text_data_rows", length(rows))
record("text_reverse_distinct_values", length(unique(reverse)))
record("text_contaminant_distinct_values", length(unique(contaminant)))
record("text_only_identified_by_site_distinct_values", length(unique(by_site)))
record("text_reverse_marked", sum(marked_reverse))
record("text_contaminant_marked", sum(marked_contaminant))
record("text_only_identified_by_site_marked", sum(by_site == "+"))
record("text_marked_reverse_and_contaminant", sum(marked_reverse & marked_contaminant))
record("text_majority_ids_empty", sum(!nzchar(majority)))
record("text_majority_ids_distinct_all_rows", length(unique(majority)))
record("text_majority_ids_duplicated_rows", length(majority) - length(unique(majority)))
record("text_majority_ids_distinct_unmarked_reverse", length(unique(majority[!marked_reverse])))
record("text_majority_ids_distinct_unmarked_reverse_and_contaminant",
    length(unique(majority[!marked_reverse & !marked_contaminant])))

majority_members <- strsplit(majority, ";", fixed = TRUE)
protein_members <- strsplit(protein_ids, ";", fixed = TRUE)
hasTag <- function(members, tag) vapply(members, function(x) any(startsWith(x, tag)), TRUE)
leadsWithTag <- function(members, tag) vapply(members, function(x) length(x) > 0L && startsWith(x[[1L]], tag), TRUE)
allTag <- function(members, tag) vapply(members, function(x) length(x) > 0L && all(startsWith(x, tag)), TRUE)

for (tag in c("CON__", "REV__")) {
    marked <- if (tag == "CON__") marked_contaminant else marked_reverse
    in_majority <- hasTag(majority_members, tag)
    in_protein_ids_only <- hasTag(protein_members, tag) & !in_majority
    leads <- leadsWithTag(majority_members, tag)
    only_tagged <- allTag(majority_members, tag)
    prefix <- paste0("text_", sub("__$", "", tag), "_")
    for (state in c("marked", "unmarked")) {
        rows_in <- if (state == "marked") marked else !marked
        record(paste0(prefix, state, "_member_in_majority"), sum(rows_in & in_majority))
        record(paste0(prefix, state, "_member_in_protein_ids_only"), sum(rows_in & in_protein_ids_only))
        record(paste0(prefix, state, "_no_member"), sum(rows_in & !in_majority & !in_protein_ids_only))
        record(paste0(prefix, state, "_majority_leads_with_member"), sum(rows_in & leads))
        record(paste0(prefix, state, "_majority_only_members"), sum(rows_in & only_tagged))
    }
}
record("text_contaminant_marked_fasta_headers_with_ECOLI", sum(marked_contaminant & grepl("_ECOLI", fasta,
    fixed = TRUE)))

removeData()
writeCounts()
