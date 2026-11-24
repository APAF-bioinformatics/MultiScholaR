# KIT-167 gate F2 (section 11, "Feasibility gates"): size and time of one mapping at the shipped settings.
#
# KPC2 (PRIDE PXD016672, BPA_SEP_KPC2_Contaminants.fasta, CC0) against the shipped K. pneumoniae MGH 78578 reference
# (UP000000265, UniProt 2026_04, CC BY 4.0, inst/extdata/reference_mapping/references/), DIAMOND v2.2.8 (GPL-3, a
# program; GitHub release asset verified by bytes and SHA-256 against section 2.2), threads 4, both directions, at the
# settings of section 5.1. Measured: wall seconds and peak RSS per direction, hit rows, whether the raw `score` field
# is an integer on every row (decides the ranking key, 5.2), the executable's SHA-256 (for the pin file, 8.6), the
# processor model and the `avx512f` flag. Stop rule (F2): over 10 minutes or 8 GB for one direction at the shipped
# settings -> returned to the owner with the numbers.
#
# Writes only aggregates to capture-out/f2/ (no hit table, no sequence, no identifier). Runs on a runner only.

out_dir <- file.path(getwd(), "capture-out", "f2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
scratch <- file.path(tempdir(), "k167-f2")
dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
reference_fasta <- normalizePath(file.path("inst", "extdata", "reference_mapping", "references", "UP000000265.fasta"))

rows <- list()
put <- function(field, value) {
    value <- paste(as.character(value), collapse = ";")
    rows[[length(rows) + 1L]] <<- data.frame(field = field, value = value, stringsAsFactors = FALSE)
    cat(sprintf("KIT167 %s %s\n", field, value))
}
finish <- function(status) {
    put("status", status)
    table <- do.call(rbind, rows)
    utils::write.table(table, file.path(out_dir, "f2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    if (!identical(status, "ok")) quit(status = 1L, save = "no")
}
sha256 <- function(path) digest::digest(file = path, algo = "sha256")

fetch <- function(url, dest) {
    for (attempt in 1:3) {
        ok <- tryCatch({
            utils::download.file(url, dest, mode = "wb", quiet = TRUE, method = "libcurl")
            TRUE
        }, error = function(e) {
            message("download attempt ", attempt, " failed: ", conditionMessage(e))
            FALSE
        })
        if (ok && file.exists(dest)) return(TRUE)
        Sys.sleep(30)
    }
    FALSE
}
options(timeout = 300)

cpu <- tryCatch(readLines("/proc/cpuinfo"), error = function(e) character())
put("processor", sub("^model name\\s*:\\s*", "", grep("^model name", cpu, value = TRUE)[1]))
put("avx512f", any(grepl("\\bavx512f\\b", grep("^flags", cpu, value = TRUE))))
put("cores", parallel::detectCores())

# The engine: v2.2.8 Linux asset, verified against section 2.2.
archive <- file.path(scratch, "diamond-linux64.tar.gz")
if (!fetch("https://github.com/bbuchfink/diamond/releases/download/v2.2.8/diamond-linux64.tar.gz", archive)) {
    put("failed_check", "download_engine")
    finish("failed")
}
put("engine_archive_bytes", file.size(archive))
put("engine_archive_sha256", sha256(archive))
if (file.size(archive) != 10221459 ||
    sha256(archive) != "33a765d192ae4f5b388f15fe2e7c23d3aaea0f9e3f6872efe8c324aec261fbc7") {
    put("failed_check", "engine_archive_mismatch")
    finish("failed")
}
utils::untar(archive, exdir = scratch)
diamond <- list.files(scratch, pattern = "^diamond$", recursive = TRUE, full.names = TRUE)[1]
Sys.chmod(diamond, "0755")
put("engine_archive_members", paste(utils::untar(archive, list = TRUE), collapse = ","))
put("engine_executable_sha256", sha256(diamond))
put("engine_version", system2(diamond, "version", stdout = TRUE))
help_text <- tryCatch(system2(diamond, "help", stdout = TRUE, stderr = TRUE), error = function(e) character())
cat(help_text, sep = "\n")

# The strain FASTA: found in the project's PRIDE FTP folder, dated by publication (year/month). The folder is taken
# from the PRIDE API's project record when it answers, else every month of 2019 to 2024 is listed; the file is the
# folder listing's entry matching the name, case-insensitive. Each step prints what it read (instrument run 2, after
# run 45334031153 found no URL).
readText <- function(url) {
    text <- tryCatch(paste(readLines(url, warn = FALSE), collapse = "\n"), error = function(e) {
        message("read failed: ", url, ": ", conditionMessage(e))
        ""
    })
    text
}
pride_url <- function(accession, name) {
    folders <- character()
    for (endpoint in c(sprintf("https://www.ebi.ac.uk/pride/ws/archive/v3/projects/%s", accession),
        sprintf("https://www.ebi.ac.uk/pride/ws/archive/v2/projects/%s", accession))) {
        text <- readText(endpoint)
        cat(sprintf("PRIDE %s: %d characters\n", endpoint, nchar(text)))
        date <- regmatches(text, regexpr("\"publicationDate\"[[:space:]]*:[[:space:]]*\"[0-9]{4}-[0-9]{2}", text))
        if (length(date)) {
            ym <- regmatches(date, regexpr("[0-9]{4}-[0-9]{2}$", date))
            folders <- c(folders, sprintf("https://ftp.pride.ebi.ac.uk/pride/data/archive/%s/%s/%s/",
                substr(ym, 1, 4), substr(ym, 6, 7), accession))
        }
    }
    folders <- unique(c(folders, as.vector(outer(sprintf("%02d/%s/", 1:12, accession), 2019:2024,
        function(tail, year) sprintf("https://ftp.pride.ebi.ac.uk/pride/data/archive/%d/%s", year, tail)))))
    for (folder in folders) {
        listing <- suppressWarnings(readText(folder))
        if (!nzchar(listing)) next
        files <- unique(regmatches(listing, gregexpr("href=\"[^\"]+\"", listing))[[1]])
        files <- sub("^href=\"(.*)\"$", "\\1", files)
        cat(sprintf("PRIDE folder %s: %d entries\n", folder, length(files)))
        hit <- files[tolower(basename(files)) == tolower(name)]
        cat("PRIDE folder FASTA entries:", paste(grep("\\.fasta", basename(files), value = TRUE,
            ignore.case = TRUE), collapse = " "), "\n")
        if (length(hit)) return(paste0(folder, basename(hit[[1]])))
    }
    NA_character_
}
strain_name <- "BPA_SEP_KPC2_Contaminants.fasta"
# Instrument run 3: run 45335506224 found PXD016672's folder (2021/07) without the file; the research note places the
# search FASTA in the DDA project of the same strains (PXD020839, K. pneumoniae species complex), tried second.
strain_url <- pride_url("PXD016672", strain_name)
if (is.na(strain_url)) strain_url <- pride_url("PXD020839", strain_name)
put("strain_url", strain_url)
strain_raw <- file.path(scratch, "strain-raw.fasta")
if (is.na(strain_url) || !fetch(strain_url, strain_raw)) {
    put("failed_check", "download_strain")
    finish("failed")
}
put("strain_bytes", file.size(strain_raw))
put("strain_sha256", sha256(strain_raw))
if (file.size(strain_raw) != 15817761) {
    put("failed_check", "strain_size_differs_from_2.2")
    finish("failed")
}

# Strain entries for F2: the `>ref|...|` class (section 8.4); every other class counted, not selected.
read_fasta <- function(path) {
    lines <- readLines(path, warn = FALSE)
    starts <- which(startsWith(lines, ">"))
    ends <- c(starts[-1L] - 1L, length(lines))
    data.frame(header = substring(lines[starts], 2L),
        sequence = vapply(seq_along(starts), function(i)
            if (ends[i] > starts[i]) paste(lines[(starts[i] + 1L):ends[i]], collapse = "") else "", character(1)),
        stringsAsFactors = FALSE)
}
strain <- read_fasta(strain_raw)
is_ref <- startsWith(strain$header, "ref|")
put("strain_entries", nrow(strain))
put("strain_ref_entries", sum(is_ref))
put("strain_uniprot_style_entries", sum(grepl("^(sp|tr)\\|", strain$header)))
put("strain_other_entries", sum(!is_ref & !grepl("^(sp|tr)\\|", strain$header)))
selected <- strain[is_ref, ]
selected$id <- sub("^ref\\|([^|]+)\\|.*$", "\\1", selected$header)
selected$sequence <- sub("\\*$", "", toupper(selected$sequence))
selected <- selected[!duplicated(selected$sequence), ]
put("strain_search_entities", nrow(selected))
writeLines(paste0(">", selected$id, "\n", selected$sequence), file.path(scratch, "strain.fasta"))
file.copy(reference_fasta, file.path(scratch, "reference.fasta"))
reference <- read_fasta(reference_fasta)
put("reference_entries", nrow(reference))
put("reference_sha256", sha256(reference_fasta))

# The search, section 5.1; relative names inside the scratch folder (v2.2.8 and spaces).
old <- setwd(scratch)
time_bin <- if (file.exists("/usr/bin/time")) "/usr/bin/time" else NA_character_
put("peak_rss_tool", if (is.na(time_bin)) "unavailable" else time_bin)
run_engine <- function(label, args) {
    log <- paste0(label, ".time")
    started <- Sys.time()
    status <- if (is.na(time_bin)) {
        system2(diamond, args, stdout = paste0(label, ".out"), stderr = paste0(label, ".err"))
    } else {
        system2(time_bin, c("-v", "-o", log, diamond, args), stdout = paste0(label, ".out"),
            stderr = paste0(label, ".err"))
    }
    seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    rss_kb <- NA_real_
    if (!is.na(time_bin) && file.exists(log)) {
        line <- grep("Maximum resident set size", readLines(log), value = TRUE)
        if (length(line)) rss_kb <- as.numeric(sub(".*:\\s*", "", line[1]))
    }
    put(paste0(label, "_exit"), status)
    put(paste0(label, "_seconds"), sprintf("%.1f", seconds))
    put(paste0(label, "_peak_rss_mb"), if (is.na(rss_kb)) "NA" else sprintf("%.1f", rss_kb / 1024))
    put(paste0(label, "_args"), paste(args, collapse = " "))
    if (!identical(as.integer(status), 0L)) {
        cat(readLines(paste0(label, ".err")), sep = "\n")
        setwd(old)
        put("failed_check", paste0(label, "_exit_status"))
        finish("failed")
    }
    seconds
}
run_engine("makedb_reference", c("makedb", "--in", "reference.fasta", "--db", "reference", "--threads", "4"))
run_engine("makedb_strain", c("makedb", "--in", "strain.fasta", "--db", "strain", "--threads", "4"))
dir.create("tmp", showWarnings = FALSE)
fields <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
    "evalue", "bitscore", "score", "qlen", "slen", "qcovhsp", "scovhsp")
search <- function(query, db, out) c("blastp", "--query", query, "--db", db, "--out", out, "--ultra-sensitive",
    "--masking", "none", "--ext", "full", "--comp-based-stats", "1", "-k0", "--max-hsps", "1", "--evalue", "1e-5",
    "--block-size", "0.4", "--index-chunks", "1", "--algo", "0", "--threads", "4", "--tmpdir", "tmp",
    "--outfmt", "6", fields)
forward <- run_engine("forward", search("strain.fasta", "reference", "forward.tsv"))
reverse <- run_engine("reverse", search("reference.fasta", "strain", "reverse.tsv"))
for (label in c("forward", "reverse")) {
    hits <- utils::read.delim(paste0(label, ".tsv"), header = FALSE, col.names = fields, colClasses = "character",
        quote = "", comment.char = "")
    put(paste0(label, "_rows"), nrow(hits))
    put(paste0(label, "_queries_with_hits"), length(unique(hits$qseqid)))
    put(paste0(label, "_score_integer"), all(grepl("^-?[0-9]+$", hits$score)))
    put(paste0(label, "_bitscore_examples"), paste(utils::head(unique(hits$bitscore), 5), collapse = ","))
    put(paste0(label, "_score_examples"), paste(utils::head(unique(hits$score), 5), collapse = ","))
    put(paste0(label, "_output_mb"), sprintf("%.2f", file.size(paste0(label, ".tsv")) / 2^20))
}
setwd(old)
put("over_10_min", max(forward, reverse) > 600)
put("pins_used", "--block-size 0.4 --index-chunks 1 --algo 0")
finish("ok")
