# The input manifest, the record phase, retrieval and phase pre (KIT-159 section 5.1).

MANIFEST_COLUMNS <- c(
    "dataset", "role", "file", "url", "accession", "licence", "kind", "redistribution", "citation",
    "sha256", "bytes", "members", "keep_archive", "optional", "commit", "retrieved", "source_note"
)
DRAFT_COLUMNS <- c("dataset", "role", "url", "members", "licence", "kind", "redistribution",
                   "citation", "source_note", "truth")
USER_AGENT <- "data-fetch/1.0"

# The one registry of verified inputs of this process (valInputPath reads it).
.val_inputs <- new.env(parent = emptyenv())

#' SHA-256 of a file, streamed
#' @noRd
valSha256 <- function(path) digest::digest(file = path, algo = "sha256")

#' Validate a manifest table (section 5.1, retrieval step 1)
#' @noRd
valValidateManifest <- function(m) {
    missing <- setdiff(MANIFEST_COLUMNS, names(m))
    bad <- function(what) valAbort("manifest: {what}", "val_manifest_invalid")
    if (length(missing) > 0L) bad(paste("missing column", paste(missing, collapse = ", ")))
    if (!all(grepl("^(https://|repo://)", m$url))) bad("a URL is neither https:// nor repo://")
    if (!all(grepl("^[0-9a-f]{64}$", m$sha256))) bad("a hash is not 64 hexadecimal characters")
    if (!all(grepl("^[0-9]+$", m$bytes)) || any(as.numeric(m$bytes) <= 0)) {
        bad("bytes is not a positive integer")
    }
    if (anyDuplicated(m$role)) bad("a duplicate role")
    if (!all(m$kind %in% c("data", "program"))) bad("kind is not data or program")
    if (!all(m$redistribution %in% c("permitted", "not_permitted"))) {
        bad("redistribution is not permitted or not_permitted")
    }
    if (any(!nzchar(trimws(m$citation)))) bad("an input row without a citation")
    for (col in c("keep_archive", "optional", "commit")) {
        m[[col]] <- toupper(as.character(m[[col]])) %in% "TRUE"
    }
    if (any(m$commit & !(m$kind == "data" & m$redistribution == "permitted"))) {
        bad("commit = TRUE for an input that is not permitted data")
    }
    m$bytes <- as.numeric(m$bytes)
    m$members[is.na(m$members)] <- ""
    m
}

#' Read and validate an input manifest (`inputs.tsv`)
#'
#' @param path Path of `inputs.tsv`.
#' @return The manifest as a data frame.
#' @export
valReadManifest <- function(path) {
    if (!file.exists(path)) valAbort("manifest {path} does not exist", "val_manifest_invalid")
    m <- utils::read.delim(path, colClasses = "character", na.strings = character(),
                           quote = "", stringsAsFactors = FALSE)
    valValidateManifest(m)
}

#' Free bytes of the file system holding `dir` (its nearest existing parent)
#' @noRd
valFreeBytes <- function(dir) {
    d <- dir
    while (!dir.exists(d)) d <- dirname(d)
    out <- system2("df", c("-B1", "--output=avail", shQuote(d)), stdout = TRUE)
    as.numeric(trimws(out[[length(out)]]))
}

#' Download with the harness's settings (3,600 s, user agent `data-fetch/1.0`, no credentials)
#' @noRd
valDownload <- function(url, destfile) {
    h <- curl::new_handle(useragent = USER_AGENT, timeout = 3600)
    curl::curl_download(url, destfile, handle = h, quiet = TRUE)
}

#' Obtain one file to `tmp`: a repo:// copy or an https download with retries; TRUE on success
#' @noRd
valObtain <- function(url, tmp, repo_root, download, attempts, backoff) {
    if (startsWith(url, "repo://")) {
        src <- file.path(repo_root, sub("^repo://", "", url))
        return(file.exists(src) && file.copy(src, tmp, overwrite = TRUE))
    }
    for (i in seq_len(attempts)) {
        ok <- tryCatch({
            download(url, tmp)
            file.exists(tmp)
        }, error = function(e) FALSE)
        if (ok) return(TRUE)
        if (i < attempts) Sys.sleep(backoff)
    }
    FALSE
}

#' Check a file's size, then its hash; on mismatch delete it and refuse
#' @noRd
valVerifyFile <- function(path, bytes, sha256, what) {
    if (file.size(path) != bytes) {
        unlink(path)
        valAbort("{what}: size differs from the manifest", "val_input_size_mismatch")
    }
    if (!identical(valSha256(path), sha256)) {
        unlink(path)
        valAbort("{what}: SHA-256 differs from the manifest", "val_input_hash_mismatch")
    }
    invisible(TRUE)
}

#' Read a members file (`members-<role>.tsv`: member, sha256, bytes)
#' @noRd
valReadMembers <- function(members_dir, name) {
    path <- file.path(members_dir, name)
    if (!file.exists(path)) valAbort("members file {name} does not exist", "val_manifest_invalid")
    mm <- utils::read.delim(path, colClasses = "character", stringsAsFactors = FALSE)
    if (!all(c("member", "sha256", "bytes") %in% names(mm)) ||
            !all(grepl("^[0-9a-f]{64}$", mm$sha256))) {
        valAbort("members file {name} is not valid", "val_manifest_invalid")
    }
    mm$bytes <- as.numeric(mm$bytes)
    mm
}

#' Extract and verify the declared members of an archive
#' @noRd
valExtractMembers <- function(archive, mm, exdir, role) {
    dir.create(exdir, recursive = TRUE, showWarnings = FALSE)
    if (grepl("\\.zip$", archive, ignore.case = TRUE)) {
        utils::unzip(archive, files = mm$member, exdir = exdir)
    } else {
        utils::untar(archive, files = mm$member, exdir = exdir, tar = "internal")
    }
    paths <- file.path(exdir, mm$member)
    for (i in seq_len(nrow(mm))) {
        if (!file.exists(paths[[i]])) {
            valAbort("{role}: member {mm$member[[i]]} is absent", "val_input_hash_mismatch")
        }
        valVerifyFile(paths[[i]], mm$bytes[[i]], mm$sha256[[i]],
                      paste0(role, ":", mm$member[[i]]))
    }
    paths
}

#' Fetch and verify the declared inputs (section 5.1, retrieval)
#'
#' Validates the manifest; refuses when 2.5 x the bytes to fetch (archives plus members) exceed the
#' free space of `dest`; downloads each file once (or copies a `repo://` file), checks size then
#' SHA-256, extracts and verifies every declared member. An optional role that cannot be fetched is
#' recorded `input_unavailable`. Verified paths are then reached only through [valInputPath()].
#'
#' @param manifest Path of `inputs.tsv`, or a manifest data frame.
#' @param dest Folder for the verified files (`$RUNNER_TEMP/val-data`).
#' @param roles Roles to fetch, or `NULL` for all.
#' @param repo_root Root that `repo://` paths are relative to.
#' @param members_dir Folder of the `members-*.tsv` files (default: the manifest's folder).
#' @param free_bytes Free space override (tests).
#' @param download Download function `(url, destfile)` (default: curl with the harness settings).
#' @param attempts,backoff Download attempts and the back-off in seconds.
#' @param truth_roles Roles only the ticket's truth code reads (refused in phase `pre`).
#' @return A data frame: `role`, `member`, `path`, `status`.
#' @export
valFetchInputs <- function(manifest, dest, roles = NULL, repo_root = ".", members_dir = NULL,
                           free_bytes = NULL, download = NULL, attempts = 3L, backoff = 30,
                           truth_roles = character()) {
    if (is.character(manifest)) {
        members_dir <- members_dir %||% dirname(manifest)
        m <- valReadManifest(manifest)
    } else {
        m <- valValidateManifest(manifest)
    }
    if (!is.null(roles)) m <- m[m$role %in% roles, , drop = FALSE]
    download <- download %||% valDownload
    members <- lapply(m$members, function(x) if (nzchar(x)) valReadMembers(members_dir, x))
    need <- sum(m$bytes) + sum(vapply(members, function(x) sum(x$bytes %||% 0), numeric(1)))
    free <- free_bytes %||% valFreeBytes(dest)
    if (2.5 * need > free) {
        valAbort("2.5 x {need} bytes to fetch exceed the {free} bytes free",
                 "val_disk_insufficient")
    }
    rows <- list()
    for (i in seq_len(nrow(m))) {
        r <- m[i, ]
        dir.create(file.path(dest, r$dataset), recursive = TRUE, showWarnings = FALSE)
        target <- file.path(dest, r$dataset, r$file)
        tmp <- paste0(target, ".part")
        if (!valObtain(r$url, tmp, repo_root, download, attempts, backoff)) {
            unlink(tmp)
            if (!r$optional) {
                valAbort("{r$role}: the input could not be fetched", "val_input_unavailable")
            }
            rows[[length(rows) + 1L]] <- data.frame(role = r$role, member = "",
                                                    path = NA_character_,
                                                    status = "input_unavailable")
            next
        }
        valVerifyFile(tmp, r$bytes, r$sha256, r$role)
        file.rename(tmp, target)
        rows[[length(rows) + 1L]] <- data.frame(role = r$role, member = "", path = target,
                                                status = "verified")
        if (!is.null(members[[i]])) {
            exdir <- file.path(dest, r$dataset, paste0(r$role, "-members"))
            paths <- valExtractMembers(target, members[[i]], exdir, r$role)
            rows[[length(rows) + 1L]] <- data.frame(role = r$role, member = members[[i]]$member,
                                                    path = paths, status = "verified")
        }
    }
    out <- do.call(rbind, rows)
    .val_inputs$paths <- out
    .val_inputs$truth_roles <- truth_roles
    out
}

#' The path of a verified input (the only way the harness reaches an input)
#'
#' @param role Declared role.
#' @param member Archive member, or `NULL` for the file itself.
#' @return The verified path.
#' @export
valInputPath <- function(role, member = NULL) {
    reg <- .val_inputs$paths
    if (identical(.val_inputs$phase, "pre") && role %in% .val_inputs$truth_roles) {
        valAbort("phase pre may not read the truth role {role}", "val_spec_invalid")
    }
    want <- member %||% ""
    hit <- if (is.null(reg)) integer() else which(reg$role == role & reg$member == want &
                                                      reg$status == "verified")
    if (length(hit) != 1L) valAbort("input {role} {want} is not declared and verified",
                                    "val_input_undeclared")
    reg$path[[hit]]
}

#' Record inputs from a draft (the record phase, section 5.1)
#'
#' Downloads each draft row, hashes it, lists the declared archive members (semicolon-separated
#' names in `members`) with their hashes, and writes `inputs.tsv` and `members-<role>.tsv` into
#' `out_dir` (hashes and sizes only).
#'
#' @param draft `draft-inputs.tsv` path or data frame.
#' @param dest Download folder (private).
#' @param out_dir Output folder.
#' @param repo_root Root of `repo://` paths.
#' @param download,backoff As [valFetchInputs()].
#' @return The recorded manifest.
#' @export
valRecordInputs <- function(draft, dest, out_dir, repo_root = ".", download = NULL, backoff = 30) {
    if (is.character(draft)) {
        draft <- utils::read.delim(draft, colClasses = "character", na.strings = character(),
                                   stringsAsFactors = FALSE)
    }
    missing <- setdiff(DRAFT_COLUMNS, names(draft))
    if (length(missing) > 0L) {
        valAbort("draft-inputs: missing column {paste(missing, collapse = ', ')}",
                 "val_manifest_invalid")
    }
    download <- download %||% valDownload
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    rows <- lapply(seq_len(nrow(draft)), function(i) {
        r <- draft[i, ]
        file <- basename(sub("[?#].*$", "", r$url))
        path <- file.path(dest, r$dataset, file)
        dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
        if (!valObtain(r$url, path, repo_root, download, 3L, backoff)) {
            valAbort("{r$role}: the input could not be fetched", "val_input_unavailable")
        }
        members_file <- ""
        if (nzchar(r$members)) {
            names_m <- trimws(strsplit(r$members, ";", fixed = TRUE)[[1L]])
            exdir <- file.path(dest, r$dataset, paste0(r$role, "-members"))
            mm <- data.frame(member = names_m, sha256 = NA_character_, bytes = NA_real_)
            valExtractMembersUnverified(path, names_m, exdir)
            mm$sha256 <- vapply(file.path(exdir, names_m), valSha256, character(1))
            mm$bytes <- file.size(file.path(exdir, names_m))
            members_file <- paste0("members-", r$role, ".tsv")
            utils::write.table(mm, file.path(out_dir, members_file), sep = "\t", quote = FALSE,
                               row.names = FALSE)
        }
        data.frame(dataset = r$dataset, role = r$role, file = file, url = r$url, accession = "",
                   licence = r$licence, kind = r$kind, redistribution = r$redistribution,
                   citation = r$citation, sha256 = valSha256(path),
                   bytes = format(file.size(path), scientific = FALSE), members = members_file,
                   keep_archive = FALSE, optional = FALSE, commit = FALSE,
                   retrieved = format(Sys.Date()), source_note = r$source_note,
                   stringsAsFactors = FALSE)
    })
    m <- do.call(rbind, rows)
    utils::write.table(m, file.path(out_dir, "inputs.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE)
    valValidateManifest(m)
}

#' Extract named members without verification (the record phase computes their hashes)
#' @noRd
valExtractMembersUnverified <- function(archive, names_m, exdir) {
    dir.create(exdir, recursive = TRUE, showWarnings = FALSE)
    if (grepl("\\.zip$", archive, ignore.case = TRUE)) {
        utils::unzip(archive, files = names_m, exdir = exdir)
    } else {
        utils::untar(archive, files = names_m, exdir = exdir, tar = "internal")
    }
    invisible(file.path(exdir, names_m))
}

#' Compare the twin records of one draft (ENG-13)
#'
#' @param dir_a,dir_b Record output folders.
#' @param listed_sizes Named sizes the source lists (role = bytes), or `NULL`.
#' @return A list: `identical` (hashes, sizes and members equal, and equal to listed sizes) and
#'   `rows` (per role).
#' @export
valCompareRecords <- function(dir_a, dir_b, listed_sizes = NULL) {
    read <- function(d) {
        utils::read.delim(file.path(d, "inputs.tsv"), colClasses = "character",
                          na.strings = character(), stringsAsFactors = FALSE)
    }
    a <- read(dir_a)
    b <- read(dir_b)
    m <- merge(a[, c("role", "sha256", "bytes", "members")],
               b[, c("role", "sha256", "bytes", "members")], by = "role", all = TRUE,
               suffixes = c("_a", "_b"))
    m$equal <- m$sha256_a %in% m$sha256_b & m$sha256_a == m$sha256_b & m$bytes_a == m$bytes_b
    m$equal[is.na(m$equal)] <- FALSE
    for (i in which(nzchar(m$members_a %||% ""))) {
        fa <- file.path(dir_a, m$members_a[[i]])
        fb <- file.path(dir_b, m$members_b[[i]])
        same <- file.exists(fa) && file.exists(fb) && identical(valSha256(fa), valSha256(fb))
        m$equal[[i]] <- m$equal[[i]] && same
    }
    if (!is.null(listed_sizes)) {
        lst <- listed_sizes[m$role]
        m$listed_equal <- is.na(lst) | as.character(lst) == m$bytes_a
        m$equal <- m$equal & m$listed_equal
    }
    list(identical = all(m$equal), rows = m)
}

#' The UniProt release from a response header (x-uniprot-release), or NA
#' @noRd
valUniprotRelease <- function() {
    url <- paste0("https://rest.uniprot.org/uniprotkb/search?",
                  "query=accession:P69905&fields=accession&size=1")
    res <- tryCatch(curl::curl_fetch_memory(url, handle = curl::new_handle(useragent = USER_AGENT)),
                    error = function(e) NULL)
    if (is.null(res)) return(NA_character_)
    h <- curl::parse_headers_list(res$headers)
    h[["x-uniprot-release"]] %||% NA_character_
}

#' Record the UniProt cache of a proteomics data set (section 5.1, G6, pool row T3)
#'
#' Runs the package's own `getUniprotAnnotationsFull` (network open) into `cache_dir`, catching the
#' classed conditions of KIT-175: `prot_annotation_uniprot_incomplete` refuses with
#' `val_annotation_incomplete`; `prot_annotation_uniprot_empty` or no rows with
#' `val_annotation_empty`. Records the UniProt release from the response header and the SHA-256 of
#' the cache table sorted by its first column.
#'
#' @param data_tbl The data set's protein table.
#' @param protein_id_column Its protein identifier column.
#' @param taxon_id The declared taxon (the query is restricted to it).
#' @param cache_dir Cache folder.
#' @param annotate Annotation function (default: the package's `getUniprotAnnotationsFull`).
#' @param release Function returning the UniProt release (default: one header request).
#' @return A list: `table`, `release`, `release_source`, `sha256`, `path`.
#' @export
valRecordUniprotCache <- function(data_tbl, protein_id_column, taxon_id, cache_dir,
                                  annotate = NULL, release = NULL) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    annotate <- annotate %||% function(data_tbl, protein_id_column, taxon_id, cache_dir) {
        fun <- getExportedValue("MultiScholaR", "getUniprotAnnotationsFull")
        fun(data_tbl = data_tbl, protein_id_column = protein_id_column, cache_dir = cache_dir,
            taxon_id = taxon_id)
    }
    seen <- character()
    tbl <- withCallingHandlers(
        annotate(data_tbl, protein_id_column, taxon_id, cache_dir),
        condition = function(cnd) {
            if (inherits(cnd, "prot_annotation_uniprot_incomplete")) seen <<- c(seen, "incomplete")
            if (inherits(cnd, "prot_annotation_uniprot_empty")) seen <<- c(seen, "empty")
        }
    )
    if ("incomplete" %in% seen) {
        valAbort("UniProt answered incompletely; the cache is not recorded",
                 "val_annotation_incomplete")
    }
    if ("empty" %in% seen || is.null(tbl) || nrow(tbl) == 0L) {
        valAbort("UniProt returned no annotation; the cache is not recorded",
                 "val_annotation_empty")
    }
    rel <- (release %||% valUniprotRelease)()
    sorted <- tbl[order(as.character(tbl[[1L]]), method = "radix"), , drop = FALSE]
    text <- utils::capture.output(utils::write.table(sorted, sep = "\t", quote = FALSE,
                                                     row.names = FALSE))
    path <- file.path(cache_dir, "uniprot_annotations.RDS")
    saveRDS(tbl, path)
    list(table = tbl, release = if (is.na(rel)) format(Sys.Date()) else rel,
         release_source = if (is.na(rel)) "retrieval date" else "x-uniprot-release",
         sha256 = digest::digest(paste(text, collapse = "\n"), algo = "sha256", serialize = FALSE),
         path = path)
}

#' Phase pre: run a ticket's pre-declaration script and admit only its listed typed tables
#'
#' The script is evaluated in an environment holding `out_dir`, with [valInputPath()] refusing
#' truth roles. Every file it writes must be a `.tsv` listed in `pre-outputs.tsv` (file, column,
#' type, labels), with exactly its columns and types, label values from the declared vocabulary and
#' at most 1,000 rows; anything else is `val_publish_forbidden`.
#'
#' @param script Path of the ticket's script.
#' @param pre_outputs Path of `pre-outputs.tsv`.
#' @param out_dir Output folder.
#' @return `TRUE` when every output is admitted.
#' @export
valRunPre <- function(script, pre_outputs, out_dir) {
    spec <- utils::read.delim(pre_outputs, colClasses = "character", na.strings = character(),
                              stringsAsFactors = FALSE)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    old <- .val_inputs$phase
    .val_inputs$phase <- "pre"
    on.exit(.val_inputs$phase <- old, add = TRUE)
    env <- new.env(parent = asNamespace("msrval"))
    env$out_dir <- out_dir
    sys.source(script, envir = env)
    files <- list.files(out_dir, recursive = TRUE)
    for (f in files) {
        s <- spec[spec$file == f, , drop = FALSE]
        if (nrow(s) == 0L || !grepl("\\.tsv$", f)) {
            valAbort("phase pre wrote an unlisted file {f}", "val_publish_forbidden")
        }
        labels <- lapply(stats::setNames(s$labels, s$column), function(x) {
            if (nzchar(x)) strsplit(x, ",", fixed = TRUE)[[1L]] else character()
        })
        valCheckResultTable(file.path(out_dir, f), f, columns = stats::setNames(s$type, s$column),
                            labels = labels[s$type == "label"], max_rows = 1000L)
    }
    TRUE
}

#' R1: record a UniProt cache from a committed accession list (10.4 R1, r1.12)
#'
#' Reads `accession` and `taxon` (one taxon) from the list, calls [valRecordUniprotCache()] into the
#' private `cache_dir`, and writes `uniprot.tsv` (taxon, accessions, rows, release,
#' release_source, sha256) into `out_dir`; the cache itself never enters `out_dir`. Refusals
#' (`val_annotation_incomplete`, `val_annotation_empty`) pass through and nothing is written.
#'
#' @param accessions Path of the accession list (`accession`, `taxon`).
#' @param out_dir The record folder.
#' @param cache_dir Private cache folder.
#' @param annotate,release As [valRecordUniprotCache()].
#' @return The record's row (a data frame), invisibly.
#' @export
valRecordUniprot <- function(accessions, out_dir, cache_dir, annotate = NULL, release = NULL) {
    acc <- utils::read.delim(accessions, colClasses = "character", stringsAsFactors = FALSE)
    taxon <- unique(acc$taxon)
    if (length(taxon) != 1L || !all(nzchar(acc$accession))) {
        valAbort("the accession list must name one taxon and no empty accession",
                 "val_manifest_invalid")
    }
    res <- valRecordUniprotCache(data.frame(Protein.Ids = acc$accession), "Protein.Ids", taxon,
                                 cache_dir, annotate = annotate, release = release)
    row <- data.frame(taxon = taxon, accessions = nrow(acc), rows = nrow(res$table),
                      release = res$release, release_source = res$release_source,
                      sha256 = res$sha256, stringsAsFactors = FALSE)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.table(row, file.path(out_dir, "uniprot.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE)
    invisible(row)
}
