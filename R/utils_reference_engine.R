# KIT-167 section 8.1, 8.4, 8.6: the sequence search engine (DIAMOND v2.2.8) as a separate program.
#   - `readDiamondPin()`: the pin file `inst/extdata/diamond-pin.tsv` (version, release assets and their SHA-256, the
#     executable's SHA-256 per platform), read by R and by the bundle build.
#   - `findDiamond()`: `MSR_DIAMOND` (set by the desktop shell), the option `MultiScholaR.diamond_path`, the path.
#   - `checkDiamond()`: the executable's SHA-256 against the pin once per session, then `diamond version`; on Windows
#     the exit statuses 0xC0000135 (DLL not found) and 0xC0000139 (entry point not found) are the missing runtime.
#   - `runDiamondCommand()`: one call with a time limit, in a scratch folder that is its working directory.
#   - `diamondSearchArguments()`: the shipped search settings of section 5.1.

#' The DIAMOND `--outfmt 6` fields every rule of the resolver uses (section 5.1)
#' @noRd
diamondOutputFields <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart",
    "send", "evalue", "bitscore", "score", "qlen", "slen", "qcovhsp", "scovhsp")

#' The determinism pins (section 5.1): fixed values passed explicitly, chosen at gate F2 on time and memory only.
#' Both proteomes are far below one block of 0.4 billion letters, so every search runs as one block.
#' @noRd
diamondPins <- c(block_size = "0.4", index_chunks = "1", algo = "0")

#' The Windows exit statuses of a program that cannot load its runtime: 0xC0000135 (a DLL not found) and 0xC0000139
#' (an entry point not found: a runtime older than the build's toolset), as signed and unsigned 32-bit values.
#' @noRd
diamondRuntimeMissingStatuses <- c(-1073741515, -1073741511, 3221225781, 3221225785)

#' Per-session record of the engine checks (the SHA-256 is computed once per executable and session)
#' @noRd
referenceEngineCache <- new.env(parent = emptyenv())

#' Forget the session's engine check (tests; a changed executable is re-checked by its size and time stamp anyway).
#' @noRd
resetDiamondCheck <- function() {
    rm(list = ls(referenceEngineCache, all.names = TRUE), envir = referenceEngineCache)
    invisible(NULL)
}

#' The platform name the pin file uses for this computer
#' @noRd
referenceEnginePlatform <- function() {
    if (identical(.Platform$OS.type, "windows")) return("windows")
    if (identical(unname(Sys.info()[["sysname"]]), "Linux")) return("linux")
    "other"
}

#' Read the engine pin file
#' @param path The pin file; the shipped one by default.
#' @return Data frame: platform, version, archive_url, archive_bytes, archive_sha256, executable_path,
#'   executable_sha256 (NA until the gate that measures it has run).
#' @noRd
readDiamondPin <- function(path = system.file("extdata", "diamond-pin.tsv", package = "MultiScholaR")) {
    if (!nzchar(path) || !file.exists(path)) {
        cli::cli_abort("the engine pin file is missing from the installation", class = "reference_mapping_internal")
    }
    utils::read.delim(path, sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings = "NA",
        colClasses = c(platform = "character", version = "character", archive_url = "character",
            archive_bytes = "numeric", archive_sha256 = "character", executable_path = "character",
            executable_sha256 = "character"))
}

#' Where the engine is: `MSR_DIAMOND`, then the option `MultiScholaR.diamond_path`, then `diamond` on the path
#' @return The path, or "" when none is set or found.
#' @noRd
findDiamond <- function() {
    from_env <- Sys.getenv("MSR_DIAMOND", "")
    if (nzchar(from_env)) return(from_env)
    from_option <- getOption("MultiScholaR.diamond_path", "")
    if (is.character(from_option) && length(from_option) == 1L && nzchar(from_option)) return(from_option)
    unname(Sys.which("diamond"))
}

#' Run one engine command with a time limit
#' @param diamond The executable.
#' @param args Its arguments (paths relative to `scratch`).
#' @param scratch The working directory of the call.
#' @param timeout Seconds; past it the call is killed and `reference_mapping_timeout` is raised.
#' @param run_fn The runner (`processx::run`; tests pass a stub).
#' @return The runner's list (`status`, `stdout`, `stderr`).
#' @noRd
runDiamondCommand <- function(diamond, args, scratch, timeout = 600, run_fn = processx::run) {
    minutes <- max(1L, as.integer(ceiling(timeout / 60)))
    result <- tryCatch(
        run_fn(diamond, args, wd = scratch, timeout = timeout, error_on_status = FALSE),
        system_command_timeout_error = function(e) list(status = NA_integer_, stdout = "", stderr = "", timeout = TRUE)
    )
    if (isTRUE(result$timeout)) referenceMappingRefuse("reference_mapping_timeout", minutes)
    status <- result$status
    if (isTRUE(as.numeric(status) %in% diamondRuntimeMissingStatuses)) {
        referenceMappingRefuse("reference_mapping_engine_runtime_missing", fields = list(status = status))
    }
    if (!identical(as.integer(status), 0L)) {
        detail <- utils::tail(strsplit(paste(result$stderr %||% "", collapse = "\n"), "\n", fixed = TRUE)[[1]], 3L)
        referenceMappingRefuse("reference_mapping_engine_failed", as.character(status),
            paste(detail, collapse = " "), fields = list(status = status))
    }
    result
}

#' Check the engine before a new mapping (section 8.4)
#' @param path The executable (`findDiamond()`).
#' @param pin The pin file's table.
#' @param run_fn The runner, passed to `runDiamondCommand`.
#' @return `list(path, version, sha256, platform)`; refuses with `reference_mapping_engine_missing`,
#'   `reference_mapping_engine_mismatch` or `reference_mapping_engine_runtime_missing`.
#' @noRd
checkDiamond <- function(path = findDiamond(), pin = readDiamondPin(), run_fn = processx::run) {
    if (!nzchar(path) || !file.exists(path)) referenceMappingRefuse("reference_mapping_engine_missing")
    platform <- referenceEnginePlatform()
    row <- pin[pin$platform == platform, , drop = FALSE]
    expected <- if (nrow(row) == 1L) row$executable_sha256 else NA_character_
    version <- if (nrow(row) == 1L) row$version else NA_character_
    info <- file.info(path)
    key <- paste(normalizePath(path, winslash = "/"), info$size, as.numeric(info$mtime), expected, sep = "|")
    if (exists(key, envir = referenceEngineCache, inherits = FALSE)) return(get(key, envir = referenceEngineCache))
    sha256 <- digest::digest(file = path, algo = "sha256")
    if (is.na(expected) || !identical(sha256, expected)) {
        referenceMappingRefuse("reference_mapping_engine_mismatch", path, version %||% "",
            fields = list(sha256 = sha256, expected = expected))
    }
    result <- runDiamondCommand(path, "version", scratch = tempdir(), timeout = 60, run_fn = run_fn)
    reported <- regmatches(result$stdout, regexpr("[0-9]+\\.[0-9]+\\.[0-9]+", result$stdout))
    if (!identical(reported, version)) {
        referenceMappingRefuse("reference_mapping_engine_mismatch", path, version,
            fields = list(sha256 = sha256, expected = expected, reported = reported))
    }
    checked <- list(path = path, version = version, sha256 = sha256, platform = platform)
    assign(key, checked, envir = referenceEngineCache)
    checked
}

#' The arguments of one search at the shipped settings (section 5.1)
#' @param query,db,out Relative names inside the scratch folder.
#' @param evalue The e-value setting.
#' @param threads Threads, or NULL to leave DIAMOND's default (results do not depend on it).
#' @return Character vector of arguments.
#' @noRd
diamondSearchArguments <- function(query, db, out, evalue, threads = NULL) {
    c("blastp", "--query", query, "--db", db, "--out", out, "--ultra-sensitive", "--masking", "none", "--ext",
        "full", "--comp-based-stats", "1", "-k0", "--max-hsps", "1", "--evalue", format(evalue),
        "--block-size", diamondPins[["block_size"]], "--index-chunks", diamondPins[["index_chunks"]],
        "--algo", diamondPins[["algo"]], if (!is.null(threads)) c("--threads", as.character(threads)),
        "--tmpdir", "tmp", "--outfmt", "6", diamondOutputFields)
}

#' Read one DIAMOND tabular output with the fields of `diamondOutputFields`
#' @param path The output file.
#' @return Data frame; identifiers as text, every statistic a double parsed once from its printed text.
#' @noRd
readDiamondHits <- function(path) {
    numeric <- setdiff(diamondOutputFields, c("qseqid", "sseqid"))
    empty <- as.data.frame(stats::setNames(c(list(character(), character()),
        rep(list(numeric()), length(numeric))), diamondOutputFields), stringsAsFactors = FALSE)
    if (!file.exists(path) || file.size(path) == 0) return(empty)
    utils::read.delim(path, header = FALSE, col.names = diamondOutputFields, quote = "", comment.char = "",
        colClasses = stats::setNames(c("character", "character", rep("numeric", length(numeric))),
            diamondOutputFields), stringsAsFactors = FALSE)
}
