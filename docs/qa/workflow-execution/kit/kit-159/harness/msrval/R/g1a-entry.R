# G1(a)'s job pieces (KIT-159 section 2; KIT-199 12.3, D7 / VAL-120 (2)): the job-environment
# recorder, the lockfile check, the DIA-with-limpa DA inputs of part (2), `gates.md`, and L1's
# committed DA tables for G13 (VAL-133 (1)). The entry script is `val-g1a.R`.

#' The environment keys a G1(a) job records, in KIT-199's `environment.tsv` order (its
#' `control.R` `writeEnvironment`), plus the lockfile's SHA-256
#' @noRd
G1A_JOB_ENV_KEYS <- c("cpu_model", "cpu_count", "image_os", "image_version", "c_library_version",
                      "libblas3_package", "liblapack3_package", "blas_library", "lapack_library",
                      "lapack_version", "r_version", "rust_toolchain_version", "pin_fastduckdb",
                      "pin_fastlimpadpc", "openblas_num_threads", "omp_num_threads",
                      "lockfile_sha256")

#' The lockfile G1(a) and D7 read
#' @noRd
G1A_LOCKFILE <- "packaging/multischolar.lock.json"

#' The first line of a command's output, or `not found`
#' @noRd
valCommandValue <- function(command, args, pattern = NULL) {
    path <- Sys.which(command)
    if (!nzchar(path)) return("not found")
    out <- tryCatch(processx::run(path, args, error_on_status = FALSE)$stdout,
                    error = function(e) "")
    out <- trimws(strsplit(out, "\n", fixed = TRUE)[[1L]][1L] %||% "")
    if (!is.null(pattern)) out <- regmatches(out, regexpr(pattern, out))
    if (!length(out) || is.na(out) || !nzchar(out)) "not found" else out
}

#' A commit pin of an internal package in DESCRIPTION (`<pkg>=...@<40 hex>`)
#' @noRd
valInternalPin <- function(pkg_root, package) {
    text <- readLines(file.path(pkg_root, "DESCRIPTION"), warn = FALSE)
    line <- grep(paste0("^\\s*", package, "="), text, value = TRUE)
    pin <- regmatches(line, regexpr("@[0-9a-f]{40}", line))
    if (length(pin)) substring(pin[[1L]], 2L) else "not found"
}

#' The G1(a) job's environment (D7: the CPU model, the image version, the ten keys of KIT-199's
#' corrected I04 text, and the lockfile's digest), derived as KIT-199's `writeEnvironment` derives
#' them so that the values compare equal with a link's `environment.tsv`
#'
#' @param pkg_root The package checkout.
#' @param cpuinfo Path of `/proc/cpuinfo`.
#' @param command `function(command, args, pattern)` reading a command's output (tests).
#' @return A named character vector over `G1A_JOB_ENV_KEYS`.
#' @export
valJobEnvironment <- function(pkg_root = ".", cpuinfo = "/proc/cpuinfo",
                              command = valCommandValue) {
    env_or <- function(name) {
        v <- Sys.getenv(name)
        if (nzchar(v)) v else "<unset>"
    }
    info <- if (file.exists(cpuinfo)) readLines(cpuinfo, warn = FALSE) else character()
    model <- sub("^model name\\s*:\\s*", "", grep("^model name", info, value = TRUE))
    rust <- Sys.which("rustc")
    if (!nzchar(rust)) rust <- file.path(Sys.getenv("HOME"), ".cargo", "bin", "rustc")
    rust_version <- if (file.exists(rust)) {
        command(rust, "--version", "[0-9]+\\.[0-9]+\\.[0-9]+")
    } else {
        "not found"
    }
    lock <- file.path(pkg_root, G1A_LOCKFILE)
    values <- c(
        cpu_model = if (length(model)) model[[1L]] else "not found",
        cpu_count = as.character(parallel::detectCores()),
        image_os = env_or("ImageOS"), image_version = env_or("ImageVersion"),
        c_library_version = command("getconf", "GNU_LIBC_VERSION", "[0-9]+\\.[0-9]+"),
        libblas3_package = command("dpkg-query", c("-W", "-f=${Version}", "libblas3"), NULL),
        liblapack3_package = command("dpkg-query", c("-W", "-f=${Version}", "liblapack3"), NULL),
        blas_library = basename(normalizePath(extSoftVersion()[["BLAS"]], mustWork = FALSE)),
        lapack_library = basename(normalizePath(La_library(), mustWork = FALSE)),
        lapack_version = La_version(),
        r_version = R.version.string,
        rust_toolchain_version = rust_version,
        pin_fastduckdb = valInternalPin(pkg_root, "fastduckdb"),
        pin_fastlimpadpc = valInternalPin(pkg_root, "fastlimpadpc"),
        openblas_num_threads = env_or("OPENBLAS_NUM_THREADS"),
        omp_num_threads = env_or("OMP_NUM_THREADS"),
        lockfile_sha256 = if (file.exists(lock)) valSha256(lock) else "not found"
    )
    values[G1A_JOB_ENV_KEYS]
}

#' Write the job environment as `environment.tsv` (key, value; KIT-199's form)
#'
#' @param values From [valJobEnvironment()].
#' @param path Output path.
#' @return `TRUE` when every value was read (none `not found`), invisibly.
#' @export
valWriteJobEnvironment <- function(values, path) {
    clean <- gsub("[\t\r\n]+", " ", unname(values))
    utils::write.table(data.frame(key = names(values), value = clean), path, sep = "\t",
                       quote = FALSE, row.names = FALSE)
    invisible(all(nzchar(clean)) && !any(clean == "not found"))
}

#' Is the lockfile unchanged between two commits?
#'
#' @param repo The checkout (full history).
#' @param since,at The two commits.
#' @return `TRUE` or `FALSE`.
#' @export
valLockfileUnchanged <- function(repo, since, at = "HEAD") {
    identical(valGit(repo, "diff", "--quiet", since, at, "--", G1A_LOCKFILE)$status, 0L)
}

#' Part (2)'s DIA-with-limpa DA inputs (D7): the live table part (1) copied, the last link's
#' `da-runner.tsv` with its control jobs, and the mode for this job
#'
#' @param repo The checkout (full history).
#' @param link_dir The link's evidence folder (KIT-199's unless a later re-freezing link
#'   supplies its own).
#' @param job From [valJobEnvironment()].
#' @param live_table The live DA table part (1) copied after step `da`.
#' @param controls The link's control-run subfolders.
#' @return `list(table, link, mode, link_commit)`.
#' @export
valG1aDaInputs <- function(repo, link_dir, job, live_table, controls = c("r1", "r2")) {
    if (!file.exists(live_table)) valAbort("no live DA table was copied", "val_spec_invalid")
    link <- valReadDaLink(file.path(repo, link_dir, "da-runner.tsv"),
                          file.path(repo, link_dir, controls, "environment.tsv"),
                          file.path(repo, link_dir, controls, "da-runner.tsv"))
    commit <- valGit(repo, "log", "-1", "--format=%H", "--",
                     file.path(link_dir, "da-runner.tsv"))$out[[1L]]
    mode <- valG1aDaMode(job, link, valLockfileUnchanged(repo, commit, "HEAD"))
    list(table = paste(readLines(live_table, warn = FALSE), collapse = "\n"), link = link,
         mode = mode, link_commit = commit)
}

#' A data frame as a markdown table
#' @noRd
valMdTable <- function(x) {
    if (is.null(x) || !nrow(x)) return("(none)")
    esc <- function(v) gsub("|", "\\|", ifelse(is.na(v), "", as.character(v)), fixed = TRUE)
    c(paste0("| ", paste(names(x), collapse = " | "), " |"),
      paste0("|", paste(rep("---", ncol(x)), collapse = "|"), "|"),
      apply(x, 1L, function(r) paste0("| ", paste(esc(r), collapse = " | "), " |")))
}

#' Write G1(a)'s `gates.md` section for a job
#'
#' @param path Output path.
#' @param at The commit under test.
#' @param env From [valJobEnvironment()].
#' @param part1 From [valG1aPart1()].
#' @param part2 From [valG1aPart2()], or `NULL` (not run).
#' @param da From [valG1aDaInputs()], or `NULL`.
#' @return `path`, invisibly.
#' @export
valG1aGatesMd <- function(path, at, env, part1, part2 = NULL, da = NULL) {
    p1_ok <- all(part1$status == "identical")
    lines <- c(
        "## G1(a) oracle reproduction (written by val-g1a.R)", "",
        sprintf("Commit under test: %s.", at), "",
        "### The job's environment (D7 keys, CPU model, lockfile digest)", "",
        valMdTable(data.frame(key = names(env), value = unname(env))), "",
        "### Part (1): installed against `load_all`", "",
        valMdTable(part1), "",
        sprintf("Part (1) result: %s.", if (p1_ok) {
            "every compared step identical on every lane run"
        } else {
            "a lane differs or did not complete; G1(a) stops at its first differing step (D3)"
        }), "",
        "### Part (2): against the oracle of the commit under test", ""
    )
    lines <- c(lines, if (is.null(part2)) "Not run in this job." else valMdTable(part2), "")
    if (!is.null(da)) {
        lines <- c(lines, "### D7: the DIA-with-limpa DA table by value", "",
                   sprintf("Link commit %s; mode `%s`%s.", da$link_commit, da$mode$mode,
                           if (length(da$mode$reasons)) {
                               paste0(" (", paste(da$mode$reasons, collapse = "; "), ")")
                           } else {
                               ""
                           }), "")
    }
    writeLines(lines, path)
    invisible(path)
}

#' L1's committed DA tables for G13 (VAL-133 (1)): each route's committed entry, its table read
#' at `at` and checked against its SHA-256
#'
#' The request's `l1_path` file holds per route `digest`, `fingerprint`, `table` (a committed
#' path), `table_sha256`, `keys`, `q_col` and `fc_col`; a missing field or a table whose hash
#' differs is `val_spec_invalid`.
#'
#' @param req The parity request.
#' @param repo The checkout (full history).
#' @param at The commit to read at.
#' @return Named list by route: `digest`, `fingerprint`, `table` (text), `keys`, `q_col`,
#'   `fc_col`, ready for [valCompareParityDa()].
#' @export
valParityL1Tables <- function(req, repo = ".", at = req[["release_commit"]]) {
    bad <- function(what) valAbort("parity request: {what}", "val_spec_invalid")
    valCheckParityL1(req, repo, at)
    committed <- jsonlite::fromJSON(valGitShow(repo, at, req[["l1_path"]]),
                                    simplifyVector = FALSE)
    lapply(stats::setNames(names(req[["l1"]]), names(req[["l1"]])), function(route) {
        w <- committed[[route]]
        need <- c("table", "table_sha256", "keys", "q_col", "fc_col")
        if (!all(need %in% names(w))) bad(paste(route, "lacks", setdiff(need, names(w))[[1L]]))
        sha <- valGitSha(repo, at, w[["table"]])
        if (!identical(sha, w[["table_sha256"]])) {
            bad(paste("the L1 table of", route, "does not have its committed SHA-256"))
        }
        list(digest = w[["digest"]], fingerprint = w[["fingerprint"]],
             table = valGitShow(repo, at, w[["table"]]), keys = unlist(w[["keys"]]),
             q_col = w[["q_col"]], fc_col = w[["fc_col"]])
    })
}
