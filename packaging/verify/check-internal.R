# KIT-134 P2 step 4 (the token build): the two internal packages in the bundle, checked from the
# moved folder with the bundle's own R, as check-library.R is run. Build-time check, not package code.
#
# Usage: Rscript --vanilla check-internal.R <bundle directory> <required|attempt>
#
# For each of fastduckdb and fastlimpadpc: bundled or not, version and commit, loads from the bundle.
# fastduckdb: its reader (read_parquet_fast, which the lipidomics and metabolomics partitioned steps
# use for tables wider than 64 columns) reads back a 100-column Parquet file written by duckdb.
# Linux: the DIA manifest MultiScholaR's configure wrote (supported, native_query, the fastduckdb
# writer), whether the DIA rollup and the direct writer admit the build, and the DPC engine the
# DIA_LIMPA route resolves. "required" fails when a package is missing or any check fails;
# "attempt" (Windows, D4) reports and fails only when a bundled package does not load.

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L || !arguments[[2L]] %in% c("required", "attempt")) {
    stop("usage: Rscript --vanilla check-internal.R <bundle directory> <required|attempt>")
}
bundle <- normalizePath(arguments[[1L]], winslash = "/", mustWork = TRUE)
required <- identical(arguments[[2L]], "required")
is_windows <- .Platform$OS.type == "windows"
canonical <- function(paths) {
    paths <- normalizePath(paths, winslash = "/", mustWork = FALSE)
    if (is_windows) tolower(paths) else paths
}
failures <- character()
check <- function(name, passed, detail = "", counts = TRUE) {
    cat(sprintf("%s %s%s\n", if (passed) "ok  " else if (counts) "FAIL" else "note", name,
        if (nzchar(detail)) paste0(": ", detail) else ""))
    if (!passed && counts) failures <<- c(failures, name)
}
message_of <- function(expr) tryCatch({ force(expr); "admitted" }, error = function(e) paste("refused:", conditionMessage(e)))

for (package in c("fastduckdb", "fastlimpadpc")) {
    path <- find.package(package, lib.loc = .libPaths(), quiet = TRUE)
    if (!length(path)) {
        check(sprintf("%s is bundled", package), FALSE, "not in the bundle's library", counts = required)
        next
    }
    description <- utils::packageDescription(package, lib.loc = .libPaths())
    loaded <- tryCatch({ loadNamespace(package); "" }, error = function(e) conditionMessage(e))
    check(sprintf("%s loads from the bundle", package),
        !nzchar(loaded) && startsWith(canonical(getNamespaceInfo(package, "path")), paste0(canonical(bundle), "/")),
        sprintf("version %s, %s/%s%s@%s%s", description$Version, description$RemoteUsername %||% "?",
            description$RemoteRepo %||% "?", if (is.null(description$RemoteSubdir)) "" else paste0("/", description$RemoteSubdir),
            description$RemoteSha %||% "?", if (nzchar(loaded)) paste(";", loaded) else ""))
}

if ("fastduckdb" %in% loadedNamespaces()) {
    columns <- stats::setNames(lapply(seq_len(99L), function(i) as.numeric(seq_len(500L)) / i), sprintf("v%03d", seq_len(99L)))
    expected <- data.frame(id = sprintf("row%04d", seq_len(500L)), columns, stringsAsFactors = FALSE)
    file <- tempfile(fileext = ".parquet")
    roundTrip <- function() {
        con <- DBI::dbConnect(duckdb::duckdb())
        on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
        duckdb::duckdb_register(con, "expected", expected)
        DBI::dbExecute(con, sprintf("COPY expected TO '%s' (FORMAT parquet)", gsub("\\\\", "/", file)))
        as.data.frame(fastduckdb::read_parquet_fast(file, lazy_strings = FALSE))
    }
    read_back <- tryCatch(roundTrip(), error = function(e) conditionMessage(e))
    same <- is.data.frame(read_back) && identical(names(read_back), names(expected)) &&
        all(vapply(names(expected), function(column) identical(as.vector(read_back[[column]]), expected[[column]]), logical(1)))
    check("fastduckdb's reader reads back a 100-column Parquet file written by duckdb, identical", same,
        if (is.data.frame(read_back)) sprintf("%d rows, %d columns", nrow(read_back), ncol(read_back)) else read_back,
        counts = required)
}

if (!is_windows) {
    manifest <- readRDS(system.file("duckdb", "linux_amd64", "manifest.rds", package = "MultiScholaR"))
    check("the DIA manifest: supported", isTRUE(manifest$supported), "", counts = required)
    check("the DIA manifest: native_query present", !is.null(manifest$native_query),
        if (is.null(manifest$native_query)) "absent (fastduckdb has no query-manifest.rds)" else manifest$native_query$api,
        counts = required)
    check("the DIA manifest: the fastduckdb writer recorded", !is.null(manifest$fastduckdb),
        if (is.null(manifest$fastduckdb)) "absent (fastduckdb exports no parquet_writer_capabilities)" else "", counts = required)
    rollup <- message_of(MultiScholaR:::diaRollupNativeInfo())
    check("the DIA precursor rollup admits this build", identical(rollup, "admitted"), rollup, counts = required)
    writer <- message_of(MultiScholaR:::workflowFastWriterInfo())
    check("the direct Parquet writer admits this build", identical(writer, "admitted"), writer, counts = FALSE)
    engine <- tryCatch(MultiScholaR:::resolveLimpaDpcEngine(), error = function(e) paste("refused:", conditionMessage(e)))
    check("the DIA_LIMPA route resolves the fastlimpadpc engine", identical(engine, "fastlimpadpc"), engine, counts = required)
}

if (length(failures)) {
    cat(sprintf("KIT-134 CHECK FAILED: %s\n", paste(failures, collapse = "; ")))
    quit(save = "no", status = 1L)
}
cat("KIT-134 CHECK PASSED: the internal packages in the bundle\n")
