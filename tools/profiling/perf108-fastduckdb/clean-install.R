# Existing PERF-108 guarded runner entry points for the published install recipe.
runPerf108CleanInstall <- function(stage, run_dir, source, library_path) {
    request <- jsonlite::read_json(file.path(source,
        "tools/profiling/perf108-fastduckdb/clean-install-request.json"))
    hash <- function(path) digest::digest(file = path, algo = "sha256")
    if (stage == "UB") {
        stopifnot(identical(hash(request$archive), request$archive_sha256),
            identical(hash(request$installer), request$installer_sha256))
        installer <- file.path(run_dir, "install-fastduckdb-0.2.0.sh")
        archive <- file.path(run_dir, "fastduckdb_0.2.0.tar.gz")
        stopifnot(file.copy(request$installer, installer), file.copy(request$archive, archive))
        destination <- file.path(run_dir, "library")
        work <- file.path(run_dir, "clean-install")
        status <- withr::with_envvar(c(R_LIBS = paste(c(library_path, .libPaths()), collapse = .Platform$path.sep)),
            system2("bash", shQuote(c(installer, archive, destination, work))))
        if (!identical(status, 0L)) stop("published clean installation failed: ", status)
        # Only bootstrap the unchanged scientific provider/test dependencies after
        # the new fastduckdb package has independently installed successfully.
        for (package in c("MultiScholaR", "duckplyr", "nanoarrow"))
            stopifnot(file.copy(file.path(library_path, package), destination, recursive = TRUE))
        receipt <- list(status = "complete", library = destination, work = work,
            archive_sha256 = hash(archive), installer_sha256 = hash(installer),
            lock_sha256 = hash(file.path(work, "fastduckdb/src/rust/Cargo.lock")),
            provenance = readLines(file.path(work, "install-provenance.txt")))
        saveRDS(receipt, file.path(run_dir, "clean-install-result.rds"))
        return(receipt)
    }
    stopifnot(identical(stage, "UT"), identical(normalizePath(find.package("fastduckdb")),
        normalizePath(file.path(library_path, "fastduckdb"))))
    build <- readRDS(file.path(dirname(library_path), "clean-install-result.rds"))
    if (!is.null(request$generic_receipt)) {
        stopifnot(identical(normalizePath(request$generic_library), normalizePath(library_path)),
            identical(hash(request$generic_receipt), request$generic_sha256))
        checks <- readRDS(request$generic_receipt)
    } else {
        checks <- testthat::test_dir(file.path(build$work, "fastduckdb/tests/testthat"),
            package = "fastduckdb", load_package = "installed", reporter = "summary",
            stop_on_failure = FALSE)
        checks <- as.data.frame(checks)
    }
    saveRDS(checks, file.path(run_dir, "generic-tests.rds"))
    stopifnot(!any(checks$failed > 0 | checks$error | checks$skipped))
    scope <- fastduckdb::open_query_fast()
    on.exit(fastduckdb::close_query_fast(scope), add = TRUE)
    fastduckdb::register_query_native_fast(scope, "MultiScholaR", "dia_ordered_sum_v1")
    actual <- fastduckdb::query_data_fast(scope, paste(
        "SELECT g, dia_r_sum(list(x ORDER BY o)) AS quantity FROM (VALUES",
        "(1,1,1e16::DOUBLE),(1,2,1::DOUBLE),(1,3,1::DOUBLE),",
        "(2,1,NULL::DOUBLE),(2,2,NULL::DOUBLE),",
        "(3,1,0::DOUBLE),(3,2,NULL::DOUBLE),(3,3,2::DOUBLE)) t(g,o,x)",
        "GROUP BY g ORDER BY g"))
    expected <- data.frame(g = 1:3, quantity = c(sum(c(1e16, 1, 1)), NA_real_, 2))
    stopifnot(identical(actual, expected))
    fastduckdb::close_query_fast(scope)
    stopifnot(!fastduckdb::query_is_valid_fast(scope))
    files <- list.files(find.package("fastduckdb"), recursive = TRUE, full.names = TRUE)
    files <- files[grepl("[.](so|rdb|rdx)$|/DESCRIPTION$|/NAMESPACE$|manifest[.]rds$", files)]
    receipt <- list(status = "complete", generic_assertions = sum(checks$passed),
        native_registration = "passed; unchanged MultiScholaR provider loaded against fresh fastduckdb",
        ordered_sum = actual, library = library_path, build = build,
        installed_hashes = stats::setNames(vapply(files, hash, character(1)), files),
        generic_receipt_reused = request$generic_receipt,
        known_limit = "all-missing literal-list provider call fails in old and new builds; grouped P14 call passes; separate open finding",
        performance_qualification_transferred = FALSE)
    saveRDS(receipt, file.path(run_dir, "clean-install-check.rds"))
    receipt
}
