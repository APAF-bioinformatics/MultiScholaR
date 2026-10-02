# Build only at package installation, never at analysis runtime.
destination <- file.path("inst", "duckdb", "linux_amd64")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
supported <- identical(Sys.info()[["sysname"]], "Linux") &&
    identical(Sys.info()[["machine"]], "x86_64") && .Machine$longdouble.digits == 64
manifest <- list(schema = 3L, supported = supported, platform = "linux_amd64",
    duckdb = "1.5.5", duckplyr = "1.2.1", api = "v1.2.0", longdouble_digits = 64L,
    runtime_release_status = "held",
    release_holds = c("scoped_unsigned_exact_sum", "pinned_private_duckdb_relation_bridge"),
    configure_invocation = '"${R_HOME}/bin/Rscript" --vanilla src/build_dia_rsum.R')
manifest$fastduckdb <- NULL
manifest$dependencies <- lapply(c("duckdb", "duckplyr"), function(package) {
    root <- system.file(package = package)
    paths <- list.files(root, "[.](rdb|rdx|so)$", recursive = TRUE)
    if (!nzchar(root) || !length(paths)) stop("exact-sum dependency artifacts are missing")
    list(package = package, files = stats::setNames(vapply(paths,
        function(path) digest::digest(file = file.path(root, path), algo = "sha256"),
        character(1)), paths))
})
if (requireNamespace("fastduckdb", quietly = TRUE) &&
    "parquet_writer_capabilities" %in% getNamespaceExports("fastduckdb")) {
    capability <- fastduckdb::parquet_writer_capabilities()
    binary <- getLoadedDLLs()[["fastduckdb"]][["path"]]
    writer_root <- system.file(package = "fastduckdb")
    writer_files <- c("R/fastduckdb.rdb", "R/fastduckdb.rdx", "NAMESPACE")
    writer_hashes <- stats::setNames(vapply(writer_files, function(path)
        digest::digest(file = file.path(writer_root, path), algo = "sha256"), character(1)), writer_files)
    manifest$fastduckdb <- list(capability = capability, files = writer_hashes,
        binary_sha256 = digest::digest(file = binary, algo = "sha256"))
}
if (supported) {
    compiler <- Sys.which("cc")
    if (!nzchar(compiler)) stop("DIA exact-sum installation requires a C compiler")
    flags <- c("-std=c11", "-O2", "-fPIC", "-shared", "-fno-fast-math", "-Wall", "-Wextra", "-Werror",
        "-DDUCKDB_EXTENSION_API_VERSION_MAJOR=1", "-DDUCKDB_EXTENSION_API_VERSION_MINOR=2",
        "-DDUCKDB_EXTENSION_API_VERSION_PATCH=0", "-I", "src")
    library <- tempfile("dia_rsum-", fileext = ".so")
    status <- system2(compiler, c(flags, "src/dia_rsum.c", "-o", shQuote(library)))
    if (status != 0L) stop("DIA exact-sum extension compilation failed")
    output <- file.path(destination, "dia_rsum.duckdb_extension")
    stopifnot(file.copy(library, output, overwrite = TRUE))
    # DuckDB C extension footer: upstream extension-ci-tools ef15a2a format.
    con <- file(output, "ab")
    writeBin(c(as.raw(c(0, 147, 4, 16)), charToRaw("duckdb_signature"), as.raw(c(128, 4))), con)
    for (field in c("", "", "", "C_STRUCT", "dia025-v1", "v1.2.0", "linux_amd64", "4")) {
        bytes <- charToRaw(field)
        writeBin(c(bytes, raw(32L - length(bytes))), con)
    }
    writeBin(raw(256L), con)
    close(con)
    unlink(library)
    manifest$binary_sha256 <- digest::digest(file = output, algo = "sha256")
    inputs <- c("src/dia_rsum.c", "src/duckdb.h", "src/duckdb_extension.h",
        "src/build_dia_rsum.R", "configure")
    manifest$source_sha256 <- stats::setNames(vapply(inputs,
        function(path) digest::digest(file = path, algo = "sha256"), character(1)), basename(inputs))
    manifest$compiler <- system2(compiler, "--version", stdout = TRUE)[[1L]]
    manifest$flags <- flags
}
saveRDS(manifest, file.path(destination, "manifest.rds"), compress = FALSE, version = 3)
