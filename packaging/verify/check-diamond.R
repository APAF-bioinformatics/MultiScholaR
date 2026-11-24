# KIT-167: the bundled DIAMOND in a bare container (network off, the bundle read-only under a path with spaces): the
# launcher's MSR_DIAMOND points at it, `diamond version` runs, and the package's own check (SHA-256 against the pin,
# then the version) accepts it. Prints `KIT167 <field> <value>` lines; exits 1 when a check fails.
arguments <- commandArgs(trailingOnly = TRUE)
bundle <- arguments[[1L]]
path <- Sys.getenv("MSR_DIAMOND")
cat(sprintf("KIT167 msr_diamond %s\n", path))
expected <- file.path(bundle, "diamond", "diamond")
ok <- identical(normalizePath(path, mustWork = FALSE), normalizePath(expected, mustWork = FALSE)) && file.exists(path)
version <- if (ok) system2(path, "version", stdout = TRUE) else character()
cat(sprintf("KIT167 diamond_version %s\n", paste(version, collapse = " ")))
checked <- tryCatch(MultiScholaR:::checkDiamond(), error = function(e) {
    cat(sprintf("KIT167 check_diamond_refused %s\n", conditionMessage(e)))
    NULL
})
if (!is.null(checked)) cat(sprintf("KIT167 check_diamond_sha256 %s\n", checked$sha256))
if (!ok || !length(version) || is.null(checked)) quit(status = 1L)
