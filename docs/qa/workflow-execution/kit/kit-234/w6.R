# KIT-234 W6: no workbook in the installed library or in the built tarball. Run by run-tests' scripts job, which
# installs package_commit with R CMD INSTALL (the installer's mechanism) and runs this script from the package root.
#   (i)  system.file("workbooks", package = "MultiScholaR") in the installed library: "" expected;
#   (ii) R CMD build --no-manual of the checkout into a temporary folder, its paths listed with utils::untar: none may
#        match ^MultiScholaR/(inst/workbooks|archive)(/|$).
# Exit 0 when both hold, 1 when either does not (a result), 2 when the instrument failed (MultiScholaR not installed,
# or R CMD build gave no tarball): no result.

installed_package <- system.file(package = "MultiScholaR")
if (!nzchar(installed_package)) {
    cat("W6: MultiScholaR is not installed in this library; no result\n")
    quit(save = "no", status = 2L)
}
cat("W6: MultiScholaR", format(utils::packageVersion("MultiScholaR")), "installed at", installed_package, "\n")

workbooks <- system.file("workbooks", package = "MultiScholaR")
cat("W6 (i): system.file(\"workbooks\", package = \"MultiScholaR\") = \"", workbooks, "\"\n", sep = "")
if (nzchar(workbooks)) {
    cat("W6 (i): files there:\n")
    writeLines(paste0("  ", list.files(workbooks, recursive = TRUE)))
}

root <- normalizePath(".", mustWork = TRUE)
build_dir <- tempfile("kit234-build-")
dir.create(build_dir)
r_binary <- file.path(R.home("bin"), "R")
build_status <- system(paste("cd", shQuote(build_dir), "&&", shQuote(r_binary), "CMD build --no-manual", shQuote(root)))
tarball <- list.files(build_dir, pattern = "^MultiScholaR_.*[.]tar[.]gz$", full.names = TRUE)
if (build_status != 0L || length(tarball) != 1L) {
    cat("W6 (ii): R CMD build exited", build_status, "and left", length(tarball), "tarballs; no result\n")
    quit(save = "no", status = 2L)
}
paths <- utils::untar(tarball, list = TRUE)
hits <- grep("^MultiScholaR/(inst/workbooks|archive)(/|$)", paths, value = TRUE)
cat("W6 (ii): ", basename(tarball), " holds ", length(paths), " paths; ", length(hits),
    " match ^MultiScholaR/(inst/workbooks|archive)(/|$)\n", sep = "")
if (length(hits) > 0L) writeLines(paste0("  ", sort(hits)))

passed <- !nzchar(workbooks) && length(hits) == 0L
cat("W6:", if (passed) "no workbook installed or built (exit 0)" else "a workbook is installed or built (exit 1)", "\n")
quit(save = "no", status = if (passed) 0L else 1L)
