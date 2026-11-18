# run-tests: the package tree made ready before its test files run (from the package tree's root).
#
# Usage: Rscript --vanilla prepare.R <library table out> [load]
# Names the DESCRIPTION dependencies (Depends, Imports, LinkingTo) the library lacks (a warning; the files that need
# them will fail), writes the library table (package, version, library), and with `load` compiles the tree's native
# code once with pkgload::load_all, outside the files' time and logs.

arguments <- commandArgs(trailingOnly = TRUE)
description <- read.dcf("DESCRIPTION")
fields <- intersect(c("Depends", "Imports", "LinkingTo"), colnames(description))
declared <- unlist(lapply(fields, function(field) strsplit(description[, field], ",")[[1L]]))
declared <- setdiff(trimws(sub("\\(.*$", "", trimws(declared))), c("", "R"))
installed <- utils::installed.packages()
missing <- setdiff(declared, rownames(installed))
if (length(missing)) cat(sprintf("::warning::the library lacks the package's dependencies: %s\n", paste(missing, collapse = ", ")))
cat(sprintf("%d DESCRIPTION dependencies, %d missing\n", length(declared), length(missing)))
utils::write.table(data.frame(package = installed[, "Package"], version = installed[, "Version"],
    library = installed[, "LibPath"], row.names = NULL), arguments[[1L]], sep = "\t", quote = FALSE, row.names = FALSE)
if (length(arguments) > 1L && identical(arguments[[2L]], "load")) {
    suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = FALSE))
    cat("MultiScholaR", as.character(utils::packageVersion("MultiScholaR")), "loads from the tree\n")
}
