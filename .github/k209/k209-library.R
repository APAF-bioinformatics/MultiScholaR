# KIT-209 (branch run-kit209-memory), gate B item B2: is each build's DESCRIPTION met by the library it runs on? For
# every Depends, Imports and LinkingTo entry of the installed build, the version the journey's library path finds and
# whether it meets the requirement. Writes <out tsv>: build, package, requirement, installed, met.
#
# Usage: Rscript --vanilla k209-library.R <build label> <build library> <out tsv>
args <- commandArgs(trailingOnly = TRUE)
build <- args[[1L]]
lib <- args[[2L]]
out <- args[[3L]]
paths <- c(lib, .libPaths())
description <- utils::packageDescription("MultiScholaR", lib.loc = lib)
fields <- unlist(lapply(c("Depends", "Imports", "LinkingTo"), function(field) {
    value <- description[[field]]
    if (is.null(value)) character() else trimws(strsplit(gsub("\n", " ", value), ",")[[1L]])
}))
fields <- fields[nzchar(fields)]
rows <- do.call(rbind, lapply(fields, function(entry) {
    name <- trimws(sub("\\(.*$", "", entry))
    requirement <- if (grepl("\\(", entry)) trimws(sub("^.*\\((.*)\\).*$", "\\1", entry)) else ""
    if (identical(name, "R")) {
        installed <- as.character(getRversion())
    } else {
        installed <- tryCatch(as.character(utils::packageVersion(name, lib.loc = paths)), error = function(e) NA_character_)
    }
    met <- if (is.na(installed)) FALSE else if (!nzchar(requirement)) TRUE else {
        op <- sub("^([<>=]+).*$", "\\1", requirement)
        version <- trimws(sub("^[<>=]+", "", requirement))
        isTRUE(do.call(op, list(package_version(installed), package_version(version))))
    }
    data.frame(build = build, package = name, requirement = requirement, installed = installed, met = met, stringsAsFactors = FALSE)
}))
utils::write.table(rows, out, sep = "\t", row.names = FALSE, quote = FALSE, append = file.exists(out), col.names = !file.exists(out))
cat(sprintf("k209 B2 %s: %d requirements, %d not met\n", build, nrow(rows), sum(!rows$met)))
if (any(!rows$met)) print(rows[!rows$met, ], row.names = FALSE)
