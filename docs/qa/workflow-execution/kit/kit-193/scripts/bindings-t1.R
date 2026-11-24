# KIT-193 section 8.2: the digests of target 1's bound objects on the installed package, for bindings.tsv rows T1
# (`digest::digest(deparse(body(fn)), algo = "sha256")`). Run by run-tests after the package is installed.
# Usage: Rscript docs/qa/workflow-execution/kit/kit-193/scripts/bindings-t1.R <out dir>
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1L]] else "kit193-out"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
ns <- asNamespace("MultiScholaR")
fns <- c("exceedsRoundedPercent", "roundedPercentNumerator", "allowedByRoundedPercent", "updateMissingValueParameters")
digests <- vapply(fns, function(name) digest::digest(deparse(body(get(name, envir = ns))), algo = "sha256"),
    character(1L))
rows <- data.frame(name = fns, digest = unname(digests), stringsAsFactors = FALSE)
utils::write.table(rows, file.path(out, "bindings-t1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(R.version.string, paste("MultiScholaR", format(utils::packageVersion("MultiScholaR"))),
    paste("digest", format(utils::packageVersion("digest")))), file.path(out, "versions.txt"))
print(rows)
