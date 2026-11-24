# KIT-193 section 8.2: the digest of row T4's bound function `processEnrichments` on the installed package after
# KIT-203's merge (`digest::digest(deparse(body(fn)), algo = "sha256")`). Run by run-tests after the package is
# installed. Usage: Rscript docs/qa/workflow-execution/kit/kit-193/scripts/bindings-t4.R <out dir>
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1L]] else "kit193-out"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
fn <- get("processEnrichments", envir = asNamespace("MultiScholaR"))
rows <- data.frame(name = "processEnrichments", value = digest::digest(deparse(body(fn)), algo = "sha256"),
    stringsAsFactors = FALSE)
utils::write.table(rows, file.path(out, "bindings-t4.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(rows)
