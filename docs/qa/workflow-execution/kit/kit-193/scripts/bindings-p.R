# KIT-193 section 8.2: the digests of the objects bound by the pairs of targets 4 (library parts) and 5 and of group R
# (all but R-I1, the lipid cases and R-N1's sweep helper), and the versions of the libraries they bind, on the installed
# package. Functions: `digest::digest(deparse(body(fn)), algo = "sha256")` (for the S4 method, the body of
# `methods::getMethod()`); constants: `digest::digest(deparse(x), algo = "sha256")`; libraries:
# `utils::packageVersion()`. Run by run-tests after the package is installed.
# Usage: Rscript docs/qa/workflow-execution/kit/kit-193/scripts/bindings-p.R <out dir>
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1L]] else "kit193-out"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
ns <- asNamespace("MultiScholaR")
fns <- c("integrationEnrichmentPermutation", "integrationEnrichmentSetT", "integrationEnrichmentFloorMinimum",
    "restoreRuvFeatureLevels", "computeRuvCovariates", "computeOmicsRuvAssayFactors", "computeProtRuvModelFactors",
    "ruvExactContrastsFit", "protRuvFirstPassControls", "omicsRuvFirstPassControls", "trxRuvDesignControls")
consts <- c("INTEGRATION_ENRICHMENT_PERMUTATION_TIE", "OMICS_RUV_FIRST_PASS_DROPPED_SHARE", "TRX_FIRST_PASS_DROPPED_SHARE")
libs <- c("qvalue", "stats", "ruv", "RUVIIIC", "limma", "edgeR")
digest_of <- function(x) digest::digest(deparse(x), algo = "sha256")
method <- methods::getMethod("ruvIII_C_Varying", "ProteinQuantitativeData", where = ns)
rows <- data.frame(
    name = c(fns, "ruvIII_C_Varying,ProteinQuantitativeData", consts, libs),
    kind = c(rep("package function", length(fns) + 1L), rep("package constant", length(consts)),
        rep("library", length(libs))),
    value = c(vapply(fns, function(name) digest_of(body(get(name, envir = ns))), character(1L)), digest_of(body(method)),
        vapply(consts, function(name) digest_of(get(name, envir = ns)), character(1L)),
        vapply(libs, function(name) format(utils::packageVersion(name)), character(1L))),
    stringsAsFactors = FALSE)
utils::write.table(rows, file.path(out, "bindings-p.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(R.version.string, paste("MultiScholaR", format(utils::packageVersion("MultiScholaR"))),
    paste("digest", format(utils::packageVersion("digest")))), file.path(out, "versions-p.txt"))
print(rows)
