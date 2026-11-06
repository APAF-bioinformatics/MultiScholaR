# KIT-099: the call profile's part (b) summarised offline from its stored Rprof file. call-profile.R's own summary
# (summaryRprof(memory = "both")) stopped on a sample line it could not parse, after the run had completed; the profile was
# collected, so it is read here, not rerun. Per candidate call: seconds with the call on the stack (samples x 0.05 s), its
# share of the run, and the largest R heap in use while it was on the stack (small + large vector heap, from the memory
# fields each sample carries; the process's resident peak is reported by the trial runner, not here).
# Usage: Rscript call-profile-summary.R <rprof.out> <call-profile.out> <output txt>
args <- commandArgs(TRUE)
lines <- readLines(args[[1L]])
interval <- as.numeric(sub(".*sample.interval=([0-9]+).*", "\\1", lines[[1L]])) / 1e6
samples <- lines[-1L]
samples <- samples[grepl("^:[0-9]+:[0-9]+:[0-9]+:[0-9]+:", samples)]
mem <- do.call(rbind, lapply(strsplit(sub("^:([0-9]+:[0-9]+:[0-9]+:[0-9]+):.*", "\\1", samples), ":"), as.numeric))
heap <- (mem[, 1L] + mem[, 2L]) * 8
stacks <- lapply(sub("^:[0-9]+:[0-9]+:[0-9]+:[0-9]+:", "", samples), function(s) gsub('"', "", strsplit(trimws(s), " ")[[1L]]))
candidates <- c("svd", "La.svd", "cancor", "ruv_cancorplot", "trxCancor", "trxGridMember", "trxPrelimOrder", "estimateDisp", "glmQLFit", "glmQLFTest",
    "RUVs", "calcNormFactors", "filterByExpr", "cpm", "trxRuvGrid", "trxRuvs", "trxDaFit", "trxDaApply", "RUVIII", "runTrxReferenceRuvs",
    "runTrxReferenceDaAnalysis", "runTrxReferenceExpressionFilter", "runTrxReferenceRuviiiExport", "runTrxImportProcessing", "trxImportPreview")
out <- c(grep("^(##|trx_)", readLines(args[[2L]]), value = TRUE),
    sprintf("## (b) summarised offline from the stored profile: %d samples at %.2f s = %.0f s of %s", length(samples), interval, length(samples) * interval, basename(args[[1L]])),
    sprintf("%-34s %9s %8s %16s", "call", "seconds", "share %", "max R heap GiB"))
rows <- lapply(candidates, function(call) {
    on <- vapply(stacks, function(s) call %in% s, logical(1))
    if (!any(on)) return(NULL)
    data.frame(call = call, seconds = sum(on) * interval, share = 100 * mean(on), heap = max(heap[on]) / 1024^3)
})
rows <- do.call(rbind, rows); rows <- rows[order(-rows$seconds), ]
out <- c(out, sprintf("%-34s %9.1f %8.1f %16.2f", rows$call, rows$seconds, rows$share, rows$heap),
    sprintf("%-34s %9s %8s %16.2f", "(whole run)", "", "", max(heap) / 1024^3))
writeLines(out, args[[3L]]); cat(out, sep = "\n")
