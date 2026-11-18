# run-tests: the counts and text of R CMD check's ERRORs, WARNINGs and NOTEs, from 00check.log.
#
# Usage: Rscript --vanilla check-summary.R <00check.log> <out file> <R CMD check exit status> <R CMD build exit status>

arguments <- commandArgs(trailingOnly = TRUE)
log_path <- arguments[[1L]]
lines <- if (file.exists(log_path)) readLines(log_path, warn = FALSE) else character()
starts <- grep("^\\* ", lines)
blocks <- list()
for (i in seq_along(starts)) {
    end <- if (i < length(starts)) starts[[i + 1L]] - 1L else length(lines)
    block <- lines[starts[[i]]:end]
    end_at <- grep("^Status:", block)
    if (length(end_at)) block <- block[seq_len(end_at[[1L]] - 1L)]
    level <- regmatches(block[[1L]], regexpr("(ERROR|WARNING|NOTE)$", block[[1L]]))
    if (!length(level)) level <- if (length(block) > 1L && any(block[-1L] %in% c("ERROR", "WARNING", "NOTE"))) {
        block[-1L][block[-1L] %in% c("ERROR", "WARNING", "NOTE")][[1L]]
    } else ""
    if (nzchar(level)) blocks[[length(blocks) + 1L]] <- list(level = level, text = block)
}
levels <- vapply(blocks, `[[`, "", "level")
status <- grep("^Status:", lines, value = TRUE)
out <- c(
    sprintf("R CMD build exit status: %s; R CMD check exit status: %s", arguments[[4L]], arguments[[3L]]),
    sprintf("ERRORs %d, WARNINGs %d, NOTEs %d", sum(levels == "ERROR"), sum(levels == "WARNING"), sum(levels == "NOTE")),
    if (length(status)) status[[length(status)]] else "no Status line in 00check.log (check did not finish)",
    "")
for (block in blocks) out <- c(out, sprintf("[%s]", block$level), block$text, "")
writeLines(out, arguments[[2L]])
# the findings for the base-versus-head diff (assemble.R)
saveRDS(list(status = arguments[[3L]], levels = levels, texts = vapply(blocks, function(b) paste(b$text, collapse = "\n"), "")),
    file.path(dirname(arguments[[2L]]), "findings.rds"))
cat(out, sep = "\n")
