# KIT-190 K190-D01: the declared metabolomics identity check, run late (audit "The decision", new finding). Compares,
# field by field, the metabolomics rows of `effect.R bound metab-check` (seeds 1 to 5 of every deciding cell, F+B and D+B,
# and their unbounded F and D where run) with the lipidomics rows of the counted stage 2 (`bound-s2/*/journeys.tsv`).
# Decides only whether "lipidomics and metabolomics" may be written beside the result. Base R only.
# Usage: Rscript metab-identity.R <metab journeys.tsv> <out dir>
args <- commandArgs(TRUE)
here <- "docs/qa/workflow-execution/kit/kit-190"
metab <- utils::read.delim(args[[1L]], stringsAsFactors = FALSE, colClasses = "character")
lipid <- do.call(rbind, lapply(list.files(file.path(here, "bound-s2"), "^journeys\\.tsv$", recursive = TRUE, full.names = TRUE),
  utils::read.delim, stringsAsFactors = FALSE, colClasses = "character"))
fields <- c("factor_group_r", "pct_u", "k_u", "percentage", "k", "n_controls", "dropped", "oracle_ok", "k_ok", "pct_ok", "calls",
  "any_call", "fdp", "power", "abs_lfc_error", "dims_below_k", "dims_at_k", "dims_above_k", "dims_below_kref", "dims_at_kref",
  "dims_above_kref", "removed_by_ruv", "separation")
key <- function(x) paste(x$cell, x$r, x$arm, x$mode)
metab <- metab[metab$cell != "E1", , drop = FALSE]
lipid <- lipid[key(lipid) %in% key(metab), , drop = FALSE]
rows <- lapply(seq_len(nrow(metab)), function(i) {
  m <- metab[i, ]
  l <- lipid[key(lipid) == key(m), , drop = FALSE]
  if (nrow(l) != 1L) return(data.frame(cell = m$cell, r = m$r, arm = m$arm, mode = m$mode, matched = FALSE, differing = NA))
  diff <- fields[vapply(fields, function(f) !identical(m[[f]], l[[f]]), logical(1L))]
  data.frame(cell = m$cell, r = m$r, arm = m$arm, mode = m$mode, matched = TRUE, differing = paste(diff, collapse = ","))
})
out <- do.call(rbind, rows)
dir.create(args[[2L]], recursive = TRUE, showWarnings = FALSE)
utils::write.table(out, file.path(args[[2L]], "metab-identity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
identical_rows <- sum(out$matched & !nzchar(out$differing))
cat(sprintf("metab identity: %d of %d rows identical in all %d fields; unmatched %d\n", identical_rows, nrow(out), length(fields),
  sum(!out$matched)))
