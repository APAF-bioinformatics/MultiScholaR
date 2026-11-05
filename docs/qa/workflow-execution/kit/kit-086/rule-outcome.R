# rule-outcome.R <proof root> <stage>...: for each assay, the chosen percentage and k and the status of every percentage, from the stored sweep product. Synthetic corpora.
`%||%` <- function(a, b) if (is.null(a)) b else a
args <- commandArgs(TRUE); P <- args[[1L]]
for (st in args[-1L]) {
  d <- file.path(P, "kit/runtime", st, "candidate-1", "project")
  dirs <- list.dirs(d, recursive = TRUE); dirs <- dirs[grepl("ruv-optimisation", basename(dirs), fixed = TRUE)]
  if (!length(dirs)) { cat(st, "no sweep product\n"); next }
  for (dir in dirs) {
    o <- readRDS(file.path(dir, "payload.rds"))$value
    for (assay in names(o)) {
      a <- o[[assay]]; t <- a$optimization_results
      cat(st, assay, "best_percentage", a$best_percentage %||% NA, "best_k", a$best_k %||% NA, "success", isTRUE(a$success))
      if (is.data.frame(t)) { s <- table(t$status); cat(" percentages", nrow(t), "statuses:", paste(names(s), s, collapse = "; ")); if ("best_k" %in% names(t)) cat(" k_range", paste(range(t$best_k, na.rm = TRUE), collapse = "-")) }
      cat("\n")
    }
  }
}
