# KIT-099: the worker operations of a journey member, in order, with their seconds (KIT-084's worker-times.R): each worker
# request's directory, from the request file's time to its last file's time. Usage: Rscript worker-times.R <proof root> <stage>...
P <- commandArgs(TRUE)[1]
for (st in commandArgs(TRUE)[-1]) {
  d <- file.path(P, "kit/runtime", st, "candidate-1", "project")
  rq <- list.files(d, pattern = "^request\\.rds$", recursive = TRUE, full.names = TRUE, all.files = TRUE)
  rows <- lapply(rq, function(f) {
    r <- tryCatch(readRDS(f), error = function(e) NULL)
    fs <- list.files(dirname(f), recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
    journal <- file.path(dirname(f), "journal.json")
    status <- if (file.exists(journal)) tryCatch(jsonlite::read_json(journal)$status, error = function(e) NA) else NA
    data.frame(op = as.character(r$operation %||% NA)[1], start = file.mtime(f),
      secs = as.numeric(difftime(max(file.mtime(fs)), file.mtime(f), units = "secs")), status = as.character(status %||% NA))
  })
  x <- do.call(rbind, rows); x <- x[order(x$start), ]
  cat("##", st, nrow(x), "worker requests\n")
  for (i in seq_len(nrow(x))) cat(format(x$start[i], "%H:%M:%S"), sprintf("%8.1f", x$secs[i]), sprintf("%-24s", x$op[i]), x$status[i], "\n")
  agg <- aggregate(secs ~ op, x, sum); agg <- agg[order(-agg$secs), ]
  cat("# total seconds per operation\n"); for (i in seq_len(nrow(agg))) cat(sprintf("%-24s %8.1f\n", agg$op[i], agg$secs[i]))
}
