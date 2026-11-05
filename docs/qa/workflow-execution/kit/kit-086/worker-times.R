# worker-times.R <proof root> <stage>...: each worker step's start and duration, from its request file and the last file written in its directory.
`%||%` <- function(a, b) if (is.null(a)) b else a
args <- commandArgs(TRUE); P <- args[[1L]]
for (st in args[-1L]) {
  d <- file.path(P, "kit/runtime", st, "candidate-1", "project")
  rq <- list.files(d, pattern = "^request\\.rds$", recursive = TRUE, full.names = TRUE, all.files = TRUE)
  rows <- lapply(rq, function(f) {
    r <- tryCatch(readRDS(f), error = function(e) NULL)
    fs <- list.files(dirname(f), recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
    op <- as.character(r$operation %||% r$operation_id %||% r$action$operation %||% r$declaration %||% paste(names(r), collapse = ","))[1]
    data.frame(op = op, start = file.mtime(f), secs = as.numeric(difftime(max(file.mtime(fs)), file.mtime(f), units = "secs")))
  })
  x <- do.call(rbind, rows); x <- x[order(x$start), ]
  cat("##", st, nrow(x), "steps\n")
  for (i in seq_len(nrow(x))) cat(format(x$start[i], "%H:%M:%S"), sprintf("%7.1f", x$secs[i]), x$op[i], "\n")
}
