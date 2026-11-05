P <- commandArgs(TRUE)[1]
for (st in c("kit-journey-084c", "kit-journey-084d")) {
  d <- file.path(P, "kit/runtime", st, "candidate-1", "project")
  rq <- list.files(d, pattern = "^request\\.rds$", recursive = TRUE, full.names = TRUE, all.files = TRUE)
  rows <- lapply(rq, function(f) {
    r <- tryCatch(readRDS(f), error = function(e) NULL)
    dir <- dirname(f)
    fs <- list.files(dir, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
    m <- file.mtime(fs)
    op <- as.character(r$operation %||% r$operation_id %||% r$action$operation %||% r$declaration %||% paste(names(r), collapse = ","))[1]
    data.frame(op = op, start = file.mtime(f), secs = as.numeric(difftime(max(m), file.mtime(f), units = "secs")))
  })
  x <- do.call(rbind, rows); x <- x[order(x$start), ]
  cat("##", st, nrow(x), "\n")
  for (i in seq_len(nrow(x))) cat(format(x$start[i], "%H:%M:%S"), sprintf("%7.1f", x$secs[i]), x$op[i], "\n")
}
