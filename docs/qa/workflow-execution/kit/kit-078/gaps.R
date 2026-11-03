# Where a member's sampler gaps fall against its worker requests. Usage: Rscript gaps.R <member dir>
m <- commandArgs(TRUE)[[1L]]
tr <- read.delim(file.path(m, "memory-trace.tsv"))
# wall = uptime + offset: the uptime stamps of the driver's phases against the trace's own clock are shared; the offset comes from the host.
offset <- as.numeric(Sys.time()) - as.numeric(strsplit(readLines("/proc/uptime", 1L), " ")[[1L]][[1L]])
d <- diff(tr$uptime); i <- which(d > 0.25)
cat("gaps above 0.25 s:\n"); for (k in i) cat(sprintf("  %.2f +%.2f s; anon %.2f GiB -> %.2f GiB\n", tr$uptime[k], d[k], tr$anon[k]/1024^3, tr$anon[k+1]/1024^3))
reqs <- list.files(m, "^request[.]rds$", recursive = TRUE, full.names = TRUE, all.files = TRUE)
reqs <- reqs[grepl("/[.]workflow-export/wfx-", reqs)]
rows <- lapply(reqs, function(p) { r <- readRDS(p); dir <- dirname(p)
  t0 <- file.mtime(p); w <- file.path(dir, "worker-result.rds"); a <- file.path(dir, "accepted-result.rds")
  data.frame(op = r$operation, requested = as.numeric(t0) - offset, result = if (file.exists(w)) as.numeric(file.mtime(w)) - offset else NA, accepted = if (file.exists(a)) as.numeric(file.mtime(a)) - offset else NA) })
rows <- do.call(rbind, rows); rows <- rows[order(rows$requested), ]
rows$worker_s <- round(rows$result - rows$requested, 1)
print(transform(rows, requested = round(requested, 1), result = round(result, 1), accepted = round(accepted, 1)), row.names = FALSE)
