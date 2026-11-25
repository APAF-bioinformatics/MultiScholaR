# KIT-190: the per-cell confidence bounds the declaration promises (section 5) and decide.R did not write.
# COMPUTED AFTER THE DECISION; DESCRIPTIVE; NON-SIMULTANEOUS; DECIDES NOTHING. Forms as the audit names them ("The
# decision", (4)): effect cells, the one-sided 95 % paired-t upper bound mean(d) + t(0.95, n - 1) sd(d) / sqrt(n),
# d = FDP(D+B) - FDP(F+B); null cells, the one-sided 95 % upper bound on Delta = p10 - p01 (any call, D+B only minus F+B
# only) found by inverting Tango's score test: the Delta0 at which Z(Delta0) = -z(0.95). Base R only.
# Usage: Rscript bounds.R <journeys.tsv> ... -- <out dir>
args <- commandArgs(TRUE)
sep <- which(args == "--")
rows <- do.call(rbind, lapply(args[seq_len(sep - 1L)], utils::read.delim, stringsAsFactors = FALSE))
out <- args[[sep + 1L]]
z <- stats::qnorm(0.95)
tangoZ <- function(b, c, n, m) {
  score <- function(q) b / (q + m) + (if (c > 0) c / q else 0) - 2 * (n - b - c) / (1 - 2 * q - m)
  lo <- max(0, -m) + 1e-12
  hi <- (1 - m) / 2 - 1e-12
  q <- if (score(hi) > 0) hi else if (score(lo) < 0) lo else stats::uniroot(score, c(lo, hi), tol = 1e-13)$root
  (b - c - n * m) / sqrt(n * (2 * q + m * (1 - m)))
}
pairs <- function(cell) {
  f <- rows[rows$cell == cell & rows$arm == "F+B", , drop = FALSE]
  d <- rows[rows$cell == cell & rows$arm == "D+B", , drop = FALSE]
  r <- intersect(f$r, d$r)
  list(f = f[match(r, f$r), ], d = d[match(r, d$r), ])
}
res <- list()
for (cell in c("N-r", "N-s", "N-w")) {
  p <- pairs(cell)
  dc <- as.logical(p$d$any_call)
  fc <- as.logical(p$f$any_call)
  n <- length(dc)
  b <- sum(dc & !fc)
  c <- sum(fc & !dc)
  est <- (b - c) / n
  upper <- stats::uniroot(function(m) tangoZ(b, c, n, m) + z, c(est + 1e-9, 0.99), tol = 1e-12)$root
  res[[cell]] <- data.frame(cell = cell, endpoint = "any_call difference", n = n, estimate = est, upper_bound_95 = upper,
    margin = 0.10, form = "inverted Tango score test")
}
for (cell in c("E-r", "E-s", "E-p", "E-g", "E-w")) {
  p <- pairs(cell)
  d <- as.numeric(p$d$fdp) - as.numeric(p$f$fdp)
  n <- length(d)
  res[[cell]] <- data.frame(cell = cell, endpoint = "mean FDP difference", n = n, estimate = mean(d),
    upper_bound_95 = mean(d) + stats::qt(0.95, n - 1) * stats::sd(d) / sqrt(n), margin = 0.05, form = "paired t")
}
res <- do.call(rbind, res)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
con <- file(file.path(out, "bounds.tsv"), "w")
writeLines("# computed after the decision; descriptive; non-simultaneous one-sided 95 % bounds; decides nothing", con)
utils::write.table(res, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)
print(res)
