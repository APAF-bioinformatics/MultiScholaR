# KIT-114 (b)/(d), rows 10, 11 and 14 of parked-20261113.md: two mechanisms of the RUV-III correction path, on toy data.
# A: RUV-III-C on uncentred log2 values shifts each protein by a constant (the controls' mean level leaks into W); the same call on
#    protein-centred values with the means restored has no offset and the same contrasts.
# B: with the groups as the replicate structure, the correction removes k of the n - g residual dimensions the DA then
#    estimates its variance from. Usage: R_LIBS=<library> Rscript ruv_why.R
suppressMessages(library(RUVIIIC))
set.seed(1)
ruv <- function(Y, M, k, ctl) RUVIII_C_Varying(k = k, Y = Y, M = M, toCorrect = colnames(Y), potentialControls = ctl, progress = FALSE)
# ---- A. Row 14: where do per-protein offsets come from? Technical replicates, level ~22, one unwanted factor.
sim <- function(units, reps, p, nctl) {
  n <- units * reps; unit <- rep(seq_len(units), each = reps)
  mu <- rnorm(p, 22, 0.7); bio <- matrix(rnorm(units * p, 0, 0.3), units, p)[unit, ]
  w <- rnorm(n); a <- rnorm(p, 0, 0.3)
  Y <- matrix(mu, n, p, byrow = TRUE) + bio + outer(w, a) + matrix(rnorm(n * p, 0, 0.15), n, p)
  dimnames(Y) <- list(paste0("s", seq_len(n)), paste0("p", seq_len(p)))
  list(Y = Y, M = model.matrix(~ factor(unit) - 1), ctl = colnames(Y)[seq_len(nctl)])
}
cat("A. offsets (technical replicates as the replicate structure)\n")
cat(sprintf("%-28s %8s %10s %12s %12s %12s %14s\n", "case", "k", "controls", "med|offset|", "max|offset|", "sd ratio", "contrast diff"))
for (case in list(c(40, 3, 400, 200, 2), c(40, 3, 400, 50, 2), c(40, 3, 400, 50, 20), c(40, 3, 400, 50, 44), c(99, 3, 400, 50, 44))) {
  d <- sim(case[1], case[2], case[3], case[4]); k <- case[5]
  out <- ruv(d$Y, d$M, k, d$ctl)
  offset <- colMeans(out) - colMeans(d$Y)
  # the same call on protein-centred input, means restored afterwards
  mu <- colMeans(d$Y); cen <- sweep(d$Y, 2, mu)
  out2 <- sweep(ruv(cen, d$M, k, d$ctl), 2, mu, "+")
  g <- rep(c(-1, 1), length.out = nrow(d$Y)) # any contrast across samples
  contrast <- max(abs(crossprod(g - mean(g), out) - crossprod(g - mean(g), out2)))
  cat(sprintf("%-28s %8d %10d %12.2f %12.2f %12.2f %14.2e   centred: max|offset| %.2e\n",
    sprintf("%d samples, %d units", nrow(d$Y), case[1]), k, case[4], median(abs(offset)), max(abs(offset)),
    median(apply(out, 2, sd) / apply(d$Y, 2, sd)), contrast, max(abs(colMeans(out2) - mu))))
}
# ---- B. Row 11: groups as the replicate structure, pure null (no group effect), what is left of the residual variance?
cat("\nB. residual variance kept after correction, groups as replicates, no true effect\n")
cat(sprintf("%4s %3s %9s %22s %18s %12s\n", "n", "k", "controls", "variance kept (median)", "(n-g-k)/(n-g)", "t inflation"))
for (case in list(c(6, 3, 44), c(6, 3, 9), c(6, 3, 5), c(12, 3, 44), c(12, 3, 12), c(12, 3, 6), c(12, 1, 44), c(24, 3, 44))) {
  n <- case[1]; k <- case[2]; nctl <- case[3]; p <- 1000; kept <- numeric()
  for (r in 1:10) {
    grp <- rep(1:2, each = n / 2); w <- rnorm(n); a <- rnorm(p, 0, 0.3)
    Y <- matrix(rnorm(p, 22, 0.7), n, p, byrow = TRUE) + outer(w, a) + matrix(rnorm(n * p, 0, 0.3), n, p)
    dimnames(Y) <- list(paste0("s", 1:n), paste0("p", 1:p)); M <- model.matrix(~ factor(grp) - 1)
    out <- ruv(Y, M, k, colnames(Y)[1:nctl])
    rv <- function(X) colSums((X - M %*% solve(crossprod(M), crossprod(M, X)))^2) / (n - 2)
    # against the variance with the one real factor removed exactly (the best any method could do)
    truth <- Y - outer(w, a)
    kept <- c(kept, median(rv(out) / rv(truth)))
  }
  cat(sprintf("%4d %3d %9d %22.2f %18.2f %12.2f\n", n, k, nctl, mean(kept), (n - 2 - k) / (n - 2), 1 / sqrt(mean(kept))))
}
