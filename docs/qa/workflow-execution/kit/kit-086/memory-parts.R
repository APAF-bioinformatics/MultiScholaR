# One assay at the large size (200 samples x 18,285 features, simulated, the declared missing pattern): the peak memory of each step
# that behaves differently with missing values. Usage: Rscript memory-parts.R <part> <missing 0|1>. Each part runs in its own process.
args <- commandArgs(TRUE); part <- args[[1]]; with_missing <- args[[2]] == "1"
suppressMessages(library(MultiScholaR))
set.seed(20261105)
n <- 200L; p <- 18285L; sid <- sprintf("S%03d", seq_len(n)); group <- rep(c("A", "B"), length.out = n)
base <- rnorm(p, 15, 1.8)
mat <- outer(base, rep(0, n), "+") + matrix(rep(rnorm(n, 0, 0.3), each = p), p, n) + matrix(rnorm(n * p, 0, 0.3), p, n)
dimnames(mat) <- list(sprintf("F%05d", seq_len(p)), sid)
if (with_missing) for (f in order(base)[seq_len(round(p * 0.3))]) mat[f, sample.int(n, round(n * 0.15))] <- NA_real_
Y <- t(mat); M <- model.matrix(~ 0 + group); rownames(M) <- sid
ctl <- colnames(Y)[order(-base)[seq_len(900L)]]
quiet <- function(x) suppressWarnings(suppressMessages(x))
t <- system.time(switch(part,
  base = NULL,
  nipals = if (with_missing) quiet(mixOmics::impute.nipals(Y, ncomp = 2L)),
  pca = quiet(mixOmics::pca(Y, ncomp = 2L)),
  ruv = quiet(RUVIIIC::RUVIII_C_Varying(k = 1L, Y = Y, M = M, toCorrect = colnames(Y), potentialControls = ctl)),
  cor = quiet(stats::cor(mat, use = "pairwise.complete.obs")),
  limma = quiet(limma::eBayes(limma::contrasts.fit(limma::lmFit(mat, stats::model.matrix(~ 0 + group)), c(1, -1))))
))[["elapsed"]]
gc_max <- sum(gc()[, 6L])
cat(sprintf("PART %s missing %s elapsed %.1f s gc_max_mb %.0f\n", part, with_missing, t, gc_max))
