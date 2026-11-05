# One assay at the large size (200 samples x 18,000 features, simulated): the time of each part of the omics RUV sweep.
# Usage: Rscript ruv-part-times.R <sweep|nipals|ruv183|ruv912>. For shares only: the parts may run side by side.
part <- commandArgs(TRUE)[[1]]
suppressMessages(library(MultiScholaR))
set.seed(20261105)
n <- 200L; p <- 18000L
group <- rep(c("A", "B"), length.out = n); sid <- sprintf("S%03d", seq_len(n))
mat <- matrix(rnorm(n * p, 15, 2), p, n, dimnames = list(sprintf("F%05d", seq_len(p)), sid)) + rep(rnorm(n, 0, 0.3), each = p)
design <- data.frame(Sample_ID = sid, group = group, row.names = sid)
Y <- t(mat); M <- model.matrix(~ 0 + group); rownames(M) <- sid
quiet <- function(x) suppressWarnings(suppressMessages(x))
stamp <- function(label, x) { t <- system.time(x)[["elapsed"]]; cat(sprintf("%s %.1f s\n", label, t)); invisible(t) }
if (part == "sweep") {
  a <- stamp("anova pct5", ctl <- quiet(MultiScholaR:::getNegCtrlProtAnovaHelper(mat, design_matrix = design, grouping_variable = "group", percentage_as_neg_ctrl = 5)))
  cat("controls", sum(ctl), "\n")
  c1 <- stamp("cancor pct5", quiet(ruv::ruv_cancorplot(Y, X = group, ctl = ctl)))
  cat(sprintf("sweep estimate, 20 percentages, one assay: %.0f s\n", 20 * (a + c1)))
} else if (part == "nipals") {
  Yn <- Y; Yn[sample(length(Yn), 0.05 * length(Yn))] <- NA
  stamp("nipals ncomp2 5pct missing", quiet(mixOmics::impute.nipals(Yn, ncomp = 2L)))
} else {
  k <- as.integer(sub("ruv", "", part)); ctl <- colnames(Y)[sample(p, k)]
  stamp(paste("RUVIII_C_Varying controls", k), quiet(RUVIIIC::RUVIII_C_Varying(k = 1L, Y = Y, M = M, toCorrect = colnames(Y), potentialControls = ctl)))
}
