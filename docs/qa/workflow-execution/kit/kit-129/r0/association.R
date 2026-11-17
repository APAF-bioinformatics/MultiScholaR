# KIT-129 R0, part B: the statistic of the factor-tied-to-groups warning (ticket section 3), shared by the four
# KIT-124 simulators in their association mode. Evidence only: nothing under R/, tests/ or tools/ changes.
#
# The mode: K129_ASSOC=1 in the environment of a simulator (docs/qa/workflow-execution/kit/kit-124/sim/s*.R; same
# arguments as before, K124_ARMS unset). The simulator then regenerates each data set from its seed, runs the RUV step
# of the declared arms on the streams the KIT-124 study used, and records the statistic of each arm's estimated
# factors W. No arm's test is run. One exception by construction: arm 5 chooses its controls from a first-pass test
# of the groups, so that first pass is run where arm 5 reads it (S2: the route's no-RUV test; S4: the no-RUV kernel
# test; S1 and S3: their own first-pass fit, as in the study).
#
# Statistic (section 3): for each estimated factor W_j, the multiple correlation with the groups, the square root of
# R^2 of lm(W_j ~ group) on the samples the test uses. The largest over the factors is `group_association` (the
# package's `ruv_group_association`); every factor's value is kept in `factor_associations` (";"-separated, W_1 first).
#
# Arms (section 3: "arms 3 or 4, 5, 7 and 8, per simulator"). Arm 3 and arm 4 share the sweep's controls and k, so
# their W is one W per simulator; it is recorded under the arm whose code computes a W:
#   S1  4 (arm 3 subtracts and keeps no W: "W is computed for this purpose only, with the same controls and k")
#   S2  4 (the route's covariates mode; arm 3, correction mode, calls computeRuvCovariates with the same arguments)
#   S3  3 (today's design keeps W; S3 has no arm 4)
#   S4  4 (as S1)
ASSOC_ARMS <- list(S1 = c("4", "5", "7", "8"), S2 = c("4", "5", "7", "8"), S3 = c("3", "5", "7", "8"),
    S4 = c("4", "5", "7", "8"))

# One row per cell, run and arm (and assay in S4, pos and neg; "all" elsewhere).
#   status  ok         the arm gave a W; the statistic is recorded
#           no_factor  the arm is no RUV on this data set (E = 0): no W, no statistic
#           refused    a designed refusal (as in the study): no W, no statistic
#           crashed    any other error
#   k                 columns of W (0 without one)
#   true_rho          the cell's correlation of the unwanted factor with the group (NA in a cell without a factor)
#   true_association  the same statistic on the true factor that carries the cell's correlation (exactly true_rho in
#                     S1, S3 and S4, where it is fixed by construction; in S2 it is rho in expectation only)
#   drift_association S4 only: the statistic on the true drift along the injection order
ASSOC_COLUMNS <- c("simulator", "cell", "cell_class", "run", "seed", "arm", "assay", "status", "k", "n_samples",
    "group_association", "factor_associations", "true_rho", "true_association", "drift_association", "seconds")

# The statistic per column of W (samples in rows, in the order of `group`).
factorAssociations <- function(W, group) {
    W <- as.matrix(W)
    stopifnot(nrow(W) == length(group), !anyNA(W))
    group <- factor(group)
    vapply(seq_len(ncol(W)), function(j) sqrt(summary(stats::lm(W[, j] ~ group))$r.squared), numeric(1L))
}

assocRow <- function(simulator, cell, cell_class, run, seed, arm, status, W = NULL, group = NULL, true_rho = NA_real_,
                     true_factor = NULL, drift = NULL, seconds = NA_real_, assay = "all") {
    values <- if (identical(status, "ok")) factorAssociations(W, group) else numeric()
    one <- function(v) if (is.null(v) || is.null(group)) NA_real_ else factorAssociations(v, group)[[1L]]
    data.frame(simulator = simulator, cell = cell, cell_class = cell_class, run = run, seed = seed,
        arm = as.character(arm), assay = assay, status = status, k = length(values),
        n_samples = if (is.null(group)) NA_integer_ else length(group),
        group_association = if (length(values)) signif(max(values), 8L) else NA_real_,
        factor_associations = if (length(values)) paste(signif(values, 8L), collapse = ";") else NA_character_,
        true_rho = true_rho, true_association = signif(one(true_factor), 8L), drift_association = signif(one(drift), 8L),
        seconds = round(seconds, 3), stringsAsFactors = FALSE)[ASSOC_COLUMNS]
}
