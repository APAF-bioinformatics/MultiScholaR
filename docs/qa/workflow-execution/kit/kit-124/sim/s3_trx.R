# KIT-124 simulator S3 (transcriptomics), ticket revision 1.2. Evidence only: nothing under R/, tests/ or tools/ changes.
#
# Usage (from the repository root, with R_LIBS holding the pinned library and kit/scratch-k124-lib for sva):
#   Rscript docs/qa/workflow-execution/kit/kit-124/sim/s3_trx.R <out dir> <runs per core cell> <runs per stress cell> <cores>
#       [cells = all | comma list, e.g. C0,S7] [first run number = 1]
# Writes <out dir>/s3-runs.tsv (one row per cell, run, arm and subset), s3-runs.meta.txt, s3-events.tsv (every refusal
# and crash with its condition class and message), s3-timing.tsv (shared preparation and total seconds per run) and
# cells/s3-<cell>.rds (the restart unit: a finished cell with the same run numbers is not run again).
#
# ---------------------------------------------------------------------------------------------------------------------
# PACKAGE FUNCTIONS CALLED (MultiScholaR through pkgload::load_all("."); arguments as a default GUI run passes them)
#   trxExpressionFilter(counts, group)           edgeR::filterByExpr defaults on TMM sizes, then the zero-variance mask
#   trxTmmLogCpm(filtered_counts, group)         TMM on the filtered counts, cpm(log = TRUE)
#   trxFilteredDge(filtered_counts, group)       the DGEList of the first pass (the calls of trxPrelimOrder, repeated here
#                                                to keep the p-values: estimateDisp, glmQLFit, glmQLFTest(coef = 2:G))
#   trxRuvGrid(filtered_counts, log_cpm_tmm, design)                 params = list(): percentage 5..30 by 5, k penalty 0.05
#   trxRuvs(filtered_counts, log_cpm_tmm, design, grid$control_genes, grid$best$best_k)
#   trxDaContrast(ruvs_counts, pdata, "B_vs_A=groupB-groupA")        arm 3: ~ 0 + group + W_1..W_k on the original counts
#   trxDaFormula(pdata); trxDaFit(counts, pdata, formula); trxDaApply(fit, "B_vs_A=groupB-groupA")
#                                                every other arm: the same fit and test halves (TMM, estimateDisp(robust =
#                                                TRUE), glmQLFit(robust = TRUE), glmQLFTest, topTags) with that arm's
#                                                covariates; tests/testthat/test-trx-kernels.R (KIT-124 W1) pins
#                                                trxDaContrast == trxDaFit == the direct edgeR calls (R-16).
#   Significance rule of this omic: edgeR's Benjamini-Hochberg FDR < 0.05 (da_q_val_thresh default 0.05; no fold-change
#   threshold, no TREAT; R/utils_workflow_trx_declarations.R trxDaFormat). A "call" below is FDR < 0.05.
#   The package's DA cannot run without W columns (its formula is "~ 0 + group + W..."); arm 1 passes "~ 0 + group" to
#   trxDaFit, which is the same fit without covariates.
#
# ---------------------------------------------------------------------------------------------------------------------
# DATA MODEL (per cell and run, generated once from set.seed(seed) and shared by all arms)
#   Genes G = 5,000. Samples j in groups A, B (and C in cell C7); the tested contrast is always B - A.
#   Baseline expected count lambda_g: log2(lambda_g) ~ N(6.5, 2^2) (log-normal, median 91 counts; a log-normal is the
#     usual approximation to the bulk RNA-seq distribution of mean counts after low-count filtering).
#   Dispersion phi_g = (0.2 + 1/sqrt(lambda_g))^2 * 40 / chisq_40: the trend and gene-wise scatter of the voom paper's
#     simulation (Law, Chen, Shi, Smyth 2014, Genome Biology 15:R29: BCV trend 0.2 + 1/sqrt(expected count), scaled
#     inverse chi-square on 40 df). BCV about 0.3 at the median gene, within edgeR's user-guide range for bulk data (0.1
#     inbred to 0.4 human). The trend is anchored at the centre depth so that sequencing depth does not change biology.
#   Library depth L_j ~ logNormal(meanlog = -0.25^2/2, sdlog = 0.25) (mean 1), divided by lib_div.
#   Group effects: round(pi * G) genes with beta = +1 or -1 log2 in group B (the share "up" exact); in C7 an independent
#     set of the same size and sign rule in group C.
#   Unwanted factors: nf columns W (samples x nf), each with sample mean 0 and sample variance 1, exactly orthogonal to
#     each other. With u the standardised score (B minus A), factor 1 is rho * u + sqrt(1 - rho^2) * z_1, z a normal
#     draw residualised on the group model and orthonormalised: the sample correlation of factor 1 with group is exactly
#     rho. Loadings alpha_gf = s_f * e_gf, e ~ N(0, 1) (S10: e = 0.7 * z(log2 lambda) + sqrt(1 - 0.49) * N(0, 1)).
#     Strength: s_f^2 = ratio_f * median_g[(phi_g + 1/mu_g) / ln(2)^2], the median residual variance of log2 counts (delta
#     method; mu_g = lambda_g / lib_div; phi_g + 0.01 in C6). ratio = 3 for one factor; 3, 1.5, 0.75 for three.
#   Counts: y_gj ~ NB(mean = L_j * lambda_g * 2^(beta_g,group(j) + sum_f alpha_gf W_jf), dispersion phi_g).
#   Technical duplicates (C6): 4 biological units per group, 2 libraries each. Unit level: a Gamma(shape 1/phi_g, mean 1)
#     multiplier (biological dispersion phi_g, BCV about 0.3); library level: NB with technical dispersion 0.01 (BCV 0.1,
#     library preparation replicates). Depth and the unwanted factor vary per library. The default path sees 16 runs.
#   Heavy tails (S7): every count independently with probability 0.05 multiplied by Uniform(1.5, 10) and rounded: the
#     "random outliers" of Soneson and Delorenzi 2013 (BMC Bioinformatics 14:91). The negative binomial is already
#     heavy-tailed on the log scale; a t-distributed log mean would only change the dispersion distribution, whereas
#     isolated gene-and-sample outliers are what edgeR's robust options exist for.
#   Known batch (S8): the unwanted factor is a step: batch 1/2 = 4/2 in group A and 2/4 in group B (correlation with
#     group 1/3 by design), standardised, strength 3. Arm 2 has it as the true factor, arm 17 as `batch` in the formula
#     (the same column space by construction).
#   Rows: gene-level parameters are built in blocks (effects first) and permuted jointly before counts are drawn; gene
#     ids are assigned after the permutation (R-19).
#
# CELLS (index = position in seed rule; n = runs per group; pi = share of effects; rho; nf; other)
#   idx id  class  n   groups pi   rho nf  other
#    0  C0  core   6   2      0.2  0   1   centre
#    1  C1  core   12  2      0.2  0   1
#    2  C2  core   6   2      0.1  0   1
#    3  C3  core   6   2      0.2  0   0   no unwanted factor
#    4  C4  core   6   2      0    0   1   global null
#    5  C5  core   6   2      0.2  0.3 1
#    6  C6  core   8   2      0.2  0   1   technical duplicates: 4 units per group, 2 libraries each
#    7  C7  core   6   3      0.2  0   1   three groups
#    8  S1  stress 3   2      0.2  0   1
#    9  S2  stress 6   2      0.4  0   1
#   10  S3  stress 6   2      0.2  0.6 1
#   11  S4  stress 6   2      0.2  0   3
#   12  S5  stress 6   2      0.2  0   1   low counts: lib_div = 10 (replaces the ticket's "0% missing", not applicable)
#   13  S6  -                              not run: "50% missing" has no counterpart for counts
#   14  S7  stress 6   2      0.2  0   1   heavy tails (outliers)
#   15  S8  stress 6   2      0.2  1/3 1   known unbalanced batch
#   16  S9  stress 6   2      0.2  0   1   80% of effects up
#   17  S10 stress 6   2      0.2  0   1   loadings tied to abundance (correlation 0.7)
#   18  S11 stress 3   2      0.2  0.3 1
#   19  S12 stress 3   2      0.2  0   3
#   20  S13 stress 3   2      0.2  0   0
#   21  S14 stress 6   2      0.2  0   0   no factor, low counts: lib_div = 10 (replaces "no factor with 50% missing")
#   Random-allocation cells (KIT-129 R0 part B supplement, decision M37; never part of "all"; run only when named):
#   22  C0R core   6   2      0.2  0   1   C0 with the factor allocated at random: z is the same normal draw, standardised
#                                          (sample mean 0, variance 1) but NOT residualised on the group model, so its
#                                          correlation with group is 0 in expectation only (as S2 does it)
#   23  C4R core   6   2      0    0   1   C4, factor allocated at random as C0R
#   24  S1R stress 3   2      0.2  0   1   S1, factor allocated at random as C0R
#   In the study mode (K129_ASSOC unset) a random-allocation cell also writes s3-assoc.tsv: one row per run and arm run,
#   the factor-tied-to-groups statistic of the W the arm tested with (columns of kit-129/r0/association.R; true_association
#   is the realised association of the true factor). Every other cell writes no such table and is computed as before.
#
# ARMS (label; random stream number in brackets where it differs)
#   1   no RUV: ~ 0 + group.
#   2   the true factor(s) in the model as W_1..W_nf (in C3, S13, S14 there is none: equal to arm 1, k = 0).
#   3   today's design: trxRuvGrid -> trxRuvs (group as replicate sets, the grid's controls and k) -> trxDaContrast.
#   5   RUVSeq::RUVg, controls = all genes but the most significant quarter of the first pass, k = today's grid's.
#   6   RUVSeq::RUVg, controls = first-pass BH q > 0.5, k = today's grid's.
#   7   RUVSeq::RUVs, all genes as controls, group as replicate sets (the vignette's RUVs usage), k = E.
#   8   RUVSeq::RUVg, 50 oracle null genes (no effect in any group), k = E.
#   9   RUVSeq::RUVg, degraded set: 40 nulls and 10 true effects of the tested contrast (50 nulls in C4), k = E.
#   10  arm 8's fit, the test on its normalizedCounts with ~ 0 + group (correct, then test).
#   11  arm 7 with k = 1.
#   R   RUVSeq::RUVr [18]: deviance residuals of a first-pass glmFit of ~ group (the vignette's), all genes, k = E.
#   14  sva::svaseq, its own number of factors; 14E [19] svaseq with n.sv = E. Surrogate variables as W in the model.
#   15  C6 only: RUVs with the technical-duplicate sets, all genes as controls, k = E, test on normalizedCounts.
#   16  C6 only: the same RUVs fit with W in the model.
#   17  S8 only: ~ 0 + group + batch.
#   22  on request (K124_ARMS): arm 8 (RUVg, the same 50 oracle genes: drawn again from arm 8's stream) with k from
#       today's separation rule on those genes: trxCancor(t(log_cpm_tmm), the group columns, the controls), then
#       trxFindBestK up to the grid's limit (3 under 15 runs, else 4), as a grid member computes it.
#   23  on request (K124_ARMS): arm 5 (RUVg, all genes but the most significant quarter of the first pass) with k = E.
#   24  on request (K124_ARMS): arm 5 with k from today's separation rule on arm 5's own controls (arm 22's rule). Arm 5
#       takes the grid's k, which the grid computed on its own winning control set; arm 24 (user, 2026-11-18) differs from
#       arm 5 in k only, and does not read the grid.
#   Arms 22 and 23 (user, 2026-11-18; 17 to 21 are taken as arm or stream numbers in one simulator or another) cross the
#   two k rules over the two control sets: 22 differs from 8, and 23 from 5, in k only. Arm 22's k is never 0 (nor 24's); arm 23
#   with E = 0 is no RUV, as arms 7 to 10.
#   K124_ARMS (environment; unset = the arms above, 22 to 24 not run): a comma list of the arm labels to run, of those
#   the cell has and 22, 23, 24. Arm 1 must be in it (refusals and E = 0 carry its result). An arm that reads another arm's
#   fit (5 and 6 the grid of 3; 10 the fit of 8; 16 of 15) is refused when that arm is left out. The data, E and every
#   stream are as without it, so an arm's rows are the same whichever other arms run.
#   Arms 5 to 11, R, 14, 15 and 16 estimate on upper-quartile normalised counts (EDASeq::betweenLaneNormalization(which =
#   "upper"), as the RUVSeq vignette does) and, 10 and 15 apart, test on the original filtered counts with W in the
#   model. Arm 3 passes the filtered counts as the package does. The first pass of arms 5 and 6 is the package's own
#   (the quasi-likelihood F-test of ~ group, not the vignette's likelihood-ratio test). The ticket's arms 4 (equal to
#   arm 3 here: today's design already puts W in the model), 12 and 13 (ruv::RUV4, RUVinv: log-scale methods) are not
#   run in S3.
#   E = sva::num.sv(log_cpm_tmm, model.matrix(~ group), method = "be", B = 20, seed = the run's seed). E = 0: the arm is
#   no RUV (status ok, k = 0, arm 1's result). The k column is the number of covariates in the fitted model.
#
# MEASURES per row (subset is always "full": counts have no missing values)
#   n_features after the expression filter; truth restricted to those genes; n_true = true effects of B - A among them;
#   n_calls = FDR < 0.05; null_p05, null_p01 = share of true nulls with p below 0.05, 0.01; fdp = false calls / calls
#   (0 with no calls); power = share of true effects called; bias = mean of sign(true) * (logFC - true) over true effects
#   (negative = attenuated); rmse = root mean square of (logFC - true); seconds = the arm's own elapsed time (the shared
#   preparation is in s3-timing.tsv). status: ok; refused (a condition of class trx_error or s3_refusal: a designed
#   refusal, scored as the run's no-RUV result, R-07); crashed (any other error; measures NA).
#   Seed rule: seed = 100000 + 10000 * cell index + run; data from set.seed(seed); arm a from set.seed(seed * 100 + a).
#
# K129_ASSOC (environment, KIT-129 R0 part B; unset = everything above): the association mode. K124_ARMS must be unset.
#   The data, the shared preparation (filter, TMM, the first pass, upper-quartile counts, E) and each arm's RUV step run
#   as above on the same streams for the arms of docs/qa/workflow-execution/kit/kit-129/r0/association.R (here 3, 5, 7,
#   8; S3 has no arm 4), and no arm is tested: an arm stops at its W and the row holds the statistic of the
#   factor-tied-to-groups warning (one row per cell, run and arm; columns in association.R). Written to s3-assoc.tsv;
#   events, timing and the sidecar as above.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L) stop("usage: Rscript s3_trx.R <out dir> <runs per core cell> <runs per stress cell> <cores> [cells = all] [first run = 1]")
out_dir <- args[[1L]]
runs_core <- as.integer(args[[2L]])
runs_stress <- as.integer(args[[3L]])
cores <- as.integer(args[[4L]])
cells_arg <- if (length(args) >= 5L) args[[5L]] else "all"
first_run <- if (length(args) >= 6L) as.integer(args[[6L]]) else 1L
arms_wanted <- Sys.getenv("K124_ARMS", "")
arms_wanted <- if (nzchar(arms_wanted)) strsplit(arms_wanted, ",", fixed = TRUE)[[1L]]
stopifnot(is.null(arms_wanted) || "1" %in% arms_wanted)
assoc_mode <- nzchar(Sys.getenv("K129_ASSOC", ""))
if (assoc_mode) {
    stopifnot(is.null(arms_wanted))
    source("docs/qa/workflow-execution/kit/kit-129/r0/association.R")
    arms_wanted <- ASSOC_ARMS$S3
}
started <- Sys.time()

suppressMessages(pkgload::load_all(".", quiet = TRUE))
for (package in c("edgeR", "RUVSeq", "EDASeq", "sva")) requireNamespace(package)

cell <- function(idx, id, class, n = 6L, groups = 2L, pi = 0.2, rho = 0, nf = 1L, tech = FALSE, lib_div = 1, outliers = FALSE,
                 batch = FALSE, up = 0.5, abund = 0, random = FALSE) {
    list(idx = idx, id = id, class = class, n = n, groups = groups, pi = pi, rho = rho, nf = nf, tech = tech, lib_div = lib_div,
        outliers = outliers, batch = batch, up = up, abund = abund, random = random)
}
CELLS <- list(
    cell(0L, "C0", "core"), cell(1L, "C1", "core", n = 12L), cell(2L, "C2", "core", pi = 0.1), cell(3L, "C3", "core", nf = 0L),
    cell(4L, "C4", "core", pi = 0), cell(5L, "C5", "core", rho = 0.3), cell(6L, "C6", "core", n = 8L, tech = TRUE),
    cell(7L, "C7", "core", groups = 3L),
    cell(8L, "S1", "stress", n = 3L), cell(9L, "S2", "stress", pi = 0.4), cell(10L, "S3", "stress", rho = 0.6),
    cell(11L, "S4", "stress", nf = 3L), cell(12L, "S5", "stress", lib_div = 10), cell(14L, "S7", "stress", outliers = TRUE),
    cell(15L, "S8", "stress", batch = TRUE), cell(16L, "S9", "stress", up = 0.8), cell(17L, "S10", "stress", abund = 0.7),
    cell(18L, "S11", "stress", n = 3L, rho = 0.3), cell(19L, "S12", "stress", n = 3L, nf = 3L),
    cell(20L, "S13", "stress", n = 3L, nf = 0L), cell(21L, "S14", "stress", nf = 0L, lib_div = 10),
    cell(22L, "C0R", "core", random = TRUE), cell(23L, "C4R", "core", pi = 0, random = TRUE),
    cell(24L, "S1R", "stress", n = 3L, random = TRUE))
names(CELLS) <- vapply(CELLS, `[[`, "", "id")

G <- 5000L
TECH_DISPERSION <- 0.01
RATIOS <- list(`1` = 3, `3` = c(3, 1.5, 0.75))
ENTRY <- "B_vs_A=groupB-groupA"
ARM_STREAM <- c(`1` = 1L, `2` = 2L, `3` = 3L, `5` = 5L, `6` = 6L, `7` = 7L, `8` = 8L, `9` = 9L, `10` = 10L, `11` = 11L, R = 18L,
    `14` = 14L, `14E` = 19L, `15` = 15L, `16` = 16L, `17` = 17L, `22` = 22L, `23` = 23L, `24` = 24L)
COLUMNS <- c("simulator", "cell", "cell_class", "run", "seed", "arm", "subset", "status", "k", "n_features", "n_true", "n_calls",
    "null_p05", "null_p01", "fdp", "power", "bias", "rmse", "seconds")

standardise <- function(v) (v - mean(v)) / stats::sd(v)

# One data set: counts (genes x runs, rows permuted), the design, the truth and the true factors.
generate <- function(cl, seed) {
    set.seed(seed)
    levels <- LETTERS[seq_len(cl$groups)]
    group <- factor(rep(levels, each = cl$n), levels = levels)
    n <- length(group)
    unit <- if (cl$tech) rep(seq_len(n / 2L), each = 2L) else seq_len(n)
    batch <- if (cl$batch) factor(c(1, 1, 1, 1, 2, 2, 1, 1, 2, 2, 2, 2)) else NULL
    log2_lambda <- stats::rnorm(G, 6.5, 2)
    lambda <- 2^log2_lambda
    phi <- (0.2 + 1 / sqrt(lambda))^2 * 40 / stats::rchisq(G, 40)
    n_eff <- round(cl$pi * G)
    signs <- rep(c(1, -1), c(round(cl$up * n_eff), n_eff - round(cl$up * n_eff)))
    beta <- matrix(0, G, cl$groups, dimnames = list(NULL, levels))
    beta[seq_len(n_eff), "B"] <- signs
    if (cl$groups == 3L) beta[sample.int(G, n_eff), "C"] <- signs
    depth <- stats::rlnorm(n, -0.25^2 / 2, 0.25) / cl$lib_div
    # the unwanted factors: sample mean 0, variance 1, factor 1's correlation with the group score exactly rho
    nf <- cl$nf
    W <- matrix(0, n, nf)
    if (cl$batch) {
        W[, 1L] <- standardise(as.numeric(batch == "2"))
    } else if (nf > 0L && cl$random) {
        W <- apply(matrix(stats::rnorm(n * nf), n, nf), 2L, standardise)
        W[, 1L] <- cl$rho * standardise((group == "B") - (group == "A")) + sqrt(1 - cl$rho^2) * W[, 1L]
    } else if (nf > 0L) {
        z <- qr.resid(qr(stats::model.matrix(~ group)), matrix(stats::rnorm(n * nf), n, nf))
        z <- qr.Q(qr(z)) * sqrt(n - 1)
        W <- z
        W[, 1L] <- cl$rho * standardise((group == "B") - (group == "A")) + sqrt(1 - cl$rho^2) * z[, 1L]
    }
    alpha <- matrix(0, G, nf)
    if (nf > 0L) {
        residual_variance <- stats::median((phi + (if (cl$tech) TECH_DISPERSION else 0) + cl$lib_div / lambda) / log(2)^2)
        e <- matrix(stats::rnorm(G * nf), G, nf)
        if (cl$abund > 0) e <- cl$abund * standardise(log2_lambda) + sqrt(1 - cl$abund^2) * e
        alpha <- sweep(e, 2L, sqrt(RATIOS[[as.character(nf)]] * residual_variance), `*`)
    }
    perm <- sample.int(G)
    lambda <- lambda[perm]; phi <- phi[perm]; beta <- beta[perm, , drop = FALSE]; alpha <- alpha[perm, , drop = FALSE]
    mu <- outer(lambda, depth) * 2^(beta[, as.integer(group), drop = FALSE] + alpha %*% t(W))
    if (cl$tech) {
        units <- unique(unit)
        bio <- matrix(stats::rgamma(G * length(units), shape = 1 / phi, scale = phi), G, length(units))
        counts <- stats::rnbinom(G * n, mu = mu * bio[, match(unit, units)], size = 1 / TECH_DISPERSION)
    } else {
        counts <- stats::rnbinom(G * n, mu = mu, size = 1 / phi)
    }
    counts <- matrix(as.double(counts), G, n)
    if (cl$outliers) {
        hit <- stats::runif(G * n) < 0.05
        counts[hit] <- round(counts[hit] * stats::runif(sum(hit), 1.5, 10))
    }
    runs <- sprintf("%s_%02d", group, stats::ave(seq_len(n), group, FUN = seq_along))
    dimnames(counts) <- list(sprintf("g%05d", seq_len(G)), runs)
    list(counts = counts, design = data.frame(Run = runs, group = group, stringsAsFactors = FALSE), unit = unit, batch = batch, W = W,
        truth = stats::setNames(beta[, "B"] - beta[, "A"], rownames(counts)),
        null_all = stats::setNames(rowSums(beta != 0) == 0, rownames(counts)))
}

refuse <- function(message) rlang::abort(message, class = "s3_refusal")

# The package's fit and test halves with an arm's covariates; the table in the order of the counts' rows.
daTest <- function(counts, design, W = NULL, formula = NULL) {
    pdata <- data.frame(design, row.names = design$Run)
    if (!is.null(W)) {
        W <- as.matrix(W)
        colnames(W) <- paste0("W_", seq_len(ncol(W)))
        pdata <- cbind(pdata, W)
    }
    if (is.null(formula)) formula <- if (is.null(W)) "~ 0 + group" else trxDaFormula(pdata)
    table <- trxDaApply(trxDaFit(counts, pdata, formula), ENTRY)$results[[1L]]$top_table
    list(table = table[rownames(counts), c("logFC", "PValue", "FDR")], k = if (is.null(W)) 0L else ncol(W))
}

measures <- function(result, truth) {
    t <- result$table
    effect <- truth != 0
    called <- t$FDR < 0.05
    error <- t$logFC[effect] - truth[effect]
    data.frame(k = result$k, n_features = length(truth), n_true = sum(effect), n_calls = sum(called),
        null_p05 = mean(t$PValue[!effect] < 0.05), null_p01 = mean(t$PValue[!effect] < 0.01),
        fdp = if (any(called)) sum(called & !effect) / sum(called) else 0,
        power = if (any(effect)) mean(called[effect]) else NA_real_,
        bias = if (any(effect)) mean(sign(truth[effect]) * error) else NA_real_,
        rmse = if (any(effect)) sqrt(mean(error^2)) else NA_real_)
}

NA_MEASURES <- data.frame(k = NA_integer_, n_features = NA_integer_, n_true = NA_integer_, n_calls = NA_integer_, null_p05 = NA_real_,
    null_p01 = NA_real_, fdp = NA_real_, power = NA_real_, bias = NA_real_, rmse = NA_real_)

# One (cell, run): every arm on the same data. Returns list(rows, events, timing).
oneRun <- function(cl, run) {
    seed <- 100000L + 10000L * cl$idx + run
    events <- list()
    rows <- list()
    side <- list()
    no_ruv <- NULL
    assocRecord <- function(arm, status, W = NULL, seconds = NA_real_) {
        rows[[length(rows) + 1L]] <<- assocRow("S3", cl$id, cl$class, run, seed, arm, status, W = W,
            group = if (is.null(prep$d)) NULL else prep$d$design$group, true_rho = if (cl$nf > 0L) cl$rho else NA_real_,
            true_factor = if (cl$nf > 0L && !is.null(prep$d)) prep$d$W[, 1L], seconds = seconds)
    }
    record <- function(arm, status, m, seconds) {
        if (assoc_mode) return(assocRecord(arm, status, seconds = seconds))
        rows[[length(rows) + 1L]] <<- cbind(data.frame(simulator = "S3", cell = cl$id, cell_class = cl$class, run = run, seed = seed,
            arm = arm, subset = "full", status = status, stringsAsFactors = FALSE), m, data.frame(seconds = seconds))[, COLUMNS]
    }
    event <- function(arm, status, condition) {
        events[[length(events) + 1L]] <<- data.frame(cell = cl$id, run = run, arm = arm, status = status,
            condition_class = class(condition)[[1L]], message = gsub("[\t\n]+", " ", conditionMessage(condition)), stringsAsFactors = FALSE)
    }
    # an arm: its own stream; a designed refusal is scored as the no-RUV result, any other error is a crash
    arm <- function(label, body) {
        if (!label %in% armsOf(cl)) return(NULL)
        t0 <- proc.time()[["elapsed"]]
        set.seed(seed * 100L + ARM_STREAM[[label]])
        result <- tryCatch(suppressWarnings(suppressMessages(body())), error = function(e) e)
        seconds <- proc.time()[["elapsed"]] - t0
        if (!assoc_mode && cl$random) {                 # the statistic of the W this arm tests with
            side_status <- if (inherits(result, "error")) {
                if (inherits(result, "trx_error") || inherits(result, "s3_refusal")) "refused" else "crashed"
            } else if (is.null(result$W)) "no_factor" else "ok"
            side[[length(side) + 1L]] <<- assocRow("S3", cl$id, cl$class, run, seed, label, side_status,
                W = if (side_status == "ok") result$W, group = prep$d$design$group,
                true_rho = if (cl$nf > 0L) cl$rho else NA_real_, true_factor = if (cl$nf > 0L) prep$d$W[, 1L], seconds = seconds)
        }
        if (assoc_mode && !inherits(result, "error")) {
            assocRecord(label, if (is.null(result$W)) "no_factor" else "ok", result$W, seconds)
            return(result)
        }
        if (inherits(result, "error")) {
            designed <- inherits(result, "trx_error") || inherits(result, "s3_refusal")
            event(label, if (designed) "refused" else "crashed", result)
            if (designed && !is.null(no_ruv)) {
                refused <- measures(no_ruv, truth)
                refused$k <- 0L
                record(label, "refused", refused, seconds)
            } else {
                record(label, if (designed) "refused" else "crashed", NA_MEASURES, seconds)
            }
            return(NULL)
        }
        record(label, "ok", measures(result, truth), seconds)
        result
    }
    noRuv <- function() list(table = no_ruv$table, k = 0L)

    t_start <- proc.time()[["elapsed"]]
    truth <- NULL
    prep <- tryCatch(suppressWarnings(suppressMessages({
        d <- generate(cl, seed)
        f <- trxExpressionFilter(d$counts, d$design$group)
        y <- f$filtered_counts
        n <- trxTmmLogCpm(y, d$design$group)
        model <- stats::model.matrix(~ group, data = d$design)
        first <- edgeR::glmQLFit(edgeR::estimateDisp(trxFilteredDge(y, d$design$group), model), model)
        first_p <- edgeR::topTags(edgeR::glmQLFTest(first, coef = 2:ncol(model)), n = nrow(y), sort.by = "none")$table$PValue
        uq <- EDASeq::betweenLaneNormalization(y, which = "upper", round = TRUE)
        E <- as.integer(sva::num.sv(n$log_cpm_tmm, model, method = "be", B = 20, seed = seed))
        list(d = d, y = y, n = n, model = model, first_p = first_p, uq = uq, E = E)
    })), error = function(e) e)
    prep_seconds <- proc.time()[["elapsed"]] - t_start
    if (inherits(prep, "error")) {
        status <- if (inherits(prep, "trx_error")) "refused" else "crashed"
        for (label in armsOf(cl)) {
            event(label, status, prep)
            record(label, status, NA_MEASURES, 0)
        }
        return(list(rows = do.call(rbind, rows), events = do.call(rbind, events),
            timing = data.frame(cell = cl$id, run = run, prep_seconds = prep_seconds, total_seconds = prep_seconds, E = NA_integer_)))
    }
    d <- prep$d; y <- prep$y; n <- prep$n; uq <- prep$uq; E <- prep$E; design <- d$design
    genes <- rownames(y)
    truth <- d$truth[genes]
    all_null <- d$null_all[genes]
    withW <- function(W) if (is.null(W) || !ncol(W)) noRuv() else if (assoc_mode) list(W = W) else
        if (cl$random) c(daTest(y, design, W), list(W = W)) else daTest(y, design, W)

    no_ruv <- arm("1", function() daTest(y, design))
    arm("2", function() withW(d$W))
    grid <- NULL
    arm("3", function() {
        grid <<- trxRuvGrid(y, n$log_cpm_tmm, design)
        s <- trxRuvs(y, n$log_cpm_tmm, design, grid$control_genes, grid$best$best_k)
        if (assoc_mode) return(list(W = s$W))
        table <- trxDaContrast(s$ruvs_counts, s$pdata, ENTRY)$top_table
        list(table = table[genes, c("logFC", "PValue", "FDR")], k = ncol(s$W))
    })
    gridK <- function() if (is.null(grid)) refuse("today's grid gave no k") else grid$best$best_k
    arm("5", function() {
        k <- gridK()
        controls <- genes[order(prep$first_p)[-seq_len(ceiling(0.25 * length(genes)))]]
        withW(RUVSeq::RUVg(uq, controls, k)$W)
    })
    arm("6", function() {
        k <- gridK()
        controls <- genes[stats::p.adjust(prep$first_p, "BH") > 0.5]
        if (length(controls) < 2L) refuse("fewer than two genes have a first-pass q above 0.5")
        withW(RUVSeq::RUVg(uq, controls, k)$W)
    })
    by_group <- RUVSeq::makeGroups(design$group)
    arm("7", function() if (E == 0L) noRuv() else withW(RUVSeq::RUVs(uq, genes, E, by_group)$W))
    oracle <- NULL
    arm("8", function() {
        if (E == 0L) return(noRuv())
        nulls <- genes[all_null]
        if (length(nulls) < 50L) refuse("fewer than 50 null genes pass the filter")
        oracle <<- RUVSeq::RUVg(uq, sample(nulls, 50L), E)
        withW(oracle$W)
    })
    arm("9", function() {
        if (E == 0L) return(noRuv())
        effects <- genes[truth != 0]
        n_bad <- min(10L, length(effects))
        controls <- c(sample(genes[all_null], 50L - n_bad), if (n_bad) sample(effects, n_bad))
        withW(RUVSeq::RUVg(uq, controls, E)$W)
    })
    arm("10", function() {
        if (E == 0L) return(noRuv())
        if (is.null(oracle)) refuse("arm 8 gave no fit")
        result <- daTest(oracle$normalizedCounts, design)
        result$k <- ncol(oracle$W)
        result
    })
    arm("11", function() withW(RUVSeq::RUVs(uq, genes, 1L, by_group)$W))
    arm("R", function() {
        if (E == 0L) return(noRuv())
        dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = y, group = design$group), method = "upperquartile")
        dge <- edgeR::estimateGLMTagwiseDisp(edgeR::estimateGLMCommonDisp(dge, prep$model), prep$model)
        residuals <- stats::residuals(edgeR::glmFit(dge, prep$model), type = "deviance")
        withW(RUVSeq::RUVr(uq, genes, E, residuals)$W)
    })
    sv <- function(n.sv) {
        utils::capture.output(fit <- sva::svaseq(uq, prep$model, prep$model[, 1L, drop = FALSE], n.sv = n.sv))
        if (fit$n.sv == 0L) noRuv() else withW(as.matrix(fit$sv))
    }
    arm("14", function() sv(NULL))
    arm("14E", function() if (E == 0L) noRuv() else sv(E))
    if (cl$tech) {
        sets <- NULL
        arm("15", function() {
            if (E == 0L) return(noRuv())
            sets <<- RUVSeq::RUVs(uq, genes, E, RUVSeq::makeGroups(d$unit))
            result <- daTest(sets$normalizedCounts, design)
            result$k <- ncol(sets$W)
            result
        })
        arm("16", function() {
            if (E == 0L) return(noRuv())
            if (is.null(sets)) refuse("arm 15 gave no fit")
            withW(sets$W)
        })
    }
    if (cl$batch) arm("17", function() {
        result <- daTest(y, data.frame(design, batch = d$batch), formula = "~ 0 + group + batch")
        result$k <- 1L
        result
    })
    separationK <- function(controls) {
        cancor_p <- trxCancor(Y = t(n$log_cpm_tmm), X = prep$model[, -1, drop = FALSE], ctl = rownames(n$log_cpm_tmm) %in% controls)
        trxFindBestK(cancor_p, max_k = if (ncol(y) < 15) 3L else 4L)
    }
    arm("22", function() {
        nulls <- genes[all_null]
        if (length(nulls) < 50L) refuse("fewer than 50 null genes pass the filter")
        set.seed(seed * 100L + ARM_STREAM[["8"]])
        controls <- sample(nulls, 50L)
        set.seed(seed * 100L + ARM_STREAM[["22"]])
        withW(RUVSeq::RUVg(uq, controls, separationK(controls))$W)
    })
    arm("23", function() {
        if (E == 0L) return(noRuv())
        controls <- genes[order(prep$first_p)[-seq_len(ceiling(0.25 * length(genes)))]]
        withW(RUVSeq::RUVg(uq, controls, E)$W)
    })
    arm("24", function() {
        controls <- genes[order(prep$first_p)[-seq_len(ceiling(0.25 * length(genes)))]]
        withW(RUVSeq::RUVg(uq, controls, separationK(controls))$W)
    })
    total <- proc.time()[["elapsed"]] - t_start
    c(list(rows = do.call(rbind, rows), events = if (length(events)) do.call(rbind, events),
        timing = data.frame(cell = cl$id, run = run, prep_seconds = prep_seconds, total_seconds = total, E = E)),
        if (length(side)) list(assoc = do.call(rbind, side)))
}

armsOf <- function(cl) {
    arms <- c("1", "2", "3", "5", "6", "7", "8", "9", "10", "11", "R", "14", "14E", if (cl$tech) c("15", "16"), if (cl$batch) "17")
    if (is.null(arms_wanted)) arms else intersect(c(arms, "22", "23", "24"), arms_wanted)
}

# A run whose worker died (or failed outside every arm) is a crash of each of its arms.
crashedRun <- function(cl, run, message) {
    seed <- 100000L + 10000L * cl$idx + run
    labels <- armsOf(cl)
    list(rows = if (assoc_mode) do.call(rbind, lapply(labels, function(label) assocRow("S3", cl$id, cl$class, run, seed, label, "crashed"))) else
            cbind(data.frame(simulator = "S3", cell = cl$id, cell_class = cl$class, run = run, seed = seed, arm = labels,
            subset = "full", status = "crashed", stringsAsFactors = FALSE), NA_MEASURES[rep(1L, length(labels)), ], seconds = NA_real_)[, COLUMNS],
        events = data.frame(cell = cl$id, run = run, arm = labels, status = "crashed", condition_class = "worker_failure",
            message = gsub("[\t\n]+", " ", message), stringsAsFactors = FALSE),
        timing = data.frame(cell = cl$id, run = run, prep_seconds = NA_real_, total_seconds = NA_real_, E = NA_integer_))
}

random_cells <- names(CELLS)[vapply(CELLS, `[[`, TRUE, "random")]
selected <- if (identical(cells_arg, "all")) setdiff(names(CELLS), random_cells) else strsplit(cells_arg, ",")[[1L]]
if (!all(selected %in% names(CELLS))) stop("unknown cell(s): ", paste(setdiff(selected, names(CELLS)), collapse = ", "))
if (!assoc_mode && any(selected %in% random_cells)) source("docs/qa/workflow-execution/kit/kit-129/r0/association.R")
dir.create(file.path(out_dir, "cells"), recursive = TRUE, showWarnings = FALSE)

for (id in selected) {
    cl <- CELLS[[id]]
    runs <- first_run - 1L + seq_len(if (cl$class == "core") runs_core else runs_stress)
    path <- file.path(out_dir, "cells", sprintf("s3-%s.rds", id))
    if (file.exists(path) && identical(readRDS(path)$runs, runs)) {
        message(sprintf("%s: present, skipped", id))
        next
    }
    t0 <- proc.time()[["elapsed"]]
    results <- parallel::mclapply(runs, function(run) tryCatch(oneRun(cl, run), error = function(e) crashedRun(cl, run, conditionMessage(e))),
        mc.cores = cores, mc.preschedule = FALSE)
    for (i in seq_along(results)) {
        if (!is.list(results[[i]]) || is.null(results[[i]]$rows)) results[[i]] <- crashedRun(cl, runs[[i]], paste(as.character(results[[i]]), collapse = " "))
    }
    part <- function(name) do.call(rbind, lapply(results, `[[`, name))
    cell_result <- list(runs = runs, rows = part("rows"), events = part("events"), timing = part("timing"),
        wall_seconds = proc.time()[["elapsed"]] - t0, cores = cores)
    if (!is.null(side <- part("assoc"))) cell_result$assoc <- side
    saveRDS(cell_result, paste0(path, ".tmp"))
    file.rename(paste0(path, ".tmp"), path)
    message(sprintf("%s: %d runs in %.0f s on %d cores", id, length(runs), cell_result$wall_seconds, cores))
}

# The tables hold every finished cell of the out dir, in cell order.
done <- Filter(function(id) file.exists(file.path(out_dir, "cells", sprintf("s3-%s.rds", id))), names(CELLS))
stored <- lapply(done, function(id) readRDS(file.path(out_dir, "cells", sprintf("s3-%s.rds", id))))
writeTable <- function(x, file) utils::write.table(x, file.path(out_dir, file), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
writeTable(do.call(rbind, lapply(stored, `[[`, "rows")), if (assoc_mode) "s3-assoc.tsv" else "s3-runs.tsv")
side <- if (!assoc_mode) do.call(rbind, lapply(stored, `[[`, "assoc"))   # the random-allocation cells' statistic
if (!is.null(side)) writeTable(side, "s3-assoc.tsv")
events <- do.call(rbind, lapply(stored, `[[`, "events"))
if (is.null(events)) events <- data.frame(cell = character(), run = integer(), arm = character(), status = character(),
    condition_class = character(), message = character())
writeTable(events, "s3-events.tsv")
writeTable(do.call(rbind, lapply(stored, `[[`, "timing")), "s3-timing.tsv")

versions <- vapply(c("edgeR", "limma", "RUVSeq", "EDASeq", "sva", "ruv", "pkgload"), function(p) as.character(utils::packageVersion(p)), "")
writeLines(c(
    "simulator: S3 (transcriptomics), KIT-124 revision 1.2",
    sprintf("commit: %s%s", system("git rev-parse HEAD", intern = TRUE),
        if (length(system("git status --porcelain -- R docs/qa/workflow-execution/kit/kit-124/sim", intern = TRUE))) " (R/ or sim/ differ from the commit)" else ""),
    sprintf("script sha256: %s", digest::digest(file = sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), algo = "sha256")),
    sprintf("command: Rscript %s", paste(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), args), collapse = " ")),
    sprintf("cells run by this command: %s; cells in the tables: %s", paste(selected, collapse = ","), paste(done, collapse = ",")),
    sprintf("runs: core %d, stress %d, first run %d; cores %d", runs_core, runs_stress, first_run, cores),
    sprintf("K124_ARMS: %s", if (assoc_mode) "unset" else if (is.null(arms_wanted)) "unset (arms 1 to 17 as the cell has them)" else paste(arms_wanted, collapse = ",")),
    if (assoc_mode) sprintf("K129_ASSOC: the association mode (KIT-129 R0 part B), arms %s, the RUV step only, no arm tested", paste(arms_wanted, collapse = ",")),
    "seed rule: seed = 100000 + 10000 * cell index + run; data set.seed(seed); arm a set.seed(seed * 100 + a) (R = 18, 14E = 19)",
    sprintf("library paths: %s", paste(.libPaths(), collapse = ":")),
    sprintf("R: %s; BLAS %s", R.version.string, extSoftVersion()[["BLAS"]]),
    sprintf("MultiScholaR (load_all): %s", as.character(utils::packageVersion("MultiScholaR"))),
    sprintf("packages: %s", paste(names(versions), versions, collapse = "; ")),
    sprintf("started %s; finished %s", format(started), format(Sys.time())),
    sprintf("wall seconds per cell: %s", paste(done, vapply(stored, function(s) sprintf("%.0f", s$wall_seconds), ""), collapse = "; "))),
    file.path(out_dir, "s3-runs.meta.txt"))
