# KIT-124 simulator S1 (ticket revision 1.2): a log2 protein matrix of 1,000 proteins through the package's classic
# sweep (findBestNegCtrlPercentage), RUV-III-C correction (ruvIII_C_Varying), RUV-III factors (computeRuvCovariates)
# and the limma kernel (runTestsContrasts), against the field's alternatives. Evidence only: nothing under R/ changes.
#
# Usage (from the repository root, R_LIBS set to the scratch library then the pinned library, one BLAS/OpenMP thread):
#   Rscript s1_classic.R <out dir> <runs per core cell> <runs per stress cell> <cores> [cells = all] [first run = 1]
# `cells` is "all" or a comma list (C0,S3); `first run` shifts the run numbers (the pilot uses 9001, so runs 1 to 1000
# stay unused for the confirming study). Writes <out dir>/s1-runs.tsv (one row per cell, run, arm and subset),
# s1-runs.meta.txt (commit, libraries, versions, command, seed rule, wall time per cell) and s1-runs.log (the message
# of every refusal and crash). Restartable per cell: a cell whose rows are complete for the requested runs is skipped;
# an incomplete cell's rows are dropped and the cell is rerun.
#
# ---------------------------------------------------------------------------------------------------------------
# PARAMETER TABLE
# ---------------------------------------------------------------------------------------------------------------
# Common to every cell (the centre, C0, is exactly this):
#   proteins            1,000 before the missing-value filter; ids PG_00001.. carry no information (the effects, the
#                       controls and the loadings are drawn at random positions, so no truth sits in the first rows)
#   abundance           mu_i ~ N(22, 2^2) log2
#   residual variance   s_i^2 = d0 * s0^2 / chisq(d0), d0 = 4, s0 = 0.25 (scaled inverse chi-square, limma's model);
#                       its median is d0 * s0^2 / qchisq(0.5, d0) = 0.0745 (the "median residual variance" below)
#   residuals           normal, independent per protein and sample
#   groups              2 ("control", "treated"), 6 samples each; contrast scored: grouptreated-groupcontrol
#   true effects        20% of proteins (exactly round(share * 1000), at random positions); |effect| = 1 log2 on the
#                       treated group; sign up for exactly round(0.5 * count), down for the rest
#   unwanted factor     1 factor; sample scores w with sample mean 0 and sample variance 1; loadings a_i ~ N(0, tau^2)
#                       independent of everything; strength = tau^2 / median residual variance = 3 (tau = 0.473)
#   factor and group    fixed by construction, not by chance: w = rho * g + sqrt(1 - rho^2) * e, where g is the
#                       contrast's group indicator (treated +1, control -1, other 0) scaled to sample variance 1 and
#                       e is a normal draw made exactly orthogonal to the whole group model (and to earlier factors,
#                       and to the batch where there is one) and scaled to sample variance 1; so the sample
#                       correlation of w with g is exactly rho. Centre: rho = 0.
#   missingness         intensity dependent: P(missing) = plogis((lod - value) / 1), lod solved per data set so that
#                       the overall missing share is the target (centre 25%); then the filter: a protein is kept when
#                       at least half of each group's samples are observed
#   technical           none (every sample its own unit; `tech_rep_group` unique per sample)
#   supplied controls   oracle set: 50 proteins without an effect in any group, drawn at random among those observed
#                       in every sample (among all such proteins when fewer than 50 are complete). Degraded set: 40 of
#                       those 50 plus 10 proteins with a true effect (drawn among the complete ones where possible).
#                       Both are drawn with the data, so arms 8, 9, 10, 12, 13 share them.
#
# Cells (index = the seed's cell index; only the settings that differ from the centre are listed):
#   idx cell class   setting
#    0  C0   core    centre
#    1  C1   core    12 per group
#    2  C2   core    10% effects
#    3  C3   core    no unwanted factor
#    4  C4   core    no true effects (global null), factor as the centre
#    5  C5   core    rho = 0.3
#    6  C6   core    technical duplicates: 4 units per group, 2 injections each (16 samples). Residual variance split
#                    75% biological (shared by a unit's 2 injections) and 25% technical (per injection), total s_i^2.
#                    The factor is per injection (rho = 0 with group). `tech_rep_group` names the unit.
#    7  C7   core    3 groups ("control", "treated", "other"), 6 each. 20% effects: half of them on the treated group,
#                    half on the other group (each +/- 1 log2 against control, signs as the centre). Contrast scored:
#                    grouptreated-groupcontrol; its true effects are the treated half (10% of proteins); the proteins
#                    changed in "other" are true nulls of the contrast but not negative controls. g is +1 treated,
#                    -1 control, 0 other.
#    8  S1   stress  3 per group
#    9  S2   stress  40% effects
#   10  S3   stress  rho = 0.6
#   11  S4   stress  3 factors, mutually orthogonal, each of strength 3; the first carries rho (0)
#   12  S5   stress  0% missing
#   13  S6   stress  50% missing
#   14  S7   stress  heavy tails: residuals t with 3 df scaled to variance 1 (t3 / sqrt(3)) times s_i
#   15  S8   stress  a known unbalanced batch added to the centre: batch b1 = 4 control + 2 treated, b2 = 2 control +
#                    4 treated (correlation of batch with group 1/3 by layout); a batch step with normal loadings of
#                    strength 3 on the standardised batch indicator; the centre's factor is kept, orthogonal to group
#                    and batch
#   16  S9   stress  80% of effects up
#   17  S10  stress  loadings tied to abundance: a_i = tau * (0.7 * z_i + sqrt(1 - 0.49) * N(0, 1)), z_i = (mu_i - 22) / 2
#   18  S11  stress  3 per group, rho = 0.3
#   19  S12  stress  3 per group, 3 factors
#   20  S13  stress  3 per group, no factor
#   21  S14  stress  no factor, 50% missing
#   22  C8   core    the exactness round (user, 2026-11-18, decision M16): C6 with no unwanted factor (technical
#                    duplicates, 4 units per group, 2 injections each; data generation as C6 with factors = 0). Run on
#                    request (a cell argument naming it); reported, judges no design.
#
# The test (every arm that tests with the package's kernel): the proteomics default path's call, as
# differentialAbundanceAnalysisHelper makes it with inst/config/config.ini [deAnalysisParameters] untouched:
#   runTestsContrasts(matrix, contrast_strings = "grouptreated-groupcontrol", design_matrix = design,
#     formula_string = "~ 0 + group" (+ covariates where the arm says), weights = NA, treat_lfc_cutoff = 0,
#     eBayes_trend = TRUE, eBayes_robust = TRUE, group_id = "group", technical_replicate_id = "tech_rep_group")
# A TREAT cutoff of 0 is not NA, so the kernel takes its TREAT branch: limma::treat(fit, lfc = 0) and topTreat, and
# treat() is called without trend or robust (review R-16: recorded, not acted on). p = raw_pvalue, q = fdr_qvalue
# (qvalue::qvalue, BH where it fails), a call is q < 0.05 (de_q_val_thresh). With technical duplicates (C6) the kernel
# runs duplicateCorrelation on the unit. The script stops if config.ini no longer holds these values.
#
# Today's sweep: findBestNegCtrlPercentage as resolveProtRuvParameters calls it on the GUI's defaults
# (percentage_range = 1:20, separation_metric = "max_difference", k_penalty_weight = 0.5, adaptive_k_penalty = TRUE,
# ruv_grouping_variable = "group"; the function's own defaults for the rest: qvalue, cutoff 0.05, 5 components).
#
# E: sva::num.sv(dat, mod = model.matrix(~ group), method = "be", B = 20, seed = stream 0) on the proteins observed
# in every sample (num.sv takes no missing values). One value per data set, shared by every arm that uses it. E = 0
# makes the arm no RUV (status ok, k = 0, the no-RUV result).
#
# Arms (W = RUV-III factors from computeRuvCovariates on the complete proteins, controls centred; "in the model" =
# covariates of the kernel's formula on the uncorrected matrix; "corrected" = ruvIII_C_Varying then the kernel):
#    1  no RUV                                                              every cell
#    2  the true factor scores in the model (and the batch in S8); in a cell without a factor this is the group
#       model, the true one                                                 every cell
#    3  today: sweep controls and k, corrected                              every cell
#    4  sweep controls and k, W in the model                                every cell
#    5  controls outside the 25% most significant of a first pass, sweep k, W in the model   every cell
#    6  controls with first-pass q > 0.5, sweep k, W in the model           every cell
#    7  all proteins as controls, k = E, W in the model                     every cell
#    8  oracle controls, k = E, W in the model                              every cell
#    9  degraded set, k = E, W in the model                                 every cell but C4 (no true effects)
#   10  oracle controls, k = E, corrected                                   every cell
#   11  all proteins as controls, k = 1, W in the model                     every cell
#   12  ruv::RUV4, oracle controls, k = E, variance_adjust p.ebayes         every cell, complete subset only
#   13  ruv::RUVinv, oracle controls, variance_adjust p.ebayes              every cell, complete subset only
#   14  sva::sva(n.sv = E) surrogate variables in the model                 every cell, complete subset only
#   15  technical replicates as the replicate structure, all proteins as controls, k = E, corrected   C6
#   16  as 15, W in the model                                               C6
#   17  the known batch in the model                                        S8
#   22  arm 8 with k from today's separation rule on the oracle controls    on request (K124_ARMS)
#   23  arm 5 with k = E                                                    on request (K124_ARMS)
#   24  arm 5 with k from today's separation rule on arm 5's own controls   on request (K124_ARMS)
# Arms 22 and 23 (user, 2026-11-18; 17 to 21 are taken as arm or stream numbers in one simulator or another) cross
# the two k rules over the two control sets: 22 differs from 8, and 23 from 5, in k only. Arm 22's k is the sweep's
# rule for one control set (its steps 3 and 4): the canonical-correlation curve of the oracle controls against all
# proteins (prepareRuvCancorInput with 5 components, ruvCancorFromInput), read up to calculateAdaptiveMaxK(samples),
# then findBestKElbow; an NA is a refusal. Arm 23 does not read the sweep, so a sweep refusal does not refuse it.
# Arm 24 (user, 2026-11-18): arm 5 takes the sweep's k, which the sweep computed on its own winning control set; arm
# 24 applies arm 22's rule to arm 5's controls (all proteins outside the 25% most significant of the first pass), so
# it differs from arm 5 in k only, and does not read the sweep either.
# The exactness round (user, 2026-11-18, decision M16; kit-124/exact-round-declaration.md), on request (K124_ARMS).
# Each arm is another arm with the named change and nothing else: it runs on that arm's random stream. Where only E
# changes (29 to 31) and the new E equals the base arm's k, the base arm's result is taken without computing it again.
#   25  arm 7 with the comparison's standard error exact per protein: the contrast is fitted as a coefficient
#       (limma::contrastAsCoef on the kernel's model matrix; the kernel's fit, fitProtDaBaseModel, on the
#       reparameterised matrix; its first coefficient taken in place of limma::contrasts.fit), then the kernel's
#       eBayes(trend, robust), treat(lfc = 0) and topTreat as before. contrasts.fit forms the standard error from one
#       correlation matrix for all proteins, which is approximate for a protein with missing values when the model
#       carries W.
#   26  arm 7 with a test that applies the settings the kernel sets: limma::treat(fit, lfc = 0, trend = TRUE,
#       robust = TRUE) in place of treat(fit, lfc = 0) (which moderates without trend and not robustly). At a cutoff
#       of 0 this equals the kernel without its TREAT step (topTable of eBayes(trend, robust)): one arm, not two (the
#       two p-values are compared under K124_EXACT_CHECK=1).
#   27  arm 1 with arm 26's test.            28  arm 3 with arm 26's test (the same sweep, controls, k and correction).
#   29  C6 and C8: arm 7 with E estimated after the biological units are removed, not only the groups:
#       sva::num.sv(dat, model.matrix(~ tech_rep_group), method = "be", B = 20, seed = stream 0).
#   30  C6 and C8: arm 16 with arm 29's E.
#   31  arm 7 with E from 200 permutations: sva::num.sv(..., B = 200, seed = stream 0), the group model.
#   32  arm 7 with all of them: arm 25's standard error, arm 26's test, B = 200, and in C6 and C8 the units removed.
#   33  C6 and C8: arm 16 with all of them (the design a declared technical-replicate design would take).
# With E = 0 arms 29 to 31 are arm 1's result, as arm 7; arms 25, 26, 32, 33 are the no-RUV fit under their own test,
# and a refusal of arms 25, 26, 28, 32, 33 is scored as that fit too (R-07 with the arm's own test).
# K124_EXACT_CHECK=1 (environment; a check, not a result): per data set, the simulator's rebuilt kernel with neither
# change must be identical() to the package's kernel (no RUV, and with arm 7's W), and treat(lfc = 0, trend, robust)
# must give topTable's p-values of eBayes(trend, robust) to 1e-12; otherwise the run crashes with the reason.
# K124_ARMS (environment; unset = the arms above, 22 to 33 not run): a comma list of the arms to run, of those the
# cell has and 22 to 33. Arm 1 must be in it (refusals carry its result). The data, E and every stream are as without
# it, so an arm's rows are the same whichever other arms run.
# First pass (arms 5, 6): limma lmFit(~ group) + eBayes, the F test over the group coefficients; q = BH.
# Replicate structure of RUV-III: arms 3 to 6 take what the package resolves from the design (the group; the unit in
# C6); arms 7 to 11 take the group in every cell (in C6 the object is declared without technical replicates for the
# RUV step only), so that arms 15 and 16 differ from arms 10 and 7 by the replicate structure. The test always has
# the true declaration. Arms 12 and 13: X = the group indicators with an intercept (Z = 1), p and the estimate of the
# treated column, variance_adjust(ebayes = TRUE, pooled = FALSE, evar = FALSE, rsvar = FALSE), q as the kernel
# (qvalue, BH where it fails).
#
# Subsets: `complete` = the proteins observed in every sample. Arms 12 to 14 are fitted and scored on them. Every
# other arm is fitted once on all proteins and scored on all (`full`) and on the complete ones (`complete`, the same
# fit's p and calls restricted to them).
#
# Status: `ok`; `refused` = a designed refusal, scored as the same data set's no-RUV result (R-07): the package's
# classed refusals (prot_ruv_insufficient_controls, prot_ruv_empty_result), the sweep's "No valid percentage found"
# stop, and the simulator's own of the same kind (class s1_refusal: no complete proteins for E; k not below the
# controls observed in every sample; k above samples minus replicate sets; no residual degree of freedom left for the
# model with k covariates; arms 4 to 6 when the sweep refused); `crashed` = any other error, measures NA.
#
# Measures per row: n_features = proteins of the subset with a finite p; n_true = true effects of the contrast in the
# subset; n_calls = q < 0.05; null_p05, null_p01 = share of the tested true nulls with p < 0.05, 0.01; fdp = false
# calls / calls (0 without calls); power = called true effects / n_true (NA when n_true = 0); bias, rmse = mean and
# root mean square of sign(effect) * (estimate - effect) over the tested true effects (negative bias = shrunk toward
# 0); seconds = the arm's elapsed time, the same on both of its subset rows (the sweep's time is in arm 3).
#
# Seeds: seed = 100000 + 10000 * cell index + run. The data (and the supplied control sets) are drawn once per (cell,
# run) under set.seed(seed) and shared by all arms. Each arm runs under set.seed(seed * 100 + arm number); E is
# stream 0. The sweep is computed once under stream 3 (inside arm 3, or by the first of arms 4 to 6 if arm 3 is
# dropped) and reused by arms 4 to 6.
# ---------------------------------------------------------------------------------------------------------------

args <- commandArgs(TRUE)
if (length(args) < 4L) stop("usage: Rscript s1_classic.R <out dir> <runs per core cell> <runs per stress cell> <cores> [cells = all] [first run = 1]")
out_dir <- args[[1L]]
runs_core <- as.integer(args[[2L]]); runs_stress <- as.integer(args[[3L]]); cores <- as.integer(args[[4L]])
cells_wanted <- if (length(args) > 4L) args[[5L]] else "all"
first_run <- if (length(args) > 5L) as.integer(args[[6L]]) else 1L
arms_wanted <- Sys.getenv("K124_ARMS", "")
arms_wanted <- if (nzchar(arms_wanted)) as.integer(strsplit(arms_wanted, ",", fixed = TRUE)[[1L]])
stopifnot(is.null(arms_wanted) || (!anyNA(arms_wanted) && 1L %in% arms_wanted))
stopifnot(first_run >= 1L, first_run + max(runs_core, runs_stress) - 1L < 10000L)

suppressMessages(pkgload::load_all(".", quiet = TRUE))
logger::log_threshold(logger::FATAL)
try(RUVIIIC::omp_set_num_threads(1L), silent = TRUE)
quietly <- function(expr) invisible(suppressWarnings(suppressMessages(utils::capture.output(expr, type = "output"))))

# ---- the default path's test settings, pinned against config.ini (R-16) -----------------------------------------
test_settings <- list(formula_string = "~ 0 + group", treat_lfc_cutoff = 0, eBayes_trend = TRUE, eBayes_robust = TRUE,
    q_threshold = 0.05, contrast = "grouptreated-groupcontrol")
local({
    lines <- readLines(file.path("inst", "config", "config.ini"))
    from <- which(trimws(lines) == "[deAnalysisParameters]")
    stopifnot(length(from) == 1L)
    block <- lines[seq(from + 1L, length(lines))]
    value <- function(key) trimws(sub("^[^=]*=", "", grep(paste0("^\\s*", key, "\\s*="), block, value = TRUE)[[1L]]))
    ok <- identical(value("formula_string"), test_settings$formula_string) &&
        as.double(value("treat_lfc_cutoff")) == test_settings$treat_lfc_cutoff &&
        identical(tolower(value("eBayes_trend")), "true") && identical(tolower(value("eBayes_robust")), "true") &&
        as.double(value("de_q_val_thresh")) == test_settings$q_threshold
    if (!ok) stop("inst/config/config.ini [deAnalysisParameters] no longer holds the settings this simulator pins")
})

# ---- cells -----------------------------------------------------------------------------------------------------
centre <- list(n = 6L, groups = 2L, share = 0.20, up = 0.5, factors = 1L, rho = 0, strength = 3, missing = 0.25,
    tails = "normal", batch = FALSE, duplicates = FALSE, abundance_tie = 0)
cell_table <- list(
    C0 = list(), C1 = list(n = 12L), C2 = list(share = 0.10), C3 = list(factors = 0L), C4 = list(share = 0),
    C5 = list(rho = 0.3), C6 = list(n = 4L, duplicates = TRUE), C7 = list(groups = 3L),
    S1 = list(n = 3L), S2 = list(share = 0.40), S3 = list(rho = 0.6), S4 = list(factors = 3L), S5 = list(missing = 0),
    S6 = list(missing = 0.50), S7 = list(tails = "t3"), S8 = list(batch = TRUE), S9 = list(up = 0.8),
    S10 = list(abundance_tie = 0.7), S11 = list(n = 3L, rho = 0.3), S12 = list(n = 3L, factors = 3L),
    S13 = list(n = 3L, factors = 0L), S14 = list(factors = 0L, missing = 0.50),
    C8 = list(n = 4L, duplicates = TRUE, factors = 0L))
cell_names <- names(cell_table)
cell_setting <- function(cell) utils::modifyList(centre, cell_table[[cell]])
cell_index <- function(cell) match(cell, cell_names) - 1L
cell_class <- function(cell) if (startsWith(cell, "C")) "core" else "stress"
arms_of <- function(cell) {
    p <- cell_setting(cell)
    arms <- 1:14
    if (p$share == 0) arms <- setdiff(arms, 9L)
    if (p$duplicates) arms <- c(arms, 15L, 16L)
    if (p$batch) arms <- c(arms, 17L)
    if (!is.null(arms_wanted)) arms <- intersect(c(arms, 22L, 23L, 24L, 25L, 26L, 27L, 28L, if (p$duplicates) c(29L, 30L),
        31L, 32L, if (p$duplicates) 33L), arms_wanted)
    arms
}
# the exactness round: the arms whose test differs from the kernel's (exact standard error, settings applied)
own_test <- list(`25` = c(TRUE, FALSE), `26` = c(FALSE, TRUE), `27` = c(FALSE, TRUE), `28` = c(FALSE, TRUE),
    `32` = c(TRUE, TRUE), `33` = c(TRUE, TRUE))
self_check <- identical(Sys.getenv("K124_EXACT_CHECK"), "1")
subsets_of <- function(arm) if (arm %in% 12:14) "complete" else c("full", "complete")

# ---- data ------------------------------------------------------------------------------------------------------
m_proteins <- 1000L; d0 <- 4; s0 <- 0.25
median_residual_variance <- d0 * s0^2 / stats::qchisq(0.5, d0)
pick <- function(pool, size) pool[sample.int(length(pool), min(size, length(pool)))]
simulate <- function(p, seed) {
    set.seed(seed)
    m <- m_proteins
    group_levels <- c("control", "treated", "other")[seq_len(p$groups)]
    per_unit <- if (p$duplicates) 2L else 1L
    units <- p$n * p$groups
    unit_of <- rep(seq_len(units), each = per_unit)
    groups <- rep(rep(group_levels, each = p$n), each = per_unit)
    n_samples <- length(groups)
    group <- factor(groups, levels = group_levels)
    mu <- stats::rnorm(m, 22, 2)
    s <- sqrt(d0 * s0^2 / stats::rchisq(m, d0))
    # true effects at random positions
    effect_treated <- numeric(m); effect_other <- numeric(m)
    n_effects <- round(p$share * m)
    if (n_effects > 0L) {
        where <- sample.int(m, n_effects)
        sign <- sample(rep(c(1, -1), c(round(p$up * n_effects), n_effects - round(p$up * n_effects))))
        on_other <- if (p$groups == 3L) sample(rep(c(FALSE, TRUE), length.out = n_effects)) else logical(n_effects)
        effect_treated[where[!on_other]] <- sign[!on_other]
        effect_other[where[on_other]] <- sign[on_other]
    }
    # the unwanted factors, their correlation with group fixed by construction
    scale1 <- function(v) { v <- v - mean(v); v / stats::sd(v) }
    g <- scale1((groups == "treated") - (groups == "control"))
    basis <- stats::model.matrix(~ group)
    batch <- NULL; unwanted <- NULL; loadings <- NULL
    tau <- sqrt(p$strength * median_residual_variance)
    if (p$batch) {
        batch <- character(n_samples)
        for (level in group_levels) {
            in_group <- which(groups == level)
            first <- if (level == "control") 4L else 2L
            batch[in_group] <- rep(c("b1", "b2"), c(first, length(in_group) - first))
        }
        step <- scale1(batch == "b2")
        basis <- cbind(basis, step)
        unwanted <- cbind(unwanted, batch_step = step)
        loadings <- cbind(loadings, stats::rnorm(m, 0, tau))
    }
    scores <- NULL
    for (f in seq_len(p$factors)) {
        e <- stats::rnorm(n_samples)
        e <- scale1(stats::lm.fit(basis, e)$residuals)
        basis <- cbind(basis, e)
        w <- if (f == 1L) p$rho * g + sqrt(1 - p$rho^2) * e else e
        a <- if (f == 1L && p$abundance_tie > 0) {
            tau * (p$abundance_tie * (mu - 22) / 2 + sqrt(1 - p$abundance_tie^2) * stats::rnorm(m))
        } else stats::rnorm(m, 0, tau)
        scores <- cbind(scores, w); loadings <- cbind(loadings, a)
    }
    if (!is.null(scores)) colnames(scores) <- paste0("F", seq_len(ncol(scores)))
    unwanted <- cbind(unwanted, scores)
    noise <- function(n) if (identical(p$tails, "t3")) stats::rt(n, 3) / sqrt(3) else stats::rnorm(n)
    x <- matrix(mu, m, n_samples) + outer(effect_treated, as.numeric(groups == "treated")) +
        outer(effect_other, as.numeric(groups == "other"))
    if (!is.null(unwanted)) x <- x + loadings %*% t(unwanted)
    if (p$duplicates) {
        biological <- matrix(noise(m * units), m, units) * s * sqrt(0.75)
        x <- x + biological[, unit_of, drop = FALSE] + matrix(noise(m * n_samples), m, n_samples) * s * sqrt(0.25)
    } else {
        x <- x + matrix(noise(m * n_samples), m, n_samples) * s
    }
    if (p$missing > 0) {
        lod <- stats::uniroot(function(l) mean(stats::plogis(l - x)) - p$missing, c(0, 60))$root
        x[stats::runif(length(x)) < stats::plogis(lod - x)] <- NA_real_
    }
    ids <- sprintf("PG_%05d", seq_len(m)); samples <- sprintf("run_%03d", seq_len(n_samples))
    dimnames(x) <- list(ids, samples)
    keep <- rep(TRUE, m)
    for (level in group_levels) {
        in_group <- groups == level
        keep <- keep & rowSums(is.finite(x[, in_group, drop = FALSE])) >= sum(in_group) / 2
    }
    x <- x[keep, , drop = FALSE]
    effect_treated <- stats::setNames(effect_treated, ids)[rownames(x)]
    any_effect <- (effect_treated != 0) | (stats::setNames(effect_other, ids)[rownames(x)] != 0)
    complete <- rowSums(!is.finite(x)) == 0
    # the supplied control sets
    null_pool <- which(!any_effect & complete); if (length(null_pool) < 50L) null_pool <- which(!any_effect)
    oracle_rows <- pick(null_pool, 50L)
    oracle <- stats::setNames(seq_len(nrow(x)) %in% oracle_rows, rownames(x))
    degraded <- NULL
    if (any(any_effect)) {
        effect_pool <- which(any_effect & complete); if (length(effect_pool) < 10L) effect_pool <- which(any_effect)
        degraded_rows <- c(utils::head(oracle_rows, max(length(oracle_rows) - 10L, 0L)), pick(effect_pool, 10L))
        degraded <- stats::setNames(seq_len(nrow(x)) %in% degraded_rows, rownames(x))
    }
    design <- data.frame(Run = samples, group = group, replicates = seq_len(n_samples),
        tech_rep_group = paste0("u", unit_of), no_tech = paste0("i", seq_len(n_samples)), stringsAsFactors = FALSE)
    if (!is.null(batch)) design$batch <- factor(batch)
    rownames(design) <- samples
    if (!is.null(scores)) rownames(scores) <- samples
    list(x = x, design = design, group = group, effect = effect_treated, complete = complete, oracle = oracle,
        degraded = degraded, scores = scores, has_batch = !is.null(batch))
}
make_object <- function(x, design, technical) {
    quant <- as.data.frame(c(list(Protein.Group = rownames(x)), as.data.frame(x)), check.names = FALSE, stringsAsFactors = FALSE)
    object_design <- design; rownames(object_design) <- NULL
    methods::new("ProteinQuantitativeData", protein_quant_table = quant, protein_id_column = "Protein.Group",
        protein_id_table = data.frame(Protein.Group = rownames(x), stringsAsFactors = FALSE), design_matrix = object_design,
        sample_id = "Run", group_id = "group", technical_replicate_id = technical,
        args = list(globalParameters = list(workflow_type = "LFQ")))
}
as_matrix <- function(object) {
    q <- object@protein_quant_table
    m <- as.matrix(q[, setdiff(names(q), object@protein_id_column), drop = FALSE]); rownames(m) <- q[[object@protein_id_column]]; m
}

# ---- the test and the measures -----------------------------------------------------------------------------------
refuse <- function(message) stop(structure(class = c("s1_refusal", "error", "condition"), list(message = message, call = NULL)))
is_refusal <- function(e) inherits(e, c("s1_refusal", "prot_ruv_insufficient_controls", "prot_ruv_empty_result")) ||
    grepl("No valid percentage found", conditionMessage(e), fixed = TRUE)
q_of <- function(p) {
    q <- rep(NA_real_, length(p)); valid <- is.finite(p)
    if (sum(valid) >= 3L) q[valid] <- tryCatch(qvalue::qvalue(p[valid])$q, error = function(e) stats::p.adjust(p[valid], "BH"))
    else q[valid] <- stats::p.adjust(p[valid], "BH")
    q
}
# The kernel as the default path calls it; `covariates` (samples by k, named) join the design and the formula.
da <- function(matrix, design, covariates = NULL, terms = character()) {
    if (!is.null(covariates)) {
        design <- cbind(design, as.data.frame(covariates[rownames(design), , drop = FALSE]))
        terms <- c(terms, colnames(covariates))
    }
    formula <- paste(c(test_settings$formula_string, terms), collapse = " + ")
    value <- NULL
    quietly(value <- runTestsContrasts(matrix, contrast_strings = test_settings$contrast, design_matrix = design,
        formula_string = formula, weights = NA, treat_lfc_cutoff = test_settings$treat_lfc_cutoff,
        eBayes_trend = test_settings$eBayes_trend, eBayes_robust = test_settings$eBayes_robust,
        group_id = "group", technical_replicate_id = "tech_rep_group"))
    tab <- value$results[[1L]]
    data.frame(id = rownames(tab), p = tab$raw_pvalue, q = tab$fdr_qvalue, estimate = tab$logFC, stringsAsFactors = FALSE)
}
# The exactness round: the kernel's steps from its own halves (fitProtDaBaseModel, applyProtDaContrastModel), so that
# the two changes can be made. `exact`: the contrast as a coefficient of a reparameterised model matrix
# (limma::contrastAsCoef), its first coefficient in place of contrasts.fit. `applied`: treat() given the kernel's
# trend and robust settings. With neither, the result is the kernel's (K124_EXACT_CHECK=1 asserts it).
# `no_treat` (the check only): topTable of the moderated fit, the kernel without its TREAT step.
da_variant <- function(matrix, design, covariates = NULL, exact = FALSE, applied = FALSE, no_treat = FALSE) {
    terms <- character()
    if (!is.null(covariates)) {
        design <- cbind(design, as.data.frame(covariates[rownames(design), , drop = FALSE]))
        terms <- colnames(covariates)
    }
    formula <- stats::as.formula(paste(c(test_settings$formula_string, terms), collapse = " + "))
    design_m <- stats::model.matrix(formula, stats::model.frame(formula, design))
    data_subset <- matrix[, rownames(design_m), drop = FALSE]
    contrast <- limma::makeContrasts(contrasts = test_settings$contrast, levels = colnames(design_m))
    fit_design <- design_m
    if (exact) {
        fit_design <- limma::contrastAsCoef(design_m, contrast)$design
        dimnames(fit_design) <- list(rownames(design_m), c(colnames(contrast), paste0("Q", seq_len(ncol(design_m))[-1L])))
    }
    silent <- function(...) invisible(NULL)
    model <- NULL
    quietly({
        base <- fitProtDaBaseModel(data_subset, fit_design, design, rownames(design_m), messageFn = silent,
            group_id = "group", technical_replicate_id = "tech_rep_group")
        model <- applyProtDaContrastModel(base, contrast, treat_lfc_cutoff = test_settings$treat_lfc_cutoff,
            eBayes_trend = test_settings$eBayes_trend, eBayes_robust = test_settings$eBayes_robust,
            contrastsFitFn = if (exact) function(fit, contrasts) fit[, 1L] else limma::contrasts.fit,
            treatFn = if (applied) function(fit, lfc) limma::treat(fit, lfc = lfc, trend = test_settings$eBayes_trend,
                robust = test_settings$eBayes_robust) else limma::treat, messageFn = silent)
    })
    tab <- if (no_treat) limma::topTable(model$moderated_fit, coef = 1L, number = Inf) else limma::topTreat(model$public_fit, coef = 1L, n = Inf)
    # q as the kernel: qvalue on the finite p-values, Benjamini-Hochberg when they are degenerate or qvalue fails
    p <- tab$P.Value; q <- rep(NA_real_, length(p)); valid <- which(!is.na(p) & is.finite(p))
    if (length(valid)) {
        v <- p[valid]
        q[valid] <- if (all(v == 0) || all(v == 1) || length(v) < 3L) stats::p.adjust(v, "BH") else
            tryCatch(qvalue::qvalue(v)$qvalues, error = function(e) stats::p.adjust(v, "BH"))
    }
    data.frame(id = rownames(tab), p = p, q = q, estimate = tab$logFC, stringsAsFactors = FALSE)
}
measures <- c("n_features", "n_true", "n_calls", "null_p05", "null_p01", "fdp", "power", "bias", "rmse")
score <- function(tab, sim, subset) {
    ids <- if (subset == "full") rownames(sim$x) else rownames(sim$x)[sim$complete]
    tab <- tab[tab$id %in% ids & is.finite(tab$p), , drop = FALSE]
    effect <- sim$effect[tab$id]
    true <- effect != 0
    called <- is.finite(tab$q) & tab$q < test_settings$q_threshold
    n_true <- sum(sim$effect[ids] != 0)
    error <- sign(effect[true]) * (tab$estimate[true] - effect[true])
    c(n_features = nrow(tab), n_true = n_true, n_calls = sum(called),
        null_p05 = if (any(!true)) mean(tab$p[!true] < 0.05) else NA_real_,
        null_p01 = if (any(!true)) mean(tab$p[!true] < 0.01) else NA_real_,
        fdp = if (any(called)) sum(called & !true) / sum(called) else 0,
        power = if (n_true > 0L) sum(called & true) / n_true else NA_real_,
        bias = if (any(true)) mean(error) else NA_real_, rmse = if (any(true)) sqrt(mean(error^2)) else NA_real_)
}

# ---- one data set, every arm --------------------------------------------------------------------------------------
one <- function(cell, run) {
    p <- cell_setting(cell)
    seed <- 100000L + 10000L * cell_index(cell) + run
    stream <- function(arm) seed * 100L + arm
    started <- proc.time()[[3L]]
    sim <- simulate(p, seed)
    x <- sim$x; design <- sim$design; complete <- sim$complete
    n_samples <- ncol(x); n_groups <- nlevels(sim$group)
    xc <- x[complete, , drop = FALSE]
    all_controls <- stats::setNames(rep(TRUE, nrow(x)), rownames(x))
    cache <- new.env()
    object <- function(technical) {                       # the full object, by replicate declaration
        key <- paste0("full_", technical)
        if (is.null(cache[[key]])) cache[[key]] <- make_object(x, design, technical)
        cache[[key]]
    }
    object_complete <- function(technical) {
        key <- paste0("complete_", technical)
        if (is.null(cache[[key]])) cache[[key]] <- make_object(xc, design, technical)
        cache[[key]]
    }
    # E, once per data set (stream 0)
    estimate <- tryCatch({
        if (nrow(xc) < 2L) refuse("no proteins observed in every sample: E is undefined")
        value <- NULL
        quietly(value <- sva::num.sv(xc, stats::model.matrix(~ group, design), method = "be", B = 20, seed = stream(0L)))
        as.integer(value)
    }, error = function(e) e)
    E <- function() if (inherits(estimate, "error")) stop(estimate) else estimate
    # the exactness round: E with 200 permutations, and with the biological units removed (stream 0, as E)
    E_variant <- function(B, units) {
        if (B == 20L && !units) return(E())
        key <- paste0("E_", B, "_", units)
        if (is.null(cache[[key]])) cache[[key]] <- tryCatch({
            if (nrow(xc) < 2L) refuse("no proteins observed in every sample: E is undefined")
            model <- if (units) stats::model.matrix(~ tech_rep_group, design) else stats::model.matrix(~ group, design)
            value <- NULL
            quietly(value <- sva::num.sv(xc, model, method = "be", B = B, seed = stream(0L)))
            as.integer(value)
        }, error = function(e) e)
        if (inherits(cache[[key]], "error")) stop(cache[[key]])
        cache[[key]]
    }
    no_ruv_under <- function(arm) {                       # the no-RUV fit under an arm's own test
        settings <- own_test[[as.character(arm)]]
        key <- paste0("none_", settings[[1L]], "_", settings[[2L]])
        if (is.null(cache[[key]])) cache[[key]] <- da_variant(x, design, exact = settings[[1L]], applied = settings[[2L]])
        cache[[key]]
    }
    # the sweep, once per data set (stream 3), its time charged to arm 3
    sweep <- function() {
        if (is.null(cache$sweep)) {
            set.seed(stream(3L))
            cache$sweep <- tryCatch({
                value <- NULL
                quietly(value <- findBestNegCtrlPercentage(normalised_protein_matrix_obj = object("tech_rep_group"),
                    percentage_range = seq(1, 20, by = 1), separation_metric = "max_difference", k_penalty_weight = 0.5,
                    adaptive_k_penalty = TRUE, ruv_grouping_variable = "group", verbose = FALSE))
                controls <- as.logical(value$best_control_genes_index); controls[is.na(controls)] <- FALSE
                list(k = validateProtNormRuvK(value$best_k, "ruv_k"), controls = stats::setNames(controls, rownames(x)))
            }, error = function(e) e)
        }
        if (inherits(cache$sweep, "error")) stop(cache$sweep)
        cache$sweep
    }
    sweep_for_reuse <- function() tryCatch(sweep(), error = function(e) {
        if (is_refusal(e)) refuse(paste("the sweep refused:", conditionMessage(e))) else stop(e)
    })
    first_pass <- function() {
        if (is.null(cache$first)) {
            fit <- limma::eBayes(limma::lmFit(x, stats::model.matrix(~ group, design)))
            cache$first <- limma::topTable(fit, coef = seq(2L, n_groups), number = Inf, sort.by = "none")$P.Value
        }
        cache$first
    }
    separation_k <- function(controls) {                  # today's rule for k, on a given control set
        plot <- NULL
        quietly(plot <- ruvCancorFromInput(prepareRuvCancorInput(object("no_tech"), 5), object("no_tech"), unname(controls), "group"))
        plot$data <- plot$data[plot$data$K <= calculateAdaptiveMaxK(n_samples), , drop = FALSE]
        k <- findBestKElbow(plot)
        if (is.na(k)) refuse("findBestKElbow returned NA for the control set")
        k
    }
    wide_controls <- function() {                         # arm 5's controls
        first <- first_pass()
        wide <- rank(first, ties.method = "first", na.last = "keep") > 0.25 * sum(is.finite(first))
        wide[is.na(wide)] <- FALSE
        stats::setNames(wide, rownames(x))
    }
    no_ruv <- list(same_as_no_ruv = TRUE, k = 0L)
    room <- function(k) if (n_samples - n_groups - k < 1L)
        refuse(sprintf("no residual degree of freedom for %d covariates with %d samples and %d groups", k, n_samples, n_groups))
    in_model <- function(controls, k, technical, test = da) {   # W of RUV-III as covariates of the test
        if (k == 0L) return(no_ruv)
        room(k)
        usable <- unname(controls[complete])
        if (sum(usable) <= k) refuse(sprintf("k = %d is not below the %d controls observed in every sample", k, sum(usable)))
        blocks <- length(unique(if (technical == "no_tech") as.character(design$group) else
            resolveProtRuvReplicateBlocks(design, "Run", "group", technical, "group")$design$.ruv_replicate_block))
        if (k > n_samples - blocks) refuse(sprintf("k = %d exceeds samples minus replicate sets (%d)", k, n_samples - blocks))
        value <- NULL
        quietly(value <- computeRuvCovariates(object_complete(technical), k = k, controls = usable, ruv_grouping_variable = "group"))
        if (self_check && identical(test, da)) {
            stopifnot("the rebuilt kernel is not the package's kernel (with W)" =
                identical(da(x, design, value$W), da_variant(x, design, value$W)))
        }
        list(tab = test(x, design, value$W), k = k)
    }
    corrected <- function(controls, k, technical, test = da) {   # RUV-III-C correction, then the test
        if (k == 0L) return(no_ruv)
        value <- NULL
        quietly(value <- ruvIII_C_Varying(object(technical), ruv_grouping_variable = "group", ruv_number_k = k, ctrl = controls))
        list(tab = test(as_matrix(value), design), k = k)
    }
    # the exactness round: `base` arm with the changes named, on the base arm's stream (E first: num.sv sets a seed)
    variant <- function(arm, base, B = 20L, units = FALSE) function() {
        settings <- own_test[[as.character(arm)]]
        test <- if (is.null(settings)) da else function(matrix, design, covariates = NULL)
            da_variant(matrix, design, covariates, exact = settings[[1L]], applied = settings[[2L]])
        k <- if (base %in% c(7L, 16L)) E_variant(B, units)
        kept <- if (is.null(settings)) cache[[paste0("result_", base)]]   # the base arm's result, when this arm's k is its k
        if (!is.null(kept) && identical(k, if (isTRUE(kept$same_as_no_ruv)) 0L else kept$k)) return(kept)
        set.seed(stream(base))
        result <- switch(as.character(base),
            `1` = list(tab = test(x, design), k = 0L),
            `3` = { s <- sweep(); corrected(s$controls, s$k, "tech_rep_group", test) },
            `7` = in_model(all_controls, k, "no_tech", test),
            `16` = in_model(all_controls, k, "tech_rep_group", test))
        if (isTRUE(result$same_as_no_ruv) && !is.null(settings)) result <- list(tab = no_ruv_under(arm), k = 0L)
        result
    }
    ruv_package <- function(method) {                     # ruv::RUV4 / ruv::RUVinv on the complete proteins
        control <- unname(sim$oracle[complete])
        if (!any(control)) refuse("no oracle control observed in every sample")
        X <- stats::model.matrix(~ group, design)[, -1L, drop = FALSE]
        if (method == "RUV4") {
            k <- E(); if (k == 0L) return(no_ruv)
            room(k)
            fit <- ruv::RUV4(t(xc), X, control, k = k, Z = 1)
        } else {
            k <- NA_integer_
            fit <- ruv::RUVinv(t(xc), X, control, Z = 1)
        }
        adjusted <- ruv::variance_adjust(fit, ebayes = TRUE, pooled = FALSE, evar = FALSE, rsvar = FALSE)
        pv <- as.numeric(adjusted$p.ebayes[1L, ])
        list(tab = data.frame(id = rownames(xc), p = pv, q = q_of(pv), estimate = as.numeric(fit$betahat[1L, ]),
            stringsAsFactors = FALSE), k = k)
    }
    arm_functions <- list(
        `1` = function() list(tab = da(x, design), k = 0L),
        `2` = function() list(tab = da(x, design, sim$scores, if (sim$has_batch) "batch" else character()),
            k = if (is.null(sim$scores)) 0L else ncol(sim$scores)),
        `3` = function() { s <- sweep(); corrected(s$controls, s$k, "tech_rep_group") },
        `4` = function() { s <- sweep_for_reuse(); in_model(s$controls, s$k, "tech_rep_group") },
        `5` = function() {
            s <- sweep_for_reuse(); first <- first_pass()
            wide <- rank(first, ties.method = "first", na.last = "keep") > 0.25 * sum(is.finite(first))
            wide[is.na(wide)] <- FALSE
            in_model(stats::setNames(wide, rownames(x)), s$k, "tech_rep_group")
        },
        `6` = function() {
            s <- sweep_for_reuse(); q <- stats::p.adjust(first_pass(), "BH")
            in_model(stats::setNames(is.finite(q) & q > 0.5, rownames(x)), s$k, "tech_rep_group")
        },
        `7` = function() in_model(all_controls, E(), "no_tech"),
        `8` = function() in_model(sim$oracle, E(), "no_tech"),
        `9` = function() in_model(sim$degraded, E(), "no_tech"),
        `10` = function() corrected(sim$oracle, E(), "no_tech"),
        `11` = function() in_model(all_controls, 1L, "no_tech"),
        `12` = function() ruv_package("RUV4"),
        `13` = function() ruv_package("RUVinv"),
        `14` = function() {
            k <- E(); if (k == 0L) return(no_ruv)
            room(k)
            value <- NULL
            quietly(value <- sva::sva(xc, stats::model.matrix(~ group, design), stats::model.matrix(~ 1, design), n.sv = k))
            sv <- value$sv; dimnames(sv) <- list(colnames(xc), paste0("SV", seq_len(ncol(sv))))
            list(tab = da(xc, design, sv), k = k)
        },
        `15` = function() corrected(all_controls, E(), "tech_rep_group"),
        `16` = function() in_model(all_controls, E(), "tech_rep_group"),
        `17` = function() list(tab = da(x, design, NULL, "batch"), k = 0L),
        `22` = function() in_model(sim$oracle, separation_k(sim$oracle), "no_tech"),
        `23` = function() in_model(wide_controls(), E(), "tech_rep_group"),
        `24` = function() in_model(wide_controls(), separation_k(wide_controls()), "tech_rep_group"),
        `25` = variant(25L, 7L), `26` = variant(26L, 7L), `27` = variant(27L, 1L), `28` = variant(28L, 3L),
        `29` = variant(29L, 7L, units = TRUE), `30` = variant(30L, 16L, units = TRUE), `31` = variant(31L, 7L, B = 200L),
        `32` = variant(32L, 7L, B = 200L, units = p$duplicates), `33` = variant(33L, 16L, B = 200L, units = TRUE))
    if (self_check) {
        stopifnot("the rebuilt kernel is not the package's kernel (no RUV)" = identical(da(x, design), da_variant(x, design)))
        with_settings <- da_variant(x, design, applied = TRUE); without_treat <- da_variant(x, design, applied = TRUE, no_treat = TRUE)
        gap <- max(abs(with_settings$p - without_treat$p[match(with_settings$id, without_treat$id)]))
        if (!is.finite(gap) || gap > 1e-12) stop(sprintf("treat(lfc = 0, trend, robust) is not eBayes(trend, robust): p differs by up to %.3g", gap))
    }
    rows <- list(); log <- character()
    reference <- NULL                                     # arm 1's table, the no-RUV result
    for (arm in arms_of(cell)) {
        set.seed(stream(arm)); t0 <- proc.time()[[3L]]
        result <- tryCatch(arm_functions[[as.character(arm)]](), error = function(e) e)
        seconds <- proc.time()[[3L]] - t0
        if (arm %in% c(7L, 16L) && !inherits(result, "error")) cache[[paste0("result_", arm)]] <- result
        status <- "ok"; k <- NA_integer_; tab <- NULL
        if (inherits(result, "error")) {
            status <- if (is_refusal(result)) "refused" else "crashed"
            log <- c(log, paste(cell, run, arm, status, gsub("[\r\n\t]+", " ", conditionMessage(result)), sep = "\t"))
            if (status == "refused") tab <- if (is.null(own_test[[as.character(arm)]])) reference else
                tryCatch(no_ruv_under(arm), error = function(e) NULL)
        } else if (isTRUE(result$same_as_no_ruv)) {
            k <- 0L; tab <- reference
        } else {
            k <- result$k; tab <- result$tab
        }
        if (arm == 1L && status == "ok") reference <- tab
        for (subset in subsets_of(arm)) {
            scored <- if (is.null(tab) || status == "crashed") stats::setNames(rep(NA_real_, length(measures)), measures) else
                tryCatch(score(tab, sim, subset), error = function(e) {
                    log <<- c(log, paste(cell, run, arm, "crashed", paste("scoring:", conditionMessage(e)), sep = "\t"))
                    stats::setNames(rep(NA_real_, length(measures)), measures)
                })
            row_status <- if (status != "crashed" && all(is.na(scored))) "crashed" else status
            rows[[length(rows) + 1L]] <- data.frame(simulator = "S1", cell = cell, cell_class = cell_class(cell), run = run,
                seed = seed, arm = arm, subset = subset, status = row_status, k = k, t(scored), seconds = round(seconds, 3),
                stringsAsFactors = FALSE)
        }
    }
    list(rows = do.call(rbind, rows), log = log, total = proc.time()[[3L]] - started)
}
crashed_run <- function(cell, run, message) {
    grid <- do.call(rbind, lapply(arms_of(cell), function(arm) data.frame(arm = arm, subset = subsets_of(arm), stringsAsFactors = FALSE)))
    blank <- as.data.frame(as.list(stats::setNames(rep(NA_real_, length(measures)), measures)))
    list(rows = data.frame(simulator = "S1", cell = cell, cell_class = cell_class(cell), run = run,
        seed = 100000L + 10000L * cell_index(cell) + run, arm = grid$arm, subset = grid$subset, status = "crashed",
        k = NA_integer_, blank, seconds = NA_real_, stringsAsFactors = FALSE),
        log = paste(cell, run, "all", "crashed", gsub("[\r\n\t]+", " ", message), sep = "\t"), total = NA_real_)
}

# ---- the study ---------------------------------------------------------------------------------------------------
columns <- c("simulator", "cell", "cell_class", "run", "seed", "arm", "subset", "status", "k", measures, "seconds")
cells <- if (identical(cells_wanted, "all")) setdiff(cell_names, "C8") else strsplit(cells_wanted, ",", fixed = TRUE)[[1L]]
stopifnot(all(cells %in% cell_names), cores >= 1L)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
runs_file <- file.path(out_dir, "s1-runs.tsv"); meta_file <- file.path(out_dir, "s1-runs.meta.txt")
log_file <- file.path(out_dir, "s1-runs.log")
runs_of <- function(cell) first_run - 1L + seq_len(if (cell_class(cell) == "core") runs_core else runs_stress)
done <- if (file.exists(runs_file)) utils::read.delim(runs_file, stringsAsFactors = FALSE) else NULL
rows_per_run <- function(cell) sum(lengths(lapply(arms_of(cell), subsets_of)))
is_complete <- function(cell) {
    if (is.null(done)) return(FALSE)
    have <- table(factor(done$run[done$cell == cell], levels = runs_of(cell)))
    all(have == rows_per_run(cell))
}
meta <- function(...) cat(..., "\n", file = meta_file, append = TRUE, sep = "")
versions <- function(packages) paste(vapply(sort(packages), function(p) paste(p, as.character(utils::packageVersion(p))), ""), collapse = "; ")
meta("== s1_classic.R started ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), " ==")
meta("commit: ", tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE), error = function(e) "unknown"),
    "; uncommitted changes in the work tree: ", tryCatch(length(system2("git", c("status", "--porcelain", "--untracked-files=no"), stdout = TRUE)) > 0L, error = function(e) NA))
meta("script md5: ", tryCatch(unname(tools::md5sum(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))), error = function(e) "unknown"))
meta("command: ", paste(commandArgs(FALSE), collapse = " "))
meta("working directory: ", getwd())
meta("R_LIBS: ", Sys.getenv("R_LIBS")); meta(".libPaths: ", paste(.libPaths(), collapse = ":"))
meta("threads: OMP_NUM_THREADS=", Sys.getenv("OMP_NUM_THREADS"), " OPENBLAS_NUM_THREADS=", Sys.getenv("OPENBLAS_NUM_THREADS"),
    "; worker processes: ", cores, "; BLAS: ", extSoftVersion()[["BLAS"]])
meta(R.version.string, "; ", R.version$platform)
meta("MultiScholaR (pkgload::load_all of the work tree): ", as.character(utils::packageVersion("MultiScholaR")))
meta("seed rule: seed = 100000 + 10000 * cell index + run (C0..C7 = 0..7, S1..S14 = 8..21, C8 = 22); data under set.seed(seed); ",
    "arm a under set.seed(seed * 100 + a) (arms 25 to 33: the stream of the arm each varies); E under stream 0; runs ", first_run, " onward")
meta("K124_ARMS: ", if (is.null(arms_wanted)) "unset (arms 1 to 17 as the cell has them)" else paste(arms_wanted, collapse = ","))
meta("cells: ", paste(cells, collapse = ","), "; runs per core cell ", runs_core, ", per stress cell ", runs_stress)
for (cell in cells) {
    if (is_complete(cell)) { meta("cell ", cell, ": complete already, skipped"); next }
    if (!is.null(done) && any(done$cell == cell)) {       # an incomplete cell is rerun whole
        done <- done[done$cell != cell, , drop = FALSE]
        utils::write.table(done[, columns], runs_file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
    started <- proc.time()[[3L]]
    runs <- runs_of(cell)
    results <- parallel::mclapply(runs, function(run) tryCatch(one(cell, run), error = function(e) crashed_run(cell, run, conditionMessage(e))),
        mc.cores = cores, mc.preschedule = FALSE)
    for (i in seq_along(results)) if (!is.list(results[[i]]) || !is.data.frame(results[[i]]$rows))
        results[[i]] <- crashed_run(cell, runs[[i]], paste("worker failed:", paste(as.character(results[[i]]), collapse = " ")))
    rows <- do.call(rbind, lapply(results, `[[`, "rows"))[, columns]
    new_file <- !file.exists(runs_file)
    utils::write.table(rows, runs_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = new_file, append = !new_file)
    messages <- unlist(lapply(results, `[[`, "log"))
    if (length(messages)) cat(messages, sep = "\n", file = log_file, append = TRUE)
    done <- rbind(if (!is.null(done)) done[, columns], rows)
    meta(sprintf("cell %s: %d runs, wall %.0f s, mean seconds per run %.1f (all arms, data and E included), refused rows %d, crashed rows %d",
        cell, length(runs), proc.time()[[3L]] - started, mean(vapply(results, function(r) r$total, 0), na.rm = TRUE),
        sum(rows$status == "refused"), sum(rows$status == "crashed")))
    cat(sprintf("%s %s done\n", format(Sys.time(), "%H:%M:%S"), cell))
}
meta("loaded namespaces: ", versions(loadedNamespaces()))
meta("== finished ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), " ==")
