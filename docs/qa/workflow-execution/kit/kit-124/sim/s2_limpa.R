#!/usr/bin/env Rscript
# KIT-124 simulator S2: RUV on the DIA_LIMPA route (ticket revision 1.2, "Simulation Setup"; review R-05, R-07,
# R-10, R-11, R-13, R-16, R-19). Evidence only: nothing under R/, tests/ or tools/ is changed or sourced apart
# from pkgload::load_all(".").
#
# Usage (from the repository root, R_LIBS as the sidecar records):
#   Rscript docs/qa/workflow-execution/kit/kit-124/sim/s2_limpa.R <out dir> <runs per core cell> \
#       <runs per stress cell> <cores> [cells = all | comma-separated, e.g. C0,S3]
#   K124_FIRST_RUN (environment, default 1) is the first run number; runs are consecutive from it.
# Writes <out dir>/s2-runs.tsv (one row per cell, run, arm and subset), s2-runs.meta.txt (commit, libraries,
# versions, command, seed rule), s2-data.tsv (one row per cell and run: what was generated and quantified, the
# estimator E, the route's sweep, and the seconds of the steps all arms share) and s2-errors.tsv (the message of
# every refusal and crash). Restartable per cell: <out dir>/cells/<cell>.rds is written when a cell completes
# and is reused by a later invocation with the same run numbers.
#
# ---- What is called (the package's own DIA_LIMPA route) ---------------------------------------------------------
# A default GUI run reaches these kernels through the saved-state store and the Shiny session; the simulator
# calls the kernels themselves with the arguments that path gives them. The shipped defaults are read from the
# package (readWorkflowShippedDefaults, resolveDIARuvCovariatesParameters) and asserted at start (R-16).
#   Quantification  proteinMissingValueImputationLimpa(<PeptideQuantitativeData>, dpc_slope = 0.8,
#                   verbose = FALSE, chunk = 1000, engine = <shipped: fastlimpadpc>, peptide_summary = <counts>)
#                   as computeDIAReferenceProteinDpcQuant calls it in its default `bounded` variant (the object
#                   holds the log2 matrix, the design and the distinct protein/precursor identities; the
#                   per-protein counts are handed over). Inside, the method estimates the DPC on the fly with
#                   <engine>::dpc(y) and quantifies with <engine>::dpcQuant(y, "protein.id", dpc = <that fit>,
#                   verbose = FALSE, chunk = 1000). The shipped engine is fastlimpadpc (the package pins it to
#                   one Rayon thread), not limpa: that is what a default GUI run uses where it is installed.
#   RUV optimiser   findBestNegCtrlPercentageCapped(<protein object>, percentage_range = seq(1, 30, by = 1),
#                   num_components_to_impute = 5, ruv_grouping_variable = "group", ruv_qval_cutoff = 0.05,
#                   ruv_fdr_method = "BH", separation_metric = "max_difference", k_penalty_weight = 0.5,
#                   adaptive_k_penalty = TRUE, verbose = FALSE), as diaRuvCovariatesOperation's kernel.
#   RUV factors     computeRuvCovariates(<protein object>, k, controls, "group"): ruv::RUVIII with the design's
#                   replicate structure (the declared technical replicates where the design has them, else the
#                   group), W from the centred controls, adjusted = Y - W alpha.
#   Covariates mode W1..Wk appended to the design, the quantified EList unchanged (the kernel's `with_w`).
#   Correction mode the adjusted matrix replaces the EList's E; standard errors and observation counts kept; no
#                   W terms (the kernel's `corrected`; here the EList only, since the test reads nothing else).
#   Test            the dpcDE branch of differentialAbundanceAnalysisHelper,ProteinQuantitativeData:
#                   design = model.matrix(as.formula("~ 0 + group" [+ W1 + ... + Wk]), design table);
#                   block = resolveProtTechnicalReplicateBlock(design, "group", "tech_rep_group")$block when
#                   the design has technical replicates, else NULL; limpa::dpcDE(elist[, samples], design,
#                   plot = FALSE, block = block) (dpcDE forwards block to voomaLmFitWithImputation, which runs
#                   duplicateCorrelation itself); convertDpcDAToStandardFormat(dpc_fit, "B_vs_A=groupB-groupA",
#                   design, eBayes_trend = TRUE, eBayes_robust = TRUE): contrasts.fit, eBayes, topTable, and
#                   adjustProtDaResultTable (qvalue::qvalue, Benjamini-Hochberg when qvalue cannot run).
#   Calls           the route's significance rule: fdr_qvalue < 0.05 (da_q_val_thresh; no fold-change cut: the
#                   dpcDE branch applies no TREAT, and the shipped treat_lfc_cutoff is 0).
# Beneath the GUI layer: the method differentialAbundanceAnalysisHelper also draws plots and builds result
# tables around the test; the simulator calls the test's lines only. The store checkpoints, the worker and the
# publication halves of the route are not run. Steps of a GUI run that are not in the simulator: the precursor
# q-value filter, the roll-up, the peptide accession cleanup, the between-sample normalisation (shipped: none),
# the peptide-count filter (a no-op here: every protein has 4 precursors) and the sample correlation filter.
#
# ---- Data --------------------------------------------------------------------------------------------------------
# limpa::simProteinDataSet (limpa 1.4.2) with n.peptides = 4000, peptides.per.protein = 4 (1,000 proteins),
# mu.range = c(2, 10), dpc.slope = 0.7, prop.missing = the cell's share, prop.de = 0. The function ignores its
# `sigma` argument (every precursor has sd 0.6), spaces the precursor means evenly over mu.range in row order,
# draws no protein-level term, and draws missingness inside. Therefore (R-11) everything below is added to the
# complete precursor matrix it returns (other$E.complete), each protein's term shared by its 4 precursors, and
# missingness is then redrawn exactly as the function draws it, with the DPC coefficients it returns:
# P(missing) = plogis(intercept + 0.7 * y, lower.tail = FALSE). Missingness, and so the standard errors of the
# quantification, depend on the group effects and on the unwanted factor.
#   Group effects   placed by the simulator in the function's own form (-0.5 log2 in group A and +0.5 in group
#                   B for an up protein, the reverse for a down protein; a third group is untouched), on a random
#                   set of proteins. The function's own placement is not used because it cannot give 80% up
#                   (S9); one placement serves every cell. True difference B - A: +1 or -1 log2.
#   Unwanted factor scores w_j = rho * g_j + sqrt(1 - rho^2) * z_j, z ~ N(0, 1), g the standardised B-versus-A
#                   indicator: the correlation with group is rho in expectation and varies run to run by chance
#                   (the realised value is in s2-data.tsv). Loadings a_p ~ N(0, strength * 0.09): `strength` is
#                   the ratio of the factor's variance to the nominal residual variance of a protein with its
#                   4 precursors complete (0.6^2 / 4 = 0.09). The realised median residual variance of the
#                   quantified proteins is larger with missing values and is recorded in s2-data.tsv.
#   Three factors   independent scores and loadings, strength 1 each (the centre's total of 3 spread over 3).
#   Heavy tails     the precursor noise is replaced by a t with 3 df scaled to sd 0.6.
#   Batch (S8)      as well as the centre's factor: two batches, 4 of A and 2 of B in the first, 2 of A and 4 of
#                   B in the second; a per-protein step N(0, 3 * 0.09) times the standardised batch indicator.
#   Abundance (S10) loadings a_p = sd * (0.7 * z(mu_p) + sqrt(1 - 0.49) * N(0, 1)), mu_p the protein's mean.
#   Duplicates (C6) 4 units per group, 2 injections each (16 samples). A per-protein unit term N(0, 0.3^2) is
#                   shared by a unit's two injections (the function has no biological variance component); the
#                   precursor noise and the unwanted factor are per injection. The design declares the unit as
#                   the technical replicate (tech_rep_group = <group>_<unit>), so the route passes the block to
#                   dpcDE and uses the units as RUV-III's replicate structure in every package arm.
#   Rows            the precursor rows are permuted before the package sees them (R-19).
# The simulated values are log2 on the function's scale (2 to 10), handed over as a logged peptide state.
#
# ---- Cells (index for the seed: C0..C7 = 0..7, S1..S14 = 11..24) ------------------------------------------------
#   cell  class   groups  samples/group        effects  up    factors (strength)  rho  missing  other
#   C0    core    2       6                    20%      50%   1 (3)               0    25%
#   C1    core    2       12                   20%      50%   1 (3)               0    25%
#   C2    core    2       6                    10%      50%   1 (3)               0    25%
#   C3    core    2       6                    20%      50%   0                   -    25%
#   C4    core    2       6                    0%       -     1 (3)               0    25%      global null
#   C5    core    2       6                    20%      50%   1 (3)               0.3  25%
#   C6    core    2       8 (4 units x 2)      20%      50%   1 (3)               0    25%      unit sd 0.3
#   C7    core    3       6                    20%      50%   1 (3)               0    25%      B - A scored
#   S1    stress  2       3                    20%      50%   1 (3)               0    25%
#   S2    stress  2       6                    40%      50%   1 (3)               0    25%
#   S3    stress  2       6                    20%      50%   1 (3)               0.6  25%
#   S4    stress  2       6                    20%      50%   3 (1 each)          0    25%
#   S5    stress  2       6                    20%      50%   1 (3)               0    0%       no missingness drawn
#   S6    stress  2       6                    20%      50%   1 (3)               0    50%
#   S7    stress  2       6                    20%      50%   1 (3)               0    25%      t3 precursor noise
#   S8    stress  2       6                    20%      50%   1 (3)               0    25%      batch 4+2 / 2+4 (3)
#   S9    stress  2       6                    20%      80%   1 (3)               0    25%
#   S10   stress  2       6                    20%      50%   1 (3)               0    25%      loadings ~ abundance 0.7
#   S11   stress  2       3                    20%      50%   1 (3)               0.3  25%
#   S12   stress  2       3                    20%      50%   3 (1 each)          0    25%
#   S13   stress  2       3                    20%      50%   0                   -    25%
#   S14   stress  2       6                    20%      50%   0                   -    50%
# Missing shares map to the DPC intercept through the function's own solution for prop.missing at slope 0.7
# (about -2.53 for 25% and -4.20 for 50%); the redraw happens after the effects and the factor are added, so the
# realised share differs slightly and is recorded.
#
# ---- Arms (ticket numbering; all on the one quantified EList of the run) ----------------------------------------
#    1  no RUV: the route's test.
#    2  the true factor scores (and the batch in S8) in the model. In a cell without a factor: arm 1.
#    3  the route's correction mode: its sweep, computeRuvCovariates, test on the corrected EList.
#    4  the route's covariates mode: the same sweep and factors, W1..Wk in the model.
#    5  controls = all proteins outside the 25% most significant of arm 1; k = the sweep's; covariates.
#    6  controls = proteins with arm 1's q above 0.5; k = the sweep's; covariates.
#    7  all proteins as controls; k = E; covariates.
#    8  oracle controls (50 true nulls); k = E; covariates.
#    9  degraded supplied set (40 true nulls and 10 true effects); k = E; covariates. Not in C4 (no effects).
#   10  arm 8's controls, subtract then test (the corrected EList).
#   11  arm 7 with k = 1.
#   12  ruv::RUV4(Y, X, ctl = arm 8's controls, k = E), p.ebayes of ruv::variance_adjust; quantified matrix.
#   13  ruv::RUVinv(Y, X, ctl = arm 8's controls), p.ebayes; k is not a parameter (NA).
#   14  sva::sva(n.sv = E) surrogate variables in limma (lmFit, contrasts.fit, eBayes(trend, robust)).
#   15  C6 only: technical replicates as the replicate structure, all proteins as controls, k = E, subtract
#       then test.
#   16  C6 only: arm 15 with the factors in the model. On this route computeRuvCovariates already takes the
#       declared technical replicates, so arm 16 repeats arm 7's computation in C6.
#   17  S8 only: the known batch in the model, no RUV.
#   22  on request (K124_ARMS): arm 8 with k from the route's separation rule on the oracle controls.
#   23  on request (K124_ARMS): arm 5 with k = E.
#   24  on request (K124_ARMS): arm 5 with k from the route's separation rule on arm 5's own controls.
# Arms 22 and 23 (user, 2026-11-18; 17 to 21 are taken as arm or stream numbers in one simulator or another) cross
# the two k rules over the two control sets: 22 differs from 8, and 23 from 5, in k only. Arm 22's k is
# findBestNegCtrlPercentageCapped's rule for one control set: the canonical-correlation curve of the oracle controls
# against all proteins (prepareRuvCancorInput, ruvCancorFromInput), capped at min(calculateAdaptiveMaxK(samples),
# controls - 1), k at the largest all-minus-control difference. Arm 23 with E = 0 is no RUV, as arms 7 to 10.
# Arm 24 (user, 2026-11-18): arm 5 takes the sweep's k, which the sweep computed on its own winning control set; arm
# 24 applies arm 22's rule to arm 5's controls (all proteins outside the 25% most significant of arm 1), so it
# differs from arm 5 in k only, and does not read the sweep.
# K124_ARMS (environment; unset = the arms above, 22 to 24 not run): a comma list of the arms to run, of those the
# cell has and 22, 23, 24. Arm 1 must be in it (arms 5, 6, 23, 24 and every refusal read its result). Without any of arms
# 3 to 6 the route's sweep is not run (sweep_status "not run" in s2-data.tsv). The data, the control sets, E and every
# stream are as without it, so an arm's rows are the same whichever other arms run.
# Arms 5 to 11, 15 and 16 use computeRuvCovariates and the route's test. Arms 12 to 14 run on the quantified
# matrix (complete on this route), ignore the technical-replicate block and the DPC standard errors, and their
# q-values come from the package's rule applied to their p-values (adjustProtDaResultTable). X is the group B
# indicator; in C7 the group C indicator is in Z (arms 12, 13) and in the design (arm 14).
# E = sva::num.sv(quantified matrix, model.matrix(~ group), method = "be", B = 20, seed = 100 * seed + 99).
# When E is 0, arms 7 to 10, 15 and 16 are no RUV (arm 1's result, k = 0); arms 12 and 14 run their own k = 0
# form. Oracle and degraded sets are drawn once per run and shared (arms 8, 10, 12, 13; arm 9).
#
# ---- Seeds, statuses, measures -----------------------------------------------------------------------------------
# seed = 100000 + 10000 * cell index + run. Data: set.seed(seed). Arm a: set.seed(100 * seed + a) (the sweep
# runs on arm 3's stream; E's seed is 100 * seed + 99). Data are generated and quantified once per (cell, run).
# status: ok; refused = a designed refusal of the package (the sweep's "No valid percentage found",
# computeRuvCovariates' "needs more than k control proteins", an invalid k), scored as arm 1's result with k NA
# (R-07); crashed = any other error, measures NA. A run whose generation or quantification fails is crashed in
# every arm. Messages are in s2-errors.tsv.
# subset is `full`: the quantified matrix has no missing values, so the complete subset of R-13 is the full set.
# null_p05, null_p01: shares of true nulls with p below 0.05, 0.01. fdp: false calls / calls (0 with no calls).
# power: called true effects / true effects among the quantified proteins. bias: mean over true effects of
# sign(truth) * (estimate - truth), negative when estimates are shrunk towards 0. rmse: of estimate - truth over
# true effects. seconds: the arm alone (the shared steps are in s2-data.tsv).

args <- commandArgs(TRUE)
if (length(args) < 4L) {
    stop("usage: s2_limpa.R <out dir> <runs per core cell> <runs per stress cell> <cores> [cells = all]")
}
out_dir <- args[[1L]]
runs_core <- as.integer(args[[2L]])
runs_stress <- as.integer(args[[3L]])
cores <- as.integer(args[[4L]])
cells_wanted <- if (length(args) > 4L && args[[5L]] != "all") strsplit(args[[5L]], ",")[[1L]] else NULL
first_run <- as.integer(Sys.getenv("K124_FIRST_RUN", "1"))
arms_wanted <- Sys.getenv("K124_ARMS", "")
arms_wanted <- if (nzchar(arms_wanted)) as.integer(strsplit(arms_wanted, ",", fixed = TRUE)[[1L]])
stopifnot(is.null(arms_wanted) || (!anyNA(arms_wanted) && 1L %in% arms_wanted))
stopifnot(runs_core >= 0L, runs_stress >= 0L, cores >= 1L, first_run >= 1L,
    first_run + max(runs_core, runs_stress) - 1L < 10000L)

suppressMessages(pkgload::load_all(".", quiet = TRUE))
logger::log_threshold(logger::FATAL)
# Namespaces are loaded before the workers fork (their load time is then in no run's seconds). The DPC engine
# itself is first called inside a worker: its thread pool must not exist in the parent of a fork.
for (package in c("limpa", "fastlimpadpc", "limma", "ruv", "sva", "qvalue")) loadNamespace(package)

N_PROTEINS <- 1000L
PRECURSORS <- 4L
MU_RANGE <- c(2, 10)
SIGMA <- 0.6
DPC_SLOPE <- 0.7
EFFECT <- 1
RESIDUAL_VAR <- SIGMA^2 / PRECURSORS
UNIT_SD <- 0.3
FORMULA <- "~ 0 + group"
CONTRAST <- "B_vs_A=groupB-groupA"
Q_CUT <- 0.05
RUN_COLUMNS <- c("simulator", "cell", "cell_class", "run", "seed", "arm", "subset", "status", "k",
    "n_features", "n_true", "n_calls", "null_p05", "null_p01", "fdp", "power", "bias", "rmse", "seconds")
MEASURES <- c("n_features", "n_true", "n_calls", "null_p05", "null_p01", "fdp", "power", "bias", "rmse")
REFUSALS <- "No valid percentage found|needs more than [0-9]+ control proteins|value is invalid"

# --- The route's settings, read from the package and pinned (R-16) ---
shipped <- readWorkflowShippedDefaults()
RUV_PARAMETERS <- resolveDIARuvCovariatesParameters(list(), grouping = "group", ruv_mode = "covariates")
QUANT <- list(engine = resolveLimpaDpcEngine(shipped$limpaDpcQuant$engine),
    dpc_slope = as.numeric(shipped$limpaDpcQuant$dpc_slope), chunk = as.numeric(shipped$limpaDpcQuant$chunk))
stopifnot(
    identical(RUV_PARAMETERS[c("percentage_min", "percentage_max", "num_components_to_impute",
        "ruv_grouping_variable", "ruv_qval_cutoff", "ruv_fdr_method", "separation_metric", "k_penalty_weight",
        "adaptive_k_penalty")], list(percentage_min = 1, percentage_max = 30, num_components_to_impute = 5,
        ruv_grouping_variable = "group", ruv_qval_cutoff = 0.05, ruv_fdr_method = "BH",
        separation_metric = "max_difference", k_penalty_weight = 0.5, adaptive_k_penalty = TRUE)),
    identical(QUANT, list(engine = "fastlimpadpc", dpc_slope = 0.8, chunk = 1000)),
    identical(trimws(shipped$deAnalysisParameters$formula_string), FORMULA),
    as.numeric(shipped$deAnalysisParameters$treat_lfc_cutoff) == 0,
    as.numeric(shipped$deAnalysisParameters$de_q_val_thresh) == Q_CUT,
    tolower(shipped$deAnalysisParameters$eBayes_trend) == "true",
    tolower(shipped$deAnalysisParameters$eBayes_robust) == "true"
)

quietly <- function(expr) {
    suppressWarnings(suppressMessages(utils::capture.output(value <- expr, type = "output")))
    invisible(NULL)
}
clock <- \() proc.time()[[3L]]

cellTable <- function() {
    centre <- list(n = 6L, groups = 2L, effects = 0.20, up = 0.5, factors = 1L, strength = 3, rho = 0,
        missing = 0.25, units = 0L, tails = "normal", batch = FALSE, abundance = 0)
    changes <- list(
        C0 = list(), C1 = list(n = 12L), C2 = list(effects = 0.10), C3 = list(factors = 0L),
        C4 = list(effects = 0), C5 = list(rho = 0.3), C6 = list(n = 8L, units = 4L), C7 = list(groups = 3L),
        S1 = list(n = 3L), S2 = list(effects = 0.40), S3 = list(rho = 0.6),
        S4 = list(factors = 3L, strength = 1), S5 = list(missing = 0), S6 = list(missing = 0.50),
        S7 = list(tails = "t3"), S8 = list(batch = TRUE), S9 = list(up = 0.8), S10 = list(abundance = 0.7),
        S11 = list(n = 3L, rho = 0.3), S12 = list(n = 3L, factors = 3L, strength = 1),
        S13 = list(n = 3L, factors = 0L), S14 = list(factors = 0L, missing = 0.50)
    )
    cells <- lapply(names(changes), \(name) {
        core <- startsWith(name, "C")
        number <- as.integer(substring(name, 2L))
        c(list(cell = name, class = if (core) "core" else "stress", index = if (core) number else 10L + number),
            utils::modifyList(centre, changes[[name]]))
    })
    stats::setNames(cells, names(changes))
}

armsFor <- function(cell) {
    arms <- 1:14
    if (cell$effects == 0) arms <- setdiff(arms, 9L)
    if (cell$units > 0L) arms <- c(arms, 15L, 16L)
    if (cell$batch) arms <- c(arms, 17L)
    if (!is.null(arms_wanted)) arms <- intersect(c(arms, 22L, 23L, 24L), arms_wanted)
    arms
}

# Precursor-level data of one (cell, run): limpa's generator, then the cell's terms on the complete matrix,
# then missingness redrawn with the generator's DPC coefficients, then the rows permuted.
simulateData <- function(cell, seed) {
    set.seed(seed)
    sim <- limpa::simProteinDataSet(n.peptides = N_PROTEINS * PRECURSORS, n.groups = cell$groups,
        samples.per.group = cell$n, peptides.per.protein = PRECURSORS, mu.range = MU_RANGE, prop.de = 0,
        dpc.slope = DPC_SLOPE, prop.missing = if (cell$missing > 0) cell$missing else 0.25)
    y <- sim$other$E.complete
    samples <- ncol(y)
    proteins <- unique(sim$genes$Protein)
    protein_of <- match(sim$genes$Protein, proteins)
    mu_precursor <- seq(from = MU_RANGE[[1L]], to = MU_RANGE[[2L]], length.out = nrow(y))
    if (cell$tails == "t3") y[] <- mu_precursor + stats::rt(length(y), df = 3) * SIGMA / sqrt(3)
    group_index <- sim$targets$Group
    contrast <- c(-0.5, 0.5, 0)[group_index]

    n_effects <- round(cell$effects * N_PROTEINS)
    effect_rows <- sample.int(N_PROTEINS, n_effects)
    n_up <- ceiling(cell$up * n_effects)
    delta <- numeric(N_PROTEINS)
    delta[effect_rows] <- EFFECT * rep(c(1, -1), c(n_up, n_effects - n_up))
    term <- outer(delta, contrast)

    standardise <- \(x) (x - mean(x)) / sqrt(mean((x - mean(x))^2))
    scores <- matrix(numeric(0), samples, 0L)
    if (cell$factors > 0L) {
        rho <- c(cell$rho, rep(0, cell$factors - 1L))
        noise <- matrix(stats::rnorm(samples * cell$factors), samples)
        scores <- sweep(noise, 2L, sqrt(1 - rho^2), "*") + outer(standardise(contrast), rho)
        colnames(scores) <- paste0("T", seq_len(cell$factors))
        mu_protein <- standardise(as.numeric(tapply(mu_precursor, protein_of, mean)))
        loadings <- sqrt(cell$strength * RESIDUAL_VAR) * (cell$abundance * mu_protein +
            sqrt(1 - cell$abundance^2) * matrix(stats::rnorm(N_PROTEINS * cell$factors), N_PROTEINS))
        term <- term + loadings %*% t(scores)
    }
    within_group <- stats::ave(group_index, group_index, FUN = seq_along)
    batch <- NULL
    if (cell$batch) {
        first <- ifelse(group_index == 1L, ceiling(2 * cell$n / 3), floor(cell$n / 3))
        batch <- as.numeric(within_group > first)
        term <- term + outer(stats::rnorm(N_PROTEINS, 0, sqrt(cell$strength * RESIDUAL_VAR)), standardise(batch))
    }
    replicates <- within_group
    if (cell$units > 0L) {
        replicates <- ceiling(within_group / (cell$n / cell$units))
        unit <- (group_index - 1L) * cell$units + replicates
        term <- term + matrix(stats::rnorm(N_PROTEINS * max(unit), 0, UNIT_SD), N_PROTEINS)[, unit]
    }
    y <- y + term[protein_of, , drop = FALSE]
    if (cell$missing > 0) {
        probability <- stats::plogis(sim$dpc[[1L]] + y * sim$dpc[[2L]], lower.tail = FALSE)
        y[stats::rbinom(length(y), prob = probability, size = 1L) == 1L] <- NA
    }
    order <- sample.int(nrow(y))
    runs <- sprintf("run_%03d", seq_len(samples))
    design <- data.frame(Run = runs, group = c("A", "B", "C")[group_index], replicates = replicates,
        stringsAsFactors = FALSE)
    design$tech_rep_group <- paste(design$group, design$replicates, sep = "_")
    if (cell$batch) design$batch <- batch
    rownames(scores) <- runs
    list(y = y[order, , drop = FALSE], protein = sim$genes$Protein[order], precursor = rownames(y)[order],
        design = design, truth = stats::setNames(delta, proteins), scores = scores,
        dpc = sim$dpc, factor_group_cor = if (cell$factors > 0L) stats::cor(scores[, 1L], contrast) else NA_real_)
}

# The log2 peptide state as the route's bounded DPC worker holds it: the matrix, the design, the distinct
# identities and the per-protein counts. Slots are set as that worker sets them (no long table exists here).
quantify <- function(sim) {
    y <- sim$y
    dimnames(y) <- list(paste(sim$protein, sim$precursor, sep = "%"), sim$design$Run)
    object <- methods::new("PeptideQuantitativeData")
    object@peptide_data <- data.frame(Protein.Ids = sim$protein, Stripped.Sequence = sim$precursor,
        stringsAsFactors = FALSE)
    object@peptide_matrix <- y
    object@protein_id_column <- "Protein.Ids"
    object@peptide_sequence_column <- "Stripped.Sequence"
    object@is_logged_data <- TRUE
    object@design_matrix <- sim$design
    object@sample_id <- "Run"
    object@group_id <- "group"
    object@technical_replicate_id <- "tech_rep_group"
    object@args <- list(globalParameters = list(workflow_type = "DIA_LIMPA"))
    counts <- tibble::tibble(Protein.Ids = unique(sim$protein), peptide_count = PRECURSORS,
        peptidoform_count = PRECURSORS)
    quantified <- NULL
    quietly(quantified <- proteinMissingValueImputationLimpa(object, dpc_slope = QUANT$dpc_slope,
        verbose = FALSE, chunk = QUANT$chunk, engine = QUANT$engine, peptide_summary = counts))
    quantified
}

# The package's q-values for a vector of p-values (qvalue, Benjamini-Hochberg when qvalue cannot run).
packageQ <- function(p) {
    adjustProtDaResultTable(data.frame(P.Value = p), "raw_pvalue", "fdr_qvalue", "fdr_value_bh_adjustment",
        messageFn = \(...) invisible(NULL))$table$fdr_qvalue
}

# The route's test: the dpcDE branch of differentialAbundanceAnalysisHelper. `covariates` (samples in rows)
# join the design table and the formula, as the W columns do there.
routeTest <- function(elist, design, covariates = NULL) {
    formula_used <- FORMULA
    if (!is.null(covariates)) {
        design <- cbind(design[setdiff(names(design), colnames(covariates))],
            as.data.frame(covariates)[rownames(design), , drop = FALSE])
        formula_used <- paste(FORMULA, "+", paste(colnames(covariates), collapse = " + "))
    }
    model <- stats::model.matrix(stats::as.formula(formula_used), design)
    samples <- rownames(model)
    replicate <- resolveProtTechnicalReplicateBlock(design[samples, , drop = FALSE], "group", "tech_rep_group")
    block <- if (isTRUE(replicate$technical)) replicate$block else NULL
    fit <- NULL
    result <- NULL
    quietly(fit <- limpa::dpcDE(elist[, samples], model, plot = FALSE, block = block))
    quietly(result <- convertDpcDAToStandardFormat(dpc_fit = fit, contrast_strings = CONTRAST,
        design_matrix = model, eBayes_trend = TRUE, eBayes_robust = TRUE))
    table <- result$results[[1L]]
    data.frame(id = rownames(table), logFC = table$logFC, p = table$raw_pvalue, q = table$fdr_qvalue,
        stringsAsFactors = FALSE)
}

measure <- function(table, truth) {
    delta <- truth[table$id]
    null <- delta == 0
    called <- !is.na(table$q) & table$q < Q_CUT
    n_true <- sum(!null)
    error <- table$logFC[!null] - delta[!null]
    list(n_features = nrow(table), n_true = n_true, n_calls = sum(called),
        null_p05 = mean(table$p[null] < 0.05), null_p01 = mean(table$p[null] < 0.01),
        fdp = if (any(called)) sum(called & null) / sum(called) else 0,
        power = if (n_true > 0L) sum(called & !null) / n_true else NA_real_,
        bias = if (n_true > 0L) mean(sign(delta[!null]) * error) else NA_real_,
        rmse = if (n_true > 0L) sqrt(mean(error^2)) else NA_real_)
}

runRow <- function(cell, run, seed, arm, status, k, measures, seconds) {
    if (is.null(measures)) measures <- stats::setNames(as.list(rep(NA_real_, length(MEASURES))), MEASURES)
    data.frame(c(list(simulator = "S2", cell = cell$cell, cell_class = cell$class, run = run, seed = seed,
        arm = arm, subset = "full", status = status, k = k), measures, list(seconds = round(seconds, 3))),
        stringsAsFactors = FALSE)[RUN_COLUMNS]
}

errorRow <- function(cell, run, arm, status, condition) {
    data.frame(cell = cell$cell, run = run, arm = arm, status = status,
        message = gsub("[\t\r\n]+", " ", conditionMessage(condition)), stringsAsFactors = FALSE)
}

crashedRun <- function(cell, run, seed, condition, data_row = NULL) {
    arms <- armsFor(cell)
    list(runs = do.call(rbind, lapply(arms, \(arm) runRow(cell, run, seed, arm, "crashed", NA, NULL, NA_real_))),
        data = data_row %||% dataRow(cell, run, seed, list(status = "crashed")),
        errors = errorRow(cell, run, 0L, "crashed", condition))
}

dataRow <- function(cell, run, seed, values) {
    fields <- list(status = "ok", n_samples = NA, n_precursors = NA, missing_share = NA, n_proteins = NA,
        dpc_intercept_set = NA, dpc_slope_set = NA, dpc_intercept_fit = NA, dpc_slope_fit = NA,
        factor_group_cor = NA, residual_var_median = NA, k_e = NA, sweep_status = NA, sweep_k = NA,
        sweep_percentage = NA, sweep_controls = NA, replicate_structure = NA, generate_seconds = NA,
        quant_seconds = NA, numsv_seconds = NA, sweep_seconds = NA, total_seconds = NA)
    stopifnot(all(names(values) %in% names(fields)))
    values <- utils::modifyList(fields, values)
    timed <- grep("_seconds$", names(values), value = TRUE)
    values[timed] <- lapply(values[timed], round, digits = 3)
    data.frame(c(list(simulator = "S2", cell = cell$cell, cell_class = cell$class, run = run, seed = seed),
        values), stringsAsFactors = FALSE)
}

runOne <- function(cell, run) {
    seed <- 100000L + 10000L * cell$index + run
    started <- clock()
    shared <- list()
    prepared <- tryCatch({
        sim <- simulateData(cell, seed)
        shared$generate_seconds <- clock() - started
        mark <- clock()
        quantified <- quantify(sim)
        shared$quant_seconds <- clock() - mark
        results <- quantified@args$limpa_dpc_quant_results
        elist <- results$quantified_elist
        stopifnot(!anyNA(elist$E), identical(colnames(elist$E), sim$design$Run),
            identical(rownames(elist$E), quantified@protein_quant_table[["Protein.Ids"]]))
        design <- sim$design
        rownames(design) <- design$Run
        true_terms <- sim$scores
        if (cell$batch) true_terms <- cbind(true_terms, batch = design$batch)
        truth_model <- cbind(stats::model.matrix(~ 0 + group, design), true_terms)
        mark <- clock()
        k_e <- NULL
        quietly(k_e <- sva::num.sv(elist$E, stats::model.matrix(~ group, design), method = "be", B = 20,
            seed = 100L * seed + 99L))
        shared$numsv_seconds <- clock() - mark
        shared <- c(shared, list(n_samples = ncol(sim$y), n_precursors = nrow(sim$y),
            missing_share = mean(is.na(sim$y)), n_proteins = nrow(elist$E),
            dpc_intercept_set = if (cell$missing > 0) sim$dpc[[1L]] else NA_real_, dpc_slope_set = sim$dpc[[2L]],
            dpc_intercept_fit = results$dpc_parameters_used[[1L]], dpc_slope_fit = results$dpc_parameters_used[[2L]],
            factor_group_cor = sim$factor_group_cor,
            residual_var_median = stats::median(limma::lmFit(elist$E, truth_model)$sigma^2),
            k_e = as.integer(k_e)))
        TRUE
    }, error = \(e) e)
    if (inherits(prepared, "condition")) return(crashedRun(cell, run, seed, prepared))

    ids <- rownames(elist$E)
    truth <- sim$truth
    null_ids <- ids[truth[ids] == 0]
    effect_ids <- ids[truth[ids] != 0]
    all_controls <- rep(TRUE, length(ids))
    set.seed(100L * seed + 8L)
    oracle <- ids %in% sample(null_ids, 50L)
    set.seed(100L * seed + 9L)
    degraded <- if (9L %in% armsFor(cell)) {
        ids %in% c(sample(null_ids, 40L), sample(effect_ids, 10L))
    }

    # The route's sweep, once per run on arm 3's stream; arms 3 to 6 take its controls or its k.
    sweeping <- any(3:6 %in% armsFor(cell))
    set.seed(100L * seed + 3L)
    mark <- clock()
    sweep <- if (!sweeping) simpleError("the sweep was not run (K124_ARMS)") else tryCatch({
        value <- NULL
        quietly(value <- findBestNegCtrlPercentageCapped(quantified,
            percentage_range = seq(RUV_PARAMETERS$percentage_min, RUV_PARAMETERS$percentage_max, by = 1),
            num_components_to_impute = RUV_PARAMETERS$num_components_to_impute,
            ruv_grouping_variable = RUV_PARAMETERS$ruv_grouping_variable,
            ruv_qval_cutoff = RUV_PARAMETERS$ruv_qval_cutoff, ruv_fdr_method = RUV_PARAMETERS$ruv_fdr_method,
            separation_metric = RUV_PARAMETERS$separation_metric,
            k_penalty_weight = RUV_PARAMETERS$k_penalty_weight,
            adaptive_k_penalty = isTRUE(RUV_PARAMETERS$adaptive_k_penalty), verbose = FALSE))
        value
    }, error = \(e) e)
    shared$sweep_seconds <- clock() - mark
    swept <- !inherits(sweep, "condition")
    shared$sweep_status <- if (!sweeping) "not run" else if (swept) "ok" else if (grepl(REFUSALS, conditionMessage(sweep))) "refused" else "crashed"
    if (swept) {
        shared <- c(shared, list(sweep_k = sweep$best_k, sweep_percentage = sweep$best_percentage,
            sweep_controls = sum(sweep$best_control_genes_index, na.rm = TRUE)))
    }
    today <- \() if (swept) sweep else stop(sweep)

    factorsOf <- function(k, controls) {
        value <- NULL
        quietly(value <- computeRuvCovariates(quantified, k, controls, RUV_PARAMETERS$ruv_grouping_variable))
        shared$replicate_structure <<- value$replicate_structure
        value
    }
    inModel <- \(k, controls) {
        factors <- factorsOf(k, controls)
        list(k = factors$k, table = routeTest(elist, design, factors$W))
    }
    subtracted <- \(k, controls) {
        factors <- factorsOf(k, controls)
        corrected <- elist
        corrected$E <- factors$adjusted[rownames(elist$E), colnames(elist$E), drop = FALSE]
        list(k = factors$k, table = routeTest(corrected, design))
    }
    separationK <- function(controls) {
        plot <- NULL
        quietly(plot <- ruvCancorFromInput(prepareRuvCancorInput(quantified, RUV_PARAMETERS$num_components_to_impute),
            quantified, controls, RUV_PARAMETERS$ruv_grouping_variable))
        curve <- plot$data[plot$data$K <= min(calculateAdaptiveMaxK(ncol(elist$E)), sum(controls) - 1), , drop = FALSE]
        of_controls <- curve$featureset == "Control"
        difference <- curve$cc[curve$featureset == "All"] - curve$cc[of_controls]
        as.integer(curve$K[of_controls][[which(difference == max(difference, na.rm = TRUE))[[1L]]]])
    }
    first <- NULL
    noRuv <- \() list(k = 0L, table = first)
    indicator <- \(level) matrix(as.numeric(design$group == level), ncol = 1L)
    other_groups <- if (cell$groups > 2L) indicator("C") else 1
    ruvTable <- \(fit) {
        p <- ruv::variance_adjust(fit, ebayes = TRUE, pooled = FALSE, evar = FALSE, rsvar = FALSE)$p.ebayes[1L, ]
        data.frame(id = ids, logFC = fit$betahat[1L, ], p = p, q = packageQ(p), stringsAsFactors = FALSE)
    }
    arm_functions <- list(
        `1` = \() list(k = 0L, table = routeTest(elist, design)),
        `2` = \() if (ncol(true_terms)) list(k = cell$factors, table = routeTest(elist, design, true_terms)) else noRuv(),
        `3` = \() subtracted(today()$best_k, today()$best_control_genes_index),
        `4` = \() inModel(today()$best_k, today()$best_control_genes_index),
        `5` = \() inModel(today()$best_k, rank(first$p, ties.method = "first") > 0.25 * nrow(first)),
        `6` = \() inModel(today()$best_k, !is.na(first$q) & first$q > 0.5),
        `7` = \() if (k_e > 0L) inModel(k_e, all_controls) else noRuv(),
        `8` = \() if (k_e > 0L) inModel(k_e, oracle) else noRuv(),
        `9` = \() if (k_e > 0L) inModel(k_e, degraded) else noRuv(),
        `10` = \() if (k_e > 0L) subtracted(k_e, oracle) else noRuv(),
        `11` = \() inModel(1L, all_controls),
        `12` = \() list(k = k_e, table = ruvTable(ruv::RUV4(t(elist$E), indicator("B"), oracle, k = k_e,
            Z = other_groups))),
        `13` = \() list(k = NA_integer_, table = ruvTable(ruv::RUVinv(t(elist$E), indicator("B"), oracle,
            Z = other_groups))),
        `14` = \() {
            model <- stats::model.matrix(stats::as.formula(FORMULA), design)
            if (k_e > 0L) {
                surrogate <- NULL
                quietly(surrogate <- sva::sva(elist$E, stats::model.matrix(~ group, design),
                    stats::model.matrix(~ 1, design), n.sv = k_e))
                model <- cbind(model, matrix(surrogate$sv, ncol = k_e, dimnames = list(NULL, paste0("SV", seq_len(k_e)))))
            }
            contrasts <- limma::makeContrasts(contrasts = "groupB-groupA", levels = colnames(model))
            fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(elist$E, model), contrasts), trend = TRUE,
                robust = TRUE)
            list(k = k_e, table = data.frame(id = ids, logFC = fit$coefficients[, 1L], p = fit$p.value[, 1L],
                q = packageQ(fit$p.value[, 1L]), stringsAsFactors = FALSE))
        },
        `15` = \() if (k_e > 0L) subtracted(k_e, all_controls) else noRuv(),
        `16` = \() if (k_e > 0L) inModel(k_e, all_controls) else noRuv(),
        `17` = \() list(k = 0L, table = routeTest(elist, design, design["batch"])),
        `22` = \() inModel(separationK(oracle), oracle),
        `23` = \() if (k_e > 0L) inModel(k_e, rank(first$p, ties.method = "first") > 0.25 * nrow(first)) else noRuv(),
        `24` = \() {
            wide <- rank(first$p, ties.method = "first") > 0.25 * nrow(first)
            inModel(separationK(wide), wide)
        }
    )

    arms <- armsFor(cell)
    rows <- vector("list", length(arms))
    errors <- list()
    first_measures <- NULL
    for (i in seq_along(arms)) {
        arm <- arms[[i]]
        set.seed(100L * seed + arm)
        mark <- clock()
        outcome <- tryCatch(arm_functions[[as.character(arm)]](), error = \(e) e)
        seconds <- clock() - mark
        if (inherits(outcome, "condition")) {
            refused <- arm != 1L && !is.null(first_measures) && grepl(REFUSALS, conditionMessage(outcome))
            status <- if (refused) "refused" else "crashed"
            errors[[length(errors) + 1L]] <- errorRow(cell, run, arm, status, outcome)
            rows[[i]] <- runRow(cell, run, seed, arm, status, NA, if (refused) first_measures, seconds)
            next
        }
        measures <- measure(outcome$table, truth)
        if (arm == 1L) {
            stopifnot(identical(outcome$table$id, ids))
            first <- outcome$table
            first_measures <- measures
        }
        rows[[i]] <- runRow(cell, run, seed, arm, "ok", outcome$k, measures, seconds)
    }
    shared$total_seconds <- clock() - started
    list(runs = do.call(rbind, rows), data = dataRow(cell, run, seed, shared),
        errors = if (length(errors)) do.call(rbind, errors))
}

writeTsv <- function(table, path) {
    utils::write.table(table, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

scriptPath <- function() {
    file_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_argument)) sub("^--file=", "", file_argument[[1L]]) else NA_character_
}

writeMeta <- function(path, started, finished = NULL, cells_run = character()) {
    git <- \(...) tryCatch(paste(system2("git", c(...), stdout = TRUE, stderr = FALSE), collapse = " "),
        error = \(e) NA_character_, warning = \(w) NA_character_)
    script <- scriptPath()
    packages <- c("MultiScholaR", "limpa", "fastlimpadpc", "limma", "ruv", "sva", "qvalue", "statmod", "pkgload")
    versions <- vapply(packages, \(p) as.character(utils::packageVersion(p)), character(1))
    lines <- c(
        "simulator: S2 (DIA_LIMPA route), KIT-124 revision 1.2",
        paste("commit:", git("rev-parse", "HEAD")),
        paste("branch:", git("rev-parse", "--abbrev-ref", "HEAD")),
        paste("script:", script),
        paste("script md5:", if (!is.na(script)) unname(tools::md5sum(script)) else NA),
        paste("script differs from commit:", nzchar(git("status", "--porcelain", "--", script))),
        paste("changes under R/, tests/, tools/:", nzchar(git("status", "--porcelain", "--", "R", "tests", "tools"))),
        paste("command: Rscript", paste(c(script, commandArgs(TRUE)), collapse = " ")),
        paste("K124_FIRST_RUN:", first_run),
        paste("K124_ARMS:", if (is.null(arms_wanted)) "unset (arms 1 to 17 as the cell has them)" else paste(arms_wanted, collapse = ",")),
        paste("working directory:", getwd()),
        paste("R_LIBS:", Sys.getenv("R_LIBS")),
        paste("library paths:", paste(.libPaths(), collapse = ":")),
        paste("R:", R.version.string),
        paste("platform:", R.version$platform),
        paste("BLAS:", extSoftVersion()[["BLAS"]]),
        paste("packages:", paste(names(versions), versions, collapse = "; ")),
        paste("DPC engine:", QUANT$engine, "(RAYON_NUM_THREADS set to 1 by the package inside the method)"),
        paste("cores:", cores),
        paste("runs per core cell:", runs_core, "; runs per stress cell:", runs_stress),
        "seed rule: seed = 100000 + 10000 * cell index + run (C0..C7 = 0..7, S1..S14 = 11..24);",
        "  data set.seed(seed); arm a set.seed(100 * seed + a); sweep on arm 3's stream; E seed 100 * seed + 99",
        paste("cells completed in this invocation:", paste(cells_run, collapse = ",")),
        paste("started:", format(started, "%Y-%m-%d %H:%M:%S %Z")),
        paste("finished:", if (is.null(finished)) "not yet" else format(finished, "%Y-%m-%d %H:%M:%S %Z"))
    )
    writeLines(lines, path)
}

# --- Main: one cell at a time, each written when complete ---
cells <- cellTable()
if (!is.null(cells_wanted)) {
    stopifnot(all(cells_wanted %in% names(cells)))
    cells <- cells[cells_wanted]
}
cell_dir <- file.path(out_dir, "cells")
dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)
meta_path <- file.path(out_dir, "s2-runs.meta.txt")
started <- Sys.time()
writeMeta(meta_path, started)
cells_run <- character()
for (cell in cells) {
    count <- if (cell$class == "core") runs_core else runs_stress
    if (count < 1L) next
    runs <- first_run + seq_len(count) - 1L
    path <- file.path(cell_dir, paste0(cell$cell, ".rds"))
    if (file.exists(path)) {
        if (!setequal(readRDS(path)$data$run, runs)) {
            stop(sprintf("%s holds other run numbers than this invocation asks for; remove it or change the runs", path))
        }
        message(sprintf("%s: reused", cell$cell))
        next
    }
    mark <- clock()
    results <- parallel::mclapply(runs, \(run) runOne(cell, run), mc.cores = cores, mc.preschedule = FALSE)
    for (i in seq_along(results)) {
        if (!is.list(results[[i]]) || is.null(results[[i]]$runs)) {
            reason <- simpleError(paste("worker failed:", paste(as.character(results[[i]]), collapse = " ")))
            results[[i]] <- crashedRun(cell, runs[[i]], 100000L + 10000L * cell$index + runs[[i]], reason)
        }
    }
    value <- list(runs = do.call(rbind, lapply(results, `[[`, "runs")),
        data = do.call(rbind, lapply(results, `[[`, "data")),
        errors = do.call(rbind, lapply(results, `[[`, "errors")))
    saveRDS(value, paste0(path, ".tmp"))
    file.rename(paste0(path, ".tmp"), path)
    cells_run <- c(cells_run, cell$cell)
    message(sprintf("%s: %d runs in %.0f s", cell$cell, count, clock() - mark))
}

done <- intersect(names(cellTable()), sub("\\.rds$", "", list.files(cell_dir, pattern = "\\.rds$")))
parts <- lapply(file.path(cell_dir, paste0(done, ".rds")), readRDS)
writeTsv(do.call(rbind, lapply(parts, `[[`, "runs")), file.path(out_dir, "s2-runs.tsv"))
writeTsv(do.call(rbind, lapply(parts, `[[`, "data")), file.path(out_dir, "s2-data.tsv"))
errors <- do.call(rbind, lapply(parts, `[[`, "errors"))
if (is.null(errors)) {
    errors <- data.frame(cell = character(), run = integer(), arm = integer(), status = character(),
        message = character())
}
writeTsv(errors, file.path(out_dir, "s2-errors.tsv"))
writeMeta(meta_path, started, Sys.time(), cells_run)
