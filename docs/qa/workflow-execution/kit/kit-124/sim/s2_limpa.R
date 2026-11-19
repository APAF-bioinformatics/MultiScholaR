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
#   C8    core    2       8 (4 units x 2)      20%      50%   0                   -    25%      unit sd 0.3; index 8
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
# The exactness round (user, 2026-11-18, decision M16; kit-124/exact-round-declaration.md), on request (K124_ARMS).
# Each arm is arm 7 with the named change and nothing else, on arm 7's random stream. Cell C8 is C6 with no unwanted
# factor (data generation as C6 with factors = 0), run when a cell argument names it; reported, judges no design.
#   25  arm 7 with the comparison's standard error exact per protein. The route's contrasts.fit forms it from one
#       correlation matrix for all proteins, which is approximate when the model carries W, because dpcDE fits with
#       per-observation weights (the quantified matrix has no missing values). Here the model matrix is
#       reparameterised so that the comparison is its first coefficient (limma::contrastAsCoef) and the route's own
#       dpcDE and convertDpcDAToStandardFormat are called with it and the contrast "B_vs_A=comparison": for a
#       contrast that is one coefficient contrasts.fit is exact.
#   29  C6 and C8: arm 7 with E estimated after the biological units are removed, not only the groups:
#       sva::num.sv(quantified matrix, model.matrix(~ tech_rep_group), method = "be", B = 20, E's seed). On this
#       route arm 7 already takes the declared technical replicates as RUV-III's replicate structure (arm 16 repeats
#       it), so this is also arm 16's variant; arm numbers 30 and 33 of S1 and S4 are not formed here.
#   31  arm 7 with E from 200 permutations: sva::num.sv(..., B = 200, E's seed), the group model.
#   32  arm 7 with all of them: arm 25's standard error, B = 200, and in C6 and C8 the units removed.
# The confirmation round (KIT-124 / decision M29; kit-124/scored/confirm-round-declaration.md), on request (K124_ARMS),
# with K124_PACKAGE=installed (the MultiScholaR installed in the job from the bound commit). The RUV step and the test
# are the package's: the declared covariates kernel, as kit-129/final-method/prove-arms.R calls it,
# diaRuvCovariatesOperation("declared")$kernel(<the quantified object>, resolveDIARuvCovariatesParameters(list(),
# "group", ruv_control_rule = <the arm's>, ruv_k_rule = <the arm's>, ruv_control_ids = <the arm's>, ruv_k_seed =
# 100 * seed + 99), <a store that keeps the checkpoints>), then differentialAbundanceAnalysis(<its final state>, the
# default path's settings). The kernel's own no-RUV (a count of 0, or what the data cannot give under an automatic rule)
# is its final state without W; a run the kernel leaves without RUV for a reason other than a count of 0 is `refused`
# (scored as arm 1). The package's classed refusals (ruv_k_refused, ruv_factor_count_unavailable,
# omics_ruv_no_residual_df, omics_da_block_correlation_failed, omics_da_no_biological_replication, and the kernel's
# workflow_execution_invalid for a control list that matches nothing) are refusals; any other error is a crash.
#   34  recommended: all proteins as controls, the estimated count (simulator arm 32; in C6 and C8 the kernel takes
#       the declared technical replicates as replicate structure, so this is also the technical-replicate design)
#   35  supplied controls: the oracle set (arm 8's), the estimated count
#   36  supplied controls: the degraded set (arm 9's), the estimated count                 every cell but C4
#   37  no RUV: differentialAbundanceAnalysis on the quantified object (simulator arm 1)
#   39  first-pass nulls (the kernel's own first pass: dpcDE, the F test over the groups), the estimated count
#   40  first-pass nulls, the separation rule (the kernel's default for this control rule: today's sweep's k)
#   41  supplied controls (arm 8's oracle set), the separation rule: from the KIT-129 sixth pass the kernel runs the rule
#       on the list (diaRuvSeparationRuleK, arm 22's separationK); before it, today's sweep's k whatever the control rule
# Arm 38 (technical replicates as their own design) is not formed: on this route it is arm 34's computation.
# K124_FAULT=<arm> (environment, the lifecycle test only): that arm stops with an unclassed error in the first run of
# the invocation, so that a crash of a shipped-kernel arm can be seen to be reported.
# Arm numbers 26 to 28 (S1: the test applies its trend and robust settings) are not formed: the dpcDE branch has no
# TREAT step and its eBayes(trend, robust) p-values are the ones reported. With E = 0 arms 29 and 31 are arm 1's
# result, as arm 7; arms 25 and 32 are the no-RUV fit with the reparameterised matrix (equal to arm 1's up to
# rounding: without W the model's coefficients are uncorrelated). A refusal is scored as arm 1's result, as before.
# K124_EXACT_CHECK=1 (environment; a check, not a result): per data set the reparameterised no-RUV fit must give
# arm 1's estimates and p-values to 1e-8; otherwise arm 25 crashes with the reason.
# K124_ARMS (environment; unset = the arms above, 22 to 32 not run): a comma list of the arms to run, of those the
# cell has and 22 to 25, 29 (C6, C8), 31, 32. Arm 1 must be in it (arms 5, 6, 23, 24 and every refusal read its result). Without any of arms
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

installed <- identical(Sys.getenv("K124_PACKAGE"), "installed")
if (installed) {
    # the confirmation round: the installed package, every object of its namespace on the search path as load_all puts it
    suppressPackageStartupMessages(library(MultiScholaR))
    attach(asNamespace("MultiScholaR"), name = "MultiScholaR:namespace", warn.conflicts = FALSE)
} else suppressMessages(pkgload::load_all(".", quiet = TRUE))
SHIPPED_ARMS <- c(34L, 35L, 36L, 37L, 39L, 40L, 41L)
if (!installed && any(SHIPPED_ARMS %in% arms_wanted)) stop("arms 34 to 41 call the installed package: set K124_PACKAGE=installed")
fault_arm <- suppressWarnings(as.integer(Sys.getenv("K124_FAULT", "")))
# The confirmation round's savings (declaration section 8; K124_SHARE, a comma list, unset = none): one call per data
# set of a package function the shipped arms each make with the same arguments. "count" shares estimateRuvFactorCount
# (L1), "sweep" findBestNegCtrlPercentageCapped (L2, S2's arms 40 and 41). The package's own function is called once;
# a later call with identical arguments (digest of the whole argument list) returns that value, and leaves the random
# number state as the call itself did (unchanged when the function restores it, else the state it left). The binding
# is replaced in the installed namespace, so the kernels call it as before; the cache is emptied at the start of every
# run (shareReset). Every shipped arm's rows must be identical with and without it (the declaration's proof).
share_wanted <- intersect(strsplit(Sys.getenv("K124_SHARE", ""), ",", fixed = TRUE)[[1L]], c("count", "sweep"))
if (length(share_wanted) && !installed) stop("K124_SHARE applies to the installed package only")
shared_cache <- new.env()
shareReset <- function() rm(list = ls(shared_cache, all.names = TRUE), envir = shared_cache)
for (shared_name in c(count = "estimateRuvFactorCount", sweep = "findBestNegCtrlPercentageCapped")[share_wanted]) local({
    name <- shared_name
    real <- get(name, envir = asNamespace("MultiScholaR"), inherits = FALSE)
    shared <- function(...) {
        key <- paste(name, digest::digest(list(...), algo = "sha256"))
        hit <- shared_cache[[key]]
        if (!is.null(hit)) {
            if (hit$changed) assign(".Random.seed", hit$after, envir = globalenv())
            return(hit$value)
        }
        before <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
        value <- real(...)
        after <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
        assign(key, list(value = value, changed = !identical(before, after) && !is.null(after), after = after), envir = shared_cache)
        value
    }
    utils::assignInNamespace(name, shared, ns = "MultiScholaR")
})
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
        S13 = list(n = 3L, factors = 0L), S14 = list(factors = 0L, missing = 0.50),
        C8 = list(n = 8L, units = 4L, factors = 0L)
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
    if (!is.null(arms_wanted)) arms <- intersect(c(arms, 22L, 23L, 24L, 25L, if (cell$units > 0L) 29L, 31L, 32L,
        if (installed) setdiff(SHIPPED_ARMS, if (cell$effects == 0) 36L)), arms_wanted)
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
routeTest <- function(elist, design, covariates = NULL, exact = FALSE) {
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
    contrast_string <- CONTRAST
    if (exact) {                                          # the exactness round: the comparison as the first coefficient
        contrast <- limma::makeContrasts(contrasts = sub("^.*=", "", CONTRAST), levels = colnames(model))
        model <- limma::contrastAsCoef(model, contrast)$design
        dimnames(model) <- list(samples, c("comparison", paste0("Q", seq_len(ncol(model))[-1L])))
        contrast_string <- sub("=.*$", "=comparison", CONTRAST)
    }
    fit <- NULL
    result <- NULL
    quietly(fit <- limpa::dpcDE(elist[, samples], model, plot = FALSE, block = block))
    quietly(result <- convertDpcDAToStandardFormat(dpc_fit = fit, contrast_strings = contrast_string,
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
    shareReset()
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
    degraded <- if (any(c(9L, 36L) %in% armsFor(cell))) {
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
    inModel <- \(k, controls, exact = FALSE) {
        factors <- factorsOf(k, controls)
        list(k = factors$k, table = routeTest(elist, design, factors$W, exact = exact))
    }
    # the exactness round: arm 7 with the changes named, on arm 7's stream (E first: num.sv sets a seed)
    counts <- list()
    seventh <- NULL
    variant <- \(exact = FALSE, B = 20L, units = FALSE) \() {
        k <- k_e
        key <- paste(B, units)
        if (B != 20L || units) {
            if (is.null(counts[[key]])) {
                model <- if (units) stats::model.matrix(~ tech_rep_group, design) else stats::model.matrix(~ group, design)
                quietly(k <- sva::num.sv(elist$E, model, method = "be", B = B, seed = 100L * seed + 99L))
                counts[[key]] <<- as.integer(k)
            }
            k <- counts[[key]]
        }
        # only E changed and it equals arm 7's k: arm 7's result, not computed again
        if (!exact && !is.null(seventh) && isTRUE(seventh$k == k)) return(seventh)
        set.seed(100L * seed + 7L)
        if (exact && identical(Sys.getenv("K124_EXACT_CHECK"), "1")) {
            none <- routeTest(elist, design, exact = TRUE)
            gap <- max(abs(none$logFC - first$logFC), abs(none$p - first$p))
            if (!identical(none$id, first$id) || !is.finite(gap) || gap > 1e-8) stop(sprintf("the reparameterised no-RUV fit is not arm 1's: differs by up to %.3g", gap))
        }
        if (k > 0L) inModel(k, all_controls, exact = exact) else if (exact) list(k = 0L, table = routeTest(elist, design, exact = TRUE)) else noRuv()
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
    # the confirmation round: the installed package's declared covariates kernel and DA entry point (arms 34 to 40)
    shippedDa <- function(state) {
        value <- NULL
        quietly(value <- differentialAbundanceAnalysis(state, contrasts_tbl = data.frame(contrasts = CONTRAST,
            stringsAsFactors = FALSE), formula_string = FORMULA, da_q_val_thresh = Q_CUT, treat_lfc_cutoff = 0,
            eBayes_trend = TRUE, eBayes_robust = TRUE, args_row_id = "Protein.Ids", qvalue_column = "fdr_qvalue",
            raw_pvalue_column = "raw_pvalue"))
        table <- value$contrasts_results$results[[1L]]
        if (is.null(table) || !nrow(table)) stop("differentialAbundanceAnalysis returned no table")
        data.frame(id = rownames(table), logFC = table$logFC, p = table$raw_pvalue, q = table$fdr_qvalue,
            stringsAsFactors = FALSE)
    }
    shipped <- function(...) {
        saved <- list()
        store <- list(checkpoint = function(value, label, from = NULL) {
            saved[[label]] <<- value
            label
        })
        parameters <- resolveDIARuvCovariatesParameters(list(), "group", ruv_k_seed = 100L * seed + 99L, ...)
        value <- NULL
        quietly(value <- diaRuvCovariatesOperation("declared")$kernel(quantified, parameters, store))
        state <- saved[[value$final]]
        if (is.null(state)) stop("the declared kernel's final state was not checkpointed")
        refused <- !isTRUE(value$ruv_applied) && !isTRUE(value$ruv_k_estimate == 0L)
        list(k = as.integer(value$best_k), table = shippedDa(state),
            refused = if (refused) value$ruv_not_applied_reason %||% "the kernel applied no RUV")
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
        },
        `25` = variant(exact = TRUE), `29` = variant(units = TRUE), `31` = variant(B = 200L),
        `32` = variant(exact = TRUE, B = 200L, units = cell$units > 0L),
        `34` = \() shipped(ruv_control_rule = "all_features"),
        `35` = \() shipped(ruv_control_rule = "supplied_list", ruv_control_ids = ids[oracle]),
        `36` = \() shipped(ruv_control_rule = "supplied_list", ruv_control_ids = ids[degraded]),
        `37` = \() list(k = 0L, table = shippedDa(quantified)),
        `39` = \() shipped(ruv_control_rule = "first_pass_drop_top_25", ruv_k_rule = "num_sv_be"),
        `40` = \() shipped(ruv_control_rule = "first_pass_drop_top_25", ruv_k_rule = "sweep_bounded"),
        `41` = \() shipped(ruv_control_rule = "supplied_list", ruv_control_ids = ids[oracle], ruv_k_rule = "sweep_bounded")
    )
    SHIPPED_REFUSALS <- c("ruv_k_refused", "ruv_factor_count_unavailable", "omics_ruv_no_residual_df",
        "omics_da_block_correlation_failed", "omics_da_no_biological_replication")

    arms <- armsFor(cell)
    rows <- vector("list", length(arms))
    errors <- list()
    first_measures <- NULL
    for (i in seq_along(arms)) {
        arm <- arms[[i]]
        set.seed(100L * seed + arm)
        mark <- clock()
        outcome <- tryCatch({
            if (identical(arm, fault_arm) && run == first_run) stop("deliberate fault (K124_FAULT, the lifecycle test)")
            arm_functions[[as.character(arm)]]()
        }, error = \(e) e)
        seconds <- clock() - mark
        if (arm %in% SHIPPED_ARMS && !inherits(outcome, "condition") && !is.null(outcome$refused)) {
            # the kernel left the run without RUV for a reason other than a count of 0: a refusal, scored as arm 1
            outcome <- structure(class = c("s2_shipped_refusal", "error", "condition"),
                list(message = paste("the kernel applied no RUV:", outcome$refused), call = NULL))
        }
        if (inherits(outcome, "condition")) {
            refused <- arm != 1L && !is.null(first_measures) && (grepl(REFUSALS, conditionMessage(outcome)) ||
                (arm %in% SHIPPED_ARMS && inherits(outcome, c("s2_shipped_refusal", SHIPPED_REFUSALS))))
            status <- if (refused) "refused" else "crashed"
            errors[[length(errors) + 1L]] <- errorRow(cell, run, arm, status, outcome)
            rows[[i]] <- runRow(cell, run, seed, arm, status, NA, if (refused) first_measures, seconds)
            next
        }
        if (arm == 7L) seventh <- outcome
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
        paste("MultiScholaR:", if (installed) paste0("installed (K124_PACKAGE=installed; built from ",
            Sys.getenv("K124_PACKAGE_COMMIT", "unknown"), ", at ", find.package("MultiScholaR"), ")") else "pkgload::load_all of the work tree"),
        paste("K124_FAULT:", if (is.na(fault_arm)) "unset" else fault_arm),
        paste("K124_SHARE:", if (length(share_wanted)) paste(share_wanted, collapse = ",") else "unset"),
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
        "seed rule: seed = 100000 + 10000 * cell index + run (C0..C8 = 0..8, S1..S14 = 11..24);",
        "  data set.seed(seed); arm a set.seed(100 * seed + a) (arms 25 to 32: arm 7's stream); sweep on arm 3's stream; E seed 100 * seed + 99",
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
} else {
    cells <- cells[setdiff(names(cells), "C8")]           # the exactness round's cell runs on request only
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
