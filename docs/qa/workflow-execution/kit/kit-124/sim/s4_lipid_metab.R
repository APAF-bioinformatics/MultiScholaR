# KIT-124 simulator S4: lipidomics and metabolomics (ticket revision 1.2, "Simulation Setup" and W3).
#
# Usage (from the repository root, R 4.6.1, OMP_NUM_THREADS=1 so that <cores> is the number of cores used):
#   OMP_NUM_THREADS=1 Rscript docs/qa/workflow-execution/kit/kit-124/sim/s4_lipid_metab.R \
#       <out dir> <runs per core cell> <runs per stress cell> <cores> [cells = all] [first run = 1]
#   cells: "all" or a comma list (C0,...,C7,S1,...,S14). Runs are numbered first run, first run + 1, ...
# Writes <out dir>/s4-runs.tsv (one row per cell, run, arm and subset; the two assays and all contrasts pooled),
#   s4-runs.meta.txt (commit, libraries, versions, command, seed rule, the lipid/metab sharing check),
#   s4-notes.tsv (every refusal and crash with its reason, apart from the measures) and s4-timing.tsv (seconds per run).
# Restartable per cell: a finished cell is kept in <out dir>/cells/ and skipped on the next call.
# Libraries: R_LIBS as the ticket names them; the package is loaded with pkgload::load_all(".").
#
# ---------------------------------------------------------------------------------------------------------------
# 1. What of the package is run, and the significance rule
# ---------------------------------------------------------------------------------------------------------------
# Default path of a GUI run, established by reading (R-16) and pinned at start-up by `checkSharedKernels()`:
#   sweep   runLipidPerAssayRuvOptimization(theObject, ruv_mode = "automatic", params = list(ruv_grouping_variable =
#           "group", percentage_min = 1, percentage_max = 20, separation_metric = "max_difference",
#           k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, manual_k = 2, manual_percentage = 5))
#           (the tab's defaults, R/mod_omic_norm_panels.R:129-146; the call of omicsRuvOperation,
#           R/utils_workflow_omics_declarations.R:305-314). Metab: runPerAssayRuvOptimization, same arguments.
#   RUV     ruvIII_C_Varying(theObject, ruv_grouping_variable = "group", ruv_number_k = extractLipidBestKPerAssay(.),
#           ctrl = extractLipidCtrlPerAssay(.)): RUVIIIC::RUVIII_C_Varying per assay, the group as the replicate
#           structure, levels not restored (convention.md 1.5).
#   test    runLipidsDA(theObject, contrasts_tbl, formula_string = "~ 0 + group", da_q_val_thresh = 0.05,
#           treat_lfc_cutoff = 0, eBayes_trend = TRUE, eBayes_robust = TRUE) (R/utils_workflow_lipid_declarations.R:78-82).
#           Both omics: a cut-off of 0 is passed on as NA (func_lipid_da_model.R:476, func_metab_da_orchestration.R:250),
#           so there is no TREAT: limma lmFit, contrasts.fit, eBayes(trend = TRUE, robust = TRUE), topTable, per assay
#           and per contrast; q = qvalue::qvalue(p)$qvalues (BH when qvalue fails or p is degenerate).
#   rule    a call is fdr_qvalue < 0.05 (and |logFC| >= 0), per assay and contrast (formatLipidDaContrastTable).
#   ITSD    normaliseUntransformedData(theObject, method = "ITSD", itsd_aggregation = "median", itsd_feature_ids = <the
#           standards per assay>) on untransformed values, then logTransformAssays(offset = 1): the tab's defaults
#           (median is the first choice; offset 1). One factor per sample and assay; the standards are removed.
# Arm 3 calls exactly these on a LipidomicsAssayData object. Every other arm tests with the kernel the orchestrator
# calls, runTestsContrastsLipidDA(data, contrast_strings, design_matrix, formula_string, "Sample_ID",
# treat_lfc_cutoff = NA, eBayes_trend = TRUE, eBayes_robust = TRUE), and the same rule; the start-up check proves the
# kernel and the orchestrator give the same p, q and calls, and that the kernel equals limma called directly.
# Lipid and metab: the start-up check runs sweep, RUV, ITSD + log2, kernel and orchestrator through both classes on
# one data set and stops unless every output is identical. So one class (lipid) is simulated.
# Not run here (not what the arms ask for): the between-sample normalisation (cyclic loess by default), the intensity
# filter and the correlation filter of the GUI route.
#
# ---------------------------------------------------------------------------------------------------------------
# 2. Data: one LC-MS run order per data set. Parameter table ("chosen" = a choice of this simulator)
# ---------------------------------------------------------------------------------------------------------------
# Feature j, injection i, on the log2 scale:  y = level_j + effect_j[group_i] + bio_j[extract_i] + tech_ji
#                                                 + b_j g(t_i) + s_j h_i + a_j f_i,   then values under the limit -> NA.
# | Parameter | Value | Source |
# |---|---|---|
# | Features | 300: two assays (pos, neg) of 150 | task |
# | Internal standards | 10: 5 per assay, beside the 300 | task (10); 5 per assay is a reading: the package treats each assay apart |
# | Injection order | position 1 a pooled QC, then a QC after every 5 study injections (every 6th position) and after the last one | task; S19 "one QC sample at every fifth injection", S20 "at least every 10th"; the first and last QC are chosen so that the QC curve is never extrapolated |
# | Order, f, batch | shared by the two assays (one extract, one sequence) | chosen |
# | level_j | N(20, 2^2); standards 22 | chosen |
# | Residual variance | sigma_j^2 = 0.25 x, x scaled inverse chi-square with 10 df, rescaled to median 1 over the candidates: median residual SD 0.5 log2 | chosen (the ceiling arm's power at the centre is then neither 0 nor 1) |
# | Technical share | 0.25 of sigma_j^2 (technical SD 0.25 at the median, QC RSD about 17%); the rest is biological, per extract | chosen; under S19's "QC RSD over 20% removed" |
# | Standards' noise | technical only, variance 0.0625 | chosen (constant amount spiked) |
# | Effects | +/-1 log2 in group B on the cell's share of each assay (rounded), half up unless the cell says 80%; 3 groups: one group, drawn per feature, is shifted | ticket (1 log2); composition chosen |
# | Drift g(t) | exp(-2u), u the position scaled to [0, 1]; mean 0, variance 1 over the study injections; b_j normal | shape chosen: S19 and S20 describe drift along the run without numbers |
# | Sample-wise factor f | mean 0, variance 1 over study injections, correlation with the group score exactly the cell's (0 at the centre); a_j normal | ticket (normal loadings, correlation by construction) |
# | Batch step h | two batches along the run, standardised; s_j normal; first batch the largest multiple of the group count under half the run | task reading of "3 factors" |
# | Strength | variance of the unwanted part = 3 x median residual variance = 0.75, split equally over the components present (centre: drift 0.375 + factor 0.375; three: 0.25 each) | ticket (3); split chosen |
# | Correlation with group | f: exact. Drift: the order is the one of 2,000 random orders whose correlation of g with the group score is nearest the cell's (0 at the centre). Batch step in the 3-factor cells: nearest 0 | ticket ("fixed by construction"). Limit: orders are balanced, not simply randomised, so chance confounding of the drift is removed |
# | Pooled QC injection | level + the mean of effect and bio over the extracts + technical noise + drift and batch at its position + a_j times its own N(0, 1) score | chosen. The factor is injection-level: if QC injections did not carry it, no QC-based method could see it |
# | Standards | drift, step and factor with their own normal loadings of the same variances; no group effect | task |
# | Contaminated standards (arm 9) | the first standard of each assay (2 of 10) gets +1 log2 in group B | task |
# | Missing | per assay one detection limit, the quantile at the cell's share (25% at the centre) of all candidate study values; study and QC values under it are NA. The standards are spiked above the limit and never missing | task (detection limit); ticket (25%); standards chosen |
# | Presence rule | a feature is kept when it is observed at least once in every group. The 150 features of an assay are the first candidates of 2,000 that pass, by kind (nulls, up, down), so the effect share is exact. The share of missing values among the kept features is therefore below the cell's; it is written per run to s4-timing.tsv with the number of complete features | found here: without the rule the package's sweep returns "No valid percentage found" (the canonical-correlation plot is invalid) and its test fails ("prior.weights contain NA values"); the GUI's intensity filter stands before both. Holding the cell's share among the kept features instead is not possible in the cells without unwanted variation (S14 keeps 0.06% of candidates) |
# | Duplicates (C6) | 4 extracts per group, 2 injections each: bio shared, technical noise, drift and f per injection | task |
# | Known batch (S8) | drift + batch step (0.375 each), no f; batch 1 has 4 of group A and 2 of B, batch 2 the reverse | task reading; shares chosen |
# | Heavy tails (S7) | bio and technical noise t with 3 df scaled to variance 1 | chosen |
# | Loadings tied to abundance (S10) | every loading = SD x (0.7 z + sqrt(0.51) e), z the standardised level | chosen |
# | Rows | permuted within each assay; ids carry no information (R-19) | ticket |
# Cells (index for the seed): C0 0 centre (6 per group, 20% effects, drift + f, correlation 0, 25% missing, 2 groups);
# C1 1 12 per group; C2 2 10% effects; C3 3 no unwanted variation; C4 4 no effects; C5 5 correlation 0.3; C6 6
# duplicates; C7 7 3 groups; S1 8 3 per group; S2 9 40% effects; S3 10 correlation 0.6; S4 11 3 factors; S5 12 0%
# missing; S6 13 50% missing; S7 14 heavy tails; S8 15 known unbalanced batch; S9 16 80% up; S10 17 loadings tied to
# abundance; S11 18 3 per group, correlation 0.3; S12 19 3 per group, 3 factors; S13 20 no unwanted variation, 3 per
# group; S14 21 no unwanted variation, 50% missing. C8 22 (the exactness round, user 2026-11-18, decision M16):
# C6 with no unwanted variation (duplicates; data generation as C6 with no drift and no f), run when the cells
# argument names it; reported, judges no design.
# Seed = 100000 + 10000 x cell index + run: the data of a (cell, run), made once and shared by every arm. Each arm
# starts its own stream at set.seed(seed * 100 + arm number) (I = 18, Qs = 19, Qm = 20, L = 21); E uses seed * 100 + 99.
#
# ---------------------------------------------------------------------------------------------------------------
# 3. Arms. Pooled QC injections are never tested. Arms 1 to 17 and I drop them before anything; Q and L use them.
# ---------------------------------------------------------------------------------------------------------------
# "W in the model": W from ruv::RUVIII(Y, M, ctl, k, return.info = TRUE) on the features (and standards) complete in
# every injection, centred, added to the formula as W1..Wk; the test is on the uncorrected values, missing ones kept.
# "Subtract": RUVIIIC::RUVIII_C_Varying(k, Y, M, toCorrect = features, potentialControls) called directly, missing
# values kept, then the plain test. The GUI offers neither all features nor the standards as controls, nor pooled QC
# or duplicates as replicates: these arms are computed here, not by a package route.
#  1  no correction.
#  2  ceiling: the true drift (and batch step) subtracted, the true f in the model (f absent in S8; arm 1 where the
#     cell has no unwanted variation).
#  3  today's design: the sweep, ruvIII_C_Varying (group as replicates), runLipidsDA on the corrected object.
#  4  today's controls and k, W in the model (M = group, as in every arm to 11).
#  5  controls = features outside the 25% most significant of arm 1's test (smallest p over the contrasts), k by
#     today's rule on those controls (ruvCancor, then findBestKElbow), W in the model.
#  6  controls = features with first-pass BH q above 0.5 in every contrast, k as arm 5, W in the model.
#  7  all features as controls, k = E, W in the model.          11  arm 7 with k = 1.
#  8  the standards as controls, k = E, W in the model.          9  arm 8 with the contaminated standards.
# 10  arm 8, subtract.
# 12  ruv::RUV4(Y, X, ctl = standards, k = E, Z = 1), ruv::variance_adjust()$p.ebayes, q by the package's rule.
# 13  ruv::RUVinv(Y, X, ctl = standards, Z = 1), p.ebayes likewise; no k.
# 14  sva::sva(dat, mod = ~ group, mod0 = ~ 1, n.sv = E), the surrogate variables in the model.
# 15  (C6) duplicates as replicates, all features as controls, k = E, subtract.     16  (C6) the same, W in the model.
# 17  (S8) the known batch in the model: "~ 0 + group + batch".
# 22  on request (K124_ARMS): arm 8 (the standards as controls, W in the model) with k by today's rule on the
#     standards (ruvCancor, then findBestKElbow, as arm 5 on its controls). The standards are not rows of the sweep's
#     object, so the rule reads an object of the features and the standards, prepared as the sweep's.
# 23  on request (K124_ARMS): arm 5's controls with k = E, W in the model.
#     Arms 22 and 23 (user, 2026-11-18; 17 to 21 are taken as arm or stream numbers in one simulator or another)
#     cross the two k rules over the two control sets: 22 differs from 8, and 23 from 5, in k only.
#     There is no arm 24 here (S1 to S3: arm 5 with today's rule on arm 5's own controls): arm 5 already is that
#     design in this simulator. A 24 in K124_ARMS is ignored.
#     The exactness round (user, 2026-11-18, decision M16; kit-124/exact-round-declaration.md), on request
#     (K124_ARMS). Each arm is another arm with the named change and nothing else, on that arm's random stream.
# 25  arm 7 with each comparison's standard error exact per feature. limma::contrasts.fit forms it from one
#     correlation matrix for all features, which is approximate for a feature with missing values when the model
#     carries W. Here, per comparison, the model matrix is reparameterised so that the comparison is its first
#     coefficient (limma::contrastAsCoef), its columns join the design table, and the same kernel is called with
#     the formula of those columns and the contrast "comparison": for a contrast that is one coefficient
#     contrasts.fit is exact. Estimates, residual variances, their moderation (trend, robust) and q are as before.
# 29  (C6, C8) arm 7 with E estimated after the extracts (the biological units) are removed, not only the groups:
#     sva::num.sv(dat, mod = ~ extract, method = "be", B = 20, E's seed).       30  (C6, C8) arm 16 with arm 29's E.
# 31  arm 7 with E from 200 permutations: sva::num.sv(dat, mod = ~ group, method = "be", B = 200, E's seed).
# 32  arm 7 with all of them: arm 25's standard error, B = 200, and in C6 and C8 the extracts removed.
# 33  (C6, C8) arm 16 with all of them.
#     Arm numbers 26 to 28 (S1: the test applies its trend and robust settings) are not formed: this kernel has no
#     TREAT step (a cut-off of 0 is passed on as NA) and reports the p-values of eBayes(trend, robust).
#     With E = 0 these arms are the no-RUV fit (for 25, 32, 33 with the reparameterised matrix: equal to arm 1's up
#     to rounding, because without W the model's coefficients are uncorrelated); a refused assay is arm 1's, as before.
#     K124_EXACT_CHECK=1 (environment; a check, not a result): per data set and view, the reparameterised no-RUV fit
#     must give arm 1's estimates and p-values to 1e-8; otherwise arm 25 crashes with the reason.
#     K124_ARMS (environment; unset = the arms of this list but 22 to 33): a comma list of the arm labels to run,
#     of those the cell has and 22, 23, 25, 29 to 33. Arm 1 must be in it. The data, E and every stream are as without it, so an
#     arm's rows are the same whichever other arms run.
#     The confirmation round (KIT-124 / decision M29; kit-124/scored/confirm-round-declaration.md), on request
#     (K124_ARMS), with K124_PACKAGE=installed (the MultiScholaR installed in the job from the bound commit). The RUV
#     step and the test are the package's new-design kernels, as kit-129/final-method/prove-arms.R calls them, on an
#     object of the study injections (the package holds no pooled-QC injection; S4's own W is formed on the features
#     complete in study and QC injections, so these arms' controls can differ from arm 32's): runOmicsRuvFactors(object,
#     <design>, ruv_k_rule = "num_sv_be", ruv_control_ids, ruv_grouping_variable, ruv_k_seed = E's seed), then
#     runLipidsDA(object, contrasts, "~ 0 + group", 0.05, 0, TRUE, TRUE, ruv_factors = its factors). The same is run
#     through the metabolomics classes (MetaboliteAssayData, runMetabolitesDA) and the arm crashes unless the two tables
#     agree (identical, or within 1e-12); the lipid table is scored. An assay the kernel refuses (no RUV, its warning)
#     is a refused assay, scored as arm 1 (R-07).
# 34  recommended: all features as controls, the groups as replicate sets, the estimated count (simulator arm 32).
# 35  supplied controls: the standards (arm 8's), the estimated count. 36  the contaminated standards (arm 9's).
# 37  no RUV: runLipidsDA on the object, no factors (simulator arm 1).
# 38  (C6, C8) technical replicates: the extract as the replicate sets (ruv_grouping_variable = "extract"), all
#     features as controls, the estimated count (simulator arm 33).
#     Not formed: 39 and 40 (first-pass nulls): the kernel has no first-pass control rule.
#     K124_FAULT=<arm> (environment, the lifecycle test only): that arm stops with an unclassed error in the first run
#     of the invocation, so that a crash of a shipped-kernel arm can be seen to be reported.
#  I  the package's ITSD normalisation and log2 alone, then the test.
#  Qs pooled QC injections as the one replicate set (M over study and QC injections), all features as controls,
#     k = E, subtract; the study injections are tested.            Qm  the same, W in the model.
#  L  QC-RLSC: per feature stats::loess of the QC values on position, subtracted from the study values (division on
#     the raw scale) and the QC median added back. degree 2 and span 0.75 (the stats::loess defaults; chosen: S19 as
#     read in convention.md gives neither) with 10 QC injections or more; with fewer (every cell here: 3 to 6),
#     degree 1 and span 1.5, a weighted straight line through all QC, because a local quadratic through 3 to 6
#     points reproduces QC noise (chosen). A feature with fewer than 3 QC values is left as it is.
# E = sva::num.sv(dat, mod = ~ group, method = "be", B = 20, seed) per assay on the features complete in every
#     injection (study samples only in dat); it may be 0, and the arm is then arm 1. Reported k: mean over assays.
# Subsets: `full` = all 300 features; `complete` = the arm run again on the features complete in every injection
#     (study and QC). Arms 12 to 14 need complete data and have the `complete` row only. Where a cell has no missing
#     values the `complete` rows repeat the `full` ones.
# Status: `refused` = the package refuses (the sweep finds no valid percentage; k above the complete controls) or
#     the arm cannot be formed (k above the controls, above the rank the replicates give, or no residual degrees of
#     freedom; RUVinv with fewer controls than residual dimensions); the refused assay is scored as arm 1 (R-07).
#     `crashed` = any other error; the row then has no measures. Reasons: s4-notes.tsv.
# Measures over the (feature, contrast) tests of both assays: n_features = tests with a p-value; n_true = true
#     effects in the subset (tested or not); n_calls; null_p05, null_p01 = share of tested nulls under 0.05, 0.01;
#     fdp (0 with no calls); power = true calls / n_true; bias = mean of sign(truth) x (estimate - truth) and rmse
#     over tested true effects. `seconds` includes the shared steps the arm uses (E, the sweep, the first pass).
# Technical replicates (KIT-129 R0, on the `kit-128` tree; this file needs that tree): the design declares the extract
#     as the technical replicate of each injection, as the design step declares `tech_rep_group`. The objects carry it
#     (`technical_replicate_id = "extract"`, read by runLipidsDA and runMetabolitesDA) and the kernel is called with
#     `group_id = "group", technical_replicate_id = "extract"`, in every cell. In C6 (4 extracts per group, 2 injections
#     each) every arm's test, the ceiling included, is therefore fitted with the extract as block
#     (limma::duplicateCorrelation, then lmFit(block, correlation); KIT-128). In every other cell each extract is one
#     injection, the package's rule finds no duplicates and the fit is the unblocked one, as before. Before KIT-128
#     the 16 injections of C6 were tested as 16 independent samples (the KIT-124 confirm rows of C6: void).
#     The RUV step does not read the declaration (the group, or the arm's own replicate matrix, as before).

if (!identical(Sys.getenv("OMP_NUM_THREADS"), "1")) stop("run with OMP_NUM_THREADS=1 (RUVIIIC is threaded)")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L) stop("usage: Rscript s4_lipid_metab.R <out dir> <runs per core cell> <runs per stress cell> <cores> [cells = all] [first run = 1]")
out_dir <- args[[1L]]
runs_core <- as.integer(args[[2L]])
runs_stress <- as.integer(args[[3L]])
cores <- as.integer(args[[4L]])
cells_arg <- if (length(args) >= 5L) args[[5L]] else "all"
first_run <- if (length(args) >= 6L) as.integer(args[[6L]]) else 1L
arms_wanted <- Sys.getenv("K124_ARMS", "")
arms_wanted <- if (nzchar(arms_wanted)) strsplit(arms_wanted, ",", fixed = TRUE)[[1L]]
stopifnot(is.null(arms_wanted) || "1" %in% arms_wanted)

installed <- identical(Sys.getenv("K124_PACKAGE"), "installed")
if (installed) {
    # the confirmation round: the installed package, every object of its namespace on the search path as load_all puts it
    suppressPackageStartupMessages(library(MultiScholaR))
    attach(asNamespace("MultiScholaR"), name = "MultiScholaR:namespace", warn.conflicts = FALSE)
} else suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
SHIPPED_ARMS <- c("34", "35", "36", "37", "38")
if (!installed && any(SHIPPED_ARMS %in% arms_wanted)) stop("arms 34 to 38 call the installed package: set K124_PACKAGE=installed")
fault_arm <- Sys.getenv("K124_FAULT", "")
logger::log_threshold(logger::FATAL)
logger::log_threshold(logger::FATAL, namespace = "MultiScholaR")
# loaded before the workers fork, so that no run pays for it
for (package in c("limma", "qvalue", "ruv", "RUVIIIC", "sva", "statmod")) loadNamespace(package)

# ---- constants -------------------------------------------------------------------------------------------------
P_ASSAY <- 150L
POOL <- 2000L
N_IS <- 5L
ASSAYS <- c("pos", "neg")
MED_RES_VAR <- 0.25
TECH_SHARE <- 0.25
STRENGTH <- 3
ID_COL <- "database_identifier"
RUV_PARAMS <- list(ruv_grouping_variable = "group", percentage_min = 1, percentage_max = 20,
    separation_metric = "max_difference", k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, manual_k = 2,
    manual_percentage = 5)
ARM_NUMBER <- c(stats::setNames(1:17, 1:17), I = 18L, Qs = 19L, Qm = 20L, L = 21L, "22" = 22L, "23" = 23L,
    "25" = 25L, "29" = 29L, "30" = 30L, "31" = 31L, "32" = 32L, "33" = 33L,
    "34" = 34L, "35" = 35L, "36" = 36L, "37" = 37L, "38" = 38L)
EXACT_STREAM <- c("25" = 7L, "29" = 7L, "30" = 16L, "31" = 7L, "32" = 7L, "33" = 16L)   # the stream of the arm each varies
COMPLETE_ONLY <- c("12", "13", "14")

cellTable <- function() {
    cell <- function(id, class, n = 6L, groups = 2L, effects = 0.2, rho = 0, unwanted = "two", missing = 0.25,
                     tails = FALSE, up = 0.5, abundance = FALSE, duplicates = FALSE) {
        list(id = id, class = class, n = n, groups = groups, effects = effects, rho = rho, unwanted = unwanted,
            missing = missing, tails = tails, up = up, abundance = abundance, duplicates = duplicates)
    }
    cells <- list(cell("C0", "core"), cell("C1", "core", n = 12L), cell("C2", "core", effects = 0.1),
        cell("C3", "core", unwanted = "none"), cell("C4", "core", effects = 0), cell("C5", "core", rho = 0.3),
        cell("C6", "core", n = 4L, duplicates = TRUE), cell("C7", "core", groups = 3L),
        cell("S1", "stress", n = 3L), cell("S2", "stress", effects = 0.4), cell("S3", "stress", rho = 0.6),
        cell("S4", "stress", unwanted = "three"), cell("S5", "stress", missing = 0), cell("S6", "stress", missing = 0.5),
        cell("S7", "stress", tails = TRUE), cell("S8", "stress", unwanted = "batch"), cell("S9", "stress", up = 0.8),
        cell("S10", "stress", abundance = TRUE), cell("S11", "stress", n = 3L, rho = 0.3),
        cell("S12", "stress", n = 3L, unwanted = "three"), cell("S13", "stress", n = 3L, unwanted = "none"),
        cell("S14", "stress", unwanted = "none", missing = 0.5),
        cell("C8", "core", n = 4L, duplicates = TRUE, unwanted = "none"))
    for (i in seq_along(cells)) cells[[i]]$index <- i - 1L
    stats::setNames(cells, vapply(cells, `[[`, "", "id"))
}

armsOfCell <- function(cell) {
    arms <- c(as.character(1:14), "I", "Qs", "Qm", "L")
    if (cell$duplicates) arms <- c(arms, "15", "16")
    if (identical(cell$unwanted, "batch")) arms <- c(arms, "17")
    if (!is.null(arms_wanted)) arms <- intersect(c(arms, "22", "23", "25", "31", "32", if (cell$duplicates) c("29", "30", "33"),
        if (installed) c("34", "35", "36", "37", if (cell$duplicates) "38")), arms_wanted)
    arms[order(ARM_NUMBER[arms])]
}

# ---- small helpers ---------------------------------------------------------------------------------------------
quiet <- function(expr) {
    con <- file(nullfile(), "w")
    sink(con)
    sink(con, type = "message")
    on.exit({
        sink(type = "message")
        sink()
        close(con)
    })
    suppressWarnings(suppressMessages(expr))
}
popStandardise <- function(x, ref = x) (x - mean(ref)) / sqrt(mean((ref - mean(ref))^2))
noise <- function(n, tails) if (tails) stats::rt(n, 3) / sqrt(3) else stats::rnorm(n)
contrastsOf <- function(groups) {
    pairs <- utils::combn(groups, 2L)
    stats::setNames(lapply(seq_len(ncol(pairs)), function(i) c(pairs[1L, i], pairs[2L, i])),
        sprintf("group%s-group%s", pairs[2L, ], pairs[1L, ]))
}
pkgQ <- function(p) {
    # The q-value rule of applyLipidDaContrasts, for the arms whose p-values are not limma's (12, 13).
    q <- rep(NA_real_, length(p))
    ok <- which(is.finite(p))
    if (!length(ok)) return(q)
    v <- p[ok]
    q[ok] <- if (all(v == 0) || all(v == 1) || length(v) < 3L) stats::p.adjust(v, "BH") else
        tryCatch(qvalue::qvalue(v)$qvalues, error = function(e) stats::p.adjust(v, "BH"))
    q
}
refusal <- function(reason) structure(list(reason = reason), class = "s4_refusal")
isRefusal <- function(x) inherits(x, "s4_refusal")

# ---- data ------------------------------------------------------------------------------------------------------
generateData <- function(cell, seed) {
    set.seed(seed)
    groups <- LETTERS[seq_len(cell$groups)]
    units <- data.frame(group = rep(groups, each = cell$n), stringsAsFactors = FALSE)
    units$extract <- sprintf("E%02d", seq_len(nrow(units)))
    inj <- if (cell$duplicates) units[rep(seq_len(nrow(units)), each = 2L), ] else units
    n <- nrow(inj)
    # the sequence: QC first, a QC after every 5 study injections and after the last
    kind <- "Q"
    for (i in seq_len(n)) kind <- c(kind, "S", if (i %% 5L == 0L || i == n) "Q")
    pos_s <- which(kind == "S")
    pos_q <- which(kind == "Q")
    g_all <- exp(-2 * (seq_along(kind) - 1) / (length(kind) - 1))
    g_s <- popStandardise(g_all[pos_s])
    g_q <- popStandardise(g_all[pos_q], g_all[pos_s])
    n_first <- max(cell$groups, cell$groups * (n %/% (2L * cell$groups)))
    if (cell$duplicates || n_first >= n) n_first <- n %/% 2L
    batch_s <- ifelse(seq_len(n) <= n_first, 1, 2)
    batch_q <- ifelse(pos_q <= pos_s[n_first], 1, 2)
    h_s <- popStandardise(batch_s)
    h_q <- popStandardise(batch_q, batch_s)
    # the order: of 2,000 random orders, the one nearest the cell's correlations
    z <- popStandardise(match(inj$group, groups))
    candidates <- if (identical(cell$unwanted, "batch")) {
        replicate(2000L, {
            first <- c(sample(which(inj$group == "A"), 4L), sample(which(inj$group == "B"), 2L))
            c(sample(first), sample(setdiff(seq_len(n), first)))
        })
    } else replicate(2000L, sample(n))
    zp <- matrix(z[candidates], nrow = n)
    miss <- abs(drop(crossprod(g_s, zp)) / n - cell$rho)
    if (identical(cell$unwanted, "three")) miss <- miss + abs(drop(crossprod(h_s, zp)) / n)
    inj <- inj[candidates[, which.min(miss)], ]
    z <- popStandardise(match(inj$group, groups))
    design <- data.frame(Sample_ID = sprintf("S%02d", seq_len(n)), group = inj$group, extract = inj$extract,
        batch = paste0("b", batch_s), pos = pos_s, stringsAsFactors = FALSE)
    qc_ids <- sprintf("QC%02d", seq_along(pos_q))
    # the sample-wise factor: exact correlation with the group score, orthogonal to the group otherwise
    e <- stats::residuals(stats::lm(stats::rnorm(n) ~ factor(design$group)))
    f_s <- cell$rho * z + sqrt(1 - cell$rho^2) * popStandardise(e)
    f_q <- stats::rnorm(length(pos_q))
    components <- switch(cell$unwanted, two = c("drift", "factor"), three = c("drift", "step", "factor"),
        batch = c("drift", "step"), none = character())
    v_each <- if (length(components)) STRENGTH * MED_RES_VAR / length(components) else 0
    extract_of <- match(design$extract, units$extract)
    contrasts <- contrastsOf(groups)

    assays <- lapply(stats::setNames(ASSAYS, ASSAYS), function(assay) {
        # candidate features; the assay keeps the first of each kind (null, up, down) that pass the presence rule
        p <- POOL
        level <- stats::rnorm(p, 20, 2)
        x <- 10 / stats::rchisq(p, 10)
        sigma2 <- MED_RES_VAR * x / stats::median(x)
        z_level <- c((level - 20) / 2, rep(1, N_IS))
        loading <- function(present) {
            if (!present) return(rep(0, p + N_IS))
            e <- stats::rnorm(p + N_IS)
            sqrt(v_each) * if (cell$abundance) 0.7 * z_level + sqrt(0.51) * e else e
        }
        b <- loading("drift" %in% components)
        s <- loading("step" %in% components)
        a <- loading("factor" %in% components)
        n_eff <- round(cell$effects * P_ASSAY)
        n_up <- round(cell$up * n_eff)
        wanted <- c("0" = P_ASSAY - n_eff, "1" = n_up, "-1" = n_eff - n_up)
        sign <- sample(c(0, 1, -1), p, replace = TRUE, prob = pmax(wanted, 1e-12))
        target <- if (length(groups) == 2L) rep(2L, p) else sample(length(groups), p, replace = TRUE)
        delta <- matrix(0, p, length(groups), dimnames = list(NULL, groups))
        delta[cbind(seq_len(p), target)] <- sign
        bio <- matrix(noise(p * nrow(units), cell$tails), p) * sqrt((1 - TECH_SHARE) * sigma2)
        unit_part <- delta[, units$group, drop = FALSE] + bio
        unwanted <- function(rows, g, h, f) outer(b[rows], g) + outer(s[rows], h) + outer(a[rows], f)
        feat <- seq_len(p)
        std <- p + seq_len(N_IS)
        y_s <- level + unit_part[, extract_of, drop = FALSE] + matrix(noise(p * n, cell$tails), p) * sqrt(TECH_SHARE * sigma2) +
            unwanted(feat, g_s, h_s, f_s)
        y_q <- level + rowMeans(unit_part) + matrix(noise(p * length(pos_q), cell$tails), p) * sqrt(TECH_SHARE * sigma2) +
            unwanted(feat, g_q, h_q, f_q)
        is_s <- 22 + matrix(noise(N_IS * n, cell$tails), N_IS) * 0.25 + unwanted(std, g_s, h_s, f_s)
        is_q <- 22 + matrix(noise(N_IS * length(pos_q), cell$tails), N_IS) * 0.25 + unwanted(std, g_q, h_q, f_q)
        # the detection limit: the cell's share of all candidate study values lies under it; then the presence rule
        limit <- if (cell$missing > 0) stats::quantile(y_s, cell$missing, names = FALSE) else -Inf
        seen <- y_s >= limit
        kept <- Reduce(`&`, lapply(groups, function(g) rowSums(seen[, design$group == g, drop = FALSE]) >= 1L))
        rows <- unlist(lapply(names(wanted), function(kind) {
            have <- which(kept & sign == as.numeric(kind))
            if (length(have) < wanted[[kind]]) stop("the candidate pool is too small for this cell")
            have[seq_len(wanted[[kind]])]
        }))
        rows <- rows[sample(length(rows))]
        censor <- function(m) {
            m[m < limit] <- NA_real_
            m
        }
        isc_s <- is_s
        isc_s[1L, design$group == "B"] <- isc_s[1L, design$group == "B"] + 1
        isc_q <- is_q
        isc_q[1L, ] <- isc_q[1L, ] + mean(units$group == "B")
        truth <- vapply(contrasts, function(pair) delta[, pair[[2L]]] - delta[, pair[[1L]]], numeric(p))
        known <- outer(b[feat], g_s) + outer(s[feat], h_s)
        ids <- sprintf("%s_f%03d", assay, seq_along(rows))
        is_ids <- sprintf("%s_ITSD_%d", assay, seq_len(N_IS))
        named <- function(m, row_ids, cols) {
            dimnames(m) <- list(row_ids, cols)
            m
        }
        list(y = named(censor(y_s)[rows, , drop = FALSE], ids, design$Sample_ID),
            qc = named(censor(y_q)[rows, , drop = FALSE], ids, qc_ids),
            is = named(is_s, is_ids, design$Sample_ID), is_qc = named(is_q, is_ids, qc_ids),
            isc = named(isc_s, is_ids, design$Sample_ID), isc_qc = named(isc_q, is_ids, qc_ids),
            truth = named(truth[rows, , drop = FALSE], ids, names(contrasts)),
            known = named(known[rows, , drop = FALSE], ids, design$Sample_ID))
    })
    list(cell = cell, seed = seed, design = design, f = f_s, has_factor = "factor" %in% components,
        pos_q = pos_q, qc_ids = qc_ids, contrasts = names(contrasts), contrast_pairs = contrasts, assays = assays)
}

completeFeatures <- function(assay) rowSums(!is.finite(cbind(assay$y, assay$qc))) == 0L

restrictData <- function(data, keep) {
    for (name in names(data$assays)) {
        rows <- keep[[name]]
        for (slot in c("y", "qc", "truth", "known")) data$assays[[name]][[slot]] <- data$assays[[name]][[slot]][rows, , drop = FALSE]
    }
    data
}

# ---- package calls ---------------------------------------------------------------------------------------------
assayTable <- function(m) {
    table <- tibble::as_tibble(as.data.frame(m, check.names = FALSE))
    table <- tibble::add_column(table, !!ID_COL := rownames(m), .before = 1L)
    tibble::add_column(table, name = rownames(m), .after = 1L)
}
lipidObject <- function(mats, design) {
    methods::new("LipidomicsAssayData", lipid_data = lapply(mats, assayTable), lipid_id_column = ID_COL,
        annotation_id_column = "name", design_matrix = design[, c("Sample_ID", "group", "extract")], sample_id = "Sample_ID",
        group_id = "group", technical_replicate_id = "extract", args = list())
}
metabObject <- function(mats, design) {
    methods::new("MetaboliteAssayData", metabolite_data = lapply(mats, assayTable), metabolite_id_column = ID_COL,
        annotation_id_column = "name", design_matrix = design[, c("Sample_ID", "group", "extract")], sample_id = "Sample_ID",
        group_id = "group", technical_replicate_id = "extract", args = list())
}
assayMatrix <- function(table, samples) {
    m <- as.matrix(as.data.frame(table[, samples], check.names = FALSE))
    rownames(m) <- table[[ID_COL]]
    m
}
longTable <- function(results) {
    do.call(rbind, lapply(names(results), function(contrast) {
        t <- results[[contrast]]
        data.frame(feature = rownames(t), contrast = contrast, est = t$logFC, p = t$P.Value, q = t$fdr_qvalue,
            bh = t$fdr_value_bh, stringsAsFactors = FALSE)
    }))
}
# The kernel of the default path: limma, trend and robust, no TREAT.
daKernel <- function(m, design, formula, contrasts, exact = FALSE) {
    if (anyNA(m)) {
        # A feature whose observed injections cannot estimate every coefficient makes the kernel fail for the whole
        # assay ("prior.weights contain NA values"): it is left out and counts as not tested.
        x <- stats::model.matrix(stats::as.formula(formula), design)
        seen <- is.finite(m)
        pattern <- drop(seen %*% 2^(seq_len(ncol(m)) - 1))
        estimable <- vapply(split(seq_len(nrow(m)), pattern), function(i) qr(x[seen[i[[1L]], ], , drop = FALSE])$rank == ncol(x), logical(1L))
        m <- m[estimable[as.character(pattern)], , drop = FALSE]
    }
    if (exact) {
        # the exactness round: per comparison, the model reparameterised so that the comparison is its first coefficient
        x <- stats::model.matrix(stats::as.formula(formula), design)
        return(do.call(rbind, lapply(contrasts, function(contrast) {
            q <- limma::contrastAsCoef(x, limma::makeContrasts(contrasts = contrast, levels = colnames(x)))$design
            colnames(q) <- c("comparison", paste0("Q", seq_len(ncol(q))[-1L]))
            table <- daKernel(m, cbind(design[setdiff(names(design), colnames(q))], q),
                paste("~ 0 +", paste(colnames(q), collapse = " + ")), "comparison")
            table$contrast <- contrast
            table
        })))
    }
    longTable(quiet(runTestsContrastsLipidDA(data = m, contrast_strings = contrasts, design_matrix = design,
        formula_string = formula, sample_id_col = "Sample_ID", treat_lfc_cutoff = NA, eBayes_trend = TRUE,
        eBayes_robust = TRUE, group_id = "group", technical_replicate_id = "extract"))$results)
}
contrastsTable <- function(contrasts) data.frame(contrasts = contrasts, friendly_names = contrasts, stringsAsFactors = FALSE)

# Lipid and metab share every kernel S4 uses, and the default path is the one the header states: checked on one
# data set, identical outputs or stop.
checkSharedKernels <- function() {
    data <- generateData(cellTable()$C0, 1L)
    mats <- lapply(data$assays, `[[`, "y")
    lipid <- lipidObject(mats, data$design)
    metab <- metabObject(mats, data$design)
    checks <- list()
    sweep_l <- quiet(runLipidPerAssayRuvOptimization(lipid, "automatic", RUV_PARAMS))
    sweep_m <- quiet(runPerAssayRuvOptimization(metab, "automatic", RUV_PARAMS))
    part <- function(x) lapply(x, `[`, c("success", "best_k", "best_percentage", "control_genes_index", "optimization_results"))
    checks$sweep <- identical(part(sweep_l), part(sweep_m))
    ruv_l <- quiet(ruvIII_C_Varying(lipid, ruv_grouping_variable = "group", ruv_number_k = extractLipidBestKPerAssay(sweep_l),
        ctrl = extractLipidCtrlPerAssay(sweep_l)))
    ruv_m <- quiet(ruvIII_C_Varying(metab, ruv_grouping_variable = "group", ruv_number_k = extractBestKPerAssay(sweep_m),
        ctrl = extractCtrlPerAssay(sweep_m)))
    checks$ruv <- identical(ruv_l@lipid_data, ruv_m@metabolite_data)
    tbl <- contrastsTable(data$contrasts)
    da_l <- quiet(runLipidsDA(ruv_l, tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE))$da_lipids_long
    da_m <- quiet(runMetabolitesDA(ruv_m, tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE))$da_metabolites_long
    cols <- c("logFC", "P.Value", "fdr_qvalue", "fdr_value_bh", "assay", "comparison", "significant")
    checks$orchestrator <- identical(da_l$lipid_id, da_m$metabolite_id) && identical(da_l[cols], da_m[cols])
    samples <- data$design$Sample_ID
    kernel_equal <- direct_equal <- kernels_equal <- TRUE
    for (assay in ASSAYS) {
        m <- assayMatrix(ruv_l@lipid_data[[assay]], samples)
        k_l <- daKernel(m, data$design, "~ 0 + group", data$contrasts)
        k_m <- longTable(quiet(runTestsContrastsMetabDA(m, data$contrasts, data$design, "~ 0 + group", "Sample_ID", NA, TRUE, TRUE,
            group_id = "group", technical_replicate_id = "extract"))$results)
        kernels_equal <- kernels_equal && identical(k_l, k_m)
        o <- da_l[da_l$assay == assay, ]
        at <- match(o$lipid_id, k_l$feature)
        kernel_equal <- kernel_equal && identical(k_l$p[at], o$P.Value) && identical(k_l$q[at], o$fdr_qvalue) &&
            identical(which(k_l$q[at] < 0.05), which(o$significant != "NS"))
        dm <- stats::model.matrix(~ 0 + group, data$design)
        fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(m, dm), limma::makeContrasts(contrasts = data$contrasts, levels = colnames(dm))),
            trend = TRUE, robust = TRUE)
        direct_equal <- direct_equal && isTRUE(all.equal(unname(fit$p.value[k_l$feature, 1L]), k_l$p, tolerance = 1e-12))
    }
    checks$kernel_lipid_equals_metab <- kernels_equal
    checks$orchestrator_equals_kernel <- kernel_equal
    checks$kernel_equals_limma_trend_robust <- direct_equal
    raw <- lapply(data$assays, function(a) 2^rbind(a$y, a$is))
    ids <- lapply(data$assays, function(a) rownames(a$is))
    itsd_l <- quiet(logTransformAssays(normaliseUntransformedData(lipidObject(raw, data$design), method = "ITSD",
        itsd_aggregation = "median", itsd_feature_ids = ids), offset = 1))
    itsd_m <- quiet(logTransformAssays(normaliseUntransformedData(metabObject(raw, data$design), method = "ITSD",
        itsd_aggregation = "median", itsd_feature_ids = ids), offset = 1))
    checks$itsd <- identical(itsd_l@lipid_data, itsd_m@metabolite_data) && nrow(itsd_l@lipid_data$pos) == P_ASSAY
    if (!all(unlist(checks))) stop("lipid/metab sharing check failed: ", paste(names(checks)[!unlist(checks)], collapse = ", "))
    checks
}

# ---- RUV building blocks ---------------------------------------------------------------------------------------
replicateMatrix <- function(blocks) {
    m <- ruv::replicate.matrix(blocks)
    rownames(m) <- names(blocks)
    m
}
groupM <- function(design) replicateMatrix(stats::setNames(design$group, design$Sample_ID))
qcM <- function(design, qc_ids) replicateMatrix(stats::setNames(c(design$Sample_ID, rep("QC", length(qc_ids))), c(design$Sample_ID, qc_ids)))

# W of RUV-III on complete columns. `y`: injections by columns, complete. Returns W for `keep` rows, centred.
ruvW <- function(y, m, ctl, k, keep = rownames(y), residual_df) {
    if (is.na(k)) return(refusal("k could not be estimated"))
    if (k == 0L) return(NULL)
    if (k > sum(ctl)) return(refusal(sprintf("k = %d above %d complete controls", k, sum(ctl))))
    if (k > nrow(y) - ncol(m)) return(refusal(sprintf("k = %d above the rank %d the replicates give", k, nrow(y) - ncol(m))))
    if (residual_df - k < 1L) return(refusal(sprintf("no residual degrees of freedom with k = %d", k)))
    # W as computeRuvCovariates forms it (KIT-041, KIT-123 R2): the controls centred, times ac' (ac ac')^-1
    ac <- ruv::RUVIII(Y = y, M = m[rownames(y), , drop = FALSE], ctl = ctl, k = k, return.info = TRUE)$fullalpha[seq_len(k), ctl, drop = FALSE]
    w <- sweep(y[, ctl, drop = FALSE], 2L, colMeans(y[, ctl, drop = FALSE])) %*% t(ac) %*% solve(ac %*% t(ac))
    rownames(w) <- rownames(y)
    w <- scale(w[keep, , drop = FALSE], center = TRUE, scale = FALSE)
    colnames(w) <- paste0("W", seq_len(k))
    w
}
# RUV-III-C called directly, missing values kept. `y`: injections by columns.
ruvSubtract <- function(y, m, controls, to_correct, k, keep = rownames(y)) {
    if (is.na(k)) return(refusal("k could not be estimated"))
    if (k == 0L) return(t(y[keep, to_correct, drop = FALSE]))
    complete <- sum(colSums(!is.finite(y[, controls, drop = FALSE])) == 0L)
    if (k > complete) return(refusal(sprintf("k = %d above %d controls observed in every sample", k, complete)))
    if (k > nrow(y) - ncol(m)) return(refusal(sprintf("k = %d above the rank %d the replicates give", k, nrow(y) - ncol(m))))
    out <- quiet(RUVIIIC::RUVIII_C_Varying(k = k, Y = y, M = m[rownames(y), , drop = FALSE], toCorrect = to_correct,
        potentialControls = controls, progress = FALSE))
    t(out[keep, to_correct, drop = FALSE])
}
testWithW <- function(view, assay, w, y = view$assays[[assay]]$y, exact = FALSE) {
    if (isRefusal(w)) return(w)
    if (is.null(w)) return(daKernel(y, view$design, "~ 0 + group", view$contrasts, exact))
    design <- cbind(view$design, as.data.frame(w[view$design$Sample_ID, , drop = FALSE]))
    daKernel(y, design, paste("~ 0 + group +", paste(colnames(w), collapse = " + ")), view$contrasts, exact)
}
residualDf <- function(view) nrow(view$design) - length(unique(view$design$group))
completeY <- function(a, extra = NULL, with_qc = FALSE) {
    # injections by complete columns; `extra` = standards rows to append when complete
    y <- if (with_qc) cbind(a$y, a$qc) else a$y
    if (!is.null(extra)) y <- rbind(y, extra[rowSums(!is.finite(extra)) == 0L, colnames(y), drop = FALSE])
    t(y[c(a$complete, rep(TRUE, nrow(y) - nrow(a$y))), , drop = FALSE])
}
controlIds <- function(index, features) {
    # a sweep's control set as feature ids: by name when the vector is named, else by position
    if (!is.null(names(index))) return(names(index)[which(index)])
    if (length(index) != length(features)) stop("the sweep's control vector cannot be aligned with the assay")
    features[which(index)]
}
elbowK <- function(view, controls, object = view$cancor_object()) {
    # today's rule for k on a given control set: ruvCancor on the sweep's prepared object, then the elbow
    plots <- quiet(ruvCancor(object, ctrl = controls, ruv_grouping_variable = "group"))
    lapply(stats::setNames(ASSAYS, ASSAYS), function(assay) {
        plot <- if (is.list(plots) && assay %in% names(plots)) plots[[assay]] else plots
        as.integer(findBestKElbow(plot))
    })
}

# ---- arms: each returns list(tables = per-assay long tables or refusals, k = per-assay k) ---------------------------
perAssay <- function(view, fn) {
    out <- lapply(stats::setNames(ASSAYS, ASSAYS), function(assay) fn(view$assays[[assay]], assay))
    list(tables = lapply(out, `[[`, "table"), k = vapply(out, function(x) as.numeric(x$k %||% NA_real_), numeric(1L)))
}
modelArm <- function(view, controls_of, k_of, m_of = function(view, a) groupM(view$design), standards = NULL, with_qc = FALSE,
                     exact = FALSE) {
    perAssay(view, function(a, assay) {
        k <- k_of(a, assay)
        extra <- if (!is.null(standards)) if (with_qc) cbind(a[[standards]], a[[paste0(standards, "_qc")]]) else a[[standards]]
        y <- completeY(a, extra, with_qc)
        ctl <- controls_of(a, assay, colnames(y))
        if (isRefusal(ctl)) return(list(table = ctl, k = 0))
        w <- ruvW(y, m_of(view, a), ctl, k, keep = view$design$Sample_ID, residual_df = residualDf(view))
        list(table = testWithW(view, assay, w, exact = exact), k = if (isRefusal(w)) 0 else k)
    })
}
subtractArm <- function(view, k_of, m, standards = NULL, with_qc = FALSE) {
    perAssay(view, function(a, assay) {
        k <- k_of(a, assay)
        y <- if (with_qc) cbind(a$y, a$qc) else a$y
        controls <- rownames(a$y)
        if (!is.null(standards)) {
            y <- rbind(y, a[[standards]])
            controls <- rownames(a[[standards]])
        }
        corrected <- ruvSubtract(t(y), m, controls, rownames(a$y), k, keep = view$design$Sample_ID)
        if (isRefusal(corrected)) return(list(table = corrected, k = 0))
        list(table = daKernel(corrected, view$design, "~ 0 + group", view$contrasts), k = k)
    })
}
allComplete <- function(a, assay, columns) rep(TRUE, length(columns))
standardsOnly <- function(a, assay, columns) grepl("_ITSD_", columns)
kE <- function(a, assay) a$E
# the exactness round: E's variants (runOne), the extracts as replicate sets, and the check of the reparameterisation
kVariant <- function(name) function(a, assay) a[[name]]
extractM <- function(view, a) replicateMatrix(stats::setNames(view$design$extract, view$design$Sample_ID))
hasDuplicates <- function(view) anyDuplicated(view$design$extract) > 0L
checkExact <- function(view) {
    if (!identical(Sys.getenv("K124_EXACT_CHECK"), "1")) return(invisible(NULL))
    for (assay in ASSAYS) {
        base <- view$base()[[assay]]
        exact <- daKernel(view$assays[[assay]]$y, view$design, "~ 0 + group", view$contrasts, exact = TRUE)
        at <- match(paste(base$feature, base$contrast), paste(exact$feature, exact$contrast))
        gap <- max(abs(exact$est[at] - base$est), abs(exact$p[at] - base$p))
        if (nrow(exact) != nrow(base) || !is.finite(gap) || gap > 1e-8) stop(sprintf("the reparameterised no-RUV fit is not arm 1's (%s): differs by up to %.3g", assay, gap))
    }
}

armFunctions <- list(
    "1" = function(view) list(tables = view$base(), k = c(NA_real_, NA_real_)),
    "2" = function(view) perAssay(view, function(a, assay) {
        y <- a$y - a$known
        if (!view$has_factor) return(list(table = daKernel(y, view$design, "~ 0 + group", view$contrasts), k = NA))
        design <- cbind(view$design, f = view$f)
        list(table = daKernel(y, design, "~ 0 + group + f", view$contrasts), k = NA)
    }),
    "3" = function(view) {
        sweep <- view$sweep()
        best_k <- extractLipidBestKPerAssay(sweep)
        object <- view$object()
        corrected <- tryCatch(quiet(ruvIII_C_Varying(object, ruv_grouping_variable = "group", ruv_number_k = best_k,
            ctrl = extractLipidCtrlPerAssay(sweep))), omics_ruv_insufficient_controls = function(e) refusal(conditionMessage(e)))
        refused <- lapply(stats::setNames(ASSAYS, ASSAYS), function(assay) {
            if (isRefusal(corrected)) return(corrected)
            if (!isTRUE(sweep[[assay]]$success)) return(refusal(sweep[[assay]]$error %||% "the sweep found no valid percentage"))
            if (!assay %in% names(corrected@lipid_data)) stop("ruvIII_C_Varying dropped assay ", assay)
            if (identical(corrected@lipid_data[[assay]], object@lipid_data[[assay]])) stop("ruvIII_C_Varying left assay ", assay, " uncorrected")
            NULL
        })
        tested <- object
        for (assay in ASSAYS) if (is.null(refused[[assay]])) tested@lipid_data[[assay]] <- corrected@lipid_data[[assay]]
        da <- quiet(runLipidsDA(tested, contrastsTable(view$contrasts), "~ 0 + group", 0.05, 0, TRUE, TRUE))$da_lipids_long
        tables <- lapply(stats::setNames(ASSAYS, ASSAYS), function(assay) {
            if (!is.null(refused[[assay]])) return(refused[[assay]])
            rows <- da[da$assay == assay, ]
            if (!nrow(rows)) stop("runLipidsDA returned no rows for assay ", assay)
            data.frame(feature = rows$lipid_id, contrast = rows$comparison, est = rows$logFC, p = rows$P.Value,
                q = rows$fdr_qvalue, bh = rows$fdr_value_bh, stringsAsFactors = FALSE)
        })
        list(tables = tables, k = vapply(ASSAYS, function(assay) if (is.null(refused[[assay]])) as.numeric(best_k[[assay]]) else 0, numeric(1L)))
    },
    "4" = function(view) {
        sweep <- view$sweep()
        modelArm(view, function(a, assay, columns) {
            if (!isTRUE(sweep[[assay]]$success)) return(refusal(sweep[[assay]]$error %||% "the sweep found no valid percentage"))
            columns %in% controlIds(sweep[[assay]]$control_genes_index, rownames(a$y))
        }, function(a, assay) if (isTRUE(sweep[[assay]]$success)) as.integer(sweep[[assay]]$best_k) else 0L)
    },
    "5" = function(view) {
        controls <- lapply(view$first_pass(), function(fp) fp$feature[rank(fp$p, na.last = "keep", ties.method = "first") > 0.25 * sum(is.finite(fp$p)) & is.finite(fp$p)])
        empiricalControlsArm(view, controls)
    },
    "6" = function(view) empiricalControlsArm(view, lapply(view$first_pass(), function(fp) fp$feature[is.finite(fp$bh) & fp$bh > 0.5])),
    "7" = function(view) modelArm(view, allComplete, kE),
    "8" = function(view) modelArm(view, standardsOnly, kE, standards = "is"),
    "9" = function(view) modelArm(view, standardsOnly, kE, standards = "isc"),
    "10" = function(view) subtractArm(view, kE, groupM(view$design), standards = "is"),
    "11" = function(view) modelArm(view, allComplete, function(a, assay) 1L),
    "12" = function(view) referenceArm(view, "RUV4"),
    "13" = function(view) referenceArm(view, "RUVinv"),
    "14" = function(view) perAssay(view, function(a, assay) {
        if (is.na(a$E)) return(list(table = refusal("k could not be estimated"), k = 0))
        if (a$E == 0L) return(list(table = testWithW(view, assay, NULL), k = 0))
        if (residualDf(view) - a$E < 1L) return(list(table = refusal(sprintf("no residual degrees of freedom with k = %d", a$E)), k = 0))
        mod <- stats::model.matrix(~ factor(view$design$group))
        sv <- quiet(sva::sva(a$y, mod, mod[, 1L, drop = FALSE], n.sv = a$E))$sv
        w <- matrix(sv, nrow = nrow(view$design), dimnames = list(view$design$Sample_ID, paste0("W", seq_len(a$E))))
        list(table = testWithW(view, assay, w), k = a$E)
    }),
    "15" = function(view) subtractArm(view, kE, replicateMatrix(stats::setNames(view$design$extract, view$design$Sample_ID))),
    "16" = function(view) modelArm(view, allComplete, kE,
        m_of = function(view, a) replicateMatrix(stats::setNames(view$design$extract, view$design$Sample_ID))),
    "17" = function(view) perAssay(view, function(a, assay) list(table = daKernel(a$y, view$design, "~ 0 + group + batch", view$contrasts), k = NA)),
    "22" = function(view) {
        with_standards <- lipidObject(lapply(view$assays, function(a) rbind(a$y, a$is)), view$design)
        k <- elbowK(view, lapply(view$assays, function(a) rownames(a$is)),
            quiet(prepareOmicsRuvCancorObject(with_standards, "lipid_data", "lipid_id_column")))
        modelArm(view, function(a, assay, columns) {
            if (is.na(k[[assay]])) return(refusal("findBestKElbow returned NA"))
            standardsOnly(a, assay, columns)
        }, function(a, assay) if (is.na(k[[assay]])) 0L else k[[assay]], standards = "is")
    },
    "23" = function(view) {
        controls <- lapply(view$first_pass(), function(fp) fp$feature[rank(fp$p, na.last = "keep", ties.method = "first") > 0.25 * sum(is.finite(fp$p)) & is.finite(fp$p)])
        usable <- all(vapply(controls, length, integer(1L)) > 0L)
        modelArm(view, function(a, assay, columns) {
            if (!usable) return(refusal("an assay has no first-pass controls"))
            columns %in% controls[[assay]]
        }, kE)
    },
    "25" = function(view) {
        checkExact(view)
        modelArm(view, allComplete, kE, exact = TRUE)
    },
    "29" = function(view) modelArm(view, allComplete, kVariant("E_units")),
    "30" = function(view) modelArm(view, allComplete, kVariant("E_units"), m_of = extractM),
    "31" = function(view) modelArm(view, allComplete, kVariant("E_200")),
    "32" = function(view) modelArm(view, allComplete, kVariant(if (hasDuplicates(view)) "E_units_200" else "E_200"), exact = TRUE),
    "33" = function(view) modelArm(view, allComplete, kVariant("E_units_200"), m_of = extractM, exact = TRUE),
    "34" = function(view) shippedArm(view, "groups_all_features"),
    "35" = function(view) shippedArm(view, "supplied_controls", standards = "is"),
    "36" = function(view) shippedArm(view, "supplied_controls", standards = "isc"),
    "37" = function(view) shippedArm(view, NULL),
    "38" = function(view) shippedArm(view, "groups_all_features", grouping = "extract"),
    "I" = function(view) {
        raw <- lapply(view$assays, function(a) 2^rbind(a$y, a$is))
        object <- quiet(logTransformAssays(normaliseUntransformedData(lipidObject(raw, view$design), method = "ITSD",
            itsd_aggregation = "median", itsd_feature_ids = lapply(view$assays, function(a) rownames(a$is))), offset = 1))
        perAssay(view, function(a, assay) {
            m <- assayMatrix(object@lipid_data[[assay]], view$design$Sample_ID)
            if (!identical(rownames(m), rownames(a$y))) stop("ITSD normalisation changed the feature set")
            list(table = daKernel(m, view$design, "~ 0 + group", view$contrasts), k = NA)
        })
    },
    "Qs" = function(view) subtractArm(view, kE, qcM(view$design, view$qc_ids), with_qc = TRUE),
    "Qm" = function(view) modelArm(view, allComplete, kE, m_of = function(view, a) qcM(view$design, view$qc_ids), with_qc = TRUE),
    "L" = function(view) perAssay(view, function(a, assay) {
        list(table = daKernel(loessCorrect(a$y, a$qc, view$design$pos, view$pos_q), view$design, "~ 0 + group", view$contrasts), k = NA)
    })
)

# The confirmation round (arms 34 to 38): the installed package's new-design kernels through both classes; the lipid
# table is scored, the metab table must agree with it. `ruv_design` NULL is no RUV.
shippedArm <- function(view, ruv_design, standards = NULL, grouping = NULL) {
    mats <- lapply(view$assays, function(a) if (is.null(standards)) a$y else rbind(a$y, a[[standards]]))
    controls <- if (!is.null(standards)) lapply(view$assays, function(a) rownames(a[[standards]]))
    contrasts <- contrastsTable(view$contrasts)
    through <- function(omic) {
        object <- if (omic == "lipid") lipidObject(mats, view$design) else metabObject(mats, view$design)
        factors <- if (!is.null(ruv_design)) quiet(runOmicsRuvFactors(object, ruv_design, ruv_k_rule = "num_sv_be",
            ruv_control_ids = controls, ruv_grouping_variable = grouping, ruv_k_seed = view$seed * 100L + 99L))$factors
        da <- quiet(if (omic == "lipid") runLipidsDA(object, contrasts, "~ 0 + group", 0.05, 0, TRUE, TRUE, ruv_factors = factors) else
            runMetabolitesDA(object, contrasts, "~ 0 + group", 0.05, 0, TRUE, TRUE, ruv_factors = factors))
        long <- if (omic == "lipid") da$da_lipids_long else da$da_metabolites_long
        id <- if (omic == "lipid") "lipid_id" else "metabolite_id"
        table <- data.frame(feature = long[[id]], contrast = long$comparison, est = long$logFC, p = long$P.Value,
            q = long$fdr_qvalue, bh = long$fdr_value_bh, assay = long$assay, stringsAsFactors = FALSE)
        list(table = table[order(table$assay, table$contrast, table$feature), , drop = FALSE], factors = factors)
    }
    lipid <- through("lipid")
    metab <- through("metab")
    a <- lipid$table
    b <- metab$table
    same <- identical(a, b) || (identical(a[c("feature", "contrast", "assay")], b[c("feature", "contrast", "assay")]) &&
        isTRUE(max(abs(c(a$est - b$est, a$p - b$p, a$q - b$q)), na.rm = TRUE) <= 1e-12) &&
        identical(is.na(c(a$est, a$p, a$q)), is.na(c(b$est, b$p, b$q))))
    if (!same) stop("the lipid and metab kernels give different tables on this data set")
    tables <- lapply(stats::setNames(ASSAYS, ASSAYS), function(assay) {
        record <- lipid$factors[[assay]]
        if (!is.null(record) && identical(record$status, "refused")) {
            return(refusal(sprintf("the package refused the assay (%s): %s", record$reason_class %||% "unclassed", record$reason %||% "")))
        }
        rows <- a[a$assay == assay, c("feature", "contrast", "est", "p", "q", "bh"), drop = FALSE]
        if (!nrow(rows)) stop("the package's DA returned no rows for assay ", assay)
        rownames(rows) <- NULL
        rows
    })
    list(tables = tables, k = vapply(ASSAYS, function(assay) as.numeric(lipid$factors[[assay]]$k %||% 0), numeric(1L)))
}

empiricalControlsArm <- function(view, controls) {
    usable <- vapply(controls, length, integer(1L)) > 0L
    k <- if (all(usable)) elbowK(view, controls) else NULL
    modelArm(view, function(a, assay, columns) {
        if (!all(usable)) return(refusal("an assay has no first-pass controls"))
        if (is.na(k[[assay]])) return(refusal("findBestKElbow returned NA"))
        columns %in% controls[[assay]]
    }, function(a, assay) if (is.null(k) || is.na(k[[assay]])) 0L else k[[assay]])
}

referenceArm <- function(view, method) {
    perAssay(view, function(a, assay) {
        k <- if (method == "RUV4") a$E else NA
        if (method == "RUV4" && is.na(k)) return(list(table = refusal("k could not be estimated"), k = 0))
        if (method == "RUV4" && k == 0L) return(list(table = testWithW(view, assay, NULL), k = 0))
        y <- completeY(a, a$is)
        ctl <- grepl("_ITSD_", colnames(y))
        n_groups <- length(unique(view$design$group))
        residual <- nrow(y) - n_groups
        if (method == "RUV4" && k > sum(ctl)) return(list(table = refusal(sprintf("k = %d above %d complete standards", k, sum(ctl))), k = 0))
        if (method == "RUV4" && residual - k < 1L) return(list(table = refusal(sprintf("no residual degrees of freedom with k = %d", k)), k = 0))
        if (method == "RUVinv" && sum(ctl) < residual) {
            return(list(table = refusal(sprintf("RUVinv needs %d controls (the residual dimensions), %d standards are complete", residual, sum(ctl))), k = NA))
        }
        table <- do.call(rbind, lapply(names(view$contrast_pairs), function(contrast) {
            pair <- view$contrast_pairs[[contrast]]
            others <- c(pair[[2L]], setdiff(unique(view$design$group), pair))
            x <- vapply(others, function(g) as.numeric(view$design$group == g), numeric(nrow(view$design)))
            fit <- if (method == "RUV4") ruv::RUV4(y, x, ctl, k, Z = 1) else ruv::RUVinv(y, x, ctl, Z = 1)
            p <- ruv::variance_adjust(fit)$p.ebayes[1L, !ctl]
            data.frame(feature = colnames(y)[!ctl], contrast = contrast, est = fit$betahat[1L, !ctl], p = p, q = pkgQ(p),
                bh = stats::p.adjust(p, "BH"), stringsAsFactors = FALSE)
        }))
        list(table = table, k = k)
    })
}

loessCorrect <- function(y, qc, pos_s, pos_q) {
    many <- length(pos_q) >= 10L
    for (j in seq_len(nrow(y))) {
        ok <- is.finite(qc[j, ])
        if (sum(ok) < 3L) next
        fit <- suppressWarnings(stats::loess(v ~ x, data.frame(x = pos_q[ok], v = qc[j, ok]), span = if (many) 0.75 else 1.5,
            degree = if (many) 2 else 1, control = stats::loess.control(surface = "direct")))
        y[j, ] <- y[j, ] - stats::predict(fit, data.frame(x = pos_s)) + stats::median(qc[j, ok])
    }
    y
}

# ---- one view (a subset of the features) with its shared, lazily computed pieces ------------------------------------
makeView <- function(data, estimates, variants = list()) {
    view <- data
    for (assay in ASSAYS) {
        view$assays[[assay]]$complete <- completeFeatures(view$assays[[assay]])
        view$assays[[assay]]$E <- estimates[[assay]]
        for (name in names(variants)) view$assays[[assay]][[name]] <- variants[[name]][[assay]]
    }
    cache <- new.env(parent = emptyenv())
    cache$seconds <- c(base = 0, sweep = 0, cancor = 0)
    lazy <- function(name, expr_fn) function() {
        if (is.null(cache[[name]])) {
            cache$seconds[[name]] <- system.time(cache[[name]] <- expr_fn())[["elapsed"]]
        }
        cache[[name]]
    }
    view$object <- function() lipidObject(lapply(view$assays, `[[`, "y"), view$design)
    view$base <- lazy("base", function() lapply(view$assays, function(a) daKernel(a$y, view$design, "~ 0 + group", view$contrasts)))
    view$first_pass <- function() lapply(view$base(), function(t) {
        # one row per feature: the smallest p and the smallest BH q over the contrasts
        stats::aggregate(cbind(p, bh) ~ feature, t, min, na.action = stats::na.pass)
    })
    view$sweep <- lazy("sweep", function() quiet(runLipidPerAssayRuvOptimization(view$object(), "automatic", RUV_PARAMS)))
    view$cancor_object <- lazy("cancor", function() quiet(prepareOmicsRuvCancorObject(view$object(), "lipid_data", "lipid_id_column")))
    view$cache <- cache
    view
}

# ---- scoring ---------------------------------------------------------------------------------------------------
scoreArm <- function(result, view) {
    tables <- result$tables
    refused <- vapply(tables, isRefusal, logical(1L))
    base <- if (any(refused)) view$base()
    rows <- do.call(rbind, lapply(ASSAYS, function(assay) {
        t <- if (refused[[assay]]) base[[assay]] else tables[[assay]]
        truth <- view$assays[[assay]]$truth
        t <- t[t$feature %in% rownames(truth), ]
        t$truth <- truth[cbind(match(t$feature, rownames(truth)), match(t$contrast, colnames(truth)))]
        t
    }))
    n_true <- sum(vapply(view$assays, function(a) sum(a$truth != 0), numeric(1L)))
    tested <- is.finite(rows$p)
    null <- tested & rows$truth == 0
    real <- tested & rows$truth != 0
    call <- tested & is.finite(rows$q) & rows$q < 0.05
    error <- rows$est[real] - rows$truth[real]
    list(status = if (any(refused)) "refused" else "ok", k = mean(result$k), n_features = sum(tested), n_true = n_true,
        n_calls = sum(call), null_p05 = if (any(null)) mean(rows$p[null] < 0.05) else NA_real_,
        null_p01 = if (any(null)) mean(rows$p[null] < 0.01) else NA_real_,
        fdp = if (any(call)) sum(call & rows$truth == 0) / sum(call) else 0,
        power = if (n_true > 0) sum(call & rows$truth != 0) / n_true else NA_real_,
        bias = if (any(real)) mean(sign(rows$truth[real]) * error) else NA_real_,
        rmse = if (any(real)) sqrt(mean(error^2)) else NA_real_,
        reasons = vapply(tables[refused], function(x) x$reason, ""))
}

COLUMNS <- c("simulator", "cell", "cell_class", "run", "seed", "arm", "subset", "status", "k", "n_features", "n_true",
    "n_calls", "null_p05", "null_p01", "fdp", "power", "bias", "rmse", "seconds")

runOne <- function(cell, run) {
    seed <- 100000L + 10000L * cell$index + run
    started <- proc.time()[["elapsed"]]
    data <- generateData(cell, seed)
    seconds_E <- system.time(estimates <- lapply(data$assays, function(a) {
        dat <- a$y[completeFeatures(a), , drop = FALSE]
        if (nrow(dat) < 20L) return(NA_integer_)
        tryCatch(as.integer(quiet(sva::num.sv(dat, stats::model.matrix(~ factor(data$design$group)), method = "be", B = 20,
            seed = seed * 100L + 99L))), error = function(e) NA_integer_)
    }))[["elapsed"]]
    # the exactness round: E with 200 permutations and with the extracts (the biological units) removed, E's seed
    estimateWith <- function(B, units) lapply(data$assays, function(a) {
        dat <- a$y[completeFeatures(a), , drop = FALSE]
        if (nrow(dat) < 20L) return(NA_integer_)
        mod <- stats::model.matrix(~ factor(if (units) data$design$extract else data$design$group))
        tryCatch(as.integer(quiet(sva::num.sv(dat, mod, method = "be", B = B, seed = seed * 100L + 99L))), error = function(e) NA_integer_)
    })
    variants <- list()
    if (any(c("31", "32") %in% armsOfCell(cell))) variants$E_200 <- estimateWith(200, FALSE)
    if (any(c("29", "30") %in% armsOfCell(cell))) variants$E_units <- estimateWith(20, TRUE)
    if (any(c("32", "33") %in% armsOfCell(cell)) && cell$duplicates) variants$E_units_200 <- estimateWith(200, TRUE)
    keep <- lapply(data$assays, completeFeatures)
    same <- all(unlist(keep))
    views <- list(full = makeView(data, estimates, variants), complete = if (!same) makeView(restrictData(data, keep), estimates, variants))
    rows <- list()
    notes <- list()
    for (subset in c("full", "complete")) {
        for (arm in armsOfCell(cell)) {
            if (subset == "full" && arm %in% COMPLETE_ONLY && !same) next
            key <- paste(arm, subset)
            if (subset == "complete" && same) {
                rows[[key]] <- rows[[paste(arm, "full")]]
                rows[[key]]$subset <- "complete"
                next
            }
            view <- views[[subset]]
            set.seed(seed * 100L + if (arm %in% names(EXACT_STREAM)) EXACT_STREAM[[arm]] else ARM_NUMBER[[arm]])
            t0 <- proc.time()[["elapsed"]]
            scored <- tryCatch({
                if (identical(arm, fault_arm) && run == first_run) stop("deliberate fault (K124_FAULT, the lifecycle test)")
                scoreArm(armFunctions[[arm]](view), view)
            }, error = function(e) list(status = "crashed", reasons = conditionMessage(e)))
            seconds <- proc.time()[["elapsed"]] - t0
            measures <- c("k", "n_features", "n_true", "n_calls", "null_p05", "null_p01", "fdp", "power", "bias", "rmse")
            values <- lapply(stats::setNames(measures, measures), function(m) if (is.null(scored[[m]])) NA_real_ else signif(scored[[m]], 6L))
            rows[[key]] <- data.frame(simulator = "S4", cell = cell$id, cell_class = cell$class, run = run, seed = seed, arm = arm,
                subset = subset, status = scored$status, values, seconds = round(seconds, 3), stringsAsFactors = FALSE)
            if (scored$status != "ok") {
                notes[[key]] <- data.frame(cell = cell$id, run = run, arm = arm, subset = subset, status = scored$status,
                    assay = if (is.null(names(scored$reasons))) "" else names(scored$reasons),
                    reason = gsub("[\t\n]+", " ", unname(scored$reasons)), stringsAsFactors = FALSE)
            }
        }
    }
    rows <- do.call(rbind, rows)
    # the shared steps an arm uses are part of its cost: they ran inside the first arm that asked for them
    for (subset in names(views)) {
        view <- views[[subset]]
        if (is.null(view)) next
        s <- view$cache$seconds
        at <- function(arm) rows$subset == subset & rows$arm == arm
        own <- function(arm) rows$seconds[at(arm)]
        if (any(at("3"))) rows$seconds[at("4")] <- own("4") + s[["sweep"]]
        if (any(at("5"))) rows$seconds[at("6")] <- own("6") + s[["cancor"]]
        for (arm in setdiff(unique(rows$arm), c("1", "2", "3", "4", "5", "6", "11", "13", "17", "22", "I", "L"))) {
            rows$seconds[at(arm)] <- own(arm) + seconds_E
        }
    }
    if (same) {
        # no missing values: the complete rows repeat the full ones at no cost; arms 12 to 14 keep their one row
        rows$seconds[rows$subset == "complete"] <- ifelse(rows$arm[rows$subset == "complete"] %in% COMPLETE_ONLY,
            rows$seconds[rows$subset == "full"], 0)
        rows <- rows[!(rows$subset == "full" & rows$arm %in% COMPLETE_ONLY), ]
    }
    rownames(rows) <- NULL
    rows$seconds <- round(rows$seconds, 3)
    list(rows = rows[COLUMNS], notes = if (length(notes)) do.call(rbind, notes),
        timing = data.frame(cell = cell$id, run = run, seconds_total = round(proc.time()[["elapsed"]] - started, 3),
            seconds_E = round(seconds_E, 3),
            missing_share = round(mean(vapply(data$assays, function(a) mean(is.na(a$y)), numeric(1L))), 4),
            complete_features = sum(unlist(keep)), stringsAsFactors = FALSE))
}

# ---- driver ----------------------------------------------------------------------------------------------------
writeTsv <- function(x, path) {
    tmp <- paste0(path, ".tmp")
    utils::write.table(x, tmp, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    file.rename(tmp, path)
}
readTsv <- function(path) utils::read.delim(path, stringsAsFactors = FALSE, colClasses = "character", quote = "")

main <- function() {
    cells <- cellTable()
    cells <- if (identical(cells_arg, "all")) cells[setdiff(names(cells), "C8")] else cells[strsplit(cells_arg, ",")[[1L]]]
    dir.create(file.path(out_dir, "cells"), recursive = TRUE, showWarnings = FALSE)
    started <- Sys.time()
    checks <- checkSharedKernels()
    for (cell in cells) {
        path <- file.path(out_dir, "cells", paste0(cell$id, ".tsv"))
        if (file.exists(path)) next
        runs <- first_run - 1L + seq_len(if (cell$class == "core") runs_core else runs_stress)
        t0 <- proc.time()[["elapsed"]]
        results <- parallel::mclapply(runs, function(run) tryCatch(runOne(cell, run), error = function(e) conditionMessage(e)),
            mc.cores = cores, mc.preschedule = FALSE)
        failed <- !vapply(results, is.list, logical(1L))
        notes <- do.call(rbind, lapply(results[!failed], `[[`, "notes"))
        if (any(failed)) {
            # a run that died outside an arm: every arm of it is a crash
            reasons <- vapply(results[failed], function(x) if (is.character(x)) x[[1L]] else "worker died", "")
            notes <- rbind(notes, data.frame(cell = cell$id, run = runs[failed], arm = "all", subset = "all", status = "crashed",
                assay = "", reason = gsub("[\t\n]+", " ", reasons), stringsAsFactors = FALSE))
        }
        rows <- do.call(rbind, c(lapply(results[!failed], `[[`, "rows"), lapply(runs[failed], function(run) {
            arms <- armsOfCell(cell)
            grid <- rbind(expand.grid(arm = setdiff(arms, COMPLETE_ONLY), subset = "full", stringsAsFactors = FALSE),
                expand.grid(arm = arms, subset = "complete", stringsAsFactors = FALSE))
            empty <- as.list(stats::setNames(rep(NA_real_, 11L), COLUMNS[9:19]))
            data.frame(simulator = "S4", cell = cell$id, cell_class = cell$class, run = run,
                seed = 100000L + 10000L * cell$index + run, arm = grid$arm, subset = grid$subset, status = "crashed", empty,
                stringsAsFactors = FALSE)[COLUMNS]
        })))
        timing <- do.call(rbind, lapply(results[!failed], `[[`, "timing"))
        if (!is.null(notes)) writeTsv(notes, file.path(out_dir, "cells", paste0(cell$id, ".notes.tsv")))
        if (!is.null(timing)) writeTsv(timing, file.path(out_dir, "cells", paste0(cell$id, ".timing.tsv")))
        writeTsv(rows, path)
        cat(sprintf("%s %s: %d runs, %.0f s wall, %d refused, %d crashed rows\n", format(Sys.time(), "%H:%M:%S"), cell$id,
            length(runs), proc.time()[["elapsed"]] - t0, sum(rows$status == "refused"), sum(rows$status == "crashed")))
    }
    gather <- function(pattern) {
        files <- file.path(out_dir, "cells", sprintf(pattern, names(cellTable())))
        do.call(rbind, lapply(files[file.exists(files)], readTsv))
    }
    writeTsv(gather("%s.tsv"), file.path(out_dir, "s4-runs.tsv"))
    notes <- gather("%s.notes.tsv")
    writeTsv(if (is.null(notes)) data.frame(cell = character(), run = character(), arm = character(), subset = character(),
        status = character(), assay = character(), reason = character()) else notes, file.path(out_dir, "s4-notes.tsv"))
    writeTsv(gather("%s.timing.tsv"), file.path(out_dir, "s4-timing.tsv"))
    git <- function(...) tryCatch(system2("git", c(...), stdout = TRUE, stderr = FALSE), error = function(e) NA_character_)
    versions <- vapply(c("limma", "qvalue", "ruv", "RUVIIIC", "sva", "statmod", "pkgload"), function(p) as.character(utils::packageVersion(p)), "")
    writeLines(c("simulator: S4 (lipidomics and metabolomics), KIT-124 revision 1.2; KIT-129 R0: the extract declared as technical replicate (kit-128 tree)",
        paste("commit:", git("rev-parse", "HEAD")),
        paste("script committed state:", if (length(git("status", "--porcelain", "--", "docs/qa/workflow-execution/kit/kit-124/sim/s4_lipid_metab.R"))) "modified or untracked at run time" else "clean"),
        paste("command:", paste(c("OMP_NUM_THREADS=1 Rscript", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), args), collapse = " ")),
        paste("cells:", paste(names(cells), collapse = ",")),
        paste("K124_ARMS:", if (is.null(arms_wanted)) "unset (the arms of the header's list but 22 and 23)" else paste(arms_wanted, collapse = ",")),
        paste("runs:", sprintf("core %d, stress %d, numbered from %d", runs_core, runs_stress, first_run)),
        paste("cores:", cores, "(OMP_NUM_THREADS=1)"),
        "seed rule: data seed = 100000 + 10000 * cell index + run (C0..C7 = 0..7, S1..S14 = 8..21); arm stream = set.seed(seed * 100 + arm number), I = 18, Qs = 19, Qm = 20, L = 21 (arms 25 to 33: the stream of the arm each varies); E = seed * 100 + 99; cell index of C8 = 22",
        paste("R:", R.version.string),
        paste("library paths:", paste(.libPaths(), collapse = ":")),
        paste("R_LIBS:", Sys.getenv("R_LIBS")),
        paste(if (installed) sprintf("MultiScholaR (installed, K124_PACKAGE=installed; built from %s, at %s):",
            Sys.getenv("K124_PACKAGE_COMMIT", "unknown"), find.package("MultiScholaR")) else "MultiScholaR (pkgload::load_all):",
            as.character(utils::packageVersion("MultiScholaR"))),
        paste("K124_FAULT:", if (nzchar(fault_arm)) fault_arm else "unset"),
        paste("packages:", paste(names(versions), versions, collapse = ", ")),
        paste("lipid/metab sharing check (identical outputs through both classes):", paste(names(checks), unlist(checks), collapse = ", ")),
        paste("started:", format(started, "%Y-%m-%d %H:%M:%S"), "finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))),
        file.path(out_dir, "s4-runs.meta.txt"))
}

main()
