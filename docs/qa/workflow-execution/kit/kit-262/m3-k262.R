# KIT-262 M3 (ticket section 6; VAL-302 O-1 (ii), O-2 (a), O-3 yes, O-4 (a)): the present rule P against candidate C1 (R5 judged at
# min(k, cap); ranking unchanged), paired on the same data set, through the product's own entry points. NOT C2 (not tested in
# M3, VAL-302 O-1 (ii)). Simulated data. Evidence code: no package change. DRAFT: written, not sent; an auditor reads it first.
#
# Usage, from the package root on a run-tests runner (package installed, full history present):
#   Rscript docs/qa/workflow-execution/kit/kit-262/m3-k262.R <out dir relative to the package root>
# One UNIT = one data set (cell, r): the P journey (the shipped Automatic RUV step, KIT-172's `kit172Journey`), then, only where
# C1's choice differs from P's, one Manual journey at C1's percentage and k (no package code: Manual applies rule (d) and no R5,
# the bound is inert there). Where C1's choice equals P's the unit re-uses P's row (identical input).
# Environment:
#   K262_M3_CELLS   cells to run, comma list among N-r,N-s,E-r,E-s,E-p,E-g,N-w,E-w,S-N,S-E (default: the eight of VAL-235)
#   K262_M3_RMAP    replicates per cell "N-r=100,N-s=100,..." (default: VAL-235's final R for the eight cells; a stress cell
#                   has NO default: its R comes from the blinded re-estimation and is declared with the counted job)
#   K262_M3_REPS    replicates for every selected cell, e.g. "9001,9002" or "1:50" (overrides K262_M3_RMAP)
#   K262_M3_E9REPS  replicates, e.g. "1,2": extra P-only units on the eight cells for E9 (rows kept like the others)
#   K262_M3_E8      "1" (default in L3 and counted runs): units with r <= K262_M3_E8_MAX (default 5) also run as metabolomics (E8)
#   K262_M3_E11MOD  replicate rule of the E11 sample (default 20: r %% 20 == 7, 5 % of the replicates)
#   K262_M3_WORKERS forked workers (default 4)
#   K262_M3_ONLY_UNITS  CONTINUATION (counted mode only): run exactly the units named (tags `cNN-rR`), once; a tag outside the plan stops
#   K262_M3_FORCE_ARMS "1" (L3): every unit whose C1 re-uses P's row also runs the Manual journey at P's own choice (E11 forced)
#   K262_M3_STAGE1  "1": blind stage 1 of the re-estimation (cells S-N, S-E, replicates 1:50): rows stay under tempdir(); only the
#                   reestimate.tsv of `reestimate-m3-k262.R` and fixed console lines leave
#   K262_M3_GATE    path (from the package root) of the checks.tsv of the scorer's self-test run earlier in the same request: when set, the script
#                   stops before any unit unless that file exists and every row passes (the runner does not stop a request after a failed script)
#   K262_LIFECYCLE  "1": row files under tempdir(); only checks.tsv, timing.tsv and hashes.tsv leave (VAL-11)
#   K262_M3_FAULT   "1" (L3 only): units 1 to 5 carry planted faults: 1 a journey stop (a data-set error row), 2 a partial match inside a
#                   journey, 3 a partial match in unit code outside a journey (both instrument rows), 4 a partial match inside an inner
#                   tryCatch(error =) in unit code (still caught by the guard: instrument row, class partial_match_guard), 5 a worker that dies with no row
#                   and no trace (an instrument row written by the parent). The checks FAULT_ROWS then read the five outcomes and that the sixth unit is clean.
# Rows: `summary.tsv` (one line per unit, appended by the unit itself when it ends), `traces/trace-<tag>.txt` for an error or
# instrument unit, `provenance.json` (at start and at end), `checks.tsv`, `timing.tsv`, `hashes.tsv`. A unit with no row and no
# trace gets an instrument row from the parent; the scorer stops on every instrument row.
suppressPackageStartupMessages(pkgload::load_all(".", export_all = TRUE, helpers = FALSE, quiet = TRUE))
logger::log_threshold(logger::FATAL)
options(warnPartialMatchDollar = TRUE)
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- the partial-match guard (lesson of KIT-188 and KIT-275) ---------------------------------------------------------------
# `$` on a list or data frame matches a name by its prefix when none is exact; that warning is an error here. It is a LOCAL calling
# handler around each journey and each unit (a stop from a GLOBAL calling handler runs with every local handler disabled and ends the
# process with nothing written). The enclosing tryCatch receives the stop; the journey's and the unit's own handlers turn it into an
# instrument row. Every read of a field that may be absent below uses `[[` (exact matching).
m3GuardMuffle <- function(w) {
  if (grepl("partial match of", tolower(conditionMessage(w)), fixed = TRUE)) {
    options(k262.guard_trace = paste(vapply(sys.calls(), function(cl) substr(paste(deparse(cl)[[1L]], collapse = ""), 1L, 80L), ""),
      collapse = " <- "))
    stop(paste("partial-match guard:", conditionMessage(w)), call. = FALSE)
  }
  invokeRestart("muffleWarning")
}
m3IsGuard <- function(e) grepl("partial-match guard", conditionMessage(e), fixed = TRUE)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1L)
out <- normalizePath(file.path(getwd(), args[[1L]]), mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
# The gate (counted requests): the scorer's self-test must have passed in this request, else no unit runs.
gate <- Sys.getenv("K262_M3_GATE")
if (nzchar(gate)) {
  g <- tryCatch(utils::read.delim(gate, stringsAsFactors = FALSE, colClasses = "character"), error = function(e) NULL)
  if (is.null(g) || !nrow(g) || !"pass" %in% names(g) || !all(g$pass %in% "TRUE")) stop("K262 M3 gate: the scorer's self-test did not pass (", gate, "); no unit was run")
}
source("docs/qa/workflow-execution/kit/kit-262/rules-k262.R")

# The test helpers of this tree in their own environment (the harness of KIT-172 and KIT-190, read-only), attached behind the global one.
helpers <- new.env(parent = globalenv())
suppressMessages(testthat::source_test_helpers(file.path("tests", "testthat"), env = helpers))
attach(helpers, name = "k262_helpers", warn.conflicts = FALSE)
# [M] The helper's own `kit172Quiet` (suppressWarnings(suppressMessages(.))) would hide a partial-match warning from the guard; inside the
# helper's enclosure it is replaced by the guard-aware quiet below. No other helper is changed.
m3Quiet <- function(expr) withCallingHandlers(expr, warning = m3GuardMuffle, message = function(m) invokeRestart("muffleMessage"))  # [M]
assign("kit172Quiet", m3Quiet, envir = helpers)  # [M]

# ---- copy: tests/testthat/test-workflow-omics-ruv-cancor-once.R:43-54 is not needed here; M3 copies the following blocks ----
# ---- copy: docs/qa/workflow-execution/kit/kit-190/effect.R e2PickCells, verbatim ----
e2PickCells <- function(allowed, groups, per_group, count) {
  cells <- which(allowed, arr.ind = TRUE)
  cells <- cells[sample.int(nrow(cells)), , drop = FALSE]
  key <- paste(cells[, 1L], groups[cells[, 2L]])
  rank <- stats::ave(seq_along(key), key, FUN = seq_along)
  cells <- cells[rank <= per_group, , drop = FALSE]
  cells[seq_len(min(count, nrow(cells))), , drop = FALSE]
}
# ---- end of the copy ----

# ---- copy: docs/qa/workflow-execution/kit/kit-190/effect.R e2Values (VAL-235's eight cells), verbatim except the [L] lines ----
# [L] the true level per feature is returned (the abundance thirds of the calibration need it); no draw moves.
e2Values <- function(seed, planted = FALSE, missing = TRUE, n = 520L, n_itsd = 5L, effect_share = 0.1,
                     factor_pattern = c("balanced", "random")) {
  factor_pattern <- match.arg(factor_pattern)
  set.seed(seed)
  design <- kit172Design()
  groups <- design$group
  n_features <- n - n_itsd
  level <- stats::runif(n_features, log2(0.01), log2(1e6))
  logs <- matrix(level, n_features, nrow(design)) + matrix(stats::rnorm(n_features * nrow(design), 0, 0.25), n_features)
  effects <- rep(0, n_features)
  loadings <- rep(0, n_features)
  if (planted) {
    changed <- sample.int(n_features, round(effect_share * n_features))
    effects[changed] <- sample(c(-1, 1), length(changed), replace = TRUE)
    logs[, groups == "B"] <- logs[, groups == "B"] + effects
    loadings <- stats::rnorm(n_features, 0, 0.5)
    pattern <- rep(c(-1, -1, -1, 1, 1, 1), 2L)
    if (identical(factor_pattern, "random")) {
      w <- stats::rnorm(nrow(design))
      pattern <- (w - mean(w)) / stats::sd(w) * sqrt(12 / 11)
    }
    logs <- logs + outer(loadings, pattern)
  }
  values <- 2^logs
  if (missing) {
    gone <- e2PickCells(matrix(TRUE, n_features, nrow(design)), groups, 3L, round(0.1 * length(values)))
    values[gone] <- NA_real_
  }
  ones <- e2PickCells(!is.na(values), groups, 2L, round(0.02 * sum(!is.na(values))))
  values[ones] <- 1
  standards <- matrix(2^(13 + stats::rnorm(n_itsd * nrow(design), 0, 0.1)), n_itsd)
  out <- rbind(values, standards)
  dimnames(out) <- list(c(sprintf("F%03d", seq_len(n_features)), sprintf("IS_%d", seq_len(n_itsd))), design$Run)
  list(values = out, effects = stats::setNames(c(effects, rep(0, n_itsd)), rownames(out)),
    loadings = stats::setNames(c(loadings, rep(0, n_itsd)), rownames(out)),
    level = stats::setNames(c(level, rep(NA_real_, n_itsd)), rownames(out)),  # [L]
    factor_group_r = if (planted) stats::cor(pattern, as.numeric(groups == "B")) else NA_real_)
}
# ---- end of the copy ----

# ---- copy: docs/qa/workflow-execution/kit/kit-190/effect.R E2_CELLS, ALIGNED_CELLS, e2Data (the seeds 190000000 + 10000 x idx + r) ----
E2_CELLS <- c("N-r", "N-s", "E-r", "E-s", "E-p", "E-g")
ALIGNED_CELLS <- c("N-w", "E-w")
e2Data <- function(cell, r) {
  idx <- match(cell, c(E2_CELLS, ALIGNED_CELLS))
  seed <- 190000000 + 10000 * idx + r
  sim <- e2Values(seed, planted = TRUE, missing = TRUE, effect_share = if (startsWith(cell, "N")) 0 else 0.1,
    factor_pattern = if (cell %in% ALIGNED_CELLS) "random" else "balanced")
  values <- sim$values
  effects <- sim$effects
  one_group <- character()
  if (identical(cell, "E-g")) {
    pool <- which(!startsWith(rownames(values), "IS_") & effects == 0)
    one_group <- rownames(values)[pool[sample.int(length(pool), 26L)]]
    values[one_group, kit172Design()$group == "A"] <- NA_real_
  }
  ord <- switch(cell,
    "N-s" = , "E-s" = order(rowMeans(log2(values), na.rm = TRUE)),
    "E-p" = order(effects == 0),
    seq_len(nrow(values)))
  list(values = values[ord, , drop = FALSE], effects = effects, one_group = one_group, seed = seed,
    factor_group_r = sim$factor_group_r, level = sim$level)  # [L] `level` added to the returned list
}
# ---- end of the copy ----

# ---- copy: docs/qa/workflow-execution/kit/kit-262/m2-k262.R m2Design (S01.. in g equal groups), verbatim ----
m2Design <- function(n, g) {
  runs <- sprintf("S%02d", seq_len(n))
  group <- rep(LETTERS[seq_len(g)], each = n / g)
  data.frame(Run = runs, group = group, batch = "B1", replicates = stats::ave(seq_len(n), group, FUN = seq_along),
    tech_rep_group = runs, stringsAsFactors = FALSE)
}
# ---- end of the copy ----

# ---- the stress generator (ticket 6.1): M2's `m2Values` (a declared generalisation of KIT-190's e2Values) in KIT-172's floor ORDER ----
# The floor order is named and fixed: `after` (the mask first, then 2 % of the OBSERVED values set to 1), KIT-172's, the ticket's and the
# order of VAL-235's eight cells (ruling, m3-plan-notes.md 1). Marked lines: [N] the order named and asserted; [H] `effect_share` (0 for S-N,
# 0.1 for S-E, as e2Values' own argument); [I] MCAR only (the stress cells' mask), the abundance-dependent branch removed; [L] the true
# level is returned; [P] fixed design n = 24, g = 2, 300 features, MCAR 10 % (M2's cell 79).
M3_STRESS_FLOOR_ORDER <- "after"
m3StressValues <- function(seed, effect_share, floor_order = M3_STRESS_FLOOR_ORDER, n = 24L, g = 2L, n_features = 300L, rate = 0.10) {  # [N] [H] [P]
  stopifnot(identical(floor_order, "after"))  # [N]
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  design <- m2Design(n, g)
  groups <- design$group
  level <- stats::runif(n_features, log2(0.01), log2(1e6))
  logs <- matrix(level, n_features, nrow(design)) + matrix(stats::rnorm(n_features * nrow(design), 0, 0.25), n_features)
  effects <- rep(0, n_features)
  loadings <- rep(0, n_features)
  changed <- sample.int(n_features, round(effect_share * n_features))  # [H]
  effects[changed] <- sample(c(-1, 1), length(changed), replace = TRUE)
  logs[, groups == LETTERS[g]] <- logs[, groups == LETTERS[g]] + effects
  loadings <- stats::rnorm(n_features, 0, 0.5)
  w <- stats::rnorm(nrow(design))
  pattern <- (w - mean(w)) / stats::sd(w) * sqrt(n / (n - 1))
  logs <- logs + outer(loadings, pattern)
  values <- 2^logs
  values[sample.int(length(values), floor(rate * length(values)))] <- NA_real_  # [I] MCAR, the mask first
  ones <- e2PickCells(!is.na(values), groups, 2L, round(0.02 * sum(!is.na(values))))  # [N] then the floor values, on observed cells
  values[ones] <- 1  # [N]
  ids <- sprintf("F%04d", seq_len(n_features))
  dimnames(values) <- list(ids, design$Run)
  list(values = values, design = design, effects = stats::setNames(effects, ids), level = stats::setNames(level, ids))  # [L]
}

# ---- the design-aware copies of the harness (ticket 6.1, P20): kit172Object and kit172StateMatrix are fixed at 12 samples ----
# ---- copy: tests/testthat/helper-kit172.R kit172Object, with the design a parameter [D] ----
m3Object <- function(values, class, design = kit172Design()) {  # [D]
  # [D] the design is the parameter above; the source's `design <- kit172Design()` is removed
  if (!is.list(values)) values <- list(LCMS_Pos = values)
  assays <- lapply(values, kit172Table, class = class)
  args <- list(deAnalysisParameters = list(formula_string = "~ 0 + group", de_q_val_thresh = 0.05, treat_lfc_cutoff = 0))
  if (identical(class, "lipid")) {
    methods::new("LipidomicsAssayData", lipid_data = assays, design_matrix = design, lipid_id_column = "LipidName",
      annotation_id_column = "LipidClass", sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
      database_identifier_type = "Unknown", internal_standard_regex = "^IS_", args = args)
  } else {
    methods::new("MetaboliteAssayData", metabolite_data = assays, design_matrix = design, metabolite_id_column = "Metabolite",
      annotation_id_column = "Ontology", sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
      database_identifier_type = "Unknown", internal_standard_regex = "^IS_", args = args)
  }
}
# ---- end of the copy ----
# ---- copy: tests/testthat/helper-kit172.R kit172StateMatrix, with the runs a parameter [D] ----
m3StateMatrix <- function(journey, state, runs, assay = "LCMS_Pos") {  # [D]
  # [D] the runs are the parameter above; the source's `runs <- kit172Design()$Run` is removed
  object <- journey$owner$getState(state)
  table <- methods::slot(object, journey$binding$data_slot)[[assay]]
  id <- methods::slot(object, journey$binding$id_slot)
  matrix(as.matrix(as.data.frame(table)[, runs, drop = FALSE]), ncol = length(runs), dimnames = list(table[[id]], runs))
}
# ---- end of the copy ----

# ---- copy: docs/qa/workflow-execution/kit/kit-190/effect.R CUTOFF, k190Oracle, k190Candidate (rule (d)'s oracle, F2), verbatim ----
CUTOFF <- 0.05
k190Oracle <- function(m, groups, pct) {
  m <- m[rowSums(!is.na(m)) > 1, , drop = FALSE]
  count <- round(nrow(m) * pct / 100)
  p <- vapply(rownames(m), function(id) {
    x <- m[id, ]
    observed <- !is.na(x)
    g <- factor(groups[observed])
    if (nlevels(g) < 2L) return(1)
    fit <- stats::aov(x[observed] ~ g)
    if (fit$df.residual == 0L) return(NA_real_)
    value <- summary(fit)[[1]][["Pr(>F)"]][1]
    if (is.na(value)) 1 else value
  }, numeric(1L))
  one_group_ids <- rownames(m)[apply(!is.na(m), 1L, function(x) length(unique(groups[x])) < 2L)]
  p <- p[!is.na(p)]
  n <- length(p)
  o <- order(p)
  q <- numeric(n)
  q[o] <- pmin(1, rev(cummin(rev(p[o] * n / seq_len(n)))))
  names(q) <- names(p)
  candidates <- names(q)[q > CUTOFF]
  ranked <- candidates[order(-q[candidates])]
  chosen <- utils::head(ranked, max(min(count, length(ranked)), 1L))
  list(selected = rownames(m)[rownames(m) %in% chosen], q = q, p = p, count = count, valid = rownames(m),
    one_group = one_group_ids)
}
k190Candidate <- function(orc, rule) {
  candidates <- names(orc$q)[orc$q > CUTOFF]
  if (rule %in% c("c", "d")) candidates <- setdiff(candidates, orc$one_group)
  ranked <- if (rule %in% c("b", "d")) candidates[order(-orc$p[candidates], candidates, method = "radix")] else
    candidates[order(-orc$q[candidates])]
  if (!length(ranked)) return(character())
  utils::head(ranked, max(min(orc$count, length(ranked)), 1L))
}
# ---- end of the copy ----

# ---- copy: docs/qa/workflow-execution/kit/kit-262/m2-k262.R m2Object, M2_RUV_PARAMS, m2Prepare, m2Sweep (E10's direct M2-form sweep), verbatim ----
m2Object <- function(values, design, omic) {
  table <- if (identical(omic, "lipid")) data.frame(LipidName = rownames(values), LipidClass = "PC", stringsAsFactors = FALSE) else
    data.frame(Metabolite = rownames(values), Ontology = "PC", stringsAsFactors = FALSE)
  for (j in seq_len(ncol(values))) table[[colnames(values)[[j]]]] <- unname(values[, j])
  assays <- list(LCMS_Pos = table)
  args <- list(deAnalysisParameters = list(formula_string = "~ 0 + group", de_q_val_thresh = 0.05, treat_lfc_cutoff = 0))
  if (identical(omic, "lipid")) {
    methods::new("LipidomicsAssayData", lipid_data = assays, design_matrix = design, lipid_id_column = "LipidName",
      annotation_id_column = "LipidClass", sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
      database_identifier_type = "Unknown", internal_standard_regex = "^IS_", args = args)
  } else {
    methods::new("MetaboliteAssayData", metabolite_data = assays, design_matrix = design, metabolite_id_column = "Metabolite",
      annotation_id_column = "Ontology", sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
      database_identifier_type = "Unknown", internal_standard_regex = "^IS_", args = args)
  }
}

M2_RUV_PARAMS <- list(percentage_min = 1, percentage_max = 20, ruv_grouping_variable = "group",
  separation_metric = "max_difference", k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, manual_k = 2L, manual_percentage = 5)

m2Prepare <- function(object, omic) {
  binding <- if (identical(omic, "lipid")) lipidomicsDeclarationBinding else metabolomicsDeclarationBinding
  parameters <- stats::setNames(list(1, 0.5, 2L, 1L), binding$intensity_parameters)
  filtered <- suppressWarnings(suppressMessages(binding$filter_intensity(object, parameters)))
  list(binding = binding, object = suppressWarnings(suppressMessages(logTransformAssays(theObject = filtered, offset = 0))),
    filtered = filtered)
}
m2Sweep <- function(prepared) {
  suppressWarnings(suppressMessages(prepared$binding$optimise_ruv(prepared$object, "automatic", M2_RUV_PARAMS)))
}
# ---- end of the copy ----

# ---- the cells (ticket 6.1) ----------------------------------------------------------------------------------------------
# idx: KIT-190's cell index for the eight (seeds 190000000 + 10000 x idx + r, VAL-235's seeds verbatim); 9 and 10 for the stress cells
# (seeds 262200000 + 10000 x idx + r). kind: null (any-call endpoint) or effect (false-discovery proportion).
M3_CELLS <- data.frame(id = c("N-r", "N-s", "E-r", "E-s", "E-p", "E-g", "N-w", "E-w", "S-N", "S-E"), idx = c(1:8, 9L, 10L),
  kind = c("null", "null", "effect", "effect", "effect", "effect", "null", "effect", "null", "effect"),
  stress = c(rep(FALSE, 8L), TRUE, TRUE), stringsAsFactors = FALSE)
M3_FINAL_R <- c(`N-r` = 100L, `N-s` = 100L, `E-r` = 200L, `E-s` = 200L, `E-p` = 300L, `E-g` = 200L, `N-w` = 250L, `E-w` = 200L)
M3_THIRDS <- log2(0.01) + (log2(1e6) - log2(0.01)) * c(0, 1 / 3, 2 / 3, 1)  # fixed equal-width thirds of the generator's level range

m3Seed <- function(idx, stress, r) if (stress) 262200000 + 10000 * idx + r else 190000000 + 10000 * idx + r
m3Tag <- function(idx, r) sprintf("c%02d-r%d", idx, r)

m3Data <- function(cell, r) {
  if (!cell$stress) {
    d <- e2Data(cell$id, r)
    return(list(values = d$values, effects = d$effects, level = d$level, design = kit172Design(), itsd = TRUE, seed = d$seed))
  }
  seed <- m3Seed(idx = cell$idx, stress = TRUE, r = r)
  sim <- m3StressValues(seed, effect_share = if (identical(cell$kind, "null")) 0 else 0.1)
  # no `IS_` rows exist in the stress cells (as M2): the ITSD step has no standard to act on and is not run (ruling item)
  list(values = sim$values, effects = sim$effects, level = sim$level, design = sim$design, itsd = FALSE, seed = seed)
}

# ---- one journey and what it yields ---------------------------------------------------------------------------------------
# A journey that stops is returned as the condition (a data-set error row, never a denominator); a partial-match guard stop is
# re-thrown by the caller (an instrument row).
m3RunJourney <- function(root, object, omic, ruv, itsd) {
  tryCatch(withCallingHandlers(kit172Journey(root, object, omic, itsd = itsd, ruv = ruv), warning = m3GuardMuffle),
    error = function(e) e)
}

m3Pick <- function(x, assay = "LCMS_Pos") {
  if (is.null(x)) return(NULL)
  if (assay %in% names(x)) x[[assay]] else if (is.list(x)) NULL else x  # a list or a named vector holding the assay
}
m3Num <- function(x) { v <- m3Pick(x); if (is.null(v) || !length(v)) NA_real_ else suppressWarnings(as.numeric(v[[1L]])) }

# The endpoints of one journey (the P journey or a Manual one): the RUV outcome by the step's own fields, the DA calls and the
# calibration counts by abundance third.
m3Endpoints <- function(run, data, runs, n_samples) {
  # L3r: the step's per-assay fields are in `run$ruv$value` (`runWorkflowOmicsReferenceRuv` lifts only `optimisation` to the top level,
  # utils_workflow_omics_declarations.R:2624-2629). A missing field or assay is an instrument stop with a named message, never a silent NULL.
  step <- run[["ruv"]][["value"]]
  if (!is.list(step)) stop("M3 instrument: the RUV step result has no `value` list")
  for (fld in c("ruv_applied", "ruv_not_applied_reason", "best_k", "k_sweep", "k_cap")) {
    if (!fld %in% names(step)) stop(sprintf("M3 instrument: the RUV step's value has no field '%s'", fld))
    if (!(is.list(step[[fld]]) || is.atomic(step[[fld]])) || !"LCMS_Pos" %in% names(step[[fld]])) stop(sprintf("M3 instrument: the RUV step's field '%s' has no assay LCMS_Pos", fld))
  }
  opt <- run[["ruv"]][["optimisation"]][["LCMS_Pos"]]
  if (!is.list(opt)) stop("M3 instrument: the RUV step result has no optimisation for assay LCMS_Pos")
  applied <- m3Pick(step[["ruv_applied"]])
  reason <- m3Pick(step[["ruv_not_applied_reason"]])
  da <- kit172DaTable(run)
  if (!is.null(da[["assay"]])) da <- da[da[["assay"]] == "LCMS_Pos", , drop = FALSE]
  truth <- data$effects[da[["feature"]]]
  truth[is.na(truth)] <- 0
  q <- da[["fdr_qvalue"]]
  called <- da[["feature"]][!is.na(q) & q < 0.05]
  planted <- names(data$effects)[data$effects != 0]
  lvl <- data$level[da[["feature"]]]
  third <- cut(lvl, breaks = M3_THIRDS, include.lowest = TRUE, labels = FALSE)
  p <- da[["raw_pvalue"]]
  isnull <- truth == 0 & is.finite(p) & !is.na(third)
  nt <- vapply(1:3, function(i) sum(isnull & third == i), numeric(1L))
  np <- vapply(1:3, function(i) sum(isnull & third == i & p < 0.05), numeric(1L))
  k_used <- m3Num(step[["best_k"]])
  list(ruv = isTRUE(applied), applied_field = applied, reason = if (is.null(reason) || is.na(reason[[1L]])) NA_character_ else as.character(reason[[1L]]),
    success = isTRUE(opt[["success"]]), pct = if (isTRUE(applied)) as.numeric(opt[["best_percentage"]] %||% NA_real_) else NA_real_,
    k_sweep = m3Num(step[["k_sweep"]]), k_used = if (isTRUE(applied)) k_used else NA_real_, cap = m3Num(step[["k_cap"]]),
    resid_df = if (isTRUE(applied)) n_samples - 2 - k_used else NA_real_,
    calls = length(called), any_call = length(called) > 0L,
    fdp = if (length(called)) sum(!called %in% planted) / length(called) else 0,
    power = if (length(planted)) sum(called %in% planted) / length(planted) else NA_real_,
    nt = nt, np = np, da = da, opt = opt, trace = opt[["optimization_results"]])
}

# The committed D rows of KIT-190 (rule (d), the same generator and seeds; E9): the sweep's percentage and k, and the bounded k used.
m3CommittedD <- function() {
  root <- "docs/qa/workflow-execution/kit/kit-190"
  rd <- function(p) utils::read.delim(p, stringsAsFactors = FALSE, colClasses = "character")
  s2 <- do.call(rbind, lapply(Sys.glob(file.path(root, "bound-s2", "*", "journeys.tsv")), function(p) {
    x <- rd(p)
    x[x[["arm"]] %in% c("D", "D+B"), c("cell", "r", "arm", "percentage", "k", "calls")]
  }))
  ed <- do.call(rbind, lapply(Sys.glob(file.path(root, "ed", "*", "journeys.tsv")), function(p) {
    x <- rd(p)
    x <- x[x[["side"]] == "head" & x[["class"]] == "lipid" & x[["kind"]] == "e2", , drop = FALSE]
    data.frame(cell = x[["cell"]], r = x[["r"]], arm = "D", percentage = x[["percentage"]], k = x[["k"]], calls = x[["calls"]],
      stringsAsFactors = FALSE)
  }))
  rbind(s2, ed)
}

# ---- the row ------------------------------------------------------------------------------------------------------------
M3_COLUMNS <- c("cell", "r", "seed", "purpose", "omic", "status", "class", "error", "n_samples", "cap", "n_features",
  "t_n_ok", "t_n_r3", "t_n_r4", "t_n_r5", "t_n_r6", "t_n_err", "t_n_skipped_no_controls", "t_n_rejected",
  "p_ruv", "p_reason", "p_pct", "p_k_sweep", "p_k_used", "p_resid_df", "p_calls", "p_any_call", "p_fdp", "p_power",
  "p_nt1", "p_nt2", "p_nt3", "p_np1", "p_np2", "p_np3",
  "c_source", "c_ruv", "c_pct", "c_k_used", "c_resid_df", "c_calls", "c_any_call", "c_fdp", "c_power",
  "c_nt1", "c_nt2", "c_nt3", "c_np1", "c_np2", "c_np3", "c_differs", "change_class",
  "e1", "e7", "e9", "e9_compared", "d_pct", "d_k", "db_k", "d_calls", "e10", "e10_direct_no_ruv", "e11", "e11_run", "e8", "f2",
  "seconds_p", "seconds_c", "seconds_other", "seconds_total")
m3EmptyRow <- function(cell, r, seed, purpose, omic) {
  row <- stats::setNames(rep(list(NA), length(M3_COLUMNS)), M3_COLUMNS)
  row[c("cell", "r", "seed", "purpose", "omic", "status", "class", "error")] <- list(cell, r, seed, purpose, omic, "ok", "", "")
  row
}
m3Line <- function(row) paste(vapply(M3_COLUMNS, function(cn) {
  v <- row[[cn]]
  if (is.null(v) || !length(v) || is.na(v[[1L]])) return("NA")
  v <- v[[1L]]
  if (is.double(v)) sprintf("%.15g", v) else gsub("[\t\r\n]", " ", as.character(v))
}, ""), collapse = "\t")

# ---- one unit ------------------------------------------------------------------------------------------------------------
m3Unit <- function(task, omic = task$omic) {
  started <- proc.time()[[3L]]
  cell <- M3_CELLS[M3_CELLS$id == task$cell, ]
  row <- m3EmptyRow(cell = task$cell, r = task$r, seed = m3Seed(idx = cell$idx, stress = cell$stress, r = task$r), purpose = task$purpose, omic = omic)
  data <- m3Data(cell, task$r)
  runs <- data$design$Run
  n_samples <- length(runs)
  row$n_samples <- n_samples
  row$n_features <- sum(!startsWith(rownames(data$values), "IS_"))
  stopped <- function(e, step) {
    if (m3IsGuard(e)) stop(e)  # an instrument defect by definition: re-thrown to the unit's handler
    row$status <<- "error"; row$class <<- class(e)[[1L]]; row$error <<- paste0(step, ": ", gsub("[\t\r\n]", " ", conditionMessage(e)))
    row$seconds_total <<- proc.time()[[3L]] - started
    row
  }
  object_for <- function() if (cell$stress) m3Object(data$values, omic, data$design) else kit172Object(data$values, omic)
  state_matrix <- function(run) {
    m <- if (cell$stress) m3StateMatrix(run, paste0(run$prefix, "_normalized"), runs) else kit172StateMatrix(run, paste0(run$prefix, "_normalized"))
    m[!is.finite(m)] <- NA_real_
    m
  }
  # ---- the P journey: the shipped Automatic step -----------------------------------------------------------------------
  planted_fault <- task$fault
  root_p <- file.path(tempdir(), sprintf("k262m3-%s-%s-P", task$tag, omic))
  on.exit(unlink(root_p, recursive = TRUE), add = TRUE)
  t0 <- proc.time()[[3L]]
  if (identical(planted_fault, "error")) {
    run_p <- simpleError("K262 planted journey stop (L3's failure check)")
  } else if (identical(planted_fault, "guard")) {
    run_p <- tryCatch(withCallingHandlers({ invisible(data.frame(abcdef = 1)$abc); simpleError("the planted partial match was not signalled") },
      warning = m3GuardMuffle), error = function(e) e)
  } else {
    run_p <- m3RunJourney(root = root_p, object = object_for(), omic = omic, ruv = kit172LegacyRuv("automatic"), itsd = data$itsd)
  }
  row$seconds_p <- proc.time()[[3L]] - t0
  if (inherits(run_p, "error")) return(stopped(run_p, "P journey"))
  on.exit(closeDIAStateAdapter(run_p$owner), add = TRUE, after = FALSE)
  if (identical(planted_fault, "unit")) invisible(data.frame(abcdef = 1)$abc)  # a partial match in unit code outside any journey
  if (identical(planted_fault, "inner_trycatch")) tryCatch(invisible(data.frame(abcdef = 1)$abc), error = function(e) NULL)  # a partial match inside an inner tryCatch(error =): the guard's stop is signalled past it (L3r2: still an instrument row, class partial_match_guard)
  p <- m3Endpoints(run = run_p, data = data, runs = runs, n_samples = n_samples)
  trace <- p$trace
  if (!is.data.frame(trace) || !nrow(trace)) stop("the P journey returned no optimisation trace")
  cap <- MultiScholaR:::calculateAdaptiveMaxK(trace$sample_size[[1L]])
  row$cap <- cap
  status <- trace$status
  row$t_n_ok <- sum(status == "ok"); row$t_n_r3 <- sum(status == "skipped_insufficient_controls")
  row$t_n_r4 <- sum(status == "invalid_cancor_plot"); row$t_n_r5 <- sum(status == K262_R5); row$t_n_r6 <- sum(status == "error_scoring")
  row$t_n_err <- sum(status %in% c("error_neg_ctrl_selection", "error_cancor"))
  row$t_n_skipped_no_controls <- sum(status == "skipped_no_controls_selected")
  row$t_n_rejected <- sum(status != "ok")  # every status but ok is a rejected percentage, `skipped_no_controls_selected` included
  # ---- the rules from the trace (no hook: C1 needs only R5's reason and the trace; O-1 (ii)) ----------------------------------
  complete <- rep(NA_integer_, nrow(trace))
  parsed <- lapply(trace$error_reason, k262ParseReason)
  is_r5 <- status == K262_R5
  complete[is_r5] <- vapply(parsed[is_r5], function(x) x$c, integer(1L))
  ru <- k262Rules(trace = trace, curves = NULL, complete = complete, cap = cap)
  wp <- ru$winner$P; w1 <- ru$winner$C1
  # E1: P's usable set, composites, winner and k used against the shipped trace and the step's own outcome
  shipped_usable <- status == "ok"
  same_winner <- if (isTRUE(p$success)) identical(wp$percentage, as.numeric(p$opt[["best_percentage"]])) &&
    identical(trace$best_k[[wp$row]], as.integer(p$opt[["best_k"]])) else is.na(wp$row)
  row$e1 <- identical(ru$rows$p_usable, shipped_usable) && same_winner &&
    isTRUE(all.equal(ru$rows$composite_recomputed[shipped_usable], trace$composite_score[shipped_usable], tolerance = 0)) &&
    all(vapply(parsed[is_r5], function(x) x$ok, logical(1L))) && identical(p$ruv, isTRUE(p$success)) && isTRUE(p$cap == cap || is.na(p$cap)) &&
    (!p$ruv || isTRUE(p$k_used == min(p$k_sweep, cap)))
  # E7 (reduced under O-1 (ii)): P within C1 row by row, C1's k used within the cap, and RUV under P implies RUV under C1
  row$e7 <- all(!ru$rows$p_usable | ru$rows$c1_usable) && (is.na(w1$row) || w1$k_used <= cap) && (is.na(wp$row) || !is.na(w1$row))
  # P's row
  row$p_ruv <- p$ruv; row$p_reason <- p$reason; row$p_pct <- p$pct; row$p_k_sweep <- p$k_sweep; row$p_k_used <- p$k_used
  row$p_resid_df <- p$resid_df; row$p_calls <- p$calls; row$p_any_call <- p$any_call; row$p_fdp <- p$fdp; row$p_power <- p$power
  for (i in 1:3) { row[[paste0("p_nt", i)]] <- p$nt[[i]]; row[[paste0("p_np", i)]] <- p$np[[i]] }
  # ---- E9: P's percentage and sweep k against KIT-190's committed D row, P's k used against the D+B row (eight cells only) ----
  if (!cell$stress) {
    d <- M3_D[M3_D$cell == task$cell & M3_D$r == as.character(task$r), , drop = FALSE]
    dd <- d[d$arm == "D", , drop = FALSE]; dbb <- d[d$arm == "D+B", , drop = FALSE]
    if (nrow(dd) == 1L && nrow(dbb) == 1L) {
      row$d_pct <- suppressWarnings(as.numeric(dd$percentage)); row$d_k <- suppressWarnings(as.numeric(dd$k))
      row$db_k <- suppressWarnings(as.numeric(dbb$k)); row$d_calls <- suppressWarnings(as.numeric(dd$calls))
      row$e9_compared <- TRUE
      same_num <- function(a, b) (is.na(a) && is.na(b)) || isTRUE(a == b)
      row$e9 <- same_num(p$pct, row$d_pct) && same_num(p$k_sweep, row$d_k) && same_num(p$k_used, row$db_k)
    } else row$e9_compared <- FALSE
  }
  if (identical(task$purpose, "e9")) {  # E9's P-only unit
    row$seconds_total <- proc.time()[[3L]] - started
    return(row)
  }
  # ---- E10 (stress cells): the direct M2-form sweep on the same data, before between-sample normalisation, against the journey's ----
  if (cell$stress) {
    t1 <- proc.time()[[3L]]
    direct <- m2Sweep(m2Prepare(m2Object(data$values, data$design, omic), omic))[[1L]]
    row$e10_direct_no_ruv <- !isTRUE(direct[["success"]])
    row$e10 <- TRUE  # the comparison itself is reported by the scorer (McNemar counts); the check here is that the sweep ran
    row$seconds_other <- proc.time()[[3L]] - t1
  }
  # ---- the C1 arm -----------------------------------------------------------------------------------------------------
  differs <- if (is.na(w1$row)) FALSE else !(p$ruv && isTRUE(w1$percentage == p$pct))
  row$c_differs <- differs
  if (is.na(w1$row)) {
    row$c_source <- "none"   # no usable percentage under C1 (hence none under P, E7): the uncorrected assay, P's row
  } else if (!differs) row$c_source <- "reused" else row$c_source <- "manual"
  copy_p <- function() {
    row$c_ruv <<- p$ruv; row$c_pct <<- p$pct; row$c_k_used <<- p$k_used; row$c_resid_df <<- p$resid_df; row$c_calls <<- p$calls
    row$c_any_call <<- p$any_call; row$c_fdp <<- p$fdp; row$c_power <<- p$power
    for (i in 1:3) { row[[paste0("c_nt", i)]] <<- p$nt[[i]]; row[[paste0("c_np", i)]] <<- p$np[[i]] }
  }
  manual_at <- function(pct, k) {
    ruv <- kit172LegacyRuv("manual")
    ruv$manual_percentage <- pct
    ruv$manual_k <- as.integer(k)
    root_c <- file.path(tempdir(), sprintf("k262m3-%s-%s-C-%s", task$tag, omic, gsub("[^0-9]", "", format(pct))))
    on.exit(unlink(root_c, recursive = TRUE), add = TRUE)
    run <- m3RunJourney(root = root_c, object = object_for(), omic = omic, ruv = ruv, itsd = data$itsd)
    if (inherits(run, "error")) return(list(error = run))
    on.exit(closeDIAStateAdapter(run$owner), add = TRUE, after = FALSE)
    list(endpoints = m3Endpoints(run = run, data = data, runs = runs, n_samples = n_samples), input = state_matrix(run))
  }
  if (identical(row$c_source, "manual")) {
    t2 <- proc.time()[[3L]]
    m <- manual_at(pct = w1$percentage, k = w1$k_used)
    row$seconds_c <- proc.time()[[3L]] - t2
    if (!is.null(m$error)) return(stopped(m$error, "C1 Manual journey"))
    c1 <- m$endpoints
    row$c_ruv <- c1$ruv; row$c_pct <- c1$pct; row$c_k_used <- c1$k_used; row$c_resid_df <- c1$resid_df; row$c_calls <- c1$calls
    row$c_any_call <- c1$any_call; row$c_fdp <- c1$fdp; row$c_power <- c1$power
    for (i in 1:3) { row[[paste0("c_nt", i)]] <- c1$nt[[i]]; row[[paste0("c_np", i)]] <- c1$np[[i]] }
    # F2: the Manual journey's percentage, k and controls are the declared ones and rule (d)'s oracle at that percentage
    orc <- k190Oracle(m$input, data$design$group, w1$percentage)
    controls <- names(c1$opt[["control_genes_index"]])[c1$opt[["control_genes_index"]]]
    row$f2 <- isTRUE(c1$ruv) && isTRUE(c1$pct == w1$percentage) && isTRUE(c1$k_used == w1$k_used) && length(controls) > 0L &&
      setequal(controls, k190Candidate(orc, "d"))
  } else {
    copy_p()
    # E11: a Manual journey at P's own percentage and k reproduces P's DA table (a declared 5 % sample of the re-used rows; every
    # re-used row in L3 under K262_M3_FORCE_ARMS)
    do_e11 <- identical(row$c_source, "reused") && (identical(Sys.getenv("K262_M3_FORCE_ARMS"), "1") || task$r %% e11_mod == 7L)
    if (do_e11) {
      t3 <- proc.time()[[3L]]
      m <- manual_at(pct = p$pct, k = p$k_used)
      row$seconds_other <- sum(row$seconds_other, proc.time()[[3L]] - t3, na.rm = TRUE)
      row$e11_run <- TRUE
      if (!is.null(m$error)) {
        row$e11 <- FALSE
      } else {
        cols <- c("feature", "logFC", "raw_pvalue", "fdr_qvalue")
        a <- p$da[order(p$da[["feature"]]), cols, drop = FALSE]; b <- m$endpoints$da[order(m$endpoints$da[["feature"]]), cols, drop = FALSE]
        rownames(a) <- NULL; rownames(b) <- NULL
        row$e11 <- isTRUE(m$endpoints$ruv) && isTRUE(m$endpoints$pct == p$pct) && isTRUE(m$endpoints$k_used == p$k_used) && identical(a, b)
      }
    } else row$e11_run <- FALSE
  }
  # the change class: rescued (P without RUV, C1 with), re-ranked (both with RUV, another percentage), same otherwise
  row$change_class <- if (!differs) "same" else if (!p$ruv) "rescued" else "reranked"
  # ---- E8: the metabolomics identity (r <= 5): the same unit through the metabolomics class gives the same row ------------------
  if (identical(omic, "lipid") && use_e8 && task$r <= e8_max && !identical(task$purpose, "e9")) {
    t4 <- proc.time()[[3L]]
    task_m <- task; task_m$fault <- NULL
    other <- m3Unit(task_m, omic = "metab")
    skip <- c("omic", "seconds_p", "seconds_c", "seconds_other", "seconds_total", "e8")
    keep <- setdiff(M3_COLUMNS, skip)
    row$e8 <- identical(vapply(keep, function(cn) as.character(other[[cn]]), ""), vapply(keep, function(cn) as.character(row[[cn]]), ""))
    row$seconds_other <- sum(row$seconds_other, proc.time()[[3L]] - t4, na.rm = TRUE)
  }
  row$seconds_total <- proc.time()[[3L]] - started
  row
}

# A unit as the child runs it: every outcome becomes one appended row; a stop from the guard or from unit code is an instrument row.
m3UnitSafe <- function(task) {
  started <- proc.time()[[3L]]
  cell <- M3_CELLS[M3_CELLS$id == task$cell, ]
  empty <- m3EmptyRow(cell = task$cell, r = task$r, seed = m3Seed(idx = cell$idx, stress = cell$stress, r = task$r), purpose = task$purpose, omic = task$omic)
  if (identical(task$fault, "kill")) tools::pskill(Sys.getpid(), tools::SIGKILL)  # the worker dies with no row and no trace
  row <- tryCatch(withCallingHandlers(m3Unit(task), warning = m3GuardMuffle), error = function(e) {
    trace <- paste(c(conditionMessage(e), getOption("k262.guard_trace", "")), collapse = "\n")
    try(writeLines(trace, file.path(out_dir, "traces", sprintf("trace-%s.txt", task$tag))), silent = TRUE)
    empty$status <- "instrument"
    empty$class <- if (m3IsGuard(e)) "partial_match_guard" else "script"
    empty$error <- gsub("[\t\r\n]", " ", conditionMessage(e))
    empty$seconds_total <- proc.time()[[3L]] - started
    empty
  })
  if (identical(row$status, "error")) try(writeLines(c(row$class, row$error), file.path(out_dir, "traces", sprintf("trace-%s.txt", task$tag))), silent = TRUE)
  cat(m3Line(row), "\n", sep = "", file = summary_file, append = TRUE)  # the unit's own record, before it returns
  cat(sprintf("unit %s: %s (%.0f s)\n", task$tag, row$status, as.numeric(row$seconds_total)))  # progress: the tag and the status only
  row$status
}

# ---- the plan --------------------------------------------------------------------------------------------------------------
parseSet <- function(text, default) {
  if (!nzchar(text)) return(default)
  unlist(lapply(strsplit(text, ",", fixed = TRUE)[[1L]], function(piece) {
    bounds <- as.integer(strsplit(piece, ":", fixed = TRUE)[[1L]])
    if (length(bounds) == 2L) seq.int(bounds[[1L]], bounds[[2L]]) else bounds
  }))
}
lifecycle <- identical(Sys.getenv("K262_LIFECYCLE"), "1")
stage1 <- identical(Sys.getenv("K262_M3_STAGE1"), "1")
fault_on <- identical(Sys.getenv("K262_M3_FAULT"), "1")
stopifnot(!fault_on || lifecycle)
workers <- as.integer(Sys.getenv("K262_M3_WORKERS", "4"))
use_e8 <- !identical(Sys.getenv("K262_M3_E8", "1"), "0")
e8_max <- as.integer(Sys.getenv("K262_M3_E8_MAX", "5"))
e11_mod <- as.integer(Sys.getenv("K262_M3_E11MOD", "20"))
cells_sel <- strsplit(Sys.getenv("K262_M3_CELLS", if (stage1) "S-N,S-E" else paste(names(M3_FINAL_R), collapse = ",")), ",", fixed = TRUE)[[1L]]
stopifnot(all(cells_sel %in% M3_CELLS$id), !stage1 || setequal(cells_sel, c("S-N", "S-E")))
rmap <- utils::modifyList(as.list(M3_FINAL_R), if (nzchar(Sys.getenv("K262_M3_RMAP"))) {
  kv <- strsplit(strsplit(Sys.getenv("K262_M3_RMAP"), ",", fixed = TRUE)[[1L]], "=", fixed = TRUE)
  stats::setNames(lapply(kv, function(x) as.integer(x[[2L]])), vapply(kv, `[[`, "", 1L))
} else list())
reps_env <- Sys.getenv("K262_M3_REPS", if (stage1) "1:50" else "")
reps_of <- function(cell_id) {
  if (nzchar(reps_env)) return(parseSet(reps_env, integer()))
  if (is.null(rmap[[cell_id]])) stop("K262 M3: no replicate count for cell ", cell_id, " (a stress cell takes its R from K262_M3_RMAP or K262_M3_REPS)")
  seq_len(rmap[[cell_id]])
}
tasks <- do.call(rbind, lapply(cells_sel, function(id) data.frame(cell = id, r = reps_of(id), purpose = "main", stringsAsFactors = FALSE)))
e9_reps <- parseSet(Sys.getenv("K262_M3_E9REPS", ""), integer())
if (length(e9_reps)) tasks <- rbind(tasks, do.call(rbind, lapply(intersect(cells_sel, names(M3_FINAL_R)), function(id)
  data.frame(cell = id, r = e9_reps, purpose = "e9", stringsAsFactors = FALSE))))
tasks$idx <- M3_CELLS$idx[match(tasks$cell, M3_CELLS$id)]
tasks$tag <- ifelse(tasks$purpose == "main", m3Tag(tasks$idx, tasks$r), paste0(m3Tag(tasks$idx, tasks$r), "-e9"))
tasks <- tasks[order(tasks$purpose, tasks$r, tasks$idx), , drop = FALSE]
stopifnot(!anyDuplicated(tasks$tag))
# CONTINUATION (counted mode): exactly the units named, once. A tag outside this job's plan stops; so does any use in L3.
only <- Sys.getenv("K262_M3_ONLY_UNITS")
if (nzchar(only)) {
  if (lifecycle) stop("K262_M3_ONLY_UNITS is a counted-mode setting (the continuation); it is refused in the lifecycle")
  keep <- trimws(strsplit(only, ",", fixed = TRUE)[[1L]])
  if (!length(keep) || !all(nzchar(keep)) || !all(keep %in% tasks$tag)) stop("K262 M3: K262_M3_ONLY_UNITS names a unit outside this job's plan: ",
    paste(setdiff(keep, tasks$tag), collapse = ", "))
  tasks <- tasks[tasks$tag %in% keep, , drop = FALSE]
}
tasks$omic <- "lipid"
tasks$fault <- NA_character_
if (fault_on) tasks$fault[seq_len(min(5L, nrow(tasks)))] <- c("error", "guard", "unit", "inner_trycatch", "kill")[seq_len(min(5L, nrow(tasks)))]
rownames(tasks) <- NULL
M3_D <- if (any(tasks$purpose == "e9" | (tasks$cell %in% names(M3_FINAL_R)))) m3CommittedD() else NULL

# The lifecycle and stage 1 keep their row files under tempdir(): no share, rate, FDP, power or k distribution leaves them (VAL-11).
out_dir <- if (lifecycle || stage1) file.path(tempdir(), "m3-private") else out
dir.create(file.path(out_dir, "traces"), recursive = TRUE, showWarnings = FALSE)
summary_file <- file.path(out_dir, "summary.tsv")
if (file.exists(summary_file)) stop("K262 M3: summary.tsv already exists in the output directory; a continuation uses its own directory")
writeLines(paste(M3_COLUMNS, collapse = "\t"), summary_file)

sha <- function(f) tryCatch(strsplit(system2("sha256sum", shQuote(f), stdout = TRUE), " ")[[1L]][[1L]], error = function(e) NA_character_)
script_dir <- "docs/qa/workflow-execution/kit/kit-262"
provenance <- function(finished) {
  versions <- vapply(c("MultiScholaR", "limma", "RUVIIIC", "ruv", "qvalue", "preprocessCore"), function(pk)
    tryCatch(as.character(utils::packageVersion(pk)), error = function(e) NA_character_), "")
  list(ticket = "KIT-262 M3", commit = tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE), error = function(e) NA_character_),
    scripts = stats::setNames(as.list(vapply(file.path(script_dir, c("m3-k262.R", "rules-k262.R", "score-m3-k262.R", "reestimate-m3-k262.R")), sha, "")),
      c("m3-k262.R", "rules-k262.R", "score-m3-k262.R", "reestimate-m3-k262.R")),
    r_version = R.version.string, packages = as.list(versions), cells = cells_sel, units_planned = nrow(tasks), main_units = sum(tasks$purpose == "main"),
    stress_floor_order = M3_STRESS_FLOOR_ORDER, seed_rule = "eight cells 190000000 + 10000 x idx + r (VAL-235's); stress 262200000 + 10000 x idx + r (idx 9, 10)",
    lifecycle = lifecycle, stage1 = stage1, continuation = nzchar(only), finished = finished, time_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
}
if (!lifecycle && !stage1) jsonlite::write_json(provenance(FALSE), file.path(out, "provenance.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)

cat(sprintf("M3: %d units planned (%d main) in %d cells; workers %d; lifecycle %s; stage1 %s; continuation %s\n", nrow(tasks), sum(tasks$purpose == "main"),
  length(unique(tasks$cell)), workers, lifecycle, stage1, nzchar(only)))
results <- parallel::mclapply(seq_len(nrow(tasks)), function(i) tryCatch(m3UnitSafe(as.list(tasks[i, ])), error = function(e) paste("instrument:", conditionMessage(e))),
  mc.cores = workers, mc.preschedule = FALSE)

# A unit with no row (a worker that did not return) gets an instrument row from the parent: with its trace if one exists, else "no result and no trace".
have <- utils::read.delim(summary_file, stringsAsFactors = FALSE, colClasses = "character")
for (i in seq_len(nrow(tasks))) {
  tag_rows <- have[have$cell == tasks$cell[[i]] & have$r == as.character(tasks$r[[i]]) & have$purpose == tasks$purpose[[i]], , drop = FALSE]
  if (nrow(tag_rows) >= 1L) next
  cell <- M3_CELLS[M3_CELLS$id == tasks$cell[[i]], ]
  row <- m3EmptyRow(cell = tasks$cell[[i]], r = tasks$r[[i]], seed = m3Seed(idx = cell$idx, stress = cell$stress, r = tasks$r[[i]]), purpose = tasks$purpose[[i]], omic = "lipid")
  trace_file <- file.path(out_dir, "traces", sprintf("trace-%s.txt", tasks$tag[[i]]))
  row$status <- "instrument"; row$class <- "no_row"
  row$error <- if (file.exists(trace_file)) paste("no row; trace:", gsub("[\t\r\n]", " ", paste(readLines(trace_file, warn = FALSE), collapse = " | "))) else
    "no result and no trace (the worker did not return)"
  cat(m3Line(row), "\n", sep = "", file = summary_file, append = TRUE)
}
final <- utils::read.delim(summary_file, stringsAsFactors = FALSE, colClasses = "character")
for (cn in c("seconds_p", "seconds_c", "seconds_other", "seconds_total")) final[[cn]] <- suppressWarnings(as.numeric(final[[cn]]))

# ---- checks (counts only; no estimate) -----------------------------------------------------------------------------------------
ok_rows <- final[final$status == "ok", , drop = FALSE]
flagcount <- function(col, rows = ok_rows) {
  v <- rows[[col]]
  data.frame(check = toupper(col), units = sum(v %in% c("TRUE", "FALSE")), passed = sum(v %in% "TRUE"), failed = sum(v %in% "FALSE"), stringsAsFactors = FALSE)
}
checks <- do.call(rbind, lapply(c("e1", "e7", "e9", "e10", "e11", "e8"), flagcount))
# S3 (audit): F2 rests on KIT-190's oracle, validated only at n = 12 with at most 3 of 6 values missing per group; in the stress cells it is OBSERVE
# (reported as F2_STRESS, not a hard check of L3: a difference goes to the coordinator). F2 on the eight cells is hard.
f2_stress <- ok_rows$cell %in% M3_CELLS$id[M3_CELLS$stress]
checks <- rbind(checks, flagcount("f2", ok_rows[!f2_stress, , drop = FALSE]), transform(flagcount("f2", ok_rows[f2_stress, , drop = FALSE]), check = "F2_STRESS"))
expected_n <- nrow(tasks)
checks <- rbind(checks,
  data.frame(check = "COMPLETION", units = expected_n, passed = nrow(final), failed = abs(expected_n - nrow(final)), stringsAsFactors = FALSE),
  data.frame(check = "INSTRUMENT_ROWS", units = nrow(final), passed = sum(final$status != "instrument"), failed = sum(final$status == "instrument"), stringsAsFactors = FALSE),
  data.frame(check = "ERROR_ROWS_REPORTED", units = nrow(final), passed = sum(final$status != "error"), failed = sum(final$status == "error"), stringsAsFactors = FALSE))
seed_ok <- vapply(seq_len(nrow(final)), function(i) {
  cell <- M3_CELLS[M3_CELLS$id == final$cell[[i]], ]
  isTRUE(as.numeric(final$seed[[i]]) == m3Seed(idx = cell$idx, stress = cell$stress, r = as.numeric(final$r[[i]])))
}, logical(1L))
checks <- rbind(checks, data.frame(check = "SEEDS", units = nrow(final), passed = sum(seed_ok), failed = sum(!seed_ok), stringsAsFactors = FALSE))
checks <- rbind(checks, data.frame(check = "FLOOR_ORDER_STRESS", units = 1L, passed = as.integer(identical(M3_STRESS_FLOOR_ORDER, "after")),
  failed = as.integer(!identical(M3_STRESS_FLOOR_ORDER, "after")), stringsAsFactors = FALSE))
# the disjointness of the stress seed range from every earlier range (E9's last clause): M2 up to 262119895, the re-measurement from 263101101
stress_seeds <- as.numeric(final$seed[final$cell %in% c("S-N", "S-E")])
disjoint <- !length(stress_seeds) || (min(stress_seeds) > 262119895 && max(stress_seeds) < 263100000)
checks <- rbind(checks, data.frame(check = "SEED_RANGE", units = length(stress_seeds), passed = if (disjoint) length(stress_seeds) else 0L,
  failed = if (disjoint) 0L else length(stress_seeds), stringsAsFactors = FALSE))
if (fault_on) {
  by_cell <- function(i) final[final$cell == tasks$cell[[i]] & final$r == as.character(tasks$r[[i]]) & final$purpose == tasks$purpose[[i]], , drop = FALSE]
  task_tag_of <- function(i) tasks$tag[[i]]
  fr <- lapply(1:5, by_cell)
  want <- list(c("error", ""), c("instrument", "partial_match_guard"), c("instrument", "partial_match_guard"), c("instrument", "partial_match_guard"), c("instrument", "no_row"))
  hit <- vapply(1:5, function(i) nrow(fr[[i]]) == 1L && identical(fr[[i]]$status, want[[i]][[1L]]) && (!nzchar(want[[i]][[2L]]) || identical(fr[[i]]$class, want[[i]][[2L]])), logical(1L))
  others_ok <- sum(final$status == "ok") == nrow(final) - 5L
  # L3r: per planted unit its id, the expected and the observed status and class (aggregate: no data value), so a failing unit is seen in the log
  for (i in 1:5) cat(sprintf("FAULT unit %d (%s, r %s): expected %s / %s, observed %s / %s -> %s\n", i, task_tag_of(i), tasks$r[[i]], want[[i]][[1L]], if (nzchar(want[[i]][[2L]])) want[[i]][[2L]] else "(any)",
    if (nrow(fr[[i]]) == 1L) fr[[i]]$status else "no row", if (nrow(fr[[i]]) == 1L) fr[[i]]$class else "", if (hit[[i]]) "ok" else "MISMATCH"))
  checks <- rbind(checks, data.frame(check = "FAULT_ROWS", units = 6L, passed = sum(hit) + as.integer(others_ok), failed = 6L - sum(hit) - as.integer(others_ok), stringsAsFactors = FALSE))
}
# ---- the lifecycle's scorer rehearsal (blind: counts only): the blind scorer twice (E6, byte-identical files, every estimand NA) and the full
# scorer once on the rows just written (the format the simulation writes is the format the scorer reads); nothing it computes leaves the run ----
if (lifecycle && !fault_on) {
  source(file.path(script_dir, "score-m3-k262.R"))  # sourced: defines the functions; its command line runs nothing when the first argument is no mode
  lc_main <- final[final$purpose == "main" & final$status == "ok", , drop = FALSE]
  lc_reps <- sort(unique(as.integer(lc_main$r)))
  lc_plan <- data.frame(cell = unique(lc_main$cell), R = length(lc_reps), stringsAsFactors = FALSE)
  ra <- tryCatch(blindM3(dirs = out_dir, plan = lc_plan, out_dir = file.path(tempdir(), "blind-a"), reps = lc_reps), error = function(e) e)
  rb <- tryCatch(blindM3(dirs = out_dir, plan = lc_plan, out_dir = file.path(tempdir(), "blind-b"), reps = lc_reps), error = function(e) e)
  same_files <- !inherits(ra, "error") && !inherits(rb, "error") && identical(list.files(file.path(tempdir(), "blind-a")), list.files(file.path(tempdir(), "blind-b"))) &&
    all(unname(tools::md5sum(file.path(tempdir(), "blind-a", list.files(file.path(tempdir(), "blind-a"))))) == unname(tools::md5sum(file.path(tempdir(), "blind-b", list.files(file.path(tempdir(), "blind-b"))))))
  na_ok <- !inherits(ra, "error") && all(is.na(ra$share_rescued)) && all(is.na(ra$share_reranked)) && all(is.na(ra$calibration)) && all(is.na(ra$power)) && all(is.na(ra$verdict))
  full_ok <- tryCatch({ fullM3(dirs = out_dir, plan = lc_plan, a_gain = TRUE, out_dir = file.path(tempdir(), "full-lc"), B = 200L, reps = lc_reps); TRUE }, error = function(e) FALSE)
  checks <- rbind(checks,
    data.frame(check = "E6_BLIND_SCORER", units = 1L, passed = as.integer(same_files), failed = as.integer(!same_files), stringsAsFactors = FALSE),
    data.frame(check = "BLIND_ESTIMANDS_NA", units = 1L, passed = as.integer(na_ok), failed = as.integer(!na_ok), stringsAsFactors = FALSE),
    data.frame(check = "FULL_SCORER_ON_ROWS", units = 1L, passed = as.integer(full_ok), failed = as.integer(!full_ok), stringsAsFactors = FALSE))
}
checks$status <- ifelse(checks$failed == 0L, "pass", "FAIL")
timing <- if (nrow(ok_rows)) stats::aggregate(cbind(seconds_total, seconds_p) ~ cell + purpose, data = final[final$status == "ok", ], FUN = mean) else data.frame()
# ---- what leaves the run -------------------------------------------------------------------------------------------------------
if (stage1) {
  # blind stage 1: only the re-estimation's own table and fixed console lines (VAL-185's form); the rows stay under tempdir()
  capped <- do.call(rbind, lapply(c("P", "C1"), function(arm) {
    pre <- if (arm == "P") "p_" else "c_"
    d <- ok_rows[ok_rows$cell %in% c("S-N", "S-E"), , drop = FALSE]
    data.frame(cell = d$cell, arm = arm, r = as.integer(d$r), fdp = suppressWarnings(as.numeric(d[[paste0(pre, "fdp")]])),
      any_call = as.logical(d[[paste0(pre, "any_call")]]), stringsAsFactors = FALSE)
  }))
  rows_rds <- file.path(tempdir(), "m3-capped.rds"); saveRDS(capped, rows_rds)
  code <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(file.path(script_dir, "reestimate-m3-k262.R")), "estimate", shQuote(rows_rds), shQuote(out)),
    stdout = file.path(tempdir(), "reestimate.log"), stderr = file.path(tempdir(), "reestimate.log"))
  for (cl in unique(final$cell)) cat(sprintf("units %s: %d\n", cl, sum(final$cell == cl)))
  cat(sprintf("instrument rows: %d; error rows: %d; seconds: %.0f\n", sum(final$status == "instrument"), sum(final$status == "error"), sum(final$seconds_total, na.rm = TRUE)))
  if (lifecycle) {
    # the lifecycle (L3) proves the stage-1 mechanics only: the estimate's schema and pair counts are shown, its values never; the file is removed
    est_file <- file.path(out, "reestimate.tsv")
    est <- NULL
    if (code == 0L && file.exists(est_file)) {
      est <- utils::read.delim(est_file, stringsAsFactors = FALSE)
      print(est[c("cell", "endpoint", "n_pairs")])
      unlink(est_file)
    } else cat("reestimate: failed\n")
    checks$status <- NULL
    checks <- rbind(checks, data.frame(check = "STAGE1_ESTIMATE_SCHEMA", units = 1L, passed = as.integer(code == 0L && identical(names(est), c("cell", "endpoint", "n_pairs", "value")) && all(est$n_pairs > 1L)),
      failed = as.integer(!(code == 0L && identical(names(est), c("cell", "endpoint", "n_pairs", "value")) && all(est$n_pairs > 1L))), stringsAsFactors = FALSE))
  } else {
    if (code == 0L && file.exists(file.path(out, "reestimate.tsv"))) cat(readLines(file.path(out, "reestimate.tsv")), sep = "\n") else cat("reestimate: failed\n")
    jsonlite::write_json(c(provenance(TRUE), list(seconds_median_per_unit = stats::median(final$seconds_total, na.rm = TRUE))), file.path(out, "provenance.json"),
      auto_unbox = TRUE, pretty = TRUE, digits = NA)
    quit(status = if (code != 0L || any(final$status == "instrument")) 1L else 0L)
  }
}
checks$status <- ifelse(checks$failed == 0L, "pass", "FAIL")
utils::write.table(checks, file.path(out, "checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(timing, file.path(out, "timing.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (!lifecycle) {
  jsonlite::write_json(provenance(TRUE), file.path(out, "provenance.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
}
hashes <- data.frame(file = "summary.tsv", md5 = unname(tools::md5sum(summary_file)))
utils::write.table(hashes, file.path(out, "hashes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(checks)
print(timing)
cat(sprintf("units: %d; ok: %d; error rows: %d; instrument rows: %d\n", nrow(final), sum(final$status == "ok"), sum(final$status == "error"), sum(final$status == "instrument")))
# the exit status reads the hard checks only; E8, E10, E11 and the error-row count are reported in checks.tsv (a lifecycle finding or a scorer input)
hard <- c("E1", "E7", "E9", "F2", "COMPLETION", "INSTRUMENT_ROWS", "SEEDS", "SEED_RANGE", "FLOOR_ORDER_STRESS", "FAULT_ROWS", "E6_BLIND_SCORER", "BLIND_ESTIMANDS_NA",
  "FULL_SCORER_ON_ROWS", "STAGE1_ESTIMATE_SCHEMA")
quit(status = if (any(checks$status[checks$check %in% hard] == "FAIL")) 1L else 0L)
