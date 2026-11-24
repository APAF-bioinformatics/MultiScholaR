# KIT-190 E (ticket section 5.3; acceptance criterion A5): the effect of choosing lipidomics and metabolomics legacy
# RUV control features by the ANOVA instead of by row order, base against head; and E3, the run time of the repaired
# selection at realistic sizes. Evidence code under docs/qa: no package code. Descriptive: no verdict.
#
# Usage, from the package root of the bound commit on a run-tests runner (full git history):
#   Rscript docs/qa/workflow-execution/kit/kit-190/effect.R <e | e3> <base commit | self> <out dir>
#   e             E1 (KIT-172 T11 and T1b fixtures) and E2 (the seeded simulation, cells N-r, N-s, E-r, E-s, E-p, E-g).
#   e3            E3: legacy Automatic and Manual (5 %) on the head tree at n features per assay, two assays, and one
#                 base Automatic journey at the largest n; seconds of the RUV step and of the journey; projection.
#   <base commit> 40 digits: checked out as a git worktree and loaded as "base"; this tree is "head". `self`: this tree
#                 is loaded as both sides (the lifecycle).
# Environment:
#   K190_LIFECYCLE=1  R = 2 per cell, E1 on the two T11 data sets, the oracle check (cell N-r again with the grouping
#                     `factor1`, counted only there), falsifiers counted and never stopping; under `e3`, n = 520 and one
#                     head journey only.
#   K190_FAULT=1      cell N-r alone; the "head" side drops one control in its first lipid journey; the stop rule on.
#   K190_CELLS        comma list among E1,N-r,N-s,E-r,E-s,E-p,E-g (default all), to split the work by cell.
#   K190_WORKERS      forked workers per side (default 2; the two sides run at once on the runner's 4 cores).
#   K190_RULE=d       the head carries rule (d) of section 2.7a (VAL-110, VAL-126): F2 compares the head's controls
#                     with that rule's selection instead of the first repair's.
#   K190_E3_N         comma list of features per assay for `e3` (default 5000,20000).
# Reported only (VAL-100, for the owner's D2): per journey, what candidate rules (b) rank by p with a feature-id
# tie-break, (c) one-group features not candidates, (d) both, would select at the chosen percentage (rule_* columns).
# Each side runs in a fresh Rscript (`--side`) with pkgload::load_all(<tree>, export_all = TRUE) and the test helpers of
# the head tree sourced read-only (tests/testthat/helper-*.R); journeys as kit172Journey (helper-kit172.R:110-130),
# RUV kit172LegacyRuv("automatic") (:149-153). parallel::mclapply, mc.preschedule = FALSE, one journey per task.
# Written to <out dir>: cells.tsv, fixtures.tsv, falsifiers.tsv, journeys.tsv, provenance.json (`e`); e3/times.tsv
# (`e3`). Outside lifecycle mode, and in the fault run, the script stops with a classed error (k190_falsifier) after
# the first cell in which a falsifier triggers, having written what it has.
#
# E2's generator is a copy of tests/testthat/helper-kit172.R:21-62 at d4e2eb2c (the branch point), verbatim except the
# names (kit172PickCells -> e2PickCells, kit172Values -> e2Values) and exactly one added argument `effect_share` (0 or
# 0.1), which replaces round(0.1 * n_features) at :44; the factor is planted whenever planted = TRUE.
GENERATOR_SOURCE <- list(file = "tests/testthat/helper-kit172.R", commit = "d4e2eb2c", lines = "21-62")

# ---- copy: helper-kit172.R:21-28 at d4e2eb2c (name changed) ----
e2PickCells <- function(allowed, groups, per_group, count) {
  cells <- which(allowed, arr.ind = TRUE)
  cells <- cells[sample.int(nrow(cells)), , drop = FALSE]
  key <- paste(cells[, 1L], groups[cells[, 2L]])
  rank <- stats::ave(seq_along(key), key, FUN = seq_along)
  cells <- cells[rank <= per_group, , drop = FALSE]
  cells[seq_len(min(count, nrow(cells))), , drop = FALSE]
}

# ---- copy: helper-kit172.R:34-62 at d4e2eb2c (name changed; `effect_share` added, :44 uses it; `factor_pattern` added
# by bound-declaration.md section 3: "random" draws the factor's sample pattern N(0, 1) per sample right after the
# loadings, centred and scaled to the balanced pattern's sample variance 12/11; "balanced" leaves the stream unchanged) ----
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
    factor_group_r = if (planted) stats::cor(pattern, as.numeric(groups == "B")) else NA_real_)
}

E2_CELLS <- c("N-r", "N-s", "E-r", "E-s", "E-p", "E-g")
ALL_CELLS <- c("E1", E2_CELLS)
# bound-declaration.md section 3: the two cells where the unwanted factor lines up with the groups by chance (cells 7, 8)
ALIGNED_CELLS <- c("N-w", "E-w")
CUTOFF <- 0.05

k190Abort <- function(class, message) stop(errorCondition(message, class = c(class, "k190_condition")))

# ---- E2 data (5.3): effects, the one-group features of E-g, the row order -------------------------------------------
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
    factor_group_r = sim$factor_group_r)
}

# ---- the rule of section 2.6, base R only (as the test file's oracle) ------------------------------------------------
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

# Reported only (VAL-100; D2's candidate rules, at the percentage the shipped sweep chose; nothing ships from this): the
# controls each candidate rule would select from the same BH family. (b) candidates q > cutoff ranked by p, highest
# first, ties broken by feature id; (c) the shipped ranking with features observed in fewer than two groups not
# candidates (they stay in the BH family); (d) both.
k190Candidate <- function(orc, rule) {
  candidates <- names(orc$q)[orc$q > CUTOFF]
  if (rule %in% c("c", "d")) candidates <- setdiff(candidates, orc$one_group)
  ranked <- if (rule %in% c("b", "d")) candidates[order(-orc$p[candidates], candidates, method = "radix")] else
    candidates[order(-orc$q[candidates])]
  if (!length(ranked)) return(character())
  utils::head(ranked, max(min(orc$count, length(ranked)), 1L))
}

# ---- one journey ----------------------------------------------------------------------------------------------------
k190Journey <- function(task, side, fault) {
  started <- proc.time()[[3L]]
  root <- file.path(tempdir(), sprintf("k190-%s-%s-%s-%s-%d", side, task$kind, gsub("[^A-Za-z0-9]", "", task$cell),
    task$class, task$r))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  object <- kit172Object(task$values, task$class)
  ruv <- kit172LegacyRuv("automatic")
  if (identical(task$grouping, "factor1")) {
    object@design_matrix$factor1 <- object@design_matrix$group
    ruv$ruv_grouping_variable <- "factor1"
  }
  run <- tryCatch(kit172Journey(root, object, task$class, ruv = ruv), error = function(e) e)
  base_row <- list(kind = task$kind, cell = task$cell, class = task$class, r = task$r, grouping = task$grouping,
    dataset = task$dataset %||% "", side = side)
  if (inherits(run, "error")) {
    return(c(base_row, list(error = conditionMessage(run), refused = inherits(run, "multischolar_ruv_grouping_unusable"),
      seconds = proc.time()[[3L]] - started)))
  }
  on.exit(closeDIAStateAdapter(run$owner), add = TRUE, after = FALSE)
  input <- kit172StateMatrix(run, paste0(run$prefix, "_normalized"))
  input[!is.finite(input)] <- NA_real_
  opt <- run$ruv$optimisation[["LCMS_Pos"]]
  groups <- kit172Design()$group
  valid <- rownames(input)[rowSums(!is.na(input)) > 1]
  pct <- opt$best_percentage %||% NA_real_
  controls <- if (isTRUE(opt$success)) names(opt$control_genes_index)[opt$control_genes_index] else character()
  if (fault && identical(side, "head") && identical(task$kind, "e2") && identical(task$class, "lipid") && task$r == 1L &&
    length(controls) > 1L)
    controls <- controls[-1L]
  first <- if (is.finite(pct)) utils::head(valid, max(round(length(valid) * pct / 100), 1L)) else character()
  orc <- if (is.finite(pct)) k190Oracle(input, groups, pct) else NULL
  tie_decided <- NA_real_
  if (!is.null(orc) && length(controls)) {
    q_cut <- min(orc$q[intersect(controls, names(orc$q))])
    block <- names(orc$q)[orc$q == q_cut & orc$q > CUTOFF]
    inside <- sum(controls %in% block)
    tie_decided <- if (length(block) > inside) inside / length(controls) else 0
  }
  observed_groups <- apply(!is.na(input), 1L, function(x) length(unique(groups[x])))
  complete <- rownames(input)[rowSums(is.na(input)) == 0L]
  one_group <- rownames(input)[observed_groups == 1L]
  da <- kit172DaTable(run)
  if (!is.null(da$assay)) da <- da[da$assay == "LCMS_Pos", , drop = FALSE]
  truth <- task$effects[da$feature]
  truth[is.na(truth)] <- 0
  called <- da$feature[!is.na(da$fdr_qvalue) & da$fdr_qvalue < 0.05]
  planted <- names(task$effects)[task$effects != 0]
  share <- function(set) if (length(controls)) sum(controls %in% set) / length(controls) else NA_real_
  c(base_row, list(error = "", refused = FALSE,
    n_v = length(valid), percentage = pct, n_controls = length(controls), k = opt$best_k %||% NA_integer_,
    dropped = !isTRUE(opt$success),
    share_first = share(first), share_planted = share(planted), share_complete = share(complete),
    share_tie_decided = tie_decided, share_one_group = share(one_group),
    share_one_group_features = if (length(valid)) sum(one_group %in% valid) / length(valid) else NA_real_,
    eg_survived = if (length(task$one_group)) sum(task$one_group %in% rownames(input)) else NA_integer_,
    calls = length(called), any_call = length(called) > 0L,
    fdp = if (length(called)) sum(!called %in% planted) / length(called) else 0,
    power = if (length(planted)) sum(called %in% planted) / length(planted) else NA_real_,
    abs_lfc_error = mean(abs(da$logFC - truth), na.rm = TRUE),
    matches_first = if (length(controls) && !is.null(orc)) setequal(controls, first) else NA,
    # K190_RULE=d (VAL-110, VAL-126): the head's lipid and metab controls are checked against rule (d) (2.7a).
    matches_oracle = if (length(controls) && !is.null(orc)) setequal(controls, if (identical(side, "head") &&
      identical(Sys.getenv("K190_RULE"), "d")) k190Candidate(orc, "d") else orc$selected) else NA,
    seconds = proc.time()[[3L]] - started,
    controls = list(controls), logfc = list(stats::setNames(da$logFC, da$feature))),
    rule_metrics(orc, controls, planted, complete, one_group))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Per candidate rule (b, c, d; reported only): its selection's size, its overlap (Jaccard) with the shipped rule's,
# and its shares of planted, one-group and complete features.
rule_metrics <- function(orc, controls, planted, complete, one_group) {
  out <- list()
  for (rule in c("b", "c", "d")) {
    sel <- if (is.null(orc)) character() else k190Candidate(orc, rule)
    shipped <- if (is.null(orc)) character() else orc$selected
    sh <- function(set) if (length(sel)) sum(sel %in% set) / length(sel) else NA_real_
    out[[paste0("rule_", rule, "_n")]] <- length(sel)
    out[[paste0("rule_", rule, "_jaccard_shipped")]] <- if (length(union(sel, shipped))) length(intersect(sel, shipped)) /
      length(union(sel, shipped)) else NA_real_
    out[[paste0("rule_", rule, "_share_planted")]] <- sh(planted)
    out[[paste0("rule_", rule, "_share_one_group")]] <- sh(one_group)
    out[[paste0("rule_", rule, "_share_complete")]] <- sh(complete)
  }
  out
}

e1Tasks <- function(lifecycle) {
  tasks <- list()
  add <- function(dataset, class, values, effects, r) tasks[[length(tasks) + 1L]] <<- list(kind = "e1", cell = "E1",
    class = class, r = r, grouping = "group", dataset = dataset, values = values, effects = effects, one_group = character())
  i <- 0L
  for (dataset in c("t11_concentration", "t11_f172")) for (class in c("lipid", "metab")) {
    values <- if (identical(dataset, "t11_f172")) kit172Values(172001L)$values else kit172ConcentrationValues()
    i <- i + 1L
    add(dataset, class, values, stats::setNames(rep(0, nrow(values)), rownames(values)), i)
  }
  if (!lifecycle) {
    seeds <- rbind(data.frame(seed = 172101:172120, class = "lipid"), data.frame(seed = 172101L, class = "metab"))
    for (j in seq_len(nrow(seeds))) for (missing in c(TRUE, FALSE)) {
      sim <- kit172Values(seeds$seed[[j]], planted = TRUE, missing = missing)
      i <- i + 1L
      add(sprintf("t1b_%d_%s", seeds$seed[[j]], if (missing) "missing" else "complete"), seeds$class[[j]], sim$values,
        sim$effects, i)
    }
  }
  tasks
}

e2Tasks <- function(cells, replicates, oracle_check) {
  tasks <- list()
  for (cell in intersect(cells, E2_CELLS)) for (class in c("lipid", "metab")) for (r in seq_len(replicates)) {
    data <- e2Data(cell, r)
    tasks[[length(tasks) + 1L]] <- list(kind = "e2", cell = cell, class = class, r = r, grouping = "group",
      values = data$values, effects = data$effects, one_group = data$one_group, seed = data$seed)
    if (oracle_check && identical(cell, "N-r")) {
      tasks[[length(tasks) + 1L]] <- list(kind = "oracle_check", cell = cell, class = class, r = r, grouping = "factor1",
        values = data$values, effects = data$effects, one_group = data$one_group, seed = data$seed)
    }
  }
  tasks
}

loadTree <- function(tree, head_tree) {
  suppressPackageStartupMessages(pkgload::load_all(tree, export_all = TRUE, helpers = FALSE, quiet = TRUE))
  logger::log_threshold(logger::FATAL)
  # The head tree's helpers in their own environment, attached behind the global one (no name here is overwritten).
  helpers <- new.env(parent = globalenv())
  suppressMessages(testthat::source_test_helpers(file.path(head_tree, "tests", "testthat"), env = helpers))
  attach(helpers, name = "k190_helpers", warn.conflicts = FALSE)
}

runSide <- function(tree, head_tree, side, work, cells, replicates, lifecycle, fault, workers) {
  loadTree(tree, head_tree)
  tasks <- c(if ("E1" %in% cells) e1Tasks(lifecycle), e2Tasks(cells, replicates, lifecycle && identical(side, "head")))
  results <- parallel::mclapply(tasks, function(task) tryCatch(
    suppressWarnings(suppressMessages(k190Journey(task, side, fault))),
    error = function(e) list(kind = task$kind, cell = task$cell, class = task$class, r = task$r, side = side,
      grouping = task$grouping, dataset = task$dataset %||% "", error = paste("instrument:", conditionMessage(e)),
      refused = FALSE)),
    mc.cores = workers, mc.preschedule = FALSE)
  commit <- system2("git", c("-C", shQuote(tree), "rev-parse", "HEAD"), stdout = TRUE)
  saveRDS(list(side = side, commit = commit, results = results), file.path(work, paste0(side, ".rds")))
  failed <- Filter(function(x) inherits(x, "try-error") || startsWith(x$error %||% "", "instrument:"), results)
  if (length(failed)) {
    cat(sprintf("side %s: %d journey(s) failed in the instrument; first: %s\n", side, length(failed),
      if (inherits(failed[[1L]], "try-error")) as.character(failed[[1L]]) else failed[[1L]]$error))
    quit(status = 1L)
  }
  cat(sprintf("side %s finished: %d journeys\n", side, length(results)))
}

# ---- E3 side: timed journeys at realistic sizes ----------------------------------------------------------------------
runE3Side <- function(tree, head_tree, side, work, sizes, modes) {
  loadTree(tree, head_tree)
  rows <- list()
  for (n in sizes) for (mode in modes) {
    values <- list(LCMS_Pos = e2Values(190900001, planted = TRUE, n = n, effect_share = 0.1)$values,
      LCMS_Neg = e2Values(190900002, planted = TRUE, n = n, effect_share = 0.1)$values)
    root <- file.path(tempdir(), sprintf("k190-e3-%s-%d-%s", side, n, mode))
    started <- proc.time()[[3L]]
    ruv_seconds <- NA_real_
    error <- ""
    outcome <- tryCatch(suppressWarnings(suppressMessages({
      run <- kit172Journey(root, kit172Object(values, "metab"), "metab", ruv = NULL, da = FALSE)
      config <- run$fixture$workflow$config_list
      ruv_seconds <- system.time(runWorkflowOmicsReferenceRuv(run$binding, run$owner, kit172LegacyRuv(mode), config))[["elapsed"]]
      runWorkflowOmicsReferenceDaAnalysis(run$binding, run$owner, getWorkflowStateReference(run$owner),
        list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), kit172Contrasts())
      closeDIAStateAdapter(run$owner)
      TRUE
    })), error = function(e) conditionMessage(e))
    if (!isTRUE(outcome)) error <- outcome
    unlink(root, recursive = TRUE)
    rows[[length(rows) + 1L]] <- data.frame(side = side, n_per_assay = n, assays = 2L, mode = mode, ruv_seconds = ruv_seconds,
      journey_seconds = proc.time()[[3L]] - started, error = gsub("[\t\n]", " ", error), stringsAsFactors = FALSE)
  }
  commit <- system2("git", c("-C", shQuote(tree), "rev-parse", "HEAD"), stdout = TRUE)
  saveRDS(list(side = side, commit = commit, rows = do.call(rbind, rows)), file.path(work, paste0("e3-", side, ".rds")))
}

# ==== bound mode (bound-declaration.md, VAL-182, VAL-185) ============================================================
# Arms: unbounded B (legacy, d4e2eb2c), F (first repair, d6d67efe), D (rule (d), a1c12d28); bounded B+B, F+B, D+B through
# the legacy Manual mode at the unbounded arm's chosen percentage with k = min(k, 2). Lipid only. Old cells (1 to 6)
# and E1 run F, D; the aligned cells N-w, E-w run B, F, D. Unbounded rows for r <= 50 in the old cells are read from
# the committed e/ (F, B) and ed/ (D) tables; otherwise the unbounded journey is run.
BOUND_TREES <- c(B = "d4e2eb2cff0302c27f98b78f90bfb51dab6aee04", F = "d6d67efea7e7a09f8ba41c28a99718dddb2ec796",
  D = "a1c12d28bfd6e8291f3e360ab75b6875458d8c88")
BOUND_CAP <- 2L
BOUND_NULL_CELLS <- c("N-r", "N-s", "N-w")
BOUND_EFFECT_CELLS <- c("E-r", "E-s", "E-p", "E-g", "E-w")
BOUND_DECIDING <- c(BOUND_NULL_CELLS, BOUND_EFFECT_CELLS)

boundKey <- function(cell, r, dataset) if (identical(cell, "E1")) paste("E1", dataset, sep = "|") else paste(cell, r, sep = "|")

boundArms <- function(cell) if (cell %in% ALIGNED_CELLS) c("B", "F", "D") else c("F", "D")

boundE1Data <- function(dataset) {
  if (identical(dataset, "t11_concentration")) {
    values <- kit172ConcentrationValues()
    return(list(values = values, effects = stats::setNames(rep(0, nrow(values)), rownames(values)), one_group = character()))
  }
  if (identical(dataset, "t11_f172")) {
    values <- kit172Values(172001L)$values
    return(list(values = values, effects = stats::setNames(rep(0, nrow(values)), rownames(values)), one_group = character()))
  }
  parts <- strsplit(dataset, "_", fixed = TRUE)[[1L]]
  sim <- kit172Values(as.integer(parts[[2L]]), planted = TRUE, missing = identical(parts[[3L]], "missing"))
  list(values = sim$values, effects = sim$effects, one_group = character())
}

boundE1Datasets <- function(all = TRUE) c("t11_concentration", "t11_f172", if (all) as.vector(outer(sprintf("t1b_%d", 172101:172120),
  c("missing", "complete"), paste, sep = "_")))

# One journey of arm `arm` (unbounded when `ruv` is Automatic, bounded when Manual), with the declaration's endpoints,
# its checks and its instrumentation. `rule` names the arm's control rule for the oracle check.
boundJourney <- function(task, arm, ruv, rule, k_ref = NA_real_) {
  started <- proc.time()[[3L]]
  root <- file.path(tempdir(), sprintf("k190b-%s-%s-%s-%d-%s", arm, gsub("[^A-Za-z0-9]", "", task$cell),
    gsub("[^A-Za-z0-9]", "", task$dataset), task$r, ruv$ruv_mode))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  data <- if (identical(task$cell, "E1")) boundE1Data(task$dataset) else e2Data(task$cell, task$r)
  object <- kit172Object(data$values, "lipid")
  head_row <- list(cell = task$cell, r = task$r, dataset = task$dataset, arm = arm, mode = ruv$ruv_mode,
    factor_group_r = data$factor_group_r %||% NA_real_)
  run <- tryCatch(kit172Journey(root, object, "lipid", ruv = ruv), error = function(e) e)
  if (inherits(run, "error")) {
    return(c(head_row, list(error = gsub("[\t\n]", " ", conditionMessage(run)), seconds = proc.time()[[3L]] - started)))
  }
  on.exit(closeDIAStateAdapter(run$owner), add = TRUE, after = FALSE)
  input <- kit172StateMatrix(run, paste0(run$prefix, "_normalized"))
  input[!is.finite(input)] <- NA_real_
  opt <- run$ruv$optimisation[["LCMS_Pos"]]
  groups <- kit172Design()$group
  valid <- rownames(input)[rowSums(!is.na(input)) > 1]
  pct <- opt$best_percentage %||% NA_real_
  k <- opt$best_k %||% NA_integer_
  controls <- if (isTRUE(opt$success)) names(opt$control_genes_index)[opt$control_genes_index] else character()
  orc <- if (is.finite(pct)) k190Oracle(input, groups, pct) else NULL
  oracle <- if (is.null(orc)) character() else switch(rule,
    B = utils::head(valid, max(round(length(valid) * pct / 100), 1L)),
    F = orc$selected,
    D = k190Candidate(orc, "d"))
  # Residual dimensions per feature: observed runs minus observed replicate sets (the groups).
  observed <- !is.na(input)
  dims <- rowSums(observed) - apply(observed, 1L, function(x) length(unique(groups[x])))
  dimCounts <- function(kk) if (is.finite(kk)) c(sum(dims < kk), sum(dims == kk), sum(dims > kk)) else rep(NA_integer_, 3L)
  corrected <- tryCatch(kit172StateMatrix(run, paste0(run$prefix, "_ruv_corrected")), error = function(e) NULL)
  removed <- if (is.null(corrected)) NA_integer_ else sum(!rownames(input) %in% rownames(corrected)[rowSums(!is.na(corrected)) > 0])
  curve <- tryCatch(.extractCancorCurveData(opt$cancor_plot), error = function(e) NULL)
  separation <- tryCatch(calculateSeparationScore(opt$cancor_plot, metric = "max_difference"), error = function(e) NA_real_)
  da <- kit172DaTable(run)
  if (!is.null(da$assay)) da <- da[da$assay == "LCMS_Pos", , drop = FALSE]
  truth <- data$effects[da$feature]
  truth[is.na(truth)] <- 0
  called <- da$feature[!is.na(da$fdr_qvalue) & da$fdr_qvalue < 0.05]
  planted <- names(data$effects)[data$effects != 0]
  at_k <- dimCounts(k)
  at_ref <- dimCounts(k_ref)
  c(head_row, list(error = "", percentage = pct, k = k, n_controls = length(controls), dropped = !isTRUE(opt$success),
    oracle_ok = if (length(controls) && !is.null(orc)) setequal(controls, oracle) else NA,
    calls = length(called), any_call = length(called) > 0L,
    fdp = if (length(called)) sum(!called %in% planted) / length(called) else 0,
    power = if (length(planted)) sum(called %in% planted) / length(planted) else NA_real_,
    abs_lfc_error = mean(abs(da$logFC - truth), na.rm = TRUE),
    dims_below_k = at_k[[1L]], dims_at_k = at_k[[2L]], dims_above_k = at_k[[3L]],
    dims_below_kref = at_ref[[1L]], dims_at_kref = at_ref[[2L]], dims_above_kref = at_ref[[3L]],
    removed_by_ruv = removed, separation = separation,
    curve = if (is.null(curve)) "" else paste(sprintf("%s:%.6g:%.6g", curve$K, curve$All, curve$Control), collapse = ";"),
    seconds = proc.time()[[3L]] - started))
}

# A side: one tree, one rule; per task the unbounded journey (unless reused) then the bounded one.
runBoundSide <- function(tree, head_tree, arm, work, tasks_file, workers, fault) {
  loadTree(tree, head_tree)
  tasks <- readRDS(tasks_file)
  results <- parallel::mclapply(seq_along(tasks), function(i) tryCatch(suppressWarnings(suppressMessages({
    task <- tasks[[i]]
    rows <- list()
    pct_u <- task$pct_u
    k_u <- task$k_u
    if (isTRUE(task$run_unbounded)) {
      u <- boundJourney(task, arm, kit172LegacyRuv("automatic"), arm)
      rows[[1L]] <- u
      if (!isTRUE(task$reuse)) {
        pct_u <- u$percentage %||% NA_real_
        k_u <- u$k %||% NA_real_
      }
    }
    if (is.finite(pct_u %||% NA_real_) && is.finite(k_u %||% NA_real_)) {
      k_b <- as.integer(min(k_u, BOUND_CAP))
      if (fault && i == 1L) k_b <- k_b + 1L
      ruv <- kit172LegacyRuv("manual")
      ruv$manual_percentage <- pct_u
      ruv$manual_k <- k_b
      b <- boundJourney(task, paste0(arm, "+B"), ruv, arm, k_ref = k_u)
      b$pct_u <- pct_u
      b$k_u <- k_u
      b$k_ok <- isTRUE(b$k == as.integer(min(k_u, BOUND_CAP)))
      b$pct_ok <- isTRUE(b$percentage == pct_u)
      rows[[length(rows) + 1L]] <- b
    }
    rows
  })), error = function(e) list(list(cell = tasks[[i]]$cell, r = tasks[[i]]$r, dataset = tasks[[i]]$dataset, arm = arm,
    error = paste("instrument:", gsub("[\t\n]", " ", conditionMessage(e)))))),
    mc.cores = workers, mc.preschedule = FALSE)
  saveRDS(list(arm = arm, commit = system2("git", c("-C", shQuote(tree), "rev-parse", "HEAD"), stdout = TRUE),
    rows = unlist(results, recursive = FALSE)), file.path(work, paste0("bound-", arm, ".rds")))
}

boundFlat <- function(rows) {
  cols <- c("cell", "r", "dataset", "arm", "mode", "error", "factor_group_r", "pct_u", "k_u", "percentage", "k", "n_controls",
    "dropped", "oracle_ok", "k_ok", "pct_ok", "calls", "any_call", "fdp", "power", "abs_lfc_error", "dims_below_k", "dims_at_k",
    "dims_above_k", "dims_below_kref", "dims_at_kref", "dims_above_kref", "removed_by_ruv", "separation", "curve", "seconds")
  do.call(rbind, lapply(rows, function(x) {
    out <- lapply(stats::setNames(cols, cols), function(cn) {
      v <- x[[cn]]
      if (is.null(v) || length(v) == 0L) NA else v[[1L]]
    })
    as.data.frame(out, stringsAsFactors = FALSE)
  }))
}

# The committed unbounded rows (lipid): F from e/ (head), B from e/ (base), D from ed/ (head).
boundCommitted <- function(root) {
  read <- function(dir) do.call(rbind, lapply(file.path(root, "docs/qa/workflow-execution/kit/kit-190", dir, c("a", "b"),
    "journeys.tsv"), function(f) if (file.exists(f)) utils::read.delim(f, stringsAsFactors = FALSE) else NULL))
  e <- read("e")
  ed <- read("ed")
  pick <- function(t, side) {
    ok <- is.na(t$error) | !nzchar(as.character(t$error))
    t <- t[t$side == side & t$class == "lipid" & t$kind %in% c("e1", "e2") & ok, , drop = FALSE]
    t$key <- vapply(seq_len(nrow(t)), function(i) boundKey(t$cell[[i]], t$r[[i]], as.character(t$dataset[[i]])), "")
    t[, c("key", "percentage", "k", "n_controls", "calls")]
  }
  list(F = pick(e, "head"), B = pick(e, "base"), D = pick(ed, "head"))
}

runBound <- function(phase, out) {
  started_at <- Sys.time()
  lifecycle <- startsWith(phase, "lc")
  blinded <- phase %in% c("stage1", "lc-stage1")
  if (blinded) options(warn = -1L) # no warning text reaches the console in a blinded run
  plant <- identical(Sys.getenv("K190_PLANT"), "1")
  fault <- identical(Sys.getenv("K190_FAULT"), "1")
  stopifnot(phase %in% c("lc", "lc-stage1", "stage1", "stage2"), !plant || identical(phase, "lc-stage1"),
    !fault || identical(phase, "lc"))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out <- normalizePath(out)
  script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))
  head_root <- normalizePath(".")
  # Private working directory: discarded with the runner, never under `out`.
  work <- file.path(tempdir(), "k190-bound-work")
  dir.create(work, showWarnings = FALSE)
  say <- function(...) if (!blinded) cat(..., "\n", sep = "")
  committed <- boundCommitted(head_root)
  plan <- NULL
  if (identical(phase, "stage2")) {
    plan_file <- Sys.getenv("K190_PLAN")
    stopifnot(nzchar(plan_file), file.exists(plan_file))
    plan <- utils::read.delim(plan_file, stringsAsFactors = FALSE)
  }
  cells_env <- Sys.getenv("K190_CELLS", paste(c(BOUND_DECIDING, if (identical(phase, "stage2") || identical(phase, "lc")) "E1"),
    collapse = ","))
  cells <- strsplit(cells_env, ",", fixed = TRUE)[[1L]]
  stopifnot(all(cells %in% c(BOUND_DECIDING, "E1")), !blinded || !("E1" %in% cells))
  R_of <- function(cell) {
    if (lifecycle) return(2L)
    if (identical(phase, "stage1")) return(50L)
    as.integer(plan$R[plan$cell == cell])
  }
  tasks <- list(B = list(), F = list(), D = list())
  for (cell in cells) {
    entries <- if (identical(cell, "E1")) lapply(boundE1Datasets(!lifecycle), function(d) list(r = 1L, dataset = d)) else
      lapply(seq_len(R_of(cell)), function(r) list(r = r, dataset = ""))
    for (arm in boundArms(cell)) for (en in entries) {
      key <- boundKey(cell, en$r, en$dataset)
      ref <- committed[[arm]][committed[[arm]]$key == key, , drop = FALSE]
      old <- !(cell %in% ALIGNED_CELLS) && nrow(ref) == 1L && is.finite(ref$percentage) && is.finite(ref$k)
      # Only the lifecycle's own check (i) re-runs reusable unbounded journeys; every other phase reuses them.
      reuse_now <- old && !identical(phase, "lc")
      task <- list(cell = cell, r = en$r, dataset = en$dataset, reuse = reuse_now, run_unbounded = !reuse_now, pct_u = if (old) ref$percentage else NA_real_, k_u = if (old) ref$k else NA_real_,
        committed = if (old) as.list(ref) else NULL)
      tasks[[arm]][[length(tasks[[arm]]) + 1L]] <- task
    }
  }
  arms <- names(tasks)[vapply(tasks, length, 1L) > 0L]
  trees <- list()
  for (arm in arms) {
    trees[[arm]] <- file.path(tempdir(), paste0("k190-tree-", arm))
    invisible(system2("git", c("worktree", "add", "--detach", shQuote(trees[[arm]]), BOUND_TREES[[arm]]), stdout = TRUE, stderr = TRUE))
    stopifnot(dir.exists(file.path(trees[[arm]], "R")))
    pkgbuild::compile_dll(trees[[arm]], quiet = TRUE)
  }
  workers <- as.integer(Sys.getenv("K190_WORKERS", if (length(arms) == 3L) "1" else "2"))
  rscript <- file.path(R.home("bin"), "Rscript")
  jobs <- lapply(arms, function(arm) {
    tf <- file.path(work, paste0("tasks-", arm, ".rds"))
    saveRDS(tasks[[arm]], tf)
    parallel::mcparallel(system2(rscript, c(shQuote(script_file), "--bound-side", shQuote(trees[[arm]]), shQuote(head_root), arm,
      shQuote(work), shQuote(tf), workers, if (fault) "1" else "0"),
      stdout = file.path(work, paste0("bound-", arm, ".log")), stderr = file.path(work, paste0("bound-", arm, ".log"))))
  })
  status <- unlist(parallel::mccollect(jobs))
  files <- file.path(work, paste0("bound-", arms, ".rds"))
  if (any(status != 0L) || !all(file.exists(files))) {
    if (!blinded) for (arm in arms) {
      log <- file.path(work, paste0("bound-", arm, ".log"))
      if (file.exists(log)) cat(paste0("[", arm, "] ", utils::tail(readLines(log, warn = FALSE), 15L)), sep = "\n")
    }
    k190Abort("k190_instrument", "a bound side failed")
  }
  rows <- boundFlat(unlist(lapply(files, function(f) readRDS(f)$rows), recursive = FALSE))
  instrument <- sum(startsWith(as.character(rows$error), "instrument:"), na.rm = TRUE)
  journey_errors <- sum(!is.na(rows$error) & nzchar(as.character(rows$error)), na.rm = TRUE) - instrument
  bounded <- rows[grepl("\\+B$", rows$arm) & !is.na(rows$mode), , drop = FALSE]
  falsifiers <- c(k_readback = sum(!bounded$k_ok, na.rm = TRUE), pct_readback = sum(!bounded$pct_ok, na.rm = TRUE),
    oracle = sum(!bounded$oracle_ok, na.rm = TRUE), instrument = instrument)
  seconds_median <- vapply(c(bounded = "\\+B$", B = "^B$", F = "^F$", D = "^D$"), function(rx) {
    x <- rows$seconds[grepl(rx, rows$arm)]
    if (length(x)) stats::median(x) else NA_real_
  }, numeric(1L))
  counts <- as.data.frame(table(cell = rows$cell, arm = rows$arm), stringsAsFactors = FALSE)
  provenance <- list(ticket = "KIT-190 bound", phase = phase, head_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    trees = as.list(BOUND_TREES[arms]), script_sha256 = strsplit(system2("sha256sum", shQuote(script_file), stdout = TRUE), " ")[[1L]][[1L]],
    seed_rule = "190000000 + 10000 * cell + r; cells N-r 1 to E-g 6, N-w 7, E-w 8", cells = cells,
    journeys = stats::setNames(as.list(counts$Freq), paste(counts$cell, counts$arm)), falsifiers = as.list(falsifiers),
    journey_errors = journey_errors,
    seconds_median_per_journey = as.list(seconds_median),
    wall_seconds = as.numeric(difftime(Sys.time(), started_at, units = "secs")))
  if (blinded) {
    # Stage 1: the rows stay in `work`; the re-estimation sees only the capped arms' endpoints.
    capped <- bounded[bounded$cell %in% BOUND_DECIDING & bounded$arm %in% c("F+B", "D+B"), c("cell", "arm", "r", "fdp", "any_call")]
    if (plant) {
      shift <- capped$cell == "E-r" & capped$arm == "D+B"
      capped$fdp[shift] <- capped$fdp[shift] + 0.30
      nr <- capped$cell == "N-r"
      capped$arm[nr] <- ifelse(capped$arm[nr] == "F+B", "D+B", "F+B")
    }
    rows_file <- file.path(work, "capped.rds")
    saveRDS(capped, rows_file)
    code <- system2(rscript, c(shQuote(file.path(dirname(script_file), "reestimate.R")), "estimate", shQuote(rows_file), shQuote(out)),
      stdout = file.path(work, "reestimate.log"), stderr = file.path(work, "reestimate.log"))
    jsonlite::write_json(provenance, file.path(out, "provenance.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
    for (i in seq_len(nrow(counts))) cat(sprintf("journeys %s %s: %d\n", counts$cell[[i]], counts$arm[[i]], counts$Freq[[i]]))
    for (nm in names(falsifiers)) cat(sprintf("falsifier %s: %d\n", nm, falsifiers[[nm]]))
    cat(sprintf("seconds: %.0f\n", provenance$wall_seconds))
    if (code == 0L && file.exists(file.path(out, "reestimate.tsv"))) {
      cat(readLines(file.path(out, "reestimate.tsv")), sep = "\n")
    } else cat("reestimate: failed\n")
    if (code != 0L || any(falsifiers > 0L)) quit(status = 1L)
    return(invisible(NULL))
  }
  utils::write.table(rows[, setdiff(names(rows), "curve")], file.path(out, "journeys.tsv"), sep = "\t", quote = FALSE,
    row.names = FALSE, na = "NA")
  utils::write.table(rows[, c("cell", "r", "dataset", "arm", "mode", "k", "separation", "curve")], file.path(out, "curves.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  checks <- list()
  if (lifecycle) {
    # (i) the unbounded journeys of the old cells reproduce the committed rows (k, percentage, controls count, calls).
    unb <- rows[rows$mode == "automatic" & rows$cell %in% c(E2_CELLS, "E1") & !is.na(rows$mode), , drop = FALSE]
    rep_ok <- vapply(seq_len(nrow(unb)), function(i) {
      key <- boundKey(unb$cell[[i]], unb$r[[i]], as.character(unb$dataset[[i]]))
      ref <- committed[[unb$arm[[i]]]]
      ref <- ref[ref$key == key, , drop = FALSE]
      nrow(ref) == 1L && isTRUE(ref$k == unb$k[[i]]) && isTRUE(ref$percentage == unb$percentage[[i]]) &&
        isTRUE(ref$n_controls == unb$n_controls[[i]]) && isTRUE(ref$calls == unb$calls[[i]])
    }, logical(1L))
    checks$i_reproduce <- c(ok = sum(rep_ok), of = length(rep_ok))
    checks$ii_k_pct_readback <- c(ok = sum(bounded$k_ok & bounded$pct_ok, na.rm = TRUE), of = nrow(bounded))
    checks$iii_oracle <- c(ok = sum(bounded$oracle_ok, na.rm = TRUE), of = nrow(bounded))
    filled <- !is.na(bounded$dims_at_k) & !is.na(bounded$removed_by_ruv) & nzchar(bounded$curve %||% "")
    checks$iv_instrumentation <- c(ok = sum(filled), of = nrow(bounded))
    secs <- seconds_median
    secs[!is.finite(secs)] <- 0
    cR <- c(`N-r` = 100, `N-s` = 100, `N-w` = 250, `E-r` = 200, `E-s` = 200, `E-p` = 200, `E-g` = 200, `E-w` = 200)
    old_cost <- function(R) 2 * secs[["bounded"]] * R + (secs[["F"]] + secs[["D"]]) * max(R - 50, 0)
    projection <- sum(vapply(names(cR), function(c) if (c %in% ALIGNED_CELLS) cR[[c]] * (3 * secs[["bounded"]] + secs[["B"]] +
      secs[["F"]] + secs[["D"]]) else old_cost(cR[[c]]), numeric(1L))) + 42 * 2 * secs[["bounded"]]
    checks$vii_projected_core_hours_option_c <- c(value = projection / 3600 * 1.25, ceiling = 100)
  }
  jsonlite::write_json(c(provenance, list(checks = lapply(checks, as.list))), file.path(out, "provenance.json"),
    auto_unbox = TRUE, pretty = TRUE, digits = NA)
  cat(sprintf("KIT-190 bound %s: %d rows; falsifiers %s\n", phase, nrow(rows), paste(names(falsifiers), falsifiers, collapse = ", ")))
  for (nm in names(checks)) cat(sprintf("check %s: %s\n", nm, paste(names(checks[[nm]]), round(checks[[nm]], 2), collapse = " ")))
  if (fault) {
    if (falsifiers[["k_readback"]] > 0L) k190Abort("k190_falsifier", "fault run: k read-back falsifier triggered (as planted)")
    k190Abort("k190_instrument", "fault run: the planted k + 1 was not detected")
  }
  if (any(falsifiers > 0L)) k190Abort("k190_falsifier", sprintf("bound %s: falsifiers %s", phase,
    paste(names(falsifiers), falsifiers, collapse = ", ")))
  invisible(NULL)
}

# ---- the side processes ---------------------------------------------------------------------------------------------
args <- commandArgs(TRUE)
if (length(args) && identical(args[[1L]], "--side")) {
  runSide(args[[2L]], args[[3L]], args[[4L]], args[[5L]], strsplit(args[[6L]], ",", fixed = TRUE)[[1L]],
    as.integer(args[[7L]]), identical(args[[8L]], "1"), identical(args[[9L]], "1"), as.integer(args[[10L]]))
  quit(status = 0L)
}
if (length(args) && identical(args[[1L]], "--e3-side")) {
  runE3Side(args[[2L]], args[[3L]], args[[4L]], args[[5L]], as.integer(strsplit(args[[6L]], ",", fixed = TRUE)[[1L]]),
    strsplit(args[[7L]], ",", fixed = TRUE)[[1L]])
  quit(status = 0L)
}

# ---- the controller -------------------------------------------------------------------------------------------------
if (length(args) && identical(args[[1L]], "--bound-side")) {
  runBoundSide(args[[2L]], args[[3L]], args[[4L]], args[[5L]], args[[6L]], as.integer(args[[7L]]), identical(args[[8L]], "1"))
  quit(status = 0L)
}
if (length(args) && identical(args[[1L]], "bound")) {
  # Usage: Rscript effect.R bound <lc | lc-stage1 | stage1 | stage2> <out dir>   (bound-declaration.md)
  stopifnot(length(args) == 3L)
  runBound(args[[2L]], args[[3L]])
  quit(status = 0L)
}
stopifnot(length(args) == 3L, args[[1L]] %in% c("e", "e3"), grepl("^([0-9a-f]{40}|self)$", args[[2L]]))
mode <- args[[1L]]
started_at <- Sys.time()
lifecycle <- identical(Sys.getenv("K190_LIFECYCLE"), "1")
fault <- identical(Sys.getenv("K190_FAULT"), "1")
stop_rule <- !lifecycle || fault
cells <- strsplit(Sys.getenv("K190_CELLS", paste(ALL_CELLS, collapse = ",")), ",", fixed = TRUE)[[1L]]
if (fault) cells <- "N-r"
stopifnot(length(cells) > 0L, all(cells %in% ALL_CELLS))
replicates <- if (lifecycle) 2L else 50L
workers <- as.integer(Sys.getenv("K190_WORKERS", "2"))
out <- args[[3L]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))
head_root <- normalizePath(".")
head_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
base_root <- head_root
if (!identical(args[[2L]], "self")) {
  base_root <- file.path(tempdir(), "k190-base")
  invisible(system2("git", c("worktree", "add", "--detach", shQuote(base_root), args[[2L]]), stdout = TRUE, stderr = TRUE))
  stopifnot(dir.exists(file.path(base_root, "R")))
}
base_commit <- system2("git", c("-C", shQuote(base_root), "rev-parse", "HEAD"), stdout = TRUE)
cat("KIT-190", mode, ": head", head_commit, "base", base_commit, "lifecycle", lifecycle, "fault", fault, "cells",
  paste(cells, collapse = ","), "R", replicates, "out", out, "\n")
for (tree in unique(c(base_root, head_root))) pkgbuild::compile_dll(tree, quiet = TRUE)
work <- file.path(tempdir(), "k190-work")
dir.create(work, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")
sha <- function(f) strsplit(system2("sha256sum", shQuote(f), stdout = TRUE), " ")[[1L]][[1L]]
environment_facts <- function() list(r_version = R.version.string, blas = extSoftVersion()[["BLAS"]], lapack = La_library(),
  runner_image = paste(Sys.getenv("ImageOS"), Sys.getenv("ImageVersion")), cores = parallel::detectCores(),
  processor = tryCatch(sub("^model name\\s*:\\s*", "", grep("^model name", readLines("/proc/cpuinfo"), value = TRUE)[[1L]]),
    error = function(e) NA), script_sha256 = sha(script_file))
showLog <- function(name) {
  log <- file.path(work, paste0(name, ".log"))
  if (file.exists(log)) cat(paste0("[", name, "] ", utils::tail(readLines(log, warn = FALSE), 20L)), sep = "\n")
}

if (identical(mode, "e3")) {
  sizes <- if (lifecycle) 520L else as.integer(strsplit(Sys.getenv("K190_E3_N", "5000,20000"), ",", fixed = TRUE)[[1L]])
  dir.create(file.path(out, "e3"), showWarnings = FALSE)
  plan <- list(list(side = "head", tree = head_root, sizes = sizes, modes = if (lifecycle) "automatic" else c("automatic", "manual")))
  if (!lifecycle) plan[[2L]] <- list(side = "base", tree = base_root, sizes = max(sizes), modes = "automatic")
  status <- integer()
  for (step in plan) {
    status[[step$side]] <- system2(rscript, c(shQuote(script_file), "--e3-side", shQuote(step$tree), shQuote(head_root),
      step$side, shQuote(work), paste(step$sizes, collapse = ","), paste(step$modes, collapse = ",")),
      stdout = file.path(work, paste0("e3-", step$side, ".log")), stderr = file.path(work, paste0("e3-", step$side, ".log")))
    showLog(paste0("e3-", step$side))
  }
  files <- file.path(work, paste0("e3-", names(status), ".rds"))
  if (any(status != 0L) || !all(file.exists(files)))
    k190Abort("k190_instrument", sprintf("an E3 side failed (exit status %s)", paste(status, collapse = ", ")))
  times <- do.call(rbind, lapply(files, function(f) readRDS(f)$rows))
  times$commit <- ifelse(times$side == "head", head_commit, base_commit)
  # Projection: the Automatic RUV step's seconds per feature at the largest measured size, times 50,000 features in all.
  auto <- times[times$side == "head" & times$mode == "automatic" & is.finite(times$ruv_seconds), , drop = FALSE]
  largest <- auto[which.max(auto$n_per_assay), , drop = FALSE]
  projected <- if (nrow(largest)) largest$ruv_seconds / (largest$n_per_assay * largest$assays) * 50000 else NA_real_
  times <- rbind(times, data.frame(side = "head", n_per_assay = 25000L, assays = 2L, mode = "automatic_projected_50000",
    ruv_seconds = projected, journey_seconds = NA_real_, error = "", commit = head_commit, stringsAsFactors = FALSE))
  times$worker_deadline_seconds <- 1800
  utils::write.table(times, file.path(out, "e3", "times.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  jsonlite::write_json(c(list(ticket = "KIT-190 E3", head_commit = head_commit, base_commit = base_commit, lifecycle = lifecycle,
    sizes = sizes, seeds = c(LCMS_Pos = 190900001, LCMS_Neg = 190900002), generator = GENERATOR_SOURCE,
    projection_rule = "head Automatic RUV seconds / (features per assay x 2) x 50000 at the largest measured size",
    wall_seconds = as.numeric(difftime(Sys.time(), started_at, units = "secs"))), environment_facts()),
    file.path(out, "e3", "provenance.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
  cat(sprintf("KIT-190 E3 written: %d rows; projected Automatic RUV at 50,000 features %.0f s (deadline 1,800 s)\n",
    nrow(times), projected))
  quit(status = 0L)
}

sides <- list(base = base_root, head = head_root)
jobs <- lapply(names(sides), function(side) parallel::mcparallel(system2(rscript, c(shQuote(script_file), "--side",
  shQuote(sides[[side]]), shQuote(head_root), side, shQuote(work), paste(cells, collapse = ","), replicates,
  if (lifecycle) "1" else "0", if (fault) "1" else "0", workers),
  stdout = file.path(work, paste0(side, ".log")), stderr = file.path(work, paste0(side, ".log")))))
status <- unlist(parallel::mccollect(jobs))
for (side in names(sides)) showLog(side)
wall_sides <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))
if (any(status != 0L) || !all(file.exists(file.path(work, c("base.rds", "head.rds")))))
  k190Abort("k190_instrument", sprintf("a side failed (exit status %s)", paste(status, collapse = ", ")))
data <- lapply(c(base = "base", head = "head"), function(side) readRDS(file.path(work, paste0(side, ".rds"))))

scalar <- c("n_v", "percentage", "n_controls", "k", "share_first", "share_planted", "share_complete", "share_tie_decided",
  "share_one_group", "share_one_group_features", "eg_survived", "calls", "fdp", "power", "abs_lfc_error", "seconds",
  as.vector(outer(c("rule_b", "rule_c", "rule_d"), c("_n", "_jaccard_shipped", "_share_planted", "_share_one_group",
    "_share_complete"), paste0)))
flags <- c("dropped", "refused", "any_call", "matches_first", "matches_oracle")
flat <- function(results) do.call(rbind, lapply(results, function(x) {
  row <- data.frame(kind = x$kind, cell = x$cell, class = x$class, r = x$r, grouping = x$grouping, dataset = x$dataset %||% "",
    side = x$side, error = gsub("[\t\n]", " ", x$error %||% ""), stringsAsFactors = FALSE)
  for (name in scalar) row[[name]] <- as.numeric(x[[name]] %||% NA_real_)
  for (name in flags) row[[name]] <- as.logical(x[[name]] %||% NA)
  row
}))
rows <- rbind(flat(data$base$results), flat(data$head$results))
utils::write.table(rows, file.path(out, "journeys.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
key <- function(x) paste(x$kind, x$cell, x$class, x$r, x$grouping, x$dataset %||% "")
controlsOf <- function(side) stats::setNames(lapply(data[[side]]$results, function(x) x$controls[[1L]] %||% character()),
  vapply(data[[side]]$results, key, ""))
logfcOf <- function(side) stats::setNames(lapply(data[[side]]$results, function(x) x$logfc[[1L]] %||% numeric()),
  vapply(data[[side]]$results, key, ""))

summarise <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  m <- if (n) mean(x) else NA_real_
  half <- if (n > 1L) 1.96 * stats::sd(x) / sqrt(n) else NA_real_
  c(mean = m, lower = m - half, upper = m + half, R = n)
}
countRow <- function(cell, class, side, metric, hits, n) {
  ci <- if (n > 0L) stats::binom.test(hits, n)$conf.int else c(NA_real_, NA_real_)
  data.frame(cell = cell, class = class, side = side, metric = metric, mean = hits, lower = ci[[1L]], upper = ci[[2L]], R = n,
    stringsAsFactors = FALSE)
}
cell_rows <- list()
falsifier_rows <- list()
seconds <- list()
stopped_at <- NULL
e2_rows <- rows[rows$kind == "e2", , drop = FALSE]
for (cell in intersect(cells, E2_CELLS)) {
  for (class in c("lipid", "metab")) {
    b <- e2_rows[e2_rows$side == "base" & e2_rows$cell == cell & e2_rows$class == class, , drop = FALSE]
    h <- e2_rows[e2_rows$side == "head" & e2_rows$cell == cell & e2_rows$class == class, , drop = FALSE]
    b <- b[order(b$r), , drop = FALSE]
    h <- h[order(h$r), , drop = FALSE]
    stopifnot(nrow(b) == replicates, nrow(h) == replicates, identical(b$r, h$r))
    seconds[[paste(cell, class)]] <- c(base = sum(b$seconds, na.rm = TRUE), head = sum(h$seconds, na.rm = TRUE))
    for (name in setdiff(scalar, "seconds")) {
      for (side in c("base", "head")) {
        s <- summarise((if (side == "base") b else h)[[name]])
        cell_rows[[length(cell_rows) + 1L]] <- data.frame(cell = cell, class = class, side = side, metric = name,
          mean = s[["mean"]], lower = s[["lower"]], upper = s[["upper"]], R = s[["R"]], stringsAsFactors = FALSE)
      }
      s <- summarise(h[[name]] - b[[name]])
      cell_rows[[length(cell_rows) + 1L]] <- data.frame(cell = cell, class = class, side = "head-base", metric = name,
        mean = s[["mean"]], lower = s[["lower"]], upper = s[["upper"]], R = s[["R"]], stringsAsFactors = FALSE)
    }
    for (side in c("base", "head")) {
      x <- if (side == "base") b else h
      cell_rows[[length(cell_rows) + 1L]] <- countRow(cell, class, side, "any_call", sum(x$any_call, na.rm = TRUE), nrow(x))
      cell_rows[[length(cell_rows) + 1L]] <- countRow(cell, class, side, "assay_dropped_or_refused",
        sum(x$dropped | x$refused | nzchar(x$error), na.rm = TRUE), nrow(x))
    }
    # F1 and F2 are evaluated on journeys that chose controls (a dropped assay has none).
    f1 <- if (lifecycle) sum(!b$matches_first, !h$matches_first, na.rm = TRUE) else sum(!b$matches_first, na.rm = TRUE)
    f2 <- sum(!h$matches_oracle, na.rm = TRUE)
    f3 <- sum(h$refused, na.rm = TRUE)
    n_eval <- c(f1 = sum(!is.na(b$matches_first)) + if (lifecycle) sum(!is.na(h$matches_first)) else 0L,
      f2 = sum(!is.na(h$matches_oracle)), f3 = nrow(h))
    for (i in 1:3) falsifier_rows[[length(falsifier_rows) + 1L]] <- data.frame(falsifier = paste0("F", i), cell = cell,
      class = class, count = c(f1, f2, f3)[[i]], evaluated = n_eval[[i]], stringsAsFactors = FALSE)
  }
  triggered <- Filter(function(row) identical(row$cell, cell) && row$count > 0L, falsifier_rows)
  if (stop_rule && length(triggered)) {
    stopped_at <- list(cell = cell, falsifiers = unique(vapply(triggered, function(row) row$falsifier, "")))
    break
  }
}
oc <- rows[rows$kind == "oracle_check", , drop = FALSE]
if (nrow(oc)) falsifier_rows[[length(falsifier_rows) + 1L]] <- data.frame(falsifier = "oracle_check", cell = "N-r",
  class = "lipid,metab", count = sum(is.na(oc$matches_oracle) | !oc$matches_oracle | nzchar(oc$error), na.rm = TRUE), evaluated = nrow(oc),
  stringsAsFactors = FALSE)

# E1: one row per data set, class and side, plus the head-against-base comparison.
fixture_rows <- list()
if ("E1" %in% cells && is.null(stopped_at)) {
  e1 <- rows[rows$kind == "e1", , drop = FALSE]
  bc <- controlsOf("base")
  hc <- controlsOf("head")
  bl <- logfcOf("base")
  hl <- logfcOf("head")
  for (k in unique(paste(e1$dataset, e1$class))) {
    b <- e1[e1$side == "base" & paste(e1$dataset, e1$class) == k, , drop = FALSE]
    h <- e1[e1$side == "head" & paste(e1$dataset, e1$class) == k, , drop = FALSE]
    id <- key(as.list(h[1L, ]))
    common <- intersect(names(bl[[id]]), names(hl[[id]]))
    union_size <- length(union(bc[[id]], hc[[id]]))
    extra <- data.frame(max_abs_dlogfc = if (length(common)) max(abs(hl[[id]][common] - bl[[id]][common]), na.rm = TRUE) else NA_real_,
      jaccard = if (union_size) length(intersect(bc[[id]], hc[[id]])) / union_size else NA_real_)
    fixture_rows[[length(fixture_rows) + 1L]] <- cbind(b, extra)
    fixture_rows[[length(fixture_rows) + 1L]] <- cbind(h, extra)
    f1 <- if (lifecycle) sum(!b$matches_first, !h$matches_first, na.rm = TRUE) else sum(!b$matches_first, na.rm = TRUE)
    falsifier_rows[[length(falsifier_rows) + 1L]] <- data.frame(falsifier = c("F1", "F2", "F3"), cell = "E1", class = k,
      count = c(f1, sum(!h$matches_oracle, na.rm = TRUE), sum(h$refused, na.rm = TRUE)),
      evaluated = c(if (lifecycle) 2L else 1L, 1L, 1L), stringsAsFactors = FALSE)
  }
}

writeRows <- function(x, name) if (length(x)) utils::write.table(do.call(rbind, x), file.path(out, name), sep = "\t",
  quote = FALSE, row.names = FALSE, na = "NA")
writeRows(cell_rows, "cells.tsv")
# The counts the lifecycle must show (ticket section 10 as ruled by VAL-100): F1, F3 and the oracle check 0; F2 every
# evaluated head journey in the cells with planted effects and in E1 (the tree is base, so the repair's change is
# absent), and 0 in the null cells N-r and N-s, where the designed selection equals the row-order one (the top q tie
# block holds every control). Outside the lifecycle every falsifier's stated count is 0.
if (length(falsifier_rows)) {
  falsifier_rows <- lapply(falsifier_rows, function(row) {
    row$stated <- ifelse(lifecycle & !fault & row$falsifier == "F2" & !row$cell %in% c("N-r", "N-s"), row$evaluated, 0L)
    row
  })
}
lifecycle_counts_as_stated <- if (length(falsifier_rows)) all(vapply(falsifier_rows, function(row) all(row$count == row$stated),
  logical(1L))) else NA
writeRows(falsifier_rows, "falsifiers.tsv")
writeRows(fixture_rows, "fixtures.tsv")
journey_seconds <- rows$seconds[is.finite(rows$seconds) & rows$kind != "oracle_check"]
full_journeys <- 46L + length(E2_CELLS) * 2L * 50L
provenance <- c(list(ticket = "KIT-190 E", head_commit = head_commit, base_commit = base_commit, base_arg = args[[2L]],
  generator = GENERATOR_SOURCE, lifecycle = lifecycle, fault = fault, stop_rule = stop_rule, cells = cells,
  replicates = replicates, workers_per_side = workers, falsifier_counts_as_stated = lifecycle_counts_as_stated,
  seed_rule = "set.seed(190000000 + 10000 * cell + r); cells N-r 1, N-s 2, E-r 3, E-s 4, E-p 5, E-g 6",
  side_exit_status = as.list(stats::setNames(status, names(sides))), journeys = nrow(rows),
  mean_journey_seconds = if (length(journey_seconds)) mean(journey_seconds) else NA_real_,
  projected_minutes_per_side_full = if (length(journey_seconds)) mean(journey_seconds) * full_journeys / workers / 60 else NA_real_,
  seconds_per_cell_class = lapply(seconds, as.list), wall_seconds_sides = wall_sides,
  wall_seconds = as.numeric(difftime(Sys.time(), started_at, units = "secs")), stopped_at = stopped_at,
  started = format(started_at, "%Y-%m-%dT%H:%M:%S%z"), finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  environment_facts())
jsonlite::write_json(provenance, file.path(out, "provenance.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat(sprintf("KIT-190 E written: %d journeys, %d cell rows, %d falsifier rows, wall %.1f s\n", nrow(rows), length(cell_rows),
  length(falsifier_rows), provenance$wall_seconds))
if (!is.null(stopped_at))
  k190Abort("k190_falsifier", sprintf("falsifier %s triggered in cell %s", paste(stopped_at$falsifiers, collapse = ","),
    stopped_at$cell))
