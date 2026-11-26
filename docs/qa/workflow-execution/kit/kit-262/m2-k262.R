# KIT-262 M2 (ticket section 5): how often the lipidomics/metabolomics automatic RUV search ends with no usable
# percentage, under the present rule P, under C1 (R5 judged at min(k, cap)) and, for coverage only (VAL-302, O-1 (ii)),
# C2 (proteomics' truncated curve); the cause of each R5 rejection; how often a candidate changes P's choice. Simulated
# data, no differential analysis, no normalisation between samples. Evidence code: no package change.
#
# Usage, from the package root on a run-tests runner:
#   Rscript docs/qa/workflow-execution/kit/kit-262/m2-k262.R <out dir relative to the package root>
# Environment:
#   K262_CELLS   cells to run, e.g. "1:36" or "3,25,50" (default: all 198)
#   K262_REPS    replicates, e.g. "1:40" (default) or "91" (the lifecycle)
#   K262_E8      "1" (default in lifecycle and counted runs): the ten E8 cells at replicates <= K262_E8_MAX also run
#                as metabolomics and compared with lipidomics (E8)
#   K262_E8_MAX  highest replicate of the E8 comparison (default 5)
#   K262_E1P     replicates at which the unhooked sweep is also run and compared (E1'; default "7,27")
#   K262_WORKERS forked workers (default 4)
#   K262_LIFECYCLE "1": row files under tempdir(); only checks.tsv, timing.tsv and hashes.tsv are written to the out dir
#   K262_FAULT   "1": the first task of the first cell signals an error after its data are generated (L2's failure
#                check: the error row must be written, never dropped, and never counted as a result)
suppressPackageStartupMessages(pkgload::load_all(".", export_all = TRUE, helpers = FALSE, quiet = TRUE))
logger::log_threshold(logger::FATAL)
args <- commandArgs(trailingOnly = TRUE)
out <- normalizePath(file.path(getwd(), args[[1L]]), mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source("docs/qa/workflow-execution/kit/kit-262/rules-k262.R")

# ---- copy: tests/testthat/test-workflow-omics-ruv-cancor-once.R:43-54, verbatim ----
ruvComparable <- function(result) {
  lapply(result, function(assay) {
    lapply(assay, function(value) {
      if (inherits(value, "ggplot")) return(value$data)
      if (is.list(value) && !is.data.frame(value)) {
        return(lapply(value, function(inner) if (inherits(inner, "ggplot")) inner$data else inner))
      }
      value
    })
  })
}
# ---- end of the copy ----

# ---- copy: docs/qa/workflow-execution/kit/kit-190/effect.R e2PickCells (helper-kit172.R:21-28 at d4e2eb2c), verbatim ----
e2PickCells <- function(allowed, groups, per_group, count) {
  cells <- which(allowed, arr.ind = TRUE)
  cells <- cells[sample.int(nrow(cells)), , drop = FALSE]
  key <- paste(cells[, 1L], groups[cells[, 2L]])
  rank <- stats::ave(seq_along(key), key, FUN = seq_along)
  cells <- cells[rank <= per_group, , drop = FALSE]
  cells[seq_len(min(count, nrow(cells))), , drop = FALSE]
}
# ---- end of the copy ----

# ---- copy: kit-196 sim-k188-69fd01bc.R:252-254 (KIT-188's mask solver), verbatim ----
solveLod <- function(values, rate) {
    stats::uniroot(function(lod) mean(stats::plogis(lod - values)) - rate, c(min(values) - 60, max(values) + 60), tol = 1e-10)$root
}
# ---- end of the copy ----

# ---- the grid (ticket section 5; O-4 (a)): cells 1 to 180 in the order n, g, missingness, F; cells 181 to 198 n = 4 ----
M2_MISSING <- data.frame(id = c("none", "MCAR1", "MCAR5", "MCAR10", "MCAR20", "AD1", "AD5", "AD10", "AD20"),
  mech = c("none", rep("MCAR", 4L), rep("AD", 4L)), rate = c(0, 0.01, 0.05, 0.10, 0.20, 0.01, 0.05, 0.10, 0.20),
  stringsAsFactors = FALSE)
m2Grid <- function() {
  main <- expand.grid(f = c(300L, 1000L), m = seq_len(nrow(M2_MISSING)), g = c(2L, 3L), n = c(6L, 12L, 24L, 48L, 96L))
  main <- main[order(main$n, main$g, main$m, main$f), c("n", "g", "m", "f")]
  small <- expand.grid(f = c(300L, 1000L), m = seq_len(nrow(M2_MISSING)))
  small <- small[order(small$m, small$f), ]
  small <- data.frame(n = 4L, g = 2L, m = small$m, f = small$f)
  grid <- rbind(main, small)
  grid$cell <- seq_len(nrow(grid))
  grid$missing_id <- M2_MISSING$id[grid$m]
  grid$mech <- M2_MISSING$mech[grid$m]
  grid$rate <- M2_MISSING$rate[grid$m]
  grid$report_only <- grid$cell > 180L
  rownames(grid) <- NULL
  grid[c("cell", "n", "g", "missing_id", "mech", "rate", "f", "report_only")]
}
E8_CELLS <- c(3L, 25L, 50L, 62L, 83L, 100L, 121L, 143L, 150L, 172L)

# ---- the generator: a declared generalisation of KIT-190's e2Values (effect.R, helper-kit172.R:34-62 at d4e2eb2c) ----
# Marked lines: [A] n samples in g equal groups of the design, no IS_ rows; [B] the effect in the LAST group, planted
# always; [C] the random factor pattern for any n (sample variance n / (n - 1)); [D] the missingness mechanism, applied
# after the 2 % values at 1; [E] no internal-standard rows.
m2Design <- function(n, g) {
  runs <- sprintf("S%02d", seq_len(n))
  group <- rep(LETTERS[seq_len(g)], each = n / g)
  data.frame(Run = runs, group = group, batch = "B1", replicates = stats::ave(seq_len(n), group, FUN = seq_along),
    tech_rep_group = runs, stringsAsFactors = FALSE)
}

m2Values <- function(seed, n, g, n_features, mech, rate, floor_order = "before") {  # [A] [G] the order option
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")  # [F] the RNG kinds named
  design <- m2Design(n, g)                                                                                    # [A]
  groups <- design$group
  level <- stats::runif(n_features, log2(0.01), log2(1e6))
  logs <- matrix(level, n_features, nrow(design)) + matrix(stats::rnorm(n_features * nrow(design), 0, 0.25), n_features)
  effects <- rep(0, n_features)
  loadings <- rep(0, n_features)
  changed <- sample.int(n_features, round(0.1 * n_features))  # [B]
  effects[changed] <- sample(c(-1, 1), length(changed), replace = TRUE)
  logs[, groups == LETTERS[g]] <- logs[, groups == LETTERS[g]] + effects                                      # [B]
  loadings <- stats::rnorm(n_features, 0, 0.5)
  w <- stats::rnorm(nrow(design))
  pattern <- (w - mean(w)) / stats::sd(w) * sqrt(n / (n - 1))                                                  # [C]
  logs <- logs + outer(loadings, pattern)
  values <- 2^logs
  ones_first <- function() {  # [G]
    ones <- e2PickCells(!is.na(values), groups, 2L, round(0.02 * sum(!is.na(values))))  # [G]
    values[ones] <<- 1  # [G]
  }  # [G]
  if (identical(floor_order, "before")) ones_first()  # [G] M2's order: floor values first, then the mask
  expected_share <- 0  # [D]
  if (identical(mech, "MCAR")) {  # [D]
    values[sample.int(length(values), floor(rate * length(values)))] <- NA_real_  # [D]
    expected_share <- floor(rate * length(values)) / length(values)  # [D]
  } else if (identical(mech, "AD")) {  # [D]
    y <- log2(values)  # [D]
    lod <- solveLod(y, rate)  # [D]
    expected_share <- mean(stats::plogis(lod - y))  # [D]
    u <- matrix(stats::runif(length(values)), n_features)  # [D]
    values[u < stats::plogis(lod - y)] <- NA_real_  # [D]
  }  # [D]
  if (identical(floor_order, "after")) ones_first()  # [G] KIT-172's order: the mask first, then 2 % of the observed values set to 1
  ids <- sprintf("F%04d", seq_len(n_features))                                                                 # [E]
  dimnames(values) <- list(ids, design$Run)  # [E]
  list(values = values, design = design, effects = stats::setNames(effects, ids), expected_share = expected_share)  # [E]
}

# The S4 object of a values matrix for one class (one assay, `LCMS_Pos`), as kit172Object builds it.
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

# The application's QC filter at its defaults (percentile 1, proportion 0.5, 2 values per group in at least 1 group:
# the binding's `filter_intensity`, the S4 method the reference step calls), the log2 step (offset 0), then the
# package's sweep through the binding's `optimise_ruv`, as the RUV step's kernel calls it.
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

M2_SUMMARY_COLUMNS <- c("cell", "r", "seed", "floor_order", "omic", "n", "g", "f", "missing_id", "report_only", "cap", "error",
  "features_after_filter", "complete_features", "share_complete", "realised_missing_share", "expected_missing_share",
  "min_values_per_feature", "n_ok", "n_r3", "n_r4", "n_r5", "n_r6", "n_error_status", "n_r5_cap_side", "n_r5_completeness",
  "p_ruv", "p_percentage", "p_k_sweep", "p_k_used", "p_resid_df", "c1_ruv", "c1_percentage", "c1_k_used", "c1_resid_df",
  "c2_ruv", "c2_percentage", "c2_k_used", "c2_resid_df", "c1_differs_from_p", "c2_differs_from_p", "c2_differs_from_c1",
  "p_no_ruv_cause", "e1", "e1prime", "e2", "e3", "e4", "e7", "e8", "finite_at_cap", "seconds_sweep", "seconds_total")

m2EmptyRow <- function(task) {
  row <- stats::setNames(rep(list(NA), length(M2_SUMMARY_COLUMNS)), M2_SUMMARY_COLUMNS)
  row[c("cell", "r", "seed", "floor_order", "omic", "n", "g", "f", "missing_id", "report_only")] <- list(task$cell, task$r, task$seed, task$floor_order, task$omic,
    task$n, task$g, task$f, task$missing_id, task$report_only)
  row$error <- ""
  row
}

m2Dataset <- function(task) {
  started <- proc.time()[[3L]]
  row <- m2EmptyRow(task)
  rows <- NULL
  outcome <- tryCatch({
    sim <- m2Values(task$seed, task$n, task$g, task$f, task$mech, task$rate, task$floor_order)
    if (identical(Sys.getenv("K262_FAULT"), "1") && isTRUE(task$fault)) stop("K262 injected fault (the lifecycle's failure check)")
    realised <- mean(is.na(sim$values))
    row$realised_missing_share <- realised
    row$expected_missing_share <- sim$expected_share
    prepared <- m2Prepare(m2Object(sim$values, sim$design, task$omic), task$omic)
    assay <- methods::slot(prepared$filtered, prepared$binding$data_slot)[[1L]]
    runs <- sim$design$Run
    matrix_after <- as.matrix(as.data.frame(assay)[, runs, drop = FALSE])
    row$features_after_filter <- nrow(matrix_after)
    row$complete_features <- sum(rowSums(is.na(matrix_after)) == 0L)
    row$share_complete <- row$complete_features / max(nrow(matrix_after), 1L)
    row$min_values_per_feature <- if (nrow(matrix_after)) min(rowSums(!is.na(matrix_after))) else NA_integer_
    # E4: MCAR within 0.5 points of r (exact up to the floor); AD: the solved lod's expected share within 0.5 points of
    # r (the realised share of a random mask is reported, not bounded); every feature entering the sweep has >= 2 values
    e4_share <- if (identical(task$mech, "none")) realised == 0 else if (identical(task$mech, "MCAR")) abs(realised - task$rate) <= 0.005 else
      abs(sim$expected_share - task$rate) <= 0.005
    row$e4 <- isTRUE(e4_share) && isTRUE(row$min_values_per_feature >= 2L)
    swept_at <- proc.time()[[3L]]
    captured <- k262WithCurveCapture(function() m2Sweep(prepared))
    row$seconds_sweep <- proc.time()[[3L]] - swept_at
    result <- captured$value
    res <- result[[1L]]
    trace <- res$optimization_results
    cap <- MultiScholaR:::calculateAdaptiveMaxK(trace$sample_size[[1L]])
    row$cap <- cap
    aligned <- k262AlignCurves(result, captured$curves)
    complete <- rep(NA_integer_, nrow(trace))
    parsed <- lapply(trace$error_reason, k262ParseReason)
    is_r5 <- trace$status == K262_R5
    complete[is_r5] <- vapply(parsed[is_r5], function(p) p$c, integer(1L))
    ru <- k262Rules(trace, aligned$curves[[1L]], complete, cap)
    cause <- k262Cause(ru$rows, cap)
    # E1: P's usable set, composites and winner against the shipped trace and result
    shipped_usable <- trace$status == "ok"
    same_winner <- if (isTRUE(res$success)) identical(ru$winner$P$percentage, res$best_percentage) &&
      identical(trace$best_k[[ru$winner$P$row]], as.integer(res$best_k)) else is.na(ru$winner$P$row)
    row$e1 <- identical(ru$rows$p_usable, shipped_usable) && same_winner &&
      isTRUE(all.equal(ru$rows$composite_recomputed[shipped_usable], trace$composite_score[shipped_usable], tolerance = 0))
    row$e2 <- all(vapply(parsed[is_r5], function(p) p$ok, logical(1L))) && aligned$matched
    reached <- !trace$status %in% K262_NOT_REACHED
    k_rec <- vapply(which(reached), function(i) as.integer(findBestKElbow(list(data = aligned$curves[[1L]][[i]]))), integer(1L))
    row$e3 <- identical(captured$record_failures, 0L) && aligned$matched && identical(k_rec, trace$best_k[reached])
    # E7: P within C1 within C2 and k_b <= min(k, cap), where the truncated curve has finite rows
    finite <- ru$rows$finite_at_cap
    sane <- is.na(finite) | finite
    row$finite_at_cap <- all(sane)
    row$e7 <- all((!ru$rows$p_usable | ru$rows$c1_usable)[sane]) && all((!ru$rows$c1_usable | ru$rows$c2_usable)[sane]) &&
      all((is.na(ru$rows$k_b) | is.na(ru$rows$k) | ru$rows$k_b <= pmin(ru$rows$k, cap))[sane])
    # E1': a sweep without the hook gives identical results (a declared sample of replicates)
    if (task$r %in% task$e1p) {
      row$e1prime <- identical(ruvComparable(result), ruvComparable(m2Sweep(prepared)))
    }
    status <- trace$status
    row$n_ok <- sum(status == "ok"); row$n_r3 <- sum(status == "skipped_insufficient_controls")
    row$n_r4 <- sum(status == "invalid_cancor_plot"); row$n_r5 <- sum(is_r5); row$n_r6 <- sum(status == "error_scoring")
    row$n_error_status <- sum(status %in% c("error_neg_ctrl_selection", "error_cancor"))
    row$n_r5_cap_side <- sum(cause == "cap_side", na.rm = TRUE)
    row$n_r5_completeness <- sum(cause == "completeness", na.rm = TRUE)
    for (nm in c("P", "C1", "C2")) {
      w <- ru$winner[[nm]]
      key <- tolower(nm)
      row[[paste0(key, "_ruv")]] <- !is.na(w$row)
      row[[paste0(key, "_percentage")]] <- w$percentage
      row[[paste0(key, "_k_used")]] <- w$k_used
      row[[paste0(key, "_resid_df")]] <- if (is.na(w$row)) NA_integer_ else as.integer(task$n - task$g - w$k_used)
    }
    row$p_k_sweep <- if (is.na(ru$winner$P$row)) NA_integer_ else trace$best_k[[ru$winner$P$row]]
    same_choice <- function(a, b) is.na(a$row) && is.na(b$row) || (!is.na(a$row) && !is.na(b$row) && a$percentage == b$percentage)
    row$c1_differs_from_p <- !same_choice(ru$winner$P, ru$winner$C1)
    row$c2_differs_from_p <- !same_choice(ru$winner$P, ru$winner$C2)
    row$c2_differs_from_c1 <- !same_choice(ru$winner$C1, ru$winner$C2)
    row$p_no_ruv_cause <- if (!is.na(ru$winner$P$row)) "" else if (!is.na(ru$winner$C1$row)) "rescued_by_c1" else
      if (row$n_r5_cap_side + row$n_r5_completeness > 0L) "completeness_only" else "no_r5"
    if (isTRUE(task$e8)) {
      other <- if (identical(task$omic, "lipid")) "metab" else "lipid"
      prepared_o <- m2Prepare(m2Object(sim$values, sim$design, other), other)
      trace_o <- m2Sweep(prepared_o)[[1L]]$optimization_results
      row$e8 <- identical(trace_o, trace)
    }
    rows <- data.frame(cell = task$cell, r = task$r, omic = task$omic, percentage = ru$rows$percentage, status = trace$status,
      k = trace$best_k, c = complete, cause = cause, p_usable = ru$rows$p_usable, c1_usable = ru$rows$c1_usable, k_b = ru$rows$k_b,
      c2_usable = ru$rows$c2_usable, stringsAsFactors = FALSE)
    ""
  }, error = function(e) gsub("[\t\n]", " ", conditionMessage(e)))
  row$error <- outcome
  row$seconds_total <- proc.time()[[3L]] - started
  list(row = row, rows = rows)
}

parseSet <- function(text, default) {
  if (!nzchar(text)) return(default)
  unlist(lapply(strsplit(text, ",", fixed = TRUE)[[1L]], function(piece) {
    bounds <- as.integer(strsplit(piece, ":", fixed = TRUE)[[1L]])
    if (length(bounds) == 2L) seq.int(bounds[[1L]], bounds[[2L]]) else bounds
  }))
}

grid <- m2Grid()
cells <- parseSet(Sys.getenv("K262_CELLS", ""), grid$cell)
reps <- parseSet(Sys.getenv("K262_REPS", ""), 1:40)
e8_max <- as.integer(Sys.getenv("K262_E8_MAX", "5"))
e1p <- parseSet(Sys.getenv("K262_E1P", ""), c(7L, 27L))
workers <- as.integer(Sys.getenv("K262_WORKERS", "4"))
use_e8 <- !identical(Sys.getenv("K262_E8", "1"), "0")
stopifnot(all(cells %in% grid$cell))
# KIT-262 re-measurement (VAL-381): the order of the 2 % floor values and the seed base are parameters; the defaults
# reproduce M2 exactly (floor values first; seeds 262100000 + 100 x cell + r).
floor_order <- Sys.getenv("K262_FLOOR_ORDER", "before")
stopifnot(floor_order %in% c("before", "after"))
seed_base <- as.numeric(Sys.getenv("K262_SEED_BASE", "262100000"))

# The lifecycle (K262_LIFECYCLE=1) keeps its row files under tempdir(): no share, rate, k distribution leaves it (VAL-11).
data_dir <- if (identical(Sys.getenv("K262_LIFECYCLE"), "1")) tempdir() else out
summary_file <- file.path(data_dir, "summary.tsv")
rows_file <- file.path(data_dir, "rows.tsv")
first <- TRUE
fault_pending <- identical(Sys.getenv("K262_FAULT"), "1")
# Cells are run in batches of about 40 data sets so that the workers are busy whatever the replicate count (a lifecycle
# with one replicate per cell would otherwise run one data set at a time). Results are written after every batch.
cells_per_batch <- max(1L, floor(40 / length(reps)))
batches <- split(cells, ceiling(seq_along(cells) / cells_per_batch))
for (batch in batches) {
  tasks <- unlist(lapply(batch, function(cell) {
    g <- grid[grid$cell == cell, ]
    lapply(reps, function(r) list(cell = cell, r = r, seed = seed_base + 100 * cell + r, floor_order = floor_order, omic = "lipid", n = g$n, g = g$g,
      f = g$f, missing_id = g$missing_id, mech = g$mech, rate = g$rate, report_only = g$report_only,
      e8 = use_e8 && cell %in% E8_CELLS && r <= e8_max, e1p = e1p, fault = FALSE))
  }), recursive = FALSE)
  if (fault_pending) { tasks[[1L]]$fault <- TRUE; fault_pending <- FALSE }
  done <- parallel::mclapply(tasks, function(task) tryCatch(m2Dataset(task), error = function(e) list(
    row = {x <- m2EmptyRow(task); x$error <- paste("instrument:", gsub("[\t\n]", " ", conditionMessage(e))); x}, rows = NULL)),
    mc.cores = workers, mc.preschedule = FALSE)
  summary_rows <- lapply(seq_along(done), function(i) {
    d <- done[[i]]
    if (!is.list(d) || is.null(d$row)) {
      x <- m2EmptyRow(tasks[[i]]); x$error <- paste("worker did not return:", paste(as.character(d), collapse = " ")); x
    } else d$row
  })
  summary_df <- do.call(rbind, lapply(summary_rows, function(x) {
    x[lengths(x) == 0L] <- NA
    as.data.frame(x, stringsAsFactors = FALSE)
  }))
  rows_df <- do.call(rbind, lapply(done, function(d) if (is.list(d)) d$rows))
  utils::write.table(summary_df, summary_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = first, append = !first)
  if (!is.null(rows_df)) utils::write.table(rows_df, rows_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = first, append = !first)
  first <- FALSE
  cat(sprintf("cells %d to %d done: %d data sets, %d errors, %.1f s of sweeps\n", min(batch), max(batch), nrow(summary_df),
    sum(nzchar(summary_df$error)), sum(as.numeric(summary_df$seconds_sweep), na.rm = TRUE)))
}

final <- utils::read.delim(summary_file, stringsAsFactors = FALSE, colClasses = c(error = "character"))
final$error[is.na(final$error)] <- ""
check_cols <- c("e1", "e1prime", "e2", "e3", "e4", "e7", "e8")
checks <- do.call(rbind, lapply(check_cols, function(cn) {
  v <- if (cn %in% names(final)) as.character(final[[cn]]) else rep(NA_character_, nrow(final))
  ok <- !nzchar(final$error)
  data.frame(check = toupper(cn), units = sum(ok & !is.na(v) & v %in% c("TRUE", "FALSE")), passed = sum(ok & v %in% "TRUE"),
    failed = sum(ok & v %in% "FALSE"), stringsAsFactors = FALSE)
}))
# completion: every task of the plan has a row (a data set that errors still has its error row); the fault run's error row
expected_n <- length(cells) * length(reps)
checks <- rbind(checks, data.frame(check = "COMPLETION", units = expected_n, passed = nrow(final), failed = abs(expected_n - nrow(final)),
  stringsAsFactors = FALSE))
if (identical(Sys.getenv("K262_FAULT"), "1")) {
  fault_ok <- sum(grepl("K262 injected fault", final$error)) == 1L && sum(nzchar(final$error)) == 1L
  checks <- rbind(checks, data.frame(check = "FAULT_ROW", units = 1L, passed = as.integer(fault_ok), failed = as.integer(!fault_ok),
    stringsAsFactors = FALSE))
}
# E6 on the lifecycle's own rows: the scorer twice, byte-identical (the outputs stay under tempdir())
if (identical(Sys.getenv("K262_LIFECYCLE"), "1")) {
  source("docs/qa/workflow-execution/kit/kit-262/score-k262.R")
  scoreM2(final, file.path(tempdir(), "e6a"), expected_reps = reps)
  scoreM2(final, file.path(tempdir(), "e6b"), expected_reps = reps)
  files <- list.files(file.path(tempdir(), "e6a"))
  same_files <- identical(files, list.files(file.path(tempdir(), "e6b"))) && length(files) > 0L &&
    all(unname(tools::md5sum(file.path(tempdir(), "e6a", files))) == unname(tools::md5sum(file.path(tempdir(), "e6b", files))))
  checks <- rbind(checks, data.frame(check = "E6", units = length(files), passed = as.integer(same_files) * length(files),
    failed = as.integer(!same_files), stringsAsFactors = FALSE))
}
# re-measurement checks: the seeds follow the declared formula; the order column equals the order requested
seed_ok <- sum(final$seed == seed_base + 100 * final$cell + final$r)
checks <- rbind(checks, data.frame(check = "SEEDS", units = nrow(final), passed = seed_ok, failed = nrow(final) - seed_ok, stringsAsFactors = FALSE))
order_ok <- sum(final$floor_order %in% floor_order)
checks <- rbind(checks, data.frame(check = "ORDERCOL", units = nrow(final), passed = order_ok, failed = nrow(final) - order_ok, stringsAsFactors = FALSE))
# K262_REPRO=1 (with the default order and seed base): the rows reproduce M2's committed summary rows, field by field
if (identical(Sys.getenv("K262_REPRO"), "1")) {
  committed <- do.call(rbind, lapply(Sys.glob("docs/qa/workflow-execution/kit/kit-262/m2/*/summary.tsv"), function(p)
    utils::read.delim(p, stringsAsFactors = FALSE, colClasses = "character")))
  mine <- utils::read.delim(summary_file, stringsAsFactors = FALSE, colClasses = "character")
  skip <- c("seconds_sweep", "seconds_total", "floor_order")
  compared <- 0L; same <- 0L
  for (i in seq_len(nrow(mine))) {
    j <- which(committed$cell == mine$cell[[i]] & committed$r == mine$r[[i]])
    if (length(j) != 1L) next
    cols <- setdiff(names(committed), skip)
    compared <- compared + 1L
    same <- same + as.integer(identical(as.character(unlist(mine[i, cols])), as.character(unlist(committed[j, cols]))))
  }
  checks <- rbind(checks, data.frame(check = "REPRO", units = nrow(mine), passed = same, failed = nrow(mine) - same, stringsAsFactors = FALSE))
}
# K262_LIFECYCLE: the two orders give different data at one AD seed (the option does something)
if (identical(Sys.getenv("K262_LIFECYCLE"), "1") && identical(Sys.getenv("K262_FAULT"), "")) {
  g1 <- grid[grid$cell == cells[[1L]], ]
  a <- m2Values(seed_base + 100 * cells[[1L]] + reps[[1L]], g1$n, g1$g, g1$f, g1$mech, g1$rate, "before")$values
  b <- m2Values(seed_base + 100 * cells[[1L]] + reps[[1L]], g1$n, g1$g, g1$f, g1$mech, g1$rate, "after")$values
  differs <- !isTRUE(all.equal(a, b)) && identical(dim(a), dim(b))
  checks <- rbind(checks, data.frame(check = "ORDER_DIFFERS", units = 1L, passed = as.integer(differs), failed = as.integer(!differs), stringsAsFactors = FALSE))
}
checks$status <- ifelse(checks$failed == 0L, "pass", "FAIL")
utils::write.table(checks, file.path(out, "checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
timing <- stats::aggregate(cbind(seconds_sweep, seconds_total) ~ n + f, data = final[!nzchar(final$error), ], FUN = mean)
utils::write.table(timing, file.path(out, "timing.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
hashes <- data.frame(file = c("summary.tsv", "rows.tsv"), md5 = unname(tools::md5sum(c(summary_file, rows_file))))
utils::write.table(hashes, file.path(out, "hashes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(checks)
print(timing)
cat(sprintf("data sets: %d; errors: %d\n", nrow(final), sum(nzchar(final$error))))
