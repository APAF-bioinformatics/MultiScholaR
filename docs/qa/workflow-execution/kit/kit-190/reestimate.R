# KIT-190 bound-declaration.md section 8: the one blinded re-estimation of R (VAL-185). Evidence code: no package code.
#
# Usage (run-tests scripts lane; base R only):
#   Rscript reestimate.R estimate <capped rows .rds> <out dir>
#       The capped arms' stage-1 rows (columns cell, arm in {F+B, D+B}, r, fdp, any_call), private to the runner.
#       Writes <out dir>/reestimate.tsv with exactly the columns cell, endpoint, n_pairs, value:
#       effect cells: fdp_difference_sd = sd(FDP(D+B) - FDP(F+B)) over the paired seeds (shift invariant);
#       null cells: any_call_discordant_share = (b + c') / n (invariant to exchanging b and c').
#       Nothing else is computed, printed or saved: no mean difference, test statistic or per-arm rate.
#   Rscript reestimate.R plan <reestimate.tsv> <provenance.json> <spent core-hours> <out dir>
#       The fixed formula: per-cell R (multiples of 50, option C's floor, at most 1,000), the target used and the
#       achieved joint probability under equality, the projected core-hours; written to <out dir>/plan.tsv and
#       <out dir>/plan-summary.tsv. If P_c(1,000) < t then R_c = 1,000.
#   Rscript reestimate.R selftest <out dir>
#       Lifecycle item (viii): planted inputs give byte-identical estimate files (a +0.30 shift of every D+B FDP in
#       an effect cell; (b, c') = (10, 0), (0, 10), (5, 5) on 50 pairs in a null cell). Exit status 1 on a failure.
EFFECT_CELLS <- c("E-r", "E-s", "E-p", "E-g", "E-w")
NULL_CELLS <- c("N-r", "N-s", "N-w")
SCHEMA <- c("cell", "endpoint", "n_pairs", "value")
FLOOR <- c(`N-r` = 100, `N-s` = 100, `N-w` = 250, `E-r` = 200, `E-s` = 200, `E-p` = 200, `E-g` = 200, `E-w` = 200)
CEILING <- 100

estimate <- function(capped) {
  out <- list()
  for (cell in intersect(c(NULL_CELLS, EFFECT_CELLS), unique(capped$cell))) {
    f <- capped[capped$cell == cell & capped$arm == "F+B", , drop = FALSE]
    d <- capped[capped$cell == cell & capped$arm == "D+B", , drop = FALSE]
    r <- intersect(f$r, d$r)
    f <- f[match(r, f$r), , drop = FALSE]
    d <- d[match(r, d$r), , drop = FALSE]
    if (cell %in% EFFECT_CELLS) {
      keep <- is.finite(f$fdp) & is.finite(d$fdp)
      diffs <- d$fdp[keep] - f$fdp[keep]
      out[[cell]] <- data.frame(cell = cell, endpoint = "fdp_difference_sd", n_pairs = length(diffs),
        value = if (length(diffs) > 1L) stats::sd(diffs) else NA_real_, stringsAsFactors = FALSE)
    } else {
      keep <- !is.na(f$any_call) & !is.na(d$any_call)
      a <- as.logical(d$any_call[keep])
      b <- as.logical(f$any_call[keep])
      out[[cell]] <- data.frame(cell = cell, endpoint = "any_call_discordant_share", n_pairs = sum(keep),
        value = if (sum(keep)) sum(a != b) / sum(keep) else NA_real_, stringsAsFactors = FALSE)
    }
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  if (!identical(names(res), SCHEMA) || !all(res$endpoint %in% c("fdp_difference_sd", "any_call_discordant_share")))
    stop("reestimate: output schema violated")
  res
}

writeEstimate <- function(res, path) utils::write.table(format(res, digits = 10), path, sep = "\t", quote = FALSE,
  row.names = FALSE)

plan <- function(est, prov, spent) {
  z <- stats::qnorm(0.95)
  P <- function(cell, R) {
    v <- est$value[est$cell == cell]
    if (!length(v) || !is.finite(v) || v <= 0) return(1)
    if (cell %in% EFFECT_CELLS) stats::pnorm(0.05 * sqrt(R) / v - z) else stats::pnorm(0.10 * sqrt(R) / sqrt(v) - z)
  }
  secs <- unlist(prov$seconds_median_per_journey)
  secs[!is.finite(secs)] <- 0
  cost <- function(R) sum(vapply(names(R), function(c) {
    if (c %in% c("N-w", "E-w")) R[[c]] * (3 * secs[["bounded"]] + secs[["B"]] + secs[["F"]] + secs[["D"]]) else
      2 * secs[["bounded"]] * R[[c]] + (secs[["F"]] + secs[["D"]]) * max(R[[c]] - 50, 0)
  }, numeric(1L))) + 42 * 2 * secs[["bounded"]]
  pick <- function(t) vapply(names(FLOOR), function(c) {
    for (R in seq(FLOOR[[c]], 1000, by = 50)) if (P(c, R) >= t) return(R)
    1000
  }, numeric(1L))
  t <- 0.80^(1 / 8)
  R <- pick(t)
  core <- function(R) spent + cost(R) / 3600 * 1.25
  while (core(R) > CEILING && t > 0) {
    t <- t - 0.0005
    R <- pick(t)
    if (all(R == FLOOR)) break
  }
  achieved <- prod(vapply(names(R), function(c) P(c, R[[c]]), numeric(1L)))
  list(plan = data.frame(cell = names(R), R = as.integer(R), stringsAsFactors = FALSE),
    summary = data.frame(target_per_cell = t, achieved_joint_probability = achieved, projected_core_hours = core(R),
      spent_core_hours = spent, ceiling = CEILING))
}

args <- commandArgs(TRUE)
mode <- args[[1L]]
if (identical(mode, "estimate")) {
  res <- estimate(readRDS(args[[2L]]))
  dir.create(args[[3L]], recursive = TRUE, showWarnings = FALSE)
  writeEstimate(res, file.path(args[[3L]], "reestimate.tsv"))
} else if (identical(mode, "plan")) {
  est <- utils::read.delim(args[[2L]], stringsAsFactors = FALSE)
  prov <- jsonlite::read_json(args[[3L]])
  res <- plan(est, prov, as.numeric(args[[4L]]))
  dir.create(args[[5L]], recursive = TRUE, showWarnings = FALSE)
  utils::write.table(res$plan, file.path(args[[5L]], "plan.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(res$summary, file.path(args[[5L]], "plan-summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  print(res$plan)
  print(res$summary)
} else if (identical(mode, "selftest")) {
  out <- args[[2L]]
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  set.seed(190888)
  base_rows <- function() {
    f <- stats::runif(50)
    d <- pmin(1, f + stats::rnorm(50, 0, 0.1))
    rbind(data.frame(cell = "E-r", arm = "F+B", r = 1:50, fdp = f, any_call = NA),
      data.frame(cell = "E-r", arm = "D+B", r = 1:50, fdp = d, any_call = NA))
  }
  nulls <- function(b, c) {
    a_d <- c(rep(TRUE, b), rep(FALSE, c), rep(FALSE, 50 - b - c))
    a_f <- c(rep(FALSE, b), rep(TRUE, c), rep(FALSE, 50 - b - c))
    rbind(data.frame(cell = "N-r", arm = "F+B", r = 1:50, fdp = NA, any_call = a_f),
      data.frame(cell = "N-r", arm = "D+B", r = 1:50, fdp = NA, any_call = a_d))
  }
  e0 <- base_rows()
  e1 <- e0
  e1$fdp[e1$arm == "D+B"] <- e1$fdp[e1$arm == "D+B"] + 0.30
  files <- list(effect_plain = e0, effect_shifted = e1, null_10_0 = nulls(10, 0), null_0_10 = nulls(0, 10),
    null_5_5 = nulls(5, 5))
  for (nm in names(files)) writeEstimate(estimate(files[[nm]]), file.path(out, paste0("selftest-", nm, ".tsv")))
  same <- function(a, b) identical(readBin(file.path(out, paste0("selftest-", a, ".tsv")), "raw", 1e6),
    readBin(file.path(out, paste0("selftest-", b, ".tsv")), "raw", 1e6))
  checks <- c(effect_shift = same("effect_plain", "effect_shifted"), null_swap = same("null_10_0", "null_0_10"),
    null_split = same("null_10_0", "null_5_5"))
  for (nm in names(checks)) cat(sprintf("selftest %s: %s\n", nm, if (checks[[nm]]) "identical" else "DIFFERENT"))
  if (!all(checks)) quit(status = 1L)
} else stop("usage: reestimate.R estimate|plan|selftest ...")
