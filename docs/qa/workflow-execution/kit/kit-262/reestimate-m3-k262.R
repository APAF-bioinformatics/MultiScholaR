# KIT-262 M3: the one blinded re-estimation of R for the two stress cells (ticket 6.1; VAL-185's form, reused). A copy of KIT-190's
# `reestimate.R` (bound-declaration.md section 8) with marked lines [R]: the arms are C1 (the candidate) and P (the reference), the
# re-estimated cells are S-N and S-E (k_cells = 2: the per-cell target is 0.80^(1/2) at one-sided alpha 0.05); N-r and N-s are NOT
# re-estimated (C2 is not tested, VAL-302 O-1 (ii)). The floors are the ticket's minima (S-N 250, S-E 200), at most 1,000, multiples of
# 50. THE CEILING LOOP OF KIT-190 IS REMOVED [R]: R is never reduced. If any R exceeds 400 (the ticket's "R doubled") the plan says so
# and the projected cost goes to the coordinator (G2); nothing in this script lowers the target to fit a budget. Evidence code only.
#
# Usage (run-tests scripts lane; base R only):
#   Rscript reestimate-m3-k262.R estimate <capped rows .rds> <out dir>
#       The stage-1 rows of the two stress cells (columns cell, arm in {P, C1}, r, fdp, any_call), private to the runner. Writes
#       <out dir>/reestimate.tsv with exactly the columns cell, endpoint, n_pairs, value: effect cell: fdp_difference_sd =
#       sd(FDP(C1) - FDP(P)) over the paired seeds; null cell: any_call_discordant_share = (b + c') / n. Nothing else is computed.
#   Rscript reestimate-m3-k262.R plan <reestimate.tsv> <provenance.json> <spent core-hours> <out dir>
#       The fixed formula; writes plan.tsv and plan-summary.tsv. provenance.json carries seconds_median_per_unit (stage 1's).
#   Rscript reestimate-m3-k262.R alter <reestimate.tsv> <copy>
#       Lifecycle: a copy with the first row's value altered by +0.01 (stage 2's byte-for-byte check must then stop).
#   Rscript reestimate-m3-k262.R selftest <out dir>
#       Planted inputs give byte-identical estimate files (a +0.30 shift of every C1 FDP in an effect cell; (b, c') = (10, 0), (0, 10),
#       (5, 5) on 50 pairs in a null cell) and the plan's R values equal the hand values (e5hand-m3-k262.py). Exit status 1 on a failure.
EFFECT_CELLS <- c("S-E")  # [R]
NULL_CELLS <- c("S-N")  # [R]
SCHEMA <- c("cell", "endpoint", "n_pairs", "value")
FLOOR <- c(`S-N` = 250, `S-E` = 200)  # [R] the ticket's minima for the stress cells
LIMIT <- 1000  # [R]
ROUTE_ABOVE <- 400  # [R] above this the projected cost goes to the coordinator
SPREAD <- 1.93  # [R] runner spread: the largest per-sweep ratio to the lifecycle seen in M2's counted jobs (audit M11); KIT-190 used 1.25

estimate <- function(capped) {
  out <- list()
  for (cell in intersect(c(NULL_CELLS, EFFECT_CELLS), unique(capped$cell))) {
    f <- capped[capped$cell == cell & capped$arm == "P", , drop = FALSE]  # [R]
    d <- capped[capped$cell == cell & capped$arm == "C1", , drop = FALSE]  # [R]
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

writeEstimate <- function(res, path) {
  # Each value printed alone, unpadded, so a job covering some cells writes the same lines as one covering all.
  res$value <- sprintf("%.10g", res$value)
  utils::write.table(res, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

# [R] plan(): KIT-190's formula with the ticket's cells and floors, k_cells = 2, and NO ceiling loop (R is never reduced).
plan <- function(est, prov, spent) {
  z <- stats::qnorm(0.95)
  P <- function(cell, R) {
    v <- est$value[est$cell == cell]
    if (!length(v) || !is.finite(v) || v <= 0) return(1)
    if (cell %in% EFFECT_CELLS) stats::pnorm(0.05 * sqrt(R) / v - z) else stats::pnorm(0.10 * sqrt(R) / sqrt(v) - z)
  }
  secs <- as.numeric(prov$seconds_median_per_unit)
  if (!length(secs) || !is.finite(secs)) stop("plan: provenance.json carries no seconds_median_per_unit")
  pick <- function(t) vapply(names(FLOOR), function(c) {
    for (R in seq(FLOOR[[c]], LIMIT, by = 50)) if (P(c, R) >= t) return(R)
    LIMIT
  }, numeric(1L))
  t <- 0.80^(1 / length(FLOOR))
  R <- pick(t)
  core <- function(R) spent + sum(R) * secs / 3600 * SPREAD
  achieved <- prod(vapply(names(R), function(c) P(c, R[[c]]), numeric(1L)))
  above <- any(R > ROUTE_ABOVE)
  list(plan = data.frame(cell = names(R), R = as.integer(R), stringsAsFactors = FALSE),
    summary = data.frame(target_per_cell = t, k_cells = length(FLOOR), achieved_joint_probability = achieved, projected_core_hours = core(R),
      spent_core_hours = spent, above_400 = above, route = if (above) "projected cost to the coordinator (G2); R is never reduced" else "none",
      stringsAsFactors = FALSE))
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
} else if (identical(mode, "alter")) {
  x <- readLines(args[[2L]])
  fields <- strsplit(x[[2L]], "\t", fixed = TRUE)[[1L]]
  fields[[4L]] <- format(as.numeric(fields[[4L]]) + 0.01, digits = 10)
  x[[2L]] <- paste(fields, collapse = "\t")
  dir.create(dirname(args[[3L]]), recursive = TRUE, showWarnings = FALSE)
  writeLines(x, args[[3L]])
} else if (identical(mode, "selftest")) {
  out <- args[[2L]]
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  set.seed(262888)
  base_rows <- function() {
    f <- stats::runif(50)
    d <- pmin(1, f + stats::rnorm(50, 0, 0.1))
    rbind(data.frame(cell = "S-E", arm = "P", r = 1:50, fdp = f, any_call = NA),
      data.frame(cell = "S-E", arm = "C1", r = 1:50, fdp = d, any_call = NA))
  }
  nulls <- function(b, c) {
    a_d <- c(rep(TRUE, b), rep(FALSE, c), rep(FALSE, 50 - b - c))
    a_f <- c(rep(FALSE, b), rep(TRUE, c), rep(FALSE, 50 - b - c))
    rbind(data.frame(cell = "S-N", arm = "P", r = 1:50, fdp = NA, any_call = a_f),
      data.frame(cell = "S-N", arm = "C1", r = 1:50, fdp = NA, any_call = a_d))
  }
  e0 <- base_rows()
  e1 <- e0
  e1$fdp[e1$arm == "C1"] <- e1$fdp[e1$arm == "C1"] + 0.30
  files <- list(effect_plain = e0, effect_shifted = e1, null_10_0 = nulls(10, 0), null_0_10 = nulls(0, 10), null_5_5 = nulls(5, 5))
  for (nm in names(files)) writeEstimate(estimate(files[[nm]]), file.path(out, paste0("selftest-", nm, ".tsv")))
  same <- function(a, b) identical(readBin(file.path(out, paste0("selftest-", a, ".tsv")), "raw", 1e6),
    readBin(file.path(out, paste0("selftest-", b, ".tsv")), "raw", 1e6))
  checks <- c(effect_shift = same("effect_plain", "effect_shifted"), null_swap = same("null_10_0", "null_0_10"),
    null_split = same("null_10_0", "null_5_5"))
  # the plan against the hand values of e5hand-m3-k262.py: (SD 0.265, discordant share 0.30) -> R 250 and 300; (SD 0.35) -> 450;
  # (SD 0.40) -> 550 (above 400: routed); (SD 0.10, share 0.01) -> the floors 200 and 250; spent 10 and 0.5 s per unit, spread 1.93
  planCase <- function(sd_value, share) plan(data.frame(cell = c("S-E", "S-N"), endpoint = c("fdp_difference_sd", "any_call_discordant_share"),
    n_pairs = 50L, value = c(sd_value, share)), list(seconds_median_per_unit = 0.5), 10)
  a <- planCase(0.265, 0.30); b <- planCase(0.35, 0.30); c4 <- planCase(0.40, 0.30); lo <- planCase(0.10, 0.01)
  r_of <- function(x, cell) x$plan$R[x$plan$cell == cell]
  checks <- c(checks, plan_265_300 = r_of(a, "S-E") == 250L && r_of(a, "S-N") == 300L && !a$summary$above_400,
    plan_035 = r_of(b, "S-E") == 450L && b$summary$above_400, plan_040 = r_of(c4, "S-E") == 550L && c4$summary$above_400,
    plan_floors = r_of(lo, "S-E") == 200L && r_of(lo, "S-N") == 250L && !lo$summary$above_400,
    plan_cost = isTRUE(all.equal(a$summary$projected_core_hours, 10 + (250 + 300) * 0.5 / 3600 * 1.93, tolerance = 1e-12)),
    plan_never_reduced = r_of(c4, "S-E") >= 550L)
  for (nm in names(checks)) cat(sprintf("selftest %s: %s\n", nm, if (checks[[nm]]) "ok" else "FAILED"))
  if (!all(checks)) quit(status = 1L)
} else stop("usage: reestimate-m3-k262.R estimate|plan|alter|selftest ...")
