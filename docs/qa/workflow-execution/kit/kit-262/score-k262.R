# KIT-262 scoring of M2 (ticket section 5): the per-cell estimands with Wilson intervals, the gain rule, the skip rule
# (O-3 (2)), the stress-cell selection (O-3 (1)), the cause table, the residual-degrees-of-freedom table and the 4-sample
# report. Base R only, no random numbers: running it twice on the same input gives byte-identical files (E6). M3's
# scoring (Tango, sign-flip, bounds, outcome table) is not here: M3 is not started. Evidence code: no package change.
#
#   Rscript score-k262.R m2 <out dir> <summary.tsv> [<summary.tsv> ...]
#   Rscript score-k262.R selftest <out dir>
#
# Declared (declaration.md): gain_C1 (a) = in at least one of the 180 cells with n >= 6 the one-sided Wilson lower limit
# at z = qnorm(1 - 0.05 / 180) of the share of data sets in which C1 has a usable percentage and P has none is >= 0.05;
# gain_C2 (b) the same for C2 against P; extra_C2 (e) for C2 against C1. The skip rule: M3 is skipped iff neither (a)
# nor (b) holds. A cell with fewer than the expected data sets without error is flagged and excluded from the gain
# rule (named in the output); a data set with an error is never a denominator.
ALPHA_FAMILY <- 0.05
N_MAIN_CELLS <- 180L
Z95 <- stats::qnorm(0.975)
Z_GAIN <- stats::qnorm(1 - ALPHA_FAMILY / N_MAIN_CELLS)
GAIN_THRESHOLD <- 0.05

wilson <- function(x, n, z) {
  if (!is.finite(n) || n == 0) return(c(NA_real_, NA_real_))
  p <- x / n
  d <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / d
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  c(centre - half, centre + half)
}

readSummary <- function(paths) {
  parts <- lapply(paths, function(p) utils::read.delim(p, stringsAsFactors = FALSE, colClasses = c(error = "character")))
  common <- Reduce(intersect, lapply(parts, names))
  d <- do.call(rbind, lapply(parts, function(x) x[common]))
  d$error[is.na(d$error)] <- ""
  d
}

flag <- function(x) isTRUE(as.logical(x))

# One cell's counts and intervals. `d` the cell's data sets without error.
cellStats <- function(d) {
  n <- nrow(d)
  p <- as.logical(d$p_ruv); c1 <- as.logical(d$c1_ruv); c2 <- as.logical(d$c2_ruv)
  out <- list(n_valid = n, p_no_ruv = sum(!p), c1_no_ruv = sum(!c1), c2_no_ruv = sum(!c2),
    c1_not_p = sum(c1 & !p), c2_not_p = sum(c2 & !p), c2_not_c1 = sum(c2 & !c1),
    p_ruv = sum(p), c1_differs_among_p = sum(p & as.logical(d$c1_differs_from_p)),
    c2_differs_among_p = sum(p & as.logical(d$c2_differs_from_p)),
    nested_violations = sum((p & !c1) | (c1 & !c2)))
  for (nm in c("p_no_ruv", "c1_no_ruv", "c2_no_ruv", "c1_not_p", "c2_not_p", "c2_not_c1")) {
    w <- wilson(out[[nm]], n, Z95)
    out[[paste0(nm, "_lo")]] <- w[[1L]]; out[[paste0(nm, "_hi")]] <- w[[2L]]
    out[[paste0(nm, "_lo_gain")]] <- wilson(out[[nm]], n, Z_GAIN)[[1L]]
  }
  for (nm in c("c1_differs_among_p", "c2_differs_among_p")) {
    w <- wilson(out[[nm]], out$p_ruv, Z95)
    out[[paste0(nm, "_lo")]] <- w[[1L]]; out[[paste0(nm, "_hi")]] <- w[[2L]]
  }
  # the cause of P's no-RUV data sets, and of the R5 rows
  cause <- d$p_no_ruv_cause
  out$no_ruv_rescued_by_c1 <- sum(cause %in% "rescued_by_c1")
  out$no_ruv_completeness_only <- sum(cause %in% "completeness_only")
  out$no_ruv_no_r5 <- sum(cause %in% "no_r5")
  out$r5_rows <- sum(d$n_r5, na.rm = TRUE)
  out$r5_cap_side <- sum(d$n_r5_cap_side, na.rm = TRUE)
  out$r5_completeness <- sum(d$n_r5_completeness, na.rm = TRUE)
  out$mean_complete_features <- mean(d$complete_features, na.rm = TRUE)
  out$mean_share_complete <- mean(d$share_complete, na.rm = TRUE)
  out
}

scoreM2 <- function(d, out_dir, expected_reps = 1:40) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  d <- d[order(d$cell, d$r), ]
  cells <- sort(unique(d$cell))
  cell_rows <- list(); df_rows <- list(); stat_cache <- list()
  for (cell in cells) {
    all_rows <- d[d$cell == cell, ]
    ok <- all_rows[!nzchar(all_rows$error), ]
    expected <- length(expected_reps)
    complete <- nrow(ok) == expected && setequal(ok$r, expected_reps) && !anyDuplicated(ok$r)
    st <- cellStats(ok)
    stat_cache[[as.character(cell)]] <- st
    head <- all_rows[1L, ]
    cell_rows[[length(cell_rows) + 1L]] <- c(list(cell = cell, n = head$n, g = head$g, f = head$f, missing_id = head$missing_id,
      report_only = as.logical(head$report_only), cap = head$cap, expected = expected, n_rows = nrow(all_rows),
      n_error = sum(nzchar(all_rows$error)), cell_complete = complete), st)
    # residual degrees of freedom and k used, per cell and rule (data sets with RUV)
    for (rule in c("p", "c1", "c2")) {
      has <- as.logical(ok[[paste0(rule, "_ruv")]])
      k_used <- ok[[paste0(rule, "_k_used")]][has]
      resid <- ok[[paste0(rule, "_resid_df")]][has]
      df_rows[[length(df_rows) + 1L]] <- data.frame(cell = cell, n = head$n, g = head$g, f = head$f, missing_id = head$missing_id,
        report_only = as.logical(head$report_only), rule = toupper(rule), with_ruv = sum(has), k_at_cap = sum(k_used == head$cap),
        k_used_min = if (length(k_used)) min(k_used) else NA, k_used_max = if (length(k_used)) max(k_used) else NA,
        resid_df_median = if (length(resid)) stats::median(resid) else NA, resid_df_le1 = sum(resid <= 1), resid_df_0 = sum(resid == 0),
        stringsAsFactors = FALSE)
    }
  }
  cells_df <- do.call(rbind, lapply(cell_rows, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
  utils::write.table(cells_df, file.path(out_dir, "m2-cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(do.call(rbind, df_rows), file.path(out_dir, "m2-resid-df.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  # the cause table (the owner's question): per cell, P's no-RUV data sets by cause, and the R5 rows by cause
  causes <- cells_df[c("cell", "n", "g", "f", "missing_id", "report_only", "n_valid", "p_no_ruv", "no_ruv_rescued_by_c1",
    "no_ruv_completeness_only", "no_ruv_no_r5", "r5_rows", "r5_cap_side", "r5_completeness", "mean_complete_features", "mean_share_complete")]
  utils::write.table(causes, file.path(out_dir, "m2-causes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  # the gain rule and the skip rule: the 180 cells with n >= 6 that are complete
  main <- cells_df[!cells_df$report_only, ]
  eligible <- main[main$cell_complete, ]
  gain <- data.frame(cell = main$cell, n = main$n, g = main$g, f = main$f, missing_id = main$missing_id, cell_complete = main$cell_complete,
    n_valid = main$n_valid, c1_not_p = main$c1_not_p, c1_not_p_lo = main$c1_not_p_lo_gain, a = main$cell_complete & main$c1_not_p_lo_gain >= GAIN_THRESHOLD,
    c2_not_p = main$c2_not_p, c2_not_p_lo = main$c2_not_p_lo_gain, b = main$cell_complete & main$c2_not_p_lo_gain >= GAIN_THRESHOLD,
    c2_not_c1 = main$c2_not_c1, c2_not_c1_lo = main$c2_not_c1_lo_gain, e = main$cell_complete & main$c2_not_c1_lo_gain >= GAIN_THRESHOLD,
    stringsAsFactors = FALSE)
  utils::write.table(gain, file.path(out_dir, "m2-gain.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  a <- any(gain$a); b <- any(gain$b); e <- any(gain$e)
  skip <- !a && !b
  best <- function(col) { i <- which.max(gain[[col]]); if (!length(i)) c(cell = NA, count = NA) else c(cell = gain$cell[[i]], count = gain[[col]][[i]]) }
  verdict <- data.frame(item = c("z_gain (one-sided, 1 - 0.05 / 180)", "threshold on the lower limit", "cells in the family (n >= 6)",
    "cells complete and eligible", "incomplete or errored cells in the family", "a: C1 gain shown (any cell)", "b: C2 gain shown (any cell)",
    "e: C2 extra gain over C1 shown (any cell)", "largest C1-not-P count (cell)", "largest C2-not-P count (cell)", "skip rule: M3 skipped (neither a nor b)"),
    value = c(sprintf("%.6f", Z_GAIN), as.character(GAIN_THRESHOLD), as.character(nrow(main)), as.character(nrow(eligible)),
      paste(main$cell[!main$cell_complete], collapse = " "), as.character(a), as.character(b), as.character(e),
      sprintf("%s (cell %s)", best("c1_not_p")[["count"]], best("c1_not_p")[["cell"]]),
      sprintf("%s (cell %s)", best("c2_not_p")[["count"]], best("c2_not_p")[["cell"]]), as.character(skip)), stringsAsFactors = FALSE)
  utils::write.table(verdict, file.path(out_dir, "m2-verdict.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  # O-3 (1): among the cells n = 24, g = 2, F = 300 with MCAR or AD at 5, 10, 20 %, the largest C1-not-P count; ties: AD first, then
  # the lower rate
  pool <- main[main$n == 24L & main$g == 2L & main$f == 300L & main$missing_id %in% c("MCAR5", "MCAR10", "MCAR20", "AD5", "AD10", "AD20"), ]
  pool$ad <- startsWith(pool$missing_id, "AD")
  pool$rate <- as.numeric(sub("^[A-Z]+", "", pool$missing_id))
  pool <- pool[order(-pool$c1_not_p, -pool$ad, pool$rate), ]
  utils::write.table(pool[c("cell", "missing_id", "n_valid", "c1_not_p", "c1_not_p_lo", "ad", "rate", "cell_complete")],
    file.path(out_dir, "m2-stress.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  list(a = a, b = b, e = e, skip = skip, stress_cell = if (nrow(pool)) pool$cell[[1L]] else NA_integer_, cells = cells_df, gain = gain)
}

# ---- self-test (E5, the M2 part): constructed inputs, outcomes as constructed ----
constructedSummary <- function(cells, c1_extra, c2_extra, n_per = 40L, report_only = FALSE) {
  rows <- list()
  for (cell in cells) {
    for (r in seq_len(n_per)) {
      p <- r > c1_extra[[as.character(cell)]] + c2_extra[[as.character(cell)]]
      c1 <- p || r <= c1_extra[[as.character(cell)]]
      c2 <- TRUE
      rows[[length(rows) + 1L]] <- data.frame(cell = cell, r = r, n = 24L, g = 2L, f = 300L,
        missing_id = if (cell == 1L) "MCAR5" else "AD10", report_only = report_only, cap = 3L, error = "", p_ruv = p, c1_ruv = c1, c2_ruv = c2,
        c1_differs_from_p = FALSE, c2_differs_from_p = FALSE, p_no_ruv_cause = if (p) "" else if (c1) "rescued_by_c1" else "completeness_only",
        n_r5 = 0L, n_r5_cap_side = 0L, n_r5_completeness = 0L, complete_features = 3, share_complete = 0.1, p_k_used = 3L, p_resid_df = 19L,
        c1_k_used = 3L, c1_resid_df = 19L, c2_k_used = 2L, c2_resid_df = 20L, stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

selfTest <- function(out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  checks <- data.frame(check = character(), pass = logical(), detail = character(), stringsAsFactors = FALSE)
  add <- function(check, pass, detail = "") checks[nrow(checks) + 1L, ] <<- list(check, isTRUE(pass), as.character(detail))
  # Wilson against hand values (e5hand-k262.py computes the same numbers independently)
  w <- wilson(7, 40, Z_GAIN); add("wilson_gain_7_of_40", w[[1L]] >= GAIN_THRESHOLD, sprintf("lower %.6f", w[[1L]]))
  w6 <- wilson(6, 40, Z_GAIN); add("wilson_gain_6_of_40", w6[[1L]] < GAIN_THRESHOLD, sprintf("lower %.6f", w6[[1L]]))
  w2 <- wilson(2, 40, Z95); add("wilson_two_sided_2_of_40", abs(w2[[1L]] - 0.0138) < 5e-4 && abs(w2[[2L]] - 0.1650) < 5e-4,
    sprintf("%.4f to %.4f", w2[[1L]], w2[[2L]]))
  add("z_gain", abs(Z_GAIN - 3.4524) < 5e-4, sprintf("%.6f", Z_GAIN))
  # a constructed table: cell 1 with 7 C1 recoveries (gain shown), cell 2 with 6 (no gain)
  d_gain <- constructedSummary(c(1L, 2L), list("1" = 7L, "2" = 6L), list("1" = 0L, "2" = 0L))
  s <- scoreM2(d_gain, file.path(out_dir, "gain-case"))
  add("gain_shown_at_7_of_40", s$a && s$gain$a[s$gain$cell == 1L] && !s$gain$a[s$gain$cell == 2L], "cell 1 yes, cell 2 no")
  add("skip_not_applied_when_gain", !s$skip, "")
  # no gain anywhere: skip applies
  d_none <- constructedSummary(c(1L, 2L), list("1" = 6L, "2" = 0L), list("1" = 0L, "2" = 0L))
  sn <- scoreM2(d_none, file.path(out_dir, "none-case"))
  add("skip_applied_without_gain", sn$skip && !sn$a && !sn$b, "")
  # C2 only: 8 recoveries by C2 (not C1)
  d_c2 <- constructedSummary(c(1L), list("1" = 0L), list("1" = 8L))
  s2 <- scoreM2(d_c2, file.path(out_dir, "c2-case"))
  add("c2_only_gain_not_skipped", !s2$a && s2$b && !s2$skip, "b true, a false: M3 not skipped by the rule as declared")
  # an errored data set is never a denominator and flags its cell
  d_err <- d_gain
  d_err$error[d_err$cell == 1L & d_err$r == 1L] <- "constructed error"
  se <- scoreM2(d_err, file.path(out_dir, "error-case"))
  add("error_row_flags_cell", !se$cells$cell_complete[se$cells$cell == 1L] && se$cells$n_valid[se$cells$cell == 1L] == 39L && !se$gain$a[se$gain$cell == 1L],
    "cell incomplete, 39 valid, excluded from the gain rule")
  # report-only cells never enter the family
  d_rep <- constructedSummary(c(181L), list("181" = 20L), list("181" = 0L), report_only = TRUE)
  sr <- scoreM2(d_rep, file.path(out_dir, "report-case"))
  add("report_only_excluded", !sr$a && nrow(sr$gain) == 0L, "")
  # the stress cell: ties on the count broken AD first, then the lower rate
  cells <- c(1L, 2L, 3L, 4L)
  tab <- data.frame(cell = cells, n = 24L, g = 2L, f = 300L, missing_id = c("MCAR10", "AD10", "AD5", "AD20"), stringsAsFactors = FALSE)
  base <- constructedSummary(cells, setNames(as.list(rep(5L, 4L)), cells), setNames(as.list(rep(0L, 4L)), cells))
  base$missing_id <- tab$missing_id[match(base$cell, tab$cell)]
  ss <- scoreM2(base, file.path(out_dir, "stress-case"))
  add("stress_tie_break", identical(ss$stress_cell, 3L), "equal counts: AD first, then the lower rate (AD5 = cell 3)")
  # E6: the scorer twice on the same input, byte-identical
  scoreM2(d_gain, file.path(out_dir, "gain-case-again"))
  files <- list.files(file.path(out_dir, "gain-case"))
  add("E6_byte_identical", identical(files, list.files(file.path(out_dir, "gain-case-again"))) && all(
    unname(tools::md5sum(file.path(out_dir, "gain-case", files))) == unname(tools::md5sum(file.path(out_dir, "gain-case-again", files)))),
    paste(length(files), "files"))
  utils::write.table(checks, file.path(out_dir, "checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  print(checks)
  quit(status = if (all(checks$pass)) 0L else 1L)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2L && identical(args[[1L]], "selftest")) {
  selfTest(normalizePath(file.path(getwd(), args[[2L]]), mustWork = FALSE))
} else if (length(args) >= 3L && identical(args[[1L]], "m2")) {
  out_dir <- normalizePath(file.path(getwd(), args[[2L]]), mustWork = FALSE)
  res <- scoreM2(readSummary(args[-(1:2)]), out_dir)
  print(utils::read.delim(file.path(out_dir, "m2-verdict.tsv")))
}
