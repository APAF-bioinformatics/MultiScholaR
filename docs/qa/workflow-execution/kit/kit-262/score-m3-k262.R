# KIT-262 M3 scoring (ticket 6.2, 6.3, 6.4, 7.1, 7.2): C1 against P on the ten cells (VAL-235's eight and the stress cells S-N, S-E),
# paired by data set. Base R only. Evidence code: no package change. DRAFT: written, not sent; an auditor reads it first.
# Random numbers: only the sign-flip test (seed 262000002) and the percentile bootstrap (seed 262000001), both with the kinds named; the
# scorer run twice on the same input gives byte-identical files (E6).
#
#   Rscript score-m3-k262.R blind <out> <plan.tsv> <source dir> [<source dir> ...]
#       The rehearsal (L3): mechanics and controls only. Every estimand is NA (no share, rate, FDP, power, k distribution, test statistic,
#       interval or verdict); per cell: units, errors, checks passed and failed, seconds. A source dir holds summary.tsv.
#   Rscript score-m3-k262.R full <out> <plan.tsv> <m2-verdict.tsv> <source dir> [<source dir> ...]
#       The counted scoring: one run over the union of the jobs and their continuations; each data set scored once.
#   Rscript score-m3-k262.R missing <plan.tsv> <source dir> [<source dir> ...]
#       The continuation's unit list: the planned units no source holds a complete line for, as K262_M3_ONLY_UNITS (no value is shown).
#   Rscript score-m3-k262.R selftest <out> <expected.tsv>
#       Constructed cases for every path; the expected values are read from e5hand-m3-expected.tsv (worked in Python, e5hand-m3-k262.py).
#
# O-1 (ii) is in force (VAL-302): only C1 is tested; the family has no Bonferroni split; alpha = 0.05 one-sided, the margins are VAL-185's
# (null cells 0.10 on the share of data sets with any call; effect cells 0.05 on the mean false-discovery proportion); O-2 (a): "no
# worse" is judged on the cell average, the changed-data-set figures (6.3) are reported beside it. The outcome table reduces to Stop,
# K1 and K4 (C2 does not exist in M3). K1 needs C1 NW in all ten cells and M2's gain `a` (read from the M2 verdict).
# A cell is determined when every planned unit has a row; a data set with an error row is never a denominator and is counted as excluded.
# STOPS (R errors with a class): an instrument row (m3_instrument), a unit complete in two sources or twice in one (m3_unit_in_two_sources,
# m3_duplicate_key), a row outside the plan's replicate range (m3_range), a cell outside the ten (m3_unknown_cell), a seed that is not the
# declared formula (m3_seed), differing headers (m3_schema). A truncated last line of a source (an interrupted unit) is dropped.
ALPHA <- 0.05
Z <- stats::qnorm(1 - ALPHA)
Z95 <- stats::qnorm(0.975)
M_NULL <- 0.10
M_EFFECT <- 0.05
BOOT_B <- 10000L
SIGNFLIP_SEED <- 262000002L
BOOT_SEED <- 262000001L
NULL_CELLS <- c("N-r", "N-s", "N-w", "S-N")
EFFECT_CELLS <- c("E-r", "E-s", "E-p", "E-g", "E-w", "S-E")
M3_IDS <- c("N-r", "N-s", "E-r", "E-s", "E-p", "E-g", "N-w", "E-w", "S-N", "S-E")
M3_IDX <- stats::setNames(c(1:8, 9L, 10L), M3_IDS)
M3_STRESS <- c("S-N", "S-E")
m3Seed <- function(idx, stress, r) if (stress) 262200000 + 10000 * idx + r else 190000000 + 10000 * idx + r
m3Tag <- function(idx, r) sprintf("c%02d-r%d", idx, r)  # the unit tag of m3-k262.R (the continuation names units by it)
# the row format of m3-k262.R (E0, e0-m3-k262.py: the vector below equals M3_COLUMNS there)
M3_COLUMNS <- c("cell", "r", "seed", "purpose", "omic", "status", "class", "error", "n_samples", "cap", "n_features",
  "t_n_ok", "t_n_r3", "t_n_r4", "t_n_r5", "t_n_r6", "t_n_err", "t_n_skipped_no_controls", "t_n_rejected",
  "p_ruv", "p_reason", "p_pct", "p_k_sweep", "p_k_used", "p_resid_df", "p_calls", "p_any_call", "p_fdp", "p_power",
  "p_nt1", "p_nt2", "p_nt3", "p_np1", "p_np2", "p_np3",
  "c_source", "c_ruv", "c_pct", "c_k_used", "c_resid_df", "c_calls", "c_any_call", "c_fdp", "c_power",
  "c_nt1", "c_nt2", "c_nt3", "c_np1", "c_np2", "c_np3", "c_differs", "change_class",
  "e1", "e7", "e9", "e9_compared", "d_pct", "d_k", "db_k", "d_calls", "e10", "e10_direct_no_ruv", "e11", "e11_run", "e8", "f2",
  "seconds_p", "seconds_c", "seconds_other", "seconds_total")

stopClass <- function(class, message) stop(errorCondition(message, class = c(class, "m3_stop")))
num <- function(x) suppressWarnings(as.numeric(x))
lgl <- function(x) ifelse(x %in% "TRUE", TRUE, ifelse(x %in% "FALSE", FALSE, NA))

# ---- copy: docs/qa/workflow-execution/kit/kit-190/decide.R tangoZ, bergerBoos, skewness, effectTest, nullTest, verbatim except [S] ----
tangoZ <- function(b, c, n, m = M_NULL) {
  score <- function(q) b / (q + m) + (if (c > 0) c / q else 0) - 2 * (n - b - c) / (1 - 2 * q - m)
  lo <- 1e-12
  hi <- (1 - m) / 2 - 1e-12
  q <- if (score(hi) > 0) hi else if (score(lo) < 0) lo else stats::uniroot(score, c(lo, hi), tol = 1e-12)$root
  (b - c - n * m) / sqrt(n * (2 * q + m * (1 - m)))
}

bergerBoos <- function(b, c, n, m = M_NULL, gamma = ALPHA / 10) {
  t_obs <- b - c
  ci <- stats::binom.test(b + c, n, conf.level = 1 - gamma)$conf.int
  lo <- max(ci[[1L]], m)
  hi <- min(ci[[2L]], 1)
  if (lo > hi) return(gamma)
  prob <- function(delta) {
    p10 <- (delta + m) / 2
    D <- 0:n
    w <- stats::dbinom(D, n, delta)
    # given D discordant pairs, b ~ Bin(D, p10 / delta); T = 2b - D <= t  <=>  b <= floor((t + D) / 2)
    sum(w * stats::pbinom(floor((t_obs + D) / 2), D, p10 / delta))
  }
  grid <- seq(lo, hi, length.out = 200)
  min(1, max(vapply(grid, prob, numeric(1L))) + gamma)
}

skewness <- function(e) {
  m2 <- mean((e - mean(e))^2)
  if (m2 <= 0) return(0)
  mean((e - mean(e))^3) / m2^1.5
}

effectTest <- function(d, m = M_EFFECT, vectors = 100000L, seed = SIGNFLIP_SEED) {  # [S] the seed is this ticket's (262000002)
  e <- d - m
  n <- length(e)
  s <- stats::sd(e)
  if (!is.finite(s) || s == 0) {
    # S4 (c): every paired difference equal; the studentised statistic is undefined; no worse iff mean(d) < margin.
    return(list(n = n, mean_d = mean(d), sd_d = 0, t = NA_real_, p_signflip = NA_real_, mc_se = NA_real_, skewness = 0,
      johnson_used = FALSE, t_johnson = NA_real_, t_upper_bound = mean(d), no_worse = mean(d) < m))
  }
  # K190-L01: the observed statistic goes through the same arithmetic as every sign vector's (the all-plus vector), so
  # the identity pattern ties exactly instead of being lost to rounding.
  tRows <- function(ex) {
    means <- rowMeans(ex)
    sds <- sqrt(rowSums((ex - means)^2) / (n - 1))
    means / (sds / sqrt(n))
  }
  t_obs <- tRows(matrix(e, 1L, n))
  set.seed(seed)
  count <- 0
  chunk <- max(1L, floor(5e6 / n))
  done <- 0L
  while (done < vectors) {
    k <- min(chunk, vectors - done)
    signs <- matrix(sample(c(-1, 1), k * n, replace = TRUE), k, n)
    tstar <- tRows(signs * rep(e, each = k))
    count <- count + sum(tstar <= t_obs)
    done <- done + k
  }
  p_flip <- (count + 1) / (vectors + 1)
  g <- skewness(e)
  mu3 <- mean((e - mean(e))^3)
  t_j <- (mean(e) + mu3 / (6 * s^2 * n) + mu3 * mean(e)^2 / (3 * s^4)) / (s / sqrt(n))
  johnson <- abs(g) > 1
  no_worse <- if (johnson) t_j < -Z else p_flip <= ALPHA
  list(n = n, mean_d = mean(d), sd_d = stats::sd(d), t = t_obs, p_signflip = p_flip, mc_se = sqrt(p_flip * (1 - p_flip) / vectors),
    skewness = g, johnson_used = johnson, t_johnson = t_j, t_upper_bound = mean(d) + stats::qt(1 - ALPHA, n - 1) * stats::sd(d) / sqrt(n),
    no_worse = no_worse)
}

nullTest <- function(d_call, f_call) {
  n <- length(d_call)
  b <- sum(d_call & !f_call)
  c <- sum(f_call & !d_call)
  z <- tangoZ(b, c, n)
  list(n = n, b = b, c = c, difference = (b - c) / n, tango_z = z, berger_boos_p = bergerBoos(b, c, n), no_worse = z < -Z)
}
# ---- end of the copy ----

# ---- copy: docs/qa/workflow-execution/kit/kit-190/bounds.R tangoZ (negative margins allowed), renamed tangoZb [S] ----
tangoZb <- function(b, c, n, m) {
  score <- function(q) b / (q + m) + (if (c > 0) c / q else 0) - 2 * (n - b - c) / (1 - 2 * q - m)
  lo <- max(0, -m) + 1e-12
  hi <- (1 - m) / 2 - 1e-12
  q <- if (score(hi) > 0) hi else if (score(lo) < 0) lo else stats::uniroot(score, c(lo, hi), tol = 1e-13)$root
  (b - c - n * m) / sqrt(n * (2 * q + m * (1 - m)))
}
# ---- end of the copy ----

# One-sided 95 % bounds on Delta = p10 - p01 by inverting Tango's score test (the upper bound is bounds.R's; the lower bound is its mirror,
# used only to label a not-NW cell "shown worse" or "inconclusive"); NA when the root is not bracketed.
tangoBound <- function(b, c, n, side) {
  est <- (b - c) / n
  f <- if (identical(side, "upper")) function(m) tangoZb(b, c, n, m) + Z else function(m) tangoZb(b, c, n, m) - Z
  lims <- if (identical(side, "upper")) c(est + 1e-9, 0.99) else c(-0.99, est - 1e-9)
  if (!all(is.finite(f(lims))) || f(lims[[1L]]) * f(lims[[2L]]) > 0) return(NA_real_)
  stats::uniroot(f, lims, tol = 1e-12)$root
}

wilson <- function(x, n, z) {
  if (!is.finite(n) || n == 0) return(c(NA_real_, NA_real_))
  p <- x / n
  d <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / d
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  c(centre - half, centre + half)
}

# Newcombe (1998), method 10, for p1 - p2 with p1 = (a + b) / n (the candidate) and p2 = (a + c) / n (the reference): a both, b candidate
# only, c reference only, d neither. The correlation correction of the paper (a zero marginal: phi 0).
newcombePaired <- function(a, b, c, d) {
  n <- a + b + c + d
  if (n == 0) return(c(NA_real_, NA_real_))
  p1 <- (a + b) / n; p2 <- (a + c) / n
  w1 <- wilson(x = a + b, n = n, z = Z95); w2 <- wilson(x = a + c, n = n, z = Z95)
  A <- (a + b) * (c + d) * (a + c) * (b + d)
  phi <- if (A == 0) 0 else {
    diff <- a * d - b * c
    C <- if (diff > n / 2) diff - n / 2 else if (diff >= 0) 0 else diff
    C / sqrt(A)
  }
  lo <- (p1 - p2) - sqrt((p1 - w1[[1L]])^2 + (w2[[2L]] - p2)^2 - 2 * phi * (p1 - w1[[1L]]) * (w2[[2L]] - p2))
  hi <- (p1 - p2) + sqrt((w1[[2L]] - p1)^2 + (p2 - w2[[1L]])^2 - 2 * phi * (w1[[2L]] - p1) * (p2 - w2[[1L]]))
  c(lo, hi)
}

# Paired-t two-sided 95 % interval and one-sided 95 % bounds of the mean of d.
tSummary <- function(d) {
  d <- d[is.finite(d)]
  n <- length(d)
  if (n < 2L) return(list(n = n, mean = if (n) mean(d) else NA_real_, half = NA_real_, upper = NA_real_, lower = NA_real_))
  s <- stats::sd(d)
  list(n = n, mean = mean(d), half = stats::qt(0.975, n - 1) * s / sqrt(n), upper = mean(d) + stats::qt(1 - ALPHA, n - 1) * s / sqrt(n),
    lower = mean(d) - stats::qt(1 - ALPHA, n - 1) * s / sqrt(n))
}

# The percentile interval of bootstrap replicates (type-7 quantiles at 0.025 and 0.975); NA when no replicate is finite.
pctInterval <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(NA_real_, NA_real_))
  stats::quantile(x, c(0.025, 0.975), type = 7, names = FALSE)
}

# ---- reading: sources, units, stops ------------------------------------------------------------------------------------------
readSource <- function(dir) {
  path <- file.path(dir, "summary.tsv")
  lines <- readLines(path, warn = FALSE)
  if (!length(lines)) stopClass("m3_schema", paste("empty summary:", path))
  header <- strsplit(lines[[1L]], "\t", fixed = TRUE)[[1L]]
  body <- lines[-1L]
  body <- body[nzchar(body)]
  fields <- strsplit(body, "\t", fixed = TRUE)
  whole <- lengths(fields) == length(header)  # a line cut short is an interrupted unit: dropped, never counted
  rows <- if (any(whole)) {
    m <- do.call(rbind, fields[whole])
    stats::setNames(as.data.frame(m, stringsAsFactors = FALSE), header)
  } else stats::setNames(as.data.frame(matrix(character(), 0L, length(header)), stringsAsFactors = FALSE), header)
  if (nrow(rows)) rows$source <- basename(dir) else rows$source <- character()
  list(rows = rows, dropped = sum(!whole), header = header)
}

# The replicates a cell is planned for: 1..R, or `reps` when given (the lifecycle's own seeds 9001 and 9002; never in a counted scoring).
repsOf <- function(plan, cell, reps = NULL) if (is.null(reps)) seq_len(plan$R[plan$cell == cell]) else as.integer(reps)

loadUnits <- function(dirs, plan, reps = NULL) {
  parts <- lapply(dirs, readSource)
  headers <- lapply(parts, `[[`, "header")
  if (!all(vapply(headers, function(h) identical(h, headers[[1L]]), logical(1L)))) stopClass("m3_schema", "the sources' summary headers differ")
  if (!identical(headers[[1L]], M3_COLUMNS)) stopClass("m3_schema", "the summary header is not the row format of m3-k262.R")
  u_all <- do.call(rbind, lapply(parts, `[[`, "rows"))
  if (is.null(u_all) || !nrow(u_all)) u_all <- stats::setNames(as.data.frame(matrix(character(), 0L, length(M3_COLUMNS) + 1L), stringsAsFactors = FALSE), c(M3_COLUMNS, "source"))
  if (any(u_all$status %in% "instrument")) stopClass("m3_instrument", sprintf("scoring stops: %d instrument row(s) (VAL-54): %s", sum(u_all$status %in% "instrument"),
    paste(utils::head(paste(u_all$cell[u_all$status %in% "instrument"], u_all$r[u_all$status %in% "instrument"], u_all$class[u_all$status %in% "instrument"]), 3L), collapse = "; ")))
  if (!all(u_all$cell %in% M3_IDS)) stopClass("m3_unknown_cell", paste("a cell outside the ten:", paste(setdiff(unique(u_all$cell), M3_IDS), collapse = ", ")))
  key <- paste(u_all$cell, u_all$r, u_all$purpose, sep = "|")
  dup <- unique(key[duplicated(key)])
  if (length(dup)) {
    two <- vapply(dup, function(k) length(unique(u_all$source[key == k])) > 1L, logical(1L))
    if (any(two)) stopClass("m3_unit_in_two_sources", paste("a unit complete in two sources:", paste(utils::head(dup[two], 3L), collapse = ", ")))
    stopClass("m3_duplicate_key", paste("a unit twice in one source:", paste(utils::head(dup, 3L), collapse = ", ")))
  }
  main <- u_all[u_all$purpose == "main", , drop = FALSE]
  in_plan <- vapply(seq_len(nrow(main)), function(i) main$cell[[i]] %in% plan$cell && num(main$r[[i]]) %in% repsOf(plan, main$cell[[i]], reps), logical(1L))
  if (!all(in_plan)) stopClass("m3_range", "a row outside the plan's replicate range 1..R (or a cell not in the plan)")
  seed_ok <- vapply(seq_len(nrow(u_all)), function(i) isTRUE(num(u_all$seed[[i]]) == m3Seed(idx = M3_IDX[[u_all$cell[[i]]]], stress = u_all$cell[[i]] %in% M3_STRESS, r = num(u_all$r[[i]]))), logical(1L))
  if (!all(seed_ok)) stopClass("m3_seed", "a seed is not the declared formula")
  u_all
}

readPlan <- function(path) {
  p <- utils::read.delim(path, stringsAsFactors = FALSE, colClasses = "character")
  data.frame(cell = p$cell, R = as.integer(p$R), stringsAsFactors = FALSE)
}

# ---- the blind scorer: mechanics and controls only; every estimand NA --------------------------------------------------------------
blindM3 <- function(dirs, plan, out_dir, reps = NULL) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  u_all <- loadUnits(dirs, plan, reps)
  main <- u_all[u_all$purpose == "main", , drop = FALSE]
  rows <- lapply(plan$cell, function(cell) {
    u <- main[main$cell == cell, , drop = FALSE]
    flag <- function(col) c(pass = sum(u[[col]] %in% "TRUE"), fail = sum(u[[col]] %in% "FALSE"))
    f <- lapply(c("e1", "e7", "e9", "e11", "e8", "f2", "e10"), flag)
    names(f) <- c("e1", "e7", "e9", "e11", "e8", "f2", "e10")
    planned <- length(repsOf(plan, cell, reps))
    data.frame(cell = cell, planned = planned, present = nrow(u), ok = sum(u$status == "ok"), error_rows = sum(u$status == "error"),
      missing = length(setdiff(repsOf(plan, cell, reps), num(u$r))), e1_fail = f$e1[["fail"]], e7_fail = f$e7[["fail"]], e9_fail = f$e9[["fail"]], e11_fail = f$e11[["fail"]],
      e8_fail = f$e8[["fail"]], f2_fail = f$f2[["fail"]], e10_fail = f$e10[["fail"]], e1_pass = f$e1[["pass"]], e9_pass = f$e9[["pass"]], e11_pass = f$e11[["pass"]],
      e8_pass = f$e8[["pass"]], f2_pass = f$f2[["pass"]], seconds_mean = mean(num(u$seconds_total), na.rm = TRUE), seconds_max = suppressWarnings(max(num(u$seconds_total), na.rm = TRUE)),
      share_rescued = NA_real_, share_reranked = NA_real_, calibration = NA_real_, power = NA_real_, verdict = NA_character_, stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  utils::write.table(tab, file.path(out_dir, "m3-blind.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  tab
}

# ---- the counted scorer -----------------------------------------------------------------------------------------------------------
verdictCell <- function(cell, v) {
  kind <- if (cell %in% NULL_CELLS) "any_call" else "fdp"
  if (nrow(v) == 0L) return(data.frame(cell = cell, endpoint = kind, n = 0L, estimate = NA_real_, statistic = NA_real_, p_value = NA_real_, p_berger_boos = NA_real_,
    skewness = NA_real_, johnson_used = NA, no_worse = NA, lower_bound = NA_real_, upper_bound = NA_real_, label = "not determined (no valid pair)", stringsAsFactors = FALSE))
  if (kind == "any_call") {
    r <- nullTest(lgl(v$c_any_call), lgl(v$p_any_call))
    lower <- tangoBound(b = r$b, c = r$c, n = r$n, side = "lower"); upper <- tangoBound(b = r$b, c = r$c, n = r$n, side = "upper")
    nw <- r$no_worse
    data.frame(cell = cell, endpoint = kind, n = r$n, estimate = r$difference, statistic = r$tango_z, p_value = NA_real_, p_berger_boos = r$berger_boos_p, skewness = NA_real_,
      johnson_used = NA, no_worse = nw, lower_bound = lower, upper_bound = upper,
      label = if (isTRUE(nw)) "NW" else if (isTRUE(lower > 0)) "shown worse" else "inconclusive", stringsAsFactors = FALSE)
  } else {
    d <- num(v$c_fdp) - num(v$p_fdp)
    r <- effectTest(d)
    ts <- tSummary(d)
    lower <- if (is.finite(ts$lower)) ts$lower else mean(d)
    upper <- if (is.finite(ts$upper)) ts$upper else mean(d)
    data.frame(cell = cell, endpoint = kind, n = r$n, estimate = r$mean_d, statistic = if (r$johnson_used) r$t_johnson else r$t, p_value = r$p_signflip, p_berger_boos = NA_real_,
      skewness = r$skewness, johnson_used = r$johnson_used, no_worse = r$no_worse, lower_bound = lower, upper_bound = upper,
      label = if (isTRUE(r$no_worse)) "NW" else if (isTRUE(lower > 0)) "shown worse" else "inconclusive", stringsAsFactors = FALSE)
  }
}

# 6.3: the data sets in which C1's choice differs from P's: effect cells the paired mean FDP difference with its paired-t 95 % interval,
# null cells b - c' over those data sets with Newcombe's interval; split into rescued and re-ranked.
changedCell <- function(cell, v) {
  one <- function(label, u) {
    n <- nrow(u)
    if (cell %in% NULL_CELLS) {
      cc <- lgl(u$c_any_call); pp <- lgl(u$p_any_call)
      a <- sum(cc & pp); b <- sum(cc & !pp); c <- sum(!cc & pp); d <- sum(!cc & !pp)
      ci <- newcombePaired(a = a, b = b, c = c, d = d)
      data.frame(cell = cell, subset = label, n = n, endpoint = "any_call", estimate = if (n) (b - c) / n else NA_real_, lo = ci[[1L]], hi = ci[[2L]], b = b, c = c, stringsAsFactors = FALSE)
    } else {
      ts <- tSummary(num(u$c_fdp) - num(u$p_fdp))
      data.frame(cell = cell, subset = label, n = n, endpoint = "fdp", estimate = ts$mean, lo = ts$mean - ts$half, hi = ts$mean + ts$half, b = NA_real_, c = NA_real_, stringsAsFactors = FALSE)
    }
  }
  ch <- v[lgl(v$c_differs) %in% TRUE, , drop = FALSE]
  rbind(one("changed", ch), one("rescued", ch[ch$change_class == "rescued", , drop = FALSE]), one("reranked", ch[ch$change_class == "reranked", , drop = FALSE]))
}

calibrationCell <- function(cell, v, B = BOOT_B) {
  n <- nrow(v)
  nt_p <- cbind(num(v$p_nt1), num(v$p_nt2), num(v$p_nt3)); np_p <- cbind(num(v$p_np1), num(v$p_np2), num(v$p_np3))
  nt_c <- cbind(num(v$c_nt1), num(v$c_nt2), num(v$c_nt3)); np_c <- cbind(num(v$c_np1), num(v$c_np2), num(v$c_np3))
  point <- function(np, nt) colSums(np) / colSums(nt)
  pt_p <- point(np_p, nt_p); pt_c <- point(np_c, nt_c)
  if (n < 2L) return(data.frame(cell = cell, third = 1:3, n_units = n, tested_p = NA_real_, alpha_p = NA_real_, lo_p = NA_real_, hi_p = NA_real_, alpha_c1 = NA_real_,
    lo_c1 = NA_real_, hi_c1 = NA_real_, diff = NA_real_, lo_diff = NA_real_, hi_diff = NA_real_, stringsAsFactors = FALSE))
  rp <- matrix(NA_real_, B, 3L); rc <- matrix(NA_real_, B, 3L)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    rp[b, ] <- colSums(np_p[idx, , drop = FALSE]) / colSums(nt_p[idx, , drop = FALSE])
    rc[b, ] <- colSums(np_c[idx, , drop = FALSE]) / colSums(nt_c[idx, , drop = FALSE])
  }
  do.call(rbind, lapply(1:3, function(i) {
    ip <- pctInterval(rp[, i]); ic <- pctInterval(rc[, i]); id <- pctInterval(rc[, i] - rp[, i])
    data.frame(cell = cell, third = i, n_units = n, tested_p = sum(nt_p[, i]), alpha_p = pt_p[[i]], lo_p = ip[[1L]], hi_p = ip[[2L]], alpha_c1 = pt_c[[i]], lo_c1 = ic[[1L]], hi_c1 = ic[[2L]],
      diff = pt_c[[i]] - pt_p[[i]], lo_diff = id[[1L]], hi_diff = id[[2L]], stringsAsFactors = FALSE)
  }))
}

powerCell <- function(cell, v) {
  d <- num(v$c_power) - num(v$p_power)
  ts <- tSummary(d)
  data.frame(cell = cell, n = ts$n, power_p = mean(num(v$p_power), na.rm = TRUE), power_c1 = mean(num(v$c_power), na.rm = TRUE), diff = ts$mean, lo = ts$mean - ts$half, hi = ts$mean + ts$half, stringsAsFactors = FALSE)
}

coverageCell <- function(cell, v) {
  n <- nrow(v)
  p <- lgl(v$p_ruv) %in% TRUE; c1 <- lgl(v$c_ruv) %in% TRUE; ch <- lgl(v$c_differs) %in% TRUE
  w <- function(x) wilson(x = x, n = n, z = Z95)
  data.frame(cell = cell, n = n, p_ruv = sum(p), p_ruv_lo = w(sum(p))[[1L]], p_ruv_hi = w(sum(p))[[2L]], c1_ruv = sum(c1), c1_ruv_lo = w(sum(c1))[[1L]], c1_ruv_hi = w(sum(c1))[[2L]],
    changed = sum(ch), changed_lo = w(sum(ch))[[1L]], changed_hi = w(sum(ch))[[2L]], rescued = sum(v$change_class %in% "rescued"), reranked = sum(v$change_class %in% "reranked"),
    p_pct_median = suppressWarnings(stats::median(num(v$p_pct), na.rm = TRUE)), c1_pct_median = suppressWarnings(stats::median(num(v$c_pct), na.rm = TRUE)),
    e1_fail = sum(v$e1 %in% "FALSE"), stringsAsFactors = FALSE)
}

# 6.4 (VAL-289): per design (samples, groups 2) and rule, the search's k used and the residual degrees of freedom left.
residTable <- function(v_all) {
  do.call(rbind, lapply(split(v_all, num(v_all$n_samples)), function(u) do.call(rbind, lapply(c("p", "c"), function(a) {
    has <- lgl(u[[paste0(a, "_ruv")]]) %in% TRUE
    k <- num(u[[paste0(a, "_k_used")]])[has]; rd <- num(u[[paste0(a, "_resid_df")]])[has]; cap <- num(u$cap)[has]
    data.frame(n_samples = num(u$n_samples[[1L]]), rule = if (a == "p") "P" else "C1", with_ruv = sum(has), k_min = if (length(k)) min(k) else NA_real_, k_max = if (length(k)) max(k) else NA_real_,
      k_at_cap = sum(k == cap), resid_df_median = if (length(rd)) stats::median(rd) else NA_real_, resid_df_le1 = sum(rd <= 1), resid_df_0 = sum(rd == 0), stringsAsFactors = FALSE)
  }))))
}

# E10 (stress cells): the journey's "no usable percentage" under P against the direct M2-form sweep before between-sample normalisation.
e10Table <- function(v_all) {
  s <- v_all[v_all$cell %in% M3_STRESS & v_all$status == "ok", , drop = FALSE]
  if (!nrow(s)) return(data.frame(cell = character(), n = integer(), both_no = integer(), journey_only = integer(), direct_only = integer(), both_ruv = integer(),
    difference = numeric(), lo = numeric(), hi = numeric(), stringsAsFactors = FALSE))
  do.call(rbind, lapply(split(s, s$cell), function(u) {
    j <- !(lgl(u$p_ruv) %in% TRUE); d <- lgl(u$e10_direct_no_ruv) %in% TRUE
    a <- sum(j & d); b <- sum(j & !d); c <- sum(!j & d); dd <- sum(!j & !d)
    ci <- newcombePaired(a = a, b = b, c = c, d = dd)
    data.frame(cell = u$cell[[1L]], n = nrow(u), both_no = a, journey_only = b, direct_only = c, both_ruv = dd, difference = (b - c) / nrow(u), lo = ci[[1L]], hi = ci[[2L]], stringsAsFactors = FALSE)
  }))
}

outcomeOf <- function(verdicts, plan, a_gain, check_fail, e11_fail, missing_cells) {
  if (check_fail) return(list(outcome = "Stop", reason = "a check failed in a counted data set (E1, E7 = F1, E9 = F3, F2)"))
  if (e11_fail) return(list(outcome = "Stop (E11)", reason = "E11 failed: the candidate rows of the re-used data sets are to be run (ticket 9.3); the coordinator decides"))
  not_nw <- verdicts$cell[!(verdicts$label %in% "NW")]
  if (!length(not_nw) && isTRUE(a_gain)) return(list(outcome = "K1", reason = "C1 NW in all deciding cells and gain shown in M2: to adopt C1 in a later implementation ticket"))
  why <- if (!length(not_nw)) "C1 is NW in every cell but M2 showed no gain (a is FALSE): keep the present rule with M2's limit" else
    paste0("not NW: ", paste(sprintf("%s (%s)", verdicts$cell[!(verdicts$label %in% "NW")], verdicts$label[!(verdicts$label %in% "NW")]), collapse = "; "),
      ". Inconclusive and not-determined cells are named; they are not shown worse unless labelled so.")
  list(outcome = "K4", reason = why)
}

fullM3 <- function(dirs, plan, a_gain, out_dir, B = BOOT_B, reps = NULL) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  u_all <- loadUnits(dirs, plan, reps)
  main <- u_all[u_all$purpose == "main", , drop = FALSE]
  e9rows <- u_all[u_all$purpose == "e9", , drop = FALSE]
  ok <- main[main$status == "ok", , drop = FALSE]
  # the counted checks: E1, E7 (F1), E9 (F3, on the main units and the E9 units), F2 are falsifiers; E11 has its own stop; E8, E10 are reported
  bad <- function(col, rows) sum(rows[[col]] %in% "FALSE")
  hard <- bad("e1", ok) + bad("e7", ok) + bad("e9", ok) + bad("e9", e9rows[e9rows$status == "ok", , drop = FALSE]) + bad("f2", ok)
  e11_bad <- bad("e11", ok)
  cells <- plan$cell
  per <- lapply(cells, function(cell) {
    u <- main[main$cell == cell, , drop = FALSE]
    v <- ok[ok$cell == cell, , drop = FALSE]
    planned <- length(repsOf(plan, cell, reps))
    missing <- length(setdiff(repsOf(plan, cell, reps), num(u$r)))
    list(cell = cell, planned = planned, present = nrow(u), valid = nrow(v), errors = sum(u$status == "error"), missing = missing, v = v)
  })
  names(per) <- cells
  det <- vapply(per, function(x) x$missing == 0L, logical(1L))
  verdicts <- do.call(rbind, lapply(per, function(x) {
    if (x$missing > 0L) return(data.frame(cell = x$cell, endpoint = if (x$cell %in% NULL_CELLS) "any_call" else "fdp", n = x$valid, estimate = NA_real_, statistic = NA_real_, p_value = NA_real_,
      p_berger_boos = NA_real_, skewness = NA_real_, johnson_used = NA, no_worse = NA, lower_bound = NA_real_, upper_bound = NA_real_, label = "not determined (did not complete)", stringsAsFactors = FALSE))
    verdictCell(cell = x$cell, v = x$v)
  }))
  rownames(verdicts) <- NULL
  verdicts$excluded <- vapply(verdicts$cell, function(cell) per[[cell]]$errors, numeric(1L))
  verdicts$missing <- vapply(verdicts$cell, function(cell) per[[cell]]$missing, numeric(1L))
  verdicts$owner_asked <- verdicts$excluded > 5
  utils::write.table(verdicts, file.path(out_dir, "m3-verdicts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  set.seed(BOOT_SEED, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  usable <- cells[det & vapply(per, function(x) x$valid > 0L, logical(1L))]
  write2 <- function(x, name) utils::write.table(if (is.null(x)) data.frame() else x, file.path(out_dir, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  write2(do.call(rbind, lapply(usable, function(cell) changedCell(cell = cell, v = per[[cell]]$v))), "m3-changed.tsv")
  write2(do.call(rbind, lapply(usable, function(cell) calibrationCell(cell = cell, v = per[[cell]]$v, B = B))), "m3-calibration.tsv")
  write2(do.call(rbind, lapply(intersect(usable, EFFECT_CELLS), function(cell) powerCell(cell = cell, v = per[[cell]]$v))), "m3-power.tsv")
  write2(do.call(rbind, lapply(usable, function(cell) coverageCell(cell = cell, v = per[[cell]]$v))), "m3-coverage.tsv")
  write2(residTable(ok), "m3-resid-df.tsv")
  write2(e10Table(ok), "m3-e10.tsv")
  counts <- data.frame(cell = cells, planned = vapply(per, `[[`, 0L, "planned"), present = vapply(per, `[[`, 0L, "present"), valid = vapply(per, `[[`, 0L, "valid"),
    error_rows = vapply(per, `[[`, 0L, "errors"), missing = vapply(per, `[[`, 0L, "missing"), stringsAsFactors = FALSE)
  write2(counts, "m3-units.tsv")
  out <- outcomeOf(verdicts = verdicts, plan = plan, a_gain = a_gain, check_fail = hard > 0L, e11_fail = e11_bad > 0L, missing_cells = cells[!det])
  write2(data.frame(outcome = out$outcome, reason = out$reason, a_gain_m2 = a_gain, falsifier_failures = hard, e11_failures = e11_bad, stringsAsFactors = FALSE), "m3-outcome.tsv")
  list(outcome = out, verdicts = verdicts, counts = counts)
}

# The CONTINUATION's unit list (ticket 9.5; declared with the job): the planned main units that no source holds a complete line for, as the
# tags `cNN-rR` that K262_M3_ONLY_UNITS takes. Counted from cell, r and purpose only: no value is read or shown. An instrument row stops it
# (an instrument row is not a missing unit). An errored data set is NOT missing (it has its row and is never re-run).
missingTags <- function(dirs, plan, reps = NULL) {
  u <- loadUnits(dirs, plan, reps)
  main <- u[u$purpose == "main", , drop = FALSE]
  have <- paste(main$cell, main$r)
  tags <- character()
  for (cell in plan$cell) for (r in repsOf(plan, cell, reps)) {
    if (!(paste(cell, r) %in% have)) tags <- c(tags, m3Tag(idx = M3_IDX[[cell]], r = r))
  }
  tags
}

readGainA <- function(path) {
  v <- utils::read.delim(path, stringsAsFactors = FALSE, colClasses = "character")
  x <- v$value[startsWith(v$item, "a: C1 gain shown")]
  if (length(x) != 1L) stopClass("m3_schema", "the M2 verdict has no unique 'a: C1 gain shown' row")
  identical(x, "TRUE")
}

# ---- the self-test (E5, M3 part): constructed inputs, outcomes as constructed, expected values from the Python hand values ----------------
mkUnit <- function(cell, r, ...) {
  idx <- M3_IDX[[cell]]
  row <- stats::setNames(rep(list("NA"), length(M3_COLUMNS)), M3_COLUMNS)
  row[c("cell", "r", "seed", "purpose", "omic", "status", "class", "error", "n_samples", "cap")] <- list(cell, as.character(r),
    sprintf("%.15g", m3Seed(idx = idx, stress = cell %in% M3_STRESS, r = r)), "main", "lipid", "ok", "", "", if (cell %in% M3_STRESS) "24" else "12", if (cell %in% M3_STRESS) "3" else "2")
  row[c("p_ruv", "c_ruv", "c_differs")] <- list("TRUE", "TRUE", "FALSE")
  row[c("p_pct", "c_pct", "p_k_used", "c_k_used", "p_resid_df", "c_resid_df", "p_calls", "c_calls")] <- list("5", "5", "2", "2", "8", "8", "0", "0")
  row[c("p_any_call", "c_any_call")] <- list("FALSE", "FALSE")
  row[c("p_fdp", "c_fdp", "p_power", "c_power")] <- list("0.1", "0.1", "0.5", "0.5")
  for (i in 1:3) { row[[paste0("p_nt", i)]] <- "10"; row[[paste0("c_nt", i)]] <- "10"; row[[paste0("p_np", i)]] <- "0"; row[[paste0("c_np", i)]] <- "0" }
  row[c("change_class", "c_source", "e1", "e7", "e9", "e11", "e8", "f2", "e10", "e10_direct_no_ruv", "seconds_total")] <- list("same", "reused", "TRUE", "TRUE", "NA", "NA", "NA", "NA", "NA", "NA", "30")
  extra <- list(...)
  for (nm in names(extra)) row[[nm]] <- as.character(extra[[nm]])
  row
}
lineOf <- function(row) paste(unlist(row[M3_COLUMNS]), collapse = "\t")
writeSource <- function(dir, rows, extra_lines = character()) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(c(paste(M3_COLUMNS, collapse = "\t"), vapply(rows, lineOf, ""), extra_lines), file.path(dir, "summary.tsv"))
}
fullPlan <- function(R) data.frame(cell = M3_IDS, R = as.integer(R), stringsAsFactors = FALSE)

selfTest <- function(out_dir, expected_path) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  exp <- utils::read.delim(expected_path, stringsAsFactors = FALSE, colClasses = "character")
  ev <- function(key) { v <- exp$value[exp$key == key]; stopifnot(length(v) == 1L); as.numeric(v) }
  checks <- data.frame(check = character(), pass = logical(), detail = character(), stringsAsFactors = FALSE)
  add <- function(check, pass, detail = "") checks[nrow(checks) + 1L, ] <<- list(check, isTRUE(pass), as.character(detail))
  near <- function(a, b, tol = 1e-6) isTRUE(all(abs(a - b) <= tol))
  stopsWith <- function(expr, class) tryCatch({ force(expr); FALSE }, m3_stop = function(e) inherits(e, class))
  # --- statistics against the Python hand values
  add("tango_b0c0", near(tangoZ(b = 0, c = 0, n = 100), ev("tango_z.b0_c0_n100_m0.10")), sprintf("%.6f", tangoZ(b = 0, c = 0, n = 100)))
  add("tango_b40c0", near(tangoZ(b = 40, c = 0, n = 100), ev("tango_z.b40_c0_n100_m0.10")), sprintf("%.6f", tangoZ(b = 40, c = 0, n = 100)))
  add("bound_upper", near(tangoBound(b = 10, c = 5, n = 100, side = "upper"), ev("bound_upper.b10_c5_n100")), sprintf("%.6f", tangoBound(b = 10, c = 5, n = 100, side = "upper")))
  add("bound_lower", near(tangoBound(b = 10, c = 5, n = 100, side = "lower"), ev("bound_lower.b10_c5_n100")), sprintf("%.6f", tangoBound(b = 10, c = 5, n = 100, side = "lower")))
  nc <- newcombePaired(a = 10, b = 5, c = 3, d = 82)
  add("newcombe_10_5_3_82", near(nc, c(ev("newcombe_lo.a10_b5_c3_d82"), ev("newcombe_hi.a10_b5_c3_d82"))), sprintf("%.6f to %.6f", nc[[1L]], nc[[2L]]))
  nz <- newcombePaired(a = 0, b = 7, c = 0, d = 3)
  add("newcombe_zero_marginal", near(nz, c(ev("newcombe_lo.a0_b7_c0_d3"), ev("newcombe_hi.a0_b7_c0_d3"))), sprintf("%.6f to %.6f", nz[[1L]], nz[[2L]]))
  w <- wilson(x = 5, n = 40, z = Z95)
  add("wilson_5_of_40", near(w, c(ev("wilson_lo.x5_n40"), ev("wilson_hi.x5_n40"))), sprintf("%.6f to %.6f", w[[1L]], w[[2L]]))
  d6 <- c(0.1, 0.2, 0.0, 0.3, 0.1, 0.2); t6 <- tSummary(d6)
  add("t_interval_and_bounds", near(c(t6$mean, t6$half, t6$upper, t6$lower), c(ev("t_mean.d6"), ev("t_half95.d6"), ev("t_upper_bound.d6"), ev("t_lower_bound.d6"))), sprintf("half %.6f", t6$half))
  ej <- effectTest(c(rep(0, 170), rep(0.20, 30)))
  add("johnson_used_and_t", ej$johnson_used && near(ej$t_johnson, ev("johnson_t.d170_30")) && near(ej$skewness, ev("johnson_skew.d170_30")) && ej$no_worse, sprintf("t %.6f", ej$t_johnson))
  es <- effectTest(c(0.03, 0.02))
  add("signflip_exact_two_pairs", abs(es$p_signflip - ev("signflip_p.n2_003_002")) <= 0.01 && !es$no_worse, sprintf("p %.4f", es$p_signflip))
  add("pct_interval_quantile7", near(pctInterval(1:101), c(ev("quantile7_lo.1to101"), ev("quantile7_hi.1to101"))) &&
    near(pctInterval(c(0, 0.1, 0.1, 0.3, 0.9)), c(ev("quantile7_lo.five"), ev("quantile7_hi.five"))) && all(is.na(pctInterval(c(NA_real_, NaN)))), "type 7; all-NA gives NA")
  # --- sources and plan
  src <- function(name) file.path(out_dir, "src", name)
  eqR <- fullPlan(rep(10L, 10L))
  # identical arms: every cell NW, K1 with a = TRUE
  all_equal <- function(cells = M3_IDS, rs = 100L) do.call(c, lapply(cells, function(cell) lapply(seq_len(rs), function(r) mkUnit(cell, r, p_fdp = sprintf("%.3f", 0.1 + 0.001 * (r %% 5)), c_fdp = sprintf("%.3f", 0.1 + 0.001 * (r %% 5))))))
  writeSource(src("identical"), all_equal())
  res <- fullM3(src("identical"), fullPlan(rep(100L, 10L)), TRUE, file.path(out_dir, "case-identical"), B = 200L)
  add("identical_arms_NW_K1", all(res$verdicts$label == "NW") && identical(res$outcome$outcome, "K1"), paste(res$verdicts$label, collapse = ","))
  res_a <- fullM3(src("identical"), fullPlan(rep(100L, 10L)), FALSE, file.path(out_dir, "case-identical-a0"), B = 200L)
  add("all_NW_but_a_false_K4", identical(res_a$outcome$outcome, "K4") && grepl("a is FALSE", res_a$outcome$reason, fixed = TRUE), res_a$outcome$reason)
  # a +0.20 shift in every pair: effect cells not NW and shown worse; null cells with 40 of 100 candidate-only calls not NW and shown worse
  shifted <- do.call(c, lapply(M3_IDS, function(cell) lapply(1:100, function(r) {
    if (cell %in% NULL_CELLS) mkUnit(cell, r, p_any_call = "FALSE", c_any_call = if (r <= 40) "TRUE" else "FALSE") else
      mkUnit(cell, r, p_fdp = sprintf("%.3f", 0.1 + 0.002 * (r %% 7)), c_fdp = sprintf("%.3f", 0.3 + 0.002 * (r %% 7)))
  })))
  writeSource(src("shifted"), shifted)
  rs <- fullM3(src("shifted"), fullPlan(rep(100L, 10L)), TRUE, file.path(out_dir, "case-shifted"), B = 200L)
  add("shift_not_NW_shown_worse", all(rs$verdicts$label == "shown worse") && identical(rs$outcome$outcome, "K4"), paste(unique(rs$verdicts$label), collapse = ","))
  add("shift_null_statistic_matches_hand", near(rs$verdicts$statistic[rs$verdicts$cell == "N-r"], ev("tango_z.b40_c0_n100_m0.10")), "")
  # inconclusive: an effect cell with mean difference 0.06 and large spread (lower bound below 0), a null cell with 15 candidate-only of 100
  # effect cell E-r: 30 data sets with FDP 0 -> 1, 24 with 1 -> 0, 46 equal: mean difference 0.06 above the margin, SD about 0.73: not NW, lower bound below 0
  inc_e2 <- lapply(1:100, function(r) if (r <= 30) mkUnit("E-r", r, p_fdp = "0", c_fdp = "1") else if (r <= 54) mkUnit("E-r", r, p_fdp = "1", c_fdp = "0") else mkUnit("E-r", r, p_fdp = "0.5", c_fdp = "0.5"))
  rest <- do.call(c, lapply(setdiff(M3_IDS, "E-r"), function(cell) lapply(1:100, function(r) mkUnit(cell, r))))
  writeSource(src("inconclusive"), c(inc_e2, rest))
  ri <- fullM3(src("inconclusive"), fullPlan(rep(100L, 10L)), TRUE, file.path(out_dir, "case-inconclusive"), B = 200L)
  lab <- ri$verdicts$label[ri$verdicts$cell == "E-r"]
  add("inconclusive_label", identical(lab, "inconclusive") && identical(ri$outcome$outcome, "K4") && grepl("E-r (inconclusive)", ri$outcome$reason, fixed = TRUE), lab)
  # missing units: a cell with 95 of 100 rows is "did not complete", never tested; a zero-valid-pair (all error rows) cell is "not determined (no valid pair)"
  miss <- c(do.call(c, lapply(setdiff(M3_IDS, c("N-r", "E-p")), function(cell) lapply(1:100, function(r) mkUnit(cell, r)))),
    lapply(1:95, function(r) mkUnit("N-r", r)), lapply(1:100, function(r) mkUnit("E-p", r, status = "error", class = "omics_ruv_failed", error = "constructed")))
  writeSource(src("missing"), miss)
  rmiss <- fullM3(src("missing"), fullPlan(rep(100L, 10L)), TRUE, file.path(out_dir, "case-missing"), B = 200L)
  add("missing_unit_not_determined", identical(rmiss$verdicts$label[rmiss$verdicts$cell == "N-r"], "not determined (did not complete)") &&
    identical(rmiss$verdicts$label[rmiss$verdicts$cell == "E-p"], "not determined (no valid pair)") && identical(rmiss$outcome$outcome, "K4") &&
    grepl("N-r (not determined (did not complete))", rmiss$outcome$reason, fixed = TRUE), paste(rmiss$outcome$reason))
  add("error_rows_never_a_denominator", rmiss$counts$valid[rmiss$counts$cell == "E-p"] == 0L && rmiss$counts$error_rows[rmiss$counts$cell == "E-p"] == 100L && rmiss$verdicts$excluded[rmiss$verdicts$cell == "E-p"] == 100, "E-p: 0 valid, 100 error rows")
  # continuation: the missing five units of N-r come from a second source; the union is scored once; the cell is then determined
  writeSource(src("cont-a"), c(do.call(c, lapply(setdiff(M3_IDS, "N-r"), function(cell) lapply(1:100, function(r) mkUnit(cell, r)))), lapply(1:95, function(r) mkUnit("N-r", r))))
  writeSource(src("cont-b"), lapply(96:100, function(r) mkUnit("N-r", r)))
  rc <- fullM3(c(src("cont-a"), src("cont-b")), fullPlan(rep(100L, 10L)), TRUE, file.path(out_dir, "case-continuation"), B = 200L)
  add("continuation_union_scored_once", rc$counts$missing[rc$counts$cell == "N-r"] == 0L && rc$counts$valid[rc$counts$cell == "N-r"] == 100L && rc$verdicts$label[rc$verdicts$cell == "N-r"] == "NW", "95 + 5 units, 100 valid")
  # a unit complete in two sources stops scoring; a unit twice in one source stops too
  writeSource(src("cont-dup"), lapply(95:100, function(r) mkUnit("N-r", r)))
  add("unit_in_two_sources_stops", stopsWith(fullM3(c(src("cont-a"), src("cont-dup")), fullPlan(rep(100L, 10L)), TRUE, file.path(out_dir, "case-dup"), B = 50L), "m3_unit_in_two_sources"), "unit 95 in both")
  writeSource(src("dup-one"), c(lapply(1:100, function(r) mkUnit("N-r", r)), list(mkUnit("N-r", 7))))
  add("unit_twice_in_one_source_stops", stopsWith(loadUnits(src("dup-one"), eqR), "m3_duplicate_key"), "")
  # an interrupted unit: a truncated last line is dropped; complete in the continuation -> counted once; only truncated -> missing
  trunc <- substr(lineOf(mkUnit("N-r", 100)), 1L, 40L)
  writeSource(src("int-a"), c(do.call(c, lapply(setdiff(M3_IDS, "N-r"), function(cell) lapply(1:100, function(r) mkUnit(cell, r)))), lapply(1:99, function(r) mkUnit("N-r", r))), extra_lines = trunc)
  writeSource(src("int-b"), list(mkUnit("N-r", 100)))
  ri2 <- fullM3(c(src("int-a"), src("int-b")), fullPlan(rep(100L, 10L)), TRUE, file.path(out_dir, "case-interrupted"), B = 200L)
  ri1 <- fullM3(src("int-a"), fullPlan(rep(100L, 10L)), TRUE, file.path(out_dir, "case-interrupted-only"), B = 200L)
  add("interrupted_unit_dropped_then_continued", ri2$counts$missing[ri2$counts$cell == "N-r"] == 0L && ri1$counts$missing[ri1$counts$cell == "N-r"] == 1L &&
    ri1$verdicts$label[ri1$verdicts$cell == "N-r"] == "not determined (did not complete)", "99 + 1 truncated")
  add("missing_tags_for_continuation", identical(missingTags(src("int-a"), fullPlan(rep(100L, 10L))), "c01-r100") &&
    length(missingTags(c(src("int-a"), src("int-b")), fullPlan(rep(100L, 10L)))) == 0L &&
    identical(missingTags(src("missing"), fullPlan(rep(100L, 10L))), sprintf("c01-r%d", 96:100)), "one interrupted unit; none after the continuation; five units of a short cell; error rows are not missing")
  # an instrument row, a row outside the plan's range, an unknown cell, a wrong seed, differing headers stop scoring
  writeSource(src("instr"), c(lapply(1:99, function(r) mkUnit("N-r", r)), list(mkUnit("N-r", 100, status = "instrument", class = "partial_match_guard", error = "constructed"))))
  add("instrument_row_stops", stopsWith(loadUnits(src("instr"), eqR), "m3_instrument"), "")
  writeSource(src("range"), c(lapply(1:100, function(r) mkUnit("N-r", r)), list(mkUnit("N-r", 9001))))
  add("row_outside_range_stops", stopsWith(loadUnits(src("range"), fullPlan(rep(100L, 10L))), "m3_range"), "r = 9001 in a plan of 100")
  unk <- mkUnit("N-r", 1); unk[["cell"]] <- "Z-z"  # (the row is built first: `cell` is also mkUnit's own formal)
  writeSource(src("unk"), list(unk))
  add("unknown_cell_stops", stopsWith(loadUnits(src("unk"), eqR), "m3_unknown_cell"), "")
  writeSource(src("seed"), list(mkUnit("N-r", 1, seed = "12345")))
  add("wrong_seed_stops", stopsWith(loadUnits(src("seed"), eqR), "m3_seed"), "")
  dir.create(src("hdr"), recursive = TRUE, showWarnings = FALSE)
  writeLines(c(paste(rev(M3_COLUMNS), collapse = "\t")), file.path(src("hdr"), "summary.tsv"))
  add("header_mismatch_stops", stopsWith(loadUnits(src("hdr"), eqR), "m3_schema"), "")
  # the falsifiers: E1, E7, E9, F2 false in a counted row -> Stop; E11 false -> Stop (E11); E8 and E10 false are reported only
  stopcase <- function(col, expect, value = "FALSE") {
    rows <- do.call(c, lapply(M3_IDS, function(cell) lapply(1:100, function(r) { u <- mkUnit(cell, r); if (cell == "N-r" && r == 1L) u[[col]] <- value; u })))
    writeSource(src(paste0("stop-", col)), rows)
    fullM3(src(paste0("stop-", col)), fullPlan(rep(100L, 10L)), TRUE, file.path(out_dir, paste0("case-stop-", col)), B = 50L)$outcome$outcome
  }
  add("E1_false_stops", identical(stopcase("e1", "Stop"), "Stop"), "")
  add("E7_false_stops", identical(stopcase("e7", "Stop"), "Stop"), "")
  add("E9_false_stops", identical(stopcase("e9", "Stop"), "Stop"), "")
  add("F2_false_stops", identical(stopcase("f2", "Stop"), "Stop"), "")
  add("E11_false_own_stop", identical(stopcase("e11", "Stop (E11)"), "Stop (E11)"), "")
  add("E8_false_reported_only", identical(stopcase("e8", "K1"), "K1"), "E8 false: outcome unchanged (reported in the blind and unit tables)")
  add("E10_false_reported_only", identical(stopcase("e10", "K1"), "K1"), "")
  # changed data sets: 5 rescued (P without RUV, C1 with) and 3 re-ranked in an effect cell: the paired mean FDP difference among the changed
  chg <- c(lapply(1:5, function(r) mkUnit("E-r", r, p_ruv = "FALSE", c_ruv = "TRUE", c_differs = "TRUE", change_class = "rescued", c_source = "manual", p_fdp = "0.4", c_fdp = "0.2")),
    lapply(6:8, function(r) mkUnit("E-r", r, c_differs = "TRUE", change_class = "reranked", c_source = "manual", p_fdp = "0.1", c_fdp = "0.2")), lapply(9:20, function(r) mkUnit("E-r", r)))
  ch <- changedCell("E-r", do.call(rbind, lapply(chg, function(x) as.data.frame(x, stringsAsFactors = FALSE))))
  add("changed_split_counts_and_means", identical(ch$n, c(8L, 5L, 3L)) && near(ch$estimate, c((5 * -0.2 + 3 * 0.1) / 8, -0.2, 0.1), 1e-9), paste(ch$estimate, collapse = ","))
  # null cell: the changed data sets' b and c' with Newcombe's interval (a = 10 both, b = 5, c = 3, d = 82 over 100 -- here over the changed set)
  nul <- c(lapply(1:10, function(r) mkUnit("N-r", r, c_differs = "TRUE", change_class = "reranked", p_any_call = "TRUE", c_any_call = "TRUE")),
    lapply(11:15, function(r) mkUnit("N-r", r, c_differs = "TRUE", change_class = "rescued", p_any_call = "FALSE", c_any_call = "TRUE")),
    lapply(16:18, function(r) mkUnit("N-r", r, c_differs = "TRUE", change_class = "reranked", p_any_call = "TRUE", c_any_call = "FALSE")),
    lapply(19:100, function(r) mkUnit("N-r", r, c_differs = "TRUE", change_class = "reranked", p_any_call = "FALSE", c_any_call = "FALSE")))
  chn <- changedCell("N-r", do.call(rbind, lapply(nul, function(x) as.data.frame(x, stringsAsFactors = FALSE))))
  add("changed_null_newcombe_matches_hand", chn$n[[1L]] == 100L && near(c(chn$lo[[1L]], chn$hi[[1L]]), c(ev("newcombe_lo.a10_b5_c3_d82"), ev("newcombe_hi.a10_b5_c3_d82"))) && near(chn$estimate[[1L]], 0.02, 1e-12), "")
  add("changed_subset_with_no_units_is_NA", chn$n[chn$subset == "rescued"] == 5L && is.finite(chn$lo[chn$subset == "rescued"]) &&
    is.na(changedCell("E-r", do.call(rbind, lapply(chg[9:12], function(x) as.data.frame(x, stringsAsFactors = FALSE))))$estimate[[1L]]), "")
  # calibration by third: degenerate (every unit alike: the interval is the point) and non-degenerate (units differ: lo < point < hi within [min, max])
  cal_rows <- function(np_p, np_c) lapply(seq_along(np_p), function(i) mkUnit("N-r", i, p_np1 = np_p[[i]][[1L]], p_np2 = np_p[[i]][[2L]], p_np3 = np_p[[i]][[3L]], c_np1 = np_c[[i]][[1L]], c_np2 = np_c[[i]][[2L]], c_np3 = np_c[[i]][[3L]]))
  df_of <- function(rows) do.call(rbind, lapply(rows, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
  same <- df_of(cal_rows(rep(list(c(1, 0, 0)), 6), rep(list(c(1, 1, 0)), 6)))
  cd <- calibrationCell("N-r", same, B = 300L)
  add("calibration_degenerate_interval_equals_point", near(cd$lo_p, cd$alpha_p, 1e-12) && near(cd$hi_p, cd$alpha_p, 1e-12) && near(cd$lo_diff, cd$diff, 1e-12) && near(cd$hi_diff, cd$diff, 1e-12), "")
  four_p <- list(c(1, 0, 0), c(0, 1, 0), c(0, 0, 0), c(1, 1, 0)); four_c <- list(c(1, 1, 0), c(0, 1, 1), c(0, 0, 1), c(1, 1, 1))
  var <- df_of(cal_rows(four_p, four_c))
  cv <- calibrationCell("N-r", var, B = 2000L)
  add("calibration_ratio_matches_hand", near(cv$alpha_p, c(ev("calib_ratio_P.third1"), ev("calib_ratio_P.third2"), ev("calib_ratio_P.third3"))) &&
    near(cv$alpha_c1, c(ev("calib_ratio_C1.third1"), ev("calib_ratio_C1.third2"), ev("calib_ratio_C1.third3"))) && near(cv$diff, c(ev("calib_diff.third1"), ev("calib_diff.third2"), ev("calib_diff.third3"))), "")
  add("calibration_nondegenerate_interval", all(cv$lo_c1 < cv$alpha_c1 + 1e-12 & cv$alpha_c1 - 1e-12 < cv$hi_c1) && all(cv$hi_c1 > cv$lo_c1) && all(cv$lo_c1 >= 0) && all(cv$hi_c1 <= 1) && all(cv$hi_diff >= cv$lo_diff) && all(cv$hi_diff[2:3] > cv$lo_diff[2:3]), "lo < point < hi, in [0, 1]")
  # residual degrees of freedom: k used and n - g - k
  rt <- residTable(df_of(c(lapply(1:4, function(r) mkUnit("N-r", r, p_k_used = "2", p_resid_df = "8", c_k_used = "2", c_resid_df = "8")), lapply(5:6, function(r) mkUnit("N-r", r, p_k_used = "3", p_resid_df = "1", c_k_used = "3", c_resid_df = "1", cap = "3")))))
  add("resid_df_table", rt$with_ruv[rt$rule == "P"] == 6L && rt$k_at_cap[rt$rule == "P"] == 4L + 2L && rt$resid_df_le1[rt$rule == "P"] == 2L && rt$resid_df_0[rt$rule == "P"] == 0L && rt$resid_df_median[rt$rule == "P"] == 8, "")
  # E10: the journey's no-RUV share against the direct sweep (McNemar counts)
  e10 <- e10Table(df_of(c(lapply(1:4, function(r) mkUnit("S-N", r, p_ruv = "FALSE", e10_direct_no_ruv = "TRUE")), lapply(5:6, function(r) mkUnit("S-N", r, p_ruv = "TRUE", e10_direct_no_ruv = "TRUE")),
    lapply(7:8, function(r) mkUnit("S-N", r, p_ruv = "TRUE", e10_direct_no_ruv = "FALSE")))))
  add("e10_mcnemar_counts", e10$both_no == 4L && e10$direct_only == 2L && e10$journey_only == 0L && e10$both_ruv == 2L, "")
  # the blind scorer: every estimand NA, mechanics only; stops on an instrument row
  writeSource(src("blind"), all_equal(rs = 100L))
  bl <- blindM3(src("blind"), fullPlan(rep(100L, 10L)), file.path(out_dir, "case-blind"))
  add("blind_all_estimands_NA", all(is.na(bl$share_rescued)) && all(is.na(bl$calibration)) && all(is.na(bl$power)) && all(is.na(bl$verdict)) && all(bl$present == 100L) && all(bl$e1_fail == 0L), "")
  add("blind_stops_on_instrument", stopsWith(blindM3(src("instr"), eqR, file.path(out_dir, "case-blind-instr")), "m3_instrument"), "")
  # an E9 unit (purpose e9) with a false E9 counts as a falsifier; its row is not a denominator
  e9u <- c(all_equal(), list(mkUnit("N-r", 1, purpose = "e9", e9 = "FALSE", e9_compared = "TRUE")))
  writeSource(src("e9"), e9u)
  add("E9_unit_false_stops", identical(fullM3(src("e9"), fullPlan(rep(100L, 10L)), TRUE, file.path(out_dir, "case-e9"), B = 50L)$outcome$outcome, "Stop"), "")
  # the lifecycle's own seeds (9001, 9002): outside 1..R they stop a counted scoring, and are read only under an explicit `reps` (the lifecycle's
  # blind and full runs); the full scorer then runs on two pairs per cell without error
  writeSource(src("lc"), do.call(c, lapply(M3_IDS, function(cell) lapply(9001:9002, function(r) mkUnit(cell, r)))))
  lcplan <- fullPlan(rep(2L, 10L))
  lc_full <- tryCatch(fullM3(src("lc"), lcplan, TRUE, file.path(out_dir, "case-lc-full"), B = 100L, reps = c(9001L, 9002L)), error = function(e) e)
  add("lifecycle_reps_override", stopsWith(loadUnits(src("lc"), lcplan), "m3_range") && nrow(loadUnits(src("lc"), lcplan, reps = c(9001L, 9002L))) == 20L &&
    all(blindM3(src("lc"), lcplan, file.path(out_dir, "case-lc-blind"), reps = c(9001L, 9002L))$missing == 0L), "")
  # with two pairs and equal arms the null cells' Tango statistic is -sqrt(2) x sqrt(0.1 / 0.9) = -0.4714, not below -1.6449: not NW, lower bound below 0
  add("full_scorer_runs_on_two_pairs_per_cell", !inherits(lc_full, "error") && identical(lc_full$outcome$outcome, "K4") && grepl("N-r (inconclusive)", lc_full$outcome$reason, fixed = TRUE) &&
    isTRUE(abs(lc_full$verdicts$statistic[lc_full$verdicts$cell == "N-r"] + sqrt(2) * sqrt(0.1 / 0.9)) < 1e-6), if (inherits(lc_full, "error")) conditionMessage(lc_full) else "K4 on two equal pairs: the null cells are inconclusive")
  # E6: the scorer twice on the same input, byte-identical
  fullM3(src("shifted"), fullPlan(rep(100L, 10L)), TRUE, file.path(out_dir, "case-shifted-again"), B = 200L)
  fa <- list.files(file.path(out_dir, "case-shifted")); fb <- list.files(file.path(out_dir, "case-shifted-again"))
  add("E6_byte_identical", identical(fa, fb) && length(fa) > 0L && all(unname(tools::md5sum(file.path(out_dir, "case-shifted", fa))) == unname(tools::md5sum(file.path(out_dir, "case-shifted-again", fb)))), paste(length(fa), "files"))
  # M2's verdict file: the unique `a` row
  vf <- file.path(out_dir, "m2-verdict-a.tsv")
  writeLines(c("item\tvalue", "a: C1 gain shown (any cell)\tTRUE", "b: C2 gain shown (any cell)\tTRUE"), vf)
  add("m2_verdict_a_read", isTRUE(readGainA(vf)), "")
  utils::write.table(checks, file.path(out_dir, "checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  print(checks)
  quit(status = if (all(checks$pass)) 0L else 1L)
}

args <- commandArgs(trailingOnly = TRUE)
# Sourced by m3-k262.R's lifecycle (its own first argument is an output directory): no command runs; the functions above are the interface.
if (!length(args) || !args[[1L]] %in% c("selftest", "blind", "full", "missing")) {
  invisible(NULL)
} else if (length(args) >= 3L && identical(args[[1L]], "selftest")) {
  selfTest(normalizePath(file.path(getwd(), args[[2L]]), mustWork = FALSE), args[[3L]])
} else if (length(args) >= 4L && identical(args[[1L]], "blind")) {
  tab <- blindM3(dirs = args[-(1:3)], plan = readPlan(args[[3L]]), out_dir = normalizePath(file.path(getwd(), args[[2L]]), mustWork = FALSE))
  print(tab[c("cell", "planned", "present", "ok", "error_rows", "missing", "e1_fail", "e9_fail", "e11_fail", "f2_fail", "seconds_mean")])
} else if (length(args) >= 5L && identical(args[[1L]], "full")) {
  res <- fullM3(dirs = args[-(1:4)], plan = readPlan(args[[3L]]), a_gain = readGainA(args[[4L]]), out_dir = normalizePath(file.path(getwd(), args[[2L]]), mustWork = FALSE))
  print(res$verdicts[c("cell", "n", "label", "excluded", "missing")])
  cat(sprintf("outcome: %s\n%s\n", res$outcome$outcome, res$outcome$reason))
} else if (length(args) >= 3L && identical(args[[1L]], "missing")) {
  tags <- missingTags(dirs = args[-(1:2)], plan = readPlan(args[[2L]]))
  cat(sprintf("units missing: %d\nK262_M3_ONLY_UNITS=%s\n", length(tags), paste(tags, collapse = ",")))
} else stop("usage: score-m3-k262.R blind|full|missing|selftest ...")
