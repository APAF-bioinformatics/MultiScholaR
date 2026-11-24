# KIT-190 bound-declaration.md sections 5 and 6: the decision of the bounded-k measurement (VAL-182 (3), VAL-185).
# Evidence code: no package code; base R only.
#
# Usage (run-tests scripts lane):
#   Rscript decide.R decide <journeys.tsv> [<journeys.tsv> ...] -- <out dir>
#       Rows of `effect.R bound stage2` (all seeds). Per deciding cell, the paired comparison of D+B with F+B:
#       null cells (N-r, N-s, N-w): Tango's score test of H0 p10 - p01 >= 0.10, the Berger-Boos exact p-value beside it;
#       effect cells (E-r, E-s, E-p, E-g, E-w): the studentised sign-flip test of H0 mean(d) >= 0.05 (100,000 sign
#       vectors, seed 19001), Johnson's statistic when |g| > 1 (R4-01), the paired t bound beside it;
#       intersection-union over the 8 cells, each one-sided at alpha = 0.05. E1-T1b (unit = seed) and T11 reported.
#       Writes verdicts.tsv, reported.tsv and decision.tsv.
#   Rscript decide.R selftest <out dir>
#       Lifecycle item (vi): constructed inputs, identical arms (must pass) and a +0.20 shift in every pair (must fail),
#       for both tests. Exit status 1 on a failure.
ALPHA <- 0.05
Z <- stats::qnorm(1 - ALPHA)
NULL_CELLS <- c("N-r", "N-s", "N-w")
EFFECT_CELLS <- c("E-r", "E-s", "E-p", "E-g", "E-w")
M_NULL <- 0.10
M_EFFECT <- 0.05

# Tango (1998): restricted MLE of p01 under p10 - p01 = m, then the score statistic.
tangoZ <- function(b, c, n, m = M_NULL) {
  score <- function(q) b / (q + m) + (if (c > 0) c / q else 0) - 2 * (n - b - c) / (1 - 2 * q - m)
  lo <- 1e-12
  hi <- (1 - m) / 2 - 1e-12
  q <- if (score(hi) > 0) hi else if (score(lo) < 0) lo else stats::uniroot(score, c(lo, hi), tol = 1e-12)$root
  (b - c - n * m) / sqrt(n * (2 * q + m * (1 - m)))
}

# Berger-Boos: statistic T = b - c, small values reject; nuisance delta = p10 + p01 restricted to its 1 - gamma
# Clopper-Pearson interval (gamma = alpha / 10), on the H0 boundary p10 - p01 = m; p = sup P(T <= t) + gamma.
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

# The studentised sign-flip test of H0 mean(d) >= m (one-sided, lower), and Johnson's statistic when |g| > 1.
effectTest <- function(d, m = M_EFFECT, vectors = 100000L, seed = 19001L) {
  e <- d - m
  n <- length(e)
  s <- stats::sd(e)
  t_obs <- mean(e) / (s / sqrt(n))
  set.seed(seed)
  count <- 0
  chunk <- max(1L, floor(5e6 / n))
  done <- 0L
  while (done < vectors) {
    k <- min(chunk, vectors - done)
    signs <- matrix(sample(c(-1, 1), k * n, replace = TRUE), k, n)
    ex <- signs * rep(e, each = k)
    means <- rowMeans(ex)
    sds <- sqrt(rowSums((ex - means)^2) / (n - 1))
    tstar <- means / (sds / sqrt(n))
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

pairs <- function(rows, cell) {
  f <- rows[rows$cell == cell & rows$arm == "F+B" & (is.na(rows$error) | !nzchar(rows$error)), , drop = FALSE]
  d <- rows[rows$cell == cell & rows$arm == "D+B" & (is.na(rows$error) | !nzchar(rows$error)), , drop = FALSE]
  key <- function(x) paste(x$r, x$dataset)
  common <- intersect(key(f), key(d))
  list(f = f[match(common, key(f)), , drop = FALSE], d = d[match(common, key(d)), , drop = FALSE],
    excluded = length(union(key(f), key(d))) - length(common))
}

decide <- function(rows) {
  verdicts <- list()
  for (cell in c(NULL_CELLS, EFFECT_CELLS)) {
    p <- pairs(rows, cell)
    if (cell %in% NULL_CELLS) {
      r <- nullTest(as.logical(p$d$any_call), as.logical(p$f$any_call))
      verdicts[[cell]] <- data.frame(cell = cell, endpoint = "any_call", n = r$n, excluded = p$excluded,
        estimate = r$difference, statistic = r$tango_z, p_value = NA, p_berger_boos = r$berger_boos_p, skewness = NA,
        johnson_used = NA, no_worse = r$no_worse, stringsAsFactors = FALSE)
    } else {
      r <- effectTest(as.numeric(p$d$fdp) - as.numeric(p$f$fdp))
      verdicts[[cell]] <- data.frame(cell = cell, endpoint = "fdp", n = r$n, excluded = p$excluded, estimate = r$mean_d,
        statistic = if (r$johnson_used) r$t_johnson else r$t, p_value = r$p_signflip, p_berger_boos = NA,
        skewness = r$skewness, johnson_used = r$johnson_used, no_worse = r$no_worse, stringsAsFactors = FALSE)
    }
  }
  verdicts <- do.call(rbind, verdicts)
  # Reported only: E1-T1b by seed (mean of a seed's two journeys), T11 fixtures, and power and k differences.
  reported <- list()
  t1b <- pairs(rows[rows$cell == "E1" & startsWith(as.character(rows$dataset), "t1b_"), , drop = FALSE], "E1")
  if (nrow(t1b$f)) {
    seed <- sub("^t1b_([0-9]+)_.*$", "\\1", t1b$f$dataset)
    dd <- tapply(as.numeric(t1b$d$fdp) - as.numeric(t1b$f$fdp), seed, mean)
    r <- effectTest(as.numeric(dd))
    reported[["T1b"]] <- data.frame(item = "E1-T1b FDP by seed", n = r$n, estimate = r$mean_d, p_value = r$p_signflip,
      note = sprintf("no worse at 0.05: %s (reported, not deciding)", r$no_worse), stringsAsFactors = FALSE)
  }
  for (fx in c("t11_concentration", "t11_f172")) for (arm in c("F+B", "D+B")) {
    x <- rows[rows$cell == "E1" & rows$dataset == fx & rows$arm == arm, , drop = FALSE]
    if (nrow(x)) reported[[paste(fx, arm)]] <- data.frame(item = paste(fx, arm), n = 1, estimate = x$calls[[1L]], p_value = NA,
      note = sprintf("calls (all false); k = %s", x$k[[1L]]), stringsAsFactors = FALSE)
  }
  for (cell in c(NULL_CELLS, EFFECT_CELLS)) {
    p <- pairs(rows, cell)
    for (v in c("k", if (cell %in% EFFECT_CELLS) c("power", "abs_lfc_error"))) {
      d <- as.numeric(p$d[[v]]) - as.numeric(p$f[[v]])
      d <- d[is.finite(d)]
      if (length(d) > 1L) {
        h <- stats::qt(0.975, length(d) - 1) * stats::sd(d) / sqrt(length(d))
        reported[[paste(cell, v)]] <- data.frame(item = paste(cell, v, "D+B minus F+B"), n = length(d), estimate = mean(d),
          p_value = NA, note = sprintf("95%% t interval [%.4f, %.4f], descriptive", mean(d) - h, mean(d) + h),
          stringsAsFactors = FALSE)
      }
    }
  }
  list(verdicts = verdicts, reported = do.call(rbind, reported),
    decision = data.frame(all_cells_no_worse = all(verdicts$no_worse), ships = if (all(verdicts$no_worse)) "D+B" else "F+B"))
}

args <- commandArgs(TRUE)
if (identical(args[[1L]], "decide")) {
  sep <- which(args == "--")
  files <- args[2:(sep - 1L)]
  out <- args[[sep + 1L]]
  rows <- do.call(rbind, lapply(files, utils::read.delim, stringsAsFactors = FALSE))
  res <- decide(rows)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(res)) utils::write.table(res[[nm]], file.path(out, paste0(if (nm == "verdicts") "verdicts" else nm, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  print(res$verdicts)
  print(res$decision)
} else if (identical(args[[1L]], "selftest")) {
  out <- args[[2L]]
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  set.seed(190777)
  n <- 200
  f <- stats::runif(n, 0, 0.5)
  same_effect <- effectTest((f + stats::rnorm(n, 0, 0.05)) - f)$no_worse
  shift_effect <- effectTest((f + 0.20 + stats::rnorm(n, 0, 0.05)) - f)$no_worse
  calls <- stats::runif(n) < 0.1
  same_null <- nullTest(calls, calls)$no_worse
  shifted <- calls
  shifted[which(!calls)[seq_len(round(0.20 * n))]] <- TRUE
  shift_null <- nullTest(shifted, calls)$no_worse
  checks <- c(effect_identical_passes = same_effect, effect_shift_fails = !shift_effect,
    null_identical_passes = same_null, null_shift_fails = !shift_null)
  utils::write.table(data.frame(check = names(checks), ok = checks), file.path(out, "decide-selftest.tsv"), sep = "\t",
    quote = FALSE, row.names = FALSE)
  for (nm in names(checks)) cat(sprintf("selftest %s: %s\n", nm, if (checks[[nm]]) "ok" else "FAILED"))
  if (!all(checks)) quit(status = 1L)
} else stop("usage: decide.R decide|selftest ...")
