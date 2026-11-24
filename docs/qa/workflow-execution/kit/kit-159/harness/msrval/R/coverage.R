# The coverage regimes T-COV (KIT-159 section 10.2): generators, the regime catalogue, the pass
# bands, the chunked runner and the chunk speed-up on simulations (f_sim, G14 (c)). Every
# simulated data set is scored at every level a consumer declares; levels add no simulations.

COVERAGE_LEVELS <- c(0.95, 0.975, 1 - 0.05 / 3, 1 - 0.05 / 4, 1 - 0.05 / 6, 1 - 0.05 / 8,
                     1 - 0.05 / 12, 1 - 0.05 / 15)
TARGET_SIMS <- 20000L

#' The harness regime cells of T-COV
#'
#' One row per simulation cell: the regime, its method, form, unit count, the data condition and
#' the family of generator. Statistic classes are set per statistic inside a family.
#'
#' @param confirmatory Also the confirmatory cell DI10K (VAL-163 (3)), after T-COV's cells and
#'   outside T-COV: C4e at g = 5 on 4,840 features with 40 changed at effect scale 1 and C4e's
#'   `accuracy` and `threshold_count` statistics (family `c4e:4840acc`), with its own index; and
#'   after it the confirmatory cell P1K (VAL-176): the rank statistic on one contrast (group C
#'   against A, 40 positives of 4,840) at the rank cells' effect rule applied to that contrast
#'   (family `c4e:4840p1`), with its own index.
#' @return A data frame of cells with a stable `cell` id.
#' @export
valCoverageCells <- function(confirmatory = FALSE) {
    cells <- list()
    add <- function(...) cells[[length(cells) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)  # nolint
    for (n in c(6L, 7L, 8L, 9L, 10L, 13L, 21L)) for (shape in c("normal", "skewed")) {
        add(regime = "C1", family = "exact_median", n = n, condition = shape)
    }
    for (n in c(4L, 5L, 6L, 8L, 13L, 21L, 100L, 800L)) for (shape in c("normal", "skewed")) {
        add(regime = "C2", family = "paired_t", n = n, condition = shape)
    }
    c3 <- valC3Configs()
    for (cfg in names(c3)) for (cls in c("accuracy", "threshold_count")) {
        for (tr in c("identity", "log")) for (shape in c("normal", "gamma")) {
            for (tau in c(0.02, 0.05)) {
                add(regime = "C3", family = paste("c3", cfg, cls, tr, sep = ":"),
                    n = c3[[cfg]]$g, condition = paste0(shape, "_sd", tau))
            }
        }
    }
    for (g in c(4L, 6L)) {
        add(regime = "C4", family = "c4:500", n = g, condition = "matrix500")
        add(regime = "C4", family = "c4:4840", n = g, condition = "matrix4840")
    }
    for (g in c(5L, 14L)) {
        add(regime = "C4e", family = "c4e:500", n = g, condition = "matrix500_block")
        add(regime = "C4e", family = "c4e:4840", n = g, condition = "matrix4840_block")
    }
    for (n in c(21L, 48L)) for (p in c(0.95, 0.8)) {
        add(regime = "C5", family = "wilson", n = n, condition = paste0("p", p))
    }
    for (n in c(21L, 48L)) for (pp in c("0.9_0.8", "0.9_0.9")) {
        add(regime = "C6", family = "newcombe_paired", n = n, condition = paste0("p", pp))
    }
    for (cls in c("accuracy", "threshold_count")) {
        add(regime = "C7", family = paste0("c7:", cls), n = 200L, condition = "iid_in_class")
    }
    for (n in c(10L, 21L, 250L)) for (p in c(0.05, 0.5, 0.95)) {
        add(regime = "C8", family = "clopper_pearson", n = n, condition = paste0("p", p))
    }
    out <- do.call(rbind, cells)
    out$cell <- sprintf("%s-%03d", out$regime, seq_len(nrow(out)))
    if (isTRUE(confirmatory)) {
        out <- rbind(out, data.frame(regime = "C4e", family = "c4e:4840acc", n = 5L,
                                     condition = "matrix4840acc_block",
                                     cell = sprintf("C4e-%03d", nrow(out) + 1L),
                                     stringsAsFactors = FALSE))
        out <- rbind(out, data.frame(regime = "C4e", family = "c4e:4840p1", n = 5L,
                                     condition = "matrix4840p1_block",
                                     cell = sprintf("C4e-%03d", nrow(out) + 1L),
                                     stringsAsFactors = FALSE))
    }
    out
}

#' Designs of regime C3 (jackknife form (i) over designs built from fewer units)
#'
#' Pairs from 5, 8 and 21 batches (KIT-163); k + k from two pools (KIT-164: 3 + 3 of 7 and 14,
#' 2 + 2 of 3 and 5, 2 + 2 of 7 and 14), all designs when at most 400, else 200 drawn once with a
#' fixed seed; and the laboratory grouping (8 laboratories of 1 to 4 batches over the 3 + 3 designs;
#' laboratories are consecutive batches, sizes 4, 3, 3, 3, 3, 2, 2, 1).
#' @noRd
valC3Configs <- function() {
    pairs <- function(n) list(designs = t(utils::combn(n, 2L)), units = as.list(seq_len(n)), g = n)
    kk <- function(k, n1, n2) {
        a <- utils::combn(n1, k, simplify = FALSE)
        b <- lapply(utils::combn(n2, k, simplify = FALSE), function(x) x + n1)
        grid <- expand.grid(i = seq_along(a), j = seq_along(b))
        if (nrow(grid) > 400L) grid <- grid[valWithSeed(159300L, sample.int(nrow(grid), 200L)), ]
        d <- t(vapply(seq_len(nrow(grid)), function(r) c(a[[grid$i[r]]], b[[grid$j[r]]]),
                      numeric(2L * k)))
        list(designs = d, units = as.list(seq_len(n1 + n2)), g = n1 + n2)
    }
    en <- kk(3L, 7L, 14L)
    lab_sizes <- c(4L, 3L, 3L, 3L, 3L, 2L, 2L, 1L)
    labs <- split(seq_len(21L), rep(seq_along(lab_sizes), lab_sizes))
    list(pairs5 = pairs(5L), pairs8 = pairs(8L), pairs21 = pairs(21L), en3x3 = en,
         r8_2x2 = kk(2L, 3L, 5L), r21_2x2 = kk(2L, 7L, 14L),
         labs8 = list(designs = en$designs, units = unname(labs), g = 8L))
}

#' The pass band of a simulated coverage (section 4.4)
#'
#' @param kind `exact`, `clopper`, `approx_run` or `approx_within`.
#' @param coverage Simulated coverage.
#' @param level Nominal level.
#' @param achieved Achieved level (exact methods).
#' @param nsim Simulations.
#' @return `pass`, `conservative`, `fail`, or `not_offered` (an exact method with no interval).
#' @export
valCoverageStatus <- function(kind, coverage, level, achieved = NA_real_, nsim = 2000L) {
    if (kind == "exact" && !is.finite(achieved)) return("not_offered")
    if (!is.finite(coverage)) return("fail")
    if (kind == "exact") {
        if (!is.finite(achieved)) return("not_offered")
        se <- sqrt(achieved * (1 - achieved) / nsim)
        return(if (abs(coverage - achieved) <= 3 * se + 1e-12) "pass" else "fail")
    }
    if (kind == "clopper") {
        se <- sqrt(level * (1 - level) / nsim)
        return(if (coverage >= level - 3 * se - 1e-12) "conservative" else "fail")
    }
    if (coverage < level - 0.02 - 1e-12) return("fail")
    # r1.13 (VAL-69): run-level `conservative` is read on the miss rate, below a quarter of the
    # nominal miss rate, so that it can be seen at every level; within-run keeps its upper band
    if (kind == "approx_run") {
        return(if (1 - coverage < (1 - level) / 4 - 1e-12) "conservative" else "pass")
    }
    if (coverage > min(level + 0.04, 1) + 1e-12) return("fail")
    "pass"
}

#' Simulated per-unit values: normal (mean = median = 0.01) or skewed lognormal of log-SD 1
#' centred at its median (`c = 1`) or its mean (`c = exp(1/2)`)
#' @noRd
valGenUnits <- function(n, shape, centre) {
    if (shape == "normal") return(0.01 + 0.02 * stats::rnorm(n))
    cc <- if (centre == "median") 1 else exp(0.5)
    0.02 * (exp(stats::rnorm(n)) - cc)
}

#' Coverage rows of one cell and statistic from per-simulation hits (`sims x levels`)
#' @noRd
valCoverageRows <- function(cell, statistic, class, hits, kind, achieved = NULL,
                            undefined = NULL, target = NA_real_, zero_width = NULL,
                            target_se = NA_real_, range = NULL, effect = NA_real_,
                            half_width = NULL, estimate = NULL) {
    nsim <- nrow(hits)
    # r1.14 (4.4) and O4 of VAL-109: the mean half-width over the simulated SD of the statistic
    # times the level's normal quantile. For a jackknife interval (critical value t(g - 1)) with an
    # exactly estimated SE this is t(g - 1) / z, not 1; `se_inflation` removes the t factor (the
    # mean jackknife SE over the simulated SD; 1.00 for an exactly estimated SE; E8 of the re-audit)
    hw_ratio <- rep(NA_real_, ncol(hits))
    se_infl <- rep(NA_real_, ncol(hits))
    if (!is.null(half_width) && !is.null(estimate) && sum(is.finite(estimate)) > 1L) {
        z <- stats::qnorm(1 - (1 - COVERAGE_LEVELS) / 2)
        sd_est <- stats::sd(estimate, na.rm = TRUE)
        hw_ratio <- colMeans(half_width, na.rm = TRUE) / (sd_est * z)
        tq <- stats::qt(1 - (1 - COVERAGE_LEVELS) / 2, as.integer(cell$n) - 1L)
        se_infl <- colMeans(half_width, na.rm = TRUE) / (sd_est * tq)
    }
    cov <- colMeans(hits)
    ach <- if (is.null(achieved)) rep(NA_real_, ncol(hits)) else colMeans(achieved, na.rm = TRUE)
    zw <- if (is.null(zero_width)) rep(NA_real_, ncol(hits)) else colMeans(zero_width)
    # r1.13 (VAL-69): `degenerate` when more than 1 % of the intervals have zero width, or the
    # target lies within one Monte Carlo SE of a bound of the statistic's range
    at_bound <- !is.null(range) && is.finite(target) && is.finite(target_se) &&
        min(abs(target - range)) <= target_se
    degenerate <- (is.finite(zw) & zw > 0.01) | at_bound
    data.frame(
        cell = cell$cell, regime = cell$regime, family = cell$family, n = cell$n,
        condition = cell$condition, statistic = statistic, class = class,
        level = COVERAGE_LEVELS, nsim = nsim, coverage = cov,
        mc_se = sqrt(cov * (1 - cov) / nsim), achieved = ach,
        undefined = if (is.null(undefined)) 0 else mean(undefined),
        target = target,
        status = vapply(seq_along(cov), function(i) {
            if (degenerate[[i]]) return("degenerate")
            valCoverageStatus(kind, cov[[i]], COVERAGE_LEVELS[[i]], ach[[i]], nsim)
        }, character(1)),
        kind = kind, zero_width = zw, effect = effect, half_width_ratio = hw_ratio,
        se_inflation = se_infl, stringsAsFactors = FALSE
    )
}

#' The cell's index (the numeric suffix of its id)
#' @noRd
valCellIndex <- function(cell) as.integer(sub("^.*-", "", cell$cell))

#' Seed of simulation `s` of a cell: fixed by the cell and `s` alone (never by the chunking)
#'
#' @param cell A row of [valCoverageCells()].
#' @param s Simulation index (1 to 99,999).
#' @return An integer seed.
#' @export
valSimSeed <- function(cell, s) 1000000000L + valCellIndex(cell) * 100000L + as.integer(s)

#' Seed of target simulation `t` of a cell (1 to 99,999)
#'
#' @inheritParams valSimSeed
#' @param t Target simulation index.
#' @return An integer seed.
#' @export
valTargetSeed <- function(cell, t) 1500000000L + valCellIndex(cell) * 100000L + as.integer(t)

#' Seed of the bootstrap inside simulation `s` (C7)
#' @noRd
valBootSeed <- function(cell, s) 2000000000L + valCellIndex(cell) * 50000L + as.integer(s)

#' The simulation spec of a harness cell: statistics, their pass band, one simulation, one target
#' simulation (or an analytic target)
#' @noRd
valCovSpec <- function(cell) {
    fam <- strsplit(cell$family, ":", fixed = TRUE)[[1L]]
    switch(fam[[1L]],
           exact_median = valSpecUnits(cell, "exact_median"),
           paired_t = valSpecUnits(cell, "paired_t"),
           wilson = valSpecBinomial(cell, "wilson"),
           clopper_pearson = valSpecBinomial(cell, "clopper_pearson"),
           newcombe_paired = valSpecNewcombe(cell),
           c3 = valSpecC3(cell, fam),
           c4 = valSpecMatrix(cell, "replicate", as.integer(fam[[2L]])),
           c4e = valSpecMatrix(cell, "block", as.integer(sub("(acc|p1)$", "", fam[[2L]])),
                               rank = fam[[2L]] %in% c("4840", "4840p1"),
                               contrast = if (identical(fam[[2L]], "4840p1")) 2L),
           c7 = valSpecBootstrap(cell, fam[[2L]]),
           valAbort("unknown coverage family {cell$family}", "val_spec_invalid"))
}

#' Interval bounds of one statistic at every level, as a list
#' @noRd
valBounds <- function(intervals) {
    list(lo = vapply(intervals, function(iv) iv$ci_low %||% NA_real_, numeric(1)),
         hi = vapply(intervals, function(iv) iv$ci_high %||% NA_real_, numeric(1)),
         achieved = vapply(intervals, function(iv) iv$achieved_level %||% NA_real_, numeric(1)),
         undefined = any(vapply(intervals, function(iv) !identical(iv$status, "ok"), logical(1))))
}

#' C1 and C2: per-unit values (analytic target: 0.01 normal; 0 skewed, its median or mean)
#' @noRd
valSpecUnits <- function(cell, method) {
    centre <- if (method == "exact_median") "median" else "mean"
    list(stats = data.frame(statistic = method, class = NA_character_,
                            kind = if (method == "exact_median") "exact" else "approx_run"),
         target = if (cell$condition == "normal") 0.01 else 0,
         sim = function(s) {
             d <- valGenUnits(cell$n, cell$condition, centre)
             b <- valBounds(lapply(COVERAGE_LEVELS, function(l) {
                 valInterval(d, method = method, conf = l)
             }))
             list(bounds = list(b), data = d)
         })
}

#' C5 and C8: binomial shares
#' @noRd
valSpecBinomial <- function(cell, method) {
    p <- as.numeric(sub("^p", "", cell$condition))
    list(stats = data.frame(statistic = method, class = NA_character_,
                            kind = if (method == "wilson") "approx_within" else "clopper"),
         target = p,
         sim = function(s) {
             x <- stats::rbinom(1L, cell$n, p)
             b <- valBounds(lapply(COVERAGE_LEVELS, function(l) {
                 valInterval(x, n = cell$n, method = method, conf = l)
             }))
             list(bounds = list(b), data = x)
         })
}

#' C6: correlated binary pairs (latent correlation 0.5)
#' @noRd
valSpecNewcombe <- function(cell) {
    pp <- as.numeric(strsplit(sub("^p", "", cell$condition), "_", fixed = TRUE)[[1L]])
    list(stats = data.frame(statistic = "newcombe_paired", class = NA_character_,
                            kind = "approx_within"),
         target = pp[[1L]] - pp[[2L]],
         sim = function(s) {
             z1 <- stats::rnorm(cell$n)
             z2 <- 0.5 * z1 + sqrt(1 - 0.25) * stats::rnorm(cell$n)
             x1 <- z1 < stats::qnorm(pp[[1L]])
             x2 <- z2 < stats::qnorm(pp[[2L]])
             tab <- c(both = sum(x1 & x2), first_only = sum(x1 & !x2),
                      second_only = sum(!x1 & x2), neither = sum(!x1 & !x2))
             b <- valBounds(lapply(COVERAGE_LEVELS, function(l) {
                 valInterval(tab, method = "newcombe_paired", conf = l)
             }))
             list(bounds = list(b), data = tab)
         })
}

#' C3: per-design values over units; theta = mean over designs, deletion j = mean over the designs
#' not using unit j
#' @noRd
valSpecC3 <- function(cell, fam) {
    cfg <- valC3Configs()[[fam[[2L]]]]
    cls <- fam[[3L]]
    tr <- fam[[4L]]
    shape <- sub("_sd.*$", "", cell$condition)
    tau <- as.numeric(sub("^.*_sd", "", cell$condition))
    n_batch <- max(cfg$designs)
    uses <- vapply(cfg$units, function(u) apply(cfg$designs, 1L, function(d) any(d %in% u)),
                   logical(nrow(cfg$designs)))
    incidence <- t(apply(cfg$designs, 1L, function(d) {
        x <- numeric(n_batch)
        x[d] <- 1 / length(d)
        x
    }))
    oneSim <- function() {
        u <- if (shape == "normal") {
            tau * stats::rnorm(n_batch)
        } else {
            tau * (stats::rexp(n_batch) - 1)
        }
        m <- drop(incidence %*% u) + 0.02 * stats::rnorm(nrow(incidence))
        if (cls == "accuracy") return(0.1 + m)
        p <- pmin(pmax(0.05 + m, 0), 1)
        stats::rbinom(length(p), 500L, p) / 500
    }
    list(stats = data.frame(statistic = paste(fam[-1L], collapse = ":"), class = cls,
                            kind = "approx_run"),
         target_one = function() mean(oneSim()),
         sim = function(s) {
             v <- oneSim()
             dels <- vapply(seq_len(ncol(uses)), function(j) mean(v[!uses[, j]]), numeric(1))
             b <- valBounds(lapply(COVERAGE_LEVELS, function(l) {
                 valJackknife(full = mean(v), deletions = dels, conf = l, form = "unit",
                              transform = if (tr == "log") "log", class = cls)
             }))
             list(bounds = list(b), data = v)
         })
}

#' C4 and C4e: jackknife forms (ii) and (iii) on simulated feature matrices; `rank` the rank
#' statistic at the calibrated scale (the 4,840-feature families), otherwise `accuracy` and
#' `threshold_count` at scale 1 (DI10K's `c4e:4840acc` too); `contrast` a single contrast column
#' for the rank statistic and its calibration (P1K's `c4e:4840p1`), otherwise all contrasts pooled
#' @noRd
valSpecMatrix <- function(cell, form, n_feat, rank = n_feat == 4840L, contrast = NULL) {
    g <- cell$n
    lay <- valCovLayout(form, g)
    # r1.13 (VAL-69): the rank cells' effect is calibrated per g into the declared pAUC band
    scale <- if (rank) valRankEffect(form, g, contrast = contrast)$effect else 1
    # I10 diagnostic (VAL-109): `fixed_symmetric` draws the feature baselines once per cell and
    # splits the changed features +1 and -1; `fixed_onedir` (DI10b, VAL-137 (1)) keeps that draw
    # with every change +1, as T2's; `redraw` (T2) draws them in every simulation
    fixed <- if (isTRUE(cell[["generator"]] %in% c("fixed_symmetric", "fixed_onedir"))) {
        valCovFixedFeatures(cell, n_feat)
    }
    statsOf <- function(x, l, truth) {
        valCovStats(valCovFit(x, l, form), truth, n_feat, form, rank = rank,
                    contrast = contrast)
    }
    names_s <- if (rank) {
        "rank.pauc_diff"
    } else if (form == "replicate") {
        c("accuracy.mae", "accuracy.mae_diff", "threshold_count.sens",
          "threshold_count.sens_diff", "threshold_count.fcs", "threshold_count.fcs_diff")
    } else {
        c("accuracy.mae", "accuracy.bias", "threshold_count.rej")
    }
    parts <- strsplit(names_s, ".", fixed = TRUE)
    list(stats = data.frame(statistic = vapply(parts, `[`, "", 2L),
                            class = vapply(parts, `[`, "", 1L), kind = "approx_run"),
         range = if (rank) c(-1, 1),
         effect = if (rank) scale else NA_real_,
         target_one = function() {
             dat <- valCovMatrixData(lay, n_feat, form, scale, fixed = fixed)
             unname(statsOf(dat$x, lay, dat$truth))
         },
         sim = function(s) {
             dat <- valCovMatrixData(lay, n_feat, form, scale, fixed = fixed)
             full <- unname(statsOf(dat$x, lay, dat$truth))
             dels <- matrix(vapply(seq_len(g), function(u) {
                 keep <- lay$unit != u
                 unname(statsOf(dat$x[, keep], lay[keep, ], dat$truth))
             }, numeric(length(names_s))), nrow = length(names_s))
             bounds <- lapply(seq_along(names_s), function(k) {
                 valBounds(lapply(COVERAGE_LEVELS, function(l) {
                     valJackknife(full = full[[k]], deletions = dels[k, ], conf = l, form = form)
                 }))
             })
             list(bounds = bounds, data = dat$x, estimate = full)
         })
}

#' C7: the paired within-run bootstrap over 200 keys in two classes, iid within class
#' @noRd
valSpecBootstrap <- function(cell, cls) {
    simData <- function() {
        cl <- rep(c("x", "y"), each = 100L)
        if (cls == "accuracy") {
            e_a <- stats::rnorm(200L, 0, ifelse(cl == "x", 0.3, 0.5))
            return(data.frame(class = cl, a = e_a, b = 0.8 * e_a + stats::rnorm(200L, 0, 0.15)))
        }
        z1 <- stats::rnorm(200L)
        z2 <- 0.5 * z1 + sqrt(0.75) * stats::rnorm(200L)
        pa <- ifelse(cl == "x", 0.7, 0.05)
        pb <- ifelse(cl == "x", 0.6, 0.04)
        data.frame(class = cl, a = as.numeric(z1 < stats::qnorm(pa)),
                   b = as.numeric(z2 < stats::qnorm(pb)))
    }
    statistic <- if (cls == "accuracy") {
        function(d) stats::median(abs(d$b)) - stats::median(abs(d$a))
    } else {
        function(d) {
            x <- d$class == "x"
            mean(d$a[x]) - mean(d$b[x])
        }
    }
    b_max <- valBootstrapMinB(max(COVERAGE_LEVELS))
    list(stats = data.frame(statistic = paste0("paired_", cls), class = cls,
                            kind = "approx_within"),
         target_one = function() statistic(simData()),
         sim = function(s) {
             d <- simData()
             bs <- valBootstrap(d, statistic, B = b_max, conf = max(COVERAGE_LEVELS),
                                seed = valBootSeed(cell, s), strata = "class",
                                keep_replicates = TRUE)
             reps <- bs$replicates[is.finite(bs$replicates)]
             q <- lapply(COVERAGE_LEVELS, function(l) {
                 tail <- (1 - l) / 2
                 stats::quantile(reps, c(tail, 1 - tail), type = 6, names = FALSE)
             })
             b <- list(lo = vapply(q, `[`, 0, 1L), hi = vapply(q, `[`, 0, 2L),
                       achieved = rep(NA_real_, length(COVERAGE_LEVELS)),
                       undefined = !identical(bs$status, "ok"))
             list(bounds = list(b), data = as.matrix(d[, c("a", "b")]))
         })
}

#' Run a range of a cell's simulations and target simulations
#'
#' Each simulation `s` runs under [valSimSeed()] and each target simulation `t` under
#' [valTargetSeed()], so any split of a cell into ranges gives the same data sets and bounds as the
#' cell run whole (coordinator's ruling of 2026-11-23, point 3).
#'
#' @param cell A row of [valCoverageCells()].
#' @param sims Simulation indices.
#' @param targets Target simulation indices (empty for a cell with an analytic target).
#' @return A list: the cell id, `sims`, bound arrays (`lo`, `hi`, `achieved`: sims x levels x
#'   statistics), `undefined` (sims x statistics), per-simulation data digests, `targets`,
#'   `target_sum`, and the seconds of each phase.
#' @export
valCoverageRange <- function(cell, sims, targets = integer()) {
    cell <- as.list(cell)
    spec <- valCovSpec(cell)
    n_stat <- nrow(spec$stats)
    n_lv <- length(COVERAGE_LEVELS)
    lo <- hi <- ach <- array(NA_real_, c(length(sims), n_lv, n_stat))
    undefined <- matrix(FALSE, length(sims), n_stat)
    estimate <- matrix(NA_real_, length(sims), n_stat)
    digests <- character(length(sims))
    t0 <- proc.time()[["elapsed"]]
    for (i in seq_along(sims)) {
        r <- valWithSeed(valSimSeed(cell, sims[[i]]), spec$sim(sims[[i]]))
        for (k in seq_len(n_stat)) {
            b <- r$bounds[[k]]
            lo[i, , k] <- b$lo
            hi[i, , k] <- b$hi
            ach[i, , k] <- b$achieved
            undefined[i, k] <- b$undefined
        }
        if (!is.null(r$estimate)) estimate[i, ] <- r$estimate
        digests[[i]] <- digest::digest(r$data, algo = "sha256")
    }
    t1 <- proc.time()[["elapsed"]]
    if (is.null(spec$target_one)) targets <- integer()
    # Each target simulation's value is kept (targets x statistics), so that the combined target is
    # summed in index order whatever the split (a sum of range sums differs in the last bits).
    target_values <- matrix(NA_real_, length(targets), n_stat)
    for (j in seq_along(targets)) {
        target_values[j, ] <- valWithSeed(valTargetSeed(cell, targets[[j]]), spec$target_one())
    }
    list(cell = cell$cell, sims = as.integer(sims), lo = lo, hi = hi, achieved = ach,
         undefined = undefined, estimate = estimate, digests = digests,
         targets = as.integer(targets),
         target_values = target_values, target_sum = colSums(target_values),
         seconds_sims = t1 - t0,
         seconds_target = proc.time()[["elapsed"]] - t1)
}

#' Combine a cell's ranges into its coverage rows
#'
#' Refuses (`val_spec_invalid`) unless the ranges hold every simulation 1..`nsim` and (for a cell
#' with a simulated target) every target simulation 1..`target_sims` exactly once.
#'
#' @param cell A row of [valCoverageCells()].
#' @param parts Results of [valCoverageRange()].
#' @param nsim,target_sims The cell's simulation counts.
#' @return Coverage rows, one per statistic and level, with attribute `data_digests`.
#' @export
valCombineCoverage <- function(cell, parts, nsim, target_sims) {
    cell <- as.list(cell)
    spec <- valCovSpec(cell)
    sims <- unlist(lapply(parts, `[[`, "sims"))
    if (!identical(sort(sims), seq_len(nsim))) {
        valAbort("the ranges of {cell$cell} do not hold each simulation exactly once",
                 "val_spec_invalid")
    }
    target <- spec[["target"]]
    if (is.null(target)) {
        targets <- unlist(lapply(parts, `[[`, "targets"))
        if (!identical(sort(targets), seq_len(target_sims))) {
            valAbort("the ranges of {cell$cell} do not hold each target simulation exactly once",
                     "val_spec_invalid")
        }
        values <- do.call(rbind, lapply(parts, `[[`, "target_values"))[order(targets), ,
                                                                       drop = FALSE]
        total <- numeric(ncol(values))
        for (j in seq_len(nrow(values))) total <- total + values[j, ]
        target <- total / target_sims
        target_se <- apply(values, 2L, stats::sd) / sqrt(target_sims)
    } else {
        target_se <- rep(NA_real_, nrow(spec$stats))
    }
    ord <- order(sims)
    bind <- function(field) {
        a <- do.call(abind3, lapply(parts, `[[`, field))
        a[ord, , , drop = FALSE]
    }
    lo <- bind("lo")
    hi <- bind("hi")
    ach <- bind("achieved")
    undefined <- do.call(rbind, lapply(parts, `[[`, "undefined"))[ord, , drop = FALSE]
    est <- do.call(rbind, lapply(parts, function(p) {
        p[["estimate"]] %||% matrix(NA_real_, length(p$sims), nrow(spec$stats))
    }))[ord, , drop = FALSE]
    out <- do.call(rbind, lapply(seq_len(nrow(spec$stats)), function(k) {
        hits <- lo[, , k, drop = FALSE] <= target[[k]] & target[[k]] <= hi[, , k, drop = FALSE]
        hits <- matrix(hits %in% TRUE, nrow = nsim)
        st <- spec$stats[k, ]
        zero <- matrix((lo[, , k] == hi[, , k]) %in% TRUE, nrow = nsim)
        # R-P1 (C66; KIT-162 6.7): misses at the lower bound (the target below the interval) and
        # at the upper bound (the target above it), as counts of simulations per level
        low <- matrix((lo[, , k] > target[[k]]) %in% TRUE, nrow = nsim)
        high <- matrix((hi[, , k] < target[[k]]) %in% TRUE, nrow = nsim)
        rows <- valCoverageRows(cell, st$statistic, st$class, hits, st$kind,
                                achieved = if (st$kind == "exact") matrix(ach[, , k], nrow = nsim),
                                undefined = undefined[, k], target = target[[k]], zero_width = zero,
                                target_se = target_se[[k]], range = spec[["range"]],
                                effect = spec[["effect"]] %||% NA_real_,
                                half_width = matrix((hi[, , k] - lo[, , k]) / 2, nrow = nsim),
                                estimate = est[, k])
        rows$miss_low <- as.integer(colSums(low))
        rows$miss_high <- as.integer(colSums(high))
        rows
    }))
    out$generator <- cell[["generator"]] %||% "redraw"
    out$seconds_sims <- sum(vapply(parts, `[[`, 0, "seconds_sims"))
    out$seconds_target <- sum(vapply(parts, `[[`, 0, "seconds_target"))
    attr(out, "data_digests") <- unlist(lapply(parts, `[[`, "digests"))[ord]
    out
}

#' Bind 3-d arrays along their first dimension
#' @noRd
abind3 <- function(...) {
    xs <- list(...)
    d <- dim(xs[[1L]])
    n <- sum(vapply(xs, function(x) dim(x)[[1L]], 0))
    out <- array(NA_real_, c(n, d[[2L]], d[[3L]]))
    at <- 0L
    for (x in xs) {
        m <- dim(x)[[1L]]
        if (m > 0L) out[at + seq_len(m), , ] <- x
        at <- at + m
    }
    out
}

#' Run one coverage regime cell whole, or a ticket regime (section 4.4)
#'
#' A harness cell is a row of [valCoverageCells()]; its simulations use the per-index seeds of
#' [valSimSeed()] and [valTargetSeed()]. A ticket regime is a list with `id`, `method`, `class`,
#' `form`, `g`, `levels`, `simulations`, a function `simulate(nsim, seed, levels)` returning a data
#' frame with `level`, `covered` (and optionally `achieved`), and `kind` of its pass band
#' (`approx_run` by default).
#'
#' @param regime A cell (one-row data frame) or a ticket regime (list).
#' @param nsim Simulations (at least 2,000 in T2 and for a ticket regime).
#' @param seed Seed of a ticket regime (harness cells take theirs from the cell).
#' @param target_sims Simulations of the full statistic for the target (20,000; fewer only when
#'   probing cost).
#' @return A data frame of coverage rows, one per level (and statistic).
#' @export
valCoverageRegime <- function(regime, nsim = 2000L, seed = NULL, target_sims = TARGET_SIMS) {
    if (is.list(regime) && !is.data.frame(regime)) return(valTicketRegime(regime, nsim, seed))
    part <- valCoverageRange(regime, seq_len(nsim), seq_len(target_sims))
    valCombineCoverage(regime, list(part), nsim, target_sims)
}

#' A ticket coverage regime (declared in `coverage_regimes`)
#' @noRd
valTicketRegime <- function(regime, nsim, seed) {
    checkmate::assert_function(regime$simulate)
    if (nsim < 2000L) {
        valAbort("a ticket regime needs at least 2,000 simulations", "val_spec_invalid")
    }
    res <- regime$simulate(nsim, seed, regime$levels)
    checkmate::assert_names(names(res), must.include = c("level", "covered"))
    kind <- regime[["kind"]] %||% "approx_run"
    do.call(rbind, lapply(regime$levels, function(lv) {
        r <- res[abs(res$level - lv) < 5e-5, , drop = FALSE]
        cov <- mean(r$covered)
        ach <- if (is.null(r[["achieved"]])) NA_real_ else mean(r[["achieved"]])
        data.frame(regime = regime$id, method = regime$method, class = regime$class,
                   form = regime$form, n = regime$g, level = lv, nsim = nrow(r), coverage = cov,
                   mc_se = sqrt(cov * (1 - cov) / nrow(r)), achieved = ach,
                   status = valCoverageStatus(kind, cov, lv, ach, nrow(r)), source = "ticket",
                   stringsAsFactors = FALSE)
    }))
}

RANK_EFFECT_GRID <- seq(0.05, 1, by = 0.05)
RANK_CACHE <- new.env(parent = emptyenv())

#' The rank cells' effect (KIT-159 10.2, r1.13, VAL-69)
#'
#' The generator's group effects are scaled by the grid value whose mean normalised pAUC at FPR 0.01
#' of the plain (moderated) arm, over `nsim` seeded calibration simulations, lies in [0.4, 0.8] and
#' is nearest 0.6; the grid is walked upwards and stops once the mean passes 0.8. Seeds depend on
#' the form, g and the simulation index only, so every process finds the same effect. No grid value
#' in the band is `val_spec_invalid`. With `contrast` (P1K, VAL-176) the same rule is applied to the
#' pAUC of that one contrast column, on the same calibration data sets.
#'
#' @param form `replicate` (C4) or `block` (C4e).
#' @param g Unit count.
#' @param nsim Calibration simulations per grid value.
#' @param contrast `NULL` (every contrast pooled) or one contrast column.
#' @return A list: `effect`, `mean_pauc`.
#' @export
valRankEffect <- function(form, g, nsim = 50L, contrast = NULL) {
    key <- paste(c(form, g, nsim, contrast), collapse = " ")
    if (!is.null(RANK_CACHE[[key]])) return(RANK_CACHE[[key]])
    lay <- valCovLayout(form, g)
    base <- 1900000000L + (if (form == "block") 50000L else 0L) + as.integer(g) * 1000L
    means <- numeric()
    for (scale in RANK_EFFECT_GRID) {
        vals <- vapply(seq_len(nsim), function(s) {
            valWithSeed(base + s, {
                dat <- valCovMatrixData(lay, 4840L, form, scale)
                fit <- valCovFit(dat$x, lay, form)
                truth <- dat$truth
                if (!is.null(contrast)) {
                    truth <- truth[, contrast]
                    fit$est <- fit$est[, contrast]
                    fit$p <- fit$p[, contrast]
                }
                sc <- data.frame(key = seq_along(truth), expected = as.vector(truth),
                                 estimate = as.vector(fit$est), p = as.vector(fit$p),
                                 s_own = as.vector(fit$p), truth_status = "scored",
                                 key_status = "ok", stringsAsFactors = FALSE)
                valPartialAuc(sc, fpr_max = 0.01)$value
            })
        }, numeric(1))
        means <- c(means, mean(vals))
        if (mean(vals) > 0.8) break
    }
    ok <- means >= 0.4 & means <= 0.8
    if (!any(ok)) {
        valAbort("no rank effect on the grid gives a mean pAUC in [0.4, 0.8]", "val_spec_invalid")
    }
    i <- which(ok)[which.min(abs(means[ok] - 0.6))]
    res <- list(effect = RANK_EFFECT_GRID[[i]], mean_pauc = means[[i]])
    RANK_CACHE[[key]] <- res
    res
}

#' Sample layout of C4 (two groups, g replicates; form `replicate`) and C4e (form `block`: g = 5
#' blocks of two samples of each of five groups, KIT-162; otherwise g blocks of two samples of each
#' of two groups)
#' @noRd
valCovLayout <- function(form, g) {
    if (form == "replicate") {
        return(data.frame(group = rep(c("A", "B"), each = g), unit = rep(seq_len(g), 2L),
                          stringsAsFactors = FALSE))
    }
    n_groups <- if (g == 5L) 5L else 2L
    lay <- expand.grid(rep = 1:2, group = LETTERS[seq_len(n_groups)], unit = seq_len(g),
                       stringsAsFactors = FALSE)
    lay[, c("group", "unit")]
}

#' A simulated log2 matrix (section 10.2, feature matrices): 20 % of 500 features (or 40 of 4,840)
#' changed by +1 log2 per step of the group effects, residual SD 0.3, per-sample offset SD 0.1, and
#' for form `block` a per-feature block effect of SD 0.1
#' @noRd
valCovMatrixData <- function(lay, n_feat, form, scale = 1, fixed = NULL) {
    n_changed <- if (n_feat == 4840L) 40L else as.integer(n_feat / 5)
    groups <- sort(unique(lay$group))
    effect <- scale * (if (length(groups) == 5L) c(0, 0.5, 1, 1.5, 2) else c(0, 1))
    changed <- seq_len(n_feat) <= n_changed
    if (!is.null(fixed)) changed <- fixed$sign
    gi <- match(lay$group, groups)
    base <- if (is.null(fixed)) stats::rnorm(n_feat, 20, 1) else fixed$baseline
    x <- matrix(base, n_feat, nrow(lay)) +
        outer(changed, effect[gi]) +
        matrix(stats::rnorm(nrow(lay), 0, 0.1), n_feat, nrow(lay), byrow = TRUE) +
        matrix(stats::rnorm(n_feat * nrow(lay), 0, 0.3), n_feat, nrow(lay))
    if (form == "block") {
        b <- matrix(stats::rnorm(n_feat * max(lay$unit), 0, 0.1), n_feat, max(lay$unit))
        x <- x + b[, lay$unit]
    }
    truth <- outer(changed, effect[-1L] - effect[[1L]])
    list(x = x, truth = truth)
}

#' The I10 diagnostic generator's feature draw (VAL-109): one seeded draw per cell of the feature
#' baselines (normal, mean 20, SD 1), and the changed features (the first fifth, or 40 of 4,840)
#' split alternately +1 and -1, so that median normalisation is not shifted by the changes and the
#' target is conditional on the features, as a fixed benchmark's is. Generator `fixed_onedir`
#' (DI10b): the same draw with every changed feature +1, as T2's and the spike-in benchmarks'
#' one-directional changes
#' @noRd
valCovFixedFeatures <- function(cell, n_feat) {
    n_changed <- if (n_feat == 4840L) 40L else as.integer(n_feat / 5)
    baseline <- valWithSeed(2100000000L + valCellIndex(cell), stats::rnorm(n_feat, 20, 1))
    step <- if (identical(cell[["generator"]], "fixed_onedir")) 1 else c(1, -1)
    sign <- c(rep_len(step, n_changed), rep(0, n_feat - n_changed))
    list(baseline = baseline, sign = sign)
}

#' A plan row's cell, with the plan's `generator` when the plan carries one (I10)
#' @noRd
valPlanCell <- function(cells, p) {
    cell <- cells[cells$cell == p$cell, ]
    if (!is.null(p[["generator"]]) && !is.na(p[["generator"]])) cell$generator <- p[["generator"]]
    cell
}

#' Two arms on one matrix: median normalisation, limma with group (and block) terms, contrasts of
#' each group against the first. `plain`: moderated p; `ordinary`: unmoderated t p
#' @noRd
valCovFit <- function(x, lay, form) {
    med <- apply(x, 2L, stats::median)
    x <- sweep(x, 2L, med) + stats::median(med)
    grp <- factor(lay$group)
    des <- if (form == "block") {
        stats::model.matrix(~ 0 + grp + factor(lay$unit))
    } else {
        stats::model.matrix(~ 0 + grp)
    }
    k <- nlevels(grp)
    cm <- matrix(0, ncol(des), k - 1L)
    cm[1L, ] <- -1
    for (j in seq_len(k - 1L)) cm[j + 1L, j] <- 1
    fit <- limma::contrasts.fit(limma::lmFit(x, des), cm)
    t_ord <- fit$coefficients / (fit$stdev.unscaled * fit$sigma)
    p_ord <- 2 * stats::pt(-abs(t_ord), fit$df.residual)
    eb <- limma::eBayes(fit)
    list(est = eb$coefficients, p = eb$p.value, p_ord = p_ord)
}

#' The statistics of C4 (accuracy, threshold_count on 500 features; rank on 4,840) and C4e
#' @noRd
valCovStats <- function(fit, truth, n_feat, form, rank = n_feat == 4840L, contrast = NULL) {
    if (!is.null(contrast)) {
        truth <- truth[, contrast]
        fit <- lapply(fit[c("est", "p", "p_ord")], function(m) m[, contrast])
    }
    bh <- function(p) apply(p, 2L, stats::p.adjust, method = "BH")
    pauc <- function(p) {
        s <- data.frame(key = seq_along(truth), expected = as.vector(truth),
                        estimate = as.vector(fit$est), p = as.vector(p), s_own = as.vector(p),
                        truth_status = "scored", key_status = "ok", stringsAsFactors = FALSE)
        valPartialAuc(s, fpr_max = 0.01)$value
    }
    if (rank) return(c(rank.pauc_diff = pauc(fit$p) - pauc(fit$p_ord)))
    pos <- truth != 0
    callp <- bh(fit$p) < 0.05
    callo <- bh(fit$p_ord) < 0.05
    if (form == "replicate") {
        return(c(
            accuracy.mae = stats::median(abs(fit$est - truth)),
            accuracy.mae_diff = stats::median(abs(fit$est - truth)) -
                stats::median(abs(0.8 * fit$est - truth)),
            threshold_count.sens = mean(callp[pos]),
            threshold_count.sens_diff = mean(callp[pos]) - mean(callo[pos]),
            threshold_count.fcs = mean(callp[!pos]),
            threshold_count.fcs_diff = mean(callp[!pos]) - mean(callo[!pos])
        ))
    }
    c(accuracy.mae = stats::median(abs(fit$est - truth)),
      accuracy.bias = stats::median(fit$est - truth),
      threshold_count.rej = mean(fit$p[!pos] < 0.05))
}

#' The regime table from coverage rows (section 10.2)
#'
#' One row per method x class x form x design x unit count x level. Over the data conditions and
#' statistics of a row: any `fail` fails it, except that a paired t count whose normal case passes
#' and skewed case fails is `flagged` (C2); otherwise any `conservative` makes it conservative;
#' `not_offered` where the method gives no interval at that count and level. Each row carries
#' `coverage_min` and `coverage_max` over its conditions and statistics and, for a row not offered
#' because a condition failed its band, `failed_condition` and `failed_side` (`low`, `high`; for
#' an exact method against its achieved level) (5.7, r1.12).
#'
#' @param detail Rows of [valCoverageRegime()].
#' @return The regime table (`coverage-table.tsv`).
#' @export
valSummariseCoverage <- function(detail) {
    method <- c(C1 = "exact_median", C2 = "paired_t", C3 = "jackknife", C4 = "jackknife",
                C4e = "jackknife", C5 = "wilson", C6 = "newcombe_paired",
                C7 = "bootstrap_within_run", C8 = "clopper_pearson")
    form <- c(C3 = "unit", C4 = "replicate", C4e = "block")
    d <- detail
    d$method <- method[d$regime]
    d$form <- unname(form[d$regime])
    d$design <- ifelse(d$regime == "C3", vapply(strsplit(d$family, ":"), `[`, "", 2L), NA)
    d$transform <- ifelse(d$regime == "C3", vapply(strsplit(d$family, ":"), function(x) {
        x[length(x)]
    }, ""), NA)
    key <- paste(d$regime, d$method, d$class, d$form, d$design, d$transform, d$n,
                 sprintf("%.6f", d$level), sep = "|")
    rows <- lapply(split(seq_len(nrow(d)), key), function(i) {
        r <- d[i, ]
        st <- r$status
        status <- if (r$regime[[1L]] == "C2" && any(st == "fail")) {
            normal_ok <- all(st[r$condition == "normal"] %in% c("pass", "conservative"))
            skewed_bad <- any(st[r$condition == "skewed"] == "fail")
            if (normal_ok && skewed_bad) "flagged" else "fail"
        } else if (any(st == "degenerate")) {
            "degenerate"
        } else if (any(st == "not_offered")) {
            "not_offered"
        } else if (any(st == "fail")) {
            "fail"
        } else if (any(st == "conservative")) {
            "conservative"
        } else {
            "pass"
        }
        w <- which.min(r$coverage)
        bad <- r$status == "fail"
        side <- ifelse(r$kind == "exact", ifelse(r$coverage < r$achieved, "low", "high"),
                       ifelse(r$coverage < r$level, "low", "high"))
        offered <- status %in% c("pass", "conservative")
        data.frame(
            regime = r$regime[[1L]], method = r$method[[1L]], class = r$class[[1L]],
            form = r$form[[1L]], design = r$design[[1L]], transform = r$transform[[1L]],
            n = r$n[[1L]], level = r$level[[1L]], status = status, coverage = r$coverage[[w]],
            mc_se = r$mc_se[[w]], source = "harness",
            statistic = paste(unique(paste(r$statistic, r$condition, sep = "@")), collapse = ", "),
            coverage_min = min(r$coverage), coverage_max = max(r$coverage),
            failed_condition = if (offered || !any(bad)) NA_character_ else
                paste(unique(r$condition[bad]), collapse = ","),
            failed_side = if (offered || !any(bad)) NA_character_ else
                paste(unique(side[bad]), collapse = ","),
            zero_width = if (is.null(r$zero_width)) NA_real_ else max(r$zero_width),
            effect = if (is.null(r$effect)) NA_real_ else r$effect[[1L]],
            half_width_ratio = if (is.null(r$half_width_ratio)) NA_real_ else
                r$half_width_ratio[[w]],
            se_inflation = if (is.null(r$se_inflation)) NA_real_ else r$se_inflation[[w]],
            generator = if (is.null(r$generator)) NA_character_ else r$generator[[1L]],
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    out[order(out$regime, out$method, out$class, out$n, out$level), , drop = FALSE]
}

#' A coverage plan: each cell split into a fixed number of consecutive ranges
#'
#' Range k of a cell holds simulations `ceiling(nsim (k - 1) / r) + 1` to `ceiling(nsim k / r)`, and
#' the same split of its target simulations; the ranges depend only on the declared `ranges`, never
#' on the chunk count.
#'
#' @param cells Rows of [valCoverageCells()].
#' @param ranges Ranges per cell (recycled).
#' @param nsim,target_sims Simulation counts.
#' @param request Request label of each cell (recycled).
#' @return The plan: `task`, `request`, `cell`, `range`, `sim_from`, `sim_to`, `target_from`,
#'   `target_to`.
#' @export
valCoveragePlan <- function(cells, ranges, nsim = 2000L, target_sims = TARGET_SIMS,
                            request = "all") {
    ranges <- rep_len(as.integer(ranges), nrow(cells))
    request <- rep_len(request, nrow(cells))
    cut <- function(n, r) as.integer(ceiling(n * (0:r) / r))
    rows <- lapply(seq_len(nrow(cells)), function(i) {
        r <- ranges[[i]]
        s <- cut(nsim, r)
        t <- cut(target_sims, r)
        data.frame(request = request[[i]], cell = cells$cell[[i]], range = seq_len(r),
                   sim_from = s[-length(s)] + 1L, sim_to = s[-1L],
                   target_from = t[-length(t)] + 1L, target_to = t[-1L],
                   stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    out$task <- sprintf("%s.r%02d", out$cell, out$range)
    out[, c("task", "request", "cell", "range", "sim_from", "sim_to", "target_from",
            "target_to")]
}

#' Run plan tasks in this process (one chunk) and save each task's range result
#'
#' @param plan_file The plan (`.tsv`).
#' @param tasks Task ids to run.
#' @param out_dir Folder for `<task>.rds`.
#' @return `TRUE`, invisibly.
#' @export
valCoverageTasks <- function(plan_file, tasks, out_dir) {
    plan <- utils::read.delim(plan_file, stringsAsFactors = FALSE)
    cells <- valCoverageCells(confirmatory = TRUE)
    for (id in tasks) {
        p <- plan[plan$task == id, ]
        cell <- valPlanCell(cells, p)
        targets <- if (p$target_to >= p$target_from) p$target_from:p$target_to else integer()
        res <- valCoverageRange(cell, p$sim_from:p$sim_to, targets)
        saveRDS(res, file.path(out_dir, paste0(id, ".rds")))
    }
    invisible(TRUE)
}

#' Start chunk processes (task groups), each an `Rscript` with the thread variables at 1; wait
#' @noRd
valRunTaskChunks <- function(groups, plan_file, work_dir) {
    procs <- lapply(seq_along(groups), function(k) {
        code <- sprintf("msrval::valCoverageTasks('%s', c(%s), '%s')", plan_file,
                        paste0("'", groups[[k]], "'", collapse = ", "), work_dir)
        valRscriptProcess(code, file.path(work_dir, sprintf("chunk-%02d.log", k)))
    })
    for (p in procs) p$wait()
    vapply(procs, function(p) p$get_exit_status(), integer(1))
}

#' Deal tasks to chunks by their expected seconds (largest first, to the least loaded chunk)
#' @noRd
valDealTasks <- function(tasks, seconds, chunks) {
    load <- numeric(chunks)
    groups <- vector("list", chunks)
    for (i in order(-seconds)) {
        k <- which.min(load)
        groups[[k]] <- c(groups[[k]], tasks[[i]])
        load[[k]] <- load[[k]] + seconds[[i]]
    }
    groups[lengths(groups) > 0L]
}

#' Run T-COV from a plan in concurrent chunks and write `coverage-detail.tsv`,
#' `coverage-table.tsv`, `coverage-tasks.tsv` and `coverage-run.txt`
#'
#' The dealing of tasks to chunks does not change any result (each simulation has its own seed);
#' a cell whose ranges are incomplete is listed as missing, never scored.
#'
#' @param plan_file The plan (`.tsv`, from [valCoveragePlan()], with `est_seconds`).
#' @param request The request label whose tasks run here (`all` for every task).
#' @param nsim,target_sims Simulation counts the plan was made for.
#' @param chunks Concurrent processes.
#' @param out_dir Output folder.
#' @return The regime table, invisibly.
#' @export
valRunCoverage <- function(plan_file, request = "all", nsim = 2000L, target_sims = TARGET_SIMS,
                           chunks = 4L, out_dir) {
    plan <- utils::read.delim(plan_file, stringsAsFactors = FALSE)
    if (request != "all") plan <- plan[plan$request == request, , drop = FALSE]
    if (nrow(plan) == 0L) valAbort("no task of request {request}", "val_spec_invalid")
    est <- if (is.null(plan$est_seconds)) rep(1, nrow(plan)) else plan$est_seconds
    work_dir <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "val-coverage")
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
    groups <- valDealTasks(plan$task, est, chunks)
    started <- Sys.time()
    status <- valRunTaskChunks(groups, normalizePath(plan_file), work_dir)
    wall <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    cells <- valCoverageCells(confirmatory = TRUE)
    detail <- list()
    tasks_out <- list()
    missing <- character()
    for (cid in unique(plan$cell)) {
        ids <- plan$task[plan$cell == cid]
        files <- file.path(work_dir, paste0(ids, ".rds"))
        if (!all(file.exists(files))) {
            missing <- c(missing, cid)
            next
        }
        parts <- lapply(files, readRDS)
        tasks_out[[cid]] <- data.frame(task = ids, seconds_sims = vapply(parts, `[[`, 0,
                                                                         "seconds_sims"),
                                       seconds_target = vapply(parts, `[[`, 0, "seconds_target"))
        detail[[cid]] <- valCombineCoverage(valPlanCell(cells, plan[plan$cell == cid, ][1L, ]),
                                            parts, nsim, target_sims)
    }
    detail <- do.call(rbind, detail)
    write <- function(x, f) {
        utils::write.table(x, file.path(out_dir, f), sep = "\t", quote = FALSE, row.names = FALSE)
    }
    write(detail, "coverage-detail.tsv")
    write(do.call(rbind, tasks_out), "coverage-tasks.tsv")
    table <- valSummariseCoverage(detail)
    write(table, "coverage-table.tsv")
    writeLines(c(sprintf("request\t%s", request), sprintf("wall_seconds\t%.1f", wall),
                 sprintf("chunks\t%d", length(groups)),
                 sprintf("chunk_exit\t%s", paste(status, collapse = ",")),
                 sprintf("cells_missing\t%s", paste(missing, collapse = ","))),
               file.path(out_dir, "coverage-run.txt"))
    invisible(table)
}

#' Time one regime in one process and in four concurrent chunks (f_sim, G14 (c), section 10.4),
#' and probe the cost of one cell of every family
#'
#' @param id Regime id whose cells are timed (whole cells).
#' @param nsim Simulations per timed cell.
#' @param out_dir Output folder: `fsim.tsv`, `cost-probe.tsv`.
#' @param chunks Concurrent chunks (4).
#' @param probe_nsim,probe_target_sims Sizes of the cost probe.
#' @param probe_cells Cell ids to probe (default: one cell per family).
#' @return f_sim, invisibly.
#' @export
valTimeRegime <- function(id, nsim, out_dir, chunks = 4L, probe_nsim = 20L,
                          probe_target_sims = 200L, probe_cells = NULL) {
    cells <- valCoverageCells()
    f_sim <- NA_real_
    if (!is.null(id) && nzchar(id)) {
        timed <- cells[cells$regime == id, ]
        work <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "val-fsim")
        dir.create(work, recursive = TRUE, showWarnings = FALSE)
        plan <- valCoveragePlan(timed, 1L, nsim = nsim)
        plan_file <- file.path(work, "plan.tsv")
        utils::write.table(plan, plan_file, sep = "\t", quote = FALSE, row.names = FALSE)
        wallOf <- function(groups) {
            t0 <- Sys.time()
            st <- valRunTaskChunks(groups, plan_file, work)
            list(wall = as.numeric(difftime(Sys.time(), t0, units = "secs")), status = st)
        }
        one <- wallOf(list(plan$task))
        four <- wallOf(split(plan$task, rep_len(seq_len(chunks), nrow(plan))))
        f_sim <- one$wall / four$wall
        utils::write.table(data.frame(regime = id, cells = nrow(timed), nsim = nsim,
                                      wall_one_process = round(one$wall, 1),
                                      wall_four_chunks = round(four$wall, 1),
                                      f_sim = round(f_sim, 3),
                                      exit_one = paste(one$status, collapse = ","),
                                      exit_four = paste(four$status, collapse = ","),
                                      cores = valCoreCount()),
                           file.path(out_dir, "fsim.tsv"), sep = "\t", quote = FALSE,
                           row.names = FALSE)
    }
    fams <- if (is.null(probe_cells)) {
        cells[!duplicated(ifelse(cells$regime == "C3", "C3", cells$family)), ]
    } else {
        cells[cells$cell %in% probe_cells, ]
    }
    probe <- do.call(rbind, lapply(seq_len(nrow(fams)), function(i) {
        r <- valCoverageRange(fams[i, ], seq_len(probe_nsim), seq_len(probe_target_sims))
        data.frame(cell = fams$cell[[i]], family = fams$family[[i]], n = fams$n[[i]],
                   probe_nsim = probe_nsim, probe_target_sims = probe_target_sims,
                   seconds_target = r$seconds_target, seconds_sims = r$seconds_sims,
                   stringsAsFactors = FALSE)
    }))
    utils::write.table(probe, file.path(out_dir, "cost-probe.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE)
    invisible(f_sim)
}


#' Cores of the runner (from /proc/cpuinfo)
#' @noRd
valCoreCount <- function() {
    length(grep("^processor", readLines("/proc/cpuinfo")))
}
