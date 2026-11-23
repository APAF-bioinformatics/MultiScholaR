# The coverage regimes T-COV (KIT-159 section 10.2): generators, the regime catalogue, the pass
# bands, the chunked runner and the chunk speed-up on simulations (f_sim, G14 (c)). Every
# simulated data set is scored at every level a consumer declares; levels add no simulations.

COVERAGE_LEVELS <- c(0.95, 0.975, 1 - 0.05 / 3, 1 - 0.05 / 4, 1 - 0.05 / 6, 1 - 0.05 / 8,
                     1 - 0.05 / 12, 1 - 0.05 / 15)
COVERAGE_SEED <- 159200L
TARGET_SIMS <- 20000L
TARGET_SEED <- 159700000L

#' The harness regime cells of T-COV
#'
#' One row per simulation cell: the regime, its method, form, unit count, the data condition and
#' the family of generator. Statistic classes are set per statistic inside a family.
#'
#' @return A data frame of cells with a stable `cell` id.
#' @export
valCoverageCells <- function() {
    cells <- list()
    add <- function(...) cells[[length(cells) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)
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
    if (coverage > min(level + 0.04, 1) + 1e-12) {
        return(if (kind == "approx_run") "conservative" else "fail")
    }
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

#' Coverage rows of one cell from per-simulation hits (`sims x levels` logical matrix)
#' @noRd
valCoverageRows <- function(cell, statistic, class, hits, kind, achieved = NULL,
                            undefined = NULL, target = NA_real_) {
    nsim <- nrow(hits)
    cov <- colMeans(hits)
    ach <- if (is.null(achieved)) rep(NA_real_, ncol(hits)) else colMeans(achieved, na.rm = TRUE)
    data.frame(
        cell = cell$cell, regime = cell$regime, family = cell$family, n = cell$n,
        condition = cell$condition, statistic = statistic, class = class,
        level = COVERAGE_LEVELS, nsim = nsim, coverage = cov, mc_se = sqrt(cov * (1 - cov) / nsim),
        achieved = ach,
        undefined = if (is.null(undefined)) 0 else mean(undefined),
        target = target,
        status = vapply(seq_along(cov), function(i) {
            valCoverageStatus(kind, cov[[i]], COVERAGE_LEVELS[[i]], ach[[i]], nsim)
        }, character(1)),
        kind = kind, stringsAsFactors = FALSE
    )
}

#' Run one coverage regime cell, or a ticket regime (section 4.4)
#'
#' A harness cell is a row of [valCoverageCells()]. A ticket regime is a list with `id`, `method`,
#' `class`, `form`, `g`, `levels`, `simulations`, a function `simulate(nsim, seed, levels)`
#' returning a data frame with `level`, `covered` (and optionally `achieved`), and `kind` of its
#' pass band (`approx_run` by default).
#'
#' @param regime A cell (one-row data frame) or a ticket regime (list).
#' @param nsim Simulations (at least 2,000 in T2 and for a ticket regime).
#' @param seed Seed of the cell.
#' @param target_sims Simulations of the full statistic for the target (20,000; fewer only when
#'   probing cost).
#' @return A data frame of coverage rows, one per level (and statistic).
#' @export
valCoverageRegime <- function(regime, nsim = 2000L, seed = NULL, target_sims = TARGET_SIMS) {
    if (is.list(regime) && !is.data.frame(regime)) return(valTicketRegime(regime, nsim, seed))
    cell <- as.list(regime)
    seed <- seed %||% (COVERAGE_SEED + as.integer(sub("^.*-", "", cell$cell)))
    fam <- strsplit(cell$family, ":", fixed = TRUE)[[1L]]
    t0 <- proc.time()[["elapsed"]]
    out <- switch(fam[[1L]],
           exact_median = valCovUnits(cell, nsim, seed, "exact_median"),
           paired_t = valCovUnits(cell, nsim, seed, "paired_t"),
           c3 = valCovC3(cell, nsim, seed, fam, target_sims),
           c4 = valCovMatrix(cell, nsim, seed, "replicate", as.integer(fam[[2L]]), target_sims),
           c4e = valCovMatrix(cell, nsim, seed, "block", as.integer(fam[[2L]]), target_sims),
           wilson = valCovBinomial(cell, nsim, seed, "wilson"),
           clopper_pearson = valCovBinomial(cell, nsim, seed, "clopper_pearson"),
           newcombe_paired = valCovNewcombe(cell, nsim, seed),
           c7 = valCovBootstrap(cell, nsim, seed, fam[[2L]], target_sims),
           valAbort("unknown coverage family {cell$family}", "val_spec_invalid"))
    if (is.null(out$seconds_sims)) {
        out <- valTimed(out, 0, proc.time()[["elapsed"]] - t0, 0L)
    }
    out
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
    kind <- regime$kind %||% "approx_run"
    do.call(rbind, lapply(regime$levels, function(lv) {
        r <- res[abs(res$level - lv) < 5e-5, , drop = FALSE]
        cov <- mean(r$covered)
        ach <- if (is.null(r$achieved)) NA_real_ else mean(r$achieved)
        data.frame(regime = regime$id, method = regime$method, class = regime$class,
                   form = regime$form, n = regime$g, level = lv, nsim = nrow(r), coverage = cov,
                   mc_se = sqrt(cov * (1 - cov) / nrow(r)), achieved = ach,
                   status = valCoverageStatus(kind, cov, lv, ach, nrow(r)), source = "ticket",
                   stringsAsFactors = FALSE)
    }))
}

#' C1 and C2: per-unit values
#' @noRd
valCovUnits <- function(cell, nsim, seed, method) {
    centre <- if (method == "exact_median") "median" else "mean"
    target <- if (cell$condition == "normal") 0.01 else 0
    hits <- matrix(FALSE, nsim, length(COVERAGE_LEVELS))
    ach <- matrix(NA_real_, nsim, length(COVERAGE_LEVELS))
    valWithSeed(seed, for (s in seq_len(nsim)) {
        d <- valGenUnits(cell$n, cell$condition, centre)
        for (j in seq_along(COVERAGE_LEVELS)) {
            iv <- valInterval(d, method = method, conf = COVERAGE_LEVELS[[j]])
            hits[s, j] <- isTRUE(iv$ci_low <= target && target <= iv$ci_high)
            ach[s, j] <- iv$achieved_level
        }
    })
    kind <- if (method == "exact_median") "exact" else "approx_run"
    out <- valCoverageRows(cell, method, NA_character_, hits, kind,
                           achieved = if (method == "exact_median") ach, target = target)
    if (method == "exact_median") {
        too_few <- is.na(colMeans(ach, na.rm = TRUE))
        out$status[too_few] <- "not_offered"
    }
    out
}

#' C5 and C8: binomial shares
#' @noRd
valCovBinomial <- function(cell, nsim, seed, method) {
    p <- as.numeric(sub("^p", "", cell$condition))
    hits <- matrix(FALSE, nsim, length(COVERAGE_LEVELS))
    valWithSeed(seed, for (s in seq_len(nsim)) {
        x <- stats::rbinom(1L, cell$n, p)
        for (j in seq_along(COVERAGE_LEVELS)) {
            iv <- valInterval(x, n = cell$n, method = method, conf = COVERAGE_LEVELS[[j]])
            hits[s, j] <- iv$ci_low <= p && p <= iv$ci_high
        }
    })
    kind <- if (method == "wilson") "approx_within" else "clopper"
    valCoverageRows(cell, method, NA_character_, hits, kind, target = p)
}

#' C6: correlated binary pairs (latent correlation 0.5)
#' @noRd
valCovNewcombe <- function(cell, nsim, seed) {
    pp <- as.numeric(strsplit(sub("^p", "", cell$condition), "_", fixed = TRUE)[[1L]])
    target <- pp[[1L]] - pp[[2L]]
    hits <- matrix(FALSE, nsim, length(COVERAGE_LEVELS))
    valWithSeed(seed, for (s in seq_len(nsim)) {
        z1 <- stats::rnorm(cell$n)
        z2 <- 0.5 * z1 + sqrt(1 - 0.25) * stats::rnorm(cell$n)
        x1 <- z1 < stats::qnorm(pp[[1L]])
        x2 <- z2 < stats::qnorm(pp[[2L]])
        tab <- c(both = sum(x1 & x2), first_only = sum(x1 & !x2), second_only = sum(!x1 & x2),
                 neither = sum(!x1 & !x2))
        for (j in seq_along(COVERAGE_LEVELS)) {
            iv <- valInterval(tab, method = "newcombe_paired", conf = COVERAGE_LEVELS[[j]])
            hits[s, j] <- iv$ci_low <= target && target <= iv$ci_high
        }
    })
    valCoverageRows(cell, "newcombe_paired", NA_character_, hits, "approx_within", target = target)
}

#' C3: per-design values over units; theta = mean over designs, deletion j = mean over the designs
#' not using unit j
#' @noRd
valCovC3 <- function(cell, nsim, seed, fam, target_sims) {
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
    t0 <- proc.time()[["elapsed"]]
    target <- valWithSeed(TARGET_SEED + seed, mean(replicate(target_sims, mean(oneSim()))))
    t1 <- proc.time()[["elapsed"]]
    hits <- matrix(FALSE, nsim, length(COVERAGE_LEVELS))
    undefined <- logical(nsim)
    valWithSeed(seed, for (s in seq_len(nsim)) {
        v <- oneSim()
        dels <- vapply(seq_len(ncol(uses)), function(j) mean(v[!uses[, j]]), numeric(1))
        for (j in seq_along(COVERAGE_LEVELS)) {
            iv <- valJackknife(full = mean(v), deletions = dels, conf = COVERAGE_LEVELS[[j]],
                               form = "unit", transform = if (tr == "log") "log",
                               class = cls)
            undefined[s] <- iv$status != "ok"
            hits[s, j] <- isTRUE(iv$ci_low <= target && target <= iv$ci_high)
        }
    })
    out <- valCoverageRows(cell, paste(fam[-1L], collapse = ":"), cls, hits, "approx_run",
                           undefined = undefined, target = target)
    valTimed(out, t1 - t0, proc.time()[["elapsed"]] - t1, target_sims)
}

#' Record the seconds of a cell's target phase and simulation phase
#' @noRd
valTimed <- function(out, seconds_target, seconds_sims, target_sims) {
    out$seconds_target <- seconds_target
    out$seconds_sims <- seconds_sims
    out$target_sims <- target_sims
    out
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
valCovMatrixData <- function(lay, n_feat, form) {
    n_changed <- if (n_feat == 4840L) 40L else as.integer(n_feat / 5)
    groups <- sort(unique(lay$group))
    effect <- if (length(groups) == 5L) c(0, 0.5, 1, 1.5, 2) else c(0, 1)
    changed <- seq_len(n_feat) <= n_changed
    gi <- match(lay$group, groups)
    x <- matrix(stats::rnorm(n_feat, 20, 1), n_feat, nrow(lay)) +
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
valCovStats <- function(fit, truth, n_feat, form) {
    bh <- function(p) apply(p, 2L, stats::p.adjust, method = "BH")
    pauc <- function(p) {
        s <- data.frame(key = seq_along(truth), expected = as.vector(truth),
                        estimate = as.vector(fit$est), p = as.vector(p), s_own = as.vector(p),
                        truth_status = "scored", key_status = "ok", stringsAsFactors = FALSE)
        valPartialAuc(s, fpr_max = 0.01)$value
    }
    if (n_feat == 4840L) return(c(rank.pauc_diff = pauc(fit$p) - pauc(fit$p_ord)))
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

#' C4 and C4e: jackknife forms (ii) and (iii) on simulated feature matrices
#' @noRd
valCovMatrix <- function(cell, nsim, seed, form, n_feat, target_sims) {
    g <- cell$n
    lay <- valCovLayout(form, g)
    full_stats <- function(dat) valCovStats(valCovFit(dat$x, lay, form), dat$truth, n_feat, form)
    t0 <- proc.time()[["elapsed"]]
    target <- valWithSeed(TARGET_SEED + seed, {
        res <- replicate(target_sims, full_stats(valCovMatrixData(lay, n_feat, form)),
                         simplify = FALSE)
        Reduce(`+`, res) / target_sims
    })
    t1 <- proc.time()[["elapsed"]]
    n_stat <- length(target)
    hits <- array(FALSE, c(nsim, length(COVERAGE_LEVELS), n_stat))
    undefined <- matrix(FALSE, nsim, n_stat)
    valWithSeed(seed, for (s in seq_len(nsim)) {
        dat <- valCovMatrixData(lay, n_feat, form)
        full <- valCovStats(valCovFit(dat$x, lay, form), dat$truth, n_feat, form)
        dels <- vapply(seq_len(g), function(u) {
            keep <- lay$unit != u
            valCovStats(valCovFit(dat$x[, keep], lay[keep, ], form), dat$truth, n_feat, form)
        }, numeric(n_stat))
        dels <- matrix(dels, nrow = n_stat)
        for (k in seq_len(n_stat)) for (j in seq_along(COVERAGE_LEVELS)) {
            iv <- valJackknife(full = full[[k]], deletions = dels[k, ], conf = COVERAGE_LEVELS[[j]],
                               form = form)
            undefined[s, k] <- iv$status != "ok"
            hits[s, j, k] <- isTRUE(iv$ci_low <= target[[k]] && target[[k]] <= iv$ci_high)
        }
    })
    t2 <- proc.time()[["elapsed"]]
    out <- do.call(rbind, lapply(seq_len(n_stat), function(k) {
        nm <- strsplit(names(target)[[k]], ".", fixed = TRUE)[[1L]]
        valCoverageRows(cell, nm[[2L]], nm[[1L]], hits[, , k, drop = TRUE], "approx_run",
                        undefined = undefined[, k], target = target[[k]])
    }))
    valTimed(out, t1 - t0, t2 - t1, target_sims)
}

#' C7: the paired within-run bootstrap over 200 keys in two classes, iid within class
#' @noRd
valCovBootstrap <- function(cell, nsim, seed, cls, target_sims) {
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
    t0 <- proc.time()[["elapsed"]]
    target <- valWithSeed(TARGET_SEED + seed, mean(replicate(target_sims, statistic(simData()))))
    t1 <- proc.time()[["elapsed"]]
    b_max <- valBootstrapMinB(max(COVERAGE_LEVELS))
    hits <- matrix(FALSE, nsim, length(COVERAGE_LEVELS))
    valWithSeed(seed, for (s in seq_len(nsim)) {
        d <- simData()
        bs <- valBootstrap(d, statistic, B = b_max, conf = max(COVERAGE_LEVELS),
                           seed = seed * 10L + s, strata = "class", keep_replicates = TRUE)
        reps <- bs$replicates[is.finite(bs$replicates)]
        for (j in seq_along(COVERAGE_LEVELS)) {
            tail <- (1 - COVERAGE_LEVELS[[j]]) / 2
            q <- stats::quantile(reps, c(tail, 1 - tail), type = 6, names = FALSE)
            hits[s, j] <- q[[1L]] <= target && target <= q[[2L]]
        }
    })
    out <- valCoverageRows(cell, paste0("paired_", cls), cls, hits, "approx_within",
                           target = target)
    valTimed(out, t1 - t0, proc.time()[["elapsed"]] - t1, target_sims)
}

#' The regime table from coverage rows (section 10.2)
#'
#' One row per method x class x form x design x unit count x level. Over the data conditions and
#' statistics of a row: any `fail` fails it, except that a paired t count whose normal case passes
#' and skewed case fails is `flagged` (C2); otherwise any `conservative` makes it conservative;
#' `not_offered` where the method gives no interval at that count and level.
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
        data.frame(
            regime = r$regime[[1L]], method = r$method[[1L]], class = r$class[[1L]],
            form = r$form[[1L]], design = r$design[[1L]], transform = r$transform[[1L]],
            n = r$n[[1L]], level = r$level[[1L]], status = status, coverage = r$coverage[[w]],
            mc_se = r$mc_se[[w]], source = "harness",
            statistic = paste(unique(paste(r$statistic, r$condition, sep = "@")), collapse = ", "),
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    out[order(out$regime, out$method, out$class, out$n, out$level), , drop = FALSE]
}

#' Run coverage cells in one process (one chunk) and save the rows
#'
#' @param cells Cell ids.
#' @param nsim Simulations per cell.
#' @param out_file RDS file for the rows.
#' @param target_sims Target simulations.
#' @return The rows, invisibly.
#' @export
valCoverageChunk <- function(cells, nsim, out_file, target_sims = TARGET_SIMS) {
    all_cells <- valCoverageCells()
    rows <- lapply(cells, function(id) {
        valCoverageRegime(all_cells[all_cells$cell == id, ], nsim = nsim, target_sims = target_sims)
    })
    out <- do.call(rbind, rows)
    saveRDS(out, out_file)
    invisible(out)
}

#' Start chunk processes, each an `Rscript` with the thread variables at 1, and wait for all
#' @noRd
valRunChunksConcurrently <- function(groups, nsim, work_dir, target_sims) {
    rscript <- file.path(R.home("bin"), "Rscript")
    env <- c("current", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    procs <- lapply(seq_along(groups), function(k) {
        out_file <- file.path(work_dir, sprintf("chunk-%02d.rds", k))
        code <- sprintf(
            "msrval::valCoverageChunk(c(%s), nsim = %dL, out_file = '%s', target_sims = %dL)",
            paste0("'", groups[[k]], "'", collapse = ", "), as.integer(nsim), out_file,
            as.integer(target_sims)
        )
        processx::process$new(rscript, c("--vanilla", "-e", code), env = env,
                              stdout = file.path(work_dir, sprintf("chunk-%02d.log", k)),
                              stderr = "2>&1")
    })
    for (p in procs) p$wait()
    status <- vapply(procs, function(p) p$get_exit_status(), integer(1))
    files <- file.path(work_dir, sprintf("chunk-%02d.rds", seq_along(groups)))
    list(status = status, files = files)
}

#' Run T-COV in concurrent chunks and write `coverage-detail.tsv` and `coverage-table.tsv`
#'
#' Cells are dealt to the chunks in order of their expected cost; every cell has its own seed, so
#' the result does not depend on the chunking.
#'
#' @param ids Regime ids (`C1`, ..., `C8`) or cell ids; `NULL` for all.
#' @param nsim Simulations per cell (2,000 in T2).
#' @param chunks Concurrent processes.
#' @param out_dir Output folder.
#' @param target_sims Target simulations.
#' @return The regime table, invisibly.
#' @export
valRunCoverage <- function(ids = NULL, nsim = 2000L, chunks = 4L, out_dir,
                           target_sims = TARGET_SIMS) {
    cells <- valCoverageCells()
    if (!is.null(ids)) cells <- cells[cells$regime %in% ids | cells$cell %in% ids, ]
    if (nrow(cells) == 0L) valAbort("no coverage cell matches {ids}", "val_spec_invalid")
    work_dir <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "val-coverage")
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
    groups <- split(cells$cell, rep_len(seq_len(chunks), nrow(cells)))
    started <- Sys.time()
    run <- valRunChunksConcurrently(groups, nsim, work_dir, target_sims)
    wall <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    ok <- file.exists(run$files)
    detail <- do.call(rbind, lapply(run$files[ok], readRDS))
    missing <- setdiff(cells$cell, detail$cell)
    utils::write.table(detail, file.path(out_dir, "coverage-detail.tsv"), sep = "\t",
                       quote = FALSE, row.names = FALSE)
    table <- valSummariseCoverage(detail)
    utils::write.table(table, file.path(out_dir, "coverage-table.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE)
    writeLines(c(sprintf("wall_seconds\t%.1f", wall), sprintf("chunks\t%d", chunks),
                 sprintf("chunk_exit\t%s", paste(run$status, collapse = ",")),
                 sprintf("cells_missing\t%s", paste(missing, collapse = ","))),
               file.path(out_dir, "coverage-run.txt"))
    invisible(table)
}

#' Time one regime in one process and in four concurrent chunks (f_sim, G14 (c), section 10.4),
#' and probe the cost of one cell of every family
#'
#' @param id Regime id whose cells are timed.
#' @param nsim Simulations per timed cell.
#' @param out_dir Output folder: `fsim.tsv`, `cost-probe.tsv`.
#' @param chunks Concurrent chunks (4).
#' @param probe_nsim,probe_target_sims Sizes of the cost probe.
#' @return f_sim, invisibly.
#' @export
valTimeRegime <- function(id, nsim, out_dir, chunks = 4L, probe_nsim = 20L,
                          probe_target_sims = 200L) {
    cells <- valCoverageCells()
    ids <- cells$cell[cells$regime == id]
    work <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "val-fsim")
    dir.create(file.path(work, "one"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(work, "four"), recursive = TRUE, showWarnings = FALSE)
    wallOf <- function(groups, dir) {
        t0 <- Sys.time()
        run <- valRunChunksConcurrently(groups, nsim, dir, TARGET_SIMS)
        list(wall = as.numeric(difftime(Sys.time(), t0, units = "secs")), status = run$status)
    }
    one <- wallOf(list(ids), file.path(work, "one"))
    four <- wallOf(split(ids, rep_len(seq_len(chunks), length(ids))), file.path(work, "four"))
    f_sim <- one$wall / four$wall
    utils::write.table(data.frame(regime = id, cells = length(ids), nsim = nsim,
                                  wall_one_process = round(one$wall, 1),
                                  wall_four_chunks = round(four$wall, 1), f_sim = round(f_sim, 3),
                                  exit_one = paste(one$status, collapse = ","),
                                  exit_four = paste(four$status, collapse = ","),
                                  cores = valCoreCount()),
                       file.path(out_dir, "fsim.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE)
    fams <- cells[!duplicated(ifelse(cells$regime == "C3", "C3", cells$family)), ]
    probe <- do.call(rbind, lapply(seq_len(nrow(fams)), function(i) {
        r <- valCoverageRegime(fams[i, ], nsim = probe_nsim, target_sims = probe_target_sims)
        data.frame(cell = fams$cell[[i]], family = fams$family[[i]], n = fams$n[[i]],
                   probe_nsim = probe_nsim, probe_target_sims = probe_target_sims,
                   seconds_target = r$seconds_target[[1L]], seconds_sims = r$seconds_sims[[1L]],
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
