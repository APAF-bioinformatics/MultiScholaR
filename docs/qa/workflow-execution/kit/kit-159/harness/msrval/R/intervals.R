# Interval methods, the classifier and the refusal fill (KIT-159 sections 4.4 to 4.7).

LABEL_WITHIN_RUN <- paste(
    "within run: conditional on the samples; not a confidence interval for a repeat of the",
    "experiment"
)
LABEL_MONTE_CARLO <- "Monte Carlo interval over draws from one pool; not between-study uncertainty"

#' Construct an interval record
#'
#' Every interval of the harness is reported in this one form (section 4.4).
#'
#' @param ci_low,ci_high Bounds (NA when there is no interval).
#' @param method Method name, or `ticket:<function>`.
#' @param n_units Units (or keys, or deletions) the interval rests on.
#' @param conf Nominal level.
#' @param scope `run_level`, `monte_carlo` or `within_run`.
#' @param class Statistic class of the metric (`accuracy`, `threshold_count`, `rank`) or NA.
#' @param form Unit form (`unit`, `replicate`, `block`, `key_block`) or NA.
#' @param status `ok` or the reason there is no usable interval.
#' @param estimate Point estimate on the full data.
#' @param ... Further fields (`se`, `critical`, `achieved_level`, `flags`, ...).
#' @return A list of class `msrval_interval`.
#' @export
valNewInterval <- function(ci_low, ci_high, method, n_units, conf = 0.95, scope = "run_level",
                           class = NA_character_, form = NA_character_, status = "ok",
                           estimate = NA_real_, ...) {
    scope <- match.arg(scope, c("run_level", "monte_carlo", "within_run"))
    label <- switch(scope, within_run = LABEL_WITHIN_RUN, monte_carlo = LABEL_MONTE_CARLO,
                    run_level = "")
    x <- list(method = method, scope = scope, estimate = estimate, ci_low = ci_low,
              ci_high = ci_high, conf = conf, achieved_level = NA_real_,
              n_units = as.integer(n_units), status = status, flags = character(),
              class = class, form = form, se = NA_real_, critical = NA_real_, label = label)
    extra <- list(...)
    x[names(extra)] <- extra
    structure(x, class = "msrval_interval")
}

#' Scope from the unit role, refusing dependent units
#' @noRd
valScopeOf <- function(unit_role, default) {
    switch(unit_role,
           dependent = valAbort("dependent units never enter an interval", "val_spec_invalid"),
           draw = "monte_carlo",
           default)
}

#' An interval by a method of section 4.4
#'
#' Per-unit values `x`; with `y` the paired difference x - y; with `fixed` the one-sample form
#' x - fixed. For `clopper_pearson` and `wilson`, `x` is the count and `n` the number of trials; for
#' `newcombe_paired`, `x` is the 2 x 2 table (both, first only, second only, neither).
#'
#' @param x Values, a count, or the paired table.
#' @param method One of `paired_t`, `exact_median`, `clopper_pearson`, `wilson`,
#'   `newcombe_paired`.
#' @param conf Level.
#' @param y Comparator values (paired form).
#' @param fixed Fixed value (one-sample form).
#' @param n Trials (proportions).
#' @param unit_role `independent`, `draw` (Monte Carlo scope) or `dependent` (refused).
#' @param class Statistic class.
#' @param scope Scope override for `clopper_pearson` over keys (`within_run`).
#' @return An interval record ([valNewInterval()]).
#' @export
valInterval <- function(x, method = c("paired_t", "exact_median", "clopper_pearson", "wilson",
                                      "newcombe_paired"),
                        conf = 0.95, y = NULL, fixed = NULL, n = NULL,
                        unit_role = c("independent", "draw", "dependent"), class = NA_character_,
                        scope = NULL) {
    method <- match.arg(method)
    unit_role <- match.arg(unit_role)
    checkmate::assert_number(conf, lower = 0.5, upper = 1)
    if (method %in% c("paired_t", "exact_median")) {
        sc <- valScopeOf(unit_role, "run_level")
        d <- if (is.null(y)) x else x - y
        if (!is.null(fixed)) d <- d - fixed
        d <- d[is.finite(d)]
        f <- if (method == "paired_t") valPairedT else valExactMedian
        return(f(d, conf, sc, class))
    }
    if (method == "newcombe_paired") return(valNewcombePaired(x, conf, class))
    checkmate::assert_count(x)
    checkmate::assert_count(n, positive = TRUE)
    if (x > n) valAbort("count {x} exceeds trials {n}", "val_spec_invalid")
    if (method == "wilson") {
        b <- valWilsonBounds(x, n, conf)
        return(valNewInterval(b[[1L]], b[[2L]], "wilson", n, conf, "within_run", class,
                              estimate = x / n))
    }
    sc <- scope %||% valScopeOf(unit_role, "run_level")
    a <- 1 - conf
    lo <- if (x == 0) 0 else stats::qbeta(a / 2, x, n - x + 1)
    hi <- if (x == n) 1 else stats::qbeta(1 - a / 2, x + 1, n - x)
    valNewInterval(lo, hi, "clopper_pearson", n, conf, sc, class, estimate = x / n,
                   flags = "conservative")
}

#' Paired t interval on per-unit differences
#' @noRd
valPairedT <- function(d, conf, scope, class) {
    n <- length(d)
    if (n < 2L) {
        return(valNewInterval(NA_real_, NA_real_, "paired_t", n, conf, scope, class,
                              status = "too_few_units"))
    }
    se <- stats::sd(d) / sqrt(n)
    crit <- valTCritical(conf, n - 1L)
    m <- mean(d)
    valNewInterval(m - crit * se, m + crit * se, "paired_t", n, conf, scope, class, estimate = m,
                   se = se, critical = crit)
}

#' Order-statistic (sign-test) interval for the median of per-unit values
#' @noRd
valExactMedian <- function(d, conf, scope, class) {
    n <- length(d)
    if (n == 0L || 1 - 2^(1 - n) < conf - 1e-12) {
        return(valNewInterval(NA_real_, NA_real_, "exact_median", n, conf, scope, class,
                              status = "too_few_units"))
    }
    ks <- seq_len(floor(n / 2))
    achieved <- 1 - 2 * stats::pbinom(ks - 1L, n, 0.5)
    k <- max(ks[achieved >= conf - 1e-12])
    xs <- sort(d)
    valNewInterval(xs[k], xs[n - k + 1L], "exact_median", n, conf, scope, class,
                   estimate = stats::median(d), achieved_level = achieved[k])
}

#' Wilson score bounds
#' @noRd
valWilsonBounds <- function(x, n, conf) {
    z <- stats::qnorm(1 - (1 - conf) / 2)
    p <- x / n
    den <- 1 + z^2 / n
    centre <- (p + z^2 / (2 * n)) / den
    half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
    c(max(0, centre - half), min(1, centre + half))
}

#' Newcombe's paired difference of proportions, method 10
#' @noRd
valNewcombePaired <- function(tab, conf, class) {
    checkmate::assert_numeric(tab, len = 4L, lower = 0, any.missing = FALSE)
    if (!is.null(names(tab))) {
        tab <- tab[c("both", "first_only", "second_only", "neither")]
    }
    a <- tab[[1L]]
    b <- tab[[2L]]
    cc <- tab[[3L]]
    d <- tab[[4L]]
    n <- a + b + cc + d
    p1 <- (a + b) / n
    p2 <- (a + cc) / n
    w1 <- valWilsonBounds(a + b, n, conf)
    w2 <- valWilsonBounds(a + cc, n, conf)
    num <- a * d - b * cc
    if (num > 0) num <- max(num - n / 2, 0)
    den <- sqrt((a + b) * (cc + d) * (a + cc) * (b + d))
    phi <- if (den == 0) 0 else num / den
    diff <- p1 - p2
    delta <- sqrt((p1 - w1[1L])^2 - 2 * phi * (p1 - w1[1L]) * (w2[2L] - p2) + (w2[2L] - p2)^2)
    eps <- sqrt((w1[2L] - p1)^2 - 2 * phi * (w1[2L] - p1) * (p2 - w2[1L]) + (p2 - w2[1L])^2)
    valNewInterval(diff - delta, diff + eps, "newcombe_paired", n, conf, "within_run", class,
                   estimate = diff, phi = phi)
}

#' Is a jackknife statistic a median or quantile (refused: the delete-one jackknife is
#' inconsistent for non-smooth statistics)?
#' @noRd
valIsNonSmooth <- function(statistic) {
    if (is.character(statistic)) return(statistic %in% c("median", "quantile"))
    identical(statistic, stats::median) || identical(statistic, stats::quantile) ||
        identical(statistic, stats::median.default)
}

#' The delete-one jackknife (section 4.4)
#'
#' Either per-unit `values` with a smooth `statistic` (delete-one-unit), or a `full` value and the
#' `deletions` of a separately run deletion set (forms `replicate`, `block`). SE^2 = (g - 1) / g *
#' sum((theta_b - mean)^2); interval theta(full) +/- c SE with c = t(g - 1) at the level, or a
#' declared larger `critical`. `transform = "log"` forms the interval for log theta and
#' back-transforms it. A non-finite (or, under log, non-positive) value gives
#' `status = "undefined_deletion"`. The bias estimate (g - 1)(mean - full) is reported; with a
#' `margin` it is flagged when its absolute value exceeds half the margin.
#'
#' @param full Statistic on the full unit set.
#' @param deletions Statistic on each deletion.
#' @param values Per-unit values (with `statistic`).
#' @param statistic A smooth function of the per-unit values (never a median or quantile).
#' @param conf Level.
#' @param form `unit`, `replicate`, `block` or `key_block` (within run).
#' @param critical Declared critical value (at least t(g - 1)).
#' @param transform `NULL` or `"log"`.
#' @param margin The metric's margin, for the bias flag.
#' @param class Statistic class.
#' @param unequal_deletion A deletion left a group whole (form `replicate`).
#' @param design For form `unit` over designs: the design family whose T-COV regime licenses it.
#' @return An interval record with `se`, `critical`, `bias`.
#' @export
valJackknife <- function(full = NULL, deletions = NULL, values = NULL, statistic = NULL,
                         conf = 0.95, form = c("unit", "replicate", "block", "key_block"),
                         critical = NULL, transform = NULL, margin = NULL,
                         class = NA_character_, unequal_deletion = FALSE, design = NA_character_) {
    form <- match.arg(form)
    if (!is.null(statistic) && valIsNonSmooth(statistic)) {
        valAbort("the jackknife is refused for a median or quantile; use exact_median",
                 "val_spec_invalid")
    }
    if (!is.null(values)) {
        if (!is.function(statistic)) {
            valAbort("{.arg statistic} must be a function of the per-unit values",
                     "val_spec_invalid")
        }
        full <- statistic(values)
        deletions <- vapply(seq_along(values), function(i) statistic(values[-i]), numeric(1))
    }
    checkmate::assert_number(full, na.ok = TRUE)
    checkmate::assert_numeric(deletions)
    scope <- if (form == "key_block") "within_run" else "run_level"
    g <- length(deletions)
    flags <- if (isTRUE(unequal_deletion)) "unequal_deletion" else character()
    blank <- function(status) {
        valNewInterval(NA_real_, NA_real_, "jackknife", g, conf, scope, class, form,
                       status = status, estimate = full, flags = flags, design = design)
    }
    if (g < 4L) return(blank("too_few_units"))
    vals <- c(full, deletions)
    log_bad <- !is.null(transform) && identical(transform, "log") && any(vals <= 0, na.rm = TRUE)
    if (any(!is.finite(vals)) || log_bad) return(blank("undefined_deletion"))
    f <- if (identical(transform, "log")) log else identity
    finv <- if (identical(transform, "log")) exp else identity
    th <- f(deletions)
    se <- sqrt((g - 1) / g * sum((th - mean(th))^2))
    t_crit <- valTCritical(conf, g - 1L)
    crit <- t_crit
    if (!is.null(critical)) {
        if (critical < t_crit - 1e-12) {
            valAbort("a declared critical value ({critical}) below t(g - 1) = {round(t_crit, 4)}",
                     "val_spec_invalid")
        }
        crit <- critical
        flags <- c(flags, "calibrated_critical")
    }
    bias <- (g - 1) * (mean(deletions) - full)
    if (!is.null(margin) && abs(bias) > margin / 2) flags <- c(flags, "depends_on_replicate_count")
    centre <- f(full)
    valNewInterval(finv(centre - crit * se), finv(centre + crit * se), "jackknife", g, conf, scope,
                   class, form, estimate = full, se = se, critical = crit, bias = bias,
                   flags = flags, transform = if (is.null(transform)) "identity" else transform,
                   design = design)
}

#' Smallest bootstrap B at a level: B >= 50 / tail probability (A07)
#'
#' @param conf Level.
#' @return The smallest admissible B.
#' @export
valBootstrapMinB <- function(conf) {
    ceiling(round(50 / ((1 - conf) / 2), 6L))
}

#' Paired percentile bootstrap over keys, within one run
#'
#' One set of resampled rows per replicate is applied to the whole (wide) table, so every arm of
#' the comparison is resampled together. Strata keep their sizes; with `cluster`, clusters are
#' resampled whole (contract row 20). More than 1 % non-finite replicates: `status = "unstable"`.
#'
#' @param data A data frame, one row per key, the arms' columns side by side.
#' @param statistic Function of a resampled data frame returning one number.
#' @param B Replicates (at least [valBootstrapMinB()]).
#' @param conf Level.
#' @param seed Seed from the declaration.
#' @param strata Column naming the stratum of each row, or `NULL`.
#' @param cluster Column naming the cluster of each row, or `NULL`.
#' @param events Column of a logical outcome for a bootstrapped proportion, or `NULL`.
#' @param min_events Fewest keys in the rarer outcome.
#' @param class Statistic class.
#' @param keep_replicates Return the replicates (`replicates`), for the coverage regimes.
#' @return An interval record (scope `within_run`).
#' @export
valBootstrap <- function(data, statistic, B, conf = 0.95, seed, strata = NULL, cluster = NULL,  # nolint
                         events = NULL, min_events = 5L, class = NA_character_,
                         keep_replicates = FALSE) {
    checkmate::assert_data_frame(data, min.rows = 1L)
    checkmate::assert_function(statistic)
    checkmate::assert_int(seed)
    if (B < valBootstrapMinB(conf)) {
        valAbort("B = {B} is below 50 / tail probability ({valBootstrapMinB(conf)})",
                 "val_spec_invalid")
    }
    estimate <- statistic(data)
    if (!is.null(events)) {
        ev <- as.logical(data[[events]])
        if (min(sum(ev), sum(!ev)) < min_events) {
            return(valNewInterval(NA_real_, NA_real_, "bootstrap_within_run", nrow(data), conf,
                                  "within_run", class, status = "too_few_events",
                                  estimate = estimate))
        }
    }
    n <- nrow(data)
    unit_id <- if (is.null(cluster)) seq_len(n) else data[[cluster]]
    units <- unique(unit_id)
    rows_by_unit <- split(seq_len(n), factor(unit_id, levels = units))
    first_row <- vapply(rows_by_unit, `[`, integer(1), 1L)
    stratum <- if (is.null(strata)) rep(1L, length(units)) else data[[strata]][first_row]
    by_stratum <- split(seq_along(units), stratum)
    reps <- valWithSeed(seed, vapply(seq_len(B), function(b) {
        pick <- unlist(lapply(by_stratum, function(u) {
            u[sample.int(length(u), length(u), replace = TRUE)]
        }), use.names = FALSE)
        statistic(data[unlist(rows_by_unit[pick], use.names = FALSE), , drop = FALSE])
    }, numeric(1)))
    na_share <- mean(!is.finite(reps))
    if (na_share > 0.01) {
        return(valNewInterval(NA_real_, NA_real_, "bootstrap_within_run", n, conf, "within_run",
                              class, status = "unstable", estimate = estimate, B = B,
                              seed = seed, na_share = na_share))
    }
    tail <- (1 - conf) / 2
    q <- stats::quantile(reps[is.finite(reps)], c(tail, 1 - tail), type = 6, names = FALSE)
    valNewInterval(q[[1L]], q[[2L]], "bootstrap_within_run", n, conf, "within_run", class,
                   estimate = estimate, B = B, seed = seed, na_share = na_share,
                   replicates = if (keep_replicates) reps)
}

#' The harness's committed coverage table (T2's `coverage-table.tsv`), or an empty table
#' @noRd
valHarnessRegimes <- function() {
    path <- system.file("extdata", "coverage-table.tsv", package = "msrval")
    if (!nzchar(path)) {
        return(data.frame(regime = character(), method = character(), class = character(),
                          form = character(), n = integer(), level = numeric(),
                          status = character(), coverage = numeric(), mc_se = numeric(),
                          source = character(), statistic = character()))
    }
    reg <- utils::read.delim(path, stringsAsFactors = FALSE)
    # O5 (VAL-109): T2's run-level rows were classified under the band as run; r1.13 reads
    # `conservative` on the miss rate (below a quarter of the nominal). Both labels are kept: the
    # second is r1.13's rule applied to the row's best condition (the rank rows of T2R were
    # classified under r1.13 already, so both labels agree there).
    approx_run <- reg$regime %in% c("C2", "C3", "C4", "C4e")
    second <- approx_run & reg$status == "pass" &
        1 - suppressWarnings(as.numeric(reg$coverage_max)) < (1 - reg$level) / 4 - 1e-12
    reg$status_r113 <- ifelse(second %in% TRUE, "conservative", reg$status)
    reg
}

#' The regime that licenses an interval (section 4.4), or a refusal
#' @noRd
valFindRegime <- function(interval, regimes, accept_flagged) {
    regimes <- regimes %||% valHarnessRegimes()
    is_ticket <- startsWith(interval$method, "ticket:")
    hit <- regimes$method == interval$method &
        regimes$n == interval$n_units &
        abs(regimes$level - interval$conf) < 5e-5 &
        (is.na(regimes$class) | regimes$class %in% interval$class) &
        (is.na(regimes$form) | regimes$form %in% interval$form)
    if (!is.null(regimes$design)) {
        hit <- hit & (is.na(regimes$design) | regimes$design %in% interval[["design"]])
    }
    if (!is.null(regimes$transform)) {
        tr <- interval[["transform"]] %||% "identity"
        hit <- hit & (is.na(regimes$transform) | regimes$transform %in% tr)
    }
    if (is_ticket) hit <- hit & regimes$source == "ticket"
    m <- regimes[hit %in% TRUE, , drop = FALSE]
    good <- m[m$status %in% c("pass", "conservative"), , drop = FALSE]
    # O4 (VAL-109): a `conservative` row with no observed miss licenses only the statistic it
    # simulated (the declared `interval$statistic` must name it)
    if (nrow(good) > 0L) {
        cmin <- suppressWarnings(as.numeric(good[["coverage_min"]] %||% good$coverage))
        zero <- good$status == "conservative" & cmin >= 1 - 1e-12
        own <- vapply(strsplit(as.character(good$statistic), ", ", fixed = TRUE), function(s) {
            stat <- interval[["statistic"]]
            !is.null(stat) && stat %in% sub("@.*$", "", s)
        }, logical(1))
        refused_zero <- good[zero %in% TRUE & !own, , drop = FALSE]
        good <- good[!(zero %in% TRUE) | own, , drop = FALSE]
        if (nrow(good) == 0L) {
            why <- paste0("regime ", refused_zero$regime[[1L]], " has no observed miss: it ",
                          "licenses only the statistic it simulated (",
                          refused_zero$statistic[[1L]], "); declare interval$statistic (O4, ",
                          "VAL-109)")
            valAbort(why, "val_spec_invalid")
        }
    }
    if (nrow(good) > 0L) return(good[1L, ])
    flagged <- m[m$status %in% "flagged", , drop = FALSE]
    if (nrow(flagged) > 0L) {
        if (is.null(accept_flagged[["label"]]) || is.null(accept_flagged[["ledger_row"]])) {
            valAbort("regime {flagged$regime[1]} is flagged; accept_flagged is not declared",
                     "val_spec_invalid")
        }
        return(flagged[1L, ])
    }
    valAbort(paste("no passing coverage regime for {interval$method}, class {interval$class},",
                   "n {interval$n_units}, level {interval$conf}"), "val_spec_invalid")
}

#' Classify an interval against a margin (section 4.5)
#'
#' Oriented so that positive is favourable: `better` lower > margin; `match` lower >= -margin and
#' upper <= margin; `worse` upper < -margin; otherwise `inconclusive`; `no_worse` lower > -margin.
#'
#' @param interval An interval record of the difference arm - comparator, or of one arm with
#'   `fixed`.
#' @param margin The declared margin (never defaulted).
#' @param higher_is_better Orientation of the metric.
#' @param fixed A fixed value for a one-arm comparison, or `NULL`.
#' @param view The view the metric was computed in (`own` only for a one-arm comparison).
#' @param feature_interval_decides Declared on this comparison: a within-run interval may decide.
#' @param regimes Regime table (default: the harness's committed coverage table).
#' @param accept_flagged `list(label, ledger_row)` admitting a flagged regime.
#' @return A list: `outcome`, `no_worse`, `status`, `regime`, `coverage`, `mc_se`, `label`.
#' @export
valClassify <- function(interval, margin, higher_is_better, fixed = NULL, view = "intersection",
                        feature_interval_decides = FALSE, regimes = NULL,
                        accept_flagged = NULL) {
    if (is.null(margin) || !checkmate::test_number(margin, lower = 0)) {
        valAbort("a margin must be declared (never defaulted)", "val_spec_invalid")
    }
    checkmate::assert_flag(higher_is_better)
    if (identical(view, "own") && is.null(fixed)) {
        valAbort("a comparison between arms in the own view", "val_spec_invalid")
    }
    if (identical(interval$scope, "within_run") && !isTRUE(feature_interval_decides)) {
        valAbort("a within-run interval decides only with feature_interval_decides",
                 "val_spec_invalid")
    }
    lo <- interval$ci_low
    hi <- interval$ci_high
    if (!identical(interval$status, "ok")) {
        if (is.finite(lo) && is.finite(hi) && lo == hi) {
            valAbort("a zero-width interval from a method with status {interval$status}",
                     "val_spec_invalid")
        }
        return(list(outcome = NA_character_, no_worse = NA, status = interval$status,
                    regime = NA_character_, coverage = NA_real_, mc_se = NA_real_,
                    label = interval$label))
    }
    reg <- valFindRegime(interval, regimes, accept_flagged)
    if (!is.null(fixed)) {
        lo <- lo - fixed
        hi <- hi - fixed
    }
    if (!higher_is_better) {
        tmp <- -hi
        hi <- -lo
        lo <- tmp
    }
    outcome <- if (lo > margin) {
        "better"
    } else if (lo >= -margin && hi <= margin) {
        "match"
    } else if (hi < -margin) {
        "worse"
    } else {
        "inconclusive"
    }
    labels <- c(interval[["label"]], if (reg$status == "flagged") accept_flagged[["label"]])
    # O1 (VAL-109): the bands stay as run; a regime whose worst-condition miss rate is at least
    # twice nominal is flagged beside the verdict.
    worst <- suppressWarnings(as.numeric(reg[["coverage_min"]] %||% reg$coverage))
    list(outcome = outcome, no_worse = lo > -margin, status = "ok", regime = reg$regime,
         regime_status = reg$status,
         regime_status_r113 = as.character(reg[["status_r113"]] %||% reg$status),
         coverage = reg$coverage, mc_se = reg$mc_se,
         coverage_min = worst, nsim = suppressWarnings(as.numeric(reg[["nsim"]] %||% NA)),
         statistic = as.character(reg[["statistic"]] %||% NA_character_),
         miss_twice_nominal = isTRUE(1 - worst >= 2 * (1 - interval$conf) - 1e-12),
         label = paste(labels[nzchar(labels)], collapse = "; "))
}

#' Classify a compound member (section 4.6)
#'
#' The member is `worse` if any non-guard component is (or, with `worse_when`, if every component
#' of some group is); `better` or `match` only if every component is; otherwise `inconclusive`;
#' `no_worse` only if every component is. A component whose control is not `worse` is
#' `not_sensitive` and leaves the conjunction; no component left: `no_verdict`. A guard never sets
#' the outcome: a guard that does not pass (or is not sensitive) makes the member at best
#' `inconclusive`. Guards are printed as `guard <component>: <outcome>`, `worse` as `contradicted`.
#'
#' @param member `list(all_of, guards, controls, worse_when)`.
#' @param outcomes Named list of [valClassify()] results, components and controls.
#' @return A list: `outcome`, `no_worse`, `components`, `guard_text`.
#' @export
valClassifyMember <- function(member, outcomes) {
    comps <- member[["all_of"]]
    guards <- vapply(member[["guards"]] %||% list(), function(g) g[["component"]], character(1))
    requires <- vapply(member[["guards"]] %||% list(), function(g) g[["require"]] %||% "no_worse",
                       character(1))
    controls <- member[["controls"]] %||% list()
    named <- c(comps, guards, unlist(controls), unlist(member[["worse_when"]]))
    if (any(unlist(member[["worse_when"]]) %in% guards)) {
        valAbort("a worse_when group names a guard (a guard never sets the outcome)",
                 "val_spec_invalid")
    }
    if (!all(named %in% names(outcomes)) || !all(guards %in% comps)) {
        valAbort("a member names a component that is not a declared comparison",
                 "val_spec_invalid")
    }
    oc <- function(c) outcomes[[c]][["outcome"]] %||% NA_character_
    not_sensitive <- vapply(comps, function(c) {
        ctl <- controls[[c]]
        !is.null(ctl) && !identical(oc(ctl), "worse")
    }, logical(1))
    role <- ifelse(comps %in% guards, "guard", "component")
    active <- comps[role == "component" & !not_sensitive]
    components <- data.frame(component = comps, role = role,
                             outcome = vapply(comps, oc, character(1)),
                             not_sensitive = not_sensitive, stringsAsFactors = FALSE,
                             row.names = NULL)
    guard_out <- vapply(guards, oc, character(1))
    guard_text <- sprintf("guard %s: %s", guards,
                          ifelse(guard_out %in% "worse", "contradicted", guard_out))
    if (length(active) == 0L) {
        return(list(outcome = "no_verdict", no_worse = NA, components = components,
                    guard_text = guard_text))
    }
    a_out <- vapply(active, oc, character(1))
    missing <- active[!a_out %in% c("better", "match", "worse", "inconclusive")]
    if (length(missing)) {
        # I16: an active component without a verdict (refused, not evaluable) leaves the member
        # not evaluable; it is never printed as a statistical word.
        return(list(outcome = "not_evaluable", no_worse = NA, components = components,
                    guard_text = guard_text, missing = missing))
    }
    a_nw <- vapply(active, function(c) isTRUE(outcomes[[c]]$no_worse), logical(1))
    is_worse <- if (is.null(member[["worse_when"]])) {
        any(a_out %in% "worse")
    } else {
        any(vapply(member[["worse_when"]], function(grp) {
            grp <- intersect(unlist(grp), active)
            length(grp) > 0L && all(vapply(grp, function(c) identical(oc(c), "worse"),
                                           logical(1)))
        }, logical(1)))
    }
    outcome <- if (is_worse) {
        "worse"
    } else if (all(a_out %in% "better")) {
        "better"
    } else if (all(a_out %in% "match")) {
        "match"
    } else {
        "inconclusive"
    }
    no_worse <- all(a_nw) && !is_worse
    guard_pass <- vapply(seq_along(guards), function(i) {
        g <- guards[[i]]
        if (requires[[i]] == "match") identical(oc(g), "match") else isTRUE(outcomes[[g]]$no_worse)
    }, logical(1))
    guard_blocked <- any(!guard_pass) || any(not_sensitive & role == "guard")
    if (guard_blocked) {
        if (outcome %in% c("better", "match")) outcome <- "inconclusive"
        no_worse <- FALSE
    }
    list(outcome = outcome, no_worse = no_worse, components = components, guard_text = guard_text)
}

#' Refusals in a run-level comparison (section 4.7; O7 ruled in VAL-109)
#'
#' Fewer than 30 units: more than one refused unit is `not_evaluable`; exactly one is dropped,
#' with the `refusal_substitute` sensitivity arm required. 30 or more: within `tolerance` (at most
#' 5 % of the units), the verdict with the refused units dropped stands only if it is unchanged when
#' each refused unit takes the least favourable observed value; beyond it, `not_evaluable`. After a
#' refusal the interval rests on n - k units: a verdict is given only where a coverage regime for
#' that count is registered (the harness's or the consuming ticket's, option (b) of O7); otherwise
#' the comparison is `not_evaluable` with the reason in `regime`. Never an abort.
#'
#' @param d Per-unit differences (or one-arm values), NA for a refused unit.
#' @param tolerance The declared `refusal_tolerance`.
#' @param higher_is_better Orientation.
#' @param margin Margin.
#' @param method Run-level method (`paired_t`, `exact_median` or `clopper_pearson`).
#' @param conf Level.
#' @param regimes Regime table for the classifier.
#' @param accept_flagged As [valClassify()].
#' @param fixed The fixed value of a one-arm comparison, or `NULL`.
#' @param unit_role `independent` or `draw` (Monte Carlo scope).
#' @param class Statistic class.
#' @param view The comparison's view.
#' @return A list: `status`, `outcome`, `no_worse`, `outcome_dropped`, `outcome_filled`,
#'   `interval_dropped`, `interval_filled`, `n_refused`, `refusal_substitute_required`,
#'   `classification`, `flags`.
#' @export
valRefusalFill <- function(d, tolerance, higher_is_better, margin, method = "paired_t",
                           conf = 0.95, regimes = NULL, accept_flagged = NULL, fixed = NULL,
                           unit_role = "independent", class = NA_character_,
                           view = "intersection") {
    checkmate::assert_numeric(d, min.len = 1L)
    checkmate::assert_count(tolerance)
    n <- length(d)
    refused <- !is.finite(d)
    n_ref <- as.integer(sum(refused))
    if (n >= 30L && tolerance > floor(0.05 * n)) {
        valAbort("refusal_tolerance {tolerance} exceeds 5 % of {n} units", "val_spec_invalid")
    }
    counted <- if (n_ref > 0L) sprintf("n_refused=%d", n_ref) else character()
    result <- list(status = "not_evaluable", outcome = NA_character_, no_worse = NA,
                   outcome_dropped = NA_character_, outcome_filled = NA_character_,
                   interval_dropped = NULL, interval_filled = NULL, n_refused = n_ref,
                   refusal_substitute_required = FALSE,
                   classification = list(regime = paste("not_evaluable:", n_ref,
                                                        "refused units exceed the rule (4.7)")),
                   flags = c("not_evaluable", counted))
    limit <- if (n < 30L) 1L else tolerance
    if (n_ref > limit) return(result)
    classify <- function(v) {
        iv <- valLevelInterval(v, method, conf, unit_role, class)
        refusedRegime <- function(e) {
            why <- if (n_ref > 0L) "not_evaluable: " else "refused: "
            after <- ""
            if (n_ref > 0L) {
                after <- paste0(" (after ", n_ref, " refused units the interval rests on ",
                                length(v), " units; no regime is registered for that count;",
                                " O7 (b), VAL-109)")
            }
            list(outcome = NA_character_, no_worse = NA, refused = TRUE,
                 regime = paste0(why, conditionMessage(e), after))
        }
        cl <- tryCatch(valClassify(iv, margin, higher_is_better, fixed = fixed, view = view,
                                   regimes = regimes, accept_flagged = accept_flagged),
                       val_spec_invalid = refusedRegime)
        list(iv = iv, cl = cl)
    }
    kept <- d[!refused]
    dropped <- classify(kept)
    result$interval_dropped <- dropped$iv
    result$outcome_dropped <- dropped$cl$outcome
    result$classification <- dropped$cl
    result$status <- "ok"
    result$outcome <- dropped$cl$outcome
    result$no_worse <- dropped$cl$no_worse
    result$flags <- counted
    if (n_ref == 0L) return(result)
    if (isTRUE(dropped$cl$refused)) {
        result$status <- "not_evaluable"
        result$flags <- c("not_evaluable", counted)
        return(result)
    }
    if (n < 30L) {
        result$refusal_substitute_required <- TRUE
        result$flags <- c(counted, "refusal_substitute_required")
        return(result)
    }
    filled_d <- d
    filled_d[refused] <- if (higher_is_better) min(kept) else max(kept)
    filled <- classify(filled_d)
    result$interval_filled <- filled$iv
    result$outcome_filled <- filled$cl$outcome
    same <- identical(filled$cl$outcome, dropped$cl$outcome) &&
        identical(filled$cl$no_worse, dropped$cl$no_worse)
    if (!same) {
        result$status <- "not_evaluable"
        result$outcome <- NA_character_
        result$no_worse <- NA
        result$flags <- c("not_evaluable", counted)
    }
    result
}
