# KIT-210: scoring of the candidate remedies (contrasts, classifications, J, safe rates, the decision rules, H-rate,
# kappa*'s fit and validation, the decision table, the terminal outcomes) and the comparison modes of the lifecycle.
# Evidence code under docs/qa, not package code. A copy of `kit-185/score-k185.R` (SHA-256
# 5a708c6e3e34feecec60b36836eaa4d090b32ded35ef315a79eea4a0593f71eb) whose helpers `summarise`, `classify`, `clopper` are
# kept (with the limits passed in) and whose `pairedContrasts` gains the feature-set tag; everything else is new and
# follows KIT-210 revision 1.3, sections 6.3 to 6.10. Run once on a runner after the counted rows are committed.
#
# Modes (from the package root):
#   Rscript score-k210.R <records dir> <out dir>                 the scoring (rows/*.tsv.gz, checks/*.tsv)
#   Rscript score-k210.R --k185 <kit-185 dir> <check dir>         E1'': this file's helpers on KIT-185's committed rows
#                                                                 reproduce kit-185/results/cells.tsv for S and C
#   Rscript score-k210.R --compare-e1 <rows 185> <rows 210> <check dir>   E1: sim-k210.R (candidates off, seed base
#                                                                 185000000) equals sim-k185.R on the same units
#   Rscript score-k210.R --compare-e10 <rows 210> <kit-185 dir> <check dir>  E10: the bound commit's inputs to RUV-III-C
#                                                                 against KIT-185's committed rows, replicates 1-3
#   Rscript score-k210.R --scale <private dir> <check dir>        the scorer on synthetic rows of the counted run's size
# The comparison and scale modes write one small table to <check dir> (pass or fail, counts, seconds; no estimand
# value) and keep everything else in <private dir> or in memory (the lifecycle's blindness, KIT-210 6.9).

args <- commandArgs(TRUE)
mode <- if (length(args) && startsWith(args[[1L]], "--")) args[[1L]] else "score"

# ---- declared constants (KIT-185 5.5 to 5.7 unchanged; KIT-210 6.3 to 6.6) ------------------------------------------
DECIDING <- c("P0", "P0r", "P1", "P2", "P3", "P4", "P5", "L1", "L2", "M1", "M2", "P4b", "L1b", "M1b")
D5_CELLS <- c("P4b", "L1b", "M1b")
BASE_MASKS <- c("vd0.1", "vd0.5", "vd2", "vd10", "vd25", "vd50", "mcar0.5", "mcar2", "mcar10")
D5_MASKS <- c("vd2", "vd10", "bb2", "bb10", "bb25")
FAMILY <- 11L * 9L * 4L + 3L * 5L * 4L # 456
DEPARTURE_FAMILY <- 14L * 4L # 56 (KIT-210 6.5)
ALPHA <- 0.05
EDGES <- data.frame(estimand = c("B", "Cstar", "fdp", "power", "A1", "null_rate", "excess_lost"),
    material = c(0.05, 0.02, 0.010, 0.02, 0.05, 0.010, 0.02), serious = c(0.25, 0.10, 0.025, 0.10, 0.25, 0.025, 0.10),
    stringsAsFactors = FALSE)
# The route classes and their deciding estimands (KIT-210 6.4, D6): excess features lost decides on the declared class
# with the material edge 0.02 (the power band's); its serious edge is taken from the same band (0.10), declared here.
ROUTE_CLASSES <- list(legacy = c("B", "Cstar", "fdp", "power"), declared = c("B", "Cstar", "excess_lost"))
K_OF <- c(P0 = 2, P0r = 2, P1 = 2, P2 = 2, P3 = 2, P4 = 2, P5 = 2, L1 = 2, L2 = 2, M1 = 2, M2 = 2, P0lo = 2, P0hi = 2,
    P4b = 4, L1b = 2, M1b = 15)
CANDIDATES <- c("C3", "C4a", "C4b", "C6", "C7")
REPORTED_CANDIDATES <- "C7r5"
K_FAMILIES <- 1L + length(CANDIDATES) # S, C3, C4a, C4b, C6, C7 (KIT-210 6.5)
MIN_ADMITTED <- 30L
CODE_CHANGE <- c(S = "none (the shipped call)", C3 = "a new estimator (directions from complete features) and a guard",
    C4a = "a new estimator (directions per pattern of observed runs) and a guard",
    C4b = "a new estimator (one W from complete features, limma's removeBatchEffect) and a guard",
    C6 = "one line before the call (centre each feature, add the mean back)",
    C7 = "a new estimator (directions from the NIPALS-imputed matrix)", `H-rate` = "a guard or a message",
    `H-kappa` = "a guard or a message (kappa-hat computed before the call)")
COMPLETE_DATA <- c(S = "none", C3 = "none", C4a = "none", C4b = "a constant per feature before level restoration",
    C6 = "a constant per feature before level restoration", C7 = "none", `H-rate` = "none", `H-kappa` = "none")

writeTsv <- function(x, file) utils::write.table(x, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
readRowsIn <- function(dir, pattern) {
    files <- list.files(file.path(dir, "rows"), pattern = pattern, full.names = TRUE)
    if (!length(files)) return(NULL)
    do.call(rbind, lapply(files, function(f) utils::read.delim(gzfile(f), stringsAsFactors = FALSE, na.strings = "NA")))
}

# ---- helpers kept from score-k185.R (the limits passed in; with the defaults they compute what that file computes) --
summarise <- function(delta, upper_alpha = ALPHA, bonf_divisor = FAMILY) {
    delta <- delta[is.finite(delta)]
    n <- length(delta)
    if (n < 2L) return(c(n = n, mean = if (n) mean(delta) else NA, sd = NA, mcse = NA, upper95 = NA, lower_bonf = NA,
        lower95_2s = NA, upper95_2s = NA))
    m <- mean(delta)
    s <- stats::sd(delta)
    se <- s / sqrt(n)
    # with the default limits the expressions are score-k185.R's own, so E1'' can be byte-identical
    t_upper <- if (identical(upper_alpha, ALPHA)) stats::qt(0.95, n - 1) else stats::qt(1 - upper_alpha, n - 1)
    c(n = n, mean = m, sd = s, mcse = se, upper95 = m + t_upper * se,
        lower_bonf = m - stats::qt(1 - ALPHA / bonf_divisor, n - 1) * se,
        lower95_2s = m - stats::qt(0.975, n - 1) * se, upper95_2s = m + stats::qt(0.975, n - 1) * se)
}
classify <- function(sm, material, serious) {
    if (!is.finite(sm[["upper95"]])) return("not determined")
    if (sm[["upper95"]] < material) return("no material effect")
    if (is.finite(sm[["lower_bonf"]]) && sm[["lower_bonf"]] >= serious) return("serious")
    if (is.finite(sm[["lower_bonf"]]) && sm[["lower_bonf"]] >= material) return("material")
    "not determined"
}
clopper <- function(x, n) if (n == 0L) c(NA, NA) else c(stats::qbeta(0.025, x, n - x + 1), stats::qbeta(0.975, x + 1, n - x))
edgeOf <- function(est) EDGES[EDGES$estimand == est, , drop = FALSE]

# ---- per-replicate paired contrasts (KIT-185's, with the feature-set tag; KIT-210 6.3) --------------------------------
# `fset_b`: the tag of the B, C* and A1 rows (NULL: no tag, KIT-185's rows); `fset_da`: the tag of the DA rows.
pick <- function(rows, corr, stratum, column, fset = NULL) {
    keep <- rows$corr == corr & rows$stratum == stratum
    if (!is.null(fset) && "fset" %in% names(rows)) keep <- keep & rows$fset %in% fset
    x <- rows[keep, c("cell", "variant", "rep", "mask", column), drop = FALSE]
    names(x)[5L] <- "v"
    x
}
pairedContrasts <- function(rows, mu, S0 = "S0", O = "O", O0 = "O0", fset_b = NULL, fset_da = NULL) {
    keys <- c("cell", "variant", "rep", "mask")
    named <- function(x, name) {
        names(x)[5L] <- name
        x
    }
    joined <- function(...) Reduce(function(p, q) merge(p, q, by = keys), list(...))
    out <- list()
    for (est in c("B", "Cstar")) {
        x <- joined(named(pick(rows, mu, "all", est, fset_b), "a"), named(pick(rows, S0, "all", est, fset_b), "s"))
        out[[est]] <- data.frame(x[keys], estimand = rep(est, nrow(x)), delta = x$a - x$s, stringsAsFactors = FALSE)
    }
    x <- joined(named(pick(rows, mu, "comp", "A1", fset_b), "a"), named(pick(rows, S0, "comp", "A1", fset_b), "s"))
    out$A1 <- data.frame(x[keys], estimand = rep("A1", nrow(x)), delta = x$a - x$s, stringsAsFactors = FALSE)
    for (est in c("fdp", "power", "null_rate")) {
        sign <- if (identical(est, "power")) -1 else 1
        x <- joined(named(pick(rows, mu, "filt", est, fset_da), "a"), named(pick(rows, O, "filt", est, fset_da), "o"),
            named(pick(rows, S0, "filt", est, fset_da), "s"), named(pick(rows, O0, "filt", est, fset_da), "o0"))
        out[[est]] <- data.frame(x[keys], estimand = rep(est, nrow(x)), delta = sign * ((x$a - x$o) - (x$s - x$o0)),
            stringsAsFactors = FALSE)
    }
    do.call(rbind, out)
}

# ---- mode --k185: E1'' --------------------------------------------------------------------------------------------
if (identical(mode, "--k185")) {
    k185 <- args[[2L]]
    check_dir <- args[[3L]]
    dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)
    t0 <- proc.time()[[3L]]
    A <- readRowsIn(k185, "^a-.*[0-9]\\.tsv\\.gz$")
    reference <- utils::read.delim(file.path(k185, "results", "cells.tsv"), colClasses = "character", na.strings = "NA")
    S <- pairedContrasts(A, "S")
    C <- pairedContrasts(A, "C")
    status <- A[A$corr == "S" & A$stratum == "all" & A$variant == "base", c("cell", "mask", "rep", "status"), drop = FALSE]
    out <- list()
    for (i in seq_len(nrow(reference))) {
        ref <- reference[i, ]
        sel <- function(x) x$delta[x$cell == ref$cell & x$mask == ref$mask & x$estimand == ref$estimand & x$variant == "base"]
        tot <- summarise(sel(S))
        dC <- summarise(sel(C))
        edge <- edgeOf(ref$estimand)
        class_S <- classify(tot, edge$material, edge$serious)
        st <- status[status$cell == ref$cell & status$mask == ref$mask, , drop = FALSE]
        n_total <- as.integer(ref$R)
        refused <- sum(st$status == "refused")
        incomplete <- n_total - sum(st$status == "ok") - refused
        if (n_total > 0 && refused / n_total > 0.05) class_S <- "refused by the route" else
            if (n_total > 0 && incomplete / n_total > 0.05) class_S <- "not determined"
        out[[i]] <- data.frame(mean = tot[["mean"]], sd = tot[["sd"]], upper95 = tot[["upper95"]], lower_bonf = tot[["lower_bonf"]],
            classification = class_S, C_mean = dC[["mean"]], C_upper95 = dC[["upper95"]], C_lower_bonf = dC[["lower_bonf"]],
            C_classification = classify(dC, edge$material, edge$serious), stringsAsFactors = FALSE)
    }
    mine <- do.call(rbind, out)
    temp <- tempfile(fileext = ".tsv")
    writeTsv(mine, temp)
    mine <- utils::read.delim(temp, colClasses = "character", na.strings = "NA")
    columns <- names(mine)
    differing <- columns[!vapply(columns, function(cn) identical(mine[[cn]], reference[[cn]]), logical(1L))]
    writeTsv(data.frame(check = "E1''", rows = nrow(mine), columns_compared = paste(columns, collapse = ","),
        pass = !length(differing), differing = paste(differing, collapse = ","), seconds = proc.time()[[3L]] - t0),
        file.path(check_dir, "check-e1pp.tsv"))
    quit(save = "no")
}

# ---- mode --compare-e1: E1 --------------------------------------------------------------------------------------------
K185_COLUMNS <- c("cell", "part", "variant", "rep", "mask", "mech", "rate", "corr", "stratum", "status", "class", "k", "n_filt",
    "n_set", "n_comp", "removed", "miss_share", "miss_run_min", "miss_run_med", "miss_run_max", "n_rd_lt", "n_rd_eq",
    "n_rd_gt", "block_share", "rank", "rho", "rho_min", "A1", "inc_share", "L", "B", "RSS", "B2", "slope", "Cstar",
    "n_tested", "n_true", "n_called", "fdp", "power", "null_rate", "qmethod", "da_warn", "da_failed", "lost",
    "lost_planted", "lost_null", "k_chosen", "n_controls", "percentage")
readAll <- function(dir) {
    files <- list.files(dir, pattern = "\\.tsv\\.gz$", full.names = TRUE, recursive = TRUE)
    if (!length(files)) return(NULL)
    do.call(rbind, lapply(files, function(f) utils::read.delim(gzfile(f), colClasses = "character", na.strings = "NA")))
}
sortedOn <- function(x, columns) {
    x <- x[, columns, drop = FALSE]
    x <- x[do.call(order, x[c("cell", "part", "variant", "rep", "mask", "corr", "stratum", "class")]), , drop = FALSE]
    rownames(x) <- NULL
    x
}
if (identical(mode, "--compare-e1")) {
    a <- readAll(args[[2L]])
    b <- readAll(args[[3L]])
    check_dir <- args[[4L]]
    dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)
    rows <- list()
    for (cellname in sort(unique(c(a$cell, b$cell)))) {
        x <- sortedOn(a[a$cell == cellname, , drop = FALSE], K185_COLUMNS)
        y <- sortedOn(b[b$cell == cellname, , drop = FALSE], K185_COLUMNS)
        same_rows <- nrow(x) == nrow(y)
        differing <- if (same_rows) K185_COLUMNS[!vapply(K185_COLUMNS, function(cn) identical(x[[cn]], y[[cn]]), logical(1L))] else "row count"
        rows[[length(rows) + 1L]] <- data.frame(check = "E1", cell = cellname, rows_185 = nrow(x), rows_210 = nrow(y),
            pass = same_rows && !length(differing), differing = paste(differing, collapse = ","), stringsAsFactors = FALSE)
    }
    writeTsv(do.call(rbind, rows), file.path(check_dir, "check-e1.tsv"))
    quit(save = "no")
}

# ---- mode --compare-e10: E10 ------------------------------------------------------------------------------------------
if (identical(mode, "--compare-e10")) {
    a <- readAll(args[[2L]])
    k185 <- args[[3L]]
    check_dir <- args[[4L]]
    dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)
    deciding_columns <- c("n_filt", "removed", "B", "Cstar", "A1")
    expected_columns <- c("fdp", "power", "n_called", "null_rate", "qmethod")
    rows <- list()
    for (cellname in sort(unique(a$cell))) {
        files <- list.files(file.path(k185, "rows"), pattern = sprintf("^a-%s-[0-9]+-[0-9]+\\.tsv\\.gz$", cellname), full.names = TRUE)
        if (!length(files)) {
            rows[[length(rows) + 1L]] <- data.frame(check = "E10", cell = cellname, rows_matched = 0L, rows_bound = NA, rows_k185 = 0L,
                pass = NA, differing = "no KIT-185 rows for this cell", differing_expected_kit201 = "", stringsAsFactors = FALSE)
            next
        }
        b <- do.call(rbind, lapply(files, function(f) utils::read.delim(gzfile(f), colClasses = "character", na.strings = "NA")))
        x <- a[a$cell == cellname & a$variant == "base", , drop = FALSE]
        b <- b[b$variant == "base" & b$rep %in% unique(x$rep), , drop = FALSE]
        keys <- c("cell", "variant", "rep", "mask", "corr", "stratum")
        m <- merge(x[, c(keys, deciding_columns, expected_columns)], b[, c(keys, deciding_columns, expected_columns)], by = keys,
            suffixes = c(".bound", ".k185"))
        same <- function(cn) identical(m[[paste0(cn, ".bound")]], m[[paste0(cn, ".k185")]])
        differing <- deciding_columns[!vapply(deciding_columns, same, logical(1L))]
        expected <- expected_columns[!vapply(expected_columns, same, logical(1L))]
        rows[[length(rows) + 1L]] <- data.frame(check = "E10", cell = cellname, rows_matched = nrow(m), rows_bound = nrow(x),
            rows_k185 = nrow(b), pass = nrow(m) == nrow(x) && nrow(m) == nrow(b) && !length(differing),
            differing = paste(differing, collapse = ","), differing_expected_kit201 = paste(expected, collapse = ","), stringsAsFactors = FALSE)
    }
    writeTsv(do.call(rbind, rows), file.path(check_dir, "check-e10.tsv"))
    quit(save = "no")
}

# ---- the scoring (KIT-210 6.3 to 6.10) ----------------------------------------------------------------------------
# One correction's per-unit summaries: Delta rows (cell, variant, rep, mask, estimand, delta) split once by unit and
# estimand.
unitKey <- function(x) paste(x$cell, x$variant, x$mask, x$estimand, sep = "|")
summariseByUnit <- function(contrasts, upper_alpha = ALPHA, bonf_divisor = FAMILY, reps = NULL) {
    if (!is.null(reps)) contrasts <- contrasts[contrasts$rep %in% reps, , drop = FALSE]
    if (is.null(contrasts) || !nrow(contrasts)) {
        return(data.frame(cell = character(), variant = character(), mask = character(), estimand = character(), n = numeric(),
            mean = numeric(), sd = numeric(), mcse = numeric(), upper95 = numeric(), lower_bonf = numeric(), lower95_2s = numeric(),
            upper95_2s = numeric(), stringsAsFactors = FALSE))
    }
    groups <- split(contrasts$delta, unitKey(contrasts))
    parts <- do.call(rbind, strsplit(names(groups), "|", fixed = TRUE))
    sm <- t(vapply(groups, summarise, numeric(8L), upper_alpha = upper_alpha, bonf_divisor = bonf_divisor))
    data.frame(cell = parts[, 1L], variant = parts[, 2L], mask = parts[, 3L], estimand = parts[, 4L], sm, stringsAsFactors = FALSE,
        row.names = NULL)
}
classifyRows <- function(sm) {
    vapply(seq_len(nrow(sm)), function(i) {
        edge <- edgeOf(sm$estimand[[i]])
        classify(unlist(sm[i, c("n", "mean", "sd", "mcse", "upper95", "lower_bonf", "lower95_2s", "upper95_2s")]), edge$material, edge$serious)
    }, character(1L))
}
isDecidingUnit <- function(cell, mask) cell %in% DECIDING & ifelse(cell %in% D5_CELLS, mask %in% D5_MASKS, mask %in% BASE_MASKS)

if (identical(mode, "score") || identical(mode, "--scale")) {
    scoreAll <- function(records, out_dir, safe_rate_file) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        A <- readRowsIn(records, "^a-.*[0-9]\\.tsv\\.gz$")
        B <- readRowsIn(records, "^b-.*[0-9]\\.tsv\\.gz$")
        Sv <- readRowsIn(records, "^s-.*[0-9]\\.tsv\\.gz$")
        R <- readRowsIn(records, "^r-.*[0-9]\\.tsv\\.gz$")
        D <- readRowsIn(records, "^d-.*[0-9]\\.tsv\\.gz$")
        base <- A[A$variant == "base", , drop = FALSE]

        # status per unit and correction: S from its "all" row (fset C3), the candidates from their DA rows (fset C3)
        statusOf <- function(rows, corr, stratum) {
            x <- rows[rows$corr == corr & rows$stratum == stratum & rows$fset %in% "C3" & rows$mask != "none", , drop = FALSE]
            x[, c("cell", "variant", "rep", "mask", "status")]
        }
        miss <- base[base$corr == "mask", c("cell", "variant", "rep", "mask", "miss_share")]
        realised <- stats::aggregate(miss_share ~ cell + mask, data = miss, FUN = mean)

        # Delta of every correction on its own feature set, and D against S on the same set (6.3)
        delta <- list(S = pairedContrasts(base, "S", fset_b = "C3", fset_da = "C3"))
        D <- list()
        for (c in CANDIDATES) {
            delta[[c]] <- pairedContrasts(base, c, fset_b = c, fset_da = "C3")
            s_on_c <- pairedContrasts(base, "S", fset_b = c, fset_da = "C3")
            keys <- c("cell", "variant", "rep", "mask", "estimand")
            m <- merge(delta[[c]], s_on_c, by = keys, suffixes = c(".c", ".s"))
            D[[c]] <- data.frame(m[keys], delta = m$delta.c - m$delta.s, stringsAsFactors = FALSE)
        }
        # excess features lost (6.3): lost(c) - lost(S) on F_m^filt, per replicate; deciding on the declared class
        lost <- base[base$stratum == "filt" & base$fset %in% "C3" & base$mask != "none", c("cell", "variant", "rep", "mask", "corr", "lost")]
        for (c in c("S", CANDIDATES)) {
            x <- merge(lost[lost$corr == c, -5L], lost[lost$corr == "S", -5L], by = c("cell", "variant", "rep", "mask"), suffixes = c(".c", ".s"))
            excess <- data.frame(x[c("cell", "variant", "rep", "mask")], estimand = "excess_lost", delta = x$lost.c - x$lost.s,
                stringsAsFactors = FALSE)
            delta[[c]] <- rbind(delta[[c]], excess)
            if (!identical(c, "S")) D[[c]] <- rbind(D[[c]], excess)
        }

        # per unit classifications and J per route class (6.4); refusals: S by the route, a candidate by its guard
        J <- list()
        cells_out <- list()
        refusalShare <- function(c) {
            st <- if (identical(c, "S")) statusOf(base, "S", "all") else statusOf(base, c, "filt")
            route <- statusOf(base, "S", "all")
            st <- merge(route, st, by = c("cell", "variant", "rep", "mask"), all.x = TRUE, suffixes = c(".s", ".c"))
            st$refused <- st$status.s == "refused" | (!identical(c, "S") & st$status.c %in% "refused")
            st$incomplete <- !st$refused & !(st$status.s == "ok" & (identical(c, "S") | st$status.c %in% "ok"))
            stats::aggregate(cbind(refused, incomplete) ~ cell + variant + mask, data = st, FUN = mean)
        }
        for (c in c("S", CANDIDATES)) {
            sm <- summariseByUnit(delta[[c]])
            sm$classification <- classifyRows(sm)
            ref <- refusalShare(c)
            sm <- merge(sm, ref, by = c("cell", "variant", "mask"), all.x = TRUE)
            sm$classification[sm$incomplete %in% NA | sm$incomplete > 0.05] <- "not determined"
            sm$candidate <- c
            sm$deciding <- isDecidingUnit(sm$cell, sm$mask) & sm$variant == "base"
            sm <- merge(sm, realised, by = c("cell", "mask"), all.x = TRUE)
            cells_out[[c]] <- sm
            for (rho in names(ROUTE_CLASSES)) {
                x <- sm[sm$estimand %in% ROUTE_CLASSES[[rho]] & sm$variant == "base", , drop = FALSE]
                units <- unique(x[, c("cell", "mask", "refused", "miss_share")])
                j <- vapply(seq_len(nrow(units)), function(i) {
                    y <- x[x$cell == units$cell[[i]] & x$mask == units$mask[[i]], , drop = FALSE]
                    if (isTRUE(units$refused[[i]] > 0.05)) return("refused")
                    if (any(y$classification == "serious")) return("harmful")
                    if (length(y$classification) == length(ROUTE_CLASSES[[rho]]) && all(y$classification == "no material effect")) return("safe")
                    "not safe"
                }, character(1L))
                J[[paste(c, rho)]] <- data.frame(candidate = c, route_class = rho, cell = units$cell, mask = units$mask,
                    mech = sub("[0-9.]+$", "", units$mask), rate = as.numeric(sub("^[a-z]+", "", units$mask)) / 100,
                    realised_share = units$miss_share, J = j, deciding = isDecidingUnit(units$cell, units$mask), stringsAsFactors = FALSE)
            }
        }
        cells_tab <- do.call(rbind, cells_out)
        J_tab <- do.call(rbind, J)
        writeTsv(cells_tab, file.path(out_dir, "cells.tsv"))
        writeTsv(J_tab, file.path(out_dir, "J.tsv"))

        # safe rates per (candidate, route class, cell, mechanism) and the more severe mechanism (6.4, VAL-143)
        safe <- list()
        for (key in unique(paste(J_tab$candidate, J_tab$route_class, J_tab$cell, sep = "|"))) {
            parts <- strsplit(key, "|", fixed = TRUE)[[1L]]
            x <- J_tab[J_tab$candidate == parts[[1L]] & J_tab$route_class == parts[[2L]] & J_tab$cell == parts[[3L]] & J_tab$deciding, , drop = FALSE]
            for (mech in unique(x$mech)) {
                y <- x[x$mech == mech, , drop = FALSE]
                y <- y[order(y$rate), , drop = FALSE]
                ok <- y$J == "safe"
                safe_rate <- if (!length(ok) || !ok[[1L]]) 0 else y$rate[max(which(cumprod(ok) == 1))]
                safe[[length(safe) + 1L]] <- data.frame(candidate = parts[[1L]], route_class = parts[[2L]], cell = parts[[3L]], mech = mech,
                    safe_rate = safe_rate, realised_at_safe = if (safe_rate > 0) y$realised_share[y$rate == safe_rate][[1L]] else 0,
                    stringsAsFactors = FALSE)
            }
        }
        safe_tab <- do.call(rbind, safe)
        severe <- do.call(rbind, lapply(split(safe_tab, paste(safe_tab$candidate, safe_tab$route_class, safe_tab$cell)), function(x) {
            x[which.min(x$safe_rate), , drop = FALSE]
        }))
        writeTsv(safe_tab, file.path(out_dir, "safe-rate.tsv"))
        writeTsv(severe, file.path(out_dir, "safe-rate-severe.tsv"))

        # complete-data departure (6.3, rule 3 (b)): e(c, none) - e(S, none), family 56
        departure <- list()
        none_rows <- base[base$mask == "none", , drop = FALSE]
        for (c in CANDIDATES) {
            for (est in c("B", "Cstar")) {
                x <- merge(pick(none_rows, c, "all", est, c), pick(none_rows, "S", "all", est, c), by = c("cell", "variant", "rep", "mask"))
                departure[[length(departure) + 1L]] <- data.frame(candidate = c, x[c("cell", "variant", "rep")], estimand = est, delta = x$v.x - x$v.y)
            }
            for (est in c("fdp", "power")) {
                sign <- if (identical(est, "power")) -1 else 1
                x <- merge(pick(none_rows, c, "filt", est, "C3"), pick(none_rows, "S", "filt", est, "C3"), by = c("cell", "variant", "rep", "mask"))
                departure[[length(departure) + 1L]] <- data.frame(candidate = c, x[c("cell", "variant", "rep")], estimand = est, delta = sign * (x$v.x - x$v.y))
            }
        }
        departure <- do.call(rbind, departure)
        departure$mask <- "none"
        dep_tab <- do.call(rbind, lapply(split(departure, paste(departure$candidate, departure$cell, departure$estimand)), function(x) {
            sm <- summarise(x$delta, bonf_divisor = DEPARTURE_FAMILY)
            edge <- edgeOf(x$estimand[[1L]])
            data.frame(candidate = x$candidate[[1L]], cell = x$cell[[1L]], estimand = x$estimand[[1L]], t(sm),
                classification = classify(sm, edge$material, edge$serious), deciding = x$cell[[1L]] %in% DECIDING, stringsAsFactors = FALSE)
        }))
        writeTsv(dep_tab, file.path(out_dir, "departure.tsv"))

        # the decision rules (6.4), unadjusted and with the across-candidate adjustment (6.5)
        verdictOf <- function(c, rho, upper_alpha, bonf_divisor) {
            est <- ROUTE_CLASSES[[rho]]
            js <- J_tab[J_tab$route_class == rho & J_tab$deciding, , drop = FALSE]
            sev <- severe[severe$route_class == rho & severe$cell %in% DECIDING, , drop = FALSE]
            rule1 <- function(x) any(sev$safe_rate[sev$candidate == x] > 0)
            if (identical(c, "S")) return(list(rule1 = rule1("S"), verdict = if (rule1("S")) "safe under stated conditions" else "not safe under stated conditions"))
            d <- D[[c]][D[[c]]$estimand %in% est & isDecidingUnit(D[[c]]$cell, D[[c]]$mask) & D[[c]]$variant == "base", , drop = FALSE]
            sd_pos <- summariseByUnit(d, upper_alpha, bonf_divisor)
            d_neg <- d
            d_neg$delta <- -d_neg$delta
            sd_neg <- summariseByUnit(d_neg, upper_alpha, bonf_divisor)
            material <- vapply(sd_pos$estimand, function(e) edgeOf(e)$material, numeric(1L))
            worse <- sd_pos[is.finite(sd_pos$lower_bonf) & sd_pos$lower_bonf >= material, , drop = FALSE]
            better <- sd_neg[is.finite(sd_neg$lower_bonf) & sd_neg$lower_bonf >= material, , drop = FALSE]
            safe_S <- js[js$candidate == "S" & js$J == "safe", c("cell", "mask")]
            in_safe <- paste(sd_pos$cell, sd_pos$mask) %in% paste(safe_S$cell, safe_S$mask)
            non_inferior <- all(is.finite(sd_pos$upper95[in_safe]) & sd_pos$upper95[in_safe] < material[in_safe])
            sev_c <- sev[sev$candidate == c, , drop = FALSE]
            sev_s <- sev[sev$candidate == "S", , drop = FALSE]
            m <- merge(sev_c, sev_s, by = "cell", suffixes = c(".c", ".s"))
            rate_ok <- nrow(m) > 0 && all(m$safe_rate.c >= m$safe_rate.s)
            r2 <- c(a = !nrow(worse), b = non_inferior, c = nrow(better) > 0, d = rate_ok)
            jc <- js[js$candidate == c, c("cell", "mask", "J")]
            js_s <- js[js$candidate == "S", c("cell", "mask", "J")]
            jj <- merge(jc, js_s, by = c("cell", "mask"), suffixes = c(".c", ".s"))
            r3a <- any(jj$J.c == "harmful" & jj$J.s == "safe")
            dep <- dep_tab[dep_tab$candidate == c & dep_tab$deciding & dep_tab$estimand %in% c("B", "Cstar", "fdp", "power"), , drop = FALSE]
            r3b <- any(dep$classification %in% c("material", "serious"))
            r3c <- any(jj$J.c == "refused" & jj$J.s == "safe")
            verdict <- if (r3a || r3b || r3c) "not adoptable" else if (all(r2)) "better than the shipped" else if (r2[["c"]]) "trade-off" else "no gain"
            list(rule1 = rule1(c), rule2 = r2, rule3 = c(a = r3a, b = r3b, c = r3c), verdict = verdict,
                worse = paste(unique(paste(worse$cell, worse$mask, worse$estimand)), collapse = "; "),
                better = paste(unique(paste(better$cell, better$mask, better$estimand)), collapse = "; "))
        }
        rules <- list()
        for (rho in names(ROUTE_CLASSES)) for (c in c("S", CANDIDATES)) {
            v <- verdictOf(c, rho, ALPHA, FAMILY)
            vk <- verdictOf(c, rho, ALPHA / K_FAMILIES, FAMILY * K_FAMILIES)
            rules[[length(rules) + 1L]] <- data.frame(candidate = c, route_class = rho, rule1 = v$rule1, rule1_robust = vk$rule1,
                rule2 = paste(names(v$rule2 %||% NULL), v$rule2 %||% NULL, sep = "=", collapse = ","),
                rule3 = paste(names(v$rule3 %||% NULL), v$rule3 %||% NULL, sep = "=", collapse = ","), verdict = v$verdict,
                verdict_robust = vk$verdict, robust = identical(v$verdict, vk$verdict) && identical(v$rule1, vk$rule1),
                materially_better = v$better %||% NA, materially_worse = v$worse %||% NA, stringsAsFactors = FALSE)
        }
        rules_tab <- do.call(rbind, rules)
        writeTsv(rules_tab, file.path(out_dir, "rules.tsv"))

        # admission (6.6, H): the data set's admission values; a lipid data set is admitted when all its assays are
        assay_rows <- base[base$corr == "assay" & base$mask != "none", c("cell", "variant", "rep", "mask", "assay", "miss_share", "kappa_hat")]
        admission <- stats::aggregate(cbind(miss_share, kappa_hat) ~ cell + variant + rep + mask, data = assay_rows, FUN = max,
            na.action = stats::na.pass)
        delta_S <- delta$S
        validateAdmitted <- function(admitted, rho, reps = NULL, upper_alpha = ALPHA) {
            est <- ROUTE_CLASSES[[rho]]
            x <- merge(delta_S[delta_S$estimand %in% est, ], admitted, by = c("cell", "variant", "rep", "mask"))
            if (!is.null(reps)) x <- x[x$rep %in% reps, , drop = FALSE]
            if (!nrow(x)) return(list(validated = NA, units = 0L, units_failing = 0L, failing = "", pooled_units = 0L, pooled_ok = NA,
                admitted = 0L))
            counts <- stats::aggregate(rep ~ cell + mask, data = unique(x[, c("cell", "mask", "rep")]), FUN = length)
            big <- counts[counts$rep >= MIN_ADMITTED, c("cell", "mask")]
            small <- counts[counts$rep < MIN_ADMITTED, c("cell", "mask")]
            perUnit <- function(y) {
                vapply(est, function(e) {
                    sm <- summarise(y$delta[y$estimand == e], upper_alpha = upper_alpha)
                    is.finite(sm[["upper95"]]) && sm[["upper95"]] < edgeOf(e)$material
                }, logical(1L))
            }
            unit_ok <- vapply(seq_len(nrow(big)), function(i) all(perUnit(x[x$cell == big$cell[[i]] & x$mask == big$mask[[i]], ])), logical(1L))
            pooled <- x[paste(x$cell, x$mask) %in% paste(small$cell, small$mask), , drop = FALSE]
            pooled_ok <- if (nrow(pooled)) all(perUnit(pooled)) else TRUE
            list(validated = all(unit_ok) && pooled_ok, units = nrow(big), units_failing = sum(!unit_ok),
                failing = paste(paste(big$cell, big$mask)[!unit_ok], collapse = "; "), pooled_units = nrow(small), pooled_ok = pooled_ok,
                admitted = nrow(unique(x[, c("cell", "variant", "rep", "mask")])))
        }
        total_sets <- nrow(unique(admission[isDecidingUnit(admission$cell, admission$mask), c("cell", "variant", "rep", "mask")]))
        # H-rate: the matched cell's realised share at its KIT-185 safe rate at the more severe mechanism
        k185_safe <- utils::read.delim(safe_rate_file, stringsAsFactors = FALSE)
        thresholds <- do.call(rbind, lapply(split(k185_safe, k185_safe$cell), function(x) x[which.min(x$safe_rate), c("cell", "realised_at_safe")]))
        adm_rate <- merge(admission, thresholds, by = "cell")
        adm_rate <- adm_rate[isDecidingUnit(adm_rate$cell, adm_rate$mask) & adm_rate$miss_share <= adm_rate$realised_at_safe, c("cell", "variant", "rep", "mask")]
        h <- list()
        for (rho in names(ROUTE_CLASSES)) {
            v <- validateAdmitted(adm_rate, rho)
            h[[length(h) + 1L]] <- data.frame(rule = "H-rate", route_class = rho, validated = v$validated, units = v$units,
                units_failing = v$units_failing, failing = v$failing, pooled_units = v$pooled_units, pooled_ok = v$pooled_ok,
                admitted_share = v$admitted / total_sets, held_share = 1 - v$admitted / total_sets, kappa_star = NA, stringsAsFactors = FALSE)
        }
        # H-kappa: kappa* fitted on replicates 1-200 (6.6), with the margin 1 - 0.05/K; validated on 201-400 and out of design
        fitKappa <- function(rho, cells_in) {
            train <- admission[admission$rep <= 200 & admission$cell %in% cells_in & isDecidingUnit(admission$cell, admission$mask), , drop = FALSE]
            thresholds <- unique(stats::quantile(train$kappa_hat, probs = seq(0.02, 1, by = 0.02), na.rm = TRUE, names = FALSE))
            thresholds <- thresholds[is.finite(thresholds)]
            best <- NA_real_
            for (t in sort(thresholds)) {
                adm <- train[is.finite(train$kappa_hat) & train$kappa_hat <= t, c("cell", "variant", "rep", "mask")]
                est <- ROUTE_CLASSES[[rho]]
                x <- merge(delta_S[delta_S$estimand %in% est, ], adm, by = c("cell", "variant", "rep", "mask"))
                if (!nrow(x)) next
                counts <- stats::aggregate(rep ~ cell + mask, data = unique(x[, c("cell", "mask", "rep")]), FUN = length)
                big <- counts[counts$rep >= MIN_ADMITTED, ]
                ok <- all(vapply(seq_len(nrow(big)), function(i) {
                    y <- x[x$cell == big$cell[[i]] & x$mask == big$mask[[i]], ]
                    all(vapply(est, function(e) {
                        sm <- summarise(y$delta[y$estimand == e], upper_alpha = ALPHA / K_FAMILIES)
                        is.finite(sm[["upper95"]]) && sm[["upper95"]] < edgeOf(e)$material
                    }, logical(1L)))
                }, logical(1L)))
                if (ok) best <- t
            }
            best
        }
        loso <- list()
        for (rho in names(ROUTE_CLASSES)) {
            kstar <- fitKappa(rho, DECIDING)
            if (is.na(kstar)) {
                h[[length(h) + 1L]] <- data.frame(rule = "H-kappa", route_class = rho, validated = NA, units = NA, units_failing = NA,
                    failing = "kappa* none: not supported", pooled_units = NA, pooled_ok = NA, admitted_share = NA, held_share = NA,
                    kappa_star = NA, stringsAsFactors = FALSE)
                next
            }
            test <- admission[admission$rep > 200 & isDecidingUnit(admission$cell, admission$mask) & is.finite(admission$kappa_hat) &
                admission$kappa_hat <= kstar, c("cell", "variant", "rep", "mask")]
            v <- validateAdmitted(test, rho)
            test_total <- nrow(unique(admission[admission$rep > 200 & isDecidingUnit(admission$cell, admission$mask), c("cell", "variant", "rep", "mask")]))
            h[[length(h) + 1L]] <- data.frame(rule = "H-kappa", route_class = rho, validated = v$validated, units = v$units,
                units_failing = v$units_failing, failing = v$failing, pooled_units = v$pooled_units, pooled_ok = v$pooled_ok,
                admitted_share = v$admitted / test_total, held_share = 1 - v$admitted / test_total, kappa_star = kstar, stringsAsFactors = FALSE)
            # out of design: P0lo, P0hi (part A), the sensitivity variants and the run-aligned masks, each with its own R
            for (part_name in c("P0lo/P0hi", "sensitivity", "run-aligned")) {
                rows <- switch(part_name, `P0lo/P0hi` = A[A$cell %in% c("P0lo", "P0hi"), , drop = FALSE], sensitivity = Sv, `run-aligned` = R)
                if (is.null(rows) || !nrow(rows)) next
                ds <- pairedContrasts(rows, "S", fset_b = "C3", fset_da = "C3")
                ar <- rows[rows$corr == "assay" & rows$mask != "none", c("cell", "variant", "rep", "mask", "kappa_hat")]
                ad <- stats::aggregate(kappa_hat ~ cell + variant + rep + mask, data = ar, FUN = max, na.action = stats::na.pass)
                ad <- ad[is.finite(ad$kappa_hat) & ad$kappa_hat <= kstar, c("cell", "variant", "rep", "mask")]
                old <- delta_S
                delta_S <- ds
                v <- validateAdmitted(ad, rho)
                delta_S <- old
                h[[length(h) + 1L]] <- data.frame(rule = paste("H-kappa out of design:", part_name), route_class = rho, validated = v$validated,
                    units = v$units, units_failing = v$units_failing, failing = v$failing, pooled_units = v$pooled_units, pooled_ok = v$pooled_ok,
                    admitted_share = NA, held_share = NA, kappa_star = kstar, stringsAsFactors = FALSE)
            }
            # leave one shape out (reported)
            for (shape in DECIDING) {
                kk <- fitKappa(rho, setdiff(DECIDING, shape))
                held <- admission[admission$cell == shape & isDecidingUnit(admission$cell, admission$mask), , drop = FALSE]
                adm <- held[is.finite(held$kappa_hat) & !is.na(kk) & held$kappa_hat <= kk, c("cell", "variant", "rep", "mask")]
                v <- if (nrow(adm)) validateAdmitted(adm, rho) else list(validated = NA, units = 0, units_failing = 0)
                loso[[length(loso) + 1L]] <- data.frame(route_class = rho, shape = shape, kappa_star = kk,
                    admitted_share = if (nrow(held)) nrow(adm) / nrow(held) else NA, admitted_units_not_safe = v$units_failing,
                    validated = v$validated, stringsAsFactors = FALSE)
            }
        }
        h_tab <- do.call(rbind, h)
        writeTsv(h_tab, file.path(out_dir, "h-admission.tsv"))
        loso_tab <- if (length(loso)) do.call(rbind, loso)
        if (!is.null(loso_tab)) {
            writeTsv(loso_tab, file.path(out_dir, "loso.tsv"))
        }

        # secondary reading (6.4, K210-P17): legacy J without FDP; FDP harm shared by C3, C4a and C4b
        fdp <- cells_tab[cells_tab$estimand == "fdp" & cells_tab$deciding, c("candidate", "cell", "mask", "classification")]
        harm <- function(c) paste(fdp$cell, fdp$mask)[fdp$candidate == c & fdp$classification %in% c("material", "serious")]
        shared_units <- Reduce(intersect, lapply(c("C3", "C4a", "C4b"), harm))
        secondary <- do.call(rbind, lapply(c("S", CANDIDATES), function(c) {
            x <- cells_tab[cells_tab$candidate == c & cells_tab$deciding & cells_tab$estimand %in% c("B", "Cstar", "power"), , drop = FALSE]
            units <- unique(x[, c("cell", "mask")])
            data.frame(candidate = c, cell = units$cell, mask = units$mask,
                J_without_fdp = vapply(seq_len(nrow(units)), function(i) {
                    y <- x$classification[x$cell == units$cell[[i]] & x$mask == units$mask[[i]]]
                    if (any(y == "serious")) "harmful" else if (all(y == "no material effect")) "safe" else "not safe"
                }, character(1L)),
                fdp_harm = paste(units$cell, units$mask) %in% harm(c),
                fdp_harm_shared_by_C3_C4a_C4b = paste(units$cell, units$mask) %in% shared_units, stringsAsFactors = FALSE)
        }))
        writeTsv(secondary, file.path(out_dir, "secondary.tsv"))

        # reported per candidate and unit: agreement with S, corrected-value change, losses, guard refusals, seconds, FDP by
        # residual-dimension stratum (6.3)
        filt <- base[base$fset %in% "C3" & base$stratum %in% c("filt", "rd_lt", "rd_eq", "rd_gt"), , drop = FALSE]
        num <- function(v) suppressWarnings(as.numeric(v))
        reported_c <- stats::aggregate(cbind(jaccard_S = num(filt$jaccard_S), mad_S = num(filt$mad_S), lost = num(filt$lost),
            guard_targets = num(filt$guard_targets), fdp = num(filt$fdp), power = num(filt$power)) ~ corr + cell + mask + stratum,
            data = filt, FUN = function(v) mean(v, na.rm = TRUE), na.action = stats::na.pass)
        seconds <- stats::aggregate(seconds ~ corr + cell, data = base[base$stratum == "all" & is.finite(num(base$seconds)), ], FUN = mean)
        writeTsv(reported_c, file.path(out_dir, "reported-candidates.tsv"))
        writeTsv(seconds, file.path(out_dir, "seconds.tsv"))

        # reported parts with every candidate (6.8)
        reportedSummary <- function(contrasts, label) {
            if (is.null(contrasts) || !nrow(contrasts)) return(NULL)
            sm <- summariseByUnit(contrasts)
            data.frame(part = label, sm, stringsAsFactors = FALSE)
        }
        reported <- list()
        for (c in c("S", CANDIDATES, REPORTED_CANDIDATES)) {
            fb <- if (identical(c, "S")) "C3" else c
            reported[[length(reported) + 1L]] <- reportedSummary(pairedContrasts(A[A$cell %in% c("P0lo", "P0hi"), ], c, fset_b = fb, fset_da = "C3"), paste("P0lo/P0hi", c))
            if (!is.null(R)) reported[[length(reported) + 1L]] <- reportedSummary(pairedContrasts(R, c, fset_b = fb, fset_da = "C3"), paste("run-aligned", c))
            if (!is.null(Sv)) reported[[length(reported) + 1L]] <- reportedSummary(pairedContrasts(Sv, c, fset_b = fb, fset_da = "C3"), paste("sensitivity", c))
            if (identical(c, "C7r5")) reported[[length(reported) + 1L]] <- reportedSummary(pairedContrasts(base, c, fset_b = c, fset_da = "C3"), "part A C7r5 (reported)")
        }
        writeTsv(do.call(rbind, reported), file.path(out_dir, "reported.tsv"))
        if (!is.null(B)) {
            oracle <- A[A$variant == "base" & A$corr %in% c("O", "O0"), , drop = FALSE]
            partb <- list()
            for (route in unique(B$variant[B$corr != "mask"])) {
                tag <- if (identical(route, "legacy")) "L" else "D"
                rows <- B[B$variant == route & !B$corr %in% c("O", "O0"), , drop = FALSE]
                oracle_route <- oracle
                oracle_route$variant <- route
                both <- rbind(rows, oracle_route[, names(rows)])
                # the route's state: KIT-185's rows (fset NA); the oracle's DA rows come from part A (fset "C3")
                partb[[length(partb) + 1L]] <- reportedSummary(pairedContrasts(both, paste0("SB", tag), S0 = paste0("SB", tag, "0"),
                    fset_b = NA, fset_da = c(NA, "C3")), paste("part B", route, "S"))
                for (c in c(CANDIDATES, REPORTED_CANDIDATES)) {
                    partb[[length(partb) + 1L]] <- reportedSummary(pairedContrasts(both, paste0(c, "B", tag), S0 = paste0("SB", tag, "0"),
                        fset_b = c, fset_da = c(c, NA, "C3")), paste("part B", route, c))
                }
            }
            writeTsv(do.call(rbind, partb), file.path(out_dir, "partb.tsv"))
        }

        # the decision table (6.10) and the terminal outcome per route class (6.9)
        table_rows <- list()
        outcome <- character()
        for (rho in names(ROUTE_CLASSES)) {
            sev <- severe[severe$route_class == rho & severe$cell %in% DECIDING, , drop = FALSE]
            for (c in c("S", CANDIDATES)) {
                r <- rules_tab[rules_tab$candidate == c & rules_tab$route_class == rho, , drop = FALSE]
                sc <- sev[sev$candidate == c & sev$safe_rate > 0, , drop = FALSE]
                js <- J_tab[J_tab$candidate == c & J_tab$route_class == rho & J_tab$deciding, , drop = FALSE]
                dep <- dep_tab[dep_tab$candidate == c & dep_tab$deciding, , drop = FALSE]
                sec <- secondary[secondary$candidate == c, , drop = FALSE]
                rc <- reported_c[reported_c$corr == c & reported_c$stratum == "filt", , drop = FALSE]
                table_rows[[length(table_rows) + 1L]] <- data.frame(candidate = c, route_class = rho,
                    verdict = paste0(if (isTRUE(r$rule1)) "safe under stated conditions; " else "", r$verdict),
                    robust = r$robust,
                    safe_under = if (nrow(sc)) paste(sprintf("%s %s up to %.3g (realised %.3g), at the true k", sc$cell, sc$mech, sc$safe_rate,
                        sc$realised_at_safe), collapse = "; ") else "nowhere",
                    units_not_safe = sum(js$J != "safe"), units_harmful = sum(js$J == "harmful"), units_refused = sum(js$J == "refused"),
                    departure_material_or_serious = if (nrow(dep)) sum(dep$classification %in% c("material", "serious")) else NA,
                    mean_lost = if (nrow(rc)) mean(rc$lost, na.rm = TRUE) else NA,
                    mean_guard_targets = if (nrow(rc)) mean(rc$guard_targets, na.rm = TRUE) else NA,
                    materially_better = r$materially_better, materially_worse = r$materially_worse,
                    legacy_fdp_harm_units = sum(sec$fdp_harm), legacy_fdp_harm_shared = sum(sec$fdp_harm & sec$fdp_harm_shared_by_C3_C4a_C4b),
                    on_complete_data = COMPLETE_DATA[[c]],
                    call_agreement_with_S = if (nrow(rc)) mean(rc$jaccard_S, na.rm = TRUE) else NA,
                    median_change_from_S = if (nrow(rc)) mean(rc$mad_S, na.rm = TRUE) else NA,
                    code_change = CODE_CHANGE[[c]], stringsAsFactors = FALSE)
            }
            for (rule in c("H-rate", "H-kappa")) {
                x <- h_tab[h_tab$rule == rule & h_tab$route_class == rho, , drop = FALSE]
                table_rows[[length(table_rows) + 1L]] <- data.frame(candidate = rule, route_class = rho,
                    verdict = if (nrow(x) && isTRUE(x$validated)) "validated" else "not validated", robust = NA,
                    safe_under = if (nrow(x)) sprintf("admitted %.3g, held %.3g%s", x$admitted_share, x$held_share,
                        if (identical(rule, "H-kappa") && is.finite(x$kappa_star)) sprintf(", kappa* %.4g", x$kappa_star) else "") else NA,
                    units_not_safe = if (nrow(x)) x$units_failing else NA, units_harmful = NA, units_refused = NA,
                    departure_material_or_serious = NA, mean_lost = NA, mean_guard_targets = NA, materially_better = NA,
                    materially_worse = NA, legacy_fdp_harm_units = NA, legacy_fdp_harm_shared = NA, on_complete_data = COMPLETE_DATA[[rule]],
                    call_agreement_with_S = NA, median_change_from_S = NA, code_change = CODE_CHANGE[[rule]], stringsAsFactors = FALSE)
            }
            s_safe <- isTRUE(rules_tab$rule1_robust[rules_tab$candidate == "S" & rules_tab$route_class == rho])
            cand <- rules_tab[rules_tab$route_class == rho & rules_tab$candidate %in% CANDIDATES, , drop = FALSE]
            t_out <- if (s_safe) {
                if (any(cand$verdict_robust == "better than the shipped")) "T1" else
                    if (any(cand$verdict %in% c("trade-off", "better than the shipped"))) "T2" else "T4a"
            } else {
                if (any(cand$rule1 & cand$verdict != "not adoptable")) "T3" else "T4"
            }
            outcome <- c(outcome, sprintf("%s: %s", rho, t_out))
        }
        writeTsv(do.call(rbind, table_rows), file.path(out_dir, "decision-table.tsv"))
        writeLines(c(paste("terminal outcome", outcome), "(T5 applies instead if a gate stopped the run; the gates are not judged here)",
            sprintf("deciding units scored per candidate: %d", length(unique(paste(J_tab$cell, J_tab$mask)[J_tab$deciding])))),
            file.path(out_dir, "outcome.txt"))

        # checks: the jobs' checks and E5 (determinism)
        check_files <- list.files(file.path(records, "checks"), pattern = "\\.tsv$", full.names = TRUE)
        checks <- if (length(check_files)) do.call(rbind, lapply(check_files, function(f) cbind(utils::read.delim(f, stringsAsFactors = FALSE,
            colClasses = "character"), source = basename(f))))
        if (!is.null(D)) {
            keys <- c("cell", "variant", "rep", "mask", "corr", "stratum", "fset")
            comparable <- function(x) {
                x <- x[, setdiff(names(x), c("seconds", "part")), drop = FALSE]
                x[] <- lapply(x, as.character)
                x <- x[do.call(order, x[keys]), , drop = FALSE]
                rownames(x) <- NULL
                x
            }
            ref <- A[A$cell %in% unique(D$cell) & A$rep %in% unique(D$rep) & A$variant == "base", , drop = FALSE]
            e5 <- data.frame(check = "E5", cell = paste(unique(D$cell), collapse = ","), rep = paste(range(D$rep), collapse = "-"),
                mask = "all", value = NA, threshold = "identical without seconds", pass = identical(comparable(D), comparable(ref)),
                note = sprintf("%d rows", nrow(D)), source = "score", stringsAsFactors = FALSE)
            checks <- rbind(checks, e5[, names(checks) %||% names(e5)])
        }
        if (!is.null(checks)) writeTsv(checks, file.path(out_dir, "checks.tsv"))
        invisible(outcome)
    }
}

# ---- mode --scale: the scorer on synthetic rows of the counted run's size (lifecycle; no estimand leaves) ------------
syntheticRows <- function(dir, reps = 400L) {
    dir.create(file.path(dir, "rows"), recursive = TRUE, showWarnings = FALSE)
    set.seed(210L)
    total <- 0
    columns <- c("cell", "part", "variant", "rep", "mask", "mech", "rate", "corr", "stratum", "status", "class", "k", "n_filt",
        "n_set", "n_comp", "removed", "miss_share", "miss_run_min", "miss_run_med", "miss_run_max", "n_rd_lt", "n_rd_eq",
        "n_rd_gt", "block_share", "rank", "rho", "rho_min", "A1", "inc_share", "L", "B", "RSS", "B2", "slope", "Cstar",
        "n_tested", "n_true", "n_called", "fdp", "power", "null_rate", "qmethod", "da_warn", "da_failed", "lost",
        "lost_planted", "lost_null", "k_chosen", "n_controls", "percentage", "seconds", "fset", "assay", "kappa_hat",
        "kappa_norm", "kappa_gap", "angle_sin", "guard_targets", "jaccard_S", "mad_S", "ncomp")
    for (cellname in c(DECIDING, "P0lo", "P0hi")) {
        masks <- c("none", if (cellname %in% D5_CELLS) D5_MASKS else BASE_MASKS)
        refs_masked <- c("O", "N", "S0", "O0", "N0")
        one <- function(mask) {
            refs <- if (identical(mask, "none")) c("O", "N") else refs_masked
            b1 <- expand.grid(corr = c("S", "C3", refs), stratum = c("all", "comp", "inc"), fset = "C3", stringsAsFactors = FALSE)
            d1 <- expand.grid(corr = c("S", CANDIDATES, if (cellname %in% c("L1", "L2", "M1", "M2", "L1b", "M1b")) "C7r5", refs),
                stratum = c("filt", "rd_lt", "rd_eq", "rd_gt"), fset = "C3", stringsAsFactors = FALSE)
            bc <- do.call(rbind, lapply(setdiff(CANDIDATES, "C3"), function(c) expand.grid(corr = c("S", c, refs),
                stratum = c("all", "comp", "inc"), fset = c, stringsAsFactors = FALSE)))
            mk <- data.frame(corr = c("mask", "assay"), stratum = "all", fset = NA, stringsAsFactors = FALSE)
            x <- rbind(mk, b1, d1, bc)
            x$mask <- mask
            x
        }
        template <- do.call(rbind, lapply(masks, one))
        n <- nrow(template)
        table <- template[rep(seq_len(n), reps), , drop = FALSE]
        table$rep <- rep(seq_len(reps), each = n)
        table$cell <- cellname
        table$part <- "A"
        table$variant <- "base"
        table$status <- "ok"
        table$mech <- sub("[0-9.]+$", "", table$mask)
        table$rate <- suppressWarnings(as.numeric(sub("^[a-z]+", "", table$mask)) / 100)
        for (cn in setdiff(columns, names(table))) table[[cn]] <- stats::runif(nrow(table))
        table$class <- NA
        table$qmethod <- "BH"
        table$assay <- "pos"
        total <- total + nrow(table)
        connection <- gzfile(file.path(dir, "rows", sprintf("a-%s-1-%d.tsv.gz", cellname, reps)), "w")
        utils::write.table(table[, columns], connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
        close(connection)
    }
    total
}
peakMemory <- function() {
    status <- tryCatch(readLines("/proc/self/status"), error = function(e) character())
    line <- grep("^VmHWM", status, value = TRUE)
    if (length(line)) as.numeric(gsub("[^0-9]", "", line)) / 1024^2 else NA_real_
}

here_dir <- tryCatch(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))),
    error = function(e) ".")
k185_safe_rate <- file.path(here_dir, "..", "kit-185", "results", "safe-rate.tsv")
if (identical(mode, "--scale")) {
    private <- args[[2L]]
    check_dir <- args[[3L]]
    dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)
    t0 <- proc.time()[[3L]]
    rows <- syntheticRows(private)
    t1 <- proc.time()[[3L]]
    result <- tryCatch({
        scoreAll(private, file.path(private, "scored"), k185_safe_rate)
        "completed"
    }, error = function(e) paste("error:", substr(gsub("[\r\n\t]+", " ", conditionMessage(e)), 1L, 300L)))
    writeTsv(data.frame(check = "scorer scale test", rows = rows, generate_seconds = t1 - t0, score_seconds = proc.time()[[3L]] - t1,
        peak_gib = peakMemory(), result = result, stringsAsFactors = FALSE), file.path(check_dir, "check-scale.tsv"))
    quit(save = "no")
}
if (identical(mode, "score")) {
    if (length(args) != 2L) stop("usage: Rscript score-k210.R <records dir> <out dir>")
    outcome <- scoreAll(args[[1L]], args[[2L]], k185_safe_rate)
    cat(sprintf("KIT-210 scoring finished: %s\n", paste(outcome, collapse = "; ")))
}
