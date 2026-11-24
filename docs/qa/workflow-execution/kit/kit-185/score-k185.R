# KIT-185: scoring (contrasts, intervals, classifications, onset, safe rates, attribution, outcome) from the committed
# rows. Run once, on a runner, after the counted rows are committed (KIT-185 6, 10). Pinned by the declaration's SHA-256.
#
# Usage (from the package root): Rscript docs/qa/workflow-execution/kit/kit-185/score-k185.R <records dir> <out dir>
#   <records dir> holds rows/*.tsv.gz (parts A, B, S, R and the determinism job D) and checks/*.tsv (the jobs' checks).
# Writes <out dir>/cells.tsv, onset.tsv, safe-rate.tsv, attribution.tsv, partb.tsv, sensitivity.tsv, reported.tsv,
#   checks.tsv, outcome.txt.

args <- commandArgs(TRUE)
if (length(args) != 2L) stop("usage: Rscript score-k185.R <records dir> <out dir>")
records <- args[[1L]]
out_dir <- args[[2L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- declared constants (KIT-185 5.5 to 5.7; the owner's D1, D3, D5 (a), VAL-67) ------------------------------------
DECIDING <- c("P0", "P0r", "P1", "P2", "P3", "P4", "P5", "L1", "L2", "M1", "M2", "P4b", "L1b", "M1b")
BASE_MASKS <- c("vd0.1", "vd0.5", "vd2", "vd10", "vd25", "vd50", "mcar0.5", "mcar2", "mcar10")
D5_MASKS <- c("vd2", "vd10", "bb2", "bb10", "bb25")
FAMILY <- 11L * 9L * 4L + 3L * 5L * 4L # 456
ALPHA <- 0.05
EDGES <- data.frame(estimand = c("B", "Cstar", "fdp", "power", "A1", "null_rate"),
    material = c(0.05, 0.02, 0.010, 0.02, 0.05, 0.010), serious = c(0.25, 0.10, 0.025, 0.10, 0.25, 0.025),
    deciding = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE), stringsAsFactors = FALSE)
K_OF <- c(P0 = 2, P0r = 2, P1 = 2, P2 = 2, P3 = 2, P4 = 2, P5 = 2, L1 = 2, L2 = 2, M1 = 2, M2 = 2, P0lo = 2, P0hi = 2,
    P4b = 4, L1b = 2, M1b = 15)

readRows <- function(pattern) {
    files <- list.files(file.path(records, "rows"), pattern = pattern, full.names = TRUE)
    if (!length(files)) return(NULL)
    do.call(rbind, lapply(files, function(f) utils::read.delim(gzfile(f), stringsAsFactors = FALSE, na.strings = "NA")))
}
A <- readRows("^a-.*[0-9]\\.tsv\\.gz$")
B <- readRows("^b-.*[0-9]\\.tsv\\.gz$")
S <- readRows("^s-.*[0-9]\\.tsv\\.gz$")
R <- readRows("^r-.*[0-9]\\.tsv\\.gz$")
D <- readRows("^d-.*[0-9]\\.tsv\\.gz$")

# ---- per-replicate paired contrasts ---------------------------------------------------------------------------
value <- function(rows, corr, stratum, column) {
    x <- rows[rows$corr == corr & rows$stratum == stratum, c("cell", "variant", "rep", "mask", column), drop = FALSE]
    names(x)[5L] <- "v"
    x
}
# Delta for one correction `mu` (S, C or Z; SBL, CBL... in part B) of every masked arm, per replicate, oriented so that
# larger is worse. B and C* are already referenced to O on the same mask: mu(m) - S0(m). FDP, power (as a loss) and the
# null rate: [mu(m) - O(m)] - [S0(m) - O0(m)].
pairedContrasts <- function(rows, mu, S0 = "S0", O = "O", O0 = "O0") {
    keys <- c("cell", "variant", "rep", "mask")
    named <- function(x, name) {
        names(x)[5L] <- name
        x
    }
    joined <- function(...) Reduce(function(p, q) merge(p, q, by = keys), list(...))
    out <- list()
    for (est in c("B", "Cstar")) {
        x <- joined(named(value(rows, mu, "all", est), "a"), named(value(rows, S0, "all", est), "s"))
        out[[est]] <- data.frame(x[keys], estimand = rep(est, nrow(x)), delta = x$a - x$s, stringsAsFactors = FALSE)
    }
    x <- joined(named(value(rows, mu, "comp", "A1"), "a"), named(value(rows, S0, "comp", "A1"), "s"))
    out$A1 <- data.frame(x[keys], estimand = rep("A1", nrow(x)), delta = x$a - x$s, stringsAsFactors = FALSE)
    for (est in c("fdp", "power", "null_rate")) {
        sign <- if (identical(est, "power")) -1 else 1
        x <- joined(named(value(rows, mu, "filt", est), "a"), named(value(rows, O, "filt", est), "o"),
            named(value(rows, S0, "filt", est), "s"), named(value(rows, O0, "filt", est), "o0"))
        out[[est]] <- data.frame(x[keys], estimand = rep(est, nrow(x)), delta = sign * ((x$a - x$o) - (x$s - x$o0)),
            stringsAsFactors = FALSE)
    }
    do.call(rbind, out)
}
unitStatus <- function(rows) {
    s <- rows[rows$corr == "S" & rows$stratum == "all", c("cell", "variant", "rep", "mask", "status", "n_comp"), drop = FALSE]
    crashed <- rows[rows$corr == "unit" & rows$status == "crashed", c("cell", "variant", "rep"), drop = FALSE]
    list(S = s, crashed = crashed)
}

# ---- classification (KIT-185 5.6, 5.7) -------------------------------------------------------------------------
summarise <- function(delta) {
    delta <- delta[is.finite(delta)]
    n <- length(delta)
    if (n < 2L) return(c(n = n, mean = if (n) mean(delta) else NA, sd = NA, mcse = NA, upper95 = NA, lower_bonf = NA,
        lower95_2s = NA, upper95_2s = NA))
    m <- mean(delta)
    s <- stats::sd(delta)
    se <- s / sqrt(n)
    c(n = n, mean = m, sd = s, mcse = se, upper95 = m + stats::qt(0.95, n - 1) * se,
        lower_bonf = m - stats::qt(1 - ALPHA / FAMILY, n - 1) * se,
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

contrastsA <- if (!is.null(A)) list(S = pairedContrasts(A, "S"), C = pairedContrasts(A, "C"),
    Z = if (any(A$corr == "Z" & A$mask != "none")) pairedContrasts(A, "Z")) else list()
statusA <- if (!is.null(A)) unitStatus(A)

cellTable <- list()
for (cellname in unique(A$cell)) {
    masks <- setdiff(unique(A$mask[A$cell == cellname]), "none")
    for (mask in masks) {
        st <- statusA$S[statusA$S$cell == cellname & statusA$S$mask == mask & statusA$S$variant == "base", , drop = FALSE]
        reps <- sort(unique(A$rep[A$cell == cellname & A$variant == "base"]))
        n_total <- length(reps)
        refused <- sum(st$status == "refused")
        incomplete <- n_total - sum(st$status == "ok") - refused
        assessable <- st$rep[st$status == "ok" & st$n_comp >= max(20, 3 * K_OF[[cellname]])]
        realised <- mean(A$miss_share[A$cell == cellname & A$mask == mask & A$corr == "mask" & A$variant == "base"], na.rm = TRUE)
        for (est in EDGES$estimand) {
            edge <- EDGES[EDGES$estimand == est, ]
            pick <- function(mu) {
                x <- contrastsA[[mu]]
                if (is.null(x)) return(numeric())
                x$delta[x$cell == cellname & x$mask == mask & x$estimand == est & x$variant == "base"]
            }
            pickRep <- function(mu) {
                x <- contrastsA[[mu]]
                x <- x[x$cell == cellname & x$mask == mask & x$estimand == est & x$variant == "base", , drop = FALSE]
                stats::setNames(x$delta, x$rep)
            }
            tot <- summarise(pick("S"))
            dC <- summarise(pick("C"))
            dS <- pickRep("S")
            dCr <- pickRep("C")
            both <- intersect(intersect(names(dS), names(dCr)), as.character(assessable))
            attr_delta <- dS[both] - dCr[both]
            att <- summarise(attr_delta)
            zr <- pickRep("Z")
            zboth <- intersect(names(zr), names(dCr))
            attZ <- summarise(zr[zboth] - dCr[zboth])
            class_S <- classify(tot, edge$material, edge$serious)
            class_C <- classify(dC, edge$material, edge$serious)
            deciding <- edge$deciding && cellname %in% DECIDING &&
                mask %in% (if (cellname %in% c("P4b", "L1b", "M1b")) D5_MASKS else BASE_MASKS)
            if (n_total > 0 && refused / n_total > 0.05) class_S <- "refused by the route"
            else if (n_total > 0 && incomplete / n_total > 0.05) class_S <- "not determined"
            ci <- clopper(refused, n_total)
            t_bonf <- if (tot[["n"]] >= 2) stats::qt(1 - ALPHA / FAMILY, tot[["n"]] - 1) else NA
            t_95 <- if (tot[["n"]] >= 2) stats::qt(0.95, tot[["n"]] - 1) else NA
            share <- if (is.finite(tot[["lower95_2s"]]) && (tot[["lower95_2s"]] > 0 || tot[["upper95_2s"]] < 0) && length(both) > 1L) {
                # the attributable share with a delta-method interval (paired ratio of means)
                a <- attr_delta
                t <- dS[both]
                ratio <- mean(a) / mean(t)
                g <- c(1 / mean(t), -mean(a) / mean(t)^2)
                v <- stats::var(cbind(a, t)) / length(both)
                se <- sqrt(drop(t(g) %*% v %*% g))
                c(ratio, ratio - 1.96 * se, ratio + 1.96 * se)
            } else c(NA, NA, NA)
            cellTable[[length(cellTable) + 1L]] <- data.frame(cell = cellname, mask = mask, mech = sub("[0-9.]+$", "", mask),
                rate = as.numeric(sub("^[a-z]+", "", mask)) / 100, realised_share = realised, estimand = est, deciding = deciding,
                R = n_total, n = tot[["n"]], refused = refused, refused_share = if (n_total) refused / n_total else NA,
                refused_lo = ci[[1L]], refused_hi = ci[[2L]], did_not_complete = incomplete,
                mean = tot[["mean"]], sd = tot[["sd"]], mcse = tot[["mcse"]], upper95 = tot[["upper95"]], lower_bonf = tot[["lower_bonf"]],
                lower95_2s = tot[["lower95_2s"]], upper95_2s = tot[["upper95_2s"]],
                resolvable_material = if (is.finite(t_bonf)) edge$material + t_bonf * tot[["mcse"]] else NA,
                resolvable_none = if (is.finite(t_95)) edge$material - t_95 * tot[["mcse"]] else NA,
                classification = class_S,
                C_mean = dC[["mean"]], C_upper95 = dC[["upper95"]], C_lower_bonf = dC[["lower_bonf"]], C_classification = class_C,
                attr_n = att[["n"]], attr_not_assessable = sum(st$status == "ok") - length(assessable), attr_mean = att[["mean"]],
                attr_lo = att[["lower95_2s"]], attr_hi = att[["upper95_2s"]], attr_share = share[[1L]], attr_share_lo = share[[2L]],
                attr_share_hi = share[[3L]], attrZ_n = attZ[["n"]], attrZ_mean = attZ[["mean"]], attrZ_lo = attZ[["lower95_2s"]],
                attrZ_hi = attZ[["upper95_2s"]], stringsAsFactors = FALSE)
        }
    }
}
cells <- do.call(rbind, cellTable)

# attribution sentence (5.7), in the declared order
cells$attribution <- NA_character_
hit <- cells$classification %in% c("material", "serious")
cells$attribution[hit] <- with(cells[hit, ], ifelse(is.finite(attr_lo) & attr_lo <= 0 & attr_hi >= 0, "not attributable",
    ifelse(attr_n < 2, "not assessable", ifelse(C_classification == "no material effect", "attributable to the zero-filling",
        ifelse(C_classification == "material" & classification == "serious", "partly", "attribution not determined")))))
utils::write.table(cells, file.path(out_dir, "cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
utils::write.table(cells[hit, c("cell", "mask", "estimand", "classification", "C_classification", "attr_mean", "attr_lo", "attr_hi",
    "attr_share", "attr_share_lo", "attr_share_hi", "attrZ_mean", "attrZ_lo", "attrZ_hi", "attribution")],
    file.path(out_dir, "attribution.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# ---- onset table and safe rates (5.7) -----------------------------------------------------------------------
decide <- cells[cells$deciding, , drop = FALSE]
onset <- list()
safe <- list()
for (cellname in unique(decide$cell)) for (mech in unique(decide$mech[decide$cell == cellname])) {
    x <- decide[decide$cell == cellname & decide$mech == mech, , drop = FALSE]
    rates <- sort(unique(x$rate))
    for (est in unique(x$estimand)) {
        y <- x[x$estimand == est, , drop = FALSE]
        y <- y[order(y$rate), ]
        material <- y$rate[y$classification %in% c("material", "serious")]
        serious <- y$rate[y$classification == "serious"]
        onset[[length(onset) + 1L]] <- data.frame(cell = cellname, mech = mech, estimand = est,
            lowest_material = if (length(material)) min(material) else NA, lowest_serious = if (length(serious)) min(serious) else NA,
            realised_at_material = if (length(material)) y$realised_share[y$rate == min(material)][[1L]] else NA,
            stringsAsFactors = FALSE)
    }
    ok_rate <- vapply(rates, function(r) all(x$classification[x$rate == r] == "no material effect"), logical(1L))
    safe_rate <- if (!ok_rate[[1L]]) 0 else rates[max(which(cumprod(ok_rate) == 1))]
    safe[[length(safe) + 1L]] <- data.frame(cell = cellname, mech = mech, safe_rate = safe_rate,
        realised_at_safe = if (safe_rate > 0) x$realised_share[x$rate == safe_rate][[1L]] else 0,
        highest_measured = max(rates), stringsAsFactors = FALSE)
}
utils::write.table(do.call(rbind, onset), file.path(out_dir, "onset.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
utils::write.table(do.call(rbind, safe), file.path(out_dir, "safe-rate.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# ---- the outcome (5.7) -------------------------------------------------------------------------------------
judged <- decide$classification[decide$classification != "refused by the route"]
outcome <- if (all(judged == "no material effect")) "O1" else if (any(judged == "serious")) "O3" else
    if (any(judged == "material")) "O2" else "O4"
held <- setdiff(DECIDING, unique(decide$cell))
writeLines(c(sprintf("outcome %s", outcome), sprintf("deciding cells scored: %s", paste(unique(decide$cell), collapse = ", ")),
    sprintf("deciding cells held (not in these rows): %s", if (length(held)) paste(held, collapse = ", ") else "none"), sprintf("deciding classifications %d (family %d)", nrow(decide), FAMILY),
    paste(names(table(decide$classification)), table(decide$classification), sep = ": ")), file.path(out_dir, "outcome.txt"))

# ---- reported parts: P0lo, P0hi; run-aligned; part B; sensitivity -------------------------------------------
reportedSummary <- function(contrasts, label) {
    if (is.null(contrasts) || !nrow(contrasts)) return(NULL)
    groups <- unique(contrasts[, c("cell", "variant", "mask", "estimand")])
    do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
        g <- groups[i, ]
        d <- contrasts$delta[contrasts$cell == g$cell & contrasts$variant == g$variant & contrasts$mask == g$mask &
            contrasts$estimand == g$estimand]
        sm <- summarise(d)
        data.frame(part = label, g, t(sm), stringsAsFactors = FALSE)
    }))
}
reported <- rbind(
    reportedSummary(contrastsA$S[contrastsA$S$cell %in% c("P0lo", "P0hi"), ], "A reported S"),
    reportedSummary(contrastsA$C[contrastsA$C$cell %in% c("P0lo", "P0hi"), ], "A reported C"),
    if (!is.null(R)) reportedSummary(pairedContrasts(R, "S"), "run-aligned S"),
    if (!is.null(R)) reportedSummary(pairedContrasts(R, "C"), "run-aligned C"))
utils::write.table(reported, file.path(out_dir, "reported.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
if (!is.null(S)) {
    sens <- rbind(reportedSummary(pairedContrasts(S, "S"), "sensitivity S"), reportedSummary(pairedContrasts(S, "C"), "sensitivity C"))
    utils::write.table(sens, file.path(out_dir, "sensitivity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
if (!is.null(B)) {
    # part B: the DA contrasts referenced to part A's O rows of the same replicate and mask (B08 (v))
    oracle <- A[A$variant == "base" & A$corr %in% c("O", "O0"), , drop = FALSE]
    partb <- list()
    for (route in unique(B$variant[B$corr != "mask"])) {
        tag <- if (identical(route, "legacy")) "L" else "D"
        rows <- B[B$variant == route & !B$corr %in% c("O", "O0"), , drop = FALSE]
        oracle_route <- oracle
        oracle_route$variant <- route
        for (mu in c(paste0("SB", tag), paste0("CB", tag))) {
            x <- pairedContrasts(rbind(rows, oracle_route[, names(rows)]), mu, S0 = paste0("SB", tag, "0"))
            partb[[length(partb) + 1L]] <- reportedSummary(x, paste("part B", route, mu))
        }
    }
    utils::write.table(do.call(rbind, partb), file.path(out_dir, "partb.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

# ---- checks: the jobs' checks, E5 (determinism), E8 (regenerated replicates) ----------------------------------
check_files <- list.files(file.path(records, "checks"), pattern = "\\.tsv$", full.names = TRUE)
checks <- do.call(rbind, lapply(check_files, function(f) cbind(utils::read.delim(f, stringsAsFactors = FALSE, colClasses = "character"),
    source = basename(f))))
drop_cols <- c("seconds")
comparable <- function(x, keys) {
    x <- x[, setdiff(names(x), c(drop_cols, "part")), drop = FALSE]
    # compared as text: tables read separately can infer different column types (an all-NA column is logical in one
    # and numeric in another), which identical() would count as a difference (dated section of the declaration)
    x[] <- lapply(x, as.character)
    x <- x[do.call(order, x[keys]), , drop = FALSE]
    rownames(x) <- NULL
    x
}
extra <- list()
if (!is.null(D)) {
    keys <- c("cell", "variant", "rep", "mask", "corr", "stratum")
    ref <- A[A$cell %in% unique(D$cell) & A$rep %in% unique(D$rep) & A$variant == "base", , drop = FALSE]
    same <- identical(comparable(D, keys), comparable(ref, keys))
    extra[[length(extra) + 1L]] <- data.frame(check = "E5", cell = paste(unique(D$cell), collapse = ","), rep = paste(range(D$rep), collapse = "-"),
        mask = "all", value = NA, threshold = "identical without seconds", pass = same, note = sprintf("%d rows", nrow(D)), source = "score")
}
if (!is.null(S)) {
    keys <- c("cell", "rep", "mask", "corr", "stratum")
    for (v in intersect(c("k1", "k4"), unique(S$variant))) {
        # the N rows of the test on F_filt (the strata of F_m depend on the variant's own S and C)
        s <- S[S$variant == v & S$corr == "N" & S$stratum == "filt", , drop = FALSE]
        a <- A[A$variant == "base" & A$corr == "N" & A$stratum == "filt" &
            paste(A$cell, A$rep, A$mask) %in% paste(s$cell, s$rep, s$mask), , drop = FALSE]
        cols <- setdiff(names(s), c("variant", "k", "seconds", "part"))
        same <- identical(comparable(s[, cols], keys), comparable(a[, cols], keys))
        extra[[length(extra) + 1L]] <- data.frame(check = "E8", cell = paste(unique(s$cell), collapse = ","), rep = "all", mask = v,
            value = NA, threshold = "N rows identical", pass = same, note = sprintf("%d rows", nrow(s)), source = "score")
    }
}
if (!is.null(R)) {
    keys <- c("cell", "rep", "mask", "corr", "stratum")
    r <- R[R$mask == "none", , drop = FALSE]
    a <- A[A$variant == "base" & A$mask == "none" & paste(A$cell, A$rep) %in% paste(r$cell, r$rep), , drop = FALSE]
    cols <- setdiff(names(r), c("variant", "seconds", "part"))
    same <- identical(comparable(r[, cols], keys), comparable(a[, cols], keys))
    extra[[length(extra) + 1L]] <- data.frame(check = "E8", cell = paste(unique(r$cell), collapse = ","), rep = "all", mask = "none (run-aligned)",
        value = NA, threshold = "none rows identical", pass = same, note = sprintf("%d rows", nrow(r)), source = "score")
}
all_checks <- rbind(checks, if (length(extra)) do.call(rbind, lapply(extra, function(x) {
    x[] <- lapply(x, as.character)
    x
})))
utils::write.table(all_checks, file.path(out_dir, "checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat(sprintf("KIT-185 scoring: outcome %s; %d deciding classifications\n", outcome, nrow(decide)))
