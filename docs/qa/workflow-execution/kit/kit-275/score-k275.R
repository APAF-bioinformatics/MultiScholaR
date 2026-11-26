# KIT-275: the scoring (sections 6, 7 of the ticket), run once on a runner over the committed rows of the counted jobs. Evidence code under
# docs/qa: no package code. Base R and stats only. Hashed in the declaration.
#
# Usage:
#   Rscript score-k275.R <rows dir> <out dir>                    scoring: the rows (*.tsv.gz) and the counted jobs' checks (checks-*.tsv) in <rows dir>
#   Rscript score-k275.R --blind <rows dir> <out dir> <job>      blind (L1, the smoke): parse success, cell coverage and E8 only; no estimand
#   Rscript score-k275.R --labels <table.tsv> <out dir>          the label rules alone, on constructed intervals (E4's boundary cases)
# Environment:
#   K275_SCORE_R   "<tier D>,<others>,<ID>" data sets per cell (default 200,100,20); E4's self-test sets smaller values
#   K275_BOOT_B    bootstrap resamples (default 10000; E4's self-test sets a smaller value on constructed degenerate data)
#   K275_G1        the G1 table (columns design, kind PC2 | NC2, r, rho): required for an outcome line (F3)
# Every table is read as TEXT (`colClasses = "character"`) and converted explicitly: `read.delim` turns a column holding only the
# letters T and F (tiers T and F) into logical (lesson 1). `$` on a list never relies on a prefix (lesson 2).
# Boundaries (a limit "equal to" a margin) are compared at 1e-12, because a margin such as 0.95 is not a binary number.
# Output (scoring): cells.tsv (every estimate with its intervals and label), decision.tsv (the outcome row and the consequence labels),
# settings.txt. A row set carrying the smoke flag is refused for an outcome line.

options(warnPartialMatchDollar = TRUE)
globalCallingHandlers(warning = function(w) {
    if (grepl("partial match of", tolower(conditionMessage(w)), fixed = TRUE)) stop(paste("partial-match guard:", conditionMessage(w)), call. = FALSE)
})
args <- commandArgs(TRUE)
blind <- length(args) >= 1L && identical(args[[1L]], "--blind")
labels_mode <- length(args) >= 1L && identical(args[[1L]], "--labels")
if (labels_mode) {
    stopifnot(length(args) == 3L)
    rows_dir <- args[[2L]]
    out_dir <- args[[3L]]
    job <- "labels"
} else if (blind) {
    stopifnot(length(args) == 4L)
    rows_dir <- args[[2L]]
    out_dir <- args[[3L]]
    job <- args[[4L]]
} else {
    stopifnot(length(args) == 2L)
    rows_dir <- args[[1L]]
    out_dir <- args[[2L]]
    job <- "score"
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- constants (section 7.1) -----------------------------------------------------------------------------------------------
M_ABS <- 0.05
M_ACT <- 0.50
M_RHO <- 0.95
M_MATCH <- 0.03
M_ANY <- 0.10
M_ALPHA <- 0.0125
TOL <- 1e-12
DESIGNS <- c("D3", "D4", "D6")
PER <- c(D3 = 3, D4 = 4, D6 = 6)
CIDX <- c(D3 = 1L, D4 = 2L, D6 = 3L)
TIER_T <- c(D = 1L, K = 2L, F = 3L, S = 4L, V = 5L, M = 6L, T = 7L, E = 8L)  # the bootstrap stream t of each tier (5.3); ID has none (E8 only); G1 uses 10 (R-4)
THETAS <- c("05", "07", "09", "095")
R_SET <- as.integer(strsplit(Sys.getenv("K275_SCORE_R", "200,100,20"), ",", fixed = TRUE)[[1L]])
BOOT_B <- as.integer(Sys.getenv("K275_BOOT_B", "10000"))
declaredR <- function(tier) if (identical(tier, "D")) R_SET[[1L]] else if (identical(tier, "ID")) R_SET[[3L]] else R_SET[[2L]]
le <- function(a, b) is.finite(a) && is.finite(b) && a <= b + TOL
gt <- function(a, b) is.finite(a) && is.finite(b) && a > b + TOL
ge <- function(a, b) is.finite(a) && is.finite(b) && a >= b - TOL
lt <- function(a, b) is.finite(a) && is.finite(b) && a < b - TOL
num <- function(x) suppressWarnings(as.numeric(x))
covered <- function(n_valid, R) n_valid >= 0.90 * R
# the label rules of section 7.1 (a limit "equal to" a margin is compared at 1e-12)
labelA <- function(lw, up, n, R) if (!covered(n, R)) "not determined (coverage)" else if (le(up, M_ABS)) "absorbed" else if (gt(lw, M_ACT)) "acts" else "mixed"
labelB <- function(lo, hi, n, R) if (!covered(n, R)) "not determined (coverage)" else if (lt(hi, M_RHO)) "understated" else if (ge(lo, M_RHO)) "not understated" else "not determined"
labelBp <- function(lo, hi, n, R, derived, est) if (!covered(n, R)) "not determined (coverage)" else if (ge(lo, derived - M_MATCH) && le(hi, derived + M_MATCH)) "as derived" else
    sprintf("differs from derived (%s)", format(round(est, 4)))
labelC <- function(lo, hi, n, R, margin) if (!covered(n, R)) "not determined (coverage)" else if (gt(lo, margin)) "raised (material)" else if (le(hi, margin)) "not raised" else "not determined"
shownC <- function(lo, n, R) covered(n, R) && gt(lo, 0)
labelPi <- function(lo, hi, n = 0, R = 0) if (n > 0 && !covered(n, R)) "not determined (coverage)" else if (gt(lo, 0.5)) "imprint retained" else if (lt(hi, 0.5)) "imprint removed" else "not determined"
labelsOf <- function(table) {
    out <- lapply(seq_len(nrow(table)), function(i) {
        x <- table[i, ]
        lo <- num(x[["lo"]]); hi <- num(x[["hi"]]); n <- num(x[["n"]]); R <- num(x[["R"]])
        label <- switch(x[["kind"]], A = labelA(lo, hi, n, R), B = labelB(lo, hi, n, R), Bp = labelBp(lo, hi, n, R, num(x[["derived"]]), num(x[["est"]])),
            C = labelC(lo, hi, n, R, num(x[["margin"]])), pi = labelPi(lo, hi), stop("unknown kind"))
        data.frame(id = x[["id"]], kind = x[["kind"]], label = label, shown = if (identical(x[["kind"]], "C")) as.character(shownC(lo, n, R)) else NA_character_,
            stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
}
TEXT_COLUMNS <- c("tier", "design", "cell", "omic", "version", "arm", "smoke", "status", "class", "step", "ruv_status_y", "ruv_status_x",
    "corr_removed_y", "corr_removed_x", "tested_sig_y", "e5_y", "e5_x", "pe_a_excluded")

# ---- reading: every table as text ---------------------------------------------------------------------------------------------
readText <- function(path) utils::read.delim(path, stringsAsFactors = FALSE, colClasses = "character", na.strings = "NA")
if (labels_mode) {
    utils::write.table(labelsOf(readText(rows_dir)), file.path(out_dir, "labels.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    quit(save = "no", status = 0L)
}
row_files <- c(list.files(rows_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE, recursive = TRUE),
    list.files(rows_dir, pattern = "^part-.*\\.tsv$", full.names = TRUE, recursive = TRUE))
rows <- do.call(rbind, lapply(row_files, function(f) readText(if (grepl("\\.gz$", f)) gzfile(f) else f)))
if (is.null(rows) || !nrow(rows)) stop("K275 scoring stops: no rows")
for (column in setdiff(names(rows), TEXT_COLUMNS)) rows[[column]] <- num(rows[[column]])
instrument_rows <- sum(rows[["status"]] %in% "instrument" | grepl("partial-match guard|partial_match_guard", rows[["class"]]))
if (!blind && instrument_rows > 0L) stop("K275 scoring stops: an instrument status or a partial-match guard stop in the rows (VAL-54)")
rows <- rows[!rows[["arm"]] %in% "unit", , drop = FALSE]
if (!blind) {
    lim <- vapply(rows[["tier"]], function(t) if (t %in% c("D", "K", "F", "S", "V", "M", "T", "E", "ID")) declaredR(t) else NA_real_, numeric(1L))
    if (any(!is.finite(rows[["r"]]) | rows[["r"]] < 1 | !is.finite(lim) | rows[["r"]] > lim)) stop("K275 scoring stops: a row outside its tier's declared data sets (RANGE-1)")
    key <- paste(rows[["tier"]], rows[["design"]], rows[["version"]], rows[["arm"]], rows[["r"]])
    if (anyDuplicated(key)) stop("K275 scoring stops: a journey row appears twice")
}
VALID <- c("ok", "dropped", "uncorrected", "no_da", "refused_after_exclusion")
rows[["valid"]] <- rows[["status"]] %in% VALID
isTrue <- function(x) as.character(x) %in% "TRUE"
isFalse <- function(x) as.character(x) %in% "FALSE"

# ---- blind mode: parse success, coverage, E8; nothing of the measurement is written -------------------------------------------
armClassOf <- function(arm) {
    if (arm %in% c("PEL2", "PEG2", "LX")) return(arm)
    sub("^(LM|LG|GK)[0-9].*$", "\\1", sub("^(NO|LA|GA|GSC)$", "\\1", arm))
}
e8Compare <- function(rows) {
    a <- rows[rows[["tier"]] %in% "D" & rows[["design"]] %in% "D4" & rows[["version"]] %in% "V1g" & rows[["arm"]] %in% c("NO", "LM2", "GK2", "LA", "GA"), , drop = FALSE]
    b <- rows[rows[["tier"]] %in% "ID", , drop = FALSE]
    skip <- c("tier", "cell", "omic", "seconds", "warnings", "smoke", "class", "step")
    key <- function(x) paste(x[["r"]], x[["version"]], x[["arm"]])
    common <- intersect(key(a), key(b))
    worst <- 0
    same <- 0L
    for (k in common) {
        x <- a[key(a) == k, , drop = FALSE][1L, ]
        y <- b[key(b) == k, , drop = FALSE][1L, ]
        deviation <- 0
        for (column in setdiff(names(x), skip)) {
            if (is.numeric(x[[column]])) {
                d <- if (is.na(x[[column]]) && is.na(y[[column]])) 0 else if (is.na(x[[column]]) || is.na(y[[column]])) Inf else abs(x[[column]] - y[[column]])
            } else d <- if (identical(as.character(x[[column]]), as.character(y[[column]]))) 0 else Inf
            deviation <- max(deviation, d)
        }
        worst <- max(worst, deviation)
        same <- same + as.integer(deviation <= 1e-10)
    }
    list(compared = length(common), same = same, worst = worst)
}
if (blind) {
    coverage <- do.call(rbind, lapply(split(rows, paste(vapply(rows[["arm"]], armClassOf, ""), rows[["design"]])), function(x) data.frame(
        item = paste("coverage", armClassOf(x[["arm"]][[1L]]), x[["design"]][[1L]]), value = sprintf("%d journeys, %d completed", nrow(x), sum(x[["valid"]])),
        pass = any(x[["valid"]]), stringsAsFactors = FALSE)))
    e8 <- e8Compare(rows)
    out <- rbind(data.frame(item = "parse", value = sprintf("%d rows, %d columns, %d row files", nrow(rows), ncol(rows), length(row_files)), pass = TRUE,
        stringsAsFactors = FALSE), coverage, data.frame(item = "E8", value = sprintf("%d compared, %d equal, max deviation %s", e8$compared, e8$same,
        format(e8$worst)), pass = e8$compared == 0L || e8$same == e8$compared, stringsAsFactors = FALSE))
    out <- rbind(out, data.frame(item = "instrument_stop", value = sprintf("%d row(s) with status instrument or a guard class", instrument_rows), pass = instrument_rows == 0L, stringsAsFactors = FALSE))
    utils::write.table(out, file.path(out_dir, sprintf("scorer-%s.tsv", job)), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    quit(save = "no", status = if (instrument_rows > 0L) 3L else 0L)
}
if (any(isTrue(rows[["smoke"]]))) stop("K275 scoring stops: the rows carry the smoke flag; a smoke row set is refused for an outcome line")

# ---- scoring -----------------------------------------------------------------------------------------------------------------------
cells <- list()
emit <- function(family, id, tier, design, version, arm, estimand, n_valid, R, estimate, lo, hi, level, sided, label = NA_character_, jk_lo = NA_real_, jk_hi = NA_real_,
    label_jk = NA_character_, note = "") {
    cells[[length(cells) + 1L]] <<- data.frame(family = family, id = id, tier = tier, design = design, version = version, arm = arm, estimand = estimand,
        n_valid = n_valid, R = R, estimate = estimate, lo = lo, hi = hi, level = level, sided = sided, label = label, jk_lo = jk_lo, jk_hi = jk_hi,
        label_jk = label_jk, note = note, stringsAsFactors = FALSE)
}
cpUpper <- function(x, n, alpha) if (n <= 0) NA_real_ else if (x >= n) 1 else stats::qbeta(1 - alpha, x + 1, n - x)
cpLower <- function(x, n, alpha) if (n <= 0) NA_real_ else if (x <= 0) 0 else stats::qbeta(alpha, x, n - x + 1)
bootIndex <- function(design, tier, R) {
    set.seed(277000000 + 100 * CIDX[[design]] + TIER_T[[tier]], kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    matrix(sample.int(R, R * BOOT_B, replace = TRUE), nrow = R)
}
idx_cache <- list()
indexOf <- function(design, tier, R) {
    key <- paste(design, tier, R)
    if (is.null(idx_cache[[key]])) idx_cache[[key]] <<- bootIndex(design, tier, R)
    idx_cache[[key]]
}
quant <- function(x, p) stats::quantile(x[is.finite(x)], probs = p, type = 7, names = FALSE)
#' The estimate of a ratio of sums over data sets (a data set with a missing numerator or denominator weighs 0), with the percentile bootstrap
#' and the delete-one jackknife-normal limits. `alpha` is the tail on each side asked for (one-sided: the level; two-sided: level / 2).
ratioStat <- function(numerator, denominator, idx, alpha, transform = identity) {
    ok <- is.finite(numerator) & is.finite(denominator)
    n0 <- ifelse(ok, numerator, 0)
    d0 <- ifelse(ok, denominator, 0)
    est <- sum(n0) / sum(d0)
    boot <- colSums(matrix(n0[idx], nrow(idx))) / colSums(matrix(d0[idx], nrow(idx)))
    jack <- vapply(which(ok), function(i) (sum(n0) - n0[[i]]) / (sum(d0) - d0[[i]]), numeric(1L))
    se <- if (length(jack) > 1L) sqrt((length(jack) - 1) / length(jack) * sum((jack - mean(jack))^2)) else NA_real_
    z <- stats::qnorm(1 - alpha)
    list(est = transform(est), lo = transform(quant(boot, alpha)), hi = transform(quant(boot, 1 - alpha)), jk_lo = transform(est - z * se),
        jk_hi = transform(est + z * se), n = sum(ok))
}
meanStat <- function(x, idx, alpha) ratioStat(x, ifelse(is.finite(x), 1, NA_real_), idx, alpha)

# the per-data-set vector of one column of one (tier, design, version, arm), indexed by data set 1..R (NA where the journey is not valid)
vec <- function(tier, design, version, arm, column, R, need_valid = TRUE) {
    x <- rep(NA_real_, R)
    s <- rows[rows[["tier"]] %in% tier & rows[["design"]] %in% design & rows[["version"]] %in% version & rows[["arm"]] %in% arm, , drop = FALSE]
    if (need_valid) s <- s[s[["valid"]], , drop = FALSE]
    s <- s[is.finite(s[["r"]]) & s[["r"]] >= 1 & s[["r"]] <= R, , drop = FALSE]
    x[s[["r"]]] <- s[[column]]
    x
}
sigVec <- function(tier, design, version, arm, R) {
    x <- rep(NA_character_, R)
    s <- rows[rows[["tier"]] %in% tier & rows[["design"]] %in% design & rows[["version"]] %in% version & rows[["arm"]] %in% arm & rows[["valid"]], , drop = FALSE]
    s <- s[is.finite(s[["r"]]) & s[["r"]] >= 1 & s[["r"]] <= R, , drop = FALSE]
    x[s[["r"]]] <- s[["tested_sig_y"]]
    x
}
REF_ARM <- c(LM1 = "NO", LM2 = "LM1", LM3 = "LM2", `LM2@20` = "LM1@20", LG2 = "LG1", GK1 = "NO", GK2 = "GK1", GK3 = "GK2", LA = "LA", GA = "GA",
    PEL2 = "LM2", PEG2 = "GK2", LX = "LX")   # R-2 (coordinator): the emulation is compared with VD at the SAME k (no factor is spent on a in VE)
refTiers <- function(tier) if (identical(tier, "V")) "V" else unique(c(tier, "D", "K"))
#' the reference (VD at k - 1) of a V1x arm: the column's vector from the first tier that holds it
refVec <- function(tier, design, arm, column, R) {
    if (!arm %in% names(REF_ARM)) return(NULL)
    ref <- REF_ARM[[arm]]
    for (t in refTiers(tier)) {
        s <- rows[rows[["tier"]] %in% t & rows[["design"]] %in% design & rows[["version"]] %in% "VD" & rows[["arm"]] %in% ref & rows[["valid"]], , drop = FALSE]
        if (nrow(s)) return(list(values = vec(t, design, "VD", ref, column, R), sig = sigVec(t, design, "VD", ref, R), tier = t, arm = ref))
    }
    NULL
}

# ---- descriptive and reported estimands, every (tier, design, version, arm) ----------------------------------------------------
combos <- unique(rows[rows[["version"]] %in% c("V0", "V1g", "V1p", "V1s", "V1t", "V2", "VD"), c("tier", "design", "version", "arm")])
for (i in seq_len(nrow(combos))) {
    tier <- combos[["tier"]][[i]]
    design <- combos[["design"]][[i]]
    version <- combos[["version"]][[i]]
    arm <- combos[["arm"]][[i]]
    if (!design %in% DESIGNS || !tier %in% names(TIER_T)) next   # tier ID: E8 only, no bootstrap cells (R-4)
    R <- declaredR(tier)
    idx <- indexOf(design, tier, R)
    sub <- rows[rows[["tier"]] %in% tier & rows[["design"]] %in% design & rows[["version"]] %in% version & rows[["arm"]] %in% arm, , drop = FALSE]
    nv <- sum(sub[["valid"]])
    # proportions with exact two-sided 95 % Clopper-Pearson limits
    for (when in c("pre", "post")) for (t in seq_along(THETAS)) {
        column <- sprintf("fa_%s_%s", when, THETAS[[t]])
        x <- sub[[column]][sub[["valid"]] & is.finite(sub[[column]])]
        if (length(x)) emit("reported", paste0("P_", when), tier, design, version, arm, sprintf("P_%s(%s)", when, c("0.5", "0.7", "0.9", "0.95")[[t]]), length(x), R,
            mean(x), cpLower(sum(x), length(x), 0.025), cpUpper(sum(x), length(x), 0.025), 0.95, "two", note = "exact Clopper-Pearson")
    }
    means <- c("cy_pre_a", "cy_post_a", "cx_pre_a", "cx_post_a", "qa_pre", "qa_post", "pr", "kappa_c", "df_y", "k_used", "k_found", "n_controls",
        "any_call_y", "null_p05_y", "any_call_x", "null_p05_x", "ffx_pre_05", "ffx_post_05", paste0("ff_pre_", THETAS), paste0("ff_post_", THETAS),
        "pe_n_excluded", "pe_excluded_x")
    for (column in means) {
        x <- vec(tier, design, version, arm, column, R)
        if (any(is.finite(x))) {
            st <- meanStat(x, idx, 0.025)
            emit("reported", column, tier, design, version, arm, sprintf("mean %s", column), st$n, R, st$est, st$lo, st$hi, 0.95, "two", jk_lo = st$jk_lo, jk_hi = st$jk_hi)
        }
    }
    # H2 and H3 (tier E): the share of data sets in which a is excluded, and the shares of the emulation's results after an exclusion
    if (identical(tier, "E")) {
        a_ex <- sub[["pe_a_excluded"]][sub[["valid"]] & sub[["pe_a_excluded"]] %in% c("TRUE", "FALSE")]
        if (length(a_ex) && identical(version, "V1g")) emit("reported", "H2", tier, design, version, arm, "share of data sets in which a is excluded", length(a_ex), R, mean(a_ex == "TRUE"),
            cpLower(sum(a_ex == "TRUE"), length(a_ex), 0.025), cpUpper(sum(a_ex == "TRUE"), length(a_ex), 0.025), 0.95, "two", note = "exact Clopper-Pearson")
        if (arm %in% c("PEL2", "PEG2") && nrow(sub)) for (st_name in c("refused_after_exclusion", "dropped", "uncorrected", "no_da")) emit("reported", paste0("H3_", st_name), tier, design, version, arm,
            sprintf("share of journeys with status %s", st_name), nrow(sub), R, mean(sub[["status"]] == st_name), cpLower(sum(sub[["status"]] == st_name), nrow(sub), 0.025),
            cpUpper(sum(sub[["status"]] == st_name), nrow(sub), 0.025), 0.95, "two", note = "exact Clopper-Pearson; denominator all rows of the arm")
    }
    # C_a quantiles
    for (column in c("cy_pre_a", "cy_post_a")) {
        x <- vec(tier, design, version, arm, column, R)
        x <- x[is.finite(x)]
        if (length(x)) emit("reported", paste0("q_", column), tier, design, version, arm, sprintf("quantiles 0.05, 0.5, 0.95 of %s", column), length(x), R, stats::median(x),
            quant(x, 0.05), quant(x, 0.95), NA_real_, "descriptive", note = "estimate = median; lo, hi = the 5th and 95th percentiles")
    }
    # paired against the clean data set (V0) in the same tier, and the references
    if (!version %in% c("V0", "VD")) {
        v0 <- NULL
        for (t in refTiers(tier)) if (is.null(v0)) {
            if (any(rows[["tier"]] %in% t & rows[["design"]] %in% design & rows[["version"]] %in% "V0" & rows[["arm"]] %in% arm & rows[["valid"]])) v0 <- t
        }
        if (!is.null(v0) && !identical(tier, "ID")) {
            for (column in c("cy_post_a", "cy_pre_a")) {
                d <- vec(tier, design, version, arm, column, R) - vec(v0, design, "V0", arm, column, R)
                if (any(is.finite(d))) {
                    st <- meanStat(d, idx, 0.025)
                    emit("reported", paste0("delta_", column), tier, design, version, arm, sprintf("%s: this version minus V0 (paired)", column), st$n, R, st$est, st$lo, st$hi, 0.95, "two")
                }
            }
            for (column in c("any_call_y", "null_p05_y")) {
                d <- vec(tier, design, version, arm, column, R) - vec(v0, design, "V0", arm, column, R)
                if (any(is.finite(d))) {
                    st <- meanStat(d, idx, 0.025)
                    emit("reported", paste0("delta_", column, "_V0"), tier, design, version, arm, sprintf("%s: this version minus V0 (paired)", column), st$n, R, st$est, st$lo, st$hi, 0.95, "two",
                        note = sprintf("V0's own rate %s; %d data sets unpaired", format(mean(vec(v0, design, "V0", arm, column, R), na.rm = TRUE)), sum(xor(is.finite(vec(tier, design, version, arm, column, R)), is.finite(vec(v0, design, "V0", arm, column, R))))))
                }
            }
        }
        ref <- refVec(tier, design, arm, "se2_mod_y", R)
        if (!is.null(ref)) {
            same <- sigVec(tier, design, version, arm, R) == ref$sig
            same[is.na(same)] <- FALSE
            dropped <- sum(!same & is.finite(vec(tier, design, version, arm, "se2_mod_y", R)) & is.finite(ref$values))
            num_mod <- ifelse(same, vec(tier, design, version, arm, "se2_mod_y", R), NA_real_)
            den_mod <- ifelse(same, ref$values, NA_real_)
            st <- ratioStat(num_mod, den_mod, idx, 0.025, sqrt)
            if (st$n > 0) emit("reported", "rho_mod", tier, design, version, arm, sprintf("rho_mod against VD %s (%s)", ref$arm, ref$tier), st$n, R, st$est, st$lo, st$hi, 0.95, "two",
                jk_lo = st$jk_lo, jk_hi = st$jk_hi, note = sprintf("%d pairs dropped for a differing tested set", dropped))
            un_v <- vec(tier, design, version, arm, "se2_unmod_y", R)
            un_r <- refVec(tier, design, arm, "se2_unmod_y", R)$values
            st <- ratioStat(ifelse(same, un_v, NA_real_), ifelse(same, un_r, NA_real_), idx, 0.025, sqrt)
            if (st$n > 0) emit("reported", "rho_unmod", tier, design, version, arm, sprintf("rho_unmod against VD %s (%s)", ref$arm, ref$tier), st$n, R, st$est, st$lo, st$hi, 0.95, "two", jk_lo = st$jk_lo, jk_hi = st$jk_hi,
                note = sprintf("%d pairs dropped for a differing tested set", dropped))
            share <- c(V1g = "g", V1s = "g", V1p = "p", V1t = "t", V2 = "two")[[version]]
            dot_v <- vec(tier, design, version, arm, paste0("dot_", share), R)
            dot_r <- refVec(tier, design, arm, paste0("dot_", share), R)$values
            den <- vec(tier, design, version, arm, paste0("den_", share), R, need_valid = FALSE)
            st <- ratioStat(ifelse(same, dot_v - dot_r, NA_real_), ifelse(same, den, NA_real_), idx, 0.025)
            if (st$n > 0) emit("reported", "pi", tier, design, version, arm, sprintf("pi, the imprint retained, against VD %s (%s)", ref$arm, ref$tier), st$n, R, st$est, st$lo, st$hi, 0.95, "two", jk_lo = st$jk_lo, jk_hi = st$jk_hi,
                label = labelPi(st$lo, st$hi, st$n, R), note = sprintf("%d pairs dropped for a differing tested set", dropped))
            for (column in c("any_call_y", "null_p05_y")) {
                d <- vec(tier, design, version, arm, column, R) - refVec(tier, design, arm, column, R)$values
                if (any(is.finite(d))) {
                    st <- meanStat(d, idx, 0.025)
                    emit("reported", paste0("delta_", column, "_VD"), tier, design, version, arm, sprintf("%s: this version minus VD %s (paired)", column, ref$arm), st$n, R, st$est, st$lo, st$hi, 0.95, "two",
                        note = sprintf("VD's own rate %s; %d data sets unpaired", format(mean(refVec(tier, design, arm, column, R)$values, na.rm = TRUE)), sum(xor(is.finite(vec(tier, design, version, arm, column, R)), is.finite(refVec(tier, design, arm, column, R)$values)))))
                }
            }
        }
    }
}

# ---- the labelled cells (section 7.1) --------------------------------------------------------------------------------------------
lab <- list()   # the labelled cells: the label by name ("A D3 LM2")
emitMissing <- function(family, id, tier, design, version, arm, estimand, R, level, sided) emit(family, id, tier, design, version, arm, estimand, 0L, R, NA_real_, NA_real_,
    NA_real_, level, sided, label = "not determined (coverage)")
pairedRatio <- function(tier, design, arm, column, R, transform = identity, alpha) {
    ref <- refVec(tier, design, arm, column, R)
    if (is.null(ref)) return(NULL)
    same <- sigVec(tier, design, "V1g", arm, R) == ref$sig
    same[is.na(same)] <- FALSE
    st <- ratioStat(ifelse(same, vec(tier, design, "V1g", arm, column, R), NA_real_), ifelse(same, ref$values, NA_real_), indexOf(design, tier, R), alpha, transform)
    c(st, list(ref = ref, dropped = sum(!same & is.finite(vec(tier, design, "V1g", arm, column, R)))))
}
for (design in DESIGNS) {
    R <- declaredR("D")
    idx <- indexOf(design, "D", R)
    # PC1: P_pre(0.5) of V1g (arm NO), lower limit above 0.90, one-sided 0.01/3
    x <- vec("D", design, "V1g", "NO", "fa_pre_05", R)
    x <- x[is.finite(x)]
    lo <- cpLower(sum(x), length(x), 0.01 / 3)
    pass <- covered(length(x), R) && gt(lo, 0.90)
    emit("control", "PC1", "D", design, "V1g", "NO", "P_pre(0.5) lower limit above 0.90", length(x), R, if (length(x)) mean(x) else NA_real_, lo, cpUpper(sum(x), length(x), 0.01 / 3),
        0.01 / 3, "one", label = if (pass) "passes" else "FAILS")
    lab[[paste("PC1", design)]] <- pass
    # NC1: P(fail(0.5)) of a in V0 on GK2, upper limit at most M_abs
    x <- vec("D", design, "V0", "GK2", "fa_post_05", R)
    x <- x[is.finite(x)]
    up <- cpUpper(sum(x), length(x), 0.01 / 3)
    pass <- covered(length(x), R) && le(up, M_ABS)
    emit("control", "NC1", "D", design, "V0", "GK2", "P(fail(0.5)) of the genuine sample after RUV, upper limit at most 0.05", length(x), R, if (length(x)) mean(x) else NA_real_,
        cpLower(sum(x), length(x), 0.01 / 3), up, 0.01 / 3, "one", label = if (pass) "passes" else "FAILS")
    lab[[paste("NC1", design)]] <- pass
    # family A: P_post(0.5), V1g, absorbed / acts / mixed at 0.025/12
    for (arm in c("LM2", "GK2", "LA", "GA")) {
        x <- vec("D", design, "V1g", arm, "fa_post_05", R)
        x <- x[is.finite(x)]
        up <- cpUpper(sum(x), length(x), 0.025 / 12)
        lw <- cpLower(sum(x), length(x), 0.025 / 12)
        label <- labelA(lw, up, length(x), R)
        emit("labelled", "A", "D", design, "V1g", arm, "P_post(0.5)", length(x), R, if (length(x)) mean(x) else NA_real_, lw, up, 0.025 / 12, "one", label = label)
        lab[[paste("A", design, arm)]] <- label
    }
    # family B: rho_mod, V1g against VD at k - 1 (LM2 against LM1, GK2 against GK1), 0.025/6
    for (arm in c("LM2", "GK2")) {
        st <- pairedRatio("D", design, arm, "se2_mod_y", R, sqrt, 0.025 / 6)
        if (is.null(st)) {
            emitMissing("labelled", "B", "D", design, "V1g", arm, "rho_mod against VD", R, 0.025 / 6, "one")
            lab[[paste("B", design, arm)]] <- "not determined (coverage)"
            next
        }
        label <- labelB(st$lo, st$hi, st$n, R)
        emit("labelled", "B", "D", design, "V1g", arm, sprintf("rho_mod against VD %s", st$ref$arm), st$n, R, st$est, st$lo, st$hi, 0.025 / 6, "one", label = label,
            jk_lo = st$jk_lo, jk_hi = st$jk_hi, label_jk = labelB(st$jk_lo, st$jk_hi, st$n, R),
            note = sprintf("%d pairs dropped for a differing tested set", st$dropped))
        lab[[paste("B", design, arm)]] <- label
    }
    # B': rho_unmod of LX against the derived value, two-sided 0.05/3
    derived <- sqrt((2 * PER[[design]] - 3) / (2 * PER[[design]] - 1))
    st <- pairedRatio("D", design, "LX", "se2_unmod_y", R, sqrt, 0.05 / 3 / 2)
    if (is.null(st)) {
        emitMissing("labelled", "B'", "D", design, "V1g", "LX", "rho_unmod of LX against the derived value", R, 0.05 / 3, "two")
        lab[[paste("B'", design)]] <- "not determined (coverage)"
    } else {
        label <- labelBp(st$lo, st$hi, st$n, R, derived, st$est)
        emit("labelled", "B'", "D", design, "V1g", "LX", sprintf("rho_unmod of LX against the derived %s", format(round(derived, 4))), st$n, R, st$est, st$lo, st$hi, 0.05 / 3, "two",
            label = label, jk_lo = st$jk_lo, jk_hi = st$jk_hi, label_jk = labelBp(st$jk_lo, st$jk_hi, st$n, R, derived, st$est),
            note = sprintf("%d pairs dropped for a differing tested set", st$dropped))
        lab[[paste("B'", design)]] <- label
    }
    # C and C': the paired rise of the any call and of the per-feature null rate, V1g against V0, 0.025/12, both readings
    for (arm in c("LM2", "GK2")) for (kind in c("C", "C'")) {
        column <- if (identical(kind, "C")) "any_call_y" else "null_p05_y"
        margin <- if (identical(kind, "C")) M_ANY else M_ALPHA
        d <- vec("D", design, "V1g", arm, column, R) - vec("D", design, "V0", arm, column, R)
        st <- ratioStat(d, ifelse(is.finite(d), 1, NA_real_), idx, 0.025 / 12)
        ref <- refVec("D", design, arm, column, R)
        label <- labelC(st$lo, st$hi, st$n, R, margin)
        emit("labelled", kind, "D", design, "V1g", arm, sprintf("%s: V1g minus V0, paired", column), st$n, R, st$est, st$lo, st$hi, 0.025 / 12, "one",
            label = label, jk_lo = st$jk_lo, jk_hi = st$jk_hi, label_jk = labelC(st$jk_lo, st$jk_hi, st$n, R, margin),
            note = sprintf("any rise shown: %s; V0 rate %s; VD rate %s", shownC(st$lo, st$n, R), format(mean(vec("D", design, "V0", arm, column, R), na.rm = TRUE)),
                if (is.null(ref)) NA else format(mean(ref$values, na.rm = TRUE))))
        lab[[paste(kind, design, arm)]] <- label
    }
    # pi, the imprint retained (V1g; LM2, GK2, LG2, LX), two-sided 95 %
    for (arm in c("LM2", "GK2", "LG2", "LX")) {
        tier <- if (identical(arm, "LG2")) "K" else "D"
        Rt <- declaredR(tier)
        ref <- refVec(tier, design, arm, "dot_g", Rt)
        if (is.null(ref)) {
            emitMissing("labelled", "pi", tier, design, "V1g", arm, "pi, the imprint retained (A7)", Rt, 0.05, "two")
            lab[[paste("pi", design, arm)]] <- "not determined (coverage)"
            next
        }
        same <- sigVec(tier, design, "V1g", arm, Rt) == ref$sig
        same[is.na(same)] <- FALSE
        st <- ratioStat(ifelse(same, vec(tier, design, "V1g", arm, "dot_g", Rt) - ref$values, NA_real_), ifelse(same, vec(tier, design, "V1g", arm, "den_g", Rt, FALSE), NA_real_),
            indexOf(design, tier, Rt), 0.025)
        emit("labelled", "pi", tier, design, "V1g", arm, "pi, the imprint retained (A7)", st$n, Rt, st$est, st$lo, st$hi, 0.05, "two", label = labelPi(st$lo, st$hi, st$n, Rt),
            jk_lo = st$jk_lo, jk_hi = st$jk_hi, label_jk = labelPi(st$jk_lo, st$jk_hi, st$n, Rt), note = sprintf("%d pairs dropped for a differing tested set", sum(!same & is.finite(vec(tier, design, "V1g", arm, "dot_g", Rt)))))
        lab[[paste("pi", design, arm)]] <- labelPi(st$lo, st$hi, st$n, Rt)
    }
}

# ---- G1 (PC2, NC2) and the falsifiers ---------------------------------------------------------------------------------------------
g1_path <- Sys.getenv("K275_G1")
if (!nzchar(g1_path) || !file.exists(g1_path)) stop("K275 scoring stops: K275_G1 (the G1 table) is not set or not found: an outcome line needs F3 judged")
g1 <- readText(g1_path)
if (any(isTrue(g1[["smoke"]]))) stop("K275 scoring stops: the G1 table carries the smoke flag")
g1[["rho"]] <- num(g1[["rho"]])
g1[["r"]] <- num(g1[["r"]])
for (design in DESIGNS) for (kind in c("PC2", "NC2")) {
    s <- g1[g1[["design"]] %in% design & g1[["kind"]] %in% kind, , drop = FALSE]
    s <- s[order(s[["r"]]), , drop = FALSE]
    n <- PER[[design]]
    target <- if (identical(kind, "PC2")) sqrt((2 * n - 3) / (2 * n - 1)) else 1
    R <- nrow(s)
    set.seed(277000000 + 100 * CIDX[[design]] + 10L, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    idx <- matrix(sample.int(max(R, 1L), max(R, 1L) * BOOT_B, replace = TRUE), nrow = max(R, 1L))
    st <- meanStat(s[["rho"]], idx, 0.01 / 3 / 2)
    pass <- R > 0 && ge(st$lo, target - M_MATCH) && le(st$hi, target + M_MATCH)
    emit("control", kind, "G1", design, NA_character_, NA_character_, sprintf("rho_unmod within 0.03 of %s", format(round(target, 4))), R, R, st$est, st$lo, st$hi, 0.01 / 3, "two",
        label = if (pass) "passes" else "FAILS")
    lab[[paste(kind, design)]] <- pass
}
# F4: an identity check fails in a counted journey (E5 from the rows; E1, E2 from the counted jobs' checks)
e5_bad <- sum(isFalse(rows[["e5_y"]]) | isFalse(rows[["e5_x"]]))
checks_files <- list.files(rows_dir, pattern = "^checks-.*\\.tsv$", full.names = TRUE, recursive = TRUE)
checks <- if (length(checks_files)) do.call(rbind, lapply(checks_files, readText))
idcheck_bad <- if (is.null(checks)) 0L else sum(checks[["check"]] %in% c("E1", "E2") & !isTrue(checks[["pass"]]))
# S11 (R-3): E8 (a finding) and the cross-job identity of tier E's pass-1 NO journeys with tier D's NO journeys (a mismatch STOPS, as F4 does)
e8 <- e8Compare(rows)
crossNO <- function() {
    skip <- c("tier", "seconds", "warnings", grep("^pe_", names(rows), value = TRUE), "valid")
    e <- rows[rows[["tier"]] %in% "E" & rows[["arm"]] %in% "NO" & rows[["r"]] <= 100, , drop = FALSE]
    d <- rows[rows[["tier"]] %in% "D" & rows[["arm"]] %in% "NO", , drop = FALSE]
    bad <- 0L
    compared <- 0L
    for (i in seq_len(nrow(e))) {
        m <- d[d[["design"]] == e[["design"]][[i]] & d[["version"]] == e[["version"]][[i]] & d[["r"]] == e[["r"]][[i]], , drop = FALSE]
        if (!nrow(m)) next
        compared <- compared + 1L
        for (column in setdiff(names(rows), skip)) {
            x <- e[[column]][[i]]
            y <- m[[column]][[1L]]
            differs <- if (is.numeric(x)) !((is.na(x) && is.na(y)) || (!is.na(x) && !is.na(y) && abs(x - y) <= 1e-10)) else !identical(as.character(x), as.character(y))
            if (differs) {
                bad <- bad + 1L
                break
            }
        }
    }
    c(compared = compared, bad = bad)
}
cross <- crossNO()
e12_bad <- if (is.null(checks)) 0L else sum(checks[["check"]] %in% "E12" & !isTrue(checks[["pass"]]))   # reported, not a stop
f1 <- !all(unlist(lab[grepl("^PC1", names(lab))]))
f2 <- !all(unlist(lab[grepl("^NC1", names(lab))]))
f3 <- !all(unlist(lab[grepl("^(PC2|NC2)", names(lab))]))
f4 <- e5_bad > 0L || idcheck_bad > 0L || cross[["bad"]] > 0L

# ---- the outcome row (7.3; the first row that holds; fixed-k cells only) ----------------------------------------------------------
cellA <- function(design, arm) lab[[paste("A", design, arm)]]
lm2 <- vapply(DESIGNS, function(d) identical(cellA(d, "LM2"), "absorbed"), logical(1L))
gk2 <- vapply(DESIGNS, function(d) identical(cellA(d, "GK2"), "absorbed"), logical(1L))
acts <- vapply(DESIGNS, function(d) identical(cellA(d, "LM2"), "acts") && identical(cellA(d, "GK2"), "acts"), logical(1L))
outcome <- if (f1) "Stop (F1)" else if (f2) "Stop (F2)" else if (f3) "Stop (F3)" else if (f4) "Stop (F4)" else if (all(lm2) && all(gk2)) "R1" else
    if (all(gk2) && !all(lm2)) "R2" else if (all(lm2) && !all(gk2)) "R3" else if (all(acts)) "R4" else "R5"
decision <- data.frame(item = "outcome", value = outcome, stringsAsFactors = FALSE)
for (name in sort(names(lab))) if (grepl("^(A|B|B'|C|C'|pi) ", name)) decision <- rbind(decision, data.frame(item = name, value = as.character(lab[[name]]), stringsAsFactors = FALSE))
decision <- rbind(decision, data.frame(item = c("F1", "F2", "F3", "F4"), value = as.character(c(f1, f2, f3, f4)), stringsAsFactors = FALSE),
    data.frame(item = c("E5 failing rows", "E1/E2 failing checks", "E12 failing checks", "E8", "crossNO"), value = c(as.character(c(e5_bad, idcheck_bad, e12_bad)),
        sprintf("%d compared, %d equal, max deviation %s", e8$compared, e8$same, format(e8$worst)), sprintf("%d compared, %d differ", cross[["compared"]], cross[["bad"]])), stringsAsFactors = FALSE))
all_cells <- do.call(rbind, cells)
utils::write.table(all_cells, file.path(out_dir, "cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
utils::write.table(decision, file.path(out_dir, "decision.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
writeLines(c(sprintf("K275_SCORE_R = %s", paste(R_SET, collapse = ",")), sprintf("K275_BOOT_B = %d", BOOT_B), sprintf("rows = %d from %d files", nrow(rows), length(row_files))),
    file.path(out_dir, "settings.txt"))
cat(sprintf("K275 scoring: %d rows, %d cells, outcome %s\n", nrow(rows), nrow(all_cells), outcome))
