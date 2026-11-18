# KIT-124 exactness round (decision M16; exact-round-declaration.md): each variant arm beside its base arm from the same
# runner job, from aggregate.R's k124-cells.tsv, and the paired per-run differences (variant minus base) from the
# per-run tables. Reads tables only; the reading per rule follows the declaration ("no worse": in every judged core
# cell, on every rule, the variant's result is not worse than the base arm's; pass < inconclusive < fail).
#
# Usage: Rscript --vanilla exact-beside.R <out dir> <tree>=<scored dir>,<S1 raw dir|->,<S2 raw dir|->,<S3 raw dir|->,<S4 raw dir|-> ...
# Writes to <out dir>: inventory.tsv, variant-beside-base.tsv, paired-differences.tsv, k-distribution.tsv,
# c-cells.tsv (C6 and C8: arms 1, 2 beside the base and variant arms of change (c)), reading.tsv (per change,
# simulator and subset: the declaration's reading over the judged core cells).
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[[1L]]; trees <- args[-1L]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pairs <- data.frame(
    variant = c("25", "26", "27", "28", "29", "30", "31", "32", "33"),
    base    = c("7",  "7",  "1",  "3",  "7",  "16", "7",  "7",  "16"),
    change  = c("a exact SE", "b trend+robust", "b trend+robust (no RUV)", "b trend+robust (legacy)",
                "c E without units", "c E without units (tech-rep design)", "d 200 permutations",
                "e all together", "e all together (tech-rep design)"), stringsAsFactors = FALSE)
lvl <- c(pass = 0, inconclusive = 1, fail = 2)
rd <- function(p) utils::read.delim(p, quote = "", stringsAsFactors = FALSE, colClasses = c(arm = "character"))
readRuns <- function(dir, sim) {
    x <- utils::read.delim(file.path(dir, sprintf("%s-runs.tsv", tolower(sim))), quote = "", colClasses = "character",
        na.strings = "NA")
    for (m in c("run", "k", "n_calls", "null_p05", "null_p01", "fdp", "power")) x[[m]] <- suppressWarnings(as.numeric(x[[m]]))
    x
}
perrun <- list(); inv <- list(); beside <- list(); paired <- list(); kdist <- list(); ccells <- list(); reading <- list()
for (tr in trees) {
    tree <- sub("=.*$", "", tr); parts <- strsplit(sub("^[^=]*=", "", tr), ",")[[1L]]
    cells <- rd(file.path(parts[[1L]], "k124-cells.tsv"))
    raw <- stats::setNames(parts[-1L], paste0("S", 1:4))
    for (sim in names(raw)[raw != "-"]) {
        runs <- readRuns(raw[[sim]], sim)
        if (sim == "S2") runs <- runs[runs$arm != "16", ]  # as aggregate.R: S2 arm 16 repeats arm 7
        # ---- inventory
        g <- split(runs, paste(runs$cell, runs$arm, runs$subset, sep = "\t"))
        inv[[length(inv) + 1L]] <- do.call(rbind, lapply(g, \(x) data.frame(tree = tree, simulator = sim,
            cell = x$cell[[1L]], arm = x$arm[[1L]], subset = x$subset[[1L]], rows = nrow(x),
            distinct_runs = length(unique(x$run)), first_run = min(x$run), last_run = max(x$run),
            duplicated_keys = sum(duplicated(x$run)), ok = sum(x$status == "ok"), refused = sum(x$status == "refused"),
            crashed = sum(x$status == "crashed"), stringsAsFactors = FALSE)))
        # ---- k (E) distribution per cell, arm, subset
        kdist[[length(kdist) + 1L]] <- do.call(rbind, lapply(g, \(x) {
            k <- x$k[x$status != "crashed"]
            tab <- table(k)
            data.frame(tree = tree, simulator = sim, cell = x$cell[[1L]], arm = x$arm[[1L]], subset = x$subset[[1L]],
                runs = length(k), mean_k = mean(k, na.rm = TRUE), se_k = stats::sd(k, na.rm = TRUE) / sqrt(sum(!is.na(k))),
                share_k_above_0 = mean(k > 0, na.rm = TRUE), share_k_above_0_se = sqrt(mean(k > 0, na.rm = TRUE) *
                    (1 - mean(k > 0, na.rm = TRUE)) / sum(!is.na(k))),
                k_table = paste(sprintf("%s:%d", names(tab), as.integer(tab)), collapse = " "), stringsAsFactors = FALSE)
        }))
        # ---- variant beside base, and paired differences
        sc <- cells[cells$simulator == sim, ]
        for (i in seq_len(nrow(pairs))) {
            v <- pairs$variant[[i]]; b <- pairs$base[[i]]
            for (cl in unique(sc$cell[sc$arm == v])) for (sub in unique(sc$subset[sc$arm == v & sc$cell == cl])) {
                rv <- sc[sc$cell == cl & sc$arm == v & sc$subset == sub, ]; rb <- sc[sc$cell == cl & sc$arm == b & sc$subset == sub, ]
                if (nrow(rv) != 1L || nrow(rb) != 1L) { message("no base for ", sim, " ", cl, " ", v, " ", sub); next }
                rules <- paste0("rule", 1:7)
                lv <- lvl[unlist(rv[rules])]; lb <- lvl[unlist(rb[rules])]
                worse <- rules[!is.na(lv) & !is.na(lb) & lv > lb]; better <- rules[!is.na(lv) & !is.na(lb) & lv < lb]
                to_fail <- rules[!is.na(lv) & !is.na(lb) & lv > lb & lv == 2]
                judged <- rv$cell_role == "core" && !rv$void
                row <- function(r, role) data.frame(tree = tree, simulator = sim, cell = cl, subset = sub,
                    change = pairs$change[[i]], pair = paste(v, "vs", b), role = role, arm = r$arm, cell_role = r$cell_role,
                    judged = judged, runs = r$runs, ok = r$ok, refused = r$refused_designed + r$refused_method,
                    crashed = r$crashed, mean_k = r$mean_k, null05 = r$null05, null05_se = r$null05_se, null01 = r$null01,
                    null01_se = r$null01_se, mean_fdp = r$mean_fdp, mean_fdp_se = r$mean_fdp_se,
                    mean_fdp_upper = r$mean_fdp_upper, fdp10_share = r$fdp10_share, fdp10_lower = r$fdp10_lower,
                    fdp10_upper = r$fdp10_upper, rule3_limit = r$rule3_limit, power = r$power, power_se = r$power_se,
                    anycall_share = r$anycall_share, r[rules], verdict = r$verdict,
                    rules_worse_than_base = if (role == "variant") paste(sub("rule", "", worse), collapse = ",") else "",
                    rules_better_than_base = if (role == "variant") paste(sub("rule", "", better), collapse = ",") else "",
                    stringsAsFactors = FALSE, check.names = FALSE)
                beside[[length(beside) + 1L]] <- rbind(row(rv, "variant"), row(rb, "base"))
                reading[[length(reading) + 1L]] <- data.frame(tree = tree, simulator = sim, cell = cl, subset = sub,
                    pair = paste(v, "vs", b), judged = judged, worse = paste(sub("rule", "", worse), collapse = ","),
                    to_fail = paste(sub("rule", "", to_fail), collapse = ","),
                    better = paste(sub("rule", "", better), collapse = ","), stringsAsFactors = FALSE)
                xv <- runs[runs$cell == cl & runs$arm == v & runs$subset == sub, ]
                xb <- runs[runs$cell == cl & runs$arm == b & runs$subset == sub, ]
                m <- match(xv$run, xb$run); stopifnot(!anyNA(m), nrow(xv) == nrow(xb)); xb <- xb[m, ]
                use <- xv$status != "crashed" & xb$status != "crashed"
                dif <- function(a, c) { d <- (a - c)[use]; d <- d[is.finite(d)]
                    c(mean = if (length(d)) mean(d) else NA_real_, se = if (length(d) > 1L) stats::sd(d) / sqrt(length(d)) else NA_real_,
                      n = length(d), differ = sum(d != 0)) }
                dfdp <- dif(xv$fdp, xb$fdp); dpow <- dif(xv$power, xb$power); dn5 <- dif(xv$null_p05, xb$null_p05)
                dn1 <- dif(xv$null_p01, xb$null_p01); d10 <- dif(as.numeric(xv$fdp > 0.10), as.numeric(xb$fdp > 0.10))
                dk <- dif(xv$k, xb$k); dkp <- dif(as.numeric(xv$k > 0), as.numeric(xb$k > 0))
                if (b == "7") perrun[[paste(tree, sim, cl, sub, v, sep = "|")]] <- data.frame(run = xv$run,
                    fdp = ifelse(use, xv$fdp - xb$fdp, NA_real_), power = ifelse(use, xv$power - xb$power, NA_real_))
                paired[[length(paired) + 1L]] <- data.frame(tree = tree, simulator = sim, cell = cl, subset = sub,
                    change = pairs$change[[i]], pair = paste(v, "vs", b), judged = judged, paired_runs = dfdp[["n"]],
                    runs_fdp_differs = dfdp[["differ"]], runs_power_differs = dpow[["differ"]],
                    d_fdp = dfdp[["mean"]], d_fdp_se = dfdp[["se"]],
                    d_power = dpow[["mean"]], d_power_se = dpow[["se"]], d_power_runs = dpow[["n"]],
                    d_null05 = dn5[["mean"]], d_null05_se = dn5[["se"]], d_null01 = dn1[["mean"]], d_null01_se = dn1[["se"]],
                    d_fdp10_share = d10[["mean"]], d_fdp10_share_se = d10[["se"]], d_k = dk[["mean"]], d_k_se = dk[["se"]],
                    runs_k_differs = dk[["differ"]], d_share_k_above_0 = dkp[["mean"]], d_share_k_above_0_se = dkp[["se"]],
                    stringsAsFactors = FALSE)
            }
        }
        # ---- C6 and C8: arms 1, 2 beside the arms of change (c)
        keep <- sc[sc$cell %in% c("C6", "C8") & sc$arm %in% c("1", "2", "7", "16", "29", "30", "31", "32", "33"), ]
        if (nrow(keep)) ccells[[length(ccells) + 1L]] <- data.frame(tree = tree, keep[, c("simulator", "cell", "arm",
            "subset", "cell_role", "runs", "mean_k", "null05", "null05_se", "null01", "null01_se", "mean_fdp",
            "mean_fdp_se", "fdp10_share", "fdp10_lower", "fdp10_upper", "rule3_limit", "power", "power_se",
            "power_gain", "power_gain_se", "anycall_share", paste0("rule", 1:7), "verdict")], stringsAsFactors = FALSE)
    }
}
w <- function(x, name) {
    x <- do.call(rbind, x); rownames(x) <- NULL
    x[] <- lapply(x, \(v) if (is.double(v)) signif(v, 6) else v)
    utils::write.table(x, file.path(out_dir, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    invisible(x)
}
w(inv, "inventory.tsv"); w(beside, "variant-beside-base.tsv"); w(paired, "paired-differences.tsv")
w(kdist, "k-distribution.tsv"); w(ccells, "c-cells.tsv")
# ---- (e) against its parts: per run, (32 - 7) minus the sum of (part - 7) over the parts run in that cell
add <- list()
for (key in grep("\\|32$", names(perrun), value = TRUE)) {
    stem <- sub("\\|32$", "", key); parts_here <- intersect(c("25", "26", "29", "31"), sub("^.*\\|", "",
        names(perrun)[startsWith(names(perrun), paste0(stem, "|"))]))
    tot <- perrun[[key]]
    for (m in c("fdp", "power")) {
        pred <- Reduce(`+`, lapply(parts_here, \(p) perrun[[paste(stem, p, sep = "|")]][[m]]))
        inter <- tot[[m]] - pred; ok <- is.finite(inter)
        f <- strsplit(stem, "|", fixed = TRUE)[[1L]]
        add[[length(add) + 1L]] <- data.frame(tree = f[[1L]], simulator = f[[2L]], cell = f[[3L]], subset = f[[4L]],
            measure = m, parts = paste(parts_here, collapse = "+"), observed_32_minus_7 = mean(tot[[m]][ok]),
            sum_of_parts = mean(pred[ok]), interaction = mean(inter[ok]),
            interaction_se = stats::sd(inter[ok]) / sqrt(sum(ok)), runs = sum(ok), stringsAsFactors = FALSE)
    }
}
w(add, "all-together-vs-parts.tsv")
rdg <- do.call(rbind, reading)
summ <- do.call(rbind, lapply(split(rdg, paste(rdg$tree, rdg$simulator, rdg$pair, rdg$subset, sep = "|")), \(x) {
    j <- x[x$judged, ]
    lst <- function(col) paste(sprintf("%s(%s)", j$cell[nzchar(j[[col]])], j[[col]][nzchar(j[[col]])]), collapse = " ")
    data.frame(tree = x$tree[[1L]], simulator = x$simulator[[1L]], pair = x$pair[[1L]], subset = x$subset[[1L]],
        judged_core_cells = paste(j$cell, collapse = ","), reported_cells = paste(x$cell[!x$judged], collapse = ","),
        worse_in_judged = lst("worse"), to_fail_in_judged = lst("to_fail"), better_in_judged = lst("better"),
        reading = if (!nrow(j)) "no judged core cell (reported only)" else if (any(nzchar(j$to_fail))) "rejected (a rule moves to fail)" else
            if (any(nzchar(j$worse))) "between: a rule moves from pass to inconclusive (user's decision)" else "no worse in every judged core cell",
        worse_in_reported = paste(sprintf("%s(%s)", x$cell[!x$judged & nzchar(x$worse)], x$worse[!x$judged & nzchar(x$worse)]), collapse = " "),
        stringsAsFactors = FALSE)
}))
rownames(summ) <- NULL
utils::write.table(summ[order(summ$tree, summ$simulator, summ$pair, summ$subset), ], file.path(out_dir, "reading.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
message("written: ", out_dir)
