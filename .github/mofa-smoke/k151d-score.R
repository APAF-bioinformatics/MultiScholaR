# KIT-151 d: assemble one task's job results and apply the declared readings (base R; run once, in the score job).
# Usage: Rscript --vanilla k151d-score.R <task: ruv|z|pin> <jobs dir> <reference dir> <out dir>
# Declarations: docs/qa/workflow-execution/kit/kit-151d/{ruv,z,pin}-declaration.md (rules there and here).

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 4L)
task <- arguments[[1L]]
jobs_dir <- arguments[[2L]]
ref_dir <- arguments[[3L]]
out_dir <- arguments[[4L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
md <- character()
say <- function(...) md <<- c(md, sprintf(...))
fmt <- function(x, digits = 3) ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "g"))
f3 <- function(x) ifelse(is.na(x), "NA", sprintf("%+.3f", x))
write_tsv <- function(x, name) utils::write.table(x, file.path(out_dir, name), sep = "\t", quote = FALSE,
    row.names = FALSE, na = "NA")
`%||%` <- function(a, b) if (is.null(a)) b else a

read_rows <- function(root) {
    files <- list.files(root, pattern = "^rows\\.tsv$", recursive = TRUE, full.names = TRUE)
    if (!length(files)) return(NULL)
    tables <- lapply(files, function(f) {
        t <- utils::read.delim(f, stringsAsFactors = FALSE)
        t$job_dir <- sub(paste0("^", normalizePath(root), "/?"), "", normalizePath(dirname(f)))
        job <- jsonlite_free_json(file.path(dirname(f), "job.json"))
        t$cpu_model <- job$cpu_model %||% NA_character_
        t
    })
    columns <- unique(unlist(lapply(tables, names)))
    do.call(rbind, lapply(tables, function(t) {
        for (col in setdiff(columns, names(t))) t[[col]] <- NA
        t[columns]
    }))
}

# job.json without jsonlite (the score job's R has base packages only): the two string fields used here.
jsonlite_free_json <- function(path) {
    if (!file.exists(path)) return(list())
    text <- readLines(path, warn = FALSE)
    get <- function(key) {
        hit <- grep(sprintf('^\\s*"%s": ', key), text, value = TRUE)
        if (!length(hit)) return(NULL)
        value <- sub(sprintf('^\\s*"%s": ', key), "", hit[[1L]])
        value <- sub(",$", "", value)
        gsub('^"|"$', "", value)
    }
    list(cpu_model = get("cpu_model"), cpu_vendor = get("cpu_vendor"), status = get("status"),
        failed_cells = get("failed_cells"), seconds_total = get("seconds_total"), condition = get("condition"),
        avx512f = get("avx512f"))
}

# Paired comparison: mean of per-seed differences with a 95% seed-bootstrap percentile interval (20,000 resamples,
# set.seed(151)); reading with margin m: interval inside [-m, m] "no practical difference"; lower bound > 0 "better"
# ("higher" without margin); upper bound < 0 "worse" ("lower"); otherwise "inconclusive".
boot_ci <- function(d, B = 20000L) {
    set.seed(151L)
    n <- length(d)
    means <- rowMeans(matrix(d[sample.int(n, B * n, replace = TRUE)], nrow = B))
    unname(stats::quantile(means, c(0.025, 0.975)))
}
read_ci <- function(lo, hi, margin) {
    if (is.na(lo)) return("NA")
    if (!is.na(margin) && lo >= -margin && hi <= margin) return("no practical difference")
    if (lo > 0) return(if (is.na(margin)) "higher" else "better")
    if (hi < 0) return(if (is.na(margin)) "lower" else "worse")
    "inconclusive"
}
paired <- function(d, margin) {
    d <- d[is.finite(d)]
    if (length(d) < 3L) return(data.frame(seeds = length(d), mean = NA, lo = NA, hi = NA, median = NA, up = NA,
        down = NA, worst = NA, reading = "NA"))
    ci <- boot_ci(d)
    data.frame(seeds = length(d), mean = mean(d), lo = ci[[1L]], hi = ci[[2L]], median = stats::median(d),
        up = sum(d > 0), down = sum(d < 0), worst = min(d), reading = read_ci(ci[[1L]], ci[[2L]], margin))
}

# ================================================================== ruv (Q1)
if (task == "ruv") {
    rows <- read_rows(jobs_dir)
    write_tsv(rows, "ruv-rows.tsv")
    checks <- list.files(jobs_dir, pattern = "^ruv-check\\.tsv$", recursive = TRUE, full.names = TRUE)
    checks <- do.call(rbind, lapply(checks, utils::read.delim))
    if (!is.null(checks)) write_tsv(checks, "ruv-check.tsv")
    say("# KIT-151 d, Q1: which RUV output feeds integration")
    say("")
    say("Cells: %d, ok %d, failed %d (failed: %s). CPU models: %s. RUVIII_C_Varying against ruv::RUVIII, largest abs difference: %s.",
        nrow(rows), sum(rows$status == "ok"), sum(rows$status != "ok"),
        paste(unique(sprintf("%s %s", rows$arm[rows$status != "ok"], rows$error[rows$status != "ok"])), collapse = "; "),
        paste(sprintf("%s (%d cells)", names(table(rows$cpu_model)), as.integer(table(rows$cpu_model))), collapse = ", "),
        if (is.null(checks)) "not run" else fmt(max(checks$max_abs_diff_vs_ruv_RUVIII)))
    ok <- rows[rows$status == "ok", ]
    ok$group_ve <- (ok$group_factor_ve_A + ok$group_factor_ve_B) / 2
    ok$neglog10_p <- -log10(pmax(ok$group_factor_p, 1e-300))
    for (f in c("S", "G", "U")) ok[[paste0("retained_", f)]] <- (ok[[paste0("retained_", f, "_A")]] + ok[[paste0("retained_", f, "_B")]]) / 2
    say("Fits at the round limit (10,000): %d; most rounds: %s.", sum(ok$stopped_by == "maxiter"), max(ok$iterations))

    # Lifecycle: the strength candidate rule (fixed before the lifecycle run).
    life <- ok[ok$arm == "ideal" & ok$scenario == "main" & ok$seed_index %in% 91:100, ]
    if (nrow(life)) {
        cand <- aggregate(r_S ~ strength + N, life, function(x) c(n = length(x), recovered = sum(x >= 0.9), median = stats::median(x)))
        cand <- data.frame(strength = cand$strength, N = cand$N, cand$r_S)
        write_tsv(cand, "ruv-lifecycle-candidates.tsv")
        say("")
        say("Lifecycle candidates (ideal arm, main scenario, data sets 91 to 100): S recovered (abs r >= 0.9) per strength and N.")
        say("")
        say("| strength | N | data sets | S recovered | median abs r S |")
        say("|---|---|---|---|---|")
        for (i in seq_len(nrow(cand))) say("| %s | %d | %d | %d | %s |", cand$strength[[i]], cand$N[[i]], cand$n[[i]],
            cand$recovered[[i]], fmt(cand$median[[i]]))
        passing <- vapply(c("C1", "C2"), function(s) {
            x <- cand[cand$strength == s, ]
            nrow(x) == 2L && all(x$recovered >= 9)
        }, logical(1))
        say("")
        say("Rule (first candidate, in order C1, C2, with >= 9 of 10 at both N): %s.",
            if (any(passing)) names(passing)[which(passing)[[1L]]] else "no candidate passes")
    }

    counted <- ok[!ok$seed_index %in% 91:100, ]
    if (nrow(counted)) {
        arms <- c("raw", "ideal", "ruv_k1", "ruv_k2", "ruv_k3", "calibrated")
        summary <- list()
        for (N in sort(unique(counted$N))) for (sc in c("main", "conf")) for (arm in arms) {
            a <- counted[counted$N == N & counted$scenario == sc & counted$arm == arm, ]
            if (!nrow(a)) next
            ci <- boot_ci(a$r_S)
            summary[[length(summary) + 1L]] <- data.frame(N = N, scenario = sc, arm = arm, seeds = nrow(a),
                S_recovered = sum(a$r_S >= 0.9), median_r_S = stats::median(a$r_S), mean_r_S = mean(a$r_S),
                mean_r_S_lo = ci[[1L]], mean_r_S_hi = ci[[2L]], G_recovered = sum(a$r_G >= 0.9),
                median_r_G = stats::median(a$r_G), median_group_ve = stats::median(a$group_ve),
                median_group_R2 = stats::median(a$group_factor_R2), median_neglog10_p = stats::median(a$neglog10_p),
                median_retained_S = stats::median(a$retained_S), median_retained_G = stats::median(a$retained_G),
                median_retained_U = stats::median(a$retained_U), median_iterations = stats::median(a$iterations),
                max_iterations = max(a$iterations), at_limit = sum(a$stopped_by == "maxiter"),
                k_A = if (arm == "calibrated") paste(names(table(a$k_A)), table(a$k_A), sep = ":", collapse = " ") else "",
                k_B = if (arm == "calibrated") paste(names(table(a$k_B)), table(a$k_B), sep = ":", collapse = " ") else "",
                k_refused = if (arm == "calibrated") sum(a$k_status_A == "refused" | a$k_status_B == "refused", na.rm = TRUE) else NA,
                stringsAsFactors = FALSE)
        }
        summary <- do.call(rbind, summary)
        write_tsv(summary, "ruv-summary.tsv")
        say("")
        say("## Per arm (median over data sets; S and G recovered = matched abs r >= 0.9)")
        say("")
        say("| N | scenario | arm | S recovered | median abs r S | mean abs r S [95%%] | G recovered | group factor VE %% | R2 vs group | -log10 p | S kept | G kept | U kept | k (calibrated, A; B) |")
        say("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for (i in seq_len(nrow(summary))) with(summary[i, ], say("| %d | %s | %s | %d/%d | %s | %s [%s, %s] | %d/%d | %s | %s | %s | %s | %s | %s | %s |",
            N, scenario, arm, S_recovered, seeds, fmt(median_r_S), fmt(mean_r_S), fmt(mean_r_S_lo), fmt(mean_r_S_hi),
            G_recovered, seeds, fmt(median_group_ve), fmt(median_group_R2), fmt(median_neglog10_p), fmt(median_retained_S),
            fmt(median_retained_G), fmt(median_retained_U), if (nzchar(k_A)) paste0(k_A, "; ", k_B) else ""))

        measures <- list(r_S = 0.02, r_G = 0.02, group_ve = NA, group_factor_R2 = NA, neglog10_p = NA, retained_S = NA,
            retained_G = NA)
        pairs <- list()
        for (N in sort(unique(counted$N))) for (sc in c("main", "conf")) for (ref in c("ideal", "raw")) for (arm in setdiff(arms, ref)) {
            a <- counted[counted$N == N & counted$scenario == sc & counted$arm == ref, ]
            b <- counted[counted$N == N & counted$scenario == sc & counted$arm == arm, ]
            m <- merge(a, b, by = "seed_index", suffixes = c(".a", ".b"))
            if (!nrow(m)) next
            for (meas in names(measures)) {
                p <- paired(m[[paste0(meas, ".b")]] - m[[paste0(meas, ".a")]], measures[[meas]])
                pairs[[length(pairs) + 1L]] <- cbind(data.frame(N = N, scenario = sc, a = ref, b = arm, measure = meas,
                    stringsAsFactors = FALSE), p)
            }
        }
        pairs <- do.call(rbind, pairs)
        write_tsv(pairs, "ruv-pairs.tsv")
        say("")
        say("## Paired against ideal (b minus ideal; mean [95%% seed-bootstrap interval], reading)")
        say("")
        say("| N | scenario | arm | abs r S | abs r G | group VE (points) | R2 | -log10 p | S kept |")
        say("|---|---|---|---|---|---|---|---|---|")
        cellf <- function(p) sprintf("%s [%s, %s] %s", f3(p$mean), f3(p$lo), f3(p$hi), p$reading)
        for (N in sort(unique(pairs$N))) for (sc in c("main", "conf")) for (arm in setdiff(arms, "ideal")) {
            p <- pairs[pairs$N == N & pairs$scenario == sc & pairs$a == "ideal" & pairs$b == arm, ]
            if (!nrow(p)) next
            g <- function(meas) cellf(p[p$measure == meas, ])
            say("| %d | %s | %s | %s | %s | %s | %s | %s | %s |", N, sc, arm, g("r_S"), g("r_G"), g("group_ve"),
                g("group_factor_R2"), g("neglog10_p"), g("retained_S"))
        }
        say("")
        say("## Declared readings")
        for (N in sort(unique(summary$N))) for (sc in c("main", "conf")) {
            base <- summary[summary$N == N & summary$scenario == sc & summary$arm == "ideal", ]
            if (!nrow(base)) next
            if (base$S_recovered < 24) {
                say("- N = %d, %s: the ideal arm recovers S in %d of %d data sets (< 24): the design cannot show removal here; measured effects are reported only.",
                    N, sc, base$S_recovered, base$seeds)
                next
            }
            for (arm in setdiff(arms, "ideal")) {
                s <- summary[summary$N == N & summary$scenario == sc & summary$arm == arm, ]
                if (!nrow(s)) next
                removal <- if (s$S_recovered <= 6) "removes the within-group factor" else if (s$S_recovered >= base$S_recovered - 3)
                    "keeps the within-group factor" else "partly removes the within-group factor"
                say("- N = %d, %s, %s: %s (S recovered %d of %d; ideal %d).", N, sc, arm, removal, s$S_recovered, s$seeds,
                    base$S_recovered)
            }
        }
    }
}

# ================================================================== z (Q2)
if (task == "z") {
    rows <- read_rows(jobs_dir)
    write_tsv(rows, "z-rows.tsv")
    say("# KIT-151 d, Q2: z-score and scale_views")
    say("")
    say("Cells: %d, ok %d, failed %d (failed: %s). CPU models: %s.", nrow(rows), sum(rows$status == "ok"),
        sum(rows$status != "ok"), paste(unique(rows$cell[rows$status != "ok"]), collapse = ", "),
        paste(sprintf("%s (%d cells)", names(table(rows$cpu_model)), as.integer(table(rows$cpu_model))), collapse = ", "))
    ok <- rows[rows$status == "ok", ]
    ok$worst <- pmin(ok$r_F1, ok$r_F2, ok$r_F3)
    say("Fits at the round limit (10,000): %d; most rounds: %s.", sum(ok$stopped_by == "maxiter"), max(ok$iterations))
    counted <- ok[!ok$seed_index %in% 91:92, ]
    summary <- do.call(rbind, lapply(split(ok, list(ok$N, ok$cell, ok$seed_index %in% 91:92), drop = TRUE), function(a) {
        data.frame(N = a$N[[1L]], lifecycle = a$seed_index[[1L]] %in% 91:92, cell = a$cell[[1L]], seeds = nrow(a),
            worst_median = stats::median(a$worst), worst_p10 = unname(stats::quantile(a$worst, 0.1)),
            F1_median = stats::median(a$r_F1), F2_median = stats::median(a$r_F2), F3_median = stats::median(a$r_F3),
            F2_ge_0.9 = sum(a$r_F2 >= 0.9), F3_ge_0.9 = sum(a$r_F3 >= 0.9), iterations_median = stats::median(a$iterations),
            iterations_max = max(a$iterations), stringsAsFactors = FALSE)
    }))
    summary <- summary[order(summary$lifecycle, summary$N, summary$cell), ]
    write_tsv(summary, "z-summary.tsv")
    say("")
    say("| N | cell | data sets | worst-of-three median (10th pct) | F1 | F2 | F3 | F2 >= 0.9 | iterations median (max) |")
    say("|---|---|---|---|---|---|---|---|---|")
    for (i in which(!summary$lifecycle)) with(summary[i, ], say("| %d | %s | %d | %s (%s) | %s | %s | %s | %d | %s (%s) |",
        N, cell, seeds, fmt(worst_median), fmt(worst_p10), fmt(F1_median), fmt(F2_median), fmt(F3_median), F2_ge_0.9,
        iterations_median, iterations_max))

    if (nrow(counted)) {
        wide <- function(N, cell) {
            a <- counted[counted$N == N & counted$cell == cell, c("seed_index", "worst", "r_F1", "r_F2", "r_F3")]
            a
        }
        comparisons <- list(
            c("Z1 z-score all vs none, missing on", "m1_orig_znone_sv1", "m1_orig_zall_sv1"),
            c("Z2 z-score all vs none, missing off", "m0_orig_znone_sv1", "m0_orig_zall_sv1"),
            c("Z4 z-score small data sets vs none, missing on", "m1_orig_znone_sv1", "m1_orig_zsmall_sv1"),
            c("Z5 z-score large data sets vs none, missing on", "m1_orig_znone_sv1", "m1_orig_zlarge_sv1"),
            c("Z6 z-score all vs small only, missing on", "m1_orig_zsmall_sv1", "m1_orig_zall_sv1"),
            c("Z7 z-score small data sets vs none, missing off", "m0_orig_znone_sv1", "m0_orig_zsmall_sv1"),
            c("Z8 z-score large data sets vs none, missing off", "m0_orig_znone_sv1", "m0_orig_zlarge_sv1"),
            c("Z9 z-score all vs none, scale_views off, missing on", "m1_orig_znone_sv0", "m1_orig_zall_sv0"),
            c("Z10 z-score all: scale_views on vs off, missing on", "m1_orig_zall_sv0", "m1_orig_zall_sv1"),
            c("Z11 z-score all vs none, scale_views off, missing off", "m0_orig_znone_sv0", "m0_orig_zall_sv0"))
        for (m in c("m1", "m0")) for (s in c("none", "orig", "lup", "ldown")) {
            comparisons[[length(comparisons) + 1L]] <- c(sprintf("S scale_views on vs off, scales %s, missing %s", s,
                if (m == "m1") "on" else "off"), sprintf("%s_%s_znone_sv0", m, s), sprintf("%s_%s_znone_sv1", m, s))
        }
        interactions <- list(
            c("Z3 z-score gain, missing on minus missing off", "m1_orig_znone_sv1", "m1_orig_zall_sv1", "m0_orig_znone_sv1", "m0_orig_zall_sv1"))
        for (m in c("m1", "m0")) for (s in c("orig", "lup", "ldown")) {
            interactions[[length(interactions) + 1L]] <- c(sprintf("SX scale_views effect at scales %s minus at none, missing %s", s,
                if (m == "m1") "on" else "off"), sprintf("%s_%s_znone_sv0", m, s), sprintf("%s_%s_znone_sv1", m, s),
                sprintf("%s_none_znone_sv0", m), sprintf("%s_none_znone_sv1", m))
        }
        measures <- c("worst", "r_F2", "r_F3", "r_F1")
        out <- list()
        for (N in sort(unique(counted$N))) {
            for (cmp in comparisons) {
                a <- wide(N, cmp[[2L]])
                b <- wide(N, cmp[[3L]])
                if (!nrow(a) || !nrow(b)) next
                m <- merge(a, b, by = "seed_index", suffixes = c(".a", ".b"))
                for (meas in measures) out[[length(out) + 1L]] <- cbind(data.frame(N = N, comparison = cmp[[1L]], a = cmp[[2L]],
                    b = cmp[[3L]], measure = meas, stringsAsFactors = FALSE), paired(m[[paste0(meas, ".b")]] - m[[paste0(meas, ".a")]], 0.02))
            }
            for (cmp in interactions) {
                parts <- lapply(cmp[2:5], function(cell) wide(N, cell))
                if (any(vapply(parts, nrow, integer(1)) == 0L)) next
                m <- Reduce(function(x, y) merge(x, y, by = "seed_index"), lapply(seq_along(parts), function(i) {
                    p <- parts[[i]]
                    names(p)[-1L] <- paste0(names(p)[-1L], ".", i)
                    p
                }))
                for (meas in measures) {
                    d <- (m[[paste0(meas, ".2")]] - m[[paste0(meas, ".1")]]) - (m[[paste0(meas, ".4")]] - m[[paste0(meas, ".3")]])
                    out[[length(out) + 1L]] <- cbind(data.frame(N = N, comparison = cmp[[1L]], a = paste(cmp[4:5], collapse = " -> "),
                        b = paste(cmp[2:3], collapse = " -> "), measure = meas, stringsAsFactors = FALSE), paired(d, 0.02))
                }
            }
        }
        out <- do.call(rbind, out)
        write_tsv(out, "z-pairs.tsv")
        say("")
        say("## Declared comparisons (b minus a; mean [95%% seed-bootstrap interval], reading with margin 0.02)")
        for (N in sort(unique(out$N))) {
            say("")
            say("N = %d", N)
            say("")
            say("| comparison | worst-of-three | F2 | F3 | F1 |")
            say("|---|---|---|---|---|")
            for (cmp in unique(out$comparison[out$N == N])) {
                g <- function(meas) {
                    p <- out[out$N == N & out$comparison == cmp & out$measure == meas, ]
                    sprintf("%s [%s, %s] %s (%d/%d)", f3(p$mean), f3(p$lo), f3(p$hi), p$reading, p$up, p$down)
                }
                say("| %s | %s | %s | %s | %s |", cmp, g("worst"), g("r_F2"), g("r_F3"), g("r_F1"))
            }
        }
        # scale_views on: the scales should not matter (one divisor per data set); differences reported.
        inv <- list()
        for (N in sort(unique(counted$N))) for (m in c("m1", "m0")) {
            cells <- sprintf("%s_%s_znone_sv1", m, c("none", "orig", "lup", "ldown"))
            parts <- lapply(cells, function(cell) wide(N, cell))
            if (any(vapply(parts, nrow, integer(1)) == 0L)) next
            for (i in 2:4) {
                mm <- merge(parts[[1L]], parts[[i]], by = "seed_index")
                d <- abs(c(mm$r_F1.x - mm$r_F1.y, mm$r_F2.x - mm$r_F2.y, mm$r_F3.x - mm$r_F3.y))
                inv[[length(inv) + 1L]] <- data.frame(N = N, missing = m, versus = cells[[i]], seeds = nrow(mm),
                    max_abs_diff = max(d), share_within_1e6 = mean(d <= 1e-6))
            }
        }
        if (length(inv)) {
            inv <- do.call(rbind, inv)
            write_tsv(inv, "z-scale-invariance.tsv")
            say("")
            say("scale_views on, scales none against each other pattern (abs r of F1 to F3; expected near 0): %s",
                paste(sprintf("N %d %s %s max %s (%.0f%% within 1e-6)", inv$N, inv$missing, inv$versus, fmt(inv$max_abs_diff),
                    100 * inv$share_within_1e6), collapse = "; "))
        }
        # tie to k151b (instrumentation): same data, selection and options as k151b's trend1000(_obs50) cells.
        ref_rows <- file.path(ref_dir, "k151b-selection-b-rows.tsv")
        ref_jobs <- file.path(ref_dir, "k151b-selection-b-jobs.tsv")
        if (file.exists(ref_rows)) {
            r <- utils::read.delim(ref_rows, stringsAsFactors = FALSE)
            rj <- utils::read.delim(ref_jobs, stringsAsFactors = FALSE)
            ties <- list()
            for (pair in list(c("m1_orig_znone_sv1", "full", "trend1000_obs50"), c("m0_orig_znone_sv1", "trend", "trend1000"))) {
                ours <- counted[counted$N == 20 & counted$cell == pair[[1L]], ]
                for (i in seq_len(nrow(ours))) for (f in c("F1", "F2", "F3")) {
                    h <- r[r$scenario == pair[[2L]] & r$arm == pair[[3L]] & r$seed_index == ours$seed_index[[i]] & r$planted == f, ]
                    if (!nrow(h)) next
                    ties[[length(ties) + 1L]] <- data.frame(cell = pair[[1L]], seed_index = ours$seed_index[[i]], factor = f,
                        ours = ours[[paste0("r_", f)]][[i]], k151b = h$abs_r[[1L]], iter_ours = ours$iterations[[i]],
                        iter_k151b = h$iterations[[1L]], cpu_ours = ours$cpu_model[[i]],
                        cpu_k151b = rj$cpu_model[rj$seed_index == ours$seed_index[[i]]][1L], stringsAsFactors = FALSE)
                }
            }
            if (length(ties)) {
                ties <- do.call(rbind, ties)
                ties$same_cpu <- ties$cpu_ours == ties$cpu_k151b
                ties$agree_1e6 <- abs(ties$ours - ties$k151b) <= 1e-6
                write_tsv(ties, "z-tie-k151b.tsv")
                say("")
                say("Tie to k151b (instrumentation, not a gate): %d of %d abs r values agree within 1e-6 (same CPU model: %d of %d); iterations equal in %d of %d; largest difference %s.",
                    sum(ties$agree_1e6), nrow(ties), sum(ties$agree_1e6 & ties$same_cpu), sum(ties$same_cpu),
                    sum(ties$iter_ours == ties$iter_k151b), nrow(ties), fmt(max(abs(ties$ours - ties$k151b))))
            }
        }
    }
}

# ================================================================== pin (Q3)
if (task == "pin") {
    greedy_pairs <- function(r) {
        rows <- rownames(r)
        cols <- colnames(r)
        out <- data.frame(a = character(), b = character(), abs_r = numeric(), stringsAsFactors = FALSE)
        r[is.na(r)] <- 0
        while (length(rows) && length(cols)) {
            sub <- abs(r[rows, cols, drop = FALSE])
            at <- which(sub == max(sub), arr.ind = TRUE)[1L, ]
            out[nrow(out) + 1L, ] <- list(rows[[at[[1L]]]], cols[[at[[2L]]]], sub[at[[1L]], at[[2L]]])
            rows <- rows[-at[[1L]]]
            cols <- cols[-at[[2L]]]
        }
        out
    }
    stack_w <- function(W) do.call(rbind, unname(W))
    compare_fits <- function(a, b) {
        r <- stats::cor(a$Z, b$Z)
        m <- greedy_pairs(r)
        wa <- stack_w(a$W)
        wb <- stack_w(b$W)
        sz <- vapply(seq_len(nrow(m)), function(i) sign(r[m$a[[i]], m$b[[i]]]), numeric(1))
        dz <- vapply(seq_len(nrow(m)), function(i) max(abs(a$Z[, m$a[[i]]] - sz[[i]] * b$Z[, m$b[[i]]])), numeric(1))
        same_shape <- identical(dim(a$Z), dim(b$Z)) && identical(dim(wa), dim(wb))
        list(min_r_Z = if (nrow(m)) min(m$abs_r) else NA_real_, max_dZ = if (nrow(m)) max(dz) else NA_real_,
            raw_dZ = if (same_shape) max(abs(a$Z - b$Z)) else NA_real_,
            raw_dW = if (same_shape) max(abs(wa - wb)) else NA_real_,
            identical = identical(unname(a$Z), unname(b$Z)) && identical(unname(wa), unname(wb)) &&
                identical(a$iterations, b$iterations))
    }
    job_dirs <- list.dirs(jobs_dir, recursive = FALSE)
    cpu <- do.call(rbind, lapply(job_dirs, function(d) {
        f <- file.path(d, "cpu.tsv")
        x <- if (file.exists(f)) utils::read.delim(f, header = FALSE, stringsAsFactors = FALSE, quote = "") else
            data.frame(V1 = c("vendor", "model", "full"), V2 = NA)
        v <- stats::setNames(x$V2, x$V1)
        data.frame(job = basename(d), vendor = v[["vendor"]], model = v[["model"]], full = identical(v[["full"]], "true"),
            stringsAsFactors = FALSE)
    }))
    cpu$intel <- cpu$vendor == "GenuineIntel"
    write_tsv(cpu, "pin-cpu.tsv")
    say("# KIT-151 d, Q3: the pin on an Intel processor")
    say("")
    say("Jobs: %d (Intel %d, AMD %d); ran fully: %d. Models: %s.", nrow(cpu), sum(cpu$intel, na.rm = TRUE),
        sum(cpu$vendor == "AuthenticAMD", na.rm = TRUE), sum(cpu$full),
        paste(sprintf("%s (%d)", names(table(cpu$model)), as.integer(table(cpu$model))), collapse = "; "))
    entries <- list()
    jobs <- list()
    for (d in job_dirs[cpu$full]) for (cond in c("U", "P", "PN")) {
        f <- file.path(d, cond, "fits.rds")
        if (!file.exists(f)) next
        x <- readRDS(f)
        jobs[[length(jobs) + 1L]] <- list(job = basename(d), condition = cond, job_record = x$job)
        for (key in names(x$fits)) entries[[length(entries) + 1L]] <- list(job = basename(d), condition = cond,
            cpu = x$job$cpu_model, vendor = x$job$cpu_vendor, case = sub("^[^|]*\\|", "", key), fit = x$fits[[key]])
    }
    job_table <- do.call(rbind, lapply(jobs, function(j) {
        blas <- j$job_record$python$blas
        ob <- Filter(function(x) identical(x$internal_api, "openblas"), blas)
        data.frame(job = j$job, condition = j$condition, cpu_model = j$job_record$cpu_model,
            vendor = j$job_record$cpu_vendor, avx512f = isTRUE(j$job_record$avx512f), status = j$job_record$status,
            failed_cells = j$job_record$failed_cells,
            openblas_threads_one = length(ob) > 0L && all(vapply(ob, function(x) identical(as.integer(x$num_threads), 1L), logical(1))),
            openblas_all_haswell = length(ob) > 0L && all(vapply(ob, function(x) identical(as.character(x$architecture), "Haswell"), logical(1))),
            openblas = paste(vapply(ob, function(x) sprintf("%s:%s:%s:threads=%s", basename(as.character(x$filepath)),
                as.character(x$version), as.character(x$architecture), as.character(x$num_threads)), character(1)), collapse = " "),
            other_blas = paste(vapply(Filter(function(x) !identical(x$internal_api, "openblas"), blas),
                function(x) as.character(x$internal_api), character(1)), collapse = " "),
            numpy_avx512_on = paste(grep("^AVX512(F|CD|_[A-Z]+)$", j$job_record$python$numpy_cpu_features_on, value = TRUE), collapse = " "),
            stringsAsFactors = FALSE)
    }))
    if (!is.null(job_table)) {
        job_table$check <- with(job_table, ifelse(status != 0 | failed_cells != 0, "cells failed",
            ifelse(condition %in% c("P", "PN") & !(openblas_threads_one & openblas_all_haswell), "pin not in force",
            ifelse(condition == "PN" & nzchar(numpy_avx512_on), "numpy AVX512 still on", "ok"))))
        write_tsv(job_table, "pin-jobs.tsv")
        say("Instrumentation: %s.", if (all(job_table$check == "ok")) "every job and condition passed its checks" else
            paste(sprintf("%s %s: %s", job_table$job, job_table$condition, job_table$check)[job_table$check != "ok"], collapse = "; "))
    }
    rows <- read_rows(jobs_dir)
    if (!is.null(rows)) write_tsv(rows, "pin-rows.tsv")
    classes <- list()
    pairs <- list()
    for (cond in c("U", "P", "PN")) for (case in unique(vapply(entries, function(e) e$case, character(1)))) {
        es <- Filter(function(e) identical(e$condition, cond) && e$case == case, entries)
        if (!length(es)) next
        class_of <- integer(length(es))
        reps <- list()
        for (i in seq_along(es)) {
            hit <- which(vapply(reps, function(r) compare_fits(r, es[[i]]$fit)$identical, logical(1)))
            if (length(hit)) class_of[[i]] <- hit[[1L]] else {
                reps[[length(reps) + 1L]] <- es[[i]]$fit
                class_of[[i]] <- length(reps)
            }
        }
        for (cl in seq_along(reps)) classes[[length(classes) + 1L]] <- data.frame(condition = cond, case = case, class = cl,
            jobs = sum(class_of == cl), cpu_models = paste(sort(unique(vapply(es[class_of == cl], function(e) e$cpu,
            character(1)))), collapse = " + "), iterations = reps[[cl]]$iterations, elbo = reps[[cl]]$elbo,
            stringsAsFactors = FALSE)
        for (i in seq_along(es)) for (k in seq_along(es)) if (i < k) {
            cmp <- compare_fits(es[[i]]$fit, es[[k]]$fit)
            pairs[[length(pairs) + 1L]] <- data.frame(condition = cond, case = case, job_a = es[[i]]$job, job_b = es[[k]]$job,
                cpu_a = es[[i]]$cpu, cpu_b = es[[k]]$cpu, vendor_a = es[[i]]$vendor, vendor_b = es[[k]]$vendor,
                same_cpu = es[[i]]$cpu == es[[k]]$cpu, cross_vendor = es[[i]]$vendor != es[[k]]$vendor,
                identical = cmp$identical, iter_a = es[[i]]$fit$iterations, iter_b = es[[k]]$fit$iterations,
                raw_max_abs_dZ = cmp$raw_dZ, raw_max_abs_dW = cmp$raw_dW, matched_min_abs_r_Z = cmp$min_r_Z,
                matched_max_abs_dZ = cmp$max_dZ, stringsAsFactors = FALSE)
        }
    }
    if (length(classes)) {
        classes <- do.call(rbind, classes)
        pairs <- if (length(pairs)) do.call(rbind, pairs) else data.frame()
        write_tsv(classes, "pin-classes.tsv")
        write_tsv(pairs, "pin-pairs.tsv")
        say("")
        say("| condition | cases | cases with one class | Intel-AMD pairs not identical | Intel-AMD pairs with different iterations | largest Intel-AMD abs dZ (raw) | lowest Intel-AMD matched abs r | same-model pairs not identical |")
        say("|---|---|---|---|---|---|---|---|")
        readings <- character()
        for (cond in c("U", "P", "PN")) {
            cc <- classes[classes$condition == cond, ]
            if (!nrow(cc)) next
            per_case <- table(cc$case)
            pp <- if (nrow(pairs)) pairs[pairs$condition == cond, ] else pairs
            cross <- if (nrow(pp)) pp[pp$cross_vendor, ] else pp
            same <- if (nrow(pp)) pp[pp$same_cpu, ] else pp
            say("| %s | %d | %d | %d/%d | %d | %s | %s | %d/%d |", cond, length(per_case), sum(per_case == 1L),
                sum(!cross$identical), nrow(cross), sum(cross$iter_a != cross$iter_b),
                if (nrow(cross)) fmt(max(cross$raw_max_abs_dZ, na.rm = TRUE)) else "NA",
                if (nrow(cross)) fmt(min(cross$matched_min_abs_r_Z, na.rm = TRUE), 6) else "NA",
                sum(!same$identical), nrow(same))
            readings[[cond]] <- if (!nrow(cross)) "no Intel-AMD pair: cannot tell" else if (all(cross$identical))
                "bit-identical across Intel and AMD" else if (all(cross$iter_a == cross$iter_b) &&
                max(cross$raw_max_abs_dZ, na.rm = TRUE) <= 1e-10) "agree to rounding across Intel and AMD, not bit-identical" else
                "Intel and AMD give different models"
        }
        say("")
        for (cond in names(readings)) say("- %s: %s", cond, readings[[cond]])
        # per case under U: Intel against AMD, planted abs r per vendor
        if (!is.null(rows)) {
            u <- rows[rows$condition == "U" & rows$status == "ok", ]
            u$vendor <- cpu$vendor[match(sub("/.*$", "", u$job_dir), cpu$job)]
            per_case <- do.call(rbind, lapply(split(u, list(u$case, u$vendor), drop = TRUE), function(a) data.frame(case = a$case[[1L]],
                vendor = a$vendor[[1L]], jobs = nrow(a), iterations = paste(sort(unique(a$iterations)), collapse = ","),
                F1 = paste(sprintf("%.4f", sort(unique(round(a$planted_F1, 4)))), collapse = ","),
                F2 = paste(sprintf("%.4f", sort(unique(round(a$planted_F2, 4)))), collapse = ","),
                F3 = paste(sprintf("%.4f", sort(unique(round(a$planted_F3, 4)))), collapse = ","), stringsAsFactors = FALSE)))
            if (!is.null(per_case)) {
                write_tsv(per_case, "pin-unpinned-by-vendor.tsv")
                say("")
                say("Unpinned, by vendor (distinct iteration counts and planted abs r): see pin-unpinned-by-vendor.tsv.")
            }
            # history (instrumentation): k151c pin rows (fast case) and k151b rows (hard cases) on the same CPU model.
            hist <- list()
            kc <- file.path(ref_dir, "k151c-pin-rows.tsv")
            kcj <- file.path(ref_dir, "k151c-pin-jobs.tsv")
            if (file.exists(kc)) {
                h <- utils::read.delim(kc, stringsAsFactors = FALSE)
                hj <- utils::read.delim(kcj, stringsAsFactors = FALSE)
                h$cpu <- hj$cpu_model[match(paste(sub("/.*$", "", h$job_dir), h$condition), paste(hj$job, hj$condition))]
                for (cond in c("U", "P")) {
                    ours <- rows[rows$condition == cond & rows$status == "ok" & startsWith(rows$case, "fast:"), ]
                    for (i in seq_len(nrow(ours))) {
                        m <- h[h$condition == cond & h$arm == sub("^fast:", "", ours$case[[i]]) &
                            (cond == "P" | h$cpu == ours$cpu_model[[i]]), ]
                        if (!nrow(m)) next
                        hist[[length(hist) + 1L]] <- data.frame(source = "k151c pin", condition = cond, case = ours$case[[i]],
                            job = ours$job_dir[[i]], cpu = ours$cpu_model[[i]], matches = sum(m$iterations == ours$iterations[[i]] &
                            abs(m$planted_F2 - ours$planted_F2[[i]]) <= 1e-9), of = nrow(m), stringsAsFactors = FALSE)
                    }
                }
            }
            kb <- file.path(ref_dir, "k151b-selection-b-rows.tsv")
            kbj <- file.path(ref_dir, "k151b-selection-b-jobs.tsv")
            if (file.exists(kb)) {
                h <- utils::read.delim(kb, stringsAsFactors = FALSE)
                hj <- utils::read.delim(kbj, stringsAsFactors = FALSE)
                for (spec in list(c("hard_prev10_none", "prev", "none", "10"), c("hard_full23_default", "full", "trend1000_obs50", "23"))) {
                    ref <- h[h$scenario == spec[[2L]] & h$arm == spec[[3L]] & h$seed_index == as.integer(spec[[4L]]) & h$planted == "F2", ]
                    ref_cpu <- hj$cpu_model[hj$seed_index == as.integer(spec[[4L]])][1L]
                    ours <- u[u$case == spec[[1L]], ]
                    for (i in seq_len(nrow(ours))) hist[[length(hist) + 1L]] <- data.frame(source = sprintf("k151b seed %s (%s)", spec[[4L]], ref_cpu),
                        condition = "U", case = spec[[1L]], job = ours$job_dir[[i]], cpu = ours$cpu_model[[i]],
                        matches = as.integer(nrow(ref) == 1L && ref$iterations == ours$iterations[[i]] &&
                            abs(ref$abs_r - ours$planted_F2[[i]]) <= 1e-9), of = nrow(ref), stringsAsFactors = FALSE)
                }
            }
            if (length(hist)) {
                hist <- do.call(rbind, hist)
                write_tsv(hist, "pin-history.tsv")
                say("")
                say("History (instrumentation, not a gate): %s", paste(vapply(split(hist, list(hist$source, hist$condition), drop = TRUE),
                    function(x) sprintf("%s %s: %d of %d rows match a committed row (same CPU model where unpinned)", x$source[[1L]],
                        x$condition[[1L]], sum(x$matches > 0), nrow(x)), character(1)), collapse = "; "))
            }
        }
    }
}

writeLines(md, file.path(out_dir, sprintf("%s-score.md", task)))
cat(md, sep = "\n")
