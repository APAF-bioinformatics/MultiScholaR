# KIT-151 c: assemble one task's job results and apply the declared readings (base R; run once, in the score job).
# Usage: Rscript --vanilla k151c-score.R <task: seed|pin|ruv> <jobs dir> <committed kit-151 dir> <out dir>
# Declarations: docs/qa/workflow-execution/kit/kit-151c/{seed,pin,ruv}-declaration.md (thresholds there and here).

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 4L)
task <- arguments[[1L]]
jobs_dir <- arguments[[2L]]
history_dir <- arguments[[3L]]
out_dir <- arguments[[4L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
md <- character()
say <- function(...) md <<- c(md, sprintf(...))
fmt <- function(x, digits = 3) ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "g"))
write_tsv <- function(x, name) utils::write.table(x, file.path(out_dir, name), sep = "\t", quote = FALSE,
    row.names = FALSE, na = "NA")

greedy_match <- function(r) {
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

# Two fits on the same features: factors matched by |r| of the factor values; sign aligned per matched pair.
compare_fits <- function(a, b) {
    r <- stats::cor(a$Z, b$Z)
    m <- greedy_match(r)
    wa <- stack_w(a$W)
    wb <- stack_w(b$W)
    same_features <- identical(rownames(wa), rownames(wb))
    sz <- vapply(seq_len(nrow(m)), function(i) sign(r[m$a[[i]], m$b[[i]]]), numeric(1))
    dz <- vapply(seq_len(nrow(m)), function(i) max(abs(a$Z[, m$a[[i]]] - sz[[i]] * b$Z[, m$b[[i]]])), numeric(1))
    rw <- if (same_features) vapply(seq_len(nrow(m)), function(i) abs(stats::cor(wa[, m$a[[i]]], wb[, m$b[[i]]])),
        numeric(1)) else rep(NA_real_, nrow(m))
    dw <- if (same_features) vapply(seq_len(nrow(m)), function(i) max(abs(wa[, m$a[[i]]] - sz[[i]] * wb[, m$b[[i]]])),
        numeric(1)) else rep(NA_real_, nrow(m))
    raw_same_shape <- identical(dim(a$Z), dim(b$Z)) && identical(dim(wa), dim(wb))
    list(match = cbind(m, abs_r_W = rw, max_abs_dZ = dz, max_abs_dW = dw),
        min_r_Z = if (nrow(m)) min(m$abs_r) else NA_real_, min_r_W = if (nrow(m)) min(rw) else NA_real_,
        max_dZ = if (nrow(m)) max(dz) else NA_real_, max_dW = if (nrow(m)) max(dw) else NA_real_,
        raw_dZ = if (raw_same_shape) max(abs(a$Z - b$Z)) else NA_real_,
        raw_dW = if (raw_same_shape) max(abs(wa - wb)) else NA_real_,
        identical = identical(unname(a$Z), unname(b$Z)) && identical(unname(wa), unname(wb)) &&
            identical(a$iterations, b$iterations),
        k_a = ncol(a$Z), k_b = ncol(b$Z))
}

read_job_dirs <- function(root) {
    files <- list.files(root, pattern = "^fits\\.rds$", recursive = TRUE, full.names = TRUE)
    lapply(files, function(f) c(readRDS(f), list(dir = dirname(f))))
}

read_rows <- function(root) {
    files <- list.files(root, pattern = "^rows\\.tsv$", recursive = TRUE, full.names = TRUE)
    tables <- lapply(files, function(f) {
        t <- utils::read.delim(f, stringsAsFactors = FALSE)
        t$job_dir <- sub(paste0("^", normalizePath(root), "/?"), "", normalizePath(dirname(f)))
        t
    })
    columns <- unique(unlist(lapply(tables, names)))
    do.call(rbind, lapply(tables, function(t) {
        for (col in setdiff(columns, names(t))) t[[col]] <- NA
        t[columns]
    }))
}

# ================================================================== seed (SCI-06)
if (task == "seed") {
    jobs <- read_job_dirs(jobs_dir)
    rows <- read_rows(jobs_dir)
    write_tsv(rows, "seed-rows.tsv")
    say("# KIT-151 c, task 1 (SCI-06): training seeds")
    say("")
    say("Jobs: %d; CPU: %s; cells: %d, failed: %d.", length(jobs),
        paste(unique(vapply(jobs, function(j) j$job$cpu_model, character(1))), collapse = ", "), nrow(rows),
        sum(rows$status != "ok"))
    fits <- do.call(c, lapply(jobs, function(j) j$fits))
    keys <- names(fits)
    parts <- do.call(rbind, strsplit(keys, "|", fixed = TRUE))
    cells <- unique(parts[, 1:3, drop = FALSE])
    pairs <- list()
    summary <- list()
    factors <- list()
    for (ci in seq_len(nrow(cells))) {
        sel <- which(parts[, 1] == cells[ci, 1] & parts[, 2] == cells[ci, 2] & parts[, 3] == cells[ci, 3])
        sel <- sel[order(as.integer(parts[sel, 4]))]
        seeds <- as.integer(parts[sel, 4])
        cf <- fits[sel]
        for (i in seq_along(cf)) for (j in seq_along(cf)) if (i < j) {
            cmp <- compare_fits(cf[[i]], cf[[j]])
            init_d <- if (!is.null(cf[[i]]$init_Z) && !is.null(cf[[j]]$init_Z)) max(abs(cf[[i]]$init_Z - cf[[j]]$init_Z))
                else NA_real_
            init_r <- if (!is.null(cf[[i]]$init_Z) && !is.null(cf[[j]]$init_Z)) {
                r_init <- stats::cor(cf[[i]]$init_Z, cf[[j]]$init_Z)
                dimnames(r_init) <- list(paste0("a", seq_len(nrow(r_init))), paste0("b", seq_len(ncol(r_init))))
                min(greedy_match(r_init)$abs_r)
            } else NA_real_
            pairs[[length(pairs) + 1L]] <- data.frame(dataset = cells[ci, 1], arm = cells[ci, 2], mode = cells[ci, 3],
                seed_a = seeds[[i]], seed_b = seeds[[j]], identical = cmp$identical, k_a = cmp$k_a, k_b = cmp$k_b,
                min_abs_r_Z = cmp$min_r_Z, min_abs_r_W = cmp$min_r_W, max_abs_dZ = cmp$max_dZ, max_abs_dW = cmp$max_dW,
                iter_a = cf[[i]]$iterations, iter_b = cf[[j]]$iterations, elbo_diff = cf[[i]]$elbo - cf[[j]]$elbo,
                init_max_abs_d = init_d, init_min_abs_r = init_r, stringsAsFactors = FALSE)
        }
        # per factor of the lowest seed (reference) against every other seed
        ref <- cf[[1L]]
        ve <- rowSums(ref$r2)
        for (f in colnames(ref$Z)) {
            rz <- rw <- numeric()
            for (j in seq_along(cf)[-1L]) {
                m <- compare_fits(ref, cf[[j]])$match
                hit <- m[m$a == f, , drop = FALSE]
                rz <- c(rz, if (nrow(hit)) hit$abs_r else 0)
                rw <- c(rw, if (nrow(hit)) hit$abs_r_W else 0)
            }
            planted_r <- abs(stats::cor(ref$planted, ref$Z[, f]))[, 1L]
            factors[[length(factors) + 1L]] <- data.frame(dataset = cells[ci, 1], arm = cells[ci, 2],
                mode = cells[ci, 3], reference_seed = seeds[[1L]], factor = f, ve_sum_views = unname(ve[[f]]),
                closest_planted = names(which.max(planted_r)), closest_planted_abs_r = max(planted_r),
                median_abs_r_Z = stats::median(rz), min_abs_r_Z = min(rz),
                median_abs_r_W = stats::median(rw), min_abs_r_W = min(rw), stringsAsFactors = FALSE)
        }
        p <- do.call(rbind, pairs[vapply(pairs, function(x) x$dataset == cells[ci, 1] && x$arm == cells[ci, 2] &&
            x$mode == cells[ci, 3], logical(1))])
        its <- vapply(cf, function(x) x$iterations, numeric(1))
        elbos <- vapply(cf, function(x) x$elbo, numeric(1))
        tier <- if (all(p$identical)) "A identical" else if (max(p$max_abs_dZ) <= 1e-8) "B machine precision" else
            if (min(p$min_abs_r_Z) >= 0.99) "C same solution" else "D different solutions"
        cr <- rows[rows$dataset == cells[ci, 1] & rows$arm == cells[ci, 2] & rows$mode == cells[ci, 3] &
            rows$status == "ok", ]
        summary[[ci]] <- data.frame(dataset = cells[ci, 1], arm = cells[ci, 2], mode = cells[ci, 3], seeds = length(cf),
            tier = tier, pairs = nrow(p), pairs_identical = sum(p$identical), pairs_min_r_below_0.99 = sum(p$min_abs_r_Z < 0.99),
            min_pair_min_abs_r_Z = min(p$min_abs_r_Z), median_pair_min_abs_r_Z = stats::median(p$min_abs_r_Z),
            min_pair_min_abs_r_W = min(p$min_abs_r_W), max_abs_dZ = max(p$max_abs_dZ), max_abs_dW = max(p$max_abs_dW),
            iterations = paste(sort(unique(its)), collapse = ","), elbo_range = max(elbos) - min(elbos),
            elbo_range_rel = (max(elbos) - min(elbos)) / abs(stats::median(elbos)),
            learned_factors = paste(sort(unique(vapply(cf, function(x) ncol(x$Z), numeric(1)))), collapse = ","),
            init_max_abs_d = max(p$init_max_abs_d), init_min_abs_r = min(p$init_min_abs_r),
            pca_solver = paste(unique(cr$pca_solver), collapse = " / "),
            F1_range = sprintf("%.4f-%.4f", min(cr$planted_F1), max(cr$planted_F1)),
            F2_range = sprintf("%.4f-%.4f", min(cr$planted_F2), max(cr$planted_F2)),
            F3_range = sprintf("%.4f-%.4f", min(cr$planted_F3), max(cr$planted_F3)),
            F2_seeds_ge_0.9 = sum(cr$planted_F2 >= 0.9), F3_seeds_ge_0.9 = sum(cr$planted_F3 >= 0.9),
            stringsAsFactors = FALSE)
    }
    pairs <- do.call(rbind, pairs)
    summary <- do.call(rbind, summary)
    factors <- do.call(rbind, factors)
    write_tsv(pairs, "seed-pairs.tsv")
    write_tsv(summary, "seed-summary.tsv")
    write_tsv(factors, "seed-factors.tsv")
    say("")
    say("| data set | arm | mode | tier | identical pairs | min pair-min abs r (Z) | min (W) | max abs dZ | iterations | ELBO range (rel) | start: max abs d | PCA solver |")
    say("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for (i in seq_len(nrow(summary))) with(summary[i, ], say("| %s | %s | %s | %s | %d/%d | %s | %s | %s | %s | %s (%s) | %s | %s |",
        dataset, arm, mode, tier, pairs_identical, pairs, fmt(min_pair_min_abs_r_Z, 6), fmt(min_pair_min_abs_r_W, 6),
        fmt(max_abs_dZ), iterations, fmt(elbo_range), fmt(elbo_range_rel), fmt(init_max_abs_d), pca_solver))
    say("")
    say("Tiers (declared): A every pair bit-identical; B max abs factor difference after matching <= 1e-8; C every matched factor |r| >= 0.99; D otherwise.")
    tiers <- unique(substr(summary$tier, 1L, 1L))
    say("")
    say("Reading: %s", if (all(tiers %in% c("A", "B"))) "seeds give the same model to machine precision in every cell: multi-seed training cannot measure stability here; resampling would be needed." else
        sprintf("seeds differ in %d of %d cells (tiers %s); magnitudes in the table and seed-factors.tsv.", sum(!substr(summary$tier, 1L, 1L) %in% c("A", "B")), nrow(summary), paste(sort(tiers), collapse = ", ")))
}

# ================================================================== pin (SCI-05)
if (task == "pin") {
    jobs <- read_job_dirs(jobs_dir)
    rows <- read_rows(jobs_dir)
    job_table <- do.call(rbind, lapply(jobs, function(j) {
        blas <- j$job$python$blas
        ob <- Filter(function(x) identical(x$internal_api, "openblas"), blas)
        threads_one <- length(ob) > 0L && all(vapply(ob, function(x) identical(as.integer(x$num_threads), 1L), logical(1)))
        haswell <- length(ob) > 0L && all(vapply(ob, function(x) identical(as.character(x$architecture), "Haswell"),
            logical(1)))
        data.frame(job = basename(dirname(j$dir)), condition = if (is.null(j$job$condition)) NA else j$job$condition,
            cpu_model = j$job$cpu_model, avx512f = isTRUE(j$job$avx512f), status = j$job$status,
            failed_cells = j$job$failed_cells, openblas_threads_one = threads_one, openblas_all_haswell = haswell,
            openblas = paste(vapply(ob, function(x) sprintf("%s:%s:%s:threads=%s", basename(as.character(x$filepath)),
                as.character(x$version), as.character(x$architecture), as.character(x$num_threads)), character(1)),
                collapse = " "),
            other_blas = paste(vapply(Filter(function(x) !identical(x$internal_api, "openblas"), blas),
                function(x) as.character(x$internal_api), character(1)), collapse = " "),
            # numpy's AVX-512 dispatch targets (AVX512F, AVX512CD, AVX512_<name>); sub-features such as AVX512BW are
            # not targets and cannot be switched off
            numpy_avx512_on = paste(grep("^AVX512(F|CD|_[A-Z]+)$", j$job$python$numpy_cpu_features_on, value = TRUE),
                collapse = " "),
            env = paste(names(unlist(j$job$python$environment)), unlist(j$job$python$environment), sep = "=", collapse = " "),
            stringsAsFactors = FALSE)
    }))
    # instrumentation checks (declared): every cell ok; P and PN: every OpenBLAS at 1 thread and Haswell; PN: no
    # numpy AVX512 feature on.
    job_table$check <- with(job_table, ifelse(status != 0 | failed_cells != 0, "cells failed",
        ifelse(condition %in% c("P", "PN") & !(openblas_threads_one & openblas_all_haswell), "pin not in force",
        ifelse(condition == "PN" & nzchar(numpy_avx512_on), "numpy AVX512 still on", "ok"))))
    write_tsv(job_table, "pin-jobs.tsv")
    write_tsv(rows, "pin-rows.tsv")
    say("# KIT-151 c, task 2 (SCI-05): pinning the maths library")
    say("")
    models <- sort(unique(job_table$cpu_model))
    say("Jobs x conditions: %d; CPU models seen: %s.", nrow(job_table), paste(sprintf("%s (%d jobs%s)", models,
        vapply(models, function(m) length(unique(job_table$job[job_table$cpu_model == m])), integer(1)),
        ifelse(vapply(models, function(m) any(job_table$avx512f[job_table$cpu_model == m]), logical(1)), ", AVX-512", "")),
        collapse = "; "))
    say("Instrumentation: %s.", if (all(job_table$check == "ok")) "every job and condition passed its checks" else
        paste(sprintf("%s %s: %s", job_table$job, job_table$condition, job_table$check)[job_table$check != "ok"], collapse = "; "))
    entries <- list()
    for (j in jobs) for (key in names(j$fits)) entries[[length(entries) + 1L]] <- list(job = basename(dirname(j$dir)),
        condition = j$job$condition, cpu = j$job$cpu_model, arm = sub("^[^|]*\\|", "", key), fit = j$fits[[key]])
    classes <- list()
    pairs <- list()
    for (cond in c("U", "P", "PN")) for (arm in unique(vapply(entries, function(e) e$arm, character(1)))) {
        es <- Filter(function(e) identical(e$condition, cond) && e$arm == arm, entries)
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
        for (cl in seq_along(reps)) classes[[length(classes) + 1L]] <- data.frame(condition = cond, arm = arm, class = cl,
            jobs = sum(class_of == cl), cpu_models = paste(sort(unique(vapply(es[class_of == cl], function(e) e$cpu,
            character(1)))), collapse = " + "), iterations = reps[[cl]]$iterations, elbo = reps[[cl]]$elbo,
            stringsAsFactors = FALSE)
        for (i in seq_along(es)) for (k in seq_along(es)) if (i < k) {
            cmp <- compare_fits(es[[i]]$fit, es[[k]]$fit)
            pairs[[length(pairs) + 1L]] <- data.frame(condition = cond, arm = arm, job_a = es[[i]]$job, job_b = es[[k]]$job,
                cpu_a = es[[i]]$cpu, cpu_b = es[[k]]$cpu, same_cpu = es[[i]]$cpu == es[[k]]$cpu, identical = cmp$identical,
                iter_a = es[[i]]$fit$iterations, iter_b = es[[k]]$fit$iterations, raw_max_abs_dZ = cmp$raw_dZ,
                raw_max_abs_dW = cmp$raw_dW, matched_min_abs_r_Z = cmp$min_r_Z, matched_max_abs_dZ = cmp$max_dZ,
                matched_max_abs_dW = cmp$max_dW, stringsAsFactors = FALSE)
        }
    }
    classes <- do.call(rbind, classes)
    pairs <- if (length(pairs)) do.call(rbind, pairs) else data.frame()
    write_tsv(classes, "pin-classes.tsv")
    write_tsv(pairs, "pin-pairs.tsv")
    say("")
    say("| condition | arms | arms with one class | most classes in an arm | same-CPU pairs not identical | cross-CPU pairs not identical | cross-CPU iterations differ | largest cross-CPU abs dZ (raw) | lowest cross-CPU matched abs r |")
    say("|---|---|---|---|---|---|---|---|---|")
    readings <- character()
    for (cond in c("U", "P", "PN")) {
        cc <- classes[classes$condition == cond, ]
        if (!nrow(cc)) next
        per_arm <- table(cc$arm)
        pp <- pairs[pairs$condition == cond, ]
        cross <- pp[!pp$same_cpu, ]
        same <- pp[pp$same_cpu, ]
        say("| %s | %d | %d | %d | %d/%d | %d/%d | %d | %s | %s |", cond, length(per_arm), sum(per_arm == 1L), max(per_arm),
            sum(!same$identical), nrow(same), sum(!cross$identical), nrow(cross), sum(cross$iter_a != cross$iter_b),
            if (nrow(cross)) fmt(max(cross$raw_max_abs_dZ, na.rm = TRUE)) else "NA",
            if (nrow(cross)) fmt(min(cross$matched_min_abs_r_Z, na.rm = TRUE), 6) else "NA")
        readings[[cond]] <- if (!nrow(cross)) "no cross-CPU pair: cannot tell" else if (all(cross$identical))
            "identical across the CPU models seen" else if (all(cross$iter_a == cross$iter_b) &&
            max(cross$raw_max_abs_dZ, na.rm = TRUE) <= 1e-10) "agree to rounding, not bit-identical" else
            "processors still give different models"
    }
    say("")
    for (cond in names(readings)) say("- %s: %s", cond, readings[[cond]])
    # the committed history for this case (first run on the Intel Xeon 6973P-C, follow-up on the AMD EPYC 7763)
    history <- list()
    first <- file.path(history_dir, "results-rows.tsv")
    follow <- file.path(history_dir, "followup-rows.tsv")
    if (file.exists(first)) {
        h <- utils::read.delim(first, stringsAsFactors = FALSE)
        h <- h[h$scenario == "main" & h$seed_index == 7 & h$planted == "F1", c("arm", "iterations", "abs_r")]
        if (nrow(h)) history[[1L]] <- data.frame(source = "first run k151-main-7", cpu = "Intel(R) Xeon(R) 6973P-C", h)
    }
    if (file.exists(follow)) {
        h <- utils::read.delim(follow, stringsAsFactors = FALSE)
        h <- h[h$scenario == "n20weak" & h$seed_index == 7 & h$mode == "fast" & h$planted == "F1", c("arm", "iterations", "abs_r")]
        if (nrow(h)) history[[2L]] <- data.frame(source = "follow-up k151f-n20weak-7", cpu = "AMD EPYC 7763 64-Core Processor", h)
    }
    history <- do.call(rbind, history)
    if (!is.null(history)) {
        u <- rows[rows$condition == "U" & rows$status == "ok", ]
        ours <- do.call(rbind, lapply(seq_len(nrow(history)), function(i) {
            cpu_jobs <- job_table$job[job_table$cpu_model == history$cpu[[i]] & job_table$condition == "U"]
            m <- u[u$arm == history$arm[[i]] & sub("/.*$", "", u$job_dir) %in% cpu_jobs, ]
            data.frame(this_run_jobs_same_cpu = nrow(m), this_run_iterations = paste(unique(m$iterations), collapse = ","),
                this_run_F1 = paste(unique(round(m$planted_F1, 6)), collapse = ","))
        }))
        history <- cbind(history, ours)
        write_tsv(history, "pin-history.tsv")
        say("")
        say("History check (instrumentation, not a gate): unpinned runs on the CPU model of a committed run; see pin-history.tsv. Rows with a same-CPU job: %d; iterations equal to history: %d.",
            sum(history$this_run_jobs_same_cpu > 0), sum(history$this_run_jobs_same_cpu > 0 &
            history$this_run_iterations == as.character(history$iterations)))
    }
}

# ================================================================== ruv (SCI-04)
if (task == "ruv") {
    rows <- read_rows(jobs_dir)
    write_tsv(rows, "ruv-rows.tsv")
    checks <- list.files(jobs_dir, pattern = "^ruv-check\\.tsv$", recursive = TRUE, full.names = TRUE)
    checks <- do.call(rbind, lapply(checks, utils::read.delim))
    if (!is.null(checks)) write_tsv(checks, "ruv-check.tsv")
    jobs <- read_job_dirs(jobs_dir)
    job_table <- do.call(rbind, lapply(jobs, function(j) data.frame(job = basename(dirname(j$dir)), cpu_model = j$job$cpu_model,
        status = j$job$status, failed_cells = j$job$failed_cells, RUVIIIC = as.character(j$job$versions$RUVIIIC %||% NA),
        ruv = as.character(j$job$versions$ruv %||% NA), seconds = j$job$seconds_total, stringsAsFactors = FALSE)))
    write_tsv(job_table, "ruv-jobs.tsv")
    ok <- rows[rows$status == "ok", ]
    ok$group_ve <- (ok$group_factor_ve_A + ok$group_factor_ve_B) / 2
    ok$neglog10_p <- -log10(ok$group_factor_p)
    ok$retained_S <- (ok$retained_S_A + ok$retained_S_B) / 2
    ok$retained_G <- (ok$retained_G_A + ok$retained_G_B) / 2
    ok$retained_U <- (ok$retained_U_A + ok$retained_U_B) / 2
    say("# KIT-151 c, task 3 (SCI-04, INT-004): group-corrected inputs")
    say("")
    say("Jobs: %d (CPU: %s); cells: %d, failed: %d. RUVIII_C_Varying against ruv::RUVIII, largest abs difference: %s.",
        nrow(job_table), paste(unique(job_table$cpu_model), collapse = ", "), nrow(rows), sum(rows$status != "ok"),
        if (is.null(checks)) "not run" else fmt(max(checks$max_abs_diff_vs_ruv_RUVIII)))
    arms <- c("raw", "oracle", "ruv_k1", "ruv_k2", "ruv_k3")
    summary <- list()
    paired <- list()
    for (N in sort(unique(ok$N))) {
        for (arm in intersect(arms, unique(ok$arm))) {
            a <- ok[ok$N == N & ok$arm == arm, ]
            summary[[length(summary) + 1L]] <- data.frame(N = N, arm = arm, seeds = nrow(a), S_recovered = sum(a$r_S >= 0.9),
                median_r_S = stats::median(a$r_S), median_best_r_S = stats::median(a$best_r_S),
                G_recovered = sum(a$r_G >= 0.9), median_r_G = stats::median(a$r_G),
                median_group_ve = stats::median(a$group_ve), median_group_R2 = stats::median(a$group_factor_R2),
                median_neglog10_p = stats::median(a$neglog10_p), median_s_factor_ve = stats::median((a$s_factor_ve_A + a$s_factor_ve_B) / 2, na.rm = TRUE),
                median_total_ve = stats::median((a$total_ve_A + a$total_ve_B) / 2),
                median_retained_S = stats::median(a$retained_S), median_retained_G = stats::median(a$retained_G),
                median_retained_U = stats::median(a$retained_U), median_iterations = stats::median(a$iterations),
                stringsAsFactors = FALSE)
            for (ref in c("raw", "oracle")) if (arm != ref) {
                b <- ok[ok$N == N & ok$arm == ref, ]
                m <- merge(a, b, by = "seed_index", suffixes = c("", ".ref"))
                if (!nrow(m)) next
                paired[[length(paired) + 1L]] <- data.frame(N = N, arm = arm, versus = ref, seeds = nrow(m),
                    median_d_group_ve = stats::median(m$group_ve - m$group_ve.ref),
                    median_ratio_group_ve = stats::median(m$group_ve / m$group_ve.ref),
                    seeds_group_ve_up = sum(m$group_ve > m$group_ve.ref),
                    median_d_R2 = stats::median(m$group_factor_R2 - m$group_factor_R2.ref),
                    seeds_R2_up = sum(m$group_factor_R2 > m$group_factor_R2.ref),
                    median_d_neglog10_p = stats::median(m$neglog10_p - m$neglog10_p.ref),
                    median_d_r_S = stats::median(m$r_S - m$r_S.ref), stringsAsFactors = FALSE)
            }
        }
    }
    summary <- do.call(rbind, summary)
    paired <- do.call(rbind, paired)
    write_tsv(summary, "ruv-summary.tsv")
    write_tsv(paired, "ruv-paired.tsv")
    say("")
    say("| N | arm | S recovered (abs r >= 0.9) | median abs r S | G recovered | group factor: median VE %% (mean of views) | median R2 vs group | median -log10 p | S kept in data | G kept | U kept |")
    say("|---|---|---|---|---|---|---|---|---|---|---|")
    for (i in seq_len(nrow(summary))) with(summary[i, ], say("| %d | %s | %d/%d | %s | %d/%d | %s | %s | %s | %s | %s | %s |", N, arm,
        S_recovered, seeds, fmt(median_r_S), G_recovered, seeds, fmt(median_group_ve), fmt(median_group_R2),
        fmt(median_neglog10_p), fmt(median_retained_S), fmt(median_retained_G), fmt(median_retained_U)))
    say("")
    say("Readings (declared rules):")
    for (N in sort(unique(summary$N))) {
        base <- summary[summary$N == N & summary$arm == "oracle", ]
        if (!nrow(base) || base$S_recovered < 15) {
            say("- N = %d: the oracle arm recovers S in %s of 20 seeds (< 15): the design cannot show removal at this N.", N,
                if (nrow(base)) base$S_recovered else "NA")
            next
        }
        for (k in c("ruv_k1", "ruv_k2", "ruv_k3")) {
            s <- summary[summary$N == N & summary$arm == k, ]
            if (!nrow(s)) next
            removal <- if (s$S_recovered <= 5) "removes the within-group factor" else if (s$S_recovered >= base$S_recovered - 2)
                "keeps the within-group factor" else "partly removes the within-group factor"
            pv <- paired[paired$N == N & paired$arm == k & paired$versus == "oracle", ]
            inflation <- if (nrow(pv) && pv$seeds_group_ve_up >= 15) sprintf("inflates the group factor's variance explained (median x%s, up in %d of %d seeds; R2 median %+.3f)",
                fmt(pv$median_ratio_group_ve), pv$seeds_group_ve_up, pv$seeds, pv$median_d_R2) else
                sprintf("no inflation of the group factor's variance explained shown (up in %s of %s seeds)",
                    if (nrow(pv)) pv$seeds_group_ve_up else "NA", if (nrow(pv)) pv$seeds else "NA")
            say("- N = %d, %s: %s (S recovered %d of %d vs oracle %d); against oracle it %s.", N, k, removal, s$S_recovered,
                s$seeds, base$S_recovered, inflation)
        }
    }
}

writeLines(md, file.path(out_dir, sprintf("%s-score.md", task)))
cat(md, sep = "\n")
