# KIT-201 M1 (ticket section 5.2, revision 1.3): base against head on simulated data with known truth, per kernel
# (L lipidomics, M metabolomics, P proteomics' legacy runTestsContrasts). Evidence code under docs/qa; no package code.
# Synthetic data only, generated here from fixed seeds.
#
# Usage, from the package root of the head commit on a run-tests runner (full git history):
#   Rscript m1-sim.R r0 <base commit> [<out dir>]
#       The blind lifecycle R0: every kernel and cell, data sets r = 9001 to 9005 (main and calibration sets), both
#       sides, plus the forced-failure unit on head. Writes ONLY r0-checks.tsv, r0-sd.tsv and the formula's output
#       precision.tsv (5.2 "Family and precision"); nothing else of R0 is written or printed.
#   Rscript m1-sim.R m1 <kernel L|M|P> <base commit> [<out dir>] [--cells=1,2,...] [--part=<name>]
#       The counted run for one kernel at the replicate counts of precision.tsv (which must exist in <out dir>).
#       Without --cells: every cell, then the summary. With --cells/--part: those cells only, raw rows to
#       raw-<kernel>-<part>.rds (a split job, 10 "When the formula sets R above 200"); the summary is then made by
#   Rscript m1-sim.R summary <kernel L|M|P> [<out dir>]
#       which reads every raw-<kernel>*.rds in <out dir>.
#   Rscript m1-sim.R family [<out dir>]
#       Writes the gated family (family.tsv) as this script enumerates it; the committed file must equal it.
#   <out dir> default docs/qa/workflow-execution/kit/kit-201/m1 (made absolute).
# Child process (started by the parent):
#   Rscript m1-sim.R --unit <side> <tree> <commit> <kernel> <cell> <set main|calib> <r from> <r to> <rds> <pv dir> <fail r>
#
# Mechanism (as KIT-196's S1):
# - Base is a detached git worktree of <base commit>; head is this tree. Each side's tree is prepared as the run-tests jobs
#   do (src/build_dia_rsum.R when the DIA manifest is absent) and its native code compiled once (pkgbuild::compile_dll),
#   before any child starts. Each unit (one kernel, cell, set and block of data sets) runs in a fresh Rscript child, four at
#   a time (parallel::mclapply in the parent, which loads no package); the child loads its side's tree with
#   pkgload::load_all(tree, export_all = TRUE, helpers = FALSE, compile = FALSE) and asserts that MultiScholaR's namespace
#   is that tree. Base runs first; each head child reads the base's p-values of the same data sets (pv dir) to compute
#   |delta p| and the calls made on one side only, and every row carries its data set's SHA-256: a base/head mismatch is
#   an instrument defect and stops the summary.
# - A data set whose kernel call fails is recorded with its message (status failed) and stops its cell's summary.
#
# Fixed before any run (ticket 5.2, recorded in the family commit):
# - Not tested: a feature with no p-value in a contrast's table (no value in any sample: KIT-196/KIT-200's untested rows,
#   AveExpr NA; or a contrast that cannot be estimated) is not called, is counted per data set (n_untested, n_na_p), and is
#   left out of every denominator (null rates overall and by stratum, power). FDP counts calls only.
# - Called: fdr_qvalue < 0.05 (the shipped significance column). Null rate: raw p < 0.05.
# - Thirds: tertiles of the true level (type-7 quantiles at 1/3 and 2/3 over the data set's features). The missing-value
#   stratum (VD20 cells): features with at least one missing value among the data set's samples.
# - z for V1 and V2: (level - generator mean) / generator SD (2.2 and 3 for L and M; 22 and 2 for P).
# - Seeds: data set r of cell c: 201000000 + 10000 c + r (M1 main r = 1..R; calibration r = 5001..5000 + R; R0
#   r = 9001..9005 for both). L and M share the generator, hence the data; L and M's two calibration thresholds share the
#   seed (the same noise at tau 0.5 and 1). RNGkind Mersenne-Twister, Inversion, Rejection.
# - Rule 2's two differences: head - base (paired) and head - 0.05; a comparison is gated only when its projected SE,
#   1.5 x (the larger of its two SDs) / sqrt(R), is at most SE*; FDP's: 1.5 x SD(head - base) / sqrt(R). A comparison
#   whose SD cannot be computed in R0 (fewer than two data sets with a defined rate) is reported, not gated.
# - Bootstrap: 2,000 resamples of data sets, set.seed(201999) before each comparison, 95 % percentile interval.
# - Cost projection (the 180.0 core-hour cap): worker-seconds = sum over cells and sides of R x (R0 seconds per data set,
#   main plus calibration) x 1.15 (child start-up); core-hours = worker-seconds / 3600 + 10 minutes x 4 cores per job,
#   jobs per kernel = ceiling(kernel's wall minutes at 4 workers / 170). Over the cap: the cell of largest projected cost
#   with R above 200 loses ceiling(5 % of its R), repeatedly; R never goes below 200.

K201_SEED_BASE <- 201000000
K201_BOOT_SEED <- 201999L
K201_BOOT_N <- 2000L
K201_DELTA <- 0.01
K201_ALPHA <- 0.05
K201_WORKERS <- 4L
K201_CHUNK <- 25L
K201_R0 <- 9001:9005
K201_R_MIN <- 200L
K201_R_MAX <- 1000L
K201_CAP <- 180.0
K201_FORCED <- list(kernel = "L", cell = 1L, r = 9003L)
K201_VARIANCES <- c("V0", "V1", "V2")
K201_MISSING <- c("none", "VD20")
K201_DESIGNS <- list(`2x3` = c(groups = 2L, per = 3L), `2x6` = c(groups = 2L, per = 6L), `4x3` = c(groups = 4L, per = 3L))
K201_KERNELS <- list(
    L = list(n = 500L, mean = 2.2, sd = 3, effects = 100L, b05 = 12L, b1 = 12L, taus = c(0.5, 1)),
    M = list(n = 500L, mean = 2.2, sd = 3, effects = 100L, b05 = 12L, b1 = 12L, taus = c(0.5, 1)),
    P = list(n = 2000L, mean = 22, sd = 2, effects = 400L, b05 = 50L, b1 = 50L, taus = c(0, 1)))
K201_STRATA <- c("overall", "third1", "third2", "third3", "missing")
K201_ROW_COLUMNS <- c("side", "commit", "kernel", "cell", "set", "r", "tau", "contrast", "status", "message", "seconds",
    "hash", "n_features", "n_untested", "n_na_p", "n_zero_df", "s2prior_distinct", "calls", "false_calls", "true_calls",
    "fdp", "power", "fdp_third1", "fdp_third2", "fdp_third3", "null_overall", "null_third1", "null_third2", "null_third3",
    "null_missing", "den_overall", "den_third1", "den_third2", "den_third3", "den_missing", "max_dp", "only_base",
    "only_head")

k201_argv <- commandArgs(trailingOnly = TRUE)
k201_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))
k201_rscript <- file.path(R.home("bin"), "Rscript")

k201Flat <- function(x, n = 500L) substr(gsub("[\r\n\t]+", " ", paste(as.character(x), collapse = " ")), 1L, n)

# ---- cells, family ----------------------------------------------------------------------------------------------------
# Cell c (1 to 18): variance model (outer), missingness, design (inner).
k201Cells <- function() {
    grid <- expand.grid(design = names(K201_DESIGNS), missing = K201_MISSING, variance = K201_VARIANCES, stringsAsFactors = FALSE)
    data.frame(cell = seq_len(nrow(grid)), variance = grid$variance, missing = grid$missing, design = grid$design,
        stringsAsFactors = FALSE)
}

k201Groups <- function(design) LETTERS[seq_len(K201_DESIGNS[[design]][["groups"]])]
k201Contrasts <- function(design) paste0("group", k201Groups(design)[-1L], "-groupA")

# The gated comparisons (5.2 "Decision rule"): per kernel, cell, tau and contrast, FDP (main sets) and the null rate
# overall and per stratum (calibration sets at tau > 0, the main sets' beta = 0 features at tau = 0).
k201Family <- function() {
    cells <- k201Cells()
    rows <- list()
    for (kernel in names(K201_KERNELS)) for (i in seq_len(nrow(cells))) {
        cell <- cells[i, ]
        for (tau in K201_KERNELS[[kernel]]$taus) for (contrast in k201Contrasts(cell$design)) {
            strata <- if (identical(cell$missing, "VD20")) K201_STRATA else K201_STRATA[-5L]
            rows[[length(rows) + 1L]] <- data.frame(kernel = kernel, cell = cell$cell, variance = cell$variance,
                missing = cell$missing, design = cell$design, tau = as.character(tau), contrast = contrast,
                rule = c("fdp", rep("null", length(strata))), stratum = c("overall", strata),
                data = c("main", rep(if (tau > 0) "calib" else "main", length(strata))), stringsAsFactors = FALSE)
        }
    }
    do.call(rbind, rows)
}

k201SeStar <- function(K) K201_DELTA / (1.96 + stats::qnorm(0.95^(1 / K)))

# ---- generator --------------------------------------------------------------------------------------------------------
k201SolveLod <- function(values, rate) {
    stats::uniroot(function(lod) mean(stats::plogis(lod - values)) - rate, c(min(values) - 60, max(values) + 60), tol = 1e-10)$root
}

# One data set: values (NA where censored), the design, the true level and, per contrast, the true beta.
k201Generate <- function(kernel, cell, set, r, tau = NULL) {
    spec <- K201_KERNELS[[kernel]]
    info <- k201Cells()[cell, ]
    set.seed(K201_SEED_BASE + 10000 * cell + r, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    n <- spec$n
    groups <- k201Groups(info$design)
    per <- K201_DESIGNS[[info$design]][["per"]]
    group <- rep(groups, each = per)
    runs <- sprintf("S%02d", seq_along(group))
    level <- stats::rnorm(n, spec$mean, spec$sd)
    s2 <- 4 * 0.25^2 / stats::rchisq(n, 4)
    z <- (level - spec$mean) / spec$sd
    if (info$variance %in% c("V1", "V2")) s2 <- s2 * exp(-0.541 * z)
    hyper <- sample.int(n, floor(0.05 * n))
    if (identical(info$variance, "V2")) s2[hyper] <- s2[hyper] * 8
    order <- sample.int(n)
    classes <- list(effect = order[seq_len(spec$effects)], b05 = order[spec$effects + seq_len(spec$b05)],
        b1 = order[spec$effects + spec$b05 + seq_len(spec$b1)])
    size <- numeric(n)
    size[classes$effect] <- stats::runif(spec$effects, 0, 2.5)
    size[classes$b05] <- 0.5
    size[classes$b1] <- 1
    sign <- numeric(n)
    for (members in classes) sign[members] <- sample(rep(c(1, -1), length.out = length(members)))
    in_group <- if (length(groups) == 2L) rep("B", n) else sample(c("B", "C", "D"), n, replace = TRUE)
    beta <- matrix(0, n, length(groups) - 1L, dimnames = list(NULL, groups[-1L]))
    if (is.null(tau)) {
        for (g in groups[-1L]) beta[, g] <- ifelse(in_group == g, sign * size, 0)
    } else {
        # Boundary-calibration set (R1-03): every feature at |beta| = tau in every non-reference group.
        calib_sign <- sample(rep(c(1, -1), length.out = n))
        for (g in groups[-1L]) beta[, g] <- calib_sign * tau
    }
    mean_matrix <- matrix(level, n, length(group))
    for (g in groups[-1L]) mean_matrix[, group == g] <- mean_matrix[, group == g] + beta[, g]
    y <- mean_matrix + matrix(stats::rnorm(n * length(group)), n, length(group)) * sqrt(s2)
    u <- matrix(stats::runif(n * length(group)), n, length(group))
    if (identical(info$missing, "VD20")) y[u < stats::plogis(k201SolveLod(y, 0.20) - y)] <- NA_real_
    ids <- sprintf("F%05d", seq_len(n))
    dimnames(y) <- list(ids, runs)
    colnames(beta) <- paste0("group", groups[-1L], "-groupA")
    rownames(beta) <- ids
    design <- data.frame(Run = runs, group = group, replicates = seq_along(group), stringsAsFactors = FALSE)
    rownames(design) <- runs
    cuts <- stats::quantile(level, c(1 / 3, 2 / 3), type = 7, names = FALSE)
    third <- 1L + (level > cuts[[1L]]) + (level > cuts[[2L]])
    list(y = y, design = design, level = level, third = third, beta = beta, contrasts = colnames(beta),
        with_missing = rowSums(is.na(y)) > 0, vd20 = identical(info$missing, "VD20"), hash = digest::digest(list(y, design, beta), algo = "sha256"))
}

# ---- kernels ----------------------------------------------------------------------------------------------------------
# The kernel's tables per contrast and its public fit, as the entry points of 5.2 call them.
k201Kernel <- function(kernel, data, tau) {
    if (kernel %in% c("L", "M")) {
        fit_fn <- get(if (identical(kernel, "L")) "fitLipidDaAssayModel" else "fitMetabDaAssayModel")
        apply_fn <- get(if (identical(kernel, "L")) "applyLipidDaContrasts" else "applyMetabDaContrasts")
        fit <- fit_fn(data$y, data$design, "~ 0 + group", sample_id_col = "Run")
        result <- apply_fn(fit, data$contrasts, tau, TRUE, TRUE)
        p_column <- "P.Value"
    } else {
        result <- runTestsContrasts(data$y, contrast_strings = data$contrasts, design_matrix = data$design,
            formula_string = "~ 0 + group", treat_lfc_cutoff = tau, eBayes_trend = TRUE, eBayes_robust = TRUE)
        p_column <- "raw_pvalue"
    }
    ids <- rownames(data$y)
    tables <- lapply(stats::setNames(data$contrasts, data$contrasts), function(contrast) {
        table <- as.data.frame(result$results[[contrast]])
        at <- match(ids, rownames(table))
        list(p = table[[p_column]][at], q = table$fdr_qvalue[at],
            untested = if ("AveExpr" %in% names(table)) sum(rownames(table) %in% ids & is.na(table$AveExpr)) else NA_integer_)
    })
    fit <- result$fit.eb
    list(tables = tables, s2prior_distinct = length(unique(fit$s2.prior)), zero_df = sum(fit$df.residual == 0, na.rm = TRUE))
}

k201Rate <- function(hit, keep) if (sum(keep)) sum(hit & keep) / sum(keep) else NA_real_

# One row per contrast for one data set, set and threshold.
k201Metrics <- function(data, set, tau, out, base_pv) {
    rows <- list()
    for (contrast in data$contrasts) {
        tab <- out$tables[[contrast]]
        tested <- !is.na(tab$p)
        called <- tested & !is.na(tab$q) & tab$q < 0.05
        beta <- abs(data$beta[, contrast])
        sig <- tested & tab$p < K201_ALPHA
        row <- list(tau = tau, contrast = contrast, n_features = length(tab$p), n_untested = tab$untested,
            n_na_p = sum(!tested), n_zero_df = out$zero_df, s2prior_distinct = out$s2prior_distinct)
        if (identical(set, "main")) {
            false <- called & beta <= tau
            true <- called & beta > tau
            row <- c(row, list(calls = sum(called), false_calls = sum(false), true_calls = sum(true),
                fdp = sum(false) / max(sum(called), 1), power = k201Rate(true, tested & beta > tau)))
            for (k in 1:3) row[[paste0("fdp_third", k)]] <- sum(false & data$third == k) / max(sum(called & data$third == k), 1)
        }
        null_set <- if (identical(set, "calib")) tested else if (tau == 0) tested & beta == 0 else NULL
        if (!is.null(null_set)) {
            strata <- list(overall = rep(TRUE, length(beta)), third1 = data$third == 1L, third2 = data$third == 2L,
                third3 = data$third == 3L, missing = if (data$vd20) data$with_missing else NULL)
            for (s in names(strata)) {
                if (is.null(strata[[s]])) next
                keep <- null_set & strata[[s]]
                row[[paste0("null_", s)]] <- k201Rate(sig, keep)
                row[[paste0("den_", s)]] <- sum(keep)
            }
        }
        if (!is.null(base_pv)) {
            other <- base_pv[[contrast]]
            both <- tested & !is.na(other$p)
            row$max_dp <- if (any(both)) max(abs(tab$p[both] - other$p[both])) else NA_real_
            other_called <- !is.na(other$p) & !is.na(other$q) & other$q < 0.05
            row$only_base <- sum(other_called & !called)
            row$only_head <- sum(called & !other_called)
        }
        rows[[length(rows) + 1L]] <- row
    }
    rows
}

# ---- child process: one unit ------------------------------------------------------------------------------------------
k201Child <- function(a) {
    stopifnot(length(a) == 12L)
    side <- a[[2L]]
    tree <- normalizePath(a[[3L]])
    commit <- a[[4L]]
    kernel <- a[[5L]]
    cell <- as.integer(a[[6L]])
    set <- a[[7L]]
    rs <- seq(as.integer(a[[8L]]), as.integer(a[[9L]]))
    rds <- a[[10L]]
    pv_dir <- a[[11L]]
    fail_r <- as.integer(a[[12L]])
    setwd(tree)
    suppressPackageStartupMessages(suppressMessages(pkgload::load_all(tree, export_all = TRUE, helpers = FALSE,
        compile = FALSE, quiet = TRUE)))
    ns_path <- normalizePath(getNamespaceInfo("MultiScholaR", "path"))
    if (!identical(ns_path, tree)) stop(sprintf("MultiScholaR is loaded from %s, not from the side's tree %s", ns_path, tree))
    for (namespace in c("global", "MultiScholaR")) logger::log_threshold(logger::ERROR, namespace = namespace)
    taus <- K201_KERNELS[[kernel]]$taus
    if (identical(set, "calib")) taus <- taus[taus > 0]
    pv_file <- file.path(pv_dir, sprintf("%s-%02d-%s-%d.rds", kernel, cell, set, rs[[1L]]))
    base_pv <- if (identical(side, "head") && file.exists(pv_file)) readRDS(pv_file) else NULL
    own_pv <- list()
    rows <- list()
    for (r in rs) for (tau in taus) {
        key <- sprintf("%d|%s", r, tau)
        started <- proc.time()[[3L]]
        data <- k201Generate(kernel, cell, set, r, if (identical(set, "calib")) tau else NULL)
        base <- list(side = side, commit = commit, kernel = kernel, cell = cell, set = set, r = r, hash = data$hash)
        outcome <- tryCatch({
            if (identical(r, fail_r)) stop("KIT-201 R0 forced failure (lifecycle failure path)")
            out <- suppressWarnings(suppressMessages(k201Kernel(kernel, data, tau)))
            seconds <- proc.time()[[3L]] - started
            own_pv[[key]] <- lapply(out$tables, function(tab) list(p = tab$p, q = tab$q))
            lapply(k201Metrics(data, set, tau, out, base_pv[[key]]), function(row)
                c(base, row, list(status = "ok", message = "", seconds = seconds)))
        }, error = function(e) list(c(base, list(tau = tau, contrast = NA_character_, status = "failed",
            message = k201Flat(conditionMessage(e)), seconds = proc.time()[[3L]] - started))))
        rows <- c(rows, outcome)
    }
    if (identical(side, "base")) saveRDS(own_pv, pv_file)
    table <- do.call(rbind, lapply(rows, function(row) {
        values <- stats::setNames(as.list(rep(NA, length(K201_ROW_COLUMNS))), K201_ROW_COLUMNS)
        given <- intersect(names(row), K201_ROW_COLUMNS)
        values[given] <- lapply(row[given], function(v) if (is.null(v) || !length(v)) NA else v[[1L]])
        as.data.frame(values, stringsAsFactors = FALSE)
    }))
    saveRDS(list(rows = table, ns_path = ns_path), rds)
}

# ---- parent -----------------------------------------------------------------------------------------------------------
k201Prepare <- function(tree) {
    log <- tempfile("kit201-prepare-", fileext = ".txt")
    if (file.exists(file.path(tree, "src", "build_dia_rsum.R")) &&
        !file.exists(file.path(tree, "inst", "duckdb", "linux_amd64", "manifest.rds"))) {
        status <- system(sprintf("cd %s && %s --vanilla src/build_dia_rsum.R > %s 2>&1", shQuote(tree), shQuote(k201_rscript), shQuote(log)))
        if (!identical(as.integer(status), 0L)) stop(sprintf("build_dia_rsum.R exit %s: %s", status, k201Flat(utils::tail(readLines(log), 5L))))
    }
    status <- system2(k201_rscript, c("-e", shQuote(sprintf("pkgbuild::compile_dll('%s', force = TRUE, quiet = TRUE)", tree))),
        stdout = log, stderr = log)
    if (!identical(as.integer(status), 0L)) stop(sprintf("compile_dll exit %s: %s", status, k201Flat(utils::tail(readLines(log), 5L))))
}

# Units: one per kernel, cell, set and block of at most K201_CHUNK data sets.
k201Units <- function(kernel, cells, R_of, r0 = FALSE) {
    units <- list()
    for (cell in cells) for (set in c("main", "calib")) {
        rs <- if (r0) K201_R0 else (if (identical(set, "main")) 0L else 5000L) + seq_len(R_of(cell))
        blocks <- split(rs, ceiling(seq_along(rs) / K201_CHUNK))
        for (b in blocks) units[[length(units) + 1L]] <- data.frame(kernel = kernel, cell = cell, set = set,
            from = min(b), to = max(b), fail = NA_integer_, stringsAsFactors = FALSE)
    }
    do.call(rbind, units)
}

k201RunSide <- function(units, side, tree, commit, work) {
    pv_dir <- file.path(work, "pv")
    dir.create(pv_dir, recursive = TRUE, showWarnings = FALSE)
    logs <- file.path(work, paste0("units-", side))
    dir.create(logs, recursive = TRUE, showWarnings = FALSE)
    runChild <- function(i) {
        u <- units[i, ]
        tag <- sprintf("%s-%02d-%s-%d%s", u$kernel, u$cell, u$set, u$from, if (is.na(u$fail)) "" else "-forced")
        rds <- file.path(logs, paste0(tag, ".rds"))
        log <- file.path(logs, paste0(tag, ".log"))
        code <- system2(k201_rscript, c(shQuote(k201_file), "--unit", side, shQuote(tree), commit, u$kernel, u$cell, u$set,
            u$from, u$to, shQuote(rds), shQuote(pv_dir), if (is.na(u$fail)) "0" else u$fail), stdout = log, stderr = log)
        res <- if (file.exists(rds)) readRDS(rds) else list(error = sprintf("child exit %s, no result: %s", code,
            k201Flat(utils::tail(readLines(log, warn = FALSE), 25L), 4000L)))
        res$forced <- !is.na(u$fail)
        res
    }
    results <- parallel::mclapply(seq_len(nrow(units)), runChild, mc.cores = K201_WORKERS, mc.preschedule = FALSE)
    failed <- vapply(results, function(res) !is.list(res) || !is.null(res$error), logical(1))
    if (any(failed)) {
        messages <- vapply(results[failed], function(res) if (is.list(res)) res$error else k201Flat(res), character(1))
        stop(sprintf("%d unit(s) on %s returned no records (instrument defect): %s", sum(failed), side,
            k201Flat(messages[[1L]], 2000L)))
    }
    if (!all(vapply(results, function(res) identical(res$ns_path, normalizePath(tree)), logical(1))))
        stop(sprintf("a child on %s did not load the side's tree", side))
    rows <- do.call(rbind, lapply(results[!vapply(results, `[[`, logical(1), "forced")], `[[`, "rows"))
    forced <- results[vapply(results, `[[`, logical(1), "forced")]
    list(rows = rows, forced = if (length(forced)) do.call(rbind, lapply(forced, `[[`, "rows")))
}

# Both sides for a set of units; base first (its p-values are the head children's comparison).
k201RunBoth <- function(units, base_commit, head_root, head_commit, forced_unit = NULL) {
    work <- file.path(tempdir(), "kit201-work")
    dir.create(work, recursive = TRUE, showWarnings = FALSE)
    base_tree <- file.path(tempdir(), "kit201-base")
    if (!dir.exists(base_tree)) {
        invisible(system2("git", c("worktree", "add", "--detach", shQuote(base_tree), base_commit), stdout = TRUE, stderr = TRUE))
    }
    if (!dir.exists(file.path(base_tree, "R"))) stop("the base worktree could not be made")
    base_tree <- normalizePath(base_tree)
    k201Prepare(base_tree)
    k201Prepare(head_root)
    base <- k201RunSide(units, "base", base_tree, base_commit, work)
    head <- k201RunSide(rbind(units, forced_unit), "head", head_root, head_commit, work)
    invisible(system2("git", c("worktree", "remove", "--force", shQuote(base_tree)), stdout = TRUE, stderr = TRUE))
    list(rows = rbind(base$rows, head$rows), forced = head$forced)
}

# ---- summary and decision ---------------------------------------------------------------------------------------------
k201Boot <- function(x) {
    set.seed(K201_BOOT_SEED, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    means <- replicate(K201_BOOT_N, mean(x[sample.int(length(x), length(x), replace = TRUE)]))
    stats::quantile(means, c(0.025, 0.975), type = 7, names = FALSE)
}

# Per data set, the two sides' values of one comparison (FDP on the main sets; the null rate on the family's sets).
k201Paired <- function(rows, kernel, cell, tau, contrast, rule, stratum, data_set) {
    keep <- rows$kernel == kernel & rows$cell == cell & rows$set == data_set & rows$status == "ok" &
        as.character(rows$tau) == tau & rows$contrast == contrast
    column <- if (identical(rule, "fdp")) "fdp" else paste0("null_", stratum)
    side <- function(name) {
        part <- rows[keep & rows$side == name, , drop = FALSE]
        stats::setNames(as.numeric(part[[column]]), part$r)
    }
    b <- side("base")
    h <- side("head")
    common <- intersect(names(b), names(h))
    ok <- common[!is.na(b[common]) & !is.na(h[common])]
    list(base = unname(b[ok]), head = unname(h[ok]))
}

# Cell status: stopped when any of its data sets failed on either side, or the data sets' hashes differ (instrument).
k201CellStatus <- function(rows, kernel, cell) {
    part <- rows[rows$kernel == kernel & rows$cell == cell, , drop = FALSE]
    failed <- unique(part[part$status != "ok", c("side", "set", "r", "message")])
    hashes <- unique(part[, c("side", "set", "r", "tau", "hash")])
    by <- split(hashes$hash, paste(hashes$set, hashes$r, hashes$tau))
    mismatch <- sum(vapply(by, function(h) length(unique(h)) > 1L, logical(1)))
    if (mismatch) return(list(status = "stopped", reason = sprintf("instrument: %d data set hash mismatch(es)", mismatch)))
    if (nrow(failed)) return(list(status = "stopped", reason = sprintf("%d failed data set(s): %s", nrow(failed),
        k201Flat(failed$message[[1L]], 300L))))
    list(status = "ok", reason = "")
}

k201Summary <- function(rows, family, precision = NULL) {
    out <- list()
    for (i in seq_len(nrow(family))) {
        f <- family[i, ]
        status <- k201CellStatus(rows, f$kernel, f$cell)
        values <- k201Paired(rows, f$kernel, f$cell, f$tau, f$contrast, f$rule, f$stratum, f$data)
        gated <- if (is.null(precision)) NA else {
            hit <- precision$kernel == f$kernel & precision$cell == f$cell & precision$tau == f$tau &
                precision$contrast == f$contrast & precision$rule == f$rule & precision$stratum == f$stratum
            isTRUE(precision$gated[hit])
        }
        row <- c(as.list(f), list(cell_status = status$status, cell_reason = status$reason, n = length(values$head),
            gated = gated))
        if (identical(status$status, "ok") && length(values$head) >= 2L) {
            d <- values$head - values$base
            ci <- k201Boot(d)
            row <- c(row, list(mean_base = mean(values$base), mean_head = mean(values$head), diff = mean(d),
                diff_lo = ci[[1L]], diff_hi = ci[[2L]], mc_se_diff = stats::sd(d) / sqrt(length(d))))
            if (identical(f$rule, "null")) {
                ci0 <- k201Boot(values$head - K201_ALPHA)
                row <- c(row, list(excess = mean(values$head) - K201_ALPHA, excess_lo = ci0[[1L]], excess_hi = ci0[[2L]],
                    mc_se_head = stats::sd(values$head) / sqrt(length(values$head))))
                pass <- row$diff_hi <= K201_DELTA || row$excess_hi <= K201_DELTA
                fail <- row$diff > K201_DELTA && row$excess > K201_DELTA
            } else {
                pass <- row$diff_hi <= K201_DELTA
                fail <- row$diff > K201_DELTA
            }
            row$verdict <- if (pass) "pass" else if (fail) "fails_materially" else "imprecise"
        } else row$verdict <- "not_judged"
        out[[i]] <- row
    }
    columns <- unique(unlist(lapply(out, names)))
    do.call(rbind, lapply(out, function(row) {
        values <- stats::setNames(as.list(rep(NA, length(columns))), columns)
        values[names(row)] <- lapply(row, function(v) if (is.null(v) || !length(v)) NA else v[[1L]])
        as.data.frame(values, stringsAsFactors = FALSE)
    }))
}

# The reported measures (not gated): power, FDP by third, call changes, |delta p|, zero-df and untested counts.
k201Reported <- function(rows) {
    ok <- rows[rows$status == "ok" & rows$set == "main", , drop = FALSE]
    measures <- c("power", "fdp_third1", "fdp_third2", "fdp_third3", "max_dp", "only_base", "only_head", "n_zero_df",
        "n_untested", "n_na_p", "calls")
    keys <- unique(ok[, c("kernel", "cell", "tau", "contrast")])
    out <- list()
    for (i in seq_len(nrow(keys))) {
        k <- keys[i, ]
        part <- ok[ok$kernel == k$kernel & ok$cell == k$cell & ok$tau == k$tau & ok$contrast == k$contrast, , drop = FALSE]
        for (m in measures) for (side in c("base", "head")) {
            v <- as.numeric(part[[m]][part$side == side])
            if (!length(v) || all(is.na(v))) next
            out[[length(out) + 1L]] <- data.frame(kernel = k$kernel, cell = k$cell, tau = as.character(k$tau),
                contrast = k$contrast, measure = m, side = side, mean = mean(v, na.rm = TRUE), max = max(v, na.rm = TRUE),
                n = sum(!is.na(v)), stringsAsFactors = FALSE)
        }
    }
    do.call(rbind, out)
}

k201Classify <- function(summary) {
    gated <- summary[summary$gated %in% TRUE, , drop = FALSE]
    if (!nrow(gated)) return("no gated comparison")
    if (any(gated$cell_status != "ok" | gated$verdict == "not_judged")) return("stopped (a gated comparison was not judged)")
    if (any(gated$verdict == "fails_materially")) return("worse")
    if (any(gated$verdict == "imprecise")) return("inconclusive")
    "no worse"
}

# ---- R0 (blind) -------------------------------------------------------------------------------------------------------
k201R0 <- function(base_commit, out) {
    head_root <- normalizePath(".")
    head_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
    family <- k201Family()
    cells <- k201Cells()$cell
    units <- do.call(rbind, lapply(names(K201_KERNELS), function(kernel) k201Units(kernel, cells, NULL, r0 = TRUE)))
    forced <- data.frame(kernel = K201_FORCED$kernel, cell = K201_FORCED$cell, set = "main", from = min(K201_R0),
        to = max(K201_R0), fail = K201_FORCED$r, stringsAsFactors = FALSE)
    run <- k201RunBoth(units, base_commit, head_root, head_commit, forced)
    rows <- run$rows
    checks <- list()
    check <- function(kernel, name, pass, detail = "") checks[[length(checks) + 1L]] <<- data.frame(kernel = kernel,
        check = name, status = if (isTRUE(pass)) "pass" else "fail", detail = detail, stringsAsFactors = FALSE)
    for (kernel in names(K201_KERNELS)) {
        part <- rows[rows$kernel == kernel, , drop = FALSE]
        complete <- all(vapply(cells, function(cell) identical(k201CellStatus(part, kernel, cell)$status, "ok"), logical(1)))
        check(kernel, "every cell completes on both sides", complete && all(c("base", "head") %in% part$side))
        hashes <- unique(part[, c("side", "set", "r", "tau", "cell", "hash")])
        by <- split(hashes$hash, paste(hashes$cell, hashes$set, hashes$r, hashes$tau))
        check(kernel, "data-set hashes agree", all(vapply(by, function(h) length(unique(h)) == 1L && length(h) == 2L, logical(1))))
        engaged <- part$status == "ok" & (as.numeric(part$tau) > 0 | identical(kernel, "P"))
        check(kernel, "distinct prior variances on head (tau > 0; P at 0 too)",
            all(part$s2prior_distinct[engaged & part$side == "head"] > 1))
        check(kernel, "one prior variance on base (tau > 0; P at 0 too)",
            all(part$s2prior_distinct[engaged & part$side == "base"] == 1))
        summary <- k201Summary(part, family[family$kernel == kernel, , drop = FALSE])
        check(kernel, "summary written (structure only)", nrow(summary) == sum(family$kernel == kernel),
            sprintf("%d rows; columns: %s", nrow(summary), paste(names(summary), collapse = ",")))
    }
    forced_rows <- run$forced
    forced_row <- forced_rows[forced_rows$r == K201_FORCED$r & forced_rows$status == "failed", , drop = FALSE]
    forced_all <- rbind(rows[!(rows$kernel == K201_FORCED$kernel & rows$cell == K201_FORCED$cell & rows$side == "head"), ],
        forced_rows)
    forced_status <- k201CellStatus(forced_all, K201_FORCED$kernel, K201_FORCED$cell)
    check(K201_FORCED$kernel, "failure path: the forced data set is recorded failed with its message",
        nrow(forced_row) >= 1L && all(grepl("forced failure", forced_row$message)))
    check(K201_FORCED$kernel, "failure path: the forced cell's summary is stopped", identical(forced_status$status, "stopped"))
    checks <- do.call(rbind, checks)
    utils::write.table(checks, file.path(out, "r0-checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    # SDs only (no mean, no difference, no level): per gated comparison, of the paired difference and of head's rate.
    sds <- list()
    for (i in seq_len(nrow(family))) {
        f <- family[i, ]
        values <- k201Paired(rows, f$kernel, f$cell, f$tau, f$contrast, f$rule, f$stratum, f$data)
        n <- length(values$head)
        sds[[i]] <- data.frame(kind = "sd", kernel = f$kernel, cell = f$cell, side = NA_character_, tau = f$tau,
            contrast = f$contrast, rule = f$rule, stratum = f$stratum, n = n,
            sd_diff = if (n >= 2L) stats::sd(values$head - values$base) else NA_real_,
            sd_head = if (n >= 2L && identical(f$rule, "null")) stats::sd(values$head) else NA_real_,
            seconds_per_dataset = NA_real_, stringsAsFactors = FALSE)
    }
    ok <- rows[rows$status == "ok", , drop = FALSE]
    # seconds per data set: one kernel call per (set, r, tau); the rows of a call's contrasts carry the same seconds.
    calls <- unique(ok[, c("kernel", "cell", "side", "set", "r", "tau", "seconds")])
    timing <- stats::aggregate(seconds ~ kernel + cell + side + set + r, calls, sum)
    timing <- stats::aggregate(seconds ~ kernel + cell + side + set, timing, mean)
    timing <- stats::aggregate(seconds ~ kernel + cell + side, timing, sum)
    times <- data.frame(kind = "time", kernel = timing$kernel, cell = timing$cell, side = timing$side, tau = NA_character_,
        contrast = NA_character_, rule = NA_character_, stratum = NA_character_, n = NA_integer_, sd_diff = NA_real_,
        sd_head = NA_real_, seconds_per_dataset = timing$seconds, stringsAsFactors = FALSE)
    r0sd <- rbind(do.call(rbind, sds), times)
    utils::write.table(r0sd, file.path(out, "r0-sd.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    k201Precision(r0sd, family, out)
    cat("KIT-201 R0: checks", paste(checks$status, collapse = " "), "\n")
    if (any(checks$status != "pass")) quit(status = 2L)
}

# The formula of 5.2 (no hand choice): R per cell from R0's SDs, capped at 1,000 and by the 180.0 core-hour cap.
k201Precision <- function(r0sd, family, out) {
    K <- nrow(family)
    se_star <- k201SeStar(K)
    sd <- r0sd[r0sd$kind == "sd", , drop = FALSE]
    sd$sd_max <- pmax(sd$sd_diff, sd$sd_head, na.rm = TRUE)
    cells <- unique(sd[, c("kernel", "cell")])
    cells$sd_max <- mapply(function(k, c) {
        v <- sd$sd_max[sd$kernel == k & sd$cell == c]
        if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
    }, cells$kernel, cells$cell)
    cells$R <- pmin(K201_R_MAX, pmax(K201_R_MIN, ceiling((1.5 * cells$sd_max / se_star)^2)))
    cells$R[is.na(cells$R)] <- K201_R_MIN
    time <- r0sd[r0sd$kind == "time", , drop = FALSE]
    cells$sec_per_dataset <- mapply(function(k, c) sum(time$seconds_per_dataset[time$kernel == k & time$cell == c]),
        cells$kernel, cells$cell)
    cost <- function(cells) {
        worker <- cells$R * cells$sec_per_dataset * 1.15
        wall_min <- tapply(worker, cells$kernel, sum) / K201_WORKERS / 60
        jobs <- sum(ceiling(wall_min / 170))
        list(total = sum(worker) / 3600 + jobs * 10 * 4 / 60, cell = worker / 3600, jobs = jobs)
    }
    projected <- cost(cells)
    while (projected$total > K201_CAP) {
        open <- which(cells$R > K201_R_MIN)
        if (!length(open)) stop("the 180.0 core-hour cap cannot be met at R = 200 (reported to the owner)")
        worst <- open[which.max(projected$cell[open])]
        cells$R[worst] <- max(K201_R_MIN, cells$R[worst] - ceiling(0.05 * cells$R[worst]))
        projected <- cost(cells)
    }
    at <- match(paste(sd$kernel, sd$cell), paste(cells$kernel, cells$cell))
    sd$R <- cells$R[at]
    sd$projected_se <- 1.5 * sd$sd_max / sqrt(sd$R)
    sd$gated <- !is.na(sd$projected_se) & sd$projected_se <= se_star
    precision <- data.frame(kernel = sd$kernel, cell = sd$cell, tau = sd$tau, contrast = sd$contrast, rule = sd$rule,
        stratum = sd$stratum, R = sd$R, projected_se = sd$projected_se, gated = sd$gated,
        cell_core_hours = projected$cell[at], stringsAsFactors = FALSE)
    utils::write.table(precision, file.path(out, "precision.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    writeLines(c(sprintf("K\t%d", K), sprintf("se_star\t%.6f", se_star), sprintf("projected_core_hours\t%.2f", projected$total),
        sprintf("jobs\t%d", projected$jobs), sprintf("gated\t%d", sum(sd$gated)), sprintf("reported_not_gated\t%d", sum(!sd$gated))),
        file.path(out, "precision-cost.tsv"))
}

# ---- M1 ---------------------------------------------------------------------------------------------------------------
k201M1 <- function(kernel, base_commit, out, cells_arg, part) {
    precision <- utils::read.delim(file.path(out, "precision.tsv"), stringsAsFactors = FALSE, colClasses = c(tau = "character"))
    family <- k201Family()
    cells <- if (is.null(cells_arg)) k201Cells()$cell else as.integer(strsplit(cells_arg, ",")[[1L]])
    R_of <- function(cell) unique(precision$R[precision$kernel == kernel & precision$cell == cell])[[1L]]
    units <- k201Units(kernel, cells, R_of)
    head_root <- normalizePath(".")
    head_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
    run <- k201RunBoth(units, base_commit, head_root, head_commit)
    saveRDS(run$rows, file.path(out, sprintf("raw-%s%s.rds", kernel, if (is.null(part)) "" else paste0("-", part))))
    if (is.null(cells_arg)) k201Summarise(kernel, out)
}

k201Summarise <- function(kernel, out) {
    files <- list.files(out, sprintf("^raw-%s(-.+)?[.]rds$", kernel), full.names = TRUE)
    rows <- do.call(rbind, lapply(files, readRDS))
    precision <- utils::read.delim(file.path(out, "precision.tsv"), stringsAsFactors = FALSE, colClasses = c(tau = "character"))
    family <- k201Family()
    family <- family[family$kernel == kernel, , drop = FALSE]
    gz <- gzfile(file.path(out, sprintf("rows-%s.tsv.gz", kernel)), "w")
    utils::write.table(rows, gz, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    close(gz)
    summary <- k201Summary(rows, family, precision)
    utils::write.table(summary, file.path(out, sprintf("summary-%s.tsv", kernel)), sep = "\t", quote = FALSE,
        row.names = FALSE, na = "NA")
    utils::write.table(k201Reported(rows), file.path(out, sprintf("reported-%s.tsv", kernel)), sep = "\t", quote = FALSE,
        row.names = FALSE, na = "NA")
    gated <- summary[summary$gated %in% TRUE, , drop = FALSE]
    v1v2 <- rows[rows$status == "ok" & rows$side == "head" & rows$kernel == kernel, , drop = FALSE]
    cells <- k201Cells()
    falsifier <- v1v2$cell %in% cells$cell[cells$variance %in% c("V1", "V2")]
    lines <- c(sprintf("kernel\t%s", kernel), sprintf("class\t%s", k201Classify(summary)),
        sprintf("gated\t%d", nrow(gated)), sprintf("pass\t%d", sum(gated$verdict == "pass")),
        sprintf("imprecise\t%d", sum(gated$verdict == "imprecise")),
        sprintf("fails_materially\t%d", sum(gated$verdict == "fails_materially")),
        sprintf("not_judged\t%d", sum(gated$verdict == "not_judged")),
        sprintf("stopped_cells\t%s", paste(unique(summary$cell[summary$cell_status != "ok"]), collapse = ",")),
        sprintf("falsifier_max_dp_V1V2\t%s", format(suppressWarnings(max(v1v2$max_dp[falsifier], na.rm = TRUE)), digits = 6)))
    writeLines(lines, file.path(out, sprintf("class-%s.tsv", kernel)))
    cat(paste("KIT-201 M1", lines), sep = "\n")
}

# ---- entry ------------------------------------------------------------------------------------------------------------
k201Main <- function(a) {
    if (length(a) && identical(a[[1L]], "--unit")) return(invisible(k201Child(a)))
    flags <- a[startsWith(a, "--")]
    positional <- a[!startsWith(a, "--")]
    flag <- function(name) {
        hit <- flags[startsWith(flags, paste0("--", name, "="))]
        if (length(hit)) sub("^[^=]*=", "", hit[[1L]]) else NULL
    }
    mode <- positional[[1L]]
    default_out <- "docs/qa/workflow-execution/kit/kit-201/m1"
    outOf <- function(i) {
        out <- if (length(positional) >= i) positional[[i]] else default_out
        dir.create(out, recursive = TRUE, showWarnings = FALSE)
        normalizePath(out)
    }
    if (identical(mode, "family")) {
        out <- outOf(2L)
        utils::write.table(k201Family(), file.path(out, "family.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    } else if (identical(mode, "r0")) {
        stopifnot(grepl("^[0-9a-f]{40}$", positional[[2L]]))
        k201R0(positional[[2L]], outOf(3L))
    } else if (identical(mode, "m1")) {
        stopifnot(positional[[2L]] %in% names(K201_KERNELS), grepl("^[0-9a-f]{40}$", positional[[3L]]))
        out <- outOf(4L)
        committed <- utils::read.delim(file.path(out, "family.tsv"), stringsAsFactors = FALSE, colClasses = "character")
        mine <- k201Family()
        mine[] <- lapply(mine, as.character)
        if (!identical(unname(as.list(committed)), unname(as.list(mine)))) stop("family.tsv differs from this script's family")
        k201M1(positional[[2L]], positional[[3L]], out, flag("cells"), flag("part"))
    } else if (identical(mode, "summary")) {
        stopifnot(positional[[2L]] %in% names(K201_KERNELS))
        k201Summarise(positional[[2L]], outOf(3L))
    } else stop("mode must be family, r0, m1 or summary")
}

k201Main(k201_argv)
