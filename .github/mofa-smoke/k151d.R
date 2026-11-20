# KIT-151 d evidence (branch run-mofa-evidence-k151d; methods ledger INT-13). Three questions, one task each.
# Declarations: docs/qa/workflow-execution/kit/kit-151d/{ruv,z,pin}-declaration.md (written after the lifecycle run,
# committed before each counted run).
#   ruv (Q1)  which RUV output should feed integration: a planted shared within-group factor S, one unwanted factor
#             per omic, a group effect; MOFA2 inputs uncorrected, ideal (the planted unwanted term removed exactly),
#             RUV-III with groups as replicate sets and every feature a control at k = 1, 2, 3
#             (RUVIIIC::RUVIII_C_Varying, as the package calls it), and the package's calibrated default design
#             (the same correction with k from the package's estimated count: sva::num.sv(method = "be", B = 200,
#             seed 20261118) on the features observed in every sample, the groups removed; R/func_general_ruv_factor_count.R,
#             R/func_omics_ruv_factors.R at e8a2ac67). Scenarios: unwanted factor unrelated to the groups (main) and
#             partly confounded with them (conf: correlation 0.6 with the group in the sample).
#   z   (Q2)  per-feature z-score and scale_views on the k151b generator (mean-variance trend; abundance-dependent
#             missing values on or off; data-set scales none, as k151b, large data set scaled up, or down), selection
#             fixed at the decided default (INT-12: 1,000 features, variance minus the fitted mean-variance trend,
#             observed in at least half the samples).
#   pin (Q3)  the k151c pin case (n20weak, data seed 151007, fast, the eight first-run arms) and two hard slow cases
#             of the k151b generator (prev@10 no selection, the case that diverged between an Intel and an AMD
#             processor; full@23 at the default selection), unpinned (U), pinned (P) and pinned with numpy's AVX-512
#             targets off (PN; the k151c fast case only).
# MOFA2 is called directly as in k151b.R / k151c.R. Exit status 1 when any cell failed (rows with status "failed").
#
# Usage:
#   Rscript --vanilla k151d.R <out> ruv <N> <scenarios,...> <strength> <seed indices,...> <arms,...>
#   Rscript --vanilla k151d.R <out> z   <N> <seed index> <cells,...>
#   Rscript --vanilla k151d.R <out> pin <condition U|P|PN> <cases,...>

arguments <- commandArgs(trailingOnly = TRUE)
out_dir <- arguments[[1L]]
task <- arguments[[2L]]
stopifnot(task %in% c("ruv", "z", "pin"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()
split_list <- function(x) strsplit(x, ",", fixed = TRUE)[[1L]]
TRAIN_SEED <- 20261119L

# ================================================================== k151c generator and preparation (pin fast case)
K151_SIZES <- c(D1 = 20000L, D2 = 2000L, D3 = 200L, D4 = 100L)
K151_ACTIVE <- list(D1 = c("F1", "F3"), D2 = "F1", D3 = c("F1", "F2"), D4 = c("F1", "F2"))
K151_LOADED_SHARE <- 0.10
K151_EFFECT <- 1.0
K151_NOISE_SD <- 1.5
K151_F1_NOISE_SD <- 0.5

generate_k151 <- function(scenario, seed_index) {
    N <- if (startsWith(scenario, "n60")) 60L else 20L
    F2_EFFECT <- if (grepl("strong$", scenario)) 1.5 else 1.0
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(151000L + seed_index)
    group <- rep(c(-1, 1), each = N / 2L)
    f1 <- group + stats::rnorm(N, 0, K151_F1_NOISE_SD)
    f2 <- stats::rnorm(N)
    f3 <- stats::rnorm(N)
    z <- scale(cbind(F1 = f1, F2 = f2, F3 = f3))
    samples <- sprintf("S%02d", seq_len(N))
    rownames(z) <- samples
    views <- list()
    for (view in names(K151_SIZES)) {
        d <- K151_SIZES[[view]]
        w <- matrix(0, d, 3L, dimnames = list(NULL, colnames(z)))
        for (k in K151_ACTIVE[[view]]) {
            loaded <- sample.int(d, round(K151_LOADED_SHARE * d))
            w[loaded, k] <- sample(c(-1, 1), length(loaded), replace = TRUE) * (if (k == "F2") F2_EFFECT else K151_EFFECT)
        }
        offset <- stats::rnorm(d, 20, 2)
        x <- w %*% t(z) + matrix(stats::rnorm(d * N, 0, K151_NOISE_SD), d, N) + offset
        dimnames(x) <- list(sprintf("%s_f%05d", view, seq_len(d)), samples)
        views[[view]] <- x
    }
    list(views = views, z = z, N = N, data_seed = 151000L + seed_index)
}

k151_arm <- function(cap = NA_integer_, zscore = FALSE, scale_views = TRUE, weight_views = FALSE) {
    list(cap = cap, zscore = zscore, scale_views = scale_views, center_groups = !zscore, weight_views = weight_views)
}
K151_ARMS <- list(
    none = k151_arm(scale_views = FALSE),
    scale_views = k151_arm(),
    cap1000 = k151_arm(cap = 1000L),
    cap2000 = k151_arm(cap = 2000L),
    cap5000 = k151_arm(cap = 5000L),
    cap10000 = k151_arm(cap = 10000L),
    weight_views = k151_arm(weight_views = TRUE),
    zscore = k151_arm(zscore = TRUE, scale_views = FALSE))

prepare_k151 <- function(views, arm) {
    out <- list()
    for (view in names(views)) {
        m <- views[[view]]
        m <- m[rowSums(!is.na(m)) >= 2L, , drop = FALSE]
        if (arm$zscore) {
            spread <- apply(m, 1L, stats::sd, na.rm = TRUE)
            m <- m[spread > 0 & !is.na(spread), , drop = FALSE]
        }
        if (!is.na(arm$cap) && nrow(m) > arm$cap) {
            variance <- apply(m, 1L, stats::var, na.rm = TRUE)
            keep <- sort(order(-variance, seq_along(variance))[seq_len(arm$cap)])
            m <- m[keep, , drop = FALSE]
        }
        if (arm$zscore) {
            centre <- rowMeans(m, na.rm = TRUE)
            spread <- apply(m, 1L, stats::sd, na.rm = TRUE)
            m <- (m - centre) / spread
        }
        rownames(m) <- paste(view, rownames(m), sep = ":")
        out[[view]] <- m
    }
    out
}

# ================================================================== k151b generator (z task, pin hard cases)
# k151b.R's generate() with N and the data-set scales as parameters; N = 20 with scales c(1, 2, 0.25, 4) is k151b's
# arithmetic exactly (same draws in the same order).
B_SIZES <- c(D1 = 20000L, D2 = 2000L, D3 = 200L, D4 = 100L)
B_ACTIVE <- list(D1 = c("F1", "F3"), D2 = "F1", D3 = c("F1", "F2"), D4 = c("F1", "F2"))
B_LOADED_SHARE <- 0.10
B_EFFECT <- 1.0
B_F2_EFFECT <- 1.5
B_NOISE_SD <- 1.5
B_F1_NOISE_SD <- 0.5
B_TREND_REF <- 20
B_TREND_DOUBLING <- 4
B_MISS_CENTRE <- 17.5
B_MISS_WIDTH <- 0.5
B_MISS_FLOOR <- 0.02
B_HALF_SHARE <- 0.5
B_TREND_SPAN <- 0.3
B_SCALES <- list(
    none = c(D1 = 1, D2 = 1, D3 = 1, D4 = 1),
    orig = c(D1 = 1, D2 = 2, D3 = 0.25, D4 = 4), # k151b
    lup = c(D1 = 4, D2 = 2, D3 = 0.5, D4 = 0.25), # large data sets scaled up, small (F2) ones down
    ldown = c(D1 = 0.25, D2 = 0.5, D3 = 2, D4 = 4)) # large data sets scaled down, small (F2) ones up

generate_b <- function(N, seed_index, trend, missing, scales) {
    data_seed <- 151000L + seed_index
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(data_seed)
    group <- rep(c(-1, 1), each = N / 2L)
    f1 <- group + stats::rnorm(N, 0, B_F1_NOISE_SD)
    f2 <- stats::rnorm(N)
    f3 <- stats::rnorm(N)
    z <- scale(cbind(F1 = f1, F2 = f2, F3 = f3))
    samples <- sprintf("S%02d", seq_len(N))
    rownames(z) <- samples
    views <- list()
    for (view in names(B_SIZES)) {
        d <- B_SIZES[[view]]
        w <- matrix(0, d, 3L, dimnames = list(NULL, colnames(z)))
        for (k in B_ACTIVE[[view]]) {
            loaded <- sample.int(d, round(B_LOADED_SHARE * d))
            w[loaded, k] <- sample(c(-1, 1), length(loaded), replace = TRUE) * (if (k == "F2") B_F2_EFFECT else B_EFFECT)
        }
        offset <- stats::rnorm(d, 20, 2)
        sd_j <- if (trend) B_NOISE_SD * 2^(-(offset - B_TREND_REF) / B_TREND_DOUBLING) else rep(B_NOISE_SD, d)
        x <- w %*% t(z) + matrix(stats::rnorm(d * N), d, N) * sd_j + offset
        dimnames(x) <- list(sprintf("%s_f%05d", view, seq_len(d)), samples)
        views[[view]] <- x
    }
    if (missing) {
        set.seed(data_seed + 400000L)
        for (view in names(B_SIZES)) {
            x <- views[[view]]
            p <- 1 - (1 - stats::plogis((B_MISS_CENTRE - x) / B_MISS_WIDTH)) * (1 - B_MISS_FLOOR)
            x[matrix(stats::runif(length(x)), nrow(x)) < p] <- NA_real_
            views[[view]] <- x
        }
    }
    for (view in names(B_SIZES)) views[[view]] <- views[[view]] * B_SCALES[[scales]][[view]]
    list(views = views, z = z, N = N, data_seed = data_seed)
}

fitted_trend <- function(abundance, y) {
    fit <- stats::lowess(abundance, y, f = B_TREND_SPAN)
    stats::approx(fit$x, fit$y, xout = abundance, ties = base::mean, rule = 2)$y
}

# As k151b.R's prepare_views; `zviews` names the data sets that get the per-feature z-score (none, all or a subset).
prepare_b <- function(views, select, cap, eligible, zviews) {
    out <- list()
    for (vi in seq_along(views)) {
        view <- names(views)[[vi]]
        m <- views[[view]]
        z_here <- view %in% zviews
        need <- if (identical(eligible, "half")) ceiling(B_HALF_SHARE * ncol(m)) else 2L
        m <- m[rowSums(!is.na(m)) >= need, , drop = FALSE]
        if (z_here) {
            spread <- apply(m, 1L, stats::sd, na.rm = TRUE)
            m <- m[spread > 0 & !is.na(spread), , drop = FALSE]
        }
        if (!is.na(cap) && nrow(m) > cap) {
            variance <- apply(m, 1L, stats::var, na.rm = TRUE)
            score <- switch(select,
                raw = variance,
                trend_diff = variance - fitted_trend(rowMeans(m, na.rm = TRUE), variance),
                stop("unknown selection ", select))
            if (!all(is.finite(score))) stop("selection score is not finite for every feature in ", view)
            keep <- sort(order(-score, seq_along(score))[seq_len(cap)])
            m <- m[keep, , drop = FALSE]
        }
        if (z_here) {
            centre <- rowMeans(m, na.rm = TRUE)
            spread <- apply(m, 1L, stats::sd, na.rm = TRUE)
            m <- (m - centre) / spread
        }
        rownames(m) <- paste(view, rownames(m), sep = ":")
        out[[view]] <- m
    }
    out
}

# ================================================================== Python side (as k151c.R): versions, BLAS, CPU
py_setup <- function() {
    reticulate::py_run_string(paste(
        "import sys, os, json, importlib.metadata as md",
        "import numpy as np, scipy, sklearn, threadpoolctl, mofapy2",
        "k151d_env = {'python': sys.version.split()[0]}",
        "for p in ('mofapy2', 'numpy', 'scipy', 'scikit-learn', 'threadpoolctl', 'h5py'):",
        "    k151d_env[p] = md.version(p)",
        "k151d_env['environment'] = {k: os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS',",
        "    'MKL_NUM_THREADS', 'OPENBLAS_CORETYPE', 'OPENBLAS_VERBOSE', 'NPY_DISABLE_CPU_FEATURES')}",
        "k151d_env['threadpool_info'] = threadpoolctl.threadpool_info()",
        "from numpy.core._multiarray_umath import __cpu_features__, __cpu_dispatch__, __cpu_baseline__",
        "k151d_env['numpy_cpu_baseline'] = list(__cpu_baseline__)",
        "k151d_env['numpy_cpu_dispatch'] = list(__cpu_dispatch__)",
        "k151d_env['numpy_cpu_features_on'] = sorted(k for k, v in __cpu_features__.items() if v)",
        "k151d_env_json = json.dumps(k151d_env, default=str, indent=1)",
        sep = "\n"))
    main <- reticulate::import_main()
    writeLines(main$k151d_env_json, file.path(out_dir, "python-env.json"))
    env <- main$k151d_env
    list(versions = env[c("python", "mofapy2", "numpy", "scipy", "scikit-learn", "threadpoolctl", "h5py")],
        environment = env$environment,
        blas = lapply(env$threadpool_info, function(x) x[intersect(names(x), c("user_api", "internal_api", "prefix",
            "version", "num_threads", "threading_layer", "architecture", "filepath"))]),
        numpy_cpu_features_on = unlist(env$numpy_cpu_features_on))
}

# ================================================================== MOFA2 training
locate_option <- function(lists, name) {
    where <- names(lists)[vapply(lists, function(l) name %in% names(l), logical(1))]
    if (length(where) != 1L) stop(sprintf("option %s found in %d MOFA2 option lists", name, length(where)))
    where
}

# opts: scale_views, center_groups, weight_views, mode, maxiter
train <- function(views, opts) {
    object <- MOFA2::create_mofa(views)
    lists <- list(data = MOFA2::get_default_data_options(object), model = MOFA2::get_default_model_options(object),
        training = MOFA2::get_default_training_options(object))
    lists$data$scale_views <- opts$scale_views
    lists$data$center_groups <- opts$center_groups
    lists$data$use_float32 <- FALSE
    lists$model$num_factors <- 5L
    lists$training$seed <- TRAIN_SEED
    lists$training$maxiter <- opts$maxiter
    lists$training$convergence_mode <- opts$mode
    lists$training$verbose <- FALSE
    weight_list <- tryCatch(locate_option(lists, "weight_views"), error = function(e) {
        if (isTRUE(opts$weight_views)) stop(conditionMessage(e))
        NA_character_
    })
    if (!is.na(weight_list)) lists[[weight_list]]$weight_views <- isTRUE(opts$weight_views)
    object <- MOFA2::prepare_mofa(object, data_options = lists$data, model_options = lists$model,
        training_options = lists$training)
    outfile <- tempfile("k151d-", fileext = ".hdf5")
    on.exit(unlink(outfile), add = TRUE)
    messages <- character()
    withCallingHandlers(
        MOFA2::run_mofa(object, outfile = outfile, save_data = TRUE, use_basilisk = FALSE),
        message = function(m) {
            messages <<- c(messages, trimws(conditionMessage(m)))
            invokeRestart("muffleMessage")
        },
        warning = function(w) {
            messages <<- c(messages, paste("WARNING:", trimws(conditionMessage(w))))
            invokeRestart("muffleWarning")
        })
    if (any(grepl("Setting use_basilisk = TRUE", messages, fixed = TRUE))) stop("MOFA2 fell back to basilisk")
    model <- withCallingHandlers(MOFA2::load_model(outfile), warning = function(w) {
        messages <<- c(messages, paste("WARNING:", trimws(conditionMessage(w))))
        invokeRestart("muffleWarning")
    })
    samples <- colnames(views[[1L]])
    factors <- MOFA2::get_factors(model)[[1L]][samples, , drop = FALSE]
    weights <- MOFA2::get_weights(model)
    r2 <- MOFA2::get_variance_explained(model)$r2_per_factor[[1L]]
    if (all(rownames(r2) %in% names(views))) r2 <- t(r2)
    stats_slot <- methods::slot(model, "training_stats")
    elbo <- as.numeric(stats_slot$elbo)
    elbo <- elbo[is.finite(elbo)]
    iterations <- length(stats_slot$time)
    list(Z = factors, W = weights, r2 = r2, iterations = iterations,
        stopped_by = if (iterations >= opts$maxiter) "maxiter" else "convergence",
        elbo = if (length(elbo)) utils::tail(elbo, 1L) else NA_real_,
        elbo_rel_change_last5 = if (length(elbo) > 6L) abs((utils::tail(elbo, 1L) - elbo[[length(elbo) - 5L]]) /
            elbo[[1L]]) * 100 else NA_real_,
        warnings = length(grep("^WARNING:", messages)))
}

greedy_match <- function(r) {
    out <- stats::setNames(rep(NA_character_, nrow(r)), rownames(r))
    value <- stats::setNames(rep(0, nrow(r)), rownames(r))
    r[is.na(r)] <- 0
    rows <- rownames(r)
    cols <- colnames(r)
    while (length(rows) && length(cols)) {
        sub <- r[rows, cols, drop = FALSE]
        at <- which(sub == max(sub), arr.ind = TRUE)[1L, ]
        out[[rows[[at[[1L]]]]]] <- cols[[at[[2L]]]]
        value[[rows[[at[[1L]]]]]] <- sub[at[[1L]], at[[2L]]]
        rows <- rows[-at[[1L]]]
        cols <- cols[-at[[2L]]]
    }
    list(factor = out, abs_r = value)
}

# ================================================================== ruv task (Q1): generator and corrections
RUV_SIZES <- c(A = 1000L, B = 500L)
RUV_SHARE <- c(G = 0.20, S = 0.20, U = 0.40)
RUV_NOISE_SD <- 1.0
# Loading weights (+/- weight on the loaded share) per strength candidate; the lifecycle run picks the first candidate,
# in this order, under which the ideal arm recovers S (matched abs r >= 0.9) in at least 9 of 10 lifecycle data sets at
# both N (rule fixed before the lifecycle run; the counted run uses the chosen one, named in its declaration).
RUV_STRENGTH <- list(C1 = c(G = 1, S = 1.5, U = 2), C2 = c(G = 1, S = 2, U = 3))
RUV_RHO <- c(main = 0, conf = 0.6) # correlation of each omic's unwanted factor with the group, exact in the sample
RUV_MODE <- "slow"
RUV_MAXITER <- 10000L
RUV_ARMS <- c("raw", "ideal", "ruv_k1", "ruv_k2", "ruv_k3", "calibrated", "probe_fail")
COUNT_SEED <- 20261118L # RUV_FACTOR_COUNT_SEED in R/func_general_ruv_factor_count.R
COUNT_PERMUTATIONS <- 200L # decision M30

generate_ruv <- function(N, seed_index, scenario, strength) {
    data_seed <- 154000L + seed_index
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(data_seed)
    group <- rep(c(-1, 1), each = N / 2L)
    within <- function(x) {
        x <- x - stats::ave(x, group)
        x / stats::sd(x)
    }
    s <- within(stats::rnorm(N))
    e_a <- within(stats::rnorm(N))
    e_b <- within(stats::rnorm(N))
    rho <- RUV_RHO[[scenario]]
    g_std <- group / stats::sd(group)
    latent <- cbind(G = group, S = s, UA = rho * g_std + sqrt(1 - rho^2) * e_a,
        UB = rho * g_std + sqrt(1 - rho^2) * e_b)
    samples <- sprintf("S%02d", seq_len(N))
    rownames(latent) <- samples
    weights <- RUV_STRENGTH[[strength]]
    views <- list()
    loadings <- list()
    for (view in names(RUV_SIZES)) {
        d <- RUV_SIZES[[view]]
        u <- if (view == "A") "UA" else "UB"
        w <- matrix(0, d, 3L, dimnames = list(NULL, c("G", "S", u)))
        for (k in c("G", "S", u)) {
            kind <- if (k %in% c("UA", "UB")) "U" else k
            idx <- sample.int(d, round(RUV_SHARE[[kind]] * d))
            w[idx, k] <- sample(c(-1, 1), length(idx), replace = TRUE) * weights[[kind]]
        }
        offset <- stats::rnorm(d, 20, 2)
        x <- w %*% t(latent[, colnames(w)]) + matrix(stats::rnorm(d * N, 0, RUV_NOISE_SD), d, N) + offset
        dimnames(x) <- list(sprintf("%s_f%04d", view, seq_len(d)), samples)
        rownames(w) <- rownames(x)
        views[[view]] <- x
        loadings[[view]] <- w
    }
    list(views = views, latent = latent, loadings = loadings, group = group, N = N, data_seed = data_seed)
}

ruv_replicate_matrix <- function(group) {
    # as getRuvIIIReplicateMatrixHelper: one 0/1 column per group, samples in rows, no intercept
    m <- stats::model.matrix(~ g - 1, data = data.frame(g = factor(group)))
    attr(m, "assign") <- NULL
    attr(m, "contrasts") <- NULL
    m
}

ruv_correct <- function(x, group, k) {
    Y <- t(x) # samples by features, as the package passes it
    M <- ruv_replicate_matrix(group)
    rownames(M) <- rownames(Y)
    out <- RUVIIIC::RUVIII_C_Varying(k = k, Y = Y, M = M, toCorrect = colnames(Y), potentialControls = colnames(Y),
        progress = FALSE)
    out <- as.matrix(out)
    if (identical(sort(rownames(out)), sort(colnames(x))) && identical(sort(colnames(out)), sort(rownames(x)))) {
        out <- t(out)
    }
    out <- out[rownames(x), colnames(x), drop = FALSE]
    if (any(!is.finite(out))) stop(sprintf("RUV-III k = %d returned non-finite values", k))
    out
}

ruv_reference <- function(x, group, k) {
    Y <- t(x)
    t(ruv::RUVIII(Y = Y, M = ruv_replicate_matrix(group), ctl = rep(TRUE, ncol(Y)), k = k))
}

# The package's estimated count (estimateRuvFactorCount, groups removed, no technical replicates): sva::num.sv on the
# features observed in every sample under the fixed RNG kind, and the limits of computeOmicsRuvAssayFactors (k below the
# controls; a residual degree of freedom beside the groups; samples minus replicate sets). k = 0 or a refused k: the
# view is not corrected (the package analyses such an assay without RUV).
calibrated_k <- function(x, group) {
    dat <- unname(x[rowSums(!is.finite(x)) == 0L, , drop = FALSE])
    mod <- stats::model.matrix(~ group, data.frame(group = as.character(group)))
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(COUNT_SEED)
    k <- as.integer(sva::num.sv(dat, mod, method = "be", B = COUNT_PERMUTATIONS, seed = COUNT_SEED))
    n <- ncol(x)
    groups <- length(unique(group))
    status <- if (k == 0L) "none" else if (k >= nrow(dat) || k > n - groups - 1L || k > n - groups) "refused" else "ok"
    list(k = k, status = status, features_used = nrow(dat))
}

# Joint least-squares coefficients of each feature on the planted G, S and the view's U; the share left of a factor's
# signal = mean abs coefficient over the features loaded on it, corrected over uncorrected.
coefficients_on <- function(x, design) {
    X <- cbind(1, design)
    t(solve(crossprod(X), crossprod(X, t(x))))[, -1L, drop = FALSE]
}

# ================================================================== run
result <- list(task = task, arguments = arguments[-1L])
result$cpu <- tryCatch(system("lscpu", intern = TRUE), error = function(e) conditionMessage(e))
result$cpu_model <- trimws(sub("^Model name:", "", grep("^Model name:", result$cpu, value = TRUE)[1L]))
result$cpu_vendor <- trimws(sub("^Vendor ID:", "", grep("^Vendor ID:", result$cpu, value = TRUE)[1L]))
flags <- grep("^Flags:", result$cpu, value = TRUE)
result$avx512f <- length(flags) == 1L && grepl("\\bavx512f\\b", flags)
result$nproc <- parallel::detectCores()
rows <- list()
fits <- list()
failed <- 0L

status <- tryCatch({
    suppressPackageStartupMessages(library(MOFA2))
    result$python <- py_setup()
    result$versions <- c(list(R = as.character(getRversion()), MOFA2 = as.character(utils::packageVersion("MOFA2")),
        reticulate = as.character(utils::packageVersion("reticulate"))), result$python$versions)
    cat(sprintf("CPU %s [%s] (avx512f %s); Python env: %s\n", result$cpu_model, result$cpu_vendor, result$avx512f,
        jsonlite::toJSON(result$python[c("environment", "blas")], auto_unbox = TRUE)))

    run_cell <- function(key, base, expr_fn, keep_fit = FALSE) {
        clock <- Sys.time()
        outcome <- tryCatch(expr_fn(), error = function(e) {
            failed <<- failed + 1L
            message("CELL FAILED ", key, ": ", conditionMessage(e))
            structure(list(message = conditionMessage(e)), class = "cell_error")
        })
        seconds <- round(as.numeric(difftime(Sys.time(), clock, units = "secs")), 2)
        if (inherits(outcome, "cell_error")) {
            row <- cbind(base, data.frame(status = "failed", error = outcome$message, stringsAsFactors = FALSE))
        } else {
            if (keep_fit) fits[[key]] <<- outcome$fit
            row <- cbind(base, data.frame(status = "ok", error = NA_character_, stringsAsFactors = FALSE), outcome$row)
        }
        row$seconds <- seconds
        rows[[key]] <<- row
        cat(sprintf("CELL %s: %s; %s s\n", key, row$status[[1L]], seconds))
    }

    if (task == "ruv") {
        N <- as.integer(arguments[[3L]])
        scenarios <- split_list(arguments[[4L]])
        strength <- arguments[[5L]]
        indices <- as.integer(split_list(arguments[[6L]]))
        arms <- split_list(arguments[[7L]])
        stopifnot(N %in% c(20L, 60L), all(scenarios %in% names(RUV_RHO)), strength %in% names(RUV_STRENGTH),
            !anyNA(indices), all(arms %in% RUV_ARMS))
        suppressPackageStartupMessages({
            library(RUVIIIC)
            library(ruv)
            library(sva)
        })
        result$versions$RUVIIIC <- as.character(utils::packageVersion("RUVIIIC"))
        result$versions$ruv <- as.character(utils::packageVersion("ruv"))
        result$versions$sva <- as.character(utils::packageVersion("sva"))
        if (!identical(result$versions$sva, "3.60.0")) stop("sva is not 3.60.0 (the version the package's count follows)")
        checks <- list()
        for (scenario in scenarios) for (index in indices) {
            data <- generate_ruv(N, index, scenario, strength)
            counts <- NULL
            for (arm_name in arms) {
                key <- sprintf("N%d|%s|%s|%d|%s", N, scenario, strength, index, arm_name)
                base <- data.frame(N = N, scenario = scenario, strength = strength, seed_index = index,
                    data_seed = data$data_seed, arm = arm_name, mode = RUV_MODE, train_seed = TRAIN_SEED,
                    stringsAsFactors = FALSE)
                run_cell(key, base, function() {
                    if (arm_name == "probe_fail") stop("probe: deliberate failure of the lifecycle cell")
                    corrected <- list()
                    k_used <- c(A = NA_integer_, B = NA_integer_)
                    k_status <- c(A = NA_character_, B = NA_character_)
                    for (view in names(RUV_SIZES)) {
                        x <- data$views[[view]]
                        u <- if (view == "A") "UA" else "UB"
                        if (arm_name == "calibrated") {
                            if (is.null(counts)) counts <<- lapply(data$views, calibrated_k, group = data$group)
                            cnt <- counts[[view]]
                            k_used[[view]] <- cnt$k
                            k_status[[view]] <- cnt$status
                        }
                        corrected[[view]] <- switch(arm_name,
                            raw = x,
                            ideal = x - outer(data$loadings[[view]][, u], data$latent[, u]),
                            calibrated = if (identical(k_status[[view]], "ok")) ruv_correct(x, data$group, k_used[[view]]) else x,
                            ruv_correct(x, data$group, as.integer(sub("ruv_k", "", arm_name))))
                        if (startsWith(arm_name, "ruv_k") && index == indices[[1L]] && scenario == scenarios[[1L]]) {
                            k <- as.integer(sub("ruv_k", "", arm_name))
                            checks[[length(checks) + 1L]] <<- data.frame(N = N, scenario = scenario, seed_index = index,
                                view = view, k = k,
                                max_abs_diff_vs_ruv_RUVIII = max(abs(corrected[[view]] - ruv_reference(x, data$group, k))),
                                stringsAsFactors = FALSE)
                        }
                    }
                    keep <- list()
                    for (view in names(RUV_SIZES)) {
                        u <- if (view == "A") "UA" else "UB"
                        design <- data$latent[, c("G", "S", u)]
                        before <- coefficients_on(data$views[[view]], design)
                        after <- coefficients_on(corrected[[view]], design)
                        w <- data$loadings[[view]]
                        for (f in c("G", "S", u)) {
                            loaded <- w[, f] != 0
                            keep[[sprintf("retained_%s_%s", if (f == u) "U" else f, view)]] <-
                                mean(abs(after[loaded, f])) / mean(abs(before[loaded, f]))
                        }
                    }
                    views <- lapply(names(RUV_SIZES), function(view) {
                        m <- corrected[[view]]
                        rownames(m) <- paste(view, rownames(m), sep = ":")
                        m
                    })
                    names(views) <- names(RUV_SIZES)
                    fit <- train(views, list(scale_views = TRUE, center_groups = TRUE, weight_views = FALSE,
                        mode = RUV_MODE, maxiter = RUV_MAXITER))
                    r <- abs(stats::cor(data$latent, fit$Z))
                    match <- greedy_match(r)
                    gf <- match$factor[["G"]]
                    if (is.na(gf)) gf <- colnames(fit$Z)[[which.max(r["G", ])]]
                    gz <- fit$Z[, gf]
                    fit_lm <- stats::lm(gz ~ factor(data$group))
                    an <- stats::anova(fit_lm)
                    sf <- match$factor[["S"]]
                    row <- data.frame(iterations = fit$iterations, stopped_by = fit$stopped_by, elbo = fit$elbo,
                        elbo_rel_change_last5 = fit$elbo_rel_change_last5, learned_factors = ncol(fit$Z),
                        r_G = match$abs_r[["G"]], r_S = match$abs_r[["S"]], r_UA = match$abs_r[["UA"]],
                        r_UB = match$abs_r[["UB"]], best_r_S = max(r["S", ]), best_r_G = max(r["G", ]),
                        group_factor = gf, group_factor_ve_A = unname(fit$r2[gf, "A"]),
                        group_factor_ve_B = unname(fit$r2[gf, "B"]), group_factor_R2 = summary(fit_lm)$r.squared,
                        group_factor_p = an[["Pr(>F)"]][[1L]], s_factor = sf,
                        s_factor_ve_A = if (is.na(sf)) NA_real_ else unname(fit$r2[sf, "A"]),
                        s_factor_ve_B = if (is.na(sf)) NA_real_ else unname(fit$r2[sf, "B"]),
                        total_ve_A = sum(fit$r2[, "A"]), total_ve_B = sum(fit$r2[, "B"]),
                        as.data.frame(keep), k_A = k_used[["A"]], k_B = k_used[["B"]], k_status_A = k_status[["A"]],
                        k_status_B = k_status[["B"]], cor_UA_group = stats::cor(data$latent[, "UA"], data$group),
                        mofa_warnings = fit$warnings, stringsAsFactors = FALSE)
                    list(fit = NULL, row = row)
                })
            }
        }
        if (length(checks)) utils::write.table(do.call(rbind, checks), file.path(out_dir, "ruv-check.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
    } else if (task == "z") {
        N <- as.integer(arguments[[3L]])
        index <- as.integer(arguments[[4L]])
        cells <- split_list(arguments[[5L]])
        stopifnot(N %in% c(20L, 60L), !is.na(index))
        # cell code: m<0|1>_<scales>_z<none|all|small|large>_sv<0|1>, or probe_fail
        parse_cell <- function(cell) {
            p <- strsplit(cell, "_", fixed = TRUE)[[1L]]
            stopifnot(length(p) == 4L, p[[1L]] %in% c("m0", "m1"), p[[2L]] %in% names(B_SCALES),
                p[[3L]] %in% c("znone", "zall", "zsmall", "zlarge"), p[[4L]] %in% c("sv0", "sv1"))
            list(missing = p[[1L]] == "m1", scales = p[[2L]], z = sub("^z", "", p[[3L]]), sv = p[[4L]] == "sv1")
        }
        Z_VIEWS <- list(none = character(), all = names(B_SIZES), small = c("D3", "D4"), large = c("D1", "D2"))
        data_sets <- list()
        result$missing_share <- list()
        for (cell in cells) {
            base <- data.frame(N = N, seed_index = index, data_seed = 151000L + index, cell = cell,
                stringsAsFactors = FALSE)
            run_cell(cell, base, function() {
                if (cell == "probe_fail") stop("probe: deliberate failure of the lifecycle cell")
                spec <- parse_cell(cell)
                dkey <- sprintf("m%d_%s", as.integer(spec$missing), spec$scales)
                if (is.null(data_sets[[dkey]])) {
                    data_sets[[dkey]] <<- generate_b(N, index, trend = TRUE, missing = spec$missing, scales = spec$scales)
                    result$missing_share[[dkey]] <<- lapply(data_sets[[dkey]]$views, function(v) round(mean(is.na(v)), 4))
                }
                data <- data_sets[[dkey]]
                views <- prepare_b(data$views, "trend_diff", 1000L, "half", Z_VIEWS[[spec$z]])
                kept <- vapply(views, nrow, integer(1))
                fit <- train(views, list(scale_views = spec$sv, center_groups = TRUE, weight_views = FALSE,
                    mode = "slow", maxiter = 10000L))
                r <- abs(stats::cor(data$z, fit$Z))
                match <- greedy_match(r)
                row <- data.frame(missing = spec$missing, scales = spec$scales, z = spec$z, scale_views = spec$sv,
                    iterations = fit$iterations, stopped_by = fit$stopped_by, elbo = fit$elbo,
                    elbo_rel_change_last5 = fit$elbo_rel_change_last5, learned_factors = ncol(fit$Z),
                    r_F1 = match$abs_r[["F1"]], r_F2 = match$abs_r[["F2"]], r_F3 = match$abs_r[["F3"]],
                    f_F1 = match$factor[["F1"]], f_F2 = match$factor[["F2"]], f_F3 = match$factor[["F3"]],
                    stringsAsFactors = FALSE)
                for (f in c("F1", "F2", "F3")) for (view in names(B_SIZES)) {
                    row[[sprintf("r2_%s_%s", f, view)]] <- if (is.na(match$factor[[f]])) NA_real_ else
                        unname(fit$r2[match$factor[[f]], view])
                }
                for (view in names(B_SIZES)) row[[paste0("kept_", view)]] <- unname(kept[[view]])
                row$mofa_warnings <- fit$warnings
                list(fit = NULL, row = row)
            })
        }
    } else {
        condition <- arguments[[3L]]
        cases <- split_list(arguments[[4L]])
        result$condition <- condition
        # case codes: fast:<k151c arm> (n20weak@7, fast, maxiter 1,000, as k151c's pin task); hard_prev10_none (k151b
        # prev@10, no selection, scale_views off); hard_full23_default (k151b full@23, the decided selection,
        # scale_views on); both slow, maxiter 10,000; probe_fail.
        for (case in cases) {
            base <- data.frame(condition = condition, case = case, stringsAsFactors = FALSE)
            run_cell(paste(condition, case, sep = "|"), base, function() {
                if (case == "probe_fail") stop("probe: deliberate failure of the lifecycle cell")
                if (startsWith(case, "fast:")) {
                    arm <- K151_ARMS[[sub("^fast:", "", case)]]
                    if (is.null(arm)) stop("unknown arm in ", case)
                    data <- generate_k151("n20weak", 7L)
                    views <- prepare_k151(data$views, arm)
                    opts <- c(arm[c("scale_views", "center_groups", "weight_views")], list(mode = "fast", maxiter = 1000L))
                } else if (case == "hard_prev10_none") {
                    data <- generate_b(20L, 10L, trend = FALSE, missing = FALSE, scales = "none")
                    views <- prepare_b(data$views, "none", NA_integer_, "two", character())
                    opts <- list(scale_views = FALSE, center_groups = TRUE, weight_views = FALSE, mode = "slow",
                        maxiter = 10000L)
                } else if (case == "hard_full23_default") {
                    data <- generate_b(20L, 23L, trend = TRUE, missing = TRUE, scales = "orig")
                    views <- prepare_b(data$views, "trend_diff", 1000L, "half", character())
                    opts <- list(scale_views = TRUE, center_groups = TRUE, weight_views = FALSE, mode = "slow",
                        maxiter = 10000L)
                } else stop("unknown case ", case)
                fit <- train(views, opts)
                match <- greedy_match(abs(stats::cor(data$z, fit$Z)))
                row <- data.frame(data_seed = data$data_seed, iterations = fit$iterations, stopped_by = fit$stopped_by,
                    elbo = fit$elbo, learned_factors = ncol(fit$Z), planted_F1 = match$abs_r[["F1"]],
                    planted_F2 = match$abs_r[["F2"]], planted_F3 = match$abs_r[["F3"]], mofa_warnings = fit$warnings,
                    stringsAsFactors = FALSE)
                list(fit = fit[c("Z", "W", "r2", "iterations", "elbo")], row = row)
            }, keep_fit = TRUE)
        }
    }
    if (failed > 0L) 1L else 0L
}, error = function(e) {
    result$error <<- conditionMessage(e)
    message("JOB ERROR: ", conditionMessage(e))
    1L
})

bind_rows <- function(rows) {
    columns <- unique(unlist(lapply(rows, names)))
    do.call(rbind, lapply(rows, function(r) {
        for (col in setdiff(columns, names(r))) r[[col]] <- NA
        r[columns]
    }))
}
if (length(rows)) {
    table <- bind_rows(rows)
    if ("error" %in% names(table)) table$error <- gsub("[\t\r\n]+", " ", table$error)
    utils::write.table(table, file.path(out_dir, "rows.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
result$failed_cells <- failed
result$status <- status
result$seconds_total <- round(as.numeric(difftime(Sys.time(), started, units = "secs")), 2)
if (task == "pin") saveRDS(list(job = result, fits = fits), file.path(out_dir, "fits.rds"))
jsonlite::write_json(result, file.path(out_dir, "job.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE,
    null = "null", na = "null")
cat(sprintf("JOB %s: %d cells, %d failed, %s s, CPU %s\n", task, length(rows), failed, result$seconds_total,
    result$cpu_model))
quit(save = "no", status = status)
