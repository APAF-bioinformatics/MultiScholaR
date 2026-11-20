# KIT-151 c evidence (branch run-mofa-evidence-k151c): three small questions about the MOFA2 model, one task per job.
# Declarations: docs/qa/workflow-execution/kit/kit-151c/{seed,pin,ruv}-declaration.md.
#   seed  (SCI-06) training seeds on fixed data sets, one job, one processor: do seeds give different models?
#   pin   (SCI-05) the divergent case (N = 20 weak, data seed 151007, fast), run under one environment condition
#         (the workflow sets the condition's variables before R, and so Python, starts).
#   ruv   (SCI-04, INT-004) planted shared within-group factor; inputs uncorrected, oracle-corrected and corrected by
#         RUV-III (groups as replicate sets, every feature a control; RUVIIIC::RUVIII_C_Varying, as the package calls
#         it) with k = 1, 2, 3.
# MOFA2 is called directly as in k151-followup.R (same view preparation and options; data sets of the seed and pin
# tasks come from the same generator and so equal the k151 runs' data).
#
# Usage:
#   Rscript --vanilla k151c.R <out> seed <datasets: scenario@index,...> <arms,...> <modes,...> <train seeds,...>
#   Rscript --vanilla k151c.R <out> pin  <condition label> <arms,...>
#   Rscript --vanilla k151c.R <out> ruv  <N> <data seed indices,...> <arms,...>
# Writes rows.tsv (one row per cell), fits.rds (factors, weights, start values per cell, the job record),
# job.json, python-env.json and, for seed, source-excerpts.txt. Exit status 1 when any cell failed.

arguments <- commandArgs(trailingOnly = TRUE)
out_dir <- arguments[[1L]]
task <- arguments[[2L]]
stopifnot(task %in% c("seed", "pin", "ruv"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()
split_list <- function(x) strsplit(x, ",", fixed = TRUE)[[1L]]

# ------------------------------------------------------------------ k151 generator (copied from k151-followup.R)
SIZES <- c(D1 = 20000L, D2 = 2000L, D3 = 200L, D4 = 100L)
ACTIVE <- list(D1 = c("F1", "F3"), D2 = "F1", D3 = c("F1", "F2"), D4 = c("F1", "F2"))
LOADED_SHARE <- 0.10
EFFECT <- 1.0
NOISE_SD <- 1.5
F1_NOISE_SD <- 0.5
MAXITER <- 1000L
PIN_TRAIN_SEED <- 20261119L

generate_k151 <- function(scenario, seed_index) {
    N <- if (startsWith(scenario, "n60")) 60L else 20L
    F2_EFFECT <- if (grepl("strong$", scenario)) 1.5 else 1.0
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(151000L + seed_index)
    group <- rep(c(-1, 1), each = N / 2L)
    f1 <- group + stats::rnorm(N, 0, F1_NOISE_SD)
    f2 <- stats::rnorm(N)
    f3 <- stats::rnorm(N)
    z <- scale(cbind(F1 = f1, F2 = f2, F3 = f3))
    samples <- sprintf("S%02d", seq_len(N))
    rownames(z) <- samples
    views <- list()
    for (view in names(SIZES)) {
        d <- SIZES[[view]]
        w <- matrix(0, d, 3L, dimnames = list(NULL, colnames(z)))
        for (k in ACTIVE[[view]]) {
            loaded <- sample.int(d, round(LOADED_SHARE * d))
            w[loaded, k] <- sample(c(-1, 1), length(loaded), replace = TRUE) * (if (k == "F2") F2_EFFECT else EFFECT)
        }
        offset <- stats::rnorm(d, 20, 2)
        x <- w %*% t(z) + matrix(stats::rnorm(d * N, 0, NOISE_SD), d, N) + offset
        dimnames(x) <- list(sprintf("%s_f%05d", view, seq_len(d)), samples)
        views[[view]] <- x
    }
    list(views = views, z = z, N = N, data_seed = 151000L + seed_index)
}

base_arm <- function(cap = NA_integer_, zscore = FALSE, scale_views = TRUE, weight_views = FALSE, fail = FALSE) {
    list(cap = cap, zscore = zscore, scale_views = scale_views, center_groups = !zscore, weight_views = weight_views,
        fail = fail)
}
K151_ARMS <- list(
    none = base_arm(scale_views = FALSE),
    scale_views = base_arm(),
    cap1000 = base_arm(cap = 1000L),
    cap2000 = base_arm(cap = 2000L),
    cap5000 = base_arm(cap = 5000L),
    cap10000 = base_arm(cap = 10000L),
    weight_views = base_arm(weight_views = TRUE),
    zscore = base_arm(zscore = TRUE, scale_views = FALSE),
    probe_fail = base_arm(fail = TRUE)) # lifecycle only

# view preparation (copied from k151-followup.R)
prepare_views <- function(views, arm) {
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

# ------------------------------------------------------------------ Python side: versions, BLAS, CPU features, hooks
py_setup <- function() {
    reticulate::py_run_string(paste(
        "import sys, os, json, io, contextlib, importlib.metadata as md",
        "import numpy as np, scipy, scipy.linalg, scipy.special, sklearn, sklearn.decomposition as skd",
        "import threadpoolctl, mofapy2",
        "import mofapy2.build_model.init_model as k151c_im",
        "k151c_env = {'python': sys.version.split()[0]}",
        "for p in ('mofapy2', 'numpy', 'scipy', 'scikit-learn', 'threadpoolctl', 'h5py'):",
        "    k151c_env[p] = md.version(p)",
        "k151c_env['environment'] = {k: os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS',",
        "    'MKL_NUM_THREADS', 'OPENBLAS_CORETYPE', 'OPENBLAS_VERBOSE', 'NPY_DISABLE_CPU_FEATURES')}",
        "k151c_env['threadpool_info'] = threadpoolctl.threadpool_info()",
        "try:",
        "    k151c_env['numpy_config'] = np.show_config(mode='dicts')",
        "except Exception as e:",
        "    k151c_env['numpy_config'] = 'show_config failed: ' + repr(e)",
        "from numpy.core._multiarray_umath import __cpu_features__, __cpu_dispatch__, __cpu_baseline__",
        "k151c_env['numpy_cpu_baseline'] = list(__cpu_baseline__)",
        "k151c_env['numpy_cpu_dispatch'] = list(__cpu_dispatch__)",
        "k151c_env['numpy_cpu_features_on'] = sorted(k for k, v in __cpu_features__.items() if v)",
        "k151c_capture = {'initZ': [], 'pca_solver': [], 'errors': []}",
        "if not getattr(k151c_im.initModel, '_k151c_wrapped', False):",
        "    _k151c_initZ = k151c_im.initModel.initZ",
        "    def _k151c_initZ_wrapped(self, *a, **kw):",
        "        out = _k151c_initZ(self, *a, **kw)",
        "        try:",
        "            k151c_capture['initZ'].append(np.array(self.nodes['Z'].getExpectation(), dtype=float, copy=True))",
        "        except Exception as e:",
        "            k151c_capture['errors'].append('initZ capture: ' + repr(e))",
        "        return out",
        "    k151c_im.initModel.initZ = _k151c_initZ_wrapped",
        "    _k151c_pca_fit = skd.PCA._fit",
        "    def _k151c_pca_fit_wrapped(self, X):",
        "        out = _k151c_pca_fit(self, X)",
        "        k151c_capture['pca_solver'].append(str(getattr(self, '_fit_svd_solver', 'unknown')) + ' ' + str(X.shape))",
        "        return out",
        "    skd.PCA._fit = _k151c_pca_fit_wrapped",
        "    k151c_im.initModel._k151c_wrapped = True",
        sep = "\n"))
    main <- reticulate::import_main()
    reticulate::py_run_string("k151c_env_json = json.dumps(k151c_env, default=str, indent=1)")
    writeLines(main$k151c_env_json, file.path(out_dir, "python-env.json"))
    env <- main$k151c_env
    list(versions = env[c("python", "mofapy2", "numpy", "scipy", "scikit-learn", "threadpoolctl", "h5py")],
        environment = env$environment,
        blas = lapply(env$threadpool_info, function(x) x[intersect(names(x), c("user_api", "internal_api", "prefix",
            "version", "num_threads", "threading_layer", "architecture", "filepath"))]),
        numpy_cpu_features_on = unlist(env$numpy_cpu_features_on))
}

# The installed mofapy2 (and scikit-learn PCA) lines that set the start and use the seed: printed and saved.
source_excerpts <- function() {
    reticulate::py_run_string(paste(
        "import pathlib, re, mofapy2, sklearn",
        "k151c_src = []",
        "def k151c_show(path, pattern, before, after, root):",
        "    lines = path.read_text(errors='replace').splitlines()",
        "    for i, line in enumerate(lines, 1):",
        "        if re.search(pattern, line):",
        "            k151c_src.append('---- ' + str(path.relative_to(root)) + ':' + str(i))",
        "            for j in range(max(1, i - before), min(len(lines), i + after) + 1):",
        "                k151c_src.append(f'{path.relative_to(root)}:{j}: {lines[j - 1]}')",
        "m = pathlib.Path(mofapy2.__file__).parent",
        "s = pathlib.Path(sklearn.__file__).parent",
        "k151c_src.append('mofapy2 installed at ' + str(m) + '; scikit-learn at ' + str(s))",
        "k151c_show(m / 'build_model' / 'build_model.py', r'init_model\\.init|seed|initModel\\(', 2, 6, m.parent)",
        "k151c_show(m / 'build_model' / 'init_model.py', r'np\\.random\\.seed|def initZ\\b|qmean == \"pca\"|def initW\\b|def initSW\\b', 2, 16, m.parent)",
        "k151c_show(m / 'run' / 'entry_point.py', r'seed', 1, 2, m.parent)",
        "k151c_show(m / 'core' / 'BayesNet.py', r'random|seed|shuffle', 1, 1, m.parent)",
        "k151c_show(s / 'decomposition' / '_pca.py', r'_fit_svd_solver = ', 3, 3, s.parent)",
        "k151c_show(s / 'decomposition' / '_pca.py', r'def _fit_truncated|random_state = check_random_state|randomized_svd\\(', 1, 10, s.parent)",
        sep = "\n"))
    lines <- unlist(reticulate::import_main()$k151c_src)
    writeLines(lines, file.path(out_dir, "source-excerpts.txt"))
    cat("SOURCE EXCERPTS (installed mofapy2 and scikit-learn)\n")
    cat(lines, sep = "\n")
    cat("END SOURCE EXCERPTS\n")
}

# ------------------------------------------------------------------ MOFA2 training (as k151-followup.R, seed a parameter)
locate_option <- function(lists, name) {
    where <- names(lists)[vapply(lists, function(l) name %in% names(l), logical(1))]
    if (length(where) != 1L) stop(sprintf("option %s found in %d MOFA2 option lists", name, length(where)))
    where
}

train <- function(views, arm, mode, seed) {
    if (isTRUE(arm$fail)) stop("probe: deliberate failure of the lifecycle arm")
    reticulate::py_run_string("k151c_capture['initZ'].clear(); k151c_capture['pca_solver'].clear()")
    object <- MOFA2::create_mofa(views)
    lists <- list(data = MOFA2::get_default_data_options(object), model = MOFA2::get_default_model_options(object),
        training = MOFA2::get_default_training_options(object))
    lists$data$scale_views <- arm$scale_views
    lists$data$center_groups <- arm$center_groups
    lists$data$use_float32 <- FALSE
    lists$model$num_factors <- 5L
    lists$training$seed <- seed
    lists$training$maxiter <- MAXITER
    lists$training$convergence_mode <- mode
    lists$training$verbose <- FALSE
    weight_list <- tryCatch(locate_option(lists, "weight_views"), error = function(e) {
        if (arm$weight_views) stop(conditionMessage(e))
        NA_character_
    })
    if (!is.na(weight_list)) lists[[weight_list]]$weight_views <- arm$weight_views
    object <- MOFA2::prepare_mofa(object, data_options = lists$data, model_options = lists$model,
        training_options = lists$training)
    outfile <- tempfile("k151c-", fileext = ".hdf5")
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
    capture <- reticulate::import_main()$k151c_capture
    samples <- colnames(views[[1L]])
    factors <- MOFA2::get_factors(model)[[1L]][samples, , drop = FALSE]
    weights <- MOFA2::get_weights(model)
    r2 <- MOFA2::get_variance_explained(model)$r2_per_factor[[1L]]
    if (all(rownames(r2) %in% names(views))) r2 <- t(r2)
    stats_slot <- methods::slot(model, "training_stats")
    elbo <- as.numeric(stats_slot$elbo)
    elbo <- elbo[is.finite(elbo)]
    init <- if (length(capture$initZ) == 1L) capture$initZ[[1L]] else NULL
    if (!is.null(init)) rownames(init) <- samples
    list(Z = factors, W = weights, r2 = r2, iterations = length(stats_slot$time),
        elbo = if (length(elbo)) utils::tail(elbo, 1L) else NA_real_, init_Z = init,
        pca_solver = paste(unlist(capture$pca_solver), collapse = "; "),
        capture_errors = paste(unlist(capture$errors), collapse = "; "),
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

# ------------------------------------------------------------------ RUV task: generator and corrections (declared)
RUV_SIZES <- c(A = 1000L, B = 500L)
RUV_SHARE <- c(G = 0.20, S = 0.20, U = 0.40)
RUV_NOISE_SD <- 1.0
RUV_TRAIN_SEED <- 20261119L
RUV_MODE <- "slow"
RUV_ARMS <- c("raw", "oracle", "ruv_k1", "ruv_k2", "ruv_k3", "probe_fail")

generate_ruv <- function(N, seed_index) {
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(153000L + seed_index)
    group <- rep(c(-1, 1), each = N / 2L)
    within <- function(x) {
        x <- x - stats::ave(x, group)
        x / stats::sd(x)
    }
    latent <- cbind(G = group, S = within(stats::rnorm(N)), UA = within(stats::rnorm(N)), UB = within(stats::rnorm(N)))
    samples <- sprintf("S%02d", seq_len(N))
    rownames(latent) <- samples
    views <- list()
    loaded <- list()
    for (view in names(RUV_SIZES)) {
        d <- RUV_SIZES[[view]]
        u <- if (view == "A") "UA" else "UB"
        w <- matrix(0, d, 3L, dimnames = list(NULL, c("G", "S", u)))
        for (k in c("G", "S", u)) {
            share <- RUV_SHARE[[if (k %in% c("UA", "UB")) "U" else k]]
            idx <- sample.int(d, round(share * d))
            w[idx, k] <- sample(c(-1, 1), length(idx), replace = TRUE)
        }
        offset <- stats::rnorm(d, 20, 2)
        x <- w %*% t(latent[, colnames(w)]) + matrix(stats::rnorm(d * N, 0, RUV_NOISE_SD), d, N) + offset
        dimnames(x) <- list(sprintf("%s_f%04d", view, seq_len(d)), samples)
        views[[view]] <- x
        loaded[[view]] <- w != 0
        colnames(loaded[[view]]) <- c("G", "S", "U")
    }
    list(views = views, latent = latent, loaded = loaded, group = group, N = N, data_seed = 153000L + seed_index)
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
    M <- ruv_replicate_matrix(group)
    t(ruv::RUVIII(Y = Y, M = M, ctl = rep(TRUE, ncol(Y)), k = k))
}

oracle_correct <- function(x, u) {
    b <- as.vector(x %*% (u - mean(u))) / sum((u - mean(u))^2)
    x - outer(b, u - mean(u))
}

retained <- function(corrected, raw, loaded_rows, covariate) {
    slope <- function(m) abs(as.vector(m[loaded_rows, , drop = FALSE] %*% (covariate - mean(covariate))))
    mean(slope(corrected)) / mean(slope(raw))
}

# ------------------------------------------------------------------ run
result <- list(task = task, arguments = arguments[-1L])
result$cpu <- tryCatch(system("lscpu", intern = TRUE), error = function(e) conditionMessage(e))
result$cpu_model <- trimws(sub("^Model name:", "", grep("^Model name:", result$cpu, value = TRUE)[1L]))
flags <- grep("^Flags:", result$cpu, value = TRUE)
result$avx512f <- length(flags) == 1L && grepl("\\bavx512f\\b", flags)
result$nproc <- parallel::detectCores()
rows <- list()
fits <- list()
failed <- 0L

cell_row <- function(base, status, error = NA_character_) {
    cbind(base, data.frame(status = status, error = error, stringsAsFactors = FALSE))
}

status <- tryCatch({
    suppressPackageStartupMessages(library(MOFA2))
    result$python <- py_setup()
    result$versions <- c(list(R = as.character(getRversion()), MOFA2 = as.character(utils::packageVersion("MOFA2")),
        reticulate = as.character(utils::packageVersion("reticulate"))), result$python$versions)
    cat(sprintf("CPU %s (avx512f %s); Python env: %s\n", result$cpu_model, result$avx512f,
        jsonlite::toJSON(result$python[c("environment", "blas")], auto_unbox = TRUE)))
    if (task == "seed") source_excerpts()

    run_cell <- function(key, base, expr_fn) {
        clock <- Sys.time()
        outcome <- tryCatch(expr_fn(), error = function(e) {
            failed <<- failed + 1L
            message("CELL FAILED ", key, ": ", conditionMessage(e))
            structure(list(message = conditionMessage(e)), class = "cell_error")
        })
        seconds <- round(as.numeric(difftime(Sys.time(), clock, units = "secs")), 2)
        if (inherits(outcome, "cell_error")) {
            row <- cell_row(base, "failed", outcome$message)
        } else {
            fits[[key]] <<- outcome$fit
            row <- cbind(cell_row(base, "ok"), outcome$row)
        }
        row$seconds <- seconds
        rows[[key]] <<- row
        cat(sprintf("CELL %s: %s; %s s\n", key, row$status, seconds))
    }

    if (task == "seed") {
        datasets <- split_list(arguments[[3L]])
        arms <- split_list(arguments[[4L]])
        modes <- split_list(arguments[[5L]])
        seeds <- as.integer(split_list(arguments[[6L]]))
        stopifnot(all(arms %in% names(K151_ARMS)), all(modes %in% c("fast", "medium", "slow")), !anyNA(seeds))
        for (ds in datasets) {
            scenario <- sub("@.*$", "", ds)
            index <- as.integer(sub("^.*@", "", ds))
            data <- generate_k151(scenario, index)
            for (arm_name in arms) {
                arm <- K151_ARMS[[arm_name]]
                views <- if (isTRUE(arm$fail)) NULL else prepare_views(data$views, arm)
                for (mode in modes) for (seed in seeds) {
                    key <- sprintf("%s|%s|%s|%d", ds, arm_name, mode, seed)
                    base <- data.frame(dataset = ds, data_seed = data$data_seed, N = data$N, arm = arm_name,
                        mode = mode, train_seed = seed, stringsAsFactors = FALSE)
                    run_cell(key, base, function() {
                        fit <- train(views, arm, mode, seed)
                        r <- abs(stats::cor(data$z, fit$Z))
                        match <- greedy_match(r)
                        row <- data.frame(iterations = fit$iterations, elbo = fit$elbo, learned_factors = ncol(fit$Z),
                            planted_F1 = match$abs_r[["F1"]], planted_F2 = match$abs_r[["F2"]],
                            planted_F3 = match$abs_r[["F3"]], init_captured = !is.null(fit$init_Z),
                            pca_solver = fit$pca_solver, capture_errors = fit$capture_errors,
                            mofa_warnings = fit$warnings, stringsAsFactors = FALSE)
                        fit$planted <- data$z
                        list(fit = fit, row = row)
                    })
                }
            }
        }
    } else if (task == "pin") {
        condition <- arguments[[3L]]
        arms <- split_list(arguments[[4L]])
        stopifnot(all(arms %in% names(K151_ARMS)))
        data <- generate_k151("n20weak", 7L)
        result$condition <- condition
        for (arm_name in arms) {
            arm <- K151_ARMS[[arm_name]]
            key <- sprintf("%s|%s", condition, arm_name)
            base <- data.frame(condition = condition, dataset = "n20weak@7", data_seed = data$data_seed, arm = arm_name,
                mode = "fast", train_seed = PIN_TRAIN_SEED, stringsAsFactors = FALSE)
            run_cell(key, base, function() {
                fit <- train(prepare_views(data$views, arm), arm, "fast", PIN_TRAIN_SEED)
                match <- greedy_match(abs(stats::cor(data$z, fit$Z)))
                row <- data.frame(iterations = fit$iterations, elbo = fit$elbo, learned_factors = ncol(fit$Z),
                    planted_F1 = match$abs_r[["F1"]], planted_F2 = match$abs_r[["F2"]],
                    planted_F3 = match$abs_r[["F3"]], mofa_warnings = fit$warnings, stringsAsFactors = FALSE)
                fit$init_Z <- NULL
                list(fit = fit, row = row)
            })
        }
    } else {
        N <- as.integer(arguments[[3L]])
        indices <- as.integer(split_list(arguments[[4L]]))
        arms <- split_list(arguments[[5L]])
        stopifnot(N %in% c(20L, 60L), !anyNA(indices), all(arms %in% RUV_ARMS))
        suppressPackageStartupMessages({
            library(RUVIIIC)
            library(ruv)
        })
        result$versions$RUVIIIC <- as.character(utils::packageVersion("RUVIIIC"))
        result$versions$ruv <- as.character(utils::packageVersion("ruv"))
        checks <- list()
        for (index in indices) {
            data <- generate_ruv(N, index)
            for (arm_name in arms) {
                key <- sprintf("N%d|%d|%s", N, index, arm_name)
                base <- data.frame(N = N, seed_index = index, data_seed = data$data_seed, arm = arm_name,
                    mode = RUV_MODE, train_seed = RUV_TRAIN_SEED, stringsAsFactors = FALSE)
                run_cell(key, base, function() {
                    if (arm_name == "probe_fail") stop("probe: deliberate failure of the lifecycle arm")
                    corrected <- list()
                    for (view in names(RUV_SIZES)) {
                        x <- data$views[[view]]
                        u <- data$latent[, if (view == "A") "UA" else "UB"]
                        corrected[[view]] <- switch(arm_name,
                            raw = x,
                            oracle = oracle_correct(x, u),
                            ruv_correct(x, data$group, as.integer(sub("ruv_k", "", arm_name))))
                        if (startsWith(arm_name, "ruv_k") && index == indices[[1L]]) {
                            k <- as.integer(sub("ruv_k", "", arm_name))
                            reference <- ruv_reference(x, data$group, k)
                            checks[[length(checks) + 1L]] <<- data.frame(N = N, seed_index = index, view = view, k = k,
                                max_abs_diff_vs_ruv_RUVIII = max(abs(corrected[[view]] - reference)),
                                max_abs_value = max(abs(x)), stringsAsFactors = FALSE)
                        }
                    }
                    keep_rows <- list()
                    for (view in names(RUV_SIZES)) {
                        l <- data$loaded[[view]]
                        u <- data$latent[, if (view == "A") "UA" else "UB"]
                        keep_rows[[paste0("retained_G_", view)]] <- retained(corrected[[view]], data$views[[view]],
                            l[, "G"], data$group)
                        keep_rows[[paste0("retained_S_", view)]] <- retained(corrected[[view]], data$views[[view]],
                            l[, "S"], data$latent[, "S"])
                        keep_rows[[paste0("retained_U_", view)]] <- retained(corrected[[view]], data$views[[view]],
                            l[, "U"], u)
                    }
                    views <- lapply(names(RUV_SIZES), function(view) {
                        m <- corrected[[view]]
                        rownames(m) <- paste(view, rownames(m), sep = ":")
                        m
                    })
                    names(views) <- names(RUV_SIZES)
                    fit <- train(views, base_arm(), RUV_MODE, RUV_TRAIN_SEED)
                    r <- abs(stats::cor(data$latent, fit$Z))
                    match <- greedy_match(r)
                    gf <- match$factor[["G"]]
                    if (is.na(gf)) gf <- colnames(fit$Z)[[which.max(r["G", ])]]
                    gz <- fit$Z[, gf]
                    fit_lm <- stats::lm(gz ~ factor(data$group))
                    an <- stats::anova(fit_lm)
                    row <- data.frame(iterations = fit$iterations, elbo = fit$elbo, learned_factors = ncol(fit$Z),
                        r_G = match$abs_r[["G"]], r_S = match$abs_r[["S"]], r_UA = match$abs_r[["UA"]],
                        r_UB = match$abs_r[["UB"]], best_r_S = max(r["S", ]), group_factor = gf,
                        group_factor_ve_A = unname(fit$r2[gf, "A"]), group_factor_ve_B = unname(fit$r2[gf, "B"]),
                        group_factor_R2 = summary(fit_lm)$r.squared, group_factor_p = an[["Pr(>F)"]][[1L]],
                        s_factor = match$factor[["S"]],
                        s_factor_ve_A = if (is.na(match$factor[["S"]])) NA_real_ else unname(fit$r2[match$factor[["S"]], "A"]),
                        s_factor_ve_B = if (is.na(match$factor[["S"]])) NA_real_ else unname(fit$r2[match$factor[["S"]], "B"]),
                        total_ve_A = sum(fit$r2[, "A"]), total_ve_B = sum(fit$r2[, "B"]),
                        as.data.frame(keep_rows), mofa_warnings = fit$warnings, stringsAsFactors = FALSE)
                    fit$init_Z <- NULL
                    fit$latent <- data$latent
                    list(fit = fit, row = row)
                })
            }
        }
        if (length(checks)) utils::write.table(do.call(rbind, checks), file.path(out_dir, "ruv-check.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
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
saveRDS(list(job = result, fits = fits), file.path(out_dir, "fits.rds"))
jsonlite::write_json(result, file.path(out_dir, "job.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE,
    null = "null", na = "null")
cat(sprintf("JOB %s: %d cells, %d failed, %s s, CPU %s\n", task, length(rows), failed, result$seconds_total,
    result$cpu_model))
quit(save = "no", status = status)
