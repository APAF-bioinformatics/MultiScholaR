# MOFA2 smoke (run-mofa-smoke; KIT-060 risks 2 and 6, KIT-134): train one MOFA2 model with the bundle's
# R library and the bundle's Python, started the way the bundle starts R (Linux: the launcher
# packaging/linux/Rscript; Windows: packaging/verify/windows-clean-run.ps1, the shell's environment).
#
# Usage: Rscript --vanilla smoke.R <bundle> <out dir> <label> <corpus dir> <factors> <maxiter> <seed>
#        <convergence mode> <mode: train|fail>
# Input: the KIT-137 tier-1 synthetic corpus (20 biological samples, six views in four omics, planted
# scores in truth_scores.tsv). Smoke choices, not methods decisions: runs of one biological sample
# are averaged per view; views are named <omic>_<assay>; features are prefixed with the view; every
# MOFA2 option not named here stays at its default. Writes <out dir>/result.json (always, also on
# failure) and exits 1 when the model did not train.

started <- Sys.time()
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 9L)
bundle <- arguments[[1L]]
out_dir <- arguments[[2L]]
label <- arguments[[3L]]
corpus <- arguments[[4L]]
num_factors <- as.integer(arguments[[5L]])
maxiter <- as.integer(arguments[[6L]])
seed <- as.integer(arguments[[7L]])
convergence <- arguments[[8L]]
mode <- arguments[[9L]]
stopifnot(mode %in% c("train", "fail"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
is_windows <- .Platform$OS.type == "windows"

result <- list(label = label, mode = mode, platform = R.version$platform, os = utils::sessionInfo()$running,
    trained = FALSE, error = NULL,
    options = list(num_factors = num_factors, maxiter = maxiter, seed = seed, convergence_mode = convergence,
        use_basilisk = FALSE))
seconds <- function(since) round(as.numeric(difftime(Sys.time(), since, units = "secs")), 2)
hex <- function(x) sprintf("%a", as.vector(x))
md5 <- function(x) {
    file <- tempfile()
    writeBin(as.vector(x), file)
    unname(tools::md5sum(file))
}
writeResult <- function() {
    result$seconds_total <- seconds(started)
    jsonlite::write_json(result, file.path(out_dir, "result.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE,
        null = "null", na = "null")
}

status <- tryCatch({
    # ------------------------------------------------------------ environment and versions
    result$environment <- as.list(Sys.getenv(c("RETICULATE_PYTHON", "PYTHONNOUSERSITE", "R_LIBS_USER", "R_LIBS_SITE",
        "PYTHONPATH", "PYTHONHOME", "RETICULATE_PYTHON_ENV", "CONDA_PREFIX", "VIRTUAL_ENV", "MSR_BUNDLE"),
        unset = NA_character_))
    result$lib_paths <- .libPaths()
    suppressPackageStartupMessages(library(MOFA2))
    r_versions <- vapply(c("MOFA2", "reticulate", "basilisk", "rhdf5", "Rhdf5lib", "HDF5Array", "DelayedArray",
        "Matrix", "jsonlite"), function(p) as.character(utils::packageVersion(p)), character(1))
    result$r_versions <- c(R = as.character(getRversion()), as.list(r_versions))
    result$run_mofa_formals <- lapply(formals(MOFA2::run_mofa), function(x) paste(deparse(x), collapse = ""))
    result$hdf5_r <- tryCatch(paste(utils::capture.output(print(rhdf5::h5version())), collapse = " "),
        error = function(e) paste("unavailable:", conditionMessage(e)))

    clock <- Sys.time()
    py_config <- reticulate::py_config()
    result$python_config <- list(python = py_config$python, libpython = py_config$libpython,
        version = py_config$version_string, forced = py_config$forced, virtualenv = py_config$virtualenv)
    expected <- normalizePath(Sys.getenv("RETICULATE_PYTHON"), winslash = "/", mustWork = FALSE)
    result$python_is_bundle <- identical(normalizePath(py_config$python, winslash = "/", mustWork = FALSE), expected)
    reticulate::py_run_string(paste(
        "import sys, importlib.metadata as md, numpy, h5py, threadpoolctl",
        "smoke_versions = {'python': sys.version.split()[0], 'executable': sys.executable, 'prefix': sys.prefix,",
        "    'numpy': numpy.__version__, 'h5py': h5py.__version__, 'hdf5_h5py': h5py.version.hdf5_version}",
        "for p in ['mofapy2', 'scipy', 'scikit-learn', 'pandas', 'dtw-python']: smoke_versions[p] = md.version(p)",
        "smoke_threads = [{k: str(v) for k, v in d.items() if k in ('user_api', 'internal_api', 'num_threads', 'version', 'filepath')}",
        "    for d in threadpoolctl.threadpool_info()]",
        "smoke_paths = [p for p in sys.path]", sep = "\n"))
    main <- reticulate::import_main()
    result$python_versions <- main$smoke_versions
    result$python_threads <- main$smoke_threads
    result$python_sys_path <- main$smoke_paths
    result$seconds_python <- seconds(clock)
    if (!isTRUE(result$python_is_bundle)) stop("reticulate's Python is not the bundle's: ", py_config$python)

    # ------------------------------------------------------------ input
    truth <- utils::read.delim(file.path(corpus, "truth_scores.tsv"), check.names = FALSE)
    samples <- truth$biological_sample
    files <- sort(list.files(corpus, recursive = TRUE, full.names = TRUE))
    result$corpus_md5 <- as.list(stats::setNames(unname(tools::md5sum(files)), substring(files, nchar(corpus) + 2L)))
    views <- list()
    for (omic in c("proteomics", "transcriptomics", "lipidomics", "metabolomics")) {
        design <- utils::read.delim(file.path(corpus, omic, "samples.tsv"), check.names = FALSE)
        for (file in sort(list.files(file.path(corpus, omic), pattern = "^matrix_.*[.]tsv$", full.names = TRUE))) {
            assay <- sub("^matrix_(.*)[.]tsv$", "\\1", basename(file))
            view <- if (identical(assay, omic)) omic else paste(omic, assay, sep = "_")
            matrix <- as.matrix(utils::read.delim(file, row.names = 1L, check.names = FALSE))
            averaged <- vapply(samples, function(sample) {
                runs <- intersect(design$Run[design$biological_sample == sample], colnames(matrix))
                if (!length(runs)) return(rep(NA_real_, nrow(matrix)))
                values <- rowMeans(matrix[, runs, drop = FALSE], na.rm = TRUE)
                values[is.nan(values)] <- NA_real_
                values
            }, numeric(nrow(matrix)))
            dimnames(averaged) <- list(paste(view, rownames(matrix), sep = ":"), samples)
            views[[view]] <- averaged
        }
    }
    result$views <- lapply(views, function(v) list(features = nrow(v), samples = ncol(v),
        missing_fraction = round(mean(is.na(v)), 4), samples_absent = sum(colSums(!is.na(v)) == 0L)))

    # ------------------------------------------------------------ create, options, prepare
    clock <- Sys.time()
    object <- create_mofa(views)
    data_options <- get_default_data_options(object)
    model_options <- get_default_model_options(object)
    model_options$num_factors <- num_factors
    training_options <- get_default_training_options(object)
    training_options$seed <- seed
    training_options$maxiter <- maxiter
    training_options$convergence_mode <- convergence
    training_options$verbose <- FALSE
    result$data_options <- data_options
    result$model_options <- model_options
    result$training_options <- training_options
    object <- prepare_mofa(object, data_options = data_options, model_options = model_options,
        training_options = training_options)
    result$seconds_prepare <- seconds(clock)
    if (identical(mode, "fail")) stop("lifecycle: deliberate failure after prepare_mofa, before run_mofa")

    # ------------------------------------------------------------ train
    model_file <- file.path(tempdir(), "mofa-smoke-model.hdf5")
    messages <- character()
    clock <- Sys.time()
    trained <- withCallingHandlers(
        run_mofa(object, outfile = model_file, save_data = TRUE, use_basilisk = FALSE),
        message = function(m) { messages <<- c(messages, conditionMessage(m)); invokeRestart("muffleMessage") },
        warning = function(w) { messages <<- c(messages, paste("WARNING:", conditionMessage(w))); invokeRestart("muffleWarning") })
    result$seconds_train <- seconds(clock)
    result$run_mofa_messages <- trimws(messages)
    result$model_file_bytes <- file.size(model_file)
    # run_mofa falls back to basilisk (a conda environment) when reticulate's Python lacks mofapy2.
    result$basilisk_fallback <- any(grepl("Setting use_basilisk = TRUE", messages, fixed = TRUE))
    if (result$basilisk_fallback) stop("run_mofa fell back to basilisk: the bundle's Python was not used")

    # ------------------------------------------------------------ load and report
    clock <- Sys.time()
    model <- load_model(model_file)
    result$seconds_load <- seconds(clock)
    factors <- get_factors(model)[[1L]]
    weights <- get_weights(model)
    r2 <- get_variance_explained(model)$r2_per_factor[[1L]]
    if (all(rownames(r2) %in% names(views))) r2 <- t(r2)
    elbo <- tryCatch(as.numeric(model@training_stats$elbo), error = function(e) NA_real_)
    result$model <- list(
        factors_requested = num_factors,
        factors_in_model = ncol(factors),
        factors_active_1pct = sum(apply(r2, 1L, max) >= 1),
        iterations = sum(!is.na(elbo)),
        elbo_last = utils::tail(elbo[!is.na(elbo)], 1L),
        elbo_last_hex = hex(utils::tail(elbo[!is.na(elbo)], 1L)),
        r2_percent = list(factors = rownames(r2), views = colnames(r2), values = unname(r2), values_hex = hex(r2)),
        r2_total_percent = as.list(get_variance_explained(model)$r2_total[[1L]]),
        factors = list(samples = rownames(factors), names = colnames(factors), values = unname(factors),
            values_hex = hex(factors), md5 = md5(factors)),
        weights = lapply(weights, function(w) list(rows = nrow(w), md5 = md5(w), values_hex = hex(w))))
    planted <- as.matrix(truth[match(rownames(factors), truth$biological_sample), c("Z1", "Z2", "Z3")])
    correlation <- stats::cor(planted, factors, use = "pairwise.complete.obs")
    result$planted <- list(
        correlation = list(planted = rownames(correlation), recovered = colnames(correlation),
            values = unname(round(correlation, 4))),
        best = lapply(stats::setNames(rownames(correlation), rownames(correlation)), function(z) {
            k <- which.max(abs(correlation[z, ]))
            list(factor = colnames(correlation)[[k]], r = round(correlation[z, k], 4))
        }))
    result$trained <- TRUE
    0L
}, error = function(e) {
    result$error <<- conditionMessage(e)
    message("SMOKE ERROR: ", conditionMessage(e))
    1L
})

# Peak memory of this process (R and the embedded Python share it), if Python came up.
result$peak_memory_mib <- tryCatch({
    if (!reticulate::py_available(initialize = FALSE)) NA_real_ else {
        reticulate::py_run_string(if (is_windows) paste(
            "import ctypes",
            "from ctypes import wintypes",
            "class _PMC(ctypes.Structure):",
            "    _fields_ = [('cb', wintypes.DWORD), ('PageFaultCount', wintypes.DWORD)] + [(n, ctypes.c_size_t) for n in",
            "        ('PeakWorkingSetSize', 'WorkingSetSize', 'QuotaPeakPagedPoolUsage', 'QuotaPagedPoolUsage',",
            "         'QuotaPeakNonPagedPoolUsage', 'QuotaNonPagedPoolUsage', 'PagefileUsage', 'PeakPagefileUsage')]",
            "_c = _PMC(); _c.cb = ctypes.sizeof(_c)",
            "ctypes.windll.kernel32.GetCurrentProcess.restype = ctypes.c_void_p",
            "ctypes.windll.psapi.GetProcessMemoryInfo.argtypes = [ctypes.c_void_p, ctypes.POINTER(_PMC), wintypes.DWORD]",
            "ctypes.windll.psapi.GetProcessMemoryInfo(ctypes.windll.kernel32.GetCurrentProcess(), ctypes.byref(_c), _c.cb)",
            "smoke_peak_mib = _c.PeakWorkingSetSize / 1048576", sep = "\n") else paste(
            "import resource",
            "smoke_peak_mib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024", sep = "\n"))
        round(reticulate::import_main()$smoke_peak_mib, 1)
    }
}, error = function(e) paste("unavailable:", conditionMessage(e)))

writeResult()
cat(sprintf("SMOKE %s (%s): trained %s; %s s in all\n", label, mode, result$trained, seconds(started)))
if (isTRUE(result$trained)) {
    cat(sprintf("MOFA2 %s, mofapy2 %s, Python %s, numpy %s, h5py %s (HDF5 %s); %d factors in the model, %d active (>= 1%%); %d iterations; train %s s\n",
        result$r_versions$MOFA2, result$python_versions$mofapy2, result$python_versions$python,
        result$python_versions$numpy, result$python_versions$h5py, result$python_versions$hdf5_h5py,
        result$model$factors_in_model, result$model$factors_active_1pct, result$model$iterations, result$seconds_train))
    for (z in names(result$planted$best)) cat(sprintf("  %s -> %s, r = %s\n", z, result$planted$best[[z]]$factor,
        result$planted$best[[z]]$r))
}
quit(save = "no", status = status)
