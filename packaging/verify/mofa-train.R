# KIT-134 release-candidate run (branch run-k134-rc): MOFA2 trained in double precision INSIDE the assembled,
# installed bundle (its R, its library's MOFA2 and reticulate, its Python with mofapy2), network closed. Copied
# unchanged from .github/mofa-smoke/diagnose.R of branch run-mofa-smoke-k062 (KIT-062), plus the MOFA lines printed
# at the end (the result stays in the job). performance-v2 has no integration model yet, so MOFA2 is called directly.
#
# Original header:
# KIT-062 diagnosis (branch run-mofa-smoke-k062): why two Linux runs of the same MOFA2 training differ.
# One job trains once, on the KIT-137 tier-1 corpus, with the views MultiScholaR's integration model builds
# (R/utils_workflow_integration_model.R at 449b9de7, provisional parameters: omics in name order, views in manifest
# order, internal standards dropped, runs averaged per biological sample (tech_rep_group), samples in C-locale order
# with missing values where a view lacks one, features observed in fewer than two samples dropped, feature names
# "<view>:<feature_id>"; MOFA2 defaults except scale_views FALSE, center_groups TRUE, 5 factors, fast convergence,
# maxiter 1000, seed 20261119), and records the machine, the BLAS libraries and thread counts, and the trained
# values as hexadecimal floats.
#
# Usage: Rscript --vanilla diagnose.R <out dir> <label> <corpus dir> <precision: float32|float64>
# float64: the option that turns float32 off is found by name in the installed MOFA2's default option lists
# (a name containing "float32"); the job stops if there is none or more than one. Its location and the mofapy2 source
# lines that read it are recorded.

arguments <- commandArgs(trailingOnly = TRUE)
# packaging/verify/windows-clean-run.ps1 passes the bundle's folder first (as smoke.R expects); it is not used here.
if (length(arguments) == 5L) arguments <- arguments[-1L]
stopifnot(length(arguments) == 4L)
out_dir <- arguments[[1L]]
label <- arguments[[2L]]
corpus <- arguments[[3L]]
precision <- arguments[[4L]]
stopifnot(precision %in% c("float32", "float64"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()
hex <- function(x) sprintf("%a", as.vector(x))
md5 <- function(x) {
    file <- tempfile()
    writeBin(as.vector(x), file)
    unname(tools::md5sum(file))
}
result <- list(label = label, precision = precision, trained = FALSE)
result$environment <- as.list(Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS", "RETICULATE_PYTHON"), unset = NA_character_))
result$platform <- R.version$platform
if (.Platform$OS.type == "windows") {
    result$cpu <- tryCatch(system2("powershell", c("-NoProfile", "-Command", "(Get-CimInstance Win32_Processor).Name"),
        stdout = TRUE), error = function(e) conditionMessage(e))
    result$cpu_model <- trimws(result$cpu[[1L]])
} else {
    result$cpu <- tryCatch(system("lscpu", intern = TRUE), error = function(e) conditionMessage(e))
    result$cpu_model <- trimws(sub("^Model name:", "", grep("^Model name:", result$cpu, value = TRUE)[1L]))
}
result$nproc <- parallel::detectCores()
result$r_blas <- list(blas = unname(extSoftVersion()["BLAS"]), lapack = La_library(), la_version = La_version())

status <- tryCatch({
    suppressPackageStartupMessages(library(MOFA2))
    result$mofa2 <- as.character(utils::packageVersion("MOFA2"))
    reticulate::py_run_string(paste(
        "import io, contextlib, sys, importlib.metadata as md, numpy, threadpoolctl",
        "_b = io.StringIO()",
        "with contextlib.redirect_stdout(_b):",
        "    numpy.show_config()",
        "diag_numpy_config = _b.getvalue()",
        "diag_versions = {'python': sys.version.split()[0], 'numpy': numpy.__version__, 'mofapy2': md.version('mofapy2'),",
        "    'scipy': md.version('scipy')}",
        sep = "\n"))
    main <- reticulate::import_main()
    result$numpy_config <- main$diag_numpy_config
    result$python_versions <- main$diag_versions
    threads <- function() {
        reticulate::py_run_string(paste("import threadpoolctl",
            "diag_threads = [{k: str(v) for k, v in d.items() if k in ('user_api', 'internal_api', 'num_threads', 'version', 'filepath', 'architecture')}",
            "    for d in threadpoolctl.threadpool_info()]", sep = "\n"))
        reticulate::import_main()$diag_threads
    }
    result$python_threads_before <- threads()
    # The mofapy2 source lines that read a float32 setting (recorded, not interpreted).
    reticulate::py_run_string(paste("import mofapy2, pathlib, re",
        "diag_float_lines = []",
        "for p in sorted(pathlib.Path(mofapy2.__file__).parent.rglob('*.py')):",
        "    for i, line in enumerate(p.read_text(errors='replace').splitlines(), 1):",
        "        if 'float32' in line: diag_float_lines.append(f'{p.name}:{i}: {line.strip()}')", sep = "\n"))
    result$mofapy2_float32_lines <- reticulate::import_main()$diag_float_lines

    # ------------------------------------------------------------ the views, as the package builds them
    omics <- sort(c("proteomics", "transcriptomics", "lipidomics", "metabolomics"), method = "radix")
    views <- list()
    for (omic in omics) {
        manifest <- jsonlite::read_json(file.path(corpus, omic, "integration_input.json"))
        samples <- utils::read.delim(file.path(corpus, omic, "samples.tsv"), colClasses = "character")
        key_column <- manifest$biological_sample_column
        for (view in manifest$views) {
            table <- utils::read.delim(file.path(corpus, omic, view$matrix_file), colClasses = "character",
                check.names = FALSE)
            features <- utils::read.delim(file.path(corpus, omic, view$features_file), colClasses = "character")
            values <- vapply(table[-1L], as.numeric, numeric(nrow(table)))
            values <- matrix(values, nrow = nrow(table), dimnames = list(table[[1L]], names(table)[-1L]))
            if (!is.null(features$is_internal_standard)) {
                values <- values[!(toupper(features$is_internal_standard) %in% "TRUE"), , drop = FALSE]
            }
            key <- samples[[key_column]][match(colnames(values), samples[[manifest$sample_id_column]])]
            levels <- sort(unique(key), method = "radix")
            averaged <- vapply(levels, function(level) {
                m <- rowMeans(values[, key == level, drop = FALSE], na.rm = TRUE)
                m[is.nan(m)] <- NA_real_
                m
            }, numeric(nrow(values)))
            averaged <- matrix(averaged, nrow = nrow(values), dimnames = list(rownames(values), levels))
            views[[view$name]] <- averaged
        }
    }
    kept <- sort(unique(unlist(lapply(views, colnames))), method = "radix")
    views <- lapply(stats::setNames(names(views), names(views)), function(view) {
        m <- matrix(NA_real_, nrow(views[[view]]), length(kept), dimnames = list(rownames(views[[view]]), kept))
        m[, colnames(views[[view]])] <- views[[view]]
        m <- m[rowSums(!is.na(m)) >= 2L, , drop = FALSE]
        rownames(m) <- paste(view, rownames(m), sep = ":")
        m
    })
    result$views <- lapply(views, function(v) list(features = nrow(v), samples = ncol(v), md5 = md5(v)))

    # ------------------------------------------------------------ options
    object <- create_mofa(views)
    data_options <- get_default_data_options(object)
    model_options <- get_default_model_options(object)
    training_options <- get_default_training_options(object)
    data_options$scale_views <- FALSE
    data_options$center_groups <- TRUE
    model_options$num_factors <- 5L
    training_options$seed <- 20261119L
    training_options$maxiter <- 1000L
    training_options$convergence_mode <- "fast"
    training_options$verbose <- FALSE
    lists <- list(data = data_options, model = model_options, training = training_options)
    found <- unlist(lapply(names(lists), function(n) paste(n, grep("float32", names(lists[[n]]), value = TRUE), sep = "$")))
    found <- found[!grepl("\\$$", found)]
    result$float32_option <- found
    result$float32_default <- if (length(found) == 1L) {
        parts <- strsplit(found, "$", fixed = TRUE)[[1L]]
        lists[[parts[[1L]]]][[parts[[2L]]]]
    } else NA
    if (identical(precision, "float64")) {
        if (length(found) != 1L) stop("no single MOFA2 option naming float32: ", paste(found, collapse = ", "))
        parts <- strsplit(found, "$", fixed = TRUE)[[1L]]
        if (parts[[1L]] == "data") data_options[[parts[[2L]]]] <- FALSE
        if (parts[[1L]] == "model") model_options[[parts[[2L]]]] <- FALSE
        if (parts[[1L]] == "training") training_options[[parts[[2L]]]] <- FALSE
    }
    result$options <- list(data = data_options, model = model_options, training = training_options)
    object <- prepare_mofa(object, data_options = data_options, model_options = model_options,
        training_options = training_options)

    # ------------------------------------------------------------ train
    outfile <- file.path(tempdir(), "diagnose.hdf5")
    clock <- Sys.time()
    messages <- utils::capture.output(run_mofa(object, outfile = outfile, save_data = TRUE, use_basilisk = FALSE),
        type = "message")
    result$seconds_train <- round(as.numeric(difftime(Sys.time(), clock, units = "secs")), 2)
    result$float_message <- grep("float", messages, value = TRUE)
    result$python_threads_after <- threads()
    model <- load_model(outfile)
    factors <- get_factors(model)[[1L]]
    weights <- get_weights(model)
    r2 <- get_variance_explained(model)$r2_per_factor[[1L]]
    if (all(rownames(r2) %in% names(views))) r2 <- t(r2)
    elbo <- as.numeric(methods::slot(model, "training_stats")$elbo)
    elbo <- elbo[is.finite(elbo)]
    truth <- utils::read.delim(file.path(corpus, "truth_scores.tsv"))
    planted <- as.matrix(truth[match(rownames(factors), truth$biological_sample), c("Z1", "Z2", "Z3")])
    correlation <- abs(stats::cor(planted, factors))
    result$model <- list(
        iterations = length(methods::slot(model, "training_stats")$time),
        elbo_last_hex = hex(utils::tail(elbo, 1L)), elbo_last = utils::tail(elbo, 1L),
        factors = list(samples = rownames(factors), values_hex = hex(factors), md5 = md5(factors), dim = dim(factors)),
        weights = lapply(weights, function(w) list(values_hex = hex(w), md5 = md5(w), dim = dim(w))),
        r2 = list(factors = rownames(r2), views = colnames(r2), values_hex = hex(r2), md5 = md5(r2)),
        dtype_note = paste(class(factors[1L, 1L]), "in R; MOFA2 stores what mofapy2 wrote"),
        planted_abs_r = as.list(apply(correlation, 1L, max)))
    result$trained <- TRUE
    0L
}, error = function(e) {
    result$error <<- conditionMessage(e)
    message("DIAGNOSE ERROR: ", conditionMessage(e))
    1L
})
result$seconds_total <- round(as.numeric(difftime(Sys.time(), started, units = "secs")), 2)
jsonlite::write_json(result, file.path(out_dir, "result.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE,
    null = "null", na = "null")
cat(sprintf("DIAGNOSE %s (%s; %s): trained %s; %s iterations; train %s s; CPU %s\n", label, precision,
    result$environment$OPENBLAS_NUM_THREADS, result$trained, result$model$iterations, result$seconds_train,
    result$cpu_model))
# KIT-134: the values for the comparison with KIT-062's double-precision runs (44466163909), printed to the log.
cat(sprintf("MOFA ENV python %s; MOFA2 %s; mofapy2 %s; numpy %s; lapack %s; float32 option %s set off: %s\n",
    paste(result$environment$RETICULATE_PYTHON), paste(result$mofa2), paste(result$python_versions$mofapy2),
    paste(result$python_versions$numpy), paste(result$r_blas$lapack), paste(result$float32_option, collapse = ","),
    identical(precision, "float64")))
cat(sprintf("MOFA LIBPATHS %s\n", paste(.libPaths(), collapse = " | ")))
cat(sprintf("MOFA PACKAGE MOFA2 from %s\n", tryCatch(find.package("MOFA2"), error = conditionMessage)))
if (isTRUE(result$trained)) {
    cat(sprintf("MOFA RESULT %s: iterations %d; ELBO %.6f (%s); factors md5 %s; planted |r| %s\n", label,
        result$model$iterations, result$model$elbo_last, result$model$elbo_last_hex, result$model$factors$md5,
        paste(sprintf("%.4f", unlist(result$model$planted_abs_r)), collapse = " ")))
    cat(sprintf("MOFA FACTORS %s %s\n", label, paste(result$model$factors$values_hex, collapse = " ")))
    cat(sprintf("MOFA R2 %s %s\n", label, paste(result$model$r2$values_hex, collapse = " ")))
} else cat(sprintf("MOFA RESULT %s: NOT TRAINED: %s\n", label, paste(result$error)))
quit(save = "no", status = status)
