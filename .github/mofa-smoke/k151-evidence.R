# KIT-151 cap evidence (branch run-mofa-evidence-k151; declaration
# docs/qa/workflow-execution/kit/kit-151/cap-evidence-declaration.md).
# One job: one scenario x one data seed, every requested arm trained in turn with MOFA2 called directly
# (create_mofa -> default options with the declared ones set -> prepare_mofa -> run_mofa(use_basilisk = FALSE)),
# as integrationTrainMofa in R/utils_workflow_integration_model.R (origin/kit-142-integration-record b20d4fd7) does.
# The view preparation follows integrationAlignViews there: features observed in fewer than two samples dropped;
# under the z-score, features without spread dropped; the cap keeps the most variable features (variance per feature
# with na.rm, ties to the earlier row, kept rows in their original order); the z-score is (x - mean) / sd with na.rm.
#
# Usage: Rscript --vanilla k151-evidence.R <out dir> <scenario: main|n60|missing> <seed index 1..10> <arms: all|a,b,...>
# Writes rows.tsv (one row per arm x planted factor), job.json, source-excerpts.txt, save-data-check.json.
# Exit status 1 when any arm failed (the failed arm is written as a row with status "failed").

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 4L)
out_dir <- arguments[[1L]]
scenario <- arguments[[2L]]
seed_index <- as.integer(arguments[[3L]])
arm_request <- arguments[[4L]]
stopifnot(scenario %in% c("main", "n60", "missing"), seed_index %in% 1:10)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()

# ------------------------------------------------------------------ fixed design (declared)
SIZES <- c(D1 = 20000L, D2 = 2000L, D3 = 200L, D4 = 100L)
ACTIVE <- list(D1 = c("F1", "F3"), D2 = "F1", D3 = c("F1", "F2"), D4 = c("F1", "F2"))
LOADED_SHARE <- 0.10
EFFECT <- 1.0
NOISE_SD <- 1.5
F1_NOISE_SD <- 0.5
MISSING_SHARE <- 0.15
DATA_SEED <- 151000L + seed_index
TRAIN_SEED <- 20261119L
N <- if (scenario == "n60") 60L else 20L
ARMS <- list(
    none = list(cap = NA_integer_, zscore = FALSE, scale_views = FALSE, center_groups = TRUE, weight_views = FALSE),
    scale_views = list(cap = NA_integer_, zscore = FALSE, scale_views = TRUE, center_groups = TRUE, weight_views = FALSE),
    cap1000 = list(cap = 1000L, zscore = FALSE, scale_views = TRUE, center_groups = TRUE, weight_views = FALSE),
    cap2000 = list(cap = 2000L, zscore = FALSE, scale_views = TRUE, center_groups = TRUE, weight_views = FALSE),
    cap5000 = list(cap = 5000L, zscore = FALSE, scale_views = TRUE, center_groups = TRUE, weight_views = FALSE),
    cap10000 = list(cap = 10000L, zscore = FALSE, scale_views = TRUE, center_groups = TRUE, weight_views = FALSE),
    weight_views = list(cap = NA_integer_, zscore = FALSE, scale_views = TRUE, center_groups = TRUE, weight_views = TRUE),
    zscore = list(cap = NA_integer_, zscore = TRUE, scale_views = FALSE, center_groups = FALSE, weight_views = FALSE),
    # Lifecycle only (never in "all"): fails on purpose, to show a failed arm is reported as failed.
    probe_fail = list(cap = NA_integer_, zscore = FALSE, scale_views = TRUE, center_groups = TRUE, weight_views = FALSE,
        fail = TRUE))
arms <- if (identical(arm_request, "all")) setdiff(names(ARMS), "probe_fail") else strsplit(arm_request, ",")[[1L]]
stopifnot(all(arms %in% names(ARMS)))

# ------------------------------------------------------------------ generator (declared)
generate <- function() {
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(DATA_SEED)
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
            w[loaded, k] <- sample(c(-1, 1), length(loaded), replace = TRUE) * EFFECT
        }
        offset <- stats::rnorm(d, 20, 2)
        x <- w %*% t(z) + matrix(stats::rnorm(d * N, 0, NOISE_SD), d, N) + offset
        dimnames(x) <- list(sprintf("%s_f%05d", view, seq_len(d)), samples)
        views[[view]] <- x
    }
    if (scenario == "missing") {
        set.seed(DATA_SEED + 100000L)
        for (view in names(views)) {
            mask <- matrix(stats::runif(length(views[[view]])) < MISSING_SHARE, nrow(views[[view]]))
            views[[view]][mask] <- NA_real_
        }
    }
    list(views = views, z = z, group = group)
}

# ------------------------------------------------------------------ view preparation (as integrationAlignViews)
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

# The MOFA2 option list that holds a name (data, model or training); stops if not exactly one.
locate_option <- function(lists, name) {
    where <- names(lists)[vapply(lists, function(l) name %in% names(l), logical(1))]
    if (length(where) != 1L) stop(sprintf("option %s found in %d MOFA2 option lists", name, length(where)))
    where
}

train <- function(views, arm, save_data = TRUE) {
    if (isTRUE(arm$fail)) stop("probe: deliberate failure of the lifecycle arm")
    object <- MOFA2::create_mofa(views)
    lists <- list(data = MOFA2::get_default_data_options(object), model = MOFA2::get_default_model_options(object),
        training = MOFA2::get_default_training_options(object))
    lists$data$scale_views <- arm$scale_views
    lists$data$center_groups <- arm$center_groups
    lists$data$use_float32 <- FALSE
    lists$model$num_factors <- 5L
    lists$training$seed <- TRAIN_SEED
    lists$training$maxiter <- 1000L
    lists$training$convergence_mode <- "fast"
    lists$training$verbose <- FALSE
    # weight_views: set where MOFA2's option lists hold it; absent, only the weight_views arm stops.
    weight_list <- tryCatch(locate_option(lists, "weight_views"), error = function(e) {
        if (arm$weight_views) stop(conditionMessage(e))
        NA_character_
    })
    if (!is.na(weight_list)) lists[[weight_list]]$weight_views <- arm$weight_views
    object <- MOFA2::prepare_mofa(object, data_options = lists$data, model_options = lists$model,
        training_options = lists$training)
    outfile <- tempfile("k151-", fileext = ".hdf5")
    on.exit(unlink(outfile), add = TRUE)
    messages <- character()
    withCallingHandlers(
        MOFA2::run_mofa(object, outfile = outfile, save_data = save_data, use_basilisk = FALSE),
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
    list(model = model, messages = messages, weight_list = weight_list,
        h5 = tryCatch(unique(rhdf5::h5ls(outfile, recursive = 1L)$name), error = function(e) conditionMessage(e)))
}

# Greedy maximum matching: the largest |r| first, one learned factor per planted factor; an unmatched planted
# factor (fewer learned factors than planted) gets |r| = 0 and no factor.
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

result <- list(scenario = scenario, seed_index = seed_index, data_seed = DATA_SEED, train_seed = TRAIN_SEED, N = N,
    arms = arms)
result$cpu <- tryCatch(system("lscpu", intern = TRUE), error = function(e) conditionMessage(e))
result$cpu_model <- trimws(sub("^Model name:", "", grep("^Model name:", result$cpu, value = TRUE)[1L]))
result$nproc <- parallel::detectCores()
rows <- list()
failed <- 0L

status <- tryCatch({
    suppressPackageStartupMessages(library(MOFA2))
    reticulate::py_run_string(paste(
        "import sys, importlib.metadata as md",
        "k151_versions = {'python': sys.version.split()[0], 'mofapy2': md.version('mofapy2'),",
        "    'numpy': md.version('numpy'), 'scipy': md.version('scipy'), 'h5py': md.version('h5py')}", sep = "\n"))
    result$versions <- c(list(R = as.character(getRversion()), MOFA2 = as.character(utils::packageVersion("MOFA2")),
        reticulate = as.character(utils::packageVersion("reticulate"))), reticulate::import_main()$k151_versions)

    # -------------------------------------------------------------- installed source excerpts (recorded, not interpreted)
    reticulate::py_run_string(paste(
        "import mofapy2, pathlib",
        "root = pathlib.Path(mofapy2.__file__).parent",
        "pats = ('scale_views', 'weight_views', 'center_groups', 'scale_groups')",
        "k151_src = ['mofapy2 installed at ' + str(root)]",
        "for p in sorted(root.rglob('*.py')):",
        "    lines = p.read_text(errors='replace').splitlines()",
        "    for i, line in enumerate(lines, 1):",
        "        if any(s in line for s in pats):",
        "            k151_src.append('---- ' + str(p.relative_to(root.parent)) + ':' + str(i))",
        "            for j in range(max(1, i - 4), min(len(lines), i + 8) + 1):",
        "                k151_src.append(f'{p.relative_to(root.parent)}:{j}: {lines[j - 1]}')",
        sep = "\n"))
    python_source <- reticulate::import_main()$k151_src
    ns <- asNamespace("MOFA2")
    r_source <- c(sprintf("MOFA2 %s installed at %s (functions deparsed from the installed namespace; line numbers are within each deparsed function)",
        utils::packageVersion("MOFA2"), find.package("MOFA2")))
    for (fn in sort(ls(ns, all.names = TRUE))) {
        f <- get(fn, envir = ns)
        if (!is.function(f)) next
        text <- deparse(f, width.cutoff = 500L)
        hit <- grep("scale_views|weight_views|center_groups|variance_explained|save_data", text)
        if (!length(hit)) next
        r_source <- c(r_source, sprintf("==== MOFA2::%s", fn))
        for (i in hit) r_source <- c(r_source, sprintf("MOFA2::%s:%d: %s", fn, i, text[[i]]))
    }
    writeLines(c("##### mofapy2 (installed .py files)", python_source, "", "##### MOFA2 R (installed namespace)",
        r_source), file.path(out_dir, "source-excerpts.txt"))

    # -------------------------------------------------------------- save_data = FALSE check (tiny model)
    tiny <- local({
        RNGkind("Mersenne-Twister", "Inversion", "Rejection")
        set.seed(151999L)
        z <- stats::rnorm(20)
        m <- function(d, v) {
            x <- outer(stats::rnorm(d), z) + matrix(stats::rnorm(d * 20), d)
            dimnames(x) <- list(sprintf("%s_%03d", v, seq_len(d)), sprintf("S%02d", 1:20))
            x
        }
        list(A = m(60L, "A"), B = m(40L, "B"))
    })
    tiny_arm <- ARMS$scale_views
    check <- list()
    for (save in c(TRUE, FALSE)) {
        key <- if (save) "save_data_true" else "save_data_false"
        check[[key]] <- tryCatch({
            fit <- train(tiny, tiny_arm, save_data = save)
            ve <- tryCatch(MOFA2::get_variance_explained(fit$model), error = function(e) e)
            list(trained = TRUE, h5_top_level = fit$h5,
                variance_explained_ok = !inherits(ve, "error"),
                variance_explained_error = if (inherits(ve, "error")) conditionMessage(ve) else NA,
                r2_per_factor = if (!inherits(ve, "error")) as.vector(ve$r2_per_factor[[1L]]) else NA,
                r2_total = if (!inherits(ve, "error")) as.vector(ve$r2_total[[1L]]) else NA,
                cache_has_variance_explained = !is.null(fit$model@cache$variance_explained),
                messages = fit$messages)
        }, error = function(e) list(trained = FALSE, error = conditionMessage(e)))
    }
    if (isTRUE(check$save_data_true$variance_explained_ok) && isTRUE(check$save_data_false$variance_explained_ok)) {
        check$max_abs_difference_r2_per_factor <- max(abs(check$save_data_true$r2_per_factor -
            check$save_data_false$r2_per_factor))
    }
    jsonlite::write_json(check, file.path(out_dir, "save-data-check.json"), auto_unbox = TRUE, digits = NA,
        pretty = TRUE, null = "null", na = "null")

    # -------------------------------------------------------------- the arms
    data <- generate()
    result$generated <- lapply(data$views, function(v) list(features = nrow(v), samples = ncol(v),
        missing = sum(is.na(v))))
    for (arm_name in arms) {
        arm <- ARMS[[arm_name]]
        clock <- Sys.time()
        base <- data.frame(scenario = scenario, seed_index = seed_index, data_seed = DATA_SEED, N = N, arm = arm_name,
            stringsAsFactors = FALSE)
        outcome <- tryCatch({
            views <- prepare_views(data$views, arm)
            kept <- vapply(views, nrow, integer(1))
            fit <- train(views, arm)
            factors <- MOFA2::get_factors(fit$model)[[1L]]
            factors <- factors[rownames(data$z), , drop = FALSE]
            r2 <- MOFA2::get_variance_explained(fit$model)$r2_per_factor[[1L]]
            if (all(rownames(r2) %in% names(views))) r2 <- t(r2)
            r <- abs(stats::cor(data$z, factors))
            match <- greedy_match(r)
            iterations <- length(methods::slot(fit$model, "training_stats")$time)
            out <- data.frame(planted = rownames(r), status = "ok", error = NA_character_,
                matched_factor = unname(match$factor), abs_r = unname(match$abs_r), stringsAsFactors = FALSE)
            for (view in names(SIZES)) {
                out[[paste0("r2_", view)]] <- vapply(match$factor, function(f) if (is.na(f)) NA_real_ else
                    unname(r2[f, view]), numeric(1))
            }
            for (view in names(SIZES)) out[[paste0("kept_", view)]] <- unname(kept[[view]])
            out$iterations <- iterations
            out$learned_factors <- ncol(factors)
            out$factors_dropped <- 5L - ncol(factors)
            out$weight_views_option_list <- fit$weight_list
            out$mofa_warnings <- length(grep("^WARNING:", fit$messages))
            out
        }, error = function(e) {
            failed <<- failed + 1L
            message("ARM FAILED ", arm_name, ": ", conditionMessage(e))
            out <- data.frame(planted = c("F1", "F2", "F3"), status = "failed", error = conditionMessage(e),
                matched_factor = NA_character_, abs_r = NA_real_, stringsAsFactors = FALSE)
            for (view in names(SIZES)) out[[paste0("r2_", view)]] <- NA_real_
            for (view in names(SIZES)) out[[paste0("kept_", view)]] <- NA_integer_
            out$iterations <- NA_integer_
            out$learned_factors <- NA_integer_
            out$factors_dropped <- NA_integer_
            out$weight_views_option_list <- NA_character_
            out$mofa_warnings <- NA_integer_
            out
        })
        outcome$seconds <- round(as.numeric(difftime(Sys.time(), clock, units = "secs")), 2)
        rows[[arm_name]] <- cbind(base[rep(1L, nrow(outcome)), , drop = FALSE], outcome)
        cat(sprintf("ARM %s %s seed %d: %s; |r| %s; %s s\n", scenario, arm_name, seed_index, outcome$status[[1L]],
            paste(sprintf("%s=%.3f", outcome$planted, outcome$abs_r), collapse = " "), outcome$seconds[[1L]]))
    }
    if (failed > 0L) 1L else 0L
}, error = function(e) {
    result$error <<- conditionMessage(e)
    message("JOB ERROR: ", conditionMessage(e))
    1L
})

if (length(rows)) {
    table <- do.call(rbind, rows)
    table$error <- gsub("[\t\r\n]+", " ", table$error)
    utils::write.table(table, file.path(out_dir, "rows.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
result$failed_arms <- failed
result$status <- status
result$seconds_total <- round(as.numeric(difftime(Sys.time(), started, units = "secs")), 2)
jsonlite::write_json(result, file.path(out_dir, "job.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE,
    null = "null", na = "null")
cat(sprintf("JOB %s seed %d: %d arms, %d failed, %s s, CPU %s\n", scenario, seed_index, length(arms), failed,
    result$seconds_total, result$cpu_model))
quit(save = "no", status = status)
