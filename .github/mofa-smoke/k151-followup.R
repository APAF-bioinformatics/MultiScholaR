# KIT-151 cap evidence, follow-up (branch run-mofa-evidence-k151; declaration
# docs/qa/workflow-execution/kit/kit-151/cap-followup-declaration.md).
# Same generator, view preparation and MOFA2 call as k151-evidence.R (the first run), with three added dimensions:
# convergence mode (fast, medium, slow), smaller caps (250, 500) and F2 strength (weak = the first run; strong = F2
# weight +/-1.5 in D3 and D4). Plus a sensitivity data variant ("noisy": 10% of D1's features loaded on no factor get
# noise sd x 2) for cap 1,000 and cap 2,000. Also records, for the capped arms, how many kept and dropped features in
# each data set are loaded on each planted factor.
#
# Usage: Rscript --vanilla k151-followup.R <out dir> <scenario: n20weak|n20strong|n60weak|n60strong> <seed index 1..10>
#   <cells: all | arm@mode,arm@mode,...> <modes for "all", comma list>
# Writes rows.tsv (one row per cell x planted factor), composition.tsv, job.json, source-excerpts.txt.
# Exit status 1 when any cell failed (the failed cell is written as rows with status "failed").

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 5L)
out_dir <- arguments[[1L]]
scenario <- arguments[[2L]]
seed_index <- as.integer(arguments[[3L]])
cell_request <- arguments[[4L]]
modes_all <- strsplit(arguments[[5L]], ",")[[1L]]
stopifnot(scenario %in% c("n20weak", "n20strong", "n60weak", "n60strong"), seed_index %in% 1:10,
    all(modes_all %in% c("fast", "medium", "slow")))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()

# ------------------------------------------------------------------ fixed design (declared)
SIZES <- c(D1 = 20000L, D2 = 2000L, D3 = 200L, D4 = 100L)
ACTIVE <- list(D1 = c("F1", "F3"), D2 = "F1", D3 = c("F1", "F2"), D4 = c("F1", "F2"))
LOADED_SHARE <- 0.10
EFFECT <- 1.0
STRENGTH <- if (grepl("strong$", scenario)) "strong" else "weak"
F2_EFFECT <- if (STRENGTH == "strong") 1.5 else 1.0
NOISE_SD <- 1.5
F1_NOISE_SD <- 0.5
INFLATED_SHARE <- 0.10 # of D1's features, drawn from those loaded on no factor
INFLATED_SD_FACTOR <- 2
DATA_SEED <- 151000L + seed_index
TRAIN_SEED <- 20261119L
MAXITER <- 1000L
N <- if (startsWith(scenario, "n60")) 60L else 20L
base_arm <- function(cap = NA_integer_, zscore = FALSE, scale_views = TRUE, center_groups = TRUE,
                     weight_views = FALSE, noisy = FALSE, fail = FALSE) {
    list(cap = cap, zscore = zscore, scale_views = scale_views, center_groups = center_groups,
        weight_views = weight_views, noisy = noisy, fail = fail)
}
ARMS <- list(
    none = base_arm(scale_views = FALSE),
    scale_views = base_arm(),
    cap250 = base_arm(cap = 250L),
    cap500 = base_arm(cap = 500L),
    cap1000 = base_arm(cap = 1000L),
    cap2000 = base_arm(cap = 2000L),
    cap5000 = base_arm(cap = 5000L),
    weight_views = base_arm(weight_views = TRUE),
    zscore = base_arm(zscore = TRUE, scale_views = FALSE, center_groups = FALSE),
    noisy_cap1000 = base_arm(cap = 1000L, noisy = TRUE),
    noisy_cap2000 = base_arm(cap = 2000L, noisy = TRUE),
    probe_fail = base_arm(fail = TRUE)) # lifecycle only
MAIN_ARMS <- c("none", "scale_views", "cap250", "cap500", "cap1000", "cap2000", "cap5000", "weight_views", "zscore")
NOISY_ARMS <- c("noisy_cap1000", "noisy_cap2000")
if (identical(cell_request, "all")) {
    arms_here <- c(MAIN_ARMS, if (N == 20L) NOISY_ARMS)
    cells <- unlist(lapply(modes_all, function(mode) paste(arms_here, mode, sep = "@")))
} else {
    cells <- strsplit(cell_request, ",")[[1L]]
}
cell_arm <- sub("@.*$", "", cells)
cell_mode <- sub("^.*@", "", cells)
stopifnot(all(cell_arm %in% names(ARMS)), all(cell_mode %in% c("fast", "medium", "slow")))

# ------------------------------------------------------------------ generator (the first run's, F2 effect a parameter)
# The random draws are those of k151-evidence.R in the same order, so weak data equal the first run's data; strong
# data differ only in F2's loadings (x 1.5); the noisy variant adds noise to chosen D1 features from its own stream.
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
    loadings <- list()
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
        rownames(w) <- rownames(x)
        views[[view]] <- x
        loadings[[view]] <- w != 0
    }
    set.seed(DATA_SEED + 200000L)
    unloaded <- which(rowSums(loadings$D1) == 0L)
    inflated <- sort(unloaded[sample.int(length(unloaded), round(INFLATED_SHARE * SIZES[["D1"]]))])
    extra_sd <- sqrt((INFLATED_SD_FACTOR * NOISE_SD)^2 - NOISE_SD^2)
    noisy <- views
    noisy$D1[inflated, ] <- noisy$D1[inflated, ] + matrix(stats::rnorm(length(inflated) * N, 0, extra_sd),
        length(inflated), N)
    inflated_flags <- lapply(views, function(v) stats::setNames(rep(FALSE, nrow(v)), rownames(v)))
    inflated_flags$D1[inflated] <- TRUE
    list(views = views, noisy = noisy, z = z, loadings = loadings, inflated = inflated_flags)
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

locate_option <- function(lists, name) {
    where <- names(lists)[vapply(lists, function(l) name %in% names(l), logical(1))]
    if (length(where) != 1L) stop(sprintf("option %s found in %d MOFA2 option lists", name, length(where)))
    where
}

train <- function(views, arm, mode) {
    if (isTRUE(arm$fail)) stop("probe: deliberate failure of the lifecycle arm")
    object <- MOFA2::create_mofa(views)
    lists <- list(data = MOFA2::get_default_data_options(object), model = MOFA2::get_default_model_options(object),
        training = MOFA2::get_default_training_options(object))
    lists$data$scale_views <- arm$scale_views
    lists$data$center_groups <- arm$center_groups
    lists$data$use_float32 <- FALSE
    lists$model$num_factors <- 5L
    lists$training$seed <- TRAIN_SEED
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
    outfile <- tempfile("k151-", fileext = ".hdf5")
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
    list(model = model, messages = messages)
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

# Kept and dropped features per data set by what they carry: loaded on F1, F2, F3 (a feature can carry two), loaded
# on none, inflated noise (the noisy variant only; a subset of "none").
composition <- function(views, data, arm_name) {
    out <- list()
    for (view in names(SIZES)) {
        kept <- sub("^[^:]*:", "", rownames(views[[view]]))
        is_kept <- rownames(data$loadings[[view]]) %in% kept
        l <- data$loadings[[view]]
        classes <- list(F1 = l[, "F1"], F2 = l[, "F2"], F3 = l[, "F3"], none = rowSums(l) == 0L,
            inflated = data$inflated[[view]] & isTRUE(ARMS[[arm_name]]$noisy))
        for (k in names(classes)) {
            out[[length(out) + 1L]] <- data.frame(view = view, carries = k, total = sum(classes[[k]]),
                kept = sum(classes[[k]] & is_kept), dropped = sum(classes[[k]] & !is_kept), stringsAsFactors = FALSE)
        }
    }
    do.call(rbind, out)
}

result <- list(scenario = scenario, strength = STRENGTH, seed_index = seed_index, data_seed = DATA_SEED,
    train_seed = TRAIN_SEED, N = N, cells = cells)
result$cpu <- tryCatch(system("lscpu", intern = TRUE), error = function(e) conditionMessage(e))
result$cpu_model <- trimws(sub("^Model name:", "", grep("^Model name:", result$cpu, value = TRUE)[1L]))
result$nproc <- parallel::detectCores()
rows <- list()
compositions <- list()
failed <- 0L

status <- tryCatch({
    suppressPackageStartupMessages(library(MOFA2))
    reticulate::py_run_string(paste(
        "import sys, importlib.metadata as md",
        "k151_versions = {'python': sys.version.split()[0], 'mofapy2': md.version('mofapy2'),",
        "    'numpy': md.version('numpy'), 'scipy': md.version('scipy'), 'h5py': md.version('h5py')}", sep = "\n"))
    result$versions <- c(list(R = as.character(getRversion()), MOFA2 = as.character(utils::packageVersion("MOFA2")),
        reticulate = as.character(utils::packageVersion("reticulate"))), reticulate::import_main()$k151_versions)
    # The installed mofapy2 lines that set the convergence criterion (recorded, not interpreted).
    reticulate::py_run_string(paste(
        "import mofapy2, pathlib",
        "root = pathlib.Path(mofapy2.__file__).parent",
        "pats = ('convergence_mode', 'tolerance', 'onverged')",
        "k151_src = ['mofapy2 installed at ' + str(root)]",
        "for p in sorted(root.rglob('*.py')):",
        "    lines = p.read_text(errors='replace').splitlines()",
        "    for i, line in enumerate(lines, 1):",
        "        if any(s in line for s in pats):",
        "            k151_src.append('---- ' + str(p.relative_to(root.parent)) + ':' + str(i))",
        "            for j in range(max(1, i - 4), min(len(lines), i + 8) + 1):",
        "                k151_src.append(f'{p.relative_to(root.parent)}:{j}: {lines[j - 1]}')",
        sep = "\n"))
    writeLines(reticulate::import_main()$k151_src, file.path(out_dir, "source-excerpts.txt"))

    data <- generate()
    prepared <- list()
    for (i in seq_along(cells)) {
        arm_name <- cell_arm[[i]]
        mode <- cell_mode[[i]]
        arm <- ARMS[[arm_name]]
        clock <- Sys.time()
        base <- data.frame(scenario = scenario, strength = STRENGTH, N = N, seed_index = seed_index,
            data_seed = DATA_SEED, arm = arm_name, mode = mode, stringsAsFactors = FALSE)
        outcome <- tryCatch({
            if (is.null(prepared[[arm_name]])) {
                prepared[[arm_name]] <- prepare_views(if (arm$noisy) data$noisy else data$views, arm)
                if (!is.na(arm$cap)) {
                    comp <- composition(prepared[[arm_name]], data, arm_name)
                    compositions[[arm_name]] <- cbind(base[rep(1L, nrow(comp)), c("scenario", "strength", "N",
                        "seed_index", "arm")], comp)
                }
            }
            views <- prepared[[arm_name]]
            kept <- vapply(views, nrow, integer(1))
            fit <- train(views, arm, mode)
            factors <- MOFA2::get_factors(fit$model)[[1L]]
            factors <- factors[rownames(data$z), , drop = FALSE]
            r2 <- MOFA2::get_variance_explained(fit$model)$r2_per_factor[[1L]]
            if (all(rownames(r2) %in% names(views))) r2 <- t(r2)
            r <- abs(stats::cor(data$z, factors))
            match <- greedy_match(r)
            stats_slot <- methods::slot(fit$model, "training_stats")
            iterations <- length(stats_slot$time)
            elbo <- as.numeric(stats_slot$elbo)
            elbo <- elbo[is.finite(elbo)]
            out <- data.frame(planted = rownames(r), status = "ok", error = NA_character_,
                matched_factor = unname(match$factor), abs_r = unname(match$abs_r), stringsAsFactors = FALSE)
            for (view in names(SIZES)) {
                out[[paste0("r2_", view)]] <- vapply(match$factor, function(f) if (is.na(f)) NA_real_ else
                    unname(r2[f, view]), numeric(1))
            }
            for (view in names(SIZES)) out[[paste0("kept_", view)]] <- unname(kept[[view]])
            out$iterations <- iterations
            out$stopped_by <- if (iterations >= MAXITER) "maxiter" else "convergence"
            out$elbo_final <- if (length(elbo)) utils::tail(elbo, 1L) else NA_real_
            out$learned_factors <- ncol(factors)
            out$factors_dropped <- 5L - ncol(factors)
            out$mofa_warnings <- length(grep("^WARNING:", fit$messages))
            out
        }, error = function(e) {
            failed <<- failed + 1L
            message("CELL FAILED ", arm_name, "@", mode, ": ", conditionMessage(e))
            out <- data.frame(planted = c("F1", "F2", "F3"), status = "failed", error = conditionMessage(e),
                matched_factor = NA_character_, abs_r = NA_real_, stringsAsFactors = FALSE)
            for (view in names(SIZES)) out[[paste0("r2_", view)]] <- NA_real_
            for (view in names(SIZES)) out[[paste0("kept_", view)]] <- NA_integer_
            out$iterations <- NA_integer_
            out$stopped_by <- NA_character_
            out$elbo_final <- NA_real_
            out$learned_factors <- NA_integer_
            out$factors_dropped <- NA_integer_
            out$mofa_warnings <- NA_integer_
            out
        })
        outcome$seconds <- round(as.numeric(difftime(Sys.time(), clock, units = "secs")), 2)
        rows[[cells[[i]]]] <- cbind(base[rep(1L, nrow(outcome)), , drop = FALSE], outcome)
        cat(sprintf("CELL %s %s@%s seed %d: %s; |r| %s; %s iterations; %s s\n", scenario, arm_name, mode, seed_index,
            outcome$status[[1L]], paste(sprintf("%s=%.3f", outcome$planted, outcome$abs_r), collapse = " "),
            outcome$iterations[[1L]], outcome$seconds[[1L]]))
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
if (length(compositions)) {
    utils::write.table(do.call(rbind, compositions), file.path(out_dir, "composition.tsv"), sep = "\t", quote = FALSE,
        row.names = FALSE, na = "NA")
}
result$failed_cells <- failed
result$status <- status
result$seconds_total <- round(as.numeric(difftime(Sys.time(), started, units = "secs")), 2)
jsonlite::write_json(result, file.path(out_dir, "job.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE,
    null = "null", na = "null")
cat(sprintf("JOB %s seed %d: %d cells, %d failed, %s s, CPU %s\n", scenario, seed_index, length(cells), failed,
    result$seconds_total, result$cpu_model))
quit(save = "no", status = status)
