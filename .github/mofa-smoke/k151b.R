# KIT-151 selection evidence, run b (branch run-mofa-evidence-k151b; declaration
# docs/qa/workflow-execution/kit/kit-151/selection-b-declaration.md). Audit SCI-01, SCI-02, SCI-03, SCI-13.
# The generator of k151-followup.R (scenario n20strong) with three switchable data features:
#   trend   : noise sd depends on the feature's abundance (its offset): sd = 1.5 * 2^(-(offset - 20) / 4), so the
#             median noise sd stays 1.5 (the earlier design) and doubles for every 4 units lower abundance;
#   missing : abundance-dependent missing values, P(missing) = 1 - (1 - plogis((17.5 - x) / 0.5)) * (1 - 0.02),
#             drawn per value on the generated (unscaled) value x, from a separate random stream;
#   scales  : each data set multiplied by its own factor (D1 1, D2 2, D3 0.25, D4 4) after the missing values,
#             so the data sets are on different scales (scale_views has something to do).
# All three off = the earlier n20strong data exactly (same random draws, same arithmetic).
# Selection arms: none; raw variance; random subset of the same size; variance above the mean-variance trend
# (difference, and log ratio); raw variance then per-feature z-score; each with eligibility "observed in >= 2 samples"
# (the current rule) or "observed in >= half the samples".
#
# Usage: Rscript --vanilla k151b.R <out dir> <seed index> <cells: all | scenario:arm,...> <maxiter> <probe: 0|1>
# Writes rows.tsv (cell x planted factor), composition.tsv (cell x data set x class), job.json, source-excerpts.txt.
# Exit status 1 when any cell failed (the failed cell is written as rows with status "failed").

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 5L)
out_dir <- arguments[[1L]]
seed_index <- as.integer(arguments[[2L]])
cell_request <- arguments[[3L]]
MAXITER <- as.integer(arguments[[4L]])
probe <- identical(arguments[[5L]], "1")
stopifnot(seed_index %in% c(1:30, 91:92), MAXITER >= 1000L)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()

# ------------------------------------------------------------------ fixed design (declared)
SIZES <- c(D1 = 20000L, D2 = 2000L, D3 = 200L, D4 = 100L)
ACTIVE <- list(D1 = c("F1", "F3"), D2 = "F1", D3 = c("F1", "F2"), D4 = c("F1", "F2"))
LOADED_SHARE <- 0.10
EFFECT <- 1.0
F2_EFFECT <- 1.5 # the follow-up's "strong" F2
NOISE_SD <- 1.5
F1_NOISE_SD <- 0.5
TREND_REF <- 20
TREND_DOUBLING <- 4
MISS_CENTRE <- 17.5
MISS_WIDTH <- 0.5
MISS_FLOOR <- 0.02
VIEW_SCALE <- c(D1 = 1, D2 = 2, D3 = 0.25, D4 = 4)
HALF_SHARE <- 0.5
TREND_SPAN <- 0.3
DATA_SEED <- 151000L + seed_index
TRAIN_SEED <- 20261119L
N <- 20L
SCENARIOS <- list(
    prev = list(trend = FALSE, missing = FALSE, scales = FALSE),
    trend = list(trend = TRUE, missing = FALSE, scales = TRUE),
    full = list(trend = TRUE, missing = TRUE, scales = TRUE))
arm <- function(select = "none", cap = NA_integer_, eligible = "two", zscore = FALSE, scale_views = TRUE,
                fail = FALSE) {
    list(select = select, cap = cap, eligible = eligible, zscore = zscore, scale_views = scale_views,
        center_groups = TRUE, fail = fail)
}
ARMS <- list(
    none = arm(scale_views = FALSE),
    scale_views = arm(),
    raw250 = arm("raw", 250L),
    raw1000 = arm("raw", 1000L),
    raw2000 = arm("raw", 2000L),
    random1000 = arm("random", 1000L),
    trend1000 = arm("trend_diff", 1000L),
    trendratio1000 = arm("trend_ratio", 1000L),
    raw1000_z = arm("raw", 1000L, zscore = TRUE))
for (a in c("raw250", "raw1000", "raw2000", "random1000", "trend1000", "trendratio1000", "raw1000_z")) {
    ARMS[[paste0(a, "_obs50")]] <- ARMS[[a]]
    ARMS[[paste0(a, "_obs50")]]$eligible <- "half"
}
ARMS$probe_fail <- arm(fail = TRUE)
BASE_ARMS <- c("none", "scale_views", "raw250", "raw1000", "raw2000", "random1000", "trend1000", "trendratio1000",
    "raw1000_z")
OBS_ARMS <- paste0(c("raw250", "raw1000", "raw2000", "random1000", "trend1000", "trendratio1000", "raw1000_z"),
    "_obs50")
if (identical(cell_request, "all")) {
    cells <- c(paste0("prev:", BASE_ARMS), paste0("trend:", BASE_ARMS), paste0("full:", c(BASE_ARMS, OBS_ARMS)))
} else {
    cells <- strsplit(cell_request, ",")[[1L]]
}
if (probe) cells <- c(cells, "full:probe_fail")
cell_scenario <- sub(":.*$", "", cells)
cell_arm <- sub("^.*:", "", cells)
stopifnot(all(cell_scenario %in% names(SCENARIOS)), all(cell_arm %in% names(ARMS)))

# ------------------------------------------------------------------ generator
# The draws of k151-followup.R in the same order; the noise is drawn as standard normals and multiplied by each
# feature's sd (with all switches off this is the same arithmetic as rnorm(d * N, 0, 1.5)).
generate <- function(sc) {
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
    noise_sd <- list()
    for (view in names(SIZES)) {
        d <- SIZES[[view]]
        w <- matrix(0, d, 3L, dimnames = list(NULL, colnames(z)))
        for (k in ACTIVE[[view]]) {
            loaded <- sample.int(d, round(LOADED_SHARE * d))
            w[loaded, k] <- sample(c(-1, 1), length(loaded), replace = TRUE) * (if (k == "F2") F2_EFFECT else EFFECT)
        }
        offset <- stats::rnorm(d, 20, 2)
        sd_j <- if (sc$trend) NOISE_SD * 2^(-(offset - TREND_REF) / TREND_DOUBLING) else rep(NOISE_SD, d)
        x <- w %*% t(z) + matrix(stats::rnorm(d * N), d, N) * sd_j + offset
        dimnames(x) <- list(sprintf("%s_f%05d", view, seq_len(d)), samples)
        rownames(w) <- rownames(x)
        views[[view]] <- x
        loadings[[view]] <- w != 0
        noise_sd[[view]] <- stats::setNames(sd_j, rownames(x))
    }
    if (sc$missing) {
        set.seed(DATA_SEED + 400000L)
        for (view in names(SIZES)) {
            x <- views[[view]]
            p <- 1 - (1 - stats::plogis((MISS_CENTRE - x) / MISS_WIDTH)) * (1 - MISS_FLOOR)
            x[matrix(stats::runif(length(x)), nrow(x)) < p] <- NA_real_
            views[[view]] <- x
        }
    }
    if (sc$scales) for (view in names(SIZES)) views[[view]] <- views[[view]] * VIEW_SCALE[[view]]
    list(views = views, z = z, loadings = loadings, noise_sd = noise_sd)
}

# ------------------------------------------------------------------ view preparation
# As integrationAlignViews (eligibility, then the z-score spread filter, then selection, then the z-score), with the
# selection statistic and the eligibility rule as arm settings.
fitted_trend <- function(abundance, y) {
    fit <- stats::lowess(abundance, y, f = TREND_SPAN)
    stats::approx(fit$x, fit$y, xout = abundance, ties = base::mean, rule = 2)$y
}

prepare_views <- function(views, a) {
    out <- list()
    for (vi in seq_along(views)) {
        view <- names(views)[[vi]]
        m <- views[[view]]
        need <- if (identical(a$eligible, "half")) ceiling(HALF_SHARE * ncol(m)) else 2L
        m <- m[rowSums(!is.na(m)) >= need, , drop = FALSE]
        if (a$zscore) {
            spread <- apply(m, 1L, stats::sd, na.rm = TRUE)
            m <- m[spread > 0 & !is.na(spread), , drop = FALSE]
        }
        if (!is.na(a$cap) && nrow(m) > a$cap) {
            variance <- apply(m, 1L, stats::var, na.rm = TRUE)
            score <- switch(a$select,
                raw = variance,
                random = {
                    set.seed(DATA_SEED + 300000L + vi)
                    stats::runif(nrow(m))
                },
                trend_diff = variance - fitted_trend(rowMeans(m, na.rm = TRUE), variance),
                trend_ratio = log(variance) - fitted_trend(rowMeans(m, na.rm = TRUE), log(variance)),
                stop("unknown selection ", a$select))
            if (!all(is.finite(score))) stop("selection score is not finite for every feature in ", view)
            keep <- sort(order(-score, seq_along(score))[seq_len(a$cap)])
            m <- m[keep, , drop = FALSE]
        }
        if (a$zscore) {
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

train <- function(views, a) {
    if (isTRUE(a$fail)) stop("probe: deliberate failure of the lifecycle cell")
    object <- MOFA2::create_mofa(views)
    lists <- list(data = MOFA2::get_default_data_options(object), model = MOFA2::get_default_model_options(object),
        training = MOFA2::get_default_training_options(object))
    lists$data$scale_views <- a$scale_views
    lists$data$center_groups <- a$center_groups
    lists$data$use_float32 <- FALSE
    lists$model$num_factors <- 5L
    lists$training$seed <- TRAIN_SEED
    lists$training$maxiter <- MAXITER
    lists$training$convergence_mode <- "slow"
    lists$training$verbose <- FALSE
    weight_list <- tryCatch(locate_option(lists, "weight_views"), error = function(e) NA_character_)
    if (!is.na(weight_list)) lists[[weight_list]]$weight_views <- FALSE
    object <- MOFA2::prepare_mofa(object, data_options = lists$data, model_options = lists$model,
        training_options = lists$training)
    outfile <- tempfile("k151b-", fileext = ".hdf5")
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

# Kept and dropped features per data set by what they carry (F1, F2, F3; a feature can carry two; none), with the
# median true noise sd and the median observed fraction of the kept features of each class.
composition <- function(views, data, raw_views) {
    out <- list()
    for (view in names(SIZES)) {
        kept <- sub("^[^:]*:", "", rownames(views[[view]]))
        is_kept <- rownames(data$loadings[[view]]) %in% kept
        l <- data$loadings[[view]]
        observed <- rowMeans(!is.na(raw_views[[view]]))
        classes <- list(F1 = l[, "F1"], F2 = l[, "F2"], F3 = l[, "F3"], none = rowSums(l) == 0L,
            all = rep(TRUE, nrow(l)))
        for (k in names(classes)) {
            sel <- classes[[k]] & is_kept
            out[[length(out) + 1L]] <- data.frame(view = view, carries = k, total = sum(classes[[k]]),
                kept = sum(sel), dropped = sum(classes[[k]] & !is_kept),
                kept_noise_sd_median = if (any(sel)) stats::median(data$noise_sd[[view]][sel]) else NA_real_,
                kept_observed_median = if (any(sel)) stats::median(observed[sel]) else NA_real_,
                stringsAsFactors = FALSE)
        }
    }
    do.call(rbind, out)
}

result <- list(seed_index = seed_index, data_seed = DATA_SEED, train_seed = TRAIN_SEED, N = N, maxiter = MAXITER,
    cells = cells)
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

    data_sets <- list()
    result$missing_share <- list()
    for (sc_name in unique(cell_scenario)) {
        data_sets[[sc_name]] <- generate(SCENARIOS[[sc_name]])
        result$missing_share[[sc_name]] <- lapply(data_sets[[sc_name]]$views, function(v) round(mean(is.na(v)), 4))
    }
    for (i in seq_along(cells)) {
        sc_name <- cell_scenario[[i]]
        arm_name <- cell_arm[[i]]
        a <- ARMS[[arm_name]]
        data <- data_sets[[sc_name]]
        clock <- Sys.time()
        base <- data.frame(scenario = sc_name, N = N, seed_index = seed_index, data_seed = DATA_SEED, arm = arm_name,
            stringsAsFactors = FALSE)
        outcome <- tryCatch({
            views <- prepare_views(data$views, a)
            comp <- composition(views, data, data$views)
            compositions[[cells[[i]]]] <- cbind(base[rep(1L, nrow(comp)), , drop = FALSE], comp)
            kept <- vapply(views, nrow, integer(1))
            fit <- train(views, a)
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
            out$elbo_rel_change_last5 <- if (length(elbo) > 6L) {
                abs((utils::tail(elbo, 1L) - elbo[[length(elbo) - 5L]]) / elbo[[1L]]) * 100
            } else NA_real_
            out$learned_factors <- ncol(factors)
            out$factors_dropped <- 5L - ncol(factors)
            out$mofa_warnings <- length(grep("^WARNING:", fit$messages))
            out
        }, error = function(e) {
            failed <<- failed + 1L
            message("CELL FAILED ", cells[[i]], ": ", conditionMessage(e))
            out <- data.frame(planted = c("F1", "F2", "F3"), status = "failed", error = conditionMessage(e),
                matched_factor = NA_character_, abs_r = NA_real_, stringsAsFactors = FALSE)
            for (view in names(SIZES)) out[[paste0("r2_", view)]] <- NA_real_
            for (view in names(SIZES)) out[[paste0("kept_", view)]] <- NA_integer_
            out$iterations <- NA_integer_
            out$stopped_by <- NA_character_
            out$elbo_final <- NA_real_
            out$elbo_rel_change_last5 <- NA_real_
            out$learned_factors <- NA_integer_
            out$factors_dropped <- NA_integer_
            out$mofa_warnings <- NA_integer_
            out
        })
        outcome$seconds <- round(as.numeric(difftime(Sys.time(), clock, units = "secs")), 2)
        rows[[cells[[i]]]] <- cbind(base[rep(1L, nrow(outcome)), , drop = FALSE], outcome)
        cat(sprintf("CELL %s seed %d: %s; |r| %s; %s iterations; %s s\n", cells[[i]], seed_index,
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
cat(sprintf("JOB seed %d: %d cells, %d failed, %s s, CPU %s\n", seed_index, length(cells), failed,
    result$seconds_total, result$cpu_model))
quit(save = "no", status = status)
