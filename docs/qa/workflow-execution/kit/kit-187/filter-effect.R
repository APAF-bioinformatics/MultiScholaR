# KIT-187 S2 (VAL-71, VAL-83 D5): the effect of the count rule on simulated data, base against head (ticket section
# 5.3; acceptance criterion A5). Evidence code under docs/qa: no package code.
#
# Usage, from the package root of the bound commit on a run-tests runner (full git history):
#   Rscript docs/qa/workflow-execution/kit/kit-187/filter-effect.R <base commit | self> <out dir>
#   <base commit>  40 digits: checked out as a git worktree and loaded as "base"; this tree is "head". `self`: this tree
#                  is loaded as both sides (the lifecycle's mechanics and fault runs).
# Environment: K187_LIFECYCLE=1 (R = 2 per cell and mask; falsifiers counted and written without stopping unless
#   K187_FAULT=1), K187_FAULT=1 (the "head" side drops one kept protein in cell K, replicate 1, mask vd5; the stop rule
#   on), K187_CELLS (comma list of cells; default A,B,C,D,E,K), K187_WORKERS (forked workers per side; default 2).
# Each side runs in a fresh Rscript (`--side`) with pkgload::load_all(<tree>, export_all = TRUE), the two sides at once.
# Per cell, mask and replicate: KIT-185's generator with strength 0 (no unwanted factor), the value-dependent mask, the
# filter through the S4 method at minimums 2 and 2 and percentile 1, then KIT-185's DA (`daTable`) on the kept set.
# Written to <out dir>: cells.tsv, falsifiers.tsv, fixtures.tsv, provenance.json. Outside lifecycle mode, and in the
# fault run, the script stops with a classed error (k187_falsifier) after the first cell in which a falsifier
# triggers, having written the rows of the cells up to and including it.
#
# The functions below marked "verbatim" are copied unchanged from KIT-185's generator,
# docs/qa/workflow-execution/kit/kit-185/sim-k185.R at commit 432f703b (branch kit-185-ruviiic-missing), at the line
# ranges named; the file itself cannot be sourced (it stops without its five arguments and attaches the installed
# package). S2's cells take idx = 200 + i, so the copied seed line gives 187000000 + 10000 i + r.

# ---- verbatim: sim-k185.R:43-52 at 432f703b ----
warning_count <- new.env()
warning_count$n <- 0L
quiet <- function(expr) {
    value <- NULL
    invisible(utils::capture.output(value <- withCallingHandlers(suppressMessages(expr), warning = function(w) {
        warning_count$n <- warning_count$n + 1L
        invokeRestart("muffleWarning")
    }), type = "output"))
    value
}

# ---- verbatim: sim-k185.R:71-71 at 432f703b ----
GROUPS <- c("control", "treated", "g3", "g4", "g5")

# ---- verbatim: sim-k185.R:72-76 at 432f703b ----
cell <- function(idx, omic, features, groups, per, k = 2L, level = 20, rho = 0, dup = FALSE, batches = 0L,
                 batch_runs = 0L, factors = 2L, assays = 1L) {
    list(idx = idx, omic = omic, features = features, groups = groups, per = per, k = k, level = level, rho = rho,
        dup = dup, batches = batches, batch_runs = batch_runs, factors = factors, assays = assays)
}

# ---- verbatim: sim-k185.R:90-91 at 432f703b ----
MED_VAR_PROT <- 4 * 0.25^2 / stats::qchisq(0.5, 4)
MED_VAR_OMIC <- 0.25

# ---- verbatim: sim-k185.R:101-108 at 432f703b ----
maskMech <- function(mask) if (identical(mask, "none")) "none" else sub("[0-9.]+$", "", mask)
maskRate <- function(mask) if (identical(mask, "none")) 0 else as.numeric(sub("^[a-z]+", "", mask)) / 100

# ---- data (KIT-185 5.2) --------------------------------------------------------------------------------------------
scale1 <- function(v) {
    v <- v - mean(v)
    v / stats::sd(v)
}

# ---- verbatim: sim-k185.R:109-199 at 432f703b ----
generateReplicate <- function(cellname, r, strength = 1.5) {
    p <- CELLS[[cellname]]
    set.seed(185000000L + 10000L * p$idx + r, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    n <- p$features
    # the runs: groups, units (P5), batches (D5 cells)
    if (p$batches > 0L) {
        per_batch <- if (identical(cellname, "P4b")) 2L else if (identical(cellname, "L1b")) 16L else 3L
        batch <- rep(seq_len(p$batches), each = p$groups * per_batch)
        group <- rep(rep(GROUPS[seq_len(p$groups)], each = per_batch), times = p$batches)
    } else {
        batch <- NULL
        group <- rep(GROUPS[seq_len(p$groups)], each = p$per * if (p$dup) 2L else 1L)
    }
    runs <- length(group)
    unit <- if (p$dup) rep(seq_len(p$groups * p$per), each = 2L) else seq_len(runs)
    run_ids <- sprintf("run_%03d", seq_len(runs))
    group_factor <- factor(group, levels = GROUPS[seq_len(p$groups)])
    median_variance <- if (identical(p$omic, "prot")) MED_VAR_PROT else MED_VAR_OMIC
    # levels
    mu <- stats::rnorm(n, p$level, 2)
    # variances
    s2 <- if (identical(p$omic, "prot")) 4 * 0.25^2 / stats::rchisq(n, 4) else {
        x <- 10 / stats::rchisq(n, 10)
        0.25 * x / stats::median(x)
    }
    # effects: exactly 20 % of the features, random positions, half up, in group 2
    n_effects <- round(0.2 * n)
    where <- sample.int(n, n_effects)
    signs <- sample(rep(c(1, -1), c(round(n_effects / 2), n_effects - round(n_effects / 2))))
    beta <- numeric(n)
    beta[where] <- signs
    # factors: centred, sample variance 1, exactly orthogonal to each other and to the group design (and batches)
    basis <- stats::model.matrix(~group_factor)
    planted <- NULL
    if (p$batches > 0L) {
        dummies <- stats::model.matrix(~ factor(batch))[, -1L, drop = FALSE]
        dummies <- scale(dummies, center = TRUE, scale = FALSE)
        basis <- cbind(basis, dummies)
        colnames(dummies) <- paste0("B", seq_len(ncol(dummies)))
        planted <- dummies
    }
    g <- scale1(as.numeric(group == "treated"))
    factor_scores <- NULL
    for (f in seq_len(p$factors)) {
        raw <- if (!identical(p$omic, "prot") && f == 1L) {
            order <- sample.int(runs)
            exp(-2 * (order - 1) / (runs - 1))
        } else stats::rnorm(runs)
        e <- scale1(stats::lm.fit(basis, raw)$residuals)
        basis <- cbind(basis, e)
        z <- if (f == 1L && p$rho > 0) p$rho * g + sqrt(1 - p$rho^2) * e else e
        factor_scores <- cbind(factor_scores, z)
    }
    if (!is.null(factor_scores)) colnames(factor_scores) <- paste0("F", seq_len(ncol(factor_scores)))
    planted <- cbind(planted, factor_scores)
    rownames(planted) <- run_ids
    # loadings
    tau <- sqrt(strength * median_variance)
    loadings <- if (p$factors > 0L) matrix(stats::rnorm(n * p$factors), n, p$factors) * tau
    offsets <- NULL
    if (p$batches > 0L) {
        offsets <- matrix(stats::rnorm(n * p$batches), n, p$batches) * sqrt(3 * median_variance)
        offsets <- offsets - rowMeans(offsets)
    }
    # residuals
    noise <- if (p$dup) {
        units <- max(unit)
        biological <- matrix(stats::rnorm(n * units), n, units) * sqrt(0.75 * s2)
        technical <- matrix(stats::rnorm(n * runs), n, runs) * sqrt(0.25 * s2)
        biological[, unit, drop = FALSE] + technical
    } else matrix(stats::rnorm(n * runs), n, runs) * sqrt(s2)
    # the uniforms of the masks, after the values (5.2); the batch-block draw W after U and V (B04)
    U <- matrix(stats::runif(n * runs), n, runs)
    V <- matrix(stats::runif(n * runs), n, runs)
    W <- if (p$batches > 0L) matrix(stats::runif(n * p$batches), n, p$batches)
    unwanted <- matrix(0, n, runs)
    if (p$factors > 0L) unwanted <- unwanted + loadings %*% t(factor_scores)
    if (p$batches > 0L) unwanted <- unwanted + offsets[, batch, drop = FALSE]
    y <- mu + outer(beta, as.numeric(group == "treated")) + unwanted + noise
    ids <- switch(p$omic, prot = sprintf("PG_%05d", seq_len(n)),
        lipid = c(sprintf("A1_%04d", seq_len(n / 2)), sprintf("A2_%04d", seq_len(n / 2))), metab = sprintf("M_%04d", seq_len(n)))
    assay <- switch(p$omic, lipid = rep(c("pos", "neg"), each = n / 2), rep("pos", n))
    dimnames(y) <- dimnames(unwanted) <- dimnames(U) <- dimnames(V) <- list(ids, run_ids)
    design <- if (identical(p$omic, "prot")) {
        data.frame(Run = run_ids, group = group, replicates = seq_len(runs), tech_rep_group = paste0("u", unit),
            stringsAsFactors = FALSE)
    } else data.frame(Sample_ID = run_ids, group = group, extract = paste0("e", seq_len(runs)), stringsAsFactors = FALSE)
    list(cell = cellname, p = p, r = r, y = y, unwanted = unwanted, beta = stats::setNames(beta, ids), U = U, V = V,
        W = W, batch = batch, group = group, unit = unit, runs = run_ids, Z = planted, design = design, assay = stats::setNames(assay, ids),
        strength = strength)
}

# ---- verbatim: sim-k185.R:202-234 at 432f703b ----
solveLod <- function(values, rate, shift = 0) {
    stats::uniroot(function(l) mean(stats::plogis(l + shift - values)) - rate, c(min(values) - 60, max(values) + 60),
        tol = 1e-10)$root
}
maskMatrix <- function(sim, mask) {
    y <- sim$y
    if (identical(mask, "none")) return(matrix(FALSE, nrow(y), ncol(y), dimnames = dimnames(y)))
    rate <- maskRate(mask)
    switch(maskMech(mask),
        vd = sim$U < stats::plogis(solveLod(y, rate) - y),
        mcar = sim$V < rate,
        bb = {
            block_mean <- vapply(seq_len(max(sim$batch)), function(b) rowMeans(y[, sim$batch == b, drop = FALSE]), numeric(nrow(y)))
            lod <- solveLod(block_mean, rate)
            block_missing <- sim$W < stats::plogis(lod - block_mean)
            out <- block_missing[, sim$batch, drop = FALSE]
            dimnames(out) <- dimnames(y)
            out
        },
        ra = {
            z1 <- sim$Z[, "F1"]
            top <- which.max(z1)
            shift_of <- function(c) matrix(c * z1, nrow(y), ncol(y), byrow = TRUE)
            share_top <- function(c) {
                lod <- solveLod(y, rate, shift_of(c))
                mean(stats::plogis(lod + c * z1[[top]] - y[, top]))
            }
            upper <- 1
            while (share_top(upper) < 2 * rate && upper < 1e3) upper <- upper * 2
            c_value <- stats::uniroot(function(c) share_top(c) - 2 * rate, c(0, upper), tol = 1e-10)$root
            sim$U < stats::plogis(solveLod(y, rate, shift_of(c_value)) + shift_of(c_value) - y)
        })
}

# ---- verbatim: sim-k185.R:237-244 at 432f703b ----
protObject <- function(mat, design) {
    quant <- as.data.frame(c(list(Protein.Group = rownames(mat)), as.data.frame(mat, check.names = FALSE)),
        check.names = FALSE, stringsAsFactors = FALSE)
    methods::new("ProteinQuantitativeData", protein_quant_table = quant, protein_id_column = "Protein.Group",
        protein_id_table = data.frame(Protein.Group = rownames(mat), stringsAsFactors = FALSE), design_matrix = design,
        sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
        args = list(globalParameters = list(workflow_type = "LFQ")))
}

# ---- verbatim: sim-k185.R:429-465 at 432f703b ----
CONTRAST <- "grouptreated-groupcontrol"
qMethod <- function(p) {
    p <- p[!is.na(p) & is.finite(p)]
    if (!length(p) || all(p == 0) || all(p == 1) || length(p) < 3L) return("BH")
    if (inherits(tryCatch(qvalue::qvalue(p), error = function(e) e), "error")) "BH" else "qvalue"
}
daTable <- function(sim, mat = NULL, object = NULL, ruv_factors = NULL) {
    omic <- sim$p$omic
    before <- warning_count$n
    if (identical(omic, "prot")) {
        object <- object %||% protObject(mat, sim$design)
        value <- quiet(differentialAbundanceAnalysis(object, contrasts_tbl = data.frame(contrasts = paste0("treated_vs_control=", CONTRAST),
            stringsAsFactors = FALSE), formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0,
            eBayes_trend = TRUE, eBayes_robust = TRUE, args_row_id = "Protein.Group", qvalue_column = "fdr_qvalue",
            raw_pvalue_column = "raw_pvalue"))
        tab <- value$contrasts_results$results[[1L]]
        if (is.null(tab) || !nrow(tab)) stop("differentialAbundanceAnalysis returned no table")
        out <- data.frame(id = rownames(tab), p = tab$raw_pvalue, q = tab$fdr_qvalue, stringsAsFactors = FALSE)
        attr(out, "qmethod") <- qMethod(out$p)
        attr(out, "failed") <- 0L
    } else {
        object <- object %||% omicObject(omic, mat, sim$assay, sim$design)
        contrasts <- data.frame(contrasts = CONTRAST, friendly_names = "treated_vs_control", stringsAsFactors = FALSE)
        value <- quiet(if (identical(omic, "lipid")) {
            runLipidsDA(object, contrasts, "~ 0 + group", 0.05, 0, TRUE, TRUE, ruv_factors = ruv_factors)
        } else runMetabolitesDA(object, contrasts, "~ 0 + group", 0.05, 0, TRUE, TRUE, ruv_factors = ruv_factors))
        long <- as.data.frame(value[[if (identical(omic, "lipid")) "da_lipids_long" else "da_metabolites_long"]])
        id_column <- if (identical(omic, "lipid")) "lipid_id" else "metabolite_id"
        out <- data.frame(id = as.character(long[[id_column]]), p = long$P.Value, q = long$fdr_qvalue, assay = long$assay,
            stringsAsFactors = FALSE)
        methods <- vapply(split(out$p, out$assay), qMethod, character(1L))
        attr(out, "qmethod") <- if (!length(methods)) NA_character_ else if (length(unique(methods)) == 1L) methods[[1L]] else "mixed"
        attr(out, "failed") <- length(value$failed_assays)
    }
    attr(out, "warnings") <- warning_count$n - before
    out
}

# ---- KIT-187 S2 (not copied) ---------------------------------------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x
KIT185_COMMIT <- "432f703b"
KIT185_LINES <- c("43-52", "71", "72-76", "90-91", "101-108", "109-199", "202-234", "237-244", "429-465")
CELLS <- list(A = cell(201L, "prot", 3000L, 2L, 3L), B = cell(202L, "prot", 3000L, 4L, 3L),
    C = cell(203L, "prot", 3000L, 5L, 3L), D = cell(204L, "prot", 3000L, 3L, 4L),
    E = cell(205L, "prot", 3000L, 2L, 6L), K = cell(206L, "prot", 3000L, 2L, 4L))
MASKS <- c("vd5", "vd10", "vd25")
MIN_REPS <- 2
MIN_GROUPS <- 2
PERCENTILE <- 1

k187Abort <- function(class, message) stop(structure(class = c(class, "error", "condition"),
    list(message = message, call = NULL)))

# Cell E: drawn as 2 groups of 6; runs 1 to 3 of `control` are dropped before the mask is drawn.
dropControlRuns <- function(sim) {
    keep <- setdiff(seq_along(sim$runs), which(sim$group == "control")[1:3])
    for (field in c("y", "U", "V")) sim[[field]] <- sim[[field]][, keep, drop = FALSE]
    sim$design <- sim$design[keep, , drop = FALSE]
    rownames(sim$design) <- NULL
    sim$group <- sim$group[keep]
    sim$runs <- sim$runs[keep]
    sim
}

# The filter as the tab runs it: updateMissingValueParameters, then the S4 method.
filterKept <- function(mat, design) {
    object <- protObject(mat, design)
    env <- new.env()
    env$config_list <- list()
    object <- quiet(updateMissingValueParameters(object, min_reps_per_group = MIN_REPS, min_groups = MIN_GROUPS, env = env))
    filtered <- quiet(removeRowsWithMissingValuesPercent(theObject = object, proteins_intensity_cutoff_percentile = PERCENTILE))
    as.character(filtered@protein_quant_table$Protein.Group)
}

# Computed here, not through the package: a value is present when finite and strictly above the 1st percentile
# (type 7) of every finite value of the masked matrix; a group passes with at most n - min(n, m) values missing and a
# protein is kept with at least g passing groups. For equal groups (cell A) this is the rule the tab states
# (falsifier 3); for cell E it is the per-group rule of D1 (b).
perGroupKept <- function(masked, group, m = MIN_REPS, g = MIN_GROUPS) {
    threshold <- stats::quantile(masked[is.finite(masked)], probs = PERCENTILE / 100, type = 7, names = FALSE)
    present <- is.finite(masked) & masked > threshold
    passing <- integer(nrow(masked))
    for (level in unique(group)) {
        columns <- group == level
        n <- sum(columns)
        passing <- passing + (rowSums(!present[, columns, drop = FALSE]) <= n - min(n, m))
    }
    rownames(masked)[passing >= g]
}

sideReplicate <- function(cellname, r, side, fault) {
    started <- proc.time()[[3L]]
    sim <- generateReplicate(cellname, r, strength = 0)
    if (identical(cellname, "E")) sim <- dropControlRuns(sim)
    ids <- rownames(sim$y)
    out <- list()
    for (mask in MASKS) {
        masked <- sim$y
        masked[maskMatrix(sim, mask)] <- NA_real_
        kept <- filterKept(masked, sim$design)
        if (fault && identical(side, "head") && identical(cellname, "K") && r == 1L && identical(mask, MASKS[[1L]]))
            kept <- kept[-1L]
        da <- daTable(sim, mat = masked[kept, , drop = FALSE])
        entry <- list(kept = match(kept, ids), calls = match(da$id[!is.na(da$q) & da$q < 0.05], ids),
            da_sha256 = digest::digest(da[c("id", "p", "q")], algo = "sha256"), qmethod = attr(da, "qmethod"))
        if (cellname %in% c("A", "E")) entry$per_group <- match(perGroupKept(masked, sim$group), ids)
        out[[mask]] <- entry
    }
    list(cell = cellname, r = r, masks = out, planted = which(sim$beta != 0), n = length(ids),
        seconds = proc.time()[[3L]] - started)
}

fixtureKept <- function(tree) {
    e2e <- file.path(tree, "tests", "testdata", "e2e")
    read <- function(...) utils::read.delim(file.path(e2e, ...), check.names = FALSE, stringsAsFactors = FALSE)
    one <- function(name, table, id_column, prefix, design_file) {
        design <- read(design_file)
        columns <- grep(prefix, names(table), value = TRUE)
        mat <- as.matrix(table[, columns, drop = FALSE])
        dimnames(mat) <- list(as.character(table[[id_column]]), sub(prefix, "", columns))
        mat <- log2(mat[, design$sample, drop = FALSE])
        kept <- filterKept(mat, data.frame(Run = design$sample, group = design$group, replicates = seq_len(nrow(design)),
            tech_rep_group = design$sample, stringsAsFactors = FALSE))
        data.frame(fixture = name, proteins = nrow(mat), kept = paste(kept, collapse = ","), n_kept = length(kept),
            stringsAsFactors = FALSE)
    }
    rbind(one("prot_lfq/proteinGroups.txt", read("prot_lfq", "proteinGroups.txt"), "Protein.IDs", "^LFQ\\.intensity\\.",
            file.path("prot_lfq", "design.tsv")),
        one("prot_tmt/pd_tmt_peptides.tsv", read("prot_tmt", "pd_tmt_peptides.tsv"), "Master.Protein.Accessions",
            "^Abundance: F[0-9]+: [0-9A-Z]+, ", file.path("prot_tmt", "design.tsv")))
}

runSide <- function(tree, side, work, cells, replicates, fault, workers) {
    suppressPackageStartupMessages(pkgload::load_all(tree, export_all = TRUE, helpers = FALSE, quiet = TRUE))
    logger::log_threshold(logger::FATAL)
    tasks <- expand.grid(r = seq_len(replicates), cell = cells, stringsAsFactors = FALSE)
    results <- parallel::mclapply(seq_len(nrow(tasks)), function(i) tryCatch(
        sideReplicate(tasks$cell[[i]], tasks$r[[i]], side, fault),
        error = function(e) list(cell = tasks$cell[[i]], r = tasks$r[[i]], error = conditionMessage(e))),
        mc.cores = workers, mc.preschedule = FALSE)
    fixtures <- tryCatch(fixtureKept(tree), error = function(e) data.frame(fixture = "error", proteins = NA,
        kept = conditionMessage(e), n_kept = NA, stringsAsFactors = FALSE))
    fixtures$side <- side
    commit <- system2("git", c("-C", shQuote(tree), "rev-parse", "HEAD"), stdout = TRUE)
    saveRDS(list(side = side, commit = commit, results = results, fixtures = fixtures, warnings = warning_count$n),
        file.path(work, paste0(side, ".rds")))
    errors <- Filter(function(x) inherits(x, "try-error") || !is.null(x$error), results)
    if (length(errors)) {
        cat(sprintf("side %s: %d replicate(s) failed; first: %s\n", side, length(errors),
            if (inherits(errors[[1L]], "try-error")) as.character(errors[[1L]]) else errors[[1L]]$error))
        quit(status = 1L)
    }
    cat(sprintf("side %s finished: %d replicates, %d warnings\n", side, length(results), warning_count$n))
}

# ---- the side process ----------------------------------------------------------------------------------------------
args <- commandArgs(TRUE)
if (length(args) && identical(args[[1L]], "--side")) {
    runSide(args[[2L]], args[[3L]], args[[4L]], strsplit(args[[5L]], ",", fixed = TRUE)[[1L]], as.integer(args[[6L]]),
        identical(args[[7L]], "1"), as.integer(args[[8L]]))
    quit(status = 0L)
}

# ---- the controller ------------------------------------------------------------------------------------------------
stopifnot(length(args) == 2L, grepl("^([0-9a-f]{40}|self)$", args[[1L]]))
started_at <- Sys.time()
lifecycle <- identical(Sys.getenv("K187_LIFECYCLE"), "1")
fault <- identical(Sys.getenv("K187_FAULT"), "1")
stop_rule <- !lifecycle || fault
cells <- strsplit(Sys.getenv("K187_CELLS", paste(names(CELLS), collapse = ",")), ",", fixed = TRUE)[[1L]]
stopifnot(length(cells) > 0L, all(cells %in% names(CELLS)))
cells <- names(CELLS)[names(CELLS) %in% cells]
replicates <- if (lifecycle) 2L else 100L
workers <- as.integer(Sys.getenv("K187_WORKERS", "2"))
out <- args[[2L]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))
head_root <- normalizePath(".")
head_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
base_root <- head_root
if (!identical(args[[1L]], "self")) {
    base_root <- file.path(tempdir(), "k187-s2-base")
    invisible(system2("git", c("worktree", "add", "--detach", shQuote(base_root), args[[1L]]), stdout = TRUE, stderr = TRUE))
    stopifnot(dir.exists(file.path(base_root, "R")))
}
base_commit <- system2("git", c("-C", shQuote(base_root), "rev-parse", "HEAD"), stdout = TRUE)
cat("KIT-187 S2: head", head_commit, "base", base_commit, "lifecycle", lifecycle, "fault", fault, "cells",
    paste(cells, collapse = ","), "R", replicates, "out", out, "\n")
# Native code compiled once per tree before the two sides start at once (one tree may serve both).
for (tree in unique(c(base_root, head_root))) pkgbuild::compile_dll(tree, quiet = TRUE)
work <- file.path(tempdir(), "k187-s2-work")
dir.create(work, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")
sides <- list(base = base_root, head = head_root)
jobs <- lapply(names(sides), function(side) parallel::mcparallel(system2(rscript, c(shQuote(script_file), "--side",
    shQuote(sides[[side]]), side, shQuote(work), paste(cells, collapse = ","), replicates, if (fault) "1" else "0", workers),
    stdout = file.path(work, paste0(side, ".log")), stderr = file.path(work, paste0(side, ".log")))))
status <- unlist(parallel::mccollect(jobs))
for (side in names(sides)) {
    log <- file.path(work, paste0(side, ".log"))
    if (file.exists(log)) cat(paste0("[", side, "] ", utils::tail(readLines(log, warn = FALSE), 20L)), sep = "\n")
}
wall_sides <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))
if (any(status != 0L) || !all(file.exists(file.path(work, c("base.rds", "head.rds")))))
    k187Abort("k187_instrument", sprintf("a side failed (exit status %s)", paste(status, collapse = ", ")))
data <- lapply(c(base = "base", head = "head"), function(side) readRDS(file.path(work, paste0(side, ".rds"))))
byCell <- function(results, cellname) {
    hits <- Filter(function(x) identical(x$cell, cellname), results)
    hits[order(vapply(hits, function(x) x$r, integer(1L)))]
}

summarise <- function(x) {
    x <- x[!is.na(x)]
    n <- length(x)
    m <- if (n) mean(x) else NA_real_
    half <- if (n > 1L) 1.96 * stats::sd(x) / sqrt(n) else NA_real_
    c(mean = m, lower = m - half, upper = m + half, R = n)
}
cell_rows <- list()
falsifier_rows <- list()
seconds <- list()
stopped_at <- NULL
addRows <- function(cellname, mask, side, metrics) for (metric in names(metrics)) {
    s <- summarise(metrics[[metric]])
    cell_rows[[length(cell_rows) + 1L]] <<- data.frame(cell = cellname, mask = mask, side = side, metric = metric,
        mean = s[["mean"]], lower = s[["lower"]], upper = s[["upper"]], R = s[["R"]], stringsAsFactors = FALSE)
}
addCount <- function(cellname, mask, side, metric, count) cell_rows[[length(cell_rows) + 1L]] <<- data.frame(cell = cellname,
    mask = mask, side = side, metric = metric, mean = count, lower = NA_real_, upper = NA_real_, R = replicates,
    stringsAsFactors = FALSE)
for (cellname in cells) {
    b_all <- byCell(data$base$results, cellname)
    h_all <- byCell(data$head$results, cellname)
    seconds[[cellname]] <- c(base = sum(vapply(b_all, function(x) x$seconds, numeric(1L))),
        head = sum(vapply(h_all, function(x) x$seconds, numeric(1L))))
    stopifnot(length(b_all) == replicates, length(h_all) == replicates)
    for (mask in MASKS) {
        m <- list()
        for (side in c("base", "head")) for (name in c("kept", "calls", "power", "fdp")) m[[side]][[name]] <- numeric(replicates)
        delta <- list(added = numeric(replicates), removed = numeric(replicates))
        added_share <- added_fdp <- rep(NA_real_, replicates)
        f1 <- f2 <- f3 <- 0L
        fate <- numeric(replicates)
        bh <- c(base = 0L, head = 0L)
        for (i in seq_len(replicates)) {
            b <- b_all[[i]]$masks[[mask]]
            h <- h_all[[i]]$masks[[mask]]
            planted <- h_all[[i]]$planted
            for (side in c("base", "head")) {
                x <- if (side == "base") b else h
                m[[side]]$kept[[i]] <- length(x$kept)
                m[[side]]$calls[[i]] <- length(x$calls)
                m[[side]]$power[[i]] <- length(intersect(x$calls, planted)) / length(planted)
                m[[side]]$fdp[[i]] <- if (length(x$calls)) length(setdiff(x$calls, planted)) / length(x$calls) else 0
                if (identical(x$qmethod, "BH")) bh[[side]] <- bh[[side]] + 1L
            }
            added <- setdiff(h$kept, b$kept)
            removed <- setdiff(b$kept, h$kept)
            delta$added[[i]] <- length(added)
            delta$removed[[i]] <- length(removed)
            if (length(added)) {
                called <- intersect(added, h$calls)
                added_share[[i]] <- length(called) / length(added)
                added_fdp[[i]] <- if (length(called)) length(setdiff(called, planted)) / length(called) else 0
            }
            if (length(removed)) f1 <- f1 + 1L
            if (cellname %in% c("K", "E") && (!identical(b$kept, h$kept) || !identical(b$da_sha256, h$da_sha256))) f2 <- f2 + 1L
            if (identical(cellname, "A") && length(union(setdiff(h$kept, h$per_group), setdiff(h$per_group, h$kept)))) f3 <- f3 + 1L
            if (identical(cellname, "E")) fate[[i]] <- length(union(setdiff(h$kept, h$per_group), setdiff(h$per_group, h$kept)))
        }
        for (side in c("base", "head")) {
            addRows(cellname, mask, side, m[[side]])
            addCount(cellname, mask, side, "bh_fallback_replicates", bh[[side]])
        }
        addRows(cellname, mask, "head-base", c(delta, list(kept = m$head$kept - m$base$kept, calls = m$head$calls - m$base$calls,
            power = m$head$power - m$base$power, fdp = m$head$fdp - m$base$fdp)))
        addRows(cellname, mask, "head", list(added_share_called = added_share, added_fdp = added_fdp))
        if (identical(cellname, "E")) addRows(cellname, mask, "head", list(per_group_rule_fate_changed = fate))
        falsifier_rows[[length(falsifier_rows) + 1L]] <- data.frame(falsifier = 1L, cell = cellname, mask = mask, count = f1)
        if (cellname %in% c("K", "E"))
            falsifier_rows[[length(falsifier_rows) + 1L]] <- data.frame(falsifier = 2L, cell = cellname, mask = mask, count = f2)
        if (identical(cellname, "A"))
            falsifier_rows[[length(falsifier_rows) + 1L]] <- data.frame(falsifier = 3L, cell = cellname, mask = mask, count = f3)
    }
    triggered <- Filter(function(row) identical(row$cell, cellname) && row$count > 0L, falsifier_rows)
    if (stop_rule && length(triggered)) {
        stopped_at <- list(cell = cellname, falsifiers = unique(vapply(triggered, function(row) row$falsifier, integer(1L))))
        break
    }
}

writeRows <- function(rows, name) utils::write.table(do.call(rbind, rows), file.path(out, name), sep = "\t", quote = FALSE,
    row.names = FALSE, na = "NA")
writeRows(cell_rows, "cells.tsv")
writeRows(falsifier_rows, "falsifiers.tsv")
utils::write.table(rbind(data$base$fixtures, data$head$fixtures), file.path(out, "fixtures.tsv"), sep = "\t", quote = FALSE,
    row.names = FALSE, na = "NA")
sha <- function(f) strsplit(system2("sha256sum", shQuote(f), stdout = TRUE), " ")[[1L]][[1L]]
provenance <- list(ticket = "KIT-187 S2", head_commit = head_commit, base_commit = base_commit, base_arg = args[[1L]],
    kit185_source = list(commit = KIT185_COMMIT, file = "docs/qa/workflow-execution/kit/kit-185/sim-k185.R", lines = KIT185_LINES),
    lifecycle = lifecycle, fault = fault, stop_rule = stop_rule, cells = cells, masks = MASKS, replicates = replicates,
    minimums = c(min_reps_per_group = MIN_REPS, min_groups = MIN_GROUPS), percentile = PERCENTILE, strength = 0,
    cell_index = vapply(CELLS[cells], `[[`, integer(1L), "idx"),
    seed_rule = "set.seed(185000000 + 10000 * idx + r, Mersenne-Twister, Inversion, Rejection); idx = 200 + i, so 187000000 + 10000 i + r",
    workers_per_side = workers, side_exit_status = as.list(stats::setNames(status, names(sides))),
    side_warnings = list(base = data$base$warnings, head = data$head$warnings),
    seconds_per_cell = lapply(seconds, as.list), wall_seconds_sides = wall_sides,
    wall_seconds = as.numeric(difftime(Sys.time(), started_at, units = "secs")), stopped_at = stopped_at,
    r_version = R.version.string, blas = extSoftVersion()[["BLAS"]], lapack = La_library(),
    runner_image = paste(Sys.getenv("ImageOS"), Sys.getenv("ImageVersion")),
    processor = tryCatch(sub("^model name\\s*:\\s*", "", grep("^model name", readLines("/proc/cpuinfo"), value = TRUE)[[1L]]),
        error = function(e) NA),
    script_sha256 = sha(script_file), started = format(started_at, "%Y-%m-%dT%H:%M:%S%z"),
    finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
jsonlite::write_json(provenance, file.path(out, "provenance.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat(sprintf("KIT-187 S2 written: %d cell rows, %d falsifier rows, wall %.1f s\n", length(cell_rows), length(falsifier_rows),
    provenance$wall_seconds))
if (!is.null(stopped_at))
    k187Abort("k187_falsifier", sprintf("falsifier %s triggered in cell %s", paste(stopped_at$falsifiers, collapse = ","),
        stopped_at$cell))
