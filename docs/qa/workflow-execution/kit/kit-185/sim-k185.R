# KIT-185 (VAL-56, VAL-67): RUV-III-C and missing values, the measurement. Evidence code under docs/qa: no package
# code. Runs on a GitHub-hosted runner with MultiScholaR installed from the bound commit (run-tests `scripts`).
#
# Usage (from the package root):
#   Rscript docs/qa/workflow-execution/kit/kit-185/sim-k185.R <part> <cells> <first> <last> <out dir>
#   part   A (part A, deciding and reported cells), B (part B, the routes choosing k and controls),
#          S (part A sensitivity panel: variants k1, k4, s1, s10, loess), R (run-aligned masks, D5)
#   cells  one cell or a comma list (P0,P0r,P1,P2,P3,P4,P5,L1,L2,M1,M2,P0lo,P0hi,P4b,L1b,M1b)
#   first, last  replicate range (S: 1 to 400; lifecycle: 401 to 403)
# Environment: K185_LIFECYCLE=1 (blind lifecycle: the masked arms' estimands are not written, only completion, timing,
#   mechanics counts and the SHA-256 of the rows; the checks E1 to E4, E6, E7 and gates G1, G2 are run and written to
#   checks-<job>.tsv), K185_FAULT=1 (E7: one deliberate error in the first unit), K185_WORKERS (default 4),
#   K185_PRESENCE (comma list of cells that take KIT-124 S4's presence rule after the shipped filter, G2; in lifecycle
#   the lipidomics and metabolomics cells are run both ways, variant "presence"), K185_VARIANTS (part S: comma list).
#
# Every random draw of a replicate is made in generateReplicate() before any package function is called; seeds:
# set.seed(185000000 + 10000 c + r, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection").
# The design, estimands and reading rules are KIT-185 revision 1.2, sections 5 and 7, and the declaration.

args <- commandArgs(TRUE)
if (length(args) != 5L) stop("usage: Rscript sim-k185.R <part> <cells> <first> <last> <out dir>")
part <- args[[1L]]
cells_wanted <- strsplit(args[[2L]], ",", fixed = TRUE)[[1L]]
first_rep <- as.integer(args[[3L]])
last_rep <- as.integer(args[[4L]])
out_dir <- args[[5L]]
stopifnot(part %in% c("A", "B", "S", "R"), !is.na(first_rep), !is.na(last_rep), first_rep <= last_rep)
lifecycle <- identical(Sys.getenv("K185_LIFECYCLE"), "1")
fault <- identical(Sys.getenv("K185_FAULT"), "1")
workers <- as.integer(Sys.getenv("K185_WORKERS", "4"))
presence_cells <- strsplit(Sys.getenv("K185_PRESENCE", ""), ",", fixed = TRUE)[[1L]]
started_at <- Sys.time()
script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))
here <- dirname(script_file)
source(file.path(here, "ruviiic-ref.R"))

suppressPackageStartupMessages(library(MultiScholaR))
attach(asNamespace("MultiScholaR"), name = "MultiScholaR:namespace", warn.conflicts = FALSE)
logger::log_threshold(logger::FATAL)
invisible(RUVIIIC::omp_set_num_threads(1L))

# ---- quiet calls, warnings counted ---------------------------------------------------------------------------------
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
refusal_classes <- c("prot_ruv_insufficient_controls", "omics_ruv_insufficient_controls", "prot_ruv_empty_result")
partb_refusal_classes <- c(refusal_classes, "prot_ruv_factors_unavailable", "ruv_factor_count_unavailable", "ruv_k_refused",
    "ruv_k_above_limit", "omics_ruv_no_residual_df", "omics_ruv_factors_refused")
attempt <- function(expr, refusals = refusal_classes) {
    t0 <- proc.time()[[3L]]
    out <- tryCatch(list(value = expr, status = "ok", class = NA_character_), error = function(e) {
        cls <- intersect(class(e), refusals)
        sweep <- grepl("No valid percentage found", conditionMessage(e), fixed = TRUE)
        if (inherits(e, "k185_not_applied")) return(list(value = NULL, status = "not_applied", class = "k185_not_applied"))
        list(value = NULL, status = if (length(cls) || sweep) "refused" else "crashed",
            class = if (length(cls)) cls[[1L]] else if (sweep) "sweep_no_valid_percentage" else
                substr(gsub("[\r\n\t]+", " ", paste(class(e)[[1L]], conditionMessage(e), sep = ": ")), 1L, 300L))
    })
    out$seconds <- proc.time()[[3L]] - t0
    out
}

# ---- cells (KIT-185 5.2) -------------------------------------------------------------------------------------------
GROUPS <- c("control", "treated", "g3", "g4", "g5")
cell <- function(idx, omic, features, groups, per, k = 2L, level = 20, rho = 0, dup = FALSE, batches = 0L,
                 batch_runs = 0L, factors = 2L, assays = 1L) {
    list(idx = idx, omic = omic, features = features, groups = groups, per = per, k = k, level = level, rho = rho,
        dup = dup, batches = batches, batch_runs = batch_runs, factors = factors, assays = assays)
}
CELLS <- list(
    P0 = cell(1L, "prot", 3000L, 2L, 6L), P0r = cell(2L, "prot", 3000L, 2L, 6L, rho = 0.3),
    P1 = cell(3L, "prot", 6000L, 2L, 4L), P2 = cell(4L, "prot", 3000L, 2L, 3L), P3 = cell(5L, "prot", 6000L, 5L, 4L),
    P4 = cell(6L, "prot", 5000L, 5L, 10L), P5 = cell(7L, "prot", 3000L, 2L, 4L, dup = TRUE),
    L1 = cell(8L, "lipid", 500L, 2L, 32L, assays = 2L), L2 = cell(9L, "lipid", 500L, 2L, 6L, assays = 2L),
    M1 = cell(10L, "metab", 150L, 4L, 45L), M2 = cell(11L, "metab", 150L, 4L, 3L),
    P0lo = cell(12L, "prot", 3000L, 2L, 6L, level = 14), P0hi = cell(13L, "prot", 3000L, 2L, 6L, level = 26),
    P4b = cell(14L, "prot", 5000L, 5L, 10L, k = 4L, batches = 5L, factors = 0L),
    L1b = cell(15L, "lipid", 500L, 2L, 32L, k = 2L, batches = 2L, factors = 1L, assays = 2L),
    M1b = cell(16L, "metab", 150L, 4L, 45L, k = 15L, batches = 15L, factors = 1L))
stopifnot(all(cells_wanted %in% names(CELLS)))
D5_CELLS <- c("P4b", "L1b", "M1b")
REF_EVERY_20 <- c("P4", "L1", "M1", "P4b", "L1b", "M1b")
MED_VAR_PROT <- 4 * 0.25^2 / stats::qchisq(0.5, 4)
MED_VAR_OMIC <- 0.25

masksOf <- function(part, cellname) {
    if (identical(part, "B")) return(if (cellname %in% c("P0", "L2", "M2")) c("none", "vd2", "vd10", "vd50") else
        c("none", "vd2", "vd10"))
    if (identical(part, "S")) return(c("none", "vd2", "vd10", "mcar2"))
    if (identical(part, "R")) return(c("none", "ra2", "ra10"))
    if (cellname %in% D5_CELLS) return(c("none", "vd2", "vd10", "bb2", "bb10", "bb25"))
    c("none", "vd0.1", "vd0.5", "vd2", "vd10", "vd25", "vd50", "mcar0.5", "mcar2", "mcar10")
}
maskMech <- function(mask) if (identical(mask, "none")) "none" else sub("[0-9.]+$", "", mask)
maskRate <- function(mask) if (identical(mask, "none")) 0 else as.numeric(sub("^[a-z]+", "", mask)) / 100

# ---- data (KIT-185 5.2) --------------------------------------------------------------------------------------------
scale1 <- function(v) {
    v <- v - mean(v)
    v / stats::sd(v)
}
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

# ---- masks (KIT-185 5.2) -------------------------------------------------------------------------------------------
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

# ---- package objects -------------------------------------------------------------------------------------------------
protObject <- function(mat, design) {
    quant <- as.data.frame(c(list(Protein.Group = rownames(mat)), as.data.frame(mat, check.names = FALSE)),
        check.names = FALSE, stringsAsFactors = FALSE)
    methods::new("ProteinQuantitativeData", protein_quant_table = quant, protein_id_column = "Protein.Group",
        protein_id_table = data.frame(Protein.Group = rownames(mat), stringsAsFactors = FALSE), design_matrix = design,
        sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
        args = list(globalParameters = list(workflow_type = "LFQ")))
}
protMatrix <- function(object) {
    q <- object@protein_quant_table
    samples <- as.character(object@design_matrix[[object@sample_id]])
    m <- as.matrix(q[, samples, drop = FALSE])
    rownames(m) <- as.character(q[[object@protein_id_column]])
    m
}
protWithTable <- function(object, mat) {
    object@protein_quant_table <- as.data.frame(c(list(Protein.Group = rownames(mat)), as.data.frame(mat, check.names = FALSE)),
        check.names = FALSE, stringsAsFactors = FALSE)
    object
}
ID_COL <- "feature_id"
assayTable <- function(m) {
    table <- tibble::as_tibble(as.data.frame(m, check.names = FALSE))
    table <- tibble::add_column(table, !!ID_COL := rownames(m), .before = 1L)
    tibble::add_column(table, name = rownames(m), .after = 1L)
}
omicSlot <- function(omic) if (identical(omic, "lipid")) "lipid_data" else "metabolite_data"
omicObject <- function(omic, mat, assay, design) {
    tables <- lapply(split(rownames(mat), factor(assay[rownames(mat)], levels = unique(assay))), function(ids) assayTable(mat[ids, , drop = FALSE]))
    if (identical(omic, "lipid")) {
        methods::new("LipidomicsAssayData", lipid_data = tables, lipid_id_column = ID_COL, annotation_id_column = "name",
            design_matrix = design, sample_id = "Sample_ID", group_id = "group", technical_replicate_id = "extract", args = list())
    } else {
        methods::new("MetaboliteAssayData", metabolite_data = tables, metabolite_id_column = ID_COL, annotation_id_column = "name",
            design_matrix = design, sample_id = "Sample_ID", group_id = "group", technical_replicate_id = "extract", args = list())
    }
}
omicMatrix <- function(object, omic, runs) {
    tables <- methods::slot(object, omicSlot(omic))
    do.call(rbind, lapply(tables, function(t) {
        m <- as.matrix(as.data.frame(t[, runs], check.names = FALSE))
        rownames(m) <- as.character(t[[ID_COL]])
        m
    }))
}
omicWithMatrix <- function(object, omic, mat, assay) {
    tables <- lapply(split(rownames(mat), factor(assay[rownames(mat)], levels = unique(assay[rownames(mat)]))),
        function(ids) assayTable(mat[ids, , drop = FALSE]))
    methods::slot(object, omicSlot(omic)) <- tables
    object
}
replicateMatrix <- function(design, sample_column, block_column) {
    getRuvIIIReplicateMatrixHelper(design, !!rlang::sym(sample_column), !!rlang::sym(block_column))
}
memStore <- function() {
    saved <- new.env()
    list(saved = saved, checkpoint = function(value, label, from = NULL) {
        assign(label, value, envir = saved)
        label
    }, product = function(value, kind, label) {
        assign(label, value, envir = saved)
        label
    })
}

# ---- the route's filter, the input of the RUV step ---------------------------------------------------------------
filterArm <- function(sim, masked, use_presence) {
    p <- sim$p
    if (identical(p$omic, "prot")) {
        object <- protObject(masked, sim$design)
        env <- new.env()
        env$config_list <- list()
        object <- quiet(updateMissingValueParameters(object, min_reps_per_group = 2, min_groups = 2, env = env))
        filtered <- quiet(removeRowsWithMissingValuesPercent(theObject = object))
        ids <- as.character(filtered@protein_quant_table$Protein.Group)
        return(list(ids = ids, object = filtered))
    }
    raw <- omicObject(p$omic, 2^masked, sim$assay, sim$design)
    filtered <- quiet(if (identical(p$omic, "lipid")) {
        lipidIntensityFiltering(raw, lipids_intensity_cutoff_percentile = 10, lipids_proportion_of_samples_below_cutoff = 0.5)
    } else {
        metaboliteIntensityFiltering(raw, metabolites_intensity_cutoff_percentile = 10, metabolites_proportion_of_samples_below_cutoff = 0.5)
    })
    ids <- unlist(lapply(methods::slot(filtered, omicSlot(p$omic)), function(t) as.character(t[[ID_COL]])), use.names = FALSE)
    ids <- rownames(masked)[rownames(masked) %in% ids]
    if (use_presence) {
        seen <- !is.na(masked[ids, , drop = FALSE])
        ok <- Reduce(`&`, lapply(unique(sim$group), function(g) rowSums(seen[, sim$group == g, drop = FALSE]) >= 1L))
        ids <- ids[ok]
    }
    list(ids = ids, object = omicObject(p$omic, masked[ids, , drop = FALSE], sim$assay, sim$design))
}

replicateSets <- function(sim) {
    if (identical(sim$p$omic, "prot")) {
        blocks <- resolveProtRuvReplicateBlocks(sim$design, "Run", "group", "tech_rep_group", "group")
        M <- replicateMatrix(blocks$design, "Run", ".ruv_replicate_block")
        list(block = blocks$design$.ruv_replicate_block, M = M[sim$runs, , drop = FALSE])
    } else {
        M <- replicateMatrix(sim$design, "Sample_ID", "group")
        list(block = sim$group, M = M[sim$runs, , drop = FALSE])
    }
}
residualDimension <- function(obs, block) {
    sets <- vapply(unique(block), function(b) rowSums(obs[, block == b, drop = FALSE]) > 0, logical(nrow(obs)))
    if (is.null(dim(sets))) sets <- matrix(sets, nrow = nrow(obs))
    rowSums(obs) - rowSums(sets)
}

buildArm <- function(sim, mask, norm, use_presence, sets) {
    t0 <- proc.time()[[3L]]
    masked <- sim$y
    missing <- maskMatrix(sim, mask)
    masked[missing] <- NA_real_
    filtered <- filterArm(sim, masked, use_presence)
    ids <- filtered$ids
    object <- filtered$object
    if (identical(sim$p$omic, "prot")) {
        object <- quiet(normaliseBetweenSamples(object, normalisation_method = norm))
        Y <- protMatrix(object)[ids, sim$runs, drop = FALSE]
    } else {
        if (!identical(norm, "none")) object <- quiet(normaliseBetweenSamples(object, normalisation_method = norm))
        Y <- omicMatrix(object, sim$p$omic, sim$runs)[ids, , drop = FALSE]
    }
    obs <- is.finite(Y)
    Y[!obs] <- NA_real_
    run_share <- colMeans(!obs)
    block_share <- if (!is.null(sim$batch)) {
        mean(vapply(seq_len(max(sim$batch)), function(b) rowSums(obs[, sim$batch == b, drop = FALSE]) == 0L, logical(nrow(obs))))
    } else NA_real_
    list(mask = mask, ids = ids, object = object, filtered = filtered$object, Y = Y, obs = obs, rd = stats::setNames(residualDimension(obs, sets$block), ids),
        removed = nrow(sim$y) - length(ids), miss_share = mean(!obs), run_share = run_share, block_share = block_share,
        seconds = proc.time()[[3L]] - t0)
}

# ---- the corrections (KIT-185 5.3) ---------------------------------------------------------------------------------
allControls <- function(arm, omic) if (identical(omic, "prot")) stats::setNames(rep(TRUE, length(arm$ids)), arm$ids) else arm$ids
fitShipped <- function(sim, arm, k, ctrl = NULL) {
    omic <- sim$p$omic
    ctrl <- ctrl %||% allControls(arm, omic)
    if (identical(omic, "prot")) {
        corrected <- quiet(ruvIII_C_Varying(arm$object, ruv_grouping_variable = "group", ruv_number_k = k, ctrl = ctrl))
        corrected <- quiet(removeRowsWithMissingValuesPercent(theObject = corrected))
        return(protMatrix(corrected))
    }
    corrected <- quiet(ruvIII_C_Varying(arm$object, ruv_grouping_variable = "group", ruv_number_k = k, ctrl = ctrl))
    before <- methods::slot(arm$object, omicSlot(omic))
    after <- methods::slot(corrected, omicSlot(omic))
    if (!identical(names(after), names(before)) || any(vapply(names(before), function(a) identical(after[[a]], before[[a]]), logical(1L)))) {
        stop(structure(class = c("k185_not_applied", "error", "condition"),
            list(message = "the legacy RUV method returned an assay uncorrected or dropped it", call = NULL)))
    }
    omicMatrix(corrected, omic, sim$runs)
}
# the shipped compiled call the S4 method makes, before its post-processing (E1 (a) raw)
fitCompiledRaw <- function(sim, arm, k, sets, assay_ids = arm$ids) {
    Y <- t(arm$Y[assay_ids, sim$runs, drop = FALSE])
    out <- RUVIIIC::RUVIII_C_Varying(k = k, Y = Y, M = sets$M, toCorrect = colnames(Y), potentialControls = colnames(Y))
    t(out)
}
fitReference <- function(sim, arm, k, sets, directions, controls = NULL, raw = FALSE) {
    omic <- sim$p$omic
    groups <- if (identical(omic, "prot")) list(all = arm$ids) else split(arm$ids, factor(sim$assay[arm$ids], levels = unique(sim$assay[arm$ids])))
    pieces <- lapply(groups, function(assay_ids) {
        Y <- t(arm$Y[assay_ids, sim$runs, drop = FALSE])
        potential <- if (is.null(controls)) assay_ids else intersect(controls, assay_ids)
        out <- ruviiicRef(k, Y, sets$M, toCorrect = assay_ids, potentialControls = potential, directions = directions)
        m <- t(out)
        attr(m, "gap_ratio") <- attr(out, "gap_ratio")
        attr(m, "pattern_errors") <- attr(out, "pattern_errors")
        m
    })
    corrected <- do.call(rbind, pieces)
    gap <- unlist(lapply(pieces, attr, "gap_ratio"))
    errors <- sum(unlist(lapply(pieces, attr, "pattern_errors")))
    if (raw) return(structure(corrected, gap_ratio = gap, pattern_errors = errors))
    out <- postProcess(sim, arm, corrected)
    structure(out, gap_ratio = gap, pattern_errors = errors)
}
# the route's post-processing of a corrected matrix (as the S4 method and the route do to S's)
postProcess <- function(sim, arm, corrected, restore = identical(sim$p$omic, "prot")) {
    corrected <- corrected[, sim$runs, drop = FALSE]
    if (restore) corrected <- restoreRuvFeatureLevels(corrected, arm$Y[rownames(corrected), sim$runs, drop = FALSE])
    corrected <- corrected[rowSums(is.finite(corrected)) > 0, , drop = FALSE]
    if (identical(sim$p$omic, "prot")) {
        object <- quiet(removeRowsWithMissingValuesPercent(theObject = protWithTable(arm$object, corrected)))
        corrected <- protMatrix(object)
    }
    corrected
}

# ---- the differential test (KIT-185 5.3, part A: the legacy routes' test) ----------------------------------------
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

# ---- estimands (KIT-185 5.5) ---------------------------------------------------------------------------------------
patternKeys <- function(obs) do.call(paste0, as.data.frame(obs * 1L))
groupDesign <- function(group) {
    G <- stats::model.matrix(~ 0 + factor(group, levels = unique(group)))
    dimnames(G) <- NULL
    G
}
# per feature: ||P_j y*_j||^2 (P_j onto the planted factors on the observed runs with the groups projected out) and
# the residual sum of squares of y*_j on the groups
projectionTerms <- function(Ystar, obs, G, Z) {
    keys <- patternKeys(obs)
    U <- RSS <- rep(NA_real_, nrow(Ystar))
    for (key in unique(keys)) {
        rows <- which(keys == key)
        runs <- which(obs[rows[[1L]], ])
        if (length(runs) < 2L) next
        Gs <- G[runs, colSums(G[runs, , drop = FALSE]) > 0, drop = FALSE]
        qG <- qr(Gs)
        Yt <- t(Ystar[rows, runs, drop = FALSE])
        RSS[rows] <- colSums(qr.resid(qG, Yt)^2)
        qZ <- qr(qr.resid(qG, Z[runs, , drop = FALSE]), tol = 1e-7)
        U[rows] <- if (qZ$rank == 0L) 0 else colSums(qr.fitted(qZ, Yt, k = qZ$rank)^2)
    }
    list(U = U, RSS = RSS)
}
slopeOf <- function(Ystar, obs, group, beta) {
    Y <- Ystar
    Y[!obs] <- NA_real_
    m1 <- rowMeans(Y[, group == "control", drop = FALSE], na.rm = TRUE)
    m2 <- rowMeans(Y[, group == "treated", drop = FALSE], na.rm = TRUE)
    ok <- beta != 0 & is.finite(m1) & is.finite(m2)
    if (!any(ok)) return(NA_real_)
    sum(beta[ok] * (m2 - m1)[ok]) / sum(beta[ok]^2)
}
recoveryComplete <- function(Ybase, Ystar, Z) {
    if (nrow(Ybase) < 2L) return(list(rank = NA_integer_, rho = rep(NA_real_, ncol(Z))))
    D <- Ybase - Ystar
    D <- D - rowMeans(D)
    s <- svd(D, nu = 0L)
    if (!length(s$d) || max(s$d) == 0) return(list(rank = 0L, rho = rep(0, ncol(Z))))
    keep <- s$d > 1e-8 * max(s$d)
    rho <- stats::cancor(s$v[, keep, drop = FALSE], Z)$cor
    list(rank = sum(keep), rho = c(rho, rep(0, ncol(Z)))[seq_len(ncol(Z))])
}
recoveryIncomplete <- function(Ybase, Ystar, obs, Z) {
    if (!nrow(Ybase)) return(NA_real_)
    keys <- patternKeys(obs)
    share <- rep(NA_real_, nrow(Ybase))
    for (key in unique(keys)) {
        rows <- which(keys == key)
        runs <- which(obs[rows[[1L]], ])
        if (length(runs) < 2L) next
        Zr <- scale(Z[runs, , drop = FALSE], center = TRUE, scale = FALSE)
        qZ <- qr(Zr, tol = 1e-7)
        d <- t(Ybase[rows, runs, drop = FALSE] - Ystar[rows, runs, drop = FALSE])
        d <- d - matrix(colMeans(d), nrow(d), ncol(d), byrow = TRUE)
        denominator <- colSums(d^2)
        numerator <- if (qZ$rank == 0L) 0 else colSums(qr.fitted(qZ, d, k = qZ$rank)^2)
        share[rows] <- ifelse(denominator > 0, numerator / denominator, NA_real_)
    }
    stats::median(share, na.rm = TRUE)
}
daMetrics <- function(tab, ids, beta) {
    if (is.null(tab)) return(c(n_tested = NA, n_true = sum(beta[ids] != 0), n_called = NA, fdp = NA, power = NA, null_rate = NA))
    planted <- beta[ids] != 0
    m <- match(ids, tab$id)
    pv <- tab$p[m]
    qv <- tab$q[m]
    called <- !is.na(qv) & is.finite(qv) & qv < 0.05
    n_true <- sum(planted)
    nulls <- !planted & !is.na(pv) & is.finite(pv)
    c(n_tested = sum(!is.na(pv) & is.finite(pv)), n_true = n_true, n_called = sum(called),
        fdp = if (any(called)) sum(called & !planted) / sum(called) else 0,
        power = if (n_true > 0L) sum(called & planted) / n_true else NA_real_,
        null_rate = if (any(nulls)) mean(pv[nulls] < 0.05) else NA_real_)
}

COLUMNS <- c("cell", "part", "variant", "rep", "mask", "mech", "rate", "corr", "stratum", "status", "class", "k", "n_filt",
    "n_set", "n_comp", "removed", "miss_share", "miss_run_min", "miss_run_med", "miss_run_max", "n_rd_lt", "n_rd_eq",
    "n_rd_gt", "block_share", "rank", "rho", "rho_min", "A1", "inc_share", "L", "B", "RSS", "B2", "slope", "Cstar",
    "n_tested", "n_true", "n_called", "fdp", "power", "null_rate", "qmethod", "da_warn", "da_failed", "lost",
    "lost_planted", "lost_null", "k_chosen", "n_controls", "percentage", "seconds")
BLIND_KEEP <- c("cell", "part", "variant", "rep", "mask", "mech", "rate", "corr", "stratum", "status", "class", "k", "n_filt",
    "removed", "miss_share", "miss_run_min", "miss_run_med", "miss_run_max", "n_rd_lt", "n_rd_eq", "n_rd_gt", "block_share",
    "qmethod", "da_warn", "da_failed", "seconds")
newRow <- function(...) {
    row <- stats::setNames(rep(list(NA), length(COLUMNS)), COLUMNS)
    values <- lapply(list(...), function(v) if (is.null(v) || !length(v)) NA else if (length(v) > 1L) paste(v, collapse = ";") else v)
    row[names(values)] <- values
    row
}

# The estimand rows of one arm. `fits`: named list, each list(matrix, status, class, seconds, da (table), obs (the
# observed runs of its values), none = TRUE for the none-valued references S0, O0, N0, SB0...). `refs`: names of the
# corrections that define L(O) and L(N) for each correction.
scoreArm <- function(sim, arm, fits, F_set, F_filt, k, common, Z, G) {
    rows <- list()
    beta <- sim$beta
    rd <- arm$rd[F_filt]
    strata_rd <- list(rd_lt = F_filt[rd < k], rd_eq = F_filt[rd == k], rd_gt = F_filt[rd > k])
    comp <- F_set[rowSums(!arm$obs[F_set, , drop = FALSE]) == 0L]
    inc <- setdiff(F_set, comp)
    strata <- list(all = F_set, comp = comp, inc = inc)
    terms <- list()
    for (name in names(fits)) {
        fit <- fits[[name]]
        if (!identical(fit$status, "ok") || is.null(fit$matrix)) next
        obs <- fit$obs[F_set, , drop = FALSE]
        terms[[name]] <- projectionTerms(fit$matrix[F_set, sim$runs, drop = FALSE], obs, G, Z)
    }
    for (name in names(fits)) {
        fit <- fits[[name]]
        base <- c(common, list(corr = name, status = fit$status, class = fit$class, k = fit$k %||% k, k_chosen = fit$k_chosen %||% NA,
            n_controls = fit$n_controls %||% NA, percentage = fit$percentage %||% NA))
        ok <- identical(fit$status, "ok") && !is.null(fit$matrix)
        for (stratum in names(strata)) {
            set <- strata[[stratum]]
            row <- do.call(newRow, c(base, list(stratum = stratum, n_set = length(set), n_comp = length(comp),
                seconds = if (identical(stratum, "all")) fit$seconds else NA)))
            if (ok && length(set)) {
                idx <- match(set, F_set)
                Ostar <- terms[[fit$ref_O]]
                Nstar <- terms[[fit$ref_N]]
                if (!is.null(terms[[name]]) && !is.null(Ostar) && !is.null(Nstar)) {
                    L_mu <- sum(terms[[name]]$U[idx])
                    L_O <- sum(Ostar$U[idx])
                    L_N <- sum(Nstar$U[idx])
                    row$L <- L_mu
                    row$B <- (L_mu - L_O) / (L_N - L_O)
                    row$RSS <- sum(terms[[name]]$RSS[idx])
                    row$B2 <- row$RSS / sum(Ostar$RSS[idx])
                }
                obs <- fit$obs[set, , drop = FALSE]
                Ystar <- fit$matrix[set, sim$runs, drop = FALSE]
                row$slope <- slopeOf(Ystar, obs, sim$group, beta[set])
                O_slope <- slopeOf(fits[[fit$ref_O]]$matrix[set, sim$runs, drop = FALSE], obs, sim$group, beta[set])
                row$Cstar <- O_slope - row$slope
                if (!fit$name_kind %in% c("O", "N")) {
                    base_values <- fit$base[set, sim$runs, drop = FALSE]
                    if (identical(stratum, "comp")) {
                        recovery <- recoveryComplete(base_values, Ystar, Z)
                        row$rank <- recovery$rank
                        row$rho <- paste(formatC(recovery$rho, digits = 10, format = "f"), collapse = ";")
                        row$rho_min <- min(recovery$rho)
                        row$A1 <- 1 - row$rho_min^2
                    }
                    if (identical(stratum, "inc")) row$inc_share <- recoveryIncomplete(base_values, Ystar, obs, Z)
                }
            }
            rows[[length(rows) + 1L]] <- row
        }
        # (d) and (f): on F_filt and by residual dimension
        for (stratum in c("filt", names(strata_rd))) {
            set <- if (identical(stratum, "filt")) F_filt else strata_rd[[stratum]]
            row <- do.call(newRow, c(base, list(stratum = stratum, n_set = length(set), qmethod = attr(fit$da, "qmethod") %||% NA,
                da_warn = attr(fit$da, "warnings") %||% NA, da_failed = attr(fit$da, "failed") %||% NA)))
            if (!is.null(fit$da) && length(set)) {
                metrics <- daMetrics(fit$da, set, beta)
                row[names(metrics)] <- as.list(metrics)
            }
            if (ok && fit$name_kind %in% c("S", "C", "Z") && length(set)) {
                finite <- set %in% rownames(fit$matrix)[rowSums(is.finite(fit$matrix)) > 0]
                planted <- beta[set] != 0
                row$lost <- mean(!finite)
                row$lost_planted <- if (any(planted)) mean(!finite[planted]) else NA_real_
                row$lost_null <- if (any(!planted)) mean(!finite[!planted]) else NA_real_
            }
            rows[[length(rows) + 1L]] <- row
        }
    }
    rows
}

# ---- one replicate of part A (and S, R) ------------------------------------------------------------------------
checks_env <- new.env()
checks_env$rows <- list()
addCheck <- function(check, cellname, r, mask, value, threshold, pass, note = "") {
    checks_env$rows[[length(checks_env$rows) + 1L]] <- data.frame(check = check, cell = cellname, rep = r, mask = mask,
        value = signif(as.numeric(value), 6), threshold = threshold, pass = pass, note = note, stringsAsFactors = FALSE)
}
maxDiff <- function(a, b) {
    common <- intersect(rownames(a), rownames(b))
    fa <- rownames(a)[rowSums(is.finite(a)) > 0]
    fb <- rownames(b)[rowSums(is.finite(b)) > 0]
    same <- setequal(fa, fb) && identical(is.finite(a[common, , drop = FALSE]), is.finite(b[common, , drop = FALSE]))
    d <- abs(a[common, , drop = FALSE] - b[common, colnames(a), drop = FALSE])
    list(same_finite = same, max = if (any(is.finite(d))) max(d, na.rm = TRUE) else 0, share_1e8 = mean(d > 1e-8, na.rm = TRUE))
}

unitA <- function(cellname, r, part, variant) {
    t_unit <- proc.time()[[3L]]
    warning_count$n <- 0L
    p <- CELLS[[cellname]]
    k <- switch(variant, k1 = 1L, k4 = 4L, p$k)
    strength <- switch(variant, s1 = 0.5, s10 = 5, 1.5)
    norm <- if (identical(variant, "loess")) "cyclicloess" else "none"
    use_presence <- identical(variant, "presence") || cellname %in% presence_cells
    sim <- generateReplicate(cellname, r, strength)
    sets <- replicateSets(sim)
    G <- groupDesign(sim$group)
    Z <- sim$Z
    masks <- masksOf(part, cellname)
    want_ref <- lifecycle || (cellname %in% REF_EVERY_20 && r %% 20L == 0L && identical(part, "A"))
    want_e2 <- lifecycle || (r %% 20L == 0L)
    timing <- c(filter = 0, S = 0, C = 0, Z = 0, DA = 0, estimands = 0)
    arms <- list()
    for (mask in masks) {
        arms[[mask]] <- buildArm(sim, mask, norm, use_presence, sets)
        timing[["filter"]] <- timing[["filter"]] + arms[[mask]]$seconds
    }
    # corrections per arm
    for (mask in masks) {
        arm <- arms[[mask]]
        O <- arm$Y - sim$unwanted[arm$ids, sim$runs, drop = FALSE]
        S <- attempt(fitShipped(sim, arm, k))
        timing[["S"]] <- timing[["S"]] + S$seconds
        C <- list(status = "not_run", class = NA_character_, value = NULL, seconds = NA)
        if (identical(S$status, "ok") && (mask != "none" || want_e2)) {
            C <- attempt({
                if (fault && identical(mask, masks[[2L]]) && r == first_rep && identical(cellname, cells_wanted[[1L]]) &&
                    identical(variant, variants_of(part)[[1L]])) stop("deliberate fault (K185_FAULT, E7)")
                fitReference(sim, arm, k, sets, "complete")
            })
            timing[["C"]] <- timing[["C"]] + C$seconds
        }
        Zf <- list(status = "not_run", class = NA_character_, value = NULL, seconds = NA)
        if (want_ref && identical(S$status, "ok")) {
            Zf <- attempt(fitReference(sim, arm, k, sets, "zero_filled"))
            timing[["Z"]] <- timing[["Z"]] + Zf$seconds
        }
        arms[[mask]]$fit <- list(S = S, C = C, Z = Zf, O = O, N = arm$Y)
        # E1 (a): the reference with zero_filled equals the compiled call, before and after post-processing
        if (lifecycle && identical(S$status, "ok") && identical(Zf$status, "ok")) {
            raw_cpp <- tryCatch({
                if (identical(p$omic, "prot")) fitCompiledRaw(sim, arm, k, sets) else
                    do.call(rbind, lapply(split(arm$ids, factor(sim$assay[arm$ids], levels = unique(sim$assay))),
                        function(a) fitCompiledRaw(sim, arm, k, sets, a)))
            }, error = function(e) NULL)
            raw_ref <- fitReference(sim, arm, k, sets, "zero_filled", raw = TRUE)
            tolerance <- if (length(sim$runs) > 20L) 1e-6 else 1e-8
            if (!is.null(raw_cpp)) {
                d <- maxDiff(raw_cpp, raw_ref[rownames(raw_cpp), , drop = FALSE])
                gap <- attr(raw_ref, "gap_ratio")
                addCheck("E1a_raw", cellname, r, mask, d$max, tolerance, d$same_finite && d$max <= tolerance,
                    sprintf("share>1e-8 %.4g; gap ratio median %.4g max %.4g; pattern errors %d", d$share_1e8,
                        if (length(gap)) stats::median(gap) else NA, if (length(gap)) max(gap) else NA, attr(raw_ref, "pattern_errors")))
            } else addCheck("E1a_raw", cellname, r, mask, NA, tolerance, FALSE, "compiled call failed")
            d <- maxDiff(S$value, Zf$value)
            addCheck("E1a_post", cellname, r, mask, d$max, tolerance, d$same_finite && d$max <= tolerance, sprintf("share>1e-8 %.4g", d$share_1e8))
        }
    }
    # E2: on none, C equals S
    none <- arms[["none"]]
    if (want_e2 && identical(none$fit$C$status, "ok")) {
        tolerance <- if (length(sim$runs) > 20L) 1e-6 else 1e-8
        d <- maxDiff(none$fit$S$value, none$fit$C$value)
        addCheck("E2", cellname, r, "none", d$max, tolerance, d$same_finite && d$max <= tolerance, variant)
    }
    # the DA calls and the estimands
    t_da <- 0
    da_errors <- character()
    da <- function(mat) {
        t0 <- proc.time()[[3L]]
        out <- tryCatch(daTable(sim, mat), error = function(e) {
            da_errors <<- c(da_errors, substr(gsub("[\r\n\t]+", " ", conditionMessage(e)), 1L, 200L))
            NULL
        })
        seconds <- proc.time()[[3L]] - t0
        t_da <<- t_da + seconds
        if (!is.null(out)) attr(out, "seconds") <- seconds
        out
    }
    S_none <- none$fit$S$value
    complete_y <- sim$y[, sim$runs, drop = FALSE]
    rows <- list()
    for (mask in masks) {
        arm <- arms[[mask]]
        fit <- arm$fit
        t_est <- proc.time()[[3L]]
        common <- list(cell = cellname, part = part, variant = variant, rep = r, mask = mask, mech = maskMech(mask), rate = maskRate(mask),
            n_filt = length(arm$ids))
        rd <- arm$rd
        rows[[length(rows) + 1L]] <- do.call(newRow, c(common, list(corr = "mask", stratum = "all", status = "ok", k = k,
            removed = arm$removed, miss_share = arm$miss_share, miss_run_min = min(arm$run_share), miss_run_med = stats::median(arm$run_share),
            miss_run_max = max(arm$run_share), n_rd_lt = sum(rd < k), n_rd_eq = sum(rd == k), n_rd_gt = sum(rd > k),
            block_share = arm$block_share, seconds = arm$seconds)))
        F_filt <- arm$ids
        finiteRows <- function(m) if (is.null(m)) character() else rownames(m)[rowSums(is.finite(m)) > 0]
        F_set <- if (identical(fit$S$status, "ok")) {
            Reduce(intersect, list(F_filt, finiteRows(fit$S$value), if (identical(fit$C$status, "ok")) finiteRows(fit$C$value) else F_filt,
                finiteRows(S_none)))
        } else character()
        obs <- arm$obs
        entry <- function(kind, matrix, status, class, seconds, base, ref_O, ref_N, obs_matrix, da_rows = NULL) {
            list(name_kind = kind, matrix = matrix, status = status, class = class, seconds = seconds, base = base, ref_O = ref_O,
                ref_N = ref_N, obs = obs_matrix, da_rows = da_rows)
        }
        fits <- list(
            S = entry("S", fit$S$value, fit$S$status, fit$S$class, fit$S$seconds, arm$Y, "O", "N", obs),
            C = entry("C", fit$C$value, fit$C$status, fit$C$class, fit$C$seconds, arm$Y, "O", "N", obs),
            Z = entry("Z", fit$Z$value, fit$Z$status, fit$Z$class, fit$Z$seconds, arm$Y, "O", "N", obs),
            O = entry("O", fit$O, "ok", NA_character_, NA, arm$Y, "O", "N", obs),
            N = entry("N", fit$N, "ok", NA_character_, NA, arm$Y, "O", "N", obs))
        if (!identical(mask, "none")) {
            none_ids <- none$ids
            O_none <- none$fit$O
            all_obs <- matrix(TRUE, length(F_filt), length(sim$runs), dimnames = list(F_filt, sim$runs))
            fits$S0 <- entry("S", S_none, none$fit$S$status, none$fit$S$class, NA, none$Y, "O0", "N0", all_obs)
            fits$O0 <- entry("O", O_none, "ok", NA_character_, NA, none$Y, "O0", "N0", all_obs)
            fits$N0 <- entry("N", none$Y, "ok", NA_character_, NA, none$Y, "O0", "N0", all_obs)
        }
        if (fit$Z$status == "not_run") fits$Z <- NULL
        if (identical(mask, "none") && fit$C$status == "not_run") fits$C <- NULL
        # the DA calls: S, C, Z, O, N on the masked values; S0 and O0 on the none values of the rows of F_filt
        for (name in names(fits)) {
            f <- fits[[name]]
            if (!identical(fit$S$status, "ok")) break
            if (!identical(f$status, "ok") || is.null(f$matrix)) next
            if (identical(mask, "none") && name %in% c("C", "N")) next
            if (identical(name, "N0")) next
            rows_tested <- if (name %in% c("S0", "O0")) intersect(F_filt, rownames(f$matrix)) else rownames(f$matrix)
            rows_tested <- rows_tested[rowSums(is.finite(f$matrix[rows_tested, , drop = FALSE])) > 0]
            tab <- da(f$matrix[rows_tested, sim$runs, drop = FALSE])
            fits[[name]]$da <- tab
            fits[[name]]$seconds <- sum(fits[[name]]$seconds, attr(tab, "seconds"), na.rm = TRUE)
        }
        if (length(F_set)) {
            rows <- c(rows, scoreArm(sim, arm, fits, F_set, F_filt, k, common, Z, G))
        } else {
            for (name in names(fits)) rows[[length(rows) + 1L]] <- do.call(newRow, c(common, list(corr = name, stratum = "all",
                status = fits[[name]]$status, class = fits[[name]]$class, k = k, n_set = 0L)))
        }
        timing[["estimands"]] <- timing[["estimands"]] + proc.time()[[3L]] - t_est
    }
    timing[["DA"]] <- t_da
    timing[["estimands"]] <- timing[["estimands"]] - t_da
    rows <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))[, COLUMNS]
    # checks of the lifecycle that read the none arm (G1) and the refusal (E7), on the first replicate
    if (lifecycle) lifecycleChecks(sim, arms, rows, k, sets, variant, da_errors)
    list(rows = rows, timing = data.frame(cell = cellname, part = part, variant = variant, rep = r, seconds = proc.time()[[3L]] - t_unit,
        t(timing), warnings = warning_count$n, da_errors = length(da_errors),
        da_error_first = if (length(da_errors)) da_errors[[1L]] else NA_character_, stringsAsFactors = FALSE))
}

# ---- one replicate of part B -----------------------------------------------------------------------------------
LEGACY_PROT <- list(norm_method = "none", ruv_mode = "automatic", auto_percentage_min = 1, auto_percentage_max = 20,
    separation_metric = "max_difference", k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, ruv_percentage = NA, ruv_k = NA,
    grouping = "group")
DECLARED_PROT <- c(LEGACY_PROT, list(ruv_design = "groups_all_features", ruv_k_rule = "num_sv_be"))
LEGACY_OMIC <- list(ruv_mode = "automatic", ruv_grouping_variable = "group", percentage_min = 1, percentage_max = 20,
    separation_metric = "max_difference", k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, manual_k = 2, manual_percentage = 5)
DECLARED_OMIC <- list(ruv_design = "groups_all_features", ruv_grouping_variable = "group", ruv_k = NA, ruv_k_rule = "num_sv_be")
omicBinding <- function(omic) if (identical(omic, "lipid")) lipidomicsDeclarationBinding else metabolomicsDeclarationBinding
# the controls a route chose, resolved as the route's RUV method resolves them: proteomics by name
# (`names(ctrl[which(ctrl)])`, R/func_prot_s4_norm_methods.R); lipidomics and metabolomics by position when the
# logical vector has one entry per feature of the assay ("Assume correct order", R/func_lipid_s4_ruv_methods.R), else
# by name
controlIds <- function(index, ids, by_name = FALSE) {
    if (is.null(index)) return(character())
    if (is.character(index)) return(intersect(index, ids))
    names_of <- names(index)
    index <- as.logical(index)
    index[is.na(index)] <- FALSE
    if (!by_name && length(index) == length(ids)) return(ids[index])
    if (!is.null(names_of)) return(intersect(names_of[index], ids))
    ids[index]
}
# run one route of part B on one arm: the shipped state, its matrix, k and controls, and the comparator beside it
routeB <- function(sim, arm, route, sets) {
    omic <- sim$p$omic
    store <- memStore()
    if (identical(omic, "prot")) {
        run <- attempt(quiet(diaNormalizationOperation(route)$kernel(arm$filtered, if (identical(route, "legacy")) LEGACY_PROT else DECLARED_PROT, store)),
            partb_refusal_classes)
        if (!identical(run$status, "ok")) return(list(status = run$status, class = run$class, seconds = run$seconds))
        value <- run$value
        state <- get(value$final, envir = store$saved)
        k <- as.integer(value$best_k)
        if (is.na(k) || k == 0L) return(list(status = "none", class = "k0", seconds = run$seconds, state = state, k = 0L,
            matrix = protMatrix(state)[, sim$runs, drop = FALSE]))
        controls <- controlIds(value$controls, arm$ids, by_name = TRUE)
        percentage <- value$percentage %||% NA
    } else {
        binding <- omicBinding(omic)
        run <- attempt(quiet(omicsRuvOperation(binding, route)$kernel(arm$filtered, if (identical(route, "legacy")) LEGACY_OMIC else DECLARED_OMIC, store)),
            partb_refusal_classes)
        if (!identical(run$status, "ok")) return(list(status = run$status, class = run$class, seconds = run$seconds))
        value <- run$value
        state <- get(value$final, envir = store$saved)
        if (identical(route, "legacy")) {
            optimisation <- get(value$optimisation, envir = store$saved)
            k_list <- binding$extract_best_k(optimisation)
            control_list <- binding$extract_controls(optimisation)
            after <- methods::slot(state, omicSlot(omic))
            before <- methods::slot(arm$filtered, omicSlot(omic))
            if (!identical(names(after), names(before)) || any(vapply(names(before), function(a) identical(after[[a]], before[[a]]), logical(1L)))) {
                return(list(status = "refused", class = "assay_not_corrected", seconds = run$seconds))
            }
        } else {
            k_list <- lapply(value$factors, `[[`, "k")
            control_list <- lapply(names(value$factors), function(a) if (identical(value$factors[[a]]$status, "ok")) sim$assay[sim$assay == a] else NULL)
            names(control_list) <- names(value$factors)
            control_list <- lapply(control_list, function(x) if (is.null(x)) NULL else names(x)[names(x) %in% arm$ids])
        }
        k <- vapply(k_list, function(x) if (is.null(x) || is.na(x)) 0L else as.integer(x), integer(1L))
        controls <- unlist(lapply(names(control_list), function(a) controlIds(control_list[[a]], arm$ids[sim$assay[arm$ids] == a])), use.names = FALSE)
        percentage <- NA
    }
    matrix <- if (identical(omic, "prot")) protMatrix(state)[, sim$runs, drop = FALSE] else omicMatrix(state, omic, sim$runs)
    # the comparator with the k and controls the route chose (per assay for lipidomics and metabolomics)
    comparator <- attempt({
        if (identical(omic, "prot")) {
            fitReference(sim, arm, k, sets, "complete", controls = controls)
        } else {
            pieces <- lapply(names(k), function(a) {
                ids <- arm$ids[sim$assay[arm$ids] == a]
                if (k[[a]] == 0L) return(arm$Y[ids, , drop = FALSE])
                sub_arm <- arm
                sub_arm$ids <- ids
                out <- fitReference(sim, sub_arm, k[[a]], sets, "complete", controls = intersect(controls, ids), raw = TRUE)
                out <- out[, sim$runs, drop = FALSE]
                if (identical(route, "declared")) out <- restoreRuvFeatureLevels(out, arm$Y[rownames(out), sim$runs, drop = FALSE])
                out[rowSums(is.finite(out)) > 0, , drop = FALSE]
            })
            do.call(rbind, pieces)
        }
    })
    list(status = "ok", class = NA_character_, seconds = run$seconds, state = state, k = k, matrix = matrix, controls = controls,
        percentage = percentage, comparator = comparator, ruv_factors = if (!identical(omic, "prot") && identical(route, "declared")) state@args$ruv_factors)
}
subsetState <- function(sim, state, ids) {
    if (identical(sim$p$omic, "prot")) {
        keep <- state@protein_quant_table[[state@protein_id_column]] %in% ids
        state@protein_quant_table <- state@protein_quant_table[keep, , drop = FALSE]
        rownames(state@protein_quant_table) <- NULL # VAL-98 (1): the DA's column_to_rownames needs default row names
        if (!is.null(state@args$ruv_uncorrected_quant)) {
            state@args$ruv_uncorrected_quant <- state@args$ruv_uncorrected_quant[rownames(state@args$ruv_uncorrected_quant) %in% ids, , drop = FALSE]
        }
        return(state)
    }
    tables <- methods::slot(state, omicSlot(sim$p$omic))
    methods::slot(state, omicSlot(sim$p$omic)) <- lapply(tables, function(t) t[t[[ID_COL]] %in% ids, , drop = FALSE])
    state
}
routeDa <- function(sim, arm, result, route, ids = NULL) {
    omic <- sim$p$omic
    if (identical(omic, "prot")) {
        state <- if (is.null(ids)) result$state else subsetState(sim, result$state, ids)
        return(daTable(sim, object = state))
    }
    if (identical(route, "legacy")) {
        state <- if (is.null(ids)) result$state else subsetState(sim, result$state, ids)
        return(daTable(sim, object = state))
    }
    # the declared design: the uncorrected values with W (utils_workflow_omics_declarations.R:801)
    tested <- arm$filtered
    if (!is.null(ids)) tested <- subsetState(sim, tested, ids)
    daTable(sim, object = tested, ruv_factors = result$ruv_factors)
}

unitB <- function(cellname, r) {
    t_unit <- proc.time()[[3L]]
    warning_count$n <- 0L
    p <- CELLS[[cellname]]
    use_presence <- cellname %in% presence_cells
    sim <- generateReplicate(cellname, r)
    sets <- replicateSets(sim)
    G <- groupDesign(sim$group)
    Z <- sim$Z
    masks <- masksOf("B", cellname)
    routes <- if (cellname %in% c("P0", "L2", "M2")) c("legacy", "declared") else "legacy"
    arms <- list()
    for (mask in masks) {
        arm <- buildArm(sim, mask, "none", use_presence, sets)
        arm$O <- arm$Y - sim$unwanted[arm$ids, sim$runs, drop = FALSE]
        arm$routes <- lapply(stats::setNames(routes, routes), function(route) routeB(sim, arm, route, sets))
        arms[[mask]] <- arm
    }
    none <- arms[["none"]]
    rows <- list()
    for (mask in masks) {
        arm <- arms[[mask]]
        common <- list(cell = cellname, part = "B", variant = "base", rep = r, mask = mask, mech = maskMech(mask), rate = maskRate(mask),
            n_filt = length(arm$ids))
        rows[[length(rows) + 1L]] <- do.call(newRow, c(common, list(corr = "mask", stratum = "all", status = "ok",
            removed = arm$removed, miss_share = arm$miss_share, miss_run_min = min(arm$run_share), miss_run_med = stats::median(arm$run_share),
            miss_run_max = max(arm$run_share), block_share = arm$block_share, seconds = arm$seconds)))
        F_filt <- arm$ids
        for (route in routes) {
            tag <- if (identical(route, "legacy")) "L" else "D"
            res <- arm$routes[[route]]
            res_none <- none$routes[[route]]
            k_route <- if (is.null(res$k)) NA_integer_ else max(res$k)
            obs <- arm$obs
            all_obs <- matrix(TRUE, length(F_filt), length(sim$runs), dimnames = list(F_filt, sim$runs))
            entry <- function(kind, matrix, status, class, seconds, base, ref_O, ref_N, obs_matrix, extra = list()) {
                c(list(name_kind = kind, matrix = matrix, status = status, class = class, seconds = seconds, base = base, ref_O = ref_O,
                    ref_N = ref_N, obs = obs_matrix), extra)
            }
            extra <- list(k = k_route, k_chosen = paste(res$k, collapse = ";"), n_controls = length(res$controls),
                percentage = res$percentage %||% NA)
            ok_matrix <- identical(res$status, "ok")
            fits <- list()
            fits[[paste0("SB", tag)]] <- entry("S", if (ok_matrix) res$matrix, if (res$status %in% c("ok", "none")) res$status else res$status,
                res$class, res$seconds, arm$Y, "O", "N", obs, extra)
            if (identical(res$status, "none")) fits[[paste0("SB", tag)]]$status <- "none"
            comparator <- res$comparator
            fits[[paste0("CB", tag)]] <- entry("C", if (ok_matrix && identical(comparator$status, "ok")) comparator$value,
                if (ok_matrix) comparator$status else "not_run", if (ok_matrix) comparator$class else NA_character_,
                if (ok_matrix) comparator$seconds else NA, arm$Y, "O", "N", obs, extra)
            fits$O <- entry("O", arm$O, "ok", NA_character_, NA, arm$Y, "O", "N", obs)
            fits$N <- entry("N", arm$Y, "ok", NA_character_, NA, arm$Y, "O", "N", obs)
            if (!identical(mask, "none")) {
                fits[[paste0("SB", tag, "0")]] <- entry("S", if (identical(res_none$status, "ok")) res_none$matrix, res_none$status,
                    res_none$class, NA, none$Y, "O0", "N0", all_obs)
                fits$O0 <- entry("O", none$O, "ok", NA_character_, NA, none$Y, "O0", "N0", all_obs)
                fits$N0 <- entry("N", none$Y, "ok", NA_character_, NA, none$Y, "O0", "N0", all_obs)
            }
            # DA: the route's own test on its state; the comparator by the legacy test (the declared DA does not read
            # the corrected matrix, E6); the none state on the rows of F_filt
            S_name <- paste0("SB", tag)
            if (res$status %in% c("ok", "none")) fits[[S_name]]$da <- tryCatch(routeDa(sim, arm, res, route), error = function(e) NULL)
            C_name <- paste0("CB", tag)
            if (identical(route, "legacy") && identical(fits[[C_name]]$status, "ok")) {
                fits[[C_name]]$da <- tryCatch(daTable(sim, fits[[C_name]]$matrix), error = function(e) NULL)
            }
            S0_name <- paste0("SB", tag, "0")
            if (!identical(mask, "none") && res_none$status %in% c("ok", "none")) {
                # VAL-98 (1): the error is recorded (in the row's class), not dropped
                fits[[S0_name]]$da <- tryCatch(routeDa(sim, none, res_none, route, ids = intersect(F_filt, none$ids)), error = function(e) {
                    fits[[S0_name]]$class <<- substr(paste("da_error:", gsub("[\r\n\t]+", " ", conditionMessage(e))), 1L, 300L)
                    NULL
                })
            }
            finiteRows <- function(m) if (is.null(m)) character() else rownames(m)[rowSums(is.finite(m)) > 0]
            F_set <- if (ok_matrix && identical(res_none$status, "ok")) Reduce(intersect, list(F_filt, finiteRows(res$matrix),
                if (identical(comparator$status, "ok")) finiteRows(comparator$value) else F_filt, finiteRows(res_none$matrix))) else character()
            route_common <- utils::modifyList(common, list(variant = route))
            if (length(F_set)) {
                rows <- c(rows, scoreArm(sim, arm, fits, F_set, F_filt, k_route, route_common, Z, G))
            } else {
                for (name in names(fits)) {
                    f <- fits[[name]]
                    row <- do.call(newRow, c(route_common, list(corr = name, stratum = "filt", status = f$status, class = f$class,
                        k = k_route, k_chosen = extra$k_chosen, n_controls = extra$n_controls, seconds = f$seconds, n_set = length(F_filt))))
                    if (!is.null(f$da)) {
                        metrics <- daMetrics(f$da, F_filt, sim$beta)
                        row[names(metrics)] <- as.list(metrics)
                    }
                    rows[[length(rows) + 1L]] <- row
                }
            }
        }
    }
    rows <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))[, COLUMNS]
    list(rows = rows, timing = data.frame(cell = cellname, part = "B", variant = "base", rep = r, seconds = proc.time()[[3L]] - t_unit,
        filter = NA, S = NA, C = NA, Z = NA, DA = NA, estimands = NA, warnings = warning_count$n, stringsAsFactors = FALSE))
}

# ---- lifecycle checks (KIT-185 7.4, gates G1 and G2) --------------------------------------------------------------
lifecycleChecks <- function(sim, arms, rows, k, sets, variant, da_errors) {
    cellname <- sim$cell
    r <- sim$r
    # G1: on none, S recovers the planted directions (rho_min >= 0.9) and leaves at most 0.10 (B on F^comp)
    g1 <- rows[rows$mask == "none" & rows$corr == "S" & rows$stratum == "comp", , drop = FALSE]
    # G1 is the gate of the part-A cells at the planted k (not of the sensitivity variants)
    if (nrow(g1) && !identical(part, "S")) addCheck("G1", cellname, r, "none", g1$rho_min[[1L]], 0.9, isTRUE(g1$rho_min[[1L]] >= 0.9 && g1$B[[1L]] <= 0.10),
        sprintf("B(S, none) on F^comp %.4f; variant %s", g1$B[[1L]], variant))
    # G2: the shipped DA of lipidomics and metabolomics on the masked arms (failed assays, any correction)
    failed <- sum(as.numeric(rows$da_failed[rows$stratum == "filt"]), na.rm = TRUE)
    addCheck(if (identical(sim$p$omic, "prot")) "DA_errors" else "G2", cellname, r, "all", failed + length(da_errors), 0,
        failed + length(da_errors) == 0, sprintf("variant %s; failed assays %d; DA calls that stopped %d%s", variant, failed,
            length(da_errors), if (length(da_errors)) paste0("; first: ", da_errors[[1L]]) else ""))
    # E1-neg: the complete-feature comparator differs from S on at least one masked arm
    differs <- vapply(names(arms)[-1L], function(mask) {
        f <- arms[[mask]]$fit
        if (!identical(f$S$status, "ok") || !identical(f$C$status, "ok")) return(NA_real_)
        maxDiff(f$S$value, f$C$value)$max
    }, numeric(1L))
    addCheck("E1neg", cellname, r, "masked", max(differs, na.rm = TRUE), 1e-6, isTRUE(max(differs, na.rm = TRUE) > 1e-6), variant)
    if (r != first_rep || !identical(part, "A") || !identical(variant, "base")) return(invisible(NULL))
    # E1 (b): upstream's R reference (installed) against the compiled call, one replicate per cell
    for (mask in intersect(c("none", "vd10", "mcar2", "bb10"), names(arms))) {
        arm <- arms[[mask]]
        ids <- if (identical(sim$p$omic, "prot")) arm$ids else arm$ids[sim$assay[arm$ids] == sim$assay[arm$ids[[1L]]]]
        fixed <- ids[sampleFixed(length(ids), min(200L, length(ids)))]
        incomplete <- setdiff(ids[rowSums(!arm$obs[ids, , drop = FALSE]) > 0], fixed)
        chosen <- c(fixed, utils::head(incomplete, 400L - length(fixed)))
        Y <- t(arm$Y[ids, sim$runs, drop = FALSE])
        result <- tryCatch({
            compiled <- RUVIIIC::RUVIII_C_Varying(k = k, Y = Y, M = sets$M, toCorrect = chosen, potentialControls = ids)
            reference <- quiet(RUVIIIC::RUVIII_C_Varying(k = k, Y = Y, M = sets$M, toCorrect = chosen, potentialControls = ids,
                version = "R", filename = NULL, progress = FALSE))
            single <- vapply(chosen, function(id) {
                runs_seen <- which(!is.na(Y[, id]))
                sum(colSums(sets$M[runs_seen, , drop = FALSE]) > 0) < 2L
            }, logical(1L))
            keep <- chosen[!single]
            a <- t(compiled[, keep, drop = FALSE])
            b <- t(reference[, keep, drop = FALSE])
            d <- maxDiff(a, b)
            list(d = d, single = sum(single))
        }, error = function(e) list(error = conditionMessage(e)))
        if (!is.null(result$error)) {
            addCheck("E1b", cellname, r, mask, NA, NA, FALSE, substr(result$error, 1L, 200L))
        } else {
            tolerance <- if (length(sim$runs) > 20L) 1e-6 else 1e-8
            addCheck("E1b", cellname, r, mask, result$d$max, tolerance, result$d$same_finite && result$d$max <= tolerance,
                sprintf("features in one replicate set excluded %d; share>1e-8 %.4g; targets %d", result$single, result$d$share_1e8, length(chosen)))
        }
    }
    # E7: a refusal by the route's own rule (k one above the controls observed in every sample)
    arm <- arms[["none"]]
    complete_ids <- arm$ids[rowSums(!arm$obs) == 0L]
    few <- unlist(lapply(split(complete_ids, sim$assay[complete_ids]), utils::head, 3L), use.names = FALSE)
    if (identical(sim$p$omic, "prot")) few <- utils::head(complete_ids, 3L)
    ctrl <- if (identical(sim$p$omic, "prot")) stats::setNames(arm$ids %in% few, arm$ids) else few
    refused <- attempt(fitShipped(sim, arm, 4L, ctrl = ctrl))
    addCheck("E7_refusal", cellname, r, "none", NA, NA, identical(refused$status, "refused"), sprintf("status %s; class %s", refused$status, refused$class))
    # E3: part A's S equals the declared kernel's corrected table (P0, L2, M2 at VD 10 %)
    if (cellname %in% c("P0", "L2", "M2") && "vd10" %in% names(arms) && identical(arms$vd10$fit$S$status, "ok")) {
        arm <- arms$vd10
        e3 <- tryCatch({
            store <- memStore()
            if (identical(sim$p$omic, "prot")) {
                filtered <- filterArm(sim, { m <- sim$y; m[maskMatrix(sim, "vd10")] <- NA_real_; m }, FALSE)$object
                parameters <- utils::modifyList(LEGACY_PROT, list(ruv_mode = "manual", ruv_percentage = 5, ruv_k = 2,
                    ruv_design = "groups_all_features", ruv_k_rule = "user"))
                value <- quiet(diaNormalizationOperation("declared")$kernel(filtered, parameters, store))
                state <- get(value$final, envir = store$saved)
                kernel <- protMatrix(state)
                S <- arm$fit$S$value
                same <- identical(unname(S), unname(kernel)) && identical(dimnames(S), dimnames(kernel))
                # E6: the declared DA does not read the corrected matrix
                swapped <- protWithTable(state, arm$fit$C$value)
                first <- daTable(sim, object = state)
                second <- daTable(sim, object = swapped)
                e6 <- identical(first$p, second$p[match(first$id, second$id)]) && identical(sort(first$id), sort(second$id))
                list(value = if (same) 0 else NA, pass = same, e6 = e6)
            } else {
                parameters <- list(ruv_design = "groups_all_features", ruv_grouping_variable = "group", ruv_k = 2, ruv_k_rule = "user")
                value <- quiet(omicsRuvOperation(omicBinding(sim$p$omic), "declared")$kernel(arm$object, parameters, store))
                state <- get(value$final, envir = store$saved)
                kernel <- omicMatrix(state, sim$p$omic, sim$runs)
                S <- arm$fit$S$value
                common <- intersect(rownames(S), rownames(kernel))
                centre <- function(m) m - rowMeans(m, na.rm = TRUE)
                d <- max(abs(centre(S[common, , drop = FALSE]) - centre(kernel[common, , drop = FALSE])), na.rm = TRUE)
                same <- setequal(rownames(S), rownames(kernel)) && d <= 1e-10
                e6 <- NA
                if (identical(sim$cell, "L2")) {
                    # the declared DA's input rule (utils_workflow_omics_declarations.R:801): per assay the uncorrected
                    # values the state carries, else the state's own assay
                    declaredInput <- function(state) {
                        assays <- methods::slot(state, omicSlot(sim$p$omic))
                        tables <- lapply(names(assays), function(a) state@args$ruv_uncorrected_assays[[a]] %||% assays[[a]])
                        names(tables) <- names(assays)
                        tested <- arm$object
                        methods::slot(tested, omicSlot(sim$p$omic)) <- tables
                        tested
                    }
                    first <- daTable(sim, object = declaredInput(state), ruv_factors = state@args$ruv_factors)
                    swapped <- omicWithMatrix(state, sim$p$omic, arm$fit$C$value, sim$assay)
                    second <- daTable(sim, object = declaredInput(swapped), ruv_factors = swapped@args$ruv_factors)
                    e6 <- identical(first$p, second$p[match(first$id, second$id)])
                }
                list(value = d, pass = same, e6 = e6)
            }
        }, error = function(e) list(value = NA, pass = FALSE, e6 = NA, error = conditionMessage(e)))
        addCheck("E3", cellname, r, "vd10", e3$value, if (identical(sim$p$omic, "prot")) "identical" else "1e-10", isTRUE(e3$pass),
            e3$error %||% "")
        if (!is.na(e3$e6 %||% NA)) addCheck("E6", cellname, r, "vd10", NA, "identical p", isTRUE(e3$e6),
            if (identical(sim$p$omic, "prot")) "differentialAbundanceAnalysis on the declared state" else
                "runLipidsDA on the state's uncorrected assays with its factors (the rule of utils_workflow_omics_declarations.R:801)")
    }
    invisible(NULL)
}
# a fixed random choice that leaves the random number state as it was
sampleFixed <- function(n, size) {
    old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()))
    set.seed(185L)
    sort(sample.int(n, size))
}
# E4: the estimand functions on noise-free constructed cases
checkE4 <- function() {
    set.seed(1L)
    runs <- 12L
    group <- rep(c("control", "treated"), each = 6L)
    G <- groupDesign(group)
    Z <- cbind(scale1(stats::lm.fit(G, stats::rnorm(runs))$residuals))
    Z <- cbind(Z, scale1(stats::lm.fit(cbind(G, Z), stats::rnorm(runs))$residuals))
    n <- 50L
    loadings <- matrix(stats::rnorm(n * 2L), n, 2L)
    U <- loadings %*% t(Z)
    beta <- rep(c(1, -1, 0, 0, 0), 10L)
    y <- 20 + outer(beta, as.numeric(group == "treated")) + U
    obs <- matrix(TRUE, n, runs)
    terms <- function(m) projectionTerms(m, obs, G, Z)$U
    B <- (sum(terms(y - 0.5 * U)) - sum(terms(y - U))) / (sum(terms(y)) - sum(terms(y - U)))
    addCheck("E4_B", "constructed", 0L, "none", B, "0.25 +- 1e-12", abs(B - 0.25) <= 1e-12)
    rho <- recoveryComplete(y, y - U, Z)$rho
    addCheck("E4_rho", "constructed", 0L, "none", min(rho), "1 +- 1e-12", all(abs(rho - 1) <= 1e-12))
    truth <- 20 + outer(beta, as.numeric(group == "treated"))
    shrunk <- 20 + outer(0.9 * beta, as.numeric(group == "treated"))
    cstar <- slopeOf(truth, obs, group, beta) - slopeOf(shrunk, obs, group, beta)
    addCheck("E4_Cstar", "constructed", 0L, "none", cstar, "0.10 +- 1e-12", abs(cstar - 0.10) <= 1e-12)
}

# ---- the job ---------------------------------------------------------------------------------------------------
variants_of <- function(part) {
    if (identical(part, "S")) {
        wanted <- Sys.getenv("K185_VARIANTS", "k1,k4,s1,s10,loess")
        return(strsplit(wanted, ",", fixed = TRUE)[[1L]])
    }
    "base"
}
units <- do.call(rbind, lapply(cells_wanted, function(cellname) {
    variants <- variants_of(part)
    if (lifecycle && identical(part, "A") && !identical(CELLS[[cellname]]$omic, "prot")) variants <- c("base", "presence")
    expand.grid(rep = seq(first_rep, last_rep), variant = variants, cell = cellname, stringsAsFactors = FALSE)
}))
job <- Sys.getenv("K185_JOB", sprintf("%s-%s-%d-%d%s", tolower(part), paste(cells_wanted, collapse = "_"), first_rep, last_rep,
    if (lifecycle) "-l1" else ""))
rows_dir <- file.path(out_dir, "rows")
dir.create(rows_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("KIT-185 %s: %d units, %d workers, lifecycle %s\n", job, nrow(units), workers, lifecycle))
if (lifecycle) checkE4()
runUnit <- function(i) {
    u <- units[i, ]
    checks_env$rows <- list()
    result <- tryCatch(if (identical(part, "B")) unitB(u$cell, u$rep) else unitA(u$cell, u$rep, part, u$variant),
        error = function(e) list(error = paste(class(e)[[1L]], conditionMessage(e))))
    if (!is.null(result$error)) {
        rows <- as.data.frame(do.call(newRow, list(cell = u$cell, part = part, variant = u$variant, rep = u$rep, mask = "all", corr = "unit",
            stratum = "all", status = "crashed", class = substr(gsub("[\r\n\t]+", " ", result$error), 1L, 300L))), stringsAsFactors = FALSE)[, COLUMNS]
        result <- list(rows = rows, timing = NULL)
    }
    rows <- result$rows
    # written to 7 significant digits (KIT-185 declaration: output precision only)
    numeric_columns <- vapply(rows, is.numeric, logical(1L))
    rows[numeric_columns] <- lapply(rows[numeric_columns], signif, 7L)
    hash <- NA_character_
    if (lifecycle) {
        temp <- tempfile(fileext = ".tsv")
        utils::write.table(rows, temp, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
        hash <- strsplit(system2("sha256sum", temp, stdout = TRUE), " ")[[1L]][[1L]]
        unlink(temp)
        blind <- rows$mask != "none" & rows$corr != "mask"
        rows[blind, setdiff(COLUMNS, BLIND_KEEP)] <- NA
    }
    part_file <- file.path(rows_dir, sprintf("part-%s-%s-%d.tsv", u$cell, u$variant, u$rep))
    utils::write.table(rows, part_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    checks <- if (length(checks_env$rows)) do.call(rbind, checks_env$rows)
    list(timing = if (!is.null(result$timing)) cbind(result$timing, rows_sha256 = hash), checks = checks)
}
results <- parallel::mclapply(seq_len(nrow(units)), runUnit, mc.cores = workers, mc.preschedule = FALSE)
for (i in seq_along(results)) if (!is.list(results[[i]])) {
    u <- units[i, ]
    rows <- as.data.frame(do.call(newRow, list(cell = u$cell, part = part, variant = u$variant, rep = u$rep, mask = "all", corr = "unit",
        stratum = "all", status = "crashed", class = substr(paste("worker failed:", paste(as.character(results[[i]]), collapse = " ")), 1L, 300L))),
        stringsAsFactors = FALSE)[, COLUMNS]
    utils::write.table(rows, file.path(rows_dir, sprintf("part-%s-%s-%d.tsv", u$cell, u$variant, u$rep)), sep = "\t", quote = FALSE,
        row.names = FALSE, na = "NA")
    results[[i]] <- list(timing = NULL, checks = NULL)
}
# assemble: one table per cell (and part), the parts removed after
for (cellname in cells_wanted) {
    files <- list.files(rows_dir, pattern = sprintf("^part-%s-.*\\.tsv$", cellname), full.names = TRUE)
    files <- files[order(as.integer(sub(".*-([0-9]+)\\.tsv$", "\\1", files)), files)]
    table <- do.call(rbind, lapply(files, utils::read.delim, stringsAsFactors = FALSE, na.strings = "NA", colClasses = "character"))
    out_file <- file.path(rows_dir, sprintf("%s-%s-%d-%d%s.tsv.gz", tolower(part), cellname, first_rep, last_rep, if (lifecycle) "-l1" else ""))
    connection <- gzfile(out_file, "w")
    utils::write.table(table, connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    close(connection)
    unlink(files)
}
timing <- do.call(rbind, lapply(results, `[[`, "timing"))
utils::write.table(timing, file.path(out_dir, sprintf("timing-%s.tsv", job)), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
checks <- do.call(rbind, lapply(results, `[[`, "checks"))
if (lifecycle) {
    checks <- rbind(do.call(rbind, checks_env$rows), checks)
    # E7 (crash): the deliberate fault was recorded crashed and the job went on
    if (fault) {
        first_file <- file.path(rows_dir, sprintf("%s-%s-%d-%d-l1.tsv.gz", tolower(part), cells_wanted[[1L]], first_rep, last_rep))
        written <- utils::read.delim(gzfile(first_file), stringsAsFactors = FALSE)
        crashed <- written[written$status == "crashed" & grepl("deliberate fault", written$class, fixed = TRUE), , drop = FALSE]
        checks <- rbind(checks, data.frame(check = "E7_crash", cell = cells_wanted[[1L]], rep = first_rep, mask = "masked", value = nrow(crashed),
            threshold = ">= 1", pass = nrow(crashed) >= 1L && nrow(timing) == nrow(units), note = sprintf("units completed %d of %d",
                nrow(timing), nrow(units)), stringsAsFactors = FALSE))
    }
}
if (!is.null(checks) && nrow(checks)) utils::write.table(checks, file.path(out_dir, sprintf("checks-%s.tsv", job)), sep = "\t", quote = FALSE,
    row.names = FALSE, na = "NA")
# provenance
sha <- function(f) strsplit(system2("sha256sum", f, stdout = TRUE), " ")[[1L]][[1L]]
tables <- list.files(out_dir, recursive = TRUE, full.names = TRUE, pattern = "\\.(tsv|gz)$")
library_table <- utils::installed.packages()[, c("Package", "Version")]
key_packages <- c("MultiScholaR", "RUVIIIC", "ruv", "RSpectra", "limma", "qvalue", "progress", "sva", "dplyr", "tidyr")
provenance <- list(job = job, part = part, cells = cells_wanted, replicates = c(first_rep, last_rep), lifecycle = lifecycle,
    fault = fault, presence_cells = presence_cells, variants = variants_of(part), workers = workers,
    commit = tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE), error = function(e) NA),
    r_version = R.version.string, blas = extSoftVersion()[["BLAS"]], lapack = La_library(),
    omp_threads = Sys.getenv("OMP_NUM_THREADS"), runner_image = paste(Sys.getenv("ImageOS"), Sys.getenv("ImageVersion")),
    processor = tryCatch(sub("^model name\\s*:\\s*", "", grep("^model name", readLines("/proc/cpuinfo"), value = TRUE)[[1L]]), error = function(e) NA),
    seed_rule = "set.seed(185000000 + 10000 * cell index + replicate, kind = 'Mersenne-Twister', normal.kind = 'Inversion', sample.kind = 'Rejection')",
    cell_index = vapply(CELLS[cells_wanted], `[[`, integer(1L), "idx"),
    key_packages = as.list(stats::setNames(library_table[match(key_packages, library_table[, "Package"]), "Version"], key_packages)),
    script_sha256 = c(sim = sha(script_file), ref = sha(file.path(here, "ruviiic-ref.R"))),
    started = format(started_at, "%Y-%m-%dT%H:%M:%S%z"), finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    units = nrow(units), units_crashed = sum(vapply(list.files(rows_dir, full.names = TRUE), function(f) {
        x <- utils::read.delim(gzfile(f), stringsAsFactors = FALSE)
        sum(x$corr == "unit" & x$status == "crashed")
    }, numeric(1L))), exit_status = 0L,
    tables_sha256 = stats::setNames(as.list(vapply(tables, sha, character(1L))), sub(paste0("^", out_dir, "/?"), "", tables)))
jsonlite::write_json(provenance, file.path(out_dir, sprintf("provenance-%s.json", job)), auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat(sprintf("KIT-185 %s finished: %d units, wall %.1f s\n", job, nrow(units), as.numeric(difftime(Sys.time(), started_at, units = "secs"))))
