# KIT-114 (b) / KIT-123 R3: how the DA should account for the classic route's RUV step. The simulator of
# d_ruv_controls.R (1,000 proteins, 20% true effects of 1 log2, one unwanted factor per sample independent of group,
# the group as the replicate structure), through the package's own sweep (findBestNegCtrlPercentage), correction
# (ruvIII_C_Varying) and DA kernel (runTestsContrasts). Arms, cells and the rule are fixed in tickets/kit/KIT-123.md
# before any run.
#   no_ruv      DA on the uncorrected values
#   today       the sweep's controls and k, RUV-III-C correction, DA on the corrected values
#   cov_sel     the sweep's controls and k, W (from the proteins complete in every sample) as covariates in the DA
#   corr_true   as today, with the same number of true null proteins as controls (diagnostic)
#   cov_true    as cov_sel, with those true null controls (diagnostic)
#   corr_wide   as today, controls = every protein outside the 25% most significant of a first pass
#   cov_wide    as cov_sel, with those wide controls
#   *_k         today, cov_sel, cov_wide at k limited to calculateAdaptiveMaxK(samples) (KIT-123 R1)
# Usage (from the repository root): Rscript <this> <out.tsv> [seeds = 20] [cores = 4]
args <- commandArgs(TRUE)
seeds <- if (length(args) > 1L) as.integer(args[[2L]]) else 20L
cores <- if (length(args) > 2L) as.integer(args[[3L]]) else 4L
suppressMessages(pkgload::load_all(".", quiet = TRUE))
logger::log_threshold(logger::FATAL)
try(RUVIIIC::omp_set_num_threads(1L), silent = TRUE)
quietly <- function(expr) suppressWarnings(suppressMessages(utils::capture.output(value <- expr, type = "output")))
simulate <- function(seed, n, missing, m = 1000L, pi1 = 0.20, effect = 1.0, d0 = 4, s0 = 0.25) {
    set.seed(seed)
    mu <- rnorm(m, 22, 2)
    sd <- sqrt(d0 * s0^2 / rchisq(m, d0))
    de <- seq_len(m) <= round(pi1 * m)
    delta <- ifelse(de, effect * sample(c(-1, 1), m, replace = TRUE), 0)
    groups <- rep(c("control", "treated"), each = n)
    w <- rnorm(2 * n)                       # one unwanted factor per sample
    alpha <- rnorm(m, 0, 0.5)               # its loading per protein
    x <- sapply(seq_len(2 * n), function(j) mu + (groups[[j]] == "treated") * delta + alpha * w[[j]] + rnorm(m, 0, sd))
    if (missing > 0) {                      # intensity-dependent missingness at the target overall share
        p_at <- function(lod) mean(stats::plogis((lod - x) / 1.0))
        lod <- stats::uniroot(function(l) p_at(l) - missing, c(0, 60))$root
        x[stats::runif(length(x)) < stats::plogis((lod - x) / 1.0)] <- NA_real_
    }
    ids <- sprintf("PG_%05d", seq_len(m)); runs <- sprintf("run_%03d", seq_len(2 * n))
    dimnames(x) <- list(ids, runs)
    keep <- rowSums(is.finite(x[, groups == "control"])) >= n / 2 & rowSums(is.finite(x[, groups == "treated"])) >= n / 2
    x <- x[keep, , drop = FALSE]
    design <- data.frame(Run = runs, group = factor(groups, levels = c("control", "treated")), replicates = seq_len(2 * n),
        tech_rep_group = paste0("s", seq_len(2 * n)), stringsAsFactors = FALSE)
    object_design <- design
    rownames(design) <- runs
    quant <- as.data.frame(c(list(Protein.Group = rownames(x)), as.data.frame(x)), check.names = FALSE, stringsAsFactors = FALSE)
    object <- methods::new("ProteinQuantitativeData", protein_quant_table = quant, protein_id_column = "Protein.Group",
        protein_id_table = data.frame(Protein.Group = rownames(x), stringsAsFactors = FALSE), design_matrix = object_design,
        sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
        args = list(globalParameters = list(workflow_type = "LFQ")))
    list(object = object, x = x, de = stats::setNames(de, ids)[rownames(x)], design = design, unwanted = w)
}
# The DA kernel; `w` (samples by k) joins the design and the formula as W1..Wk.
da <- function(matrix, design, truth, w = NULL) {
    formula <- "~ 0 + group"
    if (!is.null(w)) {
        design <- cbind(design, as.data.frame(w[rownames(design), , drop = FALSE]))
        formula <- paste(formula, "+", paste(colnames(w), collapse = " + "))
    }
    value <- NULL
    quietly(value <- runTestsContrasts(matrix, "grouptreated-groupcontrol", design, formula))
    tab <- value$results[[1L]]
    called <- !is.na(tab$fdr_qvalue) & tab$fdr_qvalue < 0.05
    t <- truth[rownames(tab)]
    c(fdp = if (sum(called)) sum(called & !t) / sum(called) else 0, power = sum(called & t) / sum(truth))
}
# W of RUV-III from the proteins observed in every sample, the controls centred (the PRPS implementation's formula).
ruv_w <- function(x, groups, controls, k) {
    complete <- rowSums(!is.finite(x)) == 0
    ctl <- controls[complete]
    if (sum(ctl) <= k) stop("too few complete controls")
    y <- t(x[complete, , drop = FALSE])
    info <- ruv::RUVIII(Y = y, M = groups, ctl = ctl, k = k, return.info = TRUE, inputcheck = FALSE)
    alpha <- info$fullalpha[seq_len(k), , drop = FALSE]
    ac <- alpha[, ctl, drop = FALSE]
    w <- sweep(y, 2, colMeans(y))[, ctl, drop = FALSE] %*% t(ac) %*% solve(ac %*% t(ac))
    dimnames(w) <- list(colnames(x), paste0("W", seq_len(k)))
    w
}
as_matrix <- function(object) {
    q <- object@protein_quant_table
    m <- as.matrix(q[, setdiff(names(q), object@protein_id_column), drop = FALSE]); rownames(m) <- q[[object@protein_id_column]]; m
}
one <- function(n, missing, seed) {
    sim <- simulate(9000L + 100L * n + round(100 * missing) + seed, n, missing)
    x <- sim$x; groups <- sim$design$group
    rows <- list()
    add <- function(arm, k, controls, usable, status, result) rows[[length(rows) + 1L]] <<- data.frame(n = n, missing = missing,
        seed = seed, arm = arm, proteins = nrow(x), status = status, k = k, controls = controls, usable = usable,
        fdp = result[["fdp"]], power = result[["power"]])
    add("no_ruv", NA, NA, NA, "ok", da(x, sim$design, sim$de))
    sweep <- NULL
    quietly(sweep <- tryCatch(findBestNegCtrlPercentage(sim$object, ruv_grouping_variable = "group", verbose = FALSE), error = function(e) NULL))
    if (is.null(sweep) || is.null(sweep$best_k) || is.na(sweep$best_k)) {
        add("today", NA, NA, NA, "refused", c(fdp = NA, power = NA))
        return(do.call(rbind, rows))
    }
    complete <- rowSums(!is.finite(x)) == 0
    selected <- as.logical(sweep$best_control_genes_index); names(selected) <- rownames(x)
    # the same number of true null proteins, taken among the complete ones where there are enough
    null_pool <- which(!sim$de & complete); if (length(null_pool) < sum(selected)) null_pool <- which(!sim$de)
    true_controls <- stats::setNames(logical(nrow(x)), rownames(x))
    true_controls[sample(null_pool, min(length(null_pool), sum(selected)))] <- TRUE
    first <- limma::eBayes(limma::lmFit(x, stats::model.matrix(~ groups)))$p.value[, 2L]
    wide <- stats::setNames(rank(first, ties.method = "first", na.last = "keep") > 0.25 * sum(is.finite(first)), rownames(x))
    wide[is.na(wide)] <- FALSE
    k_today <- as.integer(sweep$best_k)
    k_bound <- min(k_today, calculateAdaptiveMaxK(2L * n))
    correction <- function(arm, controls, k) {
        usable <- sum(controls & complete); corrected <- NULL
        quietly(corrected <- tryCatch(ruvIII_C_Varying(sim$object, ruv_grouping_variable = "group", ruv_number_k = k, ctrl = controls),
            error = function(e) NULL))
        if (is.null(corrected)) return(add(arm, k, sum(controls), usable, "failed", c(fdp = NA, power = NA)))
        add(arm, k, sum(controls), usable, "ok", da(as_matrix(corrected), sim$design, sim$de))
    }
    covariates <- function(arm, controls, k) {
        usable <- sum(controls & complete)
        w <- tryCatch(ruv_w(x, groups, controls, k), error = function(e) NULL)
        if (is.null(w)) return(add(arm, k, sum(controls), usable, "failed", c(fdp = NA, power = NA)))
        result <- tryCatch(da(x, sim$design, sim$de, w), error = function(e) NULL)
        if (is.null(result)) return(add(arm, k, sum(controls), usable, "failed", c(fdp = NA, power = NA)))
        add(arm, k, sum(controls), usable, "ok", result)
    }
    correction("today", selected, k_today); covariates("cov_sel", selected, k_today)
    correction("corr_true", true_controls, k_today); covariates("cov_true", true_controls, k_today)
    correction("corr_wide", wide, k_today); covariates("cov_wide", wide, k_today)
    correction("today_k", selected, k_bound); covariates("cov_sel_k", selected, k_bound); covariates("cov_wide_k", wide, k_bound)
    do.call(rbind, rows)
}
cells <- expand.grid(seed = seq_len(seeds), missing = c(0, 0.25, 0.50), n = c(6L, 12L))
out <- do.call(rbind, parallel::mclapply(seq_len(nrow(cells)), function(i) one(cells$n[[i]], cells$missing[[i]], cells$seed[[i]]),
    mc.cores = cores, mc.preschedule = FALSE))
utils::write.table(out, args[[1L]], sep = "\t", quote = FALSE, row.names = FALSE)
base <- out[out$arm == "no_ruv", c("n", "missing", "seed", "fdp")]; names(base)[[4L]] <- "base"
summary <- do.call(rbind, lapply(split(out, list(out$n, out$missing, out$arm), drop = TRUE), function(d) {
    ok <- merge(d[d$status == "ok", ], base, by = c("n", "missing", "seed"))
    data.frame(n = d$n[[1L]], missing = d$missing[[1L]], arm = d$arm[[1L]], runs = nrow(d), ok = nrow(ok), k = mean(ok$k),
        controls = mean(ok$controls), usable = mean(ok$usable), fdr = mean(ok$fdp), power = mean(ok$power),
        fdr_minus_no_ruv = mean(ok$fdp - ok$base), se_difference = if (nrow(ok) > 1L) stats::sd(ok$fdp - ok$base) / sqrt(nrow(ok)) else NA)
}))
summary <- summary[order(summary$n, summary$missing, summary$arm), ]
utils::write.table(summary, sub("\\.tsv$", "-summary.tsv", args[[1L]]), sep = "\t", quote = FALSE, row.names = FALSE)
print(summary, row.names = FALSE, digits = 3)
