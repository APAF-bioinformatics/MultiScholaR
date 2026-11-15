# KIT-114 (d) part 2 / KIT-118 M2: RUV-III negative controls chosen among all proteins (today) against among proteins
# complete in every sample, through the package's own sweep (findBestNegCtrlPercentage), correction (ruvIII_C_Varying)
# and DA kernel (runTestsContrasts). Parameters and the rule are fixed in tickets/kit/KIT-118.md (M2) before any run.
# The `complete` arm needs the measurement switch (`d_ruv_controls-switch.patch`, beside this file) applied to the
# working tree: the package does not carry it (the rule was not met; KIT-118 M2 result).
# Usage (from the repository root): Rscript <this> <out.tsv>
args <- commandArgs(TRUE)
suppressMessages(pkgload::load_all(".", quiet = TRUE))
logger::log_threshold(logger::FATAL)
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
    # The QC missingness filter's effect: at least half the values in each group.
    keep <- rowSums(is.finite(x[, groups == "control"])) >= n / 2 & rowSums(is.finite(x[, groups == "treated"])) >= n / 2
    x <- x[keep, , drop = FALSE]
    design <- data.frame(Run = runs, group = factor(groups, levels = c("control", "treated")), replicates = seq_len(2 * n),
        tech_rep_group = paste0("s", seq_len(2 * n)), stringsAsFactors = FALSE)
    object_design <- design               # the S4 object's design carries no row names (the package sets them itself)
    rownames(design) <- runs              # the DA kernel's design is keyed by run
    quant <- as.data.frame(c(list(Protein.Group = rownames(x)), as.data.frame(x)), check.names = FALSE, stringsAsFactors = FALSE)
    object <- methods::new("ProteinQuantitativeData", protein_quant_table = quant, protein_id_column = "Protein.Group",
        protein_id_table = data.frame(Protein.Group = rownames(x), stringsAsFactors = FALSE), design_matrix = object_design,
        sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
        args = list(globalParameters = list(workflow_type = "LFQ")))
    list(object = object, de = stats::setNames(de, ids)[rownames(x)], design = design)
}
da <- function(matrix, design, truth) {
    value <- NULL
    quietly(value <- runTestsContrasts(matrix, "grouptreated-groupcontrol", design, "~ 0 + group"))
    tab <- value$results[[1L]]
    called <- !is.na(tab$fdr_qvalue) & tab$fdr_qvalue < 0.05
    t <- truth[rownames(tab)]
    c(called = sum(called), fdp = if (sum(called)) sum(called & !t) / sum(called) else 0, power = sum(called & t) / sum(truth))
}
as_matrix <- function(object) {
    q <- object@protein_quant_table
    m <- as.matrix(q[, setdiff(names(q), object@protein_id_column), drop = FALSE]); rownames(m) <- q[[object@protein_id_column]]; m
}
rows <- list()
for (n in c(6L, 12L)) for (missing in c(0, 0.25, 0.50)) for (seed in seq_len(20L)) {
    sim <- simulate(9000L + 100L * n + round(100 * missing) + seed, n, missing)
    complete <- sum(ruvCompleteControlRows(sim$object))
    base <- da(as_matrix(sim$object), sim$design, sim$de)
    rows[[length(rows) + 1L]] <- data.frame(n = n, missing = missing, seed = seed, arm = "no_ruv", proteins = nrow(sim$object@protein_quant_table),
        complete = complete, status = "ok", percentage = NA, k = NA, controls = NA, usable = NA, fdp = base[["fdp"]], power = base[["power"]])
    for (arm in c("all", "complete")) {
        options(MultiScholaR.ruv_controls_complete_only = identical(arm, "complete"))
        sweep <- NULL; failure <- NULL
        quietly(sweep <- tryCatch(findBestNegCtrlPercentage(sim$object, ruv_grouping_variable = "group", verbose = FALSE),
            error = function(e) { failure <<- conditionMessage(e); NULL }))
        status <- if (is.null(sweep) || is.null(sweep$best_k) || is.na(sweep$best_k)) "refused" else "ok"
        result <- c(fdp = NA_real_, power = NA_real_); usable <- NA; controls <- NA
        if (identical(status, "ok")) {
            controls <- sum(sweep$best_control_genes_index)
            usable <- sum(ruvCompleteControlRows(sim$object)[which(sweep$best_control_genes_index)])
            corrected <- NULL
            quietly(corrected <- tryCatch(ruvIII_C_Varying(sim$object, ruv_grouping_variable = "group",
                ruv_number_k = sweep$best_k, ctrl = sweep$best_control_genes_index), error = function(e) { failure <<- conditionMessage(e); NULL }))
            if (is.null(corrected)) status <- "correction_failed" else {
                r <- da(as_matrix(corrected), sim$design, sim$de); result <- c(fdp = r[["fdp"]], power = r[["power"]])
            }
        }
        rows[[length(rows) + 1L]] <- data.frame(n = n, missing = missing, seed = seed, arm = arm, proteins = nrow(sim$object@protein_quant_table),
            complete = complete, status = status, percentage = if (identical(status, "refused")) NA else sweep$best_percentage,
            k = if (identical(status, "refused")) NA else sweep$best_k, controls = controls, usable = usable,
            fdp = result[["fdp"]], power = result[["power"]])
    }
    options(MultiScholaR.ruv_controls_complete_only = FALSE)
    utils::write.table(do.call(rbind, rows), args[[1L]], sep = "\t", quote = FALSE, row.names = FALSE)
}
out <- do.call(rbind, rows)
summary <- do.call(rbind, lapply(split(out, list(out$n, out$missing, out$arm), drop = TRUE), function(d) {
    ok <- d[d$status == "ok", ]
    data.frame(n = d$n[[1]], missing = d$missing[[1]], arm = d$arm[[1]], runs = nrow(d), ok = nrow(ok), refused = sum(d$status == "refused"),
        correction_failed = sum(d$status == "correction_failed"), proteins = mean(d$proteins), complete = mean(d$complete),
        k = mean(ok$k), controls = mean(ok$controls), usable = mean(ok$usable), fdr = mean(ok$fdp),
        se = if (nrow(ok) > 1) stats::sd(ok$fdp) / sqrt(nrow(ok)) else NA, power = mean(ok$power))
}))
summary <- summary[order(summary$n, summary$missing, summary$arm), ]
utils::write.table(summary, sub("\\.tsv$", "-summary.tsv", args[[1L]]), sep = "\t", quote = FALSE, row.names = FALSE)
print(summary, row.names = FALSE, digits = 3)
