# KIT-114 (e): lipid and metab zeros, kept (today: log2(x + 1), a 0 stays 0) against read as missing (0 -> NA, log2(x)),
# through the lipid DA kernel (runTestsContrastsLipidDA; the metab kernel is the same code). Parameters and the rule are
# fixed in tickets/kit/KIT-114.md (revision 1.2) before any run. Usage: Rscript <this> <out.tsv>
args <- commandArgs(TRUE)
suppressMessages(pkgload::load_all(".", quiet = TRUE))
logger::log_threshold(logger::FATAL)
quietly <- function(expr) suppressWarnings(suppressMessages(utils::capture.output(value <- expr, type = "output")))
simulate <- function(seed, n, censored, m = 500L, pi1 = 0.20, effect = 1.0, d0 = 4, s0 = 0.25) {
    set.seed(seed)
    mu <- rnorm(m, 18, 2)
    sd <- sqrt(d0 * s0^2 / rchisq(m, d0))
    de <- seq_len(m) <= round(pi1 * m)
    delta <- ifelse(de, effect * sample(c(-1, 1), m, replace = TRUE), 0)
    groups <- rep(c("A", "B"), each = n)
    logx <- sapply(seq_len(2 * n), function(j) mu + (groups[[j]] == "B") * delta + rnorm(m, 0, sd))
    # Left-censoring written as 0: P(censored) = logistic((LOD - true log2 value) / 0.5); LOD solved for the target share.
    p_at <- function(lod) mean(stats::plogis((lod - logx) / 0.5))
    lod <- stats::uniroot(function(l) p_at(l) - censored, c(0, 40))$root
    zero <- matrix(stats::runif(length(logx)) < stats::plogis((lod - logx) / 0.5), nrow(logx))
    raw <- 2^logx
    raw[zero] <- 0
    ids <- sprintf("L%04d", seq_len(m))
    dimnames(raw) <- list(ids, sprintf("S%02d", seq_len(2 * n)))
    design <- data.frame(Run = colnames(raw), group = groups, stringsAsFactors = FALSE)
    list(raw = raw, design = design, de = stats::setNames(de, ids), delta = stats::setNames(delta, ids))
}
rows <- list()
for (n in c(3L, 6L, 10L)) for (censored in c(0.10, 0.30)) for (seed in seq_len(20L)) {
    sim <- simulate(5000L + 100L * n + 10L * (censored == 0.30) + seed, n, censored)
    # Both treatments on the same features: at least 2 non-zero values per group, as the QC missingness filter leaves
    # them before DA (correction after the first run failed: limma's robust eBayes refuses features with too few values).
    groups <- sim$design$group
    keep <- rowSums(sim$raw[, groups == "A"] > 0) >= 2 & rowSums(sim$raw[, groups == "B"] > 0) >= 2
    raw <- sim$raw[keep, , drop = FALSE]
    treatments <- list(kept = log2(raw + 1), missing = { y <- raw; y[y == 0] <- NA; log2(y) })
    for (treatment in names(treatments)) {
        value <- NULL
        quietly(value <- runTestsContrastsLipidDA(treatments[[treatment]], "groupB-groupA", sim$design, "~ 0 + group"))
        tab <- value$results[[1L]]
        ids <- rownames(tab)
        truth <- sim$de[ids]
        called <- !is.na(tab$fdr_qvalue) & tab$fdr_qvalue < 0.05
        bias <- mean(tab$logFC[truth] - sim$delta[ids][truth], na.rm = TRUE)
        rel <- mean(abs(tab$logFC[truth] - sim$delta[ids][truth]) / abs(sim$delta[ids][truth]), na.rm = TRUE)
        rows[[length(rows) + 1L]] <- data.frame(n = n, censored = censored, seed = seed, treatment = treatment,
            kept_features = sum(keep), tested = sum(!is.na(tab$P.Value)), called = sum(called), false = sum(called & !truth),
            fdp = if (sum(called)) sum(called & !truth) / sum(called) else 0, power = sum(called & truth) / sum(sim$de[keep]),
            bias = bias, relative_error = rel)
    }
}
out <- do.call(rbind, rows)
utils::write.table(out, args[[1L]], sep = "\t", quote = FALSE, row.names = FALSE)
summary <- do.call(rbind, lapply(split(out, list(out$n, out$censored, out$treatment), drop = TRUE), function(d)
    data.frame(n = d$n[[1]], censored = d$censored[[1]], treatment = d$treatment[[1]], fdr = mean(d$fdp),
        se = stats::sd(d$fdp) / sqrt(nrow(d)), power = mean(d$power), log2fc_bias = mean(d$bias),
        relative_error = mean(d$relative_error), kept_features = mean(d$kept_features), tested = mean(d$tested), seeds = nrow(d))))
summary$controls <- summary$fdr <= 0.05 + 2 * summary$se
summary <- summary[order(summary$n, summary$censored, summary$treatment), ]
utils::write.table(summary, sub("\\.tsv$", "-summary.tsv", args[[1L]]), sep = "\t", quote = FALSE, row.names = FALSE)
print(summary, row.names = FALSE, digits = 3)
