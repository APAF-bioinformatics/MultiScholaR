# KIT-114 (a): qvalue() against Benjamini-Hochberg through the proteomics DA kernel (runTestsContrasts), on centred
# simulated proteins with known truth. Parameters and the decision rule are fixed in tickets/kit/KIT-114.md
# (revision 1.2) before any run. Usage (from the repository root): Rscript <this> <out.tsv>
args <- commandArgs(TRUE)
suppressMessages(pkgload::load_all(".", quiet = TRUE))
quietly <- function(expr) suppressWarnings(suppressMessages(utils::capture.output(value <- expr, type = "output")))
simulate <- function(seed, n, pi1, effect, m = 2000L, d0 = 4, s0 = 0.25) {
    set.seed(seed)
    mu <- rnorm(m, 20, 2)
    sd <- sqrt(d0 * s0^2 / rchisq(m, d0))
    de <- seq_len(m) <= round(pi1 * m)
    delta <- ifelse(de, effect * sample(c(-1, 1), m, replace = TRUE), 0)
    groups <- rep(c("A", "B"), each = n)
    x <- sapply(seq_len(2 * n), function(j) mu + (groups[[j]] == "B") * delta + rnorm(m, 0, sd))
    ids <- sprintf("P%05d", seq_len(m))
    dimnames(x) <- list(ids, sprintf("S%02d", seq_len(2 * n)))
    design <- data.frame(Run = colnames(x), group = groups, replicates = seq_len(2 * n), stringsAsFactors = FALSE)
    rownames(design) <- design$Run
    list(x = x, design = design, de = stats::setNames(de, ids))
}
rows <- list()
for (n in c(3L, 4L, 6L, 10L)) for (pi1 in c(0.05, 0.20)) for (effect in c(0.5, 1.0)) for (seed in seq_len(20L)) {
    sim <- simulate(1000L * n + 100L * (pi1 == 0.20) + 10L * (effect == 1.0) + seed, n, pi1, effect)
    value <- NULL
    quietly(value <- runTestsContrasts(sim$x, "groupB-groupA", sim$design, "~ 0 + group"))
    tab <- value$results[[1L]]
    ids <- rownames(tab)
    truth <- sim$de[ids]
    for (column in c("fdr_qvalue", "fdr_value_bh_adjustment")) {
        called <- !is.na(tab[[column]]) & tab[[column]] < 0.05
        rows[[length(rows) + 1L]] <- data.frame(n = n, pi1 = pi1, effect = effect, seed = seed, column = column,
            called = sum(called), false = sum(called & !truth), true_de = sum(truth),
            fdp = if (sum(called)) sum(called & !truth) / sum(called) else 0, power = sum(called & truth) / sum(truth))
    }
}
out <- do.call(rbind, rows)
utils::write.table(out, args[[1L]], sep = "\t", quote = FALSE, row.names = FALSE)
summary <- do.call(rbind, lapply(split(out, list(out$n, out$pi1, out$effect, out$column), drop = TRUE), function(d)
    data.frame(n = d$n[[1]], pi1 = d$pi1[[1]], effect = d$effect[[1]], column = d$column[[1]], fdr = mean(d$fdp),
        se = stats::sd(d$fdp) / sqrt(nrow(d)), power = mean(d$power), seeds = nrow(d))))
summary$controls <- summary$fdr <= 0.05 + 2 * summary$se
summary <- summary[order(summary$n, summary$pi1, summary$effect, summary$column), ]
utils::write.table(summary, sub("\\.tsv$", "-summary.tsv", args[[1L]]), sep = "\t", quote = FALSE, row.names = FALSE)
print(summary, row.names = FALSE, digits = 3)
