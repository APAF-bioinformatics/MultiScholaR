# KIT-201, for KIT-188's E2 (VAL-132 (2)): on this commit, do proteomics' two "RUV off" tests give the same table?
# Skip: the legacy test, runTestsContrasts (TREAT at 0 with trend and robust after KIT-201). Declared `none`: the declared
# designs' test, runProtRuvModelContrasts(exact = TRUE). Kernel level only (the same matrix into both); the routes'
# steps before the DA (which matrix each tests) are not compared here. Simulated data only.
# Usage, from the package root: Rscript skip-none.R <out dir>
suppressPackageStartupMessages(pkgload::load_all(".", export_all = TRUE, helpers = FALSE, quiet = TRUE))
out <- commandArgs(trailingOnly = TRUE)[[1L]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
data <- function(seed, missing, groups = c("A", "B"), per = 4L, n = 2000L) {
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    group <- rep(groups, each = per)
    level <- stats::rnorm(n, 22, 2)
    s2 <- 4 * 0.25^2 / stats::rchisq(n, 4) * exp(-0.541 * (level - 22) / 2)
    beta <- ifelse(seq_len(n) <= 200L, sample(c(-1, 1), n, replace = TRUE), 0)
    y <- matrix(level, n, length(group)) + outer(beta, as.numeric(group == "B")) +
        matrix(stats::rnorm(n * length(group)), n, length(group)) * sqrt(s2)
    if (missing) y[matrix(stats::runif(length(y)), nrow(y)) < stats::plogis(19.7 - y)] <- NA_real_
    dimnames(y) <- list(sprintf("P%05d", seq_len(n)), sprintf("S%02d", seq_along(group)))
    design <- data.frame(Run = colnames(y), group = group, replicates = seq_along(group), stringsAsFactors = FALSE)
    rownames(design) <- design$Run
    list(y = y, design = design)
}
rows <- list()
for (case in list(list(seed = 1L, missing = FALSE, per = 4L), list(seed = 2L, missing = TRUE, per = 4L),
    list(seed = 3L, missing = TRUE, per = 3L), list(seed = 4L, missing = TRUE, per = 6L))) {
    d <- data(case$seed, case$missing, per = case$per)
    contrast <- "groupB-groupA"
    skip <- suppressMessages(suppressWarnings(runTestsContrasts(d$y, contrast, d$design, "~ 0 + group", treat_lfc_cutoff = 0,
        eBayes_trend = TRUE, eBayes_robust = TRUE)))$results[[contrast]]
    none <- suppressMessages(suppressWarnings(runProtRuvModelContrasts(d$y, contrast, d$design, "~ 0 + group", treat_lfc_cutoff = 0,
        eBayes_trend = TRUE, eBayes_robust = TRUE, exact = TRUE)))$results[[contrast]]
    ids <- rownames(skip)
    at <- match(ids, rownames(none))
    diff <- function(column) {
        a <- skip[[column]]
        b <- none[[column]][at]
        if (is.null(a) || is.null(b)) return(NA_real_)
        both <- !is.na(a) & !is.na(b)
        if (!any(both)) return(NA_real_)
        max(abs(a[both] - b[both]))
    }
    rows[[length(rows) + 1L]] <- data.frame(seed = case$seed, missing = case$missing, per_group = case$per,
        rows_skip = nrow(skip), rows_none = nrow(none), same_ids = setequal(ids, rownames(none)),
        na_pattern_same = identical(is.na(skip$raw_pvalue), is.na(none$raw_pvalue[at])),
        max_logFC = diff("logFC"), max_raw_p = diff("raw_pvalue"), max_q = diff("fdr_qvalue"),
        max_bh = diff("fdr_value_bh"), columns_same = identical(sort(names(skip)), sort(names(none))))
}
table <- do.call(rbind, rows)
utils::write.table(table, file.path(out, "skip-none.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(table)
