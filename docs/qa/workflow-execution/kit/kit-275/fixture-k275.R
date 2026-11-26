# KIT-275 G1: the exactness fixtures (section 8.4 of the ticket) and the checks E6 and E10. Evidence code under docs/qa: no package code.
# Loads sim-k275.R as a library (`options(k275.library = TRUE)`): the generator, the unmoderated reading and the estimand functions are
# the ones the counted jobs use.
#
# Usage: Rscript fixture-k275.R <designs> <first> <last> <out dir>      designs: a comma list of D3, D4, D6 (the cells X3, X4, X6)
# Environment: K275_SMOKE=1 writes the smoke flag to the table (the scorer refuses such a table for an outcome line).
# Writes g1.tsv (design, kind PC2 | NC2, r, rho, smoke), checks-g1.tsv (E6, E10, and the per-data-set identity of the reading) and
# provenance-g1.json. Scenario N0 (no unwanted variation, no anomaly), one assay of 2,000 features, 20 data sets per design.

options(k275.library = TRUE)
args_fixture <- commandArgs(TRUE)
stopifnot(length(args_fixture) == 4L)
fixture_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])))
source(file.path(fixture_dir, "sim-k275.R"), local = FALSE)
designs <- strsplit(args_fixture[[1L]], ",", fixed = TRUE)[[1L]]
first_rep <- as.integer(args_fixture[[2L]])
last_rep <- as.integer(args_fixture[[3L]])
out_dir <- args_fixture[[4L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
smoke <- identical(Sys.getenv("K275_SMOKE"), "1")
cell_of <- c(D3 = "X3", D4 = "X4", D6 = "X6")
stopifnot(all(designs %in% names(cell_of)), !is.na(first_rep), !is.na(last_rep), first_rep <= last_rep)

#' The generator's per-feature residual variance (the truth): the first two draws of `generate` (the levels, then the variances' chi-squares),
#' made with the same seed.
truthVariance <- function(cellname, r) {
    p <- CELLS[[cellname]]
    set.seed(275000000 + 100000 * p$idx + 10000 * match("N0", SCENARIOS) + r, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    n_feat <- p$per_assay * length(p$assays)
    invisible(stats::rnorm(n_feat, p$level, p$sd))
    x <- 10 / stats::rchisq(n_feat, 10)
    0.25 * x / stats::median(x)
}

g1 <- list()
for (design in designs) {
    cellname <- cell_of[[design]]
    n <- CELLS[[cellname]]$per
    for (r in seq(first_rep, last_rep)) {
        sim <- generate(cellname, "N0", r)
        features <- sim$feat_ids
        sigma2 <- truthVariance(cellname, r)
        stopifnot(length(sigma2) == length(features))
        # the genuine data without the sample loading (a part of the generator that normalisation removes): the oracle of N0
        values <- sim$oracle[features, , drop = FALSE]
        a_run <- sim$runs[[length(sim$runs)]]
        # NC2: the reading on the genuine data; expected 1
        genuine <- unmoderatedReading(values, sim$design, features, NULL)
        rho_nc2 <- sqrt(genuine$sum / sum(sigma2 * (2 / n)))
        # PC2: sample a's row replaced by the mean of its mates' values; expected sqrt((2n - 3) / (2n - 1))
        made <- values
        mates <- sim$runs[sim$design$group == "B" & sim$runs != a_run]
        made[, a_run] <- rowMeans(values[, mates, drop = FALSE])
        pseudo <- unmoderatedReading(made, sim$design, features, NULL)
        rho_pc2 <- sqrt(pseudo$sum / sum(sigma2 * (1 / (n - 1) + 1 / n)))
        g1[[length(g1) + 1L]] <- data.frame(design = design, kind = c("NC2", "PC2"), r = r,
            rho = if (smoke) NA_real_ else c(rho_nc2, rho_pc2), smoke = smoke, stringsAsFactors = FALSE)   # the smoke send writes no reading
        # E6 on the first data set of each design
        if (r == first_rep) {
            fit <- quiet(fitLipidDaAssayModel(made, as.data.frame(sim$design), "~ 0 + group", "Run"))
            result <- quiet(applyLipidDaContrasts(fit, "groupB-groupA", NA, TRUE, TRUE))
            tab <- result$results[["groupB-groupA"]]
            eb <- result$fit.eb
            pos <- match(rownames(tab), rownames(eb$coefficients))
            direct <- eb$stdev.unscaled[pos, 1L] * sqrt(eb$s2.post[pos])
            moderated <- abs(tab$logFC / tab$t)
            dev_mod <- max(abs(moderated - direct), na.rm = TRUE)
            addCheck("E6", cellname, r, "moderated SE |logFC / t| equals stdev.unscaled x sqrt(s2.post)", signif(dev_mod, 4), 1e-8, dev_mod <= 1e-8)
            kernel_v <- sum(eb$sigma[pos]^2 * eb$stdev.unscaled[pos, 1L]^2)
            dev_un <- abs(kernel_v - pseudo$sum) / pseudo$sum
            addCheck("E6", cellname, r, "the script's lm.fit variance equals the kernel's sigma^2 x stdev.unscaled^2", signif(dev_un, 4), 1e-8, dev_un <= 1e-8)
            ctx <- list(tier = "G1", design = design, cell = cellname, omic = "lipid", r = r, per = n, smoke = smoke, shares = plantedShares(cellname, r))
            lx <- lxRow(ctx, list(m_out = list(LCMS_Pos = made)), sim, "V1g")
            dev_lx <- abs(lx[["se2_unmod_y"]] - pseudo$sum) / pseudo$sum
            addCheck("E6", cellname, r, "LX's reading on the hand-made pseudo-replicate matrix equals PC2's", signif(dev_lx, 4), 1e-10, dev_lx <= 1e-10)
        }
    }
}

# E10: kappa_c, PR and the dot product of pi on a constructed matrix with a known planted_delta and known controls (the values are derived by hand)
delta <- c(F1 = 1, F2 = 2, F3 = 0, F4 = 3, LBL_001 = 0)
kappa <- kappaControls(delta, c("F1", "F3"))
addCheck("E10", "constructed", 0L, "kappa_c = mean(delta^2 over the controls) / mean(delta^2 over the affected rows) = 0.5 / 2.8", signif(kappa, 10), 0.5 / 2.8,
    isTRUE(all.equal(kappa, 0.5 / 2.8, tolerance = 1e-12)))
M <- cbind(R1 = c(1, 0), R2 = c(-1, 0), R3 = c(0.5, 0.5), R4 = c(1, 0), R5 = c(-1, 0), R6 = c(4, 3))
rownames(M) <- c("F1", "F2")
design_c <- data.frame(Run = colnames(M), group = rep(c("A", "B"), each = 3L), stringsAsFactors = FALSE)
pr <- pseudoReplicateRatio(M, design_c, "R6")
addCheck("E10", "constructed", 0L, "PR = d_a / median(d_s) = sqrt(12.5) / sqrt(1.25) = sqrt(10)", signif(pr, 10), sqrt(10), isTRUE(all.equal(pr, sqrt(10), tolerance = 1e-12)))
shares_c <- list(g = list(x = c(F1 = 1, F2 = -1), delta = delta), p = list(x = c(F1 = 0, F2 = 0)), t = list(x = c(F1 = 0, F2 = 0)), two = list(x = c(F1 = 0, F2 = 0)))
numbers <- daNumbers(data.frame(feature = c("F1", "F2"), logFC = c(1.5, -0.5), t = c(3, -1), p = c(0.01, 0.2), q = c(0.02, 0.4), stringsAsFactors = FALSE), shares_c)
addCheck("E10", "constructed", 0L, "dot_g = sum(x * logFC) = 1 * 1.5 + (-1) * (-0.5) = 2; SE^2 sum = 0.25 + 0.25 = 0.5", signif(numbers$dots[["g"]], 10), 2,
    isTRUE(all.equal(numbers$dots[["g"]], 2)) && isTRUE(all.equal(numbers$se2, 0.5)) && identical(numbers$any, 1L))

g1_table <- do.call(rbind, g1)
utils::write.table(g1_table, file.path(out_dir, "g1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
checks_table <- do.call(rbind, checks_env$rows)
utils::write.table(checks_table, file.path(out_dir, "checks-g1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
jsonlite::write_json(list(job = "g1", designs = designs, data_sets = c(first_rep, last_rep), smoke = smoke, commit = tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    error = function(e) NA), r_version = R.version.string, script_sha256 = sha(file.path(fixture_dir, "sim-k275.R")),
    fixture_sha256 = sha(file.path(fixture_dir, "fixture-k275.R")), started = format(started_at, "%Y-%m-%dT%H:%M:%S%z"), finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
    file.path(out_dir, "provenance-g1.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat(sprintf("K275 G1: %d data-set readings written, %d checks, %d failed\n", nrow(g1_table), nrow(checks_table), sum(!checks_table$pass)))
