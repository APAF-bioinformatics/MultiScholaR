# KIT-101 revision 2.3, F-19 diagnosis: the prot_dia_limpa fixture's ground truth (no package code).
# Rebuilds the lane's protein matrix and precursor table with the committed generator's own functions (the file is
# parsed and its function definitions evaluated; its top-level generate_all_fixtures() call is not run), checks the
# rebuilt report against the committed report.tsv (sha256), and measures what the fixture carries per group:
# the per-sample scale terms (kitProteinMatrix sample_scale, the report's run_scale), missingness, peptides per protein.
# Then the "oracle" analyses on the truth: limma on the true protein log2 matrix, with and without median centring.
# Usage (repository root): Rscript docs/qa/workflow-execution/kit/kit-101/diagnosis/f19_fixture.R

suppressPackageStartupMessages({ library(limma) })
root <- normalizePath(".")
gen_file <- file.path(root, "tests/testdata/e2e/generate_fixtures.R")
env <- new.env(parent = globalenv())
exprs <- parse(gen_file)
for (e in exprs) {
    is_def <- is.call(e) && identical(e[[1L]], as.name("<-")) && is.call(e[[3L]]) && identical(e[[3L]][[1L]], as.name("function"))
    is_root <- is.call(e) && identical(e[[1L]], as.name("<-")) && identical(e[[2L]], as.name("FIXTURE_ROOT"))
    if (is_def) eval(e, env)
}
assign("FIXTURE_ROOT", file.path(root, "tests/testdata/e2e"), env)
manifest <- jsonlite::read_json(file.path(root, "tests/testdata/e2e/manifest.json"), simplifyVector = FALSE)
lane <- Filter(function(l) identical(l$lane_id, "prot_dia_limpa"), manifest$lanes)[[1L]]

# 1. The generator reproduces the committed report (so the truth below is the fixture's truth).
out <- file.path(tempdir(), "f19_regen"); dir.create(out, showWarnings = FALSE)
invisible(capture.output(env$generateProtDiannFixture(lane, out)))
regen <- digest::digest(file = file.path(out, "report.tsv"), algo = "sha256")
cat("regenerated report.tsv sha256:", regen, "\nmanifest sha256:            ", lane$files_sha256$report.tsv,
    "\nmatch:", identical(regen, lane$files_sha256$report.tsv), "\n\n")

# 2. The truth: protein matrix (linear), its sample scale terms, the differential set.
kit <- env$kitProteinTools()
lm_ <- env$protFixtureMatrix(lane, kit)
design <- lm_$design
recipe <- lm_$recipe
# kitProteinMatrix's own draws in order, to recover sample_scale (the matrix's per-sample shift).
truth <- withr::with_seed(recipe$seed, {
    n <- recipe$samples; p <- recipe$proteins
    group <- rep(recipe$groups, length.out = n)
    batch <- paste0("B", ((seq_len(n) - 1L) %/% 10L) + 1L)
    base <- stats::rnorm(p, mean = 20, sd = 2.2)
    differential <- seq_len(p) <= round(p * recipe$differential_fraction)
    effect <- ifelse(differential, stats::rnorm(p, 0, 1.2), 0)
    batch_effect <- stats::rnorm(length(unique(batch)), 0, 0.15)
    sample_scale <- stats::rnorm(n, 0, 0.25)
    label <- sprintf("%s_%d", group, stats::ave(seq_along(group), group, FUN = seq_along))
    list(group = group, label = label, effect = effect, differential = differential, sample_scale = sample_scale, base = base)
})
at <- match(design$sample, truth$label)
sample_scale <- truth$sample_scale[at]
names(sample_scale) <- design$sample
cat("kitProteinMatrix sample_scale (log2) by sample:\n"); print(round(sample_scale, 3))
cat(sprintf("mean KO - mean WT of sample_scale: %.3f log2\n", mean(sample_scale[design$group == "KO"]) - mean(sample_scale[design$group == "WT"])))

# The report's run_scale (Precursor.Normalised = Precursor.Quantity * run_scale): the generator's draw order.
report <- utils::read.delim(file.path(root, "tests/testdata/e2e/prot_dia_limpa/report.tsv"), stringsAsFactors = FALSE)
ratio <- tapply(log2(report$Precursor.Normalised / report$Precursor.Quantity), report$Run, median)
ratio <- ratio[design$sample]
cat("\nreport run_scale (median log2 Precursor.Normalised / Precursor.Quantity) by run:\n"); print(round(ratio, 3))
cat(sprintf("mean KO - mean WT of run_scale: %.3f log2\n", mean(ratio[design$group == "KO"]) - mean(ratio[design$group == "WT"])))

# Observed per-run median log2 Precursor.Normalised (the quantity the import reads), and per group.
med <- tapply(log2(report$Precursor.Normalised), report$Run, median)[design$sample]
cat("\nper-run median log2(Precursor.Normalised):\n"); print(round(med, 3))

# 3. Missingness and peptides per protein, by planted status.
ids <- sprintf("P%05d", seq_len(recipe$proteins))
planted <- ids[truth$differential]
pairs <- unique(report[c("Protein.Group", "Stripped.Sequence")])
pep_per_protein <- table(factor(pairs$Protein.Group, levels = ids))
cat("\npeptides per protein (observed in the report):\n"); print(table(pep_per_protein))
prec <- unique(report[c("Protein.Group", "Precursor.Id")])
n_prec <- table(factor(prec$Protein.Group, levels = ids))
cells <- table(factor(report$Protein.Group, levels = ids), factor(report$Run, levels = design$sample))
miss_frac <- 1 - cells / as.vector(n_prec)
miss_by_group <- data.frame(protein = ids, planted = ids %in% planted,
    WT = rowMeans(miss_frac[, design$group == "WT"]), KO = rowMeans(miss_frac[, design$group == "KO"]))
cat("\nprecursor missing fraction, mean over proteins, by group and planted status:\n")
print(aggregate(cbind(WT, KO) ~ planted, miss_by_group, mean))
cat("protein-run cells absent (no precursor at all):", sum(cells == 0), "of", length(cells), "\n")
cat("missing protein-run cells by group (truth matrix lm$missing):",
    sprintf("WT %d, KO %d", sum(lm_$missing[, design$group == "WT"]), sum(lm_$missing[, design$group == "KO"])), "\n")

# 4. Oracle analyses on the true protein matrix (no missing, no peptide layer).
Y <- log2(lm_$values); rownames(Y) <- ids; colnames(Y) <- design$sample
g <- factor(design$group, levels = c("WT", "KO"))
X <- model.matrix(~ g)
oracle <- function(M, label) {
    fit <- eBayes(lmFit(M, X), trend = TRUE, robust = TRUE)
    tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
    sig <- rownames(tt)[tt$adj.P.Val < 0.05]
    cat(sprintf("%-46s sig %3d  planted %2d/40  false positives %2d  median logFC of non-planted %.3f\n", label,
        length(sig), sum(sig %in% planted), sum(!sig %in% planted), median(tt$logFC[!rownames(tt) %in% planted])))
}
cat("\nOracle DA on the true protein matrix (limma, BH q < 0.05):\n")
oracle(Y, "true log2 matrix, no normalisation")
oracle(sweep(Y, 2, apply(Y, 2, median)), "true log2 matrix, median centred")
oracle(Y - outer(rep(1, nrow(Y)), sample_scale), "true log2 matrix minus sample_scale")
Yr <- Y + outer(rep(1, nrow(Y)), ratio)
oracle(Yr, "true + run_scale (what Precursor.Normalised carries)")
oracle(sweep(Yr, 2, apply(Yr, 2, median)), "true + run_scale, median centred")
