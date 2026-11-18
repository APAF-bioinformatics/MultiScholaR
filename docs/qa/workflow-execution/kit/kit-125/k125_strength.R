# KIT-125, bounded check of the unwanted variation's strength (decision M21; method: declaration.md, section
# "Bounded check of the unwanted variation's strength", written before this script was run). Evidence only.
#   Rscript <this folder>/k125_strength.R <out dir> <DIA_LIMPA data dir> <matrices per source, declared: 10>
# from the root of the pinned work tree (0b6fcb6d), R_LIBS and threads as k125_common.R says.
# Real matrices: the batch-aligned null splits 1..n of PXD024043, through k125_dia_limpa.R's own top-level code
# (its preparation, pool, batch and `quantify`; every expression but the one that sources k125_common.R and the
# closing k125Run call is evaluated as it stands). Simulated matrices: `simulateData` and `quantify` of
# kit-124/sim/s2_limpa.R on the pinned tree, evaluated as they stand, cells C0 and C5, runs 1..n under the
# simulator's seeds; the proteins without a group effect.
# Writes strength-matrices.tsv (one row per matrix), strength-summary.tsv (median and range over the matrices of a
# source) and meta.txt.

file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
here <- dirname(normalizePath(file))
source(file.path(here, "k125_common.R"))
real_code <- parse(file.path(here, "k125_dia_limpa.R"))
stopifnot(identical(real_code[[1L]][[1L]], as.name("local")), identical(real_code[[length(real_code)]][[1L]], as.name("k125Run")))
for (expression in real_code[-c(1L, length(real_code))]) eval(expression, globalenv())    # sets args, loads the tree, pool, batch, quantify
n <- args$splits
dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(args$out_dir, "strength-matrices.tsv"))) stop("the out dir already holds a result: one attempt per out dir")
started <- Sys.time()

S2_FILE <- "docs/qa/workflow-execution/kit/kit-124/sim/s2_limpa.R"
S2_NAMES <- c("N_PROTEINS", "PRECURSORS", "MU_RANGE", "SIGMA", "DPC_SLOPE", "EFFECT", "RESIDUAL_VAR", "UNIT_SD",
    "cellTable", "simulateData", "quantify")
s2 <- new.env(parent = globalenv())
for (expression in parse(S2_FILE)) {
    if (is.call(expression) && identical(expression[[1L]], as.name("<-")) && is.name(expression[[2L]]) &&
        as.character(expression[[2L]]) %in% S2_NAMES) eval(expression, s2)
}
stopifnot(all(S2_NAMES %in% ls(s2)))

# The measures of one matrix E (proteins x 12 samples), factor f (one value per sample), group = B-against-A indicator.
measure <- function(E, f, group, seed) {
    stopifnot(!anyNA(E), length(f) == ncol(E), length(group) == ncol(E))
    X <- E - rowMeans(E)
    total <- rowSums(X^2)
    X <- X[total > 0, , drop = FALSE]; total <- total[total > 0]
    shares <- function(f) {
        centred <- f - mean(f)
        explained <- as.numeric(X %*% centred)^2 / sum(centred^2)
        list(r2 = explained / total, pooled = sum(explained) / sum(total))
    }
    of_factor <- shares(f)
    set.seed(100L * seed + 98L)
    of_permuted <- shares(sample(f))
    pca <- stats::prcomp(t(X), center = FALSE, scale. = FALSE)
    variance <- pca$sdev^2 / sum(pca$sdev^2)
    q <- stats::quantile(of_factor$r2, c(0.25, 0.5, 0.75), names = FALSE)
    out <- data.frame(n_proteins = nrow(X), r2_q25 = q[[1L]], r2_median = q[[2L]], r2_q75 = q[[3L]],
        r2_above_010 = mean(of_factor$r2 > 0.10), r2_above_030 = mean(of_factor$r2 > 0.30), pooled_share = of_factor$pooled,
        permuted_r2_median = stats::median(of_permuted$r2), permuted_pooled_share = of_permuted$pooled)
    for (i in 1:5) out[[sprintf("pc%d_variance", i)]] <- variance[[i]]
    for (i in 1:5) out[[sprintf("pc%d_r2_factor", i)]] <- stats::cor(pca$x[, i], f)^2
    out$factor_group_cor <- stats::cor(f, group)
    out
}

rows <- list()
# ---- real: batch-aligned splits 1..n -----------------------------------------------------------------------------
recorded <- utils::read.delim(file.path(here, "runs", "dia_limpa", "splits.tsv"), stringsAsFactors = FALSE, quote = "")
first_batch <- sort(unique(batch))[[1L]]
for (split in seq_len(n)) {
    seed <- splitSeed("dia_limpa", "batch", split)
    units <- drawSplit(pool, batch, "batch", seed)
    was <- recorded[recorded$split_type == "batch" & recorded$split == split, ]
    stopifnot(nrow(was) == 1L, was$seed == seed, identical(paste(units$A, collapse = ","), was$group_A),
        identical(paste(units$B, collapse = ","), was$group_B))
    samples <- c(units$A, units$B)                                   # from here: the first lines of `prepare`
    y <- prepared$y[, samples, drop = FALSE]
    seen <- is.finite(y)
    keep <- rowSums(seen[, 1:6]) >= 2L & rowSums(seen[, 7:12]) >= 2L
    runs <- sprintf("run_%03d", 1:12)
    design <- data.frame(Run = runs, group = rep(c("A", "B"), each = 6L), replicates = rep(1:6, 2L), stringsAsFactors = FALSE)
    design$tech_rep_group <- paste(design$group, design$replicates, sep = "_")
    quantified <- quantify(y[keep, , drop = FALSE], prepared$protein[keep], prepared$sequence[keep], prepared$precursor[keep], design)
    E <- quantified@args$limpa_dpc_quant_results$quantified_elist$E
    stopifnot(identical(colnames(E), runs))
    f <- as.numeric(batch[match(samples, pool)] == first_batch)
    stopifnot(sum(f) == 6)
    rows[[length(rows) + 1L]] <- cbind(data.frame(source = "real, batch-aligned", matrix = split, seed = seed, factor = "culture batch"),
        measure(E, f, as.numeric(design$group == "B"), seed))
    cat(sprintf("%s real split %d done\n", format(Sys.time(), "%H:%M:%S"), split))
}
# ---- simulated: S2's cells C0 and C5, runs 1..n ------------------------------------------------------------------
cells <- s2$cellTable()
for (name in c("C0", "C5")) {
    cell <- cells[[name]]
    for (run in seq_len(n)) {
        seed <- 100000L + 10000L * cell$index + run                  # s2_limpa.R's seed rule
        sim <- s2$simulateData(cell, seed)
        quantified <- s2$quantify(sim)
        E <- quantified@args$limpa_dpc_quant_results$quantified_elist$E
        stopifnot(identical(colnames(E), sim$design$Run), all(rownames(E) %in% names(sim$truth)))
        E <- E[sim$truth[rownames(E)] == 0, , drop = FALSE]
        rows[[length(rows) + 1L]] <- cbind(data.frame(source = paste("simulated,", name), matrix = run, seed = seed, factor = "true factor score"),
            measure(E, sim$scores[sim$design$Run, 1L], as.numeric(sim$design$group == "B"), seed))
        cat(sprintf("%s simulated %s run %d done\n", format(Sys.time(), "%H:%M:%S"), name, run))
    }
}
table <- do.call(rbind, rows)
measures <- setdiff(names(table), c("source", "matrix", "seed", "factor"))
summary <- do.call(rbind, lapply(unique(table$source), function(s) do.call(rbind, lapply(measures, function(m) {
    v <- table[[m]][table$source == s]
    data.frame(source = s, measure = m, matrices = length(v), median = stats::median(v), min = min(v), max = max(v))
}))))
round4 <- function(d) { numeric <- vapply(d, is.double, logical(1L)); d[numeric] <- lapply(d[numeric], signif, digits = 4L); d }
utils::write.table(round4(table), file.path(args$out_dir, "strength-matrices.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(round4(summary), file.path(args$out_dir, "strength-summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
git <- function(...) paste(system2("git", c(...), stdout = TRUE, stderr = FALSE), collapse = " ")
sha <- function(f) unname(tools::sha256sum(f))
writeLines(c("KIT-125 bounded check of the unwanted variation's strength (decision M21)",
    sprintf("package tree commit %s; changes under R/, tests/, src/: %s", git("rev-parse", "HEAD"),
        nzchar(git("status", "--porcelain", "--untracked-files=no", "--", "R", "tests", "src"))),
    sprintf("sha256: k125_strength.R %s; k125_dia_limpa.R %s; k125_common.R %s; s2_limpa.R %s", sha(file), sha(file.path(here, "k125_dia_limpa.R")),
        sha(file.path(here, "k125_common.R")), sha(S2_FILE)),
    sprintf("input (sha256): %s %s", basename(input), prepared$input_sha256),
    sprintf("matrices per source: %d; threads: OMP_NUM_THREADS=%s OPENBLAS_NUM_THREADS=%s; one R process", n, Sys.getenv("OMP_NUM_THREADS"), Sys.getenv("OPENBLAS_NUM_THREADS")),
    sprintf("%s; packages: %s", R.version.string, paste(vapply(c("MultiScholaR", "limpa", "fastlimpadpc", "limma"), function(p) paste(p, as.character(utils::packageVersion(p))), ""), collapse = "; ")),
    sprintf("DPC engine: %s", QUANT$engine),
    sprintf("started %s; finished %s", format(started, "%Y-%m-%d %H:%M:%S %Z"), format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))),
    file.path(args$out_dir, "meta.txt"))
print(round4(summary), row.names = FALSE)
