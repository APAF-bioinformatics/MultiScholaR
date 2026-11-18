# KIT-125, transcriptomics: null splits of the wild-type replicates of Schurch et al. 2016 (ENA PRJEB5348; gene
# counts from github.com/bartongroup/profDGE48, commit 5ff79149, Preprocessed_data/WT_countdata.tar.gz) through the
# code paths of KIT-124's simulator S3 (kit-124/sim/s3_trx.R). Usage and outputs: k125_common.R.
#
# Data, in place of S3's `generate`:
#   files   <data dir>/WT_countdata.tar.gz (48 files, gene and count; one biological replicate each, the 7 lanes
#           summed by the authors) and <data dir>/exclude.lst (the authors' list of bad replicates)
#   genes   the count rows; the five summary rows of the counting tool (no_feature, ambiguous, too_low_aQual,
#           not_aligned, alignment_not_unique) are dropped
#   pool    the 48 wild-type replicates minus the 6 the authors exclude (42 biological units; no batch recorded)
#   split   12 units, 6 and 6 (k125_common.R); false groups "A" and "B"; runs renamed A_01.. B_06 as in S3
# Arms (S3's calls, s3_trx.R:220-230, 289-337, 360; the numbers are S3's):
#   1  no RUV: trxDaFit / trxDaApply with "~ 0 + group" on the filtered counts
#   3  today's design: trxRuvGrid -> trxRuvs (its controls and k, the groups as replicate sets) -> trxDaContrast
#   7  recommended: RUVSeq::RUVs on upper-quartile normalised counts, all genes as controls, the false groups as
#      replicate sets, k = E = sva::num.sv(TMM logCPM, ~ group, "be", B = 20, seed = the split's seed); W in the model
#      of the package's fit on the original filtered counts. E = 0: no RUV.
# The filter (trxExpressionFilter), TMM logCPM (trxTmmLogCpm) and the significance rule (edgeR's BH FDR < 0.05) are
# the package's, as in S3. The package kernel of arm 7 is being written on branch kit-129-trx; this is the
# simulator's own arm-7 path on the pinned tree.

local({
    file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
    source(file.path(dirname(file), "k125_common.R"))
})
args <- k125Args()
suppressMessages(pkgload::load_all(".", quiet = TRUE))
for (package in c("edgeR", "RUVSeq", "EDASeq", "sva")) requireNamespace(package)
ENTRY <- "B_vs_A=groupB-groupA"

# ---- the public counts ---------------------------------------------------------------------------------------------
inputs <- file.path(args$data_dir, c("WT_countdata.tar.gz", "exclude.lst"))
unpacked <- file.path(tempdir(), "k125-trx"); dir.create(unpacked)
utils::untar(inputs[[1L]], exdir = unpacked)
files <- sort(list.files(unpacked, pattern = "\\.gbgout$", recursive = TRUE, full.names = TRUE))
stopifnot(length(files) == 48L)
tables <- lapply(files, utils::read.delim, header = FALSE, quote = "", stringsAsFactors = FALSE, col.names = c("gene", "count"))
stopifnot(all(vapply(tables, function(t) identical(t$gene, tables[[1L]]$gene), logical(1L))))
all_counts <- vapply(tables, function(t) as.double(t$count), numeric(nrow(tables[[1L]])))
replicate_of <- sub("^(WT_rep[0-9]+)_.*$", "\\1", basename(files))
dimnames(all_counts) <- list(tables[[1L]]$gene, replicate_of)
all_counts <- all_counts[!rownames(all_counts) %in% c("no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique"), ]
stopifnot(!anyDuplicated(rownames(all_counts)), !anyNA(all_counts), all(all_counts >= 0), all(all_counts == round(all_counts)))
excluded <- sub("^(WT_rep[0-9]+)_.*$", "\\1", grep("^WT_", readLines(inputs[[2L]]), value = TRUE))
pool <- setdiff(colnames(all_counts), excluded)
stopifnot(length(excluded) == 6L, length(pool) == 42L)

refuse <- function(message) rlang::abort(message, class = "s3_refusal")
is_refusal <- function(e) inherits(e, "trx_error") || inherits(e, "s3_refusal")
# S3's daTest (s3_trx.R:220-230): the package's fit and test halves with an arm's covariates.
daTest <- function(counts, design, W = NULL) {
    pdata <- data.frame(design, row.names = design$Run)
    if (!is.null(W)) {
        W <- as.matrix(W)
        colnames(W) <- paste0("W_", seq_len(ncol(W)))
        pdata <- cbind(pdata, W)
    }
    formula <- if (is.null(W)) "~ 0 + group" else trxDaFormula(pdata)
    table <- trxDaApply(trxDaFit(counts, pdata, formula), ENTRY)$results[[1L]]$top_table[rownames(counts), ]
    list(p = table$PValue, q = table$FDR, k = if (is.null(W)) 0L else ncol(W))
}

prepare <- function(units, seed) suppressWarnings(suppressMessages({
    group <- factor(rep(c("A", "B"), each = 6L), levels = c("A", "B"))
    runs <- sprintf("%s_%02d", group, stats::ave(seq_along(group), group, FUN = seq_along))
    counts <- all_counts[, c(units$A, units$B)]
    colnames(counts) <- runs
    design <- data.frame(Run = runs, group = group, stringsAsFactors = FALSE)
    y <- trxExpressionFilter(counts, design$group)$filtered_counts
    n <- trxTmmLogCpm(y, design$group)
    model <- stats::model.matrix(~ group, data = design)
    uq <- EDASeq::betweenLaneNormalization(y, which = "upper", round = TRUE)
    E <- as.integer(sva::num.sv(n$log_cpm_tmm, model, method = "be", B = 20, seed = seed))
    genes <- rownames(y)
    quiet <- function(f) function() suppressWarnings(suppressMessages(f()))
    list(arms = lapply(list(
        `1` = function() daTest(y, design),
        `3` = function() {
            grid <- trxRuvGrid(y, n$log_cpm_tmm, design)
            s <- trxRuvs(y, n$log_cpm_tmm, design, grid$control_genes, grid$best$best_k)
            table <- trxDaContrast(s$ruvs_counts, s$pdata, ENTRY)$top_table[genes, ]
            list(p = table$PValue, q = table$FDR, k = ncol(s$W))
        },
        `7` = function() if (E == 0L) list(same_as_no_ruv = TRUE) else
            daTest(y, design, RUVSeq::RUVs(uq, genes, E, RUVSeq::makeGroups(design$group))$W)), quiet))
}))

k125Run("trx", args, pool = pool, batch = NULL, types = "random",
    designs = c(`1` = "no RUV", `3` = "today's design", `7` = "recommended"), prepare = prepare, is_refusal = is_refusal,
    inputs = inputs, packages = c("MultiScholaR", "edgeR", "limma", "RUVSeq", "EDASeq", "sva", "pkgload"),
    notes = sprintf("genes in the table: %d; excluded replicates (authors' list): %s", nrow(all_counts), paste(excluded, collapse = ",")))
