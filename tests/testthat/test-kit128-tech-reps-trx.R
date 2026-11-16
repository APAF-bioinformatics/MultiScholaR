# KIT-128 (decisions D1, D2): declared technical replicates of a transcriptomics design are summed per biological unit at
# the design step (`edgeR::sumTechReps`; edgeR User's Guide 4.9.1, limma User's Guide 18.2.5), by the proteomics rule
# (resolveProtTechnicalReplicateBlock); the route then runs on the units. The fixture is generated here: 2 groups, 4
# biological units per group, 2 libraries per unit, no true effect (KIT-124's C6 recipe: a Gamma multiplier per unit,
# biological dispersion 0.09; negative binomial libraries, technical dispersion 0.01). Size of the null-rate test: 20
# seeds by 500 genes (those the expression filter keeps); the limit is KIT-124's rule 1, the one-sided 97.5% exact
# binomial upper bound of the share of P < 0.05 at most 0.075.
library(testthat)

kit128TrxDesign <- function(replicates = c("technical", "unique")) {
  replicates <- match.arg(replicates)
  groups <- rep(c("A", "B"), each = 8L)
  unit <- rep(rep(seq_len(4L), each = 2L), times = 2L)
  runs <- sprintf("%s%d_%d", groups, unit, rep(seq_len(2L), times = 8L))
  data.frame(Run = runs, group = groups, replicates = unit,
    tech_rep_group = if (identical(replicates, "technical")) paste(groups, unit, sep = "_") else runs, stringsAsFactors = FALSE)
}

kit128TrxCounts <- function(seed, genes = 500L) {
  design <- kit128TrxDesign()
  unit <- match(design$tech_rep_group, unique(design$tech_rep_group))
  withr::with_seed(seed, {
    lambda <- 2^stats::rnorm(genes, 8, 1.5)
    depth <- stats::runif(16L, 0.7, 1.3)
    biological <- matrix(stats::rgamma(genes * 8L, shape = 1 / 0.09, scale = 0.09), nrow = genes)
    mu <- lambda * biological[, unit] * rep(depth, each = genes)
    matrix(stats::rnbinom(genes * 16L, mu = mu, size = 1 / 0.01), nrow = genes,
      dimnames = list(sprintf("G%03d", seq_len(genes)), design$Run))
  })
}

#' The design state's object, as the design step builds it from the builder's result.
kit128TrxObject <- function(counts, design) {
  table <- cbind(data.frame(GeneID = rownames(counts), GeneName = paste0("N", rownames(counts)), stringsAsFactors = FALSE),
    as.data.frame(counts, check.names = FALSE))
  rownames(table) <- NULL
  trxDesignObject(list(data_cln = list(counts = table), design_matrix = design, group_levels = c("A", "B")),
    list(lipid_id_col = "GeneID", annotation_col = "GeneName"))
}

#' The route's kernels in order on an object: the expression filter, TMM, the RUV grid, RUVs, the edgeR QL test.
kit128TrxRouteP <- function(object) {
  suppressWarnings(suppressMessages({
    x <- trxAssayInputs(object)
    y <- trxExpressionFilter(x$counts, x$group)$filtered_counts
    n <- trxTmmLogCpm(y, x$group)
    grid <- trxRuvGrid(y, n$log_cpm_tmm, x$design)
    s <- trxRuvs(y, n$log_cpm_tmm, x$design, grid$control_genes, grid$best$best_k)
    trxDaContrast(s$ruvs_counts, s$pdata, "B_vs_A = groupB-groupA")$top_table$PValue
  }))
}

test_that("KIT-128: the design object of a duplicated design is edgeR::sumTechReps of the counts with one row per unit", {
  skip_if_not_installed("edgeR")
  counts <- kit128TrxCounts(1L)
  design <- kit128TrxDesign()
  object <- kit128TrxObject(counts, design)
  x <- trxAssayInputs(object)
  units <- unique(design$tech_rep_group)
  expected <- edgeR::sumTechReps(counts, ID = design$tech_rep_group)
  expect_identical(colnames(x$counts), units)
  expect_equal(x$counts, expected[, units], tolerance = 0)
  expect_identical(as.character(x$design$Run), units)
  expect_identical(as.character(x$design$group), rep(c("A", "B"), each = 4L))
  expect_identical(x$design$replicates, rep(seq_len(4L), times = 2L))
  expect_identical(x$design$tech_rep_group, units)
  # The summed runs are recorded, unit to runs; no technical-replicate set is left.
  expect_identical(object@args$technical_replicates_summed, split(design$Run, factor(design$tech_rep_group, levels = units)))
  expect_false(anyDuplicated(object@design_matrix$tech_rep_group) > 0L)
})

test_that("KIT-128: a design without technical replicates gives the object it gave before", {
  counts <- kit128TrxCounts(1L)
  design <- kit128TrxDesign("unique")
  object <- kit128TrxObject(counts, design)
  table <- cbind(data.frame(GeneID = rownames(counts), GeneName = paste0("N", rownames(counts)), stringsAsFactors = FALSE),
    as.data.frame(counts, check.names = FALSE))
  rownames(table) <- NULL
  expect_identical(object, createTranscriptomicsAssayData(list(counts = table), design, c("A", "B"),
    roles = list(gene_id_column = "GeneID", annotation_id_column = "GeneName", technical_replicate_id = "tech_rep_group")))
  expect_null(object@args$technical_replicates_summed)
  expect_identical(colnames(trxAssayInputs(object)$counts), design$Run)
})

test_that("KIT-128: null rate at 0.05 on duplicated libraries through the route's kernels is within the limit (upper bound at most 0.075)", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("RUVSeq")
  p <- unlist(lapply(seq_len(20L), function(seed) kit128TrxRouteP(kit128TrxObject(kit128TrxCounts(seed), kit128TrxDesign()))))
  expect_gt(length(p), 9900L)
  bound <- stats::binom.test(sum(p < 0.05), length(p), alternative = "less", conf.level = 0.975)$conf.int[[2L]]
  expect_lte(bound, 0.075, label = sprintf("null rate %.4f, its 97.5%% upper bound %.4f", mean(p < 0.05), bound))
})

test_that("KIT-128: the design check counts a group's biological units, not its libraries", {
  design <- kit128TrxDesign()
  # Group A is one unit sequenced twice; group B keeps two units of two libraries.
  design <- design[c(1L, 2L, 9L, 10L, 11L, 12L), ]
  check <- trxDesignCheck(list(design_matrix = design, group_levels = c("A", "B")), list(sample_columns = design$Run))
  expect_false(check$valid)
  expect_match(check$errors, "trx_design_group_singleton", fixed = TRUE)
  # Two units per group pass, as two libraries per group without technical replicates do.
  design <- kit128TrxDesign()[c(1L, 2L, 3L, 4L, 9L, 10L, 11L, 12L), ]
  expect_true(trxDesignCheck(list(design_matrix = design, group_levels = c("A", "B")), list(sample_columns = design$Run))$valid)
  unique_runs <- kit128TrxDesign("unique")[c(1L, 2L, 9L, 10L), ]
  expect_true(trxDesignCheck(list(design_matrix = unique_runs, group_levels = c("A", "B")), list(sample_columns = unique_runs$Run))$valid)
})
