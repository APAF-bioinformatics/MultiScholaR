# The standard tables of a workflow cell (KIT-159 7.3: valReadArmTables; the stage rule A12,
# revision 1.8 item 4), on hand-made DA files and state objects shaped as the product writes them
# (proteomics `da_proteins_<contrast>_long_annot.tsv`, R/func_prot_da_results_io.R; lipidomics
# `de_<mode>_lipids_<contrast>_long_annot.tsv`, R/func_lipid_da_export.R; transcriptomics
# `de_results_<contrast>.csv`, R/utils_workflow_trx_declarations.R:903).

writeTsv <- function(x, path) {
    utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

readTsv <- function(path) {
    utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

protCell <- function() {
    dir <- withr::local_tempdir(.local_envir = parent.frame())
    da_dir <- file.path(dir, "da")
    dir.create(da_dir)
    writeTsv(data.frame(Protein.Group = c("P1", "P2", "P3"), gene_name = c("A", "B", "C"),
                        fdr_qvalue = c(0.01, 0.2, 0.04), raw_pvalue = c(0.001, 0.1, 0.01),
                        log2FC = c(1.5, -0.2, -2), comparison = "groupKO-groupWT",
                        log2norm.s1.KO = c(10, 11, 12), log2norm.s2.WT = c(12, 11, 10),
                        check.names = FALSE),
             file.path(da_dir, "da_proteins_groupKO-groupWT_long_annot.tsv"))
    quant <- function(shift) {
        data.frame(Protein.Group = c("P1", "P2", "P3", "P4"), s1 = c(10, 11, 12, 13) + shift,
                   s2 = c(12, 11, 10, 9) + shift)
    }
    states <- list(
        normalized = list(protein_quant_table = quant(0), protein_id_column = "Protein.Group"),
        correlation_filtered = list(protein_quant_table = quant(0.5)[1:3, ],
                                    protein_id_column = "Protein.Group"))
    list(paths = list(da_output_dir = da_dir), out = file.path(dir, "tables"),
         reader = function(name) states[[name]])
}

test_that("proteomics: da.tsv's columns and values, the matrices by stage, features", {
    cell <- protCell()
    paths <- valReadArmTables(NULL, cell$paths, "prot_dia", cell$out, ruv = "off",
                              state_reader = cell$reader, model_state = "correlation_filtered")
    expect_setequal(names(paths), c("da", "matrix_pre_ruv", "matrix", "features"))
    da <- readTsv(paths[["da"]])
    expect_identical(names(da), DA_COLUMNS)
    expect_identical(da$feature_id, c("P1", "P2", "P3"))
    expect_identical(unique(da$assay), "all")
    expect_identical(unique(da$contrast), "groupKO-groupWT")
    expect_identical(da$log2fc, c(1.5, -0.2, -2))
    expect_identical(da$p, c(0.001, 0.1, 0.01))
    expect_identical(da$s_own, c(0.01, 0.2, 0.04))
    expect_identical(unique(da$call_statistic), "qvalue")
    expect_identical(da$called_own, c(TRUE, FALSE, TRUE))
    # Whole numbers read back from a TSV are integers: compared by value.
    expect_equal(da$mean_expression, c(11, 11, 11))
    pre <- readTsv(paths[["matrix_pre_ruv"]])
    expect_identical(names(pre), c("feature_id", "assay", "s1", "s2"))
    expect_identical(nrow(pre), 4L)
    model <- readTsv(paths[["matrix"]])
    expect_identical(model$s1, c(10.5, 11.5, 12.5))
    feats <- readTsv(paths[["features"]])
    expect_identical(names(feats), c("feature_id", "assay", "in_matrix_pre_ruv", "in_matrix",
                                     "tested"))
    expect_identical(feats$tested, c(TRUE, TRUE, TRUE, FALSE))
    expect_false(file.exists(file.path(cell$out, "matrix_display.tsv")))
})

test_that("a fold-change threshold enters called_own; a missing DA file is refused", {
    cell <- protCell()
    paths <- valReadArmTables(NULL, cell$paths, "prot_dia", cell$out, ruv = "off",
                              lfc_cutoff = 1.8, state_reader = cell$reader,
                              model_state = "correlation_filtered")
    expect_identical(readTsv(paths[["da"]])$called_own, c(FALSE, FALSE, TRUE))
    unlink(list.files(cell$paths$da_output_dir, full.names = TRUE))
    expect_error(valReadArmTables(NULL, cell$paths, "prot_dia", cell$out, ruv = "off",
                                  state_reader = cell$reader,
                                  model_state = "correlation_filtered"),
                 class = "val_spec_invalid")
})

lipidCell <- function() {
    dir <- withr::local_tempdir(.local_envir = parent.frame())
    da_dir <- file.path(dir, "da")
    dir.create(da_dir)
    for (mode in c("posmode", "negmode")) {
        assay <- if (mode == "posmode") "LCMS_Pos" else "LCMS_Neg"
        writeTsv(data.frame(lipid_id = paste0(mode, 1:2), assay = assay, logFC = c(1, -1),
                            raw_pvalue = c(0.01, 0.5), fdr_qvalue = c(0.02, 0.5),
                            comparison = "groupKO-groupWT", intensity.s1 = c(5, 6),
                            intensity.s2 = c(7, 8), ITSD = NA, check.names = FALSE),
                 file.path(da_dir, sprintf("de_%s_lipids_groupKO-groupWT_long_annot.tsv", mode)))
    }
    assay <- function(mode, shift) {
        data.frame(lipid_id = c(paste0(mode, 1:2), "ITSD_1"), s1 = c(5, 6, 1) + shift,
                   s2 = c(7, 8, 1) + shift)
    }
    data <- function(shift) list(LCMS_Pos = assay("posmode", shift), LCMS_Neg = assay("negmode", shift))
    w <- matrix(c(0.1, -0.1), ncol = 1L, dimnames = list(c("s1", "s2"), "W1"))
    states <- list(
        lipid_normalized = list(lipid_data = data(0), lipid_id_column = "lipid_id",
                                internal_standard_regex = "^ITSD"),
        lipid_correlation_filtered = list(
            lipid_data = data(0.25), lipid_id_column = "lipid_id",
            internal_standard_regex = "^ITSD",
            args = list(ruv_uncorrected_assays = data(0),
                        ruv_factors = list(LCMS_Pos = list(W = w, k = 1L, design = "groups_all_features",
                                                           controls_n = 1L, features_tested = 2L,
                                                           fallback = "none"),
                                           LCMS_Neg = list(W = w, k = 1L,
                                                           design = "groups_all_features",
                                                           controls_n = 1L, features_tested = 2L,
                                                           fallback = "none")))))
    list(paths = list(da_output_dir = da_dir), out = file.path(dir, "tables"),
         reader = function(name) states[[name]])
}

test_that("lipidomics with factors in the model: matrix.tsv uncorrected, W, a labelled display", {
    cell <- lipidCell()
    paths <- valReadArmTables(NULL, cell$paths, "lipid_canonical", cell$out,
                              ruv = "factors_in_model", state_reader = cell$reader,
                              model_state = "lipid_correlation_filtered", itsd_applied = TRUE)
    expect_setequal(names(paths), c("da", "matrix_pre_ruv", "matrix", "features", "ruv", "itsd",
                                    "matrix_display"))
    da <- readTsv(paths[["da"]])
    expect_identical(names(da), DA_COLUMNS)
    expect_setequal(unique(da$assay), c("LCMS_Pos", "LCMS_Neg"))
    model <- readTsv(paths[["matrix"]])
    expect_equal(model$s1[model$feature_id == "posmode1"], 5)
    display <- readTsv(paths[["matrix_display"]])
    expect_identical(names(display)[1:3], c("feature_id", "assay", "label"))
    expect_identical(unique(display$label), "display; not what is tested")
    expect_identical(display$s1[display$feature_id == "posmode1"], 5.25)
    ruv <- readTsv(paths[["ruv"]])
    expect_identical(names(ruv), c("assay", "design", "k", "controls_n", "controls_share",
                                   "features_tested", "fallback", "sample", "factor", "w"))
    expect_identical(nrow(ruv), 4L)
    itsd <- readTsv(paths[["itsd"]])
    expect_identical(names(itsd), c("assay", "applied", "standards_n"))
    expect_identical(itsd$standards_n, c(1L, 1L))
})

test_that("lipidomics with RUV correcting the matrix: matrix.tsv is the corrected one, no display", {
    cell <- lipidCell()
    paths <- valReadArmTables(NULL, cell$paths, "lipid_canonical", cell$out, ruv = "correction",
                              state_reader = cell$reader,
                              model_state = "lipid_correlation_filtered", itsd_applied = TRUE)
    expect_false("matrix_display" %in% names(paths))
    model <- readTsv(paths[["matrix"]])
    expect_identical(model$s1[model$feature_id == "posmode1"], 5.25)
})

test_that("transcriptomics: da.tsv from the edgeR table, log-CPM matrices", {
    dir <- withr::local_tempdir()
    da_dir <- file.path(dir, "da")
    dir.create(da_dir)
    utils::write.csv(data.frame(GeneID = c("g1", "g2"), GeneName = c("A", "B"), logFC = c(2, 0),
                                logCPM = c(7.5, 3), PValue = c(1e-4, 0.8), FDR = c(2e-4, 0.8)),
                     file.path(da_dir, "de_results_Sera_vs_RPMI.csv"), row.names = FALSE)
    counts <- data.frame(GeneID = c("g1", "g2", "g3"), GeneName = c("A", "B", "C"),
                         s1 = c(100, 20, 5), s2 = c(120, 25, 0), s3 = c(400, 22, 3),
                         s4 = c(380, 18, 2))
    obj <- list(count_data = list(counts = counts), gene_id_column = "GeneID",
                annotation_id_column = "GeneName")
    paths <- valReadArmTables(NULL, list(da_output_dir = da_dir), "trx_screen",
                              file.path(dir, "tables"), ruv = "off",
                              state_reader = function(name) obj, model_state = "trx_qc_complete")
    da <- readTsv(paths[["da"]])
    expect_identical(names(da), DA_COLUMNS)
    expect_identical(unique(da$contrast), "Sera_vs_RPMI")
    expect_identical(da$mean_expression, c(7.5, 3))
    expect_identical(da$s_own, c(2e-4, 0.8))
    pre <- readTsv(paths[["matrix_pre_ruv"]])
    expect_identical(names(pre), c("feature_id", "assay", "s1", "s2", "s3", "s4"))
    expect_true(all(is.finite(as.matrix(pre[, -(1:2)]))))
})
