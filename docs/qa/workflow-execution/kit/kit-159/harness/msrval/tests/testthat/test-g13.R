# G13 under the SWD-083 tolerance (VAL-133 (1); KIT-199 implementation audit I05, O2, O4): the
# measured DIA-with-limpa case, OpenBLAS against the reference BLAS (KIT-199 R4, rows of
# `r4/da-runner.tsv` and `r1/da-runner.tsv`, `g1a/SOURCES.tsv`), and the Linux bundle's GUI-entry
# lane's six DIA-with-limpa DA steps excluded by name from any installer acceptance.

g13Table <- function(name) {
    t <- valReadTextTable(testthat::test_path("g1a", name))
    for (col in setdiff(names(t), c("Protein.Ids", "comparison"))) t[[col]] <- as.numeric(t[[col]])
    t
}

test_that("G13-SWD083: the measured DIA-with-limpa difference passes under the rule", {
    ref <- g13Table("k199-da-r1.tsv")
    bundle <- g13Table("k199-da-r4.tsv")
    # the measured case: two values cross the 10th significant digit (P00025 fdr_qvalue,
    # P00146 log2FC); largest relative difference 5.2e-11 (P00189 log2FC); every call identical
    f_ref <- valFingerprint(ref)
    expect_false(identical(valFingerprint(bundle)$sha256_10, f_ref$sha256_10))
    l1 <- list(digest = valMaskedDigest(ref), fingerprint = f_ref, table = ref,
               keys = c("Protein.Ids", "comparison"), q_col = "fdr_qvalue", fc_col = "log2FC")
    res <- valCompareParityDa(bundle, l1)
    expect_identical(res$status, "within_tolerance")
    expect_identical(res$comparison$calls_changed, 0L)
    expect_lt(max(res$comparison$columns$max_rel), 1e-10)
    called <- bundle
    called$fdr_qvalue[called$Protein.Ids == "P00038"] <- 0.06
    expect_identical(valCompareParityDa(called, l1)$status, "differs")
})

test_that("G13-O4: the bundle lane's six DIA-with-limpa DA steps are excluded by name", {
    ex <- msrval:::PARITY_INSTALLER_EXCLUDED
    expect_identical(unique(ex$lane), "prot_dia_limpa")
    expect_identical(ex$step, c("da", "summary_refusal", "summary_save", "summary_copy",
                                "summary_report", "summary_export"))
    rows <- data.frame(lane = c(rep("prot_dia_limpa", 3L), "prot_dia"),
                       step = c("norm", "da", "summary_save", "da"),
                       status = c("pass", "fail", "fail", "fail"))
    acc <- valInstallerAcceptanceRows(rows)
    expect_identical(acc$kept$step, c("norm", "da"))
    expect_identical(acc$kept$lane, c("prot_dia_limpa", "prot_dia"))
    expect_identical(acc$excluded$step, c("da", "summary_save"))
    expect_match(acc$excluded$reason[[1L]], "VAL-133")
})
