# T-TAB with RUV applied on a classic proteomics route (KIT-159 7.3; S0 run 45277101451): with the
# shipped arm (RUV-III automatic) LFQ commits `ruv_corrected` and no `normalized` state, so the
# matrix that entered RUV is read from the RUV state's own product. T-TAB's fixture runs use RUV off
# on this route and never reached this path.

test_that("T-TAB-RUV: prot_lfq on its generated input with the shipped arm writes every standard table", {
    root <- driverPkgRoot()
    corpus <- valWriteCorpus("prot_lfq", withr::local_tempdir(), root)
    cell <- valCorpusCell("prot_lfq", corpus)
    out <- withr::local_tempdir()
    res <- valRunWorkflow(cell$script, files = cell$files, design = cell$design,
                          arm = valArmShipped(), helpers = valGuiHelpers(root),
                          uniprot_cache = cell$uniprot_cache, annotation = cell$annotation,
                          tables_dir = out, pkg_root = root, work_dir = withr::local_tempdir())
    expect_identical(res$exit, "completed",
                     info = paste(res$steps$step, res$steps$outcome, res$steps$message,
                                  collapse = "; "))
    expect_identical(res$steps$head[res$steps$step == "norm"], "ruv_corrected")
    pre <- utils::read.delim(file.path(out, "matrix_pre_ruv.tsv"), check.names = FALSE,
                             stringsAsFactors = FALSE)
    model <- utils::read.delim(file.path(out, "matrix.tsv"), check.names = FALSE,
                               stringsAsFactors = FALSE)
    expect_identical(names(pre)[1:2], c("feature_id", "assay"))
    expect_setequal(names(pre)[-(1:2)], names(model)[-(1:2)])
    # The matrix that entered RUV is not the corrected one.
    shared <- intersect(pre$feature_id, model$feature_id)
    expect_gt(length(shared), 0L)
    s <- names(model)[[3L]]
    expect_false(isTRUE(all.equal(pre[[s]][match(shared, pre$feature_id)],
                                  model[[s]][match(shared, model$feature_id)])))
    expect_true(file.exists(file.path(out, "ruv.tsv")))
    da <- utils::read.delim(file.path(out, "da.tsv"), stringsAsFactors = FALSE)
    expect_gt(nrow(da), 0L)
})
