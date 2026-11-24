# T-TAB additions of the pre-release repair (audit I09, I12, I13): DIA_LIMPA and transcriptomics
# with RUV off on their generated inputs, and the values of matrix_pre_ruv.tsv against the product's
# own record of the matrix that entered RUV on DIA, LFQ, TMT and DIA_LIMPA.

tabRepairRun <- function(route, arm, on_complete = NULL) {
    root <- driverPkgRoot()
    skip_if(is.na(root), "the package checkout is not known (run through val-tests.R)")
    corpus <- valWriteCorpus(valCorpusFor(route), withr::local_tempdir(.local_envir = parent.frame()),
                             root)
    cell <- valCorpusCell(route, corpus)
    script <- valActionScript(route, assays = cell$assays,
                              import_tool = if (cell$omic == "metabolomics") "msdial",
                              contrasts = cell$contrasts, ruv = valArmRuvRun(arm))
    out <- withr::local_tempdir(.local_envir = parent.frame())
    res <- valRunWorkflow(script, files = cell$files, design = cell$design, arm = arm,
                          helpers = valGuiHelpers(root), uniprot_cache = cell$uniprot_cache,
                          annotation = cell$annotation, forced_inputs = cell$forced_inputs,
                          tables_dir = out, pkg_root = root,
                          work_dir = withr::local_tempdir(.local_envir = parent.frame()),
                          on_complete = on_complete)
    list(res = res, out = out)
}

for (route in c("prot_dia_limpa", "trx_screen")) {
    local({
        route <- route
        test_that(paste("T-TAB-OFF:", route, "with RUV off completes and writes no RUV factors"), {
            x <- tabRepairRun(route, valArmRuvOff())
            expect_identical(x$res$exit, "completed",
                             info = paste(x$res$steps$step, x$res$steps$outcome,
                                          x$res$steps$message, collapse = "; "))
            expect_true(file.exists(file.path(x$out, "da.tsv")))
            expect_false(file.exists(file.path(x$out, "ruv.tsv")))
            if (route == "prot_dia_limpa") {
                s <- x$res$steps
                expect_identical(s$head[s$step == "ruv_covariates"], "correlation_filtered")
            }
        })
    })
}

for (route in c("prot_dia", "prot_lfq", "prot_tmt", "prot_dia_limpa")) {
    local({
        route <- route
        test_that(paste("T-PRE-RUV:", route, "matrix_pre_ruv.tsv is the product's record of the RUV input"), {
            direct <- NULL
            read_direct <- function(ctx) {
                ns <- asNamespace("MultiScholaR")
                m <- ctx$manager
                state <- if (route == "prot_dia_limpa") "ruv_covariates" else "ruv_corrected"
                direct <<- ns$withWorkflowStateScope(m, function() {
                    rec <- ns$readComponentMetadata(m, ns$getDIAStateReference(m, state))
                    inp <- ns$protIntegrationInputRuvInput(m, rec)
                    tbl <- ns$readDIAExecutionComponents(
                        m, inp$reference, fields = c("protein_quant_table", "design_matrix")
                    )@protein_quant_table
                    msrval:::valFeatureTable(tbl, rec$data$slots$protein_id_column, "all")
                })
                invisible(NULL)
            }
            x <- tabRepairRun(route, valArmShipped(), on_complete = read_direct)
            expect_identical(x$res$exit, "completed")
            expect_false(is.null(direct))
            want <- tempfile(fileext = ".tsv")
            utils::write.table(direct, want, sep = "\t", quote = FALSE, row.names = FALSE,
                               na = "NA")
            expect_identical(sha256File(file.path(x$out, "matrix_pre_ruv.tsv")), sha256File(want))
        })
    })
}
