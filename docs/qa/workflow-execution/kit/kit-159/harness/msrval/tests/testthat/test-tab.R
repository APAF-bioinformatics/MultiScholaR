# T-TAB (KIT-159 section 10.3): per route, the standard tables' columns on the route's committed
# e2e fixture (tests/testdata/e2e), through the driver core end to end: the helpers extracted from
# this checkout (valGuiHelpers), the accepted script (valActionScript), valRunWorkflow under
# testServer, and valReadArmTables. The six-sample fixtures take the GUI-entry lanes' permissive QC
# settings as declared forced inputs and RUV off (limpa: the shipped covariates), as those lanes do;
# these runs pin columns, not science.

for (route in driverRoutes()) {
    local({
        route <- route
        test_that(paste("T-TAB:", route, "writes the standard tables with their columns"), {
            root <- driverPkgRoot()
            cell <- valFixtureCell(route, root)
            helpers <- valGuiHelpers(root)
            # The covariate routes with the shipped arm, so ruv.tsv and matrix_display.tsv are written
            # and checked on both; the others with RUV off as their GUI-entry lanes run (10.3, r1.11).
            arm <- if (route %in% c("prot_dia_limpa", "trx_screen")) valArmShipped() else valArmRuvOff()
            out <- withr::local_tempdir()
            res <- valRunWorkflow(cell$script, files = cell$files, design = cell$design, arm = arm,
                                  helpers = helpers, forced_inputs = cell$forced_inputs,
                                  uniprot_cache = cell$uniprot_cache, tables_dir = out,
                                  pkg_root = root)
            expect_identical(res$exit, "completed",
                             info = paste(res$steps$step, res$steps$outcome, res$steps$reason,
                                          res$steps$class, res$steps$message, collapse = "; "))
            # A refused diagnostic step is recorded and the cell continues (7.3, A21): on the
            # fixtures the ITSD analysis refuses, as the GUI-entry oracles show.
            scored <- res$steps$step %in% names(Filter(function(s) !identical(s$declared, "diagnostic"),
                stats::setNames(cell$script$steps, vapply(cell$script$steps, `[[`, "", "id"))))
            expect_true(all(res$steps$outcome[scored] == "ok"),
                        info = paste(res$steps$step, res$steps$outcome, collapse = "; "))
            da <- utils::read.delim(file.path(out, "da.tsv"), stringsAsFactors = FALSE)
            expect_identical(names(da), DA_COLUMNS)
            expect_gt(nrow(da), 0L)
            expect_true(all(nzchar(da$contrast)))
            if (cell$omic %in% c("lipidomics", "metabolomics")) {
                expect_setequal(unique(da$assay), cell$assays)
            } else {
                expect_identical(unique(da$assay), "all")
            }
            for (m in c("matrix_pre_ruv.tsv", "matrix.tsv")) {
                x <- utils::read.delim(file.path(out, m), stringsAsFactors = FALSE,
                                       check.names = FALSE)
                expect_identical(names(x)[1:2], c("feature_id", "assay"), info = m)
                expect_gt(ncol(x), 2L)
            }
            feats <- utils::read.delim(file.path(out, "features.tsv"), stringsAsFactors = FALSE)
            expect_identical(names(feats), c("feature_id", "assay", "in_matrix_pre_ruv", "in_matrix",
                                             "tested"))
            if (route %in% c("prot_dia_limpa", "trx_screen")) {
                expect_true(file.exists(file.path(out, "ruv.tsv")))
                expect_true(file.exists(file.path(out, "matrix_display.tsv")))
                display <- utils::read.delim(file.path(out, "matrix_display.tsv"),
                                             stringsAsFactors = FALSE, check.names = FALSE)
                expect_identical(names(display)[1:3], c("feature_id", "assay", "label"))
                expect_gt(ncol(display), 3L)
            }
            if (cell$omic %in% c("lipidomics", "metabolomics")) {
                expect_true(file.exists(file.path(out, "itsd.tsv")))
            }
            expect_true(file.exists(res$ledger_path))
        })
    })
}
