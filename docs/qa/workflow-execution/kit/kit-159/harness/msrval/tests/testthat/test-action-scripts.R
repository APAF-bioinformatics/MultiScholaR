# The accepted per-route action scripts (KIT-159 7.4; VAL-41, VAL-26): valActionScript, pinned
# against action-scripts.md and against the product's own QC spec tables.

scriptIds <- function(route, ...) {
    vapply(valActionScript(route, ...)$steps, function(s) s$id, character(1))
}

test_that("every route has a script whose steps carry an id, a declared state and a run function", {
    for (route in driverRoutes()) {
        s <- valActionScript(route)
        expect_identical(s$route, route)
        expect_true(s$omic %in% c("proteomics", "lipidomics", "metabolomics", "transcriptomics"))
        for (step in s$steps) {
            expect_true(is.character(step$id) && nzchar(step$id), info = route)
            expect_true(is.function(step$run), info = paste(route, step$id))
            expect_true(identical(step$declared, "diagnostic") ||
                            (is.list(step$declared) &&
                                 names(step$declared) %in% c("head", "tab")),
                        info = paste(route, step$id))
        }
        expect_true(valCheckActionScript(s$steps))
    }
    expect_error(valActionScript("prot_unknown"), class = "val_spec_invalid")
})

test_that("no script presses revert, skip, export, load, report, enrichment or Cancel (7.4, A17)", {
    banned <- "revert|skip|bypass|export|load_filtered|generate_report|enrich|copy_to|cancel"
    for (route in driverRoutes()) {
        presses <- unlist(lapply(valActionScript(route)$steps, function(s) s$press))
        expect_false(any(grepl(banned, presses, ignore.case = TRUE)), info = route)
    }
    with_summary <- valActionScript("prot_dia", summary = list(label = "a", description = "b"))
    presses <- unlist(lapply(with_summary$steps, function(s) s$press))
    expect_true("session_summary-save_workflow_args" %in% presses)
    expect_false(any(grepl("generate_report|copy_to|export_session", presses)))
})

test_that("the common skeleton per route: import, design, QC, normalisation, correlation, DA", {
    expect_identical(scriptIds("trx_screen"),
                     c("import", "design", "qc_filter", "qc_finalize", "norm", "da"))
    for (route in c("lipid_canonical", "metab_lc", "metab_gc", "metab_combined")) {
        expect_identical(scriptIds(route),
                         c("import", "design", "qc_intensity", "qc_duplicates", "qc_itsd",
                           "qc_finalize", "norm", "correlation", "da"), info = route)
    }
    limpa <- scriptIds("prot_dia_limpa")
    expect_identical(tail(limpa, 4L), c("norm", "correlation", "ruv_covariates", "da"))
    for (route in c("prot_dia", "prot_lfq", "prot_lfq_fragpipe", "prot_tmt")) {
        expect_identical(tail(scriptIds(route), 3L), c("norm", "correlation", "da"), info = route)
        expect_identical(head(scriptIds(route), 2L), c("import", "design"), info = route)
    }
    expect_identical(tail(scriptIds("prot_dia", summary = list(label = "a", description = "b")), 1L),
                     "summary_save")
})

test_that("the QC sub-tabs are the product's spec tables, in order (DIA, DIA_LIMPA, LFQ, TMT)", {
    ns <- asNamespace("MultiScholaR")
    qcIds <- function(route) {
        peptide <- if (route %in% c("DIA", "DIA_LIMPA")) {
            vapply(ns$getProtQcPeptideModuleSpecs(route), function(x) paste0("qc_peptide_", x$id),
                   character(1))
        }
        protein <- vapply(ns$getProtQcProteinModuleSpecs(route),
                          function(x) paste0("qc_protein_", x$id), character(1))
        c(peptide, protein)
    }
    strip <- function(ids) ids[startsWith(ids, "qc_")]
    expect_identical(strip(scriptIds("prot_dia")), qcIds("DIA"))
    expect_identical(strip(scriptIds("prot_dia_limpa")), qcIds("DIA_LIMPA"))
    expect_identical(strip(scriptIds("prot_lfq")), qcIds("LFQ"))
    expect_identical(strip(scriptIds("prot_lfq_fragpipe")), qcIds("LFQ"))
    expect_identical(strip(scriptIds("prot_tmt")), qcIds("TMT"))
    expect_length(qcIds("DIA"), 12L)
    expect_length(qcIds("DIA_LIMPA"), 7L)
    expect_length(qcIds("LFQ"), 4L)
})

test_that("each QC step declares the state its sub-tab commits; DA declares its tab complete", {
    s <- valActionScript("prot_dia")
    byId <- stats::setNames(s$steps, vapply(s$steps, function(x) x$id, character(1)))
    expect_identical(byId$qc_peptide_qvalue_filter$declared, list(head = "qvalue_filtered"))
    expect_identical(byId$qc_protein_replicate_filter$declared,
                     list(head = "protein_replicate_filtered"))
    expect_identical(byId$correlation$declared, list(head = "correlation_filtered"))
    expect_identical(byId$da$declared, list(tab = c(differential_expression = "complete")))
    expect_identical(valActionScript("prot_lfq")$steps[[2L]]$declared,
                     list(head = "protein_s4_initial"))
    lipid <- valActionScript("lipid_canonical")$steps
    expect_identical(lipid[[5L]]$declared, "diagnostic")
    expect_identical(lipid[[9L]]$declared, list(tab = c(differential_analysis = "complete")))
})

test_that("action-scripts.md lists exactly each route's step ids (the accepted scripts, VAL-41)", {
    md <- file.path(Sys.getenv("VAL_HARNESS_DIR"), "..", "action-scripts.md")
    expect_true(file.exists(md))
    lines <- grep("^- `[a-z_]+`: ", readLines(md), value = TRUE)
    pinned <- stats::setNames(
        lapply(sub("^- `[a-z_]+`: ", "", lines), function(x) trimws(strsplit(x, ",")[[1L]])),
        sub("^- `([a-z_]+)`: .*$", "\\1", lines))
    expect_setequal(names(pinned), driverRoutes())
    for (route in driverRoutes()) {
        expect_identical(scriptIds(route), pinned[[route]], info = route)
    }
})
