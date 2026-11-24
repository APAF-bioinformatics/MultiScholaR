# The cell-to-verdict path end to end (audit I01 to I08, I12, I14; register section 7 item 2): on
# each omic family's generated corpus, a ticket written to the template's contracts runs through
# valRunPhase (a blind lifecycle request on one unit, two counted requests of two units each) and
# valScorePhase, with the harness's own regime table and a ticket coverage regime. The inputs are
# shaped so that the audited faults would show: featureKey(feature_id, arm) returning a status
# (I01); a universe from the input role with each feature's assay, metabolomics compounds sharing
# their ids across two assays (I02); a ticket regime licensing the interval over n - k units after a
# refused unit (I03, I04, O7 (b)); within-run rows beside the run-level verdicts (I06, O8); the
# measured coverage and family bound in results.md (I05, O1); the run folder and provenance of 5.6
# and 5.7 (I07); the licence guard and the shell check on what production wrote (I08); the hook
# writing only under ticket/ for completed workflow cells (I14).

e2eTicket <- function(dir, route, pkg_root) {
    code <- file.path(dir, "ticket.R")
    writeLines(c(
        sprintf("ROUTE <- %s", deparse(route)),
        sprintf("PKG_ROOT <- %s", deparse(pkg_root)),
        sprintf("CONTRAST_FILE <- %s", deparse(file.path(dir, "contrast.txt"))),
        "UNITS <- c('u01', 'u02', 'u03', 'u04')",
        "LOOKUP <- NULL",
        "KEYS <- NULL",
        "e2eFeatures <- function(cell) {",
        "    if (cell$omic == 'proteomics') {",
        "        rep <- utils::read.delim(cell$files$search_results, stringsAsFactors = FALSE,",
        "                                 check.names = FALSE)",
        "        col <- intersect(c('Protein.Group', 'Protein.Ids'), names(rep))[[1L]]",
        "        return(data.frame(assay = 'all', feature_id = unique(rep[[col]]),",
        "                          stringsAsFactors = FALSE))",
        "    }",
        "    if (cell$omic == 'transcriptomics') {",
        "        x <- utils::read.delim(cell$files$counts, stringsAsFactors = FALSE)",
        "        LOOKUP <<- x[, c('GeneID', 'GeneName')]",
        "        return(data.frame(assay = 'all', feature_id = x$GeneID, stringsAsFactors = FALSE))",
        "    }",
        "    do.call(rbind, lapply(seq_along(cell$assays), function(i) {",
        "        x <- utils::read.delim(cell$files[[i]], check.names = FALSE,",
        "                               stringsAsFactors = FALSE)",
        "        id <- intersect(c('Alignment ID', 'LipidName'), names(x))[[1L]]",
        "        data.frame(assay = cell$assays[[i]], feature_id = unique(as.character(x[[id]])),",
        "                   stringsAsFactors = FALSE)",
        "    }))",
        "}",
        "prepareInputs <- function(data_dir, work_dir) {",
        "    corpus <- msrval::valWriteCorpus(msrval::valCorpusFor(ROUTE),",
        "                                     file.path(work_dir, 'corpus'), PKG_ROOT)",
        "    cell <- msrval::valCorpusCell(ROUTE, corpus)",
        "    feats <- e2eFeatures(cell)",
        "    fpath <- file.path(work_dir, 'features.tsv')",
        "    utils::write.table(feats, fpath, sep = '\\t', quote = FALSE, row.names = FALSE)",
        "    KEYS <<- unique(featureKey(feats$feature_id, NA)$key)",
        "    one <- list(files = cell$files, design = cell$design,",
        "                contrasts = vapply(cell$contrasts, paste, '', collapse = '-'),",
        "                contrast_pairs = cell$contrasts, assays = cell$assays,",
        "                import_tool = if (cell$omic == 'metabolomics') 'msdial',",
        "                uniprot_cache = cell$uniprot_cache,",
        "                accessions = cell$annotation$accessions,",
        "                forced_inputs = cell$forced_inputs, features_path = fpath)",
        "    list(units = stats::setNames(rep(list(one), length(UNITS)), UNITS), features = feats)",
        "}",
        "buildTruth <- function(data_dir) {",
        "    list(primary = data.frame(key = KEYS, class = 'all_features', expected = 0,",
        "                              truth_status = 'scored'))",
        "}",
        "featureKey <- function(feature_id, arm) {",
        "    fid <- as.character(feature_id)",
        "    if (!is.null(LOOKUP)) {",
        "        by_name <- LOOKUP$GeneID[match(fid, LOOKUP$GeneName)]",
        "        dup <- fid %in% LOOKUP$GeneName[duplicated(LOOKUP$GeneName)]",
        "        key <- ifelse(fid %in% LOOKUP$GeneID, fid, by_name)",
        "        status <- ifelse(fid %in% LOOKUP$GeneID, 'ok', ifelse(dup, 'ambiguous',",
        "                         ifelse(is.na(by_name), 'unmapped', 'ok')))",
        "        return(data.frame(key = key, status = status,",
        "                          reason = ifelse(status == 'ok', NA, status)))",
        "    }",
        "    key <- sub('^(sp|tr)[|]([^|]+)[|].*$', '\\\\2', fid)",
        "    status <- ifelse(grepl(';', fid, fixed = TRUE), 'ambiguous', 'ok')",
        "    data.frame(key = key, status = status, reason = ifelse(status == 'ok', NA, 'group'))",
        "}",
        "e2eInTaxon <- function(a) rep(TRUE, length(a))",
        "nullArm <- function(unit, inputs) {",
        "    f <- utils::read.delim(inputs$features_path, stringsAsFactors = FALSE,",
        "                           colClasses = 'character')",
        "    data.frame(feature_id = f$feature_id, assay = f$assay,",
        "               contrast = readLines(CONTRAST_FILE)[[1L]], log2fc = 0, p = 0.5,",
        "               s_own = 0.5, call_statistic = 'none', called_own = FALSE,",
        "               mean_expression = 20)",
        "}",
        "refusingArm <- function(unit, inputs) {",
        "    if (unit == 'u04') stop('this arm refuses unit u04')",
        "    nullArm(unit, inputs)",
        "}",
        "e2eHook <- function(manager, results, summary_dir, tables, unit, arm, out_dir) {",
        "    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)",
        "    utils::write.table(data.frame(arm = arm, has_results = as.character(!is.null(results))),",
        "                       file.path(out_dir, sprintf('hook-%s-%s.tsv', unit, arm)),",
        "                       sep = '\\t', quote = FALSE, row.names = FALSE)",
        "}"
    ), code)
    writeTicketRegimeGenerator(file.path(dir, "gen.R"))
    code
}

e2eDeclaration <- function(route, contrast) {
    metric <- function(id, cmp, extra = list()) {
        utils::modifyList(list(
            id = id, `function` = "valMedianAbsError", args = list(n_min = 1L),
            class = "accuracy", views = list("intersection"), truth = "primary",
            call_rule = "own", higher_is_better = FALSE, contrast = contrast, margin = 0.5,
            conf = 0.95, comparisons = list(as.list(cmp)),
            interval = list(method = "paired_t",
                            accept_flagged = list(label = "paired t over 4 units, flagged (C2)",
                                                  ledger_row = "VAL-39"))), extra)
    }
    prot <- startsWith(route, "prot_")
    list(
        ticket = "KIT-E2E", stage = "counted",
        datasets = list(list(id = "d1", route = route, contrasts = list(contrast), seed = 159500L,
                             lifecycle_units = list("u01"),
                             not_under_test = list(arms = list(), units = list()),
                             truth = list("primary"),
                             annotation = if (prot) list(expected = "full", min_share = 1,
                                                         in_taxon = "e2eInTaxon"))),
        arms = list(list(id = "wf", kind = "workflow", builder = "valArmShipped", route = route),
                    list(id = "null", kind = "custom", `function` = "nullArm"),
                    list(id = "refusing", kind = "custom", `function` = "refusingArm")),
        hooks = list(post_cell = "e2eHook", summary = FALSE),
        metrics = list(
            metric("mae_nw", c("wf", "null"), list(primary_family = "fam")),
            metric("mae_ref", c("refusing", "null"), list(interval = list(method = "paired_t"))),
            metric("mae_w", c("wf", "null"), list(within_run = list(B = 2000L, seed = 159501L)))),
        primary_family = list(fam = list(m = 1L, members = list(
            a = list(all_of = list("mae_nw:wf-null"))))),
        refusal_tolerance = list(),
        coverage_regimes = list(list(id = "e2e-pt3", method = "paired_t", class = "accuracy",
                                     g = 3L, levels = list(0.95), generator = "gen.R",
                                     simulations = 2000L, seed = 159502L)),
        outputs = list(list(path = "ticket/hook-{unit}-{arm}.tsv", key = "none", max_rows = 1L,
                            lineage = list("corpus"),
                            columns = list(list(name = "arm", type = "label"),
                                           list(name = "has_results", type = "label")),
                            labels = list(arm = list("wf"),
                                          has_results = list("TRUE", "FALSE")))),
        blind = list(seed = 159503L), ceilings = list(cell_minutes = 30)
    )
}

e2eRun <- function(route) {
    root <- driverPkgRoot()
    skip_if(is.na(root), "the package checkout is not known (run through val-tests.R)")
    dir <- withr::local_tempdir(.local_envir = parent.frame())
    code <- e2eTicket(dir, route, root)
    units <- repairUnits(4L)
    man <- repairManifest()
    writeLines("unknown", file.path(dir, "contrast.txt"))
    decl <- e2eDeclaration(route, "unknown")
    # A blind lifecycle request on u01 with the workflow arm alone: it also gives the route's
    # contrast name, read from the cell's standard DA table.
    life <- valRunPhase("lifecycle", decl, "d1", units, "life", code, dir,
                        file.path(dir, "w0"), file.path(dir, "o0"), chunks = 1L, arms = "wf",
                        manifest = man, pkg_root = root)
    da <- list.files(file.path(dir, "w0", "cells"), "^da\\.tsv$", recursive = TRUE,
                     full.names = TRUE)
    contrast <- unique(utils::read.delim(da[[1L]], stringsAsFactors = FALSE)$contrast)[[1L]]
    writeLines(contrast, file.path(dir, "contrast.txt"))
    decl <- e2eDeclaration(route, contrast)
    tc <- valTicketCoverage(decl, dir)
    ev <- file.path(dir, "ev")
    dir.create(file.path(ev, "runs"), recursive = TRUE)
    utils::write.table(tc$coverage, file.path(ev, "ticket-coverage.tsv"), sep = "\t",
                       quote = FALSE, row.names = FALSE, na = "NA")
    decl$ticket_coverage <- list(path = "ticket-coverage.tsv",
                                 sha256 = sha256File(file.path(ev, "ticket-coverage.tsv")))
    regimes <- valClassifierRegimes(decl, ev)
    for (lab in c("r1", "r2")) {
        ids <- if (lab == "r1") c("u01", "u02") else c("u03", "u04")
        o <- valRunPhase("counted", decl, "d1", units, lab, code, dir,
                         file.path(dir, paste0("w", lab)), file.path(dir, paste0("o", lab)),
                         chunks = 4L, unit_ids = ids, manifest = man, regimes = regimes,
                         ticket_coverage = tc$coverage, pkg_root = root)
        file.copy(o, file.path(ev, "runs"), recursive = TRUE)
    }
    s <- valScorePhase(decl, "d1", units, ev, c("r1", "r2"), "s1", file.path(dir, "os"), man,
                       regimes = regimes)
    list(life = life, r1 = file.path(ev, "runs", "r1"), r2 = file.path(ev, "runs", "r2"),
         score = s, route = route)
}

e2eExpectations <- function(x) {
    route <- x$route
    lm <- utils::read.delim(file.path(x$life, "metrics.tsv"), stringsAsFactors = FALSE)
    expect_true(all(lm$blinded))
    expect_false(file.exists(file.path(x$life, "scored.tsv.gz")))
    expect_false(file.exists(file.path(x$life, "results.md")))
    for (r in c(x$r1, x$r2)) {
        for (f in c("cells.tsv", "calls.tsv", "metrics.tsv", "comparisons.tsv", "coverage.tsv",
                    "tables.sha256", "not-committed.sha256", "scored.tsv.gz", "results.md",
                    "attribution.tsv", "allowed-paths.txt", "provenance.json",
                    "outputs.sha256")) {
            expect_true(file.exists(file.path(r, f)), info = paste(route, basename(r), f))
        }
        expect_true(length(list.files(r, sprintf("^ledger-%s-wf-u0[1-4]\\.tsv$", route))) == 2L,
                    info = route)
        prov <- jsonlite::read_json(file.path(r, "provenance.json"))
        expect_identical(length(prov$cells), 6L)
        nin <- prov$scoring$not_in_universe
        expect_true(all(vapply(nin, function(z) z$n == 0L, logical(1))), info = route)
        expect_true(valCheckProvenance(file.path(r, "provenance.json")))
        status <- checkEvidenceShell(r)
        if (!is.na(status)) expect_identical(status, 0L, info = paste(route, basename(r)))
        hooks <- list.files(file.path(r, "ticket"), "^hook-")
        expect_identical(length(hooks), 2L, info = route)
    }
    m1 <- utils::read.delim(file.path(x$r1, "metrics.tsv"), stringsAsFactors = FALSE)
    expect_true(all(is.finite(m1$value[m1$arm == "wf" & m1$metric == "mae_nw"])), info = route)
    sc <- utils::read.delim(file.path(x$r1, "scored.tsv.gz"), stringsAsFactors = FALSE)
    expect_identical(names(sc), msrval:::SCORED_COLUMNS)
    if (route == "metab_combined") {
        expect_setequal(unique(sc$assay), c("LCMS_Pos", "GCMS"))
        wf <- sc[sc$arm == "wf" & sc$testable, ]
        for (a in unique(wf$assay)) {
            w <- wf[wf$assay == a & wf$unit == "u01", ]
            expect_equal(w$s_bh, stats::p.adjust(w$p, "BH"), info = a)
        }
    }
    c2 <- utils::read.delim(file.path(x$r2, "cells.tsv"), stringsAsFactors = FALSE)
    expect_identical(c2$exit[c2$unit == "u04" & c2$arm == "refusing"], "refused")
    cmp <- utils::read.delim(file.path(x$score, "comparisons.tsv"), stringsAsFactors = FALSE)
    nw <- cmp[cmp$component %in% "mae_nw:wf-null", ]
    expect_identical(nw$regime, "C2")
    expect_false(is.na(nw$outcome))
    ref <- cmp[cmp$component %in% "mae_ref:refusing-null", ]
    expect_identical(ref$n_refused, 1L)
    expect_identical(ref$regime, "e2e-pt3")
    expect_false(is.na(ref$outcome))
    expect_match(ref$flags, "refusal_substitute_required", fixed = TRUE)
    w <- cmp[cmp$component %in% "mae_w:wf-null", ]
    expect_identical(lengths(regmatches(w$flags, gregexpr("within_run_reported=", w$flags))), 4L)
    expect_false(is.na(cmp$outcome[cmp$member %in% "a"]))
    txt <- readLines(file.path(x$score, "results.md"))
    expect_true(any(grepl("mae_nw: wf - null.*regime C2.*measured coverage 0.812", txt)))
    expect_true(any(grepl("FLAG: worst-condition miss rate at least twice nominal", txt,
                          fixed = TRUE)))
    expect_true(any(grepl("fam: m = 1, nominal family level 0.95", txt, fixed = TRUE)))
}

test_that("E2E-prot: prot_dia from cells to verdicts", {
    e2eExpectations(e2eRun("prot_dia"))
})

test_that("E2E-metab: metab_combined (two assays, shared ids) from cells to verdicts", {
    e2eExpectations(e2eRun("metab_combined"))
})

test_that("E2E-trx: trx_screen from cells to verdicts", {
    e2eExpectations(e2eRun("trx_screen"))
})
