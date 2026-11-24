# The G1(b) alternatives (KIT-159 r1.14; S0, then L1a): export-then-load on one route per module
# server, DIA_LIMPA without "Compute RUV-III Covariates"; their comparison with the accepted
# script's DA table (valCompareAlternative); the S0 plan and the gate that ends S0 before any route
# when the new pieces' tests did not pass (VAL-54).

altIds <- function(route, alternative = NULL) {
    vapply(valActionScript(route, alternative = alternative)$steps, function(s) s$id, character(1))
}

altPresses <- function(route, alternative = NULL) {
    unlist(lapply(valActionScript(route, alternative = alternative)$steps, function(s) s$press))
}

EXPORT_LOAD_BUTTONS <- list(
    prot_dia = c("normalization-export_filtered_session",
                 "differential_expression-load_filtered_session"),
    lipid_canonical = c("norm-export_session", "de-load_filtered_session"),
    metab_lc = c("norm-export_session", "de-load_filtered_session")
)

test_that("T-ALT-1: export-then-load: export after the last normalisation step, load before DA", {
    for (route in names(EXPORT_LOAD_BUTTONS)) {
        accepted <- altIds(route)
        alt <- altIds(route, "export_load")
        expect_identical(alt, append(accepted, c("export_session", "load_session"),
                                     after = match("da", accepted) - 1L), info = route)
        expect_setequal(setdiff(altPresses(route, "export_load"), altPresses(route)),
                        EXPORT_LOAD_BUTTONS[[route]])
        steps <- valActionScript(route, alternative = "export_load")$steps
        names(steps) <- alt
        last_norm <- steps[[match("export_session", alt) - 1L]]$declared$head
        expect_identical(steps$export_session$declared$head, last_norm, info = route)
        expect_identical(steps$load_session$declared$head, last_norm, info = route)
    }
    expect_error(valActionScript("trx_screen", alternative = "export_load"),
                 class = "val_spec_invalid")
    expect_error(valActionScript("prot_dia", alternative = "bypass"), class = "val_spec_invalid")
})

test_that("T-ALT-2: every button the alternative adds is rendered by the product's tab UI", {
    ns <- asNamespace("MultiScholaR")
    html <- function(fn, id) paste(as.character(htmltools::renderTags(fn(id))$html), collapse = "\n")
    pages <- list(
        prot_dia = c(html(ns$mod_prot_norm_ui, "normalization"),
                     html(ns$mod_prot_da_ui, "differential_expression")),
        lipid_canonical = c(html(ns$mod_lipid_norm_ui, "norm"), html(ns$mod_lipid_da_ui, "de")),
        metab_lc = c(html(ns$mod_metab_norm_ui, "norm"), html(ns$mod_metab_da_ui, "de"))
    )
    for (route in names(pages)) {
        added <- setdiff(altPresses(route, "export_load"), altPresses(route))
        expect_length(added, 2L)
        for (id in added) {
            expect_true(any(grepl(paste0("id=\"", id, "\""), pages[[route]], fixed = TRUE)),
                        info = paste(route, id))
        }
    }
})

test_that("T-ALT-3: DIA_LIMPA without the covariates button drops only that step", {
    expect_identical(altIds("prot_dia_limpa", "no_covariates"),
                     setdiff(altIds("prot_dia_limpa"), "ruv_covariates"))
    expect_false("normalization-apply_ruv_covariates" %in%
                     altPresses("prot_dia_limpa", "no_covariates"))
    expect_error(valActionScript("prot_dia", alternative = "no_covariates"),
                 class = "val_spec_invalid")
    # the accepted scripts are unchanged by the argument's default
    expect_identical(altIds("prot_dia_limpa", NULL),
                     vapply(valActionScript("prot_dia_limpa")$steps, `[[`, "", "id"))
})

daTable <- function() {
    data.frame(feature_id = c("a", "b", "c"), assay = "all", contrast = "KO-WT",
               log2fc = c(1, -0.5, 0.1), p = c(0.001, 0.02, 0.6), s_own = c(0.01, 0.04, 0.8),
               call_statistic = "qvalue", called_own = c(TRUE, TRUE, FALSE),
               mean_expression = c(20, 21, 22), stringsAsFactors = FALSE)
}

test_that("T-ALT-4: valCompareAlternative reports digest, largest difference and changed calls", {
    a <- daTable()
    same <- valCompareAlternative(a, a)
    expect_identical(same$identical, TRUE)
    expect_identical(same$digest_a, same$digest_b)
    expect_identical(same$max_abs_diff, 0)
    expect_identical(c(same$called_changed, same$keys_one_side, same$na_mismatch), c(0L, 0L, 0L))
    b <- a
    b$log2fc[[2L]] <- -1
    b$called_own[[3L]] <- TRUE
    d <- valCompareAlternative(a, b)
    expect_false(d$identical)
    expect_equal(d$max_abs_diff, 0.5)
    expect_identical(d$called_changed, 1L)
    expect_identical(c(d$rows_a, d$rows_b), c(3L, 3L))
    expect_identical(c(d$called_a, d$called_b), c(2L, 3L))
    c2 <- rbind(a, transform(a[1L, ], feature_id = "z"))
    c2$p[[3L]] <- NA
    e <- valCompareAlternative(a, c2)
    expect_identical(e$keys_one_side, 1L)
    expect_identical(e$na_mismatch, 1L)
    expect_identical(e$max_abs_diff, 0)
})

test_that("T-S0-1: the S0 plan: ten accepted routes, four alternatives, the two-species route", {
    plan <- valS0Plan()
    expect_identical(nrow(plan), 15L)
    expect_false(anyDuplicated(plan$cell) > 0L)
    accepted <- plan[is.na(plan$alternative) & plan$cell != PARITY_TWO_SPECIES, ]
    expect_setequal(accepted$route, driverRoutes())
    expect_identical(accepted$cell, accepted$route)
    alt <- plan[!is.na(plan$alternative), ]
    expect_setequal(paste(alt$route, alt$alternative),
                    c("prot_dia export_load", "lipid_canonical export_load",
                      "metab_lc export_load", "prot_dia_limpa no_covariates"))
    expect_identical(alt$compare_with, alt$route)
    two <- plan[plan$cell == PARITY_TWO_SPECIES, ]
    expect_identical(c(two$route, two$corpus, two$compare_with),
                     c("prot_dia", PARITY_TWO_SPECIES, "prot_dia"))
    expect_true(all(plan$corpus %in% valCorpusSpec()$corpus))
    expect_true(all(plan$corpus[plan$cell != PARITY_TWO_SPECIES] ==
                        vapply(plan$route[plan$cell != PARITY_TWO_SPECIES], valCorpusFor, "",
                               USE.NAMES = FALSE)))
})

test_that("T-S0-2: S0 ends before any route unless every named test file ran and passed", {
    dir <- withr::local_tempdir()
    tsv <- file.path(dir, "harness-tests.tsv")
    rows <- function(failed = 0L, error = 0L, skipped = 0L) {
        data.frame(file = c("test-corpus.R", "test-alternatives.R"), test = c("x", "y"),
                   expectations = 3L, passed = 3L - failed, failed = failed, skipped = skipped,
                   error = error, warnings = 0L, seconds = 1, stringsAsFactors = FALSE)
    }
    put <- function(x) utils::write.table(x, tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    files <- c("corpus", "alternatives")
    put(rows())
    expect_true(isTRUE(valS0TestsPassed(tsv, files)))
    put(rows(failed = c(1L, 0L)))
    expect_false(isTRUE(valS0TestsPassed(tsv, files)))
    put(rows(error = c(0L, 1L)))
    expect_false(isTRUE(valS0TestsPassed(tsv, files)))
    put(rows(skipped = c(1L, 0L)))
    expect_false(isTRUE(valS0TestsPassed(tsv, files)))
    put(rows()[1L, ])
    expect_false(isTRUE(valS0TestsPassed(tsv, files)))
    expect_false(isTRUE(valS0TestsPassed(file.path(dir, "absent.tsv"), files)))
    put(rows(failed = c(1L, 0L)))
    out <- file.path(dir, "s0")
    expect_error(valRunS0(out, driverPkgRoot(), list(path = tsv, files = files),
                          work_dir = file.path(dir, "work")),
                 class = "val_spec_invalid")
    expect_false(dir.exists(file.path(dir, "work")))
    expect_length(list.files(out, recursive = TRUE), 0L)
})
