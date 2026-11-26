# KIT-268 (owner's decision VAL-312, wording VAL-314): where Storey's pi0 cannot be estimated for a comparison and
# Benjamini-Hochberg values sit in the q-value column, the tab's notice names the comparisons (proteomics through
# warnDaQvalueFallback; lipidomics and metabolomics through completeOmicDaRun), the analysis record's fdr_fallback_reason
# carries the one added sentence, and the statistical claim of that sentence holds on the runner's qvalue (T6).
# T1, T2 and the notice names are guarded with exists() so they fail rather than error before the change.

k268Sentence <- paste("For these comparisons the q-value column holds Benjamini-Hochberg adjusted p-values, which equal",
    "Storey's q-values with the proportion of unchanged features (pi0) set to 1; because pi0 could not be estimated for these",
    "comparisons, the values are larger than Storey's q-values would be by a factor of 1/pi0 of unknown size (1 or more),",
    "so they are more conservative than, and not on the same scale as, the q-values of the other comparisons.")

k268Fn <- function(name) get0(name, asNamespace("MultiScholaR"), inherits = FALSE)

test_that("KIT-268 T1: the proteomics notice names the comparisons that fell back and none that did not", {
    notify <- k268Fn("warnDaQvalueFallback")
    expect_true(is.function(notify))
    if (!is.function(notify)) return(invisible())
    manifest <- list(
        results = list(`A-B` = list(qvalue_warnings = list(`A-B` = TRUE)), `C-D` = list(qvalue_warnings = list())),
        contrasts = data.frame(contrasts = c("A-B", "C-D"), friendly_names = c("A_vs_B", "C_vs_D"), stringsAsFactors = FALSE))
    notes <- list()
    testthat::local_mocked_bindings(showNotification = function(ui, ..., type = "default") {
        notes[[length(notes) + 1L]] <<- list(text = ui, type = type)
        invisible(NULL)
    }, .package = "shiny")
    notify(manifest)
    expect_length(notes, 1L)
    expect_identical(notes[[1L]]$type, "warning")
    expect_true(grepl("A_vs_B", notes[[1L]]$text, fixed = TRUE))
    expect_false(grepl("C_vs_D", notes[[1L]]$text, fixed = TRUE))
    notes <- list()
    manifest$results[["A-B"]]$qvalue_warnings <- list()
    notify(manifest)
    expect_length(notes, 0L)
})

test_that("KIT-268 T2: the mapping of a fallback key to the comparison's friendly name", {
    names_of <- k268Fn("daQvalueFallbackNames")
    expect_true(is.function(names_of))
    if (!is.function(names_of)) return(invisible())
    contrasts <- data.frame(contrasts = c("groupA-groupB", "X=Y"), friendly_names = c("A_vs_B", "Z_vs_W"),
        stringsAsFactors = FALSE)
    # The package's contrast form is <friendly>=<contrast>: the text before "="; an exact match on `contrasts` wins;
    # an unmatched key stays as it is.
    expect_identical(names_of(c("B_vs_A=groupB-groupA", "groupA-groupB", "X=Y", "other"), contrasts),
        c("B_vs_A", "A_vs_B", "Z_vs_W", "other"))
    expect_identical(names_of("groupA-groupB", NULL), "groupA-groupB")
})

k268Completion <- function(omic, results) {
    config <- omicGuiConfig(omic)
    if (identical(omic, "lipidomics")) {
        object <- module_ci_lipid_da_object(layout = "triple")
        workflow <- module_ci_lipid_da_workflow(object)
        tab <- module_ci_lipid_da_data_env()
        capture <- module_ci_lipid_da_capture_ui()
    } else {
        object <- module_ci_metab_da_object(layout = "combined")
        workflow <- module_ci_metab_da_workflow(object)
        tab <- module_ci_metab_da_data_env()
        capture <- module_ci_metab_da_capture_ui()
    }
    testthat::local_mocked_bindings(updateSelectInput = capture$updateSelectInput, showNotification = capture$showNotification,
        removeNotification = capture$removeNotification, .package = "shiny")
    suppressMessages(completeOmicDaRun(results = results, tab = tab, workflowData = workflow, session = "k268",
        experimentPaths = list(da_output_dir = NULL, publication_graphs_dir = NULL), da = config$da,
        binding = config$binding$binding))
    list(text = vapply(capture$notifications, function(n) as.character(n[[1L]]), ""),
        type = vapply(capture$notifications, function(n) as.character(n$type %||% "default"), ""))
}

test_that("KIT-268 T3: the lipidomics and metabolomics tab names each comparison that fell back, with its assay", {
    cases <- list(
        lipidomics = list(results = module_ci_lipid_da_results_list("mixed"), assay = "LCMS_Neg"),
        metabolomics = list(results = module_ci_metab_da_results_list("mixed"), assay = "LCMS_Pos"))
    for (omic in names(cases)) {
        results <- cases[[omic]]$results
        assay <- cases[[omic]]$assay
        # A planted fallback of one comparison in one assay.
        results$qvalue_warnings <- stats::setNames(list(list(`groupKO-groupWT` = TRUE)), assay)
        got <- k268Completion(omic, results)
        hit <- grep(paste0("KO_vs_WT in ", assay), got$text, fixed = TRUE)
        expect_length(hit, 1L)
        expect_identical(got$type[hit], "warning")
        # No fallback, no notice (transcriptomics' list is always empty).
        results$qvalue_warnings <- list()
        none <- k268Completion(omic, results)
        expect_length(grep("Benjamini-Hochberg", none$text, fixed = TRUE), 0L)
    }
    # The tied_pvalues case: the DA kernel's own list (assay -> comparison) names every comparison in every assay.
    module_ci_lipid_da_bind_orchestration_stubs("tied_pvalues")
    suppressMessages(capture.output(lipid <- runLipidsDA(theObject = module_ci_lipid_da_object(layout = "triple"),
        contrasts_tbl = module_ci_lipid_da_contrasts("multi_group"), formula_string = "~ 0 + group", da_q_val_thresh = 0.05,
        treat_lfc_cutoff = 0)))
    got <- k268Completion("lipidomics", lipid)
    expect_true(any(grepl("KO_vs_WT, RES_vs_WT in LCMS_Pos", got$text, fixed = TRUE)))
    module_ci_metab_da_bind_orchestration_stubs("tied_pvalues")
    suppressMessages(capture.output(metab <- runMetabolitesDA(theObject = module_ci_metab_da_object(layout = "combined"),
        contrasts_tbl = module_ci_metab_da_contrasts("multi_group"), formula_string = "~ 0 + group", da_q_val_thresh = 0.05,
        treat_lfc_cutoff = 0)))
    got <- k268Completion("metabolomics", metab)
    expect_true(any(grepl("KO_vs_WT, RES_vs_WT in LCMS_Pos", got$text, fixed = TRUE)))
})

test_that("KIT-268 T5: the record's reason for a fallback carries the existing text and the added sentence", {
    fallback <- workflowRecordDaFdrMethod("DA-contrast-output", list(qvalue_warnings = list(a = TRUE)))
    expect_identical(fallback$fdr_method, "BH")
    expect_match(fallback$fdr_fallback_reason, "fdr_qvalue holds Benjamini-Hochberg adjusted p-values", fixed = TRUE)
    expect_match(fallback$fdr_fallback_reason, "1/pi0", fixed = TRUE)
    expect_match(fallback$fdr_fallback_reason, "unknown size", fixed = TRUE)
    expect_match(fallback$fdr_fallback_reason, "not on the same scale", fixed = TRUE)
    plain <- workflowRecordDaFdrMethod("DA-contrast-output", list(qvalue_warnings = list()))
    expect_identical(plain, list(fdr_method = "qvalue"))
})

test_that("KIT-268 T6: Benjamini-Hochberg is Storey's q divided by pi0 at every index, ties included; the trigger", {
    skip_if_not_installed("qvalue")
    set.seed(268)
    p <- c(stats::runif(400), stats::runif(100, 0, 0.001))
    tied <- p
    tied[2:5] <- tied[1L]
    for (v in list(p, tied)) {
        qv <- qvalue::qvalue(v)
        # With pi0 = 1 the identity cannot tell q = pi0 x BH from q = BH.
        expect_lt(qv$pi0, 0.95)
        bh <- stats::p.adjust(v, "BH")
        expect_equal(bh, qv$qvalues / qv$pi0, tolerance = 1e-10)
        expect_true(all(qv$qvalues <= bh + 1e-12))
    }
    # The usual trigger of the fallback: no p-value at or above 0.95 stops qvalue()'s smoother (VAL-237, VAL-295).
    expect_error(qvalue::qvalue(stats::runif(50, 0, 0.9)))
})
