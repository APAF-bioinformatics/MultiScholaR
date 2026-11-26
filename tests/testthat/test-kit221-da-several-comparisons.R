# KIT-221 (VAL-188 (3), VAL-193): scripted proteomics DA given several comparisons in one call returned the first
# comparison only. differentialAbundanceAnalysis() fitted every comparison but built its derived tables and plots from
# contrasts_results$results[[1]] (R/func_prot_da_model_methods.R:430, :433 at 176f9d3ff); protein_deAnalysisWrapperFunction()
# read contrasts_tbl[, 1][[1]], the first row of a base data.frame (R/func_prot_da_compat.R:115). Under option (A) both
# return every comparison, each comparison's values the same as a call with it alone, read the comparisons as text and
# refuse a non-data-frame with the class prot_da_contrasts_invalid. Each W test first asserts that its call completed
# (`expect_null` on the caught message) and returns on a failure, so on the tests-only commit every failing row fails
# by an assertion and never by an error. Synthetic data only.
#
# Fixture builders, copied and renamed from tests/testthat/test-kit218-da-wrapper-tail.R at 176f9d3ff (which copies
# test-kit200-prot-da-no-value-proteins.R): kit218Quiet :15-19, kit218WithSeed :21-22, kit218Matrix :25-35,
# kit218Design :37-41, kit218Object :54-61, kit218MockPlots :86-91. The limpa object and its fake
# dpcDE are copied from test-prot-07d-da-model-stats-shared.R:670-697 (renamed, a 60 x 12 matrix, no W column). The
# RUV-factors-in-the-model builders are copied from test-prot-ruv-factors-in-model.R:8-52 (ruvModelFixture,
# ruvModelObject, ruvModelPackageState), the three groups being the fixture's twelve runs relabelled A, B, C.

kit221Quiet <- function(expr) {
    value <- NULL
    utils::capture.output(value <- suppressWarnings(suppressMessages(expr)))
    value
}

kit221WithSeed <- function(seed, code) withr::with_seed(seed, code, .rng_kind = "Mersenne-Twister",
    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

# Proteins P1..P60, levels N(22, 2^2), residual sd 0.3, P11..P16 with a group-B effect of +1 (KIT-221 section 4).
kit221Matrix <- function(sizes = c(A = 4L, B = 4L, C = 4L), seed = 200001L, proteins = 60L) {
    groups <- rep(names(sizes), sizes)
    x <- kit221WithSeed(seed, {
        level <- stats::rnorm(proteins, 22, 2)
        effect <- as.numeric(seq_len(proteins) %in% 11:16)
        level + outer(effect, as.numeric(groups == "B")) +
            matrix(stats::rnorm(proteins * length(groups), sd = 0.3), proteins)
    })
    dimnames(x) <- list(paste0("P", seq_len(proteins)), paste0("S", seq_along(groups)))
    x
}

kit221Design <- function(sizes = c(A = 4L, B = 4L, C = 4L)) {
    runs <- paste0("S", seq_len(sum(sizes)))
    data.frame(Run = runs, group = rep(names(sizes), sizes), replicates = runs, row.names = runs,
        stringsAsFactors = FALSE)
}

kit221Object <- function(x, design) {
    rownames(design) <- NULL
    ProteinQuantitativeData(protein_quant_table = data.frame(Protein.Ids = rownames(x), x, row.names = NULL,
        check.names = FALSE, stringsAsFactors = FALSE), design_matrix = design, sample_id = "Run", group_id = "group",
        technical_replicate_id = "replicates", protein_id_column = "Protein.Ids",
        protein_id_table = data.frame(Protein.Ids = rownames(x), gene_names = paste0("GENE_", rownames(x)),
            stringsAsFactors = FALSE), args = list())
}

kit221CompleteObject <- function() kit221Object(kit221Matrix(), kit221Design())

# The limpa route: an object carrying the quantified EList, as the existing test of the DA helper builds it.
kit221LimpaObject <- function() {
    x <- kit221Matrix()
    ids <- rownames(x)
    nobs <- matrix(2L, nrow(x), ncol(x), dimnames = dimnames(x))
    elist <- methods::new("EList", list(E = x,
        genes = data.frame(NPrec = rep(2L, nrow(x)), PropObs = rep(1, nrow(x)), row.names = ids),
        other = list(standard.error = matrix(0.2, nrow(x), ncol(x), dimnames = dimnames(x)), n.observations = nobs)))
    object <- kit221CompleteObject()
    object@args <- list(limpa_dpc_quant_results = list(quantified_elist = elist,
        standard_errors = elist$other$standard.error))
    object
}

kit221MockLimpa <- function(env = parent.frame()) {
    fake_dpcDE <- function(y, design, plot = TRUE, ...) limma::lmFit(y$E, design)
    testthat::local_mocked_bindings(dpcDE = fake_dpcDE, .package = "limpa", .env = env)
}

# The RUV-factors-in-the-model route: ruvModelFixture's data (twelve runs), the runs relabelled A, B, C.
kit221RuvFixture <- function(seed = 20261118L, proteins = 300L, per_group = 6L, missing = 0.25) {
    kit221WithSeed(seed, {
        groups <- rep(c("control", "treated"), each = per_group)
        n_samples <- length(groups)
        mu <- stats::rnorm(proteins, 22, 2)
        s <- sqrt(4 * 0.25^2 / stats::rchisq(proteins, 4))
        effect <- numeric(proteins)
        where <- sample.int(proteins, round(0.2 * proteins))
        effect[where] <- sample(rep(c(1, -1), length.out = length(where)))
        scale1 <- function(v) { v <- v - mean(v); v / stats::sd(v) }
        w <- scale1(stats::lm.fit(stats::model.matrix(~ factor(groups)), stats::rnorm(n_samples))$residuals)
        tau <- sqrt(3 * 4 * 0.25^2 / stats::qchisq(0.5, 4))
        x <- matrix(mu, proteins, n_samples) + outer(effect, as.numeric(groups == "treated")) +
            outer(stats::rnorm(proteins, 0, tau), w) +
            matrix(stats::rnorm(proteins * n_samples), proteins, n_samples) * s
        if (missing > 0) {
            lod <- stats::uniroot(function(l) mean(stats::plogis(l - x)) - missing, c(0, 60))$root
            x[stats::runif(length(x)) < stats::plogis(lod - x)] <- NA_real_
        }
        dimnames(x) <- list(sprintf("PG_%05d", seq_len(proteins)), sprintf("run_%03d", seq_len(n_samples)))
        keep <- rowSums(is.finite(x[, groups == "control"])) >= per_group / 2 &
            rowSums(is.finite(x[, groups == "treated"])) >= per_group / 2
        design <- data.frame(Run = colnames(x), group = rep(c("A", "B", "C"), each = n_samples / 3L),
            replicates = seq_len(n_samples), tech_rep_group = paste0("u", seq_len(n_samples)), stringsAsFactors = FALSE)
        list(x = x[keep, , drop = FALSE], design = design)
    })
}

kit221RuvState <- function(fixture, k = 1L) {
    quant <- as.data.frame(c(list(Protein.Group = rownames(fixture$x)), as.data.frame(fixture$x)), check.names = FALSE,
        stringsAsFactors = FALSE)
    normalised <- methods::new("ProteinQuantitativeData", protein_quant_table = quant, protein_id_column = "Protein.Group",
        protein_id_table = data.frame(Protein.Group = rownames(fixture$x), stringsAsFactors = FALSE),
        design_matrix = fixture$design, sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
        args = list(globalParameters = list(workflow_type = "LFQ")))
    factors <- computeProtRuvModelFactors(normalised, ruv_k = k, ruv_k_rule = "user", ruv_grouping_variable = "group")
    corrected <- suppressMessages(ruvIII_C_Varying(protRuvGroupReplicateObject(normalised), ruv_grouping_variable = "group",
        ruv_number_k = factors$k, ctrl = factors$controls))
    attachProtRuvModelFactors(corrected, normalised, factors)
}

# ---- KIT-221's own helpers -----------------------------------------------------------------------------------------------

kit221B <- "B_vs_A=groupB-groupA"
kit221C <- "C_vs_A=groupC-groupA"

# The comparisons as a base data.frame (the default), a tibble, or a data.frame whose column is a factor.
kit221Table <- function(strings, kind = c("data.frame", "tibble", "factor")) {
    switch(match.arg(kind),
        data.frame = data.frame(contrasts = strings, stringsAsFactors = FALSE),
        tibble = tibble::tibble(contrasts = strings),
        factor = data.frame(contrasts = strings, stringsAsFactors = TRUE))
}

# The plots before the DA step are not this ticket's subject (as in KIT-218 and KIT-200).
kit221MockPlots <- function(env = parent.frame()) {
    testthat::local_mocked_bindings(plotRle = function(...) ggplot2::ggplot(),
        plotPca = function(...) ggplot2::ggplot(), plotNumOfValuesNoLog = function(...) ggplot2::ggplot(),
        .package = "MultiScholaR", .env = env)
}

kit221Da <- function(object, contrasts_tbl, row_id = "Protein.Ids", trend = TRUE, robust = TRUE) {
    tryCatch(list(value = kit221Quiet(differentialAbundanceAnalysis(object, contrasts_tbl = contrasts_tbl,
        formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, eBayes_trend = trend,
        eBayes_robust = robust, args_row_id = row_id, qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue"))),
        error = function(e) list(error = conditionMessage(e), condition = e))
}

kit221Wrapper <- function(object, contrasts_tbl) {
    tryCatch(list(value = kit221Quiet(protein_deAnalysisWrapperFunction(object, contrasts_tbl = contrasts_tbl,
        formula_string = "~ 0 + group", de_q_val_thresh = 0.05, treat_lfc_cutoff = 0, eBayes_trend = TRUE,
        eBayes_robust = TRUE, args_row_id = "Protein.Ids"))),
        error = function(e) list(error = conditionMessage(e), condition = e))
}

# The call completed. On the tests-only commit a failing call fails here and prints the caught message.
kit221Completed <- function(run, label) {
    expect_null(run$error, label = label)
    is.null(run$error)
}

# Per comparison and protein, log2FC, raw_pvalue and fdr_qvalue of `got` identical to the reference call's rows
# (`refs`: a named list, comparison -> that one-comparison call's significant_rows).
kit221SameRows <- function(got, refs, label) {
    expect_identical(sort(unique(got$comparison)), sort(names(refs)), label = paste(label, "comparisons"))
    for (comparison in sort(names(refs))) {
        a <- got[got$comparison == comparison, , drop = FALSE]
        b <- refs[[comparison]][refs[[comparison]]$comparison == comparison, , drop = FALSE]
        expect_identical(nrow(a), nrow(b), info = paste(label, comparison, "rows"))
        if (nrow(a) != nrow(b)) next
        id <- intersect(c("Protein.Ids", "Protein.Group"), names(a))[[1L]]
        a <- a[order(a[[id]]), , drop = FALSE]
        b <- b[order(b[[id]]), , drop = FALSE]
        expect_identical(a[[id]], b[[id]], info = paste(label, comparison, "ids"))
        for (column in c("log2FC", "raw_pvalue", "fdr_qvalue")) {
            expect_identical(a[[column]], b[[column]], info = paste(label, comparison, column))
        }
    }
}

# W1's expectations for a result holding both comparisons of the fixture (60 proteins each).
kit221ExpectBoth <- function(value) {
    rows <- value$significant_rows
    both <- c("B_vs_A", "C_vs_A")
    expect_identical(sort(unique(rows$comparison)), both)
    expect_identical(nrow(rows), 120L)
    expect_identical(nrow(value$contrasts_results_table), 120L)
    expect_setequal(unique(value$num_sig_da_molecules$comparison), both)
    expect_setequal(unique(value$num_sig_da_molecules_first_go$table$comparison), both)
    expect_identical(nrow(value$list_of_volcano_plots), 2L)
    expect_identical(nrow(value$da_proteins_wide), 60L)
    expect_true(all(c("log2FC:B_vs_A", "log2FC:C_vs_A") %in% names(value$da_proteins_wide)))
    expect_identical(sort(unique(value$da_proteins_long$comparison)), both)
}

kit221Singles <- function(call, object, kind = "data.frame", ...) {
    list(B_vs_A = call(object, kit221Table(kit221B, kind), ...), C_vs_A = call(object, kit221Table(kit221C, kind), ...))
}

# ---- W1, W1b, W2: differentialAbundanceAnalysis(), limma route ---------------------------------------------------------

test_that("KIT-221 W1: two comparisons in a data.frame give every table and plot both comparisons (limma)", {
    kit221MockPlots()
    run <- kit221Da(kit221CompleteObject(), kit221Table(c(kit221B, kit221C)))
    if (!kit221Completed(run, "W1 error")) return(invisible())
    kit221ExpectBoth(run$value)
    cat(sprintf("\nKIT-221 W1: names(da_proteins_long) = %s\n", paste(names(run$value$da_proteins_long), collapse = ", ")))
})

test_that("KIT-221 W1b: two comparisons in a tibble give both comparisons too (M1: the base message is recorded)", {
    kit221MockPlots()
    run <- kit221Da(kit221CompleteObject(), kit221Table(c(kit221B, kit221C), "tibble"))
    cat(sprintf("\nKIT-221 M1: %s\n", if (is.null(run$error)) "the call completed" else
        paste("caught message:", run$error)))
    if (!kit221Completed(run, "W1b error")) return(invisible())
    kit221ExpectBoth(run$value)
})

test_that("KIT-221 W2: each comparison of a two-comparison call equals a call with that comparison alone (limma)", {
    kit221MockPlots()
    object <- kit221CompleteObject()
    multi <- kit221Da(object, kit221Table(c(kit221B, kit221C)))
    singles <- kit221Singles(kit221Da, object)
    if (!kit221Completed(multi, "W2 two-comparison error")) return(invisible())
    if (!kit221Completed(singles$B_vs_A, "W2 B_vs_A error") || !kit221Completed(singles$C_vs_A, "W2 C_vs_A error")) return(invisible())
    kit221SameRows(multi$value$significant_rows, lapply(singles, function(run) run$value$significant_rows), "W2")
})

# ---- W3: protein_deAnalysisWrapperFunction() ----------------------------------------------------------------------------

test_that("KIT-221 W3: the sibling wrapper given a two-row data.frame returns both comparisons, each equal to its one-row tibble call", {
    kit221MockPlots()
    object <- kit221CompleteObject()
    multi <- kit221Wrapper(object, kit221Table(c(kit221B, kit221C)))
    if (!kit221Completed(multi, "W3 error")) return(invisible())
    value <- multi$value
    expect_identical(length(value$contrasts_results$results), 2L)
    expect_identical(sort(unique(value$significant_rows$comparison)), c("B_vs_A", "C_vs_A"))
    expect_identical(nrow(value$significant_rows), 120L)
    singles <- kit221Singles(kit221Wrapper, object, "tibble")
    if (!kit221Completed(singles$B_vs_A, "W3 B_vs_A slice error") ||
        !kit221Completed(singles$C_vs_A, "W3 C_vs_A slice error")) return(invisible())
    kit221SameRows(value$significant_rows, lapply(singles, function(run) run$value$significant_rows), "W3")
})

# ---- W4: the limpa route ---------------------------------------------------------------------------------------------------

test_that("KIT-221 W4: two comparisons on the limpa route give both comparisons, each equal to its one-comparison call", {
    skip_if_not_installed("limpa")
    kit221MockPlots()
    kit221MockLimpa()
    object <- kit221LimpaObject()
    multi <- kit221Da(object, kit221Table(c(kit221B, kit221C)), trend = FALSE, robust = FALSE)
    singles <- kit221Singles(kit221Da, object, trend = FALSE, robust = FALSE)
    if (!kit221Completed(multi, "W4 two-comparison error")) return(invisible())
    if (!kit221Completed(singles$B_vs_A, "W4 B_vs_A error") || !kit221Completed(singles$C_vs_A, "W4 C_vs_A error")) return(invisible())
    expect_identical(multi$value$da_model$model, "dpcDE")
    kit221ExpectBoth(multi$value)
    kit221SameRows(multi$value$significant_rows, lapply(singles, function(run) run$value$significant_rows), "W4")
})

# ---- W5: the RUV-factors-in-the-model route (three groups, one factor) ------------------------------------------------

test_that("KIT-221 W5: two comparisons on the RUV-factors-in-the-model route give both comparisons, each equal to its one-comparison call", {
    kit221MockPlots()
    state <- kit221RuvState(kit221RuvFixture(), 1L)
    multi <- kit221Da(state, kit221Table(c(kit221B, kit221C)), row_id = "Protein.Group")
    singles <- kit221Singles(kit221Da, state, row_id = "Protein.Group")
    if (!kit221Completed(multi, "W5 two-comparison error")) return(invisible())
    if (!kit221Completed(singles$B_vs_A, "W5 B_vs_A error") || !kit221Completed(singles$C_vs_A, "W5 C_vs_A error")) return(invisible())
    expect_identical(multi$value$da_model$tested_matrix, "uncorrected")
    expect_identical(sort(unique(multi$value$significant_rows$comparison)), c("B_vs_A", "C_vs_A"))
    # A GUARD, not a failing-first row: the S4 method keeps every comparison in contrasts_results (ticket 2.1), so this
    # passes on base and on head; the derived tables below are what the fix changes.
    expect_identical(length(multi$value$contrasts_results$results), 2L)
    expect_identical(nrow(multi$value$significant_rows), 2L * nrow(singles$B_vs_A$value$significant_rows))
    kit221SameRows(multi$value$significant_rows, lapply(singles, function(run) run$value$significant_rows), "W5")
})

# ---- W7, W8: a comparisons column stored as a factor, levels in the opposite order to the rows -----------------------

test_that("KIT-221 W7: a factor column gives each comparison its own statistics (limma, differentialAbundanceAnalysis)", {
    kit221MockPlots()
    object <- kit221CompleteObject()
    run <- kit221Da(object, kit221Table(c(kit221C, kit221B), "factor"))
    singles <- kit221Singles(kit221Da, object)
    if (!kit221Completed(run, "W7 factor error")) return(invisible())
    if (!kit221Completed(singles$B_vs_A, "W7 B_vs_A error") || !kit221Completed(singles$C_vs_A, "W7 C_vs_A error")) return(invisible())
    kit221SameRows(run$value$significant_rows, lapply(singles, function(run) run$value$significant_rows), "W7")
    rows <- run$value$significant_rows
    effect <- rows[rows$Protein.Ids %in% paste0("P", 11:16), , drop = FALSE]
    cat(sprintf("\nKIT-221 W7: mean log2FC of P11..P16 under %s\n", paste(names(tapply(effect$log2FC, effect$comparison, mean)),
        sprintf("%.4f", tapply(effect$log2FC, effect$comparison, mean)), sep = " = ", collapse = "; ")))
})

test_that("KIT-221 W8: the sibling wrapper given a factor column returns both comparisons, each equal to its one-row character call", {
    kit221MockPlots()
    object <- kit221CompleteObject()
    run <- kit221Wrapper(object, kit221Table(c(kit221C, kit221B), "factor"))
    if (!kit221Completed(run, "W8 factor error")) return(invisible())
    singles <- kit221Singles(kit221Wrapper, object, "tibble")
    if (!kit221Completed(singles$B_vs_A, "W8 B_vs_A error") || !kit221Completed(singles$C_vs_A, "W8 C_vs_A error")) return(invisible())
    expect_identical(length(run$value$contrasts_results$results), 2L)
    kit221SameRows(run$value$significant_rows, lapply(singles, function(run) run$value$significant_rows), "W8")
})

# ---- W9: anything that is not a data frame is refused with a classed error ----------------------------------------------

kit221Refused <- function(code) {
    caught <- tryCatch(kit221Quiet(code), error = function(e) e)
    expect_true(inherits(caught, "prot_da_contrasts_invalid"),
        label = paste("the call was refused with class prot_da_contrasts_invalid; it gave:",
            if (inherits(caught, "condition")) conditionMessage(caught) else "no error"))
}

test_that("KIT-221 W9: a character vector is refused with the class prot_da_contrasts_invalid by both functions", {
    kit221MockPlots()
    object <- kit221CompleteObject()
    strings <- c(kit221B, kit221C)
    kit221Refused(differentialAbundanceAnalysis(object, contrasts_tbl = strings, formula_string = "~ 0 + group",
        da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, args_row_id = "Protein.Ids"))
    kit221Refused(protein_deAnalysisWrapperFunction(object, contrasts_tbl = strings, formula_string = "~ 0 + group",
        de_q_val_thresh = 0.05, treat_lfc_cutoff = 0, args_row_id = "Protein.Ids"))
})
