# KIT-218 (VAL-180, VAL-188): the exported scripted wrapper daAnalysisWrapperFunction passed its named list of
# per-contrast tables to getSignificantData and printCountDaGenesTable one list level too deep
# (R/func_prot_da_model_wrapper.R:273-276, :306-309), so every call stopped after the DA step in rownames_to_column
# ("is.data.frame(df) is not TRUE"). Each W test first asserts that its call completed (`expect_null` on the caught
# message, `conditionMessage(e)`: purrr wraps the error twice), so the tests-only commit fails by assertion at that
# message. G1 (base only) undoes the nesting with test doubles that call the package's own functions, to show that
# nothing later in the tail stops first. Synthetic data only.
#
# The fixture builders are copied, renamed, from tests/testthat/test-kit200-prot-da-no-value-proteins.R at
# 817b026a (unchanged at 3f8203c2): kit200Scatter :10, kit200Quiet :15-19, kit200WithSeed :21-22,
# kit200Matrix :25-35, kit200Design :37-41, kit200SetScatter :43-46, kit200AddRow :48-52, kit200Object :109-116.

kit218Scatter <- list(c("P2", "S1"), c("P3", "S5"), c("P3", "S9"), c("P4", "S12"))

kit218Quiet <- function(expr) {
    value <- NULL
    utils::capture.output(value <- suppressWarnings(suppressMessages(expr)))
    value
}

kit218WithSeed <- function(seed, code) withr::with_seed(seed, code, .rng_kind = "Mersenne-Twister",
    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

# Proteins P1..P60, levels N(22, 2^2), residual sd 0.3, P11..P16 with a group-B effect of +1 (KIT-218 section 4).
kit218Matrix <- function(sizes = c(A = 4L, B = 4L, C = 4L), seed = 200001L, proteins = 60L) {
    groups <- rep(names(sizes), sizes)
    x <- kit218WithSeed(seed, {
        level <- stats::rnorm(proteins, 22, 2)
        effect <- as.numeric(seq_len(proteins) %in% 11:16)
        level + outer(effect, as.numeric(groups == "B")) +
            matrix(stats::rnorm(proteins * length(groups), sd = 0.3), proteins)
    })
    dimnames(x) <- list(paste0("P", seq_len(proteins)), paste0("S", seq_along(groups)))
    x
}

kit218Design <- function(sizes = c(A = 4L, B = 4L, C = 4L)) {
    runs <- paste0("S", seq_len(sum(sizes)))
    data.frame(Run = runs, group = rep(names(sizes), sizes), replicates = runs, row.names = runs,
        stringsAsFactors = FALSE)
}

kit218SetScatter <- function(x, cells = kit218Scatter) {
    for (cell in cells) x[cell[[1L]], cell[[2L]]] <- NA_real_
    x
}

kit218AddRow <- function(x, id, values = numeric()) {
    row <- matrix(NA_real_, 1L, ncol(x), dimnames = list(id, colnames(x)))
    row[1L, names(values)] <- values
    rbind(x, row)
}

kit218Object <- function(x, design) {
    rownames(design) <- NULL
    ProteinQuantitativeData(protein_quant_table = data.frame(Protein.Ids = rownames(x), x, row.names = NULL,
        check.names = FALSE, stringsAsFactors = FALSE), design_matrix = design, sample_id = "Run", group_id = "group",
        technical_replicate_id = "replicates", protein_id_column = "Protein.Ids",
        protein_id_table = data.frame(Protein.Ids = rownames(x), gene_names = paste0("GENE_", rownames(x)),
            stringsAsFactors = FALSE), args = list())
}

# ---- KIT-218's own helpers ------------------------------------------------------------------------------------------

kit218TwoContrasts <- function() data.frame(contrasts = c("groupB-groupA", "groupC-groupA"),
    full_format = c("B_vs_A=groupB-groupA", "C_vs_A=groupC-groupA"), stringsAsFactors = FALSE)

kit218OneContrast <- function() data.frame(contrasts = "groupB-groupA", full_format = "B_vs_A=groupB-groupA",
    stringsAsFactors = FALSE)

# W4's input (OD-4 (a)): a base data.frame without full_format whose first column holds name=expression strings.
kit218NoFullFormat <- function() data.frame(contrasts = c("B_vs_A=groupB-groupA", "C_vs_A=groupC-groupA"),
    stringsAsFactors = FALSE)

kit218ReturnedNames <- c("significant_rows", "volplot_plot", "num_sig_da_molecules_first_go", "pvalhist",
    "norm_counts", "da_proteins_wide", "da_proteins_long", "list_of_volcano_plots", "num_sig_da_molecules")

# The wrapper's RLE, PCA and value-count plots precede the DA step and are not this ticket's subject (as in KIT-200 T9).
kit218MockPlots <- function(env = parent.frame()) {
    testthat::local_mocked_bindings(plotRle = function(...) ggplot2::ggplot(),
        plotPca = function(...) ggplot2::ggplot(), plotNumOfValuesNoLog = function(...) ggplot2::ggplot(),
        .package = "MultiScholaR", .env = env)
}

# G1's doubles: the package's own functions, captured before mocking, given the list of tables once (the fix's
# argument, list(contrasts_results_table)); every other argument is forwarded unchanged through `...`.
kit218MockUnnested <- function(env = parent.frame()) {
    significant <- get("getSignificantData", envir = asNamespace("MultiScholaR"))
    count <- get("printCountDaGenesTable", envir = asNamespace("MultiScholaR"))
    testthat::local_mocked_bindings(
        getSignificantData = function(list_of_da_tables, ...) significant(
            list_of_da_tables = list(list_of_da_tables[[1L]][[1L]]), ...),
        printCountDaGenesTable = function(list_of_da_tables, ...) count(
            list_of_da_tables = list(list_of_da_tables[[1L]][[1L]]), ...),
        .package = "MultiScholaR", .env = env)
}

kit218Call <- function(object, contrasts_tbl) {
    tryCatch(list(value = kit218Quiet(daAnalysisWrapperFunction(object, contrasts_tbl = contrasts_tbl,
        formula_string = "~ 0 + group", treat_lfc_cutoff = 0, eBayes_trend = TRUE, eBayes_robust = TRUE,
        args_row_id = "Protein.Ids"))), error = function(e) list(error = conditionMessage(e)))
}

# The base-side guard: the call completed. On the tests-only commit this fails and prints the caught message.
kit218Completed <- function(run, label) {
    expect_null(run$error, label = label)
    is.null(run$error)
}

# W1's value check: for each contrast and protein, log2FC, fdr_qvalue and raw_pvalue in significant_rows are identical
# to the DA step's table (logFC, fdr_qvalue, raw_pvalue), matched by row name.
kit218ExpectCarried <- function(value, ids = NULL) {
    rows <- value$significant_rows
    for (full in names(value$contrasts_results$results)) {
        ref <- value$contrasts_results$results[[full]]
        keep <- if (is.null(ids)) rownames(ref) else ids
        short <- strsplit(full, "=", fixed = TRUE)[[1L]][[1L]]
        got <- rows[rows$comparison == short, , drop = FALSE]
        expect_identical(sort(got$Protein.Ids), sort(rownames(ref)), info = full)
        at <- match(keep, got$Protein.Ids)
        expect_identical(unname(got$log2FC[at]), unname(ref[keep, "logFC"]), info = paste(full, "log2FC"))
        expect_identical(unname(got$fdr_qvalue[at]), unname(ref[keep, "fdr_qvalue"]), info = paste(full, "fdr_qvalue"))
        expect_identical(unname(got$raw_pvalue[at]), unname(ref[keep, "raw_pvalue"]), info = paste(full, "raw_pvalue"))
    }
}

kit218CompleteObject <- function() kit218Object(kit218Matrix(), kit218Design())

kit218T9Object <- function() kit218Object(kit218AddRow(kit218SetScatter(kit218Matrix()), "P_empty"), kit218Design())

# ---- W1 to W4 (failing first on base) --------------------------------------------------------------------------------

test_that("KIT-218 W1: two contrasts on complete data run through the wrapper's tail with the DA step's numbers", {
    kit218MockPlots()
    run <- kit218Call(kit218CompleteObject(), kit218TwoContrasts())
    if (!kit218Completed(run, "W1 error")) return(invisible())
    value <- run$value
    expect_true(all(kit218ReturnedNames %in% names(value)))
    rows <- value$significant_rows
    expect_identical(nrow(rows), 120L)
    expect_identical(sort(unique(rows$comparison)), c("B_vs_A", "C_vs_A"))
    expect_identical(sort(unique(rows$log_intensity)), c("groupB-groupA", "groupC-groupA"))
    kit218ExpectCarried(value)
    expect_setequal(unique(value$num_sig_da_molecules_first_go$table$comparison), c("B_vs_A", "C_vs_A"))
    expect_identical(nrow(value$da_proteins_long), 120L)
    expect_identical(nrow(value$da_proteins_wide), 60L)
    expect_identical(anyDuplicated(value$da_proteins_wide$Protein.Ids), 0L)
})

test_that("KIT-218 W2: KIT-200's T9 object (a protein without values) runs through the wrapper's tail", {
    kit218MockPlots()
    run <- kit218Call(kit218T9Object(), kit218OneContrast())
    if (!kit218Completed(run, "W2 error")) return(invisible())
    value <- run$value
    rows <- value$significant_rows
    empty <- rows[rows$Protein.Ids == "P_empty", , drop = FALSE]
    expect_identical(nrow(empty), 1L)
    expect_true(is.na(empty$log2FC))
    expect_true(is.na(empty$fdr_qvalue))
    expect_true(is.na(empty$raw_pvalue))
    expect_identical(as.character(empty$colour), "black")
    expect_identical(sum(value$da_proteins_long$Protein.Ids == "P_empty"), 1L)
    kit218ExpectCarried(value, ids = paste0("P", 1:60))
    # KIT-204 C14: printed, not asserted.
    cat("\nKIT-218 W2: num_sig_da_molecules (P_empty's status is counted here):\n")
    print(as.data.frame(value$num_sig_da_molecules))
})

test_that("KIT-218 W3: one contrast on complete data runs through the wrapper's tail", {
    kit218MockPlots()
    run <- kit218Call(kit218CompleteObject(), kit218OneContrast())
    if (!kit218Completed(run, "W3 error")) return(invisible())
    value <- run$value
    expect_identical(nrow(value$significant_rows), 60L)
    expect_identical(unique(value$significant_rows$comparison), "B_vs_A")
    kit218ExpectCarried(value)
})

test_that("KIT-218 W4: a data.frame of contrasts without full_format gives every row's comparison (OD-4 (a))", {
    kit218MockPlots()
    run <- kit218Call(kit218CompleteObject(), kit218NoFullFormat())
    if (!kit218Completed(run, "W4 error")) return(invisible())
    rows <- run$value$significant_rows
    expect_identical(sort(unique(rows$comparison)), c("B_vs_A", "C_vs_A"))
    expect_identical(nrow(rows), 120L)
})

# ---- G1: gate, base only ----------------------------------------------------------------------------------------------

test_that("KIT-218 G1: on base, with the nesting undone by test doubles, the tail runs to its end on W1 to W4's inputs", {
    skip_if(!any(grepl("nested_list", deparse(daAnalysisWrapperFunction))), "the fix is present (G1 is base only)")
    kit218MockPlots()
    kit218MockUnnested()
    cases <- list(
        W1 = list(object = kit218CompleteObject(), contrasts = kit218TwoContrasts(), rows = 120L),
        W2 = list(object = kit218T9Object(), contrasts = kit218OneContrast(), rows = 61L),
        W3 = list(object = kit218CompleteObject(), contrasts = kit218OneContrast(), rows = 60L),
        W4 = list(object = kit218CompleteObject(), contrasts = kit218NoFullFormat(), rows = NA_integer_))
    for (case in names(cases)) {
        run <- kit218Call(cases[[case]]$object, cases[[case]]$contrasts)
        if (!kit218Completed(run, paste("G1", case, "error"))) next
        rows <- run$value$significant_rows
        cat(sprintf("\nKIT-218 G1 %s: completed; significant_rows %d rows; comparisons: %s\n", case, nrow(rows),
            paste(sort(unique(rows$comparison)), collapse = ", ")))
        if (!is.na(cases[[case]]$rows)) expect_identical(nrow(rows), cases[[case]]$rows, info = case)
        expect_true(all(kit218ReturnedNames %in% names(run$value)), info = case)
    }
})

# ---- G3: guard -------------------------------------------------------------------------------------------------------

test_that("KIT-218 G3: input without any value is refused with KIT-196's class before the tail", {
    kit218MockPlots()
    x <- kit218Matrix()
    x[] <- NA_real_
    expect_error(kit218Quiet(daAnalysisWrapperFunction(kit218Object(x, kit218Design()),
        contrasts_tbl = kit218TwoContrasts(), formula_string = "~ 0 + group", treat_lfc_cutoff = 0,
        eBayes_trend = TRUE, eBayes_robust = TRUE, args_row_id = "Protein.Ids")),
        class = "omics_da_no_observed_features")
})

# ---- M1: measurement, recorded not judged ----------------------------------------------------------------------------

test_that("KIT-218 M1 (recorded, not judged): groups named 1, 2, 3 (F1)", {
    kit218MockPlots()
    sizes <- c(`1` = 4L, `2` = 4L, `3` = 4L)
    object <- kit218Object(kit218Matrix(sizes), kit218Design(sizes))
    run <- kit218Call(object, data.frame(contrasts = c("group2-group1", "group3-group1"),
        full_format = c("2_vs_1=group2-group1", "3_vs_1=group3-group1"), stringsAsFactors = FALSE))
    cat(sprintf("\nKIT-218 M1: %s\n", if (is.null(run$error)) "the call completed" else
        paste("caught message:", run$error)))
    succeed()
})

# ---- H1 (ticket section 14a): the comparisons are read as text; a non-data-frame is refused -------------------------

kit218FactorRows <- function(factors) data.frame(contrasts = c("groupC-groupA", "groupB-groupA"),
    full_format = c("C_vs_A=groupC-groupA", "B_vs_A=groupB-groupA"), stringsAsFactors = factors)

kit218FactorNoFullFormat <- function(factors) data.frame(contrasts = c("C_vs_A=groupC-groupA", "B_vs_A=groupB-groupA"),
    stringsAsFactors = factors)

# Per comparison and protein, log2FC, fdr_qvalue and raw_pvalue of `got` identical to the reference call's.
kit218ExpectSameRows <- function(got, ref) {
    expect_identical(sort(unique(got$comparison)), sort(unique(ref$comparison)))
    for (comparison in unique(ref$comparison)) {
        a <- got[got$comparison == comparison, , drop = FALSE]
        b <- ref[ref$comparison == comparison, , drop = FALSE]
        expect_identical(sort(a$Protein.Ids), sort(b$Protein.Ids), info = comparison)
        at <- match(b$Protein.Ids, a$Protein.Ids)
        for (column in c("log2FC", "fdr_qvalue", "raw_pvalue")) {
            expect_identical(a[[column]][at], b[[column]], info = paste(comparison, column))
        }
    }
}

test_that("KIT-218 W5: a contrasts_tbl that is not a data frame is refused with a classed error", {
    kit218MockPlots()
    expect_error(kit218Quiet(daAnalysisWrapperFunction(kit218CompleteObject(),
        contrasts_tbl = c("B_vs_A=groupB-groupA", "C_vs_A=groupC-groupA"), formula_string = "~ 0 + group",
        treat_lfc_cutoff = 0, eBayes_trend = TRUE, eBayes_robust = TRUE, args_row_id = "Protein.Ids")),
        class = "prot_da_contrasts_invalid")
})

test_that("KIT-218 W6: a full_format column stored as a factor gives each comparison its own statistics", {
    kit218MockPlots()
    run <- kit218Call(kit218CompleteObject(), kit218FactorRows(TRUE))
    if (!kit218Completed(run, "W6 error")) return(invisible())
    text <- kit218Call(kit218CompleteObject(), kit218FactorRows(FALSE))
    if (!kit218Completed(text, "W6 text reference error")) return(invisible())
    kit218ExpectSameRows(run$value$significant_rows, text$value$significant_rows)
})

test_that("KIT-218 W7: a first column stored as a factor, without full_format, gives both comparisons", {
    kit218MockPlots()
    run <- kit218Call(kit218CompleteObject(), kit218FactorNoFullFormat(TRUE))
    if (!kit218Completed(run, "W7 error")) return(invisible())
    expect_identical(sort(unique(run$value$significant_rows$comparison)), c("B_vs_A", "C_vs_A"))
    text <- kit218Call(kit218CompleteObject(), kit218FactorNoFullFormat(FALSE))
    if (!kit218Completed(text, "W7 text reference error")) return(invisible())
    kit218ExpectSameRows(run$value$significant_rows, text$value$significant_rows)
})

test_that("KIT-218 M2 (recorded, not judged; base only): W6's factor input through G1's doubles", {
    skip_if(!any(grepl("nested_list", deparse(daAnalysisWrapperFunction))), "the fix is present (M2 is base only)")
    kit218MockPlots()
    kit218MockUnnested()
    run <- kit218Call(kit218CompleteObject(), kit218FactorRows(TRUE))
    if (is.null(run$error)) {
        rows <- run$value$significant_rows
        effect <- rows[rows$Protein.Ids %in% paste0("P", 11:16), , drop = FALSE]
        means <- tapply(effect$log2FC, effect$comparison, mean)
        cat(sprintf("\nKIT-218 M2: mean log2FC of P11..P16 under %s\n", paste(names(means), sprintf("%.4f", means),
            sep = " = ", collapse = "; ")))
    } else {
        cat(sprintf("\nKIT-218 M2: caught message: %s\n", run$error))
    }
    succeed()
})
