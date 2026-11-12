# GUI audit R4, P-05 / A102-01 / A102-02: non-numeric text in a proteomics quantity column, through the real vroom
# reader on small temporary files. The rule (user-approved 2026-11-11), the same for all five readers:
#   - a missing-value token (`PROT_IMPORT_MISSING_TOKENS`: empty, NA, NaN, Filtered, n.q.) is read as NA, and the import
#     signals a notice (class `prot_import_missing_tokens`) with the count per text token (NaN, Filtered, n.q.; an
#     empty cell or NA is read as NA as before, without a notice);
#   - any other text that is not a number refuses the import (class `prot_import_non_numeric_quantity`), naming the
#     column and up to three example values;
#   - a column of plain numbers is read as before.
# Before the repair (Step A, recorded in the R4 report): a decimal comma was read 100 times too large by every reader
# (vroom's "number" guess drops the comma), and a token or an Excel error cell was a silent NA, or for the wide formats
# an opaque pivot error when vroom typed one column as text.

protQuantityRows <- 1000L
# vroom guesses a column's type from 100 rows spread over the file: every 10th row here. Rows 994, 996 and 998 are not
# sampled (the column is typed double and the cell is a parse problem); row 1 is (the column is typed character).
protQuantityLate <- c(994L, 996L, 998L)
protQuantitySampled <- 1L

protQuantityValues <- function() {
    withr::with_seed(20261111L, round(stats::runif(protQuantityRows, 1000, 90000), 2))
}

protQuantityText <- function(kind) {
    text <- format(protQuantityValues(), nsmall = 2, trim = TRUE, scientific = FALSE)
    switch(kind,
        clean = text,
        decimal_comma = sub(".", ",", text, fixed = TRUE),
        filtered_late = replace(text, protQuantityLate, "Filtered"),
        nan_late = replace(text, protQuantityLate, "NaN"),
        nq_late = replace(text, protQuantityLate, "n.q."),
        empty_late = replace(text, protQuantityLate, ""),
        div0_late = replace(text, protQuantityLate, "#DIV/0!"),
        filtered_sampled = replace(text, protQuantitySampled, "Filtered"),
        div0_sampled = replace(text, protQuantitySampled, "#DIV/0!"),
        hex_late = replace(text, protQuantityLate, "0x1A"),
        hex_sampled = replace(text, protQuantitySampled, "0x1A"))
}

protQuantityEdited <- function(kind) {
    if (endsWith(kind, "_late")) protQuantityLate else if (endsWith(kind, "_sampled")) protQuantitySampled else integer()
}

protQuantityWrite <- function(table, ext = ".tsv", env = parent.frame()) {
    path <- withr::local_tempfile(fileext = ext, .local_envir = env)
    utils::write.table(table, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    path
}

# One reader's file with `kind` in its first quantity column (`scope` "one") or in both (`scope` "all", wide formats),
# the reader, the edited column's name as the refusal gives it, and the edited column's values in the result.
protQuantityFixture <- function(reader, kind, scope = "one", env = parent.frame()) {
    n <- protQuantityRows
    ids <- sprintf("P%05d", seq_len(n))
    other <- format(rev(protQuantityValues()), nsmall = 2, trim = TRUE, scientific = FALSE)
    second <- if (identical(scope, "all")) protQuantityText(kind) else other
    first_run <- function(result, run) result$data[[result$column_mapping$quantity_col %||% "Abundance"]][result$data$Run == run]
    switch(reader,
        diann = list(path = protQuantityWrite(data.frame(Protein.Group = ids, Stripped.Sequence = paste0("PEPTIDE", seq_len(n)),
                Run = rep(c("S1", "S2"), n / 2), Precursor.Normalised = protQuantityText(kind), Precursor.Quantity = other,
                Q.Value = "0.001", PG.Q.Value = "0.001", check.names = FALSE), env = env),
            read = importDIANNData, column = "Precursor.Normalised", values = function(result) result$data$Precursor.Normalised),
        spectronaut = list(path = protQuantityWrite(data.frame(PG.ProteinGroups = ids, R.FileName = rep(c("S1", "S2"), n / 2),
                PG.Quantity = protQuantityText(kind), PG.Qvalue = "0.001", check.names = FALSE), env = env),
            read = importSpectronautData, column = "PG.Quantity", values = function(result) result$data$PG.Quantity),
        fragpipe = list(path = protQuantityWrite(data.frame(`Protein ID` = ids, `S1 MaxLFQ Intensity` = protQuantityText(kind),
                `S2 MaxLFQ Intensity` = second, check.names = FALSE), env = env),
            read = importFragPipeData, column = "S1 MaxLFQ Intensity", values = function(result) first_run(result, "S1")),
        maxquant = list(path = protQuantityWrite(data.frame(`Majority protein IDs` = ids, `LFQ intensity S1` = protQuantityText(kind),
                `LFQ intensity S2` = second, `Potential contaminant` = "", Reverse = "", check.names = FALSE), ".txt", env = env),
            read = importMaxQuantData, column = "LFQ intensity S1", values = function(result) first_run(result, "S1")),
        pd_tmt = list(path = protQuantityWrite(data.frame(Accession = ids, `Abundance: F1: 126, Sample_A` = protQuantityText(kind),
                `Abundance: F1: 127N, Sample_B` = second, check.names = FALSE), env = env),
            # The reader checks its abundance columns under the run names it gives them (channel and sample).
            read = importProteomeDiscovererTMTData, column = "126_Sample_A",
            values = function(result) result$data$Abundance[result$data$Run == "126_Sample_A"]))
}

# Read a fixture: the result or the error, the notices signalled and the warnings raised.
protQuantityRead <- function(fixture, read = fixture$read) {
    seen <- list(notices = list(), warnings = character())
    result <- tryCatch(withCallingHandlers(suppressMessages(withCallingHandlers(read(fixture$path),
        prot_import_missing_tokens = function(notice) seen$notices[[length(seen$notices) + 1L]] <<- notice)),
        warning = function(w) {
            seen$warnings <<- c(seen$warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        }), error = function(error) error)
    c(list(result = result), seen)
}

protQuantityReaders <- c("diann", "spectronaut", "fragpipe", "maxquant", "pd_tmt")
protQuantityWide <- c("fragpipe", "maxquant", "pd_tmt")

test_that("P-05: the missing-value tokens are one documented constant", {
    expect_identical(PROT_IMPORT_MISSING_TOKENS, c("", "NA", "NaN", "Filtered", "n.q."))
})

test_that("P-05: a column of plain numbers imports as before, with no notice, for every reader", {
    for (reader in protQuantityReaders) {
        fixture <- protQuantityFixture(reader, "clean")
        read <- protQuantityRead(fixture)
        expect_false(inherits(read$result, "error"), info = reader)
        expect_length(read$notices, 0L)
        expect_length(read$warnings, 0L)
        values <- fixture$values(read$result)
        expect_type(values, "double")
        expect_identical(values, protQuantityValues(), info = reader)
    }
})

test_that("P-05: a decimal-comma quantity column refuses the import, naming the column and example values", {
    for (reader in protQuantityReaders) {
        fixture <- protQuantityFixture(reader, "decimal_comma")
        read <- protQuantityRead(fixture)
        expect_s3_class(read$result, "prot_import_non_numeric_quantity")
        message <- conditionMessage(read$result)
        expect_match(message, sprintf("Quantity column '%s'", fixture$column), fixed = TRUE, info = reader)
        expect_match(message, sprintf("\"%s\"", protQuantityText("decimal_comma")[[1L]]), fixed = TRUE, info = reader)
    }
})

test_that("P-05: a vendor missing-value token in late rows is read as NA with a notice of the count per token", {
    cases <- c(filtered_late = "Filtered", nan_late = "NaN", nq_late = "n.q.")
    for (reader in protQuantityReaders) for (kind in names(cases)) {
        info <- paste(reader, kind)
        fixture <- protQuantityFixture(reader, kind)
        read <- protQuantityRead(fixture)
        expect_false(inherits(read$result, "error"), info = info)
        values <- fixture$values(read$result)
        expect_type(values, "double")
        expect_identical(which(is.na(values)), protQuantityLate, info = info)
        expect_false(any(is.nan(values)), info = info)
        expect_identical(values[-protQuantityLate], protQuantityValues()[-protQuantityLate], info = info)
        expect_length(read$notices, 1L)
        expect_identical(read$notices[[1L]]$counts, stats::setNames(3L, cases[[kind]]), info = info)
        expect_match(conditionMessage(read$notices[[1L]]), sprintf("3 x %s", cases[[kind]]), fixed = TRUE, info = info)
        expect_length(read$warnings, 0L)
    }
})

test_that("P-05: an empty cell is read as NA as before, with no notice (the file format's own missing value)", {
    # The frozen prot_tmt_pd fixture has empty abundance cells: a notice for them would change its entry-lane captures.
    for (reader in protQuantityReaders) {
        fixture <- protQuantityFixture(reader, "empty_late")
        read <- protQuantityRead(fixture)
        expect_false(inherits(read$result, "error"), info = reader)
        expect_identical(which(is.na(fixture$values(read$result))), protQuantityLate, info = reader)
        expect_length(read$notices, 0L)
        expect_length(read$warnings, 0L)
    }
})

test_that("P-05: an Excel error cell refuses the import, in a late row or a sampled one, in one column or in all", {
    for (reader in protQuantityReaders) for (kind in c("div0_late", "div0_sampled")) {
        for (scope in if (reader %in% protQuantityWide) c("one", "all") else "one") {
            info <- paste(reader, kind, scope)
            fixture <- protQuantityFixture(reader, kind, scope)
            read <- protQuantityRead(fixture)
            expect_s3_class(read$result, "prot_import_non_numeric_quantity")
            message <- conditionMessage(read$result)
            expect_match(message, sprintf("Quantity column '%s' contains text that is not a number: \"#DIV/0!\"",
                fixture$column), fixed = TRUE, info = info)
            if (identical(scope, "all")) expect_match(message, "1 other quantity column", fixed = TRUE, info = info)
            expect_no_match(message, "pivot", fixed = TRUE)
        }
    }
})

test_that("P-05: a token in a row vroom samples (a text column) is read as NA with the notice, in one column or in all", {
    for (reader in protQuantityReaders) for (scope in if (reader %in% protQuantityWide) c("one", "all") else "one") {
        info <- paste(reader, scope)
        fixture <- protQuantityFixture(reader, "filtered_sampled", scope)
        read <- protQuantityRead(fixture)
        expect_false(inherits(read$result, "error"), info = info)
        values <- fixture$values(read$result)
        expect_type(values, "double")
        expect_identical(which(is.na(values)), protQuantitySampled, info = info)
        expect_identical(values[-1L], protQuantityValues()[-1L], info = info)
        expect_identical(read$notices[[1L]]$counts, c(Filtered = if (identical(scope, "all")) 2L else 1L), info = info)
    }
})

test_that("P-05: importProteomicsData keeps the refusal's class and message, and passes the notice on", {
    fixture <- protQuantityFixture("fragpipe", "div0_late")
    logged <- character()
    refused <- protQuantityRead(fixture, function(path) importProteomicsData("fragpipe", path,
        logError = function(message) logged <<- c(logged, message)))
    expect_s3_class(refused$result, "prot_import_non_numeric_quantity")
    expect_match(conditionMessage(refused$result),
        "Failed to import data: Quantity column 'S1 MaxLFQ Intensity' contains text that is not a number: \"#DIV/0!\"", fixed = TRUE)
    expect_length(logged, 1L)

    noticed <- protQuantityRead(protQuantityFixture("fragpipe", "nq_late"), function(path) importProteomicsData("fragpipe", path))
    expect_identical(noticed$notices[[1L]]$counts, c(n.q. = 3L))
})

test_that("P-05: the GUI's read shows the notice as a notification and returns the import", {
    fixture <- protQuantityFixture("spectronaut", "filtered_late")
    shown <- list()
    result <- suppressMessages(readProtImportDataWithStatus("spectronaut", fixture$path, input = list(spectronaut_quantity = "pg"),
        updateStatus = function(...) invisible(NULL), logInfo = function(...) invisible(NULL),
        recordImportedData = function(...) invisible(NULL),
        showNotification = function(ui, ..., type = "default") shown[[length(shown) + 1L]] <<- list(ui = ui, type = type)))
    expect_identical(which(is.na(result$data$PG.Quantity)), protQuantityLate)
    expect_length(shown, 1L)
    expect_identical(shown[[1L]]$type, "warning")
    expect_match(shown[[1L]]$ui, "3 x Filtered", fixed = TRUE)
})

test_that("P-05: after the import check no reader leaves a quantity column for the design step to convert", {
    # The managed design step converts the columns `recordDIAImportMetadata` lists as `needs_numeric` with
    # suppressWarnings(as.numeric()) (transformDIADesignTable); a reader's quantity column is now always a plain double.
    for (reader in protQuantityReaders) for (kind in c("clean", "filtered_late", "filtered_sampled")) {
        fixture <- protQuantityFixture(reader, kind, if (reader %in% protQuantityWide) "all" else "one")
        result <- protQuantityRead(fixture)$result
        column <- result$data[[result$column_mapping$quantity_col %||% "Abundance"]]
        expect_true(is.double(column) && is.null(attributes(column)), info = paste(reader, kind))
    }
})

# KIT-115 item 8: R reads "0x1A" as the number 26. A quantity file does not write numbers that way, so it is text that
# is not a number, in a row vroom samples (a text column) as in a late one.
test_that("(KIT-115) hexadecimal text in a quantity column refuses the import, naming the value", {
    for (reader in protQuantityReaders) for (kind in c("hex_late", "hex_sampled")) {
        info <- paste(reader, kind)
        fixture <- protQuantityFixture(reader, kind)
        read <- protQuantityRead(fixture)
        expect_s3_class(read$result, "prot_import_non_numeric_quantity")
        if (inherits(read$result, "error")) expect_match(conditionMessage(read$result), "\"0x1A\"", fixed = TRUE, info = info)
    }
    # The shared rule (lipid, metab and trx sample columns): the same text is not numeric-like.
    expect_false(omicsImportValuesAreNumericLike(c("1.5", "0x1A")))
    expect_false(omicsImportValuesAreNumericLike(c("1.5", " -0X1p3 ")))
    expect_true(omicsImportValuesAreNumericLike(c("1.5", "1e5", "", NA, "0")))
})

# KIT-115 item 9: a DIA-NN parquet report goes through the same check. A quantity column stored as text is read as
# numbers, a missing-value token is counted, and other text refuses the import.
test_that("(KIT-115) a DIA-NN parquet report's quantity column is checked as a text report's is", {
    skip_if_not_installed("arrow")
    n <- protQuantityRows
    report <- function(kind) {
        path <- withr::local_tempfile(fileext = ".parquet", .local_envir = parent.frame())
        arrow::write_parquet(data.frame(Protein.Group = sprintf("P%05d", seq_len(n)),
            Stripped.Sequence = paste0("PEPTIDE", seq_len(n)), Run = rep(c("S1", "S2"), n / 2),
            Precursor.Normalised = protQuantityText(kind), Precursor.Quantity = rev(protQuantityValues()),
            Q.Value = 0.001, PG.Q.Value = 0.001), path)
        protQuantityRead(list(path = path, read = importDIANNData))
    }
    clean <- report("clean")
    expect_false(inherits(clean$result, "error"))
    expect_identical(clean$result$data$Precursor.Normalised, protQuantityValues())
    expect_length(clean$notices, 0L)
    token <- report("filtered_late")
    expect_false(inherits(token$result, "error"))
    expect_identical(which(is.na(token$result$data$Precursor.Normalised)), protQuantityLate)
    expect_identical(token$notices[[1L]]$counts, c(Filtered = 3L))
    expect_s3_class(report("decimal_comma")$result, "prot_import_non_numeric_quantity")
    expect_s3_class(report("div0_late")$result, "prot_import_non_numeric_quantity")
})


# KIT-118 (user, 2026-11-12, after the A101-07 adjudication and its audit): a DIA-NN quantity of exactly 0 is a precursor
# DIA-NN did not quantify, a missing value (as for LFQ and TMT, KIT-101 F-18/F-26). Before, a 0 was counted as observed
# and entered the technical-replicate imputation's mean. The rule covers both quantity columns whatever the
# selection: the analysis reads both (IQ from Precursor.Quantity; the filters and limpa from Precursor.Normalised). A negative or infinite quantity refuses the import.
protDiannZeroFile <- function(quantity, normalised, env = parent.frame()) {
    n <- length(quantity)
    protQuantityWrite(data.frame(Protein.Group = sprintf("P%05d", seq_len(n)),
        Stripped.Sequence = paste0("PEPTIDE", seq_len(n)), Run = rep(c("S1", "S2"), length.out = n),
        Precursor.Normalised = normalised, Precursor.Quantity = quantity, Q.Value = 0.001, PG.Q.Value = 0.001), env = env)
}

test_that("(KIT-118) DIA-NN's zero quantities are read as missing in the columns the workflow reads, with a notice", {
    n <- protQuantityRows
    zero <- c(5L, 500L)
    values <- protQuantityValues()
    other <- rev(values)
    path <- protDiannZeroFile(replace(other, zero, 0), replace(values, zero, 0))
    kept <- setdiff(seq_len(n), zero)
    selected <- protQuantityRead(list(path = path), read = function(p) importDIANNData(p, use_precursor_norm = TRUE))
    expect_false(inherits(selected$result, "error"))
    expect_identical(which(is.na(selected$result$data$Precursor.Normalised)), zero)
    expect_identical(which(is.na(selected$result$data$Precursor.Quantity)), zero)
    expect_identical(selected$result$data$Precursor.Normalised[kept], values[kept])
    expect_identical(selected$result$data$Precursor.Quantity[kept], other[kept])
    expect_length(selected$notices, 1L)
    expect_identical(selected$notices[[1L]]$counts, c(Precursor.Normalised = 2L, Precursor.Quantity = 2L))
    # Precursor.Quantity selected: the same, because the analysis reads both columns whatever the selection.
    raw <- protQuantityRead(list(path = path), read = function(p) importDIANNData(p, use_precursor_norm = FALSE))
    expect_identical(which(is.na(raw$result$data$Precursor.Quantity)), zero)
    expect_identical(which(is.na(raw$result$data$Precursor.Normalised)), zero)
    expect_identical(raw$notices[[1L]]$counts, c(Precursor.Normalised = 2L, Precursor.Quantity = 2L))
    # No zero, no notice: the file imports as before.
    clean <- protQuantityRead(list(path = protDiannZeroFile(other, values)), read = importDIANNData)
    expect_length(clean$notices, 0L)
    expect_identical(clean$result$data$Precursor.Quantity, other)
})

test_that("(KIT-118) a text Precursor.Quantity is checked and its zeros read as missing when Precursor.Normalised is selected", {
    # Before, P-05 checked only the selected column: a text Precursor.Quantity stayed text and its "0" survived.
    quantity <- c("10", "0", "Filtered", "12.5", "", "30")
    read <- protQuantityRead(list(path = protDiannZeroFile(quantity, c(1, 2, 3, 4, 5, 6))),
        read = function(p) importDIANNData(p, use_precursor_norm = TRUE))
    expect_false(inherits(read$result, "error"))
    expect_identical(read$result$data$Precursor.Quantity, c(10, NA, NA, 12.5, NA, 30))
    counts <- lapply(read$notices, `[[`, "counts")
    expect_identical(counts, list(c(Filtered = 1L), c(Precursor.Quantity = 1L)))
    refused <- protQuantityRead(list(path = protDiannZeroFile(c("10", "abc", "1", "2", "3", "4"), 1:6)),
        read = function(p) importDIANNData(p, use_precursor_norm = TRUE))
    expect_s3_class(refused$result, "prot_import_non_numeric_quantity")
})

test_that("(KIT-118) a negative or infinite DIA-NN quantity refuses the import, naming the column", {
    for (bad in c(-1, Inf)) {
        read <- protQuantityRead(list(path = protDiannZeroFile(c(10, bad, 3, 4), c(1, 2, 3, 4))),
            read = function(p) importDIANNData(p, use_precursor_norm = TRUE))
        expect_s3_class(read$result, "error")
        expect_match(conditionMessage(read$result), "Quantity column 'Precursor.Quantity' contains negative or infinite value\\(s\\) \\(e.g. .*, row 2\\)")
    }
})

test_that("(KIT-118) a DIA-NN zero beside a missing technical replicate no longer lowers the imputed value", {
    # Import, precursor rollup and technical-replicate imputation on Precursor.Quantity, as the IQ route runs them
    # (raw_quantity_column, func_pept_s4_missingness.R). Peptide PEPA in four technical replicates: 100, 200, 0 and
    # missing. Before, the 0 was observed: it stayed 0 and the missing run was imputed with mean(100, 200, 0) = 100.
    runs <- c("R1", "R2", "R3", "R4")
    report <- data.frame(Protein.Group = "P1", Protein.Ids = "P1", Stripped.Sequence = rep(c("PEPA", "PEPB"), each = 4),
        Modified.Sequence = rep(c("PEPA", "PEPB"), each = 4), Precursor.Id = rep(c("PEPA2", "PEPB2"), each = 4),
        Precursor.Charge = 2L, Run = rep(runs, 2), Precursor.Quantity = c(100, 200, 0, NA, 10, 20, 30, 40),
        Precursor.Normalised = c(100, 200, 0, NA, 10, 20, 30, 40), Q.Value = 0.001, PG.Q.Value = 0.001)
    imported <- suppressMessages(importDIANNData(protQuantityWrite(report), use_precursor_norm = TRUE))$data
    peptides <- suppressMessages(rollUpPrecursorToPeptideHelper(as.data.frame(imported), core_utilisation = NA))
    imputed <- suppressMessages(peptideMissingValueImputationHelper(as.data.frame(peptides),
        data.frame(Run = runs, technical_group = "G1"), Run, Run, technical_group, Protein.Ids, Stripped.Sequence,
        Peptide.RawQuantity, Peptide.Imputed, core_utilisation = NA))
    pepa <- imputed[imputed$Stripped.Sequence == "PEPA", ]
    pepa <- pepa[match(runs, pepa$Run), ]
    expect_identical(pepa$Peptide.Imputed, c(100, 200, 150, 150))
    expect_identical(pepa$is_imputed, c(FALSE, FALSE, TRUE, TRUE))
    pepb <- imputed[imputed$Stripped.Sequence == "PEPB", ]
    expect_identical(pepb$Peptide.Imputed[match(runs, pepb$Run)], c(10, 20, 30, 40))
})

test_that("(KIT-118) every zero in the synthetic DIA journey corpus is read as missing at import, and nothing else changes", {
    # The DIA journeys import this generator's report (tools/profiling/workflow_perf107_synthetic.R writes ~2% zeros
    # into Precursor.Quantity); here its small calibration cohort, through the reader with Precursor.Normalised selected.
    code <- new.env(parent = asNamespace("MultiScholaR"))
    sys.source(file.path(normalizePath(test_path("..", "..")), "tools/profiling/workflow_perf107_synthetic.R"), code)
    recipe <- code$perf107SyntheticRecipe("calibration", small = TRUE)
    catalogue <- code$perf107SyntheticCatalogue(recipe)
    report <- do.call(rbind, lapply(seq_len(recipe$samples), function(i) code$perf107SyntheticSample(recipe, catalogue, i)))
    zero <- !is.na(report$Precursor.Quantity) & report$Precursor.Quantity == 0
    expect_gt(sum(zero), 0L)
    read <- protQuantityRead(list(path = protQuantityWrite(report)), read = importDIANNData)
    expect_false(inherits(read$result, "error"))
    data <- read$result$data
    expect_identical(sum(is.na(data$Precursor.Quantity)), sum(is.na(report$Precursor.Quantity)) + sum(zero))
    # write.table keeps 15 significant digits: the other values equal the generator's to that precision.
    expect_equal(data$Precursor.Quantity[!zero], report$Precursor.Quantity[!zero], tolerance = 1e-14)
    expect_identical(read$notices[[1L]]$counts, c(Precursor.Quantity = sum(zero)))
})
