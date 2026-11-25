# KIT-235: the lipidomics and metabolomics report templates read the study-parameters file the writer produces ('*'
# bullets, R/func_general_filemgmt_study_parameter_writer.R). Each case writes its file with createWorkflowArgsFromConfig
# from the writer-contract fixtures (tests/testthat/test-general-filemgmt-s4-contracts.R, copied), then runs the
# template's chunks (parser, Format*, Sample Analysis, Main with its own handler) in an environment created as
# RenderReport creates it, new.env(parent = asNamespace("MultiScholaR")), with the Project Setup chunk's paths, %||%
# and a d66_log stub (ticket section 3, G0a). T6 compares the head text with the base templates' chunks
# (tests/testdata/kit235/base, verbatim at the branch point) on the older '\u2022' and '-' shapes.

k235Templates <- list(
    lipid = list(dir = "lipidomics", file = "lipidomics_report.rmd", extract = "ExtractLipidStudyParameters"),
    metab = list(dir = "metabolomics", file = "metabolomics_report.rmd", extract = "ExtractStudyParametersNew"))

k235Chunks <- c("Hierarchical Parser", "Format Sample Details", "Method Section", "Result Section", "Comment Section",
    "Opinions Section", "Differential Abundance Section", "Pathway Enrichment Section", "Sample Analysis Functions",
    "Main Execution Function")

#' A template file: the package's (head) or the base copy of its chunks
k235Path <- function(omic, side) {
    spec <- k235Templates[[omic]]
    if (identical(side, "head")) system.file("reports", spec$dir, spec$file, package = "MultiScholaR") else
        testthat::test_path("..", "testdata", "kit235", "base", spec$file)
}

#' One chunk of a template, extracted with knitr::purl and parsed
k235Exprs <- function(path, label) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    start <- grep(paste0("^```\\{r ", label, "[,}]"), lines)
    stopifnot(length(start) == 1L)
    end <- start + which(lines[-seq_len(start)] == "```")[[1L]]
    chunk <- tempfile(fileext = ".Rmd")
    writeLines(lines[start:end], chunk, useBytes = TRUE)
    code <- tempfile(fileext = ".R")
    suppressMessages(knitr::purl(chunk, output = code, quiet = TRUE, documentation = 0))
    parse(code, keep.source = FALSE, encoding = "UTF-8")
}

# --- The writer's files ----------------------------------------------------------------------------------------------

# Fixtures of tests/testthat/test-general-filemgmt-s4-contracts.R (makeMetabolomicsWorkflowObject,
# makeLipidomicsWorkflowObject), copied; `ruv` selects the RUV arguments (M, L1, L2, S of ticket section 3).
k235Metab <- function(ruv = TRUE) {
    design_matrix <- data.frame(Sample_ID = c("S1", "S2", "S3"), group = c("A", "A", "B"), stringsAsFactors = FALSE)
    assay_pos <- data.frame(database_identifier = c("M1", "M2"), metabolite_identification = c("Met1", "Met2"),
        S1 = c(10, 20), S2 = c(11, 21), S3 = c(12, 22), check.names = FALSE)
    assay_neg <- data.frame(database_identifier = c("M3", "M4", "M5"), metabolite_identification = c("Met3", "Met4", "Met5"),
        S1 = c(30, 40, 50), S2 = c(31, 41, 51), S3 = c(32, 42, 52), check.names = FALSE)
    args <- list(
        ITSDNormalization = list(applied = TRUE, method_type = "median", itsd_aggregation = "sum",
            itsd_pattern_columns = c("ISTD_1", "ISTD_2"), removed_itsd = TRUE,
            timestamp = as.POSIXct("2026-04-18 09:00:00", tz = "UTC"),
            itsd_features_per_assay = list(assay_pos = c("ISTD_A1", "ISTD_A2"), assay_neg = c("ISTD_B1")),
            itsd_counts_per_assay = list(assay_pos = 2L, assay_neg = 1L)),
        log_transformed = TRUE, log_transform_offset = 1, normalisation_method = "PQN")
    if (ruv) args <- c(args, list(ruv_grouping_variable = "group", ruv_number_k = list(assay_pos = 2L, assay_neg = 3L),
        ctrl = list(assay_pos = c(TRUE, FALSE), assay_neg = c(TRUE, TRUE, FALSE))))
    createMetaboliteAssayData(metabolite_data = list(assay_pos = assay_pos, assay_neg = assay_neg),
        design_matrix = design_matrix, metabolite_id_column = "database_identifier",
        annotation_id_column = "metabolite_identification", sample_id = "Sample_ID", group_id = "group",
        database_identifier_type = "HMDB", internal_standard_regex = "^ISTD", args = args)
}

k235Lipid <- function(ruv = c("single", "per_assay", "none")) {
    ruv <- match.arg(ruv)
    design_matrix <- data.frame(Sample_ID = c("S1", "S2", "S3"), group = c("A", "A", "B"), stringsAsFactors = FALSE)
    assay_data <- data.frame(database_identifier = c("L1", "L2", "L3"), lipid_identification = c("Lipid1", "Lipid2", "Lipid3"),
        S1 = c(100, 200, 300), S2 = c(101, 201, 301), S3 = c(102, 202, 302), check.names = FALSE)
    args <- list(ITSDNormalization = list(applied = FALSE), log_transformed = FALSE, normalisation_method = "median")
    if (identical(ruv, "single")) args <- c(args, list(ruv_grouping_variable = "group", ruv_number_k = 4L,
        ctrl = c(TRUE, FALSE, TRUE)))
    if (identical(ruv, "per_assay")) args <- c(args, list(ruv_grouping_variable = "group",
        ruv_number_k = list(assay_main = 4L), ctrl = list(assay_main = c(TRUE, FALSE, TRUE))))
    createLipidomicsAssayData(lipid_data = list(assay_main = assay_data), design_matrix = design_matrix,
        lipid_id_column = "database_identifier", annotation_id_column = "lipid_identification", sample_id = "Sample_ID",
        group_id = "group", database_identifier_type = "LipidMaps", internal_standard_regex = "^LIPID_ISTD", args = args)
}

k235Cases <- list(
    M = list(omic = "metab", workflow = "metabolomics", make = function() k235Metab(TRUE)),
    L1 = list(omic = "lipid", workflow = "lipidomics", make = function() k235Lipid("single")),
    L2 = list(omic = "lipid", workflow = "lipidomics", make = function() k235Lipid("per_assay")),
    S_lipid = list(omic = "lipid", workflow = "lipidomics", make = function() k235Lipid("none")),
    S_metab = list(omic = "metab", workflow = "metabolomics", make = function() k235Metab(FALSE)))

#' The study-parameters lines the writer produces for a case (its '*' shape)
k235WriterLines <- function(case) {
    spec <- k235Cases[[case]]
    dir <- withr::local_tempdir(pattern = "kit235-writer-")
    had_config <- exists("config_list", envir = .GlobalEnv, inherits = FALSE)
    old_config <- if (had_config) get("config_list", envir = .GlobalEnv)
    assign("config_list", list(), envir = .GlobalEnv)
    on.exit(if (had_config) assign("config_list", old_config, envir = .GlobalEnv) else rm("config_list", envir = .GlobalEnv),
        add = TRUE)
    path <- NULL
    invisible(utils::capture.output(path <- createWorkflowArgsFromConfig(workflow_name = spec$workflow,
        source_dir_path = dir, final_s4_object = spec$make(),
        contrasts_tbl = data.frame(contrast = "A_vs_B", stringsAsFactors = FALSE)), type = "output"))
    readLines(path, warn = FALSE, encoding = "UTF-8")
}

#' The same file in an older bullet shape: each leading '* ' becomes '\u2022 ' or '- '
k235Shape <- function(lines, shape = c("star", "bullet", "dash")) {
    switch(match.arg(shape), star = lines, bullet = sub("^\\* ", "\u2022 ", lines), dash = sub("^\\* ", "- ", lines))
}

#' An encoded step parameter, as the record stores it (`encodeWorkflowRecordValue`)
k235Step <- function(operation, parameters = list(), outcomes = list()) {
    encode <- function(name, value) list(name = name, value = list(type = typeof(value), value = as.list(value),
        length = length(value)))
    list(operation = operation, parameters = unname(Map(encode, names(parameters), parameters)), outcomes = outcomes)
}

#' Run one file through a template's chunks: the report content Main's handler leaves (`results_data`), the parsed
#' file, the extractor called alone and AnalyzeSamples' numbers.
k235Run <- function(omic, lines, side = "head", record = NULL) {
    spec <- k235Templates[[omic]]
    dir <- withr::local_tempdir(pattern = "kit235-")
    writeLines(enc2utf8(lines), file.path(dir, "study_parameters.txt"), useBytes = TRUE)
    writeLines(c("Run\tgroup", "S1\tA", "S2\tA", "S3\tB"), file.path(dir, "design_matrix.tab"))
    if (!is.null(record)) jsonlite::write_json(list(steps = record), file.path(dir, "analysis_record.json"),
        auto_unbox = TRUE, null = "null", digits = NA)
    env <- new.env(parent = asNamespace("MultiScholaR"))
    env$report_source_dir <- dir
    env$report_publication_tables_dir <- file.path(dir, "publication_tables")
    env$params <- list(omic_type = spec$dir, experiment_label = "kit235", workflow_name = "kit235")
    env$`%||%` <- function(a, b) if (!is.null(a)) a else b
    env$d66_log <- function(...) invisible(NULL)
    path <- k235Path(omic, side)
    for (label in k235Chunks) for (e in k235Exprs(path, label)) suppressWarnings(eval(e, env))
    study <- env$ReadStudyParamsNew(file.path(dir, "study_parameters.txt"))
    design <- readr::read_tsv(file.path(dir, "design_matrix.tab"), show_col_types = FALSE)
    list(content = env$results_data$content, extracted = env[[spec$extract]](study),
        numbers = suppressWarnings(env$AnalyzeSamples(design, study))$numbers)
}

k235Text <- function(x) paste(as.character(x), collapse = "")
k235Has <- function(text, pattern) grepl(pattern, k235Text(text), fixed = TRUE)

#' The assertions every case shares (P02): the extractor alone returns a non-empty list, so a stop seen in Main is the
#' template's and not the harness's.
k235ExpectExtractor <- function(got, info) expect_gt(length(got$extracted), 0L, label = paste("extractor length,", info))

#' Main's report reached the Methods: no handler message, the Methods heading printed
k235ExpectMethods <- function(got, info) {
    expect_false(k235Has(got$content$sample_details, "Error running Main"), info = info)
    expect_true(k235Has(got$content$method_details, "METHOD DETAILS"), info = info)
}

# --- T1 to T5 ----------------------------------------------------------------------------------------------------------

test_that("KIT-235 T1: a writer file (M metab, L2 lipid) reaches the Methods through Main", {
    for (case in c("M", "L2")) {
        got <- k235Run(k235Cases[[case]]$omic, k235WriterLines(case))
        k235ExpectExtractor(got, case)
        k235ExpectMethods(got, case)
    }
    got <- k235Run("metab", k235WriterLines("M"))
    expect_true(k235Has(got$content$method_details, "The dataset contains 2 assays: assay_pos, assay_neg."))
})

test_that("KIT-235 T2: the per-assay RUV lines are read as the assays they name (M)", {
    got <- k235Run("metab", k235WriterLines("M"))
    k235ExpectExtractor(got, "M")
    expect_identical(names(got$extracted$ruv_assay_list), c("assay_pos", "assay_neg"))
    expect_true(k235Has(got$content$method_details,
        "For assay_pos, 2 unwanted components were removed using 1 empirical negative control metabolites."))
    expect_true(k235Has(got$content$method_details,
        "For assay_neg, 3 unwanted components were removed using 2 empirical negative control metabolites."))
})

test_that("KIT-235 T3: the writer's fields are read (M, extractor, no record)", {
    got <- k235Run("metab", k235WriterLines("M"))
    k235ExpectExtractor(got, "M")
    x <- got$extracted
    expect_identical(x$num_assays, "2")
    expect_identical(x$assay_names, "assay_pos, assay_neg")
    expect_identical(x$itsd_applied, TRUE)
    expect_identical(x$itsd_method_type, "median")
    expect_identical(x$itsd_aggregation, "sum")
    expect_identical(x$itsd_removed, TRUE)
    expect_identical(x$log_transform_offset, "1")
    expect_identical(x$normalisation_method, "PQN")
    expect_identical(x$ruv_grouping_variable, "group")
})

test_that("KIT-235 T3b: the record's ITSD aggregation is kept through the overlay (M with a record)", {
    record <- list(k235Step("metab_itsd_normalisation", list(itsd_aggregation = "sum")))
    got <- k235Run("metab", k235WriterLines("M"), record = record)
    k235ExpectExtractor(got, "M with record")
    expect_identical(got$numbers$itsd_aggregation, "sum")
})

test_that("KIT-235 T4: the single RUV form (L1) reaches the Methods without a false k", {
    got <- k235Run("lipid", k235WriterLines("L1"))
    k235ExpectExtractor(got, "L1")
    k235ExpectMethods(got, "L1")
    method <- got$content$method_details
    expect_true(k235Has(method, "The dataset contains a single assay: assay_main."))
    expect_true(k235Has(method, "Data were analyzed on the original intensity scale without log transformation."))
    expect_true(k235Has(method,
        "RUV correction was applied independently to each assay using the 'group' variable to preserve biological signal."))
    expect_false(k235Has(method, "NA unwanted"))
})

test_that("KIT-235 T5: a writer file with RUV not applied says so (S lipid, S metab)", {
    for (case in c("S_lipid", "S_metab")) {
        got <- k235Run(k235Cases[[case]]$omic, k235WriterLines(case))
        k235ExpectExtractor(got, case)
        expect_true(k235Has(got$content$method_details, "RUV-III batch correction was not applied."), info = case)
    }
})

# --- T6: the older shapes print what they printed (a guard; passes on base by construction) -----------------------------

test_that("KIT-235 T6: '\u2022' and '-' files print the base text, and the same text as the '*' file (M, L1)", {
    for (case in c("M", "L1")) {
        omic <- k235Cases[[case]]$omic
        star <- k235WriterLines(case)
        head_star <- k235Text(k235Run(omic, star)$content$method_details)
        for (shape in c("bullet", "dash")) {
            lines <- k235Shape(star, shape)
            head <- k235Run(omic, lines)
            base <- k235Run(omic, lines, side = "base")
            info <- paste(case, shape)
            k235ExpectExtractor(head, info)
            expect_true(k235Has(head$content$method_details, "METHOD DETAILS"), info = info)
            expect_identical(k235Text(head$content$method_details), k235Text(base$content$method_details), info = info)
            expect_identical(k235Text(head$content$method_details), head_star, info = info)
        }
    }
})
