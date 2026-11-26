# KIT-272: the proteomics reports describe panels d to f of Figure 1 (c to d with two stages) as what they are: notched
# box plots of the samples' scores on principal components 1 and 2, by group (the owner's decisions VAL-316 (7) and
# VAL-332). Lipidomics and metabolomics get their own sentence (their box plots come from a separate PCA).
# The sentences are DRAFTS for the owner's wording pass: they are held in the constants below (one place), and a later
# wording pass changes the constants and the templates together. FormatResults and the Methods functions are called
# directly, from the template's own chunks (copies of the helpers of test-kit236-report-sentences.R, never sourced);
# the frozen base copies are tests/testdata/kit227/base.

k272Spec <- list(
    LFQ = list(file = "LFQ_report.rmd", method = "FormatMethodDetailsLFQ"),
    TMT = list(file = "TMT_report.rmd", method = "FormatMethodDetailsTMT"),
    DIANN = list(file = "DIANN_report.rmd", method = "FormatMethodDetails"),
    DIANN_limpa = list(file = "DIANN_limpa_report.rmd", method = "FormatMethodDetails"))
k272Names <- names(k272Spec)

# --- The draft sentences ---------------------------------------------------------------------------------------------

k272Variants <- list(
    three = list(box = " (panels d-f)", pca = " (panels a-c)", ref = "panels a-c"),
    two = list(box = " (panels c-d)", pca = " (panels a-b)", ref = "panels a-b"),
    unrecorded = list(box = "", pca = "", ref = "the principal component analysis panels"))

k272OldBullet <- list(
    three = "• Density plots (panels d-f): the distribution of protein intensities per sample.",
    two = "• Density plots (panels c-d): the distribution of protein intensities per sample.",
    unrecorded = "• Density plots: the distribution of protein intensities per sample.")

k272NewBullet <- lapply(k272Variants, function(v) paste0(
    "• Principal component scores by group", v$box, ": for each experimental group, a box plot of the samples' ",
    "scores on principal component 1 (upper plot) and principal component 2 (lower plot), the same scores that ", v$ref,
    " plot as points. The line in each box is the median, the box spans the first to the third quartile, the whiskers ",
    "reach the most extreme scores within 1.5 times the interquartile range of the box, and scores beyond them are drawn ",
    "as points. The notch spans the median plus and minus 1.58 times the interquartile range divided by the square root ",
    "of the number of samples in the group: as a rough guide, two groups whose notches do not overlap have medians that ",
    "differ at about the 5% level (McGill, Tukey and Larsen, 1978). With fewer than about ten samples in a group the ",
    "notch extends beyond the box and should not be read. Each stage is a separate analysis, so the scores of different ",
    "panels are not on a common scale and their signs are arbitrary."))

k272OldPca <- lapply(k272Variants, function(v) paste0("• Principal component analysis", v$pca,
    ": samples plotted on the first two principal components,\n  coloured by experimental group."))
k272NewPca <- lapply(k272Variants, function(v) paste0("• Principal component analysis", v$pca,
    ": samples plotted on the first two principal components of the log2 protein values (centred, not scaled), ",
    "computed at each stage from the tenth of the proteins with the largest ratio of standard deviation to mean of ",
    "their log2 values; coloured by experimental group."))

k272OldMethod <- "density boxplots"
k272NewMethod <- "box plots of the principal component scores"

k272OldLipidMetab <- "• Density boxplots to assess distribution consistency"
k272NewLipidMetab <- c(
    lipid = paste0("• Box plots of the samples' scores on the first two principal components of every lipid with ",
        "at least two values, a separate analysis from the PCA panels, which use the tenth of the lipids with the ",
        "largest ratio of standard deviation to mean"),
    metab = paste0("• Box plots of the samples' scores on the first two principal components of every metabolite ",
        "with at least two values, a separate analysis from the PCA panels, which use the tenth of the metabolites with ",
        "the largest ratio of standard deviation to mean"))

k272PlaceholderSentence <- paste0("In Figure 1, panel d (principal component scores by group after log2 ",
    "transformation) could not be drawn by the application; it shows a note giving the reason instead of the plot.")

# --- Helpers (copies of the KIT-236 test's) ---------------------------------------------------------------------------

#' A template file: the package's (head) or the frozen base copy (KIT-227, verbatim at cf8af245)
k272Path <- function(file, side = "head", dir = "proteomics") {
    if (identical(side, "head")) system.file("reports", dir, file, package = "MultiScholaR") else
        testthat::test_path("..", "testdata", "kit227", "base", file)
}

k272Chunks <- function(path) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    starts <- grep("^```\\{r", lines)
    lapply(starts, function(s) {
        e <- s + which(grepl("^```[[:space:]]*$", lines[-seq_len(s)]))[[1L]]
        list(label = trimws(sub("^```\\{r[[:space:],]*([^,}]*).*$", "\\1", lines[[s]])),
            code = if (e > s + 1L) lines[(s + 1L):(e - 1L)] else character())
    })
}

#' The template's functions (every `name <- function` of every chunk) in one environment, with the report's paths
k272Env <- function(name, root, side = "head") {
    env <- new.env(parent = asNamespace("MultiScholaR"))
    for (chunk in k272Chunks(k272Path(k272Spec[[name]]$file, side))) {
        for (e in parse(text = chunk$code, keep.source = FALSE)) {
            if (is.call(e) && identical(e[[1L]], as.name("<-")) && is.call(e[[3L]]) &&
                identical(e[[3L]][[1L]], as.name("function"))) eval(e, env)
        }
    }
    env$`%||%` <- function(a, b) if (!is.null(a)) a else b
    env$d66_log <- function(...) invisible(NULL)
    env$tibble <- tibble::tibble
    env$params <- list(omic_type = "proteomics", experiment_label = "kit272", workflow_name = "kit272")
    env$report_source_dir <- file.path(root, "scripts")
    env$report_results_summary_dir <- file.path(root, "summary")
    env$report_qc_figures_dir <- file.path(root, "summary", "QC_figures")
    env$report_publication_figures_dir <- file.path(root, "summary", "Publication_figures")
    env$report_publication_tables_dir <- file.path(root, "summary", "Publication_tables")
    for (d in c(env$report_source_dir, env$report_qc_figures_dir, env$report_publication_figures_dir,
        env$report_publication_tables_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
    env
}

k272Norm <- function(x) trimws(gsub("[[:space:]]+", " ", paste(as.character(x), collapse = " ")))
k272Has <- function(text, phrase) grepl(k272Norm(phrase), k272Norm(text), fixed = TRUE)
k272Count <- function(text, pattern) lengths(regmatches(text, gregexpr(pattern, text, fixed = TRUE)))

#' FormatResults as printed, for one RUV state (TRUE, FALSE, NA); `stages_known = FALSE` is the limpa route without RUV
k272Results <- function(name, ruv, side = "head", stages_known = TRUE) {
    env <- k272Env(name, withr::local_tempdir(pattern = "kit272-"), side)
    args <- list(correlation_data = list(), figures_dir = tempdir(), ruv_applied = ruv)
    if (!stages_known) args$stages_known <- FALSE
    paste(as.character(do.call(env$FormatResults, args)), collapse = "")
}

k272RuvOf <- c(three = TRUE, two = FALSE, unrecorded = NA)

#' The Methods text of the template's own function from the numbers AnalyzeSamples gives an empty study; `limpa_route`
#' sets the two numbers that select the limpa route's own paragraph
k272Method <- function(name, limpa_route = FALSE, side = "head") {
    env <- k272Env(name, withr::local_tempdir(pattern = "kit272-m-"), side)
    runs <- c("KO_1", "KO_2", "KO_3", "WT_1", "WT_2", "WT_3")
    design <- data.frame(Run = runs, group = sub("_[0-9]+$", "", runs), replicates = as.integer(sub("^.*_", "", runs)),
        tech_rep_group = runs, stringsAsFactors = FALSE)
    numbers <- suppressWarnings(env$AnalyzeSamples(design, list()))$numbers
    stopifnot(length(numbers) > 0L)
    if (limpa_route) {
        numbers$record_present <- TRUE
        numbers$dpc_engine <- "fastlimpadpc"
    }
    paste(as.character(env[[k272Spec[[name]]$method]](numbers)), collapse = "")
}

#' The lines of one chunk of a template, or NULL
k272Chunk <- function(file, label) {
    lines <- readLines(k272Path(file), warn = FALSE, encoding = "UTF-8")
    start <- grep(paste0("^```\\{r ", label, "[,}]"), lines)
    if (length(start) != 1L) return(NULL)
    lines[start:(start + which(lines[-seq_len(start)] == "```")[[1L]])]
}

#' What the placeholder chunk prints for a list written beside a combined figure (copy of KIT-259's `k259bSentence`)
k272Placeholder <- function(file, list) {
    chunk <- k272Chunk(file, "QC_metrics_placeholders")
    figures <- withr::local_tempdir()
    feature <- withr::local_tempdir()
    file.create(file.path(figures, "proteomics_composite_QC_figure.png"))
    writeLines(list, file.path(figures, "proteomics_composite_QC_figure.placeholders"))
    env <- new.env(parent = globalenv())
    env$params <- list(omic_type = "proteomics")
    env$report_qc_figures_dir <- figures
    env$current_project_paths <- list(feature_qc_dir = feature)
    env$results_data <- list(numbers = list(ruv_applied = FALSE))
    rmd <- tempfile(fileext = ".Rmd")
    writeLines(chunk, rmd, useBytes = TRUE)
    code <- tempfile(fileext = ".R")
    suppressMessages(knitr::purl(rmd, output = code, quiet = TRUE, documentation = 0))
    paste(utils::capture.output(for (e in parse(code, keep.source = FALSE, encoding = "UTF-8")) eval(e, env)),
        collapse = "\n")
}

# --- T1 to T3, T10: the Results text -----------------------------------------------------------------------------------

for (name in k272Names) {
    test_that(sprintf("KIT-272 T1 to T3, T10 (%s): the Results bullets say what panels d to f are, in every RUV state", name), {
        for (variant in names(k272Variants)) {
            text <- k272Results(name, k272RuvOf[[variant]])
            info <- paste(name, variant)
            expect_true(k272Has(text, k272NewBullet[[variant]]), info = info)
            expect_true(k272Has(text, k272NewPca[[variant]]), info = info)
            expect_false(k272Has(text, "protein intensities per sample"), info = info)
            expect_false(k272Has(text, "Density plots"), info = info)
        }
        # The limpa route without RUV (the layout is not known): the same text as RUV not recorded
        if (identical(name, "DIANN_limpa")) {
            text <- k272Results(name, FALSE, stages_known = FALSE)
            expect_true(k272Has(text, k272NewBullet$unrecorded), info = "limpa, stages not known")
            expect_true(k272Has(text, k272NewPca$unrecorded), info = "limpa, stages not known")
            expect_false(k272Has(text, "Density plots"), info = "limpa, stages not known")
        }
    })

    # --- T5: the figure's facts, in the words of the draft ---------------------------------------------------------
    test_that(sprintf("KIT-272 T5 (%s): the bullet names the scores, the group, the box, the whiskers, the notch and its limit", name), {
        for (variant in names(k272Variants)) {
            text <- k272Results(name, k272RuvOf[[variant]])
            info <- paste(name, variant)
            for (fact in c("Principal component scores by group", "principal component 1 (upper plot)",
                "principal component 2 (lower plot)", "for each experimental group",
                "the same scores that", "plot as points", "The line in each box is the median",
                "the box spans the first to the third quartile",
                "the whiskers reach the most extreme scores within 1.5 times the interquartile range of the box",
                "1.58 times the interquartile range divided by the square root of the number of samples in the group",
                "two groups whose notches do not overlap have medians that differ at about the 5% level",
                "(McGill, Tukey and Larsen, 1978)",
                "With fewer than about ten samples in a group the notch extends beyond the box and should not be read",
                "Each stage is a separate analysis", "not on a common scale", "their signs are arbitrary")) {
                expect_true(k272Has(text, fact), info = paste(info, fact))
            }
            expect_false(grepl("95%", text, fixed = TRUE), info = info)
            expect_false(k272Has(text, "confidence interval"), info = info)
        }
    })

    # --- T6: the Methods phrase ------------------------------------------------------------------------------------
    test_that(sprintf("KIT-272 T6 (%s): the Methods quality-control sentence names the box plots of the scores", name), {
        routes <- if (identical(name, "DIANN_limpa")) c(FALSE, TRUE) else FALSE
        for (route in routes) {
            text <- k272Method(name, limpa_route = route)
            info <- paste(name, "limpa route:", route)
            expect_identical(k272Count(k272Norm(text), k272NewMethod), 1L, info = info)
            expect_false(k272Has(text, "density boxplots"), info = info)
        }
    })

    # --- T7: KIT-259's placeholder label ---------------------------------------------------------------------------
    test_that(sprintf("KIT-272 T7 (%s): a replaced panel d is named by what it plots", name), {
        file <- k272Spec[[name]]$file
        expect_false(is.null(k272Chunk(file, "QC_metrics_placeholders")))
        skip_if(is.null(k272Chunk(file, "QC_metrics_placeholders")))
        expect_identical(k272Placeholder(file, "d\tpre_norm_density"), k272PlaceholderSentence)
    })

    # --- T4: no template carries the old words anywhere ------------------------------------------------------------
    test_that(sprintf("KIT-272 T4 (%s): the old description appears nowhere in the template", name), {
        lines <- readLines(k272Path(k272Spec[[name]]$file), warn = FALSE, encoding = "UTF-8")
        expect_identical(sum(grepl("protein intensities per sample", lines, fixed = TRUE)), 0L)
        expect_identical(sum(grepl("Density plots", lines, fixed = TRUE)), 0L)
    })

    # --- T8, G1': guards ------------------------------------------------------------------------------------------
    test_that(sprintf("KIT-272 T8 (%s): nothing else of the Results text changed (head, the new bullets put back, equals the base)", name), {
        for (variant in names(k272Variants)) {
            head <- k272Results(name, k272RuvOf[[variant]])
            base <- k272Results(name, k272RuvOf[[variant]], side = "base")
            back <- sub(k272NewBullet[[variant]], k272OldBullet[[variant]], head, fixed = TRUE)
            back <- sub(k272NewPca[[variant]], k272OldPca[[variant]], back, fixed = TRUE)
            expect_identical(back, base, info = paste(name, variant))
            for (kept in c("Relative log expression", "Pearson correlation",
                "Whether technical variation was reduced and biological signal retained")) {
                expect_true(k272Has(head, kept), info = paste(name, variant, kept))
            }
        }
    })

    test_that(sprintf("KIT-272 G1' (%s): the frozen base copy holds each old element exactly once (the declared edits' match assertions)", name), {
        for (variant in names(k272Variants)) {
            base <- k272Results(name, k272RuvOf[[variant]], side = "base")
            expect_identical(k272Count(base, k272OldBullet[[variant]]), 1L, info = paste(name, variant, "bullet"))
            expect_identical(k272Count(base, k272OldPca[[variant]]), 1L, info = paste(name, variant, "PCA"))
        }
        base_method <- k272Method(name, limpa_route = FALSE, side = "base")
        expect_identical(k272Count(k272Norm(base_method), k272OldMethod), 1L, info = paste(name, "Methods"))
    })

    test_that(sprintf("KIT-272 T9 (%s): the variants render as one string", name), {
        for (variant in names(k272Variants)) {
            text <- expect_no_error(k272Results(name, k272RuvOf[[variant]]))
            expect_true(is.character(text) && length(text) == 1L && nzchar(text), info = paste(name, variant))
        }
    })
}

# --- T11: lipidomics and metabolomics -----------------------------------------------------------------------------------

test_that("KIT-272 T11: the lipidomics and metabolomics bullets describe their own box plots, in the Methods and when RUV is not recorded", {
    for (omic in c("lipid", "metab")) {
        file <- if (identical(omic, "lipid")) "lipidomics_report.rmd" else "metabolomics_report.rmd"
        dir <- if (identical(omic, "lipid")) "lipidomics" else "metabolomics"
        text <- k272Norm(readLines(k272Path(file, dir = dir), warn = FALSE, encoding = "UTF-8"))
        expect_identical(k272Count(text, k272Norm(k272NewLipidMetab[[omic]])), 2L, info = omic)
        expect_identical(k272Count(text, "Density boxplots to assess distribution consistency"), 0L, info = omic)
    }
})
