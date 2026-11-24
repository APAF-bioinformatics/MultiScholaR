# KIT-231: the shipped workbooks read folders that setupDirectories() returns, and pass only arguments the package
# functions take. The DIA workbooks asked for `de_output_dir` (the set-up has returned `da_output_dir` since 79d80356),
# passed `de_output_dir =` to createDAResultsForEnrichment() (formal `da_output_dir` since a6e73542), the experienced
# workbook read never-assigned bare folder names, limpa_starter read the integration folder of a set-up without
# integration, limpa_experienced read `.x$da_proteins_long` (the DA call returns `de_proteins_long`), and RunApplet()
# saved its result to `output_dir`, which the set-up never returns. N1 to N3 read the installed workbooks and report
# templates as text; N4 evaluates the workbooks' own DA output and enrichment chunks on a real ProteinQuantitativeData
# (KIT-200 section 4, builders copied from tests/testthat/test-kit200-prot-da-no-value-proteins.R on branch
# kit-200-prot-da-proteins-without-values), with only writeInteractiveVolcanoPlotProteomics mocked; N5 runs RunApplet()
# with the testServer mocks of test-applets-direct-shared.R:172-209 (cf8af245), copied. Each test first asserts that its
# call completed, so the tests-only commit fails by assertion at the stated message.

kit231Label <- "kit231"

kit231Quiet <- function(expr) {
    value <- NULL
    utils::capture.output(utils::capture.output(value <- suppressWarnings(suppressMessages(expr)), type = "message"))
    value
}

kit231SetupDirs <- function(omic_types, label, env = parent.frame()) {
    base <- withr::local_tempdir("kit231-project-", .local_envir = env)
    kit231Quiet(setupDirectories(base_dir = base, omic_types = omic_types, label = label, force = TRUE))
}

# The workbook's R chunks: label, first code line (file line numbers), code.
kit231Chunks <- function(path) {
    lines <- readLines(path, warn = FALSE)
    is_header <- grepl("^\\s*```+\\s*\\{[rR]([ ,].*)?\\}\\s*$", lines)
    is_fence <- grepl("^\\s*```+\\s*$", lines)
    chunks <- list()
    i <- 1L
    while (i <= length(lines)) {
        if (is_header[[i]]) {
            ends <- which(is_fence & seq_along(lines) > i)
            j <- if (length(ends)) ends[[1L]] else length(lines) + 1L
            inner <- sub("\\}\\s*$", "", sub("^\\s*```+\\s*\\{[rR]\\s*,?\\s*", "", lines[[i]]))
            label <- trimws(sub(",.*$", "", inner))
            if (grepl("=", label)) label <- ""
            chunks[[length(chunks) + 1L]] <- list(label = label, code_start = i + 1L,
                code = if (j > i + 1L) lines[(i + 1L):(j - 1L)] else character())
            i <- j + 1L
        } else {
            i <- i + 1L
        }
    }
    chunks
}

kit231Chunk <- function(chunks, label) {
    hit <- Filter(function(chunk) identical(chunk$label, label), chunks)
    if (!length(hit)) stop("chunk not found: ", label)
    hit[[1L]]
}

kit231Parse <- function(chunk) {
    exprs <- parse(text = chunk$code, keep.source = TRUE)
    list(exprs = exprs, data = if (length(exprs)) utils::getParseData(exprs) else NULL)
}

kit231Literal <- function(chunks, name) {
    code <- unlist(lapply(chunks, `[[`, "code"))
    hit <- Filter(length, regmatches(code, regexec(paste0("^\\s*", name, "\\s*<-\\s*\"([^\"]*)\""), code)))
    if (length(hit)) hit[[1L]][[2L]] else NA_character_
}

# The workbook's own set-up: the uncommented arguments of its setupDirectories() call, its experiment_label and
# omic_type assignments (each DIA workbook has a commented `#omic_types = "metabolomics"` line before the real one).
kit231WorkbookSetup <- function(chunks) {
    label <- kit231Literal(chunks, "experiment_label")
    env <- new.env(parent = globalenv())
    if (!is.na(label)) env$experiment_label <- label
    omic <- kit231Literal(chunks, "omic_type")
    if (!is.na(omic)) env$omic_type <- omic
    for (chunk in chunks) {
        parsed <- kit231Parse(chunk)
        pd <- parsed$data
        if (is.null(pd)) next
        hit <- pd$id[pd$token == "SYMBOL_FUNCTION_CALL" & pd$text == "setupDirectories"]
        if (!length(hit)) next
        call_id <- pd$parent[pd$id == pd$parent[pd$id == hit[[1L]]]]
        args <- as.list(str2lang(utils::getParseText(pd, call_id)))[-1L]
        env$kit231_omic_types <- eval(args[["omic_types"]], env)
        env$kit231_label <- eval(args[["label"]], env)
        break
    }
    env
}

kit231Kids <- function(pd, id) {
    kids <- pd[pd$parent == id, , drop = FALSE]
    kids[order(kids$line1, kids$col1), , drop = FALSE]
}

# Every `<object>$<name>` read in a chunk whose object is `holder` itself (`holder$x$name`, `holder[[...]]$name` when
# `depth` is 2; `holder$name` when `depth` is 1): the read's line, its text and the object expression's text.
kit231Reads <- function(chunk, holder, depth) {
    pd <- kit231Parse(chunk)$data
    if (is.null(pd)) return(list())
    out <- list()
    for (dollar in pd$parent[pd$token == "'$'"]) {
        kids <- kit231Kids(pd, dollar)
        lhs <- kids$id[kids$token == "expr"][[1L]]
        inner <- if (depth == 2L) {
            lhs_kids <- kit231Kids(pd, lhs)
            lhs_kids$id[lhs_kids$token == "expr"][1L]
        } else {
            lhs
        }
        if (is.na(inner) || !identical(utils::getParseText(pd, inner), holder)) next
        out[[length(out) + 1L]] <- list(line = chunk$code_start - 1L + pd$line1[pd$id == dollar],
            text = utils::getParseText(pd, dollar), lhs = utils::getParseText(pd, lhs),
            name = sub("^.*\\$\\s*", "", utils::getParseText(pd, dollar)))
    }
    out
}

kit231WorkbookFiles <- function() {
    sort(list.files(system.file("workbooks", package = "MultiScholaR"), pattern = "\\.[Rr]md$", recursive = TRUE,
        full.names = TRUE))
}

kit231DiaWorkbooks <- c("DIA_workflow_starter.rmd", "DIA_workflow_limpa.rmd", "DIA_workflow_limpa_starter.rmd",
    "DIA_workflow_limpa_experienced.rmd", "DIA_workflow_experienced.rmd")

kit231DiaPath <- function(name) system.file("workbooks", "proteomics", name, package = "MultiScholaR")

test_that("N1: every directories element a workbook or report template reads exists in what setupDirectories returns", {
    misses <- character()
    seen <- 0L
    # Workbooks: each read evaluated against the workbook's own set-up (`$` partial-matches keys, as R does).
    for (path in kit231WorkbookFiles()) {
        chunks <- kit231Chunks(path)
        setup <- kit231WorkbookSetup(chunks)
        if (!is.null(setup$kit231_omic_types)) {
            setup$project_dirs <- kit231SetupDirs(setup$kit231_omic_types, setup$kit231_label)
        }
        for (chunk in chunks) {
            for (read in kit231Reads(chunk, "project_dirs", 2L)) {
                seen <- seen + 1L
                key <- tryCatch(eval(str2lang(read$lhs), setup), error = function(e) NULL)
                value <- tryCatch(eval(str2lang(read$text), setup), error = function(e) NULL)
                if (is.null(value)) {
                    misses <- c(misses, sprintf("%s:%d key %s %s", basename(path), read$line,
                        if (is.null(key)) paste(read$lhs, "(absent)") else read$lhs, read$name))
                }
            }
        }
    }
    # Report templates: `current_project_paths$<name>` against the names the set-up returns for the template's omic
    # (folder name), or for any omic (shared templates).
    omics <- c("proteomics", "metabolomics", "transcriptomics", "lipidomics", "integration")
    all_dirs <- kit231SetupDirs(omics, kit231Label)
    names_for <- lapply(stats::setNames(omics, omics), function(omic) names(all_dirs[[paste0(omic, "_", kit231Label)]]))
    reports <- sort(list.files(system.file("reports", package = "MultiScholaR"), pattern = "\\.[Rr]md$",
        recursive = TRUE, full.names = TRUE))
    for (path in reports) {
        omic <- basename(dirname(path))
        allowed <- if (omic %in% omics) names_for[[omic]] else unique(unlist(names_for))
        for (chunk in kit231Chunks(path)) {
            for (read in kit231Reads(chunk, "current_project_paths", 1L)) {
                seen <- seen + 1L
                if (!read$name %in% allowed) {
                    misses <- c(misses, sprintf("%s:%d %s", basename(path), read$line, read$name))
                }
            }
        }
    }
    expect_gt(seen, 150L)
    expect_identical(misses, character())
})

kit231DaCalls <- c("createDAResultsForEnrichment", "outputDeAnalysisResults")

test_that("N2: the DIA workbooks' createDAResultsForEnrichment and outputDeAnalysisResults calls use only their formals", {
    misses <- character()
    calls <- 0L
    for (name in kit231DiaWorkbooks) {
        for (chunk in kit231Chunks(kit231DiaPath(name))) {
            pd <- kit231Parse(chunk)$data
            if (is.null(pd)) next
            for (token in pd$id[pd$token == "SYMBOL_FUNCTION_CALL" & pd$text %in% kit231DaCalls]) {
                fn_name <- pd$text[pd$id == token]
                call_id <- pd$parent[pd$id == pd$parent[pd$id == token]]
                dummy <- function() NULL
                formals(dummy) <- formals(get(fn_name, envir = asNamespace("MultiScholaR")))
                message <- tryCatch({
                    match.call(dummy, str2lang(utils::getParseText(pd, call_id)))
                    ""
                }, error = function(e) conditionMessage(e))
                calls <- calls + 1L
                if (nzchar(message)) {
                    misses <- c(misses, sprintf("%s:%d %s: %s", name, chunk$code_start - 1L + pd$line1[pd$id == call_id],
                        fn_name, message))
                }
            }
        }
    }
    expect_identical(calls, 10L)
    expect_identical(misses, character())
})

kit231Assigned <- function(exprs) {
    out <- character()
    walk <- function(x) {
        if (!is.call(x)) return(invisible(NULL))
        if (is.name(x[[1L]]) && as.character(x[[1L]]) %in% c("<-", "=", "<<-") && length(x) >= 3L) {
            target <- x[[2L]]
            while (is.call(target) && length(target) >= 2L) target <- target[[2L]]
            if (is.name(target)) out <<- c(out, as.character(target))
        }
        parts <- as.list(x)
        for (k in seq_along(parts)) {
            if (identical(parts[[k]], quote(expr = ))) next
            walk(parts[[k]])
        }
        invisible(NULL)
    }
    for (k in seq_along(exprs)) walk(exprs[[k]])
    unique(out)
}

test_that("N3: the experienced workbook's DA output and enrichment chunks read no folder name the workbook never assigns", {
    chunks <- kit231Chunks(kit231DiaPath("DIA_workflow_experienced.rmd"))
    assigned <- unique(unlist(lapply(chunks, function(chunk) kit231Assigned(kit231Parse(chunk)$exprs))))
    free <- character()
    for (label in c("Output DE Analysis Results", "Enrichment Analysis")) {
        fn <- function() NULL
        body(fn) <- as.call(c(as.name("{"), as.list(kit231Parse(kit231Chunk(chunks, label))$exprs)))
        symbols <- codetools::findGlobals(fn, merge = FALSE)$variables
        free <- c(free, setdiff(grep("_dir$", symbols, value = TRUE), assigned))
    }
    expect_identical(sort(unique(free)), character())
})

# ---- N4: the workbooks' own chunks on a real object (KIT-200 section 4) -------------------------------------------

kit231WithSeed <- function(seed, code) withr::with_seed(seed, code, .rng_kind = "Mersenne-Twister",
    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

# Copied from kit200Matrix, kit200Design, kit200Object (test-kit200-prot-da-no-value-proteins.R).
kit231Matrix <- function(sizes = c(A = 4L, B = 4L, C = 4L), seed = 200001L, proteins = 60L) {
    groups <- rep(names(sizes), sizes)
    x <- kit231WithSeed(seed, {
        level <- stats::rnorm(proteins, 22, 2)
        effect <- as.numeric(seq_len(proteins) %in% 11:16)
        level + outer(effect, as.numeric(groups == "B")) +
            matrix(stats::rnorm(proteins * length(groups), sd = 0.3), proteins)
    })
    dimnames(x) <- list(paste0("P", seq_len(proteins)), paste0("S", seq_along(groups)))
    x
}

kit231Design <- function(sizes = c(A = 4L, B = 4L, C = 4L)) {
    runs <- paste0("S", seq_len(sum(sizes)))
    data.frame(Run = runs, group = rep(names(sizes), sizes), replicates = runs, row.names = runs,
        stringsAsFactors = FALSE)
}

kit231Object <- function(x, design) {
    rownames(design) <- NULL
    ProteinQuantitativeData(protein_quant_table = data.frame(Protein.Ids = rownames(x), x, row.names = NULL,
        check.names = FALSE, stringsAsFactors = FALSE), design_matrix = design, sample_id = "Run", group_id = "group",
        technical_replicate_id = "replicates", protein_id_column = "Protein.Ids",
        protein_id_table = data.frame(Protein.Ids = rownames(x), gene_names = paste0("GENE_", rownames(x)),
            stringsAsFactors = FALSE), args = list())
}

# Evaluate a chunk's expressions in order, as knitr does; stop at the first error. Messages are collected.
kit231Evaluate <- function(code, env) {
    messages <- character()
    error <- NULL
    utils::capture.output(error <- tryCatch(withCallingHandlers({
        for (expr in parse(text = code)) eval(expr, env)
        NULL
    }, message = function(m) {
        messages <<- c(messages, conditionMessage(m))
        invokeRestart("muffleMessage")
    }, warning = function(w) invokeRestart("muffleWarning")), error = function(e) conditionMessage(e)))
    list(error = error, messages = messages)
}

kit231DaChunks <- list(
    DIA_workflow_starter.rmd = c(out = "34 Output DE Analysis Results", enr = "35 Enrichment Analysis"),
    DIA_workflow_limpa.rmd = c(out = "34 Output DE Analysis Results", enr = "35 Enrichment Analysis"),
    DIA_workflow_limpa_starter.rmd = c(out = "20 Output DE Analysis Results", enr = "21 Enrichment Analysis"),
    DIA_workflow_limpa_experienced.rmd = c(out = "20 Output DE Analysis Results", enr = "21 Enrichment Analysis"),
    DIA_workflow_experienced.rmd = c(out = "Output DE Analysis Results", enr = "Enrichment Analysis")
)

test_that("N4: each DIA workbook's DA output chunk writes its tables and its enrichment call reads them", {
    object <- kit231Object(kit231Matrix(), kit231Design())
    contrasts_tbl <- tibble::tibble(contrasts = "B_vs_A=groupB-groupA")
    da_list <- list(B_vs_A = kit231Quiet(protein_deAnalysisWrapperFunction(object, dplyr::slice(contrasts_tbl, 1L),
        formula_string = "~ 0 + group", de_q_val_thresh = 0.05, treat_lfc_cutoff = 0, eBayes_trend = TRUE,
        eBayes_robust = TRUE, args_group_pattern = "(\\d+)", args_row_id = "Protein.Ids")))
    local_mocked_bindings(writeInteractiveVolcanoPlotProteomics = function(...) invisible(NULL))
    # The workbooks attach these with loadDependencies(); limpa_experienced's chunk calls imap_dfr() and filter() bare.
    withr::local_package("purrr")
    withr::local_package("dplyr")
    withr::local_dir(withr::local_tempdir("kit231-wd-"))

    for (name in names(kit231DaChunks)) {
        chunks <- kit231Chunks(kit231DiaPath(name))
        env <- kit231WorkbookSetup(chunks)
        env$project_dirs <- kit231SetupDirs(env$kit231_omic_types, env$kit231_label)
        env$config_list <- list(deAnalysisParameters = list(plots_format = c("pdf", "png")))
        env$ruv_normalised_for_da_analysis_obj <- object
        env$uniprot_dat_cln <- data.frame(Entry = paste0("P", 1:60), gene_names = paste0("GENE_P", 1:60))
        env$da_analysis_results_list <- da_list
        env$contrasts_tbl <- contrasts_tbl
        env$design_matrix <- object@design_matrix
        da_dir <- env$project_dirs[[paste0("proteomics_", env$kit231_label)]]$da_output_dir

        output <- kit231Evaluate(kit231Chunk(chunks, kit231DaChunks[[name]][["out"]])$code, env)
        expect_null(output$error, label = paste(name, "DA output chunk error"))
        expect_identical(grep("Error processing contrast", output$messages, value = TRUE), character(),
            label = paste(name, "messages"))
        expect_true(file.exists(file.path(da_dir, "da_proteins_B_vs_A_long_annot.tsv")), label = paste(name, "long_annot"))

        enrichment <- kit231Chunk(chunks, kit231DaChunks[[name]][["enr"]])
        first <- as.character(attr(parse(text = enrichment$code, keep.source = TRUE), "srcref")[[1L]])
        call <- kit231Evaluate(first, env)
        expect_null(call$error, label = paste(name, "createDAResultsForEnrichment error"))
        if (is.null(call$error)) {
            da_data <- env$de_results_for_enrichment@da_data
            expect_length(da_data, 1L)
            expect_s3_class(da_data[[1L]], "data.frame")
        }

        if (identical(name, "DIA_workflow_limpa_experienced.rmd")) {
            top <- kit231Evaluate(kit231Chunk(chunks, "extract-top-de-genes")$code, env)
            expect_null(top$error, label = paste(name, "extract-top-de-genes error"))
            expect_true(file.exists(file.path(da_dir, "top_de_genes.tsv")), label = paste(name, "top_de_genes.tsv"))
        }
    }
})

# ---- N5: RunApplet on a setupDirectories() project ---------------------------------------------------------------

test_that("N5: RunApplet runs the design applet on a setupDirectories() project and saves its result in source_dir", {
    skip_if_not(workflowOwnerCapable())
    report_path <- file.path(withr::local_tempdir("kit231-applet-data-"), "report.tsv")
    utils::write.table(data.frame(
        Run = rep(c("S10", "S2", "S1"), each = 2L), Protein.Group = rep(c("P1", "P2"), 3L),
        Stripped.Sequence = rep(c("AAK", "CCK"), 3L), Precursor.Quantity = c(10, 20, 30, 40, 50, 60),
        Q.Value = 0.001), report_path, sep = "\t", row.names = FALSE, quote = FALSE)

    # Copied from test-applets-direct-shared.R:172-209 (cf8af245): runApp drives the applet's own app object through
    # testServer, builds a two-sample design and clicks Save; stopApp keeps the value.
    applet_seen <- new.env(parent = emptyenv())
    setBuilderInputs <- function(session, ...) {
        values <- list(...)
        names(values) <- paste0("designBuilder-", names(values))
        do.call(session$setInputs, values)
    }
    local_mocked_bindings(
        runApp = function(appDir, ...) {
            shiny::testServer(appDir, {
                session$flushReact()
                setBuilderInputs(session, sample_to_rename = "S10", new_sample_name = "T10")
                setBuilderInputs(session, rename_sample = 1L)
                for (factor_name in c("Ctl", "Trt")) {
                    setBuilderInputs(session, new_factor = factor_name)
                    setBuilderInputs(session, add_factor = which(c("Ctl", "Trt") == factor_name))
                }
                assignments <- c(S1 = "Ctl", T10 = "Trt")
                for (index in seq_along(assignments)) {
                    setBuilderInputs(session, selected_runs = names(assignments)[[index]],
                        factor1_select = assignments[[index]], factor2_select = "", factor3_select = "",
                        replicate_start = 1)
                    setBuilderInputs(session, assign_metadata = index)
                }
                setBuilderInputs(session, contrast_group1 = "Trt", contrast_group2 = "Ctl")
                setBuilderInputs(session, add_contrast = 1L)
                setBuilderInputs(session, formula_string = "~ 0 + group")
                setBuilderInputs(session, save_results = 1L)
            })
            applet_seen$stopped
        },
        stopApp = function(returnValue = invisible()) {
            applet_seen$stopped <- returnValue
            invisible(NULL)
        },
        updateSelectizeInput = function(session, inputId, label = NULL, choices = NULL, selected = NULL, ...) {
            invisible(NULL)
        },
        showNotification = function(ui, ...) invisible(NULL),
        .package = "shiny"
    )

    for (force in c(FALSE, TRUE)) {
        applet_seen$stopped <- NULL
        project_dirs <- kit231SetupDirs("proteomics", kit231Label)
        source_dir <- project_dirs[[paste0("proteomics_", kit231Label)]]$source_dir
        data_tbl <- vroom::vroom(report_path, show_col_types = FALSE)
        config_list <- list(deAnalysisParameters = list(formula_string = "~ 0 + 1"))
        error <- tryCatch({
            kit231Quiet(RunApplet("designMatrix", "proteomics", kit231Label, force = force))
            NULL
        }, error = function(e) conditionMessage(e))
        expect_null(error, label = paste("RunApplet error, force =", force))
        expect_identical(names(applet_seen$stopped), c("design_matrix", "data_cln", "contrasts_tbl", "config_list"),
            label = paste("saved result, force =", force))
        expect_true(file.exists(file.path(source_dir, "design_matrix_results.rds")),
            label = paste("design_matrix_results.rds in source_dir, force =", force))
    }
})
