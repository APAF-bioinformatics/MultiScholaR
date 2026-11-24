# KIT-231 gate G0 (ticket section 7): what R does with the shipped workbooks' folder reads, on a runner, at the base
# commit plus this script. Records only; no assertion, no expected value is coded here.
#   P  primitives: R's behaviour on zero-length paths and the calls the workbooks make with a NULL folder element;
#   S  static scan of all eight workbooks under inst/workbooks/: every call to a MultiScholaR export matched against its
#      formals (the generic's for an S4 generic), every called function looked up, every project_dirs element read
#      evaluated against what setupDirectories() returns for the workbook's own set-up, and per chunk the symbols read
#      that no earlier chunk assigns (recorded, not judged);
#   C  the affected chunks, evaluated on a real ProteinQuantitativeData (KIT-200 section 4) and a real
#      protein_deAnalysisWrapperFunction result; only writeInteractiveVolcanoPlotProteomics is mocked;
#   F  RunApplet("designMatrix", ...) with shiny::runApp mocked, each workbook's design chunk, and the experienced
#      workbook's 'Raw Data QC' call.
# Fixture builders copied (not called) from tests/testthat/test-kit200-prot-da-no-value-proteins.R on branch
# origin/kit-200-prot-da-proteins-without-values (kit200WithSeed, kit200Matrix, kit200Design, kit200Object).
# Usage (run-tests scripts job, from the package root after R CMD INSTALL): Rscript g0.R <out dir>

options(warn = 1, cli.num_colors = 1, crayon.enabled = FALSE, width = 200)
Sys.setenv(NO_COLOR = "1")

g0Args <- commandArgs(trailingOnly = TRUE)
g0Out <- if (length(g0Args) >= 1L) g0Args[[1L]] else "docs/qa/workflow-execution/kit/kit-231/g0-out"
dir.create(g0Out, recursive = TRUE, showWarnings = FALSE)
g0Out <- normalizePath(g0Out)
g0Root <- normalizePath(getwd())
g0Scratch <- tempfile("kit231-scratch-")
dir.create(g0Scratch)
g0Scratch <- normalizePath(g0Scratch)
g0Wd <- file.path(g0Scratch, "wd")
dir.create(g0Wd)
g0Wd <- normalizePath(g0Wd)

g0Log <- function(...) {
    cat(format(Sys.time(), "%H:%M:%S"), "KIT-231 G0:", ..., "\n")
    flush(stdout())
}

# ---- Recording ---------------------------------------------------------------------------------------------------

g0Clean <- function(x, limit = 3000L) {
    x <- paste(as.character(x), collapse = " | ")
    x <- gsub("[\t\r\n]+", " ", x)
    if (nchar(x) > limit) x <- paste0(substr(x, 1L, limit), " ...[truncated]")
    x
}

g0Value <- function(v) {
    if (is.null(v)) return("NULL")
    if (is.atomic(v) && length(v) <= 6L) return(paste(deparse(v, width.cutoff = 500L), collapse = ""))
    paste0("<", paste(class(v), collapse = ","), " length ", length(v), ">")
}

g0Fields <- c("part", "file", "line", "chunk", "expression", "line_text", "value", "class", "condition_class",
    "condition_message", "condition_call", "messages", "warnings", "files_created", "note")
g0Rows <- list()
g0Add <- function(...) {
    row <- list(...)
    stopifnot(all(names(row) %in% g0Fields))
    full <- stats::setNames(as.list(rep("", length(g0Fields))), g0Fields)
    for (n in names(row)) full[[n]] <- g0Clean(row[[n]])
    g0Rows[[length(g0Rows) + 1L]] <<- full
    invisible(NULL)
}

g0Write <- function(rows, fields, name) {
    if (!length(rows)) {
        df <- as.data.frame(stats::setNames(replicate(length(fields), character(), simplify = FALSE), fields))
    } else {
        df <- do.call(rbind, lapply(rows, function(r) as.data.frame(r[fields], stringsAsFactors = FALSE)))
    }
    df[] <- lapply(df, function(col) gsub("[\t\r\n]+", " ", as.character(col)))
    utils::write.table(df, file.path(g0Out, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

g0Flush <- function() {
    g0Write(g0Rows, g0Fields, "g0.tsv")
    g0Write(g0Stale, g0StaleFields, "g0-stale.tsv")
    g0Write(g0Calls, g0CallFields, "g0-calls.tsv")
    g0Write(g0DirReads, g0DirFields, "g0-dirs.tsv")
    g0Write(g0Symbols, g0SymbolFields, "g0-symbols.tsv")
}

# Evaluate one expression; every condition is recorded, the script continues.
g0Try <- function(expr, env) {
    msgs <- character()
    warns <- character()
    res <- NULL
    out <- utils::capture.output({
        res <- tryCatch(
            withCallingHandlers(
                list(value = eval(expr, env), ok = TRUE, cond_class = "", cond_msg = "", call = ""),
                message = function(m) {
                    msgs <<- c(msgs, conditionMessage(m))
                    invokeRestart("muffleMessage")
                },
                warning = function(w) {
                    warns <<- c(warns, conditionMessage(w))
                    invokeRestart("muffleWarning")
                }
            ),
            error = function(e) list(value = NULL, ok = FALSE, cond_class = paste(class(e), collapse = ","),
                cond_msg = conditionMessage(e),
                call = paste(deparse(conditionCall(e), width.cutoff = 200L), collapse = " "))
        )
    }, type = "output")
    res$messages <- msgs
    res$warnings <- warns
    res$output <- out
    res
}

g0Snap <- function(dirs) {
    unlist(lapply(dirs, function(d) {
        if (!dir.exists(d)) return(character())
        file.path(d, list.files(d, recursive = TRUE, all.files = TRUE, include.dirs = TRUE, no.. = TRUE))
    }), use.names = FALSE)
}

g0Rel <- function(paths, root) {
    paths <- gsub(g0Wd, "<wd>", paths, fixed = TRUE)
    if (!identical(root, g0Scratch)) paths <- gsub(root, "<project>", paths, fixed = TRUE)
    gsub(g0Scratch, "<scratch>", paths, fixed = TRUE)
}

g0Created <- function(before, after, root) {
    new <- setdiff(after, before)
    if (!length(new)) return("none")
    paste(g0Rel(sort(new), root), collapse = "; ")
}

# Evaluate chunk code expression by expression, as knitr does; stops at the first error and reports its line.
g0EvalCode <- function(code, start_line, env) {
    exprs <- tryCatch(parse(text = code, keep.source = TRUE), error = function(e) e)
    if (inherits(exprs, "error")) {
        return(list(ok = FALSE, stop_line = NA, cond_class = "parse_error", cond_msg = conditionMessage(exprs),
            call = "", messages = character(), warnings = character(), done = 0L, total = NA, value = NULL))
    }
    srcrefs <- attr(exprs, "srcref")
    msgs <- character()
    warns <- character()
    value <- NULL
    for (i in seq_along(exprs)) {
        r <- g0Try(exprs[[i]], env)
        msgs <- c(msgs, r$messages)
        warns <- c(warns, r$warnings)
        if (!r$ok) {
            line <- start_line - 1L + as.integer(srcrefs[[i]])[[1L]]
            return(list(ok = FALSE, stop_line = line, cond_class = r$cond_class, cond_msg = r$cond_msg,
                call = r$call, messages = msgs, warnings = warns, done = i - 1L, total = length(exprs), value = NULL))
        }
        value <- r$value
    }
    list(ok = TRUE, stop_line = NA, cond_class = "", cond_msg = "", call = "", messages = msgs, warnings = warns,
        done = length(exprs), total = length(exprs), value = value)
}

# ---- Workbook chunks ---------------------------------------------------------------------------------------------

g0Workbooks <- c(
    "inst/workbooks/proteomics/DIA_workflow_starter.rmd",
    "inst/workbooks/proteomics/DIA_workflow_limpa.rmd",
    "inst/workbooks/proteomics/DIA_workflow_limpa_starter.rmd",
    "inst/workbooks/proteomics/DIA_workflow_limpa_experienced.rmd",
    "inst/workbooks/proteomics/DIA_workflow_experienced.rmd",
    "inst/workbooks/metabolomics/metabolomics_workflow_starter.rmd",
    "inst/workbooks/metabolomics/metabolomics_workflow_experienced.rmd",
    "inst/workbooks/integration/integration_workflow_starter.rmd"
)

g0Chunks <- function(path) {
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
            code <- if (j > i + 1L) lines[(i + 1L):(j - 1L)] else character()
            chunks[[length(chunks) + 1L]] <- list(label = label, header = i, code_start = i + 1L, code = code,
                eval_false = grepl("eval\\s*=\\s*(FALSE|F)\\b", inner))
            i <- j + 1L
        } else {
            i <- i + 1L
        }
    }
    chunks
}

g0ChunkByLabel <- function(chunks, label) {
    hit <- Filter(function(ch) identical(ch$label, label), chunks)
    if (length(hit)) hit[[1L]] else NULL
}

g0Literal <- function(code, name) {
    m <- regmatches(code, regexec(paste0("^\\s*", name, "\\s*<-\\s*\"([^\"]*)\""), code))
    m <- Filter(length, m)
    if (length(m)) m[[1L]][[2L]] else NA_character_
}

g0Target <- function(t) {
    if (is.name(t)) return(as.character(t))
    if (is.character(t)) return(t)
    if (is.call(t) && length(t) >= 2L) return(g0Target(t[[2L]]))
    character()
}

# Names a piece of code assigns (`<-`, `=`, `<<-`, `->`, assign("x", ...), for-loop variables), at any depth.
g0Assigned <- function(exprs) {
    out <- character()
    walk <- function(x) {
        if (!is.call(x)) return(invisible(NULL))
        head <- x[[1L]]
        if (is.name(head)) {
            h <- as.character(head)
            if (h %in% c("<-", "=", "<<-") && length(x) >= 3L) out <<- c(out, g0Target(x[[2L]]))
            if (identical(h, "assign") && length(x) >= 2L && is.character(x[[2L]])) out <<- c(out, x[[2L]])
            if (identical(h, "for") && length(x) >= 2L && is.name(x[[2L]])) out <<- c(out, as.character(x[[2L]]))
        }
        xs <- as.list(x)
        for (k in seq_along(xs)) {
            if (identical(xs[[k]], quote(expr = ))) next
            walk(xs[[k]])
        }
        invisible(NULL)
    }
    for (k in seq_along(exprs)) walk(exprs[[k]])
    unique(out)
}

g0Kids <- function(pd, id) {
    k <- pd[pd$parent == id, , drop = FALSE]
    k[order(k$line1, k$col1), , drop = FALSE]
}

# ---- Part P ------------------------------------------------------------------------------------------------------

g0PartP <- function() {
    g0Log("part P: primitives")
    env <- new.env(parent = globalenv())
    env$p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) + ggplot2::geom_point()
    env$png_path <- file.path(g0Scratch, "kit231-null-plot.png")
    codes <- c(
        "list(a = 1)$de_output_dir",
        "file.path(NULL)",
        "file.path(NULL, \"x\")",
        "file.exists(character(0))",
        "logical(0) && TRUE",
        "logical(0) && FALSE",
        "if (file.exists(character(0)) && !FALSE) 1",
        "if (file.exists(character(0)) && !TRUE) 1",
        "dir.create(character(0))",
        "ggplot2::ggsave(character(0), p)",
        "ggplot2::ggsave(png_path, plot = NULL)",
        "grDevices::png(character(0))",
        "grDevices::pdf(character(0))",
        "saveRDS(1, character(0))",
        "vroom::vroom_write(data.frame(a = 1), character(0))",
        "writexl::write_xlsx(data.frame(a = 1), character(0))",
        "createDAResultsForEnrichment(contrasts_tbl = NULL, design_matrix = NULL, de_output_dir = stop(\"evaluated\"))"
    )
    for (code in codes) {
        devices_before <- grDevices::dev.list()
        before <- g0Snap(c(g0Wd, g0Scratch))
        r <- g0Try(parse(text = code)[[1L]], env)
        devices_after <- grDevices::dev.list()
        opened <- devices_after[!devices_after %in% devices_before]
        note <- ""
        if (length(opened)) {
            note <- paste0("device opened: ", paste(names(opened), opened, collapse = ", "), "; closed with dev.off()")
            for (d in opened) grDevices::dev.off(d)
        }
        after <- g0Snap(c(g0Wd, g0Scratch))
        g0Add(part = "P", expression = code, value = if (r$ok) g0Value(r$value) else "",
            class = if (r$ok) paste(class(r$value), collapse = ",") else "", condition_class = r$cond_class,
            condition_message = r$cond_msg, condition_call = r$call, messages = r$messages, warnings = r$warnings,
            files_created = g0Created(before, after, g0Scratch), note = note)
    }

    # fs::dir_copy onto an existing destination that holds an edited copy (the application's call, R/app_server.R:544)
    src <- system.file("workbooks", "proteomics", package = "MultiScholaR")
    dst <- file.path(g0Scratch, "kit231-dircopy", "scripts", "proteomics")
    dir.create(dst, recursive = TRUE)
    edited <- file.path(dst, "DIA_workflow_starter.rmd")
    writeLines("USER EDIT kit231", edited)
    writeLines("the user's own notes", file.path(dst, "kit231_user_notes.txt"))
    env$src <- src
    env$dst <- dst
    r <- g0Try(quote(fs::dir_copy(src, dst, overwrite = TRUE)), env)
    first <- if (file.exists(edited)) readLines(edited, n = 1L, warn = FALSE) else "<absent>"
    src_first <- readLines(file.path(src, "DIA_workflow_starter.rmd"), n = 1L, warn = FALSE)
    g0Add(part = "P", expression = "fs::dir_copy(<system.file('workbooks','proteomics')>, <dst>, overwrite = TRUE)",
        value = paste0("dst/DIA_workflow_starter.rmd first line now: '", first, "' (shipped first line: '", src_first,
            "'); user's own file kept: ", file.exists(file.path(dst, "kit231_user_notes.txt")),
            "; nested <dst>/proteomics exists: ", dir.exists(file.path(dst, "proteomics")),
            "; files in dst: ", length(list.files(dst, recursive = TRUE)),
            "; files in src: ", length(list.files(src, recursive = TRUE))),
        condition_class = r$cond_class, condition_message = r$cond_msg, condition_call = r$call,
        messages = r$messages, warnings = r$warnings,
        note = paste("dst listing:", paste(list.files(dst, recursive = TRUE), collapse = "; ")))
    g0Flush()
}

# ---- Fixture (KIT-200 section 4; builders copied from the KIT-200 test file) -----------------------------------

kit200WithSeed <- function(seed, code) withr::with_seed(seed, code, .rng_kind = "Mersenne-Twister",
    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

kit200Matrix <- function(sizes = c(A = 4L, B = 4L, C = 4L), seed = 200001L, proteins = 60L) {
    groups <- rep(names(sizes), sizes)
    x <- kit200WithSeed(seed, {
        level <- stats::rnorm(proteins, 22, 2)
        effect <- as.numeric(seq_len(proteins) %in% 11:16)
        level + outer(effect, as.numeric(groups == "B")) +
            matrix(stats::rnorm(proteins * length(groups), sd = 0.3), proteins)
    })
    dimnames(x) <- list(paste0("P", seq_len(proteins)), paste0("S", seq_along(groups)))
    x
}

kit200Design <- function(sizes = c(A = 4L, B = 4L, C = 4L)) {
    runs <- paste0("S", seq_len(sum(sizes)))
    data.frame(Run = runs, group = rep(names(sizes), sizes), replicates = runs, row.names = runs,
        stringsAsFactors = FALSE)
}

kit200Object <- function(x, design) {
    # The KIT-200 tests run under load_all(export_all = TRUE); the installed package does not export the generator.
    # Instrument fix after run 45509141744 (VAL-54): resolved here, locally, so the search path S scans is unchanged.
    ProteinQuantitativeData <- get("ProteinQuantitativeData", envir = asNamespace("MultiScholaR"))
    rownames(design) <- NULL
    ProteinQuantitativeData(protein_quant_table = data.frame(Protein.Ids = rownames(x), x, row.names = NULL,
        check.names = FALSE, stringsAsFactors = FALSE), design_matrix = design, sample_id = "Run", group_id = "group",
        technical_replicate_id = "replicates", protein_id_column = "Protein.Ids",
        protein_id_table = data.frame(Protein.Ids = rownames(x), gene_names = paste0("GENE_", rownames(x)),
            stringsAsFactors = FALSE), args = list())
}

# ---- Part C ------------------------------------------------------------------------------------------------------

g0State <- new.env()
g0State$volcano <- 0L
g0State$runapp <- 0L

g0WithVolcanoMock <- function(code) {
    testthat::with_mocked_bindings(code,
        writeInteractiveVolcanoPlotProteomics = function(...) {
            g0State$volcano <- g0State$volcano + 1L
            invisible(NULL)
        },
        .package = "MultiScholaR")
}

g0WithRunAppMock <- function(code) {
    testthat::with_mocked_bindings(code,
        runApp = function(...) {
            g0State$runapp <- g0State$runapp + 1L
            list()
        },
        .package = "shiny")
}

g0Setup <- function(omic_types, label) {
    base <- tempfile("kit231-project-", tmpdir = g0Scratch)
    dir.create(base)
    base <- normalizePath(base)
    r <- g0Try(substitute(setupDirectories(base_dir = B, omic_types = O, label = L, force = TRUE),
        list(B = base, O = omic_types, L = label)), globalenv())
    list(base = base, dirs = r$value, ok = r$ok, msg = r$cond_msg)
}

g0LineText <- function(file, line) {
    lines <- readLines(file.path(g0Root, file), warn = FALSE)
    if (is.na(line) || line > length(lines)) return("")
    trimws(lines[[line]])
}

g0CPlan <- list(
    list(file = "inst/workbooks/proteomics/DIA_workflow_starter.rmd",
        folders = list(list(1725L, "U1", "file.path(project_dirs[[paste0(omic_type, \"_\", experiment_label)]]$de_output_dir)"),
            list(1768L, "U2", "file.path(project_dirs[[paste0(omic_type, \"_\", experiment_label)]]$de_output_dir)")),
        before = character(), out = "34 Output DE Analysis Results", enr = "35 Enrichment Analysis", after = character()),
    list(file = "inst/workbooks/proteomics/DIA_workflow_limpa.rmd",
        folders = list(list(2021L, "U3", "file.path(project_dirs[[paste0(omic_type, \"_\", experiment_label)]]$de_output_dir)"),
            list(2064L, "U3", "file.path(project_dirs[[paste0(omic_type, \"_\", experiment_label)]]$de_output_dir)")),
        before = character(), out = "34 Output DE Analysis Results", enr = "35 Enrichment Analysis", after = character()),
    list(file = "inst/workbooks/proteomics/DIA_workflow_limpa_starter.rmd",
        folders = list(list(1170L, "U5", "project_dirs[[paste0(\"integration_\", experiment_label)]]$multiomic_inputs_dir"),
            list(1333L, "U4", "file.path(project_dirs[[paste0(omic_type, \"_\", experiment_label)]]$de_output_dir)"),
            list(1378L, "U4", "file.path(project_dirs[[paste0(omic_type, \"_\", experiment_label)]]$de_output_dir)")),
        before = "17 Output Files For Audit Trail", out = "20 Output DE Analysis Results", enr = "21 Enrichment Analysis",
        after = character()),
    list(file = "inst/workbooks/proteomics/DIA_workflow_limpa_experienced.rmd",
        folders = list(list(1470L, "U6", "file.path(project_dirs[[paste0(omic_type, \"_\", experiment_label)]]$de_output_dir)"),
            list(1563L, "U7", "file.path( project_dirs[[paste0(omic_type, \"_\", experiment_label)]]$de_output_dir, \"top_de_genes.tsv\")"),
            list(1663L, "U6", "file.path(project_dirs[[paste0(omic_type, \"_\", experiment_label)]]$de_output_dir)")),
        before = character(), out = "20 Output DE Analysis Results", enr = "21 Enrichment Analysis",
        after = "extract-top-de-genes"),
    list(file = "inst/workbooks/proteomics/DIA_workflow_experienced.rmd",
        folders = list(list(1147L, "U8", "de_output_dir"), list(1148L, "U8", "publication_graphs_dir"),
            list(1175L, "U8", "de_output_dir"), list(1186L, "U8", "pathway_dir")),
        before = character(), out = "Output DE Analysis Results", enr = "Enrichment Analysis", after = character())
)

g0SetupSpec <- function(file) {
    chunks <- g0Chunks(file.path(g0Root, file))
    code <- unlist(lapply(chunks, `[[`, "code"))
    label <- g0Literal(code, "experiment_label")
    omic <- g0Literal(code, "omic_type")
    env <- new.env(parent = globalenv())
    if (!is.na(label)) env$experiment_label <- label
    spec <- NULL
    for (ch in chunks) {
        exprs <- tryCatch(parse(text = ch$code, keep.source = TRUE), error = function(e) NULL)
        if (is.null(exprs) || !length(exprs)) next
        pd <- utils::getParseData(exprs)
        hits <- pd$id[pd$token == "SYMBOL_FUNCTION_CALL" & pd$text == "setupDirectories"]
        if (!length(hits)) next
        call_id <- pd$parent[pd$id == pd$parent[pd$id == hits[[1L]]]]
        cl <- str2lang(utils::getParseText(pd, call_id))
        a <- as.list(cl)[-1L]
        spec <- list(omic_types = eval(a[["omic_types"]], env), label = eval(a[["label"]], env),
            text = g0Clean(utils::getParseText(pd, call_id)), line = ch$code_start - 1L + pd$line1[pd$id == call_id])
        break
    }
    list(chunks = chunks, label = label, omic = omic, spec = spec)
}

g0CRow <- function(file, line, chunk, expression, r, before, after, root, note = "") {
    g0Add(part = "C", file = basename(file), line = line, chunk = chunk, expression = expression,
        line_text = if (!is.na(suppressWarnings(as.integer(line)))) g0LineText(file, as.integer(line)) else "",
        value = if (isTRUE(r$ok)) g0Value(r$value) else "",
        class = if (isTRUE(r$ok)) paste(class(r$value), collapse = ",") else "",
        condition_class = r$cond_class, condition_message = r$cond_msg, condition_call = r$call,
        messages = r$messages, warnings = r$warnings, files_created = g0Created(before, after, root), note = note)
}

g0RunChunk <- function(file, chunks, label, env, root, note = "") {
    ch <- g0ChunkByLabel(chunks, label)
    if (is.null(ch)) {
        g0Add(part = "C", file = basename(file), chunk = label, note = "chunk label not found")
        return(invisible(NULL))
    }
    before <- g0Snap(c(root, g0Wd))
    volcano_before <- g0State$volcano
    r <- g0EvalCode(ch$code, ch$code_start, env)
    after <- g0Snap(c(root, g0Wd))
    g0CRow(file, if (is.na(r$stop_line)) paste0(ch$code_start, "-", ch$code_start + length(ch$code) - 1L) else r$stop_line,
        label, paste0("whole chunk (", r$done, " of ", r$total, " top-level expressions completed)"), r, before, after,
        root, note = paste0(note, if (nzchar(note)) "; " else "", "status: ", if (r$ok) "completed" else "stopped",
            "; volcano writer (mocked) calls: ", g0State$volcano - volcano_before))
    invisible(r)
}

g0PartC <- function(obj, da_list, contrasts, uniprot, config) {
    for (plan in g0CPlan) {
        file <- plan$file
        g0Log("part C:", basename(file))
        info <- g0SetupSpec(file)
        setup <- g0Setup(info$spec$omic_types, info$spec$label)
        root <- setup$base
        env <- new.env(parent = globalenv())
        env$project_dirs <- setup$dirs
        if (!is.na(info$label)) env$experiment_label <- info$label
        if (!is.na(info$omic)) env$omic_type <- info$omic
        env$config_list <- config
        env$ruv_normalised_for_da_analysis_obj <- obj
        env$ruv_normalised_filtered_results_obj <- obj
        env$uniprot_dat_cln <- uniprot
        env$da_analysis_results_list <- da_list
        env$contrasts_tbl <- contrasts
        env$design_matrix <- obj@design_matrix
        g0Add(part = "C", file = basename(file), line = info$spec$line, chunk = "set-up",
            expression = paste0("setupDirectories(base_dir = <project>, omic_types = ", deparse(info$spec$omic_types),
                ", label = '", info$spec$label, "', force = TRUE)"),
            value = paste0("keys: ", paste(names(setup$dirs), collapse = ", ")),
            condition_message = setup$msg,
            note = paste0("workbook's experiment_label: ", info$label, "; workbook's omic_type assignment: ",
                if (is.na(info$omic)) "none" else info$omic, "; workbook call: ", info$spec$text))
        for (f in plan$folders) {
            before <- g0Snap(c(root, g0Wd))
            r <- g0Try(str2lang(f[[3L]]), env)
            after <- g0Snap(c(root, g0Wd))
            g0CRow(file, f[[1L]], f[[2L]], f[[3L]], r, before, after, root, note = "folder expression")
        }
        for (label in plan$before) g0RunChunk(file, info$chunks, label, env, root)
        g0RunChunk(file, info$chunks, plan$out, env, root, note = "DA output chunk")
        # The enrichment chunk's createDAResultsForEnrichment(...) call (its first top-level expression)
        ch <- g0ChunkByLabel(info$chunks, plan$enr)
        if (is.null(ch)) {
            g0Add(part = "C", file = basename(file), chunk = plan$enr, note = "chunk label not found")
            next
        }
        exprs <- parse(text = ch$code, keep.source = TRUE)
        first <- exprs[[1L]]
        line <- ch$code_start - 1L + as.integer(attr(exprs, "srcref")[[1L]])[[1L]]
        before <- g0Snap(c(root, g0Wd))
        r <- g0Try(first, env)
        after <- g0Snap(c(root, g0Wd))
        note <- "enrichment chunk, first expression only (processEnrichments not run)"
        if (r$ok && isS4(env$de_results_for_enrichment)) {
            dd <- env$de_results_for_enrichment@da_data
            note <- paste0(note, "; da_data: ", length(dd), " element(s), non-NULL: ",
                sum(!vapply(dd, is.null, logical(1))), "; rows: ",
                paste(vapply(dd, function(t) if (is.null(t)) "NULL" else as.character(nrow(t)), character(1)), collapse = ","))
        }
        g0CRow(file, line, plan$enr, paste(deparse(first, width.cutoff = 200L), collapse = " "), r, before, after,
            root, note = note)
        for (label in plan$after) {
            g0RunChunk(file, info$chunks, label, env, root)
            if (identical(label, "extract-top-de-genes")) {
                # separately, :1563-1564, with top_de_genes as the chunk would have built it from de_proteins_long
                lines <- readLines(file.path(g0Root, file), warn = FALSE)[1563:1564]
                long <- da_list[[1L]]$de_proteins_long
                env$top_de_genes <- tryCatch(dplyr::filter(long, fdr_qvalue < 0.05), error = function(e) long)
                before <- g0Snap(c(root, g0Wd))
                r <- g0EvalCode(lines, 1563L, env)
                after <- g0Snap(c(root, g0Wd))
                g0CRow(file, if (is.na(r$stop_line)) "1563-1564" else r$stop_line, label,
                    paste(trimws(lines), collapse = " ; "), r, before, after, root,
                    note = paste0("lines 1563-1564 alone; top_de_genes set to de_proteins_long rows with fdr_qvalue < 0.05 (",
                        nrow(env$top_de_genes), " rows); status: ", if (r$ok) "completed" else "stopped"))
            }
        }
        g0Add(part = "C", file = basename(file), chunk = "project after C",
            files_created = paste(g0Rel(sort(g0Snap(root)[grepl("results", g0Snap(root))]), root), collapse = "; "),
            note = "every file and folder under <project>/results* at the end of this workbook's C rows")
        g0Flush()
    }
}

# ---- Part F ------------------------------------------------------------------------------------------------------

g0FEnv <- function(dirs, label) {
    env <- new.env(parent = globalenv())
    env$project_dirs <- dirs
    if (!is.na(label)) env$experiment_label <- label
    env$data_tbl <- data.frame(Run = paste0("S", 1:12), Precursor.Normalised = 1, stringsAsFactors = FALSE)
    env$config_list <- g0Config
    env
}

g0FRow <- function(file, line, chunk, expression, r, before, after, root, runapp_before, note = "") {
    g0Add(part = "F", file = file, line = line, chunk = chunk, expression = expression,
        value = if (isTRUE(r$ok)) g0Value(r$value) else "", condition_class = r$cond_class,
        condition_message = r$cond_msg, condition_call = r$call, messages = r$messages, warnings = r$warnings,
        files_created = g0Created(before, after, root),
        note = paste0(note, if (nzchar(note)) "; " else "", "shiny::runApp (mocked) reached: ",
            g0State$runapp - runapp_before, " time(s)"))
}

g0PartF <- function() {
    g0Log("part F: RunApplet")
    setup <- g0Setup("proteomics", "kit231")
    env <- g0FEnv(setup$dirs, "kit231")
    g0Add(part = "F", chunk = "set-up", expression = "setupDirectories(omic_types = 'proteomics', label = 'kit231')",
        value = paste0("element output_dir present: ", "output_dir" %in% names(setup$dirs[["proteomics_kit231"]]),
            "; source_dir: ", g0Rel(setup$dirs[["proteomics_kit231"]]$source_dir, setup$base)))
    for (force in c(FALSE, TRUE)) {
        code <- sprintf("RunApplet(\"designMatrix\", \"proteomics\", \"kit231\", force = %s)", force)
        before <- g0Snap(c(setup$base, g0Wd))
        n <- g0State$runapp
        r <- g0WithRunAppMock(g0Try(str2lang(code), env))
        after <- g0Snap(c(setup$base, g0Wd))
        g0FRow("", "", "RunApplet", code, r, before, after, setup$base, n)
    }
    # What :387 does with the output file RunApplet computes (the server never runs while runApp is mocked)
    code <- "saveRDS(list(design_matrix = data.frame()), file.path(project_dirs[['proteomics_kit231']]$output_dir, 'design_matrix_results.rds'))"
    before <- g0Snap(c(setup$base, g0Wd))
    r <- g0Try(str2lang(code), env)
    after <- g0Snap(c(setup$base, g0Wd))
    g0FRow("", "", "writeDesignAppletResult :387 equivalent", code, r, before, after, setup$base, g0State$runapp,
        note = "the saveRDS of R/func_general_applets.R:387 on output_file = file.path(<output_dir element>, ...)")

    # Each workbook's design chunk (the chunk holding an uncommented RunApplet call), runApp mocked
    for (file in g0Workbooks) {
        info <- g0SetupSpec(file)
        design <- NULL
        for (ch in info$chunks) {
            exprs <- tryCatch(parse(text = ch$code, keep.source = TRUE), error = function(e) NULL)
            if (is.null(exprs) || !length(exprs)) next
            pd <- utils::getParseData(exprs)
            if (any(pd$token == "SYMBOL_FUNCTION_CALL" & pd$text == "RunApplet")) {
                design <- ch
                break
            }
        }
        if (is.null(design)) {
            g0Add(part = "F", file = basename(file), chunk = "design chunk", note = "no uncommented RunApplet call in any chunk")
            next
        }
        s <- if (is.null(info$spec)) list(base = g0Scratch, dirs = NULL) else g0Setup(info$spec$omic_types, info$spec$label)
        fenv <- g0FEnv(s$dirs, info$label)
        before <- g0Snap(c(s$base, g0Wd))
        n <- g0State$runapp
        r <- g0WithRunAppMock(g0EvalCode(design$code, design$code_start, fenv))
        after <- g0Snap(c(s$base, g0Wd))
        g0FRow(basename(file), if (is.na(r$stop_line)) paste0(design$code_start, "-", design$code_start + length(design$code) - 1L) else r$stop_line,
            design$label, paste0("whole design chunk (", r$done, " of ", r$total, " expressions completed)"), r,
            before, after, s$base, n, note = paste0("status: ", if (r$ok) "completed" else "stopped",
                "; design_matrix exists in .GlobalEnv: ", exists("design_matrix", envir = .GlobalEnv)))
    }

    # The experienced workbook's 'Raw Data QC' call
    file <- "inst/workbooks/proteomics/DIA_workflow_experienced.rmd"
    info <- g0SetupSpec(file)
    ch <- g0ChunkByLabel(info$chunks, "Raw Data QC")
    s <- g0Setup(info$spec$omic_types, info$spec$label)
    fenv <- g0FEnv(s$dirs, info$label)
    fenv$data_cln <- data.frame(Run = "S1", stringsAsFactors = FALSE)
    before <- g0Snap(c(s$base, g0Wd))
    r <- g0EvalCode(ch$code, ch$code_start, fenv)
    after <- g0Snap(c(s$base, g0Wd))
    g0FRow(basename(file), if (is.na(r$stop_line)) ch$code_start else r$stop_line, "Raw Data QC",
        "whole chunk (updateProteinFiltering(...)); data_cln set to a one-row data frame", r, before, after, s$base,
        g0State$runapp)
    g0Flush()
}

# ---- Part S ------------------------------------------------------------------------------------------------------

g0StaleFields <- c("file", "line", "chunk", "kind", "name", "detail", "eval_false")
g0Stale <- list()
g0CallFields <- c("file", "line", "chunk", "name", "package", "resolves_to", "piped", "status", "eval_false")
g0Calls <- list()
g0DirFields <- c("file", "line", "chunk", "expression", "key_present", "element_present", "value", "eval_false")
g0DirReads <- list()
g0SymbolFields <- c("file", "chunk", "line", "kind", "symbols", "eval_false")
g0Symbols <- list()

g0Push <- function(store, row) {
    row <- lapply(row, g0Clean)
    assign(store, c(get(store, envir = globalenv()), list(row)), envir = globalenv())
}

g0Exports <- getNamespaceExports("MultiScholaR")

g0Where <- function(name) {
    w <- utils::find(name, mode = "function")
    if (length(w)) w[[1L]] else ""
}

g0ScanCalls <- function(file, ch, pd, workbook_assigned) {
    toks <- pd[pd$token == "SYMBOL_FUNCTION_CALL", , drop = FALSE]
    for (k in seq_len(nrow(toks))) {
        name <- toks$text[[k]]
        fexpr <- toks$parent[[k]]
        sibs <- pd[pd$parent == fexpr, , drop = FALSE]
        if (any(sibs$token %in% c("'$'", "'@'"))) next
        pkg <- sibs$text[sibs$token == "SYMBOL_PACKAGE"]
        op <- sibs$text[sibs$token %in% c("NS_GET", "NS_GET_INT")]
        call_id <- pd$parent[pd$id == fexpr]
        line <- ch$code_start - 1L + pd$line1[pd$id == call_id]
        parent_id <- pd$parent[pd$id == call_id]
        kids <- g0Kids(pd, parent_id)
        expr_kids <- kids$id[kids$token == "expr"]
        piped <- any(kids$token == "PIPE" | (kids$token == "SPECIAL" & kids$text == "%>%")) &&
            length(expr_kids) && utils::tail(expr_kids, 1L) == call_id
        base_row <- list(file = basename(file), line = line, chunk = ch$label, name = name,
            package = if (length(pkg)) paste0(pkg, op) else "", piped = piped, eval_false = ch$eval_false)
        stale <- function(kind, detail) g0Push("g0Stale", list(file = basename(file), line = line, chunk = ch$label,
            kind = kind, name = if (length(pkg)) paste0(pkg, op, name) else name, detail = detail,
            eval_false = ch$eval_false))
        if (length(pkg)) {
            if (!requireNamespace(pkg, quietly = TRUE)) {
                stale("package_not_installed", paste0("package '", pkg, "' is not installed on the runner"))
                g0Push("g0Calls", c(base_row, resolves_to = "", status = "package not installed"))
                next
            }
            if (identical(op, "::") && !name %in% getNamespaceExports(pkg)) {
                stale("not_exported", paste0("'", name, "' is not exported by ", pkg,
                    if (exists(name, envir = asNamespace(pkg), inherits = FALSE)) " (it exists as an internal object)" else
                        " (no such object in the namespace)"))
                g0Push("g0Calls", c(base_row, resolves_to = "", status = "not exported"))
                next
            }
            if (!exists(name, envir = asNamespace(pkg), inherits = FALSE)) {
                stale("not_found_in_namespace", paste0("no object '", name, "' in namespace ", pkg))
                g0Push("g0Calls", c(base_row, resolves_to = "", status = "not in namespace"))
                next
            }
            fn <- get(name, envir = asNamespace(pkg))
            resolves <- paste0("namespace:", pkg)
        } else {
            if (name %in% workbook_assigned) {
                g0Push("g0Calls", c(base_row, resolves_to = "the workbook's own assignment", status = "not checked"))
                next
            }
            if (!exists(name, envir = globalenv(), mode = "function")) {
                internal <- exists(name, envir = asNamespace("MultiScholaR"), inherits = FALSE)
                stale("function_not_found", paste0("no function '", name, "' on the search path after library(MultiScholaR) and the dependencies",
                    if (internal) "; it exists as an unexported MultiScholaR object" else ""))
                g0Push("g0Calls", c(base_row, resolves_to = "", status = "not found"))
                next
            }
            fn <- get(name, envir = globalenv(), mode = "function")
            resolves <- g0Where(name)
        }
        is_ms <- (length(pkg) && identical(pkg, "MultiScholaR")) || (!length(pkg) && name %in% g0Exports)
        if (!is_ms) {
            g0Push("g0Calls", c(base_row, resolves_to = resolves, status = "not a MultiScholaR export; formals not checked"))
            next
        }
        ms_fn <- getExportedValue("MultiScholaR", name)
        if (!length(pkg) && !identical(fn, ms_fn)) {
            stale("masked", paste0("MultiScholaR exports '", name, "' but the bare call resolves to ", resolves))
        }
        if (!is.function(ms_fn) || is.primitive(ms_fn)) {
            g0Push("g0Calls", c(base_row, resolves_to = resolves, status = "export is not a closure"))
            next
        }
        fml <- formals(ms_fn)
        if ("..." %in% names(fml)) {
            g0Push("g0Calls", c(base_row, resolves_to = resolves, status = "formals include ...; any argument accepted"))
            next
        }
        cl <- tryCatch(str2lang(utils::getParseText(pd, call_id)), error = function(e) NULL)
        if (is.null(cl)) {
            stale("match_error", "call text did not parse alone")
            next
        }
        if (piped) {
            args <- as.list(cl)[-1L]
            has_placeholder <- any(vapply(seq_along(args), function(i) {
                !identical(args[[i]], quote(expr = )) && (identical(args[[i]], as.name(".")) || identical(args[[i]], as.name("_")))
            }, logical(1)))
            if (!has_placeholder) cl <- as.call(c(cl[[1L]], list(as.name("kit231_piped_lhs")), args))
        }
        dummy <- function() NULL
        formals(dummy) <- fml
        msg <- tryCatch({
            match.call(dummy, cl)
            ""
        }, error = function(e) conditionMessage(e))
        if (nzchar(msg)) {
            stale(if (grepl("unused argument", msg)) "unused_argument" else "match_error",
                paste0(msg, if (methods::is(ms_fn, "standardGeneric")) " [S4 generic's formals]" else ""))
            g0Push("g0Calls", c(base_row, resolves_to = resolves, status = paste("match.call:", msg)))
        } else {
            g0Push("g0Calls", c(base_row, resolves_to = resolves, status = "formals match"))
        }
    }
}

g0ScanDirs <- function(file, ch, pd, dirs_env) {
    toks <- pd[pd$token == "'$'", , drop = FALSE]
    for (k in seq_len(nrow(toks))) {
        p_id <- toks$parent[[k]]
        kids <- g0Kids(pd, p_id)
        lhs_id <- kids$id[kids$token == "expr"][1L]
        if (is.na(lhs_id)) next
        lhs_kids <- g0Kids(pd, lhs_id)
        lhs_first <- lhs_kids$id[lhs_kids$token == "expr"][1L]
        if (is.na(lhs_first) || !identical(utils::getParseText(pd, lhs_first), "project_dirs")) next
        text <- utils::getParseText(pd, p_id)
        lhs_text <- utils::getParseText(pd, lhs_id)
        line <- ch$code_start - 1L + pd$line1[pd$id == p_id]
        if (!exists("project_dirs", envir = dirs_env, inherits = FALSE)) {
            key <- list(ok = FALSE, value = NULL, cond_msg = "the workbook calls no setupDirectories()")
            el <- key
        } else {
            key <- g0Try(str2lang(lhs_text), dirs_env)
            el <- g0Try(str2lang(text), dirs_env)
        }
        key_present <- isTRUE(key$ok) && !is.null(key$value)
        el_present <- isTRUE(el$ok) && !is.null(el$value)
        g0Push("g0DirReads", list(file = basename(file), line = line, chunk = ch$label, expression = text,
            key_present = key_present, element_present = el_present,
            value = if (el_present) g0Rel(as.character(el$value)[1L], dirs_env$kit231_base) else if (!isTRUE(el$ok)) paste("error:", el$cond_msg) else "NULL",
            eval_false = ch$eval_false))
        if (!el_present) {
            kind <- if (!isTRUE(el$ok)) "directory_read_error" else if (!key_present) "directory_key_missing" else "directory_element_missing"
            g0Push("g0Stale", list(file = basename(file), line = line, chunk = ch$label, kind = kind, name = text,
                detail = if (!isTRUE(el$ok)) el$cond_msg else if (!key_present) paste0("key ", lhs_text, " is NULL for this workbook's set-up") else "element is NULL",
                eval_false = ch$eval_false))
        }
    }
}

g0PartS <- function() {
    for (file in g0Workbooks) {
        g0Log("part S:", basename(file))
        info <- g0SetupSpec(file)
        dirs_env <- new.env(parent = globalenv())
        if (!is.na(info$label)) dirs_env$experiment_label <- info$label
        if (!is.na(info$omic)) dirs_env$omic_type <- info$omic
        if (!is.null(info$spec)) {
            s <- g0Setup(info$spec$omic_types, info$spec$label)
            dirs_env$project_dirs <- s$dirs
            dirs_env$kit231_base <- s$base
            g0Add(part = "S", file = basename(file), line = info$spec$line, chunk = "set-up",
                expression = info$spec$text, value = paste0("keys: ", paste(names(s$dirs), collapse = ", ")),
                note = paste0("omic_type assignment: ", if (is.na(info$omic)) "none" else info$omic,
                    "; experiment_label: ", info$label, "; chunks: ", length(info$chunks)))
        } else {
            dirs_env$kit231_base <- g0Scratch
            g0Add(part = "S", file = basename(file), chunk = "set-up", note = paste0("no setupDirectories() call; chunks: ",
                length(info$chunks)))
        }
        workbook_assigned <- unique(unlist(lapply(info$chunks, function(ch) {
            e <- tryCatch(parse(text = ch$code), error = function(e) NULL)
            if (is.null(e)) character() else g0Assigned(e)
        })))
        cumulative <- character()
        for (ch in info$chunks) {
            exprs <- tryCatch(parse(text = ch$code, keep.source = TRUE), error = function(e) e)
            if (inherits(exprs, "error")) {
                g0Push("g0Stale", list(file = basename(file), line = ch$code_start, chunk = ch$label, kind = "parse_error",
                    name = "", detail = conditionMessage(exprs), eval_false = ch$eval_false))
                next
            }
            if (!length(exprs)) next
            pd <- utils::getParseData(exprs)
            tryCatch(g0ScanCalls(file, ch, pd, workbook_assigned), error = function(e)
                g0Push("g0Stale", list(file = basename(file), line = ch$code_start, chunk = ch$label, kind = "instrument_error",
                    name = "calls", detail = conditionMessage(e), eval_false = ch$eval_false)))
            tryCatch(g0ScanDirs(file, ch, pd, dirs_env), error = function(e)
                g0Push("g0Stale", list(file = basename(file), line = ch$code_start, chunk = ch$label, kind = "instrument_error",
                    name = "dirs", detail = conditionMessage(e), eval_false = ch$eval_false)))
            g <- tryCatch({
                f <- function() NULL
                body(f) <- as.call(c(as.name("{"), as.list(exprs)))
                codetools::findGlobals(f, merge = FALSE)
            }, error = function(e) list(error = conditionMessage(e)))
            if (!is.null(g$error)) {
                g0Push("g0Symbols", list(file = basename(file), chunk = ch$label, line = ch$code_start, kind = "findGlobals error",
                    symbols = g$error, eval_false = ch$eval_false))
            } else {
                on_path <- function(s) exists(s, envir = globalenv())
                vars <- setdiff(g$variables, cumulative)
                vars <- vars[!vapply(vars, on_path, logical(1))]
                funs <- setdiff(g$functions, cumulative)
                funs <- funs[!vapply(funs, function(s) exists(s, envir = globalenv(), mode = "function"), logical(1))]
                if (length(vars)) g0Push("g0Symbols", list(file = basename(file), chunk = ch$label, line = ch$code_start,
                    kind = "variable read, not assigned by an earlier chunk, not on the search path",
                    symbols = paste(sort(vars), collapse = ", "), eval_false = ch$eval_false))
                if (length(funs)) g0Push("g0Symbols", list(file = basename(file), chunk = ch$label, line = ch$code_start,
                    kind = "function called, not assigned by an earlier chunk, not on the search path",
                    symbols = paste(sort(funs), collapse = ", "), eval_false = ch$eval_false))
            }
            cumulative <- union(cumulative, g0Assigned(exprs))
        }
        g0Flush()
    }
}

# ---- Run ---------------------------------------------------------------------------------------------------------

g0Log("start; package root", g0Root, "; out", g0Out)
suppressPackageStartupMessages(library(MultiScholaR))

# loadDependencies()'s package lists, read from its body; installed ones attached in its order; nothing installed
g0Deps <- list()
for (g0St in as.list(body(MultiScholaR::loadDependencies))[-1L]) {
    if (is.call(g0St) && identical(g0St[[1L]], as.name("<-")) && is.name(g0St[[2L]]) &&
        as.character(g0St[[2L]]) %in% c("cran_packages", "bioc_packages", "github_packages")) {
        g0Deps[[as.character(g0St[[2L]])]] <- eval(g0St[[3L]], baseenv())
    }
}
g0DepNames <- unique(c(g0Deps$cran_packages, g0Deps$bioc_packages, names(g0Deps$github_packages)))
g0Attached <- character()
g0Missing <- character()
g0AttachErrors <- character()
for (g0Pkg in g0DepNames) {
    if (requireNamespace(g0Pkg, quietly = TRUE)) {
        g0Ok <- tryCatch({
            suppressPackageStartupMessages(library(g0Pkg, character.only = TRUE))
            TRUE
        }, error = function(e) {
            g0AttachErrors <<- c(g0AttachErrors, paste0(g0Pkg, ": ", conditionMessage(e)))
            FALSE
        })
        if (g0Ok) g0Attached <- c(g0Attached, g0Pkg)
    } else {
        g0Missing <- c(g0Missing, g0Pkg)
        if (identical(g0Pkg, "tidyverse")) {
            for (g0Member in c("ggplot2", "tibble", "tidyr", "readr", "purrr", "dplyr", "stringr", "forcats", "lubridate")) {
                if (requireNamespace(g0Member, quietly = TRUE)) {
                    suppressPackageStartupMessages(library(g0Member, character.only = TRUE))
                    g0Attached <- c(g0Attached, paste0(g0Member, " (tidyverse member)"))
                }
            }
        }
    }
}
g0Log("attached", length(g0Attached), "; not installed", length(g0Missing))

g0Commit <- tryCatch(system2("git", c("-C", shQuote(g0Root), "rev-parse", "HEAD"), stdout = TRUE), error = function(e) "unknown")
writeLines(c(
    paste("commit:", g0Commit),
    paste("R:", R.version.string),
    paste("MultiScholaR:", as.character(utils::packageVersion("MultiScholaR"))),
    paste("versions:", paste(vapply(c("ggplot2", "vroom", "writexl", "fs", "shiny", "testthat", "dplyr", "codetools", "limma"),
        function(p) paste0(p, " ", tryCatch(as.character(utils::packageVersion(p)), error = function(e) "absent")),
        character(1)), collapse = "; ")),
    paste("loadDependencies packages attached (", length(g0Attached), "):", paste(g0Attached, collapse = ", ")),
    paste("loadDependencies packages not installed on the runner, not installed by G0 (", length(g0Missing), "):",
        paste(g0Missing, collapse = ", ")),
    paste("attach errors:", paste(g0AttachErrors, collapse = " | ")),
    paste("search path:", paste(search(), collapse = ", "))
), file.path(g0Out, "g0-session.txt"))

g0Stage <- function(name, code) {
    tryCatch(code, error = function(e) {
        g0Log("instrument error in", name, ":", conditionMessage(e))
        g0Add(part = name, note = "instrument error", condition_message = conditionMessage(e),
            condition_call = paste(deparse(conditionCall(e), width.cutoff = 200L), collapse = " "))
        g0Flush()
    })
}

setwd(g0Wd)
g0Stage("P", g0PartP())

g0Log("fixture and DA result")
g0Config <- tryCatch(readConfigFile(file = system.file("config", "config.ini", package = "MultiScholaR"), create_workers = FALSE),
    error = function(e) {
        g0Add(part = "C0", expression = "readConfigFile(<shipped config.ini>)", condition_message = conditionMessage(e),
            note = "fallback config_list with plots_format = c('pdf', 'png')")
        list(deAnalysisParameters = list(plots_format = c("pdf", "png")))
    })
g0Obj <- kit200Object(kit200Matrix(), kit200Design())
g0Contrasts <- tibble::tibble(contrasts = "B_vs_A=groupB-groupA")
g0Uniprot <- data.frame(Entry = paste0("P", 1:60), gene_names = paste0("GENE_P", 1:60), stringsAsFactors = FALSE)
g0Da <- g0Try(quote(list(B_vs_A = protein_deAnalysisWrapperFunction(g0Obj, dplyr::slice(g0Contrasts, 1L),
    formula_string = "~ 0 + group", de_q_val_thresh = 0.05, treat_lfc_cutoff = 0, eBayes_trend = TRUE,
    eBayes_robust = TRUE, args_group_pattern = "(\\d+)", args_row_id = g0Obj@protein_id_column))), globalenv())
g0Add(part = "C0", expression = "protein_deAnalysisWrapperFunction(<KIT-200 object>, <B_vs_A=groupB-groupA>, ...)",
    value = if (g0Da$ok) paste0("elements: ", paste(names(g0Da$value[[1L]]), collapse = ", "),
        "; de_proteins_long rows: ", NROW(g0Da$value[[1L]]$de_proteins_long),
        "; significant (fdr_qvalue < 0.05) in de_proteins_long: ",
        sum(g0Da$value[[1L]]$de_proteins_long$fdr_qvalue < 0.05, na.rm = TRUE),
        "; config plots_format: ", paste(g0Config$deAnalysisParameters$plots_format, collapse = ",")) else "",
    condition_class = g0Da$cond_class, condition_message = g0Da$cond_msg, condition_call = g0Da$call,
    warnings = g0Da$warnings, note = paste0("object: ", nrow(g0Obj@protein_quant_table), " proteins x ",
        nrow(g0Obj@design_matrix), " samples; groups ", paste(unique(g0Obj@design_matrix$group), collapse = ",")))
g0Flush()

if (isTRUE(g0Da$ok)) {
    g0Stage("C", g0WithVolcanoMock(g0PartC(g0Obj, g0Da$value, g0Contrasts, g0Uniprot, g0Config)))
} else {
    g0Add(part = "C", note = "not run: the DA result could not be built (instrument defect)")
}
g0Stage("F", g0PartF())
g0Stage("S", g0PartS())
g0Flush()
writeLines(c(format(Sys.time()), paste("rows:", length(g0Rows), "; stale:", length(g0Stale), "; calls:", length(g0Calls),
    "; directory reads:", length(g0DirReads))), file.path(g0Out, "g0-done.txt"))
g0Log("done")
