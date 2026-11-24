# KIT-213 gate G0 (ticket section 10; decides stop rules 4 and 5 before R1). Evidence code under docs/qa: no package
# code. Runs on a GitHub-hosted runner (run-tests `scripts`) from the package root of the commit under test, which holds
# the base's R/ unchanged.
#
# Usage: Rscript docs/qa/workflow-execution/kit/kit-213/g0.R <out dir> journeys|w1
#
# journeys: for each omic's journey of ticket section 4, the analysis twice, each in its own fresh project root, then
#   the Summary's Save (KIT-207's `k207Save`): lipidomics and metabolomics as KIT-207 G1 (run, Apply, DA); trx as KIT-207
#   G5 with RUV automatic, then the RUVIII export (`group`) when a k >= 1 was chosen; proteomics on the KIT-107 F-16
#   fixture (test-workflow-analysis-record.R) with the in-process correlation filter (`computeDIAReferenceCorrelation`,
#   committed as `runDIAReferenceCorrelation` commits it), the DA, its publication and the enrichment (test-mode kernel).
#   Writes per omic and project: every line of study_parameters.txt holding the store's incarnation or an entity id of
#   the store (its index, and every `entity_id` under the head's states, products, exports and the DA manifest); every
#   line E1's masks change (drop: `Timestamp:`, `^  timestamp = `, the version lines; the project root masked, with the
#   count per line); and every line differing between the two projects' files before any mask (line number, section
#   heading, record block). Also the record's steps (order, state, operation), the L lines with their block header, the
#   Save's outcome and the trx RUV k. Files: g0-lines.tsv, g0-diff.tsv, g0-steps.tsv, g0-journeys.txt, the raw files.
# w1: each case of `omicSummaryW1Cases` (tests/testthat/test-omic-summary-module.R) driven as the test drives it, with
#   every study_parameters.txt the real writer writes kept: whether it holds an L line (and which), per write, with the
#   steps whose `study_parameters` capture changed. File: g0-w1.tsv.
# Exit status 0 unless the setup itself fails; a journey or case that errors is written as such (an instrument result).

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, args[[2L]] %in% c("journeys", "w1"))
dir.create(args[[1L]], showWarnings = FALSE, recursive = TRUE)
out <- normalizePath(args[[1L]], mustWork = TRUE)
part <- args[[2L]]

suppressMessages(pkgload::load_all(".", export_all = TRUE, helpers = TRUE, quiet = TRUE))
library(testthat)
Sys.setenv(NOT_CRAN = "true")
env <- new.env(parent = asNamespace("MultiScholaR"))

#' Evaluate a test file's top-level assignments to a name (functions, constants), not its tests or loops.
k213Defs <- function(file, env) {
    for (e in parse(file.path("tests", "testthat", file), keep.source = FALSE)) {
        if (is.call(e) && identical(e[[1L]], as.name("<-")) && is.name(e[[2L]])) eval(e, env)
    }
}
k213Defs("test-kit207-summary-head-guard.R", env)
k213Defs("test-workflow-analysis-record.R", env)

evalq({
    L_PATTERN <- "^  ([^ ]+\\.)?(pairs|logcpm|grid|results) \\[input\\] = "
    tab <- function(x) gsub("[\t\r\n]", " ", x)

    #' The store's incarnation and its entity ids: the store's index, and every `entity_id` field under the head's
    #' states, products, exports and DA results, and under the DA manifest.
    k213StoreIds <- function(owner) {
        context <- owner$.context
        walk <- function(x) {
            if (!is.list(x)) return(character())
            found <- if (is.character(x$entity_id)) x$entity_id else character()
            c(found, unlist(lapply(unname(x), walk), use.names = FALSE))
        }
        snapshot <- context$snapshot
        manifest <- tryCatch(resolveWorkflowStateValue(context, snapshot$da_results), error = function(e) NULL)
        ids <- c(ls(envir = context$store$records, all.names = TRUE),
            walk(snapshot[intersect(c("states", "products", "exports", "da_results"), names(snapshot))]), walk(manifest))
        list(incarnation = context$store$identity$incarnation, ids = unique(ids[grepl("^[0-9a-f]{64}$", ids)]))
    }

    #' Per line: the heading of its section (the line before a dashes or equals rule) and the record block header.
    k213Context <- function(lines) {
        rule <- grepl("^(-{3,}|={3,})$", lines)
        heading <- block <- character(length(lines))
        current <- ""
        header <- ""
        for (i in seq_along(lines)) {
            if (i < length(lines) && rule[[i + 1L]]) { current <- lines[[i]]; header <- "" }
            if (grepl("^\\[[0-9]+\\. |^\\[[^]]+: [^]]+\\]$", lines[[i]])) header <- lines[[i]]
            heading[[i]] <- current
            block[[i]] <- header
        }
        list(heading = heading, block = block)
    }

    #' E1's masks: which lines are dropped and how many root occurrences each line holds.
    k213MaskKind <- function(lines, roots) {
        roots <- unique(c(normalizePath(roots, mustWork = FALSE), roots))
        roots <- roots[order(-nchar(roots))]
        kind <- ifelse(grepl("Timestamp:", lines, fixed = TRUE), "drop:Timestamp",
            ifelse(grepl("^  timestamp = ", lines), "drop:timestamp",
                ifelse(grepl("^(Repository|Branch|Commit|Package Version|Source):", lines), "drop:version", "")))
        count <- vapply(lines, function(line) {
            n <- 0L
            for (root in roots) {
                hits <- gregexpr(root, line, fixed = TRUE)[[1L]]
                if (hits[[1L]] > 0L) { n <- n + length(hits); line <- gsub(root, "<root>", line, fixed = TRUE) }
            }
            n
        }, 0L, USE.NAMES = FALSE)
        list(kind = kind, root_count = count)
    }

    #' One project's result: the file's lines, the store ids, the record's steps, the save's outcome, the roots.
    k213Result <- function(lines, owner, roots, rec, extra = list()) {
        record <- tryCatch(buildWorkflowAnalysisRecord(owner), error = function(e) e)
        steps <- if (inherits(record, "error")) data.frame(order = NA, state = conditionMessage(record), operation = "") else
            data.frame(order = vapply(record$steps, function(s) as.integer(s$order), 0L),
                state = vapply(record$steps, function(s) as.character(s$state), ""),
                operation = vapply(record$steps, function(s) as.character(s$operation), ""))
        analyses <- if (inherits(record, "error")) character() else vapply(names(record$analyses), function(name)
            sprintf("%s: %s source_state=%s", name, record$analyses[[name]]$operation,
                as.character(record$analyses[[name]]$source_state %||% NA)), "")
        c(list(lines = lines, ids = k213StoreIds(owner), roots = roots, steps = steps, analyses = analyses,
            done = rec$done, failed = paste(vapply(rec$failed, function(e) paste(class(e)[[1L]], conditionMessage(e)), ""),
                collapse = " | "), notes = paste(rec$notes, collapse = " | ")), extra)
    }

    k213ReadParams <- function(dirs, omic) {
        file <- file.path(dirs[[omic]]$source_dir, "study_parameters.txt")
        if (file.exists(file)) readLines(file, warn = FALSE) else "<no study_parameters.txt>"
    }

    #' Lipidomics or metabolomics: KIT-207 G1's journey (run, Apply, DA, Save) in a fresh root.
    k213FeatureJourney <- function(class) {
        root <- withr::local_tempdir()
        seen <- k194Notes()
        omic <- k194Omic(class)
        q <- k207Normalised(file.path(root, "fx"), class, da = TRUE)
        dirs <- k207Dirs(file.path(root, "project"), omic)
        rec <- k207Save(q$fx$workflow, omic, dirs)
        k213Result(k213ReadParams(dirs, omic), k194Owner(q$fx), root, rec, list(head = k194Head(q$fx)))
    }

    #' Transcriptomics: KIT-207 G5's journey with RUV automatic (manual k 2 when `manual`), the DA, the RUVIII export
    #' (`group`) when k >= 1, the Save.
    k213TrxJourney <- function(manual = FALSE) {
        seen <- new.env(parent = emptyenv())
        seen$notes <- character()
        project <- module_ci_trx_project(module_ci_trx_corpus())
        owner <- project$owner
        inputs <- if (manual) list(ruv_mode = "manual", ruv_k = 2L, ruv_percentage = 10) else list(ruv_mode = "automatic")
        k194TrxRun(project, omicNormReactiveState(), inputs, seen)
        head <- owner$current_state
        kit172Quiet(runTrxReferenceDaAnalysis(owner, getWorkflowStateReference(owner), list(da_q_val_thresh = 0.05,
            treat_lfc_cutoff = 0), trxDaContrastEntries(project$workflow$contrasts_tbl)))
        export <- tryCatch({
            kit172Quiet(runTrxReferenceRuviiiExport(owner, "group"))
            if (is.null(owner$.context$snapshot$exports$trx_ruviii_export)) "absent" else "accepted"
        }, error = function(e) paste("error:", conditionMessage(e)))
        # The export's k, as the W1 harness reads it (`omicSummaryW1Ruviii`).
        accepted <- owner$.context$snapshot$exports$trx_ruviii_export
        k <- if (is.null(accepted)) NA else tryCatch(hydrateSavedWorkflowEntity(owner$.context$store, accepted$reference)$value$k,
            error = function(e) paste("unreadable:", conditionMessage(e)))
        rec <- k207Save(project$workflow, "transcriptomics", project$project_dirs)
        base <- project$paths$base_dir
        k213Result(k213ReadParams(project$project_dirs, "transcriptomics"), owner, c(dirname(dirname(base)), base), rec,
            list(head = head, k = paste(k, collapse = ","), export = export, notes_run = paste(seen$notes, collapse = " | ")))
    }

    #' Proteomics: the KIT-107 F-16 fixture (MCI-009's object as the RUV state), the correlation filter at 0.5, the DA,
    #' its publication, the enrichment (test mode), the Save.
    k213ProtJourney <- function() {
        withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true")
        enrichment <- omicGuiProteomicsEnrichment
        testthat::local_mocked_bindings(omicGuiProteomicsEnrichment = function() recordEnrichmentConfig(enrichment()),
            .package = "MultiScholaR")
        root <- withr::local_tempdir()
        manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "kit213-prot", "DIA"), WorkflowState$new(audit_enabled = FALSE))
        withr::defer(closeDIAStateAdapter(manager))
        manager$.context$reference_states <- TRUE
        manager$setWorkflowType("DIA")
        design <- module_ci_prot_da_design()
        rownames(design) <- NULL
        quant <- module_ci_prot_da_table()
        rownames(quant) <- NULL
        setDIAStepProvenance(manager, "design_builder")
        manager$saveState("ruv_corrected", module_ci_prot_da_object(design = design, protein_quant_table = quant), list(), "fixture")
        paths <- list(source_dir = file.path(root, "source"), results_dir = file.path(root, "results"),
            da_output_dir = file.path(root, "results", "da_proteins"), pathway_dir = file.path(root, "results", "pathway_enrichment"),
            publication_graphs_dir = file.path(root, "publication_graphs"))
        for (p in paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
        manager$.context$experiment_paths <- paths
        contrasts <- module_ci_prot_da_contrasts("two_group")
        workflow <- list2env(list(state_manager = manager, contrasts_tbl = contrasts, tab_status = list(),
            qc_params = list(protein_qc = list()), config_list = list()))
        withDIAStateScope(manager, function() {
            run <- computeDIAReferenceCorrelation(manager, getDIAStateReference(manager), "group", 0.5)
            commit <- diaCorrelationCommit(manager, run$final, 0.5, FALSE)
            commitDIAReferenceState(manager, commit$state_name, commit$checkpoint, commit$config, commit$description)
            da_data <- new.env(parent = emptyenv())
            da_data$current_s4_object <- getDIAStateReference(manager)
            da_data$contrasts_available <- contrasts$contrasts
            manifest <- suppressMessages(suppressWarnings(runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group",
                da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), da_data, workflow, paths, contrasts = contrasts, publish = FALSE)))
            suppressMessages(suppressWarnings(produceDIADaPublication(manager, manifest, data.frame(), 0.05,
                paths[c("da_output_dir", "publication_graphs_dir")], save_assets = FALSE)))
            utils::capture.output(suppressMessages(suppressWarnings(runWorkflowEnrichment(workflow,
                omicGuiConfig("proteomics")$enrichment, manifest, list(selected_contrast = "B_vs_A", organism_taxid = "9606",
                    up_cutoff = 0, down_cutoff = 0, q_cutoff = 0.05, correction_method = "gSCS", enable_organism_filter = FALSE),
                paths))))
        })
        head <- manager$current_state
        dirs <- k207Dirs(file.path(root, "project"), "proteomics")
        rec <- k207Save(workflow, "proteomics", dirs)
        k213Result(k213ReadParams(dirs, "proteomics"), manager, root, rec, list(head = head))
    }

    k213Run <- function(label, journey) {
        result <- tryCatch(journey(), error = function(e) list(error = paste(class(e)[[1L]], conditionMessage(e))))
        result$label <- label
        result
    }

    #' The rows of g0-lines.tsv for one project: store-bound lines, masked lines, L lines.
    k213LineRows <- function(omic, project, result) {
        lines <- result$lines
        context <- k213Context(lines)
        masks <- k213MaskKind(lines, result$roots)
        tokens <- regmatches(lines, gregexpr("[0-9a-f]{64}", lines))
        store <- vapply(seq_along(lines), function(i) {
            what <- c(if (nzchar(result$ids$incarnation %||% "") && grepl(result$ids$incarnation, lines[[i]], fixed = TRUE))
                "incarnation", if (any(tokens[[i]] %in% result$ids$ids)) "entity_id")
            paste(what, collapse = "+")
        }, "")
        l_line <- grepl(L_PATTERN, lines)
        keep <- nzchar(store) | nzchar(masks$kind) | masks$root_count > 0L | l_line
        if (!any(keep)) return(NULL)
        i <- which(keep)
        data.frame(omic = omic, project = project, line = i, store_bound = store[i], mask = masks$kind[i],
            root_count = masks$root_count[i], l_line = l_line[i], heading = tab(context$heading[i]),
            block = tab(context$block[i]), text = tab(lines[i]))
    }

    #' The rows of g0-diff.tsv: every line differing between the two projects' files before any mask.
    k213DiffRows <- function(omic, a, b) {
        n <- max(length(a$lines), length(b$lines))
        la <- c(a$lines, rep("<no line>", n - length(a$lines)))
        lb <- c(b$lines, rep("<no line>", n - length(b$lines)))
        i <- which(la != lb)
        if (!length(i)) return(data.frame(omic = omic, line = NA, heading = "", block = "", project_1 = "<equal>",
            project_2 = sprintf("<equal: %d lines>", n)))
        context <- k213Context(la)
        data.frame(omic = omic, line = i, heading = tab(context$heading[i]), block = tab(context$block[i]),
            project_1 = tab(la[i]), project_2 = tab(lb[i]))
    }
}, env)

if (identical(part, "journeys")) evalq({
    journeys <- list(lipidomics = function() k213FeatureJourney("lipid"), metabolomics = function() k213FeatureJourney("metab"),
        transcriptomics = function() k213TrxJourney(), proteomics = k213ProtJourney)
    line_rows <- list()
    diff_rows <- list()
    step_rows <- list()
    notes <- character()
    omics <- names(journeys)
    while (length(omics)) {
        omic <- omics[[1L]]
        omics <- omics[-1L]
        pair <- lapply(1:2, function(p) k213Run(sprintf("%s-%d", omic, p), journeys[[omic]]))
        # E3's deviation (K213-R06): should RUV automatic give no RUVIII export (k = 0), manual k 2 as KIT-207 E6.
        if (identical(omic, "transcriptomics") && !identical(pair[[1L]]$export, "accepted")) {
            journeys$transcriptomics_manual <- function() k213TrxJourney(manual = TRUE)
            omics <- c(omics, "transcriptomics_manual")
        }
        for (p in 1:2) {
            r <- pair[[p]]
            if (!is.null(r$error)) {
                notes <- c(notes, sprintf("%s project %d ERROR: %s", omic, p, r$error))
                next
            }
            writeLines(r$lines, file.path(out, sprintf("g0-params-%s-%d.txt", omic, p)))
            line_rows[[length(line_rows) + 1L]] <- k213LineRows(omic, p, r)
            step_rows[[length(step_rows) + 1L]] <- cbind(omic = omic, project = p, r$steps)
            notes <- c(notes, sprintf("%s project %d: head=%s done=%s failed=%s ids=%d incarnation_lines=%d %s%s", omic, p,
                r$head %||% "", r$done, r$failed, length(r$ids$ids),
                sum(grepl(r$ids$incarnation %||% "<none>", r$lines, fixed = TRUE)),
                if (!is.null(r$k)) sprintf("k=%s export=%s ", r$k, r$export) else "",
                paste("analyses:", paste(r$analyses, collapse = "; "))), sprintf("  notes: %s", r$notes))
        }
        if (is.null(pair[[1L]]$error) && is.null(pair[[2L]]$error)) diff_rows[[length(diff_rows) + 1L]] <-
            k213DiffRows(omic, pair[[1L]], pair[[2L]])
    }
    utils::write.table(do.call(rbind, line_rows), file.path(out, "g0-lines.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(do.call(rbind, diff_rows), file.path(out, "g0-diff.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(do.call(rbind, step_rows), file.path(out, "g0-steps.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    writeLines(notes, file.path(out, "g0-journeys.txt"))
    cat(notes, sep = "\n")
}, env)

if (identical(part, "w1")) evalq({
    k213Defs("test-omic-summary-module.R", env)
    real <- get("createWorkflowArgsFromConfig", envir = asNamespace("MultiScholaR"))
    w1Case <- function(case) {
        saved <- list()
        testthat::local_mocked_bindings(createWorkflowArgsFromConfig = function(...) {
            path <- real(...)
            saved[[length(saved) + 1L]] <<- readLines(path, warn = FALSE)
            path
        }, .package = "MultiScholaR")
        captures <- omicSummaryW1Case(case)
        changed <- names(captures)[vapply(captures, function(step) !isTRUE(step$study_parameters %in% c("<absent>",
            "<unchanged>")) && !is.null(step$study_parameters), logical(1))]
        list(saved = saved, changed = changed, values = vapply(changed, function(s) captures[[s]]$study_parameters, ""))
    }
    rows <- list()
    for (case in omicSummaryW1Cases) {
        result <- tryCatch(w1Case(case), error = function(e) list(error = paste(class(e)[[1L]], conditionMessage(e))))
        if (!is.null(result$error)) {
            rows[[length(rows) + 1L]] <- data.frame(case = case, write = NA, lines = NA, l_lines = NA,
                steps_changed = "", captures = "", l_text = tab(paste("ERROR:", result$error)))
            next
        }
        steps <- paste(result$changed, collapse = ",")
        values <- paste(result$values, collapse = " | ")
        if (!length(result$saved)) {
            rows[[length(rows) + 1L]] <- data.frame(case = case, write = 0L, lines = 0L, l_lines = 0L,
                steps_changed = steps, captures = values, l_text = "<no study_parameters.txt written>")
        }
        for (w in seq_along(result$saved)) {
            lines <- result$saved[[w]]
            l <- grep(L_PATTERN, lines, value = TRUE)
            rows[[length(rows) + 1L]] <- data.frame(case = case, write = w, lines = length(lines), l_lines = length(l),
                steps_changed = steps, captures = values,
                l_text = tab(paste(gsub("[0-9a-f]{64}", "<HEX>", l), collapse = " || ")))
        }
    }
    table <- do.call(rbind, rows)
    utils::write.table(table, file.path(out, "g0-w1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    print(table[, c("case", "write", "lines", "l_lines", "steps_changed")])
}, env)
