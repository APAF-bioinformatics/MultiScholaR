# KIT-107 W1: the enrichment tab through its public entry points, written before the shared enrichment frame
# (R/mod_omic_enrich.R, R/mod_omic_enrich_panels.R; W3) and the core moves (W2) exist, and frozen on the start commit
# (8e7294b6; R/ unchanged since b9820d9e; docs/qa/workflow-execution/kit/kit-107/w1/test-first.txt).
#   (a) Identity: the HTML digest of `mod_prot_enrich_ui("enrichment_analysis")` under MULTISCHOLAR_TEST_MODE=true
#       (shiny's tabset ids masked) and its namespaced input and output ids, labels, help texts, choice lists and
#       defaults (`omicsUiInputDefaults`).
#   (b) The server in process: `mod_prot_da_server` and `mod_prot_enrich_server` share one bus in one `testServer`,
#       on KIT-106 W1's proteomics fixture shape (test-omic-da-module.R's `omicDaW1ProtSetup`): an on-disk DIA owner
#       holding MCI-009's object as `correlation_filtered` and its contrasts table. The DA runs first and publishes the
#       manifest and the DA tables; then the enrichment steps run. Cases:
#       gprofiler (taxon 9606) and clusterprofiler (taxon 999999) on the deterministic double: the tab, the contrast,
#         Run, each direction filter, the download;
#       kernel_gprofiler and kernel_clusterprofiler: the double replaced, by name, with a recorder that calls the live
#         `processEnrichments` under `install_enrichment_doubles()`, the fake `gost` and `enricher` recording their
#         query and background (the live kernel path, offline);
#       organism_filter: the filter on with a mapping found (today: the `clean_acc` failure, D11, P-07), the filter on
#         with no mapping (the 8 s warning, an unfiltered run), the mixed-species observer (the checkbox and its notice);
#       refusals: the tab opened before the DA, Run before the DA, Run with no contrast, taxon `abc`, `q_cutoff` 1.5, a
#         DA table removed before the Run (each failure that reaches the run's manifest branch ends its session: a
#         second Run on the same tab is F107-01's, which W3 fixes);
#       restore: `restoreProtEnrichResultsFromState` on a NULL holder over an owner that holds `enrichment_completed`
#         (the KIT-024 measure, at function level; the module-level Open is F107-02's, which W3 fixes);
#       prefill: the bus `taxon_id` sets `organism_taxid`.
#   Captured after each step: the registered enrichment outputs and each value as the entry suite reads it
#   (`guiEntryOutputValue`), the notifications (type, duration, id, message; each removal), the progress, the selector
#   updates, `tab_status`, the history, each state's digest (store references and paths masked), the hydrated product
#   (schema, parameters, results, contrasts, service) and the external resources, the kernel's arguments (values, and
#   digests of `da_results@da_data` and `go_annotations`; the pathway directory masked), the fake services' query and
#   background digests, the bus values (times masked), the display tables, the project's files (table digests) and
#   store entry counts, and the server's returned reactive values (names; value digests, store references masked).
#   Step provenance and the analysis record are excluded (D6).
# MULTISCHOLAR_KIT107_W1_DUMP=<file> writes the captures as R code (how the frozen values were made; two dumps in
# separate processes were byte-identical before they were pasted).

omicEnrichW1Sha <- function(x) digest::digest(serialize(x, NULL, version = 3L), algo = "sha256", serialize = FALSE)

#' Where MULTISCHOLAR_KIT107_W1_DUMP points, or NULL.
omicEnrichW1DumpPath <- function() {
    path <- Sys.getenv("MULTISCHOLAR_KIT107_W1_DUMP")
    if (nzchar(path)) path else NULL
}

#' R code for a capture (lists, named vectors and scalars): a list's elements one per line, a long vector's packed
#' into lines, indented (test-omic-da-module.R's `omicDaW1Code`).
omicEnrichW1Code <- function(x, indent = 0L) {
    pad <- strrep(" ", indent + 4L)
    label <- function(names, i) {
        if (is.null(names) || !nzchar(names[[i]])) return("")
        paste0(if (make.names(names[[i]]) == names[[i]]) names[[i]] else paste0("`", names[[i]], "`"), " = ")
    }
    if (is.list(x)) {
        if (length(x) == 0L) return("list()")
        items <- vapply(seq_along(x), function(i) paste0(pad, label(names(x), i), omicEnrichW1Code(x[[i]], indent + 4L)), "")
        return(paste0("list(\n", paste(items, collapse = ",\n"), ")"))
    }
    if (is.null(x)) return("NULL")
    scalar <- function(value) deparse(value, control = c("keepNA", "keepInteger"), width.cutoff = 500L)
    if (length(x) == 0L) return(if (is.null(names(x))) scalar(x) else sprintf("stats::setNames(%s, character())", scalar(unname(x))))
    if (length(x) == 1L && is.null(names(x))) return(scalar(x))
    items <- paste0(vapply(seq_along(x), function(i) paste0(label(names(x), i), scalar(unname(x[[i]]))), ""), ",")
    items[[length(items)]] <- sub(",$", ")", items[[length(items)]])
    if (sum(nchar(items) + 1L) + indent < 90L) return(paste0("c(", paste(items, collapse = " ")))
    lines <- character()
    for (item in items) {
        last <- length(lines)
        if (last && nchar(lines[[last]]) + nchar(item) < 110L) {
            lines[[last]] <- paste(lines[[last]], item)
        } else {
            lines <- c(lines, paste0(pad, item))
        }
    }
    paste0("c(\n", paste(lines, collapse = "\n"))
}

#' Append `name <- <value>` to the dump file.
omicEnrichW1Dump <- function(name, value) {
    cat(sprintf("%s <- %s\n\n", name, omicEnrichW1Code(value)), file = omicEnrichW1DumpPath(), append = TRUE)
}

# ------------------------------------------------------------------------------------------------
# (a) Identity: the enrichment UI, frozen at 8e7294b6.
# ------------------------------------------------------------------------------------------------

#' The enrichment UI as the shell mounts it (mod_omic_workflow.R:222; module id `enrichment_analysis`).
omicEnrichW1Ui <- function() {
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true")
    mod_prot_enrich_ui("enrichment_analysis")
}

#' The rendered HTML with shiny's random tabset ids masked.
omicEnrichW1Html <- function(ui) {
    html <- as.character(htmltools::renderTags(ui)$html)
    html <- gsub("data-tabsetid=\"[0-9]+\"", "data-tabsetid=\"ID\"", html)
    gsub("tab-[0-9]+-([0-9]+)", "tab-ID-\\1", html)
}

#' The UI's HTML digest and its parts: inputs (every namespaced `id` that is not an output, a label or a tab pane),
#' outputs, `data-testid` values, tab titles in order, each label's text by its input, the help texts in order, each
#' select's choices (`value=text`, the selected one marked) and the inputs' defaults.
omicEnrichW1UiParts <- function(ui) {
    html <- omicEnrichW1Html(ui)
    tags <- regmatches(html, gregexpr("<[a-zA-Z][^>]*>", html))[[1L]]
    attr_of <- function(name) vapply(tags, function(tag) {
        m <- regmatches(tag, regexec(sprintf("[[:space:]]%s=\"([^\"]*)\"", name), tag))[[1L]]
        if (length(m)) m[[2L]] else NA_character_
    }, "", USE.NAMES = FALSE)
    ids <- attr_of("id")
    classes <- attr_of("class")
    is_output <- !is.na(ids) & grepl("(^| )(shiny-[a-z]+-output|html-widget-output)( |$)", classes)
    inputs <- ids[!is.na(ids) & !is_output & !grepl("-label$", ids) & !startsWith(ids, "tab-ID-")]
    labels <- regmatches(html, gregexpr("<label[^>]*for=\"[^\"]+\"[^>]*>[^<]*</label>", html))[[1L]]
    label_text <- stats::setNames(sub("^<label[^>]*>([^<]*)</label>$", "\\1", labels),
        sub("^.*for=\"([^\"]+)\".*$", "\\1", labels))
    boxes <- regmatches(html, gregexpr("<input id=\"[^\"]+\" type=\"checkbox\"[^>]*/>[[:space:]]*<span>[^<]*</span>", html))[[1L]]
    box_text <- stats::setNames(sub("^.*<span>([^<]*)</span>$", "\\1", boxes), sub("^<input id=\"([^\"]+)\".*$", "\\1", boxes))
    help <- regmatches(html, gregexpr("<span class=\"help-block\">[^<]*</span>", html))[[1L]]
    selects <- regmatches(html, gregexpr("(?s)<select id=\"[^\"]+\"[^>]*>.*?</select>", html, perl = TRUE))[[1L]]
    choices <- lapply(stats::setNames(selects, sub("^<select id=\"([^\"]+)\".*$", "\\1", selects)), function(select) {
        options <- regmatches(select, gregexpr("<option value=\"[^\"]*\"( selected)?>[^<]*</option>", select))[[1L]]
        sub("^<option value=\"([^\"]*)\"( selected)?>([^<]*)</option>$", "\\1=\\3\\2", options)
    })
    defaults <- omicsUiInputDefaults(ui)
    defaults <- vapply(defaults[order(names(defaults))], function(value) paste(format(value), collapse = ","), "")
    list(html = digest::digest(html, algo = "sha256", serialize = FALSE), inputs = sort(unique(inputs)),
        outputs = sort(unique(ids[is_output])), testids = sort(unique(stats::na.omit(attr_of("data-testid")))),
        tabs = attr_of("data-value")[grepl("data-toggle=\"tab\"", tags)],
        labels = c(label_text, box_text)[order(names(c(label_text, box_text)))],
        help = sub("^<span class=\"help-block\">([^<]*)</span>$", "\\1", help),
        choices = if (length(choices)) choices[order(names(choices))] else list(),
        defaults = defaults)
}

# OMIC_ENRICH_FROZEN_UI_START
omicEnrichW1FrozenUi <- list(
    html = "3c213c43b0003c7ca879afc155c682b8b451935e7244e4f542a075fe481822ba",
    inputs = c(
        "enrichment_analysis-clusterprofiler_direction_filter", "enrichment_analysis-correction_method",
        "enrichment_analysis-down_cutoff", "enrichment_analysis-download_enrichment_results",
        "enrichment_analysis-enable_organism_filter", "enrichment_analysis-enrichment_method_tabs",
        "enrichment_analysis-gprofiler_direction_filter", "enrichment_analysis-organism_taxid",
        "enrichment_analysis-q_cutoff", "enrichment_analysis-run_enrichment_analysis",
        "enrichment_analysis-selected_contrast", "enrichment_analysis-stringdb_filter_significant",
        "enrichment_analysis-stringdb_max_results", "enrichment_analysis-stringdb_ranking_method",
        "enrichment_analysis-up_cutoff"),
    outputs = c(
        "enrichment_analysis-analysis_method_display", "enrichment_analysis-clusterprofiler_plot",
        "enrichment_analysis-clusterprofiler_results_table",
        "enrichment_analysis-clusterprofiler_summary_stats", "enrichment_analysis-contrasts_display",
        "enrichment_analysis-enrichment_status", "enrichment_analysis-gprofiler_plot",
        "enrichment_analysis-gprofiler_results_table", "enrichment_analysis-gprofiler_summary_stats",
        "enrichment_analysis-stringdb_plot", "enrichment_analysis-stringdb_results_table",
        "enrichment_analysis-stringdb_summary_stats"),
    testids = character(0),
    tabs = c("gprofiler2", "clusterprofiler", "stringdb"),
    labels = c(
        `enrichment_analysis-clusterprofiler_direction_filter` = "Direction:",
        `enrichment_analysis-correction_method` = "Correction Method:",
        `enrichment_analysis-down_cutoff` = "Down Log2FC Cutoff:",
        `enrichment_analysis-enable_organism_filter` = "Filter to target organism only",
        `enrichment_analysis-gprofiler_direction_filter` = "Direction:",
        `enrichment_analysis-organism_taxid` = "Organism Taxon ID:",
        `enrichment_analysis-q_cutoff` = "Q-value Cutoff:",
        `enrichment_analysis-selected_contrast` = "Select Contrast:",
        `enrichment_analysis-stringdb_filter_significant` = "Filter Significant Only",
        `enrichment_analysis-stringdb_max_results` = "Max Results:",
        `enrichment_analysis-stringdb_ranking_method` = "Ranking Method:",
        `enrichment_analysis-up_cutoff` = "Up Log2FC Cutoff:"),
    help = c(
        "Choose which DE contrast to analyze", "Minimum log2 fold change for up-regulated genes",
        "Minimum log2 fold change for down-regulated genes", "FDR threshold for enrichment significance",
        "NCBI taxonomy ID for organism",
        "Enable if using mixed-species FASTA (e.g., with contaminants). Filters DE results to target organism before enrichment analysis.",
        "Automatically determined based on organism support", "Contrasts from DE analysis",
        "g:SCS: More conservative, reduces false positives but may miss terms.",
        "FDR: Detects more terms but higher false positive risk."),
    choices = list(
        `enrichment_analysis-clusterprofiler_direction_filter` = c("all=All selected", "up=Up-regulated only", "down=Down-regulated only"),
        `enrichment_analysis-correction_method` = c("gSCS=g:SCS (Conservative) selected", "fdr=FDR (Less Conservative)"),
        `enrichment_analysis-gprofiler_direction_filter` = c("all=All selected", "up=Up-regulated only", "down=Down-regulated only"),
        `enrichment_analysis-selected_contrast` = character(0),
        `enrichment_analysis-stringdb_ranking_method` = c(
            "combined_score=Combined Score selected", "fdr_qvalue=FDR Q-value", "log2fc=Log2 Fold Change")),
    defaults = c(
        `enrichment_analysis-clusterprofiler_direction_filter` = "all",
        `enrichment_analysis-correction_method` = "gSCS", `enrichment_analysis-down_cutoff` = "0",
        `enrichment_analysis-enable_organism_filter` = "FALSE",
        `enrichment_analysis-gprofiler_direction_filter` = "all",
        `enrichment_analysis-organism_taxid` = "9606", `enrichment_analysis-q_cutoff` = "0.05",
        `enrichment_analysis-stringdb_filter_significant` = "FALSE",
        `enrichment_analysis-stringdb_max_results` = "50",
        `enrichment_analysis-stringdb_ranking_method` = "combined_score",
        `enrichment_analysis-up_cutoff` = "0"))
# OMIC_ENRICH_FROZEN_UI_END

test_that("(a) the enrichment UI keeps its HTML, ids, test ids, tab titles, labels, help texts, choices and defaults", {
    actual <- omicEnrichW1UiParts(omicEnrichW1Ui())
    if (!is.null(omicEnrichW1DumpPath())) omicEnrichW1Dump("omicEnrichW1FrozenUi", actual)
    expect_named(actual, names(omicEnrichW1FrozenUi))
    for (part in names(omicEnrichW1FrozenUi)) {
        expect_identical(actual[[part]], omicEnrichW1FrozenUi[[part]], info = part)
    }
})

# ------------------------------------------------------------------------------------------------
# (b) The server, driven in testServer.
# ------------------------------------------------------------------------------------------------

#' The browser, notification and progress doubles, scoped to `env`. A notification is kept as
#' "[type duration id] message" (shiny's default duration 5 when the caller gives none), a removal as "remove <id>", a
#' modal as "modal"; a selector update as "<id>: <choices> -> <selected or value>", and an update with a value is
#' queued for the browser double, which sends it back as the input. The progress is kept as "start <message> <value>"
#' and "<amount> <detail>".
omicEnrichW1InstallDoubles <- function(log, env = parent.frame()) {
    note <- function(text) log$notes <- c(log$notes, text)
    notify <- function(ui, action = NULL, duration = 5, closeButton = TRUE, id = NULL,
        type = c("default", "message", "warning", "error"), session = NULL) {
        note(sprintf("[%s %s %s] %s", type[[1L]], if (is.null(duration)) "NULL" else format(duration), id %||% "-",
            paste(as.character(ui), collapse = " ")))
        invisible(id %||% "notification")
    }
    queue <- function(session, input_id, value = NULL, choices = NULL) {
        if (is.null(session)) return(invisible(NULL))
        id <- session$ns(input_id)
        log$updates <- c(log$updates, sprintf("%s: %s -> %s", id,
            if (is.null(choices)) "<same>" else paste(as.character(unlist(choices, use.names = FALSE)), collapse = "|"),
            if (is.null(value)) "<none>" else paste(format(value), collapse = "|")))
        if (!is.null(value)) log$pending[[length(log$pending) + 1L]] <- list(id = id, value = value)
        invisible(NULL)
    }
    domain <- function() shiny::getDefaultReactiveDomain()
    bindings <- list(
        updateSelectInput = function(session = domain(), inputId, label = NULL, choices = NULL, selected = NULL) {
            queue(session, inputId, selected, choices)
        },
        updateSelectizeInput = function(session = domain(), inputId, label = NULL, choices = NULL, selected = NULL,
            options = list(), server = FALSE) queue(session, inputId, selected, choices),
        updateTextAreaInput = function(session = domain(), inputId, label = NULL, value = NULL, placeholder = NULL) {
            queue(session, inputId, value)
        },
        updateTextInput = function(session = domain(), inputId, label = NULL, value = NULL, placeholder = NULL) {
            queue(session, inputId, value)
        },
        updateNumericInput = function(session = domain(), inputId, label = NULL, value = NULL, min = NULL, max = NULL,
            step = NULL) queue(session, inputId, value),
        updateCheckboxInput = function(session = domain(), inputId, label = NULL, value = NULL) queue(session, inputId, value),
        updateTabsetPanel = function(session = domain(), inputId, selected = NULL) queue(session, inputId, selected),
        showNotification = notify,
        removeNotification = function(id, session = domain()) {
            note(paste("remove", id))
            invisible(NULL)
        },
        showModal = function(ui, session = domain()) {
            note("modal")
            invisible(NULL)
        },
        removeModal = function(session = domain()) invisible(NULL),
        withProgress = function(expr, min = 0, max = 1, value = min + (max - min) * 0.1, message = NULL, detail = NULL,
            style = NULL, session = domain(), env = parent.frame(), quoted = FALSE) {
            log$progress <- c(log$progress, paste("start", message %||% "-", format(value)))
            if (!quoted) expr <- substitute(expr)
            eval(expr, env)
        },
        incProgress = function(amount = 0.1, message = NULL, detail = NULL, session = domain()) {
            log$progress <- c(log$progress, paste(format(amount), detail %||% message %||% "-"))
            invisible(NULL)
        })
    do.call(testthat::local_mocked_bindings, c(bindings, list(.package = "shiny", .env = env)))
    imported <- intersect(names(bindings), unlist(getNamespaceImports("MultiScholaR")$shiny, use.names = FALSE))
    if (length(imported)) do.call(testthat::local_mocked_bindings, c(bindings[imported], list(.package = "MultiScholaR", .env = env)))
    invisible(log)
}

#' A short content digest (the entry suite's canonical form).
omicEnrichW1Short <- function(x) substr(guiEntryDigest(x), 1L, 16L)

#' A value with its process-local parts masked, for a digest: each saved store reference (its entity id is the
#' store's) as `<ref>`, each string with the entry suite's masks (the temporary roots, times, hexadecimal ids), each
#' time as `<time>` (the canonical form).
omicEnrichW1Clean <- function(x, roots, depth = 0L) {
    if (depth > 14L) return("<deep>")
    if (isSavedWorkflowReference(x)) return("<ref>")
    if (is.list(x) && !is.data.frame(x) && !isS4(x) && !inherits(x, "htmlwidget") && !is.object(x)) {
        return(lapply(x, omicEnrichW1Clean, roots = roots, depth = depth + 1L))
    }
    if (isS4(x) && !is.null(methods::slotNames(x)) && length(methods::slotNames(x))) {
        slots <- methods::slotNames(x)
        return(c(list(`<class>` = class(x)[[1L]]), stats::setNames(lapply(slots, function(s)
            omicEnrichW1Clean(methods::slot(x, s), roots, depth + 1L)), slots)))
    }
    if (is.character(x) && is.null(dim(x)) && !is.object(x)) {
        return(stats::setNames(vapply(x, function(value) if (is.na(value)) NA_character_ else guiEntryMask(value, roots), "",
            USE.NAMES = FALSE), names(x)))
    }
    x
}

#' The digest of a value with its process-local parts masked.
omicEnrichW1Digest <- function(x, roots) omicEnrichW1Short(omicEnrichW1Clean(x, roots))

#' A list of scalars as text (the bus's parameter lists): a time as `<time>`, a string with the masks.
omicEnrichW1Scalars <- function(x, roots) {
    if (is.null(x)) return("<NULL>")
    vapply(x, function(value) {
        if (inherits(value, c("POSIXt", "Date"))) return("<time>")
        if (is.null(value)) return("<NULL>")
        if (is.atomic(value) && length(value) <= 5L) return(guiEntryMask(paste(format(value), collapse = "|"), roots))
        paste("<digest>", omicEnrichW1Digest(value, roots))
    }, "")
}

#' The kernel's arguments as text: the scalars, the pathway directory masked, and digests of the DA tables, the
#' contrasts and design of `da_results`, and of `go_annotations`.
omicEnrichW1KernelCall <- function(name, args, roots) {
    scalars <- c("taxon_id", "up_cutoff", "down_cutoff", "q_cutoff", "exclude_iea", "protein_id_column", "contrast_names",
        "correction_method")
    shown <- vapply(scalars, function(arg) paste(format(args[[arg]]), collapse = "|"), "")
    da <- args$da_results
    c(kernel = name, shown, pathway_dir = guiEntryMask(args$pathway_dir, roots),
        da_data = omicEnrichW1Short(da@da_data), da_names = paste(names(da@da_data), collapse = "|"),
        da_rows = paste(vapply(da@da_data, function(table) if (is.null(table)) "<NULL>" else sprintf("%dx%d", nrow(table),
            ncol(table)), ""), collapse = "|"),
        da_contrasts = omicEnrichW1Short(da@contrasts), da_design = omicEnrichW1Short(da@design_matrix),
        go_annotations = if (is.null(args$go_annotations)) "<NULL>" else omicEnrichW1Short(args$go_annotations))
}

#' Recorders on the kernel the module selects (the deterministic double, looked up by name in test mode), scoped to
#' `env`. With `live`, the double is replaced by a recorder that calls the live `processEnrichments` under
#' `install_enrichment_doubles()`, with the fake `gost` and `enricher` recording their query and background.
omicEnrichW1InstallKernel <- function(acc, roots, live = FALSE, env = parent.frame()) {
    ns <- asNamespace("MultiScholaR")
    double <- get("runProtEnrichDeterministicProcessEnrichments", envir = ns)
    kernel <- get("processEnrichments", envir = ns)
    recorder <- function(...) {
        args <- list(...)
        call <- omicEnrichW1KernelCall(if (live) "processEnrichments" else "runProtEnrichDeterministicProcessEnrichments",
            args, roots)
        acc$kernel <- c(acc$kernel, paste0(names(call), "=", call, collapse = " ; "))
        if (live) kernel(...) else double(...)
    }
    testthat::local_mocked_bindings(runProtEnrichDeterministicProcessEnrichments = recorder, .package = "MultiScholaR",
        .env = env)
    if (!live) return(invisible(acc))
    install_enrichment_doubles(.local_envir = env)
    ids <- function(x) if (is.null(x)) "<NULL>" else paste(length(x), omicEnrichW1Short(as.character(x)))
    gost <- function(query = character(), organism = "hsapiens", ordered_query = FALSE, multi_query = FALSE,
        significant = TRUE, exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = FALSE,
        user_threshold = 0.05, correction_method = "g_SCS", domain_scope = "annotated", custom_bg = NULL,
        numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE) {
        acc$services <- c(acc$services, paste(c("gost", paste0("organism=", organism),
            paste0("sources=", paste(sources, collapse = "|")), paste0("user_threshold=", format(user_threshold)),
            paste0("correction_method=", correction_method), paste0("exclude_iea=", exclude_iea),
            paste0("domain_scope=", domain_scope), paste0("significant=", significant), paste0("query=", ids(query)),
            paste0("custom_bg=", ids(custom_bg))), collapse = " ; "))
        fake_gost(query = query, organism = organism, ordered_query = ordered_query, multi_query = multi_query,
            significant = significant, exclude_iea = exclude_iea, measure_underrepresentation = measure_underrepresentation,
            evcodes = evcodes, user_threshold = user_threshold, correction_method = correction_method,
            domain_scope = domain_scope, custom_bg = custom_bg, numeric_ns = numeric_ns, sources = sources,
            as_short_link = as_short_link, highlight = highlight)
    }
    enricher <- function(gene = character(), universe = NULL, TERM2GENE = NULL, TERM2NAME = NULL, pvalueCutoff = 0.05,
        pAdjustMethod = "BH", ...) {
        acc$services <- c(acc$services, paste(c("enricher", paste0("gene=", ids(gene)), paste0("universe=", ids(universe)),
            paste0("TERM2GENE=", if (is.null(TERM2GENE)) "<NULL>" else omicEnrichW1Short(TERM2GENE)),
            paste0("pvalueCutoff=", format(pvalueCutoff)), paste0("pAdjustMethod=", pAdjustMethod)), collapse = " ; "))
        fake_enricher(gene = gene, universe = universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME,
            pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, ...)
    }
    if (requireNamespace("gprofiler2", quietly = TRUE)) {
        testthat::local_mocked_bindings(gost = gost, .package = "gprofiler2", .env = env)
    }
    if (requireNamespace("clusterProfiler", quietly = TRUE)) {
        .install_ns_binding(env = asNamespace("clusterProfiler"), name = "enricher", value = enricher, .local_envir = env)
    }
    invisible(acc)
}

#' The display tables the tab holds: each contrast's saved g:Profiler and clusterProfiler tables, read back from the
#' store (rows and a digest of every row), and the in-memory tables (a digest).
omicEnrichW1Display <- function(all_results, manager) {
    if (is.null(all_results)) return("<NULL>")
    if (!length(all_results)) return("<empty>")
    out <- character()
    for (contrast in names(all_results)) {
        entry <- all_results[[contrast]]
        for (kind in c("gprofiler", "clusterprofiler", "stringdb")) {
            table <- entry[[paste0(kind, "_results")]]
            out[[paste0(contrast, "/", kind)]] <- if (is.null(table)) "<NULL>" else
                sprintf("%dx%d %s", nrow(table), ncol(table), omicEnrichW1Short(table))
            reference <- entry[[paste0(kind, "_results_table")]]
            if (!is.null(reference)) {
                out[[paste0(contrast, "/", kind, "_table")]] <- if (isWorkflowTableReference(reference)) {
                    page <- readDIADisplayTablePage(reference, limit = max(1L, as.integer(reference$metadata$shape[[1L]])),
                        manager = manager)
                    sprintf("<saved> %d %s", page$rows, omicEnrichW1Short(page$data))
                } else paste0("<", class(reference)[[1L]], ">")
            }
        }
    }
    out
}

#' The product the run saved (its reference on the tab), hydrated: schema, service, parameters, results, contrasts and
#' a digest of its summaries. `results` holds the DA tables' store entity ids, which are the store's: each is shown as
#' the manifest table it names (`<manifest da_proteins_long>`), or `<other>`.
omicEnrichW1Product <- function(reference, manager, manifest) {
    if (is.null(reference)) return("<none>")
    product <- hydrateSavedWorkflowEntity(manager$.context$store, reference)$value
    parameters <- vapply(product$parameters, function(value) paste(format(value), collapse = "|"), "")
    tables <- vapply(names(product$results), function(name) {
        table <- manifest$results[[name]]$tables$da_proteins_long$entity_id
        if (!is.null(table) && identical(unname(product$results[[name]]), table)) "<manifest da_proteins_long>" else "<other>"
    }, "")
    c(schema = product$schema, service = paste(product$service, collapse = "|"),
        stats::setNames(parameters, paste0("parameters$", names(parameters))),
        stats::setNames(tables, paste0("results$", names(product$results))),
        contrasts = paste(product$contrasts, collapse = "|"), summaries = omicEnrichW1Short(product$summaries))
}

#' The external resources the head records (the query date masked).
omicEnrichW1Resources <- function(manager) {
    resources <- manager$.context$snapshot$external_resources
    if (!length(resources)) return(character())
    vapply(resources, function(resource) paste(resource$service, resource$release, resource$release_kind,
        resource$query_summary %||% "-", sep = " ; "), "")
}

#' Each state's digest, with store references and paths masked.
omicEnrichW1States <- function(manager, roots) {
    history <- guiEntryHistory(manager)
    out <- character()
    for (name in unique(history)) {
        out[[name]] <- tryCatch(omicEnrichW1Digest(withWorkflowStateScope(manager, function() manager$getState(name)), roots),
            error = function(e) paste("<error>", conditionMessage(e)))
    }
    out
}

#' The bus values the enrichment writes: the results payload (a digest), the UI parameters and the organism filter
#' (their scalars, times masked), the annotation table (a digest).
omicEnrichW1Bus <- function(bus, roots) {
    values <- shiny::isolate(list(results = bus$enrichment_analysis_results, ui = bus$enrichment_ui_params,
        filter = bus$enrichment_organism_filter, uniprot = bus$uniprot_dat_cln))
    c(enrichment_analysis_results = if (is.null(values$results)) "<NULL>" else
            paste(paste(names(values$results), collapse = ","), omicEnrichW1Digest(values$results, roots)),
        stats::setNames(omicEnrichW1Scalars(values$ui, roots), paste0("enrichment_ui_params$", names(values$ui) %||% "")),
        stats::setNames(omicEnrichW1Scalars(values$filter, roots),
            paste0("enrichment_organism_filter$", names(values$filter) %||% "")),
        uniprot_dat_cln = if (is.null(values$uniprot)) "<NULL>" else
            sprintf("%dx%d %s", nrow(values$uniprot), ncol(values$uniprot), omicEnrichW1Short(values$uniprot)))
}

#' The server's returned reactive values: each field's digest, store references and paths masked.
omicEnrichW1Returned <- function(values, roots) {
    if (!shiny::is.reactivevalues(values)) return(paste0("<", class(values)[[1L]], ">"))
    fields <- shiny::isolate(shiny::reactiveValuesToList(values))
    fields <- fields[order(names(fields))]
    vapply(fields, function(value) {
        if (is.null(value)) return("<NULL>")
        if (isSavedWorkflowReference(value)) return("<ref>")
        if (is.atomic(value) && length(value) <= 3L && !is.object(value)) return(paste(format(value), collapse = "|"))
        paste0("<", class(value)[[1L]], "> ", omicEnrichW1Digest(value, roots))
    }, "")
}

#' A message as captured: the entry suite's masks, then without terminal styling or non-ASCII characters.
omicEnrichW1Text <- function(x, roots) {
    x <- gsub("\033\\[[0-9;]*m", "", guiEntryMask(x, roots))
    iconv(x, "UTF-8", "ASCII", sub = "")
}

#' Drive one session: the DA and enrichment servers in testServer on the fixture, the two UIs' defaults set (as the
#' browser sends them), the steps run (named functions of the context) and captured after each.
#' @param setup The session's fixture (`omicEnrichW1Setup`).
#' @param live Whether the kernel is the live `processEnrichments` on the fake services.
#' @return The captures, by step.
omicEnrichW1Drive <- function(setup, steps, live = FALSE) {
    log <- new.env(parent = emptyenv())
    log$notes <- character()
    log$updates <- character()
    log$progress <- character()
    log$pending <- list()
    omicEnrichW1InstallDoubles(log)
    acc <- new.env(parent = emptyenv())
    acc$kernel <- character()
    acc$services <- character()
    acc$captures <- list()
    acc$previous <- list()
    acc$extra <- list()
    acc$errors <- character()
    omicEnrichW1InstallKernel(acc, setup$roots, live = live)
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true", http_proxy = "http://127.0.0.1:9", https_proxy = "http://127.0.0.1:9",
        HTTP_PROXY = "http://127.0.0.1:9", HTTPS_PROXY = "http://127.0.0.1:9", no_proxy = "", NO_PROXY = "")
    withr::local_seed(20261109L)
    enrich_prefix <- "enrichment_analysis-"
    server <- local({
        mount <- setup$server(acc)
        function(input, output, session) {
            session$onUnhandledError(function(error) acc$errors <- c(acc$errors, paste("<unhandled>", conditionMessage(error))))
            mount(input, output, session)
        }
    })
    text <- function(x) omicEnrichW1Text(x, setup$roots)
    suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(server, {
        apply_updates <- function() {
            for (i in seq_len(20L)) {
                if (!length(log$pending)) return(invisible())
                pending <- log$pending
                log$pending <- list()
                values <- stats::setNames(lapply(pending, `[[`, "value"), vapply(pending, `[[`, "", "id"))
                do.call(session$setInputs, values[!duplicated(names(values), fromLast = TRUE)])
            }
        }
        put <- function(values) {
            tryCatch({
                do.call(session$setInputs, values)
                apply_updates()
            }, error = function(e) acc$errors <- c(acc$errors, conditionMessage(e)))
            invisible(NULL)
        }
        set <- function(...) {
            values <- list(...)
            names(values) <- paste0(enrich_prefix, names(values))
            put(values)
        }
        set_da <- function(...) {
            values <- list(...)
            names(values) <- paste0("differential_expression-", names(values))
            put(values)
        }
        read <- function(name) tryCatch(output[[paste0(enrich_prefix, name)]], error = function(e) {
            if (inherits(e, c("shiny.silent.error", "shiny.output.cancel"))) "<req>" else paste("<error>", text(conditionMessage(e)))
        })
        for (ui in list(list(mod_prot_da_ui, "differential_expression-"), list(mod_prot_enrich_ui, enrich_prefix))) {
            ui_defaults <- omicsUiInputDefaults(ui[[1L]]("x"))
            names(ui_defaults) <- paste0(ui[[2L]], sub("^x-", "", names(ui_defaults)))
            put(ui_defaults)
        }
        capture <- function(step) {
            session$flushReact()
            apply_updates()
            bus <- setup$bus()
            manager <- shiny::isolate(bus$state_manager)
            out_names <- sort(names(session$.__enclos_env__$private$outs))
            out_names <- sub(paste0("^", enrich_prefix), "", out_names[startsWith(out_names, enrich_prefix)])
            # A download output renders its file when read: it is recorded by name (the download step reads it).
            shown <- vapply(out_names, function(name) {
                if (startsWith(name, "download_")) return("<download>")
                guiEntryOutputValue(read(name), setup$roots)
            }, "")
            returned <- session$userData$enrich
            files <- unlist(guiEntryArtifacts(setup$files_root, setup$roots))
            record <- c(list(
                outputs = out_names,
                values = shown,
                notes = vapply(log$notes, text, "", USE.NAMES = FALSE),
                updates = vapply(log$updates, text, "", USE.NAMES = FALSE),
                progress = vapply(log$progress, text, "", USE.NAMES = FALSE),
                tab_status = vapply(shiny::isolate(bus$tab_status), as.character, ""),
                history = guiEntryHistory(manager),
                states = omicEnrichW1States(manager, setup$roots),
                product = omicEnrichW1Product(shiny::isolate(returned$enrichment_product), manager,
                    shiny::isolate(bus$da_analysis_results_list)),
                resources = omicEnrichW1Resources(manager),
                kernel = acc$kernel,
                services = acc$services,
                bus = omicEnrichW1Bus(bus, setup$roots),
                display = omicEnrichW1Display(shiny::isolate(returned$all_enrichment_results), manager),
                files = stats::setNames(as.character(files), names(files)),
                returned = omicEnrichW1Returned(returned, setup$roots),
                errors = vapply(acc$errors, text, "", USE.NAMES = FALSE)), acc$extra)
            # A field the step left as the previous step had it is held as "<unchanged>".
            for (field in c("tab_status", "history", "states", "product", "resources", "bus", "display", "files", "returned")) {
                value <- record[[field]]
                if (identical(value, acc$previous[[field]])) record[[field]] <- "<unchanged>"
                acc$previous[[field]] <- value
            }
            log$notes <- character()
            log$updates <- character()
            log$progress <- character()
            acc$kernel <- character()
            acc$services <- character()
            acc$errors <- character()
            acc$extra <- list()
            acc$captures[[step]] <- record
        }
        ctx <- list(set = set, set_da = set_da, read = read, bus = setup$bus, roots = setup$roots, tab = setup$tab,
            paths = setup$paths, returned = function() session$userData$enrich,
            keep = function(name, value) acc$extra[[name]] <- value)
        for (step in names(steps)) {
            steps[[step]](ctx)
            capture(step)
        }
    }))))
    acc$captures
}

#' The fixture: KIT-106 W1's proteomics shape. A DIA owner with the project's paths holding MCI-009's object (row names
#' dropped: the DIA reference's tables refuse them) as `correlation_filtered`, the bus as the proteomics shell keeps it
#' (mod_proteomics.R; `omicGuiProteomicsShell`), and the shell's tab selection on `qc`. The DA tab and the enrichment
#' tab are mounted under their module ids with the selection, as the shell mounts them (mod_omic_workflow.R:221-222).
#' @param bus_extra Bus fields set at the start (the mixed-species analysis, the taxon).
omicEnrichW1Setup <- function(bus_extra = list(), env = parent.frame()) {
    root <- normalizePath(withr::local_tempdir(pattern = "kit107-enrich-", .local_envir = env))
    results <- file.path(root, "results")
    paths <- list(source_dir = file.path(root, "source"), results_dir = results, da_output_dir = file.path(results, "da_proteins"),
        pathway_dir = file.path(results, "pathway_enrichment"), publication_graphs_dir = file.path(root, "publication_graphs"))
    for (path in paths) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    design <- module_ci_prot_da_design()
    rownames(design) <- NULL
    quant <- module_ci_prot_da_table()
    rownames(quant) <- NULL
    object <- module_ci_prot_da_object(design = design, protein_quant_table = quant)
    tab <- shiny::reactiveVal("qc")
    holder <- new.env(parent = emptyenv())
    server <- function(acc) {
        function(input, output, session) {
            manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "kit107-enrich", "DIA"),
                WorkflowState$new(audit_enabled = FALSE))
            session$onSessionEnded(function() closeDIAStateAdapter(manager))
            manager$.context$reference_states <- TRUE
            manager$.context$experiment_paths <- paths
            session$userData$dia_capacity_manager <- manager
            manager$saveState("correlation_filtered", object, list(), "MCI-009 fixture state")
            fields <- list(state_manager = manager, contrasts_tbl = module_ci_prot_da_contrasts("two_group"),
                tab_status = list(normalization = "complete", differential_expression = "pending",
                    enrichment_analysis = "disabled"),
                da_analysis_results_list = NULL, da_ui_params = NULL, uniprot_dat_cln = NULL, taxon_id = NULL,
                mixed_species_analysis = NULL, enrichment_analysis_results = NULL, enrichment_ui_params = NULL,
                enrichment_organism_filter = NULL, state_update_trigger = NULL)
            fields[names(bus_extra)] <- bus_extra
            holder$bus <- do.call(shiny::reactiveValues, fields)
            mod_prot_da_server("differential_expression", holder$bus, paths, "proteomics", "kit107",
                selected_tab = function() tab())
            session$userData$enrich <- mod_prot_enrich_server("enrichment_analysis", holder$bus, paths, "proteomics",
                "kit107", selected_tab = function() tab())
        }
    }
    list(server = server, bus = function() holder$bus, roots = root, files_root = root, tab = tab, paths = paths)
}

#' The steps the cases share: the DA (opened, run), the enrichment tab opened, the contrast selected, Run.
omicEnrichW1Da <- function(ctx) {
    ctx$tab("da")
    ctx$set_da(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)
    ctx$set_da(run_da_analysis = 1L)
}
omicEnrichW1Open <- function(ctx) ctx$tab("enrichment")
omicEnrichW1Run <- function(n) function(ctx) ctx$set(run_enrichment_analysis = n)

#' The download as the browser saves it: the file's name (its date masked) and each member's name and digest (the
#' summary's date line masked).
omicEnrichW1Download <- function(ctx) {
    path <- ctx$read("download_enrichment_results")
    if (!(is.character(path) && length(path) == 1L && file.exists(path))) {
        ctx$keep("download", paste("<no file>", paste(path, collapse = " ")))
        return(invisible())
    }
    members <- utils::unzip(path, list = TRUE)$Name
    out <- tempfile("kit107-download-")
    utils::unzip(path, exdir = out)
    digests <- vapply(sort(members), function(member) {
        lines <- readLines(file.path(out, member), warn = FALSE)
        lines <- sub("^Date: .*$", "Date: <TIME>", lines)
        digest::digest(guiEntryMask(lines, ctx$roots), algo = "sha256", serialize = FALSE)
    }, "")
    unlink(out, recursive = TRUE)
    ctx$keep("download", c(file = guiEntryMask(basename(path)), digests))
}

#' The steps of one backend's case: the tab, the contrast, Run, each direction filter, the download.
omicEnrichW1BackendSteps <- function(taxon, method_tab) {
    filter <- paste0(if (identical(method_tab, "gprofiler2")) "gprofiler" else "clusterprofiler", "_direction_filter")
    list(
        start = function(ctx) NULL,
        da = omicEnrichW1Da,
        open = omicEnrichW1Open,
        select = function(ctx) ctx$set(organism_taxid = taxon, enrichment_method_tabs = method_tab, selected_contrast = "B_vs_A",
            up_cutoff = 0, down_cutoff = 0, q_cutoff = 0.05),
        run = omicEnrichW1Run(1L),
        filter_up = function(ctx) do.call(ctx$set, stats::setNames(list("up"), filter)),
        filter_down = function(ctx) do.call(ctx$set, stats::setNames(list("down"), filter)),
        filter_all = function(ctx) do.call(ctx$set, stats::setNames(list("all"), filter)),
        download = omicEnrichW1Download)
}

#' Each case's sessions: a fixture, its steps and its kernel. A case with several sessions names its steps
#' `<session>/<step>`.
omicEnrichW1CaseSpec <- function(case, env) {
    one <- function(steps, live = FALSE, bus_extra = list()) list(list(setup = omicEnrichW1Setup(bus_extra, env = env),
        steps = steps, live = live))
    select <- function(taxon = "9606", ...) function(ctx) ctx$set(organism_taxid = taxon, enrichment_method_tabs = "gprofiler2",
        selected_contrast = "B_vs_A", up_cutoff = 0, down_cutoff = 0, q_cutoff = 0.05, ...)
    switch(case,
        gprofiler = one(omicEnrichW1BackendSteps("9606", "gprofiler2")),
        clusterprofiler = one(omicEnrichW1BackendSteps("999999", "clusterprofiler")),
        kernel_gprofiler = one(omicEnrichW1BackendSteps("9606", "gprofiler2")[c("da", "open", "select", "run")], live = TRUE),
        kernel_clusterprofiler = one(omicEnrichW1BackendSteps("999999", "clusterprofiler")[c("da", "open", "select", "run")],
            live = TRUE),
        organism_filter = c(
            # The import's organism mapping is empty: the warning, an unfiltered run.
            stats::setNames(one(list(da = omicEnrichW1Da, open = omicEnrichW1Open, select = select(enable_organism_filter = TRUE),
                run = omicEnrichW1Run(1L)), bus_extra = list(mixed_species_analysis = list(enabled = FALSE,
                organism_mapping = data.frame(uniprot_acc = character(), taxon_id = character())))), "no_mapping"),
            # The mixed-species observer, then a mapping found from the annotations: today's `clean_acc` failure.
            stats::setNames(one(list(da = omicEnrichW1Da, open = omicEnrichW1Open,
                mixed = function(ctx) {
                    bus <- ctx$bus()
                    bus$mixed_species_analysis <- list(enabled = TRUE, selected_organism = "Homo sapiens", organism_mapping = NULL)
                },
                select = select(),
                run = omicEnrichW1Run(1L))), "mapped")),
        refusals = c(
            stats::setNames(one(list(
                tab_before_da = omicEnrichW1Open,
                run_before_da = function(ctx) {
                    select()(ctx)
                    omicEnrichW1Run(1L)(ctx)
                },
                da = omicEnrichW1Da,
                open = omicEnrichW1Open,
                run_no_contrast = function(ctx) {
                    ctx$set(selected_contrast = "")
                    omicEnrichW1Run(2L)(ctx)
                },
                run_taxon_abc = function(ctx) {
                    select(taxon = "abc")(ctx)
                    omicEnrichW1Run(3L)(ctx)
                })), "pre_da"),
            stats::setNames(one(list(da = omicEnrichW1Da, open = omicEnrichW1Open,
                run_q_cutoff = function(ctx) {
                    select()(ctx)
                    ctx$set(q_cutoff = 1.5)
                    omicEnrichW1Run(1L)(ctx)
                })), "q_cutoff"),
            stats::setNames(one(list(da = omicEnrichW1Da, open = omicEnrichW1Open,
                run_table_removed = function(ctx) {
                    unlink(file.path(ctx$paths$da_output_dir, "da_proteins_groupB-groupA_long_annot.tsv"))
                    select()(ctx)
                    omicEnrichW1Run(1L)(ctx)
                })), "table_removed")),
        restore = one(list(da = omicEnrichW1Da, open = omicEnrichW1Open, select = select(), run = omicEnrichW1Run(1L),
            restore = function(ctx) {
                # KIT-024's measure: a holder with no results over the owner that holds `enrichment_completed`.
                holder <- shiny::reactiveValues(all_enrichment_results = NULL, enrichment_results_full = NULL,
                    analysis_complete = FALSE)
                restored <- restoreProtEnrichResultsFromState(holder, ctx$bus())
                again <- restoreProtEnrichResultsFromState(holder, ctx$bus())
                manager <- shiny::isolate(ctx$bus()$state_manager)
                ctx$keep("restore", c(restored = restored, again = again,
                    analysis_complete = shiny::isolate(holder$analysis_complete),
                    enrichment_results_full = omicEnrichW1Digest(shiny::isolate(holder$enrichment_results_full), ctx$roots),
                    omicEnrichW1Display(shiny::isolate(holder$all_enrichment_results), manager)))
            })),
        prefill = one(list(start = function(ctx) NULL,
            prefill = function(ctx) {
                bus <- ctx$bus()
                bus$taxon_id <- "10090"
            })))
}

omicEnrichW1Runs <- new.env(parent = emptyenv())

#' One case's captures, run once per process (each session has its own fixture).
omicEnrichW1Case <- function(case) {
    if (is.null(omicEnrichW1Runs[[case]])) {
        env <- environment()
        sessions <- omicEnrichW1CaseSpec(case, env)
        captures <- list()
        for (session_name in names(sessions) %||% rep("", length(sessions))) {
            spec <- if (nzchar(session_name)) sessions[[session_name]] else sessions[[1L]]
            seen <- omicEnrichW1Drive(spec$setup, spec$steps, live = spec$live)
            if (nzchar(session_name)) names(seen) <- paste0(session_name, "/", names(seen))
            captures <- c(captures, seen)
        }
        omicEnrichW1Runs[[case]] <- captures
        if (!is.null(omicEnrichW1DumpPath())) {
            cat(sprintf("omicEnrichW1Frozen[[\"%s\"]] <- %s\n\n", case, omicEnrichW1Code(captures)),
                file = omicEnrichW1DumpPath(), append = TRUE)
        }
    }
    omicEnrichW1Runs[[case]]
}

# The frozen captures, per case (the dump's entries, pasted).
# OMIC_ENRICH_FROZEN_START
# KIT-107 W3 re-freeze (once; kit-107/w3/w1-refreeze.md), from two byte-identical dumps in separate processes on the W3
# tree, of exactly the values its fixes move (every other value is as W1 froze it at 8e7294b6):
#   F107-01 (the tab keeps the DA manifest it holds across runs): `returned$da_results_data` after a Run is the
#     manifest (00c45cf144cf6614) where the start commit replaced it with `manifest$results` (8a3eda2205cf9527): the run
#     steps of gprofiler, clusterprofiler, kernel_*, organism_filter no_mapping, restore and table_removed, and the
#     failed runs' `returned` (organism_filter mapped, refusals run_taxon_abc and run_q_cutoff) now "<unchanged>";
#   F107-03 (the lead's decision: a missing DA table refuses the run before any bus write, state or product):
#     refusals table_removed/run_table_removed;
#   F107-05 (the live kernel's results made portable before the commit): kernel_gprofiler and kernel_clusterprofiler
#     run steps' history, states and store state-entry counts (the commit the start commit refused now happens).
omicEnrichW1Frozen <- list()
omicEnrichW1Frozen[["gprofiler"]] <- list(
    start = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<req>",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "No contrasts available.\nComplete differential expression\nanalysis first.",
            download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<req>",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        tab_status = c(
            normalization = "complete", differential_expression = "pending", enrichment_analysis = "disabled"),
        history = c("initial", "correlation_filtered"),
        states = c(initial = "71557d1c8bac9bbe", correlation_filtered = "bf8ffd456d2db4fb"),
        product = "<none>",
        resources = character(0),
        kernel = character(0),
        services = character(0),
        bus = c(
            enrichment_analysis_results = "<NULL>", `enrichment_ui_params$` = "<NULL>",
            `enrichment_organism_filter$` = "<NULL>", uniprot_dat_cln = "<NULL>"),
        display = "<empty>",
        files = c(
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "4",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "3"),
        returned = c(
            all_enrichment_results = "<list> 28df520b992b4a4b", analysis_complete = "FALSE",
            analysis_method = "gprofiler2", clusterprofiler_results = "<NULL>",
            contrasts_available = "<NULL>", current_contrast_results = "<list> 28df520b992b4a4b",
            current_s4_object = "<NULL>", da_results_data = "<NULL>",
            enrichment_plots = "<list> 28df520b992b4a4b", enrichment_results = "<NULL>",
            gprofiler_results = "<NULL>", organism_supported = "TRUE", stringdb_results = "<NULL>"),
        errors = character(0)),
    da = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<req>",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<req>",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = c(
            "[default NULL da_working] Running differential abundance analysis...",
            "[message 5 -] Differential abundance analysis completed successfully!"),
        updates = c(
            "differential_expression-formula_string: <same> -> ~ 0 + group",
            "differential_expression-volcano_contrast: groupB-groupA -> <none>",
            "differential_expression-heatmap_contrast: groupB-groupA -> <none>",
            "differential_expression-table_contrast: groupB-groupA -> <none>",
            "differential_expression-volcano_contrast: B_vs_A -> <none>",
            "differential_expression-heatmap_contrast: B_vs_A -> <none>",
            "differential_expression-table_contrast: B_vs_A -> <none>",
            "enrichment_analysis-selected_contrast: B_vs_A -> <none>"),
        progress = character(0),
        tab_status = c(
            normalization = "complete", differential_expression = "complete", enrichment_analysis = "pending",
            differential_abundance = "complete"),
        history = "<unchanged>",
        states = c(initial = "71557d1c8bac9bbe", correlation_filtered = "5acf333d81fb9c8a"),
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = c(
            `<store entries>..workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "15",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "8",
            `<store entries>..workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `publication_graphs/Interactive_Volcano_Plots/B_vs_A.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "8a7cae2958352a9afa06dadd17250eacf2c56a50583e8b0b1c17146d9991c109",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.png` = "<present>",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.tsv` = "fbe4b6a921c2abe20abdb68f3fae2a34f3a45e70062415d3da8acebed58e73c6",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.xlsx` = "<present>"),
        returned = c(
            all_enrichment_results = "<list> 28df520b992b4a4b", analysis_complete = "FALSE",
            analysis_method = "gprofiler2", clusterprofiler_results = "<NULL>",
            contrasts_available = "B_vs_A", current_contrast_results = "<list> 28df520b992b4a4b",
            current_s4_object = "<ref>", da_results_data = "<list> 00c45cf144cf6614",
            enrichment_plots = "<list> 28df520b992b4a4b", enrichment_results = "<NULL>",
            gprofiler_results = "<NULL>", organism_supported = "TRUE", stringdb_results = "<NULL>"),
        errors = character(0)),
    open = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<req>",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<req>",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = "enrichment_analysis-selected_contrast: B_vs_A -> <none>",
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)),
    select = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)),
    run = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<text> d9625c0c8f0d2e2d05bf97481af72158e56fcbcc0cb47a9ac59b41f286e9dd15",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[OK] Analysis Complete\n\nMethod: gprofiler2\n\nContrast: B_vs_A\n\nUp log2FC cutoff: 0.0\n\nDown log2FC cutoff: 0.0\n\nQ-value cutoff: 0.050\n\nOrganism: Homo sapiens\n\n\nResults Available:\n* gprofiler2: 10 terms\n* clusterProfileR: 0 terms\n* STRING-DB: 0 networks\n\n[OK] Results saved to workflow state",
            gprofiler_plot = "<text> d9625c0c8f0d2e2d05bf97481af72158e56fcbcc0cb47a9ac59b41f286e9dd15",
            gprofiler_results_table = "<text> 86bdcd67360fa73f5cb9cd2636c1a97294f8c85e4065ab1d82da77e3a81df11a",
            gprofiler_summary_stats = "Total enrichment terms: 10\nUp-regulated pathways: 5\nDown-regulated pathways: 5\n\nResults displayed in table below.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = c(
            "[default NULL enrichment_working] Running enrichment analysis...",
            "[message 5 -] Enrichment analysis completed successfully for contrast: B_vs_A",
            "remove enrichment_working"),
        updates = character(0),
        progress = c(
            "start Running enrichment analysis... 0", "0.2 Transforming DE data...",
            "0.3 Creating DE results S4 object...", "0.5 Running enrichment analysis...",
            "0.8 Storing results...", "1 Complete!"),
        tab_status = c(
            normalization = "complete", differential_expression = "complete",
            enrichment_analysis = "complete", differential_abundance = "complete"),
        history = c("initial", "correlation_filtered", "enrichment_completed"),
        states = c(
            initial = "71557d1c8bac9bbe", correlation_filtered = "5acf333d81fb9c8a",
            enrichment_completed = "61dae2e69108addd"),
        product = c(
            schema = "dia_enrichment_product_v1", service = "gprofiler2:gprofiler2_test_double",
            `parameters$taxon_id` = "9606", `parameters$up_cutoff` = "0", `parameters$down_cutoff` = "0",
            `parameters$q_cutoff` = "0.05", `parameters$exclude_iea` = "FALSE",
            `parameters$correction_method` = "gSCS", `parameters$protein_id_column` = "gene_name",
            `parameters$contrast_names` = "groupB-groupA",
            `results$groupB-groupA` = "<manifest da_proteins_long>", contrasts = "groupB-groupA",
            summaries = "6122a4961d0ad898"),
        resources = "gprofiler2 ; gprofiler2_test_double ; database_version ; enrichment of 1 contrasts, taxon 9606",
        kernel = "kernel=runProtEnrichDeterministicProcessEnrichments ; taxon_id=9606 ; up_cutoff=0 ; down_cutoff=0 ; q_cutoff=0.05 ; exclude_iea=FALSE ; protein_id_column=gene_name ; contrast_names=groupB-groupA ; correction_method=gSCS ; pathway_dir=<ROOT>/results/pathway_enrichment ; da_data=0656c115da47d5db ; da_names=groupB-groupA ; da_rows=6x17 ; da_contrasts=950ddbb7d830df8c ; da_design=bcd78aff6b352869 ; go_annotations=13df57c377de46dc",
        services = character(0),
        bus = c(
            enrichment_analysis_results = "gprofiler_results,clusterprofiler_results,stringdb_results,full_enrichment_results,selected_contrast,analysis_method,parameters bd71e1846eddbd60",
            `enrichment_ui_params$up_log2fc_cutoff` = "0", `enrichment_ui_params$down_log2fc_cutoff` = "0",
            `enrichment_ui_params$q_value_cutoff` = "0.05", `enrichment_ui_params$organism_selected` = "9606",
            `enrichment_ui_params$database_source` = "gprofiler2",
            `enrichment_ui_params$organism_name` = "Homo sapiens",
            `enrichment_ui_params$organism_supported` = "TRUE",
            `enrichment_ui_params$selected_contrast` = "B_vs_A", `enrichment_ui_params$timestamp` = "<time>",
            `enrichment_organism_filter$enabled` = "FALSE",
            `enrichment_organism_filter$filter_applied` = "FALSE",
            `enrichment_organism_filter$target_taxon_id` = "9606",
            `enrichment_organism_filter$proteins_before` = "0",
            `enrichment_organism_filter$proteins_after` = "0",
            `enrichment_organism_filter$proteins_removed` = "0",
            `enrichment_organism_filter$timestamp` = "<time>", uniprot_dat_cln = "5x9 13df57c377de46dc"),
        display = c(
            `groupB-groupA/gprofiler` = "10x24 81a566f72a3062f4",
            `groupB-groupA/gprofiler_table` = "<saved> 10 b2a4adf014f1d53e",
            `groupB-groupA/clusterprofiler` = "<NULL>", `groupB-groupA/stringdb` = "<NULL>"),
        files = c(
            `<store entries>..workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-enrichment-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "23",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "10",
            `<store entries>..workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `publication_graphs/Interactive_Volcano_Plots/B_vs_A.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "8a7cae2958352a9afa06dadd17250eacf2c56a50583e8b0b1c17146d9991c109",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.png` = "<present>",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.tsv` = "fbe4b6a921c2abe20abdb68f3fae2a34f3a45e70062415d3da8acebed58e73c6",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.xlsx` = "<present>",
            `results/pathway_enrichment/groupB-groupA_down_enrichment_results.tsv` = "511024ccd86710326a8c45231c97c82595eabf1732e3bab43fd98c53f3b289d7",
            `results/pathway_enrichment/groupB-groupA_up_enrichment_results.tsv` = "b9b76ea1cdcd7200ad5ca6ab181e831781c08df2df93f223d26cda697738397a"),
        returned = c(
            all_enrichment_results = "<list> 123e2a6cd1f92b90", analysis_complete = "TRUE",
            analysis_method = "gprofiler2", clusterprofiler_results = "<NULL>",
            clusterprofiler_results_table = "<NULL>", contrasts_available = "B_vs_A",
            current_contrast_results = "<list> 28df520b992b4a4b", current_s4_object = "<ref>",
            da_results_data = "<list> 00c45cf144cf6614", enrichment_plots = "<list> 28df520b992b4a4b",
            enrichment_product = "<ref>", enrichment_results = "<NULL>",
            enrichment_results_full = "<EnrichmentResults> 61dae2e69108addd",
            enrichment_service = "gprofiler2:gprofiler2_test_double",
            gprofiler_results = "<data.frame> 81a566f72a3062f4", gprofiler_results_table = "<ref>",
            organism_supported = "TRUE", stringdb_results = "<NULL>"),
        errors = character(0)),
    filter_up = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<text> d9625c0c8f0d2e2d05bf97481af72158e56fcbcc0cb47a9ac59b41f286e9dd15",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[OK] Analysis Complete\n\nMethod: gprofiler2\n\nContrast: B_vs_A\n\nUp log2FC cutoff: 0.0\n\nDown log2FC cutoff: 0.0\n\nQ-value cutoff: 0.050\n\nOrganism: Homo sapiens\n\n\nResults Available:\n* gprofiler2: 10 terms\n* clusterProfileR: 0 terms\n* STRING-DB: 0 networks\n\n[OK] Results saved to workflow state",
            gprofiler_plot = "<text> d9625c0c8f0d2e2d05bf97481af72158e56fcbcc0cb47a9ac59b41f286e9dd15",
            gprofiler_results_table = "<text> 86bdcd67360fa73f5cb9cd2636c1a97294f8c85e4065ab1d82da77e3a81df11a",
            gprofiler_summary_stats = "Showing 5 up-regulated pathways\n\nResults displayed in table below.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)),
    filter_down = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<text> d9625c0c8f0d2e2d05bf97481af72158e56fcbcc0cb47a9ac59b41f286e9dd15",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[OK] Analysis Complete\n\nMethod: gprofiler2\n\nContrast: B_vs_A\n\nUp log2FC cutoff: 0.0\n\nDown log2FC cutoff: 0.0\n\nQ-value cutoff: 0.050\n\nOrganism: Homo sapiens\n\n\nResults Available:\n* gprofiler2: 10 terms\n* clusterProfileR: 0 terms\n* STRING-DB: 0 networks\n\n[OK] Results saved to workflow state",
            gprofiler_plot = "<text> 2298435056ffd0ee14d729802221eae1cd2986822050d9ab20040d2ab90a9a5b",
            gprofiler_results_table = "<text> 86bdcd67360fa73f5cb9cd2636c1a97294f8c85e4065ab1d82da77e3a81df11a",
            gprofiler_summary_stats = "Showing 5 down-regulated pathways\n\nResults displayed in table below.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)),
    filter_all = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<text> d9625c0c8f0d2e2d05bf97481af72158e56fcbcc0cb47a9ac59b41f286e9dd15",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[OK] Analysis Complete\n\nMethod: gprofiler2\n\nContrast: B_vs_A\n\nUp log2FC cutoff: 0.0\n\nDown log2FC cutoff: 0.0\n\nQ-value cutoff: 0.050\n\nOrganism: Homo sapiens\n\n\nResults Available:\n* gprofiler2: 10 terms\n* clusterProfileR: 0 terms\n* STRING-DB: 0 networks\n\n[OK] Results saved to workflow state",
            gprofiler_plot = "<text> d9625c0c8f0d2e2d05bf97481af72158e56fcbcc0cb47a9ac59b41f286e9dd15",
            gprofiler_results_table = "<text> 86bdcd67360fa73f5cb9cd2636c1a97294f8c85e4065ab1d82da77e3a81df11a",
            gprofiler_summary_stats = "Total enrichment terms: 10\nUp-regulated pathways: 5\nDown-regulated pathways: 5\n\nResults displayed in table below.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)),
    download = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<text> d9625c0c8f0d2e2d05bf97481af72158e56fcbcc0cb47a9ac59b41f286e9dd15",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[OK] Analysis Complete\n\nMethod: gprofiler2\n\nContrast: B_vs_A\n\nUp log2FC cutoff: 0.0\n\nDown log2FC cutoff: 0.0\n\nQ-value cutoff: 0.050\n\nOrganism: Homo sapiens\n\n\nResults Available:\n* gprofiler2: 10 terms\n* clusterProfileR: 0 terms\n* STRING-DB: 0 networks\n\n[OK] Results saved to workflow state",
            gprofiler_plot = "<text> d9625c0c8f0d2e2d05bf97481af72158e56fcbcc0cb47a9ac59b41f286e9dd15",
            gprofiler_results_table = "<text> 86bdcd67360fa73f5cb9cd2636c1a97294f8c85e4065ab1d82da77e3a81df11a",
            gprofiler_summary_stats = "Total enrichment terms: 10\nUp-regulated pathways: 5\nDown-regulated pathways: 5\n\nResults displayed in table below.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0),
        download = c(
            file = "Enrichment_results_B_vs_A_<DATE>.zip",
            enrichment_analysis_summary.txt = "13ce7979ab166ba5faa78f9b2904c2446c580482a9a541e9afd429d04c59161f",
            gprofiler2_results.tsv = "4182e8063b318d00af302551c715e661a12e048f244dde5ac6a76b43aac4066d")))

omicEnrichW1Frozen[["clusterprofiler"]] <- list(
    start = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<req>",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "No contrasts available.\nComplete differential expression\nanalysis first.",
            download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<req>",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        tab_status = c(
            normalization = "complete", differential_expression = "pending", enrichment_analysis = "disabled"),
        history = c("initial", "correlation_filtered"),
        states = c(initial = "71557d1c8bac9bbe", correlation_filtered = "bf8ffd456d2db4fb"),
        product = "<none>",
        resources = character(0),
        kernel = character(0),
        services = character(0),
        bus = c(
            enrichment_analysis_results = "<NULL>", `enrichment_ui_params$` = "<NULL>",
            `enrichment_organism_filter$` = "<NULL>", uniprot_dat_cln = "<NULL>"),
        display = "<empty>",
        files = c(
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "4",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "3"),
        returned = c(
            all_enrichment_results = "<list> 28df520b992b4a4b", analysis_complete = "FALSE",
            analysis_method = "gprofiler2", clusterprofiler_results = "<NULL>",
            contrasts_available = "<NULL>", current_contrast_results = "<list> 28df520b992b4a4b",
            current_s4_object = "<NULL>", da_results_data = "<NULL>",
            enrichment_plots = "<list> 28df520b992b4a4b", enrichment_results = "<NULL>",
            gprofiler_results = "<NULL>", organism_supported = "TRUE", stringdb_results = "<NULL>"),
        errors = character(0)),
    da = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<req>",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<req>",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = c(
            "[default NULL da_working] Running differential abundance analysis...",
            "[message 5 -] Differential abundance analysis completed successfully!"),
        updates = c(
            "differential_expression-formula_string: <same> -> ~ 0 + group",
            "differential_expression-volcano_contrast: groupB-groupA -> <none>",
            "differential_expression-heatmap_contrast: groupB-groupA -> <none>",
            "differential_expression-table_contrast: groupB-groupA -> <none>",
            "differential_expression-volcano_contrast: B_vs_A -> <none>",
            "differential_expression-heatmap_contrast: B_vs_A -> <none>",
            "differential_expression-table_contrast: B_vs_A -> <none>",
            "enrichment_analysis-selected_contrast: B_vs_A -> <none>"),
        progress = character(0),
        tab_status = c(
            normalization = "complete", differential_expression = "complete", enrichment_analysis = "pending",
            differential_abundance = "complete"),
        history = "<unchanged>",
        states = c(initial = "71557d1c8bac9bbe", correlation_filtered = "5acf333d81fb9c8a"),
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = c(
            `<store entries>..workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "15",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "8",
            `<store entries>..workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `publication_graphs/Interactive_Volcano_Plots/B_vs_A.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "8a7cae2958352a9afa06dadd17250eacf2c56a50583e8b0b1c17146d9991c109",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.png` = "<present>",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.tsv` = "fbe4b6a921c2abe20abdb68f3fae2a34f3a45e70062415d3da8acebed58e73c6",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.xlsx` = "<present>"),
        returned = c(
            all_enrichment_results = "<list> 28df520b992b4a4b", analysis_complete = "FALSE",
            analysis_method = "gprofiler2", clusterprofiler_results = "<NULL>",
            contrasts_available = "B_vs_A", current_contrast_results = "<list> 28df520b992b4a4b",
            current_s4_object = "<ref>", da_results_data = "<list> 00c45cf144cf6614",
            enrichment_plots = "<list> 28df520b992b4a4b", enrichment_results = "<NULL>",
            gprofiler_results = "<NULL>", organism_supported = "TRUE", stringdb_results = "<NULL>"),
        errors = character(0)),
    open = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<req>",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<req>",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = "enrichment_analysis-selected_contrast: B_vs_A -> <none>",
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)),
    select = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[WARNING] CUSTOM ORGANISM\n Method: clusterprofiler\n Organism: Taxon ID 999999\n Using UniProt GO annotations",
            clusterprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = c(
            all_enrichment_results = "<list> 28df520b992b4a4b", analysis_complete = "FALSE",
            analysis_method = "clusterprofiler", clusterprofiler_results = "<NULL>",
            contrasts_available = "B_vs_A", current_contrast_results = "<list> 28df520b992b4a4b",
            current_s4_object = "<ref>", da_results_data = "<list> 00c45cf144cf6614",
            enrichment_plots = "<list> 28df520b992b4a4b", enrichment_results = "<NULL>",
            gprofiler_results = "<NULL>", organism_supported = "FALSE", stringdb_results = "<NULL>"),
        errors = character(0)),
    run = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[WARNING] CUSTOM ORGANISM\n Method: clusterprofiler\n Organism: Taxon ID 999999\n Using UniProt GO annotations",
            clusterprofiler_plot = "<text> a6d285083d30b88bb27483536326613cc71189d7ef8563ed2657465fe8bb6d22",
            clusterprofiler_results_table = "<text> f30b2e0be6ce4389129c8badc9e79a89c56d017116ba48b18b1f42ba73954083",
            clusterprofiler_summary_stats = "Total GO terms: 6\nUp-regulated: 3\nDown-regulated: 3\n\nResults displayed in table below.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[OK] Analysis Complete\n\nMethod: clusterprofiler\n\nContrast: B_vs_A\n\nUp log2FC cutoff: 0.0\n\nDown log2FC cutoff: 0.0\n\nQ-value cutoff: 0.050\n\nOrganism: Taxon ID 999999\n\n\nResults Available:\n* gprofiler2: 0 terms\n* clusterProfileR: 6 terms\n* STRING-DB: 0 networks\n\n[OK] Results saved to workflow state",
            gprofiler_plot = "<text> a6d285083d30b88bb27483536326613cc71189d7ef8563ed2657465fe8bb6d22",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = c(
            "[default NULL enrichment_working] Running enrichment analysis...",
            "[message 5 -] Enrichment analysis completed successfully for contrast: B_vs_A",
            "remove enrichment_working"),
        updates = character(0),
        progress = c(
            "start Running enrichment analysis... 0", "0.2 Transforming DE data...",
            "0.3 Creating DE results S4 object...", "0.5 Running enrichment analysis...",
            "0.8 Storing results...", "1 Complete!"),
        tab_status = c(
            normalization = "complete", differential_expression = "complete",
            enrichment_analysis = "complete", differential_abundance = "complete"),
        history = c("initial", "correlation_filtered", "enrichment_completed"),
        states = c(
            initial = "71557d1c8bac9bbe", correlation_filtered = "5acf333d81fb9c8a",
            enrichment_completed = "a62a41b234214817"),
        product = c(
            schema = "dia_enrichment_product_v1", service = "unavailable", `parameters$taxon_id` = "999999",
            `parameters$up_cutoff` = "0", `parameters$down_cutoff` = "0", `parameters$q_cutoff` = "0.05",
            `parameters$exclude_iea` = "FALSE", `parameters$correction_method` = "gSCS",
            `parameters$protein_id_column` = "uniprot_acc", `parameters$contrast_names` = "groupB-groupA",
            `results$groupB-groupA` = "<manifest da_proteins_long>", contrasts = "groupB-groupA",
            summaries = "5d6736b697306a94"),
        resources = "<unchanged>",
        kernel = "kernel=runProtEnrichDeterministicProcessEnrichments ; taxon_id=999999 ; up_cutoff=0 ; down_cutoff=0 ; q_cutoff=0.05 ; exclude_iea=FALSE ; protein_id_column=uniprot_acc ; contrast_names=groupB-groupA ; correction_method=gSCS ; pathway_dir=<ROOT>/results/pathway_enrichment ; da_data=0656c115da47d5db ; da_names=groupB-groupA ; da_rows=6x17 ; da_contrasts=950ddbb7d830df8c ; da_design=bcd78aff6b352869 ; go_annotations=0218200d5d349843",
        services = character(0),
        bus = c(
            enrichment_analysis_results = "gprofiler_results,clusterprofiler_results,stringdb_results,full_enrichment_results,selected_contrast,analysis_method,parameters 54dc05374647be6e",
            `enrichment_ui_params$up_log2fc_cutoff` = "0", `enrichment_ui_params$down_log2fc_cutoff` = "0",
            `enrichment_ui_params$q_value_cutoff` = "0.05",
            `enrichment_ui_params$organism_selected` = "999999",
            `enrichment_ui_params$database_source` = "clusterprofiler",
            `enrichment_ui_params$organism_name` = "Taxon ID 999999",
            `enrichment_ui_params$organism_supported` = "FALSE",
            `enrichment_ui_params$selected_contrast` = "B_vs_A", `enrichment_ui_params$timestamp` = "<time>",
            `enrichment_organism_filter$enabled` = "FALSE",
            `enrichment_organism_filter$filter_applied` = "FALSE",
            `enrichment_organism_filter$target_taxon_id` = "999999",
            `enrichment_organism_filter$proteins_before` = "0",
            `enrichment_organism_filter$proteins_after` = "0",
            `enrichment_organism_filter$proteins_removed` = "0",
            `enrichment_organism_filter$timestamp` = "<time>", uniprot_dat_cln = "5x9 0218200d5d349843"),
        display = c(
            `groupB-groupA/gprofiler` = "<NULL>", `groupB-groupA/clusterprofiler` = "6x12 6be7681575625a6d",
            `groupB-groupA/clusterprofiler_table` = "<saved> 6 6be7681575625a6d",
            `groupB-groupA/stringdb` = "<NULL>"),
        files = c(
            `<store entries>..workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-enrichment-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "22",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "10",
            `<store entries>..workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `publication_graphs/Interactive_Volcano_Plots/B_vs_A.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "8a7cae2958352a9afa06dadd17250eacf2c56a50583e8b0b1c17146d9991c109",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.png` = "<present>",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.tsv` = "fbe4b6a921c2abe20abdb68f3fae2a34f3a45e70062415d3da8acebed58e73c6",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.xlsx` = "<present>",
            `results/pathway_enrichment/groupB-groupA_down_enrichment_results.tsv` = "cb1c7be929956b9c2cfb8932d8136fcf92e6abd76b42c2a89c5b57d12b1ed8b0",
            `results/pathway_enrichment/groupB-groupA_up_enrichment_results.tsv` = "077ca85b815eecc33439f47e14ff20ceb7fb4d0ee1bd79a44d27a2f13f7d3dad"),
        returned = c(
            all_enrichment_results = "<list> f3e983fdeac04383", analysis_complete = "TRUE",
            analysis_method = "clusterprofiler", clusterprofiler_results = "<data.frame> 6be7681575625a6d",
            clusterprofiler_results_table = "<ref>", contrasts_available = "B_vs_A",
            current_contrast_results = "<list> 28df520b992b4a4b", current_s4_object = "<ref>",
            da_results_data = "<list> 00c45cf144cf6614", enrichment_plots = "<list> 28df520b992b4a4b",
            enrichment_product = "<ref>", enrichment_results = "<NULL>",
            enrichment_results_full = "<EnrichmentResults> a62a41b234214817",
            enrichment_service = "unavailable", gprofiler_results = "<NULL>",
            gprofiler_results_table = "<NULL>", organism_supported = "FALSE", stringdb_results = "<NULL>"),
        errors = character(0)),
    filter_up = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[WARNING] CUSTOM ORGANISM\n Method: clusterprofiler\n Organism: Taxon ID 999999\n Using UniProt GO annotations",
            clusterprofiler_plot = "<text> a6d285083d30b88bb27483536326613cc71189d7ef8563ed2657465fe8bb6d22",
            clusterprofiler_results_table = "<text> f30b2e0be6ce4389129c8badc9e79a89c56d017116ba48b18b1f42ba73954083",
            clusterprofiler_summary_stats = "Showing 3 up-regulated GO terms\n\nResults displayed in table below.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[OK] Analysis Complete\n\nMethod: clusterprofiler\n\nContrast: B_vs_A\n\nUp log2FC cutoff: 0.0\n\nDown log2FC cutoff: 0.0\n\nQ-value cutoff: 0.050\n\nOrganism: Taxon ID 999999\n\n\nResults Available:\n* gprofiler2: 0 terms\n* clusterProfileR: 6 terms\n* STRING-DB: 0 networks\n\n[OK] Results saved to workflow state",
            gprofiler_plot = "<text> a6d285083d30b88bb27483536326613cc71189d7ef8563ed2657465fe8bb6d22",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)),
    filter_down = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[WARNING] CUSTOM ORGANISM\n Method: clusterprofiler\n Organism: Taxon ID 999999\n Using UniProt GO annotations",
            clusterprofiler_plot = "<text> 1b2e9d6b289a56d87b7b46d31944bc7c346e0b7b1b6b4848c81a2d578a9e4188",
            clusterprofiler_results_table = "<text> f30b2e0be6ce4389129c8badc9e79a89c56d017116ba48b18b1f42ba73954083",
            clusterprofiler_summary_stats = "Showing 3 down-regulated GO terms\n\nResults displayed in table below.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[OK] Analysis Complete\n\nMethod: clusterprofiler\n\nContrast: B_vs_A\n\nUp log2FC cutoff: 0.0\n\nDown log2FC cutoff: 0.0\n\nQ-value cutoff: 0.050\n\nOrganism: Taxon ID 999999\n\n\nResults Available:\n* gprofiler2: 0 terms\n* clusterProfileR: 6 terms\n* STRING-DB: 0 networks\n\n[OK] Results saved to workflow state",
            gprofiler_plot = "<text> a6d285083d30b88bb27483536326613cc71189d7ef8563ed2657465fe8bb6d22",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)),
    filter_all = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[WARNING] CUSTOM ORGANISM\n Method: clusterprofiler\n Organism: Taxon ID 999999\n Using UniProt GO annotations",
            clusterprofiler_plot = "<text> a6d285083d30b88bb27483536326613cc71189d7ef8563ed2657465fe8bb6d22",
            clusterprofiler_results_table = "<text> f30b2e0be6ce4389129c8badc9e79a89c56d017116ba48b18b1f42ba73954083",
            clusterprofiler_summary_stats = "Total GO terms: 6\nUp-regulated: 3\nDown-regulated: 3\n\nResults displayed in table below.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[OK] Analysis Complete\n\nMethod: clusterprofiler\n\nContrast: B_vs_A\n\nUp log2FC cutoff: 0.0\n\nDown log2FC cutoff: 0.0\n\nQ-value cutoff: 0.050\n\nOrganism: Taxon ID 999999\n\n\nResults Available:\n* gprofiler2: 0 terms\n* clusterProfileR: 6 terms\n* STRING-DB: 0 networks\n\n[OK] Results saved to workflow state",
            gprofiler_plot = "<text> a6d285083d30b88bb27483536326613cc71189d7ef8563ed2657465fe8bb6d22",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)),
    download = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[WARNING] CUSTOM ORGANISM\n Method: clusterprofiler\n Organism: Taxon ID 999999\n Using UniProt GO annotations",
            clusterprofiler_plot = "<text> a6d285083d30b88bb27483536326613cc71189d7ef8563ed2657465fe8bb6d22",
            clusterprofiler_results_table = "<text> f30b2e0be6ce4389129c8badc9e79a89c56d017116ba48b18b1f42ba73954083",
            clusterprofiler_summary_stats = "Total GO terms: 6\nUp-regulated: 3\nDown-regulated: 3\n\nResults displayed in table below.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[OK] Analysis Complete\n\nMethod: clusterprofiler\n\nContrast: B_vs_A\n\nUp log2FC cutoff: 0.0\n\nDown log2FC cutoff: 0.0\n\nQ-value cutoff: 0.050\n\nOrganism: Taxon ID 999999\n\n\nResults Available:\n* gprofiler2: 0 terms\n* clusterProfileR: 6 terms\n* STRING-DB: 0 networks\n\n[OK] Results saved to workflow state",
            gprofiler_plot = "<text> a6d285083d30b88bb27483536326613cc71189d7ef8563ed2657465fe8bb6d22",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0),
        download = c(
            file = "Enrichment_results_B_vs_A_<DATE>.zip",
            clusterProfileR_results.tsv = "671ae3ed12b3c1b51526dd991311d931bca4bb1e347a492381205a28f264372b",
            enrichment_analysis_summary.txt = "e7a41994feeebd115b4192c86580047e7109cce09c1185306b314c0d71b7dd8b")))

omicEnrichW1Frozen[["kernel_gprofiler"]] <- list(
    da = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<req>",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<req>",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = c(
            "[default NULL da_working] Running differential abundance analysis...",
            "[message 5 -] Differential abundance analysis completed successfully!"),
        updates = c(
            "differential_expression-formula_string: <same> -> ~ 0 + group",
            "differential_expression-volcano_contrast: groupB-groupA -> <none>",
            "differential_expression-heatmap_contrast: groupB-groupA -> <none>",
            "differential_expression-table_contrast: groupB-groupA -> <none>",
            "differential_expression-volcano_contrast: B_vs_A -> <none>",
            "differential_expression-heatmap_contrast: B_vs_A -> <none>",
            "differential_expression-table_contrast: B_vs_A -> <none>",
            "enrichment_analysis-selected_contrast: B_vs_A -> <none>"),
        progress = character(0),
        tab_status = c(
            normalization = "complete", differential_expression = "complete", enrichment_analysis = "pending",
            differential_abundance = "complete"),
        history = c("initial", "correlation_filtered"),
        states = c(initial = "71557d1c8bac9bbe", correlation_filtered = "5acf333d81fb9c8a"),
        product = "<none>",
        resources = character(0),
        kernel = character(0),
        services = character(0),
        bus = c(
            enrichment_analysis_results = "<NULL>", `enrichment_ui_params$` = "<NULL>",
            `enrichment_organism_filter$` = "<NULL>", uniprot_dat_cln = "<NULL>"),
        display = "<empty>",
        files = c(
            `<store entries>..workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "15",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "8",
            `<store entries>..workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `publication_graphs/Interactive_Volcano_Plots/B_vs_A.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "8a7cae2958352a9afa06dadd17250eacf2c56a50583e8b0b1c17146d9991c109",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.png` = "<present>",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.tsv` = "fbe4b6a921c2abe20abdb68f3fae2a34f3a45e70062415d3da8acebed58e73c6",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.xlsx` = "<present>"),
        returned = c(
            all_enrichment_results = "<list> 28df520b992b4a4b", analysis_complete = "FALSE",
            analysis_method = "gprofiler2", clusterprofiler_results = "<NULL>",
            contrasts_available = "B_vs_A", current_contrast_results = "<list> 28df520b992b4a4b",
            current_s4_object = "<ref>", da_results_data = "<list> 00c45cf144cf6614",
            enrichment_plots = "<list> 28df520b992b4a4b", enrichment_results = "<NULL>",
            gprofiler_results = "<NULL>", organism_supported = "TRUE", stringdb_results = "<NULL>"),
        errors = character(0)),
    open = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<req>",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<req>",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = "enrichment_analysis-selected_contrast: B_vs_A -> <none>",
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)),
    select = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)),
    run = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<text> b8031f78fb587ac7a8c8072353259987c98528ec71467c76140c5c037e01a121",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[OK] Analysis Complete\n\nMethod: gprofiler2\n\nContrast: B_vs_A\n\nUp log2FC cutoff: 0.0\n\nDown log2FC cutoff: 0.0\n\nQ-value cutoff: 0.050\n\nOrganism: Homo sapiens\n\n\nResults Available:\n* gprofiler2: 10 terms\n* clusterProfileR: 0 terms\n* STRING-DB: 0 networks\n\n[OK] Results saved to workflow state",
            gprofiler_plot = "<text> b8031f78fb587ac7a8c8072353259987c98528ec71467c76140c5c037e01a121",
            gprofiler_results_table = "<text> 281a1e0abb3f2a5db8699c8ff6fe26dd252483b6a023ba73cf398fc1fa373202",
            gprofiler_summary_stats = "Total enrichment terms: 10\nUp-regulated pathways: 5\nDown-regulated pathways: 5\n\nResults displayed in table below.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = c(
            "[default NULL enrichment_working] Running enrichment analysis...",
            "[message 5 -] Enrichment analysis completed successfully for contrast: B_vs_A",
            "remove enrichment_working"),
        updates = character(0),
        progress = c(
            "start Running enrichment analysis... 0", "0.2 Transforming DE data...",
            "0.3 Creating DE results S4 object...", "0.5 Running enrichment analysis...",
            "0.8 Storing results...", "1 Complete!"),
        tab_status = c(
            normalization = "complete", differential_expression = "complete",
            enrichment_analysis = "complete", differential_abundance = "complete"),
        history = c("initial", "correlation_filtered", "enrichment_completed"),
        states = c(
            initial = "71557d1c8bac9bbe", correlation_filtered = "5acf333d81fb9c8a",
            enrichment_completed = "dbde2059c38a85f0"),
        product = c(
            schema = "dia_enrichment_product_v1", service = "unavailable", `parameters$taxon_id` = "9606",
            `parameters$up_cutoff` = "0", `parameters$down_cutoff` = "0", `parameters$q_cutoff` = "0.05",
            `parameters$exclude_iea` = "FALSE", `parameters$correction_method` = "gSCS",
            `parameters$protein_id_column` = "gene_name", `parameters$contrast_names` = "groupB-groupA",
            `results$groupB-groupA` = "<manifest da_proteins_long>", contrasts = "groupB-groupA",
            summaries = "18ee1ae7f00ccf86"),
        resources = "<unchanged>",
        kernel = "kernel=processEnrichments ; taxon_id=9606 ; up_cutoff=0 ; down_cutoff=0 ; q_cutoff=0.05 ; exclude_iea=FALSE ; protein_id_column=gene_name ; contrast_names=groupB-groupA ; correction_method=gSCS ; pathway_dir=<ROOT>/results/pathway_enrichment ; da_data=0656c115da47d5db ; da_names=groupB-groupA ; da_rows=6x17 ; da_contrasts=950ddbb7d830df8c ; da_design=bcd78aff6b352869 ; go_annotations=13df57c377de46dc",
        services = c(
            "gost ; organism=hsapiens ; sources=GO:BP|GO:MF|GO:CC|KEGG|REAC ; user_threshold=0.05 ; correction_method=gSCS ; exclude_iea=FALSE ; domain_scope=custom ; significant=TRUE ; query=3 80f1f68e42bf388e ; custom_bg=6 b6534d4a26da9fed",
            "gost ; organism=hsapiens ; sources=GO:BP|GO:MF|GO:CC|KEGG|REAC ; user_threshold=0.05 ; correction_method=gSCS ; exclude_iea=FALSE ; domain_scope=custom ; significant=TRUE ; query=1 907fd80ff1d231f2 ; custom_bg=6 b6534d4a26da9fed"),
        bus = c(
            enrichment_analysis_results = "gprofiler_results,clusterprofiler_results,stringdb_results,full_enrichment_results,selected_contrast,analysis_method,parameters e17575402cc40036",
            `enrichment_ui_params$up_log2fc_cutoff` = "0", `enrichment_ui_params$down_log2fc_cutoff` = "0",
            `enrichment_ui_params$q_value_cutoff` = "0.05", `enrichment_ui_params$organism_selected` = "9606",
            `enrichment_ui_params$database_source` = "gprofiler2",
            `enrichment_ui_params$organism_name` = "Homo sapiens",
            `enrichment_ui_params$organism_supported` = "TRUE",
            `enrichment_ui_params$selected_contrast` = "B_vs_A", `enrichment_ui_params$timestamp` = "<time>",
            `enrichment_organism_filter$enabled` = "FALSE",
            `enrichment_organism_filter$filter_applied` = "FALSE",
            `enrichment_organism_filter$target_taxon_id` = "9606",
            `enrichment_organism_filter$proteins_before` = "0",
            `enrichment_organism_filter$proteins_after` = "0",
            `enrichment_organism_filter$proteins_removed` = "0",
            `enrichment_organism_filter$timestamp` = "<time>", uniprot_dat_cln = "5x9 13df57c377de46dc"),
        display = c(
            `groupB-groupA/gprofiler` = "10x23 7bb8700062d08b4b",
            `groupB-groupA/gprofiler_table` = "<saved> 10 e4727adc23d4af92",
            `groupB-groupA/clusterprofiler` = "<NULL>", `groupB-groupA/stringdb` = "<NULL>"),
        files = c(
            `<store entries>..workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-enrichment-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "22",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "10",
            `<store entries>..workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `publication_graphs/Interactive_Volcano_Plots/B_vs_A.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "8a7cae2958352a9afa06dadd17250eacf2c56a50583e8b0b1c17146d9991c109",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.png` = "<present>",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.tsv` = "fbe4b6a921c2abe20abdb68f3fae2a34f3a45e70062415d3da8acebed58e73c6",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.xlsx` = "<present>",
            `results/pathway_enrichment/groupB-groupA_down_enrichment_plot.html` = "<present>",
            `results/pathway_enrichment/groupB-groupA_down_enrichment_plot.png` = "<present>",
            `results/pathway_enrichment/groupB-groupA_down_enrichment_results.tsv` = "51e2bf4b144de5e7f21bf42843fa1db3b9918c54b8ef8bda7e90613838efd31a",
            `results/pathway_enrichment/groupB-groupA_up_enrichment_plot.html` = "<present>",
            `results/pathway_enrichment/groupB-groupA_up_enrichment_plot.png` = "<present>",
            `results/pathway_enrichment/groupB-groupA_up_enrichment_results.tsv` = "51e2bf4b144de5e7f21bf42843fa1db3b9918c54b8ef8bda7e90613838efd31a"),
        returned = c(
            all_enrichment_results = "<list> 0cc1222b458d7a6c", analysis_complete = "TRUE",
            analysis_method = "gprofiler2", clusterprofiler_results = "<NULL>",
            clusterprofiler_results_table = "<NULL>", contrasts_available = "B_vs_A",
            current_contrast_results = "<list> 28df520b992b4a4b", current_s4_object = "<ref>",
            da_results_data = "<list> 00c45cf144cf6614", enrichment_plots = "<list> 28df520b992b4a4b",
            enrichment_product = "<ref>", enrichment_results = "<NULL>",
            enrichment_results_full = "<EnrichmentResults> dbde2059c38a85f0",
            enrichment_service = "unavailable", gprofiler_results = "<data.frame> 7bb8700062d08b4b",
            gprofiler_results_table = "<ref>", organism_supported = "TRUE", stringdb_results = "<NULL>"),
        errors = character(0)))

omicEnrichW1Frozen[["kernel_clusterprofiler"]] <- list(
    da = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<req>",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<req>",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = c(
            "[default NULL da_working] Running differential abundance analysis...",
            "[message 5 -] Differential abundance analysis completed successfully!"),
        updates = c(
            "differential_expression-formula_string: <same> -> ~ 0 + group",
            "differential_expression-volcano_contrast: groupB-groupA -> <none>",
            "differential_expression-heatmap_contrast: groupB-groupA -> <none>",
            "differential_expression-table_contrast: groupB-groupA -> <none>",
            "differential_expression-volcano_contrast: B_vs_A -> <none>",
            "differential_expression-heatmap_contrast: B_vs_A -> <none>",
            "differential_expression-table_contrast: B_vs_A -> <none>",
            "enrichment_analysis-selected_contrast: B_vs_A -> <none>"),
        progress = character(0),
        tab_status = c(
            normalization = "complete", differential_expression = "complete", enrichment_analysis = "pending",
            differential_abundance = "complete"),
        history = c("initial", "correlation_filtered"),
        states = c(initial = "71557d1c8bac9bbe", correlation_filtered = "5acf333d81fb9c8a"),
        product = "<none>",
        resources = character(0),
        kernel = character(0),
        services = character(0),
        bus = c(
            enrichment_analysis_results = "<NULL>", `enrichment_ui_params$` = "<NULL>",
            `enrichment_organism_filter$` = "<NULL>", uniprot_dat_cln = "<NULL>"),
        display = "<empty>",
        files = c(
            `<store entries>..workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "15",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "8",
            `<store entries>..workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `publication_graphs/Interactive_Volcano_Plots/B_vs_A.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "8a7cae2958352a9afa06dadd17250eacf2c56a50583e8b0b1c17146d9991c109",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.png` = "<present>",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.tsv` = "fbe4b6a921c2abe20abdb68f3fae2a34f3a45e70062415d3da8acebed58e73c6",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.xlsx` = "<present>"),
        returned = c(
            all_enrichment_results = "<list> 28df520b992b4a4b", analysis_complete = "FALSE",
            analysis_method = "gprofiler2", clusterprofiler_results = "<NULL>",
            contrasts_available = "B_vs_A", current_contrast_results = "<list> 28df520b992b4a4b",
            current_s4_object = "<ref>", da_results_data = "<list> 00c45cf144cf6614",
            enrichment_plots = "<list> 28df520b992b4a4b", enrichment_results = "<NULL>",
            gprofiler_results = "<NULL>", organism_supported = "TRUE", stringdb_results = "<NULL>"),
        errors = character(0)),
    open = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<req>",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<req>",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = "enrichment_analysis-selected_contrast: B_vs_A -> <none>",
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)),
    select = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[WARNING] CUSTOM ORGANISM\n Method: clusterprofiler\n Organism: Taxon ID 999999\n Using UniProt GO annotations",
            clusterprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = c(
            all_enrichment_results = "<list> 28df520b992b4a4b", analysis_complete = "FALSE",
            analysis_method = "clusterprofiler", clusterprofiler_results = "<NULL>",
            contrasts_available = "B_vs_A", current_contrast_results = "<list> 28df520b992b4a4b",
            current_s4_object = "<ref>", da_results_data = "<list> 00c45cf144cf6614",
            enrichment_plots = "<list> 28df520b992b4a4b", enrichment_results = "<NULL>",
            gprofiler_results = "<NULL>", organism_supported = "FALSE", stringdb_results = "<NULL>"),
        errors = character(0)),
    run = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[WARNING] CUSTOM ORGANISM\n Method: clusterprofiler\n Organism: Taxon ID 999999\n Using UniProt GO annotations",
            clusterprofiler_plot = "<text> 2f283053fa50943b4cc6916efa22c2fc73b6ba470a3db5efdb9d7e0cc49124e2",
            clusterprofiler_results_table = "<text> 4cde847cac7e9c4ab85a093c84e945c7ddf8dfa3f4810877838ac40faaca1db2",
            clusterprofiler_summary_stats = "Total GO terms: 6\nUp-regulated: 3\nDown-regulated: 3\n\nResults displayed in table below.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[OK] Analysis Complete\n\nMethod: clusterprofiler\n\nContrast: B_vs_A\n\nUp log2FC cutoff: 0.0\n\nDown log2FC cutoff: 0.0\n\nQ-value cutoff: 0.050\n\nOrganism: Taxon ID 999999\n\n\nResults Available:\n* gprofiler2: 0 terms\n* clusterProfileR: 6 terms\n* STRING-DB: 0 networks\n\n[OK] Results saved to workflow state",
            gprofiler_plot = "<text> 2f283053fa50943b4cc6916efa22c2fc73b6ba470a3db5efdb9d7e0cc49124e2",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = c(
            "[default NULL enrichment_working] Running enrichment analysis...",
            "[message 5 -] Enrichment analysis completed successfully for contrast: B_vs_A",
            "remove enrichment_working"),
        updates = character(0),
        progress = c(
            "start Running enrichment analysis... 0", "0.2 Transforming DE data...",
            "0.3 Creating DE results S4 object...", "0.5 Running enrichment analysis...",
            "0.8 Storing results...", "1 Complete!"),
        tab_status = c(
            normalization = "complete", differential_expression = "complete",
            enrichment_analysis = "complete", differential_abundance = "complete"),
        history = c("initial", "correlation_filtered", "enrichment_completed"),
        states = c(
            initial = "71557d1c8bac9bbe", correlation_filtered = "5acf333d81fb9c8a",
            enrichment_completed = "9421e4861d18927e"),
        product = c(
            schema = "dia_enrichment_product_v1", service = "unavailable", `parameters$taxon_id` = "999999",
            `parameters$up_cutoff` = "0", `parameters$down_cutoff` = "0", `parameters$q_cutoff` = "0.05",
            `parameters$exclude_iea` = "FALSE", `parameters$correction_method` = "gSCS",
            `parameters$protein_id_column` = "uniprot_acc", `parameters$contrast_names` = "groupB-groupA",
            `results$groupB-groupA` = "<manifest da_proteins_long>", contrasts = "groupB-groupA",
            summaries = "28df520b992b4a4b"),
        resources = "<unchanged>",
        kernel = "kernel=processEnrichments ; taxon_id=999999 ; up_cutoff=0 ; down_cutoff=0 ; q_cutoff=0.05 ; exclude_iea=FALSE ; protein_id_column=uniprot_acc ; contrast_names=groupB-groupA ; correction_method=gSCS ; pathway_dir=<ROOT>/results/pathway_enrichment ; da_data=0656c115da47d5db ; da_names=groupB-groupA ; da_rows=6x17 ; da_contrasts=950ddbb7d830df8c ; da_design=bcd78aff6b352869 ; go_annotations=0218200d5d349843",
        services = c(
            "enricher ; gene=3 779b8ace99d9ed37 ; universe=6 11a1f4f7bc6cf8ce ; TERM2GENE=bdb931bc333e84ee ; pvalueCutoff=0.05 ; pAdjustMethod=BH",
            "enricher ; gene=1 b3c43de36f5ed92f ; universe=6 11a1f4f7bc6cf8ce ; TERM2GENE=bdb931bc333e84ee ; pvalueCutoff=0.05 ; pAdjustMethod=BH"),
        bus = c(
            enrichment_analysis_results = "gprofiler_results,clusterprofiler_results,stringdb_results,full_enrichment_results,selected_contrast,analysis_method,parameters 882acdfd257a0a7e",
            `enrichment_ui_params$up_log2fc_cutoff` = "0", `enrichment_ui_params$down_log2fc_cutoff` = "0",
            `enrichment_ui_params$q_value_cutoff` = "0.05",
            `enrichment_ui_params$organism_selected` = "999999",
            `enrichment_ui_params$database_source` = "clusterprofiler",
            `enrichment_ui_params$organism_name` = "Taxon ID 999999",
            `enrichment_ui_params$organism_supported` = "FALSE",
            `enrichment_ui_params$selected_contrast` = "B_vs_A", `enrichment_ui_params$timestamp` = "<time>",
            `enrichment_organism_filter$enabled` = "FALSE",
            `enrichment_organism_filter$filter_applied` = "FALSE",
            `enrichment_organism_filter$target_taxon_id` = "999999",
            `enrichment_organism_filter$proteins_before` = "0",
            `enrichment_organism_filter$proteins_after` = "0",
            `enrichment_organism_filter$proteins_removed` = "0",
            `enrichment_organism_filter$timestamp` = "<time>", uniprot_dat_cln = "5x9 0218200d5d349843"),
        display = c(
            `groupB-groupA/gprofiler` = "<NULL>", `groupB-groupA/clusterprofiler` = "6x11 3fbdfddf484f8372",
            `groupB-groupA/clusterprofiler_table` = "<saved> 6 3fbdfddf484f8372",
            `groupB-groupA/stringdb` = "<NULL>"),
        files = c(
            `<store entries>..workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-enrichment-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "21",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "10",
            `<store entries>..workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `publication_graphs/Interactive_Volcano_Plots/B_vs_A.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "8a7cae2958352a9afa06dadd17250eacf2c56a50583e8b0b1c17146d9991c109",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.png` = "<present>",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.tsv` = "fbe4b6a921c2abe20abdb68f3fae2a34f3a45e70062415d3da8acebed58e73c6",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.xlsx` = "<present>",
            `results/pathway_enrichment/groupB-groupA_down_enrichment_plot.html` = "<present>",
            `results/pathway_enrichment/groupB-groupA_down_enrichment_plot.png` = "<present>",
            `results/pathway_enrichment/groupB-groupA_down_enrichment_results.tsv` = "c6983b1b6899f01b7555fad52d39a290da33e578504e83ce5b5acf1155b8d14a",
            `results/pathway_enrichment/groupB-groupA_up_enrichment_plot.html` = "<present>",
            `results/pathway_enrichment/groupB-groupA_up_enrichment_plot.png` = "<present>",
            `results/pathway_enrichment/groupB-groupA_up_enrichment_results.tsv` = "c6983b1b6899f01b7555fad52d39a290da33e578504e83ce5b5acf1155b8d14a"),
        returned = c(
            all_enrichment_results = "<list> 6961aaa3fc9d92d9", analysis_complete = "TRUE",
            analysis_method = "clusterprofiler", clusterprofiler_results = "<data.frame> 3fbdfddf484f8372",
            clusterprofiler_results_table = "<ref>", contrasts_available = "B_vs_A",
            current_contrast_results = "<list> 28df520b992b4a4b", current_s4_object = "<ref>",
            da_results_data = "<list> 00c45cf144cf6614", enrichment_plots = "<list> 28df520b992b4a4b",
            enrichment_product = "<ref>", enrichment_results = "<NULL>",
            enrichment_results_full = "<EnrichmentResults> 9421e4861d18927e",
            enrichment_service = "unavailable", gprofiler_results = "<NULL>",
            gprofiler_results_table = "<NULL>", organism_supported = "FALSE", stringdb_results = "<NULL>"),
        errors = character(0)))

omicEnrichW1Frozen[["organism_filter"]] <- list(
    `no_mapping/da` = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<req>",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<req>",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = c(
            "[default NULL da_working] Running differential abundance analysis...",
            "[message 5 -] Differential abundance analysis completed successfully!"),
        updates = c(
            "differential_expression-formula_string: <same> -> ~ 0 + group",
            "differential_expression-volcano_contrast: groupB-groupA -> <none>",
            "differential_expression-heatmap_contrast: groupB-groupA -> <none>",
            "differential_expression-table_contrast: groupB-groupA -> <none>",
            "differential_expression-volcano_contrast: B_vs_A -> <none>",
            "differential_expression-heatmap_contrast: B_vs_A -> <none>",
            "differential_expression-table_contrast: B_vs_A -> <none>",
            "enrichment_analysis-selected_contrast: B_vs_A -> <none>"),
        progress = character(0),
        tab_status = c(
            normalization = "complete", differential_expression = "complete", enrichment_analysis = "pending",
            differential_abundance = "complete"),
        history = c("initial", "correlation_filtered"),
        states = c(initial = "71557d1c8bac9bbe", correlation_filtered = "5acf333d81fb9c8a"),
        product = "<none>",
        resources = character(0),
        kernel = character(0),
        services = character(0),
        bus = c(
            enrichment_analysis_results = "<NULL>", `enrichment_ui_params$` = "<NULL>",
            `enrichment_organism_filter$` = "<NULL>", uniprot_dat_cln = "<NULL>"),
        display = "<empty>",
        files = c(
            `<store entries>..workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "15",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "8",
            `<store entries>..workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `publication_graphs/Interactive_Volcano_Plots/B_vs_A.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "8a7cae2958352a9afa06dadd17250eacf2c56a50583e8b0b1c17146d9991c109",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.png` = "<present>",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.tsv` = "fbe4b6a921c2abe20abdb68f3fae2a34f3a45e70062415d3da8acebed58e73c6",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.xlsx` = "<present>"),
        returned = c(
            all_enrichment_results = "<list> 28df520b992b4a4b", analysis_complete = "FALSE",
            analysis_method = "gprofiler2", clusterprofiler_results = "<NULL>",
            contrasts_available = "B_vs_A", current_contrast_results = "<list> 28df520b992b4a4b",
            current_s4_object = "<ref>", da_results_data = "<list> 00c45cf144cf6614",
            enrichment_plots = "<list> 28df520b992b4a4b", enrichment_results = "<NULL>",
            gprofiler_results = "<NULL>", organism_supported = "TRUE", stringdb_results = "<NULL>"),
        errors = character(0)),
    `no_mapping/open` = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<req>",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<req>",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = "enrichment_analysis-selected_contrast: B_vs_A -> <none>",
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)),
    `no_mapping/select` = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)),
    `no_mapping/run` = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<text> d9625c0c8f0d2e2d05bf97481af72158e56fcbcc0cb47a9ac59b41f286e9dd15",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[OK] Analysis Complete\n\nMethod: gprofiler2\n\nContrast: B_vs_A\n\nUp log2FC cutoff: 0.0\n\nDown log2FC cutoff: 0.0\n\nQ-value cutoff: 0.050\n\nOrganism: Homo sapiens\n\n\nResults Available:\n* gprofiler2: 10 terms\n* clusterProfileR: 0 terms\n* STRING-DB: 0 networks\n\n[OK] Results saved to workflow state",
            gprofiler_plot = "<text> d9625c0c8f0d2e2d05bf97481af72158e56fcbcc0cb47a9ac59b41f286e9dd15",
            gprofiler_results_table = "<text> 86bdcd67360fa73f5cb9cd2636c1a97294f8c85e4065ab1d82da77e3a81df11a",
            gprofiler_summary_stats = "Total enrichment terms: 10\nUp-regulated pathways: 5\nDown-regulated pathways: 5\n\nResults displayed in table below.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = c(
            "[default NULL enrichment_working] Running enrichment analysis...",
            "[warning 8 -] Multi-species filtering enabled but no organism mapping available. Proceeding without filtering.",
            "[message 5 -] Enrichment analysis completed successfully for contrast: B_vs_A",
            "remove enrichment_working"),
        updates = character(0),
        progress = c(
            "start Running enrichment analysis... 0", "0.2 Transforming DE data...",
            "0.3 Creating DE results S4 object...", "0.5 Running enrichment analysis...",
            "0.8 Storing results...", "1 Complete!"),
        tab_status = c(
            normalization = "complete", differential_expression = "complete",
            enrichment_analysis = "complete", differential_abundance = "complete"),
        history = c("initial", "correlation_filtered", "enrichment_completed"),
        states = c(
            initial = "71557d1c8bac9bbe", correlation_filtered = "5acf333d81fb9c8a",
            enrichment_completed = "61dae2e69108addd"),
        product = c(
            schema = "dia_enrichment_product_v1", service = "gprofiler2:gprofiler2_test_double",
            `parameters$taxon_id` = "9606", `parameters$up_cutoff` = "0", `parameters$down_cutoff` = "0",
            `parameters$q_cutoff` = "0.05", `parameters$exclude_iea` = "FALSE",
            `parameters$correction_method` = "gSCS", `parameters$protein_id_column` = "gene_name",
            `parameters$contrast_names` = "groupB-groupA",
            `results$groupB-groupA` = "<manifest da_proteins_long>", contrasts = "groupB-groupA",
            summaries = "6122a4961d0ad898"),
        resources = "gprofiler2 ; gprofiler2_test_double ; database_version ; enrichment of 1 contrasts, taxon 9606",
        kernel = "kernel=runProtEnrichDeterministicProcessEnrichments ; taxon_id=9606 ; up_cutoff=0 ; down_cutoff=0 ; q_cutoff=0.05 ; exclude_iea=FALSE ; protein_id_column=gene_name ; contrast_names=groupB-groupA ; correction_method=gSCS ; pathway_dir=<ROOT>/results/pathway_enrichment ; da_data=0656c115da47d5db ; da_names=groupB-groupA ; da_rows=6x17 ; da_contrasts=950ddbb7d830df8c ; da_design=bcd78aff6b352869 ; go_annotations=13df57c377de46dc",
        services = character(0),
        bus = c(
            enrichment_analysis_results = "gprofiler_results,clusterprofiler_results,stringdb_results,full_enrichment_results,selected_contrast,analysis_method,parameters bd71e1846eddbd60",
            `enrichment_ui_params$up_log2fc_cutoff` = "0", `enrichment_ui_params$down_log2fc_cutoff` = "0",
            `enrichment_ui_params$q_value_cutoff` = "0.05", `enrichment_ui_params$organism_selected` = "9606",
            `enrichment_ui_params$database_source` = "gprofiler2",
            `enrichment_ui_params$organism_name` = "Homo sapiens",
            `enrichment_ui_params$organism_supported` = "TRUE",
            `enrichment_ui_params$selected_contrast` = "B_vs_A", `enrichment_ui_params$timestamp` = "<time>",
            `enrichment_organism_filter$enabled` = "TRUE",
            `enrichment_organism_filter$filter_applied` = "FALSE",
            `enrichment_organism_filter$target_taxon_id` = "9606",
            `enrichment_organism_filter$proteins_before` = "0",
            `enrichment_organism_filter$proteins_after` = "0",
            `enrichment_organism_filter$proteins_removed` = "0",
            `enrichment_organism_filter$timestamp` = "<time>", uniprot_dat_cln = "5x9 13df57c377de46dc"),
        display = c(
            `groupB-groupA/gprofiler` = "10x24 81a566f72a3062f4",
            `groupB-groupA/gprofiler_table` = "<saved> 10 b2a4adf014f1d53e",
            `groupB-groupA/clusterprofiler` = "<NULL>", `groupB-groupA/stringdb` = "<NULL>"),
        files = c(
            `<store entries>..workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-enrichment-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "23",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "10",
            `<store entries>..workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `publication_graphs/Interactive_Volcano_Plots/B_vs_A.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "8a7cae2958352a9afa06dadd17250eacf2c56a50583e8b0b1c17146d9991c109",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.png` = "<present>",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.tsv` = "fbe4b6a921c2abe20abdb68f3fae2a34f3a45e70062415d3da8acebed58e73c6",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.xlsx` = "<present>",
            `results/pathway_enrichment/groupB-groupA_down_enrichment_results.tsv` = "511024ccd86710326a8c45231c97c82595eabf1732e3bab43fd98c53f3b289d7",
            `results/pathway_enrichment/groupB-groupA_up_enrichment_results.tsv` = "b9b76ea1cdcd7200ad5ca6ab181e831781c08df2df93f223d26cda697738397a"),
        returned = c(
            all_enrichment_results = "<list> 123e2a6cd1f92b90", analysis_complete = "TRUE",
            analysis_method = "gprofiler2", clusterprofiler_results = "<NULL>",
            clusterprofiler_results_table = "<NULL>", contrasts_available = "B_vs_A",
            current_contrast_results = "<list> 28df520b992b4a4b", current_s4_object = "<ref>",
            da_results_data = "<list> 00c45cf144cf6614", enrichment_plots = "<list> 28df520b992b4a4b",
            enrichment_product = "<ref>", enrichment_results = "<NULL>",
            enrichment_results_full = "<EnrichmentResults> 61dae2e69108addd",
            enrichment_service = "gprofiler2:gprofiler2_test_double",
            gprofiler_results = "<data.frame> 81a566f72a3062f4", gprofiler_results_table = "<ref>",
            organism_supported = "TRUE", stringdb_results = "<NULL>"),
        errors = character(0)),
    `mapped/da` = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<req>",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<req>",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = c(
            "[default NULL da_working] Running differential abundance analysis...",
            "[message 5 -] Differential abundance analysis completed successfully!"),
        updates = c(
            "differential_expression-formula_string: <same> -> ~ 0 + group",
            "differential_expression-volcano_contrast: groupB-groupA -> <none>",
            "differential_expression-heatmap_contrast: groupB-groupA -> <none>",
            "differential_expression-table_contrast: groupB-groupA -> <none>",
            "differential_expression-volcano_contrast: B_vs_A -> <none>",
            "differential_expression-heatmap_contrast: B_vs_A -> <none>",
            "differential_expression-table_contrast: B_vs_A -> <none>",
            "enrichment_analysis-selected_contrast: B_vs_A -> <none>"),
        progress = character(0),
        tab_status = c(
            normalization = "complete", differential_expression = "complete", enrichment_analysis = "pending",
            differential_abundance = "complete"),
        history = c("initial", "correlation_filtered"),
        states = c(initial = "71557d1c8bac9bbe", correlation_filtered = "5acf333d81fb9c8a"),
        product = "<none>",
        resources = character(0),
        kernel = character(0),
        services = character(0),
        bus = c(
            enrichment_analysis_results = "<NULL>", `enrichment_ui_params$` = "<NULL>",
            `enrichment_organism_filter$` = "<NULL>", uniprot_dat_cln = "<NULL>"),
        display = "<empty>",
        files = c(
            `<store entries>..workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "15",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "8",
            `<store entries>..workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `publication_graphs/Interactive_Volcano_Plots/B_vs_A.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "8a7cae2958352a9afa06dadd17250eacf2c56a50583e8b0b1c17146d9991c109",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.png` = "<present>",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.tsv` = "fbe4b6a921c2abe20abdb68f3fae2a34f3a45e70062415d3da8acebed58e73c6",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.xlsx` = "<present>"),
        returned = c(
            all_enrichment_results = "<list> 28df520b992b4a4b", analysis_complete = "FALSE",
            analysis_method = "gprofiler2", clusterprofiler_results = "<NULL>",
            contrasts_available = "B_vs_A", current_contrast_results = "<list> 28df520b992b4a4b",
            current_s4_object = "<ref>", da_results_data = "<list> 00c45cf144cf6614",
            enrichment_plots = "<list> 28df520b992b4a4b", enrichment_results = "<NULL>",
            gprofiler_results = "<NULL>", organism_supported = "TRUE", stringdb_results = "<NULL>"),
        errors = character(0)),
    `mapped/open` = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<req>",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<req>",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = "enrichment_analysis-selected_contrast: B_vs_A -> <none>",
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)),
    `mapped/mixed` = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<req>",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<req>",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = "[message 5 -] Multi-species FASTA detected. Filtering to Homo sapiens enabled.",
        updates = "enrichment_analysis-enable_organism_filter: <same> -> TRUE",
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)),
    `mapped/select` = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)),
    `mapped/run` = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = c(
            "[default NULL enrichment_working] Running enrichment analysis...",
            "[error 10 -] Error in enrichment analysis: object 'clean_acc' not found",
            "remove enrichment_working"),
        updates = character(0),
        progress = c(
            "start Running enrichment analysis... 0", "0.2 Transforming DE data...",
            "0.3 Creating DE results S4 object..."),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = c(
            enrichment_analysis_results = "<NULL>", `enrichment_ui_params$` = "<NULL>",
            `enrichment_organism_filter$` = "<NULL>", uniprot_dat_cln = "5x9 13df57c377de46dc"),
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)))

omicEnrichW1Frozen[["refusals"]] <- list(
    `pre_da/tab_before_da` = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<req>",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "No contrasts available.\nComplete differential expression\nanalysis first.",
            download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<req>",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = "[warning 5 -] Please complete the differential expression analysis before accessing enrichment analysis.",
        updates = character(0),
        progress = character(0),
        tab_status = c(
            normalization = "complete", differential_expression = "pending", enrichment_analysis = "disabled"),
        history = c("initial", "correlation_filtered"),
        states = c(initial = "71557d1c8bac9bbe", correlation_filtered = "bf8ffd456d2db4fb"),
        product = "<none>",
        resources = character(0),
        kernel = character(0),
        services = character(0),
        bus = c(
            enrichment_analysis_results = "<NULL>", `enrichment_ui_params$` = "<NULL>",
            `enrichment_organism_filter$` = "<NULL>", uniprot_dat_cln = "<NULL>"),
        display = "<empty>",
        files = c(
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "4",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "3"),
        returned = c(
            all_enrichment_results = "<list> 28df520b992b4a4b", analysis_complete = "FALSE",
            analysis_method = "gprofiler2", clusterprofiler_results = "<NULL>",
            contrasts_available = "<NULL>", current_contrast_results = "<list> 28df520b992b4a4b",
            current_s4_object = "<NULL>", da_results_data = "<NULL>",
            enrichment_plots = "<list> 28df520b992b4a4b", enrichment_results = "<NULL>",
            gprofiler_results = "<NULL>", organism_supported = "TRUE", stringdb_results = "<NULL>"),
        errors = character(0)),
    `pre_da/run_before_da` = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "No contrasts available.\nComplete differential expression\nanalysis first.",
            download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)),
    `pre_da/da` = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = c(
            "[default NULL da_working] Running differential abundance analysis...",
            "[message 5 -] Differential abundance analysis completed successfully!"),
        updates = c(
            "differential_expression-formula_string: <same> -> ~ 0 + group",
            "differential_expression-volcano_contrast: groupB-groupA -> <none>",
            "differential_expression-heatmap_contrast: groupB-groupA -> <none>",
            "differential_expression-table_contrast: groupB-groupA -> <none>",
            "differential_expression-volcano_contrast: B_vs_A -> <none>",
            "differential_expression-heatmap_contrast: B_vs_A -> <none>",
            "differential_expression-table_contrast: B_vs_A -> <none>",
            "enrichment_analysis-selected_contrast: B_vs_A -> <none>"),
        progress = character(0),
        tab_status = c(
            normalization = "complete", differential_expression = "complete", enrichment_analysis = "pending",
            differential_abundance = "complete"),
        history = "<unchanged>",
        states = c(initial = "71557d1c8bac9bbe", correlation_filtered = "5acf333d81fb9c8a"),
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = c(
            `<store entries>..workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "15",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "8",
            `<store entries>..workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `publication_graphs/Interactive_Volcano_Plots/B_vs_A.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "8a7cae2958352a9afa06dadd17250eacf2c56a50583e8b0b1c17146d9991c109",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.png` = "<present>",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.tsv` = "fbe4b6a921c2abe20abdb68f3fae2a34f3a45e70062415d3da8acebed58e73c6",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.xlsx` = "<present>"),
        returned = c(
            all_enrichment_results = "<list> 28df520b992b4a4b", analysis_complete = "FALSE",
            analysis_method = "gprofiler2", clusterprofiler_results = "<NULL>",
            contrasts_available = "B_vs_A", current_contrast_results = "<list> 28df520b992b4a4b",
            current_s4_object = "<ref>", da_results_data = "<list> 00c45cf144cf6614",
            enrichment_plots = "<list> 28df520b992b4a4b", enrichment_results = "<NULL>",
            gprofiler_results = "<NULL>", organism_supported = "TRUE", stringdb_results = "<NULL>"),
        errors = character(0)),
    `pre_da/open` = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = "enrichment_analysis-selected_contrast: B_vs_A -> <none>",
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)),
    `pre_da/run_no_contrast` = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<req>",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<req>",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)),
    `pre_da/run_taxon_abc` = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "<error> organism_taxid must be a single positive integer NCBI taxonomy ID.",
            clusterprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = c(
            "[default NULL enrichment_working] Running enrichment analysis...",
            "[error 10 -] Error in enrichment analysis: organism_taxid must be a single positive integer NCBI taxonomy ID.",
            "remove enrichment_working"),
        updates = character(0),
        progress = c(
            "start Running enrichment analysis... 0", "0.2 Transforming DE data...",
            "0.3 Creating DE results S4 object...", "0.5 Running enrichment analysis..."),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = c(
            enrichment_analysis_results = "<NULL>", `enrichment_ui_params$` = "<NULL>",
            `enrichment_organism_filter$enabled` = "FALSE",
            `enrichment_organism_filter$filter_applied` = "FALSE",
            `enrichment_organism_filter$target_taxon_id` = "abc",
            `enrichment_organism_filter$proteins_before` = "0",
            `enrichment_organism_filter$proteins_after` = "0",
            `enrichment_organism_filter$proteins_removed` = "0",
            `enrichment_organism_filter$timestamp` = "<time>", uniprot_dat_cln = "5x9 0218200d5d349843"),
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = "<unhandled> organism_taxid must be a single positive integer NCBI taxonomy ID."),
    `q_cutoff/da` = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<req>",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<req>",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = c(
            "[default NULL da_working] Running differential abundance analysis...",
            "[message 5 -] Differential abundance analysis completed successfully!"),
        updates = c(
            "differential_expression-formula_string: <same> -> ~ 0 + group",
            "differential_expression-volcano_contrast: groupB-groupA -> <none>",
            "differential_expression-heatmap_contrast: groupB-groupA -> <none>",
            "differential_expression-table_contrast: groupB-groupA -> <none>",
            "differential_expression-volcano_contrast: B_vs_A -> <none>",
            "differential_expression-heatmap_contrast: B_vs_A -> <none>",
            "differential_expression-table_contrast: B_vs_A -> <none>",
            "enrichment_analysis-selected_contrast: B_vs_A -> <none>"),
        progress = character(0),
        tab_status = c(
            normalization = "complete", differential_expression = "complete", enrichment_analysis = "pending",
            differential_abundance = "complete"),
        history = c("initial", "correlation_filtered"),
        states = c(initial = "71557d1c8bac9bbe", correlation_filtered = "5acf333d81fb9c8a"),
        product = "<none>",
        resources = character(0),
        kernel = character(0),
        services = character(0),
        bus = c(
            enrichment_analysis_results = "<NULL>", `enrichment_ui_params$` = "<NULL>",
            `enrichment_organism_filter$` = "<NULL>", uniprot_dat_cln = "<NULL>"),
        display = "<empty>",
        files = c(
            `<store entries>..workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "15",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "8",
            `<store entries>..workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `publication_graphs/Interactive_Volcano_Plots/B_vs_A.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "8a7cae2958352a9afa06dadd17250eacf2c56a50583e8b0b1c17146d9991c109",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.png` = "<present>",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.tsv` = "fbe4b6a921c2abe20abdb68f3fae2a34f3a45e70062415d3da8acebed58e73c6",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.xlsx` = "<present>"),
        returned = c(
            all_enrichment_results = "<list> 28df520b992b4a4b", analysis_complete = "FALSE",
            analysis_method = "gprofiler2", clusterprofiler_results = "<NULL>",
            contrasts_available = "B_vs_A", current_contrast_results = "<list> 28df520b992b4a4b",
            current_s4_object = "<ref>", da_results_data = "<list> 00c45cf144cf6614",
            enrichment_plots = "<list> 28df520b992b4a4b", enrichment_results = "<NULL>",
            gprofiler_results = "<NULL>", organism_supported = "TRUE", stringdb_results = "<NULL>"),
        errors = character(0)),
    `q_cutoff/open` = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<req>",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<req>",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = "enrichment_analysis-selected_contrast: B_vs_A -> <none>",
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)),
    `q_cutoff/run_q_cutoff` = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = c(
            "[default NULL enrichment_working] Running enrichment analysis...",
            "[error 10 -] Error in enrichment analysis: q_cutoff must be between 0 and 1.",
            "remove enrichment_working"),
        updates = character(0),
        progress = c(
            "start Running enrichment analysis... 0", "0.2 Transforming DE data...",
            "0.3 Creating DE results S4 object...", "0.5 Running enrichment analysis..."),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = c(
            enrichment_analysis_results = "<NULL>", `enrichment_ui_params$` = "<NULL>",
            `enrichment_organism_filter$enabled` = "FALSE",
            `enrichment_organism_filter$filter_applied` = "FALSE",
            `enrichment_organism_filter$target_taxon_id` = "9606",
            `enrichment_organism_filter$proteins_before` = "0",
            `enrichment_organism_filter$proteins_after` = "0",
            `enrichment_organism_filter$proteins_removed` = "0",
            `enrichment_organism_filter$timestamp` = "<time>", uniprot_dat_cln = "5x9 13df57c377de46dc"),
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)),
    `table_removed/da` = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<req>",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<req>",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = c(
            "[default NULL da_working] Running differential abundance analysis...",
            "[message 5 -] Differential abundance analysis completed successfully!"),
        updates = c(
            "differential_expression-formula_string: <same> -> ~ 0 + group",
            "differential_expression-volcano_contrast: groupB-groupA -> <none>",
            "differential_expression-heatmap_contrast: groupB-groupA -> <none>",
            "differential_expression-table_contrast: groupB-groupA -> <none>",
            "differential_expression-volcano_contrast: B_vs_A -> <none>",
            "differential_expression-heatmap_contrast: B_vs_A -> <none>",
            "differential_expression-table_contrast: B_vs_A -> <none>",
            "enrichment_analysis-selected_contrast: B_vs_A -> <none>"),
        progress = character(0),
        tab_status = c(
            normalization = "complete", differential_expression = "complete", enrichment_analysis = "pending",
            differential_abundance = "complete"),
        history = c("initial", "correlation_filtered"),
        states = c(initial = "71557d1c8bac9bbe", correlation_filtered = "5acf333d81fb9c8a"),
        product = "<none>",
        resources = character(0),
        kernel = character(0),
        services = character(0),
        bus = c(
            enrichment_analysis_results = "<NULL>", `enrichment_ui_params$` = "<NULL>",
            `enrichment_organism_filter$` = "<NULL>", uniprot_dat_cln = "<NULL>"),
        display = "<empty>",
        files = c(
            `<store entries>..workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "15",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "8",
            `<store entries>..workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `publication_graphs/Interactive_Volcano_Plots/B_vs_A.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "8a7cae2958352a9afa06dadd17250eacf2c56a50583e8b0b1c17146d9991c109",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.png` = "<present>",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.tsv` = "fbe4b6a921c2abe20abdb68f3fae2a34f3a45e70062415d3da8acebed58e73c6",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.xlsx` = "<present>"),
        returned = c(
            all_enrichment_results = "<list> 28df520b992b4a4b", analysis_complete = "FALSE",
            analysis_method = "gprofiler2", clusterprofiler_results = "<NULL>",
            contrasts_available = "B_vs_A", current_contrast_results = "<list> 28df520b992b4a4b",
            current_s4_object = "<ref>", da_results_data = "<list> 00c45cf144cf6614",
            enrichment_plots = "<list> 28df520b992b4a4b", enrichment_results = "<NULL>",
            gprofiler_results = "<NULL>", organism_supported = "TRUE", stringdb_results = "<NULL>"),
        errors = character(0)),
    `table_removed/open` = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<req>",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<req>",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = "enrichment_analysis-selected_contrast: B_vs_A -> <none>",
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)),
    `table_removed/run_table_removed` = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = c(
            "[default NULL enrichment_working] Running enrichment analysis...",
            "[error 10 -] Error in enrichment analysis: The DA results table is missing for contrast groupB-groupA.",
            "remove enrichment_working"),
        updates = character(0),
        progress = c(
            "start Running enrichment analysis... 0", "0.2 Transforming DE data...",
            "0.3 Creating DE results S4 object..."),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = c(
            `<store entries>..workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "15",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "8",
            `<store entries>..workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `publication_graphs/Interactive_Volcano_Plots/B_vs_A.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "8a7cae2958352a9afa06dadd17250eacf2c56a50583e8b0b1c17146d9991c109",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.png` = "<present>",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.xlsx` = "<present>"),
        returned = "<unchanged>",
        errors = character(0)))

omicEnrichW1Frozen[["restore"]] <- list(
    da = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<req>",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<req>",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = c(
            "[default NULL da_working] Running differential abundance analysis...",
            "[message 5 -] Differential abundance analysis completed successfully!"),
        updates = c(
            "differential_expression-formula_string: <same> -> ~ 0 + group",
            "differential_expression-volcano_contrast: groupB-groupA -> <none>",
            "differential_expression-heatmap_contrast: groupB-groupA -> <none>",
            "differential_expression-table_contrast: groupB-groupA -> <none>",
            "differential_expression-volcano_contrast: B_vs_A -> <none>",
            "differential_expression-heatmap_contrast: B_vs_A -> <none>",
            "differential_expression-table_contrast: B_vs_A -> <none>",
            "enrichment_analysis-selected_contrast: B_vs_A -> <none>"),
        progress = character(0),
        tab_status = c(
            normalization = "complete", differential_expression = "complete", enrichment_analysis = "pending",
            differential_abundance = "complete"),
        history = c("initial", "correlation_filtered"),
        states = c(initial = "71557d1c8bac9bbe", correlation_filtered = "5acf333d81fb9c8a"),
        product = "<none>",
        resources = character(0),
        kernel = character(0),
        services = character(0),
        bus = c(
            enrichment_analysis_results = "<NULL>", `enrichment_ui_params$` = "<NULL>",
            `enrichment_organism_filter$` = "<NULL>", uniprot_dat_cln = "<NULL>"),
        display = "<empty>",
        files = c(
            `<store entries>..workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "15",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "8",
            `<store entries>..workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `publication_graphs/Interactive_Volcano_Plots/B_vs_A.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "8a7cae2958352a9afa06dadd17250eacf2c56a50583e8b0b1c17146d9991c109",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.png` = "<present>",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.tsv` = "fbe4b6a921c2abe20abdb68f3fae2a34f3a45e70062415d3da8acebed58e73c6",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.xlsx` = "<present>"),
        returned = c(
            all_enrichment_results = "<list> 28df520b992b4a4b", analysis_complete = "FALSE",
            analysis_method = "gprofiler2", clusterprofiler_results = "<NULL>",
            contrasts_available = "B_vs_A", current_contrast_results = "<list> 28df520b992b4a4b",
            current_s4_object = "<ref>", da_results_data = "<list> 00c45cf144cf6614",
            enrichment_plots = "<list> 28df520b992b4a4b", enrichment_results = "<NULL>",
            gprofiler_results = "<NULL>", organism_supported = "TRUE", stringdb_results = "<NULL>"),
        errors = character(0)),
    open = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<req>",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<req>",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = "enrichment_analysis-selected_contrast: B_vs_A -> <none>",
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)),
    select = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<text> d4fb7d7310ead9d16c488b464913a0e5924d14ef0038fdaef46b4caa5cd76332",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)),
    run = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<text> d9625c0c8f0d2e2d05bf97481af72158e56fcbcc0cb47a9ac59b41f286e9dd15",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[OK] Analysis Complete\n\nMethod: gprofiler2\n\nContrast: B_vs_A\n\nUp log2FC cutoff: 0.0\n\nDown log2FC cutoff: 0.0\n\nQ-value cutoff: 0.050\n\nOrganism: Homo sapiens\n\n\nResults Available:\n* gprofiler2: 10 terms\n* clusterProfileR: 0 terms\n* STRING-DB: 0 networks\n\n[OK] Results saved to workflow state",
            gprofiler_plot = "<text> d9625c0c8f0d2e2d05bf97481af72158e56fcbcc0cb47a9ac59b41f286e9dd15",
            gprofiler_results_table = "<text> 86bdcd67360fa73f5cb9cd2636c1a97294f8c85e4065ab1d82da77e3a81df11a",
            gprofiler_summary_stats = "Total enrichment terms: 10\nUp-regulated pathways: 5\nDown-regulated pathways: 5\n\nResults displayed in table below.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = c(
            "[default NULL enrichment_working] Running enrichment analysis...",
            "[message 5 -] Enrichment analysis completed successfully for contrast: B_vs_A",
            "remove enrichment_working"),
        updates = character(0),
        progress = c(
            "start Running enrichment analysis... 0", "0.2 Transforming DE data...",
            "0.3 Creating DE results S4 object...", "0.5 Running enrichment analysis...",
            "0.8 Storing results...", "1 Complete!"),
        tab_status = c(
            normalization = "complete", differential_expression = "complete",
            enrichment_analysis = "complete", differential_abundance = "complete"),
        history = c("initial", "correlation_filtered", "enrichment_completed"),
        states = c(
            initial = "71557d1c8bac9bbe", correlation_filtered = "5acf333d81fb9c8a",
            enrichment_completed = "61dae2e69108addd"),
        product = c(
            schema = "dia_enrichment_product_v1", service = "gprofiler2:gprofiler2_test_double",
            `parameters$taxon_id` = "9606", `parameters$up_cutoff` = "0", `parameters$down_cutoff` = "0",
            `parameters$q_cutoff` = "0.05", `parameters$exclude_iea` = "FALSE",
            `parameters$correction_method` = "gSCS", `parameters$protein_id_column` = "gene_name",
            `parameters$contrast_names` = "groupB-groupA",
            `results$groupB-groupA` = "<manifest da_proteins_long>", contrasts = "groupB-groupA",
            summaries = "6122a4961d0ad898"),
        resources = "gprofiler2 ; gprofiler2_test_double ; database_version ; enrichment of 1 contrasts, taxon 9606",
        kernel = "kernel=runProtEnrichDeterministicProcessEnrichments ; taxon_id=9606 ; up_cutoff=0 ; down_cutoff=0 ; q_cutoff=0.05 ; exclude_iea=FALSE ; protein_id_column=gene_name ; contrast_names=groupB-groupA ; correction_method=gSCS ; pathway_dir=<ROOT>/results/pathway_enrichment ; da_data=0656c115da47d5db ; da_names=groupB-groupA ; da_rows=6x17 ; da_contrasts=950ddbb7d830df8c ; da_design=bcd78aff6b352869 ; go_annotations=13df57c377de46dc",
        services = character(0),
        bus = c(
            enrichment_analysis_results = "gprofiler_results,clusterprofiler_results,stringdb_results,full_enrichment_results,selected_contrast,analysis_method,parameters bd71e1846eddbd60",
            `enrichment_ui_params$up_log2fc_cutoff` = "0", `enrichment_ui_params$down_log2fc_cutoff` = "0",
            `enrichment_ui_params$q_value_cutoff` = "0.05", `enrichment_ui_params$organism_selected` = "9606",
            `enrichment_ui_params$database_source` = "gprofiler2",
            `enrichment_ui_params$organism_name` = "Homo sapiens",
            `enrichment_ui_params$organism_supported` = "TRUE",
            `enrichment_ui_params$selected_contrast` = "B_vs_A", `enrichment_ui_params$timestamp` = "<time>",
            `enrichment_organism_filter$enabled` = "FALSE",
            `enrichment_organism_filter$filter_applied` = "FALSE",
            `enrichment_organism_filter$target_taxon_id` = "9606",
            `enrichment_organism_filter$proteins_before` = "0",
            `enrichment_organism_filter$proteins_after` = "0",
            `enrichment_organism_filter$proteins_removed` = "0",
            `enrichment_organism_filter$timestamp` = "<time>", uniprot_dat_cln = "5x9 13df57c377de46dc"),
        display = c(
            `groupB-groupA/gprofiler` = "10x24 81a566f72a3062f4",
            `groupB-groupA/gprofiler_table` = "<saved> 10 b2a4adf014f1d53e",
            `groupB-groupA/clusterprofiler` = "<NULL>", `groupB-groupA/stringdb` = "<NULL>"),
        files = c(
            `<store entries>..workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-enrichment-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "23",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "10",
            `<store entries>..workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `publication_graphs/Interactive_Volcano_Plots/B_vs_A.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "8a7cae2958352a9afa06dadd17250eacf2c56a50583e8b0b1c17146d9991c109",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.png` = "<present>",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.tsv` = "fbe4b6a921c2abe20abdb68f3fae2a34f3a45e70062415d3da8acebed58e73c6",
            `results/da_proteins/da_proteins_groupB-groupA_long_annot.xlsx` = "<present>",
            `results/pathway_enrichment/groupB-groupA_down_enrichment_results.tsv` = "511024ccd86710326a8c45231c97c82595eabf1732e3bab43fd98c53f3b289d7",
            `results/pathway_enrichment/groupB-groupA_up_enrichment_results.tsv` = "b9b76ea1cdcd7200ad5ca6ab181e831781c08df2df93f223d26cda697738397a"),
        returned = c(
            all_enrichment_results = "<list> 123e2a6cd1f92b90", analysis_complete = "TRUE",
            analysis_method = "gprofiler2", clusterprofiler_results = "<NULL>",
            clusterprofiler_results_table = "<NULL>", contrasts_available = "B_vs_A",
            current_contrast_results = "<list> 28df520b992b4a4b", current_s4_object = "<ref>",
            da_results_data = "<list> 00c45cf144cf6614", enrichment_plots = "<list> 28df520b992b4a4b",
            enrichment_product = "<ref>", enrichment_results = "<NULL>",
            enrichment_results_full = "<EnrichmentResults> 61dae2e69108addd",
            enrichment_service = "gprofiler2:gprofiler2_test_double",
            gprofiler_results = "<data.frame> 81a566f72a3062f4", gprofiler_results_table = "<ref>",
            organism_supported = "TRUE", stringdb_results = "<NULL>"),
        errors = character(0)),
    restore = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<text> d9625c0c8f0d2e2d05bf97481af72158e56fcbcc0cb47a9ac59b41f286e9dd15",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "B_vs_A", download_enrichment_results = "<download>",
            enrichment_status = "[OK] Analysis Complete\n\nMethod: gprofiler2\n\nContrast: B_vs_A\n\nUp log2FC cutoff: 0.0\n\nDown log2FC cutoff: 0.0\n\nQ-value cutoff: 0.050\n\nOrganism: Homo sapiens\n\n\nResults Available:\n* gprofiler2: 10 terms\n* clusterProfileR: 0 terms\n* STRING-DB: 0 networks\n\n[OK] Results saved to workflow state",
            gprofiler_plot = "<text> d9625c0c8f0d2e2d05bf97481af72158e56fcbcc0cb47a9ac59b41f286e9dd15",
            gprofiler_results_table = "<text> 86bdcd67360fa73f5cb9cd2636c1a97294f8c85e4065ab1d82da77e3a81df11a",
            gprofiler_summary_stats = "Total enrichment terms: 10\nUp-regulated pathways: 5\nDown-regulated pathways: 5\n\nResults displayed in table below.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0),
        restore = c(
            restored = "TRUE", again = "FALSE", analysis_complete = "TRUE",
            enrichment_results_full = "61dae2e69108addd",
            `groupB-groupA/gprofiler` = "10x24 81a566f72a3062f4",
            `groupB-groupA/gprofiler_table` = "<saved> 10 b2a4adf014f1d53e",
            `groupB-groupA/clusterprofiler` = "<NULL>", `groupB-groupA/stringdb` = "<NULL>")))

omicEnrichW1Frozen[["prefill"]] <- list(
    start = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Homo sapiens\n All enrichment methods available",
            clusterprofiler_plot = "<req>",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "No contrasts available.\nComplete differential expression\nanalysis first.",
            download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<req>",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        tab_status = c(
            normalization = "complete", differential_expression = "pending", enrichment_analysis = "disabled"),
        history = c("initial", "correlation_filtered"),
        states = c(initial = "71557d1c8bac9bbe", correlation_filtered = "bf8ffd456d2db4fb"),
        product = "<none>",
        resources = character(0),
        kernel = character(0),
        services = character(0),
        bus = c(
            enrichment_analysis_results = "<NULL>", `enrichment_ui_params$` = "<NULL>",
            `enrichment_organism_filter$` = "<NULL>", uniprot_dat_cln = "<NULL>"),
        display = "<empty>",
        files = c(
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "4",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "3"),
        returned = c(
            all_enrichment_results = "<list> 28df520b992b4a4b", analysis_complete = "FALSE",
            analysis_method = "gprofiler2", clusterprofiler_results = "<NULL>",
            contrasts_available = "<NULL>", current_contrast_results = "<list> 28df520b992b4a4b",
            current_s4_object = "<NULL>", da_results_data = "<NULL>",
            enrichment_plots = "<list> 28df520b992b4a4b", enrichment_results = "<NULL>",
            gprofiler_results = "<NULL>", organism_supported = "TRUE", stringdb_results = "<NULL>"),
        errors = character(0)),
    prefill = list(
        outputs = c(
            "analysis_method_display", "clusterprofiler_plot", "clusterprofiler_results_table",
            "clusterprofiler_summary_stats", "contrasts_display", "download_enrichment_results",
            "enrichment_status", "gprofiler_plot", "gprofiler_results_table", "gprofiler_summary_stats",
            "stringdb_plot", "stringdb_results_table", "stringdb_summary_stats"),
        values = c(
            analysis_method_display = "[OK] SUPPORTED ORGANISM\n Method: gprofiler2\n Species: Mus musculus\n All enrichment methods available",
            clusterprofiler_plot = "<req>",
            clusterprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            clusterprofiler_summary_stats = "No clusterProfileR results available.",
            contrasts_display = "No contrasts available.\nComplete differential expression\nanalysis first.",
            download_enrichment_results = "<download>",
            enrichment_status = "[WAITING] Ready for analysis\n\n\nSteps:\n1. Select contrast from DE results\n2. Set log fold change cutoffs\n3. Set Q-value cutoff (significance threshold)\n4. Click 'Run Enrichment Analysis'\n\nMethod automatically determined by organism.",
            gprofiler_plot = "<req>",
            gprofiler_results_table = "<text> 0d4a07d6c799091f4ffb5599776586a0a562cb1b8ebe148b359df3930319390f",
            gprofiler_summary_stats = "No gprofiler2 results available.",
            stringdb_plot = "<text> b54fa6cb6edbb20c608dfc64ede6de6565c6567b1f2102fccf82f81af2280e51",
            stringdb_results_table = "<text> 6102723c1b0109d3cdf60900ead49920407f04abd755b07022ef33f47454c2ac",
            stringdb_summary_stats = "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."),
        notes = character(0),
        updates = "enrichment_analysis-organism_taxid: <same> -> 10090",
        progress = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        product = "<unchanged>",
        resources = "<unchanged>",
        kernel = character(0),
        services = character(0),
        bus = "<unchanged>",
        display = "<unchanged>",
        files = "<unchanged>",
        returned = "<unchanged>",
        errors = character(0)))
# OMIC_ENRICH_FROZEN_END

expectOmicEnrichW1Case <- function(case) {
    seen <- omicEnrichW1Case(case)
    frozen <- omicEnrichW1Frozen[[case]]
    testthat::expect_false(is.null(frozen), info = case)
    testthat::expect_identical(names(seen), names(frozen), info = case)
    for (step in names(frozen)) {
        for (field in union(names(frozen[[step]]), names(seen[[step]]))) {
            testthat::expect_identical(seen[[step]][[field]], frozen[[step]][[field]], info = paste(case, step, field))
        }
    }
}

for (case in c("gprofiler", "clusterprofiler", "kernel_gprofiler", "kernel_clusterprofiler", "organism_filter", "refusals",
    "restore", "prefill")) {
    local({
        name <- case
        test_that(sprintf("(b) in process: %s", name), expectOmicEnrichW1Case(name))
    })
}

# ------------------------------------------------------------------------------------------------
# KIT-107 W3: the defects fixed on the shared frame, each written before its fix and shown failing on the W3 start commit
# (eacd4a5b; kit-107/test-first.txt). They drive W1 (b)'s fixture and steps.
# ------------------------------------------------------------------------------------------------

#' W1 (b)'s case up to its Run (the DA, the tab, the contrast, Run), with `extra` steps after it.
omicEnrichW3Drive <- function(case, extra, env = parent.frame()) {
    w1 <- omicEnrichW1CaseSpec(case, env)[[1L]]
    omicEnrichW1Drive(w1$setup, c(w1$steps[c("da", "open", "select", "run")], extra), live = w1$live)
}

test_that("(W3, F107-01) a second Run on the same tab gives the first run's product parameters, state digest and display", {
    seen <- omicEnrichW3Drive("gprofiler", list(again = omicEnrichW1Run(2L)))
    expect_identical(seen$again$errors, character())
    expect_identical(seen$again$notes, seen$run$notes)
    expect_identical(seen$again$progress, seen$run$progress)
    # The second run commits `enrichment_completed` again (the store gains state entries; the history names each state
    # once) with the first's digest, and saves the same product and display tables (held as "<unchanged>" from the first
    # run's capture).
    state_entries <- "<store entries>..workflow-entities/entities/dia-state-<ID>"
    expect_gt(as.integer(seen$again$files[[state_entries]]), as.integer(seen$run$files[[state_entries]]))
    expect_identical(seen$again$history, "<unchanged>")
    expect_identical(seen$again$states, "<unchanged>")
    expect_identical(seen$again$product, "<unchanged>")
    expect_identical(seen$again$display, "<unchanged>")
    expect_identical(seen$again$values, seen$run$values)
})

test_that("(W3, F107-02) after an Open on an owner that holds enrichment_completed, selecting a contrast shows the restored tables", {
    reopen <- function(ctx) {
        # An Open: a fresh tab (its initial values) mounted on the bus whose owner holds the completed state.
        session <- shiny::getDefaultReactiveDomain()
        tab <- mod_prot_enrich_server("reopened", ctx$bus(), ctx$paths, "proteomics", "kit107", selected_tab = function() ctx$tab())
        defaults <- omicsUiInputDefaults(mod_prot_enrich_ui("x"))
        names(defaults) <- paste0("reopened-", sub("^x-", "", names(defaults)))
        do.call(session$setInputs, defaults)
        before <- length(shiny::isolate(tab$all_enrichment_results))
        session$setInputs(`reopened-organism_taxid` = "9606", `reopened-selected_contrast` = "B_vs_A")
        manager <- shiny::isolate(ctx$bus()$state_manager)
        shown <- vapply(c("gprofiler_summary_stats", "enrichment_status"), function(name) guiEntryOutputValue(
            tryCatch(session$output[[paste0("reopened-", name)]], error = function(e) paste("<error>", conditionMessage(e))),
            ctx$roots), "")
        ctx$keep("reopened", list(before = before, complete = shiny::isolate(tab$analysis_complete),
            display = omicEnrichW1Display(shiny::isolate(tab$all_enrichment_results), manager), shown = shown))
    }
    seen <- omicEnrichW3Drive("gprofiler", list(reopen = reopen))
    reopened <- seen$reopen$reopened
    expect_identical(reopened$before, 0L)
    expect_true(reopened$complete)
    expect_identical(reopened$display, seen$run$display)
    expect_identical(reopened$shown[["gprofiler_summary_stats"]], seen$run$values[["gprofiler_summary_stats"]])
    expect_identical(reopened$shown[["enrichment_status"]], seen$run$values[["enrichment_status"]])
})

for (case in c("kernel_gprofiler", "kernel_clusterprofiler")) {
    local({
        name <- case
        test_that(sprintf("(W3, F107-05) a Run on the live kernel (%s) commits enrichment_completed with the tab's results", name), {
            check <- function(ctx) {
                bus <- ctx$bus()
                manager <- shiny::isolate(bus$state_manager)
                held <- shiny::isolate(ctx$returned()$enrichment_results_full)
                committed <- "enrichment_completed" %in% unlist(diaStateHistory(manager))
                state <- if (committed) withWorkflowStateScope(manager, function() manager$getState("enrichment_completed"))
                ctx$keep("commit", list(committed = committed, class = class(state)[[1L]],
                    data = identical(state@enrichment_data, held@enrichment_data),
                    summaries = identical(state@enrichment_summaries, held@enrichment_summaries),
                    contrasts = identical(state@contrasts, held@contrasts),
                    plots = identical(lapply(state@enrichment_plots, function(x) lapply(x, class)),
                        lapply(held@enrichment_plots, function(x) lapply(x, class)))))
            }
            seen <- omicEnrichW3Drive(name, list(check = check))
            commit <- seen$check$commit
            expect_true(commit$committed)
            expect_identical(commit$class, "EnrichmentResults")
            expect_true(commit$data)
            expect_true(commit$summaries)
            expect_true(commit$contrasts)
            expect_true(commit$plots)
            expect_identical(seen$run$history, c("initial", "correlation_filtered", "enrichment_completed"))
        })
    })
}
