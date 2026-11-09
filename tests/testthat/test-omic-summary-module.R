# KIT-108 W1: the Summary & Export tab through its public entry points, written before the shared summary frame
# (R/mod_omic_summary.R; W3, W4) and the summary core (R/utils_workflow_summary.R, R/utils_workflow_dia_summary.R) exist,
# and frozen on the start commit (5be34ee5; R/ identical to kit_107's source 8b92d7ac;
# docs/qa/workflow-execution/kit/kit-108/w1/test-first.txt).
#   (a) Identity: the HTML digests of `mod_prot_summary_ui("session_summary")`, `mod_lipid_summary_ui("summary")`,
#       `mod_metab_summary_ui("summary")` and `mod_lipid_summary_ui("summary", profile = trxOmicProfile())` under
#       MULTISCHOLAR_TEST_MODE=true, and their namespaced input and output ids, test ids, labels, button texts,
#       placeholders, text areas, conditions, help texts and defaults (`omicsUiInputDefaults`).
#   (b) The server in process: each exported server under `testServer`, mounted as the shell mounts it
#       (`callOmicGuiComponent`'s `summary` call; trx's with the profile, beside its normalisation and DA tabs), on
#       MCI-011's (proteomics), MCI-017's (metab), MCI-023's (lipid) and MCI-032's (trx) fixtures. Cases, per omic:
#       refusals (save, copy, report and export without a label; the report before the copy; the push before the
#       report), the save (in process and managed, the record written or failing), the copy (complete, an optional
#       artifact missing, required ones missing, the proteomics export gate refusing, the `project_dirs` global absent or
#       present), the report (the template in `scripts/<omic>/`, copied from the package, the GitHub fallback with the
#       download double failing offline, the record bundle failing; proteomics' four routes and a stale requested
#       template), the session-state export and the GitHub push (disabled; enabled, with the push double); trx's RUVIII
#       export (M = group, technical replicates offered and run, the unavailable choice refused, the save before the
#       export, managed).
#   Captured after each step: the tab's registered outputs and each value as the entry suite reads it
#   (`guiEntryOutputValue`; `report_ready` is read as today's `"<error> unused arguments (self, name)"`, D2), the
#   notifications (type, duration, id, message), the progress, the input updates, the reactive values (the
#   `reactiveValues` the tab creates with `workflow_args_saved`), the calls that reach `copyToResultsSummary`,
#   `RenderReport`, `createWorkflowArgsFromConfig`, `omicsManagedProduct`, `writeWorkflowAnalysisRecord`,
#   `exportWorkflowAnalysisRecordBundle` and `pushProjectToGithubFromDirs` (recorders replaced by name; the render and
#   the push are doubles, the others call the real function), the `project_dirs` global before and after, the project's
#   files (table digests; after the first step, the step's changes), `study_parameters.txt` (masked digest), the
#   integration files (name and class or digest),
#   the session-state RDS's fields (times masked) and, for trx, the head's RUVIII export.
# MULTISCHOLAR_KIT108_W1_DUMP=<file> writes the captures as R code (how the frozen values were made; two dumps in
# separate processes were byte-identical before they were pasted).

omicSummaryW1DumpPath <- function() {
    path <- Sys.getenv("MULTISCHOLAR_KIT108_W1_DUMP")
    if (nzchar(path)) path else NULL
}

#' R code for a capture (lists, named vectors and scalars), as test-omic-enrich-module.R's `omicEnrichW1Code`.
omicSummaryW1Code <- function(x, indent = 0L) {
    pad <- strrep(" ", indent + 4L)
    label <- function(names, i) {
        if (is.null(names) || !nzchar(names[[i]])) return("")
        paste0(if (make.names(names[[i]]) == names[[i]]) names[[i]] else paste0("`", names[[i]], "`"), " = ")
    }
    if (is.list(x)) {
        if (length(x) == 0L) return("list()")
        items <- vapply(seq_along(x), function(i) paste0(pad, label(names(x), i), omicSummaryW1Code(x[[i]], indent + 4L)), "")
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

omicSummaryW1Dump <- function(name, value) {
    cat(sprintf("%s <- %s\n\n", name, omicSummaryW1Code(value)), file = omicSummaryW1DumpPath(), append = TRUE)
}

# ------------------------------------------------------------------------------------------------
# (a) Identity: the four summary UIs, frozen at 5be34ee5.
# ------------------------------------------------------------------------------------------------

#' The four summary UIs as the shell mounts them (omicGuiComponentUi; module ids `session_summary` and `summary`).
omicSummaryW1Uis <- function() {
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true")
    list(proteomics = mod_prot_summary_ui("session_summary"), lipidomics = mod_lipid_summary_ui("summary"),
        metabolomics = mod_metab_summary_ui("summary"),
        transcriptomics = mod_lipid_summary_ui("summary", profile = trxOmicProfile()))
}

#' The UI's HTML digest and its parts: inputs (every namespaced `id` that is not an output or a label), outputs,
#' `data-testid` values, each label's text by its input, each button's text by its id, the placeholders, the text areas'
#' values, the conditional panels' conditions, the help texts in order and the inputs' defaults.
omicSummaryW1UiParts <- function(ui) {
    html <- as.character(htmltools::renderTags(ui)$html)
    tags <- regmatches(html, gregexpr("<[a-zA-Z][^>]*>", html))[[1L]]
    attr_of <- function(name) vapply(tags, function(tag) {
        m <- regmatches(tag, regexec(sprintf("[[:space:]]%s=\"([^\"]*)\"", name), tag))[[1L]]
        if (length(m)) m[[2L]] else NA_character_
    }, "", USE.NAMES = FALSE)
    ids <- attr_of("id")
    classes <- attr_of("class")
    is_output <- !is.na(ids) & grepl("(^| )(shiny-[a-z]+-output|html-widget-output)( |$)", classes)
    inputs <- ids[!is.na(ids) & !is_output & !grepl("-label$", ids)]
    labels <- regmatches(html, gregexpr("<label[^>]*for=\"[^\"]+\"[^>]*>[^<]*</label>", html))[[1L]]
    label_text <- stats::setNames(sub("^<label[^>]*>([^<]*)</label>$", "\\1", labels),
        sub("^.*for=\"([^\"]+)\".*$", "\\1", labels))
    boxes <- regmatches(html, gregexpr("<input id=\"[^\"]+\" type=\"checkbox\"[^>]*/>[[:space:]]*<span>[^<]*</span>", html))[[1L]]
    box_text <- stats::setNames(sub("^.*<span>([^<]*)</span>$", "\\1", boxes), sub("^<input id=\"([^\"]+)\".*$", "\\1", boxes))
    buttons <- regmatches(html, gregexpr("(?s)<(button|a) [^>]*id=\"[^\"]+\"[^>]*>.*?</(button|a)>", html, perl = TRUE))[[1L]]
    button_text <- stats::setNames(trimws(gsub("[[:space:]]+", " ", gsub("<[^>]+>", "", buttons))),
        sub("^<[a-z]+ [^>]*id=\"([^\"]+)\".*$", "\\1", buttons))
    areas <- regmatches(html, gregexpr("(?s)<textarea id=\"[^\"]+\"[^>]*>.*?</textarea>", html, perl = TRUE))[[1L]]
    area_text <- stats::setNames(sub("(?s)^<textarea[^>]*>(.*)</textarea>$", "\\1", areas, perl = TRUE),
        sub("^<textarea id=\"([^\"]+)\".*$", "\\1", areas))
    placeholders <- attr_of("placeholder")
    help <- regmatches(html, gregexpr("(?s)<span class=\"help-block\">.*?</span>", html, perl = TRUE))[[1L]]
    defaults <- omicsUiInputDefaults(ui)
    defaults <- vapply(defaults[order(names(defaults))], function(value) paste(format(value), collapse = ","), "")
    sorted <- function(x) if (length(x)) x[order(names(x))] else x
    list(html = digest::digest(html, algo = "sha256", serialize = FALSE), inputs = sort(unique(inputs)),
        outputs = sort(unique(ids[is_output])), testids = sort(unique(stats::na.omit(attr_of("data-testid")))),
        labels = sorted(c(label_text, box_text)), buttons = sorted(button_text), textareas = sorted(area_text),
        placeholders = stats::setNames(placeholders[!is.na(placeholders)], ids[!is.na(placeholders)]),
        conditions = Filter(Negate(is.na), attr_of("data-display-if")),
        help = sub("(?s)^<span class=\"help-block\">(.*)</span>$", "\\1", help, perl = TRUE),
        defaults = defaults)
}

# OMIC_SUMMARY_FROZEN_UI_START
omicSummaryW1FrozenUi <- list(
    proteomics = list(
        html = "5dfcdf6664a0ebbe518066b2b487855d39086e78cd1a07d40b78506b93ee3b38",
        inputs = c(
            "session_summary-copy_to_publication", "session_summary-description",
            "session_summary-download_report", "session_summary-enable_github",
            "session_summary-experiment_label", "session_summary-export_session_state",
            "session_summary-generate_report", "session_summary-github_email", "session_summary-github_org",
            "session_summary-github_username", "session_summary-project_id", "session_summary-push_to_github",
            "session_summary-save_workflow_args"),
        outputs = c(
            "session_summary-copy_status", "session_summary-session_summary",
            "session_summary-template_status"),
        testids = c(
            "prot-summary-copy-publication", "prot-summary-download-report", "prot-summary-export-state",
            "prot-summary-generate-report", "prot-summary-save-args"),
        labels = c(
            `session_summary-description` = "Description:",
            `session_summary-enable_github` = "Enable GitHub Push",
            `session_summary-experiment_label` = "Experiment Label:",
            `session_summary-github_email` = "GitHub Email:",
            `session_summary-github_org` = "GitHub Organization:",
            `session_summary-github_username` = "GitHub Username:",
            `session_summary-project_id` = "Project ID:"),
        buttons = c(
            `session_summary-copy_to_publication` = "Copy to Publication Directory",
            `session_summary-download_report` = "Download Report",
            `session_summary-export_session_state` = "Export Session State (.RDS)",
            `session_summary-generate_report` = "Generate Report",
            `session_summary-push_to_github` = "Push to GitHub",
            `session_summary-save_workflow_args` = "Save Workflow Arguments"),
        textareas = c(
            `session_summary-description` = "Full protein analysis workflow with config parameters"),
        placeholders = c(`session_summary-experiment_label` = "e.g., my_proteomics_analysis"),
        conditions = c(
            "output[&#39;session_summary-report_ready&#39;]",
            "input[&#39;session_summary-enable_github&#39;]"),
        help = character(0),
        defaults = c(
            `session_summary-enable_github` = "FALSE", `session_summary-experiment_label` = "",
            `session_summary-github_email` = "", `session_summary-github_org` = "",
            `session_summary-github_username` = "", `session_summary-project_id` = "")),
    lipidomics = list(
        html = "1b542005d98447dadd3b5cc370a79882eb18b21e9536ae904981a0ffb05fe8c8",
        inputs = c(
            "summary-copy_to_publication", "summary-description", "summary-download_report",
            "summary-enable_github", "summary-experiment_label", "summary-export_session_state",
            "summary-generate_report", "summary-github_email", "summary-github_org",
            "summary-github_username", "summary-project_id", "summary-push_to_github",
            "summary-save_workflow_args"),
        outputs = c("summary-copy_status", "summary-session_summary", "summary-template_status"),
        testids = c(
            "lipid-summary-copy-publication", "lipid-summary-download-report", "lipid-summary-export-state",
            "lipid-summary-generate-report", "lipid-summary-save-args"),
        labels = c(
            `summary-description` = "Description:", `summary-enable_github` = "Enable GitHub Push",
            `summary-experiment_label` = "Experiment Label:", `summary-github_email` = "GitHub Email:",
            `summary-github_org` = "GitHub Organization:", `summary-github_username` = "GitHub Username:",
            `summary-project_id` = "Project ID:"),
        buttons = c(
            `summary-copy_to_publication` = "Copy to Publication Directory",
            `summary-download_report` = "Download Report",
            `summary-export_session_state` = "Export Session State (.RDS)",
            `summary-generate_report` = "Generate Report", `summary-push_to_github` = "Push to GitHub",
            `summary-save_workflow_args` = "Save Workflow Arguments"),
        textareas = c(
            `summary-description` = "Full lipidomics analysis workflow with normalization and DA"),
        placeholders = c(`summary-experiment_label` = "e.g., my_lipidomics_analysis"),
        conditions = c(
            "output[&#39;summary-report_ready&#39;]", "input[&#39;summary-enable_github&#39;]"),
        help = character(0),
        defaults = c(
            `summary-enable_github` = "FALSE", `summary-experiment_label` = "", `summary-github_email` = "",
            `summary-github_org` = "", `summary-github_username` = "", `summary-project_id` = "")),
    metabolomics = list(
        html = "c4cbcb157dbd109100ff2df28df0449bc25f0629ed47d8895b79121d7d24ecb1",
        inputs = c(
            "summary-copy_to_publication", "summary-description", "summary-download_report",
            "summary-enable_github", "summary-experiment_label", "summary-export_session_state",
            "summary-generate_report", "summary-github_email", "summary-github_org",
            "summary-github_username", "summary-project_id", "summary-push_to_github",
            "summary-save_workflow_args"),
        outputs = c("summary-copy_status", "summary-session_summary", "summary-template_status"),
        testids = c(
            "metab-summary-copy-publication", "metab-summary-download-report", "metab-summary-export-state",
            "metab-summary-generate-report", "metab-summary-save-args"),
        labels = c(
            `summary-description` = "Description:", `summary-enable_github` = "Enable GitHub Push",
            `summary-experiment_label` = "Experiment Label:", `summary-github_email` = "GitHub Email:",
            `summary-github_org` = "GitHub Organization:", `summary-github_username` = "GitHub Username:",
            `summary-project_id` = "Project ID:"),
        buttons = c(
            `summary-copy_to_publication` = "Copy to Publication Directory",
            `summary-download_report` = "Download Report",
            `summary-export_session_state` = "Export Session State (.RDS)",
            `summary-generate_report` = "Generate Report", `summary-push_to_github` = "Push to GitHub",
            `summary-save_workflow_args` = "Save Workflow Arguments"),
        textareas = c(
            `summary-description` = "Full metabolomics analysis workflow with normalization and DA"),
        placeholders = c(`summary-experiment_label` = "e.g., my_metabolomics_analysis"),
        conditions = c(
            "output[&#39;summary-report_ready&#39;]", "input[&#39;summary-enable_github&#39;]"),
        help = character(0),
        defaults = c(
            `summary-enable_github` = "FALSE", `summary-experiment_label` = "", `summary-github_email` = "",
            `summary-github_org` = "", `summary-github_username` = "", `summary-project_id` = "")),
    transcriptomics = list(
        html = "b73ded6526ca76899e49792601f81643e58fb2a7c8a553baa378ae7f2ed6eb81",
        inputs = c(
            "summary-copy_to_publication", "summary-description", "summary-download_report",
            "summary-enable_github", "summary-experiment_label", "summary-export_session_state",
            "summary-generate_report", "summary-github_email", "summary-github_org",
            "summary-github_username", "summary-project_id", "summary-push_to_github",
            "summary-run_ruviii_export", "summary-save_workflow_args"),
        outputs = c(
            "summary-copy_status", "summary-ruviii_before_after", "summary-ruviii_m_ui",
            "summary-ruviii_status", "summary-session_summary", "summary-template_status"),
        testids = c(
            "trx-summary-copy-publication", "trx-summary-download-report", "trx-summary-export-state",
            "trx-summary-generate-report", "trx-summary-ruviii-export", "trx-summary-save-args"),
        labels = c(
            `summary-description` = "Description:", `summary-enable_github` = "Enable GitHub Push",
            `summary-experiment_label` = "Experiment Label:", `summary-github_email` = "GitHub Email:",
            `summary-github_org` = "GitHub Organization:", `summary-github_username` = "GitHub Username:",
            `summary-project_id` = "Project ID:"),
        buttons = c(
            `summary-copy_to_publication` = "Copy to Publication Directory",
            `summary-download_report` = "Download Report",
            `summary-export_session_state` = "Export Session State (.RDS)",
            `summary-generate_report` = "Generate Report", `summary-push_to_github` = "Push to GitHub",
            `summary-run_ruviii_export` = "Run RUVIII Export",
            `summary-save_workflow_args` = "Save Workflow Arguments"),
        textareas = c(
            `summary-description` = "Full transcriptomics analysis workflow with normalization, DE and the RUVIII integration export"),
        placeholders = c(`summary-experiment_label` = "e.g., my_transcriptomics_analysis"),
        conditions = c(
            "output[&#39;summary-report_ready&#39;]", "input[&#39;summary-enable_github&#39;]"),
        help = "\n          RUVIII (ruv::RUVIII) on the TMM logCPM, with the RUVs inferred negative controls and k, for the\n          integration export only: the DA results are not changed. M is the replicate structure: the group (the notebook),\n          or the design's technical replicates when it has them.\n        ",
        defaults = c(
            `summary-enable_github` = "FALSE", `summary-experiment_label` = "", `summary-github_email` = "",
            `summary-github_org` = "", `summary-github_username` = "", `summary-project_id` = "")))
# OMIC_SUMMARY_FROZEN_UI_END

test_that("(a) the four summary UIs keep their HTML, ids, test ids, labels, buttons, placeholders, texts, conditions and defaults", {
    uis <- omicSummaryW1Uis()
    actual <- lapply(uis, omicSummaryW1UiParts)
    if (!is.null(omicSummaryW1DumpPath())) omicSummaryW1Dump("omicSummaryW1FrozenUi", actual)
    expect_named(actual, names(omicSummaryW1FrozenUi))
    for (omic in names(omicSummaryW1FrozenUi)) {
        for (part in names(omicSummaryW1FrozenUi[[omic]])) {
            expect_identical(actual[[omic]][[part]], omicSummaryW1FrozenUi[[omic]][[part]], info = paste(omic, part))
        }
    }
})

# ------------------------------------------------------------------------------------------------
# (b) The server, driven in testServer.
# ------------------------------------------------------------------------------------------------

#' The browser, notification and progress doubles and the reactive values' recorder, scoped to `env`. A notification
#' is kept as "[type duration id] message" (shiny's default duration 5 when the caller gives none), an input update as
#' "<id>: <choices> -> <value>" (queued for the browser double), the progress as "start <message> <value>" and
#' "<amount> <detail>". The `reactiveValues` the tab creates with `workflow_args_saved` is kept on `acc$values`.
omicSummaryW1InstallDoubles <- function(log, acc, env = parent.frame()) {
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
    values <- get("reactiveValues", envir = asNamespace("shiny"))
    bindings <- list(
        updateSelectInput = function(session = domain(), inputId, label = NULL, choices = NULL, selected = NULL) {
            queue(session, inputId, selected, choices)
        },
        updateRadioButtons = function(session = domain(), inputId, label = NULL, choices = NULL, selected = NULL,
            inline = FALSE, choiceNames = NULL, choiceValues = NULL) queue(session, inputId, selected, choices %||% choiceValues),
        updateTextAreaInput = function(session = domain(), inputId, label = NULL, value = NULL, placeholder = NULL) {
            queue(session, inputId, value)
        },
        updateTextInput = function(session = domain(), inputId, label = NULL, value = NULL, placeholder = NULL) {
            queue(session, inputId, value)
        },
        updateCheckboxInput = function(session = domain(), inputId, label = NULL, value = NULL) queue(session, inputId, value),
        updateNumericInput = function(session = domain(), inputId, label = NULL, value = NULL, min = NULL, max = NULL,
            step = NULL) queue(session, inputId, value),
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
        },
        reactiveValues = function(...) {
            created <- values(...)
            if ("workflow_args_saved" %in% names(list(...))) acc$values <- created
            created
        })
    do.call(testthat::local_mocked_bindings, c(bindings, list(.package = "shiny", .env = env)))
    imported <- intersect(names(bindings), unlist(getNamespaceImports("MultiScholaR")$shiny, use.names = FALSE))
    if (length(imported)) do.call(testthat::local_mocked_bindings, c(bindings[imported], list(.package = "MultiScholaR", .env = env)))
    invisible(log)
}

#' A short content digest (the entry suite's canonical form).
omicSummaryW1Short <- function(x) substr(guiEntryDigest(x), 1L, 16L)

#' A value with its process-local parts masked, for a digest: each saved store reference as `<ref>`, each string with
#' the entry suite's masks, each time as `<time>` (the canonical form), each environment by its class.
omicSummaryW1Clean <- function(x, roots, depth = 0L) {
    if (depth > 14L) return("<deep>")
    if (isSavedWorkflowReference(x)) return("<ref>")
    if (is.environment(x)) return(paste0("<", class(x)[[1L]], ">"))
    if (is.list(x) && !is.data.frame(x) && !isS4(x) && !is.object(x)) {
        return(lapply(x, omicSummaryW1Clean, roots = roots, depth = depth + 1L))
    }
    if (isS4(x) && length(methods::slotNames(x))) {
        slots <- methods::slotNames(x)
        return(c(list(`<class>` = class(x)[[1L]]), stats::setNames(lapply(slots, function(s)
            omicSummaryW1Clean(methods::slot(x, s), roots, depth + 1L)), slots)))
    }
    if (is.character(x) && is.null(dim(x)) && !is.object(x)) {
        return(stats::setNames(vapply(x, function(value) if (is.na(value)) NA_character_ else guiEntryMask(value, roots), "",
            USE.NAMES = FALSE), names(x)))
    }
    x
}

omicSummaryW1Digest <- function(x, roots) omicSummaryW1Short(omicSummaryW1Clean(x, roots))

#' A message as captured: the entry suite's masks, then without terminal styling or non-ASCII characters.
omicSummaryW1Text <- function(x, roots) {
    x <- gsub("\033\\[[0-9;]*m", "", guiEntryMask(x, roots))
    iconv(x, "UTF-8", "ASCII", sub = "")
}

#' An argument as text: NULL, a scalar string or flag (masked), a table (its rows, columns and digest), an S4 object
#' (its class and digest), a saved reference, an environment, or a digest.
omicSummaryW1Arg <- function(value, roots) {
    if (is.null(value)) return("<NULL>")
    if (isSavedWorkflowReference(value)) return("<ref>")
    if (is.environment(value) || shiny::is.reactivevalues(value)) return(paste0("<", class(value)[[1L]], ">"))
    if (is.data.frame(value)) return(sprintf("<table %dx%d %s>", nrow(value), ncol(value), omicSummaryW1Short(value)))
    if (isS4(value)) return(sprintf("<%s %s>", class(value)[[1L]], omicSummaryW1Digest(value, roots)))
    if (is.atomic(value) && length(value) <= 3L && !is.object(value)) return(omicSummaryW1Text(paste(format(value), collapse = "|"), roots))
    paste("<digest>", omicSummaryW1Digest(value, roots))
}

#' The `project_dirs` global: absent, the session's, or another value.
omicSummaryW1Global <- function(project_dirs) {
    if (!exists("project_dirs", envir = .GlobalEnv, inherits = FALSE)) return("absent")
    if (identical(get("project_dirs", envir = .GlobalEnv, inherits = FALSE), project_dirs)) "the session's" else "another"
}

#' Recorders on the collaborators the tab calls, replaced by name in the package namespace and scoped to `env`: the
#' copy, the parameters' writer, the managed product, the record and its bundle call the real functions (the record and
#' the bundle fail instead while `acc$faults` says so); the render double writes `<template>_<omic>_<label>.docx` in
#' the results summary when the project's `scripts/<omic>/<template>` exists and fails otherwise; the push double
#' records its call. With `fallback`, the package's report templates are absent (`system.file`) and the download double
#' fails offline.
omicSummaryW1InstallRecorders <- function(acc, setup, fallback = FALSE, env = parent.frame()) {
    ns <- asNamespace("MultiScholaR")
    roots <- setup$roots
    copy <- get("copyToResultsSummary", envir = ns)
    write_args <- get("createWorkflowArgsFromConfig", envir = ns)
    product <- get("omicsManagedProduct", envir = ns)
    record <- get("writeWorkflowAnalysisRecord", envir = ns)
    bundle <- get("exportWorkflowAnalysisRecordBundle", envir = ns)
    call <- function(name, ...) {
        args <- list(...)
        acc$calls <- c(acc$calls, paste0(name, "(", paste(names(args), vapply(args, omicSummaryW1Arg, "", roots = roots),
            sep = "=", collapse = ", "), ")"))
    }
    paths <- setup$paths
    bindings <- list(
        copyToResultsSummary = function(omic_type, experiment_label, contrasts_tbl = NULL, design_matrix = NULL,
            force = FALSE, current_rmd = NULL, project_dirs_object_name = "project_dirs") {
            call("copyToResultsSummary", omic_type = omic_type, experiment_label = experiment_label,
                contrasts_tbl = contrasts_tbl, design_matrix = design_matrix, force = force, current_rmd = current_rmd,
                project_dirs_object_name = project_dirs_object_name, global = omicSummaryW1Global(setup$project_dirs))
            failures <- copy(omic_type = omic_type, experiment_label = experiment_label, contrasts_tbl = contrasts_tbl,
                design_matrix = design_matrix, force = force, current_rmd = current_rmd,
                project_dirs_object_name = project_dirs_object_name)
            shown <- if (is.list(failures)) vapply(failures, function(f) as.character(f$display_name %||% "?"), "") else
                omicSummaryW1Arg(failures, roots)
            acc$calls <- c(acc$calls, paste("  copyToResultsSummary returned:", paste(shown, collapse = " | ")))
            invisible(failures)
        },
        RenderReport = function(omic_type, experiment_label, rmd_filename = "DIANN_report.rmd",
            project_dirs_object_name = "project_dirs", output_format = NULL) {
            call("RenderReport", omic_type = omic_type, experiment_label = experiment_label, rmd_filename = rmd_filename,
                project_dirs_object_name = project_dirs_object_name, output_format = output_format,
                global = omicSummaryW1Global(setup$project_dirs))
            template <- file.path(paths$base_dir, "scripts", omic_type, rmd_filename)
            if (!file.exists(template)) stop("kit108 render double: no template at ", template)
            dir.create(paths$results_summary_dir, recursive = TRUE, showWarnings = FALSE)
            rendered <- file.path(paths$results_summary_dir, sprintf("%s_%s_%s.docx", tools::file_path_sans_ext(rmd_filename),
                omic_type, experiment_label))
            writeLines("kit108 render double", rendered)
            rendered
        },
        createWorkflowArgsFromConfig = function(workflow_name, description = "", organism_name = NULL, taxon_id = NULL,
            source_dir_path = NULL, final_s4_object = NULL, contrasts_tbl = NULL, workflow_data = NULL) {
            call("createWorkflowArgsFromConfig", workflow_name = workflow_name, description = description,
                organism_name = organism_name, taxon_id = taxon_id, source_dir_path = source_dir_path,
                final_s4_object = final_s4_object, contrasts_tbl = contrasts_tbl, workflow_data = workflow_data)
            write_args(workflow_name = workflow_name, description = description, organism_name = organism_name,
                taxon_id = taxon_id, source_dir_path = source_dir_path, final_s4_object = final_s4_object,
                contrasts_tbl = contrasts_tbl, workflow_data = workflow_data)
        },
        omicsManagedProduct = function(operation, parameters, deliver, before = NULL) {
            call("omicsManagedProduct", operation = operation, parameters = paste(names(parameters),
                vapply(parameters, function(p) if (is.character(p) && grepl("^[0-9a-f-]{20,}$", p)) "<entity>" else
                    omicSummaryW1Arg(p, roots), ""), sep = ":", collapse = ";"))
            product(operation, parameters, deliver, before)
        },
        writeWorkflowAnalysisRecord = function(manager, path, now = Sys.time()) {
            call("writeWorkflowAnalysisRecord", manager = manager, path = path)
            if (isTRUE(acc$faults$record)) stop("kit108 record double fails")
            record(manager, path, now)
        },
        exportWorkflowAnalysisRecordBundle = function(manager, directory) {
            call("exportWorkflowAnalysisRecordBundle", manager = manager, directory = directory)
            if (isTRUE(acc$faults$bundle)) stop("kit108 bundle double fails")
            bundle(manager, directory)
        },
        pushProjectToGithubFromDirs = function(project_dirs, omic_type, experiment_label, project_id,
            github_org = getOption("github_org", "your_organization_here"),
            github_user_email = getOption("github_user_email", "your_email@example.com"),
            github_user_name = getOption("github_user_name", "your_github_username"),
            commit_message = "Initial project commit", private = TRUE) {
            call("pushProjectToGithubFromDirs", project_dirs = if (identical(project_dirs, setup$project_dirs))
                "the session's" else "another", omic_type = omic_type, experiment_label = experiment_label,
                project_id = project_id, github_org = github_org, github_user_email = github_user_email,
                github_user_name = github_user_name)
            invisible(TRUE)
        })
    do.call(testthat::local_mocked_bindings, c(bindings, list(.package = "MultiScholaR", .env = env)))
    if (!isTRUE(fallback)) return(invisible(acc))
    # The GitHub fallback: the package has no report template; the download fails offline.
    download <- function(url, destfile, method, quiet = FALSE, mode = "w", cacheOK = TRUE, extra = NULL, headers = NULL, ...) {
        call("download.file", url = url, destfile = destfile, quiet = quiet)
        stop("kit108 download double: offline")
    }
    testthat::local_mocked_bindings(download.file = download, .package = "MultiScholaR", .env = env)
    .install_ns_binding(asNamespace("utils"), "download.file", download, .local_envir = env)
    system_file <- get("system.file", envir = baseenv())
    absent <- function(..., package = "base", lib.loc = NULL, mustWork = FALSE) {
        parts <- c(...)
        if (identical(package, "MultiScholaR") && length(parts) && identical(parts[[1L]], "reports")) {
            call("system.file", path = paste(parts, collapse = "/"), package = package)
            return("")
        }
        system_file(..., package = package, lib.loc = lib.loc, mustWork = mustWork)
    }
    imports <- parent.env(ns)
    if (exists("system.file", envir = imports, inherits = FALSE)) {
        shim <- get("system.file", envir = imports)
        .install_ns_binding(imports, "system.file", function(..., package = "base", lib.loc = NULL, mustWork = FALSE) {
            parts <- c(...)
            if (identical(package, "MultiScholaR") && length(parts) && identical(parts[[1L]], "reports")) {
                return(absent(..., package = package, lib.loc = lib.loc, mustWork = mustWork))
            }
            shim(..., package = package, lib.loc = lib.loc, mustWork = mustWork)
        }, .local_envir = env)
    }
    .install_ns_binding(.BaseNamespaceEnv, "system.file", absent, .local_envir = env)
    invisible(acc)
}

#' The project's files under the setup's base directory (the entry suite's artifacts: table digests, others present),
#' `study_parameters.txt` (a digest of its masked text and its line count), the integration files (name; an RDS's
#' class, `<reference>` for a saved reference; a table's digest) and the session-state RDS's fields (times masked).
omicSummaryW1Files <- function(setup) {
    roots <- setup$roots
    paths <- setup$paths
    files <- unlist(guiEntryArtifacts(paths$base_dir, roots))
    study <- file.path(paths$source_dir, "study_parameters.txt")
    study_text <- if (file.exists(study)) {
        lines <- readLines(study, warn = FALSE)
        sprintf("%d lines %s", length(lines), substr(digest::digest(omicSummaryW1Text(lines, roots), algo = "sha256",
            serialize = FALSE), 1L, 16L))
    } else "<absent>"
    integration <- list.files(paths$integration_dir, full.names = TRUE)
    integration <- stats::setNames(vapply(integration, function(path) {
        if (grepl("\\.rds$", path, ignore.case = TRUE)) {
            value <- readRDS(path)
            if (is.null(value)) return("<RDS> NULL")
            if (isSavedWorkflowReference(value)) return("<RDS> <reference>")
            return(paste("<RDS>", class(value)[[1L]], omicSummaryW1Digest(value, roots)))
        }
        if (grepl("\\.(tsv|csv)$", path)) return(substr(digest::digest(guiEntryMask(readLines(path, warn = FALSE), roots),
            algo = "sha256", serialize = FALSE), 1L, 16L))
        "<present>"
    }, "", USE.NAMES = FALSE), vapply(basename(integration), guiEntryMask, "", roots = roots, USE.NAMES = FALSE))
    state <- list.files(paths$source_dir, pattern = "^session_state_.*\\.RDS$", full.names = TRUE)
    session_state <- if (length(state) == 1L) {
        value <- readRDS(state)
        vapply(value, function(field) {
            if (inherits(field, c("POSIXt", "Date"))) return("<time>")
            omicSummaryW1Arg(field, roots)
        }, "")
    } else sprintf("<%d files>", length(state))
    list(files = stats::setNames(as.character(files), names(files)), study_parameters = study_text,
        integration = integration, session_state = session_state)
}

#' Drive one session: the tab (and trx's normalisation and DA tabs) in testServer on the fixture, the UI's defaults set
#' as the browser sends them, the steps run (named functions of the context) and captured after each.
#' @param setup The session's fixture (`omicSummaryW1*Setup`).
#' @return The captures, by step.
omicSummaryW1Drive <- function(setup, steps) {
    log <- new.env(parent = emptyenv())
    log$notes <- character()
    log$updates <- character()
    log$progress <- character()
    log$pending <- list()
    acc <- new.env(parent = emptyenv())
    acc$calls <- character()
    acc$captures <- list()
    acc$previous <- list()
    acc$errors <- character()
    acc$faults <- list()
    acc$extra <- list()
    omicSummaryW1InstallDoubles(log, acc)
    omicSummaryW1InstallRecorders(acc, setup, fallback = isTRUE(setup$fallback))
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true", http_proxy = "http://127.0.0.1:9", https_proxy = "http://127.0.0.1:9",
        HTTP_PROXY = "http://127.0.0.1:9", HTTPS_PROXY = "http://127.0.0.1:9", no_proxy = "", NO_PROXY = "")
    withr::local_options(github_org = NULL, github_user_email = NULL, github_user_name = NULL)
    withr::local_seed(20261109L)
    had_global <- exists("project_dirs", envir = .GlobalEnv, inherits = FALSE)
    old_global <- if (had_global) get("project_dirs", envir = .GlobalEnv, inherits = FALSE)
    withr::defer(if (had_global) assign("project_dirs", old_global, envir = .GlobalEnv) else
        if (exists("project_dirs", envir = .GlobalEnv, inherits = FALSE)) rm("project_dirs", envir = .GlobalEnv))
    if (identical(setup$global, "absent") && exists("project_dirs", envir = .GlobalEnv, inherits = FALSE)) {
        rm("project_dirs", envir = .GlobalEnv)
    }
    if (identical(setup$global, "another")) assign("project_dirs", list(elsewhere = list()), envir = .GlobalEnv)
    if (identical(setup$global, "the session's")) assign("project_dirs", setup$project_dirs, envir = .GlobalEnv)
    prefix <- paste0(setup$ns, "-")
    server <- local({
        mount <- setup$server(acc)
        function(input, output, session) {
            session$onUnhandledError(function(error) acc$errors <- c(acc$errors, paste("<unhandled>", conditionMessage(error))))
            mount(input, output, session)
        }
    })
    text <- function(x) omicSummaryW1Text(x, setup$roots)
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
        set <- function(..., .ns = prefix) {
            values <- list(...)
            names(values) <- paste0(.ns, names(values))
            put(values)
        }
        read <- function(name) tryCatch(output[[paste0(prefix, name)]], error = function(e) {
            if (inherits(e, c("shiny.silent.error", "shiny.output.cancel"))) "<req>" else paste("<error>", text(conditionMessage(e)))
        })
        ui_defaults <- omicsUiInputDefaults(setup$ui("x"))
        names(ui_defaults) <- paste0(prefix, sub("^x-", "", names(ui_defaults)))
        put(ui_defaults)
        # The description's text area, which the browser sends too.
        set(description = setup$description)
        apply_updates()
        capture <- function(step) {
            session$flushReact()
            apply_updates()
            acc$global_after <- omicSummaryW1Global(setup$project_dirs)
            out_names <- sort(names(session$.__enclos_env__$private$outs))
            out_names <- sub(paste0("^", prefix), "", out_names[startsWith(out_names, prefix)])
            shown <- vapply(out_names, function(name) {
                if (startsWith(name, "download_")) return("<download>")
                guiEntryOutputValue(read(name), setup$roots)
            }, "")
            values <- if (shiny::is.reactivevalues(acc$values)) shiny::isolate(shiny::reactiveValuesToList(acc$values)) else list()
            values <- vapply(values[order(names(values))], function(value) omicSummaryW1Arg(value, setup$roots), "")
            files <- omicSummaryW1Files(setup)
            record <- c(list(
                outputs = out_names,
                values = shown,
                notes = vapply(log$notes, text, "", USE.NAMES = FALSE),
                updates = vapply(log$updates, text, "", USE.NAMES = FALSE),
                progress = vapply(log$progress, text, "", USE.NAMES = FALSE),
                reactive_values = values,
                # The calls' arguments are masked by `omicSummaryW1Arg`; their digests are not passed through the masks.
                calls = iconv(acc$calls, "UTF-8", "ASCII", sub = ""),
                global = c(before = acc$global_before, after = acc$global_after),
                files = files$files,
                study_parameters = files$study_parameters,
                integration = files$integration,
                session_state = files$session_state,
                errors = vapply(acc$errors, text, "", USE.NAMES = FALSE)), acc$extra)
            # The files after the first step are the step's changes: each file added or changed with its value, each one
            # removed as "<removed>".
            current <- record$files
            previous <- acc$previous_files
            acc$previous_files <- current
            if (!is.null(previous)) {
                differs <- vapply(names(current), function(name) !identical(current[[name]], previous[name][[1L]]), TRUE)
                changed <- names(current)[differs]
                removed <- setdiff(names(previous), names(current))
                record$files <- c(current[changed], stats::setNames(rep("<removed>", length(removed)), removed))
                if (!length(record$files)) record$files <- "<unchanged>"
            }
            # A field the step left as the previous step had it is held as "<unchanged>".
            for (field in c("study_parameters", "integration", "session_state")) {
                value <- record[[field]]
                if (identical(value, acc$previous[[field]])) record[[field]] <- "<unchanged>"
                acc$previous[[field]] <- value
            }
            log$notes <- character()
            log$updates <- character()
            log$progress <- character()
            acc$calls <- character()
            acc$errors <- character()
            acc$extra <- list()
            acc$captures[[step]] <- record
        }
        ctx <- list(set = set, read = read, setup = setup, bus = setup$bus, acc = acc,
            fault = function(name, value = TRUE) acc$faults[[name]] <- value,
            keep = function(name, value) acc$extra[[name]] <- value)
        for (step in names(steps)) {
            acc$global_before <- omicSummaryW1Global(setup$project_dirs)
            steps[[step]](ctx)
            capture(step)
        }
    }))))
    acc$captures
}

# ------------------------------------------------------------------------------------------------
# The fixtures.
# ------------------------------------------------------------------------------------------------

#' The masks' roots: the project's base directory and the owner's root, each as given and normalised.
omicSummaryW1Roots <- function(...) {
    roots <- unlist(list(...))
    roots <- roots[!is.na(roots) & nzchar(roots)]
    unique(c(normalizePath(roots, mustWork = FALSE), roots))
}

#' Lipid's (MCI-023) and metab's (MCI-017) fixture: the module CI workflow (a bound owner holding the lane's object up
#' to its correlation-filtered state, the contrasts, design and parameter bus fields), the module CI project paths, the
#' publication inputs (with or without the optional artifacts, or none), the project's template when `template`, and
#' the tab mounted as the shell mounts it (`server(module_id, project_dirs, omic, label, data)`), managed when
#' `managed` (the runner with the in-test worker, `omicsManagedChainMocks`).
#' @param inputs `complete`, `optional_missing` or `none` (and no contrasts or design on the bus).
#' @param global The `project_dirs` global at the start: `absent`, `another` or `the session's`.
omicSummaryW1FeatureSetup <- function(omic, inputs = "complete", template = FALSE, managed = FALSE, label = "kit108",
    global = "absent", fallback = FALSE, env = parent.frame()) {
    spec <- switch(omic,
        lipid = list(omic_type = "lipidomics", config = module_ci_lipid_summary_config("pos_neg"),
            workflow = module_ci_lipid_summary_workflow, paths = module_ci_lipid_summary_paths,
            inputs = module_ci_lipid_summary_write_publication_inputs, template = module_ci_lipid_summary_write_template,
            server = mod_lipid_summary_server, ui = mod_lipid_summary_ui,
            description = "Full lipidomics analysis workflow with normalization and DA"),
        metab = list(omic_type = "metabolomics", config = module_ci_metab_summary_config("combined"),
            workflow = module_ci_metab_summary_workflow, paths = module_ci_metab_summary_paths,
            inputs = module_ci_metab_summary_write_publication_inputs, template = module_ci_metab_summary_write_template,
            server = mod_metab_summary_server, ui = mod_metab_summary_ui,
            description = "Full metabolomics analysis workflow with normalization and DA"))
    bus <- spec$workflow(spec$config)
    project_dirs <- spec$paths()
    paths <- project_dirs[[spec$omic_type]]
    withr::defer(unlink(paths$base_dir, recursive = TRUE), envir = env)
    if (!identical(inputs, "none")) spec$inputs(project_dirs, spec$config, include_optional = identical(inputs, "complete"))
    if (identical(inputs, "none")) {
        bus$contrasts_tbl <- NULL
        bus$design_matrix <- NULL
    }
    if (isTRUE(template)) spec$template(project_dirs)
    owner_root <- dirname(bus$state_manager$.context$experiment_paths$source_dir)
    probe <- NULL
    if (managed) {
        probe <- new.env()
        omicsManagedChainMocks(list(owner = bus$state_manager, workflow = bus, paths = paths, source_dir = paths$source_dir),
            probe, env = env)
    }
    server <- function(acc) {
        function(input, output, session) {
            if (managed) {
                runner <- installWorkflowOmicsActions(session, bus, spec$omic_type, "kit108", paths)
                runner$progress <- module_ci_trx_quiet_progress
            }
            spec$server("summary", project_dirs, spec$omic_type, label, bus)
        }
    }
    list(omic = omic, server = server, ns = "summary", ui = spec$ui, bus = function() bus, project_dirs = project_dirs,
        paths = paths, roots = omicSummaryW1Roots(paths$base_dir, owner_root), global = global, fallback = fallback,
        description = spec$description, probe = probe)
}

#' Trx's fixture: MCI-032's project (`module_ci_trx_project(module_ci_trx_corpus())`, at `trx_qc_complete`, with its
#' technical replicates when `technical`) with the normalisation, DA and summary tabs under the transcriptomics profile,
#' as the shell mounts them (the summary with the profile; MCI-032's journey), managed when `managed`.
omicSummaryW1TrxSetup <- function(technical = FALSE, managed = FALSE, label = "kit108", global = "absent",
    env = parent.frame()) {
    skipWithoutTrxDeclReference()
    fixture <- module_ci_trx_project(module_ci_trx_corpus(), technical = technical, env = env)
    paths <- fixture$paths
    project <- dirname(dirname(paths$source_dir))
    probe <- NULL
    if (managed) {
        probe <- new.env()
        omicsManagedChainMocks(fixture, probe, env = env)
    }
    server <- function(acc) {
        function(input, output, session) {
            if (managed) {
                runner <- installWorkflowOmicsActions(session, fixture$workflow, "transcriptomics", "kit108", paths)
                runner$progress <- module_ci_trx_quiet_progress
            }
            profile <- trxOmicProfile()
            mod_lipid_norm_server("norm", fixture$workflow, paths, "transcriptomics", "kit108", profile = profile)
            mod_lipid_da_server("de", fixture$workflow, paths, "transcriptomics", "kit108", profile = profile)
            mod_lipid_summary_server("summary", fixture$project_dirs, "transcriptomics", label, fixture$workflow,
                profile = profile)
        }
    }
    list(omic = "trx", server = server, ns = "summary", ui = function(id) mod_lipid_summary_ui(id, profile = trxOmicProfile()),
        bus = function() fixture$workflow, project_dirs = fixture$project_dirs, paths = paths,
        roots = omicSummaryW1Roots(project, paths$base_dir), global = global, fallback = FALSE, owner = fixture$owner,
        description = "Full transcriptomics analysis workflow with normalization, DE and the RUVIII integration export",
        probe = probe)
}

#' Proteomics' fixtures (MCI-011): the module CI workflow of a route (`branch`; its in-memory state manager, the
#' config list's workflow type and requested template, the contrasts, design and parameter bus fields) with the module
#' CI project paths and publication inputs; or, with `owner`, the same bus over a DIA owner holding MCI-009's object
#' (row names dropped, as the DIA reference's tables need) as `correlation_filtered`, registered as the session's
#' active owner (`session$userData$dia_capacity_manager`), native (`session_format = "native_v1"`, with no export) when
#' `native`.
omicSummaryW1ProtSetup <- function(branch = "DIA", requested = NULL, owner = FALSE, native = FALSE, template = FALSE,
    label = "kit108", global = "absent", fallback = FALSE, env = parent.frame()) {
    config <- module_ci_prot_summary_branch_config(branch)
    bus <- module_ci_prot_summary_workflow(config)
    if (!is.null(requested)) bus$config_list$globalParameters$report_template <- requested
    project_dirs <- module_ci_prot_summary_paths()
    paths <- project_dirs$proteomics
    withr::defer(unlink(paths$base_dir, recursive = TRUE), envir = env)
    module_ci_prot_summary_write_publication_inputs(project_dirs)
    if (isTRUE(template)) module_ci_prot_summary_write_template_set(project_dirs)
    holder <- new.env(parent = emptyenv())
    server <- function(acc) {
        function(input, output, session) {
            if (owner) {
                design <- module_ci_prot_da_design()
                rownames(design) <- NULL
                quant <- module_ci_prot_da_table()
                rownames(quant) <- NULL
                object <- module_ci_prot_da_object(workflow_type = config$workflow_type, use_limpa = config$use_limpa,
                    report_template = config$report_template, design = design, protein_quant_table = quant,
                    args = module_ci_prot_summary_args(config))
                store_root <- file.path(paths$base_dir, "store")
                dir.create(store_root, showWarnings = FALSE)
                manager <- newDIAStateAdapter(newWorkflowEntityStore(store_root, "kit108-prot", config$workflow_type),
                    WorkflowState$new(audit_enabled = FALSE))
                session$onSessionEnded(function() closeDIAStateAdapter(manager))
                manager$.context$reference_states <- TRUE
                manager$.context$experiment_paths <- paths
                if (native) manager$.context$session_format <- "native_v1"
                manager$saveState("correlation_filtered", object, list(), "MCI-009 fixture state")
                session$userData$dia_capacity_manager <- manager
                bus$state_manager <- manager
                holder$manager <- manager
            }
            mod_prot_summary_server("session_summary", project_dirs, "proteomics", label, bus)
        }
    }
    list(omic = "prot", server = server, ns = "session_summary", ui = mod_prot_summary_ui, bus = function() bus,
        project_dirs = project_dirs, paths = paths, roots = omicSummaryW1Roots(paths$base_dir), global = global,
        fallback = fallback, manager = function() holder$manager,
        description = "Full protein analysis workflow with config parameters")
}

# ------------------------------------------------------------------------------------------------
# The steps and cases.
# ------------------------------------------------------------------------------------------------

omicSummaryW1Start <- function(ctx) invisible(NULL)
omicSummaryW1Label <- function(value = "kit108") function(ctx) ctx$set(experiment_label = value)
omicSummaryW1Click <- function(button, n = 1L) function(ctx) do.call(ctx$set, stats::setNames(list(n), button))
omicSummaryW1Fault <- function(name, value, then) function(ctx) {
    ctx$fault(name, value)
    then(ctx)
}
omicSummaryW1Github <- function(enable, n = 1L) function(ctx) {
    ctx$set(enable_github = enable, github_org = "kit108-org", github_email = "kit108@example.org",
        github_username = "kit108-user", project_id = "kit108-project")
    ctx$set(push_to_github = n)
}

#' The refusals: save, copy, report and export without a label (silent `req`, the report's F-11 notice); the report
#' before the copy (F-11's notice); the push before a report (silent).
omicSummaryW1RefusalSteps <- function() {
    list(start = omicSummaryW1Start, save_no_label = omicSummaryW1Click("save_workflow_args"),
        copy_no_label = omicSummaryW1Click("copy_to_publication"), report_no_label = omicSummaryW1Click("generate_report"),
        export_no_label = omicSummaryW1Click("export_session_state"),
        report_before_copy = function(ctx) {
            ctx$set(experiment_label = "kit108")
            ctx$set(generate_report = 2L)
        },
        github_before_report = omicSummaryW1Github(TRUE))
}

#' The flow: save, copy, report, export, the push disabled (silent), then enabled (the push double).
omicSummaryW1MainSteps <- function() {
    list(start = omicSummaryW1Start, save = omicSummaryW1Click("save_workflow_args"),
        copy = omicSummaryW1Click("copy_to_publication"), report = omicSummaryW1Click("generate_report"),
        export = omicSummaryW1Click("export_session_state"), github_disabled = omicSummaryW1Github(FALSE),
        github = omicSummaryW1Github(TRUE, 2L))
}

#' The variants: the record failing at the save, the copy (an optional artifact missing, or the gate), the report with
#' the template in the project and the record bundle failing.
omicSummaryW1VariantSteps <- function() {
    list(start = omicSummaryW1Start, save_record_fails = omicSummaryW1Fault("record", TRUE, omicSummaryW1Click("save_workflow_args")),
        copy = omicSummaryW1Fault("record", FALSE, omicSummaryW1Click("copy_to_publication")),
        report_bundle_fails = omicSummaryW1Fault("bundle", TRUE, omicSummaryW1Click("generate_report")))
}

#' The trx steps: the normalisation and the DA (MCI-032's journey inputs), then the tab's.
omicSummaryW1TrxNorm <- function(ctx) {
    ctx$set(ruv_mode = "automatic", .ns = "norm-")
    ctx$set(run_normalization = 1L, .ns = "norm-")
}
omicSummaryW1TrxDa <- function(ctx) {
    ctx$set(da_q_val_thresh = 0.05, treat_lfc_cutoff = 1, .ns = "de-")
    ctx$set(run_da_analysis = 1L, .ns = "de-")
}
omicSummaryW1Ruviii <- function(m_source, n) function(ctx) {
    ctx$set(ruviii_m_source = m_source)
    ctx$set(run_ruviii_export = n)
    export <- ctx$setup$owner$.context$snapshot$exports$trx_ruviii_export
    value <- if (!is.null(export)) hydrateSavedWorkflowEntity(ctx$setup$owner$.context$store, export$reference)$value
    ctx$keep("ruviii", if (is.null(value)) "<none>" else c(m_source = value$m_source, k = as.character(value$k),
        controls = as.character(value$control_count %||% sum(value$ctl_mask)),
        newY = sprintf("%dx%d %s", nrow(value$newY), ncol(value$newY), omicSummaryW1Short(value$newY))))
}

omicSummaryW1CaseSpec <- function(case, env) {
    feature <- function(omic, steps, ...) list(setup = omicSummaryW1FeatureSetup(omic, ..., env = env), steps = steps)
    prot <- function(steps, ...) list(setup = omicSummaryW1ProtSetup(..., env = env), steps = steps)
    trx_head <- list(norm = omicSummaryW1TrxNorm, da = omicSummaryW1TrxDa)
    switch(case,
        lipid_refusals = feature("lipid", omicSummaryW1RefusalSteps(), label = NULL),
        lipid_main = feature("lipid", omicSummaryW1MainSteps()),
        lipid_variants = feature("lipid", omicSummaryW1VariantSteps(), inputs = "optional_missing", template = TRUE,
            global = "the session's"),
        lipid_required_missing = feature("lipid", list(start = omicSummaryW1Start, copy = omicSummaryW1Click("copy_to_publication"),
            report = omicSummaryW1Click("generate_report")), inputs = "none"),
        lipid_github_fallback = feature("lipid", list(copy = omicSummaryW1Click("copy_to_publication"),
            report_fallback = omicSummaryW1Click("generate_report")), global = "the session's", fallback = TRUE),
        lipid_managed = feature("lipid", list(start = omicSummaryW1Start, save = omicSummaryW1Click("save_workflow_args"),
            export = omicSummaryW1Click("export_session_state")), managed = TRUE),
        # Another project's `project_dirs` global is kept (the tab assigns it only when absent): the copy reads it.
        lipid_global_another = feature("lipid", list(copy = omicSummaryW1Click("copy_to_publication")), global = "another"),
        metab_refusals =feature("metab", omicSummaryW1RefusalSteps(), label = NULL),
        metab_main = feature("metab", omicSummaryW1MainSteps()),
        metab_variants = feature("metab", omicSummaryW1VariantSteps(), inputs = "optional_missing", template = TRUE,
            global = "the session's"),
        metab_required_missing = feature("metab", list(start = omicSummaryW1Start, copy = omicSummaryW1Click("copy_to_publication"),
            report = omicSummaryW1Click("generate_report")), inputs = "none"),
        metab_github_fallback = feature("metab", list(copy = omicSummaryW1Click("copy_to_publication"),
            report_fallback = omicSummaryW1Click("generate_report")), global = "the session's", fallback = TRUE),
        metab_managed = feature("metab", list(start = omicSummaryW1Start, save = omicSummaryW1Click("save_workflow_args"),
            export = omicSummaryW1Click("export_session_state")), managed = TRUE),
        trx_main = list(setup = omicSummaryW1TrxSetup(env = env), steps = c(trx_head, list(
            report_before_copy = function(ctx) {
                ctx$set(experiment_label = "kit108")
                ctx$set(generate_report = 1L)
            },
            save_before_export = omicSummaryW1Click("save_workflow_args"),
            ruviii_unavailable = omicSummaryW1Ruviii("technical_replicate", 1L),
            ruviii = omicSummaryW1Ruviii("group", 2L),
            save = omicSummaryW1Click("save_workflow_args", 2L),
            copy = omicSummaryW1Click("copy_to_publication"), report = omicSummaryW1Click("generate_report", 2L),
            export = omicSummaryW1Click("export_session_state")))),
        trx_technical = list(setup = omicSummaryW1TrxSetup(technical = TRUE, env = env), steps = c(trx_head, list(
            ruviii_technical = omicSummaryW1Ruviii("technical_replicate", 1L)))),
        trx_managed = list(setup = omicSummaryW1TrxSetup(managed = TRUE, env = env), steps = c(trx_head, list(
            ruviii = omicSummaryW1Ruviii("group", 1L), save = omicSummaryW1Click("save_workflow_args")))),
        prot_refusals = prot(omicSummaryW1RefusalSteps(), label = NULL),
        prot_route_dia = prot(omicSummaryW1MainSteps()[c("start", "save", "copy", "report", "export")], branch = "DIA"),
        prot_route_dia_limpa = prot(omicSummaryW1MainSteps()[c("start", "save", "copy", "report", "export")], branch = "DIA_limpa"),
        prot_route_tmt = prot(omicSummaryW1MainSteps()[c("start", "save", "copy", "report", "export")], branch = "TMT"),
        prot_route_lfq = prot(omicSummaryW1MainSteps()[c("start", "save", "copy", "report", "export")], branch = "MaxQuant_LFQ"),
        prot_stale = prot(list(copy = omicSummaryW1Click("copy_to_publication"), report = omicSummaryW1Click("generate_report"),
            export = omicSummaryW1Click("export_session_state")), branch = "TMT", requested = "DIANN_report.rmd"),
        prot_owner = prot(omicSummaryW1MainSteps(), owner = TRUE, template = TRUE, global = "the session's"),
        prot_owner_variants = prot(list(start = omicSummaryW1Start,
            save_record_fails = omicSummaryW1Fault("record", TRUE, omicSummaryW1Click("save_workflow_args")),
            copy_gate_refuses = function(ctx) {
                ctx$fault("record", FALSE)
                manager <- ctx$setup$manager()
                manager$.context$export_pending <- list(label = "session-export")
                ctx$set(copy_to_publication = 1L)
            },
            copy = function(ctx) {
                manager <- ctx$setup$manager()
                manager$.context$export_pending <- NULL
                ctx$set(copy_to_publication = 2L)
            },
            report_bundle_fails = omicSummaryW1Fault("bundle", TRUE, omicSummaryW1Click("generate_report"))),
            owner = TRUE, global = "the session's"),
        prot_native = prot(list(start = omicSummaryW1Start, save = omicSummaryW1Click("save_workflow_args")), owner = TRUE,
            native = TRUE),
        prot_github_fallback = prot(list(copy = omicSummaryW1Click("copy_to_publication"),
            report_fallback = omicSummaryW1Click("generate_report")), fallback = TRUE))
}

omicSummaryW1Cases <- c("lipid_refusals", "lipid_main", "lipid_variants", "lipid_required_missing", "lipid_github_fallback",
    "lipid_managed", "lipid_global_another", "metab_refusals", "metab_main", "metab_variants", "metab_required_missing", "metab_github_fallback",
    "metab_managed", "trx_main", "trx_technical", "trx_managed", "prot_refusals", "prot_route_dia", "prot_route_dia_limpa",
    "prot_route_tmt", "prot_route_lfq", "prot_stale", "prot_owner", "prot_owner_variants", "prot_native",
    "prot_github_fallback")

#' One case's captures (each case has its own fixture and session).
omicSummaryW1Case <- function(case) {
    env <- environment()
    spec <- omicSummaryW1CaseSpec(case, env)
    captures <- omicSummaryW1Drive(spec$setup, spec$steps)
    if (!is.null(omicSummaryW1DumpPath())) {
        cat(sprintf("omicSummaryW1Frozen[[\"%s\"]] <- %s\n\n", case, omicSummaryW1Code(captures)),
            file = omicSummaryW1DumpPath(), append = TRUE)
    }
    captures
}

# The frozen captures, per case (the dump's entries, pasted).
# OMIC_SUMMARY_FROZEN_START
omicSummaryW1Frozen <- list()
omicSummaryW1Frozen[["lipid_refusals"]] <- list(
    start = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = c(
            `da_output/de_negmode_lipids_groupKO-groupWT_long_annot.tsv` = "1da4796cae289522efd567b881e22cd01afe80ebea6b5caf33a7b10ac8d0b83b",
            `da_output/de_negmode_lipids_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_lipids_groupKO-groupWT_long_annot.tsv` = "568747756defd94e3305d975ce8b0de78947e9dd8bacc0c68bad0c67b77e2ca0",
            `da_output/de_posmode_lipids_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `feature_qc/itsd_lcms_pos.png` = "<present>",
            `feature_qc/lcms_neg_ruv_corrected_rle.png` = "<present>",
            `feature_qc/lcms_pos_pre_norm_pca.png` = "<present>",
            `feature_qc/lipidomics_composite_QC_figure.pdf` = "<present>",
            `feature_qc/lipidomics_composite_QC_figure.png` = "<present>",
            `feature_qc/ruv_normalised_results.RDS` = "<present>",
            `feature_qc/ruv_normalised_results.tsv` = "ab2c445d1a3d8f5ae94c31cd216f8330033e9b029afa2b6da7bbacfaf9743984",
            `publication_graphs/Heatmap/artifact.txt` = "<present>",
            `publication_graphs/Heatmaps/artifact.txt` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/lipids_num_sig_da_molecules.tab` = "ab123422c31c34ada4ea3ab9a0c9149177c0b9cdb09acdb41b16e9d4adfeed97",
            `publication_graphs/NumSigDeMolecules/artifact.txt` = "<present>",
            `publication_graphs/NumSigDeMolecules/lipids_num_sig_de_molecules.tab` = "ab123422c31c34ada4ea3ab9a0c9149177c0b9cdb09acdb41b16e9d4adfeed97",
            `publication_graphs/Volcano_Plots/artifact.txt` = "<present>",
            `source/contrasts_tbl.tab` = "7209220541b9240edc6b7b73030b7384687b68cd7814091847ac8baf36891f72",
            `source/design_matrix.tab` = "c3435c79da809286aebfaa17f024efdea7282282a0904b2f5f1477f4790b6983",
            `source/study_parameters.txt` = "<present>",
            `time/12_correlation_filtered_combined_plots.png` = "<present>"),
        study_parameters = "1 lines 55c1b3aadde3f3b7",
        integration = stats::setNames(character(0), character()),
        session_state = "<0 files>",
        errors = "<unhandled> unused arguments (self, name)"),
    save_no_label = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    copy_no_label = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    report_no_label = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = "[warning 5 -] Please enter an experiment label before generating the report.",
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    export_no_label = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    report_before_copy = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = "[warning 5 -] Please copy the results to the publication directory before generating the report (Copy to Publication Directory).",
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    github_before_report = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)))

omicSummaryW1Frozen[["lipid_main"]] <- list(
    start = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = character(0),
        updates = "summary-experiment_label: <same> -> kit108",
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = c(
            `da_output/de_negmode_lipids_groupKO-groupWT_long_annot.tsv` = "1da4796cae289522efd567b881e22cd01afe80ebea6b5caf33a7b10ac8d0b83b",
            `da_output/de_negmode_lipids_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_lipids_groupKO-groupWT_long_annot.tsv` = "568747756defd94e3305d975ce8b0de78947e9dd8bacc0c68bad0c67b77e2ca0",
            `da_output/de_posmode_lipids_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `feature_qc/itsd_lcms_pos.png` = "<present>",
            `feature_qc/lcms_neg_ruv_corrected_rle.png` = "<present>",
            `feature_qc/lcms_pos_pre_norm_pca.png` = "<present>",
            `feature_qc/lipidomics_composite_QC_figure.pdf` = "<present>",
            `feature_qc/lipidomics_composite_QC_figure.png` = "<present>",
            `feature_qc/ruv_normalised_results.RDS` = "<present>",
            `feature_qc/ruv_normalised_results.tsv` = "ab2c445d1a3d8f5ae94c31cd216f8330033e9b029afa2b6da7bbacfaf9743984",
            `publication_graphs/Heatmap/artifact.txt` = "<present>",
            `publication_graphs/Heatmaps/artifact.txt` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/lipids_num_sig_da_molecules.tab` = "ab123422c31c34ada4ea3ab9a0c9149177c0b9cdb09acdb41b16e9d4adfeed97",
            `publication_graphs/NumSigDeMolecules/artifact.txt` = "<present>",
            `publication_graphs/NumSigDeMolecules/lipids_num_sig_de_molecules.tab` = "ab123422c31c34ada4ea3ab9a0c9149177c0b9cdb09acdb41b16e9d4adfeed97",
            `publication_graphs/Volcano_Plots/artifact.txt` = "<present>",
            `source/contrasts_tbl.tab` = "7209220541b9240edc6b7b73030b7384687b68cd7814091847ac8baf36891f72",
            `source/design_matrix.tab` = "c3435c79da809286aebfaa17f024efdea7282282a0904b2f5f1477f4790b6983",
            `source/study_parameters.txt` = "<present>",
            `time/12_correlation_filtered_combined_plots.png` = "<present>"),
        study_parameters = "1 lines 55c1b3aadde3f3b7",
        integration = stats::setNames(character(0), character()),
        session_state = "<0 files>",
        errors = "<unhandled> unused arguments (self, name)"),
    save = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Study parameters created for: kit108 \nDescription: Full lipidomics analysis workflow with normalization and DA \nTimestamp: <TIME> \nFile: <ROOT>/source/study_parameters.txt \nSource: Final S4 object @args + config_list \nIntegration Object: lipidomics_kit108_final_s4.RDS \nStatus: Parameters saved [OK]",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = c(
            "[message 5 -] Saved Integration S4 Object", "[message 5 -] Study parameters saved successfully"),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = c(
            "createWorkflowArgsFromConfig(workflow_name=kit108, description=Full lipidomics analysis workflow with normalization and DA, organism_name=<NULL>, taxon_id=<NULL>, source_dir_path=<ROOT>/source, final_s4_object=<LipidomicsAssayData 5093d18fc8ce6eca>, contrasts_tbl=<table 2x3 4323c37e5b4e992a>, workflow_data=<environment>)",
            "writeWorkflowAnalysisRecord(manager=<WorkflowStateOwner>, path=<ROOT>/source/analysis_record.json)"),
        global = c(before = "absent", after = "absent"),
        files = c(
            `integration/lipidomics_kit108_final_s4.RDS` = "<present>",
            `source/analysis_record.json` = "<present>"),
        study_parameters = "158 lines bcb087f5fa824d3d",
        integration = c(lipidomics_kit108_final_s4.RDS = "<RDS> LipidomicsAssayData 5093d18fc8ce6eca"),
        session_state = "<unchanged>",
        errors = character(0)),
    copy = list(
        outputs = c("copy_status", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full lipidomics analysis workflow with normalization and DA \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK]",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = "[message 5 -] Publication files copied",
        updates = character(0),
        progress = "start Copying files to publication directory... 0.1",
        reactive_values = c(
            files_copied = "TRUE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = c(
            "copyToResultsSummary(omic_type=lipidomics, experiment_label=kit108, contrasts_tbl=<table 2x3 4323c37e5b4e992a>, design_matrix=<table 6x5 be07121076802bfa>, force=TRUE, current_rmd=<NULL>, project_dirs_object_name=project_dirs, global=the session's)",
            "  copyToResultsSummary returned: "),
        global = c(before = "absent", after = "the session's"),
        files = c(
            `results_summary/Publication_figures/Heatmaps/artifact.txt` = "<present>",
            `results_summary/Publication_figures/NumSigDeMolecules/artifact.txt` = "<present>",
            `results_summary/Publication_figures/NumSigDeMolecules/lipids_num_sig_de_molecules.tab` = "ab123422c31c34ada4ea3ab9a0c9149177c0b9cdb09acdb41b16e9d4adfeed97",
            `results_summary/Publication_tables/DA_results_lipidomics.xlsx` = "<present>",
            `results_summary/Publication_tables/Pathway_enrichment_results_lipidomics.xlsx` = "<present>",
            `results_summary/Publication_tables/RUV_normalised_results.tsv` = "ab2c445d1a3d8f5ae94c31cd216f8330033e9b029afa2b6da7bbacfaf9743984",
            `results_summary/Publication_tables/da_lipidomics_num_sig_da_molecules.tab` = "ab123422c31c34ada4ea3ab9a0c9149177c0b9cdb09acdb41b16e9d4adfeed97",
            `results_summary/Publication_tables/da_lipidomics_num_sig_de_molecules.tab` = "ab123422c31c34ada4ea3ab9a0c9149177c0b9cdb09acdb41b16e9d4adfeed97",
            `results_summary/Publication_tables/ruv_normalised_results.RDS` = "<present>",
            `results_summary/QC_figures/12_correlation_filtered_combined_plots.png` = "<present>",
            `results_summary/QC_figures/itsd_lcms_pos.png` = "<present>",
            `results_summary/QC_figures/lcms_neg_ruv_corrected_rle.png` = "<present>",
            `results_summary/QC_figures/lcms_pos_pre_norm_pca.png` = "<present>",
            `results_summary/QC_figures/lipidomics_composite_QC_figure.pdf` = "<present>",
            `results_summary/QC_figures/lipidomics_composite_QC_figure.png` = "<present>",
            `results_summary/Study_report/contrasts_tbl.tab` = "7209220541b9240edc6b7b73030b7384687b68cd7814091847ac8baf36891f72",
            `results_summary/Study_report/design_matrix.tab` = "c3435c79da809286aebfaa17f024efdea7282282a0904b2f5f1477f4790b6983",
            `results_summary/Study_report/study_parameters.txt` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    report = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full lipidomics analysis workflow with normalization and DA \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK] \nReport location: <ROOT>/results_summary/lipidomics_report_lipidomics_kit108.docx",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = c(
            "[message 5 -] lipidomics_report.rmd template copied from package",
            "[message 5 -] Report generated successfully!"),
        updates = character(0),
        progress = c(
            "start Generating report... 0.1", "0.2 Checking for lipidomics_report.rmd template...",
            "0.5 Rendering report...", "0.9 Exporting the analysis record..."),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/lipidomics_report_lipidomics_kit108.docx",
            workflow_args_saved = "TRUE"),
        calls = c(
            "RenderReport(omic_type=lipidomics, experiment_label=kit108, rmd_filename=lipidomics_report.rmd, project_dirs_object_name=project_dirs, output_format=<NULL>, global=the session's)",
            "exportWorkflowAnalysisRecordBundle(manager=<WorkflowStateOwner>, directory=<ROOT>/results_summary)",
            "writeWorkflowAnalysisRecord(manager=<WorkflowStateOwner>, path=<ROOT>/results_summary/analysis_record.json)"),
        global = c(before = "the session's", after = "the session's"),
        files = c(
            `results_summary/analysis_record.docx` = "<present>",
            `results_summary/analysis_record.json` = "<present>",
            `results_summary/analysis_record.xlsx` = "<present>",
            `results_summary/lipidomics_report_lipidomics_kit108.docx` = "<present>",
            `scripts/lipidomics/lipidomics_report.rmd` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = "<unhandled> unused arguments (self, name)"),
    export = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full lipidomics analysis workflow with normalization and DA \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK] \nReport location: <ROOT>/results_summary/lipidomics_report_lipidomics_kit108.docx",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = "[message 5 -] Session state exported to: <ROOT>/source/session_state_<DATE>.RDS",
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/lipidomics_report_lipidomics_kit108.docx",
            workflow_args_saved = "TRUE"),
        calls = character(0),
        global = c(before = "the session's", after = "the session's"),
        files = c(`source/session_state_<DATE>.RDS` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = c(
            experiment_label = "kit108",
            description = "Full lipidomics analysis workflow with normalization and DA", timestamp = "<time>",
            omic_type = "lipidomics", workflow_args_saved = "TRUE", files_copied = "TRUE",
            report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/lipidomics_report_lipidomics_kit108.docx",
            project_dirs = "<digest> 7c59d96e1fa7a5b8", workflow_type = "lipidomics",
            report_template = "lipidomics_report.rmd", parameter_payload = "<digest> 93bc513813bbc3b2",
            assay_names = "LCMS_Pos|LCMS_Neg", feature_counts = "6|6", sample_count = "6",
            contrast_count = "2"),
        errors = character(0)),
    github_disabled = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full lipidomics analysis workflow with normalization and DA \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK] \nReport location: <ROOT>/results_summary/lipidomics_report_lipidomics_kit108.docx",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/lipidomics_report_lipidomics_kit108.docx",
            workflow_args_saved = "TRUE"),
        calls = character(0),
        global = c(before = "the session's", after = "the session's"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    github = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full lipidomics analysis workflow with normalization and DA \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK], GitHub pushed [OK]",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = "[message 5 -] Successfully pushed to GitHub",
        updates = character(0),
        progress = "start Pushing to GitHub... 0.1",
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/lipidomics_report_lipidomics_kit108.docx",
            workflow_args_saved = "TRUE"),
        calls = "pushProjectToGithubFromDirs(project_dirs=the session's, omic_type=lipidomics, experiment_label=kit108, project_id=kit108-project, github_org=kit108-org, github_user_email=kit108@example.org, github_user_name=kit108-user)",
        global = c(before = "the session's", after = "the session's"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)))

omicSummaryW1Frozen[["lipid_variants"]] <- list(
    start = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "Template: Lipidomics Report [OK]"),
        notes = character(0),
        updates = "summary-experiment_label: <same> -> kit108",
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "the session's", after = "the session's"),
        files = c(
            `da_output/de_negmode_lipids_groupKO-groupWT_long_annot.tsv` = "1da4796cae289522efd567b881e22cd01afe80ebea6b5caf33a7b10ac8d0b83b",
            `da_output/de_negmode_lipids_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_lipids_groupKO-groupWT_long_annot.tsv` = "568747756defd94e3305d975ce8b0de78947e9dd8bacc0c68bad0c67b77e2ca0",
            `da_output/de_posmode_lipids_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `feature_qc/ruv_normalised_results.RDS` = "<present>",
            `feature_qc/ruv_normalised_results.tsv` = "ab2c445d1a3d8f5ae94c31cd216f8330033e9b029afa2b6da7bbacfaf9743984",
            `scripts/lipidomics/lipidomics_report.rmd` = "<present>",
            `source/contrasts_tbl.tab` = "7209220541b9240edc6b7b73030b7384687b68cd7814091847ac8baf36891f72",
            `source/design_matrix.tab` = "c3435c79da809286aebfaa17f024efdea7282282a0904b2f5f1477f4790b6983",
            `source/study_parameters.txt` = "<present>"),
        study_parameters = "1 lines 55c1b3aadde3f3b7",
        integration = stats::setNames(character(0), character()),
        session_state = "<0 files>",
        errors = "<unhandled> unused arguments (self, name)"),
    save_record_fails = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Study parameters created for: kit108 \nDescription: Full lipidomics analysis workflow with normalization and DA \nTimestamp: <TIME> \nFile: <ROOT>/source/study_parameters.txt \nSource: Final S4 object @args + config_list \nIntegration Object: lipidomics_kit108_final_s4.RDS \nStatus: Parameters saved [OK]",
            template_status = "Template: Lipidomics Report [OK]"),
        notes = c(
            "[warning NULL -] The analysis record (analysis_record.json) was not written: kit108 record double fails",
            "[message 5 -] Saved Integration S4 Object", "[message 5 -] Study parameters saved successfully"),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = c(
            "createWorkflowArgsFromConfig(workflow_name=kit108, description=Full lipidomics analysis workflow with normalization and DA, organism_name=<NULL>, taxon_id=<NULL>, source_dir_path=<ROOT>/source, final_s4_object=<LipidomicsAssayData 5093d18fc8ce6eca>, contrasts_tbl=<table 2x3 4323c37e5b4e992a>, workflow_data=<environment>)",
            "writeWorkflowAnalysisRecord(manager=<WorkflowStateOwner>, path=<ROOT>/source/analysis_record.json)"),
        global = c(before = "the session's", after = "the session's"),
        files = c(`integration/lipidomics_kit108_final_s4.RDS` = "<present>"),
        study_parameters = "158 lines bcb087f5fa824d3d",
        integration = c(lipidomics_kit108_final_s4.RDS = "<RDS> LipidomicsAssayData 5093d18fc8ce6eca"),
        session_state = "<unchanged>",
        errors = character(0)),
    copy = list(
        outputs = c("copy_status", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory with optional missing artifacts [WARNING]",
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full lipidomics analysis workflow with normalization and DA \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK]",
            template_status = "Template: Lipidomics Report [OK]"),
        notes = "[warning 5 -] Publication files copied with optional missing artifacts",
        updates = character(0),
        progress = "start Copying files to publication directory... 0.1",
        reactive_values = c(
            files_copied = "TRUE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = c(
            "copyToResultsSummary(omic_type=lipidomics, experiment_label=kit108, contrasts_tbl=<table 2x3 4323c37e5b4e992a>, design_matrix=<table 6x5 be07121076802bfa>, force=TRUE, current_rmd=<NULL>, project_dirs_object_name=project_dirs, global=the session's)",
            "  copyToResultsSummary returned: DA Heatmaps | Num Sig DE Molecules"),
        global = c(before = "the session's", after = "the session's"),
        files = c(
            `results_summary/Publication_tables/DA_results_lipidomics.xlsx` = "<present>",
            `results_summary/Publication_tables/Pathway_enrichment_results_lipidomics.xlsx` = "<present>",
            `results_summary/Publication_tables/RUV_normalised_results.tsv` = "ab2c445d1a3d8f5ae94c31cd216f8330033e9b029afa2b6da7bbacfaf9743984",
            `results_summary/Publication_tables/ruv_normalised_results.RDS` = "<present>",
            `results_summary/Study_report/contrasts_tbl.tab` = "7209220541b9240edc6b7b73030b7384687b68cd7814091847ac8baf36891f72",
            `results_summary/Study_report/design_matrix.tab` = "c3435c79da809286aebfaa17f024efdea7282282a0904b2f5f1477f4790b6983",
            `results_summary/Study_report/study_parameters.txt` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    report_bundle_fails = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory with optional missing artifacts [WARNING]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full lipidomics analysis workflow with normalization and DA \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK] \nReport location: <ROOT>/results_summary/lipidomics_report_lipidomics_kit108.docx",
            template_status = "Template: Lipidomics Report [OK]"),
        notes = c(
            "[warning NULL -] The analysis record bundle was not exported: kit108 bundle double fails",
            "[message 5 -] Report generated successfully!"),
        updates = character(0),
        progress = c(
            "start Generating report... 0.1", "0.5 Rendering report...",
            "0.9 Exporting the analysis record..."),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/lipidomics_report_lipidomics_kit108.docx",
            workflow_args_saved = "TRUE"),
        calls = c(
            "RenderReport(omic_type=lipidomics, experiment_label=kit108, rmd_filename=lipidomics_report.rmd, project_dirs_object_name=project_dirs, output_format=<NULL>, global=the session's)",
            "exportWorkflowAnalysisRecordBundle(manager=<WorkflowStateOwner>, directory=<ROOT>/results_summary)"),
        global = c(before = "the session's", after = "the session's"),
        files = c(`results_summary/lipidomics_report_lipidomics_kit108.docx` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = "<unhandled> unused arguments (self, name)"))

omicSummaryW1Frozen[["lipid_required_missing"]] <- list(
    start = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = character(0),
        updates = "summary-experiment_label: <same> -> kit108",
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = character(0),
        study_parameters = "<absent>",
        integration = stats::setNames(character(0), character()),
        session_state = "<0 files>",
        errors = "<unhandled> unused arguments (self, name)"),
    copy = list(
        outputs = c("copy_status", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Error: Required publication artifacts missing: Contrasts Table, Design Matrix",
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = "[error 10 -] Copy error: Required publication artifacts missing: Contrasts Table, Design Matrix",
        updates = character(0),
        progress = "start Copying files to publication directory... 0.1",
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = c(
            "copyToResultsSummary(omic_type=lipidomics, experiment_label=kit108, contrasts_tbl=<NULL>, design_matrix=<NULL>, force=TRUE, current_rmd=<NULL>, project_dirs_object_name=project_dirs, global=the session's)",
            "  copyToResultsSummary returned: Contrasts Table | Design Matrix | DA Heatmaps | Num Sig DE Molecules"),
        global = c(before = "absent", after = "the session's"),
        files = c(
            `results_summary/Publication_tables/DA_results_lipidomics.xlsx` = "<present>",
            `results_summary/Publication_tables/Pathway_enrichment_results_lipidomics.xlsx` = "<present>",
            `results_summary/Study_report/study_parameters.txt` = "<present>",
            `source/study_parameters.txt` = "<present>"),
        study_parameters = "7 lines 54f68b59e9e2c5dc",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    report = list(
        outputs = c("copy_status", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Error: Required publication artifacts missing: Contrasts Table, Design Matrix",
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = "[warning 5 -] Please copy the results to the publication directory before generating the report (Copy to Publication Directory).",
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = character(0),
        global = c(before = "the session's", after = "the session's"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)))

omicSummaryW1Frozen[["lipid_github_fallback"]] <- list(
    copy = list(
        outputs = c("copy_status", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full lipidomics analysis workflow with normalization and DA \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK]",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = "[message 5 -] Publication files copied",
        updates = "summary-experiment_label: <same> -> kit108",
        progress = "start Copying files to publication directory... 0.1",
        reactive_values = c(
            files_copied = "TRUE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = c(
            "copyToResultsSummary(omic_type=lipidomics, experiment_label=kit108, contrasts_tbl=<table 2x3 4323c37e5b4e992a>, design_matrix=<table 6x5 be07121076802bfa>, force=TRUE, current_rmd=<NULL>, project_dirs_object_name=project_dirs, global=the session's)",
            "  copyToResultsSummary returned: "),
        global = c(before = "the session's", after = "the session's"),
        files = c(
            `da_output/de_negmode_lipids_groupKO-groupWT_long_annot.tsv` = "1da4796cae289522efd567b881e22cd01afe80ebea6b5caf33a7b10ac8d0b83b",
            `da_output/de_negmode_lipids_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_lipids_groupKO-groupWT_long_annot.tsv` = "568747756defd94e3305d975ce8b0de78947e9dd8bacc0c68bad0c67b77e2ca0",
            `da_output/de_posmode_lipids_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `feature_qc/itsd_lcms_pos.png` = "<present>",
            `feature_qc/lcms_neg_ruv_corrected_rle.png` = "<present>",
            `feature_qc/lcms_pos_pre_norm_pca.png` = "<present>",
            `feature_qc/lipidomics_composite_QC_figure.pdf` = "<present>",
            `feature_qc/lipidomics_composite_QC_figure.png` = "<present>",
            `feature_qc/ruv_normalised_results.RDS` = "<present>",
            `feature_qc/ruv_normalised_results.tsv` = "ab2c445d1a3d8f5ae94c31cd216f8330033e9b029afa2b6da7bbacfaf9743984",
            `publication_graphs/Heatmap/artifact.txt` = "<present>",
            `publication_graphs/Heatmaps/artifact.txt` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/lipids_num_sig_da_molecules.tab` = "ab123422c31c34ada4ea3ab9a0c9149177c0b9cdb09acdb41b16e9d4adfeed97",
            `publication_graphs/NumSigDeMolecules/artifact.txt` = "<present>",
            `publication_graphs/NumSigDeMolecules/lipids_num_sig_de_molecules.tab` = "ab123422c31c34ada4ea3ab9a0c9149177c0b9cdb09acdb41b16e9d4adfeed97",
            `publication_graphs/Volcano_Plots/artifact.txt` = "<present>",
            `results_summary/Publication_figures/Heatmaps/artifact.txt` = "<present>",
            `results_summary/Publication_figures/NumSigDeMolecules/artifact.txt` = "<present>",
            `results_summary/Publication_figures/NumSigDeMolecules/lipids_num_sig_de_molecules.tab` = "ab123422c31c34ada4ea3ab9a0c9149177c0b9cdb09acdb41b16e9d4adfeed97",
            `results_summary/Publication_tables/DA_results_lipidomics.xlsx` = "<present>",
            `results_summary/Publication_tables/Pathway_enrichment_results_lipidomics.xlsx` = "<present>",
            `results_summary/Publication_tables/RUV_normalised_results.tsv` = "ab2c445d1a3d8f5ae94c31cd216f8330033e9b029afa2b6da7bbacfaf9743984",
            `results_summary/Publication_tables/da_lipidomics_num_sig_da_molecules.tab` = "ab123422c31c34ada4ea3ab9a0c9149177c0b9cdb09acdb41b16e9d4adfeed97",
            `results_summary/Publication_tables/da_lipidomics_num_sig_de_molecules.tab` = "ab123422c31c34ada4ea3ab9a0c9149177c0b9cdb09acdb41b16e9d4adfeed97",
            `results_summary/Publication_tables/ruv_normalised_results.RDS` = "<present>",
            `results_summary/QC_figures/12_correlation_filtered_combined_plots.png` = "<present>",
            `results_summary/QC_figures/itsd_lcms_pos.png` = "<present>",
            `results_summary/QC_figures/lcms_neg_ruv_corrected_rle.png` = "<present>",
            `results_summary/QC_figures/lcms_pos_pre_norm_pca.png` = "<present>",
            `results_summary/QC_figures/lipidomics_composite_QC_figure.pdf` = "<present>",
            `results_summary/QC_figures/lipidomics_composite_QC_figure.png` = "<present>",
            `results_summary/Study_report/contrasts_tbl.tab` = "7209220541b9240edc6b7b73030b7384687b68cd7814091847ac8baf36891f72",
            `results_summary/Study_report/design_matrix.tab` = "c3435c79da809286aebfaa17f024efdea7282282a0904b2f5f1477f4790b6983",
            `results_summary/Study_report/study_parameters.txt` = "<present>",
            `source/contrasts_tbl.tab` = "7209220541b9240edc6b7b73030b7384687b68cd7814091847ac8baf36891f72",
            `source/design_matrix.tab` = "c3435c79da809286aebfaa17f024efdea7282282a0904b2f5f1477f4790b6983",
            `source/study_parameters.txt` = "<present>",
            `time/12_correlation_filtered_combined_plots.png` = "<present>"),
        study_parameters = "1 lines 55c1b3aadde3f3b7",
        integration = stats::setNames(character(0), character()),
        session_state = "<0 files>",
        errors = "<unhandled> unused arguments (self, name)"),
    report_fallback = list(
        outputs = c("copy_status", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full lipidomics analysis workflow with normalization and DA \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK]",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = "[error 10 -] Template retrieval failed: kit108 download double: offline",
        updates = character(0),
        progress = c(
            "start Generating report... 0.1", "0.2 Checking for lipidomics_report.rmd template..."),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = c(
            "system.file(path=reports/lipidomics/lipidomics_report.rmd, package=MultiScholaR)",
            "download.file(url=https://raw.githubusercontent.com/APAF-bioinformatics/MultiScholaR/main/inst/reports/lipidomics/lipidomics_report.rmd, destfile=<ROOT>/scripts/lipidomics/lipidomics_report.rmd, quiet=TRUE)"),
        global = c(before = "the session's", after = "the session's"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)))

omicSummaryW1Frozen[["lipid_managed"]] <- list(
    start = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = character(0),
        updates = "summary-experiment_label: <same> -> kit108",
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = c(
            `da_output/de_negmode_lipids_groupKO-groupWT_long_annot.tsv` = "1da4796cae289522efd567b881e22cd01afe80ebea6b5caf33a7b10ac8d0b83b",
            `da_output/de_negmode_lipids_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_lipids_groupKO-groupWT_long_annot.tsv` = "568747756defd94e3305d975ce8b0de78947e9dd8bacc0c68bad0c67b77e2ca0",
            `da_output/de_posmode_lipids_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `feature_qc/itsd_lcms_pos.png` = "<present>",
            `feature_qc/lcms_neg_ruv_corrected_rle.png` = "<present>",
            `feature_qc/lcms_pos_pre_norm_pca.png` = "<present>",
            `feature_qc/lipidomics_composite_QC_figure.pdf` = "<present>",
            `feature_qc/lipidomics_composite_QC_figure.png` = "<present>",
            `feature_qc/ruv_normalised_results.RDS` = "<present>",
            `feature_qc/ruv_normalised_results.tsv` = "ab2c445d1a3d8f5ae94c31cd216f8330033e9b029afa2b6da7bbacfaf9743984",
            `publication_graphs/Heatmap/artifact.txt` = "<present>",
            `publication_graphs/Heatmaps/artifact.txt` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/lipids_num_sig_da_molecules.tab` = "ab123422c31c34ada4ea3ab9a0c9149177c0b9cdb09acdb41b16e9d4adfeed97",
            `publication_graphs/NumSigDeMolecules/artifact.txt` = "<present>",
            `publication_graphs/NumSigDeMolecules/lipids_num_sig_de_molecules.tab` = "ab123422c31c34ada4ea3ab9a0c9149177c0b9cdb09acdb41b16e9d4adfeed97",
            `publication_graphs/Volcano_Plots/artifact.txt` = "<present>",
            `source/contrasts_tbl.tab` = "7209220541b9240edc6b7b73030b7384687b68cd7814091847ac8baf36891f72",
            `source/design_matrix.tab` = "c3435c79da809286aebfaa17f024efdea7282282a0904b2f5f1477f4790b6983",
            `source/study_parameters.txt` = "<present>",
            `time/12_correlation_filtered_combined_plots.png` = "<present>"),
        study_parameters = "1 lines 55c1b3aadde3f3b7",
        integration = stats::setNames(character(0), character()),
        session_state = "<0 files>",
        errors = "<unhandled> unused arguments (self, name)"),
    save = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Study parameters created for: kit108 \nDescription: Full lipidomics analysis workflow with normalization and DA \nTimestamp: <TIME> \nFile: <ROOT>/source/study_parameters.txt \nSource: Final S4 object @args + config_list \nIntegration Object: lipidomics_kit108_final_s4.RDS \nStatus: Parameters saved [OK]",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = c(
            "[default NULL lipidomics_action_working] Applying lipidomics step ...",
            "remove lipidomics_action_working", "[message 5 -] Saved Integration S4 Object",
            "[message 5 -] Study parameters saved successfully"),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = c(
            "createWorkflowArgsFromConfig(workflow_name=kit108, description=Full lipidomics analysis workflow with normalization and DA, organism_name=<NULL>, taxon_id=<NULL>, source_dir_path=<ROOT>/source, final_s4_object=<LipidomicsAssayData 0c6c2c4042cdbee7>, contrasts_tbl=<table 2x3 4323c37e5b4e992a>, workflow_data=<environment>)",
            "writeWorkflowAnalysisRecord(manager=<WorkflowStateOwner>, path=<ROOT>/source/analysis_record.json)",
            "omicsManagedProduct(operation=lipid_integration_export, parameters=path:<ROOT>/integration/lipidomics_kit108_final_s4.RDS;state_name:lipid_correlation_filtered)"),
        global = c(before = "absent", after = "absent"),
        files = c(
            `integration/lipidomics_kit108_final_s4.RDS` = "<present>",
            `source/analysis_record.json` = "<present>"),
        study_parameters = "158 lines bcb087f5fa824d3d",
        integration = c(lipidomics_kit108_final_s4.RDS = "<RDS> LipidomicsAssayData 5093d18fc8ce6eca"),
        session_state = "<unchanged>",
        errors = character(0)),
    export = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Study parameters created for: kit108 \nDescription: Full lipidomics analysis workflow with normalization and DA \nTimestamp: <TIME> \nFile: <ROOT>/source/study_parameters.txt \nSource: Final S4 object @args + config_list \nIntegration Object: lipidomics_kit108_final_s4.RDS \nStatus: Parameters saved [OK]",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = "[message 5 -] Session state exported to: <ROOT>/source/session_state_<DATE>.RDS",
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = c(`source/session_state_<DATE>.RDS` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = c(
            experiment_label = "kit108",
            description = "Full lipidomics analysis workflow with normalization and DA", timestamp = "<time>",
            omic_type = "lipidomics", workflow_args_saved = "TRUE", files_copied = "FALSE",
            report_generated = "FALSE", report_path = "<NULL>", project_dirs = "<digest> 7c59d96e1fa7a5b8",
            workflow_type = "lipidomics", report_template = "lipidomics_report.rmd",
            parameter_payload = "<digest> 93bc513813bbc3b2", assay_names = "LCMS_Pos|LCMS_Neg",
            feature_counts = "6|6", sample_count = "6", contrast_count = "2"),
        errors = character(0)))

omicSummaryW1Frozen[["lipid_global_another"]] <- list(
    copy = list(
        outputs = c("copy_status", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Error: Failed to get project paths: Could not find project paths for omic_type='lipidomics', experiment_label='kit108'.\nAvailable keys in 'project_dirs': elsewhere\nEnsure setupDirectories() was called with the correct omic_type and label.",
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = "[error 10 -] Copy error: Failed to get project paths: Could not find project paths for omic_type='lipidomics', experiment_label='kit108'.\nAvailable keys in 'project_dirs': elsewhere\nEnsure setupDirectories() was called with the correct omic_type and label.",
        updates = "summary-experiment_label: <same> -> kit108",
        progress = "start Copying files to publication directory... 0.1",
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = "copyToResultsSummary(omic_type=lipidomics, experiment_label=kit108, contrasts_tbl=<table 2x3 4323c37e5b4e992a>, design_matrix=<table 6x5 be07121076802bfa>, force=TRUE, current_rmd=<NULL>, project_dirs_object_name=project_dirs, global=another)",
        global = c(before = "another", after = "another"),
        files = c(
            `da_output/de_negmode_lipids_groupKO-groupWT_long_annot.tsv` = "1da4796cae289522efd567b881e22cd01afe80ebea6b5caf33a7b10ac8d0b83b",
            `da_output/de_negmode_lipids_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_lipids_groupKO-groupWT_long_annot.tsv` = "568747756defd94e3305d975ce8b0de78947e9dd8bacc0c68bad0c67b77e2ca0",
            `da_output/de_posmode_lipids_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `feature_qc/itsd_lcms_pos.png` = "<present>",
            `feature_qc/lcms_neg_ruv_corrected_rle.png` = "<present>",
            `feature_qc/lcms_pos_pre_norm_pca.png` = "<present>",
            `feature_qc/lipidomics_composite_QC_figure.pdf` = "<present>",
            `feature_qc/lipidomics_composite_QC_figure.png` = "<present>",
            `feature_qc/ruv_normalised_results.RDS` = "<present>",
            `feature_qc/ruv_normalised_results.tsv` = "ab2c445d1a3d8f5ae94c31cd216f8330033e9b029afa2b6da7bbacfaf9743984",
            `publication_graphs/Heatmap/artifact.txt` = "<present>",
            `publication_graphs/Heatmaps/artifact.txt` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/lipids_num_sig_da_molecules.tab` = "ab123422c31c34ada4ea3ab9a0c9149177c0b9cdb09acdb41b16e9d4adfeed97",
            `publication_graphs/NumSigDeMolecules/artifact.txt` = "<present>",
            `publication_graphs/NumSigDeMolecules/lipids_num_sig_de_molecules.tab` = "ab123422c31c34ada4ea3ab9a0c9149177c0b9cdb09acdb41b16e9d4adfeed97",
            `publication_graphs/Volcano_Plots/artifact.txt` = "<present>",
            `source/contrasts_tbl.tab` = "7209220541b9240edc6b7b73030b7384687b68cd7814091847ac8baf36891f72",
            `source/design_matrix.tab` = "c3435c79da809286aebfaa17f024efdea7282282a0904b2f5f1477f4790b6983",
            `source/study_parameters.txt` = "<present>",
            `time/12_correlation_filtered_combined_plots.png` = "<present>"),
        study_parameters = "1 lines 55c1b3aadde3f3b7",
        integration = stats::setNames(character(0), character()),
        session_state = "<0 files>",
        errors = "<unhandled> unused arguments (self, name)"))

omicSummaryW1Frozen[["metab_refusals"]] <- list(
    start = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = c(
            `da_output/de_gcms_metabolites_groupKO-groupWT_long_annot.tsv` = "13be76abf9da00e9d87f6e6f3e094de1961db1fdcb37f7611de9efda60992c01",
            `da_output/de_gcms_metabolites_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_metabolites_groupKO-groupWT_long_annot.tsv` = "3fe65bf0975381006239ba8531060b86c64ba5d9279f804722b1938d6cc28381",
            `da_output/de_posmode_metabolites_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `feature_qc/gcms_post_norm_density.png` = "<present>",
            `feature_qc/itsd_lcms_pos.png` = "<present>",
            `feature_qc/lcms_pos_pre_norm_pca.png` = "<present>",
            `feature_qc/lcms_pos_ruv_corrected_rle.png` = "<present>",
            `feature_qc/metabolomics_composite_QC_figure.pdf` = "<present>",
            `feature_qc/metabolomics_composite_QC_figure.png` = "<present>",
            `feature_qc/ruv_normalised_results.RDS` = "<present>",
            `feature_qc/ruv_normalised_results.tsv` = "f87484c97f41f65f89307e33c596dfaf95499e9ce230ff4f79764c6bb795a2dd",
            `pathway_enrichment/KO_vs_WT_down_enrichment_results.tsv` = "278827ebf4745b50520d3714c4389fd9c0a3c817e116638b4227257b9ed2a07a",
            `pathway_enrichment/KO_vs_WT_up_enrichment_results.tsv` = "278827ebf4745b50520d3714c4389fd9c0a3c817e116638b4227257b9ed2a07a",
            `publication_graphs/Heatmap/artifact.txt` = "<present>",
            `publication_graphs/Heatmaps/artifact.txt` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/metabolites_num_sig_da_molecules.tab` = "32c301d1f9145915cfa6a4789f1c5c7a9b8f69ac15917d31f28bb317eab608c9",
            `publication_graphs/Volcano_Plots/artifact.txt` = "<present>",
            `source/contrasts_tbl.tab` = "7209220541b9240edc6b7b73030b7384687b68cd7814091847ac8baf36891f72",
            `source/design_matrix.tab` = "c3435c79da809286aebfaa17f024efdea7282282a0904b2f5f1477f4790b6983",
            `source/study_parameters.txt` = "<present>",
            `time/12_correlation_filtered_combined_plots.png` = "<present>"),
        study_parameters = "1 lines 55c1b3aadde3f3b7",
        integration = stats::setNames(character(0), character()),
        session_state = "<0 files>",
        errors = "<unhandled> unused arguments (self, name)"),
    save_no_label = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    copy_no_label = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    report_no_label = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = "[warning 5 -] Please enter an experiment label before generating the report.",
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    export_no_label = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    report_before_copy = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = "[warning 5 -] Please copy the results to the publication directory before generating the report (Copy to Publication Directory).",
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    github_before_report = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)))

omicSummaryW1Frozen[["metab_main"]] <- list(
    start = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = character(0),
        updates = "summary-experiment_label: <same> -> kit108",
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = c(
            `da_output/de_gcms_metabolites_groupKO-groupWT_long_annot.tsv` = "13be76abf9da00e9d87f6e6f3e094de1961db1fdcb37f7611de9efda60992c01",
            `da_output/de_gcms_metabolites_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_metabolites_groupKO-groupWT_long_annot.tsv` = "3fe65bf0975381006239ba8531060b86c64ba5d9279f804722b1938d6cc28381",
            `da_output/de_posmode_metabolites_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `feature_qc/gcms_post_norm_density.png` = "<present>",
            `feature_qc/itsd_lcms_pos.png` = "<present>",
            `feature_qc/lcms_pos_pre_norm_pca.png` = "<present>",
            `feature_qc/lcms_pos_ruv_corrected_rle.png` = "<present>",
            `feature_qc/metabolomics_composite_QC_figure.pdf` = "<present>",
            `feature_qc/metabolomics_composite_QC_figure.png` = "<present>",
            `feature_qc/ruv_normalised_results.RDS` = "<present>",
            `feature_qc/ruv_normalised_results.tsv` = "f87484c97f41f65f89307e33c596dfaf95499e9ce230ff4f79764c6bb795a2dd",
            `pathway_enrichment/KO_vs_WT_down_enrichment_results.tsv` = "278827ebf4745b50520d3714c4389fd9c0a3c817e116638b4227257b9ed2a07a",
            `pathway_enrichment/KO_vs_WT_up_enrichment_results.tsv` = "278827ebf4745b50520d3714c4389fd9c0a3c817e116638b4227257b9ed2a07a",
            `publication_graphs/Heatmap/artifact.txt` = "<present>",
            `publication_graphs/Heatmaps/artifact.txt` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/metabolites_num_sig_da_molecules.tab` = "32c301d1f9145915cfa6a4789f1c5c7a9b8f69ac15917d31f28bb317eab608c9",
            `publication_graphs/Volcano_Plots/artifact.txt` = "<present>",
            `source/contrasts_tbl.tab` = "7209220541b9240edc6b7b73030b7384687b68cd7814091847ac8baf36891f72",
            `source/design_matrix.tab` = "c3435c79da809286aebfaa17f024efdea7282282a0904b2f5f1477f4790b6983",
            `source/study_parameters.txt` = "<present>",
            `time/12_correlation_filtered_combined_plots.png` = "<present>"),
        study_parameters = "1 lines 55c1b3aadde3f3b7",
        integration = stats::setNames(character(0), character()),
        session_state = "<0 files>",
        errors = "<unhandled> unused arguments (self, name)"),
    save = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Study parameters created for: kit108 \nDescription: Full metabolomics analysis workflow with normalization and DA \nTimestamp: <TIME> \nFile: <ROOT>/source/study_parameters.txt \nSource: Final S4 object @args + config_list \nIntegration Object: metabolomics_kit108_final_s4.RDS \nStatus: Parameters saved [OK]",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = c(
            "[message 5 -] Saved Integration S4 Object", "[message 5 -] Study parameters saved successfully"),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = c(
            "createWorkflowArgsFromConfig(workflow_name=kit108, description=Full metabolomics analysis workflow with normalization and DA, organism_name=<NULL>, taxon_id=<NULL>, source_dir_path=<ROOT>/source, final_s4_object=<MetaboliteAssayData 075de31ac869b80e>, contrasts_tbl=<table 2x3 4323c37e5b4e992a>, workflow_data=<environment>)",
            "writeWorkflowAnalysisRecord(manager=<WorkflowStateOwner>, path=<ROOT>/source/analysis_record.json)"),
        global = c(before = "absent", after = "absent"),
        files = c(
            `integration/metabolomics_kit108_final_s4.RDS` = "<present>",
            `source/analysis_record.json` = "<present>"),
        study_parameters = "127 lines 6f2dd5b056f531b5",
        integration = c(metabolomics_kit108_final_s4.RDS = "<RDS> MetaboliteAssayData 075de31ac869b80e"),
        session_state = "<unchanged>",
        errors = character(0)),
    copy = list(
        outputs = c("copy_status", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full metabolomics analysis workflow with normalization and DA \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK]",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = "[message 5 -] Publication files copied",
        updates = character(0),
        progress = "start Copying files to publication directory... 0.1",
        reactive_values = c(
            files_copied = "TRUE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = c(
            "copyToResultsSummary(omic_type=metabolomics, experiment_label=kit108, contrasts_tbl=<table 2x3 4323c37e5b4e992a>, design_matrix=<table 6x5 be07121076802bfa>, force=TRUE, current_rmd=<NULL>, project_dirs_object_name=project_dirs, global=the session's)",
            "  copyToResultsSummary returned: "),
        global = c(before = "absent", after = "the session's"),
        files = c(
            `results_summary/Publication_figures/Enrichment_Plots/KO_vs_WT_down_enrichment_results.tsv` = "278827ebf4745b50520d3714c4389fd9c0a3c817e116638b4227257b9ed2a07a",
            `results_summary/Publication_figures/Enrichment_Plots/KO_vs_WT_up_enrichment_results.tsv` = "278827ebf4745b50520d3714c4389fd9c0a3c817e116638b4227257b9ed2a07a",
            `results_summary/Publication_figures/Heatmaps/artifact.txt` = "<present>",
            `results_summary/Publication_tables/DA_results_metabolomics.xlsx` = "<present>",
            `results_summary/Publication_tables/Pathway_enrichment_results_metabolomics.xlsx` = "<present>",
            `results_summary/Publication_tables/RUV_normalised_results.tsv` = "f87484c97f41f65f89307e33c596dfaf95499e9ce230ff4f79764c6bb795a2dd",
            `results_summary/Publication_tables/da_metabolomics_num_sig_da_molecules.tab` = "32c301d1f9145915cfa6a4789f1c5c7a9b8f69ac15917d31f28bb317eab608c9",
            `results_summary/Publication_tables/ruv_normalised_results.RDS` = "<present>",
            `results_summary/QC_figures/12_correlation_filtered_combined_plots.png` = "<present>",
            `results_summary/QC_figures/gcms_post_norm_density.png` = "<present>",
            `results_summary/QC_figures/itsd_lcms_pos.png` = "<present>",
            `results_summary/QC_figures/lcms_pos_pre_norm_pca.png` = "<present>",
            `results_summary/QC_figures/lcms_pos_ruv_corrected_rle.png` = "<present>",
            `results_summary/QC_figures/metabolomics_composite_QC_figure.png` = "<present>",
            `results_summary/Study_report/contrasts_tbl.tab` = "7209220541b9240edc6b7b73030b7384687b68cd7814091847ac8baf36891f72",
            `results_summary/Study_report/design_matrix.tab` = "c3435c79da809286aebfaa17f024efdea7282282a0904b2f5f1477f4790b6983",
            `results_summary/Study_report/study_parameters.txt` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    report = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full metabolomics analysis workflow with normalization and DA \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK] \nReport location: <ROOT>/results_summary/metabolomics_report_metabolomics_kit108.docx",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = c(
            "[message 5 -] metabolomics_report.rmd template copied from package",
            "[message 5 -] Report generated successfully!"),
        updates = character(0),
        progress = c(
            "start Generating report... 0.1", "0.2 Checking for metabolomics_report.rmd template...",
            "0.5 Rendering report...", "0.9 Exporting the analysis record..."),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/metabolomics_report_metabolomics_kit108.docx",
            workflow_args_saved = "TRUE"),
        calls = c(
            "RenderReport(omic_type=metabolomics, experiment_label=kit108, rmd_filename=metabolomics_report.rmd, project_dirs_object_name=project_dirs, output_format=<NULL>, global=the session's)",
            "exportWorkflowAnalysisRecordBundle(manager=<WorkflowStateOwner>, directory=<ROOT>/results_summary)",
            "writeWorkflowAnalysisRecord(manager=<WorkflowStateOwner>, path=<ROOT>/results_summary/analysis_record.json)"),
        global = c(before = "the session's", after = "the session's"),
        files = c(
            `results_summary/analysis_record.docx` = "<present>",
            `results_summary/analysis_record.json` = "<present>",
            `results_summary/analysis_record.xlsx` = "<present>",
            `results_summary/metabolomics_report_metabolomics_kit108.docx` = "<present>",
            `scripts/metabolomics/metabolomics_report.rmd` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = "<unhandled> unused arguments (self, name)"),
    export = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full metabolomics analysis workflow with normalization and DA \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK] \nReport location: <ROOT>/results_summary/metabolomics_report_metabolomics_kit108.docx",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = "[message 5 -] Session state exported to: <ROOT>/source/session_state_<DATE>.RDS",
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/metabolomics_report_metabolomics_kit108.docx",
            workflow_args_saved = "TRUE"),
        calls = character(0),
        global = c(before = "the session's", after = "the session's"),
        files = c(`source/session_state_<DATE>.RDS` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = c(
            experiment_label = "kit108",
            description = "Full metabolomics analysis workflow with normalization and DA",
            timestamp = "<time>", omic_type = "metabolomics", workflow_args_saved = "TRUE",
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/metabolomics_report_metabolomics_kit108.docx",
            project_dirs = "<digest> 4701bf62c551a512", workflow_type = "metabolomics",
            report_template = "metabolomics_report.rmd", parameter_payload = "<digest> d66aa0b35bc9787b",
            assay_names = "LCMS_Pos|GCMS    ", feature_counts = "6|6", sample_count = "6",
            contrast_count = "2"),
        errors = character(0)),
    github_disabled = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full metabolomics analysis workflow with normalization and DA \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK] \nReport location: <ROOT>/results_summary/metabolomics_report_metabolomics_kit108.docx",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/metabolomics_report_metabolomics_kit108.docx",
            workflow_args_saved = "TRUE"),
        calls = character(0),
        global = c(before = "the session's", after = "the session's"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    github = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full metabolomics analysis workflow with normalization and DA \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK], GitHub pushed [OK]",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = "[message 5 -] Successfully pushed to GitHub",
        updates = character(0),
        progress = "start Pushing to GitHub... 0.1",
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/metabolomics_report_metabolomics_kit108.docx",
            workflow_args_saved = "TRUE"),
        calls = "pushProjectToGithubFromDirs(project_dirs=the session's, omic_type=metabolomics, experiment_label=kit108, project_id=kit108-project, github_org=kit108-org, github_user_email=kit108@example.org, github_user_name=kit108-user)",
        global = c(before = "the session's", after = "the session's"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)))

omicSummaryW1Frozen[["metab_variants"]] <- list(
    start = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "Template: Metabolomics Report [OK]"),
        notes = character(0),
        updates = "summary-experiment_label: <same> -> kit108",
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "the session's", after = "the session's"),
        files = c(
            `da_output/de_gcms_metabolites_groupKO-groupWT_long_annot.tsv` = "13be76abf9da00e9d87f6e6f3e094de1961db1fdcb37f7611de9efda60992c01",
            `da_output/de_gcms_metabolites_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_metabolites_groupKO-groupWT_long_annot.tsv` = "3fe65bf0975381006239ba8531060b86c64ba5d9279f804722b1938d6cc28381",
            `da_output/de_posmode_metabolites_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `feature_qc/ruv_normalised_results.RDS` = "<present>",
            `feature_qc/ruv_normalised_results.tsv` = "f87484c97f41f65f89307e33c596dfaf95499e9ce230ff4f79764c6bb795a2dd",
            `pathway_enrichment/KO_vs_WT_down_enrichment_results.tsv` = "278827ebf4745b50520d3714c4389fd9c0a3c817e116638b4227257b9ed2a07a",
            `pathway_enrichment/KO_vs_WT_up_enrichment_results.tsv` = "278827ebf4745b50520d3714c4389fd9c0a3c817e116638b4227257b9ed2a07a",
            `scripts/metabolomics/metabolomics_report.rmd` = "<present>",
            `source/contrasts_tbl.tab` = "7209220541b9240edc6b7b73030b7384687b68cd7814091847ac8baf36891f72",
            `source/design_matrix.tab` = "c3435c79da809286aebfaa17f024efdea7282282a0904b2f5f1477f4790b6983",
            `source/study_parameters.txt` = "<present>"),
        study_parameters = "1 lines 55c1b3aadde3f3b7",
        integration = stats::setNames(character(0), character()),
        session_state = "<0 files>",
        errors = "<unhandled> unused arguments (self, name)"),
    save_record_fails = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Study parameters created for: kit108 \nDescription: Full metabolomics analysis workflow with normalization and DA \nTimestamp: <TIME> \nFile: <ROOT>/source/study_parameters.txt \nSource: Final S4 object @args + config_list \nIntegration Object: metabolomics_kit108_final_s4.RDS \nStatus: Parameters saved [OK]",
            template_status = "Template: Metabolomics Report [OK]"),
        notes = c(
            "[warning NULL -] The analysis record (analysis_record.json) was not written: kit108 record double fails",
            "[message 5 -] Saved Integration S4 Object", "[message 5 -] Study parameters saved successfully"),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = c(
            "createWorkflowArgsFromConfig(workflow_name=kit108, description=Full metabolomics analysis workflow with normalization and DA, organism_name=<NULL>, taxon_id=<NULL>, source_dir_path=<ROOT>/source, final_s4_object=<MetaboliteAssayData 075de31ac869b80e>, contrasts_tbl=<table 2x3 4323c37e5b4e992a>, workflow_data=<environment>)",
            "writeWorkflowAnalysisRecord(manager=<WorkflowStateOwner>, path=<ROOT>/source/analysis_record.json)"),
        global = c(before = "the session's", after = "the session's"),
        files = c(`integration/metabolomics_kit108_final_s4.RDS` = "<present>"),
        study_parameters = "127 lines 6f2dd5b056f531b5",
        integration = c(metabolomics_kit108_final_s4.RDS = "<RDS> MetaboliteAssayData 075de31ac869b80e"),
        session_state = "<unchanged>",
        errors = character(0)),
    copy = list(
        outputs = c("copy_status", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory with optional missing artifacts [WARNING]",
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full metabolomics analysis workflow with normalization and DA \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK]",
            template_status = "Template: Metabolomics Report [OK]"),
        notes = "[warning 5 -] Publication files copied with optional missing artifacts",
        updates = character(0),
        progress = "start Copying files to publication directory... 0.1",
        reactive_values = c(
            files_copied = "TRUE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = c(
            "copyToResultsSummary(omic_type=metabolomics, experiment_label=kit108, contrasts_tbl=<table 2x3 4323c37e5b4e992a>, design_matrix=<table 6x5 be07121076802bfa>, force=TRUE, current_rmd=<NULL>, project_dirs_object_name=project_dirs, global=the session's)",
            "  copyToResultsSummary returned: Composite QC (PNG) | DA Heatmaps"),
        global = c(before = "the session's", after = "the session's"),
        files = c(
            `results_summary/Publication_figures/Enrichment_Plots/KO_vs_WT_down_enrichment_results.tsv` = "278827ebf4745b50520d3714c4389fd9c0a3c817e116638b4227257b9ed2a07a",
            `results_summary/Publication_figures/Enrichment_Plots/KO_vs_WT_up_enrichment_results.tsv` = "278827ebf4745b50520d3714c4389fd9c0a3c817e116638b4227257b9ed2a07a",
            `results_summary/Publication_tables/DA_results_metabolomics.xlsx` = "<present>",
            `results_summary/Publication_tables/Pathway_enrichment_results_metabolomics.xlsx` = "<present>",
            `results_summary/Publication_tables/RUV_normalised_results.tsv` = "f87484c97f41f65f89307e33c596dfaf95499e9ce230ff4f79764c6bb795a2dd",
            `results_summary/Publication_tables/ruv_normalised_results.RDS` = "<present>",
            `results_summary/Study_report/contrasts_tbl.tab` = "7209220541b9240edc6b7b73030b7384687b68cd7814091847ac8baf36891f72",
            `results_summary/Study_report/design_matrix.tab` = "c3435c79da809286aebfaa17f024efdea7282282a0904b2f5f1477f4790b6983",
            `results_summary/Study_report/study_parameters.txt` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    report_bundle_fails = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory with optional missing artifacts [WARNING]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full metabolomics analysis workflow with normalization and DA \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK] \nReport location: <ROOT>/results_summary/metabolomics_report_metabolomics_kit108.docx",
            template_status = "Template: Metabolomics Report [OK]"),
        notes = c(
            "[warning NULL -] The analysis record bundle was not exported: kit108 bundle double fails",
            "[message 5 -] Report generated successfully!"),
        updates = character(0),
        progress = c(
            "start Generating report... 0.1", "0.5 Rendering report...",
            "0.9 Exporting the analysis record..."),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/metabolomics_report_metabolomics_kit108.docx",
            workflow_args_saved = "TRUE"),
        calls = c(
            "RenderReport(omic_type=metabolomics, experiment_label=kit108, rmd_filename=metabolomics_report.rmd, project_dirs_object_name=project_dirs, output_format=<NULL>, global=the session's)",
            "exportWorkflowAnalysisRecordBundle(manager=<WorkflowStateOwner>, directory=<ROOT>/results_summary)"),
        global = c(before = "the session's", after = "the session's"),
        files = c(`results_summary/metabolomics_report_metabolomics_kit108.docx` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = "<unhandled> unused arguments (self, name)"))

omicSummaryW1Frozen[["metab_required_missing"]] <- list(
    start = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = character(0),
        updates = "summary-experiment_label: <same> -> kit108",
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = character(0),
        study_parameters = "<absent>",
        integration = stats::setNames(character(0), character()),
        session_state = "<0 files>",
        errors = "<unhandled> unused arguments (self, name)"),
    copy = list(
        outputs = c("copy_status", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Error: Required publication artifacts missing: Contrasts Table, Design Matrix",
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = "[error 10 -] Copy error: Required publication artifacts missing: Contrasts Table, Design Matrix",
        updates = character(0),
        progress = "start Copying files to publication directory... 0.1",
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = c(
            "copyToResultsSummary(omic_type=metabolomics, experiment_label=kit108, contrasts_tbl=<NULL>, design_matrix=<NULL>, force=TRUE, current_rmd=<NULL>, project_dirs_object_name=project_dirs, global=the session's)",
            "  copyToResultsSummary returned: Contrasts Table | Design Matrix | Composite QC (PNG) | DA Heatmaps"),
        global = c(before = "absent", after = "the session's"),
        files = c(
            `results_summary/Publication_tables/DA_results_metabolomics.xlsx` = "<present>",
            `results_summary/Publication_tables/Pathway_enrichment_results_metabolomics.xlsx` = "<present>",
            `results_summary/Study_report/study_parameters.txt` = "<present>",
            `source/study_parameters.txt` = "<present>"),
        study_parameters = "7 lines b454092b567d4211",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    report = list(
        outputs = c("copy_status", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Error: Required publication artifacts missing: Contrasts Table, Design Matrix",
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = "[warning 5 -] Please copy the results to the publication directory before generating the report (Copy to Publication Directory).",
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = character(0),
        global = c(before = "the session's", after = "the session's"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)))

omicSummaryW1Frozen[["metab_github_fallback"]] <- list(
    copy = list(
        outputs = c("copy_status", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full metabolomics analysis workflow with normalization and DA \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK]",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = "[message 5 -] Publication files copied",
        updates = "summary-experiment_label: <same> -> kit108",
        progress = "start Copying files to publication directory... 0.1",
        reactive_values = c(
            files_copied = "TRUE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = c(
            "copyToResultsSummary(omic_type=metabolomics, experiment_label=kit108, contrasts_tbl=<table 2x3 4323c37e5b4e992a>, design_matrix=<table 6x5 be07121076802bfa>, force=TRUE, current_rmd=<NULL>, project_dirs_object_name=project_dirs, global=the session's)",
            "  copyToResultsSummary returned: "),
        global = c(before = "the session's", after = "the session's"),
        files = c(
            `da_output/de_gcms_metabolites_groupKO-groupWT_long_annot.tsv` = "13be76abf9da00e9d87f6e6f3e094de1961db1fdcb37f7611de9efda60992c01",
            `da_output/de_gcms_metabolites_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_metabolites_groupKO-groupWT_long_annot.tsv` = "3fe65bf0975381006239ba8531060b86c64ba5d9279f804722b1938d6cc28381",
            `da_output/de_posmode_metabolites_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `feature_qc/gcms_post_norm_density.png` = "<present>",
            `feature_qc/itsd_lcms_pos.png` = "<present>",
            `feature_qc/lcms_pos_pre_norm_pca.png` = "<present>",
            `feature_qc/lcms_pos_ruv_corrected_rle.png` = "<present>",
            `feature_qc/metabolomics_composite_QC_figure.pdf` = "<present>",
            `feature_qc/metabolomics_composite_QC_figure.png` = "<present>",
            `feature_qc/ruv_normalised_results.RDS` = "<present>",
            `feature_qc/ruv_normalised_results.tsv` = "f87484c97f41f65f89307e33c596dfaf95499e9ce230ff4f79764c6bb795a2dd",
            `pathway_enrichment/KO_vs_WT_down_enrichment_results.tsv` = "278827ebf4745b50520d3714c4389fd9c0a3c817e116638b4227257b9ed2a07a",
            `pathway_enrichment/KO_vs_WT_up_enrichment_results.tsv` = "278827ebf4745b50520d3714c4389fd9c0a3c817e116638b4227257b9ed2a07a",
            `publication_graphs/Heatmap/artifact.txt` = "<present>",
            `publication_graphs/Heatmaps/artifact.txt` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/metabolites_num_sig_da_molecules.tab` = "32c301d1f9145915cfa6a4789f1c5c7a9b8f69ac15917d31f28bb317eab608c9",
            `publication_graphs/Volcano_Plots/artifact.txt` = "<present>",
            `results_summary/Publication_figures/Enrichment_Plots/KO_vs_WT_down_enrichment_results.tsv` = "278827ebf4745b50520d3714c4389fd9c0a3c817e116638b4227257b9ed2a07a",
            `results_summary/Publication_figures/Enrichment_Plots/KO_vs_WT_up_enrichment_results.tsv` = "278827ebf4745b50520d3714c4389fd9c0a3c817e116638b4227257b9ed2a07a",
            `results_summary/Publication_figures/Heatmaps/artifact.txt` = "<present>",
            `results_summary/Publication_tables/DA_results_metabolomics.xlsx` = "<present>",
            `results_summary/Publication_tables/Pathway_enrichment_results_metabolomics.xlsx` = "<present>",
            `results_summary/Publication_tables/RUV_normalised_results.tsv` = "f87484c97f41f65f89307e33c596dfaf95499e9ce230ff4f79764c6bb795a2dd",
            `results_summary/Publication_tables/da_metabolomics_num_sig_da_molecules.tab` = "32c301d1f9145915cfa6a4789f1c5c7a9b8f69ac15917d31f28bb317eab608c9",
            `results_summary/Publication_tables/ruv_normalised_results.RDS` = "<present>",
            `results_summary/QC_figures/12_correlation_filtered_combined_plots.png` = "<present>",
            `results_summary/QC_figures/gcms_post_norm_density.png` = "<present>",
            `results_summary/QC_figures/itsd_lcms_pos.png` = "<present>",
            `results_summary/QC_figures/lcms_pos_pre_norm_pca.png` = "<present>",
            `results_summary/QC_figures/lcms_pos_ruv_corrected_rle.png` = "<present>",
            `results_summary/QC_figures/metabolomics_composite_QC_figure.png` = "<present>",
            `results_summary/Study_report/contrasts_tbl.tab` = "7209220541b9240edc6b7b73030b7384687b68cd7814091847ac8baf36891f72",
            `results_summary/Study_report/design_matrix.tab` = "c3435c79da809286aebfaa17f024efdea7282282a0904b2f5f1477f4790b6983",
            `results_summary/Study_report/study_parameters.txt` = "<present>",
            `source/contrasts_tbl.tab` = "7209220541b9240edc6b7b73030b7384687b68cd7814091847ac8baf36891f72",
            `source/design_matrix.tab` = "c3435c79da809286aebfaa17f024efdea7282282a0904b2f5f1477f4790b6983",
            `source/study_parameters.txt` = "<present>",
            `time/12_correlation_filtered_combined_plots.png` = "<present>"),
        study_parameters = "1 lines 55c1b3aadde3f3b7",
        integration = stats::setNames(character(0), character()),
        session_state = "<0 files>",
        errors = "<unhandled> unused arguments (self, name)"),
    report_fallback = list(
        outputs = c("copy_status", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full metabolomics analysis workflow with normalization and DA \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK]",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = c(
            "[error 10 -] Template retrieval failed: kit108 download double: offline",
            "[error 15 -] Report generation failed: kit108 render double: no template at <ROOT>/scripts/metabolomics/metabolomics_report.rmd",
            "[warning 10 -] Debug info: Check R console for detailed error trace"),
        updates = character(0),
        progress = c(
            "start Generating report... 0.1", "0.2 Checking for metabolomics_report.rmd template...",
            "0.5 Rendering report..."),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = c(
            "system.file(path=reports/metabolomics/metabolomics_report.rmd, package=MultiScholaR)",
            "download.file(url=https://raw.githubusercontent.com/APAF-bioinformatics/MultiScholaR/main/inst/reports/metabolomics/metabolomics_report.rmd, destfile=<ROOT>/scripts/metabolomics/metabolomics_report.rmd, quiet=TRUE)",
            "RenderReport(omic_type=metabolomics, experiment_label=kit108, rmd_filename=metabolomics_report.rmd, project_dirs_object_name=project_dirs, output_format=<NULL>, global=the session's)"),
        global = c(before = "the session's", after = "the session's"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)))

omicSummaryW1Frozen[["metab_managed"]] <- list(
    start = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = character(0),
        updates = "summary-experiment_label: <same> -> kit108",
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = c(
            `da_output/de_gcms_metabolites_groupKO-groupWT_long_annot.tsv` = "13be76abf9da00e9d87f6e6f3e094de1961db1fdcb37f7611de9efda60992c01",
            `da_output/de_gcms_metabolites_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_metabolites_groupKO-groupWT_long_annot.tsv` = "3fe65bf0975381006239ba8531060b86c64ba5d9279f804722b1938d6cc28381",
            `da_output/de_posmode_metabolites_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `feature_qc/gcms_post_norm_density.png` = "<present>",
            `feature_qc/itsd_lcms_pos.png` = "<present>",
            `feature_qc/lcms_pos_pre_norm_pca.png` = "<present>",
            `feature_qc/lcms_pos_ruv_corrected_rle.png` = "<present>",
            `feature_qc/metabolomics_composite_QC_figure.pdf` = "<present>",
            `feature_qc/metabolomics_composite_QC_figure.png` = "<present>",
            `feature_qc/ruv_normalised_results.RDS` = "<present>",
            `feature_qc/ruv_normalised_results.tsv` = "f87484c97f41f65f89307e33c596dfaf95499e9ce230ff4f79764c6bb795a2dd",
            `pathway_enrichment/KO_vs_WT_down_enrichment_results.tsv` = "278827ebf4745b50520d3714c4389fd9c0a3c817e116638b4227257b9ed2a07a",
            `pathway_enrichment/KO_vs_WT_up_enrichment_results.tsv` = "278827ebf4745b50520d3714c4389fd9c0a3c817e116638b4227257b9ed2a07a",
            `publication_graphs/Heatmap/artifact.txt` = "<present>",
            `publication_graphs/Heatmaps/artifact.txt` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/metabolites_num_sig_da_molecules.tab` = "32c301d1f9145915cfa6a4789f1c5c7a9b8f69ac15917d31f28bb317eab608c9",
            `publication_graphs/Volcano_Plots/artifact.txt` = "<present>",
            `source/contrasts_tbl.tab` = "7209220541b9240edc6b7b73030b7384687b68cd7814091847ac8baf36891f72",
            `source/design_matrix.tab` = "c3435c79da809286aebfaa17f024efdea7282282a0904b2f5f1477f4790b6983",
            `source/study_parameters.txt` = "<present>",
            `time/12_correlation_filtered_combined_plots.png` = "<present>"),
        study_parameters = "1 lines 55c1b3aadde3f3b7",
        integration = stats::setNames(character(0), character()),
        session_state = "<0 files>",
        errors = "<unhandled> unused arguments (self, name)"),
    save = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Study parameters created for: kit108 \nDescription: Full metabolomics analysis workflow with normalization and DA \nTimestamp: <TIME> \nFile: <ROOT>/source/study_parameters.txt \nSource: Final S4 object @args + config_list \nIntegration Object: metabolomics_kit108_final_s4.RDS \nStatus: Parameters saved [OK]",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = c(
            "[default NULL metabolomics_action_working] Applying metabolomics step ...",
            "remove metabolomics_action_working", "[message 5 -] Saved Integration S4 Object",
            "[message 5 -] Study parameters saved successfully"),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = c(
            "createWorkflowArgsFromConfig(workflow_name=kit108, description=Full metabolomics analysis workflow with normalization and DA, organism_name=<NULL>, taxon_id=<NULL>, source_dir_path=<ROOT>/source, final_s4_object=<MetaboliteAssayData 9ed530ceceb0a1e7>, contrasts_tbl=<table 2x3 4323c37e5b4e992a>, workflow_data=<environment>)",
            "writeWorkflowAnalysisRecord(manager=<WorkflowStateOwner>, path=<ROOT>/source/analysis_record.json)",
            "omicsManagedProduct(operation=metab_integration_export, parameters=path:<ROOT>/integration/metabolomics_kit108_final_s4.RDS;state_name:metab_correlation_filtered)"),
        global = c(before = "absent", after = "absent"),
        files = c(
            `integration/metabolomics_kit108_final_s4.RDS` = "<present>",
            `source/analysis_record.json` = "<present>"),
        study_parameters = "127 lines 6f2dd5b056f531b5",
        integration = c(metabolomics_kit108_final_s4.RDS = "<RDS> MetaboliteAssayData 075de31ac869b80e"),
        session_state = "<unchanged>",
        errors = character(0)),
    export = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Study parameters created for: kit108 \nDescription: Full metabolomics analysis workflow with normalization and DA \nTimestamp: <TIME> \nFile: <ROOT>/source/study_parameters.txt \nSource: Final S4 object @args + config_list \nIntegration Object: metabolomics_kit108_final_s4.RDS \nStatus: Parameters saved [OK]",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = "[message 5 -] Session state exported to: <ROOT>/source/session_state_<DATE>.RDS",
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = c(`source/session_state_<DATE>.RDS` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = c(
            experiment_label = "kit108",
            description = "Full metabolomics analysis workflow with normalization and DA",
            timestamp = "<time>", omic_type = "metabolomics", workflow_args_saved = "TRUE",
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            project_dirs = "<digest> 4701bf62c551a512", workflow_type = "metabolomics",
            report_template = "metabolomics_report.rmd", parameter_payload = "<digest> d66aa0b35bc9787b",
            assay_names = "LCMS_Pos|GCMS    ", feature_counts = "6|6", sample_count = "6",
            contrast_count = "2"),
        errors = character(0)))

omicSummaryW1Frozen[["trx_main"]] <- list(
    norm = list(
        outputs = c(
            "report_ready", "ruviii_before_after", "ruviii_m_ui", "ruviii_status", "session_summary",
            "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            ruviii_before_after = "<table> 62875e4feeff8105b1284f0e99e836668bf7b7b93e25535e3c553a92a7f06583",
            ruviii_m_ui = "<ui> e7ae5dc12fe193d9706149295d232a88b3b7a9e15fced503f10f12bc405b104f",
            ruviii_status = "Not run yet.",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = "[message 5 -] Normalization complete! Proceeding to DE analysis.",
        updates = "summary-experiment_label: <same> -> kit108",
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = c(
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "23",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "5",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/partitioned-table-<ID>` = "2",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_expression_filter-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_ruv_grid-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_ruvs-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_tmm_logcpm-<ID>` = "1",
            `results/transcriptomics/gene_qc/ruv_normalised_log_values.tsv` = "83de1dbdd1018161c5ec6c1dac4dc179b0c1d0340fd2fe826d4b3c9b4ae93387",
            `results/transcriptomics/gene_qc/trx_raw_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_raw_rle.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_ruvs_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_ruvs_rle.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_tmm_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_tmm_rle.png` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_optimization_summary.xlsx` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_panels.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_panels.png` = "<present>",
            `results/transcriptomics/publication_graphs/normalization_comparison_grid.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/normalization_comparison_grid.png` = "<present>",
            `scripts/transcriptomics/design_matrix.tab` = "e3b6f40550967e975506c582126b035f5fcc887aae4ffe97e3200f3999a3584e",
            `scripts/transcriptomics/import_sources/counts.tsv` = "2c9e3b6290daf902297e7dab06425b831bfeb6261bd5b5d1819e6fd3e21dde49"),
        study_parameters = "<absent>",
        integration = stats::setNames(character(0), character()),
        session_state = "<0 files>",
        errors = "<unhandled> unused arguments (self, name)"),
    da = list(
        outputs = c(
            "report_ready", "ruviii_before_after", "ruviii_m_ui", "ruviii_status", "session_summary",
            "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            ruviii_before_after = "<table> 62875e4feeff8105b1284f0e99e836668bf7b7b93e25535e3c553a92a7f06583",
            ruviii_m_ui = "<ui> e7ae5dc12fe193d9706149295d232a88b3b7a9e15fced503f10f12bc405b104f",
            ruviii_status = "Not run yet.",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = c(
            "[default NULL da_running] Running differential expression analysis...", "remove da_running",
            "[message 5 -] Differential expression analysis complete!",
            "[message 5 -] DA results saved to disk (tables, volcano plots, heatmaps)"),
        updates = c(
            "de-volcano_contrast: Sera_vs_RPMI -> Sera_vs_RPMI",
            "de-heatmap_contrast: Sera_vs_RPMI -> Sera_vs_RPMI",
            "de-table_contrast: Sera_vs_RPMI -> Sera_vs_RPMI", "de-volcano_assay: Combined|counts -> counts",
            "de-heatmap_assay: Combined|counts -> Combined", "de-table_assay: All -> All"),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = c(
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "28",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "6",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx-da-fit-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_da_contrast-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_da_contrast-da_contrast-<ID>` = "1",
            `results/transcriptomics/da_genes/de_results_Sera_vs_RPMI.csv` = "c11158b3af965955f4284fad607de9462097852a94c2b074ffab0e78c0a4af6b",
            `results/transcriptomics/da_genes/ruv_optimization_summary.txt` = "<present>",
            `results/transcriptomics/da_genes/top_de_genes.tsv` = "10aa0f7b3fa91cdc9137a9e4f11dca3a1413f4bd4e236542637c471bf87daa21",
            `results/transcriptomics/publication_graphs/NumSigDeMolecules/summary_barplot_Sera_vs_RPMI.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/NumSigDeMolecules/summary_barplot_Sera_vs_RPMI.png` = "<present>",
            `results/transcriptomics/publication_graphs/Volcano_Plots/volcano_plot_Sera_vs_RPMI.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/Volcano_Plots/volcano_plot_Sera_vs_RPMI.png` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = c(
            "<unhandled> argument is of length zero", "<unhandled> argument is of length zero")),
    report_before_copy = list(
        outputs = c(
            "report_ready", "ruviii_before_after", "ruviii_m_ui", "ruviii_status", "session_summary",
            "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            ruviii_before_after = "<table> 62875e4feeff8105b1284f0e99e836668bf7b7b93e25535e3c553a92a7f06583",
            ruviii_m_ui = "<ui> e7ae5dc12fe193d9706149295d232a88b3b7a9e15fced503f10f12bc405b104f",
            ruviii_status = "Not run yet.",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = "[warning 5 -] Please copy the results to the publication directory before generating the report (Copy to Publication Directory).",
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    save_before_export = list(
        outputs = c(
            "report_ready", "ruviii_before_after", "ruviii_m_ui", "ruviii_status", "session_summary",
            "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            ruviii_before_after = "<table> 62875e4feeff8105b1284f0e99e836668bf7b7b93e25535e3c553a92a7f06583",
            ruviii_m_ui = "<ui> e7ae5dc12fe193d9706149295d232a88b3b7a9e15fced503f10f12bc405b104f",
            ruviii_status = "Not run yet.",
            session_summary = "Study parameters created for: kit108 \nDescription: Full transcriptomics analysis workflow with normalization, DE and the RUVIII integration export \nTimestamp: <TIME> \nFile: <ROOT>/scripts/transcriptomics/study_parameters.txt \nSource: Final S4 object @args + config_list \nIntegration Object: not written (no integration export yet) \nStatus: Parameters saved [OK]",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = c(
            "[warning 10 -] The integration export needs the RUVIII export: run it in the panel above, then save again.",
            "[message 5 -] Saved Integration S4 Object", "[message 5 -] Study parameters saved successfully"),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = c(
            "createWorkflowArgsFromConfig(workflow_name=kit108, description=Full transcriptomics analysis workflow with normalization, DE and the RUVIII integration export, organism_name=<NULL>, taxon_id=<NULL>, source_dir_path=<ROOT>/scripts/transcriptomics, final_s4_object=<ref>, contrasts_tbl=<table 1x3 ba384cec7aecd593>, workflow_data=<environment>)",
            "writeWorkflowAnalysisRecord(manager=<WorkflowStateOwner>, path=<ROOT>/scripts/transcriptomics/analysis_record.json)"),
        global = c(before = "absent", after = "absent"),
        files = c(
            `scripts/transcriptomics/analysis_record.json` = "<present>",
            `scripts/transcriptomics/study_parameters.txt` = "<present>"),
        study_parameters = "54 lines 8209e0aae98f1f65",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    ruviii_unavailable = list(
        outputs = c(
            "report_ready", "ruviii_before_after", "ruviii_m_ui", "ruviii_status", "session_summary",
            "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            ruviii_before_after = "<table> 62875e4feeff8105b1284f0e99e836668bf7b7b93e25535e3c553a92a7f06583",
            ruviii_m_ui = "<ui> e7ae5dc12fe193d9706149295d232a88b3b7a9e15fced503f10f12bc405b104f",
            ruviii_status = "Not run yet.",
            session_summary = "Study parameters created for: kit108 \nDescription: Full transcriptomics analysis workflow with normalization, DE and the RUVIII integration export \nTimestamp: <TIME> \nFile: <ROOT>/scripts/transcriptomics/study_parameters.txt \nSource: Final S4 object @args + config_list \nIntegration Object: not written (no integration export yet) \nStatus: Parameters saved [OK]",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = "[error 15 -] RUVIII export error: the technical-replicate column is absent (trx_ruviii_replicates_invalid)",
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0),
        ruviii = "<none>"),
    ruviii = list(
        outputs = c(
            "report_ready", "ruviii_before_after", "ruviii_m_ui", "ruviii_status", "session_summary",
            "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            ruviii_before_after = "<plot> alt=RUVIII before and after",
            ruviii_m_ui = "<ui> e7ae5dc12fe193d9706149295d232a88b3b7a9e15fced503f10f12bc405b104f",
            ruviii_status = "RUVIII export complete (integration only; DA results unchanged).\nM = group\nk = 2, 18 inferred negative controls",
            session_summary = "Study parameters created for: kit108 \nDescription: Full transcriptomics analysis workflow with normalization, DE and the RUVIII integration export \nTimestamp: <TIME> \nFile: <ROOT>/scripts/transcriptomics/study_parameters.txt \nSource: Final S4 object @args + config_list \nIntegration Object: not written (no integration export yet) \nStatus: Parameters saved [OK]",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = c(
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "31",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_ruviii_export-<ID>` = "1",
            `results/transcriptomics/gene_qc/trx_ruviii_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_ruviii_rle.png` = "<present>",
            `results/transcriptomics/publication_graphs/ruviii_before_and_after.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/ruviii_before_and_after.png` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0),
        ruviii = c(m_source = "group", k = "2", controls = "18", newY = "12x187 6fc806b39ca98cbc")),
    save = list(
        outputs = c(
            "report_ready", "ruviii_before_after", "ruviii_m_ui", "ruviii_status", "session_summary",
            "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            ruviii_before_after = "<plot> alt=RUVIII before and after",
            ruviii_m_ui = "<ui> e7ae5dc12fe193d9706149295d232a88b3b7a9e15fced503f10f12bc405b104f",
            ruviii_status = "RUVIII export complete (integration only; DA results unchanged).\nM = group\nk = 2, 18 inferred negative controls",
            session_summary = "Study parameters created for: kit108 \nDescription: Full transcriptomics analysis workflow with normalization, DE and the RUVIII integration export \nTimestamp: <TIME> \nFile: <ROOT>/scripts/transcriptomics/study_parameters.txt \nSource: Final S4 object @args + config_list \nIntegration Object: transcriptomics_normalized_data.tsv \nStatus: Parameters saved [OK]",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = c(
            "[message 5 -] Saved Integration S4 Object", "[message 5 -] Study parameters saved successfully"),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = c(
            "createWorkflowArgsFromConfig(workflow_name=kit108, description=Full transcriptomics analysis workflow with normalization, DE and the RUVIII integration export, organism_name=<NULL>, taxon_id=<NULL>, source_dir_path=<ROOT>/scripts/transcriptomics, final_s4_object=<ref>, contrasts_tbl=<table 1x3 ba384cec7aecd593>, workflow_data=<environment>)",
            "writeWorkflowAnalysisRecord(manager=<WorkflowStateOwner>, path=<ROOT>/scripts/transcriptomics/analysis_record.json)"),
        global = c(before = "absent", after = "absent"),
        files = c(
            `integration/transcriptomics_normalized_data.tsv` = "7b9d35dbf33331d16b6609caacbd9445879924bd8bb548b6821dfd649bfdd08e"),
        study_parameters = "57 lines 6f90092f3d590fd4",
        integration = c(transcriptomics_normalized_data.tsv = "7b9d35dbf33331d1"),
        session_state = "<unchanged>",
        errors = character(0)),
    copy = list(
        outputs = c(
            "copy_status", "report_ready", "ruviii_before_after", "ruviii_m_ui", "ruviii_status",
            "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            report_ready = "<error> unused arguments (self, name)",
            ruviii_before_after = "<plot> alt=RUVIII before and after",
            ruviii_m_ui = "<ui> e7ae5dc12fe193d9706149295d232a88b3b7a9e15fced503f10f12bc405b104f",
            ruviii_status = "RUVIII export complete (integration only; DA results unchanged).\nM = group\nk = 2, 18 inferred negative controls",
            session_summary = "Workflow args created for: kit108 \nDescription: Full transcriptomics analysis workflow with normalization, DE and the RUVIII integration export \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK]",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = "[message 5 -] Publication files copied",
        updates = character(0),
        progress = "start Copying files to publication directory... 0.1",
        reactive_values = c(
            files_copied = "TRUE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = c(
            "copyToResultsSummary(omic_type=transcriptomics, experiment_label=kit108, contrasts_tbl=<table 1x3 ba384cec7aecd593>, design_matrix=<table 12x3 440a505e731ef388>, force=TRUE, current_rmd=<NULL>, project_dirs_object_name=project_dirs, global=the session's)",
            "  copyToResultsSummary returned: "),
        global = c(before = "absent", after = "the session's"),
        files = c(
            `results_summary/transcriptomics/Publication_figures/NumSigDeMolecules/summary_barplot_Sera_vs_RPMI.pdf` = "<present>",
            `results_summary/transcriptomics/Publication_figures/NumSigDeMolecules/summary_barplot_Sera_vs_RPMI.png` = "<present>",
            `results_summary/transcriptomics/Publication_figures/Volcano_Plots/volcano_plot_Sera_vs_RPMI.pdf` = "<present>",
            `results_summary/transcriptomics/Publication_figures/Volcano_Plots/volcano_plot_Sera_vs_RPMI.png` = "<present>",
            `results_summary/transcriptomics/Publication_figures/canonical_correlation_panels.pdf` = "<present>",
            `results_summary/transcriptomics/Publication_figures/canonical_correlation_panels.png` = "<present>",
            `results_summary/transcriptomics/Publication_tables/DA_results_transcriptomics.xlsx` = "<present>",
            `results_summary/transcriptomics/Publication_tables/Pathway_enrichment_results_transcriptomics.xlsx` = "<present>",
            `results_summary/transcriptomics/Publication_tables/canonical_correlation_optimization_summary.xlsx` = "<present>",
            `results_summary/transcriptomics/Publication_tables/de_results_Sera_vs_RPMI.csv` = "c11158b3af965955f4284fad607de9462097852a94c2b074ffab0e78c0a4af6b",
            `results_summary/transcriptomics/Publication_tables/ruv_normalised_log_values.tsv` = "83de1dbdd1018161c5ec6c1dac4dc179b0c1d0340fd2fe826d4b3c9b4ae93387",
            `results_summary/transcriptomics/Publication_tables/top_de_genes.tsv` = "10aa0f7b3fa91cdc9137a9e4f11dca3a1413f4bd4e236542637c471bf87daa21",
            `results_summary/transcriptomics/QC_figures/normalization_comparison_grid.pdf` = "<present>",
            `results_summary/transcriptomics/QC_figures/normalization_comparison_grid.png` = "<present>",
            `results_summary/transcriptomics/QC_figures/ruviii_before_and_after.pdf` = "<present>",
            `results_summary/transcriptomics/QC_figures/ruviii_before_and_after.png` = "<present>",
            `results_summary/transcriptomics/Study_report/contrasts_tbl.tab` = "4c9856ae0044c3e20d9ae4722742e7eab0f2dabb11b18668162461912fd69d87",
            `results_summary/transcriptomics/Study_report/design_matrix.tab` = "e3b6f40550967e975506c582126b035f5fcc887aae4ffe97e3200f3999a3584e",
            `results_summary/transcriptomics/Study_report/study_parameters.txt` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    report = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "ruviii_before_after", "ruviii_m_ui",
            "ruviii_status", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            ruviii_before_after = "<plot> alt=RUVIII before and after",
            ruviii_m_ui = "<ui> e7ae5dc12fe193d9706149295d232a88b3b7a9e15fced503f10f12bc405b104f",
            ruviii_status = "RUVIII export complete (integration only; DA results unchanged).\nM = group\nk = 2, 18 inferred negative controls",
            session_summary = "Workflow args created for: kit108 \nDescription: Full transcriptomics analysis workflow with normalization, DE and the RUVIII integration export \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK] \nReport location: <ROOT>/results_summary/transcriptomics/transcriptomics_report_transcriptomics_kit108.docx",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = c(
            "[message 5 -] transcriptomics_report.rmd template copied from package",
            "[message 5 -] Report generated successfully!"),
        updates = character(0),
        progress = c(
            "start Generating report... 0.1", "0.2 Checking for transcriptomics_report.rmd template...",
            "0.5 Rendering report...", "0.9 Exporting the analysis record..."),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/transcriptomics/transcriptomics_report_transcriptomics_kit108.docx",
            workflow_args_saved = "TRUE"),
        calls = c(
            "RenderReport(omic_type=transcriptomics, experiment_label=kit108, rmd_filename=transcriptomics_report.rmd, project_dirs_object_name=project_dirs, output_format=<NULL>, global=the session's)",
            "exportWorkflowAnalysisRecordBundle(manager=<WorkflowStateOwner>, directory=<ROOT>/results_summary/transcriptomics/Study_report)",
            "writeWorkflowAnalysisRecord(manager=<WorkflowStateOwner>, path=<ROOT>/results_summary/transcriptomics/Study_report/analysis_record.json)"),
        global = c(before = "the session's", after = "the session's"),
        files = c(
            `results_summary/transcriptomics/Study_report/analysis_record.docx` = "<present>",
            `results_summary/transcriptomics/Study_report/analysis_record.json` = "<present>",
            `results_summary/transcriptomics/Study_report/analysis_record.xlsx` = "<present>",
            `results_summary/transcriptomics/transcriptomics_report_transcriptomics_kit108.docx` = "<present>",
            `scripts/transcriptomics/transcriptomics_report.rmd` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = "<unhandled> unused arguments (self, name)"),
    export = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "ruviii_before_after", "ruviii_m_ui",
            "ruviii_status", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            ruviii_before_after = "<plot> alt=RUVIII before and after",
            ruviii_m_ui = "<ui> e7ae5dc12fe193d9706149295d232a88b3b7a9e15fced503f10f12bc405b104f",
            ruviii_status = "RUVIII export complete (integration only; DA results unchanged).\nM = group\nk = 2, 18 inferred negative controls",
            session_summary = "Workflow args created for: kit108 \nDescription: Full transcriptomics analysis workflow with normalization, DE and the RUVIII integration export \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK] \nReport location: <ROOT>/results_summary/transcriptomics/transcriptomics_report_transcriptomics_kit108.docx",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = "[message 5 -] Session state exported to: <ROOT>/scripts/transcriptomics/session_state_<DATE>.RDS",
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/transcriptomics/transcriptomics_report_transcriptomics_kit108.docx",
            workflow_args_saved = "TRUE"),
        calls = character(0),
        global = c(before = "the session's", after = "the session's"),
        files = c(`scripts/transcriptomics/session_state_<DATE>.RDS` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = c(
            experiment_label = "kit108",
            description = "Full transcriptomics analysis workflow with normalization, DE and the RUVIII integration export",
            timestamp = "<time>", omic_type = "transcriptomics", workflow_args_saved = "TRUE",
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/transcriptomics/transcriptomics_report_transcriptomics_kit108.docx",
            project_dirs = "<digest> 9e585ab288403615", workflow_type = "transcriptomics",
            report_template = "transcriptomics_report.rmd", parameter_payload = "<digest> 4a786b50c11d43f0",
            assay_names = "counts", feature_counts = "187", sample_count = "12", contrast_count = "1"),
        errors = character(0)))

omicSummaryW1Frozen[["trx_technical"]] <- list(
    norm = list(
        outputs = c(
            "report_ready", "ruviii_before_after", "ruviii_m_ui", "ruviii_status", "session_summary",
            "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            ruviii_before_after = "<table> 62875e4feeff8105b1284f0e99e836668bf7b7b93e25535e3c553a92a7f06583",
            ruviii_m_ui = "<ui> aafb6bcfd2a27f3e561956222696cdc9fd29de5f9954d8d5cb1e2c4b4bbd0ce9",
            ruviii_status = "Not run yet.",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = "[message 5 -] Normalization complete! Proceeding to DE analysis.",
        updates = "summary-experiment_label: <same> -> kit108",
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = c(
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "23",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "5",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/partitioned-table-<ID>` = "2",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_expression_filter-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_ruv_grid-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_ruvs-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_tmm_logcpm-<ID>` = "1",
            `results/transcriptomics/gene_qc/ruv_normalised_log_values.tsv` = "83de1dbdd1018161c5ec6c1dac4dc179b0c1d0340fd2fe826d4b3c9b4ae93387",
            `results/transcriptomics/gene_qc/trx_raw_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_raw_rle.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_ruvs_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_ruvs_rle.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_tmm_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_tmm_rle.png` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_optimization_summary.xlsx` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_panels.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_panels.png` = "<present>",
            `results/transcriptomics/publication_graphs/normalization_comparison_grid.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/normalization_comparison_grid.png` = "<present>",
            `scripts/transcriptomics/design_matrix.tab` = "cae70e69474a64ceed958d2428432b57c223ca71f81dfbbf9b9570c684295481",
            `scripts/transcriptomics/import_sources/counts.tsv` = "2c9e3b6290daf902297e7dab06425b831bfeb6261bd5b5d1819e6fd3e21dde49"),
        study_parameters = "<absent>",
        integration = stats::setNames(character(0), character()),
        session_state = "<0 files>",
        errors = "<unhandled> unused arguments (self, name)"),
    da = list(
        outputs = c(
            "report_ready", "ruviii_before_after", "ruviii_m_ui", "ruviii_status", "session_summary",
            "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            ruviii_before_after = "<table> 62875e4feeff8105b1284f0e99e836668bf7b7b93e25535e3c553a92a7f06583",
            ruviii_m_ui = "<ui> aafb6bcfd2a27f3e561956222696cdc9fd29de5f9954d8d5cb1e2c4b4bbd0ce9",
            ruviii_status = "Not run yet.",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = c(
            "[default NULL da_running] Running differential expression analysis...", "remove da_running",
            "[message 5 -] Differential expression analysis complete!",
            "[message 5 -] DA results saved to disk (tables, volcano plots, heatmaps)"),
        updates = c(
            "de-volcano_contrast: Sera_vs_RPMI -> Sera_vs_RPMI",
            "de-heatmap_contrast: Sera_vs_RPMI -> Sera_vs_RPMI",
            "de-table_contrast: Sera_vs_RPMI -> Sera_vs_RPMI", "de-volcano_assay: Combined|counts -> counts",
            "de-heatmap_assay: Combined|counts -> Combined", "de-table_assay: All -> All"),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = c(
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "28",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "6",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx-da-fit-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_da_contrast-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_da_contrast-da_contrast-<ID>` = "1",
            `results/transcriptomics/da_genes/de_results_Sera_vs_RPMI.csv` = "c11158b3af965955f4284fad607de9462097852a94c2b074ffab0e78c0a4af6b",
            `results/transcriptomics/da_genes/ruv_optimization_summary.txt` = "<present>",
            `results/transcriptomics/da_genes/top_de_genes.tsv` = "10aa0f7b3fa91cdc9137a9e4f11dca3a1413f4bd4e236542637c471bf87daa21",
            `results/transcriptomics/publication_graphs/NumSigDeMolecules/summary_barplot_Sera_vs_RPMI.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/NumSigDeMolecules/summary_barplot_Sera_vs_RPMI.png` = "<present>",
            `results/transcriptomics/publication_graphs/Volcano_Plots/volcano_plot_Sera_vs_RPMI.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/Volcano_Plots/volcano_plot_Sera_vs_RPMI.png` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = c(
            "<unhandled> argument is of length zero", "<unhandled> argument is of length zero")),
    ruviii_technical = list(
        outputs = c(
            "report_ready", "ruviii_before_after", "ruviii_m_ui", "ruviii_status", "session_summary",
            "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            ruviii_before_after = "<plot> alt=RUVIII before and after",
            ruviii_m_ui = "<ui> aafb6bcfd2a27f3e561956222696cdc9fd29de5f9954d8d5cb1e2c4b4bbd0ce9",
            ruviii_status = "RUVIII export complete (integration only; DA results unchanged).\nM = technical_replicate\nk = 2, 18 inferred negative controls",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = c(
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "31",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_ruviii_export-<ID>` = "1",
            `results/transcriptomics/gene_qc/trx_ruviii_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_ruviii_rle.png` = "<present>",
            `results/transcriptomics/publication_graphs/ruviii_before_and_after.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/ruviii_before_and_after.png` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0),
        ruviii = c(
            m_source = "technical_replicate", k = "2", controls = "18", newY = "12x187 eac9dbacdef1380e")))

omicSummaryW1Frozen[["trx_managed"]] <- list(
    norm = list(
        outputs = c(
            "report_ready", "ruviii_before_after", "ruviii_m_ui", "ruviii_status", "session_summary",
            "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            ruviii_before_after = "<table> 62875e4feeff8105b1284f0e99e836668bf7b7b93e25535e3c553a92a7f06583",
            ruviii_m_ui = "<ui> e7ae5dc12fe193d9706149295d232a88b3b7a9e15fced503f10f12bc405b104f",
            ruviii_status = "Not run yet.",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = c(
            "[default NULL transcriptomics_action_working] Applying transcriptomics step ...",
            "remove transcriptomics_action_working",
            "[default NULL transcriptomics_action_working] Applying transcriptomics step ...",
            "remove transcriptomics_action_working",
            "[default NULL transcriptomics_action_working] Applying transcriptomics step ...",
            "remove transcriptomics_action_working",
            "[default NULL transcriptomics_action_working] Applying transcriptomics step ...",
            "remove transcriptomics_action_working",
            "[default NULL transcriptomics_action_working] Applying transcriptomics step ...",
            "remove transcriptomics_action_working",
            "[default NULL transcriptomics_action_working] Applying transcriptomics step ...",
            "remove transcriptomics_action_working",
            "[default NULL transcriptomics_action_working] Applying transcriptomics step ...",
            "remove transcriptomics_action_working",
            "[default NULL transcriptomics_action_working] Applying transcriptomics step ...",
            "remove transcriptomics_action_working",
            "[message 5 -] Normalization complete! Proceeding to DE analysis."),
        updates = "summary-experiment_label: <same> -> kit108",
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = c(
            "omicsManagedProduct(operation=trx_qc_figures, parameters=stage:raw)",
            "omicsManagedProduct(operation=trx_tmm_logcpm, parameters=)",
            "omicsManagedProduct(operation=trx_ruv_grid, parameters=logcpm:<entity>;ruv_mode:automatic;percentage_min:5;percentage_max:30;percentage_step:5;k_penalty_weight:0.05;manual_k:NA;manual_percentage:NA)",
            "omicsManagedProduct(operation=trx_qc_figures, parameters=stage:tmm)",
            "omicsManagedProduct(operation=trx_qc_figures, parameters=stage:ruvs)",
            "omicsManagedProduct(operation=trx_ruv_publication, parameters=)",
            "omicsManagedProduct(operation=trx_composite_figure, parameters=figure:normalization_comparison)"),
        global = c(before = "absent", after = "absent"),
        files = c(
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "36",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "5",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/partitioned-table-<ID>` = "2",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_expression_filter-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_ruv_grid-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_ruvs-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_tmm_logcpm-<ID>` = "1",
            `results/transcriptomics/gene_qc/ruv_normalised_log_values.tsv` = "83de1dbdd1018161c5ec6c1dac4dc179b0c1d0340fd2fe826d4b3c9b4ae93387",
            `results/transcriptomics/gene_qc/trx_raw_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_raw_rle.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_ruvs_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_ruvs_rle.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_tmm_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_tmm_rle.png` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_optimization_summary.xlsx` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_panels.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_panels.png` = "<present>",
            `results/transcriptomics/publication_graphs/normalization_comparison_grid.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/normalization_comparison_grid.png` = "<present>",
            `scripts/transcriptomics/design_matrix.tab` = "e3b6f40550967e975506c582126b035f5fcc887aae4ffe97e3200f3999a3584e",
            `scripts/transcriptomics/import_sources/counts.tsv` = "2c9e3b6290daf902297e7dab06425b831bfeb6261bd5b5d1819e6fd3e21dde49"),
        study_parameters = "<absent>",
        integration = stats::setNames(character(0), character()),
        session_state = "<0 files>",
        errors = "<unhandled> unused arguments (self, name)"),
    da = list(
        outputs = c(
            "report_ready", "ruviii_before_after", "ruviii_m_ui", "ruviii_status", "session_summary",
            "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            ruviii_before_after = "<table> 62875e4feeff8105b1284f0e99e836668bf7b7b93e25535e3c553a92a7f06583",
            ruviii_m_ui = "<ui> e7ae5dc12fe193d9706149295d232a88b3b7a9e15fced503f10f12bc405b104f",
            ruviii_status = "Not run yet.",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = c(
            "[default NULL da_running] Running differential expression analysis...",
            "[default NULL transcriptomics_action_working] Applying transcriptomics step ...",
            "remove transcriptomics_action_working",
            "[default NULL transcriptomics_action_working] Applying transcriptomics step ...",
            "remove transcriptomics_action_working",
            "[default NULL transcriptomics_action_working] Applying transcriptomics step ...",
            "remove transcriptomics_action_working", "remove da_running",
            "[message 5 -] Differential expression analysis complete!",
            "[message 5 -] DA results saved to disk (tables, volcano plots, heatmaps)"),
        updates = c(
            "de-volcano_contrast: Sera_vs_RPMI -> Sera_vs_RPMI",
            "de-heatmap_contrast: Sera_vs_RPMI -> Sera_vs_RPMI",
            "de-table_contrast: Sera_vs_RPMI -> Sera_vs_RPMI", "de-volcano_assay: Combined|counts -> counts",
            "de-heatmap_assay: Combined|counts -> Combined", "de-table_assay: All -> All"),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = c(
            "omicsManagedProduct(operation=trx_da_contrast, parameters=formula:~ 0 + group + W_1 + W_2;qvalue:0.05;lfc:1;contrast:Sera_vs_RPMI=groupSera-groupRPMI;friendly_name:Sera_vs_RPMI)",
            "omicsManagedProduct(operation=trx_da_tab_results, parameters=)",
            "omicsManagedProduct(operation=trx_da_publication, parameters=)"),
        global = c(before = "absent", after = "absent"),
        files = c(
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "44",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "6",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx-da-fit-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_da_contrast-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_da_contrast-da_contrast-<ID>` = "1",
            `results/transcriptomics/da_genes/de_results_Sera_vs_RPMI.csv` = "c11158b3af965955f4284fad607de9462097852a94c2b074ffab0e78c0a4af6b",
            `results/transcriptomics/da_genes/ruv_optimization_summary.txt` = "<present>",
            `results/transcriptomics/da_genes/top_de_genes.tsv` = "10aa0f7b3fa91cdc9137a9e4f11dca3a1413f4bd4e236542637c471bf87daa21",
            `results/transcriptomics/publication_graphs/NumSigDeMolecules/summary_barplot_Sera_vs_RPMI.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/NumSigDeMolecules/summary_barplot_Sera_vs_RPMI.png` = "<present>",
            `results/transcriptomics/publication_graphs/Volcano_Plots/volcano_plot_Sera_vs_RPMI.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/Volcano_Plots/volcano_plot_Sera_vs_RPMI.png` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = c(
            "<unhandled> argument is of length zero", "<unhandled> argument is of length zero")),
    ruviii = list(
        outputs = c(
            "report_ready", "ruviii_before_after", "ruviii_m_ui", "ruviii_status", "session_summary",
            "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            ruviii_before_after = "<plot> alt=RUVIII before and after",
            ruviii_m_ui = "<ui> e7ae5dc12fe193d9706149295d232a88b3b7a9e15fced503f10f12bc405b104f",
            ruviii_status = "RUVIII export complete (integration only; DA results unchanged).\nM = group\nk = 2, 18 inferred negative controls",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = c(
            "[default NULL transcriptomics_action_working] Applying transcriptomics step ...",
            "remove transcriptomics_action_working",
            "[default NULL transcriptomics_action_working] Applying transcriptomics step ...",
            "remove transcriptomics_action_working",
            "[default NULL transcriptomics_action_working] Applying transcriptomics step ...",
            "remove transcriptomics_action_working"),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = c(
            "omicsManagedProduct(operation=trx_ruviii_export, parameters=logcpm:<entity>;grid:<entity>;ruviii_m_source:group)",
            "omicsManagedProduct(operation=trx_qc_figures, parameters=stage:ruviii)",
            "omicsManagedProduct(operation=trx_composite_figure, parameters=figure:ruviii_before_and_after)"),
        global = c(before = "absent", after = "absent"),
        files = c(
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "53",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_ruviii_export-<ID>` = "1",
            `results/transcriptomics/gene_qc/trx_ruviii_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_ruviii_rle.png` = "<present>",
            `results/transcriptomics/publication_graphs/ruviii_before_and_after.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/ruviii_before_and_after.png` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0),
        ruviii = c(m_source = "group", k = "2", controls = "18", newY = "12x187 6fc806b39ca98cbc")),
    save = list(
        outputs = c(
            "report_ready", "ruviii_before_after", "ruviii_m_ui", "ruviii_status", "session_summary",
            "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            ruviii_before_after = "<plot> alt=RUVIII before and after",
            ruviii_m_ui = "<ui> e7ae5dc12fe193d9706149295d232a88b3b7a9e15fced503f10f12bc405b104f",
            ruviii_status = "RUVIII export complete (integration only; DA results unchanged).\nM = group\nk = 2, 18 inferred negative controls",
            session_summary = "Study parameters created for: kit108 \nDescription: Full transcriptomics analysis workflow with normalization, DE and the RUVIII integration export \nTimestamp: <TIME> \nFile: <ROOT>/scripts/transcriptomics/study_parameters.txt \nSource: Final S4 object @args + config_list \nIntegration Object: transcriptomics_normalized_data.tsv \nStatus: Parameters saved [OK]",
            template_status = "[WARNING] Report template will be downloaded when generating report"),
        notes = c(
            "[default NULL transcriptomics_action_working] Applying transcriptomics step ...",
            "remove transcriptomics_action_working", "[message 5 -] Saved Integration S4 Object",
            "[message 5 -] Study parameters saved successfully"),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = c(
            "createWorkflowArgsFromConfig(workflow_name=kit108, description=Full transcriptomics analysis workflow with normalization, DE and the RUVIII integration export, organism_name=<NULL>, taxon_id=<NULL>, source_dir_path=<ROOT>/scripts/transcriptomics, final_s4_object=<ref>, contrasts_tbl=<table 1x3 ba384cec7aecd593>, workflow_data=<environment>)",
            "writeWorkflowAnalysisRecord(manager=<WorkflowStateOwner>, path=<ROOT>/scripts/transcriptomics/analysis_record.json)",
            "omicsManagedProduct(operation=trx_integration_export, parameters=path:<ROOT>/integration/transcriptomics_normalized_data.tsv;state_name:trx_ruvs)"),
        global = c(before = "absent", after = "absent"),
        files = c(
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "54",
            `integration/transcriptomics_normalized_data.tsv` = "7b9d35dbf33331d16b6609caacbd9445879924bd8bb548b6821dfd649bfdd08e",
            `scripts/transcriptomics/analysis_record.json` = "<present>",
            `scripts/transcriptomics/study_parameters.txt` = "<present>"),
        study_parameters = "57 lines 6f90092f3d590fd4",
        integration = c(transcriptomics_normalized_data.tsv = "7b9d35dbf33331d1"),
        session_state = "<unchanged>",
        errors = character(0)))

omicSummaryW1Frozen[["prot_refusals"]] <- list(
    start = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = c(
            `da_output/da_groupB-groupA_long_annot.tsv` = "1ece3b87eeeedd82d2783415db679309b23b1cd9378d8d32e671989a020abb79",
            `da_output/da_groupB-groupA_long_annot.xlsx` = "<present>",
            `feature_qc/composite_QC_figure.pdf` = "<present>",
            `feature_qc/composite_QC_figure.png` = "<present>",
            `feature_qc/ruv_normalised_results_cln_with_replicates.RDS` = "<present>",
            `feature_qc/ruv_normalised_results_cln_with_replicates.tsv` = "46d9ba7fe967c79d88fad400b5998407771a19914a2d6f0997b8820e6a5074b5",
            `pathway_enrichment/B_vs_A_down_enrichment_results.tsv` = "4e168310ff617aba595f1605046a7d4782d44f232692d559db1bf56bf3e71568",
            `pathway_enrichment/B_vs_A_up_enrichment_results.tsv` = "b56e92aa35bf8f4fbd3dd96a2e0de7f4ceb9235d00eef687ce41204bf398ad68",
            `publication_graphs/Heatmap/artifact.txt` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/artifact.txt` = "<present>",
            `publication_graphs/Volcano_Plots/artifact.txt` = "<present>",
            `source/contrasts_tbl.tab` = "7e2814a29d730c668dca584b855fe542882c63a393a5c4deaba6a98f42d5d656",
            `source/design_matrix.tab` = "8812863972c9c16304791cc9cf303a962eaef1b1cb0a8fd193b98acfd414303d",
            `source/study_parameters.txt` = "<present>",
            `time/12_correlation_filtered_combined_plots.png` = "<present>"),
        study_parameters = "1 lines 55c1b3aadde3f3b7",
        integration = stats::setNames(character(0), character()),
        session_state = "<0 files>",
        errors = "<unhandled> unused arguments (self, name)"),
    save_no_label = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    copy_no_label = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    report_no_label = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = "[warning 5 -] Please enter an experiment label before generating the report.",
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    export_no_label = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    report_before_copy = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = "[warning 5 -] Please copy the results to the publication directory before generating the report (Copy to Publication Directory).",
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    github_before_report = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)))

omicSummaryW1Frozen[["prot_route_dia"]] <- list(
    start = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = character(0),
        updates = "session_summary-experiment_label: <same> -> kit108",
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = c(
            `da_output/da_groupB-groupA_long_annot.tsv` = "1ece3b87eeeedd82d2783415db679309b23b1cd9378d8d32e671989a020abb79",
            `da_output/da_groupB-groupA_long_annot.xlsx` = "<present>",
            `feature_qc/composite_QC_figure.pdf` = "<present>",
            `feature_qc/composite_QC_figure.png` = "<present>",
            `feature_qc/ruv_normalised_results_cln_with_replicates.RDS` = "<present>",
            `feature_qc/ruv_normalised_results_cln_with_replicates.tsv` = "46d9ba7fe967c79d88fad400b5998407771a19914a2d6f0997b8820e6a5074b5",
            `pathway_enrichment/B_vs_A_down_enrichment_results.tsv` = "4e168310ff617aba595f1605046a7d4782d44f232692d559db1bf56bf3e71568",
            `pathway_enrichment/B_vs_A_up_enrichment_results.tsv` = "b56e92aa35bf8f4fbd3dd96a2e0de7f4ceb9235d00eef687ce41204bf398ad68",
            `publication_graphs/Heatmap/artifact.txt` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/artifact.txt` = "<present>",
            `publication_graphs/Volcano_Plots/artifact.txt` = "<present>",
            `source/contrasts_tbl.tab` = "7e2814a29d730c668dca584b855fe542882c63a393a5c4deaba6a98f42d5d656",
            `source/design_matrix.tab` = "8812863972c9c16304791cc9cf303a962eaef1b1cb0a8fd193b98acfd414303d",
            `source/study_parameters.txt` = "<present>",
            `time/12_correlation_filtered_combined_plots.png` = "<present>"),
        study_parameters = "1 lines 55c1b3aadde3f3b7",
        integration = stats::setNames(character(0), character()),
        session_state = "<0 files>",
        errors = "<unhandled> unused arguments (self, name)"),
    save = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Study parameters created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nFile: <ROOT>/source/study_parameters.txt \nSource: Final S4 object @args + config_list \nIntegration Object: proteomics_kit108_final_s4.RDS \nStatus: Parameters saved [OK]",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = c(
            "[message 5 -] Saved Integration S4 Object", "[message 5 -] Study parameters saved successfully"),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = "createWorkflowArgsFromConfig(workflow_name=kit108, description=Full protein analysis workflow with config parameters, organism_name=<NULL>, taxon_id=<NULL>, source_dir_path=<ROOT>/source, final_s4_object=<NULL>, contrasts_tbl=<table 1x3 950ddbb7d830df8c>, workflow_data=<environment>)",
        global = c(before = "absent", after = "absent"),
        files = c(`integration/proteomics_kit108_final_s4.RDS` = "<present>"),
        study_parameters = "60 lines 7a948751f37e3fd7",
        integration = c(proteomics_kit108_final_s4.RDS = "<RDS> NULL"),
        session_state = "<unchanged>",
        errors = character(0)),
    copy = list(
        outputs = c("copy_status", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK]",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = "[message 5 -] Publication files copied",
        updates = character(0),
        progress = "start Copying files to publication directory... 0.1",
        reactive_values = c(
            files_copied = "TRUE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = c(
            "copyToResultsSummary(omic_type=proteomics, experiment_label=kit108, contrasts_tbl=<table 1x3 950ddbb7d830df8c>, design_matrix=<table 4x3 bcd78aff6b352869>, force=TRUE, current_rmd=<NULL>, project_dirs_object_name=project_dirs, global=the session's)",
            "  copyToResultsSummary returned: "),
        global = c(before = "absent", after = "the session's"),
        files = c(
            `results_summary/Publication_figures/Enrichment_Plots/B_vs_A_down_enrichment_results.tsv` = "4e168310ff617aba595f1605046a7d4782d44f232692d559db1bf56bf3e71568",
            `results_summary/Publication_figures/Enrichment_Plots/B_vs_A_up_enrichment_results.tsv` = "b56e92aa35bf8f4fbd3dd96a2e0de7f4ceb9235d00eef687ce41204bf398ad68",
            `results_summary/Publication_figures/Heatmap/artifact.txt` = "<present>",
            `results_summary/Publication_figures/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `results_summary/Publication_figures/NumSigDaMolecules/artifact.txt` = "<present>",
            `results_summary/Publication_figures/Volcano_Plots/artifact.txt` = "<present>",
            `results_summary/Publication_tables/DA_results_proteomics.xlsx` = "<present>",
            `results_summary/Publication_tables/Pathway_enrichment_results_proteomics.xlsx` = "<present>",
            `results_summary/Publication_tables/RUV_normalised_results.tsv` = "46d9ba7fe967c79d88fad400b5998407771a19914a2d6f0997b8820e6a5074b5",
            `results_summary/Publication_tables/ruv_normalised_results.RDS` = "<present>",
            `results_summary/QC_figures/12_correlation_filtered_combined_plots.png` = "<present>",
            `results_summary/QC_figures/proteomics_composite_QC_figure.pdf` = "<present>",
            `results_summary/QC_figures/proteomics_composite_QC_figure.png` = "<present>",
            `results_summary/Study_report/contrasts_tbl.tab` = "7e2814a29d730c668dca584b855fe542882c63a393a5c4deaba6a98f42d5d656",
            `results_summary/Study_report/design_matrix.tab` = "8812863972c9c16304791cc9cf303a962eaef1b1cb0a8fd193b98acfd414303d",
            `results_summary/Study_report/study_parameters.txt` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    report = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK] \nReport location: <ROOT>/results_summary/DIANN_report_proteomics_kit108.docx",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = c(
            "[message 5 -] DIANN_report.rmd template copied from package",
            "[message 5 -] Report generated successfully!"),
        updates = character(0),
        progress = c(
            "start Generating report... 0.1", "0.1 Detecting workflow type...",
            "0.2 Checking for DIANN_report.rmd template...", "0.5 Rendering report..."),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/DIANN_report_proteomics_kit108.docx",
            workflow_args_saved = "TRUE"),
        calls = "RenderReport(omic_type=proteomics, experiment_label=kit108, rmd_filename=DIANN_report.rmd, project_dirs_object_name=project_dirs, output_format=<NULL>, global=the session's)",
        global = c(before = "the session's", after = "the session's"),
        files = c(
            `results_summary/DIANN_report_proteomics_kit108.docx` = "<present>",
            `scripts/proteomics/DIANN_report.rmd` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = "<unhandled> unused arguments (self, name)"),
    export = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK] \nReport location: <ROOT>/results_summary/DIANN_report_proteomics_kit108.docx",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = "[message 5 -] Session state exported to: <ROOT>/source/session_state_<DATE>.RDS",
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/DIANN_report_proteomics_kit108.docx",
            workflow_args_saved = "TRUE"),
        calls = character(0),
        global = c(before = "the session's", after = "the session's"),
        files = c(`source/session_state_<DATE>.RDS` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = c(
            experiment_label = "kit108",
            description = "Full protein analysis workflow with config parameters", timestamp = "<time>",
            omic_type = "proteomics", workflow_args_saved = "TRUE", files_copied = "TRUE",
            report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/DIANN_report_proteomics_kit108.docx", workflow_type = "DIA",
            report_template = "DIANN_report.rmd", parameter_payload = "<digest> bef1e766f17f28fb",
            project_dirs = "<digest> 63ed7c866939a783"),
        errors = character(0)))

omicSummaryW1Frozen[["prot_route_dia_limpa"]] <- list(
    start = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = character(0),
        updates = "session_summary-experiment_label: <same> -> kit108",
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = c(
            `da_output/da_groupB-groupA_long_annot.tsv` = "1ece3b87eeeedd82d2783415db679309b23b1cd9378d8d32e671989a020abb79",
            `da_output/da_groupB-groupA_long_annot.xlsx` = "<present>",
            `feature_qc/composite_QC_figure.pdf` = "<present>",
            `feature_qc/composite_QC_figure.png` = "<present>",
            `feature_qc/ruv_normalised_results_cln_with_replicates.RDS` = "<present>",
            `feature_qc/ruv_normalised_results_cln_with_replicates.tsv` = "46d9ba7fe967c79d88fad400b5998407771a19914a2d6f0997b8820e6a5074b5",
            `pathway_enrichment/B_vs_A_down_enrichment_results.tsv` = "4e168310ff617aba595f1605046a7d4782d44f232692d559db1bf56bf3e71568",
            `pathway_enrichment/B_vs_A_up_enrichment_results.tsv` = "b56e92aa35bf8f4fbd3dd96a2e0de7f4ceb9235d00eef687ce41204bf398ad68",
            `publication_graphs/Heatmap/artifact.txt` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/artifact.txt` = "<present>",
            `publication_graphs/Volcano_Plots/artifact.txt` = "<present>",
            `source/contrasts_tbl.tab` = "7e2814a29d730c668dca584b855fe542882c63a393a5c4deaba6a98f42d5d656",
            `source/design_matrix.tab` = "8812863972c9c16304791cc9cf303a962eaef1b1cb0a8fd193b98acfd414303d",
            `source/study_parameters.txt` = "<present>",
            `time/12_correlation_filtered_combined_plots.png` = "<present>"),
        study_parameters = "1 lines 55c1b3aadde3f3b7",
        integration = stats::setNames(character(0), character()),
        session_state = "<0 files>",
        errors = "<unhandled> unused arguments (self, name)"),
    save = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Study parameters created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nFile: <ROOT>/source/study_parameters.txt \nSource: Final S4 object @args + config_list \nIntegration Object: proteomics_kit108_final_s4.RDS \nStatus: Parameters saved [OK]",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = c(
            "[message 5 -] Saved Integration S4 Object", "[message 5 -] Study parameters saved successfully"),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = "createWorkflowArgsFromConfig(workflow_name=kit108, description=Full protein analysis workflow with config parameters, organism_name=<NULL>, taxon_id=<NULL>, source_dir_path=<ROOT>/source, final_s4_object=<NULL>, contrasts_tbl=<table 1x3 950ddbb7d830df8c>, workflow_data=<environment>)",
        global = c(before = "absent", after = "absent"),
        files = c(`integration/proteomics_kit108_final_s4.RDS` = "<present>"),
        study_parameters = "60 lines 3d1983f1d34802ae",
        integration = c(proteomics_kit108_final_s4.RDS = "<RDS> NULL"),
        session_state = "<unchanged>",
        errors = character(0)),
    copy = list(
        outputs = c("copy_status", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK]",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = "[message 5 -] Publication files copied",
        updates = character(0),
        progress = "start Copying files to publication directory... 0.1",
        reactive_values = c(
            files_copied = "TRUE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = c(
            "copyToResultsSummary(omic_type=proteomics, experiment_label=kit108, contrasts_tbl=<table 1x3 950ddbb7d830df8c>, design_matrix=<table 4x3 bcd78aff6b352869>, force=TRUE, current_rmd=<NULL>, project_dirs_object_name=project_dirs, global=the session's)",
            "  copyToResultsSummary returned: "),
        global = c(before = "absent", after = "the session's"),
        files = c(
            `results_summary/Publication_figures/Enrichment_Plots/B_vs_A_down_enrichment_results.tsv` = "4e168310ff617aba595f1605046a7d4782d44f232692d559db1bf56bf3e71568",
            `results_summary/Publication_figures/Enrichment_Plots/B_vs_A_up_enrichment_results.tsv` = "b56e92aa35bf8f4fbd3dd96a2e0de7f4ceb9235d00eef687ce41204bf398ad68",
            `results_summary/Publication_figures/Heatmap/artifact.txt` = "<present>",
            `results_summary/Publication_figures/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `results_summary/Publication_figures/NumSigDaMolecules/artifact.txt` = "<present>",
            `results_summary/Publication_figures/Volcano_Plots/artifact.txt` = "<present>",
            `results_summary/Publication_tables/DA_results_proteomics.xlsx` = "<present>",
            `results_summary/Publication_tables/Pathway_enrichment_results_proteomics.xlsx` = "<present>",
            `results_summary/Publication_tables/RUV_normalised_results.tsv` = "46d9ba7fe967c79d88fad400b5998407771a19914a2d6f0997b8820e6a5074b5",
            `results_summary/Publication_tables/ruv_normalised_results.RDS` = "<present>",
            `results_summary/QC_figures/12_correlation_filtered_combined_plots.png` = "<present>",
            `results_summary/QC_figures/proteomics_composite_QC_figure.pdf` = "<present>",
            `results_summary/QC_figures/proteomics_composite_QC_figure.png` = "<present>",
            `results_summary/Study_report/contrasts_tbl.tab` = "7e2814a29d730c668dca584b855fe542882c63a393a5c4deaba6a98f42d5d656",
            `results_summary/Study_report/design_matrix.tab` = "8812863972c9c16304791cc9cf303a962eaef1b1cb0a8fd193b98acfd414303d",
            `results_summary/Study_report/study_parameters.txt` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    report = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK] \nReport location: <ROOT>/results_summary/DIANN_limpa_report_proteomics_kit108.docx",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = c(
            "[message 5 -] DIANN_limpa_report.rmd template copied from package",
            "[message 5 -] Report generated successfully!"),
        updates = character(0),
        progress = c(
            "start Generating report... 0.1", "0.1 Detecting workflow type...",
            "0.2 Checking for DIANN_limpa_report.rmd template...", "0.5 Rendering report..."),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/DIANN_limpa_report_proteomics_kit108.docx",
            workflow_args_saved = "TRUE"),
        calls = "RenderReport(omic_type=proteomics, experiment_label=kit108, rmd_filename=DIANN_limpa_report.rmd, project_dirs_object_name=project_dirs, output_format=<NULL>, global=the session's)",
        global = c(before = "the session's", after = "the session's"),
        files = c(
            `results_summary/DIANN_limpa_report_proteomics_kit108.docx` = "<present>",
            `scripts/proteomics/DIANN_limpa_report.rmd` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = "<unhandled> unused arguments (self, name)"),
    export = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK] \nReport location: <ROOT>/results_summary/DIANN_limpa_report_proteomics_kit108.docx",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = "[message 5 -] Session state exported to: <ROOT>/source/session_state_<DATE>.RDS",
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/DIANN_limpa_report_proteomics_kit108.docx",
            workflow_args_saved = "TRUE"),
        calls = character(0),
        global = c(before = "the session's", after = "the session's"),
        files = c(`source/session_state_<DATE>.RDS` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = c(
            experiment_label = "kit108",
            description = "Full protein analysis workflow with config parameters", timestamp = "<time>",
            omic_type = "proteomics", workflow_args_saved = "TRUE", files_copied = "TRUE",
            report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/DIANN_limpa_report_proteomics_kit108.docx",
            workflow_type = "DIA_LIMPA", report_template = "DIANN_limpa_report.rmd",
            parameter_payload = "<digest> 6b293857c6cae440", project_dirs = "<digest> 63ed7c866939a783"),
        errors = character(0)))

omicSummaryW1Frozen[["prot_route_tmt"]] <- list(
    start = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = character(0),
        updates = "session_summary-experiment_label: <same> -> kit108",
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = c(
            `da_output/da_groupB-groupA_long_annot.tsv` = "1ece3b87eeeedd82d2783415db679309b23b1cd9378d8d32e671989a020abb79",
            `da_output/da_groupB-groupA_long_annot.xlsx` = "<present>",
            `feature_qc/composite_QC_figure.pdf` = "<present>",
            `feature_qc/composite_QC_figure.png` = "<present>",
            `feature_qc/ruv_normalised_results_cln_with_replicates.RDS` = "<present>",
            `feature_qc/ruv_normalised_results_cln_with_replicates.tsv` = "46d9ba7fe967c79d88fad400b5998407771a19914a2d6f0997b8820e6a5074b5",
            `pathway_enrichment/B_vs_A_down_enrichment_results.tsv` = "4e168310ff617aba595f1605046a7d4782d44f232692d559db1bf56bf3e71568",
            `pathway_enrichment/B_vs_A_up_enrichment_results.tsv` = "b56e92aa35bf8f4fbd3dd96a2e0de7f4ceb9235d00eef687ce41204bf398ad68",
            `publication_graphs/Heatmap/artifact.txt` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/artifact.txt` = "<present>",
            `publication_graphs/Volcano_Plots/artifact.txt` = "<present>",
            `source/contrasts_tbl.tab` = "7e2814a29d730c668dca584b855fe542882c63a393a5c4deaba6a98f42d5d656",
            `source/design_matrix.tab` = "8812863972c9c16304791cc9cf303a962eaef1b1cb0a8fd193b98acfd414303d",
            `source/study_parameters.txt` = "<present>",
            `time/12_correlation_filtered_combined_plots.png` = "<present>"),
        study_parameters = "1 lines 55c1b3aadde3f3b7",
        integration = stats::setNames(character(0), character()),
        session_state = "<0 files>",
        errors = "<unhandled> unused arguments (self, name)"),
    save = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Study parameters created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nFile: <ROOT>/source/study_parameters.txt \nSource: Final S4 object @args + config_list \nIntegration Object: proteomics_kit108_final_s4.RDS \nStatus: Parameters saved [OK]",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = c(
            "[message 5 -] Saved Integration S4 Object", "[message 5 -] Study parameters saved successfully"),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = "createWorkflowArgsFromConfig(workflow_name=kit108, description=Full protein analysis workflow with config parameters, organism_name=<NULL>, taxon_id=<NULL>, source_dir_path=<ROOT>/source, final_s4_object=<NULL>, contrasts_tbl=<table 1x3 950ddbb7d830df8c>, workflow_data=<environment>)",
        global = c(before = "absent", after = "absent"),
        files = c(`integration/proteomics_kit108_final_s4.RDS` = "<present>"),
        study_parameters = "60 lines fe36eac93ca1aa4e",
        integration = c(proteomics_kit108_final_s4.RDS = "<RDS> NULL"),
        session_state = "<unchanged>",
        errors = character(0)),
    copy = list(
        outputs = c("copy_status", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK]",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = "[message 5 -] Publication files copied",
        updates = character(0),
        progress = "start Copying files to publication directory... 0.1",
        reactive_values = c(
            files_copied = "TRUE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = c(
            "copyToResultsSummary(omic_type=proteomics, experiment_label=kit108, contrasts_tbl=<table 1x3 950ddbb7d830df8c>, design_matrix=<table 4x3 bcd78aff6b352869>, force=TRUE, current_rmd=<NULL>, project_dirs_object_name=project_dirs, global=the session's)",
            "  copyToResultsSummary returned: "),
        global = c(before = "absent", after = "the session's"),
        files = c(
            `results_summary/Publication_figures/Enrichment_Plots/B_vs_A_down_enrichment_results.tsv` = "4e168310ff617aba595f1605046a7d4782d44f232692d559db1bf56bf3e71568",
            `results_summary/Publication_figures/Enrichment_Plots/B_vs_A_up_enrichment_results.tsv` = "b56e92aa35bf8f4fbd3dd96a2e0de7f4ceb9235d00eef687ce41204bf398ad68",
            `results_summary/Publication_figures/Heatmap/artifact.txt` = "<present>",
            `results_summary/Publication_figures/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `results_summary/Publication_figures/NumSigDaMolecules/artifact.txt` = "<present>",
            `results_summary/Publication_figures/Volcano_Plots/artifact.txt` = "<present>",
            `results_summary/Publication_tables/DA_results_proteomics.xlsx` = "<present>",
            `results_summary/Publication_tables/Pathway_enrichment_results_proteomics.xlsx` = "<present>",
            `results_summary/Publication_tables/RUV_normalised_results.tsv` = "46d9ba7fe967c79d88fad400b5998407771a19914a2d6f0997b8820e6a5074b5",
            `results_summary/Publication_tables/ruv_normalised_results.RDS` = "<present>",
            `results_summary/QC_figures/12_correlation_filtered_combined_plots.png` = "<present>",
            `results_summary/QC_figures/proteomics_composite_QC_figure.pdf` = "<present>",
            `results_summary/QC_figures/proteomics_composite_QC_figure.png` = "<present>",
            `results_summary/Study_report/contrasts_tbl.tab` = "7e2814a29d730c668dca584b855fe542882c63a393a5c4deaba6a98f42d5d656",
            `results_summary/Study_report/design_matrix.tab` = "8812863972c9c16304791cc9cf303a962eaef1b1cb0a8fd193b98acfd414303d",
            `results_summary/Study_report/study_parameters.txt` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    report = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK] \nReport location: <ROOT>/results_summary/TMT_report_proteomics_kit108.docx",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = c(
            "[message 5 -] TMT_report.rmd template copied from package",
            "[message 5 -] Report generated successfully!"),
        updates = character(0),
        progress = c(
            "start Generating report... 0.1", "0.1 Detecting workflow type...",
            "0.2 Checking for TMT_report.rmd template...", "0.5 Rendering report..."),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/TMT_report_proteomics_kit108.docx",
            workflow_args_saved = "TRUE"),
        calls = "RenderReport(omic_type=proteomics, experiment_label=kit108, rmd_filename=TMT_report.rmd, project_dirs_object_name=project_dirs, output_format=<NULL>, global=the session's)",
        global = c(before = "the session's", after = "the session's"),
        files = c(
            `results_summary/TMT_report_proteomics_kit108.docx` = "<present>",
            `scripts/proteomics/TMT_report.rmd` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = "<unhandled> unused arguments (self, name)"),
    export = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK] \nReport location: <ROOT>/results_summary/TMT_report_proteomics_kit108.docx",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = "[message 5 -] Session state exported to: <ROOT>/source/session_state_<DATE>.RDS",
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/TMT_report_proteomics_kit108.docx",
            workflow_args_saved = "TRUE"),
        calls = character(0),
        global = c(before = "the session's", after = "the session's"),
        files = c(`source/session_state_<DATE>.RDS` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = c(
            experiment_label = "kit108",
            description = "Full protein analysis workflow with config parameters", timestamp = "<time>",
            omic_type = "proteomics", workflow_args_saved = "TRUE", files_copied = "TRUE",
            report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/TMT_report_proteomics_kit108.docx",
            workflow_type = "TMT_PD", report_template = "TMT_report.rmd",
            parameter_payload = "<digest> d676c63741a9b9fc", project_dirs = "<digest> 63ed7c866939a783"),
        errors = character(0)))

omicSummaryW1Frozen[["prot_route_lfq"]] <- list(
    start = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = character(0),
        updates = "session_summary-experiment_label: <same> -> kit108",
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = c(
            `da_output/da_groupB-groupA_long_annot.tsv` = "1ece3b87eeeedd82d2783415db679309b23b1cd9378d8d32e671989a020abb79",
            `da_output/da_groupB-groupA_long_annot.xlsx` = "<present>",
            `feature_qc/composite_QC_figure.pdf` = "<present>",
            `feature_qc/composite_QC_figure.png` = "<present>",
            `feature_qc/ruv_normalised_results_cln_with_replicates.RDS` = "<present>",
            `feature_qc/ruv_normalised_results_cln_with_replicates.tsv` = "46d9ba7fe967c79d88fad400b5998407771a19914a2d6f0997b8820e6a5074b5",
            `pathway_enrichment/B_vs_A_down_enrichment_results.tsv` = "4e168310ff617aba595f1605046a7d4782d44f232692d559db1bf56bf3e71568",
            `pathway_enrichment/B_vs_A_up_enrichment_results.tsv` = "b56e92aa35bf8f4fbd3dd96a2e0de7f4ceb9235d00eef687ce41204bf398ad68",
            `publication_graphs/Heatmap/artifact.txt` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/artifact.txt` = "<present>",
            `publication_graphs/Volcano_Plots/artifact.txt` = "<present>",
            `source/contrasts_tbl.tab` = "7e2814a29d730c668dca584b855fe542882c63a393a5c4deaba6a98f42d5d656",
            `source/design_matrix.tab` = "8812863972c9c16304791cc9cf303a962eaef1b1cb0a8fd193b98acfd414303d",
            `source/study_parameters.txt` = "<present>",
            `time/12_correlation_filtered_combined_plots.png` = "<present>"),
        study_parameters = "1 lines 55c1b3aadde3f3b7",
        integration = stats::setNames(character(0), character()),
        session_state = "<0 files>",
        errors = "<unhandled> unused arguments (self, name)"),
    save = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Study parameters created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nFile: <ROOT>/source/study_parameters.txt \nSource: Final S4 object @args + config_list \nIntegration Object: proteomics_kit108_final_s4.RDS \nStatus: Parameters saved [OK]",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = c(
            "[message 5 -] Saved Integration S4 Object", "[message 5 -] Study parameters saved successfully"),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = "createWorkflowArgsFromConfig(workflow_name=kit108, description=Full protein analysis workflow with config parameters, organism_name=<NULL>, taxon_id=<NULL>, source_dir_path=<ROOT>/source, final_s4_object=<NULL>, contrasts_tbl=<table 1x3 950ddbb7d830df8c>, workflow_data=<environment>)",
        global = c(before = "absent", after = "absent"),
        files = c(`integration/proteomics_kit108_final_s4.RDS` = "<present>"),
        study_parameters = "60 lines 3497f04a3e6b5a31",
        integration = c(proteomics_kit108_final_s4.RDS = "<RDS> NULL"),
        session_state = "<unchanged>",
        errors = character(0)),
    copy = list(
        outputs = c("copy_status", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK]",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = "[message 5 -] Publication files copied",
        updates = character(0),
        progress = "start Copying files to publication directory... 0.1",
        reactive_values = c(
            files_copied = "TRUE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = c(
            "copyToResultsSummary(omic_type=proteomics, experiment_label=kit108, contrasts_tbl=<table 1x3 950ddbb7d830df8c>, design_matrix=<table 4x3 bcd78aff6b352869>, force=TRUE, current_rmd=<NULL>, project_dirs_object_name=project_dirs, global=the session's)",
            "  copyToResultsSummary returned: "),
        global = c(before = "absent", after = "the session's"),
        files = c(
            `results_summary/Publication_figures/Enrichment_Plots/B_vs_A_down_enrichment_results.tsv` = "4e168310ff617aba595f1605046a7d4782d44f232692d559db1bf56bf3e71568",
            `results_summary/Publication_figures/Enrichment_Plots/B_vs_A_up_enrichment_results.tsv` = "b56e92aa35bf8f4fbd3dd96a2e0de7f4ceb9235d00eef687ce41204bf398ad68",
            `results_summary/Publication_figures/Heatmap/artifact.txt` = "<present>",
            `results_summary/Publication_figures/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `results_summary/Publication_figures/NumSigDaMolecules/artifact.txt` = "<present>",
            `results_summary/Publication_figures/Volcano_Plots/artifact.txt` = "<present>",
            `results_summary/Publication_tables/DA_results_proteomics.xlsx` = "<present>",
            `results_summary/Publication_tables/Pathway_enrichment_results_proteomics.xlsx` = "<present>",
            `results_summary/Publication_tables/RUV_normalised_results.tsv` = "46d9ba7fe967c79d88fad400b5998407771a19914a2d6f0997b8820e6a5074b5",
            `results_summary/Publication_tables/ruv_normalised_results.RDS` = "<present>",
            `results_summary/QC_figures/12_correlation_filtered_combined_plots.png` = "<present>",
            `results_summary/QC_figures/proteomics_composite_QC_figure.pdf` = "<present>",
            `results_summary/QC_figures/proteomics_composite_QC_figure.png` = "<present>",
            `results_summary/Study_report/contrasts_tbl.tab` = "7e2814a29d730c668dca584b855fe542882c63a393a5c4deaba6a98f42d5d656",
            `results_summary/Study_report/design_matrix.tab` = "8812863972c9c16304791cc9cf303a962eaef1b1cb0a8fd193b98acfd414303d",
            `results_summary/Study_report/study_parameters.txt` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    report = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK] \nReport location: <ROOT>/results_summary/LFQ_report_proteomics_kit108.docx",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = c(
            "[message 5 -] LFQ_report.rmd template copied from package",
            "[message 5 -] Report generated successfully!"),
        updates = character(0),
        progress = c(
            "start Generating report... 0.1", "0.1 Detecting workflow type...",
            "0.2 Checking for LFQ_report.rmd template...", "0.5 Rendering report..."),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/LFQ_report_proteomics_kit108.docx",
            workflow_args_saved = "TRUE"),
        calls = "RenderReport(omic_type=proteomics, experiment_label=kit108, rmd_filename=LFQ_report.rmd, project_dirs_object_name=project_dirs, output_format=<NULL>, global=the session's)",
        global = c(before = "the session's", after = "the session's"),
        files = c(
            `results_summary/LFQ_report_proteomics_kit108.docx` = "<present>",
            `scripts/proteomics/LFQ_report.rmd` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = "<unhandled> unused arguments (self, name)"),
    export = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK] \nReport location: <ROOT>/results_summary/LFQ_report_proteomics_kit108.docx",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = "[message 5 -] Session state exported to: <ROOT>/source/session_state_<DATE>.RDS",
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/LFQ_report_proteomics_kit108.docx",
            workflow_args_saved = "TRUE"),
        calls = character(0),
        global = c(before = "the session's", after = "the session's"),
        files = c(`source/session_state_<DATE>.RDS` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = c(
            experiment_label = "kit108",
            description = "Full protein analysis workflow with config parameters", timestamp = "<time>",
            omic_type = "proteomics", workflow_args_saved = "TRUE", files_copied = "TRUE",
            report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/LFQ_report_proteomics_kit108.docx", workflow_type = "LFQ",
            report_template = "LFQ_report.rmd", parameter_payload = "<digest> 69658edd5b0dd82d",
            project_dirs = "<digest> 63ed7c866939a783"),
        errors = character(0)))

omicSummaryW1Frozen[["prot_stale"]] <- list(
    copy = list(
        outputs = c("copy_status", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK]",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = "[message 5 -] Publication files copied",
        updates = "session_summary-experiment_label: <same> -> kit108",
        progress = "start Copying files to publication directory... 0.1",
        reactive_values = c(
            files_copied = "TRUE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = c(
            "copyToResultsSummary(omic_type=proteomics, experiment_label=kit108, contrasts_tbl=<table 1x3 950ddbb7d830df8c>, design_matrix=<table 4x3 bcd78aff6b352869>, force=TRUE, current_rmd=<NULL>, project_dirs_object_name=project_dirs, global=the session's)",
            "  copyToResultsSummary returned: "),
        global = c(before = "absent", after = "the session's"),
        files = c(
            `da_output/da_groupB-groupA_long_annot.tsv` = "1ece3b87eeeedd82d2783415db679309b23b1cd9378d8d32e671989a020abb79",
            `da_output/da_groupB-groupA_long_annot.xlsx` = "<present>",
            `feature_qc/composite_QC_figure.pdf` = "<present>",
            `feature_qc/composite_QC_figure.png` = "<present>",
            `feature_qc/ruv_normalised_results_cln_with_replicates.RDS` = "<present>",
            `feature_qc/ruv_normalised_results_cln_with_replicates.tsv` = "46d9ba7fe967c79d88fad400b5998407771a19914a2d6f0997b8820e6a5074b5",
            `pathway_enrichment/B_vs_A_down_enrichment_results.tsv` = "4e168310ff617aba595f1605046a7d4782d44f232692d559db1bf56bf3e71568",
            `pathway_enrichment/B_vs_A_up_enrichment_results.tsv` = "b56e92aa35bf8f4fbd3dd96a2e0de7f4ceb9235d00eef687ce41204bf398ad68",
            `publication_graphs/Heatmap/artifact.txt` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/artifact.txt` = "<present>",
            `publication_graphs/Volcano_Plots/artifact.txt` = "<present>",
            `results_summary/Publication_figures/Enrichment_Plots/B_vs_A_down_enrichment_results.tsv` = "4e168310ff617aba595f1605046a7d4782d44f232692d559db1bf56bf3e71568",
            `results_summary/Publication_figures/Enrichment_Plots/B_vs_A_up_enrichment_results.tsv` = "b56e92aa35bf8f4fbd3dd96a2e0de7f4ceb9235d00eef687ce41204bf398ad68",
            `results_summary/Publication_figures/Heatmap/artifact.txt` = "<present>",
            `results_summary/Publication_figures/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `results_summary/Publication_figures/NumSigDaMolecules/artifact.txt` = "<present>",
            `results_summary/Publication_figures/Volcano_Plots/artifact.txt` = "<present>",
            `results_summary/Publication_tables/DA_results_proteomics.xlsx` = "<present>",
            `results_summary/Publication_tables/Pathway_enrichment_results_proteomics.xlsx` = "<present>",
            `results_summary/Publication_tables/RUV_normalised_results.tsv` = "46d9ba7fe967c79d88fad400b5998407771a19914a2d6f0997b8820e6a5074b5",
            `results_summary/Publication_tables/ruv_normalised_results.RDS` = "<present>",
            `results_summary/QC_figures/12_correlation_filtered_combined_plots.png` = "<present>",
            `results_summary/QC_figures/proteomics_composite_QC_figure.pdf` = "<present>",
            `results_summary/QC_figures/proteomics_composite_QC_figure.png` = "<present>",
            `results_summary/Study_report/contrasts_tbl.tab` = "7e2814a29d730c668dca584b855fe542882c63a393a5c4deaba6a98f42d5d656",
            `results_summary/Study_report/design_matrix.tab` = "8812863972c9c16304791cc9cf303a962eaef1b1cb0a8fd193b98acfd414303d",
            `results_summary/Study_report/study_parameters.txt` = "<present>",
            `source/contrasts_tbl.tab` = "7e2814a29d730c668dca584b855fe542882c63a393a5c4deaba6a98f42d5d656",
            `source/design_matrix.tab` = "8812863972c9c16304791cc9cf303a962eaef1b1cb0a8fd193b98acfd414303d",
            `source/study_parameters.txt` = "<present>",
            `time/12_correlation_filtered_combined_plots.png` = "<present>"),
        study_parameters = "1 lines 55c1b3aadde3f3b7",
        integration = stats::setNames(character(0), character()),
        session_state = "<0 files>",
        errors = "<unhandled> unused arguments (self, name)"),
    report = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK] \nReport location: <ROOT>/results_summary/TMT_report_proteomics_kit108.docx",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = c(
            "[message 5 -] TMT_report.rmd template copied from package",
            "[message 5 -] Report generated successfully!"),
        updates = character(0),
        progress = c(
            "start Generating report... 0.1", "0.1 Detecting workflow type...",
            "0.2 Checking for TMT_report.rmd template...", "0.5 Rendering report..."),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/TMT_report_proteomics_kit108.docx",
            workflow_args_saved = "FALSE"),
        calls = "RenderReport(omic_type=proteomics, experiment_label=kit108, rmd_filename=TMT_report.rmd, project_dirs_object_name=project_dirs, output_format=<NULL>, global=the session's)",
        global = c(before = "the session's", after = "the session's"),
        files = c(
            `results_summary/TMT_report_proteomics_kit108.docx` = "<present>",
            `scripts/proteomics/TMT_report.rmd` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = "<unhandled> unused arguments (self, name)"),
    export = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK] \nReport location: <ROOT>/results_summary/TMT_report_proteomics_kit108.docx",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = "[message 5 -] Session state exported to: <ROOT>/source/session_state_<DATE>.RDS",
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/TMT_report_proteomics_kit108.docx",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "the session's", after = "the session's"),
        files = c(`source/session_state_<DATE>.RDS` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = c(
            experiment_label = "kit108",
            description = "Full protein analysis workflow with config parameters", timestamp = "<time>",
            omic_type = "proteomics", workflow_args_saved = "FALSE", files_copied = "TRUE",
            report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/TMT_report_proteomics_kit108.docx",
            workflow_type = "TMT_PD", report_template = "TMT_report.rmd",
            parameter_payload = "<digest> b2f01fba10edf976", project_dirs = "<digest> 63ed7c866939a783"),
        errors = character(0)))

omicSummaryW1Frozen[["prot_owner"]] <- list(
    start = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "Templates: DIA-NN [OK], DIA-NN limpa [OK], TMT [OK], LFQ [OK]"),
        notes = character(0),
        updates = "session_summary-experiment_label: <same> -> kit108",
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "the session's", after = "the session's"),
        files = c(
            `<store entries>.store/.workflow-entities/entities/dia-state-<ID>` = "4",
            `<store entries>.store/.workflow-entities/entities/dia-table-<ID>` = "3",
            `da_output/da_groupB-groupA_long_annot.tsv` = "1ece3b87eeeedd82d2783415db679309b23b1cd9378d8d32e671989a020abb79",
            `da_output/da_groupB-groupA_long_annot.xlsx` = "<present>",
            `feature_qc/composite_QC_figure.pdf` = "<present>",
            `feature_qc/composite_QC_figure.png` = "<present>",
            `feature_qc/ruv_normalised_results_cln_with_replicates.RDS` = "<present>",
            `feature_qc/ruv_normalised_results_cln_with_replicates.tsv` = "46d9ba7fe967c79d88fad400b5998407771a19914a2d6f0997b8820e6a5074b5",
            `pathway_enrichment/B_vs_A_down_enrichment_results.tsv` = "4e168310ff617aba595f1605046a7d4782d44f232692d559db1bf56bf3e71568",
            `pathway_enrichment/B_vs_A_up_enrichment_results.tsv` = "b56e92aa35bf8f4fbd3dd96a2e0de7f4ceb9235d00eef687ce41204bf398ad68",
            `publication_graphs/Heatmap/artifact.txt` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/artifact.txt` = "<present>",
            `publication_graphs/Volcano_Plots/artifact.txt` = "<present>",
            `scripts/proteomics/DIANN_limpa_report.rmd` = "<present>",
            `scripts/proteomics/DIANN_report.rmd` = "<present>",
            `scripts/proteomics/LFQ_report.rmd` = "<present>",
            `scripts/proteomics/TMT_report.rmd` = "<present>",
            `source/contrasts_tbl.tab` = "7e2814a29d730c668dca584b855fe542882c63a393a5c4deaba6a98f42d5d656",
            `source/design_matrix.tab` = "8812863972c9c16304791cc9cf303a962eaef1b1cb0a8fd193b98acfd414303d",
            `source/study_parameters.txt` = "<present>",
            `time/12_correlation_filtered_combined_plots.png` = "<present>"),
        study_parameters = "1 lines 55c1b3aadde3f3b7",
        integration = stats::setNames(character(0), character()),
        session_state = "<0 files>",
        errors = "<unhandled> unused arguments (self, name)"),
    save = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Study parameters created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nFile: <ROOT>/source/study_parameters.txt \nSource: Final S4 object @args + config_list \nIntegration Object: proteomics_kit108_final_s4.RDS \nStatus: Parameters saved [OK]",
            template_status = "Templates: DIA-NN [OK], DIA-NN limpa [OK], TMT [OK], LFQ [OK]"),
        notes = c(
            "[message 5 -] Saved Integration S4 Object", "[message 5 -] Study parameters saved successfully"),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = c(
            "createWorkflowArgsFromConfig(workflow_name=kit108, description=Full protein analysis workflow with config parameters, organism_name=<NULL>, taxon_id=<NULL>, source_dir_path=<ROOT>/source, final_s4_object=<ref>, contrasts_tbl=<table 1x3 950ddbb7d830df8c>, workflow_data=<environment>)",
            "writeWorkflowAnalysisRecord(manager=<WorkflowStateOwner>, path=<ROOT>/source/analysis_record.json)"),
        global = c(before = "the session's", after = "the session's"),
        files = c(
            `integration/proteomics_kit108_final_s4.RDS` = "<present>",
            `source/analysis_record.json` = "<present>"),
        study_parameters = "102 lines db0514235ef41471",
        integration = c(proteomics_kit108_final_s4.RDS = "<RDS> <reference>"),
        session_state = "<unchanged>",
        errors = character(0)),
    copy = list(
        outputs = c("copy_status", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK]",
            template_status = "Templates: DIA-NN [OK], DIA-NN limpa [OK], TMT [OK], LFQ [OK]"),
        notes = "[message 5 -] Publication files copied",
        updates = character(0),
        progress = "start Copying files to publication directory... 0.1",
        reactive_values = c(
            files_copied = "TRUE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = c(
            "copyToResultsSummary(omic_type=proteomics, experiment_label=kit108, contrasts_tbl=<table 1x3 950ddbb7d830df8c>, design_matrix=<table 4x3 bcd78aff6b352869>, force=TRUE, current_rmd=<NULL>, project_dirs_object_name=project_dirs, global=the session's)",
            "  copyToResultsSummary returned: "),
        global = c(before = "the session's", after = "the session's"),
        files = c(
            `results_summary/Publication_figures/Enrichment_Plots/B_vs_A_down_enrichment_results.tsv` = "4e168310ff617aba595f1605046a7d4782d44f232692d559db1bf56bf3e71568",
            `results_summary/Publication_figures/Enrichment_Plots/B_vs_A_up_enrichment_results.tsv` = "b56e92aa35bf8f4fbd3dd96a2e0de7f4ceb9235d00eef687ce41204bf398ad68",
            `results_summary/Publication_figures/Heatmap/artifact.txt` = "<present>",
            `results_summary/Publication_figures/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `results_summary/Publication_figures/NumSigDaMolecules/artifact.txt` = "<present>",
            `results_summary/Publication_figures/Volcano_Plots/artifact.txt` = "<present>",
            `results_summary/Publication_tables/DA_results_proteomics.xlsx` = "<present>",
            `results_summary/Publication_tables/Pathway_enrichment_results_proteomics.xlsx` = "<present>",
            `results_summary/Publication_tables/RUV_normalised_results.tsv` = "46d9ba7fe967c79d88fad400b5998407771a19914a2d6f0997b8820e6a5074b5",
            `results_summary/Publication_tables/ruv_normalised_results.RDS` = "<present>",
            `results_summary/QC_figures/12_correlation_filtered_combined_plots.png` = "<present>",
            `results_summary/QC_figures/proteomics_composite_QC_figure.pdf` = "<present>",
            `results_summary/QC_figures/proteomics_composite_QC_figure.png` = "<present>",
            `results_summary/Study_report/contrasts_tbl.tab` = "7e2814a29d730c668dca584b855fe542882c63a393a5c4deaba6a98f42d5d656",
            `results_summary/Study_report/design_matrix.tab` = "8812863972c9c16304791cc9cf303a962eaef1b1cb0a8fd193b98acfd414303d",
            `results_summary/Study_report/study_parameters.txt` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    report = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK] \nReport location: <ROOT>/results_summary/DIANN_report_proteomics_kit108.docx",
            template_status = "Templates: DIA-NN [OK], DIA-NN limpa [OK], TMT [OK], LFQ [OK]"),
        notes = "[message 5 -] Report generated successfully!",
        updates = character(0),
        progress = c(
            "start Generating report... 0.1", "0.1 Detecting workflow type...",
            "0.2 Checking for DIANN_report.rmd template...", "0.5 Rendering report...",
            "0.9 Exporting the analysis record..."),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/DIANN_report_proteomics_kit108.docx",
            workflow_args_saved = "TRUE"),
        calls = c(
            "RenderReport(omic_type=proteomics, experiment_label=kit108, rmd_filename=DIANN_report.rmd, project_dirs_object_name=project_dirs, output_format=<NULL>, global=the session's)",
            "exportWorkflowAnalysisRecordBundle(manager=<WorkflowStateOwner>, directory=<ROOT>/results_summary)",
            "writeWorkflowAnalysisRecord(manager=<WorkflowStateOwner>, path=<ROOT>/results_summary/analysis_record.json)"),
        global = c(before = "the session's", after = "the session's"),
        files = c(
            `results_summary/DIANN_report_proteomics_kit108.docx` = "<present>",
            `results_summary/analysis_record.docx` = "<present>",
            `results_summary/analysis_record.json` = "<present>",
            `results_summary/analysis_record.xlsx` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = "<unhandled> unused arguments (self, name)"),
    export = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK] \nReport location: <ROOT>/results_summary/DIANN_report_proteomics_kit108.docx",
            template_status = "Templates: DIA-NN [OK], DIA-NN limpa [OK], TMT [OK], LFQ [OK]"),
        notes = "[message 5 -] Session state exported to: <ROOT>/source/session_state_<DATE>.RDS",
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/DIANN_report_proteomics_kit108.docx",
            workflow_args_saved = "TRUE"),
        calls = character(0),
        global = c(before = "the session's", after = "the session's"),
        files = c(`source/session_state_<DATE>.RDS` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = c(
            experiment_label = "kit108",
            description = "Full protein analysis workflow with config parameters", timestamp = "<time>",
            omic_type = "proteomics", workflow_args_saved = "TRUE", files_copied = "TRUE",
            report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/DIANN_report_proteomics_kit108.docx", workflow_type = "DIA",
            report_template = "DIANN_report.rmd", parameter_payload = "<digest> bef1e766f17f28fb",
            project_dirs = "<digest> 63ed7c866939a783"),
        errors = character(0)),
    github_disabled = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK] \nReport location: <ROOT>/results_summary/DIANN_report_proteomics_kit108.docx",
            template_status = "Templates: DIA-NN [OK], DIA-NN limpa [OK], TMT [OK], LFQ [OK]"),
        notes = character(0),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/DIANN_report_proteomics_kit108.docx",
            workflow_args_saved = "TRUE"),
        calls = character(0),
        global = c(before = "the session's", after = "the session's"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    github = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK], GitHub pushed [OK]",
            template_status = "Templates: DIA-NN [OK], DIA-NN limpa [OK], TMT [OK], LFQ [OK]"),
        notes = "[message 5 -] Successfully pushed to GitHub",
        updates = character(0),
        progress = "start Pushing to GitHub... 0.1",
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/DIANN_report_proteomics_kit108.docx",
            workflow_args_saved = "TRUE"),
        calls = "pushProjectToGithubFromDirs(project_dirs=the session's, omic_type=proteomics, experiment_label=kit108, project_id=kit108-project, github_org=kit108-org, github_user_email=kit108@example.org, github_user_name=kit108-user)",
        global = c(before = "the session's", after = "the session's"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)))

omicSummaryW1Frozen[["prot_owner_variants"]] <- list(
    start = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = character(0),
        updates = "session_summary-experiment_label: <same> -> kit108",
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "the session's", after = "the session's"),
        files = c(
            `<store entries>.store/.workflow-entities/entities/dia-state-<ID>` = "4",
            `<store entries>.store/.workflow-entities/entities/dia-table-<ID>` = "3",
            `da_output/da_groupB-groupA_long_annot.tsv` = "1ece3b87eeeedd82d2783415db679309b23b1cd9378d8d32e671989a020abb79",
            `da_output/da_groupB-groupA_long_annot.xlsx` = "<present>",
            `feature_qc/composite_QC_figure.pdf` = "<present>",
            `feature_qc/composite_QC_figure.png` = "<present>",
            `feature_qc/ruv_normalised_results_cln_with_replicates.RDS` = "<present>",
            `feature_qc/ruv_normalised_results_cln_with_replicates.tsv` = "46d9ba7fe967c79d88fad400b5998407771a19914a2d6f0997b8820e6a5074b5",
            `pathway_enrichment/B_vs_A_down_enrichment_results.tsv` = "4e168310ff617aba595f1605046a7d4782d44f232692d559db1bf56bf3e71568",
            `pathway_enrichment/B_vs_A_up_enrichment_results.tsv` = "b56e92aa35bf8f4fbd3dd96a2e0de7f4ceb9235d00eef687ce41204bf398ad68",
            `publication_graphs/Heatmap/artifact.txt` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/artifact.txt` = "<present>",
            `publication_graphs/Volcano_Plots/artifact.txt` = "<present>",
            `source/contrasts_tbl.tab` = "7e2814a29d730c668dca584b855fe542882c63a393a5c4deaba6a98f42d5d656",
            `source/design_matrix.tab` = "8812863972c9c16304791cc9cf303a962eaef1b1cb0a8fd193b98acfd414303d",
            `source/study_parameters.txt` = "<present>",
            `time/12_correlation_filtered_combined_plots.png` = "<present>"),
        study_parameters = "1 lines 55c1b3aadde3f3b7",
        integration = stats::setNames(character(0), character()),
        session_state = "<0 files>",
        errors = "<unhandled> unused arguments (self, name)"),
    save_record_fails = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Study parameters created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nFile: <ROOT>/source/study_parameters.txt \nSource: Final S4 object @args + config_list \nIntegration Object: proteomics_kit108_final_s4.RDS \nStatus: Parameters saved [OK]",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = c(
            "[warning NULL -] The analysis record (analysis_record.json) was not written: kit108 record double fails",
            "[message 5 -] Saved Integration S4 Object", "[message 5 -] Study parameters saved successfully"),
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = c(
            "createWorkflowArgsFromConfig(workflow_name=kit108, description=Full protein analysis workflow with config parameters, organism_name=<NULL>, taxon_id=<NULL>, source_dir_path=<ROOT>/source, final_s4_object=<ref>, contrasts_tbl=<table 1x3 950ddbb7d830df8c>, workflow_data=<environment>)",
            "writeWorkflowAnalysisRecord(manager=<WorkflowStateOwner>, path=<ROOT>/source/analysis_record.json)"),
        global = c(before = "the session's", after = "the session's"),
        files = c(`integration/proteomics_kit108_final_s4.RDS` = "<present>"),
        study_parameters = "102 lines db0514235ef41471",
        integration = c(proteomics_kit108_final_s4.RDS = "<RDS> <reference>"),
        session_state = "<unchanged>",
        errors = character(0)),
    copy_gate_refuses = list(
        outputs = c("copy_status", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Error: session outputs are pending; retry export before dependent actions",
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Study parameters created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nFile: <ROOT>/source/study_parameters.txt \nSource: Final S4 object @args + config_list \nIntegration Object: proteomics_kit108_final_s4.RDS \nStatus: Parameters saved [OK]",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = "[error 10 -] Copy error: session outputs are pending; retry export before dependent actions",
        updates = character(0),
        progress = "start Copying files to publication directory... 0.1",
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = character(0),
        global = c(before = "the session's", after = "the session's"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    copy = list(
        outputs = c("copy_status", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK]",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = "[message 5 -] Publication files copied",
        updates = character(0),
        progress = "start Copying files to publication directory... 0.1",
        reactive_values = c(
            files_copied = "TRUE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = c(
            "copyToResultsSummary(omic_type=proteomics, experiment_label=kit108, contrasts_tbl=<table 1x3 950ddbb7d830df8c>, design_matrix=<table 4x3 bcd78aff6b352869>, force=TRUE, current_rmd=<NULL>, project_dirs_object_name=project_dirs, global=the session's)",
            "  copyToResultsSummary returned: "),
        global = c(before = "the session's", after = "the session's"),
        files = c(
            `results_summary/Publication_figures/Enrichment_Plots/B_vs_A_down_enrichment_results.tsv` = "4e168310ff617aba595f1605046a7d4782d44f232692d559db1bf56bf3e71568",
            `results_summary/Publication_figures/Enrichment_Plots/B_vs_A_up_enrichment_results.tsv` = "b56e92aa35bf8f4fbd3dd96a2e0de7f4ceb9235d00eef687ce41204bf398ad68",
            `results_summary/Publication_figures/Heatmap/artifact.txt` = "<present>",
            `results_summary/Publication_figures/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `results_summary/Publication_figures/NumSigDaMolecules/artifact.txt` = "<present>",
            `results_summary/Publication_figures/Volcano_Plots/artifact.txt` = "<present>",
            `results_summary/Publication_tables/DA_results_proteomics.xlsx` = "<present>",
            `results_summary/Publication_tables/Pathway_enrichment_results_proteomics.xlsx` = "<present>",
            `results_summary/Publication_tables/RUV_normalised_results.tsv` = "46d9ba7fe967c79d88fad400b5998407771a19914a2d6f0997b8820e6a5074b5",
            `results_summary/Publication_tables/ruv_normalised_results.RDS` = "<present>",
            `results_summary/QC_figures/12_correlation_filtered_combined_plots.png` = "<present>",
            `results_summary/QC_figures/proteomics_composite_QC_figure.pdf` = "<present>",
            `results_summary/QC_figures/proteomics_composite_QC_figure.png` = "<present>",
            `results_summary/Study_report/contrasts_tbl.tab` = "7e2814a29d730c668dca584b855fe542882c63a393a5c4deaba6a98f42d5d656",
            `results_summary/Study_report/design_matrix.tab` = "8812863972c9c16304791cc9cf303a962eaef1b1cb0a8fd193b98acfd414303d",
            `results_summary/Study_report/study_parameters.txt` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)),
    report_bundle_fails = list(
        outputs = c(
            "copy_status", "download_report", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            download_report = "<download>", report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK] \nReport location: <ROOT>/results_summary/DIANN_report_proteomics_kit108.docx",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = c(
            "[message 5 -] DIANN_report.rmd template copied from package",
            "[message 5 -] Report generated successfully!",
            "[warning NULL -] The analysis record bundle was not exported: kit108 bundle double fails"),
        updates = character(0),
        progress = c(
            "start Generating report... 0.1", "0.1 Detecting workflow type...",
            "0.2 Checking for DIANN_report.rmd template...", "0.5 Rendering report...",
            "0.9 Exporting the analysis record..."),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "TRUE",
            report_path = "<ROOT>/results_summary/DIANN_report_proteomics_kit108.docx",
            workflow_args_saved = "TRUE"),
        calls = c(
            "RenderReport(omic_type=proteomics, experiment_label=kit108, rmd_filename=DIANN_report.rmd, project_dirs_object_name=project_dirs, output_format=<NULL>, global=the session's)",
            "exportWorkflowAnalysisRecordBundle(manager=<WorkflowStateOwner>, directory=<ROOT>/results_summary)"),
        global = c(before = "the session's", after = "the session's"),
        files = c(
            `results_summary/DIANN_report_proteomics_kit108.docx` = "<present>",
            `scripts/proteomics/DIANN_report.rmd` = "<present>"),
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = "<unhandled> unused arguments (self, name)"))

omicSummaryW1Frozen[["prot_native"]] <- list(
    start = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = character(0),
        updates = "session_summary-experiment_label: <same> -> kit108",
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = character(0),
        global = c(before = "absent", after = "absent"),
        files = c(
            `<store entries>.store/.workflow-entities/entities/dia-state-<ID>` = "4",
            `<store entries>.store/.workflow-entities/entities/dia-table-<ID>` = "3",
            `da_output/da_groupB-groupA_long_annot.tsv` = "1ece3b87eeeedd82d2783415db679309b23b1cd9378d8d32e671989a020abb79",
            `da_output/da_groupB-groupA_long_annot.xlsx` = "<present>",
            `feature_qc/composite_QC_figure.pdf` = "<present>",
            `feature_qc/composite_QC_figure.png` = "<present>",
            `feature_qc/ruv_normalised_results_cln_with_replicates.RDS` = "<present>",
            `feature_qc/ruv_normalised_results_cln_with_replicates.tsv` = "46d9ba7fe967c79d88fad400b5998407771a19914a2d6f0997b8820e6a5074b5",
            `pathway_enrichment/B_vs_A_down_enrichment_results.tsv` = "4e168310ff617aba595f1605046a7d4782d44f232692d559db1bf56bf3e71568",
            `pathway_enrichment/B_vs_A_up_enrichment_results.tsv` = "b56e92aa35bf8f4fbd3dd96a2e0de7f4ceb9235d00eef687ce41204bf398ad68",
            `publication_graphs/Heatmap/artifact.txt` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/artifact.txt` = "<present>",
            `publication_graphs/Volcano_Plots/artifact.txt` = "<present>",
            `source/contrasts_tbl.tab` = "7e2814a29d730c668dca584b855fe542882c63a393a5c4deaba6a98f42d5d656",
            `source/design_matrix.tab` = "8812863972c9c16304791cc9cf303a962eaef1b1cb0a8fd193b98acfd414303d",
            `source/study_parameters.txt` = "<present>",
            `time/12_correlation_filtered_combined_plots.png` = "<present>"),
        study_parameters = "1 lines 55c1b3aadde3f3b7",
        integration = stats::setNames(character(0), character()),
        session_state = "<0 files>",
        errors = "<unhandled> unused arguments (self, name)"),
    save = list(
        outputs = c("report_ready", "session_summary", "template_status"),
        values = c(
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Ready to save workflow parameters and generate report",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = "[warning 5 -] Study parameters saved with warnings",
        updates = character(0),
        progress = character(0),
        reactive_values = c(
            files_copied = "FALSE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "TRUE"),
        calls = c(
            "createWorkflowArgsFromConfig(workflow_name=kit108, description=Full protein analysis workflow with config parameters, organism_name=<NULL>, taxon_id=<NULL>, source_dir_path=<ROOT>/source, final_s4_object=<ref>, contrasts_tbl=<table 1x3 950ddbb7d830df8c>, workflow_data=<environment>)",
            "writeWorkflowAnalysisRecord(manager=<WorkflowStateOwner>, path=<ROOT>/source/analysis_record.json)"),
        global = c(before = "absent", after = "absent"),
        files = c(`source/analysis_record.json` = "<present>"),
        study_parameters = "102 lines db0514235ef41471",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)))

omicSummaryW1Frozen[["prot_github_fallback"]] <- list(
    copy = list(
        outputs = c("copy_status", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK]",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = "[message 5 -] Publication files copied",
        updates = "session_summary-experiment_label: <same> -> kit108",
        progress = "start Copying files to publication directory... 0.1",
        reactive_values = c(
            files_copied = "TRUE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = c(
            "copyToResultsSummary(omic_type=proteomics, experiment_label=kit108, contrasts_tbl=<table 1x3 950ddbb7d830df8c>, design_matrix=<table 4x3 bcd78aff6b352869>, force=TRUE, current_rmd=<NULL>, project_dirs_object_name=project_dirs, global=the session's)",
            "  copyToResultsSummary returned: "),
        global = c(before = "absent", after = "the session's"),
        files = c(
            `da_output/da_groupB-groupA_long_annot.tsv` = "1ece3b87eeeedd82d2783415db679309b23b1cd9378d8d32e671989a020abb79",
            `da_output/da_groupB-groupA_long_annot.xlsx` = "<present>",
            `feature_qc/composite_QC_figure.pdf` = "<present>",
            `feature_qc/composite_QC_figure.png` = "<present>",
            `feature_qc/ruv_normalised_results_cln_with_replicates.RDS` = "<present>",
            `feature_qc/ruv_normalised_results_cln_with_replicates.tsv` = "46d9ba7fe967c79d88fad400b5998407771a19914a2d6f0997b8820e6a5074b5",
            `pathway_enrichment/B_vs_A_down_enrichment_results.tsv` = "4e168310ff617aba595f1605046a7d4782d44f232692d559db1bf56bf3e71568",
            `pathway_enrichment/B_vs_A_up_enrichment_results.tsv` = "b56e92aa35bf8f4fbd3dd96a2e0de7f4ceb9235d00eef687ce41204bf398ad68",
            `publication_graphs/Heatmap/artifact.txt` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `publication_graphs/NumSigDaMolecules/artifact.txt` = "<present>",
            `publication_graphs/Volcano_Plots/artifact.txt` = "<present>",
            `results_summary/Publication_figures/Enrichment_Plots/B_vs_A_down_enrichment_results.tsv` = "4e168310ff617aba595f1605046a7d4782d44f232692d559db1bf56bf3e71568",
            `results_summary/Publication_figures/Enrichment_Plots/B_vs_A_up_enrichment_results.tsv` = "b56e92aa35bf8f4fbd3dd96a2e0de7f4ceb9235d00eef687ce41204bf398ad68",
            `results_summary/Publication_figures/Heatmap/artifact.txt` = "<present>",
            `results_summary/Publication_figures/Interactive_Volcano_Plots/artifact.txt` = "<present>",
            `results_summary/Publication_figures/NumSigDaMolecules/artifact.txt` = "<present>",
            `results_summary/Publication_figures/Volcano_Plots/artifact.txt` = "<present>",
            `results_summary/Publication_tables/DA_results_proteomics.xlsx` = "<present>",
            `results_summary/Publication_tables/Pathway_enrichment_results_proteomics.xlsx` = "<present>",
            `results_summary/Publication_tables/RUV_normalised_results.tsv` = "46d9ba7fe967c79d88fad400b5998407771a19914a2d6f0997b8820e6a5074b5",
            `results_summary/Publication_tables/ruv_normalised_results.RDS` = "<present>",
            `results_summary/QC_figures/12_correlation_filtered_combined_plots.png` = "<present>",
            `results_summary/QC_figures/proteomics_composite_QC_figure.pdf` = "<present>",
            `results_summary/QC_figures/proteomics_composite_QC_figure.png` = "<present>",
            `results_summary/Study_report/contrasts_tbl.tab` = "7e2814a29d730c668dca584b855fe542882c63a393a5c4deaba6a98f42d5d656",
            `results_summary/Study_report/design_matrix.tab` = "8812863972c9c16304791cc9cf303a962eaef1b1cb0a8fd193b98acfd414303d",
            `results_summary/Study_report/study_parameters.txt` = "<present>",
            `source/contrasts_tbl.tab` = "7e2814a29d730c668dca584b855fe542882c63a393a5c4deaba6a98f42d5d656",
            `source/design_matrix.tab` = "8812863972c9c16304791cc9cf303a962eaef1b1cb0a8fd193b98acfd414303d",
            `source/study_parameters.txt` = "<present>",
            `time/12_correlation_filtered_combined_plots.png` = "<present>"),
        study_parameters = "1 lines 55c1b3aadde3f3b7",
        integration = stats::setNames(character(0), character()),
        session_state = "<0 files>",
        errors = "<unhandled> unused arguments (self, name)"),
    report_fallback = list(
        outputs = c("copy_status", "report_ready", "session_summary", "template_status"),
        values = c(
            copy_status = "Files copied to publication directory successfully [OK]",
            report_ready = "<error> unused arguments (self, name)",
            session_summary = "Workflow args created for: kit108 \nDescription: Full protein analysis workflow with config parameters \nTimestamp: <TIME> \nStatus: Arguments saved [OK], Files copied [OK]",
            template_status = "[WARNING] Report templates will be downloaded when generating report"),
        notes = "[error 10 -] Template retrieval failed: kit108 download double: offline",
        updates = character(0),
        progress = c(
            "start Generating report... 0.1", "0.1 Detecting workflow type...",
            "0.2 Checking for DIANN_report.rmd template..."),
        reactive_values = c(
            files_copied = "TRUE", report_generated = "FALSE", report_path = "<NULL>",
            workflow_args_saved = "FALSE"),
        calls = c(
            "system.file(path=reports/proteomics/DIANN_report.rmd, package=MultiScholaR)",
            "download.file(url=https://raw.githubusercontent.com/APAF-bioinformatics/MultiScholaR/main/inst/reports/proteomics/DIANN_report.rmd, destfile=<ROOT>/scripts/proteomics/DIANN_report.rmd, quiet=TRUE)"),
        global = c(before = "the session's", after = "the session's"),
        files = "<unchanged>",
        study_parameters = "<unchanged>",
        integration = "<unchanged>",
        session_state = "<unchanged>",
        errors = character(0)))
# OMIC_SUMMARY_FROZEN_END

# MULTISCHOLAR_KIT108_W1_CASES=<case,...> runs only those cases (the dump's and the probes' scratch runs).
omicSummaryW1Selected <- if (nzchar(Sys.getenv("MULTISCHOLAR_KIT108_W1_CASES"))) {
    intersect(strsplit(Sys.getenv("MULTISCHOLAR_KIT108_W1_CASES"), ",", fixed = TRUE)[[1L]], omicSummaryW1Cases)
} else omicSummaryW1Cases

for (case in omicSummaryW1Selected) {
    local({
        name <- case
        test_that(sprintf("(b) %s: the outputs, notifications, progress, reactive values, calls, global and files of each step are today's", name), {
            actual <- omicSummaryW1Case(name)
            expected <- omicSummaryW1Frozen[[name]]
            expect_named(actual, names(expected))
            for (step in names(expected)) {
                expect_named(actual[[step]], names(expected[[step]]), info = paste(name, step))
                for (field in names(expected[[step]])) {
                    expect_identical(actual[[step]][[field]], expected[[step]][[field]], info = paste(name, step, field))
                }
            }
        })
    })
}

# KIT-108 W3, F108-03 (the lead's decision): trx's save before the RUVIII export refuses the integration write with the
# gate's warning, and shows no "Saved Integration" notice (the start commit's handler showed one, and wrote nothing).
# Written first: it fails on the start commit's handlers (kit-108/w3/test-first.txt).
test_that("(W3, F108-03) trx: the save before the RUVIII export shows the gate's warning and no \"Saved Integration\" notice", {
    spec <- omicSummaryW1CaseSpec("trx_main", environment())
    captures <- omicSummaryW1Drive(spec$setup, spec$steps[c("norm", "da", "report_before_copy", "save_before_export")])
    saved <- captures$save_before_export
    expect_identical(saved$notes, c(paste("[warning 10 -] The integration export needs the RUVIII export: run it in the",
        "panel above, then save again."), "[message 5 -] Study parameters saved successfully"))
    expect_match(saved$values[["session_summary"]], "Integration Object: not written (no integration export yet)", fixed = TRUE)
    expect_identical(saved$reactive_values[["workflow_args_saved"]], "TRUE")
    expect_identical(list.files(spec$setup$paths$integration_dir), character())
})
