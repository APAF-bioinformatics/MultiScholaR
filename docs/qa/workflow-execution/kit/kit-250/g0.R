# KIT-250 (r1.2) gate G0: base observation of the proteomics import's organism dialog ('Mixed species FASTA'), the
# application driven in a headless browser on a GitHub-hosted runner (run-tests scripts job). Observation only: no
# package code is changed; the commit under test is v0.9-staging 16112917 and this branch adds files under
# docs/qa/workflow-execution/kit/kit-250/ only (checked at the start, as KIT-233 does).
#
#   Rscript docs/qa/workflow-execution/kit/kit-250/g0.R g0 <commit under test>
#
# The driver is KIT-233's journey.R, copied byte for byte from origin/kit-233-journey cce7b04c into
# k233-journey.R (sha256 9b62d7ed2627fbdb64d2c8b03f90530c7fc6de5089b957d139252a9852574965) and sourced: launch in
# test mode on the prot_dia fixture, the network closed by a network namespace (KIT-233 5.1), k233Press, k233Js,
# k233StepHome, k233StepDesign. Only the import, the Cancel, the Confirm and the readings are this file's.
#
# Cases (KIT-250 section 10, G0; one fresh application session each):
#   A   (a)  the form's Taxonomy ID set to 10090 and the name to "Mus musculus" before ticking 'Mixed species FASTA';
#            import; the dialog read; Cancel by a real mouse event (chromote Input.dispatchMouseEvent at the button's
#            centre, not a scripted click()); the values read.
#        (a2) (r1.2, R1-01) after that Cancel: is the design builder shown, can a design be saved, and with which
#            taxon (the design step's log line "Target organism taxon ID: <id>").
#   B   (b)  the dialog open; a real mouse event at the centre of the "Design Matrix" tab link.
#   C   (c)  Confirm (the preselected organism, filter unticked).
#   D1  (d1) as (a) up to Cancel, then "Process Imported Data" again with the option still ticked, then Confirm.
#            (In its own session so that no design save lies between the Cancel and the re-import.)
#   D2  (d2) as (a) up to Cancel, then the option unticked and "Process Imported Data" again.
#
# What the browser cannot read on base without package code: workflowData$taxon_id and
# workflowData$mixed_species_analysis themselves. Read instead: the setup log (processing_log$setup_import, written
# with workflowData$taxon_id at the same points: finalizeProtImportSetupState :159/:176-188 and Confirm :172-196),
# the form's inputs (Confirm updates them, :175-176), and in (a2) the design step's own line, which prints the
# workflowData$taxon_id it was given (R/mod_prot_design_state_helpers.R:47; R/func_prot_annotation_uniprot_full_helpers.R:122).

K250_DIR <- local({
    file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    if (length(file)) dirname(normalizePath(file[[1L]])) else normalizePath("docs/qa/workflow-execution/kit/kit-250")
})
source(file.path(K250_DIR, "k233-journey.R"))
K233_RECORDS <- "docs/qa/workflow-execution/kit/kit-250"

K250_OUT <- "k250-out"
K250_BOUNDS <- list(deadline_min = 35, hard_s = 5 * 60, finding_s = 600, poll = 0.5)
K250_NS <- paste0(K233_NS, "setup_import-")
K250_TABS <- paste0(K233_NS, "workflow_tabs")
K250_FORM_TAXON <- 10090
K250_FORM_NAME <- "Mus musculus"
K250_ORGANISM_MODAL <- "Select Primary Organism"
K250_DESIGN_TAB <- sprintf("#%s a[data-value='design']", K250_TABS)
K250_CANCEL <- paste("#shiny-modal .modal-footer button[data-dismiss='modal'],",
    "#shiny-modal .modal-footer button[data-bs-dismiss='modal']")
K250_BUILDER_RUNS <- paste0(K233_NS, "design_matrix-builder-selected_runs")
K250_DATA_AVAILABLE <- paste0(K233_NS, "design_matrix-data_available")
K250_DONE_LINE <- "Data import completed successfully!"

# ---------------------------------------------------------------------------------------------------------------
# Predictions, written before any value is read (KIT-250 2.2 to 2.4, 10 and Re-check 1). "~" = a regular expression
# the observed text must match; "(recorded)" = no prediction. Written to g0-predictions.tsv before the first session.
# ---------------------------------------------------------------------------------------------------------------
k250P <- function(case, moment, field, predicted, basis) {
    list(case = case, moment = moment, field = field, predicted = predicted, basis = basis)
}
K250_OK <- "~^\\[OK\\] Data imported successfully!"
K250_PREDICTIONS <- list(
    # (a) the dialog open, before any answer
    k250P("a", "dialog_open", "modal_open", K250_ORGANISM_MODAL, "2.2 step 1"),
    k250P("a", "dialog_open", "dialog_selected", "9606", "fixture FASTA: one organism, preselected (func_prot_import_organisms.R:223)"),
    k250P("a", "dialog_open", "setup_import", "complete", "2.2 step 3: complete() before the answer"),
    k250P("a", "dialog_open", "design_matrix", "pending", "2.2 step 3"),
    k250P("a", "dialog_open", "setup_log_taxon_id", "10090", "2.2 step 2: the form's value"),
    k250P("a", "dialog_open", "setup_log_organism", K250_FORM_NAME, "2.2 step 2"),
    k250P("a", "dialog_open", "mixed_species_selection", "absent", "2.2 step 2"),
    k250P("a", "dialog_open", "import_status", K250_OK, "2.2 P02"),
    k250P("a", "dialog_open", "form_taxon_disabled", "TRUE", "2.1"),
    # (a) after Cancel by a real mouse event
    k250P("a", "after_cancel", "cancel_hit_button", "TRUE", "the Cancel button is the topmost element at its centre"),
    k250P("a", "after_cancel", "modal_open", "none", "2.2 step 4"),
    k250P("a", "after_cancel", "setup_import", "complete", "2.2 step 4; stop rule 1 if not"),
    k250P("a", "after_cancel", "design_matrix", "pending", "2.2 step 4"),
    k250P("a", "after_cancel", "later_tabs", "all disabled", "2.2 step 3"),
    k250P("a", "after_cancel", "setup_log_taxon_id", "10090", "2.2 step 4: the form's value kept"),
    k250P("a", "after_cancel", "setup_log_organism", K250_FORM_NAME, "2.2 step 4"),
    k250P("a", "after_cancel", "mixed_species_selection", "absent", "2.2 step 4"),
    k250P("a", "after_cancel", "import_status", K250_OK, "2.2 P02"),
    k250P("a", "after_cancel", "data_summary", "~Organism: Mus musculus \\(taxon: 10090 \\)", "2.2 P02 (ui_helpers.R:42-44)"),
    k250P("a", "after_cancel", "notices", "none", "2.2 step 4: Cancel runs no server code"),
    k250P("a", "after_cancel", "notice_naming_organism", "none", "2.2 step 4; stop rule 1 if not"),
    k250P("a", "after_cancel", "form_taxon_id", "10090", "a disabled input keeps its value"),
    k250P("a", "after_cancel", "form_organism_name", K250_FORM_NAME, "a disabled input keeps its value"),
    # (a2) after that Cancel: the design builder and a design save (Re-check 1, R1-01)
    k250P("a2", "design", "design_tab_active", "TRUE", "R1-01: statuses close no tab (mod_omic_workflow.R:32-48)"),
    k250P("a2", "design", "data_available", "TRUE", "R1-01: data_tbl and config_list exist (mod_omic_design.R:127-130)"),
    k250P("a2", "design", "builder_shown", "TRUE", "R1-01 (mod_omic_design.R:77)"),
    k250P("a2", "design", "design_save_outcome", "complete", "R1-01: a design can be saved after Cancel"),
    k250P("a2", "design", "design_matrix", "complete", "R1-01"),
    k250P("a2", "design", "design_log_taxon_id", "10090", "R1-01: saved with the form's organism (uniprot_full_helpers.R:122)"),
    # (b) the dialog open, a real mouse event on the Design Matrix tab link
    k250P("b", "after_click", "click_reached_tab", "FALSE", "2.3: the dialog's backdrop takes the mouse (INFERRED)"),
    k250P("b", "after_click", "design_tab_active", "FALSE", "2.3"),
    k250P("b", "after_click", "modal_open", K250_ORGANISM_MODAL, "2.3: easyClose = FALSE"),
    # (c) Confirm
    k250P("c", "after_confirm", "modal_open", "none", "2.2 step 5"),
    k250P("c", "after_confirm", "setup_import", "complete", "2.2 step 5"),
    k250P("c", "after_confirm", "design_matrix", "pending", "2.2 step 5"),
    k250P("c", "after_confirm", "setup_log_taxon_id", "9606", "2.2 step 5 (:178-196)"),
    k250P("c", "after_confirm", "setup_log_organism", "Homo sapiens", "2.2 step 5"),
    k250P("c", "after_confirm", "mixed_species_selection", "9606", "2.2 step 5"),
    k250P("c", "after_confirm", "mixed_filter_applied", "FALSE", "filter unticked by default (:31-35)"),
    k250P("c", "after_confirm", "notice_primary", "present", "2.2 step 5 (:222-227)"),
    k250P("c", "after_confirm", "notice_success", "present", "2.2 step 5 (:222-227)"),
    k250P("c", "after_confirm", "import_status", K250_OK, "ui_helpers.R:1-19"),
    k250P("c", "after_confirm", "data_summary", "~Organism: Homo sapiens \\(taxon: 9606 \\)", "the setup log updated by Confirm"),
    k250P("c", "after_confirm", "form_taxon_id", "9606", "Confirm updates the form (:175)"),
    k250P("c", "after_confirm", "form_organism_name", "Homo sapiens", "Confirm updates the form (:176)"),
    # (d1) Cancel, then a ticked re-import and Confirm
    k250P("d1", "after_cancel", "setup_import", "complete", "as (a)"),
    k250P("d1", "after_cancel", "setup_log_taxon_id", "10090", "as (a)"),
    k250P("d1", "after_reimport", "reimport_dialog_shown", "TRUE", "the option still ticked"),
    k250P("d1", "after_reimport", "modal_open", "none", "Confirm removes the dialog"),
    k250P("d1", "after_reimport", "setup_import", "complete", "stop rule 4 if refused or failing"),
    k250P("d1", "after_reimport", "design_matrix", "pending", "completeOmicImport (mod_omic_import.R:97-109)"),
    k250P("d1", "after_reimport", "later_tabs", "all disabled", "no later step run"),
    k250P("d1", "after_reimport", "setup_log_taxon_id", "9606", "Confirm"),
    k250P("d1", "after_reimport", "setup_log_organism", "Homo sapiens", "Confirm"),
    k250P("d1", "after_reimport", "mixed_species_selection", "9606", "Confirm"),
    k250P("d1", "after_reimport", "notice_primary", "present", "Confirm"),
    k250P("d1", "after_reimport", "notice_success", "present", "Confirm"),
    k250P("d1", "after_reimport", "notice_error", "none", "stop rule 4 if not"),
    # (d2) Cancel, then the option unticked and a re-import
    k250P("d2", "after_cancel", "setup_import", "complete", "as (a)"),
    k250P("d2", "after_cancel", "setup_log_taxon_id", "10090", "as (a)"),
    k250P("d2", "after_reimport", "import_completed", "TRUE", "stop rule 4 if not"),
    k250P("d2", "after_reimport", "modal_open", "Processing Data", "P01: the waiting flag stays TRUE, so the processing dialog is not removed"),
    k250P("d2", "after_reimport", "setup_import", "complete", "stop rule 4 if not"),
    k250P("d2", "after_reimport", "design_matrix", "pending", "completeOmicImport"),
    k250P("d2", "after_reimport", "later_tabs", "all disabled", "no later step run"),
    k250P("d2", "after_reimport", "setup_log_taxon_id", "10090", "the form's value, now enabled"),
    k250P("d2", "after_reimport", "setup_log_organism", K250_FORM_NAME, "the form's value"),
    k250P("d2", "after_reimport", "mixed_species_selection", "absent", "a fresh setup log without the option"),
    k250P("d2", "after_reimport", "form_mixed_species", "FALSE", "unticked by the driver"),
    k250P("d2", "after_reimport", "notice_success", "absent", "P01: the waiting branch shows nothing (:218-223)"),
    k250P("d2", "after_reimport", "notice_error", "none", "stop rule 4 if not"),
    k250P("d2", "after_reimport", "import_status", K250_OK, "processing FALSE and data_tbl present")
)

# ---------------------------------------------------------------------------------------------------------------
# Outer and inner phases (KIT-233's, with this file as the inner script).
# ---------------------------------------------------------------------------------------------------------------
k250Main <- function(args) {
    if (length(args) < 2L) stop("usage: g0.R g0 <commit under test> | journey-inner <phase file>")
    if (identical(args[[1L]], "journey-inner")) return(invisible(k250Inner(args[[2L]])))
    if (!identical(args[[1L]], "g0")) stop("unknown run mode ", args[[1L]])
    k250Outer(args[[2L]])
}

k250Outer <- function(commit) {
    started <- Sys.time()
    pkg <- normalizePath(".")
    out <- file.path(pkg, K250_OUT)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    work <- file.path(tempdir(), "k250-work")
    dir.create(work, recursive = TRUE, showWarnings = FALSE)
    launch <- k233LaunchSpec()
    env <- k233OuterEnvironment("g0", commit, pkg, launch)
    k233WriteTsv(env$rows, file.path(out, "environment.tsv"))
    k233WriteTsv(K250_PREDICTIONS, file.path(out, "g0-predictions.tsv"))
    if (!isTRUE(env$commit_ok)) {
        writeLines("stopped: the records commit differs from the commit under test outside the records folder",
            file.path(out, "stop.txt"))
        return(k233Finish(out))
    }
    phase <- list(
        run_mode = "g0", started = as.numeric(started), bounds = K250_BOUNDS, pkg = pkg, out = out, work = work,
        script = file.path(K250_DIR, "g0.R"), commit = commit, records_commit = env$head,
        inputs = list(id = "prot_dia fixture", sha = k233Sha(file.path(pkg, "tests/testdata/e2e/prot_dia/report.tsv"))),
        pandoc = list(dir = "", source = "not used (G0 renders no report)"), launch = launch,
        user_outer = Sys.info()[["effective_user"]], libpaths = .libPaths(), locale = Sys.getlocale("LC_CTYPE"),
        env = c(R_LIBS = Sys.getenv("R_LIBS"), R_LIBS_USER = Sys.getenv("R_LIBS_USER"),
            R_LIBS_SITE = Sys.getenv("R_LIBS_SITE"), LANG = Sys.getenv("LANG", "C.UTF-8"),
            GITHUB_RUN_ID = Sys.getenv("GITHUB_RUN_ID"), TZ = "UTC", OMP_NUM_THREADS = "1",
            OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", NOT_CRAN = "true")
    )
    phase$env <- phase$env[nzchar(phase$env)]
    phase_file <- file.path(work, "k250-phase.rds")
    saveRDS(phase, phase_file)
    k233RunInner(phase, phase_file)
    k233Finish(out)
}

k250Inner <- function(phase_file) {
    phase <- readRDS(phase_file)
    .libPaths(phase$libpaths)
    if (length(phase$env)) do.call(Sys.setenv, as.list(phase$env))
    Sys.unsetenv(c("MULTISCHOLAR_TEST_FEATURES", "MULTISCHOLAR_TEST_MODE", "CHROMOTE_CHROME"))
    invisible(tryCatch(Sys.setlocale("LC_ALL", phase$locale), warning = \(w) NULL, error = \(e) NULL))
    options(warn = 1, chromote.timeout = 60)
    writeLines(k233Utc(), file.path(phase$out, "inner-started.txt"))
    G <- k233Context(phase)
    G$k250 <- list()
    on.exit({
        for (S in G$sessions) k233StopSession(S)
        k233WriteAll(G)
        k250WriteValues(G)
    }, add = TRUE)
    fatal <- function(e) {
        G$fatal <- paste(conditionMessage(e), "| call:", paste(deparse(conditionCall(e)), collapse = " "))
        writeLines(G$fatal, file.path(G$out, "inner-fatal.txt"))
    }
    tryCatch(k250Plan(G), error = fatal, skip = fatal)
    invisible(G)
}

# ---------------------------------------------------------------------------------------------------------------
# The page: readings and the real mouse event.
# ---------------------------------------------------------------------------------------------------------------
# Shiny's dialog is the element #shiny-modal itself (class "modal"); shown while its computed display is not none.
K250_MODAL_TITLE_JS <- paste0("(function(){var m=document.getElementById('shiny-modal');",
    "if(!m||window.getComputedStyle(m).display==='none'){return 'none';}var t=m.querySelector('.modal-title');",
    "return t?(t.innerText||t.textContent||'').replace(/\\s+/g,' ').trim():'(untitled)';})()")

K250_POINT_JS <- r"---((function (sel, scroll) {
  var e = document.querySelector(sel);
  if (!e) { return { found: false }; }
  var b0 = e.getBoundingClientRect();
  var before = b0.bottom > 0 && b0.right > 0 && b0.top < window.innerHeight && b0.left < window.innerWidth;
  if (scroll) { window.scrollTo(0, 0); e.scrollIntoView({ block: 'center' }); }
  var b = e.getBoundingClientRect();
  var x = b.left + b.width / 2, y = b.top + b.height / 2;
  var inview = x >= 0 && y >= 0 && x <= window.innerWidth && y <= window.innerHeight;
  var top = inview ? document.elementFromPoint(x, y) : null;
  var d = function (n) {
    if (!n) { return ''; }
    var c = (typeof n.className === 'string' && n.className.trim()) ? '.' + n.className.trim().split(/\s+/).join('.') : '';
    return n.tagName.toLowerCase() + (n.id ? '#' + n.id : '') + c;
  };
  return { found: true, x: x, y: y, in_view_before: before, in_view: inview, top: d(top),
    hits: !!top && (top === e || e.contains(top)), text: (e.innerText || e.textContent || '').replace(/\s+/g, ' ').trim() };
})(%s, %s))---"

k250MouseClick <- function(S, selector, scroll = TRUE) {
    p <- k233Js(S, sprintf(K250_POINT_JS, k233Json(selector), if (scroll) "true" else "false"))
    if (!is.list(p) || !isTRUE(p$found) || !isTRUE(p$in_view)) {
        return(c(if (is.list(p)) p else list(), list(dispatched = FALSE, how = "not found or not in view")))
    }
    cs <- tryCatch(S$driver$get_chromote_session(), error = \(e) NULL)
    if (is.null(cs)) return(c(p, list(dispatched = FALSE, how = "no chromote session")))
    k233LogAction(S, "G0", selector, sprintf("mouse at %.1f,%.1f", p$x, p$y))
    done <- tryCatch({
        cs$Input$dispatchMouseEvent(type = "mouseMoved", x = p$x, y = p$y)
        cs$Input$dispatchMouseEvent(type = "mousePressed", x = p$x, y = p$y, button = "left", clickCount = 1)
        cs$Input$dispatchMouseEvent(type = "mouseReleased", x = p$x, y = p$y, button = "left", clickCount = 1)
        TRUE
    }, error = \(e) conditionMessage(e))
    c(p, list(dispatched = isTRUE(done), how = if (isTRUE(done)) "chromote Input.dispatchMouseEvent" else done))
}

k250ModalOpen <- function(S) as.character(k233Js(S, K250_MODAL_TITLE_JS) %||% "(no answer)")

k250WaitModal <- function(S, title, seconds) {
    k233WaitFor(S, sprintf("(%s) === %s", K250_MODAL_TITLE_JS, k233Json(title)), seconds)
}

k250DesignTabActive <- function(S) {
    r <- k233Js(S, sprintf(K233_TAB_JS, k233Json(K250_TABS), k233Json("design"), "false"))
    isTRUE(r$ok)
}

k250Snap <- function(S, case, moment, since, extra = list()) {
    G <- S$G
    Sys.sleep(1.5)
    ids <- paste0(K250_NS, c("import_status", "data_summary"))
    st <- k233State(S, since, ids) %||% list()
    ts <- st$dg$step_status_per_omic$proteomics
    log <- st$dg$export_paths$proteomics$setup_import
    later <- ts[setdiff(names(ts), c("setup_import", "design_matrix"))]
    notes <- st$notes %||% list()
    ntext <- vapply(notes, \(n) sprintf("[%s] %s", n$type, n$text), character(1))
    has <- function(pattern) any(vapply(notes, \(n) grepl(pattern, n$text), logical(1)))
    errors <- ntext[vapply(notes, \(n) identical(n$type, "error"), logical(1))]
    naming <- ntext[grepl("Mus musculus|10090|Homo sapiens|9606|organism", ntext, ignore.case = TRUE)]
    shot <- file.path(G$out, sprintf("g0-%s-%s.png", case, moment))
    tryCatch(S$driver$get_screenshot(shot, selector = "viewport"), error = \(e) NULL)
    jsonlite::write_json(list(tab_status = ts, setup_import = log, notes = notes, modals = st$modals, texts = st$texts,
        stepper = st$steps), file.path(G$out, sprintf("g0-%s-%s.json", case, moment)), auto_unbox = TRUE,
        pretty = TRUE, null = "null", force = TRUE)
    input <- function(id) k233InputValue(S, paste0(K250_NS, id))
    selection <- log$mixed_species_selection
    fields <- list(
        modal_open = k250ModalOpen(S),
        setup_import = ts$setup_import %||% "(absent)",
        design_matrix = ts$design_matrix %||% "(absent)",
        later_tabs = if (length(later) && all(unlist(later) == "disabled")) "all disabled" else
            paste(names(later), unlist(later), sep = "=", collapse = " "),
        tab_status_all = paste(names(ts), unlist(ts), sep = "=", collapse = " "),
        setup_log_taxon_id = log$taxon_id %||% "(absent)",
        setup_log_organism = log$organism %||% "(absent)",
        mixed_species_selection = if (is.null(selection)) "absent" else selection$selected_taxon_id %||% "(present)",
        mixed_filter_applied = if (is.null(selection)) "absent" else selection$filter_applied %||% "(present)",
        import_status = k233Text(st, ids[[1L]]),
        data_summary = k233Text(st, ids[[2L]]),
        form_taxon_id = input("taxon_id"),
        form_organism_name = input("organism_name"),
        form_mixed_species = input("mixed_species_fasta"),
        form_taxon_disabled = k233Js(S, sprintf("(function(){var e=document.getElementById(%s);return e?e.disabled:null;})()",
            k233Json(paste0(K250_NS, "taxon_id")))),
        notices = if (length(ntext)) paste(ntext, collapse = " || ") else "none",
        notice_success = if (has("^Data import successful!$")) "present" else "absent",
        notice_primary = if (has("Primary organism set to")) "present" else "absent",
        notice_error = if (length(errors)) paste(errors, collapse = " || ") else "none",
        notice_naming_organism = if (length(naming)) paste(naming, collapse = " || ") else "none",
        setup_log_timestamp = log$timestamp %||% "(absent)",
        screenshot = basename(shot)
    )
    fields <- c(fields, extra)
    for (k in names(fields)) {
        G$k250[[length(G$k250) + 1L]] <- list(case = case, moment = moment, field = k,
            observed = k233Clean(fields[[k]]), session = S$label, utc = k233Utc())
    }
    k250WriteValues(G)
    invisible(fields)
}

k250Note <- function(S, case, moment, field, value) {
    S$G$k250[[length(S$G$k250) + 1L]] <- list(case = case, moment = moment, field = field,
        observed = k233Clean(value), session = S$label, utc = k233Utc())
}

k250Judge <- function(predicted, observed) {
    if (is.na(observed)) return("not observed")
    if (startsWith(predicted, "~")) return(if (grepl(substring(predicted, 2L), observed)) "as predicted" else "DIFFERS")
    if (identical(trimws(observed), predicted)) "as predicted" else "DIFFERS"
}

k250WriteValues <- function(G) {
    obs <- if (length(G$k250)) k233Bind(G$k250) else
        data.frame(case = character(), moment = character(), field = character(), observed = character(),
            session = character(), utc = character())
    pred <- k233Bind(K250_PREDICTIONS)
    key <- function(d) paste(d$case, d$moment, d$field)
    obs$predicted <- pred$predicted[match(key(obs), key(pred))]
    obs$predicted[is.na(obs$predicted)] <- "(recorded)"
    obs$verdict <- vapply(seq_len(nrow(obs)), \(i) if (identical(obs$predicted[[i]], "(recorded)")) "recorded" else
        k250Judge(obs$predicted[[i]], obs$observed[[i]]), character(1))
    missing <- pred[!key(pred) %in% key(obs), , drop = FALSE]
    if (nrow(missing)) {
        obs <- rbind(obs, data.frame(case = missing$case, moment = missing$moment, field = missing$field,
            observed = "", session = "", utc = "", predicted = missing$predicted, verdict = "not observed"))
    }
    k233WriteTsv(obs[, c("case", "moment", "field", "predicted", "observed", "verdict", "session", "utc")],
        file.path(G$out, "g0-base-values.tsv"))
    invisible(obs)
}

# ---------------------------------------------------------------------------------------------------------------
# The application steps.
# ---------------------------------------------------------------------------------------------------------------
k250Start <- function(G, label) {
    S <- k233Launch(G, label)
    if (!isTRUE(S$answered) || is.null(S$driver)) {
        k250Note(S, label, "launch", "session", paste("not attached:", S$attach_error))
        return(NULL)
    }
    R <- k233FixtureRun(G, S$project_dir)
    R$label <- paste0("k250_", tolower(label))
    home <- k233StepHome(S, R, list(id = "S1", name = "home and project"), G$route)
    if (!identical(home$outcome, "complete")) {
        k250Note(S, label, "launch", "S1", paste(home$outcome, home$message %||% ""))
        return(NULL)
    }
    list(S = S, R = R)
}

k250LogCount <- function(S) sum(grepl(K250_DONE_LINE, k233ReadLog(S$log), fixed = TRUE))

#' The import as a user does it: the two files (first import only), the form's organism set to 10090 / Mus musculus
#' BEFORE the option is ticked (first import only), the option set, then "Process Imported Data". Waits for the
#' organism dialog (ticked) or for the import's completion line in the application's log (unticked).
k250Import <- function(S, R, case, first = TRUE, tick = TRUE) {
    k233ActivateTab(S, K250_TABS, "setup")
    if (first) {
        up <- k233Upload(S, case, paste0(K250_NS, "search_results_standard"), R$report)
        if (!up$ok) return(list(ok = FALSE, how = paste("search results upload:", up$text)))
        k233WaitFor(S, sprintf("(document.getElementById(%s)||{innerText:''}).innerText.trim().length > 0",
            k233Json(paste0(K250_NS, "format_detection"))), 60)
        up <- k233Upload(S, case, paste0(K250_NS, "fasta_file_standard"), R$fasta)
        if (!up$ok) return(list(ok = FALSE, how = paste("FASTA upload:", up$text)))
        k233Set(S, case, paste0(K250_NS, "taxon_id"), K250_FORM_TAXON)
        k233Set(S, case, paste0(K250_NS, "organism_name"), K250_FORM_NAME)
    }
    k233Set(S, case, paste0(K250_NS, "mixed_species_fasta"), tick)
    k233Idle(S, 10)
    done0 <- k250LogCount(S)
    since <- as.numeric(k233Js(S, "Date.now()") %||% 0)
    p <- k233Press(S, case, "prot-import-process", "testid")
    if (!isTRUE(p$ok)) return(list(ok = FALSE, how = paste("Process Imported Data:", p$reason), since = since))
    if (tick) {
        shown <- k250WaitModal(S, K250_ORGANISM_MODAL, 300)
        Sys.sleep(1.5)
        return(list(ok = shown, how = if (shown) "organism dialog shown" else
            paste("no organism dialog within 300 s; open:", k250ModalOpen(S)), since = since))
    }
    started <- Sys.time()
    while (k250LogCount(S) <= done0 && k233Seconds(started) < 300 && S$proc$is_alive()) Sys.sleep(0.5)
    completed <- k250LogCount(S) > done0
    list(ok = completed, completed = completed, since = since,
        how = if (completed) "the import's completion line logged" else "no completion line within 300 s")
}

k250Cancel <- function(S, case) {
    since <- as.numeric(k233Js(S, "Date.now()") %||% 0)
    click <- k250MouseClick(S, K250_CANCEL, scroll = TRUE)
    closed <- k233WaitFor(S, sprintf("(%s) === 'none'", K250_MODAL_TITLE_JS), 15)
    list(since = since, closed = closed, extra = list(
        cancel_hit_button = isTRUE(click$hits), cancel_button_text = click$text %||% "",
        cancel_top_element = click$top %||% "", cancel_dispatch = click$how %||% "",
        cancel_point = sprintf("%s,%s", format(click$x %||% NA), format(click$y %||% NA)), dialog_closed = closed))
}

k250Confirm <- function(S, case) {
    p <- k233Press(S, case, paste0(K250_NS, "confirm_organism"))
    closed <- k233WaitFor(S, sprintf("(%s) === 'none'", K250_MODAL_TITLE_JS), 30)
    list(pressed = isTRUE(p$ok), closed = closed)
}

k250DialogRead <- function(S, case, since) {
    text <- k233Js(S, "(function(){var m=document.querySelector('#shiny-modal');return m?(m.innerText||'').replace(/\\s+/g,' ').trim():'';})()")
    list(dialog_text = text %||% "", dialog_selected = k233InputValue(S, paste0(K250_NS, "selected_organism")) %||% "",
        dialog_filter = k233InputValue(S, paste0(K250_NS, "filter_to_organism")) %||% "")
}

# ---------------------------------------------------------------------------------------------------------------
# The cases.
# ---------------------------------------------------------------------------------------------------------------
k250CaseA <- function(G) {
    s <- k250Start(G, "A")
    if (is.null(s)) return(invisible(NULL))
    S <- s$S
    imp <- k250Import(S, s$R, "a")
    k250Note(S, "a", "import", "how", imp$how)
    if (!isTRUE(imp$ok)) return(k233StopSession(S))
    k250Snap(S, "a", "dialog_open", imp$since, k250DialogRead(S, "a"))
    cancel <- k250Cancel(S, "a")
    k250Snap(S, "a", "after_cancel", cancel$since, cancel$extra)
    # (a2) after that Cancel: the Design tab, its builder and a design save.
    mark <- length(k233ReadLog(S$log))
    since <- as.numeric(k233Js(S, "Date.now()") %||% 0)
    tab <- k233ActivateTab(S, K250_TABS, "design")
    builder <- k233WaitFor(S, sprintf("document.getElementById(%s) !== null", k233Json(K250_BUILDER_RUNS)), 30)
    available <- k233Js(S, sprintf("(window.Shiny&&Shiny.shinyapp&&Shiny.shinyapp.$values)?Shiny.shinyapp.$values[%s]:null",
        k233Json(K250_DATA_AVAILABLE)))
    active <- k250DesignTabActive(S)
    design <- if (builder) {
        k233StepDesign(S, s$R, list(id = "S3", name = "design after Cancel (a2)"), G$route)
    } else {
        list(outcome = "not run: the builder was not shown")
    }
    lines <- k233ReadLog(S$log)
    lines <- lines[seq_along(lines) > mark]
    taxon <- regmatches(lines, regexpr("Target organism taxon ID: [0-9]+", lines))
    k250Snap(S, "a2", "design", since, list(design_tab = tab, design_tab_active = active,
        data_available = if (is.null(available)) "(no value)" else available, builder_shown = builder,
        design_save_outcome = design$outcome, design_save_message = design$message %||% "",
        design_log_taxon_id = if (length(taxon)) paste(unique(sub(".*: ", "", taxon)), collapse = ",") else "(no line)",
        design_log_lines = paste(utils::head(grep("taxon|organism|UniProt|uniprot", lines, value = TRUE), 12L),
            collapse = " || ")))
    k233StopSession(S)
}

k250CaseB <- function(G) {
    s <- k250Start(G, "B")
    if (is.null(s)) return(invisible(NULL))
    S <- s$S
    imp <- k250Import(S, s$R, "b")
    k250Note(S, "b", "import", "how", imp$how)
    if (!isTRUE(imp$ok)) return(k233StopSession(S))
    since <- as.numeric(k233Js(S, "Date.now()") %||% 0)
    click <- k250MouseClick(S, K250_DESIGN_TAB, scroll = TRUE)
    Sys.sleep(2)
    active <- k250DesignTabActive(S)
    k250Snap(S, "b", "after_click", since, list(
        tab_in_view_before_scroll = click$in_view_before %||% NA, tab_in_view = click$in_view %||% NA,
        top_element = click$top %||% "", click_dispatch = click$how %||% "",
        click_reached_tab = isTRUE(click$hits) || active, design_tab_active = active,
        note = "the page was scrolled to the tab link by script before the mouse event (a user cannot scroll behind the dialog)"))
    k233StopSession(S)
}

k250CaseC <- function(G) {
    s <- k250Start(G, "C")
    if (is.null(s)) return(invisible(NULL))
    S <- s$S
    imp <- k250Import(S, s$R, "c")
    k250Note(S, "c", "import", "how", imp$how)
    if (!isTRUE(imp$ok)) return(k233StopSession(S))
    since <- as.numeric(k233Js(S, "Date.now()") %||% 0)
    confirm <- k250Confirm(S, "c")
    k250Snap(S, "c", "after_confirm", since, list(confirm_pressed = confirm$pressed, dialog_closed = confirm$closed))
    k233StopSession(S)
}

k250CaseD <- function(G, which) {
    label <- toupper(which)
    s <- k250Start(G, label)
    if (is.null(s)) return(invisible(NULL))
    S <- s$S
    imp <- k250Import(S, s$R, which)
    k250Note(S, which, "import", "how", imp$how)
    if (!isTRUE(imp$ok)) return(k233StopSession(S))
    cancel <- k250Cancel(S, which)
    k250Snap(S, which, "after_cancel", cancel$since, cancel$extra)
    if (identical(which, "d1")) {
        again <- k250Import(S, s$R, which, first = FALSE, tick = TRUE)
        confirm <- if (isTRUE(again$ok)) k250Confirm(S, which) else list(pressed = FALSE, closed = FALSE)
        k250Snap(S, which, "after_reimport", again$since %||% 0, list(reimport_dialog_shown = isTRUE(again$ok),
            reimport_how = again$how, confirm_pressed = confirm$pressed, dialog_closed = confirm$closed))
    } else {
        again <- k250Import(S, s$R, which, first = FALSE, tick = FALSE)
        Sys.sleep(3)
        k250Snap(S, which, "after_reimport", again$since %||% 0, list(import_completed = isTRUE(again$completed),
            reimport_how = again$how))
    }
    k233StopSession(S)
}

k250Value <- function(obs, case, moment, field) {
    v <- obs$observed[obs$case == case & obs$moment == moment & obs$field == field]
    if (length(v)) v[[1L]] else NA_character_
}

k250Verdict <- function(G) {
    obs <- k250WriteValues(G)
    v <- function(...) k250Value(obs, ...)
    sessions <- vapply(G$session_rows, \(r) isTRUE(r$attached), logical(1))
    lines <- c(
        sprintf("sessions attached: %d of %d (%s)", sum(sessions), length(sessions), paste(names(sessions), collapse = ",")),
        sprintf("predictions: %d as predicted, %d differ, %d not observed, %d recorded without prediction",
            sum(obs$verdict == "as predicted"), sum(obs$verdict == "DIFFERS"), sum(obs$verdict == "not observed"),
            sum(obs$verdict == "recorded")),
        sprintf("stop rule 1 (import NOT complete after Cancel, or a notice naming the organism): %s",
            if (is.na(v("a", "after_cancel", "setup_import"))) "NOT JUDGED (a not observed)" else
                if (!identical(v("a", "after_cancel", "setup_import"), "complete") ||
                    !identical(v("a", "after_cancel", "notice_naming_organism"), "none")) "FIRED" else "not fired"),
        sprintf("stop rule 2 (a person can press behind the dialog): %s",
            if (is.na(v("b", "after_click", "click_reached_tab"))) "NOT JUDGED (b not observed)" else
                if (identical(v("b", "after_click", "click_reached_tab"), "TRUE")) "FIRED (record the finding)" else "not fired"),
        sprintf("stop rule 4 (a second import refused or failing): %s", {
            d1 <- c(v("d1", "after_reimport", "reimport_dialog_shown"), v("d1", "after_reimport", "notice_error"))
            d2 <- c(v("d2", "after_reimport", "import_completed"), v("d2", "after_reimport", "notice_error"))
            if (anyNA(c(d1, d2))) "NOT JUDGED (d1 or d2 not observed)" else
                if (!identical(d1[[1L]], "TRUE") || !identical(d1[[2L]], "none") || !identical(d2[[1L]], "TRUE") ||
                    !identical(d2[[2L]], "none")) "FIRED" else "not fired"
        })
    )
    writeLines(lines, file.path(G$out, "verdict.txt"))
    invisible(lines)
}

k250Plan <- function(G) {
    for (case in list(k250CaseA, k250CaseB, k250CaseC, \(G) k250CaseD(G, "d1"), \(G) k250CaseD(G, "d2"))) {
        if (k233Remaining(G) <= 120) break
        tryCatch(case(G), error = \(e) k250Note(list(G = G, label = "?"), "?", "driver", "error", conditionMessage(e)))
        k250WriteValues(G)
    }
    k250Verdict(G)
}

if (sys.nframe() == 0L) k250Main(commandArgs(trailingOnly = TRUE))
