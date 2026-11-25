# KIT-250 H1 (AUDIT-KIT-250-impl I01 (b)): the fix observed in a headless browser on a runner, on the corrected head.
# Observation only. It reuses G0's driver (`g0.R`, which sources KIT-233's `k233-journey.R`) and replaces the cases and
# the predictions; the predictions below are written before any value is read.
#
#   Rscript docs/qa/workflow-execution/kit/kit-250/h1.R h1 <commit under test>
#
# Session A: the form set to 10090 / Mus musculus, the option ticked, import; the window open (W2, "Organism: not
# confirmed"); Esc and a click on the backdrop; Cancel by a real mouse event; the page and the Design tab after Cancel;
# 'Process Imported Data' again (W1 removed); a second Cancel (one W1, not two); 'Process Imported Data' again and Confirm.
# Session D2: import, Cancel, untick, import: the processing window closes, W1 is gone.

local({
    file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    dir <- if (length(file)) dirname(normalizePath(file[[1L]])) else normalizePath("docs/qa/workflow-execution/kit/kit-250")
    source(file.path(dir, "g0.R"), local = globalenv())
})

K250_OUT <- "k250h1-out"
K250_CANCEL <- paste0("#", K250_NS, "cancel_organism")
K250_W1_ID <- "shiny-notification-prot_import_organism_cancelled"
K250_SAVE_DESIGN <- paste0("#", K233_NS, "design_matrix-builder-save_results")
K250_IMPORT_DESIGN <- paste0("#", K233_NS, "design_matrix-show_import_modal")
K250_W1 <- "~^Import not finished: no organism was chosen\\. Press 'Process Imported Data' again"
K250_W1_NOTE <- "~\\[warning\\] Import not finished: no organism was chosen\\."

K250_PREDICTIONS <- list(
    # Session A: the window open
    k250P("h1a", "dialog_open", "modal_open", K250_ORGANISM_MODAL, "the window is shown"),
    k250P("h1a", "dialog_open", "setup_import", "incomplete", "fix: not complete while the window is open"),
    k250P("h1a", "dialog_open", "design_matrix", "disabled", "fix: no status written"),
    k250P("h1a", "dialog_open", "import_status", "Choose the organism in the window to finish the import.", "W2"),
    k250P("h1a", "dialog_open", "data_summary", "~Organism: not confirmed", "W3"),
    k250P("h1a", "dialog_open", "summary_names_form_organism", "FALSE", "W3: no name, no taxon"),
    k250P("h1a", "dialog_open", "setup_log_taxon_id", "10090", "the setup log is written before the window"),
    # Esc and a backdrop click with the window open (I03)
    k250P("h1a", "after_esc", "modal_open", K250_ORGANISM_MODAL, "I03: keyboard off"),
    k250P("h1a", "after_backdrop", "modal_open", K250_ORGANISM_MODAL, "I03: static backdrop"),
    k250P("h1a", "after_backdrop", "setup_import", "incomplete", "nothing ran"),
    # Cancel by a real mouse event
    k250P("h1a", "after_cancel", "cancel_hit_button", "TRUE", "the button is the topmost element at its centre"),
    k250P("h1a", "after_cancel", "modal_open", "none", "Cancel removes the window"),
    k250P("h1a", "after_cancel", "setup_import", "incomplete", "Cancel clears the import"),
    k250P("h1a", "after_cancel", "design_matrix", "disabled", "no status written"),
    k250P("h1a", "after_cancel", "later_tabs", "all disabled", "nothing ran"),
    k250P("h1a", "after_cancel", "setup_log_taxon_id", "(absent)", "the setup log cleared"),
    k250P("h1a", "after_cancel", "import_status", K250_W1, "W1 in the status panel"),
    k250P("h1a", "after_cancel", "import_status_success", "FALSE", "no success text"),
    k250P("h1a", "after_cancel", "data_summary", "", "no Data Summary (req(dataTbl))"),
    k250P("h1a", "after_cancel", "w1_count", "1", "W1 once, under its fixed id"),
    k250P("h1a", "after_cancel", "notices", K250_W1_NOTE, "W1 as a warning"),
    k250P("h1a", "after_cancel", "notice_success", "absent", "not finished"),
    # The Design tab after Cancel
    k250P("h1a", "design_after_cancel", "design_tab_active", "TRUE", "statuses close no tab"),
    k250P("h1a", "design_after_cancel", "data_available", "FALSE", "data_tbl and config_list cleared"),
    k250P("h1a", "design_after_cancel", "save_design_visible", "FALSE", "no design can be saved"),
    k250P("h1a", "design_after_cancel", "import_existing_design_visible", "TRUE", "outside data_available; own organism"),
    k250P("h1a", "design_after_cancel", "design_matrix", "disabled", "no design saved"),
    # 'Process Imported Data' again: W1 removed at the start
    k250P("h1a", "reimport_dialog_open", "modal_open", K250_ORGANISM_MODAL, "the option still ticked"),
    k250P("h1a", "reimport_dialog_open", "w1_count", "0", "I02: removed when the import starts"),
    k250P("h1a", "reimport_dialog_open", "setup_import", "incomplete", "not complete while the window is open"),
    k250P("h1a", "reimport_dialog_open", "notice_error", "none", "the re-import with the owner attached works"),
    # A second Cancel in the same session
    k250P("h1a", "after_cancel2", "modal_open", "none", "the re-bound button fires"),
    k250P("h1a", "after_cancel2", "w1_count", "1", "I02: one copy, not two"),
    k250P("h1a", "after_cancel2", "setup_import", "incomplete", "cleared again"),
    # 'Process Imported Data' again and Confirm
    k250P("h1a", "after_confirm", "modal_open", "none", "Confirm removes the window"),
    k250P("h1a", "after_confirm", "setup_import", "complete", "Confirm completes the import"),
    k250P("h1a", "after_confirm", "design_matrix", "pending", "completeOmicImport"),
    k250P("h1a", "after_confirm", "setup_log_taxon_id", "9606", "the confirmed organism"),
    k250P("h1a", "after_confirm", "setup_log_organism", "Homo sapiens", "the confirmed organism"),
    k250P("h1a", "after_confirm", "mixed_species_selection", "9606", "Confirm"),
    k250P("h1a", "after_confirm", "notice_primary", "present", "Confirm"),
    k250P("h1a", "after_confirm", "notice_success", "present", "Confirm"),
    k250P("h1a", "after_confirm", "notice_error", "none", "no error"),
    k250P("h1a", "after_confirm", "w1_count", "0", "I02"),
    k250P("h1a", "after_confirm", "import_status", K250_OK, "the import finished"),
    k250P("h1a", "after_confirm", "data_summary", "~Organism: Homo sapiens \\(taxon: 9606 \\)", "the confirmed organism"),
    k250P("h1a", "after_confirm", "design_data_available", "TRUE", "the builder can be used after Confirm"),
    # Session D2: Cancel, untick, import
    k250P("h1d2", "after_cancel", "setup_import", "incomplete", "Cancel clears the import"),
    k250P("h1d2", "after_cancel", "w1_count", "1", "W1"),
    k250P("h1d2", "after_reimport", "import_completed", "TRUE", "the unticked import completes"),
    k250P("h1d2", "after_reimport", "modal_open", "none", "P01 fixed: the processing window closes"),
    k250P("h1d2", "after_reimport", "setup_import", "complete", "completes"),
    k250P("h1d2", "after_reimport", "design_matrix", "pending", "completeOmicImport"),
    k250P("h1d2", "after_reimport", "setup_log_taxon_id", "10090", "the form's organism, now entered by the user"),
    k250P("h1d2", "after_reimport", "setup_log_organism", K250_FORM_NAME, "the form's organism"),
    k250P("h1d2", "after_reimport", "mixed_species_selection", "absent", "no window"),
    k250P("h1d2", "after_reimport", "form_mixed_species", "FALSE", "unticked"),
    k250P("h1d2", "after_reimport", "notice_success", "present", "P01 fixed"),
    k250P("h1d2", "after_reimport", "notice_error", "none", "no error"),
    k250P("h1d2", "after_reimport", "w1_count", "0", "I02"),
    k250P("h1d2", "after_reimport", "import_status", K250_OK, "the import finished")
)

k250Main <- function(args) {
    if (length(args) < 2L) stop("usage: h1.R h1 <commit under test> | journey-inner <phase file>")
    if (identical(args[[1L]], "journey-inner")) return(invisible(k250Inner(args[[2L]])))
    if (!identical(args[[1L]], "h1")) stop("unknown run mode ", args[[1L]])
    k250Outer(args[[2L]])
}

# g0.R's outer phase, with this file as the inner script.
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
        script = file.path(K250_DIR, "h1.R"), commit = commit, records_commit = env$head,
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


# ---------------------------------------------------------------------------------------------------------------
# H1's readings and cases.
# ---------------------------------------------------------------------------------------------------------------
k250W1Count <- function(S) {
    as.character(k233Js(S, sprintf("document.querySelectorAll('#%s').length", K250_W1_ID)) %||% "(no answer)")
}

k250Visible <- function(S, selector) {
    r <- k233Js(S, sprintf(K233_PRESS_JS, k233Json(selector), "false"))
    if (!isTRUE(r$found)) return("not in page")
    isTRUE(r$visible)
}

k250DataAvailable <- function(S) {
    v <- k233Js(S, sprintf("(window.Shiny&&Shiny.shinyapp&&Shiny.shinyapp.$values)?Shiny.shinyapp.$values[%s]:null",
        k233Json(K250_DATA_AVAILABLE)))
    if (is.null(v)) "(no value)" else v
}

k250H1Snap <- function(S, case, moment, since, extra = list()) {
    fields <- k250Snap(S, case, moment, since, c(list(w1_count = k250W1Count(S)), extra))
    k250Note(S, case, moment, "import_status_success",
        grepl("Data imported successfully", fields$import_status, fixed = TRUE))
    k250Note(S, case, moment, "summary_names_form_organism",
        grepl("Mus musculus|10090", fields$data_summary))
    k250WriteValues(S$G)
    invisible(fields)
}

k250Key <- function(S, key, code, keycode) {
    cs <- tryCatch(S$driver$get_chromote_session(), error = \(e) NULL)
    if (is.null(cs)) return("no chromote session")
    tryCatch({
        cs$Input$dispatchKeyEvent(type = "keyDown", key = key, code = code, windowsVirtualKeyCode = keycode,
            nativeVirtualKeyCode = keycode)
        cs$Input$dispatchKeyEvent(type = "keyUp", key = key, code = code, windowsVirtualKeyCode = keycode,
            nativeVirtualKeyCode = keycode)
        "chromote Input.dispatchKeyEvent"
    }, error = \(e) conditionMessage(e))
}

k250MouseAt <- function(S, x, y) {
    cs <- tryCatch(S$driver$get_chromote_session(), error = \(e) NULL)
    if (is.null(cs)) return("no chromote session")
    top <- k233Js(S, sprintf("(function(){var n=document.elementFromPoint(%d,%d);return n?(n.tagName.toLowerCase()+(n.id?'#'+n.id:'')):'';})()",
        x, y))
    tryCatch({
        cs$Input$dispatchMouseEvent(type = "mouseMoved", x = x, y = y)
        cs$Input$dispatchMouseEvent(type = "mousePressed", x = x, y = y, button = "left", clickCount = 1)
        cs$Input$dispatchMouseEvent(type = "mouseReleased", x = x, y = y, button = "left", clickCount = 1)
        paste("mouse at", x, y, "on", top %||% "?")
    }, error = \(e) conditionMessage(e))
}

k250H1CaseA <- function(G) {
    s <- k250Start(G, "A")
    if (is.null(s)) return(invisible(NULL))
    S <- s$S
    imp <- k250Import(S, s$R, "h1a")
    k250Note(S, "h1a", "import", "how", imp$how)
    if (!isTRUE(imp$ok)) return(k233StopSession(S))
    k250H1Snap(S, "h1a", "dialog_open", imp$since, k250DialogRead(S, "h1a"))
    since <- as.numeric(k233Js(S, "Date.now()") %||% 0)
    key <- k250Key(S, "Escape", "Escape", 27L)
    Sys.sleep(1.5)
    k250H1Snap(S, "h1a", "after_esc", since, list(key_dispatch = key))
    since <- as.numeric(k233Js(S, "Date.now()") %||% 0)
    click <- k250MouseAt(S, 10L, 10L)
    Sys.sleep(1.5)
    k250H1Snap(S, "h1a", "after_backdrop", since, list(backdrop_click = click))
    cancel <- k250Cancel(S, "h1a")
    k250H1Snap(S, "h1a", "after_cancel", cancel$since, cancel$extra)
    # The Design tab after Cancel.
    since <- as.numeric(k233Js(S, "Date.now()") %||% 0)
    tab <- k233ActivateTab(S, K250_TABS, "design")
    Sys.sleep(3)
    k250H1Snap(S, "h1a", "design_after_cancel", since, list(design_tab = tab,
        design_tab_active = k250DesignTabActive(S), data_available = k250DataAvailable(S),
        save_design_visible = k250Visible(S, K250_SAVE_DESIGN),
        import_existing_design_visible = k250Visible(S, K250_IMPORT_DESIGN)))
    # 'Process Imported Data' again (the option still ticked), then a second Cancel.
    again <- k250Import(S, s$R, "h1a", first = FALSE, tick = TRUE)
    k250H1Snap(S, "h1a", "reimport_dialog_open", again$since %||% 0, list(reimport_how = again$how))
    if (!isTRUE(again$ok)) return(k233StopSession(S))
    cancel <- k250Cancel(S, "h1a")
    k250H1Snap(S, "h1a", "after_cancel2", cancel$since, cancel$extra)
    # 'Process Imported Data' again and Confirm.
    again <- k250Import(S, s$R, "h1a", first = FALSE, tick = TRUE)
    confirm <- if (isTRUE(again$ok)) k250Confirm(S, "h1a") else list(pressed = FALSE, closed = FALSE)
    k250ActivateDesignAvailable <- function() {
        k233ActivateTab(S, K250_TABS, "design")
        Sys.sleep(2)
        v <- k250DataAvailable(S)
        k233ActivateTab(S, K250_TABS, "setup")
        Sys.sleep(1)
        v
    }
    available <- k250ActivateDesignAvailable()
    k250H1Snap(S, "h1a", "after_confirm", again$since %||% 0, list(reimport_how = again$how,
        confirm_pressed = confirm$pressed, dialog_closed = confirm$closed, design_data_available = available))
    k233StopSession(S)
}

k250H1CaseD2 <- function(G) {
    s <- k250Start(G, "D2")
    if (is.null(s)) return(invisible(NULL))
    S <- s$S
    imp <- k250Import(S, s$R, "h1d2")
    k250Note(S, "h1d2", "import", "how", imp$how)
    if (!isTRUE(imp$ok)) return(k233StopSession(S))
    cancel <- k250Cancel(S, "h1d2")
    k250H1Snap(S, "h1d2", "after_cancel", cancel$since, cancel$extra)
    again <- k250Import(S, s$R, "h1d2", first = FALSE, tick = FALSE)
    Sys.sleep(3)
    k250H1Snap(S, "h1d2", "after_reimport", again$since %||% 0, list(import_completed = isTRUE(again$completed),
        reimport_how = again$how))
    k233StopSession(S)
}

k250Plan <- function(G) {
    for (case in list(k250H1CaseA, k250H1CaseD2)) {
        if (k233Remaining(G) <= 120) break
        tryCatch(case(G), error = \(e) k250Note(list(G = G, label = "?"), "?", "driver", "error", conditionMessage(e)))
        k250WriteValues(G)
    }
    obs <- k250WriteValues(G)
    sessions <- vapply(G$session_rows, \(r) isTRUE(r$attached), logical(1))
    writeLines(c(
        sprintf("sessions attached: %d of %d (%s)", sum(sessions), length(sessions), paste(names(sessions), collapse = ",")),
        sprintf("predictions: %d as predicted, %d differ, %d not observed, %d recorded without prediction",
            sum(obs$verdict == "as predicted"), sum(obs$verdict == "DIFFERS"), sum(obs$verdict == "not observed"),
            sum(obs$verdict == "recorded"))
    ), file.path(G$out, "verdict.txt"))
}

if (sys.nframe() == 0L) k250Main(commandArgs(trailingOnly = TRUE))
