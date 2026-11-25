# KIT-265 H1: the fix observed in a headless browser on a runner (observation only; the predictions below are written
# before any value is read). It reuses KIT-250's driver (`g0.R`, which sources KIT-233's `k233-journey.R`; both copied
# byte for byte into this folder) and replaces the cases, the predictions and the records folder.
#
#   Rscript docs/qa/workflow-execution/kit/kit-265/h1.R h1 <commit under test>
#
# The prot_dia fixture (400 proteins, one organism, accessions P00001 ...) with FASTA files the driver writes:
#   A  headers without organism fields (OS=/OX= removed); the form set to Mus musculus / 10090 BEFORE ticking
#      'Mixed species FASTA'; ticked: route (i), the import continues and a persistent notice names the organism.
#   B  the option ticked and an empty FASTA file (the processing cannot read it): route (iv), refused; then, in the
#      same session, untick and press 'Process Imported Data' again (the refusal's own advice): it completes.
#   C  the option unticked with the same empty FASTA file: completes as before.
#   D  a valid FASTA that names organisms (OS=/OX= kept) and does not contain the data's proteins (accessions
#      prefixed "ZZ"): route (i)(c), refused (the owner's decision VAL-304: refuse).
# Notices and the error are read from the page's notification record; the error toast lasts 10 s, so the page's toast
# count is read at once after the refusal (`error_toast_in_dom`).

local({
    file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    dir <- if (length(file)) dirname(normalizePath(file[[1L]])) else normalizePath("docs/qa/workflow-execution/kit/kit-265")
    source(file.path(dir, "g0.R"), local = globalenv())
})
K233_RECORDS <- "docs/qa/workflow-execution/kit/kit-265"
K250_OUT <- "k265h1-out"
K250_BOUNDS <- list(deadline_min = 42, hard_s = 5 * 60, finding_s = 600, poll = 0.5)
K265_NOTICE_ID <- "shiny-notification-prot_import_organism_cancelled"

K265_NOTICE_I <- paste0("~\\[warning\\] No organism with a taxonomy ID could be matched to the proteins in the search ",
    "results .* This analysis will use the organism entered in the form: Mus musculus \\(taxon 10090\\)\\. ",
    "To use another organism, untick 'Mixed species FASTA'")
K265_ERROR_IV <- paste0("~\\[error\\] Error: The import was stopped: the FASTA file could not be read or processed, ",
    "so the organisms in it cannot be counted\\. Nothing was imported\\.")
K265_ERROR_C <- paste0("~\\[error\\] Error: The import was stopped: none of the proteins in the search results could ",
    "be found in the FASTA file, so their organisms cannot be counted\\. Nothing was imported\\. Check that this is ",
    "the FASTA file the search used")

K250_PREDICTIONS <- list(
    # A: route (i)
    k250P("h1a", "after_import", "import_completed", "TRUE", "route (i) continues"),
    k250P("h1a", "after_import", "modal_open", "none", "no organism window; the processing window closed"),
    k250P("h1a", "after_import", "setup_import", "complete", "the import completes on the form's organism"),
    k250P("h1a", "after_import", "design_matrix", "pending", "completeOmicImport"),
    k250P("h1a", "after_import", "setup_log_taxon_id", "10090", "the form's taxon"),
    k250P("h1a", "after_import", "setup_log_organism", K250_FORM_NAME, "the form's organism"),
    k250P("h1a", "after_import", "import_status", K250_OK, "the import page"),
    k250P("h1a", "after_import", "notices", K265_NOTICE_I, "N1 names Mus musculus (taxon 10090)"),
    k250P("h1a", "after_import", "notice_error", "none", "no error"),
    k250P("h1a", "after_import", "w1_count", "1", "one copy under the fixed id"),
    k250P("h1a", "after_wait", "w1_count", "1", "persistent: still shown more than 10 s later"),
    k250P("h1a", "design", "design_tab_active", "TRUE", "the Design tab opens"),
    k250P("h1a", "design", "data_available", "TRUE", "the import's table is on the bus"),
    # B: route (iv), then the advised path
    k250P("h1b", "after_refusal", "modal_open", "none", "the processing window is closed"),
    k250P("h1b", "after_refusal", "setup_import", "incomplete", "the import is cleared"),
    k250P("h1b", "after_refusal", "design_matrix", "disabled", "no status written"),
    k250P("h1b", "after_refusal", "later_tabs", "all disabled", "nothing ran"),
    k250P("h1b", "after_refusal", "setup_log_taxon_id", "(absent)", "the setup log cleared"),
    k250P("h1b", "after_refusal", "notice_error", K265_ERROR_IV, "R-iv"),
    k250P("h1b", "after_refusal", "error_toast_in_dom", "TRUE", "the 10 s error toast is on the page when read"),
    k250P("h1b", "after_refusal", "import_status_success", "FALSE", "no success text"),
    k250P("h1b", "after_refusal", "notice_success", "absent", "no success notice"),
    k250P("h1b", "after_untick", "import_completed", "TRUE", "the refusal's advice works in the same session"),
    k250P("h1b", "after_untick", "modal_open", "none", "closed"),
    k250P("h1b", "after_untick", "setup_import", "complete", "completes"),
    k250P("h1b", "after_untick", "design_matrix", "pending", "completeOmicImport"),
    k250P("h1b", "after_untick", "notice_error", "none", "no error"),
    k250P("h1b", "after_untick", "form_mixed_species", "FALSE", "unticked"),
    k250P("h1b", "after_untick", "design_data_available", "TRUE", "the Design tab is usable"),
    # C: unticked with the same file: as today
    k250P("h1c", "after_import", "import_completed", "TRUE", "the FASTA failure still logs and continues"),
    k250P("h1c", "after_import", "setup_import", "complete", "as today"),
    k250P("h1c", "after_import", "design_matrix", "pending", "as today"),
    k250P("h1c", "after_import", "notice_error", "none", "as today"),
    k250P("h1c", "after_import", "import_status", K250_OK, "as today"),
    # D: route (i)(c), refuse (VAL-304)
    k250P("h1d", "after_refusal", "modal_open", "none", "the processing window is closed"),
    k250P("h1d", "after_refusal", "setup_import", "incomplete", "the import is cleared"),
    k250P("h1d", "after_refusal", "design_matrix", "disabled", "no status written"),
    k250P("h1d", "after_refusal", "notice_error", K265_ERROR_C, "R-c"),
    k250P("h1d", "after_refusal", "error_toast_in_dom", "TRUE", "the 10 s error toast is on the page when read"),
    k250P("h1d", "after_refusal", "import_status_success", "FALSE", "no success text"),
    k250P("h1d", "after_refusal", "notice_success", "absent", "no success notice")
)

k265Main <- function(args) {
    if (length(args) < 2L) stop("usage: h1.R h1 <commit under test> | journey-inner <phase file>")
    if (identical(args[[1L]], "journey-inner")) return(invisible(k250Inner(args[[2L]])))
    if (!identical(args[[1L]], "h1")) stop("unknown run mode ", args[[1L]])
    k265Outer(args[[2L]])
}

# g0.R's outer phase, with this file as the inner script.
k265Outer <- function(commit) {
    started <- Sys.time()
    pkg <- normalizePath(".")
    out <- file.path(pkg, K250_OUT)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    work <- file.path(tempdir(), "k265-work")
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
        pandoc = list(dir = "", source = "not used (H1 renders no report)"), launch = launch,
        user_outer = Sys.info()[["effective_user"]], libpaths = .libPaths(), locale = Sys.getlocale("LC_CTYPE"),
        env = c(R_LIBS = Sys.getenv("R_LIBS"), R_LIBS_USER = Sys.getenv("R_LIBS_USER"),
            R_LIBS_SITE = Sys.getenv("R_LIBS_SITE"), LANG = Sys.getenv("LANG", "C.UTF-8"),
            GITHUB_RUN_ID = Sys.getenv("GITHUB_RUN_ID"), TZ = "UTC", OMP_NUM_THREADS = "1",
            OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", NOT_CRAN = "true")
    )
    phase$env <- phase$env[nzchar(phase$env)]
    phase_file <- file.path(work, "k265-phase.rds")
    saveRDS(phase, phase_file)
    k233RunInner(phase, phase_file)
    k233Finish(out)
}

# ---------------------------------------------------------------------------------------------------------------
# The FASTA files and the readings.
# ---------------------------------------------------------------------------------------------------------------
k265Fasta <- function(R, kind, dir) {
    path <- file.path(dir, sprintf("k265_%s.fasta", kind))
    lines <- readLines(R$fasta, warn = FALSE)
    header <- startsWith(lines, ">")
    if (identical(kind, "no_organisms")) lines[header] <- sub(" OS=.*$", "", lines[header])
    if (identical(kind, "other_proteins")) lines[header] <- sub("^(>(sp|tr)\\|)", "\\1ZZ", lines[header])
    if (identical(kind, "empty")) lines <- character()
    if (length(lines)) writeLines(lines, path) else file.create(path)
    path
}

# From KIT-250's h1.R (not part of g0.R): the design builder's `data_available` output value.
k250DataAvailable <- function(S) {
    v <- k233Js(S, sprintf("(window.Shiny&&Shiny.shinyapp&&Shiny.shinyapp.$values)?Shiny.shinyapp.$values[%s]:null",
        k233Json(K250_DATA_AVAILABLE)))
    if (is.null(v)) "(no value)" else v
}

k265ErrorToasts <- function(S) {
    as.character(k233Js(S, "document.querySelectorAll('.shiny-notification-error').length") %||% "(no answer)")
}

k265Snap <- function(S, case, moment, since, extra = list()) {
    fields <- k250Snap(S, case, moment, since, c(list(w1_count = as.character(k233Js(S,
        sprintf("document.querySelectorAll('#%s').length", K265_NOTICE_ID)) %||% "(no answer)")), extra))
    k250Note(S, case, moment, "import_status_success", grepl("Data imported successfully", fields$import_status, fixed = TRUE))
    k250WriteValues(S$G)
    invisible(fields)
}

#' 'Process Imported Data' as a user does it: the two files (first import only), the form's organism set to
#' 10090 / Mus musculus BEFORE the option is ticked (first import only), the option set, then the press. Waits for the
#' import's completion line in the application's log, or for an error toast (a refusal).
k265Import <- function(S, R, case, fasta, first = TRUE, tick = TRUE) {
    k233ActivateTab(S, K250_TABS, "setup")
    if (first) {
        up <- k233Upload(S, case, paste0(K250_NS, "search_results_standard"), R$report)
        if (!up$ok) return(list(ok = FALSE, how = paste("search results upload:", up$text)))
        k233WaitFor(S, sprintf("(document.getElementById(%s)||{innerText:''}).innerText.trim().length > 0",
            k233Json(paste0(K250_NS, "format_detection"))), 60)
        up <- k233Upload(S, case, paste0(K250_NS, "fasta_file_standard"), fasta)
        if (!up$ok) return(list(ok = FALSE, how = paste("FASTA upload:", up$text)))
        k233Set(S, case, paste0(K250_NS, "taxon_id"), K250_FORM_TAXON)
        k233Set(S, case, paste0(K250_NS, "organism_name"), K250_FORM_NAME)
    }
    k233Set(S, case, paste0(K250_NS, "mixed_species_fasta"), tick)
    k233Idle(S, 10)
    done0 <- k250LogCount(S)
    # A refusal is an error toast NEW since the press: an earlier one (10 s) may still be on the page.
    toasts0 <- suppressWarnings(as.integer(k265ErrorToasts(S)))
    since <- as.numeric(k233Js(S, "Date.now()") %||% 0)
    p <- k233Press(S, case, "prot-import-process", "testid")
    if (!isTRUE(p$ok)) return(list(ok = FALSE, how = paste("Process Imported Data:", p$reason), since = since))
    started <- Sys.time()
    repeat {
        completed <- k250LogCount(S) > done0
        toasts <- suppressWarnings(as.integer(k265ErrorToasts(S)))
        refused <- !is.na(toasts) && !is.na(toasts0) && toasts > toasts0
        if (completed || refused || k233Seconds(started) > 300 || !S$proc$is_alive()) break
        Sys.sleep(0.5)
    }
    list(ok = completed || refused, completed = completed, refused = refused, since = since,
        how = if (completed) "the import's completion line logged" else if (refused) "an error toast shown" else
            "neither a completion line nor an error within 300 s")
}

k265Reached <- function(S, case, moment, imp, extra = list()) {
    k265Snap(S, case, moment, imp$since %||% 0, c(list(import_completed = isTRUE(imp$completed),
        import_how = imp$how, error_toast_in_dom = k265ErrorToasts(S) != "0"), extra))
}

k265Design <- function(S) {
    k233ActivateTab(S, K250_TABS, "design")
    Sys.sleep(3)
    list(active = k250DesignTabActive(S), available = k250DataAvailable(S))
}

k265CaseA <- function(G) {
    s <- k250Start(G, "A")
    if (is.null(s)) return(invisible(NULL))
    S <- s$S
    fasta <- k265Fasta(s$R, "no_organisms", tempdir())
    imp <- k265Import(S, s$R, "h1a", fasta)
    k250Note(S, "h1a", "import", "how", imp$how)
    if (!isTRUE(imp$ok)) return(k233StopSession(S))
    k265Reached(S, "h1a", "after_import", imp)
    Sys.sleep(12)
    k265Snap(S, "h1a", "after_wait", imp$since, list())
    since <- as.numeric(k233Js(S, "Date.now()") %||% 0)
    d <- k265Design(S)
    k265Snap(S, "h1a", "design", since, list(design_tab_active = d$active, data_available = d$available))
    k233StopSession(S)
}

k265CaseB <- function(G) {
    s <- k250Start(G, "B")
    if (is.null(s)) return(invisible(NULL))
    S <- s$S
    fasta <- k265Fasta(s$R, "empty", tempdir())
    imp <- k265Import(S, s$R, "h1b", fasta)
    k250Note(S, "h1b", "import", "how", imp$how)
    if (!isTRUE(imp$ok)) return(k233StopSession(S))
    Sys.sleep(1)
    k265Reached(S, "h1b", "after_refusal", imp)
    again <- k265Import(S, s$R, "h1b", fasta, first = FALSE, tick = FALSE)
    Sys.sleep(3)
    available <- if (isTRUE(again$completed)) {
        d <- k265Design(S)
        k233ActivateTab(S, K250_TABS, "setup")
        d$available
    } else "(not completed)"
    k265Reached(S, "h1b", "after_untick", again, list(design_data_available = available))
    k233StopSession(S)
}

k265CaseC <- function(G) {
    s <- k250Start(G, "C")
    if (is.null(s)) return(invisible(NULL))
    S <- s$S
    fasta <- k265Fasta(s$R, "empty", tempdir())
    imp <- k265Import(S, s$R, "h1c", fasta, tick = FALSE)
    k250Note(S, "h1c", "import", "how", imp$how)
    if (!isTRUE(imp$ok)) return(k233StopSession(S))
    Sys.sleep(3)
    k265Reached(S, "h1c", "after_import", imp)
    k233StopSession(S)
}

k265CaseD <- function(G) {
    s <- k250Start(G, "D")
    if (is.null(s)) return(invisible(NULL))
    S <- s$S
    fasta <- k265Fasta(s$R, "other_proteins", tempdir())
    imp <- k265Import(S, s$R, "h1d", fasta)
    k250Note(S, "h1d", "import", "how", imp$how)
    if (!isTRUE(imp$ok)) return(k233StopSession(S))
    Sys.sleep(1)
    k265Reached(S, "h1d", "after_refusal", imp)
    k233StopSession(S)
}

k250Plan <- function(G) {
    for (case in list(k265CaseA, k265CaseB, k265CaseC, k265CaseD)) {
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

if (sys.nframe() == 0L) k265Main(commandArgs(trailingOnly = TRUE))
