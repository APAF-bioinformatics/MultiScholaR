# KIT-251 gate G0 (ticket r1.2, section 10): base values, recorded with no assertion, of ONE real mouse click on the
# proteomics Protein Replicate Filter Apply button, at v0.9-staging 16112917. Observation only: this script is the only
# file the commit under test adds to staging (under docs/qa/workflow-execution/kit/kit-251/); nothing under R/ changes.
#
# Instrument: KIT-233's journey driver (`journey.R` on branch kit-233-journey at cce7b04c, read with `git show`, sourced
# for its functions; KIT-233 N1), in a network namespace with loopback only, test mode, the 8-run prot_dia fixture,
# route S1 (home) to S5 (protein QC), the sub-tabs in the application's order up to Duplicate Removal through the
# driver's own presses. The one change to the driver: the Protein Replicate Filter Apply is pressed by ONE real mouse
# click (chromote Input.dispatchMouseEvent: mouseMoved, mousePressed, mouseReleased at the button's centre) instead of
# the driver's DOM `el.click()`, and the lane's permissive retry for that tab is removed, so that one click is all.
#
#   (a) managed: `packaging/electron/entry.R` with the plain launcher (MSR_LAUNCHER_PROBE "process"), as KIT-233 L.
#   (b) in process: `MultiScholaR::run_app(test_mode = TRUE)` in the R process (no managed application).
#
# Per session, after the click and the step's completion signal (or its end) and 15 s more:
#   - the R log's lines "Protein replicate filter apply observed via module input" and "... via root input"
#     (R/mod_prot_qc_protein_replicate.R:304, :315);
#   - the worker requests for operation protein_replicate_filter: `request.rds` files under the project folder's
#     `.workflow-export/wfx-*/` (R/utils_workflow_execution.R:1868-1869; R/utils_workflow_saved_exports.R:603-608);
#   - the page's events: click events on the button (and how many were trusted), `shiny:inputchanged` for its id;
#   - the notifications shown from the click on (text and type);
#   - the head (test-state digest) and the step's outcome.
#
# PREDICTIONS (by reading, ticket section 2; written before any value was read): (a) module input 1, root input 1,
# worker requests 1, one warning "A protein replicate filter action is already running", head
# protein_replicate_filtered; (b) module input 1, root input 1, worker requests 0, no warning, head
# protein_replicate_filtered (the second call returns already_complete silently). Stop rule 1: (a) worker requests 2.
#
# Usage (run-tests scripts job, from the package root; fetch_refs ["kit-233-journey"]):
#   Rscript docs/qa/workflow-execution/kit/kit-251/g0-base-values.R outer <commit under test> <out dir>
#   Rscript docs/qa/workflow-execution/kit/kit-251/g0-base-values.R journey-inner <phase file>   (the namespace phase)

K251_JOURNEY_COMMIT <- "cce7b04c53bcf94ef9b73223bc7d0932c6fc1394"
K251_JOURNEY_PATH <- "docs/qa/workflow-execution/kit/kit-233/journey.R"
K251_RECORDS <- "docs/qa/workflow-execution/kit/kit-251"
K251_BUTTON <- "proteomics_workflow-quality_control-protein_qc-replicate_filter-apply_protein_replicate_filter"
K251_LOG_MODULE <- "Protein replicate filter apply observed via module input"
K251_LOG_ROOT <- "Protein replicate filter apply observed via root input"
K251_WARNING <- "A protein replicate filter action is already running"

k251Main <- function(args) {
    if (length(args) >= 2L && identical(args[[1L]], "journey-inner")) return(invisible(k251Inner(args[[2L]])))
    if (length(args) < 3L || !identical(args[[1L]], "outer")) {
        stop("usage: g0-base-values.R outer <commit under test> <out dir> | journey-inner <phase file>")
    }
    k251Outer(args[[2L]], args[[3L]])
}

k251SourceJourney <- function(path) {
    source(path, local = globalenv())
    # The driver's records folder is this ticket's (its commit check: the commit under test differs only there).
    assign("K233_RECORDS", K251_RECORDS, envir = globalenv())
    k251Overrides()
}

k251Outer <- function(commit, out_dir) {
    started <- Sys.time()
    pkg <- normalizePath(".")
    out <- file.path(pkg, out_dir)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    work <- file.path(tempdir(), "k251-work")
    dir.create(work, recursive = TRUE, showWarnings = FALSE)
    journey <- file.path(work, "journey.R")
    status <- system2("git", c("show", paste0(K251_JOURNEY_COMMIT, ":", K251_JOURNEY_PATH)), stdout = journey)
    if (!identical(status, 0L) || !file.exists(journey) || file.size(journey) == 0) {
        writeLines(sprintf("instrument: journey.R could not be read at %s (git show exit %s)", K251_JOURNEY_COMMIT,
            status), file.path(out, "stop.txt"))
        return(invisible(NULL))
    }
    k251SourceJourney(journey)
    launch <- k233LaunchSpec()
    env <- k233OuterEnvironment("k251-g0", commit, pkg, launch)
    rows <- c(env$rows, list(list(key = "journey_commit", value = K251_JOURNEY_COMMIT),
        list(key = "journey_sha256", value = k233Sha(journey))))
    k233WriteTsv(rows, file.path(out, "environment.tsv"))
    if (!isTRUE(env$commit_ok)) {
        writeLines("stopped: the commit under test differs from the records commit outside the kit-251 records folder",
            file.path(out, "stop.txt"))
        return(invisible(NULL))
    }
    phase <- list(
        run_mode = "lifecycle", started = as.numeric(started),
        bounds = list(deadline_min = 35, hard_s = 5 * 60, finding_s = 600, poll = 0.5),
        pkg = pkg, out = out, work = work, script = k233ScriptPath(), journey = journey, commit = commit,
        records_commit = env$head,
        inputs = list(id = "prot_dia fixture", sha = k233Sha(file.path(pkg, "tests/testdata/e2e/prot_dia/report.tsv"))),
        pandoc = list(dir = "", source = "not needed (no report step)"), launch = launch,
        user_outer = Sys.info()[["effective_user"]], libpaths = .libPaths(), locale = Sys.getlocale("LC_CTYPE"),
        env = c(R_LIBS = Sys.getenv("R_LIBS"), R_LIBS_USER = Sys.getenv("R_LIBS_USER"),
            R_LIBS_SITE = Sys.getenv("R_LIBS_SITE"), LANG = Sys.getenv("LANG", "C.UTF-8"),
            GITHUB_RUN_ID = Sys.getenv("GITHUB_RUN_ID"), TZ = "UTC", OMP_NUM_THREADS = "1",
            OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", NOT_CRAN = "true")
    )
    phase$env <- phase$env[nzchar(phase$env)]
    phase_file <- file.path(work, "k251-phase.rds")
    saveRDS(phase, phase_file)
    k233RunInner(phase, phase_file)
    for (f in c("stop.txt", "inner-fatal.txt", "verdict.txt")) {
        if (file.exists(file.path(out, f))) cat(f, ":", readLines(file.path(out, f)), sep = "\n")
    }
    if (file.exists(file.path(out, "g0-base-values.tsv"))) {
        print(utils::read.delim(file.path(out, "g0-base-values.tsv"), stringsAsFactors = FALSE)[, 1:3], row.names = FALSE)
    }
    k233Finish(out)
}

k251Inner <- function(phase_file) {
    phase <- readRDS(phase_file)
    .libPaths(phase$libpaths)
    if (length(phase$env)) do.call(Sys.setenv, as.list(phase$env))
    Sys.unsetenv(c("MULTISCHOLAR_TEST_FEATURES", "MULTISCHOLAR_TEST_MODE", "CHROMOTE_CHROME"))
    invisible(tryCatch(Sys.setlocale("LC_ALL", phase$locale), warning = \(w) NULL, error = \(e) NULL))
    options(warn = 1, chromote.timeout = 60)
    k251SourceJourney(phase$journey)
    writeLines(k233Utc(), file.path(phase$out, "inner-started.txt"))
    G <- k233Context(phase)
    # One click is all: no permissive retry on the Protein Replicate Filter tab (the driver retries a refused step).
    for (i in seq_along(G$route$qc$protein_qc)) {
        if (identical(G$route$qc$protein_qc[[i]]$id, "replicate_filter")) G$route$qc$protein_qc[[i]]$permissive <- NULL
    }
    G$k251 <- list()
    on.exit({
        for (S in G$sessions) k233StopSession(S)
        k233WriteAll(G)
        k251WriteValues(G)
    }, add = TRUE)
    fatal <- function(e) {
        G$fatal <- paste(conditionMessage(e), "| class:", paste(class(e), collapse = ","), "| call:",
            paste(deparse(conditionCall(e)), collapse = " "))
        writeLines(G$fatal, file.path(G$out, "inner-fatal.txt"))
    }
    tryCatch(k251Plan(G), error = fatal, skip = fatal)
    invisible(G)
}

k251Plan <- function(G) {
    k251Session(G, "A", "managed (packaging/electron/entry.R, plain launcher)", entry = NULL)
    wrapper <- file.path(G$work, "k251-in-process-entry.R")
    writeLines(c(
        "# KIT-251 G0 (b): the application in this R process with no managed application (run_app(test_mode = TRUE)).",
        "arguments <- commandArgs(trailingOnly = TRUE)",
        "cat('MSR-LAUNCHER {\"launcher\":\"in-process\",\"reason\":\"KIT-251 G0 (b): run_app(test_mode = TRUE)\"}\\n')",
        "flush(stdout())",
        "setwd(arguments[[2L]])",
        "shiny::runApp(MultiScholaR::run_app(test_mode = TRUE), host = '127.0.0.1', port = as.integer(arguments[[1L]]),",
        "    launch.browser = FALSE)"
    ), wrapper)
    k251Session(G, "B", "in process (run_app(test_mode = TRUE))", entry = wrapper)
    verdict <- k251Verdict(G)
    writeLines(verdict, file.path(G$out, "verdict.txt"))
}

k251Session <- function(G, label, mode, entry) {
    assign("K251_ENTRY", entry, envir = globalenv())
    S <- k233Launch(G, label)
    row <- list(mode = mode, session = label, launcher = S$launcher, answered = S$answered, attached = !is.null(S$driver),
        attach_error = S$attach_error)
    if (!isTRUE(S$answered) || is.null(S$driver)) {
        G$k251[[label]] <- c(row, list(outcome = "instrument: the application did not answer or was not attached"))
        k233StopSession(S)
        return(invisible(NULL))
    }
    R <- k233FixtureRun(G, S$project_dir)
    R$label <- paste0("k251_", tolower(label))
    res <- k233RunRoute(S, G$route, R, upto = "S5")
    Sys.sleep(15)
    k233Idle(S, 15)
    G$k251[[label]] <- c(row, k251Observe(S, res))
    k233StopSession(S)
    k251CopyRequests(S, G$k251[[label]])
    invisible(NULL)
}

#' The values read after the click (see the header).
k251Observe <- function(S, res) {
    click <- S$k251_click
    steps <- Filter(\(r) identical(r$session, S$label) && identical(r$step, "S5.5"), S$G$steps)
    step <- if (length(steps)) steps[[1L]] else NULL
    log <- k233ReadLog(S$log)
    others <- list.files(S$run_folder, pattern = "\\.log$", recursive = TRUE, full.names = TRUE)
    other_log <- unlist(lapply(others, k233ReadLog))
    page <- k233Js(S, "window.__k251 || null")
    since <- click$t %||% NA_real_
    notes <- if (is.finite(as.numeric(since))) {
        k233Js(S, sprintf("(window.__k233 ? window.__k233.notes : []).filter(function(n){return n.t >= %.0f;})",
            as.numeric(since)))
    } else {
        list()
    }
    st <- k233State(S, 0)
    requests <- k251Requests(S)
    list(
        clicked = !is.null(click), click = click, step_outcome = step$outcome %||% "not reached",
        step_message = step$message %||% "", step_warnings = step$warnings %||% "",
        log_module = sum(grepl(K251_LOG_MODULE, log, fixed = TRUE)),
        log_root = sum(grepl(K251_LOG_ROOT, log, fixed = TRUE)),
        log_module_run_folder = sum(grepl(K251_LOG_MODULE, other_log, fixed = TRUE)),
        log_root_run_folder = sum(grepl(K251_LOG_ROOT, other_log, fixed = TRUE)),
        log_lines = paste(grep("Protein replicate filter apply observed", log, value = TRUE, fixed = TRUE), collapse = " || "),
        page = page, notes = notes, head = k233DigestState(st %||% list()),
        history = paste(as.character(unlist(st$dg$r6_state_history_per_omic$proteomics)), collapse = " > "),
        requests = requests
    )
}

#' The worker requests the controller wrote under the session's folders (request.rds per request directory).
k251Requests <- function(S) {
    files <- unique(c(
        list.files(S$project_dir, pattern = "^request\\.rds$", recursive = TRUE, all.files = TRUE, full.names = TRUE),
        list.files(S$run_folder, pattern = "^request\\.rds$", recursive = TRUE, all.files = TRUE, full.names = TRUE)))
    rows <- lapply(files, function(f) {
        value <- tryCatch(readRDS(f), error = \(e) list(operation = paste("unreadable:", conditionMessage(e))))
        list(file = f, operation = as.character(value$operation %||% ""),
            request_id = as.character(value$token$request_id %||% ""), mtime = format(file.mtime(f), tz = "UTC"))
    })
    rows
}

k251CopyRequests <- function(S, values) {
    k233WriteTsv(lapply(values$requests, \(r) {
        r$file <- sub(S$G$work, "<work>", r$file, fixed = TRUE)
        r
    }), file.path(S$G$out, sprintf("requests-%s.tsv", S$label)))
}

k251WriteValues <- function(G) {
    rows <- list()
    add <- function(session, mode, quantity, value, evidence = "") {
        rows[[length(rows) + 1L]] <<- list(session = session, mode = mode, quantity = quantity,
            value = k233Clean(value), evidence = evidence)
    }
    for (label in names(G$k251)) {
        v <- G$k251[[label]]
        m <- v$mode
        add(label, m, "launcher line", v$launcher, "MSR-LAUNCHER line of the application's console")
        if (!is.null(v$outcome)) {
            add(label, m, "instrument", v$outcome, v$attach_error)
            next
        }
        p <- v$page
        add(label, m, "driver mouse clicks on Apply (Input.dispatchMouseEvent)", if (isTRUE(v$clicked)) 1L else 0L,
            sprintf("at (%s, %s); element at point is the button: %s", v$click$x %||% "", v$click$y %||% "",
                v$click$hit %||% ""))
        add(label, m, "click events on the button in the page", p$clicks %||% NA, "capture listener on document")
        add(label, m, "trusted click events on the button", p$trusted %||% NA, "event.isTrusted")
        add(label, m, "shiny:inputchanged events for the button's id", p$changed %||% NA,
            paste("values:", paste(unlist(p$values), collapse = ",")))
        add(label, m, "R log: observed via module input", v$log_module,
            sprintf("R/mod_prot_qc_protein_replicate.R:304; run folder logs: %s", v$log_module_run_folder))
        add(label, m, "R log: observed via root input", v$log_root,
            sprintf("R/mod_prot_qc_protein_replicate.R:315; run folder logs: %s", v$log_root_run_folder))
        ops <- vapply(v$requests, \(r) r$operation, character(1))
        add(label, m, "worker requests for protein_replicate_filter", sum(ops == "protein_replicate_filter"),
            sprintf("request.rds files: %d in all (%s); requests-%s.tsv", length(ops),
                paste(names(table(ops)), table(ops), sep = "=", collapse = ", "), label))
        texts <- vapply(v$notes, \(n) sprintf("[%s] %s", n$type %||% "", n$text %||% ""), character(1))
        add(label, m, "notifications from the click on", length(texts), paste(texts, collapse = " || "))
        add(label, m, "'already running' warnings", sum(grepl(K251_WARNING, texts, fixed = TRUE)), K251_WARNING)
        add(label, m, "S5.5 outcome (driver)", v$step_outcome, v$step_message)
        add(label, m, "head after the click (digest)", v$head, v$history)
    }
    k233WriteTsv(rows, file.path(G$out, "g0-base-values.tsv"))
}

k251Verdict <- function(G) {
    a <- G$k251$A
    if (is.null(a) || !is.null(a$outcome) || !isTRUE(a$clicked)) {
        return("instrument: session A (managed) did not reach the click; no base value (VAL-54 instrument defect)")
    }
    n <- sum(vapply(a$requests, \(r) identical(r$operation, "protein_replicate_filter"), logical(1)))
    b <- G$k251$B
    b_state <- if (is.null(b) || !is.null(b$outcome) || !isTRUE(b$clicked)) "B not reached" else "B recorded"
    c(sprintf("A (managed): one click, %d worker request(s) for protein_replicate_filter; %s", n,
        if (n >= 2L) "STOP RULE 1: two worker requests from one click" else "stop rule 1 does not fire"), b_state)
}

# --- the two overrides of the driver -----------------------------------------------------------------------------

K251_COUNTER_JS <- r"---((function (sel, id) {
  var K = window.__k251 = { clicks: 0, trusted: 0, changed: 0, values: [] };
  document.addEventListener('click', function (ev) {
    var e = document.querySelector(sel);
    if (e && (ev.target === e || e.contains(ev.target))) { K.clicks++; if (ev.isTrusted) { K.trusted++; } }
  }, true);
  if (window.jQuery) {
    jQuery(document).on('shiny:inputchanged', function (ev) {
      if (ev.name === id) { K.changed++; K.values.push(String(ev.value)); }
    });
  }
  return true;
})(%s, %s))---"

K251_POINT_JS <- r"---((function (sel) {
  var e = document.querySelector(sel);
  if (!e) { return null; }
  e.scrollIntoView({ block: 'center' });
  var r = e.getBoundingClientRect();
  var x = r.left + r.width / 2, y = r.top + r.height / 2;
  var hit = document.elementFromPoint(x, y);
  return { x: x, y: y, w: r.width, h: r.height, hit: !!hit && (hit === e || e.contains(hit)), now: Date.now() };
})(%s))---"

k251Overrides <- function() {
    k233EntryPathDriver <- get("k233EntryPath", envir = globalenv())
    assign("k233EntryPath", function(launch) {
        entry <- get0("K251_ENTRY", envir = globalenv())
        if (is.character(entry) && length(entry) == 1L) entry else k233EntryPathDriver(launch)
    }, envir = globalenv())
    k233PressDriver <- get("k233Press", envir = globalenv())
    assign("k233Press", function(S, step, target, by = "id", click = TRUE) {
        if (!click || !identical(target, K251_BUTTON)) return(k233PressDriver(S, step, target, by, click))
        ready <- k233PressDriver(S, step, target, by, click = FALSE)
        if (!isTRUE(ready$ok)) return(ready)
        selector <- k233Selector(S, target, by)
        k233Js(S, sprintf(K251_COUNTER_JS, k233Json(selector), k233Json(target)))
        point <- k233Js(S, sprintf(K251_POINT_JS, k233Json(selector)))
        if (!is.list(point)) return(list(ok = FALSE, reason = paste(target, "has no box to click"), test = ready$test))
        Sys.sleep(0.3)
        point <- k233Js(S, sprintf(K251_POINT_JS, k233Json(selector)))
        t <- as.numeric(k233Js(S, "Date.now()") %||% point$now)
        cdp <- S$driver$get_chromote_session()
        sent <- tryCatch({
            cdp$Input$dispatchMouseEvent(type = "mouseMoved", x = point$x, y = point$y)
            cdp$Input$dispatchMouseEvent(type = "mousePressed", x = point$x, y = point$y, button = "left", clickCount = 1)
            cdp$Input$dispatchMouseEvent(type = "mouseReleased", x = point$x, y = point$y, button = "left", clickCount = 1)
            TRUE
        }, error = \(e) conditionMessage(e))
        k233LogAction(S, step, target, "one mouse click (Input.dispatchMouseEvent)")
        S$k251_click <- list(t = t, x = round(point$x, 1), y = round(point$y, 1), hit = isTRUE(point$hit),
            sent = sent)
        if (!isTRUE(sent)) return(list(ok = FALSE, reason = paste("mouse events not sent:", sent), test = ready$test))
        list(ok = TRUE, t = t, test = ready$test)
    }, envir = globalenv())
    invisible(NULL)
}

if (sys.nframe() == 0L) k251Main(commandArgs(trailingOnly = TRUE))
