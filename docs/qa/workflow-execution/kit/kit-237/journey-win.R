# KIT-237 (r1.2): the Windows wrapper around KIT-233's journey machinery (docs/qa/workflow-execution/kit/kit-233/
# journey.R, fetched at the commit the request names and sourced, never copied or edited). Records only: nothing under
# R/ is changed; what the journey finds is reported, never fixed.
#
# Why a wrapper (KIT-237 section 10's fallback, the smallest form): journey.R is sourceable as functions (N1), takes the
# launch as a parameter with an installed application folder and a recorded network method (N2) and runs steps from a
# route table (N3), but its entry points are the Linux run-tests ones: k233Outer checks the records commit by `git diff`
# and closes the network with `unshare`; k233Context fixes the route to the DIA table; k233Lifecycle is the DIA
# lifecycle; k233Launch starts R.home()'s Rscript in the driver's environment (so the application would inherit the
# driver's library variables); k233KillAppHook ends the parent process only (on Windows the installed R\bin\Rscript.exe
# serves through its R\bin\x64\Rscript.exe child); and the page helpers read the proteomics namespace only. This file
# replaces exactly those parts, by defining functions of its own and rebinding four of journey.R's names after sourcing
# it (k233Action, k233StopSession, k233ProjectPaths, K233_MARK_PANE_JS per session). Each rebinding is listed in
# k237Rebind().
#
# Usage (distributable.yml, stage journey, under the installed bundle's R through packaging/verify/windows-clean-run.ps1,
# so this R and the application have the environment the shell gives R):
#   Rscript --vanilla journey-win.R <request.json> <application folder> <test tools library> <out folder>
#       <source folder> <fixture folder> <checkout>
# <source folder> holds docs/qa/workflow-execution/kit/kit-233/ and kit-237/ as fetched by `git archive` from the
# request's commits; <fixture folder> holds tests/testdata/e2e/ from the commit under test (no line-ending conversion).

K237_RECORDS <- "docs/qa/workflow-execution/kit/kit-237"
K237_PROBE <- '{"launcher":"process","reason":"win32: no systemd"}'
K237_WORKER_LOG <- "(^|/)worker-[^/]*\\.log$"

# ---------------------------------------------------------------------------------------------------------------
# Small utilities (Windows forms).
# ---------------------------------------------------------------------------------------------------------------
k237Slash <- function(x) gsub("\\\\", "/", x)

k237PortAnswers <- function(port) {
    isTRUE(tryCatch({
        curl::curl_fetch_memory(sprintf("http://127.0.0.1:%d", port),
            handle = curl::new_handle(timeout = 3, connecttimeout = 3))
        TRUE
    }, error = function(e) FALSE))
}

#' Seconds until the port refuses connections, or NA when it still answers after `limit` seconds.
k237WaitRefused <- function(port, limit = 10) {
    started <- Sys.time()
    repeat {
        if (!k237PortAnswers(port)) return(round(k233Seconds(started), 1))
        if (k233Seconds(started) > limit) return(NA_real_)
        Sys.sleep(0.5)
    }
}

#' The whole process tree of an application: `taskkill /PID <pid> /T /F` as the shell does
#' (packaging/electron/main.js:216-217), then processx's marked-tree kill for anything left.
k237KillTree <- function(proc) {
    pid <- tryCatch(proc$get_pid(), error = function(e) NA_integer_)
    said <- if (is.na(pid)) "no pid" else tryCatch(paste(suppressWarnings(system2("taskkill",
        c("/PID", as.character(pid), "/T", "/F"), stdout = TRUE, stderr = TRUE)), collapse = " | "),
        error = function(e) paste("taskkill:", conditionMessage(e)))
    left <- tryCatch(length(proc$kill_tree()), error = function(e) NA_integer_)
    sprintf("taskkill: %s; kill_tree ended %s further processes", said, left)
}

#' File listing of a folder: relative path, bytes and modification time.
k237Listing <- function(dir) {
    if (is.null(dir) || !dir.exists(dir)) return(data.frame(path = character(), bytes = numeric(), mtime = numeric()))
    files <- list.files(dir, recursive = TRUE, all.files = TRUE, no.. = TRUE)
    info <- file.info(file.path(dir, files))
    data.frame(path = files, bytes = info$size, mtime = as.numeric(info$mtime), stringsAsFactors = FALSE)
}

k237Changed <- function(before, after) {
    key <- function(x) paste(x$path, x$bytes, x$mtime)
    after$path[!key(after) %in% key(before)]
}

#' The TMP folder's names to depth 2 (the application's TMP is the runner's; its tempdir is a child of it).
k237TmpNames <- function() {
    tmp <- Sys.getenv("TMP", dirname(tempdir()))
    top <- list.files(tmp, all.files = TRUE, no.. = TRUE)
    inner <- unlist(lapply(top[dir.exists(file.path(tmp, top))], function(d) file.path(d, list.files(file.path(tmp, d),
        all.files = TRUE, no.. = TRUE))))
    sort(c(top, inner))
}

# ---------------------------------------------------------------------------------------------------------------
# The launch (5.1): the installed application's entry.R with the installed R\bin\Rscript.exe, the run folder and the
# driver's pid as its three trailing arguments; the environment is this process's (the shell's, from
# windows-clean-run.ps1) plus the probe, the run folder and one thread per process. A wrapper file (KIT-233 5.1's
# recorded fallback form) sets the test-mode options and writes one probe line, the application's .libPaths(), into the
# run folder before it sources entry.R; in user mode it writes the probe line only.
# ---------------------------------------------------------------------------------------------------------------
k237Launch <- function(G, label, test_mode = TRUE) {
    S <- new.env(parent = emptyenv())
    S$G <- G
    S$label <- label
    S$test_mode <- test_mode
    S$nn <- 0L
    S$hard_s <- G$bounds$hard_s
    S$poll <- G$bounds$poll
    S$hooks <- list()
    S$console_seen <- 0L
    S$port <- 18237L + length(G$sessions)
    S$run_folder <- file.path(G$work, paste0("run-", label))
    S$project_dir <- file.path(G$work, paste0("projects-", label))
    dir.create(S$run_folder, recursive = TRUE, showWarnings = FALSE)
    dir.create(S$project_dir, recursive = TRUE, showWarnings = FALSE)
    S$log <- file.path(G$out, sprintf("app-%s.log", label))
    G$sessions[[label]] <- S
    entry <- k233EntryPath(G$phase$launch)
    rscript <- file.path(G$phase$launch$app_dir, "resources", "runtime", "R", "bin", "Rscript.exe")
    wrapper <- file.path(G$work, sprintf("k237-entry-%s.R", label))
    writeLines(c(
        if (test_mode) K233_TEST_OPTIONS,
        "writeLines(.libPaths(), file.path(commandArgs(trailingOnly = TRUE)[[2L]], 'k237-libpaths.txt'))",
        sprintf('source("%s")', k237Slash(entry))), wrapper)
    tail_args <- c(as.character(S$port), k237Slash(S$run_folder), as.character(Sys.getpid()))
    args <- c("--vanilla", k237Slash(wrapper), tail_args)
    env <- c("current", MSR_LAUNCHER_PROBE = K237_PROBE, MSR_RUN_FOLDER = k237Slash(S$run_folder),
        OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", CHROMOTE_CHROME = K233_NO_BROWSER)
    if (test_mode) env <- c(env, MULTISCHOLAR_TEST_MODE = "true")
    # The request's `test_features` (KIT-240 section 8, addendum 3): TRUE starts the application with
    # MULTISCHOLAR_TEST_FEATURES=1 (packaging/electron/entry.R:14-31); FALSE or absent leaves it unset.
    if (isTRUE(G$test_features)) env <- c(env, MULTISCHOLAR_TEST_FEATURES = "1")
    S$command <- paste(c(rscript, args), collapse = " ")
    started <- Sys.time()
    S$proc <- processx::process$new(rscript, args, env = env, wd = G$pkg, stdout = S$log, stderr = "2>&1",
        cleanup = TRUE, cleanup_tree = TRUE)
    S$url <- sprintf("http://127.0.0.1:%d", S$port)
    S$answered <- FALSE
    for (i in seq_len(240L)) {
        if (!S$proc$is_alive()) break
        S$answered <- isTRUE(tryCatch(curl::curl_fetch_memory(S$url, handle = curl::new_handle(timeout = 5))$status_code ==
            200L, error = function(e) FALSE))
        if (S$answered) break
        Sys.sleep(1)
    }
    S$answer_s <- round(k233Seconds(started), 1)
    log <- k233ReadLog(S$log)
    S$launcher <- paste(grep("^MSR-LAUNCHER ", log, value = TRUE), collapse = " ")
    S$usage_error <- any(grepl("usage: Rscript --vanilla entry.R", log, fixed = TRUE))
    S$driver <- NULL
    S$attach_error <- ""
    if (S$answered) {
        refused <- function(e) {
            S$attach_error <- conditionMessage(e)
            NULL
        }
        S$driver <- tryCatch(
            shinytest2::AppDriver$new(S$url, load_timeout = 120000, timeout = 60000, width = 1440, height = 1100),
            error = refused, skip = refused)
    }
    S$recorder <- !is.null(S$driver) && isTRUE(k233Js(S, K233_RECORDER_JS))
    G$session_rows[[label]] <- list(session = label, mode = if (test_mode) "test" else "user", port = S$port,
        pid = S$proc$get_pid(), command = S$command, answered = S$answered, answer_seconds = S$answer_s,
        launcher = S$launcher, usage_error = S$usage_error, attached = !is.null(S$driver),
        attach_error = S$attach_error, recorder = S$recorder, started_utc = k233Utc(started))
    k233WriteAll(G)
    S
}

#' Every stop of an application (between sessions, at the end, on error): the browser session, then the tree kill, then
#' the port is confirmed to refuse (WL-3's method everywhere, 6).
k237StopSession <- function(S) {
    if (isTRUE(S$stopped)) return(invisible(NULL))
    S$stopped <- TRUE
    k233CopyLogs(S)
    if (!is.null(S$driver)) tryCatch(S$driver$stop(), error = function(e) NULL)
    proc <- S$real_proc %||% S$proc
    said <- k237KillTree(proc)
    refused <- k237WaitRefused(S$port, 10)
    k233CopyLogs(S)
    G <- S$G
    if (!is.null(G$session_rows[[S$label]])) {
        G$session_rows[[S$label]]$exit_status <- tryCatch(proc$get_exit_status(), error = function(e) NA)
        G$session_rows[[S$label]]$stopped_utc <- k233Utc()
        G$session_rows[[S$label]]$stop <- said
        G$session_rows[[S$label]]$port_refused_after_s <- refused
    }
    invisible(NULL)
}

#' WL-4 (a) on metab_lc (SC-05): after S8 completes, the table's "Show:" set to "All Results" (a user's choice; with
#' no significant feature the default view is empty, WL rerun 1), its rows awaited, then the body emptied in the page
#' and the named outputs read again; a kind-3 finding naming the table is required.
k237L4aOmic <- function(S, R, route) {
    da <- sprintf("%s_workflow-de-da_results_table", route$omic)
    k233Set(S, "S8", sprintf("%s_workflow-de-table_significance", route$omic), "all")
    k233WaitFor(S, sprintf(paste0("(function(){var e=document.getElementById(%s);if(!e){return false;}",
        "return e.querySelectorAll('tbody tr').length>0 && !e.querySelector('td.dataTables_empty');})()"), k233Json(da)), 30)
    emptied <- k233Js(S, sprintf(paste0("(function(){var e=document.getElementById(%s);if(!e){return false;}",
        "var b=e.querySelectorAll('tbody');b.forEach(function(x){x.innerHTML='';});return b.length>0;})()"), k233Json(da)))
    outputs <- k233CheckOutputs(S, c(sprintf("%s_workflow-de-da_summary_stats", route$omic), da))
    empty <- vapply(Filter(function(o) !isTRUE(o$nonempty), outputs), function(o) o$id, character(1))
    for (o in empty) k233Finding(S, "S8", 3L, paste("named output empty:", o, "(WL-4 (a) on metab_lc, planted)"), planted = TRUE)
    k237Criterion(S$G, "WL-4", "(a) on metab_lc: the emptied DA results table gives a kind-3 finding naming it",
        da %in% empty, sprintf("tbody emptied: %s; empty outputs at the second check: %s", k233Clean(emptied),
            paste(empty, collapse = ", ")))
}

#' WL-3's crash: the whole tree ended while S2 waits; the port must refuse within 10 s. The processx handle is then
#' replaced by one that reports the process alive, so that the step's `crashed` can come only from the browser session
#' (its disconnection or the page no longer answering), never from the handle's exit alone (6, WL-3).
k237KillTreeHook <- function(S, el) {
    if (isTRUE(S$killed)) return(invisible(NULL))
    S$killed <- TRUE
    S$kill_at <- round(el, 1)
    S$kill_said <- k237KillTree(S$proc)
    S$port_refused_s <- k237WaitRefused(S$port, 10)
    S$real_proc <- S$proc
    real <- S$proc
    S$proc <- list(is_alive = function() TRUE, get_exit_status = function() real$get_exit_status(),
        get_pid = function() real$get_pid(), kill = function() NULL, kill_tree = function() NULL)
    invisible(NULL)
}

# ---------------------------------------------------------------------------------------------------------------
# One action (journey.R's k233Action, kept whole and called first), then this ticket's records for it: the worker units
# whose logs appear in the run folder during the action (7: `steps.tsv` names them, or `controller`); the files the
# route names for the action, each with size and SHA-256 (a missing or empty one is a kind-3 finding, 5.2), and every
# other file the action wrote into the project; the TMP listing around the Summary save (S10, D-1).
# ---------------------------------------------------------------------------------------------------------------
k237Action <- function(S, spec) {
    G <- S$G
    base <- k237Base(S)
    runs0 <- k237Listing(S$run_folder)
    proj0 <- k237Listing(base)
    summary_save <- identical(spec$step, paste0(S$R$summary_step %||% "S10", ".1"))
    if (summary_save) S$tmp_before <- k237TmpNames()
    # RA-02: a screenshot at the moment a terminal signal appears (an error notification, or a modal that is neither
    # expected nor a progress modal), before journey.R's 60-second grace lets it disappear.
    S$terminal_shot <- ""
    since_page <- as.numeric(k233Js(S, "Date.now()") %||% 0)
    inner <- spec$hook %||% S$hooks[[spec$step]]
    spec$hook <- function(S, el) {
        if (is.function(inner)) inner(S, el)
        if (nzchar(S$terminal_shot)) return(invisible(NULL))
        # Only signals of this action (WL rerun 2: earlier steps' error notes gave terminal screenshots at 0 s).
        seen <- k233Js(S, sprintf(paste0("(function(){var K=window.__k233||{notes:[],modals:[]};var t=%.0f;",
            "var e=K.notes.filter(function(n){return n.type==='error'&&n.t>=t;}).length>0;",
            "var m=K.modals.filter(function(x){return x.end===null&&!x.progress&&x.start>=t;});",
            "return {error:e, modal:m.length?m[m.length-1].title:''};})()"), since_page))
        modal <- as.character(seen$modal %||% "")
        if (isTRUE(seen$error) || (nzchar(modal) && !modal %in% spec$modals)) {
            path <- file.path(G$out, sprintf("%s-%02d-%s-terminal.png", S$label, S$nn, spec$step))
            ok <- tryCatch({
                S$driver$get_screenshot(path, selector = "viewport")
                TRUE
            }, error = function(e) FALSE)
            S$terminal_shot <- if (ok) sprintf("%s at %.1f s", basename(path), el) else sprintf("failed at %.1f s", el)
        }
        invisible(NULL)
    }
    res <- k233ActionKit233(S, spec)
    i <- length(G$steps)
    G$steps[[i]]$terminal_screenshot <- S$terminal_shot
    changed <- k237Changed(runs0, k237Listing(S$run_folder))
    workers <- grep(K237_WORKER_LOG, k237Slash(changed), value = TRUE)
    G$steps[[i]]$workers <- if (length(workers)) paste(sub("\\.log$", "", basename(workers)), collapse = " ") else
        "controller"
    G$steps[[i]]$run_folder_changes <- paste(k237Slash(changed), collapse = " ")
    last_summary <- identical(spec$step, paste0(S$R$summary_step %||% "S10", ".5"))
    if (summary_save || last_summary) k237TmpCheck(S, spec$step)
    if (identical(res$outcome, "complete")) {
        base <- k237Base(S)
        G$steps[[i]]$named_files <- tryCatch(k237FileCheck(S, spec$step, base, proj0, G$steps[[i]]),
            error = function(e) paste("driver error in the file check:", conditionMessage(e)))
    }
    k233WriteAll(G)
    res
}

k237Base <- function(S) {
    R <- S$R
    if (is.null(R) || is.null(R$live$project_dir)) return(NULL)
    file.path(R$live$project_dir, R$label)
}

k237TmpCheck <- function(S, step) {
    after <- k237TmpNames()
    before <- S$tmp_before %||% character()
    writeLines(before, file.path(S$G$out, sprintf("%s-tmp-before-summary.txt", S$label)))
    writeLines(after, file.path(S$G$out, sprintf("%s-tmp-after-%s.txt", S$label, step)))
    left <- grep("(^|/)archive-stage-", setdiff(after, before), value = TRUE)
    for (l in left) {
        k233Finding(S, step, "D-1", paste("a staging folder was left in TMP after the Summary save:", l))
    }
}

#' The named files of one action (route$files[[step]]: regular expressions on the path relative to the project base),
#' and every other file the action wrote there. Written to `<label>-<nn>-<step>-files.tsv`.
k237FileCheck <- function(S, step, base, before, row) {
    if (is.null(base) || !dir.exists(base)) return("no project folder")
    patterns <- S$R$files[[step]] %||% character()
    plant <- S$plant_delete
    if (!is.null(plant) && identical(plant$step, step)) {
        # Every copy the pattern names is deleted (WL 45605975875: one of several copies deleted left the pattern met).
        hit <- grep(plant$pattern, k237Slash(k237Listing(base)$path), value = TRUE)
        S$planted_deleted <- paste(hit, collapse = " | ")
        for (h in hit) unlink(file.path(base, h))
    }
    after <- k237Listing(base)
    rel <- k237Slash(after$path)
    rows <- list()
    named <- character()
    for (p in patterns) {
        hits <- which(grepl(p, rel))
        planted <- !is.null(plant) && identical(plant$step, step) && identical(plant$pattern, p)
        if (!length(hits)) {
            rows[[length(rows) + 1L]] <- list(pattern = p, path = "", bytes = NA, sha256 = "", status = "missing")
            k233Finding(S, step, 3L, sprintf("named file missing after the step: %s (under %s)", p, base), row, planted)
            next
        }
        for (h in hits) {
            named <- c(named, rel[[h]])
            empty <- is.na(after$bytes[[h]]) || after$bytes[[h]] <= 0
            rows[[length(rows) + 1L]] <- list(pattern = p, path = rel[[h]], bytes = after$bytes[[h]],
                sha256 = k233Sha(file.path(base, after$path[[h]])), status = if (empty) "empty" else "ok")
            if (empty) k233Finding(S, step, 3L, paste("named file empty after the step:", rel[[h]]), row, planted)
        }
    }
    for (f in setdiff(k237Slash(k237Changed(before, after)), named)) {
        full <- file.path(base, f)
        rows[[length(rows) + 1L]] <- list(pattern = "(not named; written during the step)", path = f,
            bytes = file.size(full), sha256 = k233Sha(full), status = "new")
    }
    k233WriteTsv(rows, file.path(S$G$out, sprintf("%s-%02d-%s-files.tsv", S$label, S$nn, step)))
    ok <- sum(vapply(rows, function(r) identical(r$status, "ok"), logical(1)))
    bad <- sum(vapply(rows, function(r) r$status %in% c("missing", "empty"), logical(1)))
    sprintf("%d named ok, %d missing or empty, %d other new", ok, bad, length(rows) - ok - bad)
}

#' The project folder the application was given or kept (S1 sets `R$live$project_dir`); journey.R's N1 and N3 read it.
k237ProjectPaths <- function(R) {
    dir <- if (is.environment(R$live) && !is.null(R$live$project_dir)) R$live$project_dir else R$project_dir
    base <- file.path(dir, R$label)
    omic <- R$omic %||% "proteomics"
    list(base = base, source = file.path(base, "scripts", omic),
        da = file.path(base, "results", omic, if (identical(omic, "proteomics")) "da_proteins" else "de"))
}

#' Point the page helpers that name the proteomics shell at the route's omic: the active pane for screenshots, and
#' (SC-04) the state script's stepper, `<omic>_workflow-workflow_progress` (R/mod_omic_workflow.R:40, :119).
k237SetOmic <- function(omic, tabs_id) {
    assign("K233_MARK_PANE_JS", gsub("proteomics_workflow-workflow_tabs", tabs_id, K237_MARK_PANE_JS_KIT233,
        fixed = TRUE), envir = globalenv())
    assign("K233_STATE_JS", gsub("#proteomics_workflow-workflow_progress", sprintf("#%s_workflow-workflow_progress", omic),
        K237_STATE_JS_KIT233, fixed = TRUE), envir = globalenv())
}

#' journey.R's chooser click-through, with the folders found (WL 45605975875: shinyFiles 0.9.3 marks a folder
#' `.sF-directory`, a file `.sF-file`, and opens a folder on its own double-click handler, shinyFiles.js:1138, :1155;
#' journey.R looked among `.sF-file` only, so no folder was ever found). On a miss, the names the chooser lists are kept.
k237ChooseFile <- function(S, step, button, path) {
    p <- k233Press(S, step, button)
    if (!isTRUE(p$ok)) return(list(ok = FALSE, how = paste("chooser button:", p$reason)))
    modal <- sprintf("#%s-modal", button)
    if (!k233WaitFor(S, sprintf("document.querySelectorAll(%s).length > 0", k233Json(paste(modal, ".sF-fileList > *"))),
        30)) {
        return(list(ok = FALSE, how = "the file chooser listed nothing within 30 s"))
    }
    shot <- file.path(S$G$out, sprintf("%s-%s-chooser-%s.png", S$label, step, sub(".*-", "", button)))
    tryCatch(S$driver$get_screenshot(shot, selector = paste(modal, ".modal-content")), error = function(e) NULL)
    parts <- c(strsplit(dirname(path), "/", fixed = TRUE)[[1L]][-1L], basename(path))
    rows <- function(m) sprintf("Array.prototype.slice.call(document.querySelectorAll(%s+' .sF-fileList .sF-file, '+%s+' .sF-fileList .sF-directory'))",
        k233Json(m), k233Json(m))
    find <- function(name) sprintf(paste0("(function(){var f=%s.filter(function(x){var d=x.querySelector('.sF-file-name>div');",
        "return d&&d.textContent===%s;})[0];return f||null;})()"), rows(modal), k233Json(name))
    listed <- function() paste(utils::head(as.character(unlist(k233Js(S, sprintf(paste0("%s.map(function(x){",
        "var d=x.querySelector('.sF-file-name>div');return d?d.textContent:'';})"), rows(modal))))), 60), collapse = ", ")
    for (i in seq_along(parts)) {
        if (!k233WaitFor(S, sprintf("(%s) !== null", find(parts[[i]])), 20)) {
            return(list(ok = FALSE, how = sprintf("the chooser does not list %s (it lists: %s)", parts[[i]], listed())))
        }
        event <- if (i < length(parts)) "dblclick" else "click"
        # A native mouse event with button 0, as a user's click gives: shinyFiles selects a row only when
        # event.button === 0 (shinyFiles.js:77-78); jQuery's trigger() carries no button (WL rerun 1, 45619531044).
        k233Js(S, sprintf(paste0("(function(){var f=%s;if(!f){return false;}f.dispatchEvent(new MouseEvent('%s',",
            "{bubbles:true,cancelable:true,button:0,view:window}));return true;})()"), find(parts[[i]]), event))
        k233LogAction(S, step, paste(button, "chooser"), paste(event, parts[[i]]))
    }
    p <- k233Press(S, step, paste(modal, "#sF-selectButton"), "css")
    if (!isTRUE(p$ok)) return(list(ok = FALSE, how = paste("Select:", p$reason)))
    ok <- k233WaitFor(S, sprintf("document.querySelector(%s) === null", k233Json(modal)), 20)
    list(ok = ok, how = "clicked through the chooser's folder list (folders and files)")
}

#' RA-03 and RA-05: the DA table of a contrast by an exact two-sided match. The contrast's two group names (from
#' `group<A>-group<B>` or `<A>_vs_<B>`) must each equal, in order and normalised, the two sides parsed from the file's
#' stem; the da folder is searched first, then the project (visible folders only); more than one candidate is refused.
#' The number of candidates and the path read, relative to the project, are kept as an observation.
k237ContrastPair <- function(x) {
    m <- regmatches(x, regexec("^group(.+)-group(.+)$", x))[[1L]]
    if (length(m) == 3L) return(m[2:3])
    p <- strsplit(x, "_vs_", fixed = TRUE)[[1L]]
    if (length(p) == 2L) p else NULL
}

k237DaFile <- function(R, contrast) {
    want <- k237ContrastPair(contrast)
    pick <- function(files) {
        if (is.null(want)) return(character())
        stems <- sub("^da_proteins_(.*)_long_annot\\.tsv$", "\\1", basename(files))
        keep <- vapply(stems, function(st) {
            got <- k237ContrastPair(st)
            !is.null(got) && identical(k233Norm(got), k233Norm(want))
        }, logical(1))
        files[keep]
    }
    paths <- k233ProjectPaths(R)
    hits <- pick(list.files(paths$da, "^da_proteins_.*_long_annot\\.tsv$", full.names = TRUE))
    where <- "da folder"
    if (!length(hits)) {
        hits <- pick(list.files(paths$base, "^da_proteins_.*_long_annot\\.tsv$", full.names = TRUE, recursive = TRUE))
        where <- "project"
    }
    G <- get0("K237_G", envir = globalenv())
    if (is.environment(G)) {
        rel <- substring(k237Slash(hits), nchar(k237Slash(paths$base)) + 2L)
        k233Observation(G, R$label, sprintf("N2 table for %s (RA-03, RA-05)", contrast), sprintf(
            "%d candidate(s) in the %s by exact two-sided match; read: %s", length(hits), where,
            if (length(hits) == 1L) rel else "none (refused)"))
    }
    if (length(hits) == 1L) hits[[1L]] else NA_character_
}

k237Rebind <- function() {
    env <- globalenv()
    assign("k233ActionKit233", k233Action, envir = env)
    assign("K237_MARK_PANE_JS_KIT233", K233_MARK_PANE_JS, envir = env)
    assign("K237_STATE_JS_KIT233", K233_STATE_JS, envir = env)
    if (!grepl("#proteomics_workflow-workflow_progress", K233_STATE_JS, fixed = TRUE)) {
        stop("journey.R's state script no longer names #proteomics_workflow-workflow_progress: the stepper rebinding ",
            "(SC-04) needs revising")
    }
    assign("k233Action", k237Action, envir = env)
    assign("k233StopSession", k237StopSession, envir = env)
    assign("k233ProjectPaths", k237ProjectPaths, envir = env)
    assign("k233ChooseFile", k237ChooseFile, envir = env)
    assign("k233DaFile", k237DaFile, envir = env)
    c("k233ChooseFile -> k237ChooseFile (folders are .sF-directory)",
        "k233DaFile -> k237DaFile (RA-03 exact two-sided match; RA-05 path recorded)","k233Action -> k237Action (journey.R's action, then this ticket's records)",
        "k233StopSession -> k237StopSession (tree kill, port refusal confirmed)",
        "k233ProjectPaths -> k237ProjectPaths (the folder S1 kept; the route's omic)",
        "K233_MARK_PANE_JS -> the route's workflow tabset (per session)",
        "K233_STATE_JS -> the route's omic stepper (per session, SC-04)")
}

# ---------------------------------------------------------------------------------------------------------------
# The run WL (6): WL-0 to WL-5 on the label-free route, then one journey on each further Windows route (OD-1 (b)), the
# planted refusal, the crash, and the chooser once in user mode (SD-2 answered before WL, 5.2).
# ---------------------------------------------------------------------------------------------------------------
k237Criterion <- function(G, criterion, item, result, evidence = "") k233Criterion(G, criterion, item, result,
    k233Clean(evidence))

k237RunSession <- function(G, route, R, label, test_mode = TRUE, upto = NULL, hooks = list()) {
    k237SetOmic(route$omic, route$tabs_id)
    S <- k237Launch(G, label, test_mode)
    R$live <- new.env(parent = emptyenv())
    S$R <- R
    if (!is.null(S$driver)) k233RunRoute(S, route, R, upto = upto, hooks = hooks)
    S
}

k237WL0 <- function(G, S, wf) {
    crit <- function(item, pass, evidence = "") k237Criterion(G, "WL-0", item, pass, evidence)
    crit("the bundle's commit accepted by 5.4's check (journey step)", identical(wf$commit_check, "pass"),
        paste(wf$commit_check, wf$package_commit, wf$commit_under_test))
    crit("a deliberately wrong commit under test stops at the same check: non-zero exit and val_product_tree_mismatch",
        identical(wf$planted_commit_check, "pass"), paste(wf$planted_commit_check, wf$planted_commit))
    crit("the test tools load shinytest2 and chromote under the installed Rscript",
        requireNamespace("shinytest2", quietly = TRUE) && requireNamespace("chromote", quietly = TRUE),
        paste(R.home(), paste(.libPaths(), collapse = " | ")))
    chrome <- tryCatch(as.character(chromote::find_chrome()), error = function(e) "")
    crit("Chrome found on windows-2025", nzchar(chrome) && file.exists(chrome), chrome)
    crit("entry.R starts with its three trailing arguments and answers HTTP 200",
        isTRUE(S$answered) && !isTRUE(S$usage_error), sprintf("answered=%s after %s s; %s", S$answered, S$answer_s,
            S$command))
    crit("the MSR-LAUNCHER line names process", grepl('"launcher":"process"', S$launcher, fixed = TRUE), S$launcher)
    crit("AppDriver attaches by URL", !is.null(S$driver), S$attach_error)
    crit("firewall rules exist for the installed programs", isTRUE(as.integer(wf$firewall_programs) > 0) &&
        isTRUE(as.integer(wf$firewall_unruled) == 0), sprintf("%s programs, %s without a rule",
            wf$firewall_programs, wf$firewall_unruled))
    long <- tryCatch(utils::readRegistry("SYSTEM\\CurrentControlSet\\Control\\FileSystem", "HLM")$LongPathsEnabled,
        error = function(e) NA)
    G$long_paths <- long
    crit("LongPathsEnabled reads 0 (as read by the installed R)", identical(as.integer(long), 0L), k233Clean(long))
    probe <- k233ReadLog(file.path(S$run_folder, "k237-libpaths.txt"))
    bundle <- k237Slash(normalizePath(file.path(G$phase$launch$app_dir, "resources", "runtime", "R"), mustWork = FALSE))
    inside <- startsWith(tolower(k237Slash(normalizePath(probe, mustWork = FALSE))), tolower(bundle))
    G$app_libpaths <- paste(probe, collapse = " | ")
    crit("the application's .libPaths() is the installed bundle's library only (probe line in the run folder)",
        length(probe) > 0 && all(inside), G$app_libpaths)
    log <- k233ReadLog(S$log)
    if (isTRUE(G$test_features)) {
        crit("test_features TRUE: the application starts with the test setting (one MSR-TEST-FEATURES line)",
            sum(grepl("^MSR-TEST-FEATURES", log)) == 1L, paste(grep("^MSR-", log, value = TRUE), collapse = " || "))
    } else {
        crit("MULTISCHOLAR_TEST_FEATURES unset in the application's environment (no MSR-TEST-FEATURES line)",
            !any(grepl("^MSR-TEST-FEATURES", log)) && !nzchar(Sys.getenv("MULTISCHOLAR_TEST_FEATURES")),
            paste(grep("^MSR-", log, value = TRUE), collapse = " || "))
    }
}

# ---------------------------------------------------------------------------------------------------------------
# WL' (SC-11; ticket 7 and 20): the lifecycle that proves, on the installed build it runs on, what WL proved and RA-01's
# pieces WL never met: the internal deadline (planted), a progress dialog open over 60 s (planted), and the outbound
# block (a real connection attempt from the application's R, made by the journey step; read here). Test mode.
# ---------------------------------------------------------------------------------------------------------------

#' The deadline, planted: while S2 waits, the driver's internal deadline is moved to 1 s in the past (W02: LFQ S2
#' completes 1.5 to 2.1 s after the press, so a deadline 5 s ahead never ended it); the step must end
#' `deadline`, every later step `not_reached (deadline)`, a kind-5 finding must be written, and the tables must exist.
#' The real deadline is restored when the session ends.
k237DeadlineHook <- function(S, el) {
    if (isTRUE(S$deadline_planted)) return(invisible(NULL))
    S$deadline_planted <- TRUE
    S$deadline_real <- S$G$deadline_at
    S$G$deadline_at <- as.numeric(Sys.time()) - 1
    invisible(NULL)
}

#' A progress dialog, planted in the page (a Shiny-form modal holding a `.progress` element, titled as planted), open
#' 75 s from the press and then removed; the press targets a harmless tab link and the step has no completion signal,
#' with a hard wait of 120 s. The classifier must record the dialog as a progress dialog of about 75 s and end the step
#' `timeout`, never `refused`.
K237_PLANT_MODAL_JS <- r"---((function (show) {
  var w = document.getElementById('k237-planted-modal');
  if (!show) { if (w) { w.remove(); } return true; }
  if (w || document.getElementById('shiny-modal')) { return false; }
  w = document.createElement('div');
  w.id = 'k237-planted-modal';
  w.innerHTML = '<div id="shiny-modal" class="modal fade in show" style="display:block;" tabindex="-1">' +
    '<div class="modal-dialog"><div class="modal-content"><div class="modal-header">' +
    '<h4 class="modal-title">KIT-237 planted progress dialog</h4></div><div class="modal-body">' +
    '<div class="progress"><div class="progress-bar" style="width:50%%"></div></div></div></div></div></div>';
  document.body.appendChild(w);
  return true;
})(%s))---"

k237ModalHook <- function(S, el) {
    if (!isTRUE(S$modal_shown)) {
        S$modal_shown <- isTRUE(k233Js(S, sprintf(K237_PLANT_MODAL_JS, "true")))
        S$modal_at <- el
    } else if (!isTRUE(S$modal_removed) && el - S$modal_at >= 75) {
        S$modal_removed <- isTRUE(k233Js(S, sprintf(K237_PLANT_MODAL_JS, "false")))
        S$modal_removed_at <- el
    }
    invisible(NULL)
}

k237WLprime <- function(G) {
    wf <- G$workflow
    routes <- G$routes
    lfq <- routes$prot_lfq
    guard <- function(criterion, code) tryCatch(code, error = function(e) k237Criterion(G, criterion,
        "driver error while judging", FALSE, conditionMessage(e)))
    # WL'-0: launch and attach (as WL-0), and the outbound block, from the journey step's attempt.
    R1 <- G$runs$prot_lfq("k237_wlp_lfq")
    R1$after <- list(S8 = k233L4a)
    S1 <- k237Launch(G, "WP1")
    k237WL0(G, S1, wf)
    # The outbound row, probe repaired (r1.7, WR.5): the job's GitHub host, port 443, reached by the step's shell and not
    # by the installed R in the same seconds; the DROP attributing it to the rules is recorded beside it.
    k237Criterion(G, "WL-0", "outbound block: the installed application's R cannot reach an address the runner reaches at the same time (the job's GitHub host, port 443)",
        identical(wf$outbound_control_connected, "True") && grepl("cannot be opened|cannot open the connection", wf$outbound_r_said %||% ""),
        sprintf("probe %s; shell connected: %s; R: %s", wf$outbound_probe %||% "-", wf$outbound_control_connected %||% "-",
            wf$outbound_r_said %||% "-"))
    k237Criterion(G, "WL-0", "the block attributed: a firewall DROP logged for R's attempt", NA,
        sprintf("DROP for R's attempt: %s (first found %s s after the attempts); %s", wf$outbound_drop_r %||% "-",
            wf$outbound_log_polled_seconds %||% "-", wf$outbound_method %||% ""))
    if (is.null(S1$driver)) {
        k237StopSession(S1)
        return(k237Verdict(G, "instrument not proven (WL'): launch or attach failed"))
    }
    # WL-1 (to S8) and WL-4 on LFQ, WL-5 on each route (test mode, S1 to S8).
    k237SetOmic(lfq$omic, lfq$tabs_id)
    R1$live <- new.env(parent = emptyenv())
    S1$R <- R1
    S1$plant_delete <- list(step = "S6.2", pattern = "(^|/)post_norm_pca\\.png$")
    k233RunRoute(S1, lfq, R1, upto = "S8")
    k237StopSession(S1)
    if (!any(vapply(G$criteria, function(c) identical(c$criterion, "L-4") && startsWith(c$item, "(a)"), logical(1)))) {
        k237Criterion(G, "L-4", "(a) the emptied DA results table gives a kind-3 finding naming it", FALSE,
            "S8 of WP1 did not complete, so the planted check did not run")
    }
    guard("WL-4", k233L4b(G, S1))
    planted <- Filter(function(f) identical(f$session, "WP1") && identical(f$step, "S6.2") && identical(f$kind, 3L) &&
        isTRUE(f$planted) && grepl("post_norm_pca", f$message, fixed = TRUE), G$findings)
    k237Criterion(G, "WL-4", "(c) a named normalisation figure deleted after S6 gives a kind-3 finding naming it",
        length(planted) == 1L && nzchar(S1$planted_deleted %||% ""), sprintf("deleted: %s", k233Clean(S1$planted_deleted)))
    guard("WL-5", k237L5(G, "WP1"))
    # WL-2, WL-3.
    missing <- paste0(k237Slash(G$work), "/no-such-folder")
    R2 <- G$runs$prot_lfq("k237_wlp2")
    R2$project_dir_override <- missing
    S2 <- k237RunSession(G, lfq, R2, "WP2", upto = "S2")
    k237StopSession(S2)
    s1 <- Filter(function(r) identical(r$session, "WP2") && identical(r$step, "S1"), G$steps)
    k237Criterion(G, "WL-2", "S1 refused with 'Base directory does not exist: <path>' within 60 s",
        length(s1) == 1L && identical(s1[[1L]]$outcome, "refused") &&
            grepl(paste("Base directory does not exist:", missing), s1[[1L]]$message, fixed = TRUE) &&
            as.numeric(s1[[1L]]$seconds) - 60 <= 60 &&
            grepl("-terminal.png at ", k233Clean(s1[[1L]]$terminal_screenshot), fixed = TRUE),
        sprintf("S1 %s after %s s: %s; terminal screenshot: %s", if (length(s1)) s1[[1L]]$outcome else "-",
            if (length(s1)) s1[[1L]]$seconds else "-", if (length(s1)) s1[[1L]]$message else "-",
            if (length(s1)) k233Clean(s1[[1L]]$terminal_screenshot) else "-"))
    S3 <- k237RunSession(G, lfq, G$runs$prot_lfq("k237_wlp3"), "WP3", upto = "S2", hooks = list(S2 = k237KillTreeHook))
    k237StopSession(S3)
    s2 <- Filter(function(r) identical(r$session, "WP3") && identical(r$step, "S2"), G$steps)
    msg <- if (length(s2)) k233Clean(s2[[1L]]$message) else ""
    k237Criterion(G, "WL-3", "the tree ended while S2 waits: the port refuses within 10 s and S2 is crashed from the browser session",
        length(s2) == 1L && identical(s2[[1L]]$outcome, "crashed") && !is.na(S3$port_refused_s %||% NA) &&
            grepl("browser", msg, fixed = TRUE), sprintf("S2 %s: %s; port refused after %s s",
            if (length(s2)) s2[[1L]]$outcome else "-", msg, k233Clean(S3$port_refused_s)))
    # WL'-D: the internal deadline, planted at S2.
    Sd <- k237RunSession(G, lfq, G$runs$prot_lfq("k237_wlpd"), "WPd", upto = "S3", hooks = list(S2 = k237DeadlineHook))
    if (!is.null(Sd$deadline_real)) G$deadline_at <- Sd$deadline_real
    k237StopSession(Sd)
    k233WriteAll(G)
    rows <- Filter(function(r) identical(r$session, "WPd"), G$steps)
    s2d <- Filter(function(r) identical(r$step, "S2"), rows)
    later <- Filter(function(r) !r$step %in% c("S1", "S2"), rows)
    five <- Filter(function(f) identical(f$session, "WPd") && identical(f$kind, 5L), G$findings)
    k237Criterion(G, "WL-D", "the internal deadline (planted at S2): S2 deadline, later steps not_reached (deadline), a kind-5 finding",
        # W08: the tables exist from WL-0 on, so their clause is not a test; the job-level end at a deadline is not planted.
        length(s2d) == 1L && identical(s2d[[1L]]$outcome, "deadline") && length(later) > 0 &&
            all(vapply(later, function(r) identical(r$outcome, "not_reached (deadline)"), logical(1))) && length(five) >= 1L,
        sprintf("S2 %s; later: %s; kind-5 findings %d; deadline restored: %s", if (length(s2d)) s2d[[1L]]$outcome else "-",
            paste(vapply(later, function(r) r$outcome, character(1)), collapse = ","), length(five),
            !is.null(Sd$deadline_real)))
    # WL'-M: a progress dialog open over 60 s (planted).
    Sm <- k237RunSession(G, lfq, G$runs$prot_lfq("k237_wlpm"), "WPm", upto = "S1")
    if (!is.null(Sm$driver)) {
        k233Action(Sm, list(step = "PM", name = "a planted progress dialog held open 75 s", by = "css",
            target = "#proteomics_workflow-workflow_tabs a[data-value='setup']", tabs = list(),
            tab_note = "planted (no tab change)", hard_s = 120, hook = k237ModalHook, dom = function(st) FALSE,
            planted = TRUE))
    }
    k237StopSession(Sm)
    pm <- Filter(function(r) identical(r$session, "WPm") && identical(r$step, "PM"), G$steps)
    k237Criterion(G, "WL-M", "a progress dialog open over 60 s is recorded as progress and does not end the step (timeout at 120 s, not refused)",
        length(pm) == 1L && identical(pm[[1L]]$outcome, "timeout") && grepl("KIT-237 planted progress dialog",
            k233Clean(pm[[1L]]$progress_modals), fixed = TRUE) && !nzchar(k233Clean(pm[[1L]]$other_modals)) &&
            isTRUE(Sm$modal_removed_at >= 75),
        sprintf("PM %s after %s s; progress dialogs: %s; other dialogs: %s; shown at %s s, removed at %s s",
            if (length(pm)) pm[[1L]]$outcome else "-", if (length(pm)) pm[[1L]]$seconds else "-",
            if (length(pm)) k233Clean(pm[[1L]]$progress_modals) else "-", if (length(pm)) k233Clean(pm[[1L]]$other_modals) else "-",
            k233Clean(Sm$modal_at), k233Clean(Sm$modal_removed_at)))
    # WL-5 on the three further routes (S1 to S8), and the chooser once in user mode (LFQ).
    for (id in setdiff(names(routes), "prot_lfq")) {
        if (k233Remaining(G) < 300) {
            k237Criterion(G, "WL-5", sprintf("%s: S1 to S8 in test mode", id), FALSE, "not started: less than 5 minutes left")
            next
        }
        label <- paste0("WP5-", id)
        S <- k237RunSession(G, routes[[id]], G$runs[[id]](paste0("k237_wlp_", id)), label, upto = "S8")
        k237StopSession(S)
        guard("WL-5", k237L5(G, label))
    }
    # W06 (b): the chooser in user mode on each of the four routes (as WL's WL-0u); (d): on LFQ the user-mode session
    # goes on to S8 (the only lifecycle use of user mode beyond S2; its later steps are recorded, not judged). DC-1: the
    # three judged chooser sessions run first, so a stall later in the unjudged LFQ session cannot fail them.
    for (id in c(setdiff(names(routes), "prot_lfq"), "prot_lfq")) {
        label <- paste0("WPu-", id)
        Su <- k237RunSession(G, routes[[id]], G$runs[[id]](paste0("k237_wlpu_", id)), label, test_mode = FALSE,
            upto = if (identical(id, "prot_lfq")) "S8" else "S2")
        k237StopSession(Su)
        s2u <- Filter(function(r) identical(r$session, label) && identical(r$step, "S2"), G$steps)
        wanted <- if (identical(id, "prot_lfq")) 2L else sum(vapply(Su$R$import %||% list(),
            function(a) identical(a$kind, "upload"), logical(1)))
        k237Criterion(G, "WL-0u", sprintf("the chooser selects every input file in user mode on %s", id),
            length(s2u) == 1L && identical(s2u[[1L]]$outcome, "complete") && length(Su$chooser) == wanted &&
                all(Su$chooser),
            sprintf("S2 %s; chooser selections %d of %d", if (length(s2u)) s2u[[1L]]$outcome else "-",
                sum(Su$chooser %in% TRUE), wanted))
    }
    k237Verdict(G)
}

k237WL1Criteria <- function(G, label, criterion = "WL-1", proteomics = TRUE) {
    crit <- function(item, pass, evidence = "") k237Criterion(G, criterion, item, pass, evidence)
    rows <- Filter(function(r) identical(r$session, label), G$steps)
    S <- G$sessions[[label]]
    outcome <- vapply(rows, function(r) k233Clean(r$outcome), character(1))
    step <- vapply(rows, function(r) k233Clean(r$step), character(1))
    crit(sprintf("%s: every step's row has an outcome and seconds", label), length(rows) > 0 && all(nzchar(outcome)) &&
        all(vapply(rows, function(r) !is.null(r$seconds) && !is.na(r$seconds), logical(1))),
        paste(step, outcome, collapse = "; "))
    last <- S$R$last_step
    crit(sprintf("%s: reached the last action (%s)", label, last), any(step == last & outcome == "complete"),
        paste(utils::tail(paste(step, outcome), 3L), collapse = "; "))
    workers <- vapply(rows, function(r) k233Clean(r$workers), character(1))
    crit(sprintf("%s: at least one heavy step ran in a worker process seen in the run folder's logs", label),
        any(nzchar(workers) & workers != "controller"), paste(step[workers != "controller"], workers[workers !=
            "controller"], collapse = "; "))
    complete <- Filter(function(r) identical(r$outcome, "complete") && !identical(r$step, "S10.4") &&
        !identical(r$step, "S9.d"), rows)
    listed <- unlist(lapply(complete, function(r) sub(" .*$", "", strsplit(k233Clean(r$files), " \\| ")[[1L]])))
    listed <- listed[nzchar(listed)]
    present <- file.exists(file.path(G$out, listed)) & file.size(file.path(G$out, listed)) > 0
    crit(sprintf("%s: every screenshot and JSON exists and is non-empty", label), length(listed) > 0 && all(present),
        sprintf("%d files, %d present", length(listed), sum(present)))
    named <- Filter(function(r) length(S$R$files[[r$step]]) > 0, complete)
    listed_named <- vapply(named, function(r) length(list.files(G$out, pattern = sprintf("^%s-[0-9]+-%s-files\\.tsv$",
        label, gsub(".", "\\.", r$step, fixed = TRUE)))) > 0, logical(1))
    crit(sprintf("%s: every named file listed with size and SHA-256, or a kind-3 finding", label),
        length(named) > 0 && all(listed_named), paste(vapply(named, function(r) paste(r$step, k233Clean(r$named_files)),
            character(1)), collapse = "; "))
    numbers <- Filter(function(n) identical(n$session, label), G$numbers)
    ids <- vapply(numbers, function(n) k233Clean(n$id), character(1))
    computed <- vapply(numbers, function(n) !identical(n$result, "not_computed"), logical(1))
    wanted <- if (proteomics) c("N1", "N2", "N3") else "N3"
    crit(sprintf("%s: %s computed", label, paste(wanted, collapse = ", ")), all(wanted %in% ids[computed]),
        paste(ids, vapply(numbers, function(n) k233Clean(n$result), character(1)), collapse = "; "))
    if (proteomics) {
        d <- S$enrich_download
        crit(sprintf("%s: the S9 download recorded (bytes; opens as a zip or not)", label), !is.null(d),
            if (is.null(d)) "no S9 download recorded" else sprintf("%s bytes; opens as a zip: %s; %s", d$bytes, d$opens,
                d$note))
    }
    pressed <- Filter(function(r) !startsWith(k233Clean(r$outcome), "not_reached"), rows)
    checks <- vapply(pressed, function(r) k233Clean(r$tab_checks), character(1))
    crit(sprintf("%s: the active-tab check recorded for every press", label), all(nzchar(checks)) &&
        !any(grepl("NOT ", checks, fixed = TRUE)), paste(utils::head(checks[grepl("NOT ", checks, fixed = TRUE)], 3L),
            collapse = " || "))
    terminal <- vapply(rows, function(r) identical(r$outcome, "refused") && startsWith(k233Clean(r$message), "modal '"),
        logical(1))
    crit(sprintf("%s: no progress modal classified terminal", label), !any(terminal), "")
    crit(sprintf("%s: the R log is present", label), file.exists(S$log) && file.size(S$log) > 0, basename(S$log))
}

k237L5 <- function(G, label) {
    rows <- Filter(function(r) identical(r$session, label) && identical(r$outcome, "complete"), G$steps)
    agreement <- vapply(rows, function(r) k233Clean(r$agreement), character(1))
    bad <- vapply(rows, function(r) r$step, character(1))[!agreement %in% c("agree", "dom_only")]
    k237Criterion(G, "WL-5", sprintf("%s: the digest and the DOM signals agree on completion at every step", label),
        !length(bad), if (length(bad)) paste("disagree at", paste(bad, collapse = ", ")) else
            sprintf("%d agree, %d DOM only", sum(agreement == "agree"), sum(agreement == "dom_only")))
}

k237WL <- function(G) {
    wf <- G$workflow
    routes <- G$routes
    lfq <- routes$prot_lfq
    # WL-0 and WL-1 on the label-free route, with the plants of WL-4 (a) and (c).
    R1 <- G$runs$prot_lfq("k237_wl_lfq")
    R1$after <- list(S8 = k233L4a)
    S1 <- k237Launch(G, "WL1")
    k237WL0(G, S1, wf)
    if (is.null(S1$driver)) {
        k237StopSession(S1)
        return(k237Verdict(G, "instrument not proven on Windows: WL-0 (launch or attach) failed"))
    }
    k237SetOmic(lfq$omic, lfq$tabs_id)
    R1$live <- new.env(parent = emptyenv())
    S1$R <- R1
    S1$plant_delete <- list(step = "S6.2", pattern = "(^|/)post_norm_pca\\.png$")
    k233RunRoute(S1, lfq, R1)
    k237StopSession(S1)
    k237Criterion(G, "WL-0", "the project root and its length recorded (S1 kept the rendered default)",
        !is.null(R1$live$project_dir), sprintf("%s (%s characters); base %s (%s characters)", R1$live$project_dir,
            nchar(R1$live$project_dir %||% ""), k237Base(S1), nchar(k237Base(S1) %||% "")))
    k237Criterion(G, "WL-0", "the R log, worker logs and MSR-* lines read", NA, sprintf("%d MSR lines; worker logs: %s",
        sum(grepl("^MSR-", k233ReadLog(S1$log))), paste(grep(K237_WORKER_LOG, k237Slash(k237Listing(S1$run_folder)$path),
            value = TRUE), collapse = " ")))
    guard <- function(criterion, code) tryCatch(code, error = function(e) k237Criterion(G, criterion,
        "driver error while judging", FALSE, conditionMessage(e)))
    guard("WL-1", k237WL1Criteria(G, "WL1"))
    if (!any(vapply(G$criteria, function(c) identical(c$criterion, "L-4") && startsWith(c$item, "(a)"), logical(1)))) {
        k237Criterion(G, "L-4", "(a) the emptied DA results table gives a kind-3 finding naming it", FALSE,
            "S8 of WL-1 did not complete, so the planted check did not run")
    }
    guard("WL-4", k233L4b(G, S1))
    planted <- Filter(function(f) identical(f$session, "WL1") && identical(f$step, "S6.2") && identical(f$kind, 3L) &&
        isTRUE(f$planted) && grepl("post_norm_pca", f$message, fixed = TRUE), G$findings)
    k237Criterion(G, "WL-4", "(c) a named normalisation figure deleted after S6 gives a kind-3 finding naming it",
        length(planted) == 1L && nzchar(S1$planted_deleted %||% ""), sprintf("deleted: %s; findings: %s", k233Clean(S1$planted_deleted),
            paste(vapply(planted, function(f) f$message, character(1)), collapse = " || ")))
    guard("WL-5", k237L5(G, "WL1"))
    # SC-07: the blocking WL-2 and WL-3, and the user-mode chooser sessions, before the three added journeys.
    # WL-2: a base directory that does not exist.
    missing <- paste0(k237Slash(G$work), "/no-such-folder")
    R2 <- G$runs$prot_lfq("k237_wl2")
    R2$project_dir_override <- missing
    S2 <- k237RunSession(G, lfq, R2, "WL2")
    k237StopSession(S2)
    rows <- Filter(function(r) identical(r$session, "WL2"), G$steps)
    s1 <- Filter(function(r) identical(r$step, "S1"), rows)
    later <- Filter(function(r) !identical(r$step, "S1"), rows)
    k237Criterion(G, "WL-2", "S1 refused with 'Base directory does not exist: <path>' within 60 s; later steps not_reached",
        length(s1) == 1L && identical(s1[[1L]]$outcome, "refused") &&
            grepl(paste("Base directory does not exist:", missing), s1[[1L]]$message, fixed = TRUE) &&
            as.numeric(s1[[1L]]$seconds) - 60 <= 60 && length(later) == length(lfq$steps) - 1L &&
            all(vapply(later, function(r) identical(r$outcome, "not_reached"), logical(1))),
        sprintf("S1 %s after %s s (the classifier holds a refusal 60 s): %s; later: %s",
            if (length(s1)) s1[[1L]]$outcome else "-", if (length(s1)) s1[[1L]]$seconds else "-",
            if (length(s1)) s1[[1L]]$message else "-", paste(vapply(later, function(r) r$outcome, character(1)),
                collapse = ",")))
    # WL-3: the whole tree ended while S2 waits.
    S3 <- k237RunSession(G, lfq, G$runs$prot_lfq("k237_wl3"), "WL3", upto = "S2", hooks = list(S2 = k237KillTreeHook))
    k237StopSession(S3)
    s2 <- Filter(function(r) identical(r$session, "WL3") && identical(r$step, "S2"), G$steps)
    msg <- if (length(s2)) k233Clean(s2[[1L]]$message) else ""
    k237Criterion(G, "WL-3", "the tree ended while S2 waits: the port refuses within 10 s and S2 is crashed from the browser session",
        length(s2) == 1L && identical(s2[[1L]]$outcome, "crashed") && !is.na(S3$port_refused_s %||% NA) &&
            grepl("browser", msg, fixed = TRUE),
        sprintf("S2 %s: %s; killed %s s after the press; %s; port refused after %s s", if (length(s2)) s2[[1L]]$outcome
            else "-", msg, k233Clean(S3$kill_at), k233Clean(S3$kill_said), k233Clean(S3$port_refused_s)))
    # SC-02: user mode (VAL-216): the chooser selects every input file, on each of the four routes, to S2 (blocking).
    for (id in names(routes)) {
        label <- paste0("WLu-", id)
        Su <- k237RunSession(G, routes[[id]], G$runs[[id]](paste0("k237_wlu_", id)), label, test_mode = FALSE,
            upto = "S2")
        k237StopSession(Su)
        s2u <- Filter(function(r) identical(r$session, label) && identical(r$step, "S2"), G$steps)
        wanted <- if (identical(id, "prot_lfq")) 2L else sum(vapply(Su$R$import %||% list(),
            function(a) identical(a$kind, "upload"), logical(1)))
        k237Criterion(G, "WL-0u", sprintf("the chooser selects every input file in user mode on %s", id),
            length(s2u) == 1L && identical(s2u[[1L]]$outcome, "complete") && length(Su$chooser) == wanted &&
                all(Su$chooser),
            sprintf("S2 %s: %s; chooser selections %s of %d; %s", if (length(s2u)) s2u[[1L]]$outcome else "-",
                if (length(s2u)) k233Clean(s2u[[1L]]$message) else "-", sum(Su$chooser %in% TRUE),
                wanted, paste(vapply(Filter(function(o) identical(o$session, label), G$observations),
                    function(o) paste(o$what, o$text), character(1)), collapse = " || ")))
    }
    k237Verdict(G)
    # One journey on each further Windows route (OD-1 (b), owner 2026-11-25, VAL-225); WL-4 (a) planted once on metab_lc.
    for (id in setdiff(names(routes), "prot_lfq")) {
        if (k233Remaining(G) < 600) {
            k237Criterion(G, "WL-1", sprintf("%s: journey run", id), FALSE, "not started: less than 10 minutes left")
            next
        }
        label <- paste0("WL1-", id)
        R <- G$runs[[id]](paste0("k237_wl_", id))
        if (identical(id, "metab_lc")) R$after <- list(S8 = function(S, R) k237L4aOmic(S, R, routes[[id]]))
        S <- k237RunSession(G, routes[[id]], R, label)
        k237StopSession(S)
        guard("WL-1", k237WL1Criteria(G, label, proteomics = FALSE))
        guard("WL-5", k237L5(G, label))
        if (identical(id, "metab_lc")) {
            numbers <- Filter(function(n) identical(n$session, label) && identical(n$id, "N2") && !isTRUE(n$planted),
                G$numbers)
            k237Criterion(G, "WL-1", "WL1-metab_lc: N2 computed", length(numbers) > 0 &&
                all(vapply(numbers, function(n) !identical(n$result, "not_computed"), logical(1))),
                paste(vapply(numbers, function(n) paste(n$contrast, n$shown, n$saved, n$result), character(1)),
                    collapse = "; "))
            if (!any(vapply(G$criteria, function(c) identical(c$item, "(a) on metab_lc: the emptied DA results table gives a kind-3 finding naming it"),
                logical(1)))) {
                k237Criterion(G, "WL-4", "(a) on metab_lc: the emptied DA results table gives a kind-3 finding naming it",
                    FALSE, "S8 of the metab_lc journey did not complete, so the planted check did not run")
            }
        }
    }
    k237Verdict(G)
    # Observation, while time remains: a worker ended during the longest QC step.
    qc <- Filter(function(r) identical(r$session, "WL1") && grepl("^S5\\.[0-9]+$", r$step) &&
        identical(r$outcome, "complete"), G$steps)
    if (length(qc) && k233Remaining(G) > 420) {
        longest <- qc[[which.max(vapply(qc, function(r) as.numeric(r$seconds), numeric(1)))]]
        Sw <- k237RunSession(G, lfq, G$runs$prot_lfq("k237_wl3w"), "WL3w", upto = "S5",
            hooks = stats::setNames(list(k233KillWorkerHook), longest$step))
        k237StopSession(Sw)
        shown <- Filter(function(r) identical(r$session, "WL3w") && identical(r$step, longest$step), G$steps)
        k233Observation(G, "WL3w", sprintf("worker ended during %s (%s, the longest QC step of WL-1)", longest$step,
            longest$name), sprintf("worker ended: %s at %s s; the step: %s", isTRUE(Sw$worker_killed),
                k233Clean(Sw$worker_kill_at), if (length(shown)) paste(shown[[1L]]$outcome, shown[[1L]]$message) else "-"))
    }
    k237Verdict(G)
}

k237Verdict <- function(G, stop_text = NULL) {
    for (i in seq_along(G$criteria)) G$criteria[[i]]$criterion <- sub("^L-", "WL-", G$criteria[[i]]$criterion)
    results <- vapply(G$criteria, function(c) c$result, character(1))
    names <- vapply(G$criteria, function(c) c$criterion, character(1))
    blocking <- names %in% c("WL-0", "WL-0u", "WL-1", "WL-2", "WL-3", "WL-4", "WL-D", "WL-M")
    proven <- is.null(stop_text) && !any(results[blocking] == "fail")
    # WL-0u is required on all four routes (SC-02): a missing route's row counts as not proven.
    routes_u <- sum(names == "WL-0u")
    if (!is.null(G$routes) && routes_u < length(G$routes)) proven <- FALSE
    l5 <- !any(results[names == "WL-5"] == "fail")
    writeLines(c(stop_text %||% sprintf("instrument %s on Windows (WL-0, WL-0u, WL-1 to WL-4: %d pass, %d fail, %d recorded)",
        if (proven) "proven" else "not proven", sum(results[blocking] == "pass"), sum(results[blocking] == "fail"),
        sum(results[blocking] == "recorded")),
        sprintf("WL-5 %s", if (l5) "passed" else "not passed: no user-mode WJ' may rely on the DOM signals")),
        file.path(G$out, "verdict.txt"))
    k233WriteAll(G)
}

# ---------------------------------------------------------------------------------------------------------------
# Main.
# ---------------------------------------------------------------------------------------------------------------
k237ReadWorkflowRecord <- function(out) {
    path <- file.path(out, "environment-workflow.tsv")
    if (!file.exists(path)) return(list())
    t <- utils::read.delim(path, stringsAsFactors = FALSE, colClasses = "character")
    stats::setNames(as.list(t$value), t$key)
}

k237Main <- function(args) {
    if (length(args) != 7L) stop("usage: journey-win.R <request> <application folder> <test tools> <out> <source> ",
        "<fixtures> <checkout>")
    a <- stats::setNames(as.list(args), c("request", "app_dir", "test_tools", "out", "src", "fixtures", "pkg"))
    started <- Sys.time()
    .libPaths(c(a$test_tools, .libPaths()))
    dir.create(a$out, recursive = TRUE, showWarnings = FALSE)
    out <- normalizePath(a$out, winslash = "/")
    request <- jsonlite::read_json(a$request, simplifyVector = FALSE)
    source(file.path(a$src, "docs/qa/workflow-execution/kit/kit-233/journey.R"))
    for (f in request$route_files) source(file.path(a$src, K237_RECORDS, f))
    rebound <- k237Rebind()
    # Inputs: every file each route reads, against the request's SHA-256 (the manifest's), before anything starts.
    inputs <- list()
    for (id in names(request$inputs)) {
        spec <- request$inputs[[id]]
        for (f in names(spec$files)) {
            path <- file.path(a$fixtures, spec$dir, f)
            got <- k233Sha(path)
            inputs[[length(inputs) + 1L]] <- list(route = id, file = paste(spec$dir, f, sep = "/"),
                sha256_request = spec$files[[f]], sha256 = got, equal = identical(got, spec$files[[f]]))
        }
    }
    k233WriteTsv(inputs, file.path(out, "inputs.tsv"))
    if (!length(inputs) || !all(vapply(inputs, function(r) isTRUE(r$equal), logical(1)))) {
        writeLines("stopped: an input's SHA-256 differs from the request (inputs.tsv)", file.path(out, "stop.txt"))
        return(invisible(k233Finish(out)))
    }
    work <- file.path(tempdir(), "k237-work")
    dir.create(work, recursive = TRUE, showWarnings = FALSE)
    phase <- list(run_mode = request$run, started = as.numeric(started),
        bounds = list(deadline_min = as.numeric(request$deadline_minutes), hard_s = 600, finding_s = 600, poll = 0.5),
        pkg = normalizePath(a$pkg, winslash = "/"), out = out, work = normalizePath(work, winslash = "/"),
        commit = request$commit_under_test, records_commit = request$route_commit,
        inputs = list(id = paste(names(request$inputs), collapse = "+"), sha = ""),
        launch = k233LaunchSpec(form = "installed", app_dir = normalizePath(a$app_dir, winslash = "/"),
            network = "firewall"),
        user_outer = Sys.info()[["effective_user"]])
    Sys.setenv(NOT_CRAN = "true", TZ = "UTC", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    Sys.unsetenv(c("MULTISCHOLAR_TEST_FEATURES", "MULTISCHOLAR_TEST_MODE", "CHROMOTE_CHROME"))
    options(warn = 1, chromote.timeout = 60)
    G <- k233Context(phase)
    assign("K237_G", G, envir = globalenv())
    G$workflow <- k237ReadWorkflowRecord(out)
    G$test_features <- isTRUE(request$test_features)
    G$routes <- k237Routes(request$routes)
    G$runs <- k237Runs(G, a$fixtures, request)
    file.copy(file.path(a$src, K237_RECORDS, "known.tsv"), file.path(out, "known.tsv"), overwrite = TRUE)
    on.exit({
        for (S in G$sessions) k237StopSession(S)
        k233WriteAll(G)
        k237WriteEnvironment(G, request, rebound)
        k233Finish(out)
    }, add = TRUE)
    fatal <- function(e) {
        G$fatal <- paste(conditionMessage(e), "| call:", paste(deparse(conditionCall(e)), collapse = " "))
        writeLines(G$fatal, file.path(out, "inner-fatal.txt"))
    }
    tryCatch(if (identical(request$run, "WLprime")) k237WLprime(G) else k237WL(G), error = fatal)
    invisible(G)
}

#' The routes the request names, from the route files it sourced (`K237_ROUTES`, one entry per route id).
k237Routes <- function(ids) {
    routes <- K237_ROUTES[unlist(ids)]
    if (any(vapply(routes, is.null, logical(1)))) stop("a requested route has no table: ", paste(unlist(ids), collapse = ", "))
    routes
}

#' Per route, a function of the label that gives the route's run list (`K237_RUNS`), with its inputs' paths.
k237Runs <- function(G, fixtures, request) {
    lapply(stats::setNames(names(G$routes), names(G$routes)), function(id) {
        function(label) K237_RUNS[[id]](G, file.path(fixtures, request$inputs[[id]]$dir), label)
    })
}

k237WriteEnvironment <- function(G, request, rebound) {
    wf <- G$workflow
    S <- G$sessions[["WL1"]] %||% G$sessions[["WP1"]]
    R <- if (is.null(S)) NULL else S$R
    chrome <- tryCatch(as.character(chromote::find_chrome()), error = function(e) "")
    chrome_version <- if (nzchar(chrome)) tryCatch(paste(suppressWarnings(system2(chrome, "--version", stdout = TRUE,
        stderr = TRUE)), collapse = " "), error = function(e) "?") else "not found"
    rows <- c(lapply(names(wf), function(k) list(key = k, value = wf[[k]])), list(
        list(key = "run", value = request$run),
        list(key = "commit_under_test", value = request$commit_under_test),
        list(key = "journey_commit (KIT-233 journey.R)", value = request$journey_commit),
        list(key = "route_commit (KIT-237 files)", value = request$route_commit),
        list(key = "routes", value = paste(unlist(request$routes), collapse = " ")),
        list(key = "mode", value = request$mode),
        list(key = "network", value = "firewall: outbound block rules for every .exe under the installed application"),
        list(key = "r_version", value = R.version.string),
        list(key = "r_home", value = R.home()),
        list(key = "chrome", value = chrome),
        list(key = "chrome_version", value = chrome_version),
        list(key = "shinytest2", value = tryCatch(as.character(utils::packageVersion("shinytest2")), error = function(e) "-")),
        list(key = "cpu", value = Sys.getenv("PROCESSOR_IDENTIFIER")),
        list(key = "cores", value = parallel::detectCores()),
        list(key = "launcher_line", value = if (is.null(S)) "" else S$launcher),
        list(key = "LongPathsEnabled_read_by_R", value = k233Clean(G$long_paths)),
        list(key = "project_root", value = if (is.null(R)) "" else k233Clean(R$live$project_dir)),
        list(key = "project_root_length", value = if (is.null(R)) "" else nchar(R$live$project_dir %||% "")),
        list(key = "application_libpaths", value = k233Clean(G$app_libpaths)),
        list(key = "journey_R_rebound", value = paste(rebound, collapse = "; ")),
        list(key = "test_features", value = isTRUE(G$test_features)),
        # SC-13: how the application's start here differs from the shell's (packaging/electron/main.js:125-175).
        list(key = "launch_differences_from_shell", value = paste(
            "a wrapper file sources the installed entry.R (test mode: the test-mode options first; always one probe line",
            "writing .libPaths() into the run folder); working directory = the checkout; NOT_CRAN=true, TZ=UTC and",
            "OMP/OPENBLAS/MKL_NUM_THREADS=1 inherited from the driver; MULTISCHOLAR_TEST_MODE=true in test mode;",
            "CHROMOTE_CHROME points at nothing; the shell's own variables otherwise as windows-clean-run.ps1 sets them;",
            "MSR_LAUNCHER_PROBE given as the shell gives it on Windows (process, 'win32: no systemd');",
            "the shell's pid replaced by the driver's")),
        list(key = "github_run_id", value = Sys.getenv("GITHUB_RUN_ID")),
        list(key = "finished_utc", value = k233Utc())))
    k233WriteTsv(rows, file.path(G$out, "environment.tsv"))
}

if (sys.nframe() == 0L) k237Main(commandArgs(trailingOnly = TRUE))
