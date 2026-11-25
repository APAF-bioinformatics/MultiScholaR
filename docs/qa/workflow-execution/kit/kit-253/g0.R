# KIT-253 G0 (r1.2): does a Reset or a Save Design pressed while a normalisation's figures action runs on the managed
# worker publish anything? One scripts job on a GitHub-hosted runner, bound to v0.9-staging 16112917. Records only:
# nothing under R/ is changed; the application is driven in a headless browser through the plain launcher
# (packaging/electron/entry.R, process launcher, managed worker), in test mode, with the network closed.
#
#   Rscript docs/qa/workflow-execution/kit/kit-253/g0.R outer <commit under test>
#   Rscript docs/qa/workflow-execution/kit/kit-253/g0.R inner <phase file>      (the browser phase, in the namespace)
#
# KIT-233's driver (journey.R on origin/kit-233-journey, read with git show; sourcing it runs nothing, its N1 to N3)
# launches the application, drives the proteomics route to protein QC (its S1 to S5), presses and records the page.
# The press is placed by the project store's request folders: `<store>/.workflow-export/<request id>/request.rds`
# (its `operation`) is written at dispatch, `completion.rds` (its `status`) at the job's end
# (R/utils_workflow_execution.R:1869, :2262). Each press is made after the target request appears and before its
# completion.rds exists; both times are recorded.

K253_OUT <- "k253-out"
K253_RECORDS <- "docs/qa/workflow-execution/kit/kit-253"
K253_JOURNEY <- c(ref = "origin/kit-233-journey", path = "docs/qa/workflow-execution/kit/kit-233/journey.R")
K253_POLL <- 0.2
K253_PROT_LABEL <- "k253_prot"
K253_LIPID_LABEL <- "k253_lipid"
K253_DONE <- c(proteomics = "Normali[sz]ation (and RUV correction )?completed", lipidomics = "Normalization pipeline complete!")

# ---------------------------------------------------------------------------------------------------------------
# The prediction (KIT-253 r1.2 section 9), fixed here before any value is read. P0 and L0 are the positive controls:
# the same journey with no press, which must show the target completing and a publication the instrument can see.
# ---------------------------------------------------------------------------------------------------------------
K253_SEQUENCES <- list(
    list(seq = "P0", omic = "proteomics", press = "none (positive control)", target = "normalization-qc-figures",
        after = "normalization-ruv", expect_status = "complete", expect_publication = TRUE),
    list(seq = "P-a", omic = "proteomics", press = "Reset to Pre-Normalization", target = "normalization-qc-figures",
        after = "normalization-ruv", expect_status = "stale", expect_publication = FALSE),
    list(seq = "P-c", omic = "proteomics", press = "Save Design", target = "normalization-qc-figures",
        after = "normalization-ruv", expect_status = "stale", expect_publication = FALSE),
    list(seq = "L0", omic = "lipidomics", press = "none (positive control)", target = "lipid_composite_figure",
        after = NA_character_, expect_status = "complete", expect_publication = TRUE),
    list(seq = "L-a", omic = "lipidomics", press = "Reset to Pre-Normalization", target = "lipid_composite_figure",
        after = NA_character_, expect_status = "stale", expect_publication = FALSE)
)
K253_NOTICES <- list(
    proteomics = c("Error applying proteomics step: managed action ended: stale",
        "Error in normalization: managed action ended: stale"),
    lipidomics = c("Error applying lipidomics step: managed action ended: stale", "Normalization pipeline complete!")
)
K253_RULES <- c(
    "target: proteomics = the first normalization-qc-figures request whose request.rds appears after the normalization-ruv request's completion.rds reads complete (R1-01); lipidomics = the first lipid_composite_figure request",
    "press inside = the target's request.rds exists and its completion.rds does not, checked immediately before and after the click",
    "as predicted (press sequences): target status stale; no publication request folder outside the store written after the press; no file written after the press at a path in the matching positive control's publication footprint; no completion notice (proteomics); the head moved between the press and the target's completion (the press wrote it); later actions of the omic's chain that completed after the target are recorded (later_chain_complete, later_requests), not judged",
    "positive control proven: target status complete; the completion notice seen; a non-empty publication footprint = the files outside the store written after the target appeared, less those already written by the chain between the run press and the target's appearance",
    "instrument defect (stop rule 3): no target seen, the target completed before the press, the press control not found, or the session did not reach the normalisation",
    "contradicted (stop rule 2): any press sequence not as predicted while its press was inside and its control proven",
    "notices: compared with the predicted texts and recorded; not a falsifier (section 5's falsifiers are the status and the publication)"
)

k253Prediction <- function(out) {
    rows <- lapply(K253_SEQUENCES, \(s) list(seq = s$seq, omic = s$omic, press = s$press, target = s$target,
        expected_status = s$expect_status, expected_publication = s$expect_publication,
        expected_notices = if (s$expect_publication) "completion notice" else paste(K253_NOTICES[[s$omic]], collapse = " || ")))
    utils::write.table(do.call(rbind, lapply(rows, as.data.frame)), file.path(out, "prediction.tsv"), sep = "\t",
        quote = FALSE, row.names = FALSE)
    writeLines(K253_RULES, file.path(out, "prediction-rules.txt"))
}

# ---------------------------------------------------------------------------------------------------------------
# Outer phase.
# ---------------------------------------------------------------------------------------------------------------
k253Main <- function(args) {
    if (length(args) < 2L) stop("usage: g0.R outer <commit under test> | inner <phase file>")
    if (identical(args[[1L]], "inner")) return(invisible(k253Inner(args[[2L]])))
    if (!identical(args[[1L]], "outer")) stop("unknown mode ", args[[1L]])
    k253Outer(args[[2L]])
}

k253Outer <- function(commit) {
    started <- Sys.time()
    pkg <- normalizePath(".")
    out <- file.path(pkg, K253_OUT)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    k253Prediction(out)
    work <- file.path(tempdir(), "k253-work")
    dir.create(work, recursive = TRUE, showWarnings = FALSE)
    git <- function(...) suppressWarnings(system2("git", c("-C", shQuote(pkg), ...), stdout = TRUE, stderr = TRUE))
    head <- git("rev-parse", "HEAD")
    changed <- git("diff", "--name-only", commit, "HEAD")
    diff_ok <- is.null(attr(changed, "status")) && all(startsWith(changed, paste0(K253_RECORDS, "/")))
    journey <- file.path(work, "journey.R")
    blob <- git("rev-parse", paste0(K253_JOURNEY[["ref"]], ":", K253_JOURNEY[["path"]]))
    writeLines(git("show", paste0(K253_JOURNEY[["ref"]], ":", K253_JOURNEY[["path"]])), journey)
    version <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = \(e) "not installed")
    env <- data.frame(key = c("commit_under_test", "records_commit", "diff_vs_commit_under_test", "diff_only_records",
        "journey_ref", "journey_blob", "journey_lines", "r_version", "MultiScholaR", "shinytest2", "chromote",
        "cores", "github_run_id", "started_utc"),
        value = c(commit, head, paste(changed, collapse = " "), diff_ok, K253_JOURNEY[["ref"]], paste(blob, collapse = " "),
            length(readLines(journey)), R.version.string, version("MultiScholaR"), version("shinytest2"),
            version("chromote"), parallel::detectCores(), Sys.getenv("GITHUB_RUN_ID"),
            format(started, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")))
    utils::write.table(env, file.path(out, "environment.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    if (!isTRUE(diff_ok) || length(readLines(journey)) < 100L) {
        writeLines("stopped before the browser phase: the records commit differs from the commit under test outside the records folder, or journey.R was not read",
            file.path(out, "verdict.txt"))
        return(invisible(out))
    }
    script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))
    phase <- list(run_mode = "k253", started = as.numeric(started), pkg = pkg, out = out, work = work,
        journey = journey, commit = commit, user_outer = Sys.info()[["effective_user"]], libpaths = .libPaths(),
        locale = Sys.getlocale("LC_CTYPE"), inputs = list(id = "prot_dia and lipid_canonical fixtures"),
        bounds = list(deadline_min = 68, hard_s = 300, finding_s = 600, poll = 0.5),
        launch = list(form = "tree", app_dir = NULL, network = "namespace", e_form = "e"),
        pandoc = list(dir = "", source = "not used (no report in G0)"),
        env = c(R_LIBS = Sys.getenv("R_LIBS"), R_LIBS_USER = Sys.getenv("R_LIBS_USER"),
            R_LIBS_SITE = Sys.getenv("R_LIBS_SITE"), LANG = Sys.getenv("LANG", "C.UTF-8"),
            GITHUB_RUN_ID = Sys.getenv("GITHUB_RUN_ID"), TZ = "UTC", OMP_NUM_THREADS = "1",
            OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", NOT_CRAN = "true"))
    phase$env <- phase$env[nzchar(phase$env)]
    phase_file <- file.path(work, "k253-phase.rds")
    saveRDS(phase, phase_file)
    console <- file.path(out, "inner-console.log")
    command <- paste('ip link set lo up && exec sudo -u "$0" env PATH="$1" HOME="$2"',
        "http_proxy=http://127.0.0.1:9 https_proxy=http://127.0.0.1:9 no_proxy=127.0.0.1,localhost",
        'Rscript "$3" inner "$4"')
    user <- system2("id", "-un", stdout = TRUE)
    status <- system2("sudo", c("unshare", "--net", "--", "bash", "-c", shQuote(command), shQuote(user),
        shQuote(Sys.getenv("PATH")), shQuote(Sys.getenv("HOME")), shQuote(script), shQuote(phase_file)),
        stdout = console, stderr = console)
    writeLines(sprintf("namespace command exit status %s", status), file.path(out, "namespace-status.txt"))
    for (f in c("verdict.txt", "inner-fatal.txt")) if (file.exists(file.path(out, f))) cat(readLines(file.path(out, f)), sep = "\n")
    if (file.exists(file.path(out, "g0-sequences.tsv"))) cat(readLines(file.path(out, "g0-sequences.tsv")), sep = "\n")
    invisible(out)
}

# ---------------------------------------------------------------------------------------------------------------
# Inner phase.
# ---------------------------------------------------------------------------------------------------------------
k253Inner <- function(phase_file) {
    phase <- readRDS(phase_file)
    .libPaths(phase$libpaths)
    if (length(phase$env)) do.call(Sys.setenv, as.list(phase$env))
    Sys.unsetenv(c("MULTISCHOLAR_TEST_FEATURES", "MULTISCHOLAR_TEST_MODE", "CHROMOTE_CHROME"))
    invisible(tryCatch(Sys.setlocale("LC_ALL", phase$locale), warning = \(w) NULL, error = \(e) NULL))
    options(warn = 1, chromote.timeout = 60)
    source(phase$journey, local = globalenv())
    writeLines(k233Utc(), file.path(phase$out, "inner-started.txt"))
    G <- k233Context(phase)
    fixture_root <- file.path(phase$pkg, "tests/testdata/e2e")
    assign(".e2e_fixture_root", function() fixture_root, envir = G$helpers)
    G$k253 <- list()
    on.exit({
        for (S in G$sessions) k233StopSession(S)
        k233WriteAll(G)
        k253Verdict(G)
    }, add = TRUE)
    fatal <- function(e) writeLines(paste(conditionMessage(e), "| call:", paste(deparse(conditionCall(e)), collapse = " ")),
        file.path(G$out, "inner-fatal.txt"))
    tryCatch({
        for (s in K253_SEQUENCES) {
            if (k233Remaining(G) < 240) {
                G$k253[[s$seq]] <- list(seq = s$seq, instrument = "not run: the driver's deadline was near")
                next
            }
            G$k253[[s$seq]] <- tryCatch(k253Sequence(G, s), error = \(e) list(seq = s$seq,
                instrument = paste("driver error:", conditionMessage(e))))
            k253WriteSequences(G)
        }
    }, error = fatal)
    invisible(G)
}

# --- the store watcher ------------------------------------------------------------------------------------------

k253Watcher <- function(S) {
    W <- new.env(parent = emptyenv())
    W$S <- S
    W$reqs <- list()
    W$heads <- list()
    W$head <- NA_character_
    W$stores <- character()
    W$t0 <- Sys.time()
    W
}

k253FindStores <- function(W) {
    roots <- c(W$S$project_dir, W$S$run_folder)
    hits <- unlist(lapply(roots, \(r) list.files(r, "^store\\.json$", recursive = TRUE, full.names = TRUE,
        all.files = TRUE)))
    W$stores <- unique(dirname(hits))
    W$stores
}

k253Rel <- function(W, t) round(as.numeric(t) - as.numeric(W$t0), 2)

k253Scan <- function(W) {
    now <- Sys.time()
    if (!length(W$stores)) k253FindStores(W)
    for (store in W$stores) {
        root <- file.path(store, ".workflow-export")
        ids <- if (dir.exists(root)) list.files(root) else character()
        for (id in ids) {
            dir <- file.path(root, id)
            r <- W$reqs[[id]]
            if (is.null(r)) {
                path <- file.path(dir, "request.rds")
                if (!file.exists(path)) next
                value <- tryCatch(readRDS(path), error = \(e) NULL)
                if (is.null(value)) next
                r <- list(id = id, dir = dir, operation = as.character(value$operation %||% NA_character_), seen = now,
                    request_mtime = file.mtime(path), expected_head = as.character(value$token$expected_head %||% NA),
                    completion_seen = NA, status = NA_character_, completion_mtime = NA, error = "")
            }
            if (is.na(r$status)) {
                path <- file.path(dir, "completion.rds")
                if (file.exists(path)) {
                    value <- tryCatch(readRDS(path), error = \(e) NULL)
                    if (!is.null(value)) {
                        r$completion_seen <- now
                        r$completion_mtime <- file.mtime(path)
                        r$status <- as.character(value$status %||% "?")
                        r$error <- k233Clean(value$error)
                    }
                }
            }
            W$reqs[[id]] <- r
        }
        head_path <- file.path(store, "dia-state-head.rds")
        if (file.exists(head_path)) {
            h <- tryCatch(readRDS(head_path), error = \(e) NULL)
            id <- if (is.list(h)) as.character(h$entity_id %||% NA_character_) else NA_character_
            if (!is.na(id) && !identical(id, W$head)) {
                W$head <- id
                W$heads[[length(W$heads) + 1L]] <- list(t = k253Rel(W, now), store = basename(store), entity_id = id)
            }
        }
    }
    invisible(W)
}

k253Open <- function(W) Filter(\(r) is.na(r$status), W$reqs)

k253RequestRows <- function(W) lapply(W$reqs, \(r) list(id = r$id, operation = r$operation,
    request_seen_s = k253Rel(W, r$seen), request_mtime_s = k253Rel(W, r$request_mtime),
    completion_seen_s = if (is.na(r$status)) NA else k253Rel(W, r$completion_seen),
    completion_mtime_s = if (is.na(r$status)) NA else k253Rel(W, r$completion_mtime), status = r$status,
    expected_head = r$expected_head, error = r$error))

# --- one sequence ----------------------------------------------------------------------------------------------

k253Sequence <- function(G, s) {
    label <- gsub("[^A-Za-z0-9]", "", s$seq)
    S <- k233Launch(G, label)
    rec <- list(seq = s$seq, omic = s$omic, press = s$press, target_operation = s$target)
    if (!isTRUE(S$answered) || is.null(S$driver)) {
        rec$instrument <- paste("the application did not answer or could not be attached:", S$attach_error)
        return(rec)
    }
    W <- k253Watcher(S)
    reached <- if (identical(s$omic, "proteomics")) k253ProtToNorm(G, S) else k253LipidToNorm(G, S)
    if (!isTRUE(reached$ok)) {
        rec$instrument <- paste("did not reach the normalisation:", reached$why)
        k233StopSession(S)
        return(rec)
    }
    rec$settings <- reached$settings
    k253FindStores(W)
    rec$stores <- paste(W$stores, collapse = " ")
    k253Scan(W)
    known_before <- names(W$reqs)
    press_fn <- k253PressFunction(S, s)
    since_ms <- as.numeric(k233Js(S, "Date.now()") %||% 0)
    run <- k233Press(S, s$seq, if (identical(s$omic, "proteomics")) "prot-norm-run" else "lipid-norm-run", "testid")
    W$t0 <- Sys.time()
    rec$head_at_run <- W$head
    W$heads <- list()
    if (!isTRUE(run$ok)) {
        rec$instrument <- paste("the run button:", run$reason)
        k233StopSession(S)
        return(rec)
    }
    if (is.function(press_fn$prepare)) rec$prepare <- press_fn$prepare()
    window <- k253Window(G, S, W, s, press_fn$press, known_before)
    rec <- c(rec, window)
    settle <- k253Settle(G, S, W, s, window, since_ms)
    rec <- c(rec, settle)
    k253Collect(G, S, W, s, rec, since_ms)
}

k253ProtToNorm <- function(G, S) {
    R <- k233FixtureRun(G, S$project_dir)
    R$label <- K253_PROT_LABEL
    res <- k233RunRoute(S, G$route, R, upto = "S5")
    if (!identical(res$outcome, "complete")) return(list(ok = FALSE, why = paste("route to S5:", res$outcome, res$message %||% "")))
    tab <- k233ActivateTab(S, paste0(K233_NS, "workflow_tabs"), "normalization")
    norm <- function(x) paste0(K233_NS, "normalization-", x)
    settings <- sprintf("defaults: ruv_mode=%s norm_method=%s; %s", k233Clean(k233InputValue(S, norm("ruv_mode"))),
        k233Clean(k233InputValue(S, norm("norm_method"))), tab)
    list(ok = TRUE, settings = settings)
}

k253Status <- function(S, omic, step) {
    st <- k233State(S, 0)
    value <- st$dg$step_status_per_omic[[omic]][[step]]
    if (is.null(value)) NA_character_ else as.character(value)
}

k253WaitStatus <- function(S, omic, step, seconds) {
    started <- Sys.time()
    repeat {
        if (identical(k253Status(S, omic, step), "complete")) return(TRUE)
        if (k233Seconds(started) > seconds || k233Remaining(S$G) <= 0) return(FALSE)
        Sys.sleep(1)
    }
}

k253LipidToNorm <- function(G, S) {
    lane <- jsonlite::read_json(file.path(G$pkg, "tests/testdata/e2e/manifest.json"))
    lanes <- if (!is.null(lane$lanes)) lane$lanes else lane
    lane <- Filter(\(l) identical(l$lane_id, "lipid_canonical"), lanes)
    if (!length(lane)) return(list(ok = FALSE, why = "lipid_canonical lane not in manifest.json"))
    lane <- lane[[1L]]
    H <- G$helpers
    step <- function(name, ok, why = "") {
        G$steps[[length(G$steps) + 1L]] <- list(session = S$label, step = name, name = name,
            outcome = if (isTRUE(ok)) "complete" else "refused", seconds = NA, started_utc = k233Utc(), message = why)
        k233WriteAll(G)
        isTRUE(ok)
    }
    try_step <- function(code) tryCatch({ code; "" }, error = \(e) conditionMessage(e))
    for (t in c("tile-lipidomics", "btn-start-analysis")) {
        p <- k233Press(S, "LS1", t, "testid")
        if (!isTRUE(p$ok)) return(list(ok = FALSE, why = paste("LS1", p$reason)))
        k233Idle(S, 10)
    }
    if (!k233WaitFor(S, "document.getElementById('experiment_label') !== null", 30))
        return(list(ok = FALSE, why = "LS1: the Experiment Setup dialog did not open"))
    k233Set(S, "LS1", "experiment_label", K253_LIPID_LABEL)
    k233Set(S, "LS1", "project_dir", S$project_dir)
    p <- k233Press(S, "LS1", "btn-confirm-setup", "testid")
    ok <- isTRUE(p$ok) && k233WaitFor(S, "document.querySelector(\"[data-testid='lipid-tab-import']\") !== null", 60)
    if (!step("LS1 home and project", ok, p$reason %||% "")) return(list(ok = FALSE, why = "LS1"))
    imp <- function(x) paste0("lipidomics_workflow-import-", x)
    assays <- unlist(lane$assays)
    dir <- file.path(G$pkg, "tests/testdata/e2e", lane$fixture_dir)
    k233Set(S, "LS2", imp("assay1_name"), assays[[1L]])
    k233Set(S, "LS2", imp("vendor_format"), lane$import_tool)
    up1 <- k233Upload(S, "LS2", imp("assay1_file_std"), file.path(dir, lane$seed_file), timeout_s = 120)
    k233Set(S, "LS2", imp("assay2_name"), assays[[2L]])
    up2 <- k233Upload(S, "LS2", imp("assay2_file_std"), file.path(dir, lane$assay2_file), timeout_s = 120)
    k233Set(S, "LS2", imp("sanitize_names"), FALSE)
    p <- k233Press(S, "LS2", "lipid-import-process", "testid")
    ok <- isTRUE(up1$ok) && isTRUE(up2$ok) && isTRUE(p$ok) && k253WaitStatus(S, "lipidomics", "setup_import", 180)
    if (!step("LS2 import", ok, sprintf("uploads %s / %s; press %s", up1$text, up2$text, p$reason %||% "ok")))
        return(list(ok = FALSE, why = "LS2"))
    why <- try_step(H$e2e_complete_lipid_design_from_manifest(S$driver, lane, timeout = 180000))
    ok <- !nzchar(why) && k253WaitStatus(S, "lipidomics", "design_matrix", 180)
    if (!step("LS3 design", ok, why)) return(list(ok = FALSE, why = paste("LS3", why)))
    why <- try_step(H$e2e_run_lipid_qc_finalize(S$driver, timeout = 120000))
    ok <- !nzchar(why) && k253WaitStatus(S, "lipidomics", "quality_control", 180)
    if (!step("LS4 QC finalise", ok, why)) return(list(ok = FALSE, why = paste("LS4", why)))
    norm <- function(x) paste0("lipidomics_workflow-norm-", x)
    why <- try_step({
        H$e2e_switch_workflow_tab(S$driver, "lipidomics", "norm", timeout = 60000)
        H$e2e_wait_for_selector(S$driver, "lipid-tab-norm", timeout = 60000)
    })
    defaults <- sprintf("rendered: apply_itsd=%s norm_method=%s ruv_mode=%s", k233Clean(k233InputValue(S, norm("apply_itsd"))),
        k233Clean(k233InputValue(S, norm("norm_method"))), k233Clean(k233InputValue(S, norm("ruv_mode"))))
    # The lipidomics browser lane's settings (helper-e2e-browser.R e2e_run_lipid_normalization_export): ITSD off,
    # between-sample normalisation none, RUV skipped; the composite figure is the chain's last action either way.
    k233Set(S, "LS5", norm("apply_itsd"), FALSE)
    k233Set(S, "LS5", norm("norm_method"), "none")
    k233Set(S, "LS5", norm("ruv_mode"), "skip")
    if (!step("LS5 normalisation tab", !nzchar(why), why)) return(list(ok = FALSE, why = paste("LS5", why)))
    list(ok = TRUE, settings = paste(defaults, "| set: apply_itsd=FALSE norm_method=none ruv_mode=skip (the lane's)"))
}

#' The press of a sequence: none for a control; Reset on the normalisation tab; Save Design on the design tab, which is
#' made active right after the run is pressed (before the window), so that the press itself is one click.
k253PressFunction <- function(S, s) {
    click <- function(target, by) function() {
        sel <- k233Selector(S, target, by)
        r <- k233Js(S, sprintf(K233_PRESS_JS, k233Json(sel), "true"))
        k233LogAction(S, s$seq, target, "click (G0 press)")
        list(found = isTRUE(r$found), visible = isTRUE(r$visible), disabled = isTRUE(r$disabled),
            clicked = isTRUE(r$clicked), browser_ms = r$t %||% NA, target = target)
    }
    if (startsWith(s$press, "none")) return(list(press = NULL))
    if (identical(s$press, "Save Design")) {
        return(list(press = click("prot-design-save", "testid"), prepare = function() {
            tab <- k233ActivateTab(S, paste0(K233_NS, "workflow_tabs"), "design")
            seen <- k233Press(S, s$seq, "prot-design-save", "testid", click = FALSE)
            sprintf("%s; Save Design visible and enabled before the window: %s %s", tab, isTRUE(seen$ok), seen$reason %||% "")
        }))
    }
    target <- if (identical(s$omic, "proteomics")) paste0(K233_NS, "normalization-reset_normalization") else
        "lipidomics_workflow-norm-reset_normalization"
    list(press = click(target, "id"), prepare = function() {
        seen <- k233Press(S, s$seq, target, "id", click = FALSE)
        sprintf("Reset visible and enabled before the window: %s %s", isTRUE(seen$ok), seen$reason %||% "")
    })
}

#' The window: find the target request, press inside it, follow it to its completion.
k253Window <- function(G, S, W, s, press, known_before) {
    after_t <- if (is.na(s$after)) W$t0 else NULL
    target <- NULL
    out <- list(press_made = FALSE, press_inside = NA, press_note = "")
    limit <- 420
    repeat {
        k253Scan(W)
        fresh <- W$reqs[setdiff(names(W$reqs), known_before)]
        if (is.null(after_t)) {
            ruv <- Filter(\(r) identical(r$operation, s$after) && !is.na(r$status), fresh)
            if (length(ruv)) {
                first <- ruv[[which.min(vapply(ruv, \(r) as.numeric(r$completion_seen), numeric(1)))]]
                out$after_id <- first$id
                out$after_status <- first$status
                if (!identical(first$status, "complete")) {
                    out$instrument <- sprintf("%s ended %s before the window", s$after, first$status)
                    break
                }
                after_t <- first$completion_seen
            }
        }
        if (!is.null(after_t) && is.null(target)) {
            cand <- Filter(\(r) identical(r$operation, s$target) && as.numeric(r$seen) >= as.numeric(after_t), fresh)
            if (length(cand)) {
                target <- cand[[which.min(vapply(cand, \(r) as.numeric(r$seen), numeric(1)))]]$id
                # A control's snapshot of what the chain wrote before the target: rewritten files (logs, status files)
                # are not counted in its publication footprint. A press sequence takes none here (no delay before
                # the press).
                if (!is.function(press)) W$snapshot <- k253Files(W, W$t0)$path
            }
        }
        if (!is.null(target)) {
            r <- W$reqs[[target]]
            comp <- file.path(r$dir, "completion.rds")
            if (is.function(press) && !isTRUE(out$press_tried)) {
                out$press_tried <- TRUE
                before <- file.exists(comp)
                out$head_before_press <- W$head
                if (before) {
                    out$press_note <- "the target's completion.rds existed before the press could be made"
                } else {
                    t_press <- Sys.time()
                    p <- press()
                    after <- file.exists(comp)
                    out$press_t <- t_press
                    out$press_s <- k253Rel(W, t_press)
                    out$press_made <- isTRUE(p$clicked)
                    out$press_inside <- isTRUE(p$clicked) && !before && !after
                    out$press_note <- sprintf("control %s: found=%s visible=%s disabled=%s clicked=%s", p$target, p$found,
                        p$visible, p$disabled, p$clicked)
                    out$press_browser_ms <- p$browser_ms
                }
            }
            if (!is.na(W$reqs[[target]]$status)) break
        }
        if (k233Seconds(W$t0) > limit || k233Remaining(G) <= 0 || !S$proc$is_alive()) {
            out$window_timeout <- sprintf("no %s within %d s (alive %s)", if (is.null(target)) "target request" else
                "target completion", limit, S$proc$is_alive())
            break
        }
        Sys.sleep(K253_POLL)
    }
    out$target_id <- target %||% NA_character_
    if (!is.null(target)) {
        r <- W$reqs[[target]]
        out$target_request_s <- k253Rel(W, r$seen)
        out$target_completion_s <- if (is.na(r$status)) NA else k253Rel(W, r$completion_seen)
        out$target_status <- r$status
        out$target_expected_head <- r$expected_head
        out$target_t <- r$seen
    }
    out
}

#' After the target: follow until nothing is open for 15 s and at least 25 s have passed (a control also waits for its
#' completion notice), at most 240 s.
k253Settle <- function(G, S, W, s, window, since_ms) {
    started <- Sys.time()
    idle_since <- NULL
    done <- FALSE
    repeat {
        k253Scan(W)
        if (length(k253Open(W))) idle_since <- NULL else if (is.null(idle_since)) idle_since <- Sys.time()
        if (!done && startsWith(s$press, "none")) {
            st <- k233State(S, since_ms)
            done <- k233Noted(st %||% list(notes = list()), K253_DONE[[s$omic]], "message")
        }
        quiet <- !is.null(idle_since) && k233Seconds(idle_since) >= 15
        enough <- k233Seconds(started) >= 25 && (done || !startsWith(s$press, "none"))
        if ((quiet && enough) || k233Seconds(started) > 240 || k233Remaining(G) <= 0) break
        Sys.sleep(0.5)
    }
    list(settle_s = round(k233Seconds(started), 1))
}

k253Files <- function(W, since) {
    root <- W$S$project_dir
    files <- list.files(root, recursive = TRUE, all.files = TRUE, full.names = TRUE)
    stores <- normalizePath(W$stores, mustWork = FALSE)
    full <- normalizePath(files, mustWork = FALSE)
    inside <- if (length(stores)) vapply(full, \(f) any(startsWith(f, paste0(stores, "/"))), logical(1)) else
        rep(FALSE, length(full))
    info <- file.info(files[!inside])
    keep <- !is.na(info$mtime) & as.numeric(info$mtime) >= as.numeric(since)
    data.frame(path = substring(normalizePath(rownames(info)[keep]), nchar(normalizePath(root)) + 2L),
        bytes = info$size[keep], mtime_s = k253Rel(W, info$mtime[keep]), stringsAsFactors = FALSE)
}

k253Publications <- function(W, since) {
    root <- W$S$project_dir
    dirs <- list.dirs(root, recursive = TRUE, full.names = TRUE)
    stores <- normalizePath(W$stores, mustWork = FALSE)
    roots <- dirs[basename(dirs) == ".workflow-export" & !normalizePath(dirname(dirs), mustWork = FALSE) %in% stores]
    rows <- list()
    for (r in roots) for (d in list.dirs(r, recursive = FALSE)) {
        journal <- file.path(d, "journal.json")
        mtime <- max(file.info(list.files(d, recursive = TRUE, all.files = TRUE, full.names = TRUE))$mtime, file.mtime(d),
            na.rm = TRUE)
        j <- if (file.exists(journal)) tryCatch(jsonlite::read_json(journal), error = \(e) list()) else list()
        rows[[length(rows) + 1L]] <- list(dir = substring(normalizePath(d), nchar(normalizePath(root)) + 2L),
            status = j$status %||% "", state_id = j$state_id %||% "", last_write_s = k253Rel(W, mtime),
            after = as.numeric(mtime) >= as.numeric(since))
    }
    rows
}

k253Collect <- function(G, S, W, s, rec, since_ms) {
    out <- G$out
    stem <- file.path(out, gsub("[^A-Za-z0-9]", "", s$seq))
    st <- k233State(S, since_ms)
    notes <- st$notes %||% list()
    k233WriteTsv(lapply(notes, \(n) list(browser_s = round((as.numeric(n$t) - since_ms) / 1000, 2), type = n$type,
        text = n$text)), paste0(stem, "-notices.tsv"))
    jsonlite::write_json(st$dg, paste0(stem, "-digest.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
    k233WriteTsv(unname(k253RequestRows(W)), paste0(stem, "-requests.tsv"))
    k233WriteTsv(W$heads, paste0(stem, "-heads.tsv"))
    tryCatch(S$driver$get_screenshot(paste0(stem, "-end.png"), selector = "viewport"), error = \(e) NULL)
    since <- if (!is.null(rec$press_t)) rec$press_t else rec$target_t %||% W$t0
    files <- k253Files(W, since)
    k233WriteTsv(files, paste0(stem, "-files-after.tsv"))
    pubs <- k253Publications(W, since)
    k233WriteTsv(pubs, paste0(stem, "-publications.tsv"))
    texts <- vapply(notes, \(n) as.character(n$text), character(1))
    types <- vapply(notes, \(n) as.character(n$type), character(1))
    rec$completion_notice <- any(grepl(K253_DONE[[s$omic]], texts) & types == "message")
    rec$error_notices <- paste(texts[types == "error"], collapse = " || ")
    rec$notices_as_predicted <- if (startsWith(s$press, "none")) NA else
        all(vapply(K253_NOTICES[[s$omic]], \(p) any(grepl(p, texts, fixed = TRUE)), logical(1)))
    rec$files_after <- nrow(files)
    if (startsWith(s$press, "none")) {
        rec$footprint <- setdiff(files$path, W$snapshot %||% character())
        rec$footprint_size <- length(rec$footprint)
        k233WriteTsv(data.frame(path = rec$footprint), paste0(stem, "-footprint.tsv"))
    }
    rec$files_after_paths <- paste(utils::head(files$path, 40L), collapse = " ; ")
    rec$publications_after <- sum(vapply(pubs, \(p) isTRUE(p$after), logical(1)))
    rec$head_end <- W$head
    after_press <- Filter(\(h) !is.null(rec$press_s) && h$t >= rec$press_s, W$heads)
    rec$head_after_press <- if (length(after_press)) after_press[[1L]]$entity_id else NA_character_
    rec$head_moves_after_press <- length(after_press)
    rec$head_moves_in_window <- length(Filter(\(h) !is.null(rec$press_s) && h$t >= rec$press_s &&
        h$t <= (rec$target_completion_s %||% Inf), W$heads))
    chain <- if (identical(s$omic, "proteomics")) c("normalization-ruv", "normalization-qc-figures") else
        "^lipid_"
    later_chain <- Filter(\(r) !is.null(rec$target_t) && as.numeric(r$seen) > as.numeric(rec$target_t) &&
        identical(r$status, "complete") && any(vapply(chain, \(c) grepl(c, r$operation), logical(1))), W$reqs)
    rec$later_chain_complete <- length(later_chain)
    rec$digest_state_end <- k233Clean(st$dg$r6_current_state_per_omic[[s$omic]])
    rec$digest_norm_status_end <- k233Clean(st$dg$step_status_per_omic[[s$omic]][["normalization"]])
    later <- Filter(\(r) !is.null(rec$target_t) && as.numeric(r$seen) > as.numeric(rec$target_t), W$reqs)
    rec$later_requests <- paste(vapply(later, \(r) sprintf("%s=%s", r$operation, r$status), character(1)), collapse = "; ")
    rec$files <- files
    rec$app_alive_end <- S$proc$is_alive()
    k233StopSession(S)
    rec
}

# --- judging and records ---------------------------------------------------------------------------------------

k253Judge <- function(G) {
    R <- G$k253
    control <- function(rec) {
        if (is.null(rec) || !is.null(rec$instrument) || !is.null(rec$window_timeout)) return(FALSE)
        identical(rec$target_status, "complete") && isTRUE(rec$completion_notice) && isTRUE(rec$footprint_size > 0)
    }
    for (name in names(R)) {
        rec <- R[[name]]
        spec <- Filter(\(s) identical(s$seq, name), K253_SEQUENCES)[[1L]]
        if (startsWith(spec$press, "none")) {
            rec$verdict <- if (control(rec)) "control proven" else "control NOT proven (instrument)"
        } else {
            ctrl <- R[[if (identical(spec$omic, "proteomics")) "P0" else "L0"]]
            if (!is.null(rec$instrument) || !is.null(rec$window_timeout) || !isTRUE(rec$press_inside)) {
                rec$verdict <- "instrument defect (stop rule 3)"
            } else if (!control(ctrl)) {
                rec$verdict <- "not judged: the control is not proven"
            } else {
                footprint <- intersect(rec$files$path %||% character(), ctrl$footprint %||% character())
                rec$footprint_hits <- paste(footprint, collapse = " ; ")
                checks <- c(status = identical(rec$target_status, "stale"),
                    no_publication_dir = isTRUE(rec$publications_after == 0),
                    no_footprint_file = !length(footprint),
                    no_completion_notice = !(identical(spec$omic, "proteomics") && isTRUE(rec$completion_notice)),
                    head_moved_by_press = isTRUE(rec$head_moves_in_window >= 1L))
                rec$checks <- paste(names(checks), checks, sep = "=", collapse = " ")
                rec$verdict <- if (all(checks)) "as predicted" else "CONTRADICTED (stop rule 2)"
            }
        }
        R[[name]] <- rec
    }
    G$k253 <- R
    R
}

k253WriteSequences <- function(G) {
    R <- k253Judge(G)
    cols <- c("seq", "omic", "press", "target_operation", "verdict", "checks", "target_id", "target_request_s",
        "press_s", "target_completion_s", "press_made", "press_inside", "press_note", "target_status",
        "target_expected_head", "head_before_press", "head_after_press", "head_end", "head_moves_after_press",
        "head_at_run", "head_moves_in_window", "later_chain_complete", "publications_after", "files_after",
        "footprint_size", "footprint_hits", "completion_notice", "notices_as_predicted",
        "error_notices", "later_requests", "digest_state_end", "digest_norm_status_end", "after_id", "after_status",
        "prepare", "settings", "settle_s", "window_timeout", "instrument", "app_alive_end", "files_after_paths", "stores")
    rows <- lapply(R, \(r) stats::setNames(lapply(cols, \(k) {
        v <- r[[k]]
        if (inherits(v, "POSIXct")) format(v, "%H:%M:%OS2") else v
    }), cols))
    k233WriteTsv(unname(rows), file.path(G$out, "g0-sequences.tsv"))
}

k253Verdict <- function(G) {
    k253WriteSequences(G)
    R <- G$k253
    v <- vapply(R, \(r) r$verdict %||% r$instrument %||% "?", character(1))
    line <- if (length(R) < length(K253_SEQUENCES) || any(grepl("instrument|not judged|NOT proven|\\?", v))) {
        "G0: instrument defect or incomplete (stop rule 3)"
    } else if (any(grepl("CONTRADICTED", v))) {
        "G0: the reading is CONTRADICTED (stop rule 2)"
    } else {
        "G0: observed as predicted (stop rule 1)"
    }
    writeLines(c(line, paste(names(v), v, sep = ": ")), file.path(G$out, "verdict.txt"))
}

if (sys.nframe() == 0L) k253Main(commandArgs(trailingOnly = TRUE))
