# KIT-247 gate G0-4 (ticket r1.2, section 10): the browser click, observation of the base, no package code. One e2e
# lane (prot_lfq, MaxQuant label-free) through the application in Chrome (shinytest2 AppDriver, tests/testthat/
# helper-e2e-browser.R), as test-e2e-proteomics-tmt-lfq.R drives it up to the normalisation and its export; then the
# design tab: is Save Design visible and enabled, does a click commit a new design state, and what does the stepper say.
#
# Binding: v0.9-staging 16112917 plus this folder only (coordinator's binding, r1.2).
#
# PREDICTIONS (written before any value is read; regular expressions applied to the observed value):
k247BrowserPredictions <- data.frame(stringsAsFactors = FALSE, rbind(
    c("button.visible", "after the normalisation the design tab shows Save Design", "^TRUE$"),
    c("button.enabled", "and it is enabled (nothing disables it, 2.1)", "^TRUE$"),
    c("after.head", "the click commits the design state again: the head is the design state", "^protein_s4_initial$"),
    c("after.commits_increased", "a new state record is written to the store", "^TRUE$"),
    c("after.history_unchanged", "the later states stay in the active history", "^TRUE$"),
    c("after.status_unchanged_except_design", "the stepper does not change (QC and normalisation still complete)", "^TRUE$"),
    c("after.status", "the stepper says Quality Control complete", "quality_control=complete"),
    c("after.notes", "the success text", "saved successfully")))
names(k247BrowserPredictions) <- c("key", "predicted", "pattern")
#
# Usage (package root, on a runner, NOT_CRAN=true, shinytest2 installed):
#   Rscript docs/qa/workflow-execution/kit/kit-247/g0-browser.R <out dir>

k247Args <- commandArgs(trailingOnly = TRUE)
k247Out <- if (length(k247Args) >= 1L) k247Args[[1L]] else "kit247-g0b"
dir.create(k247Out, recursive = TRUE, showWarnings = FALSE)
k247Log <- function(...) cat(sprintf("[KIT-247 G0-4 %s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)), file = stderr())
k247Rows <- list()
k247Put <- function(key, value) {
    value <- gsub("[\t\r\n]+", " ", paste(as.character(unlist(value)), collapse = " | "))
    if (!nzchar(value)) value <- "<none>"
    k247Rows[[length(k247Rows) + 1L]] <<- data.frame(key = key, value = substr(value, 1L, 2000L), stringsAsFactors = FALSE)
    k247Log(key, " = ", substr(value, 1L, 300L))
}
k247Write <- function() {
    observed <- if (length(k247Rows)) do.call(rbind, k247Rows) else data.frame(key = character(), value = character())
    utils::write.table(observed, file.path(k247Out, "g0-g04-observed.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    compare <- k247BrowserPredictions
    compare$observed <- vapply(compare$key, function(key) {
        hit <- observed$value[observed$key == key]
        if (length(hit)) hit[[length(hit)]] else "<not observed>"
    }, "")
    compare$held <- ifelse(compare$observed == "<not observed>", NA,
        mapply(function(p, o) grepl(p, o, perl = TRUE), compare$pattern, compare$observed))
    utils::write.table(compare, file.path(k247Out, "g0-g04-compare.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

Sys.setenv(MULTISCHOLAR_E2E_ARTIFACT_DIR = normalizePath(k247Out))
k247Log("loading the package and the e2e helpers")
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
library(testthat)
k247Put("not_cran", Sys.getenv("NOT_CRAN"))
k247Put("browser_dependencies", e2e_browser_dependencies_available())
k247Put("chrome", tryCatch(chromote::find_chrome(), error = function(e) paste("<error>", conditionMessage(e))))

k247Commits <- function(root) {
    files <- list.files(root, recursive = TRUE, all.files = TRUE, no.. = TRUE)
    length(grep("(^|/)\\.workflow-entities/entities/[^/]+/commit\\.json$", files))
}
k247Status <- function(digest) {
    status <- digest$step_status_per_omic$proteomics
    if (is.null(status)) return("<none>")
    paste(paste0(names(status), "=", vapply(status, function(x) paste(as.character(unlist(x)), collapse = ","), "")), collapse = "; ")
}

lane_id <- "prot_lfq"
case_id <- "KIT-247-G0-4-prot_lfq"
started <- Sys.time()
k247Main <- function() {
    manifest <- read_e2e_manifest()
    lane <- manifest[[lane_id]]
    project_root <- file.path(tempdir(), paste0(case_id, "-", Sys.getpid()))
    dir.create(project_root, recursive = TRUE, showWarnings = FALSE)
    paths <- e2e_proteomics_project_paths(project_root, case_id)
    driver <- e2e_new_app_driver(case_id = case_id, launch_surface = "run_app", timeout = 60000L, load_timeout = 180000L)
    on.exit(e2e_stop_driver(driver), add = TRUE)
    k247Log("app started")
    # VAL-54 instrument rerun 3 of 3 (coordinator's ruling): every notification ever shown is kept in the page (as
    # KIT-233's journey driver keeps them, journey.R:671-687), and on any failure the application's and the browser's
    # logs, the kept notifications and the project's files are saved (k247Evidence).
    k247Hook <- paste0("(function(){if(window.__k247){return true;}var K=window.__k247={notes:[]};",
        "setInterval(function(){var ns=document.querySelectorAll('.shiny-notification');for(var i=0;i<ns.length;i++){",
        "var n=ns[i];var t=(n.textContent||'').trim();if(!t||n.__k247===t){continue;}n.__k247=t;",
        "K.notes.push({t:Date.now(),text:t.slice(0,2000),cls:n.className||''});}},200);return true;})()")
    k247InstallHook <- function() tryCatch(driver$get_js(k247Hook), error = function(e) NULL)
    k247InstallHook()
    k247Evidence <- function(what) {
        logs <- tryCatch(as.data.frame(driver$get_logs()), error = function(e) data.frame(error = conditionMessage(e)))
        logs[] <- lapply(logs, function(x) gsub("[\t\r\n]+", " ", as.character(x)))
        utils::write.table(logs, file.path(k247Out, "g0-g04-applog.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
        notes <- tryCatch(driver$get_js("(window.__k247 ? window.__k247.notes : [])"), error = function(e) list())
        notes <- if (length(notes)) do.call(rbind, lapply(notes, function(n) data.frame(t = as.character(n$t),
            text = gsub("[\t\r\n]+", " ", as.character(n$text)), cls = as.character(n$cls)))) else data.frame(text = "<none>")
        utils::write.table(notes, file.path(k247Out, "g0-g04-notes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
        files <- list.files(project_root, recursive = TRUE, all.files = TRUE, full.names = TRUE)
        info <- file.info(files)
        utils::write.table(data.frame(path = sub(paste0("^", project_root, "/?"), "", files), size = info$size,
            mtime = format(info$mtime, "%H:%M:%OS3")), file.path(k247Out, "g0-g04-files.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
        k247Put("evidence_saved_at", what)
    }
    k247NoteSeen <- function(pattern, timeout) {
        script <- sprintf("(window.__k247 ? window.__k247.notes : []).some(function(n){return new RegExp(%s,'i').test(n.text);})",
            jsonlite::toJSON(pattern, auto_unbox = TRUE))
        e2e_wait_for_js(driver, script, timeout = timeout)
    }
    on.exit(if (!identical(k247Done, TRUE)) tryCatch(k247Evidence("failure"), error = function(e) NULL), add = TRUE, after = FALSE)
    k247Done <- FALSE

    # The lane, as test-e2e-proteomics-tmt-lfq.R runs it, up to the normalisation and its export.
    e2e_wait_for_selector(driver, "tile-proteomics", timeout = 60000L)
    e2e_select_omics(driver, "proteomics", assert_digest = TRUE, timeout = 60000L)
    e2e_complete_project_setup(driver, case_id, project_root, timeout = 60000L)
    e2e_wait_for_selector(driver, "prot-tab-setup", timeout = 60000L)
    guiEntryInstallUniprotCache(paths)
    e2e_upload_lane_inputs(driver, lane)
    e2e_click_process_import(driver, "proteomics")
    e2e_wait_for_step_status(driver, "proteomics", "setup_import", "complete", timeout = 120000L)
    e2e_complete_prot_design_from_manifest(driver, lane, timeout = 120000L)
    e2e_wait_for_step_status(driver, "proteomics", "design_matrix", "complete", timeout = 120000L)
    e2e_run_prot_protein_qc(driver, timeout = 120000L)
    e2e_wait_for_step_status(driver, "proteomics", "quality_control", "complete", timeout = 120000L)
    # VAL-54 instrument fix 1 (run 45594085728 timed out waiting for the normalisation STATUS after
    # e2e_run_prot_normalization_export, which presses Skip right after Run with no wait for the run to finish): the same
    # presses, each waited for by the head state (digest) and the export by its notice. The status is recorded, not waited on.
    k247WaitHead <- function(states, what, timeout = 600000L) {
        deadline <- Sys.time() + timeout / 1000
        repeat {
            digest <- tryCatch(e2e_state_digest(driver, timeout = 60000L), error = function(e) NULL)
            head <- if (is.null(digest)) NULL else e2e_digest_r6_current_state(digest, "proteomics")
            if (isTRUE(head %in% states)) return(invisible(digest))
            if (Sys.time() >= deadline) {
                k247Put(paste0(what, ".timeout_head"), head %||% "<NULL>")
                k247Put(paste0(what, ".timeout_notes"), vapply(e2e_collect_notifications(driver), function(n) as.character(n$text %||% ""), ""))
                stop(sprintf("Timed out waiting for the head to be one of %s after %s; last head %s", paste(states, collapse = ", "),
                    what, head %||% "<NULL>"))
            }
            Sys.sleep(1)
        }
    }
    e2e_switch_workflow_tab(driver, "proteomics", "normalization")
    e2e_wait_for_selector(driver, "prot-tab-normalization", timeout = 60000L)
    e2e_set_input_and_idle(driver, e2e_input_id("proteomics", "norm", "ruv_mode"), "skip", timeout = 60000L)
    e2e_click_testid(driver, "prot-norm-run", timeout = 60000L)
    k247WaitHead(c("normalized", "ruv_corrected"), "norm_run")
    # VAL-54 instrument fix 2 (rerun 45597423706: the head was 'normalized' while Skip, pressed then, did nothing for
    # 300 s; the state is committed before the step ends, KIT-233 A5): Skip waits for the run's own completion notice.
    # Rerun 3: KIT-233's completion signal for this step (journey.R:1555-1558), from the kept notifications, 600 s.
    k247InstallHook()
    k247NoteSeen("Normali[sz]ation (and RUV correction )?completed", 600000L)
    e2e_wait_for_idle(driver, timeout = 120000L)
    digest <- e2e_state_digest(driver, timeout = 60000L)
    k247Put("norm_run.status", k247Status(digest))
    e2e_set_input_and_idle(driver, e2e_input_id("proteomics", "norm", "norm_qc_tabs"), "Correlation Filtering", timeout = 60000L)
    e2e_click_input_id(driver, e2e_input_id("proteomics", "norm", "skip_correlation_filter"), timeout = 60000L)
    digest <- k247WaitHead("correlation_filtered", "corr_skip")
    k247Put("corr_skip.status", k247Status(digest))
    k247Put("corr_skip.notes", vapply(e2e_collect_notifications(driver), function(n) as.character(n$text %||% ""), ""))
    e2e_click_norm_export(driver, "proteomics")
    e2e_wait_for_js(driver, sprintf(
        "Array.from(document.querySelectorAll('.shiny-notification')).some(function(n){return (n.textContent||'').indexOf(%s)!==-1;})",
        jsonlite::toJSON("Filtered session data exported successfully!", auto_unbox = TRUE)), timeout = 240000L)
    e2e_wait_for_idle(driver, timeout = 120000L)
    k247Log("lane reached the normalisation export")

    before <- e2e_state_digest(driver, timeout = 60000L)
    before_history <- e2e_digest_r6_state_history(before, "proteomics")
    before_status <- before$step_status_per_omic$proteomics
    before_commits <- k247Commits(project_root)
    k247Put("before.head", e2e_digest_r6_current_state(before, "proteomics"))
    k247Put("before.history", paste(before_history, collapse = " > "))
    k247Put("before.status", k247Status(before))
    k247Put("before.commits", before_commits)

    # The design tab after the normalisation.
    builder <- function(input) e2e_input_id("proteomics", "design", paste0("builder-", input))
    e2e_switch_workflow_tab(driver, "proteomics", "design_matrix")
    e2e_wait_for_selector(driver, "prot-tab-design", timeout = 60000L)
    button <- e2e_call_driver_method(driver, "get_js", script = sprintf(paste(
        "(function(){var b=document.getElementById(%s);if(!b){return 'absent';}",
        "var visible=!!(b.offsetWidth||b.offsetHeight||b.getClientRects().length);",
        "var disabled=b.disabled===true||b.classList.contains('disabled')||b.getAttribute('aria-disabled')==='true';",
        "return JSON.stringify({visible:visible,disabled:disabled});})()"),
        jsonlite::toJSON(builder("save_results"), auto_unbox = TRUE)))
    k247Put("button.raw", button)
    state <- tryCatch(jsonlite::fromJSON(button), error = function(e) list(visible = NA, disabled = NA))
    k247Put("button.visible", isTRUE(state$visible))
    k247Put("button.enabled", identical(state$disabled, FALSE))

    # Move KO_5 to WT (replicate 6, after WT's five) and press Save Design.
    e2e_set_input_and_idle(driver, builder("main_tabs"), "Assign Metadata", timeout = 60000L)
    available <- e2e_selectize_option_values(driver, builder("selected_runs"), timeout = 60000L)
    moved <- e2e_match_available_values("KO_5", available)
    k247Put("moved_run", moved)
    e2e_set_selectize_values(driver, builder("selected_runs"), moved, timeout = 60000L)
    e2e_wait_for_input_id(driver, builder("replicate_start"), timeout = 60000L)
    e2e_set_input_and_idle(driver, builder("replicate_start"), 6L, timeout = 60000L)
    factors <- e2e_selectize_option_values(driver, builder("factor1_select"), timeout = 60000L)
    e2e_set_selectize_values(driver, builder("factor1_select"), e2e_match_available_values("WT", factors), timeout = 60000L,
        multiple = FALSE)
    e2e_click_input_id(driver, builder("assign_metadata"), timeout = 60000L)
    e2e_click_input_id(driver, builder("save_results"), timeout = 120000L)
    notes_seen <- tryCatch({
        e2e_wait_for_js(driver, paste0("Array.from(document.querySelectorAll('.shiny-notification')).some(function(n){",
            "return (n.textContent||'').indexOf('saved successfully')!==-1;})"), timeout = 120000L)
        TRUE
    }, error = function(e) FALSE)
    e2e_wait_for_idle(driver, timeout = 120000L)
    notes <- e2e_collect_notifications(driver)
    k247Put("after.notes", vapply(notes, function(n) as.character(n$text %||% ""), ""))
    k247Put("after.success_note_seen", notes_seen)

    after <- e2e_state_digest(driver, timeout = 60000L)
    after_history <- e2e_digest_r6_state_history(after, "proteomics")
    after_status <- after$step_status_per_omic$proteomics
    k247Put("after.head", e2e_digest_r6_current_state(after, "proteomics"))
    k247Put("after.history", paste(after_history, collapse = " > "))
    k247Put("after.history_unchanged", identical(after_history, before_history))
    k247Put("after.status", k247Status(after))
    fields <- setdiff(union(names(before_status), names(after_status)), "design_matrix")
    k247Put("after.status_unchanged_except_design", identical(lapply(after_status[fields], unlist), lapply(before_status[fields], unlist)))
    after_commits <- k247Commits(project_root)
    k247Put("after.commits", after_commits)
    k247Put("after.commits_increased", after_commits > before_commits)
    k247Done <- TRUE
    k247Evidence("success")
    "ok"
}
result <- tryCatch(k247Main(), error = function(e) {
    k247Put("error", conditionMessage(e))
    "error"
})
k247Put("outcome", result)
k247Put("seconds", round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1))
k247Write()
k247Log("done")
