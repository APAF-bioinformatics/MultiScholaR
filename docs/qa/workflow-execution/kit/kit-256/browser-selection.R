# KIT-256 G0 (base) and GB (the fix): what the design builder's "Select Technical Replicates" box shows after
# 'Assign Technical Replicates', in a real browser. One script for both gates (ticket section 10, r1.2).
#
# Launch (KIT-233's form, origin/kit-233-journey:docs/qa/workflow-execution/kit/kit-233/journey.R:467, :607-630): the
# harness app (the installed package's lipid builder over a 4-run, one-group table) is written to a temporary folder and
# started as a plain `Rscript --vanilla` process at warn = 1 with shiny.testmode; HTTP 200 is awaited; then
# shinytest2::AppDriver attaches by URL. AppDriver$new(<function>) is not used: its application process runs at
# warn = 2 (KIT-259 G0, run 45610170354), where a stray warning could abort the Assign observer and fake the box.
#
# Phase 0 (lifecycle, VAL-54): Chrome is found (else an error, never a skip), the app answers, the driver attaches, the
# app's exported warn is below 2, the box offers S1 to S4. Any failure there is an instrument defect: the script stops
# with an error and no value is recorded as a result.
# Each counted click is checked to have acted before the box is read (the summary changes, or a notification appears
# for a click that leaves the design as it was). A click without its effect is an instrument failure.
# The app log is kept; its warnings are counted in the values as observations, never as failures.
#
# Usage (package root, on a runner, NOT_CRAN=true, shinytest2 installed): Rscript <this file> <out dir> <g0|gb>

k256Args <- commandArgs(trailingOnly = TRUE)
k256Out <- if (length(k256Args) >= 1L) k256Args[[1L]] else "kit256-browser"
k256Gate <- if (length(k256Args) >= 2L) k256Args[[2L]] else "g0"
if (!k256Gate %in% c("g0", "gb")) stop("the gate is g0 or gb, not ", k256Gate)
dir.create(k256Out, recursive = TRUE, showWarnings = FALSE)
k256Out <- normalizePath(k256Out)
options(warn = 1, chromote.timeout = 60)

# ---------------------------------------------------------------------------------------------------------------------
# PREDICTIONS (written before any value is read; regular expressions applied to the observed value).
# By reading at v0.9-staging 63143aa3 (R/mod_omic_design_builder.R:520-541, :650-656; R/func_omics_design_shared.R:636-686):
# on base the selection observer re-sends the kept selection after every design change, so the box keeps {S1, S2}
# after Assign, the user's added S3, S4 join them, and one Assign marks all four runs as one biological sample
# (stop rule 1 if the box is empty). On the fix the box is empty after each successful Assign, the added runs are
# S3, S4 alone, and two samples are saved (stop rule 2 if the box is still filled). A refused Assign keeps the box.
# Summaries are the output's text with each line break and its indentation replaced by " | ".
# ---------------------------------------------------------------------------------------------------------------------
k256P <- function(key, gate, predicted, pattern) data.frame(key = key, gate = gate, predicted = predicted,
    pattern = pattern, stringsAsFactors = FALSE)
k256Sum1 <- "Group: [^,|]+, Biological Replicate: 1 \\| Samples: S1, S2 \\| Technical Replicates: 1, 2"
k256Predictions <- rbind(
    k256P("p0.chrome", "both", "a Chrome binary is found", "."),
    k256P("p0.answered", "both", "the app process answers HTTP 200", "^TRUE$"),
    k256P("p0.attached", "both", "AppDriver attaches by URL", "^TRUE$"),
    k256P("p0.app_warn", "both", "the app process runs below warn 2 (set to 1)", "^[01]$"),
    k256P("p0.options", "both", "the box offers S1 to S4", "^S1,S2,S3,S4$"),
    k256P("step1.groups", "both", "Assign Metadata gave the runs a group (contrast choices hold it)", "."),
    k256P("step2.set", "both", "the box holds S1, S2 before the click", '^\\["S1","S2"\\]$'),
    k256P("step2.acted", "both", "the summary changed after Assign", "^TRUE$"),
    k256P("step2.summary", "both", "one biological replicate {S1, S2}", paste0("^", k256Sum1, "$")),
    k256P("step3.box", "g0", "base: the box keeps S1, S2", '^\\["S1","S2"\\]$'),
    k256P("step3.box", "gb", "fix: the box is empty", "^\\[\\]$"),
    k256P("step4.sent", "g0", "base: the box then holds S1 to S4", "^S1,S2,S3,S4$"),
    k256P("step4.sent", "gb", "fix: the box then holds S3, S4 only", "^S3,S4$"),
    k256P("step4.acted", "both", "the summary changed after Assign", "^TRUE$"),
    k256P("step4.summary", "g0", "base: one biological replicate of four runs",
        "^Group: [^,|]+, Biological Replicate: 1 \\| Samples: S1, S2, S3, S4 \\| Technical Replicates: 1, 2, 3, 4$"),
    k256P("step4.summary", "gb", "fix: two biological replicates, {S1, S2} and {S3, S4}",
        paste0("^", k256Sum1, " \\| Group: [^,|]+, Biological Replicate: 2 \\| Samples: S3, S4 \\| ",
            "Technical Replicates: 1, 2$")),
    k256P("step4.box", "g0", "base: the box keeps the four runs", '^\\["S1","S2","S3","S4"\\]$'),
    k256P("step4.box", "gb", "fix: the box is empty", "^\\[\\]$"),
    k256P("step5.box_before", "both", "the box holds S1 before the refused click", '^\\["S1"\\]$'),
    k256P("step5.acted", "both", "a warning notification appears (too few runs)", "warning.*at least two"),
    k256P("step5.summary_unchanged", "both", "the refused Assign leaves the summary", "^TRUE$"),
    k256P("step5.box_after", "both", "the refused Assign leaves the box", '^\\["S1"\\]$'),
    k256P("step6.keys", "g0", "base: the saved design has one key of four runs", "^[^ =]+_1=\\{S1,S2,S3,S4\\}$"),
    k256P("step6.keys", "gb", "fix: the saved design has two keys, {S1, S2} and {S3, S4}",
        "^[^ =]+_1=\\{S1,S2\\} [^ =]+_2=\\{S3,S4\\}$"),
    k256P("nc.step2.summary", "gb", "second session: one biological replicate {S1, S2}", paste0("^", k256Sum1, "$")),
    k256P("nc.box_set", "gb", "second session: S1, S2 set again", '^\\["S1","S2"\\]$'),
    k256P("nc.acted", "gb", "the no-change Assign shows its notification", "Assigned 2 samples"),
    k256P("nc.summary_unchanged", "gb", "the design is unchanged (identical summary)", "^TRUE$"),
    k256P("nc.box", "gb", "fix: the no-change Assign also empties the box", "^\\[\\]$"))
k256Predictions <- k256Predictions[k256Predictions$gate %in% c("both", k256Gate), ]

# ---------------------------------------------------------------------------------------------------------------------
k256Log <- function(...) cat(sprintf("[KIT-256 %s %s] %s\n", k256Gate, format(Sys.time(), "%H:%M:%S"),
    paste0(...)), file = stderr())
k256Rows <- new.env(parent = emptyenv())
k256Rows$rows <- list()
k256Text <- function(value) {
    if (is.null(value)) return("<NULL>")
    paste(as.character(unlist(value, use.names = FALSE)), collapse = ",")
}
k256Put <- function(key, value) {
    value <- k256Text(value)
    k256Rows$rows[[length(k256Rows$rows) + 1L]] <- data.frame(key = key, value = value, stringsAsFactors = FALSE)
    k256Log(key, " = ", substr(value, 1L, 400L))
    invisible(value)
}
k256Write <- function() {
    observed <- do.call(rbind, k256Rows$rows)
    if (is.null(observed)) observed <- data.frame(key = character(), value = character())
    utils::write.table(observed, file.path(k256Out, paste0(k256Gate, "-values.tsv")), sep = "\t", quote = FALSE,
        row.names = FALSE)
    compare <- k256Predictions
    compare$observed <- vapply(compare$key, function(key) {
        hit <- observed$value[observed$key == key]
        if (length(hit)) hit[[length(hit)]] else "<not observed>"
    }, "")
    compare$held <- ifelse(compare$observed == "<not observed>", NA,
        mapply(function(p, o) grepl(p, o, perl = TRUE), compare$pattern, compare$observed))
    utils::write.table(compare, file.path(k256Out, paste0(k256Gate, "-compare.tsv")), sep = "\t", quote = FALSE,
        row.names = FALSE)
    k256Log("written: ", sum(compare$held %in% TRUE), " held, ", sum(compare$held %in% FALSE), " not held, ",
        sum(is.na(compare$held)), " not observed")
}
utils::write.table(k256Predictions, file.path(k256Out, paste0(k256Gate, "-predictions.tsv")), sep = "\t",
    quote = FALSE, row.names = FALSE)
# The end of the run, once: the drivers and the app process stopped, the app log's counts, the tables written.
k256Rows$finished <- FALSE
k256Finish <- function() {
    if (isTRUE(k256Rows$finished)) return(invisible(NULL))
    k256Rows$finished <- TRUE
    for (app in get0("k256Drivers", envir = globalenv(), ifnotfound = list())) {
        tryCatch(app$stop(), error = function(e) NULL)
    }
    proc <- get0("k256Proc", envir = globalenv())
    if (!is.null(proc)) tryCatch(proc$kill_tree(), error = function(e) tryCatch(proc$kill(), error = function(e) NULL))
    log <- get0("k256AppLog", envir = globalenv(), ifnotfound = "")
    lines <- if (nzchar(log) && file.exists(log)) readLines(log, warn = FALSE) else character()
    k256Put("applog.lines", length(lines))
    k256Put("applog.warnings", sum(grepl("^Warning", lines)))
    k256Put("applog.errors", sum(grepl("^Error", lines)))
    k256Write()
}
k256Instrument <- function(...) {
    message <- paste0(...)
    k256Put("instrument_failure", message)
    k256Finish()
    stop("KIT-256 instrument failure (no result): ", message, call. = FALSE)
}
# Any error (an instrument failure included) stops the app, writes the tables once and ends the script with status 1
# (in a non-interactive session an error handler must quit itself, or R carries on with the next expression).
options(error = function() {
    tryCatch(k256Finish(), error = function(e) NULL)
    q(save = "no", status = 1L)
})

k256Put("gate", k256Gate)
k256Put("package_version", as.character(utils::packageVersion("MultiScholaR")))
k256Put("package_sha", tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) "<unknown>"))
k256Put("shinytest2_version", as.character(utils::packageVersion("shinytest2")))
k256Put("chromote_version", as.character(utils::packageVersion("chromote")))
k256Put("shiny_version", as.character(utils::packageVersion("shiny")))

# ---------------------------------------------------------------------------------------------------------------------
# Phase 0: Chrome, the app process, the driver, the app's warn level, the box's options.
# ---------------------------------------------------------------------------------------------------------------------
k256Chrome <- tryCatch(chromote::find_chrome(), error = function(e) NULL)
cat("chromote::find_chrome(): ", if (is.null(k256Chrome)) "<none>" else k256Chrome, "\n", sep = "")
if (is.null(k256Chrome) || !any(nzchar(k256Chrome))) k256Instrument("no Chrome found by chromote::find_chrome()")
k256Put("p0.chrome", k256Chrome)

k256AppDir <- file.path(tempfile("kit256-app-"))
dir.create(k256AppDir, recursive = TRUE)
writeLines(c(
    "# KIT-256 harness app: the installed package's lipid design builder over a 4-run, one-group table.",
    "runs <- paste0('S', 1:4)",
    "values <- as.data.frame(matrix(c(11, 12, 21, 22, 31, 32, 41, 42), nrow = 2L, dimnames = list(NULL, runs)))",
    "assays <- list(LCMS_Pos = cbind(data.frame(LipidName = c('L1', 'L2'), LipidClass = 'PC'), values))",
    "mapping <- list(lipid_id_col = 'LipidName', annotation_col = 'LipidClass', sample_columns = runs,",
    "    is_pattern = NA_character_)",
    "ui <- shiny::fluidPage(MultiScholaR::mod_lipid_design_builder_ui('b'))",
    "server <- function(input, output, session) {",
    "    result <- MultiScholaR::mod_lipid_design_builder_server('b', data_tbl = function() assays,",
    "        config_list = shiny::reactive(list(deAnalysisParameters = list(formula_string = '~ 0 + group'))),",
    "        column_mapping = shiny::reactive(mapping))",
    "    shiny::exportTestValues(keys = {",
    "        saved <- result()",
    "        if (is.null(saved)) NULL else split(as.character(saved$design_matrix$Run), saved$design_matrix$tech_rep_group)",
    "    }, warn = getOption('warn'))",
    "}",
    "shiny::shinyApp(ui, server)"), file.path(k256AppDir, "app.R"))
file.copy(file.path(k256AppDir, "app.R"), file.path(k256Out, "harness-app.R"), overwrite = TRUE)

k256Port <- httpuv::randomPort()
k256AppLog <- file.path(k256Out, paste0(k256Gate, "-app.log"))
k256Expr <- sprintf(paste0(".libPaths(c(%s)); options(warn = 1, shiny.testmode = TRUE); ",
    "shiny::runApp('%s', port = %d, launch.browser = FALSE)"),
    paste(sprintf("'%s'", .libPaths()), collapse = ", "), k256AppDir, k256Port)
k256Proc <- processx::process$new(file.path(R.home("bin"), "Rscript"), c("--vanilla", "-e", k256Expr),
    stdout = k256AppLog, stderr = "2>&1", cleanup = TRUE, cleanup_tree = TRUE)
k256Url <- sprintf("http://127.0.0.1:%d", k256Port)
k256Drivers <- list()

k256Answered <- FALSE
for (i in seq_len(180L)) {
    if (!k256Proc$is_alive()) break
    k256Answered <- isTRUE(tryCatch(curl::curl_fetch_memory(k256Url,
        handle = curl::new_handle(timeout = 5))$status_code == 200L, error = function(e) FALSE))
    if (k256Answered) break
    Sys.sleep(1)
}
k256Put("p0.answered", k256Answered)
if (!k256Answered) {
    k256Instrument("the app did not answer HTTP 200 within 180 s (alive: ", k256Proc$is_alive(), ")")
}

k256Attach <- function() {
    attachError <- ""
    refused <- function(e) {
        attachError <<- conditionMessage(e)
        NULL
    }
    app <- tryCatch(shinytest2::AppDriver$new(k256Url, load_timeout = 120000, timeout = 60000, width = 1440,
        height = 1100), error = refused, skip = refused)
    if (is.null(app)) {
        k256Put("p0.attach_error", attachError)
        k256Instrument("AppDriver did not attach: ", attachError)
    }
    k256Drivers[[length(k256Drivers) + 1L]] <<- app
    # A notification recorder: every Shiny notification added to the page, with its class and text.
    app$run_js(paste0("window.k256Notes = []; new MutationObserver(function(ms) { ms.forEach(function(m) {",
        "m.addedNodes.forEach(function(n) { if (n.classList && n.classList.contains('shiny-notification')) ",
        "window.k256Notes.push(n.className + ' :: ' + n.textContent.replace(/\\s+/g, ' ').trim()); }); }); })",
        ".observe(document.body, {childList: true, subtree: true});"))
    app
}

k256Js <- function(app, js) app$get_js(js)
k256OptsJs <- function(id) sprintf(paste0("(function() { var el = document.getElementById('%s'); if (!el) return '';",
    " var v = el.selectize ? Object.keys(el.selectize.options) : Array.from(el.options).map(function(o) {",
    " return o.value; }); return v.filter(function(x) { return x !== ''; }).sort().join(','); })()"), id)
k256BoxJs <- "JSON.stringify(document.getElementById('b-tech_rep_samples').selectize.getValue())"
k256SummaryJs <- paste0("(function() { var el = document.getElementById('b-tech_rep_summary');",
    " return el ? el.textContent : ''; })()")
k256NotesJs <- "JSON.stringify(window.k256Notes || [])"
k256Box <- function(app) k256Js(app, k256BoxJs)
k256Summary <- function(app) gsub("\\s*\n\\s*", " | ", trimws(k256Js(app, k256SummaryJs)))
k256Notes <- function(app) jsonlite::fromJSON(k256Js(app, k256NotesJs))
k256Idle <- function(app) app$wait_for_idle(duration = 750, timeout = 60000)
k256Set <- function(app, ...) {
    app$set_inputs(..., wait_ = FALSE)
    k256Idle(app)
}
k256Click <- function(app, id) {
    app$click(input = id, wait_ = FALSE)
    k256Idle(app)
}
# Poll `condition()` up to `seconds`; TRUE when it held.
k256Until <- function(condition, seconds = 30) {
    deadline <- Sys.time() + seconds
    repeat {
        if (isTRUE(tryCatch(condition(), error = function(e) FALSE))) return(TRUE)
        if (Sys.time() > deadline) return(FALSE)
        Sys.sleep(0.25)
    }
}

k256App <- k256Attach()
k256Put("p0.attached", TRUE)
k256Warn <- k256App$get_value(export = "warn")
k256Put("p0.app_warn", k256Warn)
if (!is.numeric(k256Warn) || length(k256Warn) != 1L || k256Warn >= 2) {
    k256Instrument("the app process's warn is not below 2: ", k256Text(k256Warn))
}
k256Until(function() identical(k256Js(k256App, k256OptsJs("b-tech_rep_samples")), "S1,S2,S3,S4"), 60)
k256Put("p0.options", k256Js(k256App, k256OptsJs("b-tech_rep_samples")))
if (!identical(k256Js(k256App, k256OptsJs("b-tech_rep_samples")), "S1,S2,S3,S4")) {
    k256Instrument("the box does not offer S1 to S4")
}
k256Log("phase 0 passed")

# ---------------------------------------------------------------------------------------------------------------------
# Steps 1 and 2 (shared by both sessions): factor A; Assign Metadata S1 to S4 from 1; the box {S1, S2}, 'lowest', Assign.
# ---------------------------------------------------------------------------------------------------------------------
k256Steps12 <- function(app, prefix) {
    k256Set(app, `b-main_tabs` = "Factors")
    k256Set(app, `b-new_factor` = "A")
    k256Click(app, "b-add_factor")
    if (!k256Until(function() grepl("(^|,)A(,|$)", k256Js(app, k256OptsJs("b-factor1_select"))))) {
        k256Instrument(prefix, "step 1: factor A not offered after Add Factor")
    }
    k256Set(app, `b-main_tabs` = "Assign Metadata")
    k256Set(app, `b-selected_runs` = paste0("S", 1:4))
    if (!k256Until(function() isTRUE(k256Js(app, "!!document.getElementById('b-replicate_start')")))) {
        k256Instrument(prefix, "step 1: the starting replicate input did not render")
    }
    k256Set(app, `b-factor1_select` = "A", `b-replicate_start` = 1)
    k256Click(app, "b-assign_metadata")
    if (!k256Until(function() nzchar(k256Js(app, k256OptsJs("b-contrast_group1"))))) {
        k256Instrument(prefix, "step 1: Assign Metadata gave no group")
    }
    k256Put(paste0(prefix, "step1.groups"), k256Js(app, k256OptsJs("b-contrast_group1")))
    k256Set(app, `b-main_tabs` = "Technical Replicates")
    k256Until(function() nzchar(k256Summary(app)))
    before <- k256Summary(app)
    k256Put(paste0(prefix, "step1.summary"), before)
    k256Set(app, `b-tech_rep_samples` = c("S1", "S2"), `b-tech_rep_assignment_mode` = "lowest")
    k256Until(function() identical(k256Box(app), '["S1","S2"]'), 10)
    k256Put(paste0(prefix, "step2.set"), k256Box(app))
    k256Put(paste0(prefix, "step2.server_input"), app$get_value(input = "b-tech_rep_samples"))
    if (!identical(k256Box(app), '["S1","S2"]')) k256Instrument(prefix, "step 2: the box did not take S1, S2")
    k256Click(app, "b-assign_tech_reps")
    acted <- k256Until(function() !identical(k256Summary(app), before))
    k256Put(paste0(prefix, "step2.acted"), acted)
    if (!acted) k256Instrument(prefix, "step 2: Assign had no effect on the summary")
    k256Idle(app)
    k256Put(paste0(prefix, "step2.summary"), k256Summary(app))
    invisible(app)
}

k256Steps12(k256App, "")

# Step 3: the box after the successful Assign.
k256Box3 <- k256Box(k256App)
k256Put("step3.box", k256Box3)
k256Summary3 <- k256Summary(k256App)

# Step 4: the box set to what it shows together with S3 and S4 (what the user's box then holds); Assign.
k256Sent4 <- union(as.character(jsonlite::fromJSON(k256Box3)), c("S3", "S4"))
k256Put("step4.sent", k256Sent4)
k256Set(k256App, `b-tech_rep_samples` = k256Sent4)
k256Until(function() identical(k256Box(k256App), as.character(jsonlite::toJSON(k256Sent4))), 10)
k256Put("step4.set", k256Box(k256App))
if (!identical(k256Box(k256App), as.character(jsonlite::toJSON(k256Sent4)))) {
    k256Instrument("step 4: the box did not take ", k256Text(k256Sent4))
}
k256Click(k256App, "b-assign_tech_reps")
k256Acted4 <- k256Until(function() !identical(k256Summary(k256App), k256Summary3))
k256Put("step4.acted", k256Acted4)
if (!k256Acted4) k256Instrument("step 4: Assign had no effect on the summary")
k256Idle(k256App)
k256Put("step4.summary", k256Summary(k256App))
k256Put("step4.box", k256Box(k256App))

# Step 5: the box {S1}, Assign (refused: too few runs); the box before and after.
k256Summary5 <- k256Summary(k256App)
k256Set(k256App, `b-tech_rep_samples` = "S1")
k256Until(function() identical(k256Box(k256App), '["S1"]'), 10)
k256Put("step5.box_before", k256Box(k256App))
k256Notes5 <- length(k256Notes(k256App))
k256Click(k256App, "b-assign_tech_reps")
k256Acted5 <- k256Until(function() length(k256Notes(k256App)) > k256Notes5)
k256Put("step5.acted", if (k256Acted5) utils::tail(k256Notes(k256App), 1L) else "<no notification>")
if (!k256Acted5) k256Instrument("step 5: the refused Assign showed no notification")
k256Idle(k256App)
k256Put("step5.summary_unchanged", identical(k256Summary(k256App), k256Summary5))
k256Put("step5.box_after", k256Box(k256App))

# Step 6: Save; the saved design's keys (exported by the harness app).
k256Click(k256App, "b-save_results")
k256Saved <- k256Until(function() !is.null(k256App$get_value(export = "keys")))
k256Keys <- if (k256Saved) k256App$get_value(export = "keys") else NULL
k256Keys <- if (is.null(k256Keys)) "<no saved design>" else {
    k256Keys <- k256Keys[order(names(k256Keys))]
    paste(sprintf("%s={%s}", names(k256Keys), vapply(k256Keys, paste, "", collapse = ",")), collapse = " ")
}
k256Put("step6.keys", k256Keys)
k256Put("notifications", paste(k256Notes(k256App), collapse = " || "))

# ---------------------------------------------------------------------------------------------------------------------
# GB only, a second app session: the no-change path (S1, S2 again, the same mode, after step 2).
# ---------------------------------------------------------------------------------------------------------------------
if (identical(k256Gate, "gb")) {
    k256App2 <- k256Attach()
    k256Until(function() identical(k256Js(k256App2, k256OptsJs("b-tech_rep_samples")), "S1,S2,S3,S4"), 60)
    k256Put("nc.p0.options", k256Js(k256App2, k256OptsJs("b-tech_rep_samples")))
    k256Steps12(k256App2, "nc.")
    k256SummaryNc <- k256Summary(k256App2)
    k256Put("nc.box_after_step2", k256Box(k256App2))
    k256Set(k256App2, `b-tech_rep_samples` = c("S1", "S2"), `b-tech_rep_assignment_mode` = "lowest")
    k256Until(function() identical(k256Box(k256App2), '["S1","S2"]'), 10)
    k256Put("nc.box_set", k256Box(k256App2))
    if (!identical(k256Box(k256App2), '["S1","S2"]')) k256Instrument("no-change: the box did not take S1, S2")
    k256NotesNc <- length(k256Notes(k256App2))
    k256Click(k256App2, "b-assign_tech_reps")
    k256ActedNc <- k256Until(function() length(k256Notes(k256App2)) > k256NotesNc)
    k256Put("nc.acted", if (k256ActedNc) utils::tail(k256Notes(k256App2), 1L) else "<no notification>")
    if (!k256ActedNc) k256Instrument("no-change: the Assign showed no notification")
    k256Idle(k256App2)
    k256Put("nc.summary_unchanged", identical(k256Summary(k256App2), k256SummaryNc))
    k256Put("nc.box", k256Box(k256App2))
}
k256Log("done")
k256Finish()
