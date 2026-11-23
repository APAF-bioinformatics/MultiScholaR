# The default-user driver's parts (KIT-159 section 7.3): the defaults reader, the input state and
# its ledger, modals, step outcomes, forced options, helper extraction and the matrix stage rule.

LEDGER_SOURCES <- c("ui_default", "server_update", "structural", "forced", "override")
MATRIX_DISPLAY_LABEL <- "display; not what is tested"

#' Has an HTML node a class?
#' @noRd
valHasClass <- function(node, cls) {
    classes <- strsplit(xml2::xml_attr(node, "class") %||% "", "\\s+")[[1L]]
    cls %in% classes
}

#' Read rendered UI as a browser would send its inputs (A01(b))
#'
#' `input` (text, number, checkbox, slider `data-from`/`data-to`), radio groups (the checked value,
#' or `NULL`), `select` (single: the selected option, else the first; multiple: every selected
#' option, or `NULL`), `textarea` (its text) and checkbox groups (the checked values, or `NULL`).
#'
#' @param ui_html Rendered HTML (a string or a tag list).
#' @return A named list, one element per input id (a `NULL` element where the browser sends null).
#' @export
valUiDefaults <- function(ui_html) {
    html <- paste(as.character(ui_html), collapse = "\n")
    doc <- xml2::read_html(paste0("<html><body>", html, "</body></html>"))
    out <- list()
    put <- function(id, value) out[id] <<- list(value)  # nolint
    for (g in xml2::xml_find_all(doc, "//div[@id]")) {
        id <- xml2::xml_attr(g, "id")
        if (valHasClass(g, "shiny-input-checkboxgroup")) {
            v <- xml2::xml_attr(xml2::xml_find_all(g, ".//input[@checked]"), "value")
            put(id, if (length(v)) v else NULL)
        } else if (valHasClass(g, "shiny-input-radiogroup")) {
            v <- xml2::xml_attr(xml2::xml_find_all(g, ".//input[@checked]"), "value")
            put(id, if (length(v)) v[[1L]] else NULL)
        }
    }
    for (n in xml2::xml_find_all(doc, "//input[@id]")) {
        id <- xml2::xml_attr(n, "id")
        type <- tolower(xml2::xml_attr(n, "type") %||% "text")
        if (is.na(type)) type <- "text"
        if (valHasClass(n, "js-range-slider")) {
            from <- as.numeric(xml2::xml_attr(n, "data-from"))
            to <- xml2::xml_attr(n, "data-to")
            put(id, if (is.na(to)) from else c(from, as.numeric(to)))
        } else if (type == "checkbox") {
            put(id, !is.na(xml2::xml_attr(n, "checked")))
        } else if (type == "number") {
            put(id, as.numeric(xml2::xml_attr(n, "value")))
        } else if (type == "file") {
            put(id, NULL)
        } else if (type != "radio") {
            v <- xml2::xml_attr(n, "value")
            put(id, if (is.na(v)) "" else v)
        }
    }
    for (n in xml2::xml_find_all(doc, "//select[@id]")) {
        opts <- xml2::xml_find_all(n, ".//option")
        values <- xml2::xml_attr(opts, "value")
        chosen <- values[!is.na(xml2::xml_attr(opts, "selected"))]
        if (!is.na(xml2::xml_attr(n, "multiple"))) {
            put(xml2::xml_attr(n, "id"), if (length(chosen)) chosen else NULL)
        } else {
            put(xml2::xml_attr(n, "id"), if (length(chosen)) chosen[[1L]] else values[1L])
        }
    }
    for (n in xml2::xml_find_all(doc, "//textarea[@id]")) {
        put(xml2::xml_attr(n, "id"), xml2::xml_text(n))
    }
    out
}

#' A new input state: the values the driver holds, each with its source and location
#'
#' @return An environment.
#' @export
valInputState <- function() {
    st <- new.env(parent = emptyenv())
    st$values <- list()
    st$source <- character()
    st$location <- character()
    st
}

#' Store one value with its source and location
#' @noRd
valStore <- function(st, id, value, source, location) {
    st$values[id] <- list(value)
    st$source[[id]] <- source
    st$location[[id]] <- location
    invisible(st)
}

#' Apply rendered defaults: an input takes its rendered value when it first appears, and again
#' when its `uiOutput` re-renders (A33); a server-set value is never replaced otherwise
#'
#' @param st Input state.
#' @param defaults Named list from [valUiDefaults()].
#' @param location `ui`, `uiOutput:<id>` or `modal`.
#' @param rerender The `uiOutput` re-rendered.
#' @return The state, invisibly.
#' @export
valApplyRendered <- function(st, defaults, location = "ui", rerender = FALSE) {
    for (id in names(defaults)) {
        if (rerender || !id %in% names(st$values)) {
            valStore(st, id, defaults[[id]], "ui_default", location)
        }
    }
    invisible(st)
}

#' Record a server update; a choices-only update selects the first choice, as a browser does
#'
#' @param st Input state.
#' @param id Input id.
#' @param value,choices,selected As the `update*Input` call gave them.
#' @return The state, invisibly.
#' @export
valApplyServerUpdate <- function(st, id, value = NULL, choices = NULL, selected = NULL) {
    v <- selected %||% value %||% (if (length(choices)) unname(choices)[[1L]])
    loc <- st$location[id] %||% "server"
    if (is.na(loc)) loc <- "server"
    valStore(st, id, v, "server_update", loc)
}

#' Set an input with a source; a static default may never replace a server-set value (G12)
#'
#' @param st Input state.
#' @param id Input id.
#' @param value Value.
#' @param source One of the ledger sources.
#' @return The state, invisibly.
#' @export
valSetInput <- function(st, id, value, source) {
    if (!source %in% LEDGER_SOURCES) {
        valAbort("an input source outside the ledger's", "val_spec_invalid")
    }
    if (identical(source, "ui_default") && identical(unname(st$source[id]), "server_update")) {
        valAbort("a static default would overwrite the server-set value of {id}",
                 "val_spec_invalid")
    }
    loc <- st$location[id]
    valStore(st, id, value, source, if (is.na(loc %||% NA)) "ui" else loc)
}

#' Set an input from an action script: structural inputs, declared forced settings and arm
#' overrides only (section 7.3)
#'
#' @param st Input state.
#' @param id Input id.
#' @param value Value.
#' @param structural,forced,override The ids the script, the declaration and the arm may set.
#' @return The state, invisibly.
#' @export
valScriptSet <- function(st, id, value, structural, forced = character(), override = character()) {
    source <- if (id %in% structural) {
        "structural"
    } else if (id %in% forced) {
        "forced"
    } else if (id %in% override) {
        "override"
    } else {
        valAbort("the action script sets {id}, which is not structural, forced or an override",
                 "val_spec_invalid")
    }
    valSetInput(st, id, value, source)
}

#' An input's value
#' @param st Input state.
#' @param id Input id.
#' @return The value.
#' @export
valInputValue <- function(st, id) st$values[[id]]

#' An input's source
#' @param st Input state.
#' @param id Input id.
#' @return The source.
#' @export
valInputSource <- function(st, id) unname(st$source[id])

#' The ledger rows of the full input vector before an action (`ledger-<route>.tsv`)
#'
#' @param st Input state.
#' @param step,action Step id and the action about to be taken.
#' @return A data frame: step, action, input_id, value, source, location.
#' @export
valLedgerSnapshot <- function(st, step, action) {
    ids <- names(st$values) %||% character()
    value <- vapply(ids, function(id) {
        v <- st$values[[id]]
        if (is.null(v)) "NULL" else paste(format(v, digits = 15), collapse = "|")
    }, character(1))
    data.frame(step = rep(step, length(ids)), action = rep(action, length(ids)), input_id = ids,
               value = unname(value), source = unname(st$source[ids]),
               location = unname(st$location[ids]), stringsAsFactors = FALSE)
}

#' Apply a modal's rendered inputs (location `modal`, A18)
#'
#' @param st Input state.
#' @param modal_html The modal's rendered HTML.
#' @return The state, invisibly.
#' @export
valApplyModal <- function(st, modal_html) {
    valApplyRendered(st, valUiDefaults(modal_html), location = "modal", rerender = TRUE)
}

#' Press a modal button: only the primary confirming button; Cancel is never pressed
#'
#' @param modal_html The modal's rendered HTML.
#' @param button The button's id or label.
#' @return The pressed button's id.
#' @export
valModalPress <- function(modal_html, button) {
    doc <- xml2::read_html(paste0("<html><body>", paste(modal_html, collapse = "\n"),
                                  "</body></html>"))
    for (b in xml2::xml_find_all(doc, "//button")) {
        id <- xml2::xml_attr(b, "id")
        label <- trimws(xml2::xml_text(b))
        if (!identical(button, id) && !identical(button, label)) next
        dismiss <- !is.na(xml2::xml_attr(b, "data-dismiss")) ||
            !is.na(xml2::xml_attr(b, "data-bs-dismiss"))
        if (dismiss || is.na(id) || !valHasClass(b, "btn-primary")) {
            valAbort("only a modal's primary confirming button is pressed", "val_spec_invalid")
        }
        return(id)
    }
    valAbort("the modal has no button {button}", "val_spec_invalid")
}

#' Mask a message as the lanes do before hashing it (hexadecimal runs, numbers)
#' @noRd
valMaskMessage <- function(msg) {
    msg <- gsub("[0-9a-f]{10,}", "<hex>", msg)
    gsub("[0-9]+(\\.[0-9]+)?", "<n>", msg)
}

#' The outcome of one action-script step (A21)
#'
#' The declared state reached: `ok`. State and history unchanged with an error notification:
#' `refused` with the class (or `unclassed`) and the SHA-256 of the masked message. Anything else:
#' `did_not_complete (unexpected_state)`. A `diagnostic` step expects no state change.
#'
#' @param declared The declared state, or `diagnostic`.
#' @param before,after The workflow state before and after the action.
#' @param history_changed The history changed.
#' @param error_notification The error notification's text, or `NULL`.
#' @param class The refusal's class where the package exposes one.
#' @return A list: `outcome`, `class`, `message_sha256`, `reason`.
#' @export
valStepOutcome <- function(declared, before, after, history_changed, error_notification = NULL,
                           class = NULL) {
    refused <- function() {
        list(outcome = "refused", class = class %||% "unclassed",
             message_sha256 = digest::digest(valMaskMessage(error_notification), algo = "sha256",
                                             serialize = FALSE),
             reason = NA_character_)
    }
    if (identical(declared, "diagnostic")) {
        if (!is.null(error_notification)) return(refused())
        return(list(outcome = "ok", class = NA_character_, message_sha256 = NA_character_,
                    reason = NA_character_))
    }
    if (identical(after, declared)) {
        return(list(outcome = "ok", class = NA_character_, message_sha256 = NA_character_,
                    reason = NA_character_))
    }
    if (identical(after, before) && !isTRUE(history_changed) && !is.null(error_notification)) {
        return(refused())
    }
    list(outcome = "did_not_complete", class = NA_character_, message_sha256 = NA_character_,
         reason = "unexpected_state")
}

#' Run action-script steps in order with the step-outcome rule
#'
#' A refused `diagnostic` step is recorded and the cell continues; any other refusal or an
#' unexpected state ends the cell.
#'
#' @param steps List of steps: `id`, `declared`, `run()` returning `after`, `history_changed`,
#'   `error_notification`, `class`.
#' @param start The state before the first step.
#' @param on_step Called with each step id before it runs.
#' @return A data frame (step, outcome, class, message_sha256, reason) with attribute `exit`.
#' @export
valRunSteps <- function(steps, start, on_step = NULL) {
    state <- start
    rows <- list()
    exit <- "completed"
    for (s in steps) {
        if (!is.null(on_step)) on_step(s$id)
        r <- s$run()
        o <- valStepOutcome(s[["declared"]], state, r[["after"]], r[["history_changed"]],
                            r[["error_notification"]], r[["class"]])
        rows[[length(rows) + 1L]] <- data.frame(step = s$id, outcome = o$outcome, class = o$class,
                                                message_sha256 = o$message_sha256,
                                                reason = o$reason, stringsAsFactors = FALSE)
        if (o$outcome == "ok" && !identical(s$declared, "diagnostic")) state <- r$after
        if (o$outcome == "refused" && identical(s$declared, "diagnostic")) next
        if (o$outcome != "ok") {
            exit <- o$outcome
            break
        }
    }
    out <- do.call(rbind, rows)
    attr(out, "exit") <- exit
    out
}

#' Evaluate a cell's code with its forced options, scoped to the cell (KIT-163 N9)
#'
#' Records each forced option and `MultiScholaR.ruv_design_panel` as found at the cell's start; a
#' forced option already set at the start is `val_spec_invalid`; all are restored afterwards.
#'
#' @param forced Named list of forced options.
#' @param code The cell's code.
#' @param watch Further options to record.
#' @return A list: `value`, `options_at_start`.
#' @export
valWithForcedOptions <- function(forced, code, watch = "MultiScholaR.ruv_design_panel") {
    nm <- unique(c(names(forced), watch))
    at_start <- lapply(stats::setNames(nm, nm), function(o) {
        v <- getOption(o)
        if (is.null(v)) "unset" else paste(format(v), collapse = "|")
    })
    already <- names(forced)[!vapply(names(forced), function(o) is.null(getOption(o)), logical(1))]
    if (length(already) > 0L) {
        valAbort("a forced option is already set at the cell's start: {already}",
                 "val_spec_invalid")
    }
    value <- withr::with_options(forced, code)
    list(value = value, options_at_start = at_start)
}

#' Refuse an action script that presses an enrichment button (enrichment never runs under test mode)
#'
#' @param script List of steps with `press`.
#' @return `TRUE`.
#' @export
valCheckActionScript <- function(script) {
    presses <- unlist(lapply(script, function(s) s[["press"]]))
    if (any(grepl("enrich", presses, ignore.case = TRUE))) {
        valAbort("an action script presses an enrichment button", "val_enrichment_test_mode")
    }
    TRUE
}

#' The matrix stages a cell's standard tables hold (section 7.3; provisional ledger row H6)
#'
#' RUV off: the normalised matrix. RUV correcting the matrix: the corrected one. RUV factors in the
#' model (transcriptomics; DIA_LIMPA covariates; lipidomics and metabolomics where the design puts
#' the factors in the model, revision 1.8): `matrix.tsv` holds the uncorrected values the test
#' reads, `ruv.tsv` holds W, and `matrix_display.tsv` the corrected matrix, labelled display.
#'
#' @param omic The omic.
#' @param ruv `off`, `correction` or `factors_in_model`.
#' @return A list: `matrix_pre_ruv`, `matrix`, `display`, `display_label`, `ruv_w`.
#' @export
valMatrixStages <- function(omic, ruv = c("off", "correction", "factors_in_model")) {
    ruv <- match.arg(ruv)
    switch(ruv,
           off = list(matrix_pre_ruv = "normalised", matrix = "normalised", display = NULL,
                      display_label = NULL, ruv_w = FALSE),
           correction = list(matrix_pre_ruv = "normalised", matrix = "corrected", display = NULL,
                             display_label = NULL, ruv_w = TRUE),
           factors_in_model = list(matrix_pre_ruv = "normalised", matrix = "uncorrected",
                                   display = "corrected", display_label = MATRIX_DISPLAY_LABEL,
                                   ruv_w = TRUE))
}

#' Top-level function assignments of R files, with their source text and file
#' @noRd
valTopLevelFunctions <- function(files) {
    defs <- list()
    for (f in files) {
        for (e in parse(f, keep.source = FALSE)) {
            if (is.call(e) && as.character(e[[1L]]) %in% c("<-", "=") && is.name(e[[2L]]) &&
                    is.call(e[[3L]]) && identical(e[[3L]][[1L]], as.name("function"))) {
                nm <- as.character(e[[2L]])
                defs[[nm]] <- list(expr = e[[3L]], file = basename(f),
                                   text = paste(deparse(e), collapse = "\n"))
            }
        }
    }
    defs
}

#' Extract named helper functions and their closure, never running top-level code (G1(c), ENG-23)
#'
#' Only top-level `name <- function(...)` assignments are evaluated, into an environment whose
#' parent is `parent` (the MultiScholaR namespace in a job). The transitive closure over the files'
#' top-level functions is found with `codetools::findGlobals`; a function called that is defined
#' neither in the files nor reachable from `parent` is `val_driver_helpers_invalid`.
#'
#' @param files The helper files.
#' @param names The named functions.
#' @param parent The extraction environment's parent.
#' @return An environment with the functions, attribute `helpers` (file, function,
#'   sha256_of_source).
#' @export
valExtractHelpers <- function(files, names, parent = asNamespace("MultiScholaR")) {
    defs <- valTopLevelFunctions(files)
    env <- new.env(parent = parent)
    missing <- setdiff(names, base::names(defs))
    if (length(missing) > 0L) {
        valAbort("helper {missing} is not a top-level function of the files",
                 "val_driver_helpers_invalid")
    }
    need <- names
    done <- character()
    while (length(need) > 0L) {
        nm <- need[[1L]]
        need <- need[-1L]
        if (nm %in% done) next
        done <- c(done, nm)
        fun <- eval(defs[[nm]]$expr, envir = env)
        assign(nm, fun, envir = env)
        g <- codetools::findGlobals(fun, merge = FALSE)
        local_defs <- intersect(c(g$functions, g$variables), base::names(defs))
        need <- c(need, setdiff(local_defs, done))
        unresolved <- setdiff(g$functions, base::names(defs))
        unresolved <- unresolved[!vapply(unresolved, exists, logical(1), envir = parent)]
        if (length(unresolved) > 0L) {
            valAbort("helper {nm} calls {unresolved}, defined nowhere",
                     "val_driver_helpers_invalid")
        }
    }
    attr(env, "helpers") <- data.frame(
        file = vapply(done, function(n) defs[[n]]$file, character(1)),
        `function` = done,
        sha256_of_source = vapply(done, function(n) {
            digest::digest(defs[[n]]$text, algo = "sha256", serialize = FALSE)
        }, character(1)),
        check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL
    )
    env
}

# ---------------------------------------------------------------------------------------------
# The driver core: valRunWorkflow (section 7.3), the helpers' extraction from a checkout, the
# fixture cells, the standard tables (valReadArmTables) and replay (valReplay).
# ---------------------------------------------------------------------------------------------

#' The GUI-entry helpers of a checkout, extracted (G1(c)): the nine named functions and the step
#' builders, with their closure; `install_enrichment_doubles` is the one no-op (7.3)
#'
#' @param pkg_root The package checkout.
#' @return The extraction environment (attribute `helpers`).
#' @export
valGuiHelpers <- function(pkg_root) {
    files <- file.path(pkg_root, "tests", "testthat",
                       c("helper-gui-entry.R", "helper-workflow-omics-residency.R",
                         "helper-e2e-browser.R"))
    named <- c("guiEntryServerSpec", "guiEntryInstallBrowserDouble", "guiEntryInstallUniprotCache",
               "guiEntryDesignActions", "guiEntryActionIds", "guiEntryMask", "guiEntryDigest",
               "guiEntryRunLane", "guiEntryTrxSteps", "guiEntryOmicsSteps", "guiEntryProtSteps",
               "e2e_match_available_values", "guiEntryHistory")
    parent <- new.env(parent = asNamespace("MultiScholaR"))
    parent$install_enrichment_doubles <- function() invisible(NULL)
    valExtractHelpers(files, named, parent = parent)
}

#' A route's cell on its committed e2e fixture (tests/testdata/e2e/manifest.json): the script, the
#' files by role, the design, the contrast, and for the six-sample fixtures the GUI-entry lanes'
#' permissive QC settings as declared forced inputs (proteomics)
#'
#' @param route Route (lane) id.
#' @param pkg_root The package checkout.
#' @return A list: `route`, `omic`, `script`, `files`, `design`, `contrasts`, `assays`,
#'   `forced_inputs`, `uniprot_cache`.
#' @export
valFixtureCell <- function(route, pkg_root) {
    e2e <- file.path(pkg_root, "tests", "testdata", "e2e")
    manifest <- jsonlite::read_json(file.path(e2e, "manifest.json"), simplifyVector = FALSE)
    lane <- Filter(function(l) identical(l$lane_id, route), manifest$lanes)
    if (length(lane) != 1L) valAbort("no e2e fixture for route {route}", "val_spec_invalid")
    lane <- lane[[1L]]
    dir <- file.path(e2e, lane$fixture_dir)
    design <- utils::read.delim(file.path(dir, "design.tsv"), stringsAsFactors = FALSE,
                                check.names = FALSE)
    contrasts <- lapply(unlist(lane$expected_contrasts), function(x) strsplit(x, "_vs_")[[1L]])
    assays <- unlist(lane$assays)
    omic <- lane$omic_type
    files <- switch(omic,
                    proteomics = list(search_results = file.path(dir, lane$seed_file),
                                      fasta = file.path(dir, lane$fasta_file)),
                    transcriptomics = list(counts = file.path(dir, lane$seed_file)),
                    c(list(assay1 = file.path(dir, lane$seed_file)),
                      if (!is.null(lane$assay2_file)) {
                          list(assay2 = file.path(dir, lane$assay2_file))
                      }))
    import_values <- if (identical(lane$import_tool, "custom")) {
        list(`import-metabolite_id_col_custom` = "Feature.Name",
             `import-annotation_col_custom` = "",
             `import-sample_cols_pattern` = "^(WT|KO)_",
             `import-is_pattern` = "^ITSD|_d[0-9]+$")
    } else {
        list()
    }
    forced_inputs <- if (identical(omic, "proteomics")) {
        structure(list(use_strict_mode = FALSE, min_reps_per_group = 1, min_groups = 1,
                       intensity_cutoff_percentile = 0.1,
                       proteins_intensity_cutoff_percentile = 0.1, min_peptides_per_protein = 1,
                       min_peptidoforms_per_protein = 1, min_peptides_per_sample = 1,
                       proportion_missing_values = 1.0, parallel_cores = 1),
                  reason = "six-sample e2e fixture: the GUI-entry lanes' QC settings")
    } else {
        list()
    }
    script <- valActionScript(route, assays = assays, import_tool = lane$import_tool,
                              import_values = import_values, contrasts = contrasts)
    list(route = route, omic = omic, script = script, files = files, design = design,
         contrasts = contrasts, assays = assays, forced_inputs = forced_inputs,
         uniprot_cache = if (identical(omic, "proteomics")) {
             file.path(pkg_root, "tests", "testdata", "gui-entry", "uniprot",
                       "uniprot_annotations.RDS")
         })
}

#' Drop inputs from the state (a modal's inputs when it is removed)
#' @noRd
valDropInputs <- function(st, ids) {
    ids <- intersect(ids, names(st$values))
    if (!length(ids)) return(invisible(st))
    st$values[ids] <- NULL
    keep <- setdiff(names(st$source), ids)
    st$source <- st$source[keep]
    st$location <- st$location[keep]
    invisible(st)
}

#' A refusal's class, where the message names one (`<class>: ...`)
#' @noRd
valMessageClass <- function(msg) {
    if (is.null(msg)) return(NULL)
    m <- regmatches(msg, regexec("^\\s*([a-z][a-z0-9]*(_[a-z0-9]+)+):", msg))[[1L]]
    if (length(m)) m[[2L]] else NULL
}

#' The observed and declared tokens of a step's state (head state, or a tab status)
#' @noRd
valStateToken <- function(declared, before, after) {
    if (identical(declared, "diagnostic")) {
        return(list(declared = "diagnostic", before = before$head, after = after$head))
    }
    if (!is.null(declared[["head"]])) {
        allowed <- declared[["head"]]
        return(list(declared = if (after$head %in% allowed) after$head else allowed[[1L]],
                    before = before$head, after = after$head))
    }
    key <- names(declared[["tab"]])[[1L]]
    tok <- function(obs) paste0(key, "=", obs$tab_status[[key]] %||% "<none>")
    list(declared = paste0(key, "=", declared[["tab"]][[1L]]), before = tok(before),
         after = tok(after))
}

#' The ledger with the listed inputs' identifiers written as their SHA-256, one per line (5.2)
#'
#' @param ledger Ledger rows.
#' @param ids Input ids whose values are identifier lists.
#' @return The ledger.
#' @export
valMaskLedger <- function(ledger, ids) {
    hit <- ledger$input_id %in% ids
    ledger$value[hit] <- vapply(ledger$value[hit], function(v) {
        lines <- strsplit(v, "\n", fixed = TRUE)[[1L]]
        paste(vapply(lines, function(x) digest::digest(x, algo = "sha256", serialize = FALSE),
                     character(1), USE.NAMES = FALSE), collapse = "\n")
    }, character(1), USE.NAMES = FALSE)
    ledger
}

#' Extend the browser double (A01 point 3, A18): a choices-only update selects the first choice;
#' a modal's UI is kept until removeModal
#' @noRd
valInstallDoubleExtension <- function(log, env, notify = FALSE) {
    queue <- function(session, input_id, value, choices = NULL) {
        if (is.null(session)) return(invisible(NULL))
        id <- session$ns(input_id)
        if (!is.null(choices)) log$choices[[id]] <- as.character(unlist(choices, use.names = FALSE))
        if (is.null(value) && !is.null(choices) && length(log$choices[[id]])) {
            value <- log$choices[[id]][[1L]]
        }
        if (is.null(value)) return(invisible(NULL))
        log$updates[[length(log$updates) + 1L]] <- list(id = id, value = value)
        invisible(NULL)
    }
    dom <- function() shiny::getDefaultReactiveDomain()
    bindings <- list(
        updateSelectInput = function(session = dom(), inputId, label = NULL, choices = NULL,
                                     selected = NULL) {
            queue(session, inputId, selected, choices)
        },
        updateSelectizeInput = function(session = dom(), inputId, label = NULL, choices = NULL,
                                        selected = NULL, options = list(), server = FALSE) {
            queue(session, inputId, selected, choices)
        },
        updateRadioButtons = function(session = dom(), inputId, label = NULL, choices = NULL,
                                      selected = NULL, inline = FALSE, choiceNames = NULL,
                                      choiceValues = NULL) {
            queue(session, inputId, selected, choices %||% choiceValues)
        },
        showModal = function(ui, session = dom()) {
            log$modals <- c(log$modals, 1L)
            log$modal_html <- paste(as.character(htmltools::renderTags(ui)$html), collapse = "\n")
            invisible(NULL)
        },
        removeModal = function(session = dom()) {
            log$modal_html <- NULL
            invisible(NULL)
        }
    )
    if (notify) {
        bindings$showNotification <- function(ui, action = NULL, duration = 5, closeButton = TRUE,
                                              id = NULL, type = "default", session = dom()) {
            log$notes[[length(log$notes) + 1L]] <- list(type = type,
                                                        message = paste(as.character(ui),
                                                                        collapse = " "))
            invisible(id %||% "notification")
        }
        bindings$removeNotification <- function(...) invisible(NULL)
    }
    do.call(testthat::local_mocked_bindings, c(bindings, list(.package = "shiny", .env = env)))
    if (isNamespaceLoaded("MultiScholaR")) {
        imported <- intersect(names(bindings),
                              unlist(getNamespaceImports("MultiScholaR")$shiny, use.names = FALSE))
        if (length(imported)) {
            do.call(testthat::local_mocked_bindings,
                    c(bindings[imported], list(.package = "MultiScholaR", .env = env)))
        }
    }
    invisible(log)
}

#' Run one workflow cell: a route's action script on a module server under `shiny::testServer`
#' (section 7.3)
#'
#' The runtime frame of `guiEntryRunLane` (directories, test mode on for the file-input fallback
#' only, proxies at a closed port, the browser double, extended), with the thread variables at 1,
#' the cell's seed, the forced options and the arm's options scoped to the cell, and for
#' proteomics the declared UniProt cache. Before each action the full input vector goes to the
#' ledger. Inputs come only from the rendered UI (static, `uiOutput`s, modals), server updates, the
#' step's structural inputs, the declared forced inputs and the arm's overrides; anything else is
#' `val_spec_invalid`. Step outcomes as A21. With `tables_dir`, the standard tables are written
#' after the last step; `on_complete(ctx)` is then called inside the session (the post-cell hook's
#' place).
#'
#' @param script From [valActionScript()] (or a list with `route`, `omic`, `steps`, `contrasts`).
#' @param files Named list, file role to path.
#' @param design The design (sample, group, optional tech_rep_group and replicates).
#' @param arm An arm (a builder's value, or a list with `options` and `overrides`).
#' @param helpers From [valGuiHelpers()].
#' @param forced Forced options (named list).
#' @param forced_inputs Forced inputs by id or by id suffix; attribute `reason`.
#' @param seed The cell's seed.
#' @param uniprot_cache A UniProt cache file placed where the design step reads it.
#' @param annotation `list(expected, min_share, accessions, in_taxon)` for the annotation rule.
#' @param tables_dir Where the standard tables are written, or `NULL`.
#' @param ledger_path The ledger file (default in the work folder).
#' @param mask_ids Further input ids whose identifiers the ledger hashes.
#' @param roles Prepared identifier lists by role (supplied controls).
#' @param on_complete `function(ctx)` called inside the session after the tables.
#' @param pkg_root The package checkout.
#' @param work_dir The cell's work folder.
#' @param spec,observe A module spec and state reader (tests); default the route's.
#' @param project Set up project directories and the store (FALSE for a toy module).
#' @param technical_replicates The design has technical replicates (the calibrated design).
#' @param record_messages Keep a refusal's or error's masked text in `steps$message` (FALSE when any
#'   data input is `not_permitted`: then only its SHA-256, 5.2).
#' @return A list: `exit`, `reason`, `steps`, `ledger`, `ledger_path`, `ledger_sha256`,
#'   `actions`, `options_at_start`, `matched`, `fallback`, `tables`, `value`, `seconds`.
#' @export
valRunWorkflow <- function(script, files, design, arm = valArmShipped(), helpers = NULL,
                           forced = list(), forced_inputs = list(), seed = 20261107L,
                           uniprot_cache = NULL, annotation = NULL, tables_dir = NULL,
                           ledger_path = NULL, mask_ids = character(), roles = list(),
                           on_complete = NULL, pkg_root = ".",
                           work_dir = tempfile("val-cell-"), spec = NULL, observe = NULL,
                           project = TRUE, technical_replicates = FALSE,
                           record_messages = TRUE) {
    started <- Sys.time()
    valCheckActionScript(script$steps)
    arm_inputs <- if (is.null(arm$builder)) {
        list(options = arm$options %||% list(), overrides = arm$overrides %||% list(),
             mask = character(), expect = list())
    } else {
        valArmInputs(arm, script$route, roles = roles,
                     technical_replicates = technical_replicates)
    }
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
    ledger_path <- ledger_path %||% file.path(work_dir, paste0("ledger-", script$route, ".tsv"))
    cell <- function() {
        valWorkflowCell(script, files, design, arm_inputs, helpers, forced_inputs, seed,
                        uniprot_cache, annotation, tables_dir, arm, on_complete, work_dir, spec,
                        observe, project, record_messages)
    }
    run <- valWithForcedOptions(c(forced, arm_inputs$options), cell())
    res <- run$value
    ledger <- if (length(res$ledger)) do.call(rbind, res$ledger) else
        valLedgerSnapshot(valInputState(), "none", "none")
    ledger <- valMaskLedger(ledger, c(arm_inputs$mask, mask_ids))
    dir.create(dirname(ledger_path), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(ledger, ledger_path, sep = "\t", quote = TRUE, row.names = FALSE)
    actions <- if (length(res$actions)) do.call(rbind, res$actions) else
        data.frame(step = character(), kind = character(), input_id = character(),
                   value = character(), stringsAsFactors = FALSE)
    list(exit = res$exit, reason = res$reason, steps = res$steps, ledger = ledger,
         ledger_path = ledger_path,
         ledger_sha256 = digest::digest(file = ledger_path, algo = "sha256"),
         actions = actions, options_at_start = run$options_at_start, matched = res$matched,
         fallback = res$fallback, tables = res$tables, value = res$value,
         seconds = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1))
}

#' The cell itself, inside the forced options (internal to valRunWorkflow)
#' @noRd
valWorkflowCell <- function(script, files, design, arm_inputs, helpers, forced_inputs, seed,
                            uniprot_cache, annotation, tables_dir, arm, on_complete, work_dir,
                            spec, observe, project, record_messages) {
    cell_env <- environment()
    withr::local_envvar(c(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
                          MULTISCHOLAR_TEST_MODE = "true", http_proxy = "http://127.0.0.1:9",
                          https_proxy = "http://127.0.0.1:9", HTTP_PROXY = "http://127.0.0.1:9",
                          HTTPS_PROXY = "http://127.0.0.1:9", no_proxy = "", NO_PROXY = ""))
    withr::local_options(list(warn = 1))
    withr::local_seed(seed)
    omic <- script$omic
    ns <- if (project) asNamespace("MultiScholaR") else NULL
    spec <- spec %||% helpers$guiEntryServerSpec(omic)
    observe <- observe %||% function(env) {
        wd <- env$workflow_data
        m <- shiny::isolate(wd$state_manager)
        cur_head <- tryCatch(shiny::isolate(m$current_state), error = function(e) NULL)
        list(head = cur_head %||% "<none>",
             history = helpers$guiEntryHistory(m),
             tab_status = lapply(shiny::isolate(wd$tab_status), as.character))
    }
    paths <- NULL
    project_dirs <- NULL
    if (project) {
        if (!ns$workflowOwnerCapable()) {
            valAbort("the session owner is not available (qs2/duckdb)", "val_spec_invalid")
        }
        project_dirs <- suppressMessages(ns$setupDirectories(
            base_dir = file.path(work_dir, "project", script$route), omic_types = omic,  # nolint
            label = NULL, force = TRUE, reuse_existing = TRUE))
        paths <- project_dirs[[omic]]
        old_dirs <- get0("project_dirs", envir = globalenv(), inherits = FALSE)
        assign("project_dirs", project_dirs, envir = globalenv())
        withr::defer(if (is.null(old_dirs)) rm("project_dirs", envir = globalenv()) else
            assign("project_dirs", old_dirs, envir = globalenv()), envir = cell_env)  # nolint
        if (identical(omic, "proteomics") && !is.null(uniprot_cache)) {
            target <- file.path(paths$results_dir, "cache", "uniprot_annotations")
            dir.create(target, recursive = TRUE, showWarnings = FALSE)
            file.copy(uniprot_cache, file.path(target, "uniprot_annotations.RDS"), overwrite = TRUE)
        }
    }
    log <- new.env(parent = emptyenv())
    log$updates <- list()
    log$notes <- list()
    log$modals <- integer()
    log$choices <- list()
    log$modal_html <- NULL
    log$modal_seen <- NULL
    log$rendered <- list()
    log$warnings <- character()
    if (!is.null(helpers)) helpers$guiEntryInstallBrowserDouble(log, env = cell_env)
    valInstallDoubleExtension(log, cell_env, notify = is.null(helpers))
    scratch <- file.path(work_dir, "inputs")
    dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
    st <- valInputState()
    acc <- new.env(parent = emptyenv())
    acc$ledger <- list()
    acc$actions <- list()
    acc$steps <- list()
    acc$exit <- "completed"
    acc$reason <- NA_character_
    acc$matched <- NULL
    acc$fallback <- NA_character_
    acc$tables <- NULL
    acc$value <- NULL
    acc$destroyed <- FALSE
    server_args <- if (project) {
        list(id = spec$id, project_dirs = project_dirs, omic_type = omic,
             experiment_label = script$route)
    } else {
        c(list(id = spec$id), spec$args)
    }
    body <- function(session, workflow_data) {
        rootSession <- session$rootScope()
        prefix <- paste0(spec$id, "-")
        obs_env <- new.env(parent = emptyenv())
        obs_env$workflow_data <- workflow_data
        clicks <- list()
        current <- NULL
        rel <- function(ids) {
            if (!length(ids)) return(character())
            ifelse(startsWith(ids, prefix), substring(ids, nchar(prefix) + 1L), ids)
        }
        setFull <- function(values) {
            if (!length(values)) return(invisible())
            names(values) <- paste0(prefix, names(values))
            do.call(rootSession$setInputs, values)
        }
        ledgerBefore <- function(action) {
            acc$ledger[[length(acc$ledger) + 1L]] <- valLedgerSnapshot(st, current$id, action)
        }
        logAction <- function(kind, id, value) {
            acc$actions[[length(acc$actions) + 1L]] <- data.frame(
                step = current$id, kind = kind, input_id = id,  # nolint
                value = paste(format(value), collapse = "|"), stringsAsFactors = FALSE)
        }
        applyUpdates <- function() {
            changed <- FALSE
            for (i in seq_len(50L)) {
                if (!length(log$updates)) break
                pending <- log$updates
                log$updates <- list()
                values <- stats::setNames(lapply(pending, `[[`, "value"),
                                          rel(vapply(pending, `[[`, "", "id")))
                values <- values[!duplicated(names(values), fromLast = TRUE)]
                for (id in names(values)) valApplyServerUpdate(st, id, value = values[[id]])
                setFull(values)
                changed <- TRUE
            }
            changed
        }
        scanRendered <- function(module) {
            private <- rootSession$.__enclos_env__$private
            outs <- names(private$outs)
            scope <- if (is.null(module)) prefix else paste0(prefix, module, "-")
            changed <- FALSE
            for (name in outs[startsWith(outs, scope)]) {
                of <- tryCatch(attr(private$outs[[name]]$func, "outputFunc"),
                               error = function(e) NULL)
                if (!is.null(of) && !identical(of, shiny::uiOutput) &&
                        !identical(of, shiny::htmlOutput)) next
                value <- tryCatch(rootSession$getOutput(name), error = function(e) NULL)
                if (!is.list(value) || is.null(value$html)) next
                html <- paste(as.character(value$html), collapse = "\n")
                d <- digest::digest(html, algo = "sha256", serialize = FALSE)
                if (identical(log$rendered[[name]], d)) next
                rerender <- !is.null(log$rendered[[name]])
                log$rendered[[name]] <- d
                defaults <- valUiDefaults(html)
                if (!length(defaults)) next
                names(defaults) <- rel(names(defaults))
                before <- names(st$values)
                valApplyRendered(st, defaults, location = paste0("uiOutput:", rel(name)),
                                 rerender = rerender)
                new <- if (rerender) names(defaults) else setdiff(names(defaults), before)
                if (length(new)) {
                    setFull(st$values[new])
                    changed <- TRUE
                }
            }
            changed
        }
        applyModal <- function() {
            if (identical(log$modal_seen, log$modal_html)) return(FALSE)
            if (is.null(log$modal_html)) {
                valDropInputs(st, names(st$location)[st$location == "modal"])
            } else {
                defaults <- valUiDefaults(log$modal_html)
                names(defaults) <- rel(names(defaults))
                valApplyRendered(st, defaults, location = "modal", rerender = TRUE)
                setFull(st$values[names(defaults)])
            }
            log$modal_seen <- log$modal_html
            TRUE
        }
        applyForced <- function() {
            changed <- FALSE
            for (name in names(forced_inputs)) {
                ids <- names(st$values)
                hit <- ids[ids == name | endsWith(ids, paste0("-", name))]
                for (id in hit) {
                    if (identical(unname(st$source[id]), "forced") &&
                            identical(st$values[[id]], forced_inputs[[name]])) next
                    valSetInput(st, id, forced_inputs[[name]], "forced")
                    setFull(stats::setNames(list(forced_inputs[[name]]), id))
                    changed <- TRUE
                }
            }
            changed
        }
        settle <- function() {
            for (i in seq_len(10L)) {
                session$flushReact()
                updated <- applyUpdates()
                modal <- applyModal()
                rendered <- scanRendered(current$module)
                forced_now <- applyForced()
                if (!(updated || modal || rendered || forced_now)) break
            }
            invisible()
        }
        structuralId <- function(id) {
            any(vapply(current$structural, function(p) grepl(p, id), logical(1)))
        }
        setValues <- function(values) {
            for (id in names(values)) {
                if (!structuralId(id)) {
                    valAbort(paste("the action script sets {id}, which is not one of the",
                                   "step's structural inputs"), "val_spec_invalid")
                }
            }
            ledgerBefore(paste(names(values), collapse = ","))
            for (id in names(values)) {
                valSetInput(st, id, values[[id]], "structural")
                logAction("set", id, values[[id]])
            }
            setFull(values)
            settle()
        }
        pressValue <- function(id, value) {
            ledgerBefore(id)
            logAction("press", id, value)
            clicks[[id]] <<- value  # nolint
            setFull(stats::setNames(list(value), id))
            settle()
        }
        press <- function(id) pressValue(id, (clicks[[id]] %||% 0L) + 1L)
        api <- list(
            helpers = helpers,
            defaults = function(ui_fn, module = NULL) {
                html <- paste(as.character(htmltools::renderTags(ui_fn(module))$html),
                              collapse = "\n")
                defaults <- valUiDefaults(html)
                before <- names(st$values)
                valApplyRendered(st, defaults, location = "ui")
                new <- setdiff(names(defaults), before)
                if (length(new)) setFull(st$values[new])
                settle()
            },
            set = setValues,
            press = press,
            upload = function(id, role) {
                path <- files[[role]]
                if (is.null(path) || !file.exists(path)) {
                    valAbort("no declared file for role {role}", "val_spec_invalid")
                }
                if (!structuralId(id)) {
                    valAbort("{id} is not a structural file input of the step", "val_spec_invalid")
                }
                copy <- file.path(scratch, basename(path))
                file.copy(path, copy, overwrite = TRUE)
                ledgerBefore(id)
                logAction("upload", id, path)
                valSetInput(st, id, basename(path), "structural")
                setFull(stats::setNames(list(data.frame(name = basename(path),
                                                        size = file.size(copy), type = "",
                                                        datapath = copy,
                                                        stringsAsFactors = FALSE)), id))
                settle()
            },
            has_file = function(role) !is.null(files[[role]]),
            choices = function(id) log$choices[[paste0(prefix, id)]] %||% character(),
            modal_open = function() !is.null(log$modal_html),
            modal_confirm = function(button = NULL) {
                html <- log$modal_html
                if (is.null(html)) valAbort("no modal is open", "val_spec_invalid")
                id <- if (is.null(button)) valModalPrimary(html) else valModalPress(html, button)
                press(rel(id))
            },
            overrides = function() {
                for (id in names(arm_inputs$overrides)) {
                    if (!id %in% names(st$values)) {
                        valAbort("the arm overrides {id}, which the route does not render here",
                                 "val_spec_invalid")
                    }
                }
                for (id in names(arm_inputs$expect)) {
                    if (!identical(st$values[[id]], arm_inputs$expect[[id]])) {
                        valAbort("{id} is not the arm's expected value", "val_spec_invalid")
                    }
                }
                if (!length(arm_inputs$overrides)) return(invisible())
                ledgerBefore(paste(names(arm_inputs$overrides), collapse = ","))
                for (id in names(arm_inputs$overrides)) {
                    valSetInput(st, id, arm_inputs$overrides[[id]], "override")
                    logAction("set", id, arm_inputs$overrides[[id]])
                }
                setFull(arm_inputs$overrides)
                settle()
            },
            design = function(prefix_builder, tech_rep = TRUE, match = NULL) {
                contrasts <- script$contrasts
                if (!length(contrasts)) {
                    valAbort("the cell declares no contrast", "val_spec_invalid")
                }
                d <- design
                if (!is.null(match)) {
                    offered <- log$choices[[paste0(prefix, match)]] %||% character()
                    matched <- helpers$e2e_match_available_values(d$sample, offered)
                    acc$matched <- data.frame(declared = d$sample, matched = matched,
                                              stringsAsFactors = FALSE)
                    d$sample <- matched
                }
                extra <- contrasts[-1L]
                buttons <- paste0("(add_factor|assign_metadata|assign_tech_reps|add_contrast|",
                                  "save_results)$")
                dapi <- list(set = function(values) {
                    ids <- names(values)
                    formula <- paste0(prefix_builder, "formula_string")
                    if (formula %in% ids && formula %in% names(st$values)) {
                        values <- values[ids != formula]
                        ids <- names(values)
                    }
                    if (any(grepl("save_results$", ids)) && length(extra)) {
                        n <- clicks[[paste0(prefix_builder, "add_contrast")]] %||% 0L
                        for (pair in extra) {
                            pair_ids <- paste0(prefix_builder,
                                               c("contrast_group1", "contrast_group2"))
                            setValues(stats::setNames(list(pair[[1L]], pair[[2L]]), pair_ids))
                            n <- n + 1L
                            pressValue(paste0(prefix_builder, "add_contrast"), n)
                        }
                        extra <<- list()
                    }
                    is_button <- grepl(buttons, ids)
                    if (any(!is_button)) setValues(values[!is_button])
                    for (id in ids[is_button]) pressValue(id, as.integer(values[[id]]))
                })
                helpers$guiEntryDesignActions(dapi, d, contrasts[[1L]], prefix = prefix_builder,
                                              tech_rep = tech_rep)
            }
        )
        on_warning <- function(w) {
            log$warnings <- c(log$warnings, class(w)[[1L]])
            invokeRestart("muffleWarning")
        }
        for (step in script$steps) {
            current <- step
            before <- observe(obs_env)
            log$notes <- list()
            log$warnings <- character()
            err <- NULL
            tryCatch(withCallingHandlers({
                if (!is.null(step$tab)) {
                    tab <- stats::setNames(list(step$tab), spec$tabs)
                    ledgerBefore(spec$tabs)
                    valSetInput(st, spec$tabs, step$tab, "structural")
                    logAction("set", spec$tabs, step$tab)
                    setFull(tab)
                    settle()
                }
                step$run(api)
            }, warning = on_warning), error = function(e) {
                if (inherits(e, "val_error")) stop(e)
                err <<- e
            })
            after <- tryCatch(observe(obs_env), error = function(e) {
                acc$destroyed <- TRUE
                list(head = "<session destroyed>", history = before$history,
                     tab_status = before$tab_status)
            })
            errors <- Filter(function(n) identical(n$type, "error"), log$notes)
            msg <- if (length(errors)) {
                paste(vapply(errors, function(n) n$message, character(1)), collapse = "\n")
            }
            if (is.null(err) && isTRUE(acc$destroyed)) {
                # An R error in a server observer ends the session (testServer does not rethrow it):
                # the cell did not complete on an error.
                err <- structure(class = c("val_session_ended", "error", "condition"),
                                 list(message = paste("the module session ended on an error",
                                                      "in the server"),
                                      call = NULL))
            }
            if (!is.null(err)) {
                acc$steps[[length(acc$steps) + 1L]] <- data.frame(
                    step = step$id, outcome = "did_not_complete", class = class(err)[[1L]],  # nolint
                    message_sha256 = digest::digest(valMaskMessage(conditionMessage(err)),
                                                    algo = "sha256", serialize = FALSE),
                    reason = "error", head = after$head,  # nolint
                    message = valStepMessage(conditionMessage(err), record_messages),
                    stringsAsFactors = FALSE)
                acc$exit <- "did_not_complete"
                acc$reason <- "error"
                break
            }
            tok <- valStateToken(step$declared, before, after)
            o <- valStepOutcome(tok$declared, tok$before, tok$after,
                                !identical(before$history, after$history), msg,
                                valMessageClass(msg))
            acc$steps[[length(acc$steps) + 1L]] <- data.frame(
                step = step$id, outcome = o$outcome, class = o$class,  # nolint
                message_sha256 = o$message_sha256, reason = o$reason, head = after$head,
                message = valStepMessage(msg, record_messages), stringsAsFactors = FALSE)
            if (o$outcome == "ok" && identical(step$id, "design") &&
                    identical(omic, "proteomics")) {
                acc$fallback <- valAnnotationRule(shiny::isolate(workflow_data$uniprot_dat_cln),
                                                  annotation, log$warnings)
            }
            if (o$outcome == "ok") next
            if (o$outcome == "refused" && identical(step$declared, "diagnostic")) next
            acc$exit <- o$outcome
            acc$reason <- o$reason
            break
        }
        if (isTRUE(acc$destroyed)) return(invisible())
        ctx <- list(manager = if (project) shiny::isolate(workflow_data$state_manager),
                    workflow_data = workflow_data, paths = paths, inputs = st,
                    observe = function() observe(obs_env), session = session)
        if (identical(acc$exit, "completed") && !is.null(tables_dir)) {
            num <- function(suffix, default) {
                id <- names(st$values)[endsWith(names(st$values), suffix)]
                v <- if (length(id)) suppressWarnings(as.numeric(st$values[[id[[1L]]]])) else NA
                if (length(v) != 1L || is.na(v)) default else v
            }
            itsd <- names(st$values)[endsWith(names(st$values), "-apply_itsd")]
            da_head <- acc$steps[[length(acc$steps)]]$head
            itsd_applied <- if (length(itsd)) isTRUE(st$values[[itsd[[1L]]]]) else NA
            acc$tables <- valReadArmTables(ctx$manager, paths, script$route, tables_dir,
                                           ruv = valCellRuv(script$route, arm),
                                           q_thresh = num("-da_q_val_thresh", 0.05),
                                           lfc_cutoff = num("-treat_lfc_cutoff", 0),
                                           model_state = da_head,
                                           itsd_applied = itsd_applied)
            ctx$tables <- acc$tables
        }
        if (!is.null(on_complete)) acc$value <- on_complete(ctx)
        if (project) shiny::isolate(ns$closeDIAStateAdapter(workflow_data$state_manager))
    }
    shiny::testServer(spec$server, args = server_args, {
        body(session, workflow_data)  # nolint
    })
    list(exit = acc$exit, reason = acc$reason,
         steps = if (length(acc$steps)) do.call(rbind, acc$steps) else NULL,
         ledger = acc$ledger, actions = acc$actions, matched = acc$matched,
         fallback = acc$fallback, tables = acc$tables, value = acc$value)
}

#' A step's masked message text for `steps$message` (at most 300 characters), or NA
#' @noRd
valStepMessage <- function(msg, record) {
    if (!isTRUE(record) || is.null(msg) || !length(msg)) return(NA_character_)
    substr(gsub("[\t\r\n]+", " ", valMaskMessage(paste(msg, collapse = " "))), 1L, 300L)
}

#' The one primary confirming button of a modal (no dismiss), by id
#' @noRd
valModalPrimary <- function(modal_html) {
    doc <- xml2::read_html(paste0("<html><body>", modal_html, "</body></html>"))
    ids <- character()
    for (b in xml2::xml_find_all(doc, "//button")) {
        dismiss <- !is.na(xml2::xml_attr(b, "data-dismiss")) ||
            !is.na(xml2::xml_attr(b, "data-bs-dismiss"))
        id <- xml2::xml_attr(b, "id")
        if (!dismiss && !is.na(id) && valHasClass(b, "btn-primary")) ids <- c(ids, id)
    }
    if (length(ids) != 1L) {
        valAbort("the modal has no single primary confirming button", "val_spec_invalid")
    }
    ids
}

#' The annotation rule after the design step (7.3; 5.1, G6, C72)
#'
#' @param uniprot The bus's `uniprot_dat_cln`.
#' @param annotation `list(expected, min_share, accessions, in_taxon)` or `NULL` (no rule).
#' @param warning_classes Warning classes raised during the design step.
#' @return The cell's fallback flag (`NA` when none).
#' @noRd
valAnnotationRule <- function(uniprot, annotation, warning_classes) {
    if (is.null(annotation)) return(NA_character_)
    if (identical(annotation$expected, "none")) {
        return(paste(c("annotation_none", intersect(warning_classes,
                                                    c("prot_annotation_uniprot_incomplete",
                                                      "prot_annotation_uniprot_empty"))),
                     collapse = ":"))
    }
    if ("prot_annotation_uniprot_incomplete" %in% warning_classes) {
        valAbort("the UniProt answer is incomplete", "val_annotation_incomplete")
    }
    if (is.null(uniprot) || !NROW(uniprot)) {
        valAbort("the annotation table has no rows", "val_annotation_empty")
    }
    acc <- unique(annotation$accessions)
    in_taxon <- acc[annotation$in_taxon(acc)]
    share <- if (length(in_taxon)) mean(in_taxon %in% uniprot$Entry) else 0
    if (share < annotation$min_share) {
        valAbort("the annotation covers too few of the taxon's accessions", "val_annotation_empty")
    }
    NA_character_
}

# ---------------------------------------------------------------------------------------------
# Standard tables (7.3) and replay (7.3, G5).
# ---------------------------------------------------------------------------------------------

#' How a route's DA table is read (the product's own written table; G3 re-reads the maps):
#' proteomics `R/func_prot_da_results_io.R:332`, lipidomics `R/func_lipid_da_export.R:164`,
#' metabolomics `R/func_metab_da_results_output.R:166`, transcriptomics
#' `R/utils_workflow_trx_declarations.R:903` (edgeR `topTags`, BH)
#' @noRd
valDaMap <- function(omic) {
    switch(omic,
           proteomics = list(pattern = "^da_proteins_(.+)_long_annot\\.tsv$", sep = "\t",
                             feature = character(), assay = NULL, log2fc = "log2FC",
                             p = "raw_pvalue", s_own = "fdr_qvalue", contrast = "comparison",
                             mean = "^log2norm\\.", call_statistic = "qvalue"),
           lipidomics = list(pattern = "^de_(.+)_lipids_(.+)_long_annot\\.tsv$", sep = "\t",
                             feature = "lipid_id", assay = "assay", log2fc = "logFC",
                             p = "raw_pvalue", s_own = "fdr_qvalue", contrast = "comparison",
                             mean = "^intensity\\.", call_statistic = "qvalue"),
           metabolomics = list(pattern = "^de_(.+)_metabolites_(.+)_long_annot\\.tsv$",
                               sep = "\t", feature = "metabolite_id", assay = "assay",
                               log2fc = "logFC", p = "raw_pvalue", s_own = "fdr_qvalue",
                               contrast = "comparison", mean = "^intensity\\.",
                               call_statistic = "qvalue"),
           transcriptomics = list(pattern = "^de_results_(.+)\\.csv$", sep = ",",
                                  feature = character(), assay = NULL, log2fc = "logFC",
                                  p = "PValue", s_own = "FDR", contrast = NULL, mean = "logCPM",
                                  call_statistic = "BH"))
}

#' A slot of an S4 object or an element of a list
#' @noRd
valField <- function(obj, name) {
    if (isS4(obj)) {
        if (name %in% methods::slotNames(obj)) methods::slot(obj, name) else NULL
    } else {
        obj[[name]]
    }
}

#' A feature table as `feature_id, assay, <samples>` (numeric sample columns only)
#' @noRd
valFeatureTable <- function(tbl, id_col, assay) {
    tbl <- as.data.frame(tbl, stringsAsFactors = FALSE)
    samples <- setdiff(names(tbl)[vapply(tbl, is.numeric, logical(1))], id_col)
    out <- data.frame(feature_id = as.character(tbl[[id_col]]), assay = rep(assay, nrow(tbl)),
                      stringsAsFactors = FALSE)
    cbind(out, tbl[, samples, drop = FALSE])
}

#' Bind per-assay feature tables, the union of sample columns
#' @noRd
valBindAssays <- function(tables) {
    samples <- unique(unlist(lapply(tables, function(t) {
        setdiff(names(t), c("feature_id", "assay"))
    })))
    do.call(rbind, lapply(tables, function(t) {
        for (s in setdiff(samples, names(t))) t[[s]] <- NA_real_
        t[, c("feature_id", "assay", samples), drop = FALSE]
    }))
}

#' A state object's matrix by stage (`data`: the object's values; `uncorrected`: the values the
#' test reads where RUV factors enter the model, `args$ruv_uncorrected_assays`)
#' @noRd
valStateMatrix <- function(obj, omic, which = c("data", "uncorrected")) {
    which <- match.arg(which)
    if (is.null(obj)) valAbort("the state object is missing", "val_spec_invalid")
    switch(omic,
           proteomics = valFeatureTable(valField(obj, "protein_quant_table"),
                                        valField(obj, "protein_id_column"), "all"),
           transcriptomics = {
               tbl <- valField(obj, "count_data")[[1L]]
               id <- valField(obj, "gene_id_column")
               counts <- valFeatureTable(tbl, id, "all")
               m <- as.matrix(counts[, -(1:2), drop = FALSE])
               # The route's own transform (R/utils_workflow_trx_declarations.R:250-251, 278): TMM,
               # log-CPM, prior count 2.
               dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = m), method = "TMM")
               logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
               cbind(counts[, 1:2], as.data.frame(logcpm, stringsAsFactors = FALSE))
           },
           {
               slot <- if (identical(omic, "lipidomics")) "lipid_data" else "metabolite_data"
               id <- valField(obj, if (identical(omic, "lipidomics")) "lipid_id_column" else
                   "metabolite_id_column")  # nolint
               data <- if (identical(which, "uncorrected")) {
                   valField(obj, "args")$ruv_uncorrected_assays
               } else {
                   valField(obj, slot)
               }
               if (is.null(data)) valAbort("the state holds no {which} values", "val_spec_invalid")
               valBindAssays(lapply(names(data), function(a) valFeatureTable(data[[a]], id, a)))
           })
}

#' The state that enters RUV, per route (7.3's table)
#' @noRd
valPreRuvState <- function(route) {
    r <- valRoute(route)
    switch(r$omic,
           proteomics = if (identical(r$workflow, "DIA_LIMPA")) "peptide_count_filtered" else
               "normalized",
           lipidomics = "lipid_normalized",
           metabolomics = "metab_normalized",
           transcriptomics = "trx_qc_complete")
}

#' Read a route's DA tables into the standard `da.tsv` rows
#' @noRd
valReadDa <- function(da_dir, omic, q_thresh, lfc_cutoff) {
    map <- valDaMap(omic)
    files <- list.files(da_dir, pattern = map$pattern)
    if (!length(files)) {
        valAbort("no DA table of the route in its DA folder", "val_spec_invalid")
    }
    rows <- lapply(sort(files), function(f) {
        path <- file.path(da_dir, f)
        x <- if (identical(map$sep, ",")) {
            utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
        } else {
            utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
        }
        feature <- intersect(map$feature, names(x))
        feature <- if (length(feature)) feature[[1L]] else names(x)[[1L]]
        needed <- c(map$log2fc, map$p, map$s_own)
        if (!all(needed %in% names(x))) {
            valAbort("a DA table lacks {setdiff(needed, names(x))}", "val_spec_invalid")
        }
        contrast <- if (is.null(map$contrast)) {
            rep(sub(map$pattern, "\\1", f), nrow(x))
        } else {
            as.character(x[[map$contrast]])
        }
        mean_expression <- if (map$mean %in% names(x)) {
            as.numeric(x[[map$mean]])
        } else {
            cols <- grep(map$mean, names(x), value = TRUE)
            v <- if (length(cols)) rowMeans(as.matrix(x[, cols, drop = FALSE]), na.rm = TRUE) else
                rep(NA_real_, nrow(x))
            ifelse(is.nan(v), NA_real_, v)
        }
        log2fc <- as.numeric(x[[map$log2fc]])
        s_own <- as.numeric(x[[map$s_own]])
        data.frame(
            feature_id = as.character(x[[feature]]),  # nolint
            assay = if (is.null(map$assay)) rep("all", nrow(x)) else as.character(x[[map$assay]]),
            contrast = contrast, log2fc = log2fc, p = as.numeric(x[[map$p]]), s_own = s_own,
            call_statistic = rep(map$call_statistic, nrow(x)),
            called_own = !is.na(s_own) & s_own < q_thresh & !is.na(log2fc) &
                abs(log2fc) >= lfc_cutoff,
            mean_expression = mean_expression, stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}

#' The `ruv.tsv` rows of a state object: per assay the design, k, controls, features tested,
#' fallback and W (long), from `args$ruv_factors` where the state holds it
#' @noRd
valRuvRows <- function(obj, design_name) {
    factors <- valField(obj, "args")$ruv_factors
    cols <- c("assay", "design", "k", "controls_n", "controls_share", "features_tested",
              "fallback", "sample", "factor", "w")
    if (!length(factors)) {
        return(stats::setNames(data.frame("all", design_name %||% NA_character_, NA_integer_,
                                          NA_integer_, NA_real_, NA_integer_, "not_recorded",
                                          NA_character_, NA_character_, NA_real_,
                                          stringsAsFactors = FALSE), cols))
    }
    do.call(rbind, lapply(names(factors), function(a) {
        f <- factors[[a]]
        w <- f$W
        tested <- f$features_tested %||% NA_integer_
        base <- list(a, f$design %||% design_name %||% NA_character_,
                     as.integer(f$k %||% if (is.matrix(w)) ncol(w) else NA_integer_),
                     as.integer(f$controls_n %||% NA_integer_),
                     if (is.na(tested)) NA_real_ else (f$controls_n %||% NA_integer_) / tested,
                     as.integer(tested), f$fallback %||% "none")
        if (!is.matrix(w) || !length(w)) {
            return(stats::setNames(data.frame(base, NA_character_, NA_character_, NA_real_,
                                              stringsAsFactors = FALSE), cols))
        }
        grid <- expand.grid(sample = rownames(w) %||% as.character(seq_len(nrow(w))),
                            factor = colnames(w) %||% paste0("W", seq_len(ncol(w))),
                            stringsAsFactors = FALSE)
        out <- data.frame(base, stringsAsFactors = FALSE)[rep(1L, nrow(grid)), , drop = FALSE]
        out <- cbind(out, grid, as.vector(w))
        names(out) <- cols
        rownames(out) <- NULL
        out
    }))
}

#' Write a cell's standard tables (section 7.3; the stage rule A12, revision 1.8 item 4)
#'
#' `da.tsv` (`feature_id, assay, contrast, log2fc, p, s_own, call_statistic, called_own,
#' mean_expression`) from the route's own DA tables; `matrix_pre_ruv.tsv` (the normalised matrix
#' that enters RUV) and `matrix.tsv` (the matrix the DA model is fitted to: the state DA ran from;
#' where RUV factors enter the model, the uncorrected values) with an `assay` column;
#' `features.tsv`; `ruv.tsv` when RUV ran; `itsd.tsv` on lipidomics and metabolomics; and where RUV
#' factors enter the model `matrix_display.tsv`, labelled "display; not what is tested".
#'
#' @param manager The cell's state manager (used through `state_reader`).
#' @param paths The experiment paths (`da_output_dir`).
#' @param route Route id.
#' @param out_dir Where the tables are written (the cell's work folder).
#' @param ruv `off`, `correction` or `factors_in_model` ([valCellRuv()]).
#' @param q_thresh,lfc_cutoff The DA call rule's inputs (the ledger's values).
#' @param state_reader `function(state)` returning a committed state's object.
#' @param model_state The state DA ran from.
#' @param itsd_applied Whether ITSD normalisation was applied (lipidomics, metabolomics).
#' @param design_name The arm's RUV design, for `ruv.tsv`.
#' @param display_reader `function()` returning the display matrix (`feature_id, assay, <samples>`)
#'   or `NULL`; default, with a manager, the product's own corrected matrix on DIA_LIMPA and
#'   transcriptomics ([valProductDisplayReader()]).
#' @return Named character vector of the written paths.
#' @export
valReadArmTables <- function(manager, paths, route, out_dir, ruv = "off", q_thresh = 0.05,
                             lfc_cutoff = 0, state_reader = NULL, model_state = NULL,
                             itsd_applied = NA, design_name = NULL, display_reader = NULL) {
    r <- valRoute(route)
    omic <- r$omic
    if (is.null(display_reader) && !is.null(manager)) {
        display_reader <- valProductDisplayReader(manager, r)
    }
    state_reader <- state_reader %||% function(name) {
        ns <- asNamespace("MultiScholaR")
        ns$withWorkflowStateScope(manager, function() manager$getState(name))
    }
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    write <- function(x, name) {
        path <- file.path(out_dir, paste0(name, ".tsv"))
        utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
        path
    }
    da <- valReadDa(paths$da_output_dir, omic, q_thresh, lfc_cutoff)
    pre_obj <- state_reader(valPreRuvState(route))
    model_obj <- if (identical(omic, "transcriptomics")) pre_obj else
        state_reader(model_state %||% valPreRuvState(route))
    # With RUV applied on DIA, LFQ and TMT the product commits no `normalized` state; the matrix
    # that entered RUV is the RUV state's own record of it (S0 run 45277101451).
    pre <- if (is.null(pre_obj) && identical(omic, "proteomics") && !is.null(manager)) {
        valProtRuvInput(manager)
    } else {
        valStateMatrix(pre_obj, omic)
    }
    model <- if (identical(ruv, "factors_in_model") && omic %in% c("lipidomics", "metabolomics")) {
        valStateMatrix(model_obj, omic, "uncorrected")
    } else {
        valStateMatrix(model_obj, omic)
    }
    out <- c(da = write(da, "da"), matrix_pre_ruv = write(pre, "matrix_pre_ruv"),
             matrix = write(model, "matrix"))
    key <- function(x) paste(x$assay, x$feature_id, sep = "\r")
    all_keys <- unique(rbind(pre[, 1:2], model[, 1:2], da[, c("feature_id", "assay")]))
    features <- data.frame(feature_id = all_keys$feature_id, assay = all_keys$assay,
                           in_matrix_pre_ruv = key(all_keys) %in% key(pre),
                           in_matrix = key(all_keys) %in% key(model),
                           tested = key(all_keys) %in% key(da[is.finite(da$p), ]),
                           stringsAsFactors = FALSE)
    out <- c(out, features = write(features, "features"))
    if (!identical(ruv, "off")) {
        out <- c(out, ruv = write(valRuvRows(model_obj, design_name), "ruv"))
    }
    if (omic %in% c("lipidomics", "metabolomics")) {
        regex <- valField(pre_obj, "internal_standard_regex")
        assays <- unique(pre$assay)
        itsd <- data.frame(assay = assays, applied = rep(itsd_applied, length(assays)),
                           standards_n = vapply(assays, function(a) {
                               ids <- pre$feature_id[pre$assay == a]
                               if (length(regex) && nzchar(regex)) sum(grepl(regex, ids)) else 0L
                           }, integer(1), USE.NAMES = FALSE),
                           stringsAsFactors = FALSE)
        out <- c(out, itsd = write(itsd, "itsd"))
    }
    if (identical(ruv, "factors_in_model")) {
        display <- if (!is.null(display_reader)) display_reader() else
            valDisplayMatrix(model_obj, omic)
        if (!is.null(display)) {
            display <- cbind(display[, 1:2], label = MATRIX_DISPLAY_LABEL,
                             display[, -(1:2), drop = FALSE])
            out <- c(out, matrix_display = write(display, "matrix_display"))
        }
    }
    out
}

#' The display matrix where the product keeps it outside the state object: DIA_LIMPA's
#' `ruv-adjusted-science` checkpoint of the `ruv_covariates` state's product (as the Summary's
#' integration input reads it, `protIntegrationInputSource`, R/utils_workflow_dia_summary.R:96-176);
#' transcriptomics' RUVs-corrected log-CPM, the `corrected` table of the `trx_ruvs` state's product
#' (R/utils_workflow_trx_declarations.R:522, read through `trxStateRuvs` and `trxSavedMatrix`).
#' Other
#' routes: `NULL` (their display is the state's own values).
#' @noRd
valProductDisplayReader <- function(manager, r) {
    ns <- asNamespace("MultiScholaR")
    scoped <- function(f) ns$withWorkflowStateScope(manager, f)
    if (identical(r$workflow, "DIA_LIMPA")) {
        return(function() {
            src <- scoped(function() ns$protIntegrationInputSource(manager))
            if (is.null(src$corrected)) return(NULL)
            valFeatureTable(src$corrected, "feature_id", "all")
        })
    }
    if (identical(r$omic, "transcriptomics")) {
        return(function() {
            m <- scoped(function() {
                ruvs <- ns$trxStateRuvs(manager, ns$getWorkflowStateReference(manager, "trx_ruvs"))
                # No display matrix where RUV was not applied (design none, or k = 0; 7.3 r1.11), by
                # the product's own rule (trxIntegrationRuvState).
                applied <- isTRUE(ruvs$outcomes$ruv_applied) && isTRUE(as.integer(ruvs$k) > 0L)
                if (!applied || is.null(ruvs$corrected)) NULL else
                    ns$trxSavedMatrix(manager, ruvs$corrected)
            })
            if (is.null(m)) return(NULL)
            cbind(data.frame(feature_id = rownames(m), assay = "all", stringsAsFactors = FALSE),
                  as.data.frame(m, stringsAsFactors = FALSE))
        })
    }
    NULL
}

#' The normalised matrix a classic proteomics route's RUV step read, from the `ruv_corrected`
#' state's own product (`normalization_science`, its `normalized` checkpoint), as the integration
#' input reads it (`protIntegrationInputRuvInput`, R/utils_workflow_dia_summary.R:180-205)
#' @noRd
valProtRuvInput <- function(manager) {
    ns <- asNamespace("MultiScholaR")
    ns$withWorkflowStateScope(manager, function() {
        lineage <- ns$protIntegrationInputLineage(manager, until = "ruv_corrected")
        record <- lineage[["ruv_corrected"]]
        if (is.null(record)) valAbort("the state object is missing", "val_spec_invalid")
        input <- ns$protIntegrationInputRuvInput(manager, record)
        tbl <- ns$readDIAExecutionComponents(manager, input$reference,
                                             fields = c("protein_quant_table",
                                                        "design_matrix"))@protein_quant_table
        valFeatureTable(tbl, record$data$slots$protein_id_column, "all")
    })
}

#' The RUV-corrected display matrix where the state holds it (lipidomics and metabolomics: the
#' state's values, `runOmicsRuvFactors()$corrected`; proteomics and transcriptomics: an `adjusted`
#' or `corrected` element of `args` where present; otherwise `NULL`)
#' @noRd
valDisplayMatrix <- function(obj, omic) {
    if (omic %in% c("lipidomics", "metabolomics")) return(valStateMatrix(obj, omic, "data"))
    args <- valField(obj, "args")
    held <- args$adjusted %||% args$corrected
    if (is.null(held)) return(NULL)
    if (is.matrix(held)) {
        held <- data.frame(feature_id = rownames(held), held, check.names = FALSE,
                           stringsAsFactors = FALSE)
    }
    valFeatureTable(held, names(held)[[1L]], "all")
}

#' The masked digest of a table (feature-level tables compared without leaving the job)
#'
#' @param table A data frame.
#' @return SHA-256 of its tab-separated text with numbers at 15 significant digits.
#' @export
valMaskedDigest <- function(table) {
    x <- table
    for (col in names(x)) if (is.numeric(x[[col]])) x[[col]] <- ifelse(is.finite(x[[col]]),
                                                                       sprintf("%.15g", x[[col]]),
                                                                       "NA")
    text <- utils::capture.output(utils::write.table(x, sep = "\t", quote = FALSE,
                                                     row.names = FALSE))
    digest::digest(paste(text, collapse = "\n"), algo = "sha256", serialize = FALSE)
}

#' Replay a cell from its analysis record (section 7.3; statuses as G5)
#'
#' `buildWorkflowAnalysisRecord`, then `replayWorkflowAnalysisRecord` on a fresh store with the
#' same inputs; the replayed DA table is compared by masked digest, then by value.
#'
#' @param manager The cell's state manager.
#' @param original_da The cell's `da.tsv` rows.
#' @param build `function(manager)` returning the record (default the package's builder).
#' @param replay `function(record)` returning the replayed `da.tsv` rows.
#' @return A list: `status` (`identical`, `differs`, `unsupported`), `class`, `max_abs_diff`,
#'   `rows_only_one_side`.
#' @export
valReplay <- function(manager, original_da, build = NULL, replay) {
    build <- build %||% function(m) asNamespace("MultiScholaR")$buildWorkflowAnalysisRecord(m)
    res <- tryCatch(list(da = replay(build(manager))), error = function(e) e)
    if (inherits(res, "error")) {
        return(list(status = "unsupported", class = class(res)[[1L]], max_abs_diff = NA_real_,
                    rows_only_one_side = NA_integer_))
    }
    da <- res$da
    if (identical(valMaskedDigest(da), valMaskedDigest(original_da))) {
        return(list(status = "identical", class = NA_character_, max_abs_diff = 0,
                    rows_only_one_side = 0L))
    }
    k <- function(x) paste(x$assay, x$contrast, x$feature_id, sep = "\r")
    common <- intersect(k(da), k(original_da))
    only <- length(setdiff(union(k(da), k(original_da)), common))
    a <- da[match(common, k(da)), , drop = FALSE]
    b <- original_da[match(common, k(original_da)), , drop = FALSE]
    num <- intersect(names(a)[vapply(a, is.numeric, logical(1))],
                     names(b)[vapply(b, is.numeric, logical(1))])
    diffs <- unlist(lapply(num, function(col) abs(a[[col]] - b[[col]])))
    list(status = "differs", class = NA_character_,
         max_abs_diff = if (length(diffs)) max(diffs, na.rm = TRUE) else NA_real_,
         rows_only_one_side = as.integer(only))
}
