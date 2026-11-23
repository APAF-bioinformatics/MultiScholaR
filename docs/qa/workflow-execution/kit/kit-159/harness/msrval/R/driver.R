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
    put <- function(id, value) out[id] <<- list(value)
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
    ids <- names(st$values)
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
