# Part B (KIT-159 7.10; KIT-170 N12, N13): the integration through the exported workbook
# functions a user's script calls, with every setting recorded, declared departures as forced
# settings with their reasons, tool warnings verbatim and the maths-library pin.

#' The workbook functions of the MultiScholaR namespace (exported, `R/utils_integration_workbook.R`;
#' `integrationPinMathsLibrary` is internal, `R/utils_workflow_integration_model.R:1136`)
#'
#' @return Named list of functions.
#' @export
valIntegrationFunctions <- function() {
    ns <- asNamespace("MultiScholaR")
    names <- c("integrationSettings", "integrationFindInputs", "integrationOpenProject",
               "integrationMatchSamples", "integrationConfirmMatch", "integrationTrainModel",
               "integrationGetResults", "integrationReplayProject", "integrationPinMathsLibrary")
    stats::setNames(lapply(names, function(n) get(n, envir = ns)), names)
}

#' The settings table of `integrationSettings()` in long form
#' @noRd
valSettingsTable <- function(settings) {
    rows <- lapply(names(settings), function(table) {
        s <- settings[[table]]
        if (!is.data.frame(s) || !"setting" %in% names(s)) return(NULL)
        default <- if ("default" %in% names(s)) as.character(s$default) else
            rep(NA_character_, nrow(s))
        data.frame(table = table, setting = as.character(s$setting), default = default,
                   stringsAsFactors = FALSE)
    })
    do.call(rbind, Filter(Negate(is.null), rows))
}

#' Run the integration as a user's script does (section 7.10)
#'
#' `integrationSettings`, the maths-library pin, `integrationFindInputs`, `integrationOpenProject`,
#' `integrationMatchSamples` with the declared match table, `integrationConfirmMatch`,
#' `integrationTrainModel` with the declared departures only, `integrationGetResults`,
#' `integrationReplayProject`; `post_integration(project, results)` is then called outside test
#' mode.
#'
#' @param spec `project_dir`, `sample_match_file`, `store_root`, `forced` (setting =
#'   `list(value, reason)`), `input_matrix`.
#' @param fns The workbook functions ([valIntegrationFunctions()]).
#' @param post_integration The ticket's hook, or `NULL`.
#' @return A list: `settings` (table, setting, default, value, source, reason), `warnings`, `pin`,
#'   `results`, `replay`, `hook`.
#' @export
valRunIntegration <- function(spec, fns = valIntegrationFunctions(), post_integration = NULL) {
    warnings <- character()
    keep <- function(w) {
        warnings <<- c(warnings, conditionMessage(w))  # nolint
        invokeRestart("muffleWarning")
    }
    settings <- valSettingsTable(fns$integrationSettings())
    forced <- spec$forced %||% list()
    for (name in names(forced)) {
        if (!name %in% settings$setting) {
            valAbort("the departure {name} is not a workbook setting", "val_spec_invalid")
        }
        reason <- forced[[name]]$reason
        if (is.null(reason) || !nzchar(reason)) {
            valAbort("the departure {name} has no reason", "val_spec_invalid")
        }
    }
    hit <- match(settings$setting, names(forced))
    settings$value <- ifelse(is.na(hit), settings$default, vapply(seq_along(hit), function(i) {
        if (is.na(hit[[i]])) {
            NA_character_
        } else {
            paste(format(forced[[hit[[i]]]]$value), collapse = "|")
        }
    }, character(1)))
    settings$source <- ifelse(is.na(hit), "default", "forced")
    settings$reason <- ifelse(is.na(hit), NA_character_, vapply(seq_along(hit), function(i) {
        if (is.na(hit[[i]])) NA_character_ else forced[[hit[[i]]]]$reason
    }, character(1)))
    out <- withCallingHandlers({
        pin <- fns$integrationPinMathsLibrary()
        inputs <- fns$integrationFindInputs(spec$project_dir)
        project <- fns$integrationOpenProject(inputs, store_root = spec$store_root)
        match <- fns$integrationMatchSamples(project, spec$sample_match_file)
        confirmed <- fns$integrationConfirmMatch(project, match)
        model <- fns$integrationTrainModel(project, confirmed, input_matrix = spec$input_matrix,
                                           overrides = lapply(forced, function(f) f$value))
        results <- fns$integrationGetResults(project, model)
        replay <- fns$integrationReplayProject(spec$store_root %||% spec$project_dir)
        list(pin = pin, project = project, results = results, replay = replay)
    }, warning = keep)
    hook <- if (is.function(post_integration)) {
        withr::with_options(list(multischolar.test_mode = FALSE), {
            withr::with_envvar(c(MULTISCHOLAR_TEST_MODE = "false"),
                               post_integration(out$project, out$results))
        })
    }
    list(settings = settings, warnings = warnings, pin = out$pin, results = out$results,
         replay = out$replay, hook = hook)
}
