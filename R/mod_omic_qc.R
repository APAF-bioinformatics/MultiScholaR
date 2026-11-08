# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# KIT-104 W5: the shared QC shell (kit-101 design 4.6), the frame of every omic's QC tab. It renders
# the well panel, the title and the tabset of `config$qc` (the lipid orchestrator's markup,
# mod_lipid_qc.R at d8f6dde4; proteomics' named groups, mod_prot_qc.R), starts the sub-tabs once
# (on the design's trigger, or when the state's class is the binding's with `auto_init`) and
# completes the tab on `complete_state`. The sub-tabs and group plugins are looked up by name in
# the plugin registry (`omicQcComponents()`) when they are called, so a test's namespace mock
# replaces them. The four shared sub-tabs are in mod_omic_qc_tabs.R; their computation is the
# core's (func_omics_qc_shared.R, `runWorkflowOmicsQcAction`). Recorded divergences from
# conventions/R-shiny.md: the wellPanel layout and the in-process or managed actions in place of
# ExtendedTask are kept, under the existing-code scope rule of conventions/R.md.

# ------------------------------------------------------------------------------------------------
# The click guard and the shared revert.
# ------------------------------------------------------------------------------------------------

#' Whether a QC sub-tab's button holds a click made after its server started (KIT-101 F-17, F-24).
#'
#' The QC sub-module servers start at the design; a click made before then is already in the input when the server
#' starts, and its observer's first run would run it. Called when the server starts, this records the click value then
#' present; the returned function is FALSE for that value only, so a later click runs as before (`ignoreInit` would also
#' drop a first click made in the same flush as the start). Only a live input can hold an earlier click.
#' @param input The sub-module's input.
#' @param id The button's input id.
#' @return A function of no arguments: TRUE when the button's value is not the one present at the start.
#' @noRd
qcClickAfterStart <- function(input, id) {
    at_start <- if (shiny::is.reactivevalues(input)) shiny::isolate(input[[id]])
    function() !identical(input[[id]], at_start)
}

#' A QC sub-tab's Revert button: the state before the head selected without materialising it (KIT-075), the sub-tab's
#' results text saying so and its display cleared. A click made before the server started is not run (F-24). The
#' lipid intensity and duplicate reverts (mod_lipid_qc_intensity_server_helpers.R:85-124,
#' mod_lipid_qc_duplicates_server_helpers.R:488-550 at d8f6dde4) and the one trx borrows were this.
#' @param input,output The sub-tab's input and output.
#' @param button The Revert button's input id.
#' @param workflow_data The workflow bus.
#' @param results_output The output id of the sub-tab's results text.
#' @param clear Function of no arguments clearing the sub-tab's display.
#' @param log_label The log line's start ("Reverted lipid intensity filter to").
#' @return `input`, invisibly.
#' @noRd
registerOmicQcRevertObserver <- function(input, output, button, workflow_data, results_output, clear, log_label) {
    click_after_start <- qcClickAfterStart(input, button)
    shiny::observeEvent(input[[button]], {
        if (!click_after_start()) return(invisible(NULL))
        tryCatch({
            history <- workflow_data$state_manager$getHistory()
            if (length(history) <= 1) stop("No previous state to revert to.")
            previous <- history[[length(history) - 1]]
            selectWorkflowStateReference(workflow_data$state_manager, previous)
            output[[results_output]] <- shiny::renderText(paste("Reverted to previous state:", previous))
            clear()
            logger::log_info(paste(log_label, previous))
            shiny::showNotification("Reverted successfully", type = "message")
        }, error = function(e) {
            message <- paste("Error reverting:", e$message)
            logger::log_error(message)
            shiny::showNotification(message, type = "error")
        })
    })
    invisible(input)
}

# ------------------------------------------------------------------------------------------------
# The plugin registry.
# ------------------------------------------------------------------------------------------------

#' The QC plugin registry: each component a qc group or tab names (`config$qc$groups`), with its module pair (`ui`,
#' `server`: function names, looked up when called) and how the shell calls it (`call`):
#' - `shared`: the four shared sub-tabs, `ui(id, config)` and `server(id, workflow_data, omic_type, experiment_label,
#'   config)`;
#' - `plugin`: trx's expression filter, `ui(id)` and `server(id, workflow_data, omic_type, experiment_label)`;
#' - `group`: proteomics' group plugins, `ui(id, <route_arg> = route)` and `server(id, workflow_data,
#'   experiment_paths, omic_type, experiment_label)` (mod_prot_qc.R:140-141, :176-182);
#' - `subtab`: proteomics' 17 sub-tabs, which their group plugin starts.
#' @return A named list by component id.
#' @noRd
omicQcComponents <- function() {
    pair <- function(stem, call, route_arg = NULL) {
        list(ui = paste0("mod_", stem, "_ui"), server = paste0("mod_", stem, "_server"), call = call,
            route_arg = route_arg)
    }
    shared <- c("omic_qc_intensity", "omic_qc_duplicates", "omic_qc_itsd", "omic_qc_finalize")
    components <- c(stats::setNames(lapply(shared, pair, "shared"), shared),
        list(trx_qc_filter = pair("trx_qc_filter", "plugin"),
            prot_qc_peptide = pair("prot_qc_peptide", "group", "workflowType"),
            prot_qc_protein = pair("prot_qc_protein", "group", "workflow_type")))
    stems <- unique(unlist(lapply(omicGuiProteomicsQc()$groups, function(group) {
        vapply(group$tabs, `[[`, "", "component")
    })))
    c(components, stats::setNames(lapply(stems, pair, "subtab"), stems))
}

#' One registry entry by component id.
#' @noRd
omicQcComponent <- function(id) {
    components <- omicQcComponents()
    if (!is.character(id) || length(id) != 1L || !id %in% names(components)) {
        cli::cli_abort("{.val {id}} is not a QC component.", class = "omic_gui_config_invalid")
    }
    components[[id]]
}

# ------------------------------------------------------------------------------------------------
# The tabset and the start.
# ------------------------------------------------------------------------------------------------

#' The groups of `qc` shown on `route`, each with its tabs on that route; a group with none is not shown (proteomics'
#' peptide group on TMT and LFQ, mod_prot_qc.R:98-133).
#' @noRd
omicQcRouteGroups <- function(qc, route) {
    groups <- lapply(qc$groups, function(group) {
        group$tabs <- Filter(function(tab) route %in% tab$routes, group$tabs)
        group
    })
    Filter(function(group) length(group$tabs) > 0L, groups)
}

#' The QC tabset on `route`: the flat tabs of an unnamed group (lipid, metab, trx; mod_lipid_qc.R:119-125), or one
#' `tabPanel(<group title>, <group plugin UI>)` per named group (proteomics, mod_prot_qc.R:121-142).
#' @param ns The shell's namespace function.
#' @param config The omic's GUI configuration.
#' @param route The route shown.
#' @noRd
omicQcTabsetUi <- function(ns, config, route) {
    qc <- config$qc
    tabset_id <- if (length(qc$tabset_id) > 1L) qc$tabset_id[[route]] else qc$tabset_id
    panels <- lapply(omicQcRouteGroups(qc, route), function(group) {
        if (is.null(group$id)) {
            return(lapply(group$tabs, function(tab) {
                component <- omicQcComponent(tab$component)
                ui <- omicGuiFunction(component$ui)
                if (identical(component$call, "shared")) ui(ns(tab$id), config) else ui(ns(tab$id))
            }))
        }
        component <- omicQcComponent(group$component)
        ui <- omicGuiFunction(component$ui)
        list(shiny::tabPanel(group$title, do.call(ui, c(list(ns(group$id)),
            stats::setNames(list(route), component$route_arg)))))
    })
    do.call(shiny::tabsetPanel, c(list(id = ns(tabset_id)), unlist(panels, recursive = FALSE)))
}

#' Start the QC tab's sub-tab servers (or group plugins) shown on `route`, in the configuration's order: the shared
#' sub-tabs and trx's filter under their tab ids (mod_lipid_qc.R:57-85), proteomics' groups under their group ids
#' (mod_prot_qc.R:163-201). The shell calls this once.
#' @param workflow_data,experiment_paths,omic_type,experiment_label The QC tab's arguments.
#' @param config The omic's GUI configuration.
#' @param route The route the tab runs on.
#' @return NULL, invisibly.
#' @noRd
startOmicQcSubtabs <- function(workflow_data, experiment_paths, omic_type, experiment_label, config, route) {
    for (group in omicQcRouteGroups(config$qc, route)) {
        if (!is.null(group$id)) {
            server <- omicGuiFunction(omicQcComponent(group$component)$server)
            server(group$id, workflow_data, experiment_paths, omic_type, experiment_label)
            next
        }
        for (tab in group$tabs) {
            component <- omicQcComponent(tab$component)
            server <- omicGuiFunction(component$server)
            if (identical(component$call, "shared")) {
                server(tab$id, workflow_data, omic_type, experiment_label, config)
            } else {
                server(tab$id, workflow_data, omic_type, experiment_label)
            }
        }
    }
    invisible(NULL)
}

# ------------------------------------------------------------------------------------------------
# The shell.
# ------------------------------------------------------------------------------------------------

#' The QC tab's UI: the well panel with the title and the slot the server renders the tabset into (mod_lipid_qc.R
#' :38-51, mod_prot_qc.R:37-50).
#' @param id Module id (the QC tab's).
#' @param config An omic GUI configuration (`omicGuiConfig()`).
#' @noRd
mod_omic_qc_ui <- function(id, config) {
    validateOmicGuiConfig(config)
    ns <- shiny::NS(id)
    shiny::tagList(shiny::fluidRow(shiny::column(12, shiny::wellPanel(
        shiny::h3(config$qc$title),
        shiny::uiOutput(ns("dynamic_qc_tabs"))))))
}

#' The QC tab's server: the tabset once the tab is ready, the sub-tabs started once, and the completion.
#' @param id Module id.
#' @param workflow_data The workflow bus.
#' @param experiment_paths The experiment's project paths (proteomics' groups take them).
#' @param omic_type,experiment_label The tab's omic and label, passed to the sub-tabs.
#' @param qc_trigger The design's QC trigger (a reactive), TRUE once the design is saved.
#' @param config An omic GUI configuration (`omicGuiConfig()`).
#' @noRd
mod_omic_qc_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, qc_trigger = NULL,
    config) {
    validateOmicGuiConfig(config)
    shiny::moduleServer(id, function(input, output, session) {
        # The route the tab shows; proteomics' tests read it by this name (test-perf107-audit.R).
        qc_workflow_type <- shiny::reactiveVal(NULL)
        registerOmicQcShell(output, session, workflow_data, experiment_paths, omic_type, experiment_label, qc_trigger,
            config, qc_workflow_type, start = function(route) {
                startOmicQcSubtabs(workflow_data, experiment_paths, omic_type, experiment_label, config, route)
            })
    })
}

#' The shell's server inside its module session, in the orchestrators' order (mod_lipid_qc.R:87-161,
#' mod_prot_qc.R:69-218).
#' - Ready: with one route, when the binding's `read_metadata` answers (else the design alert); with several
#'   (proteomics), when the route is known (`route`, set at priority 100 from the owner on the trigger and the state
#'   update, mod_prot_qc.R:81-88).
#' - Start once: on `qc_trigger() == TRUE`, or with `auto_init` when the metadata's class is the binding's.
#' - Completion: with a `complete_state` on the route, when that state exists (mod_prot_qc.R:204-216).
#' @param route The shell's `qc_workflow_type` reactiveVal.
#' @param start `function(route)`, the sub-tabs' start.
#' @noRd
registerOmicQcShell <- function(output, session, workflow_data, experiment_paths, omic_type, experiment_label,
    qc_trigger, config, route, start) {
    ns <- session$ns
    qc <- config$qc
    read_metadata <- config$binding$read_metadata
    class <- config$binding$binding$class
    routes <- config$shell$routes
    by_route <- length(routes) > 1L
    started <- shiny::reactiveVal(FALSE)

    if (by_route) {
        # State-owner adoption does not change the scientific workflow type: the tabset's DOM is kept when only the
        # owner changes. The native WorkflowState changes its field in place during the design; the trigger notifies.
        shiny::observe({
            if (is.function(qc_trigger)) qc_trigger()
            workflow_data$state_update_trigger
            route(workflow_data$state_manager$workflow_type)
        }, priority = 100)
    } else {
        route(routes[[1L]])
    }

    output$dynamic_qc_tabs <- shiny::renderUI({
        if (by_route) {
            shiny::req(route())
        } else {
            # KIT-075: the state's metadata answers whether there is data; no object is read.
            has_data <- tryCatch({
                !is.null(workflow_data$state_manager) && !is.null(read_metadata(workflow_data$state_manager))
            }, error = function(e) FALSE)
            if (!has_data) {
                return(shiny::div(class = "alert alert-info", shiny::icon("info-circle"),
                    " Please complete the 'Design Matrix' step first. QC tabs will appear here once data is available."))
            }
        }
        omicQcTabsetUi(ns, config, route())
    })

    start_once <- function() {
        if (shiny::isolate(started())) return(invisible(NULL))
        on_route <- if (by_route) shiny::isolate(workflow_data$state_manager$workflow_type) else routes[[1L]]
        logger::log_info(paste(config$shell$label, "QC: starting the sub-modules on", on_route %||% "no route"))
        start(on_route)
        started(TRUE)
    }

    # The design's trigger starts the sub-modules once.
    shiny::observeEvent(qc_trigger(), {
        shiny::req(qc_trigger() == TRUE)
        start_once()
    }, ignoreNULL = TRUE, once = TRUE)

    # Or the state: once its class is the binding's (BUG-10: the flag is read in isolation, so this runs once per
    # state-manager appearance).
    if (isTRUE(qc$auto_init)) {
        shiny::observe({
            shiny::req(workflow_data$state_manager)
            current <- tryCatch(read_metadata(workflow_data$state_manager), error = function(e) NULL)
            if (!is.null(current) && identical(current$class, class) && !shiny::isolate(started())) start_once()
        })
    }

    # The completion: the route's final state reached (proteomics, `completeProtQcWorkflowStatus`).
    if (length(qc$complete_state) > 0L) {
        shiny::observeEvent({
            state <- unname(qc$complete_state[route() %||% ""])
            if (length(state) != 1L || is.na(state)) {
                NULL
            } else if (inherits(workflow_data$state_manager, "DIAStateAdapter")) {
                workflow_data$state_update_trigger
                if (hasDIAState(workflow_data$state_manager, state)) TRUE else NULL
            } else {
                workflow_data$state_manager$states[[state]]
            }
        }, {
            logger::log_info("Final QC state reached. Quality Control module is ready.")
            completeProtQcWorkflowStatus(workflow_data)
        })
    }

    invisible(NULL)
}
