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

# KIT-102 W4: the shared app shell of proteomics, lipidomics, metabolomics and transcriptomics, configured by
# `omicGuiConfig(<omic>)$shell` (utils_omic_gui_config.R). The public `mod_<omic>_ui/server` functions are thin wrappers
# over it. Every DOM id, test id, tabset id, tab value, module id and status key comes from the configuration and is
# today's; the markup is byte-identical to the four shells it replaces. Recorded divergence from conventions/R-shiny.md:
# the fluidRow and pill-tabset layout and the per-omic tab modules are kept (the byte-identical contract), under the
# existing-code rule of conventions/R.md.

# ------------------------------------------------------------------------------------------------
# The shell.
# ------------------------------------------------------------------------------------------------

#' The shared workflow shell's UI: the header (none for proteomics), the stepper, the extras' panels and the tabset.
#' @param id Module id (`<omic>_workflow` in the app).
#' @param config An omic GUI configuration (`omicGuiConfig()`).
#' @noRd
mod_omic_workflow_ui <- function(id, config) {
    validateOmicGuiConfig(config)
    shell <- config$shell
    ns <- shiny::NS(id)
    panels <- lapply(shell$tabs, omicWorkflowTabPanel, style = shell$tab_style, ns = ns)
    extras <- lapply(shell$extras, function(extra) omicGuiFunction(omicGuiComponent(extra)$ui)(ns))
    tabset <- testid(do.call(shiny::tabsetPanel, c(list(id = ns(shell$tabset_id), type = "pills"), panels)),
        paste0(shell$testid, "-workflow-tabs"))
    progress <- shiny::uiOutput(ns("workflow_progress"))
    if (is.null(shell$header_icon)) {
        # No header (proteomics): the stepper in its own row, the extras, then the tabs.
        return(shiny::tagList(shiny::fluidRow(shiny::column(12, progress)), extras, tabset))
    }
    header <- shiny::h2(shiny::icon(shell$header_icon), paste0(" ", shell$label, " Workflow"),
        style = "margin-bottom: 20px;")
    shiny::tagList(shiny::fluidRow(shiny::column(12, header, progress, extras, tabset)))
}

#' One workflow tab: its component's UI under the tab's test id. `icon_argument` (proteomics) gives the icon as
#' `tabPanel`'s argument and a line break; `icon_in_title` puts it in the title.
#' @noRd
omicWorkflowTabPanel <- function(tab, style, ns) {
    body <- shiny::div(`data-testid` = tab$testid, omicGuiComponentUi(tab$component, ns(tab$module_id)))
    if (identical(style, "icon_argument")) {
        return(shiny::tabPanel(tab$title, value = tab$value, icon = shiny::icon(tab$icon), shiny::br(), body))
    }
    shiny::tabPanel(title = shiny::tagList(shiny::icon(tab$icon), tab$title), value = tab$value, body)
}

#' The shared workflow shell's server: the bus and its registration, the session's state owner, the omics runner,
#' the extras, the tab servers in tab order, `active_tab`, the stepper and the status chain.
#' @param id Module id.
#' @param project_dirs Project directories by omic (`setupDirectories()`).
#' @param omic_type The omic, the key of `project_dirs`.
#' @param experiment_label The experiment label.
#' @param volumes Kept for the wrappers' signature; the tabs get `volumes = NULL`, as today.
#' @param config An omic GUI configuration (`omicGuiConfig()`).
#' @return The bus, or NULL when `project_dirs` has no entry for `omic_type`.
#' @noRd
mod_omic_workflow_server <- function(id, project_dirs, omic_type, experiment_label, volumes = NULL, config) {
    validateOmicGuiConfig(config)
    # The module function's frame holds only the names the pre-KIT-102 shells held (workflow_data, experiment_paths,
    # qc_trigger, selected_tab): `testServer` evaluates a caller's expression with this frame and its parent as a data
    # mask, so any other local would shadow the caller's own variables (the entry-point suite's `steps`, for one).
    shiny::moduleServer(id, function(input, output, session) {
        workflow_data <- do.call(shiny::reactiveValues, omicWorkflowInitialBus(config$shell))
        # The bus owns the session's progress trackers and payloads (KIT-013, KIT-030); helpers find it here.
        if (identical(config$shell$bus_registration, "userData")) {
            session$userData$workflow_data <- workflow_data
        } else {
            registerWorkflowOmicsBus(session, omic_type, workflow_data)
        }
        if (!omic_type %in% names(project_dirs)) {
            logger::log_error(paste("No directory information found for", omic_type))
            return()
        }
        experiment_paths <- project_dirs[[omic_type]]
        qc_trigger <- shiny::reactiveVal(NULL)
        selected_tab <- shiny::reactive(input[[config$shell$tabset_id]])
        startOmicWorkflow(input, output, session, config, list(workflow_data = workflow_data,
            project_dirs = project_dirs, experiment_paths = experiment_paths, omic_type = omic_type,
            experiment_label = experiment_label, qc_trigger = qc_trigger, selected_tab = selected_tab))
        workflow_data
    })
}

#' The shell's start, after the bus: the session's state owner (KIT-048, KIT-070), the omics runner under a managed
#' application (KIT-078; proteomics' owner is the DIA alias with its default omic and runs its managed steps per tab),
#' the extras, `active_tab` (for collect_state_digest), the tab servers in tab order, the stepper and the chain.
#' @param context The shell's values the components take: `workflow_data`, `project_dirs`, `experiment_paths`,
#'   `omic_type`, `experiment_label`, `qc_trigger` and `selected_tab`.
#' @noRd
startOmicWorkflow <- function(input, output, session, config, context) {
    shell <- config$shell
    workflow_data <- context$workflow_data
    if (identical(config$binding$runner, "omics_actions")) {
        workflow_data$state_manager <- config$binding$owner(workflow_data, context$experiment_paths, session,
            omic_type = context$omic_type)
        installWorkflowOmicsActions(session, workflow_data, context$omic_type, context$experiment_label,
            context$experiment_paths)
    } else {
        workflow_data$state_manager <- config$binding$owner(workflow_data, context$experiment_paths, session)
    }
    for (extra in shell$extras) omicGuiFunction(omicGuiComponent(extra)$server)(input, output, session, context)
    workflow_data$active_tab <- context$selected_tab
    for (tab in shell$tabs) omicGuiComponentServer(tab$component, tab$module_id, context)
    steps <- omicWorkflowSteps(shell$tabs)
    output$workflow_progress <- shiny::renderUI(render_workflow_stepper(steps, workflow_data$tab_status))
    installOmicWorkflowChain(workflow_data, shell)
    invisible(NULL)
}

#' The initial bus: the configured fields in their order, `tab_status` in its place with the first tab `pending` and
#' the others `disabled`.
#' @noRd
omicWorkflowInitialBus <- function(shell) {
    keys <- vapply(shell$tabs, `[[`, character(1), "status_key")
    status <- stats::setNames(as.list(rep("disabled", length(keys))), keys)
    status[[1L]] <- "pending"
    bus <- shell$bus
    bus["tab_status"] <- list(status)
    bus
}

#' The stepper's steps (`render_workflow_stepper`), one per tab.
#' @noRd
omicWorkflowSteps <- function(tabs) {
    lapply(tabs, function(tab) list(name = tab$step_name, key = tab$status_key, icon = tab$step_icon))
}

#' The status chain, one observer per link in tab order: when a tab completes, the next opens. The import opens the
#' step after it itself (`completeOmicImport`), so the chain starts at the second tab. With `chain_rules$qc_bypass`
#' (proteomics), a route in it marks the QC tab complete and opens the tab after it instead.
#' @noRd
installOmicWorkflowChain <- function(workflow_data, shell) {
    keys <- vapply(shell$tabs, `[[`, character(1), "status_key")
    bypass_routes <- shell$chain_rules$qc_bypass
    for (i in seq_along(keys)[-c(1L, length(keys))]) local({
        done <- keys[[i]]
        next_key <- keys[[i + 1L]]
        after <- if (i + 2L <= length(keys)) keys[[i + 2L]] else NULL
        bypassable <- !is.null(bypass_routes) && identical(next_key, "quality_control")
        shiny::observeEvent(workflow_data$tab_status[[done]], {
            status <- workflow_data$tab_status
            bypass <- FALSE
            if (bypassable && identical(status[[done]], "complete")) {
                workflow_type <- shiny::isolate(workflow_data$state_manager$workflow_type)
                bypass <- workflow_type %in% bypass_routes
                if (bypass) logger::log_info(sprintf("%s workflow detected, bypassing QC tab.", workflow_type))
            }
            updated <- advanceOmicWorkflowStatus(status, done, next_key, after, bypass)
            if (!identical(updated, status)) workflow_data$tab_status <- updated
        }, ignoreNULL = TRUE)
    })
    invisible(NULL)
}

#' One link of the chain: when `done` is complete, a disabled `next_key` becomes pending; with `bypass`, an incomplete
#' `next_key` becomes complete and a disabled `after` pending (mod_proteomics.R's TMT/LFQ rule before KIT-102).
#' @return The statuses, changed or not.
#' @noRd
advanceOmicWorkflowStatus <- function(status, done, next_key, after = NULL, bypass = FALSE) {
    if (!identical(status[[done]], "complete")) return(status)
    if (!bypass) {
        if (identical(status[[next_key]], "disabled")) status[[next_key]] <- "pending"
        return(status)
    }
    if (!identical(status[[next_key]], "complete")) {
        status[[next_key]] <- "complete"
        if (!is.null(after) && identical(status[[after]], "disabled")) status[[after]] <- "pending"
    }
    status
}

# ------------------------------------------------------------------------------------------------
# The component registry.
# ------------------------------------------------------------------------------------------------

#' A GUI component by id: today's public tab module (`ui`, `server`: function names, looked up when the UI is built or
#' the server starts, so a test's namespace stub replaces them), how the shell calls its server (`call`) and whether it
#' takes the transcriptomics profile. An extra (`peptide_checkpoint`) has `ui(ns)` and `server(input, output, session,
#' context)` only.
#' @param id A component id, as `shell$tabs[[i]]$component` or `shell$extras` name it.
#' @noRd
omicGuiComponent <- function(id) {
    components <- omicGuiComponents()
    if (!is.character(id) || length(id) != 1L || !id %in% names(components)) {
        cli::cli_abort("{.val {id}} is not a GUI component.", class = "omic_gui_config_invalid")
    }
    components[[id]]
}

#' The registry. Lipid and metab call their own tab modules (the pre-KIT-102 mod_lipidomics.R:207-263, the same lines
#' in mod_metabolomics.R); trx calls the lipid modules under its profile, its normalisation without `selected_tab`
#' (mod_transcriptomics.R:258-264); proteomics its own (mod_proteomics.R:328-374).
#' @noRd
omicGuiComponents <- function() {
    component <- function(module, call, profile = FALSE) {
        list(ui = paste0("mod_", module, "_ui"), server = paste0("mod_", module, "_server"), call = call,
            profile = profile)
    }
    feature <- function(prefix, module = prefix, profile = FALSE) {
        calls <- c(import = "import", design = "design", qc = "qc", norm = if (profile) "tab" else "tab_selected",
            da = "tab", summary = "summary")
        stats::setNames(lapply(names(calls), function(tab) component(paste0(module, "_", tab), calls[[tab]], profile)),
            paste0(prefix, "_", names(calls)))
    }
    prot <- list(prot_import = component("prot_import", "import"), prot_design = component("prot_design", "design"),
        prot_qc = component("prot_qc", "qc_named"), prot_norm = component("prot_norm", "tab_selected_positional"),
        prot_da = component("prot_da", "tab_selected_positional"),
        prot_enrich = component("prot_enrich", "tab_selected_positional"),
        prot_summary = component("prot_summary", "summary"))
    extras <- list(peptide_checkpoint = list(ui = "protPeptideCheckpointUi", server = "protPeptideCheckpointServer"))
    c(feature("lipid"), feature("metab"), feature("trx", "lipid", profile = TRUE), prot, extras)
}

#' A component's UI under its module id.
#' @noRd
omicGuiComponentUi <- function(id, module_id) {
    component <- omicGuiComponent(id)
    ui <- omicGuiFunction(component$ui)
    if (component$profile) ui(module_id, profile = trxOmicProfile()) else ui(module_id)
}

#' A component's server under its module id, with the arguments its call names.
#' @noRd
omicGuiComponentServer <- function(id, module_id, context) {
    component <- omicGuiComponent(id)
    server <- omicGuiFunction(component$server)
    if (component$profile) {
        return(callOmicGuiComponent(component$call, server, module_id, context, profile = trxOmicProfile()))
    }
    callOmicGuiComponent(component$call, server, module_id, context)
}

#' The tab servers' argument shapes, as the pre-KIT-102 shells passed them; `...` is trx's `profile`.
#' @noRd
callOmicGuiComponent <- function(call, server, module_id, context, ...) {
    data <- context$workflow_data
    paths <- context$experiment_paths
    omic <- context$omic_type
    label <- context$experiment_label
    switch(call,
        import = server(module_id, data, paths, volumes = NULL, ...),
        design = server(module_id, data, paths, volumes = NULL, qc_trigger = context$qc_trigger, ...),
        qc = server(module_id, data, paths, omic, label, context$qc_trigger, ...),
        qc_named = server(module_id, data, paths, omic, label, qc_trigger = context$qc_trigger, ...),
        tab = server(module_id, data, paths, omic, label, ...),
        tab_selected = server(module_id, data, paths, omic, label, selected_tab = context$selected_tab, ...),
        tab_selected_positional = server(module_id, data, paths, omic, label, context$selected_tab, ...),
        summary = server(module_id, context$project_dirs, omic, label, data, ...),
        cli::cli_abort("{.val {call}} is not a component call.", class = "omic_gui_config_invalid"))
}

#' A GUI function by name, looked up in the package namespace when it is called; a function is returned as it is.
#' @noRd
omicGuiFunction <- function(name) {
    if (is.function(name)) return(name)
    get(name, mode = "function")
}
