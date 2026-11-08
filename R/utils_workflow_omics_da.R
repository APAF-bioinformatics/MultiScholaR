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

# KIT-106 D3, D4: the shared DA tab's core for the omics that hold rows (lipidomics, metabolomics, transcriptomics). The
# interactive volcano in process is the binding's `<prefix>_volcano_page`, as the worker makes it; the session load is the
# file, the read, the assay check, the adoption and the restore. Neither shows a notification or updates an input: the
# tab does that with what they return (KIT-105 D9).

#' A volcano widget generator for an in-process session, the sibling of `workflowOmicsVolcanoPageFactory`: the binding's
#' `<prefix>_volcano_page` publication run here (`runWorkflowOmicsPublicationInProcess`), which reads the state in the
#' core and keeps the page as a store asset under its identity, shown in the managed factory's iframe. The tab holds no
#' state, widget or page bytes.
#' @param binding The omic's declaration binding.
#' @param manager The session's state owner.
#' @param experiment_paths The experiment's directories, by role.
#' @return A function with the generator's arguments (`da_results_list` is not read: the page reads the head's manifest)
#'   giving the iframe, or NULL when the widget reported nothing. A failure to draw is raised with the action's message.
#' @noRd
workflowOmicsInProcessVolcanoPageFactory <- function(binding, manager, experiment_paths) {
    function(da_results_list, selected_contrast, selected_assay, da_q_val_thresh) {
        display <- runWorkflowOmicsPublicationInProcess(binding, manager, paste0(binding$prefix, "_volcano_page"),
            list(contrast = selected_contrast, assay = selected_assay, qvalue = as.numeric(da_q_val_thresh)), experiment_paths)
        if (!is.null(display$message)) stop(display$message, call. = FALSE)
        if (is.null(display$page)) return(NULL)
        htmltools::tags$iframe(srcdoc = rawToChar(hydrateSavedWorkflowEntity(manager$.context$store, display$page)$bytes),
            style = "width:100%;height:820px;border:0;", title = paste("Interactive volcano", selected_contrast))
    }
}

#' The DA tab's session load (lipid's mod_lipid_da_session_helpers.R:1-690 and metab's mod_metab_da_observer_helpers.R
#' :446-639, :760-930 at 2269a45f): the session file in the experiment's source directory (else its export directory),
#' the owner's own export reopened or the file read, the assay check (the configuration's rule and text), a legacy session
#' adopted into a fresh owner on the bus (`ingestWorkflowOmicsLegacySession`), then the restore: the contrasts on the bus
#' and the tab, the assay names, the state name and the formula on the tab.
#' @param binding The omic's binding (`omic_type`, `class`, `data_slot`).
#' @param workflow_data The omic's bus (`state_manager`, `contrasts_tbl`).
#' @param tab The tab's values (`current_state_name`, `contrasts_available`, `assays_available`, `formula_from_s4`).
#' @param experiment_paths The experiment's directories (`source_dir`, `export_dir`).
#' @param session_load The configuration's `da$session_load`: `file`; `assay_rule`, `setequal` or `identical`;
#'   `assay_message`, a format of the session's and the object's assay names; `assay_omic_type`, when the check applies
#'   only to a session of that omic; `assay_check_first`, TRUE to check before the adoption (metab's order), else after
#'   the contrasts (lipid's).
#' @param reading `function()` called once the file is found, before it is read (the tab's loading notification).
#' @return `status`: "missing" (no file or directory; `message` is the tab's text), "error" (the read, the adoption or the
#'   check failed: `message`; what was restored before it is set and its choices returned) or "success". With the
#'   restore: `contrast_choices`, `assay_choices` and `table_assay_choices` (the selectors' choices) and `formula`.
#' @noRd
loadWorkflowOmicsDaSession <- function(binding, workflow_data, tab, experiment_paths, session_load, reading = function() NULL) {
    source_dir <- experiment_paths$source_dir
    if (is.null(source_dir) || !dir.exists(source_dir)) source_dir <- experiment_paths$export_dir
    if (is.null(source_dir) || !dir.exists(source_dir)) {
        return(list(status = "missing", message = "Could not find source directory for session data."))
    }
    path <- file.path(source_dir, session_load$file)
    if (!file.exists(path)) return(list(status = "missing", message = sprintf("Session file not found: %s", path)))
    reading()
    out <- list(status = "success")
    tryCatch({
        # KIT-071: the owner's own export reopens the owner (the head's state, no file read); any other file is read.
        session_data <- resolveWorkflowOmicsOwnedSession(workflow_data$state_manager, path)
        if (is.null(session_data)) session_data <- readRDS(path)
        logger::log_info(sprintf("   Loaded session from: %s", path))
        object <- session_data$current_s4_object
        assay_names <- session_data$assay_names
        check_assays <- function() {
            if (is.null(assay_names) || !inherits(object, binding$class)) return(invisible(NULL))
            if (!is.null(session_load$assay_omic_type) && !identical(session_data$omic_type, session_load$assay_omic_type)) {
                return(invisible(NULL))
            }
            actual <- names(methods::slot(object, binding$data_slot))
            differs <- if (identical(session_load$assay_rule, "setequal")) {
                length(assay_names) > 0L && !setequal(as.character(assay_names), as.character(actual))
            } else {
                length(actual) > 0L && !identical(as.character(assay_names), actual)
            }
            if (differs) {
                stop(sprintf(session_load$assay_message, paste(assay_names, collapse = ", "), paste(actual, collapse = ", ")),
                    call. = FALSE)
            }
            invisible(NULL)
        }
        if (isTRUE(session_load$assay_check_first)) check_assays()
        if (!is.null(object)) {
            state_name <- session_data$r6_current_state_name %||% "loaded_for_de"
            if (isTRUE(session_data$owned_generation)) {
                # KIT-071: the owner's own export reopens the owner; the head's state is the tab's and nothing is recommitted.
                logger::log_info("   Reopened the session owner at its head")
            } else {
                # KIT-071: a session from an older build adopts into a fresh owner on the bus (provenance `adopted`;
                # KIT-072: a manager that is not an owner is refused).
                ingestWorkflowOmicsLegacySession(session_data, workflow_data,
                    workflow_data$state_manager$.context$experiment_paths$source_dir %||% dirname(path), binding$omic_type,
                    state_name = state_name, path = path)
                logger::log_info("   Adopted the legacy session into a fresh state owner")
            }
            tab$current_state_name <- state_name
        }
        contrasts <- session_data$contrasts_tbl
        if (!is.null(contrasts) && nrow(contrasts) > 0) {
            workflow_data$contrasts_tbl <- contrasts
            tab$contrasts_available <- contrasts
            out$contrast_choices <- if ("friendly_names" %in% colnames(contrasts)) contrasts$friendly_names else contrasts$contrasts
            logger::log_info(sprintf("   Restored %d contrasts", nrow(contrasts)))
        }
        if (!is.null(assay_names)) {
            if (!isTRUE(session_load$assay_check_first)) check_assays()
            tab$assays_available <- assay_names
            out$assay_choices <- c("Combined", assay_names)
            out$table_assay_choices <- c("All", assay_names)
        }
        tryCatch({
            # KIT-075: the owner's own export carries the state's reference and its slots (`state_metadata`); a session
            # file from an older build carries the object.
            args <- if (!isS4(object) && !is.null(session_data$state_metadata)) session_data$state_metadata$args else object@args
            formula <- args$daAnalysisParameters$formula_string
            if (!is.null(formula) && nzchar(formula)) {
                tab$formula_from_s4 <- formula
                out$formula <- formula
            }
        }, error = function(error) logger::log_warn(sprintf("Could not extract formula from S4: %s", conditionMessage(error))))
        out
    }, error = function(error) {
        logger::log_error(sprintf("   Error loading session: %s", conditionMessage(error)))
        out$status <- "error"
        out$message <- conditionMessage(error)
        out
    })
}
