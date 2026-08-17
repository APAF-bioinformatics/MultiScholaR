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

#' @title Check if test mode is active
#' @description Returns TRUE when the app is running in test mode. Checks the
#'   R option `multischolar.test_mode` first, then falls back to the environment
#'   variable `MULTISCHOLAR_TEST_MODE`. Result is never cached so that
#'   \code{withr::with_options()} works correctly in tests.
#' @return Logical scalar.
#' @export
is_test_mode <- function() {
  opt <- getOption("multischolar.test_mode", FALSE)
  if (isTRUE(opt)) return(TRUE)
  env <- Sys.getenv("MULTISCHOLAR_TEST_MODE", "false")
  tolower(env) %in% c("true", "1")
}

#' @title Append data-testid attribute to a Shiny tag
#' @description Adds a \code{data-testid} HTML attribute to a Shiny tag object,
#'   enabling stable selectors in end-to-end tests.
#' @param tag A Shiny tag object.
#' @param id Character scalar — the test identifier value.
#' @return The modified tag with the \code{data-testid} attribute.
#' @export
testid <- function(tag, id) {
  shiny::tagAppendAttributes(tag, "data-testid" = id)
}

#' @title Text input for directory path entry in test mode
#' @description Returns a \code{shiny::textInput()} pre-configured for directory
#'   path entry with a \code{data-testid} attribute for test automation.
#' @param ns Namespace function from \code{shiny::NS()}.
#' @param input_id Character scalar — the input element ID (will be namespaced).
#' @param label Character scalar — the label displayed above the input.
#' @param value Character scalar — the default/initial value.
#' @return A Shiny tag (textInput) with a data-testid attribute.
#' @export
test_mode_dir_input <- function(ns, input_id, label, value) {
  testid(
    shiny::textInput(
      inputId = ns(input_id)
      , label = label
      , value = value
    )
    , id = input_id
  )
}

#' @title Hidden digest output for test state verification
#' @description Returns a hidden \code{verbatimTextOutput} element that exposes
#'   the current application state digest for automated test assertions. Returns
#'   NULL when test mode is disabled. When \code{ns} is NULL an identity
#'   function is used, supporting app-level usage outside a module namespace.
#' @param ns Namespace function, or NULL for app-level (no namespacing) usage.
#' @return A hidden Shiny tag, or NULL if test mode is inactive.
#' @export
test_mode_digest_ui <- function(ns = NULL) {
  if (!is_test_mode()) return(NULL)
  if (is.null(ns)) ns <- function(x) x
  shinyjs::hidden(
    shiny::div(
      shiny::verbatimTextOutput(ns("test_state_digest"))
    )
  )
}

#' @title Collect a structured digest of current application state
#' @description Assembles a named list snapshot of the app state for use in
#'   test assertions. Supports both the live reactiveValues structure returned
#'   by module servers (new) and plain-list mocks used in unit tests (legacy).
#'   All field reads are wrapped in \code{tryCatch}/\code{isolate}; missing
#'   data returns \code{NULL}.
#' @param values A list or reactiveValues-like object containing the central
#'   app state (e.g., \code{selected_omics}, \code{initialized_omics},
#'   \code{project_dirs}, \code{experiment_label}, \code{report_files}).
#' @param workflow_states A named list of per-omic workflow state objects.
#'   Each element may be a \code{reactiveValues} returned by a module server
#'   (with \code{state_manager}, \code{tab_status}, \code{active_tab}, and
#'   \code{processing_log} slots) or a plain list with legacy keys
#'   (\code{workflow_type}, \code{current_state}, \code{states}).
#' @return Named list with the following elements:
#'   \describe{
#'     \item{selected_omics}{Character vector of selected omics layers.}
#'     \item{initialized_omics}{Character vector of initialized omics modules.}
#'     \item{project_dir_keys}{Character vector of non-NULL project directory names.}
#'     \item{experiment_label}{Character scalar or NULL.}
#'     \item{workflow_type_per_omic}{Named list of workflow type strings per omic, or NULL.}
#'     \item{step_status_per_omic}{Named list of tab status objects per omic.}
#'     \item{r6_current_state_per_omic}{Named list of current WorkflowState state names per omic.}
#'     \item{r6_state_history_per_omic}{Named list of WorkflowState history vectors per omic.}
#'     \item{peptide_qc_audit_per_omic}{Named bounded summaries of persisted
#'       peptide-QC audit records per omic.}
#'     \item{active_tab_per_omic}{Named list of active tab strings per omic, or NULL.}
#'     \item{export_paths}{Named list of processing logs per omic.}
#'     \item{report_fingerprints}{Named list of MD5 hashes or NULL per report file.}
#'   }
#' @export
collect_state_digest <- function(values = list(), workflow_states = list()) {
  selected_omics <- values[["selected_omics"]] %||% character()
  initialized_omics <- values[["initialized_omics"]] %||% character()

  project_dirs <- values[["project_dirs"]]
  project_dir_keys <- if (is.null(project_dirs)) {
    character()
  } else {
    names(Filter(Negate(is.null), as.list(project_dirs)))
  }

  experiment_label <- values[["experiment_label"]]

  # workflow_type_per_omic:
  # Try state_manager$workflow_type (new module structure) then direct key (legacy)
  workflow_type_per_omic <- lapply(workflow_states, \(ws) {
    tryCatch(
      shiny::isolate({
        sm <- ws$state_manager
        if (!is.null(sm)) {
          wt <- tryCatch(sm$workflow_type, error = function(e) NULL)
          if (!is.null(wt)) return(wt)
        }
        ws[["workflow_type"]]
      })
      , error = function(e) NULL
    )
  })

  # step_status_per_omic:
  # Try tab_status list (new) then state names from states list (legacy)
  step_status_per_omic <- lapply(workflow_states, \(ws) {
    tryCatch(
      shiny::isolate({
        ts <- ws$tab_status
        if (!is.null(ts)) return(ts)
        names(ws[["states"]] %||% list())
      })
      , error = function(e) NULL
    )
  })

  r6_current_state_per_omic <- lapply(workflow_states, \(ws) {
    tryCatch(
      shiny::isolate({
        sm <- ws$state_manager
        if (!is.null(sm)) {
          current_state <- tryCatch(sm$current_state, error = function(e) NULL)
          if (!is.null(current_state)) return(current_state)
        }
        ws[["r6_current_state_name"]] %||% ws[["current_state"]]
      })
      , error = function(e) NULL
    )
  })

  r6_state_history_per_omic <- lapply(workflow_states, \(ws) {
    tryCatch(
      shiny::isolate({
        sm <- ws$state_manager
        if (!is.null(sm)) {
          history <- tryCatch(sm$getHistory(), error = function(e) NULL)
          if (is.null(history)) {
            history <- tryCatch(unlist(sm$state_history), error = function(e) NULL)
          }
          if (!is.null(history)) return(as.character(history))

          state_names <- tryCatch(names(sm$states), error = function(e) NULL)
          if (!is.null(state_names)) return(as.character(state_names))
        }
        history <- ws[["r6_state_history"]]
        if (!is.null(history)) return(as.character(history))
        names(ws[["states"]] %||% list())
      })
      , error = function(e) NULL
    )
  })

  peptide_qc_audit_per_omic <- lapply(workflow_states, \(ws) {
    tryCatch(
      shiny::isolate({
        state_manager <- ws$state_manager
        records <- if (!is.null(state_manager)) {
          state_manager$audit_records
        } else {
          ws[["audit_records"]]
        }
        records <- records %||% list()
        if (!length(records)) {
          return(list(
            status = "not_recorded",
            record_count = 0L,
            record_ids = character(),
            stage_ids = character(),
            canonical_digests = character(),
            immutable_import_digests = character(),
            all_records_complete = FALSE
          ))
        }

        record_field <- function(record, field) {
          value <- record[[field]]
          if (is.null(value) || length(value) != 1L || is.na(value)) {
            return(NA_character_)
          }
          as.character(value)
        }
        record_ids <- vapply(records, record_field, character(1L), field = "record_id")
        stage_ids <- vapply(records, record_field, character(1L), field = "stage_id")
        canonical_digests <- vapply(
          records,
          record_field,
          character(1L),
          field = "canonical_digest"
        )
        immutable_import_digests <- unique(vapply(
          records,
          record_field,
          character(1L),
          field = "immutable_import_digest"
        ))
        immutable_import_digests <- immutable_import_digests[
          !is.na(immutable_import_digests) & nzchar(immutable_import_digests)
        ]
        all_records_complete <- all(
          !is.na(record_ids) & nzchar(record_ids) &
            !is.na(stage_ids) & nzchar(stage_ids) &
            !is.na(canonical_digests) & grepl("^[0-9a-f]{64}$", canonical_digests)
        ) && length(immutable_import_digests) == 1L &&
          grepl("^[0-9a-f]{64}$", immutable_import_digests[[1L]])

        list(
          status = "recorded",
          record_count = length(records),
          record_ids = unname(record_ids),
          stage_ids = unname(stage_ids),
          canonical_digests = unname(canonical_digests),
          immutable_import_digests = unname(immutable_import_digests),
          all_records_complete = isTRUE(all_records_complete)
        )
      }),
      error = function(condition) list(
        status = "unavailable",
        record_count = 0L,
        record_ids = character(),
        stage_ids = character(),
        canonical_digests = character(),
        immutable_import_digests = character(),
        all_records_complete = FALSE
      )
    )
  })

  # active_tab_per_omic:
  # Try active_tab reactive (new, added by task-008) then current_state key (legacy)
  active_tab_per_omic <- lapply(workflow_states, \(ws) {
    tryCatch(
      shiny::isolate({
        at <- ws$active_tab %||% ws[["active_tab"]]
        if (is.function(at)) return(at())
        ws[["current_state"]]
      })
      , error = function(e) NULL
    )
  })

  # export_paths: processing_log per omic from workflow state
  export_paths <- lapply(workflow_states, \(ws) {
    tryCatch(
      shiny::isolate(ws$processing_log %||% list())
      , error = function(e) list()
    )
  })

  # report_fingerprints: md5sum of report files, NULL for missing/unreadable files
  report_files <- values[["report_files"]] %||% character()
  report_fingerprints <- lapply(report_files, \(path) {
    result <- tryCatch(
      tools::md5sum(path)
      , error = function(e) NULL
      , warning = function(w) NULL
    )
    if (is.null(result) || (length(result) == 1L && is.na(result))) NULL else result
  })
  if (length(report_files) > 0) {
    names(report_fingerprints) <- report_files
  }

  list(
    selected_omics = selected_omics
    , initialized_omics = initialized_omics
    , project_dir_keys = project_dir_keys
    , experiment_label = experiment_label
    , workflow_type_per_omic = workflow_type_per_omic
    , step_status_per_omic = step_status_per_omic
    , r6_current_state_per_omic = r6_current_state_per_omic
    , r6_state_history_per_omic = r6_state_history_per_omic
    , peptide_qc_audit_per_omic = peptide_qc_audit_per_omic
    , active_tab_per_omic = active_tab_per_omic
    , export_paths = export_paths
    , report_fingerprints = report_fingerprints
  )
}
