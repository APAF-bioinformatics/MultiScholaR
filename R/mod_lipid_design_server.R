#' @rdname lipidomicsDesignMatrixAppletModule
#' @export
#' @importFrom shiny moduleServer reactive observeEvent req renderUI showNotification removeNotification outputOptions renderText
#' @importFrom logger log_info log_error log_warn
#' @importFrom utils write.table
#' @importFrom vroom vroom
#' @importFrom shinyFiles shinyDirButton shinyDirChoose parseDirPath getVolumes
#' @importFrom jsonlite write_json read_json
mod_lipid_design_server <- function(id, workflow_data, experiment_paths, volumes = NULL, qc_trigger = NULL, profile = NULL) {
    # KIT-103 W4: the shared design tab (mod_omic_design.R) under the lipid configuration, or the
    # profile's omic (trx) under its own.
    mod_omic_design_server(id, workflow_data, experiment_paths, volumes = volumes, qc_trigger = qc_trigger,
        config = lipidDesignGuiConfig(profile))
}

#' The design tab's commit of a builder result (the shared commit's GUI side,
#' `runOmicsDesignCommit()`) under the lipid configuration, or the profile's omic (trx). KIT-103 W4
#' keeps this name for its callers (MCI-019.6, MCI-027's helper, test-omics-design-core.R) until
#' KIT-103 W5 retargets them; `createLipidomicsAssayDataFn` stays in the signature for those callers
#' and is not used: the object is the configuration's `build_object`.
#' @noRd
runLipidDesignBuilderObserverShell <- function(results, workflowData, experimentPaths, qcTrigger = NULL,
                                               createLipidomicsAssayDataFn = NULL, profile = NULL) {
    runOmicsDesignCommit(results, workflowData, experimentPaths, qcTrigger, lipidDesignGuiConfig(profile))
}
