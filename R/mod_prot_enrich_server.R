#' @rdname enrichmentAnalysisAppletModule
#' @export
#' @importFrom shiny moduleServer reactive reactiveValues observeEvent req showNotification
#' @importFrom purrr detect
#' @importFrom readr write_tsv
#' @importFrom dplyr filter select mutate case_when slice_head all_of pull
#' @importFrom stringr str_split str_detect str_extract
#' @param id,workflow_data,experiment_paths,omic_type,experiment_label,selected_tab Runtime inputs used by this function; see the usage section for accepted values.
mod_prot_enrich_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, selected_tab = NULL) {
  # KIT-107 W3 (D2): the shared enrichment tab (mod_omic_enrich.R) with proteomics' `enrichment` group; it returns the
  # tab's reactive values, as before.
  mod_omic_enrich_server(id, workflow_data, experiment_paths, omic_type, experiment_label,
    config = omicGuiConfig("proteomics")$enrichment, selected_tab = selected_tab)
}
