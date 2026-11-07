#' Lipidomics Design Matrix Builder Server Module
#'
#' @param id Module ID.
#' @param data_tbl A reactive expression that returns the data table LIST.
#' @param config_list A reactive expression that returns the main configuration list.
#' @param column_mapping A reactive expression that returns the column mapping list.
#' @param existing_design_matrix Optional reactive expression that returns an existing
#'   design matrix (e.g., from "Import Existing Design"). If provided and has valid
#'   group assignments, will be used instead of creating a fresh matrix.
#' @param existing_contrasts Optional reactive expression that returns an existing
#'   contrasts table (e.g., from "Import Existing Design").
#' @param state_owner Optional function returning the session's state owner; the builder
#'   then holds its tables by reference. Without it they stay plain values.
#' @param profile An omic profile (KIT-096): its run order for a fresh design, and the level order and removed runs on
#'   the saved result. `NULL`, the default, is lipidomics.
#'
#' @rdname lipidomicsDesignMatrixBuilderModule
#' @export
mod_lipid_design_builder_server <- function(id, data_tbl, config_list, column_mapping
                                            , existing_design_matrix = NULL
                                            , existing_contrasts = NULL
                                            , state_owner = NULL
                                            , profile = NULL) {
    # KIT-103 W4: the shared builder (mod_omic_design_builder.R) under the lipid configuration, or
    # the profile's omic (trx) under its own.
    mod_omic_design_builder_server(id, data_tbl, config_list, column_mapping,
        existing_design_matrix = existing_design_matrix, existing_contrasts = existing_contrasts,
        state_owner = state_owner, config = lipidDesignGuiConfig(profile))
}
