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

# KIT-103 W5: the metab design builder's exported names, thin wrappers over the shared builder
# (mod_omic_design_builder.R) under the metab configuration. The builder server's closure
# (mod_metab_design_builder_server.R) is deleted; its wrapper lives here, beside the UI's.

#' @title Metabolomics Design Matrix Builder Module
#'
#' @description A Shiny module to build a design matrix for metabolomics experiments.
#' This module follows the proteomics design-builder UI and server modules,
#' with the key difference that metabolomics data is stored as a LIST of assays
#' (e.g., LCMS_Pos, LCMS_Neg) where samples are COLUMN NAMES, not row values.
#'
#' @param id The module ID.
#' @param data_tbl A reactive expression that returns the data table LIST
#'   (e.g., list(LCMS_Pos = df_pos, LCMS_Neg = df_neg)).
#' @param config_list A reactive expression that returns the main configuration list.
#' @param column_mapping A reactive expression that returns the column mapping list.
#'
#' @return A reactive expression that returns a list containing the results
#'   when the "Save" button is clicked. The list includes:
#'   - `design_matrix`: The final, filtered design matrix.
#'   - `data_cln`: The data table LIST, filtered to include only assigned samples.
#'   - `contrasts_tbl`: A data frame of the defined contrasts.
#'   - `config_list`: The updated configuration list with the new formula.
#'
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @import shinydashboard
#' @import DT
#' @rawNamespace import(gtools, except = chr)
#' @import dplyr
#' @importFrom tibble tibble
#' @name metabolomicsDesignMatrixBuilderModule
NULL

#' Metabolomics Design Matrix Builder UI Module
#'
#' @param id Module ID.
#' @return A fluidPage containing the UI for the design matrix builder.
#'
#' @note **ARCHITECTURAL EXCEPTION**: This module is designed as a standalone component
#' callable from R Markdown workflows, not integrated into the main Shiny app workflow.
#' Therefore, it returns `fluidPage()` instead of `tagList()`, and uses reactive expressions
#' instead of `workflow_data` for state management. This is an intentional deviation from
#' the standard Golem module pattern.
#'
#' @rdname metabolomicsDesignMatrixBuilderModule
#' @export
mod_metab_design_builder_ui <- function(id) {
    # KIT-103 W5: the shared builder UI (mod_omic_design_builder.R) with the metab texts.
    mod_omic_design_builder_ui(id, omicGuiConfig("metabolomics"))
}

#' Metabolomics Design Matrix Builder Server Module
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
#'
#' @rdname metabolomicsDesignMatrixBuilderModule
#' @export
mod_metab_design_builder_server <- function(id, data_tbl, config_list, column_mapping
                                            , existing_design_matrix = NULL
                                            , existing_contrasts = NULL
                                            , state_owner = NULL) {
    # KIT-103 W5: the shared builder (mod_omic_design_builder.R) under the metab configuration.
    mod_omic_design_builder_server(id, data_tbl, config_list, column_mapping,
        existing_design_matrix = existing_design_matrix, existing_contrasts = existing_contrasts,
        state_owner = state_owner, config = metabDesignGuiConfig())
}
