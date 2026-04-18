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

#' @title Lipidomics Design Matrix Builder Module
#'
#' @description A Shiny module to build a design matrix for lipidomics experiments.
#' This module is 1:1 with the proteomics design builder (mod_prot_design_builder.R),
#' with the key difference that lipidomics data is stored as a LIST of assays
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
#' @import shiny
#' @import shinydashboard
#' @import DT
#' @import gtools
#' @import dplyr
#' @importFrom tibble tibble
#' @name lipidomicsDesignMatrixBuilderModule
NULL





















































