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

#' @title designMatrixBuilderModule
#'
#' @description A Shiny module to build a design matrix for proteomics experiments.
#' This module is a refactored version of the standalone design matrix applet,
#' designed to be embedded within the main MultiScholaR Shiny application.
#'
#' @param id The module ID.
#' @param data_tbl A reactive expression that returns the data table
#'   (e.g., from the setup & import tab).
#' @param config_list A reactive expression that returns the main configuration list.
#' @param column_mapping A reactive expression that returns the column mapping list.
#'
#' @return A reactive expression that returns a list containing the results
#'   when the "Save" button is clicked. The list includes:
#'   - `design_matrix`: The final, filtered design matrix.
#'   - `data_cln`: The data table, filtered to include only assigned runs.
#'   - `contrasts_tbl`: A data frame of the defined contrasts.
#'   - `config_list`: The updated configuration list with the new formula.
#'
#' @import shiny
#' @import shinydashboard
#' @import DT
#' @import gtools
#' @import dplyr
#' @importFrom tibble tibble
#' @name designMatrixBuilderModule
NULL





























































