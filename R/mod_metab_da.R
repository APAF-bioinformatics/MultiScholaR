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

# ============================================================================
# mod_metab_da.R
# ============================================================================
# Purpose: Metabolomics differential abundance analysis Shiny module
#
# This module provides per-assay differential abundance analysis using limma,
# with interactive volcano plots (Glimma), heatmaps, and result tables.
# UI and functionality match the proteomics DA UI and server modules.
# KIT-106 W4: the tab is the shared DA frame (mod_omic_da.R, mod_omic_da_panels.R)
# under metabolomics' `da` configuration; this file keeps the documentation block
# and the UI wrapper.
# ============================================================================

#' @title Metabolomics Differential Analysis Module
#' @description A Shiny module for performing per-assay differential abundance
#'              analysis on metabolomics data using limma.
#' @name mod_metab_da
NULL

#' @rdname mod_metab_da
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h3 h4 h5 p hr br selectInput numericInput actionButton uiOutput verbatimTextOutput plotOutput tags downloadButton checkboxInput textAreaInput tabsetPanel tabPanel conditionalPanel helpText
#' @importFrom DT DTOutput
#' @importFrom shinyjqui jqui_resizable
mod_metab_da_ui <- function(id) {
    # KIT-106 W4: the shared DA frame (mod_omic_da.R) under metabolomics' `da` configuration.
    mod_omic_da_ui(id, omicGuiConfig("metabolomics")$da)
}
