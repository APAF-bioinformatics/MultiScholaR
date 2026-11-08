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
# mod_lipid_da.R
# ============================================================================
# Purpose: Lipidomics differential abundance analysis Shiny module
#
# This module provides per-assay differential abundance analysis using limma,
# with interactive volcano plots (Glimma), heatmaps, and result tables.
# UI and functionality match the proteomics DA UI and server modules.
# KIT-106 W3: the tab is the shared DA frame (mod_omic_da.R, mod_omic_da_panels.R)
# under lipidomics' `da` configuration; this file keeps the documentation block.
# ============================================================================

#' @title Lipidomics Differential Analysis Module
#' @description A Shiny module for performing per-assay differential abundance
#'              analysis on lipidomics data using limma.
#' @name mod_lipid_da
NULL
