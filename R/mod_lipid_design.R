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
# mod_lipid_design.R
# ============================================================================
# Purpose: Lipidomics design matrix host module - 1:1 with proteomics
#
# This module embeds mod_lipid_design_builder.R and handles orchestration,
# including Import Existing Design, S4 object creation, and state management.
# ============================================================================

#' @title Lipidomics Design Matrix Applet Module
#'
#' @description A Shiny module that serves as the main host for the design
#' matrix creation workflow step. Embeds mod_lipid_design_builder_ui/server
#' and handles Import Existing Design, S4 object creation, and state saving.
#' Architecture is 1:1 with mod_prot_design.R.
#'
#' @param id Module ID
#' @param workflow_data A reactive values object to store workflow data.
#' @param experiment_paths A list of paths for the current experiment.
#' @param volumes A list of volumes for shinyFiles (optional).
#' @param qc_trigger A reactive trigger for QC execution (optional).
#'
#' @name lipidomicsDesignMatrixAppletModule
NULL












