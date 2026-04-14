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
# func_lipid_norm.R
# ============================================================================
# Purpose: Lipidomics normalization standalone functions (non-S4)
# 
# NOTE: S4 methods for lipidomics normalization (setMethod calls) remain in
# their source files (lipid_normalization.R, lipidVsSamplesS4Objects.R)
# to maintain proper class-method coupling.
#
# This file is a placeholder for any future standalone (non-S4) lipidomics
# normalization helper functions that may be extracted.
#
# S4 methods for lipidomics normalization:
# - logTransformAssays() -> lipid_normalization.R
# - normaliseUntransformedData() -> lipid_normalization.R
# - normaliseBetweenSamples() -> lipidVsSamplesS4Objects.R
# - getNegCtrlMetabAnova() -> lipidVsSamplesS4Objects.R
#
# Dependencies:
# - limma
# - func_general_helpers.R (for utility functions)
# ============================================================================

# ============================================================================
# Helper Functions for Lipidomics Normalization Module
# ============================================================================













#' Build Normalization Configuration Object
#'
#' @title Build normalization configuration object
#' @param input Shiny input object
#' @return List of configuration parameters for state saving
#' @noRd
buildLipidNormConfig <- function(input) {
    list(
        itsd = list(
            applied = input$apply_itsd %||% FALSE
            , method = input$itsd_method %||% "median"
        )
        , log2 = list(
            offset = input$log_offset %||% 1
        )
        , normalization = list(
            method = input$norm_method %||% "cyclicloess"
        )
        , ruv = list(
            mode = input$ruv_mode %||% "skip"
            , grouping_variable = input$ruv_grouping_variable
            , auto_percentage_min = input$auto_percentage_min
            , auto_percentage_max = input$auto_percentage_max
            , separation_metric = input$separation_metric
            , k_penalty_weight = input$k_penalty_weight
            , adaptive_k_penalty = input$adaptive_k_penalty
            , manual_k = input$ruv_k
            , manual_percentage = input$ruv_percentage
        )
        , timestamp = Sys.time()
    )
}
