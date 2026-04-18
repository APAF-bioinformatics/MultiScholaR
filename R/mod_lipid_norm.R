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
# mod_lipid_norm.R
# ============================================================================
# Purpose: Lipidomics normalization Shiny module (Proteomics UX harmonized)
#
# This module provides a unified normalization pipeline with per-assay
# visualization for lipidomics data (pos/neg mode support).
#
# Key features:
# - 3/9 column layout matching proteomics pattern
# - All options visible (single "Run" button)
# - Per-assay ITSD selection with DT::datatable
# - Independent auto-RUV optimization per assay
# - 3-column QC comparison (Post-Filter / Post-Norm / RUV-Corrected)
# - Per-assay visualization (pos row on top, neg row below)
# - Correlation filtering as final step
# ============================================================================

#' @title Lipidomics Normalization Module
#' @description A Shiny module for multi-step normalization of lipidomics data.
#'              Harmonized with proteomics UX pattern.
#' @name mod_lipid_norm
NULL


















































