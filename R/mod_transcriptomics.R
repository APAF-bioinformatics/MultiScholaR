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

# KIT-096: the transcriptomics workflow module (import, design, QC). The tabs are the shared modules under
# `omicGuiConfig("transcriptomics")` (KIT-102 to KIT-108; KIT-109 removed the transcriptomics profile the lipid modules
# took); nothing is copied. The tabs hold state references only (KIT-075) and run on the managed runner, or
# the binding's wrappers in process (KIT-078). Recorded divergence from conventions/R-shiny.md: the shinydashboard and
# wellPanel shell, the managed runner in place of ExtendedTask and the import's fileInput fallback are kept, to match the
# other omics (the user's instruction) under the existing-code scope rule of conventions/R.md.

# ---------------------------------------------------------------------------------------------------------------
# The module.
# ---------------------------------------------------------------------------------------------------------------

#' @title Transcriptomics Workflow Module
#' @description The transcriptomics workflow's tabs: import, design and QC (KIT-096), normalisation, DA and the summary
#'   (KIT-097). KIT-102: the shared shell (mod_omic_workflow.R) under the transcriptomics configuration; its tabs are the
#'   shared modules under that configuration (the component registry, `omicGuiComponents`).
#' @param id Module id (`transcriptomics_workflow` in the app).
#' @name mod_transcriptomics
#' @noRd
mod_transcriptomics_ui <- function(id) {
    mod_omic_workflow_ui(id, omicGuiConfig("transcriptomics"))
}

#' @rdname mod_transcriptomics
#' @param project_dirs,omic_type,experiment_label,volumes As `mod_lipidomics_server`.
#' @return The module's bus.
#' @noRd
mod_transcriptomics_server <- function(id, project_dirs, omic_type, experiment_label, volumes = NULL) {
    mod_omic_workflow_server(id, project_dirs, omic_type, experiment_label, volumes, omicGuiConfig("transcriptomics"))
}
