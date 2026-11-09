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

#' @title enrichmentAnalysisAppletModule
#'
#' @description A Shiny module for the Enrichment Analysis step of the proteomics
#' workflow. Handles GO, KEGG, Reactome, and String-DB pathway enrichment analysis
#' using results from differential expression analysis.
#'
#' @name enrichmentAnalysisAppletModule
NULL

#' @rdname enrichmentAnalysisAppletModule
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h4 selectInput helpText hr numericInput actionButton icon br tabsetPanel tabPanel verbatimTextOutput plotOutput downloadButton
#' @importFrom DT DTOutput renderDT
mod_prot_enrich_ui <- function(id) {
  # KIT-107 W3 (D2): the shared enrichment tab (mod_omic_enrich.R) with proteomics' `enrichment` group.
  mod_omic_enrich_ui(id, omicGuiConfig("proteomics")$enrichment)
}
