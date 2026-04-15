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
  ns <- shiny::NS(id)
  
  shiny::tagList(
    shiny::wellPanel(
      shiny::fluidRow(
    shiny::column(3,
      shiny::wellPanel(
        shiny::h4("Enrichment Analysis Settings"),
        
        # Contrast selector
        shiny::selectInput(
          ns("selected_contrast"),
          "Select Contrast:",
          choices = NULL,
          width = "100%"
        ),
        shiny::helpText("Choose which DE contrast to analyze"),
        
        shiny::hr(),
        
        # Only the parameters that ACTUALLY EXIST in processEnrichments
        shiny::h5("Enrichment Settings"),
        shiny::numericInput(
          ns("up_cutoff"),
          "Up Log2FC Cutoff:",
          value = 0,
          min = 0,
          max = 5,
          step = 0.1
        ),
        shiny::helpText("Minimum log2 fold change for up-regulated genes"),
        
        shiny::numericInput(
          ns("down_cutoff"),
          "Down Log2FC Cutoff:",
          value = 0,
          min = 0,
          max = 5,
          step = 0.1
        ),
        shiny::helpText("Minimum log2 fold change for down-regulated genes"),
        
        shiny::numericInput(
          ns("q_cutoff"),
          "Q-value Cutoff:",
          value = 0.05,
          min = 0.001,
          max = 0.2,
          step = 0.005
        ),
        shiny::helpText("FDR threshold for enrichment significance"),
        
        shiny::textInput(
          ns("organism_taxid"),
          "Organism Taxon ID:",
          value = 9606,  # Default to human, will be updated by server
          placeholder = "e.g., 9606 for human"
        ),
        shiny::helpText("NCBI taxonomy ID for organism"),
        
        # [OK] NEW: Mixed species filtering checkbox
        shiny::hr(),
        shiny::h5("Multi-Species Filtering"),
        shiny::checkboxInput(
          ns("enable_organism_filter"),
          "Filter to target organism only",
          value = FALSE
        ),
        shiny::helpText("Enable if using mixed-species FASTA (e.g., with contaminants). Filters DE results to target organism before enrichment analysis."),
        
        # [OK] NEW: Analysis method display
        shiny::h5("Analysis Method"),
        shiny::verbatimTextOutput(ns("analysis_method_display")),
        shiny::helpText("Automatically determined based on organism support"),
        
        shiny::hr(),
        
        # Available contrasts display
        shiny::h5("Available Contrasts"),
        shiny::verbatimTextOutput(ns("contrasts_display")),
        shiny::helpText("Contrasts from DE analysis"),
        
        shiny::br(),
        
        # Main action button
        shiny::actionButton(
          ns("run_enrichment_analysis"),
          "Run Enrichment Analysis",
          class = "btn-primary",
          width = "100%",
          icon = shiny::icon("play")
        ),
        
        shiny::br(),
        shiny::br(),
        
        # Download button
        shiny::downloadButton(
          ns("download_enrichment_results"),
          "Download All Results",
          class = "btn-success",
          width = "100%"
        ),
        
        shiny::br(),
        shiny::br(),
        
        # Status display
        shiny::h5("Analysis Status"),
        shiny::verbatimTextOutput(ns("enrichment_status")),
        
        shiny::conditionalPanel(
          condition = sprintf("input['%s'] == 'gprofiler2'", ns("enrichment_method_tabs")),
          shiny::selectInput(
            ns("correction_method"),
            "Correction Method:",
            choices = list(
              "g:SCS (Conservative)" = "gSCS",
              "FDR (Less Conservative)" = "fdr"
            ),
            selected = "gSCS"
          ),
          shiny::helpText("g:SCS: More conservative, reduces false positives but may miss terms."),
          shiny::helpText("FDR: Detects more terms but higher false positive risk.")
        )
      )
    ),
    
    # [OK] FIXED: Split tab structure with plot displays added
    shiny::column(9,
      shiny::tabsetPanel(
        id = ns("enrichment_method_tabs"),
        
        # [OK] Tab 1: gprofiler2 for supported organisms
        shiny::tabPanel(
          "gprofiler2 Analysis",
          value = "gprofiler2",
          icon = shiny::icon("globe"),
          shiny::br(),
          
          shiny::div(
            style = "background-color: #e8f4fd; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
            shiny::h5("gprofiler2 Enrichment", style = "color: #2c5282; margin-top: 0;"),
            shiny::p("Comprehensive functional enrichment using gprofiler2 for well-supported model organisms. 
                     Includes GO terms, KEGG pathways, Reactome pathways, and more.", 
                     style = "margin-bottom: 0; color: #2c5282;")
          ),
          
          # gprofiler2 results display
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] == 'gprofiler2'", ns("enrichment_method_tabs")),
            
            # Direction filter for viewing results
            shiny::fluidRow(
              shiny::column(6,
                shiny::selectInput(
                  ns("gprofiler_direction_filter"),
                  "Direction:",
                  choices = list(
                    "All" = "all",
                    "Up-regulated only" = "up",
                    "Down-regulated only" = "down"
                  ),
                  selected = "all"
                )
              ),
              shiny::column(6,
                shiny::div(
                  style = "padding: 20px 0; text-align: center; color: #666;",
                  shiny::strong("Select direction to filter results")
                )
              )
            ),
            
            shiny::hr(),
            
            # [OK] NEW: Add plot display for gprofiler2
            shiny::fluidRow(
              shiny::column(12,
                shiny::h5("Enrichment Plot"),
                shinyjqui::jqui_resizable(
                  plotly::plotlyOutput(ns("gprofiler_plot"), height = "500px", width = "100%")
                ),
                shiny::br()
              )
            ),
            
            # gprofiler2 results table
            DT::DTOutput(ns("gprofiler_results_table")),
            
            shiny::br(),
            
            # gprofiler2 summary stats
            shiny::wellPanel(
              shiny::h5("gprofiler2 Enrichment Summary"),
              shiny::verbatimTextOutput(ns("gprofiler_summary_stats"))
            )
          )
        ),
        
        # [OK] Tab 2: clusterProfileR for unsupported organisms  
        shiny::tabPanel(
          "clusterProfileR Analysis",
          value = "clusterprofiler",
          icon = shiny::icon("project-diagram"),
          shiny::br(),
          
          shiny::div(
            style = "background-color: #f0fff4; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
            shiny::h5("clusterProfileR Enrichment", style = "color: #22543d; margin-top: 0;"),
            shiny::p("Custom enrichment analysis using UniProt GO annotations via clusterProfileR. 
                     Ideal for organisms not supported by gprofiler2.", 
                     style = "margin-bottom: 0; color: #22543d;")
          ),
          
          # clusterProfileR results display
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] == 'clusterprofiler'", ns("enrichment_method_tabs")),
            
            # Direction filter for viewing results
            shiny::fluidRow(
              shiny::column(6,
                shiny::selectInput(
                  ns("clusterprofiler_direction_filter"),
                  "Direction:",
                  choices = list(
                    "All" = "all",
                    "Up-regulated only" = "up",
                    "Down-regulated only" = "down"
                  ),
                  selected = "all"
                )
              ),
              shiny::column(6,
                shiny::div(
                  style = "padding: 20px 0; text-align: center; color: #666;",
                  shiny::strong("Select direction to filter results")
                )
              )
            ),
            
            shiny::hr(),
            
            # [OK] NEW: Add plot display for clusterProfileR
            shiny::fluidRow(
              shiny::column(12,
                shiny::h5("Enrichment Plot"),
                shinyjqui::jqui_resizable(
                  plotly::plotlyOutput(ns("clusterprofiler_plot"), height = "500px", width = "100%")
                ),
                shiny::br()
              )
            ),
            
            # clusterProfileR results table
            DT::DTOutput(ns("clusterprofiler_results_table")),
            
            shiny::br(),
            
            # clusterProfileR summary stats
            shiny::wellPanel(
              shiny::h5("clusterProfileR Enrichment Summary"),
              shiny::verbatimTextOutput(ns("clusterprofiler_summary_stats"))
            )
          )
        ),
        
        # [OK] Tab 3: STRING-DB for protein-protein interaction networks
        shiny::tabPanel(
          "STRING-DB Networks",
          value = "stringdb",
          icon = shiny::icon("network-wired"),
          shiny::br(),
          
          shiny::div(
            style = "background-color: #fef5e7; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
            shiny::h5("STRING-DB Network Analysis", style = "color: #744210; margin-top: 0;"),
            shiny::p("Protein-protein interaction network enrichment using STRING database. 
                     Identifies functional networks and interaction clusters.", 
                     style = "margin-bottom: 0; color: #744210;")
          ),
          
          # STRING-DB results display
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] == 'stringdb'", ns("enrichment_method_tabs")),
            
            # Filter controls for STRING-DB
            shiny::fluidRow(
              shiny::column(4,
                shiny::selectInput(
                  ns("stringdb_ranking_method"),
                  "Ranking Method:",
                  choices = list(
                    "Combined Score" = "combined_score",
                    "FDR Q-value" = "fdr_qvalue",
                    "Log2 Fold Change" = "log2fc"
                  ),
                  selected = "combined_score"
                )
              ),
              shiny::column(4,
                shiny::checkboxInput(
                  ns("stringdb_filter_significant"),
                  "Filter Significant Only",
                  value = FALSE
                )
              ),
              shiny::column(4,
                shiny::numericInput(
                  ns("stringdb_max_results"),
                  "Max Results:",
                  value = 50,
                  min = 10,
                  max = 500,
                  step = 10
                )
              )
            ),
            
            shiny::hr(),
            
            # [OK] NEW: Add plot display for STRING-DB
            shiny::fluidRow(
              shiny::column(12,
                shiny::h5("Network Plot"),
                shinyjqui::jqui_resizable(
                  plotly::plotlyOutput(ns("stringdb_plot"), height = "500px", width = "100%")
                ),
                shiny::br()
              )
            ),
            
            # STRING-DB results table
            DT::DTOutput(ns("stringdb_results_table")),
            
            shiny::br(),
            
            # STRING-DB summary stats  
            shiny::wellPanel(
              shiny::h5("STRING-DB Network Summary"),
              shiny::verbatimTextOutput(ns("stringdb_summary_stats"))
            )
          )
        )
      )
    )
  )
  )
  )
}






































