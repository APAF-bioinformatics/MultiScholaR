#' @title proteinQCAppletUI
#'
#' @description UI component for protein-level quality control and filtering.
#' Extracted from qualityControlApplet.R to enable format-specific workflows.
#' This module handles all protein filtering steps after peptide-to-protein rollup.
#'
#' @param id Module ID
#' @param workflow_type Character string indicating workflow type ("TMT", "DIA", "LFQ").
#'   If NULL or LFQ/DIA, shows IQ Rollup tab. If TMT, skips IQ Rollup tab.
#' @export
#' @import shiny
#' @import shinydashboard
proteinQCAppletUI <- function(id, workflow_type = NULL) {
  ns <- NS(id)
  
  # Build tab list conditionally based on workflow type
  tab_list <- list()
  
  # Step 1: IQ Protein Rollup (chunk 17) - ONLY for LFQ/DIA workflows
  if (is.null(workflow_type) || workflow_type %in% c("LFQ", "DIA")) {
    tab_list[[length(tab_list) + 1]] <- shiny::tabPanel(
      "IQ Protein Rollup",
      shiny::br(),
      shiny::fluidRow(
        shiny::column(4,
          shiny::wellPanel(
            shiny::h4("Peptide-to-Protein Rollup"),
            shiny::p("Aggregate peptide-level data to protein-level quantification using the IQ algorithm (MaxLFQ)."),
            shiny::hr(),
            shiny::p("This step uses the IQ tool to implement the MaxLFQ algorithm for protein quantification, then automatically creates a ProteinQuantitativeData S4 object."),
            shiny::hr(),
            shiny::div(
              shiny::actionButton(ns("apply_iq_rollup"), "Run IQ Rollup & Create S4 Object", 
                class = "btn-primary", width = "48%"),
              shiny::actionButton(ns("revert_iq_rollup"), "Revert", 
                class = "btn-warning", width = "48%", style = "margin-left: 4%")
            )
          )
        ),
        shiny::column(8,
          shiny::verbatimTextOutput(ns("iq_rollup_results")),
          shiny::br(),
          shinyjqui::jqui_resizable(
            shiny::plotOutput(ns("iq_rollup_plot"), height = "800px", width = "100%")
          )
        )
      )
    )
  }
  
  # Step 2-5: Common tabs for ALL workflows (including TMT)
  tab_list <- c(tab_list, list(
    # Step 2: Protein Accession Cleanup (chunk 19)
    shiny::tabPanel(
      "Accession Cleanup",
      shiny::br(),
      shiny::fluidRow(
        shiny::column(4,
          shiny::wellPanel(
            shiny::h4("Best Protein Accession Selection"),
            shiny::p("Choose the best protein accession for ambiguous mappings using FASTA metadata."),
            shiny::hr(),
            
            shiny::textInput(ns("delimiter"), 
              "Protein ID Delimiter", 
              value = ";",
              width = "100%"
            ),
            shiny::helpText("Character separating multiple protein IDs (default: ';')"),
            
            # Detected delimiters info box
            shiny::div(
              style = "background-color: #e7f3ff; padding: 10px; border-left: 4px solid #2196F3; border-radius: 4px; margin-top: 10px;",
              shiny::icon("info-circle", style = "color: #2196F3;"),
              shiny::strong(" Detected Delimiters in Your Data:"),
              shiny::br(),
              shiny::textOutput(ns("detected_delimiters"), inline = TRUE)
            ),
            
            # Sample protein IDs display
            shiny::div(
              style = "background-color: #f5f5f5; padding: 10px; border-left: 4px solid #9E9E9E; border-radius: 4px; margin-top: 10px; font-family: monospace; font-size: 11px;",
              shiny::icon("list", style = "color: #9E9E9E;"),
              shiny::strong(" Sample Protein IDs:"),
              shiny::br(),
              shiny::verbatimTextOutput(ns("sample_protein_ids"))
            ),
            
            shiny::selectInput(ns("aggregation_method"), 
              "Aggregation Method", 
              choices = c("sum", "mean", "median"),
              selected = "mean",
              width = "100%"
            ),
            shiny::helpText("Method for combining duplicate protein measurements (S4 default: mean)"),
            
            shiny::hr(),
            shiny::div(
              shiny::actionButton(ns("apply_accession_cleanup"), "Apply Cleanup", 
                class = "btn-primary", width = "48%"),
              shiny::actionButton(ns("revert_accession_cleanup"), "Revert", 
                class = "btn-warning", width = "48%", style = "margin-left: 4%")
            )
          )
        ),
        shiny::column(8,
          shiny::verbatimTextOutput(ns("accession_cleanup_results")),
          shiny::br(),
          shinyjqui::jqui_resizable(
            shiny::plotOutput(ns("accession_cleanup_plot"), height = "800px", width = "100%")
          )
        )
      )
    ),
    
    # Step 3: Protein Intensity Filter (chunks 20+21 combined)
    shiny::tabPanel(
      "Protein Intensity Filter",
      shiny::br(),
      shiny::fluidRow(
        shiny::column(4,
          shiny::wellPanel(
            shiny::h4("Missing Value Parameters & Intensity Filter"),
            shiny::p("Configure missing value thresholds and filter proteins on intensity and missing values."),
            shiny::hr(),
            
            # Strict Mode Toggle
            shiny::checkboxInput(ns("use_strict_mode"), 
              "Strict Mode (No Missing Values Allowed)", 
              value = FALSE,
              width = "100%"
            ),
            shiny::helpText("When enabled, removes any protein with missing values in ANY sample across all groups. Overrides flexible threshold settings below."),
            
            shiny::hr(),
            
            # Flexible mode parameters (shown when strict mode is OFF)
            shiny::conditionalPanel(
              condition = "!input.use_strict_mode",
              ns = ns,
              
              # Missing value parameters (chunk 20)
              shiny::h5("Missing Value Parameters (Flexible Mode)"),
              shiny::numericInput(ns("min_reps_per_group"), 
                "Min Replicates per Group", 
                value = 2, min = 1, max = 10, step = 1,
                width = "100%"
              ),
              shiny::helpText("Minimum valid measurements required within each group (default: 2)"),
              
              shiny::numericInput(ns("min_groups"), 
                "Min Groups Required", 
                value = 2, min = 1, max = 10, step = 1,
                width = "100%"
              ),
              shiny::helpText("Minimum groups that must meet the replicate threshold (default: 2)"),
              
              shiny::hr(),
              
              # Calculated percentages (read-only display)
              shiny::h5("Calculated Thresholds"),
              shiny::p(style = "background-color: #f0f0f0; padding: 10px; border-radius: 5px;",
                shiny::strong("These values are automatically calculated from your replicate settings above:"),
                shiny::br(),
                shiny::textOutput(ns("calculated_groupwise_percent"), inline = TRUE),
                shiny::br(),
                shiny::textOutput(ns("calculated_max_groups_percent"), inline = TRUE)
              ),
              
              shiny::hr()
            ),
            
            # Intensity cutoff (always shown)
            shiny::h5("Intensity Threshold"),
            shiny::numericInput(ns("proteins_intensity_cutoff_percentile"), 
              "Protein Intensity Cutoff Percentile (%)", 
              value = 1, min = 0.1, max = 10, step = 0.1,
              width = "100%"
            ),
            shiny::helpText("Intensity threshold percentile (default: 1%)"),
            
            shiny::hr(),
            shiny::div(
              shiny::actionButton(ns("apply_protein_intensity_filter"), "Apply Filter", 
                class = "btn-primary", width = "48%"),
              shiny::actionButton(ns("revert_protein_intensity_filter"), "Revert", 
                class = "btn-warning", width = "48%", style = "margin-left: 4%")
            )
          )
        ),
        shiny::column(8,
          shiny::verbatimTextOutput(ns("protein_intensity_filter_results")),
          shiny::br(),
          shinyjqui::jqui_resizable(
            shiny::plotOutput(ns("protein_intensity_filter_plot"), height = "800px", width = "100%")
          )
        )
      )
    ),
    
    # Step 4: Duplicate Protein Removal (chunk 22)
    shiny::tabPanel(
      "Duplicate Removal",
      shiny::br(),
      shiny::fluidRow(
        shiny::column(4,
          shiny::wellPanel(
            shiny::h4("Remove Duplicate Proteins"),
            shiny::p("Aggregate duplicate protein entries by taking the mean across matching proteins."),
            shiny::hr(),
            
            shiny::selectInput(ns("duplicate_aggregation_method"), 
              "Aggregation Method", 
              choices = c("mean", "median", "max"),
              selected = "mean",
              width = "100%"
            ),
            shiny::helpText("Method for combining duplicate protein measurements (default: mean)"),
            
            shiny::hr(),
            shiny::div(
              shiny::actionButton(ns("apply_duplicate_removal"), "Remove Duplicates", 
                class = "btn-primary", width = "48%"),
              shiny::actionButton(ns("revert_duplicate_removal"), "Revert", 
                class = "btn-warning", width = "48%", style = "margin-left: 4%")
            )
          )
        ),
        shiny::column(8,
          shiny::verbatimTextOutput(ns("duplicate_removal_results")),
          shiny::br(),
          shinyjqui::jqui_resizable(
            shiny::plotOutput(ns("duplicate_removal_plot"), height = "800px", width = "100%")
          )
        )
      )
    ),
    
    # Step 5: Protein Replicate Filter (chunk 23)
    shiny::tabPanel(
      "Protein Replicate Filter",
      shiny::br(),
      shiny::fluidRow(
        shiny::column(4,
          shiny::wellPanel(
            shiny::h4("Remove Single-Replicate Proteins"),
            shiny::p("Remove proteins detected in only one replicate across experimental groups."),
            shiny::hr(),
            
            shiny::textInput(ns("protein_grouping_variable"), 
              "Grouping Variable", 
              value = "group",
              width = "100%"
            ),
            shiny::helpText("Column name for experimental grouping (default: 'group')"),
            
            shiny::numericInput(ns("parallel_cores"), 
              "Parallel Processing Cores", 
              value = 4, min = 1, max = 8, step = 1,
              width = "100%"
            ),
            shiny::helpText("Number of cores for parallel processing (default: 4)"),
            
            shiny::hr(),
            shiny::div(
              shiny::actionButton(ns("apply_protein_replicate_filter"), "Apply Filter", 
                class = "btn-primary", width = "48%"),
              shiny::actionButton(ns("revert_protein_replicate_filter"), "Revert", 
                class = "btn-warning", width = "48%", style = "margin-left: 4%")
            )
          )
        ),
        shiny::column(8,
          shiny::verbatimTextOutput(ns("protein_replicate_filter_results")),
          shiny::br(),
          shinyjqui::jqui_resizable(
            shiny::plotOutput(ns("protein_replicate_filter_plot"), height = "800px", width = "100%")
          )
        )
      )
    )
  ))
  
  # Construct tabsetPanel with conditional tab list
  do.call(shiny::tabsetPanel, c(list(id = ns("protein_filter_tabs")), tab_list))
} 