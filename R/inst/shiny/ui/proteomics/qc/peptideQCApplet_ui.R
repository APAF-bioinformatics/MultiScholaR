#' @title peptideQCAppletUI
#'
#' @description UI component for peptide-level quality control and filtering.
#' Extracted from qualityControlApplet.R to enable format-specific workflows.
#' This module handles all peptide filtering steps for DIA workflows.
#'
#' @param id Module ID
#' @export
#' @import shiny
#' @import shinydashboard
peptideQCAppletUI <- function(id) {
  ns <- NS(id)
  
  # Nested tabs for each peptide filtering step (extracted from qualityControlApplet.R)
  shiny::tabsetPanel(
    id = ns("peptide_filter_tabs"),
    
    # Step 1: Q-Value Filter (chunk 10)
    shiny::tabPanel(
      "Q-Value Filter",
      shiny::br(),
      shiny::fluidRow(
        shiny::column(4,
          shiny::wellPanel(
            shiny::h4("Q-Value & Proteotypic Peptide Filter"),
            shiny::p("Filter peptides based on statistical confidence (q-value) and ensure only proteotypic peptides are retained."),
            shiny::hr(),
            
            shiny::numericInput(ns("qvalue_threshold"), 
              "Q-Value Threshold", 
              value = 0.01, min = 0.001, max = 0.1, step = 0.001,
              width = "100%"
            ),
            shiny::helpText("Lower = stricter peptide ID confidence (default: 0.01)"),
            
            shiny::numericInput(ns("global_qvalue_threshold"), 
              "Global Q-Value Threshold", 
              value = 0.01, min = 0.001, max = 0.1, step = 0.001,
              width = "100%"
            ),
            shiny::helpText("Lower = stricter protein group ID confidence (default: 0.01)"),
            
            shiny::checkboxInput(ns("proteotypic_only"), 
              "Keep only proteotypic peptides", 
              value = TRUE
            ),
            shiny::helpText("TRUE = unique peptides only, FALSE = allow shared peptides (default: TRUE)"),
            
            shiny::hr(),
            shiny::div(
              shiny::actionButton(ns("apply_qvalue_filter"), "Apply Filter", 
                class = "btn-primary", width = "48%"),
              shiny::actionButton(ns("revert_qvalue"), "Revert", 
                class = "btn-warning", width = "48%", style = "margin-left: 4%")
            )
          )
        ),
        shiny::column(8,
          shiny::verbatimTextOutput(ns("qvalue_results")),
          shiny::br(),
          shinyjqui::jqui_resizable(
            shiny::plotOutput(ns("qvalue_plot"), height = "800px", width = "100%")
          )
        )
      )
    ),
    
    # Step 2: Precursor Rollup (chunk 11)
    shiny::tabPanel(
      "Precursor Rollup",
      shiny::br(),
      shiny::fluidRow(
        shiny::column(4,
          shiny::wellPanel(
            shiny::h4("Precursor to Peptide Rollup"),
            shiny::p("Aggregate intensity measurements from multiple precursor ions to peptide level."),
            shiny::hr(),
            shiny::p("This step has no configurable parameters - it uses statistical methods to combine precursor measurements."),
            shiny::hr(),
            shiny::div(
              shiny::actionButton(ns("apply_rollup"), "Apply Rollup", 
                class = "btn-primary", width = "48%"),
              shiny::actionButton(ns("revert_rollup"), "Revert", 
                class = "btn-warning", width = "48%", style = "margin-left: 4%")
            )
          )
        ),
        shiny::column(8,
          shiny::verbatimTextOutput(ns("rollup_results")),
          shiny::br(),
          shinyjqui::jqui_resizable(
            shiny::plotOutput(ns("rollup_plot"), height = "800px", width = "100%")
          )
        )
      )
    ),
    
    # Step 3: Intensity Filter (chunk 12)
    shiny::tabPanel(
      "Intensity Filter",
      shiny::br(),
      shiny::fluidRow(
        shiny::column(4,
          shiny::wellPanel(
            shiny::h4("Peptide Intensity Filter"),
            shiny::p("Remove peptides with low intensity values or detected in too few samples."),
            shiny::hr(),
            
            shiny::numericInput(ns("intensity_cutoff_percentile"), 
              "Intensity Cutoff Percentile (%)", 
              value = 1, min = 0.1, max = 10, step = 0.1,
              width = "100%"
            ),
            shiny::helpText("Higher = stricter intensity threshold (default: 1%)"),
            
            shiny::numericInput(ns("proportion_samples_below_cutoff"), 
              "Max Proportion Below Cutoff", 
              value = 0.5, min = 0.1, max = 1.0, step = 0.1,
              width = "100%"
            ),
            shiny::helpText("Higher = more relaxed filtering (default: 0.5)"),
            
            shiny::hr(),
            shiny::div(
              shiny::actionButton(ns("apply_intensity_filter"), "Apply Filter", 
                class = "btn-primary", width = "48%"),
              shiny::actionButton(ns("revert_intensity"), "Revert", 
                class = "btn-warning", width = "48%", style = "margin-left: 4%")
            )
          )
        ),
        shiny::column(8,
          shiny::verbatimTextOutput(ns("intensity_results")),
          shiny::br(),
          shinyjqui::jqui_resizable(
            shiny::plotOutput(ns("intensity_plot"), height = "800px", width = "100%")
          )
        )
      )
    ),
    
    # Step 4: Protein Peptide Count Filter (chunk 13)
    shiny::tabPanel(
      "Protein Peptides",
      shiny::br(),
      shiny::fluidRow(
        shiny::column(4,
          shiny::wellPanel(
            shiny::h4("Minimum Peptides per Protein"),
            shiny::p("Keep proteins only if they have sufficient peptide evidence (two-peptide rule)."),
            shiny::hr(),
            
            shiny::numericInput(ns("min_peptides_per_protein"), 
              "Min Peptides per Protein", 
              value = 2, min = 1, max = 5, step = 1,
              width = "100%"
            ),
            shiny::helpText("Higher = requires more unique peptides per protein (default: 2)"),
            
            shiny::numericInput(ns("min_peptidoforms_per_protein"), 
              "Min Peptidoforms per Protein", 
              value = 2, min = 1, max = 5, step = 1,
              width = "100%"
            ),
            shiny::helpText("Higher = requires more peptide forms per protein (default: 2)"),
            
            shiny::hr(),
            shiny::div(
              shiny::actionButton(ns("apply_protein_peptide_filter"), "Apply Filter", 
                class = "btn-primary", width = "48%"),
              shiny::actionButton(ns("revert_protein_peptide"), "Revert", 
                class = "btn-warning", width = "48%", style = "margin-left: 4%")
            )
          )
        ),
        shiny::column(8,
          shiny::verbatimTextOutput(ns("protein_peptide_results")),
          shiny::br(),
          shinyjqui::jqui_resizable(
            shiny::plotOutput(ns("protein_peptide_plot"), height = "800px", width = "100%")
          )
        )
      )
    ),
    
    # Step 5: Sample Quality Filter (chunk 14)
    shiny::tabPanel(
      "Sample Quality",
      shiny::br(),
      shiny::fluidRow(
        shiny::column(4,
          shiny::wellPanel(
            shiny::h4("Minimum Peptides per Sample"),
            shiny::p("Remove samples with insufficient peptide counts (poor sample performance)."),
            shiny::hr(),
            
            shiny::numericInput(ns("min_peptides_per_sample"), 
              "Min Peptides per Sample", 
              value = 500, min = 100, max = 5000, step = 100,
              width = "100%"
            ),
            shiny::helpText("Higher = stricter sample quality filter (default: 500)"),
            
            shiny::hr(),
            shiny::div(
              shiny::actionButton(ns("apply_sample_filter"), "Apply Filter", 
                class = "btn-primary", width = "48%"),
              shiny::actionButton(ns("revert_sample"), "Revert", 
                class = "btn-warning", width = "48%", style = "margin-left: 4%")
            )
          )
        ),
        shiny::column(8,
          shiny::verbatimTextOutput(ns("sample_results")),
          shiny::br(),
          shinyjqui::jqui_resizable(
            shiny::plotOutput(ns("sample_plot"), height = "800px", width = "100%")
          )
        )
      )
    ),
    
    # Step 6: Replicate Filter (chunk 15)
    shiny::tabPanel(
      "Replicate Filter",
      shiny::br(),
      shiny::fluidRow(
        shiny::column(4,
          shiny::wellPanel(
            shiny::h4("Remove Single-Replicate Peptides"),
            shiny::p("Remove peptides that appear in only one replicate across all groups."),
            shiny::hr(),
            
            shiny::textInput(ns("replicate_group_column"), 
              "Replicate Group Column", 
              value = "replicates",
              width = "100%"
            ),
            shiny::helpText("Column name for grouping replicates (default: 'replicates')"),
            
            shiny::hr(),
            shiny::div(
              shiny::actionButton(ns("apply_replicate_filter"), "Apply Filter", 
                class = "btn-primary", width = "48%"),
              shiny::actionButton(ns("revert_replicate"), "Revert", 
                class = "btn-warning", width = "48%", style = "margin-left: 4%")
            )
          )
        ),
        shiny::column(8,
          shiny::verbatimTextOutput(ns("replicate_results")),
          shiny::br(),
          shinyjqui::jqui_resizable(
            shiny::plotOutput(ns("replicate_plot"), height = "800px", width = "100%")
          )
        )
      )
    ),
    
    # Step 7: Missing Value Imputation (chunk 16)
    shiny::tabPanel(
      "Imputation",
      shiny::br(),
      shiny::fluidRow(
        shiny::column(4,
          shiny::wellPanel(
            shiny::h4("Missing Value Imputation"),
            shiny::p("Impute missing values using technical replicate averages."),
            shiny::hr(),
            
            shiny::numericInput(ns("proportion_missing_values"), 
              "Max Proportion Missing", 
              value = 0.5, min = 0.1, max = 1.0, step = 0.1,
              width = "100%"
            ),
            shiny::helpText("Lower = more stringent (less imputation), higher = less stringent (more imputation) (default: 0.5)"),
            
            shiny::hr(),
            shiny::div(
              shiny::actionButton(ns("apply_imputation"), "Apply Imputation", 
                class = "btn-primary", width = "48%"),
              shiny::actionButton(ns("revert_imputation"), "Revert", 
                class = "btn-warning", width = "48%", style = "margin-left: 4%")
            )
          )
        ),
        shiny::column(8,
          shiny::verbatimTextOutput(ns("imputation_results")),
          shiny::br(),
          shinyjqui::jqui_resizable(
            shiny::plotOutput(ns("imputation_plot"), height = "800px", width = "100%")
          )
        )
      )
    )
  )
} 