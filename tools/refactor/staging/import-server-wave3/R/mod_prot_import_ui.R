#' Setup & Import Applet UI
#' 
#' UI for data import and initial setup with multi-format support
#' 
#' @param id Module ID
#' 
#' @importFrom shiny NS tagList fileInput numericInput textInput actionButton uiOutput
#' @importFrom shiny fluidRow column wellPanel h3 h4 tags hr radioButtons verbatimTextOutput helpText div br conditionalPanel icon
#' @importFrom shinyjs useShinyjs disable enable
#' @export
mod_prot_import_ui <- function(id) {
  ns <- shiny::NS(id)
  
  # Check if shinyFiles is available
  use_shiny_files <- requireNamespace("shinyFiles", quietly = TRUE)
  
  shiny::tagList(
    shinyjs::useShinyjs(),  # Enable shinyjs for dynamic UI control
    shiny::fluidRow(
    shiny::column(12,
      shiny::wellPanel(
        shiny::h3("Setup & Import Data"),
        
        # Data import section
        shiny::fluidRow(
          shiny::column(6,
            shiny::h4("Step 1: Import Searched Data File"),
            if (use_shiny_files) {
              shiny::tagList(
                shinyFiles::shinyFilesButton(
                  ns("search_results"), 
                  "Select proteomics search results file",
                  "Select File",
                  multiple = FALSE,
                  icon = shiny::icon("file")
                ),
                shiny::br(),
                shiny::verbatimTextOutput(ns("search_results_path"), placeholder = TRUE)
              )
            } else {
              shiny::fileInput(ns("search_results_standard"), 
                       "Select proteomics search results file:",
                       accept = c(".tsv", ".txt", ".tab", ".csv", ".xlsx", ".zip", ".parquet"))
            },
            
            # Format detection output
            shiny::uiOutput(ns("format_detection")),
            
            # Manual format override option
            shiny::radioButtons(ns("format_override"),
                        "Override detected format:",
                        choices = c("Auto-detect" = "auto",
                                  "DIA-NN" = "diann",
                                  "Spectronaut DIA" = "spectronaut",
                                  "FragPipe LFQ" = "fragpipe",
                                  "MaxQuant LFQ" = "maxquant",
                                  "Proteome Discoverer TMT" = "pd_tmt"),
                        selected = "auto",
                        inline = TRUE),
            
            shiny::h4("Step 2: Import FASTA File"),
            if (use_shiny_files) {
              shiny::tagList(
                shinyFiles::shinyFilesButton(
                  ns("fasta_file"), 
                  "Select FASTA file",
                  "Select File",
                  multiple = FALSE,
                  icon = shiny::icon("file")
                ),
                shiny::br(),
                shiny::verbatimTextOutput(ns("fasta_file_path"), placeholder = TRUE)
              )
            } else {
              shiny::fileInput(ns("fasta_file_standard"), 
                       "Select FASTA file:",
                       accept = c(".fasta", ".fa", ".faa"))
            },
            
            # Sanitization and FASTA options
            shiny::wellPanel(
              shiny::checkboxInput(ns("sanitize_names"),
                           "Sanitize Sample Names",
                           value = TRUE),
              shiny::helpText("Clean sample IDs (e.g., '123-Sample!' -> 'x123_sample') for better compatibility with downstream analysis."),
              
              shiny::checkboxInput(ns("mixed_species_fasta"),
                           "Mixed species FASTA (analyze organism distribution)",
                           value = FALSE),
              shiny::helpText("Check this if your FASTA contains proteins from multiple species (e.g., spiked-in standards, contaminants)")
            ),
            
            shiny::h4("Step 3: Organism Information"),
            shiny::div(
              id = ns("organism_inputs_container"),
              shiny::numericInput(ns("taxon_id"), 
                          "Taxonomy ID:", 
                          value = 9606,  # Human default
                          min = 1),
              shiny::textInput(ns("organism_name"), 
                       "Organism Name:", 
                       value = "Homo sapiens")
            ),
            # Help text shown when mixed species is checked
            shiny::conditionalPanel(
              condition = paste0("input['", ns("mixed_species_fasta"), "'] == true"),
              shiny::helpText(
                shiny::icon("info-circle"),
                " Organism will be selected from the analysis popup after processing."
              ),
              style = "color: #666; font-style: italic;"
            )
          ),
          
          shiny::column(6,
            shiny::h4("Optional: ID Mapping Files"),
            if (use_shiny_files) {
              shiny::tagList(
                shinyFiles::shinyFilesButton(
                  ns("uniprot_mapping"), 
                  "UniProt ID mapping file (optional)",
                  "Select File",
                  multiple = FALSE,
                  icon = shiny::icon("file")
                ),
                shiny::verbatimTextOutput(ns("uniprot_mapping_path"), placeholder = TRUE),
                shiny::br(),
                
                shinyFiles::shinyFilesButton(
                  ns("uniparc_mapping"), 
                  "UniParc ID mapping file (optional)",
                  "Select File",
                  multiple = FALSE,
                  icon = shiny::icon("file")
                ),
                shiny::verbatimTextOutput(ns("uniparc_mapping_path"), placeholder = TRUE)
              )
            } else {
              shiny::tagList(
                shiny::fileInput(ns("uniprot_mapping_standard"), 
                         "UniProt ID mapping file (optional):",
                         accept = c(".tsv", ".txt", ".tab")),
                shiny::fileInput(ns("uniparc_mapping_standard"), 
                         "UniParc ID mapping file (optional):",
                         accept = c(".tsv", ".txt", ".tab"))
              )
            },
            
            shiny::h4("Configuration"),
            if (use_shiny_files) {
              shiny::tagList(
                shinyFiles::shinyFilesButton(
                  ns("config_file"), 
                  "Configuration file (optional)",
                  "Select File",
                  multiple = FALSE,
                  icon = shiny::icon("file")
                ),
                shiny::verbatimTextOutput(ns("config_file_path"), placeholder = TRUE)
              )
            } else {
              shiny::fileInput(ns("config_file_standard"), 
                       "Configuration file (optional):",
                       accept = c(".ini", ".cfg"))
            },
            shiny::tags$small("If not provided, default configuration will be used"),
            
            # --- TESTTHAT CHECKPOINT CAPTURE (TEMPORARY) ---
            shiny::hr(),
            shiny::checkboxInput(ns("capture_checkpoints"),
                         "Capture Test Checkpoints (Developer)",
                         value = getOption("multischolar.capture_test_checkpoints", FALSE)),
            shiny::helpText("Snapshots will be saved to tests/testdata/sepsis/proteomics/"),
            # --- END CHECKPOINT CAPTURE ---
            
            # Format-specific options will appear here
            shiny::uiOutput(ns("format_specific_options"))
          )
        ),
        
        shiny::hr(),
        
        # Process button and status
        shiny::fluidRow(
          shiny::column(12,
            shiny::actionButton(ns("process_data"), 
                        "Process Imported Data", 
                        class = "btn-success",
                        width = "100%"),
            shiny::br(),
            shiny::br(),
            shiny::uiOutput(ns("import_status")),
            shiny::uiOutput(ns("data_summary"))
          )
        )
      )
    )
  )
  )
}

