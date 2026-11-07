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
  # KIT-102 W6: the shared import frame (shinyjs, the well panel and its title) around the proteomics panel.
  mod_omic_import_ui(id, omicGuiConfig("proteomics"))
}

#' The proteomics import source panel's UI: the tags inside the frame's well panel, after its title (the frame's
#' panel contract, R/mod_omic_import.R). The choices, headings and test ids are `config$import`'s.
#' @param ns The import module's namespace function.
#' @param config `omicGuiConfig("proteomics")`.
#' @noRd
protImportPanelUi <- function(ns, config) {
  import <- config$import

  # Check if shinyFiles is available and we are not in test mode
  use_shiny_files <- requireNamespace("shinyFiles", quietly = TRUE) && !is_test_mode()

  shiny::tagList(
        # Data import section
        shiny::fluidRow(
          shiny::column(6,
            shiny::h4(import$help[["search"]]),
            shiny::div(
              `data-testid` = import$testids[["search_results_file"]],
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
              }
            ),
            
            # Format detection output
            shiny::uiOutput(ns("format_detection")),
            
            # Manual format override option
            shiny::radioButtons(ns("format_override"),
                        "Override detected format:",
                        choices = import$formats,
                        selected = import$default_format,
                        inline = TRUE),
            # KIT-049/KIT-046: the two DIA-NN workflows named at the entry point; the choice is the workflow type.
            shiny::radioButtons(ns("dia_route"),
                        "DIA-NN workflow (choose before import):",
                        choices = import$route_choice,
                        selected = import$default_route, inline = FALSE),
            shiny::helpText("Native DIA-IQ sessions keep this route fixed until a new project/import. Other formats are unchanged."),
            
            shiny::h4(import$help[["fasta"]]),
            shiny::div(
              `data-testid` = import$testids[["fasta_file"]],
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
              }
            ),
            
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
            
            shiny::h4(import$help[["organism"]]),
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
            testid(
              shiny::actionButton(ns("process_data"),
                          "Process Imported Data",
                          class = "btn-success",
                          width = "100%")
              , import$testids[["process"]]
            ),
            shiny::br(),
            shiny::br(),
            shiny::uiOutput(ns("import_status")),
            shiny::uiOutput(ns("data_summary"))
          )
        )
  )
}
