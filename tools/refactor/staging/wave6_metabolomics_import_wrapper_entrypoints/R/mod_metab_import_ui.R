#' @rdname mod_metab_import
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h3 h4 h5 p hr br radioButtons selectInput textInput actionButton uiOutput verbatimTextOutput icon tags conditionalPanel helpText div
#' @importFrom shinyjs useShinyjs disable enable
mod_metab_import_ui <- function(id) {
  ns <- shiny::NS(id)

  # Check if shinyFiles is available
  use_shiny_files <- requireNamespace("shinyFiles", quietly = TRUE)

  shiny::tagList(
    shinyjs::useShinyjs(),
    shiny::fluidRow(
      shiny::column(
        12,
        shiny::wellPanel(
          shiny::h3("Metabolomics Data Import"),
          shiny::fluidRow(
            # Left column: File import
            shiny::column(
              6,
              shiny::h4("Step 1: Select Vendor Format"),
              shiny::radioButtons(
                ns("vendor_format"),
                NULL,
                choices = c(
                  "MS-DIAL" = "msdial",
                  "Progenesis QI" = "progenesis",
                  "XCMS" = "xcms",
                  "Compound Discoverer" = "compound_discoverer",
                  "Other/Custom" = "custom"
                ),
                selected = "msdial",
                inline = TRUE
              ),
              shiny::hr(),
              shiny::h4("Step 2: Import Data Files"),
              shiny::p("Import one or more assay files (e.g., LCMS_Pos, LCMS_Neg, GCMS)")

              # Assay 1
              , shiny::wellPanel(
                style = "background-color: #f8f9fa;",
                shiny::fluidRow(
                  shiny::column(
                    4,
                    shiny::textInput(
                      ns("assay1_name"),
                      "Assay Name",
                      value = "LCMS_Pos"
                    )
                  ),
                  shiny::column(
                    8,
                    if (use_shiny_files) {
                      shiny::tagList(
                        shinyFiles::shinyFilesButton(
                          ns("assay1_file"),
                          "Select File",
                          "Choose Data File",
                          multiple = FALSE,
                          icon = shiny::icon("file")
                        ),
                        shiny::br(),
                        shiny::verbatimTextOutput(ns("assay1_path"), placeholder = TRUE)
                      )
                    } else {
                      shiny::fileInput(
                        ns("assay1_file_std"),
                        NULL,
                        accept = c(".tsv", ".tab", ".txt", ".csv", ".xlsx", ".parquet")
                      )
                    }
                  )
                )
              )

              # Assay 2 (optional)
              , shiny::wellPanel(
                style = "background-color: #f8f9fa;",
                shiny::fluidRow(
                  shiny::column(
                    4,
                    shiny::textInput(
                      ns("assay2_name"),
                      "Assay Name (Optional)",
                      value = "",
                      placeholder = "e.g., LCMS_Neg"
                    )
                  ),
                  shiny::column(
                    8,
                    if (use_shiny_files) {
                      shiny::tagList(
                        shinyFiles::shinyFilesButton(
                          ns("assay2_file"),
                          "Select File",
                          "Choose Data File",
                          multiple = FALSE,
                          icon = shiny::icon("file")
                        ),
                        shiny::br(),
                        shiny::verbatimTextOutput(ns("assay2_path"), placeholder = TRUE)
                      )
                    } else {
                      shiny::fileInput(
                        ns("assay2_file_std"),
                        NULL,
                        accept = c(".tsv", ".tab", ".txt", ".csv", ".xlsx", ".parquet")
                      )
                    }
                  )
                )
              )
            )

            # Right column: Column mapping
            , shiny::column(
              6,
              shiny::h4("Step 3: Column Mapping"),
              shiny::p("Configure column mappings for your data. Dropdowns are auto-populated based on detected format."),
              shiny::conditionalPanel(
                condition = paste0("output['", ns("file_loaded"), "']"),
                shiny::wellPanel(
                  # Format detection status (hidden for custom)
                  shiny::conditionalPanel(
                    condition = paste0("input['", ns("vendor_format"), "'] != 'custom'"),
                    shiny::uiOutput(ns("format_detection_status")),
                    shiny::hr()
                  )

                  # Custom format instructions
                  , shiny::conditionalPanel(
                    condition = paste0("input['", ns("vendor_format"), "'] == 'custom'"),
                    shiny::div(
                      class = "alert alert-info",
                      shiny::icon("edit"),
                      shiny::strong(" Custom Format"),
                      shiny::br(),
                      "Type the exact column names from your file below."
                    ),
                    shiny::hr()
                  )

                  # Metabolite ID column - Dropdown mode (non-custom)
                  , shiny::conditionalPanel(
                    condition = paste0("input['", ns("vendor_format"), "'] != 'custom'"),
                    shiny::fluidRow(
                      shiny::column(
                        6,
                        shiny::selectInput(
                          ns("metabolite_id_col"),
                          "Metabolite ID Column",
                          choices = NULL
                        )
                      ),
                      shiny::column(
                        6,
                        shiny::uiOutput(ns("metabolite_id_status"))
                      )
                    )
                  )

                  # Metabolite ID column - Text input mode (custom)
                  , shiny::conditionalPanel(
                    condition = paste0("input['", ns("vendor_format"), "'] == 'custom'"),
                    shiny::fluidRow(
                      shiny::column(
                        6,
                        shiny::textInput(
                          ns("metabolite_id_col_custom"),
                          "Metabolite ID Column",
                          value = "",
                          placeholder = "e.g., Compound_ID"
                        )
                      ),
                      shiny::column(
                        6,
                        shiny::uiOutput(ns("metabolite_id_status_custom"))
                      )
                    )
                  )

                  # Annotation column - Dropdown mode (non-custom)
                  , shiny::conditionalPanel(
                    condition = paste0("input['", ns("vendor_format"), "'] != 'custom'"),
                    shiny::fluidRow(
                      shiny::column(
                        6,
                        shiny::selectInput(
                          ns("annotation_col"),
                          "Annotation Column",
                          choices = NULL
                        )
                      ),
                      shiny::column(
                        6,
                        shiny::uiOutput(ns("annotation_status"))
                      )
                    )
                  )

                  # Annotation column - Text input mode (custom)
                  , shiny::conditionalPanel(
                    condition = paste0("input['", ns("vendor_format"), "'] == 'custom'"),
                    shiny::fluidRow(
                      shiny::column(
                        6,
                        shiny::textInput(
                          ns("annotation_col_custom"),
                          "Annotation Column (optional)",
                          value = "",
                          placeholder = "e.g., Metabolite_Name"
                        )
                      ),
                      shiny::column(
                        6,
                        shiny::uiOutput(ns("annotation_status_custom"))
                      )
                    )
                  )

                  # Sample columns pattern (custom mode only)
                  , shiny::conditionalPanel(
                    condition = paste0("input['", ns("vendor_format"), "'] == 'custom'"),
                    shiny::textInput(
                      ns("sample_cols_pattern"),
                      "Sample Column Pattern (Regex)",
                      value = "",
                      placeholder = "e.g., ^Sample_|_Rep[0-9]+$"
                    ),
                    shiny::helpText("Leave blank to auto-detect numeric columns as samples")
                  )

                  # Internal standard pattern
                  , shiny::textInput(
                    ns("is_pattern"),
                    "Internal Standard Pattern (Regex)",
                    value = "",
                    placeholder = "e.g., ^IS_|_d[0-9]+$|ISTD"
                  ),
                  shiny::helpText("Regular expression to identify internal standards"),

                  # Sanitization checkbox
                  shiny::checkboxInput(
                    ns("sanitize_names"),
                    "Sanitize Sample Names",
                    value = TRUE
                  ),
                  shiny::helpText("Clean sample IDs (e.g., '123-Sample!' -> 'x123_sample') for better compatibility with downstream analysis."),
                  shiny::hr()

                  # Sample column detection
                  , shiny::h5("Detected Sample Columns"),
                  shiny::verbatimTextOutput(ns("sample_columns_display"))

                  # Available columns helper (custom mode)
                  , shiny::conditionalPanel(
                    condition = paste0("input['", ns("vendor_format"), "'] == 'custom'"),
                    shiny::hr(),
                    shiny::h5("Available Columns in File"),
                    shiny::verbatimTextOutput(ns("available_columns_display"))
                  )
                ),
                shiny::hr()

                # Validation summary
                , shiny::h4("Step 4: Validation"),
                shiny::uiOutput(ns("validation_summary"))
              ),
              shiny::conditionalPanel(
                condition = paste0("!output['", ns("file_loaded"), "']"),
                shiny::div(
                  class = "alert alert-info",
                  shiny::icon("info-circle"),
                  " Import a data file to configure column mappings."
                )
              )
            )
          ),
          shiny::hr()

          # Process button
          , shiny::fluidRow(
            shiny::column(
              12,
              shiny::actionButton(
                ns("process_import"),
                "Process Imported Data",
                class = "btn-success",
                width = "100%",
                icon = shiny::icon("check")
              ),
              shiny::br(),
              shiny::br(),
              shiny::uiOutput(ns("import_status"))
            )
          )
        )
      )
    )
  )
}

