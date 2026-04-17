#' @rdname sessionSummaryModule
#' @export
#' @importFrom shiny NS tagList fluidPage h3 fluidRow column wellPanel h4 textInput textAreaInput actionButton icon br hr downloadButton verbatimTextOutput
mod_prot_summary_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tagList(
    shiny::fluidPage(
    shiny::h3("Session Summary & Report Generation"),
    
    shiny::fluidRow(
      # Left column: Workflow Parameters
      shiny::column(6,
        shiny::wellPanel(
          shiny::h4("Workflow Parameters"),
          shiny::textInput(ns("experiment_label"), "Experiment Label:", 
                   value = "", placeholder = "e.g., my_proteomics_analysis"),
          shiny::textAreaInput(ns("description"), "Description:", 
                       value = "Full protein analysis workflow with config parameters",
                       rows = 3, resize = "vertical"),
          shiny::br(),
          shiny::actionButton(ns("save_workflow_args"), "Save Workflow Arguments", 
                      class = "btn-primary", icon = shiny::icon("save"))
        )
      ),
      
      # Right column: File Management
      shiny::column(6,
        shiny::wellPanel(
          shiny::h4("File Management"),
          shiny::br(),
          shiny::actionButton(ns("copy_to_publication"), "Copy to Publication Directory", 
                      class = "btn-info", icon = shiny::icon("copy")),
          shiny::br(), shiny::br(),
          shiny::verbatimTextOutput(ns("copy_status"))
        )
      )
    ),
    
    shiny::fluidRow(
      # Left column: Report Generation
      shiny::column(6,
        shiny::wellPanel(
          shiny::h4("Report Generation"),
          shiny::textOutput(ns("template_status")),
          shiny::br(),
          shiny::actionButton(ns("generate_report"), "Generate Report", 
                      class = "btn-success", icon = shiny::icon("file-pdf")),
          shiny::br(), shiny::br(),
          shiny::conditionalPanel(
            condition = paste0("output['", ns("report_ready"), "']"),
            shiny::downloadButton(ns("download_report"), "Download Report", 
                          class = "btn-success")
          )
        )
      ),
      
      # Right column: GitHub Integration (Optional)
      shiny::column(6,
        shiny::wellPanel(
          shiny::h4("Version Control (Optional)"),
          shiny::checkboxInput(ns("enable_github"), "Enable GitHub Push", FALSE),
          shiny::conditionalPanel(
            condition = paste0("input['", ns("enable_github"), "']"),
            shiny::textInput(ns("github_org"), "GitHub Organization:", ""),
            shiny::textInput(ns("github_email"), "GitHub Email:", ""),
            shiny::textInput(ns("github_username"), "GitHub Username:", ""),
            shiny::textInput(ns("project_id"), "Project ID:", ""),
            shiny::br(),
            shiny::actionButton(ns("push_to_github"), "Push to GitHub", 
                        class = "btn-warning", icon = shiny::icon("github"))
          )
        )
      )
    ),
    
    shiny::fluidRow(
      shiny::column(12,
        shiny::wellPanel(
          shiny::h4("Session Summary"),
          shiny::verbatimTextOutput(ns("session_summary")),
          shiny::br(),
          shiny::actionButton(ns("export_session_state"), "Export Session State (.RDS)", 
                      class = "btn-secondary", icon = shiny::icon("download"))
        )
      )
    )
  )
  )
}

