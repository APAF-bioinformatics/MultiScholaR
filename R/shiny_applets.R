#' Create Base UI Structure
#'
#' @param title Character string for the title panel
#' @param sidebar_content Shiny UI elements for the sidebar
#' @param main_content Shiny UI elements for the main panel
#' @return A Shiny fluidPage UI
#' @keywords internal
create_base_ui <- function(title, sidebar_content, main_content) {
  fluidPage(
    titlePanel(title),
    sidebarLayout(
      sidebarPanel(sidebar_content),
      mainPanel(main_content)
    )
  )
}

#########################################################################################

#' Design Matrix Applet Components
#' @keywords internal
create_design_matrix_ui <- function(design_matrix_raw) {
  sidebar_content <- tagList(
    h4("Add New Group"),
    textInput("new_group", "Enter new group name:"),
    actionButton("add_group", "Add Group"),
    
    hr(),
    
    h4("Assign Metadata"),
    selectizeInput("selected_runs", "Select Runs:",
                  choices = gtools::mixedsort(design_matrix_raw$Run),
                  multiple = TRUE),
    
    selectInput("group_select", "Select Group:",
               choices = c("")),
    
    uiOutput("replicate_inputs"),
    
    actionButton("assign_metadata", "Assign Metadata"),
    
    hr(),
    
    downloadButton("download_data", "Download, Save and Close")
  )
  
  main_content <- DTOutput("data_table")
  
  create_base_ui(
    title = "Metadata Assignment Tool",
    sidebar_content = sidebar_content,
    main_content = main_content
  )
}

#' Design Matrix Server Logic
#' @keywords internal
design_matrix_server <- function(input, output, session) {
  # [Previous server code remains the same]
  # ... [Insert the entire server_designMatrix function content here]
}


#########################################################################################
#' Create Design Matrix UI
#'
#' Creates the user interface for the design matrix applet
#' @param design_matrix_raw A data frame containing the raw design matrix
#' @return A Shiny UI object
#' @keywords internal
#' 
#' @importFrom shiny h4 textInput actionButton hr selectizeInput selectInput uiOutput 
#' @importFrom shiny downloadButton tagList
#' @importFrom DT DTOutput
#' @importFrom gtools mixedsort
create_design_matrix_ui <- function(design_matrix_raw) {
  # Sidebar content
  sidebar_content <- tagList(
    # Group management 
    h4("Add New Group"),
    textInput("new_group", "Enter new group name:"),
    actionButton("add_group", "Add Group"),
    
    hr(),
    
    # Metadata assignment
    h4("Assign Metadata"),
    selectizeInput(
      "selected_runs", 
      "Select Runs:",
      choices = gtools::mixedsort(design_matrix_raw$Run),
      multiple = TRUE,
      options = list(
        placeholder = 'Click to select runs',
        plugins = list('remove_button'),
        closeAfterSelect = FALSE
      )
    ),
    
    selectInput(
      "group_select", 
      "Select Group:",
      choices = c(""),
      selected = NULL,
      multiple = FALSE
    ),
    
    # Dynamic input rendering
    uiOutput("replicate_inputs"),
    
    actionButton(
      "assign_metadata", 
      "Assign Metadata",
      class = "btn-primary",
      style = "margin-top: 10px; margin-bottom: 10px;"
    ),
    
    hr(),
    
    # Download section
    downloadButton(
      "download_data", 
      "Download, Save and Close",
      class = "btn-success"
    )
  )
  
  # Main panel content
  main_content <- tagList(
    DTOutput("data_table")
  )
  
  # Create UI
  create_base_ui(
    title = "Metadata Assignment Tool",
    sidebar_content = sidebar_content,
    main_content = main_content
  )
}

#########################################################################################

#' Run Various Shiny Applets
#'
#' A function to run different types of Shiny applets with consistent UI structure.
#'
#' @param applet_type The type of applet to run. Current options: "designMatrix"
#' @return Returns the result of the selected applet invisibly
#' @export
#'
#' @examples
#' \dontrun{
#' RunApplet("designMatrix")
#' }
RunApplet <- function(applet_type) {
  # Validation
  valid_types <- c("designMatrix")
  if (!applet_type %in% valid_types) {
    stop("Invalid applet type. Valid options are: ", 
         paste(valid_types, collapse = ", "))
  }
  
  # Handle different applet types
  if (applet_type == "designMatrix") {
    # Check for required objects
    if (!exists("design_matrix_raw", envir = parent.frame())) {
      stop("design_matrix_raw not found in the current environment")
    }
    if (!exists("source_dir", envir = parent.frame())) {
      stop("source_dir not found in the current environment")
    }
    
    # Get required objects
    design_matrix_raw <- get("design_matrix_raw", envir = parent.frame())
    source_dir <- get("source_dir", envir = parent.frame())
    
    # Create and run app
    ui <- create_design_matrix_ui(design_matrix_raw)
    server <- design_matrix_server
    
  } else if (applet_type == "future_applet") {
			# Insert future applet type here
    stop("Future applet type not yet implemented")
  }
  
  # Run the app
  result <- shinyApp(ui, server)
  invisible(result)
}