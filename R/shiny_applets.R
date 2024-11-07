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
  sidebar_content <- tagList(
    # Group management section
    h4("Add New Group"),
    textInput("new_group", "Enter new group name:"),
    actionButton("add_group", "Add Group"),
    
    hr(),
    
    # Metadata assignment section
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
    
    uiOutput("replicate_inputs"),
    
    actionButton(
      "assign_metadata", 
      "Assign Metadata",
      class = "btn-primary",
      style = "margin-top: 10px; margin-bottom: 10px;"
    ),
    
    hr(),
    
    # Contrast section
    h4("Add Contrasts"),
    selectInput(
      "contrast_group1",
      "Numerator Group:",
      choices = c(""),
      multiple = FALSE
    ),
    selectInput(
      "contrast_group2", 
      "Denominator Group:",
      choices = c(""),
      multiple = FALSE
    ),
    actionButton(
      "add_contrast",
      "Add Contrast",
      class = "btn-info",
      style = "margin-top: 10px; margin-bottom: 10px;"
    ),
    
    hr(),
    
    # Save button
    actionButton(
      "save_and_close",
      "Save and Close",
      class = "btn-success"
    )
  )
  
  main_content <- tagList(
    DTOutput("data_table"),
    h4("Defined Contrasts"),
    tableOutput("contrast_table")
  )
  
  create_base_ui(
    title = "Metadata Assignment Tool",
    sidebar_content = sidebar_content,
    main_content = main_content
  )
}

#' Design Matrix Server Logic
#' @keywords internal
design_matrix_server <- function(input, output, session) {
  # Reactive values for design matrix
  design_matrix <- reactiveVal(design_matrix_raw)
  
  # Reactive values for groups
  groups <- reactiveVal(character())
  
  # Reactive values for contrasts
  contrasts <- reactiveVal(data.frame(
    contrast_name = character(),
    numerator = character(),
    denominator = character(),
    stringsAsFactors = FALSE
  ))
  
  # Update group choices for contrast selection
  observe({
    current_groups <- unique(design_matrix()$Group)
    current_groups <- current_groups[!is.na(current_groups)]
    updateSelectInput(session, "contrast_group1", choices = current_groups)
    updateSelectInput(session, "contrast_group2", choices = current_groups)
    updateSelectInput(session, "group_select", choices = current_groups)
  })
  
  # Add new group handler
  observeEvent(input$add_group, {
    req(input$new_group)
    current_groups <- groups()
    if (!input$new_group %in% current_groups) {
      groups(c(current_groups, input$new_group))
    }
  })
  
  # Assign metadata handler
  observeEvent(input$assign_metadata, {
    req(input$selected_runs, input$group_select)
    current_matrix <- design_matrix()
    current_matrix$Group[current_matrix$Run %in% input$selected_runs] <- input$group_select
    design_matrix(current_matrix)
  })
  
  # Add contrast button handler
  observeEvent(input$add_contrast, {
    req(input$contrast_group1, input$contrast_group2)
    if(input$contrast_group1 != input$contrast_group2) {
      current_contrasts <- contrasts()
      contrast_name <- paste0(
        gsub(" ", "", input$contrast_group1), 
        ".minus.", 
        gsub(" ", "", input$contrast_group2)
      )
      new_contrast <- data.frame(
        contrast_name = contrast_name,
        numerator = input$contrast_group1,
        denominator = input$contrast_group2,
        stringsAsFactors = FALSE
      )
      contrasts(rbind(current_contrasts, new_contrast))
    }
  })
  
  # Display data table
  output$data_table <- renderDT({
    design_matrix()
  }, selection = 'none', options = list(pageLength = 25))
  
  # Display contrast table
  output$contrast_table <- renderTable({
    contrast_data <- contrasts()
    if(nrow(contrast_data) > 0) {
      data.frame(
        Contrast = paste0(
          contrast_data$contrast_name,
          "=",
          contrast_data$numerator,
          "-",
          contrast_data$denominator
        )
      )
    }
  })
  
  # Save and close handler
  observeEvent(input$save_and_close, {
    # Save the design matrix
    design_matrix_final <- design_matrix()
    assign("design_matrix", design_matrix_final, envir = parent.frame())
    
    # Create and save contrast strings file
    contrast_data <- contrasts()
    if(nrow(contrast_data) > 0) {
      # Create contrast strings in required format
      contrast_strings <- c(
        "contrasts",
        sapply(1:nrow(contrast_data), function(i) {
          paste0(
            contrast_data$contrast_name[i],
            "=",
            contrast_data$numerator[i],
            "-",
            contrast_data$denominator[i]
          )
        })
      )
      
      # Write contrast strings to file
      writeLines(
        contrast_strings,
        file.path(source_dir, "contrast_strings.tab")
      )
      
      # Create contrasts_tbl
      contrasts_tbl <- data.frame(
        numerator = contrast_data$numerator,
        denominator = contrast_data$denominator,
        stringsAsFactors = FALSE
      )
      assign("contrasts_tbl", contrasts_tbl, envir = parent.frame())
    }
    
    # Close the app
    stopApp(list(
      design_matrix = design_matrix_final,
      contrasts_tbl = if(exists("contrasts_tbl")) contrasts_tbl else NULL
    ))
  })
}

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
    
    # Run the app and capture results
    result <- runApp(shinyApp(ui, server))
    
    # Handle results if they exist
    if (!is.null(result)) {
      assign("design_matrix", result$design_matrix, envir = parent.frame())
      if (!is.null(result$contrasts_tbl)) {
        assign("contrasts_tbl", result$contrasts_tbl, envir = parent.frame())
      }
    }
    invisible(result)
    
  } else if (applet_type == "future_applet") {
    stop("Future applet type not yet implemented")
  }
}