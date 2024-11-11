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
    if (!exists("data_tbl", envir = parent.frame())) {
      stop("data_tbl not found in the current environment")
    }
    if (!exists("source_dir", envir = parent.frame())) {
      stop("source_dir not found in the current environment")
    }
    if (!exists("config_list", envir = parent.frame())) {
      stop("config_list not found in the current environment")
    }
    
    # Get required objects
    data_tbl <- get("data_tbl", envir = parent.frame())
    source_dir <- get("source_dir", envir = parent.frame())
    config_list <- get("config_list", envir = parent.frame())
    
    # Initialize data_cln with numeric conversions
    data_cln <- data_tbl |>
      mutate(Precursor.Normalised = as.numeric(Precursor.Normalised)) |>
      mutate(Precursor.Quantity = as.numeric(Precursor.Quantity))
    
    # Initialize design matrix with factor column
    design_matrix_raw <- tibble(
      Run = data_cln |> 
        distinct(Run) |>
        dplyr::select(Run) |>
        dplyr::pull(Run),
      group = NA_character_,
      factor = NA_character_,
      replicates = NA_integer_
    )
    
    # Create UI function
    create_design_matrix_ui <- function(design_matrix_raw) {
      # Sort the Run names initially
      sorted_runs <- gtools::mixedsort(design_matrix_raw$Run)
      
      fluidPage(
        titlePanel("Design Matrix Builder"),
        
        # Main layout
        fluidRow(
          # Management tools in tabs on the left
          column(4,
            wellPanel(
              tabsetPanel(
                # Sample renaming tab
                tabPanel("Rename Samples",
                  h4("Individual Rename"),
                  fluidRow(
                    column(6, selectizeInput("sample_to_rename", "Select Sample:",
                                       choices = sorted_runs)),
                    column(6, textInput("new_sample_name", "New Name:"))
                  ),
                  actionButton("rename_sample", "Rename"),
                  
                  hr(),
                  
                  h4("Bulk Rename"),
                  selectizeInput("samples_to_transform", "Select Samples:",
                            choices = sorted_runs,
                            multiple = TRUE),
                  radioButtons("transform_mode", "Transformation Mode:",
                          choices = c(
                            "Before underscore" = "before_underscore",
                            "After underscore" = "after_underscore",
                            "Character range" = "range"
                          )),
                  conditionalPanel(
                    condition = "input.transform_mode == 'range'",
                    numericInput("range_start", "Start Position:", 1),
                    numericInput("range_end", "End Position:", 1)
                  ),
                  actionButton("bulk_rename", "Apply Transformation")
                ),
                
                # Group management tab
                tabPanel("Groups",
                  h4("Add New Group"),
                  textInput("new_group", "New Group Name:"),
                  actionButton("add_group", "Add Group")
                ),
                
                # Factor management tab
                tabPanel("Factors",
                  h4("Add New Factor"),
                  textInput("new_factor", "New Factor Name:"),
                  actionButton("add_factor", "Add Factor")
                ),
                
                # Metadata assignment tab
                tabPanel("Assign Metadata",
                  h4("Assign Metadata"),
                  selectizeInput("selected_runs", "Select Runs:",
                            choices = sorted_runs,
                            multiple = TRUE),
                  selectInput("group_select", "Select Group:", choices = c("")),
                  selectInput("factor_select", "Select Factor:", choices = c("")),
                  uiOutput("replicate_inputs"),
                  actionButton("assign_metadata", "Assign")
                ),
                
                # Contrast tab
                tabPanel("Contrasts",
                  h4("Define Contrasts"),
                  selectInput("contrast_group1", "Group 1:", choices = c("")),
                  selectInput("contrast_group2", "Group 2:", choices = c("")),
                  actionButton("add_contrast", "Add Contrast")
                ),
                
                # Formula tab
                tabPanel("Formula",
                  h4("Model Formula"),
                  textInput("formula_string", "Formula:", 
                          value = config_list[["deAnalysisParameters"]][["formula_string"]])
                )
              )
            )
          ),
          
          # Design matrix and contrasts on the right
          column(8,
            # Design matrix
            wellPanel(
              h4("Current Design Matrix"),
              DTOutput("data_table")
            ),
            
            # Contrasts table
            wellPanel(
              h4("Defined Contrasts"),
              tableOutput("contrast_table")
            ),
            
            # Save button
            actionButton("save_and_close", "Save and Close", 
                        class = "btn-primary", 
                        style = "float: right;")
          )
        )
      )
    }
    
    # Create server function
    server <- function(input, output, session) {
      # Reactive values
      design_matrix <- reactiveVal(design_matrix_raw)
      data_cln_reactive <- reactiveVal(data_cln)
      
      # Initialize groups and factors
      initial_groups <- unique(design_matrix_raw$group)
      initial_groups <- initial_groups[!is.na(initial_groups) & initial_groups != ""]
      initial_factors <- unique(design_matrix_raw$factor)
      initial_factors <- initial_factors[!is.na(initial_factors) & initial_factors != ""]
      
      groups <- reactiveVal(initial_groups)
      factors <- reactiveVal(initial_factors)
      
      contrasts <- reactiveVal(data.frame(
        contrast_name = character(),
        numerator = character(),
        denominator = character(),
        stringsAsFactors = FALSE
      ))
      
      # Update dropdown choices
      observe({
        current_groups <- groups()
        current_factors <- factors()
        
        updateSelectInput(session, "group_select", 
                         choices = c("", current_groups),
                         selected = "")
        updateSelectInput(session, "factor_select", 
                         choices = c("", current_factors),
                         selected = "")
        updateSelectInput(session, "contrast_group1", 
                         choices = c("", current_groups),
                         selected = "")
        updateSelectInput(session, "contrast_group2", 
                         choices = c("", current_groups),
                         selected = "")
      })
      
      # Individual rename handler
      observeEvent(input$rename_sample, {
        req(input$sample_to_rename, input$new_sample_name)
        if(input$new_sample_name != "") {
          # Update design matrix
          current_matrix <- design_matrix()
          current_matrix$Run[current_matrix$Run == input$sample_to_rename] <- input$new_sample_name
          design_matrix(current_matrix)
          
          # Update data_cln
          current_data_cln <- data_cln_reactive()
          current_data_cln$Run[current_data_cln$Run == input$sample_to_rename] <- input$new_sample_name
          data_cln_reactive(current_data_cln)
          
          # Sort runs for UI updates
          sorted_runs <- gtools::mixedsort(current_matrix$Run)
          
          # Update UI elements
          updateSelectizeInput(session, "sample_to_rename",
                             choices = sorted_runs)
          updateSelectizeInput(session, "selected_runs",
                             choices = sorted_runs)
          updateSelectizeInput(session, "samples_to_transform",
                             choices = sorted_runs)
          
          updateTextInput(session, "new_sample_name", value = "")
        }
      })
      
      # Bulk rename handler
      observeEvent(input$bulk_rename, {
        req(input$samples_to_transform)
        current_matrix <- design_matrix()
        current_data_cln <- data_cln_reactive()
        
        # Create transformation function based on mode
        transform_fn <- function(sample_name) {
          if(input$transform_mode == "range") {
            req(input$range_start, input$range_end)
            extract_experiment(sample_name, 
                             mode = "range", 
                             start = input$range_start, 
                             end = input$range_end)
          } else if(input$transform_mode == "before_underscore") {
            extract_experiment(sample_name, mode = "start")
          } else if(input$transform_mode == "after_underscore") {
            extract_experiment(sample_name, mode = "end")
          }
        }
        
        # Apply transformations
        new_names <- sapply(input$samples_to_transform, transform_fn)
        
        # Update both matrices
        for(i in seq_along(input$samples_to_transform)) {
          current_matrix$Run[current_matrix$Run == input$samples_to_transform[i]] <- new_names[i]
          current_data_cln$Run[current_data_cln$Run == input$samples_to_transform[i]] <- new_names[i]
        }
        
        # Update reactive values
        design_matrix(current_matrix)
        data_cln_reactive(current_data_cln)
        
        # Sort runs for UI updates
        sorted_runs <- gtools::mixedsort(current_matrix$Run)
        
        # Update UI elements
        updateSelectizeInput(session, "sample_to_rename",
                           choices = sorted_runs)
        updateSelectizeInput(session, "selected_runs",
                           choices = sorted_runs)
        updateSelectizeInput(session, "samples_to_transform",
                           choices = sorted_runs)
      })
      
      # Add new group handler
      observeEvent(input$add_group, {
        req(input$new_group)
        if(input$new_group != "") {
          current_groups <- groups()
          if (!input$new_group %in% current_groups) {
            groups(c(current_groups, input$new_group))
          }
          updateTextInput(session, "new_group", value = "")
        }
      })
      
      # Add new factor handler
      observeEvent(input$add_factor, {
        req(input$new_factor)
        if(input$new_factor != "") {
          current_factors <- factors()
          if (!input$new_factor %in% current_factors) {
            factors(c(current_factors, input$new_factor))
          }
          updateTextInput(session, "new_factor", value = "")
        }
      })
      
      # Replicate inputs UI
      output$replicate_inputs <- renderUI({
        req(input$selected_runs)
        numericInput("replicate_start", 
                    paste("Starting replicate number for", length(input$selected_runs), "selected runs:"),
                    value = 1, 
                    min = 1)
      })
      
      # Assign metadata handler
      observeEvent(input$assign_metadata, {
        req(input$selected_runs, input$group_select, input$replicate_start)
        current_matrix <- design_matrix()
        
        replicate_numbers <- seq(input$replicate_start, 
                               length.out = length(input$selected_runs))
        
        current_matrix$group[current_matrix$Run %in% input$selected_runs] <- input$group_select
        if (!is.null(input$factor_select) && input$factor_select != "") {
          current_matrix$factor[current_matrix$Run %in% input$selected_runs] <- input$factor_select
        }
        current_matrix$replicates[match(input$selected_runs, current_matrix$Run)] <- replicate_numbers
        
        design_matrix(current_matrix)
        
        current_groups <- groups()
        if(!input$group_select %in% current_groups) {
          groups(c(current_groups, input$group_select))
        }
        
        if (!is.null(input$factor_select) && input$factor_select != "") {
          current_factors <- factors()
          if(!input$factor_select %in% current_factors) {
            factors(c(current_factors, input$factor_select))
          }
        }
      })
      
      # Add contrast handler
      observeEvent(input$add_contrast, {
        req(input$contrast_group1, input$contrast_group2)
        if(input$contrast_group1 != "" && 
           input$contrast_group2 != "" && 
           input$contrast_group1 != input$contrast_group2) {
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
        design_matrix_final <- design_matrix()
        data_cln_final <- data_cln_reactive()
        
        assign("design_matrix", design_matrix_final, envir = parent.frame())
        assign("data_cln", data_cln_final, envir = parent.frame())
        
        # Save formula to correct nested path in config list
        config_list[["deAnalysisParameters"]][["formula_string"]] <- input$formula_string
        assign("config_list", config_list, envir = parent.frame())
        
        contrast_data <- contrasts()
        if(nrow(contrast_data) > 0) {
          contrast_strings <- sapply(1:nrow(contrast_data), function(i) {
            paste0(
              contrast_data$contrast_name[i],
              "=group",
              contrast_data$numerator[i],
              "-group",
              contrast_data$denominator[i]
            )
          })
          
          writeLines(contrast_strings, file.path(source_dir, "contrast_strings.tab"))
          contrasts_tbl <- data.frame(
            contrasts = contrast_strings,
            stringsAsFactors = FALSE
          )
          assign("contrasts_tbl", contrasts_tbl, envir = parent.frame())
        }
        
        stopApp(list(
          design_matrix = design_matrix_final,
          contrasts_tbl = if(exists("contrasts_tbl")) contrasts_tbl else NULL,
          config_list = config_list,
          data_cln = data_cln_final
        ))
      })
    }
    
    # Create and run app
    ui <- create_design_matrix_ui(design_matrix_raw)
    
    # Run the app and capture results
    result <- runApp(
      shinyApp(ui, server),
      launch.browser = TRUE
    )
    
    # Handle results if they exist
    if (!is.null(result)) {
      assign("design_matrix", result$design_matrix, envir = parent.frame())
      if (!is.null(result$contrasts_tbl)) {
        assign("contrasts_tbl", result$contrasts_tbl, envir = parent.frame())
      }
      assign("config_list", result$config_list, envir = parent.frame())
      assign("data_cln", result$data_cln, envir = parent.frame())
    }
    invisible(result)
    
  } else if (applet_type == "future_applet") {
    stop("Future applet type not yet implemented")
  }
}