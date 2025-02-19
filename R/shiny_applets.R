RunApplet <- function(applet_type, force = FALSE) {
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
    
    # Check for existing files with timestamps
    existing_files <- list.files(source_dir, pattern = "^(design_matrix|data_cln).*\\.tab$")
    if (length(existing_files) > 0 && !force) {
      message("Found existing files in ", source_dir, ":")
      message(paste(" -", existing_files, collapse = "\n"))
      message("\nTo overwrite, rerun with force = TRUE")
      return(invisible(NULL))
    }
    
    # Get required objects
    data_tbl <- get("data_tbl", envir = parent.frame())
    source_dir <- get("source_dir", envir = parent.frame())
    config_list <- get("config_list", envir = parent.frame())
    
    # Initialize data_cln with numeric conversions
    data_cln <- data_tbl |>
      mutate(Precursor.Normalised = as.numeric(Precursor.Normalised)) |>
      mutate(Precursor.Quantity = as.numeric(Precursor.Quantity))
    
    # Initialize design matrix with factor columns
    design_matrix_raw <- tibble(
      Run = data_cln |> 
        distinct(Run) |>
        dplyr::select(Run) |>
        dplyr::pull(Run),
      group = NA_character_,
      factor1 = NA_character_,
      factor2 = NA_character_,
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
                  selectInput("factor1_select", "Select Factor 1:", choices = c("")),
                  selectInput("factor2_select", "Select Factor 2:", choices = c("")),
                  uiOutput("replicate_inputs"),
                  actionButton("assign_metadata", "Assign")
                ),
                
                # Contrast tab
                tabPanel("Contrasts",
                  h4("Define Contrasts"),
                  selectInput("contrast_group1", "Group 1:", choices = c("")),
                  selectInput("contrast_group2", "Group 2:", choices = c("")),
                  verbatimTextOutput("contrast_factors_info"),
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
      
      initial_factors <- unique(c(
        design_matrix_raw$factor1[!is.na(design_matrix_raw$factor1)],
        design_matrix_raw$factor2[!is.na(design_matrix_raw$factor2)]
      ))
      initial_factors <- initial_factors[initial_factors != ""]
      
      groups <- reactiveVal(initial_groups)
      factors <- reactiveVal(initial_factors)
      
      contrasts <- reactiveVal(data.frame(
        contrast_name = character(),
        numerator = character(),
        denominator = character(),
        factor1_num = character(),
        factor2_num = character(),
        factor1_den = character(),
        factor2_den = character(),
        stringsAsFactors = FALSE
      ))
      
      # Update dropdown choices
      observe({
        current_factors <- factors()
        
        updateSelectInput(session, "factor1_select", 
                         choices = c("", current_factors),
                         selected = "")
        updateSelectInput(session, "factor2_select", 
                         choices = c("", current_factors),
                         selected = "")
        updateSelectInput(session, "contrast_group1", 
                         choices = c("", groups()),
                         selected = "")
        updateSelectInput(session, "contrast_group2", 
                         choices = c("", groups()),
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
        req(input$selected_runs)
        req(input$factor1_select)  # Only require factor1
        
        current_matrix <- design_matrix()
        
        replicate_numbers <- seq(input$replicate_start, 
                               length.out = length(input$selected_runs))
        
        # Update the selected runs with factors
        current_matrix$factor1[current_matrix$Run %in% input$selected_runs] <- input$factor1_select
        current_matrix$factor2[current_matrix$Run %in% input$selected_runs] <- input$factor2_select
        
        # Create group names based on whether factor2 is empty or not
        group_name <- if (input$factor2_select == "") {
          input$factor1_select
        } else {
          paste(input$factor1_select, input$factor2_select, sep = "_")
        }
        
        current_matrix$group[current_matrix$Run %in% input$selected_runs] <- group_name
        
        # Update replicates
        current_matrix$replicates[match(input$selected_runs, current_matrix$Run)] <- replicate_numbers
        
        design_matrix(current_matrix)
        
        # Update groups
        unique_groups <- unique(current_matrix$group[!is.na(current_matrix$group)])
        groups(unique_groups)
        
        # Update factors
        current_factors <- factors()
        if(input$factor1_select != "" && !input$factor1_select %in% current_factors) {
          factors(c(current_factors, input$factor1_select))
        }
        if(input$factor2_select != "" && !input$factor2_select %in% current_factors) {
          factors(c(current_factors, input$factor2_select))
        }
      })
      
            # Add contrast handler
      observeEvent(input$add_contrast, {
        req(input$contrast_group1, input$contrast_group2)
        if(input$contrast_group1 != "" && 
           input$contrast_group2 != "" && 
           input$contrast_group1 != input$contrast_group2) {
          
          current_contrasts <- contrasts()
          current_matrix <- design_matrix()
          
          # Get the factors for each group
          group1_row <- which(current_matrix$group == input$contrast_group1)[1]
          group2_row <- which(current_matrix$group == input$contrast_group2)[1]
          
          contrast_name <- paste0(
            gsub(" ", "", input$contrast_group1), 
            ".minus.", 
            gsub(" ", "", input$contrast_group2)
          )
          
          new_contrast <- data.frame(
            contrast_name = contrast_name,
            numerator = input$contrast_group1,
            denominator = input$contrast_group2,
            factor1_num = current_matrix$factor1[group1_row],
            factor2_num = current_matrix$factor2[group1_row],
            factor1_den = current_matrix$factor1[group2_row],
            factor2_den = current_matrix$factor2[group2_row],
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
            Contrast = sapply(1:nrow(contrast_data), function(i) {
              # Create numerator group name
              num_group <- if (contrast_data$factor2_num[i] == "") {
                contrast_data$factor1_num[i]
              } else {
                paste(contrast_data$factor1_num[i], contrast_data$factor2_num[i], sep = "_")
              }
              
              # Create denominator group name
              den_group <- if (contrast_data$factor2_den[i] == "") {
                contrast_data$factor1_den[i]
              } else {
                paste(contrast_data$factor1_den[i], contrast_data$factor2_den[i], sep = "_")
              }
              
              # Create full contrast string
              paste0(
                num_group,
                ".minus.",
                den_group,
                " = group",
                num_group,
                "-group",
                den_group
              )
            })
          )
        }
      })
      
            # Save and close handler
      observeEvent(input$save_and_close, {
        design_matrix_final <- design_matrix()
        data_cln_final <- data_cln_reactive()
        
        # Filter design_matrix_final to only include rows with assigned metadata
        design_matrix_final <- design_matrix_final |>
            filter(!is.na(group))
        
        # Get assigned runs directly from filtered design matrix        
        assigned_runs <- design_matrix_final |>
            pull(Run)
        
        data_cln_final <- data_cln_final |>
            filter(Run %in% assigned_runs)
        
        assign("design_matrix", design_matrix_final, envir = parent.frame())
        assign("data_cln", data_cln_final, envir = parent.frame())
        
        # Save formula to correct nested path in config list
        config_list[["deAnalysisParameters"]][["formula_string"]] <- input$formula_string
        assign("config_list", config_list, envir = parent.frame())
        
        contrast_data <- contrasts()
        if(nrow(contrast_data) > 0) {
          # Modified contrast string generation to match required format
          contrast_strings <- sapply(1:nrow(contrast_data), function(i) {
            # Create numerator group name
            num_group <- if (contrast_data$factor2_num[i] == "") {
              contrast_data$factor1_num[i]
            } else {
              paste(contrast_data$factor1_num[i], contrast_data$factor2_num[i], sep = "_")
            }
            
            # Create denominator group name
            den_group <- if (contrast_data$factor2_den[i] == "") {
              contrast_data$factor1_den[i]
            } else {
              paste(contrast_data$factor1_den[i], contrast_data$factor2_den[i], sep = "_")
            }
            
            paste0(
              num_group,
              ".minus.",
              den_group,
              "=group",
              num_group,
              "-group",
              den_group
            )
          })
          
          writeLines(contrast_strings, file.path(source_dir, "contrast_strings.tab"))
          contrasts_tbl <- data.frame(
            contrasts = contrast_strings,
            stringsAsFactors = FALSE
          )
          assign("contrasts_tbl", contrasts_tbl, envir = parent.frame())
        }
        
        # Create base filename
        base_file <- file.path(source_dir, "design_matrix.tab")
        
        # If file exists, append timestamp
        if (file.exists(base_file)) {
          timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
          output_file <- file.path(source_dir, paste0("design_matrix_", timestamp, ".tab"))
          data_cln_output_file <- file.path(source_dir, paste0("data_cln_", timestamp, ".tab"))
        } else {
          output_file <- base_file
          data_cln_output_file <- file.path(source_dir, "data_cln.tab")
        }
        
        write.table(design_matrix_final, 
                    file = output_file,
                    sep = "\t",
                    row.names = FALSE,
                    quote = FALSE)
                    
        write.table(data_cln_final, 
                    file = data_cln_output_file,
                    sep = "\t",
                    row.names = FALSE,
                    quote = FALSE)
        
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