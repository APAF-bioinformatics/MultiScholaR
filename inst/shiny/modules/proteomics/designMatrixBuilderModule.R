#' @title designMatrixBuilderModule
#'
#' @description A Shiny module to build a design matrix for proteomics experiments.
#' This module is a refactored version of the standalone design matrix applet,
#' designed to be embedded within the main MultiScholaR Shiny application.
#'
#' @param id The module ID.
#' @param data_tbl A reactive expression that returns the data table
#'   (e.g., from the setup & import tab).
#' @param config_list A reactive expression that returns the main configuration list.
#'
#' @return A reactive expression that returns a list containing the results
#'   when the "Save" button is clicked. The list includes:
#'   - `design_matrix`: The final, filtered design matrix.
#'   - `data_cln`: The data table, filtered to include only assigned runs.
#'   - `contrasts_tbl`: A data frame of the defined contrasts.
#'   - `config_list`: The updated configuration list with the new formula.
#'
#' @import shiny
#' @import shinydashboard
#' @import DT
#' @import gtools
#' @import dplyr
#' @importFrom tibble tibble
#' @export
NULL

#' Design Matrix Builder UI Module
#'
#' @param id Module ID.
#' @return A tagList containing the UI for the design matrix builder.
#' @rdname designMatrixBuilderModule
designMatrixBuilderUI <- function(id) {
    ns <- NS(id)
    # The UI is extensive, so it is defined here. It is based on the
    # createDesignMatrixUi function from the original shiny_applets.R file.
    shiny::fluidPage(
        titlePanel("Design Matrix Builder"),

        # Main layout
        shiny::fluidRow(
            # Management tools and Info panel on the left
            shiny::column(
                4,
                shiny::wellPanel(
                    shiny::tabsetPanel(
                        id = ns("main_tabs"),
                        # Sample renaming tab
                        shiny::tabPanel(
                            "Rename Samples",
                            shiny::h4("Individual Rename"),
                            shiny::fluidRow(
                                shiny::column(6, shiny::selectizeInput(ns("sample_to_rename"), "Select Sample:",
                                    choices = NULL # Will be populated by server
                                )),
                                shiny::column(6, shiny::textInput(ns("new_sample_name"), "New Name:"))
                            ),
                            shiny::actionButton(ns("rename_sample"), "Rename"),
                            shiny::hr(),
                            shiny::h4("Bulk Rename"),
                            shiny::selectizeInput(ns("samples_to_transform"), "Select Samples:",
                                choices = NULL, # Will be populated by server
                                multiple = TRUE
                            ),
                            shiny::radioButtons(ns("transform_mode"), "Transformation Mode:",
                                choices = c(
                                    "Keep Text BEFORE First Underscore" = "before_underscore",
                                    "Keep Text AFTER Last Underscore" = "after_underscore",
                                    "Extract by Position" = "range"
                                )
                            ),
                            shiny::helpText(HTML(
                                "Select how to shorten the selected sample names using underscores as delimiters.",
                                "<br><b>Keep Text BEFORE First Underscore:</b> Retains only the part before the very first '_'. Example: 'SampleA_T1_Rep1' &rarr; 'SampleA'.",
                                "<br><b>Keep Text AFTER Last Underscore:</b> Retains only the part after the very last '_'. Example: 'SampleA_T1_Rep1' &rarr; 'Rep1'.",
                                "<br><b>Extract by Underscore Position:</b> Keeps text segments based on underscore indices (1-based).",
                                "<ul>",
                                "<li>Start = 1: Text before the 1st underscore.</li>",
                                "<li>End = 1: Text after the last underscore.</li>",
                                "<li>Start = N, End = N: Text between the (N-1)th and Nth underscore.</li>",
                                "<li>Start = N, End = M (N < M): Text between the (N-1)th and Mth underscore.</li>",
                                "</ul>",
                                "Example: 'SampleA_T1_Rep1' with Start=2, End=2 &rarr; 'T1'."
                            )),
                            shiny::conditionalPanel(
                                condition = paste0("input['", ns("transform_mode"), "'] == 'range'"),
                                shiny::numericInput(ns("range_start"), "Start Underscore Index:", 1, min = 1),
                                shiny::numericInput(ns("range_end"), "End Underscore Index:", 1, min = 1),
                                shiny::h5("Preview (First Selected Sample):"),
                                shiny::verbatimTextOutput(ns("range_preview"))
                            ),
                            shiny::actionButton(ns("bulk_rename"), "Apply Transformation")
                        ),

                        # Factor management tab
                        shiny::tabPanel(
                            "Factors",
                            shiny::h4("Add New Factor"),
                            shiny::textInput(ns("new_factor"), "New Factor Name:"),
                            shiny::actionButton(ns("add_factor"), "Add Factor")
                        ),

                        # Metadata assignment tab
                        shiny::tabPanel(
                            "Assign Metadata",
                            shiny::h4("Assign Metadata"),
                            shiny::selectizeInput(ns("selected_runs"), "Select Runs:",
                                choices = NULL, # Will be populated by server
                                multiple = TRUE
                            ),
                            shiny::selectInput(ns("factor1_select"), "Select Factor 1:", choices = c("")),
                            shiny::selectInput(ns("factor2_select"), "Select Factor 2:", choices = c("")),
                            shiny::uiOutput(ns("replicate_inputs")),
                            shiny::actionButton(ns("assign_metadata"), "Assign")
                        ),

                        # Technical replicates tab
                        shiny::tabPanel(
                            "Technical Replicates",
                            shiny::h4("Assign Technical Replicates"),
                            shiny::helpText(HTML(
                                "Select samples that are technical replicates of each other. ",
                                "These samples will be grouped under the same biological replicate number ",
                                "and assigned distinct technical replicate numbers.<br><br>",
                                "<b>Example:</b> If you select samples 3A, 3B, 3C with replicate numbers 3, 4, 5, ",
                                "they will all be assigned replicate number 3, with tech_reps 1, 2, 3 respectively."
                            )),
                            shiny::selectizeInput(ns("tech_rep_samples"), "Select Technical Replicates:",
                                choices = NULL, # Will be populated by server
                                multiple = TRUE
                            ),
                            shiny::radioButtons(ns("tech_rep_assignment_mode"), "Assignment Mode:",
                                choices = c(
                                    "Use lowest replicate number" = "lowest",
                                    "Use first selected sample's replicate" = "first",
                                    "Specify replicate number" = "manual"
                                )
                            ),
                            shiny::conditionalPanel(
                                condition = paste0("input['", ns("tech_rep_assignment_mode"), "'] == 'manual'"),
                                shiny::numericInput(ns("manual_replicate_number"), "Replicate Number:", 
                                    value = 1, min = 1)
                            ),
                            shiny::actionButton(ns("assign_tech_reps"), "Assign Technical Replicates", 
                                class = "btn-primary"),
                            shiny::hr(),
                            shiny::h5("Current Technical Replicate Assignments"),
                            shiny::verbatimTextOutput(ns("tech_rep_summary"))
                        ),

                        # Contrast tab
                        shiny::tabPanel(
                            "Contrasts",
                            shiny::h4("Define Contrasts"),
                            shiny::selectInput(ns("contrast_group1"), "Group 1:", choices = c("")),
                            shiny::selectInput(ns("contrast_group2"), "Group 2:", choices = c("")),
                            shiny::verbatimTextOutput(ns("contrast_factors_info")),
                            shiny::actionButton(ns("add_contrast"), "Add Contrast")
                        ),

                        # Formula tab
                        shiny::tabPanel(
                            "Formula",
                            shiny::h4("Model Formula"),
                            shiny::textInput(ns("formula_string"), "Formula:",
                                value = "" # Will be populated by server
                            )
                        ),

                        # Settings tab for global operations
                        shiny::tabPanel(
                            "Settings",
                            shiny::h4("Session Management"),
                            shiny::p("Reset options allow you to revert changes made during this session."),
                            shiny::div(
                                style = "margin-top: 15px;",
                                shiny::selectInput(ns("reset_scope"), "Reset Scope:",
                                    choices = c(
                                        "All Changes" = "all",
                                        "Sample Names Only" = "sample_names",
                                        "Factors Only" = "factors",
                                        "Contrasts Only" = "contrasts",
                                        "Formula Only" = "formula"
                                    ),
                                    selected = "all"
                                )
                            ),
                            shiny::div(
                                style = "margin-top: 15px;",
                                shiny::actionButton(ns("reset_changes"), "Reset to Initial State",
                                    class = "btn-warning",
                                    icon = icon("rotate-left"),
                                    width = "100%"
                                )
                            ),
                            shiny::hr(),
                            shiny::h4("Information"),
                            shiny::verbatimTextOutput(ns("session_info"))
                        )
                    )
                ),

                # Information & Assignments Panel
                shiny::wellPanel(
                    shiny::h4("Information & Assignments"),
                    shiny::helpText(HTML(
                        "<b>Factors:</b> These represent the main experimental variables or conditions ",
                        "that distinguish your samples (e.g., Treatment, Timepoint, Genotype, Condition). ",
                        "They are used to group runs for analysis. Define potential factor names in the 'Factors' tab, ",
                        "then assign the appropriate factor levels to your runs using the 'Assign Metadata' tab. ",
                        "A 'group' is automatically created based on the combination of assigned factors (e.g., 'TreatmentA_Timepoint1')."
                    )),
                    shiny::helpText(HTML(
                        "<b>Technical Replicates:</b> These are repeated measurements of the same biological sample ",
                        "(e.g., multiple injections or runs of the same sample). Use the 'Technical Replicates' tab to group ",
                        "samples that are technical replicates. They will share the same biological replicate number but have ",
                        "distinct technical replicate numbers. This information is used downstream for missing value imputation."
                    )),
                    shiny::helpText(HTML(
                        "<b>Contrasts:</b> These define the specific statistical comparisons you want to make ",
                        "between the groups you've defined (e.g., 'TreatmentA vs Control', 'Timepoint2 vs Timepoint1'). ",
                        "Contrasts are essential for differential expression/abundance analysis, allowing the statistical model ",
                        "to test specific hypotheses about differences between your experimental groups. Define them in the 'Contrasts' tab ",
                        "after assigning metadata and creating groups."
                    )),
                    shiny::hr(),
                    shiny::h5("Available Factors (Defined in 'Factors' Tab)"),
                    shiny::uiOutput(ns("available_factors_display")),
                    shiny::hr(),
                    shiny::h5("Defined Contrasts"),
                    shiny::uiOutput(ns("defined_contrasts_display"))
                )
            ),

            # Design matrix display on the right
            shiny::column(
                8,
                shiny::wellPanel(
                    shiny::h4("Current Design Matrix (All Runs) - Proteomics"),
                    DT::DTOutput(ns("data_table"))
                ),
                # The save button is now the primary action for the module to return data.
                # The parent module will handle closing/hiding this UI.
                shiny::actionButton(ns("save_results"), "Save Design",
                    class = "btn-primary",
                    icon = shiny::icon("save"),
                    style = "float: right;"
                )
            )
        )
    )
}


#' Design Matrix Builder Server Module
#' @rdname designMatrixBuilderModule
designMatrixBuilderServer <- function(id, data_tbl, config_list) {
    moduleServer(id, function(input, output, session) {

        # Reactive value to store the final results
        result_rv <- reactiveVal(NULL)

        # == Initial State Setup =================================================
        # This section sets up the initial state from the provided data.
        # It runs only once when the module is initialized.

        # Create a non-reactive list for the initial state to enable resets.
        initial_state <- reactive({
            req(data_tbl())
            req(config_list())

            df <- data_tbl()
            conf <- config_list()

            # Initialize data_cln with numeric conversions (Proteomics specific)
            data_cln <- df |>
                dplyr::mutate(Precursor.Normalised = as.numeric(Precursor.Normalised)) |>
                dplyr::mutate(Precursor.Quantity = as.numeric(Precursor.Quantity))

            # Initialize design matrix with factor columns
            design_matrix_raw <- tibble::tibble(
                Run = data_cln |>
                    dplyr::distinct(Run) |>
                    dplyr::select(Run) |>
                    dplyr::pull(Run) |>
                    gtools::mixedsort(),
                group = NA_character_,
                factor1 = NA_character_,
                factor2 = NA_character_,
                replicates = NA_integer_,
                tech_reps = NA_integer_
            )

            list(
                design_matrix = design_matrix_raw,
                data_cln = data_cln,
                groups = unique(design_matrix_raw$group[!is.na(design_matrix_raw$group) & design_matrix_raw$group != ""]),
                factors = character(0), # No factors defined initially
                formula = conf[["deAnalysisParameters"]][["formula_string"]],
                contrasts = data.frame(
                    contrast_name = character(),
                    numerator = character(),
                    denominator = character(),
                    stringsAsFactors = FALSE
                )
            )
        })

        # Reactive values for the current, mutable state of the builder
        design_matrix <- reactiveVal()
        data_cln_reactive <- reactiveVal()
        groups <- reactiveVal()
        factors <- reactiveVal()
        contrasts <- reactiveVal()

        # Observer to initialize/reset the reactive values from the initial state
        observe({
            state <- initial_state()
            req(state)
            design_matrix(state$design_matrix)
            data_cln_reactive(state$data_cln)
            groups(state$groups)
            factors(state$factors)
            contrasts(state$contrasts)
            
            # Update UI elements that depend on the initial state
            sorted_runs <- state$design_matrix$Run
            updateSelectizeInput(session, "sample_to_rename", choices = sorted_runs, selected = "")
            updateSelectizeInput(session, "selected_runs", choices = sorted_runs, selected = "")
            updateSelectizeInput(session, "samples_to_transform", choices = sorted_runs, selected = "")
            updateSelectizeInput(session, "tech_rep_samples", choices = sorted_runs, selected = "")
            updateTextInput(session, "formula_string", value = state$formula)
        })

        # Create a proxy for the main data table to allow updates without state loss
        proxy_data_table <- DT::dataTableProxy("data_table")

        # == UI Rendering and Updates =============================================

        # Update sample selection inputs when names change
        observeEvent(design_matrix(), {
            # This observer updates all sample dropdowns if the run names in the
            # design matrix change (e.g., through renaming).
            req(design_matrix())
            sorted_runs <- gtools::mixedsort(design_matrix()$Run)
            
            # Preserve existing selections if possible
            isolate({
                selected_rename <- input$sample_to_rename
                selected_meta <- input$selected_runs
                selected_transform <- input$samples_to_transform
                selected_tech <- input$tech_rep_samples
            })

            updateSelectizeInput(session, "sample_to_rename", choices = sorted_runs, selected = selected_rename)
            updateSelectizeInput(session, "selected_runs", choices = sorted_runs, selected = selected_meta)
            updateSelectizeInput(session, "samples_to_transform", choices = sorted_runs, selected = selected_transform)
            updateSelectizeInput(session, "tech_rep_samples", choices = sorted_runs, selected = selected_tech)
        }, ignoreNULL = TRUE, ignoreInit = TRUE)

        # Update factor and group dropdowns when they change
        observe({
            # This observer keeps the factor and group dropdowns in sync with the
            # reactive values `factors()` and `groups()`.
            
            # Preserve selections
            isolate({
                f1 <- input$factor1_select
                f2 <- input$factor2_select
                g1 <- input$contrast_group1
                g2 <- input$contrast_group2
            })
            
            current_factors <- factors()
            current_groups <- groups()

            updateSelectInput(session, "factor1_select", choices = c("", current_factors), selected = f1)
            updateSelectInput(session, "factor2_select", choices = c("", current_factors), selected = f2)
            updateSelectInput(session, "contrast_group1", choices = c("", current_groups), selected = g1)
            updateSelectInput(session, "contrast_group2", choices = c("", current_groups), selected = g2)
        })


        # Render the main data table
        output$data_table <- DT::renderDT({
                req(design_matrix())
                design_matrix()
            },
            selection = "none",
            options = list(pageLength = 10, scrollX = TRUE, server = FALSE)
        )

        # Use a proxy to update the table data without redrawing the whole thing
        observeEvent(design_matrix(), {
            req(proxy_data_table)
            DT::replaceData(proxy_data_table, design_matrix(), resetPaging = FALSE)
        })

        # Render UI for Available Factors Display
        output$available_factors_display <- renderUI({
            current_factors <- factors()
            if (length(current_factors) == 0) {
                return(shiny::p("No factors defined yet (use the 'Factors' tab)."))
            }
            shiny::p(paste(current_factors, collapse = ", "))
        })

        # Render UI for Defined Contrasts Display
        output$defined_contrasts_display <- renderUI({
            contrast_data <- contrasts()
            if (is.null(contrast_data) || nrow(contrast_data) == 0) {
                return(shiny::p("No contrasts defined yet."))
            }
            
            # Check if formula uses ~ 0 + group pattern
            formula_uses_group_prefix <- grepl("~ *0 *\\+ *group", input$formula_string)
            
            # Generate the actual contrast format that will be used
            contrast_info <- lapply(1:nrow(contrast_data), function(i) {
                # Create friendly name (replace dots with underscores for consistency)
                friendly_name <- gsub("\\.", "_", contrast_data$contrast_name[i])
                friendly_name <- gsub("\\.vs\\.", "_minus_", friendly_name)
                
                # Create actual contrast string based on formula
                if (formula_uses_group_prefix) {
                    contrast_string <- paste0(
                        "group", contrast_data$numerator[i],
                        "-group", contrast_data$denominator[i]
                    )
                } else {
                    contrast_string <- paste0(
                        contrast_data$numerator[i],
                        "-", contrast_data$denominator[i]
                    )
                }
                
                # Return the full format
                paste0(friendly_name, "=", contrast_string)
            })
            
            # Display the full format
            shiny::tagList(
                lapply(contrast_info, function(full_format) {
                    shiny::p(shiny::tags$code(full_format))
                })
            )
        })
        
        # Render preview for range extraction transformation
        output$range_preview <- renderText({
            req(input$samples_to_transform)
            req(input$range_start, input$range_end)
            
            # Get the first selected sample for preview
            first_sample <- input$samples_to_transform[1]
            
            # Apply the transformation
            tryCatch({
                preview_result <- extract_experiment(
                    first_sample, 
                    mode = "range", 
                    start = input$range_start, 
                    end = input$range_end
                )
                paste0("\"", first_sample, "\" → \"", preview_result, "\"")
            }, error = function(e) {
                paste("Error:", e$message)
            })
        })
        
        # Render technical replicate summary
        output$tech_rep_summary <- renderText({
            dm <- design_matrix()
            req(dm)
            
            # Find samples with technical replicate assignments
            tech_rep_groups <- dm |>
                dplyr::filter(!is.na(tech_reps)) |>
                dplyr::group_by(group, replicates) |>
                dplyr::summarise(
                    samples = paste(Run, collapse = ", "),
                    tech_rep_numbers = paste(tech_reps, collapse = ", "),
                    .groups = "drop"
                )
            
            if (nrow(tech_rep_groups) == 0) {
                return("No technical replicates assigned yet.")
            }
            
            # Format the output using vectorized approach
            output_text <- tech_rep_groups |>
                dplyr::mutate(
                    formatted = sprintf("Group: %s, Biological Replicate: %s\n  Samples: %s\n  Technical Replicates: %s",
                                      group, replicates, samples, tech_rep_numbers)
                ) |>
                dplyr::pull(formatted)
            
            paste(output_text, collapse = "\n\n")
        })
        
        # Render replicate number input UI
        output$replicate_inputs <- renderUI({
            ns <- session$ns
            req(input$selected_runs)
            shiny::numericInput(ns("replicate_start"),
                paste("Starting replicate number for", length(input$selected_runs), "selected runs:"),
                value = 1,
                min = 1
            )
        })

        # Render contrast factors info
        output$contrast_factors_info <- renderText({
            formula_uses_group_prefix <- grepl("~ *0 *\\+ *group", input$formula_string)
            if (formula_uses_group_prefix) {
                "Note: Contrasts will use 'group' prefix (e.g., groupGA_Control-groupGA_Elevated)\nbased on current formula: ~ 0 + group"
            } else {
                "Note: Contrasts will use group names as-is (e.g., GA_Control-GA_Elevated)"
            }
        })

        # == Event Handlers for UI Actions ========================================

        # Handler for individual sample renaming
        observeEvent(input$rename_sample, {
            req(input$sample_to_rename, input$new_sample_name)
            
            # Get current state
            current_matrix <- design_matrix()
            current_data_cln <- data_cln_reactive()

            # Update names in both tables
            current_matrix$Run[current_matrix$Run == input$sample_to_rename] <- input$new_sample_name
            current_data_cln$Run[current_data_cln$Run == input$sample_to_rename] <- input$new_sample_name

            # Update reactive values
            design_matrix(current_matrix)
            data_cln_reactive(current_data_cln)

            # Clear the input
            updateTextInput(session, "new_sample_name", value = "")
        })

        # Handler for bulk renaming
        observeEvent(input$bulk_rename, {
            req(input$samples_to_transform)
            
            current_matrix <- design_matrix()
            current_data_cln <- data_cln_reactive()

            # Transformation function based on selected mode
            transform_fn <- function(sample_name) {
                if (input$transform_mode == "range") {
                    req(input$range_start, input$range_end)
                    extract_experiment(sample_name, mode = "range", start = input$range_start, end = input$range_end)
                } else if (input$transform_mode == "before_underscore") {
                    extract_experiment(sample_name, mode = "start")
                } else if (input$transform_mode == "after_underscore") {
                    extract_experiment(sample_name, mode = "end")
                }
            }

            # Apply transformation to selected samples
            new_names <- sapply(input$samples_to_transform, transform_fn)

            # Update both tables
            for (i in seq_along(input$samples_to_transform)) {
                original_name <- input$samples_to_transform[i]
                new_name <- new_names[i]
                current_matrix$Run[current_matrix$Run == original_name] <- new_name
                current_data_cln$Run[current_data_cln$Run == original_name] <- new_name
            }

            # Update reactive values
            design_matrix(current_matrix)
            data_cln_reactive(current_data_cln)
        })
        
        # Handler for adding a new factor
        observeEvent(input$add_factor, {
            req(input$new_factor)
            new_factor_name <- trimws(input$new_factor)
            if (new_factor_name != "" && !new_factor_name %in% factors()) {
                factors(c(factors(), new_factor_name))
            }
            updateTextInput(session, "new_factor", value = "")
        })

        # Handler for assigning metadata to runs
        observeEvent(input$assign_metadata, {
            req(input$selected_runs, input$factor1_select)

            current_matrix <- design_matrix()
            
            # Generate replicate numbers
            replicate_numbers <- if (!is.null(input$replicate_start)) {
                seq(input$replicate_start, length.out = length(input$selected_runs))
            } else {
                NA_integer_
            }
            
            # Update the selected rows with factors and replicates
            selected_indices <- which(current_matrix$Run %in% input$selected_runs)
            current_matrix$factor1[selected_indices] <- input$factor1_select
            current_matrix$factor2[selected_indices] <- input$factor2_select
            current_matrix$replicates[selected_indices] <- replicate_numbers

            # Create group names based on factors
            group_name <- if (input$factor2_select == "") {
                input$factor1_select
            } else {
                paste(input$factor1_select, input$factor2_select, sep = "_")
            }
            current_matrix$group[selected_indices] <- group_name

            design_matrix(current_matrix)

            # Update the list of unique groups
            unique_groups <- unique(current_matrix$group[!is.na(current_matrix$group) & current_matrix$group != ""])
            groups(unique_groups)
        })

        # Handler for assigning technical replicates
        observeEvent(input$assign_tech_reps, {
            req(input$tech_rep_samples)
            current_matrix <- design_matrix()

            if (length(input$tech_rep_samples) < 2) {
                showNotification("Please select at least two samples to assign as technical replicates.", type = "warning")
                return()
            }

            # Get selected sample indices
            selected_indices <- which(current_matrix$Run %in% input$tech_rep_samples)
            
            # Check if all selected samples have the same group assignment
            selected_groups <- unique(current_matrix$group[selected_indices])
            if (length(selected_groups) > 1 || any(is.na(selected_groups))) {
                showNotification("All selected samples must belong to the same group. Please assign metadata first.", type = "warning")
                return()
            }

            # Determine the base replicate number
            if (input$tech_rep_assignment_mode == "lowest") {
                base_replicate_number <- min(current_matrix$replicates[selected_indices], na.rm = TRUE)
            } else if (input$tech_rep_assignment_mode == "first") {
                base_replicate_number <- current_matrix$replicates[current_matrix$Run == input$tech_rep_samples[1]]
            } else if (input$tech_rep_assignment_mode == "manual") {
                base_replicate_number <- input$manual_replicate_number
            }

            # Get the group of the selected samples
            sample_group <- selected_groups[1]
            
            # Find all samples in the same group
            group_indices <- which(current_matrix$group == sample_group & !is.na(current_matrix$group))
            
            # Calculate how many replicate numbers we're consolidating
            original_replicate_numbers <- unique(current_matrix$replicates[selected_indices])
            num_replicates_consolidated <- length(original_replicate_numbers) - 1
            
            # Get the highest replicate number being consolidated
            highest_consolidated_replicate <- max(original_replicate_numbers, na.rm = TRUE)
            
            # Adjust replicate numbers for samples in the same group that come after the consolidated ones
            samples_to_adjust <- which(
                current_matrix$group == sample_group & 
                !is.na(current_matrix$group) &
                current_matrix$replicates > highest_consolidated_replicate &
                !current_matrix$Run %in% input$tech_rep_samples
            )
            
            if (length(samples_to_adjust) > 0) {
                current_matrix$replicates[samples_to_adjust] <- 
                    current_matrix$replicates[samples_to_adjust] - num_replicates_consolidated
            }

            # Assign the same replicate number to all selected samples
            current_matrix$replicates[selected_indices] <- base_replicate_number
            
            # Assign technical replicate numbers (1, 2, 3, ...) using vectorized approach
            # Create a named vector mapping sample names to their technical replicate numbers
            tech_rep_mapping <- setNames(seq_along(input$tech_rep_samples), input$tech_rep_samples)
            
            # Apply the mapping to assign technical replicate numbers
            current_matrix$tech_reps[selected_indices] <- tech_rep_mapping[current_matrix$Run[selected_indices]]

            # Update reactive value
            design_matrix(current_matrix)
            
            showNotification(paste("Assigned", length(input$tech_rep_samples), 
                                 "samples as technical replicates with biological replicate number", 
                                 base_replicate_number, "and adjusted subsequent replicate numbers"), 
                           type = "message")
        })

        # Handler for adding a new contrast
        observeEvent(input$add_contrast, {
            req(input$contrast_group1, input$contrast_group2)
            g1 <- input$contrast_group1
            g2 <- input$contrast_group2

            if (g1 != "" && g2 != "" && g1 != g2) {
                contrast_name <- paste0(g1, ".vs.", g2)
                
                # Avoid adding duplicate contrasts
                if(!contrast_name %in% contrasts()$contrast_name) {
                    new_contrast <- data.frame(
                        contrast_name = contrast_name,
                        numerator = g1,
                        denominator = g2,
                        stringsAsFactors = FALSE
                    )
                    contrasts(rbind(contrasts(), new_contrast))
                }
            }
        })

        # Handler for resetting state
        observeEvent(input$reset_changes, {
            # Show confirmation modal
            showModal(modalDialog(
                title = "Confirm Reset",
                HTML(paste0("<p>This will revert <strong>", input$reset_scope, "</strong> to their initial state. This action cannot be undone.</p>")),
                footer = tagList(
                    modalButton("Cancel"),
                    actionButton(session$ns("confirm_reset"), "Reset", class = "btn-danger")
                ),
                easyClose = TRUE
            ))
        })

        # Handler for reset confirmation
        observeEvent(input$confirm_reset, {
            scope <- input$reset_scope
            state <- initial_state()

            if (scope == "all" || scope == "sample_names") {
                design_matrix(state$design_matrix)
                data_cln_reactive(state$data_cln)
            }
            if (scope == "all" || scope == "factors") {
                factors(state$factors)
                # Clear assignments if factors are reset
                current_matrix <- design_matrix()
                current_matrix$factor1 <- NA_character_
                current_matrix$factor2 <- NA_character_
                current_matrix$group <- NA_character_
                current_matrix$tech_reps <- NA_integer_
                design_matrix(current_matrix)
                groups(state$groups)
            }
            if (scope == "all" || scope == "contrasts") {
                contrasts(state$contrasts)
            }
            if (scope == "all" || scope == "formula") {
                updateTextInput(session, "formula_string", value = state$formula)
            }
            
            removeModal()
            showNotification(paste("Reset of", scope, "completed."), type = "message")
        })

        # == Final Save Action ====================================================

        observeEvent(input$save_results, {
            # This is the main action of the module. It prepares the final
            # data objects and returns them in a list via a reactive value.

            # 1. Filter design matrix to only include assigned runs
            design_matrix_final <- design_matrix() |>
                dplyr::filter(!is.na(group) & group != "")
            
            if(nrow(design_matrix_final) == 0) {
                showNotification("No samples have been assigned to groups. Please assign metadata before saving.", type="warning")
                return()
            }

            # 2. Filter main data to match the assigned runs
            assigned_runs <- design_matrix_final$Run
            data_cln_final <- data_cln_reactive() |>
                dplyr::filter(Run %in% assigned_runs)

            # 3. Prepare contrasts table with both friendly names and contrast strings
            contrast_data <- contrasts()
            contrasts_tbl <- if (!is.null(contrast_data) && nrow(contrast_data) > 0) {
                # Check if formula uses ~ 0 + group pattern (which creates groupXXX columns)
                formula_uses_group_prefix <- grepl("~ *0 *\\+ *group", input$formula_string)
                
                # Create both friendly names and actual contrast strings
                contrast_info <- lapply(1:nrow(contrast_data), function(i) {
                    # Create friendly name (replace dots with underscores for consistency)
                    friendly_name <- gsub("\\.", "_", contrast_data$contrast_name[i])
                    friendly_name <- gsub("\\.vs\\.", "_minus_", friendly_name)
                    
                    # Create actual contrast string based on formula
                    if (formula_uses_group_prefix) {
                        # If formula is ~ 0 + group, we need to prefix with "group"
                        contrast_string <- paste0(
                            "group", contrast_data$numerator[i],
                            "-group", contrast_data$denominator[i]
                        )
                    } else {
                        # Otherwise, use group names as-is
                        contrast_string <- paste0(
                            contrast_data$numerator[i],
                            "-", contrast_data$denominator[i]
                        )
                    }
                    
                    # Return in the format: friendly_name=contrast_string
                    # But for compatibility, we'll just return the contrast string
                    # The friendly name can be derived later if needed
                    list(
                        friendly_name = friendly_name,
                        contrast_string = contrast_string,
                        full_format = paste0(friendly_name, "=", contrast_string)
                    )
                })
                
                # For now, maintain compatibility by returning just the contrast strings
                # But store the full information for potential future use
                data.frame(
                    contrasts = sapply(contrast_info, function(x) x$contrast_string),
                    friendly_names = sapply(contrast_info, function(x) x$friendly_name),
                    full_format = sapply(contrast_info, function(x) x$full_format),
                    stringsAsFactors = FALSE
                )
            } else {
                NULL
            }
            
            # 4. Update config list with formula
            config_list_final <- config_list()
            config_list_final[["deAnalysisParameters"]][["formula_string"]] <- input$formula_string
            
            # 5. Set the final result list to the reactive value
            final_result <- list(
                design_matrix = design_matrix_final,
                data_cln = data_cln_final,
                contrasts_tbl = contrasts_tbl,
                config_list = config_list_final
            )
            result_rv(final_result)
            
            showNotification("Design saved successfully. You can close this builder.", type="message", duration = 5)
        })

        # The return value of the module server is the reactive containing the results
        return(result_rv)
    })
} 