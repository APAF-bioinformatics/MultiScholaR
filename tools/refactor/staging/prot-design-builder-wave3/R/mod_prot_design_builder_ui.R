#' Design Matrix Builder UI Module
#'
#' @param id Module ID.
#' @return A fluidPage containing the UI for the design matrix builder.
#' 
#' @note **ARCHITECTURAL EXCEPTION**: This module is designed as a standalone component
#' callable from R Markdown workflows, not integrated into the main Shiny app workflow.
#' Therefore, it returns `fluidPage()` instead of `tagList()`, and uses reactive expressions
#' instead of `workflow_data` for state management. This is an intentional deviation from
#' the standard Golem module pattern.
#' 
#' @rdname designMatrixBuilderModule
#' @export
mod_prot_design_builder_ui <- function(id) {
    ns <- shiny::NS(id)
    # The UI is extensive, so it is defined here. It is based on the
    # createDesignMatrixUi function from the original shiny_applets.R file.
    shiny::fluidPage(
        shiny::titlePanel("Design Matrix Builder"),

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
                            shiny::helpText(shiny::HTML(
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

                        # Remove samples tab
                        shiny::tabPanel(
                            "Remove Samples",
                            shiny::h4("Remove Samples from Analysis"),
                            shiny::helpText(shiny::HTML(
                                "Select samples you want to exclude from the analysis. ",
                                "Removed samples will not appear in the design matrix or be included in downstream analysis. ",
                                "You can restore removed samples using the 'Settings' tab reset functionality."
                            )),
                            shiny::selectizeInput(ns("samples_to_remove"), "Select Samples to Remove:",
                                choices = NULL, # Will be populated by server
                                multiple = TRUE
                            ),
                            shiny::actionButton(ns("remove_samples"), "Remove Selected Samples",
                                class = "btn-danger",
                                icon = shiny::icon("trash")
                            ),
                            shiny::hr(),
                            shiny::h5("Currently Removed Samples"),
                            shiny::verbatimTextOutput(ns("removed_samples_display"))
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
                            shiny::selectInput(ns("factor3_select"), "Select Factor 3:", choices = c("")),
                            shiny::uiOutput(ns("replicate_inputs")),
                            shiny::actionButton(ns("assign_metadata"), "Assign")
                        ),

                        # Technical replicates tab
                        shiny::tabPanel(
                            "Technical Replicates",
                            shiny::h4("Assign Technical Replicates"),
                            shiny::helpText(shiny::HTML(
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
                                        "Removed Samples Only" = "removed_samples",
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
                                    icon = shiny::icon("rotate-left"),
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
                    shiny::helpText(shiny::HTML(
                        "<b>Factors:</b> These represent the main experimental variables or conditions ",
                        "that distinguish your samples (e.g., Treatment, Timepoint, Genotype, Condition). ",
                        "They are used to group runs for analysis. Define potential factor names in the 'Factors' tab, ",
                        "then assign the appropriate factor levels to your runs using the 'Assign Metadata' tab. ",
                        "A 'group' is automatically created based on the combination of assigned factors (e.g., 'TreatmentA_Timepoint1_ConditionX'). ",
                        "You can assign up to three factors (Factor 1, Factor 2, and Factor 3), and groups will be created from their combinations."
                    )),
                    shiny::helpText(shiny::HTML(
                        "<b>Technical Replicates:</b> These are repeated measurements of the same biological sample ",
                        "(e.g., multiple injections or runs of the same sample). Use the 'Technical Replicates' tab to group ",
                        "samples that are technical replicates. They will share the same biological replicate number but have ",
                        "distinct technical replicate numbers. This information is used downstream for missing value imputation."
                    )),
                    shiny::helpText(shiny::HTML(
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

