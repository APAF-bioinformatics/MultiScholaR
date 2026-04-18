#' Lipidomics Design Matrix Builder UI Module
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
#' @rdname lipidomicsDesignMatrixBuilderModule
#' @export
mod_lipid_design_builder_ui <- function(id) {
    ns <- shiny::NS(id)
    # UI is identical to proteomics builder - same UX
    shiny::fluidPage(
        shiny::titlePanel("Design Matrix Builder")

        , shiny::fluidRow(
            # Management tools and Info panel on the left
            shiny::column(
                4
                , shiny::wellPanel(
                    shiny::tabsetPanel(
                        id = ns("main_tabs")
                        # Sample renaming tab
                        , shiny::tabPanel(
                            "Rename Samples"
                            , shiny::h4("Individual Rename")
                            , shiny::fluidRow(
                                shiny::column(6, shiny::selectizeInput(ns("sample_to_rename"), "Select Sample:"
                                    , choices = NULL
                                ))
                                , shiny::column(6, shiny::textInput(ns("new_sample_name"), "New Name:"))
                            )
                            , shiny::actionButton(ns("rename_sample"), "Rename")
                            , shiny::hr()
                            , shiny::h4("Bulk Rename")
                            , shiny::selectizeInput(ns("samples_to_transform"), "Select Samples:"
                                , choices = NULL
                                , multiple = TRUE
                            )
                            , shiny::radioButtons(ns("transform_mode"), "Transformation Mode:"
                                , choices = c(
                                    "Keep Text BEFORE First Underscore" = "before_underscore"
                                    , "Keep Text AFTER Last Underscore" = "after_underscore"
                                    , "Extract by Position" = "range"
                                )
                            )
                            , shiny::helpText(shiny::HTML(
                                "Select how to shorten the selected sample names using underscores as delimiters."
                                , "<br><b>Keep Text BEFORE First Underscore:</b> Retains only the part before the very first '_'. Example: 'SampleA_T1_Rep1' &rarr; 'SampleA'."
                                , "<br><b>Keep Text AFTER Last Underscore:</b> Retains only the part after the very last '_'. Example: 'SampleA_T1_Rep1' &rarr; 'Rep1'."
                                , "<br><b>Extract by Underscore Position:</b> Keeps text segments based on underscore indices (1-based)."
                            ))
                            , shiny::conditionalPanel(
                                condition = paste0("input['", ns("transform_mode"), "'] == 'range'")
                                , shiny::numericInput(ns("range_start"), "Start Underscore Index:", 1, min = 1)
                                , shiny::numericInput(ns("range_end"), "End Underscore Index:", 1, min = 1)
                                , shiny::h5("Preview (First Selected Sample):")
                                , shiny::verbatimTextOutput(ns("range_preview"))
                            )
                            , shiny::actionButton(ns("bulk_rename"), "Apply Transformation")
                        )

                        # Remove samples tab
                        , shiny::tabPanel(
                            "Remove Samples"
                            , shiny::h4("Remove Samples from Analysis")
                            , shiny::helpText(shiny::HTML(
                                "Select samples you want to exclude from the analysis. "
                                , "Removed samples will not appear in the design matrix or be included in downstream analysis."
                            ))
                            , shiny::selectizeInput(ns("samples_to_remove"), "Select Samples to Remove:"
                                , choices = NULL
                                , multiple = TRUE
                            )
                            , shiny::actionButton(ns("remove_samples"), "Remove Selected Samples"
                                , class = "btn-danger"
                                , icon = shiny::icon("trash")
                            )
                            , shiny::hr()
                            , shiny::h5("Currently Removed Samples")
                            , shiny::verbatimTextOutput(ns("removed_samples_display"))
                        )

                        # Factor management tab
                        , shiny::tabPanel(
                            "Factors"
                            , shiny::h4("Add New Factor")
                            , shiny::textInput(ns("new_factor"), "New Factor Name:")
                            , shiny::actionButton(ns("add_factor"), "Add Factor")
                        )

                        # Metadata assignment tab
                        , shiny::tabPanel(
                            "Assign Metadata"
                            , shiny::h4("Assign Metadata")
                            , shiny::selectizeInput(ns("selected_runs"), "Select Samples:"
                                , choices = NULL
                                , multiple = TRUE
                            )
                            , shiny::selectInput(ns("factor1_select"), "Select Factor 1:", choices = c(""))
                            , shiny::selectInput(ns("factor2_select"), "Select Factor 2:", choices = c(""))
                            , shiny::selectInput(ns("factor3_select"), "Select Factor 3:", choices = c(""))
                            , shiny::uiOutput(ns("replicate_inputs"))
                            , shiny::actionButton(ns("assign_metadata"), "Assign")
                        )

                        # Technical replicates tab
                        , shiny::tabPanel(
                            "Technical Replicates"
                            , shiny::h4("Assign Technical Replicates")
                            , shiny::helpText(shiny::HTML(
                                "Select samples that are technical replicates of each other. "
                                , "These samples will be grouped under the same biological replicate number."
                            ))
                            , shiny::selectizeInput(ns("tech_rep_samples"), "Select Technical Replicates:"
                                , choices = NULL
                                , multiple = TRUE
                            )
                            , shiny::radioButtons(ns("tech_rep_assignment_mode"), "Assignment Mode:"
                                , choices = c(
                                    "Use lowest replicate number" = "lowest"
                                    , "Use first selected sample's replicate" = "first"
                                    , "Specify replicate number" = "manual"
                                )
                            )
                            , shiny::conditionalPanel(
                                condition = paste0("input['", ns("tech_rep_assignment_mode"), "'] == 'manual'")
                                , shiny::numericInput(ns("manual_replicate_number"), "Replicate Number:",
                                    value = 1, min = 1)
                            )
                            , shiny::actionButton(ns("assign_tech_reps"), "Assign Technical Replicates"
                                , class = "btn-primary")
                            , shiny::hr()
                            , shiny::h5("Current Technical Replicate Assignments")
                            , shiny::verbatimTextOutput(ns("tech_rep_summary"))
                        )

                        # Contrast tab
                        , shiny::tabPanel(
                            "Contrasts"
                            , shiny::h4("Define Contrasts")
                            , shiny::selectInput(ns("contrast_group1"), "Group 1:", choices = c(""))
                            , shiny::selectInput(ns("contrast_group2"), "Group 2:", choices = c(""))
                            , shiny::verbatimTextOutput(ns("contrast_factors_info"))
                            , shiny::actionButton(ns("add_contrast"), "Add Contrast")
                        )

                        # Formula tab
                        , shiny::tabPanel(
                            "Formula"
                            , shiny::h4("Model Formula")
                            , shiny::textInput(ns("formula_string"), "Formula:"
                                , value = "~ 0 + group"
                            )
                            , shiny::helpText("Default: '~ 0 + group' creates groupX-groupY contrast format for limma.")
                        )

                        # Settings tab for global operations
                        , shiny::tabPanel(
                            "Settings"
                            , shiny::h4("Session Management")
                            , shiny::p("Reset options allow you to revert changes made during this session.")
                            , shiny::div(
                                style = "margin-top: 15px;"
                                , shiny::selectInput(ns("reset_scope"), "Reset Scope:"
                                    , choices = c(
                                        "All Changes" = "all"
                                        , "Sample Names Only" = "sample_names"
                                        , "Removed Samples Only" = "removed_samples"
                                        , "Factors Only" = "factors"
                                        , "Contrasts Only" = "contrasts"
                                        , "Formula Only" = "formula"
                                    )
                                    , selected = "all"
                                )
                            )
                            , shiny::div(
                                style = "margin-top: 15px;"
                                , shiny::actionButton(ns("reset_changes"), "Reset to Initial State"
                                    , class = "btn-warning"
                                    , icon = shiny::icon("rotate-left")
                                    , width = "100%"
                                )
                            )
                            , shiny::hr()
                            , shiny::h4("Information")
                            , shiny::verbatimTextOutput(ns("session_info"))
                        )
                    )
                )

                # Information & Assignments Panel
                , shiny::wellPanel(
                    shiny::h4("Information & Assignments")
                    , shiny::helpText(shiny::HTML(
                        "<b>Multi-Assay Note:</b> This lipidomics workflow supports multiple assays "
                        , "(e.g., LCMS_Pos, LCMS_Neg, GCMS). All assays share the SAME design matrix - "
                        , "sample assignments, factors, and contrasts apply uniformly across all assays."
                    ))
                    , shiny::hr()
                    , shiny::helpText(shiny::HTML(
                        "<b>Factors:</b> Experimental variables distinguishing your samples "
                        , "(e.g., Treatment, Timepoint, Genotype)."
                    ))
                    , shiny::helpText(shiny::HTML(
                        "<b>Technical Replicates:</b> Repeated measurements of the same biological sample."
                    ))
                    , shiny::helpText(shiny::HTML(
                        "<b>Contrasts:</b> Statistical comparisons between groups for differential analysis."
                    ))
                    , shiny::hr()
                    , shiny::h5("Available Factors")
                    , shiny::uiOutput(ns("available_factors_display"))
                    , shiny::hr()
                    , shiny::h5("Defined Contrasts")
                    , shiny::uiOutput(ns("defined_contrasts_display"))
                )
            )

            # Design matrix display on the right
            , shiny::column(
                8
                , shiny::wellPanel(
                    shiny::h4("Current Design Matrix (All Samples) - Lipidomics")
                    , DT::DTOutput(ns("data_table"))
                )
                , shiny::actionButton(ns("save_results"), "Save Design"
                    , class = "btn-primary"
                    , icon = shiny::icon("save")
                    , style = "float: right;"
                )
            )
        )
    )
}

