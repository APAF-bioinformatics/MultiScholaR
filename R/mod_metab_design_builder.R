#' @title Metabolomics Design Matrix Builder Module
#'
#' @description A Shiny module to build a design matrix for metabolomics experiments.
#' This module is 1:1 with the proteomics design builder (mod_prot_design_builder.R),
#' with the key difference that metabolomics data is stored as a LIST of assays
#' (e.g., LCMS_Pos, LCMS_Neg) where samples are COLUMN NAMES, not row values.
#'
#' @param id The module ID.
#' @param data_tbl A reactive expression that returns the data table LIST
#'   (e.g., list(LCMS_Pos = df_pos, LCMS_Neg = df_neg)).
#' @param config_list A reactive expression that returns the main configuration list.
#' @param column_mapping A reactive expression that returns the column mapping list.
#'
#' @return A reactive expression that returns a list containing the results
#'   when the "Save" button is clicked. The list includes:
#'   - `design_matrix`: The final, filtered design matrix.
#'   - `data_cln`: The data table LIST, filtered to include only assigned samples.
#'   - `contrasts_tbl`: A data frame of the defined contrasts.
#'   - `config_list`: The updated configuration list with the new formula.
#'
#' @import shiny
#' @import shinydashboard
#' @import DT
#' @import gtools
#' @import dplyr
#' @importFrom tibble tibble
#' @name metabolomicsDesignMatrixBuilderModule
NULL

#' Metabolomics Design Matrix Builder UI Module
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
#' @rdname metabolomicsDesignMatrixBuilderModule
#' @export
mod_metab_design_builder_ui <- function(id) {
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
                        "<b>Multi-Assay Note:</b> This metabolomics workflow supports multiple assays "
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
                    shiny::h4("Current Design Matrix (All Samples) - Metabolomics")
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


#' Metabolomics Design Matrix Builder Server Module
#'
#' @param id Module ID.
#' @param data_tbl A reactive expression that returns the data table LIST.
#' @param config_list A reactive expression that returns the main configuration list.
#' @param column_mapping A reactive expression that returns the column mapping list.
#' @param existing_design_matrix Optional reactive expression that returns an existing
#'   design matrix (e.g., from "Import Existing Design"). If provided and has valid
#'   group assignments, will be used instead of creating a fresh matrix.
#' @param existing_contrasts Optional reactive expression that returns an existing
#'   contrasts table (e.g., from "Import Existing Design").
#'
#' @rdname metabolomicsDesignMatrixBuilderModule
#' @export
mod_metab_design_builder_server <- function(id, data_tbl, config_list, column_mapping
                                            , existing_design_matrix = NULL
                                            , existing_contrasts = NULL) {
    shiny::moduleServer(id, function(input, output, session) {

        # Reactive value to store the final results
        result_rv <- shiny::reactiveVal(NULL)

        # == Helper: Extract Sample Names from Metabolomics Data List =============
        # Metabolomics data is a LIST of assays where samples are COLUMN NAMES
        # All assays share the same samples (per multiomics principle)
        get_sample_columns <- function(assay_list, col_map) {
            # DEBUG66: Defensive checks
            message("DEBUG66: get_sample_columns() called")
            message(sprintf("DEBUG66: assay_list is NULL: %s, length: %d", 
                is.null(assay_list), if(is.null(assay_list)) 0 else length(assay_list)))
            message(sprintf("DEBUG66: col_map is NULL: %s", is.null(col_map)))
            
            if (is.null(assay_list) || length(assay_list) == 0) {
                message("DEBUG66: Returning empty character(0) - assay_list is NULL or empty")
                return(character(0))
            }
            
            # Defensive check: ensure first element is a data frame
            first_assay <- assay_list[[1]]
            if (!is.data.frame(first_assay)) {
                message(sprintf("DEBUG66: WARNING - first_assay is not a data frame, class: %s", 
                    class(first_assay)[1]))
                return(character(0))
            }

            # Get the first assay to extract sample columns
            all_cols <- names(first_assay)
            message(sprintf("DEBUG66: first_assay has %d columns", length(all_cols)))

            # Defensive check: handle NULL col_map
            if (is.null(col_map)) {
                message("DEBUG66: col_map is NULL, using all columns as sample candidates")
                # Fallback: try to detect numeric columns as samples
                numeric_cols <- names(first_assay)[sapply(first_assay, is.numeric)]
                message(sprintf("DEBUG66: Found %d numeric columns", length(numeric_cols)))
                return(numeric_cols)
            }
            
            # Exclude metabolite ID and annotation columns
            exclude_cols <- c(
                col_map$metabolite_id_col
                , col_map$annotation_col
            )
            exclude_cols <- exclude_cols[!is.na(exclude_cols) & nzchar(exclude_cols)]
            message(sprintf("DEBUG66: Excluding columns: %s", paste(exclude_cols, collapse = ", ")))

            # Sample columns are the remaining columns
            sample_cols <- setdiff(all_cols, exclude_cols)

            # Also filter out known non-sample columns from MS-DIAL etc.
            non_sample_patterns <- c(
                "^Alignment.ID$"
                , "^Average.Rt"
                , "^Average.Mz"
                , "^Metabolite.name$"
                , "^Adduct.type$"
                , "^Post.curation.result$"
                , "^Fill.%$"
                , "^MS/MS.spectrum$"
                , "^Reference.RT$"
                , "^Reference.m/z$"
                , "^Formula$"
                , "^Ontology$"
                , "^INCHIKEY$"
                , "^SMILES$"
                , "^Annotation.tag"
                , "^RT.matched$"
                , "^m/z.matched$"
                , "^MS/MS.matched$"
                , "^Comment$"
                , "^Manually.modified$"
                , "^Total.score$"
                , "^RT.similarity$"
                , "^Dot.product$"
                , "^Reverse.dot.product$"
                , "^Fragment.presence"
                , "^S/N.average$"
                , "^Spectrum.reference"
                , "^MS1.isotopic.spectrum$"
                , "^MS/MS.spectrum$"
            )

            for (pattern in non_sample_patterns) {
                sample_cols <- sample_cols[!grepl(pattern, sample_cols, ignore.case = TRUE)]
            }

            # Also use the column_mapping sample_columns if available
            if (!is.null(col_map$sample_columns) && length(col_map$sample_columns) > 0) {
                sample_cols <- col_map$sample_columns
                message(sprintf("DEBUG66: Using sample_columns from col_map: %d columns", 
                    length(sample_cols)))
            }
            
            message(sprintf("DEBUG66: Returning %d sample columns", length(sample_cols)))
            return(sample_cols)
        }

        # == Initial State Setup =================================================
        initial_state <- shiny::reactive({
            shiny::req(data_tbl())
            shiny::req(config_list())

            assay_list <- data_tbl()
            conf <- config_list()
            col_map <- column_mapping()

            # Get sample names from the assay list
            sample_names <- get_sample_columns(assay_list, col_map)

            if (length(sample_names) == 0) {
                logger::log_error("No sample columns detected in metabolomics data")
                return(NULL)
            }
            
            # --- Check for existing design matrix from "Import Existing Design" ---
            # If an existing design matrix is provided with valid group assignments,
            # use it instead of creating a fresh one with NA values
            use_existing <- FALSE
            existing_dm <- NULL
            existing_ctr <- NULL
            
            if (!is.null(existing_design_matrix)) {
                existing_dm <- tryCatch(existing_design_matrix(), error = function(e) NULL)
                if (!is.null(existing_dm) && is.data.frame(existing_dm) && nrow(existing_dm) > 0) {
                    # Check if the existing design has any group assignments
                    if ("group" %in% names(existing_dm)) {
                        non_na_groups <- existing_dm$group[!is.na(existing_dm$group) & existing_dm$group != ""]
                        if (length(non_na_groups) > 0) {
                            use_existing <- TRUE
                            message(sprintf("DEBUG66: Using existing design matrix with %d assigned samples", 
                                length(non_na_groups)))
                        }
                    }
                }
            }
            
            if (!is.null(existing_contrasts)) {
                existing_ctr <- tryCatch(existing_contrasts(), error = function(e) NULL)
            }
            
            if (use_existing && !is.null(existing_dm)) {
                # Use the imported design matrix
                message("DEBUG66: initial_state using IMPORTED design matrix")
                
                # Extract groups and factors from the existing design
                groups_from_dm <- unique(existing_dm$group[!is.na(existing_dm$group) & existing_dm$group != ""])
                factors_from_dm <- character(0)
                if ("factor1" %in% names(existing_dm)) {
                    f1_vals <- unique(existing_dm$factor1[!is.na(existing_dm$factor1) & existing_dm$factor1 != ""])
                    factors_from_dm <- c(factors_from_dm, f1_vals)
                }
                if ("factor2" %in% names(existing_dm)) {
                    f2_vals <- unique(existing_dm$factor2[!is.na(existing_dm$factor2) & existing_dm$factor2 != ""])
                    factors_from_dm <- c(factors_from_dm, f2_vals)
                }
                factors_from_dm <- unique(factors_from_dm)
                
                # Convert existing contrasts to internal format
                contrasts_df <- if (!is.null(existing_ctr) && is.data.frame(existing_ctr) && nrow(existing_ctr) > 0) {
                    # Parse the contrast strings to extract numerator/denominator
                    parsed <- lapply(existing_ctr$contrasts, function(cs) {
                        # Expected format: "groupTreatment-groupControl"
                        parts <- strsplit(gsub("^group", "", cs), "-group")[[1]]
                        if (length(parts) == 2) {
                            list(
                                contrast_name = paste0(parts[1], ".vs.", parts[2])
                                , numerator = parts[1]
                                , denominator = parts[2]
                            )
                        } else {
                            # Fallback
                            list(
                                contrast_name = cs
                                , numerator = cs
                                , denominator = ""
                            )
                        }
                    })
                    data.frame(
                        contrast_name = sapply(parsed, function(x) x$contrast_name)
                        , numerator = sapply(parsed, function(x) x$numerator)
                        , denominator = sapply(parsed, function(x) x$denominator)
                        , stringsAsFactors = FALSE
                    )
                } else {
                    data.frame(
                        contrast_name = character()
                        , numerator = character()
                        , denominator = character()
                        , stringsAsFactors = FALSE
                    )
                }
                
                return(list(
                    design_matrix = existing_dm
                    , data_cln = assay_list
                    , sample_names = sample_names
                    , groups = groups_from_dm
                    , factors = factors_from_dm
                    , formula = conf[["deAnalysisParameters"]][["formula_string"]]
                    , contrasts = contrasts_df
                ))
            }

            # --- Create fresh design matrix (standard case) ---
            message("DEBUG66: initial_state creating FRESH design matrix")
            
            # Initialize design matrix with factor columns
            # Use "Run" to match proteomics structure (sample identifier)
            design_matrix_raw <- tibble::tibble(
                Run = gtools::mixedsort(sample_names)
                , group = NA_character_
                , factor1 = NA_character_
                , factor2 = NA_character_
                , factor3 = NA_character_
                , batch = NA_character_
                , replicates = NA_integer_
                , tech_reps = NA_integer_
            )

            list(
                design_matrix = design_matrix_raw
                , data_cln = assay_list
                , sample_names = sample_names
                , groups = unique(design_matrix_raw$group[!is.na(design_matrix_raw$group) & design_matrix_raw$group != ""])
                , factors = character(0)
                , formula = conf[["deAnalysisParameters"]][["formula_string"]]
                , contrasts = data.frame(
                    contrast_name = character()
                    , numerator = character()
                    , denominator = character()
                    , stringsAsFactors = FALSE
                )
            )
        })

        # Reactive values for the current, mutable state
        design_matrix <- shiny::reactiveVal()
        data_cln_reactive <- shiny::reactiveVal()
        sample_names_reactive <- shiny::reactiveVal()
        groups <- shiny::reactiveVal()
        factors <- shiny::reactiveVal()
        contrasts <- shiny::reactiveVal()
        removed_samples <- shiny::reactiveVal(character(0))

        # Observer to initialize/reset from initial state
        shiny::observe({
            state <- initial_state()
            shiny::req(state)

            design_matrix(state$design_matrix)
            data_cln_reactive(state$data_cln)
            sample_names_reactive(state$sample_names)
            groups(state$groups)
            factors(state$factors)
            contrasts(state$contrasts)
            removed_samples(character(0))

            # Update UI elements
            sorted_runs <- state$design_matrix$Run
            shiny::updateSelectizeInput(session, "sample_to_rename", choices = sorted_runs, selected = "")
            shiny::updateSelectizeInput(session, "selected_runs", choices = sorted_runs, selected = "")
            shiny::updateSelectizeInput(session, "samples_to_transform", choices = sorted_runs, selected = "")
            shiny::updateSelectizeInput(session, "tech_rep_samples", choices = sorted_runs, selected = "")
            shiny::updateSelectizeInput(session, "samples_to_remove", choices = sorted_runs, selected = "")
            shiny::updateTextInput(session, "formula_string", value = state$formula)
        })

        # Create a proxy for the main data table
        proxy_data_table <- DT::dataTableProxy("data_table")

        # == UI Rendering and Updates =============================================

        # Update sample selection inputs when names change or samples are removed
        shiny::observe({
            shiny::req(design_matrix())
            all_runs <- gtools::mixedsort(design_matrix()$Run)

            currently_removed <- removed_samples()
            available_runs <- all_runs[!all_runs %in% currently_removed]

            shiny::isolate({
                selected_rename <- input$sample_to_rename
                selected_meta <- input$selected_runs
                selected_transform <- input$samples_to_transform
                selected_tech <- input$tech_rep_samples
                selected_remove <- input$samples_to_remove

                selected_rename <- if (!is.null(selected_rename) && selected_rename %in% available_runs) selected_rename else ""
                selected_meta <- selected_meta[selected_meta %in% available_runs]
                selected_transform <- selected_transform[selected_transform %in% available_runs]
                selected_tech <- selected_tech[selected_tech %in% available_runs]
                selected_remove <- selected_remove[selected_remove %in% available_runs]
            })

            shiny::updateSelectizeInput(session, "sample_to_rename", choices = available_runs, selected = selected_rename)
            shiny::updateSelectizeInput(session, "selected_runs", choices = available_runs, selected = selected_meta)
            shiny::updateSelectizeInput(session, "samples_to_transform", choices = available_runs, selected = selected_transform)
            shiny::updateSelectizeInput(session, "tech_rep_samples", choices = available_runs, selected = selected_tech)
            shiny::updateSelectizeInput(session, "samples_to_remove", choices = available_runs, selected = selected_remove)
        })

        # Update factor and group dropdowns
        shiny::observe({
            shiny::isolate({
                f1 <- input$factor1_select
                f2 <- input$factor2_select
                f3 <- input$factor3_select
                g1 <- input$contrast_group1
                g2 <- input$contrast_group2
            })

            current_factors <- factors()
            current_groups <- groups()

            shiny::updateSelectInput(session, "factor1_select", choices = c("", current_factors), selected = f1)
            shiny::updateSelectInput(session, "factor2_select", choices = c("", current_factors), selected = f2)
            shiny::updateSelectInput(session, "factor3_select", choices = c("", current_factors), selected = f3)
            shiny::updateSelectInput(session, "contrast_group1", choices = c("", current_groups), selected = g1)
            shiny::updateSelectInput(session, "contrast_group2", choices = c("", current_groups), selected = g2)
        })

        # Render the main data table (excluding removed samples)
        output$data_table <- DT::renderDT({
                shiny::req(design_matrix())
                dm <- design_matrix()
                currently_removed <- removed_samples()
                dm |> dplyr::filter(!Run %in% currently_removed)
            }
            , selection = "none"
            , options = list(pageLength = 10, scrollX = TRUE, server = FALSE)
        )

        # Use proxy to update table data
        shiny::observe({
            shiny::req(proxy_data_table)
            shiny::req(design_matrix())
            dm <- design_matrix()
            currently_removed <- removed_samples()
            filtered_dm <- dm |> dplyr::filter(!Run %in% currently_removed)
            DT::replaceData(proxy_data_table, filtered_dm, resetPaging = FALSE)
        })

        # Render Available Factors Display
        output$available_factors_display <- shiny::renderUI({
            current_factors <- factors()
            if (length(current_factors) == 0) {
                return(shiny::p("No factors defined yet (use the 'Factors' tab)."))
            }
            shiny::p(paste(current_factors, collapse = ", "))
        })

        # Render Defined Contrasts Display
        output$defined_contrasts_display <- shiny::renderUI({
            contrast_data <- contrasts()
            if (is.null(contrast_data) || nrow(contrast_data) == 0) {
                return(shiny::p("No contrasts defined yet."))
            }

            formula_uses_group_prefix <- grepl("~ *0 *\\+ *group", input$formula_string)

            contrast_info <- lapply(1:nrow(contrast_data), function(i) {
                friendly_name <- gsub("\\.", "_", contrast_data$contrast_name[i])

                if (formula_uses_group_prefix) {
                    contrast_string <- paste0(
                        "group", contrast_data$numerator[i]
                        , "-group", contrast_data$denominator[i]
                    )
                } else {
                    contrast_string <- paste0(
                        contrast_data$numerator[i]
                        , "-", contrast_data$denominator[i]
                    )
                }

                paste0(friendly_name, "=", contrast_string)
            })

            shiny::tagList(
                lapply(contrast_info, function(full_format) {
                    shiny::p(shiny::tags$code(full_format))
                })
            )
        })

        # Render preview for range extraction
        output$range_preview <- shiny::renderText({
            shiny::req(input$samples_to_transform)
            shiny::req(input$range_start, input$range_end)

            first_sample <- input$samples_to_transform[1]

            tryCatch({
                preview_result <- extract_experiment(
                    first_sample
                    , mode = "range"
                    , start = input$range_start
                    , end = input$range_end
                )
                paste0("\"", first_sample, "\" -> \"", preview_result, "\"")
            }, error = function(e) {
                paste("Error:", e$message)
            })
        })

        # Render technical replicate summary
        output$tech_rep_summary <- shiny::renderText({
            dm <- design_matrix()
            shiny::req(dm)

            tech_rep_groups <- dm |>
                dplyr::filter(!is.na(tech_reps)) |>
                dplyr::group_by(group, replicates) |>
                dplyr::summarise(
                    samples = paste(Run, collapse = ", ")
                    , tech_rep_numbers = paste(tech_reps, collapse = ", ")
                    , .groups = "drop"
                )

            if (nrow(tech_rep_groups) == 0) {
                return("No technical replicates assigned yet.")
            }

            output_text <- tech_rep_groups |>
                dplyr::mutate(
                    formatted = sprintf("Group: %s, Biological Replicate: %s\n  Samples: %s\n  Technical Replicates: %s"
                                      , group, replicates, samples, tech_rep_numbers)
                ) |>
                dplyr::pull(formatted)

            paste(output_text, collapse = "\n\n")
        })

        # Render removed samples display
        output$removed_samples_display <- shiny::renderText({
            currently_removed <- removed_samples()
            if (length(currently_removed) == 0) {
                return("No samples have been removed.")
            }
            paste0(
                "Removed ", length(currently_removed), " sample(s):\n"
                , paste(gtools::mixedsort(currently_removed), collapse = "\n")
            )
        })

        # Render replicate number input UI
        output$replicate_inputs <- shiny::renderUI({
            ns <- session$ns
            shiny::req(input$selected_runs)
            shiny::numericInput(ns("replicate_start")
                , paste("Starting replicate number for", length(input$selected_runs), "selected samples:")
                , value = 1
                , min = 1
            )
        })

        # Render contrast factors info
        output$contrast_factors_info <- shiny::renderText({
            formula_uses_group_prefix <- grepl("~ *0 *\\+ *group", input$formula_string)
            if (formula_uses_group_prefix) {
                "Note: Contrasts will use 'group' prefix based on formula: ~ 0 + group"
            } else {
                "Note: Contrasts will use group names as-is"
            }
        })

        # == Event Handlers =======================================================

        # Handler for individual sample renaming
        shiny::observeEvent(input$rename_sample, {
            shiny::req(input$sample_to_rename, input$new_sample_name)

            old_name <- input$sample_to_rename
            new_name <- input$new_sample_name

            # Update design matrix
            current_matrix <- design_matrix()
            current_matrix$Run[current_matrix$Run == old_name] <- new_name
            design_matrix(current_matrix)

            # Update data_cln - rename columns in ALL assays
            current_data <- data_cln_reactive()
            for (assay_name in names(current_data)) {
                col_idx <- which(names(current_data[[assay_name]]) == old_name)
                if (length(col_idx) > 0) {
                    names(current_data[[assay_name]])[col_idx] <- new_name
                }
            }
            data_cln_reactive(current_data)

            # Update sample names tracking
            current_samples <- sample_names_reactive()
            current_samples[current_samples == old_name] <- new_name
            sample_names_reactive(current_samples)

            shiny::updateTextInput(session, "new_sample_name", value = "")
        })

        # Handler for bulk renaming
        shiny::observeEvent(input$bulk_rename, {
            shiny::req(input$samples_to_transform)

            current_matrix <- design_matrix()
            current_data <- data_cln_reactive()
            current_samples <- sample_names_reactive()

            transform_fn <- function(sample_name) {
                if (input$transform_mode == "range") {
                    shiny::req(input$range_start, input$range_end)
                    extract_experiment(sample_name, mode = "range", start = input$range_start, end = input$range_end)
                } else if (input$transform_mode == "before_underscore") {
                    extract_experiment(sample_name, mode = "start")
                } else if (input$transform_mode == "after_underscore") {
                    extract_experiment(sample_name, mode = "end")
                }
            }

            new_names <- sapply(input$samples_to_transform, transform_fn)

            for (i in seq_along(input$samples_to_transform)) {
                original_name <- input$samples_to_transform[i]
                new_name <- new_names[i]

                # Update design matrix
                current_matrix$Run[current_matrix$Run == original_name] <- new_name

                # Update ALL assays
                for (assay_name in names(current_data)) {
                    col_idx <- which(names(current_data[[assay_name]]) == original_name)
                    if (length(col_idx) > 0) {
                        names(current_data[[assay_name]])[col_idx] <- new_name
                    }
                }

                # Update sample names tracking
                current_samples[current_samples == original_name] <- new_name
            }

            design_matrix(current_matrix)
            data_cln_reactive(current_data)
            sample_names_reactive(current_samples)
        })

        # Handler for adding a new factor
        shiny::observeEvent(input$add_factor, {
            shiny::req(input$new_factor)
            new_factor_name <- trimws(input$new_factor)
            if (new_factor_name != "" && !new_factor_name %in% factors()) {
                factors(c(factors(), new_factor_name))
            }
            shiny::updateTextInput(session, "new_factor", value = "")
        })

        # Handler for assigning metadata
        shiny::observeEvent(input$assign_metadata, {
            shiny::req(input$selected_runs, input$factor1_select)

            current_matrix <- design_matrix()

            replicate_numbers <- if (!is.null(input$replicate_start)) {
                seq(input$replicate_start, length.out = length(input$selected_runs))
            } else {
                NA_integer_
            }

            selected_indices <- which(current_matrix$Run %in% input$selected_runs)
            current_matrix$factor1[selected_indices] <- input$factor1_select
            current_matrix$factor2[selected_indices] <- input$factor2_select
            current_matrix$factor3[selected_indices] <- if (!is.null(input$factor3_select) && input$factor3_select != "") {
                input$factor3_select
            } else {
                NA_character_
            }
            current_matrix$replicates[selected_indices] <- replicate_numbers

            # Create group names from factors
            factor_parts <- c()
            if (!is.null(input$factor1_select) && input$factor1_select != "") {
                factor_parts <- c(factor_parts, input$factor1_select)
            }
            if (!is.null(input$factor2_select) && input$factor2_select != "") {
                factor_parts <- c(factor_parts, input$factor2_select)
            }
            if (!is.null(input$factor3_select) && input$factor3_select != "") {
                factor_parts <- c(factor_parts, input$factor3_select)
            }

            group_name <- if (length(factor_parts) > 0) {
                paste(factor_parts, collapse = "_")
            } else {
                NA_character_
            }
            current_matrix$group[selected_indices] <- group_name

            design_matrix(current_matrix)

            unique_groups <- unique(current_matrix$group[!is.na(current_matrix$group) & current_matrix$group != ""])
            groups(unique_groups)
        })

        # Handler for assigning technical replicates
        shiny::observeEvent(input$assign_tech_reps, {
            shiny::req(input$tech_rep_samples)
            current_matrix <- design_matrix()

            if (length(input$tech_rep_samples) < 2) {
                shiny::showNotification("Please select at least two samples.", type = "warning")
                return()
            }

            selected_indices <- which(current_matrix$Run %in% input$tech_rep_samples)
            selected_groups <- unique(current_matrix$group[selected_indices])

            if (length(selected_groups) > 1 || any(is.na(selected_groups))) {
                shiny::showNotification("All selected samples must belong to the same group.", type = "warning")
                return()
            }

            if (input$tech_rep_assignment_mode == "lowest") {
                base_replicate_number <- min(current_matrix$replicates[selected_indices], na.rm = TRUE)
            } else if (input$tech_rep_assignment_mode == "first") {
                base_replicate_number <- current_matrix$replicates[current_matrix$Run == input$tech_rep_samples[1]]
            } else if (input$tech_rep_assignment_mode == "manual") {
                base_replicate_number <- input$manual_replicate_number
            }

            sample_group <- selected_groups[1]
            original_replicate_numbers <- unique(current_matrix$replicates[selected_indices])
            num_replicates_consolidated <- length(original_replicate_numbers) - 1
            highest_consolidated_replicate <- max(original_replicate_numbers, na.rm = TRUE)

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

            current_matrix$replicates[selected_indices] <- base_replicate_number
            tech_rep_mapping <- setNames(seq_along(input$tech_rep_samples), input$tech_rep_samples)
            current_matrix$tech_reps[selected_indices] <- tech_rep_mapping[current_matrix$Run[selected_indices]]

            design_matrix(current_matrix)

            shiny::showNotification(paste("Assigned", length(input$tech_rep_samples), "samples as technical replicates"), type = "message")
        })

        # Handler for adding a new contrast
        shiny::observeEvent(input$add_contrast, {
            shiny::req(input$contrast_group1, input$contrast_group2)
            g1 <- input$contrast_group1
            g2 <- input$contrast_group2

            if (g1 != "" && g2 != "" && g1 != g2) {
                contrast_name <- paste0(g1, ".vs.", g2)

                if(!contrast_name %in% contrasts()$contrast_name) {
                    new_contrast <- data.frame(
                        contrast_name = contrast_name
                        , numerator = g1
                        , denominator = g2
                        , stringsAsFactors = FALSE
                    )
                    contrasts(rbind(contrasts(), new_contrast))
                }
            }
        })

        # Handler for removing samples
        shiny::observeEvent(input$remove_samples, {
            shiny::req(input$samples_to_remove)

            samples_to_remove <- input$samples_to_remove
            currently_removed <- removed_samples()

            new_removed <- unique(c(currently_removed, samples_to_remove))
            removed_samples(new_removed)

            shiny::updateSelectizeInput(session, "samples_to_remove", selected = "")

            shiny::showNotification(
                paste("Removed", length(samples_to_remove), "sample(s) from analysis.")
                , type = "message"
            )
        })

        # Handler for resetting state
        shiny::observeEvent(input$reset_changes, {
            shiny::showModal(shiny::modalDialog(
                title = "Confirm Reset"
                , shiny::HTML(paste0("<p>This will revert <strong>", input$reset_scope, "</strong> to their initial state.</p>"))
                , footer = shiny::tagList(
                    shiny::modalButton("Cancel")
                    , shiny::actionButton(session$ns("confirm_reset"), "Reset", class = "btn-danger")
                )
                , easyClose = TRUE
            ))
        })

        # Handler for reset confirmation
        shiny::observeEvent(input$confirm_reset, {
            scope <- input$reset_scope
            state <- initial_state()

            if (scope == "all" || scope == "sample_names") {
                design_matrix(state$design_matrix)
                data_cln_reactive(state$data_cln)
                sample_names_reactive(state$sample_names)
            }
            if (scope == "all" || scope == "removed_samples") {
                removed_samples(character(0))
            }
            if (scope == "all" || scope == "factors") {
                factors(state$factors)
                current_matrix <- design_matrix()
                current_matrix$factor1 <- NA_character_
                current_matrix$factor2 <- NA_character_
                current_matrix$factor3 <- NA_character_
                current_matrix$group <- NA_character_
                current_matrix$tech_reps <- NA_integer_
                design_matrix(current_matrix)
                groups(state$groups)
            }
            if (scope == "all" || scope == "contrasts") {
                contrasts(state$contrasts)
            }
            if (scope == "all" || scope == "formula") {
                shiny::updateTextInput(session, "formula_string", value = state$formula)
            }

            shiny::removeModal()
            shiny::showNotification(paste("Reset of", scope, "completed."), type = "message")
        })

        # == Final Save Action ====================================================

        shiny::observeEvent(input$save_results, {
            # 1. Filter design matrix to only assigned samples, excluding removed
            currently_removed <- removed_samples()
            design_matrix_final <- design_matrix() |>
                dplyr::filter(!is.na(group) & group != "") |>
                dplyr::filter(!Run %in% currently_removed) |>
                dplyr::mutate(tech_rep_group = paste(group, replicates, sep = "_"))

            if(nrow(design_matrix_final) == 0) {
                shiny::showNotification("No samples have been assigned to groups.", type = "warning")
                return()
            }

            # 2. Filter data_cln (all assays) to match assigned samples
            assigned_samples <- design_matrix_final$Run
            current_data <- data_cln_reactive()
            col_map <- column_mapping()

            data_cln_final <- lapply(current_data, function(assay_df) {
                # Keep metabolite ID column and annotation column
                keep_cols <- c(
                    col_map$metabolite_id_col
                    , col_map$annotation_col
                )
                keep_cols <- keep_cols[!is.na(keep_cols) & nzchar(keep_cols)]
                keep_cols <- keep_cols[keep_cols %in% names(assay_df)]

                # Add assigned sample columns
                sample_cols_in_assay <- intersect(assigned_samples, names(assay_df))
                final_cols <- unique(c(keep_cols, sample_cols_in_assay))

                assay_df[, final_cols, drop = FALSE]
            })

            # 3. Prepare contrasts table
            contrast_data <- contrasts()
            contrasts_tbl <- if (!is.null(contrast_data) && nrow(contrast_data) > 0) {
                formula_uses_group_prefix <- grepl("~ *0 *\\+ *group", input$formula_string)

                contrast_info <- lapply(1:nrow(contrast_data), function(i) {
                    friendly_name <- gsub("\\.", "_", contrast_data$contrast_name[i])

                    if (formula_uses_group_prefix) {
                        contrast_string <- paste0(
                            "group", contrast_data$numerator[i]
                            , "-group", contrast_data$denominator[i]
                        )
                    } else {
                        contrast_string <- paste0(
                            contrast_data$numerator[i]
                            , "-", contrast_data$denominator[i]
                        )
                    }

                    list(
                        friendly_name = friendly_name
                        , contrast_string = contrast_string
                        , full_format = paste0(friendly_name, "=", contrast_string)
                    )
                })

                data.frame(
                    contrasts = sapply(contrast_info, function(x) x$contrast_string)
                    , friendly_names = sapply(contrast_info, function(x) x$friendly_name)
                    , full_format = sapply(contrast_info, function(x) x$full_format)
                    , stringsAsFactors = FALSE
                )
            } else {
                NULL
            }

            # 4. Update config list with formula
            config_list_final <- config_list()
            config_list_final[["deAnalysisParameters"]][["formula_string"]] <- input$formula_string

            # 5. Set the final result
            final_result <- list(
                design_matrix = design_matrix_final
                , data_cln = data_cln_final
                , contrasts_tbl = contrasts_tbl
                , config_list = config_list_final
            )
            result_rv(final_result)

            shiny::showNotification("Design saved successfully!", type = "message", duration = 5)
        })

        # Return the reactive containing results
        return(result_rv)
    })
}
