# ============================================================================
# mod_lipid_de.R
# ============================================================================
# Purpose: Lipidomics differential abundance analysis Shiny module
#
# This module provides per-assay differential abundance analysis using limma,
# with interactive volcano plots (Glimma), heatmaps, and result tables.
# UI and functionality matches the proteomics DE module (mod_prot_de.R).
# ============================================================================

#' @title Lipidomics Differential Analysis Module
#' @description A Shiny module for performing per-assay differential abundance
#'              analysis on lipidomics data using limma.
#' @name mod_lipid_de
NULL

#' @rdname mod_lipid_de
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h3 h4 h5 p hr br selectInput numericInput actionButton uiOutput verbatimTextOutput plotOutput tags downloadButton checkboxInput textAreaInput tabsetPanel tabPanel conditionalPanel helpText
#' @importFrom DT DTOutput
#' @importFrom shinyjqui jqui_resizable
mod_lipid_de_ui <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::fluidRow(
            # ================================================================
            # LEFT SIDEBAR (3/12)
            # ================================================================
            shiny::column(3
                , shiny::wellPanel(
                    shiny::h4("Differential Expression Analysis")

                    # Formula input
                    , shiny::textAreaInput(
                        ns("formula_string")
                        , "Model Formula:"
                        , value = "~ 0 + group"
                        , height = "60px"
                    )
                    , shiny::helpText("Formula fetched from S4 object @args if available")

                    , shiny::hr()

                    # Analysis Parameters
                    , shiny::h5("Analysis Parameters")

                    , shiny::numericInput(
                        ns("de_q_val_thresh")
                        , "Q-value Threshold:"
                        , value = 0.05
                        , min = 0.001
                        , max = 0.2
                        , step = 0.005
                    )

                    , shiny::numericInput(
                        ns("treat_lfc_cutoff")
                        , "Log Fold-Change Cutoff:"
                        , value = 0
                        , min = 0
                        , max = 2
                        , step = 0.1
                    )

                    , shiny::hr()

                    # Available contrasts display
                    , shiny::h5("Available Contrasts")
                    , shiny::verbatimTextOutput(ns("contrasts_display"))
                    , shiny::helpText("Contrasts defined in design matrix")

                    , shiny::hr()

                    # Load Filtered Session
                    , shiny::actionButton(
                        ns("load_filtered_session")
                        , "Load Filtered Session"
                        , class = "btn-info"
                        , width = "100%"
                        , icon = shiny::icon("upload")
                    )
                    , shiny::helpText("Load previously exported filtered data for DE analysis")

                    , shiny::br()
                    , shiny::br()

                    # Run DE Analysis
                    , shiny::actionButton(
                        ns("run_de_analysis")
                        , "Run DE Analysis"
                        , class = "btn-primary"
                        , width = "100%"
                        , icon = shiny::icon("play")
                    )

                    , shiny::br()
                    , shiny::br()

                    # Download Results
                    , shiny::downloadButton(
                        ns("download_de_results")
                        , "Download All Results"
                        , class = "btn-success"
                        , style = "width: 100%;"
                    )

                    , shiny::hr()

                    # Status display
                    , shiny::h5("Analysis Status")
                    , shiny::verbatimTextOutput(ns("de_status"))
                )
            )

            # ================================================================
            # RIGHT PANEL - TABSET (9/12)
            # ================================================================
            , shiny::column(9
                , shiny::tabsetPanel(
                    id = ns("de_results_tabs")

                    # ----------------------------------------------------------
                    # TAB 1: VOLCANO PLOT
                    # ----------------------------------------------------------
                    , shiny::tabPanel(
                        "Volcano Plot"
                        , shiny::br()
                        , shiny::fluidRow(
                            shiny::column(4
                                , shiny::selectInput(
                                    ns("volcano_contrast")
                                    , "Select Contrast:"
                                    , choices = NULL
                                    , width = "100%"
                                )
                            )
                            , shiny::column(4
                                , shiny::selectInput(
                                    ns("volcano_assay")
                                    , "Select Assay:"
                                    , choices = c("Combined")
                                    , width = "100%"
                                )
                            )
                            , shiny::column(4
                                , shiny::checkboxInput(
                                    ns("volcano_interactive")
                                    , "Interactive Plot (Glimma)"
                                    , value = TRUE
                                )
                            )
                        )

                        # Interactive volcano (Glimma)
                        , shiny::conditionalPanel(
                            condition = "input.volcano_interactive == true"
                            , ns = ns
                            , shiny::div(
                                id = ns("volcano_glimma_container")
                                , style = "height: 600px; background-color: white; border-radius: 8px; padding: 10px;"
                                , shiny::htmlOutput(ns("volcano_glimma"))
                            )
                        )

                        # Static volcano (ggplot2)
                        , shiny::conditionalPanel(
                            condition = "input.volcano_interactive == false"
                            , ns = ns
                            , shinyjqui::jqui_resizable(
                                shiny::plotOutput(ns("volcano_static"), height = "600px")
                            )
                        )
                    )

                    # ----------------------------------------------------------
                    # TAB 2: HEATMAP
                    # ----------------------------------------------------------
                    , shiny::tabPanel(
                        "Heatmap"
                        , shiny::br()

                        # Main controls (4 columns)
                        , shiny::fluidRow(
                            shiny::column(3
                                , shiny::selectInput(
                                    ns("heatmap_contrast")
                                    , "Select Contrast:"
                                    , choices = NULL
                                    , width = "100%"
                                )
                            )
                            , shiny::column(3
                                , shiny::numericInput(
                                    ns("heatmap_top_n")
                                    , "Top N Lipids:"
                                    , value = 50
                                    , min = 10
                                    , max = 500
                                    , step = 10
                                )
                            )
                            , shiny::column(3
                                , shiny::selectInput(
                                    ns("heatmap_clustering")
                                    , "Apply Clustering:"
                                    , choices = list(
                                        "Both rows & columns" = "both"
                                        , "Rows only" = "row"
                                        , "Columns only" = "column"
                                        , "None" = "none"
                                    )
                                    , selected = "both"
                                )
                            )
                            , shiny::column(3
                                , shiny::selectInput(
                                    ns("heatmap_scaling")
                                    , "Data Scaling:"
                                    , choices = list(
                                        "Row (lipid) scaling" = "row"
                                        , "Column (sample) scaling" = "column"
                                        , "Both" = "both"
                                        , "None" = "none"
                                    )
                                    , selected = "row"
                                )
                            )
                        )

                        # Advanced clustering options (conditional)
                        , shiny::conditionalPanel(
                            condition = "input.heatmap_clustering != 'none'"
                            , ns = ns
                            , shiny::wellPanel(
                                shiny::h5("Clustering Options")
                                , shiny::fluidRow(
                                    shiny::column(4
                                        , shiny::selectInput(
                                            ns("heatmap_cluster_method")
                                            , "Clustering Method:"
                                            , choices = list(
                                                "Ward (minimum variance)" = "ward.D2"
                                                , "Ward (original)" = "ward.D"
                                                , "Complete linkage" = "complete"
                                                , "Single linkage" = "single"
                                                , "Average linkage" = "average"
                                                , "McQuitty (WPGMA)" = "mcquitty"
                                            )
                                            , selected = "ward.D2"
                                        )
                                    )
                                    , shiny::column(4
                                        , shiny::selectInput(
                                            ns("heatmap_distance_method")
                                            , "Distance Metric:"
                                            , choices = list(
                                                "Euclidean" = "euclidean"
                                                , "Manhattan" = "manhattan"
                                                , "Pearson correlation" = "pearson"
                                                , "Spearman correlation" = "spearman"
                                                , "Maximum" = "maximum"
                                            )
                                            , selected = "euclidean"
                                        )
                                    )
                                    , shiny::column(4
                                        , shiny::selectInput(
                                            ns("heatmap_color_scheme")
                                            , "Color Scheme:"
                                            , choices = list(
                                                "Red-Blue" = "RdBu"
                                                , "Red-Yellow-Blue" = "RdYlBu"
                                                , "Blue-White-Red" = "coolwarm"
                                                , "Viridis" = "viridis"
                                                , "Plasma" = "plasma"
                                                , "Inferno" = "inferno"
                                            )
                                            , selected = "RdBu"
                                        )
                                    )
                                )
                                , shiny::fluidRow(
                                    shiny::column(4
                                        , shiny::selectInput(
                                            ns("heatmap_assay")
                                            , "Select Assay:"
                                            , choices = c("Combined")
                                            , width = "100%"
                                        )
                                    )
                                    , shiny::column(4
                                        , shiny::checkboxInput(
                                            ns("heatmap_show_labels")
                                            , "Show Lipid Labels"
                                            , value = FALSE
                                        )
                                    )
                                )
                            )
                        )

                        # Heatmap output
                        , shiny::div(
                            id = ns("heatmap_container")
                            , style = "height: 650px;"
                            , shinyjqui::jqui_resizable(
                                shiny::plotOutput(ns("heatmap_plot"), height = "600px")
                            )
                        )
                    )

                    # ----------------------------------------------------------
                    # TAB 3: DE RESULTS TABLE
                    # ----------------------------------------------------------
                    , shiny::tabPanel(
                        "DE Results Table"
                        , shiny::br()

                        # Controls
                        , shiny::fluidRow(
                            shiny::column(3
                                , shiny::selectInput(
                                    ns("table_contrast")
                                    , "Select Contrast:"
                                    , choices = NULL
                                    , width = "100%"
                                )
                            )
                            , shiny::column(3
                                , shiny::selectInput(
                                    ns("table_assay")
                                    , "Select Assay:"
                                    , choices = c("All")
                                    , width = "100%"
                                )
                            )
                            , shiny::column(3
                                , shiny::selectInput(
                                    ns("table_significance")
                                    , "Show:"
                                    , choices = list(
                                        "All Results" = "all"
                                        , "Significant Only" = "significant"
                                        , "Up-regulated" = "up"
                                        , "Down-regulated" = "down"
                                    )
                                    , selected = "significant"
                                )
                            )
                            , shiny::column(3
                                , shiny::numericInput(
                                    ns("table_max_rows")
                                    , "Max Rows:"
                                    , value = 1000
                                    , min = 100
                                    , max = 10000
                                    , step = 100
                                )
                            )
                        )

                        # Summary statistics
                        , shiny::fluidRow(
                            shiny::column(12
                                , shiny::wellPanel(
                                    shiny::h5("Summary Statistics")
                                    , shiny::verbatimTextOutput(ns("de_summary_stats"))
                                )
                            )
                        )

                        # Results table
                        , DT::DTOutput(ns("de_results_table"))
                    )
                )
            )
        )
    )
}

#' @rdname mod_lipid_de
#' @export
#' @importFrom shiny moduleServer reactiveValues reactive observeEvent req renderUI renderPlot renderText showNotification removeNotification updateSelectInput downloadHandler observe renderPrint
#' @importFrom DT renderDT datatable formatRound formatStyle styleEqual
#' @importFrom logger log_info log_error log_warn
mod_lipid_de_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # ================================================================
        # LOCAL REACTIVE VALUES
        # ================================================================
        de_data <- shiny::reactiveValues(
            de_results_list = NULL
            , contrasts_available = NULL
            , assays_available = NULL
            , analysis_complete = FALSE
            , current_s4_object = NULL
            , formula_from_s4 = NULL
        )

        # ================================================================
        # CONTRASTS DISPLAY
        # ================================================================
        output$contrasts_display <- shiny::renderPrint({
            contrasts_tbl <- workflow_data$contrasts_tbl

            if (is.null(contrasts_tbl) || nrow(contrasts_tbl) == 0) {
                cat("No contrasts defined.\nLoad a filtered session or define contrasts in the design tab.")
            } else {
                if ("friendly_names" %in% colnames(contrasts_tbl)) {
                    cat(paste(contrasts_tbl$friendly_names, collapse = "\n"))
                } else if ("contrasts" %in% colnames(contrasts_tbl)) {
                    cat(paste(contrasts_tbl$contrasts, collapse = "\n"))
                } else {
                    print(contrasts_tbl)
                }
            }
        })

        # ================================================================
        # STATUS DISPLAY
        # ================================================================
        output$de_status <- shiny::renderPrint({
            if (de_data$analysis_complete) {
                results <- de_data$de_results_list
                if (!is.null(results$significant_counts)) {
                    cat("Analysis Complete\n\n")
                    for (assay_name in names(results$significant_counts)) {
                        counts <- results$significant_counts[[assay_name]]
                        cat(sprintf("%s:\n  Up: %d | Down: %d | NS: %d\n"
                            , assay_name, counts$up, counts$down, counts$ns))
                    }
                } else {
                    cat("Analysis complete.")
                }
            } else if (!is.null(de_data$current_s4_object)) {
                cat("Data loaded. Ready to run analysis.")
            } else {
                cat("Waiting for data.\nClick 'Load Filtered Session' to begin.")
            }
        })

        # ================================================================
        # LOAD FILTERED SESSION
        # ================================================================
        shiny::observeEvent(input$load_filtered_session, {
            # [D66:START] ─────────────────────────
            d66_log <- function(...) message(sprintf("[D66] %s", paste0(...)))
            d66_log("─── ENTER load_filtered_session observer ───")
            .d66_start <- Sys.time()
            # [D66:END] ───────────────────────────

            logger::log_info("=== LOAD FILTERED SESSION BUTTON CLICKED ===")

            # [D66:START]
            d66_log("  STEP 1: Getting source_dir from experiment_paths")
            d66_log("    experiment_paths class: ", class(experiment_paths)[1])
            d66_log("    experiment_paths is NULL: ", is.null(experiment_paths))
            # [D66:END]

            # Look for the session file
            source_dir <- experiment_paths$source_dir

            # [D66:START]
            d66_log("  STEP 2: source_dir = ", if(is.null(source_dir)) "NULL" else source_dir)
            # [D66:END]

            if (is.null(source_dir) || !dir.exists(source_dir)) {
                # [D66:START]
                d66_log("  BRANCH: source_dir NULL or not exists, trying export_dir")
                # [D66:END]
                source_dir <- experiment_paths$export_dir
                # [D66:START]
                d66_log("  STEP 3: export_dir = ", if(is.null(source_dir)) "NULL" else source_dir)
                # [D66:END]
            }

            if (is.null(source_dir) || !dir.exists(source_dir)) {
                # [D66:START]
                d66_log("  ERROR: No valid source directory found")
                # [D66:END]
                shiny::showNotification(
                    "Could not find source directory for session data."
                    , type = "error"
                    , duration = 5
                )
                return()
            }

            session_file <- file.path(source_dir, "lipid_filtered_session_data_latest.rds")
            # [D66:START]
            d66_log("  STEP 4: session_file = ", session_file)
            d66_log("  STEP 4: file.exists = ", file.exists(session_file))
            # [D66:END]

            if (!file.exists(session_file)) {
                shiny::showNotification(
                    sprintf("Session file not found: %s", session_file)
                    , type = "error"
                    , duration = 5
                )
                return()
            }

            tryCatch({
                shiny::showNotification(
                    "Loading filtered session..."
                    , id = "loading_session"
                    , duration = NULL
                )

                # [D66:START]
                d66_log("  STEP 5: Reading RDS file...")
                # [D66:END]
                session_data <- readRDS(session_file)
                # [D66:START]
                d66_log("  STEP 5: RDS loaded successfully")
                d66_log("    names(session_data): ", paste(names(session_data), collapse = ", "))
                # [D66:END]

                logger::log_info(sprintf("   Loaded session from: %s", session_file))

                # [D66:START]
                d66_log("  STEP 6: Checking current_s4_object...")
                d66_log("    is.null(session_data$current_s4_object): ", is.null(session_data$current_s4_object))
                # [D66:END]

                # Restore S4 object
                if (!is.null(session_data$current_s4_object)) {
                    # [D66:START]
                    d66_log("  STEP 6a: S4 object class: ", class(session_data$current_s4_object)[1])
                    d66_log("  STEP 6b: Assigning to de_data$current_s4_object...")
                    # [D66:END]
                    de_data$current_s4_object <- session_data$current_s4_object

                    # [D66:START]
                    d66_log("  STEP 6c: Getting state_name...")
                    # [D66:END]
                    # Get state name with NULL-safe fallback
                    state_name <- session_data$r6_current_state_name
                    # [D66:START]
                    d66_log("    r6_current_state_name: ", if(is.null(state_name)) "NULL" else state_name)
                    # [D66:END]
                    if (is.null(state_name)) {
                        state_name <- "loaded_for_de"
                    }

                    # [D66:START]
                    d66_log("  STEP 6d: Checking workflow_data$state_manager...")
                    d66_log("    workflow_data is NULL: ", is.null(workflow_data))
                    d66_log("    workflow_data class: ", class(workflow_data)[1])
                    d66_log("    state_manager exists: ", !is.null(workflow_data$state_manager))
                    if (!is.null(workflow_data$state_manager)) {
                        d66_log("    state_manager class: ", class(workflow_data$state_manager)[1])
                    }
                    # [D66:END]

                    # [D66:START]
                    d66_log("  STEP 6e: Calling workflow_data$state_manager$saveState...")
                    # [D66:END]
                    workflow_data$state_manager$saveState(
                        state_name = state_name
                        , s4_data_object = session_data$current_s4_object
                        , config_object = list(loaded_from = session_file)
                        , description = "Loaded from filtered session for DE analysis"
                    )
                    # [D66:START]
                    d66_log("  STEP 6f: saveState completed")
                    # [D66:END]
                    logger::log_info("   Restored S4 object to state manager")
                }

                # [D66:START]
                d66_log("  STEP 7: Checking contrasts_tbl...")
                d66_log("    is.null(session_data$contrasts_tbl): ", is.null(session_data$contrasts_tbl))
                if (!is.null(session_data$contrasts_tbl)) {
                    d66_log("    nrow(contrasts_tbl): ", nrow(session_data$contrasts_tbl))
                    d66_log("    colnames: ", paste(colnames(session_data$contrasts_tbl), collapse = ", "))
                }
                # [D66:END]

                # Restore contrasts
                if (!is.null(session_data$contrasts_tbl) && nrow(session_data$contrasts_tbl) > 0) {
                    # [D66:START]
                    d66_log("  STEP 7a: Assigning contrasts to workflow_data...")
                    # [D66:END]
                    workflow_data$contrasts_tbl <- session_data$contrasts_tbl
                    de_data$contrasts_available <- session_data$contrasts_tbl

                    # [D66:START]
                    d66_log("  STEP 7b: Getting contrast_choices...")
                    # [D66:END]
                    # Get contrast choices
                    if ("friendly_names" %in% colnames(session_data$contrasts_tbl)) {
                        contrast_choices <- session_data$contrasts_tbl$friendly_names
                    } else {
                        contrast_choices <- session_data$contrasts_tbl$contrasts
                    }
                    # [D66:START]
                    d66_log("    contrast_choices: ", paste(contrast_choices, collapse = ", "))
                    # [D66:END]

                    # [D66:START]
                    d66_log("  STEP 7c: Updating volcano_contrast dropdown...")
                    # [D66:END]
                    # Update dropdowns
                    shiny::updateSelectInput(session, "volcano_contrast"
                        , choices = contrast_choices, selected = contrast_choices[1])
                    # [D66:START]
                    d66_log("  STEP 7d: Updating heatmap_contrast dropdown...")
                    # [D66:END]
                    shiny::updateSelectInput(session, "heatmap_contrast"
                        , choices = contrast_choices, selected = contrast_choices[1])
                    # [D66:START]
                    d66_log("  STEP 7e: Updating table_contrast dropdown...")
                    # [D66:END]
                    shiny::updateSelectInput(session, "table_contrast"
                        , choices = contrast_choices, selected = contrast_choices[1])

                    logger::log_info(sprintf("   Restored %d contrasts", nrow(session_data$contrasts_tbl)))
                }

                # [D66:START]
                d66_log("  STEP 8: Checking assay_names...")
                d66_log("    is.null(session_data$assay_names): ", is.null(session_data$assay_names))
                if (!is.null(session_data$assay_names)) {
                    d66_log("    assay_names: ", paste(session_data$assay_names, collapse = ", "))
                }
                # [D66:END]

                # Get assay names
                if (!is.null(session_data$assay_names)) {
                    de_data$assays_available <- session_data$assay_names
                    assay_choices <- c("Combined", session_data$assay_names)
                    # [D66:START]
                    d66_log("  STEP 8a: Updating volcano_assay dropdown...")
                    # [D66:END]
                    shiny::updateSelectInput(session, "volcano_assay", choices = assay_choices)
                    # [D66:START]
                    d66_log("  STEP 8b: Updating heatmap_assay dropdown...")
                    # [D66:END]
                    shiny::updateSelectInput(session, "heatmap_assay", choices = assay_choices)
                    # [D66:START]
                    d66_log("  STEP 8c: Updating table_assay dropdown...")
                    # [D66:END]
                    shiny::updateSelectInput(session, "table_assay", choices = c("All", session_data$assay_names))
                }

                # [D66:START]
                d66_log("  STEP 9: Attempting to extract formula from S4 args...")
                # [D66:END]
                # Try to get formula from S4 args (safely)
                tryCatch({
                    # [D66:START]
                    d66_log("  STEP 9a: Checking S4 @args slot...")
                    d66_log("    hasSlot 'args': ", methods::hasMethod("@", "LipidomicsAssayData"))
                    # [D66:END]
                    s4_args <- session_data$current_s4_object@args
                    # [D66:START]
                    d66_log("  STEP 9b: s4_args is NULL: ", is.null(s4_args))
                    if (!is.null(s4_args)) {
                        d66_log("    s4_args class: ", class(s4_args)[1])
                        d66_log("    names(s4_args): ", paste(names(s4_args), collapse = ", "))
                    }
                    # [D66:END]
                    if (!is.null(s4_args) && !is.null(s4_args$deAnalysisParameters)) {
                        formula_val <- s4_args$deAnalysisParameters$formula_string
                        if (!is.null(formula_val) && nzchar(formula_val)) {
                            de_data$formula_from_s4 <- formula_val
                            shiny::updateTextAreaInput(session, "formula_string"
                                , value = formula_val)
                        }
                    }
                }, error = function(e) {
                    # [D66:START]
                    d66_log("  STEP 9 ERROR: ", e$message)
                    # [D66:END]
                    logger::log_warn(sprintf("Could not extract formula from S4: %s", e$message))
                })

                # [D66:START]
                d66_log("  STEP 10: Removing notification and showing success...")
                # [D66:END]
                shiny::removeNotification("loading_session")
                shiny::showNotification("Session loaded successfully!", type = "message", duration = 3)

                # [D66:START]
                d66_log("─── EXIT load_filtered_session (",
                    round(difftime(Sys.time(), .d66_start, units = "secs"), 3), "s) ─── SUCCESS")
                # [D66:END]

            }, error = function(e) {
                # [D66:START]
                d66_log("  FATAL ERROR: ", e$message)
                d66_log("  ERROR CALL: ", deparse(e$call))
                d66_log("─── EXIT load_filtered_session ─── FAILED")
                # [D66:END]
                logger::log_error(sprintf("   Error loading session: %s", e$message))
                shiny::removeNotification("loading_session")
                shiny::showNotification(
                    sprintf("Error loading session: %s", e$message)
                    , type = "error"
                    , duration = 10
                )
            })
        })

        # ================================================================
        # RUN DE ANALYSIS
        # ================================================================
        shiny::observeEvent(input$run_de_analysis, {
            logger::log_info("=== RUN DE ANALYSIS BUTTON CLICKED ===")

            # Validate inputs
            current_s4 <- de_data$current_s4_object
            if (is.null(current_s4)) {
                current_s4 <- workflow_data$state_manager$getState()
            }

            if (is.null(current_s4) || !inherits(current_s4, "LipidomicsAssayData")) {
                shiny::showNotification(
                    "No lipidomics data loaded. Please load a filtered session first."
                    , type = "error"
                    , duration = 5
                )
                return()
            }

            contrasts_tbl <- workflow_data$contrasts_tbl
            if (is.null(contrasts_tbl) || nrow(contrasts_tbl) == 0) {
                shiny::showNotification(
                    "No contrasts defined. Please define contrasts in the design tab."
                    , type = "error"
                    , duration = 5
                )
                return()
            }

            shiny::showNotification(
                "Running differential expression analysis..."
                , id = "de_running"
                , duration = NULL
            )

            tryCatch({
                # Run DE analysis
                results <- runLipidsDE(
                    theObject = current_s4
                    , contrasts_tbl = contrasts_tbl
                    , formula_string = input$formula_string
                    , de_q_val_thresh = input$de_q_val_thresh
                    , treat_lfc_cutoff = input$treat_lfc_cutoff
                    , eBayes_trend = TRUE
                    , eBayes_robust = TRUE
                )

                de_data$de_results_list <- results
                de_data$analysis_complete <- TRUE

                # Update tab status - must replace entire list to trigger reactivity
                updated_status <- workflow_data$tab_status
                updated_status$differential_analysis <- "complete"
                workflow_data$tab_status <- updated_status

                shiny::removeNotification("de_running")
                shiny::showNotification(
                    "Differential expression analysis complete!"
                    , type = "message"
                    , duration = 5
                )

                logger::log_info("   DE analysis completed successfully")

                # =============================================================
                # UPDATE UI DROPDOWNS WITH RESULTS
                # =============================================================
                if (!is.null(results$de_lipids_long)) {
                    contrast_choices <- unique(results$de_lipids_long$friendly_name)
                    if (length(contrast_choices) == 0) {
                        contrast_choices <- unique(results$de_lipids_long$comparison)
                    }

                    shiny::updateSelectInput(session, "volcano_contrast"
                        , choices = contrast_choices, selected = contrast_choices[1])
                    shiny::updateSelectInput(session, "heatmap_contrast"
                        , choices = contrast_choices, selected = contrast_choices[1])
                    shiny::updateSelectInput(session, "table_contrast"
                        , choices = contrast_choices, selected = contrast_choices[1])

                    # Update assay dropdowns
                    assay_choices <- c("Combined", unique(results$de_lipids_long$assay))
                    shiny::updateSelectInput(session, "volcano_assay"
                        , choices = assay_choices, selected = "Combined")
                    shiny::updateSelectInput(session, "heatmap_assay"
                        , choices = assay_choices, selected = "Combined")

                    table_assay_choices <- c("All", unique(results$de_lipids_long$assay))
                    shiny::updateSelectInput(session, "table_assay"
                        , choices = table_assay_choices, selected = "All")

                    logger::log_info(sprintf("   Updated dropdowns: %d contrasts, %d assays"
                        , length(contrast_choices), length(assay_choices) - 1))
                }

                # =============================================================
                # WRITE DE RESULTS TO DISK
                # =============================================================
                logger::log_info("   Writing DE results to disk...")

                tryCatch({
                    de_output_dir <- experiment_paths$de_output_dir
                    publication_graphs_dir <- experiment_paths$publication_graphs_dir

                    logger::log_info(sprintf("   de_output_dir = %s", de_output_dir))
                    logger::log_info(sprintf("   publication_graphs_dir = %s", publication_graphs_dir))

                    if (is.null(de_output_dir) || is.null(publication_graphs_dir)) {
                        logger::log_warn("   Output directories not configured, skipping file output")
                    } else {
                        success <- outputMetabDeResultsAllContrasts(
                            de_results_list = results
                            , de_output_dir = de_output_dir
                            , publication_graphs_dir = publication_graphs_dir
                            , de_q_val_thresh = input$de_q_val_thresh
                            , lfc_threshold = input$treat_lfc_cutoff
                            , heatmap_top_n = 50
                            , heatmap_clustering = "both"
                            , heatmap_color_scheme = "RdBu"
                        )

                        if (success) {
                            logger::log_info("   All DE results written to disk successfully")
                            shiny::showNotification(
                                "DE results saved to disk (tables, volcano plots, heatmaps)"
                                , type = "message"
                                , duration = 5
                            )
                        }
                    }

                }, error = function(e) {
                    logger::log_warn(paste("   Could not write DE results to disk:", e$message))
                    shiny::showNotification(
                        paste("Warning: Could not save results to disk:", e$message)
                        , type = "warning"
                        , duration = 8
                    )
                })

            }, error = function(e) {
                logger::log_error(sprintf("   DE analysis error: %s", e$message))
                shiny::removeNotification("de_running")
                shiny::showNotification(
                    sprintf("Analysis error: %s", e$message)
                    , type = "error"
                    , duration = 10
                )
            })
        })

        # ================================================================
        # VOLCANO PLOT - GLIMMA
        # ================================================================
        output$volcano_glimma <- shiny::renderUI({
            shiny::req(de_data$de_results_list, input$volcano_contrast)

            # Show informative message for Combined view
            if (is.null(input$volcano_assay) || input$volcano_assay == "Combined") {
                return(shiny::div(
                    class = "alert alert-info"
                    , style = "margin: 20px;"
                    , shiny::icon("info-circle")
                    , shiny::tags$strong(" Combined View: ")
                    , "Interactive Glimma plots require a single assay selection. "
                    , "Please select a specific assay (e.g., LCMS_Pos) or uncheck 'Interactive Plot' to use the static volcano plot for combined view."
                ))
            }

            tryCatch({
                widget <- generateMetabVolcanoPlotGlimma(
                    de_results_list = de_data$de_results_list
                    , selected_contrast = input$volcano_contrast
                    , selected_assay = input$volcano_assay
                    , de_q_val_thresh = input$de_q_val_thresh
                )

                if (is.null(widget)) {
                    return(shiny::div(
                        class = "alert alert-warning"
                        , "Could not generate Glimma plot. Try the static plot option."
                    ))
                }

                widget

            }, error = function(e) {
                logger::log_error(sprintf("Glimma error: %s", e$message))
                shiny::div(
                    class = "alert alert-danger"
                    , sprintf("Error generating plot: %s", e$message)
                )
            })
        })

        # ================================================================
        # VOLCANO PLOT - STATIC
        # ================================================================
        output$volcano_static <- shiny::renderPlot({
            shiny::req(de_data$de_results_list, input$volcano_contrast)

            generateMetabVolcanoStatic(
                de_results_list = de_data$de_results_list
                , selected_contrast = input$volcano_contrast
                , selected_assay = input$volcano_assay
                , de_q_val_thresh = input$de_q_val_thresh
                , lfc_threshold = input$treat_lfc_cutoff
                , show_labels = TRUE
                , n_labels = 15
            )
        })

        # ================================================================
        # HEATMAP
        # ================================================================
        output$heatmap_plot <- shiny::renderPlot({
            shiny::req(de_data$de_results_list, input$heatmap_contrast)

            hm <- generateMetabDEHeatmap(
                de_results_list = de_data$de_results_list
                , selected_contrast = input$heatmap_contrast
                , selected_assay = input$heatmap_assay
                , top_n = input$heatmap_top_n
                , clustering_method = input$heatmap_cluster_method
                , distance_method = input$heatmap_distance_method
                , cluster_rows = input$heatmap_clustering %in% c("both", "row")
                , cluster_cols = input$heatmap_clustering %in% c("both", "column")
                , scale_data = input$heatmap_scaling
                , color_scheme = input$heatmap_color_scheme
                , show_lipid_names = input$heatmap_show_labels
                , de_q_val_thresh = input$de_q_val_thresh
            )

            if (!is.null(hm)) {
                ComplexHeatmap::draw(hm)
            }
        })

        # ================================================================
        # DE RESULTS TABLE
        # ================================================================
        output$de_summary_stats <- shiny::renderPrint({
            shiny::req(de_data$de_results_list)

            results <- de_data$de_results_list$de_lipids_long
            if (is.null(results) || nrow(results) == 0) {
                cat("No results available.")
                return()
            }

            # Filter by contrast
            if (!is.null(input$table_contrast)) {
                results <- results[results$comparison == input$table_contrast |
                                   results$friendly_name == input$table_contrast, ]
            }

            # Filter by assay
            if (!is.null(input$table_assay) && input$table_assay != "All") {
                results <- results[results$assay == input$table_assay, ]
            }

            total <- nrow(results)
            sig <- sum(results$significant != "NS", na.rm = TRUE)
            up <- sum(results$significant == "Up", na.rm = TRUE)
            down <- sum(results$significant == "Down", na.rm = TRUE)

            cat(sprintf("Total lipids: %d\n", total))
            cat(sprintf("Significant (Q < %.3f): %d (%.1f%%)\n"
                , input$de_q_val_thresh, sig, 100 * sig / max(total, 1)))
            cat(sprintf("  Up-regulated: %d\n", up))
            cat(sprintf("  Down-regulated: %d\n", down))
        })

        output$de_results_table <- DT::renderDT({
            shiny::req(de_data$de_results_list)

            results <- de_data$de_results_list$de_lipids_long
            if (is.null(results) || nrow(results) == 0) {
                return(NULL)
            }

            # Filter by contrast
            if (!is.null(input$table_contrast)) {
                results <- results[results$comparison == input$table_contrast |
                                   results$friendly_name == input$table_contrast, ]
            }

            # Filter by assay
            if (!is.null(input$table_assay) && input$table_assay != "All") {
                results <- results[results$assay == input$table_assay, ]
            }

            # Filter by significance
            if (input$table_significance == "significant") {
                results <- results[results$fdr_qvalue < input$de_q_val_thresh, ]
            } else if (input$table_significance == "up") {
                results <- results[results$fdr_qvalue < input$de_q_val_thresh &
                                   results$logFC > input$treat_lfc_cutoff, ]
            } else if (input$table_significance == "down") {
                results <- results[results$fdr_qvalue < input$de_q_val_thresh &
                                   results$logFC < -input$treat_lfc_cutoff, ]
            }

            # Limit rows
            if (nrow(results) > input$table_max_rows) {
                results <- results[1:input$table_max_rows, ]
            }

            # Select display columns
            display_cols <- c("lipid_id", "lipid_name", "assay", "logFC"
                , "raw_pvalue", "fdr_qvalue", "significant")
            display_cols <- intersect(display_cols, colnames(results))
            results <- results[, display_cols, drop = FALSE]

            DT::datatable(
                results
                , options = list(
                    pageLength = 25
                    , scrollX = TRUE
                    , dom = 'Bfrtip'
                    , buttons = c('copy', 'csv', 'excel')
                )
                , extensions = 'Buttons'
                , rownames = FALSE
                , filter = "top"
            ) |>
                DT::formatRound(columns = c("logFC"), digits = 3) |>
                DT::formatRound(columns = c("raw_pvalue", "fdr_qvalue"), digits = 6) |>
                DT::formatStyle(
                    "significant"
                    , backgroundColor = DT::styleEqual(
                        c("Up", "Down", "NS")
                        , c("#ffcccc", "#cce5ff", "white")
                    )
                )
        })

        # ================================================================
        # DOWNLOAD HANDLER
        # ================================================================
        output$download_de_results <- shiny::downloadHandler(
            filename = function() {
                paste0("lipidomics_de_results_", Sys.Date(), ".csv")
            }
            , content = function(file) {
                results <- de_data$de_results_list$de_lipids_long
                shiny::req(results)
                write.csv(results, file, row.names = FALSE)
            }
        )
    })
}
