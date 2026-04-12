#' @rdname differentialAbundanceAppletModule
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h4 h5 textAreaInput helpText hr numericInput actionButton icon verbatimTextOutput br downloadButton
mod_prot_da_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::wellPanel(
      shiny::fluidRow(
        # Sidebar Settings
        shiny::column(3, mod_prot_da_sidebar_ui(ns)),

        # Main Results Panel
        shiny::column(
          9,
          shiny::uiOutput(ns("heatmap_manual_save_warning")),
          shiny::tabsetPanel(
            id = ns("da_results_tabs"),
            mod_prot_da_volcano_tab_ui(ns),
            mod_prot_da_heatmap_tab_ui(ns),
            mod_prot_da_table_tab_ui(ns)
          )
        )
      )
    )
  )
}

mod_prot_da_sidebar_ui <- function(ns) {
  shiny::wellPanel(
    shiny::h4("DA Analysis Settings"),

    # Formula input (will be populated from S4 @args)
    shiny::textAreaInput(
      ns("formula_string"),
      "Model Formula:",
      value = "~ 0 + group",
      height = "60px"
    ),
    shiny::helpText("Formula fetched from S4 object @args"),
    shiny::hr(),

    # Analysis parameters
    shiny::h5("Analysis Parameters"),
    shiny::numericInput(
      ns("da_q_val_thresh"),
      "Q-value Threshold:",
      value = 0.05,
      min = 0.001,
      max = 0.2,
      step = 0.005
    ),
    shiny::numericInput(
      ns("treat_lfc_cutoff"),
      "Log Fold-Change Cutoff:",
      value = 0,
      min = 0,
      max = 2,
      step = 0.1
    ),
    shiny::hr(),

    # Available contrasts display
    shiny::h5("Available Contrasts"),
    shiny::verbatimTextOutput(ns("contrasts_display")),
    shiny::helpText("Contrasts defined in design matrix"),
    shiny::br(),

    # Load filtered session button
    shiny::actionButton(
      ns("load_filtered_session"),
      "Load Filtered Session",
      class = "btn-info",
      width = "100%",
      icon = shiny::icon("upload")
    ),
    shiny::helpText("Load previously exported filtered data for DA analysis"),
    shiny::br(),

    # Main action button
    shiny::actionButton(
      ns("run_da_analysis"),
      "Run DA Analysis",
      class = "btn-primary",
      width = "100%",
      icon = shiny::icon("play")
    ),
    shiny::br(),
    shiny::br(),

    # Download button
    shiny::downloadButton(
      ns("download_da_results"),
      "Download All Results",
      class = "btn-success",
      width = "100%"
    ),
    shiny::br(),
    shiny::br(),

    # Status display
    shiny::h5("Analysis Status"),
    shiny::verbatimTextOutput(ns("da_status"))
  )
}

mod_prot_da_volcano_tab_ui <- function(ns) {
  shiny::tabPanel(
    "Volcano Plot",
    icon = shiny::icon("chart-area"),
    shiny::br(),

    # Contrast selector for volcano plot
    shiny::fluidRow(
      shiny::column(
        6,
        shiny::selectInput(
          ns("volcano_contrast"),
          "Select Contrast:",
          choices = NULL,
          width = "100%"
        )
      ),
      shiny::column(
        6,
        shiny::checkboxInput(
          ns("volcano_interactive"),
          "Interactive Plot (Glimma)",
          value = TRUE
        )
      )
    ),
    shiny::hr(),

    # Volcano plot output (will switch between static and interactive)
    shiny::conditionalPanel(
      condition = "input.volcano_interactive == true",
      ns = ns,
      shiny::div(
        id = ns("volcano_glimma_container"),
        style = "height: 600px; background-color: white; border-radius: 8px; padding: 10px;",
        shiny::htmlOutput(ns("volcano_glimma"))
      )
    ),
    shiny::conditionalPanel(
      condition = "input.volcano_interactive == false",
      ns = ns,
      shinyjqui::jqui_resizable(
        shiny::plotOutput(ns("volcano_static"), height = "600px")
      )
    )
  )
}

mod_prot_da_heatmap_tab_ui <- function(ns) {
  shiny::tabPanel(
    "Heatmap",
    icon = shiny::icon("th"),
    shiny::br(),

    # Contrast selector and heatmap options
    shiny::fluidRow(
      shiny::column(
        3,
        shiny::selectInput(
          ns("heatmap_contrast"),
          "Select Contrast:",
          choices = NULL,
          width = "100%"
        )
      ),
      shiny::column(
        3,
        shiny::numericInput(
          ns("heatmap_top_n"),
          "Top N Genes:",
          value = 50,
          min = 10,
          max = 500,
          step = 10
        )
      ),
      shiny::column(
        3,
        shiny::selectInput(
          ns("heatmap_clustering"),
          "Apply Clustering:",
          choices = list(
            "Both rows & columns" = "both",
            "Rows only" = "row",
            "Columns only" = "column",
            "None" = "none"
          ),
          selected = "both"
        )
      ),
      shiny::column(
        3,
        shiny::selectInput(
          ns("heatmap_scaling"),
          "Data Scaling:",
          choices = list(
            "Row (gene) scaling" = "row",
            "Column (sample) scaling" = "column",
            "Both" = "both",
            "None" = "none"
          ),
          selected = "row"
        )
      )
    ),

    # Advanced clustering controls (collapsible)
    shiny::conditionalPanel(
      condition = "input.heatmap_clustering != 'none'",
      ns = ns,
      shiny::wellPanel(
        shiny::h5("Advanced Clustering Options"),
        shiny::fluidRow(
          shiny::column(
            3,
            shiny::selectInput(
              ns("heatmap_cluster_method"),
              "Clustering Method:",
              choices = list(
                "Ward (minimum variance)" = "ward.D2",
                "Ward (original)" = "ward.D",
                "Complete linkage" = "complete",
                "Single linkage" = "single",
                "Average linkage" = "average",
                "McQuitty (WPGMA)" = "mcquitty",
                "Median (WPGMC)" = "median",
                "Centroid (UPGMC)" = "centroid"
              ),
              selected = "ward.D2"
            )
          ),
          shiny::column(
            3,
            shiny::selectInput(
              ns("heatmap_distance_method"),
              "Distance Metric:",
              choices = list(
                "Euclidean" = "euclidean",
                "Manhattan" = "manhattan",
                "Pearson correlation" = "pearson",
                "Spearman correlation" = "spearman",
                "Maximum" = "maximum",
                "Canberra" = "canberra",
                "Binary" = "binary",
                "Minkowski" = "minkowski"
              ),
              selected = "euclidean"
            )
          ),
          shiny::column(
            3,
            shiny::selectInput(
              ns("heatmap_tree_cut_method"),
              "Tree Cutting:",
              choices = list(
                "Number of clusters" = "k_clusters",
                "Height cutoff" = "height_cutoff",
                "Dynamic tree cutting" = "dynamic",
                "No cutting" = "none"
              ),
              selected = "k_clusters"
            )
          ),
          shiny::column(
            3,
            # Dynamic UI for tree cutting parameters
            shiny::conditionalPanel(
              condition = "input.heatmap_tree_cut_method == 'k_clusters'",
              ns = ns,
              shiny::numericInput(
                ns("heatmap_n_clusters"),
                "Number of Clusters:",
                value = 4,
                min = 2,
                max = 20,
                step = 1
              )
            ),
            shiny::conditionalPanel(
              condition = "input.heatmap_tree_cut_method == 'height_cutoff'",
              ns = ns,
              shiny::numericInput(
                ns("heatmap_cut_height"),
                "Cut Height:",
                value = 0.5,
                min = 0.1,
                max = 2.0,
                step = 0.1
              )
            ),
            shiny::conditionalPanel(
              condition = "input.heatmap_tree_cut_method == 'dynamic'",
              ns = ns,
              shiny::numericInput(
                ns("heatmap_min_cluster_size"),
                "Min Cluster Size:",
                value = 3,
                min = 2,
                max = 20,
                step = 1
              )
            )
          )
        ),

        # Additional heatmap display options
        shiny::fluidRow(
          shiny::column(
            4,
            shiny::checkboxInput(
              ns("heatmap_show_dendro"),
              "Show Dendrogram",
              value = TRUE
            )
          ),
          shiny::column(
            4,
            shiny::checkboxInput(
              ns("heatmap_show_labels"),
              "Show Gene Labels",
              value = FALSE
            )
          ),
          shiny::column(
            4,
            shiny::selectInput(
              ns("heatmap_color_scheme"),
              "Color Scheme:",
              choices = list(
                "Red-Blue" = "RdBu",
                "Red-Yellow-Blue" = "RdYlBu",
                "Blue-White-Red" = "coolwarm",
                "Viridis" = "viridis",
                "Plasma" = "plasma",
                "Inferno" = "inferno"
              ),
              selected = "RdBu"
            )
          )
        )
      )
    ),
    shiny::hr(),

    # Interactive heatmap output
    shiny::div(
      id = ns("heatmap_container"),
      style = "height: 650px;",
      shinyjqui::jqui_resizable(
        shiny::plotOutput(ns("heatmap_plot"), height = "600px")
      )
    ),

    # Save Button
    shiny::div(
      style = "margin-top: 10px; margin-bottom: 20px;",
      shiny::actionButton(
        ns("save_heatmap"),
        "Save Heatmap",
        icon = shiny::icon("download"),
        class = "btn-primary"
      )
    ),

    # Cluster summary (if tree cutting is applied)
    shiny::conditionalPanel(
      condition = "input.heatmap_tree_cut_method != 'none'",
      ns = ns,
      shiny::wellPanel(
        shiny::h5("Cluster Information"),
        shiny::verbatimTextOutput(ns("cluster_summary"))
      )
    )
  )
}

mod_prot_da_table_tab_ui <- function(ns) {
  shiny::tabPanel(
    "DA Results Table",
    icon = shiny::icon("table"),
    shiny::br(),

    # Contrast selector and table options
    shiny::fluidRow(
      shiny::column(
        4,
        shiny::selectInput(
          ns("table_contrast"),
          "Select Contrast:",
          choices = NULL,
          width = "100%"
        )
      ),
      shiny::column(
        4,
        shiny::selectInput(
          ns("table_significance"),
          "Show:",
          choices = list(
            "All Results" = "all",
            "Significant Only" = "significant",
            "Up-regulated" = "up",
            "Down-regulated" = "down"
          ),
          selected = "significant"
        )
      ),
      shiny::column(
        4,
        shiny::numericInput(
          ns("table_max_rows"),
          "Max Rows:",
          value = 1000,
          min = 100,
          max = 10000,
          step = 100
        )
      )
    ),
    shiny::hr(),

    # Summary statistics
    shiny::fluidRow(
      shiny::column(
        12,
        shiny::wellPanel(
          shiny::h5("Summary Statistics"),
          shiny::verbatimTextOutput(ns("da_summary_stats"))
        )
      )
    ),

    # Interactive results table
    DT::DTOutput(ns("da_results_table"))
  )
}

