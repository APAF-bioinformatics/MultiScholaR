#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinydashboard
#' @import shinyjs
#' @noRd
app_ui <- function(request) {
  dashboardPage(
    dashboardHeader(
      title = "MultiScholaR",
      # Dark mode toggle in header
      tags$li(
        class = "dropdown",
        tags$a(
          id = "dark_mode_toggle",
          href = "#",
          onclick = "toggleDarkMode(); return false;",
          style = "padding: 10px 15px;",
          icon("moon"),
          title = "Toggle Dark Mode"
        )
      )
    ),
    
    dashboardSidebar(
      shiny::div(
        style = "text-align: center; padding: 10px;",
        shiny::tags$img(src = "www/MultiScholaR.png", style = "width: 80%; height: auto;")
      ),
      sidebarMenu(
        id = "main_menu",
        menuItem("Home", tabName = "home", icon = shiny::icon("home")),
        
        # Dynamic menu items will be added based on omics selection
        uiOutput("dynamic_menu")
      )
    ),
    
    dashboardBody(
      useShinyjs(),
      
      # CSS for better styling - use shiny:: prefix
      shiny::tags$head(
        shiny::tags$style(shiny::HTML("
        /* ========== DARK MODE STYLES ========== */
        body.dark-mode {
          background-color: #1a1a2e !important;
        }
        body.dark-mode .wrapper,
        body.dark-mode .content-wrapper,
        body.dark-mode .right-side {
          background-color: #1a1a2e !important;
        }
        body.dark-mode .main-header .navbar,
        body.dark-mode .main-header .logo {
          background-color: #16213e !important;
        }
        body.dark-mode .main-sidebar,
        body.dark-mode .left-side {
          background-color: #0f3460 !important;
        }
        body.dark-mode .sidebar-menu > li > a {
          color: #e4e4e4 !important;
        }
        body.dark-mode .sidebar-menu > li.active > a,
        body.dark-mode .sidebar-menu > li:hover > a {
          background-color: #1a1a2e !important;
          color: #fff !important;
        }
        body.dark-mode .box,
        body.dark-mode .well,
        body.dark-mode .nav-tabs-custom,
        body.dark-mode .omic-selection-box {
          background-color: #16213e !important;
          border-color: #0f3460 !important;
          color: #e4e4e4 !important;
        }
        body.dark-mode .omic-selection-box:hover {
          box-shadow: 0 4px 8px rgba(0,0,0,0.4) !important;
        }
        body.dark-mode .omic-selection-box.selected {
          background-color: #1a4a6e !important;
          border-color: #3498db !important;
        }
        body.dark-mode h1, body.dark-mode h2, body.dark-mode h3,
        body.dark-mode h4, body.dark-mode h5, body.dark-mode p,
        body.dark-mode label, body.dark-mode .box-title {
          color: #e4e4e4 !important;
        }
        body.dark-mode .form-control {
          background-color: #1a1a2e !important;
          border-color: #0f3460 !important;
          color: #e4e4e4 !important;
        }
        body.dark-mode .nav-tabs > li > a {
          color: #e4e4e4 !important;
        }
        body.dark-mode .nav-tabs > li.active > a {
          background-color: #16213e !important;
          border-color: #0f3460 !important;
        }
        body.dark-mode .dataTables_wrapper {
          color: #e4e4e4 !important;
        }
        body.dark-mode table.dataTable {
          background-color: #16213e !important;
          color: #e4e4e4 !important;
        }
        body.dark-mode table.dataTable thead th {
          background-color: #0f3460 !important;
          color: #e4e4e4 !important;
        }
        body.dark-mode table.dataTable tbody tr {
          background-color: #16213e !important;
        }
        body.dark-mode table.dataTable tbody tr:hover {
          background-color: #1a4a6e !important;
        }
        /* Modals */
        body.dark-mode .modal-content {
          background-color: #16213e !important;
          border-color: #0f3460 !important;
          color: #e4e4e4 !important;
        }
        body.dark-mode .modal-header {
          background-color: #0f3460 !important;
          border-bottom-color: #1a4a6e !important;
          color: #e4e4e4 !important;
        }
        body.dark-mode .modal-header .close {
          color: #e4e4e4 !important;
          opacity: 0.8;
        }
        body.dark-mode .modal-title {
          color: #e4e4e4 !important;
        }
        body.dark-mode .modal-body {
          background-color: #16213e !important;
          color: #e4e4e4 !important;
        }
        body.dark-mode .modal-footer {
          background-color: #16213e !important;
          border-top-color: #0f3460 !important;
        }
        body.dark-mode .modal-body .form-control,
        body.dark-mode .modal-body input,
        body.dark-mode .modal-body select,
        body.dark-mode .modal-body textarea {
          background-color: #1a1a2e !important;
          border-color: #0f3460 !important;
          color: #e4e4e4 !important;
        }
        body.dark-mode .modal-body .form-control::placeholder {
          color: #888 !important;
        }
        body.dark-mode .modal-body label {
          color: #e4e4e4 !important;
        }
        body.dark-mode .modal-body .help-block,
        body.dark-mode .modal-body .shiny-input-container .control-label {
          color: #b0b0b0 !important;
        }
        /* Selectize/select2 dropdowns in modals */
        body.dark-mode .selectize-input,
        body.dark-mode .selectize-dropdown {
          background-color: #1a1a2e !important;
          border-color: #0f3460 !important;
          color: #e4e4e4 !important;
        }
        body.dark-mode .selectize-dropdown-content .option {
          color: #e4e4e4 !important;
        }
        body.dark-mode .selectize-dropdown-content .option:hover,
        body.dark-mode .selectize-dropdown-content .option.active {
          background-color: #1a4a6e !important;
        }
        /* Smooth transition */
        body, .wrapper, .content-wrapper, .main-header, .main-sidebar,
        .box, .well, .form-control, .nav-tabs-custom, .omic-selection-box {
          transition: background-color 0.3s ease, color 0.3s ease, border-color 0.3s ease;
        }
        /* Toggle icon styling */
        #dark_mode_toggle {
          color: #fff;
          font-size: 18px;
        }
        body.dark-mode #dark_mode_toggle .fa-moon::before {
          content: '\\f185'; /* sun icon */
        }

        /* ========== LIGHT MODE (DEFAULT) ========== */
        /* Keep main content at normal position */
        .content-wrapper, .right-side {
          background-color: #f4f4f4;
        }
        .omic-selection-box {
          border: 2px solid #3c8dbc;
          border-radius: 10px;
          padding: 20px;
          margin: 10px;
          background-color: white;
          cursor: pointer;
          transition: all 0.3s;
        }
        .omic-selection-box:hover {
          transform: scale(1.05);
          box-shadow: 0 4px 8px rgba(0,0,0,0.1);
        }
        .omic-selection-box.selected {
          background-color: #d1ecf1;
          border-color: #bee5eb;
        }
        #log-terminal-wrapper {
          position: fixed;
          bottom: 10px;
          left: 10px;
          width: 230px;
          height: 250px;
          z-index: 1050; /* Above most other elements */
        }
        #log-terminal-panel {
          background-color: #2b3e50;
          border-radius: 5px;
          box-shadow: 0 4px 8px rgba(0,0,0,0.2);
          color: #ecf0f1;
          height: 100%;
          display: flex;
          flex-direction: column;
        }
        #log-terminal-header {
          padding: 5px 10px;
          background-color: #34495e;
          cursor: pointer;
          border-bottom: 1px solid #2c3e50;
          border-top-left-radius: 5px;
          border-top-right-radius: 5px;
        }
        #log-terminal-header h5 {
          margin: 0;
          color: #ecf0f1;
          font-weight: bold;
        }
        #log-terminal-content {
          padding: 10px;
          overflow-y: auto;
          flex-grow: 1;
          font-family: 'Courier New', Courier, monospace;
          font-size: 0.8em;
          white-space: pre-wrap;
          word-break: break-all;
          display: block;
        }
      "))
      ),
      
      shinydashboard::tabItems(
        # Home tab with omics selection
        shinydashboard::tabItem(
          tabName = "home",
          fluidRow(
            column(12,
                   h2("Welcome to MultiScholaR Interactive Workflow"),
                   p("Select the type(s) of omics data you want to analyze:"),
                   br()
            )
          ),
          
          fluidRow(
            column(4,
                   shiny::div(
                     id = "proteomics_box",
                     class = "omic-selection-box",
                     onclick = "toggleOmicSelection('proteomics')",
                     h3(shiny::icon("dna"), "Proteomics"),
                     p("Analyze protein abundance data from DIA-NN or similar tools")
                   )
            ),
            column(4,
                   shiny::div(
                     id = "metabolomics_box",
                     class = "omic-selection-box",
                     onclick = "toggleOmicSelection('metabolomics')",
                     h3(shiny::icon("flask"), "Metabolomics"),
                     p("Process and analyze metabolite profiling data")
                   )
            ),
            column(4,
                   shiny::div(
                     id = "transcriptomics_box",
                     class = "omic-selection-box",
                     onclick = "toggleOmicSelection('transcriptomics')",
                     h3(shiny::icon("chart-line"), "Transcriptomics"),
                     p("Analyze gene expression data from RNA-seq")
                   )
            )
          ),
          
          fluidRow(
            column(4,
                   shiny::div(
                     id = "lipidomics_box",
                     class = "omic-selection-box",
                     onclick = "toggleOmicSelection('lipidomics')",
                     h3(shiny::icon("oil-can"), "Lipidomics"),
                     p("Comprehensive lipid profiling and analysis")
                   )
            ),
            column(4,
                   shiny::div(
                     id = "integration_box",
                     class = "omic-selection-box disabled",
                     h3(shiny::icon("project-diagram"), "Integration"),
                     p("Multi-omics integration (available when 2+ omics selected)")
                   )
            )
          ),
          
          br(),
          
          fluidRow(
            column(12,
                   actionButton(
                     "start_analysis",
                     "Start Analysis",
                     class = "btn-primary btn-lg",
                     width = "100%",
                     icon = shiny::icon("play")
                   )
            )
          ),
          
          # Hidden inputs to track selections
          shinyjs::hidden(
            shiny::div(
              id = "selection_tracker",
              checkboxInput("proteomics_selected", "Proteomics", FALSE),
              checkboxInput("metabolomics_selected", "Metabolomics", FALSE),
              checkboxInput("transcriptomics_selected", "Transcriptomics", FALSE),
              checkboxInput("lipidomics_selected", "Lipidomics", FALSE),
              checkboxInput("integration_selected", "Integration", FALSE)
            )
          )
        ),
        
        # Dynamic tab items will be added here
        uiOutput("dynamic_tabs")
      ),
      
      # Log Terminal UI
      shiny::uiOutput("log_terminal_ui"),
      
      # JavaScript for handling omics selection - use shiny:: prefix
      shiny::tags$script(shiny::HTML("
      // ========== DARK MODE TOGGLE ==========
      function toggleDarkMode() {
        document.body.classList.toggle('dark-mode');
        var isDark = document.body.classList.contains('dark-mode');
        localStorage.setItem('darkMode', isDark ? 'true' : 'false');
        // Update icon
        var icon = document.querySelector('#dark_mode_toggle i');
        if (icon) {
          icon.className = isDark ? 'fa fa-sun' : 'fa fa-moon';
        }
      }

      // Apply saved preference on page load
      document.addEventListener('DOMContentLoaded', function() {
        var savedMode = localStorage.getItem('darkMode');
        // Default to dark mode if no preference saved (user prefers dark)
        if (savedMode === 'true' || savedMode === null) {
          document.body.classList.add('dark-mode');
          var icon = document.querySelector('#dark_mode_toggle i');
          if (icon) icon.className = 'fa fa-sun';
          if (savedMode === null) localStorage.setItem('darkMode', 'true');
        }
      });

      // ========== OMICS SELECTION ==========
      var selectedOmics = [];
      
      function toggleOmicSelection(omic) {
        if (omic === 'integration') return; // Integration is auto-enabled
        
        var box = document.getElementById(omic + '_box');
        var checkbox = document.getElementById(omic + '_selected');
        
        if (selectedOmics.includes(omic)) {
          selectedOmics = selectedOmics.filter(item => item !== omic);
          box.classList.remove('selected');
        } else {
          selectedOmics.push(omic);
          box.classList.add('selected');
        }
        
        // Update hidden checkboxes
        checkbox.checked = selectedOmics.includes(omic);
        checkbox.dispatchEvent(new Event('change'));
        
        // Enable/disable integration based on selection count
        var integrationBox = document.getElementById('integration_box');
        if (selectedOmics.length >= 2) {
          integrationBox.classList.remove('disabled');
          integrationBox.onclick = function() { toggleOmicSelection('integration'); };
        } else {
          integrationBox.classList.add('disabled');
          integrationBox.classList.remove('selected');
          integrationBox.onclick = null;
          selectedOmics = selectedOmics.filter(item => item !== 'integration');
          document.getElementById('integration_selected').checked = false;
          document.getElementById('integration_selected').dispatchEvent(new Event('change'));
        }
        
        Shiny.setInputValue('selected_omics', selectedOmics);
      }
    "))
    )
  )
}
