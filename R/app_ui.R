# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @import shinydashboard
#' @rawNamespace import(shinyjs, except = runExample)
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
        shiny::tags$img(src = "www/MultiScholaR_v2.png", style = "width: 80%; height: auto;")
      ),
      testid(
        sidebarMenu(
          id = "main_menu",
          menuItem("Home", tabName = "home", icon = shiny::icon("home")),

          # Dynamic menu items will be added based on omics selection
          uiOutput("dynamic_menu")
        ),
        "sidebar-menu"
      )
    ),
    dashboardBody(
      useShinyjs(),

      shiny::tags$head(
        shiny::tags$script(shiny::HTML(workflowPeptideSampleAckScript()))
      ),

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

        /* ========== WORKFLOW STEPPER ========== */
        .workflow-stepper {
          display: flex;
          align-items: flex-start;
          justify-content: center;
          padding: 24px 16px;
          margin-bottom: 20px;
          background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%);
          border-radius: 12px;
          box-shadow: 0 2px 8px rgba(0,0,0,0.06);
          overflow-x: auto;
        }
        .stepper-step {
          display: flex;
          flex-direction: column;
          align-items: center;
          position: relative;
          flex: 1;
          min-width: 80px;
          max-width: 140px;
        }
        .stepper-step:not(:last-child)::after {
          content: '';
          position: absolute;
          top: 20px;
          left: calc(50% + 24px);
          width: calc(100% - 48px);
          height: 3px;
          background-color: #dee2e6;
          z-index: 0;
          transition: background-color 0.4s ease;
        }
        .stepper-step.step-completed:not(:last-child)::after {
          background: linear-gradient(90deg, #27ae60 0%, #2ecc71 100%);
        }
        .stepper-step.step-current:not(:last-child)::after {
          background: linear-gradient(90deg, #3498db 0%, #dee2e6 100%);
        }
        .stepper-circle {
          width: 40px;
          height: 40px;
          border-radius: 50%;
          display: flex;
          align-items: center;
          justify-content: center;
          font-size: 14px;
          font-weight: 600;
          z-index: 1;
          transition: all 0.3s ease;
          border: 3px solid transparent;
        }
        .stepper-circle.completed {
          background: linear-gradient(135deg, #27ae60 0%, #2ecc71 100%);
          color: white;
          box-shadow: 0 3px 10px rgba(39, 174, 96, 0.35);
        }
        .stepper-circle.current {
          background: linear-gradient(135deg, #3498db 0%, #5dade2 100%);
          color: white;
          box-shadow: 0 3px 12px rgba(52, 152, 219, 0.45);
          animation: pulse-current 2s infinite;
        }
        .stepper-circle.pending {
          background-color: #fff;
          color: #adb5bd;
          border-color: #dee2e6;
        }
        @keyframes pulse-current {
          0%, 100% { box-shadow: 0 3px 12px rgba(52, 152, 219, 0.45); }
          50% { box-shadow: 0 3px 20px rgba(52, 152, 219, 0.7); }
        }
        .stepper-label {
          margin-top: 10px;
          font-size: 0.8em;
          font-weight: 500;
          color: #495057;
          text-align: center;
          white-space: nowrap;
          transition: color 0.3s ease;
        }
        .stepper-step.step-completed .stepper-label {
          color: #27ae60;
          font-weight: 600;
        }
        .stepper-step.step-current .stepper-label {
          color: #3498db;
          font-weight: 600;
        }
        .stepper-step.step-pending .stepper-label {
          color: #adb5bd;
        }

        /* Dark mode stepper */
        body.dark-mode .workflow-stepper {
          background: linear-gradient(135deg, #16213e 0%, #1a1a2e 100%);
          box-shadow: 0 2px 12px rgba(0,0,0,0.3);
        }
        body.dark-mode .stepper-step:not(:last-child)::after {
          background-color: #3d4f6f;
        }
        body.dark-mode .stepper-step.step-completed:not(:last-child)::after {
          background: linear-gradient(90deg, #2ecc71 0%, #27ae60 100%);
        }
        body.dark-mode .stepper-step.step-current:not(:last-child)::after {
          background: linear-gradient(90deg, #5dade2 0%, #3d4f6f 100%);
        }
        body.dark-mode .stepper-circle.pending {
          background-color: #1a1a2e;
          color: #6c7a89;
          border-color: #3d4f6f;
        }
        body.dark-mode .stepper-label {
          color: #b0b0b0;
        }
        body.dark-mode .stepper-step.step-completed .stepper-label {
          color: #2ecc71;
        }
        body.dark-mode .stepper-step.step-current .stepper-label {
          color: #5dade2;
        }
        body.dark-mode .stepper-step.step-pending .stepper-label {
          color: #6c7a89;
        }
        /* ========== RUV DESIGN SECTION (normalisation tab; mod_omic_norm_panels.R) ========== */
        .ruv-design {
          --ruv-rule: #d5d5d5;
          --ruv-muted: #555;
          --ruv-selected: #d1ecf1;
          --ruv-accent: #3c8dbc;
          --ruv-columns: 140px minmax(0, 1.7fr) 170px 110px minmax(0, 1.5fr);
          background-color: #fff;
        }
        body.dark-mode .ruv-design {
          --ruv-rule: #33466b;
          --ruv-muted: #b8c0cc;
          --ruv-selected: #1a4a6e;
          --ruv-accent: #5dade2;
        }
        .ruv-design > h4 {
          margin-top: 0;
        }
        .ruv-design .help-block {
          max-width: 110ch;
          color: var(--ruv-muted);
        }
        /* The designs: one comparison list, a row per radio choice, the columns shared with the headings. */
        .ruv-design-header,
        .ruv-design-option {
          display: grid;
          grid-template-columns: var(--ruv-columns);
          column-gap: 24px;
          align-items: start;
        }
        .ruv-design-header {
          align-items: end;
          padding: 6px 12px 6px 35px;
          border-bottom: 2px solid var(--ruv-rule);
          color: var(--ruv-muted);
          font-size: 12px;
          font-weight: 600;
        }
        /* Shiny pulls a radio group up under its label by the rows' own top margin, which these rows do not have. */
        .ruv-design .ruv-design-header + .shiny-options-group {
          margin-top: 0;
          margin-bottom: 12px;
        }
        .ruv-design-header + .shiny-options-group > .radio {
          margin: 0;
          padding: 10px 12px;
          border-bottom: 1px solid var(--ruv-rule);
          border-left: 3px solid transparent;
        }
        .ruv-design-header + .shiny-options-group > .radio:has(input:checked) {
          background-color: var(--ruv-selected);
          border-left-color: var(--ruv-accent);
        }
        .ruv-design-header + .shiny-options-group > .radio > label {
          display: block;
        }
        .ruv-design-header + .shiny-options-group > .radio:has(input:disabled) > label {
          cursor: not-allowed;
        }
        /* A design: what is done, then how on a second, lighter line. */
        .ruv-design-option > .ruv-design-label {
          font-weight: 600;
        }
        .ruv-design-description {
          display: block;
          margin-top: 2px;
          color: var(--ruv-muted);
          font-size: 13px;
          font-weight: 400;
          line-height: 1.45;
        }
        .ruv-design-header .ruv-design-where {
          font-weight: 400;
        }
        .ruv-design-fdp > strong,
        .ruv-design-power > strong {
          display: block;
          font-size: 16px;
          line-height: 20px;
          font-variant-numeric: tabular-nums;
        }
        .ruv-design-where,
        .ruv-design-note,
        .ruv-design-status {
          color: var(--ruv-muted);
          font-size: 13px;
          line-height: 1.45;
        }
        .ruv-design-where {
          display: block;
          font-size: 12px;
        }
        .ruv-design-status {
          grid-column: span 2;
        }
        /* A row's figures span the two figure columns and keep them; the figure with the other rule for k is a
           line under the two cells. */
        .ruv-design-figures {
          display: grid;
          grid-column: span 2;
          grid-template-columns: subgrid;
          align-content: start;
        }
        .ruv-design-alt {
          grid-column: span 2;
          margin-top: 6px;
          padding-top: 4px;
          border-top: 1px dotted var(--ruv-rule);
          color: var(--ruv-muted);
          font-size: 12px;
          line-height: 1.45;
        }
        .ruv-design-alt > span {
          display: block;
        }
        .ruv-design-option.is-pending .ruv-design-status {
          font-style: italic;
        }
        .ruv-design-option.is-unavailable .ruv-design-label {
          color: var(--ruv-muted);
        }
        .ruv-design-tier > strong {
          display: inline-block;
          padding: 1px 8px;
          border: 1px solid var(--ruv-muted);
          border-radius: 3px;
          color: var(--ruv-muted);
          font-size: 12px;
          line-height: 18px;
          white-space: nowrap;
        }
        .ruv-design-tier > strong.ruv-design-tier-recommended {
          border-color: #1e7b45;
          background-color: #1e7b45;
          color: #fff;
        }
        .ruv-design-tier > strong.ruv-design-tier-legacy {
          border-color: #8a5a00;
          color: #8a5a00;
        }
        body.dark-mode .ruv-design-tier > strong.ruv-design-tier-legacy {
          border-color: #f0c36a;
          color: #f0c36a;
        }
        /* The settings of the selected design: the k rules side by side, rows aligned; today's controls. */
        .ruv-design-k,
        .ruv-design-today {
          margin-top: 16px;
        }
        .ruv-design-k > h5 {
          font-weight: 600;
        }
        .ruv-design-k-choice {
          display: flex;
          flex-wrap: wrap;
          align-items: flex-start;
          column-gap: 3em;
        }
        .ruv-design-k-choice .shiny-input-container {
          width: auto !important;
        }
        .ruv-design-k-choice .shiny-options-group {
          display: flex;
          flex-wrap: wrap;
          column-gap: 2em;
          margin-bottom: 0;
        }
        .ruv-design .ruv-design-k-choice .shiny-options-group {
          margin-top: 0;
        }
        .ruv-design-k-choice .shiny-options-group > .radio {
          margin: 0;
        }
        .ruv-design-k-choice input[type='number'] {
          width: 8em;
        }
        /* The rules' explanation: a native disclosure, closed until opened. */
        .ruv-design-k-why {
          margin: 8px 0 10px;
        }
        .ruv-design-k-why > summary {
          display: list-item;
          width: fit-content;
          padding: 4px 0;
          color: var(--ruv-accent);
          font-weight: 600;
          cursor: pointer;
        }
        .ruv-design-k-why > summary:focus-visible {
          outline: 2px solid var(--ruv-accent);
          outline-offset: 2px;
        }
        .ruv-design-k-rules {
          display: grid;
          grid-auto-flow: column;
          grid-auto-columns: minmax(0, 1fr);
          grid-template-rows: repeat(5, auto);
          column-gap: 2em;
          margin: 4px 0 12px;
          padding-top: 10px;
          border-top: 1px solid var(--ruv-rule);
          max-width: 180ch;
        }
        .ruv-design-k-rule {
          display: grid;
          grid-row: span 5;
          grid-template-rows: subgrid;
        }
        .ruv-design-k-rule > p {
          margin: 0 0 8px;
          font-size: 13px;
          line-height: 1.45;
        }
        .ruv-design-k-rule > p.ruv-design-k-rule-name {
          font-size: 14px;
        }
        .ruv-design-today {
          max-width: 34em;
        }
        @media (max-width: 1279px) {
          /* Narrow window: the tier, the design and its note share the first column; the figures stay aligned. */
          .ruv-design {
            --ruv-columns: minmax(0, 1fr) 160px 110px;
          }
          .ruv-design-header > .ruv-design-tier,
          .ruv-design-header > .ruv-design-note {
            display: none;
          }
          .ruv-design-option > .ruv-design-tier,
          .ruv-design-option > .ruv-design-label,
          .ruv-design-option > .ruv-design-note {
            grid-column: 1;
          }
          .ruv-design-option > .ruv-design-tier {
            margin-bottom: 4px;
          }
          .ruv-design-option > .ruv-design-note {
            margin-top: 4px;
          }
          .ruv-design-option > .ruv-design-figures {
            grid-column: 2 / span 2;
            grid-row: 1 / span 3;
          }
          /* The k rules one under another, each rule's name beside its texts. */
          .ruv-design-k-rules {
            display: block;
            max-width: 100ch;
          }
          .ruv-design-k-rule {
            grid-template-columns: 9em minmax(0, 1fr);
            grid-template-rows: none;
            column-gap: 1.5em;
            padding-bottom: 6px;
          }
          .ruv-design-k-rule > p {
            grid-column: 2;
          }
          .ruv-design-k-rule > p.ruv-design-k-rule-name {
            grid-column: 1;
            grid-row: 1 / span 4;
          }
        }
      "))
      ),
      shinydashboard::tabItems(
        # Home tab with omics selection
        shinydashboard::tabItem(
          tabName = "home",
          fluidRow(
            column(
              12,
              h2("Welcome to MultiScholaR Interactive Workflow"),
              p("Select the type(s) of omics data you want to analyze:"),
              br()
            )
          ),
          fluidRow(
            column(
              4,
              testid(
                shiny::div(
                  id = "proteomics_box",
                  class = "omic-selection-box",
                  onclick = "toggleOmicSelection('proteomics')",
                  h3(shiny::icon("dna"), "Proteomics"),
                  p("Analyze protein abundance data from DIA-NN or similar tools")
                ),
                "tile-proteomics"
              )
            ),
            column(
              4,
              testid(
                shiny::div(
                  id = "metabolomics_box",
                  class = "omic-selection-box",
                  onclick = "toggleOmicSelection('metabolomics')",
                  h3(shiny::icon("flask"), "Metabolomics"),
                  p("Process and analyze metabolite profiling data")
                ),
                "tile-metabolomics"
              )
            ),
            column(
              4,
              testid(
                shiny::div(
                  id = "transcriptomics_box",
                  class = "omic-selection-box",
                  onclick = "toggleOmicSelection('transcriptomics')",
                  h3(shiny::icon("chart-line"), "Transcriptomics"),
                  p("Analyze gene expression data from RNA-seq")
                ),
                "tile-transcriptomics"
              )
            )
          ),
          fluidRow(
            column(
              4,
              testid(
                shiny::div(
                  id = "lipidomics_box",
                  class = "omic-selection-box",
                  onclick = "toggleOmicSelection('lipidomics')",
                  h3(shiny::icon("oil-can"), "Lipidomics"),
                  p("Comprehensive lipid profiling and analysis")
                ),
                "tile-lipidomics"
              )
            ),
            column(
              4,
              # KIT-138: with the gate `MultiScholaR.integration_gui` on (the default since KIT-240), the tile is live
              # (mod_integration_*). Set FALSE: the placeholder selector as before (finding m-2).
              if (integrationGuiEnabled()) integrationHomeTile() else testid(
                shiny::div(
                  id = "integration_box",
                  class = "omic-selection-box disabled",
                  h3(shiny::icon("project-diagram"), "Integration"),
                  p("Multi-omics integration (available when 2+ omics selected)")
                ),
                "tile-integration"
              )
            )
          ),
          br(),
          fluidRow(
            column(
              12,
              testid(
                actionButton(
                  "start_analysis",
                  "Start Analysis",
                  class = "btn-primary btn-lg",
                  width = "100%",
                  icon = shiny::icon("play")
                ),
                "btn-start-analysis"
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
        testid(uiOutput("dynamic_tabs"), "dynamic-tabs-container")
      ),

      # Log Terminal UI
      shiny::uiOutput("log_terminal_ui"),

      # App-level state digest for test automation; NULL uses identity NS (no module prefix)
      test_mode_digest_ui(NULL),

      # JavaScript for handling omics selection - use shiny:: prefix
      # KIT-138: the selection script as it is; with the integration gate on, Integration is selectable.
      shiny::tags$script(shiny::HTML(integrationHomeScript("
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
    ")))
    )
  )
}
