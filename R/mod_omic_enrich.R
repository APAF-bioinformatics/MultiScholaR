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

# KIT-107 W3 (D2, D3, D7, D11): the shared enrichment tab (kit-101 design 4.9), from proteomics' enrichment module at
# 6129b546 (mod_prot_enrich.R's markup and the server's setup, registration, observer and state layers, deleted by W3).
# `mod_omic_enrich_ui(id, config)` renders an omic's `enrichment` group: the cutoffs', correction's and taxon's values,
# the organism-filter section only when the group has `organism_filter`, and the result tabs `services` names (the
# markup is today's, byte-identical: W1 (a)). `mod_omic_enrich_server(...)` holds the tab's reactive values (today's 13
# fields and initial values), the method reactive, the taxon prefill and the mixed-species observer, the tab-selected and
# DA-arrival observers (the DA result manifest only: the legacy non-manifest branch is retired, D7), the
# selected-contrast observer with the restore, the outputs `services` names, the Run and the download, registered in
# today's order. The Run calls the core (`runWorkflowEnrichment`, utils_workflow_enrichment.R) with the manifest the tab
# holds; the tab keeps that manifest across runs (F107-01) and writes its values from the core's returned list. A DA
# table missing at the Run is refused before anything is written (F107-03, the lead's decision). Functions a test mocks
# are called by name, or read from the configuration built when the server starts. The module frame keeps no locals
# (the testServer masking rule; test-omic-design-module.R). Recorded divergences from conventions/R-shiny.md, kept under
# conventions/R.md's existing-code scope rule: the run is not an ExtendedTask (KIT-113 moves it); the wellPanel/fluidRow
# markup is kept byte-identical; the moved functions keep their `@` slot reads. Recorded: a blank taxon stops the Run
# silently at its `req()` before any notice (the lead's W3 instruction); the console traces of the retired legacy
# branch are not written.

#' The enrichment tab's UI: the layout of an omic's `enrichment` group.
#' @param id The module id.
#' @param config The omic's `enrichment` group (`omicGuiConfig(<omic>)$enrichment`).
#' @return A tag list.
#' @noRd
mod_omic_enrich_ui <- function(id, config) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::wellPanel(
      shiny::fluidRow(
        omicEnrichSettingsColumn(ns, config),
        shiny::column(9,
          do.call(shiny::tabsetPanel, c(list(id = ns("enrichment_method_tabs")),
            lapply(config$services, omicEnrichServiceTab, ns = ns)))
        )
      )
    )
  )
}

#' The settings column: the contrast, the cutoffs, the taxon, the organism filter (when configured), the method and
#' contrasts displays, Run, the download, the status and the correction method.
#' @noRd
omicEnrichSettingsColumn <- function(ns, config) {
  cutoffs <- config$cutoffs
  shiny::column(3,
    shiny::wellPanel(
      shiny::h4("Enrichment Analysis Settings"),

      # Contrast selector
      shiny::selectInput(
        ns("selected_contrast"),
        "Select Contrast:",
        choices = NULL,
        width = "100%"
      ),
      shiny::helpText("Choose which DE contrast to analyze"),

      shiny::hr(),

      # Only the parameters that ACTUALLY EXIST in processEnrichments
      shiny::h5("Enrichment Settings"),
      shiny::numericInput(
        ns("up_cutoff"),
        "Up Log2FC Cutoff:",
        value = cutoffs$up$value,
        min = cutoffs$up$min,
        max = cutoffs$up$max,
        step = cutoffs$up$step
      ),
      shiny::helpText("Minimum log2 fold change for up-regulated genes"),

      shiny::numericInput(
        ns("down_cutoff"),
        "Down Log2FC Cutoff:",
        value = cutoffs$down$value,
        min = cutoffs$down$min,
        max = cutoffs$down$max,
        step = cutoffs$down$step
      ),
      shiny::helpText("Minimum log2 fold change for down-regulated genes"),

      shiny::numericInput(
        ns("q_cutoff"),
        "Q-value Cutoff:",
        value = cutoffs$q$value,
        min = cutoffs$q$min,
        max = cutoffs$q$max,
        step = cutoffs$q$step
      ),
      shiny::helpText("FDR threshold for enrichment significance"),

      shiny::textInput(
        ns("organism_taxid"),
        "Organism Taxon ID:",
        value = config$taxon$default,  # Default to human, will be updated by server
        placeholder = "e.g., 9606 for human"
      ),
      shiny::helpText("NCBI taxonomy ID for organism"),

      # Mixed species filtering checkbox, when the omic filters by organism
      if (!is.null(config$organism_filter)) shiny::tagList(
        shiny::hr(),
        shiny::h5("Multi-Species Filtering"),
        shiny::checkboxInput(
          ns("enable_organism_filter"),
          "Filter to target organism only",
          value = FALSE
        ),
        shiny::helpText("Enable if using mixed-species FASTA (e.g., with contaminants). Filters DE results to target organism before enrichment analysis.")
      ),

      # Analysis method display
      shiny::h5("Analysis Method"),
      shiny::verbatimTextOutput(ns("analysis_method_display")),
      shiny::helpText("Automatically determined based on organism support"),

      shiny::hr(),

      # Available contrasts display
      shiny::h5("Available Contrasts"),
      shiny::verbatimTextOutput(ns("contrasts_display")),
      shiny::helpText("Contrasts from DE analysis"),

      shiny::br(),

      # Main action button
      shiny::actionButton(
        ns("run_enrichment_analysis"),
        "Run Enrichment Analysis",
        class = "btn-primary",
        width = "100%",
        icon = shiny::icon("play")
      ),

      shiny::br(),
      shiny::br(),

      # Download button
      shiny::downloadButton(
        ns("download_enrichment_results"),
        "Download All Results",
        class = "btn-success",
        width = "100%"
      ),

      shiny::br(),
      shiny::br(),

      # Status display
      shiny::h5("Analysis Status"),
      shiny::verbatimTextOutput(ns("enrichment_status")),

      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'gprofiler2'", ns("enrichment_method_tabs")),
        shiny::selectInput(
          ns("correction_method"),
          "Correction Method:",
          choices = as.list(config$correction$choices),
          selected = config$correction$selected
        ),
        shiny::helpText("g:SCS: More conservative, reduces false positives but may miss terms."),
        shiny::helpText("FDR: Detects more terms but higher false positive risk.")
      )
    )
  )
}

#' One result tab: g:Profiler, clusterProfiler or the STRING-DB placeholder (the markup at 6129b546, unchanged).
#' @noRd
omicEnrichServiceTab <- function(service, ns) {
  switch(service,
    gprofiler2 =
        shiny::tabPanel(
          "gprofiler2 Analysis",
          value = "gprofiler2",
          icon = shiny::icon("globe"),
          shiny::br(),
          
          shiny::div(
            style = "background-color: #e8f4fd; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
            shiny::h5("gprofiler2 Enrichment", style = "color: #2c5282; margin-top: 0;"),
            shiny::p("Comprehensive functional enrichment using gprofiler2 for well-supported model organisms. 
                     Includes GO terms, KEGG pathways, Reactome pathways, and more.", 
                     style = "margin-bottom: 0; color: #2c5282;")
          ),
          
          # gprofiler2 results display
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] == 'gprofiler2'", ns("enrichment_method_tabs")),
            
            # Direction filter for viewing results
            shiny::fluidRow(
              shiny::column(6,
                shiny::selectInput(
                  ns("gprofiler_direction_filter"),
                  "Direction:",
                  choices = list(
                    "All" = "all",
                    "Up-regulated only" = "up",
                    "Down-regulated only" = "down"
                  ),
                  selected = "all"
                )
              ),
              shiny::column(6,
                shiny::div(
                  style = "padding: 20px 0; text-align: center; color: #666;",
                  shiny::strong("Select direction to filter results")
                )
              )
            ),
            
            shiny::hr(),
            
            # [OK] NEW: Add plot display for gprofiler2
            shiny::fluidRow(
              shiny::column(12,
                shiny::h5("Enrichment Plot"),
                shinyjqui::jqui_resizable(
                  plotly::plotlyOutput(ns("gprofiler_plot"), height = "500px", width = "100%")
                ),
                shiny::br()
              )
            ),
            
            # gprofiler2 results table
            DT::DTOutput(ns("gprofiler_results_table")),
            
            shiny::br(),
            
            # gprofiler2 summary stats
            shiny::wellPanel(
              shiny::h5("gprofiler2 Enrichment Summary"),
              shiny::verbatimTextOutput(ns("gprofiler_summary_stats"))
            )
          )
        ),
    clusterprofiler =
        shiny::tabPanel(
          "clusterProfileR Analysis",
          value = "clusterprofiler",
          icon = shiny::icon("project-diagram"),
          shiny::br(),
          
          shiny::div(
            style = "background-color: #f0fff4; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
            shiny::h5("clusterProfileR Enrichment", style = "color: #22543d; margin-top: 0;"),
            shiny::p("Custom enrichment analysis using UniProt GO annotations via clusterProfileR. 
                     Ideal for organisms not supported by gprofiler2.", 
                     style = "margin-bottom: 0; color: #22543d;")
          ),
          
          # clusterProfileR results display
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] == 'clusterprofiler'", ns("enrichment_method_tabs")),
            
            # Direction filter for viewing results
            shiny::fluidRow(
              shiny::column(6,
                shiny::selectInput(
                  ns("clusterprofiler_direction_filter"),
                  "Direction:",
                  choices = list(
                    "All" = "all",
                    "Up-regulated only" = "up",
                    "Down-regulated only" = "down"
                  ),
                  selected = "all"
                )
              ),
              shiny::column(6,
                shiny::div(
                  style = "padding: 20px 0; text-align: center; color: #666;",
                  shiny::strong("Select direction to filter results")
                )
              )
            ),
            
            shiny::hr(),
            
            # [OK] NEW: Add plot display for clusterProfileR
            shiny::fluidRow(
              shiny::column(12,
                shiny::h5("Enrichment Plot"),
                shinyjqui::jqui_resizable(
                  plotly::plotlyOutput(ns("clusterprofiler_plot"), height = "500px", width = "100%")
                ),
                shiny::br()
              )
            ),
            
            # clusterProfileR results table
            DT::DTOutput(ns("clusterprofiler_results_table")),
            
            shiny::br(),
            
            # clusterProfileR summary stats
            shiny::wellPanel(
              shiny::h5("clusterProfileR Enrichment Summary"),
              shiny::verbatimTextOutput(ns("clusterprofiler_summary_stats"))
            )
          )
        ),
    stringdb =
        shiny::tabPanel(
          "STRING-DB Networks",
          value = "stringdb",
          icon = shiny::icon("network-wired"),
          shiny::br(),
          
          shiny::div(
            style = "background-color: #fef5e7; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
            shiny::h5("STRING-DB Network Analysis", style = "color: #744210; margin-top: 0;"),
            shiny::p("Protein-protein interaction network enrichment using STRING database. 
                     Identifies functional networks and interaction clusters.", 
                     style = "margin-bottom: 0; color: #744210;")
          ),
          
          # STRING-DB results display
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] == 'stringdb'", ns("enrichment_method_tabs")),
            
            # Filter controls for STRING-DB
            shiny::fluidRow(
              shiny::column(4,
                shiny::selectInput(
                  ns("stringdb_ranking_method"),
                  "Ranking Method:",
                  choices = list(
                    "Combined Score" = "combined_score",
                    "FDR Q-value" = "fdr_qvalue",
                    "Log2 Fold Change" = "log2fc"
                  ),
                  selected = "combined_score"
                )
              ),
              shiny::column(4,
                shiny::checkboxInput(
                  ns("stringdb_filter_significant"),
                  "Filter Significant Only",
                  value = FALSE
                )
              ),
              shiny::column(4,
                shiny::numericInput(
                  ns("stringdb_max_results"),
                  "Max Results:",
                  value = 50,
                  min = 10,
                  max = 500,
                  step = 10
                )
              )
            ),
            
            shiny::hr(),
            
            # [OK] NEW: Add plot display for STRING-DB
            shiny::fluidRow(
              shiny::column(12,
                shiny::h5("Network Plot"),
                shinyjqui::jqui_resizable(
                  plotly::plotlyOutput(ns("stringdb_plot"), height = "500px", width = "100%")
                ),
                shiny::br()
              )
            ),
            
            # STRING-DB results table
            DT::DTOutput(ns("stringdb_results_table")),
            
            shiny::br(),
            
            # STRING-DB summary stats  
            shiny::wellPanel(
              shiny::h5("STRING-DB Network Summary"),
              shiny::verbatimTextOutput(ns("stringdb_summary_stats"))
            )
          )
        ))
}

#' The enrichment tab's server.
#' @param id,workflow_data,experiment_paths,omic_type,experiment_label,selected_tab The tab's arguments, as
#'   `mod_prot_enrich_server`'s.
#' @param config The omic's `enrichment` group (`omicGuiConfig(<omic>)$enrichment`).
#' @return The tab's reactive values.
#' @noRd
mod_omic_enrich_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, config,
                                   selected_tab = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    runOmicEnrichServer(input, output, session, id, workflow_data, experiment_paths, config, selected_tab)
  })
}

#' The tab's reactive values: today's 13 fields and initial values (mod_prot_enrich_observer_helpers.R:190-206 at
#' 6129b546). The run adds the display tables, the service, the product and the full results.
#' @noRd
omicEnrichTabState <- function() {
  shiny::reactiveValues(
    enrichment_results = NULL,
    contrasts_available = NULL,
    analysis_complete = FALSE,
    current_s4_object = NULL,
    da_results_data = NULL,
    gprofiler_results = NULL,
    clusterprofiler_results = NULL,
    stringdb_results = NULL,
    analysis_method = NULL,
    organism_supported = NULL,
    all_enrichment_results = list(),
    current_contrast_results = list(),
    enrichment_plots = list()
  )
}

#' The server's body: the reactive values, the method reactive, the observers and the outputs in today's order
#' (mod_prot_enrich_server.R and its setup layer at 6129b546).
#' @noRd
runOmicEnrichServer <- function(input, output, session, id, workflowData, experimentPaths, config, selectedTab) {
  cat("--- Entering mod_prot_enrich_server ---\n")
  cat(sprintf("   mod_prot_enrich_server Arg: id = %s\n", id))
  cat(sprintf("   mod_prot_enrich_server Arg: workflow_data is NULL = %s\n", is.null(workflowData)))
  cat(sprintf("   mod_prot_enrich_server Arg: selected_tab is NULL = %s\n", is.null(selectedTab)))
  cat("=== ENRICHMENT ANALYSIS MODULE SERVER STARTED ===\n")
  cat(sprintf("Module ID: %s\n", id))

  tab <- omicEnrichTabState()
  cat("   mod_prot_enrich_server Step: Reactive values initialized\n")
  supported <- shiny::reactive(config[["organism"]]$supported())
  method <- createProtEnrichCurrentAnalysisMethodReactive(input = input, enrichmentData = tab,
    supportedOrganismsFn = supported, resolveAnalysisMethodFn = config[["organism"]]$method)

  registerOmicEnrichBusObservers(session, workflowData, config)
  registerOmicEnrichSelectedContrastObserver(input, tab, workflowData, config)
  if (!is.null(selectedTab)) {
    cat("   mod_prot_enrich_server Step: Setting up tab selection observer\n")
    shiny::observeEvent(selectedTab(), selectOmicEnrichTab(session, selectedTab(), tab, workflowData, config),
      ignoreInit = TRUE)
  } else {
    cat("   mod_prot_enrich_server Step: No selected_tab parameter provided - tab selection observer NOT set up\n")
  }
  shiny::observeEvent(workflowData$da_analysis_results_list, {
    cat("*** ENRICHMENT: DE results detected - updating contrasts ***\n")
    holdOmicEnrichManifest(session, tab, workflowData$da_analysis_results_list)
  }, ignoreInit = TRUE, ignoreNULL = TRUE)

  registerOmicEnrichStatusOutputs(output, input, tab, method)
  shiny::observeEvent(input$run_enrichment_analysis,
    runOmicEnrichAnalysis(input, tab, workflowData, experimentPaths, config))
  registerOmicEnrichResultOutputs(output, input, tab, workflowData, config)
  output$download_enrichment_results <- shiny::downloadHandler(
    filename = function() buildProtEnrichResultsDownloadFilename(input$selected_contrast),
    content = function(file) {
      writeProtEnrichResultsDownloadArchive(file = file, selectedContrast = input$selected_contrast,
        methodInfo = method(), organismTaxid = input$organism_taxid, upCutoff = input$up_cutoff,
        downCutoff = input$down_cutoff, qCutoff = input$q_cutoff, gprofilerResults = tab$gprofiler_results,
        clusterprofilerResults = tab$clusterprofiler_results, stringdbResults = tab$stringdb_results)
    })

  tab
}

# ------------------------------------------------------------------------------------------------
# The observers.
# ------------------------------------------------------------------------------------------------

#' The taxon prefill from the bus and, when the omic filters by organism, the mixed-species observer that switches the
#' filter on (mod_prot_enrich_registration_helpers.R:81-153 at 6129b546).
#' @noRd
registerOmicEnrichBusObservers <- function(session, workflowData, config) {
  shiny::observeEvent(workflowData[[config$taxon$prefill_from]], {
    taxon <- workflowData[[config$taxon$prefill_from]]
    if (!is.null(taxon)) shiny::updateTextInput(session, "organism_taxid", value = taxon)
  }, ignoreInit = TRUE, ignoreNULL = TRUE)
  if (is.null(config$organism_filter)) return(invisible(NULL))
  shiny::observeEvent(workflowData[[config$organism_filter$mixed_species_from]], {
    mixed <- workflowData[[config$organism_filter$mixed_species_from]]
    if (!is.null(mixed) && isTRUE(mixed$enabled)) {
      shiny::updateCheckboxInput(session, "enable_organism_filter", value = TRUE)
      cat("*** ENRICHMENT: Auto-enabled organism filter (mixed species FASTA detected at import) ***\n")
      shiny::showNotification(sprintf("Multi-species FASTA detected. Filtering to %s enabled.", mixed$selected_organism),
        type = "message", duration = 5)
    }
  }, ignoreInit = TRUE, ignoreNULL = TRUE)
}

#' The selected contrast: after an Open the results are restored from the completed state on first use (KIT-024,
#' F107-02); then the contrast's results are shown (mod_prot_enrich_registration_helpers.R:580-671 at 6129b546).
#' @noRd
registerOmicEnrichSelectedContrastObserver <- function(input, tab, workflowData, config) {
  shiny::observe({
    shiny::req(input$selected_contrast)
    restoreProtEnrichResultsFromState(tab, workflowData, config)
    shiny::req(tab$analysis_complete)
    shiny::req(tab$all_enrichment_results)
    cat(sprintf("*** CONTRAST CHANGED: User selected '%s' ***\n", input$selected_contrast))
    shown <- resolveProtEnrichSelectedContrastResults(input$selected_contrast, tab$all_enrichment_results,
      resolveProtEnrichContrastsTable(workflowData))
    showOmicEnrichContrast(tab, shown)
  })
}

#' Show one contrast's results (`resolveProtEnrichSelectedContrastResults`), or clear them when it has none.
#' @noRd
showOmicEnrichContrast <- function(tab, shown) {
  tab$gprofiler_results <- if (shown$found) shown$gprofilerResults
  tab$clusterprofiler_results <- if (shown$found) shown$clusterprofilerResults
  tab$stringdb_results <- if (shown$found) shown$stringdbResults
  tab$gprofiler_results_table <- if (shown$found) shown$gprofilerResultsTable
  tab$clusterprofiler_results_table <- if (shown$found) shown$clusterprofilerResultsTable
  invisible(shown)
}

#' The DA result manifest as the tab holds it (its source reference, the manifest, the contrasts), and the contrast
#' selector's choices. Anything else is not a DA result the tab can use (D7): nothing is held.
#' @return TRUE when a manifest was held.
#' @noRd
holdOmicEnrichManifest <- function(session, tab, references) {
  if (!isDIAResultManifest(references)) return(FALSE)
  tab$current_s4_object <- references$source
  tab$da_results_data <- references
  tab$contrasts_available <- references$contrasts$friendly_names
  shiny::updateSelectInput(session, "selected_contrast", choices = omicEnrichManifestChoices(references))
  TRUE
}

#' The contrast selector's choices for a DA result manifest: its contrasts' friendly names.
#' @noRd
omicEnrichManifestChoices <- function(manifest) {
  stats::setNames(manifest$contrasts$friendly_names, manifest$contrasts$friendly_names)
}

#' The shell selected a tab (mod_prot_enrich_registration_helpers.R:308-488 at 6129b546, the manifest branch): the
#' enrichment tab holds the DA manifest, or warns when the owner has no DA results yet. The legacy non-manifest
#' initialisation and its "invalid state" warning are retired (D7).
#' @noRd
selectOmicEnrichTab <- function(session, selected, tab, workflowData, config) {
  cat("--- Entering enrichment tab selection observer ---\n")
  cat(sprintf("   enrichment_tab_observer Step: Selected tab = %s\n", selected))
  if (is.null(selected) || selected != config$tab_value) {
    cat(sprintf("   enrichment_tab_observer Step: Tab '%s' is not enrichment tab, ignoring\n", selected))
  } else if (holdOmicEnrichManifest(session, tab, workflowData$da_analysis_results_list)) {
    return(invisible(TRUE))
  } else if (is.null(workflowData$state_manager)) {
    cat("*** workflow_data$state_manager is NULL - cannot check state ***\n")
  } else if (is.null(workflowData$da_analysis_results_list)) {
    cat("*** No DE analysis results found. User needs to complete differential expression analysis first. ***\n")
    shiny::showNotification("Please complete the differential expression analysis before accessing enrichment analysis.",
      type = "warning", duration = 5)
  }
  cat("--- Exiting enrichment tab selection observer ---\n")
  invisible(FALSE)
}

#' The Run (mod_prot_enrich_observer_helpers.R:1-162 at 6129b546, the manifest path): the working notice, the core's run
#' on the manifest the tab holds, inside the progress, then the tab's values from the core's returned list and the
#' completion notice; on failure the error notice. The tab keeps its manifest (F107-01); the source reference is held
#' after the run, as today's `on.exit` held it. A DA table missing at the Run is refused (F107-03).
#' KIT-113 D3: on a managed session the omic's runner (the one the normalisation shell installs; the frame installs none)
#' takes the run to the worker: a click while the runner's step runs or holds a pending publication is refused, then the
#' request half in the controller, one `enrichment` action on the manifest's source, and the accept half when it is
#' delivered. With no runner the run is in process, as before.
#' @noRd
runOmicEnrichAnalysis <- function(input, tab, workflowData, experimentPaths, config) {
  cat("=== STARTING ENRICHMENT ANALYSIS ===\n")
  shiny::req(input$selected_contrast, tab$da_results_data)
  shiny::req(input$organism_taxid)
  shiny::showNotification("Running enrichment analysis...", id = "enrichment_working", duration = NULL)
  operation <- config$operation()$name
  omic <- workflowActionSpec(operation)$omic_type
  runner <- workflowOmicsManagedActions(omic)
  if (!is.null(runner)) {
    # KIT-112's F112-02 guard: the runner would queue the action behind the running one (or its retry), and a DA run
    # completing in between would rewrite the tables the worker reads under this manifest's identity.
    if (!is.null(shiny::isolate(runner$step$active())) || !is.null(shiny::isolate(runner$step$pending()))) {
      shiny::showNotification(paste("A", omic, "step action is already running"), type = "warning")
      shiny::removeNotification("enrichment_working")
      return(invisible(FALSE))
    }
    manifest <- tab$da_results_data
    selected <- input$selected_contrast
    report <- function(e) {
      message(sprintf("*** ERROR in enrichment analysis: %s", e$message))
      cat(sprintf("*** ERROR in enrichment analysis: %s ***\n", e$message))
      shiny::showNotification(sprintf("Error in enrichment analysis: %s", e$message), type = "error", duration = 10)
    }
    progress <- runner$progress("Running enrichment analysis...")
    advance <- function(value, detail) progress$inc(value, detail = detail)
    # The tab holds the source reference after the run, success or failure, as in process.
    finish <- function(completed) {
      force(completed)
      progress$close()
      tab$current_s4_object <- manifest$source
      if (completed) {
        shiny::showNotification(paste("Enrichment analysis completed successfully for contrast:", selected),
          type = "message", duration = 5)
        cat("=== ENRICHMENT ANALYSIS COMPLETED ===\n")
      }
      shiny::removeNotification("enrichment_working")
      invisible(completed)
    }
    request <- tryCatch(requestWorkflowEnrichment(workflowData, config, manifest,
      parameters = list(selected_contrast = selected, organism_taxid = input$organism_taxid,
        up_cutoff = input$up_cutoff, down_cutoff = input$down_cutoff, q_cutoff = input$q_cutoff,
        correction_method = input$correction_method, enable_organism_filter = input$enable_organism_filter),
      paths = experimentPaths,
      notify = function(message, type, duration) shiny::showNotification(message, type = type, duration = duration),
      progress = advance), error = function(e) {
        report(e)
        NULL
      })
    if (is.null(request)) return(finish(FALSE))
    delivered <- NULL
    # The context names the controller's commit, which the handler's `commit = NULL` replaces: the accept half commits.
    runner$run(list(list(operation = operation, parameters = request$parameters, input = manifest$source,
      deliver = function(result) delivered <<- result,
      context = list(state_name = config$state_name, config = list(), description = operation))),
      # The runner calls `done` and `failed` from its step's observers: the tab and the bus are read without depending.
      done = function() shiny::isolate(finish(tryCatch({
        holdOmicEnrichRun(tab, acceptWorkflowEnrichment(workflowData, config, request$context, delivered, advance))
        TRUE
      }, error = function(e) {
        report(e)
        FALSE
      }))),
      failed = function(e) shiny::isolate({
        report(e)
        finish(FALSE)
      }))
    return(invisible(TRUE))
  }
  completed <- tryCatch({
    shiny::withProgress(message = "Running enrichment analysis...", value = 0,
      holdOmicEnrichRun(tab, runWorkflowEnrichment(workflowData, config, tab$da_results_data,
        parameters = list(selected_contrast = input$selected_contrast, organism_taxid = input$organism_taxid,
          up_cutoff = input$up_cutoff, down_cutoff = input$down_cutoff, q_cutoff = input$q_cutoff,
          correction_method = input$correction_method, enable_organism_filter = input$enable_organism_filter),
        paths = experimentPaths,
        notify = function(message, type, duration) shiny::showNotification(message, type = type, duration = duration),
        progress = function(value, detail) shiny::incProgress(value, detail = detail),
        refuse_missing_da = TRUE)))
    TRUE
  }, error = function(e) {
    message(sprintf("*** ERROR in enrichment analysis: %s", e$message))
    cat(sprintf("*** ERROR in enrichment analysis: %s ***\n", e$message))
    shiny::showNotification(sprintf("Error in enrichment analysis: %s", e$message), type = "error", duration = 10)
    FALSE
  })
  # The tab holds the source reference after the run, success or failure (the start commit's `on.exit`).
  tab$current_s4_object <- tab$da_results_data$source
  if (completed) {
    shiny::showNotification(paste("Enrichment analysis completed successfully for contrast:", input$selected_contrast),
      type = "message", duration = 5)
    cat("=== ENRICHMENT ANALYSIS COMPLETED ===\n")
  }
  shiny::removeNotification("enrichment_working")
  invisible(completed)
}

#' The tab's values from the core's returned list (`runWorkflowEnrichment`), as the start commit's run wrote them.
#' @noRd
holdOmicEnrichRun <- function(tab, run) {
  if (!is.null(run$annotation_match)) tab$annotation_match_results <- run$annotation_match
  tab$enrichment_service <- run$service
  tab$enrichment_product <- run$reference
  tab$all_enrichment_results <- run$contrasts
  if (run$selected$found) showOmicEnrichContrast(tab, run$selected)
  tab$enrichment_results_full <- run$results
  tab$analysis_complete <- TRUE
  tab$analysis_method <- run$method$method
  tab$organism_supported <- run$method$supported
  invisible(run)
}

# ------------------------------------------------------------------------------------------------
# The outputs.
# ------------------------------------------------------------------------------------------------

#' The method, contrasts and status texts (mod_prot_enrich_registration_helpers.R at 6129b546).
#' @noRd
registerOmicEnrichStatusOutputs <- function(output, input, tab, method) {
  output$analysis_method_display <- shiny::renderText(formatProtEnrichAnalysisMethodText(method()))
  output$contrasts_display <- shiny::renderText(formatProtEnrichContrastsText(tab$contrasts_available))
  output$enrichment_status <- shiny::renderText({
    formatProtEnrichStatusText(analysisComplete = tab$analysis_complete,
      methodInfo = if (tab$analysis_complete) method(), selectedContrast = input$selected_contrast,
      upCutoff = input$up_cutoff, downCutoff = input$down_cutoff, qCutoff = input$q_cutoff,
      gprofilerResults = tab$gprofiler_results, clusterprofilerResults = tab$clusterprofiler_results,
      stringdbResults = tab$stringdb_results)
  })
  invisible(output)
}

#' Each configured service's table and summary, then its plot (the STRING-DB placeholder's plot last), in today's order.
#' The reactive reads are passed as accessors and read inside the render expressions (KIT-064, roar BUG-08).
#' @noRd
registerOmicEnrichResultOutputs <- function(output, input, tab, workflowData, config) {
  services <- config$services
  if ("gprofiler2" %in% services) {
    output$gprofiler_results_table <- renderProtEnrichGprofilerResultsTable(
      gprofilerResults = function() tab$gprofiler_results,
      directionFilter = function() input$gprofiler_direction_filter,
      gprofilerResultsTable = function() tab$gprofiler_results_table)
    output$gprofiler_summary_stats <- shiny::renderText(formatProtEnrichGprofilerSummaryText(
      gprofilerResults = tab$gprofiler_results, directionFilter = input$gprofiler_direction_filter))
  }
  if ("clusterprofiler" %in% services) {
    output$clusterprofiler_results_table <- renderProtEnrichClusterProfilerResultsTable(
      clusterprofilerResults = function() tab$clusterprofiler_results,
      directionFilter = function() input$clusterprofiler_direction_filter,
      clusterprofilerResultsTable = function() tab$clusterprofiler_results_table)
    output$clusterprofiler_summary_stats <- shiny::renderText(formatProtEnrichClusterProfilerSummaryText(
      clusterprofilerResults = tab$clusterprofiler_results, directionFilter = input$clusterprofiler_direction_filter))
  }
  if ("stringdb" %in% services) {
    # KIT-064 (roar BUG-07): the threshold is the tab's declared `q_cutoff` control (default 0.05).
    output$stringdb_results_table <- renderProtEnrichStringDbResultsTable(
      stringdbResults = function() tab$stringdb_results,
      filterSignificant = function() input$stringdb_filter_significant,
      enrichmentPValThresh = function() input$q_cutoff %||% 0.05,
      maxResults = function() input$stringdb_max_results)
    output$stringdb_summary_stats <- shiny::renderText(formatProtEnrichStringDbSummaryText(
      stringdbResults = tab$stringdb_results))
  }
  raw_contrast <- createProtEnrichRawContrastNameReactive(input = input, workflowData = workflowData)
  if ("gprofiler2" %in% services) {
    output$gprofiler_plot <- renderProtEnrichGprofilerPlot(
      analysisComplete = function() tab$analysis_complete,
      enrichmentResultsFull = function() tab$enrichment_results_full,
      rawContrast = raw_contrast,
      directionFilter = function() input$gprofiler_direction_filter)
  }
  if ("clusterprofiler" %in% services) {
    output$clusterprofiler_plot <- renderProtEnrichClusterProfilerPlot(
      analysisComplete = function() tab$analysis_complete,
      enrichmentResultsFull = function() tab$enrichment_results_full,
      rawContrast = raw_contrast,
      directionFilter = function() input$clusterprofiler_direction_filter)
  }
  if ("stringdb" %in% services) output$stringdb_plot <- renderProtEnrichStringDbPlot()
  invisible(output)
}
