# Shared fixtures of the driver-core tests (KIT-159 7.3 to 7.5, 7.9, 7.10).

# The package checkout the harness sits in (val-tests.R sets VAL_HARNESS_DIR to the harness folder,
# six levels below the package root).
driverPkgRoot <- function() {
    harness <- Sys.getenv("VAL_HARNESS_DIR")
    if (!nzchar(harness)) return(NA_character_)
    normalizePath(file.path(harness, "..", "..", "..", "..", "..", ".."), mustWork = FALSE)
}

# The ten routes of section 7.4 (the GUI-entry lane ids).
driverRoutes <- function() {
    c("prot_dia", "prot_dia_limpa", "prot_lfq", "prot_lfq_fragpipe", "prot_tmt", "trx_screen",
      "lipid_canonical", "metab_lc", "metab_gc", "metab_combined")
}

# A toy module for the driver: a select the server updates with choices only, a modal with a
# primary confirming button, a re-rendering uiOutput, a refusing action, a diagnostic action that
# refuses, an erroring action and a tab status. Its state is `workflow_data` (head, history, tabs).
toyUi <- function(id) {
    ns <- shiny::NS(id)
    htmltools::tagList(
        shiny::selectInput(ns("method"), "Method", c("median", "cyclicloess")),
        shiny::selectInput(ns("mode"), "Mode", c("narrow", "wide")),
        shiny::textInput(ns("label"), "Label", "x"),
        shiny::uiOutput(ns("dynamic")),
        shiny::actionButton(ns("load"), "Load"),
        shiny::actionButton(ns("run"), "Run", class = "btn-primary")
    )
}

toyServer <- function(id) {
    shiny::moduleServer(id, function(input, output, session) {
        workflow_data <- shiny::reactiveValues(head = "initial", history = character(),
                                               tab_status = list(da = "pending"))
        output$dynamic <- shiny::renderUI({
            k <- if (identical(input$mode, "wide")) 9 else 3
            shiny::numericInput(session$ns("k"), "k", value = k)
        })
        advance <- function(state) {
            workflow_data$head <- state
            workflow_data$history <- c(shiny::isolate(workflow_data$history), state)
        }
        shiny::observeEvent(input$load, {
            shiny::updateSelectInput(session, "method", choices = c("cyclicloess", "median"))
            shiny::showModal(shiny::modalDialog(
                shiny::checkboxInput(session$ns("filter_org"), "Filter", FALSE),
                footer = htmltools::tagList(
                    shiny::modalButton("Cancel"),
                    shiny::actionButton(session$ns("confirm"), "Confirm Selection",
                                        class = "btn-primary"))))
        })
        shiny::observeEvent(input$confirm, {
            shiny::removeModal()
            advance("loaded")
        })
        shiny::observeEvent(input$run, {
            if (!identical(input$method, "cyclicloess")) {
                shiny::showNotification("toy_method_refused: wrong method", type = "error")
                return()
            }
            advance("normalised")
        })
        shiny::observeEvent(input$refuse, {
            shiny::showNotification("refused 12 rows at 0.5", type = "error")
        })
        shiny::observeEvent(input$boom, stop("toy failure"))
        shiny::observeEvent(input$finish, {
            workflow_data$tab_status <- list(da = "complete")
        })
    })
}

toySpec <- function() list(server = toyServer, id = "toy", tabs = "toy_tabs", args = list())

toyObserve <- function(env) {
    wd <- get("workflow_data", envir = env)
    shiny::isolate(list(head = wd$head, history = wd$history,
                        tab_status = lapply(wd$tab_status, as.character)))
}

toyScript <- function(steps) list(route = "toy", omic = "toy", steps = steps)

# The standard tables' columns (section 7.3).
DA_COLUMNS <- c("feature_id", "assay", "contrast", "log2fc", "p", "s_own", "call_statistic",
                "called_own", "mean_expression")
