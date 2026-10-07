# Small source-component check: no package build, omics data or scientific action.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
repo <- normalizePath(args[[1L]], mustWork = TRUE)
destination <- normalizePath(args[[2L]], mustWork = TRUE)
port <- httpuv::randomPort()
app <- callr::r_bg(function(repo, port) {
    source_env <- new.env(parent = globalenv())
    sys.source(file.path(repo, "R/mod_prot_qc_peptide_sample.R"), source_env)
    ui <- shiny::fluidPage(
        shiny::tags$head(shiny::tags$script(shiny::HTML(
            source_env$workflowPeptideSampleAckScript()))),
        shiny::div(`data-testid` = "tile-proteomics", "fixture"),
        shiny::actionButton("apply", "Apply"),
        shiny::verbatimTextOutput("result"), shiny::plotOutput("plot"),
        shiny::textOutput("acknowledged"))
    server <- function(input, output, session) {
        value <- shiny::reactiveVal(NULL)
        output$result <- shiny::renderText({ shiny::req(value()); "accepted result" })
        output$plot <- shiny::renderPlot({
            shiny::req(value())
            grid::grid.draw(grid::textGrob("accepted plot"))
        }, alt = function() { shiny::req(value()); "plot-one" })
        output$acknowledged <- shiny::renderText(input$ack$request_id)
        shiny::observeEvent(input$apply, {
            value(1L)
            session$onFlushed(function() session$sendCustomMessage(
                "multischolar-perf107-ready", list(request_id = "one",
                    content_id = "plot-one", expected_text = "accepted result",
                    plot_identity = "plot-one", text_id = "result",
                    plot_id = "plot", input_id = "ack")), once = TRUE)
        })
    }
    shiny::runApp(shiny::shinyApp(ui, server), host = "127.0.0.1",
        port = port, launch.browser = FALSE)
}, list(repo, port), stdout = file.path(destination, "app.log"),
    stderr = "2>&1", supervise = TRUE)
tryCatch({
    url <- sprintf("http://127.0.0.1:%d", port)
    deadline <- proc.time()[["elapsed"]] + 30
    repeat {
        ready <- tryCatch({ con <- base::url(url, open = "rb"); close(con); TRUE },
            error = function(e) FALSE)
        if (ready) break
        if (!app$is_alive() || proc.time()[["elapsed"]] > deadline) stop("fixture startup failed")
        Sys.sleep(0.1)
    }
    harness <- new.env(parent = globalenv())
    sys.source(file.path(repo, "tools/profiling/workflow_checkpoint_dia_driver.R"), harness)
    browser <- harness$newCheckpointBrowser(.libPaths()[[1L]], destination, managed_url = url)
    tryCatch({
        harness$stabilizeCheckpointManagedBrowser(browser$driver, url)
        session <- browser$driver$get_chromote_session()
        evaluate <- function(js) {
            value <- session$Runtime$evaluate(expression = js,
                returnByValue = TRUE, awaitPromise = TRUE)
            if (length(value$exceptionDetails)) stop("browser JavaScript failed")
            value$result$value
        }
        await <- function(js) {
            deadline <- proc.time()[["elapsed"]] + 20
            while (!isTRUE(evaluate(js))) {
                if (proc.time()[["elapsed"]] > deadline) stop("render acknowledgement timeout")
                Sys.sleep(0.05)
            }
        }
        await("!!(window.Shiny && Shiny.shinyapp && Shiny.shinyapp.$socket && Shiny.shinyapp.$socket.readyState===1 && window.MultiScholaRPerf107Ack)")
        evaluate("document.getElementById('apply').click(); true")
        await("document.getElementById('acknowledged').textContent==='one'")
        state <- evaluate("JSON.stringify(window.MultiScholaRPerf107AckState)")
        writeLines(state, file.path(destination, "render-result.json"))
        cat("real renderPlot acknowledgement passed\n")
    }, finally = browser$chrome$close())
}, finally = app$kill())
