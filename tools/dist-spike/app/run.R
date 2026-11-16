# KIT-130 piece 4: the hello-world Shiny app the Electron shell shows. Started as
#   <bundle>/.../Rscript --vanilla run.R <port>
# It reports, from inside a live session, which R is serving the page.

port <- as.integer(commandArgs(trailingOnly = TRUE)[[1L]])

ui <- shiny::fluidPage(
    shiny::titlePanel("MultiScholaR distribution spike"),
    shiny::p(id = "hello", "hello from the bundled R"),
    shiny::verbatimTextOutput("runtime")
)

server <- function(input, output, session) {
    output$runtime <- shiny::renderText({
        paste(
            c(
                sprintf("r_version=%s", as.character(getRversion())),
                sprintf("r_home=%s", normalizePath(R.home(), winslash = "/")),
                sprintf("lib_paths=%s", paste(normalizePath(.libPaths(), winslash = "/"), collapse = ";")),
                sprintf("r_pid=%d", Sys.getpid())
            ),
            collapse = "\n"
        )
    })
}

shiny::runApp(
    shiny::shinyApp(ui, server),
    host = "127.0.0.1",
    port = port,
    launch.browser = FALSE
)
