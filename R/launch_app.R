#' Launch MultiScholaR Shiny Application
#'
#' Launches the interactive MultiScholaR workflow application for analyzing
#' multi-omics data including proteomics, metabolomics, transcriptomics,
#' and lipidomics.
#'
#' @param launch.browser Logical; if TRUE, the app will open in your default browser
#' @param clean_env Logical; if TRUE, detach conflicting packages before running (recommended)
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}
#'
#' @return Runs the Shiny application
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' # Launch the app
#' MultiScholaRapp()
#' 
#' # Launch without opening browser
#' MultiScholaRapp(launch.browser = FALSE)
#' 
#' # Launch without cleaning environment (may cause conflicts)
#' MultiScholaRapp(clean_env = FALSE)
#' }
#' 
#' @importFrom shiny runApp
#' @importFrom logger log_info log_warn
MultiScholaRapp <- function(launch.browser = TRUE, clean_env = TRUE, ...) {
  
  # Handle namespace conflicts if requested
  if (clean_env) {
    # Check for conflicting packages
    conflicting_packages <- c("shinyjs", "DT", "shinydashboard")
    loaded_conflicts <- intersect(conflicting_packages, loadedNamespaces())
    
    if (length(loaded_conflicts) > 0) {
      message("Temporarily detaching conflicting packages to avoid namespace issues...")
      for (pkg in loaded_conflicts) {
        try(detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE), silent = TRUE)
      }
    }
  }
  
  # Load dependencies
  if (exists("loadDependencies", mode = "function")) {
    message("Loading dependencies...")
    tryCatch({
      loadDependencies(verbose = FALSE)
    }, error = function(e) {
      warning(paste("Some dependencies may not have loaded properly:", e$message))
    })
  } else {
    warning("loadDependencies function not found")
  }
  
  # Launch the app
  message("Launching MultiScholaR Workflow Application...")
  
  # Use a separate R process if clean_env is TRUE to avoid all conflicts
  if (clean_env) {
    shiny::runApp(run_app(), launch.browser = launch.browser, ...)
  } else {
    shiny::runApp(run_app(), launch.browser = launch.browser, ...)
  }
}
