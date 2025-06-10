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
  
  # Try to find the app directory
  app_path <- NULL
  
  # First, check if we're in development mode (package loaded with devtools::load_all)
  # This checks if the package is loaded from source
  pkg_path <- find.package("MultiScholaR", quiet = TRUE)
  
  if (length(pkg_path) > 0) {
    # Check for development location (R/shiny)
    dev_app_path <- file.path(pkg_path, "R", "shiny")
    if (dir.exists(dev_app_path) && file.exists(file.path(dev_app_path, "app.R"))) {
      app_path <- dev_app_path
      log_info("Running app from development location: {app_path}")
      
      # Source necessary files for development
      shiny_applets_path <- file.path(pkg_path, "R", "shiny_applets.R")
      if (file.exists(shiny_applets_path)) {
        source(shiny_applets_path)
        log_info("Loaded shiny_applets.R")
      }
      
      # Make sure other key functions are available
      helper_functions_path <- file.path(pkg_path, "R", "helper_functions.R")
      if (file.exists(helper_functions_path) && !exists("loadDependencies")) {
        source(helper_functions_path)
        log_info("Loaded helper_functions.R")
      }
      
      file_management_path <- file.path(pkg_path, "R", "file_management.R")
      if (file.exists(file_management_path) && !exists("setupDirectories")) {
        source(file_management_path)
        log_info("Loaded file_management.R")
      }
    }
    
    # If not found in development location, check inst/shiny (for installed package)
    if (is.null(app_path)) {
      inst_app_path <- system.file("shiny", package = "MultiScholaR")
      if (dir.exists(inst_app_path) && file.exists(file.path(inst_app_path, "app.R"))) {
        app_path <- inst_app_path
        log_info("Running app from installed package location: {app_path}")
      }
    }
  }
  
  # If still not found, error
  if (is.null(app_path) || !dir.exists(app_path)) {
    stop("Could not find the MultiScholaR Shiny app. ",
         "Make sure the package is properly installed or you're in the package directory.")
  }
  
  # Ensure key functions are available in the global environment for the app
  if (exists("loadDependencies", mode = "function")) {
    assign("loadDependencies", loadDependencies, envir = .GlobalEnv)
  } else {
    log_warn("loadDependencies function not found")
  }
  
  if (exists("setupDirectories", mode = "function")) {
    assign("setupDirectories", setupDirectories, envir = .GlobalEnv)
  } else {
    log_warn("setupDirectories function not found")
  }
  
  if (exists("RunApplet", mode = "function")) {
    assign("RunApplet", RunApplet, envir = .GlobalEnv)
  } else {
    log_warn("RunApplet function not found")
  }
  
  # Load dependencies
  if (exists("loadDependencies", mode = "function")) {
    message("Loading dependencies...")
    tryCatch({
      loadDependencies(verbose = FALSE)
    }, error = function(e) {
      log_warn("Some dependencies may not have loaded properly: {e$message}")
    })
  }
  
  # Launch the app
  message("Launching MultiScholaR Workflow Application...")
  
  # Use a separate R process if clean_env is TRUE to avoid all conflicts
  if (clean_env) {
    shiny::runApp(app_path, launch.browser = launch.browser, ...)
  } else {
    shiny::runApp(app_path, launch.browser = launch.browser, ...)
  }
} 