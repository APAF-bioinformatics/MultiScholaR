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

#' Run the Shiny Application
#'
#' @param ... arguments to pass to golem_opts. 
#' See `?golem::get_golem_options` for more details.
#' @inheritParams shiny::shinyApp
#'
#' @export
#' @importFrom shiny shinyApp addResourcePath
#' @importFrom golem with_golem_options 
run_app <- function(
  onStart = NULL,
  options = list(), 
  enableBookmarking = NULL,
  uiPattern = "/",
  ...
) {
  # Set up resource path for static files in inst/shiny/www before UI renders
  # Handle both development (pkgload) and installed package scenarios
  shiny_resource_path <- system.file("shiny/www", package = "MultiScholaR")
  if (shiny_resource_path == "" || !dir.exists(shiny_resource_path)) {
    # Development mode: try to find inst/shiny/www relative to common locations
    # Check current directory and parent directories
    possible_paths <- c(
      file.path(getwd(), "inst", "shiny", "www"),
      file.path(dirname(getwd()), "inst", "shiny", "www"),
      file.path(".", "inst", "shiny", "www")
    )
    for (path in possible_paths) {
      if (dir.exists(path)) {
        shiny_resource_path <- normalizePath(path)
        break
      }
    }
  }
  if (dir.exists(shiny_resource_path)) {
    shiny::addResourcePath("www", shiny_resource_path)
  }
  
  # Combine with user-provided onStart if any
  combined_onStart <- function() {
    if (!is.null(onStart)) {
      onStart()
    }
  }
  
  with_golem_options(
    app = shiny::shinyApp(
      ui = app_ui,
      server = app_server,
      onStart = combined_onStart,
      options = options, 
      enableBookmarking = enableBookmarking, 
      uiPattern = uiPattern
    ), 
    golem_opts = list(...)
  )
}

