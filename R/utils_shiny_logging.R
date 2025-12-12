#' @title Shiny Logger
#' @description Custom logger functions for routing log messages to a Shiny UI.
#'
#' @noRd

# This reactiveVal will hold all the log messages.
# NOTE: In a package context, this global reactiveVal is shared. 
# Ideally this should be session-specific, but preserving original logic.
log_messages <- shiny::reactiveVal("")

#' @title Appender function for Shiny
#' @description This function is a custom log appender for the 'logger' package.
#' It captures log records, formats them, and appends them to a reactiveVal
#' for display in a Shiny UI.
#' @param lines A character vector of the log message.
#' @export
appender_shiny <- function(lines) {
  # Isolate the logic to prevent errors when logging from non-reactive contexts
  shiny::isolate({
    # Get the current messages
    current_logs <- log_messages()
  
    # Collapse the new lines and add a newline character
    new_log <- paste(lines, collapse = "\n")
    
    # Append the new log message
    # We prepend the new message so it appears at the top of the log
    updated_logs <- paste(new_log, current_logs, sep = "\n")
    
    # Trim the log to the last 1000 lines to prevent performance issues
    log_lines <- strsplit(updated_logs, "\n")[[1]]
    if (length(log_lines) > 1000) {
      log_lines <- tail(log_lines, 1000)
      updated_logs <- paste(log_lines, collapse = "\n")
    }
    
    # Update the reactiveVal
    log_messages(updated_logs)
  })
}

#' @title Setup the logger for the Shiny App
#' @description This function configures the 'logger' package to use the custom
#' Shiny appender. It sets the layout to include colors for different log levels.
#' It should be called once when the Shiny app starts.
#' @export
setup_shiny_logger <- function() {
  logger::log_appender(appender_shiny)
  # Use a simple layout that does not include color codes
  logger::log_layout(logger::layout_simple)
  # Use paste formatter to avoid any interpolation issues
  logger::log_formatter(logger::formatter_paste)
  # Set the threshold to INFO to avoid overly verbose logs in production
  logger::log_threshold(logger::INFO)
}

