# ----------------------------------------------------------------------------
# saveTimeRecord
# ----------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#' @title Save Time Record
#' @description Save timing information for compute steps
#' @param compute_time_record Data frame to store timing records
#' @param step_name Name of the step being timed
#' @param toc_output Output from tictoc::toc()
#' @importFrom grid convertX convertY unit
#' @export
saveTimeRecord <- function ( compute_time_record, step_name, toc_output ) {

  if(length(which(compute_time_record[,"step_name"] == "")) == 0) {
    stop("saveTimeRecord: No more room on compute time record table.")
  }
  # find_first_empty_slot
  first_empty_slot <- which(compute_time_record[,"step_name"] == "")[1]

  print(first_empty_slot)

  compute_time_record[first_empty_slot, "step_name"] <- step_name
  compute_time_record[first_empty_slot, "time_elapsed"] <- toc_output$callback_msg

  return(compute_time_record)
}

# ----------------------------------------------------------------------------
# getFunctionName
# ----------------------------------------------------------------------------
#' @export
getFunctionName <- function() {
  calls <- sys.calls()
  current_call <- calls[[length(calls) - 1]]
  as.character(current_call[1])
}

# ----------------------------------------------------------------------------
# getFunctionNameSecondLevel
# ----------------------------------------------------------------------------
#' @export
getFunctionNameSecondLevel <- function() {
  calls <- sys.calls()
  current_call <- calls[[length(calls) - 2]]
  as.character(current_call[1])
}

#' Get Actual Process Memory (Task Manager View)
#' 
#' Returns the actual memory used by the R process as reported by the OS.
#' This is what you see in Task Manager / Activity Monitor, and includes:
#' - R heap memory
#' - Environment captures in ggplot/tidyverse objects
#' - Lazy evaluation promises
#' - Memory-mapped files
#' 
#' @return Memory usage in MB (numeric)
#' @export
getProcessMemoryMB <- function() {
  tryCatch({
    if (.Platform$OS.type == "windows") {
      # Windows: Use tasklist command (more reliable than wmic on modern Windows)
      pid <- Sys.getpid()
      cmd_output <- system(sprintf('tasklist /FI "PID eq %d" /FO CSV /NH', pid), intern = TRUE)
      # Parse CSV output: "process.exe","PID","Session","Session#","Mem Usage"
      # Mem Usage format: "1,234,567 K"
      if (length(cmd_output) > 0 && nchar(cmd_output[1]) > 0) {
        mem_str <- gsub('"', '', strsplit(cmd_output[1], ',')[[1]][5])
        mem_kb <- as.numeric(gsub('[^0-9]', '', mem_str))
        return(mem_kb / 1024)  # Convert KB to MB
      }
    } else {
      # Unix/Mac: Use ps command
      pid <- Sys.getpid()
      mem_kb <- as.numeric(system(sprintf("ps -o rss= -p %d", pid), intern = TRUE))
      return(mem_kb / 1024)  # Convert KB to MB
    }
    return(NA_real_)
  }, error = function(e) {
    return(NA_real_)
  })
}

#' Get R Heap Memory
#' 
#' Returns R's internal heap memory (what pryr::mem_used() shows).
#' NOTE: This does NOT include environment captures!
#' 
#' @return Memory usage in MB (numeric)
#' @export
getRHeapMemoryMB <- function() {
  sum(gc()[, 2])
}

# Helper function for null-coalescing
#' @export
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' environment capture issues where R heap stays flat but process memory spikes.
#' 
#' @param step_name Character string describing the current step
#' @param context Optional context string (e.g., function name)
#' @param log_level One of "DEBUG66", "INFO", "WARN"
#' @return Invisibly returns a list with r_heap_mb and process_mb
#' @export
checkMemoryBoth <- function(step_name, context = "", log_level = "DEBUG66") {
  # Disabled for performance to avoid slow system calls
  invisible(list(r_heap_mb = 0, process_mb = 0, hidden_mb = 0))
}

#' Memory Delta Reporter
#' 
#' Compares current memory to a baseline and reports the change.
#' Useful for bracketing operations to see exactly how much memory they consume.
#' 
#' @param baseline_mem List from a previous checkMemoryBoth() call
#' @param step_name Character string describing the completed step
#' @param context Optional context string
#' @return Invisibly returns the delta values
#' @export
reportMemoryDelta <- function(baseline_mem, step_name, context = "") {
  # Disabled for performance to avoid slow system calls
  invisible(list(delta_r_heap = 0, delta_process = 0))
}

