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

# ============================================================================
# func_multiomics_enrich_stringdb_api.R
# ============================================================================
# Purpose: STRING-DB submission, download, and retrieval helpers
# 
# This file owns the external STRING-DB request and result-transfer boundary
# shared by multiomics enrichment workflows.
#
# Dependencies:
# - gprofiler2, clusterProfiler, STRINGdb
# - func_general_plotting_*_helpers.R (for visualization)
# - func_general_enrichment_* (for shared enrichment functions)
# ============================================================================

submitStringDBEnrichment <- function(input_data_frame,
                                     identifier_column_name,
                                     value_column_name,
                                     caller_identity,
                                     api_key,
                                     species = "9606",
                                     ge_fdr = 0.05,
                                     ge_enrichment_rank_direction = -1) {

  # Load required packages, install if missing
  if (!requireNamespace("pacman", quietly = TRUE)) {
    install.packages("pacman")
  }
  pacman::p_load(
    char = c(
      "httr",      # For HTTP requests
      "jsonlite",  # For JSON parsing
      "dplyr",     # For data manipulation
      "readr",     # For reading TSV/CSV files
      "checkmate"  # For argument checking
    ),
    install = TRUE,
    update = FALSE
  )

  # --- Input Validation & Data Preparation ---
  checkmate::assertDataFrame(input_data_frame, min.rows = 1, .var.name = "input_data_frame")
  checkmate::assertString(identifier_column_name, min.chars = 1, .var.name = "identifier_column_name")
  checkmate::assertString(value_column_name, min.chars = 1, .var.name = "value_column_name")
  checkmate::assertChoice(identifier_column_name, choices = names(input_data_frame), .var.name = "identifier_column_name")
  checkmate::assertChoice(value_column_name, choices = names(input_data_frame), .var.name = "value_column_name")
  checkmate::assertString(caller_identity, min.chars = 1, .var.name = "caller_identity")
  checkmate::assertString(api_key, min.chars = 1, .var.name = "api_key") # Validate api_key
  # Convert numeric species to string if needed, then validate
  if (is.numeric(species)) {
    species <- as.character(species)
  }
  checkmate::assertString(species, .var.name = "species")
  checkmate::assertNumber(ge_fdr, lower = 0, upper = 1, .var.name = "ge_fdr")
  checkmate::assertChoice(as.integer(ge_enrichment_rank_direction), c(-1, 0, 1), .var.name = "ge_enrichment_rank_direction")

  if (identifier_column_name == value_column_name) {
    stop("Identifier and value column names must be different.")
  }

  ids_vector <- input_data_frame[[identifier_column_name]]
  values_vector <- input_data_frame[[value_column_name]]

  checkmate::assert(
    checkmate::checkCharacter(ids_vector, any.missing = TRUE, min.len = 1),
    checkmate::checkFactor(ids_vector, any.missing = TRUE, min.len = 1),
    .var.name = paste0("Identifier column '", identifier_column_name, "'")
  )
  checkmate::assertNumeric(values_vector, any.missing = TRUE, min.len = 1,
                           .var.name = paste0("Value column '", value_column_name, "'"))

  temp_df <- data.frame(
    ids = as.character(ids_vector),
    vals = values_vector,
    stringsAsFactors = FALSE
  )

  initial_rows <- nrow(temp_df)
  temp_df <- temp_df[!is.na(temp_df$ids) & !is.na(temp_df$vals), ]

  removed_rows_count <- initial_rows - nrow(temp_df)
  if (removed_rows_count > 0) {
    message(paste(
      removed_rows_count,
      "row(s) were removed from the input data due to NA values in the identifier or value columns."
    ))
  }

  if (nrow(temp_df) == 0) {
    stop("No valid identifier/value pairs remaining after NA removal. Cannot submit to STRING API.")
  }

  identifiers_string <- temp_df |>
    dplyr::mutate(combined_string = paste(.data$ids, .data$vals, sep = "\t")) |>
    dplyr::pull(.data$combined_string) |>
    paste(collapse = "\n")

  # --- API Configuration ---
  STRING_API_URL <- "https://version-12-0.string-db.org/api"
  OUTPUT_FORMAT  <- "json"
  METHOD_SUBMIT  <- "valuesranks_enrichment_submit"

  request_url <- paste(STRING_API_URL, OUTPUT_FORMAT, METHOD_SUBMIT, sep = "/")

  # --- Prepare Parameters for POST Request ---
  params_list <- list(
    species = species,
    caller_identity = caller_identity,
    identifiers = identifiers_string,
    api_key = api_key, # api_key is used here
    ge_fdr = ge_fdr,
    ge_enrichment_rank_direction = as.integer(ge_enrichment_rank_direction)
  )

  # --- Call STRING API ---
  response <- tryCatch({
    httr::POST(url = request_url, body = params_list, encode = "form")
  }, error = function(e) {
    message(paste("HTTP POST request failed:", e$message)) # More direct message
    # Return a list indicating failure, including the api_key for consistency
    return(
      list(
        job_id = NULL,
        api_key = api_key,
        submission_response = list(status = "error", message = paste("HTTP POST request failed:", e$message))
      )
    )
  })

  # --- Process Response ---
  response_content_text <- httr::content(response, "text", encoding = "UTF-8")

  if (httr::http_error(response)) {
    message( # More direct message
      paste0(
        "STRING API request failed with HTTP status: ", httr::status_code(response),
        "\nResponse content:\n", response_content_text
      )
    )
    return(
      list(
        job_id = NULL,
        api_key = api_key,
        submission_response = list(
          status = "error",
          message = paste("STRING API request failed with HTTP status:", httr::status_code(response)),
          details = response_content_text
        )
      )
    )
  }

  parsed_json_response <- tryCatch({
    temp_parsed <- jsonlite::fromJSON(response_content_text,
                                      simplifyDataFrame = FALSE,
                                      simplifyVector = FALSE,
                                      simplifyMatrix = FALSE)
    if (is.list(temp_parsed) && length(temp_parsed) >= 1 && is.list(temp_parsed[[1]])) {
      temp_parsed[[1]]
    } else if (is.list(temp_parsed) && !is.null(names(temp_parsed))) {
      temp_parsed
    } else {
      stop("Unexpected JSON structure from API.") # Simpler stop message
    }
  }, error = function(e) {
    message( # More direct message
      paste0(
        "Failed to parse JSON response from STRING API.",
        "\nOriginal error: ", e$message,
        "\nResponse content:\n", response_content_text
      )
    )
    return(
      list(
        job_id = NULL,
        api_key = api_key,
        submission_response = list(
          status = "error",
          message = paste("Failed to parse JSON response:", e$message),
          details = response_content_text
        )
      )
    )
  })

  # --- Output Results & Prepare Return List ---
  current_job_id <- NULL # Initialize to NULL

  if (!is.null(parsed_json_response$status) && parsed_json_response$status == "error") {
    message(paste("STRING API Error - Status:", parsed_json_response$status))
    if (!is.null(parsed_json_response$message)) {
      message(paste("Message:", parsed_json_response$message))
    }
  } else if (!is.null(parsed_json_response$job_id)) {
    current_job_id <- parsed_json_response$job_id # Assign job_id
    message(paste("Job submitted successfully to STRING API. Job ID:", current_job_id))
  } else {
    message("Warning: Unexpected API response structure. Could not find 'status' or 'job_id'.")
    message("Raw response content was:")
    message(response_content_text)
  }

  return(
    list(
      job_id = current_job_id,
      api_key = api_key, # Return the api_key used
      submission_response = parsed_json_response
    )
  )
}


# ----------------------------------------------------------------------------
# downloadStringDBGraph
# ----------------------------------------------------------------------------
#' Download STRING DB Graph Image
#'
#' @description
#' Downloads a graph image from a given URL provided by the STRING API.
#'
#' @param graph_url Character string: The direct URL to the graph image.
#'
#' @return The raw binary content of the graph image if successful (can be written
#'         to a file using `writeBin`). Returns `NULL` if the download fails.
#'         Messages about the download status are printed to the console.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # This function is typically called by retrieveStringDBEnrichmentResults
#' # example_graph_url <- "https://version-12-0.string-db.org/api/image/enrichmentfigure?networkId=bNEXfEymvDsZ"
#' # image_content <- downloadStringDBGraph(graph_url = example_graph_url)
#' # if (!is.null(image_content)) {
#' #   # Save the image to a file (e.g., as PNG if that's the format)
#' #   writeBin(image_content, "string_enrichment_graph.png")
#' #   message("Graph image downloaded and saved as string_enrichment_graph.png")
#' # }
#' }
downloadStringDBGraph <- function(graph_url) {
  checkmate::assertString(graph_url, min.chars = 1, .var.name = "graph_url")

  message(paste("Attempting to download graph image from:", graph_url))

  graph_response_http <- tryCatch({
    httr::GET(url = graph_url)
  }, error = function(e) {
    message(paste("HTTP GET request for graph image download failed:", e$message))
    return(NULL)
  })

  if (is.null(graph_response_http)) {
    message("Failed to initiate download of graph image.")
    return(NULL)
  }

  if (httr::http_error(graph_response_http)) {
    message(
      paste0(
        "STRING API returned an HTTP error during graph image download: ",
        httr::status_code(graph_response_http)
        # Avoid printing content for binary files directly unless debugging
      )
    )
    return(NULL)
  }

  # Get raw content for images
  graph_content_raw <- httr::content(graph_response_http, "raw")

  if (length(graph_content_raw) > 0) {
    message("Graph image successfully downloaded.")
    return(graph_content_raw)
  } else {
    message("Downloaded graph image content is empty.")
    return(NULL)
  }
}


# ----------------------------------------------------------------------------
# downloadStringDBResultsFile
# ----------------------------------------------------------------------------
#' Download STRING DB Results File
#'
#' @description
#' Downloads a results file (typically TSV) from a given URL provided by the STRING API
#' and parses it into an R data frame.
#'
#' @param download_url Character string: The direct URL to the results file.
#'
#' @return A data frame containing the parsed results from the URL.
#'         Returns `NULL` if the download or parsing fails.
#'         Messages about the download status are printed to the console.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # This function is typically called by retrieveStringDBEnrichmentResults
#' # but can be used standalone if you have a direct download URL.
#' # example_download_url <- "https://version-12-0.string-db.org/api/tsv/downloadenrichmentresults?networkId=bNEXfEymvDsZ"
#' # results_df <- downloadStringDBResultsFile(download_url = example_download_url)
#' # if (!is.null(results_df)) {
#' #   print(head(results_df))
#' # }
#' }
downloadStringDBResultsFile <- function(download_url) {
  checkmate::assertString(download_url, min.chars = 1, .var.name = "download_url")

  message(paste("Attempting to download results from:", download_url))

  results_response_http <- tryCatch({
    httr::GET(url = download_url)
  }, error = function(e) {
    message(paste("HTTP GET request for results download failed:", e$message))
    return(NULL)
  })

  if (is.null(results_response_http)) {
    message("Failed to initiate download of results.")
    return(NULL)
  }

  if (httr::http_error(results_response_http)) {
    message(
      paste0(
        "STRING API returned an HTTP error during results download: ",
        httr::status_code(results_response_http),
        "\nResponse content: ", httr::content(results_response_http, "text", encoding = "UTF-8")
      )
    )
    return(NULL)
  }

  results_content_text <- httr::content(results_response_http, "text", encoding = "UTF-8")

  enrichment_df <- tryCatch({
    readr::read_tsv(results_content_text, show_col_types = FALSE)
  }, error = function(e) {
    message(paste("Failed to parse TSV results:", e$message))
    message(paste("Raw TSV content (first 1000 chars):\n", substr(results_content_text, 1, 1000), "..."))
    return(NULL)
  })

  if (!is.null(enrichment_df)) {
    message("Enrichment results successfully downloaded and parsed.")
    return(enrichment_df)
  } else {
    message("Failed to create data frame from downloaded results.")
    return(NULL)
  }
}


# ----------------------------------------------------------------------------
# retrieveStringDBEnrichmentResults
# ----------------------------------------------------------------------------
#' Retrieve STRING DB Enrichment Results
#'
#' @description
#' Polls the STRING API for the status of a submitted enrichment job.
#' Upon successful completion, it obtains the download URL and then uses
#' `downloadStringDBResultsFile` to fetch and parse the results into a data frame.
#'
#' @param submission_info A list object returned by `submitStringDBEnrichment`.
#'                        This list must contain `job_id` (non-NULL) and `api_key`.
#' @param polling_interval_seconds Numeric: The number of seconds to wait between
#'                                 polling attempts for job status. Default is 10.
#' @param max_polling_attempts Numeric: The maximum number of polling attempts before
#'                             timing out. Default is 30.
#'
#' @return A list containing the following elements if the job is successful:
#'         - `enrichment_data`: A data frame of enrichment results (or `NULL` on download failure).
#'         - `page_url`: Character string, URL to the STRING results page (or `NULL` if not found).
#'         - `graph_url`: Character string, URL for the enrichment graph image (or `NULL` if not found).
#'         - `graph_image_content`: Raw vector, the binary content of the graph image (or `NULL` on download failure).
#'         - `status_details`: A data frame or list with the full status information from the last successful poll.
#'         Returns `NULL` if the job polling ultimately fails, times out, or a critical error occurs.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # This function is typically used after submitStringDBEnrichment.
#' # Assume 'submission_info_output' is the actual list returned by
#' # a successful call to submitStringDBEnrichment.
#'
#' # Example structure of submission_info_output:
#' # submission_info_output <- list(
#' #   job_id = "b3R6oioiQSRO",      # Actual job ID from a submission
#' #   api_key = "YOUR_API_KEY",   # The API key used
#' #   submission_response = list(
#' #     job_id = "b3R6oioiQSRO",
#' #     status = "submitted",
#' #     message = "Job was successfully submitted to the queue!"
#' #   )
#' # )
#'
#' # Make sure job_id is present and the submission was initially successful
#' # before calling retrieveStringDBEnrichmentResults.
#' # if (!is.null(submission_info_output$job_id) &&
#' #     !is.null(submission_info_output$submission_response$status) &&
#' #     submission_info_output$submission_response$status == "submitted") {
#' #
#' #   results_list <- retrieveStringDBEnrichmentResults(
#' #     submission_info = submission_info_output, # Pass the list directly
#' #     polling_interval_seconds = 5,
#' #     max_polling_attempts = 6
#' #   )
#' #   if (!is.null(results_list)) {
#' #     if(!is.null(results_list$enrichment_data)) {
#' #        print(head(results_list$enrichment_data))
#' #     }
#' #     if(!is.null(results_list$graph_image_content)) {
#' #        writeBin(results_list$graph_image_content, "enrichment_graph.png")
#' #        message("Graph saved to enrichment_graph.png")
#' #     }
#' #     message(paste("Page URL:", results_list$page_url))
#' #   } else {
#' #     message("Could not retrieve enrichment results package.")
#' #   }
#' # } else {
#' # message("Submission was not successful or job_id missing in submission_info_output.")
#' # print(submission_info_output) # For debugging
#' # }
#' }
retrieveStringDBEnrichmentResults <- function(submission_info,
                                              polling_interval_seconds = 10,
                                              max_polling_attempts = 30) {

  # Ensure necessary packages are loaded
  if (!all(sapply(c("httr", "jsonlite", "readr", "checkmate"), requireNamespace, quietly = TRUE))) {
    stop("One or more required packages (httr, jsonlite, readr, checkmate) are not installed/loaded.")
  }

  # --- Input Validation ---
  checkmate::assertList(submission_info, min.len = 2, .var.name = "submission_info")
  checkmate::assertString(submission_info$job_id, min.chars = 1, .var.name = "submission_info$job_id")
  checkmate::assertString(submission_info$api_key, min.chars = 1, .var.name = "submission_info$api_key")
  checkmate::assertNumber(polling_interval_seconds, lower = 1, .var.name = "polling_interval_seconds")
  checkmate::assertNumber(max_polling_attempts, lower = 1, .var.name = "max_polling_attempts")

  # --- API Configuration ---
  STRING_API_URL_BASE <- "https://version-12-0.string-db.org/api"
  STATUS_METHOD       <- "json/valuesranks_enrichment_status"

  status_request_url <- paste0(
    STRING_API_URL_BASE, "/", STATUS_METHOD,
    "?api_key=", utils::URLencode(submission_info$api_key, reserved = TRUE),
    "&job_id=", utils::URLencode(submission_info$job_id, reserved = TRUE)
  )

  page_url_from_api <- NULL
  download_url_from_api <- NULL
  graph_url_from_api <- NULL
  last_successful_status_data <- NULL
  job_final_status <- "polling"

  # --- Polling Loop for Job Status ---
  message(paste0("Polling STRING API for job status (Job ID: ", submission_info$job_id, "). Will attempt up to ", max_polling_attempts, " times."))
  for (attempt in 1:max_polling_attempts) {
    message(paste0("Attempt ", attempt, " of ", max_polling_attempts, "..."))

    status_response_http <- tryCatch({
      httr::GET(url = status_request_url)
    }, error = function(e) {
      message(paste("HTTP GET request for status failed on attempt", attempt, ":", e$message))
      return(NULL)
    })

    if (is.null(status_response_http)) {
      if (attempt < max_polling_attempts) {
        Sys.sleep(polling_interval_seconds)
        next
      } else {
        job_final_status <- "error_http"
        break
      }
    }

    status_content_text <- httr::content(status_response_http, "text", encoding = "UTF-8")

    if (httr::http_error(status_response_http)) {
      message(
        paste0(
          "STRING API returned an HTTP error on attempt ", attempt,
          ": ", httr::status_code(status_response_http),
          "\nResponse content: ", status_content_text
        )
      )
      if (attempt < max_polling_attempts) {
        Sys.sleep(polling_interval_seconds)
        next
      } else {
        job_final_status <- "error_api_http"
        break
      }
    }

    current_status_data_parsed <- tryCatch({ # Renamed to avoid conflict
      jsonlite::fromJSON(status_content_text, simplifyDataFrame = TRUE)
    }, error = function(e) {
      message(paste("Failed to parse JSON status response on attempt", attempt, ":", e$message))
      message(paste("Raw response was:", status_content_text))
      return(NULL)
    })

    if (is.null(current_status_data_parsed) || !is.data.frame(current_status_data_parsed) || nrow(current_status_data_parsed) == 0) {
      message("Parsed status data is not in the expected format (1-row data frame).")
      if (attempt < max_polling_attempts) {
        Sys.sleep(polling_interval_seconds)
        next
      } else {
        job_final_status <- "error_parsing"
        break
      }
    }

    last_successful_status_data <- current_status_data_parsed[1, ] # Store the first row (should be only one)
    current_status  <- last_successful_status_data$status
    current_message <- last_successful_status_data$message

    message(paste0("Job status: '", current_status, "'. Message: '", current_message, "'"))

    if (current_status == "success") {
      # Extract all URLs if present
      page_url_from_api     <- last_successful_status_data$page_url
      download_url_from_api <- last_successful_status_data$download_url
      graph_url_from_api    <- last_successful_status_data$graph_url

      if (is.null(download_url_from_api) || !nzchar(download_url_from_api)) {
        message("Job status is 'success', but essential download_url is missing or empty.")
        job_final_status <- "error_missing_download_url"
      } else {
        message("Job finished successfully. All relevant URLs found (or attempted to find).")
        job_final_status <- "success"
      }
      break # Exit polling loop
    } else if (current_status == "error") {
      message(paste("Job failed with error from API:", current_message))
      job_final_status <- "error_api_reported"
      break
    } else if (current_status %in% c("submitted", "queued", "running")) {
      if (attempt == max_polling_attempts) {
        message("Maximum polling attempts reached, and job is still processing.")
        job_final_status <- "timeout_processing"
        break
      }
      Sys.sleep(polling_interval_seconds)
    } else {
      message(paste("Unknown job status received:", current_status))
      if (attempt == max_polling_attempts) {
        job_final_status <- "timeout_unknown_status"
        break
      }
      Sys.sleep(polling_interval_seconds)
    }
  }

  # --- Prepare Results Package ---
  if (job_final_status == "success") {
    enrichment_df_results <- NULL
    graph_image_content_results <- NULL

    if (!is.null(download_url_from_api) && nzchar(download_url_from_api)) {
      enrichment_df_results <- downloadStringDBResultsFile(download_url = download_url_from_api)
    } else {
      message("Download URL for results was not available, skipping results table download.")
    }

    if (!is.null(graph_url_from_api) && nzchar(graph_url_from_api)) {
      graph_image_content_results <- downloadStringDBGraph(graph_url = graph_url_from_api)
    } else {
      message("Graph URL was not available, skipping graph image download.")
    }

    return(
      list(
        enrichment_data = enrichment_df_results,
        page_url = page_url_from_api, # Will be NULL if not found in API response
        graph_url = graph_url_from_api, # Will be NULL if not found
        download_url = download_url_from_api, # Will be NULL if not found
        graph_image_content = graph_image_content_results,
        status_details = last_successful_status_data
      )
    )
  } else {
    message(paste("Could not retrieve a full results package. Final job status/outcome:", job_final_status))
    if (!is.null(last_successful_status_data)) { # Changed from `exists("status_data")`
      message("Further details from last status check:")
      print(last_successful_status_data)
    }
    return(NULL) # Return NULL if polling failed or critical URLs were missing for "success"
  }
}


# ----------------------------------------------------------------------------
