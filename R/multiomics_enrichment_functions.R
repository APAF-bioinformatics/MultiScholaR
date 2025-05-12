#' Submit STRING DB Values/Ranks Enrichment Analysis
#'
#' @description
#' Submits a job to the STRING API (version 12.0) for values/ranks enrichment analysis.
#' This function reads protein identifiers from an input file and sends them
#' to the STRING database for analysis.
#'
#' @param input_data_frame A data frame containing at least two columns: one for
#'                        identifiers and one for associated numerical values.
#' @param identifier_column_name Character string: The name of the column in
#'                               `input_data_frame` that contains the protein/gene
#'                               identifiers (e.g., Ensembl IDs, gene symbols).
#' @param value_column_name Character string: The name of the column in
#'                          `input_data_frame` that contains the numerical values
#'                          (e.g., log fold change, p-value, score) associated
#'                          with each identifier. This column must be numeric.
#' @param caller_identity Character string: An identifier for your script or application
#'                        (e.g., "my_research_project_R_script").
#' @param api_key Character string: Your personal STRING API key.
#' @param species Numeric: NCBI/STRING species identifier. Default is 9606 (Homo sapiens).
#' @param ge_fdr Numeric: FDR threshold for gene expression enrichment. Default is 0.05.
#' @param ge_enrichment_rank_direction Integer: Direction for enrichment rank.
#'                                       (-1, 0, or 1). Default is -1.
#'
#' @return A list containing:
#'         - `job_id`: The job ID if submission was successful, otherwise `NULL`.
#'         - `api_key`: The API key used for the submission.
#'         - `submission_response`: The full parsed JSON response from the API.
#'         Messages about the submission status are also printed to the console.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a dummy data frame
#' example_data <- data.frame(
#'   protein_id = c("TP53", "EGFR", "BRCA1", "TNF", "IL6", "MYC", "MISSING_ID", NA_character_),
#'   expression_value = c(-0.585, 0.388, -0.079, 1.2, -2.1, NA_real_, 0.99, 0.5),
#'   other_info = letters[1:8],
#'   stringsAsFactors = FALSE
#' )
#'
#' # IMPORTANT: Replace "YOUR_API_KEY" with your actual STRING API key
#' # IMPORTANT: Replace "your_application_name" with a meaningful caller identity
#' submission_info <- submitStringDBEnrichment(
#'   input_data_frame = example_data,
#'   identifier_column_name = "protein_id",
#'   value_column_name = "expression_value",
#'   caller_identity = "my_R_enrichment_script_v3", # Updated example version
#'   api_key = "YOUR_API_KEY", # !!! REPLACE THIS !!!
#'   species = 9606, # Human
#'   ge_fdr = 0.05,
#'   ge_enrichment_rank_direction = -1
#' )
#'
#' # Check the submission result
#' # Access status and message from the submission_response element
#' if (!is.null(submission_info$job_id) &&
#'     !is.null(submission_info$submission_response$status) &&
#'     submission_info$submission_response$status == "submitted") {
#'   message(paste("Job submitted. Job ID:", submission_info$job_id))
#'
#'   # Attempt to retrieve results using the submission_info object directly
#'   enrichment_data <- retrieveStringDBEnrichmentResults(
#'     submission_info = submission_info, # Pass the whole list
#'     polling_interval_seconds = 10,
#'     max_polling_attempts = 12
#'   )
#'
#'   if (!is.null(enrichment_data)) {
#'     message("Enrichment results downloaded successfully:")
#'     print(head(enrichment_data))
#'   } else {
#'     message("Failed to retrieve enrichment results.")
#'   }
#' } else if (!is.null(submission_info$submission_response$status) &&
#'            submission_info$submission_response$status == "error") {
#'   message(paste("API Error during submission:", submission_info$submission_response$message))
#' } else {
#'   message("Job submission was not successful, job_id is NULL, or status unknown.")
#'   # It's helpful to print the whole submission_info for debugging in this case
#'   print(submission_info)
#' }
#' }


# https://version-12-0.string-db.org/api/json/valuesranks_enrichment_status?api_key=bsjXYSW0kKTt&job_id=brsuCMHhuVNz

# [{"job_id": "brsuCMHhuVNz", "creation_time": "2025-05-07 15:25:27", "string_version": "12.0", "status": "success", "message": "Job finished", "page_url": "https://version-12-0.string-db.org/cgi/globalenrichment?networkId=bNEXfEymvDsZ", "download_url": "https://version-12-0.string-db.org/api/tsv/downloadenrichmentresults?networkId=bNEXfEymvDsZ", "graph_url": "https://version-12-0.string-db.org/api/image/enrichmentfigure?networkId=bNEXfEymvDsZ"}]

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
  checkmate::assertString(species,  .var.name = "species")
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


#' Run STRING DB Rank Enrichment Analysis for MOFA Factors or Similar Data
#'
#' @description
#' This function automates the process of submitting data to the STRING database
#' for values/ranks enrichment analysis, polling for results, and saving the
#' output. It is a wrapper around `submitStringDBEnrichment` and
#' `retrieveStringDBEnrichmentResults`, tailored for analyses where inputs
#' might come from MOFA factor weights or similar ranked lists.
#'
#' Results, including the enrichment table, a URL list, and a graph image,
#' are saved to the specified `results_dir`.
#'
#' @param input_table A data frame containing the data to be submitted for enrichment.
#'   Must include columns specified by `identifier_column_name` and `value_column_name`.
#' @param identifier_column_name Character string. The name of the column in
#'   `input_table` containing the identifiers (e.g., protein IDs, gene symbols).
#'   Default: `"protein_id"`.
#' @param value_column_name Character string. The name of the column in `input_table`
#'   containing the numerical values (e.g., scores, weights, logFC) associated
#'   with each identifier. This column must be numeric. Default: `"score"`.
#' @param result_label Character string. A label used for naming output files
#'   (e.g., "MOFA_Factor1_Proteomics").
#' @param results_dir Character string. The path to the directory where enrichment
#'   results (table, URL file, graph image) will be saved. The directory will be
#'   created if it doesn't exist.
#' @param api_key Character string or NULL. Your personal STRING API key.
#'   If NULL, the API might work with limitations or require it for full access.
#'   Default: `NULL`.
#' @param species Character string. NCBI/STRING species identifier.
#'   Default: `"9606"` (Homo sapiens).
#' @param ge_fdr Numeric. FDR threshold for gene expression enrichment.
#'   Default: `0.05`.
#' @param ge_enrichment_rank_direction Integer. Direction for enrichment rank
#'   (-1, 0, or 1). Default: `-1`.
#' @param polling_interval_seconds Numeric. Seconds to wait between polling attempts
#'   for job status. Default: `10`.
#' @param max_polling_attempts Numeric. Maximum polling attempts before timing out.
#'   Default: `30`.
#'
#' @return A data frame containing the enrichment results from STRING DB
#'   (specifically `output_tbl$enrichment_data` from the retrieval step).
#'   Returns `NULL` or an empty structure if the process fails.
#'
#' @seealso \code{\link{submitStringDBEnrichment}}, \code{\link{retrieveStringDBEnrichmentResults}}
#'
#' @importFrom vroom vroom_write
#' @importFrom readr write_lines
#' @export
runOneStringDbRankEnrichment <- function( input_table
                                         ,  result_label
                                         , api_key = NULL
                                         , species = "9606"
                                         , ge_fdr = 0.05
                                         , ge_enrichment_rank_direction = -1
                                         , polling_interval_seconds = 10
                                         , max_polling_attempts = 30) {

  stringdb_input_table <-  input_table |>
    mutate( score = sign(log2FC) * -log10(fdr_qvalue)) |>
    relocate(score, .after="log2FC") |>
    arrange(desc(log2FC)) |>
    mutate( protein_id = purrr::map_chr(Protein.Ids, ~str_split(.x, ":")[[1]][1])) |>
    relocate(protein_id, .after="Protein.Ids")


  parsed_response <- submitStringDBEnrichment (input_data_frame = stringdb_input_table ,
                                               identifier_column_name = "protein_id",
                                               value_column_name = "score",
                                               caller_identity = result_label,
                                               api_key = api_key,
                                               species = species,
                                               ge_fdr = ge_fdr,
                                               ge_enrichment_rank_direction = ge_enrichment_rank_direction)


  output_tbl <- retrieveStringDBEnrichmentResults( submission_info = parsed_response,
                                                   polling_interval_seconds = polling_interval_seconds,
                                                   max_polling_attempts = max_polling_attempts)

  write_lines(c("page_url", output_tbl$page_url
                , "download_url" , output_tbl$download_url
                , "graph_url" , output_tbl$graph_url)
              , file.path( results_dir, "functional_enrichment_string_db", paste0( result_label, "_string_enrichment_page_url.txt") ))

  vroom::vroom_write( output_tbl$enrichment_data
                      , file = file.path( results_dir
                                          , "functional_enrichment_string_db"
                                          , paste0( result_label, "_string_enrichment_results.tab") ))

  dir.create( file.path( results_dir, "functional_enrichment_string_db" ), showWarnings = TRUE, recursive = TRUE)
  writeBin(output_tbl$graph_image_content
           , file.path( results_dir , "functional_enrichment_string_db", paste0( result_label, "string_enrichment_graph.png") ))

  return(output_tbl$enrichment_data)

}




#' Run STRING DB Rank Enrichment Analysis for MOFA Factors or Similar Data
#'
#' @description
#' This function automates the process of submitting data to the STRING database
#' for values/ranks enrichment analysis, polling for results, and saving the
#' output. It is a wrapper around `submitStringDBEnrichment` and
#' `retrieveStringDBEnrichmentResults`, tailored for analyses where inputs
#' might come from MOFA factor weights or similar ranked lists.
#'
#' Results, including the enrichment table, a URL list, and a graph image,
#' are saved to the specified `results_dir`.
#'
#' @param input_table A data frame containing the data to be submitted for enrichment.
#'   Must include columns specified by `identifier_column_name` and `value_column_name`.
#' @param identifier_column_name Character string. The name of the column in
#'   `input_table` containing the identifiers (e.g., protein IDs, gene symbols).
#'   Default: `"protein_id"`.
#' @param value_column_name Character string. The name of the column in `input_table`
#'   containing the numerical values (e.g., scores, weights, logFC) associated
#'   with each identifier. This column must be numeric. Default: `"score"`.
#' @param result_label Character string. A label used for naming output files
#'   (e.g., "MOFA_Factor1_Proteomics").
#' @param results_dir Character string. The path to the directory where enrichment
#'   results (table, URL file, graph image) will be saved. The directory will be
#'   created if it doesn't exist.
#' @param api_key Character string or NULL. Your personal STRING API key.
#'   If NULL, the API might work with limitations or require it for full access.
#'   Default: `NULL`.
#' @param species Character string. NCBI/STRING species identifier.
#'   Default: `"9606"` (Homo sapiens).
#' @param ge_fdr Numeric. FDR threshold for gene expression enrichment.
#'   Default: `0.05`.
#' @param ge_enrichment_rank_direction Integer. Direction for enrichment rank
#'   (-1, 0, or 1). Default: `-1`.
#' @param polling_interval_seconds Numeric. Seconds to wait between polling attempts
#'   for job status. Default: `10`.
#' @param max_polling_attempts Numeric. Maximum polling attempts before timing out.
#'   Default: `30`.
#'
#' @return A data frame containing the enrichment results from STRING DB
#'   (specifically `output_tbl$enrichment_data` from the retrieval step).
#'   Returns `NULL` or an empty structure if the process fails.
#'
#' @seealso \code{\link{submitStringDBEnrichment}}, \code{\link{retrieveStringDBEnrichmentResults}}
#'
#' @importFrom vroom vroom_write
#' @importFrom readr write_lines
#' @export
runOneStringDbRankEnrichmentMofa <- function( input_table
                                           ,   identifier_column_name = "protein_id"
                                           ,   value_column_name = "score"
                                          ,  result_label
                                          , results_dir                                                                                   
                                          , api_key = NULL
                                          , species = "9606"
                                          , ge_fdr = 0.05
                                          , ge_enrichment_rank_direction = -1
                                          , polling_interval_seconds = 10
                                          , max_polling_attempts = 30) {
  
  
  
  parsed_response <- submitStringDBEnrichment (input_data_frame = input_table ,
                                               identifier_column_name = identifier_column_name,
                                               value_column_name = value_column_name,
                                               caller_identity = result_label,
                                               api_key = api_key,
                                               species = species,
                                               ge_fdr = ge_fdr,
                                               ge_enrichment_rank_direction = ge_enrichment_rank_direction)
  
  
  output_tbl <- retrieveStringDBEnrichmentResults( submission_info = parsed_response,
                                                   polling_interval_seconds = polling_interval_seconds,
                                                   max_polling_attempts = max_polling_attempts)
  
  write_lines(c("page_url", output_tbl$page_url
                , "download_url" , output_tbl$download_url
                , "graph_url" , output_tbl$graph_url)
              , file.path( results_dir,  paste0( result_label, "_string_enrichment_page_url.txt") ))
  
  vroom::vroom_write( output_tbl$enrichment_data
                      , file = file.path( results_dir
                                          
                                          , paste0( result_label, "_string_enrichment_results.tab") ))
  
  dir.create( file.path( results_dir), showWarnings = TRUE, recursive = TRUE)
  writeBin(output_tbl$graph_image_content
           , file.path( results_dir , paste0( result_label, "string_enrichment_graph.png") ))
  
  return(output_tbl$enrichment_data)
  
}

# https://version-12-0.string-db.org/api/json/valuesranks_enrichment_status?api_key=bsjXYSW0kKTt&job_id=brsuCMHhuVNz

# [{"job_id": "brsuCMHhuVNz", "creation_time": "2025-05-07 15:25:27", "string_version": "12.0", "status": "success", "message": "Job finished", "page_url": "https://version-12-0.string-db.org/cgi/globalenrichment?networkId=bNEXfEymvDsZ", "download_url": "https://version-12-0.string-db.org/api/tsv/downloadenrichmentresults?networkId=bNEXfEymvDsZ", "graph_url": "https://version-12-0.string-db.org/api/image/enrichmentfigure?networkId=bNEXfEymvDsZ"}]

#' Generate a Bar Graph of STRING DB Functional Enrichment Results
#'
#' @description
#' This function takes a data frame of STRING DB enrichment results and
#' generates a faceted bar graph. The graph displays the enrichment score for
#' terms, with points overlaid indicating the -log10(falseDiscoveryRate)
#' (color) and the number of genes mapped (size). Results are faceted by
#' 'category' and 'comparison'.
#'
#' @param input_table A data frame containing functional enrichment results.
#'   It is expected to have the following columns:
#'   - `comparison`: Character, identifier for the comparison group.
#'   - `termDescription`: Character, description of the enriched term.
#'   - `enrichmentScore`: Numeric, the enrichment score for the term.
#'   - `falseDiscoveryRate`: Numeric, the False Discovery Rate for the term.
#'   - `genesMapped`: Numeric, the number of genes mapped to the term.
#'   - `category`: Character or Factor, the category of the enrichment (e.g., GO BP, KEGG).
#'
#' @return A `ggplot` object representing the enrichment bar graph.
#'
#' @examples
#' \dontrun{
#' # Assume 'enrichment_results_df' is a data frame structured as described above
#' # For example:
#' # enrichment_results_df <- data.frame(
#' #   comparison = rep(c("GroupA_vs_Control", "GroupB_vs_Control"), each = 2),
#' #   termDescription = c("Immune response", "Metabolic process", "Immune response", "Cell cycle"),
#' #   enrichmentScore = c(2.5, 1.8, 3.1, 2.0),
#' #   falseDiscoveryRate = c(0.01, 0.045, 0.005, 0.02),
#' #   genesMapped = c(50, 30, 65, 40),
#' #   category = rep(c("GO Biological Process", "KEGG Pathway"), times = 2)
#' # )
#' #
#' # enrichment_plot <- printStringDbFunctionalEnrichmentBarGraph(enrichment_results_df)
#' # print(enrichment_plot)
#' }
#' @export 
printStringDbFunctionalEnrichmentBarGraph <- function( input_table){
  
  
  max_enrich_score_group_tbl <- input_table  |>
    group_by(comparison, termDescription) |>
    summarise(counts = n(), max_enrichment_score = max(enrichmentScore), min_fdr = min(falseDiscoveryRate)) |>
    ungroup() 
  
  output_group_enrichment_table <- input_table  |>
    inner_join(max_enrich_score_group_tbl
               , by = join_by(comparison == comparison ,
                              termDescription == termDescription ,
                              enrichmentScore == max_enrichment_score,
                              falseDiscoveryRate == min_fdr)) |>
    mutate(termDescriptionAbbrev = str_sub(termDescription, 1, 100) ) |>
    ggplot(aes(x = termDescriptionAbbrev
               , y = enrichmentScore)) +
    geom_bar(stat = "identity", width=0.1) +
    geom_point(aes(x = termDescriptionAbbrev
                   , y = enrichmentScore
                   , colour = -log10(falseDiscoveryRate)
                   , size= (genesMapped)))  +
    theme(strip.text.y = element_text(angle = 0)) +
    coord_flip() +
    facet_grid(  category ~ comparison
                 , scales = "free_y", space="free")   
  
  output_group_enrichment_table
}
