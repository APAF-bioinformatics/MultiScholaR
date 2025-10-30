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
#' @param species Character or numeric: NCBI/STRING species identifier. Default is 9606 (Homo sapiens).
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



#' Run STRING DB Rank Enrichment Analysis
#'
#' @description
#' Performs STRING DB enrichment analysis on ranked data. This function processes
#' differential expression results by calculating a score from log fold change and
#' FDR values, extracts protein identifiers, and submits them to STRING DB for
#' functional enrichment analysis. Results are automatically saved to files.
#'
#' @param input_table A data frame containing differential expression results with
#'   columns: `log2FC` (log fold change), `fdr_qvalue` (adjusted p-value), and
#'   `Protein.Ids` (protein identifiers, potentially with isoform information).
#' @param result_label Character string: A label for the analysis results, used
#'   in output file names.
#' @param pathway_dir Character string: The pathway directory where results will be saved.
#'   All outputs are saved to pathway_dir/string_db/ subdirectory.
#' @param api_key Character string: Your personal STRING API key.
#' @param species Character or numeric: NCBI/STRING species identifier. Default is 9606 (Homo sapiens).
#' @param ge_fdr Numeric: FDR threshold for gene expression enrichment. Default is 0.05.
#' @param ge_enrichment_rank_direction Integer: Direction for enrichment rank.
#'   (-1, 0, or 1). Default is -1.
#' @param polling_interval_seconds Numeric: Seconds to wait between polling attempts.
#'   Default is 10.
#' @param max_polling_attempts Numeric: Maximum number of polling attempts. Default is 30.
#'
#' @return A data frame containing the STRING DB enrichment results.
#'
#' @details
#' This function saves enrichment results to the specified pathway directory:
#' - `{pathway_dir}/string_db/{result_label}_string_enrichment_page_url.txt` - URLs for STRING DB results
#' - `{pathway_dir}/string_db/{result_label}_string_enrichment_results.tab` - Tab-delimited enrichment results
#' - `{pathway_dir}/string_db/{result_label}_string_enrichment_graph.png` - Network graph image
#'
#' @importFrom dplyr mutate relocate arrange desc
#' @importFrom stringr str_split
#' @importFrom purrr map_chr
#' @importFrom readr write_lines
#' @importFrom vroom vroom_write
#'
#' @examples
#' \dontrun{
#' # Assume 'de_results' is a data frame with log2FC, fdr_qvalue, and Protein.Ids columns
#' enrichment_results <- runOneStringDbRankEnrichment(
#'   input_table = de_results,
#'   result_label = "treatment_vs_control",
#'   pathway_dir = "path/to/pathway_enrichment",
#'   api_key = "YOUR_API_KEY"
#' )
#' }
#'
#' @export
runOneStringDbRankEnrichment <- function( input_table
                                          ,  result_label
                                          , pathway_dir
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

  enrichment_dir <- file.path(pathway_dir, "string_db")
  dir.create(enrichment_dir, showWarnings = FALSE, recursive = TRUE)
  
  write_lines(c("page_url", output_tbl$page_url
                , "download_url" , output_tbl$download_url
                , "graph_url" , output_tbl$graph_url)
              , file.path(enrichment_dir, paste0(result_label, "_string_enrichment_page_url.txt")))

  vroom::vroom_write( output_tbl$enrichment_data
                      , file = file.path(enrichment_dir
                                          , paste0(result_label, "_string_enrichment_results.tab")))

  writeBin(output_tbl$graph_image_content
           , file.path(enrichment_dir, paste0(result_label, "_string_enrichment_graph.png")))

  return(output_tbl$enrichment_data)

}





#' Run STRING DB Rank Enrichment Analysis for MOFA
#'
#' @description
#' Performs STRING DB enrichment analysis for MOFA (Multi-Omics Factor Analysis)
#' results. This function takes pre-processed data with identifier and value columns
#' and submits them to STRING DB for functional enrichment analysis. Results are
#' automatically saved to the specified results directory.
#'
#' @param input_table A data frame containing the data to analyze with identifier
#'   and value columns as specified by the column name parameters.
#' @param identifier_column_name Character string: The name of the column in
#'   `input_table` that contains the protein/gene identifiers. Default is "protein_id".
#' @param value_column_name Character string: The name of the column in `input_table`
#'   that contains the numerical values (e.g., scores, weights) associated with
#'   each identifier. Default is "score".
#' @param result_label Character string: A label for the analysis results, used
#'   in output file names.
#' @param results_dir Character string: The directory path where results should be saved.
#' @param api_key Character string: Your personal STRING API key.
#' @param species Character or numeric: NCBI/STRING species identifier. Default is 9606 (Homo sapiens).
#' @param ge_fdr Numeric: FDR threshold for gene expression enrichment. Default is 0.05.
#' @param ge_enrichment_rank_direction Integer: Direction for enrichment rank.
#'   (-1, 0, or 1). Default is -1.
#' @param polling_interval_seconds Numeric: Seconds to wait between polling attempts.
#'   Default is 10.
#' @param max_polling_attempts Numeric: Maximum number of polling attempts. Default is 30.
#'
#' @return A data frame containing the STRING DB enrichment results.
#'
#' @details
#' This function saves results directly to the specified `results_dir`:
#' - URLs for STRING DB results page, download, and graph
#' - Tab-delimited enrichment results file
#' - PNG image of the enrichment network graph
#'
#' Unlike `runOneStringDbRankEnrichment`, this function expects the input data to
#' already be properly formatted with identifier and value columns.
#'
#' @importFrom readr write_lines
#' @importFrom vroom vroom_write
#'
#' @examples
#' \dontrun{
#' # Assume 'mofa_results' is a data frame with protein_id and score columns
#' enrichment_results <- runOneStringDbRankEnrichmentMofa(
#'   input_table = mofa_results,
#'   identifier_column_name = "protein_id",
#'   value_column_name = "score",
#'   result_label = "MOFA_factor_1",
#'   results_dir = "/path/to/results",
#'   api_key = "YOUR_API_KEY"
#' )
#' }
#'
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

printStringDbFunctionalEnrichmentBarGraph <- function (input_table, word_limit = 10)
{
  plot_data <- input_table |>
    group_by(comparison, category, termDescription) |>
    arrange( desc( enrichmentScore), falseDiscoveryRate ) |>
    #summarise(enrichmentScore = max(enrichmentScore), falseDiscoveryRate =min(falseDiscoveryRate), direction = first(direction), genesMapped = max(genesMapped)) |>
    dplyr::slice(1) |>
    ungroup() |>
    mutate(termDescriptionAbbrev = sapply(strsplit(as.character(termDescription), " "), function(x) {
      if (length(x) > word_limit) {
        paste(c(head(x, word_limit), "..."), collapse = " ")
      } else {
        paste(x, collapse = " ")
      }
    })) |>
    #mutate(facet_group = interaction(category, comparison)) |>
    # Ensure 'direction' is a factor with all possible levels
    mutate(direction = factor(direction, levels = c("top", "bottom", "both ends")))

  output_group_enrichment_table <- ggplot(plot_data,
                                          aes(y =  reorder_within(  termDescriptionAbbrev, enrichmentScore, list( category )  ) , x = enrichmentScore )) +
    geom_bar(aes(fill=direction),stat = "identity",  width = 0.1) +
    scale_fill_manual(
      name = "Direction",
      values = c("top" = "red", "bottom" = "blue", "both ends" = "grey"),
      drop = FALSE
    ) +
    geom_point(aes(y = reorder_within( termDescriptionAbbrev, enrichmentScore, list( category  )) ,
                   x = enrichmentScore, colour = -log10(falseDiscoveryRate),
                   size = (genesMapped))) + theme(strip.text.y = element_text(angle = 0)) +
    facet_grid( category ~ comparison, scales = "free_y",
                               space = "free") +
    scale_y_reordered() +
    ylab("Term Description") +
    xlab("Enrichment Score")
  output_group_enrichment_table
}



#' Run Metabolomics Enrichment Analysis for MOFA Factors
#'
#' @param weights Data frame containing MOFA weights
#' @param metabolomics_obj Metabolomics S4 object containing ChEBI IDs
#' @param mapping_table The metabolite ID mapping table from AnnotationHub
#' @param project_dirs Project directories structure
#' @param omic_type Omics type for directory access
#' @param experiment_label Experiment label for directory access
#' @param assay_name Name of the assay ("metabolome_lc" or "metabolome_gc")
#' @param kegg_species_code KEGG species code (e.g., "kpn" for Klebsiella pneumoniae)
#' @param reactome_organism Organism name to use for Reactome filtering (e.g., "Homo sapiens")
#' @return A data frame with enrichment results formatted for visualization
#' @export
runMetabolomicsEnrichmentAnalysis <- function(weights, 
                                             metabolomics_obj,
                                             mapping_table,
                                             project_dirs,
                                             omic_type,
                                             experiment_label,
                                             assay_name,
                                             kegg_species_code = "kpn",
                                             reactome_organism = NULL) {
  
  message(sprintf("--- Entering runMetabolomicsEnrichmentAnalysis ---"))
  message(sprintf("   runMetabolomicsEnrichmentAnalysis Args: assay_name = %s", assay_name))
  message(sprintf("   runMetabolomicsEnrichmentAnalysis Args: kegg_species_code = %s", kegg_species_code))
  message(sprintf("   runMetabolomicsEnrichmentAnalysis Args: reactome_organism = %s", as.character(reactome_organism)))
  
  # Use exact directory paths with paste0 construction
  data_dir <- project_dirs[[paste0(omic_type, "_", experiment_label)]]$data_dir
  results_dir <- project_dirs[[paste0(omic_type, "_", experiment_label)]]$integration_enrichment_plots_dir
  
  # Create directory structure if it doesn't exist
  dir.create(file.path(results_dir, assay_name), recursive = TRUE, showWarnings = FALSE)
  
  # 1. Extract metabolite weights from MOFA for the specific assay
  message("   runMetabolomicsEnrichmentAnalysis Step: Extracting metabolite weights for assay...")
  metabolite_weights <- weights |>
    dplyr::filter(view == assay_name & factor == "Factor1") |>
    mutate(feature = str_replace_all(feature, paste0("_", assay_name), ""))
  
  message(sprintf("   runMetabolomicsEnrichmentAnalysis: Found %d metabolite weights", nrow(metabolite_weights)))
  message("   runMetabolomicsEnrichmentAnalysis: Sample weights:")
  print(head(metabolite_weights, 5))
  
  # 2. Get ChEBI IDs from metabolomics object
  message("   runMetabolomicsEnrichmentAnalysis Step: Getting ChEBI IDs from metabolomics object...")
  assay_index <- if(assay_name == "metabolome_lc") 1 else 2
  
  # Extract ChEBI IDs and corresponding metabolite names
  chebi_mapping <- metabolomics_obj@metabolite_data[[assay_index]] |>
    dplyr::filter(!stringr::str_detect(database_identifier, "^ITSD") & 
                 !stringr::str_detect(metabolite_identification, "^ITSD")) |>
    dplyr::select(metabolite, database_identifier) |>
    # Extract just the ChEBI ID number from the database_identifier column
    dplyr::mutate(chebi_id = stringr::str_extract(database_identifier, "CHEBI:\\d+"))
  
  message(sprintf("   runMetabolomicsEnrichmentAnalysis: Extracted %d ChEBI IDs", nrow(chebi_mapping)))
  message("   runMetabolomicsEnrichmentAnalysis: Sample ChEBI mappings:")
  print(head(chebi_mapping, 5))
  
  # Check missing ChEBI IDs
  missing_chebi_count <- sum(is.na(chebi_mapping$chebi_id))
  if (missing_chebi_count > 0) {
    message(sprintf("   runMetabolomicsEnrichmentAnalysis WARNING: %d metabolites missing ChEBI IDs", missing_chebi_count))
    message("Sample entries with missing ChEBI IDs:")
    print(head(chebi_mapping[is.na(chebi_mapping$chebi_id),], 5))
  }
  
  # 3. Join metabolite weights with ChEBI IDs
  message("   runMetabolomicsEnrichmentAnalysis Step: Joining weights with ChEBI IDs...")
  metabolite_weights_with_ids <- metabolite_weights |>
    left_join(chebi_mapping, by = c("feature" = "metabolite"))
  
  message(sprintf("   runMetabolomicsEnrichmentAnalysis: Joined table has %d rows", nrow(metabolite_weights_with_ids)))
  message("   runMetabolomicsEnrichmentAnalysis: Sample joined data:")
  print(head(metabolite_weights_with_ids, 5))
  
  # 4. Filter out entries without ChEBI IDs
  message("   runMetabolomicsEnrichmentAnalysis Step: Filtering out entries without ChEBI IDs...")
  metabolite_weights_with_ids <- metabolite_weights_with_ids |>
    dplyr::filter(!is.na(chebi_id))
  
  message(sprintf("   runMetabolomicsEnrichmentAnalysis: After filtering, %d entries remain", nrow(metabolite_weights_with_ids)))
  
  # Create ranked gene list for GSEA
  message("   runMetabolomicsEnrichmentAnalysis Step: Creating ranked list for GSEA...")
  ranked_list <- metabolite_weights_with_ids |>
    dplyr::arrange(desc(value)) |>
    dplyr::pull(value, name = chebi_id)
  
  message(sprintf("   runMetabolomicsEnrichmentAnalysis: Created ranked list with %d entries", length(ranked_list)))
  message("   runMetabolomicsEnrichmentAnalysis: Top metabolites in ranked list:")
  print(head(ranked_list, 5))
  
  # Run KEGG pathway enrichment analysis
  message("   runMetabolomicsEnrichmentAnalysis Step: Running KEGG pathway enrichment...")
  kegg_results <- runKeggEnrichment(
    ranked_list = ranked_list, 
    mapping_table = mapping_table,
    project_dirs = project_dirs,
    omic_type = omic_type,
    experiment_label = experiment_label,
    assay_name = assay_name,
    kegg_species_code = kegg_species_code
  )
  
  message(sprintf("   runMetabolomicsEnrichmentAnalysis: KEGG analysis returned %d results", nrow(kegg_results)))
  
  # Run Reactome pathway enrichment analysis  
  message("   runMetabolomicsEnrichmentAnalysis Step: Running Reactome pathway enrichment...")
  reactome_results <- runReactomeEnrichment(
    ranked_list = ranked_list, 
    mapping_table = mapping_table,
    project_dirs = project_dirs,
    omic_type = omic_type,
    experiment_label = experiment_label,
    assay_name = assay_name,
    reactome_organism = reactome_organism
  )
  
  message(sprintf("   runMetabolomicsEnrichmentAnalysis: Reactome analysis returned %d results", nrow(reactome_results)))
  
  # Combine results and format for visualization
  message("   runMetabolomicsEnrichmentAnalysis Step: Combining KEGG and Reactome results...")
  combined_results <- bind_rows(
    kegg_results |> mutate(category = "KEGG"),
    reactome_results |> mutate(category = "Reactome")
  ) |> 
  distinct()  # Ensure no duplicates after combining
  
  message(sprintf("   runMetabolomicsEnrichmentAnalysis: Combined %d results total", nrow(combined_results)))
  
  # Save the combined results
  message(sprintf("   runMetabolomicsEnrichmentAnalysis Step: Saving results to %s...", 
                 file.path(results_dir, paste0(assay_name, "_enrichment_results.tab"))))
  vroom::vroom_write(combined_results, 
                    file.path(results_dir, paste0(assay_name, "_enrichment_results.tab")))
  
  message("--- Exiting runMetabolomicsEnrichmentAnalysis ---")
  return(combined_results)
}

' Run KEGG Pathway Enrichment Analysis for Metabolites
#'
#' @param ranked_list Named numeric vector with ChEBI IDs as names and values as weights
#' @param mapping_table Metabolite ID mapping table
#' @param project_dirs Project directories structure
#' @param omic_type Omics type for directory access
#' @param experiment_label Experiment label for directory access
#' @param assay_name Name of the assay
#' @param kegg_species_code KEGG species code (e.g., "kpn" for Klebsiella pneumoniae)
#' @return A data frame with KEGG pathway enrichment results
#' @export
runKeggEnrichment <- function(ranked_list, 
                             mapping_table, 
                             project_dirs,
                             omic_type,
                             experiment_label,
                             assay_name,
                             kegg_species_code = "kpn") {
  
  message(sprintf("--- Entering runKeggEnrichment ---"))
  message(sprintf("   runKeggEnrichment Args: assay_name = %s", assay_name))
  message(sprintf("   runKeggEnrichment Args: kegg_species_code = %s", kegg_species_code))
  message(sprintf("   runKeggEnrichment Args: ranked_list length = %d", length(ranked_list)))
  
  # Use exact directory paths with paste0 construction
  data_dir <- project_dirs[[paste0(omic_type, "_", experiment_label)]]$data_dir
  message(sprintf("   runKeggEnrichment: Using data_dir = %s", data_dir))
  
  # Print some debugging information
  message("KEGG Enrichment - Input ranked_list summary:")
  message(paste("Number of metabolites:", length(ranked_list)))
  message(paste("Sample ChEBI IDs:", paste(head(names(ranked_list), 3), collapse=", ")))
  
  # Load the necessary mappings with better debugging
  message("   runKeggEnrichment Step: Loading KEGG to ChEBI mappings...")
  kegg_to_chebi <- mapping_table |>
    dplyr::filter(!is.na(KEGG) & !is.na(ChEBI)) 
  
  message(paste("Number of KEGG-ChEBI mappings:", nrow(kegg_to_chebi)))
  message("   runKeggEnrichment Step: Sample of KEGG-ChEBI mappings:")
  print(head(kegg_to_chebi, 5))
  
  # Format ChEBI IDs to match our format
  message("   runKeggEnrichment Step: Formatting ChEBI IDs...")
  kegg_to_chebi <- kegg_to_chebi |>
    dplyr::select(KEGG, ChEBI) |>
    dplyr::mutate(ChEBI = paste0("CHEBI:", ChEBI))
  
  message("   runKeggEnrichment Step: Formatted KEGG-ChEBI mappings:")
  print(head(kegg_to_chebi, 5))
  
  # Load pathway to compound mapping
  message(sprintf("   runKeggEnrichment Step: Loading pathway to compound mapping from %s...", file.path(data_dir, "species_specific_pathway_to_compound_tbl.tab")))
  
  # Check if file exists
  if (!file.exists(file.path(data_dir, "species_specific_pathway_to_compound_tbl.tab"))) {
    message("   runKeggEnrichment ERROR: File 'species_specific_pathway_to_compound_tbl.tab' does not exist!")
    return(tibble(
      termDescription = character(),
      enrichmentScore = numeric(),
      falseDiscoveryRate = numeric(),
      genesMapped = integer(),
      comparison = character()
    ))
  }
  
  species_specific_pathway_to_compound_tbl <- vroom::vroom(
    file.path(data_dir, "species_specific_pathway_to_compound_tbl.tab"),
    show_col_types = FALSE
  )
  
  message(paste("Number of pathway-compound mappings:", nrow(species_specific_pathway_to_compound_tbl)))
  message("   runKeggEnrichment Step: Pathway to compound column names:")
  print(colnames(species_specific_pathway_to_compound_tbl))
  
  # Debug - show first few compounds and pathways
  message("First few compounds:")
  message(paste(head(species_specific_pathway_to_compound_tbl$compound, 5), collapse=", "))
  
  message("First few pathways:")
  message(paste(head(unique(species_specific_pathway_to_compound_tbl$pathway), 5), collapse=", "))
  
  # Directly use KEGG compound IDs from the pathway mapping
  message("   runKeggEnrichment Step: Extracting unique compound IDs...")
  compound_ids <- sort(unique(species_specific_pathway_to_compound_tbl$compound))
  message(paste("Total unique compounds in KEGG:", length(compound_ids)))
  message("Sample compound IDs in KEGG:")
  message(paste(head(compound_ids, 10), collapse=", "))
  
  # Cross-check with our ranked list - are there any direct KEGG IDs in our list?
  message("   runKeggEnrichment Step: Checking for direct KEGG IDs in ranked list...")
  kegg_ids_in_ranked_list <- intersect(compound_ids, names(ranked_list))
  message(paste("Direct KEGG IDs in ranked list:", length(kegg_ids_in_ranked_list)))
  if (length(kegg_ids_in_ranked_list) > 0) {
    message("Sample direct KEGG IDs found in ranked list:")
    message(paste(head(kegg_ids_in_ranked_list, 5), collapse=", "))
  }
  
  # If there are no direct KEGG IDs, we'll need to rely on the mapping
  if (length(kegg_ids_in_ranked_list) == 0) {
    message("   runKeggEnrichment Step: No direct KEGG IDs in ranked list, using ChEBI to KEGG mapping...")
    
    # Print the structure of the first few ChEBI IDs to check format
    chebi_ids_in_ranked_list <- names(ranked_list)[grepl("CHEBI:", names(ranked_list))]
    message(paste("ChEBI IDs in ranked list:", length(chebi_ids_in_ranked_list)))
    if (length(chebi_ids_in_ranked_list) > 0) {
      message(paste("Sample ChEBI IDs:", paste(head(chebi_ids_in_ranked_list, 3), collapse=", ")))
    }
    
    # Check how many of our ChEBI IDs are in the mapping table
    chebi_ids_in_mapping <- intersect(chebi_ids_in_ranked_list, kegg_to_chebi$ChEBI)
    message(paste("ChEBI IDs found in mapping table:", length(chebi_ids_in_mapping)))
    if (length(chebi_ids_in_mapping) > 0) {
      message("Sample ChEBI IDs found in mapping table:")
      message(paste(head(chebi_ids_in_mapping, 5), collapse=", "))
    }
    
    # Create a direct lookup from ChEBI to KEGG
    message("   runKeggEnrichment Step: Creating ChEBI to KEGG lookup...")
    chebi_to_kegg_lookup <- kegg_to_chebi |>
      dplyr::filter(ChEBI %in% chebi_ids_in_ranked_list)
    
    message(paste("ChEBI-KEGG mappings for our ranked list:", nrow(chebi_to_kegg_lookup)))
    if (nrow(chebi_to_kegg_lookup) > 0) {
      message("   runKeggEnrichment Step: Sample ChEBI-KEGG mappings:")
      print(head(chebi_to_kegg_lookup, 5))
    }
    
    # If we have mappings, create a new ranked list with KEGG IDs
    if (nrow(chebi_to_kegg_lookup) > 0) {
      # Create a new ranked list with KEGG IDs
      message("   runKeggEnrichment Step: Creating ranked list with KEGG IDs...")
      kegg_ranked_list <- numeric()
      for (i in 1:nrow(chebi_to_kegg_lookup)) {
        chebi_id <- chebi_to_kegg_lookup$ChEBI[i]
        kegg_id <- chebi_to_kegg_lookup$KEGG[i]
        if (chebi_id %in% names(ranked_list)) {
          kegg_ranked_list[kegg_id] <- ranked_list[chebi_id]
        }
      }
      
      message(paste("Created KEGG ranked list with", length(kegg_ranked_list), "entries"))
      if (length(kegg_ranked_list) > 0) {
        message("   runKeggEnrichment Step: Sample KEGG IDs in ranked list:")
        message(paste(head(names(kegg_ranked_list), 5), collapse=", "))
      }
      
      # If we have KEGG IDs, use them for the analysis
      if (length(kegg_ranked_list) > 0) {
        # Sort the ranked list in decreasing order
        kegg_ranked_list <- sort(kegg_ranked_list, decreasing = TRUE)
        
        # Diagnostic: Check for ties in the data
        message("   runKeggEnrichment DEBUG: Analyzing score distribution...")
        # Count exact duplicates in values
        duplicate_count <- sum(duplicated(kegg_ranked_list))
        duplicate_percentage <- (duplicate_count / length(kegg_ranked_list)) * 100
        message(sprintf("   runKeggEnrichment DEBUG: Found %d exact duplicate values (%.2f%% of total)",
                       duplicate_count, duplicate_percentage))
        
        # Show basic statistics
        message(sprintf("   runKeggEnrichment DEBUG: Score range: [%.4f, %.4f], Median: %.4f", 
                       min(kegg_ranked_list), max(kegg_ranked_list), median(kegg_ranked_list)))
        
        # Add small random noise to break ties in the ranked list
        message("   runKeggEnrichment Step: Adding small random noise to break ties in ranked list...")
        set.seed(42) # For reproducibility
        original_kegg_ranked_list <- kegg_ranked_list # Store original for reference
        kegg_ranked_list <- kegg_ranked_list + runif(length(kegg_ranked_list), min=-0.0001, max=0.0001)
        # Re-sort to ensure order is preserved with the new values
        kegg_ranked_list <- sort(kegg_ranked_list, decreasing = TRUE)
        message(sprintf("   runKeggEnrichment: Added jitter to %d values to prevent ties", length(kegg_ranked_list)))
        
        # CRITICAL: Verify the ordering and distribution before GSEA
        message("   runKeggEnrichment CRITICAL CHECK: Verifying ranked list ordering and distribution...")
        positive_values <- sum(kegg_ranked_list > 0)
        negative_values <- sum(kegg_ranked_list < 0)
        zero_values <- sum(kegg_ranked_list == 0)
        message(sprintf("   runKeggEnrichment: Value distribution: %d positive, %d negative, %d zero", 
                       positive_values, negative_values, zero_values))
        
        # Check that the list is actually sorted (decreasing)
        is_decreasing <- all(diff(kegg_ranked_list) <= 0)
        message(sprintf("   runKeggEnrichment: List is properly sorted in decreasing order: %s", 
                       if(is_decreasing) "YES" else "NO - THIS IS A PROBLEM!"))
        
        # Show top and bottom values to verify ordering
        message("   runKeggEnrichment: Top 5 entries in ranked list:")
        for (i in 1:min(5, length(kegg_ranked_list))) {
          message(sprintf("     %s: %.6f", names(kegg_ranked_list)[i], kegg_ranked_list[i]))
        }
        message("   runKeggEnrichment: Bottom 5 entries in ranked list:")
        for (i in (length(kegg_ranked_list) - min(4, length(kegg_ranked_list) - 1)):length(kegg_ranked_list)) {
          message(sprintf("     %s: %.6f", names(kegg_ranked_list)[i], kegg_ranked_list[i]))
        }
        
        # CRITICAL SECTION - WHERE THE ISSUE LIKELY IS
        message("   runKeggEnrichment Step: Creating TERM2GENE mapping with KEGG IDs...")
        message("   runKeggEnrichment Step: Checking for pathway-compound matches...")
        
        # Check if compounds in pathway table have a prefix that's missing in our ranked list
        pathway_compounds <- head(species_specific_pathway_to_compound_tbl$compound, 10)
        ranked_compounds <- head(names(kegg_ranked_list), 10)
        
        message("Pathway compound format examples:")
        print(pathway_compounds)
        
        message("Ranked list compound format examples:")
        print(ranked_compounds)
        
        # FIX: FORCE PREFIX ADDITION REGARDLESS OF CHECKS
        message("   runKeggEnrichment DEBUG: Adding 'cpd:' prefix to all KEGG IDs in ranked list")
        prefixed_kegg_ids <- paste0("cpd:", names(kegg_ranked_list))
        prefixed_kegg_ranked_list <- kegg_ranked_list
        names(prefixed_kegg_ranked_list) <- prefixed_kegg_ids
        
        message(sprintf("   runKeggEnrichment DEBUG: Created prefixed ranked list with %d entries", length(prefixed_kegg_ranked_list)))
        message("Sample prefixed entries:")
        print(head(prefixed_kegg_ranked_list, 5))
        
        # Use the prefixed version for TERM2GENE mapping
        term2gene <- species_specific_pathway_to_compound_tbl |>
          dplyr::filter(compound %in% names(prefixed_kegg_ranked_list)) |>
          dplyr::select(pathway, compound) |>
          dplyr::rename(term = pathway, gene = compound)
        
        message(sprintf("   runKeggEnrichment DEBUG: After prefix fix, TERM2GENE has %d entries", nrow(term2gene)))
        
        if (nrow(term2gene) == 0) {
          message("   runKeggEnrichment WARNING: Still no matches after prefix fix!")
          
          # Debug check the values directly
          message("   runKeggEnrichment DEBUG: Direct value checking...")
          sample_pathway_compounds <- head(species_specific_pathway_to_compound_tbl$compound, 5)
          sample_prefixed_ids <- head(names(prefixed_kegg_ranked_list), 10)
          
          for (compound in sample_prefixed_ids) {
            matching_idx <- which(species_specific_pathway_to_compound_tbl$compound == compound)
            message(sprintf("   runKeggEnrichment DEBUG: Compound '%s' appears %d times in pathway table", 
                          compound, length(matching_idx)))
            if (length(matching_idx) > 0) {
              sample_pathways <- head(species_specific_pathway_to_compound_tbl$pathway[matching_idx], 3)
              message(sprintf("   runKeggEnrichment DEBUG: Found in pathways: %s", 
                            paste(sample_pathways, collapse=", ")))
            }
          }
        }
        
        # Use prefixed list regardless of whether mapping succeeded
        kegg_ranked_list <- prefixed_kegg_ranked_list
        
        # Create term2name directly in this path - crucial addition
        message("   runKeggEnrichment Step: Creating TERM2NAME mapping in direct path...")
        # Get all pathways from KEGG if not already done
        if (!exists("pathways_tbl")) {
          message("   runKeggEnrichment Step: Fetching pathways from KEGG...")
          pathways_tbl <- KEGGREST::keggList("pathway")
          
          message(sprintf("   runKeggEnrichment: Retrieved %d pathways from KEGG", length(pathways_tbl)))
          if (length(pathways_tbl) > 0) {
            message("Sample pathway entries:")
            print(head(pathways_tbl, 3))
          }
        }
        
        # Extract unique terms from term2gene
        term2name_initial <- term2gene |>
          dplyr::distinct(term)
        
        message(sprintf("   runKeggEnrichment: Found %d distinct terms", nrow(term2name_initial)))
        
        # Create a mapping table from pathway IDs (e.g., "path:kpn00010") to names
        # First extract the numeric part of each term
        message("   runKeggEnrichment Step: Extracting numeric parts from KEGG pathway IDs...")
        term2name_with_numeric <- term2name_initial |>
          dplyr::mutate(
            pathway_num = stringr::str_extract(term, "\\d+"),
            # Create generic format for lookup
            map_format = paste0("map", pathway_num)
          )
        
        message("   runKeggEnrichment DEBUG: Sample pathways with numeric parts:")
        print(head(term2name_with_numeric, 3))
        
        # Create mapping from pathway numbers to pathway names
        pathway_num_to_name <- data.frame(
          map_id = names(pathways_tbl),
          pathway_name = as.character(pathways_tbl)
        ) |>
          dplyr::mutate(
            map_format = stringr::str_extract(map_id, "map\\d+")
          )
        
        message("   runKeggEnrichment DEBUG: Sample pathway mappings:")
        print(head(pathway_num_to_name, 3))
        
        # Join the terms with pathway names based on numeric part
        message("   runKeggEnrichment Step: Joining with pathway names...")
        term2name <- term2name_with_numeric |>
          dplyr::left_join(pathway_num_to_name, by = "map_format") |>
          dplyr::select(term, name = pathway_name)
        
        # Handle any missing names (fallback to term ID)
        term2name <- term2name |>
          dplyr::mutate(name = ifelse(is.na(name), term, name))
        
        message("   runKeggEnrichment: Created improved term2name mapping:")
        print(head(term2name, 3))
        
        # Try running GSEA with tryCatch to handle errors
        message("   runKeggEnrichment Step: Running GSEA analysis...")
        
        # ABSOLUTE FINAL SAFETY CHECK - verify both term2gene and term2name exist
        if (!exists("term2gene") || nrow(term2gene) == 0) {
          message("   runKeggEnrichment CRITICAL ERROR: term2gene is missing or empty before GSEA call")
          # Can't proceed without term2gene
          return(tibble(
            termDescription = character(),
            enrichmentScore = numeric(),
            falseDiscoveryRate = numeric(),
            genesMapped = integer(),
            comparison = character()
          ))
        }
        
        if (!exists("term2name") || nrow(term2name) == 0) {
          message("   runKeggEnrichment WARNING: term2name is missing or empty before GSEA call, creating emergency fallback")
          # Create emergency fallback - but try to get proper names first
          # Get all pathways from KEGG if not already done
          if (!exists("pathways_tbl")) {
            message("   runKeggEnrichment Step: Emergency fetch of pathways from KEGG...")
            pathways_tbl <- tryCatch({
              KEGGREST::keggList("pathway")
            }, error = function(e) {
              message("   runKeggEnrichment ERROR: Could not fetch pathway names:", e$message)
              return(NULL)
            })
          }
          
          if (!is.null(pathways_tbl) && length(pathways_tbl) > 0) {
            # Extract term numbers and try to map
            term2name <- term2gene |> 
              dplyr::distinct(term) |>
              dplyr::mutate(
                pathway_num = stringr::str_extract(term, "\\d+"),
                map_format = paste0("map", pathway_num)
              )
            
            # Create mapping from pathway numbers to pathway names
            pathway_num_to_name <- data.frame(
              map_id = names(pathways_tbl),
              pathway_name = as.character(pathways_tbl)
            ) |>
              dplyr::mutate(
                map_format = stringr::str_extract(map_id, "map\\d+")
              )
            
            # Join and select final columns
            term2name <- term2name |>
              dplyr::left_join(pathway_num_to_name, by = "map_format") |>
              dplyr::select(term, name = pathway_name)
              
            # Handle any missing names (fallback to term ID)
            term2name <- term2name |>
              dplyr::mutate(name = ifelse(is.na(name), term, name))
          } else {
            # Complete fallback if KEGG fetch failed
            term2name <- term2gene |> 
              dplyr::distinct(term) |>
              dplyr::mutate(name = term)
          }
        }
        
        # Debug output
        message(sprintf("   runKeggEnrichment VERIFICATION: term2gene has %d rows, term2name has %d rows", 
                       nrow(term2gene), nrow(term2name)))
        
        gsea_result <- tryCatch({
          # FIXED: Use fgsea and remove nPerm parameter
          clusterProfiler::GSEA(
            geneList = kegg_ranked_list,
            TERM2GENE = term2gene,
            TERM2NAME = term2name,
            minGSSize = 3,
            maxGSSize = 500,
            pvalueCutoff = 0.05,  # MODIFIED: Increased from 0.1 to 0.05 for more lenient filtering
            pAdjustMethod = "fdr",
            verbose = TRUE,
            seed = TRUE,      # Set seed for reproducibility
            by = "fgsea",     # FIXED: Use "fgsea" (valid option) instead of "fgseaMultilevel"
            BPPARAM = BiocParallel::SerialParam()  # Disable parallel processing
          )
        }, error = function(e) {
          message(paste("GSEA error:", e$message))
          
          # Try the enricher method instead
          message("Trying enricher instead...")
          
          # Select the top ranked genes (positive values)
          positive_genes <- names(kegg_ranked_list)[kegg_ranked_list > 0]
          message(sprintf("   runKeggEnrichment: Found %d positively ranked genes", length(positive_genes)))
          
          if (length(positive_genes) > 0) {
            tryCatch({
              message("   runKeggEnrichment Step: Running enricher analysis...")
              enricher_result <- clusterProfiler::enricher(
                gene = positive_genes,
                universe = names(kegg_ranked_list),
                TERM2GENE = term2gene,
                TERM2NAME = term2name,
                pvalueCutoff = 0.05,  # Using higher threshold (0.1) for more lenient filtering
                pAdjustMethod = "fdr",
                minGSSize = 3,
                maxGSSize = 500
              )
              
              if (!is.null(enricher_result) && nrow(enricher_result@result) > 0) {
                message(sprintf("   runKeggEnrichment: Enricher analysis successful with %d results", nrow(enricher_result@result)))
                return(enricher_result)
              } else {
                message("   runKeggEnrichment: Enricher analysis returned no significant results")
              }
            }, error = function(e2) {
              message(paste("Enricher error:", e2$message))
              return(NULL)
            })
          }
          
          return(NULL)
        })
        
        # Format results for visualization
        if (!is.null(gsea_result)) {
          message("   runKeggEnrichment Step: Formatting results...")
          
          if (class(gsea_result)[1] == "enrichResult") {
            # Result from enricher
            message("   runKeggEnrichment: Processing results from enricher")
            if (nrow(gsea_result@result) > 0) {
              message(sprintf("   runKeggEnrichment: Found %d enriched terms", nrow(gsea_result@result)))
            } else {
              message("   runKeggEnrichment: No enriched terms found")
            }
            
            kegg_results <- gsea_result@result |>
              as_tibble() |>
              dplyr::select(
                termDescription = Description,
                enrichmentScore = p.adjust,  # Use p.adjust as score
                falseDiscoveryRate = p.adjust,
                genesMapped = Count,
                mappedIDs = geneID  # Add the list of mapped IDs
              ) |>
              mutate(
                enrichmentScore = -log10(enrichmentScore),  # Convert to -log10 scale
                comparison = assay_name
              )
          } else if (class(gsea_result)[1] == "gseaResult" && nrow(gsea_result@result) > 0) {
            # Result from GSEA
            message("   runKeggEnrichment: Processing results from GSEA")
            message(sprintf("   runKeggEnrichment: Found %d enriched terms", nrow(gsea_result@result)))
            
            kegg_results <- gsea_result@result |>
              as_tibble() |>
              dplyr::select(
                termDescription = Description,
                enrichmentScore = NES, 
                falseDiscoveryRate = p.adjust,
                genesMapped = setSize,
                mappedIDs = core_enrichment  # Add the core enrichment metabolites
              ) |>
              mutate(
                enrichmentScore = abs(enrichmentScore),
                comparison = assay_name
              )
          } else {
            message("   runKeggEnrichment: No enriched terms in GSEA/enricher output, returning empty table")
            kegg_results <- tibble(
              termDescription = character(),
              enrichmentScore = numeric(),
              falseDiscoveryRate = numeric(),
              genesMapped = integer(),
              mappedIDs = character(),  # Add empty column for mapped IDs
              comparison = character()
            )
          }
        } else {
          message("   runKeggEnrichment: GSEA/enricher analysis failed or returned NULL, returning empty table")
          kegg_results <- tibble(
            termDescription = character(),
            enrichmentScore = numeric(),
            falseDiscoveryRate = numeric(),
            genesMapped = integer(),
            mappedIDs = character(),  # Add empty column for mapped IDs
            comparison = character()
          )
        }
        
        message(sprintf("--- Exiting runKeggEnrichment. Returning table with %d rows ---", nrow(kegg_results)))
        return(kegg_results)
      } else {
        message("No KEGG IDs could be mapped from ChEBI IDs. Returning empty result.")
        return(tibble(
          termDescription = character(),
          enrichmentScore = numeric(),
          falseDiscoveryRate = numeric(),
          genesMapped = integer(),
          comparison = character()
        ))
      }
    } else {
      message("No ChEBI-KEGG mappings found for our ranked list. Returning empty result.")
      return(tibble(
        termDescription = character(),
        enrichmentScore = numeric(),
        falseDiscoveryRate = numeric(),
        genesMapped = integer(),
        comparison = character()
      ))
    }
  } else {
    # Use direct KEGG IDs from the ranked list
    message("   runKeggEnrichment Step: Using direct KEGG IDs from ranked list")
    kegg_ranked_list <- ranked_list[kegg_ids_in_ranked_list]
    
    # Sort in decreasing order
    kegg_ranked_list <- sort(kegg_ranked_list, decreasing = TRUE)
    
    # Diagnostic: Check for ties in the data
    message("   runKeggEnrichment DEBUG: Analyzing score distribution...")
    # Count exact duplicates in values
    duplicate_count <- sum(duplicated(kegg_ranked_list))
    duplicate_percentage <- (duplicate_count / length(kegg_ranked_list)) * 100
    message(sprintf("   runKeggEnrichment DEBUG: Found %d exact duplicate values (%.2f%% of total)",
                   duplicate_count, duplicate_percentage))
    
    # Show basic statistics
    message(sprintf("   runKeggEnrichment DEBUG: Score range: [%.4f, %.4f], Median: %.4f", 
                   min(kegg_ranked_list), max(kegg_ranked_list), median(kegg_ranked_list)))
    
    # Add small random noise to break ties in the ranked list
    message("   runKeggEnrichment Step: Adding small random noise to break ties in ranked list...")
    set.seed(42) # For reproducibility
    original_kegg_ranked_list <- kegg_ranked_list # Store original for reference
    kegg_ranked_list <- kegg_ranked_list + runif(length(kegg_ranked_list), min=-0.0001, max=0.0001)
    # Re-sort to ensure order is preserved with the new values
    kegg_ranked_list <- sort(kegg_ranked_list, decreasing = TRUE)
    message(sprintf("   runKeggEnrichment: Added jitter to %d values to prevent ties", length(kegg_ranked_list)))
    
    # Create TERM2GENE mapping
    message("   runKeggEnrichment Step: Creating TERM2GENE mapping with direct KEGG IDs...")
    term2gene <- species_specific_pathway_to_compound_tbl |>
      dplyr::filter(compound %in% names(kegg_ranked_list)) |>
      dplyr::select(pathway, compound) |>
      dplyr::rename(term = pathway, gene = compound)
    
    message(paste("Created TERM2GENE mapping with", nrow(term2gene), "entries"))
    if (nrow(term2gene) > 0) {
      message("   runKeggEnrichment Step: Sample term2gene mappings:")
      print(head(term2gene, 5))
    }
  }
  
  # Check if we have any mappings
  if (!exists("term2gene") || nrow(term2gene) == 0) {
    message("No TERM2GENE mappings created. Returning empty result.")
    return(tibble(
      termDescription = character(),
      enrichmentScore = numeric(),
      falseDiscoveryRate = numeric(),
      genesMapped = integer(),
      comparison = character()
    ))
  }
  
  # Get pathway ID to name mapping
  # First get all pathways from KEGG
  message("   runKeggEnrichment Step: Fetching pathways from KEGG...")
  pathways_tbl <- KEGGREST::keggList("pathway")
  
  message(sprintf("   runKeggEnrichment: Retrieved %d pathways from KEGG", length(pathways_tbl)))
  if (length(pathways_tbl) > 0) {
    message("Sample pathway entries:")
    print(head(pathways_tbl, 3))
  }
  
  # Create mapping from pathway ID to name
  message("   runKeggEnrichment Step: Creating pathway ID to name mapping...")
  pathway_id_to_name <- data.frame(
    pathway_id = names(pathways_tbl), 
    pathway_name = pathways_tbl
  ) |>
    dplyr::mutate(
      hsa_pathway_id = str_replace(pathway_id, "^path:map", "") |> 
        purrr::map_chr(\(x){paste0(kegg_species_code, x)})
    )
  
  message("   runKeggEnrichment Step: Sample pathway ID to name mapping:")
  print(head(pathway_id_to_name, 3))
  
  # Create TERM2NAME mapping
  message("   runKeggEnrichment Step: Creating TERM2NAME mapping...")
  term2name_initial <- term2gene |>
    distinct(term)
  
  message(sprintf("   runKeggEnrichment: Found %d distinct terms", nrow(term2name_initial)))
  
  # Debug the term format
  message("   runKeggEnrichment DEBUG: Sample term formats:")
  print(head(term2name_initial$term))
  
  # Approach 1: Direct matching with pathway_id
  pathway_id_to_name_extended <- pathway_id_to_name |>
    # Add additional format columns for matching
    mutate(
      term_format1 = paste0("path:", kegg_species_code, str_replace(hsa_pathway_id, paste0("^", kegg_species_code), "")),
      term_format2 = paste0("path:", hsa_pathway_id)
    )
  
  message("   runKeggEnrichment DEBUG: Extended pathway_id_to_name formats:")
  print(head(pathway_id_to_name_extended))
  
  # Try to join using the extended formats
  term2name <- term2name_initial |>
    left_join(pathway_id_to_name_extended, 
              by = c("term" = "term_format1")) |>
    dplyr::select(term, pathway_name)
  
  # If no results, try the second format
  if (all(is.na(term2name$pathway_name))) {
    message("   runKeggEnrichment DEBUG: First join attempt failed, trying second format...")
    term2name <- term2name_initial |>
      left_join(pathway_id_to_name_extended, 
                by = c("term" = "term_format2")) |>
      dplyr::select(term, pathway_name)
  }
  
  # If both failed, try extracting and joining by numeric part
  if (all(is.na(term2name$pathway_name))) {
    message("   runKeggEnrichment DEBUG: Second join attempt failed, trying numeric part matching...")
    
    # Extract numeric parts from terms
    term2name <- term2name_initial |>
      mutate(
        # Extract just the numeric part from the pathway ID
        pathway_num = str_extract(term, "\\d+")
      )
    
    message("   runKeggEnrichment DEBUG: Extracted pathway numbers:")
    print(head(term2name))
    
    # Also extract numeric parts from pathway IDs
    pathway_id_to_name_numeric <- pathway_id_to_name |>
      mutate(
        pathway_num = str_extract(pathway_id, "\\d+")
      )
    
    message("   runKeggEnrichment DEBUG: Pathway ID numeric parts:")
    print(head(pathway_id_to_name_numeric))
    
    # Join on numeric parts
    term2name <- term2name |>
      left_join(pathway_id_to_name_numeric, by = "pathway_num") |>
      dplyr::select(term, name = pathway_name)
  } else {
    # We had a successful join earlier, just rename
    term2name <- term2name |>
      dplyr::rename(name = pathway_name)
  }
  
  message("   runKeggEnrichment Step: After joining with pathway names:")
  print(head(term2name))
  
  term2name <- term2name |>
    dplyr::filter(!is.na(name))  # Make sure we have names
  
  # Try running GSEA with tryCatch to handle errors
  message("   runKeggEnrichment Step: Running GSEA analysis...")
  gsea_result <- tryCatch({
    # FIXED: Use fgsea and remove nPerm parameter
    clusterProfiler::GSEA(
      geneList = kegg_ranked_list,
      TERM2GENE = term2gene,
      TERM2NAME = term2name,
      minGSSize = 3,
      maxGSSize = 500,
      pvalueCutoff = 0.05,  # MODIFIED: Increased from 0.1 to 0.05 for more lenient filtering
      pAdjustMethod = "fdr",
      verbose = TRUE,
      seed = TRUE,      # Set seed for reproducibility
      by = "fgsea",     # FIXED: Use "fgsea" (valid option) instead of "fgseaMultilevel"
      BPPARAM = BiocParallel::SerialParam()  # Disable parallel processing
    )
  }, error = function(e) {
    message(paste("GSEA error:", e$message))
    
    # Try the enricher method instead
    message("Trying enricher instead...")
    
    # Select the top ranked genes (positive values)
    positive_genes <- names(kegg_ranked_list)[kegg_ranked_list > 0]
    message(sprintf("   runKeggEnrichment: Found %d positively ranked genes", length(positive_genes)))
    
    if (length(positive_genes) > 0) {
      tryCatch({
        message("   runKeggEnrichment Step: Running enricher analysis...")
        enricher_result <- clusterProfiler::enricher(
          gene = positive_genes,
          universe = names(kegg_ranked_list),
          TERM2GENE = term2gene,
          TERM2NAME = term2name,
          pvalueCutoff = 0.05,  # Using higher threshold (0.1) for more lenient filtering
          pAdjustMethod = "fdr",
          minGSSize = 3,
          maxGSSize = 500
        )
        
        if (!is.null(enricher_result) && nrow(enricher_result@result) > 0) {
          message(sprintf("   runKeggEnrichment: Enricher analysis successful with %d results", nrow(enricher_result@result)))
          return(enricher_result)
        } else {
          message("   runKeggEnrichment: Enricher analysis returned no significant results")
        }
      }, error = function(e2) {
        message(paste("Enricher error:", e2$message))
        return(NULL)
      })
    }
    
    return(NULL)
  })
  
  # Format results for visualization
  if (!is.null(gsea_result)) {
    message("   runKeggEnrichment Step: Formatting results...")
    
    if (class(gsea_result)[1] == "enrichResult") {
      # Result from enricher
      message("   runKeggEnrichment: Processing results from enricher")
      if (nrow(gsea_result@result) > 0) {
        message(sprintf("   runKeggEnrichment: Found %d enriched terms", nrow(gsea_result@result)))
      } else {
        message("   runKeggEnrichment: No enriched terms found")
      }
      
      kegg_results <- gsea_result@result |>
        as_tibble() |>
        dplyr::select(
          termDescription = Description,
          enrichmentScore = p.adjust,  # Use p.adjust as score
          falseDiscoveryRate = p.adjust,
          genesMapped = Count,
          mappedIDs = geneID  # Add the list of mapped IDs
        ) |>
        mutate(
          enrichmentScore = -log10(enrichmentScore),  # Convert to -log10 scale
          comparison = assay_name
        )
    } else if (class(gsea_result)[1] == "gseaResult" && nrow(gsea_result@result) > 0) {
      # Result from GSEA
      message("   runKeggEnrichment: Processing results from GSEA")
      message(sprintf("   runKeggEnrichment: Found %d enriched terms", nrow(gsea_result@result)))
      
      kegg_results <- gsea_result@result |>
        as_tibble() |>
        dplyr::select(
          termDescription = Description,
          enrichmentScore = NES, 
          falseDiscoveryRate = p.adjust,
          genesMapped = setSize,
          mappedIDs = core_enrichment  # Add the core enrichment metabolites
        ) |>
        mutate(
          enrichmentScore = abs(enrichmentScore),
          comparison = assay_name
        )
    } else {
      message("   runKeggEnrichment: No enriched terms in GSEA/enricher output, returning empty table")
      kegg_results <- tibble(
        termDescription = character(),
        enrichmentScore = numeric(),
        falseDiscoveryRate = numeric(),
        genesMapped = integer(),
        mappedIDs = character(),  # Add empty column for mapped IDs
        comparison = character()
      )
    }
  } else {
    message("   runKeggEnrichment: GSEA/enricher analysis failed or returned NULL, returning empty table")
    kegg_results <- tibble(
      termDescription = character(),
      enrichmentScore = numeric(),
      falseDiscoveryRate = numeric(),
      genesMapped = integer(),
      mappedIDs = character(),  # Add empty column for mapped IDs
      comparison = character()
    )
  }
  
  message(sprintf("--- Exiting runKeggEnrichment. Returning table with %d rows ---", nrow(kegg_results)))
  return(kegg_results)
}

#' Run Reactome Pathway Enrichment Analysis for Metabolites
#'
#' @param ranked_list Named numeric vector with ChEBI IDs as names and values as weights
#' @param mapping_table Metabolite ID mapping table
#' @param project_dirs Project directories structure
#' @param omic_type Omics type for directory access
#' @param experiment_label Experiment label for directory access
#' @param assay_name Name of the assay
#' @param reactome_organism Organism to use for Reactome (e.g., "Klebsiella pneumoniae")
#' @return A data frame with Reactome pathway enrichment results
#' @export
runReactomeEnrichment <- function(ranked_list, 
                                 mapping_table, 
                                 project_dirs,
                                 omic_type,
                                 experiment_label,
                                 assay_name,
                                 reactome_organism = NULL) {
  
  message(sprintf("--- Entering runReactomeEnrichment ---"))
  message(sprintf("   runReactomeEnrichment Args: assay_name = %s", assay_name))
  message(sprintf("   runReactomeEnrichment Args: reactome_organism = %s", as.character(reactome_organism)))
  message(sprintf("   runReactomeEnrichment Args: ranked_list length = %d", length(ranked_list)))
  
  # Use exact directory paths with paste0 construction
  data_dir <- project_dirs[[paste0(omic_type, "_", experiment_label)]]$data_dir
  message(sprintf("   runReactomeEnrichment: Using data_dir = %s", data_dir))
  
  # Load ChEBI to Reactome mapping
  tryCatch({
    message(sprintf("   runReactomeEnrichment Step: Loading ChEBI to Reactome mapping from %s...", file.path(data_dir, "chebi_to_reactome.tab")))
    
    # Check if file exists
    if (!file.exists(file.path(data_dir, "chebi_to_reactome.tab"))) {
      message("   runReactomeEnrichment ERROR: File 'chebi_to_reactome.tab' does not exist!")
      return(tibble(
        termDescription = character(),
        enrichmentScore = numeric(),
        falseDiscoveryRate = numeric(),
        genesMapped = integer(),
        comparison = character()
      ))
    }
    
    chebi_to_reactome <- vroom::vroom(
      file.path(data_dir, "chebi_to_reactome.tab"),
      show_col_types = FALSE
    ) |>
      dplyr::select(chebi_id, reactome_id, reactome_term, organism)
    
    message(sprintf("   runReactomeEnrichment: Loaded %d rows from chebi_to_reactome.tab", nrow(chebi_to_reactome)))
    message("   runReactomeEnrichment Step: Sample of ChEBI to Reactome mapping:")
    print(head(chebi_to_reactome, 5))
    
    # List available organisms in the data
    unique_organisms <- unique(chebi_to_reactome$organism)
    message(paste("Available organisms in Reactome data:", paste(unique_organisms, collapse=", ")))
    
    # FIXED: Improved organism handling with flexible matching
    if (!is.null(reactome_organism)) {
      # Try exact match first
      if (reactome_organism %in% unique_organisms) {
        message(paste("Filtering Reactome data for organism:", reactome_organism))
        chebi_to_reactome <- chebi_to_reactome |> 
          dplyr::filter(organism == reactome_organism)
      } else {
        # Try partial matching if exact match fails
        message(sprintf("   runReactomeEnrichment DEBUG: Trying partial match for '%s'", reactome_organism))
        possible_matches <- unique_organisms[grepl(reactome_organism, unique_organisms, ignore.case = TRUE)]
        if (length(possible_matches) > 0) {
          selected_organism <- possible_matches[1]  # Use first match
          message(paste("Using closest matching organism:", selected_organism))
          chebi_to_reactome <- chebi_to_reactome |> 
            dplyr::filter(organism == selected_organism)
        } else {
          # Fallback to Human if no match
          message("No matching organism found. Using 'Homo sapiens' as fallback...")
          chebi_to_reactome <- chebi_to_reactome |> 
            dplyr::filter(organism == "Homo sapiens")
        }
      }
      message(sprintf("   runReactomeEnrichment: After organism filtering, %d rows remain", nrow(chebi_to_reactome)))
    } else {
      message("Using all organisms in Reactome data (no filtering)")
    }
    
    # Extract just the numeric part from ChEBI IDs to match our format
    message("   runReactomeEnrichment Step: Formatting ChEBI IDs...")
    modified_chebi_to_reactome <- chebi_to_reactome |>
      dplyr::mutate(
        numeric_id = as.numeric(str_extract(chebi_id, "\\d+")),
        new_chebi_id = paste0("CHEBI:", numeric_id)
      )
    
    message("   runReactomeEnrichment Step: Sample after numeric ID extraction:")
    print(head(modified_chebi_to_reactome, 5))
    
    chebi_to_reactome <- modified_chebi_to_reactome |>
      dplyr::select(-chebi_id) |>
      dplyr::rename(chebi_id = new_chebi_id) |>
      dplyr::filter(!is.na(numeric_id))  # Remove any that didn't convert properly
    
    message(sprintf("   runReactomeEnrichment: After ChEBI ID formatting, %d rows remain", nrow(chebi_to_reactome)))
    message("   runReactomeEnrichment Step: Sample of formatted ChEBI IDs:")
    print(head(chebi_to_reactome, 5))
    
    # Check which ChEBI IDs overlap with our ranked list
    message("   runReactomeEnrichment Step: Checking for ChEBI IDs that overlap with ranked list...")
    chebi_ids_in_ranked_list <- intersect(chebi_to_reactome$chebi_id, names(ranked_list))
    message(paste("Found", length(chebi_ids_in_ranked_list), "metabolites that overlap with Reactome pathways."))
    
    if (length(chebi_ids_in_ranked_list) > 0) {
      message("   runReactomeEnrichment Step: Sample matching ChEBI IDs:")
      message(paste(head(chebi_ids_in_ranked_list, 5), collapse=", "))
    }
    
    if (length(chebi_ids_in_ranked_list) == 0) {
      message("No overlap between ChEBI IDs in ranked list and Reactome database. Returning empty result.")
      return(tibble(
        termDescription = character(),
        enrichmentScore = numeric(),
        falseDiscoveryRate = numeric(),
        genesMapped = integer(),
        comparison = character()
      ))
    }
    
    # Create a subset of the ranked list with only IDs that are in Reactome
    message("   runReactomeEnrichment Step: Creating filtered ranked list...")
    reactome_ranked_list <- ranked_list[chebi_ids_in_ranked_list]
    
    message(sprintf("   runReactomeEnrichment: Created ranked list with %d metabolites", length(reactome_ranked_list)))
    message("   runReactomeEnrichment Step: Sample values from ranked list:")
    print(head(reactome_ranked_list))
    
    # Sort the ranked list in decreasing order 
    reactome_ranked_list <- sort(reactome_ranked_list, decreasing = TRUE)
    
    # Diagnostic: Check for ties in the data
    message("   runReactomeEnrichment DEBUG: Analyzing score distribution...")
    # Count exact duplicates in values
    duplicate_count <- sum(duplicated(reactome_ranked_list))
    duplicate_percentage <- (duplicate_count / length(reactome_ranked_list)) * 100
    message(sprintf("   runReactomeEnrichment DEBUG: Found %d exact duplicate values (%.2f%% of total)",
                   duplicate_count, duplicate_percentage))
    
    # Show basic statistics
    message(sprintf("   runReactomeEnrichment DEBUG: Score range: [%.4f, %.4f], Median: %.4f", 
                   min(reactome_ranked_list), max(reactome_ranked_list), median(reactome_ranked_list)))
    
    # Add small random noise to break ties in the ranked list
    message("   runReactomeEnrichment Step: Adding small random noise to break ties in ranked list...")
    set.seed(42) # For reproducibility
    original_reactome_ranked_list <- reactome_ranked_list # Store original for reference
    reactome_ranked_list <- reactome_ranked_list + runif(length(reactome_ranked_list), min=-0.0001, max=0.0001)
    # Re-sort to ensure order is preserved with the new values
    reactome_ranked_list <- sort(reactome_ranked_list, decreasing = TRUE)
    message(sprintf("   runReactomeEnrichment: Added jitter to %d values to prevent ties", length(reactome_ranked_list)))
    
    # CRITICAL: Verify the ordering and distribution before GSEA
    message("   runReactomeEnrichment CRITICAL CHECK: Verifying ranked list ordering and distribution...")
    positive_values <- sum(reactome_ranked_list > 0)
    negative_values <- sum(reactome_ranked_list < 0)
    zero_values <- sum(reactome_ranked_list == 0)
    message(sprintf("   runReactomeEnrichment: Value distribution: %d positive, %d negative, %d zero", 
                   positive_values, negative_values, zero_values))
    
    # Check that the list is actually sorted (decreasing)
    is_decreasing <- all(diff(reactome_ranked_list) <= 0)
    message(sprintf("   runReactomeEnrichment: List is properly sorted in decreasing order: %s", 
                   if(is_decreasing) "YES" else "NO - THIS IS A PROBLEM!"))
    
    # Show top and bottom values to verify ordering
    message("   runReactomeEnrichment: Top 5 entries in ranked list:")
    for (i in 1:min(5, length(reactome_ranked_list))) {
      message(sprintf("     %s: %.6f", names(reactome_ranked_list)[i], reactome_ranked_list[i]))
    }
    message("   runReactomeEnrichment: Bottom 5 entries in ranked list:")
    for (i in (length(reactome_ranked_list) - min(4, length(reactome_ranked_list) - 1)):length(reactome_ranked_list)) {
      message(sprintf("     %s: %.6f", names(reactome_ranked_list)[i], reactome_ranked_list[i]))
    }
    
    # Create TERM2GENE mapping
    message("   runReactomeEnrichment Step: Creating TERM2GENE mapping...")
    term2gene <- chebi_to_reactome |>
      dplyr::filter(chebi_id %in% names(reactome_ranked_list)) |>
      dplyr::select(reactome_id, chebi_id) |>
      dplyr::rename(term = reactome_id, gene = chebi_id)
    
    message(sprintf("   runReactomeEnrichment: Created TERM2GENE mapping with %d rows", nrow(term2gene)))
    message("   runReactomeEnrichment Step: Sample TERM2GENE mapping:")
    print(head(term2gene, 5))
    
    # Check if we have pathways
    if (nrow(term2gene) == 0) {
      message("   runReactomeEnrichment WARNING: No term-gene mappings could be created!")
      return(tibble(
        termDescription = character(),
        enrichmentScore = numeric(),
        falseDiscoveryRate = numeric(),
        genesMapped = integer(),
        comparison = character()
      ))
    }
    
    # Count unique terms and genes
    unique_terms <- unique(term2gene$term)
    unique_genes <- unique(term2gene$gene)
    message(sprintf("   runReactomeEnrichment: Mapping contains %d unique terms and %d unique genes", 
                   length(unique_terms), length(unique_genes)))
    
    # Create TERM2NAME mapping
    message("   runReactomeEnrichment Step: Creating TERM2NAME mapping...")
    term2name <- chebi_to_reactome |>
      dplyr::select(reactome_id, reactome_term) |>
      dplyr::distinct() |>  # FIXED: Added distinct here to avoid duplicates
      dplyr::rename(term = reactome_id, name = reactome_term)
    
    message(sprintf("   runReactomeEnrichment: Created TERM2NAME mapping with %d rows", nrow(term2name)))
    message("   runReactomeEnrichment Step: Sample TERM2NAME mapping:")
    print(head(term2name, 5))
    
    # Debug: check for duplicate pathways
    message("   runReactomeEnrichment DEBUG: Checking for duplicate pathway IDs...")
    duplicated_terms <- term2name$term[duplicated(term2name$term)]
    if (length(duplicated_terms) > 0) {
      message(sprintf("   runReactomeEnrichment WARNING: Found %d duplicate pathway IDs", length(duplicated_terms)))
      message("Sample duplicates:")
      for (dup_term in head(duplicated_terms, 3)) {
        dup_entries <- term2name |> dplyr::filter(term == dup_term)
        message(sprintf("   Term: %s has %d entries:", dup_term, nrow(dup_entries)))
        print(dup_entries)
      }
    } else {
      message("   runReactomeEnrichment DEBUG: No duplicate pathway IDs found in TERM2NAME")
    }
    
    # Add fallback in case term2name is empty
    if (!exists("term2name") || nrow(term2name) == 0) {
      message("   runReactomeEnrichment WARNING: term2name is missing or empty, creating fallback...")
      # Create a minimal term2name mapping using the term itself
      term2name <- term2gene |>
        dplyr::distinct(term) |>
        dplyr::mutate(name = term)
      
      message("   runReactomeEnrichment: Created fallback term2name mapping:")
      print(head(term2name, 3))
    }
    
    # Try to run GSEA analysis
    message("   runReactomeEnrichment Step: Running GSEA analysis...")
    
    # ABSOLUTE FINAL SAFETY CHECK - verify both term2gene and term2name exist
    if (!exists("term2gene") || nrow(term2gene) == 0) {
      message("   runReactomeEnrichment CRITICAL ERROR: term2gene is missing or empty before GSEA call")
      # Can't proceed without term2gene
      return(tibble(
        termDescription = character(),
        enrichmentScore = numeric(),
        falseDiscoveryRate = numeric(),
        genesMapped = integer(),
        comparison = character()
      ))
    }
    
    if (!exists("term2name") || nrow(term2name) == 0) {
      message("   runReactomeEnrichment WARNING: term2name is missing or empty before GSEA call, creating emergency fallback")
      # Create emergency fallback
      term2name <- term2gene |> 
        dplyr::distinct(term) |>
        dplyr::mutate(name = term)
    }
    
    # Debug output
    message(sprintf("   runReactomeEnrichment VERIFICATION: term2gene has %d rows, term2name has %d rows", 
                   nrow(term2gene), nrow(term2name)))
    
    gsea_result <- tryCatch({
      clusterProfiler::GSEA(
        geneList = reactome_ranked_list,
        TERM2GENE = term2gene,
        TERM2NAME = term2name,
        minGSSize = 3,
        maxGSSize = 500,
        pvalueCutoff = 0.05,  # MODIFIED: Increased from 0.1 to 0.05 for more lenient filtering
        pAdjustMethod = "fdr",
        verbose = TRUE,
        seed = TRUE,      # Set seed for reproducibility
        by = "fgsea",     # FIXED: Use "fgsea" (valid option) instead of "fgseaMultilevel"
        BPPARAM = BiocParallel::SerialParam()  # Disable parallel processing
      )
    }, error = function(e) {
      message(paste("GSEA error:", e$message))
      
      # Try the enricher method instead
      message("Trying enricher instead...")
      
      # Select the top ranked genes (positive values)
      positive_genes <- names(reactome_ranked_list)[reactome_ranked_list > 0]
      message(sprintf("   runReactomeEnrichment: Found %d positively ranked genes", length(positive_genes)))
      
      if (length(positive_genes) > 0) {
        tryCatch({
          message("   runReactomeEnrichment Step: Running enricher analysis...")
          message("   runReactomeEnrichment DEBUG: Top 10 genes for enricher:")
          message(paste(head(positive_genes, 10), collapse=", "))
          
          # Try with higher p-value cutoff
          message("   runReactomeEnrichment Step: Using higher p-value cutoff (0.1)...")
          enricher_result <- clusterProfiler::enricher(
            gene = positive_genes,
            universe = names(reactome_ranked_list),
            TERM2GENE = term2gene,
            TERM2NAME = term2name,
            pvalueCutoff = 0.05,  # Using higher threshold (0.1) for more enriched pathways
            pAdjustMethod = "fdr",
            minGSSize = 3,
            maxGSSize = 500
          )
          
          if (!is.null(enricher_result) && nrow(enricher_result@result) > 0) {
            message(sprintf("   runReactomeEnrichment: Enricher analysis successful with %d results", nrow(enricher_result@result)))
            return(enricher_result)
          } else {
            message("   runReactomeEnrichment: Enricher analysis returned no significant results")
          }
        }, error = function(e2) {
          message(paste("Enricher error:", e2$message))
          return(NULL)
        })
      }
      
      return(NULL)
    })
    
    # Format results for visualization
    if (!is.null(gsea_result)) {
      message("   runReactomeEnrichment Step: Formatting results...")
      
      if (class(gsea_result)[1] == "enrichResult") {
        # Result from enricher
        message("   runReactomeEnrichment: Processing results from enricher")
        if (nrow(gsea_result@result) > 0) {
          message(sprintf("   runReactomeEnrichment: Found %d enriched terms", nrow(gsea_result@result)))
        } else {
          message("   runReactomeEnrichment: No enriched terms found")
        }
        
        reactome_results <- gsea_result@result |>
          as_tibble() |>
          dplyr::select(
            termDescription = Description,
            enrichmentScore = p.adjust,  # Use p.adjust as score
            falseDiscoveryRate = p.adjust,
            genesMapped = Count,
            mappedIDs = geneID  # Add the list of mapped IDs
          ) |>
          mutate(
            enrichmentScore = -log10(enrichmentScore),  # Convert to -log10 scale
            comparison = assay_name
          ) |>
          distinct()  # FIXED: Added distinct to remove duplicates
      } else if (class(gsea_result)[1] == "gseaResult" && nrow(gsea_result@result) > 0) {
        # Result from GSEA
        message("   runReactomeEnrichment: Processing results from GSEA")
        message(sprintf("   runReactomeEnrichment: Found %d enriched terms", nrow(gsea_result@result)))
        
        reactome_results <- gsea_result@result |>
          as_tibble() |>
          dplyr::select(
            termDescription = Description,
            enrichmentScore = NES, 
            falseDiscoveryRate = p.adjust,
            genesMapped = setSize,
            mappedIDs = core_enrichment  # Add the core enrichment metabolites
          ) |>
          mutate(
            enrichmentScore = abs(enrichmentScore),
            comparison = assay_name
          ) |>
          distinct()  # FIXED: Added distinct to remove duplicates
      } else {
        message("   runReactomeEnrichment: No enriched terms in GSEA/enricher output, returning empty table")
        reactome_results <- tibble(
          termDescription = character(),
          enrichmentScore = numeric(),
          falseDiscoveryRate = numeric(),
          genesMapped = integer(),
          mappedIDs = character(),  # Add empty column for mapped IDs
          comparison = character()
        )
      }
      
      # FIXED: Check for duplicate results after processing
      if (nrow(reactome_results) > 0) {
        duplicate_count <- sum(duplicated(reactome_results$termDescription))
        if (duplicate_count > 0) {
          message(sprintf("   runReactomeEnrichment WARNING: Found %d duplicate pathway descriptions after processing", duplicate_count))
          message("   runReactomeEnrichment DEBUG: Applying additional distinct() operation...")
          reactome_results <- reactome_results |> distinct(termDescription, .keep_all = TRUE)
        }
      }
      
    } else {
      message("   runReactomeEnrichment: GSEA/enricher analysis failed or returned NULL, returning empty table")
      reactome_results <- tibble(
        termDescription = character(),
        enrichmentScore = numeric(),
        falseDiscoveryRate = numeric(),
        genesMapped = integer(),
        mappedIDs = character(),  # Add empty column for mapped IDs
        comparison = character()
      )
    }
    
    message(sprintf("--- Exiting runReactomeEnrichment. Returning table with %d rows ---", nrow(reactome_results)))
    return(reactome_results)
    
  }, error = function(e) {
    message(paste("Error processing Reactome data:", e$message))
    return(tibble(
      termDescription = character(),
      enrichmentScore = numeric(),
      falseDiscoveryRate = numeric(),
      genesMapped = integer(),
      comparison = character()
    ))
  })
}

#' Run Metabolomics Pathway Enrichment for Both Assays
#'
#' @param weights Data frame containing MOFA weights
#' @param metabolomics_obj Metabolomics S4 object containing ChEBI IDs
#' @param mapping_table The metabolite ID mapping table from AnnotationHub
#' @param project_dirs Project directories structure
#' @param omic_type Omics type for directory access
#' @param experiment_label Experiment label for directory access
#' @param kegg_species_code KEGG species code (e.g., "kpn" for Klebsiella pneumoniae)
#' @param reactome_organism Organism name to use for Reactome filtering (optional)
#' @return A combined data frame with enrichment results for both assays
#' @export
runMetabolomicsPathwayEnrichment <- function(weights, 
                                            metabolomics_obj, 
                                            mapping_table, 
                                            project_dirs,
                                            omic_type,
                                            experiment_label,
                                            kegg_species_code = "kpn",
                                            reactome_organism = NULL) {
  
  message(sprintf("--- Entering runMetabolomicsPathwayEnrichment ---"))
  message(sprintf("   runMetabolomicsPathwayEnrichment Args: kegg_species_code = %s", kegg_species_code))
  message(sprintf("   runMetabolomicsPathwayEnrichment Args: reactome_organism = %s", as.character(reactome_organism)))
  
  # Use exact directory paths with paste0 construction
  results_dir <- project_dirs[[paste0(omic_type, "_", experiment_label)]]$integration_enrichment_plots_dir
  
  # Inspect the structure of input objects
  message("   runMetabolomicsPathwayEnrichment Step: Inspecting metabolomics object structure...")
  message(sprintf("   runMetabolomicsPathwayEnrichment: Metabolomics object has %d assays", length(metabolomics_obj@metabolite_data)))
  message(sprintf("   runMetabolomicsPathwayEnrichment: Assay 1 (LC) has %d metabolites", nrow(metabolomics_obj@metabolite_data[[1]])))
  message(sprintf("   runMetabolomicsPathwayEnrichment: Assay 2 (GC) has %d metabolites", nrow(metabolomics_obj@metabolite_data[[2]])))
  
  message("   runMetabolomicsPathwayEnrichment Step: Inspecting weights data frame...")
  message(sprintf("   runMetabolomicsPathwayEnrichment: Weights data frame has %d entries", nrow(weights)))
  metabolome_lc_weights <- sum(weights$view == "metabolome_lc")
  metabolome_gc_weights <- sum(weights$view == "metabolome_gc")
  message(sprintf("   runMetabolomicsPathwayEnrichment: Found %d LC-MS weights and %d GC-MS weights", 
                 metabolome_lc_weights, metabolome_gc_weights))
  
  message("   runMetabolomicsPathwayEnrichment Step: Inspecting mapping table...")
  message(sprintf("   runMetabolomicsPathwayEnrichment: Mapping table has %d entries", nrow(mapping_table)))
  kegg_mappings <- sum(!is.na(mapping_table$KEGG))
  chebi_mappings <- sum(!is.na(mapping_table$ChEBI))
  message(sprintf("   runMetabolomicsPathwayEnrichment: Found %d KEGG IDs and %d ChEBI IDs in mapping table", 
                 kegg_mappings, chebi_mappings))
  
  # Run enrichment for LC-MS metabolomics
  message("   runMetabolomicsPathwayEnrichment Step: Running LC-MS metabolomics enrichment...")
  lc_results <- runMetabolomicsEnrichmentAnalysis(
    weights = weights,
    metabolomics_obj = metabolomics_obj,
    mapping_table = mapping_table,
    project_dirs = project_dirs,
    omic_type = omic_type,
    experiment_label = experiment_label,
    assay_name = "metabolome_lc",
    kegg_species_code = kegg_species_code,
    reactome_organism = reactome_organism
  )
  
  message(sprintf("   runMetabolomicsPathwayEnrichment: LC-MS analysis returned %d results", nrow(lc_results)))
  
  # Run enrichment for GC-MS metabolomics
  message("   runMetabolomicsPathwayEnrichment Step: Running GC-MS metabolomics enrichment...")
  gc_results <- runMetabolomicsEnrichmentAnalysis(
    weights = weights,
    metabolomics_obj = metabolomics_obj,
    mapping_table = mapping_table,
    project_dirs = project_dirs,
    omic_type = omic_type,
    experiment_label = experiment_label,
    assay_name = "metabolome_gc",
    kegg_species_code = kegg_species_code,
    reactome_organism = reactome_organism
  )
  
  message(sprintf("   runMetabolomicsPathwayEnrichment: GC-MS analysis returned %d results", nrow(gc_results)))
  
  # Combine results
  message("   runMetabolomicsPathwayEnrichment Step: Combining LC-MS and GC-MS results...")
  combined_results <- bind_rows(
    lc_results |> mutate(assay = "LC-MS"),
    gc_results |> mutate(assay = "GC-MS")
  ) |>
  mutate(comparison = "Metabolome") |> # Set overall comparison label
  distinct()  # Ensure no duplicates
  
  message(sprintf("   runMetabolomicsPathwayEnrichment: Combined results has %d entries", nrow(combined_results)))
  
  # Add original metabolite names
  message("   runMetabolomicsPathwayEnrichment Step: Adding original metabolite names...")
  
  # Create a mapping from KEGG IDs to metabolite names
  message("   runMetabolomicsPathwayEnrichment Step: Creating ID to name mapping tables...")
  
  # Extract all metabolite names and IDs from metabolomics object
  lc_names_mapping <- metabolomics_obj@metabolite_data[[1]] |>
    dplyr::select(metabolite, database_identifier) |>
    dplyr::filter(!is.na(database_identifier)) |>
    dplyr::mutate(chebi_id = stringr::str_extract(database_identifier, "CHEBI:\\d+"))
    
  gc_names_mapping <- metabolomics_obj@metabolite_data[[2]] |>
    dplyr::select(metabolite, database_identifier) |>
    dplyr::filter(!is.na(database_identifier)) |>
    dplyr::mutate(chebi_id = stringr::str_extract(database_identifier, "CHEBI:\\d+"))
    
  # Combine mappings from both assays
  all_names_mapping <- bind_rows(
    lc_names_mapping |> mutate(assay = "LC-MS"),
    gc_names_mapping |> mutate(assay = "GC-MS")
  )
  
  message(sprintf("   runMetabolomicsPathwayEnrichment: Created mapping table with %d entries", nrow(all_names_mapping)))
  
  # Also get KEGG to ChEBI mapping for translating KEGG IDs
  kegg_to_chebi <- mapping_table |>
    dplyr::filter(!is.na(KEGG) & !is.na(ChEBI)) |>
    dplyr::select(KEGG, ChEBI) |>
    dplyr::mutate(
      kegg_id = paste0("cpd:", KEGG),
      chebi_id = paste0("CHEBI:", ChEBI)
    )
  
  # Function to map IDs to names
  map_ids_to_names <- function(id_string, assay_type) {
    if (is.na(id_string) || id_string == "") {
      return(NA_character_)
    }
    
    # Split the ID string
    ids <- unlist(strsplit(id_string, split = ",|/"))
    names_vec <- character(length(ids))
    
    for (i in seq_along(ids)) {
      id <- ids[i]
      if (grepl("^cpd:", id)) {
        # For KEGG IDs, try to map to ChEBI first
        kegg_matches <- kegg_to_chebi |> dplyr::filter(kegg_id == id)
        if (nrow(kegg_matches) > 0) {
          # Found a matching ChEBI ID, now look up its name
          chebi_id <- kegg_matches$chebi_id[1]
          chebi_matches <- all_names_mapping |> 
            dplyr::filter(chebi_id == !!chebi_id & assay == !!assay_type)
          
          if (nrow(chebi_matches) > 0) {
            names_vec[i] <- chebi_matches$metabolite[1]
          } else {
            # Try without assay filter as fallback
            chebi_matches <- all_names_mapping |> dplyr::filter(chebi_id == !!chebi_id)
            if (nrow(chebi_matches) > 0) {
              names_vec[i] <- chebi_matches$metabolite[1]
            } else {
              names_vec[i] <- id  # Just use the ID if no match
            }
          }
        } else {
          names_vec[i] <- id  # Just use the ID if no match
        }
      } else if (grepl("^CHEBI:", id)) {
        # Direct ChEBI ID lookup
        chebi_matches <- all_names_mapping |> 
          dplyr::filter(chebi_id == !!id & assay == !!assay_type)
        
        if (nrow(chebi_matches) > 0) {
          names_vec[i] <- chebi_matches$metabolite[1]
        } else {
          # Try without assay filter as fallback
          chebi_matches <- all_names_mapping |> dplyr::filter(chebi_id == !!id)
          if (nrow(chebi_matches) > 0) {
            names_vec[i] <- chebi_matches$metabolite[1]
          } else {
            names_vec[i] <- id  # Just use the ID if no match
          }
        }
      } else {
        names_vec[i] <- id  # Just use the ID if no match
      }
    }
    
    # Join with commas
    paste(names_vec, collapse = ", ")
  }
  
  # Apply the mapping function to each row
  message("   runMetabolomicsPathwayEnrichment Step: Mapping IDs to names...")
  combined_results <- combined_results |>
    dplyr::rowwise() |>
    dplyr::mutate(
      mappedNames = map_ids_to_names(mappedIDs, assay)
    ) |>
    dplyr::ungroup()
  
  # Save combined results
  message(sprintf("   runMetabolomicsPathwayEnrichment Step: Saving combined results to %s...", 
                 file.path(results_dir, "combined_metabolomics_enrichment_results.tab")))
  vroom::vroom_write(combined_results, 
                    file.path(results_dir, "combined_metabolomics_enrichment_results.tab"))
  
  message("--- Exiting runMetabolomicsPathwayEnrichment ---")
  return(combined_results)
}

#' Run STRING DB Enrichment Analysis from DE Results S4 Object
#'
#' @description
#' This function extracts differential expression data from a de_results_for_enrichment S4 object,
#' applies the specified ranking method, and performs STRING DB enrichment analysis using the
#' existing runOneStringDbRankEnrichmentMofa function.
#'
#' @param de_results_for_enrichment An S4 object of class de_results_for_enrichment containing
#'   differential expression results across multiple contrasts.
#' @param contrast_name Character string. Name of the specific contrast to analyze.
#'   Must match one of the names in de_results_for_enrichment@de_data.
#'   If NULL, will use the first contrast. Default: NULL.
#' @param ranking_method Character string. Method for ranking proteins. Options:
#'   - "fdr_qvalue": Rank by FDR q-value (ascending, most significant first)
#'   - "log2fc": Rank by log2 fold change (descending, highest FC first)  
#'   - "combined_score": Use sign(log2FC) * (-log10(fdr_qvalue)) for ranking
#'   - "none": No ranking applied (proteins in original order)
#'   Default: "combined_score".
#' @param identifier_column Character string. Name of the column containing protein identifiers.
#'   Default: "Protein.Ids".
#' @param filter_significant Logical. Whether to filter to only significant proteins (fdr_qvalue < 0.05).
#'   Default: FALSE (include all proteins).
#' @param fdr_threshold Numeric. FDR threshold for filtering significant proteins when filter_significant=TRUE.
#'   Default: 0.05.
#' @param result_label Character string. A label used for naming output files.
#'   If NULL, will use the contrast name. Default: NULL.
#' @param results_dir Character string. The path to the directory where enrichment
#'   results will be saved. Default: "string_enrichment_results".
#' @param api_key Character string. Your personal STRING API key.
#'   Default: NULL.
#' @param species Character string. NCBI/STRING species identifier.
#'   Default: "9606" (Homo sapiens).
#' @param ge_fdr Numeric. FDR threshold for gene expression enrichment.
#'   Default: 0.05.
#' @param ge_enrichment_rank_direction Integer. Direction for enrichment rank
#'   (-1, 0, or 1). Default: -1.
#' @param polling_interval_seconds Numeric. Seconds to wait between polling attempts.
#'   Default: 10.
#' @param max_polling_attempts Numeric. Maximum polling attempts before timing out.
#'   Default: 30.
#'
#' @return A data frame containing the enrichment results from STRING DB.
#'   Returns NULL if the process fails.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming you have a de_results_for_enrichment object
#' enrichment_results <- runStringDbEnrichmentFromDEResults(
#'   de_results_for_enrichment = my_de_results,
#'   contrast_name = "T2.minus.MSO=groupT2-groupMSO",
#'   ranking_method = "combined_score",
#'   filter_significant = FALSE,
#'   result_label = "T2_vs_MSO_enrichment",
#'   results_dir = "string_enrichment_output",
#'   api_key = "your_api_key",
#'   species = "9606"
#' )
#' }
runStringDbEnrichmentFromDEResults <- function(de_results_for_enrichment,
                                             contrast_name = NULL,
                                             ranking_method = "combined_score",
                                             identifier_column = "Protein.Ids",
                                             filter_significant = FALSE,
                                             fdr_threshold = 0.05,
                                             result_label = NULL,
                                             results_dir = "string_enrichment_results",
                                             api_key = NULL,
                                             species = "9606",
                                             ge_fdr = 0.05,
                                             ge_enrichment_rank_direction = -1,
                                             polling_interval_seconds = 10,
                                             max_polling_attempts = 30) {
  
  # Load required packages
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.")
  }
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package 'stringr' is required but not installed.")
  }
  
  # Validate input S4 object
  if (!inherits(de_results_for_enrichment, "de_results_for_enrichment")) {
    stop("Input must be an S4 object of class 'de_results_for_enrichment'")
  }
  
  # Get available contrasts
  available_contrasts <- names(de_results_for_enrichment@de_data)
  
  if (length(available_contrasts) == 0) {
    stop("No contrasts found in de_results_for_enrichment@de_data")
  }
  
  # Select contrast
  if (is.null(contrast_name)) {
    contrast_name <- available_contrasts[1]
    message(paste("No contrast specified. Using:", contrast_name))
  } else if (!contrast_name %in% available_contrasts) {
    stop(paste("Contrast", contrast_name, "not found. Available contrasts:", 
               paste(available_contrasts, collapse = ", ")))
  }
  
  # Extract data for the specified contrast
  de_data <- de_results_for_enrichment@de_data[[contrast_name]]
  
  if (is.null(de_data) || nrow(de_data) == 0) {
    stop(paste("No data found for contrast:", contrast_name))
  }
  
  # Validate required columns
  required_cols <- c(identifier_column, "fdr_qvalue", "log2FC")
  missing_cols <- required_cols[!required_cols %in% colnames(de_data)]
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Validate ranking method
  valid_methods <- c("fdr_qvalue", "log2fc", "combined_score", "none")
  if (!ranking_method %in% valid_methods) {
    stop(paste("Invalid ranking_method. Must be one of:", paste(valid_methods, collapse = ", ")))
  }
  
  # Filter significant proteins if requested
  if (filter_significant) {
    initial_count <- nrow(de_data)
    de_data <- de_data |>
      dplyr::filter(fdr_qvalue < fdr_threshold)
    final_count <- nrow(de_data)
    message(paste("Filtered from", initial_count, "to", final_count, 
                  "proteins using FDR threshold of", fdr_threshold))
    
    if (nrow(de_data) == 0) {
      stop("No proteins remain after filtering for significance")
    }
  }
  
  # Clean protein IDs (remove anything after ":")  
  de_data_processed <- de_data |>
    dplyr::mutate(
      protein_id = purrr::map_chr(!!sym(identifier_column), 
                                  ~stringr::str_split(.x, ":")[[1]][1])
    ) |>
    dplyr::filter(!is.na(protein_id), protein_id != "")
  
  # Apply ranking method and create score column
  if (ranking_method == "fdr_qvalue") {
    # Rank by FDR (ascending - most significant first)
    # Use negative log10 for proper ordering in STRING DB
    de_data_processed <- de_data_processed |>
      dplyr::arrange(fdr_qvalue) |>
      dplyr::mutate(score = -log10(fdr_qvalue + 1e-10))  # Add small value to avoid log(0)
    
  } else if (ranking_method == "log2fc") {
    # Rank by log2FC (descending - highest FC first)
    de_data_processed <- de_data_processed |>
      dplyr::arrange(desc(abs(log2FC))) |>
      dplyr::mutate(score = log2FC)
    
  } else if (ranking_method == "combined_score") {
    # Use sign(log2FC) * (-log10(fdr_qvalue))
    de_data_processed <- de_data_processed |>
      dplyr::mutate(score = sign(log2FC) * (-log10(fdr_qvalue + 1e-10))) |>
      dplyr::arrange(desc(abs(score)))
    
  } else if (ranking_method == "none") {
    # No ranking - use original order with a neutral score
    de_data_processed <- de_data_processed |>
      dplyr::mutate(score = 1)
  }
  
  # Prepare input table for STRING DB
  string_input_table <- de_data_processed |>
    dplyr::select(protein_id, score) |>
    dplyr::filter(!is.na(score), !is.infinite(score))
  
  if (nrow(string_input_table) == 0) {
    stop("No valid protein-score pairs remain after processing")
  }
  
  # Set result label
  if (is.null(result_label)) {
    result_label <- paste0(contrast_name, "_", ranking_method)
  }
  
  message(paste("Submitting", nrow(string_input_table), "proteins to STRING DB"))
  message(paste("Ranking method:", ranking_method))
  message(paste("Result label:", result_label))
  
  # Call the existing MOFA function
  enrichment_results <- runOneStringDbRankEnrichmentMofa(
    input_table = string_input_table,
    identifier_column_name = "protein_id",
    value_column_name = "score",
    result_label = result_label,
    results_dir = results_dir,
    api_key = api_key,
    species = species,
    ge_fdr = ge_fdr,
    ge_enrichment_rank_direction = ge_enrichment_rank_direction,
    polling_interval_seconds = polling_interval_seconds,
    max_polling_attempts = max_polling_attempts
  )
  
  return(enrichment_results)
}

#' Run STRING DB Enrichment Analysis for Multiple Contrasts from DE Results
#'
#' @description
#' This function runs STRING DB enrichment analysis for all or selected contrasts
#' from a de_results_for_enrichment S4 object.
#'
#' @param de_results_for_enrichment An S4 object of class de_results_for_enrichment.
#' @param contrast_names Character vector. Names of contrasts to analyze. 
#'   If NULL, analyzes all available contrasts. Default: NULL.
#' @param ranking_method Character string. Same options as runStringDbEnrichmentFromDEResults.
#'   Default: "combined_score".
#' @param ... Additional arguments passed to runStringDbEnrichmentFromDEResults.
#'
#' @return A named list of enrichment results, one for each contrast.
#'
#' @export
runStringDbEnrichmentFromDEResultsMultiple <- function(de_results_for_enrichment,
                                                      contrast_names = NULL,
                                                      ranking_method = "combined_score",
                                                      ...) {
  
  # Get available contrasts
  available_contrasts <- names(de_results_for_enrichment@de_data)
  
  if (is.null(contrast_names)) {
    contrast_names <- available_contrasts
  } else {
    # Validate contrast names
    invalid_contrasts <- contrast_names[!contrast_names %in% available_contrasts]
    if (length(invalid_contrasts) > 0) {
      stop(paste("Invalid contrast names:", paste(invalid_contrasts, collapse = ", ")))
    }
  }
  
  # Run enrichment for each contrast
  results_list <- purrr::map(contrast_names, function(contrast) {
    message(paste("Processing contrast:", contrast))
    
    tryCatch({
      runStringDbEnrichmentFromDEResults(
        de_results_for_enrichment = de_results_for_enrichment,
        contrast_name = contrast,
        ranking_method = ranking_method,
        result_label = paste0(contrast, "_", ranking_method),
        ...
      )
    }, error = function(e) {
      message(paste("Error processing contrast", contrast, ":", e$message))
      return(NULL)
    })
  })
  
  # Name the results list
  names(results_list) <- contrast_names
  
  return(results_list)
}

#' Get Available Species from STRING DB
#'
#' @description
#' This function queries the STRING DB API to retrieve all available species
#' and their corresponding species identifiers. This is useful for finding
#' the correct species code to use in STRING DB enrichment analyses.
#'
#' @param search_term Character string. Optional search term to filter species
#'   by name (case-insensitive partial matching). If NULL, returns all species.
#'   Default: NULL.
#' @param api_key Character string. Your STRING DB API key. Default: NULL.
#'
#' @return A data frame containing available species with columns:
#'   - species_id: STRING DB species identifier (e.g., "9606", "STRG0A62HCE")
#'   - official_name: Official species name
#'   - compact_name: Compact species name
#'   - taxon_id: NCBI taxonomy ID
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Get all available species
#' all_species <- getStringDbSpecies()
#' 
#' # Search for human species
#' human_species <- getStringDbSpecies("sapiens")
#' 
#' # Search for mouse species  
#' mouse_species <- getStringDbSpecies("musculus")
#' 
#' # Search for a specific genus
#' zingiber_species <- getStringDbSpecies("zingiber")
#' }
getStringDbSpecies <- function(search_term = NULL, api_key = NULL) {
  
  # Load required packages
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required but not installed.")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required but not installed.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.")
  }
  
  # STRING DB API endpoint for species
  api_url <- "https://version-12-0.string-db.org/api/json/get_string_ids"
  
  # Prepare parameters
  params <- list(format = "json")
  if (!is.null(api_key)) {
    params$api_key <- api_key
  }
  
  message("Querying STRING DB for available species...")
  
  # Make API request
  response <- tryCatch({
    httr::GET(url = "https://version-12-0.string-db.org/api/json/species", 
              query = params)
  }, error = function(e) {
    stop(paste("Failed to connect to STRING DB API:", e$message))
  })
  
  # Check if request was successful
  if (httr::http_error(response)) {
    stop(paste("STRING DB API request failed with status:", httr::status_code(response)))
  }
  
  # Parse JSON response
  content_text <- httr::content(response, "text", encoding = "UTF-8")
  species_data <- tryCatch({
    jsonlite::fromJSON(content_text, flatten = TRUE)
  }, error = function(e) {
    stop(paste("Failed to parse JSON response from STRING DB:", e$message))
  })
  
  # Convert to data frame if it's not already
  if (!is.data.frame(species_data)) {
    stop("Unexpected response format from STRING DB species API")
  }
  
  # Standardize column names (STRING DB API might have different column names)
  # Try to identify the key columns
  col_mapping <- list()
  
  for (col in colnames(species_data)) {
    lower_col <- tolower(col)
    if (grepl("species.*id|taxon.*id", lower_col)) {
      col_mapping$species_id <- col
    } else if (grepl("official.*name|species.*name", lower_col)) {
      col_mapping$official_name <- col
    } else if (grepl("compact.*name|short.*name", lower_col)) {
      col_mapping$compact_name <- col
    } else if (grepl("^taxon$|ncbi.*id", lower_col)) {
      col_mapping$taxon_id <- col
    }
  }
  
  # Create standardized data frame
  result_df <- species_data
  
  # Rename columns if mappings were found
  for (new_name in names(col_mapping)) {
    old_name <- col_mapping[[new_name]]
    if (old_name %in% colnames(result_df)) {
      result_df <- result_df |>
        dplyr::rename(!!new_name := !!old_name)
    }
  }
  
  # Ensure we have at least species_id and one name column
  if (!"species_id" %in% colnames(result_df)) {
    # Try to find any ID column
    id_cols <- colnames(result_df)[grepl("id", colnames(result_df), ignore.case = TRUE)]
    if (length(id_cols) > 0) {
      result_df <- result_df |>
        dplyr::rename(species_id = !!id_cols[1])
    } else {
      stop("Could not identify species ID column in STRING DB response")
    }
  }
  
  # Apply search filter if provided
  if (!is.null(search_term)) {
    message(paste("Filtering species with search term:", search_term))
    
    # Create a search pattern (case-insensitive)
    search_pattern <- paste0("(?i)", search_term)
    
    # Search across all text columns
    text_cols <- colnames(result_df)[sapply(result_df, is.character)]
    
    # Create a combined search field
    if (length(text_cols) > 0) {
      result_df <- result_df |>
        dplyr::mutate(
          search_field = paste(!!!syms(text_cols), sep = " ", collapse = " ")
        ) |>
        dplyr::filter(grepl(search_pattern, search_field)) |>
        dplyr::select(-search_field)
    }
  }
  
  # Sort by species_id for consistent output
  result_df <- result_df |>
    dplyr::arrange(species_id)
  
  message(paste("Found", nrow(result_df), "species matching your criteria"))
  
  return(result_df)
}

#' Search STRING DB Species by Name
#'
#' @description
#' A convenience function to search for STRING DB species by name.
#' This is a wrapper around getStringDbSpecies with better error handling
#' and formatted output for common use cases.
#'
#' @param species_name Character string. The species name to search for
#'   (case-insensitive partial matching).
#' @param api_key Character string. Your STRING DB API key. Default: NULL.
#' @param show_top_n Integer. Maximum number of results to display. Default: 10.
#'
#' @return A data frame with the most relevant species matches.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Search for human
#' searchStringDbSpecies("human")
#' searchStringDbSpecies("homo sapiens")
#' 
#' # Search for mouse
#' searchStringDbSpecies("mouse")
#' searchStringDbSpecies("mus musculus")
#' 
#' # Search for your organism
#' searchStringDbSpecies("zingiber")
#' searchStringDbSpecies("ginger")
#' }
searchStringDbSpecies <- function(species_name, api_key = NULL, show_top_n = 10) {
  
  if (missing(species_name) || is.null(species_name) || species_name == "") {
    stop("Please provide a species name to search for")
  }
  
  # Get species data
  species_df <- getStringDbSpecies(search_term = species_name, api_key = api_key)
  
  if (nrow(species_df) == 0) {
    message(paste("No species found matching:", species_name))
    message("Try using a more general search term or check the spelling")
    return(tibble::tibble())
  }
  
  # Limit results
  if (nrow(species_df) > show_top_n) {
    message(paste("Showing top", show_top_n, "results out of", nrow(species_df), "total matches"))
    species_df <- species_df |> dplyr::slice_head(n = show_top_n)
  }
  
  # Display helpful information
  message("\nSpecies search results:")
  message("Use the 'species_id' column value as the 'species' parameter in STRING DB functions")
  
  return(species_df)
}