# ----------------------------------------------------------------------------
# importProteomeDiscovererTMTData
# ----------------------------------------------------------------------------
#' @title Import Proteome Discoverer TMT Data
#'
#' @description Imports and processes protein-level TMT data from Proteome Discoverer.
#' This function handles reading of .xlsx files, reshaping the data from wide to long format,
#' and returns a standardized list for use in the MultiScholaR workflow.
#'
#' @param filepath Path to the Proteome Discoverer exported .xlsx, .csv, or .tsv file.
#'
#' @return A list containing three elements:
#'   \item{data}{A tibble with the processed and reshaped data.}
#'   \item{data_type}{A character string, hardcoded to "protein".}
#'   \item{column_mapping}{A list mapping standard column names to the names in the processed data.}
#'
#' @importFrom readxl read_excel
#' @importFrom vroom vroom
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr rename select starts_with
#' @importFrom logger log_info
#' @importFrom purrr map_dfr
#' @importFrom tools file_ext
#' @importFrom dplyr rename_with select mutate everything
#' @importFrom stringr str_extract
#' @export
importProteomeDiscovererTMTData <- function(filepath) {
  message("========================================")
  message("=== Entering importProteomeDiscovererTMTData ===")
  message("========================================")
  message(sprintf("   Arg: filepath = %s", filepath))
  message(sprintf("   Arg: filepath class = %s", paste(class(filepath), collapse = ", ")))
  message(sprintf("   Arg: filepath exists = %s", file.exists(filepath)))
  log_info(paste("Starting Proteome Discoverer TMT import from:", filepath))
  
  process_single_file <- function(file_path, batch_name = NULL) {
    message(sprintf("   --- Entering process_single_file ---"))
    message(sprintf("      Arg: file_path = %s", file_path))
    message(sprintf("      Arg: batch_name = %s", ifelse(is.null(batch_name), "NULL", batch_name)))
    log_info(sprintf("Reading file: %s", basename(file_path)))
    
    # Read file based on extension
    message(sprintf("      Step: Checking file extension..."))
    if (endsWith(file_path, ".xlsx")) {
      message(sprintf("      Step: File is .xlsx, using readxl"))
      if (!requireNamespace("readxl", quietly = TRUE)) {
        stop("The 'readxl' package is required for .xlsx files.", call. = FALSE)
      }
      data <- tryCatch({
        message(sprintf("      Step: Reading Excel file..."))
        readxl::read_excel(file_path)
      }, error = function(e) {
        message(sprintf("      ERROR reading Excel file: %s", e$message))
        stop("Failed to read Excel file: ", e$message, call. = FALSE)
      })
    } else {
      message(sprintf("      Step: File is not .xlsx, using vroom"))
      data <- tryCatch({
        message(sprintf("      Step: Reading file with vroom..."))
        vroom::vroom(file_path, show_col_types = FALSE)
      }, error = function(e) {
        message(sprintf("      ERROR reading file: %s", e$message))
        stop("Failed to read file: ", e$message, call. = FALSE)
      })
    }
    
    message(sprintf("      Step: File read complete"))
    message(sprintf("      Data State: Dims = %d rows, %d cols", nrow(data), ncol(data)))
    message(sprintf("      Data State: Column names (first 10): %s", paste(head(names(data), 10), collapse = ", ")))
    log_info(sprintf("File contains %d rows and %d columns", nrow(data), ncol(data)))
    
    # Check for required Accession column
    if ("Accession" %in% names(data)) {
      data <- data |> dplyr::rename(Protein.Ids = "Accession")
    } else {
      stop("Required column 'Accession' not found in file: ", basename(file_path), 
           ". Available columns: ", paste(head(names(data), 10), collapse = ", "), call. = FALSE)
    }
    
    # Count abundance columns before renaming
    abundance_cols_before <- sum(grepl("^Abundance: ", names(data)))
    if (abundance_cols_before == 0) {
      stop("No 'Abundance: ' columns found in file: ", basename(file_path), 
           ". This does not appear to be a valid TMT export file.", call. = FALSE)
    }
    log_info(sprintf("Found %d 'Abundance: ' columns before renaming", abundance_cols_before))
    
    # Rename columns to be unique BEFORE pivoting
    if (!is.null(batch_name)) {
      data <- data |>
        dplyr::rename_with(
          ~ paste(
              batch_name, 
              gsub("Abundance: F[0-9]+: ([0-9]+[A-Z]?), (.+)", "\\1_\\2", .x), 
              sep = "_"
            ),
          .cols = dplyr::starts_with("Abundance: ")
        )
    } else {
        # Fallback for single file without a batch name
        data <- data |>
        dplyr::rename_with(
          ~ gsub("Abundance: F[0-9]+: ([0-9]+[A-Z]?), (.+)", "\\1_\\2", .x),
          .cols = dplyr::starts_with("Abundance: ")
        )
    }

    # Validate abundance columns exist after renaming
    # FIXED: Store the exact column names that were renamed, not a pattern match
    if (!is.null(batch_name)) {
      # Columns should start with batch_name followed by the pattern
      abundance_cols <- names(data)[grepl(paste0("^", batch_name, "_[0-9]+"), names(data))]
    } else {
      # For single files, columns should match the direct gsub pattern
      abundance_cols <- names(data)[grepl("^[0-9]+[A-Z]?_", names(data))]
    }
    
    if (length(abundance_cols) == 0) {
      message(sprintf("      ERROR: No abundance columns found after renaming"))
      message(sprintf("      Batch name: %s", ifelse(is.null(batch_name), "NULL", batch_name)))
      message(sprintf("      All column names: %s", paste(names(data), collapse = ", ")))
      stop("No abundance columns found after renaming in file: ", basename(file_path), 
           ". The column naming pattern may not match expected TMT format.", call. = FALSE)
    }
    message(sprintf("      Found %d abundance columns to pivot: %s", 
                    length(abundance_cols), paste(head(abundance_cols, 5), collapse = ", ")))
    log_info(sprintf("Found %d abundance columns to pivot", length(abundance_cols)))

    # Pivot the uniquely named columns
    message(sprintf("      Step: Pivoting %d columns to long format...", length(abundance_cols)))
    log_info("Pivoting data to long format...")
    long_data <- tryCatch({
      data |>
        tidyr::pivot_longer(
          # FIXED: Use explicit column selection instead of pattern matching
          cols = dplyr::all_of(abundance_cols),
          names_to = "Run",
          values_to = "Abundance"
        )
    }, error = function(e) {
      message(sprintf("      ERROR pivoting data: %s", e$message))
      stop("Failed to pivot data: ", e$message, call. = FALSE)
    })
    
    message(sprintf("      Step: Pivot complete. Pivoted to %d rows", nrow(long_data)))
    log_info(sprintf("Pivoted to %d rows", nrow(long_data)))
    
    if (!is.null(batch_name)) {
      long_data$Batch <- batch_name
      message(sprintf("      Step: Added Batch column with value '%s'", batch_name))
    }
    
    message(sprintf("   --- Exiting process_single_file. Returning %d rows ---", nrow(long_data)))
    return(long_data)
  }
  
  # Check if the file is a zip archive
  message(sprintf("   Step: Checking if file is ZIP..."))
  message(sprintf("   File extension: %s", tools::file_ext(filepath)))
  if (tolower(tools::file_ext(filepath)) == "zip") {
    message(sprintf("   Condition TRUE: File is ZIP archive"))
    log_info("ZIP archive detected. Unzipping and processing batch files.")
    
    message(sprintf("   Step: Creating temp directory..."))
    temp_dir <- tempfile()
    dir.create(temp_dir)
    message(sprintf("   Temp directory: %s", temp_dir))
    
    message(sprintf("   Step: Unzipping file..."))
    unzip(filepath, exdir = temp_dir)
    message(sprintf("   Step: Unzip complete"))
    
    message(sprintf("   Step: Searching for data files in extracted directory..."))
    files_to_process <- list.files(temp_dir, pattern = "\\.(xlsx|csv|tsv)$", 
                                     full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
    message(sprintf("   Step: File search complete"))
    message(sprintf("   Found %d files", length(files_to_process)))
    
    if (length(files_to_process) == 0) {
      message(sprintf("   ERROR: No data files found in ZIP"))
      stop("No data files (.xlsx, .csv, .tsv) found in the ZIP archive.")
    }
    
    message(sprintf("   Files found: %s", paste(basename(files_to_process), collapse = ", ")))
    log_info(sprintf("Found %d files in ZIP: %s", length(files_to_process), 
                     paste(basename(files_to_process), collapse = ", ")))
    
    # Use purrr::imap to iterate and process files with better error handling
    message(sprintf("   Step: Starting iteration over %d files...", length(files_to_process)))
    all_data_list <- purrr::imap(files_to_process, ~{
      # Use the index provided by imap to create a guaranteed unique batch name
      batch_name <- paste0("b", .y) # .y is the index (1, 2, 3...)
      
      message(sprintf("   >>>>>> PROCESSING FILE %d of %d <<<<<<", .y, length(files_to_process)))
      message(sprintf("   File name: %s", basename(.x)))
      message(sprintf("   Batch name: %s", batch_name))
      log_info(paste("Processing file", .y, "of", length(files_to_process), ":", 
                     basename(.x), "as Batch:", batch_name))
      
      result <- tryCatch({
        process_single_file(.x, batch_name = batch_name)
      }, error = function(e) {
        message(sprintf("   !!!!! ERROR IN FILE %d (%s) !!!!!", .y, basename(.x)))
        message(sprintf("   Error message: %s", e$message))
        stop("Error processing file ", .y, " (", basename(.x), "): ", e$message, call. = FALSE)
      })
      
      message(sprintf("   >>>>>> FILE %d COMPLETE <<<<<<", .y))
      message(sprintf("   Returned %d rows", nrow(result)))
      return(result)
    })
    
    # Combine all processed data with error handling
    message(sprintf("   Step: Combining data from all %d files...", length(all_data_list)))
    log_info("Combining data from all files...")
    all_data <- tryCatch({
      dplyr::bind_rows(all_data_list)
    }, error = function(e) {
      message(sprintf("   ERROR combining data: %s", e$message))
      stop("Error combining data from multiple files. Files may have incompatible structures: ", 
           e$message, call. = FALSE)
    })
    message(sprintf("   Step: Data combination complete. Combined %d rows", nrow(all_data)))
    
    message(sprintf("   Step: Cleaning up temp directory..."))
    unlink(temp_dir, recursive = TRUE)
    message(sprintf("   Step: Cleanup complete"))
    
    final_data <- all_data
    
  } else {
    message(sprintf("   Condition FALSE: File is NOT a ZIP archive"))
    # Process a single non-zip file
    final_data <- process_single_file(filepath)
  }
  
  message(sprintf("   Step: Preparing return value..."))
  message(sprintf("   Final data rows: %d", nrow(final_data)))
  log_info(sprintf("Total reshaped data rows: %d.", nrow(final_data)))
  
  result <- list(
    data = final_data,
    data_type = "protein",
    column_mapping = list(
      protein_col = "Protein.Ids",
      run_col = "Run",
      quantity_col = "Abundance",
      batch_col = if("Batch" %in% names(final_data)) "Batch" else NULL
    )
  )
  
  message("========================================")
  message("=== Exiting importProteomeDiscovererTMTData ===")
  message("========================================")
  
  return(result)
}

