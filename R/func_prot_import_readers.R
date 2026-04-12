# ----------------------------------------------------------------------------
# importDIANNData
# ----------------------------------------------------------------------------
#' @title Import DIA-NN Data
#' @description Imports a DIA-NN report file, optionally using precursor normalization.
#' @param filepath Path to the DIA-NN report file (tsv or parquet).
#' @param use_precursor_norm Logical, whether to use Precursor.Normalised instead of Precursor.Quantity.
#' @return A data frame containing the imported DIA-NN data.
#' @export
importDIANNData <- function(filepath, use_precursor_norm = TRUE) {
  log_info(paste("Starting DIA-NN import from:", filepath))
  
  # Check file exists
  if (!file.exists(filepath)) {
    stop("File not found: ", filepath)
  }
  
  # Read data (Parquet or TSV)
  data <- tryCatch({
    if (grepl("\\.parquet$", filepath, ignore.case = TRUE)) {
      arrow::read_parquet(filepath)
    } else {
      vroom::vroom(filepath, show_col_types = FALSE)
    }
  }, error = function(e) {
    stop("Failed to read file: ", e$message)
  })
  
  log_info(sprintf("Read %d rows and %d columns", nrow(data), ncol(data)))
  
  # Check for required columns
  required_cols <- c("Protein.Group", "Stripped.Sequence", "Run")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    log_error(paste("Missing required columns:", paste(missing_cols, collapse = ', ')))
    log_error(paste("Available columns:", paste(head(names(data), 20), collapse = ', ')))
    stop("Missing required DIA-NN columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Convert Precursor.Normalised if needed
  if (use_precursor_norm && "Precursor.Normalised" %in% names(data)) {
    log_info("Converting Precursor.Normalised to numeric")
    data$Precursor.Normalised <- as.numeric(data$Precursor.Normalised)
  }
  
  # Determine quantity column
  quantity_col <- if (use_precursor_norm && "Precursor.Normalised" %in% names(data)) {
    "Precursor.Normalised"
  } else if ("Precursor.Quantity" %in% names(data)) {
    "Precursor.Quantity"
  } else {
    stop("No suitable quantity column found (Precursor.Normalised or Precursor.Quantity)")
  }
  
  log_info(paste("Using quantity column:", quantity_col))
  
  return(list(
    data = data,
    data_type = "peptide",
    column_mapping = list(
      protein_col = "Protein.Group",
      peptide_col = "Stripped.Sequence",
      run_col = "Run",
      quantity_col = quantity_col,
      qvalue_col = "Q.Value",
      pg_qvalue_col = "PG.Q.Value"
    )
  ))
}

# ----------------------------------------------------------------------------
# importSpectronautData
# ----------------------------------------------------------------------------
#' Import Spectronaut Data
#' 
#' @param filepath Path to Spectronaut report file
#' @param quantity_type Either "pg" for protein or "pep" for peptide quantities
#' @return List with data, data_type, and column_mapping
#' @export
importSpectronautData <- function(filepath, quantity_type = "pg") {
  data <- vroom::vroom(filepath, show_col_types = FALSE)
  
  data_type <- if (quantity_type == "pg") "protein" else "peptide"
  
  return(list(
    data = data,
    data_type = data_type,
    column_mapping = list(
      protein_col = "PG.ProteinGroups",
      peptide_col = if (quantity_type == "pep") "EG.ModifiedSequence" else NULL,
      run_col = "R.FileName",
      quantity_col = if (quantity_type == "pg") "PG.Quantity" else "PEP.Quantity",
      qvalue_col = if (quantity_type == "pg") "PG.Qvalue" else "EG.Qvalue"
    )
  ))
}

# ----------------------------------------------------------------------------
# importFragPipeData
# ----------------------------------------------------------------------------
#' Import FragPipe Data
#' 
#' @param filepath Path to FragPipe output file
#' @param use_maxlfq Whether to use MaxLFQ intensities
#' @return List with data, data_type, and column_mapping
#' @export
importFragPipeData <- function(filepath, use_maxlfq = TRUE) {
  log_info(paste("Starting FragPipe LFQ import from:", filepath))
  
  # Check file exists
  if (!file.exists(filepath)) {
    stop("File not found: ", filepath)
  }
  
  # Read data
  data <- tryCatch({
    vroom::vroom(filepath, show_col_types = FALSE)
  }, error = function(e) {
    stop("Failed to read file: ", e$message)
  })
  
  log_info(sprintf("Read %d rows and %d columns", nrow(data), ncol(data)))
  
  # Find Protein ID column (case-insensitive, handle variations)
  protein_id_candidates <- c("Protein ID", "Protein.ID", "Protein_ID", "Protein", "protein id", "protein.id", "protein_id")
  protein_id_col <- NULL
  
  for (candidate in protein_id_candidates) {
    if (candidate %in% names(data)) {
      protein_id_col <- candidate
      break
    }
  }
  
  # Case-insensitive search if exact match not found
  if (is.null(protein_id_col)) {
    names_lower <- tolower(names(data))
    for (candidate in tolower(protein_id_candidates)) {
      idx <- which(names_lower == candidate)
      if (length(idx) > 0) {
        protein_id_col <- names(data)[idx[1]]
        break
      }
    }
  }
  
  if (is.null(protein_id_col)) {
    log_error(paste("Protein ID column not found. Available columns:", paste(head(names(data), 10), collapse = ", ")))
    stop("Required 'Protein ID' column not found in FragPipe file. Available columns: ", paste(head(names(data), 10), collapse = ", "))
  }
  
  log_info(sprintf("Found Protein ID column: %s", protein_id_col))
  
  # Find all columns ending with "Intensity" (case-insensitive)
  all_intensity_cols <- grep("Intensity$", names(data), value = TRUE, ignore.case = TRUE)
  
  if (length(all_intensity_cols) == 0) {
    log_error("No columns ending with 'Intensity' found in FragPipe file.")
    stop("No intensity columns found. Expected columns ending with 'Intensity' (e.g., 'WLP530_1 Intensity', 'WLP530_1 MaxLFQ Intensity')")
  }
  
  log_info(sprintf("Found %d columns ending with 'Intensity'", length(all_intensity_cols)))
  
  # Separate MaxLFQ and regular Intensity columns
  maxlfq_cols <- grep("MaxLFQ.*Intensity$", all_intensity_cols, value = TRUE, ignore.case = TRUE)
  regular_intensity_cols <- setdiff(all_intensity_cols, maxlfq_cols)
  
  # Select which intensity columns to use
  if (use_maxlfq && length(maxlfq_cols) > 0) {
    intensity_cols <- maxlfq_cols
    log_info(sprintf("Using MaxLFQ Intensity columns (%d columns)", length(intensity_cols)))
  } else {
    intensity_cols <- regular_intensity_cols
    log_info(sprintf("Using regular Intensity columns (%d columns)", length(intensity_cols)))
  }
  
  if (length(intensity_cols) == 0) {
    stop("No suitable intensity columns found. Check use_maxlfq parameter and file format.")
  }
  
  # Extract sample names by removing " Intensity" or " MaxLFQ Intensity" suffix
  sample_names <- gsub("\\s+(MaxLFQ\\s+)?Intensity$", "", intensity_cols, ignore.case = TRUE)
  
  # Select only the protein ID column and intensity columns for pivoting
  cols_to_keep <- c(protein_id_col, intensity_cols)
  data_subset <- data |> dplyr::select(dplyr::all_of(cols_to_keep))
  
  # Convert to long format
  log_info("Converting data from wide to long format...")
  long_data <- tryCatch({
    data_subset |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(intensity_cols),
        names_to = "Run",
        values_to = "Intensity"
      ) |>
      # Clean up Run column names (remove " Intensity" or " MaxLFQ Intensity" suffix)
      dplyr::mutate(
        Run = gsub("\\s+(MaxLFQ\\s+)?Intensity$", "", Run, ignore.case = TRUE)
      ) |>
      # Ensure Intensity is numeric
      dplyr::mutate(Intensity = as.numeric(Intensity)) |>
      # Rename protein ID column to standardized name
      dplyr::rename(Protein.Ids = !!rlang::sym(protein_id_col))
  }, error = function(e) {
    log_error(paste("Error converting to long format:", e$message))
    stop("Failed to convert data to long format: ", e$message)
  })
  
  log_info(sprintf("Converted to long format: %d rows", nrow(long_data)))
  log_info(sprintf("Unique proteins: %d, Unique samples: %d", 
                   length(unique(long_data$Protein.Ids)),
                   length(unique(long_data$Run))))
  
  return(list(
    data = long_data,
    data_type = "protein",
    column_mapping = list(
      protein_col = "Protein.Ids",
      run_col = "Run",
      quantity_col = "Intensity",
      qvalue_col = NULL  # FragPipe doesn't typically include q-values
    )
  ))
}

# ----------------------------------------------------------------------------
# importMaxQuantData
# ----------------------------------------------------------------------------
#' Import MaxQuant Data
#' 
#' @param filepath Path to MaxQuant proteinGroups.txt file
#' @param use_lfq Whether to use LFQ intensities
#' @param filter_contaminants Whether to filter contaminants
#' @return List with data, data_type, and column_mapping
#' @export
importMaxQuantData <- function(filepath, use_lfq = TRUE, filter_contaminants = TRUE) {
  data <- vroom::vroom(filepath, show_col_types = FALSE)
  
  # Filter contaminants and reverse hits if requested
  if (filter_contaminants) {
    if ("Potential.contaminant" %in% names(data)) {
      data <- data[data$Potential.contaminant != "+", ]
    }
    if ("Reverse" %in% names(data)) {
      data <- data[data$Reverse != "+", ]
    }
  }
  
  # Find intensity columns
  intensity_cols <- grep("^Intensity\\.", names(data), value = TRUE)
  lfq_cols <- grep("^LFQ\\.intensity\\.", names(data), value = TRUE)
  
  quantity_cols <- if (use_lfq && length(lfq_cols) > 0) {
    lfq_cols
  } else {
    intensity_cols
  }
  
  return(list(
    data = data,
    data_type = "protein",  # MaxQuant proteinGroups is protein-level
    column_mapping = list(
      protein_col = "Majority.protein.IDs",
      peptide_col = NULL,
      run_col = NULL,  # MaxQuant has wide format
      quantity_cols = quantity_cols,  # Multiple columns for wide format
      qvalue_col = "Q.value"  # If Perseus was used
    )
  ))
}

