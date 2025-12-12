# ============================================================================
# func_prot_import.R
# ============================================================================
# Purpose: Proteomics data import functions
# 
# This file contains functions for importing proteomics data from various
# platforms and formats (DIA-NN, TMT-PD, LFQ-Fragpipe, MaxQuant, Spectronaut).
# Functions in this file are used by mod_prot_import.R and related modules.
#
# Functions to extract here:
# - importDIANNData(): Import DIA-NN report files
# - importFragPipeData(): Import FragPipe output files
# - importProteomeDiscovererTMTData(): Import TMT data from Proteome Discoverer
# - importMaxQuantData(): Import MaxQuant proteinGroups.txt files
# - importSpectronautData(): Import Spectronaut output files
# - formatDIANN(): Convert data to DIA-NN format
# - formatDIANNParquet(): Convert limpa EList to DIA-NN format
# - detectProteomicsFormat(): Auto-detect proteomics data format
#
# Dependencies:
# - dplyr, tidyr, vroom, readxl
# - func_general_filemgmt.R (for file utilities)
# ============================================================================

# TODO: Extract the following functions from their current locations:

# Function 1: importDIANNData()
# Current location: R/mod_prot_import.R
# Lines: ~997-1059
# Description: Imports DIA-NN report files (.tsv or .tsv.gz)
# importDIANNData <- function(filepath, use_precursor_norm = TRUE) {
#   # Extract from R/mod_prot_import.R
# }

# Function 2: importFragPipeData()
# Current location: R/mod_prot_import.R
# Lines: ~1084-1210
# Description: Imports FragPipe output files (protein.tsv)
# importFragPipeData <- function(filepath, use_maxlfq = TRUE) {
#   # Extract from R/mod_prot_import.R
# }

# Function 3: importProteomeDiscovererTMTData()
# Current location: R/file_management.R
# Lines: ~599-829
# Description: Imports TMT data from Proteome Discoverer (.xlsx, .csv, .tsv, or .zip)
# importProteomeDiscovererTMTData <- function(filepath) {
#   # Extract from R/file_management.R
# }

# Function 4: importMaxQuantData()
# Current location: R/mod_prot_import.R
# Lines: ~1211-1252
# Description: Imports MaxQuant proteinGroups.txt files
# importMaxQuantData <- function(filepath, use_lfq = TRUE, filter_contaminants = TRUE) {
#   # Extract from R/mod_prot_import.R
# }

# Function 5: importSpectronautData()
# Current location: R/mod_prot_import.R
# Lines: ~1060-1083
# Description: Imports Spectronaut output files
# importSpectronautData <- function(filepath, quantity_type = "pg") {
#   # Extract from R/mod_prot_import.R
# }

# Function 6: formatDIANN()
# Current location: R/file_management.R
# Lines: ~857-956
# Description: Converts limpa EList object to standard DIA-NN format
# formatDIANN <- function(data_tbl_parquet_filt) {
#   # Extract from R/file_management.R
# }

# Function 7: formatDIANNParquet()
# Current location: R/file_management.R
# Lines: ~436-574
# Description: Converts limpa EList object to standard DIA-NN format (parquet version)
# formatDIANNParquet <- function(data_tbl_parquet_filt) {
#   # Extract from R/file_management.R
# }

# Function 8: detectProteomicsFormat()
# Current location: R/mod_prot_import.R
# Lines: ~884-996
# Description: Auto-detects proteomics data format from file headers
# detectProteomicsFormat <- function(headers, filename, preview_lines = NULL) {
#   # Extract from R/mod_prot_import.R
# }

# Function 9: getDefaultProteomicsConfig()
# Current location: R/mod_prot_import.R
# Lines: ~1253+
# Description: Returns default configuration for proteomics workflows
# getDefaultProteomicsConfig <- function() {
#   # Extract from R/mod_prot_import.R
# }


# ----------------------------------------------------------------------------
# detectProteomicsFormat
# ----------------------------------------------------------------------------
#' Detect Proteomics Data Format
#' 
#' Detects the format of proteomics search results based on headers and filename
#' 
#' @param headers Character vector of column headers
#' @param filename Name of the file
#' @param preview_lines First few lines of the file
#' 
#' @return List with format and confidence score
#' @export
detectProteomicsFormat <- function(headers, filename, preview_lines = NULL) {
  # Convert to lowercase for comparison
  headers_lower <- tolower(headers)
  filename_lower <- tolower(filename)
  
  # DIA-NN detection
  diann_score <- 0
  diann_markers <- c("protein.group", "protein.ids", "protein.names", 
                     "precursor.id", "modified.sequence", "stripped.sequence",
                     "precursor.charge", "q.value", "pg.q.value", "run")
  diann_found <- sum(diann_markers %in% headers_lower)
  diann_score <- diann_found / length(diann_markers)
  
  # Spectronaut detection
  spectronaut_score <- 0
  spectronaut_markers <- c("pg.proteingroups", "pg.proteinaccessions", "pg.genes",
                          "eg.precursorid", "eg.modifiedsequence", "fg.charge",
                          "eg.qvalue", "pg.qvalue", "r.filename", "pg.quantity")
  spectronaut_found <- sum(spectronaut_markers %in% headers_lower)
  spectronaut_score <- spectronaut_found / length(spectronaut_markers)
  
  # FragPipe detection
  fragpipe_score <- 0
  
  # Check for key markers (flexible matching)
  has_protein_id <- any(grepl("^protein\\s+id$|^protein\\.id$|^protein_id$", headers_lower))
  has_protein <- "protein" %in% headers_lower
  has_gene <- "gene" %in% headers_lower
  has_description <- any(grepl("description", headers_lower))
  has_spectral_count <- any(grepl("spectral.*count", headers_lower))
  
  # Strong indicator: columns ending with "Intensity" (case-insensitive)
  intensity_cols <- sum(grepl("intensity$", headers_lower))
  has_intensity_cols <- intensity_cols > 0
  
  # Check for MaxLFQ Intensity columns specifically
  has_maxlfq <- any(grepl("maxlfq.*intensity$", headers_lower))
  
  # Count basic markers found
  basic_markers_found <- sum(c(has_protein_id, has_protein, has_gene, has_description, has_spectral_count))
  
  # Weighted scoring similar to TMT detection
  # Basic markers worth 40% (8% each)
  basic_score <- (basic_markers_found / 5) * 0.4
  
  # Intensity columns presence is worth 40% (strong indicator)
  intensity_score <- if (has_intensity_cols) {
    # Bonus if multiple intensity columns (typical of FragPipe)
    min(0.4, 0.3 + (min(intensity_cols, 10) / 10) * 0.1)
  } else {
    0
  }
  
  # MaxLFQ presence is worth 20% (very specific to FragPipe)
  maxlfq_score <- if (has_maxlfq) 0.2 else 0
  
  fragpipe_score <- basic_score + intensity_score + maxlfq_score
  
  # Filename bonus (cap at 1.0)
  if (grepl("fragpipe|msfragger", filename_lower)) {
    fragpipe_score <- min(1.0, fragpipe_score + 0.1)
  }
  
  # MaxQuant detection
  maxquant_score <- 0
  maxquant_markers <- c("proteins", "majority.protein.ids", "protein.names",
                       "gene.names", "peptide.counts", "unique.peptides",
                       "intensity", "lfq.intensity", "ms.ms.count")
  maxquant_found <- sum(maxquant_markers %in% headers_lower)
  maxquant_score <- maxquant_found / length(maxquant_markers)
  if (grepl("proteingroups", filename_lower)) maxquant_score <- maxquant_score + 0.3
  
  # PD-TMT detection
  pd_tmt_score <- 0
  pd_tmt_markers <- c("protein fdr confidence", "master", "accession", "exp. q-value", "sum pep score")
  pd_tmt_found <- sum(pd_tmt_markers %in% headers_lower)
  abundance_found <- any(grepl("^abundance:", headers_lower))
  
  # Weighted scoring to prevent exceeding 100%
  # Header markers are worth 60% total (12% each)
  header_score <- (pd_tmt_found / length(pd_tmt_markers)) * 0.6
  # Abundance column presence is worth 40%
  abundance_score <- if (abundance_found) 0.4 else 0
  
  pd_tmt_score <- header_score + abundance_score
  
  # Determine best match
  scores <- c(diann = diann_score, 
              spectronaut = spectronaut_score,
              fragpipe = fragpipe_score,
              maxquant = maxquant_score,
              pd_tmt = pd_tmt_score)
  
  best_format <- names(which.max(scores))
  best_score <- max(scores)
  
  if (best_score < 0.3) {
    best_format <- "unknown"
  }
  
  return(list(
    format = best_format,
    confidence = best_score,
    all_scores = scores
  ))
}


# ----------------------------------------------------------------------------
# importDIANNData
# ----------------------------------------------------------------------------
#' Import DIA-NN Data
#' 
#' @param filepath Path to DIA-NN report file
#' @param use_precursor_norm Whether to use Precursor.Normalised values
#' @return List with data, data_type, and column_mapping
#' @export
importDIANNData <- function(filepath, use_precursor_norm = TRUE) {
  log_info(paste("Starting DIA-NN import from:", filepath))
  
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


# ----------------------------------------------------------------------------
# getDefaultProteomicsConfig
# ----------------------------------------------------------------------------
#' Get Default Proteomics Configuration
#' 
#' Returns default configuration settings for proteomics analysis
#' 
#' @return List with default configuration parameters
#' @export
getDefaultProteomicsConfig <- function() {
  list(
    generalParameters = list(
      min_peptides_per_protein = 2,
      min_peptides_per_sample = 2,
      q_value_threshold = 0.01,
      intensity_threshold = 0
    ),
    deAnalysisParameters = list(
      formula_string = "~ 0 + group",
      q_value_threshold = 0.05,
      log2_fc_threshold = 1
    ),
    normalizationParameters = list(
      normalisation_method = "cyclicloess"
    ),
    ruvParameters = list(
      percentage_as_neg_ctrl = 33,
      ruv_k = NULL
    )
  )
}


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


# ----------------------------------------------------------------------------
# formatDIANN
# ----------------------------------------------------------------------------
formatDIANN <- function(data_tbl_parquet_filt) {
  
  # Extract intensity matrix and gene annotations
  intensity_matrix <- data_tbl_parquet_filt$E
  gene_info <- data_tbl_parquet_filt$genes
  
  # Convert intensity matrix to data frame with precursor IDs
  intensity_df <- as.data.frame(intensity_matrix) |>
    rownames_to_column(var = "Precursor.Id")
  
  # Convert to long format
  intensity_long <- intensity_df |>
    pivot_longer(
      cols = -Precursor.Id,
      names_to = "Run",
      values_to = "Log2Intensity"
    ) |>
    filter(!is.na(Log2Intensity))  # Remove missing values
  
  # Join with gene information
  data_tbl_converted <- intensity_long |>
    left_join(gene_info, by = "Precursor.Id") |>
    mutate(
      # Basic required columns
      File.Name = paste0(Run, ".raw"),
      Protein.Ids = Protein.Group,  # Use Protein.Group as Protein.Ids
      
      # Extract sequence and charge from Precursor.Id
      Stripped.Sequence = gsub("\\d+$", "", Precursor.Id),  # Remove charge number
      Modified.Sequence = Stripped.Sequence,  # Assume no modifications shown
      Precursor.Charge = as.numeric(gsub(".*?(\\d+)$", "\\1", Precursor.Id)),  # Extract charge
      
      # Convert log2 intensities back to linear scale
      Precursor.Quantity = 2^Log2Intensity,
      Precursor.Normalised = 2^Log2Intensity,
      
      # Set reasonable defaults for quality metrics (since data was pre-filtered)
      Q.Value = 0.001,
      PEP = 0.001,
      Global.Q.Value = 0.001,
      Protein.Q.Value = 0.01,
      PG.Q.Value = 0.001,
      Global.PG.Q.Value = 0.001,
      GG.Q.Value = 0.001,
      Translated.Q.Value = 0,
      Lib.Q.Value = 0.001,
      Lib.PG.Q.Value = 0.001,
      
      # Protein-level quantities (same as precursor for now)
      PG.Quantity = Precursor.Quantity,
      PG.Normalised = Precursor.Normalised,
      PG.MaxLFQ = Precursor.Normalised,
      Genes.Quantity = Precursor.Quantity,
      Genes.Normalised = Precursor.Normalised,
      Genes.MaxLFQ = Precursor.Normalised,
      Genes.MaxLFQ.Unique = Precursor.Normalised,
      
      # Quality and technical columns (set defaults)
      Quantity.Quality = 1.0,
      RT = NA_real_,
      RT.Start = NA_real_,
      RT.Stop = NA_real_,
      iRT = NA_real_,
      Predicted.RT = NA_real_,
      Predicted.iRT = NA_real_,
      First.Protein.Description = "",
      Ms1.Profile.Corr = NA_real_,
      Ms1.Area = NA_real_,
      Ms1.Normalised = NA_real_,
      Normalisation.Factor = 1.0,
      Evidence = 1.0,
      Spectrum.Similarity = NA_real_,
      Averagine = NA_real_,
      Mass.Evidence = NA_real_,
      CScore = 1.0,
      Fragment.Quant.Raw = "",
      Fragment.Correlations = "",
      MS2.Scan = NA_real_,
      IM = 0,
      iIM = 0,
      Predicted.IM = 0,
      Predicted.iIM = 0
    ) |>
    # Select columns in typical DIA-NN order
    dplyr::select(File.Name, Run, Protein.Group, Protein.Ids, Protein.Names, Genes,
           PG.Quantity, PG.Normalised, PG.MaxLFQ,
           Genes.Quantity, Genes.Normalised, Genes.MaxLFQ, Genes.MaxLFQ.Unique,
           Modified.Sequence, Stripped.Sequence, Precursor.Id, Precursor.Charge,
           Q.Value, PEP, Global.Q.Value, Protein.Q.Value, PG.Q.Value,
           Global.PG.Q.Value, GG.Q.Value, Translated.Q.Value, Proteotypic,
           Precursor.Quantity, Precursor.Normalised, Quantity.Quality,
           RT, RT.Start, RT.Stop, iRT, Predicted.RT, Predicted.iRT,
           First.Protein.Description, Lib.Q.Value, Lib.PG.Q.Value,
           Ms1.Profile.Corr, Ms1.Area, Ms1.Normalised, Normalisation.Factor,
           Evidence, Spectrum.Similarity, Averagine, Mass.Evidence, CScore,
           Fragment.Quant.Raw, Fragment.Correlations, MS2.Scan,
           IM, iIM, Predicted.IM, Predicted.iIM)
  
  return(data_tbl_converted)
}


# ----------------------------------------------------------------------------
# formatDIANNParquet
# ----------------------------------------------------------------------------
#' @title Convert limpa EList object to standard DIA-NN format
#' @description Converts a limpa EList object (typically from processed proteomics data) 
#'   into the standard DIA-NN report format. This function extracts intensity matrices 
#'   and gene annotations, then restructures the data to match DIA-NN output conventions.
#' @param data_tbl_parquet_filt A limpa EList object containing:
#'   \itemize{
#'     \item \code{E}: Intensity matrix with precursor IDs as rownames
#'     \item \code{genes}: Data frame with gene/protein annotations including columns like 
#'           Protein.Group, Protein.Names, Genes, etc.
#'   }
#' @return A data frame in DIA-NN report format containing:
#'   \itemize{
#'     \item File identification: \code{File.Name}, \code{Run}
#'     \item Protein information: \code{Protein.Group}, \code{Protein.Ids}, \code{Protein.Names}, \code{Genes}
#'     \item Quantification data: \code{PG.Quantity}, \code{PG.Normalised}, \code{PG.MaxLFQ}, 
#'           \code{Genes.Quantity}, \code{Genes.Normalised}, \code{Genes.MaxLFQ}, etc.
#'     \item Precursor details: \code{Modified.Sequence}, \code{Stripped.Sequence}, 
#'           \code{Precursor.Id}, \code{Precursor.Charge}, \code{Precursor.Quantity}
#'     \item Quality metrics: \code{Q.Value}, \code{PEP}, \code{Global.Q.Value}, etc.
#'     \item Technical parameters: \code{RT}, \code{IM}, \code{Evidence}, etc.
#'   }
#' @details This function performs several key transformations:
#'   \itemize{
#'     \item Extracts log2 intensity data and converts to linear scale
#'     \item Converts wide-format intensity matrix to long format
#'     \item Joins intensity data with protein/gene annotations
#'     \item Extracts precursor sequence and charge information from precursor IDs
#'     \item Sets reasonable default values for quality metrics (since input data is pre-filtered)
#'     \item Reorders columns to match standard DIA-NN report structure
#'   }
#' @note The function assumes that precursor IDs contain charge information as a suffix 
#'   (e.g., "PEPTIDE123" where "123" is the charge). Quality metrics are set to 
#'   conservative default values since the input data is typically pre-filtered.
#' @importFrom dplyr mutate select left_join
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @export
formatDIANNParquet <- function(data_tbl_parquet_filt) {
  library(dplyr)
  library(tidyr)
  
  # Extract intensity matrix and gene annotations
  intensity_matrix <- data_tbl_parquet_filt$E
  gene_info <- data_tbl_parquet_filt$genes
  
  # Convert intensity matrix to data frame with precursor IDs
  intensity_df <- as.data.frame(intensity_matrix) %>%
    rownames_to_column(var = "Precursor.Id")
  
  # Convert to long format
  intensity_long <- intensity_df %>%
    pivot_longer(
      cols = -Precursor.Id,
      names_to = "Run",
      values_to = "Log2Intensity"
    ) %>%
    filter(!is.na(Log2Intensity))  # Remove missing values
  
  # Join with gene information
  data_tbl_converted <- intensity_long %>%
    left_join(gene_info, by = "Precursor.Id") %>%
    mutate(
      # Basic required columns
      File.Name = paste0(Run, ".raw"),
      Protein.Ids = Protein.Group,  # Use Protein.Group as Protein.Ids
      
      # Extract sequence and charge from Precursor.Id
      Stripped.Sequence = gsub("\\d+$", "", Precursor.Id),  # Remove charge number
      Modified.Sequence = Stripped.Sequence,  # Assume no modifications shown
      Precursor.Charge = as.numeric(gsub(".*?(\\d+)$", "\\1", Precursor.Id)),  # Extract charge
      
      # Convert log2 intensities back to linear scale
      Precursor.Quantity = 2^Log2Intensity,
      Precursor.Normalised = 2^Log2Intensity,
      
      # Set reasonable defaults for quality metrics (since data was pre-filtered)
      Q.Value = 0.001,
      PEP = 0.001,
      Global.Q.Value = 0.001,
      Protein.Q.Value = 0.01,
      PG.Q.Value = 0.001,
      Global.PG.Q.Value = 0.001,
      GG.Q.Value = 0.001,
      Translated.Q.Value = 0,
      Lib.Q.Value = 0.001,
      Lib.PG.Q.Value = 0.001,
      
      # Protein-level quantities (same as precursor for now)
      PG.Quantity = Precursor.Quantity,
      PG.Normalised = Precursor.Normalised,
      PG.MaxLFQ = Precursor.Normalised,
      Genes.Quantity = Precursor.Quantity,
      Genes.Normalised = Precursor.Normalised,
      Genes.MaxLFQ = Precursor.Normalised,
      Genes.MaxLFQ.Unique = Precursor.Normalised,
      
      # Quality and technical columns (set defaults)
      Quantity.Quality = 1.0,
      RT = NA_real_,
      RT.Start = NA_real_,
      RT.Stop = NA_real_,
      iRT = NA_real_,
      Predicted.RT = NA_real_,
      Predicted.iRT = NA_real_,
      First.Protein.Description = "",
      Ms1.Profile.Corr = NA_real_,
      Ms1.Area = NA_real_,
      Ms1.Normalised = NA_real_,
      Normalisation.Factor = 1.0,
      Evidence = 1.0,
      Spectrum.Similarity = NA_real_,
      Averagine = NA_real_,
      Mass.Evidence = NA_real_,
      CScore = 1.0,
      Fragment.Quant.Raw = "",
      Fragment.Correlations = "",
      MS2.Scan = NA_real_,
      IM = 0,
      iIM = 0,
      Predicted.IM = 0,
      Predicted.iIM = 0
    ) %>%
    # Select columns in typical DIA-NN order
    dplyr::select(File.Name, Run, Protein.Group, Protein.Ids, Protein.Names, Genes,
           PG.Quantity, PG.Normalised, PG.MaxLFQ,
           Genes.Quantity, Genes.Normalised, Genes.MaxLFQ, Genes.MaxLFQ.Unique,
           Modified.Sequence, Stripped.Sequence, Precursor.Id, Precursor.Charge,
           Q.Value, PEP, Global.Q.Value, Protein.Q.Value, PG.Q.Value,
           Global.PG.Q.Value, GG.Q.Value, Translated.Q.Value, Proteotypic,
           Precursor.Quantity, Precursor.Normalised, Quantity.Quality,
           RT, RT.Start, RT.Stop, iRT, Predicted.RT, Predicted.iRT,
           First.Protein.Description, Lib.Q.Value, Lib.PG.Q.Value,
           Ms1.Profile.Corr, Ms1.Area, Ms1.Normalised, Normalisation.Factor,
           Evidence, Spectrum.Similarity, Averagine, Mass.Evidence, CScore,
           Fragment.Quant.Raw, Fragment.Correlations, MS2.Scan,
           IM, iIM, Predicted.IM, Predicted.iIM)
  
  return(data_tbl_converted)
}


# ----------------------------------------------------------------------------
# PeptideQuantitativeDataDiann
# ----------------------------------------------------------------------------
#' @export
PeptideQuantitativeDataDiann <- function( peptide_data
                                          , design_matrix
                                          , sample_id = "Run"
                                          , group_id = "group"
                                          , technical_replicate_id = "replicates"
                                          , args = NA) {



  peptide_data <- new( "PeptideQuantitativeData"

    # Protein vs Sample quantitative data
    , peptide_data = peptide_data
    , protein_id_column = "Protein.Ids"
    , peptide_sequence_column = "Stripped.Sequence"
    , q_value_column = "Q.Value"
    , global_q_value_column = "Global.Q.Value"
    , proteotypic_peptide_sequence_column = "Proteotypic"
    , raw_quantity_column = "Precursor.Quantity"
    , norm_quantity_column = "Precursor.Normalised"
    , is_logged_data = FALSE

    # Design Matrix Information
    , design_matrix = design_matrix
    , sample_id= sample_id
    , group_id= group_id
    , technical_replicate_id= technical_replicate_id
    , args = args
  )

  peptide_data

}


# ----------------------------------------------------------------------------
# extractOrganismsFromFasta
# ----------------------------------------------------------------------------
#' @title Extract Organism Information from FASTA File
#' 
#' @description Parses FASTA file headers to extract organism name (OS) and 
#' taxonomy ID (OX) for each protein entry. Used for analyzing mixed-species
#' FASTA databases.
#' 
#' @param fasta_file_path Path to the FASTA file
#' 
#' @return A tibble with columns: uniprot_acc, organism_name, taxon_id
#' 
#' @importFrom seqinr read.fasta
#' @importFrom dplyr tibble mutate
#' @importFrom stringr str_detect str_replace_all str_extract
#' @importFrom logger log_info log_warn
#' @export
extractOrganismsFromFasta <- function(fasta_file_path) {
  log_info(paste("Extracting organism information from FASTA:", fasta_file_path))
  
  # Check file exists

  if (!file.exists(fasta_file_path)) {
    stop("FASTA file not found: ", fasta_file_path)
  }
  
  # Read FASTA file
  aa_seqinr <- tryCatch({
    seqinr::read.fasta(
      file = fasta_file_path
      , seqtype = "AA"
      , whole.header = TRUE
      , as.string = TRUE
    )
  }, error = function(e) {
    stop("Failed to read FASTA file: ", e$message)
  })
  
  headers <- names(aa_seqinr)
  total_entries <- length(headers)
  
  log_info(sprintf("Processing %d FASTA entries for organism information...", total_entries))
  
  # Helper function to extract field value from FASTA header
  # Replicates getFastaFields logic inline for efficiency
  extractField <- function(header_string, field_name) {
    pattern <- paste0(field_name, "=")
    if (!stringr::str_detect(header_string, pattern)) {
      return(NA_character_)
    }
    # Extract the field value (everything after = until next field or end)
    extracted <- stringr::str_replace_all(
      header_string
      , paste0("(.*)", field_name, "=(.*?)(\\s..=.*|$)")
      , "\\2"
    )
    return(trimws(extracted))
  }
  
  # Parse each header
  result_list <- vector("list", length(headers))
  
  for (i in seq_along(headers)) {
    header <- headers[i]
    
    # Extract accession from standard UniProt format: >sp|ACCESSION|ID or >tr|ACCESSION|ID
    parts <- strsplit(header, "\\|", fixed = FALSE)[[1]]
    
    if (length(parts) >= 2) {
      accession <- parts[2]
    } else {
      # Fallback: use first word
      accession <- strsplit(header, "\\s+")[[1]][1]
      accession <- sub("^>", "", accession)
    }
    
    # Extract organism fields
    organism_name <- extractField(header, "OS")
    taxon_id_str <- extractField(header, "OX")
    taxon_id <- suppressWarnings(as.integer(taxon_id_str))
    
    result_list[[i]] <- list(
      uniprot_acc = accession
      , organism_name = organism_name
      , taxon_id = taxon_id
    )
  }
  
  # Bind into tibble
  organism_tbl <- dplyr::bind_rows(result_list)
  
  # Count organisms found
  n_with_organism <- sum(!is.na(organism_tbl$organism_name))
  n_with_taxon <- sum(!is.na(organism_tbl$taxon_id))
  
  log_info(sprintf(
    "Extracted organism info: %d/%d have organism name, %d/%d have taxon ID"
    , n_with_organism, total_entries
    , n_with_taxon, total_entries
  ))
  
  return(organism_tbl)
}


# ----------------------------------------------------------------------------
# analyzeOrganismDistribution
# ----------------------------------------------------------------------------
#' @title Analyze Organism Distribution in Proteomics Data
#' 
#' @description Analyzes the distribution of proteins across different organisms
#' by matching protein IDs from search results against FASTA organism annotations.
#' Used to identify the primary organism in mixed-species databases.
#' 
#' @param protein_ids Character vector of protein IDs from search results
#' @param organism_mapping Tibble from extractOrganismsFromFasta() with 
#'   columns: uniprot_acc, organism_name, taxon_id
#' 
#' @return A tibble with columns:
#'   \item{organism_name}{Name of the organism}
#'   \item{taxon_id}{NCBI taxonomy ID}
#'   \item{protein_count}{Number of proteins matched to this organism}
#'   \item{percentage}{Percentage of total proteins}
#'   \item{rank}{Rank by protein count (1 = highest)}
#' 
#' @importFrom dplyr tibble filter group_by summarise arrange mutate n desc row_number left_join
#' @importFrom stringr str_split
#' @importFrom logger log_info
#' @export
analyzeOrganismDistribution <- function(protein_ids, organism_mapping) {
  log_info(sprintf("Analyzing organism distribution for %d protein IDs", length(protein_ids)))
  
  # Handle protein groups (semicolon-separated IDs)
  # Extract all unique accessions from protein groups
  all_accessions <- unique(unlist(stringr::str_split(protein_ids, ";")))
  all_accessions <- trimws(all_accessions)
  all_accessions <- all_accessions[all_accessions != ""]
  
  log_info(sprintf("Expanded to %d unique accessions", length(all_accessions)))
  
  # Clean accessions (remove isoform numbers for matching)
  clean_accession <- function(acc) {
    # Remove isoform suffix (e.g., P12345-2 -> P12345)
    sub("-\\d+$", "", acc)
  }
  
  # Create lookup table with cleaned accessions
  organism_mapping_clean <- organism_mapping |>
    dplyr::mutate(cleaned_acc = clean_accession(uniprot_acc))
  
  # Match accessions to organisms
  accession_tbl <- dplyr::tibble(
    accession = all_accessions
    , cleaned_acc = clean_accession(all_accessions)
  )
  
  # Join with organism mapping
  matched <- accession_tbl |>
    dplyr::left_join(
      organism_mapping_clean |>
        dplyr::select(cleaned_acc, organism_name, taxon_id) |>
        dplyr::distinct()
      , by = "cleaned_acc"
    )
  
  n_matched <- sum(!is.na(matched$organism_name))
  n_unmatched <- sum(is.na(matched$organism_name))
  
  log_info(sprintf("Matched %d accessions, %d unmatched", n_matched, n_unmatched))
  
  # Calculate distribution by organism
  distribution <- matched |>
    dplyr::filter(!is.na(organism_name)) |>
    dplyr::group_by(organism_name, taxon_id) |>
    dplyr::summarise(
      protein_count = dplyr::n()
      , .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(protein_count)) |>
    dplyr::mutate(
      percentage = round(protein_count / sum(protein_count) * 100, 2)
      , rank = dplyr::row_number()
    )
  
  # Add unmatched as a separate row if any
  if (n_unmatched > 0) {
    unmatched_row <- dplyr::tibble(
      organism_name = "[Unmatched/Unknown]"
      , taxon_id = NA_integer_
      , protein_count = n_unmatched
      , percentage = round(n_unmatched / length(all_accessions) * 100, 2)
      , rank = nrow(distribution) + 1
    )
    distribution <- dplyr::bind_rows(distribution, unmatched_row)
  }
  
  log_info(sprintf("Found %d distinct organisms in data", nrow(distribution) - ifelse(n_unmatched > 0, 1, 0)))
  
  return(distribution)
}
