# ----------------------------------------------------------------------------
# Quantity columns: missing-value tokens and non-numeric text (P-05)
# ----------------------------------------------------------------------------
#' The text a proteomics quantity column may hold for a missing value (P-05, user-approved 2026-11-11): it is read as
#' NA, and the import gives a notice with the count per text token (an empty cell and "NA" are the file format's own
#' missing value and are not counted). Case-sensitive, after trimming blanks. Any other text that is not a number
#' refuses the import (`checkProtImportQuantityColumns`).
#' @noRd
PROT_IMPORT_MISSING_TOKENS <- c("", "NA", "NaN", "Filtered", "n.q.")

#' Check a proteomics reader's quantity columns as its file read returned them, before the table is reshaped (P-05).
#' A text column: a missing-value token (`PROT_IMPORT_MISSING_TOKENS`) becomes NA and the rest must be numeric-like
#' (`omicsImportValuesAreNumericLike`, the lipid, metab and trx rule: finite numbers), and is converted to doubles.
#' A column vroom typed itself: the cells it could not parse as a number (it reads them as NA and records them) are
#' judged by the same rule, NaN becomes NA, and a column it read as a "number" (digit-grouping marks dropped, so a
#' decimal comma moves the decimal point) is refused. A column of plain numbers is returned as it was read.
#' @param data The table as the file read returned it, with its columns in file order.
#' @param columns The quantity columns' names in `data`; absent names are ignored.
#' @param path The file, read again only for the example values of a refused "number" column; NULL for none.
#' @return `data`. The text tokens read as NA (NaN, Filtered, n.q.) are signalled as a message of class
#'   `prot_import_missing_tokens` (field `counts`: the count per token); empty cells and "NA", the file format's own
#'   missing value, are not counted. Other text is an error of class `prot_import_non_numeric_quantity` naming the
#'   first such column and up to three values.
#' @noRd
checkProtImportQuantityColumns <- function(data, columns, path = NULL) {
    columns <- intersect(columns, names(data))
    from_vroom <- identical(typeof(attr(data, "problems")), "externalptr")
    # vroom's own "parsing issues" warning is replaced by the notice or the refusal below.
    quiet <- function(expr) withCallingHandlers(expr, warning = function(w) {
        if (grepl("parsing issues", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
    })
    # An empty cell or "NA" is the file format's own missing value: read as NA as before, and not counted (a file
    # with empty cells and plain numbers imports exactly as it did, with no notice).
    labels <- setdiff(PROT_IMPORT_MISSING_TOKENS, c("", "NA"))
    counts <- stats::setNames(integer(length(labels)), labels)
    count <- function(tokens) {
        found <- table(tokens[tokens %in% labels])
        counts[names(found)] <<- counts[names(found)] + as.integer(found)
    }
    refused <- list()
    missing <- list()
    unparsed <- list()
    if (from_vroom) {
        # vroom parses lazily: reading a column's missing values parses it and records the cells it could not parse.
        plain <- columns[vapply(columns, function(column) is.double(data[[column]]) && !is.object(data[[column]]), logical(1))]
        for (column in plain) missing[[column]] <- quiet(is.na(data[[column]]))
        problems <- quiet(vroom::problems(data, lazy = TRUE))
        unparsed <- split(trimws(problems$actual), names(data)[problems$col])
    }
    for (column in columns) {
        values <- data[[column]]
        if (is.logical(values) || is.integer(values)) next
        if (is.double(values) && !is.object(values)) {
            if (!from_vroom) next
            if (inherits(attr(data, "spec")$cols[[match(column, names(data))]], "collector_number")) {
                refused[[column]] <- protImportNumberColumnExamples(path, match(column, names(data)), attr(data, "spec")$delim)
                next
            }
            if (!any(missing[[column]])) next
            text <- unparsed[[column]] %||% character()
            other <- unique(text[!text %in% PROT_IMPORT_MISSING_TOKENS])
            if (length(other)) {
                refused[[column]] <- other
                next
            }
            not_a_number <- is.nan(values)
            count(text)
            counts[["NaN"]] <- counts[["NaN"]] + sum(not_a_number)
            if (any(not_a_number)) data[[column]][not_a_number] <- NA_real_
            next
        }
        text <- trimws(as.character(values))
        text[is.na(text)] <- ""
        token <- text %in% PROT_IMPORT_MISSING_TOKENS
        if (!omicsImportValuesAreNumericLike(text[!token])) {
            number <- suppressWarnings(as.numeric(text[!token]))
            refused[[column]] <- unique(text[!token][is.na(number) | !is.finite(number) |
                grepl("^[+-]?0[xX]", text[!token])])
            next
        }
        count(text[token])
        number <- rep(NA_real_, length(text))
        number[!token] <- as.numeric(text[!token])
        data[[column]] <- number
    }
    if (length(refused)) {
        examples <- utils::head(refused[[1L]], 3L)
        rlang::abort(paste0(
            sprintf("Quantity column '%s' contains text that is not a number", names(refused)[[1L]]),
            if (length(examples)) paste0(": ", paste(sprintf("\"%s\"", examples), collapse = ", ")),
            if (length(refused) > 1L) sprintf(" (and %d other quantity column(s))", length(refused) - 1L),
            ". Only empty cells, ", paste(setdiff(PROT_IMPORT_MISSING_TOKENS, ""), collapse = ", "),
            " are read as missing values, and numbers need a decimal point and no thousands separators.",
            " Correct the file and import it again."),
            class = "prot_import_non_numeric_quantity", column = names(refused)[[1L]], examples = examples, call = NULL)
    }
    counts <- counts[counts > 0L]
    if (length(counts)) {
        rlang::inform(sprintf("Missing-value text in the quantity column(s) was read as NA: %s.",
            paste(sprintf("%d x %s", counts, names(counts)), collapse = ", ")),
            class = "prot_import_missing_tokens", counts = counts)
    }
    data
}

#' Up to the first 100 values of a column vroom read as a "number" that are not plain numbers, as the file writes them
#' (the examples of `checkProtImportQuantityColumns`'s refusal).
#' @param path The file, or NULL.
#' @param index The column's position in the file.
#' @param delim The delimiter the first read used (its column specification's).
#' @return The values; none when there is no file, it cannot be read again, or its first rows are plain numbers.
#' @noRd
protImportNumberColumnExamples <- function(path, index, delim = NULL) {
    if (is.null(path)) return(character())
    # The refusal stands without examples: a failed second read only drops them.
    text <- tryCatch(vroom::vroom(path, delim = delim, col_select = dplyr::all_of(index), n_max = 100L,
        show_col_types = FALSE, col_types = vroom::cols(.default = vroom::col_character()))[[1L]],
        error = function(e) character())
    text <- trimws(text[!is.na(text)])
    unique(text[is.na(suppressWarnings(as.numeric(text)))])
}

#' Re-signal an import error under a prefix, as the readers' own `stop(prefix, message)`, keeping the class of a
#' quantity-column refusal (P-05) so a caller can tell it from a failed read.
#' @param prefix The text put before the error's message.
#' @param error The error.
#' @noRd
stopProtImport <- function(prefix, error) {
    if (inherits(error, "prot_import_non_numeric_quantity")) {
        rlang::abort(paste0(prefix, conditionMessage(error)), class = "prot_import_non_numeric_quantity",
            column = error$column, examples = error$examples, call = NULL)
    }
    stop(prefix, conditionMessage(error))
}

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
  
  # Both quantity columns are checked whatever `use_precursor_norm` says: the analysis reads both (the S4 constructor
  # fixes raw_quantity_column = Precursor.Quantity for the IQ rollup and imputation, and norm_quantity_column =
  # Precursor.Normalised for the peptide filters and limpa; func_pept_s4_core.R).
  read_columns <- intersect(c("Precursor.Normalised", "Precursor.Quantity"), names(data))
  # P-05: missing-value tokens in a quantity column are read as NA with a notice; other text refuses the import.
  data <- checkProtImportQuantityColumns(data, read_columns, filepath)

  # Convert Precursor.Normalised if needed
  if (use_precursor_norm && "Precursor.Normalised" %in% names(data)) {
    log_info("Converting Precursor.Normalised to numeric")
    data$Precursor.Normalised <- as.numeric(data$Precursor.Normalised)
  }

  # A quantity of exactly 0 is a precursor DIA-NN did not quantify: a missing value, not a quantity, as the LFQ and TMT
  # readers read theirs (KIT-101 F-18/F-26) and as iq reads it (user, 2026-11-12, after the A101-07 adjudication).
  # It is read as NA with a notice of the count per column. A negative or infinite quantity is not a missing value
  # but a damaged file: it refuses the import, as the rollup does (rollUpPrecursorToPeptideHelper).
  # Memory: a file read by vroom keeps its columns on disk until they are used. Each column is scanned in blocks, so
  # no whole column is held for the scan; a column is rebuilt in memory only when it has zeros (8 bytes per row).
  zeros <- integer()
  rows <- nrow(data)
  block <- 500000L
  starts <- if (rows > 0L) seq.int(1L, rows, by = block) else integer()
  for (column in read_columns) {
    if (!is.numeric(data[[column]])) next
    zero <- integer()
    for (first in starts) {
      at <- first:min(first + block - 1L, rows)
      part <- data[[column]][at]
      damaged <- which(part < 0 | is.infinite(part))
      if (length(damaged)) {
        rlang::abort(sprintf(paste0("Quantity column '%s' contains negative or infinite value(s) (e.g. %s, row %d). ",
            "Quantities must be zero or positive; correct the file and import it again."),
            column, paste(utils::head(part[damaged], 3L), collapse = ", "), at[[damaged[[1L]]]]), call = NULL)
      }
      zero <- c(zero, at[which(part == 0)])
    }
    if (length(zero)) {
      zeros[[column]] <- length(zero)
      values <- numeric(rows)
      for (first in starts) {
        at <- first:min(first + block - 1L, rows)
        values[at] <- data[[column]][at]
      }
      values[zero] <- NA_real_
      data[[column]] <- values
      values <- NULL
    }
  }
  if (length(zeros)) {
    rlang::inform(sprintf("Zero quantities (precursors DIA-NN did not quantify) were read as NA: %s.",
        paste(sprintf("%d x %s", zeros, names(zeros)), collapse = ", ")),
        class = "prot_import_missing_tokens", counts = zeros)
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
  # KIT-102 parity: a text-numeric quantity column is read as numbers. P-05: missing-value tokens are read as NA with
  # a notice; other text refuses the import.
  data <- checkProtImportQuantityColumns(data, if (quantity_type == "pg") "PG.Quantity" else "PEP.Quantity", filepath)
  
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
  
  # P-05: missing-value tokens in the intensity columns are read as NA with a notice; other text refuses the import.
  data <- checkProtImportQuantityColumns(data, intensity_cols, filepath)

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
      # Ensure Intensity is numeric. FragPipe writes 0 where a protein was not quantified in a sample: that is a
      # missing value, not a quantity (KIT-101 F-18; kept, the design step's log2(x + 1) made it a log2 value of 0).
      dplyr::mutate(Intensity = dplyr::na_if(as.numeric(Intensity), 0)) |>
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
#' Map MaxQuant's spaced column names to the dotted names the reader handles
#'
#' @param names Column names as read from a proteinGroups export.
#' @return The names with the known MaxQuant names in dotted form; other names unchanged.
#' @noRd
normaliseMaxQuantColumnNames <- function(names) {
  fixed <- c("Protein IDs", "Majority protein IDs", "Protein names", "Gene names", "Potential contaminant",
    "Only identified by site", "Peptide counts (all)", "Peptide counts (razor+unique)", "Peptide counts (unique)",
    "Unique peptides", "MS/MS count")
  out <- names
  at <- out %in% fixed
  out[at] <- gsub(" ", ".", out[at], fixed = TRUE)
  prefixes <- c("LFQ intensity ", "Intensity ", "iBAQ ", "MS/MS count ", "Unique peptides ", "Peptides ")
  for (prefix in prefixes) {
    at <- startsWith(out, prefix)
    out[at] <- paste0(gsub(" ", ".", trimws(prefix), fixed = TRUE), ".", substring(out[at], nchar(prefix) + 1L))
  }
  out
}

#' Import MaxQuant Data
#' 
#' @param filepath Path to MaxQuant proteinGroups.txt file
#' @param use_lfq Whether to use LFQ intensities
#' @param filter_contaminants Whether to drop rows marked "+" in the contaminant column ("Potential contaminant", or "Contaminant" in MaxQuant 1.4 and earlier) and in "Reverse"; the rows removed per column are logged
#' @return List with data, data_type, and column_mapping
#' @export
importMaxQuantData <- function(filepath, use_lfq = TRUE, filter_contaminants = TRUE) {
  # The flag columns are typed explicitly: vroom guesses a column's type from a sample of rows across the file,
  # and a sparse "+" column whose sample carries no mark is read as logical with every mark lost, so the rows it
  # flags are kept (reproduced by the KIT-021 corpus check: 50 decoys among 8,100 groups survived the import).
  header <- names(vroom::vroom(filepath, delim = "\t", n_max = 0L, show_col_types = FALSE))
  flags <- intersect(c("Potential.contaminant", "Reverse", "Potential contaminant", "Contaminant"), header)
  col_types <- do.call(vroom::cols, c(list(.default = vroom::col_guess()),
    stats::setNames(rep(list(vroom::col_character()), length(flags)), flags)))
  data <- vroom::vroom(filepath, delim = "\t", show_col_types = FALSE, col_types = col_types)
  # P-05: missing-value tokens in the intensity columns are read as NA with a notice; other text refuses the import.
  # Checked under the file's own names, before the rows are filtered (the choice of columns is the one made below).
  dotted <- normaliseMaxQuantColumnNames(names(data))
  lfq <- grepl("^LFQ\\.intensity\\.", dotted)
  data <- checkProtImportQuantityColumns(data,
    names(data)[if (use_lfq && any(lfq)) lfq else grepl("^Intensity\\.", dotted)], filepath)
  # Genuine MaxQuant output writes its column names with spaces ("LFQ intensity <experiment>", "Potential
  # contaminant", "Majority protein IDs"); the reader's tests and fixtures carry the dotted form R's readers once
  # produced. The known MaxQuant names are mapped to the dotted form (the experiment part of an intensity column is
  # left as it is), so a real proteinGroups.txt reads as the fixtures do (KIT-054 finding).
  names(data) <- normaliseMaxQuantColumnNames(names(data))
  
  # Filter contaminants and reverse hits if requested. MaxQuant 1.4 and earlier name the contaminant column
  # "Contaminant" (from 1.5.2.8 at the latest "Potential contaminant"); both are filtered (KIT-171). The whole group
  # goes when MaxQuant marked it; accession prefixes are not consulted. A row marked in two columns is counted at the
  # first, so the logged counts sum to the rows removed.
  if (filter_contaminants) {
    if (!any(c("Potential.contaminant", "Contaminant") %in% names(data))) {
      log_info("MaxQuant import: contaminant filtering requested but no contaminant column found ('Potential contaminant' or 'Contaminant'); no contaminant rows removed")
    }
    for (column in intersect(c("Potential.contaminant", "Contaminant", "Reverse"), names(data))) {
      keep <- is.na(data[[column]]) | data[[column]] != "+"
      log_info(sprintf("MaxQuant import: removed %d rows marked '+' in %s", sum(!keep), column))
      data <- data[keep, ]
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

  if (length(quantity_cols) == 0) {
    stop("No MaxQuant intensity columns found. Expected LFQ.intensity.* or Intensity.* columns.")
  }

  protein_col_candidates <- c(
    "Majority.protein.IDs",
    "Protein.IDs",
    "Protein.Ids",
    "Protein IDs",
    "Protein ID"
  )
  protein_col <- protein_col_candidates[protein_col_candidates %in% names(data)][1]
  if (is.na(protein_col)) {
    stop(
      "Required MaxQuant protein identifier column not found. Expected one of: ",
      paste(protein_col_candidates, collapse = ", ")
    )
  }

  cols_to_keep <- unique(c(
    protein_col,
    intersect(c("Gene.names", "Protein.names"), names(data)),
    quantity_cols
  ))

  long_data <- tryCatch({
    data |>
      dplyr::select(dplyr::all_of(cols_to_keep)) |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(quantity_cols),
        names_to = "Run",
        values_to = "Intensity"
      ) |>
      dplyr::mutate(
        Run = sub("^(LFQ\\.intensity\\.|Intensity\\.)", "", Run),
        # MaxQuant writes 0 where a protein group was not quantified in an experiment: a missing value, not a
        # quantity (KIT-101 F-18).
        Intensity = dplyr::na_if(as.numeric(Intensity), 0)
      ) |>
      dplyr::rename(Protein.Ids = !!rlang::sym(protein_col))
  }, error = function(e) {
    stop("Failed to convert MaxQuant proteinGroups data to long format: ", e$message)
  })
  
  return(list(
    data = long_data,
    data_type = "protein",  # MaxQuant proteinGroups is protein-level
    column_mapping = list(
      protein_col = "Protein.Ids",
      peptide_col = NULL,
      run_col = "Run",
      quantity_col = "Intensity",
      qvalue_col = "Q.value"  # If Perseus was used
    )
  ))
}

# ----------------------------------------------------------------------------
# importProteomicsData
# ----------------------------------------------------------------------------
#' Read a proteomics search-results file with the reader for its format.
#' @param format One of "diann", "spectronaut", "fragpipe", "maxquant" or "pd_tmt".
#' @param path The search-results file.
#' @param options Reader options as a list: `use_precursor_norm` (DIA-NN), `quantity_type` (Spectronaut),
#'   `use_maxlfq` (FragPipe), `use_lfq` and `filter_contaminants` (MaxQuant). A missing option takes the
#'   GUI's default.
#' @param readers The reader for each format, by format name; replaceable in tests.
#' @param logError Error logger.
#' @return The reader's result: `data`, `data_type` and `column_mapping`. A failure is re-signalled as
#'   "Failed to import data: <reason>".
#' @noRd
importProteomicsData <- function(
    format
    , path
    , options = list()
    , readers = list(
        diann = importDIANNData
        , spectronaut = importSpectronautData
        , fragpipe = importFragPipeData
        , maxquant = importMaxQuantData
        , pd_tmt = importProteomeDiscovererTMTData
      )
    , logError = logger::log_error
) {
  tryCatch({
    switch(format,
      "diann" = readers$diann(
        filepath = path,
        use_precursor_norm = options$use_precursor_norm %||% TRUE
      ),
      "spectronaut" = readers$spectronaut(
        filepath = path,
        quantity_type = options$quantity_type %||% "pg"
      ),
      "fragpipe" = readers$fragpipe(
        filepath = path,
        use_maxlfq = options$use_maxlfq %||% TRUE
      ),
      "maxquant" = readers$maxquant(
        filepath = path,
        use_lfq = options$use_lfq %||% TRUE,
        filter_contaminants = options$filter_contaminants %||% TRUE
      ),
      "pd_tmt" = readers$pd_tmt(
        filepath = path
      ),
      stop("Unsupported format: ", format)
    )
  }, error = function(e) {
    logError(paste("Error in data import function:", e$message))
    stopProtImport("Failed to import data: ", e)
  })
}
