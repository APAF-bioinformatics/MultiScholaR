resolveProteinCountColumn <- function(data, protein_id_column = NULL) {
  if (isS4(data)) {
    data_slots <- slotNames(data)
    if (is.null(protein_id_column) && "protein_id_column" %in% data_slots) {
      protein_id_column <- data@protein_id_column
    }
    if ("peptide_data" %in% data_slots) {
      return(resolveProteinCountColumn(data@peptide_data, protein_id_column))
    }
    if ("protein_quant_table" %in% data_slots) {
      return(resolveProteinCountColumn(data@protein_quant_table, protein_id_column))
    }
    stop("No protein or peptide data slot found.", call. = FALSE)
  }

  if (!is.data.frame(data)) {
    stop("Protein counting requires a data frame or supported S4 object.", call. = FALSE)
  }

  if (!is.null(protein_id_column)) {
    protein_id_column <- as.character(protein_id_column)[1]
    if (!is.na(protein_id_column) && nzchar(protein_id_column) && protein_id_column %in% names(data)) {
      return(protein_id_column)
    }
    stop("The declared protein identity column is absent from the data.", call. = FALSE)
  }

  candidates <- c("Protein.Group", "Protein.Ids", "Protein.IDs", "protein_id")
  resolved <- candidates[candidates %in% names(data)][1]
  if (is.na(resolved)) {
    stop("No protein identity column found.", call. = FALSE)
  }
  resolved
}

isProteinQuantificationTable <- function(data, protein_id_column = NULL) {
  if (isS4(data)) {
    return("protein_quant_table" %in% slotNames(data))
  }
  protein_id_column <- resolveProteinCountColumn(data, protein_id_column)
  value_columns <- setdiff(names(data), protein_id_column)
  length(value_columns) > 0L && all(vapply(data[value_columns], is.numeric, logical(1)))
}

# ----------------------------------------------------------------------------
# calcPeptidesPerProtein
# ----------------------------------------------------------------------------
#' Calculate peptides per protein
#'
#' @param data A data frame or S4 object containing protein/peptide data
#' @param protein_id_column Optional declared protein identity column. The active
#'   protein key is inferred when `NULL`.
#' @return Data frame with protein IDs and peptide counts
#' @export
calcPeptidesPerProtein <- function(data, protein_id_column = NULL) {
  # For protein quantification data, return empty data frame
  if (isS4(data)) {
    if ("protein_quant_table" %in% slotNames(data)) {
      return(data.frame(Protein.Ids = character(), 
                       n_peptides = integer()))
    }
    if ("peptide_data" %in% slotNames(data)) {
      return(calcPeptidesPerProtein(data@peptide_data, data@protein_id_column))
    }
  }

  protein_id_column <- resolveProteinCountColumn(data, protein_id_column)
  if (isProteinQuantificationTable(data, protein_id_column)) {
    return(data.frame(Protein.Ids = character(), n_peptides = integer()))
  }
  if ("Stripped.Sequence" %in% names(data)) {
    result <- data |>
      group_by(!!rlang::sym(protein_id_column)) |>
      summarise(n_peptides = n_distinct(Stripped.Sequence), .groups = "drop")
    names(result)[1] <- "Protein.Ids"
    return(result)
  }
  stop("Required peptide sequence column not found.", call. = FALSE)
}

# ----------------------------------------------------------------------------
# calcTotalPeptides
# ----------------------------------------------------------------------------
#' Calculate total unique peptides
#'
#' @param data A data frame or S4 object containing protein/peptide data
#' @param protein_id_column Optional declared protein identity column. The active
#'   protein key is inferred when `NULL`.
#' @return Integer count of unique peptide-protein combinations
#' @export
calcTotalPeptides <- function(data, protein_id_column = NULL) {
  # For protein quantification data, return NA
  if (isS4(data)) {
    if ("protein_quant_table" %in% slotNames(data)) {
      return(NA_integer_)
    }
    if ("peptide_data" %in% slotNames(data)) {
      return(calcTotalPeptides(data@peptide_data, data@protein_id_column))
    }
  }

  protein_id_column <- resolveProteinCountColumn(data, protein_id_column)
  if (isProteinQuantificationTable(data, protein_id_column)) {
    return(NA_integer_)
  }
  if ("Stripped.Sequence" %in% names(data)) {
    return(data |>
      distinct(!!rlang::sym(protein_id_column), Stripped.Sequence) |>
      nrow())
  }
  stop("Required peptide sequence column not found.", call. = FALSE)
}

# ----------------------------------------------------------------------------
# countPeptidesPerRun
# ----------------------------------------------------------------------------
#' Count peptides per run
#'
#' @param data A data frame or S4 object containing protein/peptide data
#' @param protein_id_column Optional declared protein identity column. The active
#'   protein key is inferred when `NULL`.
#' @return Data frame with run IDs and peptide counts
#' @export
countPeptidesPerRun <- function(data, protein_id_column = NULL) {
  # For protein quantification data, return empty data frame
  if (isS4(data)) {
    if ("protein_quant_table" %in% slotNames(data)) {
      return(data.frame(Run = character(), 
                       n_peptides = integer()))
    }
    if ("peptide_data" %in% slotNames(data)) {
      return(countPeptidesPerRun(data@peptide_data, data@protein_id_column))
    }
  }

  protein_id_column <- resolveProteinCountColumn(data, protein_id_column)
  if (isProteinQuantificationTable(data, protein_id_column)) {
    return(data.frame(Run = character(), n_peptides = integer()))
  }
  if (all(c("Run", "Stripped.Sequence") %in% names(data))) {
    return(data |>
      group_by(Run) |>
      summarise(n_peptides = n_distinct(Stripped.Sequence), .groups = "drop") |>
      arrange(Run))
  }
  stop("Required run or peptide sequence column not found.", call. = FALSE)
}

# ----------------------------------------------------------------------------
# count_num_peptides
# ----------------------------------------------------------------------------
# Count the number of peptides in the input table
#' @export
count_num_peptides <- function( input_table
                                , protein_id_column = Protein.Ids
                                , peptide_sequence_column = Stripped.Sequence ) {
  num_peptides <- input_table |>
    distinct( {{protein_id_column}}, {{peptide_sequence_column}}) |>
    count()

  num_peptides[[1,1]]
}

# ----------------------------------------------------------------------------
# countProteinsPerRun
# ----------------------------------------------------------------------------
#' Count proteins per run
#'
#' @param data A data frame or S4 object containing protein/peptide data
#' @param protein_id_column Optional declared protein identity column. The active
#'   protein key is inferred when `NULL`.
#' @return Data frame with run IDs and protein counts
#' @export
countProteinsPerRun <- function(data, protein_id_column = NULL) {
  if (isS4(data)) {
    if ("peptide_data" %in% slotNames(data)) {
      return(countProteinsPerRun(data@peptide_data, data@protein_id_column))
    }
    if ("protein_quant_table" %in% slotNames(data)) {
      return(countProteinsPerRun(data@protein_quant_table, data@protein_id_column))
    }
  }

  protein_id_column <- resolveProteinCountColumn(data, protein_id_column)
  if (isProteinQuantificationTable(data, protein_id_column)) {
    run_columns <- setdiff(names(data), protein_id_column)
    return(data.frame(
      Run = run_columns,
      n_proteins = unname(vapply(
        run_columns,
        function(column) sum(!is.na(data[[column]])),
        integer(1)
      )),
      stringsAsFactors = FALSE
    ) |>
      arrange(Run))
  }
  if ("Run" %in% names(data)) {
    return(data |>
      group_by(Run) |>
      summarise(
        n_proteins = n_distinct(!!rlang::sym(protein_id_column)),
        .groups = "drop"
      ) |>
      arrange(Run))
  }
  stop("Required run column not found.", call. = FALSE)
}

# ----------------------------------------------------------------------------
# countUniqueProteins
# ----------------------------------------------------------------------------
#' Count unique proteins in peptide or protein data
#'
#' @param data A data frame or S4 object containing protein/peptide data
#' @param protein_id_column Optional declared protein identity column. The active
#'   protein key is inferred when `NULL`.
#' @return Integer count of unique proteins
#' @export
countUniqueProteins <- function(data, protein_id_column = NULL) {
  if (isS4(data)) {
    if ("peptide_data" %in% slotNames(data)) {
      return(countUniqueProteins(data@peptide_data, data@protein_id_column))
    }
    if ("protein_quant_table" %in% slotNames(data)) {
      return(nrow(data@protein_quant_table))
    }
  }

  protein_id_column <- resolveProteinCountColumn(data, protein_id_column)
  if (isProteinQuantificationTable(data, protein_id_column)) {
    return(nrow(data))
  }
  data |>
    distinct(!!rlang::sym(protein_id_column)) |>
    nrow()
}

# ----------------------------------------------------------------------------
# count_num_proteins
# ----------------------------------------------------------------------------
# Count the number of peptides in the input table
#' @export
count_num_proteins <- function( input_table
                                , protein_id_column = Protein.Ids) {
  num_proteins <- input_table |>
    distinct( {{protein_id_column}}) |>
    count()

  num_proteins[[1,1]]
}

# ----------------------------------------------------------------------------
# count_num_samples
# ----------------------------------------------------------------------------
# Count the number of samples in the input table
#' @export
count_num_samples <- function( input_table
                               , sample_id_column = Run) {
  num_samples <- input_table |>
    distinct( {{sample_id_column}}) |>
    count()

  num_samples[[1,1]]
}
