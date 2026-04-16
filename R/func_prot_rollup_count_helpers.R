# ----------------------------------------------------------------------------
# calcPeptidesPerProtein
# ----------------------------------------------------------------------------
#' Calculate peptides per protein
#' 
#' @param data A data frame or S4 object containing protein/peptide data
#' @return Data frame with protein IDs and peptide counts
#' @export
calcPeptidesPerProtein <- function(data) {
  # For protein quantification data, return empty data frame
  if (isS4(data)) {
    if ("protein_quant_table" %in% slotNames(data)) {
      return(data.frame(Protein.Ids = character(), 
                       n_peptides = integer()))
    }
    if ("peptide_data" %in% slotNames(data)) {
      return(data@peptide_data |>
             group_by(Protein.Ids) |>
             summarise(n_peptides = n_distinct(Stripped.Sequence), 
                      .groups = "drop"))
    }
  }
  
  # For regular dataframes, check if it's protein quantification data
  if ("Protein.Ids" %in% names(data)) {
    if (all(sapply(data[setdiff(names(data), "Protein.Ids")], is.numeric))) {
      return(data.frame(Protein.Ids = character(), 
                       n_peptides = integer()))
    }
    
    if ("Stripped.Sequence" %in% names(data)) {
      return(data |>
             group_by(Protein.Ids) |>
             summarise(n_peptides = n_distinct(Stripped.Sequence), 
                      .groups = "drop"))
    }
  }
  stop("Required columns not found")
}

# ----------------------------------------------------------------------------
# calcTotalPeptides
# ----------------------------------------------------------------------------
#' Calculate total unique peptides
#' 
#' @param data A data frame or S4 object containing protein/peptide data
#' @return Integer count of unique peptide-protein combinations
#' @export
calcTotalPeptides <- function(data) {
  # For protein quantification data, return NA
  if (isS4(data)) {
    if ("protein_quant_table" %in% slotNames(data)) {
      return(NA_integer_)
    }
    if ("peptide_data" %in% slotNames(data)) {
      return(data@peptide_data |>
             distinct(Protein.Ids, Stripped.Sequence) |>
             nrow())
    }
  }
  
  # For regular dataframes, check if it's protein quantification data
  if ("Protein.Ids" %in% names(data)) {
    if (all(sapply(data[setdiff(names(data), "Protein.Ids")], is.numeric))) {
      return(NA_integer_)
    }
    
    if ("Stripped.Sequence" %in% names(data)) {
      return(distinct(data, Protein.Ids, Stripped.Sequence) |> nrow())
    }
  }
  stop("Required columns not found")
}

# ----------------------------------------------------------------------------
# countPeptidesPerRun
# ----------------------------------------------------------------------------
#' Count peptides per run
#' 
#' @param data A data frame or S4 object containing protein/peptide data
#' @return Data frame with run IDs and peptide counts
#' @export
countPeptidesPerRun <- function(data) {
  # For protein quantification data, return empty data frame
  if (isS4(data)) {
    if ("protein_quant_table" %in% slotNames(data)) {
      return(data.frame(Run = character(), 
                       n_peptides = integer()))
    }
    if ("peptide_data" %in% slotNames(data)) {
      return(data@peptide_data |>
             group_by(Run) |>
             summarise(n_peptides = n_distinct(Stripped.Sequence), 
                      .groups = "drop") |>
             arrange(Run))
    }
  }
  
  # For regular dataframes, check if it's protein quantification data
  if ("Protein.Ids" %in% names(data)) {
    if (all(sapply(data[setdiff(names(data), "Protein.Ids")], is.numeric))) {
      return(data.frame(Run = character(), 
                       n_peptides = integer()))
    }
    
    if (all(c("Run", "Stripped.Sequence") %in% names(data))) {
      return(data |>
             group_by(Run) |>
             summarise(n_peptides = n_distinct(Stripped.Sequence), 
                      .groups = "drop") |>
             arrange(Run))
    }
  }
  stop("Required columns not found")
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
#' @return Data frame with run IDs and protein counts
#' @export
countProteinsPerRun <- function(data) {
  if (isS4(data)) {
    if ("peptide_data" %in% slotNames(data)) {
      return(data@peptide_data |>
             group_by(Run) |>
             summarise(n_proteins = n_distinct(Protein.Ids), 
                      .groups = "drop") |>
             arrange(Run))
    }
    if ("protein_quant_table" %in% slotNames(data)) {
      data <- data@protein_quant_table
      run_cols <- setdiff(names(data), "Protein.Ids")
      
      # For each run (column), count non-NA values
      result <- data.frame(
        Run = run_cols,
        n_proteins = sapply(run_cols, function(col) {
          sum(!is.na(data[[col]]))
        })
      ) |> arrange(Run)
      
      return(result)
    }
  }
  
  # For regular dataframes
  if ("Protein.Ids" %in% names(data)) {
    # Check if it's a protein quantification table
    if (all(sapply(data[setdiff(names(data), "Protein.Ids")], is.numeric))) {
      run_cols <- setdiff(names(data), "Protein.Ids")
      
      # For each run (column), count non-NA values
      result <- data.frame(
        Run = run_cols,
        n_proteins = sapply(run_cols, function(col) {
          sum(!is.na(data[[col]]))
        })
      ) |> arrange(Run)
      
      return(result)
    }
    
    # For peptide data
    if ("Run" %in% names(data)) {
      return(data |>
             group_by(Run) |>
             summarise(n_proteins = n_distinct(Protein.Ids), 
                      .groups = "drop") |>
             arrange(Run))
    }
  }
  stop("Required columns not found")
}

# ----------------------------------------------------------------------------
# countUniqueProteins
# ----------------------------------------------------------------------------
#' Count unique proteins in peptide or protein data
#' 
#' @param data A data frame or S4 object containing protein/peptide data
#' @return Integer count of unique proteins
#' @export
countUniqueProteins <- function(data) {
  if (isS4(data)) {
    if ("peptide_data" %in% slotNames(data)) {
      return(data@peptide_data |> 
             distinct(Protein.Ids) |> 
             nrow())
    }
    if ("protein_quant_table" %in% slotNames(data)) {
      return(nrow(data@protein_quant_table))
    }
  }
  
  # For regular dataframes
  if ("Protein.Ids" %in% names(data)) {
    # Check if it's a protein quantification table
    if (all(sapply(data[setdiff(names(data), "Protein.Ids")], is.numeric))) {
      return(nrow(data))  # Each row is a unique protein
    }
    return(distinct(data, Protein.Ids) |> nrow())
  }
  stop("No Protein.Ids column found")
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

