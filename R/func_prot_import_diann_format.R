# ----------------------------------------------------------------------------
# formatDIANN
# ----------------------------------------------------------------------------
#' @title Format DIANN Data from EList
#' @description Converts a limpa EList object (from readDIANN) into a long-format data frame 
#' compatible with MultiScholaR proteomics import.
#' @param data_tbl_parquet_filt An EList object containing proteomics data.
#' @return A long-format data frame with Protein.Ids, Stripped.Sequence, Run, and Peptide.Normalised columns.
#' @importFrom dplyr mutate select filter rename distinct left_join
#' @importFrom tidyr pivot_longer separate
#' @importFrom tibble rownames_to_column
#' @export
formatDIANN <- function(data_tbl_parquet_filt) {
  
  # Extract intensity matrix and gene annotations
  intensity_matrix <- data_tbl_parquet_filt$E
  gene_info <- data_tbl_parquet_filt$genes
  
  # Convert intensity matrix to data frame with precursor IDs
  intensity_df <- as.data.frame(intensity_matrix) |>
    rownames_to_column(var = "Precursor.Id")

  # Ensure Precursor.Id is present in gene_info (needed for join)
  if (!"Precursor.Id" %in% names(gene_info)) {
    # If genes info matches intensity matrix dimensions, use intensity matrix rownames
    # which we know correspond to Precursor.Id
    if (nrow(gene_info) == nrow(intensity_matrix)) {
      gene_info$Precursor.Id <- rownames(intensity_matrix)
    }
  }

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
#' @title Format DIANN Parquet Data
#' @description Converts a limpa EList object containing parquet data into a long-format 
#' data frame compatible with MultiScholaR proteomics import.
#' @param data_tbl_parquet_filt An EList object containing proteomics data.
#' @return A long-format data frame.
#' @export
formatDIANNParquet <- function(data_tbl_parquet_filt) {
  return(formatDIANN(data_tbl_parquet_filt))
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

