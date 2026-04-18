# =========================================-
# Content from peptideVsSamplesS4Objects.R
# ==========================================
## Create S4 class for proteomics protein level abundance data
#'@exportClass PeptideQuantitativeData
PeptideQuantitativeData <- setClass("PeptideQuantitativeData"

                                     , slots = c(
                                       # Protein vs Sample quantitative data
                                       peptide_data = "data.frame"
                                       , peptide_matrix = "matrix"
                                       , protein_id_column = "character"
                                       , peptide_sequence_column = "character"
                                       , q_value_column = "character"
                                       , global_q_value_column = "character"
                                       , proteotypic_peptide_sequence_column = "character"
                                       , raw_quantity_column = "character"
                                       , norm_quantity_column = "character"
                                       , is_logged_data = "logical"

                                       # Design Matrix Information
                                       , design_matrix = "data.frame"
                                       , sample_id="character"
                                       , group_id="character"
                                       , technical_replicate_id="character"
                                       , args = "list"
                                     )

                                     , prototype = list(
                                       # Protein vs Sample quantitative data
                                       peptide_matrix = matrix(numeric(0), nrow = 0, ncol = 0)
                                       , protein_id_column = "Protein_Ids"
                                       , peptide_sequence_column = "Stripped.Sequence"
                                       , q_value_column = "Q.Value"
                                       , global_q_value_column = "Global.Q.Value"
                                       , proteotypic_peptide_sequence_column = "Proteotypic"
                                       , raw_quantity_column = "Precursor.Quantity"
                                       , norm_quantity_column = "Precursor.Normalised"
                                       , is_logged_data = FALSE

                                       # Design Matrix Information
                                       , sample_id="Sample_id"
                                       , group_id="group"
                                       , technical_replicate_id="replicates"

                                       # Parameters for methods and functions
                                       , args = list()

                                       )
                                     , validity = function(object) {
                                       if( !is.data.frame(object@peptide_data) ) {
                                         stop("peptide_data must be a data.frame")
                                       }
                                       if( !is.character(object@protein_id_column) ) {
                                         stop("protein_id_column must be a character")
                                       }

                                       if( !is.character(object@peptide_sequence_column) ) {
                                         stop("peptide_sequence_column must be a character")
                                       }
                                       if( !is.character(object@q_value_column) ) {
                                         stop("q_value_column must be a character")
                                       }
                                       if(!is.character(object@global_q_value_column) ) {
                                         stop("global_q_value_column must be a character")
                                       }
                                       if(!is.character(object@proteotypic_peptide_sequence_column) ) {
                                         stop("proteotypic_peptide_sequence_column must be a character")
                                       }
                                       if(!is.character(object@raw_quantity_column)  ) {
                                         stop("raw_quantity_column must be a character")
                                       }

                                       if(!is.character(object@norm_quantity_column) ) {
                                         stop("norm_quantity_column must be a character")
                                       }

                                       if( !is.data.frame(object@design_matrix) ) {
                                         stop("design_matrix must be a data.frame")
                                       }

                                       if( !is.character(object@sample_id) ) {
                                         stop("sample_id must be a character")
                                       }

                                       if( !is.character(object@group_id) ) {
                                         stop("group_id must be a character")
                                       }

                                       if( !is.character(object@technical_replicate_id) ) {
                                         stop("technical_replicate_id must be a character")
                                       }

                                       if( ! object@protein_id_column %in% colnames(object@peptide_data) ) {
                                         print(protein_id_column)
                                         print( colnames(object@peptide_data) )
                                         stop("Protein ID column must be in the peptide data table")
                                       }

                                       if(!object@peptide_sequence_column %in% colnames(object@peptide_data) ) {
                                         stop("Peptide sequence column must be in the peptide data table")
                                       }

                                       if(!object@q_value_column %in% colnames(object@peptide_data) ) {
                                         stop("Q value column must be in the peptide data table")
                                       }

                                       if(!object@raw_quantity_column %in% colnames(object@peptide_data) &
                                          !object@norm_quantity_column %in% colnames(object@peptide_data) ) {
                                         stop("Precursor raw quantity or normalised quantity column must be in the peptide data table")
                                       }

                                       if( ! object@sample_id %in% colnames(object@design_matrix) ) {
                                         stop("Sample ID column must be in the design matrix")
                                       }


                                       #Need to check the rows names in design matrix and the column names of the data table
                                       samples_in_peptide_data <-  object@peptide_data |> distinct(!!sym(object@sample_id)) |> dplyr::pull(!!sym(object@sample_id))
                                       samples_in_design_matrix <- object@design_matrix |> dplyr::pull( !! sym( object@sample_id ) )

                                       if( length( which( sort(samples_in_peptide_data) != sort(samples_in_design_matrix) )) > 0 ) {
                                         stop("Samples in peptide data and design matrix must be the same" )
                                       }

                                     }

)

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

#'@export
setMethod(f="calcPeptideMatrix"
          , signature="PeptideQuantitativeData"
          , definition=function( theObject ) {
            
            peptide_data <- theObject@peptide_data
            # Use the NORMALIZED quantity column for the matrix, not the raw one.
            # This is critical for all downstream stats (limpa, etc.)
            quantity_column_to_use <- theObject@norm_quantity_column 
            sample_id_column <- theObject@sample_id
            protein_id_column <- theObject@protein_id_column
            peptide_sequence_column <- theObject@peptide_sequence_column

            peptide_quant_table <- peptide_data |>
              dplyr::select(!!sym(sample_id_column)
                            , !!sym(theObject@protein_id_column)
                            , !!sym(theObject@peptide_sequence_column)
                            , !!sym(quantity_column_to_use)) |>
              mutate( !!sym(quantity_column_to_use)  := purrr::map_dbl(!!sym(quantity_column_to_use), as.numeric )) |>
              tidyr::pivot_wider(names_from = !!sym(sample_id_column)
                                 , values_from = !!sym(quantity_column_to_use)) |>
              dplyr::rename(Protein.Ids = !!sym(theObject@protein_id_column)
                            , Stripped.Sequence = !!sym(theObject@peptide_sequence_column))
            
            
            normalised_frozen_peptide_matrix_filt <- peptide_quant_table |>
              mutate(Protein.Ids = as.character(Protein.Ids)
                     , Stripped.Sequence = as.character(Stripped.Sequence) ) |>
              mutate(peptide_ids = paste(Protein.Ids, Stripped.Sequence, sep = "%")) |>
              dplyr::select(-Protein.Ids,-Stripped.Sequence) |>
              column_to_rownames("peptide_ids") |>
              as.matrix()           
            
            theObject@peptide_matrix <- normalised_frozen_peptide_matrix_filt
                                                   
            # print(head( normalised_frozen_peptide_matrix_filt))                                                   
            
            return(theObject)
            
          } )           

