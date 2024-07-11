

## Create S4 class for protomics protein level abundance data
ProteinsQuantitativeData <- setClass("ProteinsQuantitativeData"
         
         , slots = c(
                      # Protein vs Sample quantitative data
                      protein_data = "data.frame"
                      , protein_id_column = "character"

                      # Design Matrix Information
                      , design_matrix = "data.frame"
                      , sample_id="character"
                      , group_id="character"
                      , technical_replicate_id="character"

                    )
         
         , prototype = list( 
           # Protein vs Sample quantitative data
           protein_id_column = "character"

           # Design Matrix Information
           , sample_id="Sample_id"
           , group_id="group"
           , technical_replicate_id="replicates"

           
           
           )
         
         , validity = function(object) {
           if( !is.data.frame(object@protein_data) ) {
             stop("protein_data must be a data.frame")
           }
           if( !is.character(object@protein_id_column) ) {
             stop("protein_id_column must be a character")
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

            if( ! object@protein_id_column %in% colnames(object@protein_data) ) {
                stop("Protein ID column must be in the protein data table")
            }
           
           if( ! object@sample_id %in% colnames(object@design_matrix) ) {
             stop("Sample ID column must be in the design matrix")
           }
           
           
           #Need to check the rows names in design matrix and the column names of the data table
           samples_in_protein_data <- setdiff(colnames( object@protein_data), object@protein_id_column) 
           samples_in_design_matrix <- object@design_matrix |> pull( !! sym( object@sample_id ) )
  
           if( length( which( sort(samples_in_protein_data) != sort(samples_in_design_matrix) )) > 0 ) {
             stop("Samples in protein data and design matrix must be the same" )
           }
           
         }

)


setGeneric( name ="setProteinData"
            , def=function( theObject, protein_data, protein_id_column) {
                standardGeneric("setProteinData")
            })

setMethod( f ="setProteinData"
           , signature = "ProteinsQuantitativeData"
            , definition=function( theObject, protein_data, protein_id_column ) {
              theObject@protein_data <- protein_data
              theObject@protein_id_column <- protein_id_column
              
              return(theObject)
            })

setGeneric(name="proteinIntensityFiltering"
           , def=function( theObject, min_protein_intensity_percentile, proportion_samples_below_intensity_threshold, cluster) {
             standardGeneric("proteinIntensityFiltering")
           })

setMethod( f="proteinIntensityFiltering"
           , signature="ProteinsQuantitativeData"
           , definition = function( theObject, min_protein_intensity_percentile, proportion_samples_below_intensity_threshold, cluster) {
             protein_data <- theObject@protein_data

             data_long_cln <- protein_data  |>
               pivot_longer( cols=!matches(theObject@protein_id_column)
                             , names_to = "Sample_ID"
                             , values_to = "log_values")  |>
               mutate( temp = "")
             
             print(head( data_long_cln ))
             
             min_peptide_intensity_threshold <- ceiling( quantile( data_long_cln$log_values, na.rm=TRUE, probs = c(min_protein_intensity_percentile) ))[1]
             
             peptide_normalized_pif_cln <- peptideIntensityFiltering( data_long_cln
                                                                      , min_peptide_intensity_threshold = min_peptide_intensity_threshold
                                                                      , proportion_samples_below_intensity_threshold = proportion_samples_below_intensity_threshold
                                                                      , protein_id_column = !!sym( theObject@protein_id_column)
                                                                      , peptide_sequence_column = temp
                                                                      , peptide_quantity_column = log_values
                                                                      , cluster = cluster)
             
             
             theObject@protein_data <- peptide_normalized_pif_cln |>
               dplyr::select( -temp) |>
               pivot_wider( id_cols = theObject@protein_id_column , names_from = Sample_ID, values_from = log_values) 
             
             updated_object <- theObject
          return(updated_object)   
          })




setGeneric(name="removeProteinsWithOnlyOneReplicateObj"
           , def=function( theObject, cluster ) {
             standardGeneric("removeProteinsWithOnlyOneReplicateObj")
           })


setMethod(f="removeProteinsWithOnlyOneReplicateObj"
          , signature="ProteinsQuantitativeData"
          , definition=function( theObject, cluster) {
            protein_data <- theObject@protein_data
            samples_id_tbl <- theObject@design_matrix
            sample_id_tbl_sample_id_column <- theObject@sample_id
            replicate_group_column <- theObject@technical_replicate_id
            protein_id_column <- theObject@protein_id_column
            
            input_table_sample_id_column <- 'Sample_ID'
            quantity_column <- "log_values"
              
            data_long_cln <- protein_data  |>
              pivot_longer( cols=!matches(protein_id_column)
                            , names_to = input_table_sample_id_column
                            , values_to = quantity_column)  
            
            protein_data <- removeProteinsWithOnlyOneReplicate( input_table = data_long_cln
                                                                , samples_id_tbl = samples_id_tbl
                                                                , cluster = cluster
                                                                , input_table_sample_id_column = !!sym( input_table_sample_id_column )
                                                                , sample_id_tbl_sample_id_column = !!sym( sample_id_tbl_sample_id_column)
                                                                , replicate_group_column = !!sym( replicate_group_column)
                                                                , protein_id_column = !!sym( protein_id_column)
                                                                , quantity_column = !!sym( quantity_column))
            
            
            theObject@protein_data <- protein_data |>
              pivot_wider( id_cols = !!sym( protein_id_column)
                           , names_from = !!sym( input_table_sample_id_column)
                           , values_from = !!sym( quantity_column) )
            
            updated_object <- theObject
            
            return(updated_object)
          })

