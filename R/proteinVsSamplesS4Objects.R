

## Create S4 class for protomics protein level abundance data
#'@export
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


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
setGeneric( name ="setProteinData"
            , def=function( theObject, protein_data, protein_id_column) {
                standardGeneric("setProteinData")
            })

#'@export
setMethod( f ="setProteinData"
           , signature = "ProteinsQuantitativeData"
            , definition=function( theObject, protein_data, protein_id_column ) {
              theObject@protein_data <- protein_data
              theObject@protein_id_column <- protein_id_column

              return(theObject)
            })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Format the design matrix so that only metadata for samples in the protein data are retained, and also
# sort the sample IDs in the same order as the data matrix
#'@export
setGeneric(name ="cleanDesignMatrixObj"
           , def=function( theObject) {
             standardGeneric("cleanDesignMatrixObj")
           })

#'@export
setMethod( f ="cleanDesignMatrixObj"
           , signature = "ProteinsQuantitativeData"
           , definition=function( theObject ) {

            samples_id_vector <- setdiff(colnames(theObject@protein_data), theObject@sample_id )

             theObject@design_matrix <- data.frame( temp_sample_id = samples_id_vector )  |>
               inner_join( theObject@design_matrix
                          , by = join_by ( temp_sample_id == !!sym(theObject@sample_id)) ) |>
               dplyr::rename( !!sym(theObject@sample_id) := "temp_sample_id" ) |>
               dplyr::filter( !!sym( theObject@sample_id) %in% samples_id_vector )


             return(theObject)
           })
##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
setGeneric(name="proteinIntensityFilteringObj"
           , def=function( theObject, min_protein_intensity_percentile, proportion_samples_below_intensity_threshold, cluster) {
             standardGeneric("proteinIntensityFilteringObj")
           })

#'@export
setMethod( f="proteinIntensityFilteringObj"
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

             theObject <- cleanDesignMatrixObj(theObject)

             updated_object <- theObject

          return(updated_object)
          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
setGeneric(name="removeProteinsWithOnlyOneReplicateObj"
           , def=function( theObject, cluster ) {
             standardGeneric("removeProteinsWithOnlyOneReplicateObj")
           }
           , signature=c("theObject", "cluster"))

#'@export
setMethod(f="removeProteinsWithOnlyOneReplicateObj"
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

            theObject <- cleanDesignMatrixObj(theObject)

            updated_object <- theObject

            return(updated_object)
          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setGeneric(name="plotRleObj"
           , def=function( theObject, group, ylim ) {
             standardGeneric("plotRleObj")
           }
           , signature=c("theObject", "group", "ylim"))


#'@export
setMethod(f="plotRleObj"
          , definition=function( theObject, group, ylim) {
            protein_data <- theObject@protein_data
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id

            frozen_protein_matrix <- protein_data |>
              column_to_rownames(protein_id_column) |>
              as.matrix()

            design_matrix <- as.data.frame(design_matrix)
            rownames( design_matrix) <- design_matrix[,sample_id]

            # print( design_matrix)

            rowinfo_vector <- NA
            if( !is.na(group)){
              rowinfo_vector <-  design_matrix[colnames(frozen_protein_matrix), group]
            }

            # print( rowinfo_vector)

            rle_plot_before_cyclic_loess <- plotRle( t(frozen_protein_matrix)
                                                     , rowinfo = rowinfo_vector
                                                     , ylim=ylim)


            return( rle_plot_before_cyclic_loess)

          })



##----------------------------------------------------------------------------------------------------------------------------------------------------------------------



#'@export
setGeneric(name="plotPcaObj"
           , def=function( theObject, group_column, label_column, title, geom_text_size ) {
             standardGeneric("plotPcaObj")
           }
           , signature=c("theObject", "group_column", "label_column", "title", "geom_text_size"))


#'@export
setMethod(f="plotPcaObj"
          , definition=function( theObject, group_column, label_column, title, geom_text_size=8) {
            protein_data <- theObject@protein_data
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id


            frozen_protein_matrix <- protein_data |>
              column_to_rownames(protein_id_column) |>
              as.matrix()

            frozen_protein_matrix_pca <- frozen_protein_matrix
            frozen_protein_matrix_pca[!is.finite(frozen_protein_matrix_pca)] <- NA

            if( is.na(label_column) || label_column == "") {

              pca_plot_before_cyclic_loess_group <- plotPca( frozen_protein_matrix_pca
                                                             , design_matrix
                                                             , sample_id_column = !!sym( sample_id)
                                                             , group_column = !!sym( group_column)
                                                             , label_column = ""
                                                             , title = title
                                                             , geom.text.size = geom_text_size )            }
            else {

              pca_plot_before_cyclic_loess_group <- plotPca( frozen_protein_matrix_pca
                                                             , design_matrix
                                                             , sample_id_column = !!sym( sample_id)
                                                             , group_column = !!sym( group_column)
                                                             , label_column = !!sym( label_column)
                                                             , title = title
                                                             , geom.text.size = geom_text_size )

            }

            return( pca_plot_before_cyclic_loess_group)

          })




##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setGeneric(name="getPcaMatrixObj"
           , def=function( theObject) {
             standardGeneric("getPcaMatrixObj")
           }
           , signature=c("theObject"))


#'@export
setMethod(f="getPcaMatrixObj"
          , definition=function( theObject) {
            protein_data <- theObject@protein_data
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id


            frozen_protein_matrix <- protein_data |>
              column_to_rownames(protein_id_column) |>
              as.matrix()

            frozen_protein_matrix_pca <- frozen_protein_matrix
            frozen_protein_matrix_pca[!is.finite(frozen_protein_matrix_pca)] <- NA


            pca_mixomics_before_cyclic_loess <- mixOmics::pca(t(as.matrix(frozen_protein_matrix_pca)))$variates$X |>
              as.data.frame()    |>
              rownames_to_column(sample_id)  |>
              left_join(design_matrix, by = sample_id  )


            return( pca_mixomics_before_cyclic_loess)
          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------



# Calculate Pearson correlation between Tech rep 1 and 2
#'@export
setGeneric(name="proteinTechRepCorrelationObj"
           , def=function( theObject,  tech_rep_num_column, tech_rep_remove_regex) {
             standardGeneric("proteinTechRepCorrelationObj")
           }
           , signature=c("theObject", "tech_rep_num_column", "tech_rep_remove_regex"))

#'@export
setMethod( f = "proteinTechRepCorrelationObj"
          , definition=function( theObject,  tech_rep_num_column,  tech_rep_remove_regex ) {
            protein_data <- theObject@protein_data
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id
            tech_rep_column <- theObject@technical_replicate_id

            frozen_protein_matrix <- protein_data |>
              column_to_rownames(protein_id_column) |>
              as.matrix()

            frozen_protein_matrix_pca <- frozen_protein_matrix
            frozen_protein_matrix_pca[!is.finite(frozen_protein_matrix_pca)] <- NA

            protein_matrix_tech_rep <-proteinTechRepCorrelation( design_matrix, frozen_protein_matrix_pca
                                                                 , sample_id_column=sample_id
                                                                 , tech_rep_column = tech_rep_column
                                                                 , tech_rep_num_column = tech_rep_num_column
                                                                 , tech_rep_remove_regex = tech_rep_remove_regex )

            return( protein_matrix_tech_rep )
          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Normalize between Arrays

#'@export
setGeneric(name="normalizeBetweenArraysObj"
           , def=function( theObject, method) {
             standardGeneric("normalizeBetweenArraysObj")
           }
           , signature=c("theObject", "method"))


#'@export
setMethod(f="normalizeBetweenArraysObj"
          , definition=function( theObject,  method= "cyclicloess") {
            protein_data <- theObject@protein_data
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id



            frozen_protein_matrix <- protein_data |>
              column_to_rownames(protein_id_column) |>
              as.matrix()


            normalized_frozen_protein_matrix <- normalizeBetweenArrays( frozen_protein_matrix
                                                                        , method = method  )



            normalized_frozen_protein_matrix[!is.finite(normalized_frozen_protein_matrix)] <- NA

            normalized_frozen_protein_matrix_filt <- as.data.frame( normalized_frozen_protein_matrix ) |>
              #rownames_to_column("Protein.Ids") |>
              dplyr::filter( across( everything(), \(x) { !is.na(x) } ) ) |>
              #column_to_rownames("Protein.Ids") |>
              as.matrix()


            theObject@protein_data <- normalized_frozen_protein_matrix_filt |>
                      as.data.frame() |>
                      rownames_to_column(protein_id_column)

            theObject <- cleanDesignMatrixObj(theObject)

            updated_object <- theObject

            return(updated_object)

          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#'@export
setGeneric(name="pearsonCorForSamplePairsObj"
           , def=function( theObject,   tech_rep_remove_regex ) {
             standardGeneric("pearsonCorForSamplePairsObj")
           }
           , signature=c("theObject", "tech_rep_remove_regex"))

#'@export
setMethod(f="pearsonCorForSamplePairsObj"
          , definition=function( theObject, tech_rep_remove_regex ="pool" ) {
            protein_data <- theObject@protein_data
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id
            replicate_group_column <- theObject@technical_replicate_id


            frozen_mat_pca_long <- protein_data |>
              pivot_longer( cols=!matches(protein_id_column)
                            , values_to = "Protein.Normalized"
                            , names_to = sample_id) |>
              left_join( design_matrix
                         , by = join_by( !!sym(sample_id) == !!sym(sample_id))) |>
              mutate( temp = "")


            correlation_results_before_cyclic_loess <- calulatePearsonCorrelationForSamplePairs( design_matrix |>
                                                                                                   dplyr::select( !!sym(sample_id), !!sym(replicate_group_column) )
                                                                                                 , run_id_column = sample_id
                                                                                                 , replicate_group_column = replicate_group_column
                                                                                                 , frozen_mat_pca_long
                                                                                                 , num_of_cores = 1
                                                                                                 , sample_id_column = !!sym(sample_id)
                                                                                                 , protein_id_column = !!sym(protein_id_column)
                                                                                                 , peptide_sequence_column = temp
                                                                                                 , peptide_normalized_column = "Protein.Normalized")

            correlation_vec_before_cyclic_loess <- correlation_results_before_cyclic_loess |>
              dplyr::filter( !str_detect(!!sym(replicate_group_column), tech_rep_remove_regex )  )

           return( correlation_vec_before_cyclic_loess)
          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setGeneric(name="getNegCtrlProtAnovaObj"
           , def=function( theObject, ruv_group_id_column, num_neg_ctrl, q_val_thresh, fdr_method ) {
             standardGeneric("getNegCtrlProtAnovaObj")
           }
           , signature=c("theObject", "ruv_group_id_column", "num_neg_ctrl", "q_val_thresh", "fdr_method"))

#'@export
setMethod(f="getNegCtrlProtAnovaObj"
          , definition=function( theObject
                                 , ruv_group_id_column = "replicates"
                                 , num_neg_ctrl = 100
                                 , q_val_thresh = 0.05
                                 , fdr_method = "BH" ) {
            protein_data <- theObject@protein_data
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            group_id <- theObject@group_id
            sample_id <- theObject@sample_id
            replicate_group_column <- theObject@technical_replicate_id

            normalized_frozen_protein_matrix_filt <- protein_data |>
              column_to_rownames(protein_id_column) |>
              as.matrix()

            control_genes_index <- getNegCtrlProtAnova( normalized_frozen_protein_matrix_filt[,design_matrix |> pull(!!sym(sample_id)) ]
                                                        , design_matrix = design_matrix |>
                                                          column_to_rownames(sample_id) |>
                                                          dplyr::select( -!!sym(group_id))
                                                        , group_column = ruv_group_id_column
                                                        , num_neg_ctrl = num_neg_ctrl
                                                        , q_val_thresh = q_val_thresh
                                                        , fdr_method = fdr_method )

            return(control_genes_index)
          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
setGeneric(name="ruvCancorObj"
           , def=function( theObject, ctl, ncomp, ruv_group_id_column ) {
             standardGeneric("ruvCancorObj")
           }
           , signature=c("theObject", "ctl", "ncomp", "ruv_group_id_column"))

#'@export
setMethod( f = "ruvCancorObj"
           , definition=function( theObject, ctl, ncomp=2, ruv_group_id_column) {
             protein_data <- theObject@protein_data
             protein_id_column <- theObject@protein_id_column
             design_matrix <- theObject@design_matrix
             group_id <- theObject@group_id
             sample_id <- theObject@sample_id
             replicate_group_column <- theObject@technical_replicate_id

             if(! ruv_group_id_column %in% colnames(design_matrix)) {
               stop( paste0("The 'ruv_group_id_column = "
                            , ruv_group_id_column
                            , "' is not a column in the design matrix.") )
             }

             if( is.na(ncomp) || ncomp < 1) {
               stop(paste0("The ncomp = ", ncomp, " value is invalid."))
             }

             if( length( ctl) < 5 ) {
               stop(paste0( "The number of negative control molecules entered is less than 5. Please check the 'ctl' parameter."))
             }

             normalized_frozen_protein_matrix_filt <- protein_data |>
               column_to_rownames(protein_id_column) |>
               as.matrix()

             Y <- impute.nipals( t( normalized_frozen_protein_matrix_filt[,design_matrix |> pull(!!sym(sample_id))])
                                 , ncomp=ncomp)

             cancorplot_r2 <- ruv_cancorplot( Y ,
                                              X = design_matrix |>
                                                pull(!!sym(ruv_group_id_column)),
                                              ctl = ctl)
             cancorplot_r2


           })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#'@export
setGeneric(name="getRuvIIIReplicateMatrixObj"
           , def=function( theObject,  ruv_group_id_column ) {
             standardGeneric("getRuvIIIReplicateMatrixObj")
           }
           , signature=c("theObject", "ruv_group_id_column"))

#'@export
setMethod( f = "getRuvIIIReplicateMatrixObj"
           , definition=function( theObject, ruv_group_id_column) {
             protein_data <- theObject@protein_data
             protein_id_column <- theObject@protein_id_column
             design_matrix <- theObject@design_matrix
             group_id <- theObject@group_id
             sample_id <- theObject@sample_id
             replicate_group_column <- theObject@technical_replicate_id

             ruvIII_replicates_matrix <- getRuvIIIReplicateMatrix( design_matrix
                                                                   , !!sym(sample_id)
                                                                   , !!sym(ruv_group_id_column))
             return( ruvIII_replicates_matrix)
           })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#'@export
setGeneric(name="ruvIII_C_VaryingObj"
           , def=function( theObject, ruv_group_id_column, k, ctl) {
             standardGeneric("ruvIII_C_VaryingObj")
           }
           , signature=c("theObject", "ruv_group_id_column", "k", "ctl"))

#'@export
setMethod( f = "ruvIII_C_VaryingObj"
           , definition=function( theObject, ruv_group_id_column, k, ctl) {
             protein_data <- theObject@protein_data
             protein_id_column <- theObject@protein_id_column
             design_matrix <- theObject@design_matrix
             group_id <- theObject@group_id
             sample_id <- theObject@sample_id
             replicate_group_column <- theObject@technical_replicate_id

             normalized_frozen_protein_matrix_filt <- protein_data |>
               column_to_rownames(protein_id_column) |>
               as.matrix()

             Y <-  t( normalized_frozen_protein_matrix_filt[,design_matrix |> pull(!!sym(sample_id))])

             M <- getRuvIIIReplicateMatrix( design_matrix
                                            , !!sym(sample_id)
                                            , !!sym(ruv_group_id_column))

             cln_mat <- RUVIII_C_Varying( k = best_k
                                          , Y = Y
                                          , M = M
                                          , toCorrect = colnames(Y)
                                          , potentialControls = names( ctl[which(ctl)] ) )

             # Remove samples with no values
             cln_mat_2 <- cln_mat[rowSums(is.na(cln_mat) | is.nan(cln_mat)) != ncol(cln_mat),]

             # Remove proteins with no values
             cln_mat_3 <- t(cln_mat_2)
             cln_mat_4 <- cln_mat_3[rowSums(is.na(cln_mat_3) | is.nan(cln_mat_3)) != ncol(cln_mat_3),]

             ruv_normalized_results_cln <- cln_mat_4 |>
               as.data.frame() |>
               rownames_to_column(protein_id_column)

             theObject@protein_data <- ruv_normalized_results_cln

             theObject <- cleanDesignMatrixObj(theObject)

             return( theObject )
          })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
setGeneric(name="removeRowsWithMissingValuesPercentObj"
           , def=function( theObject
                           , ruv_group_id_column
                           , max_percent_miss_per_group = 60
                           , number_of_groups_missing = 2
                           , abundance_threshold = 1) {
             standardGeneric("removeRowsWithMissingValuesPercentObj")
           }
           , signature=c("theObject"
                         , "ruv_group_id_column"
                         , "max_percent_miss_per_group"
                         , "number_of_groups_missing"
                         , "abundance_threshold" ))

#'@export
setMethod( f = "removeRowsWithMissingValuesPercentObj"
           , definition=function( theObject
                                  , ruv_group_id_column
                                  , max_percent_miss_per_group = 50
                                  , number_of_groups_missing = 2
                                  , abundance_threshold = 1) {
             protein_data <- theObject@protein_data
             protein_id_column <- theObject@protein_id_column
             design_matrix <- theObject@design_matrix
             group_id <- theObject@group_id
             sample_id <- theObject@sample_id
             replicate_group_column <- theObject@technical_replicate_id

             theObject@protein_data <- removeRowsWithMissingValuesPercent( protein_data
                                                                           , !matches(protein_id_column)
                                                                           , design_matrix
                                                                           , !!sym(sample_id)
                                                                           , !!sym(protein_id_column)
                                                                           , !!sym(ruv_group_id_column)
                                                                           , max_percent_miss_per_group = max_percent_miss_per_group
                                                                           , number_of_groups_missing = number_of_groups_missing
                                                                           , abundance_threshold = abundance_threshold
                                                                           , temporary_abundance_column = "Log_Abundance")

             theObject <- cleanDesignMatrixObj(theObject)

             return(theObject)

           })





##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setGeneric(name="averageTechRepsObj"
           , def=function( theObject, design_matrix_columns ) {
             standardGeneric("averageTechRepsObj")
           }
           , signature=c("theObject", "design_matrix_columns" ))

#'@export
setMethod( f = "averageTechRepsObj"
           , definition=function( theObject, design_matrix_columns=c()  ) {

             protein_data <- theObject@protein_data
             protein_id_column <- theObject@protein_id_column
             design_matrix <- theObject@design_matrix
             group_id <- theObject@group_id
             sample_id <- theObject@sample_id
             replicate_group_column <- theObject@technical_replicate_id

             theObject@protein_data <- protein_data |>
               pivot_longer( cols = !matches( protein_id_column)
                             , names_to = sample_id
                             , values_to = "Log2.Protein.Imputed") |>
               left_join( design_matrix
                          , by = join_by( !!sym(sample_id) == !!sym(sample_id))) |>
               group_by( !!sym(protein_id_column), !!sym(replicate_group_column) )  |>
               summarise( Log2.Protein.Imputed = mean( Log2.Protein.Imputed, na.rm = TRUE)) |>
               ungroup() |>
               pivot_wider( names_from = !!sym(replicate_group_column)
                            , values_from = Log2.Protein.Imputed)

              theObject@sample_id <- theObject@technical_replicate_id

              theObject@design_matrix <- design_matrix |>
                dplyr::select(-!!sym( sample_id)) |>
                dplyr::select(one_of( unique( c( replicate_group_column,  group_id,  design_matrix_columns) ))) |>
                distinct()

              theObject <- cleanDesignMatrixObj(theObject)

              theObject

           })



##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
