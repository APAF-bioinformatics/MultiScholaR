

## Create S4 class for protomics protein level abundance data
#'@exportClass ProteinQuantitativeData
ProteinQuantitativeData <- setClass("ProteinQuantitativeData"

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
           protein_id_column = "Protein.Ids"

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
#'@export ProteinQuantitativeData

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------



##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
setGeneric( name ="setProteinData"
            , def=function( theect, protein_data, protein_id_column) {
                standardGeneric("setProteinData")
            })

#'@export
setMethod( f ="setProteinData"
           , signature = "ProteinQuantitativeData"
            , definition=function( theect, protein_data, protein_id_column ) {
              theect@protein_data <- protein_data
              theect@protein_id_column <- protein_id_column

              return(theect)
            })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Format the design matrix so that only metadata for samples in the protein data are retained, and also
# sort the sample IDs in the same order as the data matrix
#'@export
setGeneric(name ="cleanDesignMatrix"
           , def=function( theect) {
             standardGeneric("cleanDesignMatrix")
           })

#'@export
setMethod( f ="cleanDesignMatrix"
           , signature = "ProteinQuantitativeData"
           , definition=function( theect ) {

            samples_id_vector <- setdiff(colnames(theect@protein_data), theect@sample_id )

             theect@design_matrix <- data.frame( temp_sample_id = samples_id_vector )  |>
               inner_join( theect@design_matrix
                          , by = join_by ( temp_sample_id == !!sym(theect@sample_id)) ) |>
               dplyr::rename( !!sym(theect@sample_id) := "temp_sample_id" ) |>
               dplyr::filter( !!sym( theect@sample_id) %in% samples_id_vector )


             return(theect)
           })
##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
setGeneric(name="proteinIntensityFiltering"
           , def=function( theect, min_protein_intensity_percentile, proportion_samples_below_intensity_threshold, cluster) {
             standardGeneric("proteinIntensityFiltering")
           })

#'@export
setMethod( f="proteinIntensityFiltering"
           , signature="ProteinQuantitativeData"
           , definition = function( theect, min_protein_intensity_percentile, proportion_samples_below_intensity_threshold, cluster) {
             protein_data <- theect@protein_data

             data_long_cln <- protein_data  |>
               pivot_longer( cols=!matches(theect@protein_id_column)
                             , names_to = theect@sample_id
                             , values_to = "log_values")  |>
               mutate( temp = "")

             print(head( data_long_cln ))

             min_peptide_intensity_threshold <- ceiling( quantile( data_long_cln$log_values, na.rm=TRUE, probs = c(min_protein_intensity_percentile) ))[1]

             peptide_normalized_pif_cln <- peptideIntensityFilteringHelper( data_long_cln
                                                                      , min_peptide_intensity_threshold = min_peptide_intensity_threshold
                                                                      , proportion_samples_below_intensity_threshold = proportion_samples_below_intensity_threshold
                                                                      , protein_id_column = !!sym( theect@protein_id_column)
                                                                      , peptide_sequence_column = temp
                                                                      , peptide_quantity_column = log_values
                                                                      , cluster = cluster)


             theect@protein_data <- peptide_normalized_pif_cln |>
               dplyr::select( -temp) |>
               pivot_wider( id_cols = theect@protein_id_column , names_from = !!sym(theect@sample_id), values_from = log_values)

             theect <- cleanDesignMatrix(theect)

             updated_object <- theect

          return(updated_object)
          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
setGeneric(name="removeProteinsWithOnlyOneReplicate"
           , def=function( theect, cluster ) {
             standardGeneric("removeProteinsWithOnlyOneReplicate")
           }
           , signature=c("theect", "cluster"))

#'@export
setMethod(f="removeProteinsWithOnlyOneReplicate"
          , definition=function( theect, cluster) {
            protein_data <- theect@protein_data
            samples_id_tbl <- theect@design_matrix
            sample_id_tbl_sample_id_column <- theect@sample_id
            replicate_group_column <- theect@technical_replicate_id
            protein_id_column <- theect@protein_id_column

            input_table_sample_id_column <- 'Sample_ID'
            quantity_column <- "log_values"

            data_long_cln <- protein_data  |>
              pivot_longer( cols=!matches(protein_id_column)
                            , names_to = input_table_sample_id_column
                            , values_to = quantity_column)

            protein_data <- removeProteinsWithOnlyOneReplicateHelper( input_table = data_long_cln
                                                                , samples_id_tbl = samples_id_tbl
                                                                , cluster = cluster
                                                                , input_table_sample_id_column = !!sym( input_table_sample_id_column )
                                                                , sample_id_tbl_sample_id_column = !!sym( sample_id_tbl_sample_id_column)
                                                                , replicate_group_column = !!sym( replicate_group_column)
                                                                , protein_id_column = !!sym( protein_id_column)
                                                                , quantity_column = !!sym( quantity_column))


            theect@protein_data <- protein_data |>
              pivot_wider( id_cols = !!sym( protein_id_column)
                           , names_from = !!sym( input_table_sample_id_column)
                           , values_from = !!sym( quantity_column) )

            theect <- cleanDesignMatrix(theect)

            updated_object <- theect

            return(updated_object)
          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setGeneric(name="plotRle"
           , def=function( theect, group, ylim = c()) {
             standardGeneric("plotRle")
           }
           , signature=c("theect", "group", "ylim"))


#'@export
setMethod(f="plotRle"
          , signature="ProteinQuantitativeData"
          , definition=function( theect, group, ylim = c()) {
            protein_data <- theect@protein_data
            protein_id_column <- theect@protein_id_column
            design_matrix <- theect@design_matrix
            sample_id <- theect@sample_id

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
            if( length( ylim ) ==2 ) {

              rle_plot_before_cyclic_loess <- plotRleHelper( t(frozen_protein_matrix)
                                                       , rowinfo = rowinfo_vector
                                                       , ylim = ylim)
            } else {
              rle_plot_before_cyclic_loess <- plotRleHelper( t(frozen_protein_matrix)
                                                       , rowinfo = rowinfo_vector)
            }

            return( rle_plot_before_cyclic_loess)

          })



##----------------------------------------------------------------------------------------------------------------------------------------------------------------------



#'@export
setGeneric(name="plotPca"
           , def=function( theect, group_column, label_column, title, geom_text_size ) {
             standardGeneric("plotPca")
           }
           , signature=c("theect", "group_column", "label_column", "title", "geom_text_size"))


#'@export
setMethod(f="plotPca"
          , signature="ProteinQuantitativeData"
          , definition=function( theect, group_column, label_column, title, geom_text_size=8) {
            protein_data <- theect@protein_data
            protein_id_column <- theect@protein_id_column
            design_matrix <- theect@design_matrix
            sample_id <- theect@sample_id


            frozen_protein_matrix <- protein_data |>
              column_to_rownames(protein_id_column) |>
              as.matrix()

            frozen_protein_matrix_pca <- frozen_protein_matrix
            frozen_protein_matrix_pca[!is.finite(frozen_protein_matrix_pca)] <- NA

            if( is.na(label_column) || label_column == "") {

              pca_plot_before_cyclic_loess_group <- plotPcaHelper( frozen_protein_matrix_pca
                                                             , design_matrix
                                                             , sample_id_column = sample_id
                                                             , group_column =  group_column
                                                             , label_column = ""
                                                             , title = title
                                                             , geom.text.size = geom_text_size )            }
            else {

              pca_plot_before_cyclic_loess_group <- plotPcaHelper( frozen_protein_matrix_pca
                                                             , design_matrix
                                                             , sample_id_column =  sample_id
                                                             , group_column = group_column
                                                             , label_column =  label_column
                                                             , title = title
                                                             , geom.text.size = geom_text_size )

            }

            return( pca_plot_before_cyclic_loess_group)

          })




##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setGeneric(name="getPcaMatrix"
           , def=function( theect) {
             standardGeneric("getPcaMatrix")
           }
           , signature=c("theect"))


#'@export
setMethod(f="getPcaMatrix"
          , signature="ProteinQuantitativeData"
          , definition=function( theect) {
            protein_data <- theect@protein_data
            protein_id_column <- theect@protein_id_column
            design_matrix <- theect@design_matrix
            sample_id <- theect@sample_id


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
setGeneric(name="proteinTechRepCorrelation"
           , def=function( theect,  tech_rep_num_column, tech_rep_remove_regex) {
             standardGeneric("proteinTechRepCorrelation")
           }
           , signature=c("theect", "tech_rep_num_column", "tech_rep_remove_regex"))

#'@export
setMethod( f = "proteinTechRepCorrelation"
           , signature="ProteinQuantitativeData"
          , definition=function( theect,  tech_rep_num_column,  tech_rep_remove_regex ) {
            protein_data <- theect@protein_data
            protein_id_column <- theect@protein_id_column
            design_matrix <- theect@design_matrix
            sample_id <- theect@sample_id
            tech_rep_column <- theect@technical_replicate_id

            frozen_protein_matrix <- protein_data |>
              column_to_rownames(protein_id_column) |>
              as.matrix()

            frozen_protein_matrix_pca <- frozen_protein_matrix
            frozen_protein_matrix_pca[!is.finite(frozen_protein_matrix_pca)] <- NA

            protein_matrix_tech_rep <-proteinTechRepCorrelationHelper( design_matrix, frozen_protein_matrix_pca
                                                                 , sample_id_column=sample_id
                                                                 , tech_rep_column = tech_rep_column
                                                                 , tech_rep_num_column = tech_rep_num_column
                                                                 , tech_rep_remove_regex = tech_rep_remove_regex )

            return( protein_matrix_tech_rep )
          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Normalize between Arrays

#'@export
setGeneric(name="normalizeBetweenSamples"
           , def=function( theect, method) {
             standardGeneric("normalizeBetweenArrays")
           }
           , signature=c("theect", "method"))


#'@export
setMethod(f="normalizeBetweenSamples"
          , signature="ProteinQuantitativeData"
          , definition=function( theect,  method= "cyclicloess") {
            protein_data <- theect@protein_data
            protein_id_column <- theect@protein_id_column
            design_matrix <- theect@design_matrix
            sample_id <- theect@sample_id



            frozen_protein_matrix <- protein_data |>
              column_to_rownames(protein_id_column) |>
              as.matrix()

            frozen_protein_matrix[!is.finite(frozen_protein_matrix)] <- NA


            normalized_frozen_protein_matrix <- normalizeBetweenArrays( frozen_protein_matrix
                                                                        , method = method  )



            normalized_frozen_protein_matrix[!is.finite(normalized_frozen_protein_matrix)] <- NA

            normalized_frozen_protein_matrix_filt <- as.data.frame( normalized_frozen_protein_matrix ) |>
              #rownames_to_column("Protein.Ids") |>
              dplyr::filter( across( everything(), \(x) { !is.na(x) } ) ) |>
              #column_to_rownames("Protein.Ids") |>
              as.matrix()


            theect@protein_data <- normalized_frozen_protein_matrix_filt |>
                      as.data.frame() |>
                      rownames_to_column(protein_id_column)

            theect <- cleanDesignMatrix(theect)

            updated_object <- theect

            return(updated_object)

          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#'@export
setGeneric(name="pearsonCorForSamplePairs"
           , def=function( theect,   tech_rep_remove_regex ) {
             standardGeneric("pearsonCorForSamplePairs")
           }
           , signature=c("theect", "tech_rep_remove_regex"))

#'@export
setMethod(f="pearsonCorForSamplePairs"
          , signature="ProteinQuantitativeData"
          , definition=function( theect, tech_rep_remove_regex ="pool" ) {
            protein_data <- theect@protein_data
            protein_id_column <- theect@protein_id_column
            design_matrix <- theect@design_matrix
            sample_id <- theect@sample_id
            replicate_group_column <- theect@technical_replicate_id


            frozen_mat_pca_long <- protein_data |>
              pivot_longer( cols=!matches(protein_id_column)
                            , values_to = "Protein.Normalized"
                            , names_to = sample_id) |>
              left_join( design_matrix
                         , by = join_by( !!sym(sample_id) == !!sym(sample_id))) |>
              mutate( temp = "")


            correlation_results_before_cyclic_loess <- calulatePearsonCorrelationForSamplePairsHelper( design_matrix |>
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
setGeneric(name="getNegCtrlProtAnova"
           , def=function( theect, ruv_group_id_column, num_neg_ctrl, q_val_thresh, fdr_method ) {
             standardGeneric("getNegCtrlProtAnova")
           }
           , signature=c("theect", "ruv_group_id_column", "num_neg_ctrl", "q_val_thresh", "fdr_method"))

#'@export
setMethod(f="getNegCtrlProtAnova"
          , signature="ProteinQuantitativeData"
          , definition=function( theect
                                 , ruv_group_id_column = "replicates"
                                 , num_neg_ctrl = 100
                                 , q_val_thresh = 0.05
                                 , fdr_method = "BH" ) {
            protein_data <- theect@protein_data
            protein_id_column <- theect@protein_id_column
            design_matrix <- theect@design_matrix
            group_id <- theect@group_id
            sample_id <- theect@sample_id
            replicate_group_column <- theect@technical_replicate_id

            normalized_frozen_protein_matrix_filt <- protein_data |>
              column_to_rownames(protein_id_column) |>
              as.matrix()

            control_genes_index <- getNegCtrlProtAnovaHelper( normalized_frozen_protein_matrix_filt[,design_matrix |> pull(!!sym(sample_id)) ]
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

#'@description Sort proteins by their coefficient of variation and take the top N with lowest coefficient of variation
#'@export
setGeneric(name="getLowCoefficientOfVariationProteins"
           , def=function( theect, num_neg_ctrl ) {
             standardGeneric("getLowCoefficientOfVariationProteins")
           }
           , signature=c("theect", "num_neg_ctrl"))



#'@export
setMethod( f = "getLowCoefficientOfVariationProteins"
           , signature="ProteinQuantitativeData"
           , definition=function( theect
                                                  , num_neg_ctrl = 100) {

  list_of_control_genes <- theect@protein_data |>
    column_to_rownames(theect@protein_id_column) |>
    t() |>
    as.data.frame() |>
    summarise( across(everything(), ~sd(.)/mean(.))) |>
    t() |>
    as.data.frame() |>
    dplyr::rename( coefficient_of_variation = "V1") |>
    tibble::rownames_to_column(theect@protein_id_column) |>
    arrange( coefficient_of_variation) |>
    head(num_neg_ctrl)

  control_gene_index_helper <- theect@protein_data |>
    dplyr::select(theect@protein_id_column) |>
    mutate( index = row_number()) |>
    left_join( list_of_control_genes, by = theect@protein_id_column)  |>
    mutate( is_selected = case_when( is.na(coefficient_of_variation) ~ FALSE
                                     , TRUE ~  TRUE ) ) |>
    arrange( index) |>
    dplyr::select( Protein.Ids, is_selected) |>
    column_to_rownames(theect@protein_id_column) |>
    t()

  control_gene_index <- control_gene_index_helper[1,]

  control_gene_index

})

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
setGeneric(name="ruvCancor"
           , def=function( theect, ctl, ncomp, ruv_group_id_column ) {
             standardGeneric("ruvCancor")
           }
           , signature=c("theect", "ctl", "ncomp", "ruv_group_id_column"))

#'@export
setMethod( f = "ruvCancor"
           , signature="ProteinQuantitativeData"
           , definition=function( theect, ctl, ncomp=2, ruv_group_id_column) {
             protein_data <- theect@protein_data
             protein_id_column <- theect@protein_id_column
             design_matrix <- theect@design_matrix
             group_id <- theect@group_id
             sample_id <- theect@sample_id
             replicate_group_column <- theect@technical_replicate_id

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

             Y <-  t( normalized_frozen_protein_matrix_filt[,design_matrix |> pull(!!sym(sample_id))])
             if( length(which( is.na(normalized_frozen_protein_matrix_filt) )) > 0 ) {
               Y <- impute.nipals( t( normalized_frozen_protein_matrix_filt[,design_matrix |> pull(!!sym(sample_id))])
                                   , ncomp=ncomp)
             }

             cancorplot_r2 <- ruv_cancorplot( Y ,
                                              X = design_matrix |>
                                                pull(!!sym(ruv_group_id_column)),
                                              ctl = ctl)
             cancorplot_r2


           })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#'@export
setGeneric(name="getRuvIIIReplicateMatrix"
           , def=function( theect,  ruv_group_id_column ) {
             standardGeneric("getRuvIIIReplicateMatrix")
           }
           , signature=c("theect", "ruv_group_id_column"))

#'@export
setMethod( f = "getRuvIIIReplicateMatrix"
           , signature="ProteinQuantitativeData"
           , definition=function( theect, ruv_group_id_column) {
             protein_data <- theect@protein_data
             protein_id_column <- theect@protein_id_column
             design_matrix <- theect@design_matrix
             group_id <- theect@group_id
             sample_id <- theect@sample_id
             replicate_group_column <- theect@technical_replicate_id

             ruvIII_replicates_matrix <- getRuvIIIReplicateMatrixHelper( design_matrix
                                                                   , !!sym(sample_id)
                                                                   , !!sym(ruv_group_id_column))
             return( ruvIII_replicates_matrix)
           })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#'@export
setGeneric(name="ruvIII_C_Varying"
           , def=function( theect, ruv_group_id_column, k, ctl) {
             standardGeneric("ruvIII_C_Varying")
           }
           , signature=c("theect", "ruv_group_id_column", "k", "ctl"))

#'@export
setMethod( f = "ruvIII_C_Varying"
           , signature="ProteinQuantitativeData"
           , definition=function( theect, ruv_group_id_column, k, ctl) {
             protein_data <- theect@protein_data
             protein_id_column <- theect@protein_id_column
             design_matrix <- theect@design_matrix
             group_id <- theect@group_id
             sample_id <- theect@sample_id
             replicate_group_column <- theect@technical_replicate_id

             normalized_frozen_protein_matrix_filt <- protein_data |>
               column_to_rownames(protein_id_column) |>
               as.matrix()

             Y <-  t( normalized_frozen_protein_matrix_filt[,design_matrix |> pull(!!sym(sample_id))])

             M <- getRuvIIIReplicateMatrixHelper( design_matrix
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

             theect@protein_data <- ruv_normalized_results_cln

             theect <- cleanDesignMatrix(theect)

             return( theect )
          })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
setGeneric(name="removeRowsWithMissingValuesPercent"
           , def=function( theect
                           , ruv_group_id_column
                           , max_perc_below_thresh_per_group = 50
                           , max_perc_of_groups_below_thresh = 50
                           , min_protein_intensity_percentile = 1) {
             standardGeneric("removeRowsWithMissingValuesPercent")
           }
           , signature=c("theect"
                         , "ruv_group_id_column"
                         , "max_perc_below_thresh_per_group"
                         , "max_perc_of_groups_below_thresh"
                         , "min_protein_intensity_percentile" ))

#'@export
setMethod( f = "removeRowsWithMissingValuesPercent"
           , signature="ProteinQuantitativeData"
           , definition=function( theect
                                  , ruv_group_id_column
                                  , max_perc_below_thresh_per_group = 50
                                  , max_perc_of_groups_below_thresh = 50
                                  , min_protein_intensity_percentile = 1) {

             protein_data <- theect@protein_data
             protein_id_column <- theect@protein_id_column
             design_matrix <- theect@design_matrix
             group_id <- theect@group_id
             sample_id <- theect@sample_id
             replicate_group_column <- theect@technical_replicate_id


             # print(min_protein_intensity_threshold )

             theect@protein_data <- removeRowsWithMissingValuesPercentHelper( protein_data
                                                                           , cols= !matches(protein_id_column)
                                                                           , design_matrix = design_matrix
                                                                           , sample_id = !!sym(sample_id)
                                                                           , row_id = !!sym(protein_id_column)
                                                                           , group_column = !!sym(ruv_group_id_column)
                                                                           , max_perc_below_thresh_per_group = max_perc_below_thresh_per_group
                                                                           , max_perc_of_groups_below_thresh = max_perc_of_groups_below_thresh
                                                                           , min_protein_intensity_percentile = min_protein_intensity_percentile
                                                                           , temporary_abundance_column = "Log_Abundance")

             theect <- cleanDesignMatrix(theect)

             return(theect)

           })





##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setGeneric(name="averageTechReps"
           , def=function( theect, design_matrix_columns ) {
             standardGeneric("averageTechReps")
           }
           , signature=c("theect", "design_matrix_columns" ))

#'@export
setMethod( f = "averageTechReps"
           , signature="ProteinQuantitativeData"
           , definition=function( theect, design_matrix_columns=c()  ) {

             protein_data <- theect@protein_data
             protein_id_column <- theect@protein_id_column
             design_matrix <- theect@design_matrix
             group_id <- theect@group_id
             sample_id <- theect@sample_id
             replicate_group_column <- theect@technical_replicate_id

             theect@protein_data <- protein_data |>
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

              theect@sample_id <- theect@technical_replicate_id

              theect@design_matrix <- design_matrix |>
                dplyr::select(-!!sym( sample_id)) |>
                dplyr::select(one_of( unique( c( replicate_group_column,  group_id,  design_matrix_columns) ))) |>
                distinct()

              rownames(theect@design_matrix ) <- theect@design_matrix |> pull(one_of(replicate_group_column))
              theect@sample_id <- replicate_group_column
              theect@technical_replicate_id <- NA_character_

              theect <- cleanDesignMatrix(theect)

              theect

           })



##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setGeneric(name="chooseBestProteinAccession"
           , def=function( theect, delim, seqinr_obj, seqinr_accession_column ) {
             standardGeneric("chooseBestProteinAccession")
           }
           , signature=c("theect", "delim", "seqinr_obj", "seqinr_accession_column" ))

#'@export
setMethod( f = "chooseBestProteinAccession"
           , signature="ProteinQuantitativeData"
           , definition=function( theect, delim=";", seqinr_obj, seqinr_accession_column  ) {
             protein_data <- theect@protein_data
             protein_id_column <- theect@protein_id_column


             evidence_tbl_cleaned <- protein_data |>
               distinct() |>
               mutate( row_id = row_number() -1 )

             accession_gene_name_tbl <- chooseBestProteinAccessionHelper(input_tbl = evidence_tbl_cleaned,
                                                                   acc_detail_tab = seqinr_obj,
                                                                   accessions_column = !!sym( protein_id_column),
                                                                   row_id_column = !!sym( seqinr_accession_column),
                                                                   group_id = row_id,
                                                                   delim = ";")


             protein_log2_quant_cln <- evidence_tbl_cleaned |>
               left_join( accession_gene_name_tbl |>
                            dplyr::distinct( row_id, !!sym( seqinr_accession_column) )
                          , by = join_by(row_id) ) |>
               mutate( !!sym( theect@protein_id_column ) :=!!sym( seqinr_accession_column) )  |>
               dplyr::select(-row_id, -!!sym( seqinr_accession_column))

             theect@protein_data <- protein_log2_quant_cln

             return(theect)

           })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#'@export
setGeneric(name="chooseBestProteinAccessionSumDuplicates"
           , def=function( theect, delim, quant_columns_pattern, islogged ) {
             standardGeneric("chooseBestProteinAccessionSumDuplicates")
           }
           , signature=c("theect", "delim", "quant_columns_pattern", "islogged" ))

#'@export
setMethod( f = "chooseBestProteinAccessionSumDuplicates"
           , signature="ProteinQuantitativeData"
           , definition=function( theect, delim=";", quant_columns_pattern = "\\d+", islogged = TRUE ) {

             protein_data <- theect@protein_data
             protein_id_column <- theect@protein_id_column


             protein_log2_quant_cln <- protein_data |>
               mutate( !!sym(protein_id_column) := str_split_i(!!sym( protein_id_column), delim, 1 ) ) |>
               group_by( !!sym(protein_id_column) ) |>
               summarise ( across( matches(quant_columns_pattern), \(x){ if(islogged==TRUE) {log2(sum(2^x, na.rm = TRUE)) } else { sum(x, na.rm = TRUE)}    })) |>
               ungroup()

             theect@protein_data <- protein_log2_quant_cln

             theect

           })
