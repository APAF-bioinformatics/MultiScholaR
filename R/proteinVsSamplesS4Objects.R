

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
                      , technical_replicate_id="character" )

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
            , def=function( theObject, protein_data, protein_id_column) {
                standardGeneric("setProteinData")
            })

#'@export
setMethod( f ="setProteinData"
           , signature = "ProteinQuantitativeData"
            , definition=function( theObject, protein_data, protein_id_column ) {
              theObject@protein_data <- protein_data
              theObject@protein_id_column <- protein_id_column

              return(theObject)
            })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Format the design matrix so that only metadata for samples in the protein data are retained, and also
# sort the sample IDs in the same order as the data matrix
#'@export
setGeneric(name ="cleanDesignMatrix"
           , def=function( theObject) {
             standardGeneric("cleanDesignMatrix")
           })

#'@export
setMethod( f ="cleanDesignMatrix"
           , signature = "ProteinQuantitativeData"
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
setGeneric(name="proteinIntensityFiltering"
           , def=function( theObject, min_protein_intensity_percentile, proportion_samples_below_intensity_threshold, cluster) {
             standardGeneric("proteinIntensityFiltering")
           })

#'@export
setMethod( f="proteinIntensityFiltering"
           , signature="ProteinQuantitativeData"
           , definition = function( theObject, min_protein_intensity_percentile, proportion_samples_below_intensity_threshold, cluster) {
             protein_data <- theObject@protein_data

             data_long_cln <- protein_data  |>
               pivot_longer( cols=!matches(theObject@protein_id_column)
                             , names_to = theObject@sample_id
                             , values_to = "log_values")  |>
               mutate( temp = "")

             print(head( data_long_cln ))

             min_peptide_intensity_threshold <- ceiling( quantile( data_long_cln$log_values, na.rm=TRUE, probs = c(min_protein_intensity_percentile) ))[1]

             peptide_normalized_pif_cln <- peptideIntensityFilteringHelper( data_long_cln
                                                                      , min_peptide_intensity_threshold = min_peptide_intensity_threshold
                                                                      , proportion_samples_below_intensity_threshold = proportion_samples_below_intensity_threshold
                                                                      , protein_id_column = !!sym( theObject@protein_id_column)
                                                                      , peptide_sequence_column = temp
                                                                      , peptide_quantity_column = log_values
                                                                      , cluster = cluster)


             theObject@protein_data <- peptide_normalized_pif_cln |>
               dplyr::select( -temp) |>
               pivot_wider( id_cols = theObject@protein_id_column , names_from = !!sym(theObject@sample_id), values_from = log_values)

             theObject <- cleanDesignMatrix(theObject)

             updated_object <- theObject

          return(updated_object)
          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
setGeneric(name="removeProteinsWithOnlyOneReplicate"
           , def=function( theObject, cluster ) {
             standardGeneric("removeProteinsWithOnlyOneReplicate")
           }
           , signature=c("theObject", "cluster"))

#'@export
setMethod(f="removeProteinsWithOnlyOneReplicate"
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

            protein_data <- removeProteinsWithOnlyOneReplicateHelper( input_table = data_long_cln
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

            theObject <- cleanDesignMatrix(theObject)

            updated_object <- theObject

            return(updated_object)
          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setGeneric(name="plotRle"
           , def=function( theObject, group, ylim = c()) {
             standardGeneric("plotRle")
           }
           , signature=c("theObject", "group", "ylim"))


#'@export
setMethod(f="plotRle"
          , signature="ProteinQuantitativeData"
          , definition=function( theObject, group, ylim = c()) {
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

              rle_plot_before_cyclic_loess <- plotRleHelper( t(frozen_protein_matrix)
                                                       , rowinfo = rowinfo_vector
                                                       , ylim = ylim)

            return( rle_plot_before_cyclic_loess)

          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setGeneric(name="plotRleList"
           , def=function( theObject, list_of_columns, ylim = c()) {
             standardGeneric("plotRleList")
           }
           , signature=c("theObject", "list_of_columns", "ylim"))

#'@export
setMethod(f="plotRleList"
          , signature="ProteinQuantitativeData"
          , definition=function( theObject, list_of_columns, ylim = c()) {
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

            runOneRle <- function( column_name) {
              rowinfo_vector <- NA

              if ( column_name %in% colnames(design_matrix) ) {
                rowinfo_vector <- design_matrix[colnames(frozen_protein_matrix), column_name]
              }

              rle_plot_before_cyclic_loess <- plotRleHelper( t(frozen_protein_matrix)
                                                             , rowinfo = rowinfo_vector
                                                             , ylim = ylim)

              return( rle_plot_before_cyclic_loess)
            }

            list_of_rle_plots <- purrr::map( list_of_columns, runOneRle)

            names(list_of_rle_plots) <- list_of_columns

            return( list_of_rle_plots)

          })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @export
savePlotRleList <- function( input_list, prefix = "RLE", suffix = c("png", "pdf"), output_dir ) {

  list_of_filenames <- expand_grid( column=names(input_list), suffix=suffix)  |>
    mutate( filename= paste0( "RLE", "_", column , ".", suffix))  |>
    left_join( tibble( column =names( input_list)
                       ,  plots= input_list )
               , by=join_by(column ))


  purrr::walk2( list_of_filenames$plots
                , list_of_filenames$filename
                , \(.x, .y){
                  ggsave( plot=.x, filename= file.path(output_dir, .y))
                } )

  list_of_filenames

}



##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setGeneric(name="plotPca"
           , def=function( theObject, group_column, label_column, title, geom_text_size ) {
             standardGeneric("plotPca")
           }
           , signature=c("theObject", "group_column", "label_column", "title", "geom_text_size"))

#'@export
setMethod(f="plotPca"
          , signature="ProteinQuantitativeData"
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
              label_column <- ""
            }

            pca_plot <- plotPcaHelper( frozen_protein_matrix_pca
                                       , design_matrix
                                       , sample_id_column =  sample_id
                                       , group_column = group_column
                                       , label_column =  label_column
                                       , title = title
                                       , geom.text.size = geom_text_size )

            return( pca_plot)
          })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setGeneric(name="plotPcaList"
           , def=function( theObject, group_columns_list, label_column, title, geom_text_size ) {
             standardGeneric("plotPcaList")
           }
           , signature=c("theObject", "group_columns_list", "label_column", "title", "geom_text_size"))

#'@export
setMethod(f="plotPcaList"
          , signature="ProteinQuantitativeData"
          , definition=function( theObject, group_columns_list, label_column, title, geom_text_size=8) {
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
              label_column <- ""
            }

            pca_plots_list <- plotPcaListHelper( frozen_protein_matrix_pca
                                                 , design_matrix
                                                 , sample_id_column =  sample_id
                                                 , group_columns_list = group_columns_list
                                                 , label_column =  label_column
                                                 , title = title
                                                 , geom.text.size = geom_text_size )

            return( pca_plots_list)
          })


#' @export
savePlotPcaList <- function( input_list, prefix = "PCA", suffix = c("png", "pdf"), output_dir ) {

  list_of_filenames <- expand_grid( column=names(input_list), suffix=suffix)  |>
    mutate( filename= paste0( "RLE", "_", column , ".", suffix))  |>
    left_join( tibble( column =names( input_list)
                       ,  plots= input_list )
               , by=join_by(column ))


  purrr::walk2( list_of_filenames$plots
                , list_of_filenames$filename
                , \(.x, .y){
                  ggsave( plot=.x, filename= file.path(output_dir, .y))
                } )

  list_of_filenames

}

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setGeneric(name="getPcaMatrix"
           , def=function( theObject) {
             standardGeneric("getPcaMatrix")
           }
           , signature=c("theObject"))


#'@export
setMethod(f="getPcaMatrix"
          , signature="ProteinQuantitativeData"
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
setGeneric(name="proteinTechRepCorrelation"
           , def=function( theObject,  tech_rep_num_column, tech_rep_remove_regex) {
             standardGeneric("proteinTechRepCorrelation")
           }
           , signature=c("theObject", "tech_rep_num_column", "tech_rep_remove_regex"))

#'@export
setMethod( f = "proteinTechRepCorrelation"
           , signature="ProteinQuantitativeData"
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

            protein_matrix_tech_rep <-proteinTechRepCorrelationHelper( design_matrix, frozen_protein_matrix_pca
                                                                       , protein_id_column = protein_id_column
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
           , def=function( theObject, method) {
             standardGeneric("normalizeBetweenSamples")
           }
           , signature=c("theObject", "method"))


#'@export
setMethod(f="normalizeBetweenSamples"
          , signature="ProteinQuantitativeData"
          , definition=function( theObject,  method= "cyclicloess") {
            protein_data <- theObject@protein_data
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id



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


            theObject@protein_data <- normalized_frozen_protein_matrix_filt |>
                      as.data.frame() |>
                      rownames_to_column(protein_id_column)

            theObject <- cleanDesignMatrix(theObject)

            updated_object <- theObject

            return(updated_object)

          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#'@export
setGeneric(name="pearsonCorForSamplePairs"
           , def=function( theObject,   tech_rep_remove_regex ) {
             standardGeneric("pearsonCorForSamplePairs")
           }
           , signature=c("theObject", "tech_rep_remove_regex"))

#'@export
setMethod(f="pearsonCorForSamplePairs"
          , signature="ProteinQuantitativeData"
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
           , def=function( theObject, ruv_group_id_column, num_neg_ctrl, q_val_thresh, fdr_method ) {
             standardGeneric("getNegCtrlProtAnova")
           }
           , signature=c("theObject", "ruv_group_id_column", "num_neg_ctrl", "q_val_thresh", "fdr_method"))

#'@export
setMethod(f="getNegCtrlProtAnova"
          , signature="ProteinQuantitativeData"
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
           , def=function( theObject, num_neg_ctrl ) {
             standardGeneric("getLowCoefficientOfVariationProteins")
           }
           , signature=c("theObject", "num_neg_ctrl"))



#'@export
setMethod( f = "getLowCoefficientOfVariationProteins"
           , signature="ProteinQuantitativeData"
           , definition=function( theObject
                                                  , num_neg_ctrl = 100) {

  list_of_control_genes <- theObject@protein_data |>
    column_to_rownames(theObject@protein_id_column) |>
    t() |>
    as.data.frame() |>
    summarise( across(everything(), ~sd(.)/mean(.))) |>
    t() |>
    as.data.frame() |>
    dplyr::rename( coefficient_of_variation = "V1") |>
    tibble::rownames_to_column(theObject@protein_id_column) |>
    arrange( coefficient_of_variation) |>
    head(num_neg_ctrl)

  control_gene_index_helper <- theObject@protein_data |>
    dplyr::select(theObject@protein_id_column) |>
    mutate( index = row_number()) |>
    left_join( list_of_control_genes, by = theObject@protein_id_column)  |>
    mutate( is_selected = case_when( is.na(coefficient_of_variation) ~ FALSE
                                     , TRUE ~  TRUE ) ) |>
    arrange( index) |>
    dplyr::select( !!sym(theObject@protein_id_column), is_selected) |>
    column_to_rownames(theObject@protein_id_column) |>
    t()

  control_gene_index <- control_gene_index_helper[1,]

  control_gene_index

})

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
setGeneric(name="ruvCancor"
           , def=function( theObject, ctl, ncomp, ruv_group_id_column ) {
             standardGeneric("ruvCancor")
           }
           , signature=c("theObject", "ctl", "ncomp", "ruv_group_id_column"))

#'@export
setMethod( f = "ruvCancor"
           , signature="ProteinQuantitativeData"
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
           , def=function( theObject,  ruv_group_id_column ) {
             standardGeneric("getRuvIIIReplicateMatrix")
           }
           , signature=c("theObject", "ruv_group_id_column"))

#'@export
setMethod( f = "getRuvIIIReplicateMatrix"
           , signature="ProteinQuantitativeData"
           , definition=function( theObject, ruv_group_id_column) {
             protein_data <- theObject@protein_data
             protein_id_column <- theObject@protein_id_column
             design_matrix <- theObject@design_matrix
             group_id <- theObject@group_id
             sample_id <- theObject@sample_id
             replicate_group_column <- theObject@technical_replicate_id

             ruvIII_replicates_matrix <- getRuvIIIReplicateMatrixHelper( design_matrix
                                                                   , !!sym(sample_id)
                                                                   , !!sym(ruv_group_id_column))
             return( ruvIII_replicates_matrix)
           })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#'@export
setGeneric(name="ruvIII_C_Varying"
           , def=function( theObject, ruv_group_id_column, k, ctl) {
             standardGeneric("ruvIII_C_Varying")
           }
           , signature=c("theObject", "ruv_group_id_column", "k", "ctl"))

#'@export
setMethod( f = "ruvIII_C_Varying"
           , signature="ProteinQuantitativeData"
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

             theObject@protein_data <- ruv_normalized_results_cln

             theObject <- cleanDesignMatrix(theObject)

             return( theObject )
          })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
setGeneric(name="removeRowsWithMissingValuesPercent"
           , def=function( theObject
                           , ruv_group_id_column
                           , max_perc_below_thresh_per_group = 50
                           , max_perc_of_groups_below_thresh = 50
                           , min_protein_intensity_percentile = 1) {
             standardGeneric("removeRowsWithMissingValuesPercent")
           }
           , signature=c("theObject"
                         , "ruv_group_id_column"
                         , "max_perc_below_thresh_per_group"
                         , "max_perc_of_groups_below_thresh"
                         , "min_protein_intensity_percentile" ))

#'@export
setMethod( f = "removeRowsWithMissingValuesPercent"
           , signature="ProteinQuantitativeData"
           , definition=function( theObject
                                  , ruv_group_id_column
                                  , max_perc_below_thresh_per_group = 50
                                  , max_perc_of_groups_below_thresh = 50
                                  , min_protein_intensity_percentile = 1) {

             protein_data <- theObject@protein_data
             protein_id_column <- theObject@protein_id_column
             design_matrix <- theObject@design_matrix
             group_id <- theObject@group_id
             sample_id <- theObject@sample_id
             replicate_group_column <- theObject@technical_replicate_id


             # print(min_protein_intensity_threshold )

             theObject@protein_data <- removeRowsWithMissingValuesPercentHelper( protein_data
                                                                           , cols= !matches(protein_id_column)
                                                                           , design_matrix = design_matrix
                                                                           , sample_id = !!sym(sample_id)
                                                                           , row_id = !!sym(protein_id_column)
                                                                           , group_column = !!sym(ruv_group_id_column)
                                                                           , max_perc_below_thresh_per_group = max_perc_below_thresh_per_group
                                                                           , max_perc_of_groups_below_thresh = max_perc_of_groups_below_thresh
                                                                           , min_protein_intensity_percentile = min_protein_intensity_percentile
                                                                           , temporary_abundance_column = "Log_Abundance")

             theObject <- cleanDesignMatrix(theObject)

             return(theObject)

           })





##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setGeneric(name="averageTechReps"
           , def=function( theObject, design_matrix_columns ) {
             standardGeneric("averageTechReps")
           }
           , signature=c("theObject", "design_matrix_columns" ))

#'@export
setMethod( f = "averageTechReps"
           , signature="ProteinQuantitativeData"
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
                dplyr::select(all_of( unique( c( replicate_group_column,  group_id,  design_matrix_columns) ))) |>
                distinct()

              theObject@sample_id <- replicate_group_column
              theObject@technical_replicate_id <- NA_character_

              theObject <- cleanDesignMatrix(theObject)

              theObject

           })



##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setGeneric(name="chooseBestProteinAccession"
           , def=function( theObject, delim, seqinr_obj, seqinr_accession_column ) {
             standardGeneric("chooseBestProteinAccession")
           }
           , signature=c("theObject", "delim", "seqinr_obj", "seqinr_accession_column" ))

#'@export
setMethod( f = "chooseBestProteinAccession"
           , signature="ProteinQuantitativeData"
           , definition=function( theObject, delim=";", seqinr_obj, seqinr_accession_column  ) {
             protein_data <- theObject@protein_data
             protein_id_column <- theObject@protein_id_column


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
                          , by = join_by( row_id,
                                          !!sym( theObject@protein_id_column ) ==!!sym( seqinr_accession_column)) )

             protein_log2_quant_cln <- protein_log2_quant_cln |>
               dplyr::select(-row_id)

             theObject@protein_data <- protein_log2_quant_cln

             return(theObject)

           })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#'@export
setGeneric(name="chooseBestProteinAccessionSumDuplicates"
           , def=function( theObject, delim, quant_columns_pattern, islogged ) {
             standardGeneric("chooseBestProteinAccessionSumDuplicates")
           }
           , signature=c("theObject", "delim", "quant_columns_pattern", "islogged" ))

#'@export
setMethod( f = "chooseBestProteinAccessionSumDuplicates"
           , signature="ProteinQuantitativeData"
           , definition=function( theObject, delim=";", quant_columns_pattern = "\\d+", islogged = TRUE ) {

             protein_data <- theObject@protein_data
             protein_id_column <- theObject@protein_id_column


             protein_log2_quant_cln <- protein_data |>
               mutate( !!sym(protein_id_column) := str_split_i(!!sym( protein_id_column), delim, 1 ) ) |>
               group_by( !!sym(protein_id_column) ) |>
               summarise ( across( matches(quant_columns_pattern), \(x){ if(islogged==TRUE) {log2(sum(2^x, na.rm = TRUE)) } else { sum(x, na.rm = TRUE)}    })) |>
               ungroup()

             theObject@protein_data <- protein_log2_quant_cln

             theObject

           })



##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setGeneric(name="filterSamplesByProteinCorrelationThreshold"
           , def=function( theObject, pearson_correlation_per_pair, min_pearson_correlation_threshold ) {
             standardGeneric("filterSamplesByProteinCorrelationThreshold")
           }
           , signature=c("theObject", "pearson_correlation_per_pair", "min_pearson_correlation_threshold" ))

#'@export
setMethod( f = "filterSamplesByProteinCorrelationThreshold"
           , signature="ProteinQuantitativeData"
           , definition=function( theObject, pearson_correlation_per_pair, min_pearson_correlation_threshold = 0.75 ) {


             filtered_table <- filterSamplesByProteinCorrelationThresholdHelper (
               pearson_correlation_per_pair
               , protein_intensity_table = theObject@protein_data
               , min_pearson_correlation_threshold = min_pearson_correlation_threshold
               , filename_column_x = !!sym( paste0( theObject@sample_id, ".x") )
               , filename_column_y = !!sym( paste0( theObject@sample_id, ".y") )
               , protein_id_column = theObject@protein_id_column
               , correlation_column = pearson_correlation )

             theObject@protein_data <- filtered_table

             theObject <- cleanDesignMatrix(theObject)

             theObject
             })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
