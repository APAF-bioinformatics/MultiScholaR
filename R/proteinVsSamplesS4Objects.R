

## Create S4 class for protomics protein level abundance data
#'@exportClass ProteinQuantitativeData
ProteinQuantitativeData <- setClass("ProteinQuantitativeData"
         , slots = c(
                      # Protein vs Sample quantitative data
                      protein_quant_table = "data.frame"
                      , protein_id_column = "character"

                      # Design Matrix Information
                      , design_matrix = "data.frame"
                      , protein_id_table = "data.frame"
                      , sample_id="character"
                      , group_id="character"
                      , technical_replicate_id="character"
                      , args = "list")

         , prototype = list(
           # Protein vs Sample quantitative data
           protein_id_column = "Protein.Ids"
          , protein_id_table = data.frame()
           # Design Matrix Information
           , sample_id="Sample_id"
           , group_id="group"
           , technical_replicate_id="replicates"
           , args = NULL
           )

         , validity = function(object) {
           if( !is.data.frame(object@protein_quant_table) ) {
             stop("protein_quant_table must be a data.frame")
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

            if( ! object@protein_id_column %in% colnames(object@protein_quant_table) ) {
                stop("Protein ID column must be in the protein data table")
            }

           if( ! object@sample_id %in% colnames(object@design_matrix) ) {
             stop("Sample ID column must be in the design matrix")
           }


           #Need to check the rows names in design matrix and the column names of the data table
           samples_in_protein_quant_table <- setdiff(colnames( object@protein_quant_table), object@protein_id_column)
           samples_in_design_matrix <- object@design_matrix |> pull( !! sym( object@sample_id ) )

           if( length( which( sort(samples_in_protein_quant_table) != sort(samples_in_design_matrix) )) > 0 ) {
             stop("Samples in protein data and design matrix must be the same" )
           }

         }

)
#'@export ProteinQuantitativeData

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @export
getProteinQuantitativeData <- function( peptide_object, protein_quant_table) {
  protein_obj <- ProteinQuantitativeData(
    # Protein Data Matrix Information
    protein_quant_table=protein_quant_table
    , protein_id_column= peptide_object@protein_id_column

    # Design Matrix Information
    , design_matrix = peptide_object@design_matrix
    , protein_id_table = data.frame()
    , sample_id=peptide_object@sample_id
    , group_id=peptide_object@group_id
    , technical_replicate_id=peptide_object@technical_replicate_id
    , args = peptide_object@args
  )
}

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
setGeneric( name ="setProteinData"
            , def=function( theObject, protein_quant_table, protein_id_column) {
                standardGeneric("setProteinData")
            })

#'@export
setMethod( f ="setProteinData"
           , signature = "ProteinQuantitativeData"
            , definition=function( theObject, protein_quant_table, protein_id_column ) {
              theObject@protein_quant_table <- protein_quant_table
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

            samples_id_vector <- setdiff(colnames(theObject@protein_quant_table), theObject@sample_id )

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
           , def=function( theObject
                           , proteins_intensity_cutoff_percentile = NULL
                           , proteins_proportion_of_samples_below_cutoff = NULL
                           , core_utilisation = NULL) {
             standardGeneric("proteinIntensityFiltering")
           })

#'@export
setMethod( f="proteinIntensityFiltering"
           , signature="ProteinQuantitativeData"
           , definition = function( theObject
                                    , proteins_intensity_cutoff_percentile = NULL
                                    , proteins_proportion_of_samples_below_cutoff = NULL
                                    , core_utilisation = NULL) {
             protein_quant_table <- theObject@protein_quant_table

             proteins_intensity_cutoff_percentile <- checkParamsObjectFunctionSimplify( theObject
                                                                               , "proteins_intensity_cutoff_percentile"
                                                                               , NULL)
             proteins_proportion_of_samples_below_cutoff <- checkParamsObjectFunctionSimplify( theObject
                                                                               , "proteins_proportion_of_samples_below_cutoff"
                                                                               , NULL)
             core_utilisation <- checkParamsObjectFunctionSimplify( theObject
                                                                    , "core_utilisation"
                                                                    , NA)

             theObject <- updateParamInObject(theObject, "proteins_intensity_cutoff_percentile")
             theObject <- updateParamInObject(theObject, "proteins_proportion_of_samples_below_cutoff")
             theObject <- updateParamInObject(theObject, "core_utilisation")


             data_long_cln <- protein_quant_table  |>
               pivot_longer( cols=!matches(theObject@protein_id_column)
                             , names_to = theObject@sample_id
                             , values_to = "log_values")  |>
               mutate( temp = "")

             min_peptide_intensity_threshold <- ceiling( quantile( data_long_cln$log_values, na.rm=TRUE, probs = c(proteins_intensity_cutoff_percentile) ))[1]

             peptide_normalised_pif_cln <- peptideIntensityFilteringHelper( data_long_cln
                                                                      , min_peptide_intensity_threshold = min_peptide_intensity_threshold
                                                                      , proteins_proportion_of_samples_below_cutoff = proteins_proportion_of_samples_below_cutoff
                                                                      , protein_id_column = !!sym( theObject@protein_id_column)
                                                                      , peptide_sequence_column = temp
                                                                      , peptide_quantity_column = log_values
                                                                      , core_utilisation = core_utilisation)


             theObject@protein_quant_table <- peptide_normalised_pif_cln |>
               dplyr::select( -temp) |>
               pivot_wider( id_cols = theObject@protein_id_column , names_from = !!sym(theObject@sample_id), values_from = log_values)

             theObject <- cleanDesignMatrix(theObject)

             updated_object <- theObject

          return(updated_object)
          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
setGeneric(name="removeProteinsWithOnlyOneReplicate"
           , def=function( theObject, core_utilisation = NULL, grouping_variable = NULL) {
             standardGeneric("removeProteinsWithOnlyOneReplicate")
           }
           , signature=c("theObject", "core_utilisation", "grouping_variable"))

#'@export
setMethod(f="removeProteinsWithOnlyOneReplicate"
          , definition=function( theObject, core_utilisation = NULL, grouping_variable = NULL) {
            protein_quant_table <- theObject@protein_quant_table
            samples_id_tbl <- theObject@design_matrix
            sample_id_tbl_sample_id_column <- theObject@sample_id
            # replicate_group_column <- theObject@technical_replicate_id
            protein_id_column <- theObject@protein_id_column

            input_table_sample_id_column <- theObject@sample_id
            quantity_column <- "log_values"

            grouping_variable <- checkParamsObjectFunctionSimplifyAcceptNull( theObject
                                                                              , "grouping_variable"
                                                                              , NULL)

            core_utilisation <- checkParamsObjectFunctionSimplify( theObject
                                                                   , "core_utilisation"
                                                                   , NA)

            theObject <- updateParamInObject(theObject, "grouping_variable")
            theObject <- updateParamInObject(theObject, "core_utilisation")

            data_long_cln <- protein_quant_table  |>
              pivot_longer( cols=!matches(protein_id_column)
                            , names_to = input_table_sample_id_column
                            , values_to = quantity_column)

            protein_quant_table <- removeProteinsWithOnlyOneReplicateHelper( input_table = data_long_cln
                                                                , samples_id_tbl = samples_id_tbl
                                                                , input_table_sample_id_column = !!sym( input_table_sample_id_column )
                                                                , sample_id_tbl_sample_id_column = !!sym( sample_id_tbl_sample_id_column)
                                                                , replicate_group_column = !!sym(grouping_variable)
                                                                , protein_id_column = !!sym( protein_id_column)
                                                                , quantity_column = !!sym( quantity_column)
                                                                , core_utilisation = core_utilisation)


            theObject@protein_quant_table <- protein_quant_table |>
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
           , def=function( theObject, grouping_variable, yaxis_limit = c()) {
             standardGeneric("plotRle")
           }
           , signature=c("theObject", "grouping_variable", "yaxis_limit"))


#'@export
setMethod(f="plotRle"
          , signature="ProteinQuantitativeData"
          , definition=function( theObject, grouping_variable, yaxis_limit = c()) {
            protein_quant_table <- theObject@protein_quant_table
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id

            frozen_protein_matrix <- protein_quant_table |>
              column_to_rownames(protein_id_column) |>
              as.matrix()

            design_matrix <- as.data.frame(design_matrix)
            rownames( design_matrix) <- design_matrix[,sample_id]

            # print( design_matrix)

            rowinfo_vector <- NA
            if( !is.na(grouping_variable)){
              rowinfo_vector <-  design_matrix[colnames(frozen_protein_matrix), grouping_variable]
            }

              rle_plot_before_cyclic_loess <- plotRleHelper( t(frozen_protein_matrix)
                                                       , rowinfo = rowinfo_vector
                                                       , yaxis_limit = yaxis_limit)

            return( rle_plot_before_cyclic_loess)

          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setGeneric(name="plotRleList"
           , def=function( theObject, list_of_columns, yaxis_limit = c()) {
             standardGeneric("plotRleList")
           }
           , signature=c("theObject", "list_of_columns", "yaxis_limit"))

#'@export
setMethod(f="plotRleList"
          , signature="ProteinQuantitativeData"
          , definition=function( theObject, list_of_columns, yaxis_limit = c()) {
            protein_quant_table <- theObject@protein_quant_table
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id

            frozen_protein_matrix <- protein_quant_table |>
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
                                                             , yaxis_limit = yaxis_limit)

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
           , def=function( theObject, grouping_variable, shape_variable = NULL, label_column, title, font_size ) {
             standardGeneric("plotPca")
           }
           , signature=c("theObject", "grouping_variable", "shape_variable", "label_column", "title", "font_size"))

#'@export
setMethod(f="plotPca"
          , signature="ProteinQuantitativeData"
          , definition=function( theObject, grouping_variable, shape_variable = NULL, label_column, title, font_size=8) {
            # Defensive checks
            if (!is.character(grouping_variable) || length(grouping_variable) != 1) {
              stop("grouping_variable must be a single character string")
            }
            
            if (!is.null(shape_variable) && (!is.character(shape_variable) || length(shape_variable) != 1)) {
              stop("shape_variable must be NULL or a single character string")
            }
            
            if (!grouping_variable %in% colnames(theObject@design_matrix)) {
              stop(sprintf("grouping_variable '%s' not found in design matrix", grouping_variable))
            }
            
            if (!is.null(shape_variable) && !shape_variable %in% colnames(theObject@design_matrix)) {
              stop(sprintf("shape_variable '%s' not found in design matrix", shape_variable))
            }

            protein_quant_table <- theObject@protein_quant_table
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id

            frozen_protein_matrix <- protein_quant_table |>
              column_to_rownames(protein_id_column) |>
              as.matrix()

            frozen_protein_matrix_pca <- frozen_protein_matrix
            frozen_protein_matrix_pca[!is.finite(frozen_protein_matrix_pca)] <- NA

            if(is.na(label_column) || label_column == "") {
              label_column <- ""
            }

            required_cols <- c(sample_id, grouping_variable)
            if (!is.null(shape_variable)) {
              required_cols <- c(required_cols, shape_variable)
            }
            missing_cols <- setdiff(required_cols, colnames(design_matrix))
            if (length(missing_cols) > 0) {
              stop(sprintf("Missing columns in design matrix: %s", paste(missing_cols, collapse = ", ")))
            }

            tryCatch({
              pca_plot <- plotPcaHelper(frozen_protein_matrix_pca,
                                           design_matrix,
                                           sample_id_column = sample_id,
                                           grouping_variable = grouping_variable,
                                           shape_variable = shape_variable,
                                           label_column = label_column,
                                           title = title,
                                           geom.text.size = font_size)
              return(pca_plot)
            }, error = function(e) {
              stop(sprintf("Error in plotPcaHelper: %s", e$message))
            })
          })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setGeneric(name="plotPcaList"
           , def=function( theObject, grouping_variables_list, label_column, title, font_size ) {
             standardGeneric("plotPcaList")
           }
           , signature=c("theObject", "grouping_variables_list", "label_column", "title", "font_size"))

#'@export
setMethod(f="plotPcaList"
          , signature="ProteinQuantitativeData"
          , definition=function( theObject, grouping_variables_list, label_column, title, font_size=8) {
            protein_quant_table <- theObject@protein_quant_table
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id

            frozen_protein_matrix <- protein_quant_table |>
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
                                                 , grouping_variables_list = grouping_variables_list
                                                 , label_column =  label_column
                                                 , title = title
                                                 , geom.text.size = font_size )

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
            protein_quant_table <- theObject@protein_quant_table
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id


            frozen_protein_matrix <- protein_quant_table |>
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
           , def=function( theObject,  tech_rep_num_column = NULL, tech_rep_remove_regex = NULL) {
             standardGeneric("proteinTechRepCorrelation")
           }
           , signature=c("theObject", "tech_rep_num_column", "tech_rep_remove_regex"))

#'@export
setMethod( f = "proteinTechRepCorrelation"
           , signature="ProteinQuantitativeData"
           , definition=function( theObject,  tech_rep_num_column = NULL, tech_rep_remove_regex = NULL ) {
            protein_quant_table <- theObject@protein_quant_table
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id
            tech_rep_column <- theObject@technical_replicate_id

            tech_rep_num_column <- checkParamsObjectFunctionSimplifyAcceptNull(theObject, "tech_rep_num_column", NULL)
            tech_rep_remove_regex <- checkParamsObjectFunctionSimplifyAcceptNull(theObject, "tech_rep_remove_regex", NULL)

            theObject <- updateParamInObject(theObject, "tech_rep_num_column")
            theObject <- updateParamInObject(theObject, "tech_rep_remove_regex")

            frozen_protein_matrix <- protein_quant_table |>
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
# Plot Pearson Correlation
#' @export
setGeneric(name="plotPearson",
           def=function(theObject, tech_rep_remove_regex, y_limit) {
             standardGeneric("plotPearson")
           },
           signature=c("theObject", "tech_rep_remove_regex", "y_limit"))

#' @export
setMethod(f="plotPearson",
          signature="ProteinQuantitativeData",
          definition=function(theObject, tech_rep_remove_regex = "pool") {
           
            correlation_vec <- pearsonCorForSamplePairs(theObject, tech_rep_remove_regex)
            
            pearson_plot <- correlation_vec |>
              ggplot(aes(pearson_correlation)) +
              geom_histogram(breaks = seq(min(round(correlation_vec$pearson_correlation - 0.5, 2), na.rm = TRUE), 1, 0.001)) +
              scale_y_continuous(breaks = seq(0, 4, 1), limits = c(0, 4)) +
              xlab("Pearson Correlation") +
              ylab("Counts") +
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank())
            
            return(pearson_plot)
          })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create empty QC Grid
#' @export
setClass("GridPlotData",
         slots = list(
           pca_plots = "list",
           rle_plots = "list",
           pearson_plots = "list"
         ))

#' @export
setGeneric("InitialiseGrid", function(dummy = NULL) {
  standardGeneric("InitialiseGrid")
})

#' @export
setMethod("InitialiseGrid", 
          signature(dummy = "ANY"),
          function(dummy = NULL) {
            new("GridPlotData",
                pca_plots = list(),
                rle_plots = list(),
                pearson_plots = list())
          })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Create a QC composite figure

#' @export
#' @export
setGeneric(name = "createGridQC",
           def = function(theObject, pca_titles, rle_titles, pearson_titles, save_path = NULL, file_name = "pca_rle_pearson_corr_plots_merged") {
             standardGeneric("createGridQC")
           },
           signature = c("theObject", "pca_titles", "rle_titles", "pearson_titles", "save_path", "file_name"))

#' @export
setMethod(f = "createGridQC",
          signature = "GridPlotData",
          definition = function(theObject, pca_titles, rle_titles, pearson_titles, save_path = NULL, file_name = "pca_rle_pearson_corr_plots_merged") {
            
            createPcaPlot <- function(plot, title) {
              plot +
                xlim(-40, 45) + ylim(-30, 25) + ggtitle(title) +
                theme(text = element_text(size = 15),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank())
            }
            
            createRlePlot <- function(plot, title) {
              plot + ggtitle(title) +
                theme(text = element_text(size = 15),
                      axis.text.x = element_blank(),
                      axis.ticks.x = element_blank())
            }
            
            createPearsonPlot <- function(plot, title) {
              plot + ggtitle(title) +
                theme(text = element_text(size = 15))
            }
            
            created_pca_plots <- mapply(createPcaPlot, theObject@pca_plots, pca_titles, SIMPLIFY = FALSE)
            created_rle_plots <- mapply(createRlePlot, theObject@rle_plots, rle_titles, SIMPLIFY = FALSE)
            created_pearson_plots <- mapply(createPearsonPlot, theObject@pearson_plots, pearson_titles, SIMPLIFY = FALSE)
            
            combined_plot <- (
              wrap_plots(created_pca_plots, ncol = 3) /
              wrap_plots(created_rle_plots, ncol = 3) /
              wrap_plots(created_pearson_plots, ncol = 3)
            ) +
              plot_layout(guides = 'collect')

            if (!is.null(save_path)) {
              sapply(c("png", "pdf", "svg"), function(ext) {
                ggsave(
                  plot = combined_plot,
                  filename = file.path(save_path, paste0(file_name, ".", ext)),
                  width = 14,
                  height = 14
                )
              })
              message(paste("Plots saved in", save_path))
            }
            
            return(combined_plot)
          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
## normalise between Arrays

#'@export
setGeneric(name="normaliseBetweenSamples"
           , def=function( theObject, normalisation_method = NULL) {
             standardGeneric("normaliseBetweenSamples")
           }
           , signature=c("theObject", "normalisation_method"))


#'@export
#'@param theObject Object of class ProteinQuantitativeData
#'@param normalisation_method Method to use for normalisation. Options are cyclicloess, quantile, scale, none
setMethod(f="normaliseBetweenSamples"
          , signature="ProteinQuantitativeData"
          , definition=function( theObject,  normalisation_method= NULL) {
            protein_quant_table <- theObject@protein_quant_table
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id

            normalisation_method <- checkParamsObjectFunctionSimplify( theObject
                                                         , "normalisation_method"
                                                         , "cyclicloess")

            theObject <- updateParamInObject(theObject, "normalisation_method")

            frozen_protein_matrix <- protein_quant_table |>
              column_to_rownames(protein_id_column) |>
              as.matrix()

            frozen_protein_matrix[!is.finite(frozen_protein_matrix)] <- NA

            normalised_frozen_protein_matrix <- frozen_protein_matrix

            print(paste0("normalisation_method = ", normalisation_method))

            switch( normalisation_method
                    , cyclicloess = {
                      normalised_frozen_protein_matrix <- normalizeCyclicLoess( frozen_protein_matrix )
                    }
                    , quantile = {
                      normalised_frozen_protein_matrix <- normalizeQuantiles( frozen_protein_matrix  )
                    }
                    , scale = {
                      normalised_frozen_protein_matrix <- normalizeMedianAbsValues( frozen_protein_matrix  )
                    }
                    , none = {
                      normalised_frozen_protein_matrix <- frozen_protein_matrix
                    }
            )

            normalised_frozen_protein_matrix[!is.finite(normalised_frozen_protein_matrix)] <- NA

            # normalised_frozen_protein_matrix_filt <- as.data.frame( normalised_frozen_protein_matrix ) |>
            #   dplyr::filter( if_all( everything(), \(x) { !is.na(x) } ) ) |>
            #   as.matrix()

            theObject@protein_quant_table <- normalised_frozen_protein_matrix |>
                      as.data.frame() |>
                      rownames_to_column(protein_id_column)

            theObject <- cleanDesignMatrix(theObject)

            updated_object <- theObject

            return(updated_object)

          })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setGeneric(name="pearsonCorForSamplePairs"
           , def=function( theObject,   tech_rep_remove_regex = NULL ) {
             standardGeneric("pearsonCorForSamplePairs")
           }
           , signature=c("theObject", "tech_rep_remove_regex"))

#'@export
setMethod(f="pearsonCorForSamplePairs"
          , signature="ProteinQuantitativeData"
          , definition=function( theObject, tech_rep_remove_regex = NULL ) {
            protein_quant_table <- theObject@protein_quant_table
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id
            replicate_group_column <- theObject@technical_replicate_id

            tech_rep_remove_regex <- checkParamsObjectFunctionSimplifyAcceptNull(theObject, "tech_rep_remove_regex", "pool")
            theObject <- updateParamInObject(theObject, "tech_rep_remove_regex")

            frozen_mat_pca_long <- protein_quant_table |>
              pivot_longer( cols=!matches(protein_id_column)
                            , values_to = "Protein.normalised"
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
                                                                                                       , peptide_normalised_column = "Protein.normalised")

            correlation_vec_before_cyclic_loess <- correlation_results_before_cyclic_loess |>
              dplyr::filter( !str_detect(!!sym(replicate_group_column), tech_rep_remove_regex )  )

           return( correlation_vec_before_cyclic_loess)
          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setGeneric(name="getNegCtrlProtAnova"
           , def=function( theObject
                           , ruv_grouping_variable  = NULL
                           , percentage_as_neg_ctrl  = NULL
                           , num_neg_ctrl  = NULL
                           , ruv_qval_cutoff = NULL
                           , ruv_fdr_method = NULL ) {
             standardGeneric("getNegCtrlProtAnova")
           }
           , signature=c("theObject", "ruv_grouping_variable", "num_neg_ctrl", "ruv_qval_cutoff", "ruv_fdr_method"))

#'@export
setMethod(f="getNegCtrlProtAnova"
          , signature="ProteinQuantitativeData"
          , definition=function( theObject
                                 , ruv_grouping_variable = NULL
                                 , percentage_as_neg_ctrl = NULL
                                 , num_neg_ctrl = NULL
                                 , ruv_qval_cutoff = NULL
                                 , ruv_fdr_method = NULL ) {

            protein_quant_table <- theObject@protein_quant_table
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            group_id <- theObject@group_id
            sample_id <- theObject@sample_id
            replicate_group_column <- theObject@technical_replicate_id

            normalised_frozen_protein_matrix_filt <- protein_quant_table |>
              column_to_rownames(protein_id_column) |>
              as.matrix()

            ruv_grouping_variable <- checkParamsObjectFunctionSimplify( theObject, "ruv_grouping_variable", "replicates")
            percentage_as_neg_ctrl <- checkParamsObjectFunctionSimplify( theObject, "percentage_as_neg_ctrl", 10)
            num_neg_ctrl <- checkParamsObjectFunctionSimplify( theObject
                                                               , "num_neg_ctrl"
                                                               , round(nrow( theObject@protein_quant_table) * percentage_as_neg_ctrl / 100, 0))
            ruv_qval_cutoff <- checkParamsObjectFunctionSimplify( theObject, "ruv_qval_cutoff", 0.05)
            ruv_fdr_method <- checkParamsObjectFunctionSimplify( theObject, "ruv_fdr_method", "BH")

            theObject <- updateParamInObject(theObject, "ruv_grouping_variable")
            theObject <- updateParamInObject(theObject, "percentage_as_neg_ctrl")
            theObject <- updateParamInObject(theObject, "num_neg_ctrl")
            theObject <- updateParamInObject(theObject, "ruv_qval_cutoff")
            theObject <- updateParamInObject(theObject, "ruv_fdr_method")

            control_genes_index <- getNegCtrlProtAnovaHelper( normalised_frozen_protein_matrix_filt[,design_matrix |> pull(!!sym(sample_id)) ]
                                                        , design_matrix = design_matrix |>
                                                          column_to_rownames(sample_id) |>
                                                          dplyr::select( -!!sym(group_id))
                                                        , grouping_variable = ruv_grouping_variable
                                                        , percentage_as_neg_ctrl = percentage_as_neg_ctrl
                                                        , num_neg_ctrl = num_neg_ctrl
                                                        , ruv_qval_cutoff = ruv_qval_cutoff
                                                        , ruv_fdr_method = ruv_fdr_method )

            return(control_genes_index)
          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@description Sort proteins by their coefficient of variation and take the top N with lowest coefficient of variation
#'@export
setGeneric(name="getLowCoefficientOfVariationProteins"
           , def=function( theObject
                           , percentage_as_neg_ctrl = NULL
                           , num_neg_ctrl = NULL ) {
             standardGeneric("getLowCoefficientOfVariationProteins")
           }
           , signature=c("theObject", "percentage_as_neg_ctrl", "num_neg_ctrl"))



#'@export
setMethod( f = "getLowCoefficientOfVariationProteins"
           , signature="ProteinQuantitativeData"
           , definition=function( theObject
                                  , percentage_as_neg_ctrl = NULL
                                  , num_neg_ctrl = NULL) {

             percentage_as_neg_ctrl <- checkParamsObjectFunctionSimplify( theObject, "percentage_as_neg_ctrl", 10)
             num_neg_ctrl <- checkParamsObjectFunctionSimplify( theObject
                                                                , "num_neg_ctrl"
                                                                , round(nrow( theObject@protein_quant_table) * percentage_as_neg_ctrl / 100, 0))

             theObject <- updateParamInObject(theObject, "percentage_as_neg_ctrl")
             theObject <- updateParamInObject(theObject, "num_neg_ctrl")

  list_of_control_genes <- theObject@protein_quant_table |>
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

  control_gene_index_helper <- theObject@protein_quant_table |>
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
           , def=function( theObject, ctrl= NULL, num_components_to_impute=NULL, ruv_grouping_variable = NULL ) {
             standardGeneric("ruvCancor")
           }
           , signature=c("theObject", "ctrl", "num_components_to_impute", "ruv_grouping_variable"))

#'@export
setMethod( f = "ruvCancor"
           , signature="ProteinQuantitativeData"
           , definition=function( theObject, ctrl= NULL, num_components_to_impute=NULL, ruv_grouping_variable = NULL) {
             protein_quant_table <- theObject@protein_quant_table
             protein_id_column <- theObject@protein_id_column
             design_matrix <- theObject@design_matrix
             group_id <- theObject@group_id
             sample_id <- theObject@sample_id
             replicate_group_column <- theObject@technical_replicate_id

             ctrl <- checkParamsObjectFunctionSimplify( theObject, "ctrl", NULL)
             num_components_to_impute <- checkParamsObjectFunctionSimplify( theObject, "num_components_to_impute", 2)
             ruv_grouping_variable <- checkParamsObjectFunctionSimplify( theObject, "ruv_grouping_variable", NULL)

             theObject <- updateParamInObject(theObject, "ctrl")
             theObject <- updateParamInObject(theObject, "num_components_to_impute")
             theObject <- updateParamInObject(theObject, "ruv_grouping_variable")

             if(! ruv_grouping_variable %in% colnames(design_matrix)) {
               stop( paste0("The 'ruv_grouping_variable = "
                            , ruv_grouping_variable
                            , "' is not a column in the design matrix.") )
             }

             if( is.na(num_components_to_impute) || num_components_to_impute < 1) {
               stop(paste0("The num_components_to_impute = ", num_components_to_impute, " value is invalid."))
             }

             if( length( ctrl) < 5 ) {
               stop(paste0( "The number of negative control molecules entered is less than 5. Please check the 'ctl' parameter."))
             }

             normalised_frozen_protein_matrix_filt <- protein_quant_table |>
               column_to_rownames(protein_id_column) |>
               as.matrix()

             Y <-  t( normalised_frozen_protein_matrix_filt[,design_matrix |> pull(!!sym(sample_id))])
             if( length(which( is.na(normalised_frozen_protein_matrix_filt) )) > 0 ) {
               Y <- impute.nipals( t( normalised_frozen_protein_matrix_filt[,design_matrix |> pull(!!sym(sample_id))])
                                   , ncomp=num_components_to_impute)
             }

             cancorplot_r2 <- ruv_cancorplot( Y ,
                                              X = design_matrix |>
                                                pull(!!sym(ruv_grouping_variable)),
                                              ctl = ctrl)
             cancorplot_r2


           })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#'@export
setGeneric(name="getRuvIIIReplicateMatrix"
           , def=function( theObject,  ruv_grouping_variable = NULL) {
             standardGeneric("getRuvIIIReplicateMatrix")
           }
           , signature=c("theObject", "ruv_grouping_variable"))

#'@export
setMethod( f = "getRuvIIIReplicateMatrix"
           , signature="ProteinQuantitativeData"
           , definition=function( theObject, ruv_grouping_variable = NULL) {
             protein_quant_table <- theObject@protein_quant_table
             protein_id_column <- theObject@protein_id_column
             design_matrix <- theObject@design_matrix
             group_id <- theObject@group_id
             sample_id <- theObject@sample_id
             replicate_group_column <- theObject@technical_replicate_id

             ruv_grouping_variable <- checkParamsObjectFunctionSimplify( theObject, "ruv_grouping_variable", NULL)

             theObject <- updateParamInObject(theObject, "ruv_grouping_variable")

             ruvIII_replicates_matrix <- getRuvIIIReplicateMatrixHelper( design_matrix
                                                                   , !!sym(sample_id)
                                                                   , !!sym(ruv_grouping_variable))
             return( ruvIII_replicates_matrix)
           })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#'@export
setGeneric(name="ruvIII_C_Varying"
           , def=function( theObject, ruv_grouping_variable = NULL, ruv_number_k = NULL, ctrl = NULL)  {
             standardGeneric("ruvIII_C_Varying")
           }
           , signature=c("theObject", "ruv_grouping_variable", "ruv_number_k", "ctrl"))

#'@export
setMethod( f = "ruvIII_C_Varying"
           , signature="ProteinQuantitativeData"
           , definition=function( theObject, ruv_grouping_variable = NULL, ruv_number_k = NULL, ctrl = NULL) {
             protein_quant_table <- theObject@protein_quant_table
             protein_id_column <- theObject@protein_id_column
             design_matrix <- theObject@design_matrix
             group_id <- theObject@group_id
             sample_id <- theObject@sample_id
             replicate_group_column <- theObject@technical_replicate_id


             ruv_grouping_variable <- checkParamsObjectFunctionSimplify( theObject, "ruv_grouping_variable", NULL)
             k <- checkParamsObjectFunctionSimplify( theObject, "ruv_number_k", NULL)
             ctrl <- checkParamsObjectFunctionSimplify( theObject, "ctrl", NULL)

             theObject <- updateParamInObject(theObject, "ruv_grouping_variable")
             theObject <- updateParamInObject(theObject, "ruv_number_k")
             theObject <- updateParamInObject(theObject, "ctrl")

             normalised_frozen_protein_matrix_filt <- protein_quant_table |>
               column_to_rownames(protein_id_column) |>
               as.matrix()

             Y <-  t( normalised_frozen_protein_matrix_filt[,design_matrix |> pull(!!sym(sample_id))])

             M <- getRuvIIIReplicateMatrixHelper( design_matrix
                                            , !!sym(sample_id)
                                            , !!sym(ruv_grouping_variable))

             cln_mat <- RUVIII_C_Varying( k = ruv_number_k
                                          , Y = Y
                                          , M = M
                                          , toCorrect = colnames(Y)
                                          , potentialControls = names( ctrl[which(ctrl)] ) )

             # Remove samples with no values
             cln_mat_2 <- cln_mat[rowSums(is.na(cln_mat) | is.nan(cln_mat)) != ncol(cln_mat),]

             # Remove proteins with no values
             cln_mat_3 <- t(cln_mat_2)
             cln_mat_4 <- cln_mat_3[rowSums(is.na(cln_mat_3) | is.nan(cln_mat_3)) != ncol(cln_mat_3),]

             ruv_normalised_results_cln <- cln_mat_4 |>
               as.data.frame() |>
               rownames_to_column(protein_id_column)

             theObject@protein_quant_table <- ruv_normalised_results_cln

             theObject <- cleanDesignMatrix(theObject)

             return( theObject )
          })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
setGeneric(name="removeRowsWithMissingValuesPercent"
           , def=function( theObject
                           , ruv_grouping_variable = NULL
                           , groupwise_percentage_cutoff = NULL
                           , max_groups_percentage_cutoff = NULL
                           , proteins_intensity_cutoff_percentile = NULL ) {
             standardGeneric("removeRowsWithMissingValuesPercent")
           }
           , signature=c("theObject"
                         , "ruv_grouping_variable"
                         , "groupwise_percentage_cutoff"
                         , "max_groups_percentage_cutoff"
                         , "proteins_intensity_cutoff_percentile" ))

#'@export
setMethod( f = "removeRowsWithMissingValuesPercent"
           , signature="ProteinQuantitativeData"
           , definition=function( theObject
                                  , ruv_grouping_variable = NULL
                                  , groupwise_percentage_cutoff = NULL
                                  , max_groups_percentage_cutoff = NULL
                                  , proteins_intensity_cutoff_percentile = NULL) {

             protein_quant_table <- theObject@protein_quant_table
             protein_id_column <- theObject@protein_id_column
             design_matrix <- theObject@design_matrix
             group_id <- theObject@group_id
             sample_id <- theObject@sample_id
             replicate_group_column <- theObject@technical_replicate_id

             # print(groupwise_percentage_cutoff)
             # print(min_protein_intensity_threshold )

             ruv_grouping_variable <- checkParamsObjectFunctionSimplify(theObject
                                                                        , "ruv_grouping_variable"
                                                                        , NULL)
             groupwise_percentage_cutoff <- checkParamsObjectFunctionSimplify(theObject
                                                                              , "groupwise_percentage_cutoff"
                                                                              , 50)
             max_groups_percentage_cutoff <- checkParamsObjectFunctionSimplify(theObject
                                                                               , "max_groups_percentage_cutoff"
                                                                               , 50)
             proteins_intensity_cutoff_percentile <- checkParamsObjectFunctionSimplify(theObject
                                                                                   , "proteins_intensity_cutoff_percentile"
                                                                                   , 1)

             theObject <- updateParamInObject(theObject, "ruv_grouping_variable")
             theObject <- updateParamInObject(theObject, "groupwise_percentage_cutoff")
             theObject <- updateParamInObject(theObject, "max_groups_percentage_cutoff")
             theObject <- updateParamInObject(theObject, "proteins_intensity_cutoff_percentile")


             theObject@protein_quant_table <- removeRowsWithMissingValuesPercentHelper( protein_quant_table
                                                                           , cols= !matches(protein_id_column)
                                                                           , design_matrix = design_matrix
                                                                           , sample_id = !!sym(sample_id)
                                                                           , row_id = !!sym(protein_id_column)
                                                                           , grouping_variable = !!sym(ruv_grouping_variable)
                                                                           , groupwise_percentage_cutoff = groupwise_percentage_cutoff
                                                                           , max_groups_percentage_cutoff = max_groups_percentage_cutoff
                                                                           , proteins_intensity_cutoff_percentile = proteins_intensity_cutoff_percentile
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
#'@param theObject The object to be processed
#'@param design_matrix_columns The columns to be used in the design matrix
#'@param protein_id_column The column name of the protein id
#'@param sample_id The column name of the sample id
#'@param replicate_group_column The column name of the technical replicate id
setMethod( f = "averageTechReps"
           , signature="ProteinQuantitativeData"
           , definition=function( theObject, design_matrix_columns=c()  ) {

             protein_quant_table <- theObject@protein_quant_table
             protein_id_column <- theObject@protein_id_column
             design_matrix <- theObject@design_matrix
             group_id <- theObject@group_id
             sample_id <- theObject@sample_id
             replicate_group_column <- theObject@technical_replicate_id

             theObject@protein_quant_table <- protein_quant_table |>
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
setGeneric(name="preservePeptideNaValues"
           , def=function( peptide_obj, protein_obj)  {
             standardGeneric("preservePeptideNaValues")
           }
           , signature=c("peptide_obj", "protein_obj" ))

#'@export
setMethod( f = "preservePeptideNaValues"
           , signature=c( "PeptideQuantitativeData", "ProteinQuantitativeData" )
           , definition= function( peptide_obj, protein_obj) {
             preservePeptideNaValuesHelper( peptide_obj, protein_obj)
           })


preservePeptideNaValuesHelper <- function( peptide_obj, protein_obj) {

  sample_id_column <- peptide_obj@sample_id
  protein_id_column <- peptide_obj@protein_id_column

  check_peptide_value <- peptide_obj@peptide_data |>
    group_by( !!sym( sample_id_column), !!sym(protein_id_column) ) |>
    summarise( Peptide.Normalised = sum( Peptide.Normalised, na.rm=TRUE)
               , is_na = sum( is.na(Peptide.Normalised ))
               , num_values = n() ) |>
    mutate( Peptide.Normalised = if_else( is_na == num_values, NA_real_, Peptide.Normalised)) |>
    ungroup() |>
    arrange( !!sym( sample_id_column)) |>
    pivot_wider( id_cols = !!sym(protein_id_column)
                 , names_from = !!sym(sample_id_column)
                 , values_from = Peptide.Normalised
                 , values_fill = NA_real_)

  check_peptide_value_cln <- check_peptide_value[rownames(protein_obj@protein_quant_table)
                                                 , colnames(  protein_obj@protein_quant_table)]

  if( length( which (rownames(protein_obj@protein_quant_table) ==  rownames(check_peptide_value_cln))) != nrow(check_peptide_value_cln) ) {
    stop("The rows in the protein object and the peptide object do not match")
  }

  if( length( which( colnames( protein_obj@protein_quant_table) == colnames(check_peptide_value_cln) )) != ncol(check_peptide_value_cln) ) {
    stop("The columns in the protein object and the peptide object do not match")
  }

  protein_obj@protein_quant_table [is.na(check_peptide_value_cln)] <- NA

  protein_obj
}
##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#'@export
setGeneric(name="chooseBestProteinAccession"
           , def=function(theObject, delim=NULL, seqinr_obj=NULL, seqinr_accession_column=NULL, replace_zero_with_na = NULL, aggregation_method = NULL) {
             standardGeneric("chooseBestProteinAccession")
           }
           , signature=c("theObject", "delim", "seqinr_obj", "seqinr_accession_column"))

#'@export
#'@param theObject The object of class ProteinQuantitativeData
#'@param delim The delimiter used to split the protein accessions
#'@param seqinr_obj The object of class Seqinr::seqinr
#'@param seqinr_accession_column The column in the seqinr object that contains the protein accessions
#'@param replace_zero_with_na Replace zero values with NA
#'@param aggregation_method Method to aggregate protein values: "sum", "mean", or "median" (default: "sum")
setMethod(f = "chooseBestProteinAccession"
          , signature="ProteinQuantitativeData"
          , definition=function(theObject, delim=NULL, seqinr_obj=NULL
                              , seqinr_accession_column=NULL
                              , replace_zero_with_na = NULL
                              , aggregation_method = NULL) {
            
            protein_quant_table <- theObject@protein_quant_table
            protein_id_column <- theObject@protein_id_column
            
            delim <- checkParamsObjectFunctionSimplify(theObject, "delim",  default_value =  " |;|:|\\|")
            seqinr_obj <- checkParamsObjectFunctionSimplify(theObject, "seqinr_obj",  default_value = NULL)
            seqinr_accession_column <- checkParamsObjectFunctionSimplify(theObject
                                                                       , "seqinr_accession_column"
                                                                       , default_value = NULL)
            replace_zero_with_na <- checkParamsObjectFunctionSimplify(theObject
                                                                    , "replace_zero_with_na"
                                                                    , default_value = FALSE)
            aggregation_method <- checkParamsObjectFunctionSimplify(theObject
                                                                  , "aggregation_method"
                                                                  , default_value = "sum")
            
            if (!aggregation_method %in% c("sum", "mean", "median")) {
              stop("aggregation_method must be one of: 'sum', 'mean', 'median'")
            }
            
            theObject <- updateParamInObject(theObject, "delim")
            theObject <- updateParamInObject(theObject, "seqinr_obj")
            theObject <- updateParamInObject(theObject, "seqinr_accession_column")
            theObject <- updateParamInObject(theObject, "replace_zero_with_na")
            theObject <- updateParamInObject(theObject, "aggregation_method")
            
            evidence_tbl_cleaned <- protein_quant_table |>
              distinct() |>
              mutate(row_id = row_number() -1)
            
            accession_gene_name_tbl <- chooseBestProteinAccessionHelper(input_tbl = evidence_tbl_cleaned,
                                                                      acc_detail_tab = seqinr_obj,
                                                                      accessions_column = !!sym(protein_id_column),
                                                                      row_id_column = seqinr_accession_column,
                                                                      group_id = row_id,
                                                                      delim = ";")
            
            protein_log2_quant_cln <- evidence_tbl_cleaned |>
              left_join(accession_gene_name_tbl |>
                         dplyr::distinct(row_id, !!sym(as.character(seqinr_accession_column)))
                       , by = join_by(row_id)) |>
              mutate(!!sym(theObject@protein_id_column) := !!sym(as.character(seqinr_accession_column))) |>
              dplyr::select(-row_id, -!!sym(as.character(seqinr_accession_column)))
            
            protein_id_table <- evidence_tbl_cleaned |>
              left_join(accession_gene_name_tbl |>
                         dplyr::distinct(row_id, !!sym(as.character(seqinr_accession_column)))
                       , by = join_by(row_id)) |>
              distinct(uniprot_acc, !!sym(protein_id_column)) |>
              mutate(!!sym(paste0(protein_id_column, "_list")) := !!sym(protein_id_column)) |>
              mutate(!!sym(protein_id_column) := !!sym("uniprot_acc")) |>
              distinct(!!sym(protein_id_column), !!sym(paste0(protein_id_column, "_list"))) |>
              group_by(!!sym(protein_id_column)) |>
              summarise(!!sym(paste0(protein_id_column, "_list")) := paste(!!sym(paste0(protein_id_column, "_list")), collapse = ";")) |>
              ungroup() |>
              mutate(!!sym(paste0(protein_id_column, "_list")) := purrr::map_chr(!!sym(paste0(protein_id_column, "_list"))
                                                                                , \(x){ paste(unique(sort(str_split(x, ";")[[1]])), collapse=";") }))
            
            summed_data <- protein_log2_quant_cln |>
              mutate(!!sym(protein_id_column) := purrr::map_chr(!!sym(protein_id_column), \(x){ str_split(x, delim)[[1]][1] })) |>
              pivot_longer(
                cols = !matches(protein_id_column),
                names_to = "sample_id",
                values_to = "temporary_values_choose_accession"
              ) |>
              group_by(!!sym(protein_id_column), sample_id) |>
              summarise(
                is_na = sum(is.na(temporary_values_choose_accession)),
                temporary_values_choose_accession = case_when(
                  all(is.na(temporary_values_choose_accession)) ~ NA_real_,
                  aggregation_method == "sum" ~ sum(temporary_values_choose_accession, na.rm = TRUE),
                  aggregation_method == "mean" ~ mean(temporary_values_choose_accession, na.rm = TRUE),
                  aggregation_method == "median" ~ median(temporary_values_choose_accession, na.rm = TRUE)
                ),
                num_values = n()
              ) |>
              ungroup() |>
              pivot_wider(
                id_cols = !!sym(protein_id_column),
                names_from = sample_id,
                values_from = temporary_values_choose_accession,
                values_fill = NA_real_
              )
            
            if(replace_zero_with_na == TRUE) {
              summed_data[is.na(summed_data)] <- NA
            }
            
            protein_id_table <- rankProteinAccessionHelper(input_tbl = protein_id_table,
                                                         acc_detail_tab = seqinr_obj,
                                                         accessions_column = !!sym(paste0(protein_id_column, "_list")),
                                                         row_id_column = seqinr_accession_column,
                                                         group_id = !!sym(protein_id_column),
                                                         delim = ";") |>
              dplyr::rename(!!sym(paste0(protein_id_column, "_list")) := seqinr_accession_column) |>
              dplyr::select(-num_gene_names, -gene_names, -is_unique)
            
            theObject@protein_id_table <- protein_id_table
            theObject@protein_quant_table <- summed_data[, colnames(protein_quant_table)]
            
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

             protein_quant_table <- theObject@protein_quant_table
             protein_id_column <- theObject@protein_id_column

             protein_log2_quant_cln <- protein_quant_table |>
               mutate( !!sym(protein_id_column) := str_split_i(!!sym( protein_id_column), delim, 1 ) ) |>
               group_by( !!sym(protein_id_column) ) |>
               summarise ( across( matches(quant_columns_pattern)
                                   , \(x){ if(islogged==TRUE) {
                                                log2(sum(2^x, na.rm = TRUE))
                                           } else {
                                                sum(x, na.rm = TRUE)
                                           }
                                         } )) |>
               ungroup()

             theObject@protein_quant_table <- protein_log2_quant_cln

             theObject

           })



##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setGeneric(name="filterSamplesByProteinCorrelationThreshold"
           , def=function( theObject, pearson_correlation_per_pair = NULL, min_pearson_correlation_threshold = NULL ) {
             standardGeneric("filterSamplesByProteinCorrelationThreshold")
           }
           , signature=c("theObject", "pearson_correlation_per_pair", "min_pearson_correlation_threshold" ))

#'@export
setMethod( f = "filterSamplesByProteinCorrelationThreshold"
           , signature="ProteinQuantitativeData"
           , definition=function( theObject, pearson_correlation_per_pair = NULL, min_pearson_correlation_threshold = NULL  ) {

             pearson_correlation_per_pair <- checkParamsObjectFunctionSimplify( theObject
                                                                           , "pearson_correlation_per_pair"
                                                                           , default_value = NULL)
             min_pearson_correlation_threshold <- checkParamsObjectFunctionSimplify( theObject
                                                                                , "min_pearson_correlation_threshold"
                                                                                , default_value = 0.75)

             theObject <- updateParamInObject(theObject, "pearson_correlation_per_pair")
             theObject <- updateParamInObject(theObject, "min_pearson_correlation_threshold")

             filtered_table <- filterSamplesByProteinCorrelationThresholdHelper (
               pearson_correlation_per_pair
               , protein_intensity_table = theObject@protein_quant_table
               , min_pearson_correlation_threshold = min_pearson_correlation_threshold
               , filename_column_x = !!sym( paste0( theObject@sample_id, ".x") )
               , filename_column_y = !!sym( paste0( theObject@sample_id, ".y") )
               , protein_id_column = theObject@protein_id_column
               , correlation_column = pearson_correlation )

             theObject@protein_quant_table <- filtered_table

             theObject <- cleanDesignMatrix(theObject)

             theObject
             })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# I want to input two protein data objects and compare them,
# to see how the number of proteins changes and how the number of samples changed
# Use set diff or set intersect to compare the list of proteins and samples in the two objects
#' @export
compareTwoProteinDataObjects <- function( object_a, object_b) {


  object_a_proteins <- object_a@protein_quant_table |>
    distinct(!!sym(object_a@protein_id_column)) |>
    pull(!!sym(object_a@protein_id_column))

  object_b_proteins <- object_b@protein_quant_table |>
    distinct(!!sym(object_b@protein_id_column)) |>
    pull(!!sym(object_b@protein_id_column))

  object_a_samples <- object_a@design_matrix |>
    distinct(!!sym(object_a@sample_id)) |>
    pull(!!sym(object_a@sample_id))

  object_b_samples <- object_b@design_matrix |>
    distinct(!!sym(object_b@sample_id)) |>
    pull(!!sym(object_b@sample_id))


  proteins_in_a_not_b <- length( setdiff( object_a_proteins, object_b_proteins))
  proteins_intersect_a_and_b <- length( intersect( object_a_proteins, object_b_proteins))
  proteins_in_b_not_a <- length( setdiff( object_b_proteins, object_a_proteins))


  samples_in_a_not_b <- length( setdiff( object_a_samples, object_b_samples))
  samples_intersect_a_and_b <- length( intersect( object_a_samples, object_b_samples))
  samples_in_b_not_a <- length( setdiff( object_b_samples, object_a_samples))

  comparisons_list <- list( proteins = list( in_a_not_b = proteins_in_a_not_b
                                               , intersect_a_and_b = proteins_intersect_a_and_b
                                               , in_b_not_a = proteins_in_b_not_a)
                            , samples = list( in_a_not_b = samples_in_a_not_b
                                              , intersect_a_and_b = samples_intersect_a_and_b
                                              , in_b_not_a = samples_in_b_not_a)
  )

  comparison_tibble <- comparisons_list |>
    purrr::map_df( tibble::as_tibble) |>
    add_column( Levels = c( "proteins", "samples")) |>
    relocate( Levels, .before="in_a_not_b")

  comparison_tibble


}

#'@export
summariseProteinObject <- function ( theObject) {
  num_proteins <- theObject@protein_quant_table |>
    distinct(!!sym(theObject@protein_id_column)) |>
    pull(!!sym(theObject@protein_id_column))

  num_samples <- theObject@design_matrix |>
    distinct(!!sym(theObject@sample_id)) |>
    pull(!!sym(theObject@sample_id))

  summary_list <- list( num_proteins = length(num_proteins)
       , num_samples = length(num_samples))

  summary_list

}
