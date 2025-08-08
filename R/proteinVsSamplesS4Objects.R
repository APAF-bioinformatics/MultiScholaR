

#' An S4 Class for Protein-Level Quantitative Proteomics Data
#'
#' @description The `ProteinQuantitativeData` class is a container for protein-level
#' abundance data, typically derived from peptide-level data after a roll-up
#' process. It stores the protein quantification table, experimental design, and
#' related metadata.
#'
#' @slot protein_quant_table A data frame or matrix with proteins in rows and samples
#'   in columns, containing the quantitative abundance data.
#' @slot protein_id_column The name of the column in `protein_quant_table` that
#'   contains the protein identifiers.
#' @slot design_matrix A data frame describing the experimental design.
#' @slot protein_id_table A data frame mapping chosen protein IDs to all possible IDs
#'   from the original protein groups.
#' @slot sample_id The name of the column for unique sample identifiers.
#' @slot group_id The name of the column for experimental group labels.
#' @slot technical_replicate_id The name of the column for technical replicate identifiers.
#' @slot args A list for storing various parameters and settings.
#'
#' @return An object of class `ProteinQuantitativeData`.
#' @export
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
           samples_in_design_matrix <- object@design_matrix |> dplyr::pull( !! sym( object@sample_id ) )

           if( length( which( sort(samples_in_protein_quant_table) != sort(samples_in_design_matrix) )) > 0 ) {
             stop("Samples in protein data and design matrix must be the same" )
           }

         }

)
#'@export ProteinQuantitativeData

#' @title Initialize Method for ProteinQuantitativeData
#' @description This method is called upon creation of a new `ProteinQuantitativeData`
#' object. It ensures that the sample ID column in the design matrix is of type character.
#'
#' @param .Object The object being created.
#' @param ... Arguments passed to the constructor.
#'
#' @return An initialized `ProteinQuantitativeData` object.
setMethod("initialize", "ProteinQuantitativeData",
  function(.Object, ...) {
    # Capture all arguments passed to the constructor
    args_list <- list(...)

    # If design_matrix and sample_id are provided in the arguments
    if (!is.null(args_list$design_matrix) && !is.null(args_list$sample_id)) {
      if (args_list$sample_id %in% names(args_list$design_matrix)) {
        args_list$design_matrix[[args_list$sample_id]] <- as.character(args_list$design_matrix[[args_list$sample_id]])
        logger::log_debug(
          "Initialize ProteinQuantitativeData: Converted column '{args_list$sample_id}' in design_matrix to character."
        )
      } else {
        # This case should ideally be caught by validity, but good to be aware
        logger::log_warn(
          "Initialize ProteinQuantitativeData: Sample ID column '{args_list$sample_id}' not found in the provided design_matrix during initialization."
        )
      }
    }
    
    # Reconstruct the call to the default initializer with potentially modified args_list
    # This uses do.call to pass arguments correctly after modification
    .Object <- do.call(callNextMethod, c(list(.Object), args_list))
    
    # Explicitly return the object
    return(.Object)
  }
)

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @title Create a ProteinQuantitativeData Object from a PeptideDataObject
#' @description A constructor function that creates a `ProteinQuantitativeData` object,
#' typically after rolling up peptide data to the protein level. It inherits metadata
#' from a `PeptideQuantitativeData` object.
#'
#' @param peptide_object An object of class `PeptideQuantitativeData` from which to
#'   copy the design matrix and other metadata.
#' @param protein_quant_table A data frame of protein-level quantification data.
#'
#' @return An object of class `ProteinQuantitativeData`.
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
#' @title Set Protein Data
#' @description A method to update the main protein quantification data and the
#' protein ID column name in a `ProteinQuantitativeData` object.
#'
#' @param theObject An object of class `ProteinQuantitativeData`.
#' @param protein_quant_table A new data frame to replace the `protein_quant_table` slot.
#' @param protein_id_column The name of the new protein ID column.
#'
#' @return The updated `ProteinQuantitativeData` object.
#' @exportMethod setProteinData
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

#' @title Clean the Design Matrix for a ProteinData Object
#' @description This method filters the design matrix to ensure it only contains
#' metadata for samples present in the `protein_quant_table`.
#'
#' @param theObject An object of class `ProteinQuantitativeData`.
#'
#' @return The modified `ProteinQuantitativeData` object.
#' @exportMethod cleanDesignMatrix
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

#' @title Filter Proteins by Intensity
#' @description Filters out proteins that have low abundance across a large proportion of samples.
#'
#' @param theObject An object of class `ProteinQuantitativeData`.
#' @param proteins_intensity_cutoff_percentile The percentile used to calculate the minimum intensity threshold.
#' @param proteins_proportion_of_samples_below_cutoff The proportion of samples that must be above the threshold.
#' @param core_utilisation The number of cores for parallel processing.
#'
#' @return The modified `ProteinQuantitativeData` object with low-intensity proteins removed.
#' @exportMethod proteinIntensityFiltering
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
#' @title Remove Proteins with Only One Replicate in a Group
#' @description This method filters out proteins that are only observed in a single
#' replicate within their experimental group.
#'
#' @param theObject An object of class `ProteinQuantitativeData`.
#' @param core_utilisation The number of cores for parallel processing.
#' @param grouping_variable The unquoted column name in the design matrix that defines the groups.
#'
#' @return The modified `ProteinQuantitativeData` object.
#' @export
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



#' @title Plot Relative Log Expression (RLE) for Protein Data
#' @description Generates an RLE plot for the protein-level data in the object.
#'
#' @param theObject An object of class `ProteinQuantitativeData`.
#' @param grouping_variable The column name in the design matrix to use for coloring the boxplots.
#' @param yaxis_limit A numeric vector of length 2 for the y-axis limits.
#' @param sample_label An optional column name in the design matrix to use for sample labels.
#'
#' @return A ggplot object representing the RLE plot.
#' @exportMethod plotRle
setMethod(f="plotRle"
          , signature="ProteinQuantitativeData"
          , definition=function( theObject, grouping_variable, yaxis_limit = c(), sample_label = NULL) {
            protein_quant_table <- theObject@protein_quant_table
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id

            frozen_protein_matrix <- protein_quant_table |>
              column_to_rownames(protein_id_column) |>
              as.matrix()

            design_matrix <- as.data.frame(design_matrix)

            if(!is.null(sample_label)) {
              if ( sample_label %in% colnames(design_matrix)) {
                rownames( design_matrix) <- design_matrix[,sample_label]
                colnames( frozen_protein_matrix ) <- design_matrix[,sample_label]

              } } else {
                rownames( design_matrix) <- design_matrix[,sample_id]
              }

            # print( design_matrix)

            rowinfo_vector <- NA
            if( !is.na(grouping_variable)){
              rowinfo_vector <-  design_matrix[colnames(frozen_protein_matrix), grouping_variable]
            }

            print(rownames( design_matrix))
            print(colnames( frozen_protein_matrix))
            print(rowinfo_vector)
              rle_plot_before_cyclic_loess <- plotRleHelper( t(frozen_protein_matrix)
                                                       , rowinfo = rowinfo_vector
                                                       , yaxis_limit = yaxis_limit)

            return( rle_plot_before_cyclic_loess)

          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#' @title Plot a List of RLE Plots
#' @description Generates a list of RLE plots, one for each grouping variable
#' specified in `list_of_columns`.
#'
#' @param theObject An object of class `ProteinQuantitativeData`.
#' @param list_of_columns A character vector of column names from the design matrix to use for coloring each RLE plot.
#' @param yaxis_limit A numeric vector of length 2 for the y-axis limits, applied to all plots.
#'
#' @return A named list of ggplot objects, where each name is a grouping variable.
#' @exportMethod plotRleList
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

#' @title Save a List of RLE Plots to Files
#' @description A convenience function to save a list of ggplot objects (like those
#' created by `plotRleList`) to files in multiple formats.
#'
#' @param input_list A named list of ggplot objects.
#' @param prefix A string prefix for the output filenames.
#' @param suffix A character vector of file extensions (e.g., `c("png", "pdf")`).
#' @param output_dir The directory to save the plots in.
#'
#' @return A data frame summarizing the saved files.
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

#' @title Plot PCA for Protein Data
#' @description Generates a PCA plot for the protein-level data in the object.
#'
#' @param theObject An object of class `ProteinQuantitativeData`.
#' @param grouping_variable The column name in the design matrix to use for coloring points.
#' @param shape_variable Optional. The column name for point shapes.
#' @param label_column Optional. The column name for text labels.
#' @param title The title for the plot.
#' @param font_size The font size for labels.
#'
#' @return A ggplot object representing the PCA plot.
#' @exportMethod plotPca
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


#' @title Plot a List of PCA Plots
#' @description Generates a list of PCA plots, one for each grouping variable
#' specified in `grouping_variables_list`.
#'
#' @param theObject An object of class `ProteinQuantitativeData`.
#' @param grouping_variables_list A character vector of column names from the design matrix to use for coloring.
#' @param label_column Optional. The column name for text labels.
#' @param title The base title for the plots.
#' @param font_size The font size for labels.
#'
#' @return A named list of ggplot objects.
#' @exportMethod plotPcaList
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


#' @title Save a List of PCA Plots to Files
#' @description A convenience function to save a list of ggplot objects (like those
#' created by `plotPcaList`) to files in multiple formats.
#'
#' @param input_list A named list of ggplot objects.
#' @param prefix A string prefix for the output filenames.
#' @param suffix A character vector of file extensions.
#' @param output_dir The directory to save the plots in.
#'
#' @return A data frame summarizing the saved files.
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


#' @title Get PCA Matrix
#' @description A method to perform PCA and return the matrix of principal components
#' along with the corresponding sample metadata.
#'
#' @param theObject An object of class `ProteinQuantitativeData`.
#'
#' @return A data frame containing the principal components (PC1, PC2, etc.) and
#'   the joined metadata from the design matrix.
#' @exportMethod getPcaMatrix
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


#' @title Correlate Technical Replicates for Protein Data
#' @description A method to calculate the correlation between technical replicates
#' at the protein level.
#'
#' @param theObject An object of class `ProteinQuantitativeData`.
#' @param tech_rep_num_column The column name in the design matrix that identifies the replicate number.
#' @param tech_rep_remove_regex A regex string to identify samples to exclude (e.g., pools).
#'
#' @return A data frame with Pearson and Spearman correlations for each feature.
#' @exportMethod proteinTechRepCorrelation
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
#' @title Calculate Pearson Correlation for All Sample Pairs
#' @description This method calculates the Pearson correlation coefficient for all
#' pairs of samples within each group defined by `correlation_group`.
#'
#' @param theObject An object of class `ProteinQuantitativeData`.
#' @param tech_rep_remove_regex A regex string to identify samples to exclude.
#' @param correlation_group The column name in the design matrix to use for grouping samples.
#'
#' @return A data frame containing the Pearson correlation for each pair of samples.
#' @exportMethod pearsonCorForSamplePairs
setMethod(f="plotPearson",
          signature="ProteinQuantitativeData",
          definition=function(theObject, tech_rep_remove_regex = "pool", correlation_group = NA) {

            correlation_group_to_use <- correlation_group

            if( is.na( correlation_group)) {
              correlation_group_to_use <- theObject@technical_replicate_id
            }

            correlation_vec <- pearsonCorForSamplePairs(theObject
                                                        , tech_rep_remove_regex
                                                        , correlation_group = correlation_group_to_use)

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
#' @title A Class to Store Grid Plot Data
#' @description The `GridPlotData` class is a container for storing various types of
#' quality control (QC) plots and their corresponding titles, designed to facilitate
#' the creation of a composite QC figure.
#'
#' @slot pca_plots A list of PCA plots (ggplot objects).
#' @slot density_plots A list of density plots (ggplot objects).
#' @slot rle_plots A list of RLE plots (ggplot objects).
#' @slot pearson_plots A list of Pearson correlation plots (ggplot objects).
#' @slot pca_titles A list of titles for the PCA plots.
#' @slot density_titles A list of titles for the density plots.
#' @slot rle_titles A list of titles for the RLE plots.
#' @slot pearson_titles A list of titles for the Pearson correlation plots.
#'
#' @return An object of class `GridPlotData`.
#' @export
setClass("GridPlotData",
         slots = list(
           pca_plots = "list",
           density_plots = "list",
           rle_plots = "list",
           pearson_plots = "list",
           pca_titles = "list",
           density_titles = "list",
           rle_titles = "list",
           pearson_titles = "list"
         ))

#' @title Initialize an Empty GridPlotData Object
#' @description This function creates and initializes an empty `GridPlotData` object,
#' ready to be populated with QC plots.
#'
#' @param dummy A placeholder argument, not used.
#'
#' @return An empty object of class `GridPlotData`.
#' @export
setGeneric("InitialiseGrid", function(dummy = NULL) {
  standardGeneric("InitialiseGrid")
})

#' @rdname InitialiseGrid
#' @export
setMethod("InitialiseGrid", 
          signature(dummy = "ANY"),
          function(dummy = NULL) {
            new("GridPlotData",
                pca_plots = list(),
                density_plots = list(),
                rle_plots = list(),
                pearson_plots = list(),
                pca_titles = list(),
                density_titles = list(),
                rle_titles = list(),
                pearson_titles = list())
          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Create a QC composite figure

#' @title Create a Composite Quality Control (QC) Grid Plot
#' @description This method assembles a grid of various QC plots (PCA, density, RLE,
#' Pearson correlation) stored in a `GridPlotData` object into a single, composite figure.
#' It can also save the figure to a file.
#'
#' @param theObject An object of class `GridPlotData`.
#' @param pca_titles A list of titles for the PCA plots. If `NULL`, titles from the object are used.
#' @param density_titles A list of titles for the density plots. If `NULL`, titles from the object are used.
#' @param rle_titles A list of titles for the RLE plots. If `NULL`, titles from the object are used.
#' @param pearson_titles A list of titles for the Pearson correlation plots. If `NULL`, titles from the object are used.
#' @param save_path Optional. The directory path where the plot should be saved.
#' @param file_name The base name for the saved file (without extension).
#'
#' @return A composite ggplot object (created with `patchwork`).
#' @export
setGeneric(name = "createGridQC",
           def = function(theObject, pca_titles, density_titles, rle_titles, pearson_titles, save_path = NULL, file_name = "pca_density_rle_pearson_corr_plots_merged") {
             standardGeneric("createGridQC")
           },
           signature = c("theObject", "pca_titles", "density_titles", "rle_titles", "pearson_titles", "save_path", "file_name"))

#' @rdname createGridQC
#' @export
setMethod(f = "createGridQC",
          signature = "GridPlotData",
          definition = function(theObject, pca_titles = NULL, density_titles = NULL, rle_titles = NULL, pearson_titles = NULL, save_path = NULL, file_name = "pca_density_rle_pearson_corr_plots_merged") {
            
            # Use stored titles if not provided as parameters
            pca_titles <- if(is.null(pca_titles)) theObject@pca_titles else pca_titles
            density_titles <- if(is.null(density_titles)) theObject@density_titles else density_titles
            rle_titles <- if(is.null(rle_titles)) theObject@rle_titles else rle_titles
            pearson_titles <- if(is.null(pearson_titles)) theObject@pearson_titles else pearson_titles
            
            createLabelPlot <- function(title) {
              # Option 1: Use xlim to expand the plot area and position text at left edge
              ggplot() + 
                annotate("text", x = 0, y = 0.5, label = title, size = 5, hjust = 0) +
                xlim(0, 1) +  # Explicitly set the x limits
                theme_void() +
                theme(
                  plot.margin = margin(5, 5, 5, 5),
                  panel.background = element_blank()
                )
            }
            
            # Create basic plots without titles
            createPcaPlot <- function(plot) {
              plot +
                xlim(-40, 45) + ylim(-30, 25) +
                theme(text = element_text(size = 15),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank())
            }
            
            createDensityPlot <- function(plot) {
              # For all plots, just apply the theme without adding title
              if (inherits(plot, "patchwork")) {
                plot & 
                  theme(
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    text = element_text(size = 15)
                  )
              } else {
                plot +
                  theme(text = element_text(size = 15),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank())
              }
            }
            
            createRlePlot <- function(plot) {
              plot +
                theme(text = element_text(size = 15),
                      axis.text.x = element_blank(),
                      axis.ticks.x = element_blank())
            }
            
            createPearsonPlot <- function(plot) {
              plot +
                theme(text = element_text(size = 15))
            }
            
            # Create plots without titles
            created_pca_plots <- lapply(theObject@pca_plots, createPcaPlot)
            created_density_plots <- lapply(theObject@density_plots, createDensityPlot)
            created_rle_plots <- lapply(theObject@rle_plots, createRlePlot)
            created_pearson_plots <- lapply(theObject@pearson_plots, createPearsonPlot)
            
            # Create label plots
            pca_labels <- lapply(pca_titles, createLabelPlot)
            density_labels <- lapply(density_titles, createLabelPlot)
            rle_labels <- lapply(rle_titles, createLabelPlot)
            pearson_labels <- lapply(pearson_titles, createLabelPlot)
            
            # Combine with labels above each row - modified to keep legends with their plots
            combined_plot <- (
              wrap_plots(pca_labels, ncol = 3) /
              wrap_plots(created_pca_plots, ncol = 3) /
              wrap_plots(density_labels, ncol = 3) /
              wrap_plots(created_density_plots, ncol = 3) /
              wrap_plots(rle_labels, ncol = 3) /
              wrap_plots(created_rle_plots, ncol = 3) /
              wrap_plots(pearson_labels, ncol = 3) /
              wrap_plots(created_pearson_plots, ncol = 3)
            ) +
              plot_layout(heights = c(0.1, 1, 0.1, 1, 0.1, 1, 0.1, 1))

            if (!is.null(save_path)) {
              sapply(c("png", "pdf", "svg"), function(ext) {
                ggsave(
                  plot = combined_plot,
                  filename = file.path(save_path, paste0(file_name, ".", ext)),
                  width = 14,
                  height = 16 # Increased height to accommodate label rows
                )
              })
              message(paste("Plots saved in", save_path))
            }
            
            return(combined_plot)
          })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
## normalise between Arrays
#' @title Normalize Protein Abundances Between Samples
#' @description This method applies a specified normalization technique to the protein
#' quantification data to correct for systematic variations between samples.
#'
#' @param theObject An object of class `ProteinQuantitativeData`.
#' @param normalisation_method The normalization method to use. Supported options are
#'   `"cyclicloess"`, `"quantile"`, `"scale"`, and `"none"`.
#'
#' @return The `ProteinQuantitativeData` object with normalized protein abundances.
#' @export
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

#' @title Calculate Pearson Correlation for All Sample Pairs
#' @description This method calculates the Pearson correlation coefficient for all
#' pairs of samples. It can group samples before calculating correlations.
#'
#' @param theObject An object of class `ProteinQuantitativeData`.
#' @param tech_rep_remove_regex A regex string to identify samples to exclude from the analysis.
#' @param correlation_group The column name in the design matrix used to group samples.
#'   Correlations are calculated for all pairs within each group. If `NA`, the
#'   technical replicate column is used.
#'
#' @return A data frame containing the Pearson correlation for each pair of samples.
#' @export
setMethod(f="pearsonCorForSamplePairs"
          , signature="ProteinQuantitativeData"
          , definition=function( theObject, tech_rep_remove_regex = NULL, correlation_group = NA ) {
            protein_quant_table <- theObject@protein_quant_table
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id

            replicate_group_column <- theObject@technical_replicate_id
            if(!is.na( correlation_group )) {
              replicate_group_column <- correlation_group
            }

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

#' @title Get Negative Control Proteins using ANOVA
#' @description A method to identify a set of negative control proteins based on an
#' ANOVA test across experimental groups.
#'
#' @param theObject An object of class `ProteinQuantitativeData`.
#' @param ruv_grouping_variable The column name in the design matrix for the grouping variable.
#' @param percentage_as_neg_ctrl The percentage of total proteins to select as controls.
#' @param num_neg_ctrl The absolute number of control proteins to select.
#' @param ruv_qval_cutoff The q-value cutoff for considering proteins as non-significant.
#' @param ruv_fdr_method The FDR method to use (`"qvalue"` or `"BH"`).
#'
#' @return A logical vector indicating which proteins are selected as controls.
#' @export
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

            control_genes_index <- getNegCtrlProtAnovaHelper( normalised_frozen_protein_matrix_filt[,design_matrix |> dplyr::pull(!!sym(sample_id)) ]
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

#' @title Get Proteins with Low Coefficient of Variation
#' @description A method to identify a set of negative control proteins by selecting
#' those with the lowest coefficient of variation (CV) across all samples.
#'
#' @param theObject An object of class `ProteinQuantitativeData`.
#' @param percentage_as_neg_ctrl The percentage of total proteins to select.
#' @param num_neg_ctrl The absolute number of control proteins to select.
#'
#' @return A logical vector indicating which proteins are selected as controls.
#' @export
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
#' @title Plot RUV Canonical Correlation Analysis
#' @description This method performs a canonical correlation analysis as part of the RUV
#' (Remove Unwanted Variation) workflow and generates a plot of the results. It helps
#' to visualize the amount of variation explained by the factors of interest.
#'
#' @param theObject An object of class `ProteinQuantitativeData`.
#' @param ctrl A logical vector indicating the negative control proteins.
#' @param num_components_to_impute The number of principal components to use for NIPALS
#'   imputation if there are missing values.
#' @param ruv_grouping_variable The column in the design matrix that defines the groups
#'   for the RUV analysis.
#'
#' @return A ggplot object showing the canonical correlation plot.
#' @export
setMethod( f = "ruvCancor"
           , signature="ProteinQuantitativeData"
           , definition=function( theObject, ctrl= NULL, num_components_to_impute=NULL, ruv_grouping_variable = NULL) {
             protein_quant_table <- theObject@protein_quant_table
             protein_id_column <- theObject@protein_id_column
             design_matrix <- theObject@design_matrix
             group_id <- theObject@group_id
             sample_id <- theObject@sample_id

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

             Y <-  t( normalised_frozen_protein_matrix_filt[,design_matrix |> dplyr::pull(!!sym(sample_id))])
             if( length(which( is.na(normalised_frozen_protein_matrix_filt) )) > 0 ) {
               Y <- impute.nipals( t( normalised_frozen_protein_matrix_filt[,design_matrix |> dplyr::pull(!!sym(sample_id))])
                                   , ncomp=num_components_to_impute)
             }

             cancorplot_r2 <- ruv_cancorplot( Y ,
                                              X = design_matrix |>
                                                dplyr::pull(!!sym(ruv_grouping_variable)),
                                              ctl = ctrl)
             cancorplot_r2


           })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#' @title Get RUV-III Replicate Matrix
#' @description Creates a replicate matrix required for the RUV-III algorithm. This
#' matrix indicates which samples are technical replicates of each other.
#'
#' @param theObject An object of class `ProteinQuantitativeData`.
#' @param ruv_grouping_variable The column name in the design matrix that identifies
#'   technical replicates or groups of samples that should be similar.
#'
#' @return A matrix suitable for use with the RUV-III algorithm.
#' @export
setGeneric(name="getRuvIIIReplicateMatrix"
           , def=function( theObject,  ruv_grouping_variable = NULL) {
             standardGeneric("getRuvIIIReplicateMatrix")
           }
           , signature=c("theObject", "ruv_grouping_variable"))

#' @rdname getRuvIIIReplicateMatrix
#' @export
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
#' @title Apply RUV-III Normalization with Varying Controls
#' @description This method applies the RUV-III normalization algorithm, which is robust
#' to the choice of negative controls. It uses a replicate matrix to estimate and
#' remove unwanted variation.
#'
#' @param theObject An object of class `ProteinQuantitativeData`.
#' @param ruv_grouping_variable The column name in the design matrix that identifies replicates.
#' @param ruv_number_k The number of unwanted factors to estimate and remove (k).
#' @param ctrl A logical vector indicating the negative control proteins.
#'
#' @return A `ProteinQuantitativeData` object with RUV-III normalized data.
#' @export
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

             Y <-  t( normalised_frozen_protein_matrix_filt[,design_matrix |> dplyr::pull(!!sym(sample_id))])

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
#' @title Filter Proteins Based on Missing Value Percentage
#' @description Filters proteins based on the percentage of missing values within groups.
#' A protein is kept if it has a sufficient number of observations in a sufficient
#' number of groups.
#'
#' @param theObject An object of class `ProteinQuantitativeData`.
#' @param ruv_grouping_variable The column name in the design matrix for grouping samples.
#' @param groupwise_percentage_cutoff The minimum percentage of non-missing values required
#'   within a group for a protein to be considered "present" in that group.
#' @param max_groups_percentage_cutoff The minimum percentage of groups in which a protein
#'   must be "present" for it to be kept.
#' @param proteins_intensity_cutoff_percentile An intensity percentile used to define a
#'   minimum threshold. Values below this are treated as missing during filtering.
#'
#' @return A `ProteinQuantitativeData` object with filtered protein data.
#' @export
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

#' @rdname removeRowsWithMissingValuesPercent
#' @export
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
                                                                           , cols= protein_id_column
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

#' @title Average Technical Replicates
#' @description Averages the abundance values across technical replicates for each protein.
#' This collapses replicate columns into a single column representing the biological sample.
#'
#' @param theObject An object of class `ProteinQuantitativeData`.
#' @param design_matrix_columns A character vector of additional columns to retain in the
#'   updated design matrix.
#'
#' @return A `ProteinQuantitativeData` object with technical replicates averaged.
#' @export
setGeneric(name="averageTechReps"
           , def=function( theObject, design_matrix_columns ) {
             standardGeneric("averageTechReps")
           }
           , signature=c("theObject", "design_matrix_columns" ))

#' @rdname averageTechReps
#' @export
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


#' @title Preserve NA Values from Peptide Level
#' @description After protein-level data has been processed (e.g., imputed), this
#' method re-introduces NA values that were originally present at the peptide level.
#' If all peptides for a given protein in a given sample were NA, the corresponding
#' protein value is set to NA.
#'
#' @param peptide_obj An object of class `PeptideQuantitativeData` containing the
#'   original peptide-level data.
#' @param protein_obj An object of class `ProteinQuantitativeData` containing the
#'   protein-level data to be modified.
#'
#' @return The modified `ProteinQuantitativeData` object.
#' @export
setGeneric(name="preservePeptideNaValues"
           , def=function( peptide_obj, protein_obj)  {
             standardGeneric("preservePeptideNaValues")
           }
           , signature=c("peptide_obj", "protein_obj" ))

#' @rdname preservePeptideNaValues
#' @export
setMethod( f = "preservePeptideNaValues"
           , signature=c( "PeptideQuantitativeData", "ProteinQuantitativeData" )
           , definition= function( peptide_obj, protein_obj) {
             preservePeptideNaValuesHelper( peptide_obj, protein_obj)
           })

#' @title Helper to Preserve NA Values from Peptide Level
#' @description This helper function contains the logic for `preservePeptideNaValues`.
#' It checks where NAs existed at the peptide level and applies them back to the
#' protein level data.
#'
#' @param peptide_obj A `PeptideQuantitativeData` object.
#' @param protein_obj A `ProteinQuantitativeData` object.
#'
#' @return The modified `ProteinQuantitativeData` object.
#' @keywords internal
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


#' @title Choose the Best Protein Accession from Protein Groups
#' @description This method resolves protein groups (represented by delimited strings
#' of accessions) by selecting the best protein accession based on available
#' annotation, typically from a FASTA file. It then aggregates the quantitative
#' data for the selected accessions.
#'
#' @param theObject An object of class `ProteinQuantitativeData`.
#' @param delim The delimiter used to split protein accession strings.
#' @param seqinr_obj A sequence database object (e.g., from `seqinr::read.fasta`)
#'   containing protein annotations.
#' @param seqinr_accession_column The column in `seqinr_obj` that contains the
#'   accession numbers to match against.
#' @param replace_zero_with_na A boolean indicating whether to replace zero values
#'   with `NA` after aggregation.
#' @param aggregation_method The method to use for aggregating data for proteins
#'   that map to the same chosen accession. One of `"sum"`, `"mean"`, or `"median"`.
#'
#' @return A `ProteinQuantitativeData` object with resolved protein accessions.
#' @export
setGeneric(name="chooseBestProteinAccession"
           , def=function(theObject, delim=NULL, seqinr_obj=NULL, seqinr_accession_column=NULL, replace_zero_with_na = NULL, aggregation_method = NULL) {
             standardGeneric("chooseBestProteinAccession")
           }
           , signature=c("theObject", "delim", "seqinr_obj", "seqinr_accession_column"))

#' @rdname chooseBestProteinAccession
#' @export
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


#' @title Choose Best Protein Accession and Sum Duplicates
#' @description A simpler method to resolve protein groups. It takes the first
#' accession from a delimited string and then sums the abundances for any
#' resulting duplicate protein accessions.
#'
#' @param theObject An object of class `ProteinQuantitativeData`.
#' @param delim The delimiter for splitting protein accession strings.
#' @param quant_columns_pattern A regex pattern to identify the quantitative (sample) columns.
#' @param islogged A boolean indicating if the data is log-transformed. If `TRUE`,
#'   abundances are anti-logged, summed, and then log-transformed back.
#'
#' @return A `ProteinQuantitativeData` object with resolved and summed protein data.
#' @export
setGeneric(name="chooseBestProteinAccessionSumDuplicates"
           , def=function( theObject, delim, quant_columns_pattern, islogged ) {
             standardGeneric("chooseBestProteinAccessionSumDuplicates")
           }
           , signature=c("theObject", "delim", "quant_columns_pattern", "islogged" ))

#' @rdname chooseBestProteinAccessionSumDuplicates
#' @export
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

#' @title Filter Samples by Correlation Threshold
#' @description This method removes samples that have a low correlation with other
#' samples, based on a provided correlation matrix and a minimum threshold.
#'
#' @param theObject An object of class `ProteinQuantitativeData`.
#' @param pearson_correlation_per_pair A data frame containing pairwise sample correlations,
#'   typically from `pearsonCorForSamplePairs`.
#' @param min_pearson_correlation_threshold The minimum Pearson correlation required for a
#'   sample to be retained.
#'
#' @return A `ProteinQuantitativeData` object with low-quality samples removed.
#' @export
setGeneric(name="filterSamplesByProteinCorrelationThreshold"
           , def=function( theObject, pearson_correlation_per_pair = NULL, min_pearson_correlation_threshold = NULL ) {
             standardGeneric("filterSamplesByProteinCorrelationThreshold")
           }
           , signature=c("theObject", "pearson_correlation_per_pair", "min_pearson_correlation_threshold" ))

#' @rdname filterSamplesByProteinCorrelationThreshold
#' @export
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
#' @title Compare Two ProteinData Objects
#' @description Compares two `ProteinQuantitativeData` objects and reports the differences
#' in their protein and sample compositions.
#'
#' @param object_a The first `ProteinQuantitativeData` object.
#' @param object_b The second `ProteinQuantitativeData` object.
#'
#' @return A tibble summarizing the number of proteins and samples that are unique
#'   to each object or shared between them.
#' @export
compareTwoProteinDataObjects <- function( object_a, object_b) {


  object_a_proteins <- object_a@protein_quant_table |>
    distinct(!!sym(object_a@protein_id_column)) |>
    dplyr::pull(!!sym(object_a@protein_id_column))

  object_b_proteins <- object_b@protein_quant_table |>
    distinct(!!sym(object_b@protein_id_column)) |>
    dplyr::pull(!!sym(object_b@protein_id_column))

  object_a_samples <- object_a@design_matrix |>
    distinct(!!sym(object_a@sample_id)) |>
    dplyr::pull(!!sym(object_a@sample_id))

  object_b_samples <- object_b@design_matrix |>
    distinct(!!sym(object_b@sample_id)) |>
    dplyr::pull(!!sym(object_b@sample_id))


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

#' @title Summarize a ProteinData Object
#' @description Provides a simple summary of a `ProteinQuantitativeData` object,
#' reporting the total number of proteins and samples.
#'
#' @param theObject An object of class `ProteinQuantitativeData`.
#'
#' @return A list containing the number of proteins (`num_proteins`) and the
#'   number of samples (`num_samples`).
#' @export
summariseProteinObject <- function ( theObject) {
  num_proteins <- theObject@protein_quant_table |>
    distinct(!!sym(theObject@protein_id_column)) |>
    dplyr::pull(!!sym(theObject@protein_id_column))

  num_samples <- theObject@design_matrix |>
    distinct(!!sym(theObject@sample_id)) |>
    dplyr::pull(!!sym(theObject@sample_id))

  summary_list <- list( num_proteins = length(num_proteins)
       , num_samples = length(num_samples))

  summary_list

}


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @title Plot Density of Principal Components
#' @description A method to create density plots (as boxplots) for the first two
#' principal components from a PCA result stored in a ggplot object. This is a
#' wrapper for objects of class `gg`.
#'
#' @param theObject A ggplot object, expected to contain PCA data.
#' @param grouping_variable The column name in the plot data to use for grouping.
#' @param title The title for the plot.
#' @param font_size The font size for plot text.
#'
#' @return A combined ggplot object (using `patchwork`) showing boxplots for PC1 and PC2.
#' @export
setMethod(f="plotDensity"
          , signature="gg"
          , definition=function(theObject, grouping_variable, title = "", font_size = 8) {
            # For gg class objects, create a copy and change its class to ggplot
            gg_obj <- theObject
            class(gg_obj) <- "ggplot"
            
            # Then call the ggplot method
            plotDensity(gg_obj, grouping_variable, title, font_size)
          })

#' @rdname plotDensity
#' @export
setMethod(f="plotDensity"
          , signature="ggplot"
          , definition=function(theObject, grouping_variable, title = "", font_size = 8) {
            # First try to get data directly from the ggplot object's data element
            if (!is.null(theObject$data) && is.data.frame(theObject$data)) {
              pca_data <- as_tibble(theObject$data)
            } else {
              # Fall back to other extraction methods
              pca_data <- as_tibble(ggplot_build(theObject)$data[[1]])
              
              # If the data doesn't have PC1/PC2, try to extract from the plot's environment
              if (!("PC1" %in% colnames(pca_data) && "PC2" %in% colnames(pca_data))) {
                # Try to get the data from the plot's environment
                if (exists("data", envir = environment(theObject$mapping$x))) {
                  pca_data <- as_tibble(get("data", envir = environment(theObject$mapping$x)))
                } else {
                  stop("Could not extract PCA data from the ggplot object")
                }
              }
            }
            
            # Check if grouping variable exists in the data
            if (!grouping_variable %in% colnames(pca_data)) {
              stop(sprintf("grouping_variable '%s' not found in the data", grouping_variable))
            }
            
            # Create PC1 boxplot
            pc1_box <- ggplot(pca_data, aes(x = !!sym(grouping_variable), y = PC1, fill = !!sym(grouping_variable))) +
              geom_boxplot(notch = TRUE) +
              theme_bw() +
              labs(title = title,
                   x = "",
                   y = "PC1") +
              theme(
                legend.position = "none",
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                text = element_text(size = font_size),
                plot.margin = margin(b = 0, t = 5, l = 5, r = 5),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank()
              )
            
            # Create PC2 boxplot
            pc2_box <- ggplot(pca_data, aes(x = !!sym(grouping_variable), y = PC2, fill = !!sym(grouping_variable))) +
              geom_boxplot(notch = TRUE) +
              theme_bw() +
              labs(x = "",
                   y = "PC2") +
              theme(
                legend.position = "none",
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                text = element_text(size = font_size),
                plot.margin = margin(t = 0, b = 5, l = 5, r = 5),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank()
              )
            
            # Combine plots with minimal spacing
            combined_plot <- pc1_box / pc2_box + 
              plot_layout(heights = c(1, 1)) +
              plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))
            
            return(combined_plot)
          }) 

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @title Plot a List of Density Plots
#' @description Generates a list of density plots for principal components, creating
#' one plot for each grouping variable specified.
#'
#' @param theObject An object of class `ProteinQuantitativeData`.
#' @param grouping_variables_list A character vector of column names from the design
#'   matrix to use for grouping in each plot.
#' @param title The base title for the plots.
#' @param font_size The font size for plot text.
#'
#' @return A named list of combined ggplot objects (from `plotDensity`).
#' @export
setMethod(f="plotDensityList"
          , signature="ProteinQuantitativeData"
          , definition=function(theObject, grouping_variables_list, title = "", font_size = 8) {
            
            # Create a list of density plots for each grouping variable
            density_plots_list <- purrr::map(grouping_variables_list, function(group_var) {
              tryCatch({
                plotDensity(theObject, 
                           grouping_variable = group_var,
                           title = title,
                           font_size = font_size)
              }, error = function(e) {
                warning(sprintf("Error creating density plot for %s: %s", group_var, e$message))
              return(NULL)
                  })
                              })
            
            # Name the list elements with the grouping variables
            names(density_plots_list) <- grouping_variables_list
            
            # Remove any NULL elements (failed plots)
            density_plots_list <- density_plots_list[!sapply(density_plots_list, is.null)]
            
            return(density_plots_list)
          })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @title Save a List of Density Plots to Files
#' @description A convenience function to save a list of density plots (as generated
#' by `plotDensityList`) to files in various formats.
#'
#' @param input_list A named list of ggplot objects.
#' @param prefix A string prefix for the output filenames.
#' @param suffix A character vector of file extensions (e.g., `c("png", "pdf")`).
#' @param output_dir The directory where the plots will be saved.
#'
#' @return A data frame summarizing the saved files.
#' @export
savePlotDensityList <- function(input_list, prefix = "Density", suffix = c("png", "pdf"), output_dir) {
  
  list_of_filenames <- expand_grid(column = names(input_list), suffix = suffix) |>
    mutate(filename = paste0(prefix, "_", column, ".", suffix)) |>
    left_join(tibble(column = names(input_list),
              plots = input_list),
              by = join_by(column))
  
  purrr::walk2(list_of_filenames$plots,
               list_of_filenames$filename,
               \(.x, .y) {
                 ggsave(plot = .x, filename = file.path(output_dir, .y))
               })
  
  list_of_filenames
}

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
