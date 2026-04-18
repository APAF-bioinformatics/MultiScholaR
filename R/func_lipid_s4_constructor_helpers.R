#' Create LipidomicsAssayData Object
#'
#' Constructor function for the LipidomicsAssayData class.
#'
#' @param lipid_data Named list of data frames (assays).
#' @param design_matrix Experimental design data frame.
#' @param lipid_id_column Name of the **primary feature ID** column within assays (e.g., `"database_identifier"`).
#' @param annotation_id_column Name of the **annotation ID** column (e.g., `"lipid_identification"`).
#' @param sample_id Name of the sample ID column in design_matrix and assays.
#' @param group_id Name of the group column in design_matrix.
#' @param technical_replicate_id Name of the technical replicate column in design_matrix (use NA_character_ if none).
#' @param database_identifier_type Type of identifier in the `annotation_id_column` (e.g., `"Mixed_CHEBI_Unknown"`).
#' @param internal_standard_regex Regex to identify internal standards. Use `NA_character_` or `""` if none.
#' @param args List of arguments (e.g., from config).
#'
#' @return A LipidomicsAssayData object.
#' @export
#' @examples
#' \dontrun{
#' # Assuming lcms_pos_df, lcms_neg_df, gcms_df are data frames
#' # with 'Lipid' as ID column and samples as other columns
#' # Assuming design_df has 'SampleID', 'Group', 'Replicate' columns
#' assays_list <- list(
#'     LCMS_Pos = lcms_pos_df,
#'     LCMS_Neg = lcms_neg_df,
#'     GCMS = gcms_df
#' )
#' config <- list(...) # Your config list
#'
#' met_assay_obj <- createLipidomicsAssayData(
#'     lipid_data = assays_list,
#'     design_matrix = design_df,
#'     lipid_id_column = "Lipid",
#'     sample_id = "SampleID",
#'     group_id = "Group",
#'     technical_replicate_id = "Replicate",
#'     database_identifier_type = "InternalName",
#'     internal_standard_regex = "^IS_",
#'     args = config
#' )
#' }
createLipidomicsAssayData <- function(
  lipid_data,
  design_matrix,
  lipid_id_column = "database_identifier",
  annotation_id_column = "lipid_identification",
  sample_id = "Sample_ID",
  group_id = "group",
  technical_replicate_id = NA_character_,
  database_identifier_type = "Unknown",
  internal_standard_regex = NA_character_,
  args = list()
) {
    # Perform basic checks before creating the object
    stopifnot(is.list(lipid_data))
    stopifnot(all(sapply(lipid_data, is.data.frame)))
    stopifnot(is.data.frame(design_matrix))
    # Add more checks as needed...

    obj <- new("LipidomicsAssayData",
        lipid_data = lipid_data,
        lipid_id_column = lipid_id_column,
        annotation_id_column = annotation_id_column,
        database_identifier_type = database_identifier_type,
        internal_standard_regex = internal_standard_regex,
        design_matrix = design_matrix,
        sample_id = sample_id,
        group_id = group_id,
        technical_replicate_id = technical_replicate_id,
        args = args
    )
    # Validity check is automatically called by 'new'
    return(obj)
}

