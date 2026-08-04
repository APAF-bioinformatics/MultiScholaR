#' ProteinQuantitativeData S4 Class
#'
#' @description
#' An S4 class to store and manage protein-level quantitative data
#' from proteomics experiments along with experimental design metadata.
#'
#' @slot protein_quant_table A data.frame containing protein quantification data
#' @slot protein_id_column Character string naming the protein ID column
#' @slot design_matrix A data.frame containing experimental design
#' @slot protein_id_table A data.frame containing protein ID information
#' @slot sample_id Character string naming the sample ID column
#' @slot group_id Character string naming the group column
#' @slot technical_replicate_id Character string naming the replicate column
#' @slot args A list of additional arguments
#'
#' @importFrom purrr map map2 walk walk2 map_chr
#' @importFrom stringr str_extract str_detect
#' @importFrom dplyr filter select mutate distinct left_join arrange pull rename
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom viridis viridis plasma inferno
#' @importFrom ruv ruv_cancorplot
#' @importFrom patchwork wrap_plots
#' @importFrom vroom vroom_write
#' @importFrom writexl write_xlsx
#' @importFrom logger log_debug log_warn
#' @exportClass ProteinQuantitativeData
ProteinQuantitativeData <- setClass("ProteinQuantitativeData",
  slots = c(
    # Protein vs Sample quantitative data
    protein_quant_table = "data.frame",
    protein_id_column = "character"


    # Design Matrix Information

    , design_matrix = "data.frame",
    protein_id_table = "data.frame",
    sample_id = "character",
    group_id = "character",
    technical_replicate_id = "character",
    args = "list"
  ),
  prototype = list(
    # Protein vs Sample quantitative data
    protein_id_column = "Protein.Ids",
    protein_id_table = data.frame()

    # Design Matrix Information

    , sample_id = "Sample_id",
    group_id = "group",
    technical_replicate_id = "replicates",
    args = list()
  ),
  validity = function(object) {
    if (!is.data.frame(object@protein_quant_table)) {
      stop("protein_quant_table must be a data.frame")
    }

    if (!is.character(object@protein_id_column)) {
      stop("protein_id_column must be a character")
    }

    if (!is.data.frame(object@design_matrix)) {
      stop("design_matrix must be a data.frame")
    }

    if (!is.character(object@sample_id)) {
      stop("sample_id must be a character")
    }

    if (!is.character(object@group_id)) {
      stop("group_id must be a character")
    }

    if (!is.character(object@technical_replicate_id)) {
      stop("technical_replicate_id must be a character")
    }


    if (!object@protein_id_column %in% colnames(object@protein_quant_table)) {
      stop("Protein ID column must be in the protein data table")
    }

    if (ncol(object@protein_id_table) > 0L &&
        !object@protein_id_column %in% colnames(object@protein_id_table)) {
      stop("Protein ID column must be in the protein ID/provenance table")
    }


    if (!object@sample_id %in% colnames(object@design_matrix)) {
      stop("Sample ID column must be in the design matrix")
    }


    # Need to check the rows names in design matrix and the column names of the data table

    samples_in_protein_quant_table <- setdiff(colnames(object@protein_quant_table), object@protein_id_column)

    samples_in_design_matrix <- object@design_matrix |> dplyr::pull(!!sym(object@sample_id))


    if (!setequal(
      as.character(samples_in_protein_quant_table),
      as.character(samples_in_design_matrix)
    )) {
      stop("Samples in protein data and design matrix must be the same")
    }
  }
)

.proteinIdTableFromPeptideLineage <- function(peptideS4, quantifiedProteinIds) {
  proteinColumn <- peptideS4@protein_id_column
  quantifiedIds <- data.frame(
    value = as.character(quantifiedProteinIds),
    stringsAsFactors = FALSE
  )
  names(quantifiedIds) <- proteinColumn
  quantifiedIds <- quantifiedIds[!duplicated(quantifiedIds), , drop = FALSE]

  provenance <- NULL
  if (is.list(peptideS4@args)) {
    provenance <- peptideS4@args$protein_accession_provenance
  }
  requiredProvenance <- unique(c(proteinColumn, "Protein.Ids"))
  if (!is.data.frame(provenance) ||
      !all(requiredProvenance %in% names(provenance))) {
    availableProvenance <- requiredProvenance[
      requiredProvenance %in% names(peptideS4@peptide_data)
    ]
    provenance <- peptideS4@peptide_data[, availableProvenance, drop = FALSE]
  } else {
    provenance <- provenance[, requiredProvenance, drop = FALSE]
  }

  if (!proteinColumn %in% names(provenance)) {
    return(quantifiedIds)
  }
  for (column in names(provenance)) {
    provenance[[column]] <- as.character(provenance[[column]])
  }
  provenance <- provenance[!duplicated(provenance), , drop = FALSE]
  dplyr::left_join(quantifiedIds, provenance, by = proteinColumn)
}

setMethod(
  "initialize", "ProteinQuantitativeData",
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

#' @export
setMethod(
  f = "setProteinData",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, protein_quant_table, protein_id_column) {
    theObject@protein_quant_table <- protein_quant_table
    theObject@protein_id_column <- protein_id_column

    return(theObject)
  }
)
