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
                                       samples_in_peptide_data <- unique(object@peptide_data[[object@sample_id]])
                                       samples_in_design_matrix <- object@design_matrix[[object@sample_id]]

                                       invisible(NULL)

                                     }

)

#' @export
PeptideQuantitativeDataDiann <- function( peptide_data
                                          , design_matrix
                                          , sample_id = "Run"
                                          , group_id = "group"
                                          , technical_replicate_id = "replicates"
                                          , args = NA) {

  protein_id_column <- if ("Protein.Group" %in% names(peptide_data)) {
    "Protein.Group"
  } else if ("Protein.Ids" %in% names(peptide_data)) {
    "Protein.Ids"
  } else {
    stop(
      "PeptideQuantitativeDataDiann: Protein.Group or Protein.Ids is required.",
      call. = FALSE
    )
  }

  peptide_data <- new( "PeptideQuantitativeData"

    # Protein vs Sample quantitative data
    , peptide_data = peptide_data
    , protein_id_column = protein_id_column
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

# Peptide matrices need a compact row key, but that key must never become the
# source of truth for the biological identity.  Keep the reversible relation in
# @args so existing serialized PeptideQuantitativeData objects remain readable
# without changing the S4 schema.
.buildPeptideFeatureKeyMap <- function(theObject) {
  protein_column <- theObject@protein_id_column
  peptide_column <- theObject@peptide_sequence_column
  peptide_data <- theObject@peptide_data

  required_columns <- c(protein_column, peptide_column)
  missing_columns <- setdiff(required_columns, names(peptide_data))
  if (length(missing_columns) > 0L) {
    stop(
      "Cannot build peptide feature keys; missing column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  feature_pairs <- data.frame(
    active_protein_id = as.character(peptide_data[[protein_column]]),
    peptide_sequence = as.character(peptide_data[[peptide_column]]),
    stringsAsFactors = FALSE
  )
  feature_pairs <- feature_pairs[!duplicated(feature_pairs), , drop = FALSE]

  invalid_identity <- is.na(feature_pairs$active_protein_id) |
    !nzchar(feature_pairs$active_protein_id) |
    is.na(feature_pairs$peptide_sequence) |
    !nzchar(feature_pairs$peptide_sequence)
  if (any(invalid_identity)) {
    stop(
      "Peptide feature identities cannot contain missing or empty protein/peptide values.",
      call. = FALSE
    )
  }

  contains_separator <- grepl("%", feature_pairs$active_protein_id, fixed = TRUE) |
    grepl("%", feature_pairs$peptide_sequence, fixed = TRUE)
  display_keys <- paste(
    feature_pairs$active_protein_id,
    feature_pairs$peptide_sequence,
    sep = "%"
  )
  use_opaque_keys <- any(contains_separator) || anyDuplicated(display_keys) > 0L
  feature_pairs$feature_key <- if (use_opaque_keys) {
    sprintf(".multischolar_peptide_%08d", seq_len(nrow(feature_pairs)))
  } else {
    display_keys
  }

  feature_pairs <- feature_pairs[
    c("feature_key", "active_protein_id", "peptide_sequence")
  ]
  if (anyDuplicated(feature_pairs$feature_key) > 0L ||
      anyDuplicated(feature_pairs[c("active_protein_id", "peptide_sequence")]) > 0L) {
    stop("Peptide feature-key mapping is not one-to-one.", call. = FALSE)
  }

  feature_pairs
}

.peptideProteinAccessionProvenance <- function(theObject) {
  peptide_data <- theObject@peptide_data
  protein_column <- theObject@protein_id_column

  if (!protein_column %in% names(peptide_data)) {
    return(data.frame())
  }

  provenance_columns <- unique(c(protein_column, "Protein.Ids"))
  provenance_columns <- provenance_columns[provenance_columns %in% names(peptide_data)]
  provenance <- peptide_data[, provenance_columns, drop = FALSE]
  for (column in provenance_columns) {
    provenance[[column]] <- as.character(provenance[[column]])
  }
  provenance[!duplicated(provenance), , drop = FALSE]
}

.ensurePeptideFeatureKeyMap <- function(theObject, record_migration = TRUE) {
  expected_map <- .buildPeptideFeatureKeyMap(theObject)
  args <- theObject@args
  if (!is.list(args)) {
    args <- list()
  }

  existing_map <- args$peptide_feature_key_map
  metadata <- args$peptide_feature_key_metadata
  map_columns <- c("feature_key", "active_protein_id", "peptide_sequence")
  existing_valid <- is.data.frame(existing_map) &&
    all(map_columns %in% names(existing_map)) &&
    is.list(metadata) &&
    identical(metadata$protein_id_column, theObject@protein_id_column) &&
    identical(metadata$peptide_sequence_column, theObject@peptide_sequence_column) &&
    anyDuplicated(existing_map$feature_key) == 0L &&
    anyDuplicated(existing_map[c("active_protein_id", "peptide_sequence")]) == 0L

  if (existing_valid) {
    matched_map <- dplyr::left_join(
      expected_map[c("active_protein_id", "peptide_sequence")],
      existing_map[map_columns],
      by = c("active_protein_id", "peptide_sequence")
    )
    if (!anyNA(matched_map$feature_key)) {
      expected_map$feature_key <- as.character(matched_map$feature_key)
    } else {
      existing_valid <- FALSE
    }
  }

  if (!existing_valid && isTRUE(record_migration)) {
    args$peptide_feature_key_migration <- list(
      note = "Feature-key map rebuilt from declared protein and peptide columns.",
      protein_id_column = theObject@protein_id_column,
      peptide_sequence_column = theObject@peptide_sequence_column
    )
  }

  args$peptide_feature_key_map <- expected_map
  args$peptide_feature_key_metadata <- list(
    schema_version = 1L,
    protein_id_column = theObject@protein_id_column,
    peptide_sequence_column = theObject@peptide_sequence_column,
    key_is_display_only = !all(grepl("^\\.multischolar_peptide_", expected_map$feature_key))
  )
  args$protein_accession_provenance <- .peptideProteinAccessionProvenance(theObject)
  theObject@args <- args
  theObject
}

.peptideDataWithFeatureKey <- function(theObject,
                                       peptide_data = theObject@peptide_data,
                                       key_column = ".peptide_feature_key") {
  theObject <- .ensurePeptideFeatureKeyMap(theObject)
  feature_map <- theObject@args$peptide_feature_key_map
  lookup <- feature_map
  names(lookup)[names(lookup) == "active_protein_id"] <- theObject@protein_id_column
  names(lookup)[names(lookup) == "peptide_sequence"] <- theObject@peptide_sequence_column
  names(lookup)[names(lookup) == "feature_key"] <- key_column

  keyed_data <- dplyr::left_join(
    peptide_data,
    lookup,
    by = c(theObject@protein_id_column, theObject@peptide_sequence_column)
  )
  if (anyNA(keyed_data[[key_column]])) {
    stop("A peptide row could not be mapped to a feature key.", call. = FALSE)
  }

  list(theObject = theObject, data = keyed_data)
}

.peptideMatrixToIdentityLong <- function(theObject,
                                         peptide_matrix,
                                         value_column) {
  theObject <- .ensurePeptideFeatureKeyMap(theObject)
  feature_map <- theObject@args$peptide_feature_key_map
  matrix_keys <- rownames(peptide_matrix)
  if (is.null(matrix_keys) || anyDuplicated(matrix_keys) > 0L ||
      !all(matrix_keys %in% feature_map$feature_key)) {
    stop(
      "Peptide matrix row keys do not match the reversible feature-key map; recalculate the peptide matrix.",
      call. = FALSE
    )
  }

  identity_lookup <- feature_map
  names(identity_lookup)[names(identity_lookup) == "active_protein_id"] <- theObject@protein_id_column
  names(identity_lookup)[names(identity_lookup) == "peptide_sequence"] <- theObject@peptide_sequence_column

  long_data <- peptide_matrix |>
    as.data.frame() |>
    tibble::rownames_to_column("feature_key") |>
    tidyr::pivot_longer(
      cols = -feature_key,
      names_to = theObject@sample_id,
      values_to = value_column
    ) |>
    dplyr::left_join(identity_lookup, by = "feature_key") |>
    dplyr::select(-feature_key)

  identity_columns <- c(
    theObject@protein_id_column,
    theObject@peptide_sequence_column,
    theObject@sample_id,
    value_column
  )
  list(
    theObject = theObject,
    data = long_data[, identity_columns, drop = FALSE]
  )
}

#'@export
setMethod(f="calcPeptideMatrix"
          , signature="PeptideQuantitativeData"
          , definition=function( theObject ) {

            keyed <- .peptideDataWithFeatureKey(theObject)
            theObject <- keyed$theObject
            peptide_data <- keyed$data
            # Use the NORMALIZED quantity column for the matrix, not the raw one.
            # This is critical for all downstream stats (limpa, etc.)
            quantity_column_to_use <- theObject@norm_quantity_column
            sample_id_column <- theObject@sample_id
            duplicate_identity <- duplicated(peptide_data[
              c(".peptide_feature_key", sample_id_column)
            ])
            if (any(duplicate_identity)) {
              stop(
                "Peptide matrix input contains duplicate protein/peptide/sample identities.",
                call. = FALSE
              )
            }

            normalised_frozen_peptide_matrix_filt <- peptide_data |>
              dplyr::select(
                .peptide_feature_key,
                !!sym(sample_id_column),
                !!sym(quantity_column_to_use)
              ) |>
              dplyr::mutate(
                !!sym(quantity_column_to_use) := purrr::map_dbl(
                  !!sym(quantity_column_to_use),
                  as.numeric
                )
              ) |>
              tidyr::pivot_wider(
                names_from = !!sym(sample_id_column),
                values_from = !!sym(quantity_column_to_use)
              ) |>
              tibble::column_to_rownames(".peptide_feature_key") |>
              as.matrix()

            theObject@peptide_matrix <- normalised_frozen_peptide_matrix_filt

            # print(head( normalised_frozen_peptide_matrix_filt))

            return(theObject)

          } )
