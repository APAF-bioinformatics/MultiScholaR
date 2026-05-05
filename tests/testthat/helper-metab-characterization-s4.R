makeMetabCharacterizationDesign <- function(
  sample_ids = c("Sample_1", "Sample_2"),
  sample_id_col = "sample_id",
  group_id_col = "group"
) {
  groups <- rep(c("groupA", "groupB"), length.out = length(sample_ids))
  design_matrix <- data.frame(
    sample_id = as.character(sample_ids),
    group = groups,
    stringsAsFactors = FALSE
  )
  names(design_matrix) <- c(sample_id_col, group_id_col)
  design_matrix
}

makeMetabCharacterizationAssay <- function(
  metabolite_ids = c("M1", "M2"),
  sample_ids = c("Sample_1", "Sample_2"),
  metabolite_id_column = "metabolite_id",
  annotation_id_column = "metabolite_name",
  extra_columns = list()
) {
  assay <- data.frame(
    metabolite_id = as.character(metabolite_ids),
    metabolite_name = paste("Met", seq_along(metabolite_ids)),
    stringsAsFactors = FALSE
  )
  names(assay)[names(assay) == "metabolite_id"] <- metabolite_id_column
  names(assay)[names(assay) == "metabolite_name"] <- annotation_id_column

  for (sample_id in sample_ids) {
    assay[[sample_id]] <- seq_along(metabolite_ids)
  }
  for (column_name in names(extra_columns)) {
    assay[[column_name]] <- extra_columns[[column_name]]
  }

  assay
}

inferMetabCharacterizationSamples <- function(metabolite_data, metabolite_id_column) {
  if (is.list(metabolite_data) && !is.data.frame(metabolite_data) && length(metabolite_data) > 0) {
    assay <- metabolite_data[[1]]
  } else if (is.data.frame(metabolite_data)) {
    assay <- metabolite_data
  } else {
    return(c("Sample_1", "Sample_2"))
  }

  candidate_cols <- setdiff(colnames(assay), metabolite_id_column)
  numeric_cols <- candidate_cols[vapply(assay[candidate_cols], is.numeric, logical(1))]
  if (length(numeric_cols) > 0) {
    return(numeric_cols[[1]])
  }

  candidate_cols[[1]]
}

makeMetabCharacterizationObject <- function(
  metabolite_data = NULL,
  design_matrix = NULL,
  sample_ids = NULL,
  sample_id = "sample_id",
  group_id = "group",
  metabolite_id_column = "metabolite_id",
  annotation_id_column = "metabolite_name",
  database_identifier_type = "InternalName",
  internal_standard_regex = NA_character_,
  technical_replicate_id = NA_character_,
  args = list()
) {
  if (is.null(metabolite_data)) {
    metabolite_data <- list()
  }
  if (is.null(sample_ids)) {
    sample_ids <- inferMetabCharacterizationSamples(metabolite_data, metabolite_id_column)
  }
  if (is.null(design_matrix)) {
    design_matrix <- makeMetabCharacterizationDesign(
      sample_ids = sample_ids,
      sample_id_col = sample_id,
      group_id_col = group_id
    )
  }

  methods::new(
    "MetaboliteAssayData",
    metabolite_data = metabolite_data,
    metabolite_id_column = metabolite_id_column,
    annotation_id_column = annotation_id_column,
    database_identifier_type = database_identifier_type,
    internal_standard_regex = internal_standard_regex,
    design_matrix = design_matrix,
    sample_id = sample_id,
    group_id = group_id,
    technical_replicate_id = technical_replicate_id,
    args = args
  )
}

makeMetabCharacterizationObjectWithBareCounts <- function(counts_table) {
  object <- makeMetabCharacterizationObject(
    metabolite_data = list(counts_table),
    sample_ids = inferMetabCharacterizationSamples(counts_table, "metabolite_id")
  )
  object@metabolite_data <- counts_table
  object
}
