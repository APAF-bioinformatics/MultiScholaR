# ----------------------------------------------------------------------------
# DIA-NN decoy and contaminant classification
# ----------------------------------------------------------------------------

.peptideExclusionTruthy <- function(values) {
  if (is.logical(values)) {
    return(!is.na(values) & values)
  }
  if (is.numeric(values) || is.integer(values)) {
    return(!is.na(values) & values != 0)
  }

  normalised <- tolower(trimws(as.character(values)))
  !is.na(normalised) & normalised %in% c(
    "1", "true", "t", "yes", "y", "+", "decoy", "contaminant"
  )
}

.asPeptideExclusionFlag <- function(value, argument_name) {
  resolved <- .peptideExclusionTruthy(value)
  if (length(resolved) != 1L || is.na(resolved)) {
    stop(sprintf("%s must be a single TRUE/FALSE value.", argument_name), call. = FALSE)
  }
  resolved[[1]]
}

.findPeptideExclusionColumns <- function(input_names, candidates) {
  candidate_positions <- match(tolower(candidates), tolower(input_names), nomatch = 0L)
  unique(input_names[candidate_positions[candidate_positions > 0L]])
}

.normaliseProteinAccession <- function(accessions) {
  accessions <- trimws(as.character(accessions))
  accessions <- sub("^(?:REV__|DECOY_|CON__|CONTAMINANT_|cRAP-)", "", accessions,
                    ignore.case = TRUE, perl = TRUE)
  pipe_parts <- strsplit(accessions, "|", fixed = TRUE)
  accessions <- vapply(pipe_parts, function(parts) {
    if (length(parts) >= 3L && tolower(parts[[1]]) %in% c("sp", "tr")) {
      parts[[2]]
    } else {
      parts[[1]]
    }
  }, character(1))
  toupper(trimws(accessions))
}

.tokenizeProteinAccessions <- function(protein_text) {
  protein_text <- ifelse(is.na(protein_text), "", as.character(protein_text))
  lapply(strsplit(protein_text, "[;,]", perl = TRUE), function(tokens) {
    tokens <- trimws(tokens)
    tokens <- tokens[nzchar(tokens)]
    if (!length(tokens)) {
      return(character())
    }
    unique(.normaliseProteinAccession(tokens))
  })
}

#' Read a contaminant accession manifest
#'
#' @param manifest `NULL`, a character vector of accessions, a data frame with
#'   an `accession` column, or a path to a one-column/TSV/CSV manifest.
#' @return A list containing normalised accessions and provenance metadata.
#' @keywords internal
readPeptideContaminantManifest <- function(manifest = NULL) {
  empty_manifest <- list(
    accessions = character(),
    source = "none",
    version = NA_character_,
    checksum = NA_character_
  )
  if (is.null(manifest) || (is.character(manifest) &&
      length(manifest) == 1L && !nzchar(trimws(manifest)))) {
    return(empty_manifest)
  }

  manifest_data <- NULL
  source <- "user_vector"
  source_path <- NULL
  if (is.character(manifest) && length(manifest) == 1L && file.exists(manifest)) {
    source_path <- normalizePath(manifest, winslash = "/", mustWork = TRUE)
    source <- paste0("user_file:", basename(source_path))
    manifest_data <- if (grepl("\\.csv$", manifest, ignore.case = TRUE)) {
      utils::read.csv(
        manifest,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    } else {
      utils::read.delim(
        manifest,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  } else if (is.data.frame(manifest)) {
    source <- "user_data_frame"
    manifest_data <- manifest
  } else if (is.character(manifest)) {
    manifest_data <- data.frame(accession = manifest, stringsAsFactors = FALSE)
  } else {
    stop(
      "Contaminant manifest must be NULL, an accession vector, a data frame, or a readable file path.",
      call. = FALSE
    )
  }

  accession_column <- .findPeptideExclusionColumns(
    names(manifest_data),
    c("accession", "uniprot_accession", "protein_accession", "Protein.Ids")
  )
  if (!length(accession_column)) {
    if (ncol(manifest_data) == 1L) {
      accession_column <- names(manifest_data)[[1]]
    } else {
      stop(
        "Contaminant manifest needs an 'accession' column (or exactly one column).",
        call. = FALSE
      )
    }
  }

  accessions <- .normaliseProteinAccession(manifest_data[[accession_column[[1]]]])
  accessions <- sort(unique(accessions[!is.na(accessions) & nzchar(accessions)]))
  if (!length(accessions)) {
    stop("Contaminant manifest contains no usable accessions.", call. = FALSE)
  }

  list(
    accessions = accessions,
    source = source,
    version = if ("version" %in% names(manifest_data)) {
      paste(sort(unique(as.character(manifest_data$version))), collapse = ",")
    } else {
      NA_character_
    },
    checksum = if (!is.null(source_path)) {
      unname(tools::md5sum(source_path))
    } else {
      NA_character_
    }
  )
}

#' Classify DIA-NN rows as decoy or contaminant
#'
#' @param input_table A DIA-NN precursor table.
#' @param protein_id_column Active protein-group column.
#' @param protein_ids_column Protein-accession provenance column.
#' @param contaminant_manifest Optional user-supplied accession manifest.
#' @param exclude_decoys Whether confidently classified decoys are excluded from
#'   the biological analysis view.
#' @param exclude_contaminants Whether confidently classified contaminants are
#'   excluded from the biological analysis view.
#' @return A list containing the immutable input, annotated rows, biological
#'   analysis view, exclusion ledger, summary, and manifest provenance.
#' @export
classifyPeptideBiologicalExclusions <- function(
    input_table,
    protein_id_column = "Protein.Group",
    protein_ids_column = "Protein.Ids",
    contaminant_manifest = NULL,
    exclude_decoys = TRUE,
    exclude_contaminants = TRUE) {
  if (!is.data.frame(input_table)) {
    stop("input_table must be a data frame.", call. = FALSE)
  }
  if (length(protein_id_column) != 1L || !protein_id_column %in% names(input_table)) {
    stop("Active protein_id_column is missing from the DIA-NN table.", call. = FALSE)
  }
  exclude_decoys <- .asPeptideExclusionFlag(exclude_decoys, "exclude_decoys")
  exclude_contaminants <- .asPeptideExclusionFlag(
    exclude_contaminants,
    "exclude_contaminants"
  )

  manifest <- readPeptideContaminantManifest(contaminant_manifest)
  provenance_columns <- unique(c(
    protein_id_column,
    if (protein_ids_column %in% names(input_table)) protein_ids_column else character()
  ))
  combined_ids <- if (nrow(input_table)) {
    apply(
      input_table[provenance_columns],
      1L,
      function(values) paste(values[!is.na(values) & nzchar(trimws(values))], collapse = ";")
    )
  } else {
    character()
  }

  decoy_columns <- .findPeptideExclusionColumns(
    names(input_table),
    c("Decoy", "Reverse", "Is.Decoy", "IsDecoy")
  )
  contaminant_columns <- .findPeptideExclusionColumns(
    names(input_table),
    c(
      "Contaminant", "Potential.contaminant", "Potential.Contaminant",
      "Is.Contaminant", "IsContaminant"
    )
  )

  any_truthy <- function(columns) {
    if (!length(columns)) {
      return(rep(FALSE, nrow(input_table)))
    }
    Reduce(`|`, lapply(columns, function(column) {
      .peptideExclusionTruthy(input_table[[column]])
    }))
  }
  explicit_decoy <- any_truthy(decoy_columns)
  explicit_contaminant <- any_truthy(contaminant_columns)
  decoy_tag <- grepl(
    "(^|[;,[:space:]])(?:REV__|DECOY_)",
    combined_ids,
    ignore.case = TRUE,
    perl = TRUE
  )
  contaminant_tag <- grepl(
    "(^|[;,[:space:]])(?:CON__|CONTAMINANT_|cRAP-)",
    combined_ids,
    ignore.case = TRUE,
    perl = TRUE
  )

  manifest_id_text <- if (protein_ids_column %in% names(input_table)) {
    input_table[[protein_ids_column]]
  } else {
    input_table[[protein_id_column]]
  }
  tokenised <- .tokenizeProteinAccessions(manifest_id_text)
  manifest_match_count <- integer(nrow(input_table))
  manifest_token_count <- lengths(tokenised)
  if (length(manifest$accessions)) {
    manifest_match_count <- vapply(tokenised, function(tokens) {
      if (!length(tokens)) {
        return(0L)
      }
      canonical_tokens <- sub("-[0-9]+$", "", tokens)
      sum(tokens %in% manifest$accessions | canonical_tokens %in% manifest$accessions)
    }, integer(1))
  }
  manifest_any <- manifest_match_count > 0L
  manifest_all <- manifest_any & manifest_match_count == manifest_token_count
  manifest_partial <- manifest_any & !manifest_all

  is_decoy <- explicit_decoy | decoy_tag
  is_contaminant <- explicit_contaminant | contaminant_tag | manifest_all
  exclusion_reason <- dplyr::case_when(
    explicit_decoy ~ "explicit_decoy_column",
    decoy_tag ~ "decoy_identifier_tag",
    explicit_contaminant ~ "explicit_contaminant_column",
    contaminant_tag ~ "contaminant_identifier_tag",
    manifest_all ~ "all_group_accessions_in_contaminant_manifest",
    manifest_partial ~ "mixed_group_partial_contaminant_manifest_match",
    TRUE ~ NA_character_
  )
  classification_source <- dplyr::case_when(
    explicit_decoy ~ paste(decoy_columns, collapse = ","),
    decoy_tag ~ "identifier_tag",
    explicit_contaminant ~ paste(contaminant_columns, collapse = ","),
    contaminant_tag ~ "identifier_tag",
    manifest_any ~ manifest$source,
    TRUE ~ "no_reliable_signal"
  )
  exclude_from_biological_analysis <-
    (is_decoy & exclude_decoys) |
    (is_contaminant & exclude_contaminants)

  annotated <- input_table
  annotated$is_decoy <- is_decoy
  annotated$is_contaminant <- is_contaminant
  annotated$contaminant_manifest_partial_match <- manifest_partial
  annotated$exclude_from_biological_analysis <- exclude_from_biological_analysis
  annotated$exclusion_reason <- exclusion_reason
  annotated$classification_source <- classification_source

  identity_columns <- unique(c(
    ".source_row_id", "Run", "Precursor.Id", protein_id_column,
    if (protein_ids_column %in% names(annotated)) protein_ids_column else character(),
    "Stripped.Sequence", "Modified.Sequence", "is_decoy", "is_contaminant",
    "exclusion_reason", "classification_source"
  ))
  identity_columns <- identity_columns[identity_columns %in% names(annotated)]
  exclusion_ledger <- annotated[exclude_from_biological_analysis, identity_columns, drop = FALSE]

  summary <- data.frame(
    input_rows = nrow(input_table),
    decoy_rows = sum(is_decoy),
    contaminant_rows = sum(is_contaminant),
    partial_manifest_group_rows = sum(manifest_partial),
    excluded_rows = sum(exclude_from_biological_analysis),
    biological_rows = sum(!exclude_from_biological_analysis),
    manifest_accessions = length(manifest$accessions),
    classification_status = if (length(manifest$accessions) ||
        length(decoy_columns) || length(contaminant_columns) ||
        any(decoy_tag) || any(contaminant_tag)) {
      "classified_from_available_signals"
    } else {
      "not_classified_no_reliable_signal"
    },
    stringsAsFactors = FALSE
  )

  list(
    raw_data = input_table,
    classified_data = annotated,
    analysis_data = annotated[!exclude_from_biological_analysis, , drop = FALSE],
    exclusion_ledger = exclusion_ledger,
    summary = summary,
    manifest = manifest
  )
}
