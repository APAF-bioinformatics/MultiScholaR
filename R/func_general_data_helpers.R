#' @title Find Matching Column in Headers (Case-Insensitive)
#' @description Searches for a column name in headers using case-insensitive matching.
#'              Used for auto-populating column mapping inputs.
#'
#' @param headers Character vector of column headers from the data file.
#' @param candidates Character vector of candidate column names to search for.
#'
#' @return The first matching header name (preserving original case), or NULL if no match.
#'
#' @export
findMatchingColumn <- function(headers, candidates) {
    if (is.null(candidates) || length(candidates) == 0) {
        return(NULL)
    }

    headers_lower <- tolower(headers)
    candidates_lower <- tolower(candidates)

    for (cand in candidates_lower) {
        match_idx <- which(headers_lower == cand)
        if (length(match_idx) > 0) {
            return(headers[match_idx[1]])
        }
    }

    return(NULL)
}

# ----------------------------------------------------------------------------
# changeToCategorical
# ----------------------------------------------------------------------------
#' getOneContinousPalette
#' @export
#' @param metadata_tbl,column_name,num_colours Runtime inputs used by this function; see the usage section for accepted values.
changeToCategorical <- function(metadata_tbl, column_name, num_colours=9) {

  list_of_values <-  metadata_tbl |>
    dplyr::select( all_of(column_name)  ) |>
    #dplyr::filter(!is.na(!!sym(column_name))) |>
    pull()

  min_value <- min(list_of_values, na.rm =TRUE)
  max_value <-  max(list_of_values, na.rm =TRUE)

  if(min_value > 1) {
    min_value <- floor(min_value)
  }

  if(max_value > 1) {
    max_value <- ceiling(max_value)
  }

  formatted_list_of_values <- cut(list_of_values, breaks=seq( min_value, max_value, length.out=num_colours) )

  formatted_list_of_values
}

# ----------------------------------------------------------------------------
# peptidesIntensityMatrixPivotLonger
# ----------------------------------------------------------------------------
#' @title Pivot Peptide Matrix Long
#' @description
#'  Pivot peptide intensity matrix into long format table.
#' @export
#' @param input_matrix,sample_id_column,sequence_column,protein_id_column,quantity_column,unlog_data Runtime inputs used by this function; see the usage section for accepted values.
peptidesIntensityMatrixPivotLonger <- function( input_matrix
                                                , sample_id_column
                                                , sequence_column
                                                , protein_id_column
                                                , quantity_column
                                                , unlog_data = TRUE) {

  output_matrix <- input_matrix |>
    as.data.frame() |>
    rownames_to_column(protein_id_column) |>
    pivot_longer(cols=!contains(protein_id_column)
                 , names_to = sample_id_column
                 , values_to =  quantity_column ) |>
    separate( col = protein_id_column
              , into=c(protein_id_column, "Stripped.Sequence"), sep="_")

  if ( unlog_data == TRUE) {
    output_matrix <- output_matrix|>
      mutate( {{quantity_column}} := 2^(!!sym(quantity_column)))

  }

  output_matrix
}

# ----------------------------------------------------------------------------
# proteinIntensityMatrixPivotLonger
# ----------------------------------------------------------------------------
#' @title Pivot Protein Matrix Long
#' @description
#' Pivot protein intensity matrix into long format
#' @export
#' @param input_matrix,sample_id_column,protein_id_column,quantity_column Runtime inputs used by this function; see the usage section for accepted values.
proteinIntensityMatrixPivotLonger <- function( input_matrix
                                               , sample_id_column
                                               , protein_id_column
                                               , quantity_column) {

  output_matrix <- input_matrix |>
    as.data.frame() |>
    rownames_to_column(protein_id_column) |>
    pivot_longer(cols=!contains(protein_id_column)
                 , names_to = sample_id_column
                 , values_to =  quantity_column )

  output_matrix
}

# ----------------------------------------------------------------------------
# createIdToAttributeHash
# ----------------------------------------------------------------------------
## Output:
## An environment that act as a hash to convert keys to attributes.
#' @export
createIdToAttributeHash <- function(keys, attributes) {
	keys <- as.character( as.vector(keys))
	attribute <- as.vector(attributes)

	hash <- new.env(hash = TRUE, parent = parent.frame())

	if ( length(keys) != length(attributes))  {
		warning('Length of keys != Length of attributes list.')
		return(1)
	}

	mapply(function(k, a) base::assign(k, a, envir = hash), 
		   keys, attributes)

	return(hash)
}

# ----------------------------------------------------------------------------
# convertKeyToAttribute
# ----------------------------------------------------------------------------
## Output:
## A value that correspond to the query key value.
#' @export
convertKeyToAttribute <- function(key, hash) {

  if ( base::exists(key, hash) ) {
    return ( base::get(key, hash))
  }
  else {
    return (NA)
  }
}

# ----------------------------------------------------------------------------
# extract_experiment
# ----------------------------------------------------------------------------
##################################################################################################################
#' @title Extract Substrings from Underscore-Separated Strings
#'
#' @description
#' Extracts substrings from underscore-separated strings using different modes:
#' range (between positions), start (first element), or end (last element).
#'
#' @param x Character vector containing the strings to process
#' @param mode Character string specifying extraction mode:
#'   * "range": Extract elements between two underscore positions
#'   * "start": Extract from start to first underscore
#'   * "end": Extract from last underscore to end
#' @param start Integer: Starting position for range mode (1-based)
#' @param end Integer: Ending position for range mode (1-based, required for range mode)
#'
#' @return Character vector with extracted strings. Returns NA for strings where
#'   requested positions are out of bounds.
#'
#' @examples
#' x <- "20140602_ffs_expt1_r1_junk"
#' extract_experiment(x, mode = "range", start = 1, end = 3)  # "20140602_ffs_expt1"
#' extract_experiment(x, mode = "start")  # "20140602"
#' extract_experiment(x, mode = "end")    # "junk"
#'
#' # Multiple strings
#' x <- c("20140602_ffs_expt1_r1_junk", "20140603_ffs_expt2_r2_test")
#' extract_experiment(x, mode = "range", start = 2, end = 3)  # c("ffs_expt1", "ffs_expt2")
#'
#' @export
extract_experiment <- function(x, mode = "range", start = 1, end = NULL) {
  if (!mode %in% c("range", "start", "end")) {
    stop("Mode must be one of: 'range', 'start', 'end'")
  }

  process_string <- function(str) {
    parts <- unlist(strsplit(str, "_"))

    if (mode == "range") {
      if (is.null(end)) stop("End position required for range mode")
      if (start > length(parts) || end > length(parts)) {
        warning("Position out of bounds for string: ", str)

        return(NA_character_)
      }
      return(paste(parts[start:end], collapse = "_"))
    }

    else if (mode == "start") {
      return(parts[1])
    }

    else if (mode == "end") {
      return(parts[length(parts)])
    }
  }

  sapply(x, process_string)
}

# ----------------------------------------------------------------------------
# calcHtSize
# ----------------------------------------------------------------------------
#' @export
calcHtSize = function(ht, unit = "inch") {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  ComplexHeatmap::draw(ht)
  grDevices::dev.size(unit = unit)
}
