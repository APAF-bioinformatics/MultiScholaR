# ----------------------------------------------------------------------------
# detectProteomicsFormat
# ----------------------------------------------------------------------------
#' Detect Proteomics Data Format
#' 
#' Detects the format of proteomics search results based on headers and filename
#' 
#' @param headers Character vector of column headers
#' @param filename Name of the file
#' @param preview_lines First few lines of the file
#' 
#' @return List with format and confidence score
#' @export
detectProteomicsFormat <- function(headers, filename, preview_lines = NULL) {
  # Convert to lowercase for comparison
  headers_lower <- tolower(headers)
  filename_lower <- tolower(filename)
  
  # DIA-NN detection
  diann_score <- 0
  diann_markers <- c("protein.group", "protein.ids", "protein.names", 
                     "precursor.id", "modified.sequence", "stripped.sequence",
                     "precursor.charge", "q.value", "pg.q.value", "run")
  diann_found <- sum(diann_markers %in% headers_lower)
  diann_score <- diann_found / length(diann_markers)
  if (grepl("\\.parquet$", filename_lower)) diann_score <- diann_score + 0.2
  
  # Spectronaut detection
  spectronaut_score <- 0
  spectronaut_markers <- c("pg.proteingroups", "pg.proteinaccessions", "pg.genes",
                          "eg.precursorid", "eg.modifiedsequence", "fg.charge",
                          "eg.qvalue", "pg.qvalue", "r.filename", "pg.quantity")
  spectronaut_found <- sum(spectronaut_markers %in% headers_lower)
  spectronaut_score <- spectronaut_found / length(spectronaut_markers)
  
  # FragPipe detection
  fragpipe_score <- 0
  
  # Check for key markers (flexible matching)
  has_protein_id <- any(grepl("^protein\\s+id$|^protein\\.id$|^protein_id$", headers_lower))
  has_protein <- "protein" %in% headers_lower
  has_gene <- "gene" %in% headers_lower
  has_description <- any(grepl("description", headers_lower))
  has_spectral_count <- any(grepl("spectral.*count", headers_lower))
  
  # Strong indicator: columns ending with "Intensity" (case-insensitive)
  intensity_cols <- sum(grepl("intensity$", headers_lower))
  has_intensity_cols <- intensity_cols > 0
  
  # Check for MaxLFQ Intensity columns specifically
  has_maxlfq <- any(grepl("maxlfq.*intensity$", headers_lower))
  
  # Count basic markers found
  basic_markers_found <- sum(c(has_protein_id, has_protein, has_gene, has_description, has_spectral_count))
  
  # Weighted scoring similar to TMT detection
  # Basic markers worth 40% (8% each)
  basic_score <- (basic_markers_found / 5) * 0.4
  
  # Intensity columns presence is worth 40% (strong indicator)
  intensity_score <- if (has_intensity_cols) {
    # Bonus if multiple intensity columns (typical of FragPipe)
    min(0.4, 0.3 + (min(intensity_cols, 10) / 10) * 0.1)
  } else {
    0
  }
  
  # MaxLFQ presence is worth 20% (very specific to FragPipe)
  maxlfq_score <- if (has_maxlfq) 0.2 else 0
  
  fragpipe_score <- basic_score + intensity_score + maxlfq_score
  
  # Filename bonus (cap at 1.0)
  if (grepl("fragpipe|msfragger", filename_lower)) {
    fragpipe_score <- min(1.0, fragpipe_score + 0.1)
  }
  
  # MaxQuant detection
  maxquant_score <- 0
  maxquant_markers <- c("proteins", "majority.protein.ids", "protein.names",
                       "gene.names", "peptide.counts", "unique.peptides",
                       "intensity", "lfq.intensity", "ms.ms.count")
  maxquant_found <- sum(maxquant_markers %in% headers_lower)
  maxquant_score <- maxquant_found / length(maxquant_markers)
  if (grepl("proteingroups", filename_lower)) maxquant_score <- maxquant_score + 0.3
  
  # PD-TMT detection
  pd_tmt_score <- 0
  pd_tmt_markers <- c("protein fdr confidence", "master", "accession", "exp. q-value", "sum pep score")
  pd_tmt_found <- sum(pd_tmt_markers %in% headers_lower)
  abundance_found <- any(grepl("^abundance:", headers_lower))
  
  # Weighted scoring to prevent exceeding 100%
  # Header markers are worth 60% total (12% each)
  header_score <- (pd_tmt_found / length(pd_tmt_markers)) * 0.6
  # Abundance column presence is worth 40%
  abundance_score <- if (abundance_found) 0.4 else 0
  
  pd_tmt_score <- header_score + abundance_score
  
  # Determine best match
  scores <- c(diann = diann_score, 
              spectronaut = spectronaut_score,
              fragpipe = fragpipe_score,
              maxquant = maxquant_score,
              pd_tmt = pd_tmt_score)
  
  best_format <- names(which.max(scores))
  best_score <- max(scores)
  
  if (best_score < 0.3) {
    best_format <- "unknown"
  }
  
  return(list(
    format = best_format,
    confidence = best_score,
    all_scores = scores
  ))
}

# ----------------------------------------------------------------------------
# getDefaultProteomicsConfig
# ----------------------------------------------------------------------------
#' Get Default Proteomics Configuration
#' 
#' Returns default configuration settings for proteomics analysis
#' 
#' @return List with default configuration parameters
#' @export
getDefaultProteomicsConfig <- function() {
  list(
    generalParameters = list(
      min_peptides_per_protein = 2,
      min_peptides_per_sample = 2,
      q_value_threshold = 0.01,
      intensity_threshold = 0
    ),
    deAnalysisParameters = list(
      formula_string = "~ 0 + group",
      q_value_threshold = 0.05,
      log2_fc_threshold = 1
    ),
    normalizationParameters = list(
      normalisation_method = "cyclicloess"
    ),
    ruvParameters = list(
      percentage_as_neg_ctrl = 33,
      ruv_k = NULL
    )
  )
}

