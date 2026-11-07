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
  # Genuine MaxQuant output carries spaced names ("Majority protein IDs"); the markers are dotted (KIT-054).
  maxquant_found <- sum(maxquant_markers %in% unique(c(headers_lower, gsub(" ", ".", headers_lower, fixed = TRUE))))
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

#' Inspect Parquet schema without reading columns or retaining an open file.
#' @param path Parquet file being inspected.
#' @param col_select Existing internal reader callback argument; unused for headers.
#' @return Arrow schema supporting the existing names() header contract.
#' @noRd
readProtImportParquetSchema <- function(path, col_select = NULL) {
  stream <- arrow::ReadableFile$create(path)
  on.exit(stream$close(), add = TRUE)
  reader <- arrow::ParquetFileReader$create(stream)
  reader$GetSchema()
}

readProtImportHeaders <- function(
    filePath
    , unzipList = function(path) unzip(path, list = TRUE)
    , unzipFiles = unzip
    , fileExt = tools::file_ext
    , requireNamespaceFn = requireNamespace
    , readExcel = readxl::read_excel
    , readParquet = readProtImportParquetSchema
    , readLinesFn = readLines
    , tempfileFn = tempfile
    , dirCreate = dir.create
    , unlinkFn = unlink
    , basenameFn = basename
) {
  headers <- NULL

  if (tolower(fileExt(filePath)) == "zip") {
    zip_contents <- unzipList(filePath)
    first_data_file <- zip_contents$Name[grep("\\.(xlsx|tsv|csv)$", zip_contents$Name, ignore.case = TRUE)[1]]

    if (is.na(first_data_file)) {
      stop("No data files (.xlsx, .tsv, .csv) found inside the ZIP archive.")
    }

    if (tolower(fileExt(first_data_file)) == "xlsx") {
      if (!requireNamespaceFn("readxl", quietly = TRUE)) {
        stop("Package 'readxl' needed for .xlsx files.")
      }

      temp_dir <- tempfileFn()
      dirCreate(temp_dir)
      on.exit(unlinkFn(temp_dir, recursive = TRUE), add = TRUE)

      unzipFiles(filePath, files = first_data_file, exdir = temp_dir, junkpaths = TRUE)
      unzipped_file_path <- file.path(temp_dir, basenameFn(first_data_file))
      headers <- names(readExcel(unzipped_file_path, n_max = 0))
    } else {
      con <- unz(filePath, first_data_file)
      on.exit(close(con), add = TRUE)
      headers <- strsplit(readLinesFn(con, n = 1), "\t|,")[[1]]
    }
  } else if (grepl("\\.parquet$", filePath, ignore.case = TRUE)) {
    headers <- names(readParquet(filePath, col_select = NULL))
  } else {
    preview_lines <- readLinesFn(filePath, n = 10)
    headers <- strsplit(preview_lines[1], "\t|,")[[1]]
  }

  if (is.null(headers)) {
    stop("Could not read headers from the provided file.")
  }

  headers
}

#' Read configuration without eager pools only in an owned DIA GUI session.
#' @param file,file_type Existing configuration-reader arguments.
#' @return Unchanged configuration with portable workers for the DIA adapter.
#' @noRd
readProtGuiConfig <- function(file, file_type = "ini") {
  readConfigFile(file = file, file_type = file_type,
    create_workers = is.null(getDIAActiveManager()))
}

loadProtImportConfiguration <- function(
    configPath = NULL
    , experimentPaths
    , readConfig = readConfigFile
    , showNotification = shiny::showNotification
    , logInfo = logger::log_info
    , logError = logger::log_error
    , getDefaultConfig = getDefaultProteomicsConfig
    , fileExists = file.exists
    , fileCopy = file.copy
    , downloadFile = download.file
    , systemFileFn = system.file
) {
  if (!is.null(configPath)) {
    logInfo(paste("Reading configuration from", configPath))
    return(readConfig(file = configPath))
  }

  logInfo("Using default configuration")
  default_config_path <- file.path(experimentPaths$source_dir, "config.ini")

  if (fileExists(default_config_path)) {
    return(readConfig(file = default_config_path))
  }

  logInfo("config.ini not found in project. Retrieving default config.")
  tryCatch({
    pkg_config <- systemFileFn("config", "config.ini", package = "MultiScholaR")

    if (fileExists(pkg_config) && pkg_config != "") {
      logInfo(paste("Found default config.ini in package:", pkg_config))
      fileCopy(pkg_config, default_config_path)
      logInfo(paste("Default config.ini copied to:", default_config_path))
      showNotification("Copied default config.ini to scripts directory.", type = "message")
      readConfig(file = default_config_path)
    } else {
      logInfo("Default config.ini not found in package, downloading from GitHub...")
      default_config_url <- "https://raw.githubusercontent.com/APAF-bioinformatics/MultiScholaR/main/inst/config/config.ini"
      downloadFile(default_config_url, destfile = default_config_path, quiet = TRUE)
      logInfo(paste("Default config.ini downloaded to:", default_config_path))
      showNotification("Downloaded default config.ini to scripts directory.", type = "message")
      readConfig(file = default_config_path)
    }
  }, error = function(e_download) {
    msg <- paste("Failed to retrieve default config.ini:", e_download$message)
    logError(msg)
    showNotification(msg, type = "warning", duration = 10)
    logInfo("Using minimal fallback configuration")
    getDefaultConfig()
  })
}

readProtImportOptionalMapping <- function(
    mappingPath = NULL
    , mappingLabel
    , readMapping = function(path) vroom::vroom(path, show_col_types = FALSE)
    , logInfo = logger::log_info
    , logError = logger::log_error
) {
  if (is.null(mappingPath)) {
    return(NULL)
  }

  logInfo(paste("Reading", mappingLabel, "from", mappingPath))

  tryCatch({
    readMapping(mappingPath)
  }, error = function(e) {
    logError(paste("Failed to read", mappingLabel, "file:", e$message))
    stop("Failed to read ", mappingLabel, " file: ", e$message)
  })
}

loadProtImportOptionalMappings <- function(
    workflowData
    , uniprotPath = NULL
    , uniparcPath = NULL
    , readOptionalMapping = readProtImportOptionalMapping
    , logInfo = logger::log_info
    , logError = logger::log_error
) {
  if (!is.null(uniprotPath)) {
    workflowData$uniprot_mapping <- readOptionalMapping(
      mappingPath = uniprotPath
      , mappingLabel = "UniProt mapping"
      , readMapping = function(path) vroom::vroom(path, show_col_types = FALSE)
      , logInfo = logInfo
      , logError = logError
    )
  }

  if (!is.null(uniparcPath)) {
    workflowData$uniparc_mapping <- readOptionalMapping(
      mappingPath = uniparcPath
      , mappingLabel = "UniParc mapping"
      , readMapping = function(path) vroom::vroom(path, show_col_types = FALSE)
      , logInfo = logInfo
      , logError = logError
    )
  }

  invisible(list(
    uniprotMapping = workflowData$uniprot_mapping
    , uniparcMapping = workflowData$uniparc_mapping
  ))
}
