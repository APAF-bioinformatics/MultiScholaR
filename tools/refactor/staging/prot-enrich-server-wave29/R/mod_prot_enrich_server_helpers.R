resolveProtEnrichRunDependencies <- function(currentS4Object,
                                             daResultsData,
                                             workflowData,
                                             contrastsTbl = NULL,
                                             globalEnv = .GlobalEnv) {
  resolvedCurrentS4 <- currentS4Object
  s4Source <- if (!is.null(resolvedCurrentS4)) "current_s4_object" else NULL

  if (is.null(resolvedCurrentS4) &&
      !is.null(daResultsData) &&
      length(daResultsData) > 0) {
    firstResult <- daResultsData[[1]]
    firstObject <- tryCatch(firstResult$theObject, error = function(e) NULL)
    if (!is.null(firstObject)) {
      resolvedCurrentS4 <- firstObject
      s4Source <- "da_results_first_result"
    }
  }

  if (is.null(resolvedCurrentS4) && !is.null(workflowData$state_manager)) {
    currentState <- workflowData$state_manager$current_state
    stateObject <- workflowData$state_manager$getState(currentState)
    if (!is.null(stateObject)) {
      resolvedCurrentS4 <- stateObject
      s4Source <- "state_manager"
    }
  }

  designMatrix <- NULL
  designMatrixSource <- NULL
  designMatrixError <- NULL

  if (!is.null(resolvedCurrentS4)) {
    tryCatch({
      if (!is.null(resolvedCurrentS4@design_matrix)) {
        designMatrix <- resolvedCurrentS4@design_matrix
        designMatrixSource <- "s4_object"
      }
    }, error = function(e) {
      designMatrixError <<- e$message
    })
  }

  if (is.null(designMatrix) && exists("design_matrix", envir = globalEnv)) {
    designMatrix <- get("design_matrix", envir = globalEnv)
    designMatrixSource <- "global_environment"
  }

  list(
    contrastsTbl = contrastsTbl,
    currentS4 = resolvedCurrentS4,
    s4Source = s4Source,
    designMatrix = designMatrix,
    designMatrixSource = designMatrixSource,
    designMatrixError = designMatrixError
  )
}

resolveProtEnrichOutputDirectories <- function(experimentPaths,
                                               dirExistsFn = dir.exists,
                                               dirCreateFn = dir.create,
                                               filePathFn = file.path) {
  daOutputDir <- NULL
  daOutputDirSource <- NULL
  daOutputDirCreated <- FALSE

  if (!is.null(experimentPaths$da_output_dir) &&
      dirExistsFn(experimentPaths$da_output_dir)) {
    daOutputDir <- experimentPaths$da_output_dir
    daOutputDirSource <- "experiment_paths"
  } else {
    daOutputDir <- filePathFn(experimentPaths$results_dir, "da_proteins")
    daOutputDirSource <- "results_fallback"
    if (!dirExistsFn(daOutputDir)) {
      dirCreateFn(daOutputDir, recursive = TRUE)
      daOutputDirCreated <- TRUE
    }
  }

  pathwayDir <- NULL
  pathwayDirSource <- NULL
  pathwayDirCreated <- FALSE

  if (!is.null(experimentPaths$pathway_dir) &&
      dirExistsFn(experimentPaths$pathway_dir)) {
    pathwayDir <- experimentPaths$pathway_dir
    pathwayDirSource <- "experiment_paths"
  } else {
    pathwayDir <- filePathFn(experimentPaths$results_dir, "pathway_enrichment")
    pathwayDirSource <- "results_fallback"
    if (!dirExistsFn(pathwayDir)) {
      dirCreateFn(pathwayDir, recursive = TRUE)
      pathwayDirCreated <- TRUE
    }
  }

  list(
    daOutputDir = daOutputDir,
    daOutputDirSource = daOutputDirSource,
    daOutputDirCreated = daOutputDirCreated,
    pathwayDir = pathwayDir,
    pathwayDirSource = pathwayDirSource,
    pathwayDirCreated = pathwayDirCreated
  )
}

resolveProtEnrichUniprotAnnotations <- function(workflowData,
                                                experimentPaths,
                                                currentS4Object = NULL,
                                                organismTaxid = NULL,
                                                globalEnv = .GlobalEnv,
                                                fileExistsFn = file.exists,
                                                readRdsFn = readRDS,
                                                filePathFn = file.path,
                                                dirExistsFn = dir.exists,
                                                dirCreateFn = dir.create,
                                                getUniprotAnnotationsFn = getUniprotAnnotations,
                                                catFn = cat) {
  uniprotDatCln <- NULL
  annotationSource <- NULL
  sourcePath <- NULL
  cacheDir <- NULL
  cacheDirCreated <- FALSE
  loadError <- NULL
  creationError <- NULL

  if (exists("uniprot_dat_cln", envir = globalEnv)) {
    uniprotDatCln <- get("uniprot_dat_cln", envir = globalEnv)
    annotationSource <- "global_environment"
    catFn("   ENRICHMENT Step: Found uniprot_dat_cln in global environment\n")
  } else if (!is.null(workflowData$uniprot_dat_cln)) {
    uniprotDatCln <- workflowData$uniprot_dat_cln
    assign("uniprot_dat_cln", uniprotDatCln, envir = globalEnv)
    annotationSource <- "workflow_data"
    catFn("   ENRICHMENT Step: Found uniprot_dat_cln in workflow_data\n")
  } else {
    catFn("   ENRICHMENT Step: No uniprot_dat_cln found - checking source directory\n")

    sourcePath <- filePathFn(experimentPaths$source_dir, "uniprot_dat_cln.RDS")
    if (fileExistsFn(sourcePath)) {
      catFn(sprintf("   ENRICHMENT Step: Found uniprot_dat_cln.RDS at %s\n", sourcePath))
      tryCatch({
        uniprotDatCln <- readRdsFn(sourcePath)
        workflowData$uniprot_dat_cln <- uniprotDatCln
        assign("uniprot_dat_cln", uniprotDatCln, envir = globalEnv)
        annotationSource <- "source_directory"
        catFn(sprintf(
          "   ENRICHMENT Step: Successfully loaded %d UniProt annotations from source directory\n",
          nrow(uniprotDatCln)
        ))
      }, error = function(e) {
        loadError <<- e$message
        uniprotDatCln <<- NULL
        catFn(sprintf(
          "   ENRICHMENT Step: Error loading UniProt from source directory: %s\n",
          e$message
        ))
      })
    } else {
      catFn(sprintf(
        "   ENRICHMENT Step: No uniprot_dat_cln.RDS found at %s\n",
        sourcePath
      ))
    }
  }

  if (is.null(uniprotDatCln)) {
    catFn("   ENRICHMENT Step: Attempting to create UniProt annotations on-the-fly\n")

    tryCatch({
      cacheDir <- filePathFn(experimentPaths$results_dir, "cache", "uniprot_annotations")
      if (!dirExistsFn(cacheDir)) {
        dirCreateFn(cacheDir, recursive = TRUE)
        cacheDirCreated <- TRUE
      }

      if (!is.null(currentS4Object) && !is.null(currentS4Object@protein_quant_table)) {
        uniprotDatCln <- getUniprotAnnotationsFn(
          input_tbl = currentS4Object@protein_quant_table,
          cache_dir = cacheDir,
          taxon_id = as.numeric(organismTaxid)
        )

        workflowData$uniprot_dat_cln <- uniprotDatCln
        assign("uniprot_dat_cln", uniprotDatCln, envir = globalEnv)
        annotationSource <- "generated"
        catFn("   ENRICHMENT Step: Successfully created uniprot_dat_cln on-the-fly\n")
      } else {
        catFn("   ENRICHMENT Step: No protein table available for annotation creation\n")
      }
    }, error = function(e) {
      creationError <<- e$message
      catFn(sprintf("   ENRICHMENT Step: Error creating UniProt annotations: %s\n", e$message))
    })
  }

  list(
    uniprotDatCln = uniprotDatCln,
    source = annotationSource,
    sourcePath = sourcePath,
    cacheDir = cacheDir,
    cacheDirCreated = cacheDirCreated,
    loadError = loadError,
    creationError = creationError
  )
}

resolveProtEnrichAnnotationMatching <- function(uniprotDatCln,
                                                daResultsForEnrichment,
                                                currentS4Object,
                                                matchAnnotationsFn = matchAnnotations,
                                                catFn = cat) {
  if (is.null(uniprotDatCln) ||
      is.null(daResultsForEnrichment) ||
      is.null(currentS4Object)) {
    return(list(
      attempted = FALSE,
      proteinIdColumn = NULL,
      annotationMatchResults = NULL,
      matchRate = NULL,
      warning = NULL
    ))
  }

  catFn("   ENRICHMENT Step: Attempting UniProt annotation matching\n")

  proteinIdCol <- tryCatch(
    currentS4Object@protein_id_column,
    error = function(e) "uniprot_acc"
  )
  annotationMatchResults <- NULL
  matchRate <- NULL
  warningMessage <- NULL

  tryCatch({
    annotationMatchResults <- matchAnnotationsFn(
      da_results_s4 = daResultsForEnrichment,
      uniprot_annotations = uniprotDatCln,
      protein_id_column = proteinIdCol,
      uniprot_id_column = "Entry",
      gene_names_column = "gene_names"
    )

    matchRate <- annotationMatchResults$match_statistics$match_rate
    catFn(sprintf(
      "   ENRICHMENT Step: Annotation matching completed - %d%% match rate\n",
      matchRate
    ))
  }, error = function(e) {
    warningMessage <<- e$message
    annotationMatchResults <<- NULL
    catFn(sprintf("   ENRICHMENT Step: Warning in annotation matching: %s\n", e$message))
    catFn("   ENRICHMENT Step: Continuing with enrichment analysis...\n")
  })

  list(
    attempted = TRUE,
    proteinIdColumn = proteinIdCol,
    annotationMatchResults = annotationMatchResults,
    matchRate = matchRate,
    warning = warningMessage
  )
}

