metabDesignScalarString <- function(value) {
    is.character(value) && length(value) == 1L && !is.na(value) && nzchar(trimws(value))
}

metabDesignNonSamplePatterns <- function() {
    c(
        "^Alignment.ID$",
        "^Average.Rt",
        "^Average.Mz",
        "^Peak ID$",
        "^RT \\(min\\)$",
        "^Precursor m/z$",
        "^Metabolite.name$",
        "^Metabolite name$",
        "^Adduct.type$",
        "^Adduct$",
        "^Post.curation.result$",
        "^Fill.%$",
        "^MS/MS.spectrum$",
        "^Reference.RT$",
        "^Reference.m/z$",
        "^Formula$",
        "^Ontology$",
        "^INCHIKEY$",
        "^SMILES$",
        "^Annotation",
        "^Annotation.tag",
        "^RT.matched$",
        "^m/z.matched$",
        "^MS/MS.matched$",
        "^Comment$",
        "^Manually.modified$",
        "^Total.score$",
        "^Total Score$",
        "^RT.similarity$",
        "^Dot.product$",
        "^Reverse.dot.product$",
        "^Fragment.presence",
        "^S/N.average$",
        "^Spectrum.reference",
        "^MS1.isotopic.spectrum$",
        "^Ion mode$"
    )
}

metabDesignSampleColumnsForAssay <- function(assayData, columnMapping = NULL) {
    if (!is.data.frame(assayData)) {
        return(character())
    }

    if (!is.null(columnMapping$sample_columns) && length(columnMapping$sample_columns) > 0L) {
        return(as.character(columnMapping$sample_columns))
    }

    sampleCols <- names(assayData)[vapply(assayData, is.numeric, logical(1))]
    excludeCols <- c(columnMapping$metabolite_id_col, columnMapping$annotation_col)
    excludeCols <- excludeCols[!is.na(excludeCols) & nzchar(excludeCols)]
    sampleCols <- setdiff(sampleCols, excludeCols)

    for (pattern in metabDesignNonSamplePatterns()) {
        sampleCols <- sampleCols[!grepl(pattern, sampleCols, ignore.case = TRUE)]
    }

    sampleCols
}

validateMetabDesignAlignment <- function(designMatrix,
                                         assayList,
                                         columnMapping = NULL,
                                         sampleIdCol = "Run",
                                         groupCol = "group",
                                         requireGroups = TRUE,
                                         requireReplicates = TRUE) {
    errors <- character()
    warnings <- character()

    if (!is.data.frame(designMatrix) || nrow(designMatrix) == 0L) {
        errors <- c(errors, "Design matrix must be a non-empty data frame.")
    }
    if (is.null(assayList) || length(assayList) == 0L || !is.list(assayList)) {
        errors <- c(errors, "Metabolomics design requires a non-empty assay list.")
    }

    if (length(errors) > 0L) {
        return(list(valid = FALSE, errors = unique(errors), warnings = unique(warnings)))
    }

    assayNames <- names(assayList)
    if (is.null(assayNames)) {
        assayNames <- paste0("assay_", seq_along(assayList))
    } else {
        missingAssayNames <- !nzchar(assayNames)
        assayNames[missingAssayNames] <- paste0("assay_", which(missingAssayNames))
    }
    if (anyDuplicated(assayNames) > 0L) {
        errors <- c(errors, sprintf(
            "Metabolomics design assay names must be unique: %s",
            paste(unique(assayNames[duplicated(assayNames)]), collapse = ", ")
        ))
    }

    if (!sampleIdCol %in% names(designMatrix)) {
        errors <- c(errors, sprintf("Design matrix is missing sample column '%s'.", sampleIdCol))
        designSamples <- character()
    } else {
        designSamples <- as.character(designMatrix[[sampleIdCol]])
        if (any(is.na(designSamples) | !nzchar(trimws(designSamples)))) {
            errors <- c(errors, sprintf("Design matrix sample column '%s' contains blank sample IDs.", sampleIdCol))
        }
        if (anyDuplicated(designSamples) > 0L) {
            duplicatedSamples <- unique(designSamples[duplicated(designSamples)])
            errors <- c(errors, sprintf("Design matrix contains duplicate sample rows: %s", paste(duplicatedSamples, collapse = ", ")))
        }
    }

    if (isTRUE(requireGroups)) {
        if (!groupCol %in% names(designMatrix)) {
            errors <- c(errors, sprintf("Design matrix is missing group column '%s'.", groupCol))
        } else {
            groups <- as.character(designMatrix[[groupCol]])
            if (any(is.na(groups) | !nzchar(trimws(groups)))) {
                errors <- c(errors, sprintf("Design matrix group column '%s' contains blank assignments.", groupCol))
            }
        }
    }

    if (isTRUE(requireReplicates)) {
        if (!"replicates" %in% names(designMatrix)) {
            errors <- c(errors, "Design matrix is missing required 'replicates' column.")
        } else if (any(is.na(designMatrix$replicates))) {
            errors <- c(errors, "Design matrix contains missing replicate assignments.")
        }
    }

    assaySummaries <- lapply(seq_along(assayList), function(assayIdx) {
        assayName <- assayNames[[assayIdx]]
        assayData <- assayList[[assayIdx]]
        assayErrors <- character()
        assayWarnings <- character()

        if (!is.data.frame(assayData)) {
            assayErrors <- c(assayErrors, sprintf("Assay '%s' is not a data frame.", assayName))
            return(list(
                assay_name = assayName,
                sample_columns = character(),
                errors = assayErrors,
                warnings = assayWarnings
            ))
        }

        sampleCols <- metabDesignSampleColumnsForAssay(assayData, columnMapping)
        missingMappedCols <- setdiff(sampleCols, names(assayData))
        presentSampleCols <- intersect(sampleCols, names(assayData))
        inferenceMapping <- columnMapping
        if (is.null(inferenceMapping)) {
            inferenceMapping <- list()
        }
        inferenceMapping$sample_columns <- NULL
        inferredSampleCols <- metabDesignSampleColumnsForAssay(assayData, inferenceMapping)
        extraAssaySampleCols <- setdiff(inferredSampleCols, sampleCols)

        if (length(sampleCols) == 0L) {
            assayErrors <- c(assayErrors, sprintf("Assay '%s' has no sample columns available for design alignment.", assayName))
        }
        if (length(missingMappedCols) > 0L) {
            assayErrors <- c(
                assayErrors,
                sprintf(
                    "Assay '%s' is missing mapped sample column(s): %s",
                    assayName,
                    paste(missingMappedCols, collapse = ", ")
                )
            )
        }
        if (length(extraAssaySampleCols) > 0L) {
            assayErrors <- c(
                assayErrors,
                sprintf(
                    "Assay '%s' has sample-like column(s) absent from column mapping/design: %s",
                    assayName,
                    paste(extraAssaySampleCols, collapse = ", ")
                )
            )
        }
        if (anyDuplicated(presentSampleCols) > 0L) {
            duplicatedAssaySamples <- unique(presentSampleCols[duplicated(presentSampleCols)])
            assayErrors <- c(
                assayErrors,
                sprintf(
                    "Assay '%s' contains duplicate sample columns: %s",
                    assayName,
                    paste(duplicatedAssaySamples, collapse = ", ")
                )
            )
        }

        missingInDesign <- setdiff(presentSampleCols, designSamples)
        extraInDesign <- setdiff(designSamples, presentSampleCols)
        caseCandidateCols <- unique(c(presentSampleCols, extraAssaySampleCols))
        caseMismatchDesign <- designSamples[
            tolower(designSamples) %in% tolower(caseCandidateCols) &
                !designSamples %in% caseCandidateCols
        ]
        caseMismatchAssay <- caseCandidateCols[
            tolower(caseCandidateCols) %in% tolower(designSamples) &
                !caseCandidateCols %in% designSamples
        ]

        if (length(missingInDesign) > 0L) {
            assayErrors <- c(
                assayErrors,
                sprintf(
                    "Assay '%s' has sample column(s) missing from design: %s",
                    assayName,
                    paste(missingInDesign, collapse = ", ")
                )
            )
        }
        if (length(extraInDesign) > 0L) {
            assayErrors <- c(
                assayErrors,
                sprintf(
                    "Design has sample row(s) absent from assay '%s': %s",
                    assayName,
                    paste(extraInDesign, collapse = ", ")
                )
            )
        }
        if (length(caseMismatchDesign) > 0L || length(caseMismatchAssay) > 0L) {
            assayErrors <- c(
                assayErrors,
                sprintf(
                    "Assay '%s' has case-varied sample names between design and assay: design [%s], assay [%s]",
                    assayName,
                    paste(unique(caseMismatchDesign), collapse = ", "),
                    paste(unique(caseMismatchAssay), collapse = ", ")
                )
            )
        }
        if (length(presentSampleCols) > 0L && setequal(presentSampleCols, designSamples) && !identical(presentSampleCols, designSamples)) {
            assayWarnings <- c(
                assayWarnings,
                sprintf("Assay '%s' sample columns are reordered relative to the design matrix.", assayName)
            )
        }

        list(
            assay_name = assayName,
            sample_columns = presentSampleCols,
            errors = assayErrors,
            warnings = assayWarnings
        )
    })

    for (summary in assaySummaries) {
        errors <- c(errors, summary$errors)
        warnings <- c(warnings, summary$warnings)
    }

    list(
        valid = length(errors) == 0L,
        errors = unique(errors),
        warnings = unique(warnings),
        design_samples = designSamples,
        assay_summaries = assaySummaries
    )
}

metabDesignContrastExpression <- function(contrastString) {
    sub("^[^=]+=", "", trimws(as.character(contrastString)))
}

metabDesignFriendlyContrastName <- function(contrastString) {
    contrastString <- metabDesignContrastExpression(contrastString)
    cleanString <- gsub("^group", "", contrastString)
    cleanString <- gsub("-group", "-", cleanString)
    cleanString <- gsub("\\s+", "", cleanString)
    gsub("-", "_vs_", cleanString)
}

normaliseMetabDesignContrastsTable <- function(contrastsTbl, allowEmpty = FALSE) {
    if (is.null(contrastsTbl) || !is.data.frame(contrastsTbl) || nrow(contrastsTbl) == 0L) {
        if (isTRUE(allowEmpty)) {
            return(data.frame(
                contrasts = character(),
                friendly_names = character(),
                full_format = character(),
                stringsAsFactors = FALSE
            ))
        }
        stop("Metabolomics design contrasts table must be a non-empty data frame.", call. = FALSE)
    }

    rawInput <- if ("contrasts" %in% names(contrastsTbl)) {
        contrastsTbl$contrasts
    } else if ("contrast_string" %in% names(contrastsTbl)) {
        contrastsTbl$contrast_string
    } else {
        contrastsTbl[[1L]]
    }
    rawInput <- as.character(rawInput)
    rawContrasts <- vapply(rawInput, metabDesignContrastExpression, character(1))

    if (any(is.na(rawContrasts) | !nzchar(trimws(rawContrasts)))) {
        stop("Metabolomics design contrasts table contains an empty contrast.", call. = FALSE)
    }
    if (anyDuplicated(rawContrasts) > 0L) {
        stop("Metabolomics design contrasts table contains duplicate contrast definitions.", call. = FALSE)
    }

    fullFormat <- if ("full_format" %in% names(contrastsTbl)) {
        as.character(contrastsTbl$full_format)
    } else {
        rep(NA_character_, length(rawContrasts))
    }
    friendlyNames <- if ("friendly_names" %in% names(contrastsTbl)) {
        as.character(contrastsTbl$friendly_names)
    } else {
        rep(NA_character_, length(rawContrasts))
    }

    for (idx in seq_along(rawContrasts)) {
        if (is.na(friendlyNames[[idx]]) || !nzchar(trimws(friendlyNames[[idx]]))) {
            if (!is.na(fullFormat[[idx]]) && grepl("=", fullFormat[[idx]], fixed = TRUE)) {
                friendlyNames[[idx]] <- sub("=.*$", "", fullFormat[[idx]])
            } else if (grepl("=", rawInput[[idx]], fixed = TRUE)) {
                friendlyNames[[idx]] <- sub("=.*$", "", rawInput[[idx]])
            } else {
                friendlyNames[[idx]] <- metabDesignFriendlyContrastName(rawContrasts[[idx]])
            }
        }
        if (is.na(fullFormat[[idx]]) || !nzchar(trimws(fullFormat[[idx]]))) {
            fullFormat[[idx]] <- paste0(friendlyNames[[idx]], "=", rawContrasts[[idx]])
        }
    }

    if (anyDuplicated(friendlyNames) > 0L) {
        stop("Metabolomics design contrasts table contains duplicate friendly contrast labels.", call. = FALSE)
    }

    data.frame(
        contrasts = rawContrasts,
        friendly_names = friendlyNames,
        full_format = fullFormat,
        stringsAsFactors = FALSE
    )
}

metabDesignContrastTerms <- function(contrastString) {
    expression <- metabDesignContrastExpression(contrastString)
    expression <- gsub("\\s+", "", expression)
    terms <- regmatches(
        expression,
        gregexpr("[[:alpha:]_.][[:alnum:]_.]*", expression, perl = TRUE)
    )[[1]]
    unique(terms[nzchar(terms)])
}

validateMetabDesignContrasts <- function(designMatrix,
                                         contrastsTbl,
                                         formulaString = "~ 0 + group",
                                         allowEmpty = FALSE) {
    errors <- character()
    warnings <- character()
    modelMatrix <- NULL

    if (!metabDesignScalarString(formulaString)) {
        errors <- c(errors, "Metabolomics DA formula must be a non-empty scalar string.")
    }
    if (!is.data.frame(designMatrix) || nrow(designMatrix) == 0L) {
        errors <- c(errors, "Design matrix must be a non-empty data frame for contrast validation.")
    }

    normalisedContrasts <- tryCatch(
        normaliseMetabDesignContrastsTable(contrastsTbl, allowEmpty = allowEmpty),
        error = function(err) {
            errors <<- c(errors, conditionMessage(err))
            NULL
        }
    )

    if (length(errors) == 0L) {
        tryCatch(
            {
                formulaObj <- stats::as.formula(formulaString)
                modelFrame <- stats::model.frame(formulaObj, designMatrix, na.action = stats::na.pass)
                modelMatrix <- stats::model.matrix(formulaObj, modelFrame)
            },
            error = function(err) {
                errors <<- c(errors, sprintf("Metabolomics DA formula cannot produce a model matrix: %s", conditionMessage(err)))
            }
        )
    }

    if (length(errors) == 0L && !is.null(normalisedContrasts) && nrow(normalisedContrasts) > 0L) {
        modelTerms <- colnames(modelMatrix)
        for (idx in seq_len(nrow(normalisedContrasts))) {
            terms <- metabDesignContrastTerms(normalisedContrasts$contrasts[[idx]])
            missingTerms <- setdiff(terms, modelTerms)
            if (length(terms) == 0L) {
                errors <- c(
                    errors,
                    sprintf("Metabolomics DA contrast '%s' has no model terms.", normalisedContrasts$contrasts[[idx]])
                )
            }
            if (length(missingTerms) > 0L) {
                errors <- c(
                    errors,
                    sprintf(
                        "Metabolomics DA contrast '%s' references terms absent from the model matrix: %s",
                        normalisedContrasts$contrasts[[idx]],
                        paste(missingTerms, collapse = ", ")
                    )
                )
            }
        }
    }

    list(
        valid = length(errors) == 0L,
        errors = unique(errors),
        warnings = unique(warnings),
        contrasts_tbl = normalisedContrasts,
        model_terms = if (!is.null(modelMatrix)) colnames(modelMatrix) else character(),
        model_matrix = modelMatrix
    )
}

validateMetabDesignDaPreflight <- function(designMatrix,
                                           assayList,
                                           contrastsTbl = NULL,
                                           formulaString = "~ 0 + group",
                                           columnMapping = NULL,
                                           sampleIdCol = "Run",
                                           groupCol = "group",
                                           requireContrasts = TRUE) {
    alignment <- validateMetabDesignAlignment(
        designMatrix = designMatrix,
        assayList = assayList,
        columnMapping = columnMapping,
        sampleIdCol = sampleIdCol,
        groupCol = groupCol
    )

    contrastValidation <- validateMetabDesignContrasts(
        designMatrix = designMatrix,
        contrastsTbl = contrastsTbl,
        formulaString = formulaString,
        allowEmpty = !isTRUE(requireContrasts)
    )

    errors <- c(alignment$errors, contrastValidation$errors)
    warnings <- c(alignment$warnings, contrastValidation$warnings)

    list(
        valid = length(errors) == 0L,
        errors = unique(errors),
        warnings = unique(warnings),
        alignment = alignment,
        contrast_validation = contrastValidation,
        model_terms = contrastValidation$model_terms
    )
}
