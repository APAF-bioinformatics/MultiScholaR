lipidDesignScalarString <- function(value) {
    is.character(value) && length(value) == 1L && !is.na(value) && nzchar(trimws(value))
}

lipidDesignResolveFormulaString <- function(configList = NULL, default = "~ 0 + group") {
    formulaString <- NULL
    if (is.list(configList) && is.list(configList[["deAnalysisParameters"]])) {
        formulaString <- configList[["deAnalysisParameters"]][["formula_string"]]
    }
    if (lipidDesignScalarString(formulaString)) {
        formulaString
    } else {
        default
    }
}

lipidDesignNonSamplePatterns <- function() {
    c(
        "^Alignment.ID$",
        "^Alignment ID$",
        "^Peak ID$",
        "^Average.Rt",
        "^Average.Mz",
        "^RT \\(min\\)$",
        "^Precursor m/z$",
        "^Lipid.name$",
        "^Lipid name$",
        "^LipidName$",
        "^LipidClass$",
        "^Name$",
        "^Annotation$",
        "^Adduct.type$",
        "^Adduct$",
        "^Ion mode$",
        "^Post.curation.result$",
        "^Fill.%$",
        "^MS/MS.spectrum$",
        "^Reference.RT$",
        "^Reference.m/z$",
        "^Formula$",
        "^Ontology$",
        "^INCHIKEY$",
        "^SMILES$",
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
        "^Class$",
        "^Category$",
        "^Subclass$",
        "^Main class$"
    )
}

lipidDesignSampleColumnsForAssay <- function(assayData, columnMapping = NULL) {
    if (!is.data.frame(assayData)) {
        return(character())
    }

    if (!is.null(columnMapping$sample_columns) && length(columnMapping$sample_columns) > 0L) {
        return(as.character(columnMapping$sample_columns))
    }

    sampleCols <- names(assayData)[vapply(assayData, is.numeric, logical(1))]
    excludeCols <- c(columnMapping$lipid_id_col, columnMapping$annotation_col)
    excludeCols <- excludeCols[!is.na(excludeCols) & nzchar(excludeCols)]
    sampleCols <- setdiff(sampleCols, excludeCols)

    for (pattern in lipidDesignNonSamplePatterns()) {
        sampleCols <- sampleCols[!grepl(pattern, sampleCols, ignore.case = TRUE)]
    }

    sampleCols
}

validateLipidDesignAlignment <- function(designMatrix,
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
        errors <- c(errors, "Lipidomics design requires a non-empty assay list.")
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
            "Lipidomics design assay names must be unique: %s",
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

        sampleCols <- lipidDesignSampleColumnsForAssay(assayData, columnMapping)
        missingMappedCols <- setdiff(sampleCols, names(assayData))
        presentSampleCols <- intersect(sampleCols, names(assayData))
        inferenceMapping <- columnMapping
        if (is.null(inferenceMapping)) {
            inferenceMapping <- list()
        }
        inferenceMapping$sample_columns <- NULL
        inferredSampleCols <- lipidDesignSampleColumnsForAssay(assayData, inferenceMapping)
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

lipidDesignContrastExpression <- function(contrastString) {
    sub("^[^=]+=", "", trimws(as.character(contrastString)))
}

lipidDesignFriendlyContrastName <- function(contrastString) {
    contrastString <- lipidDesignContrastExpression(contrastString)
    cleanString <- gsub("^group", "", contrastString)
    cleanString <- gsub("-group", "-", cleanString)
    cleanString <- gsub("\\s+", "", cleanString)
    gsub("-", "_vs_", cleanString)
}

normaliseLipidDesignContrastsTable <- function(contrastsTbl, allowEmpty = FALSE) {
    if (is.null(contrastsTbl) || !is.data.frame(contrastsTbl) || nrow(contrastsTbl) == 0L) {
        if (isTRUE(allowEmpty)) {
            return(data.frame(
                contrasts = character(),
                friendly_names = character(),
                full_format = character(),
                stringsAsFactors = FALSE
            ))
        }
        stop("Lipidomics design contrasts table must be a non-empty data frame.", call. = FALSE)
    }

    rawInput <- if ("contrasts" %in% names(contrastsTbl)) {
        contrastsTbl$contrasts
    } else if ("contrast_string" %in% names(contrastsTbl)) {
        contrastsTbl$contrast_string
    } else {
        contrastsTbl[[1L]]
    }
    rawInput <- as.character(rawInput)
    rawContrasts <- vapply(rawInput, lipidDesignContrastExpression, character(1))

    if (any(is.na(rawContrasts) | !nzchar(trimws(rawContrasts)))) {
        stop("Lipidomics design contrasts table contains an empty contrast.", call. = FALSE)
    }
    if (anyDuplicated(rawContrasts) > 0L) {
        stop("Lipidomics design contrasts table contains duplicate contrast definitions.", call. = FALSE)
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
                friendlyNames[[idx]] <- lipidDesignFriendlyContrastName(rawContrasts[[idx]])
            }
        }
        if (is.na(fullFormat[[idx]]) || !nzchar(trimws(fullFormat[[idx]]))) {
            fullFormat[[idx]] <- paste0(friendlyNames[[idx]], "=", rawContrasts[[idx]])
        }
    }

    if (anyDuplicated(friendlyNames) > 0L) {
        stop("Lipidomics design contrasts table contains duplicate friendly contrast labels.", call. = FALSE)
    }

    data.frame(
        contrasts = rawContrasts,
        friendly_names = friendlyNames,
        full_format = fullFormat,
        stringsAsFactors = FALSE
    )
}

lipidDesignContrastTerms <- function(contrastString) {
    expression <- lipidDesignContrastExpression(contrastString)
    expression <- gsub("\\s+", "", expression)
    terms <- regmatches(
        expression,
        gregexpr("[[:alpha:]_.][[:alnum:]_.]*", expression, perl = TRUE)
    )[[1]]
    unique(terms[nzchar(terms)])
}

validateLipidDesignContrasts <- function(designMatrix,
                                         contrastsTbl,
                                         formulaString = "~ 0 + group",
                                         allowEmpty = FALSE) {
    errors <- character()
    warnings <- character()
    modelMatrix <- NULL

    if (!lipidDesignScalarString(formulaString)) {
        errors <- c(errors, "Lipidomics DA formula must be a non-empty scalar string.")
    }
    if (!is.data.frame(designMatrix) || nrow(designMatrix) == 0L) {
        errors <- c(errors, "Design matrix must be a non-empty data frame for contrast validation.")
    }

    normalisedContrasts <- tryCatch(
        normaliseLipidDesignContrastsTable(contrastsTbl, allowEmpty = allowEmpty),
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
                errors <<- c(errors, sprintf("Lipidomics DA formula cannot produce a model matrix: %s", conditionMessage(err)))
            }
        )
    }

    if (length(errors) == 0L && !is.null(normalisedContrasts) && nrow(normalisedContrasts) > 0L) {
        modelTerms <- colnames(modelMatrix)
        for (idx in seq_len(nrow(normalisedContrasts))) {
            terms <- lipidDesignContrastTerms(normalisedContrasts$contrasts[[idx]])
            missingTerms <- setdiff(terms, modelTerms)
            if (length(terms) == 0L) {
                errors <- c(
                    errors,
                    sprintf("Lipidomics DA contrast '%s' has no model terms.", normalisedContrasts$contrasts[[idx]])
                )
            }
            if (length(missingTerms) > 0L) {
                errors <- c(
                    errors,
                    sprintf(
                        "Lipidomics DA contrast '%s' references terms absent from the model matrix: %s",
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

validateLipidDesignDaPreflight <- function(designMatrix,
                                           assayList,
                                           contrastsTbl = NULL,
                                           formulaString = "~ 0 + group",
                                           columnMapping = NULL,
                                           sampleIdCol = "Run",
                                           groupCol = "group",
                                           requireContrasts = TRUE) {
    alignment <- validateLipidDesignAlignment(
        designMatrix = designMatrix,
        assayList = assayList,
        columnMapping = columnMapping,
        sampleIdCol = sampleIdCol,
        groupCol = groupCol
    )

    contrastValidation <- validateLipidDesignContrasts(
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
