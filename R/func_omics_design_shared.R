# KIT-103 W2: the design core of lipid, metab and trx (kit-103 scoping C1 to C5): the design
# validation, the formula, the sample discovery and initial design with its run order, the design
# edits and the save payload. Pure functions over tables and the omic's GUI configuration
# (`omicGuiConfig()$design`); the design module and builder (KIT-103 W4) keep the reactive state,
# the inputs and the notifications. Moved, with today's semantics, from mod_lipid_design_validation.R
# and mod_metab_design_validation.R (C1, C2), mod_lipid_design_builder_state_helpers.R and the metab
# builder server's closure (C3, C5) and mod_lipid_design_builder_action_helpers.R (C4), as at
# 70af7bb4. The commit and the design import are in R/utils_workflow_omics_design.R. The
# `DEBUG66` console messages are not kept (KIT-103 D5). Recorded divergence from conventions/R.md:
# the per-omic validation labels and non-sample patterns are held here by omic type
# (`omicsDesignRules()`), as the ticket keeps them out of the configuration; `shiny::req` is kept in
# the bulk rename's range mode, where today's builder stops its observer silently.

# ------------------------------------------------------------------------------------------------
# The per-omic rules.
# ------------------------------------------------------------------------------------------------

#' The design rules of one omic: the label of the validation messages, the column-mapping id key,
#' the validator's non-sample column patterns (both lists exactly as today's lipid and metab
#' validators) and the builder's annotation-name pattern (the one entry in which the lipid and
#' metab builders' sample discovery differs). Trx runs the lipid builder, so its builder rules are
#' lipid's; it validates its design with `trxDesignCheck()`, not with these patterns.
#' @param omic_type `lipidomics`, `metabolomics` or `transcriptomics`.
#' @return A list: `label`, `id_key`, `non_sample_patterns`, `builder_name_pattern`.
#' @noRd
omicsDesignRules <- function(omic_type) {
    lipid <- list(label = "Lipidomics", id_key = "lipid_id_col", builder_name_pattern = "^Lipid.name$",
        non_sample_patterns = c("^Alignment.ID$", "^Alignment ID$", "^Peak ID$", "^Average.Rt", "^Average.Mz",
            "^RT \\(min\\)$", "^Precursor m/z$", "^Lipid.name$", "^Lipid name$", "^LipidName$",
            "^LipidClass$", "^Name$", "^Annotation$", "^Adduct.type$", "^Adduct$", "^Ion mode$",
            "^Post.curation.result$", "^Fill.%$", "^MS/MS.spectrum$", "^Reference.RT$",
            "^Reference.m/z$", "^Formula$", "^Ontology$", "^INCHIKEY$", "^SMILES$", "^Annotation.tag",
            "^RT.matched$", "^m/z.matched$", "^MS/MS.matched$", "^Comment$", "^Manually.modified$",
            "^Total.score$", "^Total Score$", "^RT.similarity$", "^Dot.product$",
            "^Reverse.dot.product$", "^Fragment.presence", "^S/N.average$", "^Spectrum.reference",
            "^MS1.isotopic.spectrum$", "^Class$", "^Category$", "^Subclass$", "^Main class$"))
    metab <- list(label = "Metabolomics", id_key = "metabolite_id_col",
        builder_name_pattern = "^Metabolite.name$",
        non_sample_patterns = c("^Alignment.ID$", "^Average.Rt", "^Average.Mz", "^Peak ID$", "^RT \\(min\\)$",
            "^Precursor m/z$", "^Metabolite.name$", "^Metabolite name$", "^Adduct.type$", "^Adduct$",
            "^Post.curation.result$", "^Fill.%$", "^MS/MS.spectrum$", "^Reference.RT$",
            "^Reference.m/z$", "^Formula$", "^Ontology$", "^INCHIKEY$", "^SMILES$", "^Annotation",
            "^Annotation.tag", "^RT.matched$", "^m/z.matched$", "^MS/MS.matched$", "^Comment$",
            "^Manually.modified$", "^Total.score$", "^Total Score$", "^RT.similarity$",
            "^Dot.product$", "^Reverse.dot.product$", "^Fragment.presence", "^S/N.average$",
            "^Spectrum.reference", "^MS1.isotopic.spectrum$", "^Ion mode$"))
    switch(omic_type,
        lipidomics = lipid,
        metabolomics = metab,
        transcriptomics = utils::modifyList(lipid, list(label = "Transcriptomics")),
        cli::cli_abort("no design rules for omic {.val {omic_type}}.", class = "omics_design_config_invalid"))
}

#' Whether a value is one non-empty, non-missing string (blank text is empty).
#' @noRd
omicsDesignScalarString <- function(value) {
    is.character(value) && length(value) == 1L && !is.na(value) && nzchar(trimws(value))
}

# ------------------------------------------------------------------------------------------------
# C2: the formula.
# ------------------------------------------------------------------------------------------------

#' C2. The DA formula of a design: the configuration's `deAnalysisParameters$formula_string` when it
#' is one non-empty string, else `default` (lipid's `lipidDesignResolveFormulaString`; metab's inline
#' rule, mod_metab_design.R:225-228).
#' @param configList The workflow configuration list.
#' @param default The omic's `design$formula_default`.
#' @return One string.
#' @noRd
resolveOmicsDesignFormula <- function(configList = NULL, default = "~ 0 + group") {
    formulaString <- NULL
    if (is.list(configList) && is.list(configList[["deAnalysisParameters"]])) {
        formulaString <- configList[["deAnalysisParameters"]][["formula_string"]]
    }
    if (omicsDesignScalarString(formulaString)) formulaString else default
}

# ------------------------------------------------------------------------------------------------
# C1: the validation.
# ------------------------------------------------------------------------------------------------

#' The sample columns of one assay for the validator: the mapping's `sample_columns` when set, else
#' the numeric columns that are neither the id or annotation column nor a non-sample pattern.
#' @noRd
omicsDesignSampleColumnsForAssay <- function(assayData, columnMapping = NULL, rules) {
    if (!is.data.frame(assayData)) return(character())
    if (!is.null(columnMapping$sample_columns) && length(columnMapping$sample_columns) > 0L) {
        return(as.character(columnMapping$sample_columns))
    }
    sampleCols <- names(assayData)[vapply(assayData, is.numeric, logical(1))]
    excludeCols <- c(columnMapping[[rules$id_key]], columnMapping$annotation_col)
    excludeCols <- excludeCols[!is.na(excludeCols) & nzchar(excludeCols)]
    sampleCols <- setdiff(sampleCols, excludeCols)
    for (pattern in rules$non_sample_patterns) {
        sampleCols <- sampleCols[!grepl(pattern, sampleCols, ignore.case = TRUE)]
    }
    sampleCols
}

#' One assay's alignment with the design's samples: its errors, warnings and sample columns.
#' @noRd
omicsDesignAssaySummary <- function(assayName, assayData, columnMapping, designSamples, rules) {
    if (!is.data.frame(assayData)) {
        return(list(assay_name = assayName, sample_columns = character(),
            errors = sprintf("Assay '%s' is not a data frame.", assayName), warnings = character()))
    }
    errors <- character()
    warnings <- character()
    sampleCols <- omicsDesignSampleColumnsForAssay(assayData, columnMapping, rules)
    missingMappedCols <- setdiff(sampleCols, names(assayData))
    presentSampleCols <- intersect(sampleCols, names(assayData))
    inferenceMapping <- columnMapping %||% list()
    inferenceMapping$sample_columns <- NULL
    inferredSampleCols <- omicsDesignSampleColumnsForAssay(assayData, inferenceMapping, rules)
    extraAssaySampleCols <- setdiff(inferredSampleCols, sampleCols)
    listed <- function(values) paste(values, collapse = ", ")

    if (length(sampleCols) == 0L) {
        errors <- c(errors, sprintf("Assay '%s' has no sample columns available for design alignment.", assayName))
    }
    if (length(missingMappedCols) > 0L) {
        errors <- c(errors, sprintf("Assay '%s' is missing mapped sample column(s): %s", assayName,
            listed(missingMappedCols)))
    }
    if (length(extraAssaySampleCols) > 0L) {
        errors <- c(errors, sprintf("Assay '%s' has sample-like column(s) absent from column mapping/design: %s",
            assayName, listed(extraAssaySampleCols)))
    }
    if (anyDuplicated(presentSampleCols) > 0L) {
        errors <- c(errors, sprintf("Assay '%s' contains duplicate sample columns: %s", assayName,
            listed(unique(presentSampleCols[duplicated(presentSampleCols)]))))
    }
    missingInDesign <- setdiff(presentSampleCols, designSamples)
    extraInDesign <- setdiff(designSamples, presentSampleCols)
    caseCandidateCols <- unique(c(presentSampleCols, extraAssaySampleCols))
    caseMismatchDesign <- designSamples[tolower(designSamples) %in% tolower(caseCandidateCols) &
        !designSamples %in% caseCandidateCols]
    caseMismatchAssay <- caseCandidateCols[tolower(caseCandidateCols) %in% tolower(designSamples) &
        !caseCandidateCols %in% designSamples]
    if (length(missingInDesign) > 0L) {
        errors <- c(errors, sprintf("Assay '%s' has sample column(s) missing from design: %s", assayName,
            listed(missingInDesign)))
    }
    if (length(extraInDesign) > 0L) {
        errors <- c(errors, sprintf("Design has sample row(s) absent from assay '%s': %s", assayName,
            listed(extraInDesign)))
    }
    if (length(caseMismatchDesign) > 0L || length(caseMismatchAssay) > 0L) {
        errors <- c(errors, sprintf(
            "Assay '%s' has case-varied sample names between design and assay: design [%s], assay [%s]",
            assayName, listed(unique(caseMismatchDesign)), listed(unique(caseMismatchAssay))))
    }
    if (length(presentSampleCols) > 0L && setequal(presentSampleCols, designSamples) &&
        !identical(presentSampleCols, designSamples)) {
        warnings <- c(warnings, sprintf("Assay '%s' sample columns are reordered relative to the design matrix.",
            assayName))
    }
    list(assay_name = assayName, sample_columns = presentSampleCols, errors = errors, warnings = warnings)
}

#' C1. The design's alignment with the assays: a non-empty design and assay list, unique assay
#' names, unique non-blank sample ids, group and replicate assignments, and per assay the sample
#' columns against the design's samples.
#' @param rules The omic's `omicsDesignRules()`.
#' @return A list: `valid`, `errors`, `warnings` and, when both inputs are present,
#'   `design_samples` and `assay_summaries`.
#' @noRd
validateOmicsDesignAlignment <- function(designMatrix, assayList, rules, columnMapping = NULL,
                                         sampleIdCol = "Run", groupCol = "group",
                                         requireGroups = TRUE, requireReplicates = TRUE) {
    errors <- character()
    if (!is.data.frame(designMatrix) || nrow(designMatrix) == 0L) {
        errors <- c(errors, "Design matrix must be a non-empty data frame.")
    }
    if (is.null(assayList) || length(assayList) == 0L || !is.list(assayList)) {
        errors <- c(errors, paste(rules$label, "design requires a non-empty assay list."))
    }
    if (length(errors) > 0L) return(list(valid = FALSE, errors = unique(errors), warnings = character()))

    assayNames <- names(assayList)
    if (is.null(assayNames)) {
        assayNames <- paste0("assay_", seq_along(assayList))
    } else {
        unnamed <- !nzchar(assayNames)
        assayNames[unnamed] <- paste0("assay_", which(unnamed))
    }
    if (anyDuplicated(assayNames) > 0L) {
        errors <- c(errors, sprintf("%s design assay names must be unique: %s", rules$label,
            paste(unique(assayNames[duplicated(assayNames)]), collapse = ", ")))
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
            errors <- c(errors, sprintf("Design matrix contains duplicate sample rows: %s",
                paste(unique(designSamples[duplicated(designSamples)]), collapse = ", ")))
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
        omicsDesignAssaySummary(assayNames[[assayIdx]], assayList[[assayIdx]], columnMapping, designSamples, rules)
    })
    warnings <- character()
    for (summary in assaySummaries) {
        errors <- c(errors, summary$errors)
        warnings <- c(warnings, summary$warnings)
    }
    list(valid = length(errors) == 0L, errors = unique(errors), warnings = unique(warnings),
        design_samples = designSamples, assay_summaries = assaySummaries)
}

#' A contrast's expression: the text after a leading `name=`.
#' @noRd
omicsDesignContrastExpression <- function(contrastString) {
    sub("^[^=]+=", "", trimws(as.character(contrastString)))
}

#' The friendly name the validator gives a contrast without one: `group` prefixes stripped, blanks
#' removed, `-` as `_vs_`.
#' @noRd
omicsDesignFriendlyContrastName <- function(contrastString) {
    contrastString <- omicsDesignContrastExpression(contrastString)
    cleanString <- gsub("^group", "", contrastString)
    cleanString <- gsub("-group", "-", cleanString)
    cleanString <- gsub("\\s+", "", cleanString)
    gsub("-", "_vs_", cleanString)
}

#' C1. A contrasts table in the contract's columns (`contrasts`, `friendly_names`, `full_format`),
#' from any of `contrasts`, `contrast_string` or the first column; refuses an empty, duplicated or
#' ambiguously named contrast with today's message (a plain error).
#' @noRd
normaliseOmicsDesignContrastsTable <- function(contrastsTbl, rules, allowEmpty = FALSE) {
    if (is.null(contrastsTbl) || !is.data.frame(contrastsTbl) || nrow(contrastsTbl) == 0L) {
        if (isTRUE(allowEmpty)) {
            return(data.frame(contrasts = character(), friendly_names = character(), full_format = character(),
                stringsAsFactors = FALSE))
        }
        stop(paste(rules$label, "design contrasts table must be a non-empty data frame."), call. = FALSE)
    }
    rawInput <- if ("contrasts" %in% names(contrastsTbl)) {
        contrastsTbl$contrasts
    } else if ("contrast_string" %in% names(contrastsTbl)) {
        contrastsTbl$contrast_string
    } else {
        contrastsTbl[[1L]]
    }
    rawInput <- as.character(rawInput)
    rawContrasts <- vapply(rawInput, omicsDesignContrastExpression, character(1))
    if (any(is.na(rawContrasts) | !nzchar(trimws(rawContrasts)))) {
        stop(paste(rules$label, "design contrasts table contains an empty contrast."), call. = FALSE)
    }
    if (anyDuplicated(rawContrasts) > 0L) {
        stop(paste(rules$label, "design contrasts table contains duplicate contrast definitions."), call. = FALSE)
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
            friendlyNames[[idx]] <- if (!is.na(fullFormat[[idx]]) && grepl("=", fullFormat[[idx]], fixed = TRUE)) {
                sub("=.*$", "", fullFormat[[idx]])
            } else if (grepl("=", rawInput[[idx]], fixed = TRUE)) {
                sub("=.*$", "", rawInput[[idx]])
            } else {
                omicsDesignFriendlyContrastName(rawContrasts[[idx]])
            }
        }
        if (is.na(fullFormat[[idx]]) || !nzchar(trimws(fullFormat[[idx]]))) {
            fullFormat[[idx]] <- paste0(friendlyNames[[idx]], "=", rawContrasts[[idx]])
        }
    }
    if (anyDuplicated(friendlyNames) > 0L) {
        stop(paste(rules$label, "design contrasts table contains duplicate friendly contrast labels."), call. = FALSE)
    }
    data.frame(contrasts = rawContrasts, friendly_names = friendlyNames, full_format = fullFormat,
        stringsAsFactors = FALSE)
}

#' The model terms a contrast names.
#' @noRd
omicsDesignContrastTerms <- function(contrastString) {
    expression <- gsub("\\s+", "", omicsDesignContrastExpression(contrastString))
    terms <- regmatches(expression, gregexpr("[[:alpha:]_.][[:alnum:]_.]*", expression, perl = TRUE))[[1]]
    unique(terms[nzchar(terms)])
}

#' C1. The contrasts against the design's model matrix under the formula.
#' @return A list: `valid`, `errors`, `warnings`, `contrasts_tbl` (normalised, or NULL),
#'   `model_terms`, `model_matrix`.
#' @noRd
validateOmicsDesignContrasts <- function(designMatrix, contrastsTbl, rules, formulaString = "~ 0 + group",
                                         allowEmpty = FALSE) {
    errors <- character()
    modelMatrix <- NULL
    if (!omicsDesignScalarString(formulaString)) {
        errors <- c(errors, paste(rules$label, "DA formula must be a non-empty scalar string."))
    }
    if (!is.data.frame(designMatrix) || nrow(designMatrix) == 0L) {
        errors <- c(errors, "Design matrix must be a non-empty data frame for contrast validation.")
    }
    normalisedContrasts <- tryCatch(normaliseOmicsDesignContrastsTable(contrastsTbl, rules, allowEmpty = allowEmpty),
        error = function(err) {
            errors <<- c(errors, conditionMessage(err))
            NULL
        })
    if (length(errors) == 0L) {
        tryCatch({
            formulaObj <- stats::as.formula(formulaString)
            modelFrame <- stats::model.frame(formulaObj, designMatrix, na.action = stats::na.pass)
            modelMatrix <- stats::model.matrix(formulaObj, modelFrame)
        }, error = function(err) {
            errors <<- c(errors, sprintf("%s DA formula cannot produce a model matrix: %s", rules$label,
                conditionMessage(err)))
        })
    }
    if (length(errors) == 0L && !is.null(normalisedContrasts) && nrow(normalisedContrasts) > 0L) {
        modelTerms <- colnames(modelMatrix)
        for (idx in seq_len(nrow(normalisedContrasts))) {
            contrast <- normalisedContrasts$contrasts[[idx]]
            terms <- omicsDesignContrastTerms(contrast)
            missingTerms <- setdiff(terms, modelTerms)
            if (length(terms) == 0L) {
                errors <- c(errors, sprintf("%s DA contrast '%s' has no model terms.", rules$label, contrast))
            }
            if (length(missingTerms) > 0L) {
                errors <- c(errors, sprintf("%s DA contrast '%s' references terms absent from the model matrix: %s",
                    rules$label, contrast, paste(missingTerms, collapse = ", ")))
            }
        }
    }
    list(valid = length(errors) == 0L, errors = unique(errors), warnings = character(),
        contrasts_tbl = normalisedContrasts,
        model_terms = if (!is.null(modelMatrix)) colnames(modelMatrix) else character(),
        model_matrix = modelMatrix)
}

#' C1. The DA preflight: the alignment and the contrasts together.
#' @return A list: `valid`, `errors`, `warnings`, `alignment`, `contrast_validation`, `model_terms`.
#' @noRd
validateOmicsDesignDaPreflight <- function(designMatrix, assayList, rules, contrastsTbl = NULL,
                                           formulaString = "~ 0 + group", columnMapping = NULL,
                                           sampleIdCol = "Run", groupCol = "group", requireContrasts = TRUE) {
    alignment <- validateOmicsDesignAlignment(designMatrix = designMatrix, assayList = assayList, rules = rules,
        columnMapping = columnMapping, sampleIdCol = sampleIdCol, groupCol = groupCol)
    contrastValidation <- validateOmicsDesignContrasts(designMatrix = designMatrix, contrastsTbl = contrastsTbl,
        rules = rules, formulaString = formulaString, allowEmpty = !isTRUE(requireContrasts))
    errors <- c(alignment$errors, contrastValidation$errors)
    warnings <- c(alignment$warnings, contrastValidation$warnings)
    list(valid = length(errors) == 0L, errors = unique(errors), warnings = unique(warnings),
        alignment = alignment, contrast_validation = contrastValidation,
        model_terms = contrastValidation$model_terms)
}

#' C1. The save-time preflight of a builder result or an imported design (lipid's and metab's
#' commit, and lipid's import): the DA preflight under the result's formula, contrasts optional.
#' The shape of a configuration's `design$validate(results, columnMapping)`.
#' @param results A list with `design_matrix`, `data_cln`, `contrasts_tbl`, `config_list`.
#' @param omic_type `lipidomics` or `metabolomics`.
#' @noRd
validateOmicsDesignSave <- function(results, columnMapping, omic_type) {
    validateOmicsDesignDaPreflight(designMatrix = results$design_matrix, assayList = results$data_cln,
        rules = omicsDesignRules(omic_type), contrastsTbl = results$contrasts_tbl,
        formulaString = resolveOmicsDesignFormula(results$config_list), columnMapping = columnMapping,
        requireContrasts = FALSE)
}

# ------------------------------------------------------------------------------------------------
# C3: the sample discovery and the initial design.
# ------------------------------------------------------------------------------------------------

#' C3. The builder's sample columns of a wide assay list (today's lipid and metab builders): the
#' mapping's `sample_columns` when set; without a mapping, the first assay's numeric columns; else
#' its columns that are neither the id or annotation column nor a builder non-sample pattern.
#' @param config The omic's GUI configuration (its `shell$omic_type` selects the rules).
#' @return A character vector.
#' @noRd
omicsDesignSampleColumns <- function(assayList, colMap, config) {
    if (is.null(assayList) || length(assayList) == 0L) return(character(0))
    firstAssay <- assayList[[1]]
    if (!is.data.frame(firstAssay)) return(character(0))
    if (is.null(colMap)) return(names(firstAssay)[vapply(firstAssay, is.numeric, logical(1))])
    rules <- omicsDesignRules(config$shell$omic_type)
    excludeCols <- c(colMap[[rules$id_key]], colMap$annotation_col)
    excludeCols <- excludeCols[!is.na(excludeCols) & nzchar(excludeCols)]
    sampleCols <- setdiff(names(firstAssay), excludeCols)
    nonSamplePatterns <- c("^Alignment.ID$", "^Average.Rt", "^Average.Mz", rules$builder_name_pattern,
        "^Adduct.type$", "^Post.curation.result$", "^Fill.%$", "^MS/MS.spectrum$", "^Reference.RT$",
        "^Reference.m/z$", "^Formula$", "^Ontology$", "^INCHIKEY$", "^SMILES$", "^Annotation.tag",
        "^RT.matched$", "^m/z.matched$", "^MS/MS.matched$", "^Comment$", "^Manually.modified$",
        "^Total.score$", "^RT.similarity$", "^Dot.product$", "^Reverse.dot.product$",
        "^Fragment.presence", "^S/N.average$", "^Spectrum.reference", "^MS1.isotopic.spectrum$")
    for (pattern in nonSamplePatterns) {
        sampleCols <- sampleCols[!grepl(pattern, sampleCols, ignore.case = TRUE)]
    }
    if (!is.null(colMap$sample_columns) && length(colMap$sample_columns) > 0) sampleCols <- colMap$sample_columns
    sampleCols
}

#' C3. The builder's contrasts from an existing contrasts table: `groupA-groupB` is A against B
#' (`A.vs.B`); any other string is kept whole as its numerator.
#' @return A data frame: `contrast_name`, `numerator`, `denominator`.
#' @noRd
buildOmicsDesignContrastState <- function(existingContrasts) {
    if (is.null(existingContrasts) || !is.data.frame(existingContrasts) || nrow(existingContrasts) == 0) {
        return(data.frame(contrast_name = character(), numerator = character(), denominator = character(),
            stringsAsFactors = FALSE))
    }
    parsed <- lapply(existingContrasts$contrasts, function(contrastString) {
        parts <- strsplit(gsub("^group", "", contrastString), "-group")[[1]]
        if (length(parts) == 2) {
            return(list(contrast_name = paste0(parts[1], ".vs.", parts[2]), numerator = parts[1],
                denominator = parts[2]))
        }
        list(contrast_name = contrastString, numerator = contrastString, denominator = "")
    })
    data.frame(contrast_name = vapply(parsed, `[[`, character(1), "contrast_name"),
        numerator = vapply(parsed, `[[`, character(1), "numerator"),
        denominator = vapply(parsed, `[[`, character(1), "denominator"), stringsAsFactors = FALSE)
}

#' Resolve an optional existing value: a function is called (an error is NULL), else the value.
#' @noRd
omicsDesignExistingValue <- function(value) {
    if (is.null(value)) return(NULL)
    tryCatch(if (is.function(value)) value() else value, error = function(error) NULL)
}

#' C3. The builder's initial state over a wide assay list. An existing design with at least one
#' assigned group is used as it is, in its own (file) order, with its groups, the factor levels of
#' `factor1` and `factor2` and the existing contrasts. Otherwise a fresh design lists the samples in
#' the omic's run order (`config$design$run_order`, required, no default): the only place it is
#' applied.
#' @param assayList The assay tables (the builder's data).
#' @param configList The workflow configuration (its formula is the initial formula).
#' @param colMap The column mapping.
#' @param config The omic's GUI configuration.
#' @param existingDesignMatrix,existingContrasts Optional values, or functions returning them.
#' @return A list (`design_matrix`, `data_cln`, `sample_names`, `groups`, `factors`, `formula`,
#'   `contrasts`), or NULL when no sample column is found.
#' @noRd
buildOmicsDesignInitialState <- function(assayList, configList, colMap, config, existingDesignMatrix = NULL,
                                         existingContrasts = NULL) {
    runOrder <- config$design$run_order
    if (!is.function(runOrder)) {
        cli::cli_abort("the {config$shell$omic_type} design configuration has no run order.",
            class = "omics_design_config_invalid")
    }
    sampleNames <- omicsDesignSampleColumns(assayList, colMap, config)
    if (length(sampleNames) == 0) {
        logger::log_error(paste("No sample columns detected in", config$design$label_lower, "data"))
        return(NULL)
    }
    formulaString <- configList[["deAnalysisParameters"]][["formula_string"]]
    existing <- omicsDesignExistingValue(existingDesignMatrix)
    assigned <- function(values) values[!is.na(values) & values != ""]
    useExisting <- is.data.frame(existing) && nrow(existing) > 0 && "group" %in% names(existing) &&
        length(assigned(existing$group)) > 0
    # Read whenever given, as today (a reactive argument is then a dependency on either path).
    contrasts <- omicsDesignExistingValue(existingContrasts)
    if (useExisting) {
        factors <- character(0)
        for (column in intersect(c("factor1", "factor2"), names(existing))) {
            factors <- c(factors, unique(assigned(existing[[column]])))
        }
        return(list(design_matrix = existing, data_cln = assayList, sample_names = sampleNames,
            groups = unique(assigned(existing$group)), factors = unique(factors), formula = formulaString,
            contrasts = buildOmicsDesignContrastState(contrasts)))
    }
    designMatrix <- tibble::tibble(Run = runOrder(sampleNames), group = NA_character_, factor1 = NA_character_,
        factor2 = NA_character_, factor3 = NA_character_, batch = NA_character_, replicates = NA_integer_,
        tech_reps = NA_integer_)
    list(design_matrix = designMatrix, data_cln = assayList, sample_names = sampleNames,
        groups = unique(assigned(designMatrix$group)), factors = character(0), formula = formulaString,
        contrasts = buildOmicsDesignContrastState(NULL))
}

# ------------------------------------------------------------------------------------------------
# C4: the design edits (lipid's pure functions, the template).
# ------------------------------------------------------------------------------------------------

#' C4. Rename runs in the design (`renameMap`: old name -> new name), one name at a time.
#' @noRd
renameOmicsDesignRuns <- function(designMatrix, renameMap) {
    for (oldName in names(renameMap)) {
        designMatrix$Run[designMatrix$Run == oldName] <- unname(renameMap[[oldName]])
    }
    designMatrix
}

#' C4. Rename the sample columns of every assay.
#' @noRd
renameOmicsDesignAssayColumns <- function(assayList, renameMap) {
    for (assayName in names(assayList)) {
        assayColumns <- names(assayList[[assayName]])
        renamed <- renameMap[assayColumns]
        toRename <- !is.na(renamed)
        assayColumns[toRename] <- unname(renamed[toRename])
        names(assayList[[assayName]]) <- assayColumns
    }
    assayList
}

#' C4. Rename the tracked sample names.
#' @noRd
renameOmicsDesignTrackedSamples <- function(sampleNames, renameMap) {
    renamed <- renameMap[sampleNames]
    toRename <- !is.na(renamed)
    sampleNames[toRename] <- unname(renamed[toRename])
    sampleNames
}

#' C4. A rename applied to the design, the assays and the tracked names.
#' @return A list: `designMatrix`, `assayList`, `sampleNames`.
#' @noRd
applyOmicsDesignSampleRenameMap <- function(designMatrix, assayList, sampleNames, renameMap) {
    list(designMatrix = renameOmicsDesignRuns(designMatrix, renameMap),
        assayList = renameOmicsDesignAssayColumns(assayList, renameMap),
        sampleNames = renameOmicsDesignTrackedSamples(sampleNames, renameMap))
}

#' C4. The bulk rename map (`extract_experiment`): a character range, or the text before or after
#' the underscore. The range mode stops silently (`shiny::req`) without both bounds, as today.
#' @return A named character vector, old name -> new name.
#' @noRd
buildOmicsDesignBulkRenameMap <- function(samplesToTransform, transformMode, rangeStart = NULL, rangeEnd = NULL) {
    renamed <- vapply(samplesToTransform, function(sampleName) {
        if (transformMode == "range") {
            shiny::req(rangeStart, rangeEnd)
            return(extract_experiment(sampleName, mode = "range", start = rangeStart, end = rangeEnd))
        }
        if (transformMode == "before_underscore") return(extract_experiment(sampleName, mode = "start"))
        if (transformMode == "after_underscore") return(extract_experiment(sampleName, mode = "end"))
        stop("Unsupported transform mode.")
    }, character(1))
    stats::setNames(renamed, samplesToTransform)
}

#' C4. Add a factor name (trimmed) unless it is empty or already listed.
#' @noRd
appendOmicsDesignFactorName <- function(currentFactors, newFactorName) {
    trimmed <- trimws(newFactorName)
    if (trimmed == "" || trimmed %in% currentFactors) return(currentFactors)
    c(currentFactors, trimmed)
}

#' C4. The replicate numbers of an assignment: consecutive from the start, else NA.
#' @noRd
buildOmicsDesignReplicateNumbers <- function(selectedRuns, replicateStart) {
    if (!is.null(replicateStart) && !is.na(replicateStart)) {
        return(seq(replicateStart, length.out = length(selectedRuns)))
    }
    NA_integer_
}

#' C4. The group name: the chosen factor levels joined by `_`, NA when none.
#' @noRd
buildOmicsDesignGroupName <- function(factor1Selection, factor2Selection = NULL, factor3Selection = NULL) {
    chosen <- function(value) !is.null(value) && value != ""
    factorParts <- c()
    if (chosen(factor1Selection)) factorParts <- c(factorParts, factor1Selection)
    if (chosen(factor2Selection)) factorParts <- c(factorParts, factor2Selection)
    if (chosen(factor3Selection)) factorParts <- c(factorParts, factor3Selection)
    if (length(factorParts) == 0) return(NA_character_)
    paste(factorParts, collapse = "_")
}

#' C4. The assigned groups, in design order.
#' @noRd
listOmicsDesignAssignedGroups <- function(designMatrix) {
    unique(designMatrix$group[!is.na(designMatrix$group) & designMatrix$group != ""])
}

#' C4. Assign factor levels, group and replicate numbers to the selected runs.
#' @noRd
applyOmicsDesignMetadataAssignment <- function(designMatrix, selectedRuns, factor1Selection, factor2Selection = NULL,
                                               factor3Selection = NULL, replicateStart = NULL) {
    selected <- which(designMatrix$Run %in% selectedRuns)
    designMatrix$factor1[selected] <- factor1Selection
    designMatrix$factor2[selected] <- factor2Selection
    designMatrix$factor3[selected] <- if (!is.null(factor3Selection) && factor3Selection != "") {
        factor3Selection
    } else {
        NA_character_
    }
    designMatrix$replicates[selected] <- buildOmicsDesignReplicateNumbers(selectedRuns, replicateStart)
    designMatrix$group[selected] <- buildOmicsDesignGroupName(factor1Selection, factor2Selection, factor3Selection)
    designMatrix
}

#' C4. Mark runs of one group as technical replicates of one biological replicate (`lowest`,
#' `first` or a `manual` number), renumbering the group's later replicates.
#' @return A list: `ok`, `type`, `message` (the notification) and `designMatrix`.
#' @noRd
applyOmicsDesignTechRepAssignment <- function(designMatrix, techRepSamples, assignmentMode,
                                              manualReplicateNumber = NULL) {
    refuse <- function(message) list(ok = FALSE, type = "warning", message = message, designMatrix = designMatrix)
    if (length(techRepSamples) < 2) return(refuse("Please select at least two samples."))
    selected <- which(designMatrix$Run %in% techRepSamples)
    selectedGroups <- unique(designMatrix$group[selected])
    if (length(selectedGroups) > 1 || any(is.na(selectedGroups))) {
        return(refuse("All selected samples must belong to the same group."))
    }
    baseReplicateNumber <- switch(assignmentMode,
        lowest = min(designMatrix$replicates[selected], na.rm = TRUE),
        first = designMatrix$replicates[designMatrix$Run == techRepSamples[1]],
        manual = manualReplicateNumber,
        stop("Unsupported technical replicate assignment mode."))
    originalReplicateNumbers <- unique(designMatrix$replicates[selected])
    consolidated <- length(originalReplicateNumbers) - 1
    highest <- max(originalReplicateNumbers, na.rm = TRUE)
    toAdjust <- which(designMatrix$group == selectedGroups[1] & !is.na(designMatrix$group) &
        designMatrix$replicates > highest & !designMatrix$Run %in% techRepSamples)
    if (length(toAdjust) > 0) designMatrix$replicates[toAdjust] <- designMatrix$replicates[toAdjust] - consolidated
    designMatrix$replicates[selected] <- baseReplicateNumber
    techRepMapping <- stats::setNames(seq_along(techRepSamples), techRepSamples)
    designMatrix$tech_reps[selected] <- techRepMapping[designMatrix$Run[selected]]
    list(ok = TRUE, type = "message",
        message = paste("Assigned", length(techRepSamples), "samples as technical replicates"),
        designMatrix = designMatrix)
}

#' C4. Add the contrast `group1` against `group2` (`group1.vs.group2`) unless it is incomplete,
#' self-referential or already listed.
#' @noRd
appendOmicsDesignContrast <- function(currentContrasts, group1, group2) {
    if (is.null(group1) || is.null(group2) || group1 == "" || group2 == "" || group1 == group2) {
        return(currentContrasts)
    }
    contrastName <- paste0(group1, ".vs.", group2)
    if (contrastName %in% currentContrasts$contrast_name) return(currentContrasts)
    rbind(currentContrasts, data.frame(contrast_name = contrastName, numerator = group1, denominator = group2,
        stringsAsFactors = FALSE))
}

#' C4. Add runs to the removed runs.
#' @noRd
appendOmicsDesignRemovedSamples <- function(currentRemovedSamples, selectedSamples) {
    unique(c(currentRemovedSamples, selectedSamples))
}

# ------------------------------------------------------------------------------------------------
# C5: the save payload.
# ------------------------------------------------------------------------------------------------

#' C5. The saved design: the assigned runs that are not removed, in the design's order, with
#' `tech_rep_group` (`<group>_<replicates>`).
#' @noRd
buildOmicsDesignSaveDesignMatrix <- function(designMatrix, removedSamples = character(0)) {
    designMatrix |>
        dplyr::filter(!is.na(group) & group != "") |>
        dplyr::filter(!Run %in% removedSamples) |>
        dplyr::mutate(tech_rep_group = paste(group, replicates, sep = "_"))
}

#' C5. The saved assays: the id and annotation columns, then the assigned runs in the design's order
#' (the save's run order is the design's).
#' @param idKey The column mapping's id key (`design$id_mapping_key`).
#' @noRd
buildOmicsDesignSaveDataList <- function(dataCln, assignedSamples, colMap, idKey) {
    lapply(dataCln, function(assayDf) {
        keepCols <- c(colMap[[idKey]], colMap$annotation_col)
        keepCols <- keepCols[!is.na(keepCols) & nzchar(keepCols)]
        keepCols <- keepCols[keepCols %in% names(assayDf)]
        assayDf[, unique(c(keepCols, intersect(assignedSamples, names(assayDf)))), drop = FALSE]
    })
}

#' C5. The contract's contrasts table from the builder's contrasts: the friendly name is the
#' builder's `A.vs.B` with dots as `_` (`A_vs_B`); the contrast is `groupA-groupB` under a
#' `~ 0 + group` formula, else `A-B`. NULL without contrasts.
#' @return A data frame (`contrasts`, `friendly_names`, `full_format`), or NULL.
#' @noRd
buildOmicsDesignSaveContrastsTable <- function(contrastData, formulaString) {
    if (is.null(contrastData) || nrow(contrastData) == 0) return(NULL)
    groupPrefix <- grepl("~ *0 *\\+ *group", formulaString)
    friendlyNames <- gsub("\\.", "_", contrastData$contrast_name)
    contrastStrings <- if (groupPrefix) {
        paste0("group", contrastData$numerator, "-group", contrastData$denominator)
    } else {
        paste0(contrastData$numerator, "-", contrastData$denominator)
    }
    data.frame(contrasts = contrastStrings, friendly_names = friendlyNames,
        full_format = paste0(friendlyNames, "=", contrastStrings), stringsAsFactors = FALSE)
}

#' C5. The builder's save payload: the saved design, the assays subset to its runs in its order, the
#' contrasts table and the configuration with the formula. NULL when no run is assigned.
#' @return A list (`design_matrix`, `data_cln`, `contrasts_tbl`, `config_list`), or NULL.
#' @noRd
buildOmicsDesignSavePayload <- function(designMatrix, removedSamples = character(0), dataCln, colMap, contrastData,
                                        configList, formulaString, idKey) {
    designMatrixFinal <- buildOmicsDesignSaveDesignMatrix(designMatrix, removedSamples)
    if (nrow(designMatrixFinal) == 0) return(NULL)
    configList[["deAnalysisParameters"]][["formula_string"]] <- formulaString
    list(design_matrix = designMatrixFinal,
        data_cln = buildOmicsDesignSaveDataList(dataCln, designMatrixFinal$Run, colMap, idKey),
        contrasts_tbl = buildOmicsDesignSaveContrastsTable(contrastData, formulaString),
        config_list = configList)
}

# ------------------------------------------------------------------------------------------------
# The lipid wrappers (KIT-103 W4).
# ------------------------------------------------------------------------------------------------

#' D6. Lipid's validators under their names, as one-line wrappers over the shared family with the
#' lipid rules, for their callers outside the design tab (the lipid DA analysis and QC; KIT-104 and
#' KIT-106 retarget them). Signatures as mod_lipid_design_validation.R at 616e1147.
#' @noRd
validateLipidDesignAlignment <- function(designMatrix, assayList, columnMapping = NULL, sampleIdCol = "Run",
                                         groupCol = "group", requireGroups = TRUE, requireReplicates = TRUE) {
    validateOmicsDesignAlignment(designMatrix, assayList, omicsDesignRules("lipidomics"), columnMapping = columnMapping,
        sampleIdCol = sampleIdCol, groupCol = groupCol, requireGroups = requireGroups,
        requireReplicates = requireReplicates)
}

#' D6. See `validateLipidDesignAlignment`.
#' @noRd
validateLipidDesignContrasts <- function(designMatrix, contrastsTbl, formulaString = "~ 0 + group", allowEmpty = FALSE) {
    validateOmicsDesignContrasts(designMatrix, contrastsTbl, omicsDesignRules("lipidomics"),
        formulaString = formulaString, allowEmpty = allowEmpty)
}

#' D6. See `validateLipidDesignAlignment`.
#' @noRd
validateLipidDesignDaPreflight <- function(designMatrix, assayList, contrastsTbl = NULL, formulaString = "~ 0 + group",
                                           columnMapping = NULL, sampleIdCol = "Run", groupCol = "group",
                                           requireContrasts = TRUE) {
    validateOmicsDesignDaPreflight(designMatrix, assayList, omicsDesignRules("lipidomics"), contrastsTbl = contrastsTbl,
        formulaString = formulaString, columnMapping = columnMapping, sampleIdCol = sampleIdCol, groupCol = groupCol,
        requireContrasts = requireContrasts)
}

#' The configuration the lipid design wrappers run the shared design tab with: lipid's, or the
#' profile's omic (trx, `trxOmicProfile()`). Until KIT-103 W7 wires the configuration's design
#' function fields, the missing ones are filled with the core's (the values W7 wires; once they are
#' wired this changes nothing): lipid `validate` (the save preflight), `build_object`
#' (`lipidDesignObject`), the import's `reader` and `record_design_file`; trx `group_levels`,
#' `validate` and `build_object` (func_trx_design.R).
#' @param profile An omic profile, or NULL for lipidomics.
#' @return An omic GUI configuration.
#' @noRd
lipidDesignGuiConfig <- function(profile = NULL) {
    omic_type <- profile$omic_type %||% "lipidomics"
    config <- omicGuiConfig(omic_type)
    design <- config$design
    if (identical(omic_type, "transcriptomics")) {
        design$group_levels <- design$group_levels %||% trxDesignGroupLevels
        design$validate <- design$validate %||% trxDesignCheck
        design$build_object <- design$build_object %||% trxDesignObject
    } else {
        design$validate <- design$validate %||%
            function(results, columnMapping) validateOmicsDesignSave(results, columnMapping, "lipidomics")
        design$build_object <- design$build_object %||% lipidDesignObject
        design$design_import$reader <- design$design_import$reader %||% readOmicsDesignImport
        design$design_import$record_design_file <- design$design_import$record_design_file %||% TRUE
    }
    config$design <- design
    config
}

# ------------------------------------------------------------------------------------------------
# The metab wrappers (KIT-103 W5).
# ------------------------------------------------------------------------------------------------

#' D6. Metab's validators under their names, as one-line wrappers over the shared family with the
#' metab rules, for their callers outside the design tab (the metab DA observer; KIT-106 retargets
#' it). Signatures as mod_metab_design_validation.R at de2d891c.
#' @noRd
validateMetabDesignAlignment <- function(designMatrix, assayList, columnMapping = NULL, sampleIdCol = "Run",
                                         groupCol = "group", requireGroups = TRUE, requireReplicates = TRUE) {
    validateOmicsDesignAlignment(designMatrix, assayList, omicsDesignRules("metabolomics"),
        columnMapping = columnMapping, sampleIdCol = sampleIdCol, groupCol = groupCol, requireGroups = requireGroups,
        requireReplicates = requireReplicates)
}

#' D6. See `validateMetabDesignAlignment`.
#' @noRd
validateMetabDesignContrasts <- function(designMatrix, contrastsTbl, formulaString = "~ 0 + group", allowEmpty = FALSE) {
    validateOmicsDesignContrasts(designMatrix, contrastsTbl, omicsDesignRules("metabolomics"),
        formulaString = formulaString, allowEmpty = allowEmpty)
}

#' D6. See `validateMetabDesignAlignment`.
#' @noRd
validateMetabDesignDaPreflight <- function(designMatrix, assayList, contrastsTbl = NULL, formulaString = "~ 0 + group",
                                           columnMapping = NULL, sampleIdCol = "Run", groupCol = "group",
                                           requireContrasts = TRUE) {
    validateOmicsDesignDaPreflight(designMatrix, assayList, omicsDesignRules("metabolomics"),
        contrastsTbl = contrastsTbl, formulaString = formulaString, columnMapping = columnMapping,
        sampleIdCol = sampleIdCol, groupCol = groupCol, requireContrasts = requireContrasts)
}

#' The configuration the metab design wrappers run the shared design tab with. Until KIT-103 W7
#' wires the configuration's design function fields, the missing ones are filled with the core's
#' (the values W7 wires; once they are wired this changes nothing): `validate` (the save preflight),
#' `build_object` (`metabDesignObject`), the import's `reader`, `record_design_file` (FALSE: the
#' provenance takes the owner's source directory, P-14) and `check_files_first` (TRUE: metab checks
#' the import's files before it closes the modal, P-11).
#' @return An omic GUI configuration.
#' @noRd
metabDesignGuiConfig <- function() {
    config <- omicGuiConfig("metabolomics")
    design <- config$design
    design$validate <- design$validate %||%
        function(results, columnMapping) validateOmicsDesignSave(results, columnMapping, "metabolomics")
    design$build_object <- design$build_object %||% metabDesignObject
    design$design_import$reader <- design$design_import$reader %||% readOmicsDesignImport
    design$design_import$record_design_file <- design$design_import$record_design_file %||% FALSE
    design$design_import$check_files_first <- design$design_import$check_files_first %||% TRUE
    config$design <- design
    config
}
