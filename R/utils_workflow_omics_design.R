# KIT-103 W2: the design commit and the design import of lipid, metab and trx (kit-103 scoping C8,
# C10). `commitOmicsDesign()` takes a builder result to the bus, the seven artifacts, the design
# object, the owner binding, the provenance and the saved state; `checkOmicsDesignImportFiles()`,
# `readOmicsDesignImport()` and `commitOmicsDesignImport()` take an exported design directory to the
# same. Each runs today's steps in today's order, configured by the omic's GUI configuration
# (`omicGuiConfig()$design`), and signals a classed condition with today's message. The GUI keeps
# the modal, the notifications, the QC trigger and the tab status (KIT-102's split). Moved from the
# lipid builder commit and import shell (mod_lipid_design_builder_helpers.R,
# mod_lipid_design_import_helpers.R) and metab's commit and import chain (mod_metab_design.R,
# mod_metab_design_import_helpers.R) as at 70af7bb4; the `DEBUG66` console messages are not kept
# (KIT-103 D5).

#' Refuse a design with today's message, before anything is written (`omics_design_refused`).
#' @noRd
omicsDesignRefuse <- function(message) {
    rlang::abort(message, class = "omics_design_refused")
}

# ------------------------------------------------------------------------------------------------
# C8: the commit.
# ------------------------------------------------------------------------------------------------

#' C8. The design's artifacts in the source directory, in today's order: `design_matrix.tab`,
#' `contrast_strings.tab` (with contrasts), one `data_cln_<assay>.tab` per assay,
#' `assay_manifest.txt`, `column_mapping.json` (with a mapping), `manifest.json` and `config.ini`
#' (with a configuration; a failed write is a warning, as today).
#' @param results A builder result (`design_matrix`, `data_cln` as tables, `contrasts_tbl`).
#' @return The written paths, by artifact.
#' @noRd
writeOmicsDesignArtifacts <- function(sourceDir, results, columnMapping, configList) {
    paths <- list(design_matrix = file.path(sourceDir, "design_matrix.tab"))
    logger::log_info(paste("Writing design matrix to:", paths$design_matrix))
    utils::write.table(results$design_matrix, file = paths$design_matrix, sep = "\t", row.names = FALSE,
        quote = FALSE)
    if (!is.null(results$contrasts_tbl) && nrow(results$contrasts_tbl) > 0) {
        paths$contrasts <- file.path(sourceDir, "contrast_strings.tab")
        logger::log_info(paste("Writing contrasts to:", paths$contrasts))
        writeLines(results$contrasts_tbl$contrasts, paths$contrasts)
    }
    assayNames <- names(results$data_cln)
    paths$assays <- character()
    for (assayName in assayNames) {
        assayPath <- file.path(sourceDir, paste0("data_cln_", assayName, ".tab"))
        logger::log_info(paste("Writing assay data to:", assayPath))
        utils::write.table(results$data_cln[[assayName]], file = assayPath, sep = "\t", row.names = FALSE,
            quote = FALSE)
        paths$assays <- c(paths$assays, assayPath)
    }
    paths$assay_manifest <- file.path(sourceDir, "assay_manifest.txt")
    writeLines(assayNames, paths$assay_manifest)
    logger::log_info(sprintf("Saved assay manifest with %d assays: %s", length(assayNames),
        paste(assayNames, collapse = ", ")))
    if (!is.null(columnMapping)) {
        paths$column_mapping <- file.path(sourceDir, "column_mapping.json")
        jsonlite::write_json(columnMapping, paths$column_mapping, auto_unbox = TRUE)
    }
    paths$manifest <- file.path(sourceDir, "manifest.json")
    jsonlite::write_json(list(data_path = "assay_manifest.txt", design_matrix_path = "design_matrix.tab",
        contrast_strings_path = "contrast_strings.tab", column_mapping_path = "column_mapping.json"),
        paths$manifest, auto_unbox = TRUE, pretty = TRUE)
    if (!is.null(configList)) {
        configPath <- file.path(sourceDir, "config.ini")
        tryCatch({
            ini::write.ini(configList, configPath)
            paths$config <- configPath
        }, error = function(e) {
            logger::log_warn(paste("Could not save config.ini:", e$message))
        })
    }
    paths
}

#' C8. Commit a builder result, in the contract's order: the source directory, the omic's preflight
#' (`design$validate(results, columnMapping)`), the bus (design, owner binding at import, the held
#' tables, contrasts, configuration), the artifacts, the object (`design$build_object`), the owner
#' binding at design, the provenance (`design$record_provenance`) and the saved state. The QC
#' trigger, the tab status, the modal and the notifications stay in the GUI.
#' @param workflowData The workflow bus (its `column_mapping`, `data_format` and `state_manager`).
#' @param results The builder's result; its `data_cln` may be a held reference (KIT-075).
#' @param sourceDir The source directory (`experiment_paths$source_dir`).
#' @param config The omic's GUI configuration.
#' @return The object and the written paths, invisibly. A refusal (a missing source directory, a
#'   failed preflight) is `omics_design_refused`, signalled before the bus is written; a failure
#'   after it is `omics_design_save_failed`, with the original message.
#' @noRd
commitOmicsDesign <- function(workflowData, results, sourceDir, config) {
    design <- config$design
    route <- config$shell$routes[[1L]]
    logger::log_info(sprintf("Received results from %s design builder. Saving to workflow and disk.",
        design$label_lower))
    # KIT-075: the builder's result carries the tables' reference; they are read for this action.
    held <- results$data_cln
    results$data_cln <- resolveWorkflowOwnedValue(held, workflowData$state_manager)
    if (is.null(sourceDir) || !dir.exists(sourceDir)) {
        omicsDesignRefuse("Could not find source directory to save files.")
    }
    columnMapping <- workflowData$column_mapping
    preflight <- design$validate(results, columnMapping)
    if (!isTRUE(preflight$valid)) {
        omicsDesignRefuse(paste(sprintf("Invalid %s design:", design$label_lower),
            paste(preflight$errors, collapse = "; ")))
    }
    workflowData$design_matrix <- results$design_matrix
    bindWorkflowOmicsOwnerAtImport(workflowData, route, workflowData$data_format)
    setWorkflowOmicsBusTables(workflowData, "data_cln", held)
    workflowData$contrasts_tbl <- results$contrasts_tbl
    # KIT-030: contrasts_tbl and config_list live on workflow_data only (no global copies).
    workflowData$config_list <- results$config_list
    tryCatch({
        paths <- writeOmicsDesignArtifacts(sourceDir, results, columnMapping, workflowData$config_list)
        object <- design$build_object(results, columnMapping)
        # KIT-070: the session owner is bound to the omic's route here (the import wrote files only).
        bindWorkflowOwnerTypeAtDesign(workflowData, route, workflowData$data_format)
        # KIT-031/032: the design step's provenance (`design_builder`) and the workflow's inputs.
        design$record_provenance(workflowData$state_manager, object, source_dir = sourceDir)
        workflowData$state_manager$saveState(state_name = design$state_name, s4_data_object = object,
            config_object = results$config_list, description = design$state_description)
        logger::log_info(sprintf("Design save complete. Files saved to: %s", sourceDir))
        invisible(list(object = object, paths = paths))
    }, error = function(error) {
        rlang::abort(conditionMessage(error), class = "omics_design_save_failed", original = error)
    })
}

# ------------------------------------------------------------------------------------------------
# C10: the design import.
# ------------------------------------------------------------------------------------------------

#' The files of an exported design directory.
#' @noRd
omicsDesignImportFiles <- function(importPath) {
    list(design = file.path(importPath, "design_matrix.tab"),
        contrasts = file.path(importPath, "contrast_strings.tab"),
        manifest = file.path(importPath, "assay_manifest.txt"),
        column_mapping = file.path(importPath, "column_mapping.json"),
        config = file.path(importPath, "config.ini"))
}

#' C10. Whether an exported design directory has its design and assay manifest, with today's
#' message when not. The GUI shows the message; lipid closes the import modal first, metab after
#' the check (P-11 keeps both).
#' @return A list: `ok`, and `message` when not, else `files`.
#' @noRd
checkOmicsDesignImportFiles <- function(importPath) {
    files <- omicsDesignImportFiles(importPath)
    if (!file.exists(files$design)) {
        return(list(ok = FALSE, message = "Import failed: 'design_matrix.tab' not found in the selected directory."))
    }
    if (!file.exists(files$manifest)) {
        return(list(ok = FALSE, message = paste("Import failed: 'assay_manifest.txt' not found.",
            "Cannot determine which assay files to load.")))
    }
    list(ok = TRUE, files = files)
}

#' The column mapping of an import without `column_mapping.json`, inferred as today from the first
#' assay: the first id candidate present (else the first column), the first annotation candidate
#' present (else NULL), the design's runs among the columns, no internal-standard pattern.
#' @noRd
inferOmicsDesignImportMapping <- function(assayList, designMatrix, config) {
    columns <- names(assayList[[1]])
    candidates <- config$design$design_import$id_candidates
    first <- function(values) {
        for (candidate in values) if (candidate %in% columns) return(candidate)
        NULL
    }
    idColumn <- first(candidates$id) %||% columns[1]
    annotationColumn <- first(candidates$annotation)
    sampleColumns <- intersect(designMatrix$Run, columns)
    logger::log_info(sprintf("Inferred column mapping: id=%s, annotation=%s, %d samples", idColumn,
        annotationColumn, length(sampleColumns)))
    stats::setNames(list(idColumn, annotationColumn, sampleColumns, NA_character_),
        c(config$design$id_mapping_key, "annotation_col", "sample_columns", "is_pattern"))
}

#' The contrasts of an import's `contrast_strings.tab` (today's import rule): the friendly name
#' strips the `group` prefixes and writes `-` as `_vs_` (P-09 keeps it apart from the builder's).
#' @noRd
readOmicsDesignImportContrasts <- function(path) {
    if (!file.exists(path)) {
        logger::log_info("No contrast_strings.tab found.")
        return(NULL)
    }
    contrastStrings <- readLines(path)
    if (length(contrastStrings) == 0) return(NULL)
    cleanStrings <- gsub("-group", "-", gsub("^group", "", contrastStrings))
    # KIT-118 (P-09): a dot in a group name is `_`, as the builder writes it.
    friendlyNames <- gsub(".", "_", gsub("-", "_vs_", cleanStrings), fixed = TRUE)
    data.frame(contrasts = contrastStrings, friendly_names = friendlyNames,
        full_format = paste0(friendlyNames, "=", contrastStrings), stringsAsFactors = FALSE)
}

#' One table of an exported design (vroom, as today's imports), without vroom's live parser pointer
#' (the `problems` attribute): the session owner refuses to hold a process pointer (KIT-075), which
#' made today's imports fail with "payload contains a nonportable resource" (KIT-103 W2 finding).
#' The values and the column types are vroom's.
#' @noRd
readOmicsDesignTable <- function(path) {
    table <- vroom::vroom(path, show_col_types = FALSE)
    attr(table, "problems") <- NULL
    table
}

#' C10. Read an exported design directory, as today's lipid and metab imports: the configuration
#' (the directory's `config.ini`, else the source directory's, else empty), the design, the assays
#' the manifest names, the column mapping (`column_mapping.json`, else inferred with the omic's id
#' key and candidates) and the contrasts. Reads only; nothing is written.
#' @param importPath The exported design directory.
#' @param config The omic's GUI configuration.
#' @param sourceDir The source directory (its `config.ini` is the fallback).
#' @return A list: `config_list`, `design_matrix` (in file order), `assay_names`, `assay_list`,
#'   `column_mapping`, `contrasts_tbl`, `design_file`. A missing or empty input signals
#'   `omics_design_import_failed` with today's message.
#' @noRd
readOmicsDesignImport <- function(importPath, config, sourceDir = NULL) {
    files <- omicsDesignImportFiles(importPath)
    fail <- function(message) rlang::abort(message, class = "omics_design_import_failed")
    configList <- if (file.exists(files$config)) {
        logger::log_info("Loading config.ini from import directory.")
        readConfigFile(file = files$config)
    } else if (file.exists(file.path(sourceDir, "config.ini"))) {
        logger::log_info("Loaded config.ini from source_dir.")
        readConfigFile(file = file.path(sourceDir, "config.ini"))
    } else {
        logger::log_warn("No config.ini found. Using empty config.")
        list()
    }
    designMatrix <- readOmicsDesignTable(files$design)
    assayNames <- readLines(files$manifest)
    assayNames <- assayNames[nzchar(assayNames)]
    if (length(assayNames) == 0) fail("assay_manifest.txt is empty. No assays to load.")
    assayList <- list()
    for (assayName in assayNames) {
        assayFile <- file.path(importPath, paste0("data_cln_", assayName, ".tab"))
        if (!file.exists(assayFile)) fail(sprintf("Assay data file not found: %s", assayFile))
        assayList[[assayName]] <- readOmicsDesignTable(assayFile)
    }
    columnMapping <- if (file.exists(files$column_mapping)) {
        jsonlite::read_json(files$column_mapping, simplifyVector = TRUE)
    } else {
        logger::log_warn("column_mapping.json not found. Inferring from data structure.")
        inferOmicsDesignImportMapping(assayList, designMatrix, config)
    }
    list(config_list = configList, design_matrix = designMatrix, assay_names = assayNames, assay_list = assayList,
        column_mapping = columnMapping, contrasts_tbl = readOmicsDesignImportContrasts(files$contrasts),
        design_file = files$design)
}

#' C10. Commit an imported design, as today's imports: the omic's preflight when
#' `design_import$preflight` (lipid, and metab since P-11; with `design_import$preflight_contrasts`,
#' lipid, its normalised contrasts replace the file's), the bus
#' (configuration, mapping, design, contrasts, owner binding at import, the tables as `data_tbl`
#' and `data_cln`, the design's `tech_rep_group`), the object, the owner binding at design, the
#' provenance (with the imported design file when `design_import$record_design_file`, lipid; else
#' the owner's source directory, metab), the saved state and the QC baseline
#' (`design_import$qc_baseline`; a failure is a warning). The QC trigger, the tab status and the
#' notifications stay in the GUI.
#' @param imported What `readOmicsDesignImport()` returned.
#' @param sourceDir The source directory (`experiment_paths$source_dir`).
#' @return The object and the counts for the success text (`n_assays`, `n_samples`). A failed
#'   preflight is `omics_design_refused` ("Invalid <omic> design import: ..."), before the bus.
#' @noRd
commitOmicsDesignImport <- function(workflowData, imported, sourceDir, config) {
    design <- config$design
    importConfig <- design$design_import
    route <- config$shell$routes[[1L]]
    contrasts <- imported$contrasts_tbl
    if (isTRUE(importConfig$preflight)) {
        preflight <- design$validate(list(design_matrix = imported$design_matrix, data_cln = imported$assay_list,
            contrasts_tbl = contrasts, config_list = imported$config_list), imported$column_mapping)
        if (!isTRUE(preflight$valid)) {
            omicsDesignRefuse(sprintf("Invalid %s design import: %s", design$label_lower,
                paste(preflight$errors, collapse = "; ")))
        }
        normalised <- preflight$contrast_validation$contrasts_tbl
        if (isTRUE(importConfig$preflight_contrasts) && !is.null(normalised) && nrow(normalised) > 0L) contrasts <- normalised
    }
    # KIT-222 (owner's OD-2): a design file carrying technical-replicate marks is refused when a replicate number is
    # shared by an unmarked run or a key holds runs of two marking actions, before the bus.
    keyed <- as.data.frame(imported$design_matrix)
    if (all(c("group", "replicates") %in% names(keyed))) {
        keyed <- keyed[!is.na(keyed$group) & keyed$group != "", , drop = FALSE]
        keyed$tech_rep_group <- paste(keyed$group, keyed$replicates, sep = "_")
        refusal <- omicsDesignSaveRefusal(keyed, unmarkedWithoutColumn = FALSE)
        if (!is.null(refusal)) omicsDesignRefuse(sprintf("Invalid %s design import: %s", design$label_lower, refusal))
    }
    workflowData$config_list <- imported$config_list
    workflowData$column_mapping <- imported$column_mapping
    workflowData$design_matrix <- imported$design_matrix
    workflowData$contrasts_tbl <- contrasts
    # KIT-075: the owner is bound before the tables are held by reference.
    bindWorkflowOmicsOwnerAtImport(workflowData, route, workflowData$data_format)
    setWorkflowOmicsBusTables(workflowData, "data_tbl", imported$assay_list)
    setWorkflowOmicsBusTables(workflowData, "data_cln", imported$assay_list)
    workflowData$design_matrix <- workflowData$design_matrix |>
        dplyr::mutate(tech_rep_group = paste(group, replicates, sep = "_"))
    object <- design$build_object(list(data_cln = imported$assay_list, design_matrix = workflowData$design_matrix,
        config_list = workflowData$config_list), workflowData$column_mapping)
    # KIT-070: the session owner is bound to the omic's route here.
    bindWorkflowOwnerTypeAtDesign(workflowData, route, workflowData$data_format)
    manager <- workflowData$state_manager
    if (isTRUE(importConfig$record_design_file)) {
        design$record_provenance(manager, object, source_dir = sourceDir, design_file = imported$design_file)
    } else {
        design$record_provenance(manager, object, source_dir = manager$.context$experiment_paths$source_dir)
    }
    manager$saveState(state_name = design$state_name, s4_data_object = object,
        config_object = workflowData$config_list, description = importConfig$state_description)
    tryCatch({
        importConfig$qc_baseline(theObject = object, step_name = "1_Raw_Data", omics_type = design$label_lower,
            return_grid = FALSE, overwrite = TRUE)
        logger::log_info("Initialized QC progress tracking with raw data baseline")
    }, error = function(e) {
        logger::log_warn(paste("Could not initialize QC progress:", e$message))
    })
    logger::log_info(sprintf("Design import complete: '%s'", design$state_name))
    invisible(list(object = object, n_assays = length(imported$assay_list), n_samples = nrow(imported$design_matrix)))
}

# ------------------------------------------------------------------------------------------------
# KIT-247 (VAL-239, the owner's OD-1 (D)): a design saved again over later steps.
# ------------------------------------------------------------------------------------------------

#' KIT-247. A design table in a canonical form, so that two saves of one design compare equal: every column as text,
#' the columns in name order, the rows in run order, and `tech_rep_group` recomputed from `group` and `replicates` as
#' the save computes it (`buildProtDesignStateCheckpoint`, `commitOmicsDesignImport`). Every column is kept: RUV's
#' grouping, the correlation filter and the DA formula can use any of them (factors, batch).
#' @param design A design table, or NULL.
#' @return A data frame, or NULL.
#' @noRd
omicsDesignCanonicalTable <- function(design) {
    if (is.null(design)) return(NULL)
    design <- as.data.frame(design, stringsAsFactors = FALSE)
    if (all(c("group", "replicates") %in% names(design))) {
        design$tech_rep_group <- paste(design$group, design$replicates, sep = "_")
    }
    design[] <- lapply(design, omicsDesignCanonicalText)
    design <- design[, sort(names(design)), drop = FALSE]
    key <- if ("Run" %in% names(design)) "Run" else names(design)[[1L]]
    design <- design[order(design[[key]]), , drop = FALSE]
    rownames(design) <- NULL
    design
}

#' KIT-247. A column as text, an empty or blank value read as missing (a design file read back gives NA where the saved
#' design held ""). Kept from the VAL-288 commit at the withdrawal (VAL-331): without it an unchanged proteomics re-import
#' is judged changed (R1''' and R1'''b, F288b).
#' @noRd
omicsDesignCanonicalText <- function(column) {
    column <- as.character(column)
    column[!is.na(column) & !nzchar(trimws(column))] <- NA_character_
    column
}

#' KIT-247. Whether two of a design save's data inputs (a table, a list of tables, or saved references to them) hold
#' the same tables: the same saved entity, or the same content once read.
#' @noRd
omicsDesignSameData <- function(a, b, manager) {
    if (is.null(a) || is.null(b)) return(is.null(a) && is.null(b))
    if (isSavedWorkflowReference(a) && isSavedWorkflowReference(b) && identical(a$entity_id, b$entity_id)) return(TRUE)
    canonical <- function(value) {
        value <- resolveWorkflowOwnedValue(value, manager)
        tables <- if (is.data.frame(value)) list(value) else value
        lapply(tables, function(table) {
            if (!is.data.frame(table)) return(table)
            table <- as.data.frame(table, stringsAsFactors = FALSE)
            # As text: a table read back from an exported file and the saved one agree when their values do.
            table[] <- lapply(table, omicsDesignCanonicalText)
            rownames(table) <- NULL
            table
        })
    }
    identical(digest::digest(canonical(a), algo = "sha256"), digest::digest(canonical(b), algo = "sha256"))
}

#' KIT-247 (OD-1 (D); Re-check 1 R1-02, the coordinator's ruling): a pure comparison, made by the GUI frames before the
#' bus is hydrated, of a design save's result with the design the session holds. The state's side is read from the bus
#' as the last design commit left it (its design table and data input): the transcriptomics design state holds its
#' technical replicates summed into units (`trxDesignObject`), so its stored design is not comparable run by run.
#' @param manager The session owner.
#' @param results The save's result: `design_matrix` and `data` (the table, list or reference the commit reads).
#' @param stateName The omic's design state.
#' @param current The bus before hydration: `design_matrix` and `data`.
#' @return `"first"` (no design state in the active history: the first save, or no session owner), `"unchanged"` (the
#'   whole design table and the data equal: only the comparisons or the formula may differ) or `"changed"`.
#' @noRd
compareOmicsDesignResult <- function(manager, results, stateName, current) {
    if (!inherits(manager, "DIAStateAdapter") || !isTRUE(stateName %in% diaStateHistory(manager))) return("first")
    # A re-save is a design saved over an earlier saved design: with none on the bus there is nothing to compare
    # (correction pass after R2 45612382879; coordinator's ruling).
    if (is.null(current$design_matrix)) return("first")
    if (!identical(omicsDesignCanonicalTable(results$design_matrix), omicsDesignCanonicalTable(current$design_matrix))) {
        return("changed")
    }
    if (!isTRUE(omicsDesignSameData(results$data, current$data, manager))) return("changed")
    "unchanged"
}

#' KIT-247 (VAL-277, AUDIT-KIT-247-impl I02): the design and data the last SUCCESSFUL design commit left, recorded on
#' the session owner, so that a save that failed after replacing the bus cannot make the next save of the same design
#' look unchanged.
#' @noRd
rememberOmicsDesignCommit <- function(manager, design, data, config = NULL) {
    # KIT-277 (C03): the record also holds the saved formula, so that a later save of the same design is compared with it.
    if (inherits(manager, "DIAStateAdapter")) {
        manager$.context$kit247_design_compared <- list(design_matrix = design, data = data,
            formula = if (!is.null(config)) resolveOmicsDesignFormula(config))
    }
    invisible(NULL)
}

#' KIT-277 (C03). After a limited Save: the remembered record keeps its design and data and takes the saved formula.
#' @noRd
rememberOmicsDesignFormula <- function(manager, formulaString) {
    if (inherits(manager, "DIAStateAdapter")) {
        remembered <- manager$.context$kit247_design_compared
        if (!is.null(remembered)) {
            remembered$formula <- formulaString
            manager$.context$kit247_design_compared <- remembered
        }
    }
    invisible(NULL)
}

#' KIT-277. The model formula a design save would change: the saved configuration's formula against the remembered
#' record's (else the bus's), spaces removed.
#' @param manager,workflowData The session owner and the bus, before the save hydrates the bus.
#' @param config The save's configuration.
#' @return `c(from =, to =)`, or NULL when the formula is the same.
#' @noRd
omicsDesignFormulaChange <- function(manager, workflowData, config) {
    from <- omicsDesignCompared(manager, workflowData)$formula %||% resolveOmicsDesignFormula(workflowData$config_list)
    to <- resolveOmicsDesignFormula(config)
    if (identical(gsub("[[:space:]]", "", from), gsub("[[:space:]]", "", to))) return(NULL)
    c(from = from, to = to)
}

#' KIT-247 (VAL-277). What a design save is compared with: the last successful commit's record, else the bus (a session
#' opened from a file, before any commit in it).
#' @noRd
omicsDesignCompared <- function(manager, workflowData) {
    remembered <- if (inherits(manager, "DIAStateAdapter")) manager$.context$kit247_design_compared
    if (!is.null(remembered)) return(remembered)
    list(design_matrix = workflowData$design_matrix, data = workflowData$data_cln)
}

#' KIT-247. The states of the active history after the design state.
#' @noRd
omicsDesignLaterStates <- function(manager, stateName) {
    if (!inherits(manager, "DIAStateAdapter")) return(character())
    history <- diaStateHistory(manager)
    at <- match(stateName, history)
    if (is.na(at)) character() else history[-seq_len(at)]
}

#' KIT-247 (OD-1 (D)). After a changed design was committed: the active history is cut at the design state, as a
#' Quality Control Revert to the design would (`selectWorkflowStateReference`; no state record is deleted).
#' @return TRUE when states were set aside.
#' @noRd
cutOmicsDesignHistory <- function(manager, stateName) {
    if (!length(omicsDesignLaterStates(manager, stateName))) return(FALSE)
    selectWorkflowStateReference(manager, stateName)
    TRUE
}

#' KIT-247 (OD-4 (a), VAL-239): the omic's earlier integration input: its folder (`removeIntegrationInput`, as a
#' refused Summary Save removes it, VAL-148) and the top-level file the Summary's Save writes beside it
#' (`<omic>_<label>_final_s4.RDS`; transcriptomics `transcriptomics_normalized_data.tsv`).
#' @param integrationDir The project's integration directory.
#' @param omicType The omic.
#' @return The removed paths, invisibly.
#' @noRd
removeOmicIntegrationInput <- function(integrationDir, omicType) {
    if (is.null(integrationDir) || !dir.exists(integrationDir)) return(invisible(character()))
    folder <- file.path(integrationDir, omicType)
    removed <- if (dir.exists(folder)) folder else character()
    removeIntegrationInput(folder)
    files <- list.files(integrationDir, pattern = paste0("^", omicType, "_"), full.names = TRUE)
    files <- files[!dir.exists(files)]
    if (length(files) && !all(file.remove(files))) {
        abortIntegrationInput(sprintf("the earlier integration input %s could not be removed", paste(basename(files), collapse = ", ")))
    }
    invisible(c(removed, files))
}

#' KIT-247 (OD-1 (D), items 2 and 5). After the later states were set aside: Quality Control pending, every later tab
#' disabled, the omic's earlier integration input removed, and one notice (DRAFT wording, VAL-41).
#' @param workflowData The workflow bus.
#' @param omicType The omic.
#' @param experimentPaths The project paths (`integration_dir`, `base_dir`).
#' @param notify `shiny::showNotification`'s contract.
#' @noRd
completeOmicsDesignCut <- function(workflowData, omicType, experimentPaths, notify = shiny::showNotification,
                                   cutText = omicsDesignCutText) {
    status <- workflowData$tab_status
    names_status <- names(status)
    at <- match("quality_control", names_status)
    if (!is.na(at)) {
        status$quality_control <- "pending"
        for (name in names_status[seq_along(names_status) > at]) status[[name]] <- "disabled"
        workflowData$tab_status <- status
    }
    integrationDir <- experimentPaths$integration_dir %||%
        if (!is.null(experimentPaths$base_dir)) file.path(experimentPaths$base_dir, "integration")
    removed <- tryCatch(length(removeOmicIntegrationInput(integrationDir, omicType)) > 0L, error = function(e) {
        notify(cli::ansi_strip(conditionMessage(e)), type = "warning", duration = NULL)
        FALSE
    })
    notify(cutText(removed), type = "warning", duration = NULL)
    invisible(TRUE)
}

#' KIT-247. The texts (DRAFT wording, VAL-41; `docs/qa/workflow-execution/kit/kit-247/wording.md`).
#' @noRd
omicsDesignCutText <- function(removed = FALSE) {
    paste0("The design was saved again. The Quality Control, normalisation and later results made with the previous ",
        "design are set aside: apply Quality Control again, then run the normalisation and the later steps.",
        if (isTRUE(removed)) " The earlier integration input of this omic was removed." else "")
}

#' @noRd
omicsDesignResaveUnconfirmedText <- function() {
    paste0("The design was not saved: it differs from the saved design and later results (Quality Control, ",
        "normalisation and after) exist, and there is no application session to ask whether to set them aside. ",
        "Revert the Quality Control to the design state first, or save the design in the application, which asks.")
}

#' @noRd
omicsDesignResaveConfirmText <- function() {
    paste0("This design changes the samples, groups, technical replicates, factors or batch. Saving it sets aside the ",
        "Quality Control, normalisation and later results made with the previous design; you will apply them again.")
}

#' @noRd
omicsDesignUnchangedText <- function(contrastsChanged) {
    if (isTRUE(contrastsChanged)) {
        "Comparisons saved. Quality Control and normalisation are kept; run the differential analysis again for the new comparisons."
    } else {
        "Design saved; nothing changed."
    }
}

#' KIT-247. Whether a save's comparisons differ from the bus's.
#' @noRd
omicsDesignContrastsChanged <- function(new, old) {
    canonical <- function(x) if (is.null(x) || !nrow(as.data.frame(x))) NULL else {
        x <- as.data.frame(x, stringsAsFactors = FALSE)
        x[] <- lapply(x, as.character)
        rownames(x) <- NULL
        x
    }
    !identical(canonical(new), canonical(old))
}

#' KIT-247 (OD-1 (D)). A design save that changes nothing the steps after it read: the comparisons, the configuration
#' and (by the caller) the design files are updated; nothing is committed; the DA tab goes back to pending when the
#' comparisons changed (a prompt only: KIT-212's check compares tables, so the Summary still presents the earlier
#' analysis until the DA is run again); the pending design provenance, if any, is cleared.
#' @param daStatusKey The DA tab's status name (`differential_expression`, `differential_analysis`).
#' @noRd
applyOmicsDesignUnchanged <- function(workflowData, results, daStatusKey, notify = shiny::showNotification,
                                      quiet = FALSE) {
    changed <- omicsDesignContrastsChanged(results$contrasts_tbl, workflowData$contrasts_tbl)
    workflowData$contrasts_tbl <- results$contrasts_tbl
    workflowData$config_list <- results$config_list
    manager <- workflowData$state_manager
    if (inherits(manager, "DIAStateAdapter")) manager$.context$pending_provenance <- NULL
    status <- workflowData$tab_status
    if (changed && !is.null(status[[daStatusKey]]) && !identical(status[[daStatusKey]], "disabled")) {
        status[[daStatusKey]] <- "pending"
        workflowData$tab_status <- status
    }
    # KIT-277 (C04): a save that carries a limited change through says what it did, not that nothing changed.
    if (!isTRUE(quiet)) notify(omicsDesignUnchangedText(changed), type = "message")
    invisible(changed)
}

#' KIT-247 (OD-1 (D)). The GUI frames' gate before a design save hydrates the bus. With states after the design state:
#' `"unchanged"` runs `unchanged()`; `"changed"` asks first (Save runs `commit(cut = TRUE)`; Cancel changes nothing).
#' `"first"`, or nothing after the design state, runs `commit(cut = FALSE)` as before.
#' @param session The design tab's session (the confirmation's button is `confirm_design_resave` in its namespace).
#' @param verdict `compareOmicsDesignResult()`'s value.
#' @param later Whether the active history holds states after the design state.
#' @return The verdict, invisibly.
#' @noRd
runOmicsDesignResaveGate <- function(session, verdict, later, commit, unchanged,
                                     showModal = shiny::showModal, limited = character()) {
    # A first save, or a save with nothing after the design state, is committed as before (G1, G2).
    if (!isTRUE(later) || !verdict %in% c("unchanged", "changed")) {
        commit(cut = FALSE)
        return(invisible(verdict))
    }
    if (identical(verdict, "unchanged") && !length(limited)) {
        unchanged()
        return(invisible(verdict))
    }
    # KIT-277: a save whose design and data are unchanged but which changes the model formula, the organism or the FASTA
    # file after Accession Cleanup asks once; Save runs `unchanged()`, which carries the change through; Cancel changes
    # nothing.
    if (identical(verdict, "unchanged")) {
        if (is.null(session) || !is.function(session$ns)) {
            rlang::abort(omicsDesignLimitedUnconfirmedText(),
                class = c("omics_design_resave_unconfirmed", "omics_design_refused"))
        }
        registerOmicsDesignResaveConfirmation(session)
        assign(omicsDesignResaveKey(session), function() unchanged(), envir = session$userData)
        showModal(shiny::modalDialog(title = "Save the changed settings?",
            shiny::p(omicsDesignLimitedText(limited)),
            footer = shiny::tagList(shiny::modalButton("Cancel"),
                shiny::actionButton(session$ns("confirm_design_resave"), "Save", class = "btn-warning")),
            easyClose = FALSE))
        return(invisible("limited"))
    }
    # No session to ask (a scripted call): a changed design over later results cannot be confirmed, so it is refused and
    # nothing is committed (coordinator's ruling; ticket r1.4 section 19).
    if (is.null(session) || !is.function(session$ns)) {
        rlang::abort(omicsDesignResaveUnconfirmedText(), class = c("omics_design_resave_unconfirmed", "omics_design_refused"))
    }
    registerOmicsDesignResaveConfirmation(session)
    assign(omicsDesignResaveKey(session), function() commit(cut = TRUE), envir = session$userData)
    showModal(shiny::modalDialog(title = "Save the changed design?",
        shiny::p(omicsDesignResaveConfirmText()),
        footer = shiny::tagList(shiny::modalButton("Cancel"),
            shiny::actionButton(session$ns("confirm_design_resave"), "Save", class = "btn-warning")),
        easyClose = FALSE))
    invisible(verdict)
}

#' KIT-247. The pending save's key in the session's user data (the design tab's namespace: one per omic).
#' @noRd
omicsDesignResaveKey <- function(session) paste0("kit247_design_resave_", session$ns(""))

#' KIT-247. The confirmation's observer, registered once per design tab session.
#' @noRd
registerOmicsDesignResaveConfirmation <- function(session) {
    key <- omicsDesignResaveKey(session)
    flag <- paste0(key, "_observer")
    if (isTRUE(get0(flag, envir = session$userData, inherits = FALSE))) return(invisible(NULL))
    assign(flag, TRUE, envir = session$userData)
    shiny::observeEvent(session$input$confirm_design_resave, {
        pending <- get0(key, envir = session$userData, inherits = FALSE)
        assign(key, NULL, envir = session$userData)
        shiny::removeModal()
        if (is.function(pending)) pending()
    }, ignoreInit = TRUE)
    invisible(NULL)
}

# ------------------------------------------------------------------------------------------------
# KIT-277 (VAL-331, VAL-345, VAL-353): a design saved again with only the model formula, the organism or the FASTA file
# changed asks first, and Save carries the change through to every consumer. All texts are DRAFTS (VAL-41).
# ------------------------------------------------------------------------------------------------

#' KIT-277. The sentence of each kind of limited change, in the order formula, organism, FASTA.
#' @param limited A named list: `formula = c(from, to)`, `organism = c(from, to)` (labels), `fasta = list(from, to, steps)`.
#' @noRd
omicsDesignLimitedText <- function(limited) {
    kinds <- intersect(c("formula", "organism", "fasta"), names(limited))
    what <- vapply(kinds, function(kind) switch(kind,
        formula = sprintf("the model formula (from `%s` to `%s`)", limited$formula[["from"]], limited$formula[["to"]]),
        organism = sprintf("the organism (from %s to %s)", limited$organism[["from"]], limited$organism[["to"]]),
        fasta = sprintf("the FASTA file (from `%s` to `%s`)", limited$fasta$from, limited$fasta$to)), "")
    effects <- vapply(kinds, function(kind) switch(kind,
        formula = paste("Saving the formula keeps Quality Control and normalisation, replaces the formula on the",
            "Differential analysis tab (a formula typed there is replaced) and sets the differential analysis back:",
            "Copy and the Report wait until you run it again."),
        organism = paste("Saving the organism keeps Quality Control and normalisation, sets the differential analysis",
            "and enrichment back and removes this omic's earlier integration input. The protein annotation for the new",
            "organism is fetched when you next run the differential analysis or enrichment."),
        fasta = paste0("Accession Cleanup used the earlier file, so Accession Cleanup and every later step are set ",
            "aside: ", paste(limited$fasta$steps, collapse = ", "), ". Saving sets them aside and removes this omic's ",
            "earlier integration input; you apply them again.")), "")
    paste("This save changes", paste(what, collapse = ", "), "only.", paste(effects, collapse = " "))
}

#' @noRd
omicsDesignLimitedUnconfirmedText <- function() {
    paste0("The settings were not saved: the model formula, the organism or the FASTA file differs from the saved ",
        "design and later results exist, and there is no application session to ask whether to set them aside. ",
        "Save the design in the application, which asks.")
}

#' KIT-277. The notice of a limited Save, one sentence per kind.
#' @noRd
omicsDesignLimitedNotice <- function(kind, omicType = NULL, detail = NULL) {
    switch(kind,
        formula = if (identical(omicType, "transcriptomics")) {
            paste("Formula saved. The transcriptomics differential analysis models the groups and the RUV factors from",
                "the design, so this formula does not change it.")
        } else {
            "Formula saved. Quality Control and normalisation are kept; run the differential analysis again."
        },
        organism = "Organism saved. Run the differential analysis, then enrichment, again.",
        fasta_used = sprintf(paste("FASTA file saved (from `%s` to `%s`). Accession Cleanup will use it."),
            detail$from, detail$to),
        fasta_cut = sprintf(paste("FASTA file saved (from `%s` to `%s`). The history went back to the step before",
            "Accession Cleanup; apply Accession Cleanup and the later steps again."), detail$from, detail$to))
}

#' KIT-277 (4.5). One set-back for every consumer: the head's DA reference and enrichment marker are cleared and the head
#' is published once (the earlier entities stay in the store); the DA and enrichment tabs go to pending and the state
#' trigger runs the tabs' observers. `design_formula` records, for lipidomics and metabolomics, the formula the DA tab's
#' box takes until the next DA run (R01).
#' @param workflowData The bus.
#' @param daStatusKey,enrichmentStatusKey The tabs' status names (`enrichment` NULL when the omic has none).
#' @param designFormula NULL, or the formula to record on the head.
#' @noRd
setOmicsDaBack <- function(workflowData, daStatusKey, enrichmentStatusKey = NULL, designFormula = NULL) {
    manager <- workflowData$state_manager
    if (inherits(manager, "DIAStateAdapter")) {
        withWorkflowStateScope(manager, function() {
            snapshot <- manager$.context$snapshot
            snapshot$da_results <- NULL
            snapshot$enrichment_results <- NULL
            if (!is.null(designFormula)) snapshot$design_formula <- designFormula
            publishWorkflowStateHead(manager$.context, snapshot)
        })
    }
    status <- workflowData$tab_status
    for (key in c(daStatusKey, enrichmentStatusKey)) {
        if (!is.null(status[[key]]) && !identical(status[[key]], "disabled")) status[[key]] <- "pending"
    }
    workflowData$tab_status <- status
    workflowData$state_update_trigger <- Sys.time()
    invisible(NULL)
}

#' KIT-277 (4.1). On Save, a formula-only change carried through: proteomics records the formula on the head's settings
#' (the DA tab reads them whenever it is opened); lipidomics and metabolomics send the DA tab's box a bus event and keep
#' the formula on the head until the next DA run; transcriptomics does nothing the DA reads (a notice says so). The DA is
#' set back for every consumer except transcriptomics. The remembered record takes the formula (C03).
#' @return The notice's text.
#' @noRd
carryOmicsDesignFormula <- function(workflowData, omicType, formulaString, daStatusKey) {
    manager <- workflowData$state_manager
    if (!identical(omicType, "transcriptomics")) {
        if (identical(omicType, "proteomics")) {
            commitDIADesignFormula(manager, formulaString)
            setOmicsDaBack(workflowData, daStatusKey, "enrichment_analysis")
        } else {
            setOmicsDaBack(workflowData, daStatusKey, designFormula = formulaString)
            saved <- workflowData$design_formula_saved
            workflowData$design_formula_saved <- list(formula = formulaString, counter = (saved$counter %||% 0L) + 1L)
        }
    }
    rememberOmicsDesignFormula(manager, formulaString)
    omicsDesignLimitedNotice("formula", omicType)
}

#' KIT-277 (R01). The formula a design Save placed for a DA tab whose box is its own (lipidomics, metabolomics): the head's
#' `design_formula`, only while the store holds no DA (a DA run ends it; after it the box is the user's).
#' @noRd
omicDaDesignFormula <- function(workflowData) {
    manager <- workflowData$state_manager
    if (!inherits(manager, "DIAStateAdapter")) return(NULL)
    snapshot <- manager$.context$snapshot
    if (!is.null(snapshot$da_results)) return(NULL)
    snapshot$design_formula
}

#' KIT-277 (R01). One observer per DA tab: it writes the box with the formula once per Save (the bus event's counter), and
#' when the tab's server starts (a project reopened before the DA was run); never as a standing rule, so a formula typed
#' on the tab is not overwritten by a later head change.
#' @param session The DA tab's session.
#' @param workflowData The bus.
#' @noRd
registerOmicDaDesignFormula <- function(session, workflowData) {
    applied <- shiny::reactiveVal(NULL)
    shiny::observe({
        counter <- workflowData$design_formula_saved$counter
        formula <- shiny::isolate(omicDaDesignFormula(workflowData))
        key <- paste(counter %||% 0L, formula)
        if (is.null(formula) || identical(shiny::isolate(applied()), key)) return(invisible(NULL))
        applied(key)
        shiny::updateTextAreaInput(session, "formula_string", value = formula)
    })
    invisible(NULL)
}

#' KIT-277. The shared omic design's `unchanged()`: the comparisons and configuration written as KIT-247 does, and, when the
#' formula changed, the formula carried through (a failure is shown as a design save error) with its own notice.
#' @param results `contrasts_tbl` and `config_list`.
#' @param config The omic's GUI configuration.
#' @param formula `omicsDesignFormulaChange()`'s value.
#' @noRd
runOmicsDesignUnchangedSave <- function(workflowData, results, config, formula = NULL,
                                        notify = shiny::showNotification) {
    daKey <- config$da$status_key %||% "differential_analysis"
    applyOmicsDesignUnchanged(workflowData, results, daKey, notify, quiet = !is.null(formula))
    if (!is.null(formula)) {
        tryCatch(notify(carryOmicsDesignFormula(workflowData, config$shell$omic_type, formula[["to"]], daKey),
            type = "message"), error = function(e) {
                msg <- paste("Error saving design matrix results:", conditionMessage(e))
                logger::log_error(msg)
                notify(msg, type = "error", duration = 15)
            })
    }
    invisible(NULL)
}
