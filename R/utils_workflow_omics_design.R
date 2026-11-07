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
    friendlyNames <- gsub("-", "_vs_", cleanStrings)
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
#' `design_import$preflight` (lipid; its normalised contrasts replace the file's), the bus
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
        if (!is.null(normalised) && nrow(normalised) > 0L) contrasts <- normalised
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
