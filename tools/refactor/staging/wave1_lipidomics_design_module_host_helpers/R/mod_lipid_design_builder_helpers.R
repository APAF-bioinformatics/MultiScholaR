registerLipidDesignBuilderModule <- function(
    workflowData,
    moduleId = "builder",
    builderServerExists = exists("mod_lipid_design_builder_server"),
    builderServerFn = NULL,
    reactiveFn = shiny::reactive,
    reactiveValFn = shiny::reactiveVal
) {
    if (isTRUE(builderServerExists)) {
        if (is.null(builderServerFn)) {
            builderServerFn <- mod_lipid_design_builder_server
        }

        return(builderServerFn(
            moduleId,
            data_tbl = reactiveFn(workflowData$data_tbl),
            config_list = reactiveFn(workflowData$config_list),
            column_mapping = reactiveFn(workflowData$column_mapping),
            existing_design_matrix = reactiveFn(workflowData$design_matrix),
            existing_contrasts = reactiveFn(workflowData$contrasts_tbl)
        ))
    }

    reactiveValFn(NULL)
}

runLipidDesignBuilderObserverShell <- function(
    results,
    workflowData,
    experimentPaths,
    qcTrigger = NULL
) {
    shiny::showModal(shiny::modalDialog(
        title = "Processing Design Matrix"
        , shiny::div(
            style = "text-align: center; padding: 20px;"
            , shiny::icon("spinner", class = "fa-spin fa-3x")
            , shiny::br()
            , shiny::br()
            , shiny::p("Saving design matrix and preparing data...")
        )
        , footer = NULL
        , easyClose = FALSE
    ))

    logger::log_info("Received results from lipidomics design builder. Saving to workflow and disk.")

    workflowData$design_matrix <- results$design_matrix
    workflowData$data_cln <- results$data_cln
    workflowData$contrasts_tbl <- results$contrasts_tbl
    workflowData$config_list <- results$config_list

    if (!is.null(results$contrasts_tbl)) {
        assign("contrasts_tbl", results$contrasts_tbl, envir = .GlobalEnv)
        logger::log_info("Updated contrasts_tbl in global environment.")
    }

    assign("config_list", workflowData$config_list, envir = .GlobalEnv)
    logger::log_info("Updated global config_list.")

    source_dir <- experimentPaths$source_dir
    if (is.null(source_dir) || !dir.exists(source_dir)) {
        msg <- "Could not find source directory to save files."
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeModal()
        return(invisible(NULL))
    }

    tryCatch({
        design_matrix_path <- file.path(source_dir, "design_matrix.tab")
        logger::log_info(paste("Writing design matrix to:", design_matrix_path))
        utils::write.table(results$design_matrix, file = design_matrix_path
            , sep = "\t", row.names = FALSE, quote = FALSE)

        if (!is.null(results$contrasts_tbl) && nrow(results$contrasts_tbl) > 0) {
            contrast_path <- file.path(source_dir, "contrast_strings.tab")
            logger::log_info(paste("Writing contrasts to:", contrast_path))
            writeLines(results$contrasts_tbl$contrasts, contrast_path)
        }

        assay_names <- names(results$data_cln)
        for (assay_name in assay_names) {
            assay_path <- file.path(source_dir, paste0("data_cln_", assay_name, ".tab"))
            logger::log_info(paste("Writing assay data to:", assay_path))
            utils::write.table(results$data_cln[[assay_name]], file = assay_path
                , sep = "\t", row.names = FALSE, quote = FALSE)
        }

        manifest_path <- file.path(source_dir, "assay_manifest.txt")
        writeLines(assay_names, manifest_path)
        logger::log_info(sprintf("Saved assay manifest with %d assays: %s",
            length(assay_names), paste(assay_names, collapse = ", ")))

        col_map <- workflowData$column_mapping
        if (!is.null(col_map)) {
            col_map_path <- file.path(source_dir, "column_mapping.json")
            jsonlite::write_json(col_map, col_map_path, auto_unbox = TRUE)
            logger::log_info("Saved column_mapping.json")
        }

        manifest_json_path <- file.path(source_dir, "manifest.json")
        manifest_data <- list(
            data_path = "assay_manifest.txt",
            design_matrix_path = "design_matrix.tab",
            contrast_strings_path = "contrast_strings.tab",
            column_mapping_path = "column_mapping.json"
        )
        jsonlite::write_json(manifest_data, manifest_json_path, auto_unbox = TRUE, pretty = TRUE)
        logger::log_info(paste("Saved manifest.json to:", manifest_json_path))

        if (!is.null(workflowData$config_list)) {
            config_path <- file.path(source_dir, "config.ini")
            logger::log_info(paste("Writing config.ini to:", config_path))
            tryCatch({
                ini::write.ini(workflowData$config_list, config_path)
                logger::log_info("Saved config.ini.")
            }, error = function(e) {
                logger::log_warn(paste("Could not save config.ini:", e$message))
            })
        }

        s4_obj <- createLipidomicsAssayData(
            lipid_data = results$data_cln
            , design_matrix = results$design_matrix
            , lipid_id_column = col_map$lipid_id_col
            , annotation_id_column = if (!is.null(col_map$annotation_col) && !is.na(col_map$annotation_col) && nzchar(col_map$annotation_col)) {
                col_map$annotation_col
            } else {
                NA_character_
            }
            , sample_id = "Run"
            , group_id = "group"
            , technical_replicate_id = "tech_rep_group"
            , database_identifier_type = "Unknown"
            , internal_standard_regex = if (!is.null(col_map$is_pattern) && !is.na(col_map$is_pattern)) {
                col_map$is_pattern
            } else {
                NA_character_
            }
            , args = results$config_list
        )

        if (is.null(workflowData$state_manager)) {
            workflowData$state_manager <- WorkflowState$new("lipidomics")
        }

        logger::log_info("Saving LipidomicsAssayData S4 object to state manager as 'lipid_raw_data_s4'")
        workflowData$state_manager$saveState(
            state_name = "lipid_raw_data_s4"
            , s4_data_object = s4_obj
            , config_object = results$config_list
            , description = "Initial LipidomicsAssayData S4 object created after design matrix"
        )

        if (!is.null(qcTrigger)) {
            qcTrigger(TRUE)
            logger::log_info("QC trigger set to TRUE")
        }

        updated_status <- workflowData$tab_status
        updated_status$design_matrix <- "complete"
        workflowData$tab_status <- updated_status

        shiny::removeModal()
        shiny::showNotification("Design matrix and contrasts saved successfully!", type = "message")

        logger::log_info(sprintf("Design save complete. Files saved to: %s", source_dir))
        logger::log_info(sprintf("Saved: design_matrix.tab, %d assay files, assay_manifest.txt, column_mapping.json, config.ini",
            length(assay_names)))
    }, error = function(e) {
        msg <- paste("Error saving design matrix results:", e$message)
        logger::log_error(msg)
        shiny::removeModal()
        shiny::showNotification(msg, type = "error", duration = 15)
    })
}

registerLipidDesignBuilderResultsObserver <- function(
    builderResultsRv,
    workflowData,
    experimentPaths,
    qcTrigger = NULL,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    runBuilderObserverShell = runLipidDesignBuilderObserverShell
) {
    observeEventFn(builderResultsRv(), {
        results <- builderResultsRv()
        reqFn(results)

        runBuilderObserverShell(
            results = results,
            workflowData = workflowData,
            experimentPaths = experimentPaths,
            qcTrigger = qcTrigger
        )
    }, ignoreNULL = TRUE)
}

