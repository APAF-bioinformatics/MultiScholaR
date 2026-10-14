# Handler 6: Results Table Rendering
da_server_table_render_handler <- function(input, output, session, da_data) {
  # KIT-024: on the native session the results table is paged server-side from the store (the significance
  # selection applied in the read); the widget shows the selection's row count.
  pager <- new.env(parent = emptyenv())
  output$da_results_table <- DT::renderDT({
    shiny::req(input$table_contrast, da_data$da_results_list)
    pager$read <- NULL

    tryCatch(
      {
        if (isDIAResultManifest(da_data$da_results_list)) {
          manifest <- da_data$da_results_list
          contrast <- input$table_contrast
          significance <- input$table_significance %||% "all"
          qvalue <- input$da_q_val_thresh
          lfc <- input$treat_lfc_cutoff
          page_length <- as.integer(input$table_max_rows %||% 25L)
          # The owner is resolved here, in the render; the page requests run outside any reactive domain.
          manager <- getDIAActiveManager()
          first <- readDIAResultPage(manifest, contrast, significance, qvalue, lfc, limit = max(1L, min(page_length, 10000L)),
            manager = manager)
          pager$rows <- as.integer(if (identical(significance, "all")) first$full_population_rows else first$counts[[significance]])
          pager$read <- function(offset, limit) readDIAResultPage(manifest, contrast, significance, qvalue, lfc,
            offset = offset, limit = limit, manager = manager)$data
          return(DT::datatable(first$data, options = diaDisplayTableOptions(page_length, scrollX = TRUE,
            dom = "Bfrtip", buttons = c("copy", "csv", "excel")), extensions = "Buttons", rownames = FALSE) |>
            DT::formatRound(columns = c("log2FC", "raw_pvalue", "fdr_qvalue"), digits = 4))
        }
        # Get DA results from the new format
        if (!is.null(da_data$da_results_list$da_proteins_long)) {
          # Debug: Show available contrasts in data vs selected contrast
          available_contrasts <- unique(da_data$da_results_list$da_proteins_long$comparison)
          cat(sprintf("   DA TABLE: Available contrasts in data: %s\n", paste(available_contrasts, collapse = ", ")))
          cat(sprintf("   DA TABLE: Selected contrast from UI: %s\n", input$table_contrast))

          # Filter for selected contrast (input$table_contrast now matches comparison column)
          current_results <- da_data$da_results_list$da_proteins_long |>
            dplyr::filter(comparison == input$table_contrast)

          cat(sprintf("   DA TABLE: Found %d rows for contrast %s\n", nrow(current_results), input$table_contrast))

          if (nrow(current_results) > 0) {
            # Filter based on significance selection
            if (input$table_significance == "significant") {
              current_results <- current_results[current_results$fdr_qvalue < input$da_q_val_thresh, ]
            } else if (input$table_significance == "up") {
              current_results <- current_results[current_results$fdr_qvalue < input$da_q_val_thresh &
                current_results$log2FC > input$treat_lfc_cutoff, ]
            } else if (input$table_significance == "down") {
              current_results <- current_results[current_results$fdr_qvalue < input$da_q_val_thresh &
                current_results$log2FC < -input$treat_lfc_cutoff, ]
            }

            # Limit rows
            if (nrow(current_results) > input$table_max_rows) {
              current_results <- current_results[1:input$table_max_rows, ]
            }

            # Select relevant columns for display - use the correct protein ID column
            if (!is.null(da_data$da_results_list$theObject)) {
              protein_id_column <- da_data$da_results_list$theObject@protein_id_column
              cat(sprintf("   DA TABLE: Using protein ID column = %s\n", protein_id_column))
            } else {
              # Fallback: try to find protein ID column in the data
              possible_protein_cols <- c("Protein.Ids", "uniprot_acc", "protein_id", "Protein_ID")
              protein_id_column <- intersect(possible_protein_cols, names(current_results))[1]
              cat(sprintf("   DA TABLE: Using fallback protein ID column = %s\n", protein_id_column))
            }

            display_columns <- c(protein_id_column, "log2FC", "raw_pvalue", "fdr_qvalue")
            current_results <- current_results |>
              dplyr::select(any_of(display_columns))

            DT::datatable(
              current_results,
              options = list(
                pageLength = 25,
                scrollX = TRUE,
                dom = "Bfrtip",
                buttons = c("copy", "csv", "excel")
              ),
              extensions = "Buttons"
            ) |>
              DT::formatRound(columns = c("log2FC", "raw_pvalue", "fdr_qvalue"), digits = 4)
          } else {
            # No results for this contrast
            DT::datatable(data.frame(Message = "No results available for selected contrast"))
          }
        } else {
          # No DA results available
          DT::datatable(data.frame(Message = "No DE analysis results available"))
        }
      },
      error = function(e) {
        cat(paste("*** ERROR in DA results table:", e$message, "\n"))
        DT::datatable(data.frame(Message = paste("Error:", e$message)))
      }
    )
  }, server = TRUE, funcFilter = diaDisplayTableFilter(pager))

  # Summary statistics
  output$da_summary_stats <- shiny::renderText({
    shiny::req(input$table_contrast, da_data$da_results_list)

    tryCatch(
      {
        if (isDIAResultManifest(da_data$da_results_list)) {
          counts <- readDIAResultPage(da_data$da_results_list, input$table_contrast,
            qvalue = input$da_q_val_thresh, lfc = input$treat_lfc_cutoff, limit = 1L)$counts
          return(paste(sprintf("Total genes: %d", counts$total),
            sprintf("Significant (q < %.3f): %d", input$da_q_val_thresh, counts$significant),
            sprintf("Up-regulated: %d", counts$up), sprintf("Down-regulated: %d", counts$down),
            sprintf("Fold-change cutoff: %.2f", input$treat_lfc_cutoff), sep = "\n"))
        }
        if (!is.null(da_data$da_results_list$da_proteins_long)) {
          # Filter for selected contrast
          current_results <- da_data$da_results_list$da_proteins_long |>
            dplyr::filter(comparison == input$table_contrast)

          if (nrow(current_results) > 0) {
            total_genes <- nrow(current_results)
            significant <- sum(current_results$fdr_qvalue < input$da_q_val_thresh, na.rm = TRUE)
            up_reg <- sum(current_results$fdr_qvalue < input$da_q_val_thresh &
              current_results$log2FC > input$treat_lfc_cutoff, na.rm = TRUE)
            down_reg <- sum(current_results$fdr_qvalue < input$da_q_val_thresh &
              current_results$log2FC < -input$treat_lfc_cutoff, na.rm = TRUE)

            paste(
              sprintf("Total genes: %d", total_genes),
              sprintf("Significant (q < %.3f): %d", input$da_q_val_thresh, significant),
              sprintf("Up-regulated: %d", up_reg),
              sprintf("Down-regulated: %d", down_reg),
              sprintf("Fold-change cutoff: %.2f", input$treat_lfc_cutoff),
              sep = "\n"
            )
          } else {
            "No results available for selected contrast"
          }
        } else {
          "No DE analysis results available"
        }
      },
      error = function(e) {
        paste("Error calculating statistics:", e$message)
      }
    )
  })

  # Download handler for results
  output$download_da_results <- shiny::downloadHandler(
    filename = function() {
      paste0("DA_results_", Sys.Date(), ".zip")
    },
    content = function(file) {
      # Placeholder for creating downloadable zip file
      writeLines("Differential expression results would be packaged here", file)
    }
  )
}

