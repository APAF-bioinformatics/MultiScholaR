# Handler 5: Heatmap Rendering
da_server_heatmap_render_handler <- function(input, output, session, ns, da_data, experiment_paths) {
  # --- TESTTHAT CHECKPOINT CP09 (see test-prot-09-heatmap.R) ---
  shiny::observeEvent(list(input$heatmap_contrast, da_data$da_results_list), {
  shiny::req(input$heatmap_contrast, da_data$da_results_list)

  plot_data_structure <- list(
    da_proteins_long = da_data$da_results_list$da_proteins_long,
    theObject = da_data$da_results_list$theObject
  )

  .capture_checkpoint(list(
    da_results_list = plot_data_structure,
    selected_contrast = input$heatmap_contrast,
    top_n_genes = input$heatmap_top_n,
    clustering_method = input$heatmap_cluster_method,
    distance_method = input$heatmap_distance_method,
    cluster_rows = input$heatmap_clustering %in% c("both", "row"),
    cluster_cols = input$heatmap_clustering %in% c("both", "column"),
    scale_data = input$heatmap_scaling,
    color_scheme = input$heatmap_color_scheme,
    show_gene_names = input$heatmap_show_labels,
    da_q_val_thresh = input$da_q_val_thresh,
    tree_cut_method = input$heatmap_tree_cut_method,
    n_clusters = input$heatmap_n_clusters,
    cut_height = input$heatmap_cut_height,
    min_cluster_size = input$heatmap_min_cluster_size
  ), "cp09", "heatmap_input")
  }, ignoreInit = TRUE)
  # --- END CP09 ---
  output$heatmap_plot <- shiny::renderPlot({
    shiny::req(input$heatmap_contrast, da_data$da_results_list)

    tryCatch(
      {
        # Generate heatmap using new function
        # Create a compatible data structure for the plotting function
        plot_data_structure <- list(
          da_proteins_long = da_data$da_results_list$da_proteins_long,
          theObject = da_data$da_results_list$theObject
        )

        # The heatmap contrast input should now match the comparison column in da_proteins_long
        cat(sprintf("   HEATMAP: Looking for contrast = %s\n", input$heatmap_contrast))
        
        heatmap_result <- generateProtDAHeatmap(
          da_results_list = plot_data_structure,
          selected_contrast = input$heatmap_contrast,
          top_n_genes = input$heatmap_top_n,
          clustering_method = input$heatmap_cluster_method,
          distance_method = input$heatmap_distance_method,
          cluster_rows = input$heatmap_clustering %in% c("both", "row"),
          cluster_cols = input$heatmap_clustering %in% c("both", "column"),
          scale_data = input$heatmap_scaling,
          color_scheme = input$heatmap_color_scheme,
          show_gene_names = input$heatmap_show_labels,
          da_q_val_thresh = input$da_q_val_thresh,
          tree_cut_method = input$heatmap_tree_cut_method,
          n_clusters = input$heatmap_n_clusters,
          cut_height = input$heatmap_cut_height,
          min_cluster_size = input$heatmap_min_cluster_size
        )

        if (!is.null(heatmap_result)) {
          # Check if it's a list (new format) or just the plot (old format/fallback)
          if (is.list(heatmap_result) && "plot" %in% names(heatmap_result)) {
            # Store cluster info and plot object
            da_data$current_row_clusters <- heatmap_result$row_clusters
            da_data$current_col_clusters <- heatmap_result$col_clusters
            da_data$current_heatmap_plot <- heatmap_result$plot

            # Return plot for rendering
            heatmap_result$plot
          } else {
            # Fallback for old return format
            heatmap_result
          }
        } else {
          # No significant genes found
          plot(1, 1,
            type = "n", axes = FALSE, xlab = "", ylab = "",
            main = paste("No significant genes found for contrast:", input$heatmap_contrast)
          )
          text(1, 1, "Adjust significance thresholds\nor select different contrast", cex = 1.2)
        }
      },
      error = function(e) {
        cat(paste("*** ERROR in heatmap generation:", e$message, "\n"))
        plot(1, 1,
          type = "n", axes = FALSE, xlab = "", ylab = "",
          main = "Error generating heatmap"
        )
        text(1, 1, paste("Error:", e$message), cex = 1.0)
      }
    )
  })

  # Cluster Summary Output
  output$cluster_summary <- shiny::renderPrint({
    shiny::req(input$heatmap_tree_cut_method != "none")

    if (is.null(da_data$current_row_clusters)) {
      cat("No clusters defined. Enable clustering and tree cutting on the heatmap.")
      return()
    }

    clusters <- da_data$current_row_clusters
    n_clusters <- length(unique(clusters))

    cat(sprintf("Total Clusters: %d\n", n_clusters))
    cat("----------------------------------\n")

    unique_sorted_clusters <- sort(unique(clusters))
    purrr::walk(unique_sorted_clusters, function(i) {
      members <- names(clusters)[clusters == i]
      cat(sprintf("\nCluster %d (%d proteins):\n", i, length(members)))

      # If we have gene names, try to show them too
      # This requires looking up the proteins in the da_results list
      cat(paste(head(members, 20), collapse = ", "))
      if (length(members) > 20) {
        cat(sprintf(", ... and %d more", length(members) - 20))
      }
      cat("\n")
    })
  })

  # Save Heatmap Observer
  shiny::observeEvent(input$save_heatmap, {
    shiny::req(da_data$current_heatmap_plot, experiment_paths$publication_graphs_dir)

    logger::log_info("Save Heatmap button clicked")

    # Collect parameters
    params <- list(
      contrast = input$heatmap_contrast,
      top_n = input$heatmap_top_n,
      clustering_method = input$heatmap_cluster_method,
      distance_method = input$heatmap_distance_method,
      cluster_rows = input$heatmap_clustering %in% c("both", "row"),
      cluster_cols = input$heatmap_clustering %in% c("both", "column"),
      scaling = input$heatmap_scaling,
      color_scheme = input$heatmap_color_scheme,
      tree_cut_method = input$heatmap_tree_cut_method,
      n_clusters = input$heatmap_n_clusters,
      cut_height = input$heatmap_cut_height,
      min_cluster_size = input$heatmap_min_cluster_size
    )

    # File prefix
    prefix <- paste0("prot_", input$heatmap_contrast)
    # Sanitize prefix
    prefix <- gsub("[^A-Za-z0-9_]", "_", prefix)

    # Call helper function
    save_heatmap_products(
      heatmap_obj = da_data$current_heatmap_plot,
      row_clusters = da_data$current_row_clusters,
      params_list = params,
      output_dir = experiment_paths$publication_graphs_dir,
      file_prefix = prefix
    )

    shiny::showNotification(
      "Heatmap and cluster info saved to publication_graphs/Heatmap",
      type = "message",
      duration = 5
    )
  })
}

