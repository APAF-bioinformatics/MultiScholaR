# ----------------------------------------------------------------------------
# updateProteinFiltering
# ----------------------------------------------------------------------------
#' @title Update and Visualize Filtering Progress
#' @description Tracks and visualizes the impact of filtering steps on peptide
#'   and protein counts. Updates a scoped `FilteringProgress` object and optionally
#'   saves plots summarizing the changes. Handles both peptide-level and
#'   protein-level data inputs.
#'
#' @details
#' This function acts as a central hub for monitoring data reduction throughout
#' a filtering workflow. It performs the following actions:
#' \itemize{
#'   \item Initializes or retrieves an S4 object named `filtering_progress`
#'     of class `FilteringProgress`.
#'   \item Calculates key metrics (total unique proteins, proteins per run,
#'     total unique peptides, peptides per protein distribution, peptides per run)
#'     based on the input `data`. Peptide metrics are only calculated or updated
#'     if `data` is identified as peptide-level data. For protein-level data,
#'     peptide metrics from the last peptide step (if any) are carried forward or
#'     initialized as empty/NA.
#'   \item Adds or updates these metrics in the `filtering_progress` object
#'     under the specified `step_name`.
#'   \item Generates summary plots using `ggplot2`:
#'     \itemize{
#'       \item Bar plot of total unique proteins per step.
#'       \item Bar plot of total unique peptides per step (or placeholder if only protein data).
#'       \item Box plot of peptides per protein distribution per step (or placeholder).
#'       \item Line plot of proteins per run across steps.
#'       \item Line plot of peptides per run across steps (or placeholder).
#'     }
#'   \item If `omic_type` and `experiment_label` are provided and valid paths can be
#'     derived from the global `project_dirs` object, the generated plots are saved
#'     as PNG files into the derived `time_dir`. Warnings are issued if paths cannot be
#'     derived or if `project_dirs` is not found.
#'   \item If `return_grid` is `TRUE`, arranges the plots into a single grid using
#'     `gridExtra` and returns the grid object (grob). Also saves this combined grid
#'     if plot saving is enabled.
#'   \item If `return_grid` is `FALSE` (default), prints each plot individually
#'     and returns the list of plot objects invisibly.
#' }
#'
#' **Important:** By default this function modifies a global variable named
#' `filtering_progress`. Callers can provide `progress_env` for session ownership.
#' For saving plots, it can use an explicit or global `project_dirs`
#' object (expected to be populated by `setupDirectories()`) and the successful
#' derivation of `time_dir` from it using `omic_type` and `experiment_label`.
#'
#' @param data The input data object. Can be a data frame (expected to conform
#'   to typical peptide or protein quantification structures) or an S4 object
#'   containing relevant slots (e.g., inheriting from `SummarizedExperiment`).
#'   The function attempts to automatically detect if it\'s peptide or protein data.
#' @param step_name A character string uniquely identifying the current filtering
#'   step (e.g., "InitialData", "FilteredByQuality", "Normalized"). This name is
#'   used for tracking in the `filtering_progress` object and plot labels.
#' @param omic_type Optional character string. The type of omics data
#'   (e.g., "proteomics", "metabolomics"). Used with `experiment_label` to
#'   derive save paths from the global `project_dirs` object. If `NULL` (default)
#'   or `experiment_label` is `NULL`, plots are not saved.
#' @param experiment_label Optional character string. The specific experiment
#'   label (e.g., "workshop_data"). Used with `omic_type` to derive save paths
#'   from the global `project_dirs` object. If `NULL` (default) or `omic_type`
#'   is `NULL`, plots are not saved.
#' @param overwrite Logical. If `TRUE`, allows overwriting an existing entry for
#'   `step_name` in the `filtering_progress` object. If `FALSE` (default) and
#'   `step_name` already exists, the function will stop with an error.
#' @param return_grid Logical. If `TRUE`, returns a single combined plot grid
#'   object created with `gridExtra::grid.arrange()`. If `FALSE` (default), prints
#'   individual plots and returns an invisible list of the ggplot objects.
#'
#' @return If `return_grid` is `TRUE`, returns a `grob` object (a grid graphical object).
#'   If `return_grid` is `FALSE`, returns an invisible list containing the individual
#'   `ggplot` objects (`proteins_total`, `proteins_per_run`, `peptides_total`,
#'   `peptides_per_protein`, `peptides_per_run`). Has side effects: modifies the
#'   global `filtering_progress` object and potentially saves plots to disk if
#'   `omic_type` and `experiment_label` are provided and paths are valid.
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_text labs theme_minimal theme element_text element_blank geom_line geom_point scale_color_manual annotate theme_void geom_boxplot coord_cartesian ggsave
#' @importFrom dplyr bind_rows mutate group_by ungroup %>%
#' @importFrom forcats fct_reorder
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @importFrom methods slotNames new is
#' @importFrom stats quantile
#' @importFrom multidplyr partition new_cluster cluster_library
#' @importFrom future plan
#'
#' @export
#' @param formats,project_dirs Runtime inputs used by this function; see the usage section for accepted values.
#' @param progress_env Environment that owns `filtering_progress`. Defaults to the
#'   global environment for backward compatibility.
updateProteinFiltering <- function(data, step_name,
                                   omic_type = NULL, experiment_label = NULL,
                                   overwrite = FALSE, return_grid = FALSE,
                                   formats = c("png", "pdf"),
                                   project_dirs = NULL,
                                   progress_env = .GlobalEnv) {
  if (!is.environment(progress_env)) {
    stop("`progress_env` must be an environment.", call. = FALSE)
  }
  # Initialize filtering_progress if it doesn\'t exist
  if (!exists("filtering_progress", envir = progress_env)) {
    filtering_progress <- new("FilteringProgress")
    assign("filtering_progress", filtering_progress, envir = progress_env)
  }

  # Get the current filtering_progress object
  filtering_progress <- get("filtering_progress", envir = progress_env)

  # Path derivation and save_plots logic
  derived_time_dir <- NULL
  save_plots <- FALSE

  if (!is.null(omic_type) && !is.null(experiment_label)) {
    # Check if project_dirs is provided as argument or exists in Global Env
    project_dirs_to_use <- if (!is.null(project_dirs)) {
      project_dirs
    } else if (exists("project_dirs", envir = .GlobalEnv)) {
      get("project_dirs", envir = .GlobalEnv)
    } else {
      NULL
    }

    if (is.null(project_dirs_to_use)) {
      warning("Object 'project_dirs' not found (neither as argument nor in Global Env). Plots will not be saved. Ensure 'setupDirectories()' has been run.")
    } else {
      project_dirs_global <- project_dirs_to_use
      omic_project_key <- paste0(omic_type, "_", experiment_label)

      if (!omic_project_key %in% names(project_dirs_global)) {
        warning(paste0("Entry for \'", omic_project_key, "\' not found in global \'project_dirs\'. Plots will not be saved."))
      } else {
        current_project_paths <- project_dirs_global[[omic_project_key]]
        if (is.null(current_project_paths)) {
          warning(paste0("Entry for \'", omic_project_key, "\' in global \'project_dirs\' is NULL. Plots will not be saved."))
        } else {
          derived_publication_graphs_dir <- current_project_paths$publication_graphs_dir
          temp_time_dir <- current_project_paths$time_dir

          if (is.null(temp_time_dir) || !is.character(temp_time_dir) || length(temp_time_dir) != 1 ||
            is.null(derived_publication_graphs_dir) || !is.character(derived_publication_graphs_dir) || length(derived_publication_graphs_dir) != 1) {
            warning(paste0(
              "\'time_dir\' or \'publication_graphs_dir\' is missing, not a character string, or not a single path for \'", omic_project_key,
              "\' in global \'project_dirs\'. Plots will not be saved."
            ))
          } else {
            if (!dir.exists(temp_time_dir)) {
              warning(paste0(
                "The derived \'time_dir\' (", temp_time_dir, ") for \'", omic_project_key,
                "\' does not exist. Plots will not be saved. Ensure directories are created via setupDirectories()."
              ))
            } else {
              derived_time_dir <- temp_time_dir
              save_plots <- TRUE
              message(paste0("Plots will be saved to: ", derived_time_dir))
            }
          }
        }
      }
    }
  } else {
    # Message if omic_type/label are missing and saving might have been expected
    if (return_grid && (is.null(omic_type) || is.null(experiment_label))) {
      message("omic_type and/or experiment_label not provided. Plots will not be saved.")
    }
  }

  protein_id_column <- resolveProteinCountColumn(data)

  # Determine if we\'re working with protein_quant_table
  is_protein_quant <- if (methods::is(data, "S4")) {
    "protein_quant_table" %in% slotNames(data)
  } else {
    isProteinQuantificationTable(data, protein_id_column)
  }

  # Calculate protein metrics (always done)
  protein_count <- countUniqueProteins(data, protein_id_column)
  proteins_per_run <- countProteinsPerRun(data, protein_id_column)

  # Ensure consistent data types in proteins_per_run
  proteins_per_run$Run <- as.character(proteins_per_run$Run)
  proteins_per_run$n_proteins <- as.numeric(proteins_per_run$n_proteins)

  # Update filtering progress based on data type
  if (step_name %in% filtering_progress@steps) {
    if (!overwrite) {
      stop("Step name \'", step_name, "\' already exists. Use overwrite = TRUE to replace it.")
    }
    idx <- which(filtering_progress@steps == step_name)

    # Always update protein metrics
    filtering_progress@proteins[idx] <- protein_count
    filtering_progress@proteins_per_run[[idx]] <- proteins_per_run

    if (!is_protein_quant) {
      # Update peptide metrics only for peptide data
      filtering_progress@total_peptides[idx] <- calcTotalPeptides(data, protein_id_column)
      peptides_per_protein <- calcPeptidesPerProtein(data, protein_id_column)
      peptides_per_run <- countPeptidesPerRun(data, protein_id_column)

      # Ensure consistent data types
      peptides_per_protein$Protein.Ids <- as.character(peptides_per_protein$Protein.Ids)
      peptides_per_protein$n_peptides <- as.numeric(peptides_per_protein$n_peptides)

      peptides_per_run$Run <- as.character(peptides_per_run$Run)
      peptides_per_run$n_peptides <- as.numeric(peptides_per_run$n_peptides)

      filtering_progress@peptides_per_protein[[idx]] <- peptides_per_protein
      filtering_progress@peptides_per_run[[idx]] <- peptides_per_run
    }
  } else {
    filtering_progress@steps <- c(filtering_progress@steps, step_name)
    filtering_progress@proteins <- c(filtering_progress@proteins, protein_count)
    filtering_progress@proteins_per_run <- c(
      filtering_progress@proteins_per_run,
      list(proteins_per_run)
    )

    if (!is_protein_quant) {
      # Add peptide metrics only for peptide data
      filtering_progress@total_peptides <- c(
        filtering_progress@total_peptides,
        calcTotalPeptides(data, protein_id_column)
      )

      peptides_per_protein <- calcPeptidesPerProtein(data, protein_id_column)
      peptides_per_run <- countPeptidesPerRun(data, protein_id_column)

      # Ensure consistent data types
      peptides_per_protein$Protein.Ids <- as.character(peptides_per_protein$Protein.Ids)
      peptides_per_protein$n_peptides <- as.numeric(peptides_per_protein$n_peptides)

      peptides_per_run$Run <- as.character(peptides_per_run$Run)
      peptides_per_run$n_peptides <- as.numeric(peptides_per_run$n_peptides)

      filtering_progress@peptides_per_protein <- c(
        filtering_progress@peptides_per_protein,
        list(peptides_per_protein)
      )
      filtering_progress@peptides_per_run <- c(
        filtering_progress@peptides_per_run,
        list(peptides_per_run)
      )
    } else {
      # For protein data, maintain existing peptide metrics or add NA/empty entries
      if (length(filtering_progress@total_peptides) > 0) {
        filtering_progress@total_peptides <- c(
          filtering_progress@total_peptides,
          filtering_progress@total_peptides[length(filtering_progress@total_peptides)]
        )
        filtering_progress@peptides_per_protein <- c(
          filtering_progress@peptides_per_protein,
          filtering_progress@peptides_per_protein[length(filtering_progress@peptides_per_protein)]
        )
        filtering_progress@peptides_per_run <- c(
          filtering_progress@peptides_per_run,
          filtering_progress@peptides_per_run[length(filtering_progress@peptides_per_run)]
        )
      } else {
        filtering_progress@total_peptides <- c(filtering_progress@total_peptides, NA_integer_)
        filtering_progress@peptides_per_protein <- c(
          filtering_progress@peptides_per_protein,
          list(data.frame(
            Protein.Ids = character(),
            n_peptides = integer()
          ))
        )
        filtering_progress@peptides_per_run <- c(
          filtering_progress@peptides_per_run,
          list(data.frame(
            Run = character(),
            n_peptides = integer()
          ))
        )
      }
    }
  }

  # Update the filtering_progress object in its owner environment
  assign("filtering_progress", filtering_progress, envir = progress_env)

  # Create base protein count plot (always shown)
  message("   [updateProteinFiltering] Generating P1 (Protein Count Bar)...")
  p1 <- ggplot(data.frame(
    step = factor(filtering_progress@steps, levels = filtering_progress@steps),
    proteins = filtering_progress@proteins
  ), aes(x = step, y = proteins)) +
    geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
    geom_text(aes(label = proteins),
      vjust = -0.5,
      size = 4
    ) +
    labs(
      title = "Number of Proteins",
      x = "Filtering Step",
      y = "Unique Proteins"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank()
    )

  # Create proteins per run plot (always shown)
  message("   [updateProteinFiltering] Generating P4 (Proteins Per Run Line)...")
  # First ensure all data frames in the list have consistent column types
  proteins_per_run_list <- lapply(filtering_progress@proteins_per_run, function(df) {
    df$Run <- as.character(df$Run)
    df$n_proteins <- as.numeric(df$n_proteins)
    return(df)
  })

  message(sprintf("      [P4] Binding %d data frames...", length(proteins_per_run_list)))
  p4_data <- bind_rows(proteins_per_run_list, .id = "step")
  message(sprintf("      [P4] Combined data rows: %d", nrow(p4_data)))

  p4 <- p4_data |>
    mutate(step = filtering_progress@steps[as.numeric(step)]) |>
    group_by(Run) |>
    mutate(avg_proteins = mean(n_proteins)) |>
    ungroup() |>
    # Run is already character from our preprocessing
    mutate(Run = fct_reorder(Run, avg_proteins)) |>
    ggplot(aes(
      x = Run, y = n_proteins,
      group = step,
      color = factor(step, levels = filtering_progress@steps)
    )) +
    geom_line() +
    geom_point() +
    labs(
      title = "Proteins per Run",
      x = "Run ID (ordered by average protein count)",
      y = "Number of Proteins",
      color = "Step"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank()
    ) +
    scale_color_manual(values = get_color_palette(length(filtering_progress@steps), "steelblue"))

  message("   [updateProteinFiltering] P4 Generated.")

  # Initialize peptide plots
  if (is_protein_quant) {
    # For protein data, create empty placeholder plots if no peptide data exists
    if (all(is.na(filtering_progress@total_peptides))) {
      p2 <- p3 <- p5 <- ggplot() +
        annotate("text",
          x = 0.5, y = 0.5,
          label = "No peptide data available for protein quantification data"
        ) +
        theme_void()
    } else {
      # If peptide data exists from previous steps, create plots with existing data
      p2 <- ggplot(data.frame(
        step = factor(filtering_progress@steps, levels = filtering_progress@steps),
        total_peptides = filtering_progress@total_peptides
      ), aes(x = step, y = total_peptides)) +
        geom_bar(stat = "identity", fill = "forestgreen", width = 0.7) +
        geom_text(aes(label = total_peptides),
          vjust = -0.5,
          size = 4
        ) +
        labs(
          title = "Total Unique Peptides (from last peptide data)",
          x = "Filtering Step",
          y = "Unique Peptides"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major.x = element_blank()
        )

      # Ensure consistent data types in peptides_per_protein list
      peptides_per_protein_list <- lapply(filtering_progress@peptides_per_protein, function(df) {
        if (nrow(df) > 0) {
          df$Protein.Ids <- as.character(df$Protein.Ids)
          df$n_peptides <- as.numeric(df$n_peptides)
        }
        return(df)
      })

      p3 <- ggplot() +
        geom_boxplot(
          data = bind_rows(peptides_per_protein_list, .id = "step") |>
            mutate(step = filtering_progress@steps[as.numeric(step)]),
          aes(
            x = factor(step, levels = filtering_progress@steps),
            y = n_peptides
          ),
          fill = "darkred",
          alpha = 0.5,
          outlier.shape = NA
        ) +
        labs(
          title = "Peptides per Protein Distribution (from last peptide data)",
          x = "Filtering Step",
          y = "Number of Peptides"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major.x = element_blank()
        ) +
        coord_cartesian(
          ylim = c(
            0,
            quantile(bind_rows(peptides_per_protein_list)$n_peptides, 0.95)
          )
        )

      # Ensure consistent data types in peptides_per_run list
      peptides_per_run_list <- lapply(filtering_progress@peptides_per_run, function(df) {
        if (nrow(df) > 0) {
          df$Run <- as.character(df$Run)
          df$n_peptides <- as.numeric(df$n_peptides)
        }
        return(df)
      })

      p5 <- bind_rows(peptides_per_run_list, .id = "step") |>
        mutate(step = filtering_progress@steps[as.numeric(step)]) |>
        group_by(Run) |>
        mutate(avg_peptides = mean(n_peptides)) |>
        ungroup() |>
        # Run is already character from our preprocessing
        mutate(Run = fct_reorder(Run, avg_peptides)) |>
        ggplot(aes(
          x = Run, y = n_peptides,
          group = step,
          color = factor(step, levels = filtering_progress@steps)
        )) +
        geom_line() +
        geom_point() +
        labs(
          title = "Peptides per Run (from last peptide data)",
          x = "Run ID (ordered by average peptide count)",
          y = "Number of Peptides",
          color = "Step"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major.x = element_blank()
        ) +
        scale_color_manual(values = get_color_palette(length(filtering_progress@steps), "forestgreen"))
    }
  } else {
    # For peptide data, create normal plots
    p2 <- ggplot(data.frame(
      step = factor(filtering_progress@steps, levels = filtering_progress@steps),
      total_peptides = filtering_progress@total_peptides
    ), aes(x = step, y = total_peptides)) +
      geom_bar(stat = "identity", fill = "forestgreen", width = 0.7) +
      geom_text(aes(label = total_peptides),
        vjust = -0.5,
        size = 4
      ) +
      labs(
        title = "Total Unique Peptides",
        x = "Filtering Step",
        y = "Unique Peptides"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank()
      )

    # Ensure consistent data types in peptides_per_protein list
    peptides_per_protein_list <- lapply(filtering_progress@peptides_per_protein, function(df) {
      if (nrow(df) > 0) {
        df$Protein.Ids <- as.character(df$Protein.Ids)
        df$n_peptides <- as.numeric(df$n_peptides)
      }
      return(df)
    })

    p3 <- ggplot() +
      geom_boxplot(
        data = bind_rows(peptides_per_protein_list, .id = "step") |>
          mutate(step = filtering_progress@steps[as.numeric(step)]),
        aes(
          x = factor(step, levels = filtering_progress@steps),
          y = n_peptides
        ),
        fill = "darkred",
        alpha = 0.5,
        outlier.shape = NA
      ) +
      labs(
        title = "Peptides per Protein Distribution",
        x = "Filtering Step",
        y = "Number of Peptides"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank()
      ) +
      coord_cartesian(
        ylim = c(
          0,
          quantile(bind_rows(peptides_per_protein_list)$n_peptides, 0.95)
        )
      )

    # Ensure consistent data types in peptides_per_run list
    peptides_per_run_list <- lapply(filtering_progress@peptides_per_run, function(df) {
      if (nrow(df) > 0) {
        df$Run <- as.character(df$Run)
        df$n_peptides <- as.numeric(df$n_peptides)
      }
      return(df)
    })

    p5 <- bind_rows(peptides_per_run_list, .id = "step") |>
      mutate(step = filtering_progress@steps[as.numeric(step)]) |>
      group_by(Run) |>
      mutate(avg_peptides = mean(n_peptides)) |>
      ungroup() |>
      # Run is already character from our preprocessing
      mutate(Run = fct_reorder(Run, avg_peptides)) |>
      ggplot(aes(
        x = Run, y = n_peptides,
        group = step,
        color = factor(step, levels = filtering_progress@steps)
      )) +
      geom_line() +
      geom_point() +
      labs(
        title = "Peptides per Run",
        x = "Run ID (ordered by average peptide count)",
        y = "Number of Peptides",
        color = "Step"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank()
      ) +
      scale_color_manual(values = get_color_palette(length(filtering_progress@steps), "forestgreen"))
  }

  # Create plot list based on data type
  plot_list <- list(
    proteins_total = p1,
    proteins_per_run = p4,
    peptides_total = p2,
    peptides_per_protein = p3,
    peptides_per_run = p5
  )

  # Save plots if derived_time_dir is valid and save_plots is TRUE
  if (save_plots) {
    message(sprintf("   [updateProteinFiltering] Saving individual plots to %s in formats: %s...", derived_time_dir, paste(formats, collapse = ", ")))
    for (plot_name in names(plot_list)) {
      message(sprintf("      Saving %s...", plot_name))
      for (fmt in formats) {
        filename <- file.path(
          derived_time_dir,
          sprintf("%s_%s.%s", step_name, plot_name, fmt)
        )
        tryCatch(
          {
            ggsave(filename,
              plot = plot_list[[plot_name]],
              width = 10,
              height = 8,
              dpi = 300
            )
          },
          error = function(e) message(sprintf("Warning: Failed to save %s as %s: %s", plot_name, fmt, e$message))
        )
      }
    }
  }

  # Return/display plots based on return_grid parameter
  if (return_grid) {
    message("   [updateProteinFiltering] Generating final grid...")
    message(sprintf("      Memory before grid: %s", format(sum(gc()[, 2]), units = "auto")))

    tryCatch(
      {
        if (!is_protein_quant || !all(is.na(filtering_progress@total_peptides))) {
          # Create full grid with all plots if peptide data exists
          message("      Combining all 5 plots...")
          grid1 <- gridExtra::arrangeGrob(p1, p2, p3, ncol = 3)
          grid2 <- gridExtra::arrangeGrob(p4, ncol = 1)
          grid3 <- gridExtra::arrangeGrob(p5, ncol = 1)

          # Use arrangeGrob to prevent immediate drawing, which might double-render
          grid_plot <- gridExtra::arrangeGrob(
            grid1,
            grid2,
            grid3,
            heights = c(1, 1, 1)
          )
        } else {
          # For protein_quant_table without peptide data, only show protein plots
          message("      Combining protein plots (p1, p4)...")
          grid_plot <- gridExtra::arrangeGrob(
            p1,
            p4,
            ncol = 1,
            heights = c(1, 1)
          )
        }

        message(sprintf("      Memory after grid creation: %s", format(sum(gc()[, 2]), units = "auto")))

        # Save the grid if derived_time_dir is valid and save_plots is TRUE
        if (save_plots) {
          message(sprintf("      Saving combined grid plot in formats: %s...", paste(formats, collapse = ", ")))
          for (fmt in formats) {
            filename <- file.path(
              derived_time_dir,
              sprintf("%s_combined_plots.%s", step_name, fmt)
            )
            ggsave(filename,
              plot = grid_plot,
              width = 15,
              height = if (!is_protein_quant || !all(is.na(filtering_progress@total_peptides))) 18 else 12,
              dpi = 300
            )
          }
        }

        if (!is.null(grid_plot)) {
          gridExtra::grid.arrange(grid_plot)
        }
        return(invisible(grid_plot))
      },
      error = function(e) {
        message(sprintf("ERROR in grid generation: %s", e$message))
        return(NULL)
      }
    )
  } else {
    # Print each plot individually
    message("   [updateProteinFiltering] Printing individual plots...")
    for (plot_obj in plot_list) {
      print(plot_obj)
    }
    # Return the list invisibly
    invisible(plot_list)
  }
}
