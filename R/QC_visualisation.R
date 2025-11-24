#' @title Update and Visualize Filtering Progress
#' @description Tracks and visualizes the impact of filtering steps on peptide 
#'   and protein counts. Updates a global `FilteringProgress` object and optionally 
#'   saves plots summarizing the changes. Handles both peptide-level and 
#'   protein-level data inputs.
#' 
#' @details 
#' This function acts as a central hub for monitoring data reduction throughout 
#' a filtering workflow. It performs the following actions:
#' \itemize{
#'   \item Initializes or retrieves a global S4 object named `filtering_progress` 
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
#' **Important:** This function relies on and modifies a global variable named 
#' `filtering_progress`. For saving plots, it depends on the global `project_dirs` 
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
#' 
#' @export
updateProteinFiltering <- function(data, step_name, 
                                 omic_type = NULL, experiment_label = NULL,
                                 overwrite = FALSE, return_grid = FALSE) {
    
    # Initialize filtering_progress if it doesn\'t exist
    if (!exists("filtering_progress", envir = .GlobalEnv)) {
        filtering_progress <- new("FilteringProgress")
        assign("filtering_progress", filtering_progress, envir = .GlobalEnv)
    }
    
    # Get the current filtering_progress object
    filtering_progress <- get("filtering_progress", envir = .GlobalEnv)
    
    # Path derivation and save_plots logic
    derived_time_dir <- NULL
    save_plots <- FALSE

    if (!is.null(omic_type) && !is.null(experiment_label)) {
        if (!exists("project_dirs", envir = .GlobalEnv)) {
            warning("Global object \'project_dirs\' not found. Plots will not be saved. Ensure \'setupDirectories()\' has been run.")
        } else {
            project_dirs_global <- get("project_dirs", envir = .GlobalEnv)
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
                        warning(paste0("\'time_dir\' or \'publication_graphs_dir\' is missing, not a character string, or not a single path for \'", omic_project_key,
                                       "\' in global \'project_dirs\'. Plots will not be saved."))
                    } else {
                        if (!dir.exists(temp_time_dir)) {
                            warning(paste0("The derived \'time_dir\' (", temp_time_dir, ") for \'", omic_project_key,
                                           "\' does not exist. Plots will not be saved. Ensure directories are created via setupDirectories()."))
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
    
    # Determine if we\'re working with protein_quant_table
    is_protein_quant <- if (methods::is(data, "S4")) {
        "protein_quant_table" %in% slotNames(data)
    } else {
        # For data frames, check if it looks like a protein quant table
        if ("Protein.Ids" %in% names(data)) {
            all(sapply(data[setdiff(names(data), "Protein.Ids")], is.numeric))
        } else {
            FALSE
        }
    }
    
    # Calculate protein metrics (always done)
    protein_count <- countUniqueProteins(data)
    proteins_per_run <- countProteinsPerRun(data)
    
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
            filtering_progress@total_peptides[idx] <- calcTotalPeptides(data)
            peptides_per_protein <- calcPeptidesPerProtein(data)
            peptides_per_run <- countPeptidesPerRun(data)
            
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
        filtering_progress@proteins_per_run <- c(filtering_progress@proteins_per_run, 
                                               list(proteins_per_run))
        
        if (!is_protein_quant) {
            # Add peptide metrics only for peptide data
            filtering_progress@total_peptides <- c(filtering_progress@total_peptides, 
                                                 calcTotalPeptides(data))
            
            peptides_per_protein <- calcPeptidesPerProtein(data)
            peptides_per_run <- countPeptidesPerRun(data)
            
            # Ensure consistent data types
            peptides_per_protein$Protein.Ids <- as.character(peptides_per_protein$Protein.Ids)
            peptides_per_protein$n_peptides <- as.numeric(peptides_per_protein$n_peptides)
            
            peptides_per_run$Run <- as.character(peptides_per_run$Run)
            peptides_per_run$n_peptides <- as.numeric(peptides_per_run$n_peptides)
            
            filtering_progress@peptides_per_protein <- c(filtering_progress@peptides_per_protein, 
                                                       list(peptides_per_protein))
            filtering_progress@peptides_per_run <- c(filtering_progress@peptides_per_run, 
                                                   list(peptides_per_run))
        } else {
            # For protein data, maintain existing peptide metrics or add NA/empty entries
            if (length(filtering_progress@total_peptides) > 0) {
                filtering_progress@total_peptides <- c(filtering_progress@total_peptides, 
                                                     filtering_progress@total_peptides[length(filtering_progress@total_peptides)])
                filtering_progress@peptides_per_protein <- c(filtering_progress@peptides_per_protein, 
                                                           filtering_progress@peptides_per_protein[length(filtering_progress@peptides_per_protein)])
                filtering_progress@peptides_per_run <- c(filtering_progress@peptides_per_run, 
                                                       filtering_progress@peptides_per_run[length(filtering_progress@peptides_per_run)])
            } else {
                filtering_progress@total_peptides <- c(filtering_progress@total_peptides, NA_integer_)
                filtering_progress@peptides_per_protein <- c(filtering_progress@peptides_per_protein, 
                                                           list(data.frame(Protein.Ids = character(), 
                                                                         n_peptides = integer())))
                filtering_progress@peptides_per_run <- c(filtering_progress@peptides_per_run, 
                                                       list(data.frame(Run = character(), 
                                                                     n_peptides = integer())))
            }
        }
    }
    
    # Update the global filtering_progress object
    assign("filtering_progress", filtering_progress, envir = .GlobalEnv)
    
    # Create base protein count plot (always shown)
    p1 <- ggplot(data.frame(
        step = factor(filtering_progress@steps, levels = filtering_progress@steps),
        proteins = filtering_progress@proteins
    ), aes(x = step, y = proteins)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
        geom_text(aes(label = proteins), 
                  vjust = -0.5, 
                  size = 4) +
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
    # First ensure all data frames in the list have consistent column types
    proteins_per_run_list <- lapply(filtering_progress@proteins_per_run, function(df) {
        df$Run <- as.character(df$Run)
        df$n_proteins <- as.numeric(df$n_proteins)
        return(df)
    })
    
    p4 <- bind_rows(proteins_per_run_list, .id = "step") |>
        mutate(step = filtering_progress@steps[as.numeric(step)]) |>
        group_by(Run) |>
        mutate(avg_proteins = mean(n_proteins)) |>
        ungroup() |>
        # Run is already character from our preprocessing
        mutate(Run = fct_reorder(Run, avg_proteins)) |>
        ggplot(aes(x = Run, y = n_proteins, 
                  group = step, 
                  color = factor(step, levels = filtering_progress@steps))) +
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
    
    # Initialize peptide plots
    if (is_protein_quant) {
        # For protein data, create empty placeholder plots if no peptide data exists
        if (all(is.na(filtering_progress@total_peptides))) {
            p2 <- p3 <- p5 <- ggplot() + 
                annotate("text", x = 0.5, y = 0.5, 
                        label = "No peptide data available for protein quantification data") +
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
                          size = 4) +
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
                geom_boxplot(data = bind_rows(peptides_per_protein_list, .id = "step") |>
                             mutate(step = filtering_progress@steps[as.numeric(step)]),
                           aes(x = factor(step, levels = filtering_progress@steps), 
                               y = n_peptides),
                           fill = "darkred",
                           alpha = 0.5,
                           outlier.shape = NA) +
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
                    ylim = c(0, 
                             quantile(bind_rows(peptides_per_protein_list)$n_peptides, 0.95))
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
                ggplot(aes(x = Run, y = n_peptides, 
                          group = step, 
                          color = factor(step, levels = filtering_progress@steps))) +
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
                      size = 4) +
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
            geom_boxplot(data = bind_rows(peptides_per_protein_list, .id = "step") |>
                         mutate(step = filtering_progress@steps[as.numeric(step)]),
                       aes(x = factor(step, levels = filtering_progress@steps), 
                           y = n_peptides),
                       fill = "darkred",
                       alpha = 0.5,
                       outlier.shape = NA) +
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
                ylim = c(0, 
                         quantile(bind_rows(peptides_per_protein_list)$n_peptides, 0.95))
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
            ggplot(aes(x = Run, y = n_peptides, 
                      group = step, 
                      color = factor(step, levels = filtering_progress@steps))) +
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
        for (plot_name in names(plot_list)) {
            filename <- file.path(derived_time_dir,
                                sprintf("%s_%s.png", step_name, plot_name))
            ggsave(filename, 
                   plot = plot_list[[plot_name]], 
                   width = 10, 
                   height = 8, 
                   dpi = 300)
        }
    }
    
    # Return/display plots based on return_grid parameter
    if(return_grid) {
        if (!is_protein_quant || !all(is.na(filtering_progress@total_peptides))) {
            # Create full grid with all plots if peptide data exists
            grid1 <- gridExtra::arrangeGrob(p1, p2, p3, ncol = 3)
            grid2 <- gridExtra::arrangeGrob(p4, ncol = 1)
            grid3 <- gridExtra::arrangeGrob(p5, ncol = 1)
            
            grid_plot <- gridExtra::grid.arrange(
                grid1,
                grid2,
                grid3,
                heights = c(1, 1, 1)
            )
        } else {
            # For protein_quant_table without peptide data, only show protein plots
            grid_plot <- gridExtra::grid.arrange(
                p1,
                p4,
                ncol = 1,
                heights = c(1, 1)
            )
        }
        
        # Save the grid if derived_time_dir is valid and save_plots is TRUE
        if (save_plots) {
            filename <- file.path(derived_time_dir,
                                sprintf("%s_combined_plots.png", step_name))
            ggsave(filename, 
                   plot = grid_plot, 
                   width = 15, 
                   height = if (!is_protein_quant || !all(is.na(filtering_progress@total_peptides))) 18 else 12,
                   dpi = 300)
        }
        
        return(grid_plot)
    } else {
        # Print each plot individually
        for(plot_obj in plot_list) { # Changed loop variable to avoid conflict with base::plot
            print(plot_obj)
        }
        # Return the list invisibly
        invisible(plot_list)
    }
}

# ---- Metabolomics Filtering Progress Tracking ----

#' @title S4 Class Definition for Metabolomics Filtering Progress
#' @description An S4 class to store QC metrics collected at different steps
#'              of a metabolomics data processing workflow.
#'
#' @slot steps A character vector storing the names of the filtering steps.
#' @slot assay_names A list where each element corresponds to a step and contains
#'                 a character vector of the original assay names for that step.
#' @slot n_metabolites_per_assay A list where each element corresponds to a step.
#'                                Each element is a list itself, containing named
#'                                numeric vectors (assay_name = count of unique metabolites).
#' @slot n_metabolites_total A numeric vector storing the total unique metabolites
#'                            (across all assays) identified at each step.
#' @slot detected_per_sample A list where each element corresponds to a step.
#'                            Each element is a list (per assay) containing data frames
#'                            with columns for sample ID (e.g., 'Run') and the number
#'                            of detected metabolites ('n_detected').
#' @slot missingness_per_assay A list where each element corresponds to a step.
#'                               Each element is a list containing named numeric
#'                               vectors (assay_name = percentage missing).
#' @slot sum_intensity_per_sample A list where each element corresponds to a step.
#'                                 Each element is a list (per assay) containing data frames
#'                                 with columns for sample ID ('Run') and the total
#'                                 intensity ('sum_intensity').
#' @slot cv_distribution_per_assay A list where each element corresponds to a step.
#'                                   Each element is a list (per assay) containing data frames
#'                                   with columns for metabolite ID, group, and calculated
#'                                   coefficient of variation ('cv').
#' @slot is_metrics_per_assay A list where each element corresponds to a step.
#'                             Each element is a list (per assay) containing data frames
#'                             summarizing internal standard metrics (e.g., 'is_id',
#'                             'metric_name', 'value').
#'
#' @keywords internal
#' @import methods
#' @exportClass FilteringProgressMetabolomics
setClass("FilteringProgressMetabolomics",
    slots = c(
        steps = "character",
        assay_names = "list",
        n_metabolites_per_assay = "list",
        n_metabolites_total = "numeric",
        detected_per_sample = "list",
        missingness_per_assay = "list",
        sum_intensity_per_sample = "list",
        cv_distribution_per_assay = "list",
        is_metrics_per_assay = "list"
    )
)

#' @title Initialize or Retrieve Global Metabolomics Filtering Progress Object
#' @description Checks for a global object named `filtering_progress_metabolomics`
#'              of class `FilteringProgressMetabolomics`. If it doesn't exist,
#'              it creates and assigns a new one to the global environment.
#'
#' @return The global `FilteringProgressMetabolomics` object.
#' @keywords internal
#' @noRd
#' @export
getFilteringProgressMetabolomics <- function() {
    if (!exists("filtering_progress_metabolomics", envir = .GlobalEnv)) {
        filtering_progress_metabolomics <- new("FilteringProgressMetabolomics")
        assign("filtering_progress_metabolomics", filtering_progress_metabolomics, envir = .GlobalEnv)
        message("Initialized global 'filtering_progress_metabolomics' object.") # Optional message
    }
    get("filtering_progress_metabolomics", envir = .GlobalEnv)
}

#' @title Update the Global Metabolomics Filtering Progress Object
#' @description Modifies the global `filtering_progress_metabolomics` object
#'              by adding or overwriting data for a specific step.
#'
#' @param prog_met The `FilteringProgressMetabolomics` object (retrieved globally).
#' @param step_name The name of the step to add or update.
#' @param current_assay_names Character vector of assay names for this step.
#' @param metrics_list A nested list containing the calculated metrics for each assay
#'                      for the current step.
#' @param total_metabolites The total unique metabolites calculated across assays for this step.
#' @param overwrite Logical, whether to overwrite if `step_name` exists.
#'
#' @return The updated `FilteringProgressMetabolomics` object (invisibly).
#'         Has the side effect of updating the global object.
#' @keywords internal
#' @noRd
#' @export
updateFilteringProgressMetabolomics <- function(prog_met,
                                                  step_name,
                                                  current_assay_names,
                                                  metrics_list,
                                                  total_metabolites,
                                                  overwrite = FALSE) {

    if (step_name %in% prog_met@steps) {
        if (!overwrite) {
            stop("Step name '", step_name, "' already exists in filtering_progress_metabolomics. Use overwrite = TRUE to replace it.")
        }
        idx <- which(prog_met@steps == step_name)

        # Overwrite existing data
        prog_met@assay_names[[idx]] <- current_assay_names
        prog_met@n_metabolites_per_assay[[idx]] <- lapply(metrics_list, `[[`, "n_metabolites")
        prog_met@n_metabolites_total[idx] <- total_metabolites
        prog_met@detected_per_sample[[idx]] <- lapply(metrics_list, `[[`, "detected_per_sample")
        prog_met@missingness_per_assay[[idx]] <- lapply(metrics_list, `[[`, "missingness")
        prog_met@sum_intensity_per_sample[[idx]] <- lapply(metrics_list, `[[`, "sum_intensity_per_sample")
        prog_met@cv_distribution_per_assay[[idx]] <- lapply(metrics_list, `[[`, "cv_distribution")
        prog_met@is_metrics_per_assay[[idx]] <- lapply(metrics_list, `[[`, "is_metrics")

    } else {
        # Append new data
        prog_met@steps <- c(prog_met@steps, step_name)
        prog_met@assay_names <- c(prog_met@assay_names, list(current_assay_names))
        prog_met@n_metabolites_per_assay <- c(prog_met@n_metabolites_per_assay, list(lapply(metrics_list, `[[`, "n_metabolites")))
        prog_met@n_metabolites_total <- c(prog_met@n_metabolites_total, total_metabolites)
        prog_met@detected_per_sample <- c(prog_met@detected_per_sample, list(lapply(metrics_list, `[[`, "detected_per_sample")))
        prog_met@missingness_per_assay <- c(prog_met@missingness_per_assay, list(lapply(metrics_list, `[[`, "missingness")))
        prog_met@sum_intensity_per_sample <- c(prog_met@sum_intensity_per_sample, list(lapply(metrics_list, `[[`, "sum_intensity_per_sample")))
        prog_met@cv_distribution_per_assay <- c(prog_met@cv_distribution_per_assay, list(lapply(metrics_list, `[[`, "cv_distribution")))
        prog_met@is_metrics_per_assay <- c(prog_met@is_metrics_per_assay, list(lapply(metrics_list, `[[`, "is_metrics")))
    }

    # Update the global object
    assign("filtering_progress_metabolomics", prog_met, envir = .GlobalEnv)
    invisible(prog_met)
}

# ---- Metabolomics QC Metric Helper Functions ----

#' @title Extract Quantitative Data and Sample Names from Assay Tibble
#' @description Separates annotation columns from quantitative data columns
#'              in a metabolomics assay tibble.
#'
#' @param assay_data A tibble/data.frame representing one metabolomics assay,
#'                   with metabolite annotations and sample intensity columns.
#' @param sample_id_col The name of the column in the design matrix that contains
#'                      the sample identifiers (which match column names in assay_data).
#'                      ***Note: This argument is currently unused but kept for consistency.
#'                      The function infers sample columns based on numeric type.***
#'
#' @return A list containing:
#'         - `quant_data`: A data frame with only the quantitative (numeric) columns.
#'         - `sample_names`: A character vector of the sample column names.
#'         - `annotation_data`: A data frame with the non-numeric annotation columns.
#'
#' @importFrom dplyr select where
#' @keywords internal
#' @noRd
#' @export 
getMetaboliteQuantData <- function(assay_data) {
    # Identify quantitative (numeric) columns - assumes sample columns are numeric
    quant_cols <- sapply(assay_data, is.numeric)
    quant_data <- assay_data[, quant_cols, drop = FALSE]
    sample_names <- colnames(quant_data)

    # Identify annotation (non-numeric) columns
    annotation_data <- assay_data[, !quant_cols, drop = FALSE]

    return(list(
        quant_data = as.data.frame(quant_data), # Convert tibble to df if needed downstream
        sample_names = sample_names,
        annotation_data = as.data.frame(annotation_data)
    ))
}


#' @title Count Unique Metabolites in an Assay
#' @description Counts the number of unique, non-NA identifiers in the specified
#'              metabolite ID column of an assay tibble.
#'
#' @param assay_data A tibble/data.frame for one assay.
#' @param metabolite_id_col Character string, the name of the column containing
#'                         metabolite identifiers.
#'
#' @return A single numeric value: the count of unique metabolites.
#'         Returns 0 if the column doesn't exist or has no non-NA values.
#'
#' @importFrom dplyr pull distinct n
#' @keywords internal
#' @noRd
#' @export 
countUniqueMetabolites <- function(assay_data, metabolite_id_col) {
    if (!metabolite_id_col %in% colnames(assay_data)) {
        warning("Metabolite ID column '", metabolite_id_col, "' not found in assay data.")
        return(0)
    }
    
    # Count unique non-NA metabolite IDs - avoid using pipe operators
    metabolite_ids <- dplyr::pull(assay_data, .data[[metabolite_id_col]])
    metabolite_ids <- stats::na.omit(metabolite_ids)
    unique_ids <- dplyr::distinct(data.frame(id = metabolite_ids))
    return(nrow(unique_ids))
}

#' @title Count Detected Metabolites per Sample
#' @description For each sample column in an assay's quantitative data,
#'              counts the number of non-missing, non-zero intensity values.
#'
#' @param assay_data A tibble/data.frame for one assay.
#' @param sample_id_col ***Currently unused. Sample columns inferred by numeric type.***
#' @param metabolite_id_col ***Currently unused.***
#'
#' @return A data frame with columns 'Run' (sample name) and 'n_detected'.
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr group_by summarise n filter
#' @keywords internal
#' @noRd
#' @export
countMetabolitesPerSample <- function(assay_data, sample_id_col, metabolite_id_col) {
    quant_info <- getMetaboliteQuantData(assay_data)
    quant_data <- quant_info$quant_data
    sample_names <- quant_info$sample_names

    if (length(sample_names) == 0) {
        return(data.frame(Run = character(), n_detected = integer()))
    }

    # Add a temporary row identifier if no metabolite ID column exists or is non-unique
    # This ensures pivot_longer works correctly even without a unique ID col.
    quant_data$..temp_row_id.. <- seq_len(nrow(quant_data))

    quant_data |>
        tidyr::pivot_longer(
            cols = all_of(sample_names),
            names_to = "Run",
            values_to = "Intensity"
        ) |>
        # Define 'detected' as non-NA and greater than 0 (adjust if needed)
        dplyr::filter(!is.na(.data$Intensity) & .data$Intensity > 0) |>
        dplyr::group_by(.data$Run) |>
        dplyr::summarise(n_detected = dplyr::n(), .groups = "drop")
}

#' @title Calculate Overall Missingness Percentage
#' @description Calculates the percentage of missing values (NA or zero)
#'              in the quantitative portion of an assay tibble.
#'
#' @param assay_data A tibble/data.frame for one assay.
#' @param sample_id_col ***Currently unused. Sample columns inferred by numeric type.***
#'
#' @return A single numeric value: the percentage of missing data points.
#'
#' @keywords internal
#' @noRd
#' @export
calculateMissingness <- function(assay_data, sample_id_col) {
    # Get only the quantitative columns (sample columns)
    quant_info <- getMetaboliteQuantData(assay_data)
    quant_data <- quant_info$quant_data
    sample_names <- quant_info$sample_names

    # Validate input
    if (nrow(quant_data) == 0 || ncol(quant_data) == 0 || length(sample_names) == 0) {
        warning("No valid data for missingness calculation")
        return(NA_real_)
    }

    # Exclude any non-sample columns
    if ("..temp_row_id.." %in% colnames(quant_data)) {
        quant_data <- quant_data[, setdiff(colnames(quant_data), "..temp_row_id.."), drop = FALSE]
    }
    
    # Count missing (NA or zero) values in all sample columns
    total_cells <- nrow(quant_data) * length(sample_names)
    missing_values <- 0
    
    for (col in sample_names) {
        missing_values <- missing_values + sum(is.na(quant_data[[col]]) | quant_data[[col]] == 0)
    }
    
    # Calculate percentage
    missing_pct <- (missing_values / total_cells) * 100
    
    # Debug output to verify calculation
    message("DEBUG: Missing values: ", missing_values, ", Total cells: ", total_cells, 
            ", Percentage: ", missing_pct, "%")
    
    return(missing_pct)
}

#' @title Calculate Sum Intensity per Sample (TIC Proxy)
#' @description Calculates the sum of all intensities for each sample column
#'              in the quantitative portion of an assay tibble.
#'
#' @param assay_data A tibble/data.frame for one assay.
#' @param sample_id_col ***Currently unused. Sample columns inferred by numeric type.***
#'
#' @return A data frame with columns 'Run' (sample name) and 'sum_intensity'.
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr group_by summarise
#' @keywords internal
#' @noRd
#' @export
calculateSumIntensityPerSample <- function(assay_data, sample_id_col) {
    quant_info <- getMetaboliteQuantData(assay_data)
    quant_data <- quant_info$quant_data
    sample_names <- quant_info$sample_names

    if (length(sample_names) == 0) {
        return(data.frame(Run = character(), sum_intensity = numeric()))
    }

    # Add a temporary row identifier
    quant_data$..temp_row_id.. <- seq_len(nrow(quant_data))

    quant_data |>
        tidyr::pivot_longer(
            cols = all_of(sample_names),
            names_to = "Run",
            values_to = "Intensity"
        ) |>
        dplyr::group_by(.data$Run) |>
        # Sum intensities, replacing NA with 0 for the sum
        dplyr::summarise(sum_intensity = sum(.data$Intensity, na.rm = TRUE), .groups = "drop")
}

#' @title Calculate Total Unique Metabolites Across Assays
#' @description Finds all unique metabolite IDs present across a list of assay tibbles.
#'
#' @param assay_list A list of tibbles/data.frames, one for each assay.
#' @param metabolite_id_col Character string, the name of the column containing
#'                         metabolite identifiers in each assay tibble.
#'
#' @return A single numeric value: the total count of unique metabolite IDs
#'         across all provided assays.
#'
#' @importFrom purrr map
#' @importFrom dplyr distinct n pull
#' @keywords internal
#' @noRd
#' @export
calculateTotalUniqueMetabolitesAcrossAssays <- function(assay_list, metabolite_id_col) {
    if (length(assay_list) == 0) {
        return(0)
    }

    # Extract the metabolite ID column from each assay, handling missing columns
    all_ids <- purrr::map(assay_list, ~{
            if (metabolite_id_col %in% colnames(.x)) {
                dplyr::pull(.x, {{ metabolite_id_col }})
            } else {
                NULL # Return NULL if column is missing
            }
        }) |>
        unlist() # Combine all IDs into a single vector
    
    # Remove NAs
    all_ids <- stats::na.omit(all_ids)
    
    # If no valid IDs found, return 0
    if (length(all_ids) == 0) {
        return(0)
    }
    
    # Count unique IDs using a more explicit approach
    unique_ids <- dplyr::distinct(data.frame(id = all_ids))
    return(nrow(unique_ids))
}

#' @title Calculate Within-Group Metabolite CVs
#' @description Calculates the Coefficient of Variation (CV) for each metabolite
#'              across replicate samples within each experimental group.
#'
#' @param assay_data A tibble/data.frame for one assay.
#' @param design_matrix A data frame linking samples to experimental design
#'                      (must contain sample_id_col and group_id_col).
#' @param group_id_col Character string, the name of the grouping column in design_matrix.
#' @param replicate_id_col Character string, the name of the replicate identifier column
#'                         in design_matrix (used indirectly to ensure grouping works).
#' @param sample_id_col Character string, the name of the sample ID column in design_matrix
#'                      (matching column names in assay_data).
#' @param metabolite_id_col Character string, the name of the metabolite ID column in assay_data.
#'
#' @return A data frame with columns 'metabolite_id', 'group', and 'cv'.
#'         Returns an empty data frame if inputs are invalid or calculations fail.
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join group_by summarise filter select rename all_of
#' @importFrom stats sd na.omit
#' @keywords internal
#' @noRd
#' @export
calculateMetaboliteCVs <- function(assay_data,
                                   design_matrix,
                                   group_id_col,
                                   replicate_id_col, # Keep for consistency, maybe use later
                                   sample_id_col,
                                   metabolite_id_col) {

    # --- Input Validation --- #
    required_design_cols <- c(sample_id_col, group_id_col)
    if (!all(required_design_cols %in% colnames(design_matrix))) {
        warning("Design matrix missing required columns: ", paste(setdiff(required_design_cols, colnames(design_matrix)), collapse=", "))
        return(data.frame(metabolite_id = character(), group = character(), cv = numeric()))
    }
    if (!metabolite_id_col %in% colnames(assay_data)) {
        warning("Metabolite ID column '", metabolite_id_col, "' not found in assay data.")
        return(data.frame(metabolite_id = character(), group = character(), cv = numeric()))
    }

    # --- Data Preparation --- #
    quant_info <- getMetaboliteQuantData(assay_data)
    quant_data <- quant_info$quant_data
    sample_names <- quant_info$sample_names
    annotation_data <- quant_info$annotation_data

    if (length(sample_names) == 0 || nrow(quant_data) == 0) {
        warning("No sample columns or no data rows found for CV calculation")
        return(data.frame(metabolite_id = character(), group = character(), cv = numeric()))
    }

    # Ensure metabolite IDs are present for joining
    quant_data[[metabolite_id_col]] <- annotation_data[[metabolite_id_col]]

    # Select relevant columns from design matrix and ensure consistent types
    design_subset <- design_matrix |>
        dplyr::select(dplyr::all_of(required_design_cols)) |>
        dplyr::rename(Run = {{ sample_id_col }}, group = {{ group_id_col }}) |>
        dplyr::mutate(Run = as.character(Run)) # Convert Run to character to match pivoted data

    # --- Debug output to check groups --- #
    message("Groups in design matrix: ", paste(unique(design_subset$group), collapse = ", "))
    message("Number of unique samples per group:")
    for (g in unique(design_subset$group)) {
        n_samples <- sum(design_subset$group == g)
        message("  - ", g, ": ", n_samples, " samples")
    }

    # --- CV Calculation --- #
    # First pivot the data to long format
    long_data <- quant_data |>
        tidyr::pivot_longer(
            cols = all_of(sample_names),
            names_to = "Run",
            values_to = "Intensity"
        ) |>
        # Ensure Run column is character for joining
        dplyr::mutate(Run = as.character(.data$Run))
    
        # Join with design info
    long_data_with_groups <- dplyr::left_join(long_data, design_subset, by = "Run")
    
    # Check if joining worked correctly
    if (sum(is.na(long_data_with_groups$group)) > 0) {
        message("Warning: Some samples couldn't be matched to groups. Check sample IDs in design matrix.")
        message("Unmatched samples: ", paste(unique(long_data_with_groups$Run[is.na(long_data_with_groups$group)]), collapse = ", "))
    }
    
    # Remove rows where joining failed or Intensity is NA or 0
    filtered_data <- stats::na.omit(long_data_with_groups) |>
        dplyr::filter(.data$Intensity > 0) # Filter out zero values which can inflate CV
    
        # Group by metabolite and experimental group
    cv_data <- filtered_data |>
        dplyr::group_by(dplyr::across(dplyr::all_of(c(metabolite_id_col, "group")))) |>
        # Calculate SD and Mean, requiring at least 2 data points for SD
        dplyr::summarise(
            n_samples = dplyr::n(),
            mean_intensity = mean(.data$Intensity, na.rm = TRUE),
            sd_intensity = if(dplyr::n() >= 2) stats::sd(.data$Intensity, na.rm = TRUE) else NA_real_,
            .groups = "drop"
        ) |>
        # Calculate CV, handle mean close to zero or NA sd
        dplyr::mutate(
            cv = dplyr::case_when(
                n_samples < 2 ~ NA_real_, # CV is NA if fewer than 2 samples in group
                is.na(.data$sd_intensity) ~ NA_real_, # CV is NA if sd is NA
                abs(.data$mean_intensity) < .Machine$double.eps ~ NA_real_, # Avoid division by tiny number
                TRUE ~ (.data$sd_intensity / .data$mean_intensity) * 100
            ),
            # Cap CV at a reasonable maximum to avoid extreme outliers
            cv = pmin(cv, 200) # Cap at 200% to avoid extreme outliers affecting visualization
        ) |>
        # Keep only relevant columns and rename metabolite ID column back
        dplyr::select(dplyr::all_of(c(metabolite_id_col, "group", "cv", "n_samples"))) |>
        dplyr::rename(metabolite_id = {{ metabolite_id_col }})
    
    # Summary statistics for debugging
    message("CV calculation complete")
    message("CV summary statistics per group:")
    for (g in unique(cv_data$group)) {
        group_cvs <- cv_data$cv[cv_data$group == g & !is.na(cv_data$cv)]
        if (length(group_cvs) > 0) {
            group_stats <- summary(group_cvs)
            message("  - Group ", g, ":")
            message("    Min: ", round(group_stats[1], 1), "%, 1st Qu: ", round(group_stats[2], 1), 
                   "%, Median: ", round(group_stats[3], 1), "%, Mean: ", round(group_stats[4], 1), 
                   "%, 3rd Qu: ", round(group_stats[5], 1), "%, Max: ", round(group_stats[6], 1), "%")
            message("    Number of metabolites with CV: ", length(group_cvs))
        } else {
            message("  - Group ", g, ": No valid CV values")
        }
    }

    return(cv_data)
}

#' @title Calculate Internal Standard Metrics
#' @description Identifies internal standards (IS) based on a regex pattern and
#'              calculates their mean intensity and CV across all samples.
#'
#' @param assay_data A tibble/data.frame for one assay.
#' @param is_pattern Character string, regex pattern to identify IS in the metabolite_id_col.
#'                  If NULL, NA, or empty, the function returns an empty data frame.
#' @param metabolite_id_col Character string, the name of the metabolite ID column.
#' @param sample_id_col ***Currently unused. Sample columns inferred by numeric type.***
#'
#' @return A data frame with columns 'is_id', 'mean_intensity', 'cv'.
#'         Returns an empty data frame if no IS are found or pattern is invalid.
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr filter group_by summarise mutate select rename all_of
#' @importFrom stringr str_detect
#' @importFrom stats sd na.omit
#' @keywords internal
#' @noRd   
#' @export
getInternalStandardMetrics <- function(assay_data,
                                       is_pattern,
                                       metabolite_id_col,
                                       sample_id_col) {

    # --- Input Validation --- #
    if (is.null(is_pattern) || is.na(is_pattern) || is_pattern == "") {
        # message("Internal standard pattern is missing or invalid. Skipping IS metrics.")
        return(data.frame(is_id = character(), mean_intensity = numeric(), cv = numeric()))
    }
    if (!metabolite_id_col %in% colnames(assay_data)) {
        warning("Metabolite ID column '", metabolite_id_col, "' not found in assay data for IS metrics.")
        return(data.frame(is_id = character(), mean_intensity = numeric(), cv = numeric()))
    }

    # --- Data Preparation --- #
    quant_info <- getMetaboliteQuantData(assay_data)
    quant_data <- quant_info$quant_data
    sample_names <- quant_info$sample_names
    annotation_data <- quant_info$annotation_data

    if (length(sample_names) == 0 || nrow(quant_data) == 0) {
        return(data.frame(is_id = character(), mean_intensity = numeric(), cv = numeric()))
    }

    # Ensure metabolite IDs are present
    quant_data[[metabolite_id_col]] <- annotation_data[[metabolite_id_col]]

    # --- Identify and Filter IS --- #
    is_rows <- annotation_data |>
        dplyr::filter(stringr::str_detect(.data[[metabolite_id_col]], {{ is_pattern }})) |>
        dplyr::pull({{ metabolite_id_col }})

    if (length(is_rows) == 0) {
        # message("No internal standards found matching pattern: ", is_pattern)
        return(data.frame(is_id = character(), mean_intensity = numeric(), cv = numeric()))
    }

    is_data <- quant_data |>
        dplyr::filter(.data[[metabolite_id_col]] %in% is_rows)

    # --- Calculate Metrics --- #
    is_metrics <- is_data |>
        tidyr::pivot_longer(
            cols = dplyr::all_of(sample_names),
            names_to = "Run",
            values_to = "Intensity"
        ) |>
        stats::na.omit() |> # Remove missing values before calculating mean/sd
        dplyr::group_by(dplyr::across(dplyr::all_of(metabolite_id_col))) |>
        dplyr::summarise(
            mean_intensity = mean(.data$Intensity, na.rm = TRUE),
            sd_intensity = if(dplyr::n() >= 2) stats::sd(.data$Intensity, na.rm = TRUE) else NA_real_,
            n_samples = dplyr::n(), # Keep track of how many samples contributed
            .groups = "drop"
        ) |>
        dplyr::mutate(
            cv = dplyr::case_when(
                 is.na(.data$sd_intensity) ~ NA_real_,
                 abs(.data$mean_intensity) < .Machine$double.eps ~ NA_real_,
                 TRUE ~ (.data$sd_intensity / .data$mean_intensity) * 100
            )
        ) |>
        dplyr::select(dplyr::all_of(c(metabolite_id_col, "mean_intensity", "cv"))) |>
        dplyr::rename(is_id = {{ metabolite_id_col }})

    return(is_metrics)
}

#' @title Update and Visualize Metabolomics Filtering Progress
#' @description Tracks and visualizes the impact of filtering steps on metabolomics
#'              data, updating a global `FilteringProgressMetabolomics` object.
#'              Generates QC plots summarizing the changes.
#'
#' @details
#' This function serves as the primary interface for tracking metabolomics QC.
#' It performs the following:
#' \itemize{
#'   \item Initializes or retrieves the global `FilteringProgressMetabolomics` object.
#'   \item Takes a `MetaboliteAssayData` object (or similar S4 containing assays list
#'         and design matrix) as input for the current processing step.
#'   \item Extracts the list of assay data frames/tibbles and the design matrix.
#'   \item For each assay, calls helper functions (`countUniqueMetabolites`,
#'         `countMetabolitesPerSample`, `calculateMissingness`,
#'         `calculateSumIntensityPerSample`, `calculateMetaboliteCVs`,
#'         `getInternalStandardMetrics`) to calculate QC metrics.
#'   \item Calls `calculateTotalUniqueMetabolitesAcrossAssays` for an overall count.
#'   \item Updates the global `FilteringProgressMetabolomics` object with the
#'         calculated metrics for the given `step_name`.
#'   \item Generates summary plots visualizing the tracked metrics across steps.
#'   \item Optionally saves plots to disk if `publication_graphs_dir` and a time
#'         directory are available.
#'   \item Returns either a combined grid plot or an invisible list of plots.
#' }
#'
#' **Important:** Relies on and modifies the global `filtering_progress_metabolomics` object.
#' Requires helper functions (defined previously in this file) to be available.
#' For plot saving, requires either `time_dir` in the global environment or access to
#' the appropriate directory via `project_dirs$omics_type$time_dir`.
#'
#' @param theObject A S4 object (e.g., `MetaboliteAssayData`, `SummarizedExperiment`,
#'                  `MultiAssayExperiment`) containing metabolomics data. Must provide
#'                  access to a list of assays (data frames/tibbles with metabolite rows
#'                  and sample columns) and a colData/design matrix linking samples to groups.
#' @param step_name Character string uniquely identifying the current filtering step.
#' @param publication_graphs_dir Optional path for saving plots. If provided, the function
#'                              will try to find the corresponding time_dir.
#' @param omics_type Optional character string specifying the omics type (e.g., "metabolomics").
#'                  If provided and project_dirs exists in the global environment, will use
#'                  project_dirs[[omics_type]]$time_dir for plot saving.
#' @param time_dir Optional explicit path to the time directory. If provided, this overrides
#'                other methods of finding the time directory.
#' @param overwrite Logical, whether to overwrite existing data for `step_name`.
#' @param return_grid Logical, whether to return a `gridExtra` combined plot.
#' @param group_id_col Character, name of the column in `colData(theObject)` specifying groups.
#' @param sample_id_col Character, name of the sample ID column in `colData(theObject)`.
#' @param metabolite_id_col Character, name of the metabolite ID column in assay data.
#' @param is_pattern Character, regex for identifying internal standards. If not provided,
#'                   attempts to get from `theObject@internal_standard_regex` if slot exists.
#'
#' @return If `return_grid` is `TRUE`, a `grob` object. Otherwise, an invisible list
#'         containing individual `ggplot` objects.
#'
#' @importFrom methods slotNames is
#' @importFrom SummarizedExperiment assays colData assayNames
#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 ggsave
#' @importFrom purrr map map_dfr imap imap_dfr walk iwalk
#' @export
updateMetaboliteFiltering <- function(theObject,
                                      step_name,
                                      publication_graphs_dir = NULL,
                                      omics_type = NULL,
                                      time_dir = NULL,
                                      overwrite = FALSE,
                                      return_grid = FALSE,
                                      group_id_col = NULL,
                                      sample_id_col = NULL,
                                      metabolite_id_col = NULL,
                                      is_pattern = NULL) {


    prog_met <- getFilteringProgressMetabolomics()


    if (!methods::is(theObject, "S4")) {
        stop("`theObject` must be an S4 object.")
    }
    
    # Check for specific class before checking generic slots
    if (inherits(theObject, "MetaboliteAssayData")) {
        # Specific handling for MetaboliteAssayData
        design_matrix <- theObject@design_matrix
        assay_list    <- theObject@metabolite_data # Directly access the slot
        assay_names   <- names(assay_list)
        if(is.null(assay_names)) assay_names <- paste0("Assay_", seq_along(assay_list))
        names(assay_list) <- assay_names
    } else {
        # Basic checks for required methods/slots for non-MetaboliteAssayData objects
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
       stop("Package 'SummarizedExperiment' needed for this function to work.")
    }
    if (!("assays" %in% methods::slotNames(theObject) || canCoerce(theObject, "SummarizedExperiment"))) {
        stop("`theObject` must have an accessible `assays` method or slot.")
    }
     if (!("colData" %in% methods::slotNames(theObject) || canCoerce(theObject, "SummarizedExperiment"))) {
        stop("`theObject` must have an accessible `colData` method or slot.")
        }
        
        # Use SE accessors if available
        design_matrix <- SummarizedExperiment::colData(theObject)
        assay_list <- SummarizedExperiment::assays(theObject)
        # Ensure assay_list elements are data frames/tibbles if needed by helpers
        assay_list <- lapply(assay_list, as.data.frame)
        if(is.null(names(assay_list))) assay_names <- paste0("Assay_", seq_along(assay_list)) else assay_names <- names(assay_list)
        names(assay_list) <- assay_names
    }

    # Attempt to get parameters from object slots if not provided
    if (is.null(group_id_col) && "group_id" %in% slotNames(theObject)) group_id_col <- theObject@group_id
    if (is.null(sample_id_col) && "sample_id" %in% slotNames(theObject)) sample_id_col <- theObject@sample_id
    if (is.null(metabolite_id_col) && "metabolite_id_column" %in% slotNames(theObject)) metabolite_id_col <- theObject@metabolite_id_column
    if (is.null(is_pattern) && "internal_standard_regex" %in% slotNames(theObject)) is_pattern <- theObject@internal_standard_regex

    # Check if essential parameters are now available
    if (is.null(group_id_col)) stop("group_id_col must be provided or accessible via theObject@group_id")
    if (is.null(sample_id_col)) stop("sample_id_col must be provided or accessible via theObject@sample_id")
    if (is.null(metabolite_id_col)) stop("metabolite_id_col must be provided or accessible via theObject@metabolite_id_column")

    # Convert design matrix rownames to column if needed
    if (!sample_id_col %in% colnames(design_matrix)) {
        # If not, check if the rownames seem to match the sample IDs
        if(identical(rownames(design_matrix), as.character(design_matrix[[sample_id_col]]))){
             # This case is unlikely if sample_id_col isn't a column name
             warning("Sample ID column '", sample_id_col, "' not found, but rownames might match? Check object structure.")
        } else if (!is.null(rownames(design_matrix)) && sample_id_col == "Run") { # Heuristic: If rownames exist and user expects 'Run'
            message("Moving rownames of design matrix to column: ", sample_id_col)
            design_matrix <- as.data.frame(design_matrix)
            design_matrix[[sample_id_col]] <- rownames(design_matrix)
            rownames(design_matrix) <- NULL # Remove rownames after moving
        } else {
             warning("Sample ID column '", sample_id_col, "' not found in design matrix and rownames don't seem to match or weren't checked.")
        }
    }
    
    design_matrix <- as.data.frame(design_matrix) # Ensure it's a data frame


    metrics_list_this_step <- list()
    if(length(assay_list) > 0) {
        metrics_list_this_step <- purrr::imap(assay_list, function(current_assay_data, current_assay_name) {
            # Ensure current_assay_data is a data frame/tibble
            if (!is.data.frame(current_assay_data)) {
                warning("Assay ", current_assay_name, " is not a data frame. Skipping metrics calculation.")
                # Return placeholder for non-data frames
                return(list(
                    n_metabolites = 0,
                    detected_per_sample = data.frame(Run = character(), n_detected = integer()),
                    missingness = NA_real_,
                    sum_intensity_per_sample = data.frame(Run = character(), sum_intensity = numeric()),
                    cv_distribution = data.frame(metabolite_id = character(), group = character(), cv = numeric()),
                    is_metrics = data.frame(is_id = character(), mean_intensity = numeric(), cv = numeric())
                ))
            }

            list(
                n_metabolites = countUniqueMetabolites(current_assay_data, metabolite_id_col),
                detected_per_sample = countMetabolitesPerSample(current_assay_data, sample_id_col, metabolite_id_col),
                missingness = calculateMissingness(current_assay_data, sample_id_col),
                sum_intensity_per_sample = calculateSumIntensityPerSample(current_assay_data, sample_id_col),
                cv_distribution = calculateMetaboliteCVs(current_assay_data, design_matrix, group_id_col, NULL, sample_id_col, metabolite_id_col), # replicate_id_col unused for now
                is_metrics = getInternalStandardMetrics(current_assay_data, is_pattern, metabolite_id_col, sample_id_col)
            )
        })
    } else {
        warning("No assays found in theObject.")
        # Handle empty assay list case - maybe stop or proceed with empty metrics
        return(invisible(NULL)) # Or handle appropriately
    }


    total_metabolites <- calculateTotalUniqueMetabolitesAcrossAssays(assay_list, metabolite_id_col)


    updateFilteringProgressMetabolomics(prog_met, step_name, assay_names, metrics_list_this_step, total_metabolites, overwrite)


    plot_list <- generateMetaboliteFilteringPlots(getFilteringProgressMetabolomics()) # Pass updated global object

    # --- 9. Directory handling and plot saving --- #
    actual_save_dir <- NULL # Will hold the final directory path for saving

    message("--- Plot Saving Diagnostics ---")
    message(sprintf("Value of publication_graphs_dir argument in function call: %s", ifelse(is.null(publication_graphs_dir), "NULL", publication_graphs_dir)))
    # The omics_type and time_dir arguments to this function are not directly used for path construction here;
    # we rely on global omic_type, experiment_label, and project_dirs.
    message(sprintf("Value of omics_type argument in function call: %s", ifelse(is.null(omics_type), "NULL", omics_type)))
    message(sprintf("Value of time_dir argument in function call: %s", ifelse(is.null(time_dir), "NULL", time_dir)))

    message("Attempting to determine save directory using global project_dirs, omic_type, and experiment_label...")
    
    if (exists("project_dirs", envir = .GlobalEnv) &&
        exists("omic_type", envir = .GlobalEnv) &&
        exists("experiment_label", envir = .GlobalEnv)) {
        
        message("Global variables project_dirs, omic_type, and experiment_label found.")
        local_project_dirs <- get("project_dirs", envir = .GlobalEnv)
        local_omic_type <- get("omic_type", envir = .GlobalEnv)
        local_experiment_label <- get("experiment_label", envir = .GlobalEnv)
        message(sprintf("Global omic_type value: '%s', Global experiment_label value: '%s'", local_omic_type, local_experiment_label))
        
        omics_key <- paste0(local_omic_type, "_", local_experiment_label)
        message(sprintf("Constructed omics_key for project_dirs: '%s'", omics_key))
        
        if (omics_key %in% names(local_project_dirs) &&
            !is.null(local_project_dirs[[omics_key]]) && 
            "time_dir" %in% names(local_project_dirs[[omics_key]])) {
            
            message(sprintf("omics_key '%s' found in project_dirs and has a 'time_dir' entry.", omics_key))
            retrieved_time_dir <- local_project_dirs[[omics_key]]$time_dir
            message(sprintf("Retrieved time_dir from project_dirs: %s", ifelse(is.null(retrieved_time_dir), "NULL", retrieved_time_dir)))
            
            if (is.null(retrieved_time_dir) || !is.character(retrieved_time_dir) || !nzchar(retrieved_time_dir)) {
                warning(sprintf("project_dirs[['%s']]$time_dir is NULL, not a character string, or empty. Plots will not be saved.", omics_key))
                message("Reason: retrieved_time_dir is invalid.")
            } else {
                actual_save_dir <- retrieved_time_dir
                message(sprintf("Successfully set actual_save_dir for plot saving to: '%s'", actual_save_dir))
            }
        } else {
            warning(sprintf("Could not find omics_key '%s' in project_dirs, or it lacks a 'time_dir' entry. Plots will not be saved.", omics_key))
            message(sprintf("Details: omics_key '%s' in names(project_dirs): %s. project_dirs[['%s']] is NULL: %s. 'time_dir' in names(project_dirs[['%s']]): %s.", 
                          omics_key, omics_key %in% names(local_project_dirs), 
                          omics_key, is.null(local_project_dirs[[omics_key]]), 
                          omics_key, if (omics_key %in% names(local_project_dirs) && !is.null(local_project_dirs[[omics_key]])) "time_dir" %in% names(local_project_dirs[[omics_key]]) else NA))
            message("Available keys in project_dirs: ", paste(names(local_project_dirs), collapse=", "))
        }
    } else {
        warning("One or more global variables ('project_dirs', 'omic_type', 'experiment_label') not found. Plots will not be saved.")
        message(sprintf("Exists project_dirs: %s, Exists omic_type: %s, Exists experiment_label: %s", 
                      exists("project_dirs", envir = .GlobalEnv),
                      exists("omic_type", envir = .GlobalEnv),
                      exists("experiment_label", envir = .GlobalEnv)))
    }

    # Proceed with saving ONLY if actual_save_dir was successfully determined from global project_dirs
    if (!is.null(actual_save_dir)) {
        message(sprintf("Proceeding to save plots to derived directory: %s", actual_save_dir))
        if (!dir.exists(actual_save_dir)) {
            dir.create(actual_save_dir, recursive = TRUE)
            message("Created directory for QC plots: ", actual_save_dir)
        }

        # Save individual plots directly into actual_save_dir (which is the time_dir)
        purrr::iwalk(plot_list, function(plot, plot_name) {
            filename <- file.path(actual_save_dir, sprintf("%s_%s.png", step_name, plot_name))
            message(sprintf("Saving plot: %s", filename))
            ggsave(filename, 
                   plot = plot, 
                   width = 10, 
                   height = 8, 
                   dpi = 300)
        })

        # Save combined grid if return_grid is TRUE and plots exist
        if (return_grid && length(plot_list) > 0 && !is.null(plot_list[[1]]) && inherits(plot_list[[1]], "ggplot")) {
             grid_plot_obj <- do.call(gridExtra::grid.arrange, c(plot_list, ncol = 2)) 
             filename_grid <- file.path(actual_save_dir, sprintf("%s_combined_plots.png", step_name))
             message(sprintf("Saving combined grid plot: %s", filename_grid))
             ggsave(filename_grid, plot = grid_plot_obj, width = 15, height = 15, dpi = 300)
        }
        message("Metabolomics QC plots saved to: ", actual_save_dir)
    } else {
        # This block means actual_save_dir is still NULL. 
        # This implies either globals were missing or project_dirs structure was invalid for deriving time_dir.
        message("No valid save directory determined from global project_dirs. Plots will not be saved.")
        # If publication_graphs_dir was provided in the call, and we still ended up here, it means the global lookup failed.
        if(!is.null(publication_graphs_dir)){
            warning("Function was called with a publication_graphs_dir path, but plot saving still failed because a valid time_dir could not be derived from global project_dirs.")
        }
    }
    message("--- End of Plot Saving Diagnostics ---")

    if (return_grid) {
        # Combine plots into a grid
        if (length(plot_list) > 0 && !is.null(plot_list[[1]]) && inherits(plot_list[[1]], "ggplot")) {
            message("Attempting to arrange plots into a grid.")
            grid_plot_obj <- do.call(gridExtra::grid.arrange, c(plot_list, ncol = 2))
            return(grid_plot_obj)
        } else {
            message("No valid plots to arrange into a grid, or plot_list is empty. Returning NULL for the grid.")
            return(NULL)
        }
    } else {
        # Print each plot individually
        if (length(plot_list) > 0) {
            message("Printing plots individually as return_grid is FALSE or grid could not be formed.")
            purrr::walk(plot_list, function(plot) {
                if (inherits(plot, "ggplot")) {
                    print(plot)
                } else {
                    message("Encountered a non-ggplot object in plot_list when trying to print individually.")
                }
            })
        } else {
            message("Plot list is empty, nothing to print.")
        }
        # Return the list invisibly
        invisible(plot_list)
    }
}

#' @title Generate Metabolomics Filtering Progress Plots
#' @description Creates a set of quality control plots based on the metrics
#'              stored in the global `FilteringProgressMetabolomics` object.
#'
#' @details
#' Generates visualizations for key metabolomics QC metrics across processing steps:
#' \itemize{
#'   \item Total unique metabolites across assays per step (bar chart)
#'   \item Metabolites per assay per step (bar chart)
#'   \item Detected metabolites per sample per step (line chart)
#'   \item Missing value percentage per assay per step (bar chart)
#'   \item Total intensity per sample per step (line chart)
#'   \item Coefficient of variation distribution per step (box plot)
#'   \item Internal standard metrics (CV and intensity) per step (if available)
#' }
#'
#' @param prog_met A `FilteringProgressMetabolomics` object containing the tracked
#'                 metrics across various processing steps. If not provided, the
#'                 function retrieves the global object.
#'
#' @return A list of ggplot objects, one for each visualization type.
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_text labs theme_minimal theme element_text element_blank geom_line geom_point scale_color_manual annotate theme_void geom_boxplot coord_cartesian facet_wrap geom_violin scale_fill_brewer position_dodge
#' @importFrom dplyr bind_rows mutate group_by ungroup
#' @importFrom forcats fct_reorder
#' @importFrom tidyr pivot_longer
#' @keywords internal
#' @noRd
#' @export
generateMetaboliteFilteringPlots <- function(prog_met = NULL) {
    # Get the global object if not provided
    if (is.null(prog_met)) {
        prog_met <- getFilteringProgressMetabolomics()
    }
    
    # Return empty list if no steps have been tracked
    if (length(prog_met@steps) == 0) {
        message("No metabolomics filtering steps have been tracked yet.")
        return(list())
    }
    
    plot_list <- list()
    
    # --- 1. Total Metabolites Plot (All Assays Combined) --- #
    plot_list$total_metabolites <- ggplot(data.frame(
        step = factor(prog_met@steps, levels = prog_met@steps),
        total_metabolites = prog_met@n_metabolites_total
    ), aes(x = step, y = total_metabolites)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
        geom_text(aes(label = total_metabolites), 
                  vjust = -0.5, 
                  size = 4) +
        labs(
            title = "Total Unique Metabolites (All Assays)",
            x = "Filtering Step",
            y = "Unique Metabolites"
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major.x = element_blank()
        )
    
    # --- 2. Metabolites Per Assay Plot --- #
    # First create a data frame with metabolites per assay per step - use purrr
    metabolites_per_assay_df <- purrr::map_dfr(seq_along(prog_met@steps), function(step_idx) {
        step <- prog_met@steps[step_idx]
        assay_names <- prog_met@assay_names[[step_idx]]
        metabolite_counts <- unlist(prog_met@n_metabolites_per_assay[[step_idx]])
        
        if (length(metabolite_counts) > 0) {
            return(data.frame(
                step = step,
                assay = names(metabolite_counts),
                n_metabolites = as.numeric(metabolite_counts)
            ))
        }
        return(NULL)  # Return NULL if no metabolite counts (will be filtered by map_dfr)
    })
    
    if (nrow(metabolites_per_assay_df) > 0) {
        plot_list$metabolites_per_assay <- ggplot(metabolites_per_assay_df, 
                                                 aes(x = factor(step, levels = prog_met@steps), 
                                                     y = n_metabolites,
                                                     fill = assay)) +
            geom_bar(stat = "identity", position = "dodge", width = 0.7) +
            geom_text(aes(label = n_metabolites), 
                      position = position_dodge(width = 0.7),
                      vjust = -0.5, 
                      size = 3) +
            labs(
                title = "Metabolites per Assay",
                x = "Filtering Step",
                y = "Unique Metabolites",
                fill = "Assay"
            ) +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.major.x = element_blank()
            ) +
            scale_fill_brewer(palette = "Set1")
    }
    
    # --- 3. Detected Metabolites Per Sample Plot --- #
    # Create a data frame with detected metabolites per sample per step - use purrr
    detected_per_sample_df <- purrr::map_dfr(seq_along(prog_met@steps), function(step_idx) {
        step <- prog_met@steps[step_idx]
        assay_names <- prog_met@assay_names[[step_idx]]
        
        # Use purrr::map_dfr to combine results from all assays
        purrr::map_dfr(assay_names, function(assay_name) {
            if (assay_name %in% names(prog_met@detected_per_sample[[step_idx]])) {
                detected_df <- prog_met@detected_per_sample[[step_idx]][[assay_name]]
                if (nrow(detected_df) > 0) {
                    detected_df$step <- step
                    detected_df$assay <- assay_name
                    return(detected_df)
                }
            }
            return(NULL)  # Return NULL if no data (will be filtered by map_dfr)
        })
    })
    
    if (nrow(detected_per_sample_df) > 0) {
        # Ensure consistent data types
        detected_per_sample_df$Run <- as.character(detected_per_sample_df$Run)
        detected_per_sample_df$n_detected <- as.numeric(detected_per_sample_df$n_detected)
        detected_per_sample_df$step <- factor(detected_per_sample_df$step, 
                                             levels = prog_met@steps)
        
        plot_list$detected_per_sample <- detected_per_sample_df |>
            group_by(Run, assay) |>
            mutate(avg_detected = mean(n_detected)) |>
            ungroup() |>
            mutate(Run = fct_reorder(Run, avg_detected)) |>
            ggplot(aes(x = Run, y = n_detected, 
                       group = interaction(step, assay), 
                       color = step,
                       linetype = assay)) +
            geom_line() +
            geom_point() +
            labs(
                title = "Detected Metabolites per Sample",
                x = "Sample ID (ordered by average detected metabolites)",
                y = "Number of Detected Metabolites",
                color = "Step",
                linetype = "Assay"
            ) +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.major.x = element_blank()
            ) +
            scale_color_brewer(palette = "Set2")
    }
    
    # --- 4. Missingness Percentage Per Assay Plot --- #
    # Create a data frame with missingness percentage per assay per step - use purrr
    missingness_df <- purrr::map_dfr(seq_along(prog_met@steps), function(step_idx) {
        step <- prog_met@steps[step_idx]
        assay_names <- prog_met@assay_names[[step_idx]]
        
        # Use purrr::map_dfr to combine results from all assays
        purrr::map_dfr(assay_names, function(assay_name) {
            if (assay_name %in% names(prog_met@missingness_per_assay[[step_idx]])) {
                missingness <- prog_met@missingness_per_assay[[step_idx]][[assay_name]]
                if (!is.null(missingness) && !is.na(missingness)) {
                    return(data.frame(
                        step = step,
                        assay = assay_name,
                        missingness = as.numeric(missingness)
                    ))
                }
            }
            return(NULL)  # Return NULL if no data (will be filtered by map_dfr)
        })
    })
    
    if (nrow(missingness_df) > 0) {
        plot_list$missingness <- ggplot(missingness_df, 
                                       aes(x = factor(step, levels = prog_met@steps), 
                                           y = missingness,
                                           fill = assay)) +
            geom_bar(stat = "identity", position = "dodge", width = 0.7) +
            geom_text(aes(label = sprintf("%.1f%%", missingness)), 
                      position = position_dodge(width = 0.7),
                      vjust = -0.5, 
                      size = 3) +
            labs(
                title = "Missing Values Percentage per Assay",
                x = "Filtering Step",
                y = "Missing Values (%)",
                fill = "Assay"
            ) +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.major.x = element_blank()
            ) +
            scale_fill_brewer(palette = "Set1")
    }
    
    # --- 5. Sum Intensity Per Sample Plot --- #
    # Create a data frame with sum intensity per sample per step - use purrr
    sum_intensity_df <- purrr::map_dfr(seq_along(prog_met@steps), function(step_idx) {
        step <- prog_met@steps[step_idx]
        assay_names <- prog_met@assay_names[[step_idx]]
        
        # Use purrr::map_dfr to combine results from all assays
        purrr::map_dfr(assay_names, function(assay_name) {
            if (assay_name %in% names(prog_met@sum_intensity_per_sample[[step_idx]])) {
                intensity_df <- prog_met@sum_intensity_per_sample[[step_idx]][[assay_name]]
                if (nrow(intensity_df) > 0) {
                    intensity_df$step <- step
                    intensity_df$assay <- assay_name
                    return(intensity_df)
                }
            }
            return(NULL)  # Return NULL if no data (will be filtered by map_dfr)
        })
    })
    
    if (nrow(sum_intensity_df) > 0) {
        # Ensure consistent data types
        sum_intensity_df$Run <- as.character(sum_intensity_df$Run)
        sum_intensity_df$sum_intensity <- as.numeric(sum_intensity_df$sum_intensity)
        sum_intensity_df$step <- factor(sum_intensity_df$step, 
                                       levels = prog_met@steps)
        
        plot_list$sum_intensity <- sum_intensity_df |>
            group_by(Run, assay) |>
            mutate(avg_intensity = mean(sum_intensity)) |>
            ungroup() |>
            mutate(
                Run = fct_reorder(Run, avg_intensity),
                # Apply log2 transformation directly to the values
                log2_intensity = log2(pmax(sum_intensity, 1)) # Avoid log(0) with pmax
            ) |>
            ggplot(aes(x = Run, y = log2_intensity, # Plot the transformed values
                      group = interaction(step, assay), 
                      color = step,
                      linetype = assay)) +
            geom_line() +
            geom_point() +
            labs(
                title = "Total Signal Intensity per Sample",
                x = "Sample ID (ordered by average intensity)",
                y = "Sum Intensity (log2 scale)",
                color = "Step",
                linetype = "Assay"
            ) +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.major.x = element_blank()
            ) +
            scale_color_brewer(palette = "Set2")
    }
    
    # --- 6. CV Distribution Plot --- #
    # Create a data frame with CV distribution per step - use purrr
    cv_df <- purrr::map_dfr(seq_along(prog_met@steps), function(step_idx) {
        step <- prog_met@steps[step_idx]
        assay_names <- prog_met@assay_names[[step_idx]]
        
        # Use purrr::map_dfr to combine results from all assays
        purrr::map_dfr(assay_names, function(assay_name) {
            if (assay_name %in% names(prog_met@cv_distribution_per_assay[[step_idx]])) {
                step_cv_df <- prog_met@cv_distribution_per_assay[[step_idx]][[assay_name]]
                if (nrow(step_cv_df) > 0) {
                    step_cv_df$step <- step
                    step_cv_df$assay <- assay_name
                    return(step_cv_df)
                }
            }
            return(NULL)  # Return NULL if no data (will be filtered by map_dfr)
        })
    })
    
    if (nrow(cv_df) > 0) {
        # Ensure consistent data types and remove outliers
        cv_df$cv <- as.numeric(cv_df$cv)
        cv_df$step <- factor(cv_df$step, levels = prog_met@steps)
        
        # Calculate 95th percentile for y-axis limit
        q95 <- quantile(cv_df$cv, 0.95, na.rm = TRUE)
        
        plot_list$cv_distribution <- ggplot(cv_df, 
                                           aes(x = step, 
                                               y = cv,
                                               fill = assay)) +
            geom_boxplot(alpha = 0.7, outlier.shape = NA) +
            labs(
                title = "CV Distribution by Group",
                x = "Filtering Step",
                y = "Coefficient of Variation (%)",
                fill = "Assay"
            ) +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.major.x = element_blank()
            ) +
            coord_cartesian(ylim = c(0, q95)) +
            scale_fill_brewer(palette = "Set1") +
            facet_wrap(~group, scales = "free_y")
    }
    
    # --- 7. Internal Standards Metrics Plot --- #
    # Create a data frame with internal standard metrics per step - use purrr
    is_df <- purrr::map_dfr(seq_along(prog_met@steps), function(step_idx) {
        step <- prog_met@steps[step_idx]
        assay_names <- prog_met@assay_names[[step_idx]]
        
        # Use purrr::map_dfr to combine results from all assays
        purrr::map_dfr(assay_names, function(assay_name) {
            if (assay_name %in% names(prog_met@is_metrics_per_assay[[step_idx]])) {
                step_is_df <- prog_met@is_metrics_per_assay[[step_idx]][[assay_name]]
                if (nrow(step_is_df) > 0) {
                    step_is_df$step <- step
                    step_is_df$assay <- assay_name
                    return(step_is_df)
                }
            }
            return(NULL)  # Return NULL if no data (will be filtered by map_dfr)
        })
    })
    
    if (nrow(is_df) > 0) {
        # Ensure consistent data types
        is_df$cv <- as.numeric(is_df$cv)
        is_df$mean_intensity <- as.numeric(is_df$mean_intensity)
        is_df$step <- factor(is_df$step, levels = prog_met@steps)
        
        # Create CV plot for internal standards
        plot_list$is_cv <- ggplot(is_df, 
                                 aes(x = step, 
                                     y = cv,
                                     fill = assay)) +
            geom_violin(alpha = 0.7, trim = FALSE, scale = "width") +
            geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +
            labs(
                title = "Internal Standards CV",
                x = "Filtering Step",
                y = "Coefficient of Variation (%)",
                fill = "Assay"
            ) +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.major.x = element_blank()
            ) +
            scale_fill_brewer(palette = "Set1")
        
        # Create mean intensity plot for internal standards
        plot_list$is_intensity <- is_df |>
            # Apply log2 transformation directly to the data
            mutate(log2_intensity = log2(pmax(mean_intensity, 1))) |> # Avoid log(0) with pmax
            ggplot(aes(x = step, 
                       y = log2_intensity, # Plot transformed values
                       fill = assay)) +
            geom_violin(alpha = 0.7, trim = FALSE, scale = "width") +
            geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +
            labs(
                title = "Internal Standards Mean Intensity",
                x = "Filtering Step",
                y = "Mean Intensity (log2 scale)",
                fill = "Assay"
            ) +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.major.x = element_blank()
            ) +
            scale_fill_brewer(palette = "Set1")
    }
    
    return(plot_list)
}
