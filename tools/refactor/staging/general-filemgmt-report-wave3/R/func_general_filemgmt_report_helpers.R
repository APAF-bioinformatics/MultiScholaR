# ----------------------------------------------------------------------------
# savePlot
# ----------------------------------------------------------------------------
#' Save a plot in multiple formats
#'
#' This function saves a given plot in multiple specified formats and also save the ggplot object as a rds file.
#'
#' @param plot The plot object to be saved
#' @param base_path The base directory path where the plot will be saved
#' @param plot_name The name to be used for the saved plot files
#' @param formats A vector of file formats to save the plot in (default: c("pdf", "png"))
#' @param width The width of the plot (default: 7)
#' @param height The height of the plot (default: 7)
#' @param ... Additional arguments to be passed to ggsave
#'
#' @return This function is called for its side effects (saving files)
#' @export
#'
#'
savePlot <- function(plot, base_path, plot_name, formats = c("pdf", "png"), width = 7, height = 7, ...) {
    # Always save the RDS (works for both single plots and lists)
    saveRDS(plot, file.path(base_path, paste0(plot_name, ".rds")))

    # Check if plot is a list of plots
    if (is.list(plot) && !inherits(plot, "gg")) {
        # It's a list of plots - save each one individually
        plot_names <- names(plot)
        if (is.null(plot_names)) {
            plot_names <- paste0("plot_", seq_along(plot))
        }

        purrr::walk2(plot, plot_names, function(p, pname) {
            if (inherits(p, "gg")) {
                purrr::walk(formats, function(format) {
                    file_path <- file.path(base_path, paste0(plot_name, "_", pname, ".", format))
                    
                    # Use cairo_pdf for PDF format to avoid font issues on macOS
                    save_device <- format
                    if (format == "pdf") {
                        save_device <- grDevices::cairo_pdf
                    }
                    
                    ggsave(filename = file_path, plot = p, device = save_device, width = width, height = height, ...)
                })
            }
        })
    } else {
        # Single plot - original behavior
        purrr::walk(formats, \(format){
            file_path <- file.path(base_path, paste0(plot_name, ".", format))
            
            # Use cairo_pdf for PDF format to avoid font issues on macOS
            save_device <- format
            if (format == "pdf") {
                save_device <- grDevices::cairo_pdf
            }
            
            ggsave(filename = file_path, plot = plot, device = save_device, width = width, height = height, ...)
        })
    }
}

# ----------------------------------------------------------------------------
# save_plot
# ----------------------------------------------------------------------------
#' Save a plot in multiple formats
#'
#' This function saves a given plot in multiple specified formats and also save the ggplot object as a rds file.
#'
#' @param plot The plot object to be saved
#' @param base_path The base directory path where the plot will be saved
#' @param plot_name The name to be used for the saved plot files
#' @param formats A vector of file formats to save the plot in (default: c("pdf", "png"))
##' @param width The width of the plot (default: 7)
#' @param height The height of the plot (default: 7)
#' @param ... Additional arguments to be passed to ggsave
#'
#' @return This function is called for its side effects (saving files)
#' @export
#'
#'
save_plot <- function(plot, base_path, plot_name, formats = c("pdf", "png"), width = 7, height = 7, ...) {
    savePlot(plot, base_path, plot_name, formats, width, height, ...)
}

# ----------------------------------------------------------------------------
# write_results
# ----------------------------------------------------------------------------
write_results <- function(data, filename) {
    vroom::vroom_write(data, file.path(results_dir, "protein_qc", filename))
}

