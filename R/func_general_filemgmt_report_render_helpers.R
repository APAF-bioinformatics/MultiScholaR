#' @title Download Report Template from GitHub
#' @description Downloads a report template from the MultiScholaR GitHub repository
#'              and caches it locally for future use.
#' @param omic_type Character string, the omic type (e.g., "proteomics")
#' @param rmd_filename Character string, the template filename (e.g., "DIANN_limpa_report.rmd")
#' @return Character string, path to the downloaded/cached template file
#' @keywords internal
downloadReportTemplate <- function(omic_type, rmd_filename) {
    # Create cache directory using tools (base R, no extra dependencies)
    if (requireNamespace("rappdirs", quietly = TRUE)) {
        cache_base <- rappdirs::user_cache_dir("MultiScholaR")
    } else {
        # Fallback to temp directory if rappdirs not available
        cache_base <- file.path(tempdir(), "MultiScholaR_cache")
    }

    cache_dir <- file.path(cache_base, "report_templates", omic_type, "report")
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

    # Local cache file path
    cached_file <- file.path(cache_dir, rmd_filename)

    # If already cached and less than 7 days old, return it
    if (file.exists(cached_file)) {
        file_age_days <- as.numeric(difftime(Sys.time(), file.info(cached_file)$mtime, units = "days"))
        if (file_age_days < 7) {
            logger::log_info("Using cached template: {cached_file}")
            return(cached_file)
        } else {
            logger::log_info("Cached template is older than 7 days, re-downloading...")
        }
    }

    # GitHub URL (raw content from main branch)
    github_url <- sprintf(
        "https://raw.githubusercontent.com/APAF-bioinformatics/MultiScholaR/main/Workbooks/%s/report/%s",
        omic_type, rmd_filename
    )

    logger::log_info("Downloading template from GitHub: {github_url}")

    # Download
    tryCatch(
        {
            download.file(github_url, cached_file, mode = "wb", quiet = TRUE)
            logger::log_info("Successfully downloaded template to: {cached_file}")
            return(cached_file)
        },
        error = function(e) {
            rlang::abort(paste0(
                "Failed to download template '", rmd_filename, "' from GitHub.\n",
                "URL: ", github_url, "\n",
                "Error: ", e$message
            ))
        }
    )
}

##################################################################################################################
#' @export
#' @importFrom rlang abort
#' @importFrom tools file_path_sans_ext
RenderReport <- function(omic_type,
                         experiment_label,
                         rmd_filename = "DIANN_report.rmd",
                         project_dirs_object_name = "project_dirs",
                         output_format = NULL) {
    message("--- DEBUG66: Entering RenderReport ---")
    message(sprintf("   DEBUG66: omic_type = '%s'", omic_type))
    message(sprintf("   DEBUG66: experiment_label = '%s'", experiment_label))
    message(sprintf("   DEBUG66: rmd_filename = '%s'", rmd_filename))
    message(sprintf("   DEBUG66: project_dirs_object_name = '%s'", project_dirs_object_name))

    # --- Validate Inputs ---
    if (missing(omic_type) || !is.character(omic_type) || length(omic_type) != 1 || omic_type == "") {
        rlang::abort("`omic_type` must be a single non-empty character string.")
    }
    if (missing(experiment_label) || !is.character(experiment_label) || length(experiment_label) != 1 || experiment_label == "") {
        rlang::abort("`experiment_label` must be a single non-empty character string.")
    }
    if (!is.character(rmd_filename) || length(rmd_filename) != 1 || rmd_filename == "") {
        rlang::abort("`rmd_filename` must be a single non-empty character string.")
    }

    message("   DEBUG66: Input validation passed")

    # --- Retrieve Paths from Global Project Directories Object ---
    # Use the new helper function with automatic fallback
    message("   DEBUG66: Calling getProjectPaths...")
    current_paths <- tryCatch(
        {
            getProjectPaths(
                omic_type = omic_type,
                experiment_label = experiment_label,
                project_dirs_object_name = project_dirs_object_name
            )
        },
        error = function(e) {
            rlang::abort(paste0("Failed to get project paths: ", e$message))
        }
    )

    message("   DEBUG66: getProjectPaths completed")
    if (!is.null(current_paths)) {
        message(sprintf("      DEBUG66: current_paths is list: %s", is.list(current_paths)))
        message(sprintf("      DEBUG66: current_paths names: %s", paste(names(current_paths), collapse = ", ")))
        if ("base_dir" %in% names(current_paths)) {
            message(sprintf("      DEBUG66: base_dir = %s", current_paths$base_dir))
        }
        if ("results_summary_dir" %in% names(current_paths)) {
            message(sprintf("      DEBUG66: results_summary_dir = %s", current_paths$results_summary_dir))
        }
        if ("source_dir" %in% names(current_paths)) {
            message(sprintf("      DEBUG66: source_dir = %s", current_paths$source_dir))
        }
    } else {
        message("      DEBUG66: current_paths is NULL!")
    }

    if (!is.list(current_paths) ||
        is.null(current_paths$base_dir) || # Need base_dir to find the template Rmd
        is.null(current_paths$results_summary_dir)) {
        rlang::abort(paste0("Essential paths (base_dir, results_summary_dir) missing from project_dirs"))
        rlang::abort(paste0("Essential paths (base_dir, results_summary_dir) missing from project_dirs"))
    }

    message("   DEBUG66: Path validation passed")

    # --- Determine the source Rmd template directory (e.g., scripts/proteomics) ---
    # This logic mirrors part of setupDirectories to find the correct unlabelled script source leaf.
    omic_script_template_leaf <- switch(omic_type,
        proteomics = "proteomics",
        metabolomics = "metabolomics",
        transcriptomics = "transcriptomics",
        lipidomics = "lipidomics",
        integration = "integration",
        {
            rlang::abort(paste0("Internal error: Unrecognized omic_type ", sQuote(omic_type), " for Rmd template path construction."))
        }
    )

    rmd_template_dir <- file.path(current_paths$base_dir, "scripts", omic_script_template_leaf)
    rmd_input_path <- file.path(rmd_template_dir, rmd_filename)

    message(sprintf("   DEBUG66: rmd_input_path = %s", rmd_input_path))
    message(sprintf("   DEBUG66: Template file exists: %s", file.exists(rmd_input_path)))

    if (!file.exists(rmd_input_path)) {
        rlang::abort(paste0(
            "R Markdown template file not found at the expected location: ", sQuote(rmd_input_path),
            ". This should be in the general scripts/<omic_type> directory (e.g., scripts/proteomics)."
        ))
    }

    # --- Construct Output Path (in the labelled results_summary directory) ---
    message("   DEBUG66: Constructing output file path...")
    output_file_basename <- paste0(
        tools::file_path_sans_ext(rmd_filename),
        "_", omic_type,
        "_", experiment_label
    )

    output_ext <- ".docx" # Default
    if (!is.null(output_format)) {
        if (output_format == "word_document" || grepl("word", output_format, ignore.case = TRUE)) {
            output_ext <- ".docx"
        } else if (output_format == "html_document" || grepl("html", output_format, ignore.case = TRUE)) {
            output_ext <- ".html"
        } else if (grepl("pdf", output_format, ignore.case = TRUE)) {
            output_ext <- ".pdf"
        }
    }

    output_file_path <- file.path(current_paths$results_summary_dir, paste0(output_file_basename, output_ext))

    message(sprintf("   DEBUG66: output_file_path = %s", output_file_path))
    message(sprintf("   DEBUG66: output_ext = %s", output_ext))

    logger::log_info("Attempting to render report:")
    logger::log_info("- Rmd Source (Template): {rmd_input_path}")
    logger::log_info("- Output File: {output_file_path}")
    logger::log_info("- Params: omic_type=\'{omic_type}\', experiment_label=\'{experiment_label}\'")

    # Read study_parameters.txt to extract workflow_name and timestamp
    message("   DEBUG66: Reading study_parameters.txt...")
    params_path <- file.path(current_paths$source_dir, "study_parameters.txt")
    message(sprintf("      DEBUG66: params_path = %s", params_path))
    message(sprintf("      DEBUG66: params file exists: %s", file.exists(params_path)))
    if (file.exists(params_path)) {
        lines <- readLines(params_path)
        workflow_name_line <- grep("Workflow Name:", lines, value = TRUE)
        timestamp_line <- grep("Timestamp:", lines, value = TRUE)
        workflow_name <- if (length(workflow_name_line) > 0) trimws(sub("Workflow Name:", "", workflow_name_line[1])) else "Unknown Workflow"
        timestamp <- if (length(timestamp_line) > 0) trimws(sub("Timestamp:", "", timestamp_line[1])) else format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    } else {
        workflow_name <- "Unknown Workflow"
        timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    }

    message(sprintf("   DEBUG66: workflow_name = '%s'", workflow_name))
    message(sprintf("   DEBUG66: timestamp = '%s'", timestamp))

    # --- Render the Report ---
    message("   DEBUG66: About to call rmarkdown::render...")
    message("      DEBUG66: Render parameters:")
    message(sprintf("         input = %s", rmd_input_path))
    message(sprintf("         output_file = %s", output_file_path))
    message(sprintf("         output_format = %s", ifelse(is.null(output_format), "NULL (using Rmd default)", output_format)))

    rendered_path <- tryCatch(
        {
            rmarkdown::render(
                input = rmd_input_path,
                params = list(
                    omic_type = omic_type,
                    experiment_label = experiment_label,
                    workflow_name = workflow_name,
                    timestamp = timestamp
                ),
                output_file = output_file_path,
                output_format = output_format, # Pass this along; if NULL, Rmd default is used
                envir = new.env(parent = globalenv()) # Render in a clean environment
            )
        },
        error = function(e) {
            message("   DEBUG66: rmarkdown::render FAILED")
            message(sprintf("      DEBUG66: Error message: %s", e$message))
            message("      DEBUG66: Error traceback:")
            print(e)
            logger::log_error("Failed to render R Markdown report: {e$message}")
            logger::log_error("Input path: {rmd_input_path}")
            logger::log_error("Output path: {output_file_path}")
            NULL # Return NULL on failure
        }
    )

    message(sprintf("   DEBUG66: rmarkdown::render returned: %s", ifelse(is.null(rendered_path), "NULL", rendered_path)))
    if (!is.null(rendered_path)) {
        message(sprintf("   DEBUG66: Output file exists: %s", file.exists(rendered_path)))
    }

    if (!is.null(rendered_path) && file.exists(rendered_path)) {
        logger::log_info("Report successfully rendered to: {rendered_path}")
    } else {
        logger::log_warn("Report rendering failed or output file not found at expected location.")
    }

    message("--- DEBUG66: Exiting RenderReport ---")
    message(sprintf("   DEBUG66: Returning: %s", ifelse(is.null(rendered_path), "NULL", rendered_path)))

    invisible(rendered_path)
}

# ----------------------------------------------------------------------------
# saveListOfPdfs
# ----------------------------------------------------------------------------
#' @export
saveListOfPdfs <- function(list, filename) {
    # start pdf
    cairo_pdf(filename)

    # loop
    # purrr::walk( list, print)
    for (p in list) {
        print(p)
    }

    # end pdf
    dev.off()

    invisible(NULL)
}

# ----------------------------------------------------------------------------
# sourceRmdFileSimple
# ----------------------------------------------------------------------------
## Function to source Rmd files
# https://stackoverflow.com/questions/10966109/how-to-source-r-markdown-file-like-sourcemyfile-r
#' @export
sourceRmdFileSimple <- function(x, ...) {
    source(purl(x, output = tempfile()), ...)
}

# ----------------------------------------------------------------------------
# sourceRmdFile
# ----------------------------------------------------------------------------
#' https://gist.github.com/noamross/a549ee50e8a4fd68b8b1
#' Source the R code from an knitr file, optionally skipping plots
#'
#' @param file the knitr file to source
#' @param skip_plots whether to make plots. If TRUE (default) sets a null graphics device
#'
#' @return This function is called for its side effects
#' @export
sourceRmdFile <- function(file, skip_plots = TRUE) {
    temp <- tempfile(fileext = ".R")
    knitr::purl(file, output = temp)

    if (skip_plots) {
        old_dev <- getOption("device")
        options(device = function(...) {
            grDevices::pdf(NULL, ...)
        })
    }
    source(temp)
    if (skip_plots) {
        options(device = old_dev)
    }
}

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
