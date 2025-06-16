
##################################################################################################################

#' @import methods
setClass("DirectoryManager",
    slots = c(
        base_dir = "character",
        results_dir = "character",
        data_dir = "character",
        source_dir = "character",
        de_output_dir = "character",
        publication_graphs_dir = "character",
        timestamp = "character",
        qc_dir = "character",
        time_dir = "character",
        results_summary_dir = "character",
        pathway_dir = "character"
    )
)

##################################################################################################################

### Function: create_id_to_attribute_hash
### Description: Create a hash function that map keys to attributes.

## Inputs:
## keys: An array of key values
## attributes: An array of attribute values

## Output:
## An environment that act as a hash to convert keys to attributes.
#' @export
createIdToAttributeHash <- function(keys, attributes) {
	keys <- as.character( as.vector(keys))
	attribute <- as.vector(attributes)

	hash <- new.env(hash = TRUE, parent = parent.frame())

	if ( length(keys) != length(attributes))  {
		warning('Length of keys != Length of attributes list.')
		return(1)
	}

	mapply(function(k, a) base::assign(k, a, envir = hash), 
		   keys, attributes)

	return(hash)
}

##################################################################################################################

### Function: create_id_to_attribute_hash
### Description: Use a predefined hash dictionary to convert any Key to Attribute, return NA if key does not exists

## Inputs:
## key: A key value
## hash: The hash dictionary that maps keys to attributes

## Output:
## A value that correspond to the query key value.
#' @export
convertKeyToAttribute <- function(key, hash) {

  if ( base::exists(key, hash) ) {
    return ( base::get(key, hash))
  }
  else {
    return (NA)
  }
}

##################################################################################################################

#' @export
createDirectoryIfNotExists <- function(file_path, mode = "0777") {

  ### Create directory recursively if it doesn't exist
  if (! file.exists(file_path)){

    dir.create(file_path, showWarnings = TRUE, recursive = TRUE, mode = mode)

  }

}

#' @export
createDirIfNotExists  <- function(file_path, mode = "0777") {
  createDirectoryIfNotExists(file_path, mode = "0777")
}


##################################################################################################################

## Function to source Rmd files
# https://stackoverflow.com/questions/10966109/how-to-source-r-markdown-file-like-sourcemyfile-r
#' @export
sourceRmdFileSimple <- function(x, ...) {
  source(purl(x, output = tempfile()), ...)
}

##################################################################################################################

#' https://gist.github.com/noamross/a549ee50e8a4fd68b8b1
#' Source the R code from an knitr file, optionally skipping plots
#'
#' @param file the knitr file to source
#' @param skip_plots whether to make plots. If TRUE (default) sets a null graphics device
#'
#' @return This function is called for its side effects
#' @export
sourceRmdFile <- function(file, skip_plots = TRUE) {
  temp = tempfile(fileext=".R")
  knitr::purl(file, output=temp)

  if(skip_plots) {
    old_dev = getOption('device')
    options(device = function(...) {
      .Call("R_GD_nullDevice", PACKAGE = "grDevices")
    })
  }
  source(temp)
  if(skip_plots) {
    options(device = old_dev)
  }
}


##################################################################################################################
#=====================================================================================================
#' @export
createOutputDir <- function(output_dir, no_backup) {
  if (output_dir == "") {
    logerror("output_dir is an empty string")
    q()
  }
  if (dir.exists(output_dir)) {
    if (no_backup) {
      unlink(output_dir, recursive = TRUE)
    }
    else {
      backup_name <- paste(output_dir, "_prev", sep = "")
      if (dir.exists(backup_name)) { unlink(backup_name, recursive = TRUE) }
      system(paste("mv", output_dir, backup_name)) }
  }
  dir.create(output_dir, recursive = TRUE)
}


#' @export
testRequiredFiles <- function(files) {
  missing_files <- !file.exists(files)
  invisible(sapply(files[missing_files], function(file) {
    logerror("Missing required file: %s", file)
    q()
  }))
}

#' @export
testRequiredFilesWarning <- function(files) {
  missing_files <- !file.exists(files)
  invisible(sapply(files[missing_files], function(file) {
    logwarn("Missing required file: %s", file)
  }))
}

#' @export
testRequiredArguments <- function(arg_list, parameters) {
  invisible(sapply(parameters, function(par) {
    if (!par %in% names(arg_list)) {
      logerror("Missing required argument: %s", par)
      q()
    }
  }))
}

#' @export
parseType<-function (arg_list,parameters,functType){
  invisible(sapply(parameters, function(key) {
    arg_list[key]<-functType(arg_list[key])
  }))
  return (arg_list)
}

#' @export
parseString<-function (arg_list,parameters){
  invisible(sapply(parameters, function(key) {
    if(key %in% names(arg_list)) {
      arg_list[key] <- str_replace_all(arg_list[key], c("\"" = "", "\'" = ""))
    }
  }))
  return (arg_list)
}

#' @export
parseList<-function (arg_list,parameters){
  invisible(sapply(parameters, function(key) {
    items<-str_replace_all(as.character(arg_list[key])," ","")
    arg_list[key]<- base::strsplit(items,split=",")
  }))
  return (arg_list)
}

#' @export
isArgumentDefined<-function(arg_list,parameter){
  return (!is.null(arg_list[parameter]) & (parameter %in% names(arg_list)) & as.character(arg_list[parameter]) != "")
}



#' @export
setArgsDefault <- function(args, value_name, as_func, default_val=NA ) {

  if(isArgumentDefined(args, value_name))
  {
    args<-parseType(args,
                    c(value_name)
                    ,as_func)
  }else {
    logwarn(paste0( value_name, " is undefined, default value set to ", as.character(default_val), "."))
    args[[ value_name ]] <- default_val
  }

  return(args)
}

#=====================================================================================================


##################################################################################################################


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
savePlot <- function(plot, base_path, plot_name, formats = c("pdf", "png"), width=7, height=7, ... ) {
  saveRDS( plot, file.path(base_path, paste0(plot_name, ".rds")))
  purrr::walk( formats, \(format){
    file_path <- file.path(base_path, paste0(plot_name, ".", format))
    ggsave(filename = file_path, plot = plot, device = format, width=width, height=height, ...)
  })
}



#' Save a plot in multiple formats
#'
#' This function saves a given plot in multiple specified formats and also save the ggplot object as a rds file.
#'
#' @param plot The plot object to be saved
#' @param base_path The base directory path where the plot will be saved
#' @param plot_name The name to be used for the saved plot files
#' @param formats A vector of file formats to save the plot in (default: c("pdf", "png"))
#'#' @param width The width of the plot (default: 7)
#' @param height The height of the plot (default: 7)
#' @param ... Additional arguments to be passed to ggsave
#'
#' @return This function is called for its side effects (saving files)
#' @export
#'
#'
save_plot <- function(plot, base_path, plot_name, formats = c("pdf", "png"), width=7, height=7, ... ) {
  savePlot(plot, base_path, plot_name, formats, width, height, ...)
}


##################################################################################################################

#' Write results to a file
#'
#' This function writes data to a file in the protein_qc subdirectory of the results directory.
#'
#' @param data The data to be written
#' @param filename The name of the file to write the data to
#'
#' @return This function is called for its side effects (writing a file)
#' @export
#'

write_results <- function(data, filename) {
  vroom::vroom_write(data, file.path(results_dir, "protein_qc", filename))
}

##################################################################################################################

#' @export
getFunctionName <- function() {
  calls <- sys.calls()
  current_call <- calls[[length(calls) - 1]]
  as.character(current_call[1])
}


#' @export
getFunctionNameSecondLevel <- function() {
  calls <- sys.calls()
  current_call <- calls[[length(calls) - 2]]
  as.character(current_call[1])
}



#' Check the parameters in the arguments list and the function parameters to see what param applies
#' @export
checkParamsObjectFunctionSimplify <- function(theObject, param_name_string, default_value = NULL) {

  function_name <- getFunctionNameSecondLevel()

  # print(function_name)
  param_value <- dynGet(param_name_string)

  object_value <- (theObject@args)[[function_name]][[param_name_string]]

  # print(paste0("param_value = ", param_value))

  error <- paste0(function_name,  paste0(": '", param_name_string, "' is not defined.\n") )

  if( !is.null(param_value) ) {
    return( param_value)
  } else if( !is.null(object_value) ) {
    # print("use object value")
    return( object_value)
  } else if( !is.null(default_value) ) {
    return( default_value)
  } else {
    stop( error )
  }

}



#' Check the parameters in the arguments list and the function parameters to see what param applies
#' @export
checkParamsObjectFunctionSimplifyAcceptNull <- function(theObject, param_name_string, default_value = NULL) {

  function_name <- getFunctionNameSecondLevel()

  # print(function_name)
  param_value <- dynGet(param_name_string)
  object_value <- theObject@args[[function_name]][[param_name_string]]

  # print(paste0("param_value = ", param_value))

  error <- paste0(function_name, ": '", param_name_string, "' is not defined.\n")

  if( !is.null(param_value) ) {
    return( param_value)
  } else if( !is.null(object_value) ) {
    return( object_value)
  } else if ( !is.null(default_value) ) {
    return( default_value)
  } else {
    warning(error)
    return( NULL )
  }

}

#' Update the parameter in the object
#'@export
updateParamInObject <- function(theObject, param_name_string) {

  function_name <- getFunctionNameSecondLevel()

  theObject@args[[function_name]][[param_name_string]] <- dynGet(param_name_string)

  theObject

}


##################################################################################################################

#' @description Helper function to neatly print out the figures as they get produced
#' @export

summarizeQCPlot <- function(qc_figure) {
            cat("RLE Plots:\n")
            for (plot_name in names(qc_figure@rle_plots)) {
              cat(paste(" -", plot_name, "\n"))
              print(qc_figure@rle_plots[[plot_name]])
            }

            cat("\nPCA Plots:\n")
            for (plot_name in names(qc_figure@pca_plots)) {
              cat(paste(" -", plot_name, "\n"))
              print(qc_figure@pca_plots[[plot_name]])
            }

            cat("\nDensity Plots:\n")
            for (plot_name in names(qc_figure@density_plots)) {
              cat(paste(" -", plot_name, "\n"))
              print(qc_figure@density_plots[[plot_name]])
            }

            cat("\nPearson Correlation Plots:\n")
            for (plot_name in names(qc_figure@pearson_plots)) {
              cat(paste(" -", plot_name, "\n"))
              print(qc_figure@pearson_plots[[plot_name]])
            }
}

##################################################################################################################
#' @export
#' @description Read the config file and return the list of parameters
#' @param file The file path to the config file
#' @param file_type The type of the file (default: "ini")
readConfigFile <- function( file=file.path(source_dir, "config.ini")) {

  config_list <- read.config(file=file, file.type = "ini" )

  # to set the number of cores to be used in the parallel processing
  if("globalParameters" %in% names(config_list)) {
    if ( "number_of_cpus" %in% names( config_list[["globalParameters"]])  ) {

      print(paste0("Read globalParameters: number_of_cpus = "
                   , config_list$globalParameters$number_of_cpus))
      core_utilisation <- new_cluster(config_list$globalParameters$number_of_cpus)
      cluster_library(core_utilisation, c("tidyverse", "glue", "rlang", "lazyeval"))

      list_of_multithreaded_functions <- c("rollUpPrecursorToPeptide"
                                           , "peptideIntensityFiltering"
                                           , "filterMinNumPeptidesPerProtein"
                                           , "filterMinNumPeptidesPerSample"
                                           , "removePeptidesWithOnlyOneReplicate"
                                           , "peptideMissingValueImputation"
                                           , "removeProteinsWithOnlyOneReplicate")

      setCoreUtilisation <- function(config_list, function_name) {
        if (!function_name %in% names(config_list)) {
          config_list[[function_name]] <- list()
        }
        config_list[[function_name]][["core_utilisation"]] <- core_utilisation

        config_list
      }

      for( x in list_of_multithreaded_functions) {
        config_list <- setCoreUtilisation(config_list, x)
      }

      config_list[["globalParameters"]][["plots_format"]] <- str_split(config_list[["globalParameters"]][["plots_format"]], ",")[[1]]

    }}

  getConfigValue <- function (config_list, section, value) {
    config_list[[section]][[value]]
  }

  setConfigValueAsNumeric <- function (config_list, section, value) {
    config_list[[section]][[value]] <- as.numeric(config_list[[section]][[value]])
    config_list
  }

  if("srlQvalueProteotypicPeptideClean" %in% names(config_list)) {
    config_list[["srlQvalueProteotypicPeptideClean"]][["input_matrix_column_ids"]] <- str_split(config_list[["srlQvalueProteotypicPeptideClean"]][["input_matrix_column_ids"]], ",")[[1]]

    print(paste0("Read srlQvalueProteotypicPeptideClean: input_matrix_column_ids = "
                 , paste0(config_list[["srlQvalueProteotypicPeptideClean"]][["input_matrix_column_ids"]]
                          , collapse=", ")))

    config_list <- setConfigValueAsNumeric(config_list
                                           , "srlQvalueProteotypicPeptideClean"
                                           , "qvalue_threshold")
    config_list <- setConfigValueAsNumeric(config_list
                                           , "srlQvalueProteotypicPeptideClean"
                                           , "global_qvalue_threshold")
    config_list <- setConfigValueAsNumeric(config_list
                                           , "srlQvalueProteotypicPeptideClean"
                                           , "choose_only_proteotypic_peptide")

  }


  if("peptideIntensityFiltering" %in% names(config_list)) {
    config_list <- setConfigValueAsNumeric(config_list
                                           , "peptideIntensityFiltering"
                                           , "peptides_intensity_cutoff_percentile")
    config_list <- setConfigValueAsNumeric(config_list
                                           , "peptideIntensityFiltering"
                                           , "peptides_proportion_of_samples_below_cutoff")
  }


  if("filterMinNumPeptidesPerProtein" %in% names(config_list)) {
    config_list <- setConfigValueAsNumeric(config_list
                                           , "filterMinNumPeptidesPerProtein"
                                           , "peptides_per_protein_cutoff")
    config_list <- setConfigValueAsNumeric(config_list
                                           , "filterMinNumPeptidesPerProtein"
                                           , "peptidoforms_per_protein_cutoff")
    # config_list <- setConfigValueAsNumeric(config_list
    #                                        , ""
    #                                        , "")
  }

  if("filterMinNumPeptidesPerSample" %in% names(config_list)) {
    config_list <- setConfigValueAsNumeric(config_list
                                           , "filterMinNumPeptidesPerSample"
                                           , "peptides_per_sample_cutoff")

    if(!"inclusion_list" %in% names( config_list[["filterMinNumPeptidesPerSample"]])) {
      config_list[["filterMinNumPeptidesPerSample"]][["inclusion_list"]] <- ""
    }

    config_list[["filterMinNumPeptidesPerSample"]][["inclusion_list"]] <- str_split(config_list[["filterMinNumPeptidesPerSample"]][["inclusion_list"]], ",")[[1]]

  }

  if("peptideMissingValueImputation" %in% names(config_list)) {
    config_list <- setConfigValueAsNumeric(config_list
                                           , "peptideMissingValueImputation"
                                           , "proportion_missing_values")
  }

  if("removeRowsWithMissingValuesPercent" %in% names(config_list)) {
    config_list <- setConfigValueAsNumeric(config_list
                                           , "removeRowsWithMissingValuesPercent"
                                           , "groupwise_percentage_cutoff")

    config_list <- setConfigValueAsNumeric(config_list
                                           , "removeRowsWithMissingValuesPercent"
                                           , "max_groups_percentage_cutoff")

    config_list <- setConfigValueAsNumeric(config_list
                                           , "removeRowsWithMissingValuesPercent"
                                           , "proteins_intensity_cutoff_percentile")

  }


  if("ruvIII_C_Varying" %in% names(config_list)) {
    config_list <- setConfigValueAsNumeric(config_list
                                           , "ruvIII_C_Varying"
                                           , "ruv_number_k")
  }

  if("plotRle" %in% names(config_list)) {
    config_list[["plotRle"]][["yaxis_limit"]] <- str_split(config_list[["plotRle"]][["yaxis_limit"]], ",")[[1]] |>
      purrr::map_dbl( \(x) as.numeric(x) )

    print(paste0("Read plotRle: yaxis_limit = "
                 , paste0(config_list[["plotRle"]][["yaxis_limit"]], collapse=", ")))
  }

   if("deAnalysisParameters" %in% names(config_list)) {
    # Handle plots_format as array
    config_list[["deAnalysisParameters"]][["plots_format"]] <-
      str_split(config_list[["deAnalysisParameters"]][["plots_format"]], ",")[[1]]

    # Add new lfc_cutoff parameter
    config_list[["deAnalysisParameters"]][["lfc_cutoff"]] <- FALSE

    # Modify treat_lfc_cutoff to use ifelse
    config_list[["deAnalysisParameters"]][["treat_lfc_cutoff"]] <-
      ifelse(config_list[["deAnalysisParameters"]][["lfc_cutoff"]], log2(1.5), 0)

    # Handle args_group_pattern - remove quotes and fix escaping
    if("args_group_pattern" %in% names(config_list[["deAnalysisParameters"]])) {
      config_list[["deAnalysisParameters"]][["args_group_pattern"]] <-
        gsub('^"|"$', '', config_list[["deAnalysisParameters"]][["args_group_pattern"]]) |>
        gsub(pattern = "\\\\", replacement = "\\")
    }

    # Convert numeric parameters
    config_list <- setConfigValueAsNumeric(config_list,
                                         "deAnalysisParameters",
                                         "de_q_val_thresh")

    # Convert boolean parameters
    config_list[["deAnalysisParameters"]][["eBayes_trend"]] <-
      tolower(config_list[["deAnalysisParameters"]][["eBayes_trend"]]) == "true"
    config_list[["deAnalysisParameters"]][["eBayes_robust"]] <-
      tolower(config_list[["deAnalysisParameters"]][["eBayes_robust"]]) == "true"

    print(paste0("Read deAnalysisParameters: formula_string = ",
                config_list[["deAnalysisParameters"]][["formula_string"]]))
}

  config_list
}






#' @export
#' @description Read the config file and specify the section and or parameter to update the object
#' @param theObject The object to be updated
#' @param file The file path to the config file
#' @param section The section to be updated
#' @param value The parameter value to be updated
readConfigFileSection <- function( theObject
                            , file=file.path(source_dir, "config.ini")
                            , function_name
                            , parameter_name = NULL ) {

  config_list <- readConfigFile( file=file
                              , file_type = "ini" )

  if ( is.null(parameter_name) ) {
    theObject@args[[function_name]] <- config_list[[function_name]]
  } else {
    theObject@args[[function_name]][[parameter_name]] <- config_list[[function_name]][[parameter_name]]
  }

  theObject
}


##################################################################################################################
#' @title Load MultiScholaR Dependencies
#'
#' @description
#' Installs and loads all required packages for MultiScholaR This includes packages from CRAN, Bioconductor, and GitHub.
#'
#' @param verbose logical; if TRUE (default), displays progress messages during
#'   package installation and loading
#'
#' @details
#' Checks for and installs missing packages, then loads all required
#' dependencies. It handles special cases for GitHub packages (RUVIIIC and
#' MultiScholaR) and ensures all necessary packages are available for the
#' workflows.
#'
#' @return None (called for side effects)
#'
#' @examples
#' \dontrun{
#' # Load with default verbose messaging
#' loadDependencies()
#'
#' # Load silently
#' loadDependencies(verbose = FALSE)
#' }
#'
#' @importFrom devtools install_github
#' @importFrom pacman p_load
#'
#' @export
loadDependencies <- function(verbose = TRUE) {
    # --- Install Core Managers ---
    if (!requireNamespace("pacman", quietly = TRUE)) {
        if (verbose) message("Installing pacman...")
        utils::install.packages("pacman")
    }
    library(pacman)

    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        if (verbose) message("Installing BiocManager...")
        utils::install.packages("BiocManager")
    }
    # Ensure BiocManager is loaded for installation checks, but don't need library()
    # library(BiocManager) # Not strictly needed just for install()

    if (!requireNamespace("devtools", quietly = TRUE)) {
        if (verbose) message("Installing devtools...")
        utils::install.packages("devtools")
    }
    # library(devtools) # Not strictly needed just for install_github()

    # --- Define Packages by Source ---
    cran_packages <- c(
        "tidyverse", "seqinr", "lazyeval", "rlang", "glue", "GGally",
        "here", "tibble", "magrittr", "future.apply", "tictoc",
        "beepr", "furrr", "readxl", "writexl", "RColorBrewer",
        "multidplyr", "RSpectra", "progress", "Rcpp", "RcppEigen",
        "ruv", "iq", "ggrepel", "patchwork", "dplyr", "gtools",
        "shiny", "DT", "gh", "openxlsx", "plotly", "vroom",
        "gplots", "iheatmapr", "UpSetR", "gt", "gprofiler2",
        "htmltools", "rstudioapi", "flextable", "viridis", "here",
        "git2r", "fs", "logger",
        "configr", "webshot2", "shiny", "shinyjs", "shinyWidgets",
        "shinydashboard", "shinythemes", "shinycssloaders",
        # Added from Suggests:
        "testthat", "ggplot2", "ggpubr", "svglite",
        "ggraph", "reticulate", "shinyFiles"
    )

    bioc_packages <- c(
        "UniProt.ws", "mixOmics", "limma", "qvalue",
        "clusterProfiler", "GO.db", # GO.db is often a dependency, ensure it's listed
        "basilisk"
    )

    github_packages <- list(
        RUVIIIC = "cran/RUVIIIC", # Hosted on CRAN's GitHub mirror? Check source if issues
        Glimma = "APAF-bioinformatics/GlimmaV2"
    )

    # --- Installation and Loading Logic ---

    # Helper function to install/load
    install_and_load <- function(pkg, installer_func, source_name, verbose) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            if (verbose) message(sprintf("Installing %s from %s...", pkg, source_name))
            tryCatch({
                installer_func(pkg)
                # After install, load it
                pacman::p_load(char = pkg, character.only = TRUE)
            }, error = function(e) {
                warning(sprintf("Failed to install %s from %s: %s", pkg, source_name, e$message))
            })
        } else {
            if (verbose) message(sprintf("%s is already installed, loading...", pkg))
            pacman::p_load(char = pkg, character.only = TRUE)
        }
    }

    # CRAN Packages
    if (verbose) message("\n--- Processing CRAN Packages ---")
    purrr::walk(cran_packages, ~install_and_load(
        pkg = .x,
        # Use base R install.packages directly
        installer_func = function(p) utils::install.packages(p, dependencies = TRUE),
        source_name = "CRAN",
        verbose = verbose
    ))

    # Bioconductor Packages
    if (verbose) message("\n--- Processing Bioconductor Packages ---")
    purrr::walk(bioc_packages, ~install_and_load(
        pkg = .x,
        installer_func = function(p) BiocManager::install(p, update = FALSE, ask = FALSE), # Use BiocManager
        source_name = "Bioconductor",
        verbose = verbose
    ))

    # GitHub Packages
    if (verbose) message("\n--- Processing GitHub Packages ---")
    purrr::iwalk(github_packages, ~{
        pkg_name <- .y # Name of the package (e.g., "RUVIIIC")
        repo <- .x     # Repository path (e.g., "cran/RUVIIIC")
        if (!requireNamespace(pkg_name, quietly = TRUE)) {
             if (verbose) message(sprintf("Installing %s from GitHub (%s)...", pkg_name, repo))
             tryCatch({
                 # Force installation to handle potentially corrupt states
                 devtools::install_github(repo, force = TRUE)
                 pacman::p_load(char = pkg_name, character.only = TRUE)
             }, error = function(e) {
                 warning(sprintf("Failed to install %s from GitHub (%s): %s", pkg_name, repo, e$message))
             })
         } else {
             if (verbose) message(sprintf("%s is already installed, loading...", pkg_name))
             pacman::p_load(char = pkg_name, character.only = TRUE)
         }
    })

    if (verbose) message("\nAll dependencies processed successfully!")
}

##################################################################################################################
#' @title Extract Substrings from Underscore-Separated Strings
#'
#' @description
#' Extracts substrings from underscore-separated strings using different modes:
#' range (between positions), start (first element), or end (last element).
#'
#' @param x Character vector containing the strings to process
#' @param mode Character string specifying extraction mode:
#'   * "range": Extract elements between two underscore positions
#'   * "start": Extract from start to first underscore
#'   * "end": Extract from last underscore to end
#' @param start Integer: Starting position for range mode (1-based)
#' @param end Integer: Ending position for range mode (1-based, required for range mode)
#'
#' @return Character vector with extracted strings. Returns NA for strings where
#'   requested positions are out of bounds.
#'
#' @examples
#' x <- "20140602_ffs_expt1_r1_junk"
#' extract_experiment(x, mode = "range", start = 1, end = 3)  # "20140602_ffs_expt1"
#' extract_experiment(x, mode = "start")  # "20140602"
#' extract_experiment(x, mode = "end")    # "junk"
#'
#' # Multiple strings
#' x <- c("20140602_ffs_expt1_r1_junk", "20140603_ffs_expt2_r2_test")
#' extract_experiment(x, mode = "range", start = 2, end = 3)  # c("ffs_expt1", "ffs_expt2")
#'
#' @export
extract_experiment <- function(x, mode = "range", start = 1, end = NULL) {
  if (!mode %in% c("range", "start", "end")) {
    stop("Mode must be one of: 'range', 'start', 'end'")
  }

  process_string <- function(str) {
    parts <- unlist(strsplit(str, "_"))

    if (mode == "range") {
      if (is.null(end)) stop("End position required for range mode")
      if (start > length(parts) || end > length(parts)) {
        warning("Position out of bounds for string: ", str)

        return(NA_character_)
      }
      return(paste(parts[start:end], collapse = "_"))
    }

    else if (mode == "start") {
      return(parts[1])
    }

    else if (mode == "end") {
      return(parts[length(parts)])
    }
  }

  sapply(x, process_string)
}

#' @title Setup Project Directories
#' @description Creates and manages project directories with version control. If directories exist, prompts user to overwrite or reuse.
#' @param base_dir Base directory path (optional, defaults to here::here())
#' @param label Optional label to append to proteomics directory name (e.g., "proteomics_MyLabel")
#' @param force Logical; if TRUE, skips user confirmation and overwrites existing directories (default: FALSE)
#' @return List of directory paths assigned to the global environment.
#' @export
setupAndShowDirectories <- function(base_dir = here::here(), label = NULL, force = FALSE) {
    # --- Define Base Paths and Names ---
    proteomics_dirname <- if (!is.null(label)) paste0("proteomics_", substr(label, 1, 30)) else "proteomics"
    
    # Define all expected paths in one structure for easier management
    paths <- list(
        results = list(
            base = file.path(base_dir, "results", proteomics_dirname),
            subdirs = c("protein_qc", "peptide_qc", "clean_proteins", "de_proteins", 
                       "publication_graphs", "pathway_enrichment",
                       file.path("publication_graphs", "filtering_qc"))
        ),
        results_summary = list(
            base = file.path(base_dir, "results_summary", proteomics_dirname),
            subdirs = c("QC_figures", "Publication_figures", "Publication_tables", "Study_report")
        ),
        special = list(
            data = file.path(base_dir, "data"),
            scripts = file.path(base_dir, "scripts", proteomics_dirname), # Project-specific scripts
            qc_base = file.path(base_dir, "results", proteomics_dirname, "publication_graphs", "filtering_qc")
            # 'time' directory path defined later using timestamp
        )
    )
    
    results_path <- paths$results$base
    results_summary_path <- paths$results_summary$base
    scripts_path <- paths$special$scripts
    
    # Flag to track if we should skip creation/copying and just set vars
    reuse_existing <- FALSE 

    # --- Handle Existing Directories ---
    if (!force && (dir.exists(results_path) || dir.exists(results_summary_path) || dir.exists(scripts_path))) {
        cat(sprintf("\nWarning: Key directory(ies) already exist:\n"))
        if (dir.exists(results_path)) cat(sprintf("- %s\n", results_path))
        if (dir.exists(results_summary_path)) cat(sprintf("- %s\n", results_summary_path))
        if (dir.exists(scripts_path)) cat(sprintf("- %s\n", scripts_path))
        
        response_overwrite <- readline(prompt = "Do you want to overwrite these directories? (y/n): ")
        
        if (tolower(substr(response_overwrite, 1, 1)) == "y") {
            # User chose to overwrite: Delete existing directories
            message("Overwriting existing directories as requested...")
            if (dir.exists(results_path)) unlink(results_path, recursive = TRUE, force = TRUE)
            if (dir.exists(results_summary_path)) unlink(results_summary_path, recursive = TRUE, force = TRUE)
            if (dir.exists(scripts_path)) unlink(scripts_path, recursive = TRUE, force = TRUE)
            # reuse_existing remains FALSE, proceed to create new structure
        } else {
            # User chose NOT to overwrite
            response_reuse <- readline(prompt = "Do you wish to use the existing directory structure for this session? (y/n): ")
            if (tolower(substr(response_reuse, 1, 1)) == "y") {
                message("Reusing existing directory structure and setting environment variables.")
                reuse_existing <- TRUE # Set flag to skip creation/copying
            } else {
                message("Setup cancelled by user. Directory variables not set.")
                return(invisible(NULL)) # Abort if user neither overwrites nor reuses
            }
        }
    } else if (force) {
        message("Force=TRUE: Overwriting existing directories if they exist.")
        if (dir.exists(results_path)) unlink(results_path, recursive = TRUE, force = TRUE)
        if (dir.exists(results_summary_path)) unlink(results_summary_path, recursive = TRUE, force = TRUE)
        if (dir.exists(scripts_path)) unlink(scripts_path, recursive = TRUE, force = TRUE)
    }
    
    # --- Timestamp and Final Path Definitions ---
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    # Define timestamp-specific path (needed regardless of reuse)
    paths$special$time <- file.path(paths$special$qc_base, timestamp)

    # --- Directory Creation and Script Copying ---
    # Only run if NOT reusing existing structure
    if (!reuse_existing) {
        message("Creating directory structure...")
        # Create results and results_summary base and subdirs
    lapply(c("results", "results_summary"), function(type) {
        dir.create(paths[[type]]$base, recursive = TRUE, showWarnings = FALSE)
        invisible(sapply(
            file.path(paths[[type]]$base, paths[[type]]$subdirs),
            dir.create, recursive = TRUE, showWarnings = FALSE
        ))
    })
    
        # Create special directories (ensure qc_base exists for time dir)
        dir.create(paths$special$qc_base, recursive = TRUE, showWarnings = FALSE)
        dir.create(paths$special$data, recursive = TRUE, showWarnings = FALSE)
        dir.create(paths$special$scripts, recursive = TRUE, showWarnings = FALSE)
        dir.create(paths$special$time, recursive = TRUE, showWarnings = FALSE) # Timestamped dir

        # Handle scripts directory copying from template
        scripts_template_source <- file.path(base_dir, "scripts", "proteomics") # Base template location
        if (dir.exists(scripts_template_source)) {
            message("Copying script files (excluding .Rmd) from template...")
            script_files <- list.files(scripts_template_source, full.names = TRUE, recursive = TRUE)
            script_files <- script_files[!grepl("\\.[rR][mM][dD]$", script_files)] # Filter out .Rmd
        
        if (length(script_files) > 0) {
            sapply(script_files, function(f) {
                    # Calculate relative path within the source template
                    rel_path <- sub(paste0("^", tools::file_path_as_absolute(scripts_template_source), "/?"), "", tools::file_path_as_absolute(f))
                    dest_file <- file.path(paths$special$scripts, rel_path) # Destination in project scripts dir
                dir.create(dirname(dest_file), recursive = TRUE, showWarnings = FALSE)
                file.copy(f, dest_file, overwrite = TRUE)
            })
            } else {
                 message("No non-Rmd script files found in template directory.")
            }
        } else {
            message(paste("Script template directory not found at:", scripts_template_source, "- skipping script copy."))
        }
    } else {
         message("Skipping directory creation and script copying as existing structure is being reused.")
         # IMPORTANT: Ensure the timestamped directory for THIS session exists even when reusing.
         dir.create(paths$special$time, recursive = TRUE, showWarnings = FALSE)
    }

    # --- Define Final Directory Paths List for Environment ---
    # This part runs whether creating new or reusing existing structure.
    publication_graphs_dir <- file.path(paths$results$base, "publication_graphs")
    # qc_dir now refers to the timestamped directory base
    qc_dir <- paths$special$qc_base 
    
    dir_paths <- list(
        base_dir = base_dir,
        results_dir = paths$results$base,
        data_dir = paths$special$data,
        source_dir = paths$special$scripts, # This is the *project-specific* scripts dir
        de_output_dir = file.path(paths$results$base, "de_proteins"),
        publication_graphs_dir = publication_graphs_dir,
        timestamp = timestamp,
        qc_dir = qc_dir, # Base QC dir
        time_dir = paths$special$time, # Timestamped dir for current run output
        results_summary_dir = paths$results_summary$base,
        pathway_dir = file.path(paths$results$base, "pathway_enrichment"),
        protein_qc_dir = file.path(paths$results$base, "protein_qc"),
        peptide_qc_dir = file.path(paths$results$base, "peptide_qc"),
        clean_proteins_dir = file.path(paths$results$base, "clean_proteins"),
        qc_figures_dir = file.path(paths$results_summary$base, "QC_figures"),
        publication_figures_dir = file.path(paths$results_summary$base, "Publication_figures"),
        publication_tables_dir = file.path(paths$results_summary$base, "Publication_tables"),
        study_report_dir = file.path(paths$results_summary$base, "Study_report")
    )
    
    # --- Assign to Global Environment ---
    message("Assigning directory paths to global environment...")
    list2env(dir_paths, envir = .GlobalEnv)
    
    # --- Print Structure Summary ---
    cat("\nFinal Directory Structure Set in Global Environment:\n")
    print_paths <- sort(names(dir_paths)) # Sort for consistent order
    invisible(lapply(print_paths, function(name) {
        p <- dir_paths[[name]]
        # Check existence only for paths expected to be directories within the project
        is_project_dir <- startsWith(p, base_dir) && name != "base_dir" && name != "timestamp"
        
        if (is_project_dir && dir.exists(p)) {
            file_count <- length(list.files(p, recursive = TRUE, all.files = TRUE))
            cat(sprintf("%s = %s (%d files/dirs)\n", name, p, file_count))
        } else if (is.character(p)) {
             # Print non-dirs (like timestamp) or dirs outside base_dir (shouldn't happen here)
            cat(sprintf("%s = %s\n", name, p))
        } else {
            cat(sprintf("%s = %s [Non-character path]\n", name, p)) # Should not happen
        }
    }))
    
    # Return the list of paths invisibly
    invisible(dir_paths)
}

##################################################################################################################
#' @import methods
setClass("WorkflowArgs",
    slots = c(
        workflow_name = "character",
        timestamp = "character",
        args = "list",
        description = "character",
        git_info = "list",
        organism_info = "list"
    )
)

#' @title Create Workflow Arguments Container from Config
#' @param workflow_name Name of the workflow
#' @param description Optional description of the workflow
#' @param organism_name Optional organism name (defaults to value from session if available)
#' @param taxon_id Optional taxon ID (defaults to value from session if available)
#' @param source_dir_path Character string, the path to the source directory for this workflow instance.
#' @return WorkflowArgs object
#' @export
createWorkflowArgsFromConfig <- function(workflow_name, description = "", 
                                        organism_name = NULL, taxon_id = NULL,
                                        source_dir_path = NULL) { # New parameter

    # Get git information
    # Attempt to get git info, but allow failure gracefully
    git_info <- tryCatch({
        # Ensure gh is loaded if used directly. Consider adding gh:: to function calls.
        if (requireNamespace("gh", quietly = TRUE)) {
            # Check if a token is available to avoid interactive prompts if not in an interactive session
            if (nzchar(gh::gh_token())) {
                 branch_info <- gh::gh("/repos/APAF-BIOINFORMATICS/MultiScholaR/branches/main") # CHECK REPO/BRANCH
                 list(
                    commit_sha = branch_info$commit$sha,
                    branch = "main", # Consider making this dynamic or a parameter
                    repo = "MultiScholaR", 
                    timestamp = branch_info$commit$commit$author$date
                )
            } else {
                warning("GitHub token not found. Skipping Git information retrieval.", immediate. = TRUE)
                list(commit_sha = NA_character_, branch = NA_character_, repo = NA_character_, timestamp = NA_character_)
            }
        } else {
            warning("Package 'gh' not available. Skipping Git information retrieval.", immediate. = TRUE)
            list(commit_sha = NA_character_, branch = NA_character_, repo = NA_character_, timestamp = NA_character_)
        }
    }, error = function(e) {
        warning(paste("Could not retrieve Git information:", e$message), immediate. = TRUE)
        list(commit_sha = NA_character_, branch = NA_character_, repo = NA_character_, timestamp = NA_character_)
    })
    
    # Get organism information from session if not explicitly provided
    if (is.null(organism_name) && exists("organism_name", envir = .GlobalEnv)) {
        organism_name <- get("organism_name", envir = .GlobalEnv)
    }
    
    if (is.null(taxon_id) && exists("taxon_id", envir = .GlobalEnv)) {
        taxon_id <- get("taxon_id", envir = .GlobalEnv)
    }
    
    # Create organism info list
    organism_info <- list(
        organism_name = if(is.null(organism_name) || !is.character(organism_name)) NA_character_ else organism_name,
        taxon_id = if(is.null(taxon_id) || !(is.character(taxon_id) || is.numeric(taxon_id))) NA_character_ else as.character(taxon_id)
    )

    # Ensure config_list is available (it's used directly in the new() call)
    config_list_to_use <- list() # Default to empty list
    if (exists("config_list", envir = parent.frame()) && is.list(get("config_list", envir = parent.frame()))) {
        config_list_to_use <- get("config_list", envir = parent.frame())
    } else if (exists("config_list", envir = .GlobalEnv) && is.list(get("config_list", envir = .GlobalEnv))) {
        config_list_to_use <- get("config_list", envir = .GlobalEnv)
    } else {
        warning("'config_list' not found as a list in parent frame or global environment. WorkflowArgs @args slot will be empty or use prototype.", immediate. = TRUE)
    }

    new_workflow_args <- new("WorkflowArgs",
        workflow_name = workflow_name,
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        args = config_list_to_use, 
        description = description,
        git_info = git_info,
        organism_info = organism_info
    )
    
    # Store source_dir_path within the 'args' slot using a specific name
    if (!is.null(source_dir_path) && is.character(source_dir_path) && length(source_dir_path) == 1 && nzchar(source_dir_path)) {
        new_workflow_args@args$internal_workflow_source_dir <- tools::file_path_as_absolute(source_dir_path)
    } else {
        new_workflow_args@args$internal_workflow_source_dir <- NULL 
        warning("`source_dir_path` was not provided or invalid for createWorkflowArgsFromConfig. Study parameters might not be saved by show() method.", immediate. = TRUE)
    }
    
    return(new_workflow_args)
}

#' @title Format Configuration List
#' @param config_list List of configuration parameters
#' @param indent Number of spaces for indentation
#' @export
formatConfigList <- function(config_list, indent = 0) {
    output <- character()

    # Exclude internal_workflow_source_dir from printing
    names_to_process <- names(config_list)
    names_to_process <- names_to_process[names_to_process != "internal_workflow_source_dir"]


    for (name in names_to_process) {
        value <- config_list[[name]]
        # Skip core_utilisation and complex objects from display
        if (name == "core_utilisation" ||
            any(class(value) %in% c("process", "R6", "multidplyr_cluster", "cluster", "SOCKcluster"))) { # Added "cluster", "SOCKcluster"
            next
        }

        # Format the name
        name_formatted <- gsub("\\.", " ", name)
        name_formatted <- gsub("_", " ", name_formatted)
        name_formatted <- tools::toTitleCase(name_formatted)

        # Handle different value types
        if (is.list(value)) {
            if (length(value) > 0 && !is.null(names(value))) { # Only print header and recurse if list is named and not empty
                output <- c(output,
                    paste0(paste(rep(" ", indent), collapse = ""),
                          name_formatted, ":"))
                output <- c(output,
                    formatConfigList(value, indent + 2)) # Recursive call
            } else if (length(value) > 0) { # Unnamed list, print elements
                 output <- c(output,
                    paste0(paste(rep(" ", indent), collapse = ""),
                          name_formatted, ":"))
                 for(item_idx in seq_along(value)){
                     item_val <- value[[item_idx]]
                     if(is.atomic(item_val) && length(item_val) == 1){
                        output <- c(output, paste0(paste(rep(" ", indent + 2), collapse=""), "- ", as.character(item_val)))
                     } else {
                        output <- c(output, paste0(paste(rep(" ", indent + 2), collapse=""), "- [Complex List Element]"))
                     }
                 }
            } else { # Empty list
                 output <- c(output,
                    paste0(paste(rep(" ", indent), collapse = ""),
                          name_formatted, ": [Empty List]"))
            }
        } else {
            # Truncate long character vectors for display
            if (is.character(value) && length(value) > 5) {
                value_display <- paste(c(utils::head(value, 5), "..."), collapse = ", ")
            } else if (is.character(value) && length(value) > 1) {
                value_display <- paste(value, collapse = ", ")
            } else {
                value_display <- as.character(value) # Ensure it's character
                 if (length(value_display) == 0) value_display <- "[Empty/NULL]"
            }
            output <- c(output,
                paste0(paste(rep(" ", indent), collapse = ""),
                      name_formatted, ": ", value_display))
        }
    }
    return(output)
}

setMethod("show",
    "WorkflowArgs",
    function(object) {
        # Basic info
        header <- c(
            "Study Parameters",
            "================",
            "",
            "Basic Information:",
            "-----------------",
            paste("Workflow Name:", object@workflow_name),
            paste("Description:", if(nzchar(object@description)) object@description else "N/A"),
            paste("Timestamp:", object@timestamp),
            ""
        )
        
        git_info_vec <- character()
        if (!is.null(object@git_info) && is.list(object@git_info)) {
             gi <- object@git_info
             git_info_vec <- c(
                 paste("Repository:",    ifelse(!is.null(gi$repo) && !is.na(gi$repo), gi$repo, "N/A")),
                 paste("Branch:",         ifelse(!is.null(gi$branch) && !is.na(gi$branch), gi$branch, "N/A")),
                 paste("Commit:",         ifelse(!is.null(gi$commit_sha) && !is.na(gi$commit_sha), substr(gi$commit_sha,1,7), "N/A")),
                 paste("Git Timestamp:", ifelse(!is.null(gi$timestamp) && !is.na(gi$timestamp), gi$timestamp, "N/A"))
             )
        } else {
            git_info_vec <- c("Git Information: N/A")
        }
        git_section <- c("Git Information:", "---------------", git_info_vec, "")
        
        organism_info_vec <- character()
        if (!is.null(object@organism_info) && is.list(object@organism_info)) {
            oi <- object@organism_info
            if (!is.null(oi$organism_name) && !is.na(oi$organism_name) && nzchar(oi$organism_name)) {
                organism_info_vec <- c(organism_info_vec, paste("Organism Name:", oi$organism_name))
            }
            if (!is.null(oi$taxon_id) && !is.na(oi$taxon_id) && nzchar(oi$taxon_id)) {
                organism_info_vec <- c(organism_info_vec, paste("Taxon ID:", oi$taxon_id))
            }
        }
        
        organism_section <- character()
        if (length(organism_info_vec) > 0) {
            organism_section <- c(
                "Organism Information:",
                "---------------------",
                organism_info_vec,
                ""
            )
        }

        config_header <- c(
            "Configuration Parameters (from @args slot):",
            "-----------------------------------------"
        )
        
        args_to_print <- object@args
        # internal_workflow_source_dir is already excluded by formatConfigList

        params <- formatConfigList(args_to_print) 

        # Add contrasts information if it exists and is a data.frame/tibble
        contrasts_section <- character()
        # Check if contrasts_tbl exists and is not NULL and has rows
        if (!is.null(object@args$contrasts_tbl) && 
            (is.data.frame(object@args$contrasts_tbl) || tibble::is_tibble(object@args$contrasts_tbl)) && 
            nrow(object@args$contrasts_tbl) > 0) {
            
            contrasts_header <- c(
                "",
                "Contrasts (from @args$contrasts_tbl):",
                "------------------------------------"
            )
            # Ensure the 'contrasts' column exists before trying to apply
            if ("contrasts" %in% colnames(object@args$contrasts_tbl)) {
                contrasts_info <- apply(object@args$contrasts_tbl, 1, function(row) {
                    paste("  ", row["contrasts"])
                })
                contrasts_section <- c(contrasts_header, contrasts_info)
            } else {
                 contrasts_section <- c(contrasts_header, "  [Column 'contrasts' not found in contrasts_tbl]")
            }
        }

        output <- c(header, git_section, organism_section, config_header, params, contrasts_section)
        cat(paste(output, collapse = "
"), "
")

        # Save to file if internal_workflow_source_dir is defined in args and is valid
        workflow_source_dir_to_use <- object@args$internal_workflow_source_dir

        if (!is.null(workflow_source_dir_to_use) && 
            is.character(workflow_source_dir_to_use) && 
            length(workflow_source_dir_to_use) == 1 && # Ensure it's a single string
            nzchar(workflow_source_dir_to_use) && 
            dir.exists(workflow_source_dir_to_use)) {
            
            output_file <- file.path(workflow_source_dir_to_use, "study_parameters.txt")
            tryCatch({
                writeLines(output, output_file)
                cat("
Parameters saved to:", output_file, "
")
            }, error = function(e) {
                warning(paste("Failed to save study parameters to", output_file, ":", e$message), immediate. = TRUE)
            })
        } else {
            cat("
Warning: `internal_workflow_source_dir` not found, invalid, or directory does not exist in WorkflowArgs @args. Parameters not saved to file.
")
            if(!is.null(workflow_source_dir_to_use)) cat("Attempted path:", workflow_source_dir_to_use, "
")
        }
    }
)


##################################################################################################################

#' @title Copy Files to Results Summary and Show Status
#' @description Copies specified files to results summary directory and displays copy status.
#'              It now derives paths from a global `project_dirs` object using `omic_type` and `experiment_label`.
#' @param omic_type Character string specifying the omics type (e.g., "proteomics", "metabolomics").
#' @param experiment_label Character string specifying the experiment label.
#' @param contrasts_tbl A tibble containing contrast information (typically from the global environment if not passed directly).
#' @param design_matrix A data frame for the design matrix (typically from the global environment if not passed directly).
#' @param force Logical; if TRUE, skips user confirmation for backup (default: FALSE).
#' @param current_rmd Optional path to the current .Rmd file being worked on; if NULL (default),
#'                    the function will attempt to detect and save the currently active .Rmd file in RStudio.
#' @param project_dirs_object_name The name of the list object in the global environment that holds the directory structures (typically the output of setupDirectories). Defaults to "project_dirs".
#' @return Invisible list of failed copies for error handling.
#' @importFrom rlang abort
#' @importFrom tools file_path_as_absolute
#' @export
copyToResultsSummary <- function(omic_type,
                                 experiment_label,
                                 contrasts_tbl = NULL,
                                 design_matrix = NULL,
                                 force = FALSE, 
                                 current_rmd = NULL,
                                 project_dirs_object_name = "project_dirs") {

    # --- Start: Path Derivation and Validation --- #
    if (missing(omic_type) || !is.character(omic_type) || length(omic_type) != 1 || omic_type == "") {
        rlang::abort("`omic_type` must be a single non-empty character string.")
    }
    if (missing(experiment_label) || !is.character(experiment_label) || length(experiment_label) != 1 || experiment_label == "") {
        rlang::abort("`experiment_label` must be a single non-empty character string.")
    }
    if (!exists(project_dirs_object_name, envir = .GlobalEnv)) {
        rlang::abort(paste0("Global object ", sQuote(project_dirs_object_name), " not found. Run setupDirectories() first."))
    }
    
    project_dirs_global <- get(project_dirs_object_name, envir = .GlobalEnv)
    current_omic_key <- paste0(omic_type, "_", experiment_label)

    if (!current_omic_key %in% names(project_dirs_global)) {
        rlang::abort(paste0("Key ", sQuote(current_omic_key), " not found in ", sQuote(project_dirs_object_name), ". Check omic_type and experiment_label."))
    }
    current_paths <- project_dirs_global[[current_omic_key]]
    
    # Validate that current_paths is a list and contains essential directory paths
    required_paths_in_current <- c("results_dir", "results_summary_dir", "publication_graphs_dir", 
                                   "time_dir", "qc_dir", "de_output_dir", "pathway_dir", "source_dir", "feature_qc_dir")
    if (!is.list(current_paths) || !all(required_paths_in_current %in% names(current_paths))) {
        missing_req <- setdiff(required_paths_in_current, names(current_paths))
        rlang::abort(paste0("Essential paths missing from project_dirs for key ", sQuote(current_omic_key), ": ", paste(missing_req, collapse=", ")))
    }
    # --- End: Path Derivation and Validation ---

    # Track failed copies
    failed_copies <- list()
    
    cat("\nRelevant directory paths for: ", current_omic_key, "\n")
    cat(sprintf("Results Dir: %s\n", current_paths$results_dir))
    cat(sprintf("Results Summary Dir: %s\n", current_paths$results_summary_dir))
    cat(sprintf("Publication Graphs Dir: %s\n", current_paths$publication_graphs_dir))
    cat(sprintf("Time Dir (current run): %s\n", current_paths$time_dir))
    cat(sprintf("DE Output Dir: %s\n", current_paths$de_output_dir))
    cat(sprintf("Pathway Dir: %s\n", current_paths$pathway_dir))
    cat(sprintf("Source (Scripts) Dir: %s\n", current_paths$source_dir))
    cat(sprintf("Feature QC Dir: %s\n", current_paths$feature_qc_dir))
    if(!is.null(current_paths$subfeature_qc_dir)) cat(sprintf("Sub-feature QC Dir: %s\n", current_paths$subfeature_qc_dir))
    cat("\n")

    # Try to get contrasts_tbl and design_matrix from the calling environment if NULL
    if (is.null(contrasts_tbl) && exists("contrasts_tbl", envir = parent.frame())) {
        contrasts_tbl <- get("contrasts_tbl", envir = parent.frame())
        message("Using 'contrasts_tbl' from the calling environment.")
    }
    if (is.null(design_matrix) && exists("design_matrix", envir = parent.frame())) {
        design_matrix <- get("design_matrix", envir = parent.frame())
        message("Using 'design_matrix' from the calling environment.")
    }
    
    # Handle current Rmd file
    if (is.null(current_rmd) && exists("rstudioapi") && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
        context <- rstudioapi::getActiveDocumentContext()
        if (!is.null(context) && !is.null(context$path) && grepl("\\.rmd$", context$path, ignore.case = TRUE)) {
            current_rmd <- context$path
            cat(sprintf("Detected active .Rmd file: %s\n", current_rmd))
        }
    }
    
    if (!is.null(current_rmd) && file.exists(current_rmd)) {
        tryCatch({
            if (exists("rstudioapi") && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
                context <- rstudioapi::getActiveDocumentContext()
                if (!is.null(context) && !is.null(context$path) && normalizePath(tools::file_path_as_absolute(context$path)) == normalizePath(tools::file_path_as_absolute(current_rmd))) {
                    rstudioapi::documentSave(context$id)
                    cat(sprintf("Saved current Rmd file: %s\n", current_rmd))
                } else {
                    docs <- rstudioapi::getSourceEditorContexts()
                    for (doc in docs) {
                        if (!is.null(doc$path) && normalizePath(tools::file_path_as_absolute(doc$path)) == normalizePath(tools::file_path_as_absolute(current_rmd))) {
                            rstudioapi::documentSave(doc$id)
                            cat(sprintf("Saved Rmd file: %s\n", current_rmd))
                            break
                        }
                    }
                }
            }
            dest_file <- file.path(current_paths$source_dir, basename(current_rmd))
            file.copy(current_rmd, dest_file, overwrite = TRUE)
            cat(sprintf("Copied Rmd file to project scripts directory: %s\n", dest_file))
        }, error = function(e) {
            warning(sprintf("Failed to save/copy Rmd file: %s", e$message))
            failed_copies[[length(failed_copies) + 1]] <- list(type = "rmd_copy", source = current_rmd, destination = current_paths$source_dir, display_name = "Current Rmd File", error = e$message)
        })
    }
    
    # Define target directories using derived paths
    pub_tables_dir <- file.path(current_paths$results_summary_dir, "Publication_tables")
    
    # Ensure results_summary_dir exists before checking its contents (it should, from setupDirectories)
    if (!dir.exists(current_paths$results_summary_dir)) {
        dir.create(current_paths$results_summary_dir, recursive = TRUE, showWarnings = FALSE)
        logger::log_info("Created results_summary_dir as it was missing: {current_paths$results_summary_dir}")
    }
    
    # Check if the results_summary_dir has any content
    contents_of_summary_dir <- list.files(current_paths$results_summary_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE)
    
    if (length(contents_of_summary_dir) > 0) {
        logger::log_info("Results summary directory for {current_omic_key} ({current_paths$results_summary_dir}) contains existing files/folders.")
        backup_dirname <- paste0(current_omic_key, "_backup_", format(Sys.time(), "%Y%m%d_%H%M%S"))
        backup_dir <- file.path(dirname(current_paths$results_summary_dir), backup_dirname)
        
        should_proceed_with_backup <- if (!force) {
            cat(sprintf("\\nResults summary directory for %s contains content:\\n- %s\\n", current_omic_key, current_paths$results_summary_dir))
            repeat {
                response <- readline(prompt = "Do you want to backup existing directory and proceed by overwriting? (y/n): ")
                response <- tolower(substr(response, 1, 1))
                if (response %in% c("y", "n")) break
                cat("Please enter 'y' or 'n'\\n")
            }
            response == "y"
        } else {
            logger::log_info("Force mode enabled - backing up and proceeding with overwrite for {current_omic_key}...")
            TRUE
        }
        
        if (!should_proceed_with_backup) {
            logger::log_info("Overwrite of {current_paths$results_summary_dir} for {current_omic_key} cancelled by user. No backup made, original files untouched.")
            # Return a list indicating cancellation, which can be checked by the caller.
            return(invisible(list(status="cancelled", omic_key=current_omic_key, message=paste0("Backup and overwrite for ", current_omic_key, " cancelled by user."))))
        }
        
        # Proceed with backup and clearing
        if (!dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE) && !dir.exists(backup_dir)) {
             logger::log_warn("Failed to create backup directory: {backup_dir} for {current_omic_key}. Original directory will not be cleared.")
             failed_copies[[length(failed_copies) + 1]] <- list(type = "backup_dir_creation", path = backup_dir, error = "Failed to create backup directory")
             # Do not proceed with unlink if backup dir creation fails
        } else {
            items_to_backup <- list.files(current_paths$results_summary_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE)
            backup_copy_successful <- TRUE
            if (length(items_to_backup) > 0) {
                dir.create(backup_dir, showWarnings = FALSE, recursive = TRUE) # Ensure backup_dir itself exists
                copy_results <- file.copy(from = items_to_backup, to = backup_dir, overwrite = TRUE, recursive = TRUE, copy.mode = TRUE, copy.date = TRUE)
                if (!all(copy_results)) {
                    backup_copy_successful <- FALSE
                    logger::log_warn("Not all items were successfully copied from {current_paths$results_summary_dir} to backup directory {backup_dir}.")
                }
            }

            backup_has_content <- length(list.files(backup_dir, recursive=TRUE, all.files=TRUE, no..=TRUE)) > 0
            
            if (backup_copy_successful && (backup_has_content || length(contents_of_summary_dir) == 0) ) {
                logger::log_info("Successfully backed up content of {current_paths$results_summary_dir} to: {backup_dir}")
                backup_info <- data.frame(original_dir = current_paths$results_summary_dir, backup_time = Sys.time(), omic_key = current_omic_key, stringsAsFactors = FALSE)
                tryCatch(
                    write.table(backup_info, file = file.path(backup_dir, "backup_info.txt"), sep = "\\t", row.names = FALSE, quote = FALSE),
                    error = function(e) logger::log_warn("Failed to write backup_info.txt: {e$message}")
                )
                
                logger::log_info("Clearing original results_summary_dir: {current_paths$results_summary_dir}")
                unlink_success <- tryCatch({
                    unlink(current_paths$results_summary_dir, recursive = TRUE, force = TRUE)
                    TRUE # Return TRUE on success
                }, error = function(e) {
                    logger::log_error("Error unlinking {current_paths$results_summary_dir}: {e$message}")
                    FALSE # Return FALSE on error
                })

                if(unlink_success) {
                    if (!dir.create(current_paths$results_summary_dir, recursive = TRUE, showWarnings = FALSE) && !dir.exists(current_paths$results_summary_dir)) {
                        warning_msg <- sprintf("CRITICAL: Failed to recreate results_summary_dir %s after backup and unlink. Subsequent operations will likely fail.", current_paths$results_summary_dir)
                        logger::log_error(warning_msg)
                        failed_copies[[length(failed_copies) + 1]] <- list(type = "critical_dir_recreation", path = current_paths$results_summary_dir, error = warning_msg)
                        # Potentially stop execution or return an error status here
                    } else {
                         logger::log_info("Successfully cleared and recreated results_summary_dir: {current_paths$results_summary_dir}")
                    }
                } else {
                     warning_msg <- sprintf("Failed to clear original results_summary_dir %s after backup. Subsequent operations may overwrite or mix files.", current_paths$results_summary_dir)
                     logger::log_warn(warning_msg)
                     failed_copies[[length(failed_copies) + 1]] <- list(type = "dir_clear_failure", path = current_paths$results_summary_dir, error = warning_msg)
                }
            } else {
                logger::log_warn("Failed to copy all items to backup for {current_omic_key}, or backup is unexpectedly empty. Original directory {current_paths$results_summary_dir} was NOT cleared.")
                failed_copies[[length(failed_copies) + 1]] <- list(type = "backup_content_copy", source = current_paths$results_summary_dir, destination = backup_dir, error = "Failed to copy items to backup or backup empty; original not cleared")
            }
        }
    } else {
        logger::log_info("Results summary directory for {current_omic_key} ({current_paths$results_summary_dir}) is empty. No backup needed. Proceeding to create subdirectories.")
        # Ensure the main directory exists (it should if we got here and it was empty, or it was just created if missing)
        if (!dir.exists(current_paths$results_summary_dir)) {
             if(!dir.create(current_paths$results_summary_dir, recursive = TRUE, showWarnings = FALSE)) {
                 warning_msg <- sprintf("CRITICAL: Failed to create initially empty results_summary_dir %s. Subsequent operations may fail.", current_paths$results_summary_dir)
                 logger::log_error(warning_msg)
                 failed_copies[[length(failed_copies) + 1]] <- list(type = "critical_empty_dir_creation", path = current_paths$results_summary_dir, error = warning_msg)
             }
        }
    }
    
    summary_subdirs <- c("QC_figures", "Publication_figures", "Publication_tables", "Study_report")
    sapply(summary_subdirs, \(subdir) {
        dir_path <- file.path(current_paths$results_summary_dir, subdir)
        if (!dir.create(dir_path, recursive = TRUE, showWarnings = FALSE) && !dir.exists(dir_path)) {
            warning(sprintf("Failed to create directory: %s", dir_path))
            failed_copies[[length(failed_copies) + 1]] <- list(type = "directory_creation", path = dir_path, error = "Failed to create directory")
        }
    })

    # Define files to copy - make paths relative to current_paths elements
    files_to_copy <- list(
        list(source = file.path(current_paths$time_dir, "12_correlation_filtered_combined_plots.png"), dest = "QC_figures", is_dir = FALSE, display_name = "Correlation Filtered Plots"),
        list(source = file.path(current_paths$feature_qc_dir, "composite_QC_figure.pdf"), dest = "QC_figures", is_dir = FALSE, display_name = "Composite QC (PDF)", new_name = paste0(omic_type, "_composite_QC_figure.pdf")),
        list(source = file.path(current_paths$feature_qc_dir, "composite_QC_figure.png"), dest = "QC_figures", is_dir = FALSE, display_name = "Composite QC (PNG)", new_name = paste0(omic_type, "_composite_QC_figure.png")),
        list(source = file.path(current_paths$publication_graphs_dir, "Interactive_Volcano_Plots"), dest = "Publication_figures/Interactive_Volcano_Plots", is_dir = TRUE, display_name = "Interactive Volcano Plots"),
        list(source = file.path(current_paths$publication_graphs_dir, "NumSigDeMolecules"), dest = "Publication_figures/NumSigDeMolecules", is_dir = TRUE, display_name = "Num Sig DE Molecules"),
        list(source = file.path(current_paths$publication_graphs_dir, "Volcano_Plots"), dest = "Publication_figures/Volcano_Plots", is_dir = TRUE, display_name = "Volcano Plots"),
        list(source = current_paths$pathway_dir, dest = "Publication_figures/Enrichment_Plots", is_dir = TRUE, display_name = "Pathway Enrichment Plots"),
        list(source = "contrasts_tbl", dest = "Study_report", type = "object", save_as = "contrasts_tbl.tab", display_name = "Contrasts Table"),
        list(source = "design_matrix", dest = "Study_report", type = "object", save_as = "design_matrix.tab", display_name = "Design Matrix"),
        list(source = file.path(current_paths$source_dir, "study_parameters.txt"), dest = "Study_report", is_dir = FALSE, display_name = "Study Parameters")
    )

    if (omic_type == "proteomics") {
        files_to_copy <- c(files_to_copy, list(
            list(source = file.path(current_paths$feature_qc_dir, "ruv_normalised_results_cln_with_replicates.tsv"), dest = "Publication_tables", is_dir = FALSE, display_name = "RUV Normalized Results TSV", new_name = "RUV_normalised_results.tsv"),
            list(source = file.path(current_paths$feature_qc_dir, "ruv_normalised_results_cln_with_replicates.RDS"), dest = "Publication_tables", is_dir = FALSE, display_name = "RUV Normalized Results RDS", new_name = "ruv_normalised_results.RDS")
        ))
    } # Add else if for other omic-specific files if necessary
    
    # Excel files paths
    de_results_excel_path <- file.path(pub_tables_dir, paste0("DE_results_", omic_type, ".xlsx"))
    enrichment_excel_path <- file.path(pub_tables_dir, paste0("Pathway_enrichment_results_", omic_type, ".xlsx"))

    # Create combined DE workbook
    de_wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(de_wb, "DE_Results_Index")
    de_index_data <- data.frame(Sheet = character(), Description = character(), stringsAsFactors = FALSE)
    de_files <- list.files(path = current_paths$de_output_dir, pattern = paste0("de_", "\\w+", "_long_annot\\.xlsx$"), full.names = TRUE) # More generic pattern
    
    purrr::imap(de_files, \(file, idx) {
        sheet_name <- sprintf("DE_Sheet%d", idx)
        base_name <- basename(file) |> stringr::str_remove(paste0("de_", omic_type, "_|", "_long_annot\\.xlsx$"))
        de_index_data <<- rbind(de_index_data, data.frame(Sheet = sheet_name, Description = base_name, stringsAsFactors = FALSE))
        data_content <- tryCatch(openxlsx::read.xlsx(file), error = function(e) NULL)
        if (!is.null(data_content)) {
             openxlsx::addWorksheet(de_wb, sheet_name)
             openxlsx::writeData(de_wb, sheet_name, data_content)
        } else {
            warning(paste0("Failed to read DE Excel file: ", file))
            failed_copies[[length(failed_copies) + 1]] <- list(type = "de_excel_read", source = file, error = "Failed to read")
        }
    })
    openxlsx::writeData(de_wb, "DE_Results_Index", de_index_data)
    openxlsx::setColWidths(de_wb, "DE_Results_Index", cols = 1:2, widths = c(15, 50))
    openxlsx::addStyle(de_wb, "DE_Results_Index", style = openxlsx::createStyle(textDecoration = "bold"), rows = 1, cols = 1:2)

    # Create combined Enrichment workbook
    enrichment_wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(enrichment_wb, "Enrichment_Index")
    enrichment_index_data <- data.frame(Sheet = character(), Contrast = character(), Direction = character(), stringsAsFactors = FALSE)
    enrichment_files <- list.files(path = current_paths$pathway_dir, pattern = "_enrichment_results\\.tsv$", full.names = TRUE)
    
    purrr::imap(enrichment_files, \(file, idx) {
        base_name <- basename(file) |> stringr::str_remove("_enrichment_results\\.tsv$")
        direction <- ifelse(stringr::str_ends(base_name, "_up"), "up", ifelse(stringr::str_ends(base_name, "_down"), "down", "unknown"))
        contrast_label <- stringr::str_replace(base_name, "_(up|down)$", "")
        sheet_name <- sprintf("Enrich_Sheet%d", idx)
        enrichment_index_data <<- rbind(enrichment_index_data, data.frame(Sheet = sheet_name, Contrast = contrast_label, Direction = direction, stringsAsFactors = FALSE))
        data_content <- tryCatch(readr::read_tsv(file, show_col_types = FALSE), error = function(e) NULL)
         if (!is.null(data_content)) {
            openxlsx::addWorksheet(enrichment_wb, sheet_name)
            openxlsx::writeData(enrichment_wb, sheet_name, data_content)
        } else {
            warning(paste0("Failed to read Enrichment TSV file: ", file))
            failed_copies[[length(failed_copies) + 1]] <- list(type = "enrichment_tsv_read", source = file, error = "Failed to read")
        }
    })
    openxlsx::writeDataTable(enrichment_wb, "Enrichment_Index", enrichment_index_data, tableStyle = "TableStyleLight9", headerStyle = openxlsx::createStyle(textDecoration = "bold"), withFilter = TRUE)
    openxlsx::writeData(enrichment_wb, "Enrichment_Index", data.frame(Note = "Contrast represents the comparison (e.g., Group1_minus_Group2). Direction shows up-regulated or down-regulated genes."), startRow = nrow(enrichment_index_data) + 3)
    
    dir.create(pub_tables_dir, recursive = TRUE, showWarnings = FALSE)
    tryCatch({ openxlsx::saveWorkbook(de_wb, de_results_excel_path, overwrite = TRUE); cat(sprintf("Successfully saved DE results to: %s\n", de_results_excel_path)) }, error = function(e) { failed_copies[[length(failed_copies) + 1]] <- list(type = "workbook_save", path = de_results_excel_path, error = e$message); warning(sprintf("Failed to save DE workbook: %s", e$message)) })
    tryCatch({ openxlsx::saveWorkbook(enrichment_wb, enrichment_excel_path, overwrite = TRUE); cat(sprintf("Successfully saved Enrichment results to: %s\n", enrichment_excel_path)) }, error = function(e) { failed_copies[[length(failed_copies) + 1]] <- list(type = "workbook_save", path = enrichment_excel_path, error = e$message); warning(sprintf("Failed to save Enrichment workbook: %s", e$message)) })

    cat("\nCopying individual files/folders to Results Summary for ", current_omic_key, "...\n")
    cat("==============================================================\n\n")

    files_to_copy |>
        lapply(\(file_spec) {
            dest_dir_final <- file.path(current_paths$results_summary_dir, file_spec$dest)
            copy_success <- TRUE
            error_msg <- NULL
            source_display <- file_spec$source # For display purposes

            if (!is.null(file_spec$type) && file_spec$type == "object") {
                source_exists <- !is.null(get0(file_spec$source, envir = parent.frame())) # Use get0 to avoid error if object doesn't exist
                if (!source_exists) error_msg <- sprintf("Object '%s' not found in parent/global environment", file_spec$source)
            } else {
                source_exists <- if (file_spec$is_dir) dir.exists(file_spec$source) else file.exists(file_spec$source)
                if (!source_exists) error_msg <- sprintf("Source %s not found: %s", if(file_spec$is_dir) "directory" else "file", file_spec$source)
            }

            if (source_exists) {
                dir.create(dest_dir_final, recursive = TRUE, showWarnings = FALSE)
                if (!is.null(file_spec$type) && file_spec$type == "object") {
                    tryCatch({
                        obj <- get(file_spec$source, envir = parent.frame())
                        dest_path <- file.path(dest_dir_final, file_spec$save_as)
                        write.table(obj, file = dest_path, sep = "\t", row.names = FALSE, quote = FALSE)
                        if (!file.exists(dest_path) || (file.exists(dest_path) && file.size(dest_path) == 0 && nrow(obj) > 0) ) {
                            copy_success <- FALSE
                            error_msg <- "Failed to write object or file is empty"
                        }
                    }, error = function(e) { copy_success <<- FALSE; error_msg <<- sprintf("Error writing object: %s", e$message) })
                } else if (file_spec$is_dir) {
                    # For directories, copy contents. Destination is dest_dir_final itself if not nested, or specific if dest includes subdirs.
                    # The file_spec$dest might be "Publication_figures/Enrichment_Plots"
                    # In this case dest_dir_final is already .../results_summary_dir/Publication_figures/Enrichment_Plots
                    all_source_files <- list.files(file_spec$source, full.names = TRUE, recursive = TRUE)
                    if (length(all_source_files) > 0) {
                        copy_results <- sapply(all_source_files, function(f) {
                            rel_path_from_source_base <- sub(paste0("^", tools::file_path_as_absolute(file_spec$source), "/?"), "", tools::file_path_as_absolute(f))
                            dest_file_abs <- file.path(dest_dir_final, rel_path_from_source_base)
                            dir.create(dirname(dest_file_abs), showWarnings = FALSE, recursive = TRUE)
                            file.copy(f, dest_file_abs, overwrite = TRUE)
                        })
                        if (!all(copy_results)) {
                            copy_success <- FALSE
                            error_msg <- sprintf("Failed to copy some files from %s", file_spec$source)
                        }
                    } else {
                        # Source directory might be empty, which is not an error for the copy operation itself
                        message(sprintf("Source directory %s is empty. Nothing to copy.", file_spec$source))
                    }
                } else {
                    dest_path <- file.path(dest_dir_final, if (!is.null(file_spec$new_name)) file_spec$new_name else basename(file_spec$source))
                    if (!file.copy(from = file_spec$source, to = dest_path, overwrite = TRUE)) {
                        copy_success <- FALSE
                        error_msg <- "Failed to copy file"
                    } else if (file.exists(file_spec$source) && file.exists(dest_path) && file.size(file_spec$source) != file.size(dest_path)) {
                        copy_success <- FALSE
                        error_msg <- sprintf("File size mismatch: source=%d, dest=%d bytes", file.size(file_spec$source), file.size(dest_path))
                    }
                }
            }

            if (!source_exists || !copy_success) {
                failed_copies[[length(failed_copies) + 1]] <- list(type = if(!is.null(file_spec$type) && file_spec$type == "object") "object" else if(file_spec$is_dir) "directory" else "file", source = source_display, destination = dest_dir_final, display_name = file_spec$display_name, error = error_msg)
            }
            cat(sprintf("%-35s [%s -> %s] %s\n", file_spec$display_name, if(source_exists) "✓" else "✗", if(copy_success && source_exists) "✓" else "✗", if(!is.null(file_spec$type) && file_spec$type == "object") "Object" else if(file_spec$is_dir) "Directory" else "File"))
            if (!is.null(error_msg)) cat(sprintf("%35s Error: %s\n", "", error_msg))
        })

    cat("\nLegend: ✓ = exists/success, ✗ = missing/failed\n")
    cat("Arrow (->) shows source -> destination status\n")

    if (length(failed_copies) > 0) {
        cat("\nFailed Copies Summary for ", current_omic_key, ":\n")
        cat("=====================================\n")
        lapply(failed_copies, function(failure) {
            cat(sprintf("\n%s: %s\n", failure$display_name, failure$error))
            cat(sprintf("  Source: %s\n", as.character(failure$source))) # Ensure source is char
            cat(sprintf("  Destination Attempted: %s\n", failure$destination))
        })
        warning(sprintf("%d files/objects/directories failed to copy correctly for %s", length(failed_copies), current_omic_key))
    }
    cat("--- End of copyToResultsSummary for ", current_omic_key, " ---\n\n")
    invisible(failed_copies)
}


#' Update Missing Value Parameters in Configuration List
#'
#' @description
#' Automatically calculates and updates the missing value filtering parameters in the configuration list
#' based on the experimental design matrix. The function ensures at least a specified number of groups
#' have sufficient quantifiable values for analysis.
#'
#' @param design_matrix A tibble containing the experimental design. Must have a 'group' column.
#'                     Groups can have different numbers of replicates.
#' @param config_list A list containing configuration parameters. Must have a
#'                   'removeRowsWithMissingValuesPercent' nested list with 'groupwise_percentage_cutoff'
#'                   and 'max_groups_percentage_cutoff' parameters.
#' @param min_reps_per_group Integer specifying the minimum number of replicates required in each passing group.
#'                          If a group has fewer total replicates than this value, the minimum is adjusted.
#' @param min_groups Integer specifying the minimum number of groups required to have sufficient
#'                  quantifiable values. Default is 2.
#'
#' @return Updated config_list with modified missing value parameters
#'
#' @details
#' The function calculates:
#' - groupwise_percentage_cutoff: Based on minimum required replicates per group
#' - max_groups_percentage_cutoff: Based on minimum required groups
#'
#' @examples
#' \dontrun{
#' config_list <- updateMissingValueParameters(design_matrix, config_list, min_reps_per_group = 2)
#' }
#'
#' @export
updateMissingValueParameters <- function(design_matrix, config_list, min_reps_per_group = 2, min_groups = 2) {
    # Get number of replicates per group
    reps_per_group_tbl <- design_matrix |>
        group_by(group) |>
        summarise(n_reps = n()) |>
        ungroup()
    
    # Get total number of groups
    total_groups <- nrow(reps_per_group_tbl)
    
    if (min_groups > total_groups) {
        stop("min_groups cannot be larger than total number of groups")
    }
    
    # Calculate percentage missing allowed for each group
    group_thresholds <- reps_per_group_tbl |>
        mutate(
            adjusted_min_reps = pmin(n_reps, min_reps_per_group),
            max_missing = n_reps - adjusted_min_reps,
            missing_percent = round((max_missing / n_reps) * 100, 3)
        )
    
    # Use a consistent percentage threshold across all groups
    # Take the minimum percentage to ensure all groups meet minimum requirements
    groupwise_cutoff <- max(group_thresholds$missing_percent)
    
    # Calculate maximum failing groups allowed
    max_failing_groups <- total_groups - min_groups
    max_groups_cutoff <- round((max_failing_groups / total_groups) * 100, 3)
    
    # Update config_list
    config_list$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff <- groupwise_cutoff
    config_list$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff <- max_groups_cutoff
    config_list$removeRowsWithMissingValuesPercent$proteins_intensity_cutoff_percentile <- 1
    
    # Create informative message
    basic_msg <- sprintf(
        "Updated missing value parameters in config_list:
    - Requiring at least %d replicates in each passing group (varies by group size)
    - Requiring at least %d out of %d groups to pass (%.3f%% failing groups allowed)
    - groupwise_percentage_cutoff set to %.3f%%
    - max_groups_percentage_cutoff set to %.3f%%",
        min_reps_per_group,
        min_groups,
        total_groups,
        max_groups_cutoff,
        groupwise_cutoff,
        max_groups_cutoff
    )
    
    # Add details for each group
    group_detail_strings <- group_thresholds |>
        mutate(
            detail = sprintf("    Group %s: %d out of %d replicates required (%.3f%% missing allowed)",
                             group, adjusted_min_reps, n_reps, missing_percent)
        ) |>
        dplyr::pull(detail)
        
    group_details <- paste(group_detail_strings, collapse = "\n")
    
    # Print the message
    message(paste(basic_msg, "\n\nGroup details:", group_details, sep = "\n"))
    
    return(config_list)
}

##################################################################################################################
#' @export
updateRuvParameters <- function(config_list, best_k, control_genes_index, percentage_as_neg_ctrl) {
  config_list$ruvParameters$best_k <- best_k
  config_list$ruvParameters$num_neg_ctrl <- length(control_genes_index)
  config_list$ruvParameters$percentage_as_neg_ctrl <- percentage_as_neg_ctrl
  
  # Print the number of negative controls (as in the original code)
  config_list$ruvParameters$num_neg_ctrl
  
  # Return the updated config list
  return(config_list)
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

    # --- Retrieve Paths from Global Project Directories Object ---
    if (!exists(project_dirs_object_name, envir = .GlobalEnv)) {
        rlang::abort(paste0("Global object ", sQuote(project_dirs_object_name), " not found. Run setupDirectories() first."))
    }
    
    project_dirs_global <- get(project_dirs_object_name, envir = .GlobalEnv)
    current_omic_key <- paste0(omic_type, "_", experiment_label)

    if (!current_omic_key %in% names(project_dirs_global)) {
        rlang::abort(paste0("Key ", sQuote(current_omic_key), " not found in ", 
                           sQuote(project_dirs_object_name), 
                           ". Check omic_type ('", omic_type, "') and experiment_label ('", experiment_label, "')."))
    }
    current_paths <- project_dirs_global[[current_omic_key]] # This contains base_dir, results_summary_dir etc. for the *labelled* omic

    if (!is.list(current_paths) || 
        is.null(current_paths$base_dir) || # Need base_dir to find the template Rmd
        is.null(current_paths$results_summary_dir)) {
        rlang::abort(paste0("Essential paths (base_dir, results_summary_dir) missing for key ", 
                           sQuote(current_omic_key), " in ", sQuote(project_dirs_object_name), "."))
    }

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
    
    if (!file.exists(rmd_input_path)) {
        rlang::abort(paste0("R Markdown template file not found at the expected location: ", sQuote(rmd_input_path),
                           ". This should be in the general scripts/<omic_type> directory (e.g., scripts/proteomics)."))
    }

    # --- Construct Output Path (in the labelled results_summary directory) ---
    output_file_basename <- paste0(tools::file_path_sans_ext(rmd_filename), 
                                   "_", omic_type, 
                                   "_", experiment_label)
                                   
    output_ext <- ".docx" # Default
    if (!is.null(output_format)){
        if(output_format == "word_document" || grepl("word", output_format, ignore.case=TRUE)){
            output_ext <- ".docx"
        } else if (output_format == "html_document" || grepl("html", output_format, ignore.case=TRUE)){
            output_ext <- ".html"
        } else if (grepl("pdf", output_format, ignore.case=TRUE)){
            output_ext <- ".pdf"
        }
    }

    output_file_path <- file.path(current_paths$results_summary_dir, paste0(output_file_basename, output_ext))

    logger::log_info("Attempting to render report:")
    logger::log_info("- Rmd Source (Template): {rmd_input_path}")
    logger::log_info("- Output File: {output_file_path}")
    logger::log_info("- Params: omic_type=\'{omic_type}\', experiment_label=\'{experiment_label}\'")

    # --- Render the Report ---
    rendered_path <- tryCatch({
        rmarkdown::render(
            input = rmd_input_path,
            params = list(
                omic_type = omic_type,
                experiment_label = experiment_label
            ),
            output_file = output_file_path,
            output_format = output_format, # Pass this along; if NULL, Rmd default is used
            envir = new.env(parent = globalenv()) # Render in a clean environment
        )
    }, error = function(e) {
        logger::log_error("Failed to render R Markdown report: {e$message}")
        logger::log_error("Input path: {rmd_input_path}")
        logger::log_error("Output path: {output_file_path}")
        NULL # Return NULL on failure
    })

    if (!is.null(rendered_path) && file.exists(rendered_path)) {
        logger::log_info("Report successfully rendered to: {rendered_path}")
    } else {
        logger::log_warn("Report rendering failed or output file not found at expected location.")
    }
    
    invisible(rendered_path)
}

#' @title Update Parameter in S4 Object Args and Global Config List
#' @description Modifies a specific parameter within an S4 object's @args slot
#'              and also updates the corresponding value in a global list named
#'              'config_list'.
#'
#' @param theObject The S4 object whose @args slot needs updating.
#' @param function_name The name identifying the parameter section (character string,
#'                      e.g., "peptideIntensityFiltering"). Corresponds to the
#'                      first-level key in both @args and config_list.
#' @param parameter_name The specific parameter name to update (character string,
#'                       e.g., "peptides_proportion_of_samples_below_cutoff").
#'                       Corresponds to the second-level key.
#' @param new_value The new value to assign to the parameter.
#' @param config_list_name The name of the global list variable holding the
#'                         configuration (defaults to "config_list").
#' @param env The environment where the global config list resides (defaults to
#'            .GlobalEnv).
#'
#' @return The modified S4 object.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assume 'myPeptideData' is a PeptideQuantitativeData object
#' # Assume 'config_list' exists in the global environment
#'
#' # Check initial values (example)
#' # print(myPeptideData@args$peptideIntensityFiltering$peptides_proportion_of_samples_below_cutoff)
#' # print(config_list$peptideIntensityFiltering$peptides_proportion_of_samples_below_cutoff)
#'
#' # Update the parameter to 0.7
#' myPeptideData <- updateConfigParameter(
#'   theObject = myPeptideData,
#'   function_name = "peptideIntensityFiltering",
#'   parameter_name = "peptides_proportion_of_samples_below_cutoff",
#'   new_value = 0.7
#' )
#'
#' # Verify changes (example)
#' # print(myPeptideData@args$peptideIntensityFiltering$peptides_proportion_of_samples_below_cutoff) # Should be 0.7
#' # print(config_list$peptideIntensityFiltering$peptides_proportion_of_samples_below_cutoff) # Should be 0.7
#' }
updateConfigParameter <- function(theObject,
                                function_name,
                                parameter_name,
                                new_value,
                                config_list_name = "config_list",
                                env = .GlobalEnv) {

  # --- Input Validation ---
  if (!isS4(theObject)) {
    stop("'theObject' must be an S4 object.")
  }
  if (!"args" %in% methods::slotNames(theObject)) {
      stop("'theObject' must have an '@args' slot.")
  }
  if (!is.character(function_name) || length(function_name) != 1) {
    stop("'function_name' must be a single character string.")
  }
  if (!is.character(parameter_name) || length(parameter_name) != 1) {
    stop("'parameter_name' must be a single character string.")
  }
  if (!exists(config_list_name, envir = env)) {
      stop("Global config list '", config_list_name, "' not found in the specified environment.")
  }

  # Retrieve the global list safely
  current_config_list <- get(config_list_name, envir = env)

  if (!is.list(current_config_list)) {
      stop("Global variable '", config_list_name, "' is not a list.")
  }
   if (!function_name %in% names(current_config_list)) {
    warning("Function name '", function_name, "' not found in global config list '", config_list_name, "'. Adding it.")
    current_config_list[[function_name]] <- list()
  }
  if (!parameter_name %in% names(current_config_list[[function_name]])) {
       warning("Parameter '", parameter_name, "' not found under '", function_name, "' in global config list '", config_list_name, "'. Adding it.")
  }


  # --- Update S4 Object @args ---
  if (is.null(theObject@args)) {
       theObject@args <- list() # Initialize args if it's NULL
   }
  if (!is.list(theObject@args[[function_name]])) {
     # Initialize the sub-list if it doesn't exist or isn't a list
     theObject@args[[function_name]] <- list()
  }
  theObject@args[[function_name]][[parameter_name]] <- new_value
  message("Updated @args$", function_name, "$", parameter_name, " in S4 object.")


  # --- Update Global Config List ---
  current_config_list[[function_name]][[parameter_name]] <- new_value
  # Assign the modified list back to the global environment
  assign(config_list_name, current_config_list, envir = env)
  message("Updated ", config_list_name, "$", function_name, "$", parameter_name, " in global environment.")

  # --- Return the modified object ---
  return(theObject)
}