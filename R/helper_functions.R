# Author(s): Ignatius Pang
# Email: ipang@cmri.org.au
# Children's Medical Research Institute, finding cures for childhood genetic diseases

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

#TODO: Move to an Rpackage RCMRI with common functions. It may be a dependency for each CMRI's project.
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


#TODO: create an Rpackage RCMRI with common functions. It may be a dependency for each CMRI's project.
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
cmriWelcome <- function(name, autors) {
  loginfo("   ______     __    __     ______     __       ")
  loginfo('  /\\  ___\\   /\\ "-./  \\   /\\  == \\   /\\ \\      ')
  loginfo("  \\ \\ \\____  \\ \\ \\-./\\ \\  \\ \\  __<   \\ \\ \\     ")
  loginfo("   \\ \\_____\\  \\ \\_\\ \\ \\_\\  \\ \\_\\ \\_\\  \\ \\_\\    ")
  loginfo("    \\/_____/   \\/_/  \\/_/   \\/_/ /_/   \\/_/    ")
  loginfo("")
  loginfo("---- Children's Medical Research Institute ----")
  loginfo(" Finding cures for childhood genetic diseases  ")
  loginfo("")
  loginfo(" ==============================================")
  loginfo(" %s", name)
  loginfo(" Author(s): %s", paste(autors, sep = ", "))
  loginfo(" cmri-bioinformatics@cmri.org.au")
  loginfo(" ==============================================")
  loginfo("")
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
cmriFormatter <- function(record) { sprintf('CMRI Bioinformatics %s [%s] | %s', record$levelname, record$timestamp, record$msg) }


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
#' @title Load ProteomeScholaR Dependencies
#'
#' @description
#' Installs and loads all required packages for ProteomeScholaR. This includes packages from CRAN, Bioconductor, and GitHub.
#'
#' @param verbose logical; if TRUE (default), displays progress messages during
#'   package installation and loading
#'
#' @details
#' Checks for and installs missing packages, then loads all required
#' dependencies. It handles special cases for GitHub packages (RUVIIIC and
#' ProteomeScholaR) and ensures all necessary packages are available for the DIA
#' workflow.
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
        "git2r", "fs",
        "configr", "webshot2",
        # Added from Suggests:
        "testthat", "ggplot2", "ggpubr", "svglite"
    )

    bioc_packages <- c(
        "UniProt.ws", "mixOmics", "limma", "qvalue",
        "clusterProfiler", "GO.db", # GO.db is often a dependency, ensure it's listed
        # Added from Suggests:
        "EDASeq", "RUVSeq"
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
#' @return WorkflowArgs object
#' @export
createWorkflowArgsFromConfig <- function(workflow_name, description = "", 
                                        organism_name = NULL, taxon_id = NULL) {

    # Get git information
    branch_info <- gh("/repos/APAF-BIOINFORMATICS/ProteomeScholaR/branches/dev-jr")
    git_info <- list(
        commit_sha = branch_info$commit$sha,
        branch = "dev-jr",
        repo = "ProteomeScholaR",
        timestamp = branch_info$commit$commit$author$date
    )
    
    # Get organism information from session if not explicitly provided
    if (is.null(organism_name) && exists("organism_name", envir = .GlobalEnv)) {
        organism_name <- get("organism_name", envir = .GlobalEnv)
    }
    
    if (is.null(taxon_id) && exists("taxon_id", envir = .GlobalEnv)) {
        taxon_id <- get("taxon_id", envir = .GlobalEnv)
    }
    
    # Create organism info list
    organism_info <- list(
        organism_name = organism_name,
        taxon_id = taxon_id
    )

    new("WorkflowArgs",
        workflow_name = workflow_name,
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        args = config_list,  # Make sure config_list is defined somewhere
        description = description,
        git_info = git_info,
        organism_info = organism_info
    )
}

#' @title Format Configuration List
#' @param config_list List of configuration parameters
#' @param indent Number of spaces for indentation
formatConfigList <- function(config_list, indent = 0) {
    output <- character()

    for (name in names(config_list)) {
        value <- config_list[[name]]
        # Skip core_utilisation and complex objects
        if (name == "core_utilisation" ||
            any(class(value) %in% c("process", "R6", "multidplyr_cluster"))) {
            next
        }

        # Format the name
        name_formatted <- gsub("\\.", " ", name)
        name_formatted <- gsub("_", " ", name_formatted)
        name_formatted <- tools::toTitleCase(name_formatted)

        # Handle different value types
        if (is.list(value)) {
            output <- c(output,
                paste0(paste(rep(" ", indent), collapse = ""),
                      name_formatted, ":"))
            output <- c(output,
                formatConfigList(value, indent + 2))
        } else {
            output <- c(output,
                paste0(paste(rep(" ", indent), collapse = ""),
                      name_formatted, ": ", value))
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
            paste("Description:", object@description),
            paste("Timestamp:", object@timestamp),
            "",
            "Git Information:",
            "---------------",
            paste("Repository:", object@git_info$repo),
            paste("Branch:", object@git_info$branch),
            paste("Commit:", object@git_info$commit_sha),
            paste("Git Timestamp:", object@git_info$timestamp),
            ""
        )
        
        # Add organism information if available
        if (!is.null(object@organism_info) && 
            (!is.null(object@organism_info$organism_name) || 
             !is.null(object@organism_info$taxon_id))) {
            organism_header <- c(
                "Organism Information:",
                "---------------------"
            )
            
            organism_details <- character()
            if (!is.null(object@organism_info$organism_name)) {
                organism_details <- c(organism_details, 
                                     paste("Organism Name:", object@organism_info$organism_name))
            }
            if (!is.null(object@organism_info$taxon_id)) {
                organism_details <- c(organism_details, 
                                     paste("Taxon ID:", object@organism_info$taxon_id))
            }
            
            organism_info <- c(organism_header, organism_details, "")
        } else {
            organism_info <- character()
        }

        # Configuration parameters header
        config_header <- c(
            "Configuration Parameters:",
            "------------------------"
        )

        # Format configuration parameters
        params <- formatConfigList(object@args)

        # Add contrasts information if it exists
        if (!is.null(object@args$contrasts_tbl)) {
            contrasts_header <- c(
                "",
                "Contrasts:",
                "----------"
            )
            contrasts_info <- apply(object@args$contrasts_tbl, 1, function(row) {
                paste("  ", row["contrasts"])
            })
            output <- c(header, organism_info, config_header, params, contrasts_header, contrasts_info)
        } else {
            output <- c(header, organism_info, config_header, params)
        }

        # Combine and print
        cat(paste(output, collapse = "\n"), "\n")

        # Save to file if source_dir is defined
        if (exists("source_dir")) {
            output_file <- file.path(source_dir, "study_parameters.txt")
            writeLines(output, output_file)
            cat("\nParameters saved to:", output_file, "\n")
        }
    }
)


##################################################################################################################

#' @title Copy Files to Results Summary and Show Status
#' @description Copies specified files to results summary directory and displays copy status
#' @param contrasts_tbl A tibble containing contrast information
#' @param label Optional label to append to backup directory name (e.g., "backup_MyLabel")
#' @param force Logical; if TRUE, skips user confirmation (default: FALSE)
#' @param current_rmd Optional path to the current .Rmd file being worked on; if NULL (default),
#'                    the function will attempt to detect and save the currently active .Rmd file in RStudio
#' @return Invisible list of failed copies for error handling
#' @export
copyToResultsSummary <- function(contrasts_tbl, label = NULL, force = FALSE, current_rmd = NULL) {
    # Track failed copies
    failed_copies <- list()
    
    # Print current directory paths for debugging
    cat("\nCurrent directory paths:\n")
    cat(sprintf("results_dir: %s\n", results_dir))
    cat(sprintf("results_summary_dir: %s\n", results_summary_dir))
    cat(sprintf("publication_graphs_dir: %s\n", publication_graphs_dir))
    cat(sprintf("time_dir: %s\n", time_dir))
    cat(sprintf("qc_dir: %s\n", qc_dir))
    cat(sprintf("protein_qc_dir: %s\n", protein_qc_dir))
    cat(sprintf("de_output_dir: %s\n", de_output_dir))
    cat(sprintf("pathway_dir: %s\n", pathway_dir))
    cat("\n")
    
    # Try to detect the currently active .Rmd file if none provided
    if (is.null(current_rmd) && exists("rstudioapi") && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
        # Get the active document context
        context <- rstudioapi::getActiveDocumentContext()
        
        # Check if it's an .Rmd file
        if (!is.null(context) && !is.null(context$path) && grepl("\\.rmd$", context$path, ignore.case = TRUE)) {
            current_rmd <- context$path
            cat(sprintf("Detected active .Rmd file: %s\n", current_rmd))
        }
    }
    
    # Handle current Rmd file if provided or detected
    if (!is.null(current_rmd) && file.exists(current_rmd)) {
        # Attempt to save the current Rmd file to ensure latest changes are captured
        tryCatch({
            # This will trigger a save in RStudio if the file is open and has unsaved changes
            if (exists("rstudioapi") && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
                # Get all open documents
                context <- rstudioapi::getActiveDocumentContext()
                
                # Check if our target file is open and active
                if (!is.null(context) && !is.null(context$path) && 
                    normalizePath(context$path) == normalizePath(current_rmd)) {
                    # Save the document
                    rstudioapi::documentSave(context$id)
                    cat(sprintf("Saved current Rmd file: %s\n", current_rmd))
                } else {
                    # Try to find the document among all open documents
                    docs <- rstudioapi::getSourceEditorContexts()
                    for (doc in docs) {
                        if (!is.null(doc$path) && normalizePath(doc$path) == normalizePath(current_rmd)) {
                            rstudioapi::documentSave(doc$id)
                            cat(sprintf("Saved Rmd file: %s\n", current_rmd))
                            break
                        }
                    }
                }
            }
            
            # Copy the Rmd file to the scripts directory
            if (exists("source_dir") && dir.exists(source_dir)) {
                dest_file <- file.path(source_dir, basename(current_rmd))
                file.copy(current_rmd, dest_file, overwrite = TRUE)
                cat(sprintf("Copied Rmd file to scripts directory: %s\n", dest_file))
            } else {
                warning("Scripts directory not found, could not copy Rmd file")
                failed_copies[[length(failed_copies) + 1]] <- list(
                    type = "rmd_copy",
                    source = current_rmd,
                    destination = "scripts directory",
                    display_name = "Current Rmd File",
                    error = "Scripts directory not found"
                )
            }
        }, error = function(e) {
            warning(sprintf("Failed to save/copy Rmd file: %s", e$message))
            failed_copies[[length(failed_copies) + 1]] <- list(
                type = "rmd_copy",
                source = current_rmd,
                destination = "scripts directory",
                display_name = "Current Rmd File",
                error = e$message
            )
        })
    }
    
    # Define target directories
    pub_tables_dir <- file.path(results_summary_dir, "Publication_tables")
    
    # Check if results_summary directory exists
    if (dir.exists(results_summary_dir)) {
        # Create backup directory name with optional label and timestamp
        backup_dirname <- if (!is.null(label)) {
            paste0("proteomics_", substr(label, 1, 30), "_backup_", 
                  format(Sys.time(), "%Y%m%d_%H%M%S"))
        } else {
            paste0("proteomics_backup_", format(Sys.time(), "%Y%m%d_%H%M%S"))
        }
        
        # Create full backup path in the parent directory
        backup_dir <- file.path(dirname(results_summary_dir), backup_dirname)
        
        # Only prompt if force=FALSE
        should_proceed <- if (!force) {
            cat(sprintf("\nResults summary directory exists:\n- %s\n", results_summary_dir))
            repeat {
                response <- readline(prompt = "Do you want to backup existing directory and proceed? (y/n): ")
                response <- tolower(substr(response, 1, 1))
                if (response %in% c("y", "n")) {
                    break
                }
                cat("Please enter 'y' or 'n'\n")
            }
            response == "y"
        } else {
            cat("Force mode enabled - backing up and proceeding...\n")
            TRUE
        }
        
        if (!should_proceed) {
            message("Copy operation cancelled by user")
            return(invisible(NULL))
        }
        
        # Create backup directory
        dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE)
        
        # Copy contents directly to backup directory without creating nested structure
        file.copy(
            list.files(results_summary_dir, full.names = TRUE),
            backup_dir,
            recursive = TRUE
        )
        
        if (length(list.files(backup_dir)) > 0) {
            cat(sprintf("Backed up results_summary directory to: %s\n", backup_dir))
            
            # Create a backup info file
            backup_info <- data.frame(
                original_dir = results_summary_dir,
                backup_time = Sys.time(),
                label = if(!is.null(label)) label else NA_character_,
                stringsAsFactors = FALSE
            )
            
            write.table(
                backup_info,
                file = file.path(backup_dir, "backup_info.txt"),
                sep = "\t",
                row.names = FALSE,
                quote = FALSE
            )
            
            # Remove original directory after successful backup
            unlink(results_summary_dir, recursive = TRUE)
            dir.create(results_summary_dir, recursive = TRUE, showWarnings = FALSE)
        } else {
            warning("Failed to create backup, proceeding without backup")
            failed_copies[[length(failed_copies) + 1]] <- list(
                type = "backup_creation",
                path = backup_dir,
                error = "Failed to create backup"
            )
        }
    }
    
    # Check if Excel files already exist
    de_results_path <- file.path(pub_tables_dir, "DE_proteins_results.xlsx")
    enrichment_path <- file.path(pub_tables_dir, "Pathway_enrichment_results.xlsx")
    
    # Define subdirectories to create
    summary_subdirs <- c(
        "QC_figures",
        "Publication_figures",
        "Publication_tables",
        "Study_report"
    )

    # Create subdirectories in results_summary_dir
    summary_subdirs |>
        sapply(\(subdir) {
            dir_path <- file.path(results_summary_dir, subdir)
            if (!dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)) {
                if (!dir.exists(dir_path)) {  # Only warn if creation failed and dir doesn't exist
                    warning(sprintf("Failed to create directory: %s", dir_path))
                    failed_copies[[length(failed_copies) + 1]] <- list(
                        type = "directory_creation",
                        path = dir_path,
                        error = "Failed to create directory"
                    )
                }
            }
        })

    # Define files to copy with their display names
    files_to_copy <- list(
        # QC Figures - look for correlation plot in multiple locations
        list(
            source = ifelse(
                file.exists(file.path(time_dir, "12_correlation_filtered_combined_plots.png")),
                file.path(time_dir, "12_correlation_filtered_combined_plots.png"),
                # fallback to try finding it in the qc_dir if not in time_dir
                ifelse(
                    file.exists(file.path(qc_dir, "12_correlation_filtered_combined_plots.png")),
                    file.path(qc_dir, "12_correlation_filtered_combined_plots.png"),
                    # ultimate fallback - search for any correlation plots
                    ifelse(
                        length(list.files(qc_dir, pattern = "correlation.*\\.png$", recursive = TRUE, full.names = TRUE)) > 0,
                        list.files(qc_dir, pattern = "correlation.*\\.png$", recursive = TRUE, full.names = TRUE)[1],
                        file.path(time_dir, "12_correlation_filtered_combined_plots.png") # if all else fails, use original path for error reporting
                    )
                )
            ),
            dest = "QC_figures",
            is_dir = FALSE,
            display_name = "Correlation Filtered Plots"
        ),
        # Add RUV normalized results - use absolute paths to avoid confusion with results_dir
        list(
            source = file.path(protein_qc_dir, "ruv_normalised_results_cln_with_replicates.tsv"),
            dest = "Publication_tables",
            is_dir = FALSE,
            display_name = "RUV Normalized Results TSV",
            new_name = "RUV_normalised_results.tsv"
        ),
        # Add RUV RDS file (newly added)
        list(
            source = file.path(protein_qc_dir, "ruv_normalised_results_cln_with_replicates.RDS"),
            dest = "Publication_tables",
            is_dir = FALSE,
            display_name = "RUV Normalized Results RDS",
            new_name = "ruv_normalised_results.RDS"
        ),
        list(
            source = file.path(protein_qc_dir, "composite_QC_figure.pdf"),
            dest = "QC_figures", 
            is_dir = FALSE,
            display_name = "Composite QC (PDF)"
        ),
        list(
            source = file.path(protein_qc_dir, "composite_QC_figure.png"),
            dest = "QC_figures",
            is_dir = FALSE,
            display_name = "Composite QC (PNG)"
        ),
        
        # Publication Figures
        list(
            source = file.path(publication_graphs_dir, "Interactive_Volcano_Plots"),
            dest = "Publication_figures",
            is_dir = TRUE,
            display_name = "Interactive Volcano Plots"
        ),
        list(
            source = file.path(publication_graphs_dir, "NumSigDeMolecules"),
            dest = "Publication_figures",
            is_dir = TRUE,
            display_name = "Num Sig DE Molecules"
        ),
        list(
            source = file.path(publication_graphs_dir, "Volcano_Plots"),
            dest = "Publication_figures",
            is_dir = TRUE,
            display_name = "Volcano Plots"
        ),
        
        # GO Enrichment Plots (newly added)
        list(
            source = file.path(results_dir, "pathway_enrichment"),
            dest = "Publication_figures/Enrichment_Plots",
            is_dir = TRUE,
            display_name = "GO Enrichment Plots"
        ),
        
        # Study Report Tables
        list(
            source = "contrasts_tbl",
            dest = "Study_report",
            type = "object",
            save_as = "contrasts_tbl.tab",
            display_name = "Contrasts Table"
        ),
        list(
            source = "design_matrix",
            dest = "Study_report",
            type = "object",
            save_as = "design_matrix.tab",
            display_name = "Design Matrix"
        ),
        list(
            source = file.path(source_dir, "study_parameters.txt"),
            dest = "Study_report",
            is_dir = FALSE,
            display_name = "Study Parameters"
        )
    )

    # Add all DE protein files dynamically
    de_files <- list.files(
        path = de_output_dir,
        pattern = "de_proteins.*_long_annot\\.xlsx$",
        full.names = TRUE
    )

    # Create a combined workbook for DE results
    de_wb <- openxlsx::createWorkbook()
    
    # Create an index sheet first
    openxlsx::addWorksheet(de_wb, "DE_Results_Index")
    
    # Create index data frame for DE results
    de_index_data <- data.frame(
        Sheet = character(),
        Description = character(),
        stringsAsFactors = FALSE
    )
    
    # Process each DE file
    de_files |>
        purrr::imap(\(file, idx) {
            # Create simple sheet name
            sheet_name <- sprintf("DE_Sheet%d", idx)
            
            # Get base name for index
            base_name <- basename(file) |>
                stringr::str_remove("de_proteins_") |>
                stringr::str_remove("_long_annot.xlsx")
            
            # Add to index
            de_index_data <<- rbind(de_index_data, 
                               data.frame(
                                   Sheet = sheet_name,
                                   Description = base_name,
                                   stringsAsFactors = FALSE
                               ))
            
            # Add data sheet
            data <- openxlsx::read.xlsx(file)
            openxlsx::addWorksheet(de_wb, sheet_name)
            openxlsx::writeData(de_wb, sheet_name, data)
        })
    
    # Write DE index sheet
    openxlsx::writeData(de_wb, "DE_Results_Index", de_index_data)
    
    # Apply formatting to DE index sheet
    openxlsx::setColWidths(de_wb, "DE_Results_Index", cols = 1:2, widths = c(10, 50))
    openxlsx::addStyle(de_wb, "DE_Results_Index", 
                      style = openxlsx::createStyle(textDecoration = "bold"),
                      rows = 1, cols = 1:2)
    
    # Create a new workbook for pathway enrichment results
    enrichment_wb <- openxlsx::createWorkbook()
    
    # Create an index sheet for enrichment results
    openxlsx::addWorksheet(enrichment_wb, "Enrichment_Index")
    
    # Create index data frame for enrichment results
    enrichment_index_data <- data.frame(
        Sheet = character(),
        Contrast = character(),
        Direction = character(),
        stringsAsFactors = FALSE
    )
    
    # Get all enrichment result files
    enrichment_files <- list.files(
        path = pathway_dir,
        pattern = "_enrichment_results.tsv$",
        full.names = TRUE
    )
    
    # Process each enrichment file
    enrichment_files |>
        purrr::imap(\(file, idx) {
            # Extract contrast and direction from filename
            base_name <- basename(file) |>
                stringr::str_remove("_enrichment_results.tsv")
            
            # Extract direction (last part of the base_name which should be "up" or "down")
            direction <- if (stringr::str_ends(base_name, "_up")) {
                "up"
            } else if (stringr::str_ends(base_name, "_down")) {
                "down"
            } else {
                warning("Could not determine direction from filename: ", basename(file))
                "unknown"
            }
            
            # Extract contrast (everything before _up or _down)
            contrast <- stringr::str_replace(base_name, "_(up|down)$", "")
            
            sheet_name <- sprintf("Enrichment_Sheet%d", idx)
            
            # Add to index
            enrichment_index_data <<- rbind(enrichment_index_data,
                                          data.frame(
                                              Sheet = sheet_name,
                                              Contrast = contrast,
                                              Direction = direction,
                                              stringsAsFactors = FALSE
                                          ))
            
            # Add data sheet
            data <- readr::read_tsv(file, show_col_types = FALSE)
            openxlsx::addWorksheet(enrichment_wb, sheet_name)
            openxlsx::writeData(enrichment_wb, sheet_name, data)
        })
    
    # Write enrichment index sheet
    openxlsx::writeDataTable(enrichment_wb, "Enrichment_Index", 
        enrichment_index_data,
        tableStyle = "TableStyleLight9",
        headerStyle = openxlsx::createStyle(textDecoration = "bold"),
        withFilter = TRUE
    )

    # Add an explanatory note at the top of the index sheet
    openxlsx::writeData(enrichment_wb, "Enrichment_Index", 
        data.frame(Note = "Contrast represents the comparison (e.g., Group1_minus_Group2). Direction shows up-regulated or down-regulated genes."),
        startRow = nrow(enrichment_index_data) + 3
    )
    
    # Before saving, ensure the Publication_tables directory exists
    dir.create(pub_tables_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Save the workbooks with a try-catch to handle potential write errors
    tryCatch({
        openxlsx::saveWorkbook(de_wb, de_results_path, overwrite = TRUE)
        cat(sprintf("Successfully saved: %s\n", de_results_path))
    }, error = function(e) {
        failed_copies[[length(failed_copies) + 1]] <- list(
            type = "workbook_save",
            path = de_results_path,
            error = e$message
        )
        warning(sprintf("Failed to save DE results workbook: %s", e$message))
    })
    
    tryCatch({
        openxlsx::saveWorkbook(enrichment_wb, enrichment_path, overwrite = TRUE)
        cat(sprintf("Successfully saved: %s\n", enrichment_path))
    }, error = function(e) {
        failed_copies[[length(failed_copies) + 1]] <- list(
            type = "workbook_save",
            path = enrichment_path,
            error = e$message
        )
        warning(sprintf("Failed to save enrichment workbook: %s", e$message))
    })

    cat("\nCopying files to Results Summary...\n")
    cat("===================================\n\n")

    # Copy and check each file/folder
    files_to_copy |>
        lapply(\(file_spec) {
            dest_dir <- file.path(results_summary_dir, file_spec$dest)
            copy_success <- TRUE
            error_msg <- NULL

            # Get initial status and verify source exists
            if (!is.null(file_spec$type) && file_spec$type == "object") {
                source_exists <- exists(file_spec$source, envir = .GlobalEnv)
                if (!source_exists) {
                    error_msg <- sprintf("Object '%s' not found in global environment", file_spec$source)
                }
            } else {
                source_exists <- if (file_spec$is_dir) {
                    dir.exists(file_spec$source)
                } else {
                    file.exists(file_spec$source)
                }
                if (!source_exists) {
                    error_msg <- sprintf("Source %s not found: %s", 
                        if(file_spec$is_dir) "directory" else "file",
                        file_spec$source)
                }
            }

            # Perform copy operation with detailed error checking
            if (source_exists) {
                if (!is.null(file_spec$type) && file_spec$type == "object") {
                    tryCatch({
                        obj <- get(file_spec$source, envir = .GlobalEnv)
                        dest_path <- file.path(dest_dir, file_spec$save_as)
                        # Ensure destination directory exists
                        dir.create(dirname(dest_path), recursive = TRUE, showWarnings = FALSE)
                        write.table(
                            obj,
                            file = dest_path,
                            sep = "\t",
                            row.names = FALSE,
                            quote = FALSE
                        )
                        # Verify file was created and has content
                        if (!file.exists(dest_path) || file.size(dest_path) == 0) {
                            copy_success <- FALSE
                            error_msg <- "Failed to write object to file or file is empty"
                        }
                    }, error = function(e) {
                        copy_success <- FALSE
                        error_msg <<- sprintf("Error writing object: %s", e$message)
                    })
                } else if (file_spec$is_dir) {
                    source_path <- file_spec$source
                    # Check if destination path includes subdirectories
                    if (grepl("/", file_spec$dest)) {
                        # Extract the last part of the path to use as the custom destination name
                        parts <- strsplit(file_spec$dest, "/")[[1]]
                        dest_dir <- file.path(results_summary_dir, paste(parts[-length(parts)], collapse = "/"))
                        dest_path <- file.path(dest_dir, parts[length(parts)])
                    } else {
                        dest_path <- file.path(dest_dir, basename(source_path))
                    }

                    # Create destination directory
                    if (!dir.create(dest_path, showWarnings = FALSE, recursive = TRUE)) {
                        # If creation failed, check if it already exists
                        if (!dir.exists(dest_path)) {
                            copy_success <- FALSE
                            error_msg <- "Failed to create destination directory"
                        }
                    }
                    
                    if (copy_success) {
                        # Copy all files from source to destination
                        files_to_copy <- list.files(source_path, full.names = TRUE, recursive = TRUE)
                        
                        # Define these variables at the beginning
                        source_files <- list.files(source_path, recursive = TRUE)
                        dest_files <- character(0)
                        
                        if (length(files_to_copy) > 0) {
                            copy_results <- files_to_copy |>
                                sapply(\(f) {
                                    rel_path <- sub(paste0("^", source_path, "/"), "", f)
                                    dest_file <- file.path(dest_path, rel_path)
                                    # Make sure parent directory exists before copying
                                    dir.create(dirname(dest_file), showWarnings = FALSE, recursive = TRUE)
                                    file.copy(f, dest_file, overwrite = TRUE)
                                })
                            
                            if (!all(copy_results)) {
                                copy_success <- FALSE
                                failed_files <- files_to_copy[!copy_results]
                                error_msg <- sprintf("Failed to copy %d files", length(failed_files))
                            }
                            
                            # Verify file counts match
                            source_files <- list.files(source_path, recursive = TRUE)
                            dest_files <- list.files(dest_path, recursive = TRUE)
                            if (length(source_files) != length(dest_files)) {
                                copy_success <- FALSE
                                error_msg <- sprintf("File count mismatch: %d source files, %d destination files",
                                    length(source_files), length(dest_files))
                            }
                        }
                    }
                } else {
                    source_path <- file_spec$source
                    dest_path <- file.path(dest_dir,
                        if (!is.null(file_spec$new_name)) file_spec$new_name else basename(source_path)
                    )
                    
                    # Ensure destination directory exists
                    dir.create(dirname(dest_path), recursive = TRUE, showWarnings = FALSE)

                    if (!file.copy(from = source_path, to = dest_path, overwrite = TRUE)) {
                        copy_success <- FALSE
                        error_msg <- "Failed to copy file"
                    } else {
                        # Verify file sizes match
                        source_size <- file.size(source_path)
                        dest_size <- file.size(dest_path)
                        if (source_size != dest_size) {
                            copy_success <- FALSE
                            error_msg <- sprintf("File size mismatch: source=%d bytes, dest=%d bytes",
                                source_size, dest_size)
                        }
                    }
                }
            }

            # Track failed copies
            if (!source_exists || !copy_success) {
                failed_copies[[length(failed_copies) + 1]] <- list(
                    type = if(file_spec$is_dir) "directory" else "file",
                    source = file_spec$source,
                    destination = dest_dir,
                    display_name = file_spec$display_name,
                    error = error_msg
                )
            }

            # Format and display status with error messages
            source_status <- if (source_exists) "✓" else "✗"
            copy_status <- if (copy_success && source_exists) "✓" else "✗"

            cat(sprintf("%-25s [%s → %s] %s\n",
                file_spec$display_name,
                source_status,
                copy_status,
                if (!is.null(file_spec$is_dir) && file_spec$is_dir) "Directory" else "File"
            ))

            if (!is.null(error_msg)) {
                cat(sprintf("%25s Error: %s\n", "", error_msg))
            }

            # Only display file counts for directories that exist and where we've defined source_files
            if (!is.null(file_spec$is_dir) && file_spec$is_dir && source_exists) {
                # Only try to access source_files and dest_files if they're defined in this scope
                tryCatch({
                    if (exists("source_files", inherits = FALSE) && exists("dest_files", inherits = FALSE)) {
                        cat(sprintf("%25s Files: %d → %d\n", "",
                            length(source_files),
                            length(dest_files)
                        ))
                    }
                }, error = function(e) {
                    cat(sprintf("%25s Files count unavailable: %s\n", "", e$message))
                })
            }
        })

    cat("\nLegend: ✓ = exists/success, ✗ = missing/failed\n")
    cat("Arrow (→) shows source → destination status\n")

    # Report summary of failures
    if (length(failed_copies) > 0) {
        cat("\nFailed Copies Summary:\n")
        cat("=====================\n")
        lapply(failed_copies, function(failure) {
            cat(sprintf("\n%s: %s\n", failure$display_name, failure$error))
            cat(sprintf("Source: %s\n", failure$source))
            cat(sprintf("Destination: %s\n", failure$destination))
        })
        warning(sprintf("%d files/directories failed to copy correctly", length(failed_copies)))
    }

    # Return failed copies list for error handling
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
RenderReport <- function(suffix) {
  # Render this report with specific parameters
  # Using the correct path where study_parameters.txt is located in scripts/proteomics_[suffix]
  rmarkdown::render(file.path(here::here(), "scripts", "proteomics", "DIANN_report.rmd") 
                  ,params = list(suffix = suffix)
                ,output_file = file.path(here::here(), "results_summary", paste0("proteomics_", suffix), paste0("DIANN_report_", suffix, ".docx")))
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