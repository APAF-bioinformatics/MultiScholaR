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
#' This function saves a given plot in multiple specified formats.
#'
#' @param plot The plot object to be saved
#' @param base_path The base directory path where the plot will be saved
#' @param plot_name The name to be used for the saved plot files
#' @param formats A vector of file formats to save the plot in (default: c("pdf", "png"))
#'
#' @return This function is called for its side effects (saving files)
#' @export
#'
#'
save_plot <- function(plot, base_path, plot_name, formats = c("pdf", "png")) {
  for (format in formats) {
    file_path <- file.path(base_path, "protein_qc", paste0(plot_name, ".", format))
    ggsave(filename = file_path, plot = plot, device = format)
  }
}

##################################################################################################################


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
    # Define all required packages
    required_packages <- c(
        # CRAN packages
        "tidyverse", "seqinr", "lazyeval", "rlang", "glue", "GGally",
        "here", "tibble", "mixOmics", "limma", "magrittr", "future.apply",
        "tictoc", "beepr", "furrr", "readxl", "writexl", "RColorBrewer",
        "multidplyr", "RSpectra", "progress", "Rcpp", "RcppEigen",
        "qvalue", "ruv", "iq", "ggrepel", "patchwork",
        "dplyr", "gtools", "shiny", "DT", "gh", "openxlsx",
        # Additional packages
        "plotly", "vroom", "gplots", "iheatmapr",
        "UpSetR", "gt", "gprofiler2", "htmltools",
        # Git and GitHub related packages
        "git2r", "gh",
        # Bioconductor packages
        "BiocManager",
        # GitHub packages
        "UniProt.ws"
    )
    
    # Install pacman if not present
    if (!requireNamespace("pacman", quietly = TRUE)) {
        if (verbose) message("Installing pacman...")
        utils::install.packages("pacman")
    }
    
    # Install BiocManager if not present
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        if (verbose) message("Installing BiocManager...")
        utils::install.packages("BiocManager")
    }
    
    # Check and install/load packages
    required_packages |>
        purrr::map(function(pkg) {
            if (!requireNamespace(pkg, quietly = TRUE)) {
                if (verbose) message("Installing ", pkg, "...")
                pacman::p_load(char = pkg, character.only = TRUE)
            } else {
                if (verbose) message(pkg, " is already installed, loading...")
                pacman::p_load(char = pkg, character.only = TRUE)
            }
        })
    
    # Handle RUVIIIC separately as it's from GitHub
    if (!requireNamespace("RUVIIIC", quietly = TRUE)) {
        if (verbose) message("Installing RUVIIIC from GitHub...")
        tryCatch({
            devtools::install_github("cran/RUVIIIC")
        }, error = function(e) {
            warning("Failed to install RUVIIIC: ", e$message)
        })
    } else {
        if (verbose) message("RUVIIIC is already installed, loading...")
        pacman::p_load(RUVIIIC)
    }
    
    if (verbose) message("All dependencies loaded successfully!")
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
#' @description Creates and manages project directories with version control
#' @param base_dir Base directory path (optional, defaults to getwd())
#' @param label Optional label to append to proteomics directory name (e.g., "proteomics_MyLabel")
#' @param force Logical; if TRUE, skips user confirmation (default: FALSE)
#' @return List of directory paths
#' @export
setupAndShowDirectories <- function(base_dir = getwd(), label = NULL, force = FALSE) {
    # Create base paths and names
    proteomics_dirname <- if (!is.null(label)) paste0("proteomics_", substr(label, 1, 30)) else "proteomics"
    
    # Check if directories already exist
    results_path <- file.path(base_dir, "results", proteomics_dirname)
    results_summary_path <- file.path(base_dir, "results_summary", proteomics_dirname)
    
    if (!force && (dir.exists(results_path) || dir.exists(results_summary_path))) {
        cat(sprintf("\nWarning: Directory(ies) already exist:\n"))
        if (dir.exists(results_path)) cat(sprintf("- %s\n", results_path))
        if (dir.exists(results_summary_path)) cat(sprintf("- %s\n", results_summary_path))
        
        response <- readline(prompt = "Do you want to overwrite? (y/n): ")
        
        if (tolower(substr(response, 1, 1)) != "y") {
            message("Setup cancelled by user")
            return(invisible(NULL))
        }
        
        # Remove existing directories if user confirmed
        if (dir.exists(results_path)) unlink(results_path, recursive = TRUE)
        if (dir.exists(results_summary_path)) unlink(results_summary_path, recursive = TRUE)
    }
    
    # Define all paths and subdirs in one structure
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
            scripts = file.path(base_dir, "scripts", proteomics_dirname),
            time = file.path(base_dir, "results", proteomics_dirname, "publication_graphs", 
                           "filtering_qc", format(Sys.time(), "%Y%m%d_%H%M%S"))
        )
    )
    
    # Handle existing content and create directories in one pass
    lapply(c("results", "results_summary"), function(type) {
        existing <- dir.exists(file.path(base_dir, type, "proteomics"))
        if (existing) {
            # Get source path
            source_base <- file.path(base_dir, type, "proteomics")
            
            # Copy each subdirectory's contents
            sapply(paths[[type]]$subdirs, function(subdir) {
                src_dir <- file.path(source_base, subdir)
                dest_dir <- file.path(paths[[type]]$base, subdir)
                
                if (dir.exists(src_dir)) {
                    # Create destination directory
                    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
                    
                    # Copy all files maintaining structure
                    files <- list.files(src_dir, full.names = TRUE, recursive = TRUE)
                    if (length(files) > 0) {
                        sapply(files, function(f) {
                            rel_path <- sub(paste0("^", src_dir, "/"), "", f)
                            dest_file <- file.path(dest_dir, rel_path)
                            dir.create(dirname(dest_file), recursive = TRUE, showWarnings = FALSE)
                            file.copy(f, dest_file, overwrite = TRUE)
                        })
                    }
                }
            })
        }
        
        # Create any missing subdirectories
        invisible(sapply(
            file.path(paths[[type]]$base, paths[[type]]$subdirs),
            dir.create, recursive = TRUE, showWarnings = FALSE
        ))
    })
    
    # Create special directories
    invisible(sapply(paths$special, dir.create, recursive = TRUE, showWarnings = FALSE))
    
    # Build and assign global variables
    dir_paths <- c(
        list(base_dir = base_dir),
        setNames(
            lapply(paths$results$subdirs, function(d) file.path(paths$results$base, d)),
            paste0(gsub("-", "_", tolower(paths$results$subdirs)), "_dir")
        ),
        list(
            results_dir = paths$results$base,
            results_summary_dir = paths$results_summary$base,
            data_dir = paths$special$data,
            source_dir = paths$special$scripts,
            time_dir = paths$special$time
        )
    )
    
    # Assign to global environment and print status
    list2env(dir_paths, envir = .GlobalEnv)
    
    # Print simple directory structure with file counts
    cat("\nDirectory Structure:\n")
    invisible(lapply(dir_paths, function(p) {
        if (dir.exists(p)) {
            cat(sprintf("%s (%d files)\n", p, length(list.files(p, recursive = TRUE))))
        }
    }))
    
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
        git_info = "list"
    )
)

#' @title Create Workflow Arguments Container from Config
#' @param workflow_name Name of the workflow
#' @param description Optional description of the workflow
#' @return WorkflowArgs object
#' @export
createWorkflowArgsFromConfig <- function(workflow_name, description = "") {

    # Get git information
    branch_info <- gh("/repos/APAF-BIOINFORMATICS/ProteomeScholaR/branches/dev-jr")
    git_info <- list(
        commit_sha = branch_info$commit$sha,
        branch = "dev-jr",
        repo = "ProteomeScholaR",
        timestamp = branch_info$commit$commit$author$date
    )

    new("WorkflowArgs",
        workflow_name = workflow_name,
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        args = config_list,  # Make sure config_list is defined somewhere
        description = description,
        git_info = git_info
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
                    )
        } else {
            output <- c(output,
                paste0(paste(rep(" ", indent), collapse = ""),
                      name_formatted, ": ", value))
            )
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
            "",
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
            output <- c(header, params, contrasts_header, contrasts_info)
        } else {
            output <- c(header, params)
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
#' @return Invisible list of failed copies for error handling
#' @export
copyToResultsSummary <- function(contrasts_tbl) {
    # Track failed copies
    failed_copies <- list()
    
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
                warning(sprintf("Failed to create directory: %s", dir_path))
                failed_copies[[length(failed_copies) + 1]] <- list(
                    type = "directory_creation",
                    path = dir_path,
                    error = "Failed to create directory"
                )
            }
        })

    # Define files to copy with their display names
    files_to_copy <- list(
        # QC Figures
        list(
            source = file.path(time_dir, "12_correlation_filtered_combined_plots.png"),
            dest = "QC_figures",
            is_dir = FALSE,
            display_name = "Correlation Filtered Plots"
        ),
        # Add pathway directory
        list(
            source = pathway_dir,
            dest = "Publication_figures",
            is_dir = TRUE,
            display_name = "Pathway Analysis"
        ),
        # Add RUV normalized results
        list(
            source = file.path(results_dir, "protein_qc", "ruv_normalised_results_cln_with_replicates.tsv"),
            dest = "Publication_tables",
            is_dir = FALSE,
            display_name = "RUV Normalized Results",
            new_name = "RUV_normalised_results.tsv"
        ),
        list(
            source = file.path(results_dir, "protein_qc", "composite_QC_figure.pdf"),
            dest = "QC_figures",
            is_dir = FALSE,
            display_name = "Composite QC (PDF)"
        ),
        list(
            source = file.path(results_dir, "protein_qc", "composite_QC_figure.png"),
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
            file_parts <- basename(file) |>
                stringr::str_remove("_enrichment_results.tsv") |>
                stringr::str_split("_") |>
                unlist()
            
            contrast <- file_parts[1]
            direction <- file_parts[2]
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
    openxlsx::writeData(enrichment_wb, "Enrichment_Index", enrichment_index_data)
    
    # Apply formatting to enrichment index sheet
    openxlsx::setColWidths(enrichment_wb, "Enrichment_Index", cols = 1:3, widths = c(15, 30, 15))
    openxlsx::addStyle(enrichment_wb, "Enrichment_Index",
                      style = openxlsx::createStyle(textDecoration = "bold"),
                      rows = 1, cols = 1:3)
    
    # Create Publication_tables directory if it doesn't exist
    dir.create(file.path(results_summary_dir, "Publication_tables"), 
               recursive = TRUE, 
               showWarnings = FALSE)
    
    # Save the workbooks
    de_wb_path <- file.path(results_summary_dir, "Publication_tables", "DE_proteins_results.xlsx")
    enrichment_wb_path <- file.path(results_summary_dir, "Publication_tables", "Pathway_enrichment_results.xlsx")
    
    openxlsx::saveWorkbook(de_wb, de_wb_path, overwrite = TRUE)
    openxlsx::saveWorkbook(enrichment_wb, enrichment_wb_path, overwrite = TRUE)

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
                    dest_path <- file.path(dest_dir, basename(source_path))

                    # Create destination directory
                    if (!dir.create(dest_path, showWarnings = FALSE, recursive = TRUE)) {
                        copy_success <- FALSE
                        error_msg <- "Failed to create destination directory"
                    } else {
                        # Copy all files from source to destination
                        files_to_copy <- list.files(source_path, full.names = TRUE, recursive = TRUE)
                        
                        copy_results <- files_to_copy |>
                            sapply(\(f) {
                                rel_path <- sub(paste0("^", source_path, "/"), "", f)
                                dest_file <- file.path(dest_path, rel_path)
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
                } else {
                    source_path <- file_spec$source
                    dest_path <- file.path(dest_dir,
                        if (!is.null(file_spec$new_name)) file_spec$new_name else basename(source_path)
                    )

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

            if (!is.null(file_spec$is_dir) && file_spec$is_dir && source_exists && copy_success) {
                cat(sprintf("%25s Files: %d → %d\n", "",
                    length(source_files),
                    length(dest_files)
                ))
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
#' @param design_matrix A tibble containing the experimental design. Must have a 'group' column and
#'                     equal number of replicates per group.
#' @param config_list A list containing configuration parameters. Must have a 
#'                   'removeRowsWithMissingValuesPercent' nested list with 'groupwise_percentage_cutoff'
#'                   and 'max_groups_percentage_cutoff' parameters.
#' @param min_reps_per_group Integer specifying the minimum number of replicates required in each passing group
#' @param min_groups Integer specifying the minimum number of groups required to have sufficient
#'                  quantifiable values. Default is 2.
#'
#' @return Updated config_list with modified missing value parameters
#' 
#' @details 
#' The function calculates:
#' - groupwise_percentage_cutoff: Allows 1 missing value per group
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
    reps_per_group <- design_matrix |>
        group_by(group) |>
        summarise(n_reps = n()) |>
        pull(n_reps) |>
        unique()
    
    # Ensure all groups have same number of replicates
    if (length(reps_per_group) != 1) {
        stop("All groups must have the same number of replicates")
    }
    
    # Get total number of groups
    total_groups <- design_matrix |>
        pull(group) |>
        unique() |>
        length()
    
    if (min_groups > total_groups) {
        stop("min_groups cannot be larger than total number of groups")
    }
    
    if (min_reps_per_group > reps_per_group) {
        stop("min_reps_per_group cannot be larger than total replicates per group")
    }
    
    # Calculate maximum missing allowed per group
    max_missing_per_group <- reps_per_group - min_reps_per_group
    groupwise_cutoff <- round((max_missing_per_group/reps_per_group) * 100, 3)
    
    # Calculate maximum failing groups allowed
    max_failing_groups <- total_groups - min_groups
    max_groups_cutoff <- round((max_failing_groups/total_groups) * 100, 3)
    
    # Update config_list
    config_list$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff <- groupwise_cutoff
    config_list$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff <- max_groups_cutoff
    config_list$removeRowsWithMissingValuesPercent$proteins_intensity_cutoff_percentile <- 1
    
    # Print informative message
    message(sprintf("Updated missing value parameters in config_list:
    - Requiring at least %d out of %d replicates in each passing group (%.3f%% missing allowed per group)
    - Requiring at least %d out of %d groups to pass (%.3f%% failing groups allowed)
    - groupwise_percentage_cutoff set to %.3f%%
    - max_groups_percentage_cutoff set to %.3f%%", 
        min_reps_per_group,
        reps_per_group,
        groupwise_cutoff,
        min_groups,
        total_groups,
        max_groups_cutoff,
        groupwise_cutoff,
        max_groups_cutoff))
    
    return(config_list)
}
