# Author(s): Ignatius Pang
# Email: ipang@cmri.org.au
# Children’s Medical Research Institute, finding cures for childhood genetic diseases

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

	for ( i in 1:length(keys) )  {
		base::assign( keys[i], attributes[i], envir = hash)
	}

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
  loginfo("---- Children’s Medical Research Institute ----")
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
  for (file in files[missing_files]) {
    logerror("Missing required file: %s", file)
    q()
  }
}

#' @export
testRequiredFilesWarning <- function(files) {
  missing_files <- !file.exists(files)
  for (file in files[missing_files]) {
    logwarn("Missing required file: %s", file)
    #q()
  }
}

#' @export
testRequiredArguments <- function(arg_list, parameters) {
  for (par in parameters) {
    if (!par %in% names(arg_list)) {
      logerror("Missing required argument: %s", par)
      q()
    }
  }
}

#' @export
parseType<-function (arg_list,parameters,functType){
  for(key in parameters){
    arg_list[key]<-functType(arg_list[key])
  }
  return (arg_list)
}

#' @export
parseString<-function (arg_list,parameters){
  for(key in parameters){
    if(key %in% names(arg_list)) {
      arg_list[key] <- str_replace_all(arg_list[key], c("\"" = "", "\'" = ""))
    }
  }
  return (arg_list)
}

#' @export
parseList<-function (arg_list,parameters){
  for(key in parameters){
    items<-str_replace_all(as.character(arg_list[key])," ","")
    arg_list[key]<- base::strsplit(items,split=",")
  }
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
                                           , "peptideMissingValueImputation")

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
    # First ensure pacman is installed
    if (!requireNamespace("pacman", quietly = TRUE)) {
        if (verbose) message("Installing pacman...")
        install.packages("pacman")
    }
    required_packages <- c(
        # CRAN packages
        "tidyverse", "seqinr", "lazyeval", "rlang", "glue", "GGally",
        "here", "tibble", "mixOmics", "limma", "magrittr", "future.apply", 
        "tictoc", "beepr", "furrr", "readxl", "writexl", "RColorBrewer",
        "multidplyr", "RSpectra", "progress", "Rcpp", "RcppEigen",
        "qvalue", "Glimma", "ruv", "iq", "ggrepel", "patchwork",
        "dplyr", "gtools", "shiny", "DT", "gh",
        # Bioconductor packages
        "BiocManager",
        # GitHub packages
        "RUVIIIC", "ProteomeScholaR", "UniProt.ws"
    )
    library(pacman)

    # Install packages if missing
    if (!requireNamespace("RUVIIIC", quietly = TRUE)) {
        if (verbose) message("Installing RUVIIIC from GitHub...")
        devtools::install_github("cran/RUVIIIC")
    }
    if (!requireNamespace("ProteomeScholaR", quietly = TRUE)) {
        if (verbose) message("Installing ProteomeScholaR from GitHub...")
        devtools::install_github("APAF-BIOINFORMATICS/ProteomeScholaR", ref = "dev-jr")
    }
    if (verbose) message("Loading all required packages...")
    p_load(char = required_packages)
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
        results_summary_dir = "character"
    )
)

#' @title Setup Project Directories
#' @description Creates and manages project directories with version control
#' @param base_dir Base directory path (optional, defaults to here::here())
#' @return DirectoryManager object containing all directory paths
#' @export
setupDirectories <- function(base_dir = here::here()) {
    # Assign all directories to global environment
    assign("base_dir", base_dir, envir = .GlobalEnv)
    assign("results_dir", file.path(base_dir, "results", "proteomics"), envir = .GlobalEnv)
    assign("data_dir", file.path(base_dir, "data"), envir = .GlobalEnv)
    assign("source_dir", file.path(base_dir, "scripts", "proteomics"), envir = .GlobalEnv)
    assign("de_output_dir", file.path(results_dir, "de_proteins"), envir = .GlobalEnv)
    assign("publication_graphs_dir", file.path(results_dir, "publication_graphs"), envir = .GlobalEnv)
    assign("timestamp", format(Sys.time(), "%Y%m%d_%H%M%S"), envir = .GlobalEnv)
    assign("qc_dir", file.path(publication_graphs_dir, "filtering_qc"), envir = .GlobalEnv)
    assign("time_dir", file.path(qc_dir, timestamp), envir = .GlobalEnv)
    assign("results_summary_dir", file.path(base_dir, "results_summary", "proteomics"), envir = .GlobalEnv)
    
    # Directory management function with versioning
    manageDirectoryWithPrev <- function(dir_path) {
        if (dir.exists(dir_path)) {
            prev_dir <- paste0(dir_path, "_prev")
            if (dir.exists(prev_dir)) {
                unlink(prev_dir, recursive = TRUE)
            }
            if (dir.exists(dir_path) && length(list.files(dir_path)) > 0) {
                file.rename(dir_path, prev_dir)
            }
        }
        if (!dir.exists(dir_path)) {
            dir.create(dir_path, recursive = TRUE)
        }
    }
    
    # Simple directory creation without versioning
    createDirectory <- function(dir_path) {
        dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Create data and source directories (no versioning)
    c(data_dir, source_dir) |> 
        sapply(createDirectory)
    
    # Create results directories (with versioning)
    c(
        results_dir,
        de_output_dir,
        publication_graphs_dir,
        qc_dir,
        time_dir,
        results_summary_dir,
        file.path(results_dir, "clean_proteins"),
        file.path(results_dir, "protein_qc"),
        file.path(results_dir, "peptide_qc")
    ) |> 
        sapply(manageDirectoryWithPrev)
    
    # Return the DirectoryManager object
    return(new("DirectoryManager",
        base_dir = base_dir,
        results_dir = results_dir,
        data_dir = data_dir,
        source_dir = source_dir,
        de_output_dir = de_output_dir,
        publication_graphs_dir = publication_graphs_dir,
        timestamp = timestamp,
        qc_dir = qc_dir,
        time_dir = time_dir,
        results_summary_dir = results_summary_dir
    ))
}

showDirectories <- function() {
    # List of all directory variables we want to show
    dir_vars <- c(
        "base_dir", "results_dir", "data_dir", "source_dir", 
        "de_output_dir", "publication_graphs_dir", 
        "qc_dir", "time_dir", "results_summary_dir"
    )
    
    cat("Project Directory Structure:\n")
    cat("===========================\n\n")
    
    for (var_name in dir_vars) {
        # Get the path from global environment
        dir_path <- get(var_name, envir = .GlobalEnv)
        exists <- dir.exists(dir_path)
        prev_path <- paste0(dir_path, "_prev")
        has_prev <- dir.exists(prev_path)
        
        # Format the output
        status <- if (exists) "✓" else "✗"
        
        cat(sprintf("%-25s [%s] %s\n", 
            gsub("_dir$", "", var_name),  # Remove _dir suffix for display
            status,
            dir_path
        ))
        if (has_prev) {
            cat(sprintf("%25s %s\n", "", prev_path))
        }
    }
    
    cat("\nLegend: ✓ = exists, ✗ = missing\n")
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
#' @return Invisible NULL
#' @export
copyToResultsSummary <- function() {
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
            dir.create(
                file.path(results_summary_dir, subdir), 
                recursive = TRUE, 
                showWarnings = FALSE
            )
        })
    
    # Define files to copy with their display names
    files_to_copy <- list(
        # QC Figures
        list(
            source = file.path(time_dir, "correlation_filtered_combined_plots.png"),
            dest = "QC_figures",
            is_dir = FALSE,
            display_name = "Correlation Filtered Plots"
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
            source = file.path(de_output_dir, "de_proteins_long_annot.xlsx"),
            dest = "Publication_tables",
            is_dir = FALSE,
            display_name = "Proteomics Data Annotated",
            new_name = "Proteomics_data_annotated.xlsx"
        ),
        list(
            source = file.path(source_dir, "study_parameters.txt"),
            dest = "Study_report",
            is_dir = FALSE,
            display_name = "Study Parameters"
        )
    )
    
    cat("Copying files to Results Summary...\n")
    cat("===================================\n\n")
    
    # Copy and check each file/folder
    files_to_copy |> 
        lapply(\(file_spec) {
            dest_dir <- file.path(results_summary_dir, file_spec$dest)
            
            # Get initial status
            if (!is.null(file_spec$type) && file_spec$type == "object") {
                source_exists <- exists(file_spec$source, envir = .GlobalEnv)
            } else {
                source_exists <- if (file_spec$is_dir) {
                    dir.exists(file_spec$source)
                } else {
                    file.exists(file_spec$source)
                }
            }
            
            # Perform copy operation
            if (!is.null(file_spec$type) && file_spec$type == "object") {
                if (source_exists) {
                    obj <- get(file_spec$source, envir = .GlobalEnv)
                    write.table(
                        obj,
                        file = file.path(dest_dir, file_spec$save_as),
                        sep = "\t",
                        row.names = FALSE,
                        quote = FALSE
                    )
                }
                dest_path <- file.path(dest_dir, file_spec$save_as)
                dest_exists <- file.exists(dest_path)
            } else if (file_spec$is_dir) {
                source_path <- file_spec$source
                dest_path <- file.path(dest_dir, basename(source_path))
                
                if (source_exists) {
                    # Create destination directory
                    dir.create(dest_path, showWarnings = FALSE, recursive = TRUE)
                    
                    # Copy all files from source to destination
                    files_to_copy <- list.files(source_path, full.names = TRUE, recursive = TRUE)
                    
                    files_to_copy |> 
                        lapply(\(f) {
                            rel_path <- sub(paste0("^", source_path, "/"), "", f)
                            dest_file <- file.path(dest_path, rel_path)
                            dir.create(dirname(dest_file), showWarnings = FALSE, recursive = TRUE)
                            file.copy(f, dest_file, overwrite = TRUE)
                        })
                }
                dest_exists <- dir.exists(dest_path)
                
                # Get file counts for directories
                if (source_exists && dest_exists) {
                    source_files <- list.files(source_path, recursive = TRUE)
                    dest_files <- list.files(dest_path, recursive = TRUE)
                }
            } else {
                source_path <- file_spec$source
                dest_path <- file.path(dest_dir, 
                    if (!is.null(file_spec$new_name)) file_spec$new_name else basename(source_path))
                
                if (source_exists) {
                    file.copy(
                        from = source_path,
                        to = dest_path,
                        overwrite = TRUE
                    )
                }
                dest_exists <- file.exists(dest_path)
            }
            
            # Format and display status
            source_status <- if (source_exists) "✓" else "✗"
            dest_status <- if (dest_exists) "✓" else "✗"
            
            cat(sprintf("%-25s [%s → %s] %s\n",
                file_spec$display_name,
                source_status,
                dest_status,
                if (!is.null(file_spec$is_dir) && file_spec$is_dir) "Directory" else "File"
            ))
            
            if (!is.null(file_spec$is_dir) && file_spec$is_dir && source_exists && dest_exists) {
                cat(sprintf("%25s Files: %d → %d\n", "", 
                    length(source_files), 
                    length(dest_files)
                ))
            }
        })
    
    cat("\nLegend: ✓ = exists, ✗ = missing\n")
    cat("Arrow (→) shows source → destination status\n")
    
    invisible(NULL)
}