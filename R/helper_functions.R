# R/allGenerics.R
# This file contains all setGeneric definitions for the ProteomeScholaR package.
# Loading this file early ensures that generic functions are defined before
# their corresponding setMethod definitions are encountered during package loading.

# --- From proteinVsSamplesS4Objects.R ---

setGeneric(
  name = "setProteinData",
  def = function(theObject, protein_quant_table, protein_id_column) {
    standardGeneric("setProteinData")
  }
)

setGeneric(
  name = "cleanDesignMatrix",
  def = function(theObject) {
    standardGeneric("cleanDesignMatrix")
  }
)

setGeneric(
  name = "proteinIntensityFiltering",
  def = function(
    theObject,
    proteins_intensity_cutoff_percentile = NULL,
    proteins_proportion_of_samples_below_cutoff = NULL,
    core_utilisation = NULL
  ) {
    standardGeneric("proteinIntensityFiltering")
  },
  signature = c("theObject", "proteins_intensity_cutoff_percentile", "proteins_proportion_of_samples_below_cutoff", "core_utilisation")
)

setGeneric(
  name = "removeProteinsWithOnlyOneReplicate",
  def = function(theObject, core_utilisation = NULL, grouping_variable = NULL) {
    standardGeneric("removeProteinsWithOnlyOneReplicate")
  },
  signature = c("theObject", "core_utilisation", "grouping_variable")
)

setGeneric(
  name = "proteinMissingValueImputationLimpa",
  def = function(theObject,
                 dpc_results = NULL,
                 dpc_slope = 0.8,
                 quantified_protein_column = NULL,
                 verbose = TRUE,
                 chunk = 1000) {
    standardGeneric("proteinMissingValueImputationLimpa")
  },
  signature = c("theObject")
)

setGeneric(
  name = "plotRle",
  def = function(theObject, grouping_variable, yaxis_limit = c(), sample_label = NULL) {
    standardGeneric("plotRle")
  },
  signature = c("theObject", "grouping_variable", "yaxis_limit", "sample_label")
)

setGeneric(
  name = "plotRleList",
  def = function(theObject, list_of_columns, yaxis_limit = c()) {
    standardGeneric("plotRleList")
  },
  signature = c("theObject", "list_of_columns", "yaxis_limit")
)

setGeneric(
  name = "plotPca",
  def = function(theObject, grouping_variable, shape_variable = NULL, label_column, title, font_size, cv_percentile = 0.90) {
    standardGeneric("plotPca")
  },
  signature = c("theObject", "grouping_variable", "shape_variable", "label_column", "title", "font_size", "cv_percentile")
)

setGeneric(
  name = "plotPcaList",
  def = function(theObject, grouping_variables_list, label_column, title, font_size, cv_percentile = 0.90) {
    standardGeneric("plotPcaList")
  },
  signature = c("theObject", "grouping_variables_list", "label_column", "title", "font_size", "cv_percentile")
)

setGeneric(
  name = "getPcaMatrix",
  def = function(theObject) {
    standardGeneric("getPcaMatrix")
  },
  signature = c("theObject")
)

setGeneric(
  name = "proteinTechRepCorrelation",
  def = function(theObject, tech_rep_num_column = NULL, tech_rep_remove_regex = NULL) {
    standardGeneric("proteinTechRepCorrelation")
  },
  signature = c("theObject", "tech_rep_num_column", "tech_rep_remove_regex")
)

setGeneric(
  name = "plotPearson",
  def = function(theObject, tech_rep_remove_regex, correlation_group = NA) {
    standardGeneric("plotPearson")
  },
  signature = c("theObject", "tech_rep_remove_regex", "correlation_group")
)

setGeneric("InitialiseGrid", function(dummy = NULL) {
  standardGeneric("InitialiseGrid")
})

setGeneric(
  name = "createGridQC",
  def = function(theObject, pca_titles, density_titles, rle_titles, pearson_titles, save_path = NULL, file_name = "pca_density_rle_pearson_corr_plots_merged") {
    standardGeneric("createGridQC")
  },
  signature = c("theObject", "pca_titles", "density_titles", "rle_titles", "pearson_titles", "save_path", "file_name")
)

setGeneric(
  name = "normaliseBetweenSamples",
  def = function(theObject, normalisation_method = NULL) {
    standardGeneric("normaliseBetweenSamples")
  },
  signature = c("theObject", "normalisation_method")
)

setGeneric(
  name = "pearsonCorForSamplePairs",
  def = function(theObject, tech_rep_remove_regex = NULL, correlation_group = NA) {
    standardGeneric("pearsonCorForSamplePairs")
  },
  signature = c("theObject", "tech_rep_remove_regex", "correlation_group")
)

setGeneric(
  name = "getNegCtrlProtAnova",
  def = function(
    theObject,
    ruv_grouping_variable = NULL,
    percentage_as_neg_ctrl = NULL,
    num_neg_ctrl = NULL,
    ruv_qval_cutoff = NULL,
    ruv_fdr_method = NULL
  ) {
    standardGeneric("getNegCtrlProtAnova")
  },
  signature = c("theObject", "ruv_grouping_variable", "num_neg_ctrl", "ruv_qval_cutoff", "ruv_fdr_method")
)

setGeneric(
  name = "getLowCoefficientOfVariationProteins",
  def = function(
    theObject,
    percentage_as_neg_ctrl = NULL,
    num_neg_ctrl = NULL
  ) {
    standardGeneric("getLowCoefficientOfVariationProteins")
  },
  signature = c("theObject", "percentage_as_neg_ctrl", "num_neg_ctrl")
)

setGeneric(
  name = "ruvCancor",
  def = function(theObject, ctrl = NULL, num_components_to_impute = NULL, ruv_grouping_variable = NULL) {
    standardGeneric("ruvCancor")
  },
  signature = c("theObject", "ctrl", "num_components_to_impute", "ruv_grouping_variable")
)

setGeneric(
  name = "ruvCancorFast",
  def = function(theObject, ctrl = NULL, num_components_to_impute = NULL, ruv_grouping_variable = NULL, simple_imputation_method = "none") {
    standardGeneric("ruvCancorFast")
  },
  signature = c("theObject", "ctrl", "num_components_to_impute", "ruv_grouping_variable")
)

setGeneric(
  name = "getRuvIIIReplicateMatrix",
  def = function(theObject, ruv_grouping_variable = NULL) {
    standardGeneric("getRuvIIIReplicateMatrix")
  },
  signature = c("theObject", "ruv_grouping_variable")
)

setGeneric(
  name = "ruvIII_C_Varying",
  def = function(theObject, ruv_grouping_variable = NULL, ruv_number_k = NULL, ctrl = NULL) {
    standardGeneric("ruvIII_C_Varying")
  },
  signature = c("theObject", "ruv_grouping_variable", "ruv_number_k", "ctrl")
)

setGeneric(
  name = "removeRowsWithMissingValuesPercent",
  def = function(
    theObject,
    ruv_grouping_variable = NULL,
    groupwise_percentage_cutoff = NULL,
    max_groups_percentage_cutoff = NULL,
    proteins_intensity_cutoff_percentile = NULL
  ) {
    standardGeneric("removeRowsWithMissingValuesPercent")
  },
  signature = c(
    "theObject",
    "ruv_grouping_variable",
    "groupwise_percentage_cutoff",
    "max_groups_percentage_cutoff",
    "proteins_intensity_cutoff_percentile"
  )
)

setGeneric(
  name = "averageTechReps",
  def = function(theObject, design_matrix_columns, biological_replicate_column = NULL) {
    standardGeneric("averageTechReps")
  },
  signature = c("theObject", "design_matrix_columns", "biological_replicate_column")
)

setGeneric(
  name = "preservePeptideNaValues",
  def = function(peptide_obj, protein_obj) {
    standardGeneric("preservePeptideNaValues")
  },
  signature = c("peptide_obj", "protein_obj")
)

setGeneric(
  name = "chooseBestProteinAccession",
  def = function(theObject, delim = NULL, seqinr_obj = NULL, seqinr_accession_column = NULL, replace_zero_with_na = NULL, aggregation_method = NULL) {
    standardGeneric("chooseBestProteinAccession")
  },
  signature = c("theObject", "delim", "seqinr_obj", "seqinr_accession_column")
)

setGeneric(
  name = "chooseBestProteinAccessionSumDuplicates" # Note: This might be redundant or specific, verify usage
  , def = function(theObject, delim = NULL, seqinr_obj = NULL, seqinr_accession_column = NULL, replace_zero_with_na = NULL, aggregation_method = NULL) {
    standardGeneric("chooseBestProteinAccessionSumDuplicates")
  },
  signature = c("theObject", "delim", "seqinr_obj", "seqinr_accession_column")
)

setGeneric(
  name = "filterSamplesByProteinCorrelationThreshold",
  def = function(theObject, threshold = NULL, correlation_group = NULL, tech_rep_remove_regex = NULL) {
    standardGeneric("filterSamplesByProteinCorrelationThreshold")
  },
  signature = c("theObject", "threshold", "correlation_group", "tech_rep_remove_regex")
)

setGeneric(
  name = "plotDensity",
  def = function(theObject, grouping_variable, title = "", font_size = 8) {
    standardGeneric("plotDensity")
  },
  signature = c("theObject", "grouping_variable", "title", "font_size")
) # Base signature might need refinement based on methods

setGeneric(
  name = "plotDensityList",
  def = function(theObject, grouping_variables_list, title = "", font_size = 8) {
    standardGeneric("plotDensityList")
  },
  signature = c("theObject", "grouping_variables_list", "title", "font_size")
)


# --- From peptideVsSamplesS4Objects.R ---

setGeneric(
  name = "cleanDesignMatrixPeptide",
  def = function(theObject) {
    standardGeneric("cleanDesignMatrixPeptide")
  }
)

setGeneric(
  name = "srlQvalueProteotypicPeptideClean",
  def = function(theObject, qvalue_threshold = NULL, global_qvalue_threshold = NULL, choose_only_proteotypic_peptide = NULL, input_matrix_column_ids = NULL) {
    standardGeneric("srlQvalueProteotypicPeptideClean")
  },
  signature = c("theObject", "qvalue_threshold", "global_qvalue_threshold", "choose_only_proteotypic_peptide", "input_matrix_column_ids")
)

setGeneric(
  name = "rollUpPrecursorToPeptide",
  def = function(theObject, core_utilisation = NULL) {
    standardGeneric("rollUpPrecursorToPeptide")
  },
  signature = c("theObject", "core_utilisation")
)

setGeneric(
  name = "peptideIntensityFiltering",
  def = function(theObject, peptides_intensity_cutoff_percentile = NULL, peptides_proportion_of_samples_below_cutoff = NULL, core_utilisation = NULL) {
    standardGeneric("peptideIntensityFiltering")
  },
  signature = c("theObject", "peptides_intensity_cutoff_percentile", "peptides_proportion_of_samples_below_cutoff", "core_utilisation")
)

setGeneric(
  name = "removePeptidesWithMissingValuesPercent",
  def = function(
    theObject,
    grouping_variable = NULL,
    groupwise_percentage_cutoff = NULL,
    max_groups_percentage_cutoff = NULL,
    peptides_intensity_cutoff_percentile = NULL
  ) {
    standardGeneric("removePeptidesWithMissingValuesPercent")
  },
  signature = c(
    "theObject",
    "grouping_variable",
    "groupwise_percentage_cutoff",
    "max_groups_percentage_cutoff",
    "peptides_intensity_cutoff_percentile"
  )
)

setGeneric(
  name = "filterMinNumPeptidesPerProtein",
  def = function(theObject, ...) {
    standardGeneric("filterMinNumPeptidesPerProtein")
  },
  signature = c("theObject")
)

setGeneric(
  name = "filterMinNumPeptidesPerSample",
  def = function(theObject, peptides_per_sample_cutoff = NULL, core_utilisation = NULL, inclusion_list = NULL) {
    standardGeneric("filterMinNumPeptidesPerSample")
  },
  signature = c("theObject", "peptides_per_sample_cutoff", "core_utilisation", "inclusion_list")
)

setGeneric(
  name = "removePeptidesWithOnlyOneReplicate",
  def = function(theObject, replicate_group_column = NULL, core_utilisation = NULL) {
    standardGeneric("removePeptidesWithOnlyOneReplicate")
  },
  signature = c("theObject", "replicate_group_column", "core_utilisation")
)

setGeneric(
  name = "plotPeptidesProteinsCountsPerSample",
  def = function(theObject) {
    standardGeneric("plotPeptidesProteinsCountsPerSample")
  },
  signature = c("theObject")
)

setGeneric(
  name = "peptideMissingValueImputation",
  def = function(theObject, imputed_value_column = NULL, proportion_missing_values = NULL, core_utilisation = NULL) {
    standardGeneric("peptideMissingValueImputation")
  },
  signature = c("theObject", "imputed_value_column", "proportion_missing_values", "core_utilisation")
)


# --- From metabolite_qc.R ---

setGeneric(
  name = "metaboliteIntensityFiltering",
  def = function(theObject, metabolites_intensity_cutoff_percentile = NULL, metabolites_proportion_of_samples_below_cutoff = NULL) {
    standardGeneric("metaboliteIntensityFiltering")
  },
  signature = c("theObject")
) # Dispatch only on the object

setGeneric(
  name = "resolveDuplicateFeatures",
  def = function(theObject, itsd_pattern_columns = NULL) {
    standardGeneric("resolveDuplicateFeatures")
  }
)

# --- From metabolite_normalization.R ---

setGeneric(
  name = "logTransformAssays",
  def = function(theObject, offset = 1, ...) {
    standardGeneric("logTransformAssays")
  }
)

setGeneric(
  name = "normaliseUntransformedData",
  def = function(theObject, method = "ITSD", ...) {
    standardGeneric("normaliseUntransformedData")
  }
)

# --- From metaboliteVsSamplesS4Objects.R ---
# (Only generics defined uniquely in this file)

setGeneric(
  name = "getNegCtrlMetabAnova",
  def = function(theObject,
                 ruv_grouping_variable = NULL,
                 percentage_as_neg_ctrl = NULL, # Allow list/vector/single value
                 num_neg_ctrl = NULL,
                 ruv_qval_cutoff = NULL,
                 ruv_fdr_method = NULL) {
    standardGeneric("getNegCtrlMetabAnova")
  },
  signature = c("theObject")
) # Primary dispatch on object

setGeneric(
  name = "chooseBestProteinAccession",
  def = function(theObject, delim = NULL, seqinr_obj = NULL, seqinr_accession_column = NULL, replace_zero_with_na = NULL, aggregation_method = NULL) {
    standardGeneric("chooseBestProteinAccession")
  },
  signature = "theObject"
)
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
    keys <- as.character(as.vector(keys))
    attribute <- as.vector(attributes)

    hash <- new.env(hash = TRUE, parent = parent.frame())

    if (length(keys) != length(attributes)) {
        warning("Length of keys != Length of attributes list.")
        return(1)
    }

    mapply(
        function(k, a) base::assign(k, a, envir = hash),
        keys, attributes
    )

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
    if (base::exists(key, hash)) {
        return(base::get(key, hash))
    } else {
        return(NA)
    }
}

##################################################################################################################

#' @export
createDirectoryIfNotExists <- function(file_path, mode = "0777") {
    ### Create directory recursively if it doesn't exist
    if (!file.exists(file_path)) {
        dir.create(file_path, showWarnings = TRUE, recursive = TRUE, mode = mode)
    }
}

#' @export
createDirIfNotExists <- function(file_path, mode = "0777") {
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
    temp <- tempfile(fileext = ".R")
    knitr::purl(file, output = temp)

    if (skip_plots) {
        old_dev <- getOption("device")
        options(device = function(...) {
            .Call("R_GD_nullDevice", PACKAGE = "grDevices")
        })
    }
    source(temp)
    if (skip_plots) {
        options(device = old_dev)
    }
}


##################################################################################################################
# =====================================================================================================
#' @export
createOutputDir <- function(output_dir, no_backup) {
    if (output_dir == "") {
        logerror("output_dir is an empty string")
        q()
    }
    if (dir.exists(output_dir)) {
        if (no_backup) {
            unlink(output_dir, recursive = TRUE)
        } else {
            backup_name <- paste(output_dir, "_prev", sep = "")
            if (dir.exists(backup_name)) {
                unlink(backup_name, recursive = TRUE)
            }
            system(paste("mv", output_dir, backup_name))
        }
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
parseType <- function(arg_list, parameters, functType) {
    invisible(sapply(parameters, function(key) {
        arg_list[key] <- functType(arg_list[key])
    }))
    return(arg_list)
}

#' @export
parseString <- function(arg_list, parameters) {
    invisible(sapply(parameters, function(key) {
        if (key %in% names(arg_list)) {
            arg_list[key] <- str_replace_all(arg_list[key], c("\"" = "", "\'" = ""))
        }
    }))
    return(arg_list)
}

#' @export
parseList <- function(arg_list, parameters) {
    invisible(sapply(parameters, function(key) {
        items <- str_replace_all(as.character(arg_list[key]), " ", "")
        arg_list[key] <- base::strsplit(items, split = ",")
    }))
    return(arg_list)
}

#' @export
isArgumentDefined <- function(arg_list, parameter) {
    return(!is.null(arg_list[parameter]) & (parameter %in% names(arg_list)) & as.character(arg_list[parameter]) != "")
}


#' @export
setArgsDefault <- function(args, value_name, as_func, default_val = NA) {
    if (isArgumentDefined(args, value_name)) {
        args <- parseType(
            args,
            c(value_name),
            as_func
        )
    } else {
        logwarn(paste0(value_name, " is undefined, default value set to ", as.character(default_val), "."))
        args[[value_name]] <- default_val
    }

    return(args)
}

# =====================================================================================================


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
savePlot <- function(plot, base_path, plot_name, formats = c("pdf", "png", "jpg"), width = 7, height = 7, ...) {
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
                    ggsave(filename = file_path, plot = p, device = format, width = width, height = height, ...)
                })
            }
        })
    } else {
        # Single plot - original behavior
        purrr::walk(formats, \(format){
            file_path <- file.path(base_path, paste0(plot_name, ".", format))
            ggsave(filename = file_path, plot = plot, device = format, width = width, height = height, ...)
        })
    }
}


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
save_plot <- function(plot, base_path, plot_name, formats = c("pdf", "png", "jpg"), width = 7, height = 7, ...) {
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

    # Fix: Safely access nested list to avoid index errors
    object_value <- NULL
    if (!is.null(theObject@args) &&
        !is.null(theObject@args[[function_name]]) &&
        is.list(theObject@args[[function_name]])) {
        object_value <- theObject@args[[function_name]][[param_name_string]]
    }

    # print(paste0("param_value = ", param_value))

    error <- paste0(function_name, paste0(": '", param_name_string, "' is not defined.\n"))

    if (!is.null(param_value)) {
        return(param_value)
    } else if (!is.null(object_value)) {
        # print("use object value")
        return(object_value)
    } else if (!is.null(default_value)) {
        return(default_value)
    } else {
        stop(error)
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

    if (!is.null(param_value)) {
        return(param_value)
    } else if (!is.null(object_value)) {
        return(object_value)
    } else if (!is.null(default_value)) {
        return(default_value)
    } else {
        warning(error)
        return(NULL)
    }
}

#' Update the parameter in the object
#' @export
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
readConfigFile <- function(file = file.path(source_dir, "config.ini")) {
    config_list <- read.config(file = file, file.type = "ini")

    # to set the number of cores to be used in the parallel processing
    if ("globalParameters" %in% names(config_list)) {
        if ("number_of_cpus" %in% names(config_list[["globalParameters"]])) {
            print(paste0(
                "Read globalParameters: number_of_cpus = ",
                config_list$globalParameters$number_of_cpus
            ))
            core_utilisation <- new_cluster(config_list$globalParameters$number_of_cpus)
            cluster_library(core_utilisation, c("tidyverse", "glue", "rlang", "lazyeval"))

            list_of_multithreaded_functions <- c(
                "rollUpPrecursorToPeptide",
                "peptideIntensityFiltering",
                "filterMinNumPeptidesPerProtein",
                "filterMinNumPeptidesPerSample",
                "removePeptidesWithOnlyOneReplicate",
                "peptideMissingValueImputation",
                "removeProteinsWithOnlyOneReplicate"
            )

            setCoreUtilisation <- function(config_list, function_name) {
                if (!function_name %in% names(config_list)) {
                    config_list[[function_name]] <- list()
                }
                config_list[[function_name]][["core_utilisation"]] <- core_utilisation

                config_list
            }

            for (x in list_of_multithreaded_functions) {
                config_list <- setCoreUtilisation(config_list, x)
            }

            config_list[["globalParameters"]][["plots_format"]] <- str_split(config_list[["globalParameters"]][["plots_format"]], ",")[[1]]
        }
    }

    getConfigValue <- function(config_list, section, value) {
        config_list[[section]][[value]]
    }

    setConfigValueAsNumeric <- function(config_list, section, value) {
        config_list[[section]][[value]] <- as.numeric(config_list[[section]][[value]])
        config_list
    }

    if ("srlQvalueProteotypicPeptideClean" %in% names(config_list)) {
        config_list[["srlQvalueProteotypicPeptideClean"]][["input_matrix_column_ids"]] <- str_split(config_list[["srlQvalueProteotypicPeptideClean"]][["input_matrix_column_ids"]], ",")[[1]]

        print(paste0(
            "Read srlQvalueProteotypicPeptideClean: input_matrix_column_ids = ",
            paste0(config_list[["srlQvalueProteotypicPeptideClean"]][["input_matrix_column_ids"]],
                collapse = ", "
            )
        ))

        config_list <- setConfigValueAsNumeric(
            config_list,
            "srlQvalueProteotypicPeptideClean",
            "qvalue_threshold"
        )
        config_list <- setConfigValueAsNumeric(
            config_list,
            "srlQvalueProteotypicPeptideClean",
            "global_qvalue_threshold"
        )
        config_list <- setConfigValueAsNumeric(
            config_list,
            "srlQvalueProteotypicPeptideClean",
            "choose_only_proteotypic_peptide"
        )
    }


    if ("peptideIntensityFiltering" %in% names(config_list)) {
        config_list <- setConfigValueAsNumeric(
            config_list,
            "peptideIntensityFiltering",
            "peptides_intensity_cutoff_percentile"
        )
        config_list <- setConfigValueAsNumeric(
            config_list,
            "peptideIntensityFiltering",
            "peptides_proportion_of_samples_below_cutoff"
        )
    }


    if ("filterMinNumPeptidesPerProtein" %in% names(config_list)) {
        config_list <- setConfigValueAsNumeric(
            config_list,
            "filterMinNumPeptidesPerProtein",
            "peptides_per_protein_cutoff"
        )
        config_list <- setConfigValueAsNumeric(
            config_list,
            "filterMinNumPeptidesPerProtein",
            "peptidoforms_per_protein_cutoff"
        )
        # config_list <- setConfigValueAsNumeric(config_list
        #                                        , ""
        #                                        , "")
    }

    if ("filterMinNumPeptidesPerSample" %in% names(config_list)) {
        config_list <- setConfigValueAsNumeric(
            config_list,
            "filterMinNumPeptidesPerSample",
            "peptides_per_sample_cutoff"
        )

        if (!"inclusion_list" %in% names(config_list[["filterMinNumPeptidesPerSample"]])) {
            config_list[["filterMinNumPeptidesPerSample"]][["inclusion_list"]] <- ""
        }

        config_list[["filterMinNumPeptidesPerSample"]][["inclusion_list"]] <- str_split(config_list[["filterMinNumPeptidesPerSample"]][["inclusion_list"]], ",")[[1]]
    }

    if ("peptideMissingValueImputation" %in% names(config_list)) {
        config_list <- setConfigValueAsNumeric(
            config_list,
            "peptideMissingValueImputation",
            "proportion_missing_values"
        )
    }

    if ("removeRowsWithMissingValuesPercent" %in% names(config_list)) {
        config_list <- setConfigValueAsNumeric(
            config_list,
            "removeRowsWithMissingValuesPercent",
            "groupwise_percentage_cutoff"
        )

        config_list <- setConfigValueAsNumeric(
            config_list,
            "removeRowsWithMissingValuesPercent",
            "max_groups_percentage_cutoff"
        )

        config_list <- setConfigValueAsNumeric(
            config_list,
            "removeRowsWithMissingValuesPercent",
            "proteins_intensity_cutoff_percentile"
        )
    }


    if ("ruvIII_C_Varying" %in% names(config_list)) {
        config_list <- setConfigValueAsNumeric(
            config_list,
            "ruvIII_C_Varying",
            "ruv_number_k"
        )
    }

    if ("plotRle" %in% names(config_list)) {
        config_list[["plotRle"]][["yaxis_limit"]] <- str_split(config_list[["plotRle"]][["yaxis_limit"]], ",")[[1]] |>
            purrr::map_dbl(\(x) as.numeric(x))

        print(paste0(
            "Read plotRle: yaxis_limit = ",
            paste0(config_list[["plotRle"]][["yaxis_limit"]], collapse = ", ")
        ))
    }

    if ("deAnalysisParameters" %in% names(config_list)) {
        # Handle plots_format as array
        config_list[["deAnalysisParameters"]][["plots_format"]] <-
            str_split(config_list[["deAnalysisParameters"]][["plots_format"]], ",")[[1]]

        # Add new lfc_cutoff parameter
        config_list[["deAnalysisParameters"]][["lfc_cutoff"]] <- FALSE

        # Modify treat_lfc_cutoff to use ifelse
        config_list[["deAnalysisParameters"]][["treat_lfc_cutoff"]] <-
            ifelse(config_list[["deAnalysisParameters"]][["lfc_cutoff"]], log2(1.5), 0)

        # Handle args_group_pattern - remove quotes and fix escaping
        if ("args_group_pattern" %in% names(config_list[["deAnalysisParameters"]])) {
            config_list[["deAnalysisParameters"]][["args_group_pattern"]] <-
                gsub('^"|"$', "", config_list[["deAnalysisParameters"]][["args_group_pattern"]]) |>
                gsub(pattern = "\\\\", replacement = "\\")
        }

        # Convert numeric parameters
        config_list <- setConfigValueAsNumeric(
            config_list,
            "deAnalysisParameters",
            "de_q_val_thresh"
        )

        # Convert boolean parameters
        config_list[["deAnalysisParameters"]][["eBayes_trend"]] <-
            tolower(config_list[["deAnalysisParameters"]][["eBayes_trend"]]) == "true"
        config_list[["deAnalysisParameters"]][["eBayes_robust"]] <-
            tolower(config_list[["deAnalysisParameters"]][["eBayes_robust"]]) == "true"

        print(paste0(
            "Read deAnalysisParameters: formula_string = ",
            config_list[["deAnalysisParameters"]][["formula_string"]]
        ))
    }

    config_list
}


#' @export
#' @description Read the config file and specify the section and or parameter to update the object
#' @param theObject The object to be updated
#' @param file The file path to the config file
#' @param section The section to be updated
#' @param value The parameter value to be updated
readConfigFileSection <- function(
  theObject,
  file = file.path(source_dir, "config.ini"),
  function_name,
  parameter_name = NULL
) {
    config_list <- readConfigFile(
        file = file,
        file_type = "ini"
    )

    if (is.null(parameter_name)) {
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
        "git2r", "fs", "logger", "httr", "jsonlite",
        "configr", "webshot2", "shiny", "shinyjs", "shinyWidgets",
        "shinydashboard", "shinythemes", "shinycssloaders",

        # Added for BookChapter
        "conflicted", "tidytext",
        # Added from Suggests:
        "testthat", "ggplot2", "ggpubr", "svglite",
        "ggraph", "reticulate", "shinyFiles", "arrow"
    )

    bioc_packages <- c(
        "UniProt.ws", "mixOmics", "limma", "qvalue",
        "clusterProfiler", "GO.db", # GO.db is often a dependency, ensure it's listed
        "basilisk", "limpa"
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
            tryCatch(
                {
                    installer_func(pkg)
                    # After install, load it
                    pacman::p_load(char = pkg, character.only = TRUE)
                },
                error = function(e) {
                    warning(sprintf("Failed to install %s from %s: %s", pkg, source_name, e$message))
                }
            )
        } else {
            if (verbose) message(sprintf("%s is already installed, loading...", pkg))
            pacman::p_load(char = pkg, character.only = TRUE)
        }
    }

    # CRAN Packages
    if (verbose) message("\n--- Processing CRAN Packages ---")
    purrr::walk(cran_packages, ~ install_and_load(
        pkg = .x,
        # Use base R install.packages directly
        installer_func = function(p) utils::install.packages(p, dependencies = TRUE),
        source_name = "CRAN",
        verbose = verbose
    ))

    # Bioconductor Packages
    if (verbose) message("\n--- Processing Bioconductor Packages ---")
    purrr::walk(bioc_packages, ~ install_and_load(
        pkg = .x,
        installer_func = function(p) BiocManager::install(p, update = FALSE, ask = FALSE), # Use BiocManager
        source_name = "Bioconductor",
        verbose = verbose
    ))

    # GitHub Packages
    if (verbose) message("\n--- Processing GitHub Packages ---")
    purrr::iwalk(github_packages, ~ {
        pkg_name <- .y # Name of the package (e.g., "RUVIIIC")
        repo <- .x # Repository path (e.g., "cran/RUVIIIC")
        if (!requireNamespace(pkg_name, quietly = TRUE)) {
            if (verbose) message(sprintf("Installing %s from GitHub (%s)...", pkg_name, repo))
            tryCatch(
                {
                    # Force installation to handle potentially corrupt states
                    devtools::install_github(repo, force = TRUE)
                    pacman::p_load(char = pkg_name, character.only = TRUE)
                },
                error = function(e) {
                    warning(sprintf("Failed to install %s from GitHub (%s): %s", pkg_name, repo, e$message))
                }
            )
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
#' extract_experiment(x, mode = "range", start = 1, end = 3) # "20140602_ffs_expt1"
#' extract_experiment(x, mode = "start") # "20140602"
#' extract_experiment(x, mode = "end") # "junk"
#'
#' # Multiple strings
#' x <- c("20140602_ffs_expt1_r1_junk", "20140603_ffs_expt2_r2_test")
#' extract_experiment(x, mode = "range", start = 2, end = 3) # c("ffs_expt1", "ffs_expt2")
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
        } else if (mode == "start") {
            return(parts[1])
        } else if (mode == "end") {
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
            subdirs = c(
                "protein_qc", "peptide_qc", "clean_proteins", "de_proteins",
                "publication_graphs", "pathway_enrichment",
                file.path("publication_graphs", "filtering_qc")
            )
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
                dir.create,
                recursive = TRUE, showWarnings = FALSE
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
# S4 class WorkflowArgs removed - replaced with createStudyParametersFile function

# Old complex createWorkflowArgsFromConfig function removed - replaced with simple wrapper below

#' @title Format Configuration List
#' @param config_list List of configuration parameters
#' @param indent Number of spaces for indentation
#' @export
formatConfigList <- function(config_list, indent = 0) {
    message(sprintf("--- DEBUG66: Entering formatConfigList with %d items, indent=%d ---", length(config_list), indent))
    output <- character()

    # Exclude internal_workflow_source_dir from printing
    names_to_process <- names(config_list)
    names_to_process <- names_to_process[names_to_process != "internal_workflow_source_dir"]
    message(sprintf("   DEBUG66: Processing %d items after exclusions", length(names_to_process)))

    # FUNCTIONAL APPROACH - NO FOR LOOPS that cause hanging - Works in Shiny AND .rmd
    output <- if (requireNamespace("purrr", quietly = TRUE)) {
        purrr::map_chr(names_to_process, function(name) {
            message(sprintf("   DEBUG66: Processing config item '%s'", name))
            value <- config_list[[name]]
            message(sprintf("   DEBUG66: Item '%s' class: %s", name, paste(class(value), collapse = ", ")))

            # Skip core_utilisation and complex objects from display
            if (name == "core_utilisation" ||
                any(class(value) %in% c("process", "R6", "multidplyr_cluster", "cluster", "SOCKcluster"))) {
                message(sprintf("   DEBUG66: Skipping '%s' due to complex class", name))
                return("") # Return empty string instead of next
            }

            # Format the name
            name_formatted <- gsub("\\.", " ", name)
            name_formatted <- gsub("_", " ", name_formatted)
            name_formatted <- tools::toTitleCase(name_formatted)

            # Handle different value types
            if (is.list(value)) {
                message(sprintf("   DEBUG66: '%s' is a list with %d elements", name, length(value)))
                if (length(value) > 0 && !is.null(names(value))) {
                    output_lines <- paste0(paste(rep(" ", indent), collapse = ""), name_formatted, ":")
                    message(sprintf("   DEBUG66: Recursing into list '%s'", name))
                    recursive_result <- formatConfigList(value, indent + 2)
                    return(paste(c(output_lines, recursive_result), collapse = "\n"))
                } else if (length(value) > 0) { # Unnamed list, process elements functionally
                    message(sprintf("   DEBUG66: '%s' is unnamed list, processing elements", name))
                    header <- paste0(paste(rep(" ", indent), collapse = ""), name_formatted, ":")

                    # FUNCTIONAL processing of list elements - NO FOR LOOP
                    element_lines <- if (requireNamespace("purrr", quietly = TRUE)) {
                        purrr::imap_chr(value, function(item_val, item_idx) {
                            message(sprintf("   DEBUG66: Processing unnamed list item %d, class: %s", item_idx, paste(class(item_val), collapse = ", ")))
                            if (is.atomic(item_val) && length(item_val) == 1) {
                                paste0(paste(rep(" ", indent + 2), collapse = ""), "- ", as.character(item_val))
                            } else {
                                paste0(paste(rep(" ", indent + 2), collapse = ""), "- [Complex List Element]")
                            }
                        })
                    } else {
                        # Base R fallback for list processing
                        sapply(seq_along(value), function(item_idx) {
                            item_val <- value[[item_idx]]
                            if (is.atomic(item_val) && length(item_val) == 1) {
                                paste0(paste(rep(" ", indent + 2), collapse = ""), "- ", as.character(item_val))
                            } else {
                                paste0(paste(rep(" ", indent + 2), collapse = ""), "- [Complex List Element]")
                            }
                        })
                    }
                    return(paste(c(header, element_lines), collapse = "\n"))
                } else { # Empty list
                    message(sprintf("   DEBUG66: '%s' is empty list", name))
                    return(paste0(paste(rep(" ", indent), collapse = ""), name_formatted, ": [Empty List]"))
                }
            } else {
                message(sprintf("   DEBUG66: '%s' is atomic/non-list, attempting conversion", name))

                # SAFE value display formatting
                value_display <- tryCatch(
                    {
                        if (is.character(value) && length(value) > 5) {
                            paste(c(utils::head(value, 5), "..."), collapse = ", ")
                        } else if (is.character(value) && length(value) > 1) {
                            paste(value, collapse = ", ")
                        } else {
                            message(sprintf("   DEBUG66: About to convert '%s' to character", name))
                            result <- as.character(value)
                            message(sprintf("   DEBUG66: Conversion successful for '%s'", name))
                            if (length(result) == 0) "[Empty/NULL]" else result
                        }
                    },
                    error = function(e) {
                        message(sprintf("   DEBUG66: Error converting '%s': %s", name, e$message))
                        "[CONVERSION ERROR]"
                    }
                )

                return(paste0(paste(rep(" ", indent), collapse = ""), name_formatted, ": ", value_display))
            }
        }) |>
            (\(x) x[x != ""])() # Remove empty strings from skipped items
    } else {
        # Fallback using base R lapply if purrr not available
        result_list <- lapply(names_to_process, function(name) {
            value <- config_list[[name]]

            # Skip complex objects
            if (name == "core_utilisation" ||
                any(class(value) %in% c("process", "R6", "multidplyr_cluster", "cluster", "SOCKcluster"))) {
                return("")
            }

            # Format name
            name_formatted <- gsub("\\.", " ", name)
            name_formatted <- gsub("_", " ", name_formatted)
            name_formatted <- tools::toTitleCase(name_formatted)

            # Simple formatting for fallback
            value_display <- tryCatch(
                {
                    if (is.list(value)) {
                        "[List Object]"
                    } else {
                        as.character(value)[1] # Take first element only for safety
                    }
                },
                error = function(e) "[CONVERSION ERROR]"
            )

            paste0(paste(rep(" ", indent), collapse = ""), name_formatted, ": ", value_display)
        })
        unlist(result_list[result_list != ""])
    }

    # Flatten the nested strings
    output <- unlist(strsplit(output, "\n"))
    message(sprintf("--- DEBUG66: Exiting formatConfigList, returning %d lines ---", length(output)))
    return(output)
}

# S4 show method for WorkflowArgs removed - replaced with createStudyParametersFile function


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
    required_paths_in_current <- c(
        "results_dir", "results_summary_dir", "publication_graphs_dir",
        "time_dir", "qc_dir", "de_output_dir", "pathway_dir", "source_dir", "feature_qc_dir"
    )
    if (!is.list(current_paths) || !all(required_paths_in_current %in% names(current_paths))) {
        missing_req <- setdiff(required_paths_in_current, names(current_paths))
        rlang::abort(paste0("Essential paths missing from project_dirs for key ", sQuote(current_omic_key), ": ", paste(missing_req, collapse = ", ")))
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
    if (!is.null(current_paths$subfeature_qc_dir)) cat(sprintf("Sub-feature QC Dir: %s\n", current_paths$subfeature_qc_dir))
    cat("\n")

    # ROBUST: Try to get contrasts_tbl and design_matrix from environment first, then from files
    cat("Checking for required objects...\n")

    # Try contrasts_tbl from environment first
    if (is.null(contrasts_tbl)) {
        if (exists("contrasts_tbl", envir = parent.frame())) {
            contrasts_tbl <- get("contrasts_tbl", envir = parent.frame())
            cat("✓ Using 'contrasts_tbl' from calling environment\n")
        } else if (exists("contrasts_tbl", envir = .GlobalEnv)) {
            contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
            cat("✓ Using 'contrasts_tbl' from global environment\n")
        } else {
            # Fallback: try to load from file
            contrasts_file_options <- c(
                file.path(current_paths$source_dir, "contrasts_tbl.tab"),
                file.path(current_paths$source_dir, "contrast_strings.tab"),
                file.path(current_paths$source_dir, "contrasts.tab")
            )

            contrasts_loaded <- FALSE
            for (contrasts_file in contrasts_file_options) {
                if (file.exists(contrasts_file)) {
                    tryCatch(
                        {
                            contrasts_tbl <- readr::read_tsv(contrasts_file, show_col_types = FALSE)
                            cat(sprintf("✓ Loaded 'contrasts_tbl' from file: %s\n", basename(contrasts_file)))
                            contrasts_loaded <- TRUE
                            break
                        },
                        error = function(e) {
                            cat(sprintf("⚠ Failed to load contrasts from %s: %s\n", basename(contrasts_file), e$message))
                        }
                    )
                }
            }

            if (!contrasts_loaded) {
                cat("✗ 'contrasts_tbl' not found in environment or files\n")
            }
        }
    } else {
        cat("✓ Using provided 'contrasts_tbl' parameter\n")
    }

    # Try design_matrix from environment first
    if (is.null(design_matrix)) {
        if (exists("design_matrix", envir = parent.frame())) {
            design_matrix <- get("design_matrix", envir = parent.frame())
            cat("✓ Using 'design_matrix' from calling environment\n")
        } else if (exists("design_matrix", envir = .GlobalEnv)) {
            design_matrix <- get("design_matrix", envir = .GlobalEnv)
            cat("✓ Using 'design_matrix' from global environment\n")
        } else {
            # Fallback: try to load from file
            design_matrix_file <- file.path(current_paths$source_dir, "design_matrix.tab")

            if (file.exists(design_matrix_file)) {
                tryCatch(
                    {
                        design_matrix <- readr::read_tsv(design_matrix_file, show_col_types = FALSE)
                        cat(sprintf("✓ Loaded 'design_matrix' from file: %s\n", basename(design_matrix_file)))
                    },
                    error = function(e) {
                        cat(sprintf("✗ Failed to load design_matrix from %s: %s\n", basename(design_matrix_file), e$message))
                        design_matrix <- NULL
                    }
                )
            } else {
                cat(sprintf("✗ 'design_matrix' not found in environment or at expected file location: %s\n", design_matrix_file))
            }
        }
    } else {
        cat("✓ Using provided 'design_matrix' parameter\n")
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
        tryCatch(
            {
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
            },
            error = function(e) {
                warning(sprintf("Failed to save/copy Rmd file: %s", e$message))
                failed_copies[[length(failed_copies) + 1]] <- list(type = "rmd_copy", source = current_rmd, destination = current_paths$source_dir, display_name = "Current Rmd File", error = e$message)
            }
        )
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
            return(invisible(list(status = "cancelled", omic_key = current_omic_key, message = paste0("Backup and overwrite for ", current_omic_key, " cancelled by user."))))
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

            backup_has_content <- length(list.files(backup_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE)) > 0

            if (backup_copy_successful && (backup_has_content || length(contents_of_summary_dir) == 0)) {
                logger::log_info("Successfully backed up content of {current_paths$results_summary_dir} to: {backup_dir}")
                backup_info <- data.frame(original_dir = current_paths$results_summary_dir, backup_time = Sys.time(), omic_key = current_omic_key, stringsAsFactors = FALSE)
                tryCatch(
                    write.table(backup_info, file = file.path(backup_dir, "backup_info.txt"), sep = "\\t", row.names = FALSE, quote = FALSE),
                    error = function(e) logger::log_warn("Failed to write backup_info.txt: {e$message}")
                )

                logger::log_info("Clearing original results_summary_dir: {current_paths$results_summary_dir}")
                unlink_success <- tryCatch(
                    {
                        unlink(current_paths$results_summary_dir, recursive = TRUE, force = TRUE)
                        TRUE # Return TRUE on success
                    },
                    error = function(e) {
                        logger::log_error("Error unlinking {current_paths$results_summary_dir}: {e$message}")
                        FALSE # Return FALSE on error
                    }
                )

                if (unlink_success) {
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
            if (!dir.create(current_paths$results_summary_dir, recursive = TRUE, showWarnings = FALSE)) {
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
    de_files <- list.files(path = current_paths$de_output_dir, pattern = paste0("de_.+_long_annot\\.xlsx$"), full.names = TRUE) # Changed from \\w+ to .+ to allow hyphens

    purrr::imap(de_files, \(file, idx) {
        sheet_name <- sprintf("DE_Sheet%d", idx)
        base_name <- basename(file) |>
            stringr::str_remove("^de_") |>
            stringr::str_remove("_long_annot\\.xlsx$")
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
    tryCatch(
        {
            openxlsx::saveWorkbook(de_wb, de_results_excel_path, overwrite = TRUE)
            cat(paste("Successfully saved DE results to:", de_results_excel_path, "\n"))
        },
        error = function(e) {
            failed_copies[[length(failed_copies) + 1]] <- list(type = "workbook_save", path = de_results_excel_path, error = e$message)
        }
    )

    tryCatch(
        {
            openxlsx::saveWorkbook(enrichment_wb, enrichment_excel_path, overwrite = TRUE)
            cat(paste("Successfully saved Enrichment results to:", enrichment_excel_path, "\n"))
        },
        error = function(e) {
            failed_copies[[length(failed_copies) + 1]] <- list(type = "workbook_save", path = enrichment_excel_path, error = e$message)
        }
    )

    cat("\nCopying individual files/folders to Results Summary for ", current_omic_key, "...\n")
    cat("==============================================================\n\n")

    files_to_copy |>
        lapply(\(file_spec) {
            dest_dir_final <- file.path(current_paths$results_summary_dir, file_spec$dest)
            copy_success <- TRUE
            error_msg <- NULL
            source_display <- file_spec$source # For display purposes

            if (!is.null(file_spec$type) && file_spec$type == "object") {
                # ENHANCED: Use robust object sourcing - check environment first, then file
                obj <- NULL
                source_exists <- FALSE

                # Try to get object from environments first
                if (!is.null(get0(file_spec$source, envir = parent.frame()))) {
                    obj <- get(file_spec$source, envir = parent.frame())
                    source_exists <- TRUE
                    cat(sprintf("    ✓ Found '%s' in calling environment\n", file_spec$source))
                } else if (!is.null(get0(file_spec$source, envir = .GlobalEnv))) {
                    obj <- get(file_spec$source, envir = .GlobalEnv)
                    source_exists <- TRUE
                    cat(sprintf("    ✓ Found '%s' in global environment\n", file_spec$source))
                } else {
                    # Fallback: try to load from file
                    if (file_spec$source == "design_matrix") {
                        design_matrix_file <- file.path(current_paths$source_dir, "design_matrix.tab")
                        if (file.exists(design_matrix_file)) {
                            tryCatch(
                                {
                                    obj <- readr::read_tsv(design_matrix_file, show_col_types = FALSE)
                                    source_exists <- TRUE
                                    cat(sprintf("    ✓ Loaded '%s' from file: %s\n", file_spec$source, basename(design_matrix_file)))
                                },
                                error = function(e) {
                                    error_msg <- sprintf("Failed to load %s from file %s: %s", file_spec$source, basename(design_matrix_file), e$message)
                                    cat(sprintf("    ✗ %s\n", error_msg))
                                }
                            )
                        } else {
                            error_msg <- sprintf("Object '%s' not found in environment and file not found: %s", file_spec$source, design_matrix_file)
                        }
                    } else if (file_spec$source == "contrasts_tbl") {
                        contrasts_file_options <- c(
                            file.path(current_paths$source_dir, "contrasts_tbl.tab"),
                            file.path(current_paths$source_dir, "contrast_strings.tab"),
                            file.path(current_paths$source_dir, "contrasts.tab")
                        )

                        for (contrasts_file in contrasts_file_options) {
                            if (file.exists(contrasts_file)) {
                                tryCatch(
                                    {
                                        obj <- readr::read_tsv(contrasts_file, show_col_types = FALSE)
                                        source_exists <- TRUE
                                        cat(sprintf("    ✓ Loaded '%s' from file: %s\n", file_spec$source, basename(contrasts_file)))
                                        break
                                    },
                                    error = function(e) {
                                        cat(sprintf("    ⚠ Failed to load contrasts from %s: %s\n", basename(contrasts_file), e$message))
                                    }
                                )
                            }
                        }

                        if (!source_exists) {
                            error_msg <- sprintf("Object '%s' not found in environment and no readable contrasts files found", file_spec$source)
                        }
                    } else {
                        # For other objects, keep original behavior
                        error_msg <- sprintf("Object '%s' not found in parent/global environment", file_spec$source)
                    }
                }

                if (!source_exists) {
                    if (is.null(error_msg)) {
                        error_msg <- sprintf("Object '%s' not found in environment or files", file_spec$source)
                    }
                }
            } else {
                source_exists <- if (file_spec$is_dir) dir.exists(file_spec$source) else file.exists(file_spec$source)
                if (!source_exists) error_msg <- sprintf("Source %s not found: %s", if (file_spec$is_dir) "directory" else "file", file_spec$source)
            }

            if (source_exists) {
                dir.create(dest_dir_final, recursive = TRUE, showWarnings = FALSE)
                if (!is.null(file_spec$type) && file_spec$type == "object") {
                    tryCatch(
                        {
                            # Use the already-loaded object instead of getting from environment again
                            if (is.null(obj)) {
                                obj <- get(file_spec$source, envir = parent.frame()) # Fallback for other objects
                            }
                            dest_path <- file.path(dest_dir_final, file_spec$save_as)
                            write.table(obj, file = dest_path, sep = "\t", row.names = FALSE, quote = FALSE)
                            if (!file.exists(dest_path) || (file.exists(dest_path) && file.size(dest_path) == 0 && nrow(obj) > 0)) {
                                copy_success <- FALSE
                                error_msg <- "Failed to write object or file is empty"
                            }
                        },
                        error = function(e) {
                            copy_success <<- FALSE
                            error_msg <<- sprintf("Error writing object: %s", e$message)
                        }
                    )
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
                failed_copies[[length(failed_copies) + 1]] <- list(type = if (!is.null(file_spec$type) && file_spec$type == "object") "object" else if (file_spec$is_dir) "directory" else "file", source = source_display, destination = dest_dir_final, display_name = file_spec$display_name, error = error_msg)
            }
            cat(sprintf("%-35s [%s -> %s] %s\n", file_spec$display_name, if (source_exists) "✓" else "✗", if (copy_success && source_exists) "✓" else "✗", if (!is.null(file_spec$type) && file_spec$type == "object") "Object" else if (file_spec$is_dir) "Directory" else "File"))
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


#' Update Missing Value Parameters in Configuration List and S4 Object
#'
#' @description
#' Automatically calculates and updates the missing value filtering parameters in the configuration list
#' and S4 object @args based on the experimental design matrix. The function ensures at least a specified
#' number of groups have sufficient quantifiable values for analysis.
#'
#' @param theObject An S4 object containing the experimental data and design matrix.
#' @param min_reps_per_group Integer specifying the minimum number of replicates required in each passing group.
#'                          If a group has fewer total replicates than this value, the minimum is adjusted.
#' @param min_groups Integer specifying the minimum number of groups required to have sufficient
#'                  quantifiable values. Default is 2.
#' @param config_list_name The name of the global config list variable (defaults to "config_list").
#' @param env The environment where the global config list resides (defaults to .GlobalEnv).
#'
#' @return Updated S4 object with synchronized @args and global config_list
#'
#' @details
#' The function calculates:
#' - groupwise_percentage_cutoff: Based on minimum required replicates per group
#' - max_groups_percentage_cutoff: Based on minimum required groups
#'
#' Both the S4 object's @args slot and the global config_list are updated to maintain synchronization.
#'
#' @examples
#' \dontrun{
#' protein_log2_quant_cln <- updateMissingValueParameters(
#'     theObject = protein_log2_quant_cln,
#'     min_reps_per_group = 2,
#'     min_groups = 2
#' )
#' }
#'
#' @export
updateMissingValueParameters <- function(theObject,
                                         min_reps_per_group = 2,
                                         min_groups = 2,
                                         config_list_name = "config_list",
                                         env = .GlobalEnv) {
    # --- Input Validation ---
    if (!isS4(theObject)) {
        stop("'theObject' must be an S4 object.")
    }
    if (!"design_matrix" %in% methods::slotNames(theObject)) {
        stop("'theObject' must have a '@design_matrix' slot.")
    }
    if (!"args" %in% methods::slotNames(theObject)) {
        stop("'theObject' must have an '@args' slot.")
    }
    if (!exists(config_list_name, envir = env)) {
        stop("Global config list '", config_list_name, "' not found in the specified environment.")
    }

    # Extract design matrix from S4 object
    design_matrix <- theObject@design_matrix

    # Retrieve the global config list
    config_list <- get(config_list_name, envir = env)

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
    # Take the maximum percentage to ensure all groups meet minimum requirements
    groupwise_cutoff <- max(group_thresholds$missing_percent)

    # Calculate maximum failing groups allowed
    max_failing_groups <- total_groups - min_groups
    max_groups_cutoff <- round((max_failing_groups / total_groups) * 100, 3)

    # Update global config_list
    if (is.null(config_list$removeRowsWithMissingValuesPercent)) {
        config_list$removeRowsWithMissingValuesPercent <- list()
    }
    config_list$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff <- groupwise_cutoff
    config_list$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff <- max_groups_cutoff
    config_list$removeRowsWithMissingValuesPercent$proteins_intensity_cutoff_percentile <- 1

    # Assign updated config_list back to global environment
    assign(config_list_name, config_list, envir = env)

    # Update S4 object @args to match config_list
    if (is.null(theObject@args)) {
        theObject@args <- list()
    }
    if (is.null(theObject@args$removeRowsWithMissingValuesPercent)) {
        theObject@args$removeRowsWithMissingValuesPercent <- list()
    }

    theObject@args$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff <- groupwise_cutoff
    theObject@args$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff <- max_groups_cutoff
    theObject@args$removeRowsWithMissingValuesPercent$proteins_intensity_cutoff_percentile <- 1

    # Create informative message
    basic_msg <- sprintf(
        "Updated missing value parameters in both config_list and S4 object @args:
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
            detail = sprintf(
                "    Group %s: %d out of %d replicates required (%.3f%% missing allowed)",
                group, adjusted_min_reps, n_reps, missing_percent
            )
        ) |>
        dplyr::pull(detail)

    group_details <- paste(group_detail_strings, collapse = "\n")

    # Print the message
    message(paste(basic_msg, "\n\nGroup details:", group_details, sep = "\n"))
    message("✅ S4 object @args and global config_list are now synchronized")

    return(theObject)
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
        "https://raw.githubusercontent.com/APAF-bioinformatics/MultiScholaR/refs/heads/BookChapter/Workbooks/%s/report/%s",
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
                         rmd_filename = NULL,
                         project_dirs_object_name = "project_dirs",
                         output_format = NULL) {
    # --- Validate Inputs ---
    if (missing(omic_type) || !is.character(omic_type) || length(omic_type) != 1 || omic_type == "") {
        rlang::abort("`omic_type` must be a single non-empty character string.")
    }
    if (missing(experiment_label) || !is.character(experiment_label) || length(experiment_label) != 1 || experiment_label == "") {
        rlang::abort("`experiment_label` must be a single non-empty character string.")
    }
    if (!is.null(rmd_filename) && (!is.character(rmd_filename) || length(rmd_filename) != 1 || rmd_filename == "")) {
        rlang::abort("`rmd_filename` must be a single non-empty character string or NULL.")
    }

    # --- Retrieve Paths from Global Project Directories Object ---
    if (!exists(project_dirs_object_name, envir = .GlobalEnv)) {
        rlang::abort(paste0("Global object ", sQuote(project_dirs_object_name), " not found. Run setupDirectories() first."))
    }

    project_dirs_global <- get(project_dirs_object_name, envir = .GlobalEnv)
    current_omic_key <- paste0(omic_type, "_", experiment_label)

    if (!current_omic_key %in% names(project_dirs_global)) {
        rlang::abort(paste0("Key ", sQuote(current_omic_key), " not found in ", sQuote(project_dirs_object_name), ". Check omic_type and experiment_label."))
    }
    current_paths <- project_dirs_global[[current_omic_key]] # This contains base_dir, results_summary_dir etc. for the *labelled* omic

    if (!is.list(current_paths) ||
        is.null(current_paths$base_dir) || # Need base_dir to find the template Rmd
        is.null(current_paths$results_summary_dir)) {
        rlang::abort(paste0(
            "Essential paths (base_dir, results_summary_dir) missing for key ",
            sQuote(current_omic_key), " in ", sQuote(project_dirs_object_name), "."
        ))
    }

    # --- Read study_parameters.txt for workflow determination and render parameters ---
    params_path <- file.path(current_paths$source_dir, "study_parameters.txt")
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

    # --- Determine the Rmd template filename based on workflow ---
    if (is.null(rmd_filename)) {
        if (workflow_name == "DIA_limpa") {
            rmd_filename <- "DIANN_limpa_report.rmd"
        } else if (workflow_name == "DIA") {
            # TODO: Create DIA-specific report template if needed
            rmd_filename <- "DIANN_report.rmd" # Using standard template for now
        } else if (workflow_name == "TMT") {
            # TODO: Create TMT-specific report template if needed
            rmd_filename <- "TMT_report.rmd" # Assuming TMT template exists or will be created
        } else {
            rmd_filename <- "DIANN_report.rmd" # Default template
        }
    }

    # --- Locate the Rmd template file ---
    # Templates are part of the MultiScholaR package installation
    # First try package location, then download from GitHub if not found

    logger::log_info("Looking for template '{rmd_filename}' for omic type '{omic_type}'")

    # Try package installation location
    rmd_input_path <- system.file("Workbooks", omic_type, "report", rmd_filename,
        package = "MultiScholaR"
    )

    if (file.exists(rmd_input_path) && nzchar(rmd_input_path)) {
        logger::log_info("Found template in package installation: {rmd_input_path}")
    } else {
        # Download from GitHub as fallback
        logger::log_info("Template not found in package installation, attempting download from GitHub...")
        rmd_input_path <- downloadReportTemplate(omic_type, rmd_filename)
    }

    # Final validation
    if (!file.exists(rmd_input_path) || !nzchar(rmd_input_path)) {
        rlang::abort(paste0(
            "Unable to locate or download R Markdown template: ", sQuote(rmd_filename), "\n",
            "Omic type: ", omic_type, "\n",
            "Template should exist in package or be downloadable from:\n",
            "https://github.com/APAF-bioinformatics/MultiScholaR/tree/main/Workbooks/",
            omic_type, "/report/"
        ))
    }

    # --- Copy template to temporary location in project ---
    # This ensures knitr creates correct relative paths from project location
    # instead of from the package directory
    # Place temp file in BASE directory so relative paths work correctly
    temp_rmd_dir <- file.path(current_paths$base_dir, ".temp_reports")
    dir.create(temp_rmd_dir, recursive = TRUE, showWarnings = FALSE)

    temp_rmd_path <- file.path(temp_rmd_dir, basename(rmd_filename))
    file.copy(rmd_input_path, temp_rmd_path, overwrite = TRUE)

    logger::log_info("Copied template to temporary location: {temp_rmd_path}")

    # --- Construct Output Path (in the labelled results_summary directory) ---
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

    logger::log_info("Attempting to render report:")
    logger::log_info("- Rmd Source (Template): {rmd_input_path}")
    logger::log_info("- Rmd Temporary Copy: {temp_rmd_path}")
    logger::log_info("- Output File: {output_file_path}")
    logger::log_info("- Params: omic_type=\'{omic_type}\', experiment_label=\'{experiment_label}\'")


    # --- Render the Report ---
    rendered_path <- tryCatch(
        {
            result <- rmarkdown::render(
                input = temp_rmd_path, # Use temporary copy in base directory
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

            # Clean up temporary file after successful render
            if (file.exists(temp_rmd_path)) {
                unlink(temp_rmd_path)
                logger::log_info("Cleaned up temporary template file")
            }

            result
        },
        error = function(e) {
            # Clean up temporary file even on error
            if (file.exists(temp_rmd_path)) {
                unlink(temp_rmd_path)
                logger::log_info("Cleaned up temporary template file after error")
            }

            logger::log_error("Failed to render R Markdown report: {e$message}")
            logger::log_error("Input path: {rmd_input_path}")
            logger::log_error("Temporary path: {temp_rmd_path}")
            logger::log_error("Output path: {output_file_path}")
            NULL # Return NULL on failure
        }
    )

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
#'     theObject = myPeptideData,
#'     function_name = "peptideIntensityFiltering",
#'     parameter_name = "peptides_proportion_of_samples_below_cutoff",
#'     new_value = 0.7
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

#' Choose Best Protein Accession for Data Frame (S3 Version)
#'
#' @description A simplified S3 version of chooseBestProteinAccession that works
#' directly on data frames. Resolves protein groups by selecting the best
#' representative protein based on FASTA sequence data and aggregates
#' quantitative values.
#'
#' @param data_tbl Data frame containing protein quantitative data
#' @param protein_id_column Character string specifying the protein ID column name
#' @param seqinr_obj Data frame with FASTA sequence information (aa_seq_tbl_final)
#' @param seqinr_accession_column Character string specifying the accession column in seqinr_obj
#' @param delim Character string delimiter used to separate protein groups (default: ";")
#' @param replace_zero_with_na Logical, whether to replace zeros with NA (default: TRUE)
#' @param aggregation_method Character string aggregation method ("mean", "median", "sum") (default: "mean")
#'
#' @return Data frame with cleaned protein IDs and aggregated values
#'
#' @details
#' This function processes protein groups (e.g., "P12345;P67890;Q11111") by:
#' \itemize{
#'   \item Splitting groups using the specified delimiter
#'   \item Finding the best representative protein using FASTA sequence data
#'   \item Aggregating quantitative values for the chosen protein
#'   \item Returning cleaned data with single protein IDs
#' }
#'
#' The selection of "best" protein is based on:
#' \itemize{
#'   \item Presence in FASTA sequence database
#'   \item Protein length (longer sequences preferred)
#'   \item Alphabetical order as tiebreaker
#' }
#'
#' @examples
#' \dontrun{
#' cleaned_data <- chooseBestProteinAccession_s3(
#'     data_tbl = protein_data,
#'     protein_id_column = "Protein.Group",
#'     seqinr_obj = aa_seq_tbl_final,
#'     seqinr_accession_column = "uniprot_acc"
#' )
#' }
#'
#' @export
chooseBestProteinAccession_s3 <- function(data_tbl,
                                          protein_id_column,
                                          seqinr_obj,
                                          seqinr_accession_column = "uniprot_acc",
                                          delim = ";",
                                          replace_zero_with_na = TRUE,
                                          aggregation_method = "mean") {
    cat("=== Starting chooseBestProteinAccession_s3 ===\n")

    tryCatch(
        {
            # Validate inputs with detailed error messages
            cat("*** S3 CLEANUP: Validating inputs ***\n")

            if (is.null(data_tbl) || !is.data.frame(data_tbl)) {
                stop("data_tbl must be a non-null data frame")
            }

            if (nrow(data_tbl) == 0) {
                stop("data_tbl cannot be empty")
            }

            cat(sprintf("*** S3 CLEANUP: data_tbl dimensions: %d rows x %d cols ***\n", nrow(data_tbl), ncol(data_tbl)))
            cat(sprintf("*** S3 CLEANUP: data_tbl column names: %s ***\n", paste(head(names(data_tbl), 10), collapse = ", ")))

            if (!protein_id_column %in% names(data_tbl)) {
                stop(sprintf(
                    "Protein ID column '%s' not found in data. Available columns: %s",
                    protein_id_column, paste(names(data_tbl), collapse = ", ")
                ))
            }

            if (is.null(seqinr_obj) || !is.data.frame(seqinr_obj)) {
                stop("seqinr_obj must be a non-null data frame")
            }

            if (nrow(seqinr_obj) == 0) {
                stop("seqinr_obj cannot be empty")
            }

            cat(sprintf("*** S3 CLEANUP: seqinr_obj dimensions: %d rows x %d cols ***\n", nrow(seqinr_obj), ncol(seqinr_obj)))
            cat(sprintf("*** S3 CLEANUP: seqinr_obj column names: %s ***\n", paste(head(names(seqinr_obj), 10), collapse = ", ")))

            if (!seqinr_accession_column %in% names(seqinr_obj)) {
                stop(sprintf(
                    "Accession column '%s' not found in sequence data. Available columns: %s",
                    seqinr_accession_column, paste(names(seqinr_obj), collapse = ", ")
                ))
            }

            cat("*** S3 CLEANUP: Input validation passed ***\n")

            # Get non-quantitative columns (metadata columns) with error handling
            cat("*** S3 CLEANUP: Identifying column types ***\n")

            numeric_cols <- tryCatch(
                {
                    sapply(data_tbl, is.numeric)
                },
                error = function(e) {
                    stop(sprintf("Error checking column types: %s", e$message))
                }
            )

            quant_cols <- names(data_tbl)[numeric_cols]
            meta_cols <- names(data_tbl)[!numeric_cols]

            cat(sprintf(
                "*** S3 CLEANUP: Found %d quantitative columns and %d metadata columns ***\n",
                length(quant_cols), length(meta_cols)
            ))
            cat(sprintf("*** S3 CLEANUP: Quantitative columns: %s ***\n", paste(head(quant_cols, 5), collapse = ", ")))
            cat(sprintf("*** S3 CLEANUP: Metadata columns: %s ***\n", paste(head(meta_cols, 5), collapse = ", ")))

            # Replace zeros with NA if requested
            if (replace_zero_with_na && length(quant_cols) > 0) {
                cat("*** S3 CLEANUP: Replacing zeros with NA in quantitative columns ***\n")

                tryCatch(
                    {
                        data_tbl[quant_cols] <- lapply(data_tbl[quant_cols], function(x) {
                            if (is.numeric(x)) {
                                x[x == 0] <- NA
                            }
                            x
                        })
                        cat("*** S3 CLEANUP: Zero replacement completed ***\n")
                    },
                    error = function(e) {
                        stop(sprintf("Error replacing zeros with NA: %s", e$message))
                    }
                )
            }

            # Get unique protein groups
            cat("*** S3 CLEANUP: Getting unique protein groups ***\n")

            protein_groups <- tryCatch(
                {
                    unique(data_tbl[[protein_id_column]])
                },
                error = function(e) {
                    stop(sprintf("Error getting unique protein groups: %s", e$message))
                }
            )

            cat(sprintf("*** S3 CLEANUP: Processing %d unique protein groups ***\n", length(protein_groups)))
            cat(sprintf("*** S3 CLEANUP: First 5 protein groups: %s ***\n", paste(head(protein_groups, 5), collapse = ", ")))

            # Create lookup table for sequence information
            cat("*** S3 CLEANUP: Creating sequence lookup table ***\n")

            seq_lookup <- tryCatch(
                {
                    if ("length" %in% names(seqinr_obj)) {
                        cat("*** S3 CLEANUP: Using existing length column ***\n")
                        seqinr_obj |>
                            dplyr::select(dplyr::all_of(c(seqinr_accession_column, "length"))) |>
                            dplyr::distinct()
                    } else {
                        cat("*** S3 CLEANUP: Creating default length column ***\n")
                        # If no length column, create one or use a default
                        seqinr_obj |>
                            dplyr::select(dplyr::all_of(seqinr_accession_column)) |>
                            dplyr::distinct() |>
                            dplyr::mutate(length = 1000) # Default length
                    }
                },
                error = function(e) {
                    stop(sprintf("Error creating sequence lookup: %s", e$message))
                }
            )

            cat(sprintf("*** S3 CLEANUP: Created sequence lookup with %d entries ***\n", nrow(seq_lookup)))

            # Function to choose best protein from a group
            cat("*** S3 CLEANUP: Defining chooseBestFromGroup function ***\n")

            chooseBestFromGroup <- function(group_string) {
                tryCatch(
                    {
                        if (is.na(group_string) || group_string == "") {
                            return(NA_character_)
                        }

                        # Split the group
                        proteins <- stringr::str_split(group_string, delim)[[1]]
                        proteins <- stringr::str_trim(proteins) # Remove whitespace
                        proteins <- proteins[proteins != ""] # Remove empty strings

                        if (length(proteins) == 1) {
                            return(proteins[1])
                        }

                        # Find proteins that exist in sequence database
                        proteins_in_db <- proteins[proteins %in% seq_lookup[[seqinr_accession_column]]]

                        if (length(proteins_in_db) == 0) {
                            # No proteins found in database, return first one
                            return(proteins[1])
                        }

                        if (length(proteins_in_db) == 1) {
                            return(proteins_in_db[1])
                        }

                        # Multiple proteins in database - choose by length (longest first)
                        protein_info <- seq_lookup |>
                            dplyr::filter(!!sym(seqinr_accession_column) %in% proteins_in_db) |>
                            dplyr::arrange(desc(length), !!sym(seqinr_accession_column))

                        return(protein_info[[seqinr_accession_column]][1])
                    },
                    error = function(e) {
                        warning(sprintf("Error processing protein group '%s': %s", group_string, e$message))
                        return(group_string) # Return original if processing fails
                    }
                )
            }

            # Create mapping of original groups to best proteins
            cat("*** S3 CLEANUP: Creating protein group mapping ***\n")

            protein_mapping <- tryCatch(
                {
                    data.frame(
                        original_group = protein_groups,
                        best_protein = sapply(protein_groups, chooseBestFromGroup, USE.NAMES = FALSE),
                        stringsAsFactors = FALSE
                    )
                },
                error = function(e) {
                    stop(sprintf("Error creating protein mapping: %s", e$message))
                }
            )

            # Count reductions
            groups_with_multiple <- tryCatch(
                {
                    sum(stringr::str_detect(protein_mapping$original_group, delim), na.rm = TRUE)
                },
                error = function(e) {
                    cat(sprintf("*** S3 CLEANUP: Warning - could not count multi-protein groups: %s ***\n", e$message))
                    0
                }
            )

            cat(sprintf("*** S3 CLEANUP: Found %d protein groups with multiple proteins ***\n", groups_with_multiple))

            # Add mapping to data
            cat("*** S3 CLEANUP: Applying protein mapping to data ***\n")

            data_with_mapping <- tryCatch(
                {
                    data_tbl |>
                        dplyr::left_join(
                            protein_mapping,
                            by = setNames("original_group", protein_id_column)
                        )
                },
                error = function(e) {
                    stop(sprintf("Error applying protein mapping: %s", e$message))
                }
            )

            cat(sprintf("*** S3 CLEANUP: Data mapping completed. Mapped data has %d rows ***\n", nrow(data_with_mapping)))

            # Group by best protein and aggregate
            cat(sprintf("*** S3 CLEANUP: Aggregating data using method: %s ***\n", aggregation_method))

            # Define aggregation function
            agg_func <- switch(aggregation_method,
                "mean" = function(x) mean(x, na.rm = TRUE),
                "median" = function(x) median(x, na.rm = TRUE),
                "sum" = function(x) sum(x, na.rm = TRUE),
                function(x) mean(x, na.rm = TRUE) # Default to mean
            )

            # Separate quantitative and metadata columns for different aggregation
            cat("*** S3 CLEANUP: Aggregating quantitative data ***\n")

            quant_data <- tryCatch(
                {
                    if (length(quant_cols) > 0) {
                        data_with_mapping |>
                            dplyr::select(best_protein, dplyr::all_of(quant_cols)) |>
                            dplyr::group_by(best_protein) |>
                            dplyr::summarise(
                                dplyr::across(dplyr::all_of(quant_cols), agg_func),
                                .groups = "drop"
                            )
                    } else {
                        data.frame(best_protein = unique(data_with_mapping$best_protein))
                    }
                },
                error = function(e) {
                    stop(sprintf("Error aggregating quantitative data: %s", e$message))
                }
            )

            cat("*** S3 CLEANUP: Aggregating metadata ***\n")

            # For metadata columns, take the first non-NA value for each best protein
            meta_data <- tryCatch(
                {
                    if (length(meta_cols) > 0) {
                        data_with_mapping |>
                            dplyr::select(best_protein, dplyr::all_of(meta_cols)) |>
                            dplyr::group_by(best_protein) |>
                            dplyr::summarise(
                                dplyr::across(dplyr::all_of(meta_cols), ~ dplyr::first(na.omit(.))[1]),
                                .groups = "drop"
                            )
                    } else {
                        data.frame(best_protein = unique(data_with_mapping$best_protein))
                    }
                },
                error = function(e) {
                    stop(sprintf("Error aggregating metadata: %s", e$message))
                }
            )

            cat("*** S3 CLEANUP: Combining results ***\n")

            # Combine results
            result_data <- tryCatch(
                {
                    quant_data |>
                        dplyr::left_join(meta_data, by = "best_protein") |>
                        dplyr::rename(!!sym(protein_id_column) := best_protein)
                },
                error = function(e) {
                    stop(sprintf("Error combining results: %s", e$message))
                }
            )

            cat("*** S3 CLEANUP: Reordering columns ***\n")

            # Reorder columns to match original
            result_data <- tryCatch(
                {
                    result_data |>
                        dplyr::select(dplyr::all_of(names(data_tbl)))
                },
                error = function(e) {
                    stop(sprintf("Error reordering columns: %s", e$message))
                }
            )

            # Final statistics
            original_proteins <- length(unique(data_tbl[[protein_id_column]]))
            final_proteins <- length(unique(result_data[[protein_id_column]]))
            reduction_pct <- round(((original_proteins - final_proteins) / original_proteins) * 100, 1)

            cat(sprintf(
                "*** S3 CLEANUP: Protein cleanup completed: %d -> %d proteins (%.1f%% reduction) ***\n",
                original_proteins, final_proteins, reduction_pct
            ))

            return(result_data)
        },
        error = function(e) {
            cat(sprintf("*** S3 CLEANUP FATAL ERROR: %s ***\n", e$message))
            cat("*** S3 CLEANUP: Returning original data ***\n")
            return(data_tbl) # Return original data if everything fails
        }
    )
}

#' Create Study Parameters File
#'
#' @description
#' Creates a study parameters text file directly without using S4 objects.
#' This replaces the overly complex createWorkflowArgsFromConfig + WorkflowArgs show() approach.
#'
#' @param workflow_name Character string, name of the workflow
#' @param description Character string, description of the analysis
#' @param organism_name Character string, organism name (optional)
#' @param taxon_id Character string or numeric, taxon ID (optional)
#' @param source_dir_path Character string, path to save the study_parameters.txt file
#' @param contrasts_tbl Data frame, contrasts table (optional)
#' @param config_list_name Character string, name of the global config list (default: "config_list")
#' @param env Environment where the config list resides (default: .GlobalEnv)
#'
#' @return Character string path to the created study_parameters.txt file
#'
#' @examples
#' \dontrun{
#' file_path <- createStudyParametersFile(
#'     workflow_name = "proteomics_analysis",
#'     description = "DIA proteomics analysis",
#'     source_dir_path = "/path/to/scripts"
#' )
#' }
#'
#' @export
createStudyParametersFile <- function(workflow_name,
                                      description = "",
                                      organism_name = NULL,
                                      taxon_id = NULL,
                                      source_dir_path = NULL,
                                      contrasts_tbl = NULL,
                                      config_list_name = "config_list",
                                      env = .GlobalEnv) {
    # Validate required inputs
    if (missing(workflow_name) || !is.character(workflow_name) || length(workflow_name) != 1) {
        stop("workflow_name must be a single character string")
    }

    if (is.null(source_dir_path) || !is.character(source_dir_path) || length(source_dir_path) != 1) {
        stop("source_dir_path must be a single character string")
    }

    if (!dir.exists(source_dir_path)) {
        stop("source_dir_path directory does not exist: ", source_dir_path)
    }

    # Get git information
    git_info <- tryCatch(
        {
            if (requireNamespace("gh", quietly = TRUE)) {
                # Make the API call (works for public repos without token)
                branch_info <- gh::gh("/repos/APAF-BIOINFORMATICS/MultiScholaR/branches/main")
                list(
                    commit_sha = branch_info$commit$sha,
                    branch = "main",
                    repo = "MultiScholaR",
                    timestamp = branch_info$commit$commit$author$date
                )
            } else {
                list(commit_sha = NA_character_, branch = NA_character_, repo = NA_character_, timestamp = NA_character_)
            }
        },
        error = function(e) {
            message("Error fetching GitHub info: ", e$message)
            list(commit_sha = NA_character_, branch = NA_character_, repo = NA_character_, timestamp = NA_character_)
        }
    )

    # Get organism information from session if not provided
    if (is.null(organism_name) && exists("organism_name", envir = .GlobalEnv)) {
        organism_name <- get("organism_name", envir = .GlobalEnv)
    }

    if (is.null(taxon_id) && exists("taxon_id", envir = .GlobalEnv)) {
        taxon_id <- get("taxon_id", envir = .GlobalEnv)
    }

    # Get config list
    config_list <- if (exists(config_list_name, envir = env)) {
        get(config_list_name, envir = env)
    } else {
        list()
    }

    # Get contrasts_tbl if not provided
    if (is.null(contrasts_tbl)) {
        if (exists("contrasts_tbl", envir = parent.frame())) {
            contrasts_tbl <- get("contrasts_tbl", envir = parent.frame())
        } else if (exists("contrasts_tbl", envir = .GlobalEnv)) {
            contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
        }
    }

    # Build the output content
    output_lines <- c(
        "Study Parameters",
        "================",
        "",
        "Basic Information:",
        "-----------------",
        paste("Workflow Name:", workflow_name),
        paste("Description:", if (nzchar(description)) description else "N/A"),
        paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        ""
    )

    # Add git information
    if (!is.null(git_info) && is.list(git_info)) {
        git_lines <- c(
            "Git Information:",
            "---------------",
            paste("Repository:", ifelse(!is.null(git_info$repo) && !is.na(git_info$repo), git_info$repo, "N/A")),
            paste("Branch:", ifelse(!is.null(git_info$branch) && !is.na(git_info$branch), git_info$branch, "N/A")),
            paste("Commit:", ifelse(!is.null(git_info$commit_sha) && !is.na(git_info$commit_sha), substr(git_info$commit_sha, 1, 7), "N/A")),
            paste("Git Timestamp:", ifelse(!is.null(git_info$timestamp) && !is.na(git_info$timestamp), git_info$timestamp, "N/A")),
            ""
        )
        output_lines <- c(output_lines, git_lines)
    }

    # Add organism information
    if (!is.null(organism_name) && !is.na(organism_name) && nzchar(organism_name)) {
        organism_lines <- c(
            "Organism Information:",
            "---------------------",
            paste("Organism Name:", organism_name)
        )
        if (!is.null(taxon_id) && !is.na(taxon_id) && nzchar(taxon_id)) {
            organism_lines <- c(organism_lines, paste("Taxon ID:", taxon_id))
        }
        organism_lines <- c(organism_lines, "")
        output_lines <- c(output_lines, organism_lines)
    }

    # Add configuration parameters
    config_lines <- c(
        "Configuration Parameters:",
        "-------------------------"
    )

    # Clean the config list (remove problematic objects)
    clean_config <- config_list
    # Remove the internal source dir if it exists
    if (!is.null(clean_config$internal_workflow_source_dir)) {
        clean_config$internal_workflow_source_dir <- NULL
    }

    # Format the config list
    if (length(clean_config) > 0) {
        config_params <- formatConfigList(clean_config)
        config_lines <- c(config_lines, config_params)
    } else {
        config_lines <- c(config_lines, "No configuration parameters available")
    }

    output_lines <- c(output_lines, config_lines)

    # Add contrasts information
    if (!is.null(contrasts_tbl) && (is.data.frame(contrasts_tbl) || tibble::is_tibble(contrasts_tbl)) && nrow(contrasts_tbl) > 0) {
        contrasts_lines <- c(
            "",
            "Contrasts:",
            "----------"
        )

        if ("contrasts" %in% colnames(contrasts_tbl)) {
            contrasts_info <- tryCatch(
                {
                    contrasts_col <- contrasts_tbl[["contrasts"]]
                    paste("  ", as.character(contrasts_col))
                },
                error = function(e) {
                    paste("  [Error extracting contrasts:", e$message, "]")
                }
            )
            contrasts_lines <- c(contrasts_lines, contrasts_info)
        } else {
            contrasts_lines <- c(contrasts_lines, "  [Column 'contrasts' not found in contrasts_tbl]")
        }

        output_lines <- c(output_lines, contrasts_lines)
    }

    # Write to file
    output_file <- file.path(source_dir_path, "study_parameters.txt")

    tryCatch(
        {
            writeLines(output_lines, output_file)
            message("Study parameters saved to: ", output_file)
            return(output_file)
        },
        error = function(e) {
            stop("Failed to write study parameters file: ", e$message)
        }
    )
}

# S4 CLASS DEFINITION FOR RMD WORKFLOW
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

#' @title Format Configuration List
#' @param config_list List of configuration parameters
#' @param indent Number of spaces for indentation
#' @export
formatConfigList <- function(config_list, indent = 0) {
    output <- character()

    # Exclude internal_workflow_source_dir and optimization_results from printing
    names_to_process <- names(config_list)
    names_to_process <- names_to_process[!names_to_process %in% c("internal_workflow_source_dir", "optimization_results")]

    for (name in names_to_process) {
        value <- config_list[[name]]
        # Skip core_utilisation and complex objects from display
        if (name == "core_utilisation" ||
            any(class(value) %in% c("process", "R6", "multidplyr_cluster", "cluster", "SOCKcluster"))) {
            next
        }

        # Format the name
        name_formatted <- gsub("\\.", " ", name)
        name_formatted <- gsub("_", " ", name_formatted)
        name_formatted <- tools::toTitleCase(name_formatted)

        # Handle different value types
        if (is.list(value)) {
            if (length(value) > 0 && !is.null(names(value))) {
                output <- c(
                    output,
                    paste0(
                        paste(rep(" ", indent), collapse = ""),
                        name_formatted, ":"
                    )
                )
                output <- c(
                    output,
                    formatConfigList(value, indent + 2)
                )
            } else if (length(value) > 0) {
                output <- c(
                    output,
                    paste0(
                        paste(rep(" ", indent), collapse = ""),
                        name_formatted, ":"
                    )
                )
                for (item_idx in seq_along(value)) {
                    item_val <- value[[item_idx]]
                    if (is.atomic(item_val) && length(item_val) == 1) {
                        output <- c(output, paste0(paste(rep(" ", indent + 2), collapse = ""), "- ", as.character(item_val)))
                    } else {
                        output <- c(output, paste0(paste(rep(" ", indent + 2), collapse = ""), "- [Complex List Element]"))
                    }
                }
            } else {
                output <- c(
                    output,
                    paste0(
                        paste(rep(" ", indent), collapse = ""),
                        name_formatted, ": [Empty List]"
                    )
                )
            }
        } else {
            # Truncate long character vectors for display
            if (is.character(value) && length(value) > 5) {
                value_display <- paste(c(utils::head(value, 5), "..."), collapse = ", ")
            } else if (is.character(value) && length(value) > 1) {
                value_display <- paste(value, collapse = ", ")
            } else {
                value_display <- as.character(value)
                if (length(value_display) == 0) value_display <- "[Empty/NULL]"
            }
            output <- c(
                output,
                paste0(
                    paste(rep(" ", indent), collapse = ""),
                    name_formatted, ": ", value_display
                )
            )
        }
    }
    return(output)
}

setMethod(
    "show",
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
            paste("Description:", if (nzchar(object@description)) object@description else "N/A"),
            paste("Timestamp:", object@timestamp),
            ""
        )

        git_info_vec <- character()
        if (!is.null(object@git_info) && is.list(object@git_info)) {
            gi <- object@git_info
            git_info_vec <- c(
                paste("Repository:", ifelse(!is.null(gi$repo) && !is.na(gi$repo), gi$repo, "N/A")),
                paste("Branch:", ifelse(!is.null(gi$branch) && !is.na(gi$branch), gi$branch, "N/A")),
                paste("Commit:", ifelse(!is.null(gi$commit_sha) && !is.na(gi$commit_sha), substr(gi$commit_sha, 1, 7), "N/A")),
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

        # Check for optimization results
        optimization_section <- character()
        if (!is.null(object@args$optimization_results) && is.list(object@args$optimization_results)) {
            opt_res <- object@args$optimization_results

            opt_lines <- c(
                "Automatic RUV Optimization Results:",
                "-----------------------------------"
            )

            # Format optimization results
            tryCatch(
                {
                    if (!is.null(opt_res$best_percentage)) {
                        opt_lines <- c(opt_lines, paste("• Best percentage:", sprintf("%.1f%%", opt_res$best_percentage)))
                    }
                    if (!is.null(opt_res$best_k)) {
                        opt_lines <- c(opt_lines, paste("• Best k value:", opt_res$best_k))
                    }
                    if (!is.null(opt_res$best_separation_score)) {
                        opt_lines <- c(opt_lines, paste("• Separation score:", sprintf("%.4f", opt_res$best_separation_score)))
                    }
                    if (!is.null(opt_res$best_composite_score)) {
                        opt_lines <- c(opt_lines, paste("• Composite score:", sprintf("%.4f", opt_res$best_composite_score)))
                    }
                    if (!is.null(opt_res$best_control_genes_index)) {
                        control_count <- sum(opt_res$best_control_genes_index, na.rm = TRUE)
                        opt_lines <- c(opt_lines, paste("• Control genes:", control_count))
                    }
                    if (!is.null(opt_res$separation_metric_used)) {
                        opt_lines <- c(opt_lines, paste("• Separation metric:", opt_res$separation_metric_used))
                    }
                    if (!is.null(opt_res$k_penalty_weight)) {
                        opt_lines <- c(opt_lines, paste("• K penalty weight:", sprintf("%.1f", opt_res$k_penalty_weight)))
                    }
                    if (!is.null(opt_res$adaptive_k_penalty_used)) {
                        opt_lines <- c(opt_lines, paste("• Adaptive penalty:", ifelse(opt_res$adaptive_k_penalty_used, "TRUE", "FALSE")))
                    }
                    if (!is.null(opt_res$sample_size) && !is.na(opt_res$sample_size)) {
                        opt_lines <- c(opt_lines, paste("• Sample size:", opt_res$sample_size))
                    }

                    optimization_section <- c(opt_lines, "")
                },
                error = function(e) {
                    optimization_section <- c(
                        "Automatic RUV Optimization Results:",
                        "-----------------------------------",
                        paste("• [Error formatting optimization results:", e$message, "]"),
                        ""
                    )
                }
            )
        }

        config_header <- c(
            "Configuration Parameters (from @args slot):",
            "-----------------------------------------"
        )

        args_to_print <- object@args
        # internal_workflow_source_dir and optimization_results are excluded
        # (optimization_results is shown in its own section above)
        if (!is.null(args_to_print$optimization_results)) {
            args_to_print$optimization_results <- NULL
        }

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

        output <- c(header, git_section, organism_section, optimization_section, config_header, params, contrasts_section)
        cat(paste(output, collapse = "\n"), "\n")

        # Save to file if internal_workflow_source_dir is defined in args and is valid
        workflow_source_dir_to_use <- object@args$internal_workflow_source_dir

        if (!is.null(workflow_source_dir_to_use) &&
            is.character(workflow_source_dir_to_use) &&
            length(workflow_source_dir_to_use) == 1 &&
            nzchar(workflow_source_dir_to_use) &&
            dir.exists(workflow_source_dir_to_use)) {
            output_file <- file.path(workflow_source_dir_to_use, "study_parameters.txt")
            tryCatch(
                {
                    writeLines(output, output_file)
                    cat("\nParameters saved to:", output_file, "\n")
                },
                error = function(e) {
                    warning(paste("Failed to save study parameters to", output_file, ":", e$message), immediate. = TRUE)
                }
            )
        } else {
            cat("\nWarning: `internal_workflow_source_dir` not found, invalid, or directory does not exist in WorkflowArgs @args. Parameters not saved to file.\n")
            if (!is.null(workflow_source_dir_to_use)) cat("Attempted path:", workflow_source_dir_to_use, "\n")
        }
    }
)

# Updated function to extract parameters from S4 @args slot
#' @export
#' @param workflow_name Character string, name of the workflow
#' @param description Character string, description of the analysis
#' @param organism_name Character string, organism name (optional)
#' @param taxon_id Character string or numeric, taxon ID (optional)
#' @param source_dir_path Character string, path to save the study_parameters.txt file
#' @param final_s4_object S4 object containing workflow parameters in @args slot (optional)
#' @param contrasts_tbl Data frame, contrasts table (optional)
#' @param workflow_data List containing workflow state and optimization results (optional)
createWorkflowArgsFromConfig <- function(workflow_name, description = "",
                                         organism_name = NULL, taxon_id = NULL,
                                         source_dir_path = NULL,
                                         final_s4_object = NULL,
                                         contrasts_tbl = NULL,
                                         workflow_data = NULL) {
    # WORKFLOW DETECTION: If workflow_data is NULL, use original RMD logic
    if (is.null(workflow_data)) {
        cat("WORKFLOW ARGS: Detected RMD workflow (workflow_data is NULL), using original logic\n")
        return(.createWorkflowArgsFromConfig_RMD(
            workflow_name = workflow_name,
            description = description,
            organism_name = organism_name,
            taxon_id = taxon_id,
            source_dir_path = source_dir_path
        ))
    }

    # SHINY WORKFLOW: Continue with existing Shiny logic (unchanged)
    cat("WORKFLOW ARGS: Detected Shiny workflow (workflow_data provided), using Shiny logic\n")

    # Validate required inputs
    if (missing(workflow_name) || !is.character(workflow_name) || length(workflow_name) != 1) {
        stop("workflow_name must be a single character string")
    }

    if (is.null(source_dir_path) || !is.character(source_dir_path) || length(source_dir_path) != 1) {
        stop("source_dir_path must be a single character string")
    }

    if (!dir.exists(source_dir_path)) {
        stop("source_dir_path directory does not exist: ", source_dir_path)
    }

    cat("WORKFLOW ARGS: Starting parameter extraction\n")

    # Extract parameters from S4 object @args if provided
    s4_params <- list()
    # Check if final_s4_object has @args slot
    s4_has_args <- tryCatch(
        {
            !is.null(final_s4_object) && isS4(final_s4_object) && !is.null(final_s4_object@args)
        },
        error = function(e) {
            FALSE
        }
    )

    if (s4_has_args) {
        cat("WORKFLOW ARGS: Extracting parameters from S4 @args slot\n")

        tryCatch(
            {
                s4_params <- final_s4_object@args
                cat(sprintf("WORKFLOW ARGS: Found %d function groups in S4 @args\n", length(s4_params)))

                # Log the function groups for debugging
                for (func_name in names(s4_params)) {
                    param_count <- length(s4_params[[func_name]])
                    cat(sprintf("WORKFLOW ARGS: Function '%s' has %d parameters\n", func_name, param_count))
                }
            },
            error = function(e) {
                cat(sprintf("WORKFLOW ARGS: Error extracting S4 params: %s\n", e$message))
                s4_params <- list()
            }
        )
    } else {
        cat("WORKFLOW ARGS: No S4 object provided or S4 @args is NULL\n")
    }

    # Get fallback config_list from global environment
    config_list <- if (exists("config_list", envir = .GlobalEnv)) {
        get("config_list", envir = .GlobalEnv)
    } else {
        list()
    }

    cat(sprintf("WORKFLOW ARGS: Global config_list has %d sections\n", length(config_list)))

    # Extract RUV optimization results from workflow_data if available
    ruv_optimization_result <- NULL
    if (!is.null(workflow_data)) {
        cat("WORKFLOW ARGS: Checking for RUV optimization results in workflow_data\n")

        # Check multiple possible locations for RUV optimization results
        if (!is.null(workflow_data$ruv_optimization_result)) {
            ruv_optimization_result <- workflow_data$ruv_optimization_result
            cat("WORKFLOW ARGS: Found RUV optimization results in workflow_data$ruv_optimization_result\n")
        } else if (!is.null(workflow_data$state_manager)) {
            # Try to get RUV results from state manager config
            tryCatch(
                {
                    current_state_config <- workflow_data$state_manager$getStateConfig(workflow_data$state_manager$current_state)
                    if (!is.null(current_state_config$ruv_optimization_result)) {
                        ruv_optimization_result <- current_state_config$ruv_optimization_result
                        cat("WORKFLOW ARGS: Found RUV optimization results in state manager config\n")
                    }
                },
                error = function(e) {
                    cat(sprintf("WORKFLOW ARGS: Could not extract RUV results from state manager: %s\n", e$message))
                }
            )
        }

        # Also check if there's a saved RUV file
        if (is.null(ruv_optimization_result) && !is.null(source_dir_path)) {
            ruv_file <- file.path(source_dir_path, "ruv_optimization_results.RDS")
            if (file.exists(ruv_file)) {
                tryCatch(
                    {
                        ruv_optimization_result <- readRDS(ruv_file)
                        cat(sprintf("WORKFLOW ARGS: Loaded RUV optimization results from file: %s\n", ruv_file))
                    },
                    error = function(e) {
                        cat(sprintf("WORKFLOW ARGS: Could not load RUV file: %s\n", e$message))
                    }
                )
            }
        }
    }

    # Check for workflow-specific optimization results in R environment
    # These take precedence over workflow_data results
    cat("WORKFLOW ARGS: Checking for workflow-specific optimization results in R environment\n")

    if (workflow_name == "DIA_limpa") {
        cat("WORKFLOW ARGS: Workflow is DIA_limpa, checking for optimal_results_peptides\n")

        # Try to find optimal_results_peptides in R environment
        optimal_results_found <- FALSE

        # Check parent frame first, then global environment
        if (exists("optimal_results_peptides", envir = parent.frame(), inherits = FALSE)) {
            tryCatch(
                {
                    ruv_optimization_result <- get("optimal_results_peptides", envir = parent.frame())
                    optimal_results_found <- TRUE
                    cat("WORKFLOW ARGS: Found optimal_results_peptides in parent frame\n")
                },
                error = function(e) {
                    cat(sprintf("WORKFLOW ARGS: Error getting optimal_results_peptides from parent frame: %s\n", e$message))
                }
            )
        }

        if (!optimal_results_found && exists("optimal_results_peptides", envir = .GlobalEnv)) {
            tryCatch(
                {
                    ruv_optimization_result <- get("optimal_results_peptides", envir = .GlobalEnv)
                    optimal_results_found <- TRUE
                    cat("WORKFLOW ARGS: Found optimal_results_peptides in global environment\n")
                },
                error = function(e) {
                    cat(sprintf("WORKFLOW ARGS: Error getting optimal_results_peptides from global environment: %s\n", e$message))
                }
            )
        }

        if (!optimal_results_found) {
            cat("WORKFLOW ARGS: optimal_results_peptides not found in R environment for DIA_limpa workflow\n")
        }
    } else if (workflow_name == "DIA") {
        cat("WORKFLOW ARGS: Workflow is DIA\n")
        # TODO: Determine object name for DIA workflow optimization results
        # TODO: Implement logic similar to DIA_limpa once object name is determined
        # Example structure:
        # if (exists("optimal_results_OBJECTNAME", envir = .GlobalEnv)) {
        #     ruv_optimization_result <- get("optimal_results_OBJECTNAME", envir = .GlobalEnv)
        # }
    } else if (workflow_name == "TMT") {
        cat("WORKFLOW ARGS: Workflow is TMT\n")
        # TODO: Determine object name for TMT workflow optimization results
        # TODO: Implement logic similar to DIA_limpa once object name is determined
        # Example structure:
        # if (exists("optimal_results_OBJECTNAME", envir = .GlobalEnv)) {
        #     ruv_optimization_result <- get("optimal_results_OBJECTNAME", envir = .GlobalEnv)
        # }
    }

    # Merge S4 parameters with config_list (S4 takes precedence)
    merged_config <- config_list

    if (length(s4_params) > 0) {
        cat("WORKFLOW ARGS: Merging S4 parameters with config_list\n")

        # S4 parameters override config_list parameters
        for (func_name in names(s4_params)) {
            if (is.list(s4_params[[func_name]]) && length(s4_params[[func_name]]) > 0) {
                merged_config[[func_name]] <- s4_params[[func_name]]
                cat(sprintf(
                    "WORKFLOW ARGS: Updated '%s' section with %d S4 parameters\n",
                    func_name, length(s4_params[[func_name]])
                ))
            }
        }
    }

    # Get git information
    git_info <- tryCatch(
        {
            if (requireNamespace("gh", quietly = TRUE)) {
                # Make the API call (works for public repos without token)
                branch_info <- gh::gh("/repos/APAF-BIOINFORMATICS/MultiScholaR/branches/main")
                list(
                    commit_sha = branch_info$commit$sha,
                    branch = "main",
                    repo = "MultiScholaR",
                    timestamp = branch_info$commit$commit$author$date
                )
            } else {
                list(commit_sha = NA_character_, branch = NA_character_, repo = NA_character_, timestamp = NA_character_)
            }
        },
        error = function(e) {
            message("Error fetching GitHub info: ", e$message)
            list(commit_sha = NA_character_, branch = NA_character_, repo = NA_character_, timestamp = NA_character_)
        }
    )

    # Get organism information from session if not provided
    if (is.null(organism_name) && exists("organism_name", envir = .GlobalEnv)) {
        organism_name <- get("organism_name", envir = .GlobalEnv)
    }

    if (is.null(taxon_id) && exists("taxon_id", envir = .GlobalEnv)) {
        taxon_id <- get("taxon_id", envir = .GlobalEnv)
    }

    # Get contrasts_tbl if not provided
    if (is.null(contrasts_tbl)) {
        if (exists("contrasts_tbl", envir = parent.frame())) {
            contrasts_tbl <- get("contrasts_tbl", envir = parent.frame())
        } else if (exists("contrasts_tbl", envir = .GlobalEnv)) {
            contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
        }
    }

    # Build the output content
    output_lines <- c(
        "Study Parameters",
        "================",
        "",
        "Basic Information:",
        "-----------------",
        paste("Workflow Name:", workflow_name),
        paste("Description:", if (nzchar(description)) description else "N/A"),
        paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        ""
    )

    # Add git information
    if (!is.null(git_info) && is.list(git_info)) {
        git_lines <- c(
            "Git Information:",
            "---------------",
            paste("Repository:", ifelse(!is.null(git_info$repo) && !is.na(git_info$repo), git_info$repo, "N/A")),
            paste("Branch:", ifelse(!is.null(git_info$branch) && !is.na(git_info$branch), git_info$branch, "N/A")),
            paste("Commit:", ifelse(!is.null(git_info$commit_sha) && !is.na(git_info$commit_sha), substr(git_info$commit_sha, 1, 7), "N/A")),
            paste("Git Timestamp:", ifelse(!is.null(git_info$timestamp) && !is.na(git_info$timestamp), git_info$timestamp, "N/A")),
            ""
        )
        output_lines <- c(output_lines, git_lines)
    }

    # Add organism information
    if (!is.null(organism_name) && !is.na(organism_name) && nzchar(organism_name)) {
        organism_lines <- c(
            "Organism Information:",
            "---------------------",
            paste("Organism Name:", organism_name)
        )
        if (!is.null(taxon_id) && !is.na(taxon_id) && nzchar(taxon_id)) {
            organism_lines <- c(organism_lines, paste("Taxon ID:", taxon_id))
        }
        organism_lines <- c(organism_lines, "")
        output_lines <- c(output_lines, organism_lines)
    }

    # Add RUV optimization results if available
    if (!is.null(ruv_optimization_result) && is.list(ruv_optimization_result)) {
        cat("WORKFLOW ARGS: Formatting RUV optimization results\n")

        ruv_lines <- c(
            "Automatic RUV Optimization Results:",
            "-----------------------------------"
        )

        # Extract and format RUV optimization values
        tryCatch(
            {
                best_percentage <- if (!is.null(ruv_optimization_result$best_percentage)) {
                    sprintf("%.1f%%", ruv_optimization_result$best_percentage)
                } else {
                    "N/A"
                }

                best_k <- if (!is.null(ruv_optimization_result$best_k)) {
                    as.character(ruv_optimization_result$best_k)
                } else {
                    "N/A"
                }

                separation_score <- if (!is.null(ruv_optimization_result$best_separation_score)) {
                    sprintf("%.4f", ruv_optimization_result$best_separation_score)
                } else {
                    "N/A"
                }

                composite_score <- if (!is.null(ruv_optimization_result$best_composite_score)) {
                    sprintf("%.4f", ruv_optimization_result$best_composite_score)
                } else {
                    "N/A"
                }

                control_genes_count <- if (!is.null(ruv_optimization_result$best_control_genes_index)) {
                    as.character(sum(ruv_optimization_result$best_control_genes_index, na.rm = TRUE))
                } else {
                    "N/A"
                }

                separation_metric <- if (!is.null(ruv_optimization_result$separation_metric_used)) {
                    as.character(ruv_optimization_result$separation_metric_used)
                } else {
                    "N/A"
                }

                k_penalty_weight <- if (!is.null(ruv_optimization_result$k_penalty_weight)) {
                    sprintf("%.1f", ruv_optimization_result$k_penalty_weight)
                } else {
                    "N/A"
                }

                adaptive_penalty <- if (!is.null(ruv_optimization_result$adaptive_k_penalty_used)) {
                    ifelse(ruv_optimization_result$adaptive_k_penalty_used, "TRUE", "FALSE")
                } else {
                    "N/A"
                }

                sample_size <- if (!is.null(ruv_optimization_result$sample_size)) {
                    as.character(ruv_optimization_result$sample_size)
                } else {
                    "N/A"
                }

                # Also try to get RUV grouping variable from S4 parameters
                ruv_grouping_variable <- "N/A"
                if (!is.null(s4_params$ruvIII_C_Varying$ruv_grouping_variable)) {
                    ruv_grouping_variable <- s4_params$ruvIII_C_Varying$ruv_grouping_variable
                } else if (!is.null(s4_params$getNegCtrlProtAnova$ruv_grouping_variable)) {
                    ruv_grouping_variable <- s4_params$getNegCtrlProtAnova$ruv_grouping_variable
                }

                ruv_lines <- c(
                    ruv_lines,
                    paste("• Best percentage:", best_percentage),
                    paste("• Best k value:", best_k),
                    paste("• Separation score:", separation_score),
                    paste("• Composite score:", composite_score),
                    paste("• Control genes:", control_genes_count),
                    paste("• RUV grouping variable:", ruv_grouping_variable),
                    paste("• Separation metric:", separation_metric),
                    paste("• K penalty weight:", k_penalty_weight),
                    paste("• Adaptive penalty:", adaptive_penalty),
                    paste("• Sample size:", sample_size),
                    ""
                )

                cat("WORKFLOW ARGS: Successfully formatted RUV optimization results\n")
            },
            error = function(e) {
                cat(sprintf("WORKFLOW ARGS: Error formatting RUV results: %s\n", e$message))
                ruv_lines <- c(
                    ruv_lines,
                    paste("• [Error formatting RUV optimization results:", e$message, "]"),
                    ""
                )
            }
        )

        output_lines <- c(output_lines, ruv_lines)
    } else {
        cat("WORKFLOW ARGS: No RUV optimization results available\n")
    }

    # Add configuration parameters (now includes S4 parameters)
    config_lines <- c(
        "Workflow Parameters:",
        "-------------------",
        ""
    )

    # Clean the merged config (remove problematic objects)
    clean_config <- merged_config
    # Remove the internal source dir if it exists
    if (!is.null(clean_config$internal_workflow_source_dir)) {
        clean_config$internal_workflow_source_dir <- NULL
    }

    # Add special section for S4-derived parameters
    if (length(s4_params) > 0) {
        cat("WORKFLOW ARGS: About to format S4 parameters using functional approach\n")

        config_lines <- c(
            config_lines,
            "Parameters from Final S4 Object:",
            "--------------------------------"
        )

        # SAFE parameter formatting function
        formatParameterValue <- function(param_value) {
            tryCatch(
                {
                    if (is.null(param_value)) {
                        "NULL"
                    } else if (is.logical(param_value)) {
                        if (length(param_value) == 1) {
                            ifelse(param_value, "TRUE", "FALSE")
                        } else if (length(param_value) > 50) {
                            # Handle large logical vectors (like control genes index)
                            true_count <- sum(param_value, na.rm = TRUE)
                            total_count <- length(param_value)
                            sprintf(
                                "logical vector [%d TRUE, %d FALSE out of %d total]",
                                true_count, total_count - true_count, total_count
                            )
                        } else {
                            # Show first few values for smaller vectors
                            preview <- ifelse(utils::head(param_value, 5), "TRUE", "FALSE")
                            if (length(param_value) > 5) {
                                paste0("c(", paste(preview, collapse = ", "), ", ...)")
                            } else {
                                paste0("c(", paste(preview, collapse = ", "), ")")
                            }
                        }
                    } else if (is.numeric(param_value)) {
                        if (length(param_value) == 1) {
                            as.character(param_value)
                        } else if (length(param_value) > 5) {
                            sprintf(
                                "numeric vector [%d values: %s, ...]",
                                length(param_value),
                                paste(as.character(utils::head(param_value, 3)), collapse = ", ")
                            )
                        } else {
                            paste0("c(", paste(as.character(param_value), collapse = ", "), ")")
                        }
                    } else if (is.character(param_value)) {
                        if (length(param_value) == 1) {
                            param_value
                        } else if (length(param_value) > 5) {
                            sprintf(
                                "character vector [%d values: %s, ...]",
                                length(param_value),
                                paste(shQuote(utils::head(param_value, 3)), collapse = ", ")
                            )
                        } else {
                            paste0("c(", paste(shQuote(param_value), collapse = ", "), ")")
                        }
                    } else {
                        # SAFE fallback - no dput() which was causing the hang
                        paste0("[", class(param_value)[1], " object]")
                    }
                },
                error = function(e) {
                    "[SERIALIZATION ERROR]"
                }
            )
        }

        s4_sections <- if (requireNamespace("purrr", quietly = TRUE)) {
            purrr::imap(s4_params, function(func_params, func_name) {
                tryCatch(
                    {
                        if (!is.list(func_params) || length(func_params) == 0) {
                            return("")
                        }

                        cat(sprintf("WORKFLOW ARGS: Processing S4 function group '%s'\n", func_name))
                        header <- sprintf("[%s]", func_name)

                        # Use map instead of imap_chr for safer handling
                        param_lines <- purrr::imap(func_params, function(param_value, param_name) {
                            tryCatch(
                                {
                                    param_str <- formatParameterValue(param_value)
                                    sprintf("  %s = %s", param_name, param_str)
                                },
                                error = function(e) {
                                    cat(sprintf("WORKFLOW ARGS: Error formatting parameter '%s' in '%s': %s\n", param_name, func_name, e$message))
                                    sprintf("  %s = [ERROR: %s]", param_name, e$message)
                                }
                            )
                        })

                        # Convert list to character vector safely
                        param_lines_char <- unlist(param_lines)
                        paste(c(header, param_lines_char, ""), collapse = "\n")
                    },
                    error = function(e) {
                        cat(sprintf("WORKFLOW ARGS: Error processing S4 function group '%s': %s\n", func_name, e$message))
                        sprintf("[%s]\n  [ERROR: %s]\n", func_name, e$message)
                    }
                )
            })
        } else {
            # Fallback using base R if purrr not available
            lapply(names(s4_params), function(func_name) {
                func_params <- s4_params[[func_name]]
                if (!is.list(func_params) || length(func_params) == 0) {
                    return("")
                }

                header <- sprintf("[%s]", func_name)

                param_lines <- lapply(names(func_params), function(param_name) {
                    param_value <- func_params[[param_name]]
                    param_str <- formatParameterValue(param_value)
                    sprintf("  %s = %s", param_name, param_str)
                })

                paste(c(header, unlist(param_lines), ""), collapse = "\n")
            })
        }

        # Add all S4 sections at once - safely handle list from purrr::imap()
        s4_sections_char <- tryCatch(
            {
                if (is.list(s4_sections)) {
                    # Convert list to character vector
                    unlist(s4_sections)
                } else {
                    s4_sections
                }
            },
            error = function(e) {
                cat(sprintf("WORKFLOW ARGS: Error converting S4 sections: %s\n", e$message))
                "[ERROR: Could not process S4 sections]"
            }
        )

        config_lines <- c(config_lines, unlist(strsplit(paste(s4_sections_char, collapse = ""), "\n")))

        cat("WORKFLOW ARGS: S4 parameters formatted successfully\n")
    }

    # Add UI parameters from workflow_data if available (DE and Enrichment UI inputs)
    if (!is.null(workflow_data)) {
        cat("WORKFLOW ARGS: Checking for UI parameters in workflow_data\n")
        ui_sections <- c()

        # Check for DE UI parameters
        if (!is.null(workflow_data$de_ui_params)) {
            cat("WORKFLOW ARGS: Found DE UI parameters in workflow_data\n")
            de_ui_lines <- c(
                "[Differential Expression UI Parameters]",
                sprintf("  q_value_threshold = %s", ifelse(!is.null(workflow_data$de_ui_params$q_value_threshold), workflow_data$de_ui_params$q_value_threshold, "N/A")),
                sprintf("  log_fold_change_cutoff = %s", ifelse(!is.null(workflow_data$de_ui_params$log_fold_change_cutoff), workflow_data$de_ui_params$log_fold_change_cutoff, "N/A")),
                sprintf("  treat_enabled = %s", ifelse(!is.null(workflow_data$de_ui_params$treat_enabled), workflow_data$de_ui_params$treat_enabled, "N/A")),
                ""
            )
            ui_sections <- c(ui_sections, de_ui_lines)
        }

        # Check for Enrichment UI parameters
        if (!is.null(workflow_data$enrichment_ui_params)) {
            cat("WORKFLOW ARGS: Found Enrichment UI parameters in workflow_data\n")
            enrichment_ui_lines <- c(
                "[Enrichment Analysis UI Parameters]",
                sprintf("  up_log2fc_cutoff = %s", ifelse(!is.null(workflow_data$enrichment_ui_params$up_log2fc_cutoff), workflow_data$enrichment_ui_params$up_log2fc_cutoff, "N/A")),
                sprintf("  down_log2fc_cutoff = %s", ifelse(!is.null(workflow_data$enrichment_ui_params$down_log2fc_cutoff), workflow_data$enrichment_ui_params$down_log2fc_cutoff, "N/A")),
                sprintf("  q_value_cutoff = %s", ifelse(!is.null(workflow_data$enrichment_ui_params$q_value_cutoff), workflow_data$enrichment_ui_params$q_value_cutoff, "N/A")),
                sprintf("  organism_selected = %s", ifelse(!is.null(workflow_data$enrichment_ui_params$organism_selected), workflow_data$enrichment_ui_params$organism_selected, "N/A")),
                sprintf("  database_source = %s", ifelse(!is.null(workflow_data$enrichment_ui_params$database_source), workflow_data$enrichment_ui_params$database_source, "N/A")),
                ""
            )
            ui_sections <- c(ui_sections, enrichment_ui_lines)
        }

        # Add UI sections to config_lines if any were found
        if (length(ui_sections) > 0) {
            config_lines <- c(
                config_lines,
                "User Interface Parameters:",
                "-------------------------",
                ui_sections
            )
            cat("WORKFLOW ARGS: Added UI parameters to output\n")
        } else {
            cat("WORKFLOW ARGS: No UI parameters found in workflow_data\n")
        }
    }

    # Format the remaining config list
    cat("WORKFLOW ARGS: About to format remaining config list\n")
    if (length(clean_config) > 0) {
        config_lines <- c(
            config_lines,
            "Additional Configuration Parameters:",
            "-----------------------------------"
        )

        cat("WORKFLOW ARGS: Calling formatConfigList...\n")
        tryCatch(
            {
                # Check if formatConfigList exists before calling
                if (exists("formatConfigList", mode = "function")) {
                    config_params <- formatConfigList(clean_config)
                    config_lines <- c(config_lines, config_params)
                    cat("WORKFLOW ARGS: formatConfigList completed successfully\n")
                } else {
                    cat("WORKFLOW ARGS: formatConfigList function not found, using basic formatting\n")
                    # Basic fallback formatting without for loops
                    basic_config <- unlist(clean_config, recursive = TRUE)
                    config_lines <- c(config_lines, names(basic_config), " = ", as.character(basic_config))
                }
            },
            error = function(e) {
                cat(sprintf("WORKFLOW ARGS: formatConfigList failed: %s\n", e$message))
                config_lines <- c(config_lines, paste("Error formatting config:", e$message))
            }
        )
    } else {
        config_lines <- c(config_lines, "No additional configuration parameters available")
    }

    cat("WORKFLOW ARGS: Adding config lines to output\n")
    output_lines <- c(output_lines, config_lines)

    # Add contrasts information
    cat("WORKFLOW ARGS: Processing contrasts information\n")
    if (!is.null(contrasts_tbl) && (is.data.frame(contrasts_tbl) || tibble::is_tibble(contrasts_tbl)) && nrow(contrasts_tbl) > 0) {
        cat("WORKFLOW ARGS: Adding contrasts to output\n")
        contrasts_lines <- c(
            "",
            "Contrasts:",
            "----------"
        )

        if ("contrasts" %in% colnames(contrasts_tbl)) {
            contrasts_info <- tryCatch(
                {
                    contrasts_col <- contrasts_tbl[["contrasts"]]
                    paste("  ", as.character(contrasts_col))
                },
                error = function(e) {
                    paste("  [Error extracting contrasts:", e$message, "]")
                }
            )
            contrasts_lines <- c(contrasts_lines, contrasts_info)
        } else {
            contrasts_lines <- c(contrasts_lines, "  [Column 'contrasts' not found in contrasts_tbl]")
        }

        output_lines <- c(output_lines, contrasts_lines)
        cat("WORKFLOW ARGS: Contrasts added successfully\n")
    } else {
        cat("WORKFLOW ARGS: No contrasts to add\n")
    }

    # Write to file
    cat("WORKFLOW ARGS: About to write file\n")
    output_file <- file.path(source_dir_path, "study_parameters.txt")
    cat(sprintf("WORKFLOW ARGS: Target file path: %s\n", output_file))

    tryCatch(
        {
            cat("WORKFLOW ARGS: Calling writeLines...\n")
            writeLines(output_lines, output_file)
            cat(sprintf("WORKFLOW ARGS: Study parameters saved to: %s\n", output_file))
            return(output_file)
        },
        error = function(e) {
            cat(sprintf("WORKFLOW ARGS: Error writing file: %s\n", e$message))
            stop("Failed to write study parameters file: ", e$message)
        }
    )
}

# ORIGINAL RMD WORKFLOW FUNCTION (SEPARATE TO PRESERVE SHINY LOGIC)
.createWorkflowArgsFromConfig_RMD <- function(workflow_name, description = "",
                                              organism_name = NULL, taxon_id = NULL,
                                              source_dir_path = NULL) {
    cat("WORKFLOW ARGS: Using original RMD logic to create WorkflowArgs S4 object\n")

    # Get git information (original logic with token checking)
    git_info <- tryCatch(
        {
            if (requireNamespace("gh", quietly = TRUE)) {
                # Check if a token is available to avoid interactive prompts if not in an interactive session
                if (nzchar(gh::gh_token())) {
                    branch_info <- gh::gh("/repos/APAF-BIOINFORMATICS/MultiScholaR/branches/main")
                    list(
                        commit_sha = branch_info$commit$sha,
                        branch = "main",
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
        },
        error = function(e) {
            warning(paste("Could not retrieve Git information:", e$message), immediate. = TRUE)
            list(commit_sha = NA_character_, branch = NA_character_, repo = NA_character_, timestamp = NA_character_)
        }
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
        organism_name = if (is.null(organism_name) || !is.character(organism_name)) NA_character_ else organism_name,
        taxon_id = if (is.null(taxon_id) || !(is.character(taxon_id) || is.numeric(taxon_id))) NA_character_ else as.character(taxon_id)
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

    # Check for workflow-specific optimization results in R environment
    cat("WORKFLOW ARGS: Checking for workflow-specific optimization results in R environment\n")

    optimization_result <- NULL
    if (workflow_name == "DIA_limpa") {
        cat("WORKFLOW ARGS: Workflow is DIA_limpa, checking for optimal_results_peptides\n")

        # Try to find optimal_results_peptides in R environment
        if (exists("optimal_results_peptides", envir = parent.frame(), inherits = FALSE)) {
            tryCatch(
                {
                    optimization_result <- get("optimal_results_peptides", envir = parent.frame())
                    cat("WORKFLOW ARGS: Found optimal_results_peptides in parent frame\n")
                },
                error = function(e) {
                    cat(sprintf("WORKFLOW ARGS: Error getting optimal_results_peptides from parent frame: %s\n", e$message))
                }
            )
        }

        if (is.null(optimization_result) && exists("optimal_results_peptides", envir = .GlobalEnv)) {
            tryCatch(
                {
                    optimization_result <- get("optimal_results_peptides", envir = .GlobalEnv)
                    cat("WORKFLOW ARGS: Found optimal_results_peptides in global environment\n")
                },
                error = function(e) {
                    cat(sprintf("WORKFLOW ARGS: Error getting optimal_results_peptides from global environment: %s\n", e$message))
                }
            )
        }

        if (is.null(optimization_result)) {
            cat("WORKFLOW ARGS: optimal_results_peptides not found in R environment for DIA_limpa workflow\n")
        }
    } else if (workflow_name == "DIA") {
        cat("WORKFLOW ARGS: Workflow is DIA\n")
        # TODO: Determine object name for DIA workflow optimization results
        # TODO: Implement logic similar to DIA_limpa once object name is determined
    } else if (workflow_name == "TMT") {
        cat("WORKFLOW ARGS: Workflow is TMT\n")
        # TODO: Determine object name for TMT workflow optimization results
        # TODO: Implement logic similar to DIA_limpa once object name is determined
    }

    # Store optimization results in config_list if found
    if (!is.null(optimization_result)) {
        config_list_to_use$optimization_results <- optimization_result
        cat("WORKFLOW ARGS: Stored optimization results in config_list\n")
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

    cat("WORKFLOW ARGS: Created WorkflowArgs S4 object for RMD workflow\n")
    return(new_workflow_args)
}


createWorkflowArgsFromConfig <- function(workflow_name, description = "",
                                         organism_name = NULL, taxon_id = NULL,
                                         source_dir_path = NULL,
                                         final_s4_object = NULL,
                                         contrasts_tbl = NULL,
                                         workflow_data = NULL) {
    # WORKFLOW DETECTION: If workflow_data is NULL, use original RMD logic
    if (is.null(workflow_data)) {
        cat("WORKFLOW ARGS: Detected RMD workflow (workflow_data is NULL), using original logic\n")
        return(.createWorkflowArgsFromConfig_RMD(
            workflow_name = workflow_name,
            description = description,
            organism_name = organism_name,
            taxon_id = taxon_id,
            source_dir_path = source_dir_path
        ))
    }

    # SHINY WORKFLOW: Continue with existing Shiny logic (unchanged)
    cat("WORKFLOW ARGS: Detected Shiny workflow (workflow_data provided), using Shiny logic\n")

    # Validate required inputs
    if (missing(workflow_name) || !is.character(workflow_name) || length(workflow_name) != 1) {
        stop("workflow_name must be a single character string")
    }

    if (is.null(source_dir_path) || !is.character(source_dir_path) || length(source_dir_path) != 1) {
        stop("source_dir_path must be a single character string")
    }

    if (!dir.exists(source_dir_path)) {
        stop("source_dir_path directory does not exist: ", source_dir_path)
    }

    cat("WORKFLOW ARGS: Starting parameter extraction\n")

    # Extract parameters from S4 object @args if provided
    s4_params <- list()
    # Check if final_s4_object has @args slot
    s4_has_args <- tryCatch(
        {
            !is.null(final_s4_object) && isS4(final_s4_object) && !is.null(final_s4_object@args)
        },
        error = function(e) {
            FALSE
        }
    )

    if (s4_has_args) {
        cat("WORKFLOW ARGS: Extracting parameters from S4 @args slot\n")

        tryCatch(
            {
                s4_params <- final_s4_object@args
                cat(sprintf("WORKFLOW ARGS: Found %d function groups in S4 @args\n", length(s4_params)))

                # Log the function groups for debugging
                for (func_name in names(s4_params)) {
                    param_count <- length(s4_params[[func_name]])
                    cat(sprintf("WORKFLOW ARGS: Function '%s' has %d parameters\n", func_name, param_count))
                }
            },
            error = function(e) {
                cat(sprintf("WORKFLOW ARGS: Error extracting S4 params: %s\n", e$message))
                s4_params <- list()
            }
        )
    } else {
        cat("WORKFLOW ARGS: No S4 object provided or S4 @args is NULL\n")
    }

    # Get fallback config_list from global environment
    config_list <- if (exists("config_list", envir = .GlobalEnv)) {
        get("config_list", envir = .GlobalEnv)
    } else {
        list()
    }

    cat(sprintf("WORKFLOW ARGS: Global config_list has %d sections\n", length(config_list)))

    # Extract RUV optimization results from workflow_data if available
    ruv_optimization_result <- NULL
    if (!is.null(workflow_data)) {
        cat("WORKFLOW ARGS: Checking for RUV optimization results in workflow_data\n")

        # Check multiple possible locations for RUV optimization results
        if (!is.null(workflow_data$ruv_optimization_result)) {
            ruv_optimization_result <- workflow_data$ruv_optimization_result
            cat("WORKFLOW ARGS: Found RUV optimization results in workflow_data$ruv_optimization_result\n")
        } else if (!is.null(workflow_data$state_manager)) {
            # Try to get RUV results from state manager config
            tryCatch(
                {
                    current_state_config <- workflow_data$state_manager$getStateConfig(workflow_data$state_manager$current_state)
                    if (!is.null(current_state_config$ruv_optimization_result)) {
                        ruv_optimization_result <- current_state_config$ruv_optimization_result
                        cat("WORKFLOW ARGS: Found RUV optimization results in state manager config\n")
                    }
                },
                error = function(e) {
                    cat(sprintf("WORKFLOW ARGS: Could not extract RUV results from state manager: %s\n", e$message))
                }
            )
        }

        # Also check if there's a saved RUV file
        if (is.null(ruv_optimization_result) && !is.null(source_dir_path)) {
            ruv_file <- file.path(source_dir_path, "ruv_optimization_results.RDS")
            if (file.exists(ruv_file)) {
                tryCatch(
                    {
                        ruv_optimization_result <- readRDS(ruv_file)
                        cat(sprintf("WORKFLOW ARGS: Loaded RUV optimization results from file: %s\n", ruv_file))
                    },
                    error = function(e) {
                        cat(sprintf("WORKFLOW ARGS: Could not load RUV file: %s\n", e$message))
                    }
                )
            }
        }
    }

    # Check for workflow-specific optimization results in R environment
    # These take precedence over workflow_data results
    cat("WORKFLOW ARGS: Checking for workflow-specific optimization results in R environment\n")

    if (workflow_name == "DIA_limpa") {
        cat("WORKFLOW ARGS: Workflow is DIA_limpa, checking for optimal_results_peptides\n")

        # Try to find optimal_results_peptides in R environment
        optimal_results_found <- FALSE

        # Check parent frame first, then global environment
        if (exists("optimal_results_peptides", envir = parent.frame(), inherits = FALSE)) {
            tryCatch(
                {
                    ruv_optimization_result <- get("optimal_results_peptides", envir = parent.frame())
                    optimal_results_found <- TRUE
                    cat("WORKFLOW ARGS: Found optimal_results_peptides in parent frame\n")
                },
                error = function(e) {
                    cat(sprintf("WORKFLOW ARGS: Error getting optimal_results_peptides from parent frame: %s\n", e$message))
                }
            )
        }

        if (!optimal_results_found && exists("optimal_results_peptides", envir = .GlobalEnv)) {
            tryCatch(
                {
                    ruv_optimization_result <- get("optimal_results_peptides", envir = .GlobalEnv)
                    optimal_results_found <- TRUE
                    cat("WORKFLOW ARGS: Found optimal_results_peptides in global environment\n")
                },
                error = function(e) {
                    cat(sprintf("WORKFLOW ARGS: Error getting optimal_results_peptides from global environment: %s\n", e$message))
                }
            )
        }

        if (!optimal_results_found) {
            cat("WORKFLOW ARGS: optimal_results_peptides not found in R environment for DIA_limpa workflow\n")
        }
    } else if (workflow_name == "DIA") {
        cat("WORKFLOW ARGS: Workflow is DIA\n")
        # TODO: Determine object name for DIA workflow optimization results
        # TODO: Implement logic similar to DIA_limpa once object name is determined
        # Example structure:
        # if (exists("optimal_results_OBJECTNAME", envir = .GlobalEnv)) {
        #     ruv_optimization_result <- get("optimal_results_OBJECTNAME", envir = .GlobalEnv)
        # }
    } else if (workflow_name == "TMT") {
        cat("WORKFLOW ARGS: Workflow is TMT\n")
        # TODO: Determine object name for TMT workflow optimization results
        # TODO: Implement logic similar to DIA_limpa once object name is determined
        # Example structure:
        # if (exists("optimal_results_OBJECTNAME", envir = .GlobalEnv)) {
        #     ruv_optimization_result <- get("optimal_results_OBJECTNAME", envir = .GlobalEnv)
        # }
    }

    # Merge S4 parameters with config_list (S4 takes precedence)
    merged_config <- config_list

    if (length(s4_params) > 0) {
        cat("WORKFLOW ARGS: Merging S4 parameters with config_list\n")

        # S4 parameters override config_list parameters
        for (func_name in names(s4_params)) {
            if (is.list(s4_params[[func_name]]) && length(s4_params[[func_name]]) > 0) {
                merged_config[[func_name]] <- s4_params[[func_name]]
                cat(sprintf(
                    "WORKFLOW ARGS: Updated '%s' section with %d S4 parameters\n",
                    func_name, length(s4_params[[func_name]])
                ))
            }
        }
    }

    # Get git information
    git_info <- tryCatch(
        {
            if (requireNamespace("gh", quietly = TRUE)) {
                # Make the API call (works for public repos without token)
                branch_info <- gh::gh("/repos/APAF-BIOINFORMATICS/MultiScholaR/branches/main")
                list(
                    commit_sha = branch_info$commit$sha,
                    branch = "main",
                    repo = "MultiScholaR",
                    timestamp = branch_info$commit$commit$author$date
                )
            } else {
                list(commit_sha = NA_character_, branch = NA_character_, repo = NA_character_, timestamp = NA_character_)
            }
        },
        error = function(e) {
            message("Error fetching GitHub info: ", e$message)
            list(commit_sha = NA_character_, branch = NA_character_, repo = NA_character_, timestamp = NA_character_)
        }
    )

    # Get organism information from session if not provided
    if (is.null(organism_name) && exists("organism_name", envir = .GlobalEnv)) {
        organism_name <- get("organism_name", envir = .GlobalEnv)
    }

    if (is.null(taxon_id) && exists("taxon_id", envir = .GlobalEnv)) {
        taxon_id <- get("taxon_id", envir = .GlobalEnv)
    }

    # Get contrasts_tbl if not provided
    if (is.null(contrasts_tbl)) {
        if (exists("contrasts_tbl", envir = parent.frame())) {
            contrasts_tbl <- get("contrasts_tbl", envir = parent.frame())
        } else if (exists("contrasts_tbl", envir = .GlobalEnv)) {
            contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
        }
    }

    # Build the output content
    output_lines <- c(
        "Study Parameters",
        "================",
        "",
        "Basic Information:",
        "-----------------",
        paste("Workflow Name:", workflow_name),
        paste("Description:", if (nzchar(description)) description else "N/A"),
        paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        ""
    )

    # Add git information
    if (!is.null(git_info) && is.list(git_info)) {
        git_lines <- c(
            "Git Information:",
            "---------------",
            paste("Repository:", ifelse(!is.null(git_info$repo) && !is.na(git_info$repo), git_info$repo, "N/A")),
            paste("Branch:", ifelse(!is.null(git_info$branch) && !is.na(git_info$branch), git_info$branch, "N/A")),
            paste("Commit:", ifelse(!is.null(git_info$commit_sha) && !is.na(git_info$commit_sha), substr(git_info$commit_sha, 1, 7), "N/A")),
            paste("Git Timestamp:", ifelse(!is.null(git_info$timestamp) && !is.na(git_info$timestamp), git_info$timestamp, "N/A")),
            ""
        )
        output_lines <- c(output_lines, git_lines)
    }

    # Add organism information
    if (!is.null(organism_name) && !is.na(organism_name) && nzchar(organism_name)) {
        organism_lines <- c(
            "Organism Information:",
            "---------------------",
            paste("Organism Name:", organism_name)
        )
        if (!is.null(taxon_id) && !is.na(taxon_id) && nzchar(taxon_id)) {
            organism_lines <- c(organism_lines, paste("Taxon ID:", taxon_id))
        }
        organism_lines <- c(organism_lines, "")
        output_lines <- c(output_lines, organism_lines)
    }

    # Add RUV optimization results if available
    if (!is.null(ruv_optimization_result) && is.list(ruv_optimization_result)) {
        cat("WORKFLOW ARGS: Formatting RUV optimization results\n")

        ruv_lines <- c(
            "Automatic RUV Optimization Results:",
            "-----------------------------------"
        )

        # Extract and format RUV optimization values
        tryCatch(
            {
                best_percentage <- if (!is.null(ruv_optimization_result$best_percentage)) {
                    sprintf("%.1f%%", ruv_optimization_result$best_percentage)
                } else {
                    "N/A"
                }

                best_k <- if (!is.null(ruv_optimization_result$best_k)) {
                    as.character(ruv_optimization_result$best_k)
                } else {
                    "N/A"
                }

                separation_score <- if (!is.null(ruv_optimization_result$best_separation_score)) {
                    sprintf("%.4f", ruv_optimization_result$best_separation_score)
                } else {
                    "N/A"
                }

                composite_score <- if (!is.null(ruv_optimization_result$best_composite_score)) {
                    sprintf("%.4f", ruv_optimization_result$best_composite_score)
                } else {
                    "N/A"
                }

                control_genes_count <- if (!is.null(ruv_optimization_result$best_control_genes_index)) {
                    as.character(sum(ruv_optimization_result$best_control_genes_index, na.rm = TRUE))
                } else {
                    "N/A"
                }

                separation_metric <- if (!is.null(ruv_optimization_result$separation_metric_used)) {
                    as.character(ruv_optimization_result$separation_metric_used)
                } else {
                    "N/A"
                }

                k_penalty_weight <- if (!is.null(ruv_optimization_result$k_penalty_weight)) {
                    sprintf("%.1f", ruv_optimization_result$k_penalty_weight)
                } else {
                    "N/A"
                }

                adaptive_penalty <- if (!is.null(ruv_optimization_result$adaptive_k_penalty_used)) {
                    ifelse(ruv_optimization_result$adaptive_k_penalty_used, "TRUE", "FALSE")
                } else {
                    "N/A"
                }

                sample_size <- if (!is.null(ruv_optimization_result$sample_size)) {
                    as.character(ruv_optimization_result$sample_size)
                } else {
                    "N/A"
                }

                # Also try to get RUV grouping variable from S4 parameters
                ruv_grouping_variable <- "N/A"
                if (!is.null(s4_params$ruvIII_C_Varying$ruv_grouping_variable)) {
                    ruv_grouping_variable <- s4_params$ruvIII_C_Varying$ruv_grouping_variable
                } else if (!is.null(s4_params$getNegCtrlProtAnova$ruv_grouping_variable)) {
                    ruv_grouping_variable <- s4_params$getNegCtrlProtAnova$ruv_grouping_variable
                }

                ruv_lines <- c(
                    ruv_lines,
                    paste("• Best percentage:", best_percentage),
                    paste("• Best k value:", best_k),
                    paste("• Separation score:", separation_score),
                    paste("• Composite score:", composite_score),
                    paste("• Control genes:", control_genes_count),
                    paste("• RUV grouping variable:", ruv_grouping_variable),
                    paste("• Separation metric:", separation_metric),
                    paste("• K penalty weight:", k_penalty_weight),
                    paste("• Adaptive penalty:", adaptive_penalty),
                    paste("• Sample size:", sample_size),
                    ""
                )

                cat("WORKFLOW ARGS: Successfully formatted RUV optimization results\n")
            },
            error = function(e) {
                cat(sprintf("WORKFLOW ARGS: Error formatting RUV results: %s\n", e$message))
                ruv_lines <- c(
                    ruv_lines,
                    paste("• [Error formatting RUV optimization results:", e$message, "]"),
                    ""
                )
            }
        )

        output_lines <- c(output_lines, ruv_lines)
    } else {
        cat("WORKFLOW ARGS: No RUV optimization results available\n")
    }

    # Add configuration parameters (now includes S4 parameters)
    config_lines <- c(
        "Workflow Parameters:",
        "-------------------",
        ""
    )

    # Clean the merged config (remove problematic objects)
    clean_config <- merged_config
    # Remove the internal source dir if it exists
    if (!is.null(clean_config$internal_workflow_source_dir)) {
        clean_config$internal_workflow_source_dir <- NULL
    }

    # Add special section for S4-derived parameters
    if (length(s4_params) > 0) {
        cat("WORKFLOW ARGS: About to format S4 parameters using functional approach\n")

        config_lines <- c(
            config_lines,
            "Parameters from Final S4 Object:",
            "--------------------------------"
        )

        # SAFE parameter formatting function
        formatParameterValue <- function(param_value) {
            tryCatch(
                {
                    if (is.null(param_value)) {
                        "NULL"
                    } else if (is.logical(param_value)) {
                        if (length(param_value) == 1) {
                            ifelse(param_value, "TRUE", "FALSE")
                        } else if (length(param_value) > 50) {
                            # Handle large logical vectors (like control genes index)
                            true_count <- sum(param_value, na.rm = TRUE)
                            total_count <- length(param_value)
                            sprintf(
                                "logical vector [%d TRUE, %d FALSE out of %d total]",
                                true_count, total_count - true_count, total_count
                            )
                        } else {
                            # Show first few values for smaller vectors
                            preview <- ifelse(utils::head(param_value, 5), "TRUE", "FALSE")
                            if (length(param_value) > 5) {
                                paste0("c(", paste(preview, collapse = ", "), ", ...)")
                            } else {
                                paste0("c(", paste(preview, collapse = ", "), ")")
                            }
                        }
                    } else if (is.numeric(param_value)) {
                        if (length(param_value) == 1) {
                            as.character(param_value)
                        } else if (length(param_value) > 5) {
                            sprintf(
                                "numeric vector [%d values: %s, ...]",
                                length(param_value),
                                paste(as.character(utils::head(param_value, 3)), collapse = ", ")
                            )
                        } else {
                            paste0("c(", paste(as.character(param_value), collapse = ", "), ")")
                        }
                    } else if (is.character(param_value)) {
                        if (length(param_value) == 1) {
                            param_value
                        } else if (length(param_value) > 5) {
                            sprintf(
                                "character vector [%d values: %s, ...]",
                                length(param_value),
                                paste(shQuote(utils::head(param_value, 3)), collapse = ", ")
                            )
                        } else {
                            paste0("c(", paste(shQuote(param_value), collapse = ", "), ")")
                        }
                    } else {
                        # SAFE fallback - no dput() which was causing the hang
                        paste0("[", class(param_value)[1], " object]")
                    }
                },
                error = function(e) {
                    "[SERIALIZATION ERROR]"
                }
            )
        }

        s4_sections <- if (requireNamespace("purrr", quietly = TRUE)) {
            purrr::imap(s4_params, function(func_params, func_name) {
                tryCatch(
                    {
                        if (!is.list(func_params) || length(func_params) == 0) {
                            return("")
                        }

                        cat(sprintf("WORKFLOW ARGS: Processing S4 function group '%s'\n", func_name))
                        header <- sprintf("[%s]", func_name)

                        # Use map instead of imap_chr for safer handling
                        param_lines <- purrr::imap(func_params, function(param_value, param_name) {
                            tryCatch(
                                {
                                    param_str <- formatParameterValue(param_value)
                                    sprintf("  %s = %s", param_name, param_str)
                                },
                                error = function(e) {
                                    cat(sprintf("WORKFLOW ARGS: Error formatting parameter '%s' in '%s': %s\n", param_name, func_name, e$message))
                                    sprintf("  %s = [ERROR: %s]", param_name, e$message)
                                }
                            )
                        })

                        # Convert list to character vector safely
                        param_lines_char <- unlist(param_lines)
                        paste(c(header, param_lines_char, ""), collapse = "\n")
                    },
                    error = function(e) {
                        cat(sprintf("WORKFLOW ARGS: Error processing S4 function group '%s': %s\n", func_name, e$message))
                        sprintf("[%s]\n  [ERROR: %s]\n", func_name, e$message)
                    }
                )
            })
        } else {
            # Fallback using base R if purrr not available
            lapply(names(s4_params), function(func_name) {
                func_params <- s4_params[[func_name]]
                if (!is.list(func_params) || length(func_params) == 0) {
                    return("")
                }

                header <- sprintf("[%s]", func_name)

                param_lines <- lapply(names(func_params), function(param_name) {
                    param_value <- func_params[[param_name]]
                    param_str <- formatParameterValue(param_value)
                    sprintf("  %s = %s", param_name, param_str)
                })

                paste(c(header, unlist(param_lines), ""), collapse = "\n")
            })
        }

        # Add all S4 sections at once - safely handle list from purrr::imap()
        s4_sections_char <- tryCatch(
            {
                if (is.list(s4_sections)) {
                    # Convert list to character vector
                    unlist(s4_sections)
                } else {
                    s4_sections
                }
            },
            error = function(e) {
                cat(sprintf("WORKFLOW ARGS: Error converting S4 sections: %s\n", e$message))
                "[ERROR: Could not process S4 sections]"
            }
        )

        config_lines <- c(config_lines, unlist(strsplit(paste(s4_sections_char, collapse = ""), "\n")))

        cat("WORKFLOW ARGS: S4 parameters formatted successfully\n")
    }

    # Add UI parameters from workflow_data if available (DE and Enrichment UI inputs)
    if (!is.null(workflow_data)) {
        cat("WORKFLOW ARGS: Checking for UI parameters in workflow_data\n")
        ui_sections <- c()

        # Check for DE UI parameters
        if (!is.null(workflow_data$de_ui_params)) {
            cat("WORKFLOW ARGS: Found DE UI parameters in workflow_data\n")
            de_ui_lines <- c(
                "[Differential Expression UI Parameters]",
                sprintf("  q_value_threshold = %s", ifelse(!is.null(workflow_data$de_ui_params$q_value_threshold), workflow_data$de_ui_params$q_value_threshold, "N/A")),
                sprintf("  log_fold_change_cutoff = %s", ifelse(!is.null(workflow_data$de_ui_params$log_fold_change_cutoff), workflow_data$de_ui_params$log_fold_change_cutoff, "N/A")),
                sprintf("  treat_enabled = %s", ifelse(!is.null(workflow_data$de_ui_params$treat_enabled), workflow_data$de_ui_params$treat_enabled, "N/A")),
                ""
            )
            ui_sections <- c(ui_sections, de_ui_lines)
        }

        # Check for Enrichment UI parameters
        if (!is.null(workflow_data$enrichment_ui_params)) {
            cat("WORKFLOW ARGS: Found Enrichment UI parameters in workflow_data\n")
            enrichment_ui_lines <- c(
                "[Enrichment Analysis UI Parameters]",
                sprintf("  up_log2fc_cutoff = %s", ifelse(!is.null(workflow_data$enrichment_ui_params$up_log2fc_cutoff), workflow_data$enrichment_ui_params$up_log2fc_cutoff, "N/A")),
                sprintf("  down_log2fc_cutoff = %s", ifelse(!is.null(workflow_data$enrichment_ui_params$down_log2fc_cutoff), workflow_data$enrichment_ui_params$down_log2fc_cutoff, "N/A")),
                sprintf("  q_value_cutoff = %s", ifelse(!is.null(workflow_data$enrichment_ui_params$q_value_cutoff), workflow_data$enrichment_ui_params$q_value_cutoff, "N/A")),
                sprintf("  organism_selected = %s", ifelse(!is.null(workflow_data$enrichment_ui_params$organism_selected), workflow_data$enrichment_ui_params$organism_selected, "N/A")),
                sprintf("  database_source = %s", ifelse(!is.null(workflow_data$enrichment_ui_params$database_source), workflow_data$enrichment_ui_params$database_source, "N/A")),
                ""
            )
            ui_sections <- c(ui_sections, enrichment_ui_lines)
        }

        # Add UI sections to config_lines if any were found
        if (length(ui_sections) > 0) {
            config_lines <- c(
                config_lines,
                "User Interface Parameters:",
                "-------------------------",
                ui_sections
            )
            cat("WORKFLOW ARGS: Added UI parameters to output\n")
        } else {
            cat("WORKFLOW ARGS: No UI parameters found in workflow_data\n")
        }
    }

    # Format the remaining config list
    cat("WORKFLOW ARGS: About to format remaining config list\n")
    if (length(clean_config) > 0) {
        config_lines <- c(
            config_lines,
            "Additional Configuration Parameters:",
            "-----------------------------------"
        )

        cat("WORKFLOW ARGS: Calling formatConfigList...\n")
        tryCatch(
            {
                # Check if formatConfigList exists before calling
                if (exists("formatConfigList", mode = "function")) {
                    config_params <- formatConfigList(clean_config)
                    config_lines <- c(config_lines, config_params)
                    cat("WORKFLOW ARGS: formatConfigList completed successfully\n")
                } else {
                    cat("WORKFLOW ARGS: formatConfigList function not found, using basic formatting\n")
                    # Basic fallback formatting without for loops
                    basic_config <- unlist(clean_config, recursive = TRUE)
                    config_lines <- c(config_lines, names(basic_config), " = ", as.character(basic_config))
                }
            },
            error = function(e) {
                cat(sprintf("WORKFLOW ARGS: formatConfigList failed: %s\n", e$message))
                config_lines <- c(config_lines, paste("Error formatting config:", e$message))
            }
        )
    } else {
        config_lines <- c(config_lines, "No additional configuration parameters available")
    }

    cat("WORKFLOW ARGS: Adding config lines to output\n")
    output_lines <- c(output_lines, config_lines)

    # Add contrasts information
    cat("WORKFLOW ARGS: Processing contrasts information\n")
    if (!is.null(contrasts_tbl) && (is.data.frame(contrasts_tbl) || tibble::is_tibble(contrasts_tbl)) && nrow(contrasts_tbl) > 0) {
        cat("WORKFLOW ARGS: Adding contrasts to output\n")
        contrasts_lines <- c(
            "",
            "Contrasts:",
            "----------"
        )

        if ("contrasts" %in% colnames(contrasts_tbl)) {
            contrasts_info <- tryCatch(
                {
                    contrasts_col <- contrasts_tbl[["contrasts"]]
                    paste("  ", as.character(contrasts_col))
                },
                error = function(e) {
                    paste("  [Error extracting contrasts:", e$message, "]")
                }
            )
            contrasts_lines <- c(contrasts_lines, contrasts_info)
        } else {
            contrasts_lines <- c(contrasts_lines, "  [Column 'contrasts' not found in contrasts_tbl]")
        }

        output_lines <- c(output_lines, contrasts_lines)
        cat("WORKFLOW ARGS: Contrasts added successfully\n")
    } else {
        cat("WORKFLOW ARGS: No contrasts to add\n")
    }

    # Write to file
    cat("WORKFLOW ARGS: About to write file\n")
    output_file <- file.path(source_dir_path, "study_parameters.txt")
    cat(sprintf("WORKFLOW ARGS: Target file path: %s\n", output_file))

    tryCatch(
        {
            cat("WORKFLOW ARGS: Calling writeLines...\n")
            writeLines(output_lines, output_file)
            cat(sprintf("WORKFLOW ARGS: Study parameters saved to: %s\n", output_file))
            return(output_file)
        },
        error = function(e) {
            cat(sprintf("WORKFLOW ARGS: Error writing file: %s\n", e$message))
            stop("Failed to write study parameters file: ", e$message)
        }
    )
}

#' Check Missing Value Percentages in Peptide Data
#'
#' @description Calculate and report the percentage of missing values (NAs) in peptide data
#' at different levels: total dataset, per sample, and per group.
#'
#' @param peptide_obj A PeptideQuantitativeData S4 object
#' @param verbose Logical, whether to print detailed results (default: TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item total_na_percent: Overall percentage of NAs in the dataset
#'   \item per_sample_na: Data frame with NA percentages per sample
#'   \item per_group_na: Data frame with NA percentages per group
#'   \item summary_stats: Summary statistics of NA distribution
#' }
#'
#' @export
checkPeptideNAPercentages <- function(peptide_obj, verbose = TRUE) {
    # Validate input
    if (!is(peptide_obj, "PeptideQuantitativeData")) {
        stop("Input must be a PeptideQuantitativeData S4 object")
    }

    # Extract data from S4 object
    peptide_matrix <- peptide_obj@peptide_matrix
    design_matrix <- peptide_obj@design_matrix
    sample_id_col <- peptide_obj@sample_id
    group_id_col <- peptide_obj@group_id

    # Validate that matrix and design matrix are compatible
    if (ncol(peptide_matrix) != nrow(design_matrix)) {
        stop("Number of samples in peptide_matrix doesn't match design_matrix rows")
    }

    # Calculate total NA percentage
    total_values <- length(peptide_matrix)
    total_nas <- sum(is.na(peptide_matrix))
    total_na_percent <- (total_nas / total_values) * 100

    # Calculate per-sample NA percentages
    sample_na_counts <- apply(peptide_matrix, 2, function(x) sum(is.na(x)))
    sample_na_percentages <- (sample_na_counts / nrow(peptide_matrix)) * 100

    per_sample_na <- data.frame(
        sample = colnames(peptide_matrix),
        na_count = sample_na_counts,
        na_percentage = sample_na_percentages,
        stringsAsFactors = FALSE
    )

    # Add group information to per-sample results
    per_sample_na <- merge(per_sample_na, design_matrix,
        by.x = "sample", by.y = sample_id_col, all.x = TRUE
    )

    # Calculate per-group NA percentages
    per_group_na <- per_sample_na %>%
        group_by(!!sym(group_id_col)) %>%
        summarise(
            num_samples = n(),
            mean_na_percentage = mean(na_percentage, na.rm = TRUE),
            median_na_percentage = median(na_percentage, na.rm = TRUE),
            min_na_percentage = min(na_percentage, na.rm = TRUE),
            max_na_percentage = max(na_percentage, na.rm = TRUE),
            sd_na_percentage = sd(na_percentage, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        arrange(mean_na_percentage)

    # Calculate summary statistics
    summary_stats <- list(
        total_peptides = nrow(peptide_matrix),
        total_samples = ncol(peptide_matrix),
        total_groups = length(unique(design_matrix[[group_id_col]])),
        total_values = total_values,
        total_nas = total_nas,
        mean_na_per_sample = mean(sample_na_percentages),
        median_na_per_sample = median(sample_na_percentages),
        min_na_per_sample = min(sample_na_percentages),
        max_na_per_sample = max(sample_na_percentages)
    )

    # Print results if verbose
    if (verbose) {
        cat("\n=== Peptide Data Missing Value Analysis ===\n")
        cat(sprintf(
            "Dataset dimensions: %d peptides × %d samples\n",
            nrow(peptide_matrix), ncol(peptide_matrix)
        ))
        cat(sprintf("Number of groups: %d\n", summary_stats$total_groups))
        cat(sprintf(
            "Total missing values: %s out of %s (%.2f%%)\n",
            format(total_nas, big.mark = ","),
            format(total_values, big.mark = ","),
            total_na_percent
        ))

        cat("\n--- Per-Sample Missing Value Summary ---\n")
        cat(sprintf("Mean NA%% per sample: %.2f%%\n", summary_stats$mean_na_per_sample))
        cat(sprintf("Median NA%% per sample: %.2f%%\n", summary_stats$median_na_per_sample))
        cat(sprintf(
            "Range: %.2f%% - %.2f%%\n",
            summary_stats$min_na_per_sample, summary_stats$max_na_per_sample
        ))

        cat("\n--- Per-Group Missing Value Summary ---\n")
        print(per_group_na)

        cat("\n--- Samples with Highest Missing Values ---\n")
        top_missing_samples <- per_sample_na %>%
            arrange(desc(na_percentage)) %>%
            head(min(5, nrow(per_sample_na)))
        print(top_missing_samples[, c("sample", group_id_col, "na_percentage")])

        cat("\n--- Samples with Lowest Missing Values ---\n")
        bottom_missing_samples <- per_sample_na %>%
            arrange(na_percentage) %>%
            head(min(5, nrow(per_sample_na)))
        print(bottom_missing_samples[, c("sample", group_id_col, "na_percentage")])
    }

    # Return results
    results <- list(
        total_na_percent = total_na_percent,
        per_sample_na = per_sample_na,
        per_group_na = per_group_na,
        summary_stats = summary_stats
    )

    return(invisible(results))
}

#' Validate Post-Imputation Peptide Data
#'
#' @description A simple wrapper to validate peptide data after imputation,
#' specifically checking if imputation was successful (should show 0% NAs).
#'
#' @param peptide_obj A PeptideQuantitativeData S4 object (post-imputation)
#' @param expected_na_percent Expected NA percentage (default: 0 for post-imputation)
#' @param tolerance Tolerance for expected percentage (default: 0.1%)
#'
#' @return Logical indicating if validation passed, with detailed output
#'
#' @export
validatePostImputationData <- function(peptide_obj, expected_na_percent = 0, tolerance = 0.1) {
    cat("\n=== POST-IMPUTATION VALIDATION ===\n")

    # Run the full NA analysis
    na_results <- checkPeptideNAPercentages(peptide_obj, verbose = TRUE)

    # Check if imputation was successful
    actual_na_percent <- na_results$total_na_percent
    is_valid <- abs(actual_na_percent - expected_na_percent) <= tolerance

    cat("\n--- VALIDATION RESULT ---\n")
    cat(sprintf("Expected NA%%: %.2f%% (± %.2f%%)\n", expected_na_percent, tolerance))
    cat(sprintf("Actual NA%%: %.2f%%\n", actual_na_percent))

    if (is_valid) {
        cat("✓ VALIDATION PASSED: Imputation appears successful!\n")
    } else {
        cat("✗ VALIDATION FAILED: Unexpected NA percentage detected!\n")
        if (actual_na_percent > expected_na_percent + tolerance) {
            cat("  → Issue: More NAs than expected. Imputation may have failed.\n")
        } else {
            cat("  → Issue: Fewer NAs than expected. Check data integrity.\n")
        }
    }

    # Additional warnings for common issues
    if (actual_na_percent > 10) {
        cat("⚠ WARNING: High NA percentage suggests imputation problems!\n")
    }

    if (na_results$summary_stats$max_na_per_sample > actual_na_percent + 5) {
        cat("⚠ WARNING: Large variation in NA% between samples detected!\n")
    }

    cat("\n")
    return(invisible(list(
        is_valid = is_valid,
        actual_na_percent = actual_na_percent,
        expected_na_percent = expected_na_percent,
        full_results = na_results
    )))
}

#' Get Recommendations for Handling Protein-Level Missing Values
#'
#' @description Provides specific recommendations for dealing with missing values
#' in protein data based on the percentage and distribution of NAs.
#'
#' @param protein_obj A ProteinQuantitativeData S4 object
#' @param include_code Logical, whether to include example R code (default: TRUE)
#'
#' @return Prints recommendations and invisibly returns a list of strategies
#'
#' @export
getProteinNARecommendations <- function(protein_obj, include_code = TRUE) {
    # Get NA analysis
    na_results <- checkProteinNAPercentages(protein_obj, verbose = FALSE)
    na_percent <- na_results$total_na_percent

    cat("\n=== PROTEIN NA HANDLING RECOMMENDATIONS ===\n")
    cat(sprintf(
        "Your data: %.1f%% NAs across %d proteins\n\n",
        na_percent, na_results$summary_stats$total_proteins
    ))

    if (na_percent < 15) {
        cat("🎯 RECOMMENDATION: Complete Case Analysis\n")
        cat("• Your data has excellent protein coverage\n")
        cat("• Can proceed with standard analysis on proteins with complete data\n")
        if (include_code) {
            cat("\n📝 Example code:\n")
            cat("complete_proteins <- protein_obj@protein_quant_table[complete.cases(protein_obj@protein_quant_table), ]\n")
        }
    } else if (na_percent >= 15 && na_percent < 40) {
        cat("🎯 RECOMMENDATION: Consider Protein-Level Imputation\n")
        cat("• Moderate missing values - imputation could be beneficial\n")
        cat("• Options: KNN, minimum value, or mixed imputation strategies\n")
        cat("• Alternative: Filter to proteins detected in ≥X samples per group\n")
        if (include_code) {
            cat("\n📝 Example filtering code:\n")
            cat("# Keep proteins detected in ≥50% of samples per group\n")
            cat("filtered_proteins <- filterProteinsByGroupDetection(protein_obj, min_detection_rate = 0.5)\n")
        }
    } else if (na_percent >= 40 && na_percent < 60) {
        cat("🎯 RECOMMENDATION: Strict Filtering + Targeted Imputation\n")
        cat("• High missing values suggest challenging sample/detection conditions\n")
        cat("• Focus on well-detected proteins (present in majority of samples)\n")
        cat("• Consider group-wise detection requirements\n")
        if (include_code) {
            cat("\n📝 Example approach:\n")
            cat("# Keep proteins detected in ≥70% of samples in at least one group\n")
            cat("robust_proteins <- filterProteinsByGroupwise(protein_obj, min_group_detection = 0.7)\n")
        }
    } else {
        cat("⚠️  RECOMMENDATION: Review Data Quality\n")
        cat("• Very high missing values (>60%) suggest potential issues\n")
        cat("• Check: sample quality, peptide identification, rollup parameters\n")
        cat("• Consider more stringent protein identification criteria\n")
        cat("• May need to focus only on highly abundant/well-detected proteins\n")
    }

    cat("\n📚 STRATEGIES SUMMARY:\n")
    cat("1. Complete Case: Use only proteins with no NAs\n")
    cat("2. Filtering: Remove proteins with >X% missing values\n")
    cat("3. Group-wise: Require detection in ≥Y% samples per group\n")
    cat("4. Imputation: Fill NAs with estimated values (KNN, minimum, etc.)\n")
    cat("5. Hybrid: Combine filtering + imputation\n")

    cat("\n💡 TIP: Protein NAs ≠ Data Quality Issues\n")
    cat("Missing proteins often reflect:\n")
    cat("• Low abundance proteins below detection limit\n")
    cat("• Sample-specific biology (some proteins not expressed)\n")
    cat("• Normal variation in complex proteomes\n\n")

    strategies <- list(
        na_percent = na_percent,
        primary_recommendation = if (na_percent < 15) {
            "complete_case"
        } else if (na_percent < 40) {
            "imputation_or_filtering"
        } else if (na_percent < 60) {
            "strict_filtering"
        } else {
            "data_quality_review"
        },
        alternative_strategies = c("complete_case", "group_wise_filtering", "imputation", "hybrid")
    )

    return(invisible(strategies))
}

#' Check Missing Value Percentages in Protein Data
#'
#' @description Calculate and report the percentage of missing values (NAs) in protein data
#' at different levels: total dataset, per sample, and per group.
#'
#' @param protein_obj A ProteinQuantitativeData S4 object
#' @param verbose Logical, whether to print detailed results (default: TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item total_na_percent: Overall percentage of NAs in the dataset
#'   \item per_sample_na: Data frame with NA percentages per sample
#'   \item per_group_na: Data frame with NA percentages per group
#'   \item summary_stats: Summary statistics of NA distribution
#' }
#'
#' @export
checkProteinNAPercentages <- function(protein_obj, verbose = TRUE) {
    # Validate input
    if (!is(protein_obj, "ProteinQuantitativeData")) {
        stop("Input must be a ProteinQuantitativeData S4 object")
    }

    # Extract data from S4 object
    protein_quant_table <- protein_obj@protein_quant_table
    design_matrix <- protein_obj@design_matrix
    sample_id_col <- protein_obj@sample_id
    group_id_col <- protein_obj@group_id
    protein_id_col <- protein_obj@protein_id_column

    # Identify sample columns (exclude protein ID column)
    sample_columns <- setdiff(colnames(protein_quant_table), protein_id_col)

    # Validate that sample columns match design matrix
    if (length(sample_columns) != nrow(design_matrix)) {
        stop("Number of sample columns doesn't match design_matrix rows")
    }

    # Extract quantitative data matrix (samples only)
    protein_matrix <- as.matrix(protein_quant_table[, sample_columns])
    rownames(protein_matrix) <- protein_quant_table[[protein_id_col]]

    # Calculate total NA percentage
    total_values <- length(protein_matrix)
    total_nas <- sum(is.na(protein_matrix))
    total_na_percent <- (total_nas / total_values) * 100

    # Calculate per-sample NA percentages
    sample_na_counts <- apply(protein_matrix, 2, function(x) sum(is.na(x)))
    sample_na_percentages <- (sample_na_counts / nrow(protein_matrix)) * 100

    per_sample_na <- data.frame(
        sample = names(sample_na_counts),
        na_count = sample_na_counts,
        na_percentage = sample_na_percentages,
        stringsAsFactors = FALSE
    )

    # Add group information to per-sample results
    per_sample_na <- merge(per_sample_na, design_matrix,
        by.x = "sample", by.y = sample_id_col, all.x = TRUE
    )

    # Calculate per-group NA percentages
    per_group_na <- per_sample_na %>%
        group_by(!!sym(group_id_col)) %>%
        summarise(
            num_samples = n(),
            mean_na_percentage = mean(na_percentage, na.rm = TRUE),
            median_na_percentage = median(na_percentage, na.rm = TRUE),
            min_na_percentage = min(na_percentage, na.rm = TRUE),
            max_na_percentage = max(na_percentage, na.rm = TRUE),
            sd_na_percentage = sd(na_percentage, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        arrange(mean_na_percentage)

    # Calculate summary statistics
    summary_stats <- list(
        total_proteins = nrow(protein_matrix),
        total_samples = ncol(protein_matrix),
        total_groups = length(unique(design_matrix[[group_id_col]])),
        total_values = total_values,
        total_nas = total_nas,
        mean_na_per_sample = mean(sample_na_percentages),
        median_na_per_sample = median(sample_na_percentages),
        min_na_per_sample = min(sample_na_percentages),
        max_na_per_sample = max(sample_na_percentages)
    )

    # Print results if verbose
    if (verbose) {
        cat("\n=== Protein Data Missing Value Analysis ===\n")
        cat(sprintf(
            "Dataset dimensions: %d proteins × %d samples\n",
            nrow(protein_matrix), ncol(protein_matrix)
        ))
        cat(sprintf("Number of groups: %d\n", summary_stats$total_groups))
        cat(sprintf(
            "Total missing values: %s out of %s (%.2f%%)\n",
            format(total_nas, big.mark = ","),
            format(total_values, big.mark = ","),
            total_na_percent
        ))

        cat("\n--- Per-Sample Missing Value Summary ---\n")
        cat(sprintf("Mean NA%% per sample: %.2f%%\n", summary_stats$mean_na_per_sample))
        cat(sprintf("Median NA%% per sample: %.2f%%\n", summary_stats$median_na_per_sample))
        cat(sprintf(
            "Range: %.2f%% - %.2f%%\n",
            summary_stats$min_na_per_sample, summary_stats$max_na_per_sample
        ))

        cat("\n--- Per-Group Missing Value Summary ---\n")
        print(per_group_na)

        cat("\n--- Samples with Highest Missing Values ---\n")
        top_missing_samples <- per_sample_na %>%
            arrange(desc(na_percentage)) %>%
            head(min(5, nrow(per_sample_na)))
        print(top_missing_samples[, c("sample", group_id_col, "na_percentage")])

        cat("\n--- Samples with Lowest Missing Values ---\n")
        bottom_missing_samples <- per_sample_na %>%
            arrange(na_percentage) %>%
            head(min(5, nrow(per_sample_na)))
        print(bottom_missing_samples[, c("sample", group_id_col, "na_percentage")])
    }

    # Return results
    results <- list(
        total_na_percent = total_na_percent,
        per_sample_na = per_sample_na,
        per_group_na = per_group_na,
        summary_stats = summary_stats
    )

    return(invisible(results))
}

#' Validate Post-Imputation Protein Data
#'
#' @description A simple wrapper to validate protein data after imputation,
#' specifically checking if imputation was successful.
#'
#' @param protein_obj A ProteinQuantitativeData S4 object (post-imputation)
#' @param expected_na_percent Expected NA percentage (default: varies based on protein data)
#' @param tolerance Tolerance for expected percentage (default: 10%)
#'
#' @return Logical indicating if validation passed, with detailed output
#'
#' @export
validatePostImputationProteinData <- function(protein_obj, expected_na_percent = NULL, tolerance = 10) {
    cat("\n=== POST-IMPUTATION PROTEIN DATA VALIDATION ===\n")
    cat("Note: Protein-level NAs occur even after peptide imputation because:\n")
    cat("• Proteins need ≥1 detected peptide to get a quantification\n")
    cat("• Some proteins detected only in subset of samples\n")
    cat("• This is normal proteomics data behavior!\n\n")

    # Run the full NA analysis
    na_results <- checkProteinNAPercentages(protein_obj, verbose = TRUE)

    # Set expected NA percentage if not provided (proteins often have some NAs)
    if (is.null(expected_na_percent)) {
        # For protein data, NAs are very common due to missing peptides/proteins
        # Typical ranges: 20-50% depending on sample complexity and detection method
        expected_na_percent <- 35 # Realistic expectation for protein data
        cat(sprintf("Note: Using default expected NA%% of %.1f%% for protein data\n", expected_na_percent))
        cat("(Protein-level NAs are normal due to incomplete protein detection across samples)\n")
    }

    # Check if validation passes
    actual_na_percent <- na_results$total_na_percent
    is_valid <- abs(actual_na_percent - expected_na_percent) <= tolerance

    cat("\n--- VALIDATION RESULT ---\n")
    cat(sprintf("Expected NA%%: %.2f%% (± %.2f%%)\n", expected_na_percent, tolerance))
    cat(sprintf("Actual NA%%: %.2f%%\n", actual_na_percent))

    if (is_valid) {
        cat("✓ VALIDATION PASSED: Protein data NA levels are within expected range!\n")
    } else {
        cat("✗ VALIDATION FAILED: Unexpected NA percentage detected!\n")
        if (actual_na_percent > expected_na_percent + tolerance) {
            cat("  → Issue: More NAs than expected. Check for missing proteins/peptides.\n")
        } else {
            cat("  → Issue: Fewer NAs than expected. Possible over-imputation.\n")
        }
    }

    # Additional warnings for common issues
    if (actual_na_percent > 50) {
        cat("⚠ WARNING: Very high NA percentage (>50%) suggests data quality issues!\n")
    }

    if (actual_na_percent < 10) {
        cat("ℹ INFO: Very low NA percentage (<10%) - excellent protein coverage!\n")
    }

    # Educational information about protein NAs
    if (actual_na_percent > 20 && actual_na_percent < 50) {
        cat("ℹ INFO: NA percentage is typical for protein-level data\n")
        cat("  → This reflects biological reality: not all proteins detected in all samples\n")
        cat("  → Consider: protein-level imputation OR complete-case analysis\n")
    }

    if (na_results$summary_stats$max_na_per_sample > actual_na_percent + 10) {
        cat("⚠ WARNING: Large variation in NA% between samples detected!\n")
        cat("  → Some samples may have much lower protein coverage.\n")
    }

    # Check for problematic samples (>80% missing)
    high_missing_samples <- na_results$per_sample_na[na_results$per_sample_na$na_percentage > 80, ]
    if (nrow(high_missing_samples) > 0) {
        cat("⚠ WARNING: Samples with >80% missing proteins detected:\n")
        print(high_missing_samples[, c("sample", "na_percentage")])
    }

    cat("\n")
    return(invisible(list(
        is_valid = is_valid,
        actual_na_percent = actual_na_percent,
        expected_na_percent = expected_na_percent,
        full_results = na_results
    )))
}
