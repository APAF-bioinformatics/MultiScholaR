#!/usr/bin/env Rscript

# Author(s): Ignatius Pang, Pablo Galaviz
# Email: cmri-bioinformatics@cmri.org.au
# Children’s Medical Research Institute, finding cures for childhood genetic diseases


## Packages installation and loading

#Test if BioManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# load pacman package manager
if(!require(pacman)){
    install.packages("pacman")
    library(pacman)
}


## Load packages
p_load(tidyverse)
p_load(plotly)
p_load(vroom)
p_load(ggplot2)
p_load(ggpubr)
p_load(writexl)
p_load(magrittr)
p_load(knitr)
p_load(rlang)
p_load(limma)
p_load(qvalue)
p_load(ReactomeGSA)
p_load(ProteomeRiver)
p_load(configr)
p_load(logging)
p_load(optparse)
p_load(tictoc)

## Parameters
tic()

## Directories management
base_dir <- here::here()

# output_dir <- "/home/ignatius/PostDoc/2022/PYROXD1_BMP_13/Results/N155S/enrich_proteins"
# group_pattern <- "\\d+"
# de_proteins_file <- "/home/ignatius/PostDoc/2022/PYROXD1_BMP_13/Results/N155S/de_proteins/de_proteins_long.tsv"
# contrasts_file <- "/home/ignatius/PostDoc/2022/PYROXD1_BMP_13/Source/N155S/contrast_strings.tab"
# design_matrix_file <- "/home/ignatius/PostDoc/2022/PYROXD1_BMP_13/Source/N155S/design_matrix.tab"
# counts_table_file <- "/home/ignatius/PostDoc/2022/PYROXD1_BMP_13/Results/N155S/de_proteins/normalized_counts_after_ruv.tsv"
# formula_string <- "~ 0 + group "
# sample_id <- "Sample_ID"
# group_id <- "group"
# row_id <- "uniprot_acc"

command_line_options <- commandArgs(trailingOnly = TRUE)

  parser <- OptionParser(add_help_option =TRUE)

  parser <- add_option(parser, c("-c", "--config"), type = "character",  default = "/home/ubuntu/Workings/2022/Herpes_Neuropathogenesis_BMP_15/Source/config_prot.ini",
                       help = "Configuration file.  [default %default]",
                       metavar = "string")

  parser <- add_option(parser, c("-l", "--log_file"), type = "character", default = "output.log",
                       help = "Name of the logging file.  [default %default]",
                       metavar = "string")

  parser <- add_option(parser, c("-n", "--no_backup"), action = "store_true", default = FALSE,
                       help = "Deactivate backup of previous run.  [default %default]")

  parser <- add_option(parser, c("-t","--tmp_dir"), type = "character", default = "cache", dest = "tmp_dir",
                       help = "Directory path for temporary files.",
                       metavar = "string")

  parser <- add_option(parser, c( "--enrichment_method"), type="character",  dest = "enrichment_method",
                       help="Enrichment methods used: 'PADOG', 'Camera', or 'ssGSEA'",
                       metavar="string")

  parser <- add_option(parser, c( "--max_missing_values"), type="double",  dest = "max_missing_values",
                       help="Enrichment methods used: 'PADOG', 'Camera', or 'ssGSEA'",
                       metavar="string")

  parser <- add_option(parser, c( "--use_interactors"), type="logical",  dest = "use_interactors",
                       help="Indicates whether interactors from IntAct should be used to extent REACTOME's pathways in the analysis.",
                       metavar="string")

  parser <- add_option(parser, c( "--include_disease_pathways"), type="logical",  dest = "include_disease_pathways",
                       help="include_disease_pathways.",
                       metavar="string")

  parser <- add_option(parser, c( "--email"), type="character",  dest = "email",
                       help="If set to a valid e-mail address, links to the analysis result (and report) will be sent once the analysis is complete.",
                       metavar="string")

  parser <- add_option(parser, c( "--group-pattern"), type="character", dest = "group_pattern",
                       help="Regular expression pattern to identify columns with abundance values belonging to the experiment. [default %default]",
                       metavar="string")

  parser <- add_option(parser, c( "--counts"), type="character", dest = "counts_table_file",
                       help="Input file with the protein abundance values",
                       metavar="string")

  # parser <- add_option(parser, c( "--de-proteins"), type="character", dest = "de_proteins_file",
  #                      help="Input file with the list of diffierentiall expressed protein log fold-change and q-values for every contrasts.",
  #                      metavar="string")

  parser <- add_option(parser, c("--contrasts"), type="character", dest = "contrasts_file",
                       help="Input file with a table listing all comparisons to be made in string, one comparison per line (e.g. groupB.vs.group_A = groupB - groupA).",
                       metavar="string")

  parser <- add_option(parser, c("--formula"), type="character", dest = "formula_string",
                       help="A string representing the formula for input into the model.frame function. (e.g. ~ 0 + group).",
                       metavar="string")

  parser <- add_option(parser, c( "--design-matrix"), type="character", dest = "design_matrix_file",
                       help="Input file with the design matrix",
                       metavar="string")

  parser <- add_option(parser, c("-o", "--output-dir"), type="character", dest = "output_dir",
                       help="Directory path for all results files.",
                       metavar="string")

  parser <- add_option(parser, c("--sample-id"), type="character",dest = "sample_id",
                       help="A string describing the sample ID. This must be a column that exists in the design matrix.",
                       metavar="string")

  parser <- add_option(parser, c( "--group-id"), type="character", dest = "group_id",
                       help="A string describing the experimental group ID. This must be a column that exists in the design matrix.",
                       metavar="string")

  parser <- add_option(parser, c( "--row-id"), type="character", dest = "row_id",
                       help="A string describing the row id.",
                       metavar="string")

  print(commandArgs(trailingOnly = TRUE))

  args <- parse_args(parser)

  #parse and merge the configuration file options.
  if (args$config != "") {
    args <- config.list.merge(eval.config(file = args$config, config="enrich_reactome_proteins"), args)
  }

  args <- setArgsDefault(args, "output_dir", as_func=as.character, default_val="enrich_reactome_proteins" )

  createOutputDir(args$output_dir, args$no_backup)
  createDirectoryIfNotExists(args$tmp_dir)

  ## Logger configuration
  logReset()
  addHandler(writeToConsole , formatter = cmriFormatter)
  addHandler(writeToFile, file = file.path(args$output_dir, args$log_file), formatter = cmriFormatter)

  level <- ifelse(args$debug, loglevels["DEBUG"], loglevels["INFO"])
  setLevel(level = ifelse(args$silent, loglevels["ERROR"], level))

  cmriWelcome("ProteomeRiver", c("Ignatius Pang", "Pablo Galaviz"))
  loginfo("Reading configuration file %s", args$config)
  loginfo("Argument: Value")
  loginfo("----------------------------------------------------")
  for (v in names(args))
  {
    loginfo("%s : %s", v, args[v])
  }
  loginfo("----------------------------------------------------")

  print(args)

  testRequiredArguments(args, c(
    "formula_string",
    "group_pattern",
    "output_dir",
    "sample_id",
    "group_id",
    "row_id"
  ))

  testRequiredFiles(c(
    args$design_matrix_file,
    args$counts_table_file,
    # args$de_proteins_file,
    args$contrasts_file,
    args$config,
    args$tmp_dir
  ))

  args <- setArgsDefault(args, "enrichment_method", as_func=as.character, default_val="Camera" )
  args <- setArgsDefault(args, "max_missing_values", as_func=as.double, default_val=0.5 )
  args <- setArgsDefault(args, "use_interactors", as_func=as.logical, default_val=FALSE )
  args <- setArgsDefault(args, "include_disease_pathways", as_func=as.logical, default_val=TRUE )
  args <- setArgsDefault(args, "email", as_func=as.character, default_val="" )

  args<-parseString(args,
                    c(  "formula_string",
                        "group_pattern",
                        "output_dir",
                        "sample_id",
                        "group_id",
                        "row_id",
                        "tmp_dir",
                        "email",
                        "enrichment_method" ))

  args<-parseType(args,
                  c("max_missing_values"),
                  as.double)

  args<-parseType(args,
                  c("use_interactors",
                    "include_disease_pathways"),
                  as.logical)

  group_pattern <- args$group_pattern
  design_matrix_file <- args$design_matrix_file
  counts_table_file <- args$counts_table_file
  # de_proteins_file <- args$de_proteins_file
  contrasts_file <- args$contrasts_file
  output_dir <- args$output_dir
  sample_id <- args$sample_id
  group_id <-  args$group_id
  row_id <- args$row_id
  formula_string <- args$formula_string

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Create directories
createDirectoryIfNotExists( output_dir)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  loginfo("Read design matrix file.")

  captured_output <- capture.output(
    design_mat_cln <- vroom::vroom(args$design_matrix_file) %>%
      as.data.frame() %>%
      dplyr::mutate(!!rlang::sym(args$sample_id) := as.character(!!rlang::sym(args$sample_id)))
    ,type = "message"
  )
  logdebug(captured_output)

  rownames(design_mat_cln) <- design_mat_cln %>% pull(as.name(args$sample_id))

  ## Check that the sample ID and group ID name is found in the design matrix
  if (length(which(c(args$sample_id, args$group_id) %in% colnames(design_mat_cln))) != 2) {
    logerror("Sample ID and group ID are not matching to the column names used in the design matrix.")
    q()
  }

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Read the Contrasts file
contrasts_tbl <- NA
if( contrasts_file != "") {
  loginfo("Read file with lists of experimental contrasts to test %s", args$contrasts_file)
  captured_output<-capture.output(
    contrasts_tbl <- vroom::vroom( contrasts_file, delim="\t"),
    type = "message"
  )
  logdebug(captured_output)
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## List of read abundance tables after normalization with RUVIII
loginfo("Read abundance tables after normalization with RUVIII.")
captured_output <- capture.output(
norm_abundance <- vroom::vroom( counts_table_file )  %>%
  mutate( best_uniprot_acc = purrr::map_chr( uniprot_acc, ~str_split(. , ":"   )[[1]][1] )),
type = "message"
)
logdebug(captured_output)

norm_abundance_mat <- norm_abundance %>%
  dplyr::select(best_uniprot_acc, matches( group_pattern ) ) %>%
  column_to_rownames("best_uniprot_acc")


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Functions to create the replicate matrix
loginfo("Create the replicate matrix.")

ruvIII_replicates_matrix <- getRuvIIIReplicateMatrix( design_mat_cln,
                                                         !!rlang::sym(sample_id),
                                                         !!rlang::sym(group_id))

ruvIII_replicates_matrix

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Set up list of contrasts

loginfo("Create the list of contrasts.")

ff <- as.formula(  formula_string)
mod_frame <- model.frame( ff, design_mat_cln)
design_m <- model.matrix( ff, mod_frame)

contr.matrix <- makeContrasts( contrasts = contrasts_tbl %>% pull(contrasts),
                                 levels = colnames(design_m))

# ## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ## List of all proteins
# 
# loginfo("Get the list of all DE proteins.")
# captured_output <- capture.output(
# proteins_cln <- vroom::vroom(de_proteins_file) %>%
#   mutate( best_uniprot_acc = purrr::map_chr( uniprot_acc, ~str_split(. , ":"   )[[1]][1] )) %>%
#   dplyr::select( -uniprot_acc ) %>%
#   dplyr::rename( uniprot_acc = "best_uniprot_acc")
# ,type = "message"
# )
# logdebug(captured_output)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## List of contrasts
loginfo("Prepare gene ID to gene set dictionary.")

contrasts_tbl <- vroom::vroom(contrasts_file, delim = "\t")
contrast_strings <- contrasts_tbl[, 1][[1]]
  ## Make contrasts
  contr.matrix <- makeContrasts(contrasts = contrast_strings,
                                levels = colnames(design_m))

lists_of_contrasts <-   purrr::map(colnames(contr.matrix), ~contr.matrix[,.] )
names(lists_of_contrasts) <-  colnames(contr.matrix) %>% str_split( "=") %>% purrr::map_chr(1)

# lists_of_contrasts <- list( RPE_Y.vs.RPE_X=c(0, 0, -1, 1 ) ,
#                             RPE_B.vs.RPE_A=c(1, -1, 0, 0 ) ,
#                             RPE_B.vs.RPE_X=c(1, 0, -1, 0 ) ,
#                             RPE_Y.vs.RPE_A=c(0, -1, 0, 1 ) ,
#                             RPE_A.vs.RPE_X=c(0, 1, -1, 0 ) ,
#                             RPE_B.vs.RPE_Y=c(1, 0, 0, -1 )   )

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Getting available enrichment methods.")

available_methods <- get_reactome_methods(print_methods = FALSE, return_result = TRUE)

# only show the names of the available methods
available_methods$name
#> [1] "PADOG"  "Camera" "ssGSEA"
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("To get more information about a specific method, set `print_details` to `TRUE` and specify the `method`.")

method_params <- available_methods$parameters[available_methods$name == "Camera"][[1]]

paste0(method_params$name, " (", method_params$type, ", ", method_params$default, ")")


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

loginfo("Create a new request object using 'Camera' for the gene set analysis")
start_request <-ReactomeAnalysisRequest(method = args$enrichment_method)

start_request

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Setting parameters")
#To get a list of supported parameters for each method, use the ``get_reactome_methods`` function (see above).

#Parameters are simply set using the ``set_parameters`` function:

# set the maximum number of allowed missing values to 50%
request_param_defined <- set_parameters( request = start_request,
                              max_missing_values = args$max_missing_values, # 0.5 is the recommended default
                              create_reactome_visualization = TRUE,
                              create_reports = TRUE,
                              use_interactors = args$use_interactors,
                              include_disease_pathways = args$include_disease_pathways,
                              email=args$email)

request_param_defined

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Create object for storing proteomics data
input_proteomics <- new("EList")

input_proteomics$E <- norm_abundance_mat

proteomics_samples <- design_mat_cln %>% dplyr::select(-one_of(c(sample_id)) )

input_proteomics$samples <- proteomics_samples

input_proteomics$genes <- rownames( norm_abundance_mat )

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

cmriReactomeAddDataset <- function(contrast, contrast_name, input_request, elist_object, formula_string ) {

  groupA <- str_replace( names( contrast[contrast == 1] )[1], paste0("^", quo_name(enquo(group_id))), "" )
  groupB <- str_replace( names( contrast[contrast == -1] )[1], paste0("^", quo_name(enquo(group_id))), "" )

  additional_factors_for_adjustment <- formula_string %>%
    str_split( " ") %>%
    .[[1]] %>%
    discard( ~str_detect(., "\\~|\\+|\\-|^$|^\\d{1}$")  ) %>%
    discard( ~str_detect(., quo_name(enquo(group_id))))

  if(length(additional_factors_for_adjustment) == 0) {
    additional_factors_for_adjustment <- NULL
  }

  my_request <- add_dataset(request = input_request,
                            expression_values = elist_object,
                            name = contrast_name,
                            type = "microarray_norm",
                            comparison_factor = quo_name(enquo(group_id)) ,
                            comparison_group_1 = groupB,
                            comparison_group_2 = groupA,
                            additional_factors = additional_factors_for_adjustment )

  return(my_request)

}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

loginfo("The dataset can now simply be added to the request using the `add_dataset` function:")

partial_cmriReactomeAddDataset <- partial( cmriReactomeAddDataset,
                                           elist_object = input_proteomics,
                                           formula_string )

recursiveAddDataset <- function( input_contrast, input_contrast_name, input_request_status ) {

  if( length(input_contrast) !=  length( input_contrast_name) ) {
    logwarn("Function recursiveAddDataset, input_contrast should have same length as input_contrast_name")
    stop()
  }

  updated_request_status <- partial_cmriReactomeAddDataset(contrast=input_contrast[[1]],  contrast_name=input_contrast_name[1],  input_request=input_request_status )

  if( length( input_contrast) == 1) {
    return( updated_request_status)
  } else  {

    return( recursiveAddDataset(input_contrast = input_contrast[-1],
                                input_contrast_name =input_contrast_name[-1],
                                input_request_status = updated_request_status) )

  }
}

requests_list <- recursiveAddDataset( input_contrast=lists_of_contrasts,
                                     input_contrast_name=names( lists_of_contrasts),
                                     input_request_status=request_param_defined)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Submitting the request to Reactome server")

result <- perform_reactome_analysis( request = requests_list,
                                     verbose=TRUE,
                                     compress = FALSE)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

loginfo("Investigating the result")
names(result)
result_types(result)

purrr::walk( names(result), function(contrast_name) {

  proteomics_fc <- get_result(result, type = "fold_changes", name = contrast_name)
  captured_output <- capture.output(
    vroom::vroom_write( proteomics_fc, file.path(output_dir, paste0( "reactome_fold_hanges.", contrast_name ,".tab")) )
    ,type = "message" )
  logdebug(captured_output)

  captured_output <- capture.output(
    writexl::write_xlsx( proteomics_fc, file.path( output_dir, paste0( "reactome_fold_changes.", contrast_name , ".xlsx")) )
    ,type = "message" )
  logdebug(captured_output)
} )

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


loginfo("Directly merge the pathway level data for all result sets using the pathways function.")

combined_pathways <- pathways(result)
captured_output<-capture.output(
  vroom::vroom_write( combined_pathways, file.path(output_dir, paste0( "combined_pathways",".tab")) )
  ,type = "message"
)
logdebug(captured_output)

captured_output<-capture.output(
  writexl::write_xlsx( combined_pathways, file.path( output_dir, paste0( "combined_pathways",".xlsx")) )
  ,type = "message"
)
logdebug(captured_output)


captured_output<-capture.output(
  { sink( file.path( output_dir, "reactome_links.txt") )
    reactome_links(result)
    sink() }
  ,type = "message"
)
logdebug(captured_output)


if(length( names(result)) > 1) {

  captured_output <- capture.output(
    {pdf( file.path( output_dir, "plot_correlations.pdf") )
    plot_correlations(result )
    dev.off()}
    , type = "message"
  )
  logdebug(captured_output)
}

  captured_output <- capture.output(
    saveRDS(result, file = file.path( output_dir, "my_ReactomeGSA_result.rds") )
    ,type = "message"
  )
  logdebug(captured_output)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
te<-toc(quiet = TRUE)
loginfo("%f sec elapsed",te$toc-te$tic)
writeLines(capture.output(sessionInfo()), file.path(args$output_dir,"sessionInfo.txt"))





