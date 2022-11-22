#!/usr/bin/env Rscript

# Author(s): Ignatius Pang, Pablo Galaviz
# Email: cmri-bioinformatics@cmri.org.au
# Children’s Medical Research Institute, finding cures for childhood genetic diseases


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Test if BioManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install()
}

# load pacman package manager
if (!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

p_load(optparse)
p_load(tictoc)
p_load(tidyverse)

p_load(plotly)
p_load(vroom)
p_load(writexl)
p_load(ggplot2)
p_load(ggpubr)
p_load(magrittr)
p_load(knitr)
p_load(rlang)
p_load(ggrepel)

p_load(PhosR)
p_load(statmod)
p_load(limma)
p_load(qvalue)
p_load(ruv)
p_load(mixOmics)

p_load(ProteomeRiver)
p_load(configr)
p_load(logging)
p_load(svglite)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

tic()
set.seed(123456)

command_line_options <- commandArgs(trailingOnly = TRUE)
#Note: options with default values are ignored in the configuration file parsing.

parser <- OptionParser(add_help_option = TRUE)

#Note: options with default values are ignored in the configuration file parsing.
parser <- add_option(parser, c("-d", "--debug"), action = "store_true", default = FALSE,
                     help = "Print debugging output  [default %default]")

parser <- add_option(parser, c("-s", "--silent"), action = "store_true", default = FALSE,
                     help = "Only print critical information to the console.  [default %default]")

parser <- add_option(parser, c("-n", "--no_backup"), action = "store_true", default = FALSE,
                     help = "Deactivate backup of previous run.  [default %default]")

parser <- add_option(parser, c("-c", "--config"), type = "character", default = "config_prot.ini",
                     help = "Configuration file.  [default %default]",
                     metavar = "string")

parser <- add_option(parser, c("-o", "--output_dir"), type = "character",
                     help = "Directory path for all results files.",
                     metavar = "string")

parser <- add_option(parser, c("-l", "--log_file"), type = "character", default = "output.log",
                     help = "Name of the logging file.  [default %default]",
                     metavar = "string")

parser <- add_option(parser, c( "--normalization"), type = "character",
                     help = "character string specifying the normalization method to be used. Choices are none, scale, quantile or cyclicloess",
                     metavar = "string")

parser <- add_option(parser, c( "--imputation"), type = "logical",
                     help = "logical, whether imputation should be used to assign missing vlues.",
                     metavar = "double")


#Options without a default value have the following priority: configuration file < command line argument
parser <- add_option(parser, "--max_num_samples_miss_per_group", type = "integer",
                     help = "Remove protein if it exceeds this maximum number of samples with missing values per experimental group",
                     metavar = "integer")

parser <- add_option(parser, "--abundance_threshold", type = "integer",
                     help = "Abundance threshold above which the protein in the sample is accepted for analysis",
                     metavar = "integer")

parser <- add_option(parser, "--counts_table_file", type = "character",
                     help = "Input file with the protein abundance values",
                     metavar = "string")

parser <- add_option(parser, "--formula_string", type = "character",
                     help = "A string representing the formula for input into the model.frame function. (e.g. ~ 0 + group).",
                     metavar = "string")

parser <- add_option(parser, "--design_matrix_file", type = "character",
                     help = "Input file with the design matrix",
                     metavar = "string")

parser <- add_option(parser, "--sample_id", type = "character",
                     help = "A string describing the sample ID. This must be a column that exists in the design matrix.",
                     metavar = "string")

parser <- add_option(parser, "--bio_rep_col", type = "character",
                     help = "A string describing the biological replicate group ID. This must be a column that exists in the design matrix.",
                     metavar = "string")

parser <- add_option(parser, "--tech_rep_col", type = "character",
                     help = "A string describing the technical replicate group ID. This must be a column that exists in the design matrix.",
                     metavar = "string")

parser <- add_option(parser, "--group_id", type = "character",
                     help = "A string describing the experimental group ID. This must be a column that exists in the design matrix.",
                     metavar = "string")

parser <- add_option(parser, "--average_replicates_id", type = "character",
                     help = "A string describing the technical replicates that needs to be averaged together. This must be a column that exists in the design matrix.",
                     metavar = "string")

parser <- add_option(parser, "--row_id", type = "character",
                     help = "A string describing the row id.",
                     metavar = "string")

parser <- add_option(parser, "--file_prefix", type = "character",
                     help = "A string to indicate the type of analysis and is used in the file name of the output results table.",
                     metavar = "string")

parser <- add_option(parser, "--plots_format", type = "character",
                     help = "A comma separated strings to indicate the fortmat output for the plots [pdf,png,svg].",
                     metavar = "string")

#parse comand line arguments first.
args <- parse_args(parser)

#parse and merge the configuration file options.
if (args$config != "") {
  args <- config.list.merge(eval.config(file = args$config, config="impute_missing_values"), args)
}

args <- setArgsDefault(args, "output_dir", as_func=as.character, default_val="impute_missing_values" )

createOutputDir(args$output_dir, args$no_backup)

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


testRequiredArguments(args, c(
   "formula_string"
  , "sample_id"
  , "group_id"
  , "row_id"
  , "file_prefix"
  , "counts_table_file"
  , "contrasts_file"
  , "design_matrix_file"
))

testRequiredFiles(c(
  args$counts_table_file
  , args$contrasts_file
  , args$design_matrix_file
))

args <- setArgsDefault(args, "max_num_samples_miss_per_group", as_func=as.integer, default_val=NA )
args <- setArgsDefault(args, "imputation", as_func=as.logical, default_val=FALSE )

args<-parseType(args,
                c("max_num_samples_miss_per_group"
                  ,"abundance_threshold"
                )  ,as.integer)

args<-parseString(args,
                  c("formula_string"
                    ,"sample_id"
                    , "tech_rep_col"
                    , "bio_rep_col"
                    , "average_replicates_id"
                    ,"group_id"
                    ,"row_id"
                    ,"file_prefix"
                    ,"plots_format"
                  ))

if(isArgumentDefined(args,"plots_format"))
{
  args <- parseList(args,
                    c("plots_format"))
}else {
  logwarn("plots_format is undefined, default output set to pdf.")
  args$plots_format <- list("pdf")
}



if (  ! isArgumentDefined(args,"average_replicates_id") ) {
  logwarn("average_replicates_id is NA")
  args$average_replicates_id <- NA
}


## Check the normalization option
if ( ! isArgumentDefined(args, "normalization" ) ) {
  logwarn("normalization is NA, defaults to the 'scale' method")
  args$normalization <- "scale"
}

if (  ! args$normalization %in% c("none", "scale", "quantile", "cyclicloess") ) {
  logerror( paste0("the normalization method specified ", args$normalization, " is not valid, please check carefully."))
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

loginfo("Read file with counts table %s", args$counts_table_file)
captured_output<-capture.output(
  evidence_tbl_filt <- vroom::vroom(args$counts_table_file, delim = "\t") %>%
    dplyr::select(one_of(c(args$row_id)), matches(args$group_pattern))
  ,type = "message"
)
logdebug(captured_output)

if (!args$row_id %in% colnames(evidence_tbl_filt)) {
  logerror("The row ID specified in --row-id was not found in the input counts file.")
  q()
}
if (length(which(str_detect(setdiff(colnames(evidence_tbl_filt), args$row_id), args$group_pattern))) == 0) {
  logerror("The input counts file did not contain any columns that matches the string provided in --group-pattern=%s", args$group_pattern)
  q()
}


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

if ( !is.na(args$replicate_group_id) &
     args$replicate_group_id != args$group_id )  {
  if( length(which(c(args$replicate_group_id) %in% colnames(design_mat_cln))) != 1)  {
    logerror("replicate_group_id is not matching to the column names used in the design matrix.")
    q()
  }
}

cols_for_analysis <- design_mat_cln %>% pull(as.name(args$sample_id))

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



design_mat_cln <- design_mat_cln %>% 
  mutate(full_name = as.name(args$sample_id)) %>%
  mutate( biological_replicates = replicates) 

length(design_mat_cln$full_name) == length( colnames( prot_tbl)[2:ncol(prot_tbl) ])


colnames( prot_tbl)[2:ncol(prot_tbl) ] <- design_mat_cln$full_name


prot_cln <- prot_tbl  %>%
  pivot_longer( cols= design_mat_cln$full_name,
                values_to="Intensity",
                names_to="Sample_ID")   %>%
  arrange( uniprot_acc, Sample_ID) %>%
  dplyr::filter( Intensity > 0 ) %>%
  dplyr::mutate( LogIntesity = log2(Intensity))


prot_sort <- prot_cln %>%
  pivot_wider( id_cols = c("uniprot_acc"),
               values_from = "LogIntesity",
               names_from="Sample_ID")  %>%
  mutate( row_id = row_number()) %>%
  relocate( row_id, .before="uniprot_acc")


prot_sort_mat <- prot_sort %>%
  dplyr::select( -uniprot_acc ) %>%
  column_to_rownames("row_id") %>% 
  as.matrix()

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

norm_prot_sort_mat <- normalizeBetweenArrays(prot_sort_mat, method = args$normalization)

plotRle( t( prot_sort_mat ) )
plotRle( t( norm_prot_sort_mat ) )

plot(density( prot_sort_mat[,1][prot_sort_mat[,1] > 0 ], na.rm =TRUE) )

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Calculate average value from technical replicates if this is applicable
counts_rnorm.log.for.contrast <- counts_rnorm.log.quant
design_mat_updated <- design_mat_cln
if (!is.na( args$average_replicates_id)) {
  
  counts_rnorm.log.for.contrast <-  averageValuesFromReplicates(counts_rnorm.log.quant,
                                                                design_mat_cln,
                                                                args$row_id,
                                                                args$sample_id,
                                                                args$average_replicates_id)
  
  design_mat_updated <- design_mat_cln %>%
    mutate( !!rlang::sym(args$sample_id) :=  !!rlang::sym(args$average_replicates_id) ) %>%
    dplyr::select( one_of( args$sample_id, args$group_id)) %>%
    distinct()
  
  rownames( design_mat_updated) <- design_mat_updated[,args$sample_id]
}



# Averaging of technical replicates 

prot_avg <- norm_prot_sort_mat %>%
  as.data.frame %>%
  rownames_to_column("row_id") %>%
  mutate( row_id = as.integer(row_id )) %>%
  left_join( prot_sort %>%
               dplyr::select( row_id, uniprot_acc), by ="row_id") %>%
  relocate( uniprot_acc, .after="row_id") %>%
  pivot_longer( cols = matches("\\d+"),
                values_to = "NormLogIntensity",
                names_to="Sample_ID") %>%
  left_join( design_mat_cln, by = c("Sample_ID" = "full_name")) %>%
  dplyr::filter( !is.na( NormLogIntensity) ) %>%
  group_by(row_id, uniprot_acc, biological_replicates ) %>%
  summarise( AvgNormLogIntensity = mean( NormLogIntensity) ) %>%
  ungroup()


prot_avg_sort <- prot_avg %>%
  arrange( biological_replicates, uniprot_acc) %>%
  pivot_wider( id_cols = c("row_id", "uniprot_acc" ),
               values_from = "AvgNormLogIntensity",
               names_from="biological_replicates")


prot_avg_sort_mat <- prot_avg_sort%>%
  dplyr::select( -uniprot_acc  ) %>%
  column_to_rownames("row_id") %>% 
  as.matrix()

plotRle( t( prot_avg_sort_mat ) )

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------























