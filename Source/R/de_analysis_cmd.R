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

parser <- add_option(parser, c("-c", "--config"), type = "character", default = "config.ini",
                     help = "Configuration file.  [default %default]",
                     metavar = "string")

parser <- add_option(parser, c("-o", "--output_dir"), type = "character",
                     help = "Directory path for all results files.",
                     metavar = "string")

parser <- add_option(parser, c("-l", "--log_file"), type = "character", default = "output.log",
                     help = "Name of the logging file.  [default %default]",
                     metavar = "string")

parser <- add_option(parser, c( "--treat_lfc_cutoff"), type = "double",
                     help = "The minimum log2-fold-change below which changes not considered scientifically meaningful. Used in treat function of the limma library.",
                     metavar = "double")

parser <- add_option(parser, c( "--normalization"), type = "character",
                     help = "character string specifying the normalization method to be used. Choices are none, scale, quantile or cyclicloess",
                     metavar = "string")

parser <- add_option(parser, c( "--eBayes_trend"), type = "double",
                     help = "logical, should an intensity-trend be allowed for the prior variance? Default is that the prior variance is constant.",
                     metavar = "double")

parser <- add_option(parser, c( "--eBayes_robust"), type = "double",
                     help = "logical, should the estimation of df.prior and var.prior be robustified against outlier sample variances?.",
                     metavar = "double")

#Options without a default value have the following priority: configuration file < command line argument
parser <- add_option(parser, "--max_num_samples_miss_per_group", type = "integer",
                     help = "Remove protein if it exceeds this maximum number of samples with missing values per experimental group",
                     metavar = "integer")

parser <- add_option(parser, "--abundance_threshold", type = "integer",
                     help = "Abundance threshold above which the protein in the sample is accepted for analysis",
                     metavar = "integer")

parser <- add_option(parser, "--group_pattern", type = "character",
                     help = "Regular expression pattern to identify columns with abundance values belonging to the experiment.",
                     metavar = "string")

parser <- add_option(parser, "--q_val_thresh", type = "double",
                     help = "q-value threshold below which a protein has statistically significant differetial expression",
                     metavar = "double")

parser <- add_option(parser, "--control_genes_q_val_thresh", type = "double",
                     help = "q-value threshold below which a protein has statistically significant differetial expression, used for control genes",
                     metavar = "double")

parser <- add_option(parser, "--ruv_k", type = "integer",
                     help = "The number of unwanted factors to use.",
                     metavar = "integer")

parser <- add_option(parser, "--num_neg_ctrl", type = "integer",
                     help = "The number of negative control proteins to use.",
                     metavar = "number")

parser <- add_option(parser, "--ruv_method", type = "character",
                     help = "A string representing the ruv3 method to use.",
                     metavar = "character")

parser <- add_option(parser, "--counts_table_file", type = "character",
                     help = "Input file with the protein abundance values",
                     metavar = "string")

parser <- add_option(parser, "--test_pairs_file", type = "character",
                     help = "Input file with a table listing all the pairs of experimental groups to compare. First column represents group A and second column represents group B. Linear model comparisons (e.g. Contrasts) would be group B minus group A.",
                     metavar = "string")

parser <- add_option(parser, "--contrasts_file", type = "character",
                     help = "Input file with a table listing all comparisons to be made in string, one comparison per line (e.g. groupB.vs.group_A = groupB - groupA).",
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

parser <- add_option(parser, "--replicate_group_id", type = "character",
                     help = "A string describing the replicate group ID. This must be a column that exists in the design matrix.",
                     metavar = "string")

parser <- add_option(parser, "--group_id", type = "character",
                     help = "A string describing the experimental group ID. This must be a column that exists in the design matrix.",
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
  args <- config.list.merge(eval.config(file = args$config, config="de_analysis"), args)
}

args <- setArgsDefault(args, "output_dir", as_func=as.character, default_val="de_analysis" )

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
  "group_pattern"
  , "formula_string"
  , "ruv_k"
  , "q_val_thresh"
  , "control_genes_q_val_thresh"
  , "num_neg_ctrl"
  , "ruv_method"
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

args <- setArgsDefault(args, "treat_lfc_cutoff", as_func=as.double, default_val=NA )
args <- setArgsDefault(args, "eBayes_trend", as_func=as.logical, default_val=FALSE )
args <- setArgsDefault(args, "eBayes_robust", as_func=as.logical, default_val=FALSE )
args <- setArgsDefault(args, "q_val_thresh", as_func=as.double, default_val=NaN )
args <- setArgsDefault(args, "control_genes_q_val_thresh", as_func=as.double, default_val=NaN )

args<-parseType(args,
  c("ruv_k"
    ,"num_neg_ctrl"
    ,"max_num_samples_miss_per_group"
    ,"abundance_threshold"
  )  ,as.integer)

args<-parseString(args,
  c("group_pattern"
    ,"test_pairs_file"
    ,"formula_string"
    ,"sample_id"
    , "replicate_group_id"
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

if (args$group_pattern == "") {
  logwarn("Empty group pattern string, using \\d+")
  args$group_pattern <- "\\d+"
}

## Clean up replicate group ID and then take default value
if (  ! isArgumentDefined(args,"replicate_group_id") ) {
  logwarn("Replicate_group_id is NA")
  args$replicate_group_id <- NA
}

if ( is.na(args$replicate_group_id ) &
     ( !is.na(args$group_id ) |
       args$group_id == "")) {
  args$replicate_group_id <- args$group_id
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

# test_pairs <- NA
# if( limma_method == "pairs" & args$test_pairs_file != "") {
#   print("Read file with pairs of experimental groups to test.")
#   test_pairs <- vroom::vroom( args$test_pairs_file)
# }


loginfo("Read file with lists of experimental contrasts to test %s", args$contrasts_file)
captured_output<-capture.output(
contrasts_tbl <- vroom::vroom(args$contrasts_file, delim = "\t")
    ,type = "message"
)
logdebug(captured_output)

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
loginfo("Get table with columns for analysis.")

evidence_tbl_col <- evidence_tbl_filt[, c(args$row_id, cols_for_analysis) ]


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Remove empty proteins without abundance data.")
cln_dat_wide_unsorted <- ProteomeRiver::removeEmptyRows(evidence_tbl_filt,
                                                        col_pattern = args$group_pattern,
                                                        row_id = !!rlang::sym(args$row_id))

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Count the total number of missing values in total

table_value <- table(is.infinite(data.matrix(log2(cln_dat_wide_unsorted[, c(colnames(cln_dat_wide_unsorted)[1], cols_for_analysis)]%>%
                                                    column_to_rownames(args$row_id)))))

loginfo("Count the number of missing values for each sample before removing proteins with some missing values: %d",
        table_value["TRUE"])

plot_num_missing_values_before <- plotNumMissingValues(cln_dat_wide_unsorted[, cols_for_analysis])

for( format_ext in args$plots_format) {
  file_name<-file.path(args$output_dir,paste0("num_missing_values_before_filtering.",format_ext))
  captured_output<-capture.output(
    ggsave(filename = file_name, plot = plot_num_missing_values_before, limitsize = FALSE)
    ,type = "message"
  )
  logdebug(captured_output)
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cln_dat_wide_cleaned <- NA
if ("max_num_samples_miss_per_group" %in% names(args)) {
  loginfo("Remove row with any missing values")

  cln_dat_wide_cleaned <- removeRowsWithMissingValues(cln_dat_wide_unsorted,
                                                      matches(args$group_pattern),
                                                      design_mat_cln %>%
                                                            dplyr::mutate(Sample_ID = as.character(args$sample_id)),
                                                      !!rlang::sym(args$sample_id),
                                                      !!rlang::sym(args$row_id),
                                                      !!rlang::sym(args$group_id),
                                                      args$max_num_samples_miss_per_group,
                                                      args$abundance_threshold)

} else {
  cln_dat_wide_cleaned <- cln_dat_wide_unsorted
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if ( length(cols_for_analysis) != length(intersect(colnames(cln_dat_wide_cleaned)[-1], cols_for_analysis))) {
  logerror("Column names in design matrix table not consistent with column names in abundance table.")
  q()
}

# Sort the columns in the abundance table according
if (length(which(colnames(cln_dat_wide_cleaned)[-1] %in% cols_for_analysis)) != length(cols_for_analysis)) {
  logwarn("Sorting the order of the column of the abundance matrix according to the order in the design matrix.")
}

cln_dat_wide <- cln_dat_wide_cleaned[, c(colnames(cln_dat_wide_cleaned)[1], cols_for_analysis)]

loginfo("Assign the indexing row names to the data frame.")
counts_filt <- cln_dat_wide %>%
  column_to_rownames(args$row_id)

vroom::vroom_write( cln_dat_wide,
                    file.path( args$output_dir,
                               "raw_counts_after_removing_empty_rows.tsv"))

writexl::write_xlsx( cln_dat_wide,
                     file.path( args$output_dir,
                                "raw_counts_after_removing_empty_rows.xlsx"))

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ruvIII_replicates_matrix <- getRuvIIIReplicateMatrix(design_mat_cln,
                                                     !!rlang::sym(args$sample_id),
                                                     !!rlang::sym(args$replicate_group_id))

logdebug(ruvIII_replicates_matrix)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# sample_rows_lists <- get_rows_to_keep_list( cln_dat_wide,
#                                             matches(args$group_pattern),
#                                             design_mat_cln %>%
#                                               dplyr::mutate( Sample_ID = as.character(Sample_ID)),
#                                             !!rlang::sym(args$sample_id),
#                                             !!rlang::sym(args$row_id),
#                                             !!rlang::sym(args$group_id),
#                                             args$min_num_samples_per_group,
#                                             args$abundance_threshold)
#
# saveRDS(sample_rows_lists, file.path(args$output_dir, "keep_sample_rows_lists.RDS" ))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Count the total number of missing values in total
loginfo("Count the number of missing values: %d",table(is.infinite(data.matrix(log2(counts_filt))))["TRUE"])

plot_num_missing_values <- plotNumMissingValues(counts_filt[, cols_for_analysis])

for( format_ext in args$plots_format) {
  file_name<-file.path(args$output_dir,paste0("num_missing_values.",format_ext))
  captured_output<-capture.output(
ggsave(filename = file_name, plot = plot_num_missing_values,limitsize = FALSE)
    ,type = "message"
)
  logdebug(captured_output)
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

loginfo("Count the number of values available for each sample: %d",table(!is.infinite(data.matrix(log2(counts_filt))))["TRUE"])

plot_num_of_values <- plotNumOfValues(counts_filt[, cols_for_analysis])

for( format_ext in args$plots_format) {
  file_name <- file.path(args$output_dir,paste0("num_of_values.",format_ext))
  captured_output <- capture.output(
    ggsave(filename = file_name, plot = plot_num_of_values, limitsize = FALSE)
    ,type = "message" )
  logdebug(captured_output)
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Mark entries with zero values as NA, take log of values.")
na_values_marker <- (counts_filt == 0)

counts_na <- counts_filt

counts_na[na_values_marker] <- NA

counts_na.log <- log2(counts_na)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Median centering.")
counts_na.log.quant <- normalizeBetweenArrays(counts_na.log, method = args$normalization)


vroom::vroom_write(as.data.frame(counts_na.log.quant) %>%
                     rownames_to_column(args$row_id),
                   file.path(args$output_dir, "counts_after_median_scaling_before_imputation.tsv"))

writexl::write_xlsx( as.data.frame(counts_na.log.quant) %>%
                       rownames_to_column(args$row_id),
                     file.path( args$output_dir,
                                "counts_after_median_scaling_before_imputation.xlsx"))

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

loginfo("Missing value imputation")

counts_rnorm.log.quant <- counts_na.log.quant %>%
  as.data.frame %>%
  mutate_all(~{ imputePerCol(., width = 0.3, downshift = 1.8) }) %>%
  as.matrix

counts_rnorm.is_nan <- counts_rnorm.log.quant %>%
  as.data.frame %>%
  mutate_all(~{ !is.finite(.) }) %>%
  as.matrix

vroom::vroom_write( as.data.frame(counts_rnorm.log.quant) %>%
                      rownames_to_column(args$row_id),
                    file.path(args$output_dir, "counts_after_median_scaling_and_imputation.tsv"))

writexl::write_xlsx( as.data.frame(counts_rnorm.log.quant) %>%
                       rownames_to_column(args$row_id),
                     file.path( args$output_dir,
                                "counts_after_median_scaling_and_imputation.xlsx"))

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Run statistical tests without RUV.")

ID <- rownames(counts_rnorm.log.quant)

type_of_grouping <- getTypeOfGrouping(design_matrix = design_mat_cln, group_id = args$group_id, sample_id = args$sample_id)

list_rnorm.log.quant.ruv.r0 <- NA
myRes_rnorm.log.quant <- NA

list_rnorm.log.quant.ruv.r0 <- runTestsContrasts(counts_rnorm.log.quant,
                                                 contrast_strings = contrasts_tbl[, 1][[1]],
                                                 design_matrix = design_mat_cln,
                                                 formula_string = args$formula_string,
                                                 weights = NA,
                                                 treat_lfc_cutoff = as.double(args$treat_lfc_cutoff),
                                                 eBayes_trend = as.logical(args$eBayes_trend),
                                                 eBayes_robust = as.logical(args$eBayes_robust))

myRes_rnorm.log.quant <- list_rnorm.log.quant.ruv.r0$results

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# print( "Save the list of differentially expressed proteins, without using RUV.")
# save_de_protein_list(myRes_rnorm.log.quant,
#                      row_id="protein_id",
#                      sort_by_column =q.mod,
#                      results_dir = args$output_dir,
#                      file_suffix = "_without_ruv.tsv")


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Find the list of negative control genes using ANOVA.")
control_genes_index <- getNegCtrlProtAnova(counts_rnorm.log.quant,
                                           design_matrix = design_mat_cln,
                                           group_column = args$group_id,
                                           num_neg_ctrl = args$num_neg_ctrl,
                                           q_val_thresh = args$control_genes_q_val_thresh)

vroom::vroom_write(data.frame(temp_col = names(control_genes_index),
                              is_control_genes = control_genes_index) %>%
                     set_colnames(c(args$row_id, "is_control_genes")),
                   file.path(args$output_dir, "ctrl_genes_list_ruv3.tsv"),
                   delim = "\t")

writexl::write_xlsx( data.frame( temp_col = names(control_genes_index),
                                 is_control_genes = control_genes_index) %>%
                       set_colnames(c(args$row_id, "is_control_genes")),
                     file.path(args$output_dir, "ctrl_genes_list_ruv3.xlsx"))

loginfo("Total Num. Genes %d", length(control_genes_index))

loginfo("Num. Control Genes %d", length(which(control_genes_index)))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Draw canonical correlation plot.")
ruv_groups <- data.frame(temp_column = colnames(counts_rnorm.log.quant)) %>%
  dplyr::rename(!!rlang::sym(args$sample_id) := "temp_column") %>%
  left_join(design_mat_cln %>%
              dplyr::mutate(!!rlang::sym(args$sample_id) := as.character(!!rlang::sym(args$sample_id))),
            by = args$sample_id)

cancorplot_r1 <- ruv_cancorplot(t(counts_rnorm.log.quant),
                                X = ruv_groups %>% pull(!!rlang::sym(args$group_id)),
                                ctl = control_genes_index)


for( format_ext in args$plots_format) {
  file_name<-file.path(args$output_dir,paste0("cancor_plot_round_1.",format_ext))
  captured_output<-capture.output(
    ggsave(plot = cancorplot_r1, filename =  file_name, limitsize = FALSE)
    , type = "message"
  )
  logdebug(captured_output)
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Run RUVIII on the input counts table.")
counts_rnorm.log.ruvIII_v1 <- cmriRUVfit(counts_rnorm.log.quant, X = ruv_groups$group, control_genes_index, Z = 1, k = args$ruv_k,
                                         method = args$ruv_method, M = ruvIII_replicates_matrix) %>% t()

vroom::vroom_write(counts_rnorm.log.ruvIII_v1 %>%
                     as.data.frame %>%
                     rownames_to_column(args$row_id),
                   file.path(args$output_dir, "normalized_counts_after_ruv.tsv"))

writexl::write_xlsx(counts_rnorm.log.ruvIII_v1 %>%
                      as.data.frame %>%
                      rownames_to_column(args$row_id),
                    file.path(args$output_dir, "normalized_counts_after_ruv.xlsx"))

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Draw the RLE and PCA plots.")
rle_pca_plots_arranged <- rlePcaPlotList(list_of_data_matrix = list(counts_rnorm.log.quant, counts_rnorm.log.ruvIII_v1),
                                         design_matrix = design_mat_cln,
                                         sample_id_column = !!rlang::sym(args$sample_id),
                                         group_column = !!rlang::sym(args$group_id),
                                         list_of_descriptions = list("Before RUVIII", "After RUVIII"))


for( format_ext in args$plots_format) {
  file_name<-file.path(args$output_dir,paste0("rle_pca_plots.",format_ext))
  captured_output<-capture.output(
    ggsave(plot = rle_pca_plots_arranged, filename =  file_name, limitsize = FALSE, width = 14, height = 14)
    , type = "message"
  )
  logdebug(captured_output)
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### Compare more principal components
# pca.res <- pca(t(as.matrix(counts_rnorm.log.ruvIII_v1)), ncomp=5)
#
# output <- pca.res$variates$X %>%
#   as.data.frame %>%
#   rownames_to_column("Sample_ID") %>%
#   left_join(design_mat_cln, by ="Sample_ID")
#
# output
# p_load(GGally)
# ggpairs(output,          # Data frame
#         columns = 2:6,        # Columns
#         aes(color = group )  # Color by group (cat. variable)
#            )

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Compare the different experimental groups and obtain lists of differentially expressed proteins.")

list_rnorm.log.quant.ruv.r1 <- NA
myRes_rnorm.log.quant.ruv.r1 <- NA

list_rnorm.log.quant.ruv.r1 <- runTestsContrasts(counts_rnorm.log.ruvIII_v1,
                                                 contrast_strings = contrasts_tbl[, 1][[1]],
                                                 design_matrix = design_mat_cln,
                                                 formula_string = args$formula_string,
                                                 weights = NA)

myRes_rnorm.log.quant.ruv.r1 <- list_rnorm.log.quant.ruv.r1$results

saveRDS( list_rnorm.log.quant.ruv.r1$fit.eb,
         file.path(args$output_dir, "fit.eb.RDS" ) )

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# print("Save the lists of diffierentially expressed proteins (e.g. log fold-change, q-values).")
# save_de_protein_list(myRes_rnorm.log.quant.ruv.r1,
#                      row_id="protein_id",
#                      sort_by_column =q.mod,
#                      results_dir = args$output_dir,
#                      file_suffix = "_round_1_ruv.tsv")

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Prepare data for drawing the volcano plots")
selected_data <- getSignificantData(list_of_de_tables = list(myRes_rnorm.log.quant,
                                                             myRes_rnorm.log.quant.ruv.r1),
                                    list_of_descriptions = list("No RUV",
                                                                  "RUV applied"),
                                    row_id = !!sym(args$row_id),
                                    p_value_column = p.mod,
                                    q_value_column = q.mod,
                                    fdr_value_column = fdr.mod,
                                    log_q_value_column = lqm,
                                    log_fc_column = logFC,
                                    comparison_column = comparison,
                                    expression_column = expression,
                                    facet_column = analysis_type,
                                    q_val_thresh = args$q_val_thresh) %>%
  dplyr::rename(log2FC = "logFC")

## Write all the results in one single table
selected_data %>%
  dplyr:::select(-colour, -lqm) %>%
  vroom::vroom_write(file.path(args$output_dir, "lfc_qval_long.tsv"))

selected_data %>%
  dplyr:::select(-colour, -lqm) %>%
  writexl::write_xlsx(file.path(args$output_dir, "lfc_qval_long.xlsx"))

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

loginfo("Print the volcano plots")

volplot_gg.all <- plotVolcano(selected_data,
                              log_q_value_column = lqm,
                              log_fc_column = log2FC,
                              q_val_thresh = args$q_val_thresh,
                              formula_string = "analysis_type ~ comparison")


for( format_ext in args$plots_format) {
  file_name <- file.path(args$output_dir,paste0("volplot_gg_all.",format_ext))
  captured_output <- capture.output(
    ggsave(filename = file_name, plot = volplot_gg.all, width = 7.29, height = 6)
    , type = "message" )
  logdebug(captured_output)
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Count the number of up or down significnat differentially expressed proteins.")

num_sig_de_molecules <- printCountDeGenesTable(list_of_de_tables = list(myRes_rnorm.log.quant,
                                                                    myRes_rnorm.log.quant.ruv.r1),
                                           list_of_descriptions = list("No RUV",
                                                                           "RUV applied"),
                                           formula_string = "analysis_type ~ comparison")

for( format_ext in args$plots_format) {
  file_name<-file.path(args$output_dir,paste0("num_de_genes_barplot.",format_ext))
  captured_output<-capture.output(
ggsave(filename = file_name,
       plot = num_sig_de_molecules$plot,
       height = 10,
       width = 7)
    , type = "message"
  )
  logdebug(captured_output)
}

vroom::vroom_write(num_sig_de_molecules$table,
                   file.path(args$output_dir,
                             "num_significant_differentially_abundant_all.tab"))

writexl::write_xlsx(num_sig_de_molecules$table,
                    file.path(args$output_dir,
                              "num_significant_differentially_abundant_all.xlsx"))

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Print p-values distribution figure")
pvalhist <- printPValuesDistribution(selected_data,
                                     p_value_column = p.mod,
                                     formula_string = "analysis_type ~ comparison")

#logdebug(head(pvalhist))

for( format_ext in args$plots_format) {
  file_name<-file.path(args$output_dir,paste0("p_values_distn.",format_ext))
  captured_output<capture.output(
ggsave(filename = file_name,
       plot = pvalhist,
       height = 10,
       width = 7)
    , type = "message"
  )
  logdebug(captured_output)
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

loginfo("Create wide format output file")

norm_counts <- counts_rnorm.log.ruvIII_v1 %>%
  as.data.frame %>%
  set_colnames(paste0(colnames(counts_rnorm.log.ruvIII_v1), ".log2norm")) %>%
  rownames_to_column(args$row_id)


raw_counts <- counts_filt %>%
  as.data.frame %>%
  set_colnames(paste0(colnames(counts_filt), ".raw")) %>%
  rownames_to_column(args$row_id)

de_proteins_wide <- selected_data %>%
  dplyr::filter(analysis_type == "RUV applied") %>%
  dplyr::select(-lqm, -colour, -analysis_type) %>%
  pivot_wider(id_cols = c(!!sym(args$row_id)),
              names_from = c(comparison),
              names_sep = ":",
              values_from = c(log2FC, q.mod, p.mod)) %>%
  left_join(norm_counts, by = args$row_id) %>%
  left_join(raw_counts, by = args$row_id) %>%
  dplyr::arrange(across(matches("q.mod"))) %>%
  distinct()

#logdebug(head(de_proteins_wide))

vroom::vroom_write( de_proteins_wide,
                    file.path( args$output_dir,
                               paste0(args$file_prefix, "_wide.tsv")))

writexl::write_xlsx( de_proteins_wide,
                     file.path( args$output_dir,
                                paste0(args$file_prefix, "_wide.xlsx")))
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Create long format output file")

norm_counts <- counts_rnorm.log.ruvIII_v1 %>%
  as.data.frame %>%
  rownames_to_column(args$row_id) %>%
  pivot_longer(cols = matches(args$group_pattern),
               names_to = args$sample_id,
               values_to = "log2norm")  %>%
  left_join(design_mat_cln, by = args$sample_id) %>%
  mutate(temp_id = !!sym(args$sample_id)) %>%
  separate(temp_id, sep = "_", into = c("sample_number", "group_pattern")) %>%
  group_by(!!sym(args$row_id), !!sym(args$group_id)) %>%
  arrange(!!sym(args$row_id), !!sym(args$group_id), sample_number) %>%
  mutate(replicate_number = paste0("log2norm.", row_number())) %>%
  ungroup %>%
  pivot_wider(id_cols = c(!!sym(args$row_id), !!sym(args$group_id)),
              names_from = replicate_number,
              values_from = log2norm)


raw_counts <- counts_filt %>%
  as.data.frame %>%
  rownames_to_column(args$row_id) %>%
  pivot_longer(cols = matches(args$group_pattern),
               names_to = args$sample_id,
               values_to = "raw") %>%
  left_join(design_mat_cln, by = args$sample_id) %>%
  mutate(temp_id = !!sym(args$sample_id)) %>%
  separate(temp_id, sep = "_", into = c("sample_number", "group_pattern")) %>%
  group_by(!!sym(args$row_id), !!sym(args$group_id)) %>%
  arrange(!!sym(args$row_id), !!sym(args$group_id), sample_number) %>%
  mutate(replicate_number = paste0("raw.", row_number())) %>%
  ungroup %>%
  pivot_wider(id_cols = c(!!sym(args$row_id), !!sym(args$group_id)),
              names_from = replicate_number,
              values_from = raw)

left_join_columns <- rlang::set_names(c(args$row_id, args$group_id ),
                                      c(args$row_id, "left_group"))

right_join_columns <- rlang::set_names(c(args$row_id, args$group_id ),
                                       c(args$row_id, "right_group"))

de_proteins_long <- selected_data %>%
  dplyr::filter(analysis_type == "RUV applied") %>%
  dplyr::select(-lqm, -colour, -analysis_type) %>%
  dplyr::mutate(expression = str_replace_all(expression, args$group_id, "")) %>%
  separate(expression, sep = "-", into = c("left_group", "right_group")) %>%
  left_join(norm_counts, by = left_join_columns) %>%
  left_join(norm_counts, by = right_join_columns,
            suffix = c(".left", ".right")) %>%
  left_join(raw_counts, by = left_join_columns) %>%
  left_join(raw_counts, by = right_join_columns,
            suffix = c(".left", ".right")) %>%
  arrange( comparison, q.mod, log2FC) %>%
  distinct()

#head(de_proteins_long)

vroom::vroom_write( de_proteins_long,
                    file.path( args$output_dir,
                               paste0(args$file_prefix, "_long.tsv")))

writexl::write_xlsx( de_proteins_long,
                     file.path( args$output_dir,
                                paste0(args$file_prefix, "_long.xlsx")))

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
te<-toc(quiet = TRUE)
loginfo("%f sec elapsed",te$toc-te$tic)
writeLines(capture.output(sessionInfo()), file.path(args$output_dir,"sessionInfo.txt"))

