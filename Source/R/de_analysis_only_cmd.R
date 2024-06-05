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
p_load(here)

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

parser <- add_option(parser, c("-c", "--config"), type = "character", default = "/Users/ignatiuspang/Workings/2024/e_33602_CSIRO_ThomasNguyen_20231114/scripts/proteomics/config_prot.ini",
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


parser <- add_option(parser, c( "--eBayes_trend"), type = "logical",
                     help = "logical, should an intensity-trend be allowed for the prior variance? Default is that the prior variance is constant.",
                     metavar = "double")

parser <- add_option(parser, c( "--eBayes_robust"), type = "logical",
                     help = "logical, should the estimation of df.prior and var.prior be robustified against outlier sample variances?.",
                     metavar = "double")



parser <- add_option(parser, "--group_pattern", type = "character",
                     help = "Regular expression pattern to identify columns with abundance values belonging to the experiment.",
                     metavar = "string")

parser <- add_option(parser, "--q_val_thresh", type = "double",
                     help = "q-value threshold below which a protein has statistically significant differetial expression",
                     metavar = "double")



parser <- add_option(parser, "--counts_table_file", type = "character",
                     help = "Input file with the protein abundance values",
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
  , "q_val_thresh"
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


args<-parseType(args,
                c( "eBayes_trend"
                  , "eBayes_robust"
                 ) , as.logical
                )

args<-parseString(args,
  c("group_pattern"
    ,"formula_string"
    ,"sample_id"
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




## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

loginfo("Read file with lists of experimental contrasts to test %s", args$contrasts_file)
captured_output<-capture.output(
contrasts_tbl <- vroom::vroom(args$contrasts_file, delim = "\t")
    ,type = "message"
)
logdebug(captured_output)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

loginfo("Read file with counts table %s", args$counts_table_file)
captured_output<-capture.output(
evidence_tbl_filt <- vroom::vroom(args$counts_table_file, delim = "\t") |>
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
  design_mat_cln <- vroom::vroom(args$design_matrix_file) |>
    as.data.frame() |>
    dplyr::mutate(!!rlang::sym(args$sample_id) := as.character(!!rlang::sym(args$sample_id)))
  ,type = "message"
)
logdebug(captured_output)

rownames(design_mat_cln) <- design_mat_cln |>
  pull(as.name(args$sample_id))

## Check that the sample ID and group ID name is found in the design matrix
if (length(which(c(args$sample_id, args$group_id) %in% colnames(design_mat_cln))) != 2) {
  logerror("Sample ID and group ID are not matching to the column names used in the design matrix.")
  q()
}


cols_for_analysis <- design_mat_cln |>
  pull(as.name(args$sample_id))

design_mat_updated <- design_mat_cln |>
  dplyr::filter( !is.na( !!rlang::sym( args$group_id )) )


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# heatmap( counts_rnorm.log.for.imputation, Rowv=NA, Colv=NA  )

loginfo("Count the number of values ")

counts_rnorm.log.ruvIII_v1 <- evidence_tbl_filt |>
  column_to_rownames(args$row_id) |>
  log2()


plot_num_of_values <- plotNumOfValuesNoLog(counts_rnorm.log.ruvIII_v1)

for( format_ext in args$plots_format) {
  file_name <- file.path(args$output_dir,paste0("num_of_values.",format_ext))
  captured_output <- capture.output(
    ggsave(filename = file_name, plot = plot_num_of_values, limitsize = FALSE)
    ,type = "message" )
  logdebug(captured_output)
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Compare the different experimental groups and obtain lists of differentially expressed proteins.")

list_rnorm.log.quant.ruv.r1 <- NA
myRes_rnorm.log.quant.ruv.r1 <- NA


  # requires statmod library
  list_rnorm.log.quant.ruv.r1 <- runTestsContrasts(counts_rnorm.log.ruvIII_v1,
                                                   contrast_strings = contrasts_tbl[, 1][[1]],
                                                   design_matrix = design_mat_updated,
                                                   formula_string = args$formula_string,
                                                   weights = NA,
                                                   treat_lfc_cutoff = as.double(args$treat_lfc_cutoff),
                                                   eBayes_trend = as.logical(args$eBayes_trend),
                                                   eBayes_robust = as.logical(args$eBayes_robust))

myRes_rnorm.log.quant.ruv.r1 <- list_rnorm.log.quant.ruv.r1$results

## This plot is used to check the mean-variance relationship of the expression data, after fitting a linear model.
pdf(file.path(args$output_dir, "plotSA_after_ruvIII.pdf" ))
plotSA(list_rnorm.log.quant.ruv.r1$fit.eb)
dev.off()

png(file.path(args$output_dir, "plotSA_after_ruvIII.png" ))
plotSA(list_rnorm.log.quant.ruv.r1$fit.eb)
dev.off()

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
selected_data <- getSignificantData(list_of_de_tables = list(myRes_rnorm.log.quant.ruv.r1),
                                    list_of_descriptions = list("RUV applied"),
                                    row_id = !!sym(args$row_id),
                                    p_value_column = p.mod,
                                    q_value_column = q.mod,
                                    fdr_value_column = fdr.mod,
                                    log_q_value_column = lqm,
                                    log_fc_column = logFC,
                                    comparison_column = comparison,
                                    expression_column = expression,
                                    facet_column = analysis_type,
                                    q_val_thresh = args$q_val_thresh) |>
  dplyr::rename(log2FC = "logFC")

## Write all the results in one single table
selected_data |>
  dplyr:::select(-colour, -lqm) |>
  vroom::vroom_write(file.path(args$output_dir, "lfc_qval_long.tsv"))

selected_data |>
  dplyr:::select(-colour, -lqm) |>
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

num_sig_de_molecules <- printCountDeGenesTable(list_of_de_tables = list(myRes_rnorm.log.quant.ruv.r1),
                                               list_of_descriptions = list( "RUV applied"),
                                               formula_string = "analysis_type ~ comparison")

for( format_ext in args$plots_format) {
  file_name<-file.path(args$output_dir,paste0("num_sda_entities_barplot.",format_ext))
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
norm_counts <- NA

counts_table_to_use <- counts_rnorm.log.ruvIII_v1


  norm_counts <- counts_table_to_use |>
    as.data.frame() |>
    set_colnames(paste0(colnames(counts_table_to_use), ".log2norm")) |>
    rownames_to_column(args$row_id)




de_proteins_wide <- selected_data |>
  dplyr::filter(analysis_type == "RUV applied") |>
  dplyr::select(-lqm, -colour, -analysis_type) |>
  pivot_wider(id_cols = c(!!sym(args$row_id)),
              names_from = c(comparison),
              names_sep = ":",
              values_from = c(log2FC, q.mod, p.mod)) |>
  left_join(norm_counts, by = args$row_id) |>
  dplyr::arrange(across(matches("q.mod"))) |>
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


counts_table_to_use <- counts_rnorm.log.ruvIII_v1


de_proteins_long <- createDeResultsLongFormat( lfc_qval_tbl = selected_data |>
                                                 dplyr::filter(analysis_type == "RUV applied") ,
                                               norm_counts_input_tbl = counts_table_to_use,
                                               raw_counts_input_tbl = counts_table_to_use,
                                               row_id = args$row_id,
                                               sample_id = args$sample_id,
                                               group_id = args$group_id,
                                               group_pattern = args$group_pattern,
                                               design_matrix_norm = design_mat_updated,
                                               design_matrix_raw =  design_mat_updated )

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

