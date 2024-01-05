#!/usr/bin/env Rscript

# Author(s): Ignatius Pang, Pablo Galaviz
# Email: cmri-bioinformatics@cmri.org.au
# Children’s Medical Research Institute, finding cures for childhood genetic diseases


##-------------------------------------

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
p_load(GGally)
p_load(ggpubr)
p_load(magrittr)
p_load(knitr)
p_load(rlang)
p_load(ggrepel)
p_load(gridExtra)

p_load(limma)
p_load(qvalue)
p_load(ruv)
p_load(mixOmics)

p_load(ProteomeRiver)
p_load(configr)
p_load(logging)
p_load(svglite)
p_load(Glimma)

##-------------------------------------

## Input Parameters
##-------------------------------------
tic()


command_line_options <- commandArgs(trailingOnly = TRUE)
#Note: options with default values are ignored in the configuration file parsing.

parser <- OptionParser(add_help_option = TRUE)


# input_dir <- "/home/ignatius/PostDoc/2022/PYROXD1_BMP_13/Results/N155S/de_proteins/"

parser <- add_option(parser, c("-d", "--debug"), action = "store_true", default = FALSE,
                     help = "Print debugging output  [default %default]")

parser <- add_option(parser, c("-s", "--silent"), action = "store_true", default = FALSE,
                     help = "Only print critical information to the console.  [default %default]")

parser <- add_option(parser, c("-n", "--no_backup"), action = "store_true", default = FALSE,
                     help = "Deactivate backup of previous run.  [default %default]")

parser <- add_option(parser, c("-c", "--config"), type = "character", default = "config_prot.ini",
                     help = "Configuration file.  [default %default]",
                     metavar = "string")

parser <- add_option(parser, c("-l", "--log_file"), type = "character", default = "output.log",
                     help = "Name of the logging file.  [default %default]",
                     metavar = "string")

parser <- add_option(parser, c( "--input_dir"), type = "character",
                     help = "Directory path for input data files.",
                     metavar = "string")

# output_dir <- "/home/ignatius/PostDoc/2022/PYROXD1_BMP_13/Results/N155S/publication_graphs"

parser <- add_option(parser, c("-o", "--output_dir"), type = "character",
                     help = "Directory path for all results files.",
                     metavar = "string")

# args$design_matrix_file <- "/home/ignatius/PostDoc/2022/PYROXD1_BMP_13/Source/N155S/design_matrix.tab"
parser <- add_option(parser, "--design_matrix_file", type = "character",
                     help = "Input file with the design matrix",
                     metavar = "string")

parser <- add_option(parser, "--before_avg_design_matrix_file", type = "character",
                     help = "Input file with the design matrix, before technical replicates were averaged and merged into one sample.",
                     metavar = "string")

parser <- add_option(parser, "--sample_id", type = "character",
                     help = "A string describing the sample ID. This must be a column that exists in the design matrix.",
                     metavar = "string")

parser <- add_option(parser, "--group_id", type = "character",
                     help = "A string describing the replicate group ID. This must be a column that exists in the design matrix.",
                     metavar = "string")

parser <- add_option(parser, "--row_id", type = "character",
                     help = "A string describing the row id.",
                     metavar = "string")

parser <- add_option(parser, "--q_val_thresh", type = "double",
                     help = "q-value threshold below which a protein has statistically significant differetial expression",
                     metavar = "double")

# de_proteins_long_file <- "/home/ignatius/PostDoc/2022/PYROXD1_BMP_13/Results/N155S/annot_proteins/de_proteins_long_annot.tsv"
parser <- add_option(parser, "--de_proteins_long_file", type = "character",
                     help = "File path for the de_proteins_long_annot.tsv annotation file, an output from annot_prots_cmd.R",
                     metavar = "string")

parser <- add_option(parser, "--plots_format", type = "character",
                     help = "A comma separated strings to indicate the fortmat output for the plots [pdf,png,svg].",
                     metavar = "string")


# top_x_gene_name <- 5
parser <- add_option(parser, "--top_x_gene_name", type = "integer",
                     help = "Print the top X gene names f",
                     metavar = "string")

parser <- add_option(parser, "--data_type", type = "character",
                     help = "Whether the analysis is proteomics or phosphoproteomics",
                     metavar = "string")

parser <- add_option(parser, "--log2fc_column", type = "character",
                     help = "The column which contains the log2FC value",
                     metavar = "string")


parser <- add_option(parser, "--fdr_column", type = "character",
                     help = "The column which contains the FDR value",
                     metavar = "string")


##-------------------------------------


#parse comand line arguments first.
args <- parse_args(parser)

#parse and merge the configuration file options.
if (args$config != "") {
  args <- config.list.merge(eval.config(file = args$config, config="publication_graphs"), args)
}

args <- setArgsDefault(args, "output_dir", as_func=as.character, default_val="publication_graphs" )
args <- setArgsDefault(args, "before_avg_design_matrix_file", as_func=as.character, default_val=NA )

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
  "q_val_thresh"
  , "sample_id"
  , "group_id"
  , "row_id"
  , "input_dir"
  , "output_dir"
))

testRequiredFiles(c(
    args$de_proteins_long_file
  , args$design_matrix_file
))

args <- setArgsDefault(args, "q_val_thresh", as_func=as.double, default_val=0.05 )
args <- setArgsDefault(args, "data_type", as_func=as.character, default_val="proteomics" )
args <- setArgsDefault(args, "log2fc_column", as_func=as.character, default_val="log2FC" )
args <- setArgsDefault(args, "fdr_column", as_func=as.character, default_val="q.mod" )

args<-parseType(args
                , c("q_val_thresh")
                , as.double)

args<-parseString(args,
                  c("sample_id"
                    , "group_id"
                    , "row_id"
                    , "file_prefix"
                    , "plots_format"
                    , "input_dir"
                    , "output_dir"
                  ))

if(isArgumentDefined(args,"plots_format"))
{
  args <- parseList(args,
                    c("plots_format"))
}else {
  logwarn("plots_format is undefined, default output set to pdf.")
  args$plots_format <- list("pdf")
}


createOutputDir(args$output_dir, args$no_backup)

##-------------------------------------


## Before normalization

counts_before_averaging <- NA
if(!is.na(args$before_avg_design_matrix_file)) {

  counts_before_averaging_file <- file.path( args$input_dir, "counts_after_normalization_before_imputation.tsv")
  counts_before_averaging <- vroom::vroom( counts_before_averaging_file  )
}

## After averaging (before RUVIII)
counts_rnorm.log.quant <- NA

file_a <- file.path(args$input_dir, "counts_after_normalization_and_imputation.tsv")
file_b <- file.path(args$input_dir, "counts_after_normalization_before_imputation_averaged.tsv")
file_c <- file.path(args$input_dir, "counts_after_normalization_before_imputation.tsv")

counts_file <- file_c
if( file.exists(file_a)) {
  counts_file <- file_a
} else if (file.exists( file_b)) {
  counts_file <- file_b
}

if(!file.exists(counts_file)) {
  logerror( paste0("Counts after normalization file not found: ", counts_file))
}

counts_rnorm.log.quant <- vroom::vroom(  counts_file)



## After averaging and RUVIII


counts_rnorm.log.ruvIII <- NA
if( file.exists(file.path(args$input_dir, "normalized_counts_after_ruv_remove_imputed.tsv") )) {
  counts_rnorm.log.ruvIII <- vroom::vroom( file.path(args$input_dir, "normalized_counts_after_ruv_remove_imputed.tsv")  )

} else if( file.exists(file.path(args$input_dir, "normalized_counts_after_ruv.tsv") )) {

  counts_rnorm.log.ruvIII <- vroom::vroom( file.path(args$input_dir, "normalized_counts_after_ruv.tsv")  )

} else {
   logerror( "Normalized counts after RUV file not found.")
}

createDirIfNotExists(args$output_dir)

##-------------------------------------

## Read Design Matrix
##-------------------------------------

design_mat_cln <- vroom::vroom(args$design_matrix_file) %>%
    as.data.frame() %>%
    dplyr::mutate(!!rlang::sym(args$sample_id) := as.character(!!rlang::sym(args$sample_id)))


rownames(design_mat_cln) <- design_mat_cln %>% pull(as.name(args$sample_id))

## design matrix when technical replicates has been averaged and merged into one

before_avg_design_mat_cln <- NA

if(!is.na(args$before_avg_design_matrix_file )) {
  before_avg_design_mat_cln <- vroom::vroom(args$before_avg_design_matrix_file) %>%
    as.data.frame() %>%
    dplyr::mutate(!!rlang::sym(args$sample_id) := as.character(!!rlang::sym(args$sample_id)))

  rownames(before_avg_design_mat_cln) <- before_avg_design_mat_cln %>% pull(as.name(args$sample_id))
}

##-------------------------------------
## RLE plots
##-------------------------------------

counts_rnorm.log.quant_mat <- counts_rnorm.log.quant %>%
  as.data.frame() %>%
  column_to_rownames(args$row_id)

counts_rnorm.log.ruvIII_mat <- counts_rnorm.log.ruvIII %>%
  as.data.frame() %>%
  column_to_rownames(args$row_id)


before_RUVIII_rle <- plotRle(t(as.matrix(counts_rnorm.log.quant_mat)),
        rowinfo = design_mat_cln[colnames(counts_rnorm.log.quant_mat),
                                 args$group_id]) +
  theme(axis.text.x = element_text(size = 13))   +
  theme(axis.text.y = element_text(size = 13))  +
  theme(axis.title.x = element_text(size = 12))  +
  theme(axis.title.y = element_text(size = 12))  +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 12)) +
  xlab("Samples")

before_RUVIII_rle

createDirectoryIfNotExists(file.path( args$output_dir, "RLE"))

file_name_part <- file.path( args$output_dir, "RLE", "before_RUVIII_rle.")
gg_save_logging ( before_RUVIII_rle, file_name_part, args$plots_format)


##
after_RUVIII_rle <- plotRle(t(as.matrix(counts_rnorm.log.ruvIII_mat)),
        rowinfo = design_mat_cln[colnames(counts_rnorm.log.ruvIII_mat),
                                 args$group_id]) +
  theme(axis.text.x = element_text(size = 13))   +
  theme(axis.text.y = element_text(size = 13))  +
  theme(axis.title.x = element_text(size = 12))  +
  theme(axis.title.y = element_text(size = 12))  +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 12)) +
  xlab("Samples")

after_RUVIII_rle

file_name_part <- file.path( args$output_dir, "RLE", "after_RUVIII_rle.")
gg_save_logging ( after_RUVIII_rle, file_name_part, args$plots_format)


if(!is.na(args$before_avg_design_matrix_file )) {
  counts_before_averaging_mat <- counts_before_averaging %>%
    as.data.frame() %>%
    column_to_rownames("uniprot_acc")

  counts_before_averaging_rle <- plotRle(t(as.matrix(counts_before_averaging_mat)),
                              rowinfo = before_avg_design_mat_cln[colnames(counts_before_averaging_mat),
                                                       args$group_id]) +
    theme(axis.text.x = element_text(size = 13))   +
    theme(axis.text.y = element_text(size = 13))  +
    theme(axis.title.x = element_text(size = 12))  +
    theme(axis.title.y = element_text(size = 12))  +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 12)) +
    xlab("Samples")

  counts_before_averaging_rle

  file_name_part <- file.path( args$output_dir, "RLE", "counts_before_averaging_rle.")
  gg_save_logging ( counts_before_averaging_rle, file_name_part, args$plots_format)

}


##-------------------------------------
## Hierarchical clustering


before_hclust <- hclust(  dist( t(as.matrix(counts_rnorm.log.quant_mat) )  ) )
plot(before_hclust)



after_hclust <- hclust(  dist( t(as.matrix(counts_rnorm.log.ruvIII_mat) )  ), method="ward.D2" )
plot(after_hclust)


## Volcano plots
##-------------------------------------
#
# gene_names_data <- vroom::vroom(args$de_proteins_long_file) %>%
#   mutate(gene_name = str_split(gene_names, ":") %>% purrr::map_chr(1)) %>%
#   relocate( gene_name, .before="gene_names")
#
# show_gene_name <- gene_names_data %>%
#   group_by(comparison) %>%
#   arrange(comparison, !!sym(args$fdr_column)) %>%
#   mutate(row_id = row_number()) %>%
#   ungroup()  %>%
#   dplyr::filter( row_id <= args$top_x_gene_name & !!sym(args$fdr_column) < args$q_val_thresh) %>%
#   dplyr::select( comparison, !!rlang::sym(args$row_id), gene_name)

selected_data <- vroom::vroom( file.path(args$de_proteins_long_file ) )  %>% # args$input_dir, "lfc_qval_long.tsv"
    mutate( lqm = -log10(!!sym(args$fdr_column)))  %>%
    dplyr::mutate(label = case_when(abs(!!sym(args$log2fc_column)) >= 1 & !!sym(args$fdr_column) >= args$q_val_thresh ~ "Not sig., logFC >= 1",
                                   abs(!!sym(args$log2fc_column)) >= 1 & !!sym(args$fdr_column) < args$q_val_thresh ~ "Sig., logFC >= 1",
                                   abs(!!sym(args$log2fc_column)) < 1 & !!sym(args$fdr_column) < args$q_val_thresh ~ "Sig., logFC < 1",
                                   TRUE ~ "Not sig.")) %>%
    dplyr::mutate(colour = case_when(abs(!!sym(args$log2fc_column)) >= 1 & !!sym(args$fdr_column) >= args$q_val_thresh ~ "orange",
                                     abs(!!sym(args$log2fc_column)) >= 1 & !!sym(args$fdr_column) < args$q_val_thresh ~ "purple",
                                     abs(!!sym(args$log2fc_column)) < 1 & !!sym(args$fdr_column) < args$q_val_thresh ~ "blue",
                                     TRUE ~ "black")) %>%
    dplyr::mutate(colour = factor(colour, levels = c("black", "orange", "blue", "purple"))) # %>%
  # left_join( show_gene_name, by = c("comparison" = "comparison",
  #                                   "uniprot_acc" = "uniprot_acc")) %>%
  # dplyr::mutate( gene_name = case_when( !!sym(args$fdr_column) < args$q_val_thresh ~ gene_name,
  #                                            TRUE ~ NA_character_) )

#
# plotOneVolcanoWithGeneName <- function( input_data, input_title) {
#   volcano_plot <-  input_data %>%
#     ggplot(aes(y = lqm, x = !!sym(args$log2fc_column), label=gene_name)) +
#     geom_point(aes(col = colour)) +
#     scale_colour_manual(values = c(levels(selected_data$colour)),
#                         labels = c(paste0("Not sig., logFC > ",
#                                           1),
#                                    paste0("Sig,, logFC >= ",
#                                           1),
#                                    paste0("Sig., logFC <",
#                                           1),
#                                    "Not Sign,")) +
#     geom_vline(xintercept = 1, colour = "black", size = 0.2) +
#     geom_vline(xintercept = -1, colour = "black", size = 0.2) +
#     geom_hline(yintercept = -log10(args$q_val_thresh)) +
#     geom_text_repel( size  = 7, show.legend=FALSE) +
#     theme_bw() +
#     xlab("Log fold-change") +
#     ylab(expression(-log[10] ~ q ~ value)) +
#     labs(title = input_title)+  # Remove legend title
#     theme(legend.title = element_blank()) +
#     # theme(legend.position = "none")  +
#     theme(axis.text.x = element_text(size = 13))   +
#     theme(axis.text.y = element_text(size = 13))  +
#     theme(axis.title.x = element_text(size = 12))  +
#     theme(axis.title.y = element_text(size = 12))  +
#     theme(plot.title = element_text(size = 12)) +
#     theme(legend.text = element_text(size = 12)) # +
#     # theme(legend.title = element_text(size = 12))
#
#     volcano_plot
# }

createDirectoryIfNotExists(file.path( args$output_dir, "Volcano_Plots"))


list_of_volcano_plots <- selected_data %>%
  group_by( comparison) %>%
  nest() %>%
  ungroup() %>%
  mutate( title = paste( comparison)) %>%
  #mutate( data = purrr::map (data, ~{ (.) %>% mutate( !!sym(args$log2fc_column)_edited = 2^!!sym(args$log2fc_column))})) %>%
  mutate( plot = purrr:::map2( data, title, \(x,y) { plotOneVolcano(x, y,   log_fc_column = !!sym(args$log2fc_column))}) )

list_of_volcano_plots %>% pull(plot)

purrr::walk2( list_of_volcano_plots %>% pull(title),
              list_of_volcano_plots %>% pull(plot),
              ~{file_name_part <- file.path( args$output_dir, "Volcano_Plots", paste0(.x, "."))
              gg_save_logging ( .y, file_name_part, args$plots_format)} )

ggsave(
  filename = file.path(args$output_dir, "Volcano_Plots", "list_of_volcano_plots.pdf" ),
  plot = marrangeGrob( (list_of_volcano_plots  %>% pull(plot)), nrow=1, ncol=1),
  width = 7, height = 7
)


#
# list_of_volcano_plots_gene_labels <- selected_data %>%
#   group_by( analysis_type, comparison) %>%
#   nest() %>%
#   ungroup() %>%
#   mutate( title = paste( analysis_type, comparison)) %>%
#   mutate( plot = purrr:::map2( data, title, ~plotOneVolcanoWithGeneName(.x, .y)) )
#
# list_of_volcano_plots_gene_labels %>%
#   pull(plot)
#
# purrr::walk2( list_of_volcano_plots_gene_labels %>% pull(title),
#               list_of_volcano_plots_gene_labels %>% pull(plot),
#               ~{file_name_part <- file.path( args$output_dir, "Volcano_Plots", paste0(.x, "."))
#               gg_save_logging ( .y, file_name_part, args$plots_format)} )
#
# ggsave(
#   filename = file.path(args$output_dir, "Volcano_Plots", "list_of_volcano_plots_gene_labels.pdf" ),
#   plot = marrangeGrob( (list_of_volcano_plots_gene_labels  %>% pull(plot)), nrow=1, ncol=1),
#   width = 7, height = 7
# )


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


num_sig_de_molecules <- selected_data %>%
  dplyr::mutate(status = case_when( !!sym(args$fdr_column) >= 0.05 ~ "Not significant",
                                    !!sym(args$log2fc_column) >= 0 & !!sym(args$fdr_column) < 0.05 ~ "Significant and Up",
                                    !!sym(args$log2fc_column) < 0 & !!sym(args$fdr_column) < 0.05 ~ "Significant and Down",
                                   TRUE ~ "Not significant")) %>%
group_by( comparison,  status) %>% # expression, analysis_type,
  summarise(counts = n()) %>%
  ungroup()

formula_string <- ". ~ comparison"
num_sig_de_genes_barplot <- num_sig_de_molecules %>%
  dplyr::filter(status != "Not significant") %>%
  ggplot(aes(x = status, y = counts)) +
  geom_bar(stat = "identity") +
  geom_text(stat = 'identity', aes(label = counts), vjust = -0.5) +
  theme(axis.text.x = element_text(angle = 90))  +
  facet_grid(as.formula(formula_string))

num_sig_de_genes_barplot

createDirIfNotExists(file.path(args$output_dir, "NumSigDeMolecules"))

vroom::vroom_write( num_sig_de_molecules,
                    file.path(args$output_dir, "NumSigDeMolecules", "num_sig_de_molecules.tab" ) )

num_of_comparison <- num_sig_de_molecules |>
  distinct(comparison) |>
  nrow()

ggsave(filename = file.path(args$output_dir, "NumSigDeMolecules", "num_sig_de_molecules.png" ),
       plot = num_sig_de_genes_barplot,
       height = 10,
       width = (num_of_comparison + 2) *7/6 )


##-------------------------------------
## PCA plots
##-------------------------------------

counts_before_averaging_mat <- NA
if(!is.na(args$before_avg_design_matrix_file)) {
  counts_before_averaging_mat <- counts_before_averaging  %>%
    as.data.frame() %>%
    column_to_rownames(args$row_id)
}

before_ruvIII_pca <- plotPca( counts_rnorm.log.quant_mat,
                              design_matrix = design_mat_cln,
                              sample_id_column = !!rlang::sym(args$sample_id),
                              group_column = !!rlang::sym(args$group_id),
                              title = "Before RUVIII",
                              geom.text.size = 7) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 12))

before_ruvIII_pca

createDirectoryIfNotExists(file.path(args$output_dir, "PCA"))

file_name_part <- file.path( args$output_dir, "PCA", "before_ruvIII_pca." )

gg_save_logging ( before_ruvIII_pca, file_name_part, args$plots_format)

after_ruvIII_pca <- plotPca( counts_rnorm.log.ruvIII_mat,
                             design_matrix = design_mat_cln,
                             sample_id_column = !!rlang::sym(args$sample_id),
                             group_column = !!rlang::sym(args$group_id),
                             title = "After RUVIII", geom.text.size = 7) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12))   +
  theme(axis.text.y = element_text(size = 12))  +
  theme(axis.title.x = element_text(size = 12))  +
  theme(axis.title.y = element_text(size = 12))  +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 12))

after_ruvIII_pca

file_name_part <- file.path( args$output_dir, "PCA", "after_ruvIII_pca.")
gg_save_logging ( after_ruvIII_pca, file_name_part, args$plots_format)


## Labeling using the groups
before_ruvIII_pca <- plotPca( counts_rnorm.log.quant_mat,
                              design_matrix = design_mat_cln,
                              sample_id_column = !!rlang::sym(args$sample_id),
                              group_column = !!rlang::sym(args$group_id),
                              label_column = !!rlang::sym(args$group_id),
                              title = "Before RUVIII",
                              geom.text.size = 7) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 12))

before_ruvIII_pca

createDirectoryIfNotExists(file.path(args$output_dir, "PCA"))

file_name_part <- file.path( args$output_dir, "PCA", "before_ruvIII_pca_group_labels." )

gg_save_logging ( before_ruvIII_pca, file_name_part, args$plots_format)

after_ruvIII_pca <- plotPca( counts_rnorm.log.ruvIII_mat,
                             design_matrix = design_mat_cln,
                             sample_id_column = !!rlang::sym(args$sample_id),
                             group_column = !!rlang::sym(args$group_id),
                             label_column = !!rlang::sym(args$group_id),
                             title = "After RUVIII", geom.text.size = 7) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12))   +
  theme(axis.text.y = element_text(size = 12))  +
  theme(axis.title.x = element_text(size = 12))  +
  theme(axis.title.y = element_text(size = 12))  +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 12))

after_ruvIII_pca

file_name_part <- file.path( args$output_dir, "PCA", "after_ruvIII_pca_group_labels.")
gg_save_logging ( after_ruvIII_pca, file_name_part, args$plots_format)



## No sample labels
before_ruvIII_pca_no_labels <- plotPca( counts_rnorm.log.quant_mat,
                              design_matrix = design_mat_cln %>%
                                mutate( sample_labels = ""),
                              sample_id_column = !!rlang::sym(args$sample_id),
                              group_column = !!rlang::sym(args$group_id),
                              label_column = sample_labels,
                              title = "Before RUVIII",
                              geom.text.size = 7) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 12))

before_ruvIII_pca_no_labels


file_name_part <- file.path( args$output_dir, "PCA", "before_ruvIII_no_sample_labels." )
gg_save_logging ( before_ruvIII_pca_no_labels, file_name_part, args$plots_format)




after_ruvIII_pca_no_labels <- plotPca( counts_rnorm.log.ruvIII_mat,
                                       design_matrix = design_mat_cln %>%
                                         mutate( sample_labels = ""),
                                       sample_id_column = !!rlang::sym(args$sample_id),
                                       group_column = !!rlang::sym(args$group_id),
                                       label_column = sample_labels,
                                       title = "After RUVIII", geom.text.size = 7) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12))   +
  theme(axis.text.y = element_text(size = 12))  +
  theme(axis.title.x = element_text(size = 12))  +
  theme(axis.title.y = element_text(size = 12))  +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 12))

after_ruvIII_pca_no_labels

file_name_part <- file.path( args$output_dir, "PCA", "after_ruvIII_no_sample_labels.")
gg_save_logging ( after_ruvIII_pca_no_labels, file_name_part, args$plots_format)



if(!is.na(args$before_avg_design_matrix_file )) {


  counts_before_averaging_pca <- plotPca( counts_before_averaging_mat,
                                   design_matrix = before_avg_design_mat_cln,
                                   sample_id_column = !!rlang::sym(args$sample_id),
                                   group_column = !!rlang::sym(args$group_id),
                                   title = "After RUVIII", geom.text.size = 7) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12))   +
    theme(axis.text.y = element_text(size = 12))  +
    theme(axis.title.x = element_text(size = 12))  +
    theme(axis.title.y = element_text(size = 12))  +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 12))

  counts_before_averaging_pca

  file_name_part <- file.path( args$output_dir, "PCA", "counts_before_averaging_pca.")
  gg_save_logging ( counts_before_averaging_pca, file_name_part, args$plots_format)

}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Interactive Volcano Plot


if(   args$data_type  == "proteomics"
      && file.exists(file.path( args$input_dir, "fit.eb.RDS"))
      && file.exists( args$de_proteins_long_file ) ) {

  loginfo("Create proteomics interactive volcano plot")

  volcano_plot_tab <- vroom::vroom( args$de_proteins_long_file  )  %>%
    mutate( lqm = -log10(!!sym(args$fdr_column)))  |>
    dplyr::mutate(label = case_when(abs(!!sym(args$log2fc_column)) >= 1 & !!sym(args$fdr_column) >= args$q_val_thresh ~ "Not sig., logFC >= 1",
                                    abs(!!sym(args$log2fc_column)) >= 1 & !!sym(args$fdr_column) < args$q_val_thresh ~ "Sig., logFC >= 1",
                                    abs(!!sym(args$log2fc_column)) < 1 & !!sym(args$fdr_column) < args$q_val_thresh ~ "Sig., logFC < 1",
                                    TRUE ~ "Not sig.")) |>
    dplyr::mutate(colour = case_when(abs(!!sym(args$log2fc_column)) >= 1 & !!sym(args$fdr_column) >= args$q_val_thresh ~ "orange",
                                     abs(!!sym(args$log2fc_column)) >= 1 & !!sym(args$fdr_column) < args$q_val_thresh ~ "purple",
                                     abs(!!sym(args$log2fc_column)) < 1 & !!sym(args$fdr_column) < args$q_val_thresh ~ "blue",
                                     TRUE ~ "black")) |>
    dplyr::mutate(gene_name = str_split(UNIPROT_GENENAME, " |:" ) |> purrr::map_chr(1)  ) |>
    dplyr::mutate(best_uniprot_acc = str_split(!!sym(args$row_id), ":" ) |> purrr::map_chr(1)  ) |>
    dplyr::mutate(analysis_type = comparison)  |>
    dplyr::rename( PROTEIN_NAMES = "PROTEIN-NAMES") |>
    dplyr::select( best_uniprot_acc, lqm, !!sym(args$fdr_column), p.mod, !!sym(args$log2fc_column), comparison, label, colour,  gene_name, `PROTEIN_NAMES`)   |>
    dplyr::mutate( my_alpha = case_when ( gene_name !=  "" ~ 1
                                          , TRUE ~ 0.5))


  r_obj <- readRDS( file.path( args$input_dir, "fit.eb.RDS") )

  # ncol(r_obj$coefficients)
  # colnames(r_obj$coefficients)

  output_dir <- file.path( args$output_dir
                           ,  "Interactive_Volcano_Plots")

  createDirIfNotExists(output_dir)


  purrr::walk( seq_len( ncol(r_obj$coefficients))
               , \(coef) { # print(coef)
                 ProteomeRiver::getGlimmaVolcanoProteomics( r_obj
                                                            , coef = coef
                                                            , volcano_plot_tab  = volcano_plot_tab
                                                            , uniprot_column = best_uniprot_acc
                                                            , gene_name_column = gene_name
                                                            , display_columns = c( "best_uniprot_acc",  "PROTEIN_NAMES"   )
                                                            , output_dir = output_dir ) } )

}



if(   args$data_type  == "phosphoproteomics"
      && file.exists(file.path( args$input_dir, "fit.eb.RDS"))
      && file.exists( args$de_proteins_long_file ) ) {

  merge_residue_position_lists <- function(residue, position)  {

    residues_list <- str_split(residue, ";")[[1]]
    positions_list <- str_split(position, ";")[[1]]

    #print( paste(  residue, position ) )

    if( length(residues_list) != length(positions_list)) {
      stop( paste( "Length not equal", residue, position ) )
    }

    purrr::map2_chr( residues_list
                     , positions_list
                     , function(r, p){ paste0(r, p) }) |>
      paste(collapse=";")
  }
  # merge_residue_position_lists ("S;T",  "11;12")

  clean_first_positiion <- function( position) {
    first_position <- str_split(position, "\\|" ) |>
      purrr::map_chr(1) |>
      str_replace( "\\(", "") |>
      str_replace( "\\)", "")

    first_position
  }
  # clean_first_positiion("(520;526)|(562;568)|(706;712)|(724;730)|(736;742)")


  loginfo("Create phosphoproteomics interactive volcano plot")


  de_proteins_long_tbl <- vroom::vroom( args$de_proteins_long_file  )

  volcano_plot_colour_points <- de_proteins_long_tbl  %>%
    mutate( lqm = -log10(!!sym(args$fdr_column)))  |>
    dplyr::mutate(label = case_when(abs(!!sym(args$log2fc_column)) >= 1 & !!sym(args$fdr_column) >= args$q_val_thresh ~ "Not sig., logFC >= 1",
                                    abs(!!sym(args$log2fc_column)) >= 1 & !!sym(args$fdr_column) < args$q_val_thresh ~ "Sig., logFC >= 1",
                                    abs(!!sym(args$log2fc_column)) < 1 & !!sym(args$fdr_column) < args$q_val_thresh ~ "Sig., logFC < 1",
                                    TRUE ~ "Not sig.")) |>
    dplyr::mutate(colour = case_when(abs(!!sym(args$log2fc_column)) >= 1 & !!sym(args$fdr_column) >= args$q_val_thresh ~ "orange",
                                     abs(!!sym(args$log2fc_column)) >= 1 & !!sym(args$fdr_column) < args$q_val_thresh ~ "purple",
                                     abs(!!sym(args$log2fc_column)) < 1 & !!sym(args$fdr_column) < args$q_val_thresh ~ "blue",
                                     TRUE ~ "black")) |>
    dplyr::mutate(analysis_type = comparison)  |>
    dplyr::mutate(gene_name = str_split(UNIPROT_GENENAME, " |:" ) |> purrr::map_chr(1)  ) |>
    dplyr::mutate(best_uniprot_acc = str_split(!!sym(args$row_id), ":" ) |> purrr::map_chr(1)  ) |>
    mutate( first_position = purrr::map_chr(position, clean_first_positiion)) |>
    mutate( merged_sites_residues = purrr::map2_chr( residue
                                                     ,  first_position
                                                     , \(residue, position){merge_residue_position_lists(residue, position)}) ) |>
    mutate( sites_id_short = paste0( gene_name, ":", merged_sites_residues, ":", best_uniprot_acc) ) |>
    relocate( sites_id_short, .after="sites_id")

  volcano_plot_tab <- volcano_plot_colour_points |>
    dplyr::rename( PROTEIN_NAMES = "PROTEIN-NAMES") |>
    dplyr::select( sites_id, sites_id_short, best_uniprot_acc, lqm, !!sym(args$fdr_column), !!sym(args$log2fc_column), comparison, label, colour,  gene_name, sequence, `PROTEIN_NAMES`)   |>
    dplyr::mutate( my_alpha = case_when ( gene_name !=  "" ~ 1
                                          , TRUE ~ 0.5))

  r_obj <- readRDS( file.path( args$input_dir, "fit.eb.RDS") )

  # ncol(r_obj$coefficients)
  # colnames(r_obj$coefficients)

  output_dir <- file.path( args$output_dir
                           ,  "Interactive_Volcano_Plots")

  createDirIfNotExists(output_dir)

  # print(paste("nrow = ", nrow(r_obj@.Data[[1]])))
  # print(head(best_uniprot_acc))



  purrr::walk( seq_len( ncol(r_obj$coefficients))
               , \(coef) { # print(coef)
                 ProteomeRiver::getGlimmaVolcanoPhosphoproteomics( r_obj
                                                                   , coef = coef
                                                                   , volcano_plot_tab  = volcano_plot_tab
                                                                   , sites_id_column = sites_id
                                                                   , sites_id_display_column = sites_id_short
                                                                   , display_columns = c(  "sequence", "PROTEIN_NAMES" )
                                                                   , output_dir = output_dir ) } )

}





## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
te<-toc(quiet = TRUE)
loginfo("%f sec elapsed",te$toc-te$tic)
writeLines(capture.output(sessionInfo()), file.path(args$output_dir,"sessionInfo.txt"))
