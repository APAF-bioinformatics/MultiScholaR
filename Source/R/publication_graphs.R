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

parser <- add_option(parser, "--avg_design_matrix_file", type = "character",
                     help = "Input file with the design matrix, where technical replicates has already been averaged and merged into one sample.",
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

##-------------------------------------


#parse comand line arguments first.
args <- parse_args(parser)

#parse and merge the configuration file options.
if (args$config != "") {
  args <- config.list.merge(eval.config(file = args$config, config="publication_graphs"), args)
}

args <- setArgsDefault(args, "output_dir", as_func=as.character, default_val="publication_graphs" )
args <- setArgsDefault(args, "avg_design_matrix_file", as_func=as.character, default_val=NA )

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
  , "q_val_thresh"
))

testRequiredFiles(c(
    args$de_proteins_long_file
  , args$design_matrix_file
))

args <- setArgsDefault(args, "q_val_thresh", as_func=as.double, default_val=0.05 )

args<-parseType(args,
                c("q_val_thresh"
                )  ,as.double)

args<-parseString(args,
                  c("sample_id"
                    ,"group_id"
                    ,"row_id"
                    ,"file_prefix"
                    ,"plots_format"
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


counts_file <- file.path(args$input_dir, "counts_after_median_scaling_before_imputation.tsv")
if( !file.exists(counts_file)) {
  counts_file <- file.path(args$input_dir, "counts_after_normalization_before_imputation.tsv")
}

if(!file.exists(counts_file)) {
  logerror( paste0("Counts after normalization file not found: ", counts_file))
}

counts_rnorm.log.quant <- vroom::vroom(  counts_file)
counts_rnorm.log.ruvIII <- vroom::vroom( file.path(args$input_dir, "normalized_counts_after_ruv.tsv")  )
counts_rnorm.log.ruvIII.avg <- NA
if(!is.na(args$avg_design_matrix_file)) {
  counts_rnorm.log.ruvIII.avg <- vroom::vroom( file.path(args$input_dir, "normalized_counts_after_ruv_replicates_averaged.tsv")  )
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

avg_design_mat_cln <- NA

if(!is.na(args$avg_design_matrix_file )) {
  avg_design_mat_cln <- vroom::vroom(args$avg_design_matrix_file) %>%
    as.data.frame() %>%
    dplyr::mutate(!!rlang::sym(args$sample_id) := as.character(!!rlang::sym(args$sample_id)))

  rownames(avg_design_mat_cln) <- avg_design_mat_cln %>% pull(as.name(args$sample_id))
}

##-------------------------------------
## RLE plots
##-------------------------------------

counts_rnorm.log.quant_mat <- counts_rnorm.log.quant %>%
  column_to_rownames(args$row_id)

counts_rnorm.log.ruvIII_mat <- counts_rnorm.log.ruvIII %>%
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


if(!is.na(args$avg_design_matrix_file )) {


  after_RUVIII_avg_rle <- plotRle(t(as.matrix(counts_rnorm.log.ruvIII.avg_mat)),
                              rowinfo = avg_design_mat_cln[colnames(counts_rnorm.log.ruvIII.avg_mat),
                                                       args$group_id]) +
    theme(axis.text.x = element_text(size = 13))   +
    theme(axis.text.y = element_text(size = 13))  +
    theme(axis.title.x = element_text(size = 12))  +
    theme(axis.title.y = element_text(size = 12))  +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 12)) +
    xlab("Samples")

  after_RUVIII_avg_rle

  file_name_part <- file.path( args$output_dir, "RLE", "after_RUVIII_avg_rle.")
  gg_save_logging ( after_RUVIII_avg_rle, file_name_part, args$plots_format)

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
#   arrange(comparison, q.mod) %>%
#   mutate(row_id = row_number()) %>%
#   ungroup()  %>%
#   dplyr::filter( row_id <= args$top_x_gene_name & q.mod < args$q_val_thresh) %>%
#   dplyr::select( comparison, !!rlang::sym(args$row_id), gene_name)

selected_data <- vroom::vroom( file.path(args$input_dir, "lfc_qval_long.tsv") )  %>%
    mutate( lqm = -log10(q.mod))  %>%
    dplyr::mutate(label = case_when(abs(log2FC) >= 1 & q.mod >= args$q_val_thresh ~ "Not sig., logFC >= 1",
                                   abs(log2FC) >= 1 & q.mod < args$q_val_thresh ~ "Sig., logFC >= 1",
                                   abs(log2FC) < 1 & q.mod < args$q_val_thresh ~ "Sig., logFC < 1",
                                   TRUE ~ "Not sig.")) %>%
    dplyr::mutate(colour = case_when(abs(log2FC) >= 1 & q.mod >= args$q_val_thresh ~ "orange",
                                     abs(log2FC) >= 1 & q.mod < args$q_val_thresh ~ "purple",
                                     abs(log2FC) < 1 & q.mod < args$q_val_thresh ~ "blue",
                                     TRUE ~ "black")) %>%
    dplyr::mutate(colour = factor(colour, levels = c("black", "orange", "blue", "purple"))) # %>%
  # left_join( show_gene_name, by = c("comparison" = "comparison",
  #                                   "uniprot_acc" = "uniprot_acc")) %>%
  # dplyr::mutate( gene_name = case_when( q.mod < args$q_val_thresh ~ gene_name,
  #                                            TRUE ~ NA_character_) )

#
# plotOneVolcanoWithGeneName <- function( input_data, input_title) {
#   volcano_plot <-  input_data %>%
#     ggplot(aes(y = lqm, x = log2FC, label=gene_name)) +
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
  group_by( analysis_type, comparison) %>%
  nest() %>%
  ungroup() %>%
  mutate( title = paste( analysis_type, comparison)) %>%
  #mutate( data = purrr::map (data, ~{ (.) %>% mutate( log2FC_edited = 2^log2FC)})) %>%
  mutate( plot = purrr:::map2( data, title, ~plotOneVolcano(.x, .y,   log_fc_column = log2FC)) )

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
  dplyr::mutate(status = case_when( q.mod >= 0.05 ~ "Not significant",
                                    log2FC >= 0 & q.mod < 0.05 ~ "Significant and Up",
                                    log2FC < 0 & q.mod < 0.05 ~ "Significant and Down",
                                   TRUE ~ "Not significant")) %>%
group_by( comparison, expression, analysis_type, status) %>%
  summarise(counts = n()) %>%
  ungroup()

formula_string <- "analysis_type ~ comparison"
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


ggsave(filename = file.path(args$output_dir, "NumSigDeMolecules", "num_sig_de_molecules.png" ),
       plot = num_sig_de_genes_barplot,
       height = 10,
       width = 7)


##-------------------------------------
## PCA plots
##-------------------------------------

counts_rnorm.log.ruvIII.avg_mat <- NA
if(!is.na(args$avg_design_matrix_file)) {
  counts_rnorm.log.ruvIII.avg_mat <- counts_rnorm.log.ruvIII.avg  %>%
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



if(!is.na(args$avg_design_matrix_file )) {


  after_ruvIII_avg_pca <- plotPca( counts_rnorm.log.ruvIII.avg_mat,
                                   design_matrix = avg_design_mat_cln,
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

  after_ruvIII_avg_pca

  file_name_part <- file.path( args$output_dir, "PCA", "after_ruvIII_avg_pca.")
  gg_save_logging ( after_ruvIII_avg_pca, file_name_part, args$plots_format)

}



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
te<-toc(quiet = TRUE)
loginfo("%f sec elapsed",te$toc-te$tic)
writeLines(capture.output(sessionInfo()), file.path(args$output_dir,"sessionInfo.txt"))
