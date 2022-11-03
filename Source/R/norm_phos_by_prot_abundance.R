#!/usr/bin/env Rscript

# Author(s): Ignatius Pang, Pablo Galaviz
# Email: cmri-bioinformatics@cmri.org.au
# Children’s Medical Research Institute, finding cures for childhood genetic diseases


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Test if BioManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
   BiocManager::install(version = "3.12")
}

# load pacman package manager
if(!require(pacman)){
    install.packages("pacman")
    library(pacman)
}

p_load(tidyverse)
p_load(plotly)
p_load(vroom)
p_load(ggplot2)
p_load(ggpubr)
p_load(tictoc)


p_load(qvalue)
p_load(knitr)

p_load(magrittr)
p_load(optparse)
p_load(ProteomeRiver)
p_load(configr)
p_load(logging)

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
tic()


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

command_line_options <- commandArgs(trailingOnly = TRUE)
parser <- OptionParser(add_help_option = TRUE)
#Note: options with default values are ignored in the configuration file parsing.
parser <- add_option(parser, c("-d", "--debug"), action = "store_true", default = FALSE,
                     help = "Print debugging output")

parser <- add_option(parser, c("-s", "--silent"), action = "store_true", default = FALSE,
                     help = "Only print critical information to the console.")

parser <- add_option(parser, c("-n", "--no_backup"), action = "store_true", default = FALSE,
                     help = "Deactivate backup of previous run.")

parser <- add_option(parser, c("-c","--config"), type = "character", default = "config_phos_desch.ini", dest = "config",
                     help = "Configuration file.",
                     metavar = "string")

parser <- add_option(parser, c("-o","--output_dir"), type = "character", dest = "output_dir",
                     help = "Directory path for all results files.",
                     metavar = "string")

parser <- add_option(parser, c("-t","--tmp_dir"), type = "character", default = "cache", dest = "tmp_dir",
                     help = "Directory path for temporary files.",
                     metavar = "string")

parser <- add_option(parser, c("-l","--log_file"), type = "character", default = "output.log", dest = "log_file",
                     help = "Name of the logging file.",
                     metavar = "string")

#Options without a default value have the following priority: configuration file < command line argument
parser <- add_option(parser, c("--proteins_file"), type="character", dest = "proteins_file",
                     help="File with table of differentially expressed proteins.",
                     metavar="string")

parser <- add_option(parser, c("--phospho_file"), type="character", dest = "phospho_file",
                     help="File with table of differentially abundant phosphosites.",
                     metavar="string")

parser <- add_option(parser, "--plots_format", type = "character",
                     help = "A comma separated strings to indicate the fortmat output for the plots [pdf,png,svg].",
                     metavar = "string")

#parse command line arguments first.
args <- parse_args(parser)

#parse and merge the configuration file options.
if (args$config != "") {
  args <- config.list.merge(eval.config(file = args$config, config = "norm_phos_by_prot_abundance"), args)
}

args <- setArgsDefault(args, "output_dir", as_func=as.character, default_val="norm_phos_by_prot_abundance" )

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


testRequiredArguments(args, c(
  "proteins_file"
  ,"phospho_file"
))

testRequiredFiles(c(
  args$proteins_file
  ,args$phospho_file))

args<-parseString(args,c("plots_format"))

if(isArgumentDefined(args,"plots_format"))
{
  args <- parseList(args,c("plots_format"))
}else {
  logwarn("plots_format is undefined, default output set to pdf.")
  args$plots_format <- list("pdf")
}


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo ("Read phosphopeptides abundance data.")
captured_output<-capture.output(
phospho_tbl_orig <- vroom::vroom(args$phospho_file)
  ,type = "message"
)
logdebug(captured_output)


list_of_phospho_columns <- c("sites_id",
                             "uniprot_acc",
                             "gene_name",
                             "position",
                             "residue",
                             "sequence",
                             "q.mod",
                             "fdr.mod",
                             "p.mod",
                             "log2FC",
                             "comparison",
                             "maxquant_row_ids")

if( ! "fdr.mod" %in% colnames( phospho_tbl_orig )) {
  list_of_phospho_columns <- list_of_phospho_columns[list_of_phospho_columns != "fdr.mod"]
}


phospho_tbl <- phospho_tbl_orig %>%
  dplyr::select( one_of(list_of_phospho_columns)) %>%
  dplyr::mutate( phos_row_ids = sites_id )

phospho_cln <- phospho_tbl %>%
  separate_rows(uniprot_acc, position, residue, sequence, sep=":") %>%
  mutate( phospho_row_rank = row_number())


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo ("Read proteomics abundance data.")
captured_output<-capture.output(
  proteins_tbl_orig <-  vroom::vroom( args$proteins_file, delim="\t")
  ,type = "message"
)
logdebug(captured_output)


list_of_prot_columns <- c("uniprot_acc",
                          "q.mod",
                          "fdr.mod",
                          "p.mod",
                          "log2FC",
                          "comparison",
                          "maxquant_row_id")

if( ! "fdr.mod" %in% colnames( proteins_tbl_orig )) {
  list_of_prot_columns <- list_of_prot_columns[list_of_prot_columns != "fdr.mod"]
}


proteins_tbl <- proteins_tbl_orig %>%
  dplyr::select( one_of(list_of_prot_columns) ) %>%
  dplyr::rename( prot_maxquant_row_ids = "maxquant_row_id")

proteins_cln <- proteins_tbl  %>%
  separate_rows(uniprot_acc,   sep=":")  %>%
  mutate( protein_row_rank = row_number())

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

proteins_uniprot_list <- proteins_cln %>% distinct(uniprot_acc) %>% pull(uniprot_acc)
phospho_uniprot_list <- phospho_cln %>% distinct(uniprot_acc) %>% pull(uniprot_acc)
prot_phos_uniprot_list <- intersect( proteins_uniprot_list,
           phospho_uniprot_list )

logdebug(proteins_uniprot_list %>% length())
logdebug(phospho_uniprot_list %>% length())
logdebug(prot_phos_uniprot_list %>% length())

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo( "Normalisation of the phosphopeptide abundance with the protein abundance.")

# A phosphopeptide can be matched to multiple proteins and each of those protein may have their own abundance level
# Here I find distinct protein row ID and phosphosites row ID pairs
join_protein_phosopho_keys <- proteins_cln %>%
  inner_join( phospho_cln,
              by=c("uniprot_acc" = "uniprot_acc",
                   "comparison" = "comparison"),
              suffix=c(".prot", ".phos") ) %>%
  distinct ( comparison, prot_maxquant_row_ids, phos_row_ids )

# Function to find intersection between the uniprot accession of the phosphopeptide and the uniprot accession of the protein used to do fold-change normalization
intersectTwoUniprotList <-function(x, y){  paste( intersect( str_split(x, ":")[[1]], str_split(y, ":")[[1]] ), collapse=":")    }

## Which array position in the list of proteins were the selected proteins
getArrayPositionSelected <- function(x, y){   which( str_split(x, ":")[[1]]  %in% str_split(y, ":")[[1]]) }

# Once I have identified which protein was used for fold-change normalization, I have pick out the position, residue, and peptide sequence from a list that belongs to that protein
# This is the function to do that array subsetting
subsetPhosphositeDetails <- function(x,y) {    paste( str_split(x, ":")[[1]][ y], collapse=":") }

# intersectTwoUniprotList( "A:B:C", "C:D:E")

# I then have to find the unique protein row ID and phosphosites row ID pairs among the orignal
# protein groups table and phosphopeptide table, where there could be multiple possible host-proteins per phosphopeptide
basic_data_shared <- join_protein_phosopho_keys %>%
  left_join( phospho_tbl, by = c("phos_row_ids", "comparison")) %>%
  left_join( proteins_tbl, by = c("prot_maxquant_row_ids", "comparison" ),  suffix=c(".phos", ".prot") ) %>%
  mutate( uniprot_acc = purrr::map2_chr( uniprot_acc.phos, uniprot_acc.prot, intersectTwoUniprotList  )) %>%
  mutate( array_pos   = purrr::map2( uniprot_acc.phos, uniprot_acc,  getArrayPositionSelected ) ) %>%
  mutate( gene_name   = purrr::map2_chr(gene_name, array_pos, subsetPhosphositeDetails)) %>%
  mutate( position    = purrr::map2_chr(position, array_pos, subsetPhosphositeDetails)) %>%
  mutate( residue     = purrr::map2_chr(residue, array_pos, subsetPhosphositeDetails)) %>%
  mutate( residue     = purrr::map2_chr(residue, sequence, subsetPhosphositeDetails)) %>%
  dplyr::select(-uniprot_acc.phos, -uniprot_acc.prot) %>%
  distinct() %>%
  ## Calculate the new log fold-change and use the Fisher's method to calculate the updated p-value
  mutate( norm_phos_logFC = log2FC.phos - log2FC.prot) %>%
  mutate( adj_qmod.prot  = ifelse(sign(log2FC.phos)  ==  sign(log2FC.prot), 1-q.mod.prot,  q.mod.prot  )) %>%
  mutate( combined_q_mod = 1-pchisq(-2*( log(q.mod.phos) + log(adj_qmod.prot)   ), 2*2 )  ) %>%
  dplyr::mutate( status  = "Phos_and_Prot")

list_of_data_shared_columns <- c( "comparison",  "norm_phos_logFC", "combined_q_mod", "combined_fdr_mod", "sites_id",
                                  "uniprot_acc", "position", "residue", "sequence",
                                  "log2FC.phos", "q.mod.phos",
                                  "log2FC.prot", "q.mod.prot",
                                  "status", "phos_row_ids", "prot_maxquant_row_ids" )

if( "fdr.mod" %in% colnames( proteins_cln ) &
    "fdr.mod" %in% colnames( phospho_cln ) ) {
  basic_data_shared <- basic_data_shared %>%
    mutate( adj_fdrmod.prot  = ifelse(sign(log2FC.phos)  ==  sign(log2FC.prot), 1-fdr.mod.prot,  fdr.mod.prot  )) %>%
    mutate( combined_fdr_mod = 1-pchisq(-2*( log(fdr.mod.phos) + log(adj_fdrmod.prot)   ), 2*2 )  )  %>%
    dplyr::select(-adj_fdrmod.prot)

  list_of_data_shared_columns <- c( "comparison",  "norm_phos_logFC", "combined_q_mod", "combined_fdr_mod", "sites_id",
                                    "uniprot_acc", "position", "residue", "sequence",
                                    "log2FC.phos", "q.mod.phos", "fdr.mod.phos",
                                    "log2FC.prot", "q.mod.prot", "fdr.mod.prot",
                                    "status", "phos_row_ids", "prot_maxquant_row_ids" )
}

if( ! "combined_fdr_mod" %in% colnames( basic_data_shared )) {

  list_of_data_shared_columns <- list_of_data_shared_columns[ !list_of_data_shared_columns %in% c( "combined_fdr_mod", "fdr.mod.phos", "fdr.mod.prot")]
}

basic_data_shared <- basic_data_shared %>%
  dplyr::select(one_of(list_of_data_shared_columns))   %>%
  arrange(comparison, combined_q_mod, norm_phos_logFC ) %>%
  distinct() %>%
  as.data.frame


# Find all the phosphopeptides which did not have a matching protein for normalization
# Then carry over the log fold-change and p-value without changing it
basic_data_phospho_only_helper <- phospho_cln %>%
  anti_join( proteins_cln, by=c( "uniprot_acc" = "uniprot_acc" ,
                                  "comparison" = "comparison") )  %>%
  anti_join( basic_data_shared,
             by=c("sites_id" = "sites_id",
                  "comparison" = "comparison"))

basic_data_phospho_only <- phospho_tbl %>%
  inner_join( basic_data_phospho_only_helper %>%
                dplyr::select( phos_row_ids, comparison),
              by =c( "phos_row_ids", "comparison" ) ) %>%
  dplyr::mutate(norm_phos_logFC =  log2FC, combined_q_mod = q.mod) %>%
  dplyr::rename( log2FC.phos= log2FC, q.mod.phos = q.mod ) %>%
  dplyr::mutate( status  = "Phos_Only")

if( "fdr.mod" %in% colnames( proteins_cln ) &
    "fdr.mod" %in% colnames( phospho_cln ) ) {
  basic_data_phospho_only <- basic_data_phospho_only %>%
    dplyr::rename(  fdr.mod.phos = fdr.mod )
}

list_of_phospho_only_columns <- c( "comparison", "norm_phos_logFC", "combined_q_mod", "sites_id", "uniprot_acc", "position", "residue", "sequence",
"log2FC.phos", "q.mod.phos", "fdr.mod.phos", "status", "phos_row_ids" )

if( ! "combined_fdr_mod" %in% colnames( basic_data_shared )) {

  list_of_phospho_only_columns <- list_of_phospho_only_columns[ !list_of_phospho_only_columns %in% c( "combined_fdr_mod", "fdr.mod.phos")]
}

basic_data_phospho_only <- basic_data_phospho_only %>%
  dplyr::select(one_of( list_of_phospho_only_columns ))    %>%
  arrange(comparison, combined_q_mod, norm_phos_logFC ) %>%
  distinct() %>%
  as.data.frame

basic_data <- basic_data_shared  %>%
  bind_rows(basic_data_phospho_only %>%
              dplyr::mutate( position = purrr::map_chr( position, as.character)) )

 if (nrow(basic_data %>% distinct(sites_id)) != nrow(phospho_cln %>% distinct(sites_id))){
   logwarn("nrow(basic_data) != nrow(phospho_cln)")
 }

vroom::vroom_write( basic_data, file.path( args$output_dir,  "norm_phosphosite_lfc_minus_protein_lfc_basic.tsv"))

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Join normalized phosphopeptide abundance table with phosphosite annotations")
annotation_from_phospho_tbl <- phospho_tbl_orig %>%
  dplyr::mutate( phos_row_ids = sites_id ) %>%
  dplyr::select(-q.mod, -p.mod, -log2FC, -uniprot_acc,  -position, -residue, -sequence)

annotated_phos_tbl <- basic_data %>%
  left_join( annotation_from_phospho_tbl, by=c("sites_id" = "sites_id",
                                               "comparison" = "comparison",
                                               "phos_row_ids" = "phos_row_ids")) %>%
  dplyr::select( !matches( "log2norm\\.\\d+\\.(left|right)") &
                 !matches("raw\\.\\d+\\.(left|right)")) %>%
  arrange(comparison, combined_q_mod, norm_phos_logFC ) %>%
  distinct()

vroom::vroom_write( annotated_phos_tbl,
                    file.path( args$output_dir,
                               "norm_phosphosite_lfc_minus_protein_lfc_annotated.tsv"))

list_of_long_columns <- intersect( colnames(annotated_phos_tbl),
                                   c("protein_names",
                                     "ENSEMBL",
                                     "PROTEIN-NAMES",
                                     "KEYWORDS",
                                     "GO-ID",
                                     "go_biological_process",
                                     "go_cellular_compartment",
                                     "go_molecular_function",
                                     "reactome_term",
                                     "majority_protein_ids") )


writexl::write_xlsx(annotated_phos_tbl %>%
                      mutate_at( list_of_long_columns, ~substr(., 1, 32760) ),
                    file.path (args$output_dir,
                               "norm_phosphosite_lfc_minus_protein_lfc_annotated.xlsx")  )

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Compare before and after normalization with protein abundance")
before_prot_norm <- phospho_cln %>%
  dplyr::filter( q.mod < 0.05) %>%
  dplyr::distinct(comparison, sites_id) %>%
  dplyr::mutate(  Normalization =  "Before")


after_prot_norm <- basic_data %>%
  dplyr::filter( combined_q_mod < 0.05) %>%
  dplyr::distinct(comparison, sites_id) %>%
  dplyr::mutate(  Normalization =  "After")

comparisons_order <- before_prot_norm %>%
  distinct( comparison) %>%
  pull(comparison)

compare_before_and_after <- before_prot_norm %>%
  full_join( after_prot_norm, by=c("comparison", "sites_id"), suffix=c(".before", ".after")) %>%
  mutate( Normalization = case_when (  is.na(Normalization.before) & !is.na(Normalization.after) ~ "After Only",
          !is.na(Normalization.before) & !is.na(Normalization.after) ~ "Before & After",
          !is.na(Normalization.before) & is.na(Normalization.after) ~ "Before Only",
          TRUE ~ NA_character_ ) ) %>%
  mutate(Normalization = factor( Normalization,
                                 levels=c("Before Only",
                                          "Before & After",
                                          "After Only")) ) %>%
  mutate( comparison = factor(comparison, levels = comparisons_order)) %>%
  group_by( comparison,  Normalization) %>%
  summarise( Counts = n()) %>%
  ungroup()

 vroom::vroom_write( compare_before_and_after,
                     file.path(args$output_dir, "compare_before_and_after_norm_by_prot_abundance.tsv"))


 cmp_before_after_plot <- compare_before_and_after %>%
   ggplot( aes( Normalization, Counts)) +
   geom_col() +
   facet_grid( . ~ comparison ) +
   geom_text(stat='identity', aes(label= Counts), vjust=-0.5) +
   theme(axis.text.x = element_text(angle = 90))



for( format_ext in args$plots_format) {
  file_name<-file.path(args$output_dir,paste0("compare_before_and_after_norm_by_prot_abundance.",format_ext))
  captured_output<capture.output(
   ggsave( plot=cmp_before_after_plot, file_name, width = 14, height=10)
    , type = "message"
  )
  logdebug(captured_output)
}

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Count num. differentiall abundant phosphosites after normalization by protein abundance.")

output_counts_tbl <- annotated_phos_tbl %>%
   dplyr::select( comparison, norm_phos_logFC, combined_q_mod) %>%
   group_by(comparison) %>%
   nest( data = c( norm_phos_logFC, combined_q_mod) ) %>%
   ungroup %>%
   mutate( counts = purrr::map( data, function(x) {    countStatDeGenes(x, lfc_thresh = 0,
                                                                        q_val_thresh = 0.05,
                                                                        log_fc_column = norm_phos_logFC,
                                                                        q_value_column = combined_q_mod )}   )  ) %>%
  dplyr::select(-data) %>%
  unnest(counts)

 vroom::vroom_write(output_counts_tbl,
                    file.path(args$output_dir,
                              "num_sig_diff_abundant_norm_phosphosites.tab"))
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Create Volcano Plot")
qm.threshold <- 0.05
logFC.threshold <- 1

basic_data_volcano_plot_data <- basic_data %>%
  left_join( annotated_phos_tbl %>%
               dplyr::distinct( uniprot_acc, gene_name),
             by=c("uniprot_acc" = "uniprot_acc")) %>%
  dplyr::mutate( colour= case_when ( abs(norm_phos_logFC) >= logFC.threshold & combined_q_mod >= qm.threshold ~ "orange",
                                     abs(norm_phos_logFC) >= logFC.threshold & combined_q_mod < qm.threshold ~ "purple",
                                     abs(norm_phos_logFC) < logFC.threshold & combined_q_mod  < qm.threshold ~ "blue",
                                     TRUE ~ "black" ))  %>%
  dplyr::mutate( colour = factor( colour, levels=c("orange", "purple", "blue" ,"black")))

basic_data_volcano_plot <- basic_data_volcano_plot_data %>%
  ggplot( aes( norm_phos_logFC, -log10(combined_q_mod), col=colour, key=gene_name )) +
  geom_point() +
  facet_grid( . ~ comparison)  +
            scale_colour_manual(values = c("orange", "purple", "blue" ,"black"),
                                labels=c(paste0("Not significant, logFC > ",
                                  logFC.threshold),
                                  paste0("Significant, logFC >= ",
                                         logFC.threshold),
                                  paste0("Significant, logFC <",
                                         logFC.threshold),
                                  "Not Significant"))

for( format_ext in args$plots_format) {
  file_name<-file.path(args$output_dir,paste0("volplot_gg_phos_vs_prot_all.",format_ext))
  captured_output<capture.output(
  ggsave(  filename=file_name, plot=basic_data_volcano_plot, width=15, height=6   )
    , type = "message"
  )
  logdebug(captured_output)
}

captured_output<capture.output(
  ggplotly(basic_data_volcano_plot)
    , type = "message"
  )
logdebug(captured_output)

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
te<-toc(quiet = TRUE)
loginfo("%f sec elapsed",te$toc-te$tic)
writeLines(capture.output(sessionInfo()), file.path(args$output_dir,"sessionInfo.txt"))

