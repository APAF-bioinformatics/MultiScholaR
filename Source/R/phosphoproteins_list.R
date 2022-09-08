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
p_load(rlang)


p_load(qvalue)
p_load(knitr)
p_load(GO.db)
p_load(clusterProfiler)

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

parser <- add_option(parser, c("-c","--config"), type = "character", default = "config_phos.ini", dest = "config",
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

parser <- add_option(parser, "--norm_phos_logfc_file", type="character",
                     help="Results table in which the phosphorylation log fold-change is normalized by protein log fold-change.",
                     metavar="string")

parser <- add_option(parser, "--log_fc_column_name", type="character",
                     help="Column name in the input file that contains the log fold-change values of the phosphosites.",
                     metavar="string")

parser <- add_option(parser, "--fdr_column_name", type="character",
                     help="Column name in the input file that contains the false discovery rate values of the phosphosites.",
                     metavar="string")

parser <- add_option(parser, c("--proteins_file"), type="character", dest = "proteins_file",
                     help="File with table of differentially expressed proteins.",
                     metavar="string")

parser <- add_option(parser, c("--background"), type="character", dest = "background",
                     help="phosphoproteins or proteins_and_phosphoproteins.",
                     metavar="string")

parser <- add_option(parser, c("--annotation_file"), type="character", dest = "annotation_file",
                     help="File with protein accession and functional annotation data.",
                     metavar="string")

parser <- add_option(parser, c("--protein_id"), type="character", dest = "protein_id",
                     help="The protein id in the annotation file.",
                     metavar="string")

parser <- add_option(parser, c("--annotation_id"), type="character", dest = "annotation_id",
                     help="The annotatio id in the annotation file.",
                     metavar="string")

parser <- add_option(parser, c("--annotation_column"), type="character", dest = "annotation_column",
                     help="Column with the long name and biological details of the functional annotation.",
                     metavar="string")


parser <- add_option(parser, c("--annotation_type"), type="character", dest = "annotation_type",
                     help="Name of the functional annotation data, added to the output file name.",
                     metavar="string")

parser <- add_option(parser, c("--aspect_column"), type="character", dest = "aspect_column",
                     help="The aspect of the GO term.",
                     metavar="string")

parser <- add_option(parser, c("--dictionary_file"), type="character", dest = "dictionary_file",
                     help="File that converts annotation ID to annotation name.",
                     metavar="string")

parser <- add_option(parser, c("--max_gene_set_size"), type="character", dest = "max_gene_set_size",
                     help="The maximum number of genes associaed with each gene set.",
                     metavar="string")

parser <- add_option(parser, c("--min_gene_set_size"), type="character", dest = "min_gene_set_size",
                     help="The minimum number of genes associaed with each gene set.",
                     metavar="string")

parser <- add_option(parser, "--site_p_val_thresh", type = "double",
                     help = "p-value threshold below which a phosphosite is significantly enriched",
                     metavar = "double")

parser <- add_option(parser, "--p_val_thresh", type = "double",
                     help = "p-value threshold below which a GO term is significantly enriched",
                     metavar = "double")

parser <- add_option(parser, "--uniprot_to_gene_symbol_file", type = "character",
                     help = "A file that contains the dictionary to convert uniprot accession to gene symbol. Uses the column specified in 'protein_id_lookup_column' flag for protein ID. Uses the column specified in 'gene_symbol_column' for the gene symbol column.",
                     metavar = "string")

parser <- add_option(parser, "--protein_id_lookup_column", type = "character",
                     help = "The  name of the column that contained the protein ID to convert into gene symobl.",
                     metavar = "string")

parser <- add_option(parser, "--gene_symbol_column", type = "character",
                     help = "The  name of the column that contained the gene symobl.",
                     metavar = "string")

#parse comand line arguments first.
args <- parse_args(parser)

#parse and merge the configuration file options.
if (args$config != "") {
  args <- config.list.merge(eval.config(file = args$config, config = "phosphoproteins_list"), args)
}

args <- setArgsDefault(args, "output_dir", as_func=as.character, default_val="phosphoproteins_list" )

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

args <- setArgsDefault(args, "log_fc_column_name", as_func=as.character, default_val="norm_phos_logFC" )
args <- setArgsDefault(args, "fdr_column_name", as_func=as.character, default_val="combined_q_mod" )
args <- setArgsDefault(args, "p_val_thresh", as_func=as.double, default_val=0.05 )
args <- setArgsDefault(args, "site_p_val_thresh", as_func=as.double, default_val=0.05 )

args <- setArgsDefault(args, "max_gene_set_size", as_func=as.character, default_val="200" )
args <- setArgsDefault(args, "min_gene_set_size", as_func=as.character, default_val="4" )
args <- setArgsDefault(args, "protein_id", as_func=as.character, default_val="uniprot_acc" )
args <- setArgsDefault(args, "annotation_id", as_func=as.character, default_val="go_id" )
args <- setArgsDefault(args, "aspect_column", as_func=as.character, default_val=NULL )
args <- setArgsDefault(args, "background", as_func=as.character, default_val="proteins_and_phosphoproteins" )



if( isArgumentDefined(args, "aspect_column" )) {
  if( args$aspect_column == "NULL" ) {
    args$aspect_column <- NULL
  }
}
args <- setArgsDefault(args, "annotation_column", as_func=as.character, default_val="term" )

args <- setArgsDefault(args, "annotation_type", as_func=as.character, default_val="annotation" )

testRequiredArguments(args, c(
  "proteins_file",
  "norm_phos_logfc_file"
))

if( isArgumentDefined( args, "annotation_file" )) {

  testRequiredArguments(args, c(
    "annotation_file",
    "annotation_type"
  ))

  testRequiredFiles(c(
    args$annotation_file))

  if (! isArgumentDefined(args, "aspect_column")) {

    testRequiredArguments(args, c(
      "dictionary_file",
      "annotation_column"
    ))

    testRequiredFiles(c(
      args$dictionary_file))

  }

}

testRequiredFiles(c(
  args$proteins_file,
  args$norm_phos_logfc_file))

args<-parseString(args,c("plots_format",
                         "log_fc_column_name",
                         "fdr_column_name",
                         "max_gene_set_size",
                         "min_gene_set_size",
                         "protein_id",
                         "annotation_id",
                         "aspect_column",
                         "annotation_column",
                         "annotation_type"))

args<-parseType(args,
                c("p_val_thresh",
                  "site_p_val_thresh"),
                as.double)

if(isArgumentDefined(args,"plots_format")) {
  args <- parseList(args,c("plots_format"))
} else {
  logwarn("plots_format is undefined, default output set to pdf.")
  args$plots_format <- list("pdf")
}

if(isArgumentDefined(args, "uniprot_to_gene_symbol_file")) {
  testRequiredFiles(c(
    args$uniprot_to_gene_symbol_file))

  args <- setArgsDefault(args, "protein_id_lookup_column", as_func=as.character, default_val="Entry" )
  args <- setArgsDefault(args, "gene_symbol_column", as_func=as.character, default_val="Gene names" )

}

## ---------------------------------------------------------------------------------------------------------------------------------------------------

loginfo("Read phosphosite log fold-change normalized by protein log fold-change file")
captured_output<-capture.output(
  de_phos <- vroom::vroom( args$norm_phos_logfc_file )
  ,type = "message"
)
logdebug(captured_output)

## ---------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Check log fold-change and FDR column names exists in input file.")


if ( ! args$log_fc_column_name %in% colnames(de_phos ) ) {
  logerror("Column '%s' is not found in the input table.",args$log_fc_column_name )
}

if ( ! args$fdr_column_name %in% colnames(de_phos ) ) {
  logerror("Column '%s' is not found in the input table.",args$fdr_column_name)
}

if(  ! args$log_fc_column_name %in% colnames(de_phos )  |
     ! args$fdr_column_name %in% colnames(de_phos )  ) {
  stop()
}
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

positive_phosphoproteins <- de_phos %>%
  dplyr::filter( !!rlang::sym(args$fdr_column_name) < args$site_p_val_thresh & !!rlang::sym(args$log_fc_column_name) > 0 ) %>%
  mutate( uniprot_acc_first = purrr::map_chr( uniprot_acc, ~str_split(., ":") %>% map_chr(1)))  %>%
  mutate( uniprot_acc_first = str_replace_all( uniprot_acc_first, "-\\d+$", ""))  %>% # Strip away isoform information
  mutate( gene_name_first = purrr::map_chr( gene_name, ~str_split(., ":") %>% map_chr(1)))  %>%
  mutate( protein_name_first = purrr::map_chr( `PROTEIN-NAMES`, ~str_split(., ":") %>% map_chr(1)))  %>%
  group_by(comparison, uniprot_acc_first, gene_name_first, protein_name_first) %>%
  summarise( max_norm_phos_logFC = max(!!rlang::sym(args$log_fc_column_name))) %>%
  ungroup() %>%
  arrange( comparison, desc(max_norm_phos_logFC  ) )

vroom::vroom_write( positive_phosphoproteins,
                    file.path( args$output_dir,
                               "all_phosphoproteins_with_positive_logFC_sites.tab" ))

 list_of_comparisons <- positive_phosphoproteins %>% distinct( comparison) %>% pull( comparison)

 purrr::walk( list_of_comparisons, ~createDirIfNotExists( file.path(args$output_dir,. )) )

 purrr::walk( list_of_comparisons, function( input_comparison){

   positive_phosphoproteins %>%
     dplyr::filter( comparison == input_comparison) %>%
     dplyr::select( uniprot_acc_first) %>%
     vroom::vroom_write( file.path( args$output_dir,
                                    input_comparison,
                                    "all_phosphoproteins_with_positive_logFC_sites.tab" ),
                         col_names=FALSE)

 } )

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# norm_phos_logFC

negative_phosphoproteins <- de_phos %>%
  dplyr::filter( !!rlang::sym(args$fdr_column_name) < args$site_p_val_thresh & !!rlang::sym(args$log_fc_column_name)  < 0 ) %>%
  mutate( uniprot_acc_first = purrr::map_chr( uniprot_acc, ~str_split(., ":") %>% map_chr(1)))  %>%
  mutate( uniprot_acc_first = str_replace_all( uniprot_acc_first, "-\\d+$", ""))  %>% # Strip away isoform information
  mutate( gene_name_first = purrr::map_chr( gene_name, ~str_split(., ":") %>% map_chr(1)))  %>%
  mutate( protein_name_first = purrr::map_chr( `PROTEIN-NAMES`, ~str_split(., ":") %>% map_chr(1)))  %>%
  group_by(comparison, uniprot_acc_first, gene_name_first, protein_name_first) %>%
  summarise( min_norm_phos_logFC = min(!!rlang::sym(args$log_fc_column_name))) %>%
  ungroup() %>%
  arrange( comparison, min_norm_phos_logFC  )

vroom::vroom_write( negative_phosphoproteins,
                    file.path( args$output_dir,
                               "all_phosphoproteins_with_negative_logFC_sites.tab" ))

list_of_comparisons <- negative_phosphoproteins %>% distinct( comparison) %>% pull( comparison)

purrr::walk( list_of_comparisons, ~createDirIfNotExists( file.path(args$output_dir,. )) )

purrr::walk( list_of_comparisons, function( input_comparison){

  negative_phosphoproteins %>%
    dplyr::filter( comparison == input_comparison) %>%
    dplyr::select( uniprot_acc_first) %>%
    vroom::vroom_write( file.path( args$output_dir,
                                   input_comparison,
                                   "all_phosphoproteins_with_negative_logFC_sites.tab" ),
                        col_names=FALSE)

} )

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

positive_only_phosphoproteins <- positive_phosphoproteins %>%
  anti_join( negative_phosphoproteins, by=c("uniprot_acc_first" = "uniprot_acc_first"))

vroom::vroom_write( positive_only_phosphoproteins,
                    file.path(args$output_dir,
                              "phosphoproteins_with_only_positive_logFC_sites.tab" ))


list_of_comparisons <- positive_only_phosphoproteins %>% distinct( comparison) %>% pull( comparison)

purrr::walk( list_of_comparisons, ~createDirIfNotExists( file.path(args$output_dir,. )) )

purrr::walk( list_of_comparisons, function( input_comparison){

  positive_only_phosphoproteins %>%
    dplyr::filter( comparison == input_comparison) %>%
    dplyr::select( uniprot_acc_first) %>%
    vroom::vroom_write( file.path( args$output_dir,
                                   input_comparison,
                                   "phosphoproteins_with_only_positive_logFC_sites.tab" ),
                        col_names=FALSE)

} )

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

negative_only_phosphoproteins <- negative_phosphoproteins  %>%
  anti_join( positive_phosphoproteins, by=c("uniprot_acc_first" = "uniprot_acc_first"))

vroom::vroom_write( negative_only_phosphoproteins,
                    file.path(args$output_dir,
                              "phosphoproteins_with_only_negative_logFC_sites.tab" ) )

list_of_comparisons <- negative_only_phosphoproteins %>% distinct( comparison) %>% pull( comparison)

purrr::walk( list_of_comparisons, ~createDirIfNotExists( file.path(args$output_dir,. )) )

purrr::walk( list_of_comparisons, function( input_comparison){

  negative_only_phosphoproteins %>%
    dplyr::filter( comparison == input_comparison) %>%
    dplyr::select( uniprot_acc_first) %>%
    vroom::vroom_write( file.path( args$output_dir,
                                   input_comparison,
                                   "phosphoproteins_with_only_negative_logFC_sites.tab" ),
                        col_names=FALSE)

} )



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

overlapping_phosphoproteins <- positive_phosphoproteins %>%
  inner_join( negative_phosphoproteins, by=c("uniprot_acc_first" = "uniprot_acc_first",
                                             "comparison" = "comparison",
                                             "gene_name_first" = "gene_name_first",
                                             "protein_name_first" = "protein_name_first"))

vroom::vroom_write( overlapping_phosphoproteins,
                    file.path(args$output_dir,
                              "phosphoproteins_with_positive_and_negative_logFC_sites.tab" ) )


list_of_comparisons <- overlapping_phosphoproteins %>% distinct( comparison) %>% pull( comparison)

purrr::walk( list_of_comparisons, ~createDirIfNotExists( file.path(args$output_dir,. )) )

purrr::walk( list_of_comparisons, function( input_comparison){

  overlapping_phosphoproteins %>%
    dplyr::filter( comparison == input_comparison) %>%
    dplyr::select( uniprot_acc_first) %>%
    vroom::vroom_write( file.path( args$output_dir,
                                   input_comparison,
                                   "phosphoproteins_with_positive_and_negative_logFC_sites.tab" ),
                        col_names=FALSE)

} )




## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

all_phosphoproteins_with_significant_da_sites <- positive_phosphoproteins %>%
  bind_rows( negative_phosphoproteins) %>%
  distinct()

vroom::vroom_write( all_phosphoproteins_with_significant_da_sites,
                    file.path(args$output_dir,
                              "all_phosphoproteins_with_significant_da_sites.tab" ))

list_of_comparisons <- all_phosphoproteins_with_significant_da_sites %>% distinct( comparison) %>% pull( comparison)

purrr::walk( list_of_comparisons, ~createDirIfNotExists( file.path(args$output_dir,. )) )

purrr::walk( list_of_comparisons, function( input_comparison){

  all_phosphoproteins_with_significant_da_sites %>%
    dplyr::filter( comparison == input_comparison) %>%
    dplyr::select( uniprot_acc_first) %>%
    vroom::vroom_write( file.path( args$output_dir,
                                   input_comparison,
                                   "all_phosphoproteins_with_significant_da_sites.tab" ),
                        col_names=FALSE)

} )



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Classify phosphoproteins as up and down by total log2 FC of significant phosphosites

group_phosphoproteins_by_phosphosites_lfc_total <- de_phos  %>%
  mutate( uniprot_acc_first = purrr::map_chr( uniprot_acc, ~str_split(., ":") %>% map_chr(1)))  %>%
  mutate( uniprot_acc_first = str_replace_all( uniprot_acc_first, "-\\d+$", ""))  %>%
  dplyr::filter(  !!rlang::sym(args$fdr_column_name) < args$site_p_val_thresh ) %>%
  group_by( comparison, uniprot_acc_first ) %>%
  summarise( total_log2FC = sum( !!rlang::sym(args$log_fc_column_name)) )  %>%
  ungroup() %>%
  mutate( direction = case_when( total_log2FC > 0 ~ "up",
                                 total_log2FC < 0 ~ "down",
                                 TRUE ~ "no_change")) %>%
  dplyr::filter( direction != "no_change")

gene_names_list <- de_phos %>%
  mutate( uniprot_acc_first = purrr::map_chr( uniprot_acc, ~str_split(., ":") %>% map_chr(1)))  %>%
  mutate( uniprot_acc_first = str_replace_all( uniprot_acc_first, "-\\d+$", ""))  %>% # Strip away isoform information
  mutate( gene_name_first = purrr::map_chr( gene_name, ~str_split(., ":") %>% map_chr(1)))  %>%
  mutate( protein_name_first = purrr::map_chr( `PROTEIN-NAMES`, ~str_split(., ":") %>% map_chr(1)))  %>%
  distinct(comparison, uniprot_acc_first, gene_name_first, protein_name_first)


group_phosphoproteins_by_phosphosites_lfc_total_up <- group_phosphoproteins_by_phosphosites_lfc_total %>%
  dplyr::filter( direction == "up") %>%
  left_join ( gene_names_list, by =c( "uniprot_acc_first" = "uniprot_acc_first",
                                      "comparison" = "comparison"))


vroom::vroom_write( group_phosphoproteins_by_phosphosites_lfc_total_up,
                    file.path(args$output_dir,
                              "group_phosphoproteins_by_phosphosites_lfc_total_up.tab" ))

list_of_comparisons <- group_phosphoproteins_by_phosphosites_lfc_total_up %>% distinct( comparison) %>% pull( comparison)

purrr::walk( list_of_comparisons, ~createDirIfNotExists( file.path(args$output_dir,. )) )

purrr::walk( list_of_comparisons, function( input_comparison){

  group_phosphoproteins_by_phosphosites_lfc_total_up %>%
    dplyr::filter( comparison == input_comparison) %>%
    dplyr::select( uniprot_acc_first) %>%
    vroom::vroom_write( file.path( args$output_dir,
                                   input_comparison,
                                   "group_phosphoproteins_by_phosphosites_lfc_total_up.tab" ),
                        col_names=FALSE)

} )


group_phosphoproteins_by_phosphosites_lfc_total_down <- group_phosphoproteins_by_phosphosites_lfc_total %>%
  dplyr::filter( direction == "down") %>%
  left_join ( gene_names_list, by =c( "uniprot_acc_first" = "uniprot_acc",
                                      "comparison" = "comparison"))



vroom::vroom_write( group_phosphoproteins_by_phosphosites_lfc_total_down,
                    file.path(args$output_dir,
                              "group_phosphoproteins_by_phosphosites_lfc_total_down.tab" ))


list_of_comparisons <- group_phosphoproteins_by_phosphosites_lfc_total_down %>% distinct( comparison) %>% pull( comparison)

purrr::walk( list_of_comparisons, ~createDirIfNotExists( file.path(args$output_dir,. )) )

purrr::walk( list_of_comparisons, function( input_comparison){

  group_phosphoproteins_by_phosphosites_lfc_total_down %>%
    dplyr::filter( comparison == input_comparison) %>%
    dplyr::select( uniprot_acc_first) %>%
    vroom::vroom_write( file.path( args$output_dir,
                                   input_comparison,
                                   "group_phosphoproteins_by_phosphosites_lfc_total_down.tab" ),
                        col_names=FALSE)

} )




## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo ("Read proteomics abundance data.")
captured_output<-capture.output(
  proteins_tbl_orig <-  vroom::vroom( args$proteins_file, delim="\t")
  ,type = "message"
)
logdebug(captured_output)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

background_phosphoproteins  <- de_phos %>%
  mutate( uniprot_acc_first = purrr::map_chr( uniprot_acc, ~str_split(., ":") %>% map_chr(1)))  %>%
  mutate( uniprot_acc_first = str_replace_all( uniprot_acc_first, "-\\d+$", ""))  %>% # Strip away isoform information
  mutate( gene_name_first = purrr::map_chr( gene_name, ~str_split(., ":") %>% map_chr(1)))  %>%
  mutate( protein_name_first = purrr::map_chr( `PROTEIN-NAMES`, ~str_split(., ":") %>% map_chr(1)))  %>%
  distinct( uniprot_acc_first) %>%
  arrange( uniprot_acc_first )

vroom::vroom_write( background_phosphoproteins,
                    file.path(args$output_dir, "background_phosphoproteins.tab" ),
                    col_names=FALSE )

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

background_proteins <- proteins_tbl_orig %>%
  distinct(uniprot_acc) %>%
  mutate( uniprot_acc_first = purrr::map_chr( uniprot_acc, ~str_split(., ":") %>% map_chr(1)))  %>%
  mutate( uniprot_acc_first = str_replace_all( uniprot_acc_first, "-\\d+$", ""))  %>%
  distinct(uniprot_acc)

vroom::vroom_write( background_proteins,
                    file.path(args$output_dir, "background_proteins.tab" ),
                    col_names=FALSE )

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

background_proteins_phosphoproteins <- background_proteins %>%
  bind_rows( background_phosphoproteins) %>%
  distinct( uniprot_acc_first)

vroom::vroom_write( background_proteins_phosphoproteins,
                    file.path(args$output_dir, "background_proteins_phosphoproteins.tab" ),
                    col_names=FALSE )

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


if (  !is.null( args$annotation_file )) {


  ## Tidy up the annotation ID to annotation term name dictionary

  dictionary <- vroom::vroom( args$dictionary_file )

  id_to_annotation_dictionary <- buildAnnotationIdToAnnotationNameDictionary( input_table=dictionary,
                                                                annotation_column = !!rlang::sym(args$annotation_column),
                                                                annotation_id_column = !!rlang::sym(args$annotation_id))

  ## Create gene sets list




  ## preparing the enrichment test
  go_annot <- vroom::vroom(  args$annotation_file   )

  background_list <- background_phosphoproteins
  if( args$background == "phosphoproteins") {
    background_list <- background_phosphoproteins
  } else if ( args$background == "proteins_and_phosphoproteins") {
    background_list <- background_proteins_phosphoproteins
  }

  #print( args$min_gene_set_size)
  min_gene_set_size_list <- parseNumList(args$min_gene_set_size)
  max_gene_set_size_list <- parseNumList(args$max_gene_set_size)

 list_of_comparisons <- all_phosphoproteins_with_significant_da_sites %>%
   distinct(comparison) %>%
   pull(comparison)

 # Tidy up GO aspect list, marked as null if not using GO terms
 go_aspect_list <- NA
 if(!is.null(args$aspect_column)) {
   go_aspect_list <- go_annot %>%
     dplyr::filter( !is.na(!!rlang::sym( args$aspect_column) )) %>%
     distinct( !!rlang::sym( args$aspect_column) ) %>%
     pull( !!rlang::sym( args$aspect_column) ) # c("C", "F", "P")
 } else {
   go_aspect_list <- NA
 }

 # print( paste("is.na(go_aspect_list) =", is.na(go_aspect_list)) )

 list_of_genes_list <- list( all_significant=all_phosphoproteins_with_significant_da_sites,
                             overlap_only=overlapping_phosphoproteins,
                             negative_only=negative_only_phosphoproteins,
                             positive_only=positive_only_phosphoproteins,
                             negative_plus_overlap=negative_phosphoproteins,
                             positive_plus_overlap=positive_phosphoproteins,
                             negative_sum_sig_phosphosites=group_phosphoproteins_by_phosphosites_lfc_total_down,
                             positive_sum_sig_phosphosites=group_phosphoproteins_by_phosphosites_lfc_total_up)

 input_params <- cross( list(
   names_of_genes_list = names( list_of_genes_list),
   go_aspect=go_aspect_list,
   input_comparison = list_of_comparisons,
   min_size = min_gene_set_size_list,
   max_size = max_gene_set_size_list) )

 input_params_updated <- purrr::map( input_params,
                                     function(x) {
                                       x$input_table <- list_of_genes_list[[x$names_of_genes_list]]
                                     return(x)})

 runOneGoEnrichmentInOutFunctionPartial <- purrr::partial ( runOneGoEnrichmentInOutFunction,
                  comparison_column = comparison,
                  protein_id_column = uniprot_acc_first,
                  go_annot = go_annot,
                  background_list = background_list,
                  id_to_annotation_dictionary=id_to_annotation_dictionary,
                  annotation_id=!!rlang::sym(args$annotation_id),
                  protein_id=!!rlang::sym(args$protein_id),
                  aspect_column=args$aspect_column,
                  p_val_thresh=args$p_val_thresh)

 enrichment_result <- purrr::map( input_params_updated,
                                  ~runOneGoEnrichmentInOutFunctionPartial(
                                    names_of_genes_list = .$names_of_genes_list,
                                    input_table=.$input_table,
                                    go_aspect=.$go_aspect,
                                    input_comparison=.$input_comparison,
                                    min_gene_set_size=.$min_size,
                                    max_gene_set_size=.$max_size)) %>%
   bind_rows()

 if(is.null(enrichment_result) |
    nrow(enrichment_result) == 0 ) {
   warnings("No enriched terms were identified.")
 } else {


   enrichment_result_add_gene_symbol <- NA
   ## Convert Uniprot accession to gene names
   if(isArgumentDefined(args, "uniprot_to_gene_symbol_file")) {
     # args$uniprot_to_gene_symbol_file <- "/home/ubuntu/Workings/2021/ALPK1_BMP_06/Data/UniProt/data.tab"

     # args$protein_id_lookup_column <- "Entry"
     # args$gene_symbol_column <- "Gene names"

     # Clean up protein ID to gene symbol table
     uniprot_tab_delimited_tbl <- vroom::vroom( file.path( args$uniprot_to_gene_symbol_file))

     uniprot_to_gene_symbol_dict <- getUniprotAccToGeneSymbolDictionary( uniprot_tab_delimited_tbl,
                                                                         !!rlang::sym(args$protein_id_lookup_column),
                                                                         !!rlang::sym(args$gene_symbol_column),
                                                                         !!rlang::sym(args$protein_id) )

     convertProteinAccToGeneSymbolPartial <- purrr::partial( convertProteinAccToGeneSymbol,
                                                             dictionary = uniprot_to_gene_symbol_dict)

     print( head(enrichment_result))

     enrichment_result_add_gene_symbol <- enrichment_result %>%
       mutate( gene_id_list = str_split( geneID, "/") ) %>%
       mutate( gene_symbol = purrr::map_chr( gene_id_list,
                                             convertProteinAccToGeneSymbolPartial)) %>%
       dplyr::select(-gene_id_list)

   } else {
     enrichment_result_add_gene_symbol <- enrichment_result
   }



   ## generate the output files
   purrr::walk(list_of_comparisons, function(input_comparison) {
     output_file <- paste0( args$annotation_type, "_table_", input_comparison, ".tab" )

     vroom::vroom_write(enrichment_result_add_gene_symbol %>%
                          dplyr::filter( comparison == input_comparison),
                        file=file.path( args$output_dir,
                                        input_comparison,
                                        output_file ) ) })

 }

}



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
te<-toc(quiet = TRUE)
loginfo("%f sec elapsed",te$toc-te$tic)
writeLines(capture.output(sessionInfo()), file.path(args$output_dir,"sessionInfo.txt"))

