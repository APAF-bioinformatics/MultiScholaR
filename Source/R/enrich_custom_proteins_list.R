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

parser <- add_option(parser, c("-c","--config"), type = "character", default = "config_prot.ini", dest = "config",
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

parser <- add_option(parser, c("--query_proteins_file"), type="character", dest = "query_proteins_file",
                     help="File with table of differentially expressed proteins.",
                     metavar="string")

parser <- add_option(parser, "--qp_protein_id_column", type="character",
                     help="Column name in the input file that contains the protein ID.",
                     metavar="string")

parser <- add_option(parser, "--qp_protein_sets_column", type="character",
                     help="Column name in the input file that contains the protein sets each protein is part of.",
                     metavar="string")

parser <- add_option(parser, c("--background_file"), type="character", dest = "background_file",
                     help="List of proteins in the background set",
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
  args <- config.list.merge(eval.config(file = args$config, config = "custom_proteins_list"), args)
}

args <- setArgsDefault(args, "output_dir", as_func=as.character, default_val="custom_proteins_list" )

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

args <- setArgsDefault(args, "qp_protein_id_column", as_func=as.character, default_val="uniprot_acc" )
args <- setArgsDefault(args, "qp_protein_sets_column", as_func=as.character, default_val="cluster" )
args <- setArgsDefault(args, "p_val_thresh", as_func=as.double, default_val=0.05 )

args <- setArgsDefault(args, "max_gene_set_size", as_func=as.character, default_val="200" )
args <- setArgsDefault(args, "min_gene_set_size", as_func=as.character, default_val="4" )
args <- setArgsDefault(args, "protein_id", as_func=as.character, default_val="uniprot_acc" )
args <- setArgsDefault(args, "annotation_id", as_func=as.character, default_val="go_id" )
args <- setArgsDefault(args, "aspect_column", as_func=as.character, default_val=NULL )

if( isArgumentDefined(args, "aspect_column" )) {
  if( args$aspect_column == "NULL" ) {
    args$aspect_column <- NULL
  }
}
args <- setArgsDefault(args, "annotation_column", as_func=as.character, default_val="term" )
args <- setArgsDefault(args, "annotation_type", as_func=as.character, default_val="annotation" )

testRequiredArguments(args, c(
  "query_proteins_file",
  "background_file"
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

if(isArgumentDefined(args, "uniprot_to_gene_symbol_file")) {
  testRequiredFiles(c(
    args$uniprot_to_gene_symbol_file))

  args <- setArgsDefault(args, "protein_id_lookup_column", as_func=as.character, default_val="Entry" )
  args <- setArgsDefault(args, "gene_symbol_column", as_func=as.character, default_val="Gene names" )

}


testRequiredFiles(c(
  args$query_proteins_file,
  args$background_file))

args<-parseString(args,c("plots_format",
                         "qp_protein_id_column",
                         "qp_protein_sets_column",
                         "max_gene_set_size",
                         "min_gene_set_size",
                         "protein_id",
                         "annotation_id",
                         "aspect_column",
                         "annotation_column",
                         "annotation_type"))

args<-parseType(args,
                c("p_val_thresh"),
                as.double)

if(isArgumentDefined(args,"plots_format")) {
  args <- parseList(args,c("plots_format"))
} else {
  logwarn("plots_format is undefined, default output set to pdf.")
  args$plots_format <- list("pdf")
}



## ---------------------------------------------------------------------------------------------------------------------------------------------------

loginfo("Read background proteins file")
captured_output<-capture.output(
  background_list <- vroom::vroom( args$background_file, col_names=FALSE, delim="\t" )
  ,type = "message"
)
logdebug(captured_output)

## ---------------------------------------------------------------------------------------------------------------------------------------------------
captured_output<-capture.output(
  query_proteins <- vroom::vroom( args$query_proteins_file )  %>%
    mutate( comparison = "Analysis_Results")
  ,type = "message"
)
logdebug(captured_output)

## ---------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Check protein ID and protein sets label columns exists in input file.")

if ( ! args$qp_protein_id_column %in% colnames(query_proteins ) ) {
  logerror("Column '%s' is not found in the input table.",args$qp_protein_id_column )
}

if ( ! args$qp_protein_sets_column %in% colnames(query_proteins ) ) {
  logerror("Column '%s' is not found in the input table.",args$qp_protein_sets_column)
}

if(  ! args$qp_protein_id_column %in% colnames(query_proteins )  |
     ! args$qp_protein_sets_column %in% colnames(query_proteins )  ) {
  stop()

}


  # sets_list <- c(rep( "Set1",3), rep("Set2", 3), rep("Set3",3) )
  # genes_list <- LETTERS[1:9]
  #
  # query_proteins <- data.frame( uniprot_acc = genes_list, cluster = sets_list)
  #
  # args$qp_protein_sets_column <- "cluster"
  # args$qp_protein_id_column <- "uniprot_acc"

  temp_proteins_list <- query_proteins %>%
    group_nest( !!rlang::sym( args$qp_protein_sets_column), .key="data" ) %>%
    ungroup

  list_of_genes_list <-temp_proteins_list %>%pull(data )

  names(list_of_genes_list) <- temp_proteins_list %>%pull(!!rlang::sym(args$qp_protein_sets_column) )

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

  #print( args$min_gene_set_size)
  min_gene_set_size_list <- parseNumList(args$min_gene_set_size)
  max_gene_set_size_list <- parseNumList(args$max_gene_set_size)

 list_of_comparisons <- query_proteins %>%
   distinct(comparison) %>%
   pull(comparison)
 purrr::walk( list_of_comparisons, ~createDirIfNotExists( file.path(args$output_dir,. )) )


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

 input_params <- expand_grid(
   names_of_genes_list = names( list_of_genes_list),
   go_aspect=go_aspect_list,
   input_comparison = list_of_comparisons,
   min_size = min_gene_set_size_list,
   max_size = max_gene_set_size_list)


 input_params_updated <- input_params %>%
   mutate( input_table = purrr::map(names_of_genes_list,  \(x) list_of_genes_list[[x]] ))


 runOneGoEnrichmentInOutFunctionPartial <- purrr::partial ( runOneGoEnrichmentInOutFunction,
                  comparison_column = comparison,
                  protein_id_column = !!rlang::sym(args$qp_protein_id_column),
                  go_annot = go_annot,
                  background_list = background_list,
                  id_to_annotation_dictionary=id_to_annotation_dictionary,
                  annotation_id=!!rlang::sym(args$annotation_id),
                  protein_id=!!rlang::sym(args$protein_id),
                  aspect_column=args$aspect_column,
                  p_val_thresh=args$p_val_thresh)

 enrichment_result <- input_params_updated %>%
   mutate( enrichment_results = purrr::pmap( list( names_of_genes_list
                                                   , input_table
                                                   , go_aspect
                                                   , input_comparison
                                                   , min_size
                                                   , max_size) ,
                                             \(names_of_genes_list
                                               ,input_table
                                               , go_aspect
                                               , input_comparison
                                               , min_size
                                               , max_size) {runOneGoEnrichmentInOutFunctionPartial(
                                                 names_of_genes_list = names_of_genes_list,
                                                 input_table=input_table,
                                                 go_aspect=go_aspect,
                                                 input_comparison=input_comparison,
                                                 min_gene_set_size=min_size,
                                                 max_gene_set_size=max_size)} ) ) %>%
   dplyr::select( enrichment_results ) %>%
   unnest(cols=enrichment_results)

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

