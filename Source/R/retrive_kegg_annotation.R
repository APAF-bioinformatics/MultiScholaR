#!/usr/bin/env Rscript


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

p_load(tidyverse)
p_load(plotly)
p_load(vroom)
p_load(writexl)
p_load(ggplot2)

p_load(magrittr)
p_load(xml2)
p_load(rvest)
p_load(ProteomeRiver)
p_load(httr)
p_load(stringi)
p_load(rvest)

p_load(tictoc)
p_load(configr)
p_load(logging)
p_load(optparse)

p_load("KEGGREST")
p_load("EnrichmentBrowser")


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
base_dir <-  here::here()
data_dir <- file.path( base_dir, "Data")
results_dir <- file.path(base_dir, "Results")
source_dir <- file.path(base_dir, "Source")



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Wen, Ruhui. (2021). Re: How i can get a list of KEGG pathways and its list of genes? . Retrieved from: https://www.researchgate.net/post/How_i_can_get_a_list_of_KEGG_pathways_and_its_list_of_genes/5ff5ff79bffe6552672ea3fc/citation/download.

#You can retrieve the gene set from KEGG pathways of an organism using the following codes in R. Here I took Staphylcoccus aureus MRSA252 (kegg code: sar) as an example.
#step 1: install related packages and open them in R:





tic()
set.seed(123456)

command_line_options <- commandArgs(trailingOnly = TRUE)
#Note: options with default values are ignored in the configuration file parsing.

parser <- OptionParser(add_help_option = TRUE)

#Note: options with default values are ignored in the configuration file parsing.
parser <- add_option(parser, c( "--species"), type = "character",
                     help = "Species for which to download annotation. Use the KEGG three letter species code.")

parser <- add_option(parser, c( "--gene_id_file"), type = "character",
                     help = "A file with a list of gene ID for which to retrieve KEGG annotation for.")

parser <- add_option(parser, c( "--results_dir"), type = "character",
                     help = "Results Directory")

#parse comand line arguments first.
args <- parse_args(parser)

# args <- setArgsDefault(args, "species", as_func=as.character, default_val=NA )
# args <- setArgsDefault(args, "results_dir", as_func=as.character, default_val="./Results" )
# args <- setArgsDefault(args, "gene_id_file", as_func=as.character, default_val="gene_ids_table_python_all.tab" )

args <- setArgsDefault(args, "species", as_func=as.character, default_val="hsa" )
args <- setArgsDefault(args, "results_dir", as_func=as.character, default_val="/home/ubuntu/Workings/2023/rett_pbmcs_wendy_gold_bmp_12_20230201/Data/KEGG/" )
args <- setArgsDefault(args, "gene_id_file", as_func=as.character, default_val="/home/ubuntu/Workings/2023/rett_pbmcs_wendy_gold_bmp_12_20230201/Results/UniProt/gene_ids_table_python_all.tab" )

testRequiredFiles(c(
  args$gene_id_file
))

createDirIfNotExists(file.path( args$results_dir ) )

gene_id_to_uniprot_acc_file <- file.path( args$gene_id_file)

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#step2: check and obtain a list of entry identifiers (in this case: sar) and associated definition for a given database or a given set of database entries.
identifiers <- keggList("pathway", args$species)
saveRDS( identifiers, file.path( args$results_dir,  paste0("identifiers_kegg_", args$species,".RDS") ) )

#step 3: download the pathways of that organism:
# /home/ignatius/.cache/R/EnrichmentBrowser
# pathways_download <- downloadPathways(args$species, out.dir= file.path( data_dir, "KEGG")) ## only use this line if you want to download all the XML files locally
pathways <- downloadPathways(args$species)
# pathways$mmu05416

# purrr::map( pathways, \(x) {x@pathwayInfo@name})
# purrr::map( pathways, \(x) {x@pathwayInfo@title})

#
# parsePathway <- function(pathway) {
#
#   pathway@pathwayInfo@name
#
#   pathway@pathwayInfo@org
#
#   pathway@pathwayInfo@number
#
#   purrr::map( pathway@nodes, parseNode) }
#
# parseNode <- function(node) {node@name}
#
# purrr::map( pathways, parsePathway)

#step 4: retrieve gene sets for an organism from databases such as GO and KEGG:
gene_sets <- getGenesets(org = args$species, db = "kegg", cache = TRUE, return.type="list")
saveRDS( gene_sets, file.path( args$results_dir,  paste0("gene_sets_kegg_", args$species,".RDS") ) )

gene_sets <- readRDS(file.path( args$results_dir,  paste0("gene_sets_kegg_", args$species,".RDS")) )

#step5: Parse and write the gene sets to a flat text file in GMT format for other pathway enrichment analysis programs (e.g., GSEA):
writeGMT(gene_sets, gmt.file = file.path( args$results_dir,  paste0("20210106_kegg_", args$species, "_gmt") ) )


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pathways_list <- data.frame( pathway_id = names(pathways),
            pathway_name = purrr::map_chr( names(pathways), ~pathways[[.]]@pathwayInfo@title))

vroom::vroom_write( pathways_list, file.path( args$results_dir,  "pathways_list.tab") )

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# need a gene_id column and pathway_id column
uniprot_acc_to_gene_id <- vroom::vroom( gene_id_to_uniprot_acc_file, col_types ="cc" ) %>%
  dplyr::filter( !is.na(gene_id ))

uniprot_acc_to_gene_id


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
gene_sets_tbl <- data.frame( pathway_id_name= names ( gene_sets )  ) %>%
  mutate( gene_id = purrr::map( pathway_id_name, ~gene_sets[[.]]) ) %>%
  unnest(gene_id) %>%
  left_join( uniprot_acc_to_gene_id, by="gene_id") %>%
  mutate( pathway_id =  str_split( pathway_id_name, pattern="_") %>% purrr::map_chr(1)) %>%
  left_join( pathways_list, by = "pathway_id") %>%
  dplyr::select(-pathway_id_name, - gene_id) %>%
  dplyr::filter(!is.na(uniprot_acc))

vroom::vroom_write( gene_sets_tbl, file.path( args$results_dir,  "gene_sets.tab" ))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
te<-toc(quiet = TRUE)
loginfo("%f sec elapsed",te$toc-te$tic)
writeLines(capture.output(sessionInfo()), file.path(args$results_dir,"sessionInfo.txt"))

