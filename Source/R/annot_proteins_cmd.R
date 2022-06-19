#!/usr/bin/env Rscript

# Author(s): Ignatius Pang, Pablo Galaviz
# Email: cmri-bioinformatics@cmri.org.au
# Children’s Medical Research Institute, finding cures for childhood genetic diseases


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Test if BioManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
   BiocManager::install()
}

# load pacman package manager
if(!require(pacman)){
    install.packages("pacman")
    library(pacman)
}
p_load(optparse)
p_load(tictoc)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
p_load(tidyverse)
p_load(vroom)
p_load(magrittr)
p_load(rlang)
p_load( janitor)
p_load(UniProt.ws)
p_load(biomaRt)
p_load(GO.db)
p_load(ProteomeRiver)
p_load(configr)
p_load(logging)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
tic()


command_line_options <- commandArgs(trailingOnly = TRUE)
#Note: options with default values are ignored in the configuration file parsing.

parser <- OptionParser(add_help_option =TRUE)

#Note: options with default values are ignored in the configuration file parsing.
parser <- add_option(parser, c("-d", "--debug"), action = "store_true", default = FALSE,
                     help = "Print debugging output")

parser <- add_option(parser, c("-s", "--silent"), action = "store_true", default = FALSE,
                     help = "Only print critical information to the console.")

parser <- add_option(parser, c("-n", "--no_backup"), action = "store_true", default = FALSE,
                     help = "Deactivate backup of previous run.")

parser <- add_option(parser, c("-c", "--config"), type = "character", default = "config_prot.ini", dest = "config",
                     help = "Configuration file.",
                     metavar = "string")

parser <- add_option(parser, c("-o", "--output_dir"), type = "character", dest = "output_dir",
                     help = "Directory path for all results files.",
                     metavar = "string")

parser <- add_option(parser, c("-t","--tmp_dir"), type = "character", default = "cache", dest = "tmp_dir",
                     help = "Directory path for temporary files.",
                     metavar = "string")

parser <- add_option(parser, c("-l", "--log_file"), type = "character", default = "output.log", dest = "log_file",
                     help = "Name of the logging file.",
                     metavar = "string")

#Options without a default value have the following priority: configuration file < command line argument
parser <- add_option(parser,  "--taxonomy_id", type="integer", dest = "taxonomy_id",
                       help="The NCBI taxonomy ID of the organism being investigated (e.g. M. musculus=10090, H. sapien=9606).",
                       metavar="integer")

parser <- add_option(parser, "--input_wide_file", type="character", dest = "input_wide_file",
                       help="Results table with values in wider format.",
                       metavar="string")

parser <- add_option(parser, "--input_long_file", type="character",  dest = "input_long_file",
                     help="Results table with values in longer format.",
                     metavar="string")

parser <- add_option(parser, "--ids_file", type="character", dest = "ids_file",
                     help="File to link the cleaned list of accessions to the original list of protein groups from MaxQuant file. Also, contains the MaxQuant output row ID.",
                     metavar="string")

parser <- add_option(parser, "--raw_counts_file", type="character", dest = "raw_counts_file",
                     help="Input file with the protein abundance data.",
                     metavar="string")

parser <- add_option(parser,  "--output_wide_file", type="character",  dest = "output_wide_file",
                     help="Results table with values in wider format.",
                     metavar="string")

parser <- add_option(parser,  "--output_long_file", type="character",  dest = "output_long_file",
                     help="Results table with values in longer format.",
                     metavar="string")

parser <- add_option(parser,  "--reactome_file", type="character",  dest = "reactome_file",
                     help="Name of the reactome data (Download and save if it does not exists). ",
                     metavar="string")

parser <- add_option(parser,  "--uniprot_file", type="character",  dest = "uniprot_file",
                     help="Name of the uniprot data (Download and save if it does not exists). ",
                     metavar="string")

#parse comand line arguments first.
args <- parse_args(parser)

#parse and merge the configuration file options.
if (args$config != "") {
   args <- config.list.merge(eval.config(file = args$config, config = "annot_proteins"), args)
}

args <- setArgsDefault(args, "output_dir", as_func=as.character, default_val="annot_proteins" )


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
  "taxonomy_id"
  ,"output_wide_file"
  ,"output_long_file"
  ,"ids_file"
  ,"input_wide_file"
  ,"input_long_file"
  ,"raw_counts_file"
  ,"reactome_file"
  ,"uniprot_file"
))


testRequiredFiles(c(
  args$input_wide_file
  , args$input_long_file
  , args$raw_counts_file
  ,args$ids_file
))

#unused
#args<-parseType(args,
#  c("some_double_par")
#  ,as.double)

args<-parseType(args,
  c("taxonomy_id")
                ,as.integer)



loginfo("Read file with results table with values in wider format %s", args$input_wide_file)
captured_output<-capture.output(
  de_proteins_wider <- vroom::vroom( args$input_wide_file )
    , type = "message"
)
logdebug(captured_output)

loginfo("Read file with results table with values in longer format %s", args$input_long_file)
captured_output<-capture.output(
  de_proteins_longer <- vroom::vroom( args$input_long_file )
    , type = "message"
)
logdebug(captured_output)

loginfo("Read file to link the cleaned list of accessions to the original list of protein groups from MaxQuant file %s", args$ids_file)
captured_output<-capture.output(
  ids_tbl <- vroom::vroom( args$ids_file )
    , type = "message"
)
logdebug(captured_output)

loginfo("Read file with the protein abundance data %s", args$raw_counts_file)
captured_output<-capture.output(
  dat_tbl <- vroom::vroom( args$raw_counts_file )
    , type = "message"
)
logdebug(captured_output)

dat_cln <-  janitor::clean_names(dat_tbl) %>%
  dplyr::rename( gene_names_maxquant = "gene_names")
colnames(dat_cln) <- str_replace(colnames(dat_cln), "_i_ds", "_ids" )



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


reactome_file<-file.path(args$tmp_dir,args$reactome_file)
if(!file.exists(reactome_file))
{
  logwarn("Download Reactome UniProt to pathways file.")
  status <- download.file(url="https://reactome.org/download/current/UniProt2Reactome.txt", destfile=reactome_file)
  loginfo(status)
}
loginfo("Reading Reactome UniProt to pathways file.")
captured_output<-capture.output(
  reactome_map <- vroom::vroom( reactome_file ,
                              col_names = c("uniprot_acc", "reactome_id", "url", "reactome_term", "evidence", "organism") )
    , type = "message"
)
logdebug(captured_output)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

loginfo("Get the best UniProt accession per row.")

uniprot_acc_tbl <- de_proteins_wider %>%
  mutate( uniprot_acc_copy = uniprot_acc ) %>%
  separate_rows(uniprot_acc_copy, sep=":" ) %>%
  mutate( join_uniprot_acc = cleanIsoformNumber(uniprot_acc_copy)) %>%
  dplyr::distinct( uniprot_acc, join_uniprot_acc) %>%
  group_by( uniprot_acc) %>%
  mutate( acc_order_id = row_number()) %>%
  ungroup




## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Download information from UniProt.")
uniprot_file<-file.path(args$tmp_dir,args$uniprot_file)
if( ! file.exists( uniprot_file )) {

  up <- UniProt.ws(taxId=args$taxonomy_id )
  list_of_sp_columns <- c("EXISTENCE"
                          , "SCORE"
                          , "REVIEWED"
                          , "GENENAME"
                          , "PROTEIN-NAMES"
                          , "LENGTH"
                          , "ENSEMBL"
                          , "GO-ID"
                          , "KEYWORDS"
  )
  up_cls<-unlist(columns(up))
  list_intersect<-intersect(list_of_sp_columns,up_cls)
  if(length(setdiff( list_of_sp_columns,list_intersect)) > 0)
  {
    logerror("UniProt fields not found: %s",paste(list_of_sp_columns[,list_intersect],sep=", "))
  }
  uniprot_dat <- batchQueryEvidence(uniprot_acc_tbl, join_uniprot_acc, uniprot_handle=up,
                                uniprot_columns = list_of_sp_columns)
  saveRDS( uniprot_dat, uniprot_file)

}

loginfo("Reading UniProt data from %s.",uniprot_file)
uniprot_dat <- readRDS( uniprot_file )



## ----eval=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##  # uniprot_dat <- batch_query_evidence(uniprot_acc_tbl %>% head(100), best_uniprot_acc, uniprot_handle=up,
##  #                                      uniprot_columns = list_of_sp_columns)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

loginfo("Merge with Gene Ontology terms.")
goterms <- Term(GOTERM)
gotypes <- Ontology(GOTERM)


uniprot_dat_cln <- uniprotGoIdToTerm(uniprot_dat, sep="; ", goterms, gotypes  )


uniprot_dat_multiple_acc <- uniprot_acc_tbl %>%
  left_join( uniprot_dat_cln, by=c("join_uniprot_acc" = "UNIPROTKB") )   %>%
  arrange( uniprot_acc, acc_order_id) %>%
  group_by(uniprot_acc ) %>%
  summarise( across( .cols=setdiff( colnames( uniprot_dat_cln), "UNIPROTKB")   , ~paste(., collapse=":"))   ) %>%
  ungroup()   %>%
  dplyr::rename( UNIPROT_GENENAME = "GENENAME")


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Add reactome pathways annotation.")
reactome_term_tbl <- uniprot_acc_tbl %>%
  left_join( reactome_map, by=c("join_uniprot_acc" = "uniprot_acc") )   %>%
  dplyr::filter(reactome_term != "NA" ) %>%
  group_by(uniprot_acc, join_uniprot_acc) %>%
  summarise( reactome_term = paste(reactome_term, collapse="; ") ) %>%
  ungroup()     %>%
  mutate(reactome_term = str_replace_all( reactome_term , ":", "-")) %>%
  group_by(uniprot_acc ) %>%
  summarise( reactome_term = paste(reactome_term, collapse=":") ) %>%
  ungroup()

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# print("Output wider format results table with protein annotation.")
# de_proteins_wider_annot <- de_proteins_wider %>%
#   left_join( ids_tbl, by=c("uniprot_acc" = "uniprot_acc") ) %>%
#   left_join( uniprot_dat_multiple_acc, by = c("uniprot_acc" = "uniprot_acc") ) %>%
#   left_join( reactome_term_tbl, by = c("uniprot_acc" = "uniprot_acc"))  %>%
#   left_join( dat_cln, by=c("maxquant_row_id" = "id",
#                            "protein_ids" = "protein_ids")) %>%
#   arrange( comparison, q.mod, log2FC) %>%
#   distinct()
#
# head( de_proteins_wider_annot )
#
# vroom::vroom_write(de_proteins_wider_annot, path=file.path(args$output_dir,args$output_wide_file ))

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Output longer format results table with protein annotation.")
de_proteins_longer_annot <- de_proteins_longer %>%
  left_join( ids_tbl, by=c("uniprot_acc" = "uniprot_acc") ) %>%
  left_join( uniprot_dat_multiple_acc, by = c("uniprot_acc" = "uniprot_acc") ) %>%
  left_join( reactome_term_tbl, by = c("uniprot_acc" = "uniprot_acc"))  %>%
  mutate( maxquant_row_id = as.character(maxquant_row_id)) %>%
  left_join( dat_cln %>%
             mutate( id = as.character(id)) , by=c("maxquant_row_id" = "id",
                           "protein_ids" = "protein_ids")) %>%
  arrange( comparison, q.mod, log2FC) %>%
  distinct() %>%
  relocate( UNIPROT_GENENAME,	`PROTEIN-NAMES`, .before="left_group" )

vroom::vroom_write(de_proteins_longer_annot, file.path(args$output_dir,args$output_long_file ) )

list_of_long_columns <- intersect(colnames(de_proteins_longer_annot), c( "uniprot_acc",
                                                                         "protein_ids",
                                                                         "protein_names",
                                                                         "ENSEMBL",
                                                                         "PROTEIN-NAMES",
                                                                         "KEYWORDS",
                                                                         "GO-ID",
                                                                         "go_biological_process",
                                                                         "go_cellular_compartment",
                                                                         "go_molecular_function",
                                                                         "reactome_term",
                                                                         "majority_protein_ids",
                                                                         "fasta_headers",
                                                                         "evidence_ids",
                                                                         "ms_ms_ids") )

# list_of_long_columns <- colnames(de_proteins_longer_annot)

writexl::write_xlsx(de_proteins_longer_annot %>%
                      mutate_at( list_of_long_columns, ~substr(., 1, 32760) ),
                    file.path(args$output_dir, str_replace(args$output_long_file, "\\..*", ".xlsx")  ) )


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
te<-toc(quiet = TRUE)
loginfo("%f sec elapsed",te$toc-te$tic)
writeLines(capture.output(sessionInfo()), file.path(args$output_dir,"sessionInfo.txt"))

