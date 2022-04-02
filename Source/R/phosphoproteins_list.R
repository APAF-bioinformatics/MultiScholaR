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

parser <- add_option(parser, c("-c","--config"), type = "character", default = "", dest = "config",
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


testRequiredArguments(args, c(
  "norm_phos_logfc_file"
))

testRequiredFiles(c(
  args$norm_phos_logfc_file))

args<-parseString(args,c("plots_format"))

if(isArgumentDefined(args,"plots_format"))
{
  args <- parseList(args,c("plots_format"))
}else {
  logwarn("plots_format is undefined, default output set to pdf.")
  args$plots_format <- list("pdf")
}

args <- setArgsDefault(args, "log_fc_column_name", as_func=as.character, default_val="norm_phos_logFC" )
args <- setArgsDefault(args, "fdr_column_name", as_func=as.character, default_val="combined_q_mod" )

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
  dplyr::filter( !!rlang::sym(args$combined_q_mod) < 0.05 &   !!rlang::sym(args$log_fc_column_name) > 0 ) %>%
  mutate( uniprto_acc_first = purrr::map_chr( uniprot_acc, ~str_split(., ":") %>% map_chr(1)))  %>%
  mutate( gene_name_first = purrr::map_chr( gene_name, ~str_split(., ":") %>% map_chr(1)))  %>%
  mutate( protein_name_first = purrr::map_chr( `PROTEIN-NAMES`, ~str_split(., ":") %>% map_chr(1)))  %>%
  group_by(comparison, uniprto_acc_first, gene_name_first, protein_name_first) %>%
  summarise( max_norm_phos_logFC = max(!!rlang::sym(args$log_fc_column_name))) %>%
  ungroup() %>%
  arrange( comparison, desc(max_norm_phos_logFC  ) ) 

vroom::vroom_write( positive_phosphoproteins, file.path( args$output_dir, "all_phosphoproteins_with_positive_logFC_sites.tab" ))

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# norm_phos_logFC

negative_phosphoproteins <- de_phos %>%
  dplyr::filter( !!rlang::sym(args$combined_q_mod) < 0.05 & !!rlang::sym(args$log_fc_column_name)  < 0 ) %>%
  mutate( uniprto_acc_first = purrr::map_chr( uniprot_acc, ~str_split(., ":") %>% map_chr(1)))  %>%
  mutate( gene_name_first = purrr::map_chr( gene_name, ~str_split(., ":") %>% map_chr(1)))  %>%
  mutate( protein_name_first = purrr::map_chr( `PROTEIN-NAMES`, ~str_split(., ":") %>% map_chr(1)))  %>%  
  group_by(comparison, uniprto_acc_first, gene_name_first, protein_name_first) %>%
  summarise( min_norm_phos_logFC = min(!!rlang::sym(args$log_fc_column_name))) %>%
  ungroup() %>%
  arrange( comparison, min_norm_phos_logFC  )

vroom::vroom_write( negative_phosphoproteins, file.path(args$output_dir, "all_phosphoproteins_with_negative_logFC_sites.tab" ) )

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

positive_only_phosphoproteins <- positive_phosphoproteins %>%
  anti_join( negative_phosphoproteins, by=c("uniprto_acc_first" = "uniprto_acc_first"))

vroom::vroom_write( positive_only_phosphoproteins, file.path(args$output_dir, "phosphoproteins_with_only_positive_logFC_sites.tab" ) )

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

negative_only_phosphoproteins <- negative_phosphoproteins  %>%
  anti_join( positive_phosphoproteins, by=c("uniprto_acc_first" = "uniprto_acc_first"))

vroom::vroom_write( negative_only_phosphoproteins, file.path(args$output_dir, "phosphoproteins_with_only_negative_logFC_sites.tab" ) )

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

all_phosphoproteins_with_significant_da_sites <- positive_phosphoproteins %>%
  bind_rows( negative_phosphoproteins) %>%
  distinct()

vroom::vroom_write( all_phosphoproteins_with_significant_da_sites, 
                    file.path(args$output_dir, "all_phosphoproteins_with_significant_da_sites.tab" ) )

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

background_phosphoproteins  <- de_phos %>%
  mutate( uniprto_acc_first = purrr::map_chr( uniprot_acc, ~str_split(., ":") %>% map_chr(1)))  %>%
  mutate( gene_name_first = purrr::map_chr( gene_name, ~str_split(., ":") %>% map_chr(1)))  %>%
  mutate( protein_name_first = purrr::map_chr( `PROTEIN-NAMES`, ~str_split(., ":") %>% map_chr(1)))  %>%  
  distinct(comparison, uniprto_acc_first, gene_name_first, protein_name_first) %>%
  arrange( comparison, gene_name_first  )

vroom::vroom_write( background_phosphoproteins, file.path(args$output_dir, "background_phosphoproteins.tab" ) )

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
te<-toc(quiet = TRUE)
loginfo("%f sec elapsed",te$toc-te$tic)
writeLines(capture.output(sessionInfo()), file.path(args$output_dir,"sessionInfo.txt"))

