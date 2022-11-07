#!/usr/bin/env Rscript

#Test if BioManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# load pacman package manager
if(!require(pacman)){
  install.packages("pacman")
  library(pacman)
}


p_load(vroom)
p_load(here)
p_load(ProteomeRiver)

base_dir <- here::here()
results_dir <- file.path(base_dir, "Results")

createDirIfNotExists(file.path(results_dir, "Reactome"))

reactome_file <- file.path(results_dir, "Reactome", "reactome_data.txt" )
if(!file.exists(reactome_file))
{
  status <- download.file(url="https://reactome.org/download/current/UniProt2Reactome.txt", destfile=reactome_file)
}
reactome_map <- vroom::vroom( reactome_file ,
                              col_names = c("uniprot_acc", "reactome_id", "url", "reactome_term", "evidence", "organism") )



vroom::vroom_write( reactome_map, file.path(results_dir, "Reactome", "reactome_data_with_header.txt"))
