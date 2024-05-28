#!/usr/bin/env Rscript

# Author(s): Ignatius Pang, Pablo Galaviz
# Email: cmri-bioinformatics@cmri.org.au
# Children’s Medical Research Institute, finding cures for childhood genetic diseases


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



#Test if BioManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager",repos = "https://cran.csiro.au/")
  BiocManager::install()
}

# load pacman package manager
if (!require(pacman)) {
  install.packages("pacman",repos = "https://cran.csiro.au/")
  library(pacman)
}




## needed
p_load(mixOmics)
p_load(tidyverse)
p_load(devtools)
p_load(dplyr)
p_load(furrr)
p_load(GGally)
p_load(ggplot2)
p_load(ggpubr)
p_load(ggrepel)
p_load(Glimma)
p_load(glue)
p_load(GO.db)
p_load(gridExtra)
p_load(here)
p_load(knitr)
p_load(limma)
p_load(magrittr)
p_load(plotly)
p_load(purrr)
p_load(qvalue)
p_load(readxl)
p_load(rlang)
p_load(ruv)
p_load(statmod)
p_load(stringi)
p_load(tibble)
p_load(tidyr)
p_load(tidyselect)
p_load(vroom)
p_load(writexl)



## For parsing FASTA file for MaxQuant file parsing
# p_load(seqinr)
# p_load(janitor)

## Only needed for missing values imputation in command line scripts
# p_load(PhosR)

## For enrichment analysis
# p_load(clusterProfiler)

## For annotation scripts and Functional enrichment scripts
# p_load(UniProt.ws)
# p_load(httr)
# p_load(biomaRt)
# p_load(xml2)
# p_load(rvest)
# p_load(igraph)


## For running command line scripts
# p_load(logging)
# p_load(optparse)
# p_load(tictoc)
# p_load(svglite)
# p_load(configr)



devtools::document()


