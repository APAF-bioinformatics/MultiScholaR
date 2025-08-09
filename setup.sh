#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

# Update package lists
sudo apt-get update

# Install system dependencies
sudo apt-get install -y r-base r-base-dev libcurl4-openssl-dev libxml2-dev libssl-dev pandoc

# Install CRAN packages
R -e "install.packages(c('BiocManager', 'remotes', 'testthat', 'ggplot2', 'plotly', 'seqinr', 'UpSetR', 'tidyverse', 'GGally', 'logger', 'ggpubr', 'ggrepel', 'patchwork', 'svglite', 'iheatmapr', 'shiny', 'gprofiler2', 'webshot2', 'beepr', 'configr', 'devtools', 'DT', 'dplyr', 'flextable', 'furrr', 'future.apply', 'fs', 'gh', 'git2r', 'glue', 'gplots', 'gridExtra', 'gtools', 'gt', 'here', 'htmltools', 'httr', 'igraph', 'iq', 'janitor', 'lazyeval', 'logging', 'magrittr', 'multidplyr', 'openxlsx', 'optparse', 'pacman', 'progress', 'purrr', 'Rcpp', 'RcppEigen', 'RColorBrewer', 'readxl', 'rlang', 'RSpectra', 'rstudioapi', 'ruv', 'rvest', 'statmod', 'stringi', 'tibble', 'tictoc', 'tidyr', 'tidyselect', 'viridis', 'vroom', 'writexl', 'xml2', 'ggraph', 'reticulate', 'forcats', 'later', 'htmlwidgets', 'readr'), repos='http://cran.rstudio.com/')"

# Install Bioconductor packages
R -e "BiocManager::install(c('MOFA2', 'clusterProfiler', 'GO.db', 'GlimmaV2'))"

echo "Setup complete!"
