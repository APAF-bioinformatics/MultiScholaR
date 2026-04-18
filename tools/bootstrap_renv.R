#!/usr/bin/env Rscript

repo_root <- normalizePath(".", mustWork = TRUE)

Sys.setenv(
  RENV_PATHS_ROOT = file.path(repo_root, ".renv-paths"),
  RENV_PATHS_CACHE = file.path(repo_root, ".renv-paths", "cache"),
  RENV_PATHS_LIBRARY = file.path(repo_root, "renv", "library"),
  RENV_CONFIG_SANDBOX_ENABLED = "FALSE"
)

source(file.path(repo_root, "renv", "activate.R"))

suppressPackageStartupMessages({
  library(renv)
})

base_packages <- c(
  "base", "compiler", "datasets", "graphics", "grDevices", "grid",
  "methods", "parallel", "splines", "stats", "stats4", "tcltk",
  "tools", "utils"
)

bioc_packages <- c(
  "ComplexHeatmap", "Glimma", "GlimmaV2", "limma", "qvalue"
)

remote_packages <- c(
  RUVIIIC = "cran/RUVIIIC",
  GlimmaV2 = "APAF-bioinformatics/GlimmaV2"
)

read_declared_packages <- function(path = file.path(repo_root, "DESCRIPTION")) {
  desc <- read.dcf(path)[1, ]
  fields <- intersect(c("Imports", "Suggests"), colnames(read.dcf(path)))
  raw_values <- unlist(strsplit(paste(desc[fields], collapse = ","), ","))
  packages <- trimws(gsub("[[:space:]]+", " ", raw_values))
  packages <- sub(" .*", "", packages)
  packages <- sort(unique(packages[nzchar(packages)]))
  setdiff(packages, base_packages)
}

declared <- read_declared_packages()
default_user_lib <- file.path(
  path.expand("~"),
  "R",
  sprintf("%s-library", R.version$platform),
  sprintf("%s.%s", R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1])
)
global_libs <- unique(c(
  default_user_lib,
  .Library.site,
  .Library
))
global_libs <- global_libs[nzchar(global_libs) & dir.exists(global_libs)]

installed_global <- character()
for (lib in global_libs) {
  installed_global <- union(installed_global, rownames(installed.packages(lib.loc = lib)))
}

hydrate_pkgs <- intersect(declared, installed_global)
install_pkgs <- setdiff(declared, hydrate_pkgs)

cat("Hydrating ", length(hydrate_pkgs), " package(s) from existing libraries.\n", sep = "")
if (length(hydrate_pkgs)) {
  renv::hydrate(packages = hydrate_pkgs, prompt = FALSE)
}

cran_pkgs <- setdiff(install_pkgs, bioc_packages)
bioc_pkgs <- intersect(install_pkgs, bioc_packages)
remote_pkgs <- intersect(install_pkgs, names(remote_packages))
cran_pkgs <- setdiff(cran_pkgs, remote_pkgs)
bioc_pkgs <- setdiff(bioc_pkgs, remote_pkgs)

if (length(cran_pkgs)) {
  cat("Installing ", length(cran_pkgs), " CRAN package(s).\n", sep = "")
  renv::install(cran_pkgs, prompt = FALSE, transactional = FALSE)
}

if (length(bioc_pkgs)) {
  cat("Installing ", length(bioc_pkgs), " Bioconductor package(s).\n", sep = "")
  renv::install(paste0("bioc::", bioc_pkgs), prompt = FALSE, transactional = FALSE)
}

if (length(remote_pkgs)) {
  cat("Installing ", length(remote_pkgs), " remote package(s).\n", sep = "")
  renv::install(unname(remote_packages[remote_pkgs]), prompt = FALSE, transactional = FALSE)
}

cat("Snapshotting lockfile.\n")
renv::snapshot(prompt = FALSE)
