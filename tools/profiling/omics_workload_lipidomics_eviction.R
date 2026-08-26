args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]))
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
source(file.path(
    repo_root,
    "tools",
    "profiling",
    "omics_workload_assay_eviction.R"
))
omicsAssayEvictionMain(
    args,
    omic = "lipidomics",
    schema = "multischolar.lipidomics_eviction_measurement",
    script_path = script_path,
    repo_root = repo_root
)
