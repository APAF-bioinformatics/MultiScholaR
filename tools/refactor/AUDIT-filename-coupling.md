# Refactor Coupling Audit

This report flags filename-coupled references to `R/*.R` that can break wave-based splitting.

## Summary

- `high` / `dev`: 38
- `low` / `docs`: 3

## Findings

- `high` `dev` [dev/test_lipid_app.R:6](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_app.R:6): source("R/allGenerics.R")
- `high` `dev` [dev/test_lipid_app.R:7](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_app.R:7): source("R/func_lipid_s4_objects.R")
- `high` `dev` [dev/test_lipid_app.R:8](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_app.R:8): source("R/func_general_plotting.R")
- `high` `dev` [dev/test_lipid_app.R:9](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_app.R:9): if (file.exists("R/func_general_s4_generics.R")) source("R/func_general_s4_generics.R")
- `high` `dev` [dev/test_lipid_app.R:12](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_app.R:12): source("R/mod_lipid_import.R")
- `high` `dev` [dev/test_lipid_app.R:13](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_app.R:13): source("R/mod_lipid_design.R")
- `high` `dev` [dev/test_lipid_app.R:14](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_app.R:14): source("R/mod_lipid_qc.R")
- `high` `dev` [dev/test_lipid_app.R:15](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_app.R:15): source("R/mod_lipid_norm.R")
- `high` `dev` [dev/test_lipid_app.R:16](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_app.R:16): source("R/mod_lipid_de.R")
- `high` `dev` [dev/test_lipid_app.R:20](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_app.R:20): source("R/mod_lipid_summary.R")
- `high` `dev` [dev/test_lipid_app.R:21](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_app.R:21): source("R/mod_lipidomics.R")
- `high` `dev` [dev/test_lipid_core.R:7](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_core.R:7): source("R/allGenerics.R")
- `high` `dev` [dev/test_lipid_core.R:10](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_core.R:10): source("R/func_lipid_s4_objects.R")
- `high` `dev` [dev/test_lipid_core.R:14](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_core.R:14): source("R/mod_lipid_import.R")
- `high` `dev` [dev/test_lipid_core.R:15](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_core.R:15): source("R/func_lipid_import.R")
- `high` `dev` [dev/test_lipid_core.R:16](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_core.R:16): source("R/mod_lipid_design.R")
- `high` `dev` [dev/test_lipid_core.R:17](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_core.R:17): source("R/mod_lipid_design_builder.R")
- `high` `dev` [dev/test_lipid_de.R:6](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_de.R:6): source("R/allGenerics.R")
- `high` `dev` [dev/test_lipid_de.R:7](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_de.R:7): source("R/func_lipid_s4_objects.R")
- `high` `dev` [dev/test_lipid_de.R:8](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_de.R:8): source("R/func_general_plotting.R")
- `high` `dev` [dev/test_lipid_de.R:9](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_de.R:9): if (file.exists("R/func_general_s4_generics.R")) source("R/func_general_s4_generics.R")
- `high` `dev` [dev/test_lipid_de.R:13](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_de.R:13): source("R/func_lipid_de.R")
- `high` `dev` [dev/test_lipid_de.R:14](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_de.R:14): source("R/mod_lipid_de.R")
- `high` `dev` [dev/test_lipid_qc.R:6](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_qc.R:6): source("R/allGenerics.R")
- `high` `dev` [dev/test_lipid_qc.R:7](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_qc.R:7): source("R/func_lipid_s4_objects.R")
- `high` `dev` [dev/test_lipid_qc.R:8](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_qc.R:8): source("R/func_general_plotting.R")
- `high` `dev` [dev/test_lipid_qc.R:10](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_qc.R:10): if (file.exists("R/func_general_s4_generics.R")) source("R/func_general_s4_generics.R")
- `high` `dev` [dev/test_lipid_qc.R:14](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_qc.R:14): source("R/func_lipid_qc.R")
- `high` `dev` [dev/test_lipid_qc.R:15](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_qc.R:15): source("R/func_lipid_norm.R")
- `high` `dev` [dev/test_lipid_qc.R:16](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_qc.R:16): source("R/mod_lipid_qc.R")
- `high` `dev` [dev/test_lipid_qc.R:17](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_qc.R:17): source("R/mod_lipid_qc_duplicates.R")
- `high` `dev` [dev/test_lipid_qc.R:18](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_qc.R:18): source("R/mod_lipid_qc_intensity.R")
- `high` `dev` [dev/test_lipid_qc.R:19](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_qc.R:19): source("R/mod_lipid_qc_itsd.R")
- `high` `dev` [dev/test_lipid_qc.R:20](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_qc.R:20): source("R/mod_lipid_qc_s4.R")
- `high` `dev` [dev/test_lipid_qc.R:21](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_qc.R:21): source("R/mod_lipid_norm.R")
- `high` `dev` [dev/test_lipid_s4.R:7](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_s4.R:7): if (file.exists("R/allGenerics.R")) source("R/allGenerics.R")
- `high` `dev` [dev/test_lipid_s4.R:9](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_s4.R:9): if (file.exists("R/func_general_s4_generics.R")) source("R/func_general_s4_generics.R")
- `high` `dev` [dev/test_lipid_s4.R:11](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_s4.R:11): source("R/func_lipid_s4_objects.R")
- `low` `docs` [docs/proteomics_gui_testthat_strategy.md:377](/home/doktersmol/Documents/MultiScholaR/docs/proteomics_gui_testthat_strategy.md:377): source_files = c("R/func_prot_import.R", "R/func_prot_qc.R",
- `low` `docs` [docs/proteomics_gui_testthat_strategy.md:378](/home/doktersmol/Documents/MultiScholaR/docs/proteomics_gui_testthat_strategy.md:378): "R/func_prot_rollup.R", "R/func_prot_norm.R",
- `low` `docs` [docs/proteomics_gui_testthat_strategy.md:379](/home/doktersmol/Documents/MultiScholaR/docs/proteomics_gui_testthat_strategy.md:379): "R/func_prot_da.R", "R/func_prot_annotation.R"),
