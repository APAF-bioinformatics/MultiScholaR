# Integration GUI: every user-facing string (KIT-138, KIT-139)

Drafts for the user to set the final wording. All of them live in one table, `integrationGuiTexts()` in
`R/utils_integration_gui.R`; change a string there and nowhere else. `%s` and `%d` are filled in the order given.
Shown only with `options(MultiScholaR.integration_gui = TRUE)`.

Rules followed: plain words, exact numbers, none of "manifest", "owner", "declaration", "view", "alignment" on screen
(a "view" is called a "data set"; "alignment" is "sample matching").

## Home page

| Key | Text | Where |
|---|---|---|
| home_tile_title | Integration | the tile's heading (unchanged) |
| home_tile_text | Combine two or more omics already processed in this project | under the heading; replaces "Multi-omics integration (available when 2+ omics selected)" |

## Shell

| Key | Text | Where |
|---|---|---|
| tab_import | Import | tab title |
| step_import | Import | stepper |
| workflow_label | Integration | the shell's label |

## Import tab

| Key | Text | Fills |
|---|---|---|
| import_heading | Import processed omics | |
| import_intro | MultiScholaR looks in this project for the omics saved for integration at their Summary tab. Choose the data sets to combine; two or more omics are needed. | |
| import_folder | Folder searched: %s | the folder |
| import_rescan | Search again | button |
| import_found | Omics ready to combine: %d. Not ready: %d. | counts |
| import_none | No omic saved for integration was found in this folder. Open each omic's Summary tab, save its results for integration, then search again. | |
| omics_heading | Omics found | |
| col_omic, col_data_type, col_values, col_correction, col_saved, col_status | Omic; Data type; Values; Correction for unwanted variation (RUV); Saved; Status | table headers |
| values_corrected | corrected for unwanted variation | |
| values_normalised | normalised, not corrected | |
| correction_none | not applied | |
| correction_design | applied: %s | the RUV approach (below) |
| correction_legacy | applied: the approach of earlier versions (%s; replicate sets: %s) | method, replicate sets |
| correction_k | %s factor(s) | k |
| correction_view_applied | %s: %s factor(s) | assay, k (lipidomics, metabolomics) |
| correction_view_not_applied | %s: not applied | assay |
| correction_unknown | applied (approach not recorded) | |
| method_ruviii_c, method_ruviii | RUV-III-C; RUV-III | |
| sets_group, sets_technical_replicates | the RUV panel's texts (`omicGuiRuvRecordTexts`): "the experimental groups"; "the declared technical replicates" | |
| design_groups_all_features ... design_none | the RUV panel's texts (`omicGuiRuvRecordTexts`), e.g. "all features as controls, the experimental groups as replicates (recommended)" | under the RUV wording review already |
| saved_at | %s UTC | date and time |
| status_ready | ready | |
| status_changed | changed since it was imported: import again | |
| status_not_saved | not saved for integration: save it at that omic's Summary tab | |
| status_refused | cannot be used | |
| refused_heading | Not ready | box heading |
| refused_reason | %s: %s | omic, reason |
| views_heading | Data sets found | |
| col_view, col_features, col_runs, col_biological_samples, col_groups | Data set; Features; Runs; Biological samples; Groups | table headers |
| view_choice_label | Data sets to combine | checkbox group |
| view_choice_item | %s (%s; %d features, %d biological samples) | data set, omic, counts |
| import_button | Import selected data | button |
| import_need_two_omics | Select data sets from two or more omics; %d selected. | |
| import_need_two_views | Select two or more data sets; %d selected. | |
| import_need_ready | Only omics marked ready can be imported. | |
| import_refused | The import was refused: %s | the check's message |
| import_done | Imported %d omics: %d data sets, %d runs, %d biological samples. | |
| store_unavailable | Integration cannot keep its work in this project: %s | the reason |

## Shown verbatim, not from the table (for the user's decision)

- **The contract reader's refusals** (`readIntegrationInput`, KIT-061) appear as the reason of a "cannot be used" omic
  and after "The import was refused:". The brief asked for the reader's message; several of its sentences use
  internal words, for example "no integration input manifest at <path>", "the file matrix_x.tsv does not match its
  sha256", "the manifest's schema must be multischolar.integration_input.v1". Translating them means a second table
  keyed by the reader's checks; not done here.
- **The data type** column shows the omic's own route name as saved (`LFQ`, `DIA`, `TMT`, `LIPIDOMICS`,
  `METABOLOMICS`, `TRANSCRIPTOMICS`).
- **Data set names** are the saved names (`lipidomics_LCMS_Pos`, INT-4).
