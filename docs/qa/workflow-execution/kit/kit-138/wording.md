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
| tab_alignment | Sample Matching | tab title (KIT-139) |
| step_alignment | Samples | stepper (KIT-139) |

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

## Sample Matching tab (KIT-139)

| Key | Text | Fills |
|---|---|---|
| align_heading | Match samples across omics | |
| align_intro | Each row is one biological sample, matched across omics by the sample column each omic saved. Check that every row is one biological sample, then confirm the table. | |
| align_not_imported | Import two or more omics at the Import tab first. | |
| align_key_heading | How samples are matched | |
| align_key | %s: by the column %s | omic, column |
| align_caution | Samples are matched by name: a sample must carry the same name in every omic. A sample named differently in two omics appears as two rows. | |
| align_rules_heading | What happens before the model | |
| align_provisional | These rules are provisional: the integration methods may still change them. | |
| fact_samples | %d biological samples in all; %d measured in every data set. | |
| fact_missing_keep_missing | %d biological samples are missing from at least one data set; they are kept, with that data set's values missing. | D4 |
| fact_missing_drop_incomplete | %d biological samples missing from at least one data set are left out. | the other choice of the parameter |
| fact_replicates_average | %s: %d runs from %d biological samples; the runs of the %d samples measured more than once are averaged. | INT-2 |
| fact_replicates_single | %s: %d runs, one per biological sample; nothing to average. | |
| fact_standards_drop | %s: %d internal standards removed. | INT-3 |
| fact_standards_keep | %s: %d internal standards kept. | the other choice |
| fact_sparse | %s: %d features measured in fewer than two samples are left out. | the model's rule |
| few_samples | Only %d biological samples are measured in every data set (fewer than %d). A factor model needs more samples; its results may be unstable. | D7 (a warning) |
| view_table_heading | Data sets | |
| col_standards_removed, col_sparse, col_features_used | Internal standards removed; Features left out (fewer than two samples); Features used | table headers (with col_view, col_runs, col_biological_samples) |
| sample_table_heading | Biological samples | |
| col_sample, col_group, col_views_present | Biological sample; Group; Data sets measured | table headers; then one column per data set |
| cell_missing | missing | a data set without that sample |
| confirm_button | Confirm this sample table | button |
| confirmed | Sample table confirmed at %s and saved to %s. | time, file |
| not_confirmed | Not confirmed yet. | |
| align_refused | The samples cannot be matched: %s | reason |
| refuse_no_shared | These data sets share no biological sample (%s). Check that each omic names its samples the same way. | data sets |
| refuse_too_few_views | Two or more data sets are needed; %d selected. | |
| refuse_changed | A data set has changed since it was imported: %s | the check's message |

## Shown verbatim, not from the table (for the user's decision)

- **The contract reader's refusals** (`readIntegrationInput`, KIT-061) appear as the reason of a "cannot be used" omic
  and after "The import was refused:". The brief asked for the reader's message; several of its sentences use
  internal words, for example "no integration input manifest at <path>", "the file matrix_x.tsv does not match its
  sha256", "the manifest's schema must be multischolar.integration_input.v1". Translating them means a second table
  keyed by the reader's checks; not done here.
- **The data type** column shows the omic's own route name as saved (`LFQ`, `DIA`, `TMT`, `LIPIDOMICS`,
  `METABOLOMICS`, `TRANSCRIPTOMICS`).
- **Data set names** are the saved names (`lipidomics_LCMS_Pos`, INT-4).
- **The model layer's other refusals** of the sample matching (a biological sample in two groups, a run its omic's
  design does not list, a data set left without features) and the changed-input detail after `refuse_changed` are
  KIT-062's messages, shown as they are; they use "view", "biological-sample column" and "manifest's sha256".
