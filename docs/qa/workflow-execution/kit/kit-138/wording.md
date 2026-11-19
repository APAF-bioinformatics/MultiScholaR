# Integration GUI: every user-facing string (KIT-138, KIT-139, KIT-140)

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
| tab_model | Model | tab title (KIT-140) |
| step_model | Model | stepper (KIT-140) |

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
| refuse_changed | A data set has changed since it was imported. Import the data again at the Import tab. | KIT-140: the check's message moved under "Technical detail" |

## Shown verbatim, not from the table (for the user's decision)

- **The data type** column shows the omic's own route name as saved (`LFQ`, `DIA`, `TMT`, `LIPIDOMICS`,
  `METABOLOMICS`, `TRANSCRIPTOMICS`).
- **Data set names** are the saved names (`lipidomics_LCMS_Pos`, INT-4).
- **The checks' own messages** (the contract reader's, KIT-061; the model layer's, KIT-062) are no longer the
  sentence on screen (KIT-140): each refusal shows its class's plain sentence (below), and the check's message only
  under a collapsed "Technical detail" line. Those messages still use internal words ("manifest", "sha256", "view");
  they are kept unchanged because tests and records rely on them.
- **Parameter values in the settings' names**: a parameter the workbook adds before it has a label here shows its own
  name (for example `likelihood`), and a value without a label shows the value.

## Refusals in plain words (KIT-140)

One sentence per refusal class, `integrationRefusalText()`; the check's own message under "Technical detail". The
tab's own sentences (Import's "Select data sets from two or more omics", the settings' "enter a whole number") are
shown as written.

| Key | Class | Text |
|---|---|---|
| technical_detail | | Technical detail |
| refusal_input_invalid | integration_input_invalid (the contract reader; a run its omic's design does not list; no biological-sample column) | This omic's saved files cannot be used as they are: they are incomplete, or were changed after they were saved. Save the omic again at its Summary tab, then search again. |
| refusal_input_stale | integration_input_stale | A data set has changed since it was imported. Import the data again at the Import tab. |
| refusal_input_duplicate | integration_input_duplicate | Two inputs are of the same omic, or two data sets have the same name. Keep one input per omic. |
| refusal_too_few_omics | integration_too_few_omics | Two or more omics are needed. |
| refusal_too_few_views | integration_too_few_views | Two or more data sets are needed. |
| refusal_inputs_missing | integration_inputs_missing | Nothing is imported yet. Import two or more omics at the Import tab first. |
| refusal_no_shared | integration_no_shared_samples | The data sets share no biological sample. Check that each omic names its samples the same way. |
| refusal_group_conflict | integration_group_conflict | A biological sample is in different groups in different omics. Check the groups in each omic's design. |
| refusal_view_empty | integration_view_empty | A data set has no feature measured in two or more samples, so nothing of it is left for the model. |
| refusal_parameters_invalid | integration_parameters_invalid | A model setting is not valid. Check the values entered. |
| refusal_mofa_unavailable | integration_mofa_unavailable (and the software check) | The model cannot be trained on this computer: the MOFA2 software, or the Python it uses, is not installed. |
| refusal_other | any other error | This step could not be completed; the technical detail below says why. |

## Model tab (KIT-140)

The settings are built from the model's parameter table (`integrationMofaParameterSpec()`): a label below per
parameter (`param_<name>`), per choice (`choice_<name>_<value>`) and optional help (`help_<name>`). The defaults are
the table's, never these strings'.

| Key | Text | Fills |
|---|---|---|
| model_heading | Train the factor model (MOFA2) | |
| model_intro | MOFA2 looks for a small number of factors that explain the variation the data sets share. Check the settings, then train the model. | |
| model_not_confirmed | Confirm the sample table at the Sample Matching tab first. | |
| readiness_heading | Software check | |
| readiness_ok | MOFA2 and the Python it uses are installed: the model can be trained on this computer. | |
| readiness_recheck | Check again | button |
| settings_heading | Model settings | |
| settings_matching_heading | Set by the confirmed sample table | the three matching rules, shown, not offered |
| settings_provisional | Provisional: the integration methods may still change this default. Values marked so are the current defaults, not final ones. | |
| mark_provisional | provisional default | after each provisional setting's label |
| mark_decided | decided | after precision's label (INT-6) |
| setting_fixed | %s: %s | label, value |
| setting_not_whole | %s: enter a whole number. | label |
| param_technical_replicates; choice _average | Runs of the same biological sample; averaged | INT-2 |
| param_internal_standards; choice _drop, _keep | Internal standards; removed; kept | INT-3 |
| param_missing_samples; choice _keep_missing, _drop_incomplete | Samples missing from a data set; kept, with that data set's values missing; left out | D4 |
| param_feature_selection; choice _none, _top_variance | Features used from each data set; all features; only the most variable features | |
| param_features_per_view; help | Number of most variable features kept per data set; Used only when only the most variable features are used. | |
| param_scaling; choice _none, _feature_zscore | Scaling of each feature; none: the values as saved; each feature scaled to mean 0 and standard deviation 1 | |
| param_scale_views | Give every data set the same total variance | tick box |
| param_center_groups | Centre each feature on its mean | tick box |
| param_num_factors | Number of factors | |
| param_convergence_mode; choice _fast, _medium, _slow | How closely training must settle before it stops; fast; medium; slow (the closest) | |
| param_maxiter | Maximum number of training rounds | |
| param_seed; help | Random seed; The same data, settings and seed give the same model. | |
| param_precision; choice _double, _single | Numeric precision; double: the same model on every computer; single: less memory for very large studies | INT-6 |
| help_precision | Double precision is the chosen setting: it gives the same model on every computer. Single precision uses less memory for very large studies, but the model can differ slightly from one computer to another. | |
| run_button | Train the model | button |
| run_working | Training the model. This can take several minutes. | progress and status |
| run_refused | The model was not trained. %s | the refusal's sentence |
| run_none | No model trained yet. | |
| run_done | Model trained at %s. | time |
| run_reused | These data and settings were trained before; the saved model is shown. | same settings again |
| run_settings_changed | The settings have changed since this model was trained. Train again to apply them. | |
| run_count | Models trained since the import: %d. Each one is kept in the project; the latest is shown. | count |
| result_heading | The trained model | |
| result_converged | Training settled after %d rounds (of at most %d). | rounds, maximum |
| result_not_converged | Training stopped at the maximum of %d rounds without settling. Allow more rounds or a less close setting, then train again. | maximum |
| result_factors | Factors in the model: %d (asked for: %d). | |
| result_samples | Biological samples in the model: %d; measured in every data set: %d. | |
| result_variance_heading | Variance explained per data set (all factors together) | |
| col_view, col_model_features, col_variance | Data set; Features in the model; Variance explained (%) | table headers |
| result_saved | The model file is saved in the project; its SHA-256 checksum is %s. | checksum |
| few_samples | (KIT-139's) | shown again under a model trained with fewer than 16 shared samples |
