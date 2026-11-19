# Integration GUI: every user-facing string (KIT-138, KIT-139, KIT-140, KIT-063, KIT-142)

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

## Sample Matching tab: correcting the match with an edited table (KIT-063 item A)

| Key | Text | Where |
|---|---|---|
| match_heading | Correct the match (optional) | heading under the sample table |
| match_intro | If two omics name the same biological sample differently, download the sample table, change the biological_sample column so that each run names its biological sample, and upload the edited table. Runs you leave unchanged keep the name their omic saved. The new match is shown above before you confirm it. | names the file's column as it is written |
| match_download | Download the sample table | button; the file is integration_sample_match.tsv (columns omic, run, saved_biological_sample, biological_sample) |
| match_upload | Upload the edited table (tab-separated) | file input |
| match_clear | Use the names as saved | button |
| match_active | Your uploaded table is in use. | |
| match_changed | Your table gives %d runs a different biological sample. | runs |
| match_unchanged | Your table gives no run a different biological sample; samples are matched as saved. | |
| match_unknown_omics, _unknown_runs, _unknown_names, _rows_empty | %d rows name an omic that is not imported; they are ignored. / %d rows name a run that is not in the data; they are ignored. / %d rows name a sample name that is not in the data; they are ignored. / %d rows name neither a run nor a sample name; they are ignored. | rows |
| match_refused | The uploaded table cannot be used. %s | the refusal's sentence |
| refusal_sample_match_invalid | Check that the table has the columns omic, run and biological_sample, that every run listed has a biological sample, and that each run appears only once. | the refusal's sentence; the check's own message under "Technical detail" |
| setting_sample_match | Sample match: your uploaded table, which gives %d runs a different biological sample | Model tab, the settings the confirmed table fixed |

## Results tab (KIT-063)

| Key | Text | Where |
|---|---|---|
| tab_results; step_results | " Results"; Results | tab title; stepper |
| results_heading | Results of the factor model | |
| results_intro | What the trained model found: how much of each data set's variation each factor explains, how the factors relate to the groups, which features weigh most on each factor, and how the factors relate to each other. | |
| results_not_trained | Train a model at the Model tab first. | before a model |
| results_model | Results of the model trained at %s. | time |
| results_refused | The results cannot be shown. %s | the refusal's sentence |
| results_provisional | How many features are listed, which factors are shown, how the factors are related to the groups and how the figures look are provisional: the integration methods may still change them. | above the save button |
| advisories_heading | Notes from the model | MOFA2's warnings |
| advisory_overall_level | %s rises and falls with each sample's overall level in at least one data set. Such a factor often reflects differences in overall level between samples, for example from normalisation, rather than biology. Read it with care. | factors; MOFA2's text under "Technical detail" |
| advisory_other | MOFA2 gave a warning while training; its text is under Technical detail. | any other MOFA2 warning |
| variance_heading; variance_intro | Variance explained; The percentage of each data set's variation that each factor explains. A factor that explains variation in several data sets is shared between them; one that explains variation in only one data set is specific to it. | |
| variance_total_heading | All factors together, per data set | |
| col_factor | Factor | table header |
| association_heading | Factors and the groups | |
| association_intro | For each factor: the share of its variation that lies between the groups (one-way analysis of variance of the factor's values on the group), the test's p-value, and the p-value adjusted for the number of factors tested (%s). | the adjustment's name |
| padjust_BH, _holm, _bonferroni, _none | Benjamini-Hochberg; Holm; Bonferroni; not adjusted | |
| claim_heading | How to read this | the claim limit (KIT-060 D10), beside the factors-against-groups figure |
| claim_descriptive | The factors were found without using the groups. How they relate to the groups describes the model; it is not a test of the experiment. | always |
| claim_group_corrected | In %s, unwanted variation was removed with the groups as replicate sets, which keeps the differences between groups by construction. Factors can then separate the groups more than the data alone would, so no p-value is given for the groups. | data sets; only then |
| claim_correction | %s: correction for unwanted variation %s | per omic, the Import tab's words |
| col_samples_tested, col_groups_n, col_group_r2, col_f, col_df, col_p, col_p_adjusted | Samples; Groups; Variation between groups (%); F; Degrees of freedom; p-value; Adjusted p-value | table headers |
| cell_not_tested; cell_cannot_test | not given (see How to read this); cannot be computed (fewer than two groups or too few samples) | table cells |
| factor_values_heading | Factor values per biological sample | |
| weights_heading | Features with the largest weights | |
| weights_intro | The %d features with the largest weights, in absolute value, on the chosen factor in the chosen data set. A feature's weight says how strongly it follows the factor; its sign says in which direction. The scaled weight divides by the largest weight of that factor in that data set, so 1 or -1 marks the strongest feature. | number of features |
| weights_factor; weights_view | Factor; Data set | the two choices |
| col_rank, col_feature, col_feature_id, col_weight, col_weight_scaled | Rank; Feature; Identifier; Weight; Scaled weight | table headers |
| correlation_heading; correlation_intro | How the factors relate to each other; The correlation (Pearson) of the factors' values across the biological samples. Each factor is meant to capture its own pattern; a correlation far from 0 means two factors share part of theirs. | |
| plot_dataset, plot_factor, plot_variance, plot_group, plot_factor_value, plot_weight, plot_weight_scaled, plot_correlation | Data set; Factor; Variance explained (%); Group; Factor value; Weight; Scaled weight; Correlation | figure axes and legends |
| plot_weights_title | %s, %s | factor, data set |
| save_heading; save_button | Save the results; Save tables and figures to the project | |
| save_done | Saved at %s: %d tables to %s and %d figures to %s. | time, counts, folders |
| save_refused | The results were not saved. %s | the refusal's sentence |

## Record tab, study parameters, methods paragraph, report and replay (KIT-142)

These strings live in their own table, `integrationRecordTexts()` in `R/utils_integration_record.R` (the tab's title and
step are `integrationGuiTexts()`'s). The same words are used by the Record tab, the saved `study_parameters.txt`, the
record's workbook (sheets `integration_*`) and the report (`inst/reports/integration/integration_report.rmd`). The
methods paragraph is a draft for the user: each sentence is filled from the record, nothing is hard-coded.
Words that appear: "checksum (SHA-256)" (the record's identity of a file), "provisional default", "data set".

| Key | Text | Where |
|---|---|---|
| tab_record; step_record | " Record"; Record | tab title; stepper |
| record_heading; record_intro | Record and report; Save the integration's analysis record (what was combined, how samples were matched, every setting, the software versions and the model's checksum), its study parameters and a report with a draft methods paragraph. | Record tab |
| record_button; record_working | Save the record and the report; Writing the record and the report. | button; progress |
| record_not_ready | Train a model at the Model tab and save its results at the Results tab first. | before a model |
| record_done | Saved at %s: %d files in %s. | time, files, folder |
| record_refused | The record was not saved. %s | the refusal's sentence |
| record_methods_heading | Methods paragraph (draft, from the record) | Record tab |
| record_no_model | There is no trained model to record. Train a model at the Model tab first. | refusal |
| status_decided, status_provisional_default, status_user_set, status_input | decided; provisional default; set by the user; from the data | a setting's status |
| legend_decided, legend_provisional_default, legend_user_set, legend_input | a decision recorded in the methods ledger; not expected to change / the current default; the integration methods may still change it / changed by the user from the default / taken from the imported data | report legend |
| report_title; report_intro | Multi-omics integration report; This report was written by MultiScholaR from the integration's analysis record. Every value below is read from the record. | report, study parameters |
| section_methods; methods_draft_note | Methods (draft); A draft for the methods section, written from the record: check the wording before using it. | |
| section_inputs, section_matching, section_samples, section_settings, section_model, section_variance, section_results, section_notes, section_legend | Data sets combined; Sample matching; Biological samples in the model; Settings; The trained model; Variance explained per data set (all factors together); Saved results; Notes from the model; What the status of a setting means | headings |
| none_recorded; not_recorded; yes; no | None recorded.; not recorded; yes; no | |
| col_omic, col_datasets_in_model, col_values, col_correction, col_group_corrected, col_checksum | Omic; Data sets in the model; Values; Correction for unwanted variation (RUV); Corrected with the groups as replicate sets; Checksum of the saved data (SHA-256) | inputs table |
| col_field, col_value, col_setting, col_name, col_status, col_settled_by, col_applies_to, col_file, col_view, col_variance, col_note, col_status_meaning | Item; Value; Setting; Internal name; Status; Settled by; Applies to; File; Data set; Variance explained (%); Note; Meaning | headers |
| applies_model; applies_results | model; results | settings table |
| values_corrected; values_normalised | corrected for unwanted variation; normalised, not corrected | inputs table |
| match_rule, match_rule_saved, match_rule_user | How samples were matched; by the biological-sample column each omic saved; by the biological-sample column each omic saved, corrected by an uploaded table (%d runs moved) | runs |
| match_table_checksum, match_samples, match_shared, match_confirmed_at, match_design_file, match_design_checksum, match_file | Checksum of the uploaded match (SHA-256); Biological samples in the model; Measured in every data set; Sample table confirmed at; Sample table saved as; Checksum of the saved sample table (SHA-256); Uploaded match saved as | matching table |
| model_mofa2 ... model_identity | MOFA2; mofapy2; Python; R; reticulate; Numeric precision; Training rounds; Training settled; Factors (kept, asked for); Checksum of the model file (SHA-256); Size of the model file (bytes); Model's identifier in the project; Checksum of the model's identity (data and settings) | model table |
| model_factors_value | %d of %d | kept, asked for |
| note_few_samples | Only %d biological samples are measured in every data set (fewer than %d). | samples, threshold |
| methods_inputs; methods_input_item | %d data sets from %d omics were combined: %s.; %s (%s; %s) | counts, list; omic, data sets, values |
| methods_input_corrected; methods_input_normalised | corrected for unwanted variation, %s; normalised, not corrected for unwanted variation | the approach in the Import tab's words |
| methods_software | Shared factors were learned with multi-omics factor analysis, MOFA2 %s (mofapy2 %s, Python %s, R %s), in %s precision. | versions, precision |
| methods_matching_saved; methods_matching_user | Samples were matched across omics by the biological-sample column each omic saved.; ... corrected by a user-supplied match table that assigned %d runs to another biological sample (SHA-256 %s). | runs, checksum |
| methods_replicates_average | Runs of the same biological sample were averaged (the mean of the observed values). | INT-2 |
| methods_standards_drop; methods_standards_keep | Internal standards were removed (%d features).; Internal standards were kept. | INT-3 |
| methods_missing_keep_missing; methods_missing_drop_incomplete | Biological samples missing from a data set were kept, with that data set's values missing.; Only the biological samples measured in every data set were kept. | D4 |
| methods_sparse | Features measured in fewer than two biological samples were left out. | |
| methods_selection_none; methods_selection_top_variance; methods_scaling_none; methods_scaling_feature_zscore | All remaining features were used; The %d most variable features of each data set were used; without per-feature scaling.; each scaled to mean 0 and standard deviation 1. | joined with ", " |
| methods_model | The model was asked for %d factors (convergence mode %s, at most %d training rounds, seed %d; each data set scaled to the same total variance: %s; features centred: %s). | |
| methods_trained_settled; methods_trained_not_settled | Training settled after %d rounds and kept %d factors; %d biological samples entered the model, %d of them measured in every data set.; Training stopped at the maximum of %d rounds without settling and kept %d factors; ... | |
| methods_association | Each factor was related to the experimental groups after training: a one-way analysis of variance of the factor's values on the group, with the share of the factor's variance that lies between the groups (R-squared) as the effect size, and p-values adjusted across the factors tested (%s). | only when results were saved |
| methods_claim_limit; methods_claim_none | In %s, unwanted variation had been removed with the groups as replicate sets, which keeps the differences between groups by construction; p-values for the association of the factors with the groups are therefore not reported, and the factors are described, not tested.; The factors were learned without the groups; their association with the groups describes the model. | D10 |
| methods_provisional; methods_user_set | Settings that were provisional defaults at the time of the analysis: %s.; Settings changed from their defaults: %s. | settings in the Model tab's words |
| replay_no_model | This record holds no integration model, so there is nothing to repeat. | replay refusal |
| replay_input_missing | The saved data of %s is no longer at %s, so the analysis cannot be repeated. Give its new location, or save it again at that omic's Summary tab. | omic, file |
| replay_input_changed | The saved data of %s has changed since this analysis was recorded (its checksum differs), so the analysis cannot be repeated on it. | omic |
| replay_mofa_unavailable | The analysis cannot be repeated on this computer: the MOFA2 software, or the Python it uses, is not installed. | |
| replay_versions_differ; replay_version_item | The analysis was recorded with other software versions (%s); it is repeated and compared with the recorded model numerically.; %s %s, here %s | warning |
