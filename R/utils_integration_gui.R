# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# KIT-138, KIT-139: the integration GUI's configuration, its wording table and the helpers its tabs call (the project
# scan, the import request, the sample-matching preview). The model layer is KIT-062's
# (utils_workflow_integration_model.R); nothing here computes beyond reading what it returns. The whole surface sits
# behind the option `MultiScholaR.integration_gui` (off by default): with it off, the home page and the app are today's,
# byte for byte. The single-omic `validateOmicGuiConfig` does not fit (it requires design, qc and da groups); the shell
# is validated by `validateOmicGuiShell` only (KIT-060 audit section 5.3).

#' Whether the integration GUI is offered (KIT-138; the gate in the manner of `MultiScholaR.ruv_design_panel`).
#' `options(MultiScholaR.integration_gui = TRUE)` shows it.
#' @return One logical.
#' @noRd
integrationGuiEnabled <- function() {
    isTRUE(getOption("MultiScholaR.integration_gui", FALSE))
}

#' Every user-facing string of the integration GUI (KIT-138, KIT-139): one table, drafts the user sets the final
#' wording of (docs/qa/workflow-execution/kit/kit-138/wording.md). The RUV approaches' names are the RUV panel's
#' (`omicGuiRuvRecordTexts`), so that one approach reads the same everywhere. `%s`/`%d` are filled in the order the
#' comment beside each entry gives.
#' @return A named character vector.
#' @noRd
integrationGuiTexts <- function() {
    ruv <- omicGuiRuvRecordTexts()
    c(
        # The home page (gate on).
        home_tile_title = "Integration",
        home_tile_text = "Combine two or more omics already processed in this project",
        # The tabs and the stepper.
        tab_import = " Import",
        step_import = "Import",
        tab_alignment = " Sample Matching",
        step_alignment = "Samples",
        workflow_label = "Integration",
        # The Import tab.
        import_heading = "Import processed omics",
        import_intro = paste("MultiScholaR looks in this project for the omics saved for integration at their Summary",
            "tab. Choose the data sets to combine; two or more omics are needed."),
        # %s: the folder.
        import_folder = "Folder searched: %s",
        import_rescan = "Search again",
        # %d: ready, then not ready.
        import_found = "Omics ready to combine: %d. Not ready: %d.",
        import_none = paste("No omic saved for integration was found in this folder. Open each omic's Summary tab, save",
            "its results for integration, then search again."),
        omics_heading = "Omics found",
        col_omic = "Omic",
        col_data_type = "Data type",
        col_values = "Values",
        col_correction = "Correction for unwanted variation (RUV)",
        col_saved = "Saved",
        col_status = "Status",
        values_corrected = "corrected for unwanted variation",
        values_normalised = "normalised, not corrected",
        # KIT-153: an omic that applied RUV saves both.
        values_both = "normalised, and corrected for unwanted variation (both saved)",
        correction_none = "not applied",
        # %s: the approach.
        correction_design = "applied: %s",
        # %s: the method, then the replicate sets.
        correction_legacy = "applied: the approach of earlier versions (%s; replicate sets: %s)",
        # %s: the number of factors.
        correction_k = "%s factor(s)",
        # %s: the data set, then the number of factors.
        correction_view_applied = "%s: %s factor(s)",
        # %s: the data set.
        correction_view_not_applied = "%s: not applied",
        correction_unknown = "applied (approach not recorded)",
        method_ruviii_c = "RUV-III-C",
        method_ruviii = "RUV-III",
        sets_group = ruv[["sets_group"]],
        sets_technical_replicates = ruv[["sets_technical_replicates"]],
        design_groups_all_features = ruv[["design_groups_all_features"]],
        design_technical_replicates = ruv[["design_technical_replicates"]],
        design_supplied_controls = ruv[["design_supplied_controls"]],
        design_first_pass_nulls = ruv[["design_first_pass_nulls"]],
        design_none = ruv[["design_none"]],
        # %s: the time of saving.
        saved_at = "%s UTC",
        status_ready = "ready",
        status_changed = "changed since it was imported: import again",
        status_not_saved = "not saved for integration: save it at that omic's Summary tab",
        status_refused = "cannot be used",
        refused_heading = "Not ready",
        # %s: the omic, then the reason (the check's own message).
        refused_reason = "%s: %s",
        views_heading = "Data sets found",
        col_view = "Data set",
        col_features = "Features",
        col_runs = "Runs",
        col_biological_samples = "Biological samples",
        col_groups = "Groups",
        view_choice_label = "Data sets to combine",
        # %s: the data set, the omic; %d: features, biological samples.
        view_choice_item = "%s (%s; %d features, %d biological samples)",
        import_button = "Import selected data",
        # %d: the omics selected.
        import_need_two_omics = "Select data sets from two or more omics; %d selected.",
        # %d: the data sets selected.
        import_need_two_views = "Select two or more data sets; %d selected.",
        import_need_ready = "Only omics marked ready can be imported.",
        # %s: the reason.
        import_refused = "The import was refused: %s",
        # %d: omics, data sets, runs, biological samples.
        import_done = "Imported %d omics: %d data sets, %d runs, %d biological samples.",
        # %s: the reason.
        store_unavailable = "Integration cannot keep its work in this project: %s",
        # The Sample Matching tab.
        align_heading = "Match samples across omics",
        align_intro = paste("Each row is one biological sample, matched across omics by the sample column each omic",
            "saved. Check that every row is one biological sample, then confirm the table."),
        align_not_imported = "Import two or more omics at the Import tab first.",
        align_key_heading = "How samples are matched",
        # %s: the omic, then its column.
        align_key = "%s: by the column %s",
        align_caution = paste("Samples are matched by name: a sample must carry the same name in every omic. A sample",
            "named differently in two omics appears as two rows."),
        align_rules_heading = "What happens before the model",
        align_provisional = "These rules are provisional: the integration methods may still change them.",
        # %d: all biological samples, those measured in every data set.
        fact_samples = "%d biological samples in all; %d measured in every data set.",
        # %d: samples missing from at least one data set.
        fact_missing_keep_missing = paste("%d biological samples are missing from at least one data set; they are kept,",
            "with that data set's values missing."),
        fact_missing_drop_incomplete = "%d biological samples missing from at least one data set are left out.",
        # %s: the data set; %d: runs, biological samples, samples measured more than once.
        fact_replicates_average = paste("%s: %d runs from %d biological samples; the runs of the %d samples measured more",
            "than once are averaged."),
        # %s: the data set; %d: runs. Every biological sample measured once.
        fact_replicates_single = "%s: %d runs, one per biological sample; nothing to average.",
        # %s: the data set; %d: internal standards.
        fact_standards_drop = "%s: %d internal standards removed.",
        fact_standards_keep = "%s: %d internal standards kept.",
        # %s: the data set; %d: features.
        fact_sparse = "%s: %d features measured in fewer than two samples are left out.",
        # %d: samples in every data set, the threshold.
        few_samples = paste("Only %d biological samples are measured in every data set (fewer than %d). A factor model",
            "needs more samples; its results may be unstable."),
        view_table_heading = "Data sets",
        col_standards_removed = "Internal standards removed",
        col_sparse = "Features left out (fewer than two samples)",
        col_features_used = "Features used",
        sample_table_heading = "Biological samples",
        col_sample = "Biological sample",
        col_group = "Group",
        col_views_present = "Data sets measured",
        cell_missing = "missing",
        confirm_button = "Confirm this sample table",
        # %s: the time, then the file.
        confirmed = "Sample table confirmed at %s and saved to %s.",
        not_confirmed = "Not confirmed yet.",
        # %s: the reason.
        align_refused = "The samples cannot be matched: %s",
        # %s: the data sets.
        refuse_no_shared = "These data sets share no biological sample (%s). Check that each omic names its samples the same way.",
        # %d: the data sets selected.
        refuse_too_few_views = "Two or more data sets are needed; %d selected.",
        # KIT-140: the check's own message is shown under "Technical detail", not in the sentence.
        refuse_changed = "A data set has changed since it was imported. Import the data again at the Import tab.",
        # KIT-140: the refusals in plain words, by the refusal's condition class (`integrationRefusalText`); the check's
        # own message is kept under a collapsed "Technical detail" line.
        technical_detail = "Technical detail",
        refusal_input_invalid = paste("This omic's saved files cannot be used as they are: they are incomplete, or were",
            "changed after they were saved. Save the omic again at its Summary tab, then search again."),
        refusal_input_stale = "A data set has changed since it was imported. Import the data again at the Import tab.",
        refusal_input_duplicate = paste("Two inputs are of the same omic, or two data sets have the same name. Keep one",
            "input per omic."),
        refusal_too_few_omics = "Two or more omics are needed.",
        refusal_too_few_views = "Two or more data sets are needed.",
        refusal_inputs_missing = "Nothing is imported yet. Import two or more omics at the Import tab first.",
        refusal_no_shared = paste("The data sets share no biological sample. Check that each omic names its samples the",
            "same way."),
        refusal_group_conflict = paste("A biological sample is in different groups in different omics. Check the groups",
            "in each omic's design."),
        refusal_view_empty = "A data set has no feature measured in two or more samples, so nothing of it is left for the model.",
        # KIT-153.
        refusal_input_matrix_unavailable = paste("The corrected values were chosen, but at least one data set in the model",
            "has none: its omic did not apply the correction for unwanted variation (RUV). Choose the normalised values, or",
            "leave that data set out."),
        refusal_parameters_invalid = "A model setting is not valid. Check the values entered.",
        refusal_mofa_unavailable = paste("The model cannot be trained on this computer: the MOFA2 software, or the Python",
            "it uses, is not installed."),
        refusal_other = "This step could not be completed; the technical detail below says why.",
        # The Model tab (KIT-140).
        tab_model = " Model",
        step_model = "Model",
        model_heading = "Train the factor model (MOFA2)",
        model_intro = paste("MOFA2 looks for a small number of factors that explain the variation the data sets share.",
            "Check the settings, then train the model."),
        model_not_confirmed = "Confirm the sample table at the Sample Matching tab first.",
        readiness_heading = "Software check",
        readiness_ok = "MOFA2 and the Python it uses are installed: the model can be trained on this computer.",
        readiness_recheck = "Check again",
        settings_heading = "Model settings",
        settings_matching_heading = "Set by the confirmed sample table",
        settings_provisional = paste("Provisional: the integration methods may still change this default. Values",
            "marked so are the current defaults, not final ones."),
        mark_provisional = "provisional default",
        mark_decided = "decided",
        # %s: the label, then the value in words.
        setting_fixed = "%s: %s",
        # %s: the setting's label.
        setting_not_whole = "%s: enter a whole number.",
        # KIT-152. %s: the setting's label.
        setting_not_number = "%s: enter a number.",
        # The settings, by the parameter's name (`param_<name>`, `help_<name>`, `choice_<name>_<value>`). A setting
        # without an entry here shows its name: a parameter the workbook adds appears with no edit to the tab.
        # KIT-153 (INT-23): which values each data set enters the model with; the numbers are the simulation of
        # docs/qa/workflow-execution/kit/kit-151d/ruv-results.md (branch run-mofa-evidence-k151d).
        param_input_matrix = "Values each data set enters the model with",
        choice_input_matrix_uncorrected = "normalised values, before correction for unwanted variation (RUV)",
        choice_input_matrix_ruv_corrected = "values corrected for unwanted variation (RUV)",
        help_input_matrix = paste("The normalised values keep biological differences between samples of the same group,",
            "which the correction can remove: in a simulation, the default correction removed a biological pattern shared",
            "by two omics that varied within the groups (found in 0 of 30 data sets after correction, 30 of 30 before).",
            "The corrected values protect comparisons between groups when unwanted variation is correlated with the",
            "groups: with 20 samples and unwanted variation correlated 0.6 with the groups, the factor for the groups was",
            "found in 6 of 30 data sets before correction and 30 of 30 after. The corrected values can be chosen only when",
            "every data set in the model has them."),
        param_technical_replicates = "Runs of the same biological sample",
        choice_technical_replicates_average = "averaged",
        param_internal_standards = "Internal standards",
        choice_internal_standards_drop = "removed",
        choice_internal_standards_keep = "kept",
        param_missing_samples = "Samples missing from a data set",
        choice_missing_samples_keep_missing = "kept, with that data set's values missing",
        choice_missing_samples_drop_incomplete = "left out",
        param_feature_selection = "Features used from each data set",
        choice_feature_selection_variance_above_trend = paste("the features that vary most for their abundance (variance",
            "above the trend of variance with mean abundance)"),
        choice_feature_selection_top_variance = "the most variable features (variance alone)",
        choice_feature_selection_none = "all features",
        help_feature_selection = paste("Selecting the same number of features from each data set keeps a large data set",
            "from outweighing the small ones. Less abundant features are often noisier; ranking by variance above the",
            "trend keeps that noise from filling the selection."),
        param_features_per_view = "Number of features kept per data set",
        help_features_per_view = paste("A data set with more features keeps this many; a smaller data set keeps all of",
            "its features. Not used when all features are used."),
        param_min_observed_fraction = "Share of biological samples a feature must be measured in",
        help_min_observed_fraction = paste("Features measured in fewer biological samples are left out before the",
            "selection (0.5: at least half of them, and never fewer than two). Not used when all features are used."),
        param_scaling = "Scaling of each feature",
        choice_scaling_none = "none: the values as saved",
        choice_scaling_feature_zscore = "each feature scaled to mean 0 and standard deviation 1",
        # KIT-151 (INT-24). KIT-239: the simulated pattern's size and the condition in numbers (KIT-151 d, z-summary.tsv:
        # F2 |r| >= 0.9 in 26 against 12 of 30 with missing values, 29 against 26 without, 30 of 30 at 60 samples). The
        # last sentence is AUDIT-KIT-239-impl I02's (KIT-246 option (a)): a draft pending the owner's answer.
        help_scaling = paste("Each selected feature is scaled to mean 0 and standard deviation 1, after the features are",
            "chosen (the choice uses the values as saved). In simulations with 20 samples in which about a fifth of the",
            "values were missing, mostly low ones, this recovered closely (a correlation of at least 0.9 with the true",
            "pattern) a pattern carried by a tenth of the features of the two smallest data sets, each changed by about as",
            "much as its measurement noise, in 26 rather than 12 of 30 simulated data sets. Without missing values the",
            "difference could not be told from none (29 rather than 26 of 30); with 60 samples the pattern was recovered",
            "closely in 30 of 30 either way."),
        # KIT-151 (audit SCI-12): what scale_views does.
        param_scale_views = "Divide each data set by its overall standard deviation",
        help_scale_views = paste("Each data set is divided by one number, the standard deviation of all its values, so that",
            "data sets measured on different scales become comparable. It does not balance data sets with different",
            "numbers of features; the feature selection does that."),
        param_center_groups = "Centre each feature on its mean",
        help_center_groups = paste("Required unless each feature is scaled to mean 0 and standard deviation 1, which",
            "centres it already."),
        param_num_factors = "Number of factors",
        param_convergence_mode = "How closely training must settle before it stops",
        choice_convergence_mode_fast = "fast",
        choice_convergence_mode_medium = "medium",
        choice_convergence_mode_slow = "slow (the closest)",
        param_maxiter = "Maximum number of training rounds",
        help_maxiter = paste("Training stops by its convergence rule or at this limit, whichever comes first; the result",
            "says which."),
        param_seed = "Random seed",
        # KIT-151 (INT-17): what the seed promises, and what it does not.
        help_seed = paste("On the same computer and software, the same data, settings and seed give the same model. On",
            "another computer the model usually agrees to about 12 decimal places, but a factor near the limit of what",
            "the data can show can come out differently."),
        param_training_seeds = "Number of models trained (one per seed)",
        help_training_seeds = paste("The models are trained from the seed above and the next whole numbers. The results",
            "shown are the first model's; how well each factor agrees across the models is shown with the results. More",
            "models take longer to train."),
        param_precision = "Numeric precision",
        choice_precision_double = "double: the most exact (the chosen setting)",
        choice_precision_single = "single: less memory for very large studies",
        help_precision = paste("Double precision is the chosen setting: on the same computer and software it gives the",
            "same model every time, and on another computer it usually agrees to about 12 decimal places. Single",
            "precision uses less memory for very large studies, but the model differs more from one computer to another."),
        run_button = "Train the model",
        run_working = "Training the model. This can take several minutes.",
        # %s: the reason.
        run_refused = "The model was not trained. %s",
        run_none = "No model trained yet.",
        # %s: the time.
        run_done = "Model trained at %s.",
        run_reused = "These data and settings were trained before; the saved model is shown.",
        run_settings_changed = "The settings have changed since this model was trained. Train again to apply them.",
        # %d: models trained since the import.
        run_count = "Models trained since the import: %d. Each one is kept in the project; the latest is shown.",
        result_heading = "The trained model",
        # %d: rounds, then the limit. KIT-151 (INT-15): said only when training stopped by its convergence rule.
        result_converged = "Training settled after %d rounds: it stopped by its convergence rule (the limit was %d rounds).",
        # %d: the limit.
        result_not_converged = paste("Training did not settle: it stopped at the limit of %d rounds before its",
            "convergence rule was met. Allow more rounds, then train again."),
        # %d: factors kept, factors asked for.
        result_factors = "Factors in the model: %d (asked for: %d).",
        # %d: biological samples, those measured in every data set.
        result_samples = "Biological samples in the model: %d; measured in every data set: %d.",
        # KIT-151 (audit SCI-12): variance explained refers to the features given to the model.
        result_variance_heading = "Variance explained in each data set's selected features (all factors together)",
        col_variance = "Variance explained in the selected features (%)",
        col_model_features = "Features in the model",
        # %s: the checksum.
        result_saved = "The model file is saved in the project; its SHA-256 checksum is %s.",
        # KIT-063 item A: correcting the sample match with an edited table (Sample Matching tab).
        match_heading = "Correct the match (optional)",
        match_intro = paste("If two omics name the same biological sample differently, download the sample table, change",
            "the biological_sample column so that each run names its biological sample, and upload the edited table. Runs",
            "you leave unchanged keep the name their omic saved. The new match is shown above before you confirm it."),
        match_download = "Download the sample table",
        match_upload = "Upload the edited table (tab-separated)",
        match_clear = "Use the names as saved",
        match_active = "Your uploaded table is in use.",
        # %d: runs.
        match_changed = "Your table gives %d runs a different biological sample.",
        match_unchanged = "Your table gives no run a different biological sample; samples are matched as saved.",
        # %d: rows.
        match_unknown_omics = "%d rows name an omic that is not imported; they are ignored.",
        match_unknown_runs = "%d rows name a run that is not in the data; they are ignored.",
        match_unknown_names = "%d rows name a sample name that is not in the data; they are ignored.",
        match_rows_empty = "%d rows name neither a run nor a sample name; they are ignored.",
        # %s: the reason.
        match_refused = "The uploaded table cannot be used. %s",
        refusal_sample_match_invalid = paste("Check that the table has the columns omic, run and biological_sample, that",
            "every run listed has a biological sample, and that each run appears only once."),
        # %d: runs (Model tab, the settings the confirmed table fixed).
        setting_sample_match = "Sample match: your uploaded table, which gives %d runs a different biological sample",
        # KIT-063: the Results tab.
        tab_results = " Results",
        step_results = "Results",
        results_heading = "Results of the factor model",
        results_intro = paste("What the trained model found: how much of the variation of each data set's selected",
            "features each factor explains,",
            "how the factors relate to the groups, which features weigh most on each factor, and how the factors relate to",
            "each other."),
        results_not_trained = "Train a model at the Model tab first.",
        # %s: the time the model was trained.
        results_model = "Results of the model trained at %s.",
        # %s: the reason.
        results_refused = "The results cannot be shown. %s",
        results_provisional = paste("How many features are listed, which factors are shown, how the factors are related to",
            "the groups and how the figures look are provisional: the integration methods may still change them."),
        advisories_heading = "Notes from the model",
        # %s: the factors (Factor1, Factor2).
        advisory_overall_level = paste("%s rises and falls with each sample's overall level in at least one data set.",
            "Such a factor often reflects differences in overall level between samples, for example from normalisation,",
            "rather than biology. Read it with care."),
        advisory_other = "MOFA2 gave a warning while training; its text is under Technical detail.",
        variance_heading = "Variance explained",
        variance_intro = paste("The percentage of the variation of each data set's selected features (the features given",
            "to the model, not all features measured) that each factor explains. A factor that explains variation in data",
            "sets of different omics is shared across those omics; one that explains variation in only one data set is",
            "specific to it. Two data sets of one omic (for example positive and negative mode) measure related",
            "molecules, so a factor active in both is not by itself shared across omics."),
        variance_total_heading = "All factors together, per data set (its selected features)",
        # KIT-151 (INT-14): each factor's agreement across the training seeds.
        stability_heading = "How stable each factor is across training seeds",
        # %d: models; %s: the seeds; %d: the seed of the results shown.
        stability_intro = paste("The model was trained %d times, each from a different random seed (%s); the results",
            "shown are those of seed %d. For each factor: its agreement with the most similar factor of each other model,",
            "that is the correlation of the two factors' values across the biological samples (sign ignored), each factor",
            "paired with one factor only, the most similar pair first. Near 1, the factor is found whatever the seed; well",
            "below 1, it depends on the seed, so read it with care."),
        stability_single = paste("The model was trained from one seed only, so how stable the factors are across seeds",
            "was not measured."),
        col_agreement_lowest = "Lowest agreement across seeds",
        col_agreement_median = "Median agreement across seeds",
        col_factor = "Factor",
        association_heading = "Factors and the groups",
        # %s: the adjustment's name.
        association_intro = paste("For each factor: the share of its variation that lies between the groups (one-way",
            "analysis of variance of the factor's values on the group), the test's p-value, and the p-value adjusted for",
            "all the factors of the model, whichever are shown (%s)."),
        padjust_BH = "Benjamini-Hochberg",
        padjust_holm = "Holm",
        padjust_bonferroni = "Bonferroni",
        padjust_none = "not adjusted",
        claim_heading = "How to read this",
        claim_descriptive = paste("The factors were found without using the groups. How they relate to the groups describes",
            "the model; it is not a test of the experiment."),
        # %s: the data sets.
        claim_group_corrected = paste("In %s, unwanted variation was removed with the groups as replicate sets, which keeps",
            "the differences between groups by construction. Factors can then separate the groups more than the data",
            "alone would, so no p-value is given for the groups."),
        # KIT-153: the model used the normalised values before correction.
        claim_uncorrected = paste("The data sets entered the model before correction for unwanted variation (RUV). Unwanted",
            "variation that is correlated with the groups has not been removed, so it can make a factor follow the groups;",
            "read the p-values with that in mind."),
        # %s: the omic, then its correction in words.
        claim_correction = "%s: correction for unwanted variation %s",
        # KIT-239 (INT-28): the limitations of the model, on the Results tab after the notes from the model and in the
        # integration report after its Methods draft (`integrationLimitationsText`). The `limits_fixed_*` sentences hold
        # for every run; the others are filled from the model's sample table or the record. Sources of every number:
        # docs/qa/workflow-execution/kit/kit-239/wording.md.
        limits_heading = "Limitations of this model",
        limits_lead = "Read these before drawing conclusions from the factors.",
        # %d: the biological samples in the model, then those measured in every data set.
        limits_samples = "This model has %d biological samples; %d of them are measured in every data set.",
        limits_fixed_samples = paste("How small a pattern the model can find depends on the number of samples. In",
            "simulations made for this software with 20 samples (two groups of 10), a pattern carried only by the two",
            "smallest of four data sets, in a tenth of their features, each changed by about two-thirds of its measurement",
            "noise, was recovered closely (a correlation of at least 0.9 with the true pattern) in at most 2 of 10",
            "simulated data sets, whichever setting was tried; the settings tried did not include this software's current",
            "default feature selection. With 60 samples and the slow training this software uses, the same pattern was",
            "recovered in 7 to 10 of 10, and with this software's default settings a pattern as large as the noise was",
            "recovered in 26 to 29 of 30 simulated data sets with 20 samples. MOFA2's documentation says that factor models",
            "are useful only with large numbers of samples, at least more than 15. A pattern this model does not show may",
            "still be in the data: treat the factors as descriptive."),
        # %d: the biological samples missing from at least one data set, then the biological samples in the model.
        limits_missing_samples = paste("Biological samples missing from at least one data set: %d of %d; the model",
            "leaves those values out."),
        limits_fixed_missing = paste("Where a value is missing, the model leaves it out (MOFA2 does not impute it), which",
            "assumes that whether a value is missing does not depend on the value itself (missing at random). Values",
            "missing because they were too low to measure are not modelled as such, so the model's estimates for",
            "low-abundance features can be biased. MOFA2 fills missing values with each feature's mean only to compute",
            "where training starts."),
        # %s: the data set; %d: its biological samples with two or more runs, its biological samples, the most runs of
        # one biological sample.
        limits_averaged_item = "%s: the runs of %d of its %d biological samples were averaged (at most %d runs for one sample).",
        limits_averaged_none = "No biological sample has more than one run in any data set, so no runs were averaged.",
        # Said only when runs were averaged.
        limits_fixed_averaged = paste("Each average uses only the runs in which the value was measured. When a value is",
            "missing from some of a sample's runs because it was too low to measure, the average of the other runs is too",
            "high; this affects low-abundance features most."),
        limits_fixed_model = paste("The model treats all samples as one group (the groups are not given to it, and MOFA2's",
            "option for several groups of samples is not used), assumes normally distributed (Gaussian) variation in every",
            "data set, and takes one input per omic (one input can hold several data sets, such as positive and negative",
            "ion mode). The integration fills in no missing value itself; values that an omic's own workflow estimated",
            "before saving are treated as measured."),
        limits_fixed_check = paste("Given the settings of the analysis script it was adapted from, this software trains the",
            "identical model; this was checked on simulated data only, because the script's own data were not available.",
            "The check is about the software, not about any biological result, and this software's default settings",
            "differ from the script's. The defaults for the input matrix, the feature selection, the scaling and the",
            "convergence rule were chosen from simulations of one design, not from real data; the number of factors comes",
            "from the analysis script."),
        # The Results tab (it comes before the Enrichment tab).
        limits_fixed_services = paste("Enrichment, at the Enrichment tab, needs external services: STRING, queried over the",
            "internet (by default its version 12.5), for proteins and genes, and tables downloaded from KEGG or Reactome",
            "for metabolites and lipids. Its results depend on the release used, which is recorded with each result."),
        # The report. %s: the databases and releases the analysis's enrichment used.
        limits_services_report = paste("The enrichment in this analysis used %s; its results depend on the release of",
            "each external database named."),
        limits_services_none = "No enrichment was run in this analysis.",
        # The report, when the record's integration block cannot be read.
        limits_not_computed = "This analysis's own numbers could not be read from its record.",
        # The Model tab, after training.
        model_limits_pointer = "Read the limitations of this model on the Results tab before drawing conclusions from it.",
        col_samples_tested = "Samples",
        col_groups_n = "Groups",
        col_group_r2 = "Variation between groups (%)",
        col_f = "F",
        col_df = "Degrees of freedom",
        col_p = "p-value",
        col_p_adjusted = "Adjusted p-value",
        # KIT-151 (plan item E): the Kruskal-Wallis columns.
        col_epsilon2 = "Effect size, epsilon-squared (%)",
        col_h = "Kruskal-Wallis H",
        # %s: the adjustment's name.
        association_intro_kruskal_wallis = paste("For each factor: the effect size epsilon-squared (the Kruskal-Wallis",
            "statistic H divided by the number of samples minus one, as a percentage: how much of the ordering of the",
            "factor's values lies between the groups), the Kruskal-Wallis test's p-value, and the p-value adjusted for all",
            "the factors of the model, whichever are shown (%s)."),
        cell_not_tested = "not given (see How to read this)",
        cell_cannot_test = "cannot be computed (fewer than two groups or too few samples)",
        factor_values_heading = "Factor values per biological sample",
        weights_heading = "Features with the largest weights",
        # %d: the number of features.
        weights_intro = paste("The %d features with the largest weights, in absolute value, on the chosen factor in the",
            "chosen data set. A feature's weight says how strongly it follows the factor; its sign says in which",
            "direction. The scaled weight divides by the largest weight of that factor in that data set, so 1 or -1 marks",
            "the strongest feature.",
            # KIT-152 (audit SCI-16): the sign and the list's standing.
            "A factor's sign is arbitrary, so a direction has meaning only against the other features of the same factor.",
            "The list is exploratory: its order was not tested for stability."),
        weights_factor = "Factor",
        weights_view = "Data set",
        col_rank = "Rank",
        col_feature = "Feature",
        col_feature_id = "Identifier",
        col_weight = "Weight",
        col_weight_scaled = "Scaled weight",
        correlation_heading = "How the factors relate to each other",
        correlation_intro = paste("The correlation (Pearson) of the factors' values across the biological samples. Each",
            "factor is meant to capture its own pattern; a correlation far from 0 means two factors share part of theirs."),
        plot_dataset = "Data set",
        plot_factor = "Factor",
        plot_variance = "Variance explained in the selected features (%)",
        plot_group = "Group",
        plot_factor_value = "Factor value",
        plot_weight = "Weight",
        plot_weight_scaled = "Scaled weight",
        plot_correlation = "Correlation",
        # %s: the factor, then the data set.
        plot_weights_title = "%s, %s",
        save_heading = "Save the results",
        save_button = "Save tables and figures to the project",
        # %s: the time; %d: tables; %s: their folder; %d: figures; %s: their folder.
        save_done = "Saved at %s: %d tables to %s and %d figures to %s.",
        # %s: the reason.
        save_refused = "The results were not saved. %s",
        # KIT-152 correction. %d: the factors in the model.
        factors_shown_reset = paste("Some factors picked under \"Factors shown\" are not in this model, which has %d factors,",
            "so all factors are shown."),
        # KIT-152: the Results tab's settings (`integrationResultsParameterSpec()`), by the parameter's name as above.
        results_settings_heading = "Settings of the results",
        param_factors_shown = "Factors shown",
        choice_factors_shown_all = "all factors",
        param_top_features = "Number of features listed per factor and data set",
        param_weight_scale = "Weights shown",
        choice_weight_scale_scaled = "divided by the largest weight, so the strongest feature is 1 or -1",
        choice_weight_scale_raw = "as estimated",
        param_design_covariate = "Groups the factors are compared with",
        choice_design_covariate_group = "the group of each biological sample",
        param_association_test = "Test of each factor against the groups",
        choice_association_test_anova = "one-way analysis of variance",
        # KIT-151 (plan item E).
        choice_association_test_kruskal_wallis = paste("Kruskal-Wallis rank test (does not assume the factor's values are",
            "normally distributed)"),
        help_association_test = paste("Choose the test before looking at the results. Analysis of variance is the default;",
            "the Kruskal-Wallis test compares the ranks of the factor's values instead."),
        param_p_adjust_method = "Adjustment of the p-values",
        choice_p_adjust_method_BH = "Benjamini-Hochberg",
        choice_p_adjust_method_holm = "Holm",
        choice_p_adjust_method_bonferroni = "Bonferroni",
        choice_p_adjust_method_none = "not adjusted",
        param_figure_format = "File type of the saved figures",
        choice_figure_format_png = "PNG image",
        choice_figure_format_pdf = "PDF",
        param_figure_width = "Width of the saved figures (inches)",
        param_figure_height = "Height of the saved figures (inches)",
        param_figure_dpi = "Resolution of the saved figures (dots per inch)",
        param_feature_values = "Values shown for the top features",
        choice_feature_values_model_input = "the values given to the model",
        help_feature_values = paste("The values each data set gave the model: the runs of a biological sample averaged, only",
            "the features the model used, and each feature scaled if that scaling was chosen; before MOFA2 centres and",
            "scales the data sets itself."),
        param_heatmap_scaling = "Heatmap colours",
        choice_heatmap_scaling_centre = "each feature centred on its mean",
        choice_heatmap_scaling_centre_scale = "each feature centred on its mean and scaled to standard deviation 1",
        help_heatmap_scaling = "For display only: the values shown in the scatter plots and saved in the table do not change.",
        param_umap = "Draw a UMAP map of the factors",
        param_umap_n_neighbors = "UMAP: number of neighbours",
        help_umap_n_neighbors = "At least 2, and fewer than the number of biological samples.",
        param_tsne = "Draw a t-SNE map of the factors",
        param_tsne_perplexity = "t-SNE: perplexity",
        help_tsne_perplexity = paste("Three times the perplexity must be less than the number of biological samples minus 1",
            "(below 6.33 with 20 samples)."),
        param_map_seed = "Random seed for the maps",
        help_map_seed = paste("Used by the UMAP map: on the same computer, the same factors, settings and seed give the",
            "same map. The t-SNE map starts from the factors' first two principal components and does not depend on the seed."),
        # KIT-152: the factor pair plot.
        pairs_heading = "Factors against each other",
        pairs_intro = paste("Each panel shows two factors; each point is one biological sample, coloured by its group. A",
            "factor's sign is arbitrary: a factor and its mirror image describe the same pattern."),
        plot_pairs_one = "Choose two or more factors in the settings to compare them in pairs.",
        # KIT-152: the top features' values against the factor (scatter) and as a heatmap.
        scatter_heading = "Values of the top features against the factor",
        # %d: the number of features.
        scatter_intro = paste("For the chosen factor and data set, each panel is one of the %d features with the largest",
            "weights in absolute value: its value in each biological sample against that sample's factor value. The",
            "features were chosen from this model, so a pattern is expected; the ranking is exploratory and was not tested",
            "for stability. A factor's sign is arbitrary: a feature with a negative weight falls as the factor rises."),
        # %s: the values in words.
        scatter_values = "Values shown: %s. Missing values are left out.",
        values_log2 = "the values given to the model, on the log2 scale the omic saved",
        values_log2_cpm = "the values given to the model, in log2 counts per million as the omic saved them",
        values_feature_zscore = "the values given to the model, each feature scaled to mean 0 and standard deviation 1",
        values_unknown = "the values given to the model, on the scale the omic saved",
        plot_values_log2 = "Value given to the model (log2)",
        plot_values_log2_cpm = "Value given to the model (log2 counts per million)",
        plot_values_feature_zscore = "Value given to the model (scaled per feature)",
        plot_values_unknown = "Value given to the model",
        # %s: the factor.
        plot_scatter_x = "%s value",
        # %s: the factor, then the data set.
        plot_scatter_title = "%s, %s: features with the largest weights",
        plot_scatter_caption = "Panels in order of the weight's absolute value; missing values left out.",
        heatmap_heading = "Heatmap of the top features",
        heatmap_intro_centre = paste("The same features and values, each feature centred on its mean across the biological",
            "samples, for display only; samples in order of the factor's value, features in order of rank."),
        heatmap_intro_centre_scale = paste("The same features and values, each feature centred on its mean and scaled to",
            "standard deviation 1 across the biological samples, for display only; samples in order of the factor's value,",
            "features in order of rank."),
        plot_heatmap_centre = "Value given to the model,\ncentred per feature",
        plot_heatmap_centre_scale = "Value given to the model,\ncentred and scaled per feature",
        # %s: the factor.
        plot_heatmap_x = "Biological samples, in order of %s",
        plot_heatmap_caption = "For display only; the saved table holds the values given to the model.",
        # KIT-152: the UMAP and t-SNE maps of the factors.
        maps_heading = "Maps of the factors (UMAP, t-SNE)",
        maps_off = paste("No map is drawn. Tick \"Draw a UMAP map of the factors\" or \"Draw a t-SNE map of the factors\"",
            "in the settings above to draw one."),
        # %d: the biological samples on the map.
        maps_caution = paste("These maps are a picture of the factors, not evidence of groups. With %d biological samples",
            "they are unstable, and the distances between points and the sizes of clusters on them have no meaning. On",
            "another computer the same map can come out slightly different."),
        # KIT-152 correction: what each map depends on, shown with that map only.
        maps_caution_umap = "The UMAP map changes with the random seed and with the number of neighbours.",
        maps_caution_tsne = paste("The t-SNE map starts from the factors' first two principal components, so it does not",
            "depend on the random seed; it changes with the perplexity."),
        # %s: the data sets.
        maps_caution_group_corrected = paste("In %s, unwanted variation was removed with the groups as replicate sets, so",
            "the groups separate on these maps by construction."),
        map_umap = "UMAP",
        map_tsne = "t-SNE",
        # %s: the map; %d: 1 or 2.
        plot_map_axis = "%s %d",
        # %s: the map; %d: factors, biological samples; %s: the settings.
        plot_map_title = "%s of %d factors, %d biological samples (%s)",
        # %d: neighbours, seed.
        map_settings_umap = "%d neighbours, seed %d",
        # %s: perplexity; %d: seed.
        map_settings_tsne = "perplexity %s, seed %d",
        # %d: the factors in the model.
        map_refused_factors = "The map was not drawn: a map needs two or more factors, and the model has %d.",
        # %d: neighbours, biological samples, the largest number allowed.
        map_refused_neighbors = paste("The UMAP map was not drawn: the number of neighbours is %d, but with %d biological",
            "samples it must be at least 2 and at most %d."),
        # %s: perplexity; %d: biological samples; %s: the limit.
        map_refused_perplexity = paste("The t-SNE map was not drawn: the perplexity is %s, but with %d biological samples it",
            "must be above 0 and below %s (three times the perplexity must be less than the number of samples minus 1)."),
        # %s: the map, then the software.
        map_unavailable = "The %s map cannot be drawn on this computer: the software %s is not installed.",
        # KIT-142: the Record tab (its own words are `integrationRecordTexts()`).
        tab_record = " Record",
        step_record = "Record",
        # KIT-141: the Enrichment tab.
        tab_enrichment = " Enrichment",
        step_enrichment = "Enrichment",
        enrich_heading = "Enrichment of a factor",
        enrich_intro = paste("Find the biological processes, pathways or compound classes that follow a factor. Every",
            "measured feature of the chosen data set, not only those given to the model, is ranked by how closely its",
            "values follow the factor across the biological samples (a correlation from -1 to 1). Proteins and genes are",
            "tested with STRING; metabolites and lipids with KEGG or Reactome pathways."),
        enrich_optional = "This step is optional: the Record tab is open without it.",
        enrich_not_trained = "Train a model at the Model tab first.",
        enrich_factor = "Factor",
        enrich_view = "Data set",
        enrich_database = "Database",
        enrich_database_string = "STRING (proteins and genes)",
        enrich_database_kegg = "KEGG pathways (metabolites and lipids)",
        enrich_database_reactome = "Reactome pathways (metabolites and lipids)",
        enrich_id_column = "Column of the feature table holding the compound identifiers",
        enrich_id_column_kegg = "KEGG compound identifiers (such as C00031) are needed.",
        enrich_id_column_reactome = "ChEBI identifiers (such as CHEBI:17234 or 17234) are needed.",
        enrich_settings_heading = "Settings",
        enrich_correlation_method = "Correlation",
        enrich_correlation_pearson = "Pearson (straight-line relation)",
        enrich_correlation_spearman = "Spearman (rank order)",
        enrich_min_samples = "Smallest number of biological samples a feature needs to be ranked",
        enrich_string_species = "STRING species number (9606 is human, 10090 mouse)",
        enrich_string_fdr = "STRING: show terms up to this false discovery rate",
        enrich_kegg_organism = "KEGG organism code (hsa is human, mmu mouse)",
        enrich_reactome_organism = "Reactome organism (such as Homo sapiens)",
        enrich_min_set = "Smallest pathway tested (measured compounds in it)",
        enrich_max_set = "Largest pathway tested (measured compounds in it)",
        enrich_seed = "Seed of the pre-ranked and permutation tests (the correlation-adjusted tests use none)",
        enrich_key = "Your STRING key",
        enrich_key_note = paste("STRING needs a personal key (free from string-db.org). It is used only for this query",
            "and is never saved, recorded or written to a file. Enter it each time."),
        enrich_kegg_terms = paste("KEGG's tables are free for academic use only; other users need a licence from KEGG",
            "(Pathway Solutions, www.pathway.jp) before downloading them."),
        enrich_download = "Download the pathway tables",
        enrich_downloading = "Downloading the pathway tables.",
        # %s: database, release; %d: rows.
        enrich_downloaded = "%s tables downloaded: %s, %d rows.",
        enrich_run = "Run the enrichment",
        enrich_running = "Running the enrichment.",
        # %s: the reason.
        enrich_refused = "The enrichment was not run. %s",
        # %s: factor, data set, database, release.
        enrich_done = "%s in %s against %s (%s).",
        # %d: measured, ranked, too few samples, no spread, internal standards.
        enrich_counts = paste("Features: %d measured, %d ranked, %d left out for too few samples, %d left out for no",
            "spread, %d internal standards left out."),
        # %d: matched, not matched, several on one identifier.
        enrich_mapping = "Matched to the database: %d; not matched: %d; set aside as a second feature on one identifier: %d.",
        # %s: the feature, its data set; %d: +1 or -1.
        enrich_direction = paste("Direction: the factor is oriented so that its largest weight (%s in %s) is positive",
            "(multiplied by %d). Terms at the top of the list follow the factor; terms at the bottom go against it."),
        # %d: terms or pathways.
        enrich_tests = paste("%d terms or pathways in the table. P-values are adjusted within this factor, data set and",
            "database only, not across factors, data sets or databases."),
        enrich_none = "No term or pathway passed the threshold.",
        enrich_terms_heading = "Enriched terms and pathways",
        # KIT-183 (VAL-51, VAL-53; proposed wording, VAL-41: the owner edits before release). The numbers in brackets
        # are the simulation's measured rates (docs/qa/workflow-execution/kit/kit-183/results.md), written in as constants.
        # VAL-115 (4): the permutation notes take the number of reorderings chosen (%s; their literal % doubled), and every
        # offered test's note gives the caution about paired, repeated or batched samples.
        # KIT-191 (VAL-85 (2), VAL-97): four tests, each note per kind of data set (metabolites, shown also for lipids;
        # genes and proteins, through Gene Ontology), its figures from KIT-183's counted simulation (run 45225854502)
        # as `INTEGRATION_ENRICHMENT_TEST_EVIDENCE` holds them; regenerated by kit-191/evidence-k191.R (U14).
        enrich_test = "Test",
        enrich_test_camera = "Correlation-adjusted test (camera), the default",
        enrich_test_camera_floor = "Correlation-adjusted test, negative correlations counted as zero (camera)",
        enrich_test_permutation = "Permutation of the factor's values (PCGSE)",
        enrich_test_gsea = "Pre-ranked test (GSEA)",
        enrich_test_camera_pr = "Correlation-adjusted test with an assumed correlation (cameraPR)",
        enrich_test_camera_note_genes_proteins = paste("Each feature's values are regressed on the factor across the",
            "biological samples. A term is called when its features follow the factor more closely, in one direction, than",
            "the other features, allowing for how strongly the term's features vary together, as measured in these data",
            "(limma camera; Wu and Smyth 2012). Only features measured in every biological sample of this data set can be",
            "used. The factor was learned from these same data, which the test does not allow for, so read the result as a",
            "list of terms worth a closer look, not as proof. In a check on simulated data with 12 biological samples, terms",
            "with no link to the factor were called at the 5 % level in 0.4 to 6.6 % of cases, depending on the kind of",
            "term; terms made only of the most variable features, the ones the model selects, were called in up to 10.1 %",
            "even when the factor was replaced by values unrelated to the data, and in up to 12.8 % on factors weaker than",
            "the main effect. In the same check, terms of 10 genes or proteins in which half the features truly followed",
            "the main factor were found in 60 to 62 % of cases, and terms of 30 such features in 99 to 100 % (in result",
            "tables of about 120 to 150 terms, about 1 in 6 truly linked).",
            "It treats the biological samples as independent, which samples taken in pairs, repeatedly from one subject,",
            "or in batches are not."),
        enrich_test_camera_note_metabolites = paste("Each feature's values are regressed on the factor across the biological",
            "samples. A term is called when its features follow the factor more closely, in one direction, than the other",
            "features, allowing for how strongly the term's features vary together, as measured in these data (limma",
            "camera; Wu and Smyth 2012). Only features measured in every biological sample of this data set can be used.",
            "The factor was learned from these same data, which the test does not allow for, so read the result as a list",
            "of terms worth a closer look, not as proof. In a check on simulated metabolite data with 12 biological",
            "samples, terms with no link to the factor were called at the 5 % level in 0.4 to 4.5 % of cases, depending on",
            "the kind of term; terms made only of the most variable features, the ones the model selects, were called in",
            "up to 4.5 % even when the factor was replaced by values unrelated to the data, and in up to 6.7 % on factors",
            "weaker than the main effect. In the same check, terms of 5 metabolites, 2 of them truly following the main",
            "factor, were found in 2.0 to 2.2 % of cases (in result tables of about 90 terms, about 1 in 9 truly linked).",
            "It treats the biological samples as independent, which samples taken in pairs, repeatedly from one subject,",
            "or in batches are not."),
        enrich_test_camera_floor_note_genes_proteins = paste("The correlation-adjusted test, with any negative correlation",
            "measured between a term's features counted as zero (limma's own default). This never makes a p-value smaller,",
            "so it never calls a term that the default test does not call. In a check on simulated data with 12 biological",
            "samples, terms with no link were called at the 5 % level in 1.1 to 5.1 % of cases against values unrelated to",
            "the data, depending on the kind of term, and terms made only of the most variable features in up to 5.5 % (the",
            "default test: up to 10.1 %); terms of 10 genes or proteins in which half the features truly followed the main",
            "factor were found in 47 to 50 % of cases (the default test: 60 to 62 %; 30 features: 99 %) (in result tables",
            "of about 120 to 150 terms, about 1 in 6 truly linked). With 6 biological samples it found at most 12 % of such",
            "terms of 10 or 30 features. On factors weaker than the main effect, terms of features the model used were",
            "called in up to 12.8 %, as with the default test. The factor was learned from these same data, which the test",
            "does not allow for.",
            "It treats the biological samples as independent, which samples taken in pairs, repeatedly from one subject,",
            "or in batches are not."),
        enrich_test_camera_floor_note_metabolites = paste("The correlation-adjusted test, with any negative correlation",
            "measured between a term's features counted as zero (limma's own default). This never makes a p-value smaller,",
            "so it never calls a term that the default test does not call. In a check on simulated metabolite data with 12",
            "biological samples, terms with no link to the factor were called at the 5 % level in 1.7 to 3.6 % of cases",
            "against values unrelated to the data, depending on the kind of term (the default test: 3.4 to 4.5 %); terms of",
            "5 metabolites, 2 of them truly following the main factor, were found in 0.1 to 0.2 % of cases (the default",
            "test: 2.0 to 2.2 %) (in result tables of about 90 terms, about 1 in 9 truly linked). On factors weaker than the",
            "main effect, terms of metabolites the model used were called in up to 6.7 %, as with the default test. No lipid",
            "data were simulated. The factor was learned from these same data, which the test does not allow for.",
            "It treats the biological samples as independent, which samples taken in pairs, repeatedly from one subject,",
            "or in batches are not."),
        enrich_test_permutation_note_genes_proteins = paste("Each feature's correlation with the factor is compared with",
            "what %s random reorderings of the factor's values across the biological samples give, so features that",
            "vary together keep varying together (the permutation test of PCGSE, Frost et al. 2015). The reorderings stand",
            "for a factor that follows no feature; this factor was learned from these same data and follows the features",
            "it was learned from, which no reordering reproduces, so the test does not allow for that. It treats the",
            "biological samples as interchangeable, which samples taken in pairs, repeatedly from one subject, or in",
            "batches are not. In a check on simulated data with 12 biological samples, terms with no link were called at",
            "the 5 %% level in 4.8 to 5.3 %% of cases against values unrelated to the data, whatever the kind of term (the",
            "most variable features: 4.9 to 5.3 %%), and in 0.4 to 4.4 %% on the main factor; terms of 10 genes or proteins",
            "in which half the features truly followed the main factor were found in 39 to 52 %% of cases (the default",
            "test: 60 to 62 %%), of 30 in 80 to 96 %% (in result tables of about 120 to 150 terms, about 1 in 6 truly",
            "linked). With 6 biological samples it found 0.2 to 3.2 %% of such terms (10 to 100 features). On factors",
            "weaker than the main effect, terms of features the model used were called in up to 12.7 %%."),
        enrich_test_permutation_note_metabolites = paste("Each feature's correlation with the factor is compared with what",
            "%s random reorderings of the factor's values across the biological samples give, so features that vary",
            "together keep varying together (the permutation test of PCGSE, Frost et al. 2015). The reorderings stand for a",
            "factor that follows no feature; this factor was learned from these same data and follows the features it was",
            "learned from, which no reordering reproduces, so the test does not allow for that. It treats the biological",
            "samples as interchangeable, which samples taken in pairs, repeatedly from one subject, or in batches are not.",
            "In a check on simulated metabolite data with 12 biological samples, terms with no link to the factor were",
            "called at the 5 %% level in 4.2 to 5.2 %% of cases against values unrelated to the data, whatever the kind of",
            "term, and in 0.6 to 4.7 %% on the main factor; terms of 5 metabolites, 2 of them truly following the main",
            "factor, were found in 1.0 %% of cases (the default test: 2.0 to 2.2 %%) (in result tables of about 90 terms,",
            "about 1 in 9 truly linked). On factors weaker than the main effect, terms of metabolites the model used were",
            "called in up to 8.6 %% (the default test: 6.7 %%). No lipid data were simulated."),
        enrich_test_gsea_note_genes_proteins = paste("Every feature is ranked by its correlation with the factor, and a",
            "term is called when its features sit nearer the top or the bottom of the list than chance (GSEA). The test",
            "shuffles features as if each were independent. Features that vary together, such as the members of one",
            "pathway or one complex, make its p-values too small: in a check on simulated data with 12 biological samples,",
            "terms with no link to the factor whose features were correlated at 0.3 were called at the 5 % level in 40.2",
            "to 43.0 % of cases (terms of 10 to 100 genes or proteins), against 5 % expected. The factor was learned from",
            "these same data. Use it to explore; it also uses features with missing values, whose correlations, computed",
            "on fewer samples, are noisier and sit more often near the ends of the list. In the same check, terms of 10",
            "genes or proteins in which half the features truly followed the main factor were found in 63 to 72 % of",
            "cases, a figure that comes with the excess on correlated terms above (in result tables of about 120 to 150",
            "terms, about 1 in 6 truly linked).",
            "It treats the biological samples as independent, which samples taken in pairs, repeatedly from one subject,",
            "or in batches are not."),
        enrich_test_gsea_note_metabolites = paste("Every feature is ranked by its correlation with the factor, and a term",
            "is called when its features sit nearer the top or the bottom of the list than chance (GSEA). The test shuffles",
            "features as if each were independent. Features that vary together, such as the members of one pathway, make",
            "its p-values too small: in a check on simulated metabolite data with 12 biological samples, terms with no",
            "link to the factor whose features were correlated at 0.3 were called at the 5 % level in 13.5 % of cases",
            "(terms of 5 to 10 metabolites), against 5 % expected. The factor was learned from these same data. Use it to",
            "explore; it also uses features with missing values, whose correlations, computed on fewer samples, are",
            "noisier and sit more often near the ends of the list. In the same check, terms of 5 metabolites, 2 of them",
            "truly following the main factor, were found in 1.0 to 1.3 % of cases (in result tables of about 90 terms,",
            "about 1 in 9 truly linked).",
            "It treats the biological samples as independent, which samples taken in pairs, repeatedly from one subject,",
            "or in batches are not."),
        enrich_test_camera_pr_note = paste("Like the correlation-adjusted test, but with the correlation between a",
            "term's features assumed to be 0.01 instead of measured (not offered in this release: the correlation-adjusted",
            "test with the correlation measured held the error rate in the check on simulated data). Terms whose features vary together more strongly get p-values that are too",
            "small: 29.9 % of terms with no link to the factor, correlated at 0.3, were called at the 5 % level.",
            "The factor was learned from these same data. Read the result as a list of terms worth a closer look, not as",
            "proof."),
        # KIT-184 (VAL-58 (1)): for the tests that have no design (the pre-ranked and permutation tests) on a model trained
        # on RUV-corrected values.
        enrich_values_ruv_note = paste("This model was trained on RUV-corrected values, but this result uses the",
            "normalised values before RUV correction without the unwanted factors the correction removed; variation the",
            "correction removed stays in them, which lowers power and, if it follows the factor, can bias the result."),
        # KIT-184 (VAL-53 (1), VAL-58, VAL-62, VAL-160; ticket 5.5). DRAFTS for the owner's wording pass (VAL-41), listed
        # verbatim in docs/qa/workflow-execution/kit/kit-184/wording.md. The setting is off unless chosen (VAL-160). The
        # numbers are those of the counted simulation (run 45353248796; kit-184/results.md): what the adjustment kept,
        # and its limit, in the words the results audit proposed (K184-R12). The setting, the label of the form with
        # every result of a test that has a design, and the sentences of each form.
        enrich_adjust_ruv = "Adjust for removed RUV factors",
        enrich_adjust_ruv_help = paste("Off unless you switch it on. On: each feature is tested for what it shares with",
            "the factor beyond the unwanted factors that the RUV correction removed from this data set. In a check on",
            "simulated data this kept false calls within the checked limit (at most 5.4 %) for the kinds of sets with no",
            "link to what was tested that the check decided on; sets made only of the most variable features were called",
            "more often (up to 9.4 %, 6.8 % without it), and so were sets in which half the features truly followed the",
            "main factor, tested against an unrelated covariate (up to 8.6 %, up to 7.0 % without it); it found more true",
            "terms than off (with as many factors removed as were planted; with four removed where two were planted it",
            "found fewer). Because the unwanted factors are estimated, sets of features that share them can still be",
            "called when that variation follows the factor: in a check on simulated data where it followed the groups,",
            "such sets were called on the model's main factor in 8.1 to 29.6 % of cases with the adjustment and in 2.0 to",
            "11.2 % without it. Off: the values before RUV are tested without those factors. On needs a data set to which",
            "the omic applied RUV. The setting does not apply to the pre-ranked test, the permutation test or STRING."),
        # %d: the removed factors in the design.
        enrich_adjust_label_yes = "Adjusted for removed RUV factors: yes (%d factor(s)).",
        enrich_adjust_label_no = "Adjusted for removed RUV factors: no.",
        # %d: the removed factors in the design, the same, the residual degrees of freedom. These three notes are filled
        # with sprintf(), so a per cent sign in them is written doubled.
        enrich_values_ruv_covariates_note = paste("This model was trained on values corrected for unwanted variation",
            "(RUV). The test reads the values before that correction and puts the %d unwanted factors the correction",
            "removed from this data set into its model, so a feature counts only for what it shares with the factor",
            "beyond those unwanted factors. For features measured in every run the regression is exactly that of the same",
            "model (the factor and these %d unwanted factors) fitted to the corrected values, with %d residual degrees of",
            "freedom; the test's variance trend reads each feature's average before the correction, so its p-values can",
            "differ where the correction moved a feature's level. It is not the test of the corrected values on the factor",
            "alone. The unwanted factors were estimated from these same data and are treated as known. Because the",
            "unwanted factors are estimated, sets of features that share them can still be called when that variation",
            "follows the factor: in a check on simulated data where it followed the groups, such sets were called on the",
            "model's main factor in 8.1 to 29.6 %% of cases with the adjustment and in 2.0 to 11.2 %% without it."),
        # An approximate set of factors (VAL-58 (6)): no sentence claims an exact match. %d: the removed factors in the
        # design, the residual degrees of freedom.
        enrich_values_ruv_covariates_note_approximate = paste("This model was trained on values corrected for unwanted",
            "variation (RUV). The test reads the values before that correction and puts the %d unwanted factors the",
            "correction removed from this data set into its model, so a feature counts only for what it shares with the",
            "factor beyond those unwanted factors. The model has %d residual degrees of freedom. It is not the test of the",
            "corrected values on the factor alone. The unwanted factors were estimated from these same data and are",
            "treated as known. Because the unwanted factors are estimated, sets of features that share them can still be",
            "called when that variation follows the factor: in a check on simulated data where it followed the groups,",
            "such sets were called on the model's main factor in 8.1 to 29.6 %% of cases with the adjustment and in 2.0 to",
            "11.2 %% without it."),
        # %d: the factors the correction declared; %s: the share of what it removed that lies within them.
        enrich_removed_approximate = paste("Runs were removed after the correction, or it removed variation outside these",
            "%d factors: %s of what it removed lies within them, so the match with the corrected values is approximate for",
            "every feature."),
        # %d: the removed factors in the design.
        enrich_values_adjusted_default_note = paste("This model was trained on the values before RUV correction. For this",
            "result the adjustment for removed RUV factors was switched on: the test puts the %d unwanted factors that the",
            "omic's RUV correction removed from this data set into its model, so a feature counts only for what it shares",
            "with the factor beyond those unwanted factors. If the factor itself follows them, little is left to test (see",
            "the share below). The unwanted factors were estimated from these same data and are treated as known. Because",
            "the unwanted factors are estimated, sets of features that share them can still be called when that variation",
            "follows the factor: in a check on simulated data where it followed the groups, such sets were called on the",
            "model's main factor in 8.1 to 29.6 %% of cases with the adjustment and in 2.0 to 11.2 %% without it."),
        # The form without the adjustment on a model trained on RUV-corrected values: the default (VAL-160). Its numbers
        # are the ranges over the three simulated data sets where the unwanted variation followed the groups.
        enrich_values_ruv_unadjusted_note = paste("This model was trained on RUV-corrected values. This result is the form",
            "without the adjustment for removed RUV factors, which is the default (the adjustment is optional): it tests",
            "the normalised values before RUV correction without the unwanted factors the correction removed. Its degrees",
            "of freedom are counted correctly, but variation the correction removed stays in the values, which lowers",
            "power. In a check on simulated data where the unwanted variation followed the groups, sets of features",
            "sharing it were called on the model's main factor in 2.0 to 11.2 % of cases in this form and in 8.1 to 29.6 %",
            "with the adjustment."),
        # The counts of a result with the removed factors, in three sentences (the second is left out for an approximate
        # set of factors). %d: the removed factors in the design, biological samples, residual degrees of freedom.
        enrich_test_covariates_design = paste("Design: intercept, factor and %d removed unwanted factor(s); %d biological",
            "samples, %d residual degrees of freedom."),
        # %d: the features used that were measured in every run, the features used.
        enrich_test_covariates_runs = paste("%d of the %d features used were measured in every run; for the others the",
            "correction worked feature by feature and the match with the corrected values is approximate."),
        # %s: the share of the factor's variance the removed factors explain, the share expected by chance.
        enrich_test_covariates_share = paste("The removed unwanted factors explain %s of this factor's variation across",
            "these samples (about %s would be expected for an unrelated factor)."),
        # The three composed (%d, %d, %d; %d, %d; %s, %s), as shown for an exact set of factors.
        enrich_test_covariates_counts = paste("Design: intercept, factor and %d removed unwanted factor(s); %d biological",
            "samples, %d residual degrees of freedom. %d of the %d features used were measured in every run; for the others",
            "the correction worked feature by feature and the match with the corrected values is approximate. The removed",
            "unwanted factors explain %s of this factor's variation across these samples (about %s would be expected for an",
            "unrelated factor)."),
        enrich_factor_removed_overlap = paste("More than half of this factor's variation across the samples of this data",
            "set follows the unwanted factors removed from it. Little is left to test here, and in this data set the factor",
            "may reflect variation the correction treated as unwanted."),
        # %d: the removed factors, the residual degrees of freedom, the smallest number at which the rate was checked.
        enrich_test_few_df = paste("With %d removed factor(s) in the model this test has %d residual degrees of freedom,",
            "fewer than the smallest number (%d) at which its error rate was checked."),
        # camera_floor with the removed factors (KIT-191 8.2 item 4; D3 (a): its measured power, arm 3 of the counted
        # simulation, sets of 10 on the main factor, the gene and protein data sets of both deciding cells).
        enrich_floor_covariates_note = paste("False calls are at most those of the correlation-adjusted test with the same",
            "unwanted factors (counting negative correlations as zero never makes a p-value smaller). In a check on",
            "simulated data with 12 biological samples and 2 removed factors, terms of 10 features in which half the",
            "features truly followed the main factor were found in 13 to 32 % of cases with these factors in the model (the",
            "correlation-adjusted test: 30 to 48 %)."),
        enrich_caution = paste("Exploratory: the factor was learned from these same data, so these terms are worth a",
            "closer look, not proof. P-values are adjusted within this factor, data set and database only."),
        # %s: the factor's variance explained, the largest factor's (r2_percent of the data set's model features), the
        # kind of simulated data, the kind of term, the entry's test's figure (row E of the evidence table), the test.
        enrich_factor_strength = paste("In this data set, this factor explains %s of the variance of the features the",
            "model was trained on (the factor that explains most explains %s). In a check on simulated %s, on factors",
            "weaker than the main effect, %s were called in up to %s of cases by the %s, although they had no link to the",
            "simulated biology; no error rate is claimed for such factors."),
        # VAL-115 (4): the kind of simulated data, the lipid case in one bracket after it.
        enrich_strength_kind_genes_proteins = "gene and protein data",
        enrich_strength_kind_metabolites = "metabolite data",
        enrich_strength_kind_lipids = "metabolite data (no lipid data were simulated)",
        enrich_strength_terms_genes_proteins = paste("terms made of features the model used, or of features that vary",
            "together with them"),
        enrich_strength_terms_metabolites = "terms of metabolites the model used",
        # The tests as a sentence names them.
        enrich_test_phrase_camera = "correlation-adjusted test",
        enrich_test_phrase_camera_floor = "correlation-adjusted test with negative correlations counted as zero",
        enrich_test_phrase_permutation = "permutation test",
        enrich_test_phrase_gsea = "pre-ranked test",
        enrich_test_phrase_camera_pr = "correlation-adjusted test with an assumed correlation",
        # %d: the data set's biological samples (fewer than 11, the smallest number at which the error rates were decided).
        enrich_test_few_samples_camera_genes_proteins = paste("This data set has %d biological samples, fewer than the",
            "smallest number (11) at which the tests' error rates were decided; at 6 samples, terms with no link to the",
            "factor were called at the 5 %% level in up to 5.0 %% of cases (terms of the most variable features: up to",
            "9.9 %%), and terms of 30 features half following the factor were found in 23 to 44 %%."),
        enrich_test_few_samples_camera_floor_genes_proteins = paste("This data set has %d biological samples, fewer than",
            "the smallest number (11) at which the tests' error rates were decided; at 6 samples, terms with no link to the",
            "factor were called at the 5 %% level in up to 4.1 %% of cases (terms of the most variable features: up to",
            "2.7 %%), and terms of 30 features half following the factor were found in 0 to 12 %%."),
        enrich_test_few_samples_permutation_genes_proteins = paste("This data set has %d biological samples, fewer than",
            "the smallest number (11) at which the tests' error rates were decided; at 6 samples, terms with no link to the",
            "factor were called at the 5 %% level in up to 5.3 %% of cases (terms of the most variable features: up to",
            "5.3 %%), and terms of 30 features half following the factor were found in 0.7 to 1.9 %%."),
        enrich_test_few_samples_gsea_genes_proteins = paste("This data set has %d biological samples, fewer than the",
            "smallest number (11) at which the tests' error rates were decided; at 6 samples, terms with no link to the",
            "factor were called at the 5 %% level in up to 5.2 %% of cases (64.7 %% when the term's features vary",
            "together; terms of the most variable features: up to 6.2 %%), and terms of 30 features half following the",
            "factor were found in 95 to 99 %%."),
        enrich_test_few_samples_camera_metabolites = paste("This data set has %d biological samples, fewer than the",
            "smallest number (11) at which the tests' error rates were decided; at 6 samples, terms with no link to the",
            "factor were called at the 5 %% level in up to 3.7 %% of cases (terms of the most variable features: up to",
            "3.8 %%), and terms of 5 metabolites, 2 following the factor, were found in 0.1 to 0.2 %%."),
        enrich_test_few_samples_camera_floor_metabolites = paste("This data set has %d biological samples, fewer than",
            "the smallest number (11) at which the tests' error rates were decided; at 6 samples, terms with no link to the",
            "factor were called at the 5 %% level in up to 2.2 %% of cases (terms of the most variable features: up to",
            "0.3 %%), and terms of 5 metabolites, 2 following the factor, were found in 0 %%."),
        enrich_test_few_samples_permutation_metabolites = paste("This data set has %d biological samples, fewer than",
            "the smallest number (11) at which the tests' error rates were decided; at 6 samples, terms with no link to the",
            "factor were called at the 5 %% level in up to 5.0 %% of cases (terms of the most variable features: up to",
            "5.0 %%), and terms of 5 metabolites, 2 following the factor, were found in 0.1 to 0.3 %%."),
        enrich_test_few_samples_gsea_metabolites = paste("This data set has %d biological samples, fewer than the",
            "smallest number (11) at which the tests' error rates were decided; at 6 samples, terms with no link to the",
            "factor were called at the 5 %% level in up to 4.6 %% of cases (33.4 %% when the term's features vary",
            "together; terms of the most variable features: up to 4.2 %%), and terms of 5 metabolites, 2 following the",
            "factor, were found in 0.7 %%."),
        # %d: features used, biological samples, ranked features with a missing value, second features on one identifier.
        enrich_test_counts = paste("Used by the correlation-adjusted test: %d features measured in all %d biological",
            "samples (%d ranked features with an identifier had a missing value and were left out; %d set aside as a second feature on one",
            "identifier)."),
        enrich_test_counts_permutation = paste("Used by the permutation test: %d features measured in all %d biological",
            "samples (%d ranked features with an identifier had a missing value and were left out; %d set aside as a",
            "second feature on one identifier)."),
        # %d: reorderings, biological samples; %s: the smallest p-value; %d: terms tested, the smallest number of terms
        # called together (5.2 step 8, as corrected by AUDIT-KIT-193 A05: BH calls k* terms with p_(k*) <= 0.05 k* / m,
        # and p_(k*) >= p_min, so k* >= m p_min / 0.05).
        enrich_test_permutation_floor = paste("With %d reorderings of %d biological samples no p-value can be below %s;",
            "among the %d terms tested, none can be called at a 5 %% false discovery rate unless at least %d are called",
            "together."),
        # %d: reorderings, biological samples; %s: the smallest p-value, above 0.05 (AUDIT-KIT-191-impl I01).
        enrich_test_permutation_floor_none = paste("With %d reorderings of %d biological samples no p-value can be below",
            "%s, which is above 0.05, so no term can be called at a 5 %% false discovery rate."),
        # %s: seconds; %d: features, biological samples, terms, reorderings (the scripted run, D4 (c)).
        enrich_test_permutation_time = paste("The permutation test is predicted to take about %s seconds (%d features, %d",
            "biological samples, %d terms, %d reorderings); the prediction was measured on a test computer, and this",
            "computer may be slower or faster."),
        enrich_permutations = "Number of random reorderings (permutation test)",
        enrich_permutations_1000 = "1,000 (the number checked on simulated data)",
        enrich_permutations_10000 = "10,000 (smaller p-values possible; not checked separately)",
        # The evidence table beside the selector (D6 (a)); %s: the kind of simulated data.
        enrich_evidence_heading = paste("How each test behaved in a check on simulated %s data with 12 biological",
            "samples (5 %% level)"),
        enrich_evidence_lipids = "No lipid data were simulated; these are the figures for metabolites.",
        enrich_evidence_measure = "Kind of term",
        enrich_evidence_A = "Terms with no link, factor replaced by values unrelated to the data",
        enrich_evidence_B = paste("Terms made only of the most variable features (those the model selects), values",
            "unrelated to the data"),
        enrich_evidence_C = "Terms with no link, main factor",
        enrich_evidence_D = "Terms with no link whose features vary together (correlation 0.3), main factor",
        enrich_evidence_E_genes_proteins = paste("Weaker factors: terms of features the model used, or of features that",
            "vary together with them, largest"),
        enrich_evidence_E_metabolites = "Weaker factors: terms of metabolites the model used, largest",
        enrich_evidence_F = "At least one false call among a result's unrelated terms",
        enrich_evidence_C3 = "Terms of compounds that carry 2 or 3 measured features, no link",
        enrich_evidence_G = "Terms of 10 genes or proteins, half truly following the main factor: found",
        enrich_evidence_H = "The same, 30 features: found",
        enrich_evidence_I = "Terms of 5 metabolites, 2 of them following the main factor: found",
        enrich_evidence_kind_genes_proteins = "gene and protein",
        enrich_evidence_kind_metabolites = "metabolite",
        # The first sentence of `enrich_direction` (%s: the feature, its data set; %d: +1 or -1), for camera results.
        enrich_direction_oriented = paste("Direction: the factor is oriented so that its largest weight (%s in %s) is",
            "positive (multiplied by %d)."),
        enrich_direction_camera = paste("Up: the term's features rise with the factor more than the other features do",
            "on average; Down: they rise less, or fall more, than the other features do on average.")
    )
}

#' The integration GUI's configuration (KIT-138): the shell group only, in the single-omic shell's form, validated by
#' `validateOmicGuiShell`. KIT-140 adds its tab to `tabs`, KIT-063 the Results tab, KIT-142 the Record tab.
#' @return `list(shell, texts)`.
#' @noRd
integrationGuiConfig <- function() {
    texts <- integrationGuiTexts()
    tab <- function(status_key, value, title, icon, step) {
        omicGuiTab(status_key, value, value, paste0("integration_", value), texts[[title]], icon,
            paste0("integration-tab-", value), texts[[step]])
    }
    tabs <- list(
        tab("integration_import", "import", "tab_import", "file-import", "step_import"),
        tab("integration_alignment", "alignment", "tab_alignment", "th", "step_alignment"),
        tab("integration_model", "model", "tab_model", "cubes", "step_model"),
        # KIT-063: the Results tab.
        tab("integration_results", "results", "tab_results", "chart-bar", "step_results"),
        # KIT-141: the Enrichment tab (optional; it opens the Record tab as the Results tab completes).
        tab("integration_enrichment", "enrichment", "tab_enrichment", "project-diagram", "step_enrichment"),
        # KIT-142: the Record tab (the analysis record, the study parameters and the report).
        tab("integration_record", "record", "tab_record", "file-alt", "step_record"))
    shell <- list(omic_type = "integration", prefix = "integration", testid = "integration",
        label = texts[["workflow_label"]], label_lower = "integration", header_icon = NULL,
        tabset_id = "integration_tabs", tab_style = "icon_in_title", routes = "integration", tabs = tabs,
        bus = list(state_manager = NULL, tab_status = NULL, integration_scan = NULL, integration_import = NULL,
            integration_alignment = NULL, integration_model = NULL, integration_results = NULL, integration_enrichment = NULL,
            integration_record = NULL),
        bus_registration = "userData")
    validateOmicGuiShell(shell, omicGuiConfigSchema())
    list(shell = shell, texts = texts)
}

#' The home page's Integration tile with the gate on: selectable like an omic (the omics it combines were processed
#' earlier, so no other omic needs to be selected in this session).
#' @noRd
integrationHomeTile <- function(texts = integrationGuiTexts()) {
    testid(
        shiny::div(
            id = "integration_box",
            class = "omic-selection-box",
            onclick = "toggleOmicSelection('integration')",
            shiny::h3(shiny::icon("project-diagram"), texts[["home_tile_title"]]),
            shiny::p(texts[["home_tile_text"]])
        ),
        "tile-integration"
    )
}

#' The home page's selection script: as it is with the gate off; with the gate on, Integration is selected like an
#' omic and its tile is never disabled. Refused when the script no longer holds the lines it changes.
#' @param script The home page's script.
#' @return The script.
#' @noRd
integrationHomeScript <- function(script) {
    if (!integrationGuiEnabled()) return(script)
    edits <- c("if (omic === 'integration') return; // Integration is auto-enabled" =
            "// Integration is selected like an omic (KIT-138)",
        "if (selectedOmics.length >= 2) {" = "if (true) { // KIT-138: the Integration tile stays enabled")
    for (from in names(edits)) {
        if (!grepl(from, script, fixed = TRUE)) {
            cli::cli_abort("the home page script no longer holds {.val {from}}", class = "integration_gui_invalid")
        }
        script <- sub(from, edits[[from]], script, fixed = TRUE)
    }
    script
}

#' The "module not available" panel app_server shows for an omic without a module (the integration module with the
#' gate off shows the same panel, so the app is as before).
#' @noRd
integrationGuiUnavailablePanel <- function() {
    shiny::div(
        class = "alert alert-warning",
        shiny::h4("Module not yet implemented"),
        shiny::p(paste("The", "integration", "workflow module is under development."))
    )
}

# ------------------------------------------------------------------------------------------------
# The Import tab: the project scan (KIT-138).
# ------------------------------------------------------------------------------------------------

INTEGRATION_GUI_OMICS <- c("lipidomics", "metabolomics", "proteomics", "transcriptomics")

#' One value of a manifest field read by `jsonlite::read_json(simplifyVector = FALSE)`: NA when absent, null or empty.
#' @noRd
integrationGuiScalar <- function(value) {
    if (is.null(value) || !length(value)) return(NA)
    value <- unlist(value)
    if (!length(value)) NA else value[[1L]]
}

#' The correction for unwanted variation an input names, in words (KIT-138): not applied for a normalised matrix; for
#' a corrected one the declared RUV approach (KIT-129's `ruv_design`), else the approach of earlier versions (its method
#' and replicate sets), then the number of factors, overall or per data set (lipidomics, metabolomics).
#' @param manifest The input's manifest.
#' @param texts `integrationGuiTexts()`.
#' @return One string.
#' @noRd
integrationCorrectionText <- function(manifest, texts = integrationGuiTexts()) {
    correction <- manifest$correction
    # KIT-153: the correction an input names is the omic's (its corrected matrices), whichever matrix a model uses.
    if (!is.list(correction)) return(texts[["correction_none"]])
    design <- integrationGuiScalar(correction$ruv_design)
    method <- integrationGuiScalar(correction$method)
    sets <- integrationGuiScalar(correction$replicate_structure)
    key <- paste0("design_", design)
    text <- if (!is.na(design) && key %in% names(texts)) {
        sprintf(texts[["correction_design"]], texts[[key]])
    } else if (!is.na(method)) {
        method_text <- switch(as.character(method), RUVIII_C_Varying = texts[["method_ruviii_c"]],
            RUVIII = texts[["method_ruviii"]], as.character(method))
        sets_text <- if (identical(sets, "technical_replicates")) texts[["sets_technical_replicates"]] else texts[["sets_group"]]
        sprintf(texts[["correction_legacy"]], method_text, sets_text)
    } else {
        texts[["correction_unknown"]]
    }
    if (is.list(correction$views) && length(correction$views)) {
        per_view <- vapply(correction$views, function(view) {
            assay <- as.character(integrationGuiScalar(view$assay))
            if (isTRUE(as.logical(integrationGuiScalar(view$applied)))) {
                sprintf(texts[["correction_view_applied"]], assay, as.character(integrationGuiScalar(view$k)))
            } else {
                sprintf(texts[["correction_view_not_applied"]], assay)
            }
        }, character(1))
        return(paste0(text, "; ", paste(per_view, collapse = "; ")))
    }
    k <- integrationGuiScalar(correction$ruv_k)
    if (is.na(k)) k <- integrationGuiScalar(correction$k)
    if (!is.na(k)) text <- paste0(text, "; ", sprintf(texts[["correction_k"]], as.character(k)))
    text
}

#' Scan a project's integration folder for the omics' integration inputs (KIT-138). Each omic's folder
#' (`<project>/integration/<omic>/`) is read by the contract's reader (`readIntegrationInput`): a folder without the
#' input is not saved; one the reader refuses is listed with the reader's message; one whose manifest has changed
#' since the import is marked changed. Nothing is written.
#' @param integration_dir The project's integration folder.
#' @param imported The imported manifests' sha256, named by omic (NULL before an import).
#' @param texts `integrationGuiTexts()`.
#' @return `list(dir, omics, views)`: one row per omic folder found (`omic`, `dir`, `status`: `ready`, `changed`,
#'   `not_saved` or `refused`; `reason`, `data_type`, `values`, `correction`, `saved`, `manifest_sha256`) and one row
#'   per data set of a usable omic (`view`, `omic`, `features`, `runs`, `biological_samples`, `groups`).
#' @noRd
integrationScanProject <- function(integration_dir, imported = NULL, texts = integrationGuiTexts()) {
    folders <- file.path(integration_dir, INTEGRATION_GUI_OMICS)
    folders <- folders[dir.exists(folders)]
    omics <- list()
    views <- list()
    for (dir in folders) {
        omic <- basename(dir)
        row <- list(omic = omic, dir = normalizePath(dir, winslash = "/", mustWork = TRUE), status = "ready",
            reason = NA_character_, data_type = NA_character_, values = NA_character_, correction = NA_character_,
            saved = NA_character_, manifest_sha256 = NA_character_)
        manifest_path <- file.path(dir, INTEGRATION_INPUT_MANIFEST)
        if (!file.exists(manifest_path)) {
            row$status <- "not_saved"
            row$reason <- texts[["status_not_saved"]]
            omics[[omic]] <- row
            next
        }
        input <- tryCatch(readIntegrationInput(dir), integration_input_invalid = function(error) error)
        if (inherits(input, "error")) {
            row$status <- "refused"
            row$reason <- conditionMessage(input)
            omics[[omic]] <- row
            next
        }
        manifest <- input$manifest
        row$manifest_sha256 <- digest::digest(file = manifest_path, algo = "sha256")
        if (!is.null(imported) && omic %in% names(imported) && !identical(unname(imported[[omic]]), row$manifest_sha256)) {
            row$status <- "changed"
            row$reason <- texts[["status_changed"]]
        }
        row$data_type <- as.character(integrationGuiScalar(manifest$workflow))
        # KIT-153: every input holds the normalised values; one whose omic applied RUV holds the corrected ones too.
        row$values <- texts[[if (any(vapply(input$views, function(view) isTRUE(view$corrected), logical(1)))) "values_both" else
            "values_normalised"]]
        row$correction <- integrationCorrectionText(manifest, texts)
        written <- as.character(integrationGuiScalar(manifest$written_at))
        row$saved <- if (is.na(written)) NA_character_ else
            sprintf(texts[["saved_at"]], sub("^(\\d{4}-\\d{2}-\\d{2})T(\\d{2}:\\d{2}).*$", "\\1 \\2", written))
        omics[[omic]] <- row
        samples <- input$samples
        biological <- manifest$biological_sample_column %||% integrationInputBiologicalSampleColumn(samples,
            integrationGuiScalar(manifest$technical_replicate_column), manifest$sample_id_column)
        biological_count <- if (is.character(biological) && biological %in% names(samples))
            length(unique(samples[[biological]])) else NA_integer_
        group_count <- length(unique(samples[[manifest$group_column]]))
        for (view in names(input$views)) {
            # KIT-147: the reader gives each data set's size; no matrix is read here.
            views[[view]] <- data.frame(view = view, omic = omic, features = input$views[[view]]$rows,
                runs = input$views[[view]]$columns, biological_samples = biological_count, groups = group_count,
                stringsAsFactors = FALSE)
        }
    }
    omic_table <- if (length(omics)) do.call(rbind, lapply(omics, as.data.frame, stringsAsFactors = FALSE)) else
        data.frame(omic = character(), dir = character(), status = character(), reason = character(),
            data_type = character(), values = character(), correction = character(), saved = character(),
            manifest_sha256 = character(), stringsAsFactors = FALSE)
    rownames(omic_table) <- NULL
    view_table <- if (length(views)) do.call(rbind, views) else data.frame(view = character(), omic = character(),
        features = integer(), runs = integer(), biological_samples = integer(), groups = integer(),
        stringsAsFactors = FALSE)
    rownames(view_table) <- NULL
    list(dir = integration_dir, omics = omic_table, views = view_table)
}

#' The scan's omics as the Import tab shows them: one row per omic folder, the headers in words.
#' @noRd
integrationScanDisplay <- function(scan, texts = integrationGuiTexts()) {
    omics <- scan$omics
    status <- vapply(seq_len(nrow(omics)), function(i) switch(omics$status[[i]], ready = texts[["status_ready"]],
        changed = texts[["status_changed"]], not_saved = texts[["status_not_saved"]], texts[["status_refused"]]),
        character(1))
    table <- data.frame(tools::toTitleCase(omics$omic), omics$data_type, omics$values, omics$correction, omics$saved,
        status, stringsAsFactors = FALSE)
    names(table) <- unname(texts[c("col_omic", "col_data_type", "col_values", "col_correction", "col_saved", "col_status")])
    table
}

#' The data sets of the omics an import can take: those of the omics marked ready, or changed since the import (a new
#' import is the remedy).
#' @noRd
integrationScanReadyViews <- function(scan) {
    ready <- scan$omics$omic[scan$omics$status %in% c("ready", "changed")]
    scan$views[scan$views$omic %in% ready, , drop = FALSE]
}

#' The import a selection asks for (KIT-138): the folders of the omics with a selected data set and the selected data
#' sets in the scan's order. Refused in words (`integration_refused`): a data set of an omic not ready, fewer than two
#' omics (D7), fewer than two data sets.
#' @param scan `integrationScanProject()`'s value.
#' @param selected The selected data sets.
#' @return `list(dirs, omics, views)`.
#' @noRd
integrationImportRequest <- function(scan, selected, texts = integrationGuiTexts()) {
    selected <- unique(as.character(selected))
    ready <- integrationScanReadyViews(scan)
    refuse <- function(message, class) integrationRefuse(message, c(class, "integration_gui_refusal"))
    if (!all(selected %in% ready$view)) refuse(texts[["import_need_ready"]], "integration_input_not_ready")
    views <- ready$view[ready$view %in% selected]
    omics <- unique(ready$omic[ready$view %in% views])
    if (length(omics) < 2L) refuse(sprintf(texts[["import_need_two_omics"]], length(omics)), "integration_too_few_omics")
    if (length(views) < 2L) refuse(sprintf(texts[["import_need_two_views"]], length(views)), "integration_too_few_views")
    list(dirs = scan$omics$dir[match(omics, scan$omics$omic)], omics = omics, views = views)
}

#' Commit an import on the integration store (KIT-138): `importIntegrationInputs` on the folders, then the summary the
#' tabs read (the state's reference, the omics, the selected data sets in the state's order, the manifests' sha256 and
#' the counts).
#' @param manager The integration store (`newIntegrationStateOwner`).
#' @param request `integrationImportRequest()`'s value.
#' @return The import's summary.
#' @noRd
integrationCommitImport <- function(manager, request) {
    reference <- importIntegrationInputs(manager, request$dirs)
    # KIT-149 (audit ENG-05): the names, the design and the manifests only; no data set is read.
    object <- integrationImportedObject(manager, fields = c("design_matrix", "inputs"))
    views <- names(object@views)[names(object@views) %in% request$views]
    design <- object@design_matrix
    view_omics <- unlist(unname(lapply(object@inputs, function(entry) stats::setNames(rep(entry$omic, length(entry$views)),
        entry$views))))
    used <- design$omic %in% unique(view_omics[views])
    list(reference = reference, omics = names(object@inputs), views = views, manifests = integrationStateManifests(object),
        runs = sum(used), biological_samples = length(unique(design$biological_sample[used])), committed_at = Sys.time())
}

#' The imported inputs state, every field hydrated, or (KIT-149, audit ENG-05) only the named fields: the data sets
#' left unread are zero-row prototypes that keep their names and columns.
#' @param manager The integration store.
#' @param fields TRUE (every field), or the fields to read (`readWorkflowExecutionComponents`).
#' @noRd
integrationImportedObject <- function(manager, fields = TRUE) {
    withWorkflowStateScope(manager, function() {
        source <- getWorkflowStateReference(manager, INTEGRATION_INPUTS_STATE)
        readWorkflowExecutionComponents(manager, source, fields = fields)
    })
}

# The fields the Sample Matching tab's preview reads (KIT-149, audit ENG-05): every field but the data sets, whose
# values the preview reads in row pages (`integrationAlignmentPreview(manager = )`).
INTEGRATION_PREVIEW_FIELDS <- c("features", "design_matrix", "inputs")

# ------------------------------------------------------------------------------------------------
# The Sample Matching tab (KIT-139).
# ------------------------------------------------------------------------------------------------

#' The sample-matching preview of the imported data sets (KIT-139): the model layer's alignment
#' (`integrationAlignViews`) under the parameters as resolved for the selected data sets (the provisional defaults of
#' `integrationMofaParameterSpec`, never set here), with its facts counted for the tab. The few-samples warning (D7)
#' is kept as a note, not raised. Refused (`integration_refused`): an input changed since the import, fewer than two
#' data sets, no shared biological sample, a sample in two groups.
#' @param object The imported inputs state (`integrationImportedObject`).
#' @param views The selected data sets.
#' @param match NULL, or the user's sample match table as uploaded (KIT-063 item A): checked and applied as the model
#'   applies it; a table that gives no run another biological sample changes nothing.
#' @param manager KIT-149 (audit ENG-05): the integration store, whose data sets are then read in row pages
#'   (`integrationStateViewPages`; `object` needs only `INTEGRATION_PREVIEW_FIELDS`); NULL reads them from `object`,
#'   hydrated. Either way the alignment is the model's without its values (`integrationAlignViewsPartitioned(values =
#'   FALSE)`): the counts and the sample table of `integrationAlignViews`, `aligned$views` NULL.
#' @return `list(parameters, aligned, few_samples, keys, views, replicates, all_samples, match)`; `match` is the
#'   match's check (`integrationResolveSampleMatch`) or NULL.
#' @noRd
integrationAlignmentPreview <- function(object, views, match = NULL, manager = NULL) {
    integrationAssertInputsCurrent(object)
    parameters <- integrationResolveMofaParameters(object, if (is.null(match)) list(views = views) else
        list(views = views, sample_match_table = match))
    checked <- if (!is.null(match)) integrationResolveSampleMatch(object, match)
    few <- FALSE
    pages <- if (is.null(manager)) integrationObjectViewPages(object) else
        integrationStateViewPages(manager, withWorkflowStateScope(manager, function() {
            getWorkflowStateReference(manager, INTEGRATION_INPUTS_STATE)
        }))
    aligned <- withCallingHandlers(integrationAlignViewsPartitioned(object, parameters, pages, values = FALSE),
        integration_few_samples = function(warning) {
        few <<- TRUE
        invokeRestart("muffleWarning")
    })
    design <- integrationApplySampleMatch(object@design_matrix, parameters$sample_match_table)
    view_omic <- unlist(unname(lapply(object@inputs, function(entry) stats::setNames(rep(entry$omic, length(entry$views)),
        entry$views))))
    all_samples <- character()
    replicates <- do.call(rbind, lapply(parameters$views, function(view) {
        runs <- names(object@views[[view]])[-1L]
        rows <- design[design$omic == view_omic[[view]] & design$run %in% runs, , drop = FALSE]
        all_samples <<- union(all_samples, rows$biological_sample)
        per_sample <- table(rows$biological_sample)
        standards <- sum(integrationInternalStandards(object@features[[view]]))
        data.frame(view = view, omic = view_omic[[view]], runs = length(runs), biological_samples = length(per_sample),
            repeated = sum(per_sample > 1L), standards = standards, stringsAsFactors = FALSE)
    }))
    keys <- data.frame(omic = names(object@inputs), column = vapply(object@inputs, function(entry)
        as.character(entry$biological_sample_column), character(1)), stringsAsFactors = FALSE, row.names = NULL)
    out <- list(parameters = parameters, aligned = aligned, few_samples = few, keys = keys, views = parameters$views,
        replicates = replicates, all_samples = length(all_samples))
    # KIT-260 (VAL-241): the marked rows named per data set (OD-4 (a)), and those the stored pattern matched through the
    # alternatives the corrected MS-DIAL default removed (OD-3 (a)).
    out$standard_names <- integrationInternalStandardNames(object, parameters$views)
    out$standard_legacy <- integrationInternalStandardLegacy(object, parameters$views)
    if (!is.null(checked)) out$match <- checked
    out
}

#' Read a sample match table the user uploaded (KIT-063 item A): tab-separated, as the tab's download writes it, every
#' column read as text. Refused in words (`integration_sample_match_invalid`) when the file cannot be read as a table.
#' @param path The uploaded file.
#' @return A data frame.
#' @noRd
integrationReadSampleMatch <- function(path) {
    tryCatch(utils::read.delim(path, colClasses = "character", check.names = FALSE, na.strings = character(),
        stringsAsFactors = FALSE, quote = "", comment.char = ""),
        error = function(error) integrationRefuse(sprintf("the sample match file cannot be read as a tab-separated table: %s",
            conditionMessage(error)), "integration_sample_match_invalid"))
}

#' What an uploaded sample match does, in words (KIT-063 item A): runs given another biological sample, and the rows
#' that name nothing the import has (ignored), with their counts.
#' @param report `integrationResolveSampleMatch()$report`.
#' @return A character vector.
#' @noRd
integrationSampleMatchFacts <- function(report, texts = integrationGuiTexts()) {
    facts <- if (report$runs_changed > 0L) sprintf(texts[["match_changed"]], report$runs_changed) else texts[["match_unchanged"]]
    ignored <- c(match_unknown_omics = report$unknown_omics, match_unknown_runs = report$unknown_runs,
        match_unknown_names = report$unknown_names, match_rows_empty = report$rows_empty)
    for (key in names(ignored)) if (ignored[[key]] > 0L) facts <- c(facts, sprintf(texts[[key]], ignored[[key]]))
    unname(facts)
}

#' Save the confirmed sample match beside the integration design (KIT-063 item A): the runs it gives another
#' biological sample, tab-separated.
#' @param match The canonical match.
#' @param dir The folder.
#' @return The file's path.
#' @noRd
integrationSaveSampleMatch <- function(match, dir) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    path <- file.path(dir, "integration_sample_match.tsv")
    vroom::vroom_write(as.data.frame(match, stringsAsFactors = FALSE), path, delim = "\t", na = "NA", progress = FALSE)
    path
}

#' The facts the Sample Matching tab states, in words, each from the preview's counts and the parameter's value (so a
#' changed default changes the sentence, not a hard-coded rule).
#' @param preview `integrationAlignmentPreview()`'s value.
#' @return A character vector.
#' @noRd
integrationAlignmentFacts <- function(preview, texts = integrationGuiTexts()) {
    parameters <- preview$parameters
    summary <- preview$aligned$view_summary
    text <- function(prefix, value) texts[[paste0(prefix, value)]]
    # Counted over the union of the data sets' biological samples, whichever rule keeps them.
    total <- preview$all_samples
    missing <- total - preview$aligned$shared_samples
    facts <- c(sprintf(texts[["fact_samples"]], total, preview$aligned$shared_samples),
        sprintf(text("fact_missing_", parameters$missing_samples), missing))
    for (i in seq_len(nrow(preview$replicates))) {
        row <- preview$replicates[i, ]
        facts <- c(facts, if (row$repeated == 0L) sprintf(texts[["fact_replicates_single"]], row$view, row$runs) else
            sprintf(text("fact_replicates_", parameters$technical_replicates), row$view, row$runs,
                row$biological_samples, row$repeated))
        if (row$standards > 0L) {
            facts <- c(facts, sprintf(text("fact_standards_", parameters$internal_standards), row$view, row$standards))
            # KIT-260 (OD-4 (a)): which rows, "display name (identifier)".
            named <- preview$standard_names[[row$view]]
            if (length(named)) facts <- c(facts, sprintf(omicGuiItsdTexts()[["integration_names"]], row$view, paste(named, collapse = "; ")))
        }
        sparse <- summary$sparse_dropped[summary$view == row$view]
        if (length(sparse) && sparse > 0L) facts <- c(facts, sprintf(texts[["fact_sparse"]], row$view, sparse))
    }
    # KIT-260 (OD-3 (a)): the note, only where marked rows were matched through the removed alternatives.
    note <- omicsItsdLegacyNote(preview$standard_legacy, unit = "data set")
    if (!is.null(note)) facts <- c(facts, note)
    unname(facts)
}

#' The Sample Matching tab's table of data sets: runs, biological samples, internal standards removed, features left
#' out and used.
#' @noRd
integrationAlignmentViewTable <- function(preview, texts = integrationGuiTexts()) {
    summary <- preview$aligned$view_summary
    replicates <- preview$replicates[match(summary$view, preview$replicates$view), , drop = FALSE]
    table <- data.frame(summary$view, replicates$runs, replicates$biological_samples, summary$internal_standards_dropped,
        summary$sparse_dropped, summary$features_model, stringsAsFactors = FALSE)
    names(table) <- unname(texts[c("col_view", "col_runs", "col_biological_samples", "col_standards_removed",
        "col_sparse", "col_features_used")])
    table
}

#' The Sample Matching tab's sample table: one row per biological sample, its group, the number of data sets it is
#' measured in, and per data set the runs it comes from or "missing" (D4).
#' @noRd
integrationAlignmentSampleTable <- function(preview, texts = integrationGuiTexts()) {
    samples <- preview$aligned$samples
    table <- data.frame(samples$biological_sample, samples$group, samples$views_present, stringsAsFactors = FALSE)
    names(table) <- unname(texts[c("col_sample", "col_group", "col_views_present")])
    for (view in preview$views) {
        runs <- samples[[paste0(view, "_runs")]]
        table[[view]] <- ifelse(nzchar(runs), gsub(";", ", ", runs, fixed = TRUE), texts[["cell_missing"]])
    }
    table
}

#' A refusal of the sample matching in words: the tab's own sentence for the refusals it names (KIT-139: no shared
#' sample, fewer than two data sets, a changed input), the plain sentence of the refusal's class for the others
#' (KIT-140, `integrationRefusalText`). The check's own message is the refusal's detail, shown apart.
#' @param condition The `integration_refused` condition.
#' @param views The selected data sets.
#' @noRd
integrationAlignmentRefusal <- function(condition, views, texts = integrationGuiTexts()) {
    reason <- if (inherits(condition, "integration_no_shared_samples")) {
        sprintf(texts[["refuse_no_shared"]], paste(views, collapse = ", "))
    } else if (inherits(condition, "integration_too_few_views")) {
        sprintf(texts[["refuse_too_few_views"]], length(views))
    } else if (inherits(condition, "integration_input_stale")) {
        texts[["refuse_changed"]]
    } else {
        integrationRefusalText(condition, texts)$sentence
    }
    sprintf(texts[["align_refused"]], reason)
}

# ------------------------------------------------------------------------------------------------
# Refusals in plain words (KIT-140): the GUI's sentence by the refusal's condition class. The classes and the messages
# the checks raise are unchanged (tests and records rely on them); the message is kept as the refusal's detail.
# ------------------------------------------------------------------------------------------------

#' The plain sentence of each refusal class the integration tabs can meet, in the order they are looked up (a refusal
#' carries its own class and `integration_refused`).
#' @noRd
INTEGRATION_GUI_REFUSALS <- c(
    integration_input_invalid = "refusal_input_invalid",
    integration_input_stale = "refusal_input_stale",
    integration_input_duplicate = "refusal_input_duplicate",
    integration_too_few_omics = "refusal_too_few_omics",
    integration_too_few_views = "refusal_too_few_views",
    integration_inputs_missing = "refusal_inputs_missing",
    integration_no_shared_samples = "refusal_no_shared",
    integration_group_conflict = "refusal_group_conflict",
    integration_view_empty = "refusal_view_empty",
    integration_input_matrix_unavailable = "refusal_input_matrix_unavailable",
    integration_parameters_invalid = "refusal_parameters_invalid",
    integration_mofa_unavailable = "refusal_mofa_unavailable",
    integration_sample_match_invalid = "refusal_sample_match_invalid")

#' A refusal in words: a sentence the tab wrote itself (`integration_gui_refusal`) as it is; any other refusal as the
#' plain sentence of its class, the check's own message as its detail.
#' @param condition The condition.
#' @param texts `integrationGuiTexts()`.
#' @return `list(sentence, detail)`; `detail` NULL for the tab's own sentences.
#' @noRd
integrationRefusalText <- function(condition, texts = integrationGuiTexts()) {
    if (inherits(condition, "integration_gui_refusal")) return(list(sentence = conditionMessage(condition), detail = NULL))
    known <- names(INTEGRATION_GUI_REFUSALS)[vapply(names(INTEGRATION_GUI_REFUSALS), function(class) inherits(condition, class),
        logical(1))]
    key <- if (length(known)) INTEGRATION_GUI_REFUSALS[[known[[1L]]]] else "refusal_other"
    list(sentence = texts[[key]], detail = conditionMessage(condition))
}

#' Refuse with one of the tab's own sentences (KIT-140): a plain condition, so that the sentence reaches the screen
#' exactly as written (no console formatting); classed `integration_gui_refusal` and `integration_refused`.
#' @param sentence The sentence.
#' @param class The refusal's own class.
#' @noRd
integrationGuiRefuse <- function(sentence, class) {
    stop(structure(class = c(class, "integration_gui_refusal", "integration_refused", "error", "condition"),
        list(message = sentence, call = NULL)))
}

#' A refusal on screen: the sentence, and the check's message under a collapsed "Technical detail" line.
#' @param sentence The sentence.
#' @param detail The check's message, or NULL.
#' @param class The box's class.
#' @noRd
integrationRefusalUi <- function(sentence, detail = NULL, class = "alert alert-danger", texts = integrationGuiTexts()) {
    shiny::div(class = class, shiny::p(sentence),
        if (!is.null(detail) && length(detail) && nzchar(detail)) {
            shiny::tags$details(shiny::tags$summary(texts[["technical_detail"]]), shiny::tags$code(detail))
        })
}

#' Save the confirmed sample table as the integration design (KIT-139): the aligned sample table, tab-separated, in the
#' integration results' inputs folder.
#' @param samples The aligned sample table (`integrationAlignViews()$samples`).
#' @param dir The folder.
#' @return The file's path.
#' @noRd
integrationSaveDesign <- function(samples, dir) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    path <- file.path(dir, "integration_design.tsv")
    vroom::vroom_write(samples, path, delim = "\t", na = "NA", progress = FALSE)
    path
}

# ------------------------------------------------------------------------------------------------
# The Model tab (KIT-140).
# ------------------------------------------------------------------------------------------------

#' The parameters the confirmed sample table fixes (KIT-139 shows their effect and confirms it): shown on the Model
#' tab as set, not offered again, so that the model is trained on the table the user confirmed.
#' @noRd
INTEGRATION_GUI_MATCHING_PARAMETERS <- c("technical_replicates", "internal_standards", "missing_samples")

#' Whether a parameter's default is provisional: the spec's `settled_by` names the workbook as the one that settles or
#' may override it (every parameter but `precision`, decided by INT-6, today).
#' @param entry One entry of `integrationMofaParameterSpec()`.
#' @noRd
integrationParameterProvisional <- function(entry) {
    grepl("provisional|may override", entry$settled_by)
}

#' A setting's label, its choices' labels and its help, from the wording table; a parameter or value without an entry
#' shows its own name.
#' @noRd
integrationParameterText <- function(texts, prefix, name, value = NULL) {
    key <- paste(c(prefix, name, value), collapse = "_")
    if (key %in% names(texts)) texts[[key]] else if (is.null(value)) name else as.character(value)
}

#' A parameter's input type (KIT-152): the spec entry's `type`, else read from the entry as the Model tab always did
#' (`choice` with choices, `flag` with a logical default, `whole` otherwise). The other types: `number` (any number) and
#' `factors` (`all` or factor names, KIT-063's `factors_shown`).
#' @param entry One entry of a parameter spec.
#' @noRd
integrationParameterType <- function(entry) {
    entry$type %||% (if (!is.null(entry$choices)) "choice" else if (is.logical(entry$default)) "flag" else "whole")
}

#' A tab's settings, built from its parameter spec (KIT-140; KIT-152 for any tab and a type per parameter): one input
#' per parameter, its default the spec's, so a default the workbook changes changes the screen with no edit here; a
#' choice as buttons, a flag as a tick box, a whole number or a number as a number box, factors as a list to pick from
#' (its choices filled by the tab once a model is shown). Each label carries "provisional default" or "decided".
#' @param ns The module's namespace function.
#' @param spec `integrationMofaParameterSpec()` or `integrationResultsParameterSpec()`.
#' @param exclude Parameters not offered (the Model tab: those the confirmed sample table fixed).
#' @param testid The prefix of each input's test id.
#' @param note The sentence above the inputs.
#' @return A tag list.
#' @noRd
integrationModelSettingsUi <- function(ns, spec = integrationMofaParameterSpec(), texts = integrationGuiTexts(),
    exclude = INTEGRATION_GUI_MATCHING_PARAMETERS, testid = "integration-model-param-", note = texts[["settings_provisional"]]) {
    inputs <- lapply(setdiff(names(spec), exclude), function(name) {
        entry <- spec[[name]]
        id <- ns(paste0("param_", name))
        mark <- texts[[if (integrationParameterProvisional(entry)) "mark_provisional" else "mark_decided"]]
        label <- shiny::tagList(integrationParameterText(texts, "param", name), " ",
            shiny::tags$small(class = "text-muted", sprintf("(%s)", mark)))
        input <- switch(integrationParameterType(entry),
            choice = shiny::radioButtons(id, label, choiceNames = lapply(entry$choices, function(value)
                integrationParameterText(texts, "choice", name, value)), choiceValues = entry$choices, selected = entry$default),
            flag = shiny::checkboxInput(id, label, value = isTRUE(entry$default)),
            number = shiny::numericInput(id, label, value = entry$default),
            factors = shiny::selectInput(id, label, choices = stats::setNames(entry$default,
                integrationParameterText(texts, "choice", name, entry$default)), selected = entry$default, multiple = TRUE),
            shiny::numericInput(id, label, value = entry$default, step = 1))
        help <- paste0("help_", name)
        shiny::div(`data-testid` = paste0(testid, name), input,
            if (help %in% names(texts)) shiny::helpText(texts[[help]]))
    })
    shiny::tagList(shiny::p(shiny::em(note)), inputs)
}

#' The settings the confirmed sample table fixed, in words.
#' @param parameters The confirmed parameters (`workflow_data$integration_alignment$parameters`).
#' @return A character vector.
#' @noRd
integrationModelMatchingSettings <- function(parameters, spec = integrationMofaParameterSpec(), texts = integrationGuiTexts()) {
    names <- intersect(INTEGRATION_GUI_MATCHING_PARAMETERS, names(spec))
    settings <- vapply(names, function(name) sprintf(texts[["setting_fixed"]], integrationParameterText(texts, "param", name),
        integrationParameterText(texts, "choice", name, parameters[[name]])), character(1), USE.NAMES = FALSE)
    # KIT-063 item A: a confirmed table with the user's sample match says so.
    match <- parameters[["sample_match_table", exact = TRUE]]
    if (!is.null(match)) settings <- c(settings, sprintf(texts[["setting_sample_match"]], length(match$run)))
    settings
}

#' A tab's settings as its inputs give them (KIT-140; KIT-152 for any tab), in the spec's forms: a choice as text, a
#' flag as TRUE or FALSE, a whole number as an integer, a number as a double, factors as `all` or the names picked (`all`
#' among them, or none picked, is `all`). The number of most variable features is set only when that selection is
#' chosen. Refused in words (`integration_parameters_invalid`): a number box left empty, or not a whole number where one
#' is needed.
#' @param values The inputs' values, named by parameter.
#' @param spec,exclude As `integrationModelSettingsUi()`.
#' @return The parameters, named.
#' @noRd
integrationModelSettings <- function(values, spec = integrationMofaParameterSpec(), texts = integrationGuiTexts(),
    exclude = INTEGRATION_GUI_MATCHING_PARAMETERS) {
    refuse <- function(name) integrationGuiRefuse(sprintf(texts[["setting_not_whole"]], integrationParameterText(texts, "param", name)),
        "integration_parameters_invalid")
    out <- list()
    for (name in setdiff(names(spec), exclude)) {
        entry <- spec[[name]]
        value <- values[[name]]
        type <- integrationParameterType(entry)
        out[name] <- list(if (type == "choice") {
            if (is.null(value)) entry$default else as.character(value)
        } else if (type == "flag") {
            if (is.null(value)) entry$default else isTRUE(value)
        } else if (type == "factors") {
            if (!length(value) || "all" %in% value) entry$default else as.character(value)
        } else if (is.null(value)) {
            # An input not yet drawn takes the table's default, as the screen shows it.
            entry$default
        } else if (length(value) != 1L || is.na(value)) {
            if (type == "number") NA_real_ else NA_integer_
        } else if (type == "number") {
            if (!is.numeric(value)) integrationGuiRefuse(sprintf(texts[["setting_not_number"]],
                integrationParameterText(texts, "param", name)), "integration_parameters_invalid")
            as.numeric(value)
        } else {
            if (!is.numeric(value) || value != round(value)) refuse(name)
            as.integer(value)
        })
    }
    # KIT-151 (INT-12): the number kept and the observed share are set only when features are selected.
    unused <- character()
    if ("feature_selection" %in% names(out) && identical(out$feature_selection, "none")) {
        unused <- intersect(c("features_per_view", "min_observed_fraction"), names(out))
        out[unused] <- list(features_per_view = NA_integer_, min_observed_fraction = NA_real_)[unused]
    }
    for (name in setdiff(names(out), unused)) {
        if (is.double(out[[name]]) && length(out[[name]]) == 1L && is.na(out[[name]])) {
            integrationGuiRefuse(sprintf(texts[["setting_not_number"]], integrationParameterText(texts, "param", name)),
                "integration_parameters_invalid")
        }
        if (is.integer(out[[name]]) && length(out[[name]]) == 1L && is.na(out[[name]])) refuse(name)
    }
    out
}

#' The model's parameters for a run: the confirmed table's (inputs, data sets and the matching rules) with the tab's
#' method settings.
#' @noRd
integrationModelParameters <- function(confirmed, settings) {
    parameters <- confirmed
    parameters[names(settings)] <- settings
    parameters
}

#' The checks a model run makes before any training (KIT-140; `runIntegrationModel`'s, made in the app's process before
#' a run is handed to the managed worker): the inputs imported and at the head, every input unchanged on disk, the
#' parameters resolved against the inputs, MOFA2 and its Python available. Refused in the model's words and classes.
#' @param manager The integration store.
#' @param parameters The run's parameters.
#' @param probe The capability probe.
#' @return `list(source, parameters, request)`: the inputs state's reference, the resolved parameters (accepted on the
#'   head with the model), and the worker request's parameters: the resolved ones with a sample match kept in the store
#'   and named by its table's reference (KIT-148; `integrationStoreSampleMatch`), so that a match of any size stays
#'   within the request's 64 KiB.
#' @noRd
integrationModelPrecheck <- function(manager, parameters, probe = integrationMofaCapability) {
    withWorkflowStateScope(manager, function() {
        source <- getWorkflowStateReference(manager, INTEGRATION_INPUTS_STATE)
        if (is.null(source) || !identical(manager$current_state, INTEGRATION_INPUTS_STATE)) {
            integrationRefuse("the integration inputs are not imported; import them first", "integration_inputs_missing")
        }
        object <- readWorkflowExecutionComponents(manager, source, fields = "inputs")
        integrationAssertInputsCurrent(object)
        resolved <- integrationResolveMofaParameters(object, parameters)
        capability <- probe()
        if (!isTRUE(capability$available)) {
            integrationRefuse(sprintf("the integration model cannot be trained here: %s", capability$reason),
                "integration_mofa_unavailable")
        }
        request <- resolved
        if (!is.null(resolved$sample_match_table)) {
            request$sample_match_table <- integrationStoreSampleMatch(manager, resolved$sample_match_table)
        }
        list(source = source, parameters = resolved, request = request)
    })
}

#' What the Model tab shows of a trained model (KIT-140; only this, anything richer is KIT-063's): whether training
#' settled (rounds against the maximum), the factors kept against those asked for, the samples, the variance explained
#' per data set as a compact table and the model file's checksum.
#' @param value The model's product (`integrationModelProduct`).
#' @param parameters The run's resolved parameters.
#' @return `list(converged, convergence, factors, samples, variance, saved)`.
#' @noRd
integrationModelSummary <- function(value, parameters, texts = integrationGuiTexts()) {
    iterations <- as.integer(value$elbo$iterations)
    maxiter <- as.integer(parameters$maxiter)
    # KIT-151 (INT-15): how training stopped, as the product keeps it (a product without it: fewer rounds than the limit).
    stopped_by <- value$elbo$stopped_by %||% NA_character_
    converged <- if (!is.na(stopped_by)) identical(stopped_by, "convergence_rule") else isTRUE(iterations < maxiter)
    total <- value$variance_total
    features <- value$view_summary$features_model[match(total$view, value$view_summary$view)]
    variance <- data.frame(total$view, as.integer(features), sprintf("%.2f", total$r2_percent), stringsAsFactors = FALSE)
    names(variance) <- unname(texts[c("col_view", "col_model_features", "col_variance")])
    sha256 <- value$model_file$sha256
    list(converged = converged,
        convergence = if (converged) sprintf(texts[["result_converged"]], iterations, maxiter) else
            sprintf(texts[["result_not_converged"]], maxiter),
        factors = sprintf(texts[["result_factors"]], ncol(value$factors) - 1L, as.integer(parameters$num_factors)),
        samples = sprintf(texts[["result_samples"]], nrow(value$samples), as.integer(value$shared_samples)),
        variance = variance,
        saved = if (is.character(sha256) && length(sha256) == 1L && !is.na(sha256)) sprintf(texts[["result_saved"]], sha256))
}
