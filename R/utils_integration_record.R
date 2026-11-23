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

# KIT-142: the integration analysis record, its replay and its report. The record is the omics' analysis record
# (`buildWorkflowAnalysisRecord`, derived from the integration store's head); for an integration store its model
# analysis carries one more block, `integration`, with what the integration adds: every input (omic, data sets,
# manifest checksum, matrix role, correction), the confirmed sample match (with an uploaded match table), every model
# setting with its status (decided, provisional default, set by the user, from the data), the software versions, the
# model's identity and checksum, the saved results' settings and files, MOFA2's warnings and the claim limit (KIT-060
# D10). Replay re-imports the recorded inputs (each checked by its checksum) into a new store, trains the model with
# the recorded settings and compares it with the recorded model numerically (INT-6). The study parameters, the methods
# paragraph and the report read every value from the record; nothing here is hard-coded. One wording table,
# `integrationRecordTexts()`.

#' Every string of the integration record, its study parameters, methods paragraph, report and Record tab (KIT-142):
#' one table, drafts the user sets the final wording of (docs/qa/workflow-execution/kit/kit-138/wording.md, section
#' KIT-142). `%s`/`%d` are filled in the order the comment beside each entry gives.
#' @return A named character vector.
#' @noRd
integrationRecordTexts <- function() {
    c(
        # Status of a setting in the record (and its legend).
        status_decided = "decided",
        status_provisional_default = "provisional default",
        status_user_set = "set by the user",
        status_input = "from the data",
        # The study parameters' and report's sections.
        report_title = "Multi-omics integration report",
        report_intro = paste("This report was written by MultiScholaR from the integration's analysis record. Every value",
            "below is read from the record."),
        section_methods = "Methods (draft)",
        methods_draft_note = paste("A draft for the methods section, written from the record: check the wording before",
            "using it."),
        section_inputs = "Data sets combined",
        section_matching = "Sample matching",
        section_samples = "Biological samples in the model",
        section_settings = "Settings",
        section_model = "The trained model",
        section_variance = "Variance explained in each data set's selected features (all factors together)",
        section_results = "Saved results",
        section_notes = "Notes from the model",
        section_legend = "What the status of a setting means",
        none_recorded = "None recorded.",
        # Table headers.
        col_omic = "Omic",
        col_datasets_in_model = "Data sets in the model",
        col_values = "Values",
        col_correction = "Correction for unwanted variation (RUV)",
        col_group_corrected = "Corrected with the groups as replicate sets",
        col_checksum = "Checksum of the saved data (SHA-256)",
        col_field = "Item",
        col_value = "Value",
        col_setting = "Setting",
        col_name = "Internal name",
        col_status = "Status",
        col_settled_by = "Settled by",
        col_applies_to = "Applies to",
        col_file = "File",
        col_view = "Data set",
        col_variance = "Variance explained in the selected features (%)",
        col_note = "Note",
        col_status_meaning = "Meaning",
        legend_decided = "a decision recorded in the methods ledger; not expected to change",
        legend_provisional_default = "the current default; the integration methods may still change it",
        legend_user_set = "changed by the user from the default",
        legend_input = "taken from the imported data",
        applies_model = "model",
        applies_results = "results",
        yes = "yes",
        no = "no",
        not_recorded = "not recorded",
        values_corrected = "corrected for unwanted variation",
        values_normalised = "normalised, not corrected",
        # KIT-153: the omic saved corrected values too; the model used the normalised ones.
        values_uncorrected_available = "normalised, not corrected (the omic's corrected values were not used)",
        # The sample-matching facts.
        match_rule = "How samples were matched",
        match_rule_saved = "by the biological-sample column each omic saved",
        # %d: runs.
        match_rule_user = "by the biological-sample column each omic saved, corrected by an uploaded table (%d runs moved)",
        match_table_checksum = "Checksum of the uploaded match (SHA-256)",
        match_samples = "Biological samples in the model",
        match_shared = "Measured in every data set",
        match_confirmed_at = "Sample table confirmed at",
        match_design_file = "Sample table saved as",
        match_design_checksum = "Checksum of the saved sample table (SHA-256)",
        match_file = "Uploaded match saved as",
        # The model facts.
        model_mofa2 = "MOFA2",
        model_mofapy2 = "mofapy2",
        model_python = "Python",
        model_r = "R",
        model_reticulate = "reticulate",
        model_precision = "Numeric precision",
        model_iterations = "Training rounds",
        model_settled = "Training settled (stopped by its convergence rule, not at the limit)",
        model_factors = "Factors (kept, asked for)",
        model_checksum = "Checksum of the model file (SHA-256)",
        model_bytes = "Size of the model file (bytes)",
        model_entity = "Model's identifier in the project",
        model_identity = "Checksum of the model's identity (data and settings)",
        # KIT-151 (INT-14).
        model_seeds = "Models trained (one per seed)",
        # %d: models; %s: the seeds.
        model_seeds_value = "%d (seeds %s)",
        # %d: kept, asked for.
        model_factors_value = "%d of %d",
        # The few-samples note (D7). %d: samples in every data set, the threshold.
        note_few_samples = "Only %d biological samples are measured in every data set (fewer than %d).",
        # The methods paragraph (KIT-142), each sentence from the record's values.
        # %d: data sets, omics; %s: the list of data sets.
        methods_inputs = "%d data sets from %d omics were combined: %s.",
        # %s: the omic, its data sets, its values in words.
        methods_input_item = "%s (%s; %s)",
        # %s: the approach in words.
        methods_input_corrected = "corrected for unwanted variation, %s",
        methods_input_normalised = "normalised, not corrected for unwanted variation",
        # KIT-153. %s: the omic's correction in words (not used).
        methods_input_uncorrected_available = paste("normalised, before correction for unwanted variation; the omic's",
            "corrected values (%s) were not used"),
        # %s: MOFA2, mofapy2, Python, R versions, the precision.
        methods_software = paste("Shared factors were learned with multi-omics factor analysis, MOFA2 %s (mofapy2 %s,",
            "Python %s, R %s), in %s precision."),
        methods_matching_saved = "Samples were matched across omics by the biological-sample column each omic saved.",
        # %d: runs; %s: the checksum.
        methods_matching_user = paste("Samples were matched across omics by the biological-sample column each omic saved,",
            "corrected by a user-supplied match table that assigned %d runs to another biological sample (SHA-256 %s)."),
        methods_replicates_average = "Runs of the same biological sample were averaged (the mean of the observed values).",
        # %d: internal standards removed.
        methods_standards_drop = "Internal standards were removed (%d features).",
        methods_standards_keep = "Internal standards were kept.",
        methods_missing_keep_missing = "Biological samples missing from a data set were kept, with that data set's values missing.",
        methods_missing_drop_incomplete = "Only the biological samples measured in every data set were kept.",
        methods_sparse = "Features measured in fewer than two biological samples were left out.",
        # KIT-151 (INT-12). %d: the fewest samples, the samples; %s: the share.
        methods_sparse_selected = "Features measured in fewer than %d of the %d biological samples (a share of %s) were left out.",
        methods_selection_none = "All remaining features were used",
        # %d: features per data set, twice.
        methods_selection_top_variance = paste("In each data set with more than %d features, the %d with the largest",
            "variance were used (a smaller data set kept all its features)"),
        # KIT-151 (INT-12). %d: features per data set, twice.
        methods_selection_variance_above_trend = paste("In each data set with more than %d features, the %d whose",
            "variance lies furthest above the data set's trend of variance with mean abundance (a lowess fit of each",
            "feature's variance on its mean, span 0.3) were used, a feature selection that also makes the data sets more",
            "equal in size (a smaller data set kept all its features)"),
        methods_scaling_none = "without per-feature scaling.",
        # KIT-151 (INT-24).
        methods_scaling_feature_zscore = paste("each scaled to mean 0 and standard deviation 1 after the selection (which",
            "used the unscaled values); in simulation, with few samples and missing values this scaling improved the",
            "recovery of weak factors, and otherwise made no practical difference."),
        # %d: factors asked for; %s: convergence mode; %d: maximum rounds, seed; %s: views scaled, features centred.
        methods_model = paste("The model was asked for %d factors (convergence mode %s, at most %d training rounds, seed",
            "%d; each data set divided by its overall standard deviation: %s; features centred: %s). Variance explained",
            "refers to the selected features of each data set, not to every feature measured."),
        # %d: rounds, factors kept, samples, samples in every data set.
        methods_trained_settled = paste("Training settled after %d rounds, stopped by its convergence rule, and kept %d",
            "factors; %d biological samples entered the model, %d of them measured in every data set."),
        # %d: maximum rounds, factors kept, samples, samples in every data set.
        methods_trained_not_settled = paste("Training did not settle: it stopped at the limit of %d rounds before its",
            "convergence rule was met, and kept %d factors; %d biological samples entered the model, %d of them measured in",
            "every data set."),
        # KIT-151 (INT-14). %d: models, first seed, last seed, the reported model's seed; %s: each factor's lowest agreement.
        methods_seeds = paste("The model was trained %d times, from seeds %d to %d (the declared seed and the next whole",
            "numbers), and the model reported is that of seed %d. Each factor's agreement across seeds is the absolute",
            "correlation of its values with the most similar factor of each other model (factors paired one to one, the",
            "most similar pair first); the lowest agreement per factor was %s."),
        # %d: the seed.
        methods_seeds_single = "The model was trained once, from seed %d; the factors' agreement across seeds was not measured.",
        # KIT-151 (INT-17). %s: the variables set.
        methods_pin_applied = paste("Training ran with the maths library on one thread and one fixed code path (%s), under",
            "which models were bit-identical across the processor models tested (three AMD models; an Intel processor not",
            "yet tested); the same data, settings and seed give the identical model on the same software."),
        # %s: the reason.
        methods_pin_not_applied = paste("The maths library could not be fixed to one thread and one code path for this",
            "training (%s); the model is identical on the same computer and software, agrees to about 1e-12 on other",
            "processors, and a factor near the detection limit can differ between processor models."),
        # %s: the adjustment's name.
        methods_association = paste("Each factor was related to the experimental groups after training: a one-way analysis",
            "of variance of the factor's values on the group, with the share of the factor's variance that lies between",
            "the groups (R-squared) as the effect size, and p-values adjusted across all the factors of the model (%s)."),
        # KIT-151 (plan item E). %s: the adjustment's name.
        methods_association_kruskal_wallis = paste("Each factor was related to the experimental groups after training: a",
            "Kruskal-Wallis rank test of the factor's values across the groups, with epsilon-squared (H divided by the",
            "number of samples minus one) as the effect size, and p-values adjusted across all the factors of the model",
            "(%s)."),
        # %s: the data sets.
        methods_claim_limit = paste("In %s, unwanted variation had been removed with the groups as replicate sets, which",
            "keeps the differences between groups by construction; p-values for the association of the factors with the",
            "groups are therefore not reported, and the factors are described, not tested."),
        methods_claim_none = paste("The factors were learned without the groups; their association with the groups",
            "describes the model."),
        # KIT-153: the model used the normalised values before correction for unwanted variation.
        methods_claim_uncorrected = paste("The data sets entered the model as their normalised values before correction for",
            "unwanted variation (RUV), so unwanted variation that is correlated with the groups was not removed and can",
            "contribute to a factor's association with the groups."),
        # %s: the settings.
        methods_provisional = "Settings that were provisional defaults at the time of the analysis: %s.",
        # %s: the settings.
        methods_user_set = "Settings changed from their defaults: %s.",
        # KIT-152: the further figures, in the methods paragraph (audits SCI-15 to SCI-17). %d: the features listed.
        methods_top_features = paste("Within each factor and data set, features were ranked by the absolute value of their",
            "weight and the %d highest were listed; since the sign of a factor is arbitrary, a feature's direction has",
            "meaning only against the other features of the same factor. These lists are exploratory: their order was not",
            "tested for stability."),
        # %s: the scale of each data set's values.
        methods_feature_values = paste("Scatter plots show these features' values as given to the model (the runs of a",
            "biological sample averaged, after feature selection and any per-feature scaling, before MOFA2 centres and",
            "scales each data set; %s) against the factor's values, with missing values left out; the features were chosen",
            "from the same model."),
        methods_heatmap_centre = paste("Heatmaps show the same values with each feature centred on its mean, for display",
            "only, the samples in order of the factor's values."),
        methods_heatmap_centre_scale = paste("Heatmaps show the same values with each feature centred on its mean and scaled",
            "to standard deviation 1, for display only, the samples in order of the factor's values."),
        # %d: factors; %s: uwot's version; %d: neighbours; %s: minimum distance; %d: training rounds; %s: the start;
        # %d: seed.
        methods_umap = paste("A UMAP map of the %d factors was drawn with uwot %s (%d neighbours, minimum distance %s, %d",
            "training rounds, %s start, one thread, seed %d)."),
        # %d: factors; %s: Rtsne's version, the perplexity; %d: iterations, seed.
        methods_tsne = paste("A t-SNE map of the %d factors was drawn with Rtsne %s (perplexity %s, exact computation with",
            "theta 0, %d iterations, started from the first two principal components scaled to standard deviation 0.0001,",
            "one thread, seed %d)."),
        # %d: biological samples.
        methods_maps_caution = paste("These maps are pictures of the factors at %d biological samples and are unstable at",
            "this number: the distances and cluster sizes on them are not interpretable, they are not evidence of groups,",
            "and they are reproduced exactly only on the same computer."),
        # KIT-152 correction: what each map drawn depends on.
        methods_maps_umap = "The UMAP map changes with the random seed and the number of neighbours.",
        methods_maps_tsne = paste("The t-SNE map, started from the principal components, does not depend on the random seed",
            "but changes with the perplexity."),
        # %s: the data sets.
        methods_maps_group_corrected = paste("In %s, unwanted variation had been removed with the groups as replicate sets,",
            "so groups separate on these maps by construction."),
        # KIT-152: the further figures in the record's tables, study parameters and report.
        section_outputs = "Further figures",
        col_output = "Figure",
        col_shows = "What it shows",
        col_settings = "Settings",
        out_pairs = "Factors against each other",
        out_pairs_shows = "the factor values of each biological sample, two factors per panel, coloured by group",
        out_top = "Features with the largest weights",
        # %d: the features listed.
        out_top_shows = paste("the %d features per factor and data set with the largest weight in absolute value; a",
            "factor's sign is arbitrary; exploratory, not tested for stability"),
        out_scatter = "Top features against the factor",
        # %s: the scale of each data set's values.
        out_scatter_shows = "the values given to the model (%s) against the factor's values; missing values left out",
        out_heatmap = "Heatmap of the top features",
        out_heatmap_shows_centre = "the same values, each feature centred on its mean, for display only; samples in order of the factor",
        out_heatmap_shows_centre_scale = paste("the same values, each feature centred on its mean and scaled to standard",
            "deviation 1, for display only; samples in order of the factor"),
        out_map_umap = "UMAP map of the factors",
        out_map_tsne = "t-SNE map of the factors",
        # %d: factors, biological samples.
        out_map_shows = "the %d factors of %d biological samples in two dimensions; a picture, not evidence of groups",
        # %s: version; %d: neighbours; %s: minimum distance; %d: training rounds; %s: the start; %d: seed.
        out_umap_settings = "uwot %s; %d neighbours; minimum distance %s; %d training rounds; %s start; one thread; seed %d",
        # %s: version, perplexity; %d: iterations, seed.
        out_tsne_settings = paste("Rtsne %s; perplexity %s; exact (theta 0); %d iterations; started from the first two",
            "principal components; one thread; seed %d"),
        out_caution = "Caution about the maps",
        # %s: the data set, then its scale in words.
        scale_item = "%s: %s",
        scale_log2 = "log2",
        scale_log2_cpm = "log2 counts per million",
        scale_feature_zscore = "scaled per feature",
        scale_unknown = "the scale the omic saved",
        # The Record tab (its title and step are `integrationGuiTexts()`'s `tab_record`, `step_record`).
        record_heading = "Record and report",
        record_intro = paste("Save the integration's analysis record (what was combined, how samples were matched, every",
            "setting, the software versions and the model's checksum), its study parameters and a report with a draft",
            "methods paragraph."),
        record_not_ready = "Train a model at the Model tab and save its results at the Results tab first.",
        record_button = "Save the record and the report",
        record_working = "Writing the record and the report.",
        # %s: the time; %d: files; %s: the folder.
        record_done = "Saved at %s: %d files in %s.",
        # %s: the reason.
        record_refused = "The record was not saved. %s",
        record_methods_heading = "Methods paragraph (draft, from the record)",
        record_no_model = "There is no trained model to record. Train a model at the Model tab first.",
        # Replay refusals and warnings, in plain words.
        replay_no_model = "This record holds no integration model, so there is nothing to repeat.",
        # %s: the omic, then the file.
        replay_input_missing = paste("The saved data of %s is no longer at %s, so the analysis cannot be repeated. Give",
            "its new location, or save it again at that omic's Summary tab."),
        # %s: the omic.
        replay_input_changed = paste("The saved data of %s has changed since this analysis was recorded (its checksum",
            "differs), so the analysis cannot be repeated on it."),
        replay_mofa_unavailable = paste("The analysis cannot be repeated on this computer: the MOFA2 software, or the",
            "Python it uses, is not installed."),
        # %s: the differences.
        replay_versions_differ = paste("The analysis was recorded with other software versions (%s); it is repeated and",
            "compared with the recorded model numerically."),
        # %s: the software, recorded, here.
        replay_version_item = "%s %s, here %s",
        # KIT-141: the enrichment of the factors, in the record's tables, study parameters, report and methods paragraph.
        section_enrichment = "Enrichment of the factors",
        section_enrichment_terms = "Enriched terms and pathways",
        # The enrichment sheets' first header. Without it the header was NA, and one NA string empties the workbook's
        # whole shared-string table as openxlsx writes it (run-tests 44799395381).
        col_factor = "Factor",
        col_enrich_database = "Database",
        col_enrich_release = "Release",
        col_enrich_ranking = "How features were ranked",
        col_enrich_direction = "Direction of the factor",
        col_enrich_features = "Features (measured, ranked, matched)",
        col_enrich_tests = "Terms or pathways in the table",
        col_enrich_family = "Corrected together for multiple testing",
        col_enrich_term = "Term or pathway",
        col_enrich_description = "Description",
        col_enrich_category = "Category",
        col_enrich_score = "Enrichment score (NES) or effect (camera)",
        col_enrich_fdr = "False discovery rate (adjusted p-value)",
        # %s: STRING's version, its address.
        enrich_string_release = "STRING %s (%s)",
        # %s: the correlation; %d: the smallest number of samples.
        enrich_ranking_value = "%s correlation with the factor, a feature needing values in at least %d biological samples",
        # %s: the feature, its data set; %d: +1 or -1.
        enrich_direction_value = "oriented so that the largest weight (%s in %s) is positive (factor multiplied by %d)",
        # %d: measured, ranked, matched.
        enrich_counts_value = "%d measured, %d ranked, %d matched to the database",
        # %s: the correlation (pearson or spearman); %d: the smallest number of samples.
        methods_enrich_ranking = paste("For enrichment, every measured feature of a data set (not only those given to the",
            "model) was ranked by the %s correlation of its values, the runs of each biological sample averaged, with the",
            "factor's values over the biological samples where both were observed; a feature needed values in at least %d",
            "biological samples. Since the sign of a factor is arbitrary, each factor was first oriented so that its largest",
            "weight over all data sets was positive."),
        # %s: factor, data set; %d: identifiers; %s: STRING's version, the species, the date.
        methods_enrich_string = paste("%s in %s: the complete ranked list (%d identifiers) was submitted to STRING's",
            "ranked enrichment (values/ranks), STRING %s, species %s, on %s."),
        # %s: factor, data set, database; %d: compounds; %s: release, organism; %d: smallest, largest pathway size,
        # pathways tested, seed.
        methods_enrich_pathway = paste("%s in %s: %s pathways were tested by ranked gene-set enrichment (clusterProfiler",
            "GSEA) on %d compounds matched by the chosen identifier column (several features on one compound: the one",
            "with the largest absolute correlation kept), %s, organism %s, pathways with %d to %d measured compounds",
            "(%d pathways), every pathway kept, no random tie-breaking, seed %d."),
        methods_enrich_family = paste("P-values were adjusted with the Benjamini-Hochberg method within one factor, one",
            "data set and one database (STRING: within each of its term categories); they were not adjusted across",
            "factors, data sets or databases. The results depend on external databases at the release stated."),
        # KIT-183 (VAL-51, VAL-53; proposed wording, VAL-41). The numbers in brackets are the simulation's measured rates
        # (docs/qa/workflow-execution/kit/kit-183/results.md), written in as constants.
        col_enrich_strength = "Strength of the factor in this data set",
        enrich_factor_strength = integrationGuiTexts()[["enrich_factor_strength"]],
        enrich_string_test_note = "STRING runs its own ranked test; the test setting does not apply.",
        # %d: biological samples, features used.
        enrich_ranking_value_camera = paste("regression of each feature on the factor over %d biological samples",
            "(features measured in all of them: %d); correlation-adjusted test (camera, correlation estimated)"),
        enrich_ranking_value_camera_pr = paste("regression of each feature on the factor over %d biological samples",
            "(features measured in all of them: %d); pre-ranked test with an assumed correlation of 0.01 (cameraPR)"),
        methods_test_name_camera = paste("a competitive gene-set test that allows for the correlation between a term's",
            "features (limma camera, the correlation estimated for each term from the residuals; Wu and Smyth 2012),",
            "applied to the regression of each feature's values on the oriented factor"),
        methods_test_name_camera_pr = paste("a competitive gene-set test with the correlation between a term's features",
            "assumed to be 0.01 (limma cameraPR on the moderated t of the regression of each feature's values on the",
            "oriented factor)"),
        methods_test_name_gsea = paste("ranked gene-set enrichment (clusterProfiler GSEA) on each feature's correlation",
            "with the oriented factor"),
        # %d: biological samples, identifiers used, matched identifiers.
        methods_test_settings_camera = paste("the correlation estimated for each term (negative estimates kept),",
            "moderated variances with an abundance trend, an intercept and no other covariate, features measured in",
            "every one of %d biological samples (%d of %d matched identifiers; several features on one identifier: the",
            "one with the largest mean value kept)"),
        methods_test_settings_camera_pr = paste("the correlation assumed to be 0.01, moderated variances with an",
            "abundance trend, an intercept and no other covariate, features measured in every one of %d biological",
            "samples (%d of %d matched identifiers; several features on one identifier: the one with the largest mean",
            "value kept)"),
        # %d: the seed.
        methods_test_settings_gsea = "no random tie-breaking, seed %d",
        # %s: factor, data set, database, the test's name, its settings, release, organism; %d: smallest, largest
        # pathway size, pathways tested.
        methods_enrich_camera = paste("%s in %s: %s pathways were tested with %s, %s, %s, organism %s, pathways with %d",
            "to %d of these identifiers (%d pathways), every pathway kept."),
        methods_enrich_camera_pr = paste("%s in %s: %s pathways were tested with %s, %s, %s, organism %s, pathways with",
            "%d to %d of these identifiers (%d pathways), every pathway kept."),
        methods_enrich_limits_tests = paste("The factor was learned from these same data, which neither test allows for,",
            "so the enrichment results are exploratory: features the model was trained on shaped the factor and tend to be",
            "more strongly associated with it than other features even without a biological link, the more so for factors",
            "weaker than the main effect. In a check on simulated data shaped like this study (12 biological samples, 800",
            "simulated studies), the correlation-adjusted test called sets with no link to the factor at the 5 % level in",
            "0.4 to 6.6 % of cases, depending on the kind of set, on the main factor and on values unrelated to the data",
            "(largest: sets of low-abundance proteins against values unrelated to the data, 6.6 %, 95 % interval 5.8 to",
            "7.4 %, above the nominal 5 %), within the programme's tolerance of an upper 95 % limit of at most 7.5 %; sets",
            "made only of the most variable features, those the model selects, were called in up to 10.1 % (95 % interval",
            "9.5 to 10.7 %) even against values unrelated to the data. On the main factor, sets with no link whose features",
            "were correlated at 0.3 were called by the pre-ranked test in 13.5 % (sets of 5 to 10 metabolites) to 43.0 %",
            "(sets of 10 to 100 genes; 95 % interval 41.9 to 44.1 %) of cases. The correlation-adjusted test uses only features measured in every",
            "biological sample, which favours abundant features and leaves out features that fall below detection in",
            "some samples, often those that change most. Both tests use each data set's normalised values before RUV",
            "correction, the values the model is trained on by default; for a model trained on RUV-corrected values,",
            "variation the correction removed stays in the values tested. P-values were adjusted within each factor, data",
            "set and database (and Gene Ontology aspect), not across them.")
    )
}

#' The record's sheets for an integration record (KIT-142), in the order the workbook writes them.
#' @noRd
INTEGRATION_RECORD_SHEETS <- c("integration_inputs", "integration_matching", "integration_samples",
    "integration_settings", "integration_model", "integration_variance", "integration_results", "integration_outputs",
    "integration_notes", "integration_enrichment", "integration_enrichment_terms")

# ------------------------------------------------------------------------------------------------
# Small readers that take the record as built (lists, data frames) or as read back from JSON
# (`jsonlite::read_json(simplifyVector = FALSE)`: lists of rows, scalars as length-one vectors).
# ------------------------------------------------------------------------------------------------

#' One value of the record: NA when absent, null or empty; the first element otherwise.
#' @noRd
integrationRecordScalar <- function(value) {
    value <- unlist(value, use.names = FALSE)
    if (!length(value)) NA else value[[1L]]
}

#' A table of the record as a data frame: a data frame as it is; a list of rows (JSON) bound by the union of their
#' names, an absent or null field read as NA.
#' @noRd
integrationRecordFrame <- function(x) {
    if (is.data.frame(x)) return(x)
    if (!is.list(x) || !length(x)) return(data.frame())
    columns <- unique(unlist(lapply(x, names)))
    out <- lapply(stats::setNames(columns, columns), \(column) unlist(lapply(x, \(row) {
        value <- row[[column]]
        if (is.null(value) || !length(value)) NA else unlist(value)[[1L]]
    })))
    as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}

#' A recorded file's name (its folder left out), or the words for "not recorded".
#' @noRd
integrationRecordBasename <- function(path, texts = integrationRecordTexts()) {
    path <- integrationRecordScalar(path)
    if (is.na(path) || !nzchar(as.character(path))) texts[["not_recorded"]] else basename(as.character(path))
}

#' One entry of a wording table, or `fallback` when the table has no such entry.
#' @noRd
integrationTextOr <- function(texts, key, fallback) if (length(key) == 1L && key %in% names(texts)) texts[[key]] else fallback

#' The integration model's analysis of a record, or NULL.
#' @noRd
integrationRecordAnalysisOf <- function(record) record$analyses[[INTEGRATION_MODEL_OPERATION]]

#' One recorded model parameter's value (decoded), or NULL.
#' @param parameters The recorded parameters (a list of `list(name, value, ...)`).
#' @noRd
integrationRecordParameter <- function(parameters, name) {
    hit <- Filter(\(parameter) identical(integrationRecordScalar(parameter$name), name), parameters %||% list())
    if (!length(hit)) return(NULL)
    decodeWorkflowRecordValue(hit[[1L]]$value)
}

# ------------------------------------------------------------------------------------------------
# What the record holds (KIT-142 item 1).
# ------------------------------------------------------------------------------------------------

#' The status of one setting in the record (KIT-142): `input` for what the data give (the inputs' checksums, the data
#' sets), `user_set` for the uploaded sample match and for any value that differs from its default, else `decided` or
#' `provisional_default` as the parameter table marks its default (`integrationParameterProvisional`, as the Model tab
#' shows it).
#' @param name The parameter.
#' @param value Its recorded value.
#' @param spec `integrationMofaParameterSpec()` or `integrationResultsParameterSpec()`.
#' @param factors For the results' `factors_shown`: the model's factors (every factor shown is the default `all`).
#' @return One of `input`, `user_set`, `decided`, `provisional_default`.
#' @noRd
integrationParameterStatus <- function(name, value, spec, factors = NULL) {
    if (name %in% c("input_manifests", "views")) return("input")
    if (name %in% c("sample_match", "sample_match_table")) return("user_set")
    entry <- spec[[name]]
    if (is.null(entry)) return("user_set")
    default <- entry$default
    if (identical(name, "factors_shown") && identical(default, "all") && !is.null(factors) &&
        identical(as.character(value), as.character(factors))) value <- "all"
    same <- length(value) == length(default) && identical(is.na(value), is.na(default)) &&
        identical(as.character(value[!is.na(value)]), as.character(default[!is.na(default)]))
    if (!same) return("user_set")
    if (integrationParameterProvisional(entry)) "provisional_default" else "decided"
}

#' The confirmation of the sample table, as the record keeps it (KIT-142): when, the saved sample table and the saved
#' match, each with its checksum.
#' @param confirmed The Sample Matching tab's confirmation (`workflow_data$integration_alignment`), or NULL.
#' @return NULL, or `list(confirmed_at, design_file, design_sha256, sample_match_file, sample_match_sha256, samples)`.
#' @noRd
integrationConfirmationRecord <- function(confirmed) {
    if (is.null(confirmed)) return(NULL)
    checksum <- function(path) if (is.character(path) && length(path) == 1L && file.exists(path))
        digest::digest(file = path, algo = "sha256") else NA_character_
    time <- confirmed$confirmed_at
    list(confirmed_at = if (inherits(time, "POSIXt")) format(time, "%Y-%m-%dT%H:%M:%S%z") else NA_character_,
        design_file = confirmed$file %||% NA_character_, design_sha256 = checksum(confirmed$file),
        sample_match_file = confirmed$sample_match_file %||% NA_character_,
        sample_match_sha256 = checksum(confirmed$sample_match_file),
        samples = if (is.data.frame(confirmed$samples)) nrow(confirmed$samples) else NA_integer_)
}

#' Keep a saved Results publication with the model it was made from (KIT-142): the publication's label and identity,
#' the Results settings, the time, every table and figure with its checksum, the factors' relation to the groups and
#' the claim limit's facts, under the model's entry of `snapshot$exports` (the record reads them there). A publication
#' made from a model that is no longer the latest is not kept. KIT-152 (audit ENG-15): with them what the further
#' figures show and every map's settings, as plain lists (`outputs`), and the top features' values and each map's
#' coordinates as tables of the integration store, only their references on the head (`output_tables`), for the replay's
#' comparison.
#' @param manager The integration owner.
#' @param model_reference The model product's reference.
#' @param results `integrationResults()`'s value.
#' @param files The saved files by kind (`tables`, `figures`: their paths).
#' @param identity The publication's identity.
#' @return TRUE when kept, invisibly.
#' @noRd
integrationRecordResultsPublication <- function(manager, model_reference, results, files, identity) {
    withWorkflowStateScope(manager, function() {
        context <- manager$.context
        snapshot <- context$snapshot
        export <- snapshot$exports[[INTEGRATION_MODEL_OPERATION]]
        if (is.null(export) || !identical(export$reference$entity_id, model_reference$entity_id)) return(invisible(FALSE))
        # Plain lists and vectors only: a reopened store refuses any object (a data frame) among the head's roots.
        listed <- function(paths) list(file = basename(paths),
            sha256 = vapply(paths, \(path) digest::digest(file = path, algo = "sha256"), character(1), USE.NAMES = FALSE))
        export$results <- list(label = INTEGRATION_RESULTS_LABEL, identity = identity, model_entity_id = model_reference$entity_id,
            parameters = results$parameters, factors = results$factors,
            saved_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), tables = listed(files$tables), figures = listed(files$figures),
            association = as.list(results$association),
            claim = list(group_corrected = results$claim$group_corrected, corrections = as.list(results$claim$corrections),
                uncorrected = isTRUE(results$claim$uncorrected)),
            outputs = integrationOutputsRecord(results), output_tables = integrationStoreOutputTables(context, results))
        snapshot$exports[[INTEGRATION_MODEL_OPERATION]] <- export
        publishWorkflowStateHead(context, snapshot)
        invisible(TRUE)
    })
}

#' What the further figures of a saved publication show (KIT-152; audits ENG-15, SCI-15 to SCI-17), as the head keeps it:
#' plain lists of plain values only.
#' @param results `integrationResults()`'s value.
#' @return `list(feature_values, top_features, value_scales = list(view, scale), heatmap_scaling, maps)`: `maps` per map
#'   drawn its settings and the software's warnings (`notes`).
#' @noRd
integrationOutputsRecord <- function(results) {
    p <- results$parameters
    maps <- lapply(Filter(function(map) is.null(map$refusal), results$maps), function(map)
        c(map$settings, list(notes = as.character(map$notes))))
    list(feature_values = p$feature_values, top_features = p$top_features,
        value_scales = list(view = as.character(results$value_scales$view), scale = as.character(results$value_scales$scale)),
        heatmap_scaling = p$heatmap_scaling, maps = maps)
}

#' Keep the further figures' tables in the integration store (KIT-152; ENG-15): the top features' values and each
#' map's coordinates, each a table of the store (`saveWorkflowOwnedTable`; the same table is the same entity).
#' @param context The owner's context (in its state scope).
#' @param results `integrationResults()`'s value.
#' @return Named list of table references (`feature_values`, `map_umap`, `map_tsne`).
#' @noRd
integrationStoreOutputTables <- function(context, results) {
    out <- list(feature_values = saveWorkflowOwnedTable(context, results$feature_values))
    for (method in names(results$maps)) {
        map <- results$maps[[method]]
        if (is.null(map$refusal)) out[[paste0("map_", method)]] <- saveWorkflowOwnedTable(context, map$coordinates)
    }
    out
}

#' The further figures' tables of the latest saved publication, read from the store (KIT-152; for the replay's
#' comparison, before the store is closed).
#' @param manager The integration owner.
#' @return Named list of data frames (`feature_values`, `map_umap`, `map_tsne`), or NULL when none were kept.
#' @noRd
integrationRecordedOutputs <- function(manager) {
    tables <- manager$.context$snapshot$exports[[INTEGRATION_MODEL_OPERATION]]$results$output_tables
    if (!length(tables)) return(NULL)
    withWorkflowStateScope(manager, function() lapply(tables, function(reference) materializeWorkflowTable(manager$.context,
        reference)))
}

#' The integration block of the model's analysis in the record (KIT-142), and each model setting's status: called by
#' `buildWorkflowAnalysisRecord` for an integration store only, so the omics' records are as before.
#' @param manager The integration owner.
#' @param analysis The model's analysis as `workflowRecordExportAnalysis` lists it.
#' @param export The accepted model (`snapshot$exports$integration_mofa_model`).
#' @return The analysis with `status` and `settled_by` on each parameter (the uploaded match's provenance `user_set`)
#'   and the block `integration`.
#' @noRd
integrationRecordAnalysis <- function(manager, analysis, export) {
    spec <- integrationMofaParameterSpec()
    parameters <- export$parameters
    analysis$parameters <- lapply(analysis$parameters, function(parameter) {
        name <- parameter$name
        # KIT-063 item A's match is the user's own table; it has no default to compare with (it was "unknown").
        if (name %in% c("sample_match", "sample_match_table")) parameter$provenance <- "user_set"
        c(parameter, list(status = integrationParameterStatus(name, parameters[[name]], spec),
            settled_by = spec[[name]]$settled_by %||% NA_character_))
    })
    context <- manager$.context
    value <- tryCatch(hydrateSavedWorkflowEntity(context$store, export$reference)$value, error = function(error) NULL)
    object <- tryCatch(withWorkflowStateScope(manager, function() readWorkflowExecutionComponents(manager, export$source,
        fields = "inputs")), error = function(error) NULL)
    if (!is.list(value) || is.null(object)) {
        analysis$integration <- list(schema = "multischolar.integration_record.v1", complete = FALSE)
        return(analysis)
    }
    views <- as.character(parameters$views)
    # KIT-153: the matrix the model used (`input_matrix`) decides the claim limit and each input's values.
    input_matrix <- parameters$input_matrix
    corrected <- integrationGroupCorrectedViews(object, views, input_matrix)
    available <- integrationCorrectedViewNames(object)
    gui <- integrationGuiTexts()
    inputs <- lapply(unname(object@inputs), function(entry) {
        manifest <- entry$manifest
        list(omic = entry$omic, views = as.character(entry$views), views_in_model = intersect(as.character(entry$views), views),
            manifest_path = entry$path, manifest_sha256 = entry$manifest_sha256,
            workflow = integrationGuiScalar(manifest$workflow), matrix_role = integrationGuiScalar(manifest$matrix_role),
            matrix_used = if (is.null(input_matrix)) NA_character_ else input_matrix,
            corrected_available = any(as.character(entry$views) %in% available),
            scale = integrationGuiScalar(manifest$scale), source_state = integrationGuiScalar(manifest$source_state$name),
            biological_sample_column = entry$biological_sample_column, correction = manifest$correction,
            correction_text = integrationCorrectionText(manifest, gui),
            group_corrected = any(corrected[intersect(as.character(entry$views), views)]))
    })
    match <- parameters$sample_match_table
    iterations <- as.integer(value$elbo$iterations)
    stopped_by <- value$elbo$stopped_by %||% NA_character_
    if (is.na(stopped_by)) stopped_by <- integrationStoppedBy(iterations, parameters$maxiter)
    results <- export$results
    if (!is.null(results)) {
        # KIT-152: the further figures' tables stay in the store; the record keeps what they show and their settings.
        results$output_tables <- NULL
        for (name in c("tables", "figures", "association")) {
            results[[name]] <- as.data.frame(results[[name]], stringsAsFactors = FALSE)
        }
        results_spec <- integrationResultsParameterSpec()
        results$parameters <- lapply(names(results$parameters), function(name) list(name = name,
            value = encodeWorkflowRecordValue(results$parameters[[name]]),
            status = integrationParameterStatus(name, results$parameters[[name]], results_spec, results$factors),
            settled_by = results_spec[[name]]$settled_by %||% NA_character_))
    }
    advisories <- lapply(integrationModelAdvisories(value, gui), function(advisory) advisory[c("kind", "sentence", "detail")])
    if (length(value$notes)) advisories <- c(advisories, list(list(kind = "few_samples",
        sentence = sprintf(integrationRecordTexts()[["note_few_samples"]], as.integer(value$shared_samples),
            INTEGRATION_FEW_SAMPLES), detail = paste(value$notes, collapse = " "))))
    claim_facts <- list(group_corrected = names(corrected)[corrected],
        corrections = vapply(object@inputs, function(entry) integrationCorrectionText(entry$manifest, gui), character(1)),
        uncorrected = identical(input_matrix, "uncorrected"))
    analysis$integration <- list(schema = "multischolar.integration_record.v1", complete = TRUE,
        inputs = inputs,
        matching = list(rule = if (is.null(match)) "saved_names" else "user_table",
            sample_match_sha256 = parameters$sample_match %||% NA_character_,
            runs_changed = if (is.null(match)) 0L else length(match$run), table = match,
            technical_replicates = parameters$technical_replicates, internal_standards = parameters$internal_standards,
            missing_samples = parameters$missing_samples, samples = value$samples, view_summary = value$view_summary,
            shared_samples = as.integer(value$shared_samples), confirmation = export$confirmation),
        model = list(entity_id = export$reference$entity_id,
            identity_sha256 = digest::digest(export$reference$identity, algo = "sha256"),
            model_file_sha256 = value$model_file$sha256 %||% NA_character_, model_file_bytes = value$model_file$size %||% NA_real_,
            precision = parameters$precision, iterations = iterations, maxiter = as.integer(parameters$maxiter),
            # KIT-151 (INT-15): stopped by the convergence rule or at the limit, as the product keeps it.
            stopped_by = stopped_by, converged = identical(stopped_by, "convergence_rule"),
            # KIT-151 (INT-14): the models trained, one per seed, and each factor's agreement across them.
            seeds = value$seeds$seeds, seed_models = value$seeds$models, stability = value$seeds$stability,
            factors_requested = as.integer(parameters$num_factors), factors_kept = ncol(value$factors) - 1L,
            elbo_last = value$elbo$last, variance_total = value$variance_total, versions = value$versions,
            # MOFA2's data, model and training options as trained (every one, the defaults included), encoded exactly.
            mofa2_options = encodeWorkflowRecordValue(value$options)),
        results = results,
        advisories = advisories,
        claim_limit = list(applies = any(corrected), group_corrected = claim_facts$group_corrected,
            uncorrected = claim_facts$uncorrected, text = integrationClaimLimitText(list(claim = claim_facts), gui)),
        # KIT-141: the enrichments of the factors (references stay in the store; the record keeps their facts and terms).
        enrichment = integrationRecordEnrichmentBlock(context, export))
    analysis
}

# ------------------------------------------------------------------------------------------------
# The record in words (KIT-142 item 3): the tables, the study parameters, the methods paragraph.
# ------------------------------------------------------------------------------------------------

#' The integration record as the tables people read (KIT-142): one source for the workbook's integration sheets and the
#' report. Every header is from `integrationRecordTexts()`.
#' @param record An analysis record of an integration store (built or read from JSON).
#' @return Named list of data frames (`INTEGRATION_RECORD_SHEETS`); empty list without an integration model.
#' @noRd
integrationRecordTables <- function(record, texts = integrationRecordTexts()) {
    analysis <- integrationRecordAnalysisOf(record)
    block <- analysis$integration
    if (is.null(block) || !isTRUE(integrationRecordScalar(block$complete))) return(list())
    s <- integrationRecordScalar
    word <- function(value) if (isTRUE(as.logical(s(value)))) texts[["yes"]] else texts[["no"]]
    text <- function(value) { value <- s(value); if (is.na(value)) texts[["not_recorded"]] else as.character(value) }
    framed <- function(names, rows) stats::setNames(as.data.frame(rows, stringsAsFactors = FALSE), names)
    inputs <- lapply(block$inputs, function(input) data.frame(omic = s(input$omic),
        views = paste(unlist(input$views_in_model), collapse = ", "),
        values = texts[[integrationRecordValuesKey(input)]],
        correction = text(input$correction_text), group_corrected = word(input$group_corrected),
        checksum = text(input$manifest_sha256), stringsAsFactors = FALSE))
    inputs <- do.call(rbind, inputs)
    names(inputs) <- unname(texts[c("col_omic", "col_datasets_in_model", "col_values", "col_correction",
        "col_group_corrected", "col_checksum")])
    m <- block$matching
    confirmation <- m$confirmation
    matching <- framed(unname(texts[c("col_field", "col_value")]), list(
        c(texts[["match_rule"]], texts[["match_table_checksum"]], texts[["match_samples"]], texts[["match_shared"]],
            texts[["match_confirmed_at"]], texts[["match_design_file"]], texts[["match_design_checksum"]], texts[["match_file"]]),
        c(if (identical(s(m$rule), "user_table")) sprintf(texts[["match_rule_user"]], as.integer(s(m$runs_changed))) else
            texts[["match_rule_saved"]], text(m$sample_match_sha256), as.character(nrow(integrationRecordFrame(m$samples))),
            text(m$shared_samples), text(confirmation$confirmed_at), integrationRecordBasename(confirmation$design_file, texts),
            text(confirmation$design_sha256), integrationRecordBasename(confirmation$sample_match_file, texts))))
    samples <- integrationRecordFrame(m$samples)
    status_word <- function(status) integrationTextOr(texts, paste0("status_", s(status)), text(status))
    model_labels <- workflowParameterLabels[[INTEGRATION_MODEL_OPERATION]]
    setting <- function(applies, parameter, label) data.frame(applies = applies, setting = label,
        name = as.character(s(parameter$name)), value = workflowRecordValueText(parameter$value),
        status = status_word(parameter$status), settled_by = text(parameter$settled_by), stringsAsFactors = FALSE)
    settings <- lapply(analysis$parameters, function(parameter) {
        name <- as.character(s(parameter$name))
        if (identical(name, "sample_match_table")) return(NULL)
        setting(texts[["applies_model"]], parameter, model_labels[[name]] %||% name)
    })
    results <- block$results
    for (parameter in results$parameters %||% list()) {
        settings[[length(settings) + 1L]] <- setting(texts[["applies_results"]], parameter, as.character(s(parameter$name)))
    }
    settings <- do.call(rbind, Filter(Negate(is.null), settings))
    names(settings) <- unname(texts[c("col_applies_to", "col_setting", "col_name", "col_value", "col_status", "col_settled_by")])
    model <- block$model
    versions <- model$versions
    model_table <- framed(unname(texts[c("col_field", "col_value")]), list(
        unname(texts[c("model_mofa2", "model_mofapy2", "model_python", "model_r", "model_reticulate", "model_precision",
            "model_iterations", "model_settled", "model_factors", "model_seeds", "model_checksum", "model_bytes",
            "model_entity", "model_identity")]),
        c(text(versions$MOFA2), text(versions$mofapy2), text(versions$python), text(versions$R), text(versions$reticulate),
            text(model$precision), sprintf("%s (%s)", text(model$iterations), text(model$maxiter)), word(model$converged),
            sprintf(texts[["model_factors_value"]], as.integer(s(model$factors_kept)), as.integer(s(model$factors_requested))),
            if (length(unlist(model$seeds))) sprintf(texts[["model_seeds_value"]], length(unlist(model$seeds)),
                paste(unlist(model$seeds), collapse = ", ")) else texts[["not_recorded"]],
            text(model$model_file_sha256), text(model$model_file_bytes), text(model$entity_id), text(model$identity_sha256))))
    variance <- integrationRecordFrame(model$variance_total)
    variance <- if (nrow(variance)) data.frame(variance$view, sprintf("%.2f", as.numeric(variance$r2_percent)),
        stringsAsFactors = FALSE) else data.frame(character(), character())
    names(variance) <- unname(texts[c("col_view", "col_variance")])
    listed <- function(x) { x <- integrationRecordFrame(x); if (nrow(x)) data.frame(x$file, x$sha256, stringsAsFactors = FALSE) else NULL }
    files <- rbind(listed(results$tables), listed(results$figures)) %||% data.frame(character(), character())
    names(files) <- unname(texts[c("col_file", "col_checksum")])
    notes <- c(vapply(block$advisories %||% list(), \(advisory) as.character(s(advisory$sentence)), character(1)),
        as.character(unlist(block$claim_limit$text)))
    notes <- framed(texts[["col_note"]], list(notes))
    enrichment <- integrationRecordEnrichmentTables(block, texts)
    out <- list(integration_inputs = inputs, integration_matching = matching, integration_samples = samples,
        integration_settings = settings, integration_model = model_table, integration_variance = variance,
        integration_results = files, integration_outputs = integrationRecordOutputsTable(block, texts),
        integration_notes = notes, integration_enrichment = enrichment$summary,
        integration_enrichment_terms = enrichment$terms)
    out[INTEGRATION_RECORD_SHEETS]
}

#' Which values an input gave the model, as the key of its words (KIT-153): `values_corrected` (the RUV-corrected
#' matrix), `values_uncorrected_available` (the normalised matrix, the omic's corrected one not used) or
#' `values_normalised`. A record made before KIT-153 (no `matrix_used`) is read by the manifest's role, as before.
#' @param input One input of the record's integration block.
#' @noRd
integrationRecordValuesKey <- function(input) {
    s <- integrationRecordScalar
    used <- s(input$matrix_used)
    if (is.na(used)) used <- if (identical(s(input$matrix_role), "corrected")) "ruv_corrected" else "uncorrected"
    if (identical(used, "ruv_corrected")) return("values_corrected")
    if (isTRUE(as.logical(s(input$corrected_available)))) "values_uncorrected_available" else "values_normalised"
}

#' The scales of the values given to the model, in words, from the record (KIT-152): "data set: scale" per data set.
#' @param outputs The record's `results$outputs`.
#' @noRd
integrationRecordScaleText <- function(outputs, texts = integrationRecordTexts()) {
    views <- as.character(unlist(outputs$value_scales$view))
    scales <- as.character(unlist(outputs$value_scales$scale))
    words <- vapply(scales, function(scale) integrationTextOr(texts, paste0("scale_", scale), scale), character(1))
    paste(sprintf(texts[["scale_item"]], views, words), collapse = "; ")
}

#' The further figures of the saved results as the record's table (KIT-152; ENG-15): one row per figure, what it
#' shows and its settings; the maps only when drawn, then the caution about them. No row without saved results, or for
#' results saved before KIT-152.
#' @param block The record's integration block.
#' @return `data.frame` with the headers `col_output`, `col_shows`, `col_settings`.
#' @noRd
integrationRecordOutputsTable <- function(block, texts = integrationRecordTexts()) {
    s <- integrationRecordScalar
    outputs <- block$results$outputs
    rows <- list()
    if (!is.null(outputs)) {
        heatmap <- as.character(s(outputs$heatmap_scaling))
        rows <- list(c(texts[["out_pairs"]], texts[["out_pairs_shows"]], ""),
            c(texts[["out_top"]], sprintf(texts[["out_top_shows"]], as.integer(s(outputs$top_features))), ""),
            c(texts[["out_scatter"]], sprintf(texts[["out_scatter_shows"]], integrationRecordScaleText(outputs, texts)), ""),
            c(texts[["out_heatmap"]], integrationTextOr(texts, paste0("out_heatmap_shows_", heatmap), heatmap), ""))
        maps <- outputs$maps %||% list()
        for (method in intersect(c("umap", "tsne"), names(maps))) {
            m <- maps[[method]]
            settings <- if (method == "umap") sprintf(texts[["out_umap_settings"]], s(m$version), as.integer(s(m$n_neighbors)),
                format(s(m$min_dist)), as.integer(s(m$n_epochs)), s(m$init), as.integer(s(m$seed))) else
                sprintf(texts[["out_tsne_settings"]], s(m$version), format(s(m$perplexity)), as.integer(s(m$max_iter)),
                    as.integer(s(m$seed)))
            rows[[length(rows) + 1L]] <- c(texts[[paste0("out_map_", method)]],
                sprintf(texts[["out_map_shows"]], as.integer(s(m$factors)), as.integer(s(m$samples))), settings)
        }
        if (length(maps)) {
            samples <- as.integer(s(maps[[1L]]$samples))
            corrected <- as.character(unlist(block$claim_limit$group_corrected))
            caution <- integrationMapCaution(list(claim = list(group_corrected = corrected)), samples, names(maps))
            rows[[length(rows) + 1L]] <- c(texts[["out_caution"]], paste(caution, collapse = " "), "")
        }
    }
    table <- if (length(rows)) as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE) else
        data.frame(character(), character(), character(), stringsAsFactors = FALSE)
    names(table) <- unname(texts[c("col_output", "col_shows", "col_settings")])
    rownames(table) <- NULL
    table
}

#' The methods paragraph of an integration record (KIT-142): a draft for the user, stating the facts the record holds
#' (the data sets and their corrections, the software and its precision, how samples were matched, what happens before
#' the model, the model's settings and training, the factor-group statistic, the claim limit, and which settings were
#' provisional defaults or changed). Every value is read from the record.
#' @param record An analysis record of an integration store (built or read from JSON).
#' @return One string; NA without an integration model.
#' @noRd
integrationMethodsText <- function(record, texts = integrationRecordTexts()) {
    analysis <- integrationRecordAnalysisOf(record)
    block <- analysis$integration
    if (is.null(block) || !isTRUE(integrationRecordScalar(block$complete))) return(NA_character_)
    s <- integrationRecordScalar
    p <- function(name) integrationRecordParameter(analysis$parameters, name)
    used <- Filter(\(input) length(unlist(input$views_in_model)) > 0L, block$inputs)
    # KIT-153: the values each input gave the model (`integrationRecordValuesKey`).
    items <- vapply(used, function(input) sprintf(texts[["methods_input_item"]], s(input$omic),
        paste(unlist(input$views_in_model), collapse = ", "),
        switch(integrationRecordValuesKey(input),
            values_corrected = sprintf(texts[["methods_input_corrected"]], s(input$correction_text)),
            values_uncorrected_available = sprintf(texts[["methods_input_uncorrected_available"]], s(input$correction_text)),
            texts[["methods_input_normalised"]])), character(1))
    views <- unlist(lapply(used, \(input) unlist(input$views_in_model)))
    model <- block$model
    versions <- model$versions
    version <- function(x) { x <- s(x); if (is.na(x)) texts[["not_recorded"]] else sub(" .*$", "", as.character(x)) }
    m <- block$matching
    sentences <- c(
        sprintf(texts[["methods_inputs"]], length(views), length(used), paste(items, collapse = "; ")),
        sprintf(texts[["methods_software"]], version(versions$MOFA2), version(versions$mofapy2), version(versions$python),
            version(versions$R), as.character(s(model$precision))),
        if (identical(s(m$rule), "user_table")) sprintf(texts[["methods_matching_user"]], as.integer(s(m$runs_changed)),
            s(m$sample_match_sha256)) else texts[["methods_matching_saved"]])
    if (identical(p("technical_replicates"), "average")) sentences <- c(sentences, texts[["methods_replicates_average"]])
    summary <- integrationRecordFrame(m$view_summary)
    standards <- if (nrow(summary)) sum(as.integer(summary$internal_standards_dropped), na.rm = TRUE) else 0L
    # KIT-151 (INT-12): with a feature selection, the observed share and the selection rule.
    selection <- as.character(p("feature_selection"))
    selected <- selection %in% c("top_variance", "variance_above_trend")
    samples_n <- nrow(integrationRecordFrame(m$samples))
    fraction <- as.numeric(p("min_observed_fraction"))
    sparse <- if (selected && !is.na(fraction)) sprintf(texts[["methods_sparse_selected"]],
        integrationMinObserved(list(feature_selection = selection, min_observed_fraction = fraction), samples_n), samples_n,
        format(fraction)) else texts[["methods_sparse"]]
    sentences <- c(sentences,
        if (identical(p("internal_standards"), "drop")) sprintf(texts[["methods_standards_drop"]], standards) else
            texts[["methods_standards_keep"]],
        integrationTextOr(texts, paste0("methods_missing_", p("missing_samples")), NULL),
        sparse,
        paste(if (selected) sprintf(texts[[paste0("methods_selection_", selection)]], as.integer(p("features_per_view")),
            as.integer(p("features_per_view"))) else texts[["methods_selection_none"]],
            integrationTextOr(texts, paste0("methods_scaling_", p("scaling")), as.character(p("scaling"))), sep = ", "),
        sprintf(texts[["methods_model"]], as.integer(p("num_factors")), as.character(p("convergence_mode")),
            as.integer(p("maxiter")), as.integer(p("seed")), if (isTRUE(p("scale_views"))) texts[["yes"]] else texts[["no"]],
            if (isTRUE(p("center_groups"))) texts[["yes"]] else texts[["no"]]),
        if (isTRUE(as.logical(s(model$converged)))) sprintf(texts[["methods_trained_settled"]], as.integer(s(model$iterations)),
            as.integer(s(model$factors_kept)), nrow(integrationRecordFrame(m$samples)), as.integer(s(m$shared_samples))) else
            sprintf(texts[["methods_trained_not_settled"]], as.integer(s(model$maxiter)), as.integer(s(model$factors_kept)),
                nrow(integrationRecordFrame(m$samples)), as.integer(s(m$shared_samples))))
    # KIT-151 (INT-17): whether training ran with the maths library pinned (nothing said for a model without the fact).
    pin <- versions$maths_pin
    if (!is.null(pin) && !is.na(s(pin$applied))) {
        sentences <- c(sentences, if (isTRUE(as.logical(s(pin$applied)))) sprintf(texts[["methods_pin_applied"]],
            paste(sprintf("%s=%s", names(pin$variables), vapply(pin$variables, \(x) as.character(s(x)), character(1))),
                collapse = ", ")) else sprintf(texts[["methods_pin_not_applied"]], as.character(s(pin$reason))))
    }
    # KIT-151 (INT-14, audit SCI-20): how many models were trained, from which seeds, and each factor's agreement.
    seeds <- as.integer(unlist(model$seeds))
    stability <- integrationRecordFrame(model$stability)
    if (length(seeds) > 1L && nrow(stability)) {
        sentences <- c(sentences, sprintf(texts[["methods_seeds"]], length(seeds), seeds[[1L]], seeds[[length(seeds)]],
            seeds[[1L]], paste(sprintf("%s %.3f", stability$factor, as.numeric(stability$agreement_lowest)), collapse = ", ")))
    } else if (length(seeds) == 1L) {
        sentences <- c(sentences, sprintf(texts[["methods_seeds_single"]], seeds[[1L]]))
    }
    results <- block$results
    if (!is.null(results)) {
        adjust <- integrationRecordParameter(results$parameters, "p_adjust_method")
        test <- integrationRecordParameter(results$parameters, "association_test")
        sentences <- c(sentences, sprintf(texts[[if (identical(test, "kruskal_wallis")) "methods_association_kruskal_wallis" else
            "methods_association"]],
            integrationTextOr(integrationGuiTexts(), paste0("padjust_", adjust), as.character(adjust))))
        sentences <- c(sentences, integrationMethodsOutputs(results$outputs, block$claim_limit, texts))
    }
    # KIT-141: the enrichment of the factors.
    sentences <- c(sentences, integrationMethodsEnrichment(block, texts))
    claim <- block$claim_limit
    sentences <- c(sentences, if (isTRUE(as.logical(s(claim$applies))))
        sprintf(texts[["methods_claim_limit"]], paste(unlist(claim$group_corrected), collapse = ", ")) else texts[["methods_claim_none"]])
    # KIT-153: the uncorrected matrices' caution.
    if (isTRUE(as.logical(s(claim$uncorrected)))) sentences <- c(sentences, texts[["methods_claim_uncorrected"]])
    # The settings in the Model tab's words: its label, and the chosen value's words where the tab has them.
    gui <- integrationGuiTexts()
    by_status <- function(status) {
        hits <- Filter(\(parameter) identical(s(parameter$status), status), analysis$parameters)
        names <- vapply(hits, \(parameter) as.character(s(parameter$name)), character(1))
        names <- names[!names %in% c("sample_match", "sample_match_table") & !vapply(names, \(name) all(is.na(p(name))), logical(1))]
        vapply(names, function(name) {
            value <- p(name)
            words <- if (is.character(value) && length(value) == 1L) integrationParameterText(gui, "choice", name, value) else
                if (is.logical(value)) (if (isTRUE(value)) texts[["yes"]] else texts[["no"]]) else paste(format(value), collapse = ", ")
            sprintf("%s (%s)", integrationParameterText(gui, "param", name), words)
        }, character(1), USE.NAMES = FALSE)
    }
    provisional <- by_status("provisional_default")
    if (length(provisional)) sentences <- c(sentences, sprintf(texts[["methods_provisional"]], paste(provisional, collapse = "; ")))
    changed <- by_status("user_set")
    if (length(changed)) sentences <- c(sentences, sprintf(texts[["methods_user_set"]], paste(changed, collapse = "; ")))
    paste(sentences[nzchar(sentences)], collapse = " ")
}

#' The methods paragraph's sentences on the further figures (KIT-152; audits SCI-15 to SCI-17), from the record: how the
#' top features are ranked and what their lists claim; which values the scatter and heatmap show; each map drawn with
#' its software and settings, then the caution (and the group-correction sentence where the claim limit applies). A map
#' not drawn is not mentioned; results saved before KIT-152 give none.
#' @param outputs The record's `results$outputs`.
#' @param claim The record's `claim_limit`.
#' @return A character vector.
#' @noRd
integrationMethodsOutputs <- function(outputs, claim, texts = integrationRecordTexts()) {
    if (is.null(outputs)) return(character())
    s <- integrationRecordScalar
    heatmap <- as.character(s(outputs$heatmap_scaling))
    sentences <- c(sprintf(texts[["methods_top_features"]], as.integer(s(outputs$top_features))),
        sprintf(texts[["methods_feature_values"]], integrationRecordScaleText(outputs, texts)),
        integrationTextOr(texts, paste0("methods_heatmap_", heatmap), ""))
    maps <- outputs$maps %||% list()
    for (method in intersect(c("umap", "tsne"), names(maps))) {
        m <- maps[[method]]
        sentences <- c(sentences, if (method == "umap") sprintf(texts[["methods_umap"]], as.integer(s(m$factors)), s(m$version),
            as.integer(s(m$n_neighbors)), format(s(m$min_dist)), as.integer(s(m$n_epochs)), s(m$init), as.integer(s(m$seed))) else
            sprintf(texts[["methods_tsne"]], as.integer(s(m$factors)), s(m$version), format(s(m$perplexity)),
                as.integer(s(m$max_iter)), as.integer(s(m$seed))))
    }
    if (length(maps)) {
        sentences <- c(sentences, sprintf(texts[["methods_maps_caution"]], as.integer(s(maps[[1L]]$samples))),
            unname(texts[paste0("methods_maps_", intersect(c("umap", "tsne"), names(maps)))]))
        corrected <- as.character(unlist(claim$group_corrected))
        if (isTRUE(as.logical(s(claim$applies))) && length(corrected)) {
            sentences <- c(sentences, sprintf(texts[["methods_maps_group_corrected"]], paste(corrected, collapse = ", ")))
        }
    }
    sentences
}

#' The integration's study parameters (KIT-142), in the form of the omics' `study_parameters.txt`: one section per
#' table of `integrationRecordTables()` with "* item: value" lines, then the methods paragraph.
#' @param record An analysis record of an integration store.
#' @return A character vector (the file's lines).
#' @noRd
integrationStudyParameterLines <- function(record, texts = integrationRecordTexts()) {
    tables <- integrationRecordTables(record, texts)
    if (!length(tables)) return(character())
    section <- function(heading, lines) c(heading, strrep("-", nchar(heading)), lines, "")
    rows <- function(table) if (!nrow(table)) sprintf("* %s", texts[["none_recorded"]]) else
        vapply(seq_len(nrow(table)), \(i) paste0("* ", paste(sprintf("%s: %s", names(table), unlist(table[i, ], use.names = FALSE)),
            collapse = "; ")), character(1))
    pairs <- function(table) sprintf("* %s: %s", table[[1L]], table[[2L]])
    # KIT-141: an enrichment sheet without an enrichment holds one row saying so.
    enriched <- function(table) if (identical(as.character(table[[1L]]), unname(texts[["none_recorded"]])))
        sprintf("* %s", texts[["none_recorded"]]) else rows(table)
    c(section(texts[["report_title"]], texts[["report_intro"]]),
        section(texts[["section_inputs"]], rows(tables$integration_inputs)),
        section(texts[["section_matching"]], pairs(tables$integration_matching)),
        section(texts[["section_settings"]], rows(tables$integration_settings)),
        section(texts[["section_model"]], pairs(tables$integration_model)),
        section(texts[["section_variance"]], pairs(tables$integration_variance)),
        section(texts[["section_results"]], rows(tables$integration_results)),
        # KIT-152: the further figures.
        section(texts[["section_outputs"]], rows(tables$integration_outputs)),
        # KIT-141: the enrichments and their terms.
        section(texts[["section_enrichment"]], enriched(tables$integration_enrichment)),
        section(texts[["section_enrichment_terms"]], enriched(tables$integration_enrichment_terms)),
        section(texts[["section_notes"]], if (nrow(tables$integration_notes)) sprintf("* %s", tables$integration_notes[[1L]]) else
            sprintf("* %s", texts[["none_recorded"]])),
        section(texts[["section_methods"]], c(texts[["methods_draft_note"]], "", integrationMethodsText(record, texts))))
}

#' Render the integration report (KIT-142) from a record file with the package's template
#' (`inst/reports/integration/integration_report.rmd`), in a temporary folder of its own.
#' @param record_path The record (`analysis_record.json`).
#' @param output_dir Where the report is written.
#' @param format `docx` (the omics' reports) or `html`.
#' @return The report's path.
#' @noRd
renderIntegrationReport <- function(record_path, output_dir, format = c("docx", "html")) {
    format <- match.arg(format)
    record_path <- normalizePath(record_path, winslash = "/", mustWork = TRUE)
    output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
    template <- system.file("reports", "integration", "integration_report.rmd", package = "MultiScholaR")
    if (!nzchar(template)) cli::cli_abort("the integration report template is not installed", class = "integration_report_missing")
    work <- tempfile("ir-")
    dir.create(work)
    on.exit(unlink(work, recursive = TRUE), add = TRUE)
    file.copy(template, work)
    rmarkdown::render(file.path(work, basename(template)), output_format = if (format == "docx") "word_document" else "html_document",
        output_file = paste0("integration_report.", format), output_dir = output_dir, intermediates_dir = work,
        params = list(record_path = record_path), envir = new.env(parent = globalenv()), quiet = TRUE)
}

#' Save the integration's record and report into a folder (KIT-142): `analysis_record.json`, the record's workbook
#' (`analysis_record.xlsx`, with the integration sheets), `study_parameters.txt` and the rendered report; KIT-148: with
#' them the model file MOFA2 wrote (`model.hdf5`, the trained model without the training data; its SHA-256 is the
#' record's) and the uploaded sample match (`sample_match.tsv`: the runs it gives another biological sample, `omic`,
#' `run`, `biological_sample`; the header alone when none was uploaded). Everything is
#' written in a short temporary folder and copied into `directory` with every file checked
#' (`withWorkflowShortFigureFolder`, KIT-052): a path too long for Windows, or a file the renderer did not write, stops
#' the save by name instead of losing the file silently.
#' @param manager The integration owner with an accepted model.
#' @param directory The folder (the project's study report folder).
#' @param format The report's format (`docx` or `html`).
#' @param renderFn `function(record_path, output_dir, format)` writing the report (`renderIntegrationReport`).
#' @param now The record's time.
#' @return `list(directory, files, methods)`: the folder, the files written there, the methods paragraph.
#' @noRd
integrationWriteRecordBundle <- function(manager, directory, format = "docx", renderFn = renderIntegrationReport, now = Sys.time()) {
    texts <- integrationRecordTexts()
    record <- buildWorkflowAnalysisRecord(manager, now)
    if (is.null(integrationRecordAnalysisOf(record))) {
        integrationGuiRefuse(texts[["record_no_model"]], "integration_record_no_model")
    }
    export <- manager$.context$snapshot$exports[[INTEGRATION_MODEL_OPERATION]]
    value <- withWorkflowStateScope(manager, function() hydrateSavedWorkflowEntity(manager$.context$store, export$reference)$value)
    match <- export$parameters$sample_match_table
    match <- if (is.null(match)) data.frame(omic = character(), run = character(), biological_sample = character()) else
        data.frame(omic = match$omic, run = match$run, biological_sample = match$biological_sample, stringsAsFactors = FALSE)
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    out <- withWorkflowShortFigureFolder(directory, function(scratch) {
        writeBin(integrationModelFileBytes(manager, value), file.path(scratch, "model.hdf5"))
        utils::write.table(match, file.path(scratch, "sample_match.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
        json <- file.path(scratch, "analysis_record.json")
        jsonlite::write_json(record, json, auto_unbox = TRUE, digits = I(17), null = "null", na = "null", pretty = TRUE)
        read <- jsonlite::read_json(json, simplifyVector = FALSE)
        writeWorkflowAnalysisRecordWorkbook(read, file.path(scratch, "analysis_record.xlsx"))
        writeLines(integrationStudyParameterLines(read, texts), file.path(scratch, "study_parameters.txt"), useBytes = TRUE)
        renderFn(json, scratch, format)
        expected <- c("analysis_record.json", "analysis_record.xlsx", "study_parameters.txt", paste0("integration_report.", format),
            "model.hdf5", "sample_match.tsv")
        checkWorkflowFigureFiles(file.path(scratch, expected), "the integration record or report was not written")
        list(files = expected, methods = integrationMethodsText(read, texts))
    }, message = "the integration record and report could not be copied into the project")
    list(directory = directory, files = file.path(directory, out$files), methods = out$methods)
}

# ------------------------------------------------------------------------------------------------
# Replay (KIT-142 item 2): the record alone, with its inputs checked by checksum, gives the model again.
# ------------------------------------------------------------------------------------------------

#' Refuse a replay in plain words, classed `integration_replay_refused` and `integration_refused` with its own class.
#' @noRd
integrationReplayRefuse <- function(sentence, class) {
    integrationGuiRefuse(sentence, c(class, "integration_replay_refused"))
}

#' Compare the recorded software versions with this computer's (KIT-142): MOFA2, mofapy2, Python and R.
#' @param recorded The record's versions.
#' @param current This computer's (`integrationMofaVersions()`).
#' @return `data.frame(software, recorded, current, same)`.
#' @noRd
integrationReplayVersions <- function(recorded, current) {
    names <- c("MOFA2", "mofapy2", "python", "R")
    one <- function(x) { x <- integrationRecordScalar(x); if (is.na(x)) NA_character_ else sub(" .*$", "", as.character(x)) }
    out <- data.frame(software = names, recorded = vapply(names, \(n) one(recorded[[n]]), character(1)),
        current = vapply(names, \(n) one(current[[n]]), character(1)), stringsAsFactors = FALSE, row.names = NULL)
    out$same <- !is.na(out$recorded) & !is.na(out$current) & out$recorded == out$current
    out
}

#' Compare a replayed model with the recorded one numerically (KIT-142; INT-6 as corrected by INT-17, KIT-151: the
#' same inputs, settings and seed give the identical model on the same computer and software; on another processor
#' they usually agree to about 1e-12, the maths library pinned or not, but a factor near the detection limit can
#' differ, so a difference beyond `tolerance` on another computer is not by itself an error). Per table (factors, weights, variance explained per factor and in all, the sample table, the data set
#' summary): the largest absolute difference and the largest relative difference (largest absolute over the largest
#' absolute value of the recorded table). Tables of another shape, other labels or missing values in other places
#' differ outright (Inf).
#' @param recorded,replayed The model products as held outside their stores (`integrationModelComparable()`: the weights
#'   as one table; KIT-148).
#' @param tolerance The relative tolerance.
#' @return `list(tables, largest_relative, equal, identical)`.
#' @noRd
integrationCompareModels <- function(recorded, replayed, tolerance = 1e-5) {
    tables <- do.call(rbind, lapply(c("factors", "weights", "variance_explained", "variance_total", "samples", "view_summary"),
        \(name) integrationCompareTable(name, recorded[[name]], replayed[[name]], tolerance)))
    list(tables = tables, largest_relative = max(tables$largest_relative), equal = all(tables$within_tolerance),
        identical = all(tables$identical))
}

#' Compare one recorded table with its repeat (KIT-142's rule, KIT-152 shared with the further figures): the largest
#' absolute and relative difference of the numeric columns; another shape, other labels or missing values in other
#' places differ outright (Inf).
#' @param name The table's name.
#' @param a,b The recorded and the repeated table.
#' @param tolerance The relative tolerance.
#' @return One row: `table, rows, largest_absolute, largest_relative, within_tolerance, identical`.
#' @noRd
integrationCompareTable <- function(name, a, b, tolerance = 1e-5) {
    a <- as.data.frame(a)
    b <- as.data.frame(b)
    same_shape <- identical(dim(a), dim(b)) && identical(names(a), names(b))
    numeric <- if (same_shape) names(a)[vapply(a, is.numeric, logical(1))] else character()
    labels_same <- same_shape && identical(a[setdiff(names(a), numeric)], b[setdiff(names(b), numeric)])
    x <- as.matrix(a[numeric])
    y <- as.matrix(b[numeric])
    na_same <- labels_same && identical(is.na(x), is.na(y))
    absolute <- if (na_same && length(x)) max(c(0, abs(x - y)), na.rm = TRUE) else if (na_same) 0 else Inf
    scale <- if (na_same && length(x)) max(c(0, abs(x)), na.rm = TRUE) else 0
    relative <- if (!is.finite(absolute)) Inf else if (absolute == 0) 0 else absolute / max(scale, .Machine$double.xmin)
    data.frame(table = name, rows = nrow(b), largest_absolute = absolute, largest_relative = relative,
        within_tolerance = relative <= tolerance, identical = identical(a, b), stringsAsFactors = FALSE)
}

#' Compare the further figures' recorded tables with those of the repeated analysis (KIT-152; ENG-15, ENG-22): the top
#' features' values (`top_feature_values`; the values given to the model, read from the store, so equal wherever the
#' model is) and each map's coordinates (`map_umap`, `map_tsne`). The maps' starts depend on the computer's numerical
#' libraries, so their coordinates are expected identical only on the same computer (`same_computer_only`).
#' @param recorded `integrationRecordedOutputs()`'s value.
#' @param results The repeated analysis's `integrationResults()`.
#' @param tolerance The relative tolerance.
#' @return `list(tables, equal, identical)`; `tables` with `same_computer_only`.
#' @noRd
integrationCompareOutputs <- function(recorded, results, tolerance = 1e-5) {
    names <- names(recorded)
    replayed <- lapply(stats::setNames(names, names), function(name) {
        if (identical(name, "feature_values")) return(results$feature_values)
        results$maps[[sub("^map_", "", name)]]$coordinates
    })
    # The columns compared as they are; a table read back from the store may carry other attributes.
    plain <- function(x) if (is.null(x)) data.frame() else data.frame(as.list(x), check.names = FALSE, stringsAsFactors = FALSE)
    tables <- do.call(rbind, lapply(names, function(name) integrationCompareTable(if (name == "feature_values")
        "top_feature_values" else name, plain(recorded[[name]]), plain(replayed[[name]]), tolerance)))
    tables$same_computer_only <- startsWith(names, "map_")
    list(tables = tables, equal = all(tables$within_tolerance), identical = all(tables$identical))
}

#' Repeat an integration analysis from its record (KIT-142). Before anything is trained, in plain words: every
#' recorded input must be found (at its recorded place, or where `manifests` says) with its recorded checksum
#' (`integration_replay_input_missing`, `integration_replay_input_changed`); MOFA2 and its Python must be available
#' (`integration_mofa_unavailable`). Software versions other than the recorded ones are a warning
#' (`integration_replay_versions_differ`), and the comparison is recorded. The inputs are imported into a new store under
#' `scratch`, the model trained with the recorded settings (the uploaded sample match included), the saved Results
#' settings applied, and the model compared with `recorded` when given.
#' @param record The analysis record (built, or read from JSON).
#' @param scratch A folder for the replay's store.
#' @param recorded The recorded model product as held outside its store (`integrationModelComparable()`), or NULL (no
#'   comparison).
#' @param manifests Named by omic: the input manifests' current paths (absent ones at their recorded place).
#' @param probe The MOFA2 capability probe.
#' @param trainer NULL for MOFA2, or a stand-in (tests).
#' @param current_versions This computer's versions; NULL takes `integrationMofaVersions()` when MOFA2 trains.
#' @param tolerance The relative tolerance of the comparison.
#' @param recorded_outputs The further figures' recorded tables (`integrationRecordedOutputs()`), or NULL (KIT-152).
#' @param recorded_enrichment KIT-141: the recorded enrichments read from the store (`integrationRecordedEnrichment()`),
#'   or NULL; each is repeated on the replayed model and compared (`integrationReplayEnrichment`).
#' @param enrichment_mode `recorded` (the recorded STRING answer and KEGG/Reactome tables; no network) or `query_again`.
#' @param enrichment_api_key The STRING key, for `query_again` only; never kept.
#' @param enrichment_query,enrichment_download The STRING client and the downloader (tests pass stand-ins).
#' @return `list(replayed, results, comparison, versions, parameters, outputs_comparison)`; `replayed` is the repeated
#'   model as held outside its store (`integrationModelComparable()`), since the replay's store is closed on return;
#'   `outputs_comparison` is `integrationCompareOutputs()`'s value when recorded outputs and saved Results settings are
#'   given; KIT-141: `enrichment_comparison` is `integrationReplayEnrichment()`'s value when recorded enrichments are given.
#' @noRd
replayIntegrationRecord <- function(record, scratch, recorded = NULL, manifests = NULL, probe = integrationMofaCapability,
    trainer = NULL, current_versions = NULL, tolerance = 1e-5, recorded_outputs = NULL, recorded_enrichment = NULL,
    enrichment_mode = "recorded", enrichment_api_key = NULL, enrichment_query = integrationStringRankedEnrichment,
    enrichment_download = integrationEnrichmentDownload) {
    texts <- integrationRecordTexts()
    analysis <- integrationRecordAnalysisOf(record)
    block <- analysis$integration
    if (is.null(analysis) || is.null(block) || !isTRUE(integrationRecordScalar(block$complete))) {
        integrationReplayRefuse(texts[["replay_no_model"]], "integration_replay_no_model")
    }
    s <- integrationRecordScalar
    paths <- character()
    for (input in block$inputs) {
        omic <- as.character(s(input$omic))
        path <- if (!is.null(manifests) && omic %in% names(manifests)) manifests[[omic]] else as.character(s(input$manifest_path))
        if (is.na(path) || !file.exists(path)) {
            integrationReplayRefuse(sprintf(texts[["replay_input_missing"]], omic, path), "integration_replay_input_missing")
        }
        if (!identical(digest::digest(file = path, algo = "sha256"), as.character(s(input$manifest_sha256)))) {
            integrationReplayRefuse(sprintf(texts[["replay_input_changed"]], omic), "integration_replay_input_changed")
        }
        paths[[omic]] <- path
    }
    if (is.null(trainer)) {
        capability <- probe()
        if (!isTRUE(capability$available)) {
            stop(structure(class = c("integration_mofa_unavailable", "integration_replay_refused", "integration_gui_refusal",
                "integration_refused", "error", "condition"), list(message = texts[["replay_mofa_unavailable"]], call = NULL,
                reason = capability$reason)))
        }
        current_versions <- current_versions %||% integrationMofaVersions()
    }
    versions <- NULL
    if (!is.null(current_versions)) {
        versions <- integrationReplayVersions(block$model$versions, current_versions)
        if (!all(versions$same)) {
            differ <- versions[!versions$same, , drop = FALSE]
            cli::cli_warn(sprintf(texts[["replay_versions_differ"]], paste(sprintf(texts[["replay_version_item"]], differ$software,
                ifelse(is.na(differ$recorded), texts[["not_recorded"]], differ$recorded),
                ifelse(is.na(differ$current), texts[["not_recorded"]], differ$current)), collapse = "; ")),
                class = "integration_replay_versions_differ")
        }
    }
    dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
    manager <- newIntegrationStateOwner(tempfile("replay-", tmpdir = scratch))
    on.exit(closeWorkflowStateOwner(manager), add = TRUE)
    importIntegrationInputs(manager, dirname(unname(paths)), description = "Integration inputs imported (replay)")
    entry <- workflowOperationEntry(INTEGRATION_MODEL_OPERATION)
    parameters <- workflowReplayParameters(analysis$parameters, c(entry$parameters, "sample_match", "sample_match_table"))
    run <- runIntegrationModel(manager, parameters, probe = probe, trainer = trainer)
    results <- NULL
    if (!is.null(block$results)) {
        settings <- list()
        for (parameter in block$results$parameters) settings[[as.character(s(parameter$name))]] <- decodeWorkflowRecordValue(parameter$value)
        results <- integrationResults(run$value, integrationImportedObject(manager), settings, manager = manager)
    }
    # KIT-141: the recorded enrichments, repeated on the replayed model.
    enrichment_comparison <- if (length(recorded_enrichment)) integrationReplayEnrichment(manager, recorded_enrichment,
        enrichment_mode, enrichment_api_key, enrichment_query, enrichment_download, tolerance)
    enrichment_api_key <- NULL
    # KIT-148: read from the replay's store before it is closed.
    replayed <- integrationModelComparable(manager, run$value)
    list(replayed = replayed, results = results, enrichment_comparison = enrichment_comparison,
        comparison = if (!is.null(recorded)) integrationCompareModels(recorded, replayed, tolerance),
        versions = versions, parameters = manager$.context$snapshot$exports[[INTEGRATION_MODEL_OPERATION]]$parameters,
        outputs_comparison = if (!is.null(recorded_outputs) && !is.null(results))
            integrationCompareOutputs(recorded_outputs, results, tolerance))
}

#' Reopen a saved integration project and repeat its analysis from its record (KIT-142): the project's store is
#' reopened, its record built and its latest model read, then `replayIntegrationRecord()` compares the repeated model
#' with it.
#' @param project_root The integration store's folder (the one holding `.workflow-entities`).
#' @param scratch A folder for the replay's store.
#' @param ... Passed to `replayIntegrationRecord()`.
#' @return `replayIntegrationRecord()`'s value with the record (`record`).
#' @noRd
replayIntegrationProject <- function(project_root, scratch, ...) {
    identity <- readSavedStoreIdentity(file.path(project_root, ".workflow-entities", "store.json"))
    manager <- newWorkflowStateOwner(newWorkflowEntityStore(project_root, identity$project_id, identity$workflow_id), reopen = TRUE)
    opened <- TRUE
    on.exit(if (opened) closeWorkflowStateOwner(manager), add = TRUE)
    record <- buildWorkflowAnalysisRecord(manager)
    export <- manager$.context$snapshot$exports[[INTEGRATION_MODEL_OPERATION]]
    # KIT-148: the recorded model read from its store before the store is closed.
    recorded <- if (!is.null(export)) integrationModelComparable(manager, withWorkflowStateScope(manager, function()
        hydrateSavedWorkflowEntity(manager$.context$store, export$reference)$value))
    # KIT-152: the further figures' recorded tables, read before the store is closed.
    recorded_outputs <- integrationRecordedOutputs(manager)
    # KIT-141: the recorded enrichments (products, rankings, STRING answers, KEGG/Reactome tables), read likewise.
    recorded_enrichment <- integrationRecordedEnrichment(manager)
    closeWorkflowStateOwner(manager)
    opened <- FALSE
    c(replayIntegrationRecord(record, scratch, recorded = recorded, recorded_outputs = recorded_outputs,
        recorded_enrichment = recorded_enrichment, ...), list(record = record))
}
