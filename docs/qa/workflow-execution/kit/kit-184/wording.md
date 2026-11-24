# KIT-184 wording list: every user-facing text of this ticket, verbatim (DRAFTS for the owner's wording pass)

**State (2026-11-24, after the owner's decision VAL-160).** The adjustment for removed RUV factors ships with corrected
claims and a stated limit, and it is off unless the user switches it on. The simulation's nine numbers are written into
the package's texts; no placeholder is left (a test finds none: `test-integration-removed-factors.R`, "no placeholder
of the simulation's numbers is left in a shipped text"). Every text below is a draft (VAL-41: implemented as proposed;
the owner edits before release). The texts marked *changed for VAL-160* were rewritten after the decision, using the
results audit's replacement wordings (`tickets/kit/audits/AUDIT-KIT-184-declaration-20261124.md`, K184-R12) and adapted
only as far as the default being off requires; the others are as first proposed. `%d` and `%s` stand for the numbers and
names the code fills in. Extracted from the sources by a script, not retyped.

## Where each of the nine numbers stands (`results.md`, `s/numbers.tsv`)

| Number | Value | In which text |
|---|---|---|
| [R3], [I3] | 5.4 %; 4.5 to 6.4 % | `methods_enrich_limits_ruv` ("in at most 5.4 % (95 % interval of that largest figure 4.5 to 6.4 %)"); 5.4 % also in `enrich_adjust_ruv_help`, the starter workbook and NEWS |
| [R5], [I5] | 29.6 %; 27.8 to 31.3 % | `methods_enrich_limits_ruv` (genes, with the adjustment); 29.6 % as the top of the range "8.1 to 29.6 %" in `enrich_values_ruv_unadjusted_note`, the starter and NEWS |
| [R6], [I6] | 11.2 %; 9.8 to 12.6 % | `methods_enrich_limits_ruv` (proteins, without the adjustment: "7.3 to 11.2 %, 95 % interval of the larger figure 9.8 to 12.6 %"); 11.2 % as the top of the range "2.0 to 11.2 %" in the note, the starter and NEWS |
| [R5m] | 13.3 % | `methods_enrich_limits_ruv` (proteins, with the adjustment: "11.2 to 13.3 %") |
| [PF], [PC] | 13 to 32 %; 30 to 48 % | `enrich_floor_covariates_note` |

Other figures of `results.md` used in the limit, as VAL-160 words it: 4.1 % (the largest rate with the adjustment where
the unwanted variation was unrelated to the groups); with the adjustment where it followed the groups 29.6 % (genes),
11.2 and 13.3 % (proteins, by the set sizes counted), 8.1 % (metabolites); without it 2.0 %, 7.3 and 11.2 %, 10.1 %; with
the true unwanted factors 1.6 to 2.8 %.

## For the owner's wording pass: what to look at first

- `methods_enrich_limits_ruv` is long because it carries every number once. It says: the adjustment is optional and off
  unless chosen; what it kept (at most 5.4 %, and more true terms); the limit with its condition and the numbers per
  data set.
- The limit sentence ("Because the unwanted factors are estimated, sets of features that share them can still be
  called when that variation follows the factor (see Methods).") stands in four texts: the control's help and the note of
  every adjusted result (a model trained on corrected values, exact and approximate; a model trained on the values
  before RUV). The results audit named two of them; I added it to the other two notes of an adjusted result so that no
  adjusted result is shown without the limit.
- `enrich_values_ruv_unadjusted_note` now names the form without the adjustment as the default and gives both ranges; it
  no longer says "a sensitivity analysis" or "can bias a set".
- The data sets are called genes, proteins and metabolites in the Methods sentence (the simulation's transcriptomics,
  proteomics and metabolomics data sets).
- `enrich_values_ruv_note` (the pre-ranked and permutation tests on a model trained on RUV-corrected values) still says
  the remaining variation "can bias the result"; the results audit did not ask for a change there.

## The texts

### Tab (`R/utils_integration_gui.R`)

`enrich_adjust_ruv`

> Adjust for removed RUV factors

`enrich_adjust_ruv_help` *(changed for VAL-160)*

> Off unless you switch it on. On: each feature is tested for what it shares with the factor beyond the unwanted factors that the RUV correction removed from this data set. In a check on simulated data this kept false calls within the checked limit for sets with no link to what was tested (at most 5.4 %) and found more true terms than off. Because the unwanted factors are estimated, sets of features that share them can still be called when that variation follows the factor (see Methods). Off: the values before RUV are tested without those factors. On needs a data set to which the omic applied RUV. The setting does not apply to the pre-ranked test, the permutation test or STRING.

`enrich_adjust_label_yes`

> Adjusted for removed RUV factors: yes (%d factors).

`enrich_adjust_label_no`

> Adjusted for removed RUV factors: no.

`enrich_values_ruv_covariates_note` *(changed for VAL-160)*

> This model was trained on values corrected for unwanted variation (RUV). The test reads the values before that correction and puts the %d unwanted factors the correction removed from this data set into its model, so a feature counts only for what it shares with the factor beyond those unwanted factors. For features measured in every run this gives exactly the result of fitting the same model (the factor and these %d unwanted factors) to the corrected values, with %d residual degrees of freedom. It is not the test of the corrected values on the factor alone. The unwanted factors were estimated from these same data and are treated as known. Because the unwanted factors are estimated, sets of features that share them can still be called when that variation follows the factor (see Methods).

`enrich_values_ruv_covariates_note_approximate` *(changed for VAL-160)*

> This model was trained on values corrected for unwanted variation (RUV). The test reads the values before that correction and puts the %d unwanted factors the correction removed from this data set into its model, so a feature counts only for what it shares with the factor beyond those unwanted factors. The model has %d residual degrees of freedom. It is not the test of the corrected values on the factor alone. The unwanted factors were estimated from these same data and are treated as known. Because the unwanted factors are estimated, sets of features that share them can still be called when that variation follows the factor (see Methods).

`enrich_removed_approximate`

> Runs were removed after the correction, or it removed variation outside these %d factors: %s of what it removed lies within them, so the match with the corrected values is approximate for every feature.

`enrich_values_adjusted_default_note` *(changed for VAL-160)*

> This model was trained on the values before RUV correction. For this result the adjustment for removed RUV factors was switched on: the test puts the %d unwanted factors that the omic's RUV correction removed from this data set into its model, so a feature counts only for what it shares with the factor beyond those unwanted factors. If the factor itself follows them, little is left to test (see the share below). The unwanted factors were estimated from these same data and are treated as known. Because the unwanted factors are estimated, sets of features that share them can still be called when that variation follows the factor (see Methods).

`enrich_values_ruv_unadjusted_note` *(changed for VAL-160)*

> This model was trained on RUV-corrected values. This result is the form without the adjustment for removed RUV factors, which is the default (the adjustment is optional): it tests the normalised values before RUV correction without the unwanted factors the correction removed. Its degrees of freedom are counted correctly, but variation the correction removed stays in the values, which lowers power. In a check on simulated data where the unwanted variation followed the groups, sets of features sharing it were called on the model's main factor in 2.0 to 11.2 % of cases in this form and in 8.1 to 29.6 % with the adjustment.

`enrich_test_covariates_design`

> Design: intercept, factor and %d removed unwanted factors; %d biological samples, %d residual degrees of freedom.

`enrich_test_covariates_runs`

> %d of the %d features used were measured in every run; for the others the correction worked feature by feature and the match with the corrected values is approximate.

`enrich_test_covariates_share`

> The removed unwanted factors explain %s of this factor's variation across these samples (about %s would be expected for an unrelated factor).

`enrich_test_covariates_counts`

> Design: intercept, factor and %d removed unwanted factors; %d biological samples, %d residual degrees of freedom. %d of the %d features used were measured in every run; for the others the correction worked feature by feature and the match with the corrected values is approximate. The removed unwanted factors explain %s of this factor's variation across these samples (about %s would be expected for an unrelated factor).

`enrich_factor_removed_overlap`

> More than half of this factor's variation across the samples of this data set follows the unwanted factors removed from it. Little is left to test here, and in this data set the factor may reflect variation the correction treated as unwanted.

`enrich_test_few_df`

> With %d removed factors in the model this test has %d residual degrees of freedom, fewer than the smallest number (%d) at which its error rate was checked.

`enrich_floor_covariates_note` *(changed for VAL-160)*

> False calls are at most those of the correlation-adjusted test with the same unwanted factors (counting negative correlations as zero never makes a p-value smaller). In a check on simulated data with 12 biological samples and 2 removed factors, terms of 10 features in which half the features truly followed the main factor were found in 13 to 32 % of cases with these factors in the model (the correlation-adjusted test: 30 to 48 %).

`enrich_values_ruv_note`

> This model was trained on RUV-corrected values, but this result uses the normalised values before RUV correction without the unwanted factors the correction removed; variation the correction removed stays in them, which lowers power and, if it follows the factor, can bias the result.

### Record and Methods (`R/utils_integration_record.R`)

`enrich_ranking_value_camera_ruv`

> regression of each feature on the factor and %d removed unwanted factors over %d biological samples (%d residual degrees of freedom; features measured in all samples: %d); correlation-adjusted test (camera, correlation estimated)

`enrich_ranking_value_camera_pr_ruv`

> regression of each feature on the factor and %d removed unwanted factors over %d biological samples (%d residual degrees of freedom; features measured in all samples: %d); pre-ranked test with an assumed correlation of 0.01 (cameraPR)

`enrich_ranking_value_camera_floor_ruv`

> regression of each feature on the factor and %d removed unwanted factors over %d biological samples (%d residual degrees of freedom; features measured in all samples: %d); correlation-adjusted test, negative correlations counted as zero (camera, correlation estimated, allow.neg.cor = FALSE)

`methods_test_covariates_ruv`

> an intercept and the %d unwanted factors that the RUV correction removed from this data set (a basis of the difference between the normalised and the corrected values over the %d features measured in every run, averaged per biological sample; SHA-256 %s; %d residual degrees of freedom)

`methods_test_covariates_ruv_off` *(changed for VAL-160)*

> an intercept and no other covariate (without the adjustment for the removed RUV factors, which is optional: variation the correction removed stays in the values tested)

`methods_test_covariates_model_uncorrected`

> ; the model was trained on the values before RUV, and the adjustment was switched on for this result

`input_removed_factors`

> ; %d unwanted factors recorded as removed from %s (%s; SHA-256 %s)

`methods_enrich_limits_ruv` *(changed for VAL-160)*

> For a model trained on RUV-corrected values the correlation-adjusted tests can add to their model the unwanted factors the correction removed from each data set (the adjustment for removed RUV factors, a setting that is off unless chosen), which for features measured in every run is identical to fitting that model (the factor and the removed factors) to the corrected values; those factors were estimated from the same data and are treated as known. Each result names the form used: a result without the adjustment tests the values before RUV as they are (its degrees of freedom are counted correctly; variation the correction removed stays in the values tested), and a result with it on a model trained on the values before RUV counts a feature for what it shares with the factor beyond the removed factors. The pre-ranked and permutation tests do not use them: variation the correction removed stays in the values they rank. In a check on simulated data with 12 biological samples and 2 removed factors, with the removed factors in the model, sets with no link to what was tested (any set against a covariate drawn independently of the data; sets of features without unwanted variation against the model's main factor) were called at the 5 % level in at most 5.4 % (95 % interval of that largest figure 4.5 to 6.4 %) of cases, and true terms were found more often than without the adjustment. Sets of features that share unwanted variation and have no link to the simulated biology were called on the model's main factor in at most 4.1 % of cases with the adjustment when that variation was unrelated to the groups, and, when it followed the groups (correlation 0.6), in 8.1 to 29.6 % with the adjustment (genes 29.6 %, 95 % interval 27.8 to 31.3 %; proteins 11.2 to 13.3 %; metabolites 8.1 %; more often for larger sets), 2.0 to 11.2 % without it (genes 2.0 %; proteins 7.3 to 11.2 %, 95 % interval of the larger figure 9.8 to 12.6 %; metabolites 10.1 %), and 1.6 to 2.8 % had the true unwanted factors been known.

`methods_enrich_limits_family`

> P-values were adjusted within each factor, data set and database (and Gene Ontology aspect), not across them.

`methods_enrich_limits_end`

> The correlation-adjusted and permutation tests use only features measured in every biological sample, which favours abundant features and leaves out features that fall below detection in some samples, often those that change most. The tests used read each data set's normalised values before RUV correction, the values the model is trained on by default.

### Refusals (`R/utils_integration_enrichment.R`, ticket 8.4)

`integration_enrichment_removed_factors_missing` *(changed for VAL-160)*

> 'Adjust for removed RUV factors' was switched on, and it needs the unwanted factors that the correction removed from %s. The saved data of %s does not hold them (it was saved by an earlier version). Save the %s Summary again, import the data and train the model again. Switched off (the default), the test runs as it is.

`integration_enrichment_removed_factors_missing_undetermined` *(changed for VAL-160)*

> 'Adjust for removed RUV factors' was switched on, and it needs the unwanted factors that the correction removed from %s. They could not be determined (%s), so the test cannot be run with the adjustment on this data set. Switched off (the default), the test runs as it is.

`integration_enrichment_removed_factors_missing_no_ruv`

> No correction for unwanted variation (RUV) was applied to %s in %s, so there are no removed factors to adjust for. Switch 'Adjust for removed RUV factors' off.

`integration_enrichment_too_few_samples_for_test_covariates`

> With %d removed unwanted factors in its model, the correlation-adjusted test needs at least %d biological samples in this data set; it has %d.

`integration_enrichment_factor_collinear`

> This factor is a combination of the unwanted factors removed from %s, so nothing is left to test once they are in the model.

Two refusal sentences were removed with the default (`..._missing_default`, `..._missing_default_undetermined`): a
refusal for missing factors now arises only when the setting was switched on, and reads the same on every model.

### The setting's line in `integrationSettings()$enrichment` *(changed for VAL-160)*

> no (optional: yes adds the unwanted factors the RUV correction removed to the test's model)

### `integrationEnrichFactor()`'s documentation of the setting (`R/utils_integration_workbook.R`, `man/integrationEnrichFactor.Rd`) *(changed for VAL-160)*

> `adjust_ruv_factors` (`TRUE` or `FALSE`), "Adjust for removed RUV factors", applies to `"camera"` and `"camera_floor"`: with `TRUE` the test puts the unwanted factors that the omic's RUV correction removed from the data set into its model, so a feature counts for what it shares with the factor beyond them. Left out, it is `FALSE` whatever the model was trained on: the adjustment is optional. Because the removed factors are estimated, the adjustment does not protect pathways whose features share unwanted variation that follows the factor (the Methods text of a record gives the numbers). `TRUE` needs a data set to which the omic applied RUV, saved by this version or later. Every result says which form it is. The pre-ranked test, the permutation test and STRING do not use the setting.

### The starter workbook's paragraph (`inst/workbooks/integration/integration_workflow_starter.rmd`) *(changed for VAL-160)*

> **Adjust for removed RUV factors (optional, off unless you choose it).** The correlation-adjusted tests can put the unwanted factors that an omic's RUV correction removed from a data set into their model, so a feature counts for what it shares with the factor beyond those unwanted factors: `enrichment_overrides <- list(adjust_ruv_factors = TRUE)`. For a model trained on RUV-corrected values (`input_matrix <- "ruv_corrected"` above) and features measured in every run, that is the same as fitting the test's model to the corrected values. In a check on simulated data the adjusted test kept false calls within the checked limit for sets with no link to what was tested (at most 5.4 %) and found more true terms. It does not protect pathways whose features share unwanted variation that follows the groups: such sets were called in 8.1 to 29.6 % of cases with the adjustment and in 2.0 to 11.2 % without it. Each result says whether it was adjusted. Without the setting every model is tested on the values before RUV without those factors. An omic saved before this version must be saved again at its Summary step before the adjusted test can run. The pre-ranked test, the permutation test and STRING do not use the setting.

### NEWS, under the unreleased heading, New features *(changed for VAL-160; it replaces the entry under Scientific changes, since no result changes unless the setting is switched on)*

> Integration enrichment, an optional adjustment for removed RUV factors: the setting "Adjust for removed RUV factors" (`adjust_ruv_factors`), in the Enrichment tab and in `integrationEnrichFactor(overrides = )`, off unless chosen. Switched on, the correlation-adjusted tests (the default test and the form with negative correlations counted as zero) put the unwanted factors the RUV correction removed from the data set into their model (intercept, factor, removed factors), so a feature counts for what it shares with the factor beyond those factors; for a model trained on RUV-corrected values and features measured in every run this is identical to fitting that model to the corrected values. In a check on simulated data (12 biological samples, 2 removed factors) the adjusted test called sets with no link to what was tested in at most 5.4 % of cases (95 % interval of that largest figure 4.5 to 6.4 %) and found more true terms than the test without it. It does not protect sets of features that share unwanted variation when that variation follows the groups: such sets were called on the model's main factor in 8.1 to 29.6 % of cases with the adjustment, against 2.0 to 11.2 % without it and 1.6 to 2.8 % had the true unwanted factors been known. No result changes unless the setting is switched on: every model, one trained on RUV-corrected values included, is tested as before on the values before RUV without the factors. Every result of a correlation-adjusted test now carries the line "Adjusted for removed RUV factors: yes (n factors)." or "... no." in the tab, in `integrationEnrichFactor()`'s summaries and in records. The adjusted test is refused when fewer than max(4, the smallest number of samples per feature) biological samples are left after one per removed factor, when the factor is a combination of the removed factors, and for a data set without usable removed factors. The pre-ranked test, the permutation test and STRING do not use the removed factors. A recorded result replays in the form it was saved in; a result recorded before this version replays without the adjustment (VAL-53, VAL-58, VAL-62, VAL-160; KIT-184).
