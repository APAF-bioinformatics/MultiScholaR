# KIT-184 wording list: every user-facing text of this ticket, verbatim (DRAFTS for the owner's wording pass)

**State (2026-11-25, after the owner's decisions VAL-160 and VAL-178).** The adjustment for removed RUV factors ships
with a narrowed claim and a stated limit, and it is off unless the user switches it on. The simulation's numbers are in
the package's texts; no placeholder is left (a test finds none). Every text below is a draft (VAL-41: implemented as
proposed; the owner edits before release). `%d` and `%s` stand for the numbers and names the code fills in. Extracted
from the sources by a script, not retyped; in three notes that the code fills with `sprintf()` (marked I06 (f)) the
source writes each per cent sign doubled, and they are shown here as the user sees them.

## Owner sets final wording: the items of the implementation audit (AUDIT-KIT-184-impl-20261125)

The facts below are required by VAL-160 and VAL-178; the words are the owner's. Each is applied as the audit register
recommends, or, where it gives no text, in the smallest form that states the fact.

| Item | What was applied | Where |
|---|---|---|
| I01, VAL-178 (1): the claim no wider than what was checked | the register's replacement sentences: "of the kinds the check decided on" with the kinds named; sets made only of the most variable features up to 9.4 % (6.8 % without the adjustment, 9.8 % with the true factors). The register's sentences do not hold the planted sets' 8.6 % that VAL-178 (1) asks to be printed: added in my words ("sets in which half the features truly followed the main factor ... up to 8.6 %", up to 7.0 % without the adjustment, up to 8.5 % with the true factors) | `methods_enrich_limits_ruv`, `enrich_adjust_ruv_help`, the starter workbook, NEWS |
| I02: the limit in the record of any model with an adjusted entry | the register's opening of the clause ("The correlation-adjusted tests can add to their model ..."); the clause is printed when the model was trained on corrected values or any KEGG or Reactome entry is adjusted | `methods_enrich_limits_ruv` |
| I03: an approximate set of factors | the register's fragment, with the share; and in the clause "where no run was removed after the correction". The register words that second change against the old opening; with I02's new opening it reads "for a model trained on RUV-corrected values and features measured in every run where no run was removed after the correction, this is identical to ..." | `methods_test_covariates_ruv_approximate`, `methods_enrich_limits_ruv` |
| I04: what is exact is the regression | the register's sentence | `enrich_values_ruv_covariates_note` |
| I05: where it found more true terms | the register's parenthesis; joined as "; it found more true terms than off (...)" | `enrich_adjust_ruv_help`, the starter workbook |
| I06 (a): two compositions, not a range | "proteins 11.2 % (sets of 10 and 30 features) and 13.3 % (sets of 100 included)"; "proteins 7.3 % (sets of 10 and 30 features) and 11.2 % (sets of 100 included)" | `methods_enrich_limits_ruv` |
| I06 (b): the correction that was simulated | the register's parenthesis | `methods_enrich_limits_ruv` |
| I06 (c): "1 factors" | "factor(s)" in the label and in the phrases that count the removed factors ("%d removed unwanted factor(s)", "%d removed factor(s)"). **Not changed:** the running sentences that say "the %d unwanted factors" (four notes, `methods_test_covariates_ruv`, `input_removed_factors`, `enrich_removed_approximate`): at one removed factor they still read "the 1 unwanted factors". A true singular needs a small change in the code that fills the texts; the owner may prefer that | `enrich_adjust_label_yes`, `enrich_test_covariates_design`, `enrich_test_covariates_counts`, `enrich_test_few_df`, the three `enrich_ranking_value_*_ruv`, the too-few-samples refusal, NEWS |
| I06 (d): the limit as a finding of that check | "In that check it did not protect ..." | the starter workbook, NEWS |
| I06 (e): both NEWS entries listed | below | this file |
| I06 (f): "(see Methods)" in the tab | replaced by the limit's numbers, my words: "... can still be called when that variation follows the factor: in a check on simulated data where it followed the groups, such sets were called on the model's main factor in 8.1 to 29.6 % of cases with the adjustment and in 2.0 to 11.2 % without it." | `enrich_adjust_ruv_help` and the three notes of an adjusted result |

Still as first proposed, and worth the owner's eye: `enrich_values_ruv_note` (the pre-ranked and permutation tests on
a model trained on RUV-corrected values) says the remaining variation "can bias the result"; the manual's sentence for
`adjust_ruv_factors` says "the adjustment does not protect pathways ..." as a general fact (the audit's I06 (d) names
NEWS and the starter only).

## Where each of the simulation's numbers stands (`results.md`, `s/numbers.tsv`, `s/error-rates.tsv`)

| Number | Value | In which text |
|---|---|---|
| [R3], [I3] | 5.4 %; 4.5 to 6.4 % | `methods_enrich_limits_ruv`; 5.4 % also in `enrich_adjust_ruv_help`, the starter workbook and NEWS |
| [R5], [I5] | 29.6 %; 27.8 to 31.3 % | `methods_enrich_limits_ruv` (genes, with the adjustment); 29.6 % as the top of "8.1 to 29.6 %" in the notes, the help, the starter and NEWS |
| [R6], [I6] | 11.2 %; 9.8 to 12.6 % | `methods_enrich_limits_ruv` (proteins, without the adjustment, sets of 100 included); 11.2 % as the top of "2.0 to 11.2 %" elsewhere |
| [R5m] | 13.3 % | `methods_enrich_limits_ruv` (proteins, with the adjustment, sets of 100 included) |
| [PF], [PC] | 13 to 32 %; 30 to 48 % | `enrich_floor_covariates_note` |
| other null sets against the independent covariate (VAL-178 (1)) | sets made only of the most variable features: up to 9.4 % (6.8 % without the adjustment, 9.8 % with the true factors); sets in which half the features truly followed the main factor: up to 8.6 % (up to 7.0 %, up to 8.5 %) | `methods_enrich_limits_ruv`; the first two figures of each also in the help, the starter and NEWS |
| the limit, per data set | with the adjustment 29.6 % (genes), 11.2 % and 13.3 % (proteins), 8.1 % (metabolites); without it 2.0 %, 7.3 % and 11.2 %, 10.1 %; with the true factors 1.6 to 2.8 %; 4.1 % where the unwanted variation was unrelated to the groups | `methods_enrich_limits_ruv` |

## The texts

### Tab (`R/utils_integration_gui.R`)

`enrich_adjust_ruv`

> Adjust for removed RUV factors

`enrich_adjust_ruv_help` *(changed after the implementation audit, I01, I05, I06 (f))*

> Off unless you switch it on. On: each feature is tested for what it shares with the factor beyond the unwanted factors that the RUV correction removed from this data set. In a check on simulated data this kept false calls within the checked limit (at most 5.4 %) for the kinds of sets with no link to what was tested that the check decided on; sets made only of the most variable features were called more often (up to 9.4 %, 6.8 % without it), and so were sets in which half the features truly followed the main factor, tested against an unrelated covariate (up to 8.6 %, up to 7.0 % without it); it found more true terms than off (with as many factors removed as were planted; with four removed where two were planted it found fewer). Because the unwanted factors are estimated, sets of features that share them can still be called when that variation follows the factor: in a check on simulated data where it followed the groups, such sets were called on the model's main factor in 8.1 to 29.6 % of cases with the adjustment and in 2.0 to 11.2 % without it. Off: the values before RUV are tested without those factors. On needs a data set to which the omic applied RUV. The setting does not apply to the pre-ranked test, the permutation test or STRING.

`enrich_adjust_label_yes` *(changed after the implementation audit, I06 (c))*

> Adjusted for removed RUV factors: yes (%d factor(s)).

`enrich_adjust_label_no`

> Adjusted for removed RUV factors: no.

`enrich_values_ruv_covariates_note` *(changed after the implementation audit, I04, I06 (f))*

> This model was trained on values corrected for unwanted variation (RUV). The test reads the values before that correction and puts the %d unwanted factors the correction removed from this data set into its model, so a feature counts only for what it shares with the factor beyond those unwanted factors. For features measured in every run the regression is exactly that of the same model (the factor and these %d unwanted factors) fitted to the corrected values, with %d residual degrees of freedom; the test's variance trend reads each feature's average before the correction, so its p-values can differ where the correction moved a feature's level. It is not the test of the corrected values on the factor alone. The unwanted factors were estimated from these same data and are treated as known. Because the unwanted factors are estimated, sets of features that share them can still be called when that variation follows the factor: in a check on simulated data where it followed the groups, such sets were called on the model's main factor in 8.1 to 29.6 % of cases with the adjustment and in 2.0 to 11.2 % without it.

`enrich_values_ruv_covariates_note_approximate` *(changed after the implementation audit, I06 (f))*

> This model was trained on values corrected for unwanted variation (RUV). The test reads the values before that correction and puts the %d unwanted factors the correction removed from this data set into its model, so a feature counts only for what it shares with the factor beyond those unwanted factors. The model has %d residual degrees of freedom. It is not the test of the corrected values on the factor alone. The unwanted factors were estimated from these same data and are treated as known. Because the unwanted factors are estimated, sets of features that share them can still be called when that variation follows the factor: in a check on simulated data where it followed the groups, such sets were called on the model's main factor in 8.1 to 29.6 % of cases with the adjustment and in 2.0 to 11.2 % without it.

`enrich_removed_approximate`

> Runs were removed after the correction, or it removed variation outside these %d factors: %s of what it removed lies within them, so the match with the corrected values is approximate for every feature.

`enrich_values_adjusted_default_note` *(changed after the implementation audit, I06 (f))*

> This model was trained on the values before RUV correction. For this result the adjustment for removed RUV factors was switched on: the test puts the %d unwanted factors that the omic's RUV correction removed from this data set into its model, so a feature counts only for what it shares with the factor beyond those unwanted factors. If the factor itself follows them, little is left to test (see the share below). The unwanted factors were estimated from these same data and are treated as known. Because the unwanted factors are estimated, sets of features that share them can still be called when that variation follows the factor: in a check on simulated data where it followed the groups, such sets were called on the model's main factor in 8.1 to 29.6 % of cases with the adjustment and in 2.0 to 11.2 % without it.

`enrich_values_ruv_unadjusted_note` *(changed for VAL-160)*

> This model was trained on RUV-corrected values. This result is the form without the adjustment for removed RUV factors, which is the default (the adjustment is optional): it tests the normalised values before RUV correction without the unwanted factors the correction removed. Its degrees of freedom are counted correctly, but variation the correction removed stays in the values, which lowers power. In a check on simulated data where the unwanted variation followed the groups, sets of features sharing it were called on the model's main factor in 2.0 to 11.2 % of cases in this form and in 8.1 to 29.6 % with the adjustment.

`enrich_test_covariates_design` *(changed after the implementation audit, I06 (c))*

> Design: intercept, factor and %d removed unwanted factor(s); %d biological samples, %d residual degrees of freedom.

`enrich_test_covariates_runs`

> %d of the %d features used were measured in every run; for the others the correction worked feature by feature and the match with the corrected values is approximate.

`enrich_test_covariates_share`

> The removed unwanted factors explain %s of this factor's variation across these samples (about %s would be expected for an unrelated factor).

`enrich_test_covariates_counts` *(changed after the implementation audit, I06 (c))*

> Design: intercept, factor and %d removed unwanted factor(s); %d biological samples, %d residual degrees of freedom. %d of the %d features used were measured in every run; for the others the correction worked feature by feature and the match with the corrected values is approximate. The removed unwanted factors explain %s of this factor's variation across these samples (about %s would be expected for an unrelated factor).

`enrich_factor_removed_overlap`

> More than half of this factor's variation across the samples of this data set follows the unwanted factors removed from it. Little is left to test here, and in this data set the factor may reflect variation the correction treated as unwanted.

`enrich_test_few_df` *(changed after the implementation audit, I06 (c))*

> With %d removed factor(s) in the model this test has %d residual degrees of freedom, fewer than the smallest number (%d) at which its error rate was checked.

`enrich_floor_covariates_note` *(changed for VAL-160)*

> False calls are at most those of the correlation-adjusted test with the same unwanted factors (counting negative correlations as zero never makes a p-value smaller). In a check on simulated data with 12 biological samples and 2 removed factors, terms of 10 features in which half the features truly followed the main factor were found in 13 to 32 % of cases with these factors in the model (the correlation-adjusted test: 30 to 48 %).

`enrich_values_ruv_note`

> This model was trained on RUV-corrected values, but this result uses the normalised values before RUV correction without the unwanted factors the correction removed; variation the correction removed stays in them, which lowers power and, if it follows the factor, can bias the result.

### Record and Methods (`R/utils_integration_record.R`)

`enrich_ranking_value_camera_ruv` *(changed after the implementation audit, I06 (c))*

> regression of each feature on the factor and %d removed unwanted factor(s) over %d biological samples (%d residual degrees of freedom; features measured in all samples: %d); correlation-adjusted test (camera, correlation estimated)

`enrich_ranking_value_camera_pr_ruv` *(changed after the implementation audit, I06 (c))*

> regression of each feature on the factor and %d removed unwanted factor(s) over %d biological samples (%d residual degrees of freedom; features measured in all samples: %d); pre-ranked test with an assumed correlation of 0.01 (cameraPR)

`enrich_ranking_value_camera_floor_ruv` *(changed after the implementation audit, I06 (c))*

> regression of each feature on the factor and %d removed unwanted factor(s) over %d biological samples (%d residual degrees of freedom; features measured in all samples: %d); correlation-adjusted test, negative correlations counted as zero (camera, correlation estimated, allow.neg.cor = FALSE)

`methods_test_covariates_ruv`

> an intercept and the %d unwanted factors that the RUV correction removed from this data set (a basis of the difference between the normalised and the corrected values over the %d features measured in every run, averaged per biological sample; SHA-256 %s; %d residual degrees of freedom)

`methods_test_covariates_ruv_approximate` *(changed after the implementation audit, I03, new)*

> ; approximate: %s of what the correction removed lies within these factors (runs were removed after the correction), so they do not match the corrected values exactly

`methods_test_covariates_ruv_off` *(changed for VAL-160)*

> an intercept and no other covariate (without the adjustment for the removed RUV factors, which is optional: variation the correction removed stays in the values tested)

`methods_test_covariates_model_uncorrected`

> ; the model was trained on the values before RUV, and the adjustment was switched on for this result

`input_removed_factors`

> ; %d unwanted factors recorded as removed from %s (%s; SHA-256 %s)

`methods_enrich_limits_ruv` *(changed after the implementation audit, I01, I02, I03, I06 (a), (b))*

> The correlation-adjusted tests can add to their model the unwanted factors the RUV correction removed from each data set (the adjustment for removed RUV factors, a setting that is off unless chosen); for a model trained on RUV-corrected values and features measured in every run where no run was removed after the correction, this is identical to fitting that model (the factor and the removed factors) to the corrected values; those factors were estimated from the same data and are treated as known. Each result names the form used: a result without the adjustment tests the values before RUV as they are (its degrees of freedom are counted correctly; variation the correction removed stays in the values tested), and a result with it on a model trained on the values before RUV counts a feature for what it shares with the factor beyond the removed factors. The pre-ranked and permutation tests do not use them: variation the correction removed stays in the values they rank. In a check on simulated data with 12 biological samples and 2 removed factors (RUV-III with the groups as replicate sets and every feature a control; the legacy routes' choice of controls and of the number of factors was not simulated), with the removed factors in the model, sets with no link to what was tested, of the kinds the check decided on (random sets, correlated blocks of features, sets of low-abundance features and sets of features sharing unwanted variation, against a covariate drawn independently of the data; the first three also against the model's main factor), were called at the 5 % level in at most 5.4 % (95 % interval of that largest figure 4.5 to 6.4 %) of cases; against such a covariate, sets made only of the most variable features, those the model selects, were called in up to 9.4 % (6.8 % without the adjustment, 9.8 % had the true unwanted factors been known), a property of the correlation-adjusted test on such sets, and sets in which half the features truly followed the main factor in up to 8.6 % (up to 7.0 % without the adjustment, up to 8.5 % had the true unwanted factors been known); true terms were found more often than without the adjustment. Sets of features that share unwanted variation and have no link to the simulated biology were called on the model's main factor in at most 4.1 % of cases with the adjustment when that variation was unrelated to the groups, and, when it followed the groups (correlation 0.6), in 8.1 to 29.6 % with the adjustment (genes 29.6 %, 95 % interval 27.8 to 31.3 %; proteins 11.2 % (sets of 10 and 30 features) and 13.3 % (sets of 100 included); metabolites 8.1 %; more often for larger sets), 2.0 to 11.2 % without it (genes 2.0 %; proteins 7.3 % (sets of 10 and 30 features) and 11.2 % (sets of 100 included), 95 % interval of the larger figure 9.8 to 12.6 %; metabolites 10.1 %), and 1.6 to 2.8 % had the true unwanted factors been known.

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

`integration_enrichment_too_few_samples_for_test_covariates` *(changed after the implementation audit, I06 (c))*

> With %d removed unwanted factor(s) in its model, the correlation-adjusted test needs at least %d biological samples in this data set; it has %d.

`integration_enrichment_factor_collinear`

> This factor is a combination of the unwanted factors removed from %s, so nothing is left to test once they are in the model.

### The setting's line in `integrationSettings()$enrichment`

> no (optional: yes adds the unwanted factors the RUV correction removed to the test's model)

### `integrationEnrichFactor()`'s documentation of the setting (`R/utils_integration_workbook.R`, `man/integrationEnrichFactor.Rd`)

> `adjust_ruv_factors` (`TRUE` or `FALSE`), "Adjust for removed RUV factors", applies to `"camera"` and `"camera_floor"`: with `TRUE` the test puts the unwanted factors that the omic's RUV correction removed from the data set into its model, so a feature counts for what it shares with the factor beyond them. Left out, it is `FALSE` whatever the model was trained on: the adjustment is optional. Because the removed factors are estimated, the adjustment does not protect pathways whose features share unwanted variation that follows the factor (the Methods text of a record gives the numbers). `TRUE` needs a data set to which the omic applied RUV, saved by this version or later. Every result says which form it is. The pre-ranked test, the permutation test and STRING do not use the setting.

### The starter workbook's paragraph (`inst/workbooks/integration/integration_workflow_starter.rmd`) *(changed after the implementation audit, I01, I05, I06 (d))*

> **Adjust for removed RUV factors (optional, off unless you choose it).** The correlation-adjusted tests can put the unwanted factors that an omic's RUV correction removed from a data set into their model, so a feature counts for what it shares with the factor beyond those unwanted factors: `enrichment_overrides <- list(adjust_ruv_factors = TRUE)`. For a model trained on RUV-corrected values (`input_matrix <- "ruv_corrected"` above) and features measured in every run, that is the same as fitting the test's model to the corrected values. In a check on simulated data the adjusted test kept false calls within the checked limit (at most 5.4 %) for the kinds of sets with no link to what was tested that the check decided on; sets made only of the most variable features were called more often (up to 9.4 %, 6.8 % without it), and so were sets in which half the features truly followed the main factor, tested against an unrelated covariate (up to 8.6 %, up to 7.0 % without it); it found more true terms (with as many factors removed as were planted; with four removed where two were planted it found fewer). In that check it did not protect pathways whose features share unwanted variation that follows the groups: such sets were called in 8.1 to 29.6 % of cases with the adjustment and in 2.0 to 11.2 % without it. Each result says whether it was adjusted. Without the setting every model is tested on the values before RUV without those factors. An omic saved before this version must be saved again at its Summary step before the adjusted test can run. The pre-ranked test, the permutation test and STRING do not use the setting.

### NEWS, under the unreleased heading, New features: the first entry

> Integration inputs record the unwanted factors a correction for unwanted variation (RUV) removed: where an omic applied RUV to a data set, saving its Summary writes one more file, `removed_<data set>.parquet`, beside the normalised and the corrected matrix (a basis of the difference between the two over the features measured in every run, checked against the correction's number of factors; contract v2, an addition that earlier versions read as before). An input saved earlier must be saved again at the omic's Summary step before the correlation-adjusted enrichment test can run with the adjustment for removed RUV factors (KIT-184).

### NEWS, under the unreleased heading, New features: the second entry *(changed after the implementation audit, I01, I06 (c), (d))*

> Integration enrichment, an optional adjustment for removed RUV factors: the setting "Adjust for removed RUV factors" (`adjust_ruv_factors`), in the Enrichment tab and in `integrationEnrichFactor(overrides = )`, off unless chosen. Switched on, the correlation-adjusted tests (the default test and the form with negative correlations counted as zero) put the unwanted factors the RUV correction removed from the data set into their model (intercept, factor, removed factors), so a feature counts for what it shares with the factor beyond those factors; for a model trained on RUV-corrected values and features measured in every run this is identical to fitting that model to the corrected values. In a check on simulated data (12 biological samples, 2 removed factors) the adjusted test called sets with no link to what was tested, of the kinds the check decided on, in at most 5.4 % of cases (95 % interval of that largest figure 4.5 to 6.4 %; sets made only of the most variable features: up to 9.4 %, against 6.8 % without the adjustment; sets in which half the features truly followed the main factor, tested against an unrelated covariate: up to 8.6 %, against up to 7.0 %) and found more true terms than the test without it. In that check it did not protect sets of features that share unwanted variation when that variation follows the groups: such sets were called on the model's main factor in 8.1 to 29.6 % of cases with the adjustment, against 2.0 to 11.2 % without it and 1.6 to 2.8 % had the true unwanted factors been known. No result changes unless the setting is switched on: every model, one trained on RUV-corrected values included, is tested as before on the values before RUV without the factors. Every result of a correlation-adjusted test now carries the line "Adjusted for removed RUV factors: yes (n factor(s))." or "... no." in the tab, in `integrationEnrichFactor()`'s summaries and in records. The adjusted test is refused when fewer than max(4, the smallest number of samples per feature) biological samples are left after one per removed factor, when the factor is a combination of the removed factors, and for a data set without usable removed factors. The pre-ranked test, the permutation test and STRING do not use the removed factors. A recorded result replays in the form it was saved in; a result recorded before this version replays without the adjustment (VAL-53, VAL-58, VAL-62, VAL-160, VAL-178; KIT-184).
