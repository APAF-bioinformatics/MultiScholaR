# KIT-184 wording list (VAL-41: implemented as proposed; the owner edits before release)

**State.** The package's texts still hold the placeholders [R3], [I3], [R5], [I5], [R5m], [R6], [I6], and
`enrich_floor_covariates_note` still says its power "has not been measured yet". They are written into the code, with
the tests that read them, only after the coordinator's word (the result was reported first; the full suite has not
been run). No build may be released or merged with the placeholders in place. Below: the keys, the numbers, the three
sentences as they would read with the numbers in, and what the result means for two of them.

## The keys

Tab (`R/utils_integration_gui.R`): `enrich_adjust_ruv`, `enrich_adjust_ruv_help`, `enrich_adjust_label_yes`,
`enrich_adjust_label_no`, `enrich_values_ruv_covariates_note`, `enrich_values_ruv_covariates_note_approximate`,
`enrich_removed_approximate`, `enrich_values_adjusted_default_note`, `enrich_values_ruv_unadjusted_note`,
`enrich_test_covariates_design`, `enrich_test_covariates_runs`, `enrich_test_covariates_share`,
`enrich_test_covariates_counts`, `enrich_factor_removed_overlap`, `enrich_test_few_df`,
`enrich_floor_covariates_note`, `enrich_values_ruv_note` (reworded: now for the pre-ranked and the permutation test
only).

Record and Methods (`R/utils_integration_record.R`): `enrich_ranking_value_camera_ruv`,
`enrich_ranking_value_camera_pr_ruv`, `enrich_ranking_value_camera_floor_ruv`, `methods_test_covariates_ruv`,
`methods_test_covariates_ruv_off`, `methods_test_covariates_model_uncorrected`, `input_removed_factors`,
`methods_enrich_limits_ruv`, `methods_enrich_limits_family`, `methods_enrich_limits_end` (shortened).

Refusals (`R/utils_integration_enrichment.R`, ticket 8.4): `integration_enrichment_removed_factors_missing` and its
forms `_undetermined`, `_no_ruv`, `_default`, `_default_undetermined`;
`integration_enrichment_too_few_samples_for_test_covariates`; `integration_enrichment_factor_collinear`.

## The numbers from the counted simulation (`results.md`, `s/numbers.tsv`)

[R3] 5.4 %, [I3] 4.5 to 6.4 %; [R5] 29.6 %, [I5] 27.8 to 31.3 %; [R6] 11.2 %, [I6] 9.8 to 12.6 %; [R5m] 13.3 %; [PF]
13 to 32 %; [PC] 30 to 48 %.

## The three sentences with the numbers in

**`methods_enrich_limits_ruv`** (its last two sentences): "In a check on simulated data with 12 biological samples and 2
removed factors, with the removed factors in the model, sets with no link to what was tested (any set against a
covariate drawn independently of the data; sets of features without unwanted variation against the model's main
factor) were called at the 5 % level in 5.4 % (95 % interval 4.5 to 6.4 %) of cases. Sets of features that share
unwanted variation and have no link to the simulated biology were called in 29.6 % (27.8 to 31.3 %) of cases on the
model's main factor, and in 11.2 % (9.8 to 12.6 %) of cases without the adjustment."

**`enrich_values_ruv_unadjusted_note`** (its last sentence): "In a check on simulated data, sets of features that share
unwanted variation and have no link to the simulated biology were called on the model's main factor in 11.2 % of cases
in this form and in 13.3 % with the adjustment, in the same simulated data set."

**`enrich_floor_covariates_note`**: "False calls are at most those of the correlation-adjusted test with the same
unwanted factors (counting negative correlations as zero never makes a p-value smaller). In a check on simulated data
with 12 biological samples and 2 removed factors, terms of 10 features in which half the features truly followed the
main factor were found in 13 to 32 % of cases with these factors in the model (the correlation-adjusted test: 30 to
48 %)."

## For the owner: what the result does to two of these sentences

The numbers are filled in exactly as the ticket defines them, and they are true. But both sentences were drafted on
the expectation that the adjustment protects sets of features that share unwanted variation, and the simulation found
the opposite where it matters:

- Where the unwanted factor follows the groups (the simulation's second cell), such sets were called on the main
  factor in 8.1 to 29.6 % of cases **with** the adjustment and in 2.0 to 11.2 % **without** it. Where it does not
  follow the groups, 0.9 to 4.1 % with and 0.0 % without. With the planted factors themselves in the model: 1.3 to
  2.8 %.
- So the sentence before the numbers in `enrich_values_ruv_unadjusted_note` ("... can bias a set where that variation
  lines up with the factor") describes the form without the adjustment as the one at risk, and the numbers that follow
  it show the adjusted form at least as exposed. The same holds for the last sentence of
  `methods_enrich_limits_ruv`, and for `enrich_values_ruv_covariates_note` and `enrich_adjust_ruv_help`, which say what
  the adjustment does without saying what it does not do.
- What the adjustment did deliver in the simulation: the error rate within the limit on every deciding class (at most
  5.4 %), and more power than the form without it (sets of 10: 45 to 48 % against 40 % where the unwanted factor is
  independent of the groups, 30 to 34 % against 17 % where it follows them).
- By ticket 6.5 (R7, and VAL-34) this does not hold the merge and the Methods prints the number. The wording of the
  sentences around the numbers is the owner's.
