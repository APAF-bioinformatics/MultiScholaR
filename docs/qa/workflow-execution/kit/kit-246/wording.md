# KIT-246: wording record

One shipped text changed: `integrationRecordTexts()[["methods_scaling_feature_zscore"]]`
(`R/utils_integration_record.R`), the clause that ends the selection sentence of the integration record's Methods
paragraph (report "Methods (draft)", `study_parameters.txt`, Record tab) when the model used the default per-feature
scaling (INT-24).

## Decision

- VAL-251 (owner, 2026-11-25): the Model tab's help for per-feature scaling (`help_scaling`), last sentence option (a).
- VAL-293 (owner, 2026-11-26): the Methods sentence is option (a), "the same audited passage in full" as that help,
  without the effect size (option (c) not taken).

## Old (base 127011429 and earlier)

> each scaled to mean 0 and standard deviation 1 after the selection (which used the unscaled values); in simulation,
> with few samples and missing values this scaling improved the recovery of weak factors, and otherwise made no
> practical difference.

## New (fix c74f2a95)

> each scaled to mean 0 and standard deviation 1 after the selection (which used the unscaled values). In simulations
> with 20 samples in which about a fifth of the values were missing, mostly low ones, this recovered closely (a
> correlation of at least 0.9 with the true pattern) a pattern carried by a tenth of the features of the two smallest
> data sets, each changed by about as much as its measurement noise, in 26 rather than 12 of 30 simulated data sets.
> Without missing values the difference could not be told from none (29 rather than 26 of 30); with 60 samples the
> pattern was recovered closely in 30 of 30 either way.

The first sentence keeps the paragraph's own description of the step (it follows the selection clause after ", ", so
it starts in lower case). Everything from "In simulations" is the help's passage, character for character.

## How the two texts are kept the same

The evidence passage is written once, in `integrationScalingEvidenceText()` (`R/utils_integration_gui.R`, `@noRd`).
`help_scaling` is `paste(<its first sentence>, integrationScalingEvidenceText())`; `methods_scaling_feature_zscore`
is `paste(<the paragraph's first sentence>, integrationScalingEvidenceText())`. A change to the passage therefore
changes both. The test (E5, `tests/testthat/test-integration-record.R`, KIT-246 block) asserts that both texts end
with `integrationScalingEvidenceText()`.

`help_scaling` is unchanged in content: its text on the fix commit is identical to staging 127011429's (checked by
string comparison of the two `paste()` calls before the commit; R2's `test-integration-results.R` L7 and
`test-integration-gui.R` read it).

## Differences from the ticket's section 7 draft (a)

The ticket's draft (a) read "In simulations made for this software with 20 samples (two groups of 10), ...,
this scaling recovered ... of the two smallest of four data sets ...". The owner's decision (VAL-293) is the help's
passage in full (its last sentence decided in VAL-251; the sentence before it is KIT-239's text, VAL-255), which reads "In simulations with 20 samples in which ...", "this recovered" and "of the
two smallest data sets". The decided words are used; every item of the ticket's "What must survive" list holds (no
"weak"; criterion and threshold; condition; "26 rather than 12 of 30"; the size in numbers; the other two results
kept apart and neither called small; "simulations"). E2's and E5's term lists of option (a) are all present.

## Whole-table scans

The new text contains none of "every computer", "any computer", "on every machine", "same total variance", "unit
variance" (`test-integration-model.R`), no `[R3]`-style placeholder and not "has not been measured yet"
(`test-integration-removed-factors.R`).
