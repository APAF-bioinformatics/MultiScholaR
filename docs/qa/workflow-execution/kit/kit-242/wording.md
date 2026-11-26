# KIT-242 wording drafts (VAL-41: the owner sets the final wording)

- Between-sample normalisation, an assay with duplicate feature identifiers (a malformed assay; it stops):
  "Between-sample normalisation was not run: assay '<assay>' has <n> duplicated feature identifiers (<first three>).
  Run 'Resolve Duplicates' in the Quality Control tab, then run the normalisation again."
- Between-sample normalisation, a sample with no value in an assay (log line and notification):
  "Assay <assay> was normalised without samples <samples>, which have no value in <assay>; they are left empty in
  <assay> only." A sample with no value in any assay: "Samples <samples> have no value in any assay and were removed
  from the design."
- ITSD (log line): "Internal standards used per assay: <assay> <n>; ..." Not applied (log line and notification):
  "Internal-standard normalisation was not applied to assay <assay>: no internal standard matched the pattern
  '<pattern>' in column '<column>'." (or: "the import named no annotation column in which to look for internal
  standards"; "no internal standard was selected for it in the ITSD selection table"). A sample without signal:
  "Assay <assay>: samples <samples> have no internal-standard signal; their values in <assay> are set missing, in this
  assay only."
- Session summary: "- ITSD applied: A1 yes (5 standards); A2 no (no internal standard matched)".
- Proteomics normalisation text: "Warning: samples <samples> have no value and were removed from the design."

## Part B (legacy RUV-III-C and the declared design), drafts (VAL-41)

- An assay kept without RUV, no elbow (log line "Warning: ..." and notification): "RUV-III was not applied to assay <A>:
  no number of unwanted factors was supported by the data (the automatic search found no elbow in the canonical
  correlations at any control percentage). The assay keeps its values without RUV."
- An assay kept without RUV, no usable controls: "RUV-III was not applied to assay <A>: the automatic search found too
  few control features observed in every sample: <the search's refusal>. The assay keeps its values without RUV." plus,
  when a sample has one value: " Samples <s> have one value each in this assay, so no control is observed in every
  sample." **Wording defects to fix (text only; strings in R not changed):** (1) when every percentage had fewer than
  five controls selected, the text states the wrong cause ("observed in every sample"; the cause is too few controls
  selected); (2) with KIT-190's refusal embedded the text reads "... or run normalisation without RUV.. The assay keeps
  its values without RUV." (a double full stop, and advice that contradicts the action taken).
- An assay kept without RUV by `ruvIII_C_Varying` (no control feature among its features with values): "RUV-III was not
  applied to assay <A>: none of its features with values is a control feature. The assay keeps its values without RUV."
- Samples left out (log line and notification): "RUV-III for assay <A> was computed without samples <s>, which have no
  value there; they stay empty in that assay." Features left out (log line): "Assay <A>: <n> features with no value
  were left out of RUV-III and stay empty."
- Closing log line: "RUV-III correction applied" (every assay corrected); "RUV-III correction applied to <A1>; not
  applied to <A2> (see the warnings above)"; "RUV-III not applied: no assay could be corrected (see the warnings above);
  the data continue without RUV".
- Stops (`omics_ruv_failed`): "RUV-III was not run: the number of unwanted factors for assay '<A>' must be one whole
  number of 0 or more (given: <k>)."; "RUV-III was not run: assay '<A>' needs the feature id column '<col>', at least one
  feature and at least two samples of the design."; "RUV-III was not run: assay '<A>' has values in fewer than two
  samples or in no feature."; "RUV-III was not run: the control features for assay '<A>' <are positions outside the
  assay | are an unnamed logical vector of another length than the assay | must be feature ids, positions or a logical
  vector>."; "RUV-III-C could not correct assay <A> (an error, no corrected matrix, or no finite value; see the log).
  Choose another k or control set, or run without RUV."; "RUV-III was not run: the search for assay <A> failed: <error>".
- The Manual search with no control set for an assay it could not read: "RUV-III was not run: the control-feature search
  found no control features for assay <A>, because it could not read the assay. The assay has <n> duplicated feature
  identifiers (<first three>). Run 'Resolve Duplicates' in the Quality Control tab, then run the normalisation again."
  (the last two sentences only with duplicates). Known limit (audit I02): for a group left with fewer than two samples
  holding values the first sentence states a wrong cause.
- NEWS known-limit clauses (drafts): the RUV bullet's "Known limit, not yet resolved: when a group is left with fewer than
  two samples holding values in an assay ..., the RUV step stops, and its message does not yet give the right cause (in
  Automatic mode it names no assay)."; the normalisation bullet's "Known limit, not yet resolved: on the DIA route that
  normalises and applies RUV in one step, the removed sample is not yet named in the tab or recorded in the analysis
  record; only the smaller design shows it."

## VAL-373 (Q-C) texts, drafts

- An assay in which a group is left with fewer than two samples holding values (log line "Warning: ..." and
  notification; also the RUV summary): "RUV-III was not applied to assay <A>: the groups <groups> have fewer than two
  samples with values in this assay, so control features cannot be chosen. The assay keeps its values without RUV."
- Proteomics, the DIA normalisation step (log and warning notification): "Warning: samples <samples> have no value and
  were removed from the design." (the limpa route's text).
- NEWS: the two known-limit clauses are replaced by what the code now does (the RUV bullet's group sentence; the
  normalisation bullet's DIA-route sentence).
