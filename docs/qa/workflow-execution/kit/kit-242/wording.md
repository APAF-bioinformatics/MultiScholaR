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
