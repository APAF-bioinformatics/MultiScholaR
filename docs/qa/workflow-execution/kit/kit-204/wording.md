# KIT-204 wording for the owner (VAL-41): drafts, the owner sets the final text

1. Summary line (KIT-172's wording, kept for parity, S4): `Not tested (no value in a compared group): <n>`.
2. Count-table status: `Not tested`.
3. NEWS entry (`### Fixes`), the implementation audit's draft (AUDIT-KIT-204-impl-20261125, I05), as committed:

    * Proteomics: proteins that could not be tested in a comparison (no value in a compared group, or none in the
      analysed samples) are now counted as 'Not tested'. The differential-abundance summary always shows the line
      'Not tested (no value in a compared group): n' (0 when every protein was tested); before, such proteins were
      counted only in its total. In the application and the `differentialAbundanceAnalysis()` method, the table of
      significant counts (`*_num_sig_da_molecules.tab` and `.xlsx`) lists them as 'Not tested' instead of 'Not
      significant' when there are any, and the bar chart of every status that `outputDaAnalysisResults()` writes shows
      them as their own bar. This change alters no other result. Analyses saved by an earlier version keep their table of
      significant counts until the analysis is run again (KIT-204).

   Optional sentence (the audit's I11, for scripted users), not committed: "The exported wrappers
   `protein_deAnalysisWrapperFunction()` and `daAnalysisWrapperFunction()` still count them as 'Not significant'."
4. Open with the owner's wording pass: "Total genes" counts proteins (F6).
