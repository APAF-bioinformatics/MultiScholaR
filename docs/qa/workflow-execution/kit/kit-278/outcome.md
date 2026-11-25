# KIT-278 outcome: test-omics-module-globals.R allows mod_omic_norm.R

Declared test edit (coordinator ruling 2026-11-26): KIT-261's `shown_idle <<- idle` (R/mod_omic_norm.R:282, the only `<<-` in the file, assigned in the enclosing `registerOmicNormRunButtonState` at :276) is a local closure of the kind the test allows. Branch `kit-278` from v0.9-staging 75acd421; head (the edit) 07e07e67; no product code changed.

## Run 45722796710 (file alone, base 75acd421 vs head 07e07e67, no check)

| File | Predicted base | Observed base | Predicted head | Observed head | Verdict |
|---|---|---|---|---|---|
| test-omics-module-globals.R | 13 passed, 1 failed (setdiff = "mod_omic_norm.R") | 13 / 1, the same failure | 14 / 0 | 14 / 0 | fixed |

Job minutes: about 17 (plan 0.4, head 7.0, base 9.1, results 0.5), under the 30-minute ceiling.

## Not checked

R CMD check (test-only; rides the merge train). Other files. Windows.

## For the change report

Not applicable (test only; no product behaviour).

## Lesson

A ticket that edits R/mod_omic_*.R or the lipid/metab code the test scans must include tests/testthat/test-omics-module-globals.R in its scope.
