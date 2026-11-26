# KIT-278 outcome: test-omics-module-globals.R allows mod_omic_norm.R

Declared test edit (coordinator ruling 2026-11-26): KIT-261's `shown_idle <<- idle` (R/mod_omic_norm.R:282, the only `<<-` in the file, assigned in the enclosing `registerOmicNormRunButtonState` at :278) is a local closure of the kind the test allows. Branch `kit-278` from v0.9-staging 75acd421; head (the edit) 07e07e67; no product code changed.

## Run 45722796710 (file alone, base 75acd421 vs head 07e07e67, no check)

| File | Predicted base | Observed base | Predicted head | Observed head | Verdict |
|---|---|---|---|---|---|
| test-omics-module-globals.R | 13 passed, 1 failed (setdiff = "mod_omic_norm.R") | 13 / 1, the same failure | 14 / 0 | 14 / 0 | fixed |

Job minutes: about 17 (plan 0.4, head 7.0, base 9.1, results 0.5), under the 30-minute ceiling.

## Correction 1 (coordinator's ruling on audit K01): the line, not the file

The file entry (07e07e67) is replaced by the helper's own `allow` argument for exactly that line:
`omics_global_hits(sources, "<<-", allow = "shown_idle <<- idle")`, with `allowed_files` as on staging. A future `<<-`
elsewhere in `R/mod_omic_norm.R` still fails. Tests only; nothing under `R/`. Head 84941aba (predictions committed and
pushed first).

Run 45726680674 (file alone, base 75acd421 vs head 84941aba, no check): prediction held (`corr1-compare.tsv`,
`corr1-failures-base.txt`).

| File | Predicted base | Observed base | Predicted head | Observed head | Verdict |
|---|---|---|---|---|---|
| test-omics-module-globals.R | 13 passed, 1 failed (setdiff = "mod_omic_norm.R") | 13 / 1, the same failure | 14 / 0 | 14 / 0 | fixed |

Job minutes 17.83 (plan 0.38, head 6.65, base 10.31, results 0.49) = 1.19 core-hours, under the 40-minute ceiling. The
first run (45722796710, for the file entry) is superseded for the final head and kept as history.

## Not checked

R CMD check (test-only; rides the merge train). Other files. Windows.

## For the change report

Not applicable (test only; no product behaviour).

## Lesson

A ticket that edits R/mod_omic_*.R or the lipid/metab code the test scans must include tests/testthat/test-omics-module-globals.R in its scope.
