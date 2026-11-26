# KIT-221 step 0: the ticket's cited lines on the merged base

Base: `origin/v0.9-staging` 176f9d3ff (ledger VAL-329), branch `kit-221` in worktree
`/home/doktersmol/Projects/APAF/MultiScholaR-kit-221`. Read only; no R run. The ticket (r1.1) cites lines at
0cf0e440 / 9a0bf9f0 / d048898d. The register's Re-check 1 (R01 to R04) governs where the ticket is unclear; the
owner has answered OD-1 = (A) (VAL-193 (1)); OD-4 moved to KIT-230 (parked, not this ticket).

## What moved (ticket line -> base line)

`R/func_prot_da_model_methods.R`
- `contrast_strings_to_use <- contrasts_tbl[, 1]`: :295 -> :295 (unchanged).
- the three routes: `convertDpcDAToStandardFormat` :357-364 -> :362 (limpa, whole `dpcDE` block :326-373 holds);
  `runProtRuvModelContrasts` :397-401 -> :398; `runTestsContrasts` :403-412 -> :403.
- KIT-201 added `:415-416` (the `da_model` TREAT-moderation record; ticket's "three lines after :413").
- first-table pattern: `da_results_table <- contrasts_results$results[[1]]` :427 -> :430;
  `contrast_name <- names(...)[1]` :430 -> :433; data-frame check :433-435 -> :436; row-id and comparison
  columns :439-455 -> :442-444; group-name mapping :444-455 -> :447-458; `return_list$contrasts_results` :457 -> :460.
- `nested_list` for `getSignificantData` :474-477 -> :477-480; the count call :505-509 -> :508-514, now
  with `q_val_thresh` (KIT-220) and `p_value_column` (KIT-204); the second count table (`num_sig_da_molecules`
  from `significant_rows`) now has the `Not tested` status first (KIT-204) and the two barplot filters leave out
  `Not significant` and `Not tested`. Wide, long, volcano list, p-value histogram: same code, shifted by +3.
- `exists("group_mapping")` stays :447 (F4, out of scope).

`R/func_prot_da_compat.R`
- `contrast_strings = contrasts_tbl[, 1][[1]]` :115 -> :115 (unchanged); `contrasts_tbl` read from the object
  at :32 (`checkParamsObjectFunctionSimplify`), so the new check goes after it; `getSignificantData` :132-133
  and `printCountDaGenesTable` :162-163 -> :133 and :163-167 (KIT-220 added `q_val_thresh`, as SD-3 said).
  The alias at :601-603. NAMESPACE :110 and :398 (functions), :533-534 (methods): unchanged, no diff since 0cf0e440.

`R/func_prot_da_model_stats.R` changed (96 lines: KIT-200, KIT-201, KIT-249, base-fit helpers). The facts the
ticket rests on still hold: `runTestsContrasts` makes one table per element of `contrast_strings` with
`purrr::map(contrast_strings, ... topTreat(t.fit, coef = contrast) / topTable(...))` (:634-640, :759-765) and
`names(result_tables) <- contrast_strings` (:883, was :823); `makeContrasts(contrasts = contrast_strings)` :523 ->
:487/:569-ish; q-values are per comparison (KIT-249: the Benjamini-Hochberg fallback now assigns); the
moderation (`contrasts.fit`, `eBayes`, `treat`) is one fit for all columns (`applyProtDaContrastModel` :291).
`R/utils_workflow_dia_da.R` changed (402 lines) but the application's kernel still builds a one-row
`data.frame(contrasts = parameters$contrast, stringsAsFactors = FALSE)` (:535) and calls
`differentialAbundanceAnalysis()` once (:543); the only call of the method in `R/` (grep). Workbooks: the
ticket's addendum says all eight are archived by KIT-234 (merged, train 4); the alias has no caller in `R/`
except its own wrapper at :603.

## Does the design apply unchanged?

Yes. Specifically for KIT-204: the count status `Not tested` is made inside `countStatDaGenes` for each table of
the list it is given (`R/func_prot_da_results_summary.R`, `countStatDaGenesHelper` maps `countStatDaGenes` over the
tables with `purrr::imap`, then `printCountDaGenesTable` row-binds), and the second count table groups
`significant_rows` by `comparison, status`; so passing the list of all comparison tables (the sibling's shape)
gives each comparison its own `Not tested` count; no cited behaviour changes so that an acceptance criterion
cannot be met. KIT-249 (finite q-value on every tested row), KIT-200 (untested rows appended to every table)
and KIT-220 (threshold passed) act per table. KIT-218 is merged (its file `test-kit218-da-wrapper-tail.R`,
guard G4, is in the base; the wrapper already carries the data-frame check and class `prot_da_contrasts_invalid`).

## Open for other tickets

- Open ticket branches touching the two files (`git diff --name-only origin/v0.9-staging...<branch>` over every
  local and remote branch): none among remote branches; the local branches `merge-rehearsal`, `run-k243-corr-base`,
  `sim-k124-06` and `sim-k124-confirm` differ from staging in these files but are not ticket branches for an open
  edit (a rehearsal and run/simulation bases). No other open ticket edits them.
- Register R01 (stop rule 6) is adopted: if W7 passes on base (the rows labelled `C_vs_A` equal the one-row call) or
  W8 completes on base, `purrr::map` passed labels and the factor premise is false; recorded, not a stop for the fix.

## Test fixtures (ticket section 4)

- limma and wrapper fixture: copied from `tests/testthat/test-kit218-da-wrapper-tail.R:10-52, :99-` (builders) which
  copy `test-kit200-prot-da-no-value-proteins.R`.
- limpa: fake `dpcDE`, `tests/testthat/test-prot-07d-da-model-stats-shared.R:670-697` (the object at :672-689).
- RUV factors in the model: `tests/testthat/test-prot-ruv-factors-in-model.R:8-52` (`ruvModelFixture`, `ruvModelObject`,
  `ruvModelPackageState`) and :54-62 (`ruvModelDa`).
