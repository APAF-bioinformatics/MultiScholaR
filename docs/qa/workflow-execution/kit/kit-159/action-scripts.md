# KIT-159 action scripts (section 7.4; VAL-41, VAL-26)

The accepted per-route scripts of the default-user driver, as `valActionScript(route)` builds them
(`harness/msrval/R/scripts.R`). This file is outside the harness pin; the test
`test-action-scripts.R` pins the step ids below against the code, so a change to either fails the
harness tests until both agree. The rules and reasons are section 7.4's; they are not repeated here.

## What every script does

- **Import**: uploads the declared files (structural), chooses the route's format (and on DIA the
  route choice `DIA` or `DIA_LIMPA`); every other import input stays at its rendered default
  (lipidomics and metabolomics `sanitize_names` stays TRUE). Custom metabolomics column mappings are
  declared structural import values of the data set. The proteomics organism modal, when the app
  shows it, is confirmed with its primary button ("Confirm Selection") at its rendered values;
  `mixed_species_fasta = TRUE` is a declared forced input of a mixed-species data set, never part of
  the script.
- **Design**: the design builder through `guiEntryDesignActions` (extracted from the checkout): the
  declared groups, technical-replicate sets where the route takes them, the contrasts (a second and
  later contrast added before Save), sample names matched to the builder's offered runs with
  `e2e_match_available_values` on every route (the match is recorded; on routes whose import may
  rename samples, such as `sanitize_names`, the builder's offered names are the ones used). The builder's formula is left
  at its rendered value when the UI renders one (the helper's `~ 0 + group` is dropped then).
- **QC**: each sub-tab of the route's spec tables in order, its first `-apply_` button at rendered
  defaults (declared forced inputs applied where the sub-tab renders them).
- **Normalisation**: rendered or server-set defaults, then the arm's overrides, then "Run
  normalisation".
- **Correlation filter** (where rendered): its tab, then "Apply Correlation Filter" at the rendered
  threshold.
- **DIA_LIMPA**: "Compute RUV-III Covariates" after the correlation filter.
- **DA**: the DA tab, rendered defaults, "Run"; a confirmation modal the run raises is answered with
  its one primary button.
- **Summary** (only with `hooks.summary`): label and description, "Save Workflow Arguments".
- Never pressed: revert, skip/bypass, export, load filtered session, report, publication copy,
  enrichment, Cancel.

## Declared states

Each step declares the head state it must reach (any of a listed set), a tab status, or
`diagnostic`. QC sub-tab states are those of the GUI-entry oracles (`tests/testdata/gui-entry`);
normalisation admits the RUV-off and RUV-corrected states of the route; import and DA declare their
tab `complete` (neither changes the head). G1(b) records the states actually reached per route.

## Step ids (pinned)

- `prot_dia`: import, design, qc_peptide_qvalue_filter, qc_peptide_protein_peptide_filter, qc_peptide_rollup, qc_peptide_intensity_filter, qc_peptide_sample_filter, qc_peptide_replicate_filter, qc_peptide_imputation, qc_protein_rollup, qc_protein_cleanup, qc_protein_intensity_filter, qc_protein_duplicate_removal, qc_protein_replicate_filter, norm, correlation, da
- `prot_dia_limpa`: import, design, qc_peptide_qvalue_filter, qc_peptide_rollup, qc_peptide_log2_matrix, qc_peptide_accession_cleanup, qc_protein_dpc_quant, qc_protein_normalisation, qc_protein_peptide_count, norm, correlation, ruv_covariates, da
- `prot_lfq`: import, design, qc_protein_cleanup, qc_protein_intensity_filter, qc_protein_duplicate_removal, qc_protein_replicate_filter, norm, correlation, da
- `prot_lfq_fragpipe`: import, design, qc_protein_cleanup, qc_protein_intensity_filter, qc_protein_duplicate_removal, qc_protein_replicate_filter, norm, correlation, da
- `prot_tmt`: import, design, qc_protein_cleanup, qc_protein_intensity_filter, qc_protein_duplicate_removal, qc_protein_replicate_filter, norm, correlation, da
- `trx_screen`: import, design, qc_filter, qc_finalize, norm, da
- `lipid_canonical`: import, design, qc_intensity, qc_duplicates, qc_itsd, qc_finalize, norm, correlation, da
- `metab_lc`: import, design, qc_intensity, qc_duplicates, qc_itsd, qc_finalize, norm, correlation, da
- `metab_gc`: import, design, qc_intensity, qc_duplicates, qc_itsd, qc_finalize, norm, correlation, da
- `metab_combined`: import, design, qc_intensity, qc_duplicates, qc_itsd, qc_finalize, norm, correlation, da

With `summary`, each script ends with `summary_save`.

## Open points for G1(b) and G12 (recorded before any lifecycle run)

1. The DA confirmation modal (proteomics q-value warning) is assumed to carry one `btn-primary`
   button; if it does not, the cell stops `val_spec_invalid` and G12 flags it.
2. Forced inputs are matched by id or by id suffix (`-<name>`), so a declared QC setting applies in
   every sub-tab that renders it, as the GUI-entry lanes apply theirs.
3. The display matrix of DIA_LIMPA and transcriptomics cells (`matrix_display.tsv`) is read from
   the state object's `args` (`adjusted` or `corrected`) where present; the product may hold it
   elsewhere (the DIA_LIMPA adjusted matrix is a result value, `R/mod_prot_norm_workflow_helpers.R`
   :290). T-TAB on `prot_dia_limpa` shows whether it is found.
4. The design builder's formula: dropped from the helper's sets when rendered; if the builder
   renders none, the helper's `~ 0 + group` is set as part of the structural design.

## G1(b) alternatives (r1.14; S0, then L1a; no consumer runs them)

`valActionScript(route, alternative = ...)`, pinned by `test-alternatives.R`:

- `export_load` (any route but `trx_screen`, which has no load button; S0 and L1a run it on
  `prot_dia`, `lipid_canonical`, `metab_lc`): after the last normalisation step, `export_session`
  presses the tab's export button (`normalization-export_filtered_session`; lipidomics and
  metabolomics `norm-export_session`), then `load_session` opens the DA tab and presses "Load
  Filtered Session" (`differential_expression-load_filtered_session`; `de-load_filtered_session`)
  with no file chosen, which loads the latest export; both declare the head unchanged; then `da`.
- `no_covariates` (`prot_dia_limpa` only): the accepted script without `ruv_covariates`.

The two-species route (`prot_dia_two_species`) is `prot_dia`'s accepted script on the two-species
generated input with `mixed_species_fasta = TRUE` declared forced; the organism modal is confirmed at
its rendered values, as in every script.
