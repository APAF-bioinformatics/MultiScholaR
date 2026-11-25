# KIT-247 outcome: gate G0 (observation of the base)

**State: G0 judged.** G0-1 to G0-3 as below; G0-4 ended in an application observation (the browser normalisation fails at a QC figure), stop rule 2 recorded as not decided by observation (coordinator's ruling). R1 and the fix follow (the fix was drafted before R1's verdict). The owner decided OD-1 to OD-4 on 2026-11-25 (VAL-239). E-5 (should the Summary's Save refuse when the DA's contrasts differ from the bus's?) stays open
for the coordinator. No package code was changed: `git diff --name-only 16112917 <each commit below>` lists only this
folder.

## Binding

- Package: v0.9-staging **16112917** plus this folder (coordinator's binding, ticket r1.2). The ticket had bound G0
  to staging after merge train 4; Re-check 1 of the plan audit allows 16112917 because KIT-222 changes no lane's
  design (every lane's design is merge-free), so nothing G0 observes depends on it.
- Commits on `kit-247-g0`: `56ad3187` (the scripts, predictions written before any run), `19caa137` (G0-4 instrument
  fix 1), `dbef1d85` (trx probe instrument fix). Runs and minutes: `runs.tsv`.
- Scripts: `g0-lanes.R` (GUI-entry driver, `tests/testthat/helper-gui-entry.R`, each omic's module server under
  `shiny::testServer`; each case a fresh lane run up to a named lane step, then the re-save and what follows) and
  `g0-browser.R` (one e2e lane in Chrome through shinytest2). Spies only observe and call the wrapped function.
- Tables per run: `g0x-<run id>/` (`*-predictions.tsv`, `*-observed.tsv`, `*-compare.tsv`: each prediction with its
  pattern, the observed value and `held`).

## Step 0 (by reading at 16112917)

- Every cited line re-read; shifts for the three files the r1.1 header omitted are in the ticket's r1.2 header.
- Lanes: each GUI-entry lane saves the design once (`helper-gui-entry.R:464`); `:498`, `:582`, `:738` are refused
  presses before a valid design. No lane saves over later steps or presses Reset after a re-save.
- KIT-222 is not on 16112917 (`omicsDesignReplicateKeyCollisionText` absent); its text and refusal position were read
  at kit-222 75501687 for the ticket (r1.1, P05) and are rechecked at the branch point after train 4.
- R2 file list at 16112917 (the ticket's grep over `tests/testthat`, helpers excluded): test-kit129-run-path,
  test-kit207-summary-head-guard, test-kit212-da-current, test-lipid-01-qc-filtering-helpers,
  test-lipid-12a-norm-module-characterization, test-lipid-12-norm-module-contracts,
  test-metab-03b-norm-module-helper-characterization, test-metab-03c-norm-module-seams,
  test-metab-03d-norm-module-characterization, test-module-ci-lipid-design, test-module-ci-lipid-qc,
  test-module-ci-metab-design, test-module-ci-metab-qc, test-omic-design-module, test-omic-gui-config,
  test-omic-norm-module, test-omic-norm-ruv-design, test-omic-qc-module, test-omics-design-core,
  test-omics-design-import-refusal, test-omics-norm-core, test-omics-norm-itsd-selection, test-omics-norm-state-guards,
  test-prot-04a-design-module-characterization, test-prot-04-design, test-prot-05b-norm-module-contracts,
  test-prot-05d-norm-module-characterization, test-prot-annotation-uniprot-empty, test-prot-design-core,
  test-prot-design-log2, test-prot-design-save-preflight, test-prot-norm-core, test-prot-qc-coordinator-direct-shared,
  test-trx-da-without-ruv, test-trx-norm-core, test-trx-norm-module, test-workflow-checkpoint-dia-state,
  test-workflow-native-session, test-workflow-neutral-names; plus test-kit225-prot-norm-reset-source,
  test-prot-norm-input-current, test-gui-entry-*.R and the new file. `test-kit222-replicate-keys.R` arrives with
  KIT-222 (absent at 16112917).

## Observed against predicted

| Group | Run | Predictions held / not held / not observed |
|---|---|---|
| G0-1 prot_lfq (main, b to h) | 45594073040 | 78 / 0 / 0 |
| G0-2 prot_dia_limpa (main, e) | 45594077082 | 18 / 0 / 0 |
| G0-3 lipid, metab, trx | 45594081675 | 130 / 10 / 0 (6 observations below; 4 trx keys with no value: probe defect) |
| G0-3 trx probe rerun (VAL-54, 1 of 1) | 45598375255 | no error; the four trx group keys still `<absent>`: not observed (rerun budget spent) |
| G0-4 browser | 45594085728 | instrument defect, nothing observed |
| G0-4 browser rerun 1 (VAL-54) | 45597423706 | no observation: after Run the head was `normalized` (status pending); Skip then had no effect for 300 s (head stayed `normalized`, no notification) |
| G0-4 browser rerun 3 (VAL-54, last; KIT-233's capture and signal, logs saved) | 45604291093 | **application observation**: the normalisation fails in the browser with an error notification (below) |
| G0-4 browser rerun 2 (VAL-54) | 45600644799 | no observation: the head reached `normalized`, but the run's completion notice ('Normalization completed') never appeared within 300 s |

### What the base does (OBSERVED, every route run)

- **The re-save** (every route): head = the design state; the active history and every later state's entity
  unchanged; no tab status changed except as noted in observation 1; only "saved successfully" texts.
- **Path 1, Run with no Reset and no QC** (label-free; lipid, metab, trx): the run reads the NEW design state
  (`resave:<design state>`), the moved sample in its new group, with Quality Control still shown complete; the
  previous design's QC states stay between the design state and the run in the history and in the analysis record
  (a replay would apply them). Label-free: the KIT-182 notice fires ("This normalisation uses your current
  quality-control result. ..."), false here; the tab-entry warning stays silent. DIA with limpa: refused with a
  technical message, "Error in normalization: invalid protein normalization QC product request" (no numbers).
- **Path 2, Reset then Run**: label-free: Reset returns the previous design's Protein Replicate Filter table (notice:
  "Normalisation reset: the current table is again the one the normalisation started from (Protein Replicate
  Filter)."); DIA with limpa: the previous design's Peptide Count Filter table; lipid, metab, trx: `<prefix>_reset`
  re-commits the recorded start. The run, the correlation filter (and on limpa the RUV covariates) and the DA then
  run with no message, and the DA model's design holds the moved sample in its OLD group (KO) on label-free, limpa,
  lipid and metab (trx: see the rerun).
- **The record contradicts itself** after path 2 (P08): its `design_builder` step is the new design's entity
  (`resave:...`), its design rows put the moved sample in the old group (label-free and limpa observed).
- **(b) a renamed group** (label-free, WT to Ctrl): normalisation and the correlation filter run, then the DA
  refuses: "Error in DA analysis: DA contrast 'groupKO-groupCtrl' references terms absent from the model matrix:
  groupCtrl".
- **(c) DA pressed directly after the re-save** (label-free): refused, "no filtered protein data is loaded" (with the
  DA tab's entry warning "Please complete the normalization ... before accessing differential abundance analysis.").
- **(d) Accession Cleanup re-applied, then a run** (label-free): the record lists, in order, the re-applied Accession
  Cleanup, then the previous design's Intensity, Duplicate and Replicate filters, then the new normalisation and the
  previous correlation filter.
- **(f) Summary Save after the re-save, no run**: label-free with no DA: not refused by any guard; it takes the
  previous design's `correlation_filtered` (moved sample in KO); see observation 5. Label-free with a DA before the
  re-save: refused as stale ("The differential analysis was run on an earlier table ...") and the omic's integration
  folder removed (R1-08's prediction held). Lipid, metab, trx: refused at the design head in both copies.
- **(g) comparison-only re-save**: history, head (the design state) and every later state kept; the bus holds both
  comparisons. **(h) unchanged re-save**: the same, and every table identity of the design state identical
  (label-free: `protein_quant_table`, `design_matrix`, `protein_id_table` all TRUE; lipid, metab, trx all TRUE).

### Observations that differ from the predictions or were not predicted (OBSERVED; cause by reading)

1. **DA status reset to "pending" by a re-save when a DA had been run** (lipid, metab, trx, f1; prediction "no
   status changes" NOT HELD). Cause: the omics DA tab's head observer (`R/mod_omic_da.R:135`) calls
   `omicDaFollowHead`, which sets the DA status "pending" when `workflowDaResultStatus` is "stale"
   (`R/mod_omic_da.R:145-157`); after the re-save the head is the design state, so the earlier DA is stale.
2. **A refused Save removes `<integration>/<omic>` but leaves a top-level file** (lipid, metab, trx, f1; prediction
   "integration input removed" NOT HELD as worded): `removeIntegrationInput` removes the omic's folder only
   (`R/utils_workflow_integration_input.R:592-597`, called at `R/utils_workflow_summary.R:169-170`); left behind:
   `lipidomics_lipid_canonical_final_s4.RDS`, `metabolomics_metab_lc_final_s4.RDS`,
   `transcriptomics_normalized_data.tsv`. For OD-4 (a).
3. **The refusal text is misleading** (lipid, metab, trx, f0 and f1): "The current table is not the result of the
   last normalisation: Quality Control was reverted, or the normalisation was reset, after it ran." after a re-saved
   design (2.8 predicted this).
4. **trx's design rows hold units, not runs**: trx sums technical replicates into one sample per biological unit
   named by `tech_rep_group` (`R/func_trx_design.R:126-149`, `design$Run <- units` at `:137`), so the moved runs
   `Sera_03_r1`/`_r2` are absent from every trx table's design (the record's `<absent>`); the unit is `Sera_03`
   before the move and `RPMI_04` after. The first probe looked the runs up only, and could not read a trx design
   alone ("counts must be numeric": rebuilding a trx object for its design fails validation): instrument defect,
   rerun once (ruling).
5. **Label-free (f0): Save not refused, but "Study parameters saved with warnings"** and no integration input
   written; the integration export had started on `correlation_filtered` (moved sample in KO). The error behind the
   warning was not captured (the script records notifications, not the Summary's log); not traced.
6. **G0-4**: see below.

### G0-4, the browser

- First run (45594085728, instrument defect): the lane timed out waiting for the normalisation status to become
  complete after `e2e_run_prot_normalization_export` (`tests/testthat/helper-e2e-browser.R:1411-1428`), which presses
  Run, then Skip, then Export with no wait between them. By reading the application does mark it complete after a
  completed Skip (`completeProtNormCorrelationWorkflow`, `R/mod_prot_norm_correlation_helpers.R:396` →
  `finalizeProtNormCorrelationWorkflowState`, `:80-82`), so the likely cause is a Skip pressed before the run ended
  (INFERRED; no application log was captured). These e2e lanes had only been skipped on runners before.
- Rerun 1 of 3 (VAL-54; 19caa137): each press waited for the head state (`normalized`/`ruv_corrected`, then
  `correlation_filtered`) and the export for its notice; the status recorded after the run and after the Skip.
- Rerun 2 of 3 (633feb90): Skip waited for the run's completion notice; that notice never appeared in 300 s although
  the head was `normalized`.
- What is known (OBSERVED): in the browser on a runner, the label-free normalisation (RUV skipped) commits
  `normalized`, but neither a completion notice nor a working Skip follows within 300 s; the stepper's normalisation
  status stays pending. In process (G0-1, testServer) the same lane completes and shows "Normalization completed (RUV
  skipped)!". The browser's application log was not captured, so whether the normalisation never finishes in the
  managed browser session (an application observation, as KIT-233 A5 saw the state committed about 80 s before the
  step ended) or the driver misses its end (an instrument defect) was not established by rerun 2.
- **Rerun 3 of 3 (eedb59d0; coordinator's ruling) settles it: an APPLICATION OBSERVATION, not an instrument defect.**
  The driver kept every notification shown (as KIT-233's journey driver does, `journey.R:671-687`) and saved the
  application's log (`g04d-45604291093/applog-excerpt.txt`, `g0-g04-notes.tsv`). On entering Normalisation the
  pre-normalisation QC fails ("Error generating pre-normalization QC: Problem while converting geom to grob."), and
  'Run Normalization' then fails at the post-normalisation Pearson correlation plot with the error notification
  "Error in normalization: Problem while converting geom to grob." (log: "Error in normalization workflow: Problem
  while converting geom to grob.", after "[pearsonCorForSamplePairs] Final pair count: 20" on 400 proteins x 10
  samples). `normalized` is committed before the figures, so the head moves while the step fails; no completion notice
  follows, so Skip and the stepper never advance. OBSERVED on v0.9-staging 16112917, MaxQuant label-free e2e fixture,
  in the application launched by shinytest2 (`run_app(test_mode = TRUE)`, in process). The same lane under testServer
  (G0-1) completes, its figures not being drawn there. KIT-233's run L (45577106977, packaged entry, plain process
  launcher, DIA-NN 8-run fixture) completed this step (S6.1, 111 s, "Normalisation ... completed"): a different route,
  fixture and launch, so it neither confirms nor contradicts this. Not traced: which geom fails and why only here.
  For the coordinator to file (not this ticket's).
- Driver comparison (ruling): KIT-233 waits for the step's notification from a kept list with a long hard wait and
  launches the packaged entry with the plain process launcher; this driver kept AppDriver's `run_app` launch and took
  KIT-233's notification capture and completion signal.

## Stop rules

- **1** (both readings wrong for a route): not triggered. Path 1 read the new design's table with no QC on every
  route that runs it (label-free, lipid, metab, trx; limpa refuses); path 2 returned the previous design's table on
  every route, and the next run's source carried the old group where observed (label-free, limpa, lipid, metab;
  trx: the groups were not observed; Reset returned `trx_reset`, the run read it, path 1 read the new `trx_raw_data_s4`).
- **2**: not decided by observation; by reading nothing disables the tab or the button (coordinator's ruling). The browser lane could not complete a normalisation (the observation above); the lanes have shown the re-save path programmatically on every omic.
- **5** (an unchanged re-save gives different table identities): not triggered on any omic; item 1 compares the
  store's identities as written.

## Not checked

- TMT and DIA (iq) (no G0 lane; the ticket reads them as label-free and as refused on class).
- Import Existing Design as the re-save (E7).
- The cause of label-free (f0)'s "saved with warnings".
- G0-4 (the browser click and the stepper after a re-save); trx's moved-sample groups (G0-3).
- Real data; Windows; the integration tab's reading of a stale input.
- The DA outcome for a batch- or factor-only change (E22; by reading).
