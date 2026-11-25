# KIT-207 outcome: the lipidomics, metabolomics and transcriptomics Summary refuses a table that is not the normalisation's result

Ticket `tickets/kit/KIT-207.md` r1.2 (owner decisions SD-1 to SD-6 taken as proposed, methods ledger VAL-148; pool X140,
X141). Branch `kit-207-summary-guard` from `v0.9-staging` 3f066385. No R was run on the laptop; every run below is a
`run-tests` run on GitHub-hosted runners (R 4.6.1, reference BLAS).

## Commits

| Commit | What |
|---|---|
| ce0cfbef | tests only: new `tests/testthat/test-kit207-summary-head-guard.R`; E10's change to `test-integration-input-lipid-metab.R` (declared exception) |
| 2d49258d | tests only, instrument fix: the lipid/metab rows install KIT-194's Shiny notification and progress doubles (`k194Notes()`) |
| 8f0c2eda | the fix: `R/utils_workflow_summary.R`, `R/mod_omic_norm.R`, `R/utils_workflow_trx_declarations.R`; `NEWS.md` |
| bb5c0c22 | tests only, instrument: G1 prints its masked study parameters |
| 5b9d1beb | records: this folder; the ledger row `VAL-158` |
| fe82f331 | tests only (audit I03): E10 asserts the `summary_none` text; E12 |
| (this commit) | audit repairs I01, I02, I04 to I06 (NEWS, ledger, records) |

## Step 0 survey (section 10)

Scope `grep` (section 10) at 3f066385: 23 files plus the five named and the new file (29). Read for a Summary action on a
`DIAStateAdapter` at a head outside the accepted set: only `test-integration-input-lipid-metab.R:221-231` (the declared
E10: a never-normalised metabolomics owner at `metab_raw_data_s4`). Every other lipid/metab/trx Summary call in the
scope is at a normalised head (`test-integration-input-*`: `_normalized`, `_ruv_corrected`, `trx_ruvs`; `test-kit172-downstream.R`:
after a journey; `test-workflow-omics-application-reads.R`: after Apply; `test-module-ci-trx-summary.R`: after DA; the
W1 and core captures: the module-CI owners, G8) or on proteomics, or on a manager that is not a session owner
(`test-workflow-figure-paths.R:881` NULL; `test-omics-summary-core.R:601` an environment). No other hit; R2 confirmed it.

## Runs

See `runs.tsv`. Core-hours used: 96.45 (R1 1.07, R1b 0.70, R2 15.15, R2-G1 1.85, R3 77.68) against 91.3 expected and a 449.4 ceiling; no correction pass was used. Observed:

- **R1** (45344983182, ce0cfbef): the trx rows E6 and E9, E10 and E11 failed at the base values stated; every lipid/metab row
  errored before the Save ("'session' is not a ShinySession object" in `withProgress` of the run): the new file did not
  install KIT-194's notification/progress doubles. Instrument defect (VAL-54): no Save was reached, no expectation judged.
- **R1b** (45347712087, 2d49258d, the fix of the mechanics only; no expected value changed): 155 passed, 101 failed,
  0 errors, 2 skipped (G7, G8: the guard does not exist on this commit). Every failure is a failing-first expectation of
  E1 to E9, E11 at its stated base value (`r1b-base-values.txt`), every G row passes:
  - E1 lipid/metab: "Retrieved DATA S4 object from state '<prefix>_qc_complete'"; `saveRDSFn` received
    `LipidomicsAssayData`/`MetaboliteAssayData` with LCMS_Pos maximum 1,373,458.887 (KIT-194 G0: 1,373,458.9); parameters
    and record written; `done` called.
  - E2 (head `<prefix>_log2`): retrieved `<prefix>_qc_complete` (the history fallback). E3 (head
    `<prefix>_intensity_filtered`): no state retrieved, `saveRDSFn` called with NULL, `done` called.
  - E4 (head `<prefix>_reset`): retrieved the earlier run's `<prefix>_correlation_filtered` (maximum 20.11).
  - E5 (loaded `<prefix>_qc_complete`): retrieved `<prefix>_qc_complete`, maximum 1,373,458.887.
  - E6 (trx head `trx_qc_complete`): "final transcriptomics state 'trx_qc_complete'"; `transcriptomics_normalized_data.tsv`
    written; "Study parameters saved successfully".
  - E7 (managed, held): one queued press, `state_name` `<prefix>_qc_complete`.
  - E8: the report rendered and the record bundle exported; `session_state_<date>.RDS` written (`feature_counts` 514).
  - E9: the TSV written; only "Study parameters saved with warnings"; `failed` with `integration_input_invalid`.
  - E10 (in R1): no error; `metabolomics_kit135_final_s4.RDS` written; the earlier input removed (that expectation passes
    on base, as stated).
  - E11 (the tab mounted under `testServer`): Save succeeded ("Study parameters saved successfully"), Copy copied
    ("Publication files copied with optional missing artifacts"), `study_parameters.txt` present.
  Stop rule 1 did not fire: no failing-first row passed on base.
- **R2** (45350537647, head 8f0c2eda vs base 2d49258d; section 10's 29 files, 3 shards per side, check): 27 files
  `same`; `test-kit207-summary-head-guard.R` `fixed` (base 155/101, head 264 passed, 0 failed, 0 skipped);
  `test-integration-input-lipid-metab.R` `fixed` (base 130/2, head 132/0: exactly E10's two failing-first expectations).
  No `new-failure`, no `message-changed`. Files failing identically on base and head (stored laptop references, KIT-156
  record): `test-integration-input-transcriptomics.R` 3, `test-omic-summary-module.R` 28, `test-omics-summary-core.R` 25,
  `test-module-ci-trx-summary.R` 3, `test-workflow-checkpoint-harness.R` 1 error, `test-workflow-kit-harness.R` 1,
  `test-gui-entry-transcriptomics.R` 13 (`r2-compare.tsv`, `r2-compare-failures.txt`). R CMD check base and head:
  0 ERROR, 3 WARNING, 3 NOTE; nothing only on head or only on base (`r2-check-compare.txt`). The GUI-entry lanes and the
  W1 captures are `same`; stop rules 2 and 3 did not fire.
- **R2-G1** (45359339101, head bb5c0c22 vs base a88c8fde = 2d49258d plus the same print; the new file only): an instrument
  run for G1 (below). The file is again base 155/101, head 264/0.

### G1 digests (A2)

| Omic | Run | object (sha256) base = head | assays base = head | record steps base = head | study parameters |
|---|---|---|---|---|---|
| lipidomics | R2 | 4a781cc6… = 4a781cc6… | 866bb980… = 866bb980… | equal | d2d2317e… ≠ 55016f56… |
| metabolomics | R2 | 925d4551… = 925d4551… | 8a270033… = 8a270033… | equal | b2a6d5b6… ≠ 3580bc58… |
| lipidomics | R2-G1 | 4a781cc6… = 4a781cc6… | equal | equal | d5712145… ≠ 070becb1… |
| metabolomics | R2-G1 | 925d4551… = 925d4551… | equal | equal | 1658a1f8… ≠ b0d5b21c… |

Record steps on both sides: `initial, <prefix>_raw_data_s4, _intensity_filtered, _qc_complete, _log2, _normalized,
_correlation_filtered`. The study-parameters digest is not reproducible between runs of one commit: at 2d49258d it was
85fd6a39… (R1b) and d2d2317e… (R2 base) for lipidomics. R2-G1 printed the masked text: of 208 lines, base and head differ
in one line per omic, line 93, `pairs [input] = <sha256>` of the correlation step (the entity id of its input pairs
product, which embeds the store's per-session incarnation), and in nothing else (`r2g1-params-diff.txt`). **Observed:** the default-path
save hands the RDS writer the same object (byte-identical serialisation digest), writes the same record steps, and writes
the same study parameters apart from that run-varying identity line.

### G8 (W1 premise)

Head only (R2-G1, `r2g1-g1-g8-lines.txt`): the guard returns TRUE on `module_ci_lipid_summary_workflow()` and
`module_ci_metab_summary_workflow()` with the default configuration (head `<prefix>_correlation_filtered`) and with RUV
skipped (head `<prefix>_normalized`); the head's record components are `identical()` to `_normalized`'s in all four.

## R3

**R3** (45362773892; head e158c6eb = branch bb5c0c22 merged with `v0.9-staging` 42dd1751 (branch `run-k207-r3`), base 42dd1751; the full suite, 20 shards per side, check off): 463 files; 461 `same`; `test-kit207-summary-head-guard.R` `only-on-head` (264 passed, 0 failed); `test-integration-input-lipid-metab.R` `count-changed` (base 131 passed, 0 failed with the old block; head 132 passed, 0 failed: E10's block, three expectations in place of two). No `new-failure`, no `message-changed`, no timeout or crash; 39 files fail identically on base and head (stored laptop references, KIT-156 record) (`r3-full-suite-compare.tsv`, `r3-full-suite-compare-failures.txt`). 40 jobs, 1,165.2 job-minutes, longest job 53.2 min. KIT-206 and KIT-202 had not merged, so R2b is not due and no other ticket's new file needed adding (E-6).

## Acceptance criteria

| # | Criterion | Met | Evidence |
|---|---|---|---|
| A1 | failing-first tests fail on the tests-only commit at section 10's values; guards pass on both; all pass on the fix | met | R1 (trx, E10, E11) and R1b (lipid/metab); R2 head 264/0 and 132/0 |
| A2 | regression scope base against head `same` except the new file and E10's block; check: nothing only on head; G1's digests equal | not met as declared for G1's parameter digest, explained: line 93 prints the correlation step's `pairs` input as an owner-scoped entity id (store incarnation drawn per session, `R/utils_workflow_saved_entities.R:79-82`, `:352-358`), so the digest differs between any two runs, also of one commit (R1b, R2 base, R2-G1 base); every other G1 comparison is equal (object, assays, record steps, 207 of 208 masked lines); all other A2 clauses met | R2: 27 `same`, 2 `fixed`, check identical; G1 object, assays and record steps equal; study parameters equal except one run-varying identity line (R2-G1) |
| A3 | full suite on the merge with staging: no difference other than the new file and E10's block | met | R3 |
| A4 | NEWS, ledger row; no version, `NAMESPACE`, `man/` change; outcome with "For the change report" | met | `NEWS.md` (### Fixes), ledger `VAL-158`; `git diff 3f066385 -- DESCRIPTION NAMESPACE man` empty |

E11 was built and executed (not INFERRED): R1b base values above; head passes in R2.

## Deviations and recorded divergences

- E6, E9 and G5 run RUVs with `ruv_mode = "manual"`, k 2 (KIT-194's harness), not "automatic": the rows need k ≥ 1,
  which manual k 2 guarantees; E9's second run uses k 1 (as KIT-194 T4) so that its RUVs state is a new entity.
- G3 runs for lipidomics only: KIT-194 runs RUV on the F172 corpus for lipidomics only; metabolomics' RUV on F172 is not
  established there.
- G1's study-parameters digest masks, beyond the 'Timestamp:' line, the version block (Repository, Branch, Commit,
  Package Version, Source: they name the commit, which differs between base and head by construction) and the temporary
  root. It still varied between runs; R2-G1 shows why (above).
- Two tests-only commits instead of one (2d49258d: an instrument fix after R1), and one tests-only print (bb5c0c22).
- `summaryFailed`'s removal of the earlier integration input catches any error (an error inside the Save's error handler would
  escape the observer) and shows it as a warning; it is not logged (conventions/R.md:118-119). Recorded, as the Report's
  catch-all.
- Recorded divergences from `conventions/`: the refusals are plain words, not "state the offending value" (as KIT-194);
  the Report call site catches any error from the guard (as KIT-194 A05); new lines follow the surrounding code's width
  (about 130 characters, as KIT-194's texts), not `R.md`'s 100.

## User-facing texts for the wording pass (DRAFTS, VAL-41)

As implemented in `omicNormGuardTexts()` (`R/mod_omic_norm.R`), shown after "The session could not be saved: " (Save),
"Report generation failed: " (Report) and "Export failed: " (Session State):

- `summary`: "The current table is not the result of the last normalisation: Quality Control was reverted, or the
  normalisation was reset, after it ran. Run the normalisation (and the correlation filter or its Skip) again, or, if you
  have only just pressed Reset, press Revert once in a Quality Control sub-tab; then press this button again."
- `summary_trx`: "The current table is not the result of the last normalisation: Quality Control was reverted, or the
  normalisation was reset, after it ran. Run the normalisation again, then press this button again."
- `summary_none`: "No normalisation has been run in this session, so there is no final table to save. Run the
  normalisation (and the correlation filter or its Skip), then press this button again."
- `summary_loaded`: "This session file was saved after Quality Control was reverted past the normalisation, so its table
  ("<state>") was not normalised. Import the data again, run quality control and the normalisation, then press this
  button again."

Also shown by the Save (SD-4, existing text of `trxIntegrationInputSource`): "The session could not be saved: the accepted
RUVIII export was made from an earlier RUV step: run the export again, then save again".

- After a refused save whose removal failed (existing text of `removeIntegrationInput`, now shown as a persistent
  warning): "the earlier integration input at <folder> could not be removed".

The ticket's note stands: the
owner may prefer "The workflow arguments were not saved: " for these omics.

## Implementation audit repairs (`AUDIT-KIT-207-impl-20261124.md`)

I01 and I02 (NEWS, ledger, this file), I04 (`r1b-base-values.txt` regenerated whole from R1b's artifact: 23 lines, E1
carrying `qc_complete_max`; R3's identical failures are 39), I05 (A2's G1 clause), I06 (divergence and text recorded)
applied as the register gives them. I03 (the coordinator's ruling: add both tests): commit fe82f331, tests only.
E10 now also asserts the `summary_none` text; E12 (lipid, metab) runs SD-5 (b)'s main case: a valid save, a Revert to
`<prefix>_qc_complete`, a refused save into the same integration folder. Run 45376347870 (head fe82f331, the two files, 1 job, 18.7 job-min, 1.25 core-hours, from the correction row): `test-kit207-summary-head-guard.R` 292 passed, 0 failed; `test-integration-input-lipid-metab.R` 132 passed, 0 failed. E12 lipid and metab: the valid save wrote the input; the refused save showed the `summary` text, wrote nothing and removed the earlier input; E10's `summary_none` text asserted. The transcriptomics `summary_none` branch, the trx Report and Session State refusals, an `_itsd_norm` head, the A01 head at the Summary and trx Reset stay by reading (the register's I03 (c)).

## Hand-offs and with the owner

- KIT-206 (shared `summaryFailed`): `integration_input_invalid` is now shown at every omic's Save, proteomics included;
  KIT-206's tests of proteomics Save notifications must expect "The session could not be saved: <reason>" wherever they
  reach that class, and whichever of the two merges second writes the union of both branches (E-2).
- Reproducibility candidate (plan section 13, coordinator): "The analysis record and study_parameters.txt name a product
  input (the lipid/metab correlation step's `pairs`) by its owner-scoped entity id, which differs between sessions of an
  identical analysis and reads like a content hash; record the product's content sha256 (`payload$sha256`) or label the
  value session-scoped; check that replay of a lipid/metab record resolves it."
- With the owner (K207-I07, pre-existing, batched): a refused integration-input write (SD-4's refusal and every other
  `abortIntegrationInput` reason of the lipid, metab and trx input writers) leaves the omic's earlier input in
  `integration/<omic>/` while the same Save writes the head's parameters and record. Candidate text: "A save whose
  integration input is refused (`integration_input_invalid` from the omic's own input writer) removes the omic's earlier
  input, as K061-3 and KIT-207 SD-5 (b) do."

## Observed and inferred

Observed (runs above): the refusals and their texts at every E row; no file, RDS, record, queued product, report render,
record bundle or session-state file at a refused head; the earlier integration input removed (E10); the mounted tab's
Save refusal followed by Copy's 'Required publication artifacts missing: Study Parameters' (E11); the G rows' saves; the
default path's object, record steps and parameters (G1); the module-CI owners (G8); no regression in the 29 files.

Inferred (by reading, not executed): that the managed refusal also holds on a real worker (stub runner only); that the
Report and Session State refusals behave the same for transcriptomics (E8 ran lipid and metab); that a removal failure
of the earlier input is shown as a warning (no test makes `unlink` fail).

## Not checked

- A real managed worker; a browser; Windows; real data.
- Transcriptomics' Report and Session State refusals (by reading: the same guard call).
- A failing removal of the earlier integration input.
- Why the correlation step's `pairs` input identity differs between runs of one commit (outside this ticket; it moves the
  study parameters' text between runs, not between base and head).
- The integration tab's behaviour when an omic's input folder is absent after a refusal (SD-5 (b)).
- KIT-206 (shared `summaryFailed`) has not merged; its classes are not in the handler (the second to merge adds them, E-2).

## For the change report

**What was wrong.** In lipidomics and metabolomics, after Quality Control was reverted past the normalisation (or the
normalisation was reset), the Summary tab's 'Save Workflow Arguments' wrote the workflow's "final" object
(`<omic>_<label>_final_s4.RDS`), `study_parameters.txt` and `analysis_record.json` from the quality-control table (not
normalised), from an empty object (an RDS holding NULL), or from an earlier run's table, while the differential analysis
shown and on disk came from the normalised data; 'Generate Report' exported the record bundle of that table and 'Export
Session State' recorded its counts and arguments. A session file named after a QC state was saved the same way. In
transcriptomics the same Save wrote parameters and a record of a QC head and wrote an earlier run's RUVIII matrix as
`transcriptomics_normalized_data.tsv`; after a second normalisation without a new RUVIII export it wrote the earlier
export and showed only 'Study parameters saved with warnings'. Affected: any lipidomics, metabolomics or
transcriptomics user who saved the Summary after such a Revert or Reset, or after loading a QC-named session file, and
transcriptomics users who saved after a second normalisation without a new RUVIII export (the integration file).

**Code references.** Base 3f066385: `R/utils_workflow_summary.R:94-145` (`saveWorkflowSummaryParameters`, no guard;
`workflowSummaryFinalState` at `:47-73` takes the head or the history fallback, `R/utils_workflow_dia_sessions.R:1532-1541`),
`:304-355` (report), `:396-411` (session state); `R/utils_workflow_trx_declarations.R:942-946` (the export written without
checking its source). Fix 8f0c2eda: `R/utils_workflow_summary.R:86-108` (`assertWorkflowSummaryHeadFinal`), `:129-148`
(`summaryFailed`: refusal, removal of the earlier input, `integration_input_invalid` shown; to `:155`), `:157` (Save), `:359-366`
(Report), `:455` (Session State); `R/mod_omic_norm.R:692-703` (the four texts); `R/utils_workflow_trx_declarations.R:944-948`.

**What changed in results.** No completed analysis changes. Summary saves, reports and session exports at a head that
is not the normalisation's result are refused before anything is written; on a refused Save the omic's earlier
integration input (`integration/<omic>/`) is removed; a transcriptomics integration file is not written from an RUVIII
export made from an earlier RUVs state; a refused integration input is now named to the user at every omic's Save (proteomics included; not tested for proteomics). Saves on the default path
are unchanged (G1: same object and record steps; parameters equal apart from a run-varying identity line).

**Evidence.** R1 45344983182 and R1b 45347712087 (base values), R2 45350537647 (29 files, check), R2-G1 45359339101
(G1, G8), R3 45362773892. Synthetic data only: KIT-172's F172 corpus (520 features, seed 172001) for lipid/metab, the
module-CI `trx_screen` corpus for transcriptomics; in-process and stub-runner managed paths.


## Merge record (coordinator)

Merged into `v0.9-staging` on 2026-11-24 from 6b96f6a0 after the implementation audit and its re-check 1 (register `AUDIT-KIT-207-impl-20261124`: mergeable; I01 to I06 closed; I07 with the owner). The ledger row is VAL-158; the owner's decisions it implements are VAL-148. The full suite R3 (run 45362773892, on the merge with staging 42dd1751) stands: since then staging gained only KIT-203's enrichment code, which shares no path with the Summary guard; no further run was made on this merge. Corrections to the records above, made at the merge: R3 had 39 files failing identically on both sides (not 40, as `runs.tsv`'s R3 row still says); the I03 run (45376347870) took 1.25 core-hours from the correction row, so the total is 97.70 of 449.4 core-hours. KIT-206 (proteomics export and Summary) merges after this ticket: it writes the combined `summaryFailed` handler and adds `test-kit207-summary-head-guard.R` and `test-integration-input-lipid-metab.R` to its merge check. For the owner's wording pass (VAL-41): the four refusal texts and the removal warning listed above; the NEWS sentence "A refused save writes nothing and removes…" describes the head refusal, not the transcriptomics refusal after a second run without a new export, which writes the parameters and removes nothing.
