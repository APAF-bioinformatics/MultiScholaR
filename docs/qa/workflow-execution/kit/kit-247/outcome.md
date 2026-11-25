# KIT-247 outcome: a design saved again after Quality Control or normalisation

**State: done and ready for the closure audit.** The second and LAST correction (the owner's VAL-277 for AUDIT I01 and I02, and VAL-288) holds every prediction in the final comparison (45698414166). Earlier: implementation audited (AUDIT-KIT-247-impl-20261126): NOT ready for a merge train. R2 stopped at stop rule 3 (a regression in `test-workflow-native-session.R`); the one correction pass (R1', R2') holds every prediction. Two code defects found by the audit, I01 (major) and I02 (minor), are OPEN (below): the correction pass is used, and they wait for the owner's word on an extension (parked row 47). This commit is records only.

- Ticket r1.4 (owner's decisions VAL-239: OD-1 (D) on the whole design table, OD-2 (a), OD-3 (a), OD-4 (a) including the
  top-level integration file). E-5 stays open for the coordinator.
- Branch `kit-247` from v0.9-staging **e0747fc1** (KIT-222 merged). Commits: `526ec2aa` tests only; `6bd2c194` fix;
  `0da49542` fix 2 (departure 3); `de744acf` departure 3's test (its tests-first pair: `1bbe5c40` on `kit-247-r1b`);
  `1d412824` E13's mechanics (its rerun: `5909d0e8` on `kit-247-r1c`); records. The fix was written before R1 was sent: the local commit
  e888e439 (22:23, rebuilt as 6bd2c194) predates R1's request (d5ec5163, 22:26:07), and 6bd2c194 and 0da49542 were
  pushed at 22:25:29 and 22:25:57, before that request and before G0-4's third rerun returned (request 3fd8cca4,
  22:24:58). R1 ran on the tests-only commit 526ec2aa, so its result is unaffected; G0-4 no longer held R1 back by the
  coordinator's ruling of 2026-11-25 (DISPATCH-STATE). (AUDIT I03.)
- G0 (the base observation, v0.9-staging 16112917): records on branch `kit-247-g0`,
  `docs/qa/workflow-execution/kit/kit-247/outcome.md` there. G0-1 78/78, G0-2 18/18, G0-3 130/140 (three observations,
  trx groups NOT OBSERVED), G0-4 an application observation (the browser normalisation fails at a QC figure, "Problem
  while converting geom to grob", `g04d-45604291093/` there; another ticket's, not touched here). Stop rule 2: not
  decided by observation; by reading nothing disables the tab or the button (coordinator's ruling). KIT-259's plan audit
  (K259-P17): the earlier browser failure was a shinytest2 `warn = 2` artefact with two real findings owned by KIT-259;
  after KIT-259's A1 no removed-rows warning is raised on the label-free fixture, so this lane passes that point at
  `warn = 2` (any other warning still fails it there, F4): KIT-259 removed this instance, and KIT-247 does not count it
  as its own fix.

## Runs (`runs.tsv`)

### R1, failing first (45604461406, 526ec2aa, 35.9 min): every E row fails at its G0 value; G2 passes

| Row | Fails on base at | G0 value it matches |
|---|---|---|
| E1, E2, E12, E16 (label-free) | history keeps the later states (tail `correlation_filtered`, 8 states); QC "complete", normalisation "complete"; no "saved again" notice; Reset returns `protein_replicate_filtered` | G0-1 main: history unchanged, statuses unchanged, Reset P4 |
| E4, E5, E6 (TMT, DIA (iq), DIA with limpa) | the same; Reset returns `protein_replicate_filtered` (TMT, DIA) and `peptide_count_filtered` (limpa) | G0-2 main (limpa L7); TMT and DIA (iq) by reading (no G0 lane) |
| E8, E10, E11, E9 (lipid, metab, trx) | history kept; statuses complete; Reset commits `lipid_reset`, `metab_reset`, `trx_reset` | G0-3 main |
| E15 (lipid run after the re-save) | QC states between the design state and the run | G0-3 (e) |
| E3, E14 (label-free run after the re-save) | history keeps the QC states before `normalized`; the `prot_norm_input_current` notice fires once | G0-1 (e) |
| E17, E18 (label-free), E17 (lipid) | head moves to the design state (new entity); no "Comparisons saved" / "nothing changed" text | G0-1 (g), (h); G0-3 (g) |
| E19 (Cancel) | head moves to the design state and the bus design is replaced (no question) | G0-1 (b) to (h): the save commits at once |
| E21 (OD-4 (a)) | 4 proteomics integration files left | by reading (no re-save cut on base) |
| E22 (R1-03) | error: `compareOmicsDesignResult` does not exist | by construction |
| E23 (refused Save) | `lipidomics_lipid_canonical_final_s4.RDS`, `metabolomics_metab_lc_final_s4.RDS`, `transcriptomics_normalized_data.tsv` left | G0-3 f1 (observation 2) |
| E13 | INSTRUMENT: one Revert press reached `lipid_normalized` only (Revert selects the state before the head, `R/mod_omic_qc.R:59-79`), not below the run's start; Reset then committed `lipid_reset` | mechanics fixed (`1d412824`: Revert pressed until the design state), rerun R1c |
| G2 | passes (as base) | — |

Note: the first expectation of the cut rows (a modal was shown) also passes on base (the processing modal); the rows
fail on their other expectations. Not changed after R1.

### R1b (45608591579, 1bbe5c40): departure 3's tests-first pair

`test-kit247-import-config.R` on 6bd2c194 + the test: D3 fails (the configuration is the imported one, "imported", not
"session"); its guard passes. Fix 2 (`0da49542`) is the change it tests.

### R1c (45612214974, 5909d0e8): E13's mechanics, VAL-54 instrument rerun 1 of 3

With Revert pressed until the design state, E13 fails at its base value (Reset re-commits `lipid_reset`; no "no longer
current" text). Every other row fails as in R1.

### R2 (45612382879, base 526ec2aa, head 3ec7bcf8; declared in `r2-declaration.md`, files in `r2-files.txt`)

As it happened (`r2/compare.tsv`): `test-kit247-design-resave.R` fixed (base 35 passed, 77 failed, 1 error; head 117
passed, 0 failed); `test-kit247-import-config.R` only on head (2/2); 63 files `same`; check findings identical (0
ERROR, 3 WARNING, 3 NOTE on both; nothing only on head). **One new failure: `test-workflow-native-session.R`**, test
"KIT-026: the head records inputs by role and external resources, ...": passes on base, on head errors with
"Error in `paste0(\"kit247_design_resave_\", session$ns(\"\"))`: attempt to apply non-function". Cause (by reading):
the test calls `runProtDesignImportConfirmationFlow` with `session = NULL` on an owner whose history holds
`raw_data_s4` and later states; the imported and bus designs are both NULL, so `compareOmicsDesignResult` returned
"changed" and the gate tried to ask with no session. Stop rule 3: stopped and reported. The coordinator ruled it a
result keeping ONE correction pass (R1', R2').

### Correction pass (the one allowed; `correction-declaration.md`, predictions committed before sending)

- Merged origin/v0.9-staging 42846c0e (b88752e8): NEWS and ledger appends kept both, staging's first; the ledger's
  duplicated context rows VAL-230 to VAL-233 removed; no other conflict.
- Tests first `e6871939` (`test-kit247-resave-no-session.R`, N1 to N3); fix `0ce87bc2`: no earlier design on the bus
  is "first" (nothing to compare); a changed design over later results with no session is refused
  (`omics_design_resave_unconfirmed`), nothing committed (departure 4, ticket r1.4 section 19). The KIT-026 test is not
  edited.
- R1' (45629316617, e6871939, 21 min): as predicted, N1 and N2 fail, N3 passes, KIT-026 errors as in R2.
- R2' (45629323979, base staging 42846c0e, head d3a95721, 69 files, check, 545 min): as predicted, 66 files `same`
  (`test-workflow-native-session.R` passes again), `test-kit247-design-resave.R` (117), `test-kit247-import-config.R`
  (2), `test-kit247-resave-no-session.R` (7) only on head and passing; check findings identical (0 / 3 / 3), nothing only
  on head.

## Acceptance criteria

| # | Criterion | Evidence |
|---|---|---|
| A1 | G0 before code | `kit-247-g0` records (G0-4: application observation; stop rule 2 by reading, coordinator's ruling) |
| A2 | tests-only commit fails at G0's values; passes on the fix | R1, R1b, R1c; R2 head |
| A3 | every file `same` except the new ones; check: nothing only on head | R2: not met (one new failure); **R2' met**: 66 files same, the three new files only on head and passing, check identical |
| A4 | merge train's full suite | the coordinator's train |
| A5 | NEWS, ledger, records; `man/`, `NAMESPACE` unchanged | diff: `NEWS.md`, `docs/methods-decision-ledger.md` (`VAL-<merge>`), this folder; no `man/`, `NAMESPACE`, `DESCRIPTION` change |

### G rows (AUDIT I05 (f))

G1 (first saves): the GUI-entry lanes and the design module files, R2' same. G2: the new file (passes on base and
head). G3 (numbers unchanged when the steps are redone): no row (AUDIT I07). G4: `test-kit225-prot-norm-reset-source.R`,
`test-prot-norm-input-current.R`, R2' same. G5: `test-omic-norm-module.R`, `test-trx-norm-core.R`,
`test-omics-norm-state-guards.R`, R2' same. G6: the GUI-entry lanes, R2' same with identical pre-existing failures.

## OPEN defects (AUDIT-KIT-247-impl-20261126; not repaired: the correction pass is used; parked row 47)

- **I01 (major).** The confirmed proteomics Import Existing Design runs its commit outside the import's error handler,
  so an import that fails after the user presses Save (a KIT-222 refusal, a column-mapping refusal) ends the session
  instead of showing "Error during import: ..."; "Importing design files..." is left after Cancel (proteomics and
  lipid/metab); E7 is feasible and not written; E20's closure is wrong for this entry. The register's repair (quoted,
  not applied): "(a) In `runProtDesignImportConfirmationFlow`, make `flow` run `runProtDesignImportConfirmationFlowNow(...)`
  inside the shell's handler: `tryCatch(<call>, error = function(e) { msg <- paste("Error during import:",
  conditionMessage(e)); logger::log_error(msg); shiny::showNotification(msg, type = "error", duration = 15);
  shiny::removeNotification(successNotificationId); invisible(NULL) })` (and `shiny::removeModal()` for any processing
  modal it opened). (b) In both import frames remove the `importing_design` notification before the question is shown
  (or on the gate's "changed" branch). Test (one row, E7): on `prot_lfq` after `correlation`, copy `source_dir`'s three
  files into a new folder, move one sample in its `design_matrix.tab`, press Import Existing Design with
  `import_dir_text` = that folder: the question is shown (assert the modal title "Save the changed design?"); Save → the
  cut, the statuses and the notice as E1; a second copy whose `design_matrix.tab` carries a KIT-222 collision → Save →
  "Error during import" shown, history and head unchanged, the test continues (the session is alive)."
- **I02 (minor).** Departure 1 reads the comparison's state side from the live bus; a gated save that fails after
  hydrating the bus makes the next save of the same design "unchanged": the previous design's later steps kept, with a
  "kept" message. The register's repair (quoted, not applied): "Compare against what the last *successful* design
  commit left, not the live bus: record, at the end of each successful commit (after `saveState` in
  `commitOmicsDesign` / `commitOmicsDesignImport`; after `buildProtDesignStateCheckpoint` in both proteomics flows),
  `manager$.context$kit247_design_compared <- list(design_matrix = <the committed design>, data = <the committed data
  reference>)`, and pass it as `current` (the bus only when it is absent). Alternative: in each gated commit closure,
  snapshot the four bus fields (`design_matrix`, `data_cln`, `contrasts_tbl`, `config_list`) before hydration and
  restore them when the commit returns FALSE or errors. Test: a function-level row: a gated commit made to fail after
  the bus write, then the same design saved again → `"changed"`." Reachability (INFERRED): lipid/metab/trx rare (a
  failure after the bus write); proteomics builder blocked first by U-01; proteomics import ended first by I01.

## Engineering departures (recorded, with reasons)

1. The comparison's state side is read from the bus as the last design commit left it (design table and data input):
   the trx design state holds technical replicates summed into units (`R/func_trx_design.R:126-149`), so its stored
   design is not comparable run by run. In normal use the bus and the design state are written by the same commit.
   The comparison is of content in a canonical form (every column as text, columns in name order, rows in `Run` order,
   `tech_rep_group` recomputed from `group` and `replicates`; the data by saved entity, else by the SHA-256 of its
   content; `omicsDesignCanonicalTable`, `omicsDesignSameData`, `:331-367`), not of the store's table identities named
   in item 1 (kit-247-g0's outcome says 'item 1 compares the store's identities as written'). (AUDIT I05 (c).) Its
   limit: I02 (OPEN).
2. A save with nothing after the design state (the first save, or before QC) is committed as before even when it is
   unchanged, so G2 stays identical to base.
3. Proteomics Import Existing Design read the project's `config.ini` (`<source_dir>/config.ini`, copied from the
   package default when missing; `loadProtDesignImportedConfigAndTables`, `R/mod_prot_design_import_helpers.R:181`,
   `:209`) into the bus before the gate; left alone, Cancel
   would have kept the imported configuration (QC and normalisation defaults, the formula), which later steps read and
   which could change a result. Fixed (`0da49542`): the session's configuration is kept while the question is open and
   the imported one applies only on Save. It changes behaviour beyond the audited plan; its own tests-first pair is R1b.
4. A changed design over later steps reached with no application session (a scripted call) is refused with
   `omics_design_resave_unconfirmed`; nothing is committed (`runOmicsDesignResaveGate`,
   `R/utils_workflow_omics_design.R:536-540`).
5. With no earlier design on the bus the comparison returns "first" and the save is committed as before, with no
   question and no cut, even when later steps exist (`:380`; coordinator's ruling (i)); in the application this needs a
   session file without a design. **This is a limit (AUDIT I11):** such a save leaves the later steps standing after the
   re-committed design state, the state the defect starts from; reachable only from a session file without a design
   matrix, or a script.

Also: the core commits do not take a `cut` argument; the frames call `cutOmicsDesignHistory` right after a confirmed
commit (same effect, fewer signatures changed). The refused Summary Save's removal now uses
`removeOmicIntegrationInput` (folder and top-level file; coordinator's ruling on VAL-239 (iv), ticket r1.4).

## Rule matters (reported plainly)

- I attempted to amend fix commit `6bd2c194` and force-push `kit-247`; the permission system denied it. A force push is
  forbidden by the brief whether or not it is denied. Nothing on origin was rewritten; the change went in as a new
  commit (`0da49542`).
- Local history rewriting of unpushed commits: before anything on `kit-247` was pushed, I reset the branch
  (`git reset --mixed e0747fc1`, working tree kept) to rebuild it in order. It discarded two local, unpushed commits,
  `61c21f28` (tests only, without E13 and E23) and `e888e439` (the fix without the refused-Save removal), and recommitted
  their content with the additions as `526ec2aa` and `6bd2c194`.

## Notes from the implementation audit (recorded; no change made)

- **I07.** The original symptom (the DA on the old groups; the record's design rows against its `design_builder` step)
  is not asserted after the fix, and G3 has no row; by construction both now follow the new design (INFERRED).
- **I08.** The confirmation was never observed in a browser and its content is asserted nowhere (E1's modal count also
  passes on base, the processing modal). One browser observation of the confirmation is owed after KIT-259 merges and
  before the release tag (label-free: normalise, move one sample, Save Design: the question's title, text, buttons; Save:
  QC pending, later steps disabled, the notice; a second re-save answered Cancel: nothing changes; Import Existing
  Design's Cancel once I01 is repaired), recorded in the ledger row's status.
- **I10.** Lipid/metab/trx with a DA run: the DA tab's head observer (`omicDaFollowHead`, `R/mod_omic_da.R:135`,
  `:145-157`) sets DA back to "pending" after the cut set it "disabled"; a DA press is refused at the design head; no
  number changes.
- **I13.** Departure 4 is not a breaking change for users (no frame exported; application callers pass the session);
  through `runProtDesignImportObserverShell` with `session = NULL` the refusal escapes before the configuration is
  restored, leaving the project `config.ini`'s configuration on the bus (scripted only).
- **I14.** The confirmed proteomics builder commit runs outside `withDIAStateScope`; an unchanged design with different
  data entities is digested in full; memory not measured.
- **I16.** E-5 stays open. After the merge with KIT-243, Report and Copy are refused after a CUT until 'Save Workflow
  Arguments' is pressed again (INFERRED), not in the kept path.

## The merge train's scope (AUDIT I12, declared before sending)

R2's 69 files plus `test-module-ci-trx-design.R`, `test-module-ci-trx-qc.R`, `test-module-ci-lipid-norm.R`,
`test-module-ci-metab-norm.R` (helper-mediated: never run in R2/R2'), `test-kit243-summary-current-outputs.R`,
`test-kit244-limpa-pass-through.R`, `test-kit256-tech-rep-selection-cleared.R` (plus any file the I01/I02 pass adds);
`check: true`; base the staging head at the train. Predictions: every file present on base `same` (identical failure
messages where it fails on base); the three kit247 files only on head and passing (117, 2, 7, plus the new rows);
check findings identical, nothing only on head.

## Not checked

- E7: not written; feasible (each lane's source_dir holds the exported files; test mode takes the folder as text): see
  AUDIT-KIT-247-impl I01. E20: holds for both builders and the lipid/metab import; on proteomics Import Existing Design
  the refusal fires after the question, after hydration (AUDIT I01).
- The DA and the record after the fix (I07); the confirmation in a browser (I08); the unchanged re-save on DIA, DIA with
  limpa, TMT, metabolomics and transcriptomics (only label-free "same" and lipid "contrast" are tested); both import
  gates (I01).
- The trx moved-sample groups (G0: not observed).
- G0-4's browser click (the browser normalisation fails first; another ticket's).
- TMT and DIA (iq) base values were not observed in G0 (R1 shows them).
- Real data; Windows.

## For the change report

- **What was wrong.** Saving the design again after Quality Control or normalisation set nothing aside. Pressing
  'Run Normalization' then normalised the NEW design's table with NO Quality Control on label-free, TMT, lipidomics,
  metabolomics and transcriptomics, while the Quality Control tab still said complete, proteomics said the run used the
  current quality-control result, and the analysis record listed the earlier design's QC steps (a replay would apply
  them). Pressing Reset first brought back the PREVIOUS design's table on every route (DIA with limpa through KIT-225),
  and the normalisation, correlation filter and differential analysis ran on the previous design's samples and groups
  with the new comparisons, silently when group names were kept. Who: any user who changes the design after
  normalising (a mislabelled sample, an outlier, technical replicates, or one more comparison, which needs a re-save).
  Direction (by reasoning, not measured): a moved sample analysed under its old label attenuates fold changes and
  inflates residual variance; RUV-III on old groups removes real signal; runs of one sample analysed as separate
  samples is pseudo-replication (anti-conservative).
- **Legacy.** On `main` for proteomics, lipidomics and metabolomics; new to 0.9 for transcriptomics and the DIA with
  limpa numbers path.
- **Code.** Base e0747fc1 / fix `6bd2c194`, `0da49542`, `0ce87bc2` (lines at 55c3c607: `R/utils_workflow_omics_design.R:331-568`,
  `R/mod_omic_design.R:239-293`, `:380-433`, `R/mod_prot_design_builder_helpers.R:53-199`,
  `R/mod_prot_design_import_helpers.R:336-446`, `:482`, `:501`, `R/mod_omic_qc.R:307-308`,
  `R/utils_workflow_dia_products.R:965-968`, `R/mod_omic_norm.R:493-533`, `R/utils_workflow_trx_norm.R:169-188`,
  `R/utils_workflow_summary.R:170`, `:186`, `R/func_omics_design_shared.R:796-797`; base lines as the ticket's 2.2 to
  2.8): `R/utils_workflow_omics_design.R` (`compareOmicsDesignResult`,
  `runOmicsDesignResaveGate`, `cutOmicsDesignHistory`, `completeOmicsDesignCut`, `removeOmicIntegrationInput`,
  `applyOmicsDesignUnchanged`); `R/mod_omic_design.R` (`runOmicsDesignCommit`, `runOmicsDesignImport`);
  `R/mod_prot_design_builder_helpers.R` (`runProtDesignBuilderObserverShell`, `runProtDesignBuilderSaveFlow`);
  `R/mod_prot_design_import_helpers.R`; `R/mod_omic_qc.R` (the QC observer reads the active history);
  `R/utils_workflow_dia_products.R` (the notice's condition); `R/mod_omic_norm.R`, `R/utils_workflow_trx_norm.R`
  (Reset only while current); `R/utils_workflow_summary.R` (refused Save removes the top-level file too);
  `R/func_omics_design_shared.R` (KIT-222's text).
- **What changes in results.** For affected users only: after a changed re-save the later steps are set aside and must
  be redone (the numbers are then the new design's, with its Quality Control); an unchanged or comparison-only re-save
  keeps everything. No stored reference moves (R2).
- **Evidence.** G0 runs 45594073040, 45594077082, 45594081675, 45598375255, 45594085728, 45597423706, 45600644799,
  45604291093 (base); R1 45604461406, R1b 45608591579, R1c 45612214974, R2 45612382879 (one regression); R1' 45629316617, R2' 45629323979. Fixtures: the GUI-entry e2e lanes
  (small synthetic data); no simulation; no real data.

## Second correction (the owner's VAL-277 and VAL-288; the last correction for KIT-247)

- **Scope.** AUDIT-KIT-247-impl-20261126 I01 (the confirmed proteomics Import Existing Design runs inside the import's
  error handling: an error after Save shows "Error during import: ..." and the bus is put back; the "Importing design
  files..." notice is closed before the question) and I02 (a design save is compared with what the last successful design
  commit left, not the live bus); I07 and G3 rows; VAL-288: (a) a formula-only change asks and on Save sets the
  differential analysis back to be run again (both builders and both imports); (b) on the proteomics import an organism,
  taxonomy id or FASTA-only change asks, names what changed, and on Save applies it and sets only enrichment back. The
  confirmation's wording stays a draft (wording.md). Commits: merges of staging abfc99ac (0f9196a5) and 12701142
  (32c0ec57), tests fbf729ef and 97fd133f, fixes 40948dd8 and 1d0364e0, test mechanics 740a712f, records.
- **Declarations before each run:** `second-correction-declaration.md`, `val288-declaration.md` (with the I07
  classification rule and E15's gap).
- **Runs.** R1'' 45666094374 (as predicted except I07 lipid); R2'' 45666101179 (recorded as it is: 73 same, check
  identical, two head failures); R1''' 45676491133 (instrument: I02's planted failure, set with
  `local_mocked_bindings(.env)` and `deferred_run()`, outlived its row and failed every later lipidomics lane's design
  save; it also explains I07 lipid in R1'' and the two head failures of R2''); R1'''b 45688025752 (VAL-54 instrument
  rerun 1 of 3, the planted failure scoped with `with_mocked_bindings`: everything as predicted; the VAL-288 rows fail
  first as their tests-first proof).
- **I07 on lipidomics, classified by the rule written before the run:** neither a product result nor a test-design
  error. After the confirmed re-save, Reset brings nothing back and says why; the run from the design state completes
  ("Normalization pipeline complete!"), the correlation filter, export, load and DA run, the DA is on the NEW group
  (WT) and the record's design_builder step is the design state's entity (`r1pppb/kit247-log.txt`). I07b (Quality
  Control re-applied) and E15b (a run after the cut commits or explains) hold.
- **Final comparison 45698414166** (base staging 12701142, head 166ef4f5, the 76 files of `r2pp-files.txt`, R CMD
  check, 615 job-minutes): every prediction held: 73 files `same`; `test-kit247-design-resave.R` (188 passed),
  `test-kit247-import-config.R` (2), `test-kit247-resave-no-session.R` (7) only on head and passing; check findings
  identical (0 / 3 / 3), nothing only on head (`final/`).
- **Budget.** The owner raised the ticket's ceiling from 148.67 to 190 core-hours (VAL-309). Used: 90.8 before the
  second correction, 48.9 for R1'' to R1'''b, 41.0 for the final comparison: **180.7 of 190**.
- **Not checked (added):** the VAL-288 (b) FASTA path through a FASTA file (only an `aa_seq_tbl_final.RDS` in the import
  folder is compared; the organism path is tested); a formula-only change on the import paths (by reading); the
  confirmation in a browser (I08, owed after KIT-259 and before the release tag).
