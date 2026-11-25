# KIT-250: outcome

Ticket r1.3. The owner's decisions of 2026-11-25 (item 7, "take all recommendations"; methods ledger VAL-240): OD-1 = (A) as corrected by Re-check 1, OD-2 = before the release commit, O3 = the single-organism limit stated to the user (draft in `wording.md`).

Branch `kit-250` from `v0.9-staging` e0747fc1. The import, design and test files the ticket cites are byte-identical between 16112917 (where G0 ran) and e0747fc1.

| Commit | What |
|---|---|
| c101eeb9 | G0 records brought over from `kit-250-g0` (scripts and observations only) |
| d6bbcafb | tests only: `tests/testthat/test-kit250-organism-dialog.R` (E1a to E7, G1 to G3) and the declared edit `test-prot-01c-import-module-contracts.R:2099-2108`; nothing under `R/` |
| 7b674539 | the fix (5 files under `R/`) |
| da0593f7 | NEWS, ledger row `VAL-<merge>`, `wording.md` |
| f842a8ec | records: R1 and R2 tables, this file, `runs.tsv` |
| 8239aae1 | merge of `v0.9-staging` 63143aa3 (ledger: staging's rows kept, KIT-250's row after VAL-259; nothing else conflicted) |
| 43c1a304 | correction pass, tests only: E8; I02's assertions in E4 and E5; I03's in E3 |
| 77296fcf | correction pass, fix: I02 (W1's fixed id, removed when the next import starts) |
| cf908a7d | correction pass: `correction-predictions.md` (runs declared before sending) and `h1.R` |
| 6377cde6 | correction pass, test mechanics only (VAL-54 instrument rerun 1 of 3): E8's seams replace the caller's arguments of the same name |
| (records) | the correction pass's records, this file |

Core-hours: 299.8 job-minutes = 19.99 core-hours against the ticket's expected 24.0 and ceiling 102.0:
- first pass: G0 0.99, R1 0.57, R2 5.95, the refused request 0.01;
- correction pass: R2' 5.55, R2'' 6.04, H1 0.89.

The correction pass is the ticket's one correction (AUDIT-KIT-250-impl-20261125). R2' was an instrument defect (VAL-54); R2'' is its rerun 1 of 3.

Commit order (AUDIT I09): in the first pass the fix 7b674539 was committed 31 s after the tests-only commit d6bbcafb and 13 s before R1 was sent. R1 ran the tests-only commit (`package_commit` d6bbcafb), so the falsification stands. In the correction pass the tests-only commit 43c1a304 precedes the fix 77296fcf.

Conventions (AUDIT I10): the new comments and a few parameter lines exceed R.md's 100 characters, as the module's existing KIT comments do. The test file uses `function(...)` for its fakes, as the module's tests do.

# Gate G0 (base observation, run 45594289742)

- **Commit under test:** `v0.9-staging` 16112917329a63ab062655f04db5bf1dd7df9729 (its head on 2026-11-25).
- **Records commit:** `kit-250-g0` 804a69bd. It differs from 16112917 only in this folder (`environment.tsv`: `diff_only_records_folder TRUE`).
- **Run:** 45594289742, request `kit250-g0-804a69bd-2`, one scripts job, 14.8 job-minutes (0.99 core-hours). The ceiling was 60 job-minutes (4.0 core-hours).
- **Instrument:** `g0.R` sources KIT-233's driver, copied byte for byte into `k233-journey.R` (from `origin/kit-233-journey` cce7b04c, sha256 9b62d7ed…4965). Each session runs in test mode on the prot_dia fixture (one organism, Homo sapiens, 400 proteins), in a network namespace with loopback only (`environment-inner.tsv`). There were five fresh sessions (A, B, C, D1, D2) and all five attached.
- **Predictions:** written in `g0.R` before any value was read, and written to `g0-predictions.tsv` before the first session started.
- **Result:** 72 of 72 predictions as predicted; none differs and none went unobserved. Another 175 values were recorded without a prediction. Table: `g0-45594289742/g0-base-values.tsv`. Verdict: `g0-45594289742/verdict.txt`.

A request file `kit250-g0-804a69bd` was pushed as the blank template by mistake. Run 45594230854 refused it at `plan` and no job ran. It was not an attempt. It was replaced by `-2`, which was read back before it was pushed.

## Observed (base, 16112917)

Before the import, the form's Taxonomy ID was set to 10090 and its name to "Mus musculus". The option was then ticked. The FASTA's organism is 9606.

| Case | What was observed |
|---|---|
| (a) dialog open, before any answer | "Select Primary Organism" was shown, with 9606 preselected. The import was already `setup_import` = complete and the design step pending. The setup log held taxon 10090 / Mus musculus. The page behind the dialog said "[OK] Data imported successfully!". The form's taxon input was disabled. |
| (a) after Cancel by a real mouse event (chromote `Input.dispatchMouseEvent`; the topmost element at the point was the Cancel button) | The dialog closed. The statuses were unchanged (complete / pending, all later tabs disabled). The setup log kept 10090 / Mus musculus, with no `mixed_species_selection`. The page said "[OK] Data imported successfully!" and the Data Summary said "Organism: Mus musculus (taxon: 10090 )". **No notification at all** appeared after the Cancel. |
| (a2) after that Cancel (Re-check 1, R1-01) | The Design tab opened, `data_available` was TRUE and the builder was shown. A design was built and **saved** (design_matrix complete; QC opened). The design step logged **"Target organism taxon ID: 10090"**, so `workflowData$taxon_id` was the form's value. With the network closed, the UniProt warning said that annotation was not retrieved; it did not name the organism. |
| (b) dialog open, a real mouse event on the "Design Matrix" tab link | The tab link was **not in view** before the script scrolled the page. A user cannot scroll the page behind the dialog. After scrolling, the topmost element at the link's centre was the dialog's own `div.modal-footer`. The click did not reach the tab: the Design tab was not active and the dialog stayed open. |
| (c) Confirm | The dialog closed and the statuses were complete / pending. The setup log, `mixed_species_selection` and the form all held 9606 / Homo sapiens, with the filter FALSE. Notices: "Primary organism set to: Homo sapiens (Taxon: 9606)" and "Data import successful!". The Data Summary said "Organism: Homo sapiens (taxon: 9606 )". A third message, "Multi-species FASTA detected. Filtering to Homo sapiens enabled.", also appeared; by reading, it comes from the enrichment tab's auto-enable (`R/mod_omic_enrich.R:538-546`). |
| (d1) Cancel, then "Process Imported Data" again with the option ticked, then Confirm | The dialog was shown again. After Confirm the import was complete, the setup log held 9606, both notices appeared and there was no error. The second import on the attached session **succeeds**. |
| (d2) Cancel, then option unticked and import again | The import completed (its completion line was logged). The setup log held 10090 / Mus musculus and the page said "[OK] Data imported successfully!". The **"Processing Data" dialog stayed open** and "Data import successful!" was **not** shown (P01 confirmed). There was no error. |

Every tab status in every case is in `tab_status_all`. In all (d) states the tabs after the Design tab stayed disabled.

## Inferred (not executed)

- `workflowData$taxon_id` cannot be read from the browser on base. It was read directly only in (a2), from the design step's log line. In (a), (c) and (d) it is inferred from the setup log and the form, which the code writes in the same calls (`finalizeProtImportSetupState` :159/:176-188; Confirm :172-196).
- `workflowData$mixed_species_analysis` is inferred in the same way. Its proxy is the setup log's `mixed_species_selection`, which Confirm writes in the same call.
- (d2)'s stuck dialog is consistent with the waiting flag staying TRUE after Cancel. The flag itself was not read.

## Stop rules and acceptance criterion A1

- **Stop rule 1** (import not complete after Cancel, or a notice naming the organism): **not fired.** The reading stands: Cancel keeps the import complete, with the form's organism and no notice.
- **Stop rule 2** (a person can press a button behind the dialog): **not fired.** Limitation: at the point used, the dialog box itself covered the link, so this shows the link is unreachable there. It does not separately test a point covered only by the backdrop. Before scrolling, the link was also out of view.
- **Stop rule 4** (a second import refused or failing): **not fired.** Both re-imports complete. (d2) leaves "Processing Data" open, which is P01 as predicted, not a refusal. W1's advice is therefore usable.
- **A1:** met for G0. (a) to (d) are recorded with values and screenshots (`g0-45594289742/g0-*.png`, `g0-*.json`, `g0-base-values.tsv`).
- G0 replaces "by reading" in section 10's base column: every base value was observed as stated. Re-check 1's premise is confirmed on base: after Cancel a design can be saved with the form's organism (a2).

# Implementation

## The fix (7b674539)

| Ticket outcome | Where (fix commit) |
|---|---|
| 1. While the dialog is open the import is not complete: `complete()` is not called, `setup_import` is written `incomplete`, no status is written for the design tab | `R/mod_prot_import_processing_helpers.R:207-214` (`completeProtImportSuccessState`) |
| 1. The status panel shows W2 while the dialog is open and W1 after Cancel, never "[OK] Data imported successfully!"; the Data Summary's organism line reads "Organism: not confirmed" while either flag is set | `R/mod_prot_import_ui_helpers.R:1-24`, `:59-63`; wiring `R/mod_prot_import_server_helpers.R:166-179` |
| 2. Confirm completes the import (`complete = ...completeOmicImport(...)`, the default form of `completeProtImportSuccessState`), after the optional filter and before its two notices; clears `organism_unconfirmed` | `R/mod_prot_import_organism_helpers.R:155-156`, `:222-226` |
| 3. Cancel is the module's own `actionButton(ns("cancel_organism"), "Cancel")`; `cancelProtImportOrganismSelection` clears the waiting flag, sets `organism_unconfirmed`, removes the dialog, clears the import as the error path does and shows W1 as a warning that stays until closed (`duration = NULL`) | `R/mod_prot_import_organism_helpers.R:41-43`, `:244-259`; observer `R/mod_prot_import_server_helpers.R:447-452` |
| The error path's clearing, shared (no change in what it clears) | `clearProtImportWorkflowState`, `R/mod_prot_import_processing_helpers.R:266-282`, called by `resetProtImportWorkflowStateOnError` and by Cancel |
| 5. The flags are reset when an import starts and when an import fails | `R/mod_prot_import_orchestration_helpers.R:104-106`; `R/mod_prot_import_processing_helpers.R:255-257` |
| The new flag's initial value | `R/mod_prot_import_server_helpers.R:17` |

Unchanged (outcome 4): the form, the option's default, the radio's preselection, the values Confirm writes, and the "could not extract organism information" path (G3). `man/` and `NAMESPACE` are unchanged, because the new functions are internal (`@noRd`).

## R1: failing first (run 45602311219, tests-only commit d6bbcafb)

Declared before sending: E1a to E7 and the declared edit fail by assertion, G1 to G3 pass; ceiling 45 job-minutes.

Observed (`r1-45602311219/tests.tsv`, `failures.txt`):
- The new file: 31 passed, 32 failed, 0 errors. Every failure is in E1a to E7 and none in G1 to G3.
- `test-prot-01c`: 434 passed, 3 failures, all in "KIT-102 W6: completeProtImportSuccessState writes the statuses…" (the declared edit).

Each E test failed at the base value the ticket states:

| Test | Failed at |
|---|---|
| E1a | `complete` called ("complete:notify=FALSE"); setup_import "complete"; design "pending" |
| E1b | setup_import "complete"; complete called |
| E2 | complete not called (0) |
| E3 | no `import-cancel_organism`; `data-dismiss` and `data-bs-dismiss` present |
| E4 | `cancelProtImportOrganismSelection` absent |
| E5 | flags TRUE after the start and after the error path; events "complete:notify=FALSE" only (no removeModal, no notice) |
| E6 | "[OK] Data imported successfully!" in both states; "Mus musculus" and "10090" in the summary |
| E7 | function absent |

## R2: base against head (run 45603816096)

Base d6bbcafb, head da0593f7, three shards per side, plus R CMD check on both sides. Declared before sending: every file `same` except the new file and the declared edit; ceiling 600 job-minutes.

Observed (`r2-45603816096/compare.tsv`, `check-compare.txt`):

| File | Base | Head | Verdict |
|---|---|---|---|
| `test-kit250-organism-dialog.R` | 31 passed / 32 failed | 77 / 0 | fixed |
| `test-prot-01c-import-module-contracts.R` | 434 / 3 | 437 / 0 | fixed (the declared edit) |
| `test-prot-01d-import-module-characterization.R` | 79 / 0 | 79 / 0 | same |
| `test-module-ci-prot-import.R` | 107 / 0 | 107 / 0 | same |
| `test-module-ci-trx-app.R` | 50 / 0 | 50 / 0 | same |
| `test-workflow-gui-memory-lifetimes.R` | 81 / 0 (11 skipped) | 81 / 0 (11 skipped) | same |
| `test-workflow-server-shared.R` | 226 / 0 | 226 / 0 | same |

R CMD check: 0 errors, 3 warnings and 3 notes on both sides; the findings are identical on both sides (none only on head, none only on base). No kit-250 path appears in the portable-file-names note.

The confirm guards (`test-prot-01c` :1792-1895 and :1905-1934; `test-prot-01d` :626-645) pass on head, as predicted by reading.

## Correction pass (AUDIT-KIT-250-impl-20261125, I01 to I12)

What was done:
- **I01 (a), E8.** A `shiny::testServer` test of the module's observers with a real DIA owner (pattern of `test-workflow-gui-memory-lifetimes.R`). In one session it presses:
  - a ticked import;
  - the module's own Cancel (`input$cancel_organism`, `R/mod_prot_import_server_helpers.R:446-452`);
  - 'Process Imported Data' again;
  - Confirm;
  - a third ticked import, Cancel, then untick and import.

  `MockShinySession` carried it. Cancel's and Confirm's modal and notification seams go through `cancelOrganismSelectionFn` and `confirmOrganismSelectionFn`.
- **I02, the fix (77296fcf).**
  - W1 is shown under the fixed id `prot_import_organism_cancelled` (`R/mod_prot_import_organism_helpers.R:256-257`; constant in `R/mod_prot_import_ui_helpers.R:8`). A second Cancel therefore replaces it.
  - `startProtImportProcessing` removes W1 beside the flag reset (`R/mod_prot_import_orchestration_helpers.R:95`, `:109`), as KIT-175's notice is removed.
  - The seam's default is a no-op outside a Shiny session. A plain `shiny::removeNotification` default would break the three existing function-level tests that call `startProtImportProcessing` without a session: `test-prot-01c` :1469 and the F-11 test, and `test-prot-01d` :664. That is a deviation from the register's literal "`removeNotification = shiny::removeNotification`", made for this reason.
  - Nothing else in `R/` changed.
- **I02 tests.** E4 asserts the id. E5 asserts the removal at the start.
- **I03 tests.** E3 asserts on the whole window:
  - `data-backdrop="static"`, `data-bs-backdrop="static"`, `data-keyboard="false"` and `data-bs-keyboard="false"` (the installed Shiny 1.14.0 renders all four);
  - no close button in the header.
- **I12** (neither observer checks that the window is waiting): recorded under "Not checked"; no code change.

R2' (run 45611281361, base 43c1a304, head 77296fcf):
- **Instrument defect (VAL-54).** E8's `completeSuccessStateFn` fake added `messageFn`, which `finalizeProtImportProcessing` also passes. On both sides the import observer failed with "formal argument "messageFn" matched by multiple actual arguments", so E8 was not judged.
- **What was judged:**
  - E4 and E5 failed on base as predicted;
  - the six other files were `same`;
  - the check findings were identical.
- **Repair (6377cde6, test mechanics only):** the fakes now replace the caller's arguments (`utils::modifyList`). No expectation was changed.
- **Rerun as R2''.** Its base 34f3f606 is 43c1a304 with the same patch, on run branch `kit-250-r2c-base`.

R2'' (run 45614051035, base 34f3f606, head 6377cde6; `R/` as 77296fcf). Every prediction in `correction-predictions.md` held:

| File | Base | Head | Verdict |
|---|---|---|---|
| `test-kit250-organism-dialog.R` | 114 / 5 | 119 / 0 | fixed |
| `test-prot-01c-import-module-contracts.R` | 437 / 0 | 437 / 0 | same |
| `test-prot-01d-import-module-characterization.R` | 79 / 0 | 79 / 0 | same |
| `test-module-ci-prot-import.R` | 107 / 0 | 107 / 0 | same |
| `test-module-ci-trx-app.R` | 50 / 0 | 50 / 0 | same |
| `test-workflow-gui-memory-lifetimes.R` | 81 / 0 (11 skipped) | 81 / 0 (11 skipped) | same |
| `test-workflow-server-shared.R` | 226 / 0 | 226 / 0 | same |

Base failures in the new file: exactly the five declared, with 0 errors. They are E4 (the id), E5 (the removal at the start) and three in E8 (W1's id at Cancel; the removal at the second and fourth start).

E8 otherwise passes on both sides. This is observed in a Shiny session:
- Cancel through the module's observer clears the import (`data_tbl`, `config_list` NULL; setup incomplete);
- with the owner still attached, the re-import shows the window again;
- Confirm completes the import on 9606 with its table;
- Cancel, then untick and import, completes and removes the processing window.

E3's new assertions pass on both sides (guards). R CMD check: identical findings on both sides (0 errors, 3 warnings, 3 notes).

H1 (run 45611275368, I01 (b)): the corrected head 77296fcf in a headless browser on a runner. Records commit cf908a7d, which adds only `h1.R` and `correction-predictions.md` under this folder (`environment.tsv`: `diff_only_records_folder TRUE`). Test mode, network namespace. **61 of 61 predictions as predicted** (`h1-45611275368/g0-base-values.tsv`, `verdict.txt`).

Session A: the form was set to 10090 / Mus musculus, then the option ticked.
- **Window open:**
  - setup_import incomplete, design disabled;
  - the status panel shows W2;
  - the Data Summary says "Organism: not confirmed" and names neither Mus musculus nor 10090.
- **Esc** (`Input.dispatchKeyEvent`) and **a mouse click on the backdrop** at (10, 10), which landed on `div#shiny-modal`: the window stays open.
- **Cancel by a real mouse click** (the topmost element was the button's label):
  - the window closes, setup stays incomplete and the setup log is gone;
  - the status panel shows W1, with no success text, and the Data Summary is empty;
  - exactly one W1 notification (warning) appears.
- **The Design tab after Cancel:**
  - `data_available` FALSE;
  - Save Design not visible;
  - 'Import Existing Design' visible (it asks for its own organism);
  - the design status stays disabled.
- **'Process Imported Data' again:**
  - the window opens, W1 is gone (count 0) and there is no error;
  - so the re-import over a cleared bus with the owner attached works.
- **A second Cancel:** one W1, not two.
- **Again, then Confirm:**
  - complete; 9606 / Homo sapiens; both notices; W1 count 0;
  - "[OK] Data imported successfully!" and "Organism: Homo sapiens (taxon: 9606 )";
  - the builder's `data_available` TRUE.

Session D2: Cancel, untick, import.
- It completes, the "Processing Data" window closes and "Data import successful!" is shown.
- W1 is gone.
- The setup log holds 10090 / Mus musculus, which the user now entered in the enabled form.

To tell the owner, and notes (not changed here):
- **I04 (owner item D-1).** With 'Mixed species FASTA' ticked, four routes complete the import with the disabled form's organism and no window:
  - FASTA headers without organisms: the warning does not name the organism;
  - an error while counting organisms: the warning names no organism;
  - no usable protein column: silent;
  - FASTA processing failed: silent.

  This is pre-existing and outside VAL-240's Cancel scope. The audit recommends (a): refuse in the three silent or erroring cases, and name the organism in the warning.
- **I06 (owner item D-2).** Where the O3 limit is stated. The corrected draft and the open question are in `wording.md`.
- **I07.** A re-import over a completed analysis replaces the earlier table, FASTA table, configuration and setup log before the window appears. A Cancel then leaves the earlier analysis's design, later states and tab statuses over a cleared import, with nothing on screen saying so. The session owner names the cancelled import's input files. A failed re-import leaves the same state today. This is for a follow-up (coordinator).
- **I08.** What Cancel does not clear: `taxon_id`, `organism_name`, `fasta_metadata`, the mappings, `organism_distribution`, the owner's snapshot and input hashes, and the project files `aa_seq_tbl_final.RDS`, `fasta_metadata.RDS`, `config.ini` and the cache. None lets the builder save a design. 'Import Existing Design' asks for its own organism.
- **I11 (merge note for KIT-178).** KIT-178's count line "after 'Data import successful!'" must, for a mixed import, also sit in `confirmProtImportOrganismSelection`. Confirm receives no `dataImportResult`, so the counts must be carried there, with a test of the mixed case. KIT-233's J' answers the window with Confirm and is compatible; A4 in its `known.tsv` closes at the merge.

## Acceptance criteria

| # | Verdict | Evidence |
|---|---|---|
| A1 | met | G0 above; `g0-45594289742/` |
| A2 | met | R1: E1a to E7 and the declared edit fail by assertion, G1 to G3 pass. R2: all of them pass on head. Correction pass: R2'' shows the new assertions (E4, E5, E8) failing on the tests-only commit and passing on the fix |
| A3 | met | R2: every file `same` except the new file and the declared edit (both `fixed`); check findings identical. R2'': the same for the correction pass |
| A4 | met | NEWS under `### Fixes` (wording corrected per I05 (1)); ledger row `VAL-<merge>` after VAL-259 with all five cells (I05 (3)); this file with "For the change report" |
| I01 (b) | met | H1, run 45611275368: 61 of 61 predictions in a browser on the corrected head |

Stop rule 3 (a G row failing on head, or an unpredicted R2 verdict): not fired.

## Not checked

- I12: neither the Confirm nor the Cancel observer checks that the window is waiting. A Confirm click queued while the import observer still runs after the window was sent, followed by an error in that observer, would complete an import with no table. Very narrow and harmless, since the builder stays hidden; no code change by the coordinator's ruling.
- H1 covered DIA-NN in test mode only. It did not cover other routes, user mode with the file chooser, or a FASTA with more than one organism.
- I07's mixed state in practice (QC, normalisation or DA over a cleared import).
- `workflowData$taxon_id` outside (a2), and `mixed_species_analysis` on base: proxies only (see above).
- The UniProt cache (P04). It is not in this ticket and has its own ticket; in G0 the network was closed, so nothing was cached.
- Enrichment with a partially annotated mixed data set.
- A point covered only by the backdrop in G0 (b).
- User mode (G0 ran in test mode with uploads, not the file chooser).
- A FASTA with more than one organism (the fixture has one).
- Routes other than DIA-NN.
- The DIA session owner after a Cancel keeps the cancelled import's store and input hashes (I08). The re-import after it works (E8; H1).
- The full suite: none of its own; the merge train's.

## For the change report

**What was wrong, and who was affected.** This affected proteomics users on every route (DIA-NN, MaxQuant, FragPipe, Spectronaut, TMT) who ticked "Mixed species FASTA" and closed the "Select Primary Organism" window with Cancel. The option is off by default.

The import counted as finished before the window was answered. Cancel was a plain dismiss button that ran no code. As a result:
- The analysis kept the organism of the disabled form: Homo sapiens (9606) unless the user had typed another before ticking.
- The page said "[OK] Data imported successfully!" and showed that organism, with no message.
- A design could be saved under it, so the UniProt annotation queried that organism. In a mixed data set, where some accessions match, the user's own organism's proteins got no gene names and nothing said so.
- Enrichment's organism and the study parameters' Taxon ID were the form's.
- A second import in the same session with the option unticked left the "Processing Data" window open, so the user had to reload.

Quality control, normalisation and differential-analysis numbers were not affected.

**Code references.**
- Legacy, `main` 436cdbe8, `R/mod_prot_import.R`:
  - `:114`: the default 9606.
  - `:824`: the waiting flag.
  - `:868`: `modalButton("Cancel")`.
  - `:1059-1073`: completion before the answer.
- Before the fix, `v0.9-staging` e0747fc1 (identical to 16112917 for these files):
  - `R/mod_prot_import_organism_helpers.R:41`: `modalButton`.
  - `R/mod_prot_import_processing_helpers.R:206`: `complete()` before the answer.
  - `:218-223`: the waiting branch.
  - `R/mod_prot_import_organism_helpers.R:219`: the only place the waiting flag was reset.
- After the fix, `kit-250` 7b674539 and 77296fcf: see the table under "The fix" and the correction pass.
- **Since when:** first on `main` in 179eb0d2 (2025-12-04, "process multi-organism fasta files"); in the tags v0.3.6.3, v0.5.0-rc1 and v0.9.0 (`git log 436cdbe8 -S'waiting_for_organism_selection'`; `git tag --contains`).

**What changed in results.**
- With the option ticked and the 'Select Primary Organism' window shown, the Design Matrix Builder cannot save a design until an organism is confirmed: Cancel clears the import, and the page says "Import not finished: no organism was chosen. …" (draft wording). 'Import Existing Design' stays available and asks for its own organism.
- Not changed: with the option ticked but no window (FASTA headers without organisms, an error while counting organisms, FASTA processing that failed, no protein column) the import completes with the form's organism as before (AUDIT-KIT-250-impl I04). After a re-import over a completed analysis, Cancel leaves that analysis's later steps in place over a cleared import (I07).
- After Confirm, the organism, the annotation, the enrichment organism and the study parameters are the confirmed one, exactly as before.
- Imports with the option unticked are unchanged.
- A re-import after Cancel no longer leaves a window open.
- No numeric result changes for any import that was confirmed.
- Not repaired by this change:
  - A project that already annotated under a wrong organism keeps that annotation, because the UniProt cache is keyed by neither organism nor proteins (P04, its own ticket).
  - With a mixed-species FASTA, one organism is annotated and enriched (O3, a stated limit; `wording.md`).

**Evidence.**
- G0, run 45594289742: base behaviour in a headless browser on a runner, with a real mouse Cancel. The fixture was prot_dia (Homo sapiens, 400 proteins); the form was set to Mus musculus / 10090; test mode; network closed. All 72 predictions held. After Cancel the design saved with "Target organism taxon ID: 10090".
- R1, run 45602311219: the defect shown at function level on the tests-only commit.
- R2, run 45603816096: the fix, base against head, with R CMD check.
- Correction pass:
  - R2'', run 45614051035: E8 shows the module's Cancel observer and the re-import in a Shiny session; W1's id and its removal are shown;
  - H1, run 45611275368: the fix in a browser, 61 of 61 predictions.
- Limits: the fixture holds one organism; DIA-NN route only; test mode.

**Rollback** (AUDIT I05 (4)): revert the merge commit, which takes the fix and the tests-only commit together; reverting the fix alone leaves E1a to E7 and the declared edit failing. No stored format, record field or saved state changes (`organism_unconfirmed` is a session field).
