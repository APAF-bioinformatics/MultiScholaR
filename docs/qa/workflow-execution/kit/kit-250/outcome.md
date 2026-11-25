# KIT-250: outcome

Ticket r1.3. The owner's decisions of 2026-11-25 (item 7, "take all recommendations"; methods ledger VAL-240): OD-1 = (A) as corrected by Re-check 1, OD-2 = before the release commit, O3 = the single-organism limit stated to the user (draft in `wording.md`).

Branch `kit-250` from `v0.9-staging` e0747fc1. The import, design and test files the ticket cites are byte-identical between 16112917 (where G0 ran) and e0747fc1.

| Commit | What |
|---|---|
| c101eeb9 | G0 records brought over from `kit-250-g0` (scripts and observations only) |
| d6bbcafb | tests only: `tests/testthat/test-kit250-organism-dialog.R` (E1a to E7, G1 to G3) and the declared edit `test-prot-01c-import-module-contracts.R:2099-2108`; nothing under `R/` |
| 7b674539 | the fix (5 files under `R/`) |
| da0593f7 | NEWS, ledger row `VAL-<merge>`, `wording.md` |
| (records) | R1 and R2 tables, this file, `runs.tsv` |

Core-hours: 112.7 job-minutes = 7.51 core-hours (G0 0.99, R1 0.57, R2 5.95, refused request 0.01), against the ticket's expected 24.0 and ceiling 102.0. No correction pass and no instrument rerun were used.

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

## Acceptance criteria

| # | Verdict | Evidence |
|---|---|---|
| A1 | met | G0 above; `g0-45594289742/` |
| A2 | met | R1: E1a to E7 and the declared edit fail by assertion, G1 to G3 pass. R2: all of them pass on head |
| A3 | met | R2: every file `same` except the new file and the declared edit (both `fixed`); check findings identical |
| A4 | met | NEWS under `### Fixes`; ledger row `VAL-<merge>` after VAL-233; this file with "For the change report" |

Stop rule 3 (a G row failing on head, or an unpredicted R2 verdict): not fired.

## Not checked

- The fix in a browser. G0 ran on base only, and the ticket declares no head browser run. The head is judged by the unit tests E1a to E7 and G1 to G3.
- `workflowData$taxon_id` outside (a2), and `mixed_species_analysis` on base: proxies only (see above).
- The UniProt cache (P04). It is not in this ticket and has its own ticket; in G0 the network was closed, so nothing was cached.
- Enrichment with a partially annotated mixed data set.
- A point covered only by the backdrop in G0 (b).
- User mode (G0 ran in test mode with uploads, not the file chooser).
- A FASTA with more than one organism (the fixture has one).
- Routes other than DIA-NN.
- Interaction with the DIA session owner after a Cancel: the owner keeps the first import's store, as after a failed import today. G0 (d1) showed a re-import works on base with the data still in place; on head the state after Cancel is the error path's state.
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
- After the fix, `kit-250` 7b674539: see the table under "The fix".

**What changed in results.**
- With the option ticked, nothing downstream can run until an organism is confirmed. After Cancel the import is cleared and the page says "Import not finished: no organism was chosen. …" (draft wording).
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
- Limits: the fixture holds one organism; DIA-NN route only; test mode; no browser run on head.
