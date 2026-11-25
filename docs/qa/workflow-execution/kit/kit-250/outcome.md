# KIT-250 gate G0: outcome (base observation only)

Ticket r1.2. G0 only: no package code, no tests commit. The tests-only commit (R1) and the fix wait for the owner's OD-1 and OD-2 (both OPEN; O3 OPEN).

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
- A2 to A4 are not started: they wait for OD-1 and OD-2.

## For the next step (when the owner has answered)

- G0 replaces "by reading" in section 10's base column: every base value was observed as stated.
- Re-check 1's premise is confirmed on base: after Cancel a design can be saved with the form's organism (a2).
- The (a) and (a2) rows give R1's failing-first values.

## Not checked

- workflowData$taxon_id outside (a2), and mixed_species_analysis (proxies only, see above).
- The UniProt cache (P04): the network was closed, so nothing was cached.
- Enrichment with a partially annotated mixed data set.
- A point covered only by the backdrop in (b).
- User mode (G0 ran in test mode with uploads, not the file chooser).
- A FASTA with more than one organism (the fixture has one).
- Routes other than DIA-NN.
