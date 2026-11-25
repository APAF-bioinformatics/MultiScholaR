# KIT-265 outcome: with 'Mixed species FASTA' ticked, the import never silently keeps an organism nobody chose

Ticket r1.3. Branch `kit-265` from v0.9-staging 82a5df02. Owner decisions: VAL-291 (D-1 (a), D-2) and VAL-304 (O-1:
route (i)(c), the FASTA names organisms and none of the data's proteins is in it, is REFUSED). Result: done; every
prediction of every run held; ready for its implementation audit. Nothing was run on this machine.

## What was done

| Commit | What |
|---|---|
| 2269d03ad | Tests first: `test-kit265-organism-routes.R` (E1, E1b, E1c, E1d, E2, E2c, E2d, E3, E3a, E3b, E4, E5, E6, E7, G3, G4, G6), the three declared edits (KIT-250 G3; two tests of test-prot-01c), R1's predictions. Nothing under R/ |
| c40435efe | The fix (3 files under R/) and the README cell |
| 4a8915ac9 | R1 result, R2 predictions, wording drafts, NEWS, ledger row `VAL-<merge>` |
| 01485d525 | R2 result and tables, runs.tsv, the H1 driver and its predictions (records only) |
| f8c144fbb | H1 run 1 tables; two driver defects fixed (mechanics only) |
| (the commit with this file) | H1 rerun result, this file |

The fix, in the ticket's words: route (i), causes (a) and (b) (no organism field, or no taxon ID on the matched
entries): the import continues and a persistent warning, under the fixed id `prot_import_organism_cancelled`, names
the form's organism and taxon (or says the form has none); route (ii) (an error counting organisms), (iii) (no
protein column: absent from the table on the non-owner path, `protein_ids` NULL on the owner path, no name in the
mapping), (iv) (FASTA processing failed: ticked with no FASTA table) and route (i)(c) (VAL-304: the FASTA has an entry
with organism and taxon ID and none of the data's accessions, split at ";", trimmed, isoform suffix removed, both
sides, is among the FASTA's) raise `rlang::abort(..., class = "prot_import_organism_unavailable")`, which the import's
existing error path shows as "Error: <message>" and uses to clear the import. The (ii) handler passes a refusal through
unchanged. The organism window carries the single-organism sentence; the README cell is qualified. Unticked imports,
the dialog, Confirm and Cancel are not changed.

`R/mod_prot_import_organism_helpers.R` (new: `PROT_IMPORT_ORGANISM_LIMIT_TEXT`, `abortProtImportOrganismUnavailable`,
`buildProtImportNoOrganismNotice`, `protImportFastaLacksDataProteins`; changed: the modal builder, the analysis
function: `organismName` and `taxonId` arguments, the protein-column check, the no-window branch, the error handler),
`R/mod_prot_import_orchestration_helpers.R` (`runProtImportMixedSpeciesAnalysisIfNeeded`: the split condition, the (iv)
abort, the two arguments passed on), `R/mod_prot_import_server_helpers.R` (the call passes `input$organism_name` and
`input$taxon_id`). No `man/` or `NAMESPACE` change (internal functions with `@noRd`).

## Runs

| Run | Binding | Result |
|---|---|---|
| 45698567707 (R1) | tests only 2269d03a, base none | as predicted: new file 67 passed / 52 failed / 0 errors; kit250 118 / 1 (G3); 01c 436 / 2 (the two declared edits); G3 and G4 pass. Stop rules 1 and 2 not fired. 13.0 job-minutes |
| 45700788798 (R2) | base 2269d03a, head 4a8915ac, 10 files, check both | as predicted, per file below. 56.5 job-minutes |
| 45703724770 (H1 run 1) | records 01485d52, commit under test 4a8915ac | 4 of 4 sessions; 37 as predicted, 2 differ, 2 not observed: two instrument defects (below); 13.5 job-minutes |
| 45706090379 (H1 rerun 1 of 3) | records f8c144fb, commit under test 4a8915ac | 41 of 41 as predicted. 13.0 job-minutes |

R2 per file (base -> head): `test-kit265-organism-routes.R` 67/52 -> 119/0 `fixed`; `test-kit250-organism-dialog.R`
118/1 -> 119/0 `fixed`; `test-prot-01c-import-module-contracts.R` 436/2 -> 438/0 `fixed`; `same`: e2e-test-mode-integration
61/0, module-ci-prot-import 107/0, prot-01d 79/0, prot-import-core 38/0, workflow-gui-memory-lifetimes 81/0 (11
skipped, both sides), workflow-server-shared 226/0, workflow-ui-shared 178/0. `R CMD check`: 0 errors, 3 warnings, 3
notes on both sides, 6 findings identical, none only on head or base (`check-compare.txt`).

Instrument defects (VAL-54, H1 only; no result was produced for them; no expected value changed; the ceiling is
unchanged; rerun 1 of 3): (1) session A's Design-tab readings did not run: "could not find function k250DataAvailable"
(a helper in KIT-250's `h1.R`, not in `g0.R`); (2) session B's `import_completed` was read FALSE because the driver's
wait returned on the refusal's own 10 s error toast, still on the page; the same session's other readings already showed
the completed import. Both fixed in `h1.R` (the helper defined; a refusal is a toast new since the press) and re-observed
against the same predictions. The first run's tables are kept in `h1-45703724770/`.

## Acceptance criteria

| # | Criterion | Met | Evidence |
|---|---|---|---|
| A1 | R1 sent after O-1 answered; every E row (E2c for the owner's choice) and the three declared edits fail at their base values by assertion; every G row passes | met | R1: the failing set equals the prediction, 0 errors; O-1 = (a) refuse (VAL-304) |
| A2 | The fix: every E, G row and the declared edits pass | met | R2 head: 119/0, 119/0, 438/0 |
| A3 | R2 base against head: every file `same` except the new file and the two files with declared edits; `R CMD check` identical | met | `compare.tsv`, `check-compare.txt` |
| A4 | H1: the (i) warning names the organism (Mus musculus, 10090) and the import completes; the (iv) refusal shows the error, closes the processing window, shows no success, and after unticking the same session completes | met | H1 rerun 41/41 (session A, B; C unticked as today; D the (c) refusal) |
| A5 | NEWS, ledger row, wording rows, README cell; outcome.md with "For the change report" and "Not checked"; plan rows | met here; plan rows are the coordinator's | NEWS (Fixes), ledger row `VAL-<merge>`, `wording.md`, README.md, this file |

OBSERVED (in the runs above): the stated old text, 5 s default and no id on the code before the fix; the stated new
texts, persistence and ids after; the clearing of a refused import (E5, E6, G6: error shown, `data_tbl` NULL,
`setup_import` incomplete, `design_matrix` still disabled, flags FALSE) on the owner and non-owner paths; the advised
untick completing in the same session (E5 and H1 B); an unticked import with a failing FASTA step completing as today
(G3, H1 C); a persistent notice still on the page 12 s after the import (H1 A, `w1_count` 1). INFERRED: that no number
of an import that completes changes (the R2 `same` verdicts over nine existing files and the unchanged code path for
them are the evidence; no result file was compared).

One addition to the ticket's list of rows: E1d (the observer passes the form's organism and taxon to the notice), a
falsifier for the plumbing in `mod_prot_import_server_helpers.R`; it failed on the code before the fix and passes on it.

## Budget

Core-hours = job-minutes x 4 / 60: R1 13.0 (0.87); R2 56.5 (3.77); H1 26.5 (1.77) = 96.0 job-minutes, **6.40
core-hours** against 9.37 expected and the 64.33 ceiling. R2 came in under its estimate because the ten files ran in
10 to 16 s each (the ticket assumed 1.75 minutes in all, plus install) and the checks took 13 to 16 minutes. The
correction pass was not used; one H1 instrument rerun (1 of 3).

## Divergences and notes

- User messages are sentence-case GUI text, not R.md's lowercase condition style (as KIT-250 recorded). The refusal
  uses `rlang::abort` with a class (not `cli::cli_abort`: the messages embed an error's text and the form's name).
- The (ii) handler's first line still logs "Error analyzing organism distribution: <message>" for a refusal passed
  through from (iii) or (c) (the ticket: "logs as today and then re-signals"); only the log line, not the user's text.
- `tools/refactor/staging/import-server-wave2/` holds an old copy of the organism helpers with the old texts; untouched
  (not part of the package build).
- Unrelated defects: none found.

## Not checked

- No R on this machine; everything above ran on GitHub runners. No real data. No Windows.
- Why the FASTA step failed in (iv): the reason goes only to the log; the message is neutral. H1 used an empty file.
- Whether an unticked import without the protein column fails later at the design step (INFERRED, not tested).
- Whether other organisms' proteins carry gene names in the DA results (the window sentence avoids the claim).
- Whether Shiny sends a notification raised inside the observer before the observer returns (affects when the toast
  appears, not its text); H1 read the toast after the import.
- E3b's base behaviour is a prediction (INFERRED); R1 saw the same failing pattern as predicted (one warning shown, no
  refusal).
- The (c) detection compares accessions; a FASTA whose accessions are written differently from the search engine's is
  also refused with the "could be found" text (for the wording pass).
- A stale `mixed_species_analysis` from an earlier confirmed import after a later no-window import, and a re-import over
  a finished analysis (I07): separate candidates, not made worse; the refusals move I07's "completed with other data"
  form to its "error" form (ticket 3.2).
- The user guide and the paper's sentence (the manuscripts repository) are owed elsewhere and not written here.
- All message and sentence wordings are drafts for the owner's pass (`wording.md`).

## For the change report

**What was wrong and who was affected.** In the proteomics import, with 'Mixed species FASTA' ticked, the user is
told that the organism will be chosen in a window after processing. In four cases no window appears (no organism in the
FASTA file can be matched to the search results' proteins, the organisms could not be counted, the search results have
no protein column, the FASTA file could not be read or processed) and the import finished with the organism typed in
the form, Homo sapiens by default, with at most a short warning that did not name it. That organism then set the
UniProt annotation query, the Taxon ID in the study parameters and the enrichment organism; quality control,
normalisation and differential analysis were not affected. A fifth case (the FASTA names organisms and none of the data's
proteins is in it, usually a FASTA other than the one the search used) behaved the same, with the warning's text false for
it. Probably rare (not known how many users import mixed-species FASTA files).

**Code references.** Base (82a5df02): `R/mod_prot_import_organism_helpers.R:72-76` (no protein column: silent
`return(FALSE)`), `:86-98` (the warning, no organism named, 5 s), `:104-111` (counting error: a warning);
`R/mod_prot_import_orchestration_helpers.R:358-360` (no FASTA table: silent `return(FALSE)`); the FASTA failure itself is
logged at `R/func_prot_annotation_fasta_processing_helpers.R:341-352`. Fix (c40435efe):
`R/mod_prot_import_organism_helpers.R:154-163` (protein column), `:178-197` (route (i)(c) refusal and the persistent
notice), `:204-213` (the handler), `:50-110` (the constants and helpers); `R/mod_prot_import_orchestration_helpers.R:358-372`
(the (iv) refusal); `R/mod_prot_import_server_helpers.R` (the call passes the form's organism and taxon). Present on
`main` (`git show main:R/mod_prot_import.R`, around `:785` to `:883`) since 2025-12-04 (179eb0d2); in the tags v0.3.6.3,
v0.5.0-rc1 and v0.9.0.

**What changed in results.** No number of an import that completes changes. What changes is what the user is told and,
for the four stopped cases, whether the import completes: it now stops, with a message saying what failed and what to do
(untick 'Mixed species FASTA' and enter the organism, or fix the file), and clears what it wrote. For the case with no
organism information in the FASTA headers the import still completes, now with a notice that stays until closed and
names the organism used.

**Evidence.** Runs 45698567707 (R1, failing first), 45700788798 (R2, base against head, check identical), 45703724770
and 45706090379 (H1, a headless browser on a runner on the prot_dia fixture: 400 proteins, one organism, synthetic FASTA
files written by the driver; the rerun 41 of 41 as predicted). Limits: synthetic data only; the empty FASTA file is one
way for FASTA processing to fail; Linux only.
