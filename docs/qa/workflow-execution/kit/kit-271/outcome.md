# KIT-271 outcome: DIA with limpa, the correlation filter's table is no longer named as an RUV-normalised table

Branch `kit-271` from origin/v0.9-staging 894b9cc821a70b03baab1e1fa47b3d47440c7ed5. Tests-only commit
4a80353ccce5452ed8096303901646ed92a1faee (after 52f09f10, its first form); fix commit 23d0d15881c6dabc2e0607042fa40416f29c2c1c.
Owner's decision: ledger VAL-317 (Q1 = N1). Decision record, predictions and step 0: `step0.md`. Runs: `runs.tsv`. Wording drafts: `wording.md`.
OBSERVED = seen in a returned run table or log (run id given). INFERRED = follows from what was read.

## What was done

- `R/utils_workflow_dia_norm.R:779` `resolveProtNormCorrelationTableNames(manager, ...)`: on `workflow_type == "DIA_LIMPA"` it calls the
  pinned resolver once with `list(ruv_skipped = TRUE)` (names `normalised_...`, console line silenced), writes its own console line and
  returns the resolver's list; on every other route it returns the resolver's value unchanged (one resolver call on every route). The
  resolver `:742-767` and its pinned digest (`test-prot-norm-core.R:31`) are not edited.
- Callers: `utils_workflow_dia_products.R:1976` (the correlation's writer), `mod_prot_norm_observer_helpers.R:226` (S4 compatibility
  download), `utils_workflow_dia_sessions.R:1925` (export metadata). The in-memory writer and the legacy ingest are untouched.
- Older projects (R06 placement): `utils_workflow_dia_sessions.R:1684-1692` (`collectDIANativeGuiFiles`): on a DIA-with-limpa owner the
  export's file list takes the `normalised_` file if it is on disk, else the `ruv_` one, else the name the metadata gives; covers the
  export, the Summary's archive refresh and metadata from builds before KIT-206.
- Copy, one branch: `func_general_filemgmt_results_summary_helpers.R:380-393`: on a DIA-with-limpa active owner, the `normalised_` table
  if on disk, else the `ruv_` one, is delivered as `Publication_tables/normalised_results.tsv` ("Normalized Results TSV"; the RDS entry follows).
- NEWS (`### Fixes`, the single unreleased heading) and ledger row `VAL-<merge>` (citing VAL-317). No version change. No stored value edited.

## Acceptance table

| # | Criterion | Result | Evidence |
|---|---|---|---|
| A0 | Step 0: the eight stored places classified | met | `step0.md`: the six W1 limpa `files` and `limpa_one_contrast files` already fail on runners; the gui-entry limpa lane passed (the ticket expected it to fail already) ; run 45650328815, fetched once |
| A1 | Tests-only commit; failing-first at the stated values | met | run 45712056817 (4a80353c): 56 expectations, 37 passed, 19 failed, 2 skipped. Failed: T1 (2, both paths), T2 name row (2), T3 (2), T4 (4), T5 (3), T5b results_file and destination rows (4), T8 Copy (2): each at the stated value (`present` = `ruv_normalised_...tsv`; `"normalised_results.tsv" %in% copy$tables` FALSE; `results_file` = `ruv_...`). Passed: T6 a, T7, T8 export, the T2 after==before / content / adjusted-matrix rows, T4's digest row, T5/T5b/T8 export error rows. Skipped: T6 b, T8b (function absent). Stop rule 1 did not fire (the base name is `ruv_`) |
| A2 | The fix: all rows pass, in process and worker | met | run 45714327005 head: `test-kit271-limpa-table-name.R` 78 passed, 0 failed, 0 skipped (T6 b and T8b run) |
| A3 | R2 base against head: every file `same` except the declared ones; check identical; head values dumped | met | below |
| A4 | The name is the owner's answered one | met | N1, VAL-317 |
| A5 | The table's values identical before and after | met | the table's SHA-256 is `156903a4c09b186926465be29477acf60b37f0c96b51c69cdfadba25cc7e2462` in R1 (name `ruv_`, both paths, and before/after the RUV action in both modes) and in R2 (name `normalised_`); the delivered copy has the same digest on both sides (T4) |
| A6 | NEWS, ledger row, this file | met | NEWS.md, `docs/methods-decision-ledger.md` (`VAL-<merge>`), this file |

## R2: base (tests-only 4a80353c) against head (fix 23d0d158), run 45714327005, 32 files, 3 shards per side, check on both

Prediction (`step0.md`) against the returned `compare.tsv`:

| File | Predicted | Observed |
|---|---|---|
| `test-kit271-limpa-table-name.R` | `fixed` | `fixed` (19 failures to 0; 37 to 78 passed) |
| `test-prot-norm-core.R` | `message-changed`, counts same, only the six `files` | `message-changed`, 1548/306 both sides; per-expectation diff: only `limpa_covariates_correlation`, `_ruv`, `_ruv_repeat`, `limpa_correction_correlation`, `_ruv`, `_ruv_repeat` `files`, in process and worker (12); `calls`, the `:31` definition block and the pass-through captures unchanged |
| `test-workflow-saved-da.R` | `message-changed`, `limpa_one_contrast files` | `message-changed`, 755/196 both sides; only `limpa_one_contrast files` (in process, worker) |
| `test-gui-entry-proteomics.R` | `new-failure`, ten steps of prot_dia_limpa, 168/25 to 158/35 | `new-failure`, 168/25 to 158/35; lane prot_dia_limpa 38 to 28 passed; the ten steps are correlation, ruv_covariates, session_export, da_load, da, summary_refusal, summary_save, summary_copy, summary_report, summary_export; the differences are the key names only (`normalised_results_cln_with_replicates.tsv` for `ruv_...`, `normalised_results.tsv` for `RUV_normalised_results.tsv` at three steps); no other lane |
| `test-kit244-limpa-pass-through.R`, `test-kit257-norm-own-files.R`, `test-prot-norm-worker-figures.R` | `same` | `same` |
| the other 26 files | `same` | `same` |
| R CMD check | identical | base and head: 0 ERRORs, 3 WARNINGs, 3 NOTEs; findings identical on both (6), none only on head, none only on base |

Failure counts: base 682 failing expectations, head 673 = 682 + 10 (gui-entry) - 19 (new file), and nothing else.

Head values of the declared stored places (from the dumps in the head job logs; runner values, which differ from the laptop-recorded
frozen values on base and head alike, so they are not the values to record at the end-of-project re-record; they show what moved):

| Place | Base (runner) | Head (runner) |
|---|---|---|
| W1 `limpa_*_correlation` `files` (`:862`, `:992`) | `926838a97f9b...` | `fe6b0c95b1e9...` |
| W1 `limpa_*_ruv`, `_ruv_repeat` `files` (`:897`, `:933`, `:1027`, `:1063`) | `19bc6774b76d...` | `6c52f9ea998f...` |
| W1 `limpa_*_pass_through` `files` (not declared) | `8ab794fecc83...` | `8ab794fecc83...` (unchanged) |
| `limpa_one_contrast files` (`test-workflow-saved-da.R:1136`) | `8fc2b7310ced...` | `d9db4b793014...` |
| `prot_dia_limpa.json`, 13 lines | `ruv_normalised_results_cln_with_replicates.tsv` (ten steps), `.../Publication_tables/RUV_normalised_results.tsv` (three) | `normalised_results_cln_with_replicates.tsv`, `normalised_results.tsv`; the digest `1a61f3bd...` on the lines is unchanged |

## By mode: what a user receives against what was analysed (DIA with limpa; INFERRED from the code order, A5 and T2 observed the table's identity)

| RUV mode | What the analysis used | The delivered table, before | The delivered table, now |
|---|---|---|---|
| Skip | the filtered table; design without RUV terms | that table, as `RUV_normalised_results.tsv` | that table, as `normalised_results.tsv` |
| Factors as covariates (default) | the filtered table with RUV-III's W factors in the design | that table, as `RUV_normalised_results.tsv`; the W columns (`ruv_covariates_design.tsv`) are not delivered | that table, as `normalised_results.tsv`; the W columns still not delivered (KIT-273) |
| Correction | RUV-III's corrected matrix (`ruv_adjusted_protein_matrix.tsv`) | the UNCORRECTED filtered table as `RUV_normalised_results.tsv`; the matrix the analysis used is delivered nowhere | the same uncorrected table, now correctly named `normalised_results.tsv`; the analysed matrix is still not delivered (KIT-273) |

OBSERVED: the table's SHA-256 does not change across the RUV action in either mode (T2, R1 and R2), and `ruv_adjusted_protein_matrix.tsv`
is written as a separate file (T2 guard).

## Deviations and instrument records

- Instrument defect (brief rule 3, VAL-54), run 45709521655: the export helper of the test file wrote only the session summary, so the
  fixture's native export failed in staging ("saved_entity_corrupt: saved table changed during read": `qc_params.RDS` named but not written);
  no result was judged. Fixed in 4a80353c by calling the app's own `saveProtNormExportMetadataFiles` in the prepare step; no expected value,
  threshold or comparison changed; one repair, ceiling unchanged. The six affected rows (T5b x4, T8 export x2) were predicted to pass on base and did after the repair.
- Request mistake, run 45711760556: refused at `plan` (the `requested_by` text was over 300 characters); no test ran; re-sent as
  `kit271-r1-3-4a80353c` (45712056817). Not a test result and not counted against the instrument allowance.
- The step 0 finding that the gui-entry limpa lane passes on runners (the ticket expected it to fail already) changed that file's declared
  verdict from `message-changed` to `new-failure`; the prediction was committed with that verdict before R1 and R2.
- The GitHub API rate limit (403) cut off the first watch of R2; the run was not affected; a new single watch was started after the reset.
- One working-tree detail: `test-kit271-limpa-table-name.R` first form (52f09f10) and its repair (4a80353c) are both on the branch; the R2 base is 4a80353c.

## Core-hours

Expected 30.93, ceiling 103.0. Used (job minutes x 4 / 60, plan, test, check and results jobs): R1 first run 13.4 min, refused request 0.1,
R1 re-run 14.5 (together 28.0 min = 1.86 core-hours); R2 397.4 min (26.49). Total about 425 min = 28.4 core-hours. Step 0 fetched an existing run (no cost).
No correction pass was used.

## Not checked

- No R on this machine; every statement about behaviour rests on the returned runner tables.
- The real application (Shiny session, Summary tab, a restored v0.9.0 archive): the older-project routes are tested on a fixture (the old file
  planted on a native-flagged owner), not through a real restore of an archive written by v0.9.0 or by a build between. The reopen (not restore) route.
- The Summary's archive refresh (`refreshDIANativeSummaryArchive`) with stored metadata was not run; it reuses the changed file list by reading.
- The report templates (`DIANN_limpa_report.rmd:1404` reads `ruv_normalised_results.RDS`, a miss in native sessions) were not run and are unchanged.
- The end-of-project values of the eight stored places (runner values differ from the laptop-recorded ones; the re-record is VAL-279).
- The label text on Copy's list and the console lines are drafts (VAL-41), not seen in the GUI. Windows.
- Whether Copy's `design_matrix.tab` carries the W terms in the covariates mode (ticket F6); what else Copy delivers (KIT-273).

## For the change report

What was wrong and who was affected. On DIA with limpa (the route that computes RUV-III's factors after the correlation filter) the table
the filter writes is the protein table before any RUV step, in the none, covariates and correction modes alike. It was written as
`protein_qc/ruv_normalised_results_cln_with_replicates.tsv` and delivered by Copy as `Publication_tables/RUV_normalised_results.tsv`,
a name that says RUV was applied. In the correction mode it also implied that the delivered table was the matrix the differential
analysis used, which it was not (the corrected matrix, `ruv_adjusted_protein_matrix.tsv`, is delivered nowhere). Everyone who ran DIA with limpa
to Copy or Summary in tag v0.9.0, staging and every build since was affected; nobody on `main`, which has no limpa GUI route (`git grep DIA_LIMPA main -- R`
finds nothing; read with `git show main:<path>` only); nobody on the DIA, label-free and TMT routes, where the name follows the RUV choice.

Cause and code references. The resolver `resolveProtNormCorrelationResultFilenames` treats an absent RUV result as "applied": base 894b9cc8
`R/utils_workflow_dia_norm.R:742-767` (`:748-760`). The limpa pass-through clears the result (`:716`) and RUV runs after the filter
(`R/mod_prot_norm_workflow_helpers.R:79-82`), so the result is always absent at the filter. Callers at base: `R/utils_workflow_dia_products.R:1976`,
`R/mod_prot_norm_observer_helpers.R:226`, `R/utils_workflow_dia_sessions.R:1917`; Copy `R/func_general_filemgmt_results_summary_helpers.R:373-405`.
At tag v0.9.0 the same logic is at `R/utils_workflow_dia_norm.R:635-658` (resolver) and `:590-` (pass-through). On `main` only the limpa workbook wrote the
`ruv_` name, after peptide-level RUV, where it was true (inst/workbooks; not touched).

What changed. At fix 23d0d158: `R/utils_workflow_dia_norm.R:769-798` (the new function; the resolver is unchanged), the three callers
(`utils_workflow_dia_products.R:1976`, `mod_prot_norm_observer_helpers.R:226`, `utils_workflow_dia_sessions.R:1925`), the export's file list
(`utils_workflow_dia_sessions.R:1684-1692`) and Copy (`func_general_filemgmt_results_summary_helpers.R:380-393`). The table is stored as
`normalised_results_cln_with_replicates.tsv` and delivered as `normalised_results.tsv` (label "Normalized Results TSV"), the name every other
proteomics route gives a table RUV did not correct. A project written by an earlier build needs no re-run: the export lists the file on disk and Copy delivers it under the new name.

Results. No value changed: the table's SHA-256 is identical before and after (`156903a4...`, runs 45712056817 and 45714327005, both paths, all
modes). The correction mode's analysed matrix and the covariates mode's W columns are still not delivered (KIT-273).

Evidence. Run 45712056817 (R1, 19 failing-first rows at the base values), run 45714327005 (R2, 78 of 78 pass at head; base-against-head table above; check identical).
Limits: fixture projects of 12 samples x 64 proteins on runners; the app, a real restore and the report were not run. The new tests ran in process and through the worker.
