# KIT-273 outcome

Branch kit-273 from staging 7b4ae7dad (KIT-271 merged). Records: step0.md, predictions-r1.md, predictions-r2.md, runs.tsv, wording.md.

## What was run

| Run | Commit | Result |
|---|---|---|
| R1 45736368781 | a5c66eb7 (tests only) | 84 passed, 13 failed, 1 error. 11 failures are the predicted rows at the stated values. 2 failures and the 1 error were mechanics of the test file, no result judged: T9 compared Copy's listing before and after the DA (the DA's own files differ), T10's native-session fixture made Copy read an archive pointer no export had written. |
| R1b 45738622121 | 2d46e4a8 (tests only, instrument run 1 of 3) | 93 passed, 11 failed, 0 errors: exactly the 11 predicted failures (T1, T7 covariates and T2, T7 correction, in process and worker; T6; T11; T12). T3, T4, T5, T9, T10, T13 pass. The repair removed one expectation (T9's listing comparison across the DA) and set the T10 fixture non-native; no expected value, threshold or comparison otherwise changed. |
| R2 45740193086 | c9caf84d against 2d46e4a8 | 46 files and R CMD check; see below. |

Core-hours: R1 1.19, R1b 0.77, R2 25.6: 27.6 against expected 31.08 and ceiling 103.0 (no correction pass, one instrument run).

## Acceptance criteria

| # | Evidence |
|---|---|
| A1 met | R1b: the 11 failing rows fail by value (file absent / no `ruv_delivery` record); guards pass; the base Copy listings are in the logs (`KIT-273 T1/T2/T6/T10` lines). |
| A2 met | R2 head: new file 130 passed, 0 failed, in process and worker. |
| A3 met | R2: 44 files same, new file fixed, test-gui-entry-proteomics.R message-changed as declared; standing failures identical (test-prot-norm-core.R 1547/307, test-workflow-saved-da.R 755/196, test-gui-entry-proteomics.R 92/101); R CMD check identical (ERRORs 0, WARNINGs 3, NOTEs 3 on both; no finding only on head). |
| A4 met | the constants are the owner's names (VAL-332). |
| A5 met | T2: protein ids and samples of the delivered matrix equal the DA's E and values equal on the delivered text (correction); T8: byte-identical to the protein-QC files; T1: the delivered design's W terms equal `da_model$w_terms` and the DA source's design; T2b: `normalised_results.tsv` equals the DA's E (passes on base and head, no finding about KIT-271); T5, T6, T13 (no stale file, nothing without a DA). |
| A6 | NEWS, ledger row `VAL-<merge>`, this file. |

Observed: the gui-entry lane records, base against head, for prot_dia_limpa differ at summary_copy, summary_report and summary_export only, by the key `Study_report/design_matrix_with_ruv_factors.tab` (and the length and chunk lines it moves). With that key, the length lines, chunk headers and quoted path context removed, base and head records are equal at all three steps (script, Python -I, residue none). The four KIT-271 file-name key lines are identical on both sides. Every other lane record is identical; the design_matrix.tab digest does not move.

## Stored values that move (none edited; for the end-of-project re-record, VAL-279)

`tests/testdata/gui-entry/prot_dia_limpa.json`: three Copy-step listings gain one key each (`results_summary/proteomics/Study_report/design_matrix_with_ruv_factors.tab`, covariates mode), near `:4989-4994` (summary_copy), `:5324-5332` (summary_report), `:5667-5675` (summary_export) as the ticket read them (before KIT-271 shifted nothing in line count; re-read at the re-record). The digest of each new key is to be the lane's stored `protein_qc/ruv_covariates_design.tsv` digest if byte-identical; T8 shows byte identity in the fixture.

## Implementation notes

- The failure record is made inside `writeDIARuvDeliveryTables` (it names the file it was about to write), not by a `tryCatch` around the call in Copy as the ticket's wording had it; Copy appends the `failed_copies` entry (type `ruv_delivery`, source, destination, display name, error) from the returned value. Observable behaviour is the ticket's.
- Value comparisons in the tests are on the text the package writes (same writer), tolerance 0 on that representation.
- T4's "declared design that applied no RUV (k = 0)" is covered through the Skip case: the DA record is the same (`ruv_mode` none, no W terms), so the decision branch is the same; a k = 0 simulation fixture needs sva and the rules gate and was not built.
- The Skip-RUV state (KIT-276) is untouched; NEWS says what is true in it.
- Scope rule: no R/ file matching the omic-module pattern was edited; test-omics-module-globals.R not run.

## I04 (routed by AUDIT-KIT-271-impl)

The ticket's scope as written does not cover it (Copy's store-based delivery; the S4 download is `protNormS4DownloadHandler`, `R/mod_prot_norm_observer_helpers.R`, outside the ticket's files). Not changed. Candidate, smallest repair: when the head state is `ruv_corrected`, write the canonical file under `ruv_normalised_results_cln_with_replicates.RDS` instead of `names$rdsFilename`.

## Not checked

The tar delivery with the new call on a native export (the call does not touch it; not exercised); the k = 0 fixture; a restored v0.9.0 archive; Windows; real data; byte equality of in-process and worker files beyond T1/T2/T8 digests (both paths passed those).

## For the change report

- What was wrong: on DIA with limpa the results package held neither the matrix the analysis fitted in the RUV-III correction mode (only the table before correction) nor, in the covariates mode, the design columns W1..Wk the model contained. Every DIA-with-limpa project in 0.9 using a RUV mode was affected; not on `main`.
- Code: base `R/func_general_filemgmt_results_summary_helpers.R:371-427` (Copy's proteomics block, no `ruv_*` file delivered) and the writers `R/utils_workflow_dia_products.R:1802,1805` (protein_qc only); fix c9caf84d adds `writeDIARuvDeliveryTables` after `publishDIAReferenceRuvCovariates` and one call after the native-session artifact in Copy.
- Results: two optional files are added to the package in the stated modes; no value of any existing table moved.
- Evidence: runs 45736368781, 45738622121, 45740193086; the standard errors and imputed flags remain only in the session archive (candidate).
