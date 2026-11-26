# KIT-273 outcome

Branch kit-273 from staging 7b4ae7dad (KIT-271 merged). Records: step0.md, predictions-r1.md, predictions-r2.md, runs.tsv, wording.md.

## What was run

| Run | Commit | Result |
|---|---|---|
| R1 45736368781 | a5c66eb7 (tests only) | 84 passed, 13 failed, 1 error. 11 failures are the predicted rows at the stated values. 2 failures and the 1 error: T9 compared Copy's listing before and after the DA (the DA's own count table, `da_proteomics_num_sig_da_molecules.tab`, differs; a result, see the ruling below), T10's native-session fixture made Copy read an archive pointer no export had written (an instrument defect, an error before any expectation). |
| R1b 45738622121 | 2d46e4a8 (tests only, instrument run 1 of 3) | 93 passed, 11 failed, 0 errors: exactly the 11 predicted failures (T1, T7 covariates and T2, T7 correction, in process and worker; T6; T11; T12). T3, T4, T5, T9, T10, T13 pass. R1b removed one expectation (T9's listing comparison across the DA: a result, ruled I01, below) and set the T10 fixture non-native (instrument run 1 of 3). Corrected 2026-11-26 (I01): the earlier wording "no expected value changed" was not accurate; the removed T9 expectation was an expected value. |
| R2 45740193086 | c9caf84d against 2d46e4a8 | 46 files and R CMD check; see below. |

Ruling on the R1b removal, quoted (AUDIT-KIT-273-impl-20261126 I01):

> The removal of the T9 expectation after R1 was the correction of a mis-specified expectation, not a repair of test mechanics: the expectation compared Copy's listing before the differential analysis with the listing after it, and the analysis's own count table makes them differ on the code before and after the fix alike. R1b therefore counts as the ticket's one correction pass as well as instrument run 1 of 3 (T10's fixture); the correction pass is spent. Stop rule 1 named T9 and was passed over without a ruling; the removal was disclosed. What T9 was written to guard (existing deliveries unchanged by the fix) is shown by R2: the new file's listings and the limpa lane's masked comparison. (Coordinator, 2026-11-26, on AUDIT-KIT-273-impl-20261126 I01.)

Core-hours (jobs API, audit 2026-11-26): R1 17.77 job-minutes = 1.18, R1b 11.32 = 0.75, R2 383.78 = 25.59: 412.87 job-minutes = 27.52 core-hours against expected 31.08 and ceiling 103.0. The ticket's correction pass is spent (R1b, ruling above); one instrument run of 3 used.

Run tables are committed beside this file: `r1-tests.tsv`, `r1-failures.txt`, `r1b-tests.tsv`, `r1b-failures.txt`, `r2-compare.tsv`, `r2-tests-base.tsv`, `r2-tests-head.tsv`, `r2-check-compare.txt`, `r2-compare-failures.txt` (copied unchanged from the runs' `results` artifacts).

## Acceptance criteria

| # | Evidence |
|---|---|
| A1 met | R1b: the 11 failing rows fail by value (file absent / no `ruv_delivery` record); guards pass; the base Copy listings are in the logs (`KIT-273 T1/T2/T6/T10` lines); the base listings print file names only, the files' digests are in the lane's comparison. |
| A2 met | R2 head: new file 130 passed, 0 failed, in process and worker. R2 predicted 104/0 and observed 130/0: the 26 expectations inside `if (file.exists(...))` run only on head; the verdict `fixed` was predicted. |
| A3 met | R2: 44 files same, new file fixed, test-gui-entry-proteomics.R message-changed as declared; standing failures identical (test-prot-norm-core.R 1547/307, test-workflow-saved-da.R 755/196, test-gui-entry-proteomics.R 92/101); R CMD check identical (ERRORs 0, WARNINGs 3, NOTEs 3 on both; no finding only on head). |
| A4 met | the constants are the owner's names (VAL-332). |
| A5 met | T2: protein ids and samples of the delivered matrix equal the DA's E and values equal on the delivered text (correction); T8: byte-identical to the protein-QC files (digests in R2's head, both paths: correction `5057296d916c6b55a7cdaa7ee6baebfd3d926bee84384b38439bf45e96e41947`, covariates `b65f0f63761679b3be8da828b2783323c32b7c16c8100d4af8d67158d52bdcba`); T1: the delivered design's W terms equal `da_model$w_terms` and the DA source's design; T2b: `normalised_results.tsv` equals the DA's E (passes on base and head, no finding about KIT-271); T5, T6, T13 (no stale file, nothing without a DA). |
| A6 | NEWS, ledger row `VAL-<merge>`, this file. |

Observed: the gui-entry lane records, base against head, for prot_dia_limpa differ at summary_copy, summary_report and summary_export only, by the key `Study_report/design_matrix_with_ruv_factors.tab` (and the length and chunk lines it moves). With that key, the length lines, chunk headers and quoted path context removed, base and head records are equal at all three steps (script, Python -I, residue none). The four KIT-271 file-name key lines are identical on both sides. Every other lane record is identical; the design_matrix.tab digest does not move.

## Stored values that move (none edited; for the end-of-project re-record, VAL-279)

`tests/testdata/gui-entry/prot_dia_limpa.json`: three Copy-step listings gain one key each (`results_summary/proteomics/Study_report/design_matrix_with_ruv_factors.tab`, covariates mode), near `:4989-4994` (summary_copy), `:5324-5332` (summary_report), `:5667-5675` (summary_export) as the ticket read them (before KIT-271 shifted nothing in line count; re-read at the re-record). The digest of each new key is `31bee7fa8d6bf07efd85dac297b1a9148868c8cdf07fdad4fdf759bd4593aeb0` (R2), the stored `protein_qc/ruv_covariates_design.tsv` digest of the lane; T8 shows byte identity in the fixture.

## Implementation notes

- The failure record is made inside `writeDIARuvDeliveryTables` (it names the file it was about to write), not by a `tryCatch` around the call in Copy as the ticket's wording had it; Copy appends the `failed_copies` entry (type `ruv_delivery`, source, destination, display name, error) from the returned value. Observable behaviour is the ticket's.
- Value comparisons in the tests are on the text the package writes (same writer), tolerance 0 on that representation.
- T4's "declared design that applied no RUV (k = 0)" is covered through the Skip case: the DA record is the same (`ruv_mode` none, no W terms), so the decision branch is the same; a k = 0 simulation fixture needs the rules gate and a declared-rule fixture (sva is installed on the runners) and was not built; covered by reading (AUDIT-KIT-273-impl 3.4).
- The Skip-RUV state (KIT-276) is untouched; NEWS says what is true in it.
- Conventions divergence (AUDIT-KIT-273-impl I10): the new `tryCatch` catches every `error` (R.md `:118` asks for specific conditions); this is the form the ticket prescribed, and it logs the failure.
- Scope rule: no R/ file matching the omic-module pattern was edited; test-omics-module-globals.R not run.

## I04 (routed by AUDIT-KIT-271-impl)

The ticket's scope as written does not cover it (Copy's store-based delivery; the S4 download is `protNormS4DownloadHandler`, `R/mod_prot_norm_observer_helpers.R`, outside the ticket's files). Not changed. Candidate, smallest repair: when the head state is `ruv_corrected`, write the canonical file under `ruv_normalised_results_cln_with_replicates.RDS` instead of `names$rdsFilename`.

## Not checked

Correction mode on a native session (the new test file's limpa fixtures are not native; the GUI-entry limpa lane is native and covers the covariates mode: R2 head adds the key beside the unchanged `filtered_session.msdia.tar` key); the k = 0 fixture (by reading); a restored v0.9.0 archive; Windows; real data; byte equality of in-process and worker files beyond T1/T2/T8 digests (both paths passed those).

## For the change report

- What was wrong: on DIA with limpa the results package held neither the matrix the analysis fitted in the RUV-III correction mode (only the table before correction) nor, in the covariates mode, the design columns W1..Wk the model contained. Every DIA-with-limpa project in 0.9 using a RUV mode was affected; not on `main`.
- Code: base `R/func_general_filemgmt_results_summary_helpers.R:371-427` (Copy's proteomics block, no `ruv_*` file delivered) and the writers `R/utils_workflow_dia_products.R:1802,1805` (protein_qc only); fix c9caf84d adds `writeDIARuvDeliveryTables` after `publishDIAReferenceRuvCovariates` and one call after the native-session artifact in Copy.
- Results: two optional files are added to the package in the stated modes; no value of any existing table moved.
- Evidence: runs 45736368781, 45738622121, 45740193086; the standard errors and imputed flags remain only in the session archive (candidate).
- Tag v0.9.0 carries the same gap: Copy's proteomics block delivers no `ruv_*` file (`v0.9.0:R/func_general_filemgmt_results_summary_helpers.R:373-410`), the two tables are written to `protein_qc` only (`v0.9.0:R/utils_workflow_dia_products.R:1548`, `:1551`), and the correction mode fits the corrected E (`:1046`). `main` (0.4.1.2) has no limpa GUI route. Fix c9caf84d: `writeDIARuvDeliveryTables` (`R/utils_workflow_dia_products.R:1832-1870`) and Copy's call (`R/func_general_filemgmt_results_summary_helpers.R:427-438`).
