# Merge train 26 of 2026-11-26: KIT-271

**This train is NOT stacked on other open trains.** It was branched from `origin/v0.9-staging` at
`31a4d684b6382261a14fd5da58d93fe24af5e914` (fetched before building; it holds trains 15 to 18, 19 to 22, 23 and 24 and
the ledger row VAL-344). Its base is therefore the audit's staging 7f17547f plus ledger-only commits.

## Train

Base: `31a4d684b6382261a14fd5da58d93fe24af5e914` (staging).

| Order | Ticket | Branch head | Merge commit | Ledger rows | Audit |
|---|---|---|---|---|---|
| 1 | KIT-271 | `origin/kit-271` 403603f5 (its base 894b9cc8) | 2c3a9c20 | `VAL-348` (the fix) and `VAL-349` (stored references), numbered by the coordinator at landing | AUDIT-KIT-271-impl-20261126 (ready once I01 and I02 land: they are in the merge commit) |

Head (bound by the run): `2c3a9c20cdb08f8decdf36e0fb3fb93d5a9a84b1`. This file is in the declaration commit after it, under
`docs/` only.

### Conflicts met and their resolutions

- **`NEWS.md`**: staging's `### Fixes` bullets kept first, then KIT-271's bullet, with the four precision points of the
  audit's I06 worded plainly (Skip computes no RUV-III factors; in the correction mode the delivered
  `normalised_results.tsv` is the uncorrected table; the old file stays until the normalisation or the correlation filter
  is run again; with the limpa default "normalised" means "before RUV"). They are drafts for the owner's wording pass
  (VAL-41), marked so in `kit-271/wording.md`.
- **`docs/methods-decision-ledger.md`**: every staging row kept (VAL-319, VAL-341 to VAL-344); then KIT-271's row as
  `VAL-348` and the stored-reference row as `VAL-349`, the two exact texts of the audit's section 6 (I02),
  with `<t>` = 26; `<date>` and `<id>` are left for landing. KIT-271's branch row (VAL-<merge>) is replaced by them.
- No conflict under `R/`, `tests/`, `inst/`, `man/`, `DESCRIPTION` or `NAMESPACE`.

### Merged code against the base

`git diff --stat origin/v0.9-staging 2c3a9c20 -- R tests inst man DESCRIPTION NAMESPACE`: six files, 606 insertions, 7
deletions, exactly KIT-271's own change against its base 894b9cc8 (the same six files, 606/7):
`R/func_general_filemgmt_results_summary_helpers.R` (14), `R/mod_prot_norm_observer_helpers.R` (2),
`R/utils_workflow_dia_norm.R` (30), `R/utils_workflow_dia_products.R` (2), `R/utils_workflow_dia_sessions.R` (16) (64
changed lines in R) and the new `tests/testthat/test-kit271-limpa-table-name.R` (549). Nothing under `man/`, `inst/`,
`DESCRIPTION` or `NAMESPACE`; no other test file.

### Record fixes in the merge commit (documents only)

- I01: ledger row and stored-reference row say that on staging the limpa GUI-entry lane already fails at 21 of its 24
  steps since merge train 15, so the file's verdict is message-changed with counts unchanged; one line added to
  `kit-271/outcome.md` under R2 (and a pointer on the step 0 line there).
- I02: the ledger row gains the labels for R1's first send (45709521655, instrument defect, VAL-54) and the refused
  request (45711760556), the VAL-41 drafts clause and the limits; the stored-reference row in the VAL-282/VAL-301 form.
- I06: the four NEWS points (above).
- I08: `outcome.md` ":590-" corrected to ":597 (pass-through; :609 clears the result)"; `runs.tsv` R1 first send
  "0.8" to "1.15" plan/results minutes; `kit-271/r2/compare-failures.txt` (161 KB, 75-byte path with the `MultiScholaR/`
  prefix) committed from run 45714327005's artifact with one `fetch.sh` (it was still available).

## Declaration (written before the run is sent)

**Purpose.** KIT-271's merge check on the merged tree: the KIT-271 scope on staging, base against head, so that the
declared message changes are confirmed on the landing base. The full suite is not repeated.

**Scope recount.** The ticket's grep (`git grep -lE "ruv_normalised_results|RUV_normalised_results|normalised_results_cln_with_replicates|resolveProtNormCorrelationResultFilenames|saveProtNormCorrelationResults|publishDIAReferenceCorrelation|runDIAReferenceLimpaNormalisationPassThrough|results_file|copyToResultsSummary" -- tests/testthat ':!tests/testthat/helper-*' ':!data/benchmark'`)
gives 30 files on the base (staging 31a4d684, counted again here; the audit counted 30 on 7f17547f) and 31 on the merged tree (the 30 and the
new file, which holds the strings). Plus `test-gui-entry-proteomics.R` and `test-prot-norm-worker-figures.R` (neither
matches the grep): **33 files**. Against R2's 32: the same set plus `test-kit259b-second-correction.R` (KIT-259, landed),
minus nothing (set comparison, no other difference).

**Settings.** As KIT-271's R2 request: `env` `MULTISCHOLAR_KIT111_W1_DUMP=/dev/stdout`, `MULTISCHOLAR_KIT112_W1_DUMP=1`;
`blas` reference; `check` true; no extra packages; `file_timeout_minutes` 45 (so `test-prot-norm-core.R` as R2 had it),
`file_timeouts` `test-gui-entry-proteomics.R` 60; **4 shards per side** (R2's 3-shard largest was 72.3 minutes).
Base = `31a4d684b6382261a14fd5da58d93fe24af5e914`, head = `2c3a9c20cdb08f8decdf36e0fb3fb93d5a9a84b1`; one run.

**Ceiling.** 820 job-minutes = 54.67 core-hours (core-hours = job-minutes x 4 / 60); 8 test jobs x 80 + 2 check jobs x 90.
Expected 420 to 450 job-minutes (28.0 to 30.0 core-hours).

### Predictions

Base counts verified against the newest committed tables: train 23's `merge-train-20261126p/compare.tsv` (run
45716988827, head 460253ba, the code of staging 7f17547f), and `merge-train-20261126o` and `q` where train 23 did not
run the file. Counts are passed/failed/errors/skipped.

| File | Base | Head | Verdict |
|---|---|---|---|
| test-general-filemgmt-report-helpers-shared.R | 85/0/0/1 ok (committed table) | as base | same |
| test-gui-entry-proteomics.R | 92/101 (train 23 head, run 45716988827) | 92/101 | message-changed: prot_dia_limpa at its ten steps only, each differing only by the four file-name key lines (see Verification of the masked comparison) |
| test-kit206-prot-norm-export-guard.R | 116/0/0/0 ok (committed table) | as base | same |
| test-kit212-da-current.R | 309/0/0/0 ok (committed table) | as base | same |
| test-kit225-prot-norm-reset-source.R | 56/0/0/0 ok (committed table) | as base | same |
| test-kit244-limpa-pass-through.R | 204/0/0/0 ok (committed table) | as base | same |
| test-kit248-correlation-skip-record.R | 206/0/0/0 ok (committed table) | as base | same |
| test-kit257-norm-own-files.R | 133/0/0/0 ok (committed table) | as base | same |
| test-kit259b-second-correction.R | 90/0 (train 23 head, run 45716988827) | 90/0 | same |
| test-kit271-limpa-table-name.R | absent | 78/0/0/0 | only-on-head |
| test-metab-05-summary-module-characterization.R | 0/0/0/0 ok (committed table) | as base | same |
| test-module-ci-lipid-summary.R | 14/0/0/6 ok (committed table) | as base | same |
| test-module-ci-metab-da.R | 138/0/0/0 ok (committed table) | as base | same |
| test-module-ci-metab-summary.R | 12/0/0/6 ok (committed table) | as base | same |
| test-module-ci-prot-enrich.R | 79/0/0/2 ok (committed table) | as base | same |
| test-module-ci-prot-summary.R | 20/0/0/5 ok (committed table) | as base | same |
| test-omic-da-module.R | read from the run | as base | same |
| test-omics-summary-core.R | 676/25/0/2 failed (committed table) | as base | same |
| test-omic-summary-module.R | 1675/28/0/0 failed (committed table) | as base | same |
| test-prot-01c-import-module-contracts.R | 438/0/0/0 ok (committed table) | as base | same |
| test-prot-05b-norm-module-contracts.R | 440/0/0/22 ok (committed table) | as base | same |
| test-prot-05d-norm-module-characterization.R | 82/0/0/10 ok (committed table) | as base | same |
| test-prot-12-summary-module-contracts.R | 0/0/0/0 ok (committed table) | as base | same |
| test-prot-norm-core.R | 1547/307 (train 23 head, run 45716988827) | 1547/307 | message-changed: exactly the 12 W1 (a) and W2 limpa `files` records (actual 926838a9... -> fe6b0c95... for the two correlation captures, 19bc6774... -> 6c52f9ea... for the four RUV captures, each in process and through the worker); KIT-259 call-count failure (VAL-343) and every other record identical |
| test-prot-norm-input-current.R | 321/0/0/0 ok (committed table) | as base | same |
| test-prot-norm-worker-figures.R | 71/0/0/0 ok (committed table) | as base | same |
| test-workflow-checkpoint-harness.R | 628/0/1/8 failed (committed table) | as base | same |
| test-workflow-integration-input.R | 241/0/0/0 ok (committed table) | as base | same |
| test-workflow-memory-mechanism-member.R | 0/0/0/13 ok (committed table) | as base | same |
| test-workflow-native-session.R | 734/0/0/2 ok (committed table) | as base | same |
| test-workflow-omics-figures.R | 109/0/0/0 ok (committed table) | as base | same |
| test-workflow-ruv-pearson-diagnostic.R | read from the run | as base | same |
| test-workflow-saved-da.R | 755/196 (train 23 head, run 45716988827) | 755/196 | message-changed: exactly the 2 `limpa_one_contrast files` records (8fc2b731... -> d9db4b79...) |
| R CMD check | 0 E / 3 W / 3 N (train 23: 6 findings) | identical | identical findings, nothing only on head |

Files without a committed count are read from the run; for them the prediction is the verdict `same`.

### Verification of the masked comparison (planned)

For `test-gui-entry-proteomics.R` the committed `failures-base.txt` and `failures.txt` records of the run are compared by
script (Python, `-I`, the script kept outside the downloaded folder). Step 1: every record of the file outside
prot_dia_limpa's ten steps (correlation, ruv_covariates, session_export, da_load, da, summary_refusal, summary_save,
summary_copy, summary_report, summary_export) must be identical on both sides after the run's own path masking. Step 2:
for the ten records, delete every line that names one of the four keys
(`results/proteomics/protein_qc/normalised_results_cln_with_replicates.tsv` and its `ruv_` form,
`results_summary/proteomics/Publication_tables/normalised_results.tsv` and `RUV_normalised_results.tsv`), whether as a
`+` or `-` line or as an "is a character vector" or "is absent" line, the chunk header lines
(`names(actual$artifacts)[a:b] vs ...`, whose ranges move with the name), the indented quoted path lines of the
chunks (context or re-sorted copies of the same strings) and blank lines; the rest of each of the ten base records must
equal the rest of the head record. Any residue is reported as a missed prediction. The same script reports, for
`test-prot-norm-core.R` and `test-workflow-saved-da.R`, the records that differ (expected 12 and 2, each by the
`actual` digest only) and that all others are identical.

### Decision and attempts

- Every prediction holds: TRAIN GREEN for KIT-271.
- A missed prediction is a result: recorded, handed back, not repaired and not re-sent. A mechanics failure with no
  result: handed back with the evidence; no re-send without the coordinator's word. An HTTP 500/502 at dispatch with no
  run created: started once by hand with `gh workflow run run-tests.yml --ref run-tests -f request=<name>`.

## What the coordinator resolves at landing

- The two ledger numbers: `VAL-348` (KIT-271) and `VAL-349` (stored references), both rows; `<date>` and
  the train's run id in both status cells (`<id>`).
- VAL-317's status cell: "decided; (1) implemented by KIT-271 (VAL-<n>), merged by merge train 26; (2) KIT-273; (3)
  recorded in VAL-<n+1>"; REC-083 and LIM-104 status cells (audit section 6).
- Placeholders in the branch: `docs/methods-decision-ledger.md` `VAL-348` (the KIT-271 row's id cell and the
  "Stored values ... VAL-348" and "VAL-349" references inside both rows), and in the two rows `<date>`
  and `<id>`; this record (the table above and this list).
- `NEWS.md`: KIT-271's bullet in landing order (it stands after the bullets of staging's `### Fixes` that
  the conflict held, before "Proteomics, normalisation QC"); the wording is the owner's (VAL-41).
- Plan and ticket records of the audit's section 6 (rows `:126`, `:121`, `:128`, I03 and I04 candidates).

## Result

Run 45728201286 (request `train26-compare-2c3a9c20`, request commit facbdfb3), dispatched by the dispatcher without error
(no hand dispatch); green (infrastructure). Head 2c3a9c20 vs base 31a4d684. Tables:
`merge-train-20261126s-compare.tsv`, `merge-train-20261126s-check-compare.txt`,
`merge-train-20261126s-masked-comparison.txt` (same directory).

**Job-minutes 392.5 = 26.17 core-hours** (job-minutes x 4 / 60) against the ceiling of 820 job-minutes (54.67
core-hours); expected 420 to 450 (so under it): tests head 26.9/57.1/31.3/53.1, tests base 56.2/53.0/53.4/30.5, check head
15.3, check base 14.7, plan 0.4, results 0.7 (job start to end from the jobs API).

| File | Predicted | Observed base | Observed head | Verdict |
|---|---|---|---|---|
| test-kit271-limpa-table-name.R | only-on-head, 78/0/0/0 | absent | 78/0/0/0 | only-on-head: HELD |
| test-prot-norm-core.R | 1547/307 both, message-changed in exactly the 12 limpa `files` records | 1547/307 | 1547/307 | message-changed: HELD (12 records on each side; only the `actual` digest differs: 926838a9 -> fe6b0c95 for the two correlation captures and the two W2 repeats, 19bc6774 -> 6c52f9ea for the eight RUV captures; expected 34e64dc1 and b1817922 unchanged; after masking the digest, 0 records differ) |
| test-workflow-saved-da.R | 755/196 both, message-changed in the 2 `limpa_one_contrast files` records | 755/196 | 755/196 | message-changed: HELD (8fc2b731 -> d9db4b79; expected 2ce30818 unchanged; 0 other differences) |
| test-gui-entry-proteomics.R | 92/101 both, message-changed only in prot_dia_limpa's ten steps, key lines only | 92/101 | 92/101 | message-changed: HELD (see the masked comparison) |
| test-kit259b-second-correction.R | 90/0 same | 90/0 | 90/0 | same: HELD |
| the other 28 files | same (counts as the committed tables, or read from the run) | as listed | as listed | same: HELD |
| R CMD check | identical | 0 E / 3 W / 3 N | 0 E / 3 W / 3 N | identical (6 findings on both; only on head: none; only on base: none): HELD |

Totals: 33 files, 29 `same`, 3 `message-changed`, 1 `only-on-head`; no `new-failure`, no `count-changed`, no `fixed`.
Tests head: passed 10938, failed 738, errors 3, skipped 82; base: passed 10860, failed 738, errors 3, skipped 82 (the 78
are the new file). Counts read from the run for the two files without a committed count: test-omic-da-module.R 1288/81
and test-workflow-ruv-pearson-diagnostic.R 301/0 with 2 errors and 5 skipped, both the same on both sides.

**Masked comparison** (script: Python `-I`, kept outside the downloaded folder; output in
`merge-train-20261126s-masked-comparison.txt`). 741 failure records on each side; 717 identical after temporary-path
masking; 24 only on base and 24 only on head, all in the three message-changed files: 12 + 2 + 10. Test-gui-entry file:
every record outside prot_dia_limpa's ten steps is identical on both sides; for the ten steps (correlation,
ruv_covariates, session_export, da_load, da, summary_refusal, summary_save, summary_copy, summary_report,
summary_export) the head records contain 6 key lines (12 at summary_copy, summary_report, summary_export) and the base
records 0; with the key lines, the chunk headers and the quoted path lines removed, **10 of 10 steps are equal to their
base records, no residue**. So the change is the four file-name keys and nothing else, as predicted.

**Decision.** Every prediction held. TRAIN GREEN for KIT-271: the new file passes 78/0 on head; the declared limpa
messages changed in exactly the declared records and by the declared digests; no count and no other file moved; check
identical.

## Rule slips

None found. The run-tests worktree held other agents' untracked folders, untouched. One `fetch.sh` per run (45714327005
for the I08 record, 45728201286 for this result), one `gh run watch`; no polling of the API. The first `gh run list` calls
after the push (three) found the run not yet created, then found it: the dispatcher needs about a minute.

## Merge record (coordinator)

Landed on v0.9-staging on 2026-11-26 by a no-fast-forward merge of this branch at fb49493a onto 7bbc344e (the train was built on 31a4d684; staging gained only ledger rows VAL-345 to VAL-347 since).

- Verdict: green. Run 45728201286 (33 files and R CMD check; 392.5 of 820 job-minutes, 26.17 core-hours): the new file 78/0 only on head; test-prot-norm-core.R (1547/307), test-workflow-saved-da.R (755/196) and test-gui-entry-proteomics.R (92/101) message-changed at unchanged counts, exactly as declared (12, 2 and 10 records; with the digests and the four file-name key lines masked, no record differs); 29 files same; check identical.
- Conflicts at landing: the ledger only. NEWS merged by itself.
- Ledger at landing: staging's rows kept through VAL-347; KIT-271 = VAL-348; the stored-reference row = VAL-349; the placeholders in this record and in kit-271/outcome.md are replaced. Status cells updated: VAL-317; REC-083 (file-name clause closed; the correction-mode order stays open, REC-082); LIM-104 (the DIA_LIMPA file-name item closed).
- Open after this landing, not blocking it (AUDIT-KIT-271-impl-20261126): candidate I03 (a restore refusal treating the two table names as one path, before the release commit); I04 (the S4 download's correction-mode file name, routed to KIT-273); I07 (the limpa GUI lane has no green guard on runners: trains that run test-gui-entry-proteomics.R compare its messages with the key lines masked).
- Unblocked by this landing: KIT-266 and KIT-273 may branch.
