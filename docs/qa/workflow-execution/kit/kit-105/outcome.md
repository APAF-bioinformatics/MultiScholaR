# KIT-105 outcome: the normalisation core moves and one shared lipid/metab normalisation frame (2026-11-08)

**Done.** Lipid and metab now run their normalisation tab on one shared frame (`mod_omic_norm.R`) and panel library (`mod_omic_norm_panels.R`). Each omic is configured by the `norm` group in `omicGuiConfig`.
- **Core moves:** the lipid/metab normalisation computation is in the core (`utils_workflow_omics_{figures,actions}.R`, `utils_workflow_saved_exports.R`). Proteomics' 13 back-referenced normalisation functions, their six callees and the DIA_LIMPA pass-through moved unchanged into `utils_workflow_dia_norm.R`.
- **Back-references:** 25 → **6**, the ticket's target. Neither stop rule fired.
- **Code:** lipid's norm modules went from 2,665 lines to 44 and metab's from 3,852 to 30, with 978 lines of shared frame. Across `R/`, 4,873 lines fewer (2,362 added, 7,235 deleted, against `746fde45`; [lines.tsv](lines.tsv)).
- **Trx** still mounts through the lipid wrappers (KIT-110). **Proteomics** keeps its own frame (P-20, KIT-111).
- **F-20's dead workaround is removed.**
- **Behaviour:** identical for every omic, apart from two display defects found and fixed (F105-01, F105-02).
- **Build:** `kit_105` = `core/browser-smoke.KA61uA` (source `98da5939`).

## Work items

| Item | Commit(s) | What |
|---|---|---|
| W1 | `aa79fe3f` | The `norm` group and its validator. `test-omic-norm-module.R` written first: (a) UI identity, (b) server captures in process, (c) the twin differences. |
| F105-01 | `66838358`, `a2b9b41d` | Lipid/trx normalisation images and labels all showed assay 2 (lazy promises in a loop). Fixed with `force()`. lipid_canonical and trx_screen re-frozen once ([lane-changes.md](lane-changes.md)). |
| W2 | `f8d4c5fc` | The lipid/metab core: one composite and composite-from-files, one final PCA, one in-process sequence, one session collector and saver. Back-references #23 to #28 re-pointed. |
| W3 | `49aa97ce`, `28afdcd4` | Back-references #10 to #22, their callees and the pass-through moved unchanged into `utils_workflow_dia_norm.R`. F-20 removed and its block retargeted. Collate updated. |
| W4 | `120efa7f`, `54bef4ab` | The shared frame and panel library with lipid on it. Trx's 11 helper calls re-pointed. Six lipid helper files deleted. In-process correlation twin added. F105-02 fixed for lipid. |
| W5 | `fcc16c6c`, `f1f1dc3a`, `98da5939` | Metab on the frame: five helper files and the 640-line hand-written UI deleted. Metab's export progress, summary time and correlation-failure notice are config data. F105-02 fixed for metab; F105-03. A live core-only block keeps the quarantined comparisons' core-side expectations. |
| W6 | `8d58c9db`, `96f80fd8` | The in-process twins and the export collector and saver are now required config fields (the test failed first). Quarantine rows. Merge checks. |

## Proof

| Check | Result |
|---|---|
| Merge checks ([merge-checks.txt](merge-checks.txt)) | `Collate` 361/361. NAMESPACE is byte-identical; `man/`, `app_ui.R`, `app_server.R`, `mod_omic_workflow.R` and `tools/` have no diff. Back-references are 6. The dependents grep, computed against `746fde45`, found no live caller of the deleted names or the 11 deleted files. |
| Simulation Setup ([entry-suite.md](entry-suite.md)) | `metab_combined` alone gives 35/0 and goes through `mod_omic_norm_server`, the in-process normalisation, the correlation, the collector and the saver once each. `prot_dia_limpa` alone gives 37/0 and goes through the moved pass-through once. The first attempt was a harness failure (wrong working directory), recorded there. |
| Run 1, targeted ([targeted.txt](targeted.txt)) | 46 files, 7,132 passed, 0 failed, 0 errors, under strict skips. Every skip is recorded. |
| Run 2, entry suite ([entry-suite.md](entry-suite.md)) | All 10 lanes identical to their oracles, 0 fail/error/skip. |
| Run 3, DIA session and memory set ([dia-memory-set.txt](dia-memory-set.txt)) | 733/2, 191/1, 38/1, 174/4, 80/11: each equal to its baseline |
| Run 4, `R CMD check` ([check.log](check.log)) | 3 WARNINGs, 3 NOTEs. The sections and the code NOTE are identical to KIT-104's. |
| Run 6, journeys on `kit_105` ([member-results.txt](member-results.txt)) | 105a lipid, 105b metab, 105c trx, 105m DIA_LIMPA, 105l LFQ, 105t TMT: all accepted under `exactness_and_completion` |
| Normalisation-phase memory ([norm-peaks.txt](norm-peaks.txt)) | Every phase within ±1.6% of KIT-104 |
| Run 7, e2e lanes ([e2e-final.txt](e2e-final.txt)) | lipid 88, metab 190, trx 14, prot DIA 42, DIA_LIMPA 79, enrichment 98, TMT/LFQ 153: 0 failures, 0 skips, each equal to KIT-104 |
| Run 8, skip comparison ([skip-compare.md](skip-compare.md)) | Against KIT-104, each change is one of: a recorded quarantine (Q-86 to Q-89, Q-95, Q-96, Q-105 to Q-111); one of three new files; or the retargeted F-20 block (05b, 504 → 502). The new `test-prot-norm-core.R` has 3 check-mode errors, the known source-tree read of the shared DIA checkpoint fixture (its source file has the same 14 at KIT-104). Check mode: 18 fewer failures, 1 fewer error. |

## Findings, decisions and parked items

- **F105-01** (fixed): see the work items. Display only.
- **F105-02** (fixed): an error in skip-correlation's store commit ended the user's session. The frame now catches it, shows an error notification and keeps the state. Only the two W1 (c) unit captures were re-frozen; no entry lane changed.
- **F105-03** (fixed): a clock flake in a test's compared value (test only).
- **Divergences:**
  - Metab logs a few console lines it did not before: the correlation start, the export collection, and a refused skip commit.
  - Some edge cases the tab cannot reach now follow lipid's handling (`w5/notes-for-w6.txt`).
  - No capture, state, file, log line in the tab, or notification differs.
- **Parked from this ticket:** P-18 (the export's per-assay sample count), P-19 (the in-process vs worker ITSD resolution), and P-20 (a lead decision, for review: proteomics keeps its own frame until KIT-111).
- **Hand-offs:**
  - **KIT-110:** trx's `norm` values on the frame. The lipid wrappers' `profile` branches and `mod_trx_norm.R:217-227` are its start.
  - **KIT-111:** proteomics' norm core is in `utils_workflow_dia_norm.R`. `generateProtNormCompositeFromFiles` can use `omicsNormCompositeFromFiles`.
  - **KIT-106:** back-references #29 to #32.
  - **KIT-109:** an impact-map rule for `mod_omic_norm*.R`.

## Rollback

`git revert` of the KIT-105 commits, in reverse order. Restore `tests/testdata/gui-entry/{lipid_canonical,trx_screen}.json` from `a2b9b41d^` only if F105-01 itself is reverted. The exported module functions keep their signatures. The quarantined tests are restored by removing their Q-86 onwards skip lines.
