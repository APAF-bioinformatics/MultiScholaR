# KIT-110 outcome: transcriptomics' normalisation on the shared frame (2026-11-09)

**Done.** Transcriptomics' normalisation tab now runs on KIT-105's shared frame, configured by `omicGuiConfig("transcriptomics")$norm`.
- **What moved:**
  - The TMM/logCPM and RUV-grid option sections joined the panel library.
  - The frame registers its outputs by the omic's tabs and panels, and dispatches to trx's own run and reset handlers.
  - Trx's run is now core code in `utils_workflow_trx_norm.R`. It publishes through one in-process publication runner that every omic can use, `runWorkflowOmicsPublicationInProcess`.
  - `R/mod_trx_norm.R` is deleted. The lipid wrappers' `profile` argument now only selects the omic's configuration, and its signature and docs are unchanged.
- **Behaviour:** identical for every omic. No oracle changed, nothing was re-frozen, and no test was quarantined.
- **Code:** across `R/`, 439 lines added and 323 deleted against `d5244ffe` ([lines.tsv](lines.tsv)). The trx module file's 235 lines became 161 lines of core code plus the frame's gating and two panel sections.
- **Build:** `kit_110` = `core/browser-smoke.r5uUT0` (source `c1079dc7`; R/ and tests as `46ec9e01`).

## Work items

| Item | Commit(s) | What |
|---|---|---|
| W1 | `88b14014` | The trx tab's server captures, written first: in process and managed, with notification durations. 89 expectations, frozen on `d5244ffe`. |
| W2 | `2ea15095` | The trx normalisation core; the generic in-process publication runner; trx's server on the core through `notify`/`setStatus` callbacks. |
| W3 | `46ec9e01`, `23e75cf8` | Trx on the frame: the two library sections, the gated registrations, run and reset dispatch, trx's `norm` config and validator, `profile` reduced, `mod_trx_norm.R` deleted. |
| W4 | `ddcb7f42` onwards | Merge checks and the proof |

## Proof

| Check | Result |
|---|---|
| Merge checks ([merge-checks.txt](merge-checks.txt)) | `Collate` 361/361. NAMESPACE is byte-identical; `man/`, `app_ui.R`, `app_server.R`, `mod_omic_workflow.R` and `tools/` are unchanged. Back-references are still 6. The dependents grep finds no live caller of the 3 deleted names or of `mod_trx_norm.R`. |
| Run 1, targeted ([targeted.txt](targeted.txt)) | 45 files: 5,881 passed, 0 failed, 0 errors, under strict skips |
| Run 2, Simulation Setup ([entry-suite.md](entry-suite.md)) | trx_screen alone is identical. Traced calls: the frame 1, `runTrxNormalisation` 2, `runTrxNormInProcess` 2, the publications 1, the generic runner 8. These are exactly the ticket's expectations. |
| Run 3, entry suite ([entry-suite.md](entry-suite.md)) | trx_screen, lipid_canonical, metab_lc, metab_gc and metab_combined: all identical |
| Run 4, DIA session and memory set ([dia-memory-set.txt](dia-memory-set.txt)) | 733/2, 191/1, 38/1, 174/4, 80/11: each equal to its baseline |
| Run 5, `R CMD check` ([check.log](check.log)) | 3 WARNINGs, 3 NOTEs. The sections and the code NOTE are identical to KIT-105's. |
| Run 7, journeys on `kit_110` ([member-results.txt](member-results.txt)) | 110a to 110t: all accepted under `exactness_and_completion`. 110c ran trx's managed path through the frame. |
| Normalisation-phase memory ([norm-peaks.txt](norm-peaks.txt)) | Every phase within ±1.7% of KIT-105 (110c normalisation 0.575 GiB against 0.574) |
| Run 8, e2e lanes ([e2e-final.txt](e2e-final.txt)) | lipid 88, metab 190, trx 14, DIA 42, DIA_LIMPA 79, enrichment 98, TMT/LFQ 153: 0 failures and 0 skips, each equal to KIT-105 |
| Run 9, skip comparison ([skip-compare.md](skip-compare.md)) | In both modes, only the two new files differ from KIT-105. The trx MCI check-mode skips are as declared (3, 3, 4). |

## Findings, decisions and parked items

- **No product defect found.** See [findings.md](findings.md).
- **Instrumentation:**
  - The lifecycle harness's `ListReporter` version, written by the lead after KIT-105's lanes, was broken and had never been run. W2 repaired it (`be0a2601`), and `kit-105/entry-suite.md` now says so.
  - Run 9's first check-mode attempt exited before any test ran (a variable was not exported). It was fixed and rerun once.
- **Not verified:** under the in-process test worker, a managed trx run retried after a refused run in the same session fails with "native reference roots contain a live payload". W1 gives the refusal its own session. A real worker process was not checked, so whether the app shows this is unknown. No journey or e2e lane drives a refusal followed by a retry.
- **Recorded divergences:**
  - Trx's tab now writes the frame's two console start lines.
  - Its outputs and observers register in the frame's order.
  - Nothing captured changed.
- **Parked:**
  - P-21 (a lead decision, for review): trx keeps registering the outputs its tab does not show.
  - P-22: the harness step label `mod_trx_norm`.
- **Hand-offs:**
  - **KIT-109:** `profile` can be removed from the lipid wrappers' signatures. The impact-map rule for `mod_omic_norm*.R` and `utils_workflow_trx_norm.R`.
  - **KIT-111:** proteomics joins the frame. The gating by tabs and panels (`norm_by_tab`, `norm_by_panel`) and the handler path (`run$handler`) are the hooks it can use.

## Rollback

`git revert` of the KIT-110 commits, in reverse order. This restores `R/mod_trx_norm.R`, the `profile` branches and trx's `norm = NULL`. No oracle or quarantine needs restoring.
