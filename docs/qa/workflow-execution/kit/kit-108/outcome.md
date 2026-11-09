# KIT-108 outcome: summary and report on one shared frame; F104-01; the report-step memory (2026-11-10)

**Done.** All four omics now run Summary and Report on one shared frame (`mod_omic_summary.R`), configured by a `summary` group in `omicGuiConfig`. The work runs through one core (`utils_workflow_summary.R`), with proteomics' route pieces in `utils_workflow_dia_summary.R`.
- **Transcriptomics' RUVIII export** is a configured export step, run managed and in process through `runWorkflowSummaryExportStep`. Its acceptance, record and figures are unchanged.
- **The six exported summary wrappers** keep their names, signatures and roxygen blocks. NAMESPACE and `man/` are byte-identical.
- **Ten helper files deleted.** The summary set went from 4,348 lines to 1,696; `R/` is 2,349 lines fewer (1,915 added, 4,264 deleted, against `5be34ee5`; [lines.tsv](lines.tsv)).
- **Two defects fixed, neither changing a method or a result:**
  - **F104-01:** every proteomics route's export summary said "Report Template: NA". It now names the template the report step uses.
  - **F108-03:** trx showed "Saved Integration S4 Object" after a save that wrote nothing.
- **The report-step memory** (the user's question, 2026-11-07) is diagnosed and partly fixed:
  - **Cause:** "Copy to publication" re-reads each large DA workbook with openxlsx (about 2.1 GB of transient heap). R frees it, but glibc keeps the pages.
  - **Route A** (build the workbook one sheet at a time) cut the application's report-step peak by 0.52 GiB (metab) and 0.27 GiB (lipid), and its level after the report by 0.60 and 0.36 GiB. Content is identical.
  - **D8's threshold is not met:** the peak misses it by 0.68 and 0.79 GiB, the teardown level by 0.57 and 0.62 GiB. Closing the rest needs P-30 (route B, `malloc_trim`; or route C, the copy on the worker), which is the user's decision.
- **Build:** `kit_108` = `core/browser-smoke.X6BzNV` (`7b52f535`).

## Work items

| Item | Commit(s) | What |
|---|---|---|
| W1 | `28d1e29a`, `031a6ca0`, `0ce7e60e` | The `summary` group and validator. The summary captures frozen on the start commit. Probes; F108-01/02/03 recorded. |
| W2 | `308dadad`, `5e9a4b65`, `b529f4cf`, `fa4e6181` (members `24853737`, `b97ad9e8`, `f4e0fe2d`) | The memory instruments I1 to I3 and their lifecycle proof (M1); the baselines `108p`/`108q` (M2); the replica (M3); the diagnosis ([memory/report-memory-diagnosis.md](memory/report-memory-diagnosis.md)) |
| W3 | `c2702c21`, `1b3b9205`, `00a8c712`, `d92fa6d7`, `91192835`, `1cb2a79b` | The core and the frame. Lipid, metab and trx on them; six helper files deleted. F108-03 fixed, with one attributed W1 re-freeze ([w3/w1-refreeze.md](w3/w1-refreeze.md)). Q-170 to Q-175. |
| W4 | `86a5c18b`, `f6fbbb68`, `30c769a8`, `1eead052` | Proteomics on the frame; four helper files deleted. F104-01 fixed (test 1 fails first). D11's retargets. Q-180 to Q-183. |
| W5 | `f5131643`, `98f89f8f`, `4bc2b239` | Route A, with its test (fails at W4, passes at W5). The replica at the W4 and W5 trees; content identical ([memory/m4-replica.md](memory/m4-replica.md)). P-30 raised. |
| W6 | `9a50afcf`, `14b5b20c`, `58e03b8a`, `f718c5ea`, `7b52f535`, then the run commits | The lane and e2e F104-01 checks; the required function fields; merge checks; the proof |

## Proof

| Check | Result |
|---|---|
| Merge checks ([merge-checks.txt](merge-checks.txt)) | `Collate` 342/342. NAMESPACE byte-identical. `man/`, `inst/`, and the app and shell files unchanged. `tools/` is D11's three lines. Back-references 2. The dependents grep found 1 live hit: a string in a refusal test, not a call. |
| D5 test 2 and the Simulation Setup ([entry-suite.md](entry-suite.md)) | F104-01's lane check fails on the start tree ("NA" vs "DIANN_report") and passes on the merged tree. `trx_screen` and `prot_dia` alone: identical, with the expected traces. |
| D5's lane diffs ([lane-changes.md](lane-changes.md)) | None; no oracle was re-frozen |
| Run 1, targeted ([targeted.txt](targeted.txt)) | 31 files, 9,034 passed, 0 failed, 0 errors |
| Run 3, entry suite | All 10 lanes identical: trx 29, lipid 35, metab 105, proteomics 193 (188 plus F104-01's five) |
| Run 4, DIA session and memory set ([dia-memory-set.txt](dia-memory-set.txt)) | Each equal to its baseline |
| Run 5, `R CMD check` ([check.log](check.log)) | 3 WARNINGs, 3 NOTEs; every section identical to KIT-107's, apart from paths and time |
| Run 7, journeys ([member-results.txt](member-results.txt)) | 108a, b, c, m, l, t all accepted under `exactness_and_completion`, with dispatches identical to KIT-107. Each runs the shared report step. |
| Run 8, e2e ([e2e-final.txt](e2e-final.txt)) | lipid 88, metab 190, trx 17, DIA 42, DIA_LIMPA 80 (79 plus F104-01's check), enrichment 98, TMT/LFQ 153: 0 failures, 0 skips |
| Run 9, skip comparison ([skip-compare.md](skip-compare.md)) | Every change is a recorded quarantine, a new file, or a check-mode-only skip or source-tree read. The D8 test's check-mode error was a defect in the test, fixed as the run's correction pass. |
| M4, memory ([memory-outcome.md](memory-outcome.md)) | `108v`/`108w` valid, exact, accepted, with their baselines' dispatches. Peak and teardown lower, but D8 not met. No slice or worker peak worse by more than 0.1 GiB. |

## Findings, decisions and parked items

- **Fixed:** F104-01 and F108-03 (each with a test that fails first), and the run-9 D8 test defect ([findings.md](findings.md)).
- **To KIT-109:**
  - **F108-01:** `copyToResultsSummary`'s failure records are lost inside closures.
  - **F108-04:** the core test's `study_parameters` capture depends on whether it runs in a worktree. KIT-109's full suite runs that file, so this needs settling there.
  - **Dead code:** `trxSummaryIntegrationReady` and `buildLipidSummarySessionState` have no callers.
  - **Profile hooks:** `trxOmicProfile()`'s remaining base fields and the `profile` argument of `mod_lipid_summary_*`.
  - **Impact-map rules** for the three new files.
  - **P-31.**
- **Parked:**
  - **P-30:** the memory routes B and C, updated with M4's numbers.
  - **P-31:** the branch map.
- **Divergences, accepted:**
  - W1's three report-flow logicals and proteomics' detection text, kept as today's behaviour;
  - the proteomics integration notice returned by the route writer.
- **Process, recorded:**
  - W3's host-rule breach during `108p`; the attribution stands;
  - the lane check's first version looked in the wrong folder (an instrumentation failure, rerun once);
  - run 9's two rerun slips (wrong library, then wrong checkout), which led to F108-04.

## Rollback

`git revert` of the KIT-108 commits, in reverse order. Then:
- remove the Q-170 onwards skip lines;
- regenerate the bindings after reverting the generator block (`b495ae89`).

No entry oracle needs restoring. Route A is its own commit (`f5131643`) and reverts alone.
