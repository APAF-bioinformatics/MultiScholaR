# KIT-107 outcome: the enrichment module as a shared, omic-agnostic module; F-16 (2026-11-09)

**Done.** Enrichment runs on one shared frame (`mod_omic_enrich.R`, `mod_omic_enrich_panels.R`), configured by an `enrichment` group in `omicGuiConfig`.
- **Who uses it:** proteomics is the only consumer; the group is NULL for lipid, metab and trx. Another omic is added later as configuration for g:Profiler organisms, or with a new kernel for lipid and metab compound enrichment (D13).
- **The core:** the computation moved out of the module into `utils_workflow_enrichment.R` and `utils_workflow_dia_enrichment.R`, behind two entry points, `runWorkflowEnrichment` and `restoreWorkflowEnrichment`. The kernels are unchanged, and nothing is renamed.
- **F-16** (the user, 2026-11-07): enrichment now records its 9 declared parameters as step provenance. The analysis record is complete, and replay supports the step. No captured value changed.
- **Four enrichment defects fixed**, none changing a method or result:
  - **F107-05, the important one:** a live enrichment run never saved its result. The kernel's static plots captured their working frame, and the store refused the lazy bindings in it.
  - **F107-01:** a second Run on the same tab failed.
  - **F107-02:** results were never restored after an Open.
  - **F107-03:** a missing DA table still produced a product.
- **The legacy non-manifest branch is retired.** The app could not reach it.
- **Code:** the proteomics enrichment module went from 4,882 lines to 47, with 1,393 lines of frame and 1,787 of core. `R/` is 1,509 lines fewer (2,274 added, 3,783 deleted, against `8e7294b6`; [lines.tsv](lines.tsv)).
- **Build:** `kit_107` = `core/browser-smoke.66sPmE` (`8b92d7ac`).

## Work items

| Item | Commit(s) | What |
|---|---|---|
| W1 | `eed336ea` | The `enrichment` group and validator. `test-omic-enrich-module.R`: 947 captures frozen on the start commit (both backends, the live kernels offline, the organism filter, refusals, restore, prefill). F107-01/02 confirmed by probe. |
| W2 | `037698f7`, `cc65db89`, `116b4c6d`, `eacd4a5b` | The core moves; F-16's provenance and replay; the equivalence of the core and the old run, which held. F107-04 fixed. The lead wired the config. |
| W3 | `951682e2`, `e8a8b2fc`, `6f6ed201` | The frame and panels with proteomics on them; the wrappers; 8 helper files deleted; D7 (the legacy retirement); F107-01/02/03/05 fixed. W1 (b) re-frozen once, with attribution. |
| W4 | `2b1e02c8`, `ade3c654`, `aeb1d3fd` | Required config fields; D6 test 2 in the `prot_dia` lane; quarantine rows Q-163 to Q-166; merge checks; the lane checks |

## Proof

| Check | Result |
|---|---|
| Merge checks ([merge-checks.txt](merge-checks.txt)) | `Collate` 349/349. NAMESPACE byte-identical; `man/`, the app and shell files, `mod_proteomics.R`, the kernels and `tools/` unchanged. Back-references 2. The dependents grep found 0 live hits. |
| D6 test 2 and the Simulation Setup ([entry-suite.md](entry-suite.md)) | The F-16 lane check fails on the start tree and passes on the merged tree. `prot_dia` alone: identical, with every traced function called (`restoreWorkflowEnrichment` 1 rather than the ticket's 0, because of F107-02's guard; it returns NULL). |
| F-16's record diff ([lane-changes.md](lane-changes.md)) | No step value differs, so no oracle was re-frozen |
| Run 1, targeted ([targeted.txt](targeted.txt)) | 20 files, 7,476 passed, 0 failed, 0 errors |
| Run 3, entry suite | All 10 lanes identical: trx 29, lipid 35, metab 105, proteomics 188 (187 plus the F-16 check) |
| Run 4, DIA session and memory set ([dia-memory-set.txt](dia-memory-set.txt)) | Each equal to its baseline |
| Run 5, `R CMD check` ([check.log](check.log)) | 3 WARNINGs, 3 NOTEs; sections and the code NOTE identical to KIT-106 |
| Run 7, journeys ([member-results.txt](member-results.txt), [dispatches.txt](dispatches.txt)) | 107a to 107t all accepted under `exactness_and_completion`, with dispatches identical to KIT-106. None reaches enrichment. |
| Run 8, e2e ([e2e-final.txt](e2e-final.txt)) | lipid 88, metab 190, trx 17, DIA 42, DIA_LIMPA 79, enrichment 98, TMT/LFQ 153: 0 failures, 0 skips |
| Run 9, skip comparison ([skip-compare.md](skip-compare.md)) | Every change is a recorded quarantine, a new file, or (check mode) the known source-tree read |

## Findings, decisions and parked items

- **F107-01, 02, 03 and 05:** fixed, each with a test that fails first ([findings.md](findings.md)). F107-05's triage is in [w3/f107-05-triage.md](w3/f107-05-triage.md).
- **F107-04:** fixed in W2. Saved tables are not content-addressed, so replay compares them by content.
- **Lead decisions:**
  - F107-03's refusal;
  - the W1 (b) re-freeze of the values the four fixes move, including F107-01's held-manifest value, which was not pre-authorised and was accepted as the fix's direct effect ([w3/w1-refreeze.md](w3/w1-refreeze.md));
  - the required config fields.
- **Divergences:**
  - a blank taxon now stops silently (before, it wrote to the bus and showed an empty error);
  - the console-only log lines;
  - live enrichment states now carry the plots' captured frames (size: KIT-113's memory note).
- **Parked from this ticket:** P-28 (declare the organism filter; decide with P-07) and P-29 (the step-map row). P-07's `clean_acc` failure stays frozen as today's behaviour.
- **Hand-offs:**
  - **KIT-113:** the enrichment run to the worker, measured; the state size with the plots' frames.
  - **KIT-109:** the impact-map rules for `mod_omic_enrich*.R`.
  - **KIT-063:** the multi-omics enrichment helpers are untouched.

## Rollback

`git revert` of the KIT-107 commits, in reverse order. Remove the Q-163 onwards skip lines. No entry oracle needs restoring.
