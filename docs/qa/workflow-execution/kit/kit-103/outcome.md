# KIT-103 outcome: one shared design tab and builder for every omic (2026-11-08)

**Done.** The four omics now share one design tab (`mod_omic_design.R`) and one design builder (`mod_omic_design_builder.R`), configured by each omic's `design` group in `omicGuiConfig`.
- **Lipid, metab and trx** go through one design core (`func_omics_design_shared.R`, `utils_workflow_omics_design.R`).
- **Proteomics** goes through the same builder server by way of its `dia_long` adapter, with its own builder UI kept as the exported plugin.
- **Code:** 4,532 lines fewer across `R/` (3,211 added, 7,743 deleted, against `3c657c14`).
- **Behaviour:** identical for every omic. No entry lane changed and nothing was re-frozen.
- **Defect fixed:** the lipid and metab "Import Existing Design" failed in the app (F103-01).
- **Build:** `kit_103` = `core/browser-smoke.yODoyT` (source `63d71960`).

## Work items

| Item | Commit(s) | What |
|---|---|---|
| W1 | `e82f72a9` | The `design` group of the contract, its validation, and the DIA adaptor's design entries |
| W2 | `fe1a9742` | The design core for lipid, metab and trx (C1 to C5, C7, C8, C10). F103-01 fixed: vroom's live parser pointer was rejected by the owner after a partial bus write. |
| W3 | `9b8b7d62`, `a22a5e94` | The proteomics design core. Back-reference #3 moved first (31 → 29); the UniProt and FASTA annotation step became core functions. |
| W4 | `5d077268` | The shared design tab and builder, with lipid and trx on them. Seven lipid design files deleted (3,309 → 115 lines). |
| W5 | `c56d1c64` | Metab and trx on the shared modules. The 866-line metab builder closure and the metab import and validation files deleted (2,694 → 166). F-01 retargeted with identical expectations. |
| W6 | `adad5fc0` | Proteomics on the shared builder server through the `dia_long` adapter (3,136 → 1,800). The stop rule was not needed. |
| W7 | `f0e90318`, `f0c36114`, `63d71960` | The config's design functions wired and required; the stop-gaps and `runLipidDesignBuilderObserverShell` deleted; Q-46 and Q-47; [run-order.md](run-order.md); merge checks |

## Proof

| Check | Result |
|---|---|
| Merge checks ([w7/merge-checks.txt](w7/merge-checks.txt)) | `Collate` 374/374. The dependents grep found no live caller of the 137 deleted names or the 11 paths. Back-references are 29. |
| Targeted tests ([w7/targeted.txt](w7/targeted.txt)) | 46 files, 4,459 passed, 0 failed, 0 errors, under strict skips |
| Entry suite ([w7/entry-suite.txt](w7/entry-suite.txt)) | All 10 lanes identical, nothing re-frozen |
| DIA session and memory set ([w7/memory-session.txt](w7/memory-session.txt)) | native-session 733/2, perf107-audit 191/1, perf107-checkpoints 38/1, checkpoint-dia-state 174/4, memory-lifetimes 80/11: each equal to its baseline |
| `R CMD check` ([check.log](check.log)) | 3 WARNINGs, 3 NOTEs. The code NOTE is identical to KIT-102's. |
| Journeys on `kit_103` ([member-results.txt](member-results.txt)) | 103a lipid, 103b metab, 103c trx, 103m DIA_LIMPA, 103l LFQ, 103t TMT: all accepted and exact |
| Import and design phase memory ([design-peaks.txt](design-peaks.txt)) | Every route within ±1.2% of KIT-102 |
| e2e lanes ([e2e-final.txt](e2e-final.txt)) | lipid 88, metab 190, trx 14, prot DIA 42, DIA_LIMPA 12 plus 1 skip by design, enrichment 98, TMT/LFQ 153: 0 failures, as KIT-102 |
| Skip comparison ([skips-load_all.tsv](skips-load_all.tsv), [skips-check.tsv](skips-check.tsv)) | Against KIT-102, every change is a recorded quarantine (Q-31 to Q-34, Q-36, Q-46, Q-47) or one of three new files. 0 new failures in either mode; check mode has 2 fewer errors. |

## Findings, divergences and parked items

- **F103-01** (fixed): see W2.
- **F103-02** (accepted): metab's design import now reads everything before writing the bus, so an error path no longer leaves half-written fields.
- **F103-03** (parked as P-15): the standalone design applet cannot save without a DIA owner. This has been true since KIT-048. The applet acceptance item "returns a save result" is unmet as written.
- **Kept per omic as config flags:** lipid records the design file on import and metab does not (P-14); metab checks the import files before closing the modal (P-11); proteomics keeps its integer replicates.
- **Console messages** from the old module servers were dropped.
- **Parked from this ticket:** P-08 to P-15.

## Rollback

`git revert` of the KIT-103 commits. The exported module functions keep their signatures, and the quarantined tests are restorable from `quarantined-tests.md`.
