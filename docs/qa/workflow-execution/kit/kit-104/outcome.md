# KIT-104 outcome: one shared QC shell and sub-tabs for every omic (2026-11-08)

**Done.** The four omics now run their QC tab on one shared shell (`mod_omic_qc.R`). Lipid, metab and trx also share one set of four sub-tabs (`mod_omic_qc_tabs.R`: intensity, duplicates, ITSD and finalise), all configured by each omic's `qc` group in `omicGuiConfig`.
- **Proteomics** keeps its peptide and protein orchestrators as group plugins on the same shell.
- **Code:** 3,173 lines fewer across `R/` (2,090 added, 5,263 deleted, against `b423b764`).
- **Parity change** (the user's 2026-11-06 decision): the QC-finalise design-alignment check now runs for every feature omic. It landed without changing any frozen lane.
- **P-16** (the lead's decision under the standing authority): DIA_LIMPA's QC completes at `peptide_count_filtered`. The prot_dia_limpa lane was re-frozen once for it, with attribution.
- **F-25** is fixed: the DIA_LIMPA e2e lane now runs to the end.
- **Build:** `kit_104` = `core/browser-smoke.xocrNB`.

## Work items

| Item | Commit(s) | What |
|---|---|---|
| W1 | `5f8cb531` | The `qc` group and its validation. Four tests written first: parity (a), UI identity (b), F-24's display-only cases (c), the ITSD capture (d). |
| W2 | `7b99e7e8` | The QC core for lipid, metab and trx (counts, duplicates, ITSD, finalise validation, `runWorkflowOmicsQcAction`). Back-reference #9 moved. F104-03 fixed. |
| W3 | `3cfad270` | Proteomics back-references #6 to #8 moved with their roxygen tags; NAMESPACE byte-identical. |
| W4 | `b7df8291`, `ce642947` | F-25: the DIA_LIMPA e2e lane drives the limpa QC, normalisation and native export and runs to the end (79/0/0/0). F104-02 is a test-only fix. |
| W5 | `cf539d67` | The shared shell and sub-tabs, with lipid on them (2,444 → 411). F-24's guard is on Detect duplicates and Analyze internal standards. |
| W6 | `b3d4a5b0` | Metab and trx on the shell (metab 2,667 → 199). The parity check is in place (metab D11, trx D10 on the F104-04 structural shell). The dead lipid profile code is removed. |
| W7 | `65010461` | Proteomics QC on the shell through wrappers; P-16. The stop rule was not needed. |
| W8 | `8b7f3eae`, `735ad485` | The finalise fields are required; clean-ups; merge checks; the prot_dia_limpa re-freeze; [parity.md](parity.md) |

## Proof

| Check | Result |
|---|---|
| Merge checks ([w8/merge-checks.txt](w8/merge-checks.txt)) | `Collate` 369/369. NAMESPACE is byte-identical (roxygen was not run). The dependents grep found no live caller of the 102 deleted names or the 7 paths. Back-references are 25. |
| Targeted tests ([w8/targeted.txt](w8/targeted.txt)) | 92 files, 8,980 passed, 0 failed, 0 errors, under strict skips |
| Entry suite ([w8/entry-suite.txt](w8/entry-suite.txt)) | 9 lanes identical. prot_dia_limpa was re-frozen once at `8b7f3eae` for P-16 (32 captures across 14 steps: QC and normalization status and the stepper; no state, file or notification). Only that oracle file changed. |
| DIA session and memory set ([w8/memory-session.txt](w8/memory-session.txt)) | 733/2, 191/1, 38/1, 174/4, 80/11: each equal to its baseline |
| `R CMD check` ([check.log](check.log)) | 3 WARNINGs, 3 NOTEs. The code NOTE is identical to KIT-103's. |
| Journeys on `kit_104` ([member-results.txt](member-results.txt)) | 104a lipid, 104b metab, 104c trx, 104m DIA_LIMPA, 104l LFQ, 104t TMT: all accepted and exact |
| QC-phase memory ([qc-peaks.txt](qc-peaks.txt)) | Every route within ±2.4% of KIT-103 |
| e2e lanes ([e2e-final.txt](e2e-final.txt)) | lipid 88, metab 190, trx 14, prot DIA 42, **DIA_LIMPA 79 (no skip; was 12 plus 1 skip)**, enrichment 98, TMT/LFQ 153: 0 failures |
| Skip comparison ([skips-load_all.tsv](skips-load_all.tsv), [skips-check.tsv](skips-check.tsv)) | Against KIT-103, every change is a recorded quarantine (Q-56 to Q-76, Q-80) or one of three new files. The check layout first showed 18 failures in the new files, from source-file checks the installed package cannot answer. Those checks now run under `load_all` only (`ca35c462`), and the new files' remaining check-mode errors are the known source-tree read. |

## Findings, decisions and parked items

- **F104-01** (queued to KIT-108): the native export summary shows "Report Template: NA" on every proteomics route. Display only.
- **F104-02** (fixed): stale DIA_LIMPA e2e study_parameters checks.
- **F104-03** (fixed): lipid's ITSD analysis errored on a NULL pattern (test sessions only).
- **F104-04** (lead decision, implemented): trx states cannot be read as a validated shell, so trx's finalise uses a structural shell. No class contract changed.
- **P-16** was fixed here by the lead's decision, and is listed in the parked file for the user's review. **P-17** (no proteomics QC-finalise) is parked.
- **F-20's** workaround (`mod_prot_norm_workflow_helpers.R:518-524`) is now dead in the app; it is left for KIT-105.
- **Divergences (logs and edge cases only):**
  - with standards found but no readable intensities, the ITSD intensity figure is cleared, as metab did;
  - the shell's start and completion log lines are new;
  - metab's zero-assay duplicate text branch is gone.
- **Hazard recorded:** never regenerate NAMESPACE with the installed roxygen2 8.0.0.

## Rollback

`git revert` of the KIT-104 commits; restore `tests/testdata/gui-entry/prot_dia_limpa.json` from `8b7f3eae^`. The exported module functions keep their signatures, and the quarantined tests are restorable from `quarantined-tests.md`.
