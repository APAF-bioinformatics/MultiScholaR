# KIT-102 outcome: one shared shell and import for every omic (2026-11-08)

**Done.** The four omics now run on one shared app shell (`mod_omic_workflow.R`) and one import frame (`mod_omic_import.R`), each configured by its omic's entry in the validated contract (`omicGuiConfig`, `utils_omic_gui_config.R`).
- **Lipid, metab and trx** import through one shared feature-table panel (`mod_omic_import_feature_table.R`) over one import core (`utils_workflow_omics_import.R`).
- **Proteomics** imports through its panel on the same frame, over its core (`importProteomicsData`).
- **Deleted:** the lipid and metab import copies (9 files) and the trx import section.
- **Net change across `R/`:** 3,071 lines fewer (2,749 added, 5,820 deleted, against `c1041c56`).
- **Behaviour:** every omic is identical, apart from the one unification this ticket owned: **metab's import now opens the Design step**, as lipid's and trx's do.
- **Build:** `kit_102` = `core/browser-smoke.DGpv4v` (source `8bd0b4a1`).

## Work items

| Item | Commit(s) | What |
|---|---|---|
| W1 | `165ef779` | The configuration contract, its validator (`omic_gui_config_invalid`) and the DIA binding adaptor (`proteomicsDiaBinding`) |
| W2 | `d469cea3` | The feature-table import core (I1 to I13), the shared numeric-like conversion (`coerceOmicsNumericLikeColumns`), and trx's count check moved into the core |
| W3 | `0a3141c4`, `8a20c98f` | The proteomics import core. Back-references #1 and #2 moved into `utils_workflow_dia_sessions.R` (32 → 30). Spectronaut and PD-TMT reader parity. |
| W4 | `db6ffd27` | The shared shell (component registry, one status chain, the proteomics QC bypass) and the import frame, with byte-identical UIs for all four omics |
| W5 | `55671a5a`, `1c81c8c4` | The feature-table panel; lipid, metab and trx on it; the deletions; quarantines Q-25 to Q-30 |
| W6 | `d26aa08c` | Proteomics import on the frame. The process observer's order, including `releaseDIAWorkingReference`, is unchanged, and the UI md5 is unchanged. |
| W7 | `1b12a7bf`, `f8e4a942`, `8bd0b4a1` | The configs carry each omic's import preview, process and check (required); W5's stop-gap removed; the metab lanes re-frozen |

## Proof

| Check | Result |
|---|---|
| Merge checks ([w7/merge-checks.txt](w7/merge-checks.txt)) | `Collate` 381/381. The dependents grep covers 121 deleted names and 9 paths: the only hits are in quarantined files and `tools/refactor/`. Back-references 30. |
| Targeted tests ([w7/targeted.txt](w7/targeted.txt)) | 38 files, 4,980 passed, 0 failed, 0 errors, under strict skips |
| Entry suite ([entry-suite.md](entry-suite.md)) | The 7 non-metab lanes identical. The 3 metab lanes differ only in `tab_status$design_matrix` (`disabled` → `pending`) at `import` and `design_refusal`, re-frozen once; the freeze and verify runs are identical. |
| `R CMD check` ([check.log](check.log)) | 3 WARNINGs, 3 NOTEs. The code-problems NOTE is identical to KIT-101's. |
| Journeys on `kit_102` ([member-results.txt](member-results.txt)) | 102a lipid, 102b metab, 102c trx, 102m DIA_LIMPA, 102l LFQ, 102t TMT: all accepted and exact |
| Import-phase memory ([import-peaks.txt](import-peaks.txt)) | Every journey within 1% of KIT-101, so no assay copy was added |
| e2e lanes ([e2e-final.txt](e2e-final.txt)) | lipid 88, metab 190, trx 14, prot DIA 42, DIA_LIMPA 12 plus 1 skip by design, enrichment 98, TMT/LFQ 153: 0 failures, exactly as KIT-101 |
| Skip comparison, `load_all` ([skips-load_all.tsv](skips-load_all.tsv)) | 384 files, 33,602 passed, 0 failed. The only changes are the Q-25 to Q-30 quarantines (file skips and hidden blocks), Q-23, three new files, and the 2 known diagnostics. |
| Skip comparison, check layout ([skips-check.tsv](skips-check.tsv)) | 22,567 passed, 149 failed, 409 errors, against the baseline's 21,748, 152 and 410 (improved). File skips are new only where intended: F-12's two guards and Q-29. Five files each gained one error; each is a block added by a KIT-101 fix (F-08, F-09, F-11, F-15, F-27). It fails with the same source-tree read (`cannot open the connection`) as the existing blocks in those files under check. The aggregate's header step halted on `git rev-parse` inside the check directory; every per-file row was written. |

## Findings and parked items

- **F102-01:** metab's import validation warning is now computed after the numeric conversion. Only the warning text changes, and it never shows in the GUI. Accepted ([findings.md](findings.md)).
- **Divergences recorded by W5:** lipid's "sanitized" notice now appears only after a successful import, as metab's did; console messages are trimmed.
- **Parked from this ticket:** P-05 (non-numeric vendor tokens; refined by W3 for PD-TMT), P-06 (trx's `lipid_id_col` key), P-07 (the undefined `clean_acc` breaks the mixed-species organism filter).

## Rollback

`git revert` of the KIT-102 commits. The exported module functions keep their signatures throughout, and the quarantined tests are restorable from `quarantined-tests.md`.
