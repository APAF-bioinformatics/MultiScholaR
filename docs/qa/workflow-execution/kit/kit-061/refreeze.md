# KIT-061 merge gate: re-freeze of the five proteomics GUI-entry oracles

Date: 2026-11-17. Branch `kit-061-merge` (from `performance-v2` at `26751813`, `kit-061` at `a2ccc5dc` merged
without fast-forward in `5d381272`). The new oracles record `source_sha` `5d38127250d0122ebbe09d52b0477295c870112c`.

Cause: KIT-061 replaces the proteomics integration archive `integration/proteomics_<lane>_final.msdia.tar` with the
integration input at `integration/proteomics/` (`integration_input.json`, `matrix_proteomics.tsv`, `samples.tsv`,
`features_proteomics.tsv`).

## Method

The existing freeze mechanism (`tests/testthat/helper-gui-entry.R`, `guiEntryExpectLane`):
`MULTISCHOLAR_GUI_ENTRY_FREEZE=1` with `MULTISCHOLAR_GUI_ENTRY_ORACLE_DIR` pointed at a scratch folder, running
`test-gui-entry-proteomics.R` in one process (`pkgload::load_all`, proof library with fastduckdb 0.2.0, `nice -n 19`,
BLAS/OpenMP threads 1). Each scratch oracle was compared with the committed one leaf by leaf (JSON flattened to paths)
before it replaced the committed file.

Attribution rule: only these rows may move: the integration file rows in the four Summary steps (the `.msdia.tar`
row removed, the integration input's files added), the Summary text line `Integration Object: ...`, and the
notice `Saved Integration artifact` -> `Saved Integration input`. Any other difference stops the re-freeze.

## Result per oracle

All five: 22 moved leaves, all attributed, 0 unattributed. Header fields: `source_sha` and `frozen` only
(`r_version`, `load_mode`, `lane` unchanged). The 22 leaves are the same in every lane:

- Steps `summary_copy`, `summary_export`, `summary_report`, `summary_save` (4 steps), `artifacts`:
  `integration/proteomics_<lane>_final.msdia.tar` `<present>` -> absent, and 4 rows added:
  `integration/proteomics/integration_input.json` (`<present>`), `matrix_proteomics.tsv`, `samples.tsv`,
  `features_proteomics.tsv` (sha256, identical across the 4 steps). 4 x 5 = 20 leaves.
- `summary_save` / `notifications` / 0 / `message`: `Saved Integration artifact` -> `Saved Integration input`.
- `summary_save` / `outputs` / `proteomics_workflow-session_summary-session_summary`: one line only,
  `Integration Object: proteomics_<lane>_final.msdia.tar` -> `Integration Object: proteomics/integration_input.json`.

| Oracle | Previous `source_sha` (frozen) | Leaves moved | Unattributed | `matrix_proteomics.tsv` | `samples.tsv` | `features_proteomics.tsv` |
|---|---|---|---|---|---|---|
| `prot_dia` | `bdbc0286` (2026-11-15) | 22 | 0 | `350b5642...` | `a352e685...` | `4318f385...` |
| `prot_dia_limpa` | `5e33c65f` (2026-11-16) | 22 | 0 | `23eb5201...` | `31bee7fa...` | `7a9d96b6...` |
| `prot_lfq` | `dc2a05b2` (2026-11-12) | 22 | 0 | `4b0e04fc...` | `2ccfb360...` | `d6a9a448...` |
| `prot_lfq_fragpipe` | `dc2a05b2` (2026-11-12) | 22 | 0 | `6341bdce...` | `2ccfb360...` | `d6a9a448...` |
| `prot_tmt` | `dc2a05b2` (2026-11-12) | 22 | 0 | `78e411e2...` | `45cdd557...` | `d6a9a448...` |

The LFQ, FragPipe and TMT fixtures share one protein set (same `features_proteomics.tsv`); LFQ and FragPipe share
one design (same `samples.tsv`). Their matrices differ.

## Tests on `kit-061-merge` with the re-frozen oracles

One file per process, in order, `nice -n 19`, threads 1, proof library (fastduckdb 0.2.0), Tctl 65 C throughout:

| File | Passed | Failed | Skipped | Seconds |
|---|---|---|---|---|
| `test-workflow-integration-input.R` | 76 | 0 | 0 | 100 |
| `test-omics-summary-core.R` | 701 | 0 | 2 (quarantines Q-175, Q-183) | 521 |
| `test-omic-summary-module.R` | 1703 | 0 | 0 | 559 |
| `test-workflow-native-session.R` | 734 | 0 | 2 (owner/proof-root guards) | 1186 |
| `test-omic-gui-config.R` | 1906 | 0 | 0 | 75 |
| `test-omic-norm-module.R` | 639 | 0 | 0 | 547 |
| `test-gui-entry-proteomics.R` | 193 | 0 | 0 | 2309 |

## Audit fixes K061-1..6 (2026-11-18): no oracle movement

After the independent audit's fixes (the final state and the role from the head's own lineage, the backup-rename-restore
publish, the removal of a stale input, the accession's first member without its isoform suffix, the reader's checks),
`test-gui-entry-proteomics.R` passes against the oracles above unchanged (193/0, 2519 s): no `features_proteomics.tsv`
row moved (no fixture id carries `;` or a hyphen), so no re-freeze was needed. The seven files above, rerun one per
process on the final code: integration input 118/0, summary core 701/0 (2 skipped), summary module 1703/0, native
session 734/0 (2 skipped), gui config 1906/0, norm module 639/0, gui-entry proteomics 193/0.

## Notes

- Freeze run: 76 passed, 1 failed, 2888 s. The failure is the lane wall-time ceiling (`record$seconds < 900`) on
  `prot_dia`: 1209.8 s, under a 12-core simulation load with the R process in the idle scheduling class. It is a
  timing assertion, not oracle content.
- Not changed (harness revision, with KIT-129): `tools/profiling/run_workflow_checkpoint_proof.R:1905` still expects
  `integration/proteomics_<label>_final.msdia.tar` (`native_first`) or `..._final_s4.RDS`.
