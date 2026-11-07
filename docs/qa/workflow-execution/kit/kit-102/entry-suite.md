# KIT-102: the entry-point suite (budget runs 2 and 3)

## Declaration (the ticket's experiment budget, runs 2 and 3)

- **Build.** Source tree `1b12a7bf` (origin/performance-v2 `1c81c8c4` plus the W7 wiring commit), `pkgload::load_all`
  (export_all, helpers) through `kit/scripts-078/run1.R`, one R process per test file; R 4.6.1;
  `R_LIBS` = `core/browser-smoke.rvol2q/library`; `MULTISCHOLAR_STRICT_SKIPS=1`; each run's records written with
  `MULTISCHOLAR_GUI_ENTRY_RECORD_DIR` to scratch; under `systemd-inhibit --what=idle:sleep`.
- **Run 2 (verify).** The four `test-gui-entry-*.R` files in parallel against the frozen KIT-101 oracles
  (`tests/testdata/gui-entry/*.json` at `1c81c8c4`). Ceilings 900 s per lane (asserted), 4,500 s per file.
  Expected: the 7 non-metab lanes identical; `metab_lc`, `metab_gc` and `metab_combined` differ only at the steps
  `import` and `design_refusal`, in `tab_status$design_matrix` (`disabled` to `pending`), and in the `workflow_progress`
  output. Any other difference is a defect.
- **Run 3 (re-freeze, once).** `test-gui-entry-metabolomics.R` with `MULTISCHOLAR_GUI_ENTRY_FREEZE=1` (writes the three
  metab oracles), then once more under strict skips against the new oracles. The two records must be identical. No other
  lane is re-frozen.
- Comparison: each record's steps against the oracle's (or the other record's) after the JSON round trip
  `guiEntryExpectLane` applies, path by path (scratch `compare.R`).

## Run 2: verify against the KIT-101 oracles (2026-11-08 03:11-03:24)

| File | Pass | Fail | Error | Skip | Elapsed |
|---|---|---|---|---|---|
| test-gui-entry-proteomics.R (5 lanes) | 187 | 0 | 0 | 0 | 715 s |
| test-gui-entry-lipidomics.R | 35 | 0 | 0 | 0 | 82 s |
| test-gui-entry-transcriptomics.R | 29 | 0 | 0 | 0 | 48 s |
| test-gui-entry-metabolomics.R (3 lanes) | 99 | 6 | 0 | 0 | 181 s |

Path-by-path against the oracles:

- `prot_dia`, `prot_dia_limpa`, `prot_tmt`, `prot_lfq`, `prot_lfq_fragpipe`, `lipid_canonical`, `trx_screen`:
  0 differing paths each (identical).
- `metab_lc`, `metab_gc`, `metab_combined`: 2 differing paths each, and the 6 failures are these 6 step comparisons:
  - `<lane>/import/tab_status/design_matrix`: oracle `"disabled"`, record `"pending"`;
  - `<lane>/design_refusal/tab_status/design_matrix`: oracle `"disabled"`, record `"pending"`.

Within the expected set; nothing else differs. The `workflow_progress` capture (the stepper's rendered UI hash,
`metabolomics_workflow-workflow_progress`) is identical to the oracle: `render_workflow_stepper`
(func_general_applets.R) shows a `disabled` and a `pending` step the same way, as W5 found (`w5/entry-suite.txt`).
The ticket allowed it to change; it did not, so it is not among the changed captures.

## Run 3: the metab re-freeze (2026-11-08 03:24-03:29)

- **Freeze run:** `test-gui-entry-metabolomics.R`, `MULTISCHOLAR_GUI_ENTRY_FREEZE=1`: pass 45, fail 0, error 0, skip 0,
  153 s. It rewrote `tests/testdata/gui-entry/metab_{lc,gc,combined}.json` (`source_sha` `1b12a7bf`, frozen
  2026-11-08).
- **Verify run** against the new oracles, strict skips: pass 105, fail 0, error 0, skip 0, 155 s.
- **Determinism:** the verify run's records against the freeze run's records: 0 differing paths in each of the three
  lanes; against the new oracles: 0. The freeze run's records against run 2's metab records: 0.

### Every changed capture, attributed to the metab-opens-Design unification (KIT-102 D2; test MCI-012.7)

| Oracle | Step | Capture | Before (KIT-101, `ddf42372`) | After (`1b12a7bf`) |
|---|---|---|---|---|
| metab_lc.json | import | `tab_status$design_matrix` | `disabled` | `pending` |
| metab_lc.json | design_refusal | `tab_status$design_matrix` | `disabled` | `pending` |
| metab_gc.json | import | `tab_status$design_matrix` | `disabled` | `pending` |
| metab_gc.json | design_refusal | `tab_status$design_matrix` | `disabled` | `pending` |
| metab_combined.json | import | `tab_status$design_matrix` | `disabled` | `pending` |
| metab_combined.json | design_refusal | `tab_status$design_matrix` | `disabled` | `pending` |

Besides these, only each oracle's header changed (`source_sha`, `frozen`). Every other capture of the three lanes
(heads, histories, state digests, artifacts and their sha256, notifications, modals and every output, including
`workflow_progress`) is unchanged. No other lane was re-frozen.
