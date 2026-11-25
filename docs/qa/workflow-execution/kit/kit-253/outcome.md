# KIT-253 G0: a Reset or Save Design pressed while the normalisation's figures run

Ticket: `tickets/kit/KIT-253.md` r1.2 (Re-check 1 of `AUDIT-KIT-253-plan-20261125`: may start). Commit under test:
v0.9-staging `16112917329a63ab062655f04db5bf1dd7df9729`. This branch (`kit-253-g0`) adds only this folder; no package
code.

## Declared before the run

- **Purpose.** Observe P-a (proteomics DIA, Reset to Pre-Normalization), P-c (proteomics DIA, Save Design) and L-a
  (lipidomics, Reset), each pressed while a named managed action runs, against the ticket's prediction (section 9);
  with two positive controls, P0 and L0 (the same journeys, no press), showing that the instrument sees the target
  completing and a publication.
- **Binding.** The package at 16112917 (the script commit differs from it only in this folder; the script checks
  this and stops otherwise). Driver: KIT-233's `journey.R` read with `git show origin/kit-233-journey:...` (its blob
  hash is recorded in `environment.tsv`). Plain launcher `packaging/electron/entry.R`, process launcher, managed
  worker, test mode, network namespace with loopback only. Fixtures `tests/testdata/e2e/prot_dia` and
  `tests/testdata/e2e/lipid_canonical` (the lipidomics browser lane's settings: ITSD off, normalisation none, RUV
  skipped; proteomics at the rendered defaults, RUV Automatic).
- **Ceiling.** One job, 90 job-minutes (script limit 75 minutes) = 6.0 core-hours; up to three VAL-54 instrument
  reruns, 18.0 core-hours; total 24.0. Expected 45 minutes = 3.0 core-hours.
- **What it decides.** All three press sequences as predicted with both controls proven: stop rule 1, the ticket
  closes with no code and the section 13 candidate row is written. Any press sequence not as predicted (with its press
  inside and its control proven): stop rule 2, stop and report. Target not seen, press not inside the action, or a
  control not proven: stop rule 3, an instrument defect, rerun under VAL-54.
- **Prediction and judging rules**, fixed in `g0.R` (`K253_SEQUENCES`, `K253_NOTICES`, `K253_RULES`) before any value
  is read, and written by the run to `prediction.tsv` and `prediction-rules.txt` before the browser phase.

## Runs

`runs.tsv`. One run, 45594223920, on commit 85df7ca1 (binds 16112917; `environment.tsv`: the diff is this folder
only, journey.R blob 7cc74bcf). Scripts job 27.4 minutes, whole run 28.2 job-minutes = 1.88 core-hours (expected 3.0,
ceiling 24.0). No instrument rerun. Returned tables: `g0-sequences.tsv` and `g0-45594223920/`.

## Observed (from the run's files)

Times are seconds after the press of the run button. "Head" is `<store>/dia-state-head.rds`'s entity id.

| Seq | Press | Target request appeared | Press | Target completion | Target status | Head moved by the press | Publication after the press | Completion notice |
|---|---|---|---|---|---|---|---|---|
| P0 (control) | none | 26.15 (`normalization-qc-figures`, after `normalization-ruv` complete) | - | 33.42 | complete | - | yes: 1 new action-output folder, footprint of 62 files (normalised matrices, post-normalisation and RUV figures) | yes |
| P-a | Reset to Pre-Normalization | 25.11 (after `normalization-ruv` complete at 25.11) | 25.11, inside (completion.rds absent before and after the click) | 31.55 | **stale** | yes, 26.33 (one move) | none: 0 files, 0 folders | no |
| P-c | Save Design | 26.73 (after `normalization-ruv` complete at 26.53) | 26.73, inside | 34.15 | **stale** | yes, 28.32 and 30.93 (the design save and its annotation commit) | none: 0 folders; 5 files written, all the design save's own (`design_matrix.tab`, `contrast_strings.tab`, `data_cln.tab`, `config.ini`, `manifest.json` in `scripts/proteomics`), none in P0's footprint | no |
| L0 (control) | none | 21.28 (`lipid_composite_figure`) | - | 37.15 | complete | - | yes: 1 new action-output folder, footprint of 8 files (`lipidomics_composite_QC_figure.pdf/png/rds` and their staging) | yes |
| L-a | Reset to Pre-Normalization | 22.49 | 22.49, inside | 38.23 | **stale** | yes, 23.57 (one move) | none: 0 files, 0 folders | **yes**, see below |

- In each press sequence the target's token `expected_head` equals the head just before the press, and the head at
  the end is the one written by the press (P-a, L-a; P-c: the last of its two moves). No request of the normalisation
  chain was dispatched after the stale target (`later_requests` empty). The tab status for normalisation stays
  `pending` in all five.
- Head state at the end (digest): P-a `protein_replicate_filtered` (the normalisation's source); P-c `raw_data_s4`
  (the re-saved design state, as KIT-247 reads it); L-a `lipid_reset`.
- Notices (`*-notices.tsv`): P-a and P-c: "Error applying proteomics step: managed action ended: stale" and "Error in
  normalization: managed action ended: stale" (the second at the same instant). P-a also shows KIT-225's warning
  "Normalisation reset: the current table is again the one the normalisation started from (Protein Replicate
  Filter)."; P-c "Design saved successfully. You can close this builder." and "Design matrix and contrasts saved
  successfully!". L-a: "Reset to pre-normalization state" (at 22.88), then at 38.28 "Error applying lipidomics step:
  managed action ended: stale" followed by "Normalization pipeline complete!": the contradictory message of the
  ticket's 2.4, observed.
- Settings: proteomics at the rendered defaults (RUV automatic, cyclicloess); lipidomics with the browser lane's
  settings (ITSD off, normalisation none, RUV skip; rendered defaults were ITSD on, cyclicloess, RUV automatic).

## Inferred

- The stale check at acceptance (`R/utils_workflow_execution.R:2058-2063`) ended the running action in all three
  press sequences and no publication followed, as the ticket read: no file, record or flag of the superseded table was
  written. The proteomics publication (`publishDIAReferenceNormalization`) did not run: no action-output folder, no
  footprint file and no completion notice after the press.

## Acceptance criteria

| # | Criterion | Met | Evidence |
|---|---|---|---|
| A1 | G0's table for P-a, P-c, L-a with the values of section 10 | met | `g0-sequences.tsv`; per sequence `g0-45594223920/<seq>-requests.tsv`, `-heads.tsv`, `-notices.tsv`, `-files-after.tsv`, `-publications.tsv`, `-digest.json`, `-end.png` |
| A2 | The decision applied and recorded | met | outcome (1), stop rule 1: the ticket closes with no code; the section 13 candidate row is handed to the coordinator (text in the implementer's report); no ledger row |

Verdict written by the script (`g0-45594223920/verdict.txt`): "G0: observed as predicted (stop rule 1)"; P0 and L0
"control proven"; P-a, P-c, L-a "as predicted".

## Not checked

- Metabolomics (L-a stands for it by the shared runner, not run); transcriptomics; a QC Revert (P-b); Export (P-d).
- A press during an earlier figures action of the lipid chain, or during the proteomics composite (only the first
  figures action after `normalization-ruv` was targeted).
- The tab's `normalization_complete` / `ruv_complete` values directly (not in the digest); read only through the
  completion notice, the publication's absence and the tab status.
- Real data, the installed bundle, Windows, user mode; the window's length on larger tables.
- Whether a benign in-process head write can make a running action spuriously stale (ticket section 13).

## For the change report

No behaviour changed. The window between the normalisation's state commit and its publication (ticket 2.1) is handled
by the existing stale check at acceptance (`R/utils_workflow_execution.R:1846-1849`, `:2058-2063` at 16112917):
observed on a runner (run 45594223920, prot_dia and lipid_canonical fixtures) for Reset and Save Design (proteomics)
and Reset (lipidomics); nothing of the superseded table is published. Users see internal codes ('managed action
ended: stale', KIT-252) and, on lipidomics (and by the shared runner metabolomics) when the press hits the composite
figure, a contradictory 'Normalization pipeline complete!' (`R/utils_workflow_omics_actions.R:429-432`,
`R/mod_omic_norm.R:440-444`): a message only, no result affected.
