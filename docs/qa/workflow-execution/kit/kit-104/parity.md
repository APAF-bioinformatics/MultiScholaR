# KIT-104 parity record

Two behaviour changes in KIT-104 can change an entry-lane capture. This file records both, the cases each now refuses
or changes, and the lanes they touched.

## 1. The QC-finalise design-alignment check, for every feature omic (W6, `b3d4a5b0`)

**What changed.** Finalise now runs the omic's design-alignment check for lipid, metab and trx. It runs in
`validateWorkflowOmicsQcFinalize` after the existing checks (the class, the assays under `finalize_zero_rows`, the
design rows), through `omicGuiConfig(<omic>)$qc$finalize_alignment`. Since W8 that field is required for every omic
with a finalise sub-tab (lipid, metab, trx). Proteomics has no finalise step (P-17, kept as it is).

| Omic | Check | Before KIT-104 |
|---|---|---|
| Lipid | `omicsQcFinalizeAlignment(shell, "lipidomics")`, the shared alignment (`validateOmicsDesignAlignment` under `omicsDesignRules("lipidomics")`, the design's runs as the sample columns, `requireReplicates = FALSE`) | the same check (the control) |
| Metab (D11) | `omicsQcFinalizeAlignment(shell, "metabolomics")`, the same check with lipid's arguments | class, assays and design rows only |
| Trx (D10) | `trxDesignAlignment` (`validateTrxDesign` on the count table's sample columns, the design and `args$group_levels`), the check every other trx step applies | class, row counts and design rows only |

**Trx's shell (F104-04, the lead's option 1).** `TranscriptomicsAssayData`'s validity refuses a count table with no
rows, so the generic shell read (`readWorkflowOmicsShell`, `measures = "features"`) would refuse every trx state. Trx's
finalise reads a structural shell instead, through `qc$finalize_shell` (optional, trx only): the stored slots and
design, the count table as its zero-row prototype (`trxReadTemplate`, no validity run) and the stored row counts as
the attribute `assay_rows`. No class contract changed. The memory boundary holds: the finalise still reads no count
rows.

**The cases now refused.** A refusal commits nothing, leaves QC `pending` and shows
`Error finalizing QC: <text>`.
- **Metab:** a design with a blank group assignment. Text: `Current MetaboliteAssayData object has assay/design
  mismatch: Design matrix group column 'group' contains blank assignments.` Before, it committed
  `metab_qc_complete`. A design run absent from the assays was already refused, by the class's validity on the shell
  read, and still is.
- **Trx:** a design run absent from the count table (for example `S_absent`), and a group outside the declared
  levels. Before, both committed `trx_qc_complete`.
- **Lipid:** unchanged (it refused both before).

Tests: `test-omic-qc-module.R` (a), failing first for metab and trx (`w6/test-first.txt`), passing after
(`w6/targeted.txt`, `w8/targeted.txt`).

**Lanes and corpora.** No lane changed because of this check, and no design that a lane or fixture declares valid was
refused, so the park rule did not fire.
- Entry lanes: `trx_screen`, `lipid_canonical`, `metab_lc`, `metab_gc` and `metab_combined` are identical to their
  oracles after W6 (`w6/entry-suite.txt`) and at W8 (`w8/entry-suite.txt`). The proteomics lanes have no finalise
  step.
- Corpus and module-CI fixtures that run the finalise all finalise as before: MCI-014.5, the metab fixtures of
  `test-omic-qc-module.R`, the KIT-082 application-reads journeys, the `trx_screen` corpus (MCI-028.1 and .3, the
  W6 F104-04 block) and the KIT-040 `trx_edge_12` reference design (MCI-028.2).
- The journeys (104a/b/c) and the e2e lanes are recorded at runs 6 and 7.

Nothing was re-frozen for this change.

## 2. P-16: DIA_LIMPA's QC completes at `peptide_count_filtered` (W7, `65010461`; the lead's decision of 2026-11-08)

**What changed.** `omicGuiConfig("proteomics")$qc$complete_state` gains `DIA_LIMPA = "peptide_count_filtered"`. The
shared shell completes the route's QC tab when that state is committed, so normalisation opens as on the other routes.
Before, the DIA_LIMPA route left `quality_control` `pending` and `normalization` `disabled` until the export. This is
GUI status only: no state, history, provenance, artifact, notification, refusal or computed value changes.

**F-20.** The workaround in `mod_prot_norm_workflow_helpers.R:518-524` is dead in the app: when the normalisation tab
opens on `peptide_count_filtered`, `quality_control` is already `complete`, so the route-specific
`limpa_qc_complete` term decides nothing. Only its unit test (`test-prot-05b-norm-module-contracts.R`, "KIT-101
F-20") still reaches it, with a pending status built by hand. It stays until KIT-105, as the ticket says.

**The re-freeze of `prot_dia_limpa`.** The lane differed from its oracle (ddf42372) only in these captures (`w7/p16-diff.txt`,
reproduced at W8 in `w8/entry-suite.txt`). 32 captures across 14 steps, from `qc_protein_peptide_count` on:
- `tab_status$quality_control`, `"pending"` to `"complete"`: all 14 steps (`qc_protein_peptide_count`,
  `da_refusal`, `norm_refusal`, `norm`, `correlation`, `ruv_covariates`, `session_export`, `da_load`, `da`,
  `summary_refusal`, `summary_save`, `summary_copy`, `summary_report`, `summary_export`).
- `tab_status$normalization`, `"disabled"` to `"pending"`: 4 steps (`qc_protein_peptide_count`, `da_refusal`,
  `norm_refusal`, `norm`). From `correlation` on, the normalisation tab sets its own status, as before.
- `outputs$proteomics_workflow-workflow_progress`, the stepper that renders `tab_status`: all 14 steps (the UI
  digest changes with the two statuses above).

The 10 steps from `import_refusal` to `qc_protein_normalisation` are identical, and so is the step order.

**The re-freeze (W8, once).** Run 2 of the budget reproduced the difference exactly: the `prot_dia_limpa` record against
the ddf42372 oracle is line for line the same diff as `w7/p16-diff.txt`. The other nine lanes were identical. Only
that lane was then re-frozen, at `8b7f3eae` (the W8 code commit). The freeze ran
the one `test_that` block (`testthat::test_file(desc = "KIT-101 entry lane prot_dia_limpa: ...")`) under
`MULTISCHOLAR_GUI_ENTRY_FREEZE=1`, load_all, R 4.6.1, so no other oracle file was rewritten. The sha256 of the other
nine oracle files is unchanged. The new oracle against the old one differs in exactly the 32 captures above and nothing
else, including its step order. The run 2 record is identical to the new oracle (the lane is deterministic). The
verify run against the new oracles, under strict skips, is in `w8/entry-suite.txt`.

No other lane changed because of P-16: `prot_dia`, `prot_tmt`, `prot_lfq` and `prot_lfq_fragpipe` are identical to
their oracles (`w7/entry-suite.txt`, `w8/entry-suite.txt`), and their oracle files are unchanged.
