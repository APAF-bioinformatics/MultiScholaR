# KIT-220 outcome (r1.3; stopped before run D, control C, the re-freeze and R3, as ordered)

Branch `kit-220` from `origin/v0.9-staging` d048898d (before KIT-204, at the coordinator's order). Commits: 0804ed93
tests first; bb092685 and 19b9a4dc instrument fixes 1 and 2 (VAL-54; tests and G0 only, no expected value changed);
c788f04a the fix; c69e4df0 the r1.3 test first (I05); d0b8672c the r1.3 repair (lipid/metab up/down, I01; E-7 note, I04);
then the records and the I04 rework (below). Runs: `runs.tsv`. About 53.0 core-hours used (27.1 to R2, 25.9 for R2b).

## Observed

**G0** (`g0.tsv` base at 19b9a4dc, `g0-head.tsv` on the fix; `g0-w1.tsv`, `g0-versions.tsv`, `g0-feature*.tsv`).
- Synthetic proteomics fixture (400 proteins): 17 rows with 0.01 <= q < 0.05, 10 with 0.05 <= q < 0.10. Base: the
  first count table, the published table and the report sentence say 65 up / 62 down at every threshold (q < 0.05);
  the tab says 54/56 at 0.01 and 72/65 at 0.10. Head: all four agree at each threshold (54/56; 72/65; 30/30 with TREAT
  0.5), and the table carries `q_threshold`.
- trx_edge_12 (Sera_vs_RPMI): 155 genes with 0.01 <= FDR < 0.05. Base: published volcano line at -log10 0.05 and bar
  320/299 at every threshold; tab 229/235 at 0.01. Head: line at 2.0 (0.01) and 1.0 (0.10), bar = tab.
- KIT-112 W1 fixture at q 0.01, TREAT 0.5: 3 rows per contrast between 0.01 and 0.05 (so `dia_qvalue_only` moves).
- edgeR 4.10.1 has no `decideTestsDGE` (ticket E-3's name); `decideTests.DGEExact` is exported and is Glimma 2.15.2's
  documented default. The fix uses `edgeR::decideTests.DGEExact(qlf, adjust.method = "BH", p.value = t)`.

**R1** (45501071657, base code): every failing-first row failed at its stated base value (`r1-failures.txt`): T-P1
(65/62 against 54/56 and 72/65), T-P2, T-P3 (no column; 65/62; no phrase), T-PAR prot and trx, T-P5 (rows with
|log2FC| < 1 uncoloured), T-T1 (line 1.30, not 2.00), T-T2 (no `status`), G5, E-7 rows. Stop rules 1 and 2: not fired.

**R2** (45503395846, c788f04a against 19b9a4dc, 64 files, check): 60 files `same`; check findings identical (6 on both).
- `test-kit220-da-threshold-consumers.R` `fixed`: 73 passed, 4 failed on head, all four the lipidomics/metabolomics
  T-PAR rows, failing identically on base (finding F6 below).
- `test-omic-da-module.R`, `test-omic-enrich-module.R` (`new-failure`, `message-changed`): the only moved value is the
  `files` leaf `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` (section 6, RC1-01;
  `8a7cae29…` -> `2477b5a9…`, `34e7cdf6…` -> `a5562559…`). Also moved, not listed by RC1-01: `prot_load`'s run step
  and the KIT-112 W3 managed comparison, which compares against `prot_main`'s frozen value. Enrichment messages that
  differ were failing on base and differ only by the added NumSig line.
- `test-workflow-saved-da.R` (`message-changed`; 171 failures on both, laptop references): by value, the moved leaves
  are `dia_qvalue_only` products, dependencies, manifest, head_manifest, display, files (predicted); `files` of
  dia_one_contrast, dia_two_contrasts, dia_repeat, dia_reopened_rerun, dia_non_head, dia_changed_arguments, dia_retry,
  limpa_one_contrast (predicted, SD-2 (a)); the same in W2's worker cases; and, NOT predicted, `dia_refused` products,
  dependencies, head_manifest, display, files (W1 and W2): the refused click captures the head left by
  `dia_qvalue_only`. **Stop rule 3 fires on these five leaves**; reported, nothing edited.

## Inferred
- The dia_refused moves follow from dia_qvalue_only's moved products (nothing runs in that case); no new behaviour.

## Acceptance
- A1 met (G0). A2 met (R1 failing-first; R2b: the new file passes 79 of 79 on head). A3 not met as worded: only
  `dia_refused`'s five leaves are outside section 6 (`prot_load` and the W3 comparison are inside it; RC1-01 mislabelled
  `:3348`, which is `prot_load`). The owner has since allowed them (below). A4 (R3), run D,
  control C: not run (stopped as ordered). A5: NEWS, ledger row, no version or NAMESPACE change; this record.

## Findings
- F6 (folded into KIT-220 at r1.3 by the coordinator): the lipidomics and metabolomics published count tables'
  `up_regulated` and `down_regulated` were always 0:
  `summarise(significant = sum(significant != "NS"), up_regulated = sum(significant == "Up"), ...)` reads the new
  numeric `significant` (`R/func_lipid_da_export.R:311-315`; the metabolomics twin). Shown by `g0-feature.tsv` (0/0
  where the tab says 2/2).
- E-3's `decideTestsDGE` does not exist in edgeR 4.10.1 (above).

## Not checked
- R3 (full suite), control C (GUI-entry keys), run D (re-dumps), the rendered reports, Windows, real data.
- The I04 rework after R2b (the transcriptomics note's wording moved from the frozen GUI configuration into the module,
  `R/mod_omic_da.R`, so that `test-omic-gui-config.R`'s KIT-106 W1 value does not move): not yet run; the merge run covers it.
- The workbook route's interactive volcano still colours at q 0.05 and |log2FC| >= 1 (KIT-230 N3; routed by the coordinator).
- `daAnalysisWrapperFunction` (left for after KIT-218, ticket 3 item 4).

## r1.3 (implementation audit AUDIT-KIT-220-impl-20261125; coordinator's rulings)

**R1 for the lipid/metab repair.** R1b (45501071657) and R2 (45503395846) already showed the four lipid/metab T-PAR
expectations failing on base at published Up/Down 0 / 0 against the tab's 2 / 2, contrary to the ticket's base value
("prot and trx differ"); not a stop rule then; reported as F6 after R2. r1.3 declares them failing-first rows after that
observation; no new tests-only run was made. The I05 top-genes row is new: it fails on base by construction
(`fdr_thresh` is not a formal there) and R2b's base (f925f652) shows it failing; head passes.

**R2b** (45513008193, d0b8672c against f925f652; 67 files; check identical, 6 findings each side). Base shard 3 was lost
at the library install (GitHub HTTP 500; infrastructure); its seven files are compared with R2's base (same package code)
and their head failures are identical to R2's head. Verdicts: 54 `same`; `fixed` for the new file (79 passed, 0 failed)
and for the six files of the lost shard; the moves below. Expectations that now fail ONLY because a stored value awaits
its re-freeze, by leaf (every other line of those messages is identical base and head, MECHANICAL):
- (VAL-191 (2)) `test-omic-da-module.R` `files` leaf `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab`:
  `prot_load` run; `prot_main` run, select, tree_k, tree_height, tree_dynamic, tree_none, dendrogram_off, save_heatmap
  (`8a7cae29…` -> `2477b5a9…`); `prot_no_contrasts` run (`34e7cdf6…` -> `a5562559…`); the KIT-112 W3 comparison reads
  `prot_main`'s run leaf. `test-omic-enrich-module.R`: the same leaf, 18 values (R2; unchanged in R2b's head).
- (i) lipid/metab repair: `test-omic-da-module.R` `files` leaf `publication_graphs/NumSigDeMolecules/lipids_num_sig_de_molecules.tab`
  of `lipid_main` run, select, save_heatmap; `lipid_no_friendly` run, table; `lipid_managed` run, select; and
  `.../metabolites_num_sig_de_molecules.tab` of `metab_main` run, select, save_heatmap; `metab_no_friendly` run, table;
  `metab_managed` run, select: 14 leaves. `test-omics-da-core.R` `omicsDaCoreExportFiles` (`:343`), the one value used
  for lipid and metab (`a36fe768…` -> `8f71c903…`). `tests/testdata/gui-entry/lipid_canonical.json`: 12 keys
  (`results/lipidomics/publication_graphs/NumSigDeMolecules/lipids_num_sig_de_molecules.tab` x6,
  `results_summary/lipidomics/Publication_figures/NumSigDeMolecules/lipids_num_sig_de_molecules.tab` x3,
  `results_summary/lipidomics/Publication_tables/da_lipidomics_num_sig_de_molecules.tab` x3; `7590496b…` -> `6a21cfa6…`).
  The metabolomics GUI-entry captures do not move (`test-gui-entry-metabolomics.R` `same`).
- (ii) `test-workflow-saved-da.R` `savedDaW1Frozen()[["dia_refused"]]`: products, dependencies, head_manifest, display,
  files (W1 and the W2 twin), each equal to `dia_qvalue_only`'s new value (R2; R2b head identical).
- (iii) the `files` leaf of dia_one_contrast, dia_two_contrasts, dia_repeat, dia_reopened_rerun, dia_non_head,
  dia_changed_arguments, dia_retry, limpa_one_contrast (W1; W2 twins), besides `dia_qvalue_only`'s six leaves.
- Outside every list: `test-omic-gui-config.R` KIT-106 W1 (`summary` of transcriptomics gained `box_consumers`). Not
  re-frozen; the wording was moved into the module instead (Not checked).

**Owner's word (owner, 2026-11-25, "allow all three"):** (i), (ii), (iii), each only through the control; order unchanged
(control C, run D and the re-freeze commit after the merge train and a merge of staging, on the coordinator's word).

**Conduct.** The fix was drafted (local commit b29faa34) before the valid failing-first run R1b; R1b still ran on base
code. 5978ce97 was pushed under the instrument-fix title while carrying the whole fix (no run bound it); two
force-pushes with lease on `kit-220` (5978ce97 -> bb092685; 0adde685 -> c788f04a) followed. `git stash list` was run
once (read-only; the stash log was last written before the branch existed).

## Control C and run D (coordinator's word, 2026-11-25; staging 185ebf75 merged as 437e1d9a)

Script `redump.R` (each file in its own process with the dump switches and the GUI-entry record directory). First send
(six runs) wrote no dump (relative directory; VAL-54 instrument fix 1, four runs cancelled); re-run 1: base 2db79d20
(staging + script), head da5dc27b twice. Leaf tables: `control-cd.tsv`, `control-c-preexisting.tsv`.
- **Control C** (45522656187 base; 45522658880, 45522663108 heads): every key that differs base/head is an allowed
  NumSig key: `prot_dia.json` 13, `prot_dia_limpa.json`, `prot_lfq.json`, `prot_lfq_fragpipe.json`, `prot_tmt.json`
  12 each (VAL-191 (2)), `lipid_canonical.json` 12 (owner, (i)); the metabolomics lanes do not move. Head 1 = head 2 in
  every lane. At every changed key the base record equals the stored oracle. Pre-existing: `prot_dia.json` differs
  base/stored at 24 other keys (integration files; left untouched).
- **Run D, module captures** (45522667790 base; 45522671732, 45522675562 heads; head 1 = head 2 byte for byte):
  `test-omic-da-module.R` KIT-106 dump: only the 10 proteomics and 14 lipid/metab NumSig `files` values change, and each
  base value equals the stored one; `test-omic-enrich-module.R` KIT-107 dump: only the 18 NumSig values, base = stored.
  `test-omics-da-core.R` (no dump hook): base passes with the stored `a36fe768…`; head gives `8f71c903…` (R2b).
- **Run D, W1** (`test-workflow-saved-da.R` KIT-112 dump): exactly the 19 allowed leaves change (dia_qvalue_only 6,
  dia_refused 5, the eight other `files`); head 1 = head 2. **But at every one of the 19 the base dump differs from
  the stored value**: the W1 object holds laptop-recorded values, and 88 W1 leaves differ base/stored on the runner
  (the 171 known failures). Under the coordinator's rule ("a changed leaf whose base dump differs from the stored value
  stops the ticket") this is a stop. No re-freeze commit was made, and R3 was not sent.

## Re-freeze and the W1 object (coordinator's ruling, 2026-11-25)

- Re-frozen in 2c7e2959, exactly the leaves where control C / run D showed a change with base equal to stored:
  GUI-entry NumSig keys (`prot_dia.json` 13, `prot_dia_limpa.json`, `prot_lfq.json`, `prot_lfq_fragpipe.json`,
  `prot_tmt.json` 12 each, `lipid_canonical.json` 12); `test-omic-da-module.R` 10 proteomics + 14 lipid/metab NumSig
  values; `test-omic-enrich-module.R` 18; `test-omics-da-core.R:343`. `prot_dia.json`'s 24 pre-existing keys untouched.
- **Not re-frozen (end-of-project re-record list):** `tests/testthat/test-workflow-saved-da.R`, `savedDaW1Frozen()`
  (laptop-recorded; read also by W2): `dia_qvalue_only` products, dependencies, manifest, head_manifest, display,
  files; `dia_refused` products, dependencies, head_manifest, display, files; `files` of dia_one_contrast,
  dia_two_contrasts, dia_repeat, dia_reopened_rerun, dia_non_head, dia_changed_arguments, dia_retry, limpa_one_contrast
  (19 leaves). Runner values: `control-cd.tsv`.
- **A3/A4 restated.** A3: R2-scope files `same` base against head except the new file, the re-frozen files (now passing
  their re-frozen leaves) and `test-workflow-saved-da.R`, whose runner failures are the known ones plus message
  changes at those 19 leaves only. **Expected R3 verdict** (staging dd769671 against the merge): every file `same`,
  except `test-kit220-da-threshold-consumers.R` (`only-on-head`, all passing) and `test-workflow-saved-da.R`
  (`message-changed`, the same failure count, messages differing only at the 19 W1 leaves); `test-omic-gui-config.R`
  `same`; the re-frozen GUI-entry, module and core files `same` (base passes the old values, head the new); check
  findings identical. Anything else is a stop.
- Staging dd769671 merged (8bda2aac; ledger: staging's rows to VAL-211 first, `VAL-<merge>` after); NEWS: the
  workbook clause removed (VAL-209), 5d327b2f.

## For the change report
- **What was wrong.** Proteomics (every workflow and route): the published count table (`.tab`/`.xlsx`), its NumSig
  bar chart (the report's Figure 2), the `Publication_tables` copy and the report's per-comparison count sentence
  counted at q < 0.05 whatever DA threshold was set, while the volcano plots and the methods sentence used the
  threshold set; the interactive volcano highlighted only |log2FC| >= 1 at every setting. Transcriptomics: the
  published volcano (line, colours, labels), significance bar chart, `top_de_genes.tsv` and its report caption, the RUV
  summary line and the interactive page used FDR 0.05. Lipidomics and metabolomics: the published count table's
  `up_regulated` and `down_regulated` were 0 in every analysis (since 2026-01-09 metabolomics, 2026-01-18 lipidomics;
  on main), its `significant` and bar chart right.
- **Who was affected.** Proteomics and transcriptomics users who set a threshold other than 0.05 and used the published
  files or report; every proteomics user of the interactive volcano; every lipidomics and metabolomics user of the count
  table's up and down columns. No statistic (log2FC, p, q, FDR) was affected.
- **What changed in results.** At t != 0.05 the proteomics and transcriptomics published counts, figures and lists
  now equal the tab's at t (e.g. synthetic fixture at 0.01: 65/62 -> 54/56; trx_edge_12 at 0.01: bar 320/299 ->
  229/235). The proteomics count table gains `q_threshold` and the report sentence "(q-value < t)" at every threshold;
  the interactive page colours by q and the run's cutoff, at the default too. Lipid/metab up/down give the counts.
- **Before the application.** The workbook's proteomics table (`_num_significant_differentially_abundant_all.tab`)
  counted at 0.05 since before the application; it now counts at the workbook's threshold.
- **Code.** At d048898d: `R/func_prot_da_model_methods.R:508-512`; `R/func_prot_da_results_summary.R:154-176`;
  `R/func_prot_da_results_methods.R:340-350`; `R/func_prot_da_compat.R:162-166`; `R/func_prot_da_volcano_glimma.R:216-220`;
  `R/utils_workflow_dia_da.R:876-880`; `R/utils_workflow_trx_declarations.R:902-916`, `:1100`; `R/func_trx_da.R:105`;
  `R/func_trx_publication.R:253`; `inst/reports/proteomics/{DIANN,LFQ,TMT}_report.rmd` FormatDESummaryTable;
  `inst/reports/transcriptomics/transcriptomics_report.rmd:152`; `R/func_lipid_da_export.R:311-315`;
  `R/func_metab_da_results_output.R:313-317`. At c788f04a / d0b8672c: the same functions (`printCountDaGenesTable`
  `q_val_thresh`; `q_threshold` after the bind; `Status` with `lfc_thresh`; identity `lfc`; `fdr_thresh` in the trx
  calls; `decideTests.DGEExact` status; up/down before `significant` with `relocate`).
- **Evidence and its limits.** G0 (45503391942 base; R2 45503395846 head), R1b 45501071657, R2, R2b 45513008193. A
  synthetic proteomics fixture, trx_edge_12 and an 8-row lipid/metab table only; no real data; GUI-entry lanes not yet
  run for proteomics (control C).

## Past results
- Wrong when the threshold set was not 0.05: proteomics `NumSigDaMolecules/da_proteins_num_sig_da_molecules.*`,
  `Publication_tables/da_proteomics_num_sig_da_molecules.tab`, the report's count sentence and Figure 2 (DIA-NN, LFQ,
  TMT), the workbooks' table; transcriptomics published volcano, bar chart, `top_de_genes.tsv` (captioned 0.05,
  honestly), interactive page.
- Different at any threshold: the proteomics interactive volcano (|log2FC| >= 1).
- Every lipid/metab count table ever published shows 0 up and 0 down; its total and bar chart are right.
- After the fix, a transcriptomics report rendered from a project whose DA ran before it at t < 0.05 captions the
  earlier 0.05 top-genes list "FDR (BH) < t"; a direction with fewer than 10 genes below t then shows rows above t.
- The workbook route's interactive volcano still colours at q 0.05 and |log2FC| >= 1 (KIT-230 N3).
- How a user can tell: `analysis_record.json` `DA$qvalue` != 0.05 (proteomics, trx); its `software.sha`/`build` tell
  whether the files predate the fix; recount from `DA_results_*.xlsx` or `de_results_*.csv` at the stated threshold.
