# KIT-220 outcome (r1.2; stopped before run D, control C and R3, as ordered)

Branch `kit-220` from `origin/v0.9-staging` d048898d (before KIT-204, at the coordinator's order). Commits: 0804ed93
tests first; bb092685 and 19b9a4dc instrument fixes 1 and 2 (VAL-54; tests and G0 only, no expected value changed);
c788f04a the fix. Runs: `runs.tsv`. About 27.1 core-hours used (expected to this point 3.0 + 2.0 + 29.5 = 34.5).

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
- A1 met (G0). A2 met (R1 failing-first; R2 head passes except F6's pre-existing rows). A3 not met as worded: the
  dia_refused leaves, prot_load and the W3 comparison are outside section 6's list (stop rule 3). A4 (R3), run D,
  control C: not run (stopped as ordered). A5: NEWS, ledger row, no version or NAMESPACE change; this record.

## Findings (not in scope)
- F6: the lipidomics and metabolomics published count tables' `up_regulated` and `down_regulated` are always 0:
  `summarise(significant = sum(significant != "NS"), up_regulated = sum(significant == "Up"), ...)` reads the new
  numeric `significant` (`R/func_lipid_da_export.R:311-315`; the metabolomics twin). Shown by `g0-feature.tsv` (0/0
  where the tab says 2/2).
- E-3's `decideTestsDGE` does not exist in edgeR 4.10.1 (above).

## Not checked
- R3 (full suite), control C (GUI-entry keys), run D (re-dumps), the rendered reports, Windows.
- `daAnalysisWrapperFunction` (left for after KIT-218, ticket 3 item 4).

## For the change report
- What was wrong: proteomics published count table, NumSig figure and report count sentence counted at q < 0.05
  whatever threshold was set; the proteomics interactive volcano highlighted only |log2FC| >= 1; transcriptomics
  published volcano, bar chart, top-genes list, summary line and interactive page used FDR 0.05.
- Code at the base (d048898d): `R/func_prot_da_model_methods.R:508-512`; `R/func_prot_da_results_summary.R:154-176`;
  `R/func_prot_da_compat.R:162-166`; `R/func_prot_da_volcano_glimma.R:216-220`; `R/utils_workflow_trx_declarations.R:902-916`,
  `:1100`; `R/func_trx_da.R:105`; `R/func_trx_publication.R:253`. At c788f04a: the same functions.
- Evidence: G0 (base against head values above), R1 45501071657, R2 45503395846. Synthetic fixture and trx_edge_12 only.
