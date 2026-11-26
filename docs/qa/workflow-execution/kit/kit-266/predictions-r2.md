# R2 declaration (before the run), restated on the merged tree

Purpose: regressions of the fix against the tests-only base, same runner image, 57 files, 5 shards per side, R CMD check on both
sides. Base 5e3de9c0 (staging 36e3180b with the branch's tests: the new file, the declared edits, the G3 edit, the corrected E12 (c) and
caption rows). Head 41035e75 (tree equal to 9b3481f7: the fix on that base, KIT-272's sentences merged in the limpa Results
branches). The files: the ticket's grep at the merged tree (52 test files, including the new file, kit259b, kit271's and kit272's
neighbours) plus test-gui-entry-proteomics.R, test-kit252-managed-failure-text.R, test-kit259-pearson-histogram-figure-failure.R,
test-workflow-operation-registry.R and test-kit272-prot-box-panel-text.R (57). test-omics-module-globals.R is not in scope: no R/ file
named (mod|func)_(lipid|metab|omic|omics)* is edited (14/0 expected if run). No transcriptomics file is in scope, so
MULTISCHOLAR_TRX_MANIFEST_NOTE is not needed. Ceiling 980 job-minutes (65.33 core-hours), expected 608.6 (40.57).
Decision: A3, A4, A6.

## Predictions (verdicts of compare.tsv, base against head)

- test-kit266-limpa-combined-qc-figure.R: `fixed` (base about 100 failed, head 0 failed, 0 errors; G rows pass on both).
- test-kit259b-second-correction.R: `fixed`; base-only failures 5 (C01-G1 2 expectations, the limpa (ii)-E3 iteration 3).
- test-kit227-report-ruv-sentences.R: `fixed`; base-only failure 1 (E7 `skipped`).
- test-prot-norm-core.R: `fixed` by the declared edits only: base-only failures on the W3 limpa block (`:1640-1672`: the operations lists and the
  `composite_files` set-aside, both modes) and the A111-01 limpa block (length 10, 16 operations, input index 14); stored-value failures
  identical on both sides (1547/307 on staging); a base-only failure on any other line is a stop. The ticket predicted 2 base-only
  assertion groups; the exact count is not predicted beyond that.
- test-prot-norm-worker-figures.R: `fixed`; base-only failures 3 (`:136`, `:137`, `:144-145`).
- test-gui-entry-proteomics.R: head failures = base failures (92/101 on staging) plus the limpa steps from `norm` on that pass on base
  (they gain `proteomics_composite_QC_figure.*`, `pre_norm_*` keys and the changed outputs/store entries, ticket section 8); compared with the
  four file-name key lines masked (the KIT-271 lines, merge train 26's method); no classic lane (prot_dia, prot_lfq, prot_tmt,
  prot_lfq_fragpipe) moves.
- test-kit272-prot-box-panel-text.R: `same` (373/0 on both sides).
- test-omic-norm-module.R: `same` (631/8, standing).
- Every other file: `same`, with the pre-existing failure lists identical on both sides (standing: test-workflow-saved-da.R 755/196, test-kit271
  stored digests as on staging).
- R CMD check: identical on both sides (no `man/`, `NAMESPACE` or DESCRIPTION change).
Stop rule 3: any unpredicted verdict or a check difference; one correction pass (R2, R3, M1 again). Budget: running total 8.74 core-hours
before the sends; M1 (re-send) ceiling 4.00, R3 ceiling 4.00, R2 ceiling 65.33 give a worst case of 82.07 plus the correction pass 76.33 =
158.40 against 163.46; no further run is sent that could pass 163.46.

# R3 declaration (render, scripts job on the head)

Binding: head 41035e75; script docs/qa/workflow-execution/kit/kit-266/render.R; ceiling 60 job-minutes (4.00), expected 30. Predictions are
in the script (the checks of V0, V1, V2). V0 reference restated: LFQ, TMT and DIA texts equal (normText) to the same lane rendered with staging
36e3180b's templates (which hold KIT-272's bullets), not to the committed KIT-236 render (old bullets); drawings 3, 3, 3 and 2 for the limpa lane.
Limpa V0: the new caption line ("first two columns ... a third column, when the figure has one"), the three-stage covariates sentence (method
none) with KIT-272's two bullets, no "Density plots", QC_figures holds the composite PNG and PDF, PNG 1950x3900. V1: one added sentence, panels
a, f, m with limpa words. V2: PNG 1500x3300 and the two-stage sentence. Decision: any difference outside the limpa lane, or a letter or caption that
contradicts the figure, is a defect; a title clipped in the PNG (R3 opens it) is shortened (cyclic loess title not rendered by this lane).
