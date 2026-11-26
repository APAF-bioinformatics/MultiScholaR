# KIT-268 outcome: the notice names the comparisons whose q-values fell back to Benjamini-Hochberg (VAL-312, VAL-314, VAL-351)

Branch kit-268. Fix commits: 28b88665 (the notice, the sentences, NEWS), f5caefb8 (the by-position repair of `daQvalueFallbackText`),
dedc4c36 (the test-only repair of the A112-03 test, VAL-351). Run table: `runs.tsv`. Wording drafts: `wording.md`. Core-hours: 36.6 of 56.33.

## What was run (observed)
| run | what | result |
|---|---|---|
| 45705978519 | G0 attempt 1 | instrument defect (relative output path); no observation; fixed under VAL-54 (1 of 3) |
| 45708195790 | G0 attempt 2 | positive control live (6 keys); the observation below |
| 45712438665 | R1 tests-only | every prediction held exactly (new file 12/11, KIT-236 file 307/6, da-module 1286/83 against base 1288/81); T6 passed |
| 45715855490 | re-freeze of lipid_canonical on the runner | NOT committed: the diff holds KIT-217's state digests (16 steps) besides the header and the one warning |
| 45721176971 | R2, 39 files, check | stopped on four misses (T2b; new file 29/1; KIT-236 313; enrich 1005/65); records corrected in runs.tsv |
| 45726799791 | R3, correction pass | new file 30/0, KIT-236 313/0, GUI-entry as predicted, check identical; one miss (T2b); enrich 1011/59 (R2's six not reproduced) |
| 45736314028 | R4, VAL-351, test-omic-da-module.R only | prediction held: head 1288/81/0/0 = base 1288/81/0/0, verdict same, the A112-03 test 5/0 |

## G0 (observed, run 45708195790, base 4d588ba3)
The positive control (tied_pvalues lipid through `completeOmicDaRun`) printed 6 fallback keys. Lipidomics `lipid_canonical`: one fallback,
`groupKO-groupWT` in LCMS_Pos and in LCMS_Neg (an exact match on the contrasts table; friendly name KO_vs_WT). `metab_lc`, `metab_gc`,
`metab_combined`: the completion ran, 0 keys. `test-omic-summary-module.R`: 4 calls (transcriptomics), 0 keys. `test-omic-da-module.R`: 15
calls, 0 keys. Stop rule 3 (mapping) not triggered on the observed routes.

## Who gains a notice (by omic)
- Proteomics (DIA iq, LFQ, TMT, DIA with limpa; in process and managed): the existing generic notice now names the comparisons. On `main` (436cdbe8, `git show main:R/mod_prot_da.R`) a dialog named the contrasts; staging's generic notification had lost the names (GUI audit A112-03).
- Lipidomics and metabolomics: a notice they never had (on staging or `main`), with the assay, via `completeOmicDaRun` (managed and in process).
- Transcriptomics: none (edgeR's FDR is Benjamini-Hochberg throughout; the list is always empty).
- Users with no fallback comparison: no change.
- The report sentence follows the record's per-comparison method: for proteomics now; for lipidomics and metabolomics only once KIT-242's record hunk is merged (until then their report stays in the rule branch).

## Acceptance criteria
| # | Criterion | Evidence |
|---|---|---|
| A1 | proteomics notice names exactly the fallback comparisons; none when none (both routes) | T1 (R3: 30/0); T2b repaired (R4: the A112-03 test 5/0 on both routes, one notice naming B_vs_A, none without a fallback). T2b's earlier failure was the test's (it asked for A_vs_B in a one-contrast fixture), plus the defect of the fix repaired in f5caefb8 |
| A2 | lipid/metab tab shows the notice with the assay; transcriptomics none | T3 (R3); the GUI-entry lipid lane (R2, R3) shows "KO_vs_WT in LCMS_Pos; KO_vs_WT in LCMS_Neg"; managed delivery by reading (`omicsDaTabResults` keeps `qvalue_warnings`) |
| A3 | six reports: existing sentence then the added one; none or no facts unchanged | T4 (R3: KIT-236 file 313/0) |
| A4 | record reason holds the existing text and the sentence; "qvalue" otherwise | T5 (R3: new file 30/0); lipid/metab half by composition with KIT-242's E10 once merged |
| A5 | BH = q/pi0 at every index, ties included; q <= BH; the trigger | T6 (R1, observed on the runner's qvalue 2.44.0) |
| A6 | record and report sentences agree | one text (VAL-314 verbatim) in both; T4, T5 |
| A7 | no q-value, p-value, count, table or file changes | by reading (only notice, report text, one record string); R2: the other 33 files same |
| A8 | R CMD check identical to base | R3: identical (0 errors, 3 warnings, 3 notes, nothing only on either side) once the records' paths were shortened |
| A9 | stored GUI-entry steps | NOT executed (re-capture); see below |
| A10 | NEWS, ledger row, outcome.md with "For the change report" and "Not checked" | NEWS under "Fixes"; ledger row `VAL-368`; this file |

## Stored value that waits for the end-of-project re-record (A9, not edited)
Only `tests/testdata/gui-entry/lipid_canonical.json`, step `da` (line 2402; notifications 2625-2638): one added `warning` notification at
position 3 (text in wording.md and the ledger row). The owner's condition for re-capturing it (VAL-314 (3)) could not be met: the file
already differs in 16 steps' state digests since KIT-217 (VAL-324), and by the coordinator's ruling it is not re-frozen here. The evidence
(`refreeze/`) shows the notification's exact text and that nothing else of this ticket moves. The three metabolomics oracles and the W1 and
summary captures do not move. The six `stringdb_plot` failures of R2 in `test-omic-enrich-module.R` are not a stored value of this ticket (below).

## Known flake (audit I06)
`test-omic-enrich-module.R`, "(b) in process: refusals", element `stringdb_plot` (a constant-text plotly widget): three digests across
processes (b54fa6cb frozen; 62479f92 in kit-201 r3 on the base side; a3f2703c in R2 on the head side), each seen on one side only, never
reproduced (R3: 1011/59 on head, same as base). Likely mechanism (inferred): plotly's per-plot hex id escaping the capture's mask.

## Convention notes (audit I08; code left as it is)
38 added lines exceed 120 characters (max 140) against R.md's 100; `function(result)` is used for short lambdas in new code where R.md has
`\(x)`; 2-space indentation follows each file's existing style (the existing-code scope rule).

## Not checked
- The notice's wording and the report/record sentence through the owner's pass (VAL-41); the notice is the ticket's draft.
- The proteomics routes through a real GUI session (only the W3 fixture and T1).
- That a lipid/metab record carries the method per comparison (waits for KIT-242).
- The six report templates rendered end to end with a fallback comparison beyond T4's harness.
- The full 39-file scope after the last fix (R2 ran it at the earlier fix; the merge train runs it); `test-omics-module-globals.R` and `test-kit196-da-no-value-features.R` (audit I03) were not run.

## For the change report
Who was affected: users whose differential abundance has a comparison where Storey's pi0 cannot be estimated (small tables, or tables with no p-value at or above 0.95), so Benjamini-Hochberg values sit in the 'q-value' column.
What was wrong: the label hid a different scale. On staging the proteomics notice named nothing (`main`'s dialog had), and lipidomics and metabolomics said nothing on either branch; the reports named the comparisons but not the scale.
Code at the base (4d588ba3): `R/mod_prot_da_handlers_run_analysis.R:13-17` (generic text), `R/mod_omic_da.R:461-480` (no notice), `R/utils_workflow_analysis_record.R:373-378`, `ReportFdrSentence` in the six templates. At the fix: `R/func_omics_da_shared.R` (the shared text), `warnDaQvalueFallback`, `completeOmicDaRun`, the record string, the templates.
What changed in results: nothing numeric; only the notice, one report sentence and one record string.
Evidence: R1 run 45712438665 (failing first by value), R3 45726799791, R4 45736314028; T6 on the runner (BH = q/pi0 at every index, qvalue 2.44.0); G0 45708195790.
