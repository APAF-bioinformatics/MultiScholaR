# KIT-280 predictions (committed before any run)

Base counts are from merge train 37's full-suite tables (`merge-train-20261126aa-tests-head.tsv`, code of merge head b65f67706 = staging e4a6d4aea). Work item: rule B in the proteomics missing-value filter, the replay marker `threshold_values_in_filter`, the flag per comparison, the texts. Exit: A1 to A12 of the ticket.

## R1: tests-only commit (package = this branch's tests commit + records, old R/; base "none"; 1 job; ceiling 75 min = 5.00 core-h; expected 38.6 min = 2.57)

Files: `test-kit280-prot-threshold-values.R`, `test-prot-missing-value-counts.R`, `test-kit227-report-ruv-sentences.R`, `test-prot-summary-report-templates.R`. Each new test must fail on the old code at the stated value (an expectation, not an error); T7 and T9 pass.

| Test | Old code | Reason |
|---|---|---|
| T1 | fails | kept 39 (P_low1 removed), protein_count 39, not 40 |
| T2 | fails | P_low1 absent from the DA table: nrow 0, not 1 |
| T3 | fails | no `low_intensity_value`, no `low_intensity_groups` column |
| T3b (both chains) | fails | the route completes; the column is absent |
| T3c | fails | no flag column; Pmin removed (passes that expectation) |
| T3d | fails | no flag column (the DA with `1h`/`2h` completes) |
| T4 | fails | `low_intensity_value` NULL, not FALSE / TRUE |
| T5 | fails | marker row absent (NULL), declarations lack the entry, the record of (a) lacks the marker; (b) and (c) expectations pass (39 kept, no flag) |
| T5b | fails | the normalisation step lacks the marker; the filter-step, flag and replay expectations pass |
| T6 | fails | strict keeps 38, the clean-up removes P_low1 |
| T7 | passes | routes and defaults today |
| T8 | fails | no marked-count line; setdiff of the column names lacks `low_intensity_value` |
| T8b | fails for token "present" (no flag column in the page); "missing" passes |
| T9 | passes | static check of the three limpa compute functions |
| T10 | fails | no new sentence, no neutral sentence (the earlier sentence prints without a record), DIA reads the peptide percentile |
| T10b | fails | label and help are the earlier ones |

`test-prot-missing-value-counts.R`: 3 failures, 64 passed (the three edited expectations: T4's second case, T5 strict `protein_count`, T5 strict kept set; the file's 4 warnings unchanged). `test-kit227-report-ruv-sentences.R`: 3 failures (G1 `file_only`, LFQ, TMT, DIANN: the base's sentence is substituted by the neutral sentence, the old head still prints the earlier one), 183 passed. `test-prot-summary-report-templates.R`: 123 passed, 0 failed.

## R2: base against head (head = the fix commit; base = staging e4a6d4aea)

Main request: the 105 files other than `test-kit247-design-resave.R`, 5 shards x 2 sides (67.3 min each expected, 75 ceiling) and the check on both sides (40 expected, 90 ceiling). Second request: `test-kit247-design-resave.R` alone, one job per side (expected 79.7 min, ceiling 95), timeout 120 min. Expected 912.4 job-min = 60.8 core-h; declared ceiling 1,120 min = 74.67 core-h.

Verdicts, `compare.tsv`:
- `test-kit280-prot-threshold-values.R`: `only-on-head`; every test passes on head.
- `test-prot-missing-value-counts.R`: `same` (67 passed / 0 failed on both sides; the edited expectations pass on head, the file on base is the unedited one).
- `test-kit227-report-ruv-sentences.R`: `same` (186 / 0; G1's `file_only` substitution passes on head).
- `test-omic-qc-module.R`: `new-failure`: base 156 passed / 4 failed, head 153 / 7. The three new failures are W5 (:420-421), DIA, TMT, LFQ (the HTML digests of `omicQcFrozen$tabs`); the two identity expectations (:271-272, `tabs` and `subtabs` as wholes) fail on both sides with changed messages. `prot_DIA_LIMPA` unchanged.
- `test-prot-norm-core.R`: counts unchanged (base 1,555 / 307, head 1,555 / 307); the 16 `protNormW1Frozen` values (:579-580, :606-607, :634-635, :665-666, :698-699, :732-733, :766-767, :800-801) already fail on base and keep failing with a changed message: `message-changed`.
- `test-gui-entry-proteomics.R` (base 92 / 101): any change of verdict is confined to the lane steps named in the ticket: `prot_dia` import to summary_export (28 of 29 steps; the `dynamic_qc_tabs` UI digest from import on, and from the filter on the states digests, the filter text and the DA long-table artifact `da_proteins_groupKO-groupWT_long_annot.tsv`), and `prot_lfq`, `prot_lfq_fragpipe`, `prot_tmt` filter to summary_export (15 steps each): 73 steps. A step already failing on base stays failing (message changes); a step passing on base may turn red only inside these 73. A failure outside them is unpredicted.
- `test-kit247-design-resave.R`: `same` (188 / 0). It drives the lanes through the filter and checks that a re-save is faithful; it holds no stored digest of the filter's result (INFERRED from the table: 188 passed on base). A red there is unpredicted.
- every other file in the scope (103): `same`, passed/failed counts equal to the base's.
- check: `check-compare.txt` shows no finding only on head (same notes, warnings and errors both sides).

Any verdict not listed stops the ticket and goes back to the coordinator before anything is changed.

## Stored values that move (listed for the end-of-project re-record; not edited)

- `tests/testdata/gui-entry/prot_dia.json`, `prot_lfq.json`, `prot_lfq_fragpipe.json`, `prot_tmt.json`: the 73 steps above.
- `tests/testthat/test-prot-norm-core.R` `protNormW1Frozen()`: the 16 values above (already failing).
- `tests/testthat/test-omic-qc-module.R`: `omicQcFrozen$tabs` `prot_DIA`, `prot_TMT`, `prot_LFQ` (:163, :167, :168) and `$subtabs` `prot_protein_DIA`, `prot_protein_TMT`, `prot_protein_LFQ` (:182, :184, :185): three green expectations turn red by declaration (W5), two identity messages change.
