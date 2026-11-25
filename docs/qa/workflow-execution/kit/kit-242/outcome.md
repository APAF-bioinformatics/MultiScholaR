# KIT-242 gate G0: outcome

Ticket r1.2 section 6. Branch `kit-242-g0` (script and records only; never merged). Base `v0.9-staging`
3eeaf100110e661944f50c427d266d9b3f6a21df. Script `tools/kit242/g0.R`.

- Run 45584784285 at 33f11d9e (3eeaf100 + script): every case; 21.3 job-minutes.
- Run 45588916720 at 2942615f (VAL-54 instrument rerun 1 of 3, C9 only): the first run's proteomics arms of C9
  failed in the script's own object builder ("`.data` must be a data frame without row names", the protein table
  kept the matrix's row names); no result was produced. Fix: row names dropped. No expectation, prediction or stop
  rule changed. 11.6 job-minutes.
- C10 reuses run 45579225993 (obs-itsd, package commit 4a82ddff = 3eeaf100 + its script), tables in
  `g0-c10-45579225993/`.
- Core-hours: (21.3 + 11.6) x 4 / 60 = 2.19 (G0 expected 4.0, ceiling 8.0; instrument line 1 of 3 used).
- Seeds: case k draws with `set.seed(242000 + k, "Mersenne-Twister", "Inversion", "Rejection")`: C1, C8, C11
  242001; C2 242002; C3, C3m 242003; C4 242004; C6 242006; C7 242007; C9 242009; C12 242012; C14 242014. C12 (c)
  uses KIT-188's `generate("L12", "N1", 401)` (its own seed), sim-k188.R at e2d18a79 (SHA-256 in
  `g0-45584784285/provenance.txt`).

Tables: `observations.tsv` (case, class, arm, item, value), `conditions.tsv` (every warning, error, kept message and
logger line), `tab-log.tsv`, `notifications.tsv`; lane records under `lanes/`.

## Observed, case by case (S = sample B_6; A1, A2 the two assays)

| Case | Observed | Against the prediction |
|---|---|---|
| C1 | lipid and metab, declaration: warning "Assay 'A2': At least one sample column contains only NA values. Skipping normalization."; A2 equal to its input (all.equal, tolerance 0), A1 normalised; kernel `applied` TRUE; S4 args `normalisation_method` cyclicloess; record step parameter `normalisation_method=cyclicloess`, outcomes only rows/columns; DA tables for A1 and A2. Worker (`omicsWorkerRun`): handoff status `complete`, no error, no warning text and no `applied` field in the handoff. Managed tab: log "Applying between-sample normalization (method: cyclicloess )", "Between-sample normalization complete"; `normalization_complete` TRUE; notifications only "Applying ... step ..." and "Normalization pipeline complete!" | as predicted |
| C2 | quantile, A2's S with one value: "Assay 'A2': Error during 'quantile' normalization: need at least two non-NA values to interpolate. Returning unnormalized data for this assay."; A2 equal to input; `applied` TRUE. Planted cyclic-loess error: the same warning with "planted" for A1 and A2 (instrument worked); both assays equal to input; `applied` TRUE | as predicted |
| C3 | regex matching nothing in A2: "Assay 'A2': No ITSD features identified using regex '^IS_' in columns 'LipidClass' (metab: 'Ontology'). Skipping normalization."; A2 unchanged with its 5 `IS_` rows kept; A1 normalised, standards removed; args `applied` TRUE, `method_type` average_centered, features block A1 only; record `itsd_applied` TRUE. Both assays unmatched: the warning twice, both unchanged, no features block, `itsd_applied` TRUE. Managed tab: "ITSD normalization complete", `normalization_complete` TRUE, no notice. **Not predicted:** the S4 args already hold `itsd_counts_per_assay` at base (lipid `:261`, `:375`; metab `:490`, `:604`), for matched assays only ({"A1":5}); the record does not carry it | as predicted, plus the args field |
| C3m | manual ids for A1 only: "Assay 'A2': No manual ITSD features provided in itsd_feature_ids. Skipping normalization."; A2 unchanged; record `itsd_applied` TRUE; tab log names only A1 ("Assay A1 : 5 ITSD features selected") | as predicted |
| C4 | ITSD: "Assay 'A2': Aggregate ITSD signal (SampleNormFactor) is NA or zero for samples: B_6. Normalized values will be NA for these samples."; S has 0 finite values in A2 after ITSD and log2; between-sample: C1's warning, A2 equal to input; record as C1. **Legacy RUV-III-C (Automatic): A2 was dropped whole, not S's column.** At every percentage the search's cancor step logged "Assay 'A2': Error during NIPALS imputation: Some rows or columns are entirely missing: remove them first.. Skipping." (`R/func_lipid_s4_ruv_methods.R:158`), then "RUV: no valid percentage for assay: A2", then "Assay 'A2': Invalid K value resolved (NA) ... Skipping." (`:481-484`); A2 absent from the RUV state; record `ruv_applied` TRUE, `k_per_assay` {A1 1, A2 NA}; nothing at `:713-726` for A2 (A1: 40x12 -> 40x12, no feature removed). The design kept S (A1 is first and holds S). Correlation filter at 0.5: S kept in the design (only A1 left). DA: a table for A1 only | **contradicts the prediction** (stop rule 0): the empty column does not reach `:719-722`; it makes the automatic search fail and the whole assay is dropped silently; the correlation filter's handling of an empty column was not reached |
| C5 | lanes as committed. lipid_canonical: regex `^IS_\|ISTD\|Internal\|_IS` in `LipidClass` matched 0 rows; metab_lc: annotation column NA, so ITSD skipped at the "pattern columns" row (args `itsd_pattern_columns` NA). Both: tab log "ITSD normalization complete"; record `itsd_applied` TRUE; study_parameters.txt "ITSD Normalization: * Applied: Yes"; session summary "- ITSD applied: Yes"; no notification at the step | as predicted (metab_lc reaches the skip by the "pattern columns" row, not the regex row) |
| C6 | guard on, A2 of 3 features: in process `omics_cyclic_loess_too_few_features` (assay A2, sample A_1, fewest 3, minimum 10); `omicsWorkerRun`: handoff `failed` with the guard's message | as predicted |
| C7 | DA: `qvalue_warnings` for (groupB-groupA, A2) and also for A1; record v2 entries hold only contrast, assay, entity_id, rows (no `fdr_method`). A1's fallback came from `qvalue()` itself failing on 40 null features ("missing or infinite values in inputs are not allowed, using BH") | as predicted, and reachable at 40 features |
| C8 | C1's record valid; replay reproduces the warning once, states equal the history, no mismatch | as predicted |
| C9 | proteomics, empty column: cyclic loess completes, S empty, the others normalised; quantile stops ("'x' and 'y' lengths differ"); **scale completes with every value of the table missing (0 of 720 finite)**; one-value quantile stops ("need at least two non-NA values to interpolate"). Transcriptomics TMM with a planted `calcNormFactors` error stops ("planted"); the unplanted control completes | as predicted; Part C stands |
| C10 | run 45579225993: declarations and KIT-172 fixtures, "No ITSD features identified using regex '^IS_' in columns 'LipidClass'/'Ontology'" per assay, `IS_` rows still present after ITSD, features block absent, record `itsd_applied` TRUE | as predicted |
| C11 | guard on, C1's object: `omics_cyclic_loess_too_few_features`, sample B_6, fewest 0, raised before any warning | as predicted |
| C12 | (a) S4, controls found in A1 only: A2 absent, logger "Assay 'A2': No control features identified after resolution. Skipping." (`:667`), no `warning()`. (b) found in neither: "No assays were successfully processed by RUV-III. Returning original object.", both assays equal to input. (a-decl) A2 without a complete feature, legacy Automatic through the declaration and the managed tab: A2 dropped, "RUV: no valid percentage", invalid-k skip; record `ruv_applied` TRUE, `k_per_assay` A2 NA; tab `ruv_complete` TRUE; DA table for A1 only. (c) KIT-188 L12 N1 401 vd20 LA: status `uncorrected`, both assays k NA, optimisation "No valid percentage found", final skip `:481-484` for both; inside the search, LCMS_Pos failed at NIPALS (entirely missing row or column, `:158`) and LCMS_Neg at "Some assay features not found in named logical 'ctrl' vector" (`:220`) at every percentage | as predicted; (c) adds the search's own failure points |
| C13 | lipid_canonical: 0 regex rows and 0 candidates per assay (the lane's data hold no standard-like row). metab_lc: no annotation column (NA), the regex count is not applicable, and the ITSD Selection builder errors with that column ("replacement has 0 rows, data has 3") | recorded only; metab_lc's builder error not predicted |
| C14 | MS-DIAL default pattern on `Metabolite name`: the 5 'w/o MS2:' rows used as A1's standards and removed; the 2 '(d7)' rows kept as features; A1's other features changed; features block and args `itsd_counts_per_assay` A1 = 5; record `itsd_applied` TRUE. Integration input with ITSD unticked: the 5 rows marked `is_internal_standard` | as predicted |

## Inferred (not executed)

- OD-1 (c) and OD-3 (d): on the shipped default RUV route, an assay that keeps an empty sample column fails the
  automatic search (NIPALS refuses an entirely missing column) and is dropped whole; so "nothing lost beyond that
  sample in that assay" holds only if the search and the correction exclude the empty column (C4).
- The "no k and no controls" cause of KIT-188's failures is the search's own cancor step failing at every percentage
  (NIPALS on an entirely missing row or column; a control-name mismatch), not the absence of usable controls (C12 c).

## Not checked

- What the correlation filter does with an empty column in an assay that is kept (A2 was dropped before it, C4).
- Why LCMS_Neg's named control vector misses assay features (C12 c, `:220`); recorded as a new adjacent item.
- A separately started worker's log; DT pre-selection in a browser; real data of any kind.

## C4b (ticket r1.3, optional; run 45594264982 at 30e89881 = 3eeaf100 + script; 8.6 job-minutes, 0.57 core-hours)

C4's object (seed 242004; S's standards missing in A2), ITSD median, log2, cyclic loess, **RUV skipped**, then the
correlation filter at 0.5, then the DA; lipid and metab. Prediction, written in the script and emitted before any
value was read: "B_6 removed from A2 and from the design (so from every assay's DA)". Records: `g0-c4b-45594264982/`.

| Item | Observed (lipid and metab alike unless stated) | Against the prediction |
|---|---|---|
| Up to between-sample | as C4: S has 0 finite values in A2; "At least one sample column contains only NA values. Skipping normalization." for A2 | as before |
| Correlation filter | A2's pairs naming B_6: 5 rows, 0 finite correlations (A1's: 5 finite). After the filter: design 11 samples, B_6 not in it; A2 has 11 samples, B_6 absent | **as predicted** |
| A1 after the filter | A1 still holds B_6's column (40 finite values), so the assays hold different samples and A1 holds a sample the design lacks | not predicted |
| DA | lipid: **refused**, "Assay 'A1' has sample column(s) missing from design: B_6"; metab: completed, tables for A1 and A2 (40 rows each), B_6 absent from the design | not predicted; lipid and metab differ |

Observed: a sample emptied in one assay is removed from the whole design by the correlation filter, as predicted.
Inferred: under OD-1 (b) or "Skip RUV" with the filter applied, the sample leaves every assay's DA (metab), or the
lipid DA stops with a message naming the other assay (lipid). Not checked: why the metab DA accepts an assay column
the design lacks; whether the GUI's DA tab reaches the same refusal.

G0 core-hours now 2.19 + 0.57 = 2.76 of 8.0.
