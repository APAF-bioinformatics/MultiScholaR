# KIT-224 outcome (r1.2; R2's three unlisted files accepted by the owner by file, 2026-11-25)

Branch `kit-224` from `origin/v0.9-staging` cf8af245. Commits: cea72da6 (tests first, nothing under `R/`), 12205bd4 (the
fix, NEWS, ledger), 5f261d89 (stored references accepted by file, OD-4), the provenance re-dump of
`test-omic-norm-module.R` and these records. Decisions: owner, 2026-11-25 ("yes to all, and extra arm"; OD-0 to OD-5).

## Observed

**G0 and R1** (run 45508928818, tests-only commit cea72da6; `g0/`, `r1/`). Every E row failed on base at its stated
value; E4 (ii), G1, G2 passed (stop rule 1 not triggered):

| Row | Base value |
|---|---|
| E1 lipid, complete | largest level change 10.23 (pos), 11.59 (neg) log2 |
| E1m lipid, 5 % missing | 7.56, 5.98 |
| E2 metab | 10.23, 11.59 (identical to lipid) |
| E3 lipid_ruv Automatic (KIT-172 seed 172101) | 22.37; record lists no marker |
| E4 (i), (iii), (iv) | marker absent; the two rules' corrected tables identical |
| E5 transcriptomics | file and integration input: largest gene-mean change 2.431; `correction$levels` absent |
| E6 (T1b seed 172102, x 2^10) | 22 calls of 516 change, max abs(delta t) 0.4205, max abs(delta logFC) 6.4e-13 (equal to KIT-172 r2) |

G0's site E5 stopped in the KIT-040 corpus manifest check (laptop-recorded digest; the runner switch
`MULTISCHOLAR_TRX_MANIFEST_NOTE=1` was not set). R1, on the same commit, gives E5's base values (above), so G0 was not
re-sent.

**R2** (run 45510911807; 12205bd4, tests unedited, against cea72da6; `r2/`). R CMD check: identical findings on both
sides (0 only on head). The new file: 69 pass, 0 fail (`fixed`). `test-kit129-run-path.R`: `same` (711 pass).
`test-prot-ruv-levels.R` (G3), `test-trx-kernels.R` (G4): `same`. Moves:

| File | Verdict | What moved | Approved (OD-4)? |
|---|---|---|---|
| test-lipid-12c-ruv-s4-characterization.R | new-failure 10 | W1 (8), forwarding (2) | yes |
| test-metab-03a-ruviii-c-varying-characterization.R | new-failure 10 | W1 (8), forwarding (2) | yes |
| test-kit127-omics-ruv-factors.R | new-failure 4 | :241-242 | yes |
| test-omic-norm-module.R | new-failure 12 | the 12 `provenance` fields only (the dumps differ in exactly those 12 lines, `norm/base-head-diff.txt`) | yes |
| test-integration-input-transcriptomics.R | new-failure 2 | :80, :120 | yes |
| test-workflow-trx-declarations.R | message-changed | the export file's SHA-256 at two paths | yes |
| test-module-ci-trx-summary.R | message-changed | the export file's SHA-256 | yes |
| test-workflow-omics-declarations.R (KIT-073) | 22 -> 23 failing, message-changed | registry `fixed`; states lipid/metab `_ruv_corrected` and `_correlation_filtered`; correlations; DA results. No earlier state | deferred; R1-01 holds |
| test-workflow-omics-binding-steps.R (KIT-094) | message-changed | `records`: the RUV step's parameters gain the marker (9 -> 10) | deferred; R1-01 holds |
| **test-prot-norm-core.R** | message-changed (306 failing on both) | the pinned digests of the `lipid_ruv` and `metab_ruv` declaration results | **no** |
| **test-omics-summary-core.R** | message-changed (25 failing on both) | trx_main and trx_managed: the digest of `integration/transcriptomics_normalized_data.tsv` | **no** |
| **test-omic-summary-module.R** | message-changed (28 failing on both) | the same file digest, trx_main and trx_managed | **no** |

Every other file: `same`.

**S1 units, criterion (vi)** (`s1u/`). Passed. Base: 26 of 84 rows change calls (up to 22; max abs(delta t) 1.247), the
same as KIT-172 T1b. Head: 0 of 84 (max abs(delta t) 2.9e-11, max abs(delta logFC) 1.1e-10).

**S1 synthetic** (`s1/`). Instrument defect, no verdict. The declared-path part (v) raised
`omics_ruv_correction_failed` on every `k185` object on both sides. It shared one tryCatch with the journey, so all 50
`k185` results were lost, and with them (iv) and (v). On law `level` (49 objects): (i), (ii) and (iii) passed in every
row. At cutoff 0, 36 of 49 objects changed at least one call; at cutoff 0.5, none did. The correlation filter kept the
same samples in all 49. Object 16 (`level`) failed on the head only, downstream of RUV:
`Assay 'LCMS_Neg' has sample column(s) missing from design: A_6`. The correlation filter removed a sample that was
missing in one assay only (see the unrelated finding below).

## Inferred

- The three unapproved moves are the decided change seen through stored references the ticket did not list. The
  lipid/metab RUV result digests and the transcriptomics export file digest move because the values were restored.
  All three files already fail on runners on base, so only the messages change.
- The fix is unit-invariant on complete data: (vi), E6.

## Acceptance criteria

- A1: met (G0 plus R1 for E5).
- A2: met. R1 failed at the stated values; R2 head passes.
- A3: met (R2: (vi); R2c: (i) to (v), with (a) and (b) reported).
- A4: met for every file except three: `test-prot-norm-core.R`, `test-omics-summary-core.R` and `test-omic-summary-module.R` were accepted by the owner by file on 2026-11-25 (deferred to the end-of-project re-record). The edited tests pass (R2c).
- A5: met except one file outside the lists, test-gui-entry-transcriptomics.R (R3 above), accepted by the owner by file on 2026-11-25 (VAL-236); E5 of the new file was not reached in R3 (switch not set) and passes in R2 and R2c.
- A6: NEWS, ledger row and M8's status are done; `man/`, `NAMESPACE` and `DESCRIPTION` are unchanged.

## Coordinator rulings of 2026-11-25 and what was done

### The three files not in the owner's approved list (to the owner by file; the merge waits; nothing edited or re-frozen)

| File | Leaves whose message changed (R2, run 45510911807) | Why the decided fix moves them | Base |
|---|---|---|---|
| `test-prot-norm-core.R` | test "the lipid and metab RUV declaration give the recorded result through the moved owner stripper": `expect_identical(found, pinned[[key]], info = key)` (`:73`) for `key = "lipid_ruv"` and `key = "metab_ruv"`; the actual digest moves from 4ed0ad71... to 9ada5791... and from 5065333a... to cc4276ca... (2 leaves) | the digest is taken over the legacy `lipid_ruv` / `metab_ruv` declaration's result, whose corrected assays now carry restored feature levels (M8, OD-0) | already fails on base (306 failing on both; the pinned values are laptop-recorded) |
| `test-omics-summary-core.R` | tests "(W3, core side) trx_main ..." and "(W3, core side) trx_managed ...": the `save` step's `files` capture, leaf `integration/transcriptomics_normalized_data.tsv` (digest e0e8d896... on base, 7b9d35db... expected; head's value differs from both) (2 leaves) | the file is the transcriptomics RUVIII export, which under OD-2 T-a' keeps each gene's TMM logCPM level, so its bytes change | already fails on base (25 failing on both) |
| `test-omic-summary-module.R` | tests "(b) trx_main ..." and "(b) trx_managed ...": the `save` step's `files` capture, the same leaf (2 leaves) | as above (T-a') | already fails on base (28 failing on both) |

**Owner's decision, 2026-11-25 ("yes to all", item 3 of the coordinator's list): accepted by file.** Nothing is edited for them now; they join the end-of-project re-record list beside the KIT-073 and KIT-094 fixtures. The merge no longer waits on them.
Correction to E-3a (the KIT-225 shared-file list): `test-prot-norm-core.R` **does** change under KIT-224 (the two
leaves above). KIT-224's declared differences in the shared files are therefore `test-omic-norm-module.R` (12
`provenance` fields) and `test-prot-norm-core.R` (`lipid_ruv`, `metab_ruv` result digests); `test-kit129-run-path.R`
stays `same`.

### Merge of staging 67569959 (merge commit ac312e59)

Three files conflicted; everything else merged without conflict.
- `R/utils_workflow_analysis_record.R`: in `workflowReplayRuleMarkers`, KIT-202's row `intensity_filter_samples` and
  KIT-224's row `ruv_feature_levels` are both kept, KIT-224's last. The roxygen header names both.
- The `fixed` lines did not clash. KIT-202 fixes `intensity_filter_samples = "design"` on the intensity-filter
  declaration (`R/utils_workflow_omics_declarations.R:204` after the merge), while KIT-224 fixes
  `ruv_feature_levels = "restored"` on `omicsRuvOperation` (`:444`). Both behaviours are kept unchanged.
- `NEWS.md`: both entries are kept under `### Scientific changes`, KIT-224's first.
- `docs/methods-decision-ledger.md`: staging's rows to VAL-210 are kept. KIT-224's row stays `VAL-260` after them
  and now cites VAL-206, the coordinator's record of the owner's KIT-224 decisions.

### S1 instrument rerun 1 of 3 (VAL-54)

**Failure.** In R2, the declared-path part (v) raised `omics_ruv_correction_failed` on both sides. It ran on a KIT-172
object, whose design carries `tech_rep_group`. The cause is not established. It shared the journey's error trap, so it
discarded every `k185` journey.

**Fix** (commit 7110c7ab, script mechanics only):
- The declared object is built as test-kit224's G2 builds it.
- The declared part runs in its own error trap.
- A failure after the RUV step is recorded in `downstream-errors.tsv`, and the RUV state is kept.

No expected value, threshold, criterion or stop rule changed. The run binds the merged head 7110c7ab against staging
67569959 (R2c, run 45521043460).

### R2c (run 45521043460): edited tests and S1 instrument rerun 1, on 7110c7ab against staging 67569959

- **Edited stored-reference tests.** Every one passes on the head.
  - `test-omic-norm-module.R` is `same`, 639 passed and 0 failed on both sides, after the re-dump.
  - The new file passes 69 of 69.
  - The two characterization files, `test-kit127`, `test-integration-input-transcriptomics.R`,
    `test-workflow-trx-declarations.R` and `test-module-ci-trx-summary.R` gain passes. Their remaining failures are
    laptop-pinned digests that already fail on base: `trxDeclSha` products, and `kit127LegacyDigests`.
- **S1 synthetic** (`s1-r1/`; the first attempt's tables stay in `s1/`). All 748 criterion rows pass: (i) 200, (ii) 198,
  (iii) 200, (iv) 50, (v) 100. Object 16 under the `level` law gives the head-side downstream error described below.
  It is reported in `downstream-errors.tsv`; its (i) and (iii) rows pass.
- **Calls changed at q < 0.05** (cutoff 0): `k185` 31 of 50 objects, `level` 36 of 49. No object changed a call at
  cutoff 0.5, as expected: `treat` has no trend. The correlation filter kept the same samples in every judged object.
- **Rates (b)**, Monte Carlo 95 % intervals of the paired difference, head minus base:

  | Law | FDP | Power | Null share with any call |
  |---|---|---|---|
  | `k185` | +0.0003 (-0.003 to 0.004) | +0.0007 (-0.001 to 0.003) | 0 |
  | `level` | -0.0003 (-0.007 to 0.007) | -0.0003 (-0.004 to 0.004) | -0.04 (-0.12 to 0.04) |

  On these synthetic objects the repair changes which features are called, but the false discovery proportion and power
  stay the same within Monte Carlo error.
- **The legacy route itself, both sides (descriptive; not this ticket's; synthetic, not real data).** Design: two groups of six runs, two assays of 300 features, one unwanted factor orthogonal to the groups, feature variances from KIT-185's law (scaled inverse chi-square, 4 df) or that law with the SD doubling every 4 log2 below 20, 5 % value-dependent missingness in the even objects; route intensity filter, log2, cyclic loess, legacy RUV-III-C Automatic (percentage 1 to 20, no upper bound on k), Pearson filter at 0.5, limma on the corrected values with 10 residual df, eBayes(trend, robust), q < 0.05 at cutoff 0. Automatic chose k = 1 to 7 (`k185`) and 1 to 10 (`level`) per assay; some assay had k ≥ 3 in 30 of 50 and 39 of 49 objects. Null objects (complete null: the false discovery rate is the probability of any call): a call in 15 of 25 (0.60; Monte Carlo 95 % 0.40 to 0.80) and 16 of 24 (0.67; 0.47 to 0.86) against a nominal 0.05. Effect objects: mean false discovery proportion 0.36 (0.28 to 0.44) and 0.42 (0.36 to 0.49), power 0.91 and 0.88. By the chosen k (observational strata, small): null objects with k ≤ 2 in both assays had a call in 5 of 14 and 2 of 7 (mean under one call), with k ≥ 3 in some assay in 10 of 11 and 14 of 17 (mean 20 and 27 calls; 217 in one object with k = 10, the residual rank); effect objects' FDP 0.33 (6) and 0.22 (3) at k ≤ 2, 0.37 (19) and 0.45 (22) at k ≥ 3, and 0.24 (3) and 0.26 (2) with k = 1 in both assays, the true number of factors. The unbounded k (VAL-182) accounts for most of the null objects' false calls but not for the false discovery proportion, which stays at 4 to 7 times nominal at k ≤ 2. The cause of that remainder is not established by S1; testing corrected values with full residual df (VAL-53 (1)) is one candidate, not shown. The repair leaves all of this unchanged (paired differences above). To KIT-188 (truth-referenced calibration of the default route: arm LA against LCm at equal k separates the residual-df part; report FDP by chosen k), KIT-190 (a bound reduces the null objects' false calls but on S1 would not bring the FDP near nominal) and KIT-210 (k ≥ 3 with missing values gave the highest FDP, 0.43 and 0.49).

  **This ticket's repair is a correctness and unit-invariance fix, not a calibration fix**: it removes the unit
  dependence and leaves the route's error rates as they were.

### R3 (run 45573139074): full suite, 2871125b against staging 67569959

R3 was sent before staging moved to dd769671 and is kept as sent.

**Overall result**
- R CMD check: identical on both sides (0 findings only on head).
- 468 files `same`; the new file is `only-on-head` with 65 passed and 1 failed: E5 stopped at the KIT-040 corpus manifest check (`trxDeclReference`, `helper-workflow-trx-declarations.R:34`) because R3 ran without `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` (R1, R2 and R2c set it), so its assertions were not reached (they pass in R2 and R2c, 69 of 69). For the same reason R3 did not reach the assertions after `trxDeclReference` in the transcriptomics files that call it; R2 did.
- `test-omic-norm-module.R` and `test-kit129-run-path.R` are `same`.

**Files that changed, all within the approved lists and the three owner-accepted files**
- `count-changed`: the two characterization files, `test-kit127`, `test-integration-input-transcriptomics.R`. These are the accepted edits.
- `message-changed`: `test-module-ci-trx-summary.R` and `test-workflow-trx-declarations.R` (accepted edits; their remaining failures are the base's laptop-pinned digests).
- `message-changed`: `test-prot-norm-core.R`, `test-omics-summary-core.R`, `test-omic-summary-module.R`. These are the owner-accepted leaves, the same as in R2.
- `message-changed`: the KIT-073 fixture. Only `_ruv_corrected` and `_correlation_filtered` (plus registry, correlations and DA) differ. Every earlier state's message is identical on base and head, so R1-01 holds.
- `message-changed`: the KIT-094 fixture (`records`).

**Outside every list at R3 (stopped and reported; accepted by the owner by file, 2026-11-25, VAL-236):** `test-gui-entry-transcriptomics.R`, test "KIT-101 entry lane trx_screen: the KIT-040 count corpus, 12 samples, 300 genes", leaf `artifacts$integration/transcriptomics_normalized_data.tsv`.
- Four expectations change. The digest is 4645738d... on base and 7b1698df... on head; the stored value is 0f500346.... Every other line of the file's 13 failures is identical on both sides.
- The file already fails on base (16 passed, 13 failed on both sides).
- The cause is the same as for the two summary files: under T-a', the exported transcriptomics file keeps each gene's TMM logCPM level, so its bytes change.
- It was not in R2's scope (the ticket's grep does not match it), so it was first seen here.
- Owner's decision, evening of 2026-11-25 ("take all recommendations", item 3; VAL-236): the file may change as KIT-224's fix requires. Nothing is edited; it joins the end-of-project re-record list.

### The unrelated correlation-filter observation (for a ticket to reproduce)

**Seen in R2's S1** (run 45510911807), on the head side (12205bd4) only.

**Inputs.** Object r = 16 under the `level` variance law of `levels-control.R`:
- seed 224000016;
- 5 % of values missing, because r is even;
- no planted effect;
- two assays, `LCMS_Pos` and `LCMS_Neg`, 300 features each, runs A_1..A_6 and B_1..B_6.

**Route.**
1. Intensity filter (percentile 1, proportion 0.5).
2. log2 with offset 0.
3. Cyclic loess.
4. Legacy RUV Automatic (`kit172LegacyRuv("automatic")`).
5. Correlation filter (`group`, threshold 0.5).
6. DA.

**Result.** The DA step failed with `Assay 'LCMS_Neg' has sample column(s) missing from design: A_6`, raised at
`R/func_omics_design_shared.R:142`. The correlation filter (`filterSamplesByLipidCorrelationThreshold` through
`lipidomicsDeclarationBinding$filter_correlation`) removed A_6 from the design while `LCMS_Neg` kept that column.

**Reproduce.** Generate the inputs with the `values(16, "level")` function in `levels-control.R` and run the route
above. The data this repair restores made the filter drop a sample in one assay only; the defect is the filter's
handling of a sample dropped in one assay, not this change. Not fixed here.

## Hand-offs

KIT-210: if KIT-224 is in its bound commit, R3's post-processing restores lipid/metab levels and E10's lipid/metab DA columns are compared under `MultiScholaR.ruv_feature_levels = "shifted"` or dropped (ticket 2.5). KIT-193: VAL-96's statement of the legacy RUV-III-C option gains the restoration. KIT-129: U8 implemented here (E-5). KIT-190: its VAL-185 run states whether it binds a commit with KIT-224; the second of the two to merge resolves the marker list and `omicsRuvOperation`'s `fixed` by hand and re-dumps the 12 `provenance` fields.

## Note on the declared path (N2)

R2's S1 part (v) raised `omics_ruv_correction_failed` on both sides on a KIT-172-style object (`technical_replicate_id = "tech_rep_group"`, one run per block). This may be a correct refusal (such blocks leave no residual if the declared design takes technical replicates as the replicate structure); not checked. Recorded for the declared-path owner.

## Hand-off to KIT-188 (OD-5 (c))

The un-repaired legacy behaviour is reached by running the legacy route (Automatic or Manual) under
`options(MultiScholaR.ruv_feature_levels = "shifted")` (the replay marker's earlier token, `workflowReplayRuleMarkers$ruv_feature_levels$earlier`). Under it, `ruvIII_C_Varying` for lipid and metabolite data returns RUVIIIC's values exactly (E4 (ii); S1 (iv) once rerun). The option must be set in the process that computes the RUV step. Unset, or `"restored"`, gives the repaired path.

## Unrelated finding (not fixed)

Routed to the coordinator for its own ticket (I07), with the reproduction recipe in "The unrelated correlation-filter observation" above (S1 object r = 16, seed 224000016, `level` law, `levels-control.R`, run 45510911807 head side). It is on the shipped default lipid route.


The lipid correlation filter (`filterSamplesByLipidCorrelationThreshold`, through `binding$filter_correlation`) can remove
a sample from the design while one assay keeps its column. The object then fails validity at the DA
(`R/func_omics_design_shared.R:142`). This was seen once in S1, on the head; it does not come from this change.

## Not checked

KIT-224's evidence covers the legacy route without internal-standard normalisation: E3, E6 and S1 (vi) ran KIT-172 T1b journeys whose 'ITSD on' step did nothing (annotation `ISTD`, regex `^IS_`; AUDIT-KIT-242-plan-20261125 P01), so the five `IS_` standards stayed in the data as features; S1's synthetic objects run no ITSD step. With ITSD acting, the repaired route's unit invariance follows by algebra (median ITSD rescales each sample by a ratio of standard aggregates, in which a change of unit cancels), not by observation.


Windows; real data; R3 against staging dd769671 (R3 bound 67569959).

## For the change report

- **What was wrong.** Legacy lipidomics and metabolomics RUV-III-C (Automatic, the shipped default, and Manual) shifted
  every feature by a constant, s_j = colMeans(W)' alpha_j. That constant depends on the unit of measurement. limma-trend's
  prior therefore read shifted mean abundances, so moderated t, p, q and calls depended on the unit. Fold changes did not.
- **Who was affected.** Every lipid or metabolite analysis that applied RUV-III-C through Automatic (the default) or Manual, on `main` (metabolomics since d375d9d4, 2025-05-06; lipidomics since ce640728, 2026-01-18; at `main` 436cdbe8 no restoration after `corrected_matrix <- t(cln_mat)`, `R/func_lipid_s4_objects.R:2887`, `R/func_metab_s4_objects.R:2887`; DA defaults `eBayes_trend = TRUE`, `eBayes_robust = TRUE`, `R/func_lipid_da.R:138-139`, `:229`, `R/func_metab_da.R:139-140`, `:230`; RUV `automatic` selected by default, `R/mod_lipid_norm.R:149`, `R/mod_metab_norm.R:149`) and on staging. Staging only (`main` has no transcriptomics workflow): the transcriptomics RUVIII export file and the RUV-corrected integration input carried the same kind of shift; the transcriptomics DA did not.
- **Code.** Base cf8af245: `R/func_lipid_s4_ruv_methods.R:711`, `R/func_metab_s4_ruv_methods.R:711` (no restoration);
  `R/utils_workflow_trx_declarations.R:949-950`, `:1014`. Fix 12205bd4: the same lines plus `restoreRuvFeatureLevels`
  under the marker `ruv_feature_levels`; `R/func_omics_ruv_factors.R:361-364` (the declared path restores once);
  `trxRuviiiExportRestoredValues`.
- **Evidence.** KIT-172 T1b and S1 (vi): a change of unit changed calls in 26 of 84 complete-data rows before the fix
  and in 0 of 84 after (run 45510911807). E6 (seed 172102): 22 calls before, 0 after. (KIT-172 T1b journeys: legacy Manual, k 2, 5 % controls, complete data, the 'ITSD on' step a no-op, AUDIT-KIT-242-plan P01). The repair is not shown to change error rates (S1: paired changes at most 0.0007).
- **Not a calibration fix.** This repair is a correctness and unit-invariance fix, not a calibration fix: on S1 the legacy route was anticonservative on both sides of it (complete-null objects with a call 0.60 and 0.67, FDP 0.36 and 0.42 at nominal 0.05; cause not established; KIT-188, KIT-190, KIT-210, the owner).
- **Limits.** Missing values with zero-filling, and the Automatic k chosen from the uncentred cancor plot, still depend
  on the unit (KIT-185, KIT-188, KIT-210).
