# KIT-201 outcome (route (a); ticket r1.7)

State: see the end of this file ("Current state").

## Commits

| Commit | Content |
|---|---|
| ed56c7a6 | tests first: T-A1 to T-A7, the rewritten legacy-contract tests, `m1-sim.R`, `m1/family.tsv` (K = 990); nothing under R/ |
| 2c399d0c | T-A5's protein table without row names (R1 instrument defect); nothing under R/ |
| fa1f9909 | the fix, `man/runTestsContrasts.Rd`, NEWS, two ledger rows |
| 39f4c2b4 | T-196 (guard; KIT-196 audit K196-I04) |
| b0d7828a | merge of staging dd769671 (KIT-196, KIT-200, KIT-202 and others): marker list keeps every row; NEWS and ledger kept both |
| 65da7eab | T-196 (b) (guard, proteomics, after KIT-200) |
| a87b01a6 | M1 declaration, alone: R0's blind outputs, `m1/precision.tsv`, hashes, seeds, rule (`m1/declaration.md`) |
| c25e3e1d | `skip-none.R` (for KIT-188's E2) |

## Observed

**R1 (run 45328062374, ed56c7a6, tests only, no base).** `r1/tests.tsv`, `r1/failures.txt`.
- T-A1: lipid and metab `length(unique(s2.prior))` = 1 (stated value); s2.prior differs from the trended robust one;
  `df.prior` length 1.
- T-A2: `length(unique(s2.prior))` = 1 (stated, asserted first); p-values differ from `eBayes(trend, robust)`
  (e.g. 0.1727 against 0.1078).
- T-A3: proteomics at 1 and both S4 methods at 0: `length(unique(s2.prior))` = 1.
- T-A4 (guards): pass.
- T-A5: `workflowReplayRuleMarkers$treat_moderation` is NULL (stated value); then an error in the test's own fixture
  (`column_to_rownames`: the protein table had row names) before any DA ran: an instrument defect (VAL-54), repaired in
  2c399d0c with no change to any expected value, re-run in R1-2.
- T-A6: the second assertion fails with the strings equal (stated).
- T-A7: `fdp` 0.385 (stated); the note lacks "10.1%" and the reworded sentence.
- `test-prot-ruv-factors-in-model.R`: the legacy model record lacks `treat_moderation` (the new field, 1 failure); the
  rewritten legacy test fails at `length(unique(s2.prior))` = 1 and unequal p-values for entry, entry_technical and kernel
  (6 failures). `test-omic-norm-ruv-design.R`: the classic Legacy figure 0.385 against 0.389 (1 failure).
- `m1-sim.R family` on the runner wrote a `family.tsv` identical to the committed one (diff empty).

**R1-2 (run 45329779064; base 2c399d0c, head fa1f9909; one job per side).** `r1-2/`.
- Base (tests-only, T-A5 repaired): T-A5 now runs through and fails at the marker NULL (stated), then at (b)'s record
  and replayed model lacking the marker, "(a) and (b) differ" (base computes one prior under both rules) and (c)'s model,
  for proteomics and lipidomics alike; its replay-identity expectations pass on base. The other base failures are R1's.
- Head (fix): every test passes: 60 expectations in the new file, 312 in `test-prot-ruv-factors-in-model.R`, 1,153 in
  `test-omic-norm-ruv-design.R`; `compare.tsv` verdict `fixed` for all three files. A1 is met for this commit; R2 repeats
  it on the merged base.

**R0 (run 45573413958; base dd769671, head 65da7eab; one scripts job).** Blind: only `m1/r0-checks.tsv`,
`m1/r0-sd.tsv` and the formula's `m1/precision.tsv` / `precision-cost.tsv` were written and read. All 17 checks pass on
the first attempt: every cell completes on both sides in L, M and P; the data-set hashes agree; `length(unique(s2.prior))`
is above 1 on head and 1 on base wherever TREAT runs; the summary has 330 rows per kernel; the forced data set (L, cell 1,
r 9003) is recorded failed with its message and stops its cell's summary. The formula set R between 200 and 1,000 per
cell (982 comparisons gated, 8 reported, not gated: FDP overall and the lowest third's null rate in two L and two M
cells at R = 1,000); projected M1 cost 11.63 core-hours, under the 180.0 cap.

**M1 (counted; runs 45577233409 L, 45577238407 M, 45577242529 P; head a87b01a6, base dd769671).** Declared alone at
a87b01a6 (`m1/declaration.md`). Every cell completed on both sides with matching data-set hashes; no stopped cell.
`m1/class-*.tsv`, `m1/summary-*.tsv`, `m1/reported-*.tsv`, `m1/rows-*.tsv.gz`.

| Kernel | Class | Gated | Pass | Imprecise | Fails materially | Reported, not gated | Falsifier: max |delta p| (V1, V2) |
|---|---|---|---|---|---|---|---|
| L | **no worse** | 326 | 326 | 0 | 0 | 4 (all pass) | 0.286 |
| M | **no worse** | 326 | 326 | 0 | 0 | 4 (all pass) | 0.286 |
| P | **no worse** | 330 | 330 | 0 | 0 | 0 | 0.454 |

Sizes (means over data sets, head against base; ranges over cells, thresholds and contrasts):
- Null rate, lowest third of the level: base 0.043 to 0.096 → head 0.043 to 0.058 (P); 0.043 to 0.081 → 0.041 to
  0.052 (L, M). Largest decrease 0.040 (P: cell 10, V1, VD20, 2 x 3, tau 0) and 0.031 (L, M: cell 10, tau 0.5).
- Null rate, highest third: base 0.030 to 0.051 → head 0.047 to 0.051 (P); 0.034 to 0.051 → 0.048 to 0.051 (L, M):
  head - base up to +0.020 (P), +0.016 (L, M), with head within 0.01 of 0.05 (passes by rule 2's second difference).
- Null rate among features with a missing value (VD20): base up to 0.081 → head at most 0.057 (P); 0.074 → 0.051 (L, M).
- Null rate overall: head 0.046 to 0.053 (P), 0.045 to 0.051 (L, M).
- FDP: head - base -0.013 to +0.007 (P), -0.002 to +0.003 (L, M); upper limits at most 0.008 and 0.004.
- Power (reported): head - base -0.013 to +0.030 (P), -0.009 to +0.030 (L, M).
- Calls made on one side only: mean per data set up to 17 (P), 3 (L, M); zero-df and untested counts identical on
  both sides.
- V0 (no trend; KIT-124 arm 27's law), P at 0: FDP head - base -0.0003 to +0.0066 (reported beside arm 27's
  -0.0003 to +0.0005; this simulation has 3 to 6 per group and missing values, arm 27 had 6 and 12).
- **L and M are one experiment (AUDIT-KIT-201-impl I04).** They received byte-identical data (the seed has no kernel
  term) and run code identical up to names; `rows-L.tsv.gz` and `rows-M.tsv.gz` differ only in `kernel` and `seconds`.
  Their two "no worse" classes are one result for lipidomics and metabolomics, not a replication.
- At threshold above 0 (L, M; P at 1) the null rate is read on the boundary-calibration sets: features whose true change
  equals the threshold, not unchanged features (I02).
- With 20 % abundance-dependent missing values, proteomics at threshold 0 stays above the 5 % target: lowest-third null
  rate up to 0.0575, with 46 of 75 intervals wholly above 0.05; FDP up to 0.057. Without an abundance trend (V0, with
  missing values) the new test is more liberal than the old: FDP head - base +0.004 to +0.007, intervals excluding 0
  (I03; cause inferred: the trend's covariate is the mean of observed values, biased by censoring, ticket 2.4).
- Falsifier of the premise (5.2): not met; V1 and V2 differ (max |delta p| 0.29 to 0.45).

**R2 (run 45577246802; base dd769671, head c25e3e1d; 3 shards per side, check).** `r2/`. Of 101 files: 95 `same`,
the new file `only-on-head` (89 expectations pass, T-196 and T-196 (b) included), `test-prot-ruv-factors-in-model.R`
`fixed` (the old digest test is gone; base fails it on runners), and five files move: `test-gui-entry-proteomics.R` (the
four proteomics lanes from the DA step on: 6.2 (iii)), `test-workflow-saved-da.R` (W1/W2 DIA DA captures: 6.2 (i),
(iii)), `test-omic-da-module.R` (`prot_main`, named in 6.2 (iii); `prot_load`, `prot_no_contrasts` and the KIT-112 W3
managed run: the DA manifest digest, moved by the model field, 6.2 (i)), and `test-omic-enrich-module.R` (the enrichment
module's frozen captures: their `da_results_data` is the proteomics DA manifest and their runs the enrichment of its
legacy-test results: downstream of 6.2 (i) and (iii); **not named in 6.2**, listed here for plan section 11). No
`compare-error`. `R CMD check`: findings identical on both sides (0 only on head).

**For KIT-188's E2 (VAL-132 (2)): Skip and declared `none` on c25e3e1d.** `r2/skip-none.tsv`: on the same matrices
(2,000 proteins, 3, 4 and 6 per group, with and without missing values), `runTestsContrasts` (Skip's test) and
`runProtRuvModelContrasts(exact = TRUE)` (`none`'s) give the same rows, the same missing pattern and the same columns;
largest differences logFC 3.9e-14, raw p 3.2e-13, q 3.6e-13: **within 1e-10**, not bit-identical. Kernel level only:
whether the two routes hand the DA the same matrix is not compared here (KIT-188's E2 does).

## Decisions and deviations

- **Workbooks (owner, VAL-209, 2026-11-25).** No R Markdown workbook ships; nothing under `inst/workbooks/` is edited.
  Route (a) edited none; the "proteomics workbooks' DA prose" items belong to route (b), not taken (ticket r1.7).

- **Marker location (coordinator's ruling; ticket r1.4 7.1).** The DA is not a record step; the marker is the DA model
  record's `treat_moderation = "as_main_test"`, written while the current rule is in force by
  `acceptWorkflowOmicsDaManifest` and `omicsDaOperation` (lipidomics, metabolomics; legacy and declared lipid/metab
  designs) and by the proteomics legacy route's `da_model` (`differentialAbundanceAnalysisHelper`, `ruv_model` NULL).
  The DA replay (`replayWorkflowAnalysisRecord`) sets `MultiScholaR.treat_moderation_rule` from it: absent means the
  earlier rule (`"constant"`).
- **Declared proteomics designs and DIA with limpa** record no marker. Their computation does not read the option
  (`runProtRuvModelContrasts` passes the settings itself; limpa runs no TREAT), so their replay is unchanged
  whatever the option holds (T-A4 (c), (d) assert it under both rules).
- **Not rewritten:** `test-workflow-lipid-declarations.R:49` and `test-workflow-metab-declarations.R:49` (no
  declaration's `fixed` changes); `test-prot-ruv-factors-in-model.R:332-333` (both sides of that identity are
  `runTestsContrasts`, so it holds on head unchanged).
- **Rewritten** (6.2 (a)(ii)): `test-prot-ruv-factors-in-model.R` legacy digests (now: the legacy DA at 0 equals
  `eBayes(trend, robust)`, to 1e-10, for entry, entry with technical replicates and kernel) and its `da_model`
  expectation (gains the marker); `test-omic-norm-ruv-design.R:17` (0.385 to 0.389).
- **Continued session (T-A5 (c)).** The operation identity carries the build (`context$software`), not the rule; in the
  application a new build never serves an earlier build's product. Under `load_all` the build string is equal across
  commits, so T-A5 (c) models the new build by changing `context$software` before the re-run.
- **KIT-190 at merge.** Both branches append a row to `workflowReplayRuleMarkers`: a textual conflict at the list's
  closing line; keep both rows. KIT-190's `completeWorkflowAction` change filters its own marker from a declaration's
  `fixed`; this ticket's marker is in no declaration's `fixed`, so the two do not interact.
- **Old-record fixtures.** One committed record of a proteomics legacy DA made before the change exists:
  `tests/testdata/go-propagation/replay-base/record.json` (DIA, `DA-contrast-output`, cutoff 0; its model has no
  `treat_moderation`). Its project store is not committed; KIT-173's T7(b) replays it only inside that ticket's script
  `t7-replay.R` (skipped in the suite). Under this change it would replay with the earlier rule, so its DA, and the
  enrichment built on it, are expected unchanged (inferred; not run here). No committed lipid/metab record with a cutoff
  above 0 exists. T-A5 builds its records in process under each rule.

- **Wording drafts for the owner** (VAL-41): `wording.md` (NEWS, ledger, DAS-015's wording, the Legacy row, the
  Methods, and the correction of KIT-177's Methods limitation sentence, AUDIT-KIT-201-impl I11).

## Moved stored references (R3, run 45592673648; base dd769671, head c25e3e1d)

Of 482 files: 475 `same`; `test-kit201-treat-moderation.R` only on head; `test-prot-ruv-factors-in-model.R` `fixed`
(rewritten here, 6.2 (a)(ii)); six files move. None is edited (D5); each goes to plan section 11.

| File | Reference | Expectations failing, base -> head | 6.2 |
|---|---|---|---|
| `tests/testthat/test-gui-entry-proteomics.R` | `tests/testdata/gui-entry/prot_dia.json`, `prot_lfq.json`, `prot_lfq_fragpipe.json`, `prot_tmt.json`, from the DA step on (DIA lane `message-changed`; the three protein lanes `new-failure`) | 4 -> 25 | (iii), named |
| `tests/testthat/test-workflow-saved-da.R` | `savedDaW1Frozen` (`:764` on), W1 (a) and W2 DIA cases | 171 -> 175 (`message-changed`) | (i), (iii), named |
| `tests/testthat/test-omic-da-module.R` | `omicDaW1Frozen[["prot_main"]]` (`:3366`): DA files and manifest | 53 -> 80 (with the next two rows) | (iii), named |
| same | `omicDaW1Frozen[["prot_load"]]` (`:3254`), `[["prot_no_contrasts"]]` (`:4052`), and KIT-112 W3's managed run compared with `prot_main`: the DA manifest digest (model record) | | (i), by category |
| `tests/testthat/test-omic-enrich-module.R` | `omicEnrichW1Frozen` (`:836` on): `da_results_data` (the proteomics DA manifest) and the enrichment of the legacy test's results, in gprofiler, kernel_gprofiler, organism_filter, refusals, restore (`new-failure`), clusterprofiler and kernel_clusterprofiler (`message-changed`) | 22 -> 59 | **not named**: downstream of (i) and (iii) |
| `tests/testthat/test-omics-enrich-core.R` | W2's core-side comparison with W1 (b)'s frozen run step (`omicEnrichW1Frozen`, read from the module file): `da_data` digest of the run step | 9 -> 20 | **not named**: downstream of (i) and (iii) |

Unchanged as predicted: `prot_dia_limpa`, `trx_screen` and the lipid and metabolite GUI-entry lanes (`same`).

**Stop rule (8): resolved by the implementation audit (AUDIT-KIT-201-impl-20261125, I01).** Every leaf that moves in
`test-omic-enrich-module.R` (40 leaves; 22 → 59 failing expectations, 3 of base's being a flaky leaf, I13) and
`test-omics-enrich-core.R` (11 leaves; 9 → 20) holds one of three values of the proteomics DA that their fixture runs
on the legacy route at a cutoff of 0: the DA long table's digest (`fbe4b6a9…` → `67bdf842…`, the value
`omicDaW1Frozen[["prot_main"]]` also pins, 6.2 (iii)), the DA manifest's digest (`da_results_data`, `00c45cf1…` →
`936933ed…`, 6.2 (i) and (iii)), and that digest in the enrichment kernel's call (`da_data=0656c115…` →
`8098947f…`). No enrichment result moved: result files, display tables, the services' queries and the result objects
are equal on base and head. The two files were missing from 6.2's list, not from its prediction; they go to plan
section 11 (D5) and none is edited.

Coordinator's ruling (2026-11-25): resolved, a declaration omission with no consequence; A3 and A5 are met.

## Figures measured on the earlier TREAT moderation (D9; labelled, none re-measured)

- KIT-124 S1 arms 1 and 3 and S1's stress cells (legacy proteomics test, TREAT at 0, one prior variance).
- KIT-177 R0 (run 45222712083; `kit-177/r0/effect.tsv`) and its Methods limitation.
- VAL-108's simulation (KIT-187 S2).
- KIT-185's legacy proteomics DA.
- KIT-197 (the pattern it explains was seen on the earlier test).
- The RUV panel's classic Legacy row before this ticket (KIT-124 arm 3: 0.385, 9.8 %); restated here from arm 28.

## For KIT-188: Skip against declared `none`

Kernel level, measured in R2 (`r2/skip-none.tsv`): see 'Observed' above (largest differences 3.9e-14 in logFC, 3.2e-13 in
p, 3.6e-13 in q; t and BH not compared, AUDIT-KIT-201-impl I07). The route-level check is KIT-188's E2 on its bound
commit.

**Gaps of `skip-none.R` (I07; not re-run).** It does not compare `t`; its BH column name is wrong (`fdr_value_bh`; the
table's column is `fdr_value_bh_adjustment`), so `max_bh` is NA in every row; it covers two-group designs at cutoff 0
only. The 1e-10 tolerance is the owner's, declared before the run (VAL-132 (2)). KIT-188's E2 compares logFC, t, p and
q at route level, with the tables' own column names; the kernel result does not stand in for it (hand-off note in
`KIT-188.md`).

## Not checked

- No R was run on this machine; every result above is a runner's.
- The managed worker path's handling of the marker beyond R2 and R3's worker tests (KIT-112 W2 moved as predicted).
- Skip against `none` through the two whole routes (only the kernels, `r2/skip-none.tsv`).
- The replay of the committed KIT-173 record (its store is not committed).
- Windows; real data.

## Current state

R0, M1 (no worse, all three kernels), R2 and R3 done. Stop resolved (I01); records repaired; ready for a merge train.

## Acceptance criteria

| # | State | Evidence |
|---|---|---|
| A1 | met | R1 and R1-2 (each failing-first value; T-A5 after one instrument repair); R2 (the new file passes on the merged head: 89 expectations, T-196 and T-196 (b) included) |
| A2 | met | R0 (17 of 17 checks, blind); declaration a87b01a6 committed alone before M1; M1: L, M, P no worse; sizes above |
| A3 | met (stop resolved, AUDIT-KIT-201-impl I01) | R2: no `compare-error`; the `new-failure` verdicts are the moved references above |
| A4 | met | R2 check: findings identical on base and head |
| A5 | met (stop resolved, AUDIT-KIT-201-impl I01) | R3: the list above, with file, line and size |
| A6 | met | NEWS under the unreleased heading; two `VAL-<merge>` rows (one closes DAS-016); no version change; `man/runTestsContrasts.Rd` |

## Core-hours

R1 1.7, R1-2 1.9, R0 1.5, M1 15.7 (L 6.0, M 5.9, P 3.9), R2 35.2, R3 82.8: 139.1, against 107.3 expected,
524.4 ceiling (812.4 cap). R2 and R3 ran above their expected figures (16.8 and 65.1), under their ceilings.

## For the change report

**What was wrong.** limma's `treat()` re-estimates the variances itself and, unless told otherwise, uses one prior
variance for all features with no robust step (limma 3.68.4 `R/treat.R:3`, `:25-38`). The package ran
`eBayes(trend = TRUE, robust = TRUE)` and then called `treat(fit, lfc)` without those settings. So the reported t, p
and q came from a single, non-robust prior, while the configuration and the reports said "trended and robust".

**Where, in the legacy pipeline (`main` 436cdbe8).**
- Proteomics: every analysis, because TREAT ran at the default cutoff of 0 (`R/func_prot_da.R:2573-2580`;
  `inst/config/config.ini:110-112`).
- Lipidomics and metabolomics: whenever a fold-change cutoff above 0 was set (`R/func_lipid_da.R:229-236`,
  `R/func_metab_da.R:230-237`).
- The report sentence: `inst/reports/*/*_report.rmd` (for example `lipidomics_report.rmd:865-871`).

**Where, in development builds before KIT-201 (staging dd769671).**
- The same calls: `R/func_prot_da_model_stats.R:333`, `:618`; `R/func_lipid_da_model.R:135`;
  `R/func_metab_da_core_engine.R:135`; and the lipid/metab exported S4 methods through `runTestsContrasts`.
- The analysis record said `eBayes_trend = TRUE, eBayes_robust = TRUE`. In proteomics it also said
  `treat_applied = FALSE` at a cutoff of 0, although TREAT ran. (`treat_applied` and the analysis record exist only in
  development builds; `main` has neither.)

**Who was affected.**
- Every proteomics analysis of DIA without limpa, label-free and TMT data with Skip or a legacy RUV option.
- Lipidomics and metabolomics analyses with a fold-change cutoff above 0, and the exported S4 methods at any cutoff.
- Not affected: DIA with limpa, the declared RUV designs and transcriptomics.

**What changed (KIT-201, 50bdfed0).**
- TREAT takes the moderated test's trend and robust settings (`R/func_prot_da_model_stats.R:607`, `:625-626`;
  `R/func_lipid_da_model.R:137-139`; `R/func_metab_da_core_engine.R:137-139`). At a cutoff of 0 the proteomics TREAT
  p-values now equal the moderated test's up to rounding.
- Analyses saved before the change replay as they were made: replay rule marker `treat_moderation`
  (`R/utils_workflow_analysis_record.R:665-693`, `:932-952`).
- In scripts the earlier statistics are obtained with `eBayes_trend = FALSE, eBayes_robust = FALSE`.

**Evidence.**
- M1 (runs 45577233409, 45577238407, 45577242529; simulated data only; sizes as in NEWS).
- R1 and R1-2 (the old computation shown by failing tests).
- R2 45577246802 and R3 45592673648 (only DA statistics and DA record digests of the affected routes move).

**Limits.** Synthetic data: normal residuals, independent features, 2 to 4 groups of 3 to 6, 500 or 2,000
features, 0 or 20 % abundance-dependent missing values, no RUV, no blocking. Real-data figures are KIT-160 to KIT-166's.
