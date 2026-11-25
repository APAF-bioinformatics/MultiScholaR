# KIT-248 G0 and G0b: outcome (base observation only)

Gate G0 of KIT-248 r1.2 (section 10), the first part of the ticket. It is an observation of base and runs no package code
change. Branch `kit-248-g0` is scripts-only: it carries this directory and nothing under `R/`, and it is never merged.

- **Binding.** `v0.9-staging` at 16112917 plus the script, commit bfc76de9.
- **Run.** 45592410449 (scripts job, 22.8 job-minutes, ceiling 90). R 4.6.1, reference BLAS, ubuntu-22.04, limpa 1.4.2,
  mixOmics 6.36.0, `workflowOwnerCapable()` TRUE.
- **Files.** The script is `g0-base-values.R`. The predictions, written before any value was read, are in its header and
  in `g0-predictions.tsv`. The observed values are in `g0-base-values.tsv` (one row per case, route, path and
  observable). `g0-best-partner-*.tsv` is supplementary.
- **Instrument errors.** None.

## Step 0 (by reading, at 16112917)

**Ticket text.** r1.2 applies Re-check 1:
- R1-01: the replay's `pearson-skip` branch, and the (a)/(a′) branches, set provenance before the commit; E2 asserts it.
- R1-02: the line shifts are corrected.
- R1-03 and R1-04 are written in as notes.

**Line shifts.** Re-reading at 16112917 confirmed them:
- `R/mod_prot_norm_workflow_helpers.R` moved by +10: `:628`, `:693`, `:724`, and `:743-744` for the Skip's call.
- `R/utils_omic_gui_config.R:2517` and `:2744` are at those lines at 16112917; they were `:2500` and `:2727` at 3eeaf100.
- The replay's `legacy` abort is at `R/utils_workflow_analysis_record.R:809-812`, `pearson-selection` at `:867`, and the
  correlation commit at `:880-881`.
- `R/utils_workflow_dia_products.R:1810-1812` and `:1886-1897` and `R/utils_workflow_dia_state.R:724-737`, `:778` and
  `:799-803`, `:834-843` are as the ticket cites them.

**R2 grep** (the ticket's pattern, without helpers, at 16112917): 73 files. The list belongs at the implementation's branch
point, after KIT-244's merge, and is re-taken there.

**Stored references.**
- The `tests/testdata/gui-entry/prot_*.json` oracles hold entity paths only (`.../entities/pearson-selection-<ID>`). They
  hold no record text and no study-parameters text.
- There is no `_snaps` directory.
- The one stored programme record, `tests/testdata/go-propagation/replay-base/record.json`, has `correlation_filtered` as
  a `design_builder` precondition. It has no `legacy` correlation step.

## Fixture and cases (as declared in the script header)

**Classic owners.** The KIT-182 object: "small" recipe, 64 proteins × 12 samples, groups control/treated. It is committed
as `protein_replicate_filtered` with `design_builder` provenance on a native owner of workflow type LFQ, TMT or DIA (iq).
Normalisation is quantile with RUV skipped (`runDIAReferenceNormalization`).

**Limpa owner.** `protNormW1LimpaOwner`'s commits, then the pass-through. RUV covariates run before the DA.

**Cases.** Each runs on a fresh owner, followed by a DA, the record, the study parameters, and a replay into a
`kitPrepareReplayProject` copy:
- **B:** Skip.
- **A:** Apply at t, then Skip. Run in process, and with the Apply through the worker (`protNormW2Dispatch`).
- **A2:** Apply at t, normalise again, then Skip.

## Observed against predicted

**SCAN.** Thresholds 0.5, 0.9 and 0.95 each kept all 12 samples, on every route. At 0.99 the Apply raised "Samples in
protein data and design matrix must be the same" on every route (INFERRED: it removed every sample).
- So no listed threshold removes some but not all samples. **Stop rule 2 applies.** The (A) cases ran at t = 0.5, which
  removed nothing.
- Predicted: all kept at 0.5 and 0.9 on the classic routes, because r is near 0.97. Most or all removed at 0.99. 0.95 not
  predicted. **Held.**
- Supplementary: the best within-group partner correlations ranged 0.9737 to 0.9836 on the classic routes (identical on
  LFQ, TMT and DIA) and 0.9626 to 0.9721 on limpa.

**B (Skip, no Apply). Every route, held.**

| Observable | Predicted | Observed (LFQ, TMT, DIA, DIA_LIMPA) |
|---|---|---|
| `step_provenance$correlation_filtered$operation` | legacy | legacy (all 4) |
| record step operation | legacy | legacy (all 4) |
| `incomplete` | has "legacy step: correlation_filtered" | has it (all 4); the record is `valid` |
| record `outcomes` | skipped TRUE | `{"rows":64,"columns":13,"skipped":true}` (all 4) |
| study parameters, correlation block | none | `<none>` (all 4) |
| live config | threshold 0, skipped TRUE | `min_pearson_correlation_threshold=0;skipped=TRUE` (all 4) |
| live samples | all 12 | 12 (all 4) |
| replay | `workflow_replay_unsupported`, "operation legacy" | `error [workflow_replay_unsupported/rlang_error/error/condition] step "correlation_filtered" cannot be replayed: operation "legacy"` (all 4) |

**A (Apply at 0.5, then Skip). In process and through the worker, every route: held** (record observables; numbers under
stop rule 2).

| Observable | Predicted | Observed (all 8 route × path cells) |
|---|---|---|
| provenance operation | pearson-selection, threshold t, preceded_by pearson-pairs | `pearson-selection`, values `{"pairs":"<id>","threshold":0.5}`, preceded_by `pearson-pairs` |
| record step | pearson-selection, threshold t, pairs `preceded_by` | `pearson-selection`; `preceded_by` `pearson-pairs(group, tech_rep_remove_regex)` |
| record `outcomes` | skipped TRUE | `{"rows":64,"columns":13,"skipped":true}` |
| `incomplete` | no legacy correlation step | none listed |
| study parameters | a pearson-selection block with threshold t; nothing about the skip | `[4. correlation_filtered: pearson-selection] \| pearson-pairs.group … = group \| pearson-pairs.tech_rep_remove_regex [default] = pool \| pairs [input] = <id> \| threshold [no_shipped_default] = 0.5` (limpa: step 5) |
| live config / samples | threshold 0, skipped TRUE / 12 | the same / 12 |
| replay | completes; selection at t; the samples kept at t | completed; replayed provenance `pearson-selection`; replayed config `min_pearson_correlation_threshold=0.5`; 12 samples (0.5 kept all), equal to live |
| DA equal | TRUE when t removes none | TRUE (live and replay digests equal: classic 2f75288ca73ae1a2, limpa 14c0a413ec587404); limpa `ruv_design_equal` TRUE |

**A worker.** The managed completion writes `pearson-selection`, as INFERRED. Every value is identical to the in-process
A.

**A2 (Apply, normalise again, Skip): held, as A, on every route.**
- On the classic routes the second normalisation left the head at `normalized`, and the Skip then re-committed
  `correlation_filtered`, keeping the Apply's `pearson-selection` entry. This is the P07 route, observed.
- On DIA with limpa the second pass-through left the head at `correlation_filtered`: it commits no state. The Skip then
  replaced the selection on the head with the same outcome as A.

**Stop rule 1** (the record after Apply then Skip naming a Skip, or replay committing the unfiltered table): **not
triggered.**
- The record names `pearson-selection` on every route and path.
- Replay re-runs the selection: its committed config is `min_pearson_correlation_threshold=0.5`, not the Skip's.
- At t = 0.5 the selection's table happens to hold every sample.

**Stop rule 2: triggered by G0's listed thresholds; resolved by G0b** (coordinator's ruling, ticket r1.3). E3 and E6's numbers rows take their base values from G0b (below), not from G0.

## G0b (run 45597013061): the (A) cases at t_hi

**Binding.** Staging 16112917 plus the script at ad0622df (mode `g0b`). Scripts job 20.2 job-minutes. No instrument errors.
Predictions were written in the script header and `g0b-predictions.tsv` before any G0b value was read; observed values
are in `g0b-base-values.tsv`.

**Rule** (coordinator's ruling; KIT-244's audited G0 rule): per route, t_hi = (m + m2) / 2 over G0's best-partner table.

| Route | m (sample) | m2 | t_hi |
|---|---|---|---|
| LFQ, TMT, DIA (iq) | 0.97371953882200679 (run_0008) | 0.97431187393916308 | 0.97401570638058499 |
| DIA with limpa | 0.96257597326065381 (run_0001) | 0.96464941373759827 | 0.96361269349912604 |

**Observed against predicted.** Every prediction held, on every cell.

| Case, routes | Observable | Predicted | Observed |
|---|---|---|---|
| all | Apply at t_hi keeps | 11 of 12 (classic: run_0008 out; limpa: run_0001 out) | 11; classic without run_0008, limpa without run_0001 (all 12 cells). **Held** |
| A in process and worker; A2 classic (9 + 2 limpa A cells) | live `correlation_filtered` | 12 samples; config threshold 0, skipped TRUE | 12; `min_pearson_correlation_threshold=0;skipped=TRUE`. **Held** |
| same | provenance and record step | `pearson-selection`, threshold t_hi, `preceded_by` `pearson-pairs`; outcomes skipped TRUE | `pearson-selection`, values threshold = t_hi, outcomes `{"rows":64,"columns":13,"skipped":true}`. **Held** |
| same | study parameters | a pearson-selection block with threshold = t_hi; nothing about the skip | `[4. correlation_filtered: pearson-selection] ... threshold [no_shipped_default] = 0.97401570638058499` (limpa: step 5, 0.96361269349912604). **Held** |
| same | replay | completes; selection at t_hi; 11 samples | completed; replayed provenance `pearson-selection`; config threshold t_hi; 11 samples; `replay_removed_vs_live` = run_0008 (classic) / run_0001 (limpa). **Held** |
| same | DA, live against replay | differ (numeric difference predicted) | `da_equal` FALSE on all 11 cells. Classic: live 2f75288ca73ae1a2, replay 98a76f0e9f445f3b. Limpa: live 14c0a413ec587404, replay c29d808643a24874. **Held** |
| A limpa, in process and worker | RUV covariates design, live against replay | differs | `ruv_design_equal` FALSE. **Held** |
| A2 limpa (KIT-244's base sequence) | live, record, replay | live 11 (C_hi, run_0001 absent); record `pearson-selection` at t_hi; replay 11, equal to live; DA and RUV design equal; no numeric difference | live 11; record `pearson-selection` at t_hi, outcomes `{"rows":64,"columns":12,"skipped":true}`; replay 11, `replay_removed_vs_live` none; `da_equal` TRUE (c29d808643a24874 both); `ruv_design_equal` TRUE. **Held** |

**What it shows.** At base, on every proteomics route, in process and managed, a user who applies the filter (removing
one sample) and then skips it gets results on all 12 samples. The analysis record and the study parameters state the
filter at t_hi, and a replay of that record removes the sample again, so the replayed DA table differs from the one
the user reported (and on limpa the RUV covariates differ too). On limpa's "Apply, press again, Skip" the live run
keeps the filtered table (KIT-244's defect, not this ticket's), its config says skipped, and the replay agrees with it.

**Base values for the implementation's tests.** E3 (and the numbers part of E6): base = the replayed DA table differs
from the live one at t_hi (digests above; limpa RUV design differs). E2's base = 11 samples (run_0008 / run_0001
absent). E1, E7, E8, E12 base = `pearson-selection` at t_hi with the study-parameters threshold line.

## Inferred (not observed)

- G0 alone could not show (A)'s numeric consequence (its listed thresholds removed no sample or every sample); G0b
  shows it at t_hi on every route (above). Only one sample removed was observed; larger removals are inferred to
  behave the same way (the replay re-runs the selection at the recorded threshold).
- The classic routes give identical values on LFQ, TMT and DIA (iq). On this fixture the workflow type does not enter
  the correlation step or the DA.

## Unrelated discovery (for the coordinator; not fixed)

At 0.99 the Apply fails with "Samples in protein data and design matrix must be the same". By reading, this is the
all-samples-removed case: the helper keeps only samples in a within-group pair with r >= threshold
(`R/func_prot_qc_correlation_helpers.R:389` onward); with no sample column left `cleanDesignMatrix` returns the object
unchanged (`R/func_prot_s4_qc_filter_methods.R:460`), and the ProteinQuantitativeData validity check then fails
(`R/func_prot_s4_core.R:108`). The user sees a technical error rather than "this threshold removes every sample".

## Verdict

**G0 and G0b passed as observations.**
- Every prediction held in both runs. There were no instrument errors.
- The readings of (B) and (A) are confirmed on all four routes, in process and managed.
- Stop rule 1 is not triggered. Stop rule 2 was triggered by G0's listed thresholds and is resolved by G0b: (A)'s
  numeric difference between live and replay is shown on every route.
- (B) stays in scope, because it is not replayable at base.

Code waits for the owner's OD-1 and OD-3 and for KIT-244's merge.

## Not checked

- E10's base (KIT-244's sequence on a tree holding KIT-244).
- Any threshold other than G0's list and G0b's t_hi.
- The rendering of the supplement (`analysis_record_SI.rmd`).
- The GUI observers themselves: the runners were called directly, as the tests call them.
- FragPipe and Spectronaut imports: their label-free route was represented by the LFQ type.
- Windows; real data.

## Implementation (branch kit-248 from staging 828dd105, after KIT-244's merge; owner's decision VAL-238)

### R1 predictions (tests-only commit on the base code; written and committed before R1 is sent)

The tests-only commit adds `tests/testthat/test-kit248-correlation-skip-record.R` and the two declared list entries
(`tests/testthat/test-workflow-operation-registry.R:6-18`, `:97-117`; VAL-238 (ii), by file); nothing under `R/`.
Base = staging 828dd105 (KIT-244 merged). t_hi by the G0b rule, computed in the test: classic about 0.974016 (removes
run_0008), limpa about 0.963613 (removes run_0001).

| Row | Predicted on the base code (each fails at this value) |
|---|---|
| E1, E7 (LFQ, TMT, DIA) | entry `pearson-selection` (not `pearson-skip`); no fixed constant; one `preceded_by`; record step `pearson-selection` with `threshold`; no `skipped` parameter; study block `[4. correlation_filtered: pearson-selection]` with a threshold line |
| G2 (inside E1) | passes: config threshold 0, skipped TRUE; components = the normalised source's |
| G3 (inside E1) | after DA the entry is `pearson-selection` (fails; G3's base) |
| E2 | replay completes; replayed entry `pearson-selection`; 11 samples against live 12; table differs |
| E3 | replayed DA table differs from live (G0b: classic 2f75288c… live, 98a76f0e… replay) |
| E4 | entry and record step `legacy`; incomplete lists "legacy step: correlation_filtered"; no study block |
| E5 | replay refused: `workflow_replay_unsupported`, 'step "correlation_filtered" cannot be replayed: operation "legacy"' |
| E8 | as E1 (record part) |
| E12 | as E1 (record part); live 12 samples (passes) |
| E11 | first assertion: entry `legacy` (LFQ and limpa); the other assertions pass |
| E6 | as E1 (record part); replay entry `pearson-selection`; 11 samples; covariates design and DA differ |
| E10 (I06: base on a tree holding KIT-244) | live components = L7's, 12 samples (passes); record step `pearson-selection`; replay 11 samples, table differs; in process and worker |
| E9 | the constructed pre-fix record is `pearson-selection` with outcome skipped (passes); replay completes; entry `pearson-selection`; config threshold t_hi; 11 samples |
| E13 | the constructed pre-fix step is `legacy` with outcome skipped (passes); replay refused ("operation legacy"); the earlier-legacy refusal names `normalized` (passes) |
| G1 | passes |
| E10b (G row) | passes: record `pearson-selection` at 0.5; replay 12 samples equal to live |
| registry lists | the two exact-list expectations fail on base (no `pearson-skip`), declared |

### R1 (run 45645244865, tests-only 2630bf57 on staging 828dd105): observed against predicted

`r1/tests.tsv`, `r1/failures.txt`. New file: 90 passed, 96 failed, 1 skipped (E5's replay after its refusal); 18.5
job-minutes. **Every failure is at its predicted value; G1, G2 and E10b pass.**
- E1/E7 (LFQ, TMT, DIA): entry and record step `pearson-selection`, `threshold` present, no `skipped` parameter, study
  block `[4. correlation_filtered: pearson-selection]` with a threshold line; G3: `pearson-selection` after DA.
- E2/E3: replay completes, entry `pearson-selection`, 11 samples against live 12 (run_0008 removed), DA table differs.
- E4/E5: `legacy`, incomplete lists "legacy step: correlation_filtered", no study block; replay refused
  ('step "correlation_filtered" cannot be replayed: operation "legacy"').
- E8, E12: as E1. E11: entry `legacy` (LFQ and limpa), the other assertions pass.
- E6 (limpa): as E1; replay 11 samples (run_0001 removed); RUV covariates design and DA differ.
- **E10 (I06 re-derivation, on a tree holding KIT-244): live components = L7's, 12 samples (passes); record step
  `pearson-selection`; replay 11 samples, table differs; in process and through the worker.** As predicted in r1.4.
- E9: constructed pre-fix record `pearson-selection` with outcome skipped (passes); replay entry `pearson-selection`,
  config threshold t_hi only, 11 samples. E13: constructed `legacy` step with outcome skipped (passes); replay refused;
  the earlier-legacy refusal names `normalized` (passes).
- Registry: the two declared lists fail on base (`pearson-skip` absent) and also on `integration_mofa_model`, a
  difference that exists at staging independently of this ticket (see R2).

### The fix (29072595)

- `R/utils_workflow_dia_products.R`: `runDIAReferenceCorrelation` sets `pearson-skip` (cleared pending first, fixed
  skipped = TRUE) only when the Skip commits; new `diaPearsonSkipOperation()` (store unit, no parameters, no-op kernel).
- `R/utils_workflow_operations.R`: `pearson-skip` registered in the in-process loop beside `enrichment`; empty defaults.
- `R/utils_workflow_analysis_record.R`: replay maps a proteomics `correlation_filtered` step recorded as
  `pearson-selection` or `legacy` whose outcome says skipped to `pearson-skip` (OD-1 (a), (a′)); a `pearson-skip` branch
  commits the correlation's source unfiltered (`normalization_final %||%` the head) with the Skip's config, description
  and provenance (R1-01).
- NEWS (### Fixes), ledger `VAL-<merge>`, `wording.md`. `man/`, `NAMESPACE`, `DESCRIPTION` unchanged.

### R2 (run 45649121588; base 2630bf57, head 29072595; 80 files, 4 shards per side, check): observed

`r2/compare.tsv`, `r2/check-compare.txt`, `r2/tests*.tsv`; 793.4 job-minutes (over the ticket's R2 row ceiling of 600;
the scope recounted at 828dd105 is 80 files, about twice the files the ceiling was sized for).
- 76 files `same`.
- `test-kit248-correlation-skip-record.R`: fixed, 206 passed, 0 failed (every E and G row passes on head).
- `test-workflow-operation-registry.R`: message-changed. The declared `pearson-skip` lines of both exact lists are
  resolved on head (A3's declared edit); both expectations still fail on base and head, identically, on
  `integration_mofa_model` (registered, absent from the stored lists at staging; not this ticket's). Counts +1/+2/+1 in
  three loops over the registry (the new operation).
- `test-workflow-partitions.R`: count-changed, passed 95 -> 96 (its loop over every operation; no failure).
- **`test-prot-norm-core.R`: message-changed (not predicted).** KIT-111 W1 (a) and W2 compare each click's capture with
  frozen digests (`tests/testthat/test-prot-norm-core.R:505-517`, frozen values from `:559`, `dia_correlation_skip` at
  `:718`). Those frozen values already fail on runners on base and head (laptop-recorded); on head the
  `dia_correlation_skip` provenance digest and one record digest (and the provenance digests of the later cases on the
  same owner) take new values (`r2/prot-norm-core-w1-*.txt`, `-w2-*.txt`). This is the intended change (the Skip's
  step is now `pearson-skip`) reaching a stored reference the ticket did not declare.
- R CMD check: identical on both sides (0 errors, 3 warnings, 3 notes; nothing only on head).

**Stop rule 3 applies** ("any file other than A3's moving: stop and report"): `test-prot-norm-core.R` (frozen KIT-111
W1/W2 digests) and, benignly, `test-workflow-partitions.R` (+1 passed) moved beyond A3's declared pair. Work stopped
here; no correction run was sent. Decision needed (owner, by file, as VAL-246 did for the same file): whether the
`dia_correlation_skip` and the following cases' provenance and record digests in `test-prot-norm-core.R` may be
re-pinned from a runner dump; the partitions count change needs no decision but is recorded.

## For the change report (draft; final after the implementation audit)

- **What was wrong (development branch only; `main` has no record or replay).** Proteomics, every route (label-free,
  TMT, DIA with iq, DIA with limpa), in process and managed. "Skip Filtering & Proceed to DE" committed the unfiltered
  table but set no step provenance. A Skip with no Apply before it was recorded as `legacy`: the record was marked
  incomplete, the study parameters omitted the step and a replay refused the whole record (G0 45592410449, all four
  routes). An Apply followed by a Skip kept the Apply's entry (`R/utils_workflow_dia_state.R:838`, the DA's keep rule):
  the record and the study parameters said the filter was applied at its threshold, and a replay removed the sample
  again, so the replayed DA table differed from the reported one (G0b 45597013061: one sample of 12 at t_hi on every
  route; 11 of 11 cells differ; limpa RUV covariates differ too). Live results were correct throughout.
- **Code.** Base 828dd105: `R/utils_workflow_dia_products.R:1820-1821` (the Skip's commit, no provenance),
  `R/utils_workflow_analysis_record.R:856-859` (`legacy` refused), `:914-928` (`pearson-selection` re-filters). Fix
  29072595: the same function's commit block and `diaPearsonSkipOperation`; `R/utils_workflow_operations.R` (the
  registration); `R/utils_workflow_analysis_record.R` (the mapping before the entry lookup and the `pearson-skip` branch).
- **What changed in results.** No live result. The record and the study parameters now state the skip
  ('[n. correlation_filtered: pearson-skip]', 'skipped [fixed] = TRUE', no threshold); a replay keeps every sample and
  its DA equals the live one (R2: E1 to E13 pass on head, including limpa's covariates and DA and KIT-244's second-press
  sequence, E10). Records made before the fix replay as the skips they were (E9, E13).
- **Evidence and limits.** Runs G0 45592410449, G0b 45597013061, R1 45645244865, R2 45649121588. Synthetic "small"
  recipe (64 proteins x 12 samples); one sample removed at t_hi; no real data; Linux runners only.
- **Parity (R1-04).** Lipidomics and metabolomics skip steps still print nothing in the study parameters (KIT-254).
- **Overlap.** KIT-257 (branch `kit-257`, not merged; AUDIT-KIT-257-impl I01) changes which files a DIA-with-limpa second
  press removes; this ticket's E10/E10b press twice and assert the store and the record, not those files.

## Not checked (implementation)

- The full suite (the merge train's); `test-prot-norm-core.R` re-pinned values (awaiting the decision).
- The supplement's rendering of a `pearson-skip` step; the GUI observers themselves (runners called directly).
- KIT-257's merged code with this ticket; Windows; real data.
