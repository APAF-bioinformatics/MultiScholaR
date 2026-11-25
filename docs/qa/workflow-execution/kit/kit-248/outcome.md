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
