# KIT-248 G0: outcome (base observation only)

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

**Stop rule 2: triggered.**
- (A)'s numbers rows (E3, and the numbers part of E6) keep their record assertions only.
- The change report must say "numbers not shown to differ on the fixture".

## Inferred (not observed)

- (A)'s numeric consequence, that a replay removes samples the live run kept, is not shown on this fixture. The
  mechanism is shown: the replayed state is the selection at t, not the skip, so its config differs from the live one.
  It would show with a t that removes a sample.
- The best-partner table suggests such thresholds exist: between 0.9737 and 0.9836 on the classic routes, and between
  0.9626 and 0.9721 on limpa. The ticket's threshold list does not include one, and G0 chose none (stop rule 2 as
  written).
- The classic routes give identical values on LFQ, TMT and DIA (iq). On this fixture the workflow type does not enter
  the correlation step or the DA.

## Verdict

**G0 passed as an observation.**
- Every prediction held. There were no instrument errors.
- The readings of (B) and (A) are confirmed on all four routes, in process and managed.
- Stop rule 1 is not triggered; stop rule 2 is triggered.
- (B) stays in scope, because it is not replayable at base.

Code waits for the owner's OD-1 and OD-3 and for KIT-244's merge.

## Not checked

- E10's base (KIT-244's sequence on a tree holding KIT-244).
- Any threshold outside the ticket's list.
- The rendering of the supplement (`analysis_record_SI.rmd`).
- The GUI observers themselves: the runners were called directly, as the tests call them.
- FragPipe and Spectronaut imports: their label-free route was represented by the LFQ type.
- Windows; real data.
