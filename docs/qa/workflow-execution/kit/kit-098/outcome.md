# KIT-098 outcome: the transcriptomics kit route (protocol v1.8), 2026-11-06

Ticket revision 1.1. User decisions of 2026-11-06:
- "Approve v1.8 as proposed";
- "Second notebook reference", on the run-order question.

Plan: [plan.md](plan.md). Protocol: `protocol.md` section 10. Build: `kit_098` = `$P/core/browser-smoke.eoJzd6`, from `c783e97a`. Nothing in `R/` changed.

**Result.**
- **Lifecycle pair:** observed (completion and failure).
- **Journeys:** `098a` and `098b` accepted, exact against their in-process oracles, with every listed worker step on a completed worker journal.
- **Notebook check (both journeys):** every compared product from the journey store is `identical()` to the notebook reference in the GUI's run order. For 098b, the RUVIII product is identical to KIT-093's direct `ruv::RUVIII` call.
- **Attempts:** one correction pass was used, on the first lifecycle member, for a demonstrated driver defect.

## Members (recorded, uncounted; one attempt each, 1,800 s)

| Member | Attempt | Decision | Evidence |
|---|---|---|---|
| `kit:lifecycle:trx:complete` | 1 (`e1c131cd`) | **invalid** (`verifier_accepted`) | `driver_failed` wait in the design step (below); kept as `member.attempt1` and [kit098-kit-lifecycle-trx-attempt1-stage.json](kit098-kit-lifecycle-trx-attempt1-stage.json) |
| `kit:lifecycle:trx:complete` | correction rerun (`3bd9bbc3`) | **lifecycle_complete_observed** | 189 s; 0 oracle mismatches; 8 steps judged; history equal; `in_memory = FALSE`; populations complete; [stage](kit098-kit-lifecycle-trx-stage.json) |
| `kit:lifecycle:trx:fail` | 1 (`7e58122c`) | **lifecycle_failure_observed** | 56 s; `normalisation` failed after the guard injected after `qc_finalize`; **`head_unchanged = TRUE`**, head `trx_qc_complete`; `in_memory = FALSE`; [stage](kit098-kit-lifecycle-trx-fail-stage.json) |
| `kit:journey-098a:candidate:1` | 1 (`b5d2bae9`) | **accepted**, exact | 188 s, peak 1.95 GiB; 0 mismatches over 8 judged steps; 6 of 6 worker steps with a completed journal; [stage](kit098-kit-journey-098a-stage.json) |
| `kit:journey-098b:candidate:1` | 1 (`f0e4e611`) | **accepted**, exact | 188 s, peak 1.94 GiB; 0 mismatches over 8 judged steps; 6 of 6 worker steps; [stage](kit098-kit-journey-098b-stage.json) |

**Judged steps:** `design`, `expression_filter`, `qc_finalize`, `normalisation`, `da`, `ruviii_export`, `report`, `session_export`.

**History, by name:** `initial`, `trx_raw_data_s4`, `trx_expression_filtered`, `trx_qc_complete`, `trx_ruvs`, equal to each oracle's.

**Worker steps**, one request and one complete journal each in both journeys:
- `trx_expression_filter`;
- `trx_tmm_logcpm`;
- `trx_ruv_grid`;
- `trx_ruvs`;
- `trx_da_contrast`;
- `trx_ruviii_export`.

**Power and host**, per member: [power-state.txt](power-state.txt). All members ran on mains, with no other managed unit and 35 GiB available. The sleep lock was held while each member ran.

### The correction pass

**Failure.** Attempt 1 of the complete member timed out in the driver's wait for the design builder's technical-replicate summary (`diagnostic browser condition timed out`, design step). The driver expected the pair's runs in `design.tsv`'s order, `RPMI_03_r1, RPMI_03_r2`.

**Cause.** The builder lists them in the design matrix's row order (`formatLipidDesignTechRepSummary`). On this route that is the counts file's order, `RPMI_03_r2` first. The application had done the assignment; only the driver's expectation was wrong. On the lipid and metab routes the design order is `mixedsort`, equal to `design.tsv`'s, which is why the wait had always matched there.

**Fix (test first).** `kitTechRepSummaryMatches`, `0fd03bd6`, see [correction-test.txt](correction-test.txt):
- the driver waits for the pair's header;
- it then matches the samples and replicate numbers as sets.

The bindings were regenerated for the driver's hash (`3bd9bbc3`). The member was rerun through `run_workflow_kit.sh rerun`, which keeps attempt 1. No case binding changed.

## Run order and the notebook reference

**(a) is false** ([run-order.txt](run-order.txt)): aligning runs by name gives 13 of 21 kernel products identical.
- RUVs' SVD flips W_1's sign under the permutation.
- logCPM, newY and the edgeR table move in the last digits.

**(b) would need `R/`.** The lipid and metab builders use `mixedsort`, which is not the notebook's order either.

**The user chose a second notebook reference** ([reference-gui-order.txt](reference-gui-order.txt)):
- **Render.** The unchanged notebook (`44d5f643…`) through KIT-040's driver (`8d1a2181…`, 900 s ceiling):
  - one attempt, exit 0, 48 s, 95 products;
  - `git status` unchanged.
- **Input.** `trx_screen` with `design_matrix.tab`'s rows in the counts file's order:
  - the permuted design has sha256 `f8763aaf…`;
  - every other file is byte-identical.
- **Confirmed before use.** KIT-093's kernels in that order are identical on 21 of 21 products.
- **Stores.**
  - Full products: `$P/kit/trx-reference/reference-trx_screen-gui-order.rds` (`8ff1202a…`).
  - Fingerprint store: [reference-trx_screen_gui_order.rds](reference-trx_screen_gui_order.rds).
- KIT-040's reference is unchanged.

## Notebook comparison (from the journey stores; `identical()`, no alignment, no tolerance)

| Journey | Compared | Identical | Excluded with a reason | Unmapped | Files |
|---|---|---|---|---|---|
| 098a | 32 | 32 | 66 | 0 | [txt](notebook-compare-098a.txt), [tsv](notebook-compare-098a.tsv) |
| 098b | 34 | 34 | 66 | 0 | [txt](notebook-compare-098b.txt), [tsv](notebook-compare-098b.tsv) |

**Compared** (plan section 6):
- the filtered counts;
- the filter's effective library sizes, minimum sample size, median library size and its three stage counts;
- the TMM samples table and logCPM;
- the grid's order, k cap, members, CCA data, control sets, summary, best result and controls;
- best k;
- W;
- the RUVs counts and the corrected logCPM;
- the DA design, fit, formula, contrast, contrast matrix, `qlf`, its table and topTags;
- the RUVIII control mask;
- for 098a, newY and M against the notebook.

**098b's RUVIII rows:**
- the journey design's technical-replicate partition equals `design.tsv`'s;
- `newY` is `identical` to `ruv::RUVIII(Y = t(ref log_cpm_tmm), M = ruv::replicate.matrix(sets), ctl = ref mask, k = ref best_k)$newY`;
- `M` is `identical` to `ruv::replicate.matrix(sets)`;
- `m_source` is `technical_replicate`.

**Excluded, with their declared reasons:**
- the masks as vectors, raw logCPM, and the notebook intermediates (not stored by the journey);
- the renderer and PCA/RLE statistics outputs (display values outside the head's products; their input matrices are compared).

**Dry run.** The same comparator on the two recorder stores gave the same result: [notebook-compare-dryrun-trx98a.txt](notebook-compare-dryrun-trx98a.txt), [trx98b](notebook-compare-dryrun-trx98b.txt).

## Oracles and chain check

**Recordings.** `trx98a` (M = group) and `trx98b` (M = technical replicates), made with [oracle-record.R](oracle-record.R) on `kit_098`:
- 33 s each against the 1,800 s ceiling;
- records valid;
- settings read from the bindings;
- each records the `workflow_kit.R` sha256, R 4.6.1, reference BLAS/LAPACK 3.12.0, the variant per step (`whole`), protocol v1.8 and the history.

They were bound by `8930613a`. See [oracle-trx98a.txt](oracle-trx98a.txt) and [oracle-trx98b.txt](oracle-trx98b.txt).

**README chain check** ([chain-check.txt](chain-check.txt)):
- **Items 1-3 and 5:** KIT-095's registry and worker-harness tests on `kit_098`:
  - `test-workflow-operation-registry.R` 1007 pass;
  - `test-workflow-trx-declarations.R` 308 pass;
  - 0 fail, 0 skip in both.
- **Item 4:** not applicable (no operation member).
- **Item 6:** [chain-check.R](chain-check.R) on the bound bindings:
  - `kitCases` accepts all four trx stages;
  - both recordings re-fingerprint with 0 mismatches over 8 judged steps;
  - a positive control (one oracle product and the history altered) is caught.

**`trx_*_science` in `kitOperationIdentityProducts`** ([science-compare.txt](science-compare.txt)):
- **Not run-specific in process.** KIT-031's cause (run-specific bytes) is not shown between the two in-process recordings: the science products are equal.
- **Aliases only.** Each `trx_*_science` reference is the same entity as an attached `trx.*:<sha12>` product, which is judged (equal sha256).
- **Effect.** Listing them as identity products drops an alias only; no content goes unjudged.
- **Next build.** If the next build should judge them directly, that is a one-line regex change for the coordinator.

## Harness (test first)

**Test first.** [test-first.txt](test-first.txt): against the harness at `c4839d25` the new tests gave 3 failures and 6 errors. After `c783e97a`: `test-workflow-kit-harness.R` 688 pass. After the correction: 694 pass, 0 fail.

**Sites** (plan section 7):

| Site | Change |
|---|---|
| `kitJourneySteps("transcriptomics")` | the trx step list |
| `kitMemberSeconds` | `lifecycle_trx`, at most 1,800 s |
| `kitBindings` | carries `lifecycle_trx` |
| `kitCases` | the trx step-map guard; an unknown route refused (no DIA fallback); the trx lifecycle block; every trx case binds oracles, is never in memory, and has an oracle per state step and a RUVIII product (`oracle_missing`, also in the v1.7 refusal message) |
| `kitSelectorStage` | the trx lifecycle selectors |
| `kitLifecycleDecision` | the stage by route; an unknown route refused |
| `kitOperationIdentityProducts` | `trx_*_science` |
| `kitSingleJourneyDecision` | `worker_steps`, using `kitWorkerDispatches`, shared with `kitGuiMemberDecision` (`worker_dispatch_missing`) |
| `kitHistoryMismatches` | new: the history check |
| `verifyKitMember` | the history check; a journey mismatch raised as "oracle mismatch: …" |

**Driver:** `runKitTrxImportActions` and `runKitTrxStep`. `runKitOmicsImportActions` and `runKitOmicsStep` refuse other routes.

**Coordinator:**
- the stage regex;
- trx takes its own normalisation, DA and M settings;
- the member directory is passed to the single-journey decision.

**Runner:** usage text and stage regex.

**Generator:**
- the `trx_screen` corpus (its manifest asserted equal to KIT-040's reference input) and the second reference (asserted derived from it);
- `step_map_trx`;
- `lifecycle_trx`;
- keys `098a` and `098b`;
- oracle ingestion with `<label>_history`.

## Checks

- **Dependents** ([dependents.txt](dependents.txt), kit_098 library):

  | Test file | Pass |
  |---|---|
  | `test-workflow-lipid-declarations.R` | 259 |
  | `test-workflow-metab-declarations.R` | 259 |
  | `test-workflow-analysis-record.R` | 243 |
  | `test-workflow-kit-harness.R` | 688 |
  | `test-workflow-trx-declarations.R` | 308 |

  0 fail and 0 skip in each.

  The grep finds no lipid-or-metab binary left on another route's path. The one remaining `"lipid" else "metab"` (`workflow_kit.R:39`) sits inside the lipid/metab branch of `kitJourneySteps`.
- **Backward compatibility** ([bindings-diff.txt](bindings-diff.txt)):
  - **The first bind** changed only:
    - `builds/kit_098`, `corpora/trx_screen`, `harness/step_map_trx*`, `lifecycle_trx`, `tickets/098a`, `tickets/098b`;
    - the protocol hash;
    - the five harness file hashes.
  - **The later binds** added only the `trx98a_*`/`trx98b_*` populations, oracles and histories, then the driver's hash.
  - The existing route tests pass. Earlier receipts stand under their own protocol hash.
- **R CMD check** ([check.log](check.log), `--no-tests --no-manual --no-build-vignettes`, HEAD `bacb0138`, kit_098 library): **3 WARNINGs, 3 NOTEs**, the same sections as KIT-097's baseline.
- **Install:** `kit_098` from a pinned worktree at `c783e97a`. The RUVSeq/EDASeq family was carried from kit_097's library. The package code equals `kit_097`'s (`35f855bd`): no `R/` change since.

## Divergences and notes

- **Mismatch classification, all routes.** `verifyKitMember` now raises a journey oracle mismatch as "oracle mismatch: …". Before, the bare assertion read as a harness refusal, so a single journey's genuine mismatch was decided `invalid`, not `rejected`. This applies to every route; recorded decisions stand.
- **One GUI action, three worker steps.** The normalisation button runs `trx_tmm_logcpm`, `trx_ruv_grid` and `trx_ruvs`, so they are one journey step. Their products are judged as attached products, and every history state is judged by name (ticket's alternative).
- **Driver checks.** The driver checks completion texts only. It does not check the DA tab's displays as the lipid route's KIT-082 checks do. The interactive volcano was left off (KIT-097's defect).
- **Publication products.** In these members the head's products carried no additional publication products (0 additional), so every product the journey's head held was judged.
- **Runtime.** The members' runtime directories and the recorder stores are under `$P/kit/runtime/` (`kit-lifecycle-trx`, `-fail`, `kit-journey-098a`, `-098b`, `oracles/trx98a`, `trx98b`). The comparator read reflinked copies under `$P/kit/scratch-098-compare/`.

## Claims and limits

- **Journey against the same-build oracle:** path equivalence only.
- **Journey store against the notebook:** the independent check, on one synthetic corpus.
- **No memory claim, no cross-build claim, no real-data claim.**

## Hand-offs

**KIT-099:**
- the trx route;
- `lifecycle_trx`;
- the oracles `trx98a` and `trx98b`;
- the comparator;
- v1.8's operation-member rules;
- the driver's set-wise technical-replicate wait.

## Rollback

Revert this ticket's commits, from `9b54b6ab` to the evidence commits. Evidence is never deleted.
