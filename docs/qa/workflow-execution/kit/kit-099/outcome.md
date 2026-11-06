# KIT-099 outcome: transcriptomics at size (protocol v1.8), 2026-11-06

**Synthetic corpora only. Every run is recorded and uncounted.**
- Nothing here is claimed about real data, about a safe admission limit, or about the validity of the RUV procedure (the unwanted factor is synthetic).
- TRX-007's representative-workload qualification is not discharged.
- Runs were governed by the KIT plan's run ledger, not the superseded TRX allowance.
- No default, package code or behaviour changed. The default-variant decision is KIT-100's, by the user.

**Result.**
- **Sizes:** the four KIT-088 sizes are fixed at 60,000 raw genes (about 35% retained): small 24 samples, medium 100, large 200, ultra 600. 800 samples is the limit: it stopped at the 1,800 s step limit in normalisation.
- **Journeys:** the large journey and the screen small test are both accepted and exact.
- **Variant pairs at large:** both are accepted and exact.
  - Filter: `partitioned` 0.941 x `whole` on the kernel peak, a reduction.
  - TMM/logCPM: 0.979, not worse, no reduction.
- **Peaks at large:** the application stays under 1.9 GiB. The whole slice peaks at 5.01 GiB in normalisation, set by the composite QC figure, not by RUV.
- **Time:** time grows fastest in normalisation, which is what ends the curve at 800 samples.

Ticket revision 1.1. Declaration: [declaration.md](declaration.md), committed `9b884c27` before any run. The user's go: "you have my approval for large sizes." Build `kit_098` = `$P/core/browser-smoke.eoJzd6` (`c783e97a`), reused: no package change since.

## First step (KIT-098 hand-off)

**Change.** The `trx_*_science` products left `kitOperationIdentityProducts`, so they are judged directly. Test first ([test-first.txt](test-first.txt)): the harness test went from 694 pass / 3 fail to 697 / 0. Commit `2f32efcb`.

**098a and 098b re-verified offline** on reflinked copies, the runtime directories only read ([reverify.txt](reverify.txt)):
- both **accepted, exact**;
- 0 oracle mismatches;
- 12 `trx_*_science` entries now judged in each.

**Correction.** The first pass used KIT-092's `reverify.R` unchanged ([reverify-attempt1.txt](reverify-attempt1.txt)). That script predates the `member_dir` argument KIT-098 added to `kitSingleJourneyDecision`, so the decision read `worker_dispatch_missing`. KIT-099's copy passes the member directory. This was the step's one correction, and it was for a script defect.

## Members and runs (one attempt each, in the declared order)

| # | Run | Kind | Result | Time | Peak | Evidence |
|---|---|---|---|---|---|---|
| 1 | `kit-journey-099s` (screen, large configuration, vs `trx98a`) | single journey | **accepted, exact**; 6 of 6 worker steps with a completed journal | 208 s | slice anon 1.93 GiB, application 0.72 GiB | [stage](kit099-kit-journey-099s-stage.json) |
| 2 | forced failure | reused | KIT-098's `kit-lifecycle-trx-fail` on `kit_098` (`7e58122c`): failure observed, head unchanged; same build, fault and launch path (declaration 5.2) | | | `kit-098/` |
| 3 | `trx99sf`, `trx99sq` (screen pair populations) | in process | valid; states equal to `trx98a`'s (checked at bind) | 11 s, 14 s | | [oracle-runs.txt](oracle-runs.txt) |
| 4 | `kit-099a` (screen pairs) | 4 operation members | **accepted**; all computed (no cache hit), all exact | 14 to 16 s each | 1.07 to 1.08 GiB | [decision](operation-decision-kit-099a.json), [receipts](verification-099a.txt) |
| 5 | stepped trial | 6 in-process points | 24 to 600 inside the limits; 800 stopped at the step limit | 58 s to 2,160 s | RSS 1.68 to 5.75 GiB | [sizes.md](sizes.md), [size-trial.txt](size-trial.txt) |
| 6 | `trx99l` (large oracle) | in process | valid record, 9 states, history equal to `trx98a`'s | 260 s | RSS 3.0 GiB | [oracle-trx99l.txt](oracle-trx99l.txt) |
| 6 | `trx99lf`, `trx99lq` (large pair populations) | in process | valid; states equal to `trx99l`'s | 50 s, 81 s | | [oracle-runs.txt](oracle-runs.txt) |
| 6 | call profile | in process | complete (summary made offline, below) | 303 s | | [call-profile.txt](call-profile.txt) |
| 7 | `kit-099b` (large pairs) | 4 operation members | **accepted**; all computed, all exact | 23 to 34 s each | 1.32 to 1.68 GiB | [decision](operation-decision-kit-099b.json), [receipts](verification-099b.txt) |
| 8 | `kit-journey-099l` (large, vs `trx99l`) | single journey | **accepted, exact**; 6 of 6 worker steps | 706 s | slice anon 5.01 GiB, application 1.83 GiB, resident end 1.80 GiB | [stage](kit099-kit-journey-099l-stage.json), [steps](journey-steps.txt) |

**Host:** [power-state.txt](power-state.txt). Every run was on mains, in `low-power`, with 24 to 36 GiB available and the sleep lock held. Host swap use (6.7 to 8.8 GiB) came from other work. No member's slice used swap.

**Machine time:** about 2 h 40 min, from the first member at 21:12 to the last at 23:52, plus 45 min of binds. The declaration expected about 6 h; the trial ran faster than estimated.

## Sizes (KIT-088 hand-off)

| Size | Corpus key | Samples | Raw genes | Retained genes | Raw cells | Retained cells | Trial point | Finished | Manifest sha256 |
|---|---|---|---|---|---|---|---|---|---|
| small | `trx_60k_024` | 24 | 60,000 | 21,079 | 1,440,000 | 505,896 | 24 | yes (58 s) | `651f0a4b15f7…` |
| medium | `trx_60k_100` | 100 | 60,000 | 21,133 | 6,000,000 | 2,113,300 | 100 | yes (132 s) | `6d1de529daaf…` |
| large | `trx_60k_200` | 200 | 60,000 | 21,247 | 12,000,000 | 4,249,400 | 200 | yes (259 s) | `c233aaa00786…` |
| ultra | `trx_60k_600` | 600 | 60,000 | 21,280 | 36,000,000 | 12,768,000 | 600 | yes (2,056 s) | `eb5e8e7746c3…` |
| limit | `trx_60k_800` | 800 | 60,000 | 21,273 | 48,000,000 | 17,018,400 | 800 | **no**: step limit in normalisation | `d767df400708…` |

The 400 point finished too (780 s). Full sha256s and recipes are in the bindings (`corpora/trx_60k_*`).

**Recipe:**
- `kit_omics_corpus.R` `trx_screen`, seed 20261106, low-count fraction 0.6;
- two groups, two technical replicates, three batches;
- biological samples per group 6, 25, 50 and 150 for the four sizes (100 and 200 for the 400 and 800 points).

**Comparison on KIT-088's axis:** lipid and metab large are 8.0 M cells (2 x 20,000 x 200). The transcriptomics large has 12.0 M raw and 4.25 M retained cells at the same 200 samples.

**For KIT-088:**
- No `kit_086` control exists for this route, so transcriptomics has a single-build curve.
- Count corpora carry no missing values.
- A KIT-099 member stands as KIT-088 evidence only where build (`kit_098`), protocol (v1.8, sha in the bindings) and corpus hash all match.
- Only the large size was measured here as journeys. Small, medium and ultra exist as trial points, and a trial point is not a measurement.

## Variant pairs (KIT-100 hand-off)

**`kitOperationDecision`, verbatim** ([099b](operation-decision-kit-099b.json)):
- `status: accepted`, `failing_rule: null`, `promotion_authorized: false`.
- Ratios are candidate (`partitioned`) over control (`whole`), one run each, control first.

| Stage | Operation | Kernel peak (control, candidate) | Peak ratio | Slice anon ratio | Elapsed ratio | Not worse (≤ 1.02) | Reduction (≤ 0.95) |
|---|---|---|---|---|---|---|---|
| `kit-099b` (large) | `trx_expression_filter` | 1.682, 1.583 GiB | **0.941** | 0.931 | 0.861 | yes | **yes** |
| `kit-099b` (large) | `trx_tmm_logcpm` | 1.343, 1.316 GiB | 0.979 | 0.972 | 0.818 | yes | no |
| `kit-099a` (screen) | `trx_expression_filter` | 1.079, 1.078 GiB | 0.999 | 1.000 | 1.001 | yes | no |
| `kit-099a` (screen) | `trx_tmm_logcpm` | 1.074, 1.072 GiB | 0.998 | 0.998 | 1.017 | yes | no |

**Exactness:**
- The filter members reproduce the oracle's filtered state: class, rows 21,247 at large, components, support and summary.
- The logCPM members are product-only. Each created exactly one `trx.tmm_logcpm` product whose content hash equals the oracle's (large `87306ead0f10`, screen `c88cf6fabfd1`); their head state is unchanged.
- Every member's output was a new entity, not one its population store held (`kitOperationCacheHit` FALSE).

**Filter pre-pass, reported separately** (F07; the call profile, in process, [call-profile.txt](call-profile.txt) part a). The members record no phases for this operation. Peaks are this process's RSS, its baseline included.

| Declaration | Variant | Pre-pass | Rest (batched) |
|---|---|---|---|
| filter | `whole` | 13.7 s, 1.11 GiB | 8.9 s, 1.15 GiB |
| filter | `partitioned` | 11.8 s, **1.22 GiB** | 2.5 s, 1.14 GiB |
| logCPM | `whole` | 5.6 s, 1.21 GiB | 7.5 s, 1.20 GiB |
| logCPM | `partitioned` | 5.9 s, **1.28 GiB** | 1.2 s, 1.15 GiB |

In the partitioned variants the pre-pass sets the peak. It reads the whole count matrix (`trxSourceInputs`) for the TMM-effective library sizes, which is the declared pre-pass. The batched phase stays below it.

**KIT-080's selection rule, for reference:** accepted with a reduction. On these numbers it would select the filter and not logCPM. It was written for real data, and these pairs are synthetic, so the choice is the user's in KIT-100. Time ratios may favour the candidate through the file cache.

## Per step and per call (F08)

**Large journey (`099l`, [journey-steps.txt](journey-steps.txt), [peak-attribution-099l.txt](peak-attribution-099l.txt), [worker-times.txt](worker-times.txt)), anonymous GiB and seconds:**

| Step | Application | Whole slice | Seconds | Worker operations in the step (seconds, slice peak) |
|---|---|---|---|---|
| import | 1.08 | 1.42 | 18 | |
| design | 1.08 | 1.43 | 77 | |
| expression filter | 1.20 | 3.01 | 62 | filter 25 s 2.34; QC figures (raw) 30 s 3.01 |
| QC finalise | **1.83** | 2.59 | 14 | filtering metrics 6 s 2.59 |
| normalisation | 1.68 | **5.01** | **316** | QC figures 29 s 3.71; logCPM 18 s 2.47; **grid 67 s** 2.72; **RUVs 45 s** 3.57; QC figures 25 + 25 s 3.43; RUV publication 10 s; **composite figure 91 s 5.01** |
| DA | 1.49 | 2.57 | 81 | DA contrast 55 s 2.57; tab results 11 s; DA publication 12 s |
| RUVIII export | 1.45 | 4.24 | 105 | export 12 s 2.45; QC figures 28 s 3.23; composite figure 62 s **4.24** |
| report | 1.49 | 2.24 | 11 | |
| session export | 1.44 | 1.80 | 0 | |

**What sets the peak:**
- **Whole slice:** the normalisation step's composite figure (`trx_composite_figure`, the normalisation comparison grid), at 5.01 GiB. The RUVIII export step's composite figure is next at 4.24 GiB, then the QC figures at 3.0 to 3.7 GiB and RUVs at 3.57 GiB. The RUV grid peaks at only 2.72 GiB.
- **Application:** 1.83 GiB at QC finalise, and 1.4 to 1.7 GiB from there on. Unlike lipid and metab at KIT-084's large (3.33 GiB), the report step does not raise the application (1.49 GiB).

**What sets the time:**
- The normalisation step, 316 of 706 s. Within it: the composite figure 91 s, the grid 67 s, the QC figures about 54 s and RUVs 45 s.
- Over the journey, the figures take the most worker time: composite figures 153 s and QC figures 136 s.

**Within the calls, at large** ([call-profile.txt](call-profile.txt) part b; `Rprof` over the in-process recorder):
- The grid is 58 of the 96 s of in-process normalisation:
  - `ruv_cancorplot` (`trxCancor`) 31 s: its full `svd` and the `cancor` loop over k = 1..n;
  - the preliminary QL ranking (`trxPrelimOrder`, `estimateDisp` and the QL fit) 26 s.
- RUVs 14 s, of which `RUVSeq::RUVs` 6 s.
- DA 49 s, almost all `estimateDisp` in the fit.
- The R heap stayed at or below 1.06 GiB throughout.

**Across sizes** (trial, in process, [size-trial.txt](size-trial.txt)). Normalisation sets the time at every size:

| Samples | Normalisation | DA | Peak RSS |
|---|---|---|---|
| 24 | 26 s | 11 s | 1.68 GiB |
| 100 | 57 s | 29 s | 2.48 GiB |
| 200 | 123 s | 51 s | 2.83 GiB |
| 400 | 499 s | 106 s | 4.14 GiB |
| 600 | 1,635 s | 139 s | 5.75 GiB |
| 800 | over 1,800 s (stopped) | | 5.30 GiB at the stop |

Normalisation grows about as n^2.5 from 200 to 600 samples, consistent with the grid's `cancor` loops (k = 1..n per member, twice per member, six members). The in-process trial makes no figures and no report. The call split was profiled at large only.

**KIT-091 hand-off.** The transcriptomics RUV grid is the ceiling: normalisation passes 1,800 s between 600 and 800 samples at about 21,000 retained genes, in process. In a journey the grid is its own worker step, so the managed limit may sit a little higher; that was not run. At large, the grid's `ruv_cancorplot` (full SVD and `cancor` over k = 1..n) and its preliminary QL fit are the calls.

## Harness, generator and bindings (test first)

- **[test-first.txt](test-first.txt):**
  - step 1: 694/3 → 697/0;
  - step 2: 703 pass with 3 errors → 719/0. Added:
    - the corpus low-count fraction, with the default corpus regenerating byte for byte (the KIT-040 manifest sha256);
    - the trx operation stage (two pairs, at most 1,800 s);
    - `kitOracleProductMismatches`;
    - `kitOperationCacheHit`.
- **`workflow_kit.R`:**
  - the identity set without `trx`;
  - in `verifyKitOperationMember`, a product-oracle branch for product-only operations, and the cache-hit check (recorded for every member; it refuses a transcriptomics member).
- **`kit_omics_corpus.R`:** the trx recipe's `low_count_fraction`, with a CLI fifth argument for `trx_screen`.
- **`bindings-generator.R`:**
  - `operation_pair` and `operation_stage` (KIT-084's `pair_084`/`stage_084` generalised, F16);
  - the six trx size corpora;
  - `099s`, `099l`, `099a`, `099b`;
  - the `trx99l` oracle and history;
  - the product oracles `trx98a_logcpm_product` and `trx99l_logcpm_product`;
  - the four pair populations, each bound only when its state's fingerprint equals the oracle's.
- **[bindings-diff.txt](bindings-diff.txt):** the existing keys are unchanged except the two harness file hashes. KIT-084's `084p`, `084q` and `084r` are byte-identical.
- **Dependents ([dependents.txt](dependents.txt)):** every test file that sources the harness or the generator passes with 0 fail and 0 error:
  - `test-trx-design.R` 70;
  - `test-trx-kernels.R` 373;
  - `test-trx-publication.R` 174;
  - `test-workflow-analysis-record.R` 239;
  - `test-workflow-lipid-declarations.R` 259;
  - `test-workflow-metab-declarations.R` 259;
  - `test-workflow-trx-declarations.R` 308;
  - `test-workflow-kit-harness.R` 719.
- **`R CMD check`** ([check.log](check.log), `--no-tests --no-manual --no-build-vignettes`, kit_098 library): **3 WARNINGs, 3 NOTEs**, the same sections as KIT-098's baseline.
- **Full suite checkpoint** ([suite.txt](suite.txt); one process per file, 8 in parallel, kit_098 library, after every member): **376 files, 33,314 pass, 0 fail, 2 error, 314 skip**.
  - The two errors are the known diagnostics: `test-workflow-perf101-diagnostic.R` and `test-workflow-ruv-pearson-diagnostic.R`.
  - The skips are KIT-095's 313 plus one in `test-e2e-trx-journey.R`, an e2e file KIT-097 added after that checkpoint. Every other file's skip count is unchanged.

## Deviations and notes

1. **Pair populations (declared in advance, declaration section 5).** They come from in-process records stopped at the input step, not from the full oracle's stores. The scientific product identity includes the build's software digest, and a member's store is a full copy of its population store, so a same-build member on the oracle's store would have read the oracle's own product. Every population equals the oracle's state by fingerprint, and every member's output was a new entity.
   - **Unrelated discovery, recorded and not investigated:** KIT-081's `081a` ran same-build members on its own build's recordings, so its members may have read cached products.
2. **The forced failure was reused** from KIT-098 (declared; same build, fault and launch path; fixed stage names).
3. **Call profile summary.** `call-profile.R`'s own summary (`summaryRprof(memory = "both")`) stopped on a sample line it could not parse, after the run had completed. The stored profile was summarised offline by [call-profile-summary.R](call-profile-summary.R), with no rerun. That script reports the time per call and the largest R heap while the call was on the stack, rather than `summaryRprof`'s allocation column.
4. **Peak attribution script.** Its first pass found no worker requests (Python's `glob` skips dot-directories), and it was fixed before any result was read. This is analysis only.
5. **The ledger row** was added to `tickets/kit/00-KIT-PLAN.md`. That directory is git-ignored, so the edit is local, as earlier plan edits were.
6. **No new install.** `kit_098` is the one recorded install: no `R/`, `DESCRIPTION`, `NAMESPACE`, `inst/` or `src/` change since `c783e97a`.
7. **The trial's step limit is a conservative proxy.** The in-process normalisation step holds three worker steps.

## Claims and limits

- **Path equivalence:** the journeys and pairs reproduce the same build's in-process path. This is not a cross-build claim.
- **No memory claim against a control:** none exists for this route. The pair ratios are recorded, uncounted and synthetic.
- **No real-data claim.** No safe-limit claim: 800 samples is where this synthetic corpus met the in-process step limit, no more.

## Hand-offs

- **KIT-088:** the size table above, and [sizes.md](sizes.md).
- **KIT-100:** the two large pairs and their decisions; the user decides per step.
- **KIT-091:** the transcriptomics RUV grid ceiling (normalisation passes 1,800 s between 600 and 800 samples), with the per-call split at large.

## Rollback

`git revert` of this ticket's commits, from `2f32efcb` through the evidence commits. Evidence is never deleted, and the corpora stay outside the repository.
