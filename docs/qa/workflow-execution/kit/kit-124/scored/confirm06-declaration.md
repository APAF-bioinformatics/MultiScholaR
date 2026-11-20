# RUV options at association 0.6 (extra-hard cell S3): declaration

Status: **approved 2026-11-21, committed and pushed before any job of this measurement.** Authority: the user's
approval of 2026-11-21 ("measure the 0.6 scenario for every option, like for like, pre-declared, on GitHub runners";
then "go" on this sizing). Owning tickets: KIT-124 (simulation), KIT-129 (figures). Choices not carried over from
`confirm-round-declaration.md` (the confirmation round, decisions M34, M38, M40, M41) were marked "[for the user]" in
the draft; all six were taken as proposed (2026-11-21) as engineering choices: bind 44b5f219; add S1 arm 37; S2 in two
jobs of 250 runs; branch `sim-k124-06` from f67b510e; the reported, non-gating check of arms 1 and 2 (section 4); a
ceiling of 50 core-hours. The marks are kept below as the record of what was approved.

On this commit `.github/k124-sim/params.json` and `score.json` hold inert placeholders that both workflows refuse at
their first step, so that creating the branch starts nothing; the lifecycle and the counted run are later pushes.
Paths below without a folder are under `docs/qa/workflow-execution/kit/kit-124/`. `final-figures.tsv` is the copy on
`performance-v2` (e8a2ac67); this branch does not carry `kit-129/`, so the successor is added here as
`docs/qa/workflow-execution/kit/kit-129/final-figures-06.tsv`.

## 1. Purpose

Fill the "not measured" figure at association 0.6 in the in-app check (the false discovery proportion at a 5 % target
in the extra-hard cell S3; `ruvGroupAssociationConsequences()`, `R/func_general_ruv_factor_count.R` lines 459 to 500,
column `fdp_06`; rows `consequence_fdp_06` of `docs/qa/workflow-execution/kit/kit-129/final-figures.tsv`). Two of the 12
rows are measured (transcriptomics and lipid/metab recommended, arm 34, `confirm-p2x`); 10 are not.

**The 0.6 scenario.** Extra-hard cell **S3** of each simulator: the centre cell C0 (2 groups, 6 per group, 20 %
true effects, one unwanted factor of strength 3, 25 % missing where the simulator has missingness) with the factor's
sample scores built as `w = rho * g + sqrt(1 - rho^2) * e`, `g` the contrast's group indicator scaled to variance 1,
`e` exactly orthogonal to the group model, so the sample correlation of the factor with the groups is exactly
**rho = 0.6** (C0: 0; C5: 0.3). Defined in each simulator's header (`sim/s1_classic.R` cell index 10, `s2_limpa.R`
index 13, `s3_trx.R` index 10, `s4_lipid_metab.R` index 10). Read: mean FDP at the 5 % target, subset `full`, 500 runs.

## 2. Bindings

| Item | Value | Why |
|---|---|---|
| Package under test (`package_commit` = `tools_ref`) | **44b5f2195e355bb6e6efe1dae21905e30b81c26c** (`kit-129-integration`, the confirmation round's binding) | Like for like with every other figure in the table. `performance-v2` (e8a2ac67) contains it; `git diff 44b5f219 HEAD` over the arms' entry points and the DA/RUV files changes only debug printouts (KIT-146), publication staging (KIT-052), record/integration fields and the figures table itself, so no numeric change is expected, but HEAD also changes the library (MOFA2-related Imports). Binding HEAD would need its own identity proof; not proposed. **[for the user to approve]** |
| Internal packages | as pinned by that commit: fastlimpadpc d4f9cc1d, fastduckdb d7851150 (`DESCRIPTION` Remotes, `packaging/multischolar.lock.json`) | unchanged |
| Simulators, gate, scripts | branch `sim-k124-06` **[name for the user]** cut from `origin/sim-k124-confirm` f67b510e; only `.github/k124-sim/params*.json` (and later `score.json`) change; `sim/*.R`, `gate.R`, `check-*.R`, workflow unchanged | same instrument |
| Scoring | `aggregate.R` of `sim-k124-exact` 2ba5cf87, sha256 `089d852c1f3b08d70096b1218ef26550f393cffd213e4bdd9980453075fd2a66`, named in `score.json` before any result is read; workflow `k124-score.yml` | unchanged |
| Runners | GitHub-hosted `ubuntu-22.04`, 4 cores, `timeout_minutes` 300 (job limit 345); savings L1 and L2 on (`share: "count,sweep"`, proven identical, declaration section 7) | unchanged |

## 3. Cells, arms, runs, seeds

Cell **S3** only, runs **1 to 500** (the confirmation round's rule: 500 per extra-hard cell, as `confirm-p2x`). Seed
rule unchanged: data seed = 100000 + 10000 x cell index + run; arm a on its own stream (seed x 100 + a); the count's
seed as each simulator's E seed. Arms 1 and 2 always (gate, rules 3 to 5).

| Sim | Arms (beyond 1, 2) | Rows filled | Chunk |
|---|---|---|---|
| S1 classic | 34, 35, 40, **37** | classic recommended, supplied list, first-pass nulls; **37 = shipped "No RUV" (M30's test) for the classic "without correction at 0.6" field, also "not measured" today [for the user to approve]** | 500 (one job) |
| S2 DIA limpa | 34, 35, 40 | DIA recommended, supplied list, first-pass nulls | 250 (two jobs) **[for the user]** |
| S3 transcriptomics | 35, 40 | supplied list, first-pass nulls (34 and 1 already measured, `confirm-p2x`) | 500 |
| S4 lipid/metab | 35, 40 | supplied list (standards), first-pass nulls (one table for both omics) | 500 |

Launch file (counted), label `ruv06`:

```json
{ "label": "ruv06", "package": "installed",
  "package_commit": "44b5f2195e355bb6e6efe1dae21905e30b81c26c", "tools_ref": "44b5f2195e355bb6e6efe1dae21905e30b81c26c",
  "runs_core": 0, "runs_stress": 500, "cores": 4, "timeout_minutes": 300, "chunk": { "s2": 250 },
  "arms": { "s1": "1,2,34,35,37,40", "s2": "1,2,34,35,40", "s3": "1,2,35,40", "s4": "1,2,35,40" },
  "cells": { "s1": "S3", "s2": "S3", "s3": "S3", "s4": "S3" }, "share": "count,sweep" }
```

## 4. Gate and checks (carried over; nothing decides from them but as stated)

- Gate job unchanged: S4's S3 is gated against the committed laptop rows (`s4-runs.laptop-extra.tsv.gz`, arms 1, 2,
  runs 1 to 500). S1, S2 and S3's cell S3 have no committed reference: reported "no reference" (as S3's S3 in
  `confirm-p2x`). A gated cell that fails is not scored and is reported.
- Reported, not gating: arms 1 and 2 cell means (M19 values) against the earlier scored rows of the same seeds:
  S3 and S4 against `scored/confirm-gha/s3`, `s4` (`confirm-p2x`, same binding and runners: expected identical); S1
  against `final-20261118/s1` (laptop rows, cross-platform), S2 against `final-20261118/s2` (runner run 44201244793).
  A difference is reported to the user before any shipped arm's figure is read. **[for the user to approve]**

## 5. Attempts and budget

| Step | What | Size | Core-hours (expected / upper) | Decision |
|---|---|---|---|---|
| Lifecycle (one) | label `ruv06-lc`: cells C0 and S3, runs 1 and 2 (`runs_core` 2, `runs_stress` 2), the counted arms and settings, S2 `chunk` 1 (proves the range path); then one scoring of it (label `ruv06-lc-score`, proves the scoring path on these cells) | 10 jobs (S2: 4 ranges), about 12 min each (mostly the 6 to 8 min library build) | about 6 to 8 | Pass = every simulator job success, every shipped arm `ok` or a classed refusal, nothing crashed, artifacts uploaded after the output scan, gate status 0 (C0 gated; S3 "no reference" except S4), the scoring run writes `k124-cells.tsv` with the S3 rows, and the C0 rows (runs 1, 2) of arms 1, 2, 35, 37, 40 (S1), 1, 2, 35, 40 (S2, S3, S4) equal `confirm-p2c`'s rows of the same seeds and arms column for column (the smaller arm set does not shift the streams; read by a text comparison of the downloaded result tables, no R locally; S2 bias and rmse may differ at about 1e-15 between machines, as declared in section 7 of the confirmation round). Fail = stop and report; no counted launch |
| Counted (one) | label `ruv06` (section 3), launched only after the lifecycle passed | 5 simulator jobs + plan + gate | about 23 / 36 | scored once |
| Instrumentation correction (at most one) | rerun of only the jobs that failed for an instrument reason (library build, token, upload, assembly), same parameters | the failed jobs | within the ceiling | an unfavourable result, a refusal or a timeout is not an instrument failure and is not rerun; a timeout goes to the user (as section 11 of the confirmation round) |
| Scoring (one) | `k124-score.yml`, `score.json` label `confirm06-gha`, the counted run, sims s1 to s4 | 1 job | under 1 | |

**Ceiling: 50 core-hours in all** (lifecycle, counted, correction, scoring). Reaching it stops the work and is reported.
Runtime ceiling per job: 300 minutes (none is expected above about 2 hours).

## 6. Scoring and reading (identical to the confirmation round)

Pass rules 1 to 7 unchanged through `aggregate.R` 089d852c; extra-hard cells list failing rules without failing a
design; primary subset `full`; refusals scored as arm 1 (R-07) and reported with their counts. The figure per row:
`mean_fdp` of the shipped arm in cell S3, subset `full`, from `k124-cells.tsv`, with its Monte Carlo standard error.
No threshold, rule, workload, run count or seed is changed after any result is read.

## 7. What each outcome means

- **A row is filled:** its `fdp_06` (and, for classic, `uncorrected_06` from arm 37) becomes the measured value with
  its source. Nothing in `R/` changes as part of this run; the figure update of `ruvGroupAssociationConsequences()`
  and its tests is a separate, later change on the user's word.
- **A result that contradicts a current recommendation or app text** (for example: a recommended design above the
  rule-2 bound of 0.075 at 0.6 where the app says it holds; an alternative clearly below the recommended design; a
  figure that conflicts with a preset note or caution, which today quote 0.3 or older 0.6 figures) is reported to
  the user as a decision. Nothing is changed automatically.
- **A refusal share** (lipid/metab supplied list under M32; first-pass with no residual degree of freedom) is reported
  beside the figure, never hidden.
- **A cell not run** (timeout, failed job after the one correction) stays "not measured" and is reported.

## 8. Where results are recorded

- `docs/qa/workflow-execution/kit/kit-124/scored/confirm06-gha/` (layout of `confirm-gha/`: `s1`..`s4/k124-*.tsv`,
  `inventory.tsv`, `jobs.tsv`, `not-committed.sha256`, `gate/`); raw per-run tables not committed (sha256 only).
- `scored/confirm06-results.md`: the results in plain words, written once after the single scoring, with the
  refusal share beside each figure; this declaration is not edited after the lifecycle except to record run ids.
- `docs/qa/workflow-execution/kit/kit-129/final-figures-06.tsv`: the successor of `final-figures.tsv`, the same 156
  rows with only the 13 "not measured" rows of association 0.6 (10 `consequence_fdp_06`, 3 classic
  `consequence_uncorrected_06`) filled with value, source, arm and cell; every other row byte-identical.
- Methods ledger: one new row (next M number, proposed **M42**), recorded on `performance-v2` by the main session.

## 9. Risks

- **DEPS_READ_TOKEN:** used only by the library action at `tools_ref`; never printed; nothing built with it cached or
  uploaded; artifacts hold the simulator out dir only, scanned by `check-output.R` before upload. Result tables may be
  artifacts (user, 2026-11-18). An expired or rotated token fails every job at the library step: an instrument
  failure (one correction), not a result.
- **S2 depends on the internal fastlimpadpc** (pinned d4f9cc1d through the binding) for the DPC quantification; its
  bias and rmse are not bit-reproducible across runner machines (1e-15; rule-bearing values are). S2 is the largest
  and most machine-variable cost (61 to 108 s per run in Part 2), hence two ranges.
- **Lipid/metab arm 35 at 0.6:** if the estimated count reaches the 5 standards, M32 refuses the assay and the run is
  scored as no RUV; the figure is then partly a no-RUV figure and must carry the refusal share.
- **Gate coverage:** only S4's S3 is gated; the others rest on section 4's reported check.
- **Earlier artifacts expire** about 2026-12-18/19 (30 days): the lifecycle's C0 comparison with `confirm-p2c` and any
  S2 row comparison with run 44201244793 need them before then.
- **Runners:** public repository, 4-core hosted runners, 20 concurrent jobs, 6 h per job; the counted run uses 5.
