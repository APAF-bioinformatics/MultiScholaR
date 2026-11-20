# KIT-150 declaration: memory and time of the integration module's on-disk rework

Fixed before any run (lifecycle, probe or counted). Follows KIT-150 revision 1.1 ("Revision 2026-11-20 (audit
repairs)"), with the user's decisions of 2026-11-21 written in below. Lane: branch `run-mofa-measure-k150`,
workflow `.github/workflows/k150-measure.yml`, scripts `.github/mofa-smoke/k150-*`. Measurement only: no package code
is changed.

## 1. Purpose

Measure, per mechanism, how much memory and time the integration module's rework (KIT-147, KIT-148, KIT-149) saves
against the merged skeleton, stage by stage, at fixed synthetic study sizes; and state separately what MOFA2 and
mofapy2 themselves need, which the package cannot reduce. The measurement is of the mechanisms; it is not of the
final default model (the feature-selection statistic is changing in another ticket, and the cap value used here is a
fixed 1,000).

## 2. Decisions since revision 1.1 (user, 2026-11-21)

- INT-22: the direct mofapy2 hand-off is being removed. No build measures it and no keep rule is decided here; the
  revision's "B-head, switch on" row and "The direct hand-off's keep rule" section are dropped. In builds B3 and C the
  internal switch `MultiScholaR.integration_mofa_handoff` is set to `run_mofa` (off) explicitly.
- INT-15: convergence mode `slow`, round limit 10,000. Timings are under `slow`. The uncapped "MOFA2 floor" rows keep
  the fixed 10 rounds revision 1.1 proposes, and report no time.
- Builds are the five commits of section 3, all with the same model settings passed explicitly (section 6).
- Feature selection: top variance, cap 1,000 per data set, for the capped claim.
- The synthetic studies are generated inside the job by `k150-stage.R` (seeded) and written through each build's own
  input writer (`writeIntegrationInput`: TSV for base, Parquet for the others). No data files are committed.

## 3. Bindings

| Build | Commit | Mechanism it isolates (against the row above) | Stage it changes |
|---|---|---|---|
| base | 2bc0c5f5adb5ecc56d2a4a07411c54d429a2fdb4 | none (merged skeleton: TSV inputs, whole matrices) | |
| A | 0f49c7502659526a7ee97f03c9fbc66aa9af6802 | Parquet inputs, streamed import (KIT-147) | input write, import |
| B1 | 953e0996eeae3a0e4118d8539ed761f6e413b6ae | `save_data = FALSE`, one model load, fewer copies (KIT-148 item 1) | training |
| B3 | 58ddd20b8044dcb21c8f239624daf05e19376d88 | product v2 (KIT-148 item 3; head of `kit-148-integration-handoff`, item 2's switch off) | training (product write), results read |
| C | 324e14b44626d25bc247fd88c2cee9084f114c4c | partitioned alignment (KIT-149, default variant `partitioned`, passed explicitly) | alignment (and training, which aligns first) |

"Final configuration" = C. "Before" for the final comparison = base.

Environment: GitHub-hosted `ubuntu-22.04` runner (4 vCPU, 16 GB; image version, CPU model and kernel recorded per job in
`runner.txt`); R 4.6.1 (Posit r-builds); the KIT-134 snapshot library (`packaging/multischolar.lock.json`, identical in
all five commits) installed from build C's checkout by `.github/actions/snapshot-library`; MOFA2 1.22.0 and its closure
from the same snapshot (`install-mofa-library.R`); the bundle's Python (python-build-standalone CPython 3.12.16, release
20261110; mofapy2 0.7.3, numpy 1.26.4, h5py 3.10.0; hashes checked by `packaging/linux/build.sh python`); R started
through the bundle's launcher (`packaging/linux/Rscript`, RETICULATE_PYTHON = the bundle's Python). Each build is
installed with `R CMD INSTALL` into its own library (`$RUNNER_TEMP/builds/<build>`), placed first on the library path
of each stage process; the stage refuses to run unless MultiScholaR loads from that library and the library's marker
names the build's commit. No thread variable is set (the in-process app's default: OpenBLAS uses the runner's cores).

## 4. Terms and how each is measured

| Term | Definition |
|---|---|
| Stage | One of: input write, import, alignment, training, results read (section 5). Each stage runs in a fresh `Rscript` process, in process (no worker), alone in its own memory cgroup (cgroup v2, `memory.max` 14.5 GiB, `memory.swap.max` 0). |
| Baseline | Taken in that process after loading MultiScholaR and MOFA2 and importing mofapy2 (and numpy) through reticulate (and, for the input write only, after the synthetic study is generated in memory), then `gc(full = TRUE)` and `gc(full = TRUE, reset = TRUE)`: R heap in use (Ncells + Vcells, MiB, from `gc()`) and `VmRSS`. Then `VmHWM` is reset to the current RSS (`/proc/self/clear_refs`, value 5); the reset is checked (VmHWM after reset at most VmRSS + 16 MiB). |
| R heap peak | `gc()` "max used" (Ncells + Vcells, MiB) at stage end minus the baseline's R heap in use: an increase, not an absolute. It does not see Arrow, DuckDB, HDF5 or Python allocations. |
| Process peak | `VmHWM` of the stage process at stage end (MiB): the peak RSS since the baseline, absolute (it includes the baseline RSS and the embedded Python, which runs in the same process). The cgroup's `memory.peak` is recorded as a cross-check only (it also counts page cache, so file writes raise it). |
| Python peak | `tracemalloc` peak inside the embedded interpreter, started at the trainer's entry (immediately before MOFA2's `create_mofa`, `prepare_mofa` and the hand-off) and read when training returns; numpy's buffers are traced. tracemalloc is on in every training stage of every build, so its overhead is the same on both sides of every comparison; absolute training times include it. |
| Wall time | Elapsed seconds from the baseline to the stage's end (package loading excluded; the load time is recorded apart). |
| Marks | Inside the training stage: time and VmHWM at the trainer's entry (alignment done), at training's end, and after the product is written; inside the alignment stage, at the trainer's entry. No `gc()` is called at a mark. |
| Largest selected data set | Rows x samples x 8 bytes of the largest aligned data set handed to the trainer (from the trainer's arguments). |
| Dense floor | Sum over data sets of selected rows x samples x 8 bytes: what MOFA2/mofapy2 must hold as dense doubles. |
| Noise floor | Of a figure, at a size, for a build: the range (max - min) of its 3 repeats. Of a comparison: the larger of base's and C's noise floors for that stage, figure and size; at a size without repeats, the nearest smaller size's (labelled). |
| Tolerance | Memory: max(64 MiB, 5 % of the before value, 2 x noise floor). Time: max(1 s, 5 % of the before value, 2 x noise floor). |
| Value compared | The median of the completed, valid repeats (one value for A, B1, B3). |

## 5. Stages

Every stage is the build's own code path, called as the module calls it; the harness adds only measurement.

| Stage | What runs (in the build's namespace) |
|---|---|
| input write | The study generated (before the baseline), then `writeIntegrationInput()` for each omic (base: TSV, contract v1; A to C: Parquet, contract v2). |
| import | `newIntegrationStateOwner(<new folder>)`, `importIntegrationInputs(manager, <the omics' folders>)`. |
| alignment | The store reopened; `runIntegrationModel(manager, <capped parameters>, trainer = <stand-in>)`: the build's own model path (kernel or partitioned variant, its reads, its alignment and hand-off naming) up to the trainer's call; the stand-in records the data sets' sizes and stops there. |
| training | The store reopened; `runIntegrationModel(manager, <parameters>, trainer = <measuring trainer>)`, the measuring trainer calling the build's own `integrationTrainMofa(views, options)` unchanged. The stage therefore includes alignment, hand-off and MOFA2 training, and product write and acceptance; the marks split it. |
| results read | The store reopened; the accepted product hydrated and `integrationResults(value, integrationImportedObject(manager), list(), integrationGuiTexts()[, manager = manager])`, as `R/mod_integration_results.R` of each build calls it. |
| floor (uncapped) | As training, without feature selection and with 10 rounds (base and C only). |

The trainer argument also skips `integrationMofaCapability()`'s probe in every build (it checks only that MOFA2 and
mofapy2 are importable, which the baseline already did).

## 6. Model settings (passed explicitly in every build)

`technical_replicates = "average"`, `internal_standards = "drop"`, `missing_samples = "keep_missing"`,
`feature_selection = "top_variance"`, `features_per_view = 1000`, `scaling = "none"`, `scale_views = TRUE`,
`center_groups = TRUE`, `num_factors = 5`, `convergence_mode = "slow"`, `maxiter = 10000`, `seed = 20261119`,
`precision = "double"`; build C also `variant = "partitioned"`. Floor rows: `feature_selection = "none"`,
`features_per_view = NA`, `maxiter = 10`, the rest the same.

## 7. Sizes and the synthetic studies

| Size | Samples x features x data sets | Values (doubles) | Role |
|---|---|---|---|
| S | 20 x 5,000 x 4 | 0.4 M (3 MB) | lifecycle and noise floor; no claim |
| M | 200 x 30,000 x 4 | 24 M (192 MB) | claim size |
| L | 500 x 60,000 x 4 | 120 M (960 MB) | claim size if the feasibility probe passes |
| L' (`Lp`) | 300 x 40,000 x 4 | 48 M (384 MB) | replaces L if base does not complete L in the probe |
| P1 / P2 / P3 | 200 x 30,000 / 200 x 60,000 / 100 x 30,000, one measured data set plus a 1,000-feature companion (the import needs two omics) | | import pair (claim 3), build C |

Generator (`k150_study`, seeded per size: S 15001, M 15002, L 15003, L' 15004, P1 15011, P2 15012, P3 15013; the same
data in every build and repeat): four omics (proteomics, transcriptomics, lipidomics, metabolomics; one data set each,
named as INT-4 names them); one run per biological sample; two groups (A, B) alternating; five latent factors (factor
1 shifted by group, present in every omic; the others in two or one omics); per data set, loadings on a random 20 %
of features per active factor, feature means N(20, 2), feature noise SD U(0.3, 1.5), 3 % of values missing at random.

## 8. Comparisons

Per mechanism, each build against the row above it in section 3, on the stages that mechanism changes; and C against
base on every stage. For each of R heap peak, process peak and wall time: reduced (change below minus the tolerance),
no measurable change, or increased (above the tolerance). Reported, not pass/fail. Noise floors come from base and C,
3 repeats each, at S and M (and at L or L' only if the job's time allows, section 11).

## 9. Success statement

| # | Claim | Where | Pass when |
|---|---|---|---|
| 1 | Training R heap is bounded by the selected data | C, capped, at M and at L (or L') | In every repeat, the training stage's R heap peak <= 2 x largest selected data set + 256 MiB. The training stage includes alignment and the product write, so the test is conservative. |
| 2 | No regression in total memory | C against base, capped, at every size base completes | Pipeline process peak (per repeat, the largest process peak among the five stages; median over repeats) of C <= base's + tolerance (noise floor: the larger range of the two sets of repeats; without repeats, the training stage's noise floor of the nearest smaller size). |
| 3 | Import memory is bounded by one row group, not by file size | C, P1 / P2 / P3, 3 repeats each | (a) every import's R heap peak <= 4 x one row group (10,000 x runs x 8 bytes) + the feature and sample tables (their R object sizes as `readIntegrationInput` returns them, measured after the stage) + 64 MiB; (b) doubling features at fixed runs (P2 against P1) raises the median by no more than the tolerance; (c) doubling runs at fixed features (P1 against P3) raises it by no more than in proportion (P3's value) plus the tolerance. Methods wording if it passes: "import memory is bounded by one row group plus the feature and sample tables; it grows with the number of runs, not the number of features". |
| 4 | Python | every size, capped and floor, base and C | Python peak and process peak reported against the dense floor; not claimed reducible. |

MOFA2's own in-memory floor (the floor rows, and the capped Python peak) is stated separately as not reducible by the
package.

## 10. Terminal outcomes

Per row (claim, comparison or stage): pass; fail; base did not complete (stop reason recorded; no comparison made;
never used as a "before" or as a time denominator); after did not complete; instrumentation invalid (a figure not
observed through completion, the VmHWM reset not observed, or an induced failure's stop reason not recorded). A stage
stopped by the memory limit (cgroup `oom_kill` > 0), the time limit (`timeout`, exit 124) or the job's time limit is
"did not complete" with that reason; an R error is "did not complete (error)" with its message. An instrumentation
failure is kept apart from an unfavourable valid result; both are reported.

## 11. App-path rows (ENG-05)

At M and at L (or L'), base and C, in fresh processes, without a browser, by the functions the tabs call:
the Import action (`integrationCommitImport(manager, request)` on a new store, the request as
`integrationImportRequest` gives it for every data set; the tab's folder scan is not part of the action and is not
measured), then the Sample Matching preview (store reopened; base:
`integrationAlignmentPreview(integrationImportedObject(manager), views)`; C:
`integrationAlignmentPreview(integrationImportedObject(manager, INTEGRATION_PREVIEW_FIELDS), views, manager = manager)`,
as each build's `R/mod_integration_alignment.R`). Same definitions; reported without a claim. (Build C already
projects both calls, KIT-149; the Results tab still hydrates the inputs whole in C, which the results-read stage
measures.)

## 12. Runs, ceilings and attempts

| Step | Content | Job time limit |
|---|---|---|
| Lifecycle | One job: size S, all five builds (all stages; floor and app-path rows for base and C), then three induced failures on build C: a stage that allocates until a 4 GiB cgroup limit stops it, a stage that sleeps past a 240 s time limit, a stage whose import errors. Passes when every stage completes with valid figures and each induced failure is recorded after its baseline with its stop reason (`lifecycle.txt`). | 60 min |
| Feasibility probe | One job: base at L, capped pipeline and floor row. If any capped stage or the floor is killed or stopped, L is reported "base did not complete" and L' becomes the largest before/after size. | 300 min |
| Counted run | One workflow run, two jobs: `counted-small` (S and M: five builds once, base and C 3 times interleaved; floor rows for base and C; app-path rows at M; P1/P2/P3 on C 3 times) and `counted-large` (L, or L' per the probe: five builds once, floor and app-path rows for base and C; base and C repeated twice more only if the remaining time holds twice their first pipelines' time plus 10 minutes). If the probe's times project `counted-small` or `counted-large` beyond its limit, the repeats move to a third job of the same run (noise floors then include runner-to-runner variation, stated). | 240 min / 330 min |
| Correction pass | At most one, of the instrumentation, rerunning only the invalid cases, inside the remaining budget. | as needed |

Ceilings: each stage has its own time limit (write 30, import 60, alignment 60, training 150, results 30, floor 90,
app-path 60 minutes), cut to the job's deadline (job limit minus 12 minutes); a stage does not start with less than 3
minutes left. Disk: before each stage the measurement folder is checked; the job stops and reports when it passes 10 GB;
each run's folder (inputs and stores) is removed after its stages by a guarded remover. Budget: at most 66 core-hours
(990 runner-minutes on 4 cores) in total, counted from the jobs' actual durations; the job limits above sum to 930
minutes, leaving 60 for a correction pass.

## 13. Not measured or not claimed

- The app as a whole (Shiny session, worker process, the folder scan); only the two app-path rows above.
- The KC datasets or any real data.
- The final default model: the cap is fixed at 1,000 and the selection statistic is the current one (variance).
- Training time of the floor rows (10 fixed rounds).
