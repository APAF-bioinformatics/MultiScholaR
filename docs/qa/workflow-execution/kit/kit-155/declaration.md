# KIT-155 declaration: re-measurement of the alignment and training stages after pages sized by bytes

Fixed before the lifecycle run and the counted run. Lane: branch `run-mofa-measure-k155` (from `run-mofa-measure-k150`
at a1a19dfd), workflow `.github/workflows/k150-measure.yml` (trigger moved to this branch), scripts
`.github/mofa-smoke/k150-*` (two additive changes, below) and `k155-builds.tsv`, `k155-summarise.py`. Everything not
stated here is as in `docs/qa/workflow-execution/kit/kit-150/declaration.md` sections 3 to 10 (terms, figures, cgroup,
baseline, VmHWM reset, sizes and generator, tolerances, terminal outcomes).

## 1. Purpose

KIT-150 (INT-27) found claim 1 failing: the training stage's R heap grows with the number of samples because each
alignment page held 10,000 features across every run. KIT-155 sizes the pages by a byte budget (8 MiB of values: rows =
floor(8 MiB / (8 bytes x columns read)), 1 to 10,000). This run measures, at KIT-150's sizes, whether claim 1 now holds
and what the change costs in time. Only the alignment and training stages are measured and compared; write and import
run as set-up (their rows are recorded, not compared).

## 2. Bindings

| Build | Commit | Content |
|---|---|---|
| before | 12d68ca78442a50b46f52e33525dd85cf86df5e7 | performance-v2 (pages of 10,000 rows) |
| after | f263bc6db75a63dd7abf2281b2f711cf98b6c457 | kit-155-alignment-pages (pages of 8 MiB) |

Both with `variant = "partitioned"` passed explicitly. Environment, library, MOFA2 closure, bundle Python, launcher,
one library per build, the build-commit marker: as KIT-150 section 3 (MultiScholaR of the snapshot library at build
after). Model settings: KIT-150 section 6, plus `training_seeds = 1` (KIT-151 added training seeds, default 5, after
KIT-150's builds; one seed keeps the stage the same single training KIT-150 measured). The other parameters added
since KIT-150 (`min_observed_fraction` 0.5, `input_matrix` uncorrected) take their defaults, the same in both builds.

## 3. Stages and sizes

Per build and repeat: write, import, alignment (`runIntegrationModel` up to the trainer, stand-in stops there),
training (`runIntegrationModel` with the measuring trainer), each in a fresh process in its own cgroup (14.5 GiB, no
swap), time limits as KIT-150 section 12. Sizes S (20 x 5,000 x 4), M (200 x 30,000 x 4), L (500 x 60,000 x 4), same
seeds and generator. Order: for each size S, M, L, for repeat 1 to 3, before then after (interleaved).

## 4. Claim 1 and how it is read

| Form | Pass when (every repeat) |
|---|---|
| As declared (KIT-150 section 9) | training R heap peak <= 2 x largest selected data set + 256 MiB |
| Restated (INT-27, KIT-150 section 15) | training R heap peak - overhead <= the declared bound + 2 x noise; overhead = the same build's training R heap peak at S (median of 3); noise = the larger range of the 3 repeats at S and at the size |

Evaluated at M and L for both builds (before is the comparison; after is the claim). Expected from KIT-150's figures
(not a result): the declared form fails in both builds (the overhead alone is about 262 MiB against bounds of 259 and
264 MiB); the restated form fails for before at L. The 8 MiB budget was chosen from KIT-150's figures (training R heap
grew by about 13 MiB per MiB of page between M and L), not tuned on this run, and is not changed after it.

Also reported (no claim): alignment stage figures; process peaks; Python peak against the dense floor; per stage the
comparison before against after (R heap peak, process peak, wall time; median; tolerance max(64 MiB or 1 s, 5 %, 2 x
noise floor); reduced / no measurable change / increased). Time cost = the alignment and training wall times, after
against before. Model sanity: the factors (rounded to 1e-8, SHA-256) and iterations the same in both builds at every
size; a difference is a defect of KIT-155, reported.

## 5. Runs, ceilings, attempts

| Step | Content | Job limit |
|---|---|---|
| Lifecycle (`k155-lifecycle`) | Size S, both builds' four stages once; then on build after a stage stopped by a 4 GiB memory limit, one stopped by a 240 s time limit, one whose import errors. Passes when every stage completes with valid figures and each induced failure is recorded after its baseline with its stop reason (`lifecycle.txt`). | 60 min |
| Counted (`k155-counted`) | One job: S, M, L, both builds, 3 repeats each. Estimated from KIT-150's stage times about 160 minutes. | 330 min |
| Correction pass | At most one, of the instrumentation only, rerunning only invalid cases. No thresholds, sizes, budget or builds change. | 120 min |

Budget: at most 510 runner-minutes (34 core-hours) in all. A stage or the job stopped by a limit is "did not complete"
with its reason; instrumentation invalid is kept apart from an unfavourable valid result.

## 6. What each outcome means

- After passes the restated claim at M and L: pages sized by bytes bound the training R heap above its fixed overhead
  at KIT-150's sizes; recorded in INT-29. The declared form is reported as it falls.
- After fails the restated claim: a valid unfavourable result, reported with its figures; no other budget is tried
  under this ticket (a smaller budget is a new decision).
- Before does not complete at a size: no comparison at that size.
- Factors or iterations differ between builds: KIT-155 changed results; reported as a defect, the claim is not stated.
- Time: the change in alignment and training wall time is reported as measured, whatever its sign.

## 7. Record (2026-11-22; terms, thresholds, sizes and budget unchanged)

- Lifecycle run 44952799013: PASS (both builds' four stages at S valid; memory-limit, time-limit and error failures
  recorded after the baseline). 19.6 runner-minutes.
- Counted run 44954546145 (runner Intel Xeon 6973P-C): 72 stage rows, 0 instrumentation invalid; S and M complete
  (3 repeats each build); at L, before 2 and after 1 repeats. The other L cases were "not run: the disk ceiling": the
  guarded remover (`k150-remove.py`) refused the `before`/`after` run folders, so no run folder was removed. An
  instrumentation failure of this lane copy, not a result. 88.3 runner-minutes.
- The one correction pass, run 44962192602 (runner Intel Xeon Platinum 8573C): the remover's name rule extended to
  the two build names; only the three L cases not run were run (after r2, before r3, after r3), all completed and valid.
  65.2 runner-minutes. The L repeats therefore span two runners: R heap figures agree to 1.2 MiB; wall times and the
  models' iterations differ by runner (below). The combined tables (`summary.md`, `rows.tsv`, `claims.tsv`,
  `comparisons.tsv`) are `k155-summarise.py final` over both jobs' artifacts; `counted/` and `correction/` keep each
  job's own. (One duplicate row: L after r2 write completed in both jobs; it enters only the write stage, not compared.)
- Build after was f263bc6d. Its later commits on the branch (d1fd2597: a page holding only internal standards gives an
  empty block, reached only at tiny pages with standards; 5ed06132: a test's expected type) do not change any path these
  synthetic studies take (they have no internal standards); not re-measured.
- Total 173.1 runner-minutes (about 11.5 core-hours of 34).

Results (R heap peaks have a range of 0.0 to 1.2 MiB over repeats):

| Claim 1 | Build | M | L |
|---|---|---|---|
| as declared (<= 2 x largest selected + 256 MiB) | before | fail (541.8 vs 259.0) | fail (848.2 vs 263.6) |
| as declared | after | fail (463.9 vs 259.0) | fail (558.8 to 560.0 vs 263.6) |
| restated (growth above the S overhead) | before | fail (280.0 vs 259.0) | fail (586.4 vs 263.6) |
| restated | after | pass (208.6 vs 259.0) | fail (303.5 to 304.7 vs 263.6 + 2.4) |

Overheads (training R heap at S): before 261.8, after 255.3 MiB. Before against after (median): alignment R heap -77.9
MiB at M, -289.3 at L; training R heap -77.9 at M, -289.3 at L; process peak, alignment -107.6 at M and -383.8 at L,
training no measurable change at M and -293.1 at L. Time, compared on one runner (before r1/r2 and after r1 of the
counted run; before r3 and after r2/r3 of the correction run): at L the alignment stage +16.0 s (58.5 to 74.5) and
+15.4 s (66.8 to 82.2), about +27 %; the training stage (which aligns first) +32.7 s (649.7 to 682.4) and +32.0 s
(806.8 to 838.8), about +5 %; at M and S no measurable change. Models: on each runner the factors (rounded to 1e-8)
and iterations are the same in both builds at S, M and L; at L the two runners' processors give different models
(4,786 and 4,361 iterations) in both builds alike (INT-25: unpinned maths library).

Reading (section 6): the restated claim passes at M and fails at L by about 40 MiB, a valid unfavourable result. Pages
of 8 MiB lowered the training R heap at L by 289 MiB (before failed the restated bound by 323 MiB, after by 41 MiB); what is left grows with the samples by about 96 MiB between
M and L (training R heap 463.9 to 558.9 MiB with pages of the same 8 MiB), so it is not the page alone. No other budget
is tried under this ticket.
