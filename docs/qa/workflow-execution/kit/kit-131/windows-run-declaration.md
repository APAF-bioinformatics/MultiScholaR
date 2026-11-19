# Windows verification run for KIT-131 and KIT-052 (declared 2026-11-19, before launch)

## Purpose

One run on a GitHub-hosted `windows-2025` runner, with the installed application, for the Windows criteria of two
tickets that are still open:

- **KIT-131** (the second launcher): step 7 on Windows, a GUI-entry lane with every managed step on a plain worker
  (launcher `"process"`), compared with the same lane in process on the same runner (design A1) and with its oracle.
- **KIT-052** (Windows paths): the Summary repair `767ad7d6` (the session archive staged under `tempdir()`; a failed
  copy into `files/<sha256>` named; the save's failure shown as an error), so far verified on Linux only.
- **Kept:** the four browser e2e lanes without DIA still pass with the laptop's counts (88, 190, 16, 153).

## Bindings

| | |
|---|---|
| Branch | `run-win-k131` (run-only), cut from `origin/performance-v2` at `791cc5f3` (contains kit-134, kit-052-win-paths to `0aecabd1` with `767ad7d6`, and KIT-131's launcher) |
| Commit run | the commit after this declaration's. It changes only `packaging/verify/` and this branch's copy of `.github/workflows/distributable.yml`. Nothing under `R/`, `tests/`, `DESCRIPTION` or `NAMESPACE`; no oracle |
| Workflow | `KIT-134 distributable`, dispatched on `run-win-k131`, stage `full`, platform `windows` (the Linux jobs are not needed as a reference and are skipped) |
| Runner | `windows-2025` (image as GitHub serves it; the image version is read from the log) |
| R | 4.6.1 (the workflow's) |
| Rust | pinned at 1.99.1 on the GNU target, as kit-134 `075a2640` set it (1.100.0 does not build the internal packages on Windows, run 44368107499) |
| Read token | handled exactly as the existing workflow does: given to the one install step only; no cache; no artifact (the build job uploads nothing) |

### Instrumentation (packaging/verify and the workflow copy only)

- `packaging/verify/lane-worker.R` (new): the KIT-131 step 7 driver of the Linux run
  (`kit-131/step8/lane_worker.R`), adapted to the installed application's R. The module's session gets a plain managed
  application, built as `workflowManagedMain` builds it for the plain launcher; after each input the driver waits until
  the omic's runner is idle. It prints the expectations, the comparison with the oracle and with the in-process record
  of the same job (step by step and field by field, the managed step's own "Applying ... step ..." notice set aside),
  the number of managed waits and of distinct worker processes seen as children of the controller, and whether the
  application closed with no worker left.
- `packaging/verify/figure-trace.R`: also traces `omicSummaryNotify` (every Summary notification as shown to the user,
  `SUMMARY` lines) and `rmarkdown::render` (the file each render wrote, `RENDER` lines). The existing `RAISED` lines
  already catch the classed `workflow_figure_unwritten` refusal that `767ad7d6` raises for a failed `files/<sha256>`
  copy.
- `packaging/verify/figure-paths.R`: after the `prot_tmt_lane` item, the lane record's Summary steps (`REPORT` lines:
  the `.docx` among the artifacts, the notifications), when the lane wrote a record.
- `packaging/verify/windows-figure-cases.ps1`: prints the new line types.
- The workflow's Windows job: the GUI-entry lanes step runs `lipid_canonical` only (the in-process reference for
  step 7); a new step runs `lane-worker.R` on `lipid_canonical`; the figure steps run the Summary cases below instead of
  the earlier figure cases (those were confirmed in 44379470655 and are not repeated); one more row in the results table.

## Cases

The project root's length is set through `TMPDIR` as before (prot_tmt root = TMPDIR + 59 to 60 in 44275615462 and
44379470655). Roots below are targets; the run states the root reached.

| # | Ticket | Long paths | Root (TMPDIR) | Item | Counts as pass |
|---|---|---|---|---|---|
| W1 | KIT-131 | runner default (on) | runner's temp | `lipid_canonical` on plain workers | the lane runs to its end with at least one managed wait and a worker process seen; every step identical to the in-process lane of the same job after the working notice is set aside, except store-internal entity counts (the difference accepted on Linux at step 8); no worker left after the application closes |
| W2 | KIT-131 | on | runner's temp | `lipid_canonical` in process (lanes step) | the reference for W1; its own differences from the oracle are recorded, not judged here |
| S1 | KIT-052 | off | about 91 (31) | prot_tmt entry lane | the Summary `.docx` is rendered (`RENDER` with a `.docx`; "Report generated successfully!"); nothing refused, nothing lost silently |
| S2 | KIT-052 | off | about 96 (36) | prot_tmt entry lane | as S1 (the highest root below the Q-value staged copy's calculated 99) |
| S3 | KIT-052 | off | about 120 (60) | prot_tmt entry lane | see the note: expected to stop loudly at QC (as at 117 in 44379470655); pass = nothing lost silently. If the lane reaches the Summary, pass = the `.docx` rendered |
| S4 | KIT-052 | off | about 130 (70) | prot_tmt entry lane | if the save is reached: the save shows the error naming the `files/<sha256>` copy (`SUMMARY error` and `RAISED` naming the file), and copy and report stay refused (no `RENDER` of a `.docx`). If the lane stops at QC first: not reachable, recorded as such |
| S5 | KIT-052 | on | about 130 (70) | prot_tmt entry lane | nothing else changes with long paths on: every figure present (as at 126 to 131 in 44379470655) and the Summary `.docx` rendered |
| S6 | KIT-052 | off | 63 (20) | `test-workflow-figure-paths.R` | its Summary tests (767ad7d6's, with copy doubles) pass on Windows; the fixture errors of 44379470655 (folders over 248) are expected and recorded |
| E | kept | on | runner's temp | the four browser e2e lanes | 88, 190, 16, 153 |

**Note on S3 and S4 (written before the run).** The ticket asks for the Summary at about 120 with long paths off and
for the refusal above 125. Run 44379470655 showed the prot_tmt lane stopping loudly at QC at root 117 with long paths
off (the Q-value staged copy, calculated limit 99), before normalisation, DA and the Summary. If that holds, the Summary
cannot be reached above about 99 with long paths off on this lane, and with long paths on a copy over 259 does not fail,
so the `files/<sha256>` refusal cannot be shown on the real lane without a package change (shortening the Q-value staged
copy's path). S3 and S4 are run anyway, to show what happens; S6 is the only direct evidence of the refusal's text, with
a copy double.

### KIT-131 criteria not exercised by this run

- **Screen journeys under the worker:** the kit's screen journeys are Linux journeys under systemd
  (`tools/profiling`, a journey record, DIA routes); under the plain launcher a journey record is refused by design
  (L2.8) and no DIA route runs on Windows (D3, KIT-134 D4). Not run here; they need the laptop.
- **Audit F10 (the stop error is swallowed):** a change to `R/utils_workflow_execution.R` (record the message), package
  code; not in this run's scope. Reported, not exercised.

## Runtime ceiling and attempts

- The Windows job's own ceiling (300 minutes); the worker step 60 minutes; the figure steps 120 minutes each as before.
- **One run plus one bounded correction pass**, the pass used only for an instrumentation fault (a script error, a
  build fault outside the package, a case that never started). An unfavourable valid result is a finding and is not
  rerun.

## Decision after the run

Each KIT-131 criterion and each KIT-052 case is stated met or not met with its evidence in the results section below,
in `kit-134/outcome.md` and as dated Status lines in the tickets. Nothing is ticked; the user closes tickets.

## Results

**Run [44457268709](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44457268709)**, dispatched
2026-11-19 on `2c3a9b21` (declaration `eb72cdde`), Windows job only, 04:17 to 05:09 UTC. Image `windows-2025-vs2026`
20261103.298.1, rustc 1.99.1. Both internal packages built; library, installer, install and lifecycle passed. The
run's conclusion is "failure" only because the in-process lanes step fails on `lipid_canonical`'s known differences
from its oracle (27 passed, 8 failed: `norm_log` and `session_summary` text, as in the Windows runs before). 0
artifacts; no cache. **The correction pass is not used**: no step had an instrumentation fault; the results below are
valid, some unfavourable.

### KIT-131

| Criterion | Result | Evidence |
|---|---|---|
| The installed application chooses the second launcher on Windows | **met** | lifecycle test: probe and launcher `process`, "win32: no systemd"; network cut |
| W1: a lane with every managed step on a plain worker runs to its end | **met** | `lipid_canonical`, 20 steps in order; 9 managed waits (130 s on workers); 38 distinct worker processes seen as children of the controller; the application closed with 0 workers left; lane 177 s |
| W1: identical to the in-process lane of the same runner except store-internal counts | **not met** | identical at 6 of 20 steps (5 raw, 6 with the 19 "Applying lipidomics step ..." notices set aside). Beyond the store entry counts (`qc_intensity` 30 against 31, as on Linux; 101 against 144 at the end), from `norm_refusal` on: the lipid QC figures under `results/lipidomics/lipid_qc/` (`lcms_{pos,neg}_{pre,post}_norm_{correlation,density,pca,rle}.png`) are **absent** on workers and present in process, and the normalisation tab's plot outputs are a `<table>` digest on workers where in process they are `<plot> alt=...`. **No notification differs**: the figures are lost without a message (this instrumentation printed at most 6 keys per field). Against the oracle the same 14 steps differ the same way. On Linux (step 8) the same lane on plain workers was identical apart from the store count, so this is Windows-specific. Cause not read here (the worker logs and the full record stayed in the job) |
| Screen journeys under the worker | **not exercised** | Linux systemd journeys (journey record, DIA routes); refused by design under the plain launcher (L2.8); they need the laptop |
| Audit F10 (the stop error is swallowed) | **not met, not exercisable** | `R/utils_workflow_execution.R`, `workflowProcessStop`: `tryCatch(..., error = function(error) NULL)` around `ps_kill_tree` and `kill`; a package change |

### KIT-052 (prot_tmt entry lane in process; roots as reached)

| Case | Long paths (as R read it) | Root | Reached | Summary | Figures |
|---|---|---|---|---|---|
| S1 | off (0) | 90 | the end | **`.docx` written** (`TMT_report_proteomics_prot_tmt.docx`, "Report generated successfully!"); save "Study parameters saved successfully"; copy "Publication files copied" | 121 of 121 final ok (longest final 250) |
| S2 | off | 96 | the end | **`.docx` written**, same notifications | 121 of 121 (256) |
| S3 | off | 120 | QC | **not reached**: the lane stops at QC; 4 classed refusals naming the 48 Q-value staged copies (`workflow_figure_unwritten/saved_entity_outputs_pending`), then the lane's own read fails ("cannot open the connection") | 61 drawn, 13 final ok; every missing final named in a refusal |
| S4 | off | 130 | the import | **not reached**: the lane errors after 33 s ("cannot open the connection"), 9 figures drawn, **no classed refusal raised**; what the user would see is not captured | 9 of 9 |
| S5 | on (1) | 130 | the end | **`.docx` written**, same notifications as S1 | 121 of 121 (290) |
| S6 | off | 48 (fixtures) | `test-workflow-figure-paths.R` | 92 passed, **0 failed**, 7 errors, 4 skipped (44379470655: 84 passed, 7 errors). The 767ad7d6 refusal is raised on Windows with its text: "a file of the session archive could not be copied into its staging folder ... the longest of these paths is 323 characters: the project folder's path is too long for Windows with long paths disabled ..." (class `workflow_figure_unwritten/saved_entity_write`) | the errors printed are fixtures making folders over 248 characters, as before; 3 of the 7 are not printed |

- **767ad7d6 confirmed on Windows where it can be reached:** at 90 and 96 with long paths off the Summary report is
  written (not at 91 in 44379470655), and with long paths on at 130 nothing else changes.
- **The ">125 with long paths off" refusal cannot be shown on the real lane:** with long paths off the lane stops at
  QC from about 99 (S3), or earlier at 130 (S4); with long paths on a copy over 259 succeeds. Only S6 (a copy that
  really fails at 323 characters in the test) shows the refusal's text on Windows. Reaching it in the lane needs a
  package change (the Q-value staged copy's path), or S6 is accepted as the evidence.
- S1's 2 extra failures (33/4): the session archive's name `native-export-<hex>` has 9 hex digits in the oracle, which
  the recorder's mask (10 or more) leaves, and 10 or more here: the tempfile name's length, not a product difference.
- **S4 is new:** with long paths off at root 130 the lane stops early with an unclassed file error, before any figure
  check. Package code; not this run's to repair.

### Kept

Browser e2e lanes from the installed application: lipidomics-canonical 88, metabolomics-lc-gc 190,
transcriptomics-screen 16, proteomics-tmt-lfq 153: the laptop's counts, 0 failures.

### What remains

1. KIT-131: the Windows worker route loses the lipid normalisation QC figures silently. Next bounded item: one
   diagnostic run printing the worker logs and the full step difference for `norm_refusal` and `norm`, then a package
   repair (likely the same family as KIT-052's worker-route copies). Then step 7 on Windows again.
2. KIT-131: the screen journeys under systemd on the laptop; F10 (record the stop error), a package change.
3. KIT-052: the decision whether S6 stands for the ">125, long paths off" refusal; S4 (an unclassed early error at
   130 with long paths off) to be read.

## Diagnostic run (declared 2026-11-19, before launch)

**Purpose.** Read the cause of the silent loss of the lipid normalisation QC figures on Windows workers (W1 above) and
of S4's unclassed early stop, before any package change. Instrumentation only; nothing under `R/`, `tests/`,
`DESCRIPTION`, `NAMESPACE`, oracles or the harness-bound folder names.

**Hypothesis to falsify (from the code, not yet observed).** The plain launcher (`workflowLaunchProcess`) gives the
worker `TMPDIR = <request folder>/worker-scratch`, deep in the project. R on Windows takes its temporary folder from
`TMPDIR`, so `withWorkflowShortFigureFolder`'s `tempfile("wfx-")` is not short in the worker: a QC figure is drawn at
a path over 259 characters, ragg writes nothing and says nothing, `writeWorkflowOmicsFigure` returns NULL for a file
that is not there, and `generateLipidQcPlots` logs and continues. On Linux nothing limits the path, so the same lane is
identical. Falsified if the probe's ggsave in the worker's `wfx-` folder writes the file, or the replayed figures are
drawn at paths of 259 characters or fewer and are still missing.

**Bindings.** Branch `run-win-k131` at the instrumentation commit after this one (on `32cbe69e`); workflow
`KIT-134 distributable`, stage `full`, platform `windows`; same runner image, R 4.6.1, Rust 1.99.1, token handling, no
cache, no artifact.

**Instrumentation.**
- `packaging/verify/lane-worker.R`: traces `workflowLaunchProcess` (LAUNCH lines: operation, the worker's `TMPDIR` and
  its length). Before the first worker, the first `lipid_qc_figures` and the first `lipid_composite_figure` worker, it
  starts `packaging/verify/worker-probe.R` through the same function with the same request (the worker's environment,
  `TMPDIR`, working directory) and prints its log: environment, `tempdir()`, `capabilities()`, `bitmapType`, ragg/Cairo/png,
  ggsave's device, ggsave/agg_png/png/writeLines of `lcms_pos_pre_norm_correlation` in the worker's `wfx-` folder and in
  the controller's temporary folder; for the `qc_figures` request the handler replayed with the figure store switched
  off (nothing saved) and every figure traced (drawn path, final path, lengths, existence), then the replay's publication
  folder removed so the real worker starts as before. After the lane: every worker's log (in full for the norm-tab
  publications), its result and its request publication folder; the `norm_refusal` and `norm` steps' difference with
  the in-process record in full.
- `packaging/verify/figure-paths.R`: for `prot_tmt_lane`, every failing `file()`/`gzfile()` call with its path, length,
  folder existence and calls (CONN lines), and the testthat trace of each error (TRACE lines).
- Workflow copy: in-process `lipid_canonical` (reference) and the worker lane as before; e2e lanes and the long-paths-on
  figure step not run; with long paths off only `off-tmt130=70:prot_tmt_lane`.

**Ceiling and attempts.** The job's own ceilings; this is the one diagnostic run. An instrumentation fault (script error,
probe that never ran) is recorded as such and not rerun within this item's budget.

**Decision.** The cause stated with its evidence, then the package repair on `kit-131-win-worker`.

### Diagnostic run result (run 44464483568, 2c3a9b21 + eb1f9bd4)

- **Cause read (valid):** the `lipid_qc_figures` worker's own log ends "Generated 0 QC plots for stage: post_filter"
  after drawing every plot (plotRle, notch messages), with no error or warning from the figure code; all 16 QC PNGs and
  the composite are absent on workers and present in process (`norm_refusal`, `norm` difference printed in full; the
  normalisation tab's 16 plot outputs are `<table>` digests because their files are missing). The worker's `TMPDIR` is
  `<project>/scripts/lipidomics/.dia-native-<14>/.workflow-entities/.workflow-export/wfx-<32>/worker-scratch`, so its
  `tempfile("wfx-")` figure folder gives about 209 characters after a project root of at least 66 on this runner:
  every QC file name is drawn at 266 to 275 characters, where ragg writes nothing (KIT-052, run 44241609555) and
  `writeWorkflowOmicsFigure` returned NULL for the absent file. On Linux nothing limits the path.
- **Instrumentation faults (recorded, not rerun as declared):** (1) the LAUNCH and PROBE lines were printed while the
  lane ran and the lane's harness captured them, so the probe's environment, device and replay lines are lost;
  (2) the long-paths-off step failed before its case started (`$OUT\figures` is created by the skipped long-paths-on
  step), so S4's trace (item 3) was not taken. Both are fixed for the confirmation run below.

## Confirmation run (declared 2026-11-19, before launch)

**Purpose.** KIT-131 step 7 on Windows again, with the repair `kit-131-win-worker` (d1fa3ada: in a worker the short
figure folder is made in `TMP`/`TEMP`, not under the request-owned `TMPDIR`; a figure drawn but not written stops the
step with `workflow_figure_unwritten`; F10), verified on Linux in run-tests 44468417325 (no new failure, check identical).

**Bindings.** Branch `run-win-k131` = this declaration's commit, then the instrumentation fixes, then a merge of
`kit-131-win-worker` at d1fa3ada; workflow `KIT-134 distributable`, stage `full`, platform `windows`; same image,
R 4.6.1, Rust 1.99.1, token handling, no cache, no artifact.

**Instrumentation changes.** `lane-worker.R` keeps the LAUNCH and PROBE lines and prints them after the lane;
`windows-figure-cases.ps1` creates its evidence folder; the e2e step runs again; the long-paths-off step runs
`prot_tmt_lane` at about 90, 120 and 130 (TMPDIR 31, 60, 70) with CONN and TRACE lines.

**Cases and what counts as a pass.**

| # | Item | Pass |
|---|---|---|
| C1 | `lipid_canonical` on plain workers (long paths as the runner has them) | runs to its end; every step identical to the in-process lane of the same runner after the working notice is set aside, except the store-internal entity counts (accepted at Linux step 8); any other difference explained; no worker left; the probe's replay draws the QC figures under `TMP` at fewer than 260 characters |
| C2 | the four browser e2e lanes | 88, 190, 16, 153 |
| C3 | prot_tmt lane, long paths off, roots about 90, 120, 130 | 90: Summary `.docx` written as in 44457268709; 120: stops loudly at QC (classed refusal), as before; 130: the CONN/TRACE lines name the file, its length and the call that raised "cannot open the connection" (item 3 read, not repaired in this run) |

**Ceiling and attempts.** The job's own ceilings; one run. The one bounded correction pass of this work item is kept
for a repair after this run if C1 fails for a package reason.

**Decision.** C1 met or not met per step; item 3's cause stated with its evidence; then the tickets' Status lines.

### Confirmation run result (run 44472789530 on 6b08bace = 3330b3ae + merge of kit-131-win-worker d1fa3ada)

The run's conclusion is "failure" only through the in-process lanes step's known oracle differences (27/8, as before).
No instrumentation fault.

- **C1 (met, with one difference explained):** `lipid_canonical` on plain workers ran to its end (9 managed waits,
  38 worker processes, 0 left after close). Against the in-process lane of the same runner: no figure, plot output,
  file or notification differs any more (the 16 QC PNGs, the composite and the 16 normalisation plot outputs are now
  present and identical). What remains: (a) the store's entity counts at every step (accepted at Linux step 8);
  (b) the normalisation log text (`norm-norm_log`, 8 steps from `session_export`) and `session_summary` (2 steps):
  the same keys also differ between the in-process lane and the oracle, and the in-process `norm_log` digest itself
  changed between two runs of the same lane on the same image (8c3a5ec9... in 44464483568, 97c7723c... here), so it is
  run-dependent text the recorder does not mask, not a worker difference.
  Probe evidence: the worker's `TMPDIR` is 244 characters and its `tempdir()` 255; before the repair a ggsave or
  `ragg::agg_png` in its `wfx-` folder (317 to 318 characters) wrote nothing and said nothing, while
  `grDevices::png` wrote; the worker has `TMP`/`TEMP` = `C:\Users\RUNNER~1\AppData\Local\Temp` (36). With the repair
  the replayed `lipid_qc_figures` handler drew all 8 figures under `TMP` at 91 to 99 characters and copied them into
  the request folder (259 to 267 characters, long paths on).
- **C2 (met):** e2e lanes 88, 190, 16, 153, 0 failures.
- **C3:** root 89: the Summary `.docx` written, 0 of 346 figures missing (35/2, the two known mask failures).
  Roots 118 and 129: the lane stops loudly at QC with the four classed refusals naming the 48 Q-value staged copies
  (`workflow_figure_unwritten/saved_entity_outputs_pending`). **Item 3:** the unclassed "cannot open the connection"
  is raised by the lane's own test step, `helper-gui-entry.R:794` (`summary_report`: `readLines` of
  `filtered_session_summary.txt`, which the refused workflow never wrote), not by package code; it is not a path-length
  failure of the package. The earlier 33-second stop at root 130 (44457268709) was not reproduced (root reached 129
  here, same TMPDIR); the CONN lines were again captured by the lane harness.
