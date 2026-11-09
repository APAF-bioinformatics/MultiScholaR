# KIT-108 W2: the memory instruments, their lifecycle proofs so far, and the lead's commands for M1 to M3

Written 2026-11-09 by the W2 agent. No member was started. Nothing under `R/`, `inst/`, `tools/` or `tests/` was changed.

## The instruments

All scripts are in `docs/qa/workflow-execution/kit/kit-108/memory/`. They are read-only towards the application.

| Id | File | What it reads | What it writes |
|---|---|---|---|
| I1 | `analyse-traces.R` | a member's `memory-trace.tsv`, `journey-phases.tsv` and `driver-phases.jsonl` | the per-window table, plus a per-member summary (`<out>-summary.tsv`) |
| I2 | `sampler.sh` | `systemctl --user list-units` and `show -p ControlGroup`; each unit's `memory.stat` and `cgroup.procs`; per PID `/proc/<pid>/{comm,status,smaps_rollup}` | one TSV row per PID per tick, one per unit cgroup per tick, and event rows. Every row carries epoch and uptime. |
| I2 (reader) | `analyse-sampler.R` | I2's TSV and the member's harness files | `<prefix>-check.tsv` (the M1 criteria), `-tracking.tsv`, `-windows.tsv`, `-children.tsv` |
| I3 | `replica.R` | a copy of a member's project tree, through the build's library (or `--load-all <tree>`) | `replica-rows.tsv`, `rprof-{copy,render}-by-{total,self,line}.tsv`, `object-sizes.tsv` |

**I2's design points:**
- **Discovery:** units are found through `systemctl --user` only (no `pgrep -f`). Every PID comes from the unit's cgroup (nested cgroups included); a `.slice` unit contributes its cgroup row only.
- **No ptrace:** `smaps_rollup` and `status` of the user's own processes are readable without ptrace; the probe read them.
- **Clocks:** each row carries both clocks (`date +%s.%N` and `/proc/uptime`), so the driver's epoch marks map onto the harness's uptime trace exactly. I1 can only bracket that mapping (M0: ±0.01 to ±4.09 s).
- **Ending.** It never fails on the unit's state:
  - `--pattern` waits up to `--wait` for a match, then samples until no match has been active for `--grace` seconds (default 30). That gives `unit_gone`, exit 0. A shorter gap is recorded as `unit_gap` and sampling continues.
  - `--unit` gives `unit_gone`, exit 0, as soon as the unit is not active.
  - `never_seen` and `max_reached` also exit 0.

**I3's design points:**
- **Setup:** `project_dirs` comes from the application's own "reuse existing" call, `setupDirectories(base_dir, omic, label = NULL, reuse_existing = TRUE)` (`app_server.R:363-369`). It is assigned to `.GlobalEnv` only when absent, as the handler does. `time_dir`/`timestamp` are re-pointed at the member's own `filtering_qc/<timestamp>` directory, because `setupDirectories` mints a new one and the copy reads `time_dir`.
- **Inputs:**
  - The contrasts and design are read from `source_dir` (`contrasts_tbl.tab`, `design_matrix.tab`), which is the handler's fallback.
  - On the lipid screen project `contrasts_tbl.tab` is absent (the tree has `contrast_strings.tab`), so `contrasts_tbl = NULL` reaches `copyToResultsSummary`, which loads `contrast_strings.tab` itself. In the application the bus supplies it.
  - The label is the driver request's `journey$experiment_label`.
- **The calls:** `copyToResultsSummary(omic, label, contrasts, design, force = TRUE)`; then the handler's template step (`scripts/<omic>/` first, else the package's `reports/<dir>/`); then `RenderReport(omic, label, rmd)`.
- **Per call:**
  - `gc(reset = TRUE)` and a row before;
  - after it, a row with no gc (`after_<call>_raw`, the state the application is left in), then `gc()` and a row (`after_<call>`: used, max used since the reset, `/proc/self/status`, `mallinfo2`);
  - an `Rprof(memory.profiling = TRUE, interval = 0.05, line.profiling = TRUE)` profile summarised by `summaryRprof(memory = "both")`.
- **At the end:** `lobstr::obj_size` of every `.GlobalEnv` object (the replica's own are tagged `replica-own`) and of both return values; then `final_gc`; then `malloc_trim(0)` with `/proc/self/status` after it (`after_malloc_trim`).
- **The probe:** `mallinfo2`/`malloc_trim` go through `tools/profiling/workflow_checkpoint_gui_diagnostic.R:41-60`'s Rcpp code, compiled in the replica under `TMPDIR`.
- **Rows and failure:** each row is appended as it is taken. `main()`'s `on.exit` appends a final `failed` row with the phase and the error. The exit status is 0 on completion and 1 on failure.
- **The work copy** (`<out>/work`) is made with `cp -a`; on `/home` this is a cheap reflink copy (135 MB in 0.03 s). It is removed after a completed run and kept after a failure; remove it by hand.

## Lifecycle proofs run here (outputs in `lifecycle/`)

### I2 (`lifecycle/i2/`, log `lifecycle.txt`)

| Case | What | Recorded | Exit |
|---|---|---|---|
| completion | `--pattern 'kit108-sampler-probe*'`, a scratch unit with two processes (`R` holding about 230 MB anon, and `sleep`), stopped after 8 s | `start`, `unit_seen`, 8 cgroup rows, 16 process rows (R and sleep each tick), then `unit_gone` "after sampling"; all rows have 19 fields; max tick gap 1.045 s; last tick R `RssAnon` 237,420 kB against the cgroup's anon 237,536 kB (ratio 0.9995) | 0 |
| gap | two scratch units in sequence with a 2 s gap, `--grace 5` | `unit_seen`, `unit_gap` "no match for 2s", then `unit_gone` | 0 |
| failure (ticket) | `--unit kit108-sampler-probe.service` (`sleep 20`), stopped after 5 s | `unit_seen`, process rows, `unit_gone` "after sampling" | 0 |
| failure (already stopped) | `--unit kit108-sampler-stopped.service`, a unit that ran `true` and is `inactive` | `unit_gone` "not active at start" | 0 |
| failure (never seen) | `--pattern 'kit108-nothing*' --wait 3` | `never_seen` | 0 |

### I3 on a finished screen member's project (`107a`: `lipid_screen_missing` on `kit_107`; R 4.6.1; no KC data)

- **Completion** (`lifecycle/i3-107a-complete/`): exit 0, all 16 rows through `end complete`. The copy wrote `DA_results_lipidomics.xlsx` and `Pathway_enrichment_results_lipidomics.xlsx`; the render wrote `lipidomics_report_lipidomics_dia_baseline_smoke.docx`. The profile tables are non-empty and the object-size table is written. What it recorded (kB, MB):

  | Row | `RssAnon` kB | gc used Vcells MB | gc max used Ncells / Vcells MB | mallinfo2 in-use / free / mmap bytes |
  |---|---|---|---|---|
  | `package_loaded` | 368,456 | — | — | 349.2 M / 14.1 M / 16.6 M |
  | `before_copy` | 368,828 | 38.4 | 161.3 / 38.4 | 340.8 M / 23.3 M / 16.6 M |
  | `after_copy_raw` (0.31 s) | 368,916 | — | — | 349.6 M / 14.4 M / 16.6 M |
  | `after_copy` | 368,916 | 39.4 | 255.7 / 52.1 | 316.1 M / 47.9 M / 16.6 M |
  | `before_render` | 368,916 | 39.6 | 164.8 / 39.6 | 316.2 M / 47.8 M / 16.6 M |
  | `after_render_raw` (1.19 s) | 458,696 | — | — | 337.3 M / 28.8 M / 105.5 M |
  | `after_render` | 371,856 | 43.4 | 255.7 / 131.9 | 334.2 M / 31.9 M / 16.6 M |
  | `after_malloc_trim` (`trim_return` 1) | 350,420 | — | — | 334.3 M / 29.0 M / 16.6 M |

  At this size:
  - the copy leaves `RssAnon` unchanged (+0.1 MB) with a transient of +94 MB Ncells and +14 MB Vcells;
  - the render's +88 MB of raw `RssAnon` is mmapped R vectors (mallinfo2 mmap 16.6 → 105.5 MB) that `gc()` returns (+3 MB left);
  - `malloc_trim(0)` returned 21 MB;
  - `.GlobalEnv` holds only `project_dirs` (6.5 kB) and `.Random.seed`;
  - the return values are 48 B and 264 B;
  - the profile's top allocators are `openxlsx::writeData` in the copy (77 MB allocated) and `knitr::knit` in the render (195 MB allocated).
- **Failure** (`lifecycle/i3-107a-fail/`): the run used `--remove-template --template-dir <missing dir>`.
  - It exited 1 and wrote every row up to the failure: the copy rows, then `template_missing` (warn: "not in scripts/ and not at …/no-such-reports-dir/lipidomics_report.rmd"), then `before_render`.
  - The final row is `end failed` "in render: R Markdown template file not found at the expected location …", written by `on.exit`.
  - The copy's profile tables were written. The render's raw `rprof-render.out` holds only its header (56 bytes, no samples), and no render summary was written.
- **`--load-all`** (`lifecycle/i3-107a-load-all/`): on an exported tree of HEAD `50ae3054` under `TMPDIR`, it completed (exit 0), and the profile names file:line:
  - `func_general_filemgmt_results_summary_helpers.R#665` (the `writeData` into the combined workbook) and `#702` in the copy;
  - `func_general_filemgmt_report_render_helpers.R#212` (`rmarkdown::render`) in the render.

  An installed build has no srcrefs, so on `kit_107` the profile names functions, not lines (`<no location>`).

**Still to prove (M1, the lead):** I2 on the GUI member `108d` against D7's criteria, and I3 complete and forced-fail on `108d`'s own project.

## Commands for the lead

These run under the host rule: nothing else heavy alongside, and no other `multischolar*` unit active. Use `TMPDIR` on disk.

```
R=/home/doktersmol/Projects/APAF/MultiScholaR; P=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5
K=$R/docs/qa/workflow-execution/kit/kit-108/memory; L=$P/core/browser-smoke.66sPmE/library   # kit_107
O=/home/doktersmol/.cache/k108-mem; mkdir -p $O; export TMPDIR=/home/doktersmol/.cache/k108w2
systemctl --user list-units --state=active --no-legend 'multischolar*'      # must print nothing
```

### M1: the lifecycle member `108d`, with I2 attached

```
bash $K/sampler.sh $O/108d-i2.tsv --pattern 'multischolar*' --wait 900 --grace 30 --max 3600 > $O/108d-i2.log 2>&1 &
bash $R/docs/qa/workflow-execution/kit/kit-100/scripts/member.sh kit-journey-108d $L 2400
wait            # the sampler ends by itself on unit_gone (exit 0) about 30 s after the slice stops
D=$P/kit/runtime/kit-journey-108d/candidate-1
R_LIBS=$L Rscript $K/analyse-traces.R $O/108d-i1.tsv $D                   # I1: save, copy, render rows present
R_LIBS=$L Rscript $K/analyse-sampler.R $O/108d-i2.tsv $D $O/108d-i2         # I2: the M1 checks and windows
# I3 on 108d's project, once the member has ended:
R_LIBS=$L Rscript $K/replica.R --project $D --out $O/108d-i3 > $O/108d-i3.log 2>&1; echo "exit $?"   # expect 0
R_LIBS=$L Rscript $K/replica.R --project $D --out $O/108d-i3-fail --remove-template \
    --template-dir $O/no-such-reports-dir > $O/108d-i3-fail.log 2>&1; echo "exit $?"                  # expect 1
tail -1 $O/108d-i3-fail/replica-rows.tsv      # expect phase end, status failed
rm -rf $O/108d-i3-fail/work
```

**D7's M1 criteria, read from the outputs:**
- **I1:** `$O/108d-i1-summary.tsv` has save, copy and render with positive seconds.
- **I2:** `$O/108d-i2-check.tsv` must show:
  - `max tick gap inside the report step` ≤ 2;
  - `controller R PIDs` names one PID;
  - `rows within 10% of controller_anon` = all rows;
  - `events` ending in `unit_gone`.
- **I3:** `108d-i3/replica-rows.tsv` ends `end complete`, `rprof-copy-by-total.tsv` and `rprof-render-by-total.tsv` are non-empty, and `object-sizes.tsv` exists. The forced-fail run ends with an `end failed` row.

If an instrument fails, correct it once (an instrumentation failure, D7) and rerun the case.

### M2: the baselines `108p` and `108q`, with I2 attached

Run them one after the other, the same way:

```
for k in 108p 108q; do
  bash $K/sampler.sh $O/$k-i2.tsv --pattern 'multischolar*' --wait 900 --grace 30 --max 16000 > $O/$k-i2.log 2>&1 &
  bash $R/docs/qa/workflow-execution/kit/kit-100/scripts/member.sh kit-journey-$k $L 7200
  wait
  D=$P/kit/runtime/kit-journey-$k/candidate-1
  R_LIBS=$L Rscript $K/analyse-traces.R $O/$k-i1.tsv $D
  R_LIBS=$L Rscript $K/analyse-sampler.R $O/$k-i2.tsv $D $O/$k-i2
done
```

`member.sh`'s third argument is its observation wait, not the member's ceiling; the binding's 14,400 s ceiling governs. If the wait ends first, inspect the live job and do not restart it. I2's `--max 16000` stays above the ceiling.

### M3: the replica, twice on each baseline and once on `099l`

After M2, with no member live:

```
for k in 108p 108q; do for r in 1 2; do
  R_LIBS=$L Rscript $K/replica.R --project $P/kit/runtime/kit-journey-$k/candidate-1 --out $O/$k-i3-r$r > $O/$k-i3-r$r.log 2>&1; echo "$k r$r exit $?"
done; done
R_LIBS=$L Rscript $K/replica.R --project $P/kit/runtime/kit-journey-099l/candidate-1 --out $O/099l-i3 > $O/099l-i3.log 2>&1; echo "099l exit $?"
```

- **Optional, for file:line:** add `--load-all <tree>`, where the tree is a `git archive` export of the start commit under `TMPDIR` (`R/` there is identical to `kit_107`'s source, per the ticket's start condition). The process baseline is higher under `load_all` (`RssAnon` 508,948 kB against 368,456 kB on `107a`), so compare rises, not levels.
- **The two repetitions must agree within 10%** (D7) on the copy's and the render's raw `RssAnon` rise.
- **The validity rule (D7).** Compute the replica's rise and the GUI member's controller rise across the same two windows, and require them to be within ±25%:
  - replica: `after_render_raw.rss_anon_kb − before_copy.rss_anon_kb`;
  - GUI member: the controller's anon at the render's exit minus at the copy's entry, in `$O/<k>-i2-windows.tsv` (exact clocks) or `$O/<k>-i1-summary.tsv` (bracketed).

  ```
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} $3=="before_copy"{a=$h["rss_anon_kb"]} $3=="after_render_raw"{b=$h["rss_anon_kb"]} END{printf "replica rise %.3f GiB\n", (b-a)/1048576}' $O/108q-i3-r1/replica-rows.tsv
  ```

- **Disk:** each large replica copies about 2.2 to 3.0 GB with `cp -a`, a reflink on `/home`, and removes it on completion.

## What each candidate cause should look like in M2 and M3

These are the readings that confirm a cause under D7's attribution method. Each names the rows to read.

| Cause (D6) | I2 in M2 (GUI member, exact clocks) | I3 in M3 (replica) | Fix class (D6's table) |
|---|---|---|---|
| **H3: a live reference held after the step** | The R PID's `RssAnon` rises in the copy (or render) window and stays to teardown. No child accounts for it. | `after_copy`'s gc **used** (Ncells/Vcells MB) stays up by about the raw `RssAnon` rise. `final_gc` keeps it. `object-sizes.tsv` names the holder: a `.GlobalEnv` object or a large return value. `malloc_trim` returns little. If the replica shows no such hold but the GUI member does, the holder is application-side (a reactive, output or bus value the replica does not have). Say so; the replica then fails the validity rule for this cause. | drop or scope the reference: allowed |
| **H1 + H2: a transient in the copy, kept by the allocator** | The R PID's `RssAnon` climbs through the copy (stepwise, as in M0's 084d shape) and stays flat after it. No child. The worker units' anon is unchanged. | High **max used** in `after_copy` (Vcells/Ncells well above `before_copy`) and low **used** after `gc()`. `after_copy_raw` ≈ `after_copy` `RssAnon`: gc frees the heap but RSS stays. mallinfo2 `free_arena_bytes` ≈ the retained rise (or `VmData` well above arena + mmap if non-main arenas hold it). `after_malloc_trim` drops `RssAnon` by most of the rise (`trim_return` 1). The copy profile's top allocators are `openxlsx::read.xlsx`/`writeData`/`saveWorkbook` (with `--load-all`: `func_general_filemgmt_results_summary_helpers.R:646-702`). | reduce the transient with identical output (one sheet at a time, release, `gc()`): allowed if the `results_summary` files are identical. If only returning memory helps (`malloc_trim`), it is parked (P-30). |
| **H4: render-time memory held after `RenderReport`** | The R PID rises in the render window and stays. | `after_render` (after gc) `RssAnon` and gc used stay above `before_render` by ≥ 0.1 GiB. `object-sizes.tsv` shows objects left in `.GlobalEnv` by the render (e.g. the template's tables), or the rise sits in mallinfo free bytes (then it is H2 for the render). | clear the render environment, `gc()`: allowed |
| **A child process (pandoc) counted in the controller's cgroup** | `-children.tsv` shows a non-R PID (`pandoc`) in the controller cgroup. `-windows.tsv` shows `children_anon_peak` in the render window, and the cgroup's anon minus the R PID's `RssAnon` ≈ that child's. After the child exits, the cgroup's anon falls back. The R PID's own `RssAnon` does not hold the rise. | not visible (the replica records only its own process) | none: a measurement artefact; the ticket says so with I2's rows |
| **Unattributed** | The R PID's `RssAnon` rise is matched by none of the above. | Neither gc used, nor mallinfo free, nor `malloc_trim` accounts for it. `VmData` − (arena + mmap) shows how much lies outside the main arena that mallinfo2 reports. | none: the diagnosis says so, with its evidence |

**From M0, before any of this:** the rise is in the copy, at least +1.03 GiB on every large member at any alignment, and it is kept to teardown. So H1 + H2 and H3 are the live candidates for the copy, H4 for the smaller render share (≤ +0.39 GiB), and the save window's share (−0.19 to +0.55 GiB) is unresolved. I2's exact clocks decide that share in M2.

## Limits the diagnosis must state

- **mallinfo2 covers the main arena** (plus `hblkhd`, the mmapped chunks). Memory in other glibc arenas (allocations from other threads) is not in it. Compare with `VmData`/`RssAnon`.
- **The replica's starting level is the package alone** (about 0.37 GiB), not the application's 1.2 to 1.4 GiB of workflow state. D7's validity rule compares rises, not levels.
- **The first call compiles closures** (`cmp`/`cmpCall` frames in the copy profile); that is profile noise on a 0.3 s copy, and small against a 47 to 106 s copy.
- **`Rprof`'s `mem.total` is allocation volume,** not the peak. The peak is gc's max used, and RSS is `/proc/self/status`.
