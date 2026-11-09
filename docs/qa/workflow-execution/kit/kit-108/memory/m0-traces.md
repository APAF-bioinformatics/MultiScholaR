# KIT-108 W2 · M0: the existing controller traces by action (I1)

Run 2026-11-09 by the W2 agent. No member was started. The run took seconds and read files only.

## What was run

```
P=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5
K=docs/qa/workflow-execution/kit/kit-108/memory
R_LIBS=$P/core/browser-smoke.66sPmE/library Rscript $K/analyse-traces.R $K/m0-windows.tsv \
  $P/kit/runtime/kit-journey-{084c,084d,086p,086q,086u,086v,099l,107a}/candidate-1
```

- **Inputs, per member:** `memory-trace.tsv` (the harness's 50 ms trace; `controller_anon` is the controller service's cgroup anon and `anon` is the whole slice's), `journey-phases.tsv` (step marks, uptime clock) and `driver-phases.jsonl` (the driver's enter/exit marks, epoch clock).
- **Outputs:**
  - `m0-windows.tsv`: every window, with the application's anon at entry, peak and exit, the slice's anon, the largest sample gap and the alignment error.
  - `m0-windows-summary.tsv`: one row per member, the source of the tables below.
- **Members:** the ticket's seven, plus `107a` (`lipid_screen_missing` on `kit_107`, the diagnosis build) as the screen reference on the start build.

## How the windows are cut, and how good the alignment is

- **The windows:**
  - **save:** from the report step's start to the copy's `enter` mark. It covers opening the tab, the label, the description and the save click. The save click has no driver mark.
  - **copy:** from `enter:action:<ns>summary-copy_to_publication` to its `exit` mark.
  - **render:** from the copy's `exit` mark to the report step's `complete`. It covers the `generate_report` click, which has no mark.
- **The clocks.** The step marks and the trace use uptime; the driver's marks use epoch. I1 brackets the offset: every driver mark must fall inside its own step. The bracket is `[max(last mark − step complete), min(first mark − step start)]` over all steps.
  - The midpoint is used, and the half-width is reported as the alignment error.
  - Every value that depends on a driver mark is also given at both ends of the bracket (`*_range` columns).
- **A check on a known clock.** For `107a`, which ran in the current boot, the bracket is 1794036493.317 to 1794036493.845. `/proc/stat` btime is 1794036493 (whole seconds), which is consistent with it.
- **The copy's duration is exact.** It is the difference of two driver marks on one clock. The bracket only decides where the copy sits inside the report step, which is what splits growth between save and copy, and between copy and render.

## The per-member table (application = controller cgroup anon, GiB)

"After save" is the anon at the copy's enter, and "after copy" is the anon at the copy's exit. Each is given at the midpoint, with its range over the alignment bracket in brackets.

| Member | Build | Corpus | Align ± s | DA-step peak | Before save (report start) | After save | After copy | After render (report complete) | Report-step peak | Teardown (app) | Slice `resident_end` |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `084c` | kit_084 | metab_large | 0.77 | 1.556 | 1.404 | 1.819 [1.415–1.821] | 3.115 [2.993–3.292] | 3.292 | 3.328 | 3.293 | 3.652 |
| `084d` | kit_084 | lipid_large | 4.09 | 1.539 | 1.244 | 1.767 [1.196–1.794] | 3.281 [3.116–3.281] | 3.282 | 3.328 | 3.282 | 3.641 |
| `086p` | kit_087 | metab_large_missing | 0.64 | 1.480 | 1.381 | 1.716 [1.393–1.786] | 3.002 [2.832–3.225] | 3.225 | 3.262 | 3.229 | 3.588 |
| `086q` | kit_087 | lipid_large_missing | 1.40 | 1.479 | 1.330 | 1.645 [1.141–1.841] | 2.871 [2.865–2.871] | 2.871 | 2.917 | 2.871 | 3.229 |
| `086u` | kit_087 | lipid_screen_missing | 0.26 | 0.870 | 0.870 | 0.879 [0.876–0.881] | 0.881 [0.879–0.880] | 0.975 | 1.008 | 0.975 | 1.321 |
| `086v` | kit_087 | metab_screen_missing | 0.14 | 0.810 | 0.810 | 0.823 [0.819–0.823] | 0.824 [0.823–0.826] | 0.840 | 0.862 | 0.830 | 1.176 |
| `099l` | kit_098 | trx_60k_200 | 0.01 | 1.491 | 1.443 | 1.443 | 1.443 | 1.443 | 1.486 | 1.443 | 1.798 |
| `107a` | kit_107 | lipid_screen_missing | 0.26 | 0.791 | 0.791 | 0.791 | 0.791 [0.791–0.792] | 0.913 | 0.939 | 0.915 | 1.261 |

**Durations (s):**

| Member | Save | Copy (exact) | Render |
|---|---|---|---|
| `084c` | 7.8 [7.1–8.6] | 50.1 | 0.8 [0.0–1.5] |
| `084d` | 15.1 [11.1–19.2] | 105.7 | 5.4 [1.3–9.5] |
| `086p` | 7.3 [6.6–7.9] | 47.3 | 0.6 [0.0–1.3] |
| `086q` | 8.5 [7.1–9.9] | 47.1 | 1.4 [0.0–2.8] |
| `086u` | 5.4 | 0.27 | 2.0 |
| `086v` | 5.3 | 0.26 | 1.6 |
| `099l` | 6.9 | 1.8 | 2.6 |
| `107a` | 5.1 | 0.26 | 1.8 |

**The rise per window (GiB): the minimum and maximum over the alignment bracket.**

| Member | Save | Copy | Render | Whole report step |
|---|---|---|---|---|
| `084c` | +0.01 to +0.42 | **+1.30 to +1.58** | 0.00 to +0.30 | +1.89 |
| `084d` | −0.05 to +0.55 | **+1.49 to +1.92** | 0.00 to +0.17 | +2.04 |
| `086p` | +0.01 to +0.41 | **+1.29 to +1.44** | 0.00 to +0.39 | +1.84 |
| `086q` | −0.19 to +0.51 | **+1.03 to +1.72** | 0.00 to +0.01 | +1.54 |
| `086u` | +0.01 | 0.00 | +0.09 to +0.10 | +0.11 |
| `086v` | +0.01 | 0.00 | +0.01 to +0.02 | +0.03 |
| `099l` | 0 | 0 | 0 (transient peak +0.04, released) | 0 |
| `107a` | 0 | 0 | +0.12 | +0.12 |

The table matches the ticket's D6 figures: the report-step peaks 3.33/3.33/2.92/3.26, the teardowns 3.28/3.29/2.87/3.23, and the `resident_end` values 3.64/3.65/3.23/3.59.

## What M0 already shows

1. **The growth is in the report step and is never released.**
   - On the four large lipid and metab members, the application holds 1.24 to 1.40 GiB at the report step's start and 2.87 to 3.29 GiB at teardown: +1.54 to +2.04 GiB.
   - The `post` window, from report complete to teardown, is flat on every member: 25 to 47 s with no fall (084d: 3.282 → 3.282 over 47 s).
   - So whatever the cause, the memory is kept (H2 or H3), not merely a transient peak.
2. **The copy carries most of it** (confirmed). Its rise is at least +1.03 GiB on every large member, at any alignment inside the bracket, and the copy runs 47 to 106 s.
   - The render adds at most +0.01 to +0.39 GiB, depending on alignment.
   - The save window's share is not resolved by these traces. It ranges from −0.19 to +0.55 GiB across the bracket.
3. **The first jump sits at the save/copy boundary.**
   - The 1 s shape of the trace shows one step of about +0.4 to +0.5 GiB about 8 to 12 s into the report step:
     - 084d at t = 12 s: 1.20 → 1.68 GiB;
     - 084c at t = 8 s: 1.41 → 1.82;
     - 086q at t = 8 s: 1.14 → 1.61.
   - In each case that falls inside the bracket of the copy's enter mark. The existing traces cannot say whether it is the save click's end or the copy's start.
   - I2's rows carry both clocks, so M2 settles it without a bracket. The ticket's "it shows none in 084d: 1.24 to 1.18–1.20 GiB" holds for the first 11 s of the step: 1.24 to 1.19 GiB.
4. **The copy's growth is stepwise, with partial returns.**
   - 084d rises 1.68 → 2.78 → 2.65 → 2.34 → 2.22 → 2.55 → … → 3.11 GiB between t = 12 and 95 s, then stays flat at 3.11 GiB from t = 95 to 116 s.
   - Over 88 to 116 s the slice's page cache (`file`) rises from 1.43 to 1.86 GiB, consistent with the workbook being written. That is an observation, not an attribution.
   - The render then adds 3.11 → 3.28 GiB at t = 117 to 119 s.
5. **Screen members.** On `086u`, `086v` and `107a` the copy is flat (≤ 0.004 GiB). The render adds +0.01 to +0.12 GiB, which stays to teardown. This is the "+0.11 to +0.15 GiB" of D8's threshold.
6. **The transcriptomics contrast (`099l`, 60,000 genes × 200 samples).**
   - The report step is flat: 1.443 GiB at its start, after the copy (1.8 s), after the render and at teardown.
   - The render's transient peak (1.486) is released. The RUVIII export, its own step before the report, is flat too: 1.443 → 1.443, peak 1.448.
   - So transcriptomics at a comparable sample count does not show the growth.
7. **The inputs differ in kind and size** (file listings of the members' projects, not an attribution):
   - The lipid and metab DA directories hold two `*_long_annot.xlsx` files of 44 to 46 MB each, next to 71 to 73 MB TSVs. The copy reads those xlsx files with `openxlsx::read.xlsx` and writes them into one workbook. The resulting `Publication_tables/DA_results_<omic>.xlsx` is 85.7 to 87.3 MB.
   - Transcriptomics reads one 2.6 MB `de_results_*.csv`, and its `DA_results_transcriptomics.xlsx` is 1.7 MB.
   - This fits H1 as a hypothesis. Only I3's profile and I2's per-process rows can attribute it.
8. **The slice above the application at teardown is about 0.36 GiB** on every member (`resident_end` − app teardown: 0.35 to 0.36 GiB). That is the other processes in the slice (the monitor, and workers if any are alive). It is about the same (0.35 to 0.36 GiB) on screen, large and trx members.

## What M0 cannot see

- The save/copy split of the first +0.4 to +0.5 GiB, because of the alignment.
- Whether the kept memory is live R heap (H3) or freed memory the allocator keeps (H2).
- Whether a child process (pandoc) contributes. The controller cgroup's anon includes children, and the render is short and adds little here.

M2 (I2's per-process rows, on exact clocks) and M3 (I3's gc, mallinfo2, malloc_trim and object sizes) decide these.
