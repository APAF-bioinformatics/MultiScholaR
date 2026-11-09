# KIT-108: the report-step application memory, diagnosed (2026-11-09)

**Question** (the user, 2026-11-07): "The report step is unexplained. Lipid/metab app memory sits at 3.33 GiB at the report." The user asked for a thorough evaluation and diagnosis, followed by a fix.

**Answer.** The rise happens inside "Copy to publication", not the report render. For each assay, the copy re-reads the large DA results workbook with `openxlsx::read.xlsx` and writes it into a new workbook with `writeData`. That builds a transient R heap of about 2.1 GB. R frees it afterwards, but glibc's allocator keeps the freed pages instead of returning them to the system, so the application stays about 1.6 to 1.9 GiB higher until it exits. Only about 0.1 GiB stays live. Transcriptomics copies one 2.6 MB csv, so it stays flat.

## Sources

**Runs.** The lead ran M1 to M3 on `kit_107` (R 4.6.1, the start commit's build), through `$P/kit/scripts-108/memory-m123.sh` and the instruments in [instruments.md](instruments.md):
- I1, the harness's 50 ms controller trace, split by the driver's marks;
- I2, the per-process side sampler;
- I3, the in-process replica of the copy and the render.

**Raw outputs** are in `/home/doktersmol/.cache/k108-mem`. Reduced copies are in [m1-m3/](m1-m3/); the raw I2 hashes are in `m1-m3/i2-raw-sha256.txt`; validity is in [m3-validity.tsv](m3-validity.tsv), from `validity.R`. The members' evidence was committed at `b97ad9e8` (`108p`) and `f4e0fe2d` (`108q`), and the analysis at `b529f4cf`. The analysis was done by the W2 agent, and this file was written by the lead from its report.

## 1. The instruments' lifecycle (M1): passed

- **I1** found all three windows on `108d` (save 5.63 s, copy 0.26 s, render 1.77 s).
- **I2** on `108d`:
  - largest sample gap 1.12 s (1.11 s inside the report step);
  - one controller R PID;
  - 173/173 rows within 10% of the harness's controller anon (max relative error 0.0058);
  - events: start, unit_seen, unit_gone.

  On `108p` and `108q`: gaps of at most 1.19 s, 1,136/1,136 rows within 10%, one R PID each.
- **I3** on `108d`'s project completed with every row, a non-empty profile and the object-size table. The forced failure ended with an `end failed` row ("template not found").
- **D7's one instrumentation correction.** The I2 reader paired the R process with its cgroup by interpolation, which read up to 15% off at step changes. It now uses the same tick's cgroup row. Neither the windows nor the D7 check changed.

## 2. The replica's validity (D7)

| Member | Repeats r1 vs r2 | Replica rise vs GUI rise (copy entry to render exit) | Valid |
|---|---|---|---|
| `108p` metab_large_missing | copy 0.000%, render 0% | 1.660 vs 1.829 GiB (−9.3%; the raw rises without the replica's gc are within +2.1%) | yes |
| `108q` lipid_large_missing | copy 0.021%, render 0% | 1.931 vs 1.246 GiB (+55%) | **no**: not attribution evidence |

The lipid diagnosis therefore rests on I1 and I2. It shows the same mechanism, but its size is not validated by the replica. A possible reason, which was not measured: the replica starts from an empty heap, while the application can reuse free memory it already holds.

## 3. Where the rise happens (I2, the R process's RssAnon, GiB, exact clocks)

| Window | `108p` (metab) | `108q` (lipid) |
|---|---|---|
| Before the report | 1.308 | 1.398 |
| Save | 1.308 → 1.115 | → 1.423 |
| **Copy** | **→ 2.879** (49.6 s) | **→ 2.591**, peak 2.788 (43.0 s) |
| Render | +0.065 (→ 2.944) | +0.078 (→ 2.669) |
| After the report, to teardown | flat for 30 s, 2.945 | flat for 27 s, 2.669 |
| DA-step peak / report-step peak | 1.540 / 2.980 | 1.398 / 2.815 |
| Teardown / slice `resident_end` | 2.945 / 3.304 | 2.669 / 3.029 |

- The save adds nothing. M0's first jump at the save/copy boundary is the copy starting.
- No child process carries memory: only `sh`, `bash` and `sync` appear, each under 0.3 MB, and none during the render.

## 4. The mechanism (I3)

| Measure | `108p` | `108q` (not validated) | `099l` trx |
|---|---|---|---|
| gc used before / max during / after the copy (MB) | 199.7 / 2,124.5 / 299.4 | 199.7 / 2,115.2 / 299.4 | 199.7 / 332.9 / 207.0 |
| RssAnon raw (GiB) | 0.352 → 2.217 | 0.352 → 2.282 | 0.352 → 0.384 |
| RssAnon after `gc()` | 2.010 | 2.282 | 0.384 |
| mallinfo2 in use / free before the render (GiB) | 0.340 / 1.596 | 0.340 / 1.868 | 0.310 / 0.061 |
| `malloc_trim(0)`: final gc → after trim (GiB) | 2.012 → 0.453 | 2.284 → 0.453 | 0.384 → 0.335 |
| Left above the level before the copy (GiB) | +0.101 | +0.101 | −0.017 |
| Copy time | 52.9 s | 52.8 s | 1.3 s |

**The allocations, from Rprof by function** (the installed build gives function names):
- `openxlsx::read.xlsx`: 2,456 MB in total, of which `getCellInfo` 1,352 MB;
- `writeData`: 1,663 MB;
- the rest: `unlist`, `as.matrix`/`t`, `write_worksheet_xml`, `read_workbook`.

All of it sits under the DA-workbook loop of `copyToResultsSummary`, in `R/func_general_filemgmt_results_summary_helpers.R` (`:653`; `read.xlsx` at `:661`, `writeData` at `:665`, `saveWorkbook` at `:702`). A `--load-all` run named `:665` and `:702`.

## 5. The causes, ranked

1. **H1 + H2, dominant: the copy's transient heap, kept by the allocator.** Gc max-used reaches 2.1 GB while R holds only 0.3 GB afterwards. glibc's main arena then holds 1.6 to 1.9 GiB free, and `malloc_trim(0)` returns 1.56 GiB of `108p`'s 1.66 GiB rise (94%). The copy reads and rewrites two 44 to 46 MB xlsx files per assay into an 85 to 87 MB workbook. Transcriptomics reads one 2.6 MB csv (heap peak 333 MB).
2. **H3, a small live residue: +0.10 GiB** of gc-used memory remains after the copy and the final `gc()`. `.GlobalEnv` holds only `project_dirs` (6.5 kB), and the return values are 48 and 264 bytes. The holder was not identified.
3. **H4 or the allocator: the render adds +0.065 / +0.078 GiB,** in the R process. That is below the screen members' render rise (+0.09 to +0.12 GiB).

## 6. The fix routes (the ticket's decision table)

| Route | What | Class | Expected effect |
|---|---|---|---|
| **A** | Build the DA workbook one sheet at a time: the `purrr::imap` becomes a `for` loop with the same index, sheet names and failure records. `rm(data_content); gc(FALSE)` after each `writeData`, and `rm(da_wb, enrichment_wb); gc(FALSE)` after the saves. Order, content and defaults stay the same. | **Allowed** (W5's one mechanism) | Not measured, and partial. It removes at most one sheet's read garbage overlapping the next. The workbook's own cell store stays until the save. It is unlikely on its own to meet D8's thresholds: peak ≤ DA peak + 0.25 GiB, about 1.79 / 1.65 GiB; teardown ≤ pre-report + 0.25 GiB, about 1.56 / 1.65 GiB. |
| B | `malloc_trim(0)` after the copy (native code) | Parked (P-30) | In the replica it returns 94% of `108p`'s rise, so teardown would be fixed. The report-step peak stays about 2.9 GiB. |
| C | The copy on the worker (a new worker action) | Parked (P-30) | Fixes both peak and teardown. A `callr` child would not help, because it stays in the controller's cgroup. |
| D | A lighter xlsx reader | Not allowed as written | `readxl` would be a new dependency. Not measured. |

**Route A's test and guards** (W5):
- **The test:** in `test-omics-summary-core.R`, on a synthetic DA directory (2 assays × 2,000 features × 200 samples, as xlsx), gc's max used over the copy stays below the cost of holding both sheets' data at once, and nothing larger than one sheet survives the call. It fails at W4 and passes at W5.
- **The guards:**
  - every `results_summary` file is identical sheet by sheet (`sha256` per sheet) between the W4 and W5 trees on `108p`, `108q` and `099l`;
  - the report text is identical, with the timestamp masked;
  - W1 (a)/(b) and the ten lanes are identical;
  - M4's states and products equal the baselines'.

**Recommendation (adopted by the lead):** W5 tries route A once and measures it with the replica at W4 and W5, then with `108v`/`108w` (M4). If it misses D8's threshold, the stop rule applies and P-30 is raised with routes B and C and these numbers. No second mechanism is tried.

## Note on the run conditions

While `108p` (the metab baseline) was live, a W3 agent ran a scratch test for about 3.5 minutes that grew to about 1.9 GB, against the host rule ([../findings.md](../findings.md)). The figures above are the application's own process (I2 per PID, cross-checked against the harness's controller trace), and `108p`'s replica reproduces its copy rise within 10%. So the attribution stands, but `108p`'s wall-clock timings were not taken alone.
