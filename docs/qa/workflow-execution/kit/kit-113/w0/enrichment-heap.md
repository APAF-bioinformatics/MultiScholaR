# KIT-113 W0 (M0p): the enrichment reachability probe

**Decision (D8): "better" is UNREACHABLE at W0.**

On the decision case, the managed controller's figure is higher than the in-process figure. The live `processEnrichments` over the e2e doubles at 16,000 proteins × 300 samples × 3 contrasts, pass 2 (warm), gives:

| Figure | GiB |
|---|---|
| in process | 1.0751 |
| managed controller | 1.7319 |
| difference (in process − managed) | **−0.6568** |

The threshold is +0.25 GiB. Pass 1 (cold) decides the same way: 1.2321 − 2.1787 = −0.9466 GiB. Contingent P-43 is raised (the lead's row).

| Measure | Rule | Result |
|---|---|---|
| Better is reachable | in-process figure − managed-controller figure ≥ 0.25 GiB on the decision case | −0.6568 GiB (warm); −0.9466 GiB (cold): **not met** |
| Otherwise | "better" recorded **unreachable at W0**; contingent P-43 raised; the move proceeds and is judged on the no-worse rules only | **applies** |

**Why the move cannot lower the peak:** the live kernel's results carry plots, so the managed controller would hold more than the in-process one.

- **The results do not travel small.** The portable results serialize to 186.8 MiB for 3 contrasts and 64.2 MiB for 1. Hydrated from the `rds` state value, they are a 535 MiB object in the controller (`lobstr::obj_size`). The in-process live object is 16 MiB, because it shares its environments with the loaded packages.
- **The managed controller's peak is the commit (A8):** about +1.17 GiB over a heap already holding +602 MiB of hydrated results.
- **The in-process peak is the same commit (G9e).** `portableWorkflowEnrichmentResults` forces the live object's lazy bindings, then the same value is serialized.
- **What does move is smaller than the commit:** the DA tables (155 MiB in R at 16,000 × 3, 491 MB on disk), the annotation match (+571 MiB over its start) and the kernel (+667 MiB).

## The run

- **Counted attempt:** started 15:18:45, finished 15:47:45 (2026-11-10). 1,744 s inside the script, 29 min 4 s wall. Exit 0.
- **Output:**
  - `enrichment-heap.txt` (the report);
  - `enrichment-heap-groups.tsv` (every group row, written by the script);
  - the kernel's and the DA's console output stayed in the scratch log.
- **Completed:**
  - 16,000 × 300 × 3, live (the decision case) and double;
  - 4,000 × 300 × 3, live and double;
  - 16,000 × 300 × 1, live and double.
- **Not started:** 4,000 × 300 × 1. The ceiling guard stopped it at 1,744 s ("CEILING: 1744 s elapsed; case 4000 x 300 x 1 not started").
- **A launch before the attempt:** the first launch command wrapped the run in `/usr/bin/time`, which this host lacks. The shell exited 127 at 15:18:25 before `Rscript` started: no R process ran and no output was written. The same command without the wrapper was then launched at 15:18:41 (the R script started at 15:18:45). That second launch is the one counted attempt.
- **The DA setup** (in-process DA, then its publication):

  | Owner | DA | Publication | Published DA tables |
  |---|---|---|---|
  | 16,000 × 300 × 3 | 218.7 s | 34.6 s | 3 × 122.6 MB |
  | 4,000 × 300 × 3 | 252.2 s | 10.5 s | 3 × 30.7 MB |
  | 16,000 × 300 × 1 | 110.6 s | 13.9 s | 182.2 MB |

- **Every sequence committed** `enrichment_completed` (saved = TRUE).

## The two sequence figures, decision case (16,000 × 300 × 3, live kernel over the doubles)

Each figure is the sequence's highest "max used" minus the live heap at its start (Ncells × 56 + Vcells × 8 bytes). The managed sequence's worker rows (W1-W5) run in this process but are left out of its figure.

| Sequence | Live heap at start | Pass 1 (cold) | Pass 2 (warm) | Peak group |
|---|---|---|---|---|
| in process (G1-G9f, today's order) | 882.9 / 845.3 MiB | 1261.7 MiB = 1.2321 GiB | 1100.9 MiB = **1.0751 GiB** | G9e, the state commit |
| managed controller (R1-R6, then A1-A9) | 846.9 / 847.5 MiB | 2231.0 MiB = 2.1787 GiB | 1773.4 MiB = **1.7319 GiB** | A8, the state commit (A1, the hydration, is next: 1386.2 / 1833.7 MiB) |
| difference | | −0.9466 GiB | **−0.6568 GiB** | |

## The per-group table, decision case

"Peak over start" is the group's "max used" minus the sequence's live heap at its start. "Over the group start" is the same "max used" minus the heap live just before the group: it is the group's own transient. Group rows cannot be added.

**In process** (`runWorkflowEnrichment`'s body, `utils_workflow_enrichment.R:822-925`, copied verbatim; `[moves]` marks D8's "moves" groups):

| | Group | Peak over start, warm (MiB) | Over the group start, warm | Live after, warm | Value (MiB) | Seconds, warm | Peak over start, cold |
|---|---|---|---|---|---|---|---|
| G1 | context (:822-845) | 45.6 | 45.2 | 0.7 | 4.52 | 0.1 | 33.7 |
| G2 | DA tables + refusal (:846-852) [moves] | 295.2 | 294.5 | 153.6 | 154.85 | 0.9 | 293.1 |
| G3 | annotation resolve (:853-854) | 167.0 | 13.5 | 153.6 | 3.05 | 0.0 | 156.1 |
| G4 | annotation match (:855-856) [moves] | 724.8 | 571.2 | 155.2 | 3.42 | 0.9 | 734.4 |
| G5 | organism filter, bus metadata (:857-876) [filter off] | 155.2 | 0.0 | 155.2 | 0.00 | 0.0 | 157.7 |
| G6 | progress 0.5, method (:877-880) | 159.3 | 4.1 | 155.2 | 0.00 | 0.0 | 158.3 |
| G7 | id column, kernel args, kernel, cp10 (:881-889) | 170.4 | 15.2 | 155.2 | 155.95 | 0.0 | 157.8 |
| G8 | config$run = runDIAReferenceEnrichment (:890) [moves] | 821.7 | 666.5 | 162.3 | 15.99 | 2.8 | 839.1 |
| G9a | display tables (:891-900) | 225.2 | 62.9 | 162.4 | 0.02 | 4.4 | 160.4 |
| G9b | bus payload (:901-907) | 162.4 | 0.0 | 162.4 | – | 0.0 | 123.2 |
| G9c | results' @args stamped (:908-911) | 170.1 | 7.8 | 162.4 | 15.98 | 0.0 | 123.2 |
| G9d | UI parameters (:912-915) | 171.3 | 8.9 | 162.4 | 2.26 | 0.0 | 123.3 |
| G9e | **state commit** (:916-917) | **1100.9** | 938.5 | 162.7 | 0.00 | 6.2 | **1261.7** |
| G9f | tab_status, progress, the returned list (:918-925) | 162.6 | 0.0 | 162.6 | 19.32 | 0.0 | 123.6 |

**Managed controller** (D3's groupings over today's statements and functions; W rows are the emulated worker, not counted):

| | Group | Counted | Peak over start, warm (MiB) | Over the group start, warm | Live after, warm | Value (MiB) | Seconds, warm | Peak over start, cold |
|---|---|---|---|---|---|---|---|---|
| R1 | context (:822-845) | yes | 29.0 | 29.0 | 0.3 | 4.52 | 0.1 | 29.0 |
| R2 | annotation resolve (:853-854) | yes | 0.3 | 0.0 | 0.3 | 3.05 | 0.0 | 0.3 |
| R3 | mapping and its warning [filter off] | yes | 0.3 | 0.0 | 0.3 | 0.00 | 0.0 | 0.3 |
| R4 | progress 0.5, method (:877-880) | yes | 0.8 | 0.5 | 0.3 | 0.00 | 0.0 | 0.8 |
| R5 | the kernel's name | yes | 0.9 | 0.6 | 0.3 | 0.00 | 0.0 | 0.9 |
| R6 | transfer values: the annotation table as a saved display table | yes | 5.9 | 5.7 | 0.3 | 0.14 | 0.6 | 18.4 |
| W1 | source components, annotations materialized | worker | 47.4 | 47.2 | 1.7 | 5.31 | 0.1 | 47.5 |
| W2 | DA tables + refusal | worker | 290.6 | 288.9 | 154.5 | 154.85 | 0.9 | 290.6 |
| W3 | annotation match | worker | 1184.7 | 1030.2 | 156.0 | 3.42 | 0.8 | 722.9 |
| W4 | kernel inputs, run, product, roots (computeDIAReferenceEnrichment, `utils_workflow_dia_da.R:657-676`) | worker | 1561.1 | 1405.1 | 162.6 | 15.99 | 2.6 | 851.9 |
| W5 | portable results saved as an `rds` state value | worker | 1773.4 | 1610.8 | 162.8 | 19.31 | 4.4 | 1099.4 |
| — | the moved values removed and collected: live heap +2.1 MiB over start; delivered value 3.44 MiB | | | | | | | |
| A1 | **results hydrated** from the `rds` state value | yes | 1386.2 | 1384.1 | **602.4** | **535.42** | 3.4 | 1833.7 |
| A2 | provenance cleared; external-resource records | yes | 602.4 | 0.0 | 602.4 | – | 0.0 | 602.4 |
| A3 | bus filter metadata (:874-876) | yes | 602.4 | 0.0 | 602.4 | 0.00 | 0.0 | 602.4 |
| A4 | display tables (:891-900) | yes | 635.0 | 32.6 | 602.4 | 0.02 | 4.7 | 634.7 |
| A5 | bus payload (:901-907) | yes | 602.4 | 0.0 | 602.4 | – | 0.0 | 602.4 |
| A6 | results' @args stamped (:908-911) | yes | 602.4 | 0.0 | 602.4 | 535.42 | 0.0 | 602.4 |
| A7 | UI parameters (:912-915) | yes | 602.5 | 0.1 | 602.4 | 2.26 | 0.0 | 602.5 |
| A8 | **state commit** (:916-917) | yes | **1773.4** | 1171.0 | 602.4 | 0.00 | 7.5 | **2231.0** |
| A9 | tab_status, progress, the returned list (:918-925) | yes | 602.4 | 0.0 | 602.4 | 538.87 | 0.0 | 602.4 |

## Every case (reported, not judged)

| Case | Kernel | In process, cold / warm (GiB) | Managed controller, cold / warm (GiB) | Difference, cold / warm (GiB) | Portable results, serialized | In-process peak group | Managed peak group |
|---|---|---|---|---|---|---|---|
| **16,000 × 300 × 3** | **live over the doubles** | 1.2321 / **1.0751** | 2.1787 / **1.7319** | −0.9466 / **−0.6568** | 186.82 MiB | G9e commit | A8 commit |
| 16,000 × 300 × 3 | deterministic double | 1.1556 / 1.1556 | 0.0774 / 0.0784 | 1.0782 / 1.0772 | 0.06 MiB | G4 match | A8 commit |
| 4,000 × 300 × 3 | live over the doubles | 1.7170 / 1.7338 | 1.7335 / 1.7329 | −0.0165 / 0.0009 | 186.82 MiB | G9e commit | A8 commit |
| 4,000 × 300 × 3 | deterministic double | 0.2931 / 0.2931 | 0.0761 / 0.0771 | 0.2170 / 0.2160 | 0.06 MiB | G4 match | A8 commit |
| 16,000 × 300 × 1 | live over the doubles | 1.9092 / 1.6142 | 1.6141 / 1.6139 | 0.2951 / 0.0004 | 64.17 MiB | G9e commit | A8 commit |
| 16,000 × 300 × 1 | deterministic double | 1.0800 / 1.0800 | 0.0495 / 0.0503 | 1.0305 / 1.0297 | 0.02 MiB | G4 match | A8 commit |
| 4,000 × 300 × 1 | — | not started (ceiling) | | | | | |

**What the cases show:**

- **The portable results are not data-scaled in proteins.** They serialize to 186.82 MiB for 3 contrasts at 400 (dry run), 4,000 and 16,000 proteins, and to 64.17 MiB for 1 contrast. They scale with the plots, about 31 MiB per direction and contrast:
  - one diagnostic at 400 × 300 × 1 found the whole size in `enrichment_plots` and `enrichment_plotly`;
  - each `ggplotly` widget's `x` is about 31 MiB;
  - each static ggplot's `plot_env` also holds that widget.

  A `ggplotly` of a 5-row data frame outside the package already serializes to 5.9 MiB, almost all of it `x$visdat`. So this is not a `load_all` artefact: the package namespace serializes as a 68-byte reference. The doubles return 5 fixed terms; a real g:Profiler result size was not probed (D8's limit).
- **The data-scaled parts are the moved ones.**
  - The DA tables: 155 MiB in R at 16,000 × 3, 39 MiB at 4,000 × 3.
  - The annotation match: +571 to +1,030 MiB transient at 16,000; +262 MiB at 4,000.
  - The kernel: +667 to +1,405 MiB transient at 16,000.

  With the double, whose results hold no ggplot, the in-process peak is the annotation match (G4). The managed controller's peak is under 0.08 GiB, so moving the enrichment would lower the peak by 1.03 to 1.08 GiB at 16,000 proteins and by 0.22 GiB at 4,000. With the live kernel, the commit of a ~187 MiB value dominates both paths, and the managed path additionally holds the 535 MiB hydrated copy.
- **The managed figure is never lower than the in-process one by the threshold with the live kernel, warm:** −0.66, +0.001 and +0.0004 GiB. One cold figure (16,000 × 300 × 1) is +0.295 GiB, against +0.0004 warm. That pass's start heap included 78 MiB the previous case left behind, which was freed during the sequence. This is not the decision case.
- **An order effect, stated:** the decision case's in-process commit (G9e) ran first in the process and peaked lower than the same commit later in the process.

  | In-process commit (G9e) | Transient over the group start |
  |---|---|
  | decision case, warm / cold | +938 / +1,138 MiB |
  | 4,000 × 300 × 3 | +1,730 MiB |
  | 16,000 × 300 × 1 | +1,574 MiB |

  In the later live cases the in-process and managed figures converge to the same value (1.73 and 1.61 GiB): both serialize the same value, and the in-process commit first forces the live object to the hydrated copy's size. Had the decision case's in-process figure been at that later level, the difference would be about 0, still below +0.25. The decision does not depend on this.
- **Run-to-run variance:** the annotation match's transient at 16,000 × 300 × 3 was +571 MiB in G4 and +1,030 MiB in W3 (warm) on the same inputs.

## The network

- **The proxy:** the lane's closed-port proxy was set (`http(s)_proxy` and upper case = `http://127.0.0.1:9`, `no_proxy` empty; `helper-gui-entry.R:275-276`).
- **The doubles:** `helper-e2e-enrichment-doubles.R` was sourced, and `install_enrichment_doubles()` replaced `gprofiler2::gost` and `clusterProfiler::enricher`.
- **The tracers:** tracers counted every call into the HTTP entry points:
  - curl: `curl_fetch_memory`, `curl_fetch_disk`, `curl_fetch_stream`, `curl_fetch_multi`, `curl_download`, `curl`;
  - RCurl: `getURL`, `getURLContent`, `postForm`, `curlPerform`;
  - httr2: `req_perform` (httr has none of these and calls curl);
  - `utils::download.file` and `base::url`;
  - the UniProt fallback `getUniprotAnnotations`.
- **Result: no network call.** Every HTTP and UniProt counter stayed at zero over the whole probe. The only counter that moved was the offline double `gprofiler2::gost`, called 56 times: 6 per live sequence with 3 contrasts, 2 per sequence with 1 contrast. The `enricher` double was not called: taxon 9606 takes the g:Profiler path. The annotation table came from the bus, so the UniProt fallback was not reached.

## The probe's limits (as D8 words them, then the probe's own)

- **D8:**
  - it measures R heap in one process at fixture scale, not a managed controller's anonymous memory;
  - g:Profiler's real result size (terms, not samples) is not probed.
- **This probe's own limits:**
  - **The halves are emulated.** The new halves do not exist yet, so the managed sequence is today's statements in D3's groupings.
  - **The worker shares the process.** It runs in this process: the store, the loaded namespaces and the heap's fragmentation are shared. The worker rows are kept out of the controller's figure.
  - **A collection before every group.** `gc(reset = TRUE)` runs between groups, so a real process without those collections can peak higher on both paths.
  - **Heap only.** Memory outside the R heap is not counted: vroom and arrow buffers, graphics devices, Glimma and widget writers.
  - **The organism filter was off** (as in the replay), so its application, a moved group, was measured as a no-op.
  - **One annotation table:** a synthetic per-protein table on the bus (16,000 rows, the test-mode table's columns).
  - **One taxon:** 9606, the g:Profiler path. The clusterProfiler path over the `enricher` double was not run.
  - **No technical replicates.** The fixture's design has one technical-replicate group per sample, so the DA fits without `duplicateCorrelation` (see the declaration). Enrichment reads only the published DA tables.
  - **Cold passes start with remnants.** A cold pass's start heap can include the previous case's remnants. The decision is taken warm.

## Declared before the counted attempt (written 2026-11-10, before the run)

- **One attempt.** `Rscript docs/qa/workflow-execution/kit/kit-113/w0/enrichment-heap.R . docs/qa/workflow-execution/kit/kit-113/w0/enrichment-heap.txt counted`.
  - **Build:** HEAD 2bad7519, whose `R/` equals the start commit 6dce329d and `kit_112`'s source 5ecc9fd1. Loaded with `pkgload::load_all`.
  - **Environment:** R 4.6.1, `TMPDIR` on disk (`…/kit/tmp-113`). `R_LIBS` puts the browser-smoke library first, for `fastduckdb`, which the DA publication's Glimma read needs. The managed launcher's thread and duckplyr variables are set (`OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, `DUCKPLYR_FALLBACK_AUTOUPLOAD=0`, `DUCKPLYR_FALLBACK_COLLECT=0`; `utils_workflow_execution.R:1164-1165`).
  - **Ceiling:** 30 min wall, under `timeout 1800`. The script starts no new case or kernel after 1,740 s.
- **Order:** the decision case comes first: the live `processEnrichments` over the e2e doubles at 16,000 × 300 × 3. Then the deterministic double on the same owner. Then 4,000 × 300 × 3, 16,000 × 300 × 1 and 4,000 × 300 × 1, each with the live kernel and then the double.
- **Sequences per kernel:**
  - the in-process sequence twice: pass 1 is the first Run in the process ("cold"), pass 2 is a second Run ("warm");
  - then the managed-controller sequence twice.
- **Decision figure:** pass 2 (warm) of each sequence, because by then the namespaces the kernel loads are already in the heap at the start of both sequences. Pass 1 is reported too. If pass 1 would decide differently, that is stated.
- **Rule (D8):** "better is reachable" iff (in-process figure − managed-controller figure) ≥ 0.25 GiB on the decision case. GiB = 2^30 bytes.
- **Dry runs before the counted attempt (throw-away, at 400 proteins × 300 × 3):** three probe runs and one DA profile.
  - Dry run 1 failed: `fastduckdb` was missing from the default library.
  - Dry run 2 completed (533 s). Its DA took 155 s at 400 proteins with KIT-112 W1's technical-replicate pairs.
  - Dry run 3 completed (381 s) with one technical-replicate group per sample (no `duplicateCorrelation`). Its DA took 145 s.
  - The profile (1 contrast, 41 s) showed the time going to the store's table hashing through a subprocess pipe (`hashWorkflowSerializedValue`), not to the model.
  - The fixture keeps one technical-replicate group per sample: `duplicateCorrelation`'s per-protein REML (about 10 s per 400 proteins) would cost minutes at 16,000 inside the ceiling. The DA is setup, not measured.
  - No dry run used the declared sizes.
- **After the attempt, two throw-away diagnostics** (not probe reruns; they explain the portable size):
  - 400 × 300 × 1 with the live kernel: the serialized size by slot and by plot element;
  - a 5-row `ggplotly` outside the package.
