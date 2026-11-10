# KIT-112 outcome: proteomics' DA kernel on the worker, measured (2026-11-10)

**Done.** On a managed session, proteomics' DA contrast computation (`DA-contrast-output`) now runs on the worker's `omics` role, one dispatch per contrast. The publication (`da_publication`) follows as before. Without a managed application, everything runs in process exactly as before.

**Results and records unchanged.**
- Every state, product, manifest, record, file and entry lane is identical.
- Nothing was re-frozen.

**The memory claim is made.** On the 300-sample DIA corpus, the application's `da` peak fell from 2.208 to 1.575 GiB ("better" needed ≤ 1.958). The kernel's +0.50 GiB rise in the application is gone: +0.015 on the worker path. Every "no worse" rule holds, and the DA step is faster (100.9 → 82.9 s) ([memory-outcome.md](memory-outcome.md)).

- **How it is built** (KIT-111's pattern):
  - `runDIAReferenceDaAnalysis` keeps its name and signature, and composes the request half (`diaDaRequest`, no matrix read), the compute half per contrast (`computeDIAReferenceDaContrast`, worker-safe), the accept half (`acceptDIADaManifest`) and today's publication, in today's order;
  - the worker's handler is `computeWorkflowDiaContrastAction`, registered with `named_input`;
  - the DIA binding adaptor gains a `da` entry;
  - the tab's run plugin dispatches through the existing omics runner on a managed session, behind two click guards (F112-02).
- **Defect fixed, not changing a method or a result:** F112-03, a DA rerun that met a saved Glimma page crashed the publication. It existed before this ticket and would have hit the worker path too. Fixed test first.
- **Build:** `kit_112` = `core/browser-smoke.OftnRb` (`5ecc9fd1`). `R/` +242/−89 lines in 4 files ([lines.tsv](lines.tsv)):
  - `utils_workflow_dia_da.R` +91 net (expected +50 to +110);
  - `mod_prot_da_handlers_run_analysis.R` +51 net (77 → 128), one line above the expected +20 to +50; it holds the managed dispatch and F112-02's two guards;
  - the others small.

## Work items

| Item | Commit(s) | What |
|---|---|---|
| W0 | `64e7c2ae` | Start condition met (`R/` = `kit_111`'s source). Large baseline = 111q (no 112p). Reachability: `da` 2.208 against `session_export` 1.708, so +0.50, and "better" is reachable ([memory/baseline-large.md](memory/baseline-large.md)) |
| W1 | `b388d082`, `525cbd78` | The in-process captures frozen on the start commit (`test-workflow-saved-da.R`, 11 cases). D11 shown failing first. Probes P1 to P4. The request half reads no matrix ([w1/](w1/)). Findings F112-01 to F112-04 |
| W2 | `8085d7d0`, `9e2ee21d`, `d2928858`, `6f41a5dd` | The split, the handler, the registry entry, the adaptor. F112-03 fixed. F112-05 and the stale KIT-019 comment recorded for KIT-109 |
| W3 | `47e779eb`, `f6cdf762` | The run plugin dispatches on a managed session; F112-02's two guards. `test-omic-da-module.R` 1339/0 |
| W4 | `8a97ee22`, `a53f9060` | The bindings block; merge checks ([merge-checks.txt](merge-checks.txt)); F112-06 recorded for KIT-109 |
| W5 | `5ecc9fd1` onward | The proof and the measurement, below |

## Proof

| Check | Result |
|---|---|
| Merge checks | Protected paths unchanged (`NAMESPACE`, `man/`, `DESCRIPTION`, `inst/`, `tools/`); `Collate` equals `R/*.R`; `diaContrastOperation`'s body is unchanged |
| Simulation Setup ([proof/chain-a.md](proof/chain-a.md)) | `prot_dia` and `prot_dia_limpa` identical to their oracles. The runner, request, compute and accept halves were each called once; the worker handler 0 times (in-process session) |
| Run 1, targeted ([targeted.txt](targeted.txt)) | 54 files: 16,994 passed, 0 failed, 0 errors |
| Run 3, entry suite ([entry-suite.md](entry-suite.md)) | All 10 lanes identical, 0 step differences, no re-freeze: trx 29, lipid 35, metab 105, proteomics 193 |
| Run 4, DIA session and memory set ([dia-memory-set.txt](dia-memory-set.txt)) | 733/2, 191/1, 38/1, 174/4, 80/11, each equal |
| Run 5, `R CMD check` ([check.log](check.log)) | 3 WARNINGs, 3 NOTEs; the R-code NOTE is identical to KIT-111's |
| Run 7, journeys ([proof/chain-b.md](proof/chain-b.md)) | 112l, m, t, a, b, c accepted and exact. 112l/m/t add only `DA-contrast-output` 1/1; 112a/b/c dispatch as 111a/b/c |
| M2, lifecycle | Completion as above. 112f fails at `da` as declared: `DA-contrast-output` 1/0, no `da_publication`, head `correlation_filtered`, `da_results` NULL |
| M3, large pair ([memory-outcome.md](memory-outcome.md)) | Null check 111p against 111q: 0 different. 111q against 112q: 16 states, 105 comparisons with the `da_results` manifest, 0 different. Same 18 verifier mismatches. "Better" met; every "no worse" rule met |
| Run 8, e2e ([e2e-final.txt](e2e-final.txt)) | lipid 88, metab 190, trx 17, DIA 42, DIA_LIMPA 80, enrichment 98, TMT/LFQ 153: 0 failures, 0 skips |
| Run 9, skips ([skip-compare.md](skip-compare.md)) | `load_all` has no change. Check mode has one: `test-workflow-saved-da.R`'s 8 new KIT-112 blocks error, all on the known source-tree read of the checkpoint fixture |

## Budget use

- **GUI proof:**
  - the correction pass was not used;
  - one listing rerun, `test-workflow-saved-da.R` in check mode, recorded each block's message (run 9).
- **Memory:**
  - the correction pass was not used;
  - no instrumentation rerun was needed.
- **Host:**
  - every member ran under `systemd-inhibit` with no other `multischolar*` unit live ([power-state.txt](power-state.txt));
  - no agents ran during the measured windows (KIT-111's lesson).

## Hand-offs

- **To KIT-109:**
  - **F112-01:** recorded with P4's result ([findings.md](findings.md)).
  - **F112-04:** validation accepts a formula term on an all-NA design column, and the kernel then refuses it.
  - **F112-05:** the managed path leaves `da_publication` as the controller's pending step provenance (since KIT-019/025).
  - **F112-06:** the runner's own pending retry is not guarded.
  - The stale KIT-019 comment at `utils_workflow_operations.R:683-685`.
- **To KIT-113:** reuse the request/compute/accept split, `named_input`, the click guards and `compare-stores.R`. Check reachability first, as here.
- **Parked for the user:**
  - **P-37** (information): proteomics' DA tab reads its saved tables in the application.
  - **P-39** (information): the two click guards (F112-02).
  - P-36 and P-38 were not raised.
