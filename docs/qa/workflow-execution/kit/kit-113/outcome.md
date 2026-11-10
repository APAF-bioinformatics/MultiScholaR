# KIT-113 outcome: proteomics' enrichment on the worker (2026-11-10)

**Done.** On a managed session, proteomics' enrichment computation now runs on the worker's `omics` role, one action per Run, dispatched by the shared enrichment frame through the proteomics runner. That covers the DA tables' read, the annotation match, the organism filter's application, the kernel and the product. The controller keeps:
- the annotation lookup, the mapping and the method;
- the display and the state commit with its provenance;
- the external-resource records.

Without a managed application, everything runs in process exactly as before.

**Results and records unchanged.**
- **In process:** every lane, capture and count is identical, with nothing re-frozen.
- **Worker against in process:** KIT-107's eight enrichment cases give identical states, products, provenance, records, files and bus values, at core and module level. The only differences are the ticket's recorded managed-path exclusions, each asserted exactly.

**No memory claim, and a projected cost.**
- No journey reaches enrichment, so no managed session ran one (P-40).
- D8's W0 probe found "better" unreachable. It projects the application **about 0.66 GiB higher** at the enrichment step on the managed path (1.73 against 1.08 GiB). The reason: the results carry interactive plots, and the controller must hydrate the worker's copy (535 MiB), where in process the live object shares memory with the loaded packages.
- As D8 declared, the move was made and judged on the "no worse" rules, which all three journeys meet. **P-43** asks the user to keep or revert it; the lead recommends reverting ([memory-outcome.md](memory-outcome.md)).

- **How it is built** (KIT-111's and KIT-112's pattern):
  - `runWorkflowEnrichment` and `runDIAReferenceEnrichment` keep their names and signatures and call helpers in today's order;
  - the managed halves: `requestWorkflowEnrichment` (reads no DA table), `computeWorkflowEnrichment` (worker) and `acceptWorkflowEnrichment`;
  - the handler `computeWorkflowDiaEnrichmentAction` (`commit = NULL`);
  - the registry entry moves to the `omics` loop with `named_input`, and the in-process loop is gone;
  - the frame's managed branch sits behind the F112-02-style guard, and the transfer values travel as `rds` state values.
- **Defect fixed, test first:** F113-02. The annotation match results would have made the worker's value exceed its 8 MiB bound above about 20,000 proteins; they now travel by reference. No result changed.
- **Build:** `kit_113` = `core/browser-smoke.cmBMvS` (`95600989`). `R/` +388/−54 lines in 4 files ([lines.tsv](lines.tsv)), all above the ticket's estimates:
  - `utils_workflow_enrichment.R` +235 (expected +70 to +140): roxygen for 15 internal helpers, no duplicated logic (F113-05);
  - `utils_workflow_dia_da.R` +38 (expected +25 to +50);
  - `mod_omic_enrich.R` +67 (expected +25 to +50): the managed branch inline, because the work item allows one function only;
  - `utils_workflow_operations.R` −6.

## Work items

| Item | Commit(s) | What |
|---|---|---|
| W0 | `2bad7519`, `132e5445` | Start condition met; M0; the line map (no stop); the reachability probe: "better" unreachable ([w0/](w0/)) |
| W1 | `9f9cbccc`, `a4b871e7`, `9a5bd23c`, `dcdc4a3f` | Guards on the start commit; D11 shown failing first; one new frozen capture block (two byte-identical dumps); probes P1 to P5, no difference; codec `rds`; F113-02 |
| W2 | `2984636c`, `c2db9053`, `b11be207`, `1b5d7f44` | The splits, the halves, the handler, the registry; F113-02 fixed test first; core equivalence on all eight cases |
| W3 | `752b2ad3`, `0b770216` | The frame dispatches through the runner on a managed session; the unmanaged frame byte-identical |
| W4 | `ae3861d4`, `49625d63` | The bindings block; merge checks ([merge-checks.txt](merge-checks.txt)) |
| W5 | `95600989` onward | The proof, below |

## Proof

| Check | Result |
|---|---|
| Merge checks | `Collate` 342 = `R/*.R`; NAMESPACE identical; protected paths unchanged; `R/` touches only the four files in scope; `diaEnrichmentOperation`'s body unchanged |
| Simulation Setup ([proof/chain-a.md](proof/chain-a.md)) | `prot_dia` identical to its oracle. The in-process runner, its helpers and the two `runDIAReferenceEnrichment` halves were each called once; the managed request half and the worker handler 0 times |
| Run 1, targeted ([targeted.txt](targeted.txt)) | 59 files: 19,924 passed, 0 failed, 0 errors |
| Run 3, entry suite ([entry-suite.md](entry-suite.md)) | All 10 lanes identical, 0 step differences, no re-freeze: trx 29, lipid 35, metab 105, proteomics 193 |
| Run 4, DIA session and memory set ([dia-memory-set.txt](dia-memory-set.txt)) | 733/2, 191/1, 38/1, 174/4, 80/11, each equal |
| Run 5, `R CMD check` ([check.log](check.log)) | 3 WARNINGs, 3 NOTEs; the R-code NOTE identical to KIT-112's |
| Run 7 / M1, journeys ([proof/chain-b.md](proof/chain-b.md)) | 113l, m, t accepted and exact; dispatches identical to 112l/m/t; "no worse" met |
| M0p, the probe ([w0/enrichment-heap.md](w0/enrichment-heap.md)) | "Better" unreachable: in process 1.075, managed controller 1.732 GiB |
| Run 8, e2e ([e2e-final.txt](e2e-final.txt)) | lipid 88, metab 190, trx 17, DIA 42, DIA_LIMPA 80, enrichment 98, TMT/LFQ 153: 0 failures, 0 skips |
| Run 9, skips ([skip-compare.md](skip-compare.md)) | RUN9 |

## Budget use

- **GUI proof:**
  - the correction pass was not used;
  - W3 used one in-item correction (a lazily evaluated accept half, caught by its own test before the commit).
- **Memory:**
  - one probe attempt; a launch that failed before R started (a missing `/usr/bin/time`) is recorded in the probe's report;
  - no correction pass.
- **Host:** every member ran under `systemd-inhibit` with no other `multischolar*` unit live ([power-state.txt](power-state.txt)). No agents ran during the measured windows.

## Hand-offs

- **To KIT-109:** the convention divergences recorded in W2 and W3 (parameter counts, the inline managed branch). The findings are otherwise closed here ([findings.md](findings.md)).
- **Parked for the user:**
  - **P-40:** no managed session ran an enrichment; accept the probe, or a harness ticket.
  - **P-41:** information, what the worker does not inherit.
  - **P-42:** information, the worker writes the pathway files directly.
  - **P-43:** keep or revert the placement; the lead recommends reverting.
  - **F113-07**, with P-39.
