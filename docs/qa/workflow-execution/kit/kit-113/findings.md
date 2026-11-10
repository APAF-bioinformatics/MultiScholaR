# KIT-113 findings and lead decisions

## Lead decision at W0 (2026-11-10)

**The line map** (`w0/lines-at-start.txt`) does not stop W0. Every named function exists. Every named body is byte-identical to `525cbd78` apart from moving. Only KIT-112's anticipated changes touch the referenced lines.

**D8's reachability probe** (`w0/enrichment-heap.md`; one counted attempt, 29 min, inside the 30-min ceiling) gives **"better" UNREACHABLE at W0.**
- **The figures:** on the decision case (the live kernel over the offline doubles, 16,000 × 300 × 3), the in-process figure is 1.075 GiB and the managed controller's is 1.732 GiB. The difference is −0.657 GiB against the +0.25 GiB rule. The cold pass gives −0.947 GiB.
- **The move is projected to make the application worse, not just "not better".** The kernel's results carry `ggplotly` widgets (about 31 MiB each). Serialized portable results are 187 MiB for 3 contrasts. Hydrated in the controller, they are 535 MiB, against 16 MiB for the live object, which shares its environments with the loaded packages. The state commit then adds about +1.2 GiB on top, on both paths. The figure is not data-scaled: it is the same at 400, 4,000 and 16,000 proteins.
- **The deterministic double stores no plots.** With it, the move would lower the peak by 1.03 to 1.08 GiB at 16,000 proteins.

**Decision: the frozen protocol is followed.** D8's declared outcome for "unreachable at W0" is that the move is made, as the user's placement decision ("Join + move backend to worker", KIT-111 to KIT-113). It is then judged on the "no worse" rules only, and contingent **P-43** is raised beside P-40. P-43's wording anticipated this cause ("the results carry the plots' captured frames").
- The outcome will say plainly that the probe projects a higher application peak on the managed path, about +0.66 GiB at the enrichment step.
- **The journeys cannot see it,** because none reaches enrichment.
- P-43 recommends the user's choice between keeping the move and reverting it, with this evidence.
- **Reverting is cheap:** KIT-113's commits touch only the files in its `files_to_modify`.

**Pre-run declaration by the W0 agent** (D8 is silent on it): the warm pass decides. Both passes decide the same way.

## Lead decisions at W1 (2026-11-10)

- **Merged** fast-forward (`9f9cbccc`, `a4b871e7`, `9a5bd23c`). P1, P2 and P5 show no difference that D3 cannot absorb, so no stop. P3 is confirmed by reading, and P-41 stands.
- **The transfer codec is `rds`** (D3, "the transfer values"), from P1. The display-table route returns NULL for a zero-row annotation table. The `full` corpus's table is empty, so that case is real.
- **The coverage check added one frozen block,** the D6 fields W1 (b) lacked. It was frozen from two byte-identical dumps in separate processes.
- **The `.html` plot pages** differ between two in-process Runs too (generated widget ids). D6 compares names, and sha256 for TSVs only, so this is not a difference the move introduces.
- **F113-02 is accepted into W2's design.**

## Lead decisions at W2 (2026-11-10)

- **Merged** fast-forward (`2984636c`, `c2db9053`, `b11be207`). No stop rule fired. Equal to the start:
  - the in-process guards;
  - the five proteomics lanes (193/0/0, no re-freeze);
  - the DIA set (733/2, 191/1, 38/1, 174/4, 80/11);
  - the other omics' worker tests (289/0).
- **D11's blocks fail first and now pass** (registry 1044/1044). **F113-02 was fixed test first.** The core equivalence holds on all eight W1 (b) cases, each difference being one of D6's recorded exclusions, asserted exactly.
- **F113-03 is information.** It corrects F113-02's mechanism; the fix is unchanged.
- **F113-04 is accepted** as one more of D6's recorded exclusions. It is the same kind as D6's worker-failure case: a notice-identical refusal on the managed path whose bus metadata write is the accept half's, so it does not happen.
- **F113-05 (lines) is recorded.** The excess is roxygen for 15 internal functions, and no logic is duplicated.
- **The convention divergences are recorded:** two helpers carry 11 and 9 parameters, keeping today's statements.

## Lead decisions at W3 and W4 (2026-11-10)

- **W3 was cherry-picked** onto the generator commit (`752b2ad3`, `0b770216`), as KIT-112's W3 was.
  - `mod_omic_enrich.R` +67 lines, insertions only. The unmanaged frame is byte-identical; the stop rule did not fire.
  - The new blocks fail first on the start frame (101 failed, 2 errors) and pass after (2296/2296).
  - Every targeted count equals the start; the five lanes are 193/0/0, with no re-freeze.
- **F113-06 is accepted** as one more of D6's recorded exclusions, the same kind as the missing-table divergence. On the managed path an R7 refusal has already written `uniprot_dat_cln` and shown the no-mapping warning; nothing is saved or dispatched. In process there is no R7 refusal at all.
- **F113-07 is information,** carried with P-39 (the click guards, F112-02). The pending-publication guard refuses enrichment while the shared proteomics runner holds any pending publication.
- **W4's merge checks pass** (`merge-checks.txt`): `Collate` 342 = `R/*.R`; NAMESPACE identical; the protected paths unchanged; `R/` touches only the four files in scope; `diaEnrichmentOperation`'s body is unchanged.
- **The targeted set** is KIT-112's 54 files plus five from the dependents grep (59).

## Findings

| Id | Where | What | Decision |
|---|---|---|---|
| F113-01 | — | (used and dropped at the ticket's review: the proxy variables; no difference by the code) | — |
| F113-02 | the handler's value (D3, R3) | **The value carries the annotation match results, about 416 B per DA protein.** At 16,000 proteins they are 6.66 MB of the 6.76 MB value, against the 8 MiB acceptance limit (`utils_workflow_execution.R:1771-1774`). By linear scaling, a Run on about 19,900 or more proteins would be refused at acceptance (`invalid_result`) after the worker has computed, although the same Run completes in process today. Found by P4. | **Fixed in W2, inside D3, test first:** the match results are saved as an `rds` state value, as the results are, and the value sends the reference. A value over the limit is shown refused before the fix and accepted after, with the returned list identical. No methodology change. |
| F113-03 | `writeWorkflowExecutionMessage` (`utils_workflow_execution.R:1300-1302`) | F113-02's refusal happens earlier than stated: the worker's own handoff refuses a value over 8 MiB ("execution handoff must be bounded immutable metadata"), so the controller sees a worker failure, not an `invalid_result` acceptance. | Information. F113-02's fix covers it. |
| F113-04 | the managed request half, the taxon refusal | On the managed path a Run refused at the method (an invalid taxon) does not write the bus's `enrichment_organism_filter` metadata. In process that write comes before the refusal. The notice is today's. | Accepted as one more of D6's recorded exclusions, the same kind as D6's worker-failure case. |
| F113-05 | lines | `utils_workflow_enrichment.R` 1,134 → 1,369 (+235, against the expected +70 to +140): roxygen for 15 internal helpers and halves, no duplicated logic. `utils_workflow_dia_da.R` +38; `utils_workflow_operations.R` −6. | Recorded in the outcome. |
| F113-06 | the managed request half (W2 order), R7 | The R7 identity check comes after the annotation resolve and the mapping. So a refused Run on the managed path has already written `uniprot_dat_cln` to the bus (when absent) and shown the no-mapping warning (filter on, no mapping). Nothing is saved or dispatched. In process there is no R7 refusal: a stale manifest is run and committed. | Accepted as a D6 recorded exclusion. |
| F113-07 | the frame's guard (D3) | The pending-publication guard refuses enrichment while the shared proteomics runner holds any pending publication, including another proteomics tab's, which only that tab's button retries. The notice says "already running". | Information, carried with P-39 for the user. |
