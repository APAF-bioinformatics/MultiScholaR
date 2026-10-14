# KIT-022: kit debt before closure

Status 2026-10-14: STATUS_LINE. Authority: `tickets/kit/00-KIT-PLAN.md` (decision row 2026-10-14; run ledger KIT-022: none, one install recorded); protocol v1.3. Control `kit_018`; this build is KIT-017's candidate.

## The four items

1. **DA source reference** (`R/utils_workflow_dia_da.R`, `runDIAReferenceDaAnalysis`). The DA UI arguments are committed onto the current state before the contrast loop; the manifest's `source` is then the head. Identity check on the KIT-008 screen project ([da-identity-check.R](da-identity-check.R); `before` from a clean worktree at `befc0e3f`, `after` from this ticket): before, the manifest source was not the head and the module holder kept the pre-commit reference; after, the source and the holder are the head, a second run with the same settings hits every contrast product (same entity ids), the contrast identity computed on the pre-commit source and on the head in one process is equal (the DA defaults the commit writes, `differentialAbundanceAnalysis` and its helper, are now stripped from the identity like the UI block: they duplicate the operation's explicit parameters), and the result table is identical by content (4,790 rows, same digest) with identical products ([da-identity-before.json](da-identity-before.json), [da-identity-after.json](da-identity-after.json)). Cross-run identity hashes are not comparable under `load_all` (the software identity falls back to a function body with source references), which the record notes. The KIT-019 worker handler compares the manifest's source with the head by id again.
2. **The active owner from the state scope** (`R/utils_workflow_dia_state.R`). `getDIAActiveManager()` resolves the innermost open `withDIAStateScope` owner first and the Shiny session's owner second; the scope pushes its owner on entry and pops it on exit (errors included). The kit-level readers keep their signatures (the publication kernel calls `readDIAContrastOutput` itself, so a threaded argument could not have reached it without touching the kernel: ticket revision 1.1). The KIT-019 worker runs its publication inside the scope with no mock session. Test: the scope resolves without any session, nests and restores, pops on error; a reader pages a result inside the scope with no session; the session still serves the modules outside a scope.
3. **Q-value admission** (`computeDIAReferenceQvalue`; KIT-014 post-hoc audit PH-02). The family's admission estimate before any work: one streamed page (10,000 rows of the retained columns, or fewer within 64 MiB) times the row bytes plus 2,048 declared working bytes per row, plus the evidence lookup's 64 MiB bound; refused with `workflow_execution_unsupported` above the cap, recorded in `execution$admission` (policy `qvalue-streamed-page-v1`). Test: the estimate is recorded on the fixture and a cap just above the overhead refuses.
4. **Legacy-global adapters recorded** (`docs/qa/workflow-execution/kit/README.md`, "Kernel contracts the kit adapts to"): the publication phases' local global bindings for the legacy figure writers, the kernels' call-stack parameter resolution, and the active-owner resolution, named as adapters for the KIT-017 checklist. No code change.

## Tests (working tree, proof library on `R_LIBS`)

| Suite | Result |
|---|---|
| test-workflow-dia-products.R | 100 pass (new: the active owner from the scope: no session, nesting, pop on error, a reader inside the scope; the KIT-019 worker tests unchanged) |
| test-workflow-partitions.R | 299 pass (new: the Q-value admission recorded; a cap just above the overhead refuses) |
| test-module-ci-prot-da.R / test-prot-07b-da-handlers-compat.R | 245 / 55 pass |
| test-workflow-execution.R | 64 pass (10 skips) |
| test-workflow-native-session.R | 500 pass, 1 pre-existing failure+error (the process-loss test under load_all; baseline) |
| test-workflow-checkpoint-dia-state.R | 236 pass |
| KIT-019 probe rerun on this code (worker without the mock session) | reproduces ([da-publication-probe-kit022.json](da-publication-probe-kit022.json)): Glimma absent from the application-side process, loaded in the worker, +0.052 GiB over the DA step, the page asset shown |

## Claim limits

- Item 1 changes the order of two existing commits; the contrast identity strips the committed arguments, and the identity check shows the result tables and products unchanged. The DA state after the step carries the same arguments as before.
- Item 3's estimate is an admission bound, not a measurement; the cgroup guard remains authoritative (as for the family).
- No counted member; KIT-017's journey pair measures the stage with this build as the candidate.
