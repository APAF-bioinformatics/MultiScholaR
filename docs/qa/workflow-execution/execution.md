# WFX-003 guarded execution boundary

Current decision: WFX-003 is verified; [final closeout](wfx-003-closeout-20261001.md)
records 192 installed lifecycle assertions, 1,135 affected regression assertions
and predecessor readback. The [audit](wfx-003-audit-20261001.md) remains immutable;
its [four repairs](wfx-003-remediation-20261001.md) and the subsequent recovery
harness correction are complete. No GUI scientific route or default has changed.
[WFX-004](../../../tickets/openai/WFX-004.md) owns the scientific adapters.
The [frozen contract](contract.md) and [programme audit](../wfx-ticket-audit-20260930.md)
remain authoritative.

## Runtime and ownership

The internal `runWorkflowAction(operation, input_references, parameters,
resource_policy)` receives a leased controller as its server-owned resource policy.
It returns a request ID after monitor launch. `pollWorkflowAction()` advances
admission/completion; `resumeWorkflowAction()` obtains the existing bounded
reference/summary/dependency result or retries its required publication. It never
recomputes accepted science. The entry point is finite: WFX-003 registers only
`wfx_fixture` in diagnostic mode. Unregistered production operations refuse until
WFX-004 supplies their audited adapters. Arbitrary function/command dispatch and
full-assay fallback are not supported.

Each managed application uses a systemd user slice:

```text
Application A: aggregate kernel memory limit, zero swap
  +-- controller service C: Shiny R process and project mutation leases
  +-- monitor service M: current/anonymous/file memory and host headroom
  +-- worker service W: independent hard cap; includes all descendants
```

The service manager binds W to both C and M. Loss of either stops W independently
of the GUI. `OOMPolicy=kill` sets the worker's kernel `memory.oom.group=1`. All
numeric libraries use one internal thread; the server policy retains an eight-thread
compute budget for the qualified adapter. An OS admission lock admits one heavy
action per application, including across controller processes.

A has 16 GiB hard memory, zero swap, 15.5 GiB early stop, 4 GiB startup headroom
and 2 GiB continuous host headroom. The startup value follows the
[user amendment of 2026-10-02](startup-policy-20261002.md); historical receipts
retain their original policies. Before W launches, reserve observed non-worker
residency plus 1 GiB growth; W receives the remainder of the aggregate envelope.
The monitor targets 50 ms observations. A gap above 100 ms invalidates performance
measurement. Kernel hard caps remain the backstop; unqualified GUI/ancestor OOM is
not a successful qualification. Browser/client and operator-launcher memory are
outside the server-side A measurement.

Lifecycle tests use the explicitly smaller diagnostic policy: 3 GiB A, 2.9 GiB
early stop, 128 MiB growth reserve and 2 GiB host headroom. These synthetic fixtures
are not scientific capacity claims. Replacement-controller crash cases have their
own independent application slice; their caps were not raised to accommodate an
extra test controller in the same application.

Project mutation uses filelock's OS advisory lease in the existing store.
Coordination metadata lives alongside that store, while execution journals reuse
`workflowExportPaths()` and `writeWorkflowExportJournal()`. No competing scientific
store, checkpoint codec or cache was added. Cooperating current-version writers
obey the lease; stop managed controllers before rolling back to an older binary
that predates that lease contract.

Dispatch freezes request ID, project incarnation, input checkpoint IDs, expected
head, controller generation, parameter hash and installed build identity. Workers
receive only those references and bounded parameters. Stateful head/commit entry
points refuse worker callers. The controller validates the token under its retained
lease and rechecks the durable head at publication. Revert/cancellation advances
generation before acknowledgement. Closed-controller or superseded-generation
completions cannot deliver a GUI update; an older completion cannot cancel a newer
active job.

Before head advancement, the accepted-result descriptor and proposed head are
durable in the request journal. Recovery checks whether the old or proposed head
is current, finishes that one recorded transition, and refuses conflicting later
heads. Duplicate messages return the recorded result, not another checkpoint.
Required publication retains its existing independent pending/completed journal.

## Measurement and qualification

The action timer begins before request validation/admission. An FD-local
`memory.peak` window opens before those allocations, is read once after required
publication/teardown, and then closes. Retained GUI residency is included; parent
and child peaks are not added together. Receipts separately retain whole-service
peak, sampled anonymous/file maxima, worker CPU/I/O properties, application
accounting deltas, staged bytes and cleanup outcomes. Missing monitor evidence or
slow observations do not become invented zero counters.

A failed/publication-pending attempt is not a completed performance comparison.
Output-only recovery does not relabel its incomplete timing as a fresh successful
operation. Later verifier/oracle work is outside the closed action window and can
still affect the separately reported service peak.

`evaluateWorkflowQualification()` is a pure decision function for verified,
source-bound comparison records. It requires typed flags, two complete pairs,
independently frozen operation/input/parameter/runtime/boundary/policy/output
bindings, and finite positive counters. It rejects open runtime holds, stale
bindings, malformed/missing evidence and absolute-limit violations. All 59 frozen
numeric cases pass. Publication/science failures never supply capacity denominators.
The existing WFX-002 private-API/unsigned-runtime and real-corpus release holds
remain; a numeric decision is not permission to activate GUI defaults.

## Installed evidence

Evidence root: existing `MultiScholaR-dia-proof-20260911-VhC1n5/wfx-003/`.
The audited predecessor library is `browser-smoke.ZB8MQ2/library`; the repaired
library is `browser-smoke.iR0a6x/library` (192 installed assertions passed).
The new R file is registered in DESCRIPTION Collate. Filelock 1.0.3 is installed
and declared; parallel is an R-recommended dependency. The repository has no
renv.lock; no project-wide dependency migration was introduced.

| Suite | Evidence and result |
|---|---|
| `protocol-1` | Preserved failure: namespace assertion spelling, unsupported MemoryOOMGroup property, and nested fault-controller admission refusal |
| `protocol-2` | Corrected installed protocol checks passed, including both controller-death windows |
| `resources-1` | Core runtime/resource cases passed; regression wrapper lacked package context and exposed inherited standalone-source assumptions |
| `resources-2` | 151 assertions passed on final installed code, including 59 frozen numeric cases, real allocation descendant membership, Shiny session/HTTP liveness, publication-only retry, stale/duplicate/cancel paths, monitor loss and controller recovery |

The last suite ran affected regressions against the installed package:
resources 293, identity 73, projection 406 and exports 255 assertions passed.
Standalone recovery tests retained 81 passes, 17 failed assertions and one error
in their two existing source-loaded subprocess cases. Their hash function is
unchanged from the predecessor, and evaluating that predecessor in the same
standalone environment reproduces missing `workflowTableFileIdentity`. The raw
failures remain in `affected-regressions.rds`; they are not counted as passes or
silently repaired. New controller crash/recovery tests use the installed namespace.

The correction preserved the 3 GiB diagnostic cap and all runtime ceilings.
The unsupported systemd property was replaced with documented OOMPolicy, and
replacement-controller fixtures were separated into independent applications.
The final correction also fixed regression package context and completion delivery
from closed/invalidated controllers. Both affected rechecks are consumed.

## Operator entry points

Legacy `run_app()` formals and ordinary execution remain supported. Managed GUI
launch is explicit, using the installed package:

```r
MultiScholaR:::launchWorkflowApplication(existing_dedicated_directory,
    mode = "gui", port = 3838L)
```

The managed path reuses run_app/app_server, attaches application-level ownership
to sessions and closes project execution controllers when the last session ends.
It does not migrate existing GUI operation handlers. To stop an owned application,
stop its returned controller/slice units through the user service manager.

The existing profiling entry point now supports the two installed suites:

```sh
bash tools/profiling/run_workflow_checkpoint_proof.sh --workflow-lifecycle \
  "$installed_library" "$fresh_evidence_directory" "$repository" protocol
bash tools/profiling/run_workflow_checkpoint_proof.sh --workflow-lifecycle \
  "$installed_library" "$fresh_evidence_directory" "$repository" resources
```

These commands do not authorize retries. Source test inputs are hash-bound before
launch; the package implementation is loaded from its installed namespace, never
injected from helpers. The resources suite repeats the protocol checks on the same
final build and runs the affected regressions afterward.

Rollback restores only this ticket's source changes from `wfx-003/predecessor.json`
and the WFX-002 repaired installed library; accepted files and journals are retained.
The explicitly authorized five-minute predecessor readback passed in
`wfx-003/browser-smoke.FqyIDO`: the old installed library read matching table data
and history and recognized pending publication. Original input head bytes were
unchanged. The four audit findings subsequently closed in the repair suite. The standalone recovery regression gate subsequently passed in the
[final closeout](wfx-003-closeout-20261001.md).
