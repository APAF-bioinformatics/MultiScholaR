# WFX execution and qualification contract

Version 2, owned by [WFX-001](../../../tickets/openai/WFX-001.md).
The [qualification manifest](../../../tests/testdata/workflow-execution/qualification.json)
is the machine-readable successor binding to the
[programme v1.1](../../../tickets/openai/02-WORKFLOW-EXECUTION-PLAN.md).
The [original audit](../wfx-ticket-audit-20260930.md) remains immutable;
[remediation](../wfx-ticket-remediation-20260930.md) describes the runtime tests
still required. This contract freezes acceptance inputs; it does not qualify a
runtime or enable production selection.

## Authorities and applicability

The manifest binds XsaPvT application, fastduckdb and exact-sum binaries by SHA256,
its retained source files, native scientific oracle sources, dependency receipts,
and accepted lxQvhv/0bVIHx outputs. Resolve paths only under the declared evidence
root. Verify identities before use; a directory name or latest-state filename is
not an authority. Missing files or changed hashes hold the consuming check.
Canonical descriptor.rds and the existing serialized Parquet codec metadata own
R types, factors, attributes and row identity. A JSON descriptor is an index,
not a replacement for that typed descriptor.

Primary is the BpSyxB **protein_peptide_filtered** checkpoint: 24 million
precursors, 300 samples. Its provenance is synthetic OQADra (seed 902002),
through cuurUw and BpSyxB. Wide is synthetic seed 902005: 3.84 million precursors,
48 samples and 2,048-byte annotation strings. Select all complete samples in
recorded design order. Never substitute BpSyxB's later replicate_filtered head.
The manifest binds payload IDs/hashes and the heads used to locate their history.
Payload hashes are inherited from accepted receipts; WFX-001 reads only small
metadata and does not claim a fresh full-payload integrity scan.

The independent native receipts are 7X1Wpw/B (primary) and EYBmKx/B (wide).
P14's accepted output hashes and summaries match these controls. The 60v50g/P3
accepted reference is additionally retained. These historical receipts establish
reusable operation science, not new GUI timing, new release qualification or
universal capacity. Frozen A's memory failure is never a speedup denominator.
Original figures, checkpoints, failures and pending-publication records remain
owned by their existing stores.

**corpus_status = held.** No legally usable irregular real DIA case with frozen
source identity and permission has been established here. Both retained large
cases are labelled synthetic. Component implementation may proceed; WFX-005
production qualification and WFX-006 default activation cannot pass this hold.
Source-bound real data must fit WFX-005's existing small-check allocation, or a
budget amendment must precede execution. Do not quietly replace a capacity case.
The user's existing local experiment authorization permits reuse of generated
fixtures; it is not an assertion of external redistribution rights.

## Finite execution boundary

The proposed internal call is
`runWorkflowAction(operation, input_references, parameters, resource_policy)`.
It accepts only registered, qualified operation/domain/runtime/platform tuples.
Initially these are precursor rollup and Minimum Peptides per Sample. It returns
existing saved result references, summary and dependency structures through the
current writer/checkpoint/publication machinery. This is neither an arbitrary
function wrapper nor a new public R/S4 API or checkpoint format.

Workers receive immutable references and bounded parameters/metadata. They must
not capture an assay, S4 workflow, reactive environment or live connection.
Server-owned resource policy cannot be overridden by GUI arguments. Unsupported
schema, types, environment or oversized indivisible group refuses before
scientific execution; refusal is not successful migration. Legacy/small paths
remain available only inside their declared envelope. Never silently collect
an unsupported large query into R.

WFX-005 supplies the minimal rollup adoption bridge before WFX-006: reuse existing
references, or boundedly persist resident/imported and restored legacy state with
the existing codec, verifying restorability before releasing owners. Preserve
history. Refuse before copying when adoption cannot fit. Capability selection
must not flip the global reference_states flag to dispatch other methods.
Serialize/adopt within the measured action. Queue unqualified heavy actions while
a worker runs; test the next consumer's bounded compatibility path or explicit
refusal before full hydration. Broader import/design work belongs to WFX-007.

## Scientific adapters

Rollup preserves the native linear double-quantity domain, duplicate identity and
key validation, finite/nonnegative checks, frozen identification evidence,
distinct modified-sequence counts and native accession tokenization on bounded
distinct annotations. Partial missing groups sum observed values; all missing
groups retain missing results. The exact-sum extension consumes values in immutable
source ordinal order, including the `1e16 + 1 + 1` threshold case; tolerances and
native kernels remain unchanged. Global validation precedes batching.

Freeze P14: eight compute threads, three complete-sample batches, scoped
DuckDB 1.5.5 / duckplyr 1.2.1, insertion-order preservation **enabled**, and
fastduckdb 0.1.0.9001 `fast_r` writer. A batch may materialize an engine result;
10,000-row / 64 MiB R transfers do not prove upstream operators stream. Namespace
qualified duckplyr uses the scoped connection, stingy fallback policy and disabled
telemetry. Restore canonical sample/protein/peptide order and native codec output
metadata explicitly. Do not restart P15–P19 tuning or substitute the Arrow writer.
Required outputs are the accepted table, exact summary/dependencies, durable
checkpoint and the 12 required rollup figures represented in the bound receipts.

For the sample filter, C(s) is the distinct count of (Protein.Group,
Stripped.Sequence) over **all** valid identities for observed sample s, regardless
of zero or missing quantity. Retain iff C(s) >= k or s belongs to the trimmed,
deduplicated, nonblank/nonmissing inclusion list, with positive integral k.
An absent inclusion entry is reported but cannot create a sample. Reject missing
keys, duplicate/blank/missing design IDs and unmapped observed runs. Extra valid
design runs are permitted by the native validator and removed during native
`cleanDesignMatrixPeptide`; do not incorrectly require exact input/design sets.

Preserve source row order, sorted support/removal ledgers, inclusion reporting,
method argument metadata and native design cleanup in first-retained-data-appearance order,
even when the supplied design has a different ordering. The native join intentionally
converts factor Run to character; exact compatibility does not mean blindly
restoring the original factor. Other metadata follows the codec and native result.
The native S4 method and helper are independently source-bound, not rewritten as
the candidate's oracle.

`makeWorkflowExecutionFilterFixture(repeats=1)` has four heterogeneous samples:
A has two identities carried only by NA/zero quantities, B has one identity and
passes by inclusion, C has one and is removed, D has four and supplies skew.
Revision 2 keeps ten source rows and nine survivors. D now has four distinct
peptides on one protein, while A retains the same peptide on two proteins.
Protein-only and peptide-only count mutants therefore fail in opposite samples. Repeats=1000 is the
fixed 10,000-row, 4,000-sample WFX-004 comparison population, with 3,000 retained
samples. This is a contrasting semantic proof, not a large-study filter capacity
claim. Future tests vary reader boundaries 1/2/3/7 without changing scientific
membership or source ordinal. Dropping inclusions or counting only observed
positive quantities must fail the hand-authored expectations. The manifest also
specifies empty/all/none cases, invalid thresholds/keys/design, unused levels,
wide annotation and skew cases. Arbitrary floating-point reorder invariance is
not part of either contract.

## Controller, checkpoint and publication ownership

WFX-003 implements an OS-held exclusive project mutation lease released on process
death. Other controllers may read immutable artifacts but cannot mutate that
project. Revert, restore and legacy setters obey the same lease.

Dispatch freezes request ID, project incarnation, input checkpoint/entity IDs,
expected head, controller generation, parameter hash and execution identity.
Workers resolve these references, never latest. They create only request-owned
staged/immutable results; they cannot invoke commitDIAReferenceState,
publishDIAStateHead, finishDIAReferencePeptideStep or alter reactive state.
WFX-004 separates computation from existing stateful finish helpers.

The controller validates artifacts and the full token, rechecks the durable head
inside the lease critical section, then advances it. Cancel/Revert increments
generation before acknowledgement. The existing owned request journal records
request ID and accepted checkpoint. Duplicate delivery/restart returns that
outcome, never a second checkpoint. Reconcile head and journal before accepting
new mutations after a crash.

Required publication follows scientific commit using existing pending/acknowledged
states. Lost acknowledgement or failed output publication retries outputs only.
A valid stale artifact can remain request-owned but never advances the head;
cleanup never deletes accepted results. WFX-003/005 fault tests cover worker and
controller death, Revert/project replacement, duplicate completion, crashes on
both sides of head advancement and lost publication acknowledgement.

## Accounting and admission

Linux x86_64 is the initial supported platform; WFX-020/021 own other platforms.
Probe delegated cgroup-v2 accounting/limits and managed service lifecycle before
mutable work. An arbitrary existing R process is not automatically qualified.

Application cgroup A includes Shiny/controller, monitor and every worker subtree.
A has 16 GiB memory.max and zero swap. Browser/client is excluded and responsiveness
is measured separately. Worker W is independently killable. Reserve non-worker
residency plus 1 GiB growth and set W's cap to 16 GiB minus that reservation.
Admit only if the adapter's width/cardinality/skew/global-state estimate plus
native/R/writer headroom fits. Compressed bytes or row count alone cannot decide.

DuckDB's 12 GiB limit is an upper engine limit, not a reservation or whole-process
cap. Spill stays disabled. Preserve eight active compute threads with no nested
BLAS/pool multiplication. Allow one heavy action per application. Pending writer
thresholds are flush triggers, not memory guarantees; account for expanded strings
before conversion and refuse an oversized indivisible scientific group.

Start only with 16 GiB host headroom. The monitor outside W but inside A observes
at least every 100 ms, terminating W at 15.5 GiB combined usage, less than 2 GiB
host headroom or non-worker reservation overrun. Controller/monitor loss causes
fail-closed worker termination under managed lifetime/deadline. Qualified worker
overshoot must leave GUI usable; an ancestor/host OOM or unqualified GUI allocation
is infrastructure failure. Polling is not a hard limit and kernel caps can be
transiently exceeded; preserve events.

P is A's action-window kernel charged peak including anonymous and file cache,
without adding overlapping parent/child peaks or subtracting GUI residency.
Reliable action-window peak reset/read is mandatory for qualification. T spans
request snapshot/adoption through publication acknowledgement and teardown,
including admission/serialization/startup. Independent readback/oracle runs after
that window under a guard; report its service-wide cost separately.

Qualification is serial and uncontended. Both arms use metadata/footer preparation
only, no deliberate payload warming or eviction, recording natural cache state
and counterbalanced order. Missing counters or interruptions hold the measurement;
do not reconstruct time or grant free retries.

## Qualification, reuse and release

WFX-003 implements one pure evaluator against the manifest's hand-authored numeric
fixtures. Evaluate domain, evidence eligibility, science/lifecycle/integrity/outputs
and absolute guards before ratios. Invalid science/output/lifecycle is rejected;
missing runtime/data/bindings or invalid measurements is held. Neither can become
a capacity pass. All goals first require the complete frozen member set: two pairs, each with a
control and candidate. Missing/empty/incomplete pairs are held before any capacity
or ratio calculation; empty evidence cannot pass through vacuous truth.
Raw denominators and candidate counters must be positive finite
numbers; JSON null and explicit NaN/Inf sentinels are invalid.

For each of two paired repetitions, rP=Pc/Pb and rT=Tc/Tb. Benefit requires both
rP<=0.90 and rT<=1.20. Integration non-regression requires both rP<=1.05 and
rT<=1.20 plus bound predecessor benefit and absolute limits. One pair passing
and one failing is inconclusive; both failing is rejected. Medians are descriptive.
Capacity-only requires a valid control failing its predeclared memory/time guard
and exact independently verified candidate completion inside its limits; calculate
no ratio for that failed control. Other control failures cannot establish capacity.
Unknown goals and open runtime holds are held. Unsupported domain safely refuses.

The manifest fixes WFX-004's two cells and eight members; WFX-005's primary/wide
GUI cells and eight members; WFX-006's twelve explicit native/qualified-relational
crossover members and four supplementary
activation members. Use counterbalanced control/candidate then candidate/control
for repeated cells. Every future build/adapter/boundary/selector hash must be
filled by its owning ticket **before** timing. Unbuilt components are explicit
holds, not fabricated identifiers. Control and candidate use the same action
boundary, required outputs and absolute limits. Existing P14 timings cannot be
compared to an expanded GUI timer as though boundaries matched.

Calibrate routing from already bound historical measurements/analysis, then
freeze its hash before crossover validation. Insufficient calibration holds the
selector; validation outcomes cannot retune it. Hold out width/skew/boundary cases.
The final four WFX-006 members require matching WFX-005 two-pair core/runtime/
resource/output/boundary identities; only a separately tested router may differ.
Each supplementary pair must meet 1.05/1.20. Core changes require separate
requalification authority, never an automatic four-to-eight-member expansion.

A qualification receipt must bind schema/version, operation/input/parameters,
all component/runtime/build hashes, boundary/metrics, goal, members/ratios,
science/lifecycle/output checks, corpus applicability and open holds. Release
requires accepted WFX-002 runtime eligibility, WFX-004/005 qualification and
matching WFX-006 routing evidence with no holds. The current unsigned/private-API
runtime and missing real corpus cannot activate defaults. A checked ticket is
not a substitute for these receipts.

| Owner | Required proof |
|---|---|
| WFX-001 | Source bindings, finite contracts, typed metadata, small native fixtures and numeric acceptance data |
| WFX-002 | Reproducible native distribution, ABI/API/loading policy, supported platform and types |
| WFX-003 | Enforced topology, instrumentation completion/failure, pure evaluator and controller lifecycle |
| WFX-004 | Exact rollup/filter adapters, ledgers/metadata/mutation cases and paired action comparisons |
| WFX-005 | Installed GUI, three entry states, real corpus, downstream/Revert/session/publication recovery and repeated GUI comparisons |
| WFX-006 | Frozen selector, crossover/held-out tests, matching reuse bindings and operation-specific activation |

No performance members are authorized by WFX-001. Its budget is two five-minute
fixture/metadata checks, without assay hydration. The parent permits at most two
preflight correction passes and two affected small rechecks; no automatic full
retries. Preserve all failures. Later workflow outlines need their own bounded
plans before execution.

## Rollback and verification record

Rollback of WFX-001 removes only its two new deliverables and the appended fixture
helper, restoring the preceding helper header. It changes no application/default,
store, head, public method or publication format. Do not reset the branch or delete
unrelated work. The manifest binds the predecessor build/source/config; compare
those bindings and retained heads/descriptors after rollback. Runtime rollback
of future tickets must preserve accepted history and publication-pending states;
WFX-003/005 own the fault-injection proof of that behavior.

Historical revision 1 preparation completed 2026-09-30 (its overbroad completion
claim is superseded by the audit and repair below). Two bounded R checks (300-second
ceilings, both exited 0) read the accepted typed footers/descriptors/heads and ran
50 native-fixture assertions. Accepted output metadata and checkpoint/history
snapshots remain readable; removing the appended helper leaves the predecessor
fixture definitions loadable. No store or publication record was modified.
The accepted stores had no pending journal entries; this documentation-only
rollback does not claim a new pending-publication fault-injection pass. That
runtime proof remains assigned to WFX-003/005.

The static check verified 48 source/evidence bindings, all 44 authored numeric
outcomes and local links. Primary/wide P14 summaries and 1,200/96 fingerprints
match their independent native receipts. Original input head hashes and native
sources remain unchanged. No assay was hydrated, benchmark run, dependency
installed or GUI route enabled. Package startup emitted a harmless timedatectl
sandbox warning; both checks passed without substitution or correction.

Preparation scripts and receipts are preserved under the existing external root
at `wfx-001-20260930/`; their hashes are in the manifest's `verification.receipts`.
The commands were `timeout 300 Rscript --vanilla wfx001-metadata.R <repo>` and
`timeout 300 Rscript --vanilla wfx001-fixtures.R <repo>`, followed by
`python3 wfx001-validate.py` from the repository root. The latter is an independent
arithmetic sanity check of acceptance data, not the WFX-003 production evaluator.
Reader boundary reconstruction here proves fixture expectations only; actual
streaming adapter and mutation tests remain WFX-004's responsibility.

## Revision 2 audit repair

The [WFX-001 audit](wfx-001-audit-20260930.md) is immutable. The
[repair record](wfx-001-remediation-20260930.md) closes WFX001-A01–A05 against
successor evidence; the [predecessor snapshot](wfx-001-remediation-predecessor-20260930.json)
preserves exact audited bytes. Revision 1 receipts remain under
`historical_verification_v1` in the successor manifest and are not relabelled as
proof of revision 2.

The explicit crossover allocation is native / qualified relational / qualified
relational / native at each size. Automatic routing belongs to the separate
small/final checks. The same twelve crossover and four final members remain.
The authored decision set adds missing/empty/incomplete comparisons for all
three goals. `tools/profiling/validate_workflow_qualification.py` is the repository
fixture checker; the production evaluator remains WFX-003 work.

`tests/testthat/test-workflow-execution-fixtures.R` retains the native assertions
against the current package and adds counting-only mutants on the original data.
It uses no external evidence path or historical package. Run it with the normal
repository test environment after loading the current package and fixture helper.
The helper population remains 10,000 rows / 4,000 samples for repeats=1000;
its scientific discriminator changed under manifest version 2 before timing.

This ticket's rollback means restoring the byte-bound predecessor documentation
and fixture helper; it does not alter production readers, writers, checkpoints or
publication journals. The earlier checked claim of a demonstrated pending-output
rollback is withdrawn. Runtime recovery/rollback remains an unchecked obligation
in WFX-003/005. Removing only the entire helper would revert the initial WFX-001
addition; reverting this repair alone restores its saved version-1 bytes.

Production eligibility remains held for real-corpus and native-runtime
qualification. No new performance or release claim is made.
