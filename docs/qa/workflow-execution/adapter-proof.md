# WFX-004 adapter proof

**Current decision:** the eight planned comparisons are complete. The 300-sample
rollup satisfies both paired integration non-regression limits; sample filtering
is scientifically correct and memory-compliant, with inconclusive timing.
WFX-004 remains open for full wide-annotation handoff qualification. No GUI default
or production activation changed.

Authority: [WFX-004](../../../tickets/openai/WFX-004.md), [frozen contract](contract.md),
[programme audit](../wfx-ticket-audit-20260930.md), and the user-authorized
[4 GiB startup amendment](startup-policy-20261002.md).
[Bound measurements and ratios](wfx-004-comparison-20261002.json).

## Comparison

Both arms use the same guarded execution, controller commit and required-output
publication boundary. Controls are explicit P14 computation and the unchanged
native S4 sample filter. P14 retains eight threads, three complete-sample batches,
ordered exact summation, insertion-order preservation and fastduckdb fast_r writing.

| Cell | Control seconds (rep 1 / 2) | Candidate seconds (rep 1 / 2) | Control peak GiB (rep 1 / 2) | Candidate peak GiB (rep 1 / 2) | Paired decision |
|---|---:|---:|---:|---:|---|
| Rollup: 24 million precursors, 300 samples | 118.060 / 90.613 | 105.277 / 99.881 | 8.151 / 8.103 | 8.184 / 8.302 | Integration non-regression |
| Filter: 10,000 rows, 4,000 samples | 19.019 / 18.508 | 18.489 / 28.084 | 2.468 / 2.432 | 2.412 / 2.410 | Inconclusive timing |

The rule is peak ratio <=1.05 and time ratio <=1.20 in **both** pairs.
Rollup ratios are 1.0040/0.8917 and 1.0246/1.1023 (peak/time).
Filter ratios are 0.9776/0.9721 and 0.9913/1.5174. One filter pair fails the time
limit; it was not retuned or rerun. Median rollup time is 102.579 s candidate versus
104.3365 s control; median peak is 8.243 versus 8.127 GiB. These are descriptive
medians, not an additional benefit claim or a comparison with the old resident R
baseline. No failed control supplies a finite denominator.

All eight valid members completed, passed scientific and required-output checks,
preserved original inputs, and had valid action-window measurements. Each primary
output matched all 1,200 independent native fingerprints and the native summary.
Each member produced its 12 required figures. Worker and monitor teardown completed.
All runs retained the 16 GiB application cap, 15.5 GiB early stop, zero swap, 2 GiB
continuous host reserve and eight-thread budget. Only startup headroom was amended
by the user from 16 to 4 GiB, before valid timing began.

## Implementation and verification

The compute-only adapters return existing references, summaries and dependencies.
Only the leased controller finishes scientific state and publishes outputs.
Sample identity counts include zero/NA quantities; inclusions, sorted ledgers,
source order, native container types and design cleanup are preserved. Admission
uses observed cardinality, width and complete-group estimates against the actual
worker cap; the kernel guard remains authoritative. No arbitrary-function plugin
or silent full-assay candidate fallback was added.

Retained evidence under `MultiScholaR-dia-proof-20260911-VhC1n5/wfx-004/`:

- `xvjBRz`: 54 small scientific assertions.
- `adapter-lifecycle-1`: 22 managed lifecycle/publication assertions, including
  output-only retry and duplicate completion.
- `wftnMg`: 27 domain/admission assertions and 1,079 affected native/saved/reference-
  state assertions, all passing. Reader boundaries 1/2/3/7 and incorrect counting
  controls discriminate the declared science; oversized groups refuse.
- `7jRwxZ`: actual historical checkpoint selection passes with assay materialization
  and reads forbidden. A preceding assertion-only failure `zbsGlW` is retained.
- `aZdseq`: predecessor readback passed on owned copies of both new adapter outputs,
  including history, matching bounded table reads and pending-publication recognition.
  Original accepted heads remained unchanged; no science was recomputed.

Final timing libraries: candidate `eZAHKm`, control `d0H056`.
`comparison-plan-reference-selection.rds/json` binds the primary harness successor;
completed filter members retain their original `comparison-plan-host4` bindings.
The filter path already selected its target state and did not execute the corrected
historical-selection branch. [Source predecessor](wfx-004-predecessor.json) supports
scoped source rollback; never reset the worktree or remove accepted artifacts.

## Preserved failures and scope limits

Two startup refusals preceded the user's 4 GiB amendment. The first primary setup
called public Revert and hydrated the historical S4 outside the action window.
It was stopped and preserved as `primary-control-1`; its 12.48 GiB service peak is
**not** a P14 operation measurement. A worker appeared after the initial observation
and before stop completed. The user explicitly authorized selector validation and
a replacement control. The verified selector changes only the owned replay head by
references; it does not alter public Revert or GUI behavior.

Initial packaging lacked isolated dependencies; copies were subsequently hash-
verified. Scientific preflights exposed a missing fixture q-value role and row-name
comparison/ledger issues; their failures and bounded corrections remain retained.
The synthetic filter wrapper adds Q.Value=0 and declares Quantity roles without
changing frozen identity/counting inputs. Its raw inclusion arguments are preserved
under a finite 256 KiB operation-specific parameter bound; process requests remain
capped at 1 MiB.

Charged cgroup peak and complete-operation elapsed time are measured. Worker CPU,
anonymous/file observations, staged bytes and raw I/O counters are retained in the
machine-readable report. Layered-device I/O totals can count both mapper and
underlying devices and must not be described as unique physical disk bytes.
Spill is disabled; attempted-spill counts are not independently instrumented and
remain explicit NA rather than invented zeros. Costs not separately attributable
remain within the complete-operation measurement.

**Remaining WFX-004 gate:** historical wide numerical evidence does not establish
that the new worker handoff supports its larger summaries. The worker currently
returns filtering counts inline under an 8 MiB message bound. The retained wide
fixture has 10,000 distinct approximately 2 KiB accession labels, so its per-protein
summary deserves explicit reference-backed handoff verification before claiming
that whole domain supported. No new full-wide success is claimed here.

Production remains held for the private/unsigned runtime and missing irregular
real-data corpus. WFX-005 owns the installed GUI journey and downstream/Revert/
session checks; WFX-006 owns automatic selection. These successful primary
comparisons do not authorize either default activation or all-workflow rollout.
