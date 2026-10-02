# DIA policy handoff findings — 2026-09-18

Mode: audit and documentation repair, authorized by the user-approved amendment
plan. Published before repairing audited documents. This register is immutable;
closure is recorded separately in policy-remediation-20260918.json.

Scope: the existing DIA plan/index, DIA-009/013–016/018–024 and consumer/gate
handoffs. No application, dependency, scientific, runner or benchmark changes.
The dirty worktree predates this amendment and is preserved. Exact predecessor
contents/hashes are in [policy-predecessor-20260918.json](policy-predecessor-20260918.json).
Budget: document/schema/link/graph/hash checks only; no empirical runs.

| ID / severity | Observed evidence | Claim and impact | Smallest repair | Closure evidence / blocked state |
| --- | --- | --- | --- | --- |
| DIA-NEXT-01 / high | DIA-015 revision 1.2 and consumer-results.md omit the later qs2 and reader outcomes; its required history scope is standalone-rds-history-v1 | Stale handoff can select the wrong treatment or silently inherit an incompatible gate | Refresh the current evidence summary; add a distinct initially held qs2 history scope, leaving old gates intact | Successor references, unchanged old gate records and no promotion; decision |
| DIA-NEXT-02 / high | R/utils_workflow_saved_table_codec.R:6 requires plain data.frame, includes all row names under 8 MiB; the fast reader/class/name bridge exists only in the profiling driver | Successful native I/O experiment is not an accepted package-store adapter | Bound a DIA-013 continuation around the existing store, restoration, read capabilities and action lifecycle | Pending continuation with unchecked tests, exact write scope and legacy compatibility; implementation |
| DIA-NEXT-03 / high | projection-decision.md:189 qualifies only five actions on the same 960000-row table versus whole-table RDS; one-off is held and smaller shapes rejected | Does not establish benefit over resident GUI state, changing QC generations or projected-only peak | Require actual artifact reuse and complete producer/consumer/output boundaries in DIA-015; keep small-table paths | Explicit comparators, natural-use mapping and claim limits; promotion |
| DIA-NEXT-04 / high | checkpoint-results.md:220 is one matched disk-RDS/qs2 pair, separate from the earlier resident build; required history matrix and seven pairs remain absent | History acceptance still blocks composition and GUI integration | DIA-009 qualifies selected qs2 versus matched resident ownership, sharing/aliases, forward/revert/reset and complete export | New held scope with unchanged numerical authority and missing-evidence requirements; execution |
| DIA-NEXT-05 / medium | measurement.md:2045 attributes about 81.1% sampled CPU of one Apply to audit; :2062 identifies the duplicate substantive digest | Reader improvements do not remove the measured active-operation cost | Add a separately gated DIA-019 subtask to reuse only that identical digest within one record | Named owner, exact audit contract, separate Apply comparison and no canonicalization rewrite; implementation |

Evidence identities checked in the preceding read-only assessment:

- Frozen protocol: 517edccb0906165698f3c1f165dfe7d961a17910d93363a60e843fa466c62379.
- Reader qualification: d8caba96bd256095ab6b1ea85e9717d0431777d5263d6ef11a74710e8c5b11ac.
- Reader screen: 662ec5e005bd1e54959ea31c850cf7dd91b8fdbb076a1cbbc742d9b290ac12a0.

Documentation closure will not establish implementation, scientific agreement,
performance acceptance or permission to integrate. Existing PQC contracts and
all historical failures remain unchanged.
