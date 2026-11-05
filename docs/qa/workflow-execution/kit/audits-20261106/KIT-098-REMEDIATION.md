# KIT-098 remediation, 2026-11-06

Mode `audit_and_repair`, authorised by the user ("fix all recommendations"). This record closes the findings of [KIT-098-AUDIT.md](KIT-098-AUDIT.md), which is left unchanged.

The only file changed is `tickets/kit/KIT-098.md`, now at revision 1.1. No code, protocol, plan, index or other ticket was touched.

Closure here is at the ticket-text layer. Each repair becomes a requirement or acceptance criterion of KIT-098. It is verified when KIT-098 executes, with the closure evidence named in the audit.

## Changes for the index (coordinator to sync)

- **Dependencies:** `["KIT-097"]` becomes `["KIT-097", "KIT-092"]`.
- **Title:** "… protocol v1.7 …" becomes "The transcriptomics kit route: harness, in-process recorder, protocol v1.8, lifecycle and screen journeys exact".
- **Description:** rewritten for v1.8, the §7.1 oracle and the direct `ruv::RUVIII` oracle for technical-replicate RUVIII.
- **`ticket_revision`:** `1.0` becomes `1.1`.
- **`files_to_modify`:** `docs/qa/workflow-execution/kit/bindings.json` added.

**Validation.** `next-openai-ticket.sh KIT-098` reports one error, "frontmatter dependencies do not match index dependencies". It is expected and clears when the index is synced.

## Finding closure

| ID | Severity | Disposition | Repair in KIT-098 rev 1.1 (or owner) |
|---|---|---|---|
| F01 | blocking | Closed in ticket | The comparator `kit-098/notebook-compare.R` reads the journey member's store and maps each KIT-040 product to a step, state or product, with its projection. It uses `identical` unless a tolerance is declared first, and must be specified in `plan.md` before any run. Both journeys bind the notebook's design levels, contrast, RUV range, k cap and automatic best k. Two journeys are named: `098a` with `M` from group, compared with the notebook, and `098b` with `M` from technical replicates, whose RUVIII product is compared with KIT-093's direct `ruv::RUVIII` oracle (user decision 4). The RUVIII export is a journey step after DA that writes a product only. |
| F02 | blocking | Routed to KIT-092 (v1.7), with assertions kept here | Two harness defects go to KIT-092 (coordinator decision 5): a step with no bound oracle is skipped, and a non-DIA member without oracles is treated as in-memory (exact = NA). KIT-098 keeps three assertions: every trx case binds `oracles`, the lifecycle pair included; every trx step with a state, and the RUVIII product, has an oracle, and a missing one makes the member `invalid` (`oracle_missing`); acceptance requires `in_memory = FALSE` and `head_unchanged = TRUE` on the failure member. KIT-092 is added as a dependency. |
| F03 | should-fix | Closed in ticket | Every route-specific site is listed with its line: kit (`:15`, `:111`, `:178`, `:194-207`, `:251-268`, `:482-486`, `:622`, `:729`), driver (`:196`, `:506`, `:538`), coordinator (`:18`, `:298-316`), runner (`:6-8`, `:85`) and the generator. The three silent fallbacks (DIA step map, metab driving, the 600 s lifecycle ceiling) become refusals or trx entries. Trx gets its own lifecycle spec and fault, each covered by a harness test. |
| F04 | should-fix | Closed in ticket | Each KIT-095 worker step is its own journey step. If KIT-097's GUI runs several in one action, the trx verifier compares every history state by state name instead. The comparator also reads the intermediate products from the journey store. |
| F05 | should-fix | Closed in ticket | v1.8 names the oracle under §7.1, not as a control. The notebook check is defined as the second, implementation-independent check, read from the journey store. The ticket states there is no memory, cross-build or real-data claim, and the scientific contract says what each check proves. |
| F06 | should-fix | Closed by numbering, plus ticket text | KIT-092 takes v1.7 and KIT-098 takes v1.8 (user decision 5). v1.8 is decided after KIT-092 and before any KIT-098 member is bound. The evidence for backward compatibility is now a `bindings.json` diff limited to the trx keys, the protocol hash and the harness hashes, plus the existing route tests. Earlier receipts stand under their own protocol hash. |
| F07 | should-fix | Closed in ticket | v1.8 lists the worker steps. The trx decision checks one completed worker journal per listed step, reusing the `kitGuiMemberDecision` scan (`:967-985`). A missing journal makes the member `invalid` (`worker_dispatch_missing`). It is an acceptance criterion, with a harness test. |
| F08 | should-fix | Closed in ticket | The budget now names: all members recorded and uncounted; the selectors (`kit:lifecycle:trx:*`, `kit:journey-098a`/`098b:candidate:1`); two in-process recordings of at most 1,800 s each; the comparator as a read-only step, not a run; one attempt each at 1,800 s; the terminal outcomes (a rejection is recorded and never retried; the correction pass only for a demonstrated instrumentation failure); and a stop if the user declines v1.8. The run-ledger row is the coordinator's to sync in the plan. |
| F09 | should-fix | Closed in ticket (KIT-098 owns it) | KIT-098 binds the trx screen corpus. The generator asserts each file's sha256 equals the KIT-040 reference's recorded inputs. The coordinator gets the trx corpus fields. Trx uses its own bound normalisation and DA settings, not the shared `journey$normalisation`. |
| F10 | minor | Closed in ticket | The recorder reads its settings from the bindings. It records the `workflow_kit.R` sha256, the R/BLAS/LAPACK identity, the per-step variants and the protocol in force. The oracle labels name the protocol in force. |
| F11 | minor | Closed in ticket | The README chain check (items 1-6) is a Simulation Setup item, recorded in process before the lifecycle pair, reusing KIT-095's tests where they cover an item. |
| F12 | minor | Closed in ticket; one item out of scope | "793 lines" is correct (no trailing newline) and the shared contract text is left as is. `bindings.json` is added to `files_to_modify`. The order of DA and the session export follows KIT-097's GUI, cited in `plan.md`. **Out of scope, routed to the coordinator:** the plan table row for KIT-088 (`00-KIT-PLAN.md:240`) has stale dependencies. |
| F13 | minor | Closed in ticket | `kitOperationIdentityProducts` stays a harness constant (coordinator routing, together with KIT094-F07). KIT-098 extends it for the trx `*_science` identity products only, with the cause stated as KIT-031 did. Every other trx product is judged, and a harness test shows one is. |

## Counts

| | Count |
|---|---|
| Findings | 13 |
| Closed in the ticket text | 11 (F01, F03 to F11, F13) |
| Routed to KIT-092, with residual assertions kept in KIT-098 | 1 (F02) |
| Closed with one out-of-scope item for the plan | 1 (F12: the KIT-088 dependencies in the plan table) |
| Not closable here | 0 |

Execution-layer closure, meaning the tests, receipts and comparator output, waits for KIT-098's own run, after KIT-092 and KIT-097.
