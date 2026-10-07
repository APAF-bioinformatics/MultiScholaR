# CORE-001 evidence (2026-10-08)

Ticket: `tickets/openai/CORE-001.md` revision 1.1. Source base: `9df060f00a5d8e59b712d05c86d31c0165980035`. No counted runs were used.

Test environment: `pkgload::load_all()` from source with `.libPaths(c(<control jEmE2F library>, <default library>))`, as recorded in [bindings.json](bindings.json), because the default library lacks `fastduckdb`, `duckplyr` and `nanoarrow`.

## Bindings

`bindings.json` records the PERF-107 paused commit and consistent `blocked` status, the control library and build `0c38dce1…626bd`, the `full$p14` (24,000,000 rows), `screen$p14` (960,000) and `full$primary` (6,225,000) population identities, the P14 oracle (1,200 full and 48 screen page hashes) and the primary filter oracle. All 45 runtime files bound by the control build match their recorded SHA-256.

## Repairs and falsifying tests

Each test was run against the unrepaired code first and failed, then passed after the repair.

| Repair | Test | Before | After |
|---|---|---|---|
| `releaseWorkflowTableReadScope` releases everything, raises the first error last | resources: "a failing cursor close still releases the scope, its action and connection" | 6 failures, 1 error (connection open, action lease held) | 9/9 pass |
| `workflowRegisterDIAOrderedSum` retains the engine handle after success | bulk-relational: "the exact-sum provider registers on repeated native scopes in one process" | not reproducible from R (see below) | 9/9 pass |
| `publishWorkflowLegacyExport` resumes from `latest.replacement` | exports: "latest publication resumes from latest.replacement after an interrupted rename" | error "published latest session differs" | 6/6 pass |
| Staging and completion journals written atomically; a failed write removes its temporary | exports: "an interrupted journal write never leaves an unparseable export journal" | journal left as `{` (parse error) | 11/11 pass, including no orphaned temporary journal |
| `publishDIAActionOutputs` clears unjournaled staging before `produce()` | exports: "a retried action publication never publishes files from a failed attempt" | stale file published | 4/4 pass |

**Engine handle.** The premature `dlclose` is latent: `fastduckdb` keeps the engine library loaded, so the defect cannot be forced from R without C-level fault injection. Retention after success, the engine-identity check before any pointer is written, and pointer reset on a first-registration failure are verified by code review of the diff, not by a test. The kernel callback body remains byte-identical to `src/dia_rsum.c` (build-time check).

Additional exact-sum evidence: `test-perf108-native-null.R` with `MULTISCHOLAR_PERF108_NATIVE_NULL=true`: 34/34 pass on the rebuilt provider (missingness, constant and grouped lists, ordered sums, refusal of negative, infinite and NaN values).

## Regression checks (full files, from source)

| File | Result |
|---|---|
| `test-workflow-saved-resources.R` | 302 pass, 0 fail |
| `test-workflow-saved-exports.R` | 276 pass, 0 fail, 4 skip (re-run after the temporary-file cleanup) |
| `test-workflow-bulk-relational.R` | 9 pass, 0 fail, 3 skip (harness-gated preflights) |
| `test-workflow-saved-recovery.R` | 58 pass, 4 fail, 2 error |

The recovery failures are identical on the unmodified HEAD (`9df060f0`) run in a detached worktree with the same runner: the two `callr` subprocess cases cannot find `MultiScholaR.so` under a source load, the pre-existing limitation recorded by WFX-003. They are not counted as passes and not repaired here.

`R CMD INSTALL` of a clean snapshot (tracked files plus the edits) into a scratch library: exit 0; `workflow_dia_exact_sum.c` compiled with the package rule `-std=c11 -O2 -fPIC -fno-fast-math -Wall -Wextra -Werror`; `configure` regenerated `manifest.rds`, which runs the kernel-body equality check.

## Open item

The ticket's acceptance names "the existing five exact P14 cases". Those live in the harness-gated relational preflight (`gui-bulk-rollup:relational-preflight`), which needs an installed library and a proof-root receipt through the guarded runner. CORE-001's budget allows no counted runs, so they were not executed. The change does not touch arithmetic, ordering or SQL; the end-to-end exact rollup is next exercised by CORE-004's lifecycle members against `p14_oracle$screen_pages`.

## One-off guarded preflight (user-authorised 2026-10-08)

Terms agreed before launch: one attempt, 5-minute ceiling, existing guards; a harness or setup failure is recorded and falls back to CORE-004 evidence with no harness repair or retry; an exact-case mismatch would block CORE-001 and revert the C change. The plan precondition was lowered by the user to `MemAvailable` ≥ 12 GiB for this run only.

| Step | Result |
|---|---|
| Install (`install-corrected`, member `core/browser-smoke.zPIbOm`) | exit 0; source `9df060f0` plus a runtime patch of exactly the three repaired files; kernel source identical to the working tree; schema-3 manifest with native query binding; library holds MultiScholaR, fastduckdb, duckplyr, nanoarrow |
| Preflight (`reader-check … gui-bulk-rollup:writer-r-preflight`) | exit 76 before launch: `run_workflow_checkpoint_proof.sh` refuses bulk-rollup selectors when host headroom is below 16 GiB (host `MemAvailable` ≈ 15.8 GiB). No member directory or systemd unit was created. |

Outcome: setup refusal under the runner's own 16 GiB startup guard, not a scientific result. Per the agreed terms it is not retried and the runner is not changed. The five exact P14 cases fall back to CORE-004, whose acceptance now requires the installed rollup on the repaired provider to match the bound oracle.
