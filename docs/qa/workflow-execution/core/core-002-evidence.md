# CORE-002 evidence (2026-10-08)

Ticket: `tickets/openai/CORE-002.md` revision 1.1. Base: `2e4887df`. No counted runs. Test runner as recorded in [bindings.json](bindings.json); baseline failures per [test-baselines.json](test-baselines.json) (none in this ticket's files).

## Outcome per item

| Item | Outcome | Evidence |
|---|---|---|
| 1. Undeclared execution-view slots | **Criterion not implementable as written; contract pinned instead (user decision).** Kernels read slots with `@`, which R cannot intercept without changing the S4 classes; narrow-view callers legitimately read the zero-row prototype's column schema (for example `names(object@peptide_data)` in the rollup adapter); `checkpointDIAExecutionResult` already uses `dia_execution_fields` so unloaded slots never replace stored components. All 19 `readDIAExecutionComponents` call sites were read (about 12 during the ticket, the remaining 7 in the [post-hoc audit](../../../../tickets/openai/CORE-002-POSTHOC-AUDIT-20261008.md)). None computes scientific results on an undeclared slot. Three intentionally run a native kernel on zero-row prototypes to resolve parameters and output types (`resolveWorkflowActionParameters`, `computeDIAReferenceSampleFilter`, and the probe in the filtering summary), so a refusal would break working paths. | New test "narrow execution views carry schema-only prototypes that never replace stored components" (8/8); roxygen of `readDIAExecutionComponents` documents the contract. No behaviour change. |
| 2. Non-order-exact DuckDB `sum()` arm | `computeSavedWorkflowPrecursorRollup(..., non_exact = FALSE)` refuses the persist-less arm with class `dia_rollup_non_exact`; with `non_exact = TRUE` it returns `execution$aggregation_frontend = "duckdb_sum_non_exact"`. Its one production caller, the experimental capacity trial (`computeDIACapacityRollup`), passes the flag explicitly. | Failing first, then 7/7: refusal, marker, released leases, and that `runDIAReferencePrecursorRollup`, `computeDIAReferencePrecursorRollup` and `computeDIARelationalRollup` never reference `non_exact`. Two existing tests now pass the flag explicitly. |
| 3. Blank factor keys | **Review premise false; no code change.** The native kernel flags blanks only in character columns (`.rollupAssertPopulatedKey`: `is.character(values) & !nzchar(trimws(values))`), so blank factor labels are accepted by design, and both SQL validators already matched it. A cast added in this ticket made both routes stricter than the kernel; it was caught before commit and reverted. | Parity test "blank keys follow the native kernel": blank factor label accepted by kernel, diagnostic arm (identical output) and `validateDIARelationalInput`; blank character key refused by all three (7/7). |
| 4. Group order under non-alphabetic levels | Indexed arm already matches the kernel; no change. | Test with reversed, ordered and unordered `Run` levels: output identical to the kernel (4/4), passing from the first run. |

## Regression checks (full files, from source)

| File | Result |
|---|---|
| `test-workflow-saved-qc-rollup.R` | 76 pass, 1 skip (receipt-gated frozen-kernel test) |
| `test-workflow-checkpoint-dia-state.R` | 238 pass |
| `test-workflow-saved-qc.R` | 108 pass |
| `test-workflow-saved-consumers.R` | 36 pass |

## Observation for the user (not changed)

Whether a blank factor label should be a valid rollup key is a scientific question about the native kernel's contract. It is preserved here because every route must match the kernel; changing it would be a separate, kernel-level decision.

## Post-hoc audit corrections (2026-10-08)

The item 1 call-site statement above was corrected (finding CORE-002-P01). The item 3 closure was ratified by the user on 2026-10-08 (CORE-002-P03).

## Observation verified (2026-10-08)

The post-hoc audit's out-of-ticket observation was checked by running the kernels on the small fixture with one blank `Protein.Group` label.

| Blank protein label | Native rollup helper | Native sample-filter method | SQL sample filter (physical) |
|---|---|---|---|
| Factor | accepted | refused | refused |
| Character | refused | refused | refused |

The SQL sample filter matches its native kernel, so there is no parity defect in the SQL route. The inconsistency is between the two native kernels: `.rollupAssertPopulatedKey` tests blanks only when `is.character(values)`, while `filterMinNumPeptidesPerSampleHelper` applies `as.character()` first. A blank factor label therefore passes rollup and is refused, with a clear error, at sample filtering. Changing either kernel is a scientific decision outside CORE, which preserves native kernels; it is recorded for the owner of the peptide QC kernels.

