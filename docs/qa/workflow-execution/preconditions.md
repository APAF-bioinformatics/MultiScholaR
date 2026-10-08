# Managed DIA execution: platform preconditions

What the managed DIA workflow (P14 rollup, guarded actions, entity store) needs from the host, and what a user sees when a precondition is missing. Checked against the code at the CORE-015 commit.

| Precondition | Why it is needed | Where it is enforced | Refusal |
|---|---|---|---|
| Linux on x86-64 with an 80-bit `long double` | The exact rollup sum is defined by sequential long-double accumulation, the same as R's `sum()` on this platform | `diaRollupNativeInfo()` | `dia_rollup_unavailable`: "native DIA rollup requires Linux x86_64 with 80-bit long double" |
| DuckDB R package 1.5.5 (`duckdb`) | The native exact-sum provider is built against this engine ABI | `diaRollupNativeInfo()` compares the installed package with the build manifest | `dia_rollup_unavailable`: "native DIA rollup requires the bound fastduckdb query runtime" |
| fastduckdb 0.2.0 at the pinned commit | Native query owner and incremental Parquet writer | `diaRollupNativeInfo()` (query manifest and binary hashes) and `workflowFastWriterInfo()` (writer capability) | `dia_rollup_unavailable` (query runtime) or `saved_entity_unsupported` (writer) |
| systemd user manager able to run `systemd-run --user` with `MemoryMax`, `MemorySwapMax` and `BindsTo` | Every action runs in its own memory-capped unit | No version check; the first `systemd-run` or `systemctl` call fails | `workflow_execution_unsupported`: "managed execution command failed: …" (the systemd output is included) |
| cgroup v2 with the memory controller delegated to the user manager | Memory counters and caps come from the unit's cgroup | No explicit check; a missing counter fails when first read | `workflow_execution_unsupported`: "managed execution command failed" or "action peak counter unavailable" |
| Linux 6.12 or later (writable `memory.peak`) | Each action resets the peak counter so its own peak can be measured | `openWorkflowPeakWindow()` | `workflow_execution_unsupported`: "cannot reset the cgroup memory.peak counter; managed execution requires cgroup v2 on Linux 6.12 or later" |
| `/proc/meminfo` with `MemAvailable` (Linux 3.14 or later) | Admission compares free memory with the job's needs (CORE-010) | `workflowHostAvailable()` | `workflow_execution_unsupported`: "host memory unavailable: /proc/meminfo has no usable MemAvailable" |
| GNU coreutils `sync` that accepts file arguments | Durable writes sync each written file and its directory | `syncWorkflowExportPaths()` | `saved_entity_unsupported`: "durable export sync is unavailable" |

Since CORE-015, `duckplyr` is not used by the package code; it moved from Imports to Suggests and is used only by one installed test and the historical harness driver.

Every refusal above is a classed condition, so callers and the GUI can show the reason instead of a generic failure.
