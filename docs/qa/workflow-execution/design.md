# Managed DIA execution: current design

A short orientation for someone who has not read the qualification logs. It describes the code as of the CORE series (CORE-015). For the platform requirements see [preconditions.md](preconditions.md).

## The pieces

```
 GUI session (Shiny, one R process)
   │  runWorkflowAction(...)             pollWorkflowAction(...) every 50 ms
   ▼
 Controller ── journal + control files ──► project store (.workflow-entities/)
   │                                             entities/<kind>-<id>/  payload + descriptor
   │  systemd-run --user                         dia-state-head.rds      current project head
   ▼
 ┌──────────── action slice (16 GiB cap, no swap) ────────────┐
 │  monitor unit: samples memory every 50 ms, stops the action │
 │                at 15.5 GiB or below 2 GiB of host memory    │
 │  worker unit:  MemoryMax = admitted cap; runs the science   │
 └─────────────────────────────────────────────────────────────┘
```

- **Controller** (in the GUI process). It admits an action, launches the monitor and the worker, polls for the result, and publishes it as the new project head. Every step is recorded in a request journal so a crash can be recovered (`utils_workflow_execution.R`).
- **Admission.** The worker's memory cap is the smaller of what the 16 GiB slice leaves and what the host has free (minus a 2 GiB floor and a 1 GiB growth reserve), so a busy machine gets a smaller cap rather than a crash (`workflowAdmissionCap`, CORE-010).
- **Entity store.** Every saved table or value is an immutable entity, written to `.pending/` and then renamed into `entities/`, with a descriptor that records its hash. The project head points at the current state snapshot, which references every entity the project needs. Save archives exactly that closure; Open extracts it into a fresh private folder (`utils_workflow_saved_entities.R`, `utils_workflow_dia_sessions.R`).
- **Table codec.** Tables are Parquet files read in pages of at most 10,000 rows or 64 MiB through fastduckdb, and written by an incremental native writer (`utils_workflow_saved_table_codec.R`).

## The precursor-to-peptide rollup (P14)

1. **Validate** the whole input: identity keys, finite quantities, duplicate identities and conflicting evidence. The grouping-heavy checks run once per batch of complete samples (CORE-017).
2. **Admit.** Estimate the working memory for 3, 4, … batches and take the fewest that fit the cap (CORE-010). With plenty of free memory this is three.
3. **Aggregate** each batch of complete samples in DuckDB on the native fastduckdb connection. The exact ordered sum `dia_r_sum(list(x ORDER BY source order))` reproduces R's `sum()` bit for bit (CORE-012 kept this over faster compensated sums).
4. **Stage** each batch's ordered result to a private Parquet file, close the engine, then stream the files into the output writer in 10,000-row pages.

The peak memory is in step 3, in DuckDB's grouping of the ordered value lists. That state cannot spill to disk (CORE-011), so batching by sample is the lever that reduces it.

## Where measurements come from

- **Per action:** the monitor's 50 ms samples and the kernel's `memory.peak` for the slice, written to the action receipt.
- **Qualification runs:** the CORE harness (`tools/profiling/workflow_core.R`, `workflow_perf106_repair.R`) adds phase timings, phase peaks (cache-inclusive and anonymous) and an independent verifier that checks every output page against frozen oracles. The rules are in `core/protocol.md`.
