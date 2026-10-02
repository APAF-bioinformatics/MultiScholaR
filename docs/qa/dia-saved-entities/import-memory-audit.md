# DIA import and active-owner memory audit

Date: 2026-09-23. Owning ticket: DIA-025. Mode: **audit only**.
No runtime, harness, scientific-method, threshold or ticket-status changes are
authorized by this audit. This report and its JSON companion are additive findings,
not a new repair series or performance acceptance. The worktree was already dirty.

## Conclusion

**A concrete shared lifetime defect is confirmed.** The experimental adapter saves
large objects to disk, but also pins every saved/resolved record in an in-memory
working map until the entire Shiny flush finishes. Replacing a GUI owner does not
evict its old working-map entry. A nested action scope does not release that entry
while the outer flush is active. A compound import/design update can therefore
accumulate retained versions before reaching its release callback.

This is not evidence that every record occupies a separate full-size allocation:
R columns can share memory. It is evidence that obsolete owners remain reachable,
and repeated whole-record serialization creates distinct disk files/cache pages.
The exact fraction of each large-run peak attributable to this mechanism is not
measured. No claim that one repair solves the entire workflow is supported yet.

There is also an important architectural limit: **import is not an out-of-core
Parquet/DuckDB pipeline**. The DIA reader returns a complete R table; import/design
owners use whole-record qs2 storage. The Parquet/DuckDB projection is prepared
later for precursor roll-up. Saving a whole table is not a selective metadata lookup.

## Evidence and applicability

- Installed candidate: xz9cb6, RDB SHA-256
  `74a7dc9b50821e49e1a6d8c9a4366bdcc879997d61d06a9876acd55387404a08`.
  Seven key audited source files match the installed source snapshot by SHA-256;
  the hashes are in the companion JSON. Both stress runs loaded this package path.
- 500 samples: iGvvnA, prepared input 570cC2, 40 million rows, report
  5,792,259,275 B. Early guard observed 16,916,217,856 B, zero swap; controlled
  stop after 99.292 seconds. One committed `owner:data_tbl` payload, no FASTA outputs.
- 300 samples: 8WgDVe, input OQADra, 24 million rows, report 3,475,378,370 B.
  Early guard observed 16,684,773,376 B, zero swap; controlled stop after 145.151
  seconds. FASTA/configuration outputs and design-initialization payloads exist.
  **This run reached design initialization, despite the last flushed GUI receipt
  still saying `initial`.** It did not reach `raw_data_s4` or QC.
- Last complete pre-stop sample: 500 has 9,901,445,120 B anonymous and
  6,681,321,472 B file memory; 300 has 11,560,730,624 B anonymous and 5,069,545,472 B
  file memory. These asynchronous counters are not an exact additive peak partition.
  Anonymous memory covers processes, not just the R heap. No per-file cache or
  per-function allocation trace was collected in these unprofiled runs.
- The FASTA is the same 710,000-byte, 10,000-sequence input in both cases. The
  retained 300-sample annotation table has `object.size` 4,762,520 B and its metadata
  2,472 B. These are R object-size measurements, not transient/RSS peaks.
- The already-running clean 96-sample history verifier GjvReR subsequently passed
  all 16 GMUgZW/v6Blcs states at the existing 1e-12 tolerance and exact attribute/
  type/missingness checks, peak 6,113,869,824 B and zero swap. X9hWYX's earlier
  failed process exit remains retained. This does not qualify the larger failed runs.

All run artifacts are under
`/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/`.
Large qs2 payloads were **not decoded during this audit**. Descriptor inventories
and small retained FASTA objects were inspected. No new full journey was launched.

## What the actual path does

```text
Complete DIA-NN R table
  -> materialize parser diagnostics for portability
  -> save original owner:data_tbl
  -> sanitize Run and save another owner:data_tbl
  -> FASTA/configuration and import summary
  -> save design-initial bundle containing the whole cleaned table
  -> save design-cleaned-data containing that table again
  -> metadata/status consumers and other reactive work
  -> onFlushed: release the working map

Saved records remain reachable in the working map across the middle steps.
The original import result also remains live through the import callback.
```

The 300-sample committed descriptor inventory contains:

| Role | Codec | Payload bytes on disk |
| --- | --- | ---: |
| owner:data_tbl, first version | qs2 | 515,239,720 |
| owner:data_tbl, second version | qs2 | 517,622,295 |
| owner:design-initial | qs2 | 516,806,416 |
| owner:design-cleaned-data | qs2 | 517,621,779 |

These are four **logical whole-table records**, totaling 2,067,290,210 compressed
bytes on disk, not a claim of four independent full-table RAM allocations.

## Findings and closure requirements

### DIA-IMPORT-01 — High: release boundary is too broad

Evidence: `R/utils_workflow_dia_state.R:47` retains saved qs2 values;
`:83` reuses/retains decoded values; `:585` replaces the owner without evicting its
predecessor; `:603` does the same for saved reactive values; `:676` sets a scope
around the entire flush and `:679` clears it only after output delivery.
`withDIAStateScope` at `:368` cannot empty the map while the enclosing depth is one.
The optional store cache ceiling does not bound this separate working map.

Small deterministic control-flow check: actual owner/scope/callback helpers,
two-row tables, persistence/validation doubled to avoid disk/data workloads.
After owner replacement: two working records. After another save inside a nested
action scope: three records, outer depth one, original record still reachable.
After the actual flush callback: zero records. This reproduces lifetime behavior,
not large-data performance or unique allocated bytes.

Smallest repair: define consumed-owner release at the end of its actual use,
preserving committed references and other live aliases. Do not merely add another
nested scope, clear the optional cache, run gc(), or delete durable checkpoints.
Closure: active-action tests must observe superseded references released before
subsequent heavy work, preserve failed-save rollback and aliases, and keep
revert/native restore intact. Blocked state: implementation/performance promotion.

### DIA-IMPORT-02 — High: the import front door still materializes the full table

Evidence: `R/func_prot_import_readers.R:18` uses unrestricted `vroom()` or
`arrow::read_parquet()`, returning the full data table. `R/utils_workflow_saved_entities.R:798`
calls `readr::problems(table)` to make the parser pointer portable. Installed
readr 2.2.0 delegates to vroom 1.7.1; its `problems(x, lazy=FALSE)` calls
`vroom_materialize(x, replace=FALSE)`. The save boundary thus forces lazy parser
columns to materialize; it is not a cheap reference-only flush.

Also, `R/mod_prot_import_detection_helpers.R:43` reads an entire Parquet table
just to obtain column names. That second issue is **not exercised by these TSV
stress runs**, but means changing their input extension alone is not a solution.

The existing query projection starts at
`R/utils_workflow_dia_state.R:699`, after the protein/peptide QC result, not import.
Smallest repair: separate metadata/schema discovery from full materialization;
if import still exceeds the limit after lifetime repair, introduce bounded
ingestion through the existing table codec and reduced metadata queries.
Keep standalone public reader behavior compatible. Do not skip parser-error
validation merely to retain ALTREP laziness. Closure: original parser diagnostics,
column types/attributes/order and missing/nonfinite semantics plus complete raw
input preservation and complete-Apply memory/time evidence. Blocked: implementation.

### DIA-IMPORT-03 — High: original and sanitized owners overlap unnecessarily

Evidence: `R/mod_prot_import_processing_helpers.R:114` saves the original before
`:119` sanitizes. `:86` resolves it and `:95` replaces the whole Run column before
`:96` saves again. Sanitization defaults TRUE (`R/mod_prot_import_ui.R:93`).
`R/mod_prot_import_server_helpers.R:292` retains `data_import_result`, uses it for
later logging/finalization and returns it at `:399`, after FASTA/configuration.
The 300-sample store confirms both serialized import versions.

The named lookup used for Run also attaches names to the replacement vector,
including when visible run labels are unchanged; a tiny base-R example confirms
that distinction. A no-op shortcut must not silently change those attributes.

Smallest repair: compute the required small metadata and sanitization once before
the necessary durable handoff; release the original full result when no consumer
needs it. Retain exact sanitization collisions, order, attributes and diagnostics.
Closure: same returned table and setup log, fewer unnecessary whole-table writes,
and bounded active owners on success/failure. Blocked: implementation.

### DIA-IMPORT-04 — High: design metadata carries and resaves the full assay

Evidence: `R/mod_prot_design_builder_server_helpers.R:134` builds and saves a
`design-initial` value. `R/mod_prot_design_builder_state_helpers.R:344` embeds
`data_cln` in that bundle. The sync observer resolves that bundle and at
`R/mod_prot_design_builder_server_helpers.R:54` saves its entire table again as
`design-cleaned-data`. Both large files are present in the failed 300-sample run.
The initial design matrix itself needs only distinct run IDs and optional Batch
metadata, plus small group/factor/contrast fields.

Smallest repair: keep that metadata separate and reuse an immutable table
reference until a real data edit requires a successor. Register dependencies
through the existing store/export lifecycle; do not create opaque orphan links.
Closure: initialization/reset/rename/removal/assignment/save/revert preserve
exact results; unchanged data is not reserialized simply to populate another
holder, and all referenced data survives portable session restore. Blocked: implementation.

### DIA-IMPORT-05 — Medium: status and metadata consumers ask for complete payloads

Evidence: `R/mod_prot_import_server_helpers.R:165` and `:173` call
`getDIAOwnedPayload`. Status uses only existence; summary uses the already-small
setup log (`R/mod_prot_import_ui_helpers.R:1` and `:22`). Design availability and
initialization likewise evaluate `dataTbl()` for requirements/events, not only
actual table transformations. The qs2 getter returns a complete record, not a
projection. Same-flush reads of one identity may reuse the working map: **one
getter call does not necessarily mean one decode**. Across flushes, the map clears.
The processing-status short-circuit can avoid forcing the lazy argument while busy.

Smallest repair: use reference presence and bounded metadata/sample summaries for
these consumers, and hydrate only actual scientific/compatibility consumers.
Closure: status/summary/design-choice invalidations perform zero full table decodes,
while their displayed contents and dependency invalidation stay correct.
Blocked: implementation.

### DIA-IMPORT-06 — High: checkpoint I/O has a separate file-residency cost

Evidence: `R/utils_workflow_saved_entities.R:714` writes whole qs2 objects;
`:994` fingerprints the newly written file and `:1017` synchronizes it.
Recovery/validation at `:863` rereads payloads for SHA-256 and synchronizes them.
These paths do not use the bounded/cache-advised I/O added specifically to native
archive export. The large runs stop with approximately 4.72/6.22 GiB of file
memory in their last complete samples, before native export is reached.

Smallest repair: first eliminate unnecessary files/reads, then measure whether
file-specific drop-behind is appropriate for consumed checkpoint I/O. Preserve
durability, genuine external-import integrity and recovery inputs. Do not globally
flush caches or remove file memory from accounting. Closure: separately report
process/file memory and complete elapsed/physical-read tradeoffs; measured gain,
not inferred gain from smaller R objects. Per-file peak attribution remains open.
Blocked: performance decision.

### DIA-IMPORT-07 — Low here: FASTA processing has waste, but not the measured scale

Evidence: `R/func_prot_annotation_fasta_processing_helpers.R:192` reads the file
with vroom to inspect one header, then reads it again with seqinr. Parsing builds
per-entry lists and tabular outputs; `fasta_meta_file` is written rather than
checked for reusable cached content. `R/mod_prot_import_processing_helpers.R:271`
adds a global annotation alias and writes another RDS plus metadata. Non-standard
FASTA also writes a TSV. Global aliases can extend lifetime; they are not by
themselves proof of duplicate physical allocations.

Smallest repair: bounded first-header inspection and, only if repeated parsing
is measured, cache reuse keyed by source bytes/parser/version/mapping inputs.
Preserve sequence/header/accession/organism semantics and required files.
Closure: parser/annotation parity across supported FASTA variants and measured
FASTA-only benefit. The current 710 KB synthetic FASTA/4.54 MiB annotation result
does not explain a multi-GiB, sample-scaled spike. Large real FASTAs remain untested.
Blocked: optimization justification, not primary capacity diagnosis.

### DIA-IMPORT-08 — High: post-flush counters do not observe the failing interval

Evidence: `attachDIAStateSession` writes its ownership receipt after clearing the
working map (`R/utils_workflow_dia_state.R:679`). Both failures retain an old
`initial` receipt with zero working records, while their committed owner files
prove subsequent heavy work. Existing lifecycle assertions after delivery can
pass despite an unsafe peak within the flush. The cgroup peak/guard evidence is
still valid; the active-owner attribution is missing, not zero.

Smallest repair: add bounded scalar observations at import/save/sanitize/design
boundaries in the existing diagnostic machinery; no deep heap walks or continuous
allocation logs. Closure: tests reproduce intra-flush replacement and failure,
count distinct live references during the action, and retain pre-stop attribution.
Blocked: attributing exact peak shares and claiming complete architecture coverage.

## Reach beyond import

The same setter is used by normalization and RUV. RUV stores the corrected S4 at
`R/mod_prot_norm_ruv_helpers.R:295`, replaces it after cleanup at `:375`, then
saves the scientific state at `:380`. That is another concrete multi-owner path
exposed to the lifetime rule; no per-stage savings are claimed from inspection.

There is already a working local precedent: Q-value and protein/peptide evidence
Apply explicitly call `releaseDIAWorkingReference` after a successful successor
commit (`R/mod_prot_qc_peptide_qvalue.R:219`,
`R/mod_prot_qc_peptide_protein.R:157`). Their tests observe release during the
action, not only afterward. Other active matrices/model buffers and the separate
native-export file-cache issue are not explained away by this finding.

The importer remains synchronous in the GUI observer. Moving the unchanged full
table to another process alone would not prove lower combined memory. Existing
GUI configuration correctly uses `create_workers=FALSE` for the active adapter;
there is no evidence here that the previous eager-worker-pool defect returned.
Developer checkpoint capture is explicitly disabled by the harness; it is not
the source of these duplicate owner records.

## Bounded repair priority, not authorization to implement

1. Correct the consumed-owner lifetime and remove redundant import/design
   whole-record handoffs, using the existing store and action boundaries.
2. Make status/design metadata consumers reference-based. Preserve source table
   contents and scientific populations; no arbitrary row/column dropping.
3. Verify active-step owners and total peak on a small real GUI lifecycle, then
   repeat an approved retained workload with the same guards and outputs.
4. Only if active import itself remains the bottleneck, implement bounded
   columnar ingestion/selective lookup. Audit FASTA separately when its cost is
   measured. No new backend is justified by the current evidence.

All repairs must preserve exact identifiers/order, factor and other attributes,
sanitization mapping, quantity-column choice, missing/nonfinite values, parser
errors, annotations, design/exclusion denominators, and complete recoverable
checkpoints/portable exports. No FDR/model/population or file-format change follows
from this audit. Keep implementation checks, peak attribution, science and full
workflow performance as separate gates.

## Audit artifact checks

JSON parsing, required top-level/nested fields, critical types/enums, unique
finding IDs and audit-only/no-repair flags pass explicit checks. A full JSON
Schema validator is not installed (`jsonschema` CLI/Python module and R
`jsonvalidate` unavailable); none was installed for this audit. Whitespace checks
pass. No application/harness changes or new long runs were performed; only this
report and its machine-readable companion were created.
