# CORE-005 evidence: hashing memoisation and writer metadata cache

Ticket: `tickets/openai/CORE-005.md` revision 1.1. Commit `c0de46ea`. No counted runs. The installed build is recorded as `core_005` in `bindings.json`; its non-regression is CORE-006's control pair, by ticket design.

## Premise checks (before editing)

| Item | Premise | Verified |
|---|---|---|
| 1 | The fast writer, extension and native capability checks hash binaries and dependency files on every open | Yes. `workflowFastWriterInfo` hashes three package files and the DLL; `diaRollupExtensionInfo` hashes the extension and the duckdb/duckplyr files; `diaRollupNativeInfo` hashes the query manifest, the query binary, the duckdb files and the provider DLL. Note: `workflowTableFileIdentity` records path, size, mtime, **ctime** and mode, not an inode as the ticket says. ctime changes on any write or rename-over and cannot be set back, so it is the stronger key. |
| 2 | Input files are re-hashed at every commit | Yes. `describeDIANativeFile(path)` is called on `iq_output_file` and `peptide_input_file` in `commitDIAReferenceState` and `saveDIAStateRecord`. With `log_path = NULL` it is `digest::digest(file = …)`. |
| 3 | Writer metadata is re-serialised on every append | Yes. `schema_metadata` was rebuilt per batch and kept only under `optimization = "metadata"`, which only profile P15 passes. On `fast_r` it is used only when the writer opens, so later batches discarded it. A v1 table's metadata carries every row name, so the cost grows quadratically. |

## Changes

1. `workflowFileSha256(path)` in the codec keeps a per-process environment keyed by normalised path, storing the identity and the SHA-256. It returns the cached digest only while the identity is unchanged. It caches only when the identity is the same before and after hashing. It falls back to an uncached hash when the identity cannot be taken (for example, a symlinked path), so capability checks fail exactly as before. Only digests are cached: every capability comparison runs on every call, so a failed check is re-evaluated. All eleven file digests in the three capability checks use it.
2. `describeDIAInputFile(path)` in `dia_state` returns the same `list(bytes, sha256)` as `describeDIANativeFile`, with the same directory refusal and the same changed-during-validation abort, and uses the cache. Only the two commit sites use it. `describeDIANativeFile` is unchanged: the archive and copy validations in `utils_workflow_dia_sessions.R` (outside this ticket's files) keep hashing every time and keep their IO log.
3. `appendWorkflowTableWriter` stores the encoded metadata on the writer after the first batch. The `metadata` value is removed from `openWorkflowTableWriter` and from the P15 profile switch; P15 now maps to the default. `batch` (P16) is untouched.

## Tests (source runner, control `jEmE2F` library first)

The tests were written after the code, against the failing-test-first rule. To compensate, each was checked against a mutant restoring the old behaviour:

| Test | Result | Mutant (old behaviour) |
|---|---|---|
| file hashes are computed once per unchanged file identity (same size and restored mtime still rehash, through ctime) | pass | no cache: fails (2 digests instead of 1) |
| two fast writer capability checks hash each runtime file once (4 digests, then 0) | pass | no cache: fails |
| a writer encodes its restoration metadata once across batches (v1 table, 50,000 row names, five 10,000-row batches) | pass (1 encoding) | old `metadata` guard: fails (5 encodings) |
| a commit-time input descriptor matches the archive descriptor and is hashed once | pass (identical to `describeDIANativeFile`; 0 reads on the second call; changed content gives the new descriptor) | no cache: fails |

The native and extension capability checks call the same `workflowFileSha256`; no separate test opens `newDIARelationalScope` twice.

**Byte identity (acceptance item 3):** the same v1 table (50,000 row names, a factor) and v3 table (30,000 rows, NA) were written with `fast_r` in 10,000-row batches. The installed control build `jEmE2F` (old code) and the source tree (new code) produced identical file hashes: v1 `76a3f9c9…`, v3 `b128d911…`.

## Regressions (source)

| File | Result |
|---|---|
| `test-workflow-saved-resources.R` | 16 tests, 313 passed |
| `test-workflow-saved-identity.R` | 11 tests, 77 passed |
| `test-workflow-bulk-relational.R` | 4 tests, 9 passed, 3 installed-only skipped (including the edited writer case, which now asserts cached metadata for every optimization) |
| `test-workflow-saved-installed.R` | 1 test, installed-only skip |

## Install (recorded, not counted)

`core/browser-smoke.UgUYyi`: source `c0de46ea`, empty patch, exit 0. Build `693afda9…`; exact-sum identity `4f0a90eb…`. The provider `.so` is recompiled per install, and its hash is part of the identity, so this differs from `core_004` although the C source is unchanged. Bound in `bindings.json` as `core_005` with outcome `pending`. `bindings-generator.R` now takes `key=dir` pairs, so earlier builds are kept.

## Deferred by ticket design (not done here)

- **Exact P14 cases on the installed candidate, and the outcome record.** CORE-005 has no counted run. CORE-006's control pair runs full P14 on `core_005` with oracle verification. It records CORE-005 as `accepted` (peak ≤ 1.00 and elapsed ≤ 1.05 against the PERF-108 `full_p14` receipts) or `reverted`. Until then `core_005` is not an accepted control.
- **Measurement.** No measurement claim is made; caching is expected to save hashing time only.

## Update (2026-10-08)

CORE-006 closed without change (see `core-006/outcome.md`), so CORE-007's control pair on `core_005` qualifies this build and records its outcome in place of CORE-006's.
