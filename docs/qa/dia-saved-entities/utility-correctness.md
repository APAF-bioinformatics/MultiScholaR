# DIA IQ saved-entity utility qualification

Status: `G1:iq` accepted for the exact standalone utility build and capabilities
recorded here. This is DIA-008 evidence, not a GUI integration, scientific-method
acceptance or performance result.

## Bound build and prerequisite

The prerequisite is `G0:iq` scope `reduced-development-v2` in
[gates.json](gates.json). It permits standalone DIA IQ saved-entity work only.
Frozen-v1 G0 and Limpa remain held.

The qualified candidate is MultiScholaR 0.5.0 installed from the current
`performance-v2` worktree into
`/tmp/multischolar-dia008.APxrbx/library`. The Git parent is janitor
`c7c12851ea4d9e96f91df0e1f9e7b91c2eb51017`; the worktree is deliberately dirty
with the ticket series, so that commit alone is not its identity. Relevant
bindings are:

| Artifact | SHA-256 |
| --- | --- |
| Installed `R/MultiScholaR.rdb` | `a5ef768af6bc73e84952652fe359f4e18b5c5c818354a6dd4bdcf1ee96339b9b` |
| Installed `R/MultiScholaR.rdx` | `64c3abddce1c264da1cab2e9e3c7a887bb48ffe13fc27999a36ccb7a126dacf8` |
| Installed `NAMESPACE` | `9dc606a759229234e12f5a0b55879af5b7230a42f1177a4d86fa1e088879ccc5` |
| Saved-entity utility source | `f0d4b0bfbf2e580f603135ddbdec2a4e852f976aaf55862b26b133c0bfa904c9` |
| Installed qualification test | `85efe759e6fe93c018d4d17d2db39193abdb4a90b323753a6e986ad5331185bf` |
| Installed qualification receipt | `ffed6798365dc97bfda148e9c2d4dd3cd24e2e0542c31dfccb71bea986c55742` |

The comparator for public package compatibility is the preserved pristine B0
installation at
`/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/b0-install-RNROup/library`.
It is not used as a runtime performance denominator. The external installed
receipt at `/tmp/multischolar-dia008.APxrbx/installed-receipt.json` binds both
package paths and the complete snapshots.

## Installed runtime qualification

One fresh candidate installation completed under a 180-second command ceiling.
The installed-only test then ran with testthat's stop-on-failure reporter in fresh
callr processes and exited zero. It compared candidate and pristine B0 and found:

- exactly 624 exported names in each package;
- exactly the same 188 loaded namespaces after the same package-load probe;
- identical WorkflowState public methods and fields;
- identical `runTestsContrasts` and `run_app` formals;
- identical PeptideQuantitativeData and ProteinQuantitativeData slot signatures;
- identical `normaliseBetweenSamples` and `differentialAbundanceAnalysis`
  generic formals and registered method signatures; and
- no MultiScholaR native library in either build.

The package declares no native routines of its own. The utility's required
`digest::digest` hashing path did load and exercise a real `NativeSymbolInfo`
(`digest_impl`) while saving the installed checkpoint. Thus absence of a
MultiScholaR DLL is asserted rather than silently treated as an untested routine.
Candidate and B0 `Depends`, `Imports`, `Suggests` and `LinkingTo` fields are
identical. The candidate's only metadata difference relevant here is the expected
DIA-006 system requirement for GNU coreutils file-specific sync; it adds no R
dependency or eager namespace.

The installed probe constructed and validated a real ProteinQuantitativeData,
dispatched `normaliseBetweenSamples(..., "none")`, and preserved identifiers and
nonfinite-to-missing behavior. It exercised the exported WorkflowState with two
states, external `states` mutation, revert, and a complete legacy RDS round trip.
It then prepared, saved, hydrated, cached, released, closed, reopened and hydrated
the complete history. Serialized values were identical. A foreign project
reference and a live connection were rejected with `saved_entity_scope` and
`saved_entity_unsupported` respectively.

The positive cache held one 19,048-byte entry, recorded one decode/read and one
hit, stayed below 64 MiB and reported zero open resources. A newly opened default
store reported a zero-byte cache ceiling. The request-owned test artifact had six
files totalling 11,226 bytes. This is a footprint check, not a memory benchmark.

The first installed test invocation exposed a missing-argument error in the new
test's compatibility-snapshot `lapply`; no candidate utility code had run. The
bounded correction removed those test-only empty arguments. The next invocation
passed all 22 test expectations and wrote the receipt, but its shell post-check
called the non-exported `testthat::is_failure` and exited after the passing test.
Both instrumentation failures are retained. The supported StopReporter invocation
then exited zero without changing or reinstalling the tested package.

## Source exactness, lifecycle and failure evidence

Command:

```text
Rscript --vanilla -e 'devtools::test(filter = "workflow-(checkpoint|saved)")'
```

Result: 1,102 passed, zero failures or warnings, two declared skips. Coverage was
99 checkpoint-fixture assertions, 474 harness assertions, 107 complete checkpoint
assertions, 73 identity assertions, 99 recovery assertions and 250 resource
assertions. The recovery group includes real child-process termination and
fresh-process reopen. The checkpoint and identity groups cover factors, types,
missing/nonfinite values, method-valid S4 objects, immutable identities, corruption,
conflicts and explicit unsupported inputs. The resource group covers zero/oversized/
64 MiB caches, eviction, aliases, mutable-reference bypass, close/reopen and
twenty ten-round lifecycle cases.

The installed-only source member skipped because its libraries are intentionally
bound only in the separate fresh-installed command; that command passed above.
The other skip is the explicitly enabled disk-backed cgroup fixture preflight.
It belongs to DIA-003's external service harness, not the standalone G1 utility
capability. The skip is retained and prevents treating this result as large
performance or service-boundary acceptance.

## Existing workflow compatibility

Command:

```text
Rscript --vanilla -e 'devtools::test(filter = "(dia-baseline-compatibility|prot-02b-qc-peptide-module-contracts|prot-02q-qc-peptide-orchestrator-shared|prot-05b-norm-module-contracts)")'
```

Result: 1,371 passed with zero failures, warnings or skips: 42 DIA
WorkflowState/real-S4 compatibility assertions, 772 peptide module assertions,
15 peptide-orchestrator assertions and 542 protein-normalization assertions.
These exercise the existing resident/default workflow. No GUI module is wired to
the saved store, default cache remains zero, and no other-omic default changed.

## Decision and limits

`G1:iq` is accepted only for this installed build and the standalone capabilities
above: immutable metadata lookup, complete RDS checkpoint save/hydrate, durable
acceptance/recovery, bounded optional residency and explicit resource release.
This permits the next standalone IQ experiments. It does not establish lower
memory, faster CPU/elapsed time, scientific agreement across treatments, GUI
completion, Limpa support or G2/G3/G4. Large performance acceptance remains held.

Code conformance, exactness and runtime completion are separate observations:
the implementation remains one internal flat-R utility with no export or new R
dependency; source/installed exactness checks passed; and the fresh installed
process completed. None is substituted for a performance result.

Rollback is ticket-scoped: preserve committed entity directories and `.pending`
recovery inputs, then remove only the saved-entity internal implementation and
tests in reverse ticket order. Existing public WorkflowState/S4/session files stay
readable because no public class, slot, export or file format was replaced. A
rollback must reopen and hydrate accepted records before disabling the internal
store if any later integration has begun; deleting the branch, records or failed
receipts is not a valid rollback.
