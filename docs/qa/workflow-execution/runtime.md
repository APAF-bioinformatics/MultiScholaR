# WFX-002 P14 runtime packaging

Design overview: [design.md](design.md); platform requirements: [preconditions.md](preconditions.md).

**Current: PERF-108's retained native runtime is accepted with an explicit primary
timing exception; fastduckdb 0.2.0 is published.** See the
[exact-artifact qualification](perf-108/closeout/qualification.json).
Historical WFX-002 unsigned/private-bridge holds below remain unchanged for those
older artifacts. GUI/default activation remains gated by PERF-107.
The [clean-install packaging repair](perf-108/packaging-repair/remediation.json)
is published and verified. The separate [native literal-list missing-result defect](perf-108/packaging-repair/native-literal-finding.md)
is now [repaired](perf-108/native-null-repair/remediation.json) in staged build
`browser-smoke.dWAQuA`. Compact list inputs are normalized with public C APIs;
ordered arithmetic remains unchanged. All 408 installed assertions pass. This
closes the correctness finding, while original performance measurements remain
bound to their original binaries. PERF-107 must qualify the corrected deployment
build through the actual GUI journey before default activation.
Owner: [WFX-002](../../../tickets/openai/WFX-002.md). Authority:
[execution contract](contract.md), [original audit](../wfx-ticket-audit-20260930.md)
and [remediation requirements](../wfx-ticket-remediation-20260930.md).

The current machine-readable closeout is `wfx-002-repair/wfx002-repair-closeout.json` under the existing
`MultiScholaR-dia-proof-20260911-VhC1n5` evidence root. It binds source files,
installed artifacts, dependency hashes, tests, failures and protocol exceptions.
No GUI default, public MultiScholaR R/S4 API or checkpoint/session format changed.

## Implementation and supported component

The [successor patch](../../../tools/profiling/fastduckdb-production-writer.patch)
applies directly to fastduckdb upstream
`d78511504fbdb89fbf34bc26b22b3ff1524b0555`, reproducing all 95 prepared files.
The original [experimental patch](../../../tools/profiling/fastduckdb-direct-writer.patch)
and accepted P14 artifacts remain unchanged. The successor identifies itself as
fastduckdb 0.1.0.9003, capability `incremental-parquet-v2`; it is an isolated Linux
x86_64 component, not a published project-wide dependency upgrade.

The incremental interface checks named equal-length primitive/factor columns,
rejects unsupported temporal/nested/matrix/classed inputs and counts expanded
UTF-8 payload bytes, including repeated strings, before native conversion.
It retains the 10,000-row / 64 MiB transfer limit, existing row-group flush policy,
encoder and eight encoder threads. Legacy `write_parquet_fast()` keeps its type
support. Constructor failure removes its private temporary; append/flush/close
errors abort owned unfinished output. Lifecycle functions and the capability query
are exported; generated FFI wrappers remain internal.

The retained FKrMqv Cargo.lock includes the experimental Arrow FFI graph absent
from the pristine upstream lock. Source installation uses `--locked --offline -j 8`
and records R/rustc/Cargo versions and source hashes. Locked cache and vendor
sources are preserved in `wfx-002/`; no project renv migration occurred.

The application installs artifact-manifest schema 3 with exact-sum source/header,
builder/configure identities, compiler flags, pinned DuckDB/duckplyr artifact hashes
and the installed fastduckdb capability, DLL and RDB/RDX/NAMESPACE hashes. The writer now checks this API
and build binding instead of experimental version equality alone. This internal
artifact manifest does not change saved workflow formats. Ordinary Arrow callers,
the exact numerical C kernel, three P14 batches and insertion-order setting are
unchanged. There is no bit-identical rebuild claim.

## Predecessor evidence (preserved history)

All members below are inside the evidence root's `wfx-002/` directory.

| Member | Outcome |
|---|---|
| `nsCiLa` | Clean native writer build passed; about 155 s, 4.54 GiB peak, zero swap |
| `8Mjpcm` | Clean application installation passed; about 69 s, 1.26 GiB peak, zero swap |
| `fUFZGc` | Installed preflight: 92 assertions, five exact P14 cases and 12 required figures; about 43 s, 1.15 GiB peak, zero swap |
| `nHENyY` | Projection 406/406 and resources 293/293 passed; known legacy harness failures preserved separately |
| `GjkaeY` | Final readback: 11 assertions passed; about 16 s, 1.66 GiB peak, zero swap; no scientific recomputation |

The exact cases cover ordinary, empty, missing/factor/annotation and both ordered
sum threshold fixtures. Installed checks cover non-rollup codec restoration,
projected labels, UTF-8, NA/NaN, attributes, partial/empty batches, expanded strings,
unsupported types, wrong/missing runtime, pending owners through GC, caller mutation,
constructor/append/flush/finish failures, abort/finalizer cleanup and publication
failure/retry. WFX qualification remains restricted to P14; P15–P19 were not reopened.

The predecessor XsaPvT package independently read the newly written saved table,
matched its current-reader hash and history, and recognized pending publication
on an owned copy. The accepted source head remained unchanged. Rollback can restore
the source files in `wfx-002/predecessor.json` and the retained XsaPvT library;
never reset the branch or delete accepted references/history/pending outputs.

The regression run retained three failed assertions and one error: two
`fixed_reclaim` capture expectations and two obsolete no-native-DLL expectations.
These match the documented unrelated baseline failures; they are not counted as
passes or repaired here. Public exports, class/generic/method compatibility and
legacy checkpoint checks executed before those stale DLL assertions.

## Preserved failures and authorized exceptions

The closeout retains the initial headroom refusals, upstream-lock resolution
failure, incomplete baseline-library admission and `e3cKLR` source snapshot refusal
caused by local compiler products. A runner copy excluded `.o`/`.so` files from the
clean source snapshot without deleting workspace files.

`dJmb0T` failed on a test-only namespace mock; its corrected preflight passed.
A test variable wrote the completed `fUFZGc` JSON to
`wfx-runtime/finish.parquet`. The original remains intact and hash-bound. Its
complete status and exact cases were independently verified in `GjkaeY`; the
variable was corrected without repeating scientific computation.

`9cRbAM` lacked DESCRIPTION in its isolated test input. `sZFOjL` then exposed the
test's incorrect metadata-only scope for row reads. Both failures remain; normal
bounded read scopes and cleanup order were corrected before the passing `GjkaeY`.

The user explicitly authorized startup below 16 GiB for this non-performance work.
The one-off component runner retains the 16 GiB hard cap, 15.5 GiB early stop,
zero swap, 2 GiB continuous host headroom and original deadlines. The regression
runner used its existing stricter 8 GiB cap and a five-minute deadline. The repository
runner and performance-qualification policy were not changed.

Two additional five-minute readback attempts were individually authorized after
the original recheck allowance: `sZFOjL` failed, `GjkaeY` passed. This history does
not grant further automatic retries. No full-size or timed performance member ran.

## Release eligibility

`runtime_release_status=held` remains explicit. Two compatibility exceptions
prevent production qualification and default activation:

- The fixed hash-verified exact-sum extension is loaded unsigned in a scoped
  connection; arbitrary external extension loading is not supported.
- DuckDB's private `rel_from_altrep_df` and `rel_to_view` bridge remains pinned to
  recorded DuckDB 1.5.5 and duckplyr 1.2.1 artifacts.

WFX-003/004 may consume this accepted component evidence. WFX-005 production
qualification and WFX-006 activation must reject these holds. Other platforms
remain unqualified. The next bounded work item is WFX-003; no more WFX-002 runs
are needed.

## Audit repair and executable clean-install recipe

The [WFX-002 audit](wfx-002-audit-20260930.md) remains immutable. Its four findings
are closed in the [repair record](wfx-002-remediation-20260930.md). Version 0.1.0.9003 advertises
`incremental-parquet-v2`, with **only `r_columns` qualified**. Arrow-stream input
is refused before import, including primitive streams; `owner`/`logical_bytes`
are retained signature parameters for compatibility, not trusted accounting inputs.
The old experimental Arrow interface can only be replayed with its retained build.
P14 continues to use fast_r. Incremental close uses same-filesystem atomic hard-link
creation (no replacement); legacy close retains replacing rename. Native abort
only removes its private temporary and cannot delete the final destination.
The application manifest is schema 3 and also verifies fastduckdb RDB/RDX/NAMESPACE.

The commands below reproduce the retained Linux source-install procedure in a
fresh owned directory. They do not authorize extra builds or performance runs.
Use the current source-bound patch and the pinned archive/cache from the evidence
root; installation does not download dependencies. The source request contains
hashes for every prepared file, checked by the existing runner before building.
The preserved runner copies have an explicit startup-headroom exception from this
session; all other hard/swap/continuous guards remain. A new qualification protocol
must explicitly decide whether it permits that exception.

```bash
set -euo pipefail
repo_root=/home/doktersmol/Projects/APAF/MultiScholaR
proof_root=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5
recipe_root=$(mktemp -d "$proof_root/wfx-runtime-install.XXXXXX")
mkdir "$recipe_root/candidate" "$recipe_root/bulk-rollup-library"
printf '%s  %s\n' \
  e2a13cc96c0da2389e8bbbbef17b29756a8a702823c828d88bd0e401f5fcbd68 \
  "$proof_root/wfx-002/wfx002-upstream.tar.gz" | sha256sum -c -
tar -xzf "$proof_root/wfx-002/wfx002-upstream.tar.gz" \
  --strip-components=1 -C "$recipe_root/candidate"
patch -p1 -d "$recipe_root/candidate" \
  < "$repo_root/tools/profiling/fastduckdb-production-writer.patch"
cp -a "$proof_root/wfx-002/bulk-rollup-library/." "$recipe_root/bulk-rollup-library/"
python3 - "$recipe_root" "$proof_root/wfx-002/cargo" <<'PY'
import pathlib, hashlib, json, sys
root = pathlib.Path(sys.argv[1])
source = root / 'candidate'
request = {'source': str(source), 'cargo_cache': sys.argv[2], 'files': {
    str(p.relative_to(source)): hashlib.sha256(p.read_bytes()).hexdigest()
    for p in source.rglob('*') if p.is_file()}}
(root / 'fast-writer-source.json').write_text(json.dumps(request, indent=2)+'\n')
PY
bash "$proof_root/wfx-002/authorized-build-runner.sh" --service reader-check \
  "$proof_root/wfx-002/browser-smoke.8Mjpcm/library" "$recipe_root" \
  "$repo_root" gui-bulk-rollup:writer-build
```

Record the exact `browser-smoke.<id>` emitted by that command and require its
`fast-writer-build.json` status and terminal exit to show success. Do not select
an arbitrary latest directory. The build helper copies the request's Cargo cache
into the package's private `src/.cargo` **before** invoking R CMD INSTALL. This is
necessary because Makevars overrides CARGO_HOME and uses `--locked --offline`.
The existing user Cargo cache alone does not satisfy this recipe.

```bash
# Set this to the exact successful member printed above.
writer_member="$recipe_root/browser-smoke.REPLACE_WITH_BUILD_ID"
test -f "$writer_member/fast-writer-build.json"
# The source-install runner requires a complete baseline application alongside
# the new writer; this is only a bootstrap copy, not the newly installed package.
cp -a "$proof_root/wfx-002/browser-smoke.8Mjpcm/library/MultiScholaR" \
  "$writer_member/library/"
bash "$proof_root/wfx-002/wfx002-source-install.sh" --service install-corrected \
  "$writer_member/library" "$recipe_root"
# Set this to the exact successful application member printed above.
application_member="$recipe_root/browser-smoke.REPLACE_WITH_APPLICATION_ID"
test -f "$application_member/library/MultiScholaR/duckdb/linux_amd64/manifest.rds"
bash "$proof_root/wfx-002/wfx002-component-runner.sh" --service reader-check \
  "$application_member/library" "$recipe_root" "$repo_root" \
  gui-bulk-rollup:writer-r-preflight
```

The source-install runner snapshots the current repository, excludes existing
compiler products, copies the new writer and isolated duckplyr/nanoarrow dependencies,
and installs MultiScholaR into a clean library. It binds the actual dependency
artifacts in that library. Other already-installed R dependencies must match the
source-bound environment recorded in the build receipts; missing dependencies
stop the procedure rather than auto-installing. Expected outputs are the native
build receipt, successful application terminal receipt, installed manifest schema 3
with complete writer identities, and `relational-preflight.json` plus
`reader-tests.rds` from the installed check. Private-API/unsigned-loading release
holds remain even when every command succeeds.

## Current repair decision

WFX002-A01–A04 are verified closed. The new build and checks are under
`wfx-002-repair/`: native `6i8DYV`, application `veas6R`, installed preflight
`28IVRZ`, and affected regressions/readback `a1r69d`. Preflight passed 109
assertions and five exact P14 cases; the second check passed projection 406,
resources 293 and predecessor readback 11. Readback also confirms preserved
history, recognized pending publication and an unchanged original checkpoint head.

The isolated second-check input initially failed to parse because an extracted
comment lost its marker (`UHreaj`). That input was corrected and all selected
files were parsed before the single permitted correction recheck. Its failure
remains preserved. No native build correction or performance member was needed.
Both builds and both successful checks had zero swap. The old closeout, audit,
source snapshot and original P14 evidence remain unchanged.

The lifecycle guarantee now distinguishes private temporary files from another
producer's final destination. Collision and competing-producer tests preserve
existing bytes; the legacy writer still replaces when requested. All Arrow-stream
inputs are outside the qualified incremental interface and fail before import.
An altered installed R artifact is rejected with the DLL/capability unchanged.
No production/default or broader workflow qualification follows from this repair.

## PERF-108 approved local successor

Current owner: [PERF-108](../../../tickets/openai/PERF-108.md). The accepted plan
retains fastduckdb's existing readers and `fast_r` writer, eight threads, three
complete-sample P14 batches and the exact source-ordered numerical kernel.
The application call sites now explicitly choose public SQL arm C; historical
experimental arm D retains its private bridge and identity. Results record
`aggregation_frontend` so old D evidence cannot be attributed to the new route.
Automatic selection and production release remain held.

The [local extension bundle](../../../tools/profiling/perf108-extension/README.md)
uses the official community C-API example's pinned template and build tools.
Its kernel is byte-identical to the accepted source. A local build is unsigned;
this is not a trusted-loading solution by itself. The package's existing held
manifest and scoped experimental loading remain intact until a compatible
signed artifact exists. No dependency upgrade, public R/S4 API or checkpoint
format change is included. No upstream upload has occurred.

Local build/check evidence and spent attempts are recorded under
`perf-108/` in the retained proof root and in
[the local protocol](perf-108/local-protocol.json). The full-size Stage 2 budget
is conditional on trusted artifact availability. The PERF-106 P14 exception is
not inherited by this runtime. PERF-107 still owns GUI action/default eligibility.

The user subsequently rejected the external community-publication direction.
PERF-108 now investigates a native query scope owned by fastduckdb using public
DuckDB C API functions from the pinned installed engine. The current fastduckdb
package is an Arrow/Parquet storage implementation and does not already own a
query connection. No private R connection pointer, unsigned extension or second
engine is authorized as a substitute. The initial native feasibility build hit
an R macro naming collision before executing any query; the narrow correction
is prepared, but the no-retry probe allowance is spent. See PERF-108's current
internal-route evidence; earlier successful SQL/P14 tests do not prove this new
native ownership interface. No external artifact has been published.

The authorized corrected feasibility member `wCua8m` passed; `4IVkLl` then
confirmed actual streaming over 100,001 rows and three independent scope cycles.
The [current decision](perf-108/fastduckdb-feasibility-decision.json) recommends
internal fastduckdb integration, with no community publication. The public
streaming call is deprecated but available in the pinned engine; isolate it and
requalify any future replacement. This native feasibility result does not yet
qualify installed fastduckdb handles, Arrow/R conversion, application outputs or
large-data performance. PERF-108 remains in progress, not release eligible.

The subsequent [PERF-108 cache qualification](perf-108/cache-qualification/outcome.md)
supersedes that feasibility-only status. Installed checks and all 18 timed
members passed exact science and measurement validation. Candidate memory passed
every paired comparison (approximately 7.6 GiB rollup, 2.9 GiB primary filter,
3.5 GiB wide filter). Rollup and wide-filter timing passed; primary-filter timing
failed the fixed 1.05 paired limit. Production release remains held, with no
push/tag/release or GUI default activation. Historical failures and the explicit
setup-recheck authorization remain bound in the owning evidence location.

The [structural experiment](perf-108/structural/outcome.md) is also closed without
promotion. Paired aggregation did not meet the gain threshold; per-sample
validation was substantially slower; the persistent reader retains an unresolved
application cleanup defect. Owned native output passed exact science, but full
300-sample medians were 7.39 GiB / 106.45s versus 7.64 GiB / 99.21s for valid
accepted controls. One separate current control failed the frozen monitoring-gap
rule. All four small and six full rollup outputs verified; the gate prevented
consumer confirmations. Retain P14 and the original release hold. No defaults,
dependency versions, checkpoint formats or release status were changed.

The [ordered-group/reader follow-up](perf-108/ordered-groups/outcome.md) repaired
and verified the reader cleanup defect. Reader component gains remained below
the advancement threshold. Ordered reduction preserved exact science, but on
300 samples used median 6.58 GiB / 138.99s versus 7.08 GiB / 87.24s for matched
current P14: about 7% less RAM and 59% more time. All ten complete operations
verified with valid monitoring and zero swap. The decision is
`no_justified_rollup_gain`; consumer confirmations were gated off. Keep current
P14 and the correctness repair; experimental routes and release remain disabled.

The [retained-runtime closeout](perf-108/closeout/outcome.md) qualifies the exact
installed P14 route scientifically and preserves all required outputs. Its reused
full rollup pairs pass the original peak/time gate (median 7.08 GiB / 87.24s).
All eight new primary/wide filter members passed independent science, memory,
monitoring and cleanup; installed API/S4/history/rollback checks and the source
package check also passed. One primary-filter pair is 5.91% slower (47.51s versus
44.86s), exceeding the fixed 5% allowance by 0.407s. The original numerical decision
remains held; the user explicitly accepted the bound narrow timing exception.
The retained runtime is accepted through `closeout/qualification.json`, and
fastduckdb 0.2.0 is published from commit `942c56caf17a45eeb41b81984bed010b5ed7a09a`.
No retry or threshold change occurred.
The [inherited-clause reconciliation](perf-108/closeout/acceptance.json) keeps
science/compatibility separate from that release decision. GUI/default rollout
remains PERF-107.
