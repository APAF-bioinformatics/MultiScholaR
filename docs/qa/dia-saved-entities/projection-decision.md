# DIA selective-table projection decision

## Revision 1.6 — installed adapter decision, 2026-09-19

Frozen before performance execution: seven independent blocks for each of three
single-action cases, with native RDS, current store DuckDB and store fastduckdb
members in alternating forward/reverse order. This is 63 timed members, each
with independent output verification, plus one preparation. No repeats beyond
the frozen seven-block ceiling. Keep 8 GiB/zero swap/900 s, >=10 GiB launch
headroom, 7.5 GiB job/<2 GiB host early stops, physical storage and serial lease.
No profiling or forced filesystem-cache conditioning; fresh process does not
mean cold OS cache. No proof job was live at entry.

| Case / unchanged source shape | Request | Actual consumption |
| --- | --- | --- |
| Small DA, 8,020 x 57 | All 8,020 rows; Protein.Group/logFC/raw_pvalue/fdr_qvalue/comparison/uniprot_acc | Existing volcano filter and result-page helpers, full adjusted population |
| Wide DA, 10,000 x 105 | All 10,000 rows; same six columns | Same helpers, no estimator refit or population truncation |
| Recorded precursor, 960,000 x 14 | Rows 480001..489600; Run/Precursor.Id/Precursor.Quantity | Materialized typed range consumed and exported; prospective cold-table range use, not an existing GUI speed claim |

Reuse installed oSQQvG (`9f6077d53e3c710e30e3c7db2163d64c5f43bbe286db2866a904231f1170e4f2`
RDB) and recorded SxedsH source (original session
`4671bdfab031cd880ab2bd2b414c4625cc30902a247e7d60f1c8dc11411a50c3`). The application
adapter/optional native binary and numeric protocol are unchanged. Preparation
uses the real metadata-v2 store writer and native RDS writer, measuring each
write/validation/persistence cost. Each store action opens the real package
scope, validates once, restores through the package, and closes its lease before
output delivery. fastduckdb footer clearing remains per-read. No native bypass,
increased collection limit or unsupported ranked/key/full-precursor fallback.

Every member includes installed startup, input integrity/decode, the named
consumer, required exact-result file/fingerprint, cleanup and exit. The native
RDS comparator decodes once, uses direct column selection or bounded native row
indices, and drops its owned full table before output delivery. Its post-decode
selection/consumption timing is a separate resident-action reference—not evidence
that an already-resident GUI pays the cold-file cost. Do not manufacture five
uses or infer that previews are repeated scientific QC generations.

Decision rule: component advancement only if point and 95% one-sided paired
log-ratio upper elapsed/CPU ratios are <=1 and peak <1 against **both** cold
comparators. Point failure rejects; uncertainty-only failure holds, with no
adoption or favorable reruns. Report source/write costs and producer-plus-consumer
point estimates separately; a single producer observation cannot qualify a full
workflow confidence claim. Unsupported consumers and the history gates stay held.

Initial small setup check: 892 assertions passed, one exactness failure exposed
the native column-subset control dropping provenance. The sole bounded setup
correction restores original table attributes in that control/oracle; values,
rows, classes and exact comparisons are not relaxed. No performance run had
started. Historical failures below are retained.

Setup correction passed 240 projection/lifecycle assertions with no failures,
warnings or skips; the unaffected harness selection checks passed 666 assertions
in the preceding run. Preparation `browser-smoke.uuB3ie` completed at
1,550,540,800 B peak and zero swap. Its manifest is
`caa34a2ce01607c378631010e0eda0ee3c11c866366cd1188d8d426eac7883b3`.
Pre-measurement driver/runner/shell SHA-256 bindings are
`31bbbe644bb5044e989607552cd08d3c288071b8a38979a10e3d3664718abd96`,
`696433086d934c8f4d15251f88bd948d5cc70e56a26592cb08e8cb0074f6da79`,
`c3a137bcdcc2f7936df7107a784697dcfedc2ca1ced0cb2e62b37b87470da30b`.
Each output file is fingerprinted and independently compared as an exact R
object against the original source, including attributes. Attribute serialization
order is not substituted for scientific/object equality. All seven blocks are
declared before reading any performance outcome; they are not retries of a
screening result. The existing reader-screen receipt durably records every
member/verifier and block, including partial completion or failure.

### Final decision: fastduckdb not adopted for these consumers

All **63 timed members and 63 independent verifiers** completed. Every result
matched the original R values, classes, attributes and requested population;
owned resources were released and every swap peak was zero. Maximum timed-member
peak was 1,081,032,704 B (1.007 GiB). No benchmark failed, was replaced or exceeded
the fixed budget. No proof services remain running.

Complete receipt:
`/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/reader-screen.2d29555578cb1d/screen.json`,
SHA-256 `856e956b2df3158e8c31cff40a4bd2d61e456e3a30b58d22615602de28cbbcd9`.
It binds all blocks, member/verifier directories, invocation identities, inputs,
outputs and resource/timing metrics. Runtime source hashes above were unchanged.

Ratios are paired geometric means with one-sided 95% upper bounds in parentheses;
below 1 favors fastduckdb. These include complete cold-read service boundaries.

| Case | Comparator | Elapsed ratio (upper) | CPU ratio (upper) | Peak ratio (upper) | Decision |
| --- | --- | ---: | ---: | ---: | --- |
| Small DA | RDS | 1.0520 (1.1688) | 1.0526 (1.1647) | 1.0006 (1.0013) | Not adopted |
| Small DA | Store DuckDB | 0.9754 (1.0061) | 0.9760 (1.0048) | 1.0007 (1.0012) | Not adopted |
| Wide DA | RDS | 1.0137 (1.0231) | 1.0148 (1.0234) | 0.9967 (0.9971) | Not adopted |
| Wide DA | Store DuckDB | 0.9765 (0.9846) | 0.9770 (0.9839) | 1.0006 (1.0010) | Not adopted |
| Precursor range | RDS | 0.7999 (0.8067) | 0.7987 (0.8046) | 0.9016 (0.9023) | Read-only benefit vs RDS |
| Precursor range | Store DuckDB | 0.9805 (0.9831) | 0.9802 (0.9821) | 1.0048 (1.0055) | Not adopted |

Keep resident/RDS DA access: neither DA population benefits overall. For large
already-cold ranges, avoiding whole-RDS decoding helps, but the optional reader
does not improve both time and peak against the existing DuckDB reader. Its
~1.95% elapsed advantage has ~0.48% higher peak. No post-result tolerance is added.

A **secondary comparison of the controls** shows existing DuckDB versus RDS at
elapsed 0.8158 (upper 0.8219), CPU 0.8148 (0.8204) and peak 0.8973 (0.8977) for
the precursor range. This supports narrow cold selective access, not a new GUI,
persistence-format or general workflow qualification.

Creation cost matters: precursor RDS persistence took 2.352 s and committed
Parquet 3.737 s. Median read-service times were 7.701 s RDS, 6.297 s DuckDB and
6.162 s fastduckdb. Adding that single measured write gives descriptive totals
10.053/10.034/9.899 s: most of the first-use elapsed benefit disappears. These
sums are not measured composed workflows or producer confidence intervals.
Already-resident selection/consumption in the RDS arms took median 0.004 s for
DA and 0.003 s for the range, excluding decode. Those are different boundaries,
not a valid whole-service speedup ratio; flushing hot tables is not justified.

**DIA-013 optional-reader evaluation is complete, not adopted for default/GUI
use.** Keep the tested capability opt-in and preserve metadata-v2 readability,
but do not activate it or repeat unchanged cases seeking a pass. The historical
native five-action win remains bound to its older experiment. History and
G2/G3/G4 gates are not promoted; the complete DIA GUI target remains unproven.
The architectural recommendation is to preserve hot working data, retain
resident/RDS DA results, and consider existing bounded Parquet/DuckDB reads only
for genuinely cold large tables. Further work should target measured workflow
computation/copy/output reuse, not another checkpoint-codec or reader search.

## Revision 1.5 — package bridge, 2026-09-18

Current bounded item: implement the existing table-store bridge and prove small
exactness/resource/recovery plus fresh installed compatibility. Resident history,
GUI defaults and all scientific methods remain unchanged. This is not a new
reader search or performance acceptance. The historical native-I/O results below
remain bound to their original inputs/builds/boundaries.

The optional, explicitly selected `fastduckdb` reader is declared in Suggests
and pinned in Remotes to the already-tested commit
`d78511504fbdb89fbf34bc26b22b3ff1524b0555`. It uses two threads and eager strings;
it creates no DuckDB connection. It exposes ordered column projection and
contiguous source-row ranges, not key joins, ranked heatmaps, arbitrary SQL,
full-object hydration or scientific recomputation. Default readers are unchanged.

Restoration metadata version 2 is explicit and opt-in. Version 1 remains the
writer default and stays readable. Version 2 admits plain data frames and
ordinary tibbles containing plain atomic/factor/ordered-factor columns; untested
special column classes fail explicitly. Table/column attributes remain bounded;
explicit row names and element names are ordinary hidden Parquet columns, not
row-sized metadata. The 8 MiB restoration, registry, 64 MiB payload-cache and
existing collection limits are not increased. DuckDB can also read version 2,
so disabling the optional reader does not strand the new payloads. Payload and
descriptor codec identities, Arrow/Snappy writer, publication/recovery and
public R/S4 exports remain unchanged.

The pinned native reader caches 64 footers by entry count, without a byte cap.
The application adapter clears that process-wide metadata cache after each native
read, including failures. R restoration metadata/integrity remain action-scoped
and are reused once per action; no table payload cache is added. This eviction
can affect other native-reader cache hits in the same process, but not their
files or results. It is an explicit resource/performance tradeoff. The old
scope-cached five-action performance result cannot qualify this adapter unchanged.

Initial focused verification had two new test-order setup failures (helpers were
defined later); 372 existing assertions passed. After fixing the setup, 457
projection/resource assertions passed without warnings/skips. Expanded checkpoint,
export, projection, recovery and resource checks then passed 726 assertions;
the original failure is preserved here rather than erased. Fresh installed checks
also cover the final explicit rejection of untested special column classes.

Next verification is one build and one installed reader-check through the existing
guarded runner: 8 GiB, zero swap, 900 seconds each, >=10 GiB launch headroom and
the unchanged early stops. No new Rust build: copy only the already-pinned optional
reader (binary SHA-256 `3989e6fdc620e90cc1d04e9c44d373d2cd3ad564ae8aeea7ae497a29c40999dd`)
from `/tmp/multischolar-dia013-fastduckdb.NxbxFU/library`. Its unchanged application
is the installed compatibility comparator. The installer includes the four
existing untracked R files already named in Collate; omitting them would not
test this working build. No long performance member runs in this bounded item.

Pre-edit SHA-256 bindings: store
`fe689bb76be748a53fac89b8d66048bcf17a79860ab2e7e1ac39b3baf7a64624`, codec
`f494d75fbbc2aa7975b9e79e3327f6818b26743dc0b67ed7166ee172936c001a`, DESCRIPTION
`400099e5a2f9e86d2cf58b7d850d7c7822853a5c35973901ea1c3318c6091e29`.
NAMESPACE remains `9dc606a759229234e12f5a0b55879af5b7230a42f1177a4d86fa1e088879ccc5`.
The established package uses DESCRIPTION/Remotes (renv.lock is ignored); this
bridge does not introduce an environment-manager migration. New reader calls are
qualified and lazy; no eager optional namespace import is added.

### Installed bridge closeout

Build `browser-smoke.oSQQvG` completed in 49.630157 s, peak 1,128,108,032 B,
zero swap. The candidate library is
`/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/browser-smoke.oSQQvG/library`.
The existing `reader-check` service `browser-smoke.CJihRE` then passed **1,149
assertions across 73 tests**, with zero failures/warnings/errors/skips, in
130.490595 s. Peak was 1,949,036,544 B and swap zero. It checked the new installed
adapter, the unchanged native trial, action/resource contracts, and fresh-process
public exports, loaded namespaces, S4 class/generic dispatch, R6 workflow behavior
and legacy RDS round-trip against the bound predecessor. All services stopped.

These are implementation/compatibility checks, not performance members. The
next required DIA-013 item is an explicitly bound adapter comparison for intended
requests within the unchanged limits (for example the recorded preview/range
requests), including validation, restoration, creation and cleanup costs. The
old full-precursor native-I/O bundle bypassed store collection ceilings and
cannot be relabelled as this adapter's performance. No benchmark or GUI was run
in this implementation item; DIA-013 remains in progress for that decision.

Installed RDB/RDX hashes:
`9f6077d53e3c710e30e3c7db2163d64c5f43bbe286db2866a904231f1170e4f2` /
`6e60c5d25b68bca7ec93d612c40bdbf0b616aa71f86aa60c14be743e7b6f71b0`.
NAMESPACE and the copied native reader retain their predecessor hashes above.
`CJihRE/installed-compatibility.json` SHA-256 is
`6c80491ee5a452602808214edffedf17d76da7713688a5dc5b57fdd1d47bfcc4`;
`reader-tests.rds` is `efc952c1b229476e0c12825907df2e712390a9099664dac2d913b668e2811fe0`;
`reader-operation.json` is
`294299a228ce3295f9303c16408bd42097670f06739291615e277ba08cad57cc`.
Exact commands, source patch/untracked-file bindings, library origin and service
accounting are retained in the request-owned directories under the proof root.

Final source hashes: store
`cc26c074b0d86a342f4b76a91de39b9f077d47265048d0af356087605417eafd`, codec
`65993a267dcb9af31db2f2d439bce7220fdbc085d19546576c15a6de38a7ba07`, DESCRIPTION
`3a959832ee7e47825e4ecad6432e61d844affbac10c3a9d4f0a52a67e286a4d2`, driver
`8e2ad0376d154e8a3279ad4c6b6054f3cded6a96b346dcc529a2a98b7b3e2d46`, shell
`fd3fed566f11d99acf39678216210f136324ea0d8054cd47b0f6ac13c510dd88`, projection tests
`861fed59ea3a1eea1574c4391a167e5597f8c4efa4c5188764eb6cb42fe5ed45`.

Rollback policy: stop selecting `reader = "fastduckdb"` and keep default DuckDB
reads; stop opting new writes into metadata version 2 if desired. **Keep the
version-2 decoder** for any already-written v2 tables. Do not remove those readers,
reset the dirty worktree, rewrite accepted payloads or erase failed receipts.
The optional native library can remain isolated; no global dependency was installed.

Historical decision follows unchanged.

Current decision: **DIA-013 revision 1.4 experiment complete**. The repeated
precursor-reader bundle qualifies to advance to DIA-015; one-off precursor access
is held on timing uncertainty, and small/wide bundles are not adopted. No GUI or
G2/G3/G4 gate is promoted. Revisions 1.2 and 1.3 remain
completed, not-adopted experiments. Their
results below are not superseded by a new dependency build or correctness pass.

## Revision 1.4 — frozen reader experiment, 2026-09-18

Authority: user-approved implementation plan, DIA-013 revision 1.4. No application
code/dependencies/defaults, checkpoints, BLOBs, writer replacement or GUI changes.
One pinned fastduckdb source: d78511504fbdb89fbf34bc26b22b3ff1524b0555, archive
SHA-256 e2a13cc96c0da2389e8bbbbef17b29756a8a702823c828d88bd0e401f5fcbd68,
Cargo.lock SHA-256 0dfffbc58e67f7c3eff793b6b7fdbe17ebe8733838e52297aebeef02dc02e073.
Isolated library: /tmp/multischolar-dia013-fastduckdb.NxbxFU/library. Its copied
MultiScholaR build is byte-identical to the earlier qs2-pilot installation;
checkpoint results are not rerun or reinterpreted by this experiment.

Build pWCnQ5 succeeded in 450.119854 seconds, peak 3,005,743,104 bytes, zero
swap, two Rust build jobs. Native reader binary SHA-256:
3989e6fdc620e90cc1d04e9c44d373d2cd3ad564ae8aeea7ae497a29c40999dd.

Fixed shapes: existing seeded 8020 x 57 and 10000 x 105 DA tables, and only the
960000 x 14 raw_data_s4 table from recorded preparation SxedsH. The latter is
native table I/O, not a new saved-table codec or proof that production metadata/
collection ceilings admit a complete precursor table. Restore its compact native
row names in a trial-only sidecar; do not change application metadata limits.

Each member is a fresh installed-package process, one reader only (RDS, Arrow,
DuckDB or fastduckdb). Parquet readers use identical Arrow-written Snappy files.
Two requested decoder threads, eager strings, zero optional payload cache.
One action reads full, projected, and range/page data; the RDS comparator decodes
once and reuses that value within the action. DA projection consumes the existing
full-population volcano and sorted-page helpers, with thresholds 0.05 and 1;
precursor projection consumes Run/Precursor.Id/Precursor.Quantity. Pages select
rows 1..25; the precursor range selects 480001..489600. Whole-table comparators
are explicitly not bounded store reads. Unsupported arbitrary-key/ranked disk
selection is not emulated. Exact first-action outputs are published as uncompressed
RDS and checked in a separate service; output work remains in measured completion.

Measure one action and five successive actions separately. Reader/action timing,
producer write+fingerprint time and service startup-to-exit elapsed/CPU/combined
peak are separate fields. Independent preparation is not folded into read time;
producer-plus-consumer totals explicitly include the corresponding persistence
cost. Do not sum nested read/action intervals. No profiler or host cache-dropping;
fresh process does not imply cold OS cache. File-cache charges remain in peak.

Screening member order is RDS, Arrow, DuckDB, fastduckdb for the one-action case,
reversed for five actions, serially for small/wide/precursor. A valid unfavorable
screen is retained, not rerun. Only promising exact cases advance to seven blocks,
alternating paired order and using the frozen paired log-ratio/95% upper-bound
method. Elapsed and CPU must not regress against RDS/current access; memory must
decrease. Report direct DuckDB/Arrow mechanism comparisons separately. No pooled
shapes, no finite speedup denominator from a failed/OOM control, no G2/G3/G4 claim.
One bounded correction and one independently justified invalid-run replacement
are the entire allowance, including instrumentation failures. Guards remain
8 GiB/zero swap/900 seconds, launch headroom 10 GiB, early stop 7.5 GiB or host
availability below 2 GiB. The existing runner owns every service and serial lease.

Historical results follow unchanged.

Revision-1.4 setup evidence: ubcv8G passed all four readers' small exactness and
965 focused harness/projection/resource/installed-compatibility assertions, zero
failures/warnings. Peak 2,013,728,768 bytes, zero swap. Preparation dH9ROw then
failed closed at the recorded precursor's non-plain container class; its partial
files and failure receipt remain preserved (peak 1,299,943,424 bytes, zero swap).
The single bounded correction is allocated to preserving that native container,
not relaxing the application's plain-data-frame codec. Guarded schema diagnostic
GCRXyT confirmed `tbl_df/tbl/data.frame`, only normal table attributes, and a named
Run vector. In the trial-only native I/O case, container class/compact row names
remain bounded metadata and element names travel as ordinary Parquet columns.
Exact tibble/named-vector small checks are added before replacement preparation.
This is not a change to `R/utils_workflow_saved_table_codec.R` or GUI storage.

The corrected check t5UFIe failed before any performance member: the Arrow arm
did not preserve per-element column names in the tibble fixture. It also emitted
two deprecated-tibble-row-name warnings. Peak 998,215,680 bytes, zero swap,
9.049626 seconds main-process interval; its failure receipt is preserved. A
read-only diagnostic over the retained 300-row fixture, explicitly loading tibble
before comparison, isolates the difference: original and RDS result have 300
element names on Protein.Group, Arrow result has zero. Top-level provenance is
preserved in that diagnostic. This is an unresolved trial adapter/round-trip
defect, not evidence that fastduckdb itself is slower or scientifically incorrect.

The one bounded correction allowance has now been exercised and verification
still fails. No replacement preparation, screening member, seven-block series or
GUI activation is authorized by this result. The screen command is implemented
but has not been qualified by executing a complete member/verifier sequence.
Producer-plus-consumer summaries likewise remain unexecuted; do not claim those
interfaces are accepted evidence. No dependency, limit or scientific contract
was changed to manufacture a pass. All owned services are stopped.

Next bounded action, requiring a further correction allowance: repair per-element
name restoration across collected data-frame/tibble containers, preserve source
row-name semantics without deprecated assignments, rerun the new tiny exactness
checks and affected regressions, then retry the frozen preparation only after
they pass. Do not add more readers, upstream changes or workflow refactors.

That additional bounded correction was authorized by user "please continue" on
2026-09-18. Base-R reproduction isolates the cause: `names(df$x) <- labels` invokes
data.frame column replacement and strips the names, whereas assignment on an
unclassed column list preserves them. The trial now restores names on that list,
sets row-name attributes without deprecated tibble assignment, then reinstates
the original class. Previous failures above remain unchanged; no timing evidence
is invalidated because screening had not begun.

Correction check Lbl5gS passed 966 assertions with no failures, warnings or skips,
including exact tibble/vector names across all four readers and fresh installed
compatibility. Service elapsed 71.292429 seconds, peak 1,952,350,208 bytes, zero
swap. Corrected driver SHA-256:
6f8ea4ccbf016a3ce260fbffe99f49c62a0f02448a0bab6872f9794c1d450917.
A tiny complete member/verifier/late-consumer-failure test is also required before
screening, because reader equality alone does not validate its completion receipts.

That lifecycle test and all affected checks passed in 6NuFoK: 976 assertions,
zero failures/warnings/skips, complete RDS and fastduckdb members independently
verified, and an injected late consumer error left no success receipt or open
reader. Service elapsed 145.512455 seconds, peak 1,949,601,792 bytes, zero swap.
The reader build is unchanged; its previous successful build is reused.

Replacement preparation u4YWEo succeeded: all three frozen shapes, service
28.966763 seconds, peak 1,368,440,832 bytes, zero swap. Input manifest SHA-256:
097fceb1060615a5325f51f95081026022865536cc37903b14921d7dfb6fa7dd.
Screening uses this immutable preparation only. Source extraction and independent
correctness preparation are not represented as reader speed improvements.

Reproduction uses the existing approved runner banner:

```bash
bash tools/profiling/run_workflow_checkpoint_proof.sh --service reader-check \
  /tmp/multischolar-dia013-fastduckdb.NxbxFU/library \
  /home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5 \
  /home/doktersmol/Projects/APAF/MultiScholaR
```

At the historical correctness hold, driver SHA-256 was
4bf95c353cb346f99b985aaed48d466dea00c9c6576123b566b832c107bcf127.
The original passing driver was
82137ae5521c847f18509606346c8d6b9e7aef8046a970e3965cca3f3d641be6;
its 965 assertions do not qualify the subsequent tibble changes. Installed
MultiScholaR R database/index hashes remain those in checkpoint-results.md's
matched qs2 pilot. The fastduckdb installation is isolated and never imported by
the application; removing/disabling reader trial calls leaves original RDS,
Parquet and qs2 files/readers untouched. Do not delete evidence or reset the
already-dirty worktree when rolling back this profiling-only experiment.

### Revision-1.4 screening outcome

All 24 members and their independent verifiers completed under the guard, with
zero job swap and exact first-action tables/attributes plus all consumption
summaries. Screening receipt: reader-screen.1edd0814f7260f/screen.json under the
physical proof root. These are single pairs, not uncertainty-qualified claims.

| Bundle | Actions | RDS service s | fastduckdb service s | RDS peak bytes | fastduckdb peak bytes | Decision |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| 8020 x 57 | 1 | 17.674727 | 10.898343 | 992100352 | 1003839488 | not adopted; peak higher |
| 8020 x 57 | 5 | 12.337824 | 19.245021 | 991363072 | 1003859968 | not adopted; time/CPU/peak higher |
| 10000 x 105 | 1 | 12.872066 | 17.786046 | 1007169536 | 1024475136 | not adopted; time/CPU/peak higher |
| 10000 x 105 | 5 | 20.464053 | 11.265306 | 1008001024 | 1019953152 | not adopted; peak higher |
| 960000 x 14 | 1 | 26.309608 | 12.238108 | 1665359872 | 1483075584 | eligible for seven pairs |
| 960000 x 14 | 5 | 55.199811 | 21.690337 | 1666039808 | 1504108544 | eligible for seven pairs |

The measured bundle includes full-table, projection and contiguous-range access;
its peak is not a separate projected-read-only peak. Service elapsed minus the
action intervals includes namespace startup, GC, receipt/exit work and other
unattributed costs; startup has not been isolated from that remainder. Full
RDS decoding occurs at action open, so its inner selection timings alone cannot
be compared with a Parquet decoder's read timings. The screen retains all four
engines' individual events/CPU/peak metrics; do not pool shapes or select favorable
startup variation to claim a universal speedup.

The producer-plus-consumer arithmetic adds the measured format write/fingerprint
cost to service completion. Common fixture/metadata construction remains in the
separate preparation receipt: this arithmetic is not a measured composed GUI
pipeline. Cgroup peak includes native allocations, overlapping input/output and
serialization buffers, kernel/file-cache charges; allocation ownership and these
contributors were not separately profiled. Do not sum object sizes or process RSS.

Only the two eligible precursor bundles proceed to six additional independent
pairs, reusing their matching screen pair as block 1. Member code/build/inputs/
settings/boundaries remain unchanged; order alternates by block. The coordinator
calls the existing frozen runner and writes a recoverable block map at
reader-qualification.1ff3914f43aa20/qualification.json (orchestration source is
preserved alongside it). Small/wide failures are not retried. No G2/G3/G4 or
GUI promotion follows from this component experiment.

### Revision-1.4 seven-pair outcome and handoff

Both precursor cases completed seven independent matched pairs, with alternating
order and unchanged input/build/member code. There were 48 unique performance
members overall (24 screening + 24 additional precursor members), each with a
successful independent verifier. Every checked first-action output preserves
values, ordering, attributes, class and element names; every action's consumption
summary agrees. Open reader resources and job swap are zero. Maximum measured
performance-member peak is 1,780,273,152 bytes; maximum elapsed is 55.199811 s.
All proof services are stopped. No additional favorable reruns were performed.

Ratios below are fastduckdb / whole-table RDS. Points are geometric means of the
seven within-pair ratios; upper bounds are the frozen one-sided 95% Student-t
bounds on paired log ratios (df = 6). The protocol's approximate-normal-log-ratio
assumption remains an assumption, not a claim proved by seven observations.

| Case | Metric | Point ratio | Upper bound | Interpretation |
| --- | --- | ---: | ---: | --- |
| One action | Service elapsed | 0.872857 | 1.125679 | held; cannot rule out regression |
| One action | Service CPU | 0.879592 | 1.122535 | held; cannot rule out regression |
| One action | Service combined peak | 0.894626 | 0.902608 | lower peak |
| Five actions | Service elapsed | 0.660676 | 0.851433 | 33.9% lower point estimate |
| Five actions | Service CPU | 0.670779 | 0.855590 | 32.9% lower point estimate |
| Five actions | Service combined peak | 0.897689 | 0.900645 | 10.2% lower point estimate |

Adding the separately measured format write/fingerprint cost gives elapsed ratios
0.839098 (upper 1.059173) for one action and 0.652787 (upper 0.830947) for five.
CPU equivalents are 0.847400 (upper 1.059565) and 0.663503 (upper 0.836403).
These remain arithmetic producer/read accounting, not a measured composed
workflow. Parquet is larger here: 29,289,907 bytes versus 22,501,089 bytes for RDS.

**Decision:** the five-action case meets this trial's no-elapsed/no-CPU-regression
and lower-peak component advancement rule, including its uncertainty bounds.
Only that recorded, repeated precursor-reader bundle advances to DIA-015.
One-off access stays held and small/wide tables retain RDS/resident policy.
Do not infer an automatic size cutoff, projected-read-only peak reduction,
resident-GUI speedup, or another omics/normalization route's qualification.

The full GUI target remains unproven. Even within this component experiment,
the five-action peak upper ratio 0.900645 does **not** meet a strict <=0.90
confidence-bound threshold. Its point reduction of 10.2% must not be promoted
as proof of the final 10%-elapsed-and-10%-peak GUI gate. DIA-015 still needs its
history prerequisites and a matched composed producer/consumer/output boundary.
No application code, dependencies, defaults, S4 APIs or checkpoint codecs changed.

Raw paired observations (RDS then fastduckdb):

| Actions | Block | Elapsed s, RDS | Elapsed s, fast | CPU s, RDS | CPU s, fast | Peak bytes, RDS | Peak bytes, fast |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1 | 26.309608 | 12.238108 | 25.132393 | 12.118392 | 1665359872 | 1483075584 |
| 1 | 2 | 31.242964 | 26.371180 | 29.442008 | 24.986686 | 1665622016 | 1483403264 |
| 1 | 3 | 17.608796 | 14.349396 | 16.892462 | 13.959138 | 1641230336 | 1509064704 |
| 1 | 4 | 19.131570 | 24.666808 | 18.360772 | 23.540626 | 1665343488 | 1483796480 |
| 1 | 5 | 14.692745 | 19.074743 | 14.290125 | 18.395798 | 1665470464 | 1483718656 |
| 1 | 6 | 30.649958 | 27.535584 | 28.957371 | 26.170484 | 1664856064 | 1483214848 |
| 1 | 7 | 30.131941 | 24.174879 | 28.507770 | 23.022517 | 1667457024 | 1483046912 |
| 5 | 1 | 55.199811 | 21.690337 | 52.072342 | 21.336100 | 1666039808 | 1504108544 |
| 5 | 2 | 39.824082 | 20.447021 | 38.046614 | 20.023251 | 1665228800 | 1489252352 |
| 5 | 3 | 26.282183 | 24.367995 | 25.554854 | 23.777637 | 1665265664 | 1489580032 |
| 5 | 4 | 44.425621 | 22.587981 | 42.462692 | 21.969078 | 1665466368 | 1488781312 |
| 5 | 5 | 39.567411 | 33.979033 | 37.948734 | 32.798083 | 1665081344 | 1494421504 |
| 5 | 6 | 53.633834 | 38.439881 | 50.677101 | 36.836401 | 1665073152 | 1504555008 |
| 5 | 7 | 40.873843 | 38.364769 | 39.232944 | 36.759032 | 1665597440 | 1494433792 |

The recoverable qualification receipt contains every member/verifier directory,
unit identity, metric and pairing. SHA-256 bindings:

- qualification.json: d8caba96bd256095ab6b1ea85e9717d0431777d5263d6ef11a74710e8c5b11ac
- screen.json: 662ec5e005bd1e54959ea31c850cf7dd91b8fdbb076a1cbbc742d9b290ac12a0
- installed-member driver: 6f8ea4ccbf016a3ce260fbffe99f49c62a0f02448a0bab6872f9794c1d450917
- R runner: beebb869fd9b8c0bec79d7fd6bbed11cb25f14e7bbc4af32504bb50c38cdad40
- shell runner: f8e72ee49397efb71ed483fc477750c1b92dc6fdadb2bb51918c307489d1c75b

All 48 member receipts match the same input manifest, driver and native reader
binary. The earlier failed setup runs remain available and are excluded from
performance denominators. Disabling the profiling treatment leaves original
readers and all saved application files usable; the private dependency remains
only in the isolated proof library. No cleanup of evidence or branch reset is
part of this closeout.

## Selected treatment

Use one internal codec: Arrow writes a single Parquet payload with bounded
serialized R type metadata in the Parquet schema; DuckDB performs parameterized
column/row selection directly against that owned file. DBI owns the connection
for one request and closes it before returning. There is no arbitrary SQL API,
long-lived database, second table backend or checkpoint/S4 migration.

Installed versions at selection are Arrow 24.0.0 and DuckDB 1.5.5. Arrow is
already a package import. DIA-013 adds DBI and duckdb as explicit package
dependencies for qualified internal calls. Existing complete RDS remains the
whole-object comparator and mandatory compatibility codec.

## Named consumers and populations

The qualified reduced DIA IQ journey produced 8,020 adjusted protein rows, 24
ordered samples and a 57-column DA export. The selective requests are frozen as:

| Consumer | Statistical population | Requested disk data | Avoidance target |
| --- | --- | --- | --- |
| Results page | All already adjusted rows; no page-specific test/FDR | One ordered 25-row page and displayed columns | Avoid collecting the other 7,995 rows and unused sample columns |
| Volcano | All already adjusted rows | Protein key, contrast, log2FC, raw p-value, FDR and requested annotation | Avoid the raw/normalized assay columns while retaining all test rows |
| Heatmap selector | All already adjusted rows select top 50; complete tested matrix remains the authority | Selection columns plus top-50 keys and their ordered normalized sample columns | Avoid collecting nonselected assay rows; never fit or adjust on top 50 |
| QC preview/table | Existing calculated QC product population | Requested preview columns/rows only | Avoid source S4/assay hydration |

Limpa model input is explicitly not a selective-table consumer. `dpcDE` requires
the complete aligned quantified EList—E, genes and every required `other` support
member such as standard.error and n.observations, with exact feature/sample order.
An abundance-only projection must return `saved_entity_unsupported`. Limpa
science/performance remains held for DIA-017.

## Correctness and safety contract

- Validate owned paths, Parquet/schema fingerprints, codec version, key column,
  column order and bounded serialized metadata before querying.
- Accept only validated existing columns. Quote identifiers through DBI and bind
  all paths, row keys, offsets and limits as parameters.
- Preserve requested key order and duplicates through an ordinal selector table;
  report missing keys rather than silently dropping them.
- Restore exact column order, factors/ordered levels, row names, bounded safe
  attributes, missing values, NaN and infinities after collection.
- Enforce declared row, column and collected-byte ceilings. Close result sets,
  unregister request resources and disconnect DuckDB before returning.
- Corrupt/stale/foreign tables cannot be trusted hits. RDS then slice is measured
  as full hydration and never described as projection.

## Predeclared comparison

For one-off and repeated page, volcano and heatmap-shaped reads, record initial
table construction, Parquet/RDS persistence and fingerprinting, request elapsed/
CPU, file bytes, selected rows/columns, collected object bytes, full-RDS decoded
bytes, cache counters and open resources. Use actual DIA-shaped types and order,
including factor/attribute/nonfinite edge cases.

Correctness and resource closure are mandatory. The adapter is eligible only for
consumers where it demonstrably collects fewer bytes than full RDS and avoids
source hydration. Report elapsed regressions directly; do not call reduced bytes
a speedup. If no declared consumer shows a useful net tradeoff after construction
and realistic reuse, record the treatment as not adopted and leave RDS/resident
paths unchanged. GUI activation is never authorized by DIA-013 alone.

## Implementation and correctness result

The implemented codec writes Snappy-compressed Parquet with a bounded serialized
R restoration contract embedded in schema metadata. Zstd was rejected during the
first focused run because this installed Arrow build does not include that codec;
Snappy is installed and became the single fixed compression path. DuckDB queries
bind the owned Parquet path, selector ordinals/values, limit and offset. Requested
identifiers are validated and DBI-quoted. No caller supplies SQL.

Focused `workflow-saved-projection|workflow-saved-resources` passes 291
assertions with zero failures, warnings or skips. It proves selection before
collection; safe ordered duplicate/missing selectors; exact columns/order,
factors, attributes, NA/NaN/infinities; corruption/scope/limit rejection; and
zero open resources after delivery. RDS remains explicitly full hydration.
Abundance-only Limpa projection is rejected, while equal E with changed
uncertainty or support order produces distinct identities.

All existing checkpoint/product/recovery paths pass 533 assertions after the
descriptor/recovery capability extension. NAMESPACE and public APIs remain
unchanged. DBI and duckdb are newly declared package dependencies; Arrow was
already imported.

## Frozen access comparison

Fresh installed build `/tmp/multischolar-dia013.MMf4vF/library` ran five page,
volcano and heatmap-shaped accesses against a deterministic 8,020-row, 57-column
DIA reduced-control shape: 24 normalized and 24 raw sample columns plus complete
statistics/identifiers. Receipt:

`/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/consumers.943656aeee8a8/projection/projection-consumers.json`

Receipt SHA-256 is
`3c18451f7cc644071dd2ebbfa440a09a094f3ff44bec48fb17e47d9e765b84a6`;
combined consumer receipt SHA-256 is
`97d6ae34ac44b1b3d846a999f492cd8193e44ff8756f0c12317f074e1c9637e3`.

| Observation | Whole RDS | Parquet/DuckDB | Outcome |
| --- | ---: | ---: | --- |
| Persistence + fingerprint | 0.055 s | 0.144 s | Parquet 2.62× slower |
| File bytes | 3,081,472 | 4,610,392 | Parquet 1.50× larger |
| Median requested access | 0.014 s | 0.055 s | Parquet 3.93× slower |
| Five repeated page/volcano/heatmap decoded/collected bytes | 85,310,160 | 16,025,400 | Parquet collects 81.2% fewer bytes |
| Open resources after delivery | — | 0 | Passed |
| Result agreement | Reference | Exact | Passed |

The first Parquet page paid a 0.140-second cold connection/query cost; later pages
were approximately 0.048–0.055 seconds. Full-RDS reads were approximately
0.012–0.023 seconds. Volcano projection was 0.051–0.061 seconds versus
0.012–0.015; two-stage heatmap access was 0.099–0.114 versus 0.012–0.016.
The driver included table construction, persistence, hashing, misses and every
requested read. The preceding request `consumers.940667b179e57` completed its
queries but failed receipt assembly on an unqualified tooling-only
`object.size`; it is retained and was not treated as evidence. The bounded fix
qualified `utils::object.size`; no codec or scientific work changed.

## Decision

The capability is valid and materially bounds collected data, but its declared
default performance benefit fails. Outcome: **not adopted** for GUI/default DIA
pages, volcanoes, heatmaps, checkpoints or small products. Existing resident and
RDS paths remain authoritative. The codec stays internal and unactivated as
evidence for a possible future explicitly memory-constrained consumer; reduced
bytes must never be reported as a speedup.

No result is inherited by Limpa, other omics, history integration or G2/G3/G4.
DIA-014 proceeds with complete compatibility export hydration using accepted RDS.
Rollback removes the table codec, DBI/duckdb declarations and projection scenario;
existing Parquet import support and all saved RDS records remain untouched.

## DIA-013 revision 1.3: bounded action-scope correction, 2026-09-16

The measured defect was repeated full validation, schema restoration and DuckDB
startup inside one page/volcano/heatmap action. The internal codec now has an
action-scoped owner: it fully validates the accepted descriptor and file once,
opens one connection, performs bounded reads, checks project/file identity
between reads, and closes before delivery. The default one-call interface still
fully validates; no cross-action trust cache or GUI route was added. A fixed
single-query heatmap read matched R full-population ordering in tests including
ties, NA, NaN and infinities. The Parquet format and legacy RDS are unchanged.

The focused projection/resource suite passed; the affected saved-consumer, DA,
export and recovery suite also passed. A fresh package install loaded and
registered the existing public exports and S4 classes. One installed test remains
skipped because its separate DIA-008 baseline-library bindings were not supplied.

An initial diagnostic retained a table alias through `constructed$value`, so its
memory observation was invalid. The driver removed that owner before the final
compound comparisons. Both final trials used the same installed package build,
five serial repetitions and exact results, with zero open resources after each
action. They ran on `/tmp` (tmpfs), not qualified physical-disk service storage;
therefore they are diagnostic elapsed/CPU evidence, **not** combined-peak or G2
acceptance evidence.

| DIA-shaped table | One RDS decode + three views | One scoped Parquet action | Elapsed ratio | CPU ratio |
| --- | ---: | ---: | ---: | ---: |
| 8,020 rows × 24 samples | 0.021 s | 0.101 s | 4.81× slower | 4.95× slower |
| 10,000 rows × 48 samples | 0.035 s | 0.141 s | 4.03× slower | 4.00× slower |

The 8,020-row [receipt](projection-scope-small-20260916.json)
has SHA-256 `38deccaf72bc2e2f41b1a3dd2b28ba55d59edd428b3155d73434ebfdea3918a4`;
the 10,000-row [receipt](projection-scope-representative-20260916.json)
has SHA-256 `33023d5a0229593bfd425c35431bd9184d98656c7ab49da3b9e0a414e5a15d1a`.
Both receipts are preserved with this decision record; their original run
directories were temporary and do not qualify as physical-disk evidence.
The final installed R database SHA-256 is
`570ca456b280cdeb201ef5ebdf22dc8f80092e2c0c68e33c7f73340e420f3759`.

Decision: **not adopted**. Neither tested consumer shape meets the user's strict
no-elapsed/no-CPU-regression prerequisite, so no positive memory gate or seven-pair
promotion run is justified. The lack of qualified physical-disk combined-peak
evidence remains explicit; no claim is made about larger scaling or GUI benefit.
The resident/RDS GUI policy and DIA-015/G2 hold remain unchanged. A future
positive attempt would require a newly scoped, measured mechanism rather than
relabeling these diagnostic runs as acceptance.
