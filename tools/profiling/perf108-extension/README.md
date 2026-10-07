# PERF-108 local exact-sum extension bundle

Local review bundle only. No publication, signing or production qualification.
Run `bash build.sh /absolute/owned/output` inside the existing guarded runner.
The offline build uses CMake and the unmodified upstream metadata writer; it
never downloads a dependency. `Makefile` retains the community C-API build entry.

Template: duckdb/extension-template-c at
`7f71365c5ce61b2b346717af07c9d448cfc9d3c3`, the revision declared by the official
community `capi_quack` descriptor inspected on 2026-10-05. Build tools:
`5ae521a4c9709bb3cb35dc6d79e044f67bc58bf0`. Repository source URLs:
https://github.com/duckdb/extension-template-c
https://github.com/duckdb/extension-ci-tools
https://github.com/duckdb/community-extensions/tree/main/extensions/capi_quack

`src/dia_rsum.c` is byte-identical to MultiScholaR's accepted kernel; it inherits
MultiScholaR's LGPL-3.0-or-later package declaration. LICENSE, COPYING and
COPYING.LESSER are included. DuckDB headers retain their upstream copyright and
MIT notice in src/duckdb-LICENSE; vendored tools retain their own LICENSE.
No contributor copyright or relicensing authority is invented here.

Only Linux x86_64 with 64-bit long-double mantissa is admitted. The compile-time
kernel checks remain. API v1.2.0 is a C ABI target, not the R package version;
qualification remains pinned to the installed DuckDB 1.5.5 artifacts. Both sums
must receive values in immutable source-ordinal order. NULL/empty lists return
NULL; finite nonnegative doubles are accepted; invalid values fail. No R callbacks.

The output contains the standard placeholder signature and is UNSIGNED.
A successful local build cannot resolve the production loading hold. A future
community descriptor must bind an approved published repository commit and
limit platforms appropriately; neither repository nor commit exists yet.
Do not submit a descriptor with invented coordinates. Signed artifact availability
for the pinned engine must be established before Stage 2 begins.

`test/sql/dia_rsum.test` provides portable extension fixtures. The local stage
executes the corresponding kernel checks through the pinned R DuckDB driver;
the upstream SQL test runner is not installed or claimed as executed here.
