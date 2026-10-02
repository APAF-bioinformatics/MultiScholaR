# DIA-025 full-GUI capacity trial

## Current decision — select P14; candidate tuning closed, 2026-09-30

The user chose to end marginal optimization and select from the completed
evidence, retaining eight threads and the requirement for actual fastduckdb use.
Select **P14**, frozen application build XsaPvT with fastduckdb 0.1.0.9001:
three sequential complete-sample batches, exact DuckDB/duckplyr aggregation and
the direct fastduckdb R-column writer. No additional intermediate assay files.
This is the selected integration candidate, not a claim of production readiness.
Defaults and production dependency pins remain unchanged.

The 300-sample / 24-million-precursor primary result lxQvhv is **94.866 seconds,
6.705 GiB peak**, producing 12 million peptide rows. The wide-annotation result
0bVIHx is **39.659 seconds, 4.057 GiB peak**, producing 960,000 peptide rows from
3.84 million precursors / 48 samples. Both match the accepted outputs exactly
(1,200 / 96 bounded fingerprints), all summaries and 12 required figures, with
zero swap. These remain one full pair per input, not repeated qualification.

| Candidate or family | Completed evidence | Closeout decision |
|---|---|---|
| A: frozen eight-worker R | Full operation stopped by memory guard | Preserve failure; no finite speedup denominator |
| B: bulk native R | Primary 577.847 s / 14.59 GiB | Excess memory and elapsed time |
| C: direct exact SQL | Correct screen; no full confirmation | No basis to displace the qualified relational path |
| D / P0 / P1 / P2 / P5 / P6: unbatched relational or streaming/staging changes | Full primary peaks approximately 11.13–12.15 GiB | Higher memory than complete-sample batching; historical timing gap remains unresolved |
| P3: three batches, Arrow writer | Matched P14 control 106.450 s / 6.771 GiB; both inputs exact | Superseded by direct fastduckdb writer |
| P4: six batches, Arrow writer | Primary 126.708 s / 4.928 GiB | Lowest observed full-size peak; retain as memory alternative, not the selected fastduckdb implementation |
| P7: bounded evidence lookup | Primary 104.390 s / 6.680 GiB | Small change, no consistent contrast benefit |
| P8/P9/P11/P12: added layout or result stages | Exact screens; no full confirmation of these individual profiles | No compelling screen benefit; do not add staging |
| P10: native cohort staging | Primary 109.458 s / 8.381 GiB; contrast also worse than its control | Reject added intermediate I/O |
| P13: Arrow C-stream to native writer | Exact screen, slightly slower than P14; no full confirmation | Prefer direct R-column path |
| **P14: direct fastduckdb writer** | **Exact full outputs on both retained inputs** | **Selected** |
| P15/P16: metadata cache / batch reuse | Exact preflight and twelve completed screen measurements | Marginal memory changes; not promoted or full-size qualified |
| P17/P18 | Attribution did not justify native ownership/conversion changes | Not implemented |
| P19 combination | Planned only | Cancelled when the user ended marginal tuning |

P4 is the actual lowest-memory measured primary candidate. P14 is selected
because it meets the explicit fastduckdb requirement, preserves a substantial
memory reduction relative to bulk R, and has full evidence on both input shapes.
P4 and P14 came from separate rounds: their reported time/peak differences are
observations, not a fresh paired causal comparison. No six-batch/fastduckdb hybrid
has been tested or silently substituted for either implementation.

The resumed P15/P16 preflight **0ItJ9q** passes 15 native-science fixture cases,
publication failure/recovery and the bounded non-rollup consumer. The earlier
Co6YJP mismatch was in the verifier's projected-versus-complete table boundary:
the corrected fixture independently checks source labels/ordinals and uses the
existing completeWorkflowTableBatch codec for exact container restoration.
Unoptimized fast_r is included as a fixture control. No reader or scientific
comparison was weakened.

| Screen median | Primary seconds / GiB | Wide annotation seconds / GiB |
|---|---:|---:|
| P14 | 17.0125 / 2.4392 | 14.1265 / 1.8521 |
| P15 | 17.2960 / 2.4691 | 10.9360 / 1.9323 |
| P16 | 16.8445 / 2.4257 | 11.1915 / 1.8513 |

Primary screen: 960,000 precursors / 12 complete samples. Wide screen: 480,000
precursors / six complete samples. P16's primary memory reduction is about 0.55%;
its wide reduction is about 0.04%. P15 increases memory. The apparent wide time
advantage is not attributed to these tiny writer changes: P14 itself varies from
17.488 to 10.765 seconds across an interrupted queue. No full-size gain is claimed.

All twelve operations retain exact fingerprints, summaries, output completeness
and zero swap. **CHnC1c** completed its operation and own output readback, but its
oracle verifier rejected a JSON-simplified descriptor. **1wXXDD** independently
verifies that unchanged output against the canonical native RDS descriptor;
the original failed receipt remains. Its verified selected-native fingerprints
are reused for later contrast checks, avoiding full oracle construction in each
measured service. **a6M1rH** completed during a turn interruption and was recovered
from its successful service receipt. Neither computation was repeated. The
finished queue and preserved failure/interruption snapshots reside in
browser-smoke.U0ZRfy/p14-writer-screen*.json.

**Closed scope:** no further P15–P19 tuning, combination run or new eight-member
full confirmation series. The original P14 confirmation remains the selection
evidence. No experiment service is running. Original failures and inputs remain
retained; P15/P16 remain explicit experimental profiles only.

**Next discussion:** package the selected incremental fastduckdb API reproducibly
with supported platforms, expanded-payload limits, ownership/failure contracts
and capability checks; resolve exact-sum extension loading/private DuckDB API
dependencies; then qualify the installed rollup GUI action, bounded downstream
consumers, Revert, saved-session restoration and publication recovery. Automatic
large-input selection still needs its own qualification. Shared-writer staging,
GUI integration and default activation have not been completed by this experiment.

## Retained attribution and verification history

**Prior P14 optimization stop, 2026-09-29 — superseded by the closeout above.**
The user approved the four profile-gated candidates and shared-writer staging,
with eight full confirmation members maximum and one correction pass. No timed
candidate or full confirmation member has run; no staged package or default
activation is claimed.

Diagnostic-only native build 1N8X9R retains P14's converter/encoder and adds scalar
timers/counters. Its source manifest is fast-writer-opt-source.json in the existing
proof root; the delta from the preceding native source is retained in
tools/profiling/fastduckdb-writer-diagnostic.patch. Build peak 4,822,364,160 B,
zero swap. Preflight u8VfvQ proves both completed writes and schema-failure cleanup.
Fbgcgi failed before science because the new harness omitted selection of the
protein_peptide_filtered history entry. The one correction pass restored the
established selection and design-derived sample IDs. Replacement tnlBkE and the
original contrast member cACKRn then completed, with untouched source heads,
required figures, zero swap and no failed writer appends.

| Attributed writer cost | Primary screen: 480,000 output rows | Wide screen: 120,000 output rows |
|---|---:|---:|
| Metadata serialization/base64 | 0.003283 s / 49 calls | 0.000854 s / 13 calls |
| R storage slicing | 0.018618 s / 77,852,840 accounted bytes | 0.005579 s / 145,464,584 accounted bytes |
| Native conversion | 0.029304 s | 0.047154 s |
| Concatenation worker elapsed, summed | 0.007592 s / 46,569,375 buffer bytes | 0 s / 0 calls |
| Native flush wall time, including encoding/concat | 0.044106 s | 0.179770 s |
| Peak pending Arrow charge | 16,846,800 B | 21,723,900 B |
| Whole-service peak | 2,648,576,000 B | 1,911,304,192 B |

These are profiled observations, not timing results or additive memory accounts.
R object.size charges include shared objects; concatenation worker times overlap
and are included in flush time. Profvis samples attribute 0.08/0.24 seconds to
appendWorkflowTableWriter; profiles are retained in each member. Other operation
costs dominate. P17/P18 were not implemented: the observed native conversion/concat
cost does not justify changing native buffering/worker ownership in this pass.

P15 (metadata once) and P16 (whole-batch column reuse) are implemented as explicit
private profiles in application build U0ZRfy, using unchanged native fastduckdb.
Their scientific query, three sample batches and eight-thread configuration are
unchanged. Co6YJP fails the first new non-rollup consumer fixture: projected
readback has character automatic row names, while the original data frame has
integer automatic row names. No value differences were reported. It then hits
the existing testthat muffle_expectation error while reporting that mismatch.
The later scientific/lifecycle cases did not run. The unchanged projected reader
uses character source row labels for automatic plain-data-frame row names, so
this may be a fixture/read-boundary mismatch rather than a writer regression;
no repair or additional comparison was attempted after the exhausted correction
pass. Co6YJP is preserved (1,019,355,136 B peak, zero swap).

At that stop P15/P16 were unqualified and the comparison harness unexecuted.
The authorized continuation, completed screen and final scope decision are
recorded above. Earlier P14 evidence below remains valid; all previous failures
and the original input artifacts are retained.

**Direct writer replacement series — 2026-09-29.** Following the user's
clarification, the new candidates replace P3's existing writer directly and add
no cohort/result files. P13 hands the same bounded Arrow batches to fastduckdb's
native parallel writer; P14 hands bounded R columns to its converters/writer.
The pinned Rust core already supported incremental row groups; the isolated
0.1.0.9001 build exposes that lifecycle to R, including caller codec metadata,
schema checks and abort/finalizer cleanup. Project dependency pins are unchanged.
Reproducible source delta: tools/profiling/fastduckdb-direct-writer.patch against
d78511504fbdb89fbf34bc26b22b3ff1524b0555. FKrMqv contains the frozen
native source/build receipt; XsaPvT is the installed application/extension library.

Preflight hYbi3e passes both backends' GC-with-pending-data, named-column/factor/
attribute restoration, incompatible-schema rejection, abort/finalizer cleanup,
ten exact native rollup fixtures and both publication-failure/recovery paths.
A test-local path variable misplaced the JSON receipt; the intact complete
receipt was copied to relational-preflight.json with original bytes retained and
provenance in preflight-receipt-location.json. No science or timing was changed.

All six forward/reverse screens match accepted fingerprints, summaries and
required outputs, zero swap. P3 median 17.6515 s; P13 17.7505; P14 17.333.
P14 was selected without claiming a meaningful screen improvement. Both full
pairs now pass exact ordered fingerprints (1200 primary / 96 contrast), all
scientific/filtering summaries, the 12 required figures and zero swap:

| Input | P3 Arrow writer | P14 direct fastduckdb R-column writer | Members |
|---|---|---|---|
| 24m primary | 106.450 s / 6.771 GiB | 94.866 s / 6.705 GiB | W0DRDp / lxQvhv |
| 3.84m wide-annotation contrast | 42.860 s / 3.928 GiB | 39.659 s / 4.057 GiB | GA84Dk / 0bVIHx |

P14 lowers complete-operation elapsed time by 10.9% and 7.5%, respectively.
Primary peak is slightly lower; contrast peak increases by 0.129 GiB (3.3%).
No extra cohort/result files are written. Final Parquet sizes are
476,007,933 -> 425,579,513 B (primary) and 719,611,815 -> 704,388,586 B (contrast).
These are logical scientific matches, not claims of identical physical encoding.
The primary aggregate/order/summary/write phase falls from 94.626 to 83.003 s;
its internal costs remain inseparable in these unprofiled runs.

**Recommendation: advance P14 to integrated application qualification.** This
is an actual writer replacement using the existing fastduckdb encoder, not an
added staging layer. DuckDB/duckplyr still perform the exact scientific aggregation.
The observations are one full pair per input, not formal repeated qualification.
The isolated API extension, its supported platforms and GUI/default routing still
need the agreed integration qualification; no automatic/default activation follows.
Application R database SHA256:
`32b55cf8ff138ed7c62ad8afff0bdd889aaa74fec5606d83c872cc196af56a6b`.
Source patch SHA256:
`ffc1cc0a1c970f12f3f5120da645e6760e7f4551c321a7870a63cfb0476a2bf6`.
Native fastduckdb shared-library SHA256:
`3a2758bbdd5c7b94b2608bf5ea6b5d3450f14d38d28d56447232b55b971bc171`.
Eight encoder threads, Snappy, existing codec metadata and the single durable
output/publication boundary are retained. Queue receipts are in
browser-smoke.XsaPvT/p3-direct-writers-*.json.

Affected regression result OxPWzQ: all 406 projection/codec and 293 resource
assertions pass. The broader collection is not clean: two launcher expectations
omit the existing fixed_reclaim binding, and the installed-compatibility check
has one failure and a follow-on error because it still requires zero native DLLs.
Both predecessor 7wnUsR and candidate XsaPvT load one native library; the reclaim
launch option predates this series and is documented in the retained history.
These unrelated expectations were left unchanged. Preserve the full 1496-check
receipt with three assertion failures and one error; do not relabel it as a pass.
The preceding zaDtIE runner dispatch error occurred before tests; its NULL-selector
guard was corrected for OxPWzQ. This and the receipt filename repair are verification
bookkeeping changes only; no successful timed member was repeated.
The staged hybrids below remain a separate unfavorable comparison, with all
evidence preserved. Defaults remain unchanged; direct writer experiments are complete.

**Completed staged P3 / fastduckdb mechanism experiment — 2026-09-29.** The user explicitly
requested real P3 variants using fastduckdb. Reader qualification is now complete:
w0yrMj restores the recorded class/attributes/order in the harness and verifies
every actual prepared request plus failure cleanup. The new oracle references
the unchanged QnPzsh files and preserves Rh2Ssg. Paired vJ2lKS/raxrDm/iW63oM/FUDFid
members pass exactly at eight threads: combined median DuckDB 2.953 seconds,
fastduckdb 2.202; index 2.2055 versus 1.090; wide 0.7075 versus 1.0735. Thus reader
choice is workload-specific, and these are not rollup speedup measurements.

Build 7wnUsR freezes P8–P12 and unchanged P3 as a matched control. Installed
preflight 9I6dT8 passes 25 native equality fixtures, actual fastduckdb-call checks,
staging cancellation cleanup and five publication-failure/recovery lifecycles.
R database SHA256:
`1adfb1530a08348f150f07a7fd3d4e1546e643d245baab7b3f638a20e93b9dc6`.
The exact-sum binary remains unchanged. Native cohort row-group size was probed
as 262144 and matched in the DuckDB staging control. All new profiles retain
three complete-sample batches, eight threads per sequential engine phase,
12 GiB DuckDB/no spill and the common 16 GiB process guard.

| Profile | Changed boundary | Screen median, s | Median combined peak |
|---|---|---:|---:|
| P3 | Existing control | 19.6815 | 2.462 GiB |
| P8 | Projected input aligned by sample rank/source ordinal | 19.5595 | 2.521 GiB |
| P9 | P8 plus DuckDB cohort Parquet staging | 19.5675 | 2.567 GiB |
| P10 | P9 with Rust-native fastduckdb cohort filtering/writing | 19.262 | 2.617 GiB |
| P11 | P3 plus final batch Parquet and bounded DuckDB reads | 19.3395 | 2.410 GiB |
| P12 | P11 with bounded fastduckdb result reads | 19.8235 | 2.493 GiB |

All twelve forward/reverse screen members match the accepted fingerprints,
scientific/filtering summaries and 12 required figures, with zero swap and no
remaining private stage files. P10 makes three native fastduckdb streaming calls;
P12 makes 48 bounded fastduckdb reads. P10 is the selected fast finalist, without
claiming a meaningful screen speedup: the observed difference is small and P3
varies between repetitions. The declared full pairs are complete, all exact,
zero swap, complete figures and clean private staging:

| Input | P3 time / combined peak | P10 time / combined peak | Members |
|---|---|---|---|
| 24m primary | 106.326 s / 6.830 GiB | 109.458 s / 8.381 GiB | UiONBQ / 2IsEZN |
| 3.84m wide-annotation contrast | 48.812 s / 3.885 GiB | 49.650 s / 4.097 GiB | RxLXHi / okC3qK |

P10's primary operation writes 1,403,939,618 additional intermediate bytes,
with 924,721,631 peak live intermediate bytes. Layout preparation takes 2.659 s
and its three native cohort stages total 1.963 s, all included in operation time.
At the largest sampled current-memory points, P3 anon/file counters were
4.748/1.104 GiB and P10 5.093/1.957 GiB; these are attribution observations,
not components summed to manufacture cgroup peak. Whole-service CPU was
304.804/326.541 s and physical writes 994,582,528/3,291,299,840 B, respectively;
service accounting includes startup and independent readback outside operation time.

**Decision: retain P3; do not promote this added-staging design.** P10 is slower
and uses more memory on both full inputs. P12 also loses to its matched result-read
control in the screen. These are staged hybrids: P10 replaces P9's DuckDB cohort
filter/writer, and P12 replaces P11's DuckDB staged-file reader. Neither is a
direct replacement of the unchanged P3 path. The user clarified that the next
series must target such direct component replacements. The pinned fastduckdb R
API lacks incremental writing with caller-supplied codec metadata and a query/scan
interface, so inspect the Rust source before designing the smallest necessary
extension. Do not add staging merely to include fastduckdb or present these
results as a test of a direct replacement. All layout, staging, summary/writer,
publication and cleanup costs were counted. Queue evidence is in
browser-smoke.7wnUsR/p3-fast-mechanisms-*.json. Default activation remains off.

**P3 adaptation continuation, 2026-09-29: resumed on explicit user instruction.**
Publication-context repair TqcXDE passes the pending/failure/retry/actual Apply
lifecycle (1,195,642,880 B peak, zero swap); djt8ku republishes the three affected
selected populations while preserving their assay references and all old inputs.
Build bvT5nP freezes P6 (remove the redundant all-sample filter) and P7 (P3 plus
one bounded evidence lookup). Installed preflight zOmRA3 passes ten exact fixture
comparisons, both publication-failure/recovery paths and resource cleanup.
Installed R database SHA256 is
`7c8d5182c37ef40833172d618ccda91f7ff84cd6950c3e2aeacd34f87dfc57e4`;
the exact-sum extension remains
`48d54592464571b7e159de1a7dafeb3af0e4f39d70d7d73f7d5acb984f3f542b`.
The installed relational source matches the working-tree file exactly.
Eight screen and four full members completed with exact output hashes, summaries,
12 figures and zero swap. No default change.

| Profile | Fresh screen median, s | Fresh full operation, s | Full combined peak | Full member |
|---|---:|---:|---:|---|
| P0 packaged control | 16.741 | 104.006 | 11.390 GiB | 3ftOV2 |
| P6 omit redundant all-sample selector | 16.8645 | 104.452 | 12.011 GiB | ecnA9D |
| P3 three complete-sample batches | 17.797 | 105.112 | 6.995 GiB | t9tE1t |
| P7 P3 plus bounded evidence lookup | 17.741 | 104.390 | 6.680 GiB | BueO9n |

P3's substantial memory reduction repeats; its earlier time advantage does not.
P7 provides a modest 0.315 GiB lower observed peak, without a meaningful throughput
gain in this single full pair. P6 offers no measured benefit. The matched P0 was
added before any full execution because the fresh screen control changed markedly;
four full members ran, without retries. The old 72.290-second result remains
unexplained; the redundant selector does not explain that gap in this comparison.
Machine-readable queue records are in browser-smoke.bvT5nP/p3-adaptations-*.json.

Untimed Nu4Fn9 confirms all 200 Parquet row groups span run_0001–run_0300 and
overlap every P3 batch; each sample has exactly 80,000 rows. The aggregation plans
show projection/filter pushdown, but sample statistics cannot prune row groups.
No within-input size skew justifies a weighted batching policy on this workload.
Future layout changes must count preparation costs; no input rewrite was made.

I/O preflight dqjmpt failed at hashing arithmetic-created versus restored R NA
payloads before reaching fastduckdb. Bounded diagnosis found identical values,
types, attributes and missingness. The verifier-only correction canonicalizes
non-NaN NA payloads and implicit block row names on both sides, preserving all
finite bits and NA/NaN distinction. Bw0qqS passes both eight-thread readers and
four injected-failure cleanups; QnPzsh prepares the fixed index/wide inputs.
Rollup science, tolerances and its existing 1200-block verifier are unchanged.
Crossover is complete with exact native/P0 agreement in all twelve fresh members:

| Precursor rows | Native B median, s | Single-batch P0 median, s |
|---|---:|---:|
| 9,600 | 7.658 | 7.9435 |
| 96,000 | 9.074 | 8.606 |
| 960,000 | 21.6345 | 13.876 |

Native retains the small-input advantage; a clear relational benefit appears at
960,000 rows in this schema. This is not a universal row-count routing threshold.
The earlier valid 4608-row pair is retained, not combined across timing periods
to interpolate a crossover.

The contrasting input (3.84 million precursors, four per peptide, 2048-byte
accession annotations) also passes exact agreement across all 96 output blocks,
scientific/filtering summaries and 12 figures, with zero swap:

| Arm | Complete operation, s | Combined peak | Member |
|---|---:|---:|---|
| Native B, single process | 102.765 | 6.248 GiB | EYBmKx |
| P3, eight threads | 41.756 | 3.797 GiB | OlkPVn |
| P7, eight threads | 42.293 | 3.822 GiB | aFOgPr |

P3 lowers elapsed time by 59.4% and peak memory by 39.2% on this contrasting
input. P7's small primary-input memory advantage does not generalize here.
Recommend retaining P3 as the production-qualification candidate; neither P6 nor
P7 demonstrates a consistent improvement sufficient to replace it. Weighted
batch admission remains uncalibrated, and changing ordered aggregation state
would require separate scientific qualification. Do not add those mechanisms
without the relevant evidence. Reuse existing ownership/codec/writer/publication
helpers and operation-specific scientific contracts; no generic executor was
extracted without a measured second transformation consumer.

**Historical I/O qualification gate, resolved above.** First timed reader member Rh2Ssg stopped
at its DuckDB index verifier, before fastduckdb ran. Bounded diagnosis confirms
identical values but a different top-level attribute order introduced by the
new harness's manual column assembly. Applying the existing recorded attribute
order reproduces the prepared reference SHA256 exactly for the inspected page.
The small preflight missed this because its oracle used ordinary data-frame
subsetting rather than the actual index buffer. The failed receipt and original
prepared inputs remain intact; no I/O performance result is accepted. Additional
verifier-only repair authority has been requested under the repository's
correction-budget rule. The completed rollup comparisons are unaffected.
R/shell syntax and diff checks pass; all experiment services have stopped.

**Historical preflight failures, preserved.**
The user authorized the staged adaptation/experiment plan recorded in DIA-025.
The attempted fixture-only repair rebuilds selected-population progress through
finishDIAReferencePeptideStep and uses the existing publication/failure/retry
machinery, preserving Ny7QzY and vwKtL9. It has not passed verification.

Preflight yJ4Kdn failed an expected pending-output assertion before any publication
or scientific execution (combined peak 992,391,168 B, zero swap). The single
correction enabled reference_states on the copied adapter. Recheck lXK7vs still
failed at that assertion (992,272,384 B, zero swap): the copy also lacks
experiment_paths$time_dir, so resumeDIAQvalueOutput returns NULL. This is a new
preflight setup defect, not a P3 scientific/performance failure. The surrounding
test error handler additionally reports the testthat muffle_expectation restart;
both original logs/receipts are retained. Neither run reached its intended
completion/failure publication lifecycle, so neither qualifies instrumentation.

That correction allowance was spent and work stopped until the user's continuation.
Initialize the fixture's ordinary publication context before checking readiness
in the next authorized repair; preserve the pending-output assertion and all
production gates. P6/P7, scan attribution, remaining crossover/contrast/fastduckdb,
memory-weighted policy and shared-executor decisions were not executed then. Existing
P0–P5 scientific/performance evidence below is unchanged; defaults remain off.

**Prior production qualification gate — 2026-09-29, superseded by the continuation above.**
The package implementation and P0–P5 primary-dataset comparisons are complete.
Automatic GUI/default routing remains off pending the agreed results review and
remaining qualification. Crossover preparation has a publication-contract defect;
fastduckdb comparisons and the contrasting-input timings have not run. This is
not completion of the approved production rollout plan.

| Profile | Screen median, s | Full operation, s | Full combined peak | Member |
|---|---:|---:|---:|---|
| P0 packaged control | 20.363 | 120.383 | 11.477 GiB | CLXTzy |
| P1 incidental insertion order off | 20.8145 | 119.334 | 11.396 GiB | 9K0pJV |
| P2 public DBI/nanoarrow streaming | 20.2145 | 114.747 | 11.893 GiB | q4paZu |
| P3 three complete-sample batches | 21.6255 | 108.490 | 6.822 GiB | 60v50g |
| P4 six complete-sample batches | 22.920 | 126.708 | 4.928 GiB | td5UV6 |
| P5 aggregate-Parquet intermediate | 20.3085 | 121.673 | 12.147 GiB | QPGvHL |

All completed screen/full results have exact matching fingerprints, scientific and
filtering summaries and required figures; all swap peaks are zero. P3 uses its one
permitted affected large replacement: original TdS3ji was scientifically correct
but its reported 1710.573 s operation timing is invalid. System journal confirms
suspend at 2026-09-28 20:39:58 UTC and resume at 21:06:30 UTC; memory samples have
a 1591.969 s gap, while systemd's monotonic service interval is 179.216 s.
The original receipt and `timing-validity.json` remain retained. No sleep duration
was subtracted to manufacture a corrected timing. The existing verified sleep
inhibitor now protects every bulk-rollup reader member; oOBKgM passed its lifecycle
check before 60v50g. All jobs and inhibitors have ended.

Interpretation remains observational, with no preset performance trade-off target.
Complete-sample batching materially lowers peak RAM on the primary input. True
streaming and staged Parquet did not reduce peak there; P5 also wrote a
245,041,391-byte intermediate. The packaged control's 120.383 s is slower than the
prior 1wDSnG result of 72.290 s. This regression is unresolved, not silently
attributed to packaging or dismissed as noise. The current control adds sample
selection/batch plumbing; the diagnostic plan shows an additional all-sample IN
filter/mark join. Its contribution versus other runtime differences has not been
isolated. No existing large execution allowance remains for repeating that audit.

Package work: relational functions moved from test helpers into
`R/utils_workflow_dia_relational_rollup.R`; the existing saved-rollup entry point
accepts an explicit qualification profile while its default remains unchanged.
Configure builds a separate exact-sum extension from pinned headers/source and
writes compiler/flags/source/binary metadata. The fixed-path loader validates
Linux x86-64, long-double precision, DuckDB 1.5.5, duckplyr 1.2.1 and binary SHA.
Nanoarrow 0.9.0 is isolated and declared for the true-streaming candidate. Public
S4 interfaces, checkpoint formats, writer and publication/recovery machinery are
preserved. No runtime dependency installer or arbitrary extension path is added.
The public-streaming candidate uses DBI's streaming result type; the ordinary
Arrow-batch control materializes the engine result before draining.

Build HZ0QIR and installed P0–P4 preflight so7Gxt passed. Diagnostic preflight
si563u failed because JSON profiling was enabled after its .json output path;
the wrapper also lacked setup-error cleanup. Correction pass 1 fixed both.
Sh9YNU passed completion/failure instrumentation. Full diagnostic 74wKXS matched
all 1200 accepted fingerprints and summaries. Main-query RSS was about 1.06 GB
before execution, 8.43 GB after execution, and 1.88 GB after draining. Its plan
contains aggregation, two joins and sorting; this justified P5. Engine profiling
reported buffer counters exceeding measured RSS; these are not summed or treated
as cgroup measurements. Build DNMnfd froze P0–P5; P5-only preflight UtzYHe passed.
Installed R database SHA256:
`46a8926d632cc117b789e1f50fd954975d6c2c1fdbfcae8496c6d15a5b034335`.
The installed extension remains
`48d54592464571b7e159de1a7dafeb3af0e4f39d70d7d73f7d5acb984f3f542b`.

Screen members, forward then reverse: QGWkjJ, ao9UCg, VQ4MK1, j3kDRZ, AHxNJu,
xTgx05, 9dlW6C, iy64u5, JpxJNW, wz7juo, FkgfdC, y8EQbX. All 12 passed.
Each used 960,000 input rows, 480,000 outputs and 12 figures. Full members used
24 million input rows, 12 million outputs and 12 figures. Guards remain eight
threads, 12 GiB engine, 16 GiB hard / 15.5 GiB early / zero swap / 16 GiB startup /
2 GiB host; five-minute screens and 30-minute full members. Profiling was disabled
in timed members. All earlier experimental failures remain preserved below.

**Qualification inputs and limiting failure.** Preparation preflight 3FOYCI
passed creation and injected failure cleanup. Ny7QzY prepared immutable 4608,
9600, 96000 and 960000-row cases, plus the declared 3.84-million-row contrasting
input with four precursors per peptide and 2048-byte accession annotations.
The completed 4608-row crossover has exact native/relational agreement:
ZNOhpH B1 8.044 s; Zkqtva P0-1 8.564 s; YHsY6t P0-2 10.924 s;
kt7RNz B2 9.923 s. Medians B 8.9835 s and P0 9.744 s; this single tiny point
is insufficient to choose a crossover.

vwKtL9 (B, 9600 rows, repetition 1) then failed **before scientific execution**
with `QC outputs are pending; retry Apply to publish them`. The preparation
copied the predecessor's `qc_progress` product into a new selected-population
checkpoint, which has a different reference and no corresponding publication
acknowledgement. The product is also bound to the old population. The production
gate correctly refuses to proceed. Do not clear the flag, acknowledge the old
product, or weaken the gate. The next bounded repair must construct correct
selected-population progress and publish it through the existing machinery,
retaining the old preparation and failed member. Verify a prepared checkpoint
through actual Apply/publication before resuming affected crossover members.
The original 4608-row cases and passed primary-data measurements can be reused.

Both correction passes are consumed (diagnostic setup; sleep-inhibitor routing),
with two affected small checks and the single affected large rerun used. The queue
stopped at vwKtL9 and no further code correction or measurement was made. Renewed
correction authority is needed to repair fixture publication and finish the
remaining crossover/contrast work. Existing-API fastduckdb I/O draft work is not
executed or qualified; no fastduckdb speed/memory benefit is claimed. Automatic
routing thresholds, ordinary-GUI activation, selected-profile confirmations and
the fresh integrated GUI journey remain pending. R/shell syntax and diff checks
pass; no other omics/scientific transform has been migrated.

**Current — DuckDB and duckplyr comparison complete (2026-09-28).** Both
relational adapters now preserve the native numerical contract using the isolated
`dia_r_sum` extension. Recommend the **duckplyr relational candidate (D)** for
subsequent integration planning. The selected full-size member **1wDSnG**
completed in **72.290 s**, versus native bulk B's **577.847 s**: **87.5% lower
complete-operation time (7.99x)** in this exploratory confirmation. Peak combined
memory fell from 14.59 to **11.13 GiB**. Exact full scientific/output verification
passes. No GUI/default activation or formal repeated qualification is implied.

| Arm | 12-sample screen median, complete operation | 24-million-row confirmation | Combined peak at full size |
|---|---:|---|---:|
| A — frozen R, eight workers | 89.8345 s | Early memory stop; no finite baseline time | 15.51 GiB observed at stop |
| B — current bulk R, one process | 28.277 s | 577.847 s, verified | 14.59 GiB |
| C — DuckDB + exact-sum extension | 16.3665 s | Not selected for full-size execution | — |
| D — duckplyr relational API + exact-sum extension | 16.6185 s | **72.290 s, verified** | **11.13 GiB** |

The fixed new C–D–D–C screen passed: AHuDlS C1 16.167 s / 2,859,511,808 B;
NlH6P3 D1 16.885 s / 2,677,637,120 B; ocYGlW D2 16.352 s / 2,676,219,904 B;
QY8v1L C2 16.566 s / 2,860,290,048 B. Zero swap throughout. All 48 ordered
fingerprints, scientific/filtering summaries and the 12 declared sample selectors
match accepted native B exactly. Source checkpoint, installed package and dependency
identities match the native comparison; C/D additionally freeze the extension and
new adapter sources. Experimental harness source identities differ from the earlier
native runs because C/D adapters were added; native scientific/output machinery
is unchanged. Reused native evidence is not relabeled as newly executed.

C and D times differ by less than 10%; D wins the frozen lower-peak rule. Its
screen median is 41.2% below B, exceeding the 10% finalist requirement. Only one
new full-size member ran, D, with no retry. Full direct DuckDB performance remains
unmeasured; do not extrapolate its near-equal screen time into a full-size result.

1wDSnG processed all **24,000,000 precursors -> 12,000,000 peptide rows** and
published the durable checkpoint and all **12 required PNG/PDF figures**. Its
1200 bounded readback fingerprints and complete scientific/filtering summaries
match 7X1Wpw, whose independent comparison against BpSyxB's accepted complete-sample
native calculations already passed. This reuses the valid native reference instead
of repeating its computation. Verification receipt:
`browser-smoke.1wDSnG/bulk-rollup-independent-verification.json`, which hashes
both compared receipts and embeds the prior native verification binding.
Input/head identities remain unchanged. Payload is 476,007,933 B; figures total
1,814,168 B. The Parquet bytes need not hash identically: logical values, order,
types and attributes match exactly after only implicit block row-name normalization.

Full operation: validation **1.729 s**; aggregation, ordering, summary accumulation
and streamed writing **62.613 s**; enclosing Apply/checkpoint **68.030 s**;
publication **4.160 s**. Apply overlaps the preceding phases and is not additive.
Whole service including fresh-process preparation and bounded readback was
**103.185222 s**. Whole-cgroup CPU **253.405921 s**, physical read/write counters
**7,127,179,264 / 990,142,464 B**. Peak **11,945,205,760 B**, zero swap.
Sampled maximum anonymous memory **11,896,307,712 B**, file cache **1,142,308,864 B**;
these maxima are noncoincident and cannot be added. DuckDB's internal aggregation,
sort and output costs remain inseparable; no per-engine-thread CPU attribution is
claimed. A separate cumulative staged-byte counter was not retained. Final staged
payload sizes and physical I/O are reported instead.

All timed relational members used the existing runner, eight DuckDB threads,
**12 GiB engine memory**, spill disabled, **16 GiB hard / 15.5 GiB early / zero swap /
16 GiB startup / 2 GiB continuous host stop**. Screens had five-minute ceilings;
full D had 30 minutes. No spill error or spill file was observed; final pending
staging is empty and scopes/readers/connections close. No service or worker remains
live. Limits and numerical tolerances were not changed to obtain these results.

Implementation: `helper-workflow-relational-rollup.R` holds experimental C/D
validation, query construction and streaming; shared native/output helpers moved
unchanged in substance into `helper-workflow-bulk-rollup.R` so all tests load them
before execution. C uses one SQL aggregation; D uses duckplyr's exported relational
expression API with pinned DuckDB bridges to extract/register the lazy relation.
It is **not stock duckplyr sum()**: both use `dia_r_sum(list(quantity ORDER BY
source_ordinal))`. `read_sql_duckdb(..., con = scope$connection, prudence = "stingy")`
uses the same scoped engine; unsupported construction fails and R fallback/telemetry
remain disabled. The original ordinary/compensated-sum failures remain valid.

SQL input/group/evidence checks follow the existing validation rules; factors,
attributes and canonical order are restored through the existing codec. Distinct
annotations are streamed and collapsed with `collapseRollupProteinIds()`. Quantities
are restricted to the declared double schema; integer/other types fail explicitly.
Only lazy views and the bounded annotation lookup are created. `dbSendQuery(arrow =
TRUE)` plus `duckdb_fetch_record_batch(..., chunk_size = 10000)` avoids DBI's full-R
resultset behavior. Each batch goes directly to the existing writer and filtering
accumulator; no full-assay temporary SQL table or full-result R collection is used.
The scoped action owns extension connection teardown. Production package APIs,
checkpoints, dependency files and defaults are unchanged.

**Numerical and lifecycle proof:** source-ordered long-double accumulation matches
this frozen R build's 64-bit mantissa, including both earlier mismatches. Direct
DuckDB and duckplyr both pass **507 groups**, including 500 seeded randomized
positive groups, missing/zero values, subnormals, overflow and an order-sensitive
extended-precision case. `tools/profiling/dia_rsum.c` has compile-time precision
checks and never calls R from DuckDB threads. The numerical approach follows
[R's accumulation behavior](https://github.com/wch/r-source/blob/trunk/src/main/summary.c)
and the [DuckDB C extension interface](https://github.com/duckdb/extension-template-c).

Full adapter preflight **dd7pcH passed** exact C/D ordinary, empty, missing,
factor/unused-level, annotation and numerical-boundary fixtures; invalid quantity,
logged/type, duplicate/evidence rejection; cancellation cleanup; both publication
failure/recovery paths and 12 figures each. Peak 1,395,404,800 B, zero swap.
Preserved failures: the first numerical binary advertised unsupported stable C API
v1.5.5; corrected to stable v1.2.0. The next duckplyr expression duplicated its
grouping output; renewed authority removed that duplicate and both numeric gates
passed. Adapter preflight **3wsKXu** failed only because the tiny test readback used
an unsupported Arrow selective-reader call; the bounded correction uses the existing
DuckDB selective reader for that tiny oracle. Scientific C/D streaming is unchanged.

Evidence root remains `MultiScholaR-dia-proof-20260911-VhC1n5`. The
`relational-rollup-exact-sum/` directory retains source, pinned v1.5.5 headers,
metadata script, binary, failed check, corrected check and the initial failure receipt.
That initial receipt is historical; dd7pcH and the completed timed members supersede
its held status. C source SHA256
`7f76be1f3f80ef2fd70d00ba0d36d50b5b19c53572ecf46e446a0008384abddd`;
extension SHA256
`48d54592464571b7e159de1a7dafeb3af0e4f39d70d7d73f7d5acb984f3f542b`.
The binary is built with `cc -std=c11 -O2 -fPIC -shared -fno-fast-math -Wall -Wextra
-Werror`, all three `DUCKDB_EXTENSION_API_VERSION_{MAJOR,MINOR,PATCH}` defines set
to 1/2/0, and the retained pinned headers. The retained upstream footer script
uses platform `linux_amd64`, API `v1.2.0`, extension version `dia025-v1`, ABI
`C_STRUCT`. Only these experimental connections allow this local unsigned extension.
R/shell parse checks and actual guarded routing pass. All C/D source and extension
hashes match between screen and full D. This bounded comparison ends here;
production integration, broader omics and repeated qualification are subsequent work.

**Prior native-only comparison — completed before relational repair (2026-09-28).** Recommend
**bulk current native R (B)** as the experimentally supported implementation for
this precursor-to-peptide boundary. It preserves the full accepted scientific
result and durable outputs within the fixed guard, and reduces small-screen
complete-operation time by 68.5%. Production integration/default activation is
not performed or authorized by this result. Original-speed parity and a matched
full-size speedup remain unproved because the frozen eight-worker baseline
cannot complete within the guard. Preserve that limiting failure; never use its
cutoff time as a speedup denominator.

| Arm | Small-screen median complete operation | Full 24-million-row outcome | Full combined peak | Decision |
|---|---:|---|---:|---|
| A — frozen x3Dkg2, eight R workers | 89.8345 s | naNmMJ: early memory stop during kernel; no complete output | 16,654,983,168 B observed at stop (15.51 GiB) | Retain limiting baseline failure |
| B — current native R, one process | 28.277 s | 7X1Wpw: complete operation 577.847 s; independent exact verification passed | 15,665,553,408 B (14.59 GiB) | Recommended experimental candidate |
| C — DuckDB 1.5.5, ordered compensated sum | Disqualified before timing | Numerical mismatch | Not measured | Do not adopt |
| D — duckplyr 1.2.1 translated sum | Disqualified before timing | Numerical mismatch | Not measured | Do not adopt |

Scope is precursor-to-peptide only. No IQ/MaxLFQ, imputation, GUI migration,
public API, checkpoint-format or production-default changes. Native workers use
one internal numerical thread. C/D use the existing scoped connection, eight
DuckDB threads, 12 GiB engine memory and the existing no-spill settings.

The fixed A–B–B–A screen used samples 1–4, 149–152 and 297–300 in declared design
order: 960,000 precursor rows -> 480,000 peptides. All four members had identical
source/build/dependency bindings, selectors, all 48 ordered output fingerprints,
scientific summaries and filtering summaries, and each published 12 figures.
A1 **2rhuRz**: 83.957 s, peak 12,171,870,208 B; B1 **Sh06Vb**: 28.130 s,
2,970,763,264 B; B2 **mHCMDZ**: 28.424 s, 2,968,113,152 B; A2 **6dM49o**:
95.712 s, 11,205,124,096 B. All zero swap. Whole-service intervals were
99.181855 / 43.312715 / 43.740761 / 111.196256 s respectively. These include
fresh-process preparation and bounded readback outside the operation timer.
B passed the frozen >=10% improvement rule and was the sole finalist.

Full A then B consumed exactly the two full-size attempts, with no retry:
30 minutes each, 16 GiB hard / 15.5 GiB early / zero swap / 16 GiB startup
headroom / 2 GiB continuous host stop. A's read/assembly completed in 441.482 s;
the controller then stopped it during the unmodified frozen kernel. B completed
in 577.847 s (9m 37.847s); whole service including preparation/readback was
624.206674 s. B phases: read/assembly **402.816 s**, native kernel **128.533 s**,
summary plus staged writing **36.189 s**, required figure publication **6.514 s**,
teardown **0.001 s**. The enclosing Apply phase is 571.124 s and overlaps the
read/kernel/write phases; do not add it to them. Input read/assembly accounts
for 69.7% of B's complete operation. No further optimization is included here.

**Independent full verification passed.** B produced 12,000,000 peptide rows.
All 1,200 bounded 10,000-row block fingerprints match BpSyxB's retained native
complete-sample rollup calculations exactly, including ordering, types,
attributes, values and missingness; only implicit row-name representation is
normalized at block boundaries. The complete scientific summary is identical:
842,978 partial-missing outputs and 17,956 all-missing outputs for each quantity,
frozen identification evidence and accession provenance preserved. This reuses
the accepted complete-sample native reference because full A failed; it does not
substitute that earlier run's timing for the failed baseline. Source input and
head identities match. The result table is 476,007,933 B; all 12 required PNG/PDF
figures total 1,814,168 B. Checkpoint, publication, bounded readback and teardown
completed. Candidate and stopped-baseline pending directories are empty; no
services/workers remain live. Original source head and validated artifacts are
unchanged. No full-result R collection was used for independent verification.

Resource attribution: B whole-cgroup CPU is 615.704051 s, disk read/write
counters 7,177,830,400 / 1,012,654,080 B. Its sampled maximum anonymous memory
is 14,590,287,872 B and file cache 1,243,447,296 B. A's stop observation records
606.386482 s cgroup CPU and read/write counters 4,178,780,160 / 9,109,504 B;
sampled maximum anonymous memory is 16,532,131,840 B and file cache 968,445,952 B.
Sample maxima are noncoincident and must not be added to each other or to the
combined peak. Kernel-internal phases and per-worker CPU are **not separately
attributed**: process PIDs/residency and total cgroup CPU are retained, but no
per-worker CPU total was captured. Cgroup I/O counts physical traffic, not all
logical cached reads. Final payload/figure sizes are known; a separate cumulative
staged-byte counter was not retained. Native scientific execution has no spill
engine; SQL preflight explicitly prohibits spill. These measurement limitations
remain explicit rather than inferred as zero or reconstructed as exact costs.

Evidence root: `/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5`.
Each member retains `bulk-rollup.json`, logs, controller bindings and available
terminal/memory/process receipts. For stopped A, use `guard-stop.txt`: its
phase receipt still says running after external termination and is not a live
job or completion certificate. Independent verification is retained in
`browser-smoke.7X1Wpw/bulk-rollup-independent-verification.json` and
`verify-dia-bulk-rollup.R` (SHA256
`32962b9a6dface9a098b481766597b9d3a7e099d73867b00ecb58853ba95230d`).
Reference result: `ae00dcaa952f5838892afef219fb47407959e0e50ffb7bce0fe570426483c76c`;
candidate: `5bb6e4fb12d3ecb4b831d107a0a5e9c494ab0d3238ac384eb092f0b644a6069f`.

Relational rejection is exact, without tolerance changes or R repairs. For
`c(1e16, 1, 1e16, 2)`, native R gives `20000000000000004`, while ordered
`fsum(quantity ORDER BY ordinal)` gives `20000000000000000`. Ordered compensated
sum passes `c(1e16, 1, 1)`; duckplyr's translated `sum()` gives
`10000000000000000` instead of native `10000000000000002` on that case.
JSON rounds these numbers at this magnitude; exact `identical()` eligibility
booleans and retained fixture inputs establish rejection. Duckplyr enters via
`duckplyr::read_sql_duckdb(..., con = scope$connection, prudence = "stingy")`;
its registered methods use namespace-qualified dplyr generics. Automatic
materialization/fallback and fallback collection/uploads are disabled. Only the
two-row diagnostic result is collected. Since these scientific gates fail,
complete C/D streaming/publication adapters and timed members were not pursued.
This does not reject every conceivable relational aggregation algorithm.

Implementation is confined to opt-in helpers/tests in
`tests/testthat/test-workflow-memory-mechanism-member.R` and exact guard/selector
routing in the existing shell runner, R controller and DIA driver. The current
buffer, unmodified native kernels, staged writer, accumulator, checkpoint and
publisher are reused. A loads only frozen rollup/private helpers from x3Dkg2,
with current shared I/O/publication and dependency binaries; historical S4
methods are not registered. Kernel object ownership is unchanged. R parse,
shell syntax and helper-definition-order checks pass. No project dependency
files were changed for the isolated duckplyr installation.

**Preserved preflight/harness failures and accepted corrections:**
- **NN2fyx**: fixture directory missing; 6.969099 s main interval,
  802,988,032 B peak, zero swap. Created the directory in the allowed correction.
- **BQpBPJ**: full sequential callback has no `dia_source_order` attribute;
  15.611342 s, 2,591,449,088 B peak, zero swap. User renewed authority; derive
  full-stream ordinals from offset while preserving selected-row ordinals.
- **6cbEnD**: session-dependent publication wrapper had no active manager;
  40.672468 s, 8,959,352,832 B peak, zero swap. Renewed correction calls the same
  underlying publisher with the explicit owned manager and unchanged arguments.
- **dcsHuz**: JSON cannot serialize DIAFilteringSummary; 47.370681 s,
  9,069,445,120 B peak, zero swap. Its last JSON says running because error-flush
  serialization also failed; terminal/test logs establish failure. Renewed
  correction unclasses only the receipt copy, leaving the scientific product intact.
- **AswwYh passed**: complete A/B lifecycle, exact results, both sets of 12
  figures, publication failure/recovery, cancellation lease cleanup, injected
  worker failure and stopped workers. Native fixtures include empty, zeros,
  partial/all missing, invalid quantities/logged data, duplicate identities,
  conflicting evidence, factor levels, annotations and both numerical boundaries.
  Peak 9,058,512,896 B, zero swap; this valid proof was reused.
- **a2liEU**: A1 test ran before copyDIAReplayCheckpoint was defined;
  16.665525 s, 796,442,624 B peak, zero swap, no scientific work. User renewed
  authority; moved the test below existing helpers and reran only A1 as 2rhuRz.
  Four successful screen members plus this failed launch remain retained.

All source/package hashes match across successful screen and both full members.
Duckplyr 1.2.1 is isolated in `bulk-rollup-library`, using existing DuckDB 1.5.5.
Duckplyr source archive SHA256:
`7cae0967602aea89554c5acc52c729f1bade09eaf9a021735feef62247070400`;
frozen x3Dkg2 rollup source:
`681bb93b536c5ccb66bdb1bb95ad5ce556e27eb90aab4defac51a3651fb415e5`;
installed IsMJrU R database:
`2ef17e48f5368212ec633959a7681db8c6368621a8f63b7058ce0a37a7db9aff`.
Member receipts record dplyr/Arrow/multidplyr versions and binary hashes too.
The bounded experiment ends here. GUI integration, promotion, repeat
qualification, extra concurrency points and further reader optimization require
a subsequent bounded plan; DIA-025's broader GUI objective is not marked complete.

**Current — scientific workers execute correctly; original-speed parity unproven.**
Four-worker preflight **lhElgo passed**, including four distinct scientific
worker PIDs, exact serial-native result checks for all four jobs and injected
worker-error cleanup. Full input is the retained 12-million-row peptide table,
all 300 design observations, fixed 9,600-row tasks comprising complete proteins.
The coordinator holds one input cohort and streams ordered results; each task's
input was at most 522,832 B and scientific output at most 945,240 B.

| Member | R scientific workers | Input read | Worker startup + imputation + ordered writes + teardown | Whole service | Combined peak | Swap |
|---|---:|---:|---:|---:|---:|---:|
| yInxfd | 4 | 89.109 s | 271.585 s | 375.552229 s | 7,166,885,888 B (6.67 GiB) | 0 |
| tlOfJc | 8 | 58.761 s | 152.297 s | 219.844046 s | 11,420,143,616 B (10.64 GiB) | 0 |

Both passed: 1,250 jobs, 12 million data rows, 6 million support rows, selected
independent serial-native comparisons, clean worker/resource teardown, unchanged
source head. Every eight-worker scientific job fingerprint (data, support and
summary) matches four workers. Eight workers reduced the measured worker phase
by 43.9% in this exploratory pair; read-time variation is reported separately.
This is actual parallel R scientific execution, unlike earlier reader-thread tests.
It is not original implementation speed parity, complete GUI integration or
production qualification. No full journey was restarted.

User challenges the overhead relative to previously near-instant native steps.
That comparison remains missing. Repeated execution of the R method on 1,250
small tasks may repeat setup/grouping/validation, but its contribution versus
intrinsic calculation, six-million-row support generation and I/O is not yet
measured. The next necessary audit is a matched single native R call on the same
input/settings, with those costs separated; use the frozen pre-refactor method
as an additional authority where its science/output contract matches. Do not
claim parallel scaling repairs the end-to-end regression. Duckplyr was assessed
from official documentation but not installed or benchmarked; no backend or
dependency change has been made.

**Current — full peptide table fits one cohort; speed benefit not established.**
User requested the 768 MiB peptide cohort experiment. Guard-routing klZlzg passed.
Both members use IsMJrU, two Arrow CPU threads, one I/O thread, the same frozen
BpSyxB replicate checkpoint, 16 GiB hard / 15.5 GiB early stop / zero swap,
16 GiB startup headroom, 2 GiB continuous host stop and a 900 s ceiling.
Only the diagnostic reader copy changes the three cohort-budget checks.

| Member | Budget | Cohorts / total rows | Read + assembly | Whole service | Combined peak | Swap |
|---|---:|---|---:|---:|---:|---:|
| kBOWuI control | 256 MiB | 3 / 12,000,000 | 67.286 s | 89.913453 s | 2,840,535,040 B (2.65 GiB) | 0 |
| fqJrBY candidate | 768 MiB | 1 / 12,000,000 | 75.164 s | 94.708569 s | 4,264,083,456 B (3.97 GiB) | 0 |

Both passed. The candidate's assembled table occupies 627,221,104 B by R's
object-size accounting (not its total residency). All 1,250 ordered 9,600-row
block fingerprints match, including values, nonfinite encodings, column types,
attributes and source ordinals. Only implicit row-name representation is
normalized at comparison block boundaries; explicit row identifiers are retained.
Original source head is unchanged and reader resources close cleanly.

The candidate was 11.7% slower for read/assembly in this single pair and used
more memory. This demonstrates feasibility for the seven-column peptide
projection, not imputation/rollup kernel capacity, multiple workers or an
end-to-end speedup. Do not promote the larger budget as a speed improvement.
Production/default cohort limits remain unchanged. The deferred eight-thread
point, R-worker scaling and whole scientific-stage validation remain unexecuted.

**Current — cohort copy correction verified before concurrency experiments.**
Initial ownership-transfer build oVpwmz failed focused check **yd8ME0**:
values/attributes remained correct, but later batches still copied the traced
destinations. A 20-element reproduction showed temporary local bindings retained
an additional reference after return. Explicitly clearing those temporary
owners fixed the issue without per-batch GC or native mutation code.
Corrected **IsMJrU** passed core/frozen-science check **snINlQ** (49.950124 s,
peak 1,487,441,920 B, zero swap).

Same-real-cohort candidate **U8O7Od passed** (62.184934 s complete service,
peak 2,889,039,872 B, zero swap). Its exact result fingerprint matches ki2V0N.
Assay callback assembly fell 167.015 -> 7.548 s; traced quantity copies fell
449 -> 1 and position copies 448 -> 0. The two traced vectors' cumulative copy
lower bound fell 24,085,862,400 -> 35,788,800 B. Complete measured cohort work
fell 188.239 -> 44.314 s. Index/decode time varied (assay decode 14.451 ->
25.667 s), so retain phase costs and do not present this one diagnostic pair as
formal qualification or whole-workflow speedup. No limits or science changed.
Affected integration is pending; no concurrency result exists yet.

**Current — measured cohort copying defect before parallel scaling.**
User accepted the isolated scaling proposal. Serial attribution **ki2V0N passed**
on one real 4,473,600-row / 3,728-complete-protein cohort from BpSyxB, using
FIDXq2 (198.488 s service, peak 2,915,852,288 B, zero swap). Ordinal-index
stream: 4.405 s, including 0.086 s callback assembly. Assay stream: 181.466 s,
including 167.015 s callback assembly and 14.451 s decoding/restoration.
Tracing one double destination and the integer positions vector recorded
449 and 448 copies, respectively: at least 24,085,862,400 B of copying for
those two vectors alone. This is cumulative allocation/copy work, not peak
resident memory. The exact cohort fingerprint is
f4e9951c1ab59eb9e541aca2c33c3302ea2c2aea7730de9a51c6e1b74868d92a.
Original head was unchanged; no full diagnostic object was retained.

The bounded correction transfers private vector ownership out of the container
before subassignment, shared by existing ordered/grouped assembly consumers.
No batch/cohort/cache limit, scientific kernel, dependency or backend changed.
Build oVpwmz installed; focused parity and then the same real-cohort comparison
are required before starting any concurrency ladder. No parallel result exists.

**Current — BpSyxB explicitly stopped by the user (2026-09-28).**
The unit is inactive. Last observed progress: imputation 6,832/10,000 protein
groups; accepted state replicate_filtered. Completed phase times including
commit/output work: qvalue 384.108 s, evidence 277.758 s, precursor rollup
1536.139 s, peptide intensity 842.410 s, sample filtering 90.314 s, peptide
replicates 949.655 s. Last observed combined memory peak 4,529,827,840 B
(4.22 GiB), zero observed swap. This is not a final peak certificate or a
completed capacity journey. Controller exit 76/missing_memory_observation
occurred after the explicit service stop removed the cgroup; no memory guard
breach or timeout was observed. Preserve stale running receipts as historical
artifacts; this explicit stop record is the current decision.

User now requests a proposal for isolated parallel scaling on the retained large
dataset, followed by reintegration into the guarded journey. No parallel
implementation or new benchmark has been launched. The completed serial phases
are valid evidence; IQ/protein QC/normalization/DA/Summary capacity is unproven.

**Live 300 replay: BpSyxB (2026-09-28).** Guard routing **uxtNYV passed**
(29.728003 s, peak 1,421,897,728 B, zero swap). BpSyxB is running the bound
FIDXq2/cuurUw raw-checkpoint route, with
24,000,000 input rows. Actual systemd limits were checked: 16 GiB hard,
zero swap, 90 minutes. The existing controller retains 15.5 GiB early stop,
16 GiB launch headroom and 2 GiB continuous host-memory stop. First active phase
is qvalue. This is not yet a successful replay or fresh GUI journey. No source,
test or runner edits are permitted during the live member; inspect/resume it
before starting another job. The first approved large replay allowance is used.

**Current — small lifecycle and retained-output gates passed (2026-09-28).**
FIDXq2 native-session suite **i91nq8 passed** (447.458976 s, peak 3,046,178,816 B,
zero swap), including ten normalize/revert/style/export cycles, scientific
comparisons, resource/GUI ownership and native/legacy compatibility. The replay
now supplies an explicit IQ QC PDF path and performs Summary's third session-state
export. Prior replay IQ wrote its default PDF into the test working directory;
that unbound file was not used as verification evidence.

Output-only recovery/verification **ZqcqbU passed** (153.026965 s, peak
3,593,404,416 B, zero swap). It restored Vyr2MC's native archive into independent
storage, regenerated figures from saved summaries/results and IQ QC plots from
the saved preprocessing input, published Summary and exported session-state
metadata. It did not rerun MaxLFQ or create another scientific checkpoint.
Recovery forbids DA refit/full-S4 reconstruction; the subsequent independent
verifier explicitly reconstructs one required historical record at a time and
uses the existing DA kernel as an oracle. IQ/QC inventory, additional protein/
RUV/Pearson progress grids, Summary copies, all history states, final quantities
and DA TSV/workbook comparisons passed. Original visible source files were
fingerprinted before/after and unchanged. Vyr2MC remains a failed member; this
is recovered output/science evidence, not a successful timed GUI journey.

Next bounded member: one cuurUw raw-checkpoint-to-Summary 300 replay on FIDXq2,
90-minute ceiling; exact installed RDB/native-library and source-head bindings
are in the existing diagnostic specification. Preserve 16 GiB hard/15.5 GiB
early stop/zero swap/16 GiB startup/2 GiB continuous host limits. Profiling stays
off. Validate the exact guard route before launch. The large correction budget
remains one additional affected member after a bounded correction; no large
member has run in this continuation yet.

**Current — focused publication/teardown repair verified (2026-09-28).**
Build Eimw8a isolates Glimma dependency generation in an owned temporary
directory, copies only self-contained HTML to publication staging, and releases
temporary assets on return/error. Focused retained-output check **uTLNAN passed**
(9.664465 s, peak 1,062,776,832 B, zero swap): nonempty HTML only, no external
asset-directory references, no DA refit or full-S4 reconstruction, unchanged
scientific state/source head, successful closed-session finalization and an
injected finalization failure recorded correctly. Vyr2MC remains a failed member;
this output check does not retroactively qualify its incomplete output inventory.

First affected integration **w2vOXE failed** (99.274078 s, peak 1,593,044,992 B,
zero swap), exposing a real Glimma compatibility defect: annotation preparation
hardcoded Protein.Ids despite args_row_id declaring Protein.Group. The bounded
correction uses `.data[[args_row_id]]`. Final build **FIDXq2** passed affected
integration **qoZqBu** (105.480063 s, peak 1,608,376,320 B, zero swap), including
reference DA publication, native-equivalent tables/projected views, resource
ownership and affected normalization/export consumers. The added test checks
one nonempty interactive HTML output per contrast. No full scientific journey
was repeated for this repair. All earlier failures remain retained.

This bounded repair is complete. Remaining gates: independent complete output
verification, final native lifecycle checks, bound 90-minute 300 raw-checkpoint
replay and, conditionally, the 120-minute fresh GUI run. New build bindings must
be recorded before large execution. No 300 member has launched in this repair;
no full-workflow performance or capacity claim follows. The requested final
playbook remains due after the journey outcome.

**Earlier focused output repair — correction authority was required.**
User authorized finalization/output repairs and checks against retained Vyr2MC
results without repeating its science. Replay captures ownership before session
end and finalizes from the session-end receipt; finalization errors now write a
failed receipt. The output audit found production gaps: normalization/Pearson
omitted publication_graphs_dir when calling filtering plot output, and reference
DA omitted interactive-volcano files expected by Summary. Path arguments are
fixed; output request identities are versioned so older journals cannot suppress
the added files. DA stages linked-expression HTML via the existing publication
journal and requires its output files before acknowledgement.

Build cL1lmY / focused member **wgnOMa failed** (9.172754 s, peak 1,003,298,816 B,
zero swap): portable Glimma dependencies include two empty DataTables CSS files,
which violate the publication journal's nonempty-file contract. Bounded
correction switched to self-contained HTML in build **32yLMg**. Member
**331wIE failed** (10.312758 s, peak 1,368,834,048 B, zero swap): the valid
1,354,370-byte HTML was produced, but saveWidget also retained the asset directory
containing those empty CSS files. Both checks passed no-DA-refit, no-full-S4,
unchanged scientific state/source head, closed-store finalization and injected
finalization-failure assertions; neither is an overall output pass. Original
Vyr2MC files were not modified. The focused initial attempt and one correction
are consumed. Next concrete repair: use an owned temporary asset directory
outside publication staging, publish only self-contained HTML, then rerun the
focused check and affected integration. Do not weaken journal integrity or
launch large execution. Full output/science verification remains outstanding.

**Latest 2026-09-28 — authorized additional preflight failed after Summary.**
User authorized one Pearson callback correction and one additional small replay.
The test driver now passes the callback rather than its evaluated list. Member
**Vyr2MC** on unchanged DqtvOc completed all scientific phases, Pearson,
native-session export, DA views and Summary (237.262695 s service elapsed,
peak 2,895,101,952 B, zero swap). It then failed because the replay's final
`diaStateOwnership()` call ran after `testServer()` had ended the session and
closed the store. The existing session-end receipt records zero full-S4 reads,
zero working records and zero open resources. Its replay receipt remains stale
`running`; the failed service result is authoritative, not a completion claim.
Publication also warned that Interactive_Volcano_Plots was absent: replay renders
the Glimma widget but does not invoke its file export. Two filtering-plot warnings
report incomplete publication path arguments. Therefore Summary's boolean flags
are insufficient output evidence. Preserve all outputs and verify the actual
inventory before accepting this replay. The additional attempt is consumed;
next scope decision is a bounded replay-finalization/output-wiring correction
with focused checks and reuse of retained outputs where valid, not another
automatic whole-science replay. No 300-sample run has started.

**Current 2026-09-28 — small replay correction allowance exhausted.**
Corrected build DqtvOc passed affected integration **xPelq6** (89.354758 s,
peak 1,634,783,232 B, zero swap), including streaming-progress regression.
Small replay retry **eB9LHF failed** (163.304514 s, peak 2,917,416,960 B,
zero swap) after all peptide/IQ/protein QC and normalization/RUV phases completed.
Pearson scientific commit was reached, but its publication failed because the
new replay caller supplies `aesthetics()` (a list) instead of the `aesthetics`
callback required by `runDIAReferenceCorrelation()`. This is a replay-driver
API mismatch, not evidence of a numerical or capacity regression. Preserve the
failed output journal and receipts; do not call this a completed replay.
The small work item's initial attempt and one correction retry are consumed.
Per AGENTS.md, stop before further edits/runs and request authority for the
one-line driver correction and one additional small replay. The 300 replay,
fresh GUI run, independent output verification, final native lifecycle check
and requested playbook remain outstanding. No large-run allowance was consumed.

**Earlier today — native summation and replay preflight evidence.**
Final native-summation build **BlBrjE** passed core/frozen parity **yqsVz8**
(66.132736 s, peak 1,599,483,904 B), retained-selector comparison **VObGAY**
(50.112136 s, peak 2,926,710,784 B), and routing/Design checks **CKLXpI**
(19.119299 s, peak 1,469,702,144 B); all zero swap. These do not prove full
journey performance. Small raw-checkpoint replay **rLWt8p failed** after qvalue
publication, during evidence filtering: progress reporting treated a streaming
schema prototype's zero rows as the final output total. The strict observer
correctly rejected completed rows greater than that total. Scientific outputs
were not compared after this stop. The failure receipt records its active phase.
The bounded correction reports an unknown total for schema-v3 streaming writes;
footer validation still establishes the actual final row count. A regression
test covers writing nonempty batches from an empty prototype. Corrected build
**DqtvOc** installed successfully (48.097780 s, peak 1,138,495,488 B, zero swap).
The affected integration and retry outcomes are recorded above; no 300-sample
replay or fresh GUI run has launched.

**Earlier resumption — native summation verification pending at that time.**
Progress-enabled integration check 3VR5Wd passed on jRSc5A (166.326669 s,
peak 1,616,400,384 B, zero swap). Final review identified that rowsum accumulation
can differ from native sum at a threshold boundary (1e16 + 1 + 1). The kernel now
retains native R sum within the shared grouping; frozen tests require identical
quantity vectors. Installation perQSz was interrupted during byte compilation:
its package directory is empty and no valid terminal build receipt exists. It is
not a successful build or scientific failure. Preserve it and install afresh.
No job was live at resumption. The subsequent tool-host reset occurred before a
replacement launched. Comparisons now resolve the same frozen index artifacts so
instrumentation-only changes cannot rebuild indexes in original evidence.
Checkpoint replay, final GUI proof and the requested playbook remain pending.

**Reader repair verified; kernel repair under verification.**
Build 0u0XNw's core check xWVGDd passed (90.645996 s, peak 1,480,019,968 B,
zero swap). GyA198's earlier failure is retained: the old comparison mock violated
the selected-row callback contract, and a new assertion expected original rather
than projected row-name storage. Correcting those tests did not change the build.

Unprofiled frozen-selector control **50i8o8 passed** (229.331774 s, peak
3,551,653,888 B, zero swap); reader-only candidate **fsiuJm passed** (75.390051 s,
peak 2,964,635,648 B, zero swap). Across nine complete-sample measurements, median
gather elapsed fell 20.197 -> 3.262 s; assembly callback elapsed fell 10.183 ->
0.057 s; decoding/restoration fell 9.691 -> 2.713 s. Kernels were unchanged in this
comparison. All sample results matched frozen y0oLNN oracles at the existing
1e-12 tolerance, with exact summary equality. These are bounded mechanism results,
not full-journey speedups. The first control nm8VvI failed when repetition two hit
an existing output filename; retain it. Unique per-repetition paths fixed that
instrumentation error. Completed control/candidate source checks run after cleanup
and redirect all diagnostic I/O logs outside the original evidence store.

The remaining contiguous-group path uses one forward cursor, verifies actual group
identities, and delivers complete groups even across read-batch boundaries. Build
6Jo9qv / core check **MYrQyr passed** (91.996030 s, peak 1,492,279,296 B, zero swap).
No native mutation routine or new dependency was required.

Kernel build wgnO2p / check **Q1pvAj failed** (80.422994 s, peak 1,460,940,800 B,
zero swap). Frozen-kernel ordinary/empty/integer/missing/factor/annotation/rounding
comparisons passed, but the reference adapter still called the removed local
annotation helper. Correct the call to the shared helper; one correction is under
verification. The optimized kernel shares group construction across numerical and
missingness statistics, vectorizes R aggregation, reuses uniqueness validation,
and orders once. No full GUI or checkpoint replay has been launched in this repair.
Corrected kernel build **Q2q4D2** passed core/frozen-kernel check **zI0oCy**
(89.998527 s, peak 1,486,012,416 B, zero swap). Combined retained-selector check
**Nb9lTF passed** (84.804558 s, peak 2,974,515,200 B, zero swap). Compared with
frozen control 50i8o8, median kernel time fell 3.105 -> 0.435 s and median gather
20.197 -> 4.885 s. All frozen sample outputs and summaries passed existing parity
checks. Gather time varied from the reader-only member; no universal timing claim
or full-workflow extrapolation follows from these measurements.

Bounded progress and Design-driver wait reduction are under verification on
**jRSc5A**. The common writer receives an optional observer from its owning store;
it does not discover DIA session state. Receipts distinguish I/O/group delivery
from scientific commits, are throttled and atomically replace an <=8 KiB latest
record. Design assignments set existing controls and await the exact production
summary and idle acknowledgement once per pair; final saved-design validation is
retained. No checkpoint replay or new full GUI journey has launched.

The requested reusable playbook is a final handoff deliverable; the broader
sequential audit will be planned afterward, not applied speculatively here.

**Shared execution repair in progress — measured attribution before optimization.**
Approved continuation is recorded in DIA-025. No large GUI job is live. Diagnostic
preflight 2u28ce failed bounded JSON serialization of the S3 summary; ck7avq passed
its correction, including completion and injected-failure receipt checks.

Frozen SW6Naw retained-input diagnostic **y0oLNN completed**: 148.322575 s main
interval, peak 3,312,828,416 B, zero swap. First/middle/last complete sample gathers
(80,000 rows each) took 36.393 / 35.966 / 37.933 s; kernels took 5.993 / 5.314 /
5.008 s. Middle-sample assay streaming split into 18.354 s callback assembly and
17.302 s decoding/restoration. Its summary took 0.069 s and write 0.140 s.
The profiled kernel spent substantial time in per-group dplyr data-pronoun
resolution; no full-workflow speedup is inferred. Frozen sample output oracles and
bounded phase receipts are retained in browser-smoke.y0oLNN.

Attribution limitation: its source-size/mtime assertion preceded deferred scope
cleanup. Cleanup appended operational cache-release diagnostics to the retained
store's checkpoint-io.log; assay payload bytes were unchanged. Future diagnostic
cleanup logging is redirected to its request directory and assertions run after
scope release. Do not characterize y0oLNN as preserving every source log byte.
A bounded tracemem reproduction confirms closure-based nested column assignment
replaces the target vector on successive incoming batches. Sparse Arrow fragments
therefore amplify copying. First repair: coalesce sparse Arrow batches before R
restoration, use compact selected-row maps, and retain one cursor for identity-
ordered reads, within existing limits. No new native mutation routine or dependency
is added. Build 0u0XNw is under preparation; correctness and timing are pending.

**Terminal outcome — cuurUw timed out during precursor rollup.**
SW6Naw/OQADra, invocation `292b81ef22154d0788b4f9c0f5772ece`: authoritative
service timeout, main interval 7200.075596 s, peak 10,408,431,616 B, zero swap;
no engine spill files observed. Last delivered state `protein_peptide_filtered`.
The smoke receipt remains stale `running`; the service is failed/inactive.
No retry started. The complete 300-sample capacity/performance objective is unmet.

Post-run read-only audit: checkpoint observations were raw_data_s4 at 14:18:18,
qvalue_filtered at 14:32:05, and protein_peptide_filtered at 14:55:20 (local time).
The rollup sample-group index finished at 15:10:22; its staged Parquet was
223,676,830 B and last modified at 16:02:03, near the timeout. Thus there is evidence
of ongoing rollup output work, not merely a live process or unchanged GUI state.
Cgroup CPU total was 15,291.212383 seconds (all job processes); this does not
attribute CPU to reading, R kernels, summary accumulation, rendering or writing.

Confirmed remaining mechanisms: repeated cohort scanners/inverse maps in the shared
group reader; complete-sample native rollup calls with separate missingness/sum
groupings, duplicate/evidence checks, joins and repeated ordering; per-batch
fastduckdb reopening in identity-ordered map reads. No phase profile exists to rank
these costs. The new producer-fed display summaries removed their post-write
traversals but did not eliminate this underlying scientific/group execution cost.
The same group reader is used by downstream peptide intensity/sample/replicate/
imputation and uncached IQ lineage reads; ordered maps/ledgers share the ordered
reader. Protein QC, normalization/RUV/Pearson and DA mainly use different matrix/
product paths and have not been exercised at 300 samples by this member. They are
unqualified at this scale, not shown to have this same bottleneck.

**Current work — authorized end-to-end summary reuse and explicit 120-minute journey.**
The user requested a provenance plan, removal of redundant summary work across the
included DIA IQ journey, regression verification, and a complete GUI launch with a
120-minute allowance. The ordered plan and exact boundaries are recorded in
[DIA-025](../../../tickets/openai/DIA-025.md#authorized-end-to-end-summary-repair--2026-09-27).

**CueTuk superseded and stopped, not passed or timed out.** Exact-build preflight
jGsetl passed. CueTuk used 2IZpDk/OQADra with reclamation and unchanged 60-minute
guards, invocation `f5d7c71640ec47f09432f96ee738f3db`. It completed Design, evidence
counting and accepted-table/ledger writes, but not first-QC summary/figure delivery.
Last sampled interval: 2840.826076 s; sampled peak 10,431,299,584 B. The operator
stop for the newly authorized repair left systemd inactive/dead. The controller
exited 74 when its sampler encountered the removed cgroup; the smoke receipt's
`running` status is stale. Preserve all files and do not classify this as a complete
capacity result. No replacement baseline attribution diagnostic was run.

Consumer inventory and repair contract:

| Journey stage | Summary source and required change |
| --- | --- |
| Import / Design | Existing bounded import catalogue and component metadata; preserve metadata-only reuse. Serial GUI assignment waits are a separate cost. |
| First QC | Accumulate exact progress while final accepted rows are written; eliminate both post-write grouped summary traversals. |
| Peptide evidence / intensity / sample / replicate / imputation | Feed the same accumulator from final transformed/generated output batches. No-op components reuse a saved count product. |
| Precursor rollup | Accumulate from already-required complete-sample kernel outputs; retain the native scientific rollup summary separately. |
| IQ / protein QC | IQ's numerical work remains necessary. Register progress from the already-required accepted protein object once; reuse its product for output work. |
| Normalization / RUV / Pearson | Cache exact presence counts by physical assay, declared design/sample population and software. On a miss, query up to 64 sample columns together rather than one query per column. Skip/style/retry must reuse valid counts. |
| DA | Existing saved selection products already contain full-population counts; preserve their cache, FDR population and no-refit boundaries. |
| Summary / native publication | Existing bounded state metadata and output products remain sources. Archive dependency closure includes new count products and physical inputs. Legacy report compatibility stays explicit. |

Peptide accumulator workspace is capped at 256 MiB; dictionaries at 32 MiB and
final products at 64 MiB. Identity codes and raw presence vectors implement exact
distinct counts, including missing identifiers; no approximate counts or altered
scientific denominators. Producer summaries include accepted row counts, checked
against the committed table. Compatibility misses get one projected sequential
read; valid hits read no assay. Unchanged public counters remain test oracles.
Build 1etITd / core check **GJG5SE failed** (92.923429 s main interval,
peak 1,561,444,352 B, zero swap). Existing scientific assertions and the new
all-peptide-actions no-post-production-scan assertion passed. New checks exposed
an out-of-file test helper dependency and protein-summary row-name inconsistency
between native and cache-miss production. The single correction makes the fixture
self-contained and retains native named sample-count row names. Retain GJG5SE as failed.
Corrected build **SW6Naw** passed core check **qdTYiH** (91.013719 s main interval,
peak 1,508,433,920 B, zero swap). All covered native scientific outputs and exact
summary populations pass, including all peptide actions forbidding post-production
summary scans, missing/empty/factor cases and no-assay-read product hits.
Runtime RDB: `c55102a392a3384021f979af0536d734a2863350a5321472ffecad9cbfd7f7c8`;
native library: `def04822fb78abb395b099c218075b389ff6b18dd94b65dd58857d63362af15c`.
The explicit `--runtime-120-minutes` route is separately bound to this build and
only the reclamation-enabled native-first complete 300 journey. Existing 900/3600
routes remain unchanged. Routing/completion/failure check **2yGajr passed**
(31.033997 s main interval, peak 1,545,981,952 B, zero swap), including shell limit
propagation, unchanged memory stops and rejection of an unapproved extended build.
Affected GUI/DA integration **Cdexba passed** (162.018865 s main interval,
peak 1,599,651,840 B, zero swap). Native lifecycle/archive **aEfyj3 passed** (630.169947 s main interval,
peak 3,052,453,888 B, zero swap), including all ten lifecycle cycles and the covered
native/legacy compatibility, scientific, ownership and output checks. Runtime R
sources match the installed SW6Naw source copy exactly. These are correctness and
lifecycle results, not proof of no elapsed-time regression: this native suite took
longer than earlier suites, and no matched timing comparison was conducted.

Launch decision: one complete SW6Naw/OQADra 300-sample journey, native-first,
reclamation enabled, profiling disabled, `--runtime-120-minutes` explicitly selected.
The unchanged memory/swap/headroom/host guards remain. No automatic retry is
allowed. Full journey/required-output success is pending the running member. The 120-minute launch is now authorized by the passed checks. Historical
3600-second limits/results remain unchanged, and a future
7200-second result cannot be represented as a matched performance comparison.

**Current outcome — zASKAX timed out; bounded reader correction verified on small checks.**
The 300-sample candidate failed at the unchanged one-hour limit during first QC.
Authoritative service result: `timeout`, inactive, main interval 3600.104531 s,
peak 10,428,100,608 B (9.712 GiB), zero swap. Design completed and the last durable
state was `raw_data_s4`; no accepted first-QC result was published. The stale smoke
receipt is not a running status. Preserve terminal-properties.txt and
terminal-interval.json in browser-smoke.zASKAX. No 300-sample restart was launched.

The user authorized a bounded correction of this throughput defect. The group
reader materialized the entire projected assay once per cohort and reopened the
index reader every 10,000 rows, permitting repeated row-group decompression.
These are confirmed code mechanisms, not a profile assigning the entire timeout.
The replacement reads each cohort's index in one Arrow stream, maps physical to
scientific group order, and pushes assay ordinal selections into the scanner.
Contiguous ranges use balanced range predicates; scattered selections use a
bounded membership predicate. Existing complete-group/cohort/read limits remain.
Shared group and ordered-row consumers receive this correction through the codec.
No scientific kernel, global population, threshold, backend or dependency changes.

Verification uses existing core/scientific fixtures plus a 98,304-row subset of
the existing representative recipe, four interleaved cohorts, exact output and
read-volume checks, and an isolated full-materialization/selection comparison.
The comparison holds the index reader fixed; it is not a frozen-build workflow
benchmark or a measurement of physical disk bytes. Row groups shared by cohorts
can still be read more than once.

Build XtNl8N / check **2mq1xa failed**: the membership-expression branch omitted
Arrow's required skip_nulls option; two new test assertions also used base
subsetting that dropped the units attribute which the codec correctly preserved.
The bounded correction supplies skip_nulls=FALSE and repairs those assertions.
Main interval 90.484662 s, peak 1,557,573,632 B, zero swap.

Corrected build **2IZpDk**, RDB
`8a460a1d76c9c1c7be5692c419dca86e5d8297432a994b4dbf5ec2227b246212`, passed
core check **mXIHOh**: 97.236964 s main interval, peak 1,501,409,280 B, zero swap.
Native first-QC/peptide/rollup oracles passed, as did exact cohort order/coverage,
attribute restoration, empty/sparse selections, cancellation and existing limits.
The 98,304-row/four-cohort mechanism comparison measured 2.341 s full-assay
materialization versus 2.028 s pushed selection. Total index-plus-assay rows
materialized in R fell from 491,520 to 196,608 (assay alone: 393,216 to 98,304).
The candidate timing additionally includes exactness assertions. This is one small
mechanism observation with the same new index reader on both sides; it does not
quantify the old repeated-index-decompression cost, physical disk I/O, or 300-sample
runtime. Affected GUI/DA, normalization/RUV, ownership and publication integration
check **pYoITS passed**: 168.996153 s main interval, peak 1,620,713,472 B, zero swap.
The bounded implementation/check exit is satisfied. Full 300 capacity remains
unmet; no additional large attempt has been launched. Before another complete
journey, the next capacity decision should establish the corrected reader's
throughput on retained 300-sample data under the unchanged guard; small timing is
insufficient to promise completion. No production/default activation is claimed.

**Historical launch — zASKAX (2026-09-27).**
ZQ8Tbj's corrected native-session suite **5WhxhF passed**: all ten lifecycle cycles,
logical binding/physical dependency archive closure, restored S4 attributes,
missing-component rejection, owned DA handoff and scientific/output assertions.
Elapsed 391.780209 s, peak 3,218,989,056 B, zero swap. Preserve the earlier failed
uISw91 result; it was not reclassified.

Final build: RDB `1dbc963d9eedfc8eba2b168df928cf470b14e04a4e2928ccf781dba7f2814e7b`,
native library `87f875e6a44653851882fc35759e6c44b2d6df323cff47d79b276cec4b89fa21`.
BDNuFm failed a test-local uninitialized spec lookup when checking this newly bound
build; corrected that fixture, and **hwD6jQ passed** exact-build routing plus browser
completion/disconnection. The application build did not change for that fixture fix.

Launched the approved single complete 300 candidate **zASKAX**, invocation
`5d1fa00f77a943e1b64bb6ea5a787547`, on retained OQADra with reclamation enabled,
DuckDB spill disabled, profiling off and existing passive process sampling. All
full-size guards remain unchanged. This consumes the planned candidate member;
there is no automatic restart/correction member. Outputs, Summary and teardown
must complete before acceptance; independent required-output checks follow a
successful run.

The failed original-baseline diagnostic Snxb8W remains failed. Its namespacing fix
passed a module-scoped small check, but permission for a replacement diagnostic is
still pending and no replacement was run. Candidate reclamation receipts provide
new-build observations; they do not retroactively decompose mU50nZ's 6.26 GiB.
No production activation or matched performance qualification is claimed.

**Shared-fix final verification in progress — 2026-09-27.**
Managed connections now explicitly disable the DuckDB temporary directory as well
as setting its quota to zero. The initial zero-quota-only check kbrurK found an
engine temporary-file entry despite the quota; retain that failure. VHV9ng's
**nCY90x passed** the corrected forced-allocation failure/no-spill/accepted-data
checks and the affected native scientific oracles. **pXQiy6 passed** the broader
GUI, DA selection-index, resource and publication regressions.

The compiled reclamation receipt now records bounded process/allocator observations
before GC, after its existing GC, and after its existing trim, plus post-GC live R
bytes. It adds neither extra collections nor allocation profiling. In 6uOiX0's
uISw91 lifecycle, all ten cycles and the logical-binding archive/restore check
completed. Cycle 2 measured anonymous memory 1,493,458,944 B before/after GC and
1,042,194,432 B after trim, with live R 695,020,640 B. These are overlapping small-
fixture observations, not the original 300-sample baseline.

The overall **uISw91 suite failed** one owned-export integration check: low-level
store checksum validation invoked describeDIANativeFile, which that check forbids
for owned DA handoff to prevent archive rehashing. The actual path was a metadata
record checksum, not an archive read. Corrected the layer boundary by extracting
shared hashWorkflowFile and making both store validation and archive description
call it independently. Build **ZQ8Tbj** installed; its native-session verification
is running. Preserve uISw91 (378.739902 s, peak 3,159,601,152 B, zero swap) as failed.

The failed baseline diagnostic and pending explicit request for one replacement
remain unchanged. No full 300-sample candidate has run for this implementation.
No production activation or complete capacity/performance claim is made.

**Shared fix implementation progress — 2026-09-27 (not final qualification).**
The approved baseline diagnostic **Snxb8W** stopped on a namespaced Shiny request
acknowledgement timeout, before Design setup/Save. Preserve this instrumentation
failure: 222.430530 s, peak 10,685,767,680 B, zero swap. Its final process sample
attributes about 4.675 GiB anonymous to application R, 0.176 GiB to controller R
and much smaller browser allocations; there is no post-import GC/trim observation
and no decomposition of the original 6.26 GiB Design baseline. The request now
reads the root-session input; the corrected module-scoped fixture passed **xAAIBT**.
An explicit asynchronous request for one diagnostic replacement is pending; none
has been launched and the final candidate allowance remains one.

Implemented and checked so far:
- Logical `WorkflowTableComponent` records use existing RDS and reference immutable
  Parquet, with validated class/attribute/column-projection metadata. Component
  checkpoints now emit v2; native archives emit v3 and retain older readers. Native
  closure expands logical table sources without assay hydration.
- Design reuses the physical assay for unchanged rows/values. Genuine changes use
  ordered batch writers. **JSwvd5/pyydxp passed** the no-assay-scan/hash/rewrite test,
  shared readers and native peptide-state fixtures.
- Complete-group processing uses reusable narrow group/ordinal indexes and bounded
  cohorts. **swKCHC/ZrSRnr** matched science but failed an obsolete expected SQL-stream
  counter (expected 1; actual 0). The expectation now expresses the new no-global-
  cursor contract; failure retained.
- Change maps are committed without a copy, ordered with bounded integer/key
  indexes, and merged sequentially. Bounded ordinal range reads replace ordinary
  SQL OFFSET paging. **QUFb19/TAlZIF passed** all affected core checks.
- Registered file-specific prefix advice synchronizes completed staged prefixes
  every 64 MiB. Store hashing reuses the existing drop-behind hash pipeline.
  Original report/FASTA leases release at idle import delivery; import and Design
  now participate in opt-in allocator reclamation. **BAsxcr/4bNsDY passed** resource,
  byte-preservation, source-change, idle-owner and affected GUI/scientific checks.
- First QC uses the existing R confidence/exclusion/evidence rules over complete
  protein groups and sequentially writes accepted rows/ledgers. Progress counters
  use native R counters over complete indexed groups. **DzGa8q/ZJJnNB passed**.
- Persisted rollup uses the native kernel on complete samples, with independent
  global protein evidence/annotation validation. Annotation and lineage reads are
  projected/bounded. **fpGb7G/UprAYT passed**, including existing rollup oracles.

Current build bwDT2G adds `max_temp_directory_size = '0B'` on managed database
connections and cached DA selection indexes. Its zero-spill/core check is running.
The full native-session lifecycle, archive-version compatibility, remaining
consumer regressions, final no-op/index review and 300-sample candidate are still
required. No production-default activation or complete scaling claim is made.

**Approved shared-fix implementation — baseline attribution first (2026-09-27).**
The user approved the comprehensive shared production/ownership plan, reusable
narrow checkpoint indexes, and the sequence of one 300-sample import/design
attribution diagnostic followed by one complete fixed candidate. Reclamation and
new managed policies remain candidate-only. Existing scientific contracts, 16 GiB
hard / 15.5 GiB early / zero swap / 60-minute / 16 GiB launch headroom / 2 GiB host
stop remain unchanged. No historical budget or failed receipt is reset.

The existing runner now supports `--diagnostic-design-memory`: frozen PCPabv,
retained OQADra, original parser and original design assignments; stop before Save
Design. Capture startup, completed import and completed Design setup, then full GC
and allocator trim at the idle Design boundary. No first QC or known failing
Design rewrite executes. Use existing per-process/cgroup sampling and bounded
allocator/owner helpers; no assay snapshots. Diagnostic overhead is explicit and
this member is not a timed performance comparison. aFMNxH passed the small stop/
normal-route, active-owner/idle-owner acknowledgement, duplicate-request and
browser completion/disconnection checks. Restrict traces to import and adapter
initialization to avoid per-read event floods; a recording-limit event makes the
diagnostic acknowledgement incomplete rather than silently passing.

Next bounded implementation after attribution: immutable physical table plus
logical metadata binding, unchanged-Design reuse and direct streamed transforms;
then index/group consumers and bounded cache release. One consolidated correction
pass applies; no unapproved extra large reruns. All final gates in the approved
plan remain required, not presumed complete by this instrumentation check.

**Architecture-wide read-only follow-up (2026-09-27).**
The same execution-boundedness exposure exists beyond the measured Design failure.
A bounded R fetch is not a bound on the upstream SQL sort, aggregation, temporary
relation or OS cache. These are source-inspection findings, not additional measured
300-sample failures:

- First QC (`utils_workflow_dia_state.R:1593–1670`): distinct classification keys,
  global evidence aggregation, then ordered row numbering and COPY for accepted
  data/ledger. Evidence counts are scientifically required; the persistence/order
  mechanism still has the Design path's global-sort/spill exposure.
- Grouped peptide QC/imputation (`utils_workflow_saved_table_codec.R:1156`): the
  iterator issues a whole-input grouping/join/order query before bounded group
  fetches. Complete-group numerical processing does not itself bound that SQL work.
- Rollup (`utils_workflow_dia_qc_rollup.R:280`): the persisted path materializes
  `qc_aggregated` in DuckDB before ordered final Parquet publication. Aggregation is
  necessary; temporary-output overlap and spill remain possible. The alternative
  full R collection at line 144 is a compatibility path, not the reference GUI's
  main rollup persistence path.
- Group result publication (`utils_workflow_dia_qc_rollup.R:583–639` and codec
  `:1268–1350`): unchanged data can already reuse its source, which is beneficial.
  However, identity checking uses count(DISTINCT) over the full change map; metadata
  differences prevent reuse; change maps and support products still undergo global
  ordered queries before merging/writing. Maps are narrower than assays, so this
  is a smaller exposure, not proof of equal peak cost.
- DuckDB paging and ranked selection (codec `:846`, `:929`, DIA DA `:148`, `:198`):
  LIMIT bounds the returned page, not necessarily the scanned/ordered population.
  OFFSET paging may repeat work. Ranking/selection can be required, and smaller DA
  tables/top-N optimizations mean these are lower-priority risks, not confirmed
  redundant whole-table sorts in every plan.
- Shared action cleanup (`utils_workflow_saved_entities.R:888`): spill removal and
  validated-file cache release happen at action completion. That is too late to
  limit an expensive action's peak; original import files are outside store-owned
  release coverage. The measured Design failure demonstrates this gap.

Existing Arrow sequential scanning with pre-buffering disabled, direct staged
batch writers, source/change-map merging and unchanged-assay reuse are real
improvements and should be reused. Do not replace them with another backend or
claim every SQL ORDER BY is redundant. The fix must cover the shared production,
ordering/identity and resource-lifetime contracts, prioritizing Design then
first-QC/grouped processing/rollup. Normalization/RUV/DA scientific matrix peaks
are a separate residency class. No implementation or new execution performed for
this audit question.

**Diagnosed cause of mU50nZ stop — Design rewrite/spill and cached files (2026-09-26).**
Read-only audit; no workload rerun, runtime patch, cache eviction or failure-file
cleanup performed. The final `operation_enter` is at 1790393179.45. Just before it,
combined/anonymous/file memory was 10,447,491,072 / 6,717,657,088 / 3,665,477,632 B.
The final pre-stop sample 29.168 s later was 16,602,152,960 / 7,044,444,160 /
9,469,444,096 B. File memory accounts for **94.3% of this measured growth**:
+5,803,966,464 B file versus +326,787,072 B anonymous. This is a paired-sample
comparison, not an exact decomposition of the later guard peak.

The retained head/import metadata and Parquet footer prove: 24,000,000 rows,
14 columns, `spec_tbl_df` class, parser `spec`/`problems` attributes, automatic row
names, quantities already double, `needs_numeric = character(0)`, and the same
300 sample names in the imported metadata and declared design. The current and
initial data descriptors still name the original imported Parquet. The browser
script assigns groups/replicates to all samples, performs no rename/removal, then
clicks Save Design.

The failing path is `buildProtDesignSaveResults`' saved-table transform (see
`R/mod_prot_design_builder_state_helpers.R:63`): it passes `strip_parser = TRUE`
because a numeric quantity column exists, even though no conversion is needed.
`transformDIADesignTable` changes class/attributes, so its exact-metadata no-op
check fails. It then applies ordered row_number/string row-key reconstruction and
ORDER BY to the complete assay, and `saveWorkflowTableQuery` runs DuckDB COPY.
A metadata-only compatibility change therefore causes a full 24-million-row
rewrite/order operation. This is a missed native-reference boundary, not an IQ,
RUV or S4-history computation.

Post-stop `fincore` (page-residency inspection, without reading/warming payload
contents) found all these pages resident:

| Files | Logical bytes | Approximate GiB |
| --- | ---: | ---: |
| Original OQADra report.tsv | 3,475,378,370 | 3.24 |
| Imported payload.parquet | 1,012,586,394 | 0.94 |
| Seven DuckDB spill files in the active table-spill directory | 4,804,411,392 | 4.47 |

Their page-aligned residency totals 9,292,382,208 B (~8.65 GiB), close to the
9,469,444,096 B (~8.82 GiB) last file-memory observation. This post-mortem inventory
supports ownership attribution but is not a simultaneous per-file peak census.
Spill files were written during the final query interval; source/report pages
had been cold at launch but were resident again after ingestion. The store has
scope-end resource/file cleanup, but this operation was stopped before returning.
The 1 GB DuckDB memory setting is not a cap on the job's combined R/browser/native
allocations and OS file cache. The fastduckdb projected-read path is not used by
this DuckDB COPY/ordered-transform path; qs2 cannot prevent its spill/cache cost.

Reclamation was enabled only for completed QC/normalization publications. Design
never reached one, and malloc_trim concerns free allocator pages, not these cached
file pages. Thus enabling the option could not address this observed trigger.
Anonymous memory was already ~6.26 GiB before the query and rose only ~0.30 GiB in
that interval. Its live-R/browser/native/free-page split remains unresolved in this
timed member: no process/heap census was recorded. Do not relabel it all as a leak
or as allocator retention.

Necessary next fix scope: prevent metadata-only and unchanged-row Design operations
from rewriting/reordering an assay, preserve row identity/order through the existing
codec contract, and manage consumed source-file cache at its actual ownership
boundary. Any genuine remaining large transform needs bounded spill/cache behavior;
merely changing reader names, increasing the guard or moving malloc_trim is not a
scientific/scaling solution. Earlier allocator release may reduce the baseline but
cannot be quantified from this run. No implementation or retry is authorized by
this diagnostic request; the 300-sample objective remains unmet.

**300-sample result — mU50nZ stopped by memory guard (2026-09-26).**
The authorized reclamation-enabled attempt did not complete. The controller exited
76 with `reason=early_memory_stop`; the service is now inactive. Guard measurement:
16,668,934,144 B (15.524 GiB), zero swap, host available 15,449,362,432 B. Thus the
job's unchanged 15.5 GiB early limit triggered, not the host-memory stop or timeout.
Approximately 484 s elapsed from the first sample to the final sampled interval.

The last pre-stop breakdown measured 16,602,152,960 B total, including
7,044,444,160 B anonymous (6.561 GiB) and 9,469,444,096 B file-backed/cache
(8.819 GiB). This breakdown precedes the exact guard sample and is not an exact
peak decomposition. File-backed/cache is not all guaranteed immediately reclaimable.
The latest completed state receipt remained `initial` with zero working records;
reclamation was enabled but `pending` and `last` were null. No completed release
was recorded before the stop. No first scientific checkpoint or full journey is
claimed. The stale `smoke-receipt.json` status `running` reflects interrupted
acknowledgement, not a live job; guard-stop.txt plus controller exit 76 and inactive
service are authoritative. No retry started. This is a retained limiting failure,
not achievement of the 300-sample objective.

**300-sample run launched with reclamation — mU50nZ (2026-09-26).**
The user's renewed “get it going” explicitly authorizes the prepared stale-test
correction and one preflight rerun. Applied that exact correction; **U80tqn passed**
completion/disconnection, exact-build route and reclamation propagation checks.
Launched **mU50nZ**, invocation `8ce199ac55ab4f5094634f28262cadb0`, using PCPabv,
retained OQADra and `--stress-300 --native-session --reclaim-memory`. Guards and
profiling-off settings are unchanged. This is the one requested candidate-only
attempt. Confirm actual initialized-session reclamation from its live receipt;
launch intent alone is not proof of a completed release. The preceding preflight
failures remain below as historical evidence, not current launch status.

**300 launch not started — routing preflight stopped (2026-09-26).**
XXS7pf failed because the spec validator still required an empty 300 allowlist.
The bounded correction pinned that validator to the explicitly authorized PCPabv
RDB and added an end-to-end parser assertion for `:reclaim`. utjoyR then passed
those checks but failed the older blanket rejection expectation at
`test-workflow-ruv-pearson-diagnostic.R:350`: it still expects any candidate absent
from the 96/first-QC allowlists to be forbidden at both 96 and 300, ignoring the
new explicitly authorized 300 allowlist. This is a stale test expectation, not a
300 execution or scientific/memory failure. Both failed receipts are retained.

Per AGENTS.md's one-correction limit, no additional test mutation/rerun or 300 job
has been performed. A concrete proposed patch is `/tmp/dia-300-approved-route-test.patch`:
check each sample size's existing exact-build allowlist before asserting rejection.
The positive PCPabv/300/reclamation assertion remains. Renewed authority is needed
for this narrow correction and one preflight rerun before the already authorized
300 launch. No 300 capacity attempt has been consumed. Independently recovered
puNFDj confirms the completed UbvCXM small journey's 16 scientific states match
FKZX69; full output verification remains pending.

**User-authorized immediate 300-sample attempt with reclamation — 2026-09-26.**
The user explicitly requested “do it for 300 sample and run it” with reclamation
on. This authorizes one candidate-only capacity attempt now, superseding the
previous ordering prerequisite of a fresh matched 96 pair for this attempt only.
It does not pass that prerequisite, replenish historical correction/replacement
budgets, authorize a control rerun or activate a production default. No automatic
retry is authorized by this request.

Use installed PCPabv (RDB `fd69f611d68392b6583187ef0004021f4a589dbb0cac895d7bae165cb42f99db`,
native library `69ddf7237299ec580fc4f70cf4befd2379a4a837eee3bb6eaf0d8d89fd42ca76`)
and retained OQADra input: 300 samples, 10,000 proteins, 24,000,000 precursor rows.
The explicit `--reclaim-memory` runner option passes the native-state launch flag;
default-off behavior and invalid routes are checked before launch. The 300 build
allowlist is extended only for this candidate. All guards remain 16 GiB hard,
15.5 GiB early stop, zero swap, 60 minutes, 16 GiB startup headroom and 2 GiB
continuous host stop. Profiling is disabled. Completion includes automatic
outputs, Summary and teardown; a retained limiting failure remains a valid capacity
observation, not achievement of the performance objective.

Prior work: 132VKF passed the corrected automatic-RUV failure-owner regressions;
c05IQR passed browser completion/disconnection checks. UbvCXM completed the full
24-sample GUI through Summary with reclamation off: 329.480765 s, peak
3,341,148,160 B, zero swap. Its final science/output verification was interrupted
by the user's steering and is not silently promoted to passed. No memory saving
is inferred by subtracting its later anonymous memory from its earlier peak.

**Active follow-up — automatic-RUV failure owners and Summary delivery (2026-09-26).**
The user's renewed continuation authorizes these remaining boundaries, using the
existing fixtures and guarded runner. Installed reproducer **AEwtHp failed** on
GJzhWZ: genuine automatic RUV reached both injected failure points (parameter
validation and correction), and each left raw RUV results in both GUI holders plus
a replaced global plot. All state-head/scope/resource assertions passed. Six
ownership assertions failed; 232.269627 s, peak 1,525,657,600 B, zero swap. The raw
failed receipt is retained. Its rich test diffs were excessive; assertions now
report only boolean identity failures rather than traversing plot internals.

The single bounded correction keeps parameter-discovery results in temporary local
holders and restores the caller's global plot on every exit. GUI delivery still
uses the existing converted grob/scientific-product boundary after successful
computation. No numerical method, parameter, population or output requirement is
changed. Verification: affected installed tests plus existing browser completion/
failure preflight, then one retained-24-sample native GUI journey through automatic
RUV, Pearson, DA and all three Summary actions if focused checks pass. Bind the
candidate RDB/native bytes only for the small route; no 96/300 allowlist changes.
Use unchanged small guards and preserve failures; no automatic retry beyond this
bounded correction. Independently verify all 16 scientific states and required
native/output artifacts against retained valid evidence. This is diagnostic
correctness/ownership verification, not an OFF/ON performance comparison.

**Current outcome — volcano owner repaired; ten opt-in cycles passed (2026-09-26).**
The user's next continuation authorizes the remaining included display-owner audit
and ten existing small lifecycle cycles. Read-only inspection identified static
volcano mappings created alongside the input DA result. Installed audit **WBFHgE**
failed its source-collectability assertion for the static volcano; the heatmap
source marker collected. A second failure was the audit fixture missing Glimma's
required `Protein.Ids` annotation, not a Glimma ownership finding. Preserve the
failed receipt (171.781205 s, peak 1,516,863,488 B, zero swap).

The bounded correction moves the unchanged static-volcano plotting statements into
a helper receiving only prepared display rows and scalar settings, forces all
arguments to release caller promises, and restores the previous global plot at the
public boundary. The fixture now supplies the required annotation. Candidate
**GJzhWZ** installed successfully: RDB
`7d1a935e76c5c5e7a7bab409a38cf27ad9297dd2c636c1ec75af7e9dd0d2aae6`,
native library `02bd6ac15f3258c4d372aa8f57873ab3a1d7eeffad63f9ed72e1ec3d2d6d284b`.

The existing ten-cycle fixture now uses separate Normalize and Pearson observers,
matching separate GUI actions. It opts in, checks the normalization release receipt
at delivery, checks duplicate flushes do not repeat release, and retains the native
science/archive/output assertions. Bounded global/final-QC owner summaries accompany
each cycle. Installed native-session selection **186Snt passed** in
691.542573 s, peak 3,081,101,312 B, zero swap, under 8 GiB hard / 7.5 GiB early /
900 s and unchanged host guards. This consumed the one bounded correction check;
there was no second repair/retry or large execution.

All ten normalization publications had distinct request IDs and `released` results;
extra flushes did not repeat release. Every post-cycle ownership receipt had zero
working records, native records, cached payload bytes and open resources; metadata
peaked at 2,049,664 B against the unchanged 64 MiB ceiling. The 37,218-byte cycle
receipt contains bounded observations, not allocation traces or captured objects.
All four normalization choices, manual/skipped RUV, Pearson Apply/Skip/reapply,
Revert, style rendering, native archive restoration and the existing `1e-12`
scientific comparisons passed. Explicit public-S4 oracle reads remain in the test;
the GUI operations still forbid materializing a complete S4 state.

The static-volcano source-result sentinel now collects while the returned plot
remains usable, and plotted effect sizes/q-values match the original native DA
results at `1e-12`. Heatmap and Glimma source containers also collect while their
returned displays remain usable. Glimma deliberately retains expression values
needed for linked expression display. No annotation, testing population or
statistical algorithm was changed. Affected tracked-file whitespace checks pass.

| Included owner | Evidence and boundary |
| --- | --- |
| Static volcano | Reproduced retained source result in WBFHgE; corrected data-only renderer passes finalizer and layer-value checks in 186Snt. |
| Heatmap | Source-result finalizer passes; selected-matrix equality and cached clustering/style-hit guards pass. |
| Glimma | Corrected fixture supplies its required annotation; source-result finalizer passes while the widget remains usable. Full linked expression is required display data. |
| Pearson final QC | Ten cycles inspect plot/mapping frames holding 12×9 display data (5,488 B), with existing no-refit reapply/scientific checks passing. |
| Manual/skipped RUV delivery | Post-normalization inspected global fields contain filtering summaries (4,168 B for manual RUV) or 30×4 count tables (3,584 B), not a source S4/matrix. The retained canonical display is a grob. This is bounded field inspection, not exhaustive object-graph traversal. |
| Summary/export | Code inspection routes native archives and bounded workflow metadata; ten native exports/restores and required archive dependency checks pass. Summary metadata includes the reusable RUV display/result. Full Summary-browser publication is not qualified by this fixture. |
| Revert/reset/replacement | Ten reference Revert/reset cycles pass; reuse qvuKyZ's unchanged failure/cancel/replacement and real-scope guard evidence. |

The suite retains warnings from native helper fallback publication paths, missing
plot rows and upstream RUV `aes_string()` deprecation. Required files checked by
the existing native-archive assertions pass; this is not a claim that the fixture
exercises every full-GUI output destination. Automatic-RUV transient/failure-path
plot ownership and full Summary-browser delivery still require explicit closure
before declaring the whole owner inventory complete. No new defect is inferred
from their absence in this manual/skip lifecycle fixture.

Reclamation remains off by default. No elapsed/peak reduction is claimed from this
suite, and no matched 24/96 comparison, 300-sample capacity allowance or production
activation follows. The next bounded item is closure of those remaining owner
boundaries, then same-binary OFF/ON whole-journey measurement with the previously
recorded budget reconciliation. Every historical failed/stopped member is retained.


**Current implementation result — opt-in allocator release (2026-09-26).**
Build **Xk1XZp** installed successfully. RDB SHA-256
`e9a8ef0c9d57dc714bdb627c1e3356a131554a156df6244b52563d51c79e83a7`;
package native library SHA-256
`3259f795721e076038bb57b020f4b358a60259a9f5953328240dafc95da43efb`.
The registered C interface probes support or calls Linux/glibc malloc_trim(0);
unsupported platforms return a bounded unsupported status. The reusable internal
R helper collects before release and distinguishes unsupported, released,
no-pages-released and failed outcomes. Source installation now requires a C
compiler; there is no new R dependency or runtime compilation.

Newly completed QC and normalization publications can request one post-delivery
release using their existing request/state/incarnation identities. The policy is
**off by default**, enabled only by the explicit experimental native-state
`reclaim_memory = TRUE` setting. Failed/cancelled actions, replaced/closed owners,
superseded states and active scopes, leases, working payloads or worker groups
cannot trigger release. onFlushed no longer resets a live scope counter to claim
idleness. Optional release failure does not invalidate science or published files.
The snapshot/installed identity includes package native bytes, and the guarded
installer includes the C sources. No scientific algorithm or cache limit changed.

Installed affected checks **qvuKyZ passed**: 190.482642 s, combined peak
1,545,949,184 B, zero swap. These include real native registration/release with
byte-identical live scientific data; injected unsupported/no-pages/failure paths;
actual callback scope guards; worker/resource guards; cancellation, replacement,
closure and duplicate request handling; and failed publication/recovery without
repeating science, rendering or reclamation. Existing saved QC, DA, state, export
and GUI ownership regressions pass. The warning for the injected optional-release
failure is expected. No implementation correction was needed for this suite.

The subsequent harness review normalized missing/empty native inventories for
historical pure-R receipts; compiled builds still require matching native hashes.
Installed harness/reference checks **nPMvCq passed**: 333.166867 s, combined
peak 1,741,750,272 B, zero swap. These cover the native-byte inventory, missing
compiled-library rejection, pure-R compatibility, reference-state projection and
first-QC regressions. Both installed suites completed under the existing
8 GiB hard / 7.5 GiB early / 900 s and host-memory guards. Shell syntax and
affected tracked-file whitespace checks also passed.
No capacity allowance, large allowlist or production-default activation changed.
Remaining owner audit, ten-cycle lifecycle, complete OFF/ON scientific/output
comparisons and conditional 96/300 execution remain gates. This result proves
safe optional release behavior, not a measured peak or elapsed-time improvement.

**Planning update (2026-09-26).** The user requested a plan for remaining fixes.
The proposed allocator-reclamation, ownership and capacity sequence is recorded in
[DIA-025](../../../tickets/openai/DIA-025.md#proposed-remaining-memory-and-capacity-plan--2026-09-26).
It identifies the necessary compiled-interface/build-identity extension, selective
completed-action integration, small correctness/repetition checks, matched OFF/ON
24- and 96-sample comparisons, and conditional 300-sample candidate. It proposes
an explicit supplemental execution envelope without resetting historical budgets.
This update changes documentation only; implementation, runs and large allowlists
remain unchanged. The verified current outcome follows.

**Current outcome: stock QC scientific-input retention repaired and verified
(2026-09-26).** The user's renewed continuation authorized the existing data-only
reconstruction approach. Build gKWNUR has RDB SHA-256
`01a212085dee7be6f46286f800196be906f0c4e50fafb28c396f41a5accf229e`.

`buildProtNormQcPlotProduct` now preserves the caller's last plot, extracts the
existing native producer's calculated display data and bounded metadata, discards
the native plot, and invokes `rebuildWorkflowDiaQcPlot` for the four stock QC types.
That reconstruction frame never receives the original S4, scientific matrix or fit.
The existing renderer retains native PCA tick formatting and RLE text sizing.
Numerical producers/populations, density's distinct PCA population, public S4 methods
and explicitly injected producer semantics remain unchanged. Custom injected
producers are intentionally preserved and are outside the stock-producer
reclamation guarantee; their external expressions are not silently stripped.
No generic prototype cloning remains, and no storage format or dependency changed.

Installed check **WcK017 passed**, in 100.596345 seconds with 1,541,525,504 B combined
peak and zero swap, under the existing 8 GiB / 7.5 GiB early / 900 s and host guards.
It covers the affected integration selection and existing saved-QC regression file:

- Stock PCA, RLE, density and correlation restore the prior last-plot reference;
  generation failure restores it too.
- Built layer data matches native producers at **1e-12**, and trained tick labels
  match exactly, including density's nested plots.
- The original scientific-input collection marker fires for **all four stock plot
  types while their returned product remains live and usable**. This demonstrates
  input collectability, not merely empty store counters.
- The ownership receipt shows no leftover global plot for any stock type. Returned
  environments contain calculated display dataframes/seed plots, rather than the
  original S4, scientific matrices or PCA fit. These display products are legitimate
  retained objects, not targets for blanket clearing.
- Existing scientific-product and rendered-asset hit/miss tests pass, including a
  real PNG render and no repeated scientific production on valid hits. Constant and
  all-missing input behavior, native scientific configuration and affected legacy/
  publication checks pass. External-expression callback compatibility also passes.

Evidence is retained in `browser-smoke.WcK017/reader-tests.rds`,
`qc-global-owner-audit.json`, `service.log` and its terminal receipt. This resumed
repair passed its first installed check; no additional correction or measurement
run was needed. The earlier failed cloning candidates below remain withdrawn and
preserved. The verified composite last-plot repair remains intact.

**Limits:** no MiB reduction, complete latest-build workflow peak, or 300-sample
capacity is inferred from this small ownership proof. Generic injected callback
lifetimes and other plot families have not been certified by these checks.
Allocator free-page retention is still a separate requirement: this repair releases
live references but adds no production GC/trim policy or native binding. No full GUI,
96/300 trial, formal qualification or production-default activation was performed.
The bounded QC ownership repair is complete; further reclamation work requires a
new measured target and scoped check rather than reopening this successful repair.

### Previous ownership audit and withdrawn attempt (retained)

**Previous ownership audit: additional scientific roots confirmed; attempted repair
failed and was withdrawn (2026-09-25).** The user's latest request authorized auditing
and reclaiming unneeded memory. Audit Ysd0N1 reproduced an exact original
`ProteinQuantitativeData` retained by the global correlation plot. Audit 166NBD
extended this to returned plot mappings and formatter environments:

- PCA's global and returned aesthetic mappings retain the full input matrix,
  CV-filtered matrix and PCA fit.
- Density's mappings retain an upstream PCA plot.
- Correlation's global plot and returned mapping retain the exact original S4 input.
- Returned `plot_env` is already baseenv; that shallow check did **not** establish
  that the rest of the plot's environment graph was free of scientific inputs.

Both audits used the existing 60-feature/8-sample saved-QC fixture, logging only
classes, dimensions and identity flags. No unbounded object dump was taken.
The existing integration selection now also includes the affected saved-QC file.
Saved QC products already omit the native plot on disk; their reference-path display
is rebuilt from bounded data. The new findings concern live producer/global and
returned plotting references, not a newly discovered full-assay file cache.

A proposed cleanup preserved the prior last plot and detached data-only mappings/
formatters while cloning affected ggproto inheritance. Initial build iOSrBz installed
(RDB `6dc990517d3894ff8001487316f822e8375a9ba7639803695ff1c227c4af0cdc`), but hiTvF0
failed native plot building: `.class`/`.inherit` were not ggplot2's constructor
argument names and inherited `clone()` was lost. The bounded correction switched
to positional arguments and baseenv construction to avoid capturing cleanup frames.

Corrected build fOjk2m installed, but **gyN6c1 still failed**. Prototype fields include
recorded language calls; `do.call` evaluated one as code instead of preserving it as
data, causing `scale_x_continuous` lookup errors. The stronger collection-marker
test also failed: the retained correlation product did not permit collection of its
input marker. Fixing constructor spelling was therefore not sufficient to establish
reclamation or rendering parity. No numeric tolerance or output criterion was relaxed.

| Retained member | Outcome | Combined peak bytes | Swap |
| --- | --- | ---: | ---: |
| Ysd0N1 | Observation audit passed; global S4 owner reproduced | 1,538,826,240 | 0 |
| 166NBD | Extended observation audit passed; mapping/fit owners reproduced | 1,523,851,264 | 0 |
| hiTvF0 | Failed first cleanup verification | 1,563,721,728 | 0 |
| gyN6c1 | Failed corrected cleanup verification | 1,555,648,512 | 0 |

The failed candidate sources and logs remain in their installed proof snapshots.
The new regression tests remain to express the missing ownership guarantees.
**The unaccepted runtime cleanup was removed**, restoring
`R/mod_prot_norm_qc_generation_helpers.R` byte-for-byte to verified lYu5vC. The earlier
verified composite last-plot fix remains intact. No new allocator policy, native
package code, full GUI or 96/300 run was introduced. All members kept the existing
8 GiB / 7.5 GiB early / zero swap / 900 s and host guards.

Per AGENTS.md, work stops after the one bounded correction still fails. The next
repair should prefer the existing data-only reconstruction path over general
prototype/environment cloning, while preserving native plot data/tick labels,
external-expression semantics and public S4 dispatch. That requires renewed scope
for a concrete repair/check pass. This audit discovers real additional roots but
**does not claim those roots have been reclaimed**. Prior allocator-cost and
composite-fix evidence below remain valid and unchanged.

### Prior allocator-release comparison (retained)

**Completed-boundary release comparison complete (2026-09-25).** On unchanged
lYu5vC and retained FKZX69 images, one diagnostic control/treatment pair measures a
**367.270 MiB anonymous-memory reduction in a 6 ms release call** between renders.
The following render shows no observed slowdown in this pair, but the sample is
not sufficient for a performance nonregression claim. Both processes return to
nearly the same residency after the next render, and their process high-water marks
are essentially unchanged. This supports transient between-action headroom, not a
solution to rendering's active peak or 300-sample capacity. No production policy,
native package interface, dependency, process topology or allocator tuning changed.

### Design and measured outcome

The existing gui-diagnostic fixture uses two fresh, sequential callr children,
control first then treatment, with the same installed build and original images.
Each calls the production composite builder/saver twice, including automatic
PNG/PDF/RDS output. Both get identical baseline GC/trim and post-render full GC.
Only the operation between render 1 and render 2 differs: allocator counters alone
versus counters plus `malloc_trim(0)`. Timers exclude compilation, package loading,
verification and the separate GC observations; they include complete render/export
calls. No allocation or CPU profiler runs. Millisecond timer resolution limits the
precision of the short release-call measurement.

| Observation | Control | One between-action trim |
| --- | ---: | ---: |
| First render/export elapsed | 26.680 s | 21.181 s |
| Between-action operation elapsed | 0.001 s | 0.006 s |
| Anonymous RSS before operation | 1,257,213,952 B | 1,257,164,800 B |
| Anonymous RSS after operation | 1,257,213,952 B | 872,054,784 B |
| Immediate second render/export elapsed | 18.709 s | 17.857 s |
| Anonymous RSS after second render/full GC | 1,653,485,568 B | 1,653,424,128 B |
| Live R after second render/full GC | 582,171,184 B | 582,171,176 B |
| Observed child process high-water RSS | 1886.262 MiB | 1886.574 MiB |

The treatment's within-process reduction is 385,110,016 B. The release call's
increment over the counter-only control is approximately 5 ms. Do **not** label
the faster second-render observation a speedup caused by trimming: the first
renders already differ substantially *before* the intervention, and there is only
one sequential pair. This tests immediate reuse of the same image workload, not
DA/RUV or a different next scientific action. The post-second-render convergence
also demonstrates that a one-off release does not prevent residency accumulating
again. Repeated completed-boundary use and whole-workflow cost remain unmeasured.

### Output verification and retained verifier failure

Raw member **0dj9B1 remains failed verification**. Its first oracle incorrectly
required both iterations' RDS bytes to equal a fresh process's first-render
reference. Both child processes and all four renders completed and their bounded
measurements were saved before that assertion failed. The overall failed-member
interval is not accepted as a workflow timing result; the completed per-action
measurements above remain explicit diagnostic observations.

All four PNGs are byte-identical to 6EqMOx's retained PNG. Each first-render RDS is
also byte-identical to the retained RDS. The two second-render RDS files match each
other exactly but differ from the first-render RDS, including in the **control**.
No claim about the internal cause of that iteration-dependent serialization is
made here; it is not attributed to allocator release. The strict relevant oracle
is matching execution context: compare corresponding iterations byte-for-byte,
retain both first-render reference checks, and compare every PNG with the retained
reference. No numerical tolerance or required output was relaxed.

Read-only verifier **iMrdOH passed** those checks against the existing artifacts and
confirmed that the original receipt remained unchanged (SHA-256
`10835792a7ed3dbfa75cd7615d79ab6bce8e0011a38410988b9eb2e9a3aa3309`). No measured arm
was rerun. Both source-image hash inventories remain unchanged, and required PDFs
are present/nonempty. Exact matched PNG/RDS hashes:

- All PNGs: `d41c5063d2246cb60b4c6c57327a8f3279cf73d7da22a4c1ba9dcbb9b365c7de`.
- First-render RDS, both arms and reference:
  `c24a682183650d56f2576cae9de44a65794b1126efdb97fc4abbb88876b398cb`.
- Second-render RDS, both arms:
  `fec3400961d7606b933be2cb8e0481c91bf1f88aab64b41de1430653611f9129`.

Evidence: `browser-smoke.0dj9B1/allocator-release-comparison.json`, both child phase
logs/output directories, and `browser-smoke.iMrdOH/allocator-output-verification.json`.
The existing lifecycle/failure fixtures precede the comparison. 0dj9B1's combined
suite peak is 4,375,695,360 B, zero swap; it includes worker/browser preflights, not
just a rendering child. iMrdOH completed in 14.931501 s at 1,951,322,112 B peak,
zero swap. Both kept the existing 8 GiB hard / 7.5 GiB early / 900 s and host guards.
The one declared correction was used for the output oracle; every failure is kept.

This bounded cost question is answered: explicit release can recover substantial
between-action residency cheaply in this diagnostic, but does not remove the next
action's live demand or establish a general production policy. No fresh full GUI,
96/300, formal qualification or default activation was performed. The next decision
is whether to develop and qualify an explicitly scoped completed-action release
interface; a diagnostic Rcpp compilation is not a production implementation.

### Prior allocator-retention audit (retained)

**Allocator-retention audit complete (2026-09-25): the measured residual is
principally free interior pages in the main heap.** No production code or allocator
policy changed. The repaired lYu5vC binary and original FKZX69 image bytes were
reused. One existing guarded `gui-diagnostic` member, b8UwdY, passed without retry;
its native-counter/lifetime/failure checks precede the retained production-image
replay. No new full-GUI or 96/300 execution occurred.

The replay still has **zero last-plot retention**, byte-identical PNG output and
unchanged source hashes. At completion it has none of `arrow`, `duckdb`,
`fastduckdb` or `qs2` loaded. Thus this graphics-path residual can be reproduced
without those reader/serialization engines; this is not a claim about every
allocation elsewhere in the full workflow.

### Arena and resident-mapping observations

| After complete figure export and repeated full GC | Observed |
| --- | ---: |
| Live R heap | 582,856,184 B |
| Total anonymous RSS | 1,181,704,192 B |
| glibc free arena bytes (`mallinfo2`) | 351,272,320 B (334.999 MiB) |
| Main arena free bytes (`malloc_info`) | 350,172,264 B (333.950 MiB) |
| Free space at the main heap's top (`keepcost`) | 269,152 B (0.257 MiB) |
| Arena count | 2 |
| Secondary arena total / free bytes | 1,130,496 / 52,226 B |
| Main-heap mapping resident bytes | 1,171,939,328 B |

Diagnostic-only trim reduces anonymous RSS by **332,697,600 B (317.285 MiB)**.
Of that change, **332,693,504 B (317.281 MiB)** is observed in the `[heap]` mapping.
Its virtual allocation shrinks by only **266,240 B (0.254 MiB)**; other anonymous
mappings change only by KiB. Live R changes by -1,456 B. This establishes that the
large residency drop is chiefly page reclamation *inside* the existing main-heap
mapping, not disappearance of a large live R object, a large mmapped allocation,
or removal of many thread arenas. The secondary arena is only 1.078 MiB total.
Reducing the arena count cannot explain away this ~317 MiB observed main-heap cost.

Interpretation: freed blocks lie below live allocations rather than forming one
large free heap tip. Ordinary freeing/top trimming and explicit whole-free-page
reclamation are different mechanisms. Current upstream glibc's
[`_int_free_maybe_trim` and `mtrim`](https://github.com/bminor/glibc/blob/master/malloc/malloc.c)
show the top-of-heap versus interior-page distinction; the installed replay's
mapping/allocator measurements above supply the actual evidence here. A lower
trim threshold alone is not demonstrated as a remedy for these interior pages.
The audit does not identify which individual allocating function created each hole,
prove an allocation-rate model, or establish a cost-free production release policy.

Measurement bounds: allocator XML uses a fixed 1-MiB buffer and logs only per-arena
summaries; `/proc/self/smaps` is read in 256-line batches with a 100,000-line limit,
logging four aggregate mapping categories without addresses or file names. XML
capture therefore perturbs allocation by up to its bounded scratch buffer; its free
bytes need not exactly equal the earlier `mallinfo2` snapshot. Resident snapshots
are also non-atomic. The native arena-byte total exactly matches the summed XML
arena totals at every recorded boundary. The installed manual contains a main-only
`mallinfo2` caveat, but these runtime cross-checks and current upstream source show
all-arena aggregation for these observations; do not treat the two counters as
independent memory pools or add them to RSS.

### Existing API and cleanup audit

- Pinned fastduckdb 0.1.0 exports `clear_parquet_cache`; its retained Rust source
  clears the 64-entry metadata cache (`PLAN_CACHE.lock().unwrap().clear()`), without
  allocator page-release logic. Its retained Cargo.lock SHA matches the runner's
  pinned binding `0dfffbc58e67f7c3eff793b6b7fdbe17ebe8733838e52297aebeef02dc02e073`.
  `readSavedWorkflowTable` already invokes that API in `finally` after fast reads.
- `releaseWorkflowTableReadScope` shuts down DuckDB on release of the last scope;
  `endWorkflowTableAction` clears the connection/spill/validation ownership.
  `closeDIAStateAdapter` clears its working roots, workers and callbacks, then closes
  the store. Code inspection found no missing ordinary close at these boundaries;
  it is not a general proof that all native allocations in every path are freed.
- Arrow writer completion closes writer/sink and drops pending tables; GC controls
  finalization of remaining private wrappers. The pinned Arrow 24 R memory-pool API
  exposes backend/live/peak counters, not a general page-release operation.
- qs2's exported `qs_cache(clear=TRUE)` removes cached files. It is not a memory
  reclamation API and was **not invoked**. Deleting those files would not target the
  measured graphics-only residual.
- `finalizeProtNormWorkflowState` already calls GC after clearing its GUI composite
  holder; the last-plot fix now removes the separate global root as well. In this
  replay, another full GC did not lower the residual RSS; explicit page reclamation
  did. Adding another GC alone is not a demonstrated solution to this residue.

Result/receipt: `browser-smoke.b8UwdY/anonymous-memory-composite/attribution.json`
and its bounded phase log. The whole diagnostic suite took 282.581934 seconds,
4,371,091,456 B combined peak, zero swap, under unchanged 8 GiB / 7.5 GiB early /
900 s and host guards. That peak includes browser/worker fixtures and is not the
isolated graphics peak. Production source remains identical to lYu5vC; the GUI
allocator opt-in remains false. Earlier failed evidence and experiment budgets
are preserved. This completes the requested bounded discovery/audit; a prospective
production release strategy still needs an explicitly bounded cost/benefit test at
an existing completed-action boundary. Whole-workflow peak and 300-sample capacity
remain unqualified.

### Previously verified composite-lifetime repair

**Current outcome: composite last-plot retention fixed and verified (2026-09-25).**
The user's “Continue” authorized the scoped repair of the measured plot owner.
`generateProtNormCompositeQcFigure()` now saves the caller's ggplot2 last-plot
reference and restores it with `on.exit`, including generation/save failure paths.
This removes the export's separate package-global root while preserving a prior
user plot or NULL. It adds no production GC, allocator trim, native code or dependency.
The only runtime source difference from BYmpMh is four added lines in
`R/mod_prot_norm_workflow_helpers.R`; no other R source file changed in this build.

Current build lYu5vC RDB SHA-256:
`e3cd18670bc0af614f3f8db936e835493aeb019db07881d39d66f8cf3d14ba26`.

| Verification | Outcome | Service-tree seconds | Combined peak bytes |
| --- | --- | ---: | ---: |
| EO6Gvu, affected integration | Pass: normalization/publication consumers; preserves prior NULL/plot on success, generation failure and save failure | 87.371827 | 1,576,062,976 |
| oe2rfb, retained production-image replay and diagnostic fixtures | Pass: no composite last-plot retention after return; unchanged images, required files, source hashes, allocator/lifecycle/failure observations | 290.585305 | 4,394,061,824 |

Both use existing small guards, zero swap, and the established runner. The
instrumentation's completion/failure and lifecycle fixtures precede the replay.
No correction/retry was needed. Diagnostic-suite peak includes worker/browser
fixtures and must not be called the isolated composite peak.

The replay uses the same FKZX69 images and production builder/saver as pre-fix
6EqMOx. Its last-plot retained size is **0 bytes**, versus **572,298,400 B** before.
Immediately after return and full GC, live R heap is **582,026,712 B (555.064 MiB)**,
versus **1,150,943,120 B (1097.625 MiB)** before: **568,916,408 B (542.561 MiB)** less.
The reduction occurs before the diagnostic's manual last-plot clear or trim.
Anonymous RSS at that boundary is 1,206,546,432 B versus 1,591,009,280 B, a
366.652 MiB observed reduction. This is isolated post-return residency evidence,
not a whole-GUI timing/peak or 300-sample claim.

Required composite outputs are preserved:

- PNG: byte-identical, 192,000 B, SHA-256
  `d41c5063d2246cb60b4c6c57327a8f3279cf73d7da22a4c1ba9dcbb9b365c7de`.
- RDS: byte-identical, 24,077,550 B, SHA-256
  `c24a682183650d56f2576cae9de44a65794b1126efdb97fc4abbb88876b398cb`.
- PDF: both 430,331 B; all differing bytes are in one compressed object stream.
  Decompressing that stream shows only the expected `/CreationDate` change;
  page structure, dimensions and remaining bytes are unchanged.
- Original source image hashes remain unchanged. The independently saved output
  remains available even though the temporary global display reference is restored.

This bounded repair's exit conditions are met. **Allocator residency remains a
separate unresolved performance requirement.** The diagnostic still observes freed
pages resident after GC; no blanket trim policy is introduced. Rendering's active
peak, a complete latest-build journey, matched performance qualification, 96/300
capacity and production-default activation remain unqualified/held. No new full GUI
or capacity member was run for this four-line repair, and no other-omics migration
or generic plotting framework was added. The previous measured failures and
diagnostic outcomes below remain historical evidence.

### Diagnosis preceding the plot-lifetime repair (retained)

**Continued diagnosis complete: a missed live-plot owner and reclaimable glibc
pages are both demonstrated (2026-09-25).** The user explicitly required continuing
past the isolated first-grid replay. Production source and installed BYmpMh binary
remain unchanged; only the existing diagnostic harness/fixtures were extended.
The opt-in is disabled again. No production allocator policy, scientific method,
dependency, 96/300 trial or default activation changed.

### Identified live owner: ggplot2's global last-plot store

`generateProtNormCompositeQcFigure()` constructs the composite through
`generateProtNormCompositeFromFiles()` and saves it through `savePlot()`. The
composite contains raster panels loaded by `loadProtNormImageAsPlot()`. The local
composite function returns and `finalizeProtNormWorkflowState()` clears
`normData$QC_composite_figure`, but **ggplot2's namespace-global `.store` still holds
the last plot**. `ggplot2::get_last_plot()` exposes this separate root; clearing
only the GUI holder does not release it. The earlier GUI-holder inventory missed
this package-global owner.

The real GUI probe 0ZcHHS observes 731.58 MiB live R after protein QC/full GC,
1281.24 MiB after normalization/full GC, and a later fall exceeding 500 MiB.
Focused replay 6EqMOx runs the exact production composite builder and saver with
FKZX69's retained source images, writing required PNG/PDF/RDS to a fresh directory.
All source-image hashes remain unchanged. After the production function returns,
`get_last_plot()` retains a **572,298,400-byte (545.786 MiB) patchwork object**.
Clearing only that root in the disposable diagnostic child and collecting releases
**568,761,232 B (542.413 MiB) of live R heap**. This identifies a retaining owner,
not merely a correlation with RSS.

| Production composite replay boundary | Live R MiB | Anonymous RSS MiB |
| --- | ---: | ---: |
| Before construction, GC and baseline trim | 540.802 | 767.352 |
| Builder/saver returned, full GC | 1097.625 | 1517.305 |
| Only last-plot root cleared, full GC | 555.212 | 1187.738 |
| Following diagnostic allocator trim | 555.212 | 801.211 |

The root clear releases 329.566 MiB anonymous RSS immediately. The following trim
returns another **405,303,296 B (386.527 MiB)** with live R changing by only 416 B.
Thus object lifetime and allocator residency are distinct problems. The residual
post-trim increase over the loaded baseline is ~14.410 MiB live R and ~33.859 MiB
anonymous memory; it is not assigned to a specific cache or asserted to be a leak.
The replay's source hashes and owner class/size are in
`browser-smoke.6EqMOx/anonymous-memory-composite/attribution.json`; its bounded
phase log records the interventions and measurements. No heap/object dump was taken.

### Identified native residency: free glibc arena pages

In the real GUI process (0ZcHHS), the first Q-value figure action increases anonymous
RSS from 1,123,741,696 to 1,370,730,496 B, while live R after full GC changes from
677,360,656 to 701,866,752 B (+23.371 MiB). At delivery, glibc reports **326,570,800 B
(311.442 MiB) free arena space**. A diagnostic-only `malloc_trim(0)` is followed by
anonymous RSS of 1,153,933,312 B: a **216,797,184 B (206.754 MiB)** reduction before
another scientific action. Glibc in-use arena bytes are unchanged across the trim
(979,758,288 B). Arrow live pool allocation is zero at these delivery observations.
This demonstrates reclaimable glibc pages in the real application. The reclaimed
pages need not all have originated in this one action: baseline free arenas already
existed, so the release must not be equated with an exact allocation census of the
235.547 MiB action-local RSS rise.

Later QC actions refill arenas and residency rises again. One successful trim does
not establish a lasting or cost-free production fix. The composite replay above
also shows why GC alone is insufficient: first a live root must expire; subsequently
freed allocator pages can still remain resident.

Instrumentation limitation preserved: in 0ZcHHS, R's lazy evaluation caused the
trim row's RSS to be sampled before its native-counter argument executed. The next
unchanged-state `flush:owners` row, **35.17 ms later**, records the post-trim RSS.
Both raw rows remain retained. The helper now evaluates the native operation
before memory recording; J8NVBC verifies that sequencing and passive defaults.
No full GUI run was discarded or repeated for this correction.

### Verification, budgets and conclusion

| Member | Purpose/outcome | Service-tree seconds | Combined peak bytes |
| --- | --- | ---: | ---: |
| s8lAHm | Pass: real allocator/GC/lifetime and intentional-failure preflight | 237.162078 | 4,377,583,616 |
| 0ZcHHS | Complete retained 24-sample GUI with bounded allocator observations and one Q-value-delivery trim | 443.391412 | 3,242,962,944 |
| 6EqMOx | Pass: retained production composite replay identifies last-plot root; immutable source hashes | 320.893890 | 4,470,566,912 |
| J8NVBC | Pass: passive defaults and native-release-before-memory-record sequencing | 33.227649 | 2,007,339,008 |
| 20GaMT | Pass: all 16 saved scientific states unchanged versus FKZX69 at 1e-12, read-only comparison | 51.286727 | 1,347,858,432 |

All members retained 8 GiB hard / 7.5 GiB early / zero swap / 900 seconds and existing
host guards. Peaks for diagnostic suites include their worker/browser fixtures,
not just the isolated plot. GUI uses unchanged BYmpMh RDB SHA-256
`a91ed1e91c743d8d4734ead86573827e5ef619f84574c5e76246f07dc9e49ef4` and retained ykIZwA
24-sample bytes. The extra composite replay was justified by the newly measured
~550 MiB live-R increase; it was not a renamed retry or another GUI journey.
One instrumentation ordering correction was used. The small-route probe is now
disabled, and all prior failures/allowances remain recorded.

**Diagnosis:** the sustained memory has a concrete missed plot owner plus allocator
retention; it cannot be summarized as only live assays, only file cache or only
freed pages. These measurements explain major observed residency, but do not
attribute every byte, prove which earlier change altered allocator layout, or
qualify 300-sample computational peaks. A future fix should manage the composite's
last-plot lifetime while preserving any unrelated prior plot, then evaluate native
residency separately. No unconditional production trim, new backend, dependency,
process topology or further optimization is authorized by this diagnosis alone.

### Earlier isolated-grid diagnosis (retained)

**Anonymous-memory diagnosis, 2026-09-25: one source demonstrated; full GUI gap
not yet attributed.** User explicitly selected persistent anonymous memory as the
next diagnostic target. No production code, allocator policy or scientific
execution was changed. No full workflow or 96/300 member was rerun.

Retained full-journey evidence locates the first substantial rise around Q-value
figure publication: anonymous memory is 1,078.0 -> 1,317.0 MiB in FKZX69
(+239.0 MiB), versus 1,143.9 -> 1,178.6 MiB in lORe8g (+34.7 MiB).
These are recorded boundary observations; they do not alone identify allocations.

One guarded existing `gui-diagnostic` member cOBFQP passed its completion/failure,
live-heap and browser fixtures, then replayed the first retained FKZX69 display in
an isolated child. Source: 20,760-byte `.workflow-export/qvalue-a277baf727685fa6432f382a889f20fdbb42cbcd/display.rds`,
SHA-256 `4d31ab4789b0774a0c8204a0f0a02014b0c0c3c0e34c8084368cd0428a476551`.
The source hash was unchanged after replay. PNG/PDF export used the original
15 x 18 inch / 300 dpi settings. ggplot2/ragg namespaces were loaded before the
baseline. The test-only counter uses existing Rcpp and glibc `mallinfo2`; no package
dependency or production helper was introduced. Each observation is bounded
scalar metadata, not an allocation log or object capture.

| Replay boundary | Live R MiB | RSS MiB | Anonymous MiB | glibc in-use arena MiB | glibc free arena MiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| Loaded display, after baseline GC/trim | 68.735 | 131.723 | 107.445 | 100.308 | 2.867 |
| PNG complete, full GC | 69.905 | 177.156 | 147.773 | 138.442 | 3.328 |
| PDF complete, full GC | 70.224 | 177.285 | 147.773 | 108.096 | 33.674 |
| Display released, full GC | 69.233 | 177.285 | 147.773 | 108.052 | 33.718 |
| Diagnostic-only allocator trim | 69.234 | 145.516 | 116.004 | 107.307 | 33.880 |

**Established:** after figure release and full GC, glibc has 35,355,408 bytes of
free arena space while anonymous RSS remains elevated. Diagnostic-only
`malloc_trim(0)` releases 33,312,768 resident bytes (31.769 MiB), with live R heap
changing by only 936 bytes. Thus reclaimable glibc allocator pages, not live R
objects, account for that measured part of the replay's residency. Free arena
address space can remain recorded after its physical pages are reclaimed; the
free-arena counter is not itself resident memory. The post-trim anonymous excess
over baseline is still 8.559 MiB and is not assigned to a specific owner.

**Limit:** this replays the saved combined grid, not the five individual plots,
plot construction/publication, prior scientific allocations or the full Shiny
process. It reproduces 40.328 MiB of post-GC anonymous growth, not the complete
239 MiB full-GUI rise. It proves allocator retention is one real source, not that
all excess RSS is glibc, that there is no leak, or that routine trimming is the fix.
The full-journey Arrow observations still exclude large live Arrow pool buffers
at the composite boundaries, but do not identify other native allocations.

Receipt: `browser-smoke.cOBFQP/anonymous-memory-figure/attribution.json` plus its
bounded phase log. The whole diagnostic suite completed in 254.819702 seconds,
4,395,782,144 B combined peak, zero swap, under existing 8 GiB/7.5 GiB early/900 s
and host guards. That suite peak includes worker/browser fixtures and is **not**
the isolated figure's peak. One member was run; no correction/retry was needed.

This bounded diagnosis is complete. The next discriminating measurement is the
same live-R/glibc/RSS accounting in the real process immediately before/after the
complete Q-value figure publication action, with GC and a diagnostic-only release
at a completed boundary. It must precede any production allocator change or claim
that the whole residual-memory problem is explained. SQL late-load instrumentation
is unrelated to this allocation measurement and remains untouched.

**Current outcome (2026-09-25): publication fixed; complete small science/output
verification passes; performance objective remains unmet.** The user's “CONTINUE”
explicitly resumed the stopped correction. `snapshotWorkflowFile` now distinguishes
an absent destination from a non-NA symlink target. Existing files and dangling
symlinks remain rejected; failed replacement preserves the prior completed output.
No dependency, backend, numerical rule, public S4 API or required output was changed.

Current build BYmpMh has RDB SHA-256
`a91ed1e91c743d8d4734ead86573827e5ef619f84574c5e76246f07dc9e49ef4`.
All members use the existing guarded runner and retained inputs. The small journey
FKZX69 used ykIZwA report/FASTA/design (24 samples, 192,000 precursor rows; report
SHA-256 `ce7868fb0873e32df845694e6b10f2ca421495d2714f87a6d49529ecf7037e51`).
Limits remained 8 GiB hard, 7.5 GiB early stop, zero swap, 900 seconds, 10 GiB startup
headroom and 2 GiB continuous host stop. Every completed member below had zero swap.
No 96/300 run, formal qualification or production-default activation was performed.

| Evidence | Outcome | Service-tree seconds | Combined peak bytes |
| --- | --- | ---: | ---: |
| 4TCw0Z, f9hubU core | Pass: exact codec/scientific parity, complete peptide QC, identity maps, row groups, cancellation; code unchanged by destination fix | 81.642120 | 1,626,009,600 |
| h4rSOX, diagnostic lifecycle | Pass: completion/failure, blocking browser action, pool fields, SQL dispatch with an already-loaded namespace; instrumentation unchanged | 244.826988 | 4,394,725,376 |
| 8QrUXA, BYmpMh integration | Pass: publication/recovery, source independence, fresh/existing/dangling destinations, affected scientific/legacy consumers | 154.952359 | 1,553,178,624 |
| 8xhb2A, routing | Pass: real browser completion/failure and rejection of BYmpMh for 96/300 | 38.551532 | 1,362,190,336 |
| FKZX69, full small GUI | Complete: all 16 states, native publication, DA, Summary and teardown | 519.632806 | 3,257,507,840 |
| bjtM3d, independent history | Pass: all 16 states against lORe8g at 1e-12, exact structural/provenance contract, read-only source comparison | 79.197261 | 1,454,596,096 |
| 6od1mR, independent required outputs | Pass: native restore/integrity, IQ/QC/DA/Summary inventory and real RUV/Pearson execution; producer library BYmpMh | 57.543456 | 1,714,065,408 |

FKZX69 has zero full-S4 GUI reads, zero open store resources at delivery, and no
resident assay in the inventoried GUI holders. This does not certify arbitrary
unobserved closures or replace the independent process-memory measurements.

### Resource outcome and attribution limits

| Complete small member | Seconds | Combined peak bytes | IOReadBytes | IOWriteBytes |
| --- | ---: | ---: | ---: | ---: |
| lORe8g, frozen native-first control | 448.139854 | 2,971,516,928 | 2,260,189,184 | 2,142,322,688 |
| kihu2x, prior shared-I/O repair | 422.919151 | 3,317,284,864 | 2,222,022,656 | 2,106,556,416 |
| FKZX69, audit repairs | 519.632806 | 3,257,507,840 | 1,728,258,048 | 1,236,918,272 |

Against kihu2x, observed reads/writes fall 22.22%/41.28%, peak falls only 1.80%,
and elapsed rises 22.87%. Against lORe8g, peak is still 9.62% higher and elapsed
15.95% higher. This is **not performance acceptance**. FKZX69 adds bounded Arrow
pool/query diagnostic fields; earlier members lack those fields. There is no fresh
matched timing pair or statistical qualification, and no causal attribution of the
whole elapsed change. Profilers were disabled in all three diagnostic members.

Concrete mechanism evidence from retained descriptors, Parquet footers and phases:

- Rollup duration is 34.211 -> 2.910 seconds against kihu2x. Annotation production
  uses one projected DISTINCT-pair cursor; store query-stream count increases from
  6 to 7 because this previously uncounted operation is now a cursor. This counter
  is not the total SQL-statement count.
- Precursor rollup and intensity/sample/replicate states now share exactly one
  assay entity (`82b3014dc9b2a5cfab9a62d312a7aa670e30417190e827d88054eedf41cee2a3`),
  2,449,517 bytes. The prior repair wrote three distinct 3,528,081-byte copies.
- Evidence-table row groups fall 20 -> 2, but that file grows 9,298,920 -> 9,479,653
  bytes: larger row groups are not a universal compression win. Imputation falls
  10 -> 1 row groups and 4,469,865 -> 4,255,985 bytes; it remains larger than the
  original control's 3,219,362-byte output. Types/values remain independently equal.
- Retained logical file lengths fall 859,341,891 -> 800,574,512 bytes; archive lengths
  fall 467,077,120 -> 424,304,640 bytes. These are logical lengths, not physical
  allocation on a copy-on-write filesystem or process I/O counters.
- Store staged/output bytes fall 13,408,325/141,472,084 -> 9,594,422/125,486,906.
  Group-buffer row/byte accounting is now included, so old materialization counters
  must not be treated as comparable total-allocation measurements.

At composite entry, main-process RSS is 1,598,464,000 B; at composite `savePlot`
entry it is 2,246,414,336 B. Arrow live allocation is **zero at both boundaries**,
with a lifetime pool peak of 53,888,192 B (mimalloc). The largest sampled cgroup
current has anon/file 2,927,906,816/195,145,728 B, versus prior repair
3,009,347,584/189,161,472 B. These are non-atomic samples, not exact peak ownership.
Existing grid `gc()` used-heap observations start at 662.8 MiB (control 653.6,
prior repair 662.2) and end at 744.2 MiB (734.9/743.8). Live Arrow buffers alone
cannot explain the residual RSS. R/native allocator retention and other native
allocations remain unresolved; these observations do not identify an allocator
or prove a retained whole assay. Normalization and composite assembly elapsed
increase from 97.238/8.284 to 140.893/16.912 seconds against the prior repair.

Instrumentation gap: DuckDB loaded after `installCheckpointGuiTraces`, so SQL
tracing was not installed in FKZX69. Its zero SQL fields are unavailable data,
not zero SQL execution. h4rSOX covered an already-loaded DuckDB namespace only.
Arrow-pool, phase, owner and independent cgroup evidence remain valid. No whole
member was discarded or repeated to obtain a preferred outcome.

The resumed bounded correction and remaining small verification are complete.
Do not start another optimization or large trial from these results alone. The next
scope decision concerns the measured residual anonymous memory and normalization/
plotting cost, plus SQL trace loading-order coverage. Native-first scientific
correctness and reduced I/O do not establish the requested 300-sample capacity.

### Retained failed audit-repair attempts and authority

The initial “okay - fix it” authorized row-group/cancellation, repeated rollup scan,
confidence-window, no-op table and native-copy repairs, plus bounded diagnostics.
Build ln15FB installed; core aWaeb4 failed before table production because manual
`.unsafe_delete()` calls are unsupported through ArrowTabular's `$` accessor.
Its peak was 1,255,985,152 B, zero swap. The declared correction removed manual
Table/RecordBatch disposal; writer references are dropped normally and private
reader cancellation retains its supported explicit disposal.

Corrected build f9hubU (`c4e6fd6e8270c80f8b54fc64080ccf2c82b1e7b63ed0ca260b2b260c659e05f4`)
passed core but integration NiwoEp failed two copy checks: `Sys.readlink()` returns
NA for absent destinations and `nzchar(NA)` is TRUE by default. NiwoEp recorded
131.875439 seconds, 1,529,733,120 B peak, zero swap. Work stopped at the declared
correction limit. The user's subsequent “CONTINUE” authorized the one path-check
correction and affected checks, then the previously unexecuted small journey.
No historical failure or large-run allowance was erased or reset. BYmpMh was added
only to the small native-session build list; 96/300 allowlists remain unchanged.

**User-directed implementation-first continuation.** The user explicitly requires
finishing the entire DIA IQ reference-backed refactor before any more testing.
This supersedes the slice-B implementation stop and inter-slice test ordering.
No tests, rebuilds, profiling or capacity runs are to start until C–E and all
included consumers are implemented. Preserve numerical/public S4 compatibility,
all previous failed evidence, unchanged resource limits and the existing large-run
budget; this instruction does not erase or replenish experiment allowances.
The complete included refactor is now implemented and the full small journey
passes strict science and required-output checks. UoWcXR remains a retained
first-QC capacity failure, stopped by `early_memory_stop` at 16,748,548,096 B,
zero swap, before an assay component was committed. No full-size performance
qualification or production activation is implied.

Previous outcome (DIA-025 shared I/O repair, dbDWgV): **scientific/output verification
passes; performance acceptance fails.** The complete 24-sample repaired journey
matches all 16 scientific states and passes independent native/IQ/QC/DA/Summary
output and real RUV/Pearson execution checks. Core ownership/integrity tests, ten
normalize/Revert/style/export cycles, fit/product hits, and publication failure/
retry checks are covered by the receipts below. Zero full-S4 reads occur in both
GUI journeys; explicit compatibility reconstruction remains supported.

The declared fresh small pair observes 448.139854 -> 422.919151 s (-5.63%),
2,971,516,928 -> 3,317,284,864 B combined peak (**+11.64%**), and only ~1.7%
lower recorded read/write counters. This is not the requested substantial
complete-workflow I/O/memory improvement. The repaired process enters composite
QC plotting with ~299 MiB more RSS; the extra residency is already present after
early Q-value/evidence QC. File-backed memory is lower at the peak, while
anonymous memory is higher. The responsible R/native allocation is unproved.

This completes the authorized small verification attempt, not the performance
objective. Stop optimization here. No new 96/300 member, formal qualification,
other-omics migration or production activation is authorized. A further repair
requires a bounded scope/budget decision around the newly measured anonymous
residency; do not discard these results or repeat them for a favorable outcome.

The implemented shared helpers remove whole-input group staging/indexes and
per-batch Parquet/database round trips, reuse action-local validation without
routine recovery flushing, provide bounded Arrow writers/cursors and supported
fastduckdb projections, and reuse verified publication copies. Restoration v3
uses footer row counts; v1/v2 and public S4 compatibility remain tested. The helpers
and their current limits are documented in contract.md. All work remains on
performance-v2 and the experimental route; no dependency/backend/topology change
or default promotion is made. Historical failures and the user's resumed authority
for the row-name correction remain recorded below.

1z6Yvy is crash-interrupted; its one replacement 1viZKk was explicitly stopped by
the user during the regression investigation, after reaching precursor rollup.
Neither has complete timing/science/output acceptance. 1viZKk's observed combined
peak is 16,214,118,400 B (15.10 GiB), versus completed rRdoPJ's 9,947,648,000 B.
At the largest sampled currents, control anon/file are 8,427,675,648/1,421,168,640 B;
candidate anon/file are 5,765,808,128/9,704,394,752 B. These non-atomic samples
identify file-backed growth, not exact allocation ownership. Candidate durable
entities total 1,194,088,609 B after stopping; no temporary spill remains. The
partial memory regression remains valid adverse evidence. The previous large-run
allowances are not replenished by this small-only repair.

Initial shared-I/O verification declaration (historical): implementation and
targeted R syntax/whitespace checks were complete. Next:  existing isolated build, focused reference-core and affected
native/product/recovery checks, then a complete retained small journey only after
those pass. Small service ceilings remain 8 GiB hard / 7.5 GiB early / zero swap /
900 seconds with the existing host gates. One bounded correction pass is available
for this repair. Runtime/fixtures/runner freeze during each live check. No numerical
or performance acceptance is claimed from syntax checks or code inspection.

LukoSX installs successfully (1,268,334,592 B peak, zero swap). Focused BytFl0
fails at the incremental Arrow constructor: unnamed list input is not accepted;
the constructor errors and the obsolete v3-rejection assertion are retained in its service log (18.712468 s main interval,
1,278,291,968 B peak, zero swap). The single bounded correction pass wraps the
bounded column list as a data frame without another whole-table allocation and
packs complete scientific groups into <=10,000-row batches, avoiding a per-protein
R-kernel-call regression found during review. Exact science and I/O checks remain
unchanged. No complete small journey has been launched for this repair.

Corrected build QSTh25 installs (107.598209 s main interval, 1,141,035,008 B peak,
zero swap), RDB SHA-256
`a3b57715fb5b75bd00ac3b19a702fcf36f52fcea3e89ffb6ecfc658356bec25f`.
Core member yw7u7e is a retained valid failed check: 93.552796 s main interval,
1,660,465,152 B peak, zero swap. Its sole failed assertion is the integer-versus-
character projected row-name storage in the new shared group/map fixture; complete
native peptide S4/config/text comparisons pass, including generated imputation
rows. The earlier BytFl0 constructor failure remains retained. The correction also
updates the obsolete unsupported-v3 assertion to reject v4, packs complete groups
into bounded batches, and lets only the workflow codec own restoration metadata.
The test-service peak is not a comparable full-journey performance result.

Resumed row-name correction dbDWgV installs successfully (103.929865 s main
interval, 1,140,269,056 B peak, zero swap), RDB
`058687245763b7ad5927d7c1e825eb900531eb293e6ce0ec64eb51dd2a7fc660`.
hlNtmv passes the unchanged core assertions (86.345726 s main interval,
1,612,443,648 B peak, zero swap). One shared restoration helper preserves declared
integer row names in pages, streams and ranked heatmap projections; original
character/tibble behavior remains tested. Existing gui-native-session verification
now runs the real ten-cycle fixture and saved QC/DA product checks on this build.

aReBsh completes all ten real reference lifecycle cycles, all normalization choices,
manual RUV/Pearson Apply/Skip/reapply, native restore, IQ relocation/reapply and
QC/DA cache/scientific checks. After each cycle, working records, open resources
and payload cache are zero; ordinary reads record zero recovery flushes. The suite
retains one failed assertion: the old copy fault fixture expected an unchanged,
already verified destination to be copied again. Runtime behavior correctly skips
that copy under this repair's output-reuse contract. Preserve the failed receipt:
754.638617 s main interval, 3,423,453,184 B peak, zero swap. This is not a complete-
journey performance observation.

Keep runtime dbDWgV unchanged. Align that assertion with unchanged-copy reuse and
add the required changed-destination failure case to the existing saved-exports
file, including exact preservation of previous bytes and successful retry. Run
these checks with the already-required affected integration selector. Reuse the
passing ten-cycle/native/scientific evidence rather than repeat those cycles;
do not relabel aReBsh's overall failed receipt as passed. This fixture adaptation
verifies the explicitly approved reuse behavior and does not reset any experiment
allowance or relax scientific/output requirements.

lM6HHZ passes affected DA/publication/legacy integration, including the aligned
unchanged-copy assertion and required changed-destination failure/retry test,
on unchanged dbDWgV (120.336332 s main interval, 1,736,048,640 B peak, zero swap).
Reuse the passing native/scientific and ten-cycle evidence from aReBsh; its overall
failed receipt remains retained with the obsolete assertion clearly identified.

Release only a fresh matched SMALL pair after the existing routing/completion/
failure preflight: Qb1Btf (unfixed reference implementation) then dbDWgV (repair),
retained ykIZwA 24-sample / 192,000-row report SHA-256
`ce7868fb0873e32df845694e6b10f2ca421495d2714f87a6d49529ecf7037e51`.
Build RDB hashes are `c859941275a7ff2e1986aea235ff0be622a9fdd49d26558979da99c67d5be3de`
and `058687245763b7ad5927d7c1e825eb900531eb293e6ce0ec64eb51dd2a7fc660`.
One planned member each, serial, same small automatic-RUV/Pearson/native-first
route, 8 GiB hard / 7.5 GiB early / zero swap / 900 s and unchanged host gates.
No allocation or CPU profiler; both use the same bounded diagnostic traces and
include startup, required outputs, publication and teardown. Report this as a
small diagnostic comparison, not formal performance qualification. The fresh pair
avoids comparing elapsed time across the intervening editor crash/session break;
hpek2k/R4h60K remain valid scientific/output evidence, not discarded results.
After completion, independently compare all states at 1e-12 and required native/
DA/Summary outputs. Existing 96/300 allowlists remain unchanged; the repair is not
admitted to them. No new large member is authorized.

wuBvLx passes the existing real-browser completion/failure and routing preflight:
27.075197 s main interval, 1,351,585,792 B peak, zero swap. The repaired build is
explicitly rejected by both full-size routes. The declared small pair can proceed;
no runtime, fixture or runner changes occur while either member is live.

The declared small pair completes: lORe8g on Qb1Btf and kihu2x on dbDWgV,
both retained ykIZwA inputs, all 16 GUI states, zero full-S4 reads and zero swap.
Both GUI receipts are `passed_summary_ui_only`; independent output checks follow.

| Small diagnostic service-tree measurement | Unfixed lORe8g | Repaired kihu2x | Change |
| --- | ---: | ---: | ---: |
| Elapsed seconds, startup through tree exit | 448.139854 | 422.919151 | -5.63% |
| Combined peak bytes | 2,971,516,928 | 3,317,284,864 | **+11.64%** |
| Systemd IOReadBytes | 2,260,189,184 | 2,222,022,656 | -1.69% |
| Systemd IOWriteBytes | 2,142,322,688 | 2,106,556,416 | -1.67% |
| CPU seconds | 527.938299 | 496.846773 | -5.89% |

This fails the no-peak-regression acceptance condition and does not demonstrate a
substantial complete-journey I/O reduction at this size. Do not promote or launch
96/300. Finish the already-planned independent science/output verification, then
record the bounded outcome without adding optimization or favorable reruns.
These are one diagnostic pair's observed counters, not statistical qualification.
Device/OS cache effects and filesystem write amplification are not isolated.

84Kict independently passes all 16 states against the fresh control at unchanged
1e-12, exact identifiers/order/types/attributes/provenance and unchanged source
heads (57.127032 s main interval, 1,710,559,232 B peak, zero swap).

Both sampled peaks occur during `savePlot` for `proteomics_composite_QC_figure`.
At maximum sampled current, control anon/file/kernel are
2,602,962,944 / 267,874,304 / 25,006,080 B, with 75,599,872 B residual;
repair anon/file/kernel are 3,009,347,584 / 189,161,472 / 25,321,472 B,
with 92,909,568 B residual. These components are non-atomic observations and may
overlap; the residual is not a measured allocation category. Independent maxima
file bytes are 397,742,080 vs 384,425,984, and kernel 32,145,408 vs 30,318,592.
The remaining regression is anonymous memory, not a repeat of the 96-sample
file-cache-dominated attribution.

Main-process RSS already differs before plotting: before Q-value 1344.7 vs
1213.2 MiB; before protein-evidence 1282.8 vs 1423.9; before rollup 1274.1 vs
1577.2; before normalization 1280.3 vs 1748.6; entering the composite 1323.0 vs
1621.6. Plot assembly then adds a similar ~658/~653 MB before savePlot. This
locates persistent extra residency after early QC but does not distinguish live
R owners, unreachable R allocations, Arrow/DuckDB buffers or allocator retention.
No causal allocation claim is inferred from RSS alone.

Repair delivery counters: zero payload-cache bytes/open resources; 200 metadata
entries / 1,263,200 B; 329,543,554 payload-hash bytes; zero routine recovery flushes;
13,408,325 staged-stream bytes; 4,157,960 maximum sampled scratch bytes;
141,472,084 accepted store-output bytes; six query streams, 106 sequential streams
and four fastduckdb reads. Logical returned-table accounting is 122,492,760 B
versus control 113,113,240 B; decoded bytes are not physical I/O or retained heap.
These show implemented mechanisms, not a complete-workflow memory win.

65eHxE retains a verifier-invocation rejection before any verification work:
the supplied dbDWgV library did not match control lORe8g's frozen Qb1Btf library.
The runner correctly refused it (1.558725 s main interval, 66,215,936 B peak,
zero swap). Correct only that command binding and verify the control with Qb1Btf;
verify the candidate with dbDWgV. No benchmark member or runtime is replaced, and
no scientific/output comparison or resource limit is relaxed.

Final independent required-output verifiers pass with matching producer libraries:
79a3Zr for control lORe8g (65.734181 s main interval, 1,709,084,672 B peak,
zero swap) and Jpj6oJ for repaired kihu2x (52.359737 s, 2,157,694,976 B peak,
zero swap). Both report `passed_ruv_pearson_output_checks`, complete native archive
integrity and true RUV/Pearson trace coverage. Their verification costs are outside
the compared GUI journey intervals. No proof service remains live after completion.

Counter limitation: `table_rows`/`table_bytes` cover the existing page/stream
accounting and do not include the new group-buffer assembly path. They are not
total decoded workload or physical I/O and must not be used as a before/after
materialization-volume claim. The cgroup service counters and phase RSS evidence
used for the unfavorable performance decision are independent of this limitation.
Required archive/TSV/figure output, filesystem amplification, allocator retention,
process overlap and residual/unattributed memory remain measured or stated costs;
no zero-I/O or constant-total-memory claim follows from bounded R reads.

Implementation pass now covers the included C–E path, prior to any new test run:
- Existing state/table codec: column-blocked assembly, bounded batch production,
  component commits, explicit dependencies, GUI role roots and reference-only Revert.
- Peptide QC: complete independence groups, full-design denominators, original
  population intensity quantile, direct persisted SQL precursor rollup, native
  scientific kernels, durable support tables and bounded progress/output products.
- Genuine IQ: streamed aliases/compatibility values, request-owned scientific
  files and separate publication; protein kernels consume scoped components and
  deliver only references, summaries and display products.
- Normalization/RUV/Pearson: reusable science and render identities, original
  Pearson population, metadata-only Skip, file-backed figures and durable output
  journals; retries reuse committed science.
- DA: metadata initialization, reusable native base/contrast fits, reference
  result manifests, projected table/volcano reads, linked Glimma expression and
  selected heatmap matrices with separate clustering/rendering.
- Native archive v2 with complete declared roots and v1 readers/conversion,
  portable output acknowledgements, explicit guarded S4 exports, reference-backed
  Summary and passive enrichment. Active enrichment/report migration stays held.

At completion of the implementation pass it was **unverified**. No test/build/benchmark was run while these
adapters were being implemented. The following bounded item was syntax/build and focused
scientific, ownership, publication and old/new-session checks using the existing
fixtures/runner, then the complete small journey. No new large-run allowance or
performance claim is implied; all historical failures remain retained below.

Full implementation build xC9twh installs successfully. Focused 9TKJKR fails
(38.462954 s main interval, 1,696,313,344 B peak, zero swap): the native GUI's
`rollupFn` provenance alias was omitted; the imputation summary attempted to
construct an invalid empty S4; and equivalent table attribute layouts failed the
component-content cache key. Peptide numerical/table/config/text comparisons
through replicate filtering otherwise pass. Retain this failed member.
The correction preserves the original alias, reads imputation summary metadata
directly and fingerprints canonical codec metadata plus individual columns. The
latter also avoids serializing an entire assay solely to identify a component;
no memory improvement is claimed without measurements. Cache consumer routing
now uses the existing exact typed factors, and the legacy getState() seam keeps
its original no-argument call. A corrected build and affected checks follow.

Corrected build Lpgjvx passes bXU4m0 core checks: 38.294544 s main interval,
1,562,898,432 B peak, zero swap. Full Q-value/evidence/rollup/intensity/sample/
replicate/imputation S4 science, configuration and text match the native oracle,
including absent-observation imputation ordering. Wide reconstruction and
metadata-only component reuse pass. The strict oracle was retained. The next
existing `gui-native-session` check also includes the existing QC-product and
DA-fit files; its ten-cycle fixture now calls the complete reference GUI paths
and forbids full-state reads during the actions. No capacity member is released.

8VuP7T retains the first D/E check failure (47.159383 s main interval,
1,973,440,512 B peak, zero swap). QC science-product and DA base/contrast cache
checks pass. The lifecycle renderer supplied NULL shape metadata to the existing
strict render-key validator, and IQ tried to sync an export-stage file through
the store-only sync boundary. Omit absent render-only shape fields and use the
existing export-path sync helper. The old IQ-relocation assertion is updated to
expect a reference handoff and verify its explicit public S4 reconstruction.
Numerical assertions remain unchanged; the lifecycle/recovery checks are rerun.

qYBl0Y/JVrYMM completes all ten reference normalization/Revert/style/export cycles,
including all four offered normalization methods, manual RUV, Pearson Apply/Skip/
reapply, full native archive import and genuine IQ reapply after relocation.
QC product hits and DA base/contrast fit hits pass their no-hydration/no-refit spies.
The suite remains failed only because the new DA GUI fixture passed a matrix's
rownames alongside an ID column to the unchanged native oracle. Remove only that
fixture row-name duplication, then run the affected DA file with the required
legacy/module regressions on the unchanged build. Retain JVrYMM: 352.251466 s main
interval, 3,118,002,176 B peak, zero swap. Matching passing lifecycle evidence is
reused; the ten cycles are not repeated for this fixture-only correction.

Ls2KKv retains the affected-integration failure (63.588215 s main interval,
1,494,646,784 B peak, zero swap). Two legacy correlation-observer injected seams
require checking the workflow container before reading its manager. The new DA
publication path reached native results, then its count-summary product retained
a plot environment; persist only the count table used by the output consumer.
IQ lineage lookup is kept transient, preserving existing argument provenance.
The stricter IQ fixture now also compares complete scientific arguments.

uzfJL8 passes the legacy/module regressions and injected action-publication
process-loss recovery. Its new DA fixture then exposes the existing writer's
non-NULL annotation requirement (63.565613 s, 1,500,803,072 B peak, zero swap).
The reference consumer now supplies an empty annotation table when annotations
are absent, selecting the writer's existing no-annotation branch. Preserve the
native count-table type while removing only unused plot owners. This failure is
retained; science and output assertions are not relaxed.

wVpW3z/osnk6B passes affected integration checks (65.106784 s main interval,
1,509,249,024 B peak, zero swap): exact native DA tables, complete DA publication,
projected pages/counts, selected-matrix heatmap equality, no reclustering or assay
reads on colour/label hits, no repeated DA science, legacy/module APIs and durable
output-only recovery. Core science and the matching ten-cycle evidence are reused.
Freeze wVpW3z RDB `95d28c206dd50163f7e7fc755ca06219e610136e91abdb4358c8089fca869747`
for one complete small GUI proof on retained ykIZwA (24 samples/192,000 rows).
Purpose: all included QC/IQ/normalization/automatic RUV/Pearson/DA/Summary actions,
automatic native-first outputs and teardown, followed by independent comparison
against R4h60K. Profiling off, existing 8 GiB hard/7.5 GiB early/zero swap/900 s
small guard. This binding is admitted for the small native route only; the stress
allowlist is unchanged. No new large member or production promotion is authorized
by these focused checks.

Complete small GUI dXa6av finishes on wVpW3z/retained ykIZwA: 238.044472 s
main interval, 2,987,560,960 B peak, zero swap, all included QC/IQ/normalization/
automatic RUV/Pearson/DA/Summary actions and teardown. Status is
`passed_summary_ui_only`: independent science and output checks are still required.
There are 16 states (15 component scientific states plus initial), zero full S4
reads, and no resident assays in the inventoried workflow holders. The all-state
R4h60K comparison and required-file/native archive verifier are next. No scaling
or performance conclusion follows from GUI completion alone.

BpsG34 retains the first full-history comparison failure (18.303017 s,
1,141,551,104 B peak, zero swap). Peptide states match; protein_s4_created differs
only in the two declared IQ file locations. Extend the existing IQ-replay path
normalizer to the component full-journey verifier, requiring matching SHA-256 for
oracle file, immutable scientific file and published candidate file, plus exact
owned-directory/basename checks. All other record fields and 1e-12 tolerance
remain unchanged. This is a file-layout comparison adaptation, not a science
waiver or a replacement GUI member.

pmz11T and bounded all-state diagnostic HBKmAR preserve the remaining strict
comparison failures. Both IQ files are byte-identical, including their published
copies. No numerical/table population differences are reported. The differences
are IQ/protein table attribute ordering and the original GUI's dynamic RUV/Pearson
parameter aliases (`ruvIII_C_VaryingFn`, cleanup and filter aliases). Record bounded
attribute-name order in the existing v2 codec (old v2 metadata remains readable),
restore it exactly, and invoke the unchanged kernels through the original aliases.
The comparator stays strict. HBKmAR: 25.061783 s, 1,112,518,656 B peak, zero swap;
pmz11T: 18.266444 s, 1,147,609,088 B peak, zero swap. A fresh small member will be
required after these compatibility corrections; dXa6av is not relabeled as passed
science/output qualification.

AEn1eg independently validates dXa6av's complete native archive, 16 states,
766 proteins/24 samples, required IQ/QC/DA/Summary files and DA calculations,
then fails the old wrapper-based trace requirement (26.097417 s,
1,705,414,656 B peak, zero swap). Keep this failed receipt. Trace the new reference
normalization/correlation entry points while retaining all genuine RUV/Pearson
numerical-kernel labels; legacy routes retain their old wrapper requirements.
No trace is invented for dXa6av. The already-required fresh small GUI after the
attribute/provenance fixes must provide complete new-route instrumentation.

xHO6KQ retains the attribute-restoration check failure (38.391977 s,
1,563,570,176 B peak, zero swap): full peptide parity passes, but reassigning the
ordered attribute list changed R's compact automatic row names into explicit row
names and invalidated component reuse. A minimal R check confirms `NA,-3` became
`NA,3`. Preserve `.row_names_info(..., 0L)` while restoring attribute order; keep
the metadata-only no-rewrite assertion unchanged.

3DXRdV/bPPRKk passes the strict codec, wide-table, metadata-only reuse and complete
peptide science/provenance checks: 38.626941 s, 1,563,774,976 B peak, zero swap.
Freeze RDB `2121fbae6a9a9a61a5fa5576c1359abdced672b8343002adb589f97bbe7babca` for the corrected complete small GUI on unchanged ykIZwA,
with the same 8 GiB/7.5 GiB/zero-swap/900 s guard and new-route bounded traces.
The failed prior members and comparisons remain retained; stress admission remains
held pending this member's independent all-state/output checks.

Corrected full GUI 6wetoP completes on 3DXRdV: 239.363177 s main interval,
3,024,310,272 B peak, zero swap, zero full-S4 reads, all included stages and
outputs delivered. Strict 5xPHgI matches all states through RUV and reports only
the insertion order of differentialAbundanceAnalysis versus daAnalysisUI in the
final state (27.045344 s, 1,147,863,040 B peak, zero swap). Preserve the original
GUI ordering, add the absent-DA-section shape to its existing fixture, and retain
both receipts. Numerical and attribute comparisons remain strict.

MXf1VI passes complete native/required-output/DA/Summary checks and all new-route
RUV/Pearson trace requirements for 6wetoP (25.759621 s, 1,702,129,664 B peak,
zero swap). Qb1Btf/f0NTee passes the DA provenance-order correction and affected
integration suite (65.460767 s, 1,502,941,184 B peak, zero swap).
Freeze Qb1Btf RDB `c859941275a7ff2e1986aea235ff0be622a9fdd49d26558979da99c67d5be3de` for the final complete small proof on unchanged ykIZwA.
Only the existing small native route is released, under unchanged guards;
independent all-state/output checks still gate the full-size pair.

Final small hpek2k on Qb1Btf/ykIZwA completes the whole included GUI journey:
243.524684 s main interval, 2,950,512,640 B combined peak, zero swap, 16 states,
zero full-S4 reads and no resident assays in inventoried GUI holders. mzisVb
passes all 16 states against R4h60K at unchanged 1e-12 with exact ordering,
types, attributes and provenance (28.519745 s, 1,159,610,368 B). IQ path relocation
is normalized only after three-way file hash and ownership checks. 0IvoLU passes
native archive integrity, required IQ/QC/DA/Summary files and genuine automatic
RUV/Pearson trace coverage (27.718314 s, 1,700,556,800 B, zero swap).

The planned full 96 pair binds jfR7Ma RDB
`f529ddc96a9bc22883b2f984ba97edb156bd32284cf1add528cbbadfe0725399`
and Qb1Btf RDB
`c859941275a7ff2e1986aea235ff0be622a9fdd49d26558979da99c67d5be3de`
to nDqUl6 report SHA-256
`3d670aa2567263e9ae90e2966b570e4d65200e9d5adac7214b11eb705a01032a`.
One planned member per build, serial, profiling off; startup, all automatic
native-first outputs, publication and teardown included. Limits remain 16 GiB
hard, 15.5 GiB early stop, zero swap, 3600 s, 16 GiB startup headroom and 2 GiB
continuous host-memory stop. The existing frozen historical bindings stay intact;
neither new build is admitted to 300 until the 96 decision. No failed 300 control
repeat is authorized.

8Wwf9J passes the existing small completion/failure, real browser delivery and
exact build/input/guard routing checks: 16.740435 s main interval,
1,359,560,704 B peak, zero swap. The native-first control rRdoPJ starts on the
declared jfR7Ma/nDqUl6 binding. Runtime, tests and runner remain frozen while live.

Hydration-audit closure for the included route: hpek2k records reference-backed
import/design/DA workflow holders, zero full-S4 reads, zero open resources and a
bounded 200-entry metadata cache after delivery. Core/DA fixtures additionally
inspect normalization/correlation and result-manifest holders. JVrYMM's ten real
normalize/Revert/style/export cycles forbid GUI materialization and check owners,
resource/cache ceilings and complete native restores after each cycle. Product
and fit-hit spies forbid assay access or repeated fitting. Native/publication
fixtures cover missing/corrupt components, replacement and durable process-loss
output-only recovery. mzisVb verifies the explicit public compatibility boundary
against all original scientific records. These establish the tested ownership
boundaries, not an attribution of anonymous memory or a capacity claim.

Planned full-96 control rRdoPJ completes on jfR7Ma/nDqUl6:
758.557128 s service-tree interval (758.293066 s main interval),
9,947,648,000 B combined peak, zero swap, all 16 states and included GUI outputs/
teardown. GUI status is `passed_summary_ui_only`; independent checks remain.
Bounded samples: maximum anonymous 8,427,675,648 B, file 2,394,624,000 B,
kernel 43,413,504 B, dirty 73,486,336 B, writeback 91,779,072 B, each at its own
time. At maximum sampled current 9,881,645,056 B: anon 8,427,675,648 B + file
1,421,168,640 B + kernel 39,514,112 B; residual -6,713,344 B reflects non-atomic
observations/overlap and must not be treated as an allocation decomposition.
I/O: 34,244,825,088 B read, 17,701,863,424 B write. Candidate 1z6Yvy begins on
the declared Qb1Btf/nDqUl6 binding, same route/settings/outputs/guards, profiling off.

1z6Yvy is interrupted during precursor rollup when the user reports VS Code
crashed. On resumption, its service is inactive/dead (MainPID 0), controller and
terminal receipts are absent, and the GUI receipt remains `running` with four
states through protein-evidence filtering, zero full-S4 reads. The journal records
an explicit stop at 01:18:56; the controller's existing EXIT trap stops its service
on interruption. Preserve [the journal](1z6Yvy-interrupted-service.log) and all
partial member files. No OOM/early-memory failure is recorded. This is an
incomplete interrupted member, not a completed timing or scientific comparison.
Its observed peak 13,072,822,272 B already exceeds the control's 9,947,648,000 B;
the interruption does not erase that unfavorable measurement.

Consume the single documented instrumentation-invalid replacement for this
interrupted candidate, with exactly Qb1Btf/nDqUl6, unchanged native-first outputs,
settings, profiling-off route and 16 GiB/15.5 GiB/zero-swap/3600 s/16 GiB startup/
2 GiB host-stop limits. Reuse rRdoPJ; do not rerun the valid control. Matching
8Wwf9J completion/failure instrumentation evidence remains valid (no runner/runtime
change). The purpose is to finish the scientific/output and complete elapsed/peak
result. The 300 candidate stays held while the measured 96 memory regression
remains; the replacement is not permission for another optimization pass.

### Retained prior implementation and experiment evidence

**Resumed by explicit user instruction: "please continue all the way through".**
This supersedes the exhausted focused-correction stop, authorizing necessary
in-scope implementation/verification through the remaining slices. It does not
raise resource limits, change scientific acceptance, erase failures, authorize a
failed 300 control repeat or activate production defaults. The original large
experiment sequence and allowances remain the experiment boundary.
Current bounded item remains B (DIA-018/019/025). XCgkv0/kixnQf repairs
empty tokenization and legacy metadata lookup. Complete reconstructed Q-value S4
science/provenance, text and counts pass on the edge fixture; one page expectation
incorrectly drops the custom table attribute preserved by the established reader.
Retain that failed check (59.997798 s main interval, 1,272,442,880 B peak, zero swap).
Correct its oracle, finish unchanged-assay reuse and durable progress/output-only
recovery using the existing adapter/codec/journal, then the small browser proof.
PAx5Vo/GoGFsV passes first-QC scientific parity, table reuse and failed-publication
reopen/retry checks: 74.983578 s main interval, 1,400,377,344 B peak, zero swap.
The retry fixture forbids assay reads, new scientific commits and repeated renders
once figures are staged; prior output bytes remain retained. Stream/schema readers
and the first-QC branch of the existing browser driver are being completed before
its small pilot. sAExGg/vq9YGl passes the streaming (including exact TSV bytes and interrupted
write), science, recovery, ownership and old/new guard-routing assertions; its
added Chrome input-delivery fixture failed because the selector supplied an empty
browser path. Retain this failed check (363.097442 s main interval, 1,490,329,600 B
peak, zero swap). Repair only the launcher selector and rerun the affected
browser/routing file via the existing reader-check dispatcher. Runtime is unchanged;
matching passing checks are reused. Freeze first-QC-only candidate sAExGg RDB
SHA-256 `07c3792c18bd04ee221e1c488d0c84752ee38ba3de2f56f89ccbed89e03bc375`; it is not approved for the full journey.
FEiQSU passes the repaired browser/routing checks on unchanged sAExGg
(37.437862 s main interval, 1,268,285,440 B peak, zero swap). Small first-QC pilot
z9XUPG stops at Q-value output production (62.949146 s main interval,
1,913,151,488 B peak, zero swap). Import/design retain references and record zero
full S4 reads; the staged Q-value table lacks the physical element-name column
required by its Run restoration metadata. Its bounded footer proves this mismatch.
Carry that existing named Run column through selection/production and through the
small design-cleaning input; add the actual named-vector shape to the parity fixture.
This is a correctness failure, not a successful first-QC or a large-run replacement.
The named-column repair vzoyaI/aiTgEM passes the affected state/codec/Q-value and
publication-recovery checks: 42.451951 s main interval, 1,425,166,336 B peak,
zero swap. Freeze vzoyaI RDB `5aec3d13c76fcb7993a392a55a5f7f9ba2e3fdea3a5ac834242437020ffc3f3a` for the corrected small first-QC pilot;
replace the failed sAExGg first-QC launch binding, retaining its failure receipts.
The existing independent history verifier now has an explicit first-QC subset:
initial/raw/Q-value records plus six PNGs, compared with unchanged tolerance and
attribute checks against retained R4h60K (small) or gW09Ul (96). This does not
qualify later states or turn the retained old-output pair into a new comparison.
Corrected small GUI 5XmnS4 completes first QC, all 12 required figures and Revert:
68.690912 s main interval, 2,861,514,752 B peak, zero swap, zero full S4 reads and
no resident assays in the inventoried holders. Independent EKG269 (13.247256 s,
1,148,698,624 B peak, zero swap) passes raw-state comparison, then finds the
candidate decoder drops the empty Run `names` attribute in the zero-row exclusion
ledger. The stored metadata and legacy checkpoint retain it; a small qs2 check
confirms serialization preserves it. Restore this attribute only in the zero-row
materialization path and add an empty named-table regression. Keep the strict
comparison and the failed receipt. aFoxZd/TmewTD passes the repaired zero-row named-table and affected core checks:
68.585385 s main interval, 1,417,134,080 B peak, zero swap. Freeze aFoxZd RDB
`ce4cd922073370236f61e7efee8437a18411beca74895a7e1344259a6411843b` for the final small pilot and comparison; the repair changes
compatibility decoding, so the pilot and verifier will use the same new build.
ZqHJWb completes the repaired small GUI (96.099717 s main interval,
2,899,329,024 B peak, zero swap). Independent 8jZZQi passes initial/raw/Q-value
records with unchanged 1e-12 tolerance, exact attributes/types/missingness and all
six identical PNGs; both store heads remain unchanged (27.350838 s main interval,
1,163,345,920 B peak, zero swap).
Next authorized large member: one 96-sample first-QC candidate on aFoxZd and
retained nDqUl6 (7,680,000 rows, report SHA-256
`3d670aa2567263e9ae90e2966b570e4d65200e9d5adac7214b11eb705a01032a`).
Purpose: prove Import/Design/Q-value/automatic figures/Revert with references at
96 samples. Profiling off; 16 GiB hard, 15.5 GiB early, zero swap, 3600 s,
16 GiB startup headroom and 2 GiB continuous host stop. One planned candidate;
independent first-QC science/PNG comparison against retained gW09Ul follows.
No failed 300 control is repeated and this member is not a speedup comparison.
The first admission refused below 16 GiB, with no service/member started. After the
user requested retry, host available was 17.23 GiB and no proof job was live.
96 first-QC member **5TcekA failed in Design** on aFoxZd/nDqUl6.
Import persisted 7,680,000 rows. The Design COPY failed when DuckDB attempted
`/tmp/RtmplaBR8b/duckdb/temp`, an expired temporary root, during disk spill.
Observed peak was 8,120,528,896 B; no memory limit was reached. The application
ended with the initial head, zero component states and zero full S4 reads.
The timed browser wait did not notice session loss; its idle service was stopped
explicitly after preserving application logs. The outer guard then recorded
`missing_memory_observation` because the service had been removed. That secondary
receipt is not the scientific failure cause and is not an instrumentation-only
replacement. Retained stderr/stdout accompany this record as `5TcekA-shiny-*.log`.
Bounded correction: the existing table scope owns a fresh spill directory under
its store and deletes it after disconnect; the ordinary browser completion wait
rejects a disconnected Shiny socket. Focused checks must force real disk spilling
and verify cleanup on completion/failure, plus a real browser completion/disconnect
lifecycle and the small first-QC proof. Only then use the one large correction
pass for an affected 96 first-QC repeat (at most three additional affected members;
none started yet). The instrumentation-invalid replacement remains unused.
Correction build 6VEW4I RDB SHA-256
`0eb8cefabb09f6848f642407f098c3e3b56a0a87c75d65d95ee3a34909c76de4`
is frozen. OQxOqY passes reference-state science, publication recovery, ownership
and the real timed browser completion/disconnect lifecycle, but fails the existing
NULL-reader-selector harness assertion and the new spill fixture's COPY at its
artificial 32 MB setting (376.060094 s main interval, 1,799,442,432 B peak, zero swap).
Disk spilling itself succeeds and cleanup runs; the Parquet writer needs more
than that fixture-only budget. Retain this failure. Make the dispatcher predicate
scalar and restore the actual unchanged 1 GB DIA query budget after forcing spill
in the fixture; rerun the affected focused suite and current build binding.
Neither correction raises an application or job limit.
d3xIve passes the affected focused suite and frozen launch bindings on unchanged
6VEW4I: 341.286587 s main interval, 1,812,348,928 B peak, zero swap. Real disk
spilling, scoped cleanup after writer failure, exact state/codec parity, publication
recovery, browser completion/disconnect and legacy dispatcher checks pass.
Next: one small first-QC GUI on retained ykIZwA with the existing 8 GiB/7.5 GiB
limits and 900 s, followed by the strict R4h60K first-QC comparison. A successful
small proof releases the affected 96 repeat, with all original full-size guards.
LOWoIb passes the complete small first-QC GUI/outputs/Revert (91.812061 s main
interval, 3,074,260,992 B peak, zero swap). It records zero full S4 reads, no
resident assays in all three inventoried holders and zero open store resources.
iQyTjB independently passes initial/raw/Q-value equality and all six exact PNGs
against R4h60K, preserving original heads (26.652319 s, 1,078,272,000 B peak,
zero swap). The next large member is the corrected 96 first-QC repeat on 6VEW4I
and unchanged nDqUl6. This uses the one large correction pass and its first of at
most three additional affected members. Profiling stays off, with 16 GiB hard,
15.5 GiB early, zero swap, 3600 s, 16 GiB startup and 2 GiB continuous host stop.
Corrected 96 member **VxuMT5 completed** on this binding: 340.708626 s main
interval, 15,187,730,432 B combined peak (14.14 GiB), zero swap. First QC, all
12 required figures and Revert pass; zero full S4 reads, zero open resources,
no resident assay in the three inventoried holders and no remaining scoped spill
directory. EsHLi4 independently passes the gW09Ul first-QC science/output comparison
(93.041721 s main interval, 6,902,845,440 B peak, zero swap), with unchanged
1e-12 tolerance, exact attributes/identifiers/missingness, six identical PNGs and
unmodified source heads. This releases the conditional 300 first-QC candidate.
Measured cost remains material: sampled maxima are 8,479,125,504 B anonymous,
6,733,123,584 B file cache, 58,351,616 B kernel, 123,768,832 B dirty and
58,699,776 B writeback. These independent maxima must not be added. At the
largest observed current (15,187,292,160 B), anonymous/file/kernel fields were
7,934,877,696 / 6,649,524,224 / 58,208,256 B; the 544,681,984 B residual is
unattributed by these non-atomic observations. Service I/O records
19,686,547,456 B written. No allocation profiler was enabled and transient phase
attribution is incomplete. This is a first-QC capacity proof, not a performance
comparison; native-first and reference-only delivery do not establish lower
active computational peaks.
Next planned capacity member: 6VEW4I with retained OQADra (300 samples,
24,000,000 rows, 3,475,378,370-byte report). Purpose: Import/Design/Q-value,
durable outputs and reference Revert at 300 samples. One planned first-QC
candidate; unchanged 16 GiB hard/15.5 GiB early/zero swap/3600 s/16 GiB startup/
2 GiB continuous host limits. Profiling disabled. The failed 300 control is not
repeated. The large correction pass has already been used for VxuMT5 (one extra
affected member); instrumentation-invalid replacement remains unused. A valid
limiting failure is retained and blocks later adapters under slice B's exit rule.
The full 96/300 sequence remains unstarted.
Required files add the existing saved-export helpers and first-QC branch in the
existing runner/driver; no new harness, dependency or service. Exit remains the
small science/output/ownership checks followed by 96/conditional 300 first QC.

`UCTKOQ` is stopped (service inactive/dead, MainPID 0). Its retained guard receipt
reports `early_memory_stop`, current 16,790,175,744 B, peak 16,790,544,384 B,
zero swap and host available 7,956,357,120 B. The last committed scientific state
is `raw_data_s4`; the final owner event enters first QC. Q-value completion is
not evidenced. Keep the verified 4BVSnf/gW09Ul 96 pair and all failed runs intact.
The old pair used native archives plus standalone S4 outputs; it is not the
native-first comparator for this continuation.

Completed bounded item A (DIA-025/022): correlation output policy in
`mod_prot_norm_correlation_helpers.R`; archive inventory/copy helpers in
`utils_workflow_dia_sessions.R`; Summary integration output in
`mod_prot_summary_session_helpers.R`; publication inventory in
`func_general_filemgmt_results_summary_helpers.R`; corresponding native-session,
Summary and verifier tests; existing `run_workflow_checkpoint_proof.R` verifier.
No scientific/state execution changes in the output-policy-only control. Freeze
that control from x3Dkg2 before implementing B. Exit: focused native/legacy output
checks pass and output-only source/build identity is retained.

Native-first required inventory: cleaned-data/IQ/QC/normalized/final TSVs,
existing QC figures/composites and DA TSV/workbooks/volcano/heatmap artifacts,
metadata RDS files that contain no assay S4, one complete native session archive
and its completion/latest receipts. Summary copies the complete native archive
into Integration and Publication_tables and retains study/session metadata,
figures and TSVs. Automatic `*_cln_with_replicates.RDS`, Summary `*_final_s4.RDS`
and normalized publication S4 RDS are compatibility exports only. Descriptors
must never be written under those S4 filenames. Legacy report/workbook consumers
that require these files need explicit compatibility export; their migration is
held. Publication failure retains recoverable inputs and previous completed data.

Hydration audit / closure criteria (not yet satisfied): adapter `getState`, public
Revert and history currently decode whole records; import/design still build S4;
Q-value/remaining QC operate on full S4; rollup collects saved output; norm/RUV/
Pearson access source S4 before products; DA initializes/materializes and retains
collections; passive enrichment and Summary resolve S4; native export prerequisite
can hydrate. Replace GUI boundaries with state/component references, bounded
metadata/projected reads, one required scientific matrix and persisted products.
Public full reconstruction remains explicit. Verify each inventoried holder has
no retained assay/S4 after delivery, component closure packages without assay
hydration, valid product/base-fit hits never refit, and output retries create no
new scientific checkpoint. Ten existing lifecycle cycles must pass. No silent
new-format fallback, altered scientific populations or raised registry/read caps.

Budget: sequential A–E, one correction per slice; 96 first-QC, conditional 300
first-QC, matched native-first 96 pair, conditional complete 300. Across large
experiments one correction pass capped at three affected members and one
instrumentation-invalid replacement. No member started for this continuation.
Full-size guards unchanged: 16 GiB hard / 15.5 GiB early / zero swap / 3600 s /
16 GiB startup / 2 GiB continuous host stop. Formal gates/default activation held.

Slice A verification: DGp7xx is the output-only build from frozen x3Dkg2;
comparison of all runtime R files before installation finds only correlation,
Summary, native-session output/copy helpers and the declared output-policy field.
Native suite zCuKiG passes under the existing guard: 111.433056 s main-process interval,
2,159,968,256 B peak, zero swap. All ten cycles complete with zero working/native
records, cache bytes and open resources; this is the existing lifecycle boundary,
not yet the stronger reference-backed holder proof. Native restore, process-loss
recovery, real IQ reapply, omitted automatic S4 and verified portable Summary copy
checks pass. Source-only contract review then preserves the old non-native
notification text before control rebuild N7eHpR and the affected output regressions.
N7eHpR/HevpZs affected output, Summary, normalization, verifier and runner
regressions pass (181.891525 s, 1,490,014,208 B peak, zero swap; three deliberate
copy-failure warnings retained). The one slice-A correction pass also closes the
legacy-import standalone S4 write, with an explicit compatibility-export argument;
its existing restore fixture now checks no automatic S4 and exact explicit export.
Publication inventory and native-first verifier branches receive focused checks.
Final corrected control jfR7Ma is installed; RDB SHA-256
`f529ddc96a9bc22883b2f984ba97edb156bd32284cf1add528cbbadfe0725399`.
Its source differs from x3Dkg2 only in the four named output R files and DESCRIPTION.
Correction check FqBod5 passes under 8 GiB/900 s: 265.214770 s,
2,229,837,824 B peak, zero swap. Native restore/recovery/ten cycles, explicit S4
compatibility, Summary publication and verifier policy checks pass. Slice A is
complete; jfR7Ma is the frozen native-first output-only control.

Current bounded item B (DIA-018/019/025): extend
`R/utils_workflow_saved_table_codec.R` and `R/utils_workflow_dia_state.R` for
streamed table production and component checkpoints; connect existing import
processing, design builder/state helpers, peptide Q-value/preview/GUI Revert,
filtering-progress helpers and native dependency closure. Extend the existing
projection/state/Q-value tests and guarded browser runner for first-QC completion.
Exit: exact small-fixture values/types/attributes/order/provenance, reference-only
metadata/Revert/preview and released owners, then one 96 first-QC and conditional
300 first-QC under unchanged limits. Remaining adapters must wait for this gate.
Initial B build tH7W6B/check 243EJ2 remains failed: 30.080535 s main-process
interval, 1,271,394,304 B peak, zero swap. Findings were a missing global-config
fixture binding, an incorrect page-reader row-name expectation, and unsupported
readr col_spec fields. The one correction pass repaired those and reviewed ledger
argument order, integer row names, missing identifiers and absent design rows.
M0S2X8/eTxeXq are retained intermediate build snapshots, with no recheck on either;
the interruption/resumption did not reset the allowance.

**Corrected nFxCX2/EdTu02 fails:** 25.686161 s main-process interval (not a complete
workflow timing), 1,268,011,008 B combined service peak, zero swap. One error and
three failed assertions remain; five missing-parameter fixture warnings are retained.
The new query-producer and readr parser assertions now pass, but:

- Reference Q-value execution calls the existing R method on a zero-row schema
  input. `.tokenizeProteinAccessions()` converts that empty input to logical(0)
  through `ifelse()` and then fails in `strsplit()` with "non-character argument".
  The real reference-backed scientific result and its parity checks were not reached.
- Adding `row_names_type` collides with legacy `$row_names` partial matching when
  the latter field is absent. Three existing version-2 metadata assertions now fail.
  This is a compatibility regression; do not change those assertions to accept it.

Exact receipts: `browser-smoke.EdTu02/{reader-tests.rds,reader-operation.json,
service.log,terminal-properties.txt,terminal-interval.json}` under the retained
`MultiScholaR-dia-proof-20260911-VhC1n5` root. Candidate RDB SHA-256:
`79bbb634990aeb400325c0ef8a8f1632cafc7b2ff1622dbefc9a4be4adb111a4`.
Driver SHA-256: `3a770060f27f19c9119ac7dfacc46968194dee06bc6ecfbbe0aed3e0622d99f2`.
The service terminated with exit 1; no live proof job remains.

Historical decision before the renewed user instruction: stop under the
one-correction limit. A next authorized correction would cover the empty-input tokenizer contract,
a scalar metadata name that does not change legacy partial lookup, and the same
focused recheck under 8 GiB/900 s/zero swap. Keep both failures and every assertion.
That decision would not reset the large-run allowance or qualify B by itself.

B still requires its complete browser Import/Design/Q-value/output/Revert/selected
read proof, explicit unchanged-assay reuse for metadata-only changes, durable progress
products and output-only recovery, remaining failure/ownership cases, then the 96 and
conditional 300 first-QC members. The full remaining QC/IQ, normalization/RUV/Pearson,
DA/native-format/downstream migration, matched 96 pair and final 300 journey remain
unimplemented or unverified. No large experimental member, large correction member
or instrumentation-invalid replacement has been consumed. Production activation,
formal qualification, other omics and active enrichment/report migration stay held.

## Historical decision — approved active-owner repairs, revision 1.42

User explicitly approves **16 GiB startup headroom** for the final 96/300 queue.
The **16 GiB job cap, 15.5 GiB early stop, zero job swap, 2 GiB live host stop and
3600-second ceiling remain unchanged**. No extra timed attempts or application
changes. Check the altered startup predicate, then launch the frozen control.

**Final timed control 4BVSnf completed**: installed xz9cb6, retained nDqUl6
(96 samples), profiling disabled, native session and complete GUI/DA/Summary
boundary. Hsv9xZ passed the startup-only exception check. Exact build/driver/spec
hashes and unchanged limits are frozen in DIA-025. Service lifetime 1647.638822 s,
combined peak 9,379,393,536 B (8.735 GiB), zero swap; all 16 states and included
DA/Summary actions completed. Independent native/output verifier **Gn9xL6 passes**
(8,950,530,048 B peak/zero swap, 203.939460 s; separate from timed GUI cost).
**Candidate gW09Ul completed** on frozen x3Dkg2 and the same nDqUl6 inputs,
settings, native outputs and 16 GiB timed boundary: **934.165324 s, 8,951,353,344 B
peak (8.337 GiB), zero swap**. All 16 states and included DA/Summary/teardown
complete. Against fresh 4BVSnf this is **43.30% lower elapsed and 4.56% lower
combined peak**: no regression, but not a 10% memory reduction against this
control or broader adoption approval. Independent output verifier **1aeDgR
passes** (8,467,513,344 B peak/zero swap, 120.923285 s). Strict all-16-state
verifier **DrkV5R passes** against 4BVSnf with unchanged 1e-12 tolerance and exact
attributes/types/identifiers/missingness; original store heads remain unchanged
(6,190,518,272 B peak/zero swap, 135.911294 s). The 96 no-regression, scientific
and output conditions are satisfied. **300-sample UCTKOQ stopped at the early memory guard during first QC**, frozen
x3Dkg2/retained OQADra: 24,000,000 rows, 3,475,378,370-byte report, native outputs,
profiling off, unchanged 16 GiB hard/15.5 GiB early/zero swap/3600 s. Member three
of three; instrumentation-invalid replacement unused. Capacity only: no finite
large control exists and the earlier guard failure is not a speedup denominator.
Loaded R/dependency versions and
native paths/sizes agree except callr's expected per-process temporary path;
native binaries were not separately hashed by those inventories.

Import/design boundary observation (startup through the first raw_data_s4
completed-state receipt, not parser-only): cumulative sampled peak falls from
7,735,402,496 to 6,341,308,416 B (18.02%), and payload hydrations from 12 to 3.
Candidate metadata setup records zero hydrations and the working-table release
before FASTA. Total boundary memory is not uniformly lower: last pre-receipt
current memory is 4.350 versus 4.701 GiB, anonymous 4.123 versus 2.504 GiB and
file cache 0.195 versus 2.147 GiB (control/candidate). The later 8.22 GiB peak
first occurs during protein-evidence QC. Do not infer an import leak or a
whole-workflow proportional saving from these phase observations alone.
The older same-build GMUgZW used an 18 GiB boundary and remains historical,
not the comparator: fresh elapsed is 31.90% higher and CPU time 30.95% higher.
Inputs/settings/loaded namespace versions match; memory-limit events are zero.
The reason for that slowdown is unproven. Proceed with the user-approved frozen
candidate comparison, not another repair or diagnostic run.

Final routing preflight **CqlBmD passes**, 1,960,251,392 B peak/zero swap.
The user's immediate-launch request was attempted on the unchanged 16 GiB
96-sample control route; the launcher refused admission because host memory
available was below 18 GiB. No large job started and the final member budget
remains unused. Small qualification is complete; resource admission remains held.

Current continuation: **x3Dkg2/J72ScQ passes** the focused design-attribute
correction and affected regressions, 2,437,623,808 B peak/zero swap. Two tibble
schema-lookup warnings are retained. The previous cEmHvm history mismatch remains
a failed qualification; one corrected 24-sample GUI plus independent output/
history checks is next. No large member has run.

**R4h60K/x3Dkg2 completes the corrected small GUI** (2,889,367,552 B peak,
zero swap, 360.603 s). Its diagnostic timing is slower across normalization,
correlation, export and DA than cEmHvm; the cause is not established and no
performance improvement is claimed. **2xuGFv passes all 16 saved states against
vOSqsa**, with strict attributes/types/identifiers/missingness and unchanged
scientific tolerance (1,067,405,312 B peak/zero swap). Corrected-member output/
native-session verification is next; final large members remain unstarted.

**VgMr8I passes corrected-member native/DA/Summary output verification**, peak
1,340,977,152 B/zero swap. Small science/output/lifecycle qualification is complete.
Next is the fixed 16 GiB final routing preflight and freeze; the three-member
96 control/96 candidate/conditional 300 budget remains unused. No performance
acceptance or default promotion yet.

User authorizes the exact missing identity_root fixture argument and one
gui-capacity recheck on unchanged StLiHl. All preceding failures remain; no
application changes, limit changes or new large-run allowance.

**StLiHl/BshfMC passes**, 2,206,388,224 B peak/zero swap, 365.100 s. Combined
with unchanged passing z3ohQR component tests, all six implementation slices
are verified. Proceeding to native lifecycle/ten-cycle and small full-GUI
qualification. This is implementation evidence, not large-workflow performance
acceptance; the 96/300 queue has not started.

**StLiHl/CKQ8r4 native lifecycle passes**, 2,021,834,752 B peak/zero swap,
82.799 s. All ten normalization/revert/style/export cycles leave zero working
records/cache payload bytes/open resources; recovery and source-relocated real
IQ reapply also pass. Proceeding to the retained 24-sample ykIZwA full GUI and
independent science/output checks before the final large queue.

**cEmHvm/StLiHl completes the retained 24-sample full GUI**, 16 states and
766 proteins, 2,917,175,296 B combined peak/zero swap. No full-history hydration,
working/native payload records, cache bytes or open resources at delivery.
This diagnostic small member is not a performance comparison. Independent
native/DA/Summary output and saved-history agreement checks are still required.

User approves the exact four-failure test/harness correction and one affected
gui-capacity recheck on unchanged StLiHl. No application or scientific changes,
raised limits or additional large-run attempts. Resume the existing integrated
qualification/final queue only after this check passes.

**StLiHl/3lNIfv fixes the four recorded failures but exits failed on one new
fixture omission**: the added native-launch assertion lacks identity_root.
Peak 2,199,879,680 B/zero swap, 347.900 s. Schema-only read/stream cleanup,
ordinary-reader error routing and installed/S4 compatibility pass. The exact
remaining one-argument test correction has been requested asynchronously; no
approval is assumed. Application code/build is unchanged and all jobs are stopped.

Current state: slices 1–5 verified; slice 6 implemented but held after its
correction allowance. StLiHl/KsYGhR fails on an Arrow-spy environment conversion,
the omitted-selector reader-harness path and two stale closure expectations.
Peak 2,196,250,624 B/zero swap; fresh installed/S4 compatibility passes 22 checks.
The prior StLiHl/z3ohQR checkpoint I/O fixture retains exact payload equality
and zero resident file pages after save/read for an 8,389,076-byte payload;
FASTA golden outputs and other selected scientific regressions passed there.
DIA-025 records the narrow test/harness-only continuation needed. No new small
full GUI or 96/300 member has run, no job remains live and no whole-workflow
performance claim is made.

User authorizes the bounded slice-4 test-only continuation: repair the remaining
mock cleanup lifetimes and make the passive fixture exercise normal Shiny stop,
retaining actual stage receipts. One focused recheck on unchanged 7icXAI under
the existing small-check limits; no new large-run allowance or runtime change.

**7icXAI/N56UPL passes** after that test-only correction: 2,383,831,040 B peak,
zero swap, 150.687 s. This includes actual design files, project replacement,
numeric conversion/native reference restoration and saved browser enter/stop
receipts. Slice 4 is complete; proceeding to downstream operation ownership.
No new full-workflow performance claim follows from this focused result.

**TepvoW/UQB2Vr passes slice 5**, 2,351,288,320 B peak/zero swap, 199.706 s:
downstream operation scopes, consumed-input releases, IQ staging and norm/style
handoffs pass with affected protein-QC/normalization/DA regressions. Its first
Kqo6KA run remains failed due to a source-only fixture missing package helpers;
the correction uses the actual installed observer. Proceeding to final named-file
I/O and bounded-header reads, not large performance testing yet.

User approves six bounded repairs from the immutable import-memory audit, with
focused checks as refactoring proceeds, no streamed-ingestion rewrite and **16 GiB
for both final workloads**. DIA-025 records the exact slices, compatibility
contracts and budgets. First work: intra-action lifetime tests and bounded scalar
observations, followed by operation-scoped ownership. No new large member yet.

The previous 300/500 jobs are stopped at their guards. Clean verifier **GjvReR
passes all 16 GMUgZW/v6Blcs states**, 6,113,869,824 B peak/zero swap; the previous
X9hWYX failed exit remains retained. No application performance acceptance or
default adoption follows from this repair authorization.

Implementation progress: **iBnD2I/NpkJgf passes** the shared ownership correction
and active-owner diagnostics (2,355,789,824 B peak/zero swap). The expected
pre-repair failure CufMze remains, including its corrected cancellation-mock bug.

The single-handoff import is implemented in **X7Y50A**, but its corrected check
**OL0bkG is failed**, peak 2,354,327,552 B/zero swap. Actual release-before-FASTA,
exact table/summary, one imported payload and unchanged-history assertions pass;
the only failure is a scalar callback-observation flag isolated from updates by
Shiny testServer's cloned test mask. WDaM1C's preceding partial `$data` matching
failure is also retained; explicit `data = NULL` fixed that handoff ambiguity.

The separately approved fixture-only recheck **FOoJQC/X7Y50A passes**, peak
2,360,307,712 B/zero swap, with exactly one FASTA invocation and every previous
release/parity assertion retained. Slice 3 is complete; design/reference repair
is next. No new 96/300 run or large-workload performance claim yet. Historical
failures above remain retained.

Slice 4 first check **dpwmtE/8KPdAr failed**, 2,309,709,824 B peak/zero swap.
Its one correction addresses leaking historical test mocks and native export's
rejection of the bounded summary embedded directly in reference roots. The
summary now uses an existing RDS metadata record; no source assay is added to
the resident metadata cache. Same focused selector is required before advancing.

Correction **7icXAI/YmalER also failed**, 2,405,232,640 B peak/zero swap,
131.794 s. The numeric-conversion/native-reference archive test now passes, as
do import and core ownership tests. Additional historical mocks in the design
characterization tests still leak fake S4 constructors and the artifact writer
into later tests (five errors); a separate passive-browser receipt assertion
fails. No demonstrated numeric mismatch, but the suite is not qualified and
the new project-replacement assertion is not reached. DIA-025 records exact
remaining test-only work. Paused at the agreed correction limit; no new large
run, no running service and no new performance claim.

## Retained requested 300-sample attempt, revision 1.36

User explicitly requests **300 samples** at the same 16 GiB allowance, using
pre-approved commands while AFK. Keep xz9cb6 and other data/scientific settings:
24,000,000 precursor rows. One guarded attempt after the existing small route/
chunk checks and bounded generation; no automatic retry, further resizing or
cap increase. Exact scope and stopping conditions are in DIA-025 revision 1.36.

**05uICS passes** both stress-route safety/chunk checks (2,169,593,856 B peak,
zero swap). **OQADra** prepares 24,000,000 rows in chunks <=199,200 rows, report
**3,475,378,370 B**, peak 3,739,131,904 B/zero swap. Exact input/build/harness hashes
are frozen in DIA-025. **8WgDVe / xz9cb6 is running** using pre-approved commands;
no application package edits or change to the prior 500-sample failure.

**500-sample iGvvnA stops at the memory guard during import/setup**, before
raw_data_s4 or any QC checkpoint. Guard observation: **16,916,217,856 B
(15.755 GiB), zero job swap**, host available 13,413,482,496 B. Journal confirms
the controller stopped it after 99.292 seconds; no completed-workflow timing or
speedup denominator. No RUV/Pearson/export/DA result. Last complete sample has
9,901,445,120 B anonymous and 6,681,321,472 B file memory; asynchronous components
are not an exact peak partition. Exact import-function attribution is unavailable
without profiling. All inputs/artifacts remain; the 300 case does not replace it.

## Retained requested 500-sample stress, revision 1.35

User requests one 500-sample run with a 16 GiB allowance. Keep xz9cb6, 10,000
proteins, four peptides/two forms, seed and scientific settings; this is
40 million precursor rows and approximately 5.8 GB report, not numeric-matrix
bytes. Actual input size/hashes will be frozen after bounded generation. The
live 96-sample output verifier DQtNGk and saved-history checks finish first.

Add only explicit 500-case routing/guards to the existing runner and fixture
manifest, with small checks before generation. Full-job limits: 16 GiB hard,
15.5 GiB early, zero swap, 3600 seconds, >=18 GiB start headroom and >=2 GiB host
stop. One candidate, no larger resident control or automatic rerun. This is a
capacity stress, not a matched performance claim. DIA-025 records the exact scope.

User prioritizes the stress launch over remaining verification cleanup. **S2bEmG
passes** the 500-only safety/routing/chunk checks (1,980,743,680 B peak/zero swap;
one existing upstream ruv plotting deprecation warning). **570cC2** prepares the
complete input: **5,792,259,275-byte report, 40,000,000 rows, 500 samples** in
at most 200,000-row chunks; peak 6,027,583,488 B/zero swap. Exact input and harness
hashes are recorded in DIA-025. No application code or installed build changed.

**iGvvnA / xz9cb6 is running**. Live systemd limits independently read back as
16 GiB MemoryMax, zero MemorySwapMax and one-hour RuntimeMax. Report, FASTA and
design each have zero resident file-cache bytes at entry. No harness edits while
the job is live; retain a valid guard stop without repeating or shrinking input.

96-sample verification status: **DQtNGk passes** native restore/DA/heatmap/Summary,
8,949,923,840 B peak/zero swap. **pCp5qY passes** the retained small history pair.
**X9hWYX** writes a passing 16-state GMUgZW/v6Blcs comparison but subsequently
exits with a parsing error caused by editing its shared runner while it was live
(6,213,025,792 B peak/zero swap). Retain both comparison and failed exit; clean
saved-history recheck is deferred at the user's request, not silently certified.

## Retained corrected full-GUI recheck, revision 1.34

User authorizes the remaining full journey. Run one xz9cb6 candidate on unchanged
nDqUl6 (1.11 GB report, 96 samples), then verify outputs/restore and compare all
historical states with v6Blcs. No runtime or harness edits. Reuse zi8nIV lifecycle
and Lb1gdW/TNjX32 mechanism evidence. This is the third full timed member overall,
not release of the old conditional second pair. No automatic repeat or adoption.

Initial host available memory is 22,622,236,672 B; no proof job is running.
Keep 18 GiB hard/17.5 GiB early/3600 seconds/zero swap, >=20 GiB launch and >=2 GiB
host stop. Same frozen input/cache/settings/boundary and profiling disabled.
Corrected RDB and unchanged runner/driver bindings are recorded below and in
DIA-025. Full target remains <=1254.110384 s and <=10,530,879,897.6 B against
qualified v6Blcs. Complete results and remaining QC limits must be reported;
an export-only gain does not waive the full-workflow gate.

Initial launcher **sOMHCc** rejects the corrected RDB before package/input load:
1.278155 s, 66,011,136 B peak, zero swap. Retain it as an infrastructure-invalid
launch. The single binding correction appends only xz9cb6's verified fingerprint
to the native build allowlist in the workload JSON; all other fields and the
shell/R runner/driver remain unchanged. One replacement launch follows successful
binding/unchanged-field validation. No application-performance result is inferred.
Validation passes: reverting only that allowlist field recovers the exact prior
manifest hash; both declared builds match and an unknown hash is rejected. New
manifest SHA-256: 5088f056b487737c26247fd2d628de144e60d5b022e3bc55a98f86306857169e.
Replacement **GMUgZW / xz9cb6 completes** through DA/Summary and teardown:
**1249.193656 s / 10,466,906,112 B (9.748 GiB) peak / zero swap**. Against v6Blcs,
elapsed is 10.353% lower and peak 10.547% lower; margins inside the 10% thresholds
are only 4.917 seconds and 61.010 MiB. Against the failed native ARdJpI, peak is
25.555% lower but elapsed is 5.778% higher. Independent output/history verification
is pending; do not yet label science or repeated-pair acceptance passed. Final
ownership reports 16 states, no full-history reads, zero working/native/cache
payloads and zero open resources. sOMHCc remains retained separately.

## Retained export cache-overlap correction, revision 1.33

First native candidate **ARdJpI / 9SWlGu** completes at **1180.963283 s /
14,060,003,328 B peak / zero swap**: 15.249% faster but **20.161% higher peak**
than v6Blcs. This fails the memory gate; no second pair or default adoption.
Independent candidate output/history verification remains pending. During native
packaging anonymous memory stays about 6.89 GiB while file memory rises from
2.54 to 5.28 GiB; the final 13.094 GiB peak follows during downstream work.
The archive is 2,075,473,920 B, staged dependencies 571,104,348 B and checkpoint
members 1,503,941,584 B. These file sizes are not additive RAM attribution; no
per-file trace at the full-workflow peak exists. Earlier QC also exceeds the
whole-workflow target, so an export correction alone cannot establish acceptance.

The export-only correction is implemented and verified on **xz9cb6**. It uses
synchronous 4 MiB archive/copy writes with named-file cache-release advice and
drop-behind SHA-256 reads, retaining all files, full integrity checks and recovery.
No global flush or excluded memory. GNU dd's cache advice concerns persisted
pages, not dirty or RAM-backed files ([GNU coreutils documentation](https://www.gnu.org/s/coreutils/manual/html_node/dd-invocation.html)).
This is a Linux/disk-backed export result, not a portable cache-size guarantee.

**zi8nIV passes** the full native-session file, including ten successive cycles,
all normalization choices, manual RUV, Pearson Apply/Skip/reapply, portable IQ,
corruption rejection and process-loss recovery. Peak 2,014,674,944 B, zero swap.
The 33,554,555-byte non-aligned source, staged copy and 33,761,280-byte archive
each have zero resident file pages at the post-export check; exact bytes restore.
Cycles 1/10 RSS: 994,693,120 / 997,130,240 B; maximum metadata 604,392 B;
all completed-cycle working/native/cache payload and open-resource counts zero.

One representative recheck **Lb1gdW / xz9cb6** uses unchanged IVL1gM input/head,
dependencies, operation boundary, guard and cold-source setup (252 files,
1,410,622,924 B, zero resident bytes at entry). Profiling is off. Compared with
the preceding native implementation, not the full GUI:

| Native export | Whole-service seconds | Combined peak bytes | Operation seconds | Copy/reopen seconds | Unattributed seconds |
| --- | ---: | ---: | ---: | ---: | ---: |
| Before: lSCYZW / 9SWlGu | 51.282888 | 4,747,264,000 | 34.293 | 1.661 | 15.328888 |
| Corrected: Lb1gdW / xz9cb6 | 20.420501 | 3,606,306,816 | 12.668 | 1.865 | 5.887501 |

Observed **24.034% lower peak and 60.181% lower elapsed**, both zero swap.
The complete 765,726,720-byte archive is the same size, not a reduced inventory.
Independent **TNjX32 passes all 16 restored historical states** at 1e-12 with
exact type/attribute/missingness checks and unchanged source heads; verifier peak
4,579,127,296 B/zero swap, outside performance timing. CPU falls from
47,324,432,000 to 16,364,486,000 ns. Physical reads increase from 2,662,907,904 to
5,405,999,104 B; writes are 3,665,428,480 / 3,680,870,400 B. Releasing cached pages
trades reuse for additional disk reads; this result does not establish speed on
slower disks. No representative per-phase memory trace was added to the frozen
isolated mode. Small-file residency evidence and whole-service peak are distinct.

Bindings: corrected RDB 74a7dc9b50821e49e1a6d8c9a4366bdcc879997d61d06a9876acd55387404a08;
native helper fccf3e64979ed09ee6ddef19123b6b7654d853e07ceb3e80f011921e30d2c93c;
native test f29ccd5453368b52b6fa929aee58c2a94fad6a6092949a3cbaca9b227fa25212;
mechanism test ba1d070f5a3f3434967f977cfcd2c76c889ce5f8cec6ecf7e708de6d45f48626.
Runner/driver/spec hashes below are unchanged. The mechanism test changes only
the exact admitted corrected-build identity, not its workload or thresholds.

**Bounded correction complete. Full-GUI memory acceptance remains failed/unproven
for this new build; no full rerun, second pair or default activation occurred.**
The next full-journey decision must retain the earlier QC limit and independently
qualify the changed build; it cannot reuse ARdJpI as if unchanged.

Retained initial corrected-build check GmLOj9/VcLDXX fails at 1,836,204,032 B/zero swap:
the new cache fixture was on tmpfs /tmp, and GNU dd's version banner lacks the
word GNU expected by activation. Byte round-trip and I/O failure propagation
pass. The single bounded correction moved only the cache fixture onto the proof
disk and fixed the banner check. GmLOj9 and all earlier failures remain retained.

## Retained first matched full-journey pair, revision 1.32

Host headroom is 20.29 GiB after the user closes applications, above the fixed
20 GiB launch requirement. First timed full control **v6Blcs / B4Errr** has
started on frozen nDqUl6 (1.11 GB report, 96 samples). Next is integrated 9SWlGu;
no silent removal of the measured roll-up regression. Profiling is off, required
GUI/DA/Summary/output work and teardown remain included, and the exact saved design
is checked before QC. Limits: 18 GiB hard/17.5 GiB early, 3600 seconds, zero job
swap, 2 GiB host safety stop. Second pair is conditional on the first passing
science and both >=10% gates. No default adoption is claimed.

v6Blcs completes at **1,393.455982 s / 11,700,977,664 B peak / zero swap**.
Independent ZEEPON passes complete output checks: 16 states, 8,152 proteins,
96 samples; verification peak 9,460,670,464 B, zero swap. Candidate **ARdJpI /
9SWlGu** started with the same frozen boundary after a 21.74 GiB headroom
check. The first pair requires candidate time <=1,254.110384 s and peak
<=10,530,879,897.6 B, with science/output/guard checks. Two timed members are
used; the completed result above fails the memory gate, so no second pair follows.

Same file-specific cold-data-cache setup is applied to the exact report/FASTA/
design inputs before each timed member. No file-content writes, payload warming
or global flush. This adds only a three-file selector to the tested cache helper;
application builds and driver are unchanged. Full-pair hashes:

- shell: 1fbf833409a9aea82a8a300742a7d3de67a203d815aaad7c42d13b4c7de4cd7a
- R runner: 27f9711a4eb7aed0109c3b16bc54002d896502a3b7abe17a78117b1fca491bdc
- driver: 860c071cdb8c436d1e29958ee03da4ed0f0eb8759c4139a15d2fed41bea7d3e7
- workload spec: e9dbe37fc163002e55bb7d2340194e701f301eeff084eb99c20a8d24f6aa22ab

### Completed isolated comparisons — not full-GUI acceptance

Both pairs used the same frozen IVL1gM 556 MB / 48-sample history, matched cold
data pages (252 files, 1,410,622,924 B; zero resident bytes at entry), identical
recorded R/namespace versions, fresh installed services, no profiling and zero
swap. All four outputs pass independent scientific verification.

| Mechanism/member | Whole-service seconds | Peak bytes | Operation seconds | Copy/reopen seconds | Remaining service seconds |
| --- | ---: | ---: | ---: | ---: | ---: |
| Roll-up control RCfz3J / B4Errr | 35.880867 | 5,900,009,472 | 24.412 | 1.674 | 9.794867 |
| Roll-up candidate He1p1l / SJLZhY | 49.694495 | 5,844,054,016 | 32.637 | 1.826 | 15.231495 |
| Export control TNKt0l / B4Errr | 172.479533 | 5,770,969,088 | 154.094 | 1.996 | 16.389533 |
| Export candidate lSCYZW / 9SWlGu | 51.282888 | 4,747,264,000 | 34.293 | 1.661 | 15.328888 |

Roll-up is **38.5% slower for only 0.95% lower peak**; its operation time also
regresses, so startup alone cannot explain it. Native export uses **70.3% less
elapsed and 17.7% less peak memory** in this isolated case. Remaining service
time includes startup, bookkeeping and teardown and stays unattributed; it is
not discarded. Peak includes cgroup-charged file/kernel memory, not only R heap.
Do not add these independent peaks or sum the timings into a claimed GUI result.

PN0Lp1/220L67 verify roll-up science; H143Vz/9dflLm verify all 16 export-history
states at 1e-12 plus type/attribute/missingness checks. Oracle peaks are separately
1,420,124,160 / 1,407,205,376 / 6,225,272,832 / 4,899,766,272 B, all zero swap;
these are correctness processes, not production restore benchmarks. Native archive
765,726,720 B; legacy timestamped and latest files 609,339,926 B each. All earlier
failures remain. Outstanding fresh-project browser coverage is not waived.

## Retained revision 1.31 — isolated performance testing

The approved warning-only fixture correction is complete. **yGfILk passes** all
native-session checks, including ten successive cycles and archived PNG equality,
portable real IQ reapply and three process-loss recoveries. Peak 1,991,540,736 B,
zero swap, no failures/errors/warnings/skips. Cycle 1/10 RSS is
990,232,576 / 993,783,808 B; maximum metadata 604,392 B; all working/native/cache
payload and open-resource counts are zero after every completed cycle. These are
small fixture observations, not large-workflow gains.

The next measured item is one serial control/candidate pair each for roll-up
and session export on the frozen IVL1gM representative history (556,028,391-byte
report, 3,840,000 rows, 48 samples), under 8 GiB/900 seconds/zero swap. Control
B4Errr; roll-up-only candidate SJLZhY; session candidate 9SWlGu. Scope, source
head/build bindings, complete service boundary and the established file-specific
cold-cache setup are recorded in DIA-025. No profiling; all actual reads, owned
copy/reopen, preparation and publication remain measured. Roll-up's isolated
boundary excludes GUI rendering, which remains mandatory in full-journey timing.

Adapter preflight p1U8h5 fails before operation execution (test registration order)
and is retained. Its single bounded correction is test ordering only. **PNB4hg
passes** both routes and small scientific/native readback at 1,620,017,152 B peak,
zero swap. Frozen executable harness hashes before representative members:

- shell runner: eacb18aac306a9ebd7fd4bb0b4e26e1ea71a07e2222c3ac43ca6e92e682c4b09
- R runner: 802be4460648f167fd08b6188c43faf6acdb8b6f60c12a6d129d1f3f8a40dbba
- driver: 860c071cdb8c436d1e29958ee03da4ed0f0eb8759c4139a15d2fed41bea7d3e7
- mechanism tests: e9ff1a808645fc5bc2794b4eac337f12bdc0a64fdc2c79613d54236980274244

GNU tar 1.35; member receipts retain installed build/R/dependency identities.
Scientific oracle work is separate. No percentage gain or default adoption is
claimed before matched results. Full 96-sample and canonical guards remain as
approved, and remaining fresh-project GUI coverage is not waived.

## Retained revision 1.30 — native lifecycle checks

Continuation uses unchanged installed 9SWlGu and the existing 8 GiB/900-second/
zero-swap runner. 6BFgf2 checks ten Shiny-owned normalize/revert/style/native-
export/restore cycles, all four normalization choices, Pearson Apply/Skip/reapply,
manual RUV, real IQ reapply after source-project removal, and journal recovery
after fixture-owned child-process loss at three publication boundaries. The
existing small automatic-RUV journey/history evidence is reused. Per-cycle RSS
and resource observations are written to native-cycle-checks.json inside the
guarded member. No production changes or performance/default-adoption claim.
6BFgf2 is retained as failed at 2,204,499,968 B peak, zero swap: the new cycle
fixture wrongly required a warning for scale normalization, and the runner's
native selector lacked the browser binding required by an existing browser test.
The one fixture correction permits absence of the known warning and adds the
existing browser setup; no test is removed and production remains unchanged.

**IKol3q is also retained as failed**, peak 2,207,150,080 B (2.056 GiB), zero swap.
Browser setup, real IQ reapply after fixture source removal and three process-loss
recoveries pass. Nine normalize/revert/style/native-export/restore cycles complete.
Full-history reads, working/native payloads, payload-cache bytes and open resources
are zero after each; maximum metadata is 543,336 B against the unchanged 64 MiB
ceiling. Process RSS after cycles 1 and 9 is 1,409,597,440 / 1,409,986,560 B. These
are small-fixture test-process observations, not large-workflow gains.

The tenth cycle passes real manual-RUV numerical comparison and Pearson reapply,
then fails the fixture's still-too-narrow "Removed 1 row" plotting assertion:
the existing geom_bar histogram emits "Removed 2 rows". That expectation failure
destroys the mock session before export/restore; ten successful cycles are not
established. The fixture-correction allowance is spent. DIA-025 records the exact
test-only correction/recheck decision needed; no scientific or production change
is proposed. Helper-level portable restore is not proof of fresh-project browser
reapply. Remaining lifecycle, performance and default-adoption gates stay open;
all services are stopped and no large performance member has run.

## Retained revision 1.29 — small history comparison passes

User approves normalization of only data@args$daAnalysisUI$timestamp's wall-clock
value in comparison copies. Stored provenance is unchanged; field presence/type,
other timestamps and scientific parameters remain protected. Focused **wAYztS
passes**, peak 893,272,064 B, zero swap, no failures/errors/warnings; four explicit
historical/opt-in skips. No production/package or GUI change.

**odghtX passes all 16 saved historical states** from WD1nXG against vOSqsa,
including unchanged 1e-12 scientific tolerance and exact table/column type,
attribute and missingness checks. Peak 1,163,366,400 B, zero swap. The receipt
records the exact additional comparison-only clock normalization. Reuse verified
WD1nXG/dvY8S2 outputs; no new GUI/scientific execution, correction retry or large
performance member. Retain jpRHRI and every earlier failure.

The small native output/restore and cross-build history-science checks are now
clear. Remaining native lifecycle/ten-cycle coverage and isolated/full-workflow
performance gates are still open. No performance gain or default adoption is
claimed; DIA-025 stays in_progress. All services are stopped.

## Retained revision 1.28 — design-driver correction

User approves continuation after the specifically requested design-driver repair.
Send selected inputs before the Assign event, wait for the production technical-
replicate summary, and verify saved sample/group/replicate assignments before QC.
Use unchanged 9SWlGu/ykIZwA, one small recheck and its independent output/history
verification. No production scientific edits, guard relaxation or large-run use.

New event-order and design-mismatch checks pass in bRCXp9, including rejection of
retained SfeXyz and acceptance of vOSqsa's design. The broader check fails on an
old passive-browser test double that lacks the launcher's session_format argument.
Retain bRCXp9 (peak 1,441,439,744 B, zero swap). The single bounded fixture
correction adds that argument and checks the unchanged legacy default. Recheck
**xc4ikk passes**, including the real passive-browser lifecycle, new design
guards and affected RUV/export/recovery checks. Peak 2,088,407,040 B, zero swap;
four historical skips and the known upstream ruv warning.

**WD1nXG completes / dvY8S2 independently verifies outputs.** The corrected
driver confirms all 24 sample/group/biological and technical replicate assignments
before QC. GUI peak 2,923,327,488 B (2.723 GiB), zero swap, service-tree interval
288.601472 s. Verification peak 1,264,246,784 B, zero swap; 16 restored states,
766 proteins/24 samples, archive integrity, TSV/RDS, DA/heatmap, Summary and
RUV/Pearson traces pass. Full-history reads remain zero, hydrations 41, and final
working/native/cache payloads and open resources are zero.

Three scientific TSVs match retained vOSqsa byte-for-byte:

- imputed peptide: 2ef8dd544e2390d0996b8b5f45bfca63f1f3e66b103dade9e3ebc3214f0a5911
- normalized pre-RUV: 257394d670500333447883045b242714d65a5240e5ef28d2c839574cd6054b29
- final post-RUV/Pearson: 746f2e302e4476a65663e36a1225578608856a289462a33f8855a152d89ae0b3

**Formal history comparison remains held.** jpRHRI fails after preceding states
pass, reporting only args$daAnalysisUI$timestamp in the last state (peak
1,162,928,128 B, zero swap). Production records Sys.time() there; the current
comparator normalizes record/audit clock metadata but not this nested UI field.
This is not a demonstrated numeric difference, but it must not be silently
ignored. The used correction budget requires the narrow additional verifier-only
decision recorded in DIA-025. Reuse these saved files; no GUI/science rerun is
needed to check that correction. No large member has run, no performance gain is
claimed, and the production default remains off. All earlier failures remain.

## Retained revision 1.27 — native-session and roll-up implementation

User approves the native single-archive/import-only-legacy plan and conditional
DIA-IQ default adoption after both >=10% elapsed/combined-peak gains. Preserve
B4Errr as control and nDqUl6 inputs. Workload-specific limits and the complete
bounded implementation/run budgets are in DIA-025. No default has changed.

First roll-up slice **SJLZhY** installs (peak 1,235,095,552 B, zero swap).
Focused **ryg3c9 passes**, peak 1,546,641,408 B, zero swap, no failure/error/warning;
two historical withdrawn-IQ fixtures skip explicitly. Existing science/state/
resource/recovery checks and new typed collection, empty/oversized result,
envelope identity/corruption, reopen/revert and no-full-predecessor-read assertions
pass. The 200,000-row allocation fixture records two 1,600,048-byte allocations
(driver plus destination), not a full-column allocation per batch. This is
correctness/allocation evidence, not the pending large performance decision.

SJLZhY RDB SHA-256:
2aead2cd7912c9303331e2a5a0bb67f599583b40871745643e0b88fdaf0c8f7c.
Only the roll-up collector/lifetimes, generation-bound envelope/state binding and
its Apply caller differ from B4Errr. Envelope/projection preparation is persisted
as auxiliary metadata, not an additional scientific state. Audit-on and older
preparation retain full-state fallback. No correction used, large members zero.
Subsequent native archive helpers TIJS4W/QIWEaM pass reference-only packaging,
fresh restore after the source store is unavailable, retained branches/revert,
external integrity rejection and journaled output-only retries. Peak 1,430,130,688 B,
zero swap. GUI integration r675bH/Ppy22Y then fails because its output inventory
used a singular RUV metadata filename instead of the producer's plural filename;
the failure is retained (peak 1,440,022,528 B, zero swap).

The one bounded correction, installed **9SWlGu / PQlDmp**, passes the corrected
inventory, native GUI export/owned DA load, artifact restore/conflict rejection,
old RDS import and explicit DIA-IQ activation checks. Peak 1,437,110,272 B, zero
swap, no failure/error/warning; two historical withdrawn-IQ tests skip. RDB:
b6c7ef9aff9546d2dbb758931be6d7095277dc613919ca4f50e7dd9190f6c187.
Expanded **3Wzdlu passes** the format-aware harness/verifier and affected module
contracts: peak 1,489,928,192 B, zero swap, two historical skips and one known
upstream ruv plotting warning. Runtime 9SWlGu is not rebuilt.

**SfeXyz completes / 0N61pX independently verifies native outputs**, on unchanged
ykIZwA files: 16 restored states, 757 proteins/24 samples, archive integrity,
TSV/RDS, DA/heatmap, Summary/publication and real RUV/Pearson traces. GUI peak
2,969,505,792 B (2.766 GiB), zero swap, service-tree interval 367.413806 s;
verifier peak 1,249,292,288 B, zero swap. Final adapter has zero full-history
reads, 30 hydrations, no retained working/native/cache payloads or open resources.

**Cross-build science is held:** saved-history comparison **tem6Xz fails** at the
first nonempty state's design metadata (peak 1,178,021,888 B, zero swap). SfeXyz's
saved design differs from the declared input/control in the first four samples:
incorrect biological replicate numbers for runs 1/2/4, and missing technical
replicate numbers for runs 3/4. Later design rows match. This precedes QC and
explains why the 757/24 result cannot establish agreement with vOSqsa's 766/24.
Driver input/assignment synchronization is the suspected source; not yet proven.
Do not attribute this mismatch to altered scientific calculations or dismiss it
as numerical tolerance. Nominal elapsed exceeds vOSqsa's 244.917067 s, but the
executed designs are not matched and these diagnostics are not performance gates.

All members are retained. No large performance member ran; default remains off.
The implementation correction allowance is exhausted. DIA-025 records the exact
additional decision needed: bounded design-driver synchronization/pre-QC design
validation, focused checks and one small recheck, without changing runtime
science, inputs, limits or the remaining large-run budget. Remaining lifecycle,
ten-cycle and performance acceptance are explicitly incomplete.

## Retained revision 1.26 — 96-sample diagnostic execution

User authorizes the proposed approximately 1 GB / 96-sample case after the small
gate passes. Frozen B4Errr and the verified harness/settings are reused without
code edits. One serial GUI member under 18 GiB hard/17.5 GiB early/3600 seconds/
zero swap and existing host safeguards; preparation retains 8 GiB/900 seconds and
bounded streaming. No repeat after a valid failure and no performance/scaling
claim. The 192-sample case remains conditional. Preflight confirms matching
build/harness hashes and no live jobs before execution.

Preparation **nDqUl6 completes** under the stricter preparation limits, peak
1,367,638,016 B (1.274 GiB), zero swap. Actual report size 1,112,120,153 B;
7,680,000 precursor rows, 10,000 proteins, 96 samples, four peptides/two forms.
Largest generator chunk 199,680 rows. FASTA is 710,000 B and design 5,710 B.
The unfiltered numeric protein matrix would contain 7,680,000 B; this is a
1.11 GB report, not a 1.11 GB numeric matrix. Surviving dimensions are measured
during execution. Input hashes:

- report: 3d670aa2567263e9ae90e2966b570e4d65200e9d5adac7214b11eb705a01032a
- FASTA: 352b56704c330f3e59646101df2c52a469e5ebe17da652fdd0d3ea85e672f6d4
- design: a8ec88469f582dcfd12b1715ccf20bbcd49a83ead3d84bfc05543a25e446040f

**Z8rE5w / B4Errr completes and independent NWub6d passes.** The full 96-sample
GUI journey includes automatic RUV, Pearson Apply, automatic files, complete
legacy export, DA/interactive views/heatmap, Summary/publication and teardown.
Service-tree elapsed is 1,228.357202 s; whole-job peak 11,900,006,400 B
(11.083 GiB); zero swap. This is an instrumented stress diagnostic, not a timed
performance member. It exceeds the original 8 GiB/900-second acceptance limits;
no matched resident control, percentage improvement or broader promotion follows.

Independent **NWub6d** verifies 16 historical states, 8,152 final proteins and
96 samples, complete legacy restore, TSV/RDS and native DA/heatmap agreement,
Summary/publication files and complete RUV/Pearson trace coverage. Verification
alone peaks at 9,708,539,904 B (9.042 GiB), zero swap, separate from application
completion. Both complete legacy files contain 1,218,285,069 B and share SHA-256
7c419822a24362469c86906198ffcf9ddd433413a4080aa794b8db489f01e09a.

### Where memory rose

These are sampled **whole-job current-memory maxima during each phase**, not
incremental allocations or amounts to add together. Boundaries use successive
Apply helper entries, then the recorded normalization/correlation/export and
DA/Summary action boundaries. Adjacent rendering/persistence work is included.

| Phase | Sampled peak GiB |
| --- | ---: |
| Startup, import and design | 8.073 |
| Peptide q-value | 8.572 |
| Protein/peptide filtering | 9.084 |
| Precursor roll-up | 9.544 |
| Peptide intensity / sample / replicate | 7.060 / 7.028 / 7.077 |
| Peptide imputation | 8.531 |
| IQ and protein QC | 8.632 |
| Normalization, automatic RUV and composite | 7.499 |
| Pearson Apply | 7.505 |
| Complete legacy export | 11.082 |
| DA, Summary and teardown | 10.113 |

1. **Early wide-table QC is a measured application cost.** The raw 7,680,000-row
   table occupies 936,895,096 B by object.size; q-value output has the same row
   count, 22 columns and 1,213,376,312 B; protein/peptide output has 24 columns and
   1,274,816,600 B. These sizes are not unique resident heap. Across
   runPeptideQvalueApplyStep, application RSS rises 2.499 GiB; across
   runProteinPeptideApplyStep it rises 3.515 GiB. Those wrappers include their
   surrounding persistence/preparation, not only the scientific filter.
2. **Roll-up has two distinct costs.** The combined forward peak occurs at
   1790112652.313124 inside computeDIACapacityRollup, before it returns. Sampled
   current is 10,247,720,960 B, including 2,786,115,584 B (2.595 GiB) of file
   memory. After the SQL result returns at 1790112654.5563, the caller hydrates
   the whole protein_peptide_filtered S4 again. Across that hydration, process
   HWM rises from 7,425,404,928 to 7,923,150,848 B (+474.69 MiB) while the rolled
   result is already present. Avoiding that reload could remove working overlap,
   but alone is not proven to lower the earlier combined peak.
3. **RUVIII/Pearson did not reproduce the reported explosion on this case.**
   Automatic selection performs 20 native cancor calls on Y=96x8,222, then native
   RUVIII uses k=1, M=96x2 and 82 controls. Y is 6,314,496 numeric bytes. RUVIII
   entry/exit application RSS is 5,223,714,816 / 5,280,022,528 B (+53.70 MiB),
   with unchanged prior HWM. Post-RUV/Apply Pearson uses 8,152x96 (6,260,736
   numeric bytes); all four cor calls have unchanged entry/exit RSS. This is not
   zero allocation or proof about larger/different matrices. The composite's
   entry/exit RSS increase is 7.07 MiB here, unlike the small pilot's 699.64 MiB;
   already-resident allocator space can absorb work. No forced heap census was
   performed, so do not treat this as eliminated image allocation.
4. **Full-history legacy export creates the largest overlap.**
   collectProtNormExportSessionData/materializeDIAStates raises application RSS
   from 5,345,599,488 to 9,078,247,424 B (+3.476 GiB). The highest process HWM,
   9,287,225,344 B (8.649 GiB), first appears while hydrating the expanded
   protein_peptide_filtered history record. stageWorkflowLegacyExport then keeps
   that history alive through serialization, hashing and the durable serialized
   checkpoint. The whole-job peak is near the latter boundary. After the child
   scope/working owners release and existing collection runs, application RSS
   falls from 9,074,307,072 to 5,280,374,784 B (-3.533 GiB), before publication
   copying. Moving the same operation to another R process alone would not
   establish lower combined memory.

### Attribution limits and ownership

At the largest sampled cgroup current (11,899,535,360 B, time 1790113417.215148),
memory.stat reports anonymous 9,574,391,808 B, file 2,268,684,288 B (2.113 GiB)
and kernel 49,217,536 B. Near it, process sampling at 1790113416.698062 records
application R PID 3533986 RSS 9,074,307,072 / PSS 9,036,449,792 B; harness R PID
3532948 RSS 223,363,072 / PSS 192,430,080 B. Browser/helpers account for additional
memory. RSS/PSS and cgroup counters share mappings and are asynchronous; they are
not an exact additive partition. Unmeasured within-function allocations, allocator
retention, shared pages and the unsampled peak difference remain unattributed.
The harness R process's largest sampled RSS is 273,395,712 B, not the multi-GiB
application history. No continuous allocation/CPU profiler or diagnostic GC runs.

There are 8,948 cgroup samples, maximum gap 0.550351 s; process-snapshot maximum
gap is 1.68689 s. Cgroup sampling logs no error. Preserve two process-read gaps:
one denied harness smaps_rollup read and one exited PID 3536493 status read.
The registered sleep inhibitor covers this GUI service. The 64 MiB trace limit
was not reached and all required scientific entry/exit counts match.

Final adapter delivery has 16 states, one full-history read, 57 hydrates,
zero working/native/cache payloads and zero open store resources. This describes
adapter ownership, not every object retained by GUI panels or the R allocator.
The existing scientific/storage methods and limits were unchanged. B4Errr runtime
and all harness hashes match the verified bindings recorded above/below.

Decision: the 96-sample diagnostic and verification are complete. All services
are stopped; one large input/member used, no large repeat or runtime change. The
192-sample approximately 2 GB case is still unstarted. User now asks for review
of roll-up and export improvements; review does not authorize implementation or
a session-format change.

### Follow-up review requested by the user (not implementation authority)

The existing roll-up already executes its scientific aggregation in DuckDB.
collectWorkflowRollupQuery nevertheless retains every 10,000-row result batch
and then rbinds them; computeSavedWorkflowPrecursorRollup holds aggregated and
output representations through R joins/sorting, and closes the DB scope only on
return. These are concrete lifetime/collection candidates, not yet separately
measured allocation causes. The caller then demonstrably hydrates the full
predecessor despite replacing its peptide_data with the already-computed result.
An owned non-assay state description could avoid that assay reload on the default
audit-off path, while retaining complete slots/attributes/provenance and the
explicit audit-on behavior. That requires tests, not dropping validation or
changing any grouping, abundance sum, missingness, evidence or row order.

Legacy export already streams through saveRDS and reuses identical whole tables.
The observed residency increase is in history assembly, before serialization.
Earlier broader column-sharing changes were rejected on whole-GUI performance;
do not silently reintroduce them or claim their benefit here. Component-identity
reuse before decoding and shorter history/publication overlap are candidates
under the unchanged legacy contract, with no guaranteed peak reduction yet.
A self-contained native project bundle could instead copy the existing immutable
components/history manifest without hydrating all states, reconstructing ordinary
S4 values on demand. It can preserve scientific/R6 semantics but is not the same
legacy RDS format. Making it replace any required automatic legacy output needs
an explicit format/behavior decision. Offloading the same full hydration to a
worker alone does not establish lower combined memory.

## Retained revision 1.25 — Summary-verifier continuation

User authorizes the remaining verifier correction and recheck of already-complete
vOSqsa / B4Errr. checkCheckpointSummaryFiles now explicitly selects the existing
case-sensitive applied-RUV publication names; its caller propagates the already
validated diagnostic profile. The default Skip names and all other checks are
unchanged. Tests exercise both branches on retained small members, reject the
wrong route and reject an empty required table. No production edits, package
build, GUI rerun or additional large attempt.

**Small-run gate passed.** Guarded checks **IxY74y pass**, including both Summary
routes and the eight new routing/rejection assertions, peak 2,086,428,672 B,
zero swap. The four historical withdrawn-mechanism skips and known upstream ruv
deprecation warning remain explicit. No new failure or skip affects this repair.

Independent **IvJ02b passes** the same completed vOSqsa / B4Errr output, with
status passed_ruv_pearson_output_checks: complete 16-state legacy restore,
766 proteins/24 samples, TSV/RDS and native DA/heatmap agreement, required IQ/QC
and Summary/publication files, and complete scientific trace coverage. Coverage
records 20 paired native cancor calls, one paired native RUVIII call, four paired
Pearson calls and paired automatic-selection/Apply helpers. Verification peak
1,263,636,480 B, zero swap. The original GUI result remains 3,012,567,040 B
(2.806 GiB), zero swap, 244.917067 s diagnostic service interval; it was not rerun.

This completes the approved repair/reverification item. R/ still matches B4Errr
exactly; no additional runtime change was made. Prior Y9atKE/0saoyZ/Rdd0dy/SMWoRA
failures remain retained below. All services are stopped, and both larger inputs/
runs remain unstarted. Next is the existing 96-sample approximately 1 GB case;
the 192-sample case remains conditional on its completion. Limits and attempt
counts are unchanged. This is small-workload correctness/instrumentation evidence,
not performance acceptance or a large-data/scaling claim.

Frozen verifier: 500d1b7be55fee2838b5d5e8b722ce15611952a31b0140b3bc34940186d160a1.
Diagnostic test: 5b4d125b2a3df9f764c4bb71f26fe752fe42f8707db51d02ffe205f5a45e8e81.
B4Errr runtime and original vOSqsa measurement bindings below remain unchanged.

## Retained revision 1.24 — approved export-validation repair

User authorizes the narrow repair and one corrected small GUI recheck. Modify
only missing-argument handling in checkSavedWorkflowPayload, retaining complete
plot contents and resource rejection. Use fresh installed checks, native RDS/qs2
plot round-trips, rendered equivalence and affected checkpoint/export failures
before the GUI recheck. One bounded correction pass; larger attempt counts and
all limits remain unchanged. WxHeXq remains the retained control, and the repaired
runtime receives a separate frozen identity.

**B4Errr** installs successfully (peak 1,210,585,088 B, zero swap). Against frozen
WxHeXq, the only runtime difference is checkSavedWorkflowPayload: recognize the
missing sentinel and use non-dispatching environment extraction after the same
active/lazy-binding checks. No plot/session field is removed or replaced. The
implementation follows R's missing-sentinel/subsetting semantics; see
[R subsetting source](https://github.com/wch/r-source/blob/trunk/src/main/subset.c)
and [rlang missing-argument handling](https://rlang.r-lib.org/reference/missing_arg.html).

**MIQbCU passes** the guarded installed checks, peak 2,149,654,528 B, zero swap.
Both previously failing plot validations pass. The complete retained RUV result
round-trips through the actual RDS and qs2 checkpoint codecs: non-plot metadata
is identical, plotted data and classes match, and PNG SHA-256 is identical after
each round-trip. Missing values in environments/formals/lists are accepted;
connections, pointers, weak references, active bindings and lazy bindings remain
rejected, with no evaluation side effects. Scientific trace agreement and affected
checkpoint/export/recovery checks pass. Four withdrawn-mechanism checks skip
explicitly and the known upstream ruv aes_string warning remains. This is repair
correctness evidence, not performance acceptance; no correction pass was used.

The one small corrected GUI member **vOSqsa / B4Errr completes** on unchanged
ykIZwA inputs and automatic-RUV/Pearson-Apply settings. It retains the stricter
8 GiB/900-second/zero-swap small-only guard. Peak 3,012,567,040 B (2.806 GiB),
zero swap, service interval 244.917067 s. Legacy export, DA tables/static and
interactive volcano/heatmap, and all three Summary actions complete. Final
delivery records 16 states, one full-history read, 57 hydrates, zero native/
working/cache payloads and zero open store resources. Larger members remain unused.

Independent **Rdd0dy fails verification**, peak 978,141,184 B, zero swap: its
old required inventory and final-read paths still expect the RUV-Skip names
normalised_results_cln_with_replicates.tsv/.RDS. The real applied-RUV files exist
as ruv_normalised_results_cln_with_replicates.tsv/.RDS, as specified by unchanged
resolveProtNormCorrelationResultFilenames. This is verifier wiring failure, not
missing application output. Spend the one correction pass on that explicit
real-RUV branch, additionally requiring its five corrected QC PNGs; the old
Skip route remains unchanged. New checks test correct inventory selection and
rejection of a missing RUV figure. Reverify saved vOSqsa; do not repeat the GUI,
production science or package build. B4Errr remains unchanged.

**cdTci5 passes** those focused/affected checks, peak 2,129,915,904 B, zero swap,
with the same upstream warning and four historical skips. The corrected verifier
is 46932eaad2136e63ea168624bb329e02abffab318a2fab08fa94191f0aebe645; diagnostic
test 756e9591a5a87bafedfd399d19baa2d54cd313feefd51b0577817e8f1c9b4106.

**SMWoRA remains a failed full verification**, peak 1,263,706,112 B, zero swap.
Its IQ/QC inventory, complete legacy restore, final S4/TSV and native DA/heatmap
checks pass: 16 history entries, 766 proteins, 24 samples. Failure moves to the
Summary verifier, which still hardcodes normalised_results.tsv/.RDS. Actual
unchanged production copy logic (func_general_filemgmt_results_summary_helpers.R)
uses RUV_normalised_results.tsv and ruv_normalised_results.RDS when RUV applies.
Both files exist. Independent SHA-256 checks establish exact source/copy equality
for both files; this does not turn the incomplete whole verifier into a pass.
TSV/source-copy SHA-256: 746f2e302e4476a65663e36a1225578608856a289462a33f8855a152d89ae0b3.
RDS/source-copy SHA-256: 2a26b6018fd3ca60b7e1cd06fc9836d75289147a7ad54d1b9f01535da77692aa.

The first filename correction missed this downstream consumer. The one bounded
correction pass is now used. **Stop at this gate:** no further repair/recheck or
large input preparation/execution without renewed authority. The remaining narrow
action is explicit RUV filename routing in checkCheckpointSummaryFiles plus its
test, then verification of the same completed vOSqsa output—no application rebuild,
scientific change or GUI rerun. The original export-validator defect is repaired
and its codec/science tests pass; complete journey verification is still held.
All services are stopped, large attempts remain zero, and all failures are kept.

Repaired runtime RDB SHA-256:
8769e614e8b8126c3ed1ccaa6d491fbf857ffad65c988fa595e3cd33f2301a76.
Validator source: 7a7a15470f79b9c0d7ae51b5a6db0dda1073664ebd9bde4acc0cd3bec1a791e3.
GUI bindings: unchanged shell 1f6c3008d6e6d0b6031d14b2fdcc72eff688c3b4247dd8438c65a77162c26e6a;
R runner 27cf780be04af1ca09a559896f0c25aae9fcdabb8b389830e8215e440119783a;
driver 624a70d40bf6cba1e7384adf55fa51df731c33c68836808979aaa5a6325d9ae6;
spec 3ea949a48809598612552eaa89e194ae96d943dff815f77efd5a41c7c40f87a0.
Only the approved repaired-build hash was added to the diagnostic specification;
all recipe/settings/limit values and the original WxHeXq hash are preserved.
The affected test selector now also includes checkpoint/export/recovery checks.

## Retained revision 1.23 — sample-scaled RUVIII/Pearson diagnosis

User authorizes WxHeXq diagnostics at 96/192 samples and approximately 1/2 GB
report files, with the final 18 GiB total cap superseding 24/28 GiB. Runtime remains
unchanged. The original harness skipped RUV and correlation Apply; it did run
Pearson QC plots. This new diagnostic explicitly drives real automatic RUVIII and
Pearson Apply, and records actual downstream matrix dimensions rather than calling
a 1 GB report a 1 GB protein matrix.

The bounded surface, recipes, ordered checks/runs and stopping rule are in DIA-025
revision 1.23. Memory attribution separates application R RSS, other processes and
file-cache/kernel components, without forced diagnostic GC or retaining payloads.
The small pilot reaches RUVIII and Pearson Apply but fails session export. The
large runs remain held; no large-input conclusion or gain is asserted. All earlier
failed work and the restored WxHeXq result remain below.

Initial focused check **GEXm6K** is retained: native RUV/Pearson result comparisons
agree, but the new trace shape assertion fails. The observer was evaluated by
do.call during trace installation instead of scientific execution, so those
new trace events are not valid attribution evidence. Correct by quoting tracer
language arguments; assert that installation emits no scientific enter events,
then require actual call shapes and preserve the small trace artifact. This is
the one bounded correction pass. The upstream ruv aes_string deprecation warning
is retained. No small GUI or large input diagnostic has run.

Corrected **oGRIqt** passes 48 assertions, zero failures/errors and one upstream
ruv deprecation warning. Actual native-call events now show Pearson x=64x12,
RUV Y=12x64 and M=12x2; native RUV/Pearson/canonical-correlation results agree
with untraced calls at 1e-12, including missing-value input. **QXVwGF** then passes
55 assertions covering the added real-RUV legacy verifier and trace summarizer,
with the same single upstream warning. WxHeXq R/ remains unchanged.

Small preparation **ykIZwA** produces 192,000 rows / 24 samples / 1,000 proteins:
report 27,806,083 bytes, FASTA 71,000 bytes, design 1,450 bytes. Generation is
bounded and does not hydrate a full large report. Host available memory is about
19.1 GiB, below the 20 GiB launch headroom for an 18 GiB diagnostic. The small
pilot therefore uses the stricter existing 8 GiB/900-second guard, with an exact
small-recipe-only selector; large 18 GiB guards are not weakened. Guard-routing
check **qHim9A passes** (peak 2,112,888,832 B, zero swap, the same single upstream
deprecation warning). Large execution still requires sufficient host headroom.

### Small pilot result and export blocker

**Y9atKE / WxHeXq fails** at export, not a memory guard: whole-job peak
2,936,303,616 B (2.735 GiB), zero swap. The receipt's completed elapsed is null;
247.967451 s is only the failed main-process interval, not a speed denominator.
This uses ykIZwA's unchanged 27,806,083-byte report and the stricter 8 GiB/900-second
small-only guard. Twenty native canonical-correlation calls, automatic selection,
native RUVIII with k=1, and Pearson Apply at 0.5 execute. RUV receives Y=24x773
(148,416 numeric bytes), M=24x2 and 15 selected control indices. Post-RUV/Apply
Pearson receives x=766x24 (147,072 numeric bytes). All 24 samples remain.

The adapter records 16 committed states, including ruv_corrected and
correlation_filtered; zero working/native/cache payload records and zero open
store resources at final delivery. It records one full-history read during the
failed export. This is not proof that the application owns no other data:
workflowData still holds ruv_optimization_result as a list and, after Apply,
ruv_normalised_for_da_analysis_obj as a resident protein S4. Shallow owner records
are not a unique retained-heap census. Legacy session export, DA and Summary do
not complete, so no complete-output verifier runs.

The failure is `argument "..." is missing, with no default` while
serializeDIAOwnedSessionExport calls saveDIAStateValue for export-metadata.
The already-written native ruv_optimization_results.RDS remains intact: 567,070 B,
SHA-256 368b9db51a3ee57fdbfbd2e9c60076318d88218ba17d9d7e4e33e7115b26ed40.

Read-only reproduction **0saoyZ / WxHeXq** uses the existing guarded reader-check
runner, 8 GiB/900 seconds/zero swap, without another GUI or production edit. The
untraced payload validator rejects both that saved best_cancor_plot and a newly
constructed ordinary three-point ggplot. Every other saved RUV result field
validates. Both errors traverse plot layer `super` -> function environment ->
`env` -> `call_env` -> missing `...`. This reproduces a runtime payload-validation
defect independently of the new scientific tracing. The two portability
assertions fail and remain visible; this run is not described as a passing test.
Peak 2,112,786,432 B, zero swap. Scalar failure paths and the input hash are in
0saoyZ/ruv-export-validation.json; no call frames or payloads are captured there.

### Small-pilot memory attribution (not large-workload evidence)

| Observed boundary | Application RSS before -> after | What this establishes |
| --- | --- | --- |
| Automatic negative-control selection | 1,616,031,744 -> 1,616,400,384 B | No new observed process HWM in this small case |
| Native RUVIII | 1,616,457,728 -> 1,617,473,536 B | No new observed process HWM; 24x773 input only |
| Composite QC PNG loading | 1,617,866,752 -> 2,293,223,424 B at its highest image-load exit | A substantial rise during image decoding/assembly |
| Entire composite construction/save | 1,617,866,752 -> 2,351,493,120 B | Net RSS increase 733,626,368 B (699.64 MiB); HWM reaches 2,436,833,280 B |
| Pearson Apply native cor() | 2,356,699,136 -> 2,356,699,136 B | It runs after the rise, without a new observed HWM |

The composite interval is 1790108160.95183–1790108194.4698. Fifteen
loadProtNormImageAsPlot invocations decode/assemble its plot slots, followed by
savePlot. Source R/mod_prot_norm_qc_support_helpers.R reads PNGs as arrays, wraps
them in raster grobs and retains the panels together in patchwork. This names an
observed application-memory cost; it does not quantify unique live image bytes,
allocator retention or a proposed fix's benefit. No plot dimensions/formats,
science or production collection behavior were changed.

Peak sampled cgroup current is 2,935,734,272 B at 1790108191.928847, inside
composite savePlot: anonymous 2,740,662,272 B, file 163,487,744 B, kernel
23,896,064 B. Near that time (1790108192.051580), application R PID 3326309 has
RSS 2,350,952,448 / PSS 2,310,274,048 B; harness R PID 3326037 has RSS 220,360,704 /
PSS 189,749,248 B. Browser/helper processes are recorded separately. These are
asynchronous samples with shared mappings, not an exact additive decomposition;
the gap to the service peak and remaining process/cache/allocator costs are not
assigned to a specific R object. The scientific tracing is passive, with no added
GC or continuous allocation/CPU profiler.

Decision: stop scaling at the failed small lifecycle gate. The one correction
pass was used to correct trace installation. A narrowly scoped production
metadata-validation repair and one small GUI recheck require renewed authority;
the question has been put to the user. Do not skip plot/session contents, bypass
resource validation, repeat the unchanged pilot or start the 1/2 GB members.
Large attempts used: zero; large inputs generated: zero. R/ and NAMESPACE still
match frozen WxHeXq, and no service remains running. Even after repair, the
18 GiB diagnostic requires >=20 GiB host-available memory at launch.

Evidence paths above are relative to the existing proof root
/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/browser-smoke.<member>/.
Pilot measurement bindings: shell
1f6c3008d6e6d0b6031d14b2fdcc72eff688c3b4247dd8438c65a77162c26e6a;
R runner 97bde21f91a1145ab3b2c791652031b3e1538b82b04ed9b7358ca055f5d18309;
driver c1e980a455a71bcb82ea62c98b221a1bb74ca22e105e8c0a5b09de827b0a483c;
diagnostic helper 9d5a77d6ea6e05b2709a28bec4e74d90d5c6e9d3043d83c75ddd7c981ac1e377;
recipe 2b8fd057297871eaf06549ddf321794c62388ecf83e0b3e44b9ead1b0aef0e2d.
The subsequent reproduction changes only the diagnostic test (SHA-256
38c8c8b3e759adc7513d067dbc67db34f0d71f6bb1e698605d2ca65643ae3ef2),
not those bindings or installed runtime.

## Retained revision 1.22 — verifier continuation and rejected candidate

**Bounded trial finished; candidate not accepted and runtime changes withdrawn.**
The fresh current control HqtCut/WxHeXq completes and verifies. Candidate
Lr0SpH/IqtAsc hits the unchanged early-memory guard before legacy serialization,
DA and Summary complete. The >=10% whole-workflow memory gate fails. Stop after
two of four allowed full-size members; do not run the conditional second pair or
repeat an unsuccessful member. All services are stopped. R/ and NAMESPACE again
match the previously verified WxHeXq source exactly; earlier repairs and user
changes are preserved. IqtAsc source/runtime, all fixtures and failed evidence
remain available. Verifier-only corrections remain; they change no application
calculation. Wider DIA-025/performance/adoption gates remain open.

| Matched full GUI member | Result | Combined peak | Complete elapsed |
| --- | --- | --- | --- |
| HqtCut / WxHeXq | Complete; independent 9p6i0J passes | 7,365,005,312 B (6.86 GiB) | 842.570410 s |
| Lr0SpH / IqtAsc | Guard exit 76; incomplete | 8,163,643,392 B observed (7.60 GiB) | Not available; not a speed denominator |

Both have zero recorded job swap. Candidate is a controlled early stop, not a
machine OOM. Its last delivery records correlation_filtered and all 16 committed
states, zero working payloads, and 32 hydrates before export. The export request
directory contains no serialized session/journal, and no public latest/completion
file exists. The last sampled span is 618.442 s, not completed elapsed time.
No full-output verifier is run against this incomplete candidate.

The target was first exceeded at 1790073994.828580, between qvalue_filtered
delivery and protein_peptide_filtered delivery, before either modified heavy
operation (IQ/export). Within that interval, maximum sampled current is
6,703,603,712 B for control versus 7,312,908,288 B for candidate. Associated
anonymous-memory samples are 5,026,619,392 / 5,160,042,496 B and file-memory
samples 1,638,903,808 / 2,121,666,560 B. These are asynchronous cgroup snapshots,
not exact additive decompositions or proof that the refactor caused the extra
memory. At the final pre-stop sample, file memory is 2,654,773,248 B. The tested
export/IQ lifetime changes do not establish control over the earlier whole-job
peak or file-cache pressure. Do not subtract file cache or weaken the target to
declare a pass.

Outcome: isolated export sharing reduces the observed peak only modestly and
costs additional CPU; isolated IQ R-process HWM is essentially unchanged. Small
science/PNG/session checks pass, but complete large performance acceptance fails.
Further optimization requires a separately bounded decision about earlier QC
and I/O/ownership costs, not another unmeasured exporter variant. No larger
dataset or other-omic claim follows from this trial.

### Completed continuation and preserved intermediate evidence

User approves completing verification-only wiring within DIA-025. IqtAsc and all
scientific/performance settings remain frozen. Correct only checked relocation of
the two IQ files for the replay layout, and remove the unsupported reader-route
16 GiB flag (full-GUI verify-da's existing 16 GiB allowance is unchanged). Reuse
completed outputs; do not repeat a scientific or performance member.

IQ control **zCcQIY / WxHeXq** completes at 29.741180 s service elapsed, 15.085 s
Apply, peak 3,900,661,760 B and zero swap. **q23wvg** independently confirms
byte-identical peptide_values_imputed.tsv and iq_output_file.txt, then rejects
only config$peptide_input_file and config$iq_output_file: the replay uses outputs/
instead of the original GUI QC directories. This verifier failure is retained
(peak 3,854,069,760 B, zero swap), not treated as a scientific disagreement or a
passed verification. Explicit owned-path and checksum checks must precede the
two-field mapping. Candidate IQ, small GUI pair and all four conditional full-size
GUI members are still unused. No new performance acceptance is claimed.

Verifier-only **JQunyO** passes 26 assertions, zero failures/errors/warnings,
three explicit guarded-member skips; peak 1,020,772,352 B, zero swap. It covers
the prior comparison checks, copied-store ownership and rejecting an IQ config
path that differs from the specifically verified artifact. Runtime R/ still
matches installed IqtAsc exactly. Recheck saved zCcQIY next; no IQ repeat.

Only verifier code/selection and removal of the unsupported reader 16 GiB option
changed; measurement bodies and installed libraries did not. Updated bindings:
shell fbbf0a70bfe92c9abdf19ad2198c76a99dfdcaff9ecfa06eecfd19d1d4529bec;
driver 0328738be392cf89543d8614fb0ccc9920cd8d372e909c903a0351d71d3380a3;
mechanism fixture 784b9c28d253c3827a3b4557fbde5c401269ad00631e1148e57d050d4acb8ce8.
The R runner hash and both runtime hashes remain unchanged. Earlier valid export
verification is reusable: this correction is conditional on IQ replay only.

Small resident **gFWuGu / IqtAsc** completes the original full journey with the
l6OHWn report/FASTA/design, uninstrumented and sleep-inhibited: 251.310102 s
service-tree elapsed, peak 6,516,121,600 B, zero swap. **gJ4dFT** independently
passes complete legacy/DA/Summary checks (peak 1,050,775,552 B, zero swap).
The saved-state small member on the identical build/input is next. These small
members are correctness checks, not the large memory/time acceptance comparison.

Small saved-state **tCCZKy / IqtAsc** completes: 238.779630 s service-tree elapsed,
peak 2,425,114,624 B, zero swap. **nTVluk** independently passes complete legacy,
DA/Summary and the matched small-pair comparison: all 15 scientific states at
1e-12, exact types/attributes/missingness, and identical scientific PNG contents.
Verifier peak is 1,088,442,368 B, zero swap. Final candidate ownership is 16 states,
one full-history read, 59 hydrates and zero native/working/cache payloads/open
resources. Both small members have verified installed load/exit identity.

The full-size queue is now released on unchanged representative/report.tsv and
tJ6U1c FASTA/design: WxHeXq, IqtAsc, then IqtAsc, WxHeXq only if the first pair
passes correctness and peak <=90% / elapsed <=120%. Maximum four GUI members,
fresh installed processes, serial, profiling disabled, existing sleep inhibitor,
8 GiB hard / 7.5 GiB early / 900 s / zero swap and unchanged host safeguards.
Startup, complete automatic outputs, DA/Summary and teardown remain included.
All runtime and latest verifier/driver hashes above are frozen. The first current
member completes as **HqtCut / WxHeXq**: 842.570410 s service-tree elapsed,
peak 7,365,005,312 B (6.86 GiB), zero swap. All 16 states and GUI actions complete,
installed load/exit identity verifies, and delivery reports one full-history read,
59 hydrates and zero native/working/cache payloads or open resources. The sampler
has no gap over two seconds; the sleep inhibitor was registered. Independent
legacy/DA/Summary verification **9p6i0J passes** using the existing 16 GiB verify-da
allowance, separate from the application's unchanged 8 GiB limit: 16 historical
states, 8,139 proteins, 48 samples, peak 5,190,971,392 B, zero swap.

First-pair candidate target is <=6,628,504,780.8 B peak. The <=120% time ratio
would allow 1,011.084492 s, but the unchanged 900-second service ceiling remains
binding. This is a newly matched current control, not a comparison against an
old interrupted/OOM member. IqtAsc candidate member 2/4 is Lr0SpH; its retained
guard failure and withdrawal decision are recorded at the top of this document.
No second pair runs and no remaining attempt is silently repurposed.

Saved control IQ recheck **gAM1ov passes** (peak 3,623,555,072 B, zero swap),
without repeating zCcQIY. Candidate **To7T3v / IqtAsc** completes: 29.605827 s
service, 14.778 s Apply, peak 3,236,192,256 B, zero swap. Independent **KDyQ9D**
passes exact IQ input/output TSV hashes and the protein S4/provenance comparison
at 1e-12 after the two checked path relocations (peak 3,625,254,912 B, zero swap).

Do not attribute the combined IQ peak difference entirely to R ownership:
application HWM is essentially unchanged (control 1,914,892,288 B; candidate
1,913,987,072 B). File-cache/I/O attribution differs substantially (control device
read bytes 1,145,520,128; candidate 65,536). Both members finish with two hydrates
including the head, zero working/cache payloads and zero open resources. These
isolated observations do not establish the >=10% whole-GUI target.

Both mechanism pairs and their scientific verification are finished. Start the
approved small IqtAsc resident/saved GUI pair using unchanged l6OHWn inputs, with
the existing sleep inhibitor and 8 GiB/7.5 GiB early/900 s/zero-swap guards. Full
GUI timing includes startup, scientific work, automatic files, DA/Summary and
teardown; the still-unused large comparison remains gated on this small pair.

## Retained revision 1.21 — comparator correction and export observations

User requests continuation after the revision 1.20 stop. Reinstate its two runtime
files from frozen aAlj3D and repair only attribute-name/value traversal; retain the
attribute-order fix and unchanged IQ computation. A focused installed recheck
must pass before the remaining approved mechanism/GUI comparisons. WxHeXq remains
the frozen control. No extra performance attempt or limit change is introduced.
Installed correction **IqtAsc** passes **tLXxfQ**: 1,856 assertions, zero failures,
errors or warnings, six explicit historical/guarded-member skips. All 116 new
sharing/IQ/copy-lifecycle assertions pass. The existing ten small cycles, export
recovery and scientific/module checks in that selector pass. Peak is 1,168,818,176 B,
zero swap. This is correctness evidence, not a performance improvement.

Frozen candidate RDB: 36ef62e517ac8f06e879af7a3f0047134adcf52e6adc77eaddf7f7ec235ea0d0.
Source R/ matches its installed snapshot. The renewed correction is used; failed
gvgew7/lv04yQ remain unchanged. The pre-existing unrelated density-dispatch failure
is still retained, not repaired or represented as passing.

Mechanism/GUI bindings: shell da549c185e9c2c4a54f16718097d7cc51f12d100219875db0f8dec92a7957b89;
R runner 62020a11ec782bfb72bb7057ffb9581087470cf391f24d07c423bf71c4ee02b3;
driver 8055e62f51f164c74463450a13742a2e44cf27cff00bc9d776b004ec596b3512;
mechanism fixture 6d118f064c6b83c09614e4f2ba72d975a1c4629c89d72ddc0d35b809a5907815.
No further runtime/harness changes are planned during comparisons. Begin the
one current/candidate export pair on IVL1gM's complete copied store, then IQ pair;
startup/copy/reopen/preparation/publication are measured, independent readback is
separate. Guards stay 8 GiB hard/7.5 GiB early/900 s/zero swap. No performance
GUI comparison member has started.

Current export replay **8jA7BF / WxHeXq** completes: peak 5,830,955,008 B, zero
swap, 24 whole-table reuses, 18 hydrates (head/metadata plus 16 history records),
one serialization and no scientific checkpoint. Final working/native/cache
payloads and open store resources are zero. Two complete 609,339,926-byte legacy
files have matching SHA-256 003b6828a10fdfec0e1678b8f8f379e1364081cc76d79e7ececca538f6c470b3.
Independent scientific verification subsequently passes as aCPn6X below.

**Timing is infrastructure-invalid, not a speed denominator:** system suspend
occurred 19:14:08–19:15:24 AEST and the sampler has a 75.6355-second gap. Service
active monotonic interval is 115.523319 s, whereas operation wall elapsed is
179.636 s; neither is presented as an uninterrupted end-to-end timing. The
mechanism reader route does not acquire the existing browser-only sleep inhibitor.
Retain the observation/limitation and do not repeat this member or infer an export
speedup. Remaining mechanism observations are diagnostic only; the unchanged full
GUI route already acquires the inhibitor and remains the performance decision.

Independent verifier launch **gSHrdH** is retained as an instrumentation failure:
the shell requested the allowed 16 GiB verifier ceiling, but the reader entry
point requires 8 GiB and rejected it before loading the package or data (peak
65,003,520 B, exit 1). No export or scientific computation was repeated. Resume
the same verification as **aCPn6X** using that route's existing stricter 8 GiB
cap, unchanged data/comparisons and no code change. This is not a performance
member or an application resource increase. **aCPn6X passes** complete legacy
restore and all 16 historical states at 1e-12, with exact table/column attributes,
types and missingness. Peak 7,562,170,368 B, zero swap under the 8 GiB verifier
limit. No control repeat is needed.

Candidate export replay **rMcRSb / IqtAsc** completes: 166.446410 s active service
interval, 151.120 s operation, peak 5,532,545,024 B, zero swap. Application HWM
is 3,639,173,120 B versus the control's 3,776,487,424 B. It reuses 48 tables in
total, including 25 support-table instances, and 57 columns. Both members decode
18 payloads including head/metadata, with one complete-history collection and
zero working/native/cache payloads or open resources after publication. Candidate
produces two complete 609,339,927-byte legacy files with matching SHA-256
ef3f1f8bb651672acd43ed6aa3f236d89da72ff8d01381496361400235af04c9.

The observed combined-peak difference is only 298,409,984 B (about 0.28 GiB);
it is not whole-GUI acceptance. Candidate service CPU is 162.169374 s versus
113.953470 s for the interrupted control; record the extra work, not a speedup.
Control suspension prevents a clean paired elapsed comparison. Independent
candidate verification **2zPaDd passes** complete legacy restore and all 16
historical states under the same 8 GiB route; peak 7,932,071,936 B, zero swap.
The verifier uses throwing checks and a bound passed_history_science receipt;
testthat's "empty test" label reflects the absence of expect_* calls in that
standalone verifier body, not skipped execution of its scientific comparisons.
Both exports are now independently verified. Proceed to the unused IQ pair.

## Retained revision 1.20 — export sharing and IQ lifetimes

**Stopped at correctness; combined candidate withdrawn.** The approved correction
allowance is exhausted. No large mechanism or full-GUI comparison ran; there is
no new memory or throughput gain. R/ and the three profiling entrypoints match
their pre-pass verified WxHeXq/current hashes again. Earlier repairs stay intact.
The frozen control is WxHeXq, with the same 556,028,391-byte report/tJ6U1c companions.
Its source R/ exactly matches the pre-change working tree. Runtime RDB SHA-256:
153db129f2911ff008dbeda2a77b8d229c0816c726d2d9e368e71b7be9e32044.

Two mechanisms only: export-local immutable column/support-table sharing, and
IQ TSV-preparation/post-commit object lifetimes. Scientific calculations and full
legacy files stay fixed. Small saved-session inspection finds exact support-table
duplicates across 8–11 states; this is feasibility evidence, not measured large
memory relief. IVL1gM's earlier sampled peak is 7,153,307,648 B during the interval
between imputed and protein_s4_created, versus 7,894,335,488 B overall. It does not
identify IQ arithmetic as the cause; export alone has little room to deliver the
agreed 10% whole-run reduction if that earlier peak remains unchanged.

The exact ordered checks, four-member maximum, one correction allowance and
<=90% peak/<=120% elapsed paired decision are recorded in DIA-025 revision 1.20.
Original application guards remain; independent large verification is separate.
No old OOM control repeat, unbounded search or general scaling claim is authorized.

Initial installed candidate cbMUYn/check gvgew7 retains two failures: tibble
column replacement changed attribute pairlist order (values/attributes were not
dropped), and the previously recorded MCI-007.6 density S4 dispatch error. The
single correction restores exact attributes including compact automatic row names.
The unrelated density implementation remains unchanged; its known failing module
artifact test is retained, not repaired or counted as passed. The recheck selector
keeps the affected IQ module contracts and new native/saved ownership fixtures,
but omits that broader unrelated artifact file. No large members have run.

| Candidate / check | Correctness result | Guarded check peak / swap |
| --- | --- | --- |
| cbMUYn / gvgew7 | 1,909 passed assertions; one failed attribute-order assertion; one known unrelated density-dispatch error; zero warnings; three skips | 1,231,495,168 B / 0 |
| aAlj3D / lv04yQ | 1,646 passed assertions; six export errors; zero failed assertions/warnings; five skips | 1,166,553,088 B / 0 |

The corrected comparator recursively passes `attributes(value)` back into itself.
A named list's attributes are another named list, eventually reproducing the same
`names` structure indefinitely. This causes node-stack overflow; existing export
guards then reject the incomplete history. This is an introduced implementation
defect, not an OOM, scientific-method problem or valid negative performance result.
The earlier tibble failure was attribute pairlist ordering, not dropped values.
The independent IQ release/TSV checks and tiny copied-store lifecycle checks pass,
but cannot justify retaining the combined candidate without the required gates.

Failed runtimes and complete R sources remain in browser-smoke.cbMUYn and
browser-smoke.aAlj3D; check artifacts remain in gvgew7 and lv04yQ. Their RDB hashes:

- cbMUYn: 528106edad9d4c92f68e5011645ab1bd857355662f93039a82668e22fff5fde1.
- aAlj3D: 2e7d0298bbbb801ecff1262c6a27a00a24ce58158a46e108ecb8c85acd2c9b57.

Only this pass's two runtime changes and unfinished harness routing/comparison
changes were withdrawn. Tests/fixture source remains, capability-gated for the
withdrawn candidate functions; skipped trials are not counted as passing. No
additional install, test recheck or performance run follows the spent correction.
The remaining targeted correction is direct comparison of attribute names plus
recursion through individual attribute values, retaining compact row names and
attribute order. It requires a renewed bounded correction allowance. Application
capacity remains the previously verified WxHeXq/IVL1gM/TJ8UMW result, not this trial.

## Retained revision 1.19 — targeted export-owner release

**Targeted repair verified.** Installed WxHeXq completes the frozen 556 MB full
GUI journey under the unchanged 8 GiB/7.5 GiB early/900-second/zero-swap limits.
Independent large verification passes. No further run or repair is needed for
this bounded item; all services are stopped. This is a same-input capacity
result, not a qualified 10% speed/memory gain or a scaling claim. Broader DIA-025
performance/adoption gates remain open.

User authorizes fixing the observed serialization/copy memory overlap and running
the harness again. The immediate reference is zDUwrx; failed 2CYM6M and rejected
GC trial JZBPMb/bek8Wq remain immutable evidence. Only the existing export/state
ownership boundary changes: collect/serialize in a child call, release history
and its working-map roots, then copy/fingerprint/publish using metadata only.
One post-release collection is part of this explicit boundary, not a revival of
the rejected pre-hydration collection or a blanket GC policy.

The existing v1 journal records durable serialization before copying; its latest
checksum is NULL until copying and independent fingerprinting succeed. Retries
must use that staged source without hydrating history or serializing again.
Complete legacy session RDS, output inventory, science and public APIs stay fixed.
Small ownership/failure fixtures and affected installed checks precede a small
matched GUI pair/comparison, then one uninstrumented large candidate and verifier.
One bounded defect correction/recheck; no unchanged failed large repeat. GUI runs
retain 8 GiB/7.5 GiB early/900 seconds/zero swap and existing host guards. Large
independent verification has the separate 16 GiB/900-second ceiling. No larger
recipe, worker/backend change or new diagnostic series.

Implementation is installed as **WxHeXq**; source R/ matches its snapshot and
NAMESPACE matches reference zDUwrx. RDB SHA-256:
153db129f2911ff008dbeda2a77b8d229c0816c726d2d9e368e71b7be9e32044.
Affected **6wH4xC** passes 1,872 assertions, zero failures/errors/warnings and
three explicit historical rejected-experiment skips; peak 2,364,252,160 B, zero
swap. All 46 new assertions pass: weak-reference proof for both history and
working-map roots before copying, partial-copy output-only retry, staged-source
corruption rejection, serialized-checkpoint interruption, and serialization/
incomplete-history failure safety. Existing round trips and ten small cycles pass.

Small resident **Qjb1su** completes on WxHeXq/l6OHWn under 8 GiB/900-second limits:
150.992473 s, 6,249,705,472 B peak, zero swap. Independent **ykUtkA** passes legacy,
DA and required-output checks (peak 1,049,595,904 B, zero swap). Small candidate
**DZ2AwQ** completes at 147.001649 s, 2,359,570,432 B peak, zero swap. Independent
**IB48Pj** passes all historical science at 1e-12, exact types/attributes/missingness,
identical PNG contents and full legacy/DA/required-output checks (peak 1,080,565,760 B,
zero swap). These small observations are not formal performance acceptance.

Uninstrumented large **IVL1gM** completes WxHeXq on the unchanged report/FASTA/design:
**443.926401 s** service-tree elapsed, **7,894,335,488 B (7.35 GiB)** peak, zero
job swap, within 8 GiB hard / 7.5 GiB early / 900 seconds and existing host guards.
It includes startup, all included QC/IQ/normalization/correlation actions, complete
session export, DA views, Summary and teardown. Input and shell/R-runner/driver
hashes match the preceding frozen bindings. It has only about 0.15 GiB headroom
below the early guard; do not extrapolate to the unrun 1–2 GB recipes.

The publication journal is complete: one serialization, zero added scientific
checkpoints, two 609,340,738-byte complete legacy files with matching SHA-256
f9b9ff5a0e1cb563a3408321f402d0e0a193abf1a30cd5abe819b6eb2ab32ace.
Delivery reports 16 states, one full-history read, 59 hydrates and zero native/
working/cache payloads or open store resources. Compared with the previous
guarded failure, the new run reaches publication and every downstream consumer;
the aborted predecessor is not a finite timing or percentage-gain denominator.

Independent **TJ8UMW** passes full legacy restore, DA and required Summary/output
checks: 16 history entries, 8,139 proteins and 48 samples. Peak 5,220,052,992 B,
zero swap under the separately declared 16 GiB/900-second verifier ceiling.
This continuation used one large GUI attempt and no correction. Earlier failures
remain below. The only collection added is after the serialized-history frame
and its working roots have expired; no diagnostic/pre-hydration collection runs.

## Retained revision 1.18 — bounded repair completed; capacity failed

Revision 1.18's four-run allowance is exhausted. The lifetime release is
implemented and correctness-checked, but **no useful full-size performance gain
is demonstrated and the 8 GiB workflow target fails**. The pre-export GC trial
worsened memory and is withdrawn. All services are stopped. DIA-025 and broader
performance/adoption gates remain held; no further run or adjacent refactor is
authorized by this continuation.

### Frozen executions and outcome

All large members use the same 556,028,391-byte report and original tJ6U1c
FASTA/design. Scientific methods, settings, codec and output inventory are
unchanged. Passive diagnostics disable diagnostic GC and retain bounded scalar
boundaries plus process/cgroup sampling; no continuous profiler or cache flushing.

| Member / installed build | Boundary | Peak bytes (GiB) | Service elapsed | Decision |
| --- | --- | --- | --- | --- |
| NwNqwZ / U5J0lV | Passive baseline, 16 GiB | 8,289,714,176 (7.72) | 422.722159 s | Complete diagnostic, not acceptance |
| zldOEd / zDUwrx | Passive lifetime-only, 16 GiB | 8,242,167,808 (7.68) | 851.856699 s | Insufficient memory benefit; substantially slower |
| bek8Wq / JZBPMb | Passive + one application pre-export collection, 16 GiB | 8,619,323,392 (8.03) | 544.444968 s | Memory regression; GC withdrawn |
| 2CYM6M / zDUwrx | Uninstrumented, 8 GiB | 8,186,499,072 (7.62) | Incomplete; guard exit 76 | Final capacity failure |

Every member has zero recorded job swap. Diagnostics had 15.5 GiB early stops
and 3600-second ceilings. Final 2CYM6M retained 7.5 GiB early stop, 900 seconds,
zero swap and unchanged host safeguards. Sampling can overshoot the early
threshold before the stop; 7.62 GiB is the observed guard peak, not a changed cap.
Do not calculate a completion speedup from the stopped member.

Application HWM is 5,692,030,976 B in NwNqwZ, 5,738,033,152 B in zldOEd and
5,876,330,496 B in bek8Wq. The tiny lifetime-only cgroup difference is not a
demonstrated application-memory improvement. Near-peak anonymous/file samples
are approximately 6.007/1.463 GiB in baseline versus 6.058/1.574 GiB in
lifetime-only; asynchronous accounting is not an exact additive decomposition.
Elapsed increases also occur outside the changed paths: design post-checkpoint
2.729 -> 6.247 s, Q-value Apply 9.561 -> 20.317 s, protein evidence
8.562 -> 18.618 s, hydration 17.083 -> 38.337 s, RDS staging 41.778 -> 96.105 s.
The cause is not isolated; neither attribute all of it to the code change nor
dismiss it as proven host noise.

### Final failure is localized to the latest-session staging copy

2CYM6M committed all 16 states through correlation_filtered. Its last completed
flush reports zero working/native/cache payloads and open resources. That flush
precedes export; its full_reads=0 does not mean export never materialized history.

The retained request directory is
`browser-smoke.2CYM6M/dia_baseline_smoke/scripts/proteomics/.workflow-export/dia-export-28900b46587a1a/`.
Its session.rds is 609,343,884 B, while latest.rds is only 576,978,944 B.
There is no publication journal/completion manifest in this request. Together
with the existing serialize/hash/copy order, this localizes the stop to copying
the completed serialization into the staged latest file, before publication.
Do not treat either the incomplete latest file or the initial running smoke
receipt as completed output. All pending inputs/files are preserved.

The last cgroup sample records 8,039,649,280 B current, 5,915,635,712 B anonymous
and 1,864,888,320 B file memory; the guard then records 8,186,499,072 B.
The sampled span before termination is 545.842183 s, not elapsed completion.
This was a controlled early stop, not an observed machine OOM. Final DA/Summary
were not reached.

Code inspection identifies the next bounded ownership target:
stageWorkflowLegacyExport retains its session_data argument during hashing and
latest-file copying, while its callers also retain the full history. A compatible
next proposal is to finish serialization, release full-history owners while
retaining metadata/current data needed by remaining outputs, then copy/fingerprint/
publish. Removing one inner alias or collecting before hydration does not achieve
that separation. This is **not implemented or authorized as an extra repair**;
it needs a scoped change and verification without dropping files/history or
weakening recovery. A new storage format is not established as necessary.

### Retained implementation and correctness evidence

Retained runtime changes are only releaseDIAWorkingReference in the existing state
adapter and its two post-commit Q-value/protein-evidence call sites. The exact old
binding and local input are released before downstream preparation; active output,
durable history and files are preserved. Native managers are unchanged. No added
production GC remains. R/ matches installed zDUwrx byte-for-byte.

- plEmIb: 11 passive instrumentation assertions pass; no failures/errors/warnings/
  skips. Peak 2,032,898,048 B, zero swap. The long-wait transport proof was reused.
- tL0NnS: 1,826 assertions pass, no failures/errors/warnings, two historical audit
  experiment skips. All 91 lifetime assertions pass: ownership, active/foreign
  references, nested scope, reapply, audit on/off, post-commit ordering, failed
  payload/head writes, revert/reset and no extra current-output decode.
  Peak 2,250,248,192 B, zero swap.
- Small matched zDUwrx/l6OHWn GUI pair: Qa8ZoM resident 280.848823 s /
  6,236,020,736 B; DChluN candidate 258.386212 s / 2,345,369,600 B; both zero swap.
  3VgJDk verifies resident outputs (peak 1,049,792,512 B). 4J1gBV verifies the
  complete pair (peak 1,079,934,976 B): all historical scientific states at 1e-12,
  exact types/attributes/missingness, identical PNGs, full legacy/DA/Summary outputs.
  Existing execution-metadata normalization is unchanged. This is not a G4 pass.
- 9ChEZl: JZBPMb trial checks pass 1,841 assertions, no failures/errors/warnings,
  two historical skips; all 15 collection-order/retry assertions pass. Peak
  2,180,616,192 B, zero swap. That fixture remains capability-gated to the rejected
  trial build; passing correctness did not qualify its memory regression.
- gzgzmt independently verifies the **completed zldOEd diagnostic**, not failed
  2CYM6M: full legacy restore, required outputs, DA and Summary pass; 16 history
  entries, 8,139 proteins, 48 samples. Peak 6,314,049,536 B, zero swap under the
  separately declared 16 GiB/900-second verification ceiling. This cannot convert
  the final guard failure into a capacity pass.

Keep every failed small check: GWoEY1 omitted browser propagation in the new
selector (peak 1,139,675,136 B); 3s2qhk had seven passing assertions and one
unsupported wait_for_value(value=) error (peak 2,040,643,584 B); kYOfHj stopped
on an extra closing parenthesis in the new export fixture (peak 2,295,660,544 B).
All had zero swap. Corrections were fixture/runner-only and retained under the
user's explicit continuation, not hidden retries.

### Immutable bindings and disposition

Evidence root is the repository sibling
`MultiScholaR-dia-proof-20260911-VhC1n5/`; each named member has its own directory.
Retain receipts, guard stops, phase/process/cgroup samples, failed logs and
recoverable checkpoints. Full-size attempts: **4 of 4**. No larger data or old
resident OOM control was rerun.

RDB SHA-256:
- U5J0lV: 8347b064e340c8d79c1d10e6abc0dc8d9bafa363ac67fba75846bff14085429b
- zDUwrx (retained): 5d7a25a5f8d11534f0a018c173107a44ff7cc0b6bfbff3e8e563cc17baa038f1
- JZBPMb (rejected GC): 8a84dd18a03ae36fa0069deb45d0d11ca97ee3d619a63d9b7c3c3695b9bf8b4b

Passive driver: 138a91c1b37d8baea963219f1ea54dae0b9f836291a6994a04ddf5b53510f9a8.
R runner: 62020a11ec782bfb72bb7057ffb9581087470cf391f24d07c423bf71c4ee02b3.
Shell: 33ff31b3d931e16a345a113cdb6f2de4fc05edd64c0d5d71d669ac70058abf97.
All input hashes match the earlier frozen bindings. Science/correctness,
diagnostic observations, final capacity and formal performance acceptance remain
separate. The bounded plan ends at its measured-limit stopping condition.

## Retained revision 1.16 — raised-cap export-memory diagnosis

Uninstrumented **4Aqm6u** stopped at the 7.5 GiB early guard, exit 76:
current 8,095,657,984 B; peak 8,095,748,096 B (**7.54 GiB**); zero job swap.
The journal records 392.160 s wall time; no normal completion receipt exists.
This is a controlled guard failure, not an observed OOM or a performance pass.
The stale initial smoke receipt says running; the service is stopped.

All 16 checkpoint states through correlation_filtered were committed. The
request-owned export directory is empty: no session RDS, latest RDS, journal or
export metadata/correlation descriptor was produced. The last completed flush
has zero native/working/cache payloads; it predates export preparation. The final
sample reports approximately 5.02 GiB anonymous and 2.13 GiB file memory; sampled
accounting is not an exact additive decomposition. Export preparation/history
hydration is implicated, but this uninstrumented run cannot identify its exact
internal allocation. DA and Summary were not completed; no large verifier ran.

User now authorizes **one** unchanged U5J0lV/556 MB full-path memory diagnostic:
16 GiB hard, 15.5 GiB early stop, 3600 seconds, zero swap, existing host guards.
Extend only the harness trace boundaries for history hydration, immutable-table
sharing, serialization and publication; verify them first on small fixtures.
Application RSS/HWM, per-process attribution and cgroup anonymous/file counters
remain separate. Existing completed-action GC is diagnostic-only; no continuous
allocation profiler, production GC or scientific change. No acceptance-target
relaxation, larger input, old OOM-control rerun or unchanged large retry.
Instrumentation check **Wpf2ky**: 59 assertions passed, two failed, zero errors/
warnings/skips. The new test invoked pre-trace copies from the test environment,
so installed entry/exit events were absent for its direct calls. Correct only
those test calls to resolve from the instrumented namespace. The one allowed
correction/recheck **mgqM0y** passes all 61 assertions with zero failures/errors/
warnings/skips (peak 4,473,888,768 B, zero swap). Preserve Wpf2ky. The helper itself and
all runtime code are unchanged by this test correction.

The authorized full-path diagnostic **Pt8flT** completed on the fixed installed
candidate/input under the declared 16 GiB/60-minute diagnostic limits. All proof
jobs are stopped. The small check verifies actual hydration, serialization,
publication and injected serialization-failure boundaries, plus existing browser
survival, observer cleanup and capped-log checks. No further correction/retry is
authorized in this continuation.

Frozen diagnostic helper SHA-256:
3cd83a3058cce0919338e32da5488efd1868c5ae6675a08428e939e8ec4ad188.
U5J0lV RDB, input hashes, shell/R runner and driver match the previous release.

### Pt8flT outcome and memory attribution

The installed GUI route completed import, all included peptide/protein QC,
cyclicloess normalization, RUV Skip, correlation preview/Skip, legacy export,
owned-session DA loading, DA tables/static and interactive volcano/heatmap,
Summary publication and teardown. Service-tree interval **397.727409 s**;
peak **7,050,833,920 B (6.5666 GiB)**; zero job swap and zero recorded OOM events.
This is diagnostic execution completion, not scientific/performance acceptance.
The receipt is `passed_summary_ui_only`; no new independent large-session
oracle was run. Existing small-pair scientific verification remains separate.

Application PID 2268559; driver PID 2268242. Only these two R PIDs appear in the
process samples. No `startDIAWorkerCluster` or peptide-audit emitter trace occurs.
The prior replicated R-worker memory surge is absent on this route/build.

| Observed operation | Application RSS high-water / sampled phase maximum | Finding |
| --- | --- | --- |
| Q-value Apply | 3.633 GiB phase maximum | Large early R allocations remain with provenance audit disabled. |
| Protein/peptide-support Apply | 4.375 GiB phase maximum | Largest early peptide-QC application peak; post-GC live heap 1.232 GiB versus 3.921 GiB RSS. |
| Protein replicate Apply | 4.108 GiB phase maximum | No step worker pool; no old multi-process surge. |
| Correlation plot preparation | 4.107 GiB first-call phase maximum | No new application high-water mark on this fixed preview/Skip route; does not cover correlation Apply/reapply. |
| Complete-history reconstruction | HWM rises from 4.374 to 5.173 GiB | Largest application high-water increase is within saved-history hydration, before RDS serialization. |
| RDS staging/publication | Approximately 4.993 GiB sampled application RSS | No further application HWM increase; whole-job peak occurs around completion/post-export collection, with file-cache memory present. |
| Owned-session DA load / model fit | 4.119 / 4.077 GiB phase maximum | Does not rematerialize complete history; no new high-water mark. |

Phase maxima use application RSS boundary records and external PID samples;
`VmHWM` separately captures a transient peak between samples. They are not
incremental allocation sizes, live-object sizes, or additive phase costs.

**Named export path:** `collectProtNormExportSessionData()` calls
`materializeDIAStates()` (`R/utils_workflow_dia_state.R:97`), which resolves each
record through `resolveDIAStateValue()` and `hydrateSavedWorkflowEntity()`
(`R/utils_workflow_saved_entities.R:1081`). During this one export, working
records grow from 1 to 16. The high-water increases first appear on hydration
exit for raw_data_s4 (**4.682 GiB**), qvalue_filtered (**4.846 GiB**), then
protein_peptide_filtered (**5.173 GiB**). This boundary includes decoding and
decoded-object validation; it does not isolate a specific native qs2 allocation.
The three compressed payloads are each about 84 MB; their saved S4 size estimates
are approximately 476, 614 and 645 MB respectively, not unique live-heap sizes.

Exactly 24 immutable tables are shared, but sharing occurs **after** each
complete record is hydrated. It therefore does not prevent that hydration's
transient allocation. History-materialization entry/exit RSS is 4.106/4.927 GiB;
the transient HWM is higher. The journal reports one RDS serialization, two
609,344,805-byte legacy files with matching recorded SHA-256, no added scientific
checkpoint. Final manager ownership has zero native/working records, zero payload
cache and zero open store resources; complete-history read count is one.

**What the remaining memory is—and is not:** near the whole-job maximum, the
application has approximately 4.955 GiB PSS, driver 0.175 GiB PSS, and sampled
browser processes 0.371 GiB PSS. The cgroup sample reports 5.317 GiB anonymous,
1.203 GiB file and 0.034 GiB kernel memory. These views overlap and are sampled
at slightly different times; do not add them into an exact peak decomposition.
Application high-water is 5,554,733,056 B. The main application—not an R worker
pool or multi-gigabyte profiler log—is the dominant process memory owner.

At completed export, diagnostic GC changes application RSS from 5,360,558,080 to
4,345,507,840 B; live R heap after collection is 935,398,056 B (0.871 GiB).
The remaining approximately 3.18 GiB RSS/live-heap difference is **unattributed**:
R free heap capacity, allocator retention, native allocations and mappings are
not separated by this tracker. It is not evidence that all of it is live S4
history, nor proof that R-process memory is harmless. Do not prescribe blanket
production GC from this observation alone.

**Why this does not erase 4Aqm6u:** that run's last sample was 7.491 GiB current,
5.017 GiB anonymous and 2.134 GiB file memory before its 7.54 GiB guard receipt.
Pt8flT's peak sample has *more*, not less, anonymous memory (5.317 GiB), but lower
file charge (1.203 GiB). File-cache conditions/reclamation and diagnostic GC both
differ; their causal contributions are not isolated. No new runtime fix occurred.
This is evidence locating export hydration and earlier peptide-QC memory costs,
not a reproduced uninstrumented failure or an accepted 8 GiB solution.

Evidence root: `MultiScholaR-dia-proof-20260911-VhC1n5/browser-smoke.Pt8flT/`
(sibling of the repository). Retain application-phases.jsonl, process-samples.tsv,
memory-samples.tsv, smoke-receipt.json, terminal-interval.json and terminal-properties.txt.
Application phase log is 2,660,175 B; combined phase/process/memory/console logs
are under 5 MB. No Rprof/Rprofmem; no truncated trace coverage, process-sampling
gaps or memory-sample errors; maximum cgroup sample gap 0.232 seconds. Sampling
still cannot attribute every native allocation or exact simultaneous peak.

Bounded diagnostic is finished; no further jobs or runtime changes are launched.
Remaining decision: isolate the uninstrumented export failure with collection
disabled and controlled/recorded file-cache conditions before claiming a repair;
then constrain any change to the measured hydration/ownership mechanism. That is
not authorization for a new backend, scientific change, further trial or blanket
optimization. DIA-025 and the original performance gates remain open.

## Retained revision 1.15 verification continuation

User approved the verifier-only alias correction after kukNbG. Add only the
legacy removeProteinsWithOnlyOneReplicateFn/core_utilisation location already
covered semantically by the scheduling allowance, with regression checks that
grouping parameters and scientific values still fail comparison. Recheck the
saved pair after the focused verifier test; no runtime rebuild or GUI repeats.
The 1e-12 tolerance, immutable input/build bindings, output inventory and resource
limits remain unchanged. The observed small elapsed regression remains evidence.

Alias check **SEgZIG** passes 13 assertions with zero failures/errors/warnings/skips,
including rejection of changed grouping parameters, scientific values, attributes
and missingness. Re-verifier **BP7hfo** passes the saved WExzW7/tshxAi comparison:
all historical scientific states at 1e-12, exact types/missingness, identical project
PNGs, complete legacy restore and downstream output checks. Neither GUI was rerun.
Original kukNbG remains retained. This clears correctness/output release only;
the small elapsed regression still prevents an overall performance-pass claim.

At release, final uninstrumented candidate **4Aqm6u** was launched on unchanged U5J0lV and
the original 556 MB report/FASTA/design. Limits: 8 GiB hard, 7.5 GiB early stop,
900 seconds, zero swap and existing host guards. Startup, automatic outputs,
full session export, owned DA load, all included views, Summary and teardown
remain inside completion. No Rprof/Rprofmem, diagnostic trace or diagnostic GC.
The updated verifier driver SHA-256 is
c757b3459b63c529ad33af0117bb9b98927f97c5f294c4ff48c5ce58571d7dac;
runtime, inputs, shell/R-runner bindings and scientific tolerances remain unchanged.
If this member completes, its independent verifier also retains 8 GiB/900 seconds.

The final fixture-only recheck **2voFmA** passes **2,330 assertions**, zero
failures/errors/warnings/skips; peak 1,141,780,480 B, zero swap. The protein UI
fixtures now supply their own non-live worker token and assert requested counts;
real worker lifecycle/science evidence remains separate. Runtime U5J0lV is unchanged.

GUI release: one matched small resident/candidate pair using existing l6OHWn
inputs, each followed by the existing independent verifier; then the original
556 MB candidate only if small science/output comparison passes. Ordinary members
and verifiers use 8 GiB/900 seconds, 7.5 GiB early stop, zero swap and existing
host guards. No profiler or diagnostic GC; temporary sleep inhibition stays
job-owned. No repeated component trials, old large resident control or larger data.
Both small GUI members have completed; comparison is held at the verifier as
described below. No large candidate has been launched.

Small resident **WExzW7**: 203.745303 s service interval, peak 6,235,004,928 B,
zero swap. Independent verifier **mf7g8L** passes complete outputs/legacy restore.
Small candidate **tshxAi**: 238.150978 s, peak 2,374,557,696 B, zero swap. Both use
the frozen l6OHWn input hashes and installed U5J0lV, without profiling/diagnostic GC.
This single pair has substantially lower candidate memory but **higher elapsed
time** (approximately 16.9%); it is not a performance-gate pass. Candidate delivery
has 16 states, zero working/native payloads/cache bytes/database leases and one
complete-history materialization, with 47 indexed descriptors.

Candidate verifier **kukNbG** reaches the historical comparison after its output/
DA checks, then fails at protein_replicate_filtered on only the reported worker
backend metadata difference: `args$removeProteinsWithOnlyOneReplicateFn$core_utilisation`
is logical NA versus DIAWorkerSpec. The existing comparison explicitly permits
this execution-only difference under the canonical method name, but omits the
legacy Fn alias produced by updateParamInObject/getFunctionNameSecondLevel.
It has not completed comparison of all subsequent states/PNGs, so no complete
pairwise scientific pass is claimed. Receipt retained; verifier peak
1,078,956,032 B, zero swap. All services stopped.

A verifier-only alias correction, negative regression checks for real scientific
parameter differences, and re-verification of the existing saved pair have been
requested. No runtime change, tolerance change, GUI repeat or larger run is
authorized by this observation. Reuse the valid GUI executions if that correction
is approved; keep the original failed verifier unchanged. A documentation patch-
context mismatch made no file change before this corrected record.

Execution bindings: U5J0lV RDB
8347b064e340c8d79c1d10e6abc0dc8d9bafa363ac67fba75846bff14085429b;
shell e23f04e1180d957c46ea5d4771cf0aacf89e5766a7e74df37b56ea502de09495;
R runner 92659d2d3e55c4356fe25b7946e1cc9ddef21cf45ead55c3815c037a41ece56b;
driver e0f6278159dc71765f3e57ba9946c68bcef8253657c5afb1f19f7d1445b7830a.
Input fingerprints are captured inside the service and checked against frozen
preparation receipts. Limits and formal G4 thresholds are unchanged.

Renewed user continuation follows **uu1Qc3**: 2,301 passed assertions, ten failed
assertions/two errors confined to two older protein UI fixtures, zero warnings or
skips; peak 1,170,247,680 B, zero swap. Import mocking and mock lifetime are fixed.
Removing the leaked new_cluster mock exposed the fixtures' accidental real-worker
creation around a fake S4 without args. Their worker factory is now explicitly
mocked and requested counts asserted; real lifecycle/component checks are not
replaced. One affected recheck is authorized, followed by the existing GUI queue
only if clean. No runtime source change/rebuild or additional diagnostic.
The new local-worker test passes, as do owned restore/retry/sharing checks and
all 110 cycle assertions (including the ten explicitly expected plot warnings).

User "okay well - let's queue it up!" authorizes the requested one fixture-only
correction/recheck and then the already planned GUI sequence if clean. The
interruption occurred before edits/jobs; no run was duplicated. U5J0lV remains
the immutable runtime candidate. Four test files change: import reader mock,
fake protein-ID slot, test-owned mock cleanup and explicit expected cycle warning.
No runtime science or limits change. Reuse completed component results and retain
the previous NLnJla failure and spent host-suspend replacement. Current stage:
corrected affected batch, before GUI release. No extra correction is implied.

## Active continuation — full-path memory repair, 2026-09-22

The user approved the diagnosis/repair plan and changed GUI scheduling authority:
local protein counting and deferred worker creation are allowed; standalone APIs
remain compatible. Audit-off baseline is oax7SK. Revision 1.12 of DIA-025 records
the fixed workload, scope, sequence and one total correction allowance. Attribution
and runtime implementation are now complete; the candidate is not released for
full-GUI verification because affected fixture checks failed. All jobs are stopped.

### Current handoff — runtime changes present, affected batch not clean

Installed **U5J0lV** matches the runtime source. Affected batch **NLnJla** completed
with 2,296 passed assertions, seven failed assertions, three errors and ten plot
warnings; peak 1,326,436,352 B, zero swap. The new deferred-config and unused-pool
checks pass. Export hydration/sharing, serialize-once/owned-generation reuse,
publication retry, external restoration and all 100 ten-cycle assertions pass.
The cycle renders retain the existing histogram's clipped-row warning ten times;
no numeric or plot-limit change was made to hide it.

Remaining failures are localized, not waived: the import server fixture still
overrides readConfigFile instead of the new readProtGuiConfig boundary. The older
protein-replicate fake S4 lacks protein_id_column, which the unchanged count code
requires. Its local_mocked_bindings calls use the package namespace as cleanup
lifetime instead of the owning test frame, leaving the intentional `mock replicate
failure` installed and contaminating the later new local-count test. These findings
are supported by the test source, failure trace and testthat's installed cleanup
implementation; no failed case was silently removed or rerun.

The original correction allowance was used for the independently documented host
suspend. An additional **fixture-only correction and affected-batch recheck** was
requested asynchronously; it has not been assumed granted. The small GUI pair,
final uninstrumented 8 GiB/900-second member and its verifier remain unrun.
Do not claim ticket closure, full-workflow savings or production promotion.

### Fresh-process mechanism comparisons

Each pair used identical input/scientific hashes, separate installed processes,
the original 8 GiB/900-second/zero-swap guard, and no profiling or prewarming.
Controls use oax7SK; candidates use U5J0lV. The failed fixture batch does not
contaminate these independently started R processes.

| Mechanism | Control / candidate | Whole-service peak bytes | Operation plus readback seconds | Service seconds |
| --- | --- | ---: | ---: | ---: |
| Protein count | JIjGkX / 9DEERz | 4,487,421,952 -> 988,413,952 | 20.913 -> 2.731 | 36.183497 -> 18.625679 |
| Full-history export | F9q5BU / Ft6eDd | 1,352,896,512 -> 1,334,988,800 | 10.440 -> 10.419 | 29.485458 -> 28.915745 |

Worker main-R vector high-water is unchanged at 137.6 MiB; the removed processes
account for the clear service-memory benefit. This is a small protein-count result,
not a whole-GUI percentage. Its input SHA-256 is
ab975799c44b551d539eb8f1967adaa48e971219c035594b5b275e30ecdf21a9;
scientific table/design SHA-256 is
36511dad20732a646db993c72c866c08f8094ae49e66242bbcfc94fd897860b9.

Export vector high-water falls 255.6 -> 239.2 MiB, serializations 2 -> 1, and
37 identical immutable tables are shared. Its elapsed result is effectively flat,
not a meaningful speedup claim; service peak improves only about 1.3%. Input SHA-256
4430bbcf258501f1e7d0e14ec01a01708915790ee39519ebafa753ba8646227c
and all 14 historical scientific-state hashes match. Required copies/fingerprints
and readback remain measured. Both mechanisms deliver zero working/native payload
records and zero database resources. All four members completed with zero swap.

Next bounded work, only with the requested allowance: fix the reader mock and fake
S4/cleanup lifetime in the two affected fixture files, explicitly assert the known
cycle-render warning, rerun the affected batch, then resume the already planned
GUI verification if clean. Runtime science and resource limits must not change.

### Findings before implementation

| ID / severity | Evidence and impact | Smallest repair / closure evidence | Blocked layer |
| --- | --- | --- | --- |
| MEM-01 / high | tHy0Cw: four active protein workers approximately 4.1 GiB plus eight existing workers approximately 1.05 GiB. readConfigFile starts a pool eagerly; several peptide helpers explicitly ignore it. | Defer creation, use existing local count branch, preserve exact result/order and prove owned process cleanup. | Runtime acceptance |
| MEM-02 / high | Main R remains approximately 6.1 GiB after peptide QC despite empty state working/cache owners. Live-object versus allocator residency is unresolved. | Diagnostic-only collection and bounded owner evidence; repair only reproduced obsolete references/preparation. | Diagnosis |
| MEM-03 / high | collectProtNormExportSessionData materializes complete history; saveProtNormExportArtifacts serializes it twice. | Existing export-scoped decode ownership and serialize-once verified-copy publication; full legacy round-trip, one serialization, failure recovery. | Implementation |
| MEM-04 / high | DA load reads the complete legacy session and republishes all history even for the unchanged session just exported. | Reuse only a verified same-project/committed-generation export; changed/external files keep full validation and restore. | Implementation |
| MEM-05 / medium | Current lightweight trace cannot distinguish live heap from RSS or identify worker namespace/payload costs. Earlier continuous profiling produced a 7.56 GB log. | Bounded completed-action heap/owner observations and small worker probes; prove detection, release, error cleanup and recording limits. | Evidence |

Findings above remain unchanged; append closure evidence separately. No full GUI
gain is inferred from a helper result or from the historical OOM control.

### Diagnostic progress (not repair acceptance)

Replacement **qBXdTj** completed the full audit-off route on unchanged oax7SK:
13,238,194,176 B peak (12.329 GiB), zero swap, 925.588646 s service interval.
This includes completed-action diagnostic GC and is not timing acceptance.
All input hashes match the original report/FASTA/design. Sixteen states, zero
audits, zero working/native payloads and zero cache bytes/leases at final delivery.
The temporary sleep inhibitor was released; no job processes remain.

Post-collection live heap (bytes) versus RSS: after design 632,037,744 versus
2,482,642,944; Q-value 1,277,648,984 versus 3,959,001,088; protein evidence
1,330,678,184 versus 4,549,337,088; precursor rollup 842,091,696 versus
4,365,361,152; imputation 963,715,056 versus 4,480,901,120; IQ 822,552,264
versus 4,562,612,224. This does not show cumulative retention of every earlier
assay. RSS includes allocator residency and cannot be called a live-history leak.
Export RSS reaches 5,972,799,488 then falls to 4,504,211,456 after collection
(943,095,392 live bytes); restore reaches 6,399,602,688 then 4,506,357,760
(943,325,992 live bytes). No blanket production GC or speculative owner cleanup
is implemented from those observations.

Candidate **U5J0lV** implements deferred unused pools, local experimental protein
counting, one full-session serialization, output-only retry and same-generation
DA reuse with a streamed selected-file fingerprint. Export materialization also
shares exactly identical immutable assay/design/identifier tables and replaces
duplicate scoped-cache records, preserving attributes and copy-on-modify isolation.
Mutable tables and differing values/attributes are never merged. No component
schema or backend is added. This targets the measured complete-history rebuild,
not a guessed general leak. Candidate RDB SHA-256:
8347b064e340c8d79c1d10e6abc0dc8d9bafa363ac67fba75846bff14085429b.
Runtime source matches the installed candidate; affected checks are in progress.

Fixed isolated members: worker uses the existing small protein fixture; export
uses the existing representative **protein-only** fixture (10,000 proteins,
48 samples, 14 states), with unchanged tables between each fifth small value
change. No new large precursor input is generated. Run one fresh installed
control/candidate process per mechanism, serially at 8 GiB/900 seconds/zero swap.
No prewarming or profiling; preparation/startup are included in service metrics,
and the operation/readback timer is reported separately. These are component
observations, not full-GUI performance acceptance.

qughH6 completed: 50 passed assertions, zero failures/errors/warnings/skips,
peak 4,466,380,800 B, zero swap. Full audit-off attribution **LLiVXb** then failed
in design setup: browser transport EOF/no response, peak 5,326,340,096 B,
zero swap and zero cgroup OOM events. It is incomplete infrastructure evidence,
not a capacity failure. The retained host journal shows suspend requested at
08:01:32 AEST, suspend entry at 08:01:36 and resume at 08:05:43; the controller
has a matching 247.439-second sampling gap. Import-only collection reduced main
RSS 2,327,191,552 -> 1,373,900,800 B, with 627,372,752 B live heap after GC;
later QC ownership remains unmeasured.

Use the single correction/replacement allowance for this documented suspend
interruption. The existing browser runner gets a process-owned, verified
sleep/idle inhibitor, released on ordinary or failed exit. No persistent host
power setting changes. Test that lifecycle separately, reuse the passing
diagnostic batch, then replace only the invalid large diagnostic on unchanged
installed oax7SK. No additional correction allowance remains afterward.
Two documentation patch-context mismatches made no file changes before this
corrected record; neither was an empirical run or a runtime correction.

Small attribution batch **qughH6**, installed oax7SK: the completed-heap fixture
detects a deliberately retained 16 MB numeric allocation and its release. The
protein-worker probe uses a 13,072-byte table. Both local and four-worker results
have identical protein tables and design. Each worker reports MultiScholaR loaded,
194 namespaces and approximately 877,000 KiB RSS; serialized dispatched expressions
peak at 524,251 B. This identifies heavyweight package loading for a tiny grouped
count, not simply large assay shards. The batch still includes its existing
real-browser blocking-action lifecycle check before large attribution is released.
No runtime scheduling change or full-path audit-off claim follows yet.

Attribution helper binding:
a5e6d1bbf14ae83cc68a2a7558b30630e09cccf3cf3b6d6a6473ba98f6379c26;
driver fc45331f3ffd8e29c14e7a6447231f02067e3437ddd0101d2bc16ada1bb4487e;
R runner 30edf467564dfdf785b1649250c7620259004a00bc97511c6b8db157db6ef889;
shell bdf6f92cff835640ce6e3f4b691f5405f2bf3c1f248b13fce1193325ec8167bb.
Only the explicit diagnostic uses completed-action GC. It reports pre-collection
RSS and post-collection live-cell bytes/RSS, not a fabricated pre-GC live counter.
All owner descriptions are shallow class/address/shape values, never retained
producer frames or scientific payloads. No full-journey allocation trace is added.

## Current decision — audit-enabled journey verified; GUI audit switch verified

The unchanged 556,028,391-byte report (3,840,000 rows, 10,000 input proteins,
48 samples) completed the installed **1YRrHo** candidate route in **tHy0Cw**.
The service-tree interval was 1616.403990 seconds (main process 1616.143990),
whole-job peak **14,234,578,944 B / 13.257 GiB**, and job swap zero. Startup,
workers, automatic outputs, session export/reload, DA views, Summary and teardown
are inside that service interval. It is not outer-controller resource accounting.
The 16 GiB/3600-second diagnostic succeeded, not the original 8 GiB/900-second
performance gate. No large resident control was repeated or used as a denominator.

Independent verifier **Nz1IpR** passed `passed_capacity_output_checks` under the
existing 16 GiB/900-second guard, peaking at 6,330,843,136 B with zero swap. It
restored the complete legacy session (16 history entries, 8,139 final proteins,
48 samples), checked installed DA results and heatmap/output/summary consumers.
The IQ/QC inventory subcheck is explicitly presence-only; complete historical
scientific equivalence comes from the matched small pair, not an unrun large
resident comparison. All three large input hashes match the frozen records.

Coverage: all seven peptide QC actions, real IQ, protein QC, cyclicloess,
RUV Skip, correlation Skip (including pre/post Pearson plots), full session
export/reload, DA fit, complete table, static and interactive volcano, heatmap
save and Summary actions. The interactive view contained 8,139 points and 48
expression points and passed selection behavior. Other normalization/RUV choices,
correlation Apply/reapply, larger recipes and consolidated acceptance remain held.

### Memory-first attribution of the completed audit-enabled run

These are sampled concurrent service currents within named helper intervals,
not independently additive allocation costs. Main-R values are RSS at the helper
boundaries, not live-object heap measurements.

| Interval | Sampled service maximum, GiB | Main-R RSS entry -> exit, GiB |
| --- | ---: | ---: |
| Q-value Apply | 11.540 | 2.797 -> 8.381 |
| Protein-evidence Apply | 10.873 | 7.348 -> 8.764 |
| Precursor rollup | 9.716 | 7.813 -> 6.013 |
| Peptide intensity | 8.490 | 6.023 -> 6.052 |
| Peptide sample | 8.497 | 6.052 -> 6.052 |
| Peptide replicate | 8.567 | 6.052 -> 6.055 |
| Peptide imputation | 8.620 | 6.055 -> 6.082 |
| IQ | 8.754 | 6.082 -> 6.112 |
| Protein replicate | 13.253 | 6.105 -> 6.105 |
| Normalization | 8.387 | 6.106 -> 6.092 |
| Complete session export | 9.471 | 6.070 -> 6.367 |
| DA views and Summary | 9.524 | 6.904 at final Summary exit |

- **Application provenance audit remains the major early allocation path.**
  At the Q-value maximum, main-R PSS was about 8.64 GiB, the original eight
  workers about 1.05 GiB and the driver about 0.17 GiB; nearest cgroup file
  charges were about 1.43 GiB. The whole-table canonicalization/digest was
  executing then. Main-R cumulative RSS high-water reached 9,787,936,768 B
  (9.116 GiB) during the subsequent protein-evidence data digest. This is
  application work, not the harness memory logger. Early audit peak attribution
  must not be confused with the later whole-job maximum.
- **Active protein-replicate workers now set the whole-job maximum.** At the
  nearest process snapshot to sampled maximum (1789989247.978579), main-R PSS
  was 6.061 GiB, four step workers 3.880 GiB, original eight workers about
  1.047 GiB and driver 0.175 GiB. The four step workers rose to about 4.095 GiB
  PSS just before exit, but were last observed at 1789989253.898896 and did not
  survive the helper return. Service current fell to 8,873,451,520 B by
  1789989255.036981. The ownership repair works; active worker replication is
  still a cost, not permission to change worker counts or science.
- **Late consumers still have material residency.** Export increases main-R
  RSS from 6.070 to 6.367 GiB; restore reaches 6.885 GiB. The largest sampled
  DA renderer PSS is 689,034,240 B (0.642 GiB), versus the earlier incomplete
  run's approximately 1.73 GiB observation. These are different complete/incomplete
  runs, not a qualified percentage memory gain. No claim that all GUI roots or
  every live allocation have been identified follows from these observations.

The adapter delivered 16 states with zero working/native payload records,
zero native audits, zero payload-cache bytes and zero open database resources;
metadata was 378,592 B. That establishes those owners are cold, not that R has
returned all memory to the OS. Live GUI objects, allocator-retained heap and
temporary overlap are not individually resolved by RSS. PSS, RSS, file charges
and cgroup accounting overlap and must not be summed into an exact breakdown.

Boundary/console/driver logs total 2,657,052 B. No Rprof/Rprofmem was enabled.
There are 11,288 cgroup samples (largest gap 0.393031 s) and 1,129 process
snapshots (largest gap 1.5258 s). The one retained memory-sample error is the
terminal cgroup disappearing after successful exit; process-gap log is empty.
The exact service peak is slightly above the sampled current maximum; no
sub-sample function attribution or precise live-heap decomposition is claimed.

### Separately authorized GUI audit-off change

The user requested disabling the audit while the installed harness was running.
Only the new proteomics GUI state-manager constructor changes to
`WorkflowState$new(audit_enabled = FALSE)`; the cold adapter inherits that flag.
No hot patch or restart altered tHy0Cw. This skips the additional audit emitter,
canonical hashes, audit summaries and audit-ledger construction. Scientific
methods' own counts/exclusion records, ordinary logging, checkpoint persistence
and storage integrity checks remain. Standalone WorkflowState defaults and
explicit `audit_enabled = TRUE` still work; imported audit history is not erased.

This provenance mechanism originated in janitor commit **ccce54f**, 2026-08-04,
under **PQC-008**, before the performance branches. Its ticket explicitly permits
disabling additive audit capture without changing scientific results. The costly
whole-table canonical sort/serialization is not required to calculate the filter.

Fresh installed audit-off build **oax7SK** has RDB SHA-256
`867d724b3217fbcebf08fddbe2e3420861890c6bc18ee8cf7d0011713f78629e`.
Affected check **K0kvsX** passes **265 assertions**, zero failures/errors/warnings,
with two explicitly skipped withdrawn audit experiments. It verifies the GUI
flag, cold-adapter propagation, no emitter/hash call on the disabled save path,
exact scientific S4 preservation, checkpoint/revert/RDS restore and explicit
audit opt-in. Other affected state/module regressions pass. Build and checks ran
serially after the large verifier under 8 GiB/900 seconds/zero swap; no correction
was needed. Runtime source matches oax7SK. No audit-off full journey or performance
measurement has been run, so none of the audit-enabled timings/peaks quantify
the switch's benefit. DIA-025 and broader adoption/performance gates remain open.

## Retained execution protocol — post-repair full harness authorization

The user has released the full-GUI pause. Run the matched small resident/candidate
pair and independent verifiers on installed **1YRrHo**, then one unchanged 556 MB
candidate diagnostic if they pass. The older small control belongs to Q7BMJh and
cannot satisfy the comparer's matched-package requirement. No large OOM control
or larger recipe will be repeated. The exact allowance and stopping rules are in
DIA-025 revision 1.10. No job survived the interrupted turn; runtime/build matching
and host headroom checks pass.

Small post-repair pair passed: resident **lm4iE5** completed in 249.284228 s with
6,260,658,176 B peak; independent verifier **WgZxcU** passed. Candidate **JL02bz**
completed in 239.852992 s with 6,362,570,752 B peak; verifier **Y6aqMD** passed,
including all historical scientific states at the existing tolerance, full legacy
restore and identical project PNGs against lm4iE5. Both members had zero job swap
and matched the l6OHWn input hashes. These small results release the large
diagnostic, not performance acceptance. The candidate is slightly faster but has
a slightly higher small-case peak; no universal improvement is inferred.

At launch the unchanged 556 MB diagnostic was next. If complete, its independent
large-session verifier will use the existing 16 GiB/900-second verifier option,
zero swap and the same early/host safety stops; the 60-minute allowance applies
to the GUI diagnostic, not silently to the verifier.

Frozen run bindings: RDB
`129ac7b8cf44b19c202859287283fe72f17e4cb8b151a4ce42e131ac01f56c18`;
driver `4f6b08eeb54d65c17b39396a3f0291ad93dcfbb5520887f1602a3de17276d454`;
R runner `12e6317cfc79c8422dad0433a748f950df1421c84948a922d311a69d3d656f92`;
shell `98b573020178a3df497b499bb2ae0eed1a89e391fc15a65a42fbbfb99ef6a5e7`;
lightweight helper `ff0271deb881cdc049190a07044f9dde6a40b6b0fe38e0a4853051e655730048`.
Input fingerprints are captured inside the existing service boundary and checked
against the retained preparation/frozen input records. Small runs retain 8 GiB/
900 s; large diagnostic retains 16 GiB/3600 s, zero job swap and the existing
15.5 GiB/host-headroom safety stops. No continuous allocation profiler is enabled.

## Retained repair handoff — pre-harness repairs implemented; performance held

User requested the repair plan end to end, then a pause before rerunning the GUI
harness. The implementation and focused verification of steps 1–3 and 5–6 are
complete. The six retained changes cover worker lifetime, unused Glimma preparation,
image ownership, retired-session references, browser table rows and audit-key
encoding. Step 4's full-GUI run remains explicitly paused. The current installed
candidate is **1YRrHo**. No small or large full-GUI journey was rerun.

The audit-key change now has matched fresh-process evidence of lower memory and
identical input/scientific/audit hashes. **Elapsed time did not improve in that
pair**. This is a memory-repair candidate, not throughput or full-workflow
performance acceptance; the original 10% journey gates are not passed or weakened.

Current verification: **ttScvL** passes 55 new exact-encoding assertions, 19
key/canonical/hash assertions and the affected science/state checks. It also
retains the known predecessor density-dispatch error and a stale image-reference
test assumption. **uATYVF** corrects only that reference assumption and passes all
82 focused assertions, with zero failures/errors/warnings/skips. Runtime code is
unchanged between those checks. All services are stopped. DIA-025 stays held for
the paused full-GUI/performance work and its broader unfinished acceptance checks.

Retained earlier affected check **LNcnvQ** has 1,971 passed assertions, zero failed
assertions, one known predecessor density-dispatch error, one existing tidyselect
warning and three skips (unreadable cp08 plus two explicitly held audit-experiment
checks). All retained repair-specific checks pass. This is not a clean broad-suite
pass, ticket closure or scientific/performance promotion. All services are stopped.

### Retained implementation

| Item | Change and evidence | Remaining limitation |
| --- | --- | --- |
| Protein-replicate workers | The GUI action stops its own pool on exit, including calculation/output failure. Results save the existing portable worker specification. The S4 method can recreate that exact pool for restored execution and stops only pools it creates. Seventeen focused assertions pass, including live execution, terminated owned workers, surviving borrowed workers, restored S4 calculation and reversion. | The four-worker 4.10 GiB observation belongs to the prior large run; whole-journey savings have not been remeasured. Native/shared callers retain their scientific behavior; ordinary history storage stays unchanged. |
| Glimma preparation | Removed discarded annotation/display preparation and released the intermediate counts data frame. The complete widget payload, dependencies and prepended content compare identically with Q7BMJh. On the 2,048-row/12-sample fixture, starting vector heap is 151.9 MiB for both; vector-heap high-water is 163.3 -> 158.3 MiB; recorded allocations >=10 kB are 4,921,768 -> 4,440,056 B. | Isolated preparation evidence, not browser or combined-process peak. Single instrumented elapsed observations 0.295 -> 0.112 s are not performance acceptance. |
| Composite image ownership | Reused the existing plot-owner stripping helper. The fixture proves the former plot environment retained a 1,573,088-byte decoded image; the new plot keeps the raster but not that frame. PNG pixels are identical and plot RDS remains printable. | No extrapolated large-GUI memory percentage. |
| Retired session references | Closing an adapter clears its obsolete commit callback; the replacement initializer releases its previous-manager variable. Closed callback removal and active replacement callback/manager behavior pass the existing state tests. | This removes concrete reference links; retained GiB were not quantified. |
| Glimma browser rows | The existing DataTables deferred-row option is set on the Glimma table only, using an htmlwidgets render hook. LNcnvQ observes one table with all 1,024 records and 53 columns, but only 10 initial row nodes versus 1,024 in the native control. Full-row search, selection/clear, identical widget data and unchanged global defaults pass. | This proves avoided initial table-row construction, not a percentage reduction in whole-browser or workflow RSS. Total DOM counters include prior navigated documents and must not be used as an additive memory comparison. |
| Audit-key encoding | One helper avoids ifelse's full-length branch copies for ordinary atomic columns and skips replacement when no values are missing. Classed, nested and matrix columns keep the original path. Exact byte/encoding, empty/all-missing/nonfinite, factor/date, matrix/list and input-preservation tests pass; full Apply result hashes match. | Fresh-process peaks are lower, but this pair's elapsed times are higher. No digest/schema/sort/population change, row chunking, new cache or explicit GC is retained. Full-GUI time/memory gates remain open. |

### Final audit qualification and implementation — 2026-09-21

The user's renewed completion request authorized a diagnostic, one repair
comparison, its bounded correction and affected verification. The frozen
55,296-row/12-sample input and existing 8 GiB/900-second/zero-swap guard were kept.

**Measurement correction, not erased evidence.** QGxAHT ran unchanged audit code
twice. Both executed one emitter, five digests and two substantive digests, and
produced the same full result. Yet vector-heap high-water rose 241.5 -> 290.3 MiB.
The emitter body SHA-256 was identical. Thus the earlier sequential observations
are real observations, but attributing their differences to the candidate code
was not qualified: R process history alone can produce a substantial difference.
The earlier candidates remain withdrawn and all failed evidence remains intact.
The fixture's confidence q-values are all below 0.01, so it also was not a strong
performance probe of large removed-row ledgers for the duplicate-digest proposal.

The retained bounded allocation logs in QGxAHT show 142,432,000 B in stacks
containing `ifelse`, out of 388,897,368 B recorded at the >=10 kB threshold.
Canonicalization stacks account for 165,485,712 B; these inclusive totals overlap
and must not be added. The selected fix replaces only the avoidable branch-vector
construction, preserving byte-length prefixes, separators, missing-value tokens,
global ordering and digest serialization. It is used by both canonicalization and
identity-key preparation. The old path remains for special column types.

**Matched members.** Control **9JGTF2** uses ZnV7XD; initial candidate **7Neko5**
uses RaRw1Q. The initial candidate reduced memory but Apply elapsed was 5.189 s
versus 4.866 s. The one correction skips even empty subassignment on NA-free
ordinary columns via `anyNA`; final candidate **oj1ozA** uses **1YRrHo**. The
unchanged control was reused, not rerun. Each member has its own fresh installed
R process; no profiler or prewarming is enabled. The existing service includes
startup, fixture preparation, Apply, checkpoint publication, downstream readback,
result fingerprinting and teardown. The Apply sub-timer is reported separately.

| Observation | Control 9JGTF2 | Final candidate oj1ozA |
| --- | ---: | ---: |
| Vector-heap high-water | 255.0 MiB | 221.6 MiB (13.1% lower) |
| Main-R RSS high-water | 1,213,612,032 B | 1,146,298,368 B (5.5% lower) |
| Whole-service cgroup peak | 1,171,193,856 B | 1,130,356,736 B (3.5% lower) |
| Apply helper + downstream read elapsed | 4.866 s | 5.219 s (7.3% higher) |
| Full service main lifetime | 20.672937 s | 22.623661 s (9.4% higher) |
| Job swap | 0 B | 0 B |

These memory counters are different accounting views, not additive quantities.
This bounded helper pair is neither the full GUI journey nor a 1–2 GB scaling
claim. The elapsed tradeoff is retained explicitly; it is not called a speedup
or excused by the memory improvement.

Both members have input SHA-256
`06fe18ee1ae235b9bba3cd5f48311f974fb21030b0fa51ba6fa2a2afb507cf42`
and complete S4/science/audit SHA-256
`ff583ff5cdef6baa28d05c59a7909beb10206390575c74eaf30ffe4abc3e1147`.
Both have exactly the initial/raw/qvalue states, one audit record and no remaining
working records or database leases. The fresh member test source hash is
`9a3558019a8d00bf83b60695ff765e41758ede0f16f6f03715a89171cca3523a`.

**Verification closeout.** ttScvL has 2,025 passed assertions, one failed assertion,
two errors, one warning and three skips. Besides the known density error, its
image comparison incorrectly expected the ZnV7XD reference to retain the old
decoded image even though that build already contains the image repair. The test
now explicitly distinguishes an already-fixed reference; candidate owner, pixel
and plot-RDS checks are unchanged. uATYVF reruns that affected file and the existing
narrow verification checks: 82 passed, no failures/errors/warnings/skips. All new
encoding, scientific-hash, worker, session and browser checks pass using these
matching receipts. The broad run's baseline error/warning/held-prototype skips are
not relabelled as a clean suite pass. Peaks: ttScvL 3,243,016,192 B; uATYVF
3,100,061,696 B; both zero swap.

Final installed RDB SHA-256 (1YRrHo):
`129ac7b8cf44b19c202859287283fe72f17e4cb8b151a4ce42e131ac01f56c18`.
Final audit runtime source SHA-256:
`703fbadc7c09739db555657b6c25d815b527123d55095ce9719922bb01fb79f2`.
All six pre-harness repairs are implemented and focused verification is complete.
No further optimization or run is started here. The next authorized-work boundary
is the user's explicit pause before the full GUI harness; full performance
acceptance and the unrelated existing density-dispatch issue remain open.

### Retained preceding completion continuation — 2026-09-21

The user's instruction to continue through all repair items released the remaining
bounded browser diagnosis, audit comparison/correction and final affected checks.
No service limits, scientific thresholds or input populations were relaxed.

Browser diagnosis **OGVhog** captured the actual render exception: the synthetic
annotation supplied `gene`, which Glimma also adds from the count row names. Its
duplicate tooltip key caused a Vega expression parse error before table creation.
This was a fixture defect, not a memory stop. Renaming that synthetic annotation
to `Protein_Ids` preserves its values, width and 1,024-row/12-sample population.
The first rendering comparison **ZT1zRW** completed but its diagnostic directly
invoked renderValue while the normal callback was pending, producing two tables.
That run is retained but is not the accepted one-render comparison.

The final fixture requests a diagnostic direct render only after a readiness
failure, requires one table, and uses the **installed application helper**, not a
global DataTables default override. `deferProtGlimmaTableRows` adds a per-table
data attribute between Glimma's synchronous table creation and its asynchronous
jQuery-ready initialization. DataTables reads that attribute normally. No new
dependency, backend, prototype patch, model change, data truncation or lookup
adapter is introduced. LNcnvQ required no diagnostic direct render. Its historical
receipt still has the old `application_changed: false` label from the earlier
read-only trial; the test now identifies the installed application helper instead.
The observations themselves are retained unchanged.

Two later audit candidates were also rejected:

| Candidate / check | Exact result | Decision |
| --- | --- | --- |
| Reuse the duplicate substantive digest; W8ePkB / **4xihXe** | Same complete S4/science/audit SHA-256. On 55,296 rows and 12 samples, recorded allocations >=10 kB rise 386,376,792 -> 390,153,768 B; vector-heap high-water rises 237.4 -> 244.8 MiB; elapsed falls 3.469 -> 2.831 s. | Faster in this observation, but not an accepted memory repair. |
| Same reuse plus a single large first-import collection; sYeEyp / **N9nEHj** | Same complete S4/science/audit SHA-256. Allocations remain 386,376,792 -> 390,153,768 B; high-water is 237.4 -> 244.9 MiB; elapsed 3.669 -> 3.529 s. | Correction did not solve the observed memory regression; withdrawn too. |

The comparison includes the production Q-value Apply helper, audit, durable
checkpoint publication and downstream state readback; it excludes GUI plots and
startup. Before/after result SHA-256 in both trials is
`ff583ff5cdef6baa28d05c59a7909beb10206390575c74eaf30ffe4abc3e1147`.
These are valid negative observations at that boundary, not proof that no audit
optimization can work. No failed result was silently rerun with a changed workload.
Both experimental implementations remain in their installed build/source artifacts.
The current audit runtime is restored exactly to
`0ccc3e8a3f7686d53f1aca72d96d16e91e2c46d1d3a7b67d3dc172d0fcccf1ae`.
The two new audit-experiment tests require their explicit guarded selector; they
are held, not reported as passing production features. A further audit change
requires another evidence-directed scope, not automatic acceptance of a faster
but higher-memory result.

Build/check records (same evidence root as below):

- **ZnV7XD** is the final installed repair build. Its source contains the five
  retained repairs and the original audit code. Build peak 1,124,802,560 B, swap 0.
  RDB SHA-256: `43c582ab47356e0a47c1d8c4d59b6e3013c26b40593a737a57997548a5d2863d`.
- **LNcnvQ** final affected checks: peak 3,285,143,552 B, swap 0; all retained
  repair checks pass. One unchanged density dispatch failure remains, as confirmed
  in fresh predecessor/candidate processes. The check service reports failure;
  its receipt has not been relabelled successful.
- Final Glimma preparation comparison in LNcnvQ preserves the full scientific
  widget payload: vector-heap high-water 162.2 -> 157.2 MiB; allocations >=10 kB
  4,921,768 -> 4,440,056 B. These remain isolated observations, not workflow gates.
- Browser native/deferred row-node counts are 1,024/10, both with 1,024 records,
  53 columns and global deferred rendering false. Search, selection and clear
  pass. Browser-wide DOM totals span 1 versus 4 documents after navigation; no
  total-DOM or combined-peak saving is inferred from them.

All planned attempts and the correction for this continuation are used. No
full-GUI rerun occurred. The audit item is still unresolved, so the repair plan
and DIA-025 remain held rather than being marked complete merely because the
other repairs passed. The next decision is whether to authorize a further bounded
audit investigation or evaluate the five verified repairs in the paused harness;
neither is started implicitly.

### Retained rejected audit experiment

The 20,000-row temporary-key chunking experiment preserved exact keys,
canonicalized tables and SHA-256 digests across 40,003-row edge fixtures, but
failed the isolated memory decision. With the same 151.8 MiB starting vector
heap, high-water increased **207.1 -> 211.2 MiB**. Recorded allocations >=10 kB
increased **51,535,880 -> 55,865,416 B**. The elapsed observations were
0.345/0.337 s; this does not rescue the memory regression. It is neither a
complete-Apply measurement nor evidence that all chunking approaches must fail.

The runtime audit file was restored byte-for-byte to its pre-repair contents:
SHA-256 `0ccc3e8a3f7686d53f1aca72d96d16e91e2c46d1d3a7b67d3dc172d0fcccf1ae`.
No threshold was moved and no larger/easier recipe was substituted. The failed
implementation remains recoverable in installed trial lb0zYc's source and library.
The audit-memory bottleneck remains unresolved; another mechanism requires a
bounded follow-up, not an unrecorded third attempt.

### Verification, failed evidence and exact pause

All artifacts below are under
`/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5`.

- Installed repair trial **lb0zYc** built successfully under the ordinary guard;
  peak 1,172,803,584 B, zero swap. Its predecessor is Q7BMJh. RDB SHA-256
  `ad35878f44800a5d66c67f88bb42c54800adeabe897dbac2ab4b97b1d787659d`.
- Initial focused check **KtrIRg**: 1,935 passed assertions, two failed assertions,
  three errors, one warning and one skipped unreadable cp08 snapshot; peak
  3,594,637,312 B, zero swap. The old golden-master test reloaded development code,
  so later tests in that run do not establish installed-package authority.
- The single correction check **zKKMo0**: 1,945 passed assertions, one failed
  assertion, two errors, one warning and the same skip; peak 3,066,089,472 B,
  zero swap. The explicit installed-namespace identity assertion passes after
  removing the development reload. Stale import-environment injection was also
  removed: the static volcano implementation already qualifies `stringr::fixed`.
- Remaining new verification defects: the golden-master source-text assertion
  does not match the compiled native-pipe body; the bounded 1,024-row Glimma replay
  failed its 30-second initialization check. Its first malformed counts fixture
  and the subsequent timeout are both retained. No browser-memory reduction or
  successful deferred-rendering trial is claimed, and no application browser
  delivery code/dependency was changed.
- The `plotDensity(ProteinQuantitativeData)` -> `ggplot2::ggplot` S4 dispatch error
  occurs identically in separate fresh R processes on Q7BMJh and lb0zYc
  (`repair-existing-density-dispatch.json`). It is not introduced by these repairs
  and remains untouched. The tidyselect warning and unreadable cp08 skip are not
  converted to passes.

Each run retains `reader-tests.rds`, `reader-operation.json`, service logs and
terminal receipts; zKKMo0 also holds the `repair-*.json` observations. The initial
memory helper retained its first result during the second measurement and used
a different allocation threshold; use the corrected zKKMo0 observations, not a
ratio spliced across runs. Source backups before all repairs remain at
`/tmp/dia025-memory-repair-before.rm9UIE/R`.

The preceding repair budget was exhausted. The audit
rollback happened after zKKMo0; **lb0zYc is not the final working-tree build** and
must not be reused for the eventual GUI run. No rebuild or extra run was launched
in that turn after withdrawal.

### Verification-only continuation — 2026-09-21

The user then requested "continue please". DIA-025 records one focused
verification-only check and one correction/recheck, retaining the full-GUI pause.
No runtime R code, audit mechanism, application dependency or browser delivery
implementation changed. The existing reader-check service gained a fixed
`gui-memory-verification` selector for only the golden-master and isolated-browser
tests. The browser test moved into its own file so the already valid worker/QC/
normalization evidence need not be rerun. The broader repair selector still
includes it; the known predecessor density failure was not masked or reclassified.

The annotation test now captures the actual `glimmaXY` arguments and verifies
identifier/order, character annotations, missing-value cleaning, excluded math
columns and unchanged x values. It passes on the installed package, independently
of native-pipe deparsing. Its namespace remains installed throughout the check.
The Glimma runtime helper source in lb0zYc matches the current file exactly
(`0b0ee9d68bf654764f67d79a4f0aa9cedbc0708f4162cbac6a073ce56bf613c0`),
so this narrow evidence can reuse that library despite its rejected audit code.
It does not qualify lb0zYc for any audit or full-GUI run.

- **ZDGQdj**: 16 assertions passed, zero failed assertions/warnings/skips, one
  browser initialization error. Peak 1,185,300,480 B, zero swap. The failed page,
  dependency files and bounded observations are now retained instead of discarded.
- **wTPUSC** (the one correction): 16 assertions passed, zero failed assertions/
  warnings/skips, one browser initialization error. Peak 1,190,686,720 B, zero swap.
  The page is complete; jQuery and DataTables are loaded; one Glimma binding,
  widget, bound marker and plot container exist; **zero table elements or registered
  DataTables exist**. No JS error was captured. This does not establish that no
  error occurred. The conditional explicit static-render call was not invoked
  because the widget already had its bound marker. A fixture-only alert capture
  avoids an unobserved modal stall but reported no alert.

Both runs used the unchanged 1,024-row/12-sample/48-annotation fixture and the
30-second readiness limit under 8 GiB/900-second/zero-swap service guards.
Neither reached the deferred-rendering comparison. These are initialization
failures, not memory-limit failures or evidence against deferred rendering.
The retained HTML widget and JSON script identifiers agree; the remaining gap is
inside binding/render execution before table creation, not a missing dependency
registration. The exact cause is still unproved. Inspecting the retained binding's
`renderValue` directly with captured return/exception details is the next bounded
diagnostic; do not substitute a new full workflow or an easier data population.

All jobs are stopped and this continuation's correction allowance is spent.
The browser verification remains incomplete and the audit-memory bottleneck is
unresolved. DIA-025 and broader adoption gates stay held. The full GUI harness
remains paused as requested; no automatic further attempt is authorized here.

Final verification-only source bindings: golden-master test
`9c8ddd5efd3198b4736a5a5c37a0881217c227027b7f14e80e4278eb11308b7d`;
browser test `88c741ec5d3b99b68072e4addb2476a4232f9652552b9ec1aa530b33d2a734ff`;
driver `9d5af450d2617d9f1a9058e0b4cf4b147fbdcb51eb61814f448d1b9c7192a1a9`;
shell `fc2624f7ac47bf385dbf155b992b75aa5a76e2c2d02fd0f1f483e7967f309f25`.

## Retained diagnostic decision — lightweight audit reached DA views, then the memory guard

Replacement bWzeTc completed import/design, all seven peptide QC steps, real IQ,
protein QC, normalization including both Pearson plots, correlation Skip, full
legacy-session export, session reload and DA fitting. It then hit the 15.5 GiB
early guard while switching to the DA all-results table: current/peak
16,694,382,592 B (15.548 GiB), zero job swap, host available 9,635,487,744 B.
This was a controlled stop below the 16 GiB hard cap, not a completed journey or
an OOM. The explicit table/volcano/heatmap/Summary checks remain incomplete. No
independent large output verifier was run against an incomplete GUI receipt.

The instrumentation now works through long operations: no Rprof/Rprofmem files,
2,421,562 B total app/driver boundary logs, no recording-limit marker. This is
not the previous multi-GiB allocation-log problem. Required application output
file cache remains counted. The unit and the recorded app/driver/browser/worker
PIDs were confirmed stopped. All exported files/checkpoints and both lightweight
attempts remain intact; no additional run or runtime optimization was performed.

### Memory-first function audit (coverage is partial, not every heap owner)

| Function / stage | Actual memory observation in bWzeTc | Repair category / limitation |
| --- | --- | --- |
| `runPeptideQvalueApplyStep` and its audit helpers | Main-R RSS 2.85 -> 7.19 GiB; high-water rises 3.45 -> 8.00 GiB. Sampled group peak within the helper 9.54 GiB. | Full-table canonical strings, identity keys, sorting and serialization are confirmed allocation paths. Bound/reuse these preparations without changing digests, ledgers or scientific populations. |
| `runProteinPeptideApplyStep` / `.peptideQcCanonicalise` | Main-R high-water rises again to 9.25 GiB; sampled group peak within helper 10.59 GiB. Canonicalization exit records the new main-R high-water. | Another active audit peak, not cumulative resident checkpoint payloads. Exact live-object versus allocator ownership is still incomplete. |
| Rollup and subsequent peptide QC | Rollup main-R RSS falls 8.08 -> 6.27 GiB; later peptide helpers stay near 6.3 GiB without a new main-R high-water. | Historical manager/working payload counters are zero between actions. This does not prove every GUI closure or native allocation is released. |
| `runProteinReplicateFilterApplyStep` | Group sampled peak rises to 12.77 GiB while app RSS stays near 6.28 GiB. Four new workers persist through export, reload and DA. | Concrete resource-lifetime target: stop the step-created cluster after its last consumer, preserving configured execution counts and never killing a borrowed pool. |
| `buildProtNormCorrelationPlot` (pre/post) | Main RSS stays near 6.28 GiB; neither call raises the main-process high-water. Calls are too short for a reliable per-call cgroup peak in this sample stream. | Pearson is not demonstrated to cause a new resident peak for this 48-sample fixture. Allocations may reuse existing heap; wider data and correlation Apply/reapply are not covered. |
| Composite QC / normalization | Main RSS about 6.26 -> 6.29 GiB in this run; sampled group peak about 12.14 GiB for normalization. | Image decoding/raster construction is known from the small run; do not extrapolate its earlier RSS increment to this already-large heap. |
| `collectProtNormExportSessionData`, `saveProtNormExportArtifacts` | Main RSS 6.29 -> 6.32 GiB during collection, then steady during writing. Sampled group peaks 13.14 and 13.72 GiB respectively. | Full self-contained history was hydrated and written. Output file-cache costs are real; no reference-only substitution or hidden copy was used. |
| `restoreDIAStateSession` / `runTestsContrasts` | Main RSS 6.46 -> 6.52 GiB during restore and about 6.55 GiB during fitting; sampled group peaks 14.18 / 14.22 GiB. | Full external RDS load/validation precedes restore. This run did not show a new main-R high-water here. |
| DA browser views | One renderer's PSS rises from 202,427,392 B near DA fit completion to 1,859,792,896 B about 19 s later. The guard occurs during the table-significance switch. | Confirmed browser-side growth. Default Glimma is enabled and starts with the initial Volcano tab; the exact widget allocation needs isolation. Do NOT claim the table ships every row: its code selects four columns, limits rows and uses DT's normal paging. |

Numbers in the table are inclusive boundary observations, not additive allocation
budgets. VmHWM is cumulative main-process high-water; an unchanged value does not
mean zero allocations. Cgroup interval maxima are sampled and include concurrent
workers, browser and file cache. Maximum memory-sample gap was 0.4325831 s.

At the final concurrent process sample (epoch 1789973098.414104), app R PSS was
7,009,723,392 B (6.53 GiB); the four protein workers totaled 4,397,181,952 B
(4.10 GiB); the original eight workers approximately 1.12 GB (1.05 GiB); the large
browser renderer 1.73 GiB; and the driver R process 0.18 GiB. The nearest cgroup
sample has 2,009,759,744 B file charges (1.87 GiB). PSS/shared mappings and cgroup
file charges overlap and use different accounting; these must not be added into
a purported exact breakdown. At the last delivery the store had 16 states, zero
native/working/cache payloads, 62 metadata entries (378,592 B), and zero DB leases.

The complete audit list still lacks the remaining large DA interactions, heatmap
save and Summary. Completing them under the same limit now requires a measured
application-memory repair before another full run. The narrowest next candidate
is the retained step-created protein worker pool (~4.10 GiB at the guard), not a
new storage backend or a weaker memory threshold. Exact cleanup ownership,
reapply/revert/session compatibility and science must be checked before claiming
that headroom. This is outside the completed diagnostic-only allowance. Following
the user's request for memory-only cleanup and a written plan, the ordered scope,
checks and proposed verification budget are now in
[DIA-025 revision 1.5](../../../tickets/openai/DIA-025.md#memory-only-repair-work-plan).
That update is documentation only: no runtime cleanup or additional run has
started. DIA-025 remains held; the previous diagnostic allowance stays exhausted.

### Lightweight audit authority and verification history

User explicitly requests the complete GUI path before runtime refactors and
approves up to 60 minutes for the diagnostic. DIA-025 revision 1.4 removes the
unbounded allocation trace and deep object-size scanning. Scalar function/phase
records and external process/cgroup measurements remain. Same Q7BMJh package,
556 MB input and science, 16 GiB hard / 15.5 GiB early stop / zero job swap,
unchanged host-headroom guards. Only the explicit diagnostic flag gains time;
ordinary acceptance remains 900 seconds. No active proof job existed at entry.

First lightweight large member OLyg5I committed qvalue_filtered, but the driver
reported WebSocket EOF at 15:41:03 (about 180 seconds after the Q-value click).
The Shiny session subsequently ended at epoch 1789969466.84462, so advancing the
same GUI was no longer possible. It was explicitly stopped, not retried in place.
Pre-stop peak 10,384,080,896 B, zero job swap; all services are inactive. The
controller's subsequent missing-counter exit 76 follows explicit stop and is not
a memory-guard failure. The exact socket-close cause is not established.

Use this scope's one correction for the observed long-wait failure mode: the
memory diagnostic now polls existing DOM conditions using short DevTools commands
and rejects an absent transport response or disconnected Shiny session. Ordinary
waits remain unchanged. A real browser with a 190-second blocking Shiny action
must pass before the replacement full-size run. All limits/science remain fixed;
the original disconnected member and its incomplete coverage stay retained.

Correction check Cux4Uz passed, including the real 190-second blocking Shiny
action, transport/session-loss rejection and unchanged ordinary-wait delegation.
Peak 1,604,272,128 B, zero swap, main lifetime 210.847628 s. This validates the
replacement waiting strategy, not the exact cause of OLyg5I's disconnection.
The unchanged short small-GUI scientific/output evidence is retained; this test
specifically exercises the changed long-wait boundary. The one correction is used.
Replacement large member is now starting. No further retry allowance is assumed.
Corrected driver SHA-256:
`a9b94448d8a56aa4b72b5a979622196fcf5844d9ce2d205f090c4181af108447`.

One small verification then one large full-path attempt; one bounded correction
if needed, with all prior failures and hashes retained. No package/algorithm,
backend or worker-count change. The output is a memory-first audit with coverage
and uncertainty stated, not a CPU ranking or performance acceptance claim.

Lightweight checks fc0C7S (38 assertions) and final 9wq9Qj (40 assertions) pass
without failures/errors/warnings/skips. The latter includes bounded logging,
the real browser handshake/action/shutdown, expanded installed function names,
driver action labels and isolation of the one-hour allowance from ordinary runs.
Final-check peak 1,316,397,056 B, zero swap, main lifetime 18.517073 s. Small
full-GUI member flBseC is running; no large member has been launched in this
scope and no correction allowance has been used.

Small flBseC completed the full GUI with 16 states and a final diagnostic:stop
marker. Boundary logs total 2,061,344 B, with no truncation marker or Rprof/Rprofmem
files. Peak 7,041,634,304 B, zero swap; main lifetime 260.996498 s. Independent
verifier vDLr79 passes historical science, full session/output checks and identical
project PNGs against native 6Mc16q (peak 1,090,867,200 B, zero swap). These are
instrumentation/correctness checks, not qualified performance ratios. The approved
same-input large memory member is now starting; no runtime source changed.

Frozen helper SHA-256:
`ff0271deb881cdc049190a07044f9dde6a40b6b0fe38e0a4853051e655730048`;
driver `697db8c2ae685746479ad9a2fe28fd73dc55168f50f501cd789908edfc4e63a5`;
R runner `12e6317cfc79c8422dad0433a748f950df1421c84948a922d311a69d3d656f92`;
shell `813ced6c92a2aca3f4a7a7ea4f487b466a1c5b33ec4e6bf27d42e0aa7db0572a`.
Application and driver boundary logs have explicit 64 MiB bounds. No Rprof or
Rprofmem output is produced. Shallow holder metadata does not traverse closures
or compute object.size. A log cap reports incomplete coverage rather than hiding it.

### Retained allocation-profiled timeout and findings

Large diagnostic ckolpO terminated at the unchanged 900-second service limit
(main lifetime 900.070553 s; Result=timeout, SIGTERM). It did not OOM or trip the
memory guard. Exact cgroup peak 15,661,907,968 B (14.586 GiB), CPU 1,940.114389 s,
zero job swap. There is no finite completion-time denominator. All proof services
and their workers were confirmed inactive after controller cleanup.

Q-value and Protein Evidence committed; precursor SQL completed but its audit
had not committed when timeout arrived. Four durable states exist: initial,
raw_data_s4, qvalue_filtered and protein_peptide_filtered. At the last completed
delivery, native payload history/working map/cache/DB leases were zero. This is
not a census of all GUI roots. Large IQ, later QC, normalization/Pearson, legacy
export/load, DA and Summary were not reached. No large science/output acceptance
or runtime memory fix is claimed. No 1/2 GB case or extra large retry ran.

### Named costs and interpretation

| Measured path | Evidence | Meaning / proposed repair boundary |
| --- | --- | --- |
| Peptide audit: `.peptideQcCanonicalise`, `.peptideQcDigest`, summary distinct counts and removal-ledger keys | 559.66 s of 829.26 s represented by sampled main-R stack records contain `.emitPeptideQcAuditRecord` (about 67.5%). Nested digest/canonicalization 315.76/275.82 s, summary 107.58 s and removal ledger 136.32 s. GC appears in 204.98 s, overlapping these paths. | Main early bottleneck. Bound/reuse audit preparation while preserving exact keys, digests, ledgers and populations; do not remove audit or alter science. |
| Whole-object hydration | 74.68 s of inclusive sampled main-R CPU; 7,188,491,616 B of cumulative >=5 MB allocations in the checkpoint-hydration stack prefix | Repeated complete reads still cost time and allocation. This is cumulative allocation, not retained memory or unique bytes. Metadata-only callers and per-action reuse need targeted review. |
| Protein replicate worker lifetime (small full journey) | Four additional workers at about 0.860 GB PSS each remain after their step, alongside the original eight approximately 0.140 GB workers | `runProteinReplicateFilterApplyStep` creates a new cluster; current candidate retains the handles as borrowed. Bound its lifecycle without changing the configured execution population/count. Not measured at large size. |
| Composite rendering (small full journey) | 31.32 s inclusive main-R CPU; main-app sampled PSS grows from roughly 1.39 to 1.99 GB during composite generation/rendering | `loadProtNormImageAsPlot` / `generateProtNormCompositeFromFiles` decode PNGs into raster grobs and rebuild/save the composite. Separate graphics allocation from Pearson calculation and worker residency. Not measured at large size. |

Large traced helper intervals (profiled, not whole-button performance claims):
Q-value Apply 344.64068 s; Protein Evidence Apply 251.92004 s, including 5.30020 s
to persist the rollup projection. The saved precursor calculation took 10.37978 s
and then the unchanged audit required full-state hydration (4.29747 s) and resumed
large table summary work. Moving that relational calculation did not remove the
surrounding audit bottleneck. No statistical calculation or population was bypassed.

Largest sampled main-app PSS was 9,472,355,328 B (8.822 GiB), versus imported S4
476,131,952 B. The largest sampled R-heap counters correspond to 6,656,309,104 B
inside audit canonicalization of Protein Evidence. Conversion follows installed
utils::summaryRprof: vector counters * 8 plus node bytes. These maxima occur at
different instants and are not subtracted or summed. Allocation logs include a
single 898,561,008-byte serialization allocation under `.peptideQcDigest`.
At completed Q-value delivery the only inventoried workflow-bus field above 1 MB
was the 4,762,520-byte FASTA table; other reactive/closure/native/allocator ownership
was not exhaustively measured. Empty manager counters do not imply low whole-R RSS.

### Diagnostic overhead and remaining uncertainty

The allocation trace itself grew to 7,556,935,680 B; Rprof was 85,987,612 B. At the
largest current-memory sample, anonymous charges were 10,076,209,152 B, file charges
5,470,896,128 B and kernel charges 56,307,712 B. Counters are asynchronous and not
an exact additive peak decomposition. Independent fincore observations confirmed
several GiB of allocation-log residency; not every file charge is uniquely assigned
to that file. The diagnostic's object.size hooks also appear in allocation records.
Therefore the 14.586 GiB peak and 900-second timeout do NOT qualify unprofiled
16 GiB workflow capacity or performance. The separate original unprofiled 8 GiB
failure remains valid. Future full-journey attribution should bound allocation
trace windows/volume while retaining whole-run CPU/process/phase observations;
this is a proposed follow-up, not an unrecorded change or retry here.

The 41,463 large Rprof sample-prefix records have 40,025 continuation lines from
multiline anonymous call labels. Stock summaryRprof could not parse memory rows;
raw files are retained. Reported function counts use quoted names on sample-prefix
lines and can omit ancestors in continuation lines. They are inclusive, diagnostic
main-R samples, not exclusive wall times or a reconciliation of the combined
1,940.114389 CPU seconds. Worker/thread/browser/controller-related residual costs
remain unattributed. Process RSS/PSS and cgroup values are not interchangeable;
per-process maxima are never added. Maximum memory-sample gap was 0.406034 s;
the process sampler recorded no file-access gaps, but is coarser than peak accounting.
Zero job swap does not imply the already-running host had no swap usage.

Next bounded implementation proposal: audit preparation in
`R/func_prot_qc_peptide_audit_helpers.R`, starting with bounded canonical/identity
key temporaries and exact reference equivalence fixtures. Preserve all substantive
provenance/digests and the existing full action boundary. Worker lifecycle and
graphics work remain separate measured follow-ups. The diagnostic-only allowance
is used; no further runtime correction or expensive rerun is inferred. DIA-025
and the production/performance gates remain held pending an approved repair scope.

### Verified small full journey

Small diagnostic 6BLNvj completed through Summary and teardown. Independent
verifier 2eJh7w passes the complete legacy history, required output inventory,
independent DA and the preserved native comparison: all historical science and
project PNGs agree. Diagnostic service-tree elapsed 303.735044 s, CPU 672.894853 s,
peak 7,795,953,664 B, zero swap. Profiling perturbs these values; they are not a
performance-qualified comparison. Verifier main lifetime 19.127192 s, peak
1,090,539,520 B, zero swap.

Completed-small observations, not extrapolated large results: the original eight
configuration workers persist; protein replicate filtering starts four more
workers (first seen at epoch 1789965714.831331), each subsequently about 0.860 GB
PSS, retained through DA/Summary. The manager reports borrowed groups 8 and 4.
The main app's largest sampled PSS is 1,990,971,392 B during composite rendering;
the service peak also includes those worker pools, browser and profiling-file
charges. Do not sum individual process maxima or call the whole peak plot memory.
Composite rendering has 31.32 s of inclusive main-R CPU samples; GC has 90.62 s
across the small journey (overlapping categories). Function/phase observations
identify costs, not a disjoint allocation/retention census.

Small raw Rprof contains multiline anonymous-function call labels, so summaryRprof's
memory parser failed on those continuation lines. Preserve that limitation and
the raw file. Direct inspection counts 11,354 sample-prefix records plus 589
continuation lines; quoted-function first-line counts are bounded, potentially
incomplete inclusive observations. No replacement run or rewritten raw profile.

### Startup-message correction authority and bindings

User "yes", 2026-09-21, authorizes one further diagnostic-only correction/recheck:
leave Shiny stderr untouched, verify a real browser handshake/action/shutdown,
repeat the small full-GUI diagnostic, then execute the unused same-input large
diagnostic. Existing Q7BMJh library, science, inputs and 16 GiB diagnostic limits
remain unchanged. Reuse unaffected aEn94b harness checks; rerun only the affected
instrumentation tests. No active proof job existed at continuation entry.

The approved message-stream correction passed focused service 1Uc5Cv: 22
assertions, no failures/errors/warnings/skips; actual browser connection, action
and shutdown were observed with diagnostics enabled. Main lifetime 21.231150 s,
combined peak 1,343,057,920 B, zero swap. The extra correction is now implemented;
no additional retry is assumed.

Successor helper SHA-256:
`af03ce750ae25c04c5dfd8c56610933e65a399a209cdd123208d58298ce4ffbe`;
driver `3ec98c224b888a55988b6a94e712140fcb13b84bc233f578b1aa7fae6f5d5b61`;
shell `8d398918c4024d60e4a9a5055b130fc464cc5b441697f2636e11b837144e6a71`.
The R runner and installed package are unchanged. Stdout is mirrored; Shiny's
stderr remains with its existing owner. Direct phase/profile files survive an
interruption, but an abrupt kill can still prevent the final browser-log bundle.

### Retained startup-instrumentation failure

Small diagnostic MZumM9 was explicitly stopped after identifying a harness
startup deadlock. `startCheckpointGuiDiagnostic()` redirected the R message
connection with `sink(..., type = "message")`. Installed shinytest2's
`app_start_shiny()` waits for `Listening on http` in the child's original stderr
file before connecting the browser. That line appeared only in the new
application-console.log. The app was listening, but no browser connected and no
import/QC/DA action ran. This is an instrumentation failure, not a new workflow
memory failure. No full-workflow culprit attribution or runtime fix is claimed.

Before the explicit stop: cgroup current 971,280,384 B, observed peak 981,008,384 B,
zero job swap, CPU 17.068189 s. These are pre-stop observations, not terminal or
completion metrics. The controller subsequently exited 76 with missing counters
because the unit had been deliberately removed; that secondary guard result is
not the primary cause. All proof services were confirmed inactive. Startup
profiles, process samples, phase records, logs and the original running receipt
are retained without rewriting the attempt as successful.

The one correction was already spent on the fixture issues below. Under the
declared budget this continuation is blocked. Next authority needed: one
diagnostic-only correction/recheck that preserves Shiny's original message
stream, tests the actual browser startup handshake, and repeats the small
diagnostic before the already-authorized 550 MB / 16 GiB member. Leave memory
limits, science, package build, input and earlier failures unchanged. The large
diagnostic allowance has not been used; the 1/2 GB cases remain held.

### Approved diagnostic scope and verification history

2026-09-21 user authority: same 556,028,391-byte report, full existing DIA IQ
journey, 16 GiB hard cap / 15.5 GiB early stop / zero job swap / 900 seconds.
Require >=18 GiB host available at launch and stop below 2 GiB. Reuse Q7BMJh
unchanged for the first diagnostic, with harness-only profiling. No large resident
control, no 1/2 GB run, no gate promotion. Existing 8 GiB failure remains valid.
No live proof services existed at entry. Host available was approximately 23 GB;
the runner must recheck at launch. Existing host swap is not job swap.

Bounded allowance and file surfaces are in DIA-025 revision 1.2: instrumentation
lifecycle checks, one small GUI diagnostic, one same-input large diagnostic,
one evidence-directed correction pass and affected small/large recheck. Preserve
every attempt. Name observed functions/owners/process costs; do not equate total
allocations with retained bytes or sum per-process RSS into a combined peak.
Any full-workflow stages not reached must remain explicitly unmeasured.

Instrumentation check SltEeo stopped with two failed assertions and one cleanup
error (peak 1,509,912,576 B, zero swap). The dynamic fixture function appears as
`<Anonymous>` in Rprofmem, although its 16,000,048-byte allocation is recorded;
the test incorrectly required the function label. The other fixture scheduled
diagnostic shutdown after its temporary directory was removed. Correct the
assertion to require the actual allocation and order cleanup using withr::defer.
The new allowance's one correction is consumed by this verification repair;
preserve the failed test receipt. No large member has run yet. Runtime science
and the installed Q7BMJh library remain unchanged.

Corrected check aEn94b passed all 792 assertions with no failures, errors,
warnings or skips (247.305008 s main lifetime, peak 1,300,201,472 B, zero swap).
This verifies instrumented success/failure, installed trace names, actual large
allocation capture, process RSS/PSS sampling and affected harness routing. It
does not certify a GUI journey. Small diagnostic GUI check is the next member.

Frozen diagnostic helper SHA-256:
`d8500146db0568c88c81d9e8e5598323e56bd14db983314fc672d1f9b1b1f2b2`;
driver `879ea1b6e54f02c83e9b71cc17ea2b3e32b0fdcbb305cfc6e28443c60f7edac5`;
R runner `3a1d08437fdddb23eff8bd9289620e93e3277e4d6835f068d33f5b875a07b900`;
shell `4622e1f2b160e684c5780b61d9188f1ecc06f0e826c970111efa818c2dc2c0dd`.
The shell's final change records a vanished per-process list as an attribution
gap; mandatory cgroup accounting still fails closed. Protocol/scenarios are unchanged.

### Retained 8 GiB negative capacity result

Candidate `browser-smoke.4rIoEg` on the exact preserved 556,028,391-byte report
was stopped by the 7.5 GiB early guard. Controller exit 76; observed current
8,075,481,088 bytes, peak 8,075,743,232 bytes (7.52 GiB), zero job swap. Host
MemAvailable was 15,409,541,120 bytes (14.35 GiB). Journal evidence records a
controlled stop, not an OOM kill; the unit and all proof jobs are inactive.
This is a valid unfavorable capacity outcome, not an infrastructure replacement
or finite completion-time denominator. The 1/2 GB cases were not launched.

The last delivered state was raw_data_s4 (3,840,000 rows, 14 columns,
476,131,952 bytes object.size), with only initial + raw state references and no
completed QC checkpoint. At that delivery, native history and action working
maps were empty and no DB lease was held. Those are between-action observations;
they do not establish the active working set at the later guard instant.
Pearson, the DuckDB rollup and full legacy export had not been reached.

Near the stop, the last sample recorded 7,576,559,616 bytes anonymous memory,
410,292,224 bytes file charges and 36,184,064 bytes kernel charges. Counters are
asynchronous and not an exact peak decomposition, but this was predominantly
process-owned memory, not file-cache growth. In the final half-second anonymous
memory continued growing while file charges stayed flat. No retained app console
log identifies the exact allocation; report the interval as after imported-state
delivery and before first QC completion, not a proven exact Q-value hotspot.

Important implementation limit: whole-object hydration is still used for current
R operations. Imported data, design-initial, design-cleaned-data, design-result,
workflow data_cln and raw S4 exist as distinct saved records; scoped reuse is by
entity identity, not shared underlying component. Multiple active reads can
therefore coexist and lose allocation sharing. Disk descriptors do not prove
those representations remained resident at the guard. GUI/worker overhead,
other retained references and transformation/audit/serialization temporaries
remain unattributed; zero history-manager entries cannot be equated with a
bounded whole-GUI heap.

Conclusion: the actual small GUI path is verified, but the same architecture has
not achieved the required larger capacity under the unchanged limits. There is
no approved correction budget left. Stop scaling. The next necessary work is
bounded attribution of this early live-memory window, including per-process
ownership and R allocation/retained-root evidence, before proposing another
repair. Do not change backends, sample populations, formats or limits by inference.

### Verified small pair and scale-up history

The approved verifier-only correction passed. Focused installed test bT6AcX has
11 passing assertions and no failures/errors/warnings/skips, peak 1,067,491,328
bytes and zero swap. Corrected verifier nugUNX passes the original small GUI pair:
all 16 historical states including initial, scientific table types/missingness,
substantive metadata and every project PNG. The existing session/output/DA checks
also pass. It recomputes DA independently; neither GUI journey was rerun.
Lifetime 19.196603 seconds, peak 1,110,835,200 bytes, zero swap. Receipt SHA-256:
`4dcb8dd48a6a439e0cdb1dd34a3c27ba757874315fd749b036e0fb2930377683`.

Continue only the previously unexecuted candidate cases, in the approved order,
through existing guarded services. Keep the original stopped outer receipt and
all failures unchanged; this section binds continuation instead of replaying
completed members. No correction allowance remains; the infrastructure-invalid
replacement remains unused. The installed Q7BMJh package is unchanged.

Reference binding Kmj9K8 passed, matching the original report/FASTA/design hashes:
556,028,391 report bytes, 3,840,000 rows, 10,000 proteins and 48 samples. No new
large input or S4 was generated. Its peak was 623,624,192 bytes, zero swap.
The candidate-only GUI run on those exact files used the existing
8 GiB/zero-swap/900-second guards and subsequently stopped as recorded above.
Neither 1/2 GB case has run.

Current verifier-only driver SHA-256:
`3d46093993c6d6b1b926a7ec8d72970c3d894bd15c7b9470190e2bd9e138522f`;
shell `81b8ae86a71b048eea656154a8561cc6ffb6e735ac3edd38d29a1817f5a25820`.
Only comparison namespacing/metadata preservation and the focused test selector
changed. The package, scientific methods, GUI action boundary, input recipes,
limits and R runner are unchanged; earlier timing evidence is not relabeled.

### Retained verifier stop and approval

Continuation authority (user "yes", 2026-09-21): one additional verifier-only
correction/recheck is approved. Qualify utils calls and directly test the paired
comparer in its isolated environment, including NULL/attribute preservation.
Use the existing guarded reader-check service with a fixed focused selector.
No package rebuild, GUI rerun, scientific relaxation or prior-budget reset.
The completed small pair remains authoritative; larger cases stay gated on its
comparison. The historical failure and stopped outer receipt remain unchanged.

User approved the full plan on 2026-09-21. Actual installed Shiny, IQ only;
1–2 GB denotes on-disk DIA-NN reports. Whole qs2 checkpoints plus the existing
DuckDB precursor rollup; no other estimator/storage change. Production defaults
and all historical performance gates stay unchanged.

Both small actual GUI journeys completed through Summary and teardown. The
candidate's separate verifier passed the existing full-session/output/DA checks,
then failed in the new cross-member comparison on an unqualified `object.size`
call. That driver is deliberately sourced with baseenv() as parent; the two
calls in compareCheckpointGuiCapacityMembers require utils::object.size.
The cross-run comparison is incomplete, not a scientific mismatch or a pass.
The one correction allowance is exhausted, so no repair/retry was performed.
No retained-OOM-input or 1/2 GB candidate was prepared or launched. All services
and the outer dispatcher are stopped. DIA-025 is blocked; production defaults
and historical gates remain unchanged.

| Small actual GUI member | Full service elapsed | Combined CPU | Combined peak | Job swap |
| --- | ---: | ---: | ---: | ---: |
| Resident 6Mc16q; verifier E6cS13 passed | 229.069107 s | 482.568236 s | 6,941,270,016 B | 0 |
| Candidate gk3vOF; verifier 12Yl3P stopped | 235.079248 s | 500.732380 s | 7,011,180,544 B | 0 |

Candidate point values are 2.62% slower, 3.76% higher CPU and 1.01% higher peak.
These are a single small correctness pair (3,557,875-byte report), not qualified
performance inference, a capacity improvement, or evidence about 1–2 GB reports.
No performance gain is established. Both members verified installed-package
loading; the candidate's explicit activation was observed with 16 states,
zero resident native-history records, zero working records after delivery,
zero optional cache and zero database leases. These counters cover named owners,
not every possible GUI/global/plot allocation.

The candidate verifier retained successful checks of 108 IQ/QC files, 13 Summary
files, full 16-state legacy content, a 256-row/12-sample final table, independently
recomputed 256-row DA and heatmap/required publication copies. It then stopped
while recording sizes for the first scientific history comparison. Remaining
cross-member states and PNG comparisons were not executed. Its lifetime was
14.192099 seconds, peak 1,060,081,664 bytes, zero swap; this is a verifier defect,
not an OOM or a performance-member failure.

Next authority needed: one verifier-only correction/recheck, preserving all
existing members. Qualify the two utils calls and add a focused test that actually
executes this paired comparer in its base-environment context, preserving NULLs
and attributes. Reverify the unchanged small pair; do not rebuild the package or
rerun either GUI member for a verifier-only change. Only a passing comparison may
release the already-planned candidate-only larger cases. Do not reset the failed
pilot receipt or the correction allowance by creating a new attempt directory.

The outer receipt is `gui-capacity-21cd5929b2aed616/capacity.json` under the proof
root below; it correctly records `stopped`. Input preparation l6OHWn, native
6Mc16q/E6cS13 and candidate gk3vOF/12Yl3P all use the browser-smoke prefix.
Member/verifier SHA-256 values are:
`faabcc7efd44ce4fa6e98b2eb2249de851dcabc6426c362a28b1a25d6a0e0456`,
`6deba42519ea1c7bf3a0fcb9ef3012da05cafac1eba0fba82e7540675969e2bc`,
`ee9a288afe10639dafb5ec3b70ca6ceef572aabecc789d2c9b24eb88187aadc7`.

### Implementation/check history

No live proof jobs existed at entry; loaded failed units are historical evidence.
No large run or large control has been launched. One correction and one
infrastructure-invalid replacement was allowed for this entire experiment.
The correction was used and passed; the infrastructure replacement remains unused.

Installed first-slice build `browser-smoke.fPsDaH` passed in 46.986127 seconds,
peak 1,124,929,536 bytes, zero swap. State/installed compatibility verification
`browser-smoke.bifLjr` passed (81 assertions, no errors/warnings/skips), 29.887753
seconds, peak 1,889,484,800 bytes, zero swap. Coverage includes scoped single
hydration/release, ordinary public history replacement, inactive audit branches,
failed payload/head publication, corrupt payloads and project incarnation changes.
Both services stopped normally. No correction/replacement allowance was used.
Subsequent GUI wiring changes require their own affected recheck and installed
build; this first-slice receipt does not certify them.

### Wiring check failure and bounded correction

Build `browser-smoke.TzL0oj` installed successfully. Affected verification
`browser-smoke.cvgaqF` has 2,528 passing assertions, two errors, no assertion
failures/warnings/skips. Service lifetime 42.986539 seconds, peak 2,005,598,208
bytes, zero swap. Errors: export's new scope accessed an atomic test-only
workflow stub; the new real-QC test omitted the working fixture's required
column/scientific configuration, so native Q-value execution lost Run before
any SQL comparison. These are not OOM or performance results.

Use the one correction to preserve the injected native export contract and
bind the same actual scientific fixture configuration/local update environment
as the earlier working sequence. Also finish the planned same-action reuse at
the GUI write/read boundary: a just-written qs2 record remains in the required
action scope until delivery, avoiding immediate duplicate hydration. The cache
still clears at delivery and the optional cross-action cache remains zero.
Rebuild and rerun affected checks once. If verification still fails, stop under
the exhausted correction allowance. No browser or larger candidate has run.

Correction build `browser-smoke.Hsk6Cy` and check `browser-smoke.OJ7MrA` passed:
2,547 assertions, no errors/warnings/skips. Check lifetime 46.331838 seconds,
peak 2,035,716,096 bytes, zero swap; service stopped normally. The real three-QC
native/candidate values and audit records agree. The native injected export
contract, same-action write reuse, cold GUI-holder contract, atomic external
history ingestion and complete native RDS export checks pass. No correction
allowance remains. Continue only the planned harness/remaining lifecycle wiring;
any later verification failure is retained and stops work without another repair.

### Remaining harness implementation

The existing `--mode gui` entry now recognizes only the separate
`tests/testdata/workflow-checkpoints/gui-capacity.json` scope for this experiment.
It dispatches the five declared selections serially through existing services;
preparation, every member and every independent verifier are recorded before
advancing. It refuses to repeat an existing attempt directory for the same build.
The ordinary paired evaluator and frozen scenarios/protocol are unchanged.
Large resident capacity execution is rejected even at the internal browser entry.

Small-pair comparison covers all legacy historical states (existing worker/parser
portability, explicit timestamp fields and member-relative metadata paths), with
1e-12 scientific tolerance, table types/missingness and scientific PNG contents.
The normal output verifier still checks the full required output/session/DA
surface. This is not an independent IQ estimator or a universal scientific gate.
The generated reports increase distinct peptides at fixed proteins/design, never
padding or repeated aliases; generation streams <=200,000 rows per chunk and
does not pre-hydrate a full report/S4 outside the journey.

Final candidate build `browser-smoke.Q7BMJh` installed successfully, zero swap.
Expanded check `browser-smoke.IQFpkf` passed 3,325 assertions, with no failures,
errors, warnings or skips. Service interval 315.297096 seconds, peak
2,534,563,840 bytes, CPU 300.240722 seconds, zero swap. This includes the existing
harness, actual three-step scientific/output regressions, explicit capacity
routing/stop checks, project replacement and affected module/installed contracts.

The three real sampler/accounting lifecycle probes passed: success `nICf0g`,
intentional main exit 7 `sPfMLu`, and retained child `DauOfP` (all browser-smoke
prefix under the same proof root). Peaks were 399,417,344 / 399,699,968 /
399,773,696 bytes, zero swap. Exact terminal evidence and successful owned-unit
cleanup were retained; the deliberate failure is not an infrastructure retry.
The child case observes service-tree exit after main-process exit. No allowance
was reset or additional correction used.

Historical pre-run sequence: one final installed build; affected state/module/harness checks;
three bounded accounting lifecycles (success/failure/retained child) through the
same sampling path. Only on success execute the declared GUI capacity sequence.
No automatic repair or larger run follows a verification failure. Keep every
terminal service receipt. Source/driver/build/input hashes are frozen by the
installed build, preparation manifests and capacity receipt before each member.

Execution freeze: installed RDB
`21cd5929b2aed616afb94d36d535a556ef5a73ad75b551d7a2594710a04251e5`;
driver `f96700eaf9e5f0f93b53ce49bf19614282aed00c05e7ec72a9f483f507dbc359`;
R runner `3e79ef945c9aee5de46aa5d29e49c2d457d6afe4da9401eb8451688711209c16`;
shell `6dedc6677f6203acf49b53151dfdb2e27fd3d638d64f197cf87bdcf8f36008f6`;
capacity contract `f6d9514fa56e1d415c33dfc43a950a3b9c7bc191ce6ea044f3ddbb988c7dce23`.

Unverified lifecycle observation from source review, not an adoption claim:
close clears registered callbacks/working inputs/workers, but on_commit still
references the workflow bus. Check release of that closure before claiming
complete old-project heap release. This was not changed as an adjacent repair;
the passing replacement test proves scope rejection, not exhaustive heap release.

Run order: focused/lifecycle checks; one small installed resident/candidate
correctness pair; candidate on the exact original failed input; only on success,
candidate on each larger recipe. Each member requires independent verification.
Stop scaling on a valid failure. All services retain 8 GiB/zero swap/900 seconds,
7.5 GiB early stop and existing host-headroom/serial guards. No 16 GiB inheritance.

## Pre-change bindings

Branch performance-v2; source parent janitor
`c7c12851ea4d9e96f91df0e1f9e7b91c2eb51017`. Dirty tracked runtime patch SHA-256
`5e14d61e181e504252fa776592ab16a5d7fda0fcbeab31b322dd63265983fd1d` does not include
untracked package helpers; it is not by itself a candidate-build fingerprint.
Prior complete installed prototype is `browser-smoke.2Oee9H/library` beneath
`/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5`.

| Existing surface | SHA-256 before DIA-025 |
| --- | --- |
| shell runner | 4ed3d5a91ef322753f43d26706848d01f640ddac2d3424c3755b0478c37abf0c |
| R runner | 368a2fffbb063ff37c885d1c475f2cc15d7de39b49bd7e2fdbe8fadf06220fc0 |
| driver | 251515ade9047f34528309575a845f29a05643652ff226b4d431a4be6a8ab76b |
| WorkflowState | a89d321cf23c32fbdaf2779c28bc4e2304b2724b187e7247ff2921bb180636e5 |
| saved-entity core | cc26c074b0d86a342f4b76a91de39b9f077d47265048d0af356087605417eafd |
| DuckDB rollup | 3f25c4252b3e69fa62b9a280df43d82e2a9afdc4818098b9c4a40d8360aeed6d |
| frozen scenarios | cfdf7da2fc5a6123da269f85011550d8ffd2e8a611de071e0cd1a52a6f35d152 |
| frozen protocol | 517edccb0906165698f3c1f165dfe7d961a17910d93363a60e843fa466c62379 |

The exact 3,840,000-row/48-sample OOM input is the original `representative`
report (556,028,391 bytes), with the baseline.md companion bindings and failed
members bQWqlM/zqsj4M. They remain negative-capacity fixtures, not finite time or
peak denominators. Prior three-step results remain in checkpoint-results.md;
they are not full-GUI evidence.
