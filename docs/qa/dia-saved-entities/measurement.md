# DIA installed-package measurement harness

DIA-003, in progress. This is measurement tooling, not runtime optimization.
No G0/scientific/performance acceptance is claimed.

## Durable controller completion evidence — 2026-09-14

The diagnostic shell runner now retains complete service stdout/stderr and the
systemd-run terminal summary in request-owned service.log. It writes unit,
operation, service_exit and log_exit to controller-status.txt before releasing
its serial lease. service_exit is the systemd-run --wait exit status; launch/
signal/cause interpretation still requires the accompanying service log. The
service's nonzero status is preserved; a successful service with a failed log
pipeline returns 74 instead of reporting success. No resource limits changed.

Permanent shell-orchestration tests use command doubles, not scientific/service
acceptance substitutes. They verify success, service exit 7, simulated log exit 1,
retained logs/status, overlap refusal and lease release. The old ad-hoc lease
check failed once (78350) and passed when traced (29077); its flock probes can
race acquisition. A readiness-marker test now verifies overlap without acquiring
the lock being tested. Neither the initial failure nor the old fixture was erased.
Harness runs passed 148 then 153 assertions, zero failures/warnings/skips
(64118 and 63739 respectively).

Real service VdtVmt (40411) rechecked the retained bad Vpqc21 export and failed
for its expected protein-count inconsistency. Both the failure log/terminal
summary and service_exit=1, log_exit=0 are on disk. Runtime 13.740 seconds, CPU
12.740 seconds, reported peak 863.9M, zero swap; terminal I/O reported 12.3M read.
This is failure/accounting evidence, never a usable timing denominator.

Real passing verifier yCfKWB (78184) retained its full log and service_exit=0,
log_exit=0. Session/DA/heatmap checks passed. Runtime 15.984 seconds, CPU 14.796
seconds, reported peak 926.5M, zero swap. Terminal summaries remain human-readable
and rounded, not substitutes for precise cgroup peak accounting.

The log-copying tee process runs in the outer controller scope, not the measured
service cgroup. Its buffers/file-cache charges and controller-loss recovery are
not yet accounted for in a qualified completion boundary. Current operations are
diagnostics; baseline timing dispatch still refuses execution. Full controller/
process-tree accounting, exact terminal reconciliation and recovery remain required
before performance acceptance. No package runtime source changed.

## Resource-accounting continuation — 2026-09-14

Guarded receipts now retain cpu.stat, memory.stat, memory.events, memory.current
and device-scoped io.stat snapshots at entry/exit. Missing or malformed optional
counters have explicit unavailable/error states, not fabricated zero values.
memory.stat mixes byte components and event counters; its components can overlap
and cannot be summed into memory.peak or interpreted as a simultaneous peak
breakdown. Device I/O lines remain separate rather than summing possibly stacked
devices. Exit snapshots precede final receipt writing/process teardown, so they
are not authoritative terminal service totals or exclusive phase attribution.

First parser test run 23237 failed because unqualified setNames was unavailable
in the isolated test environment (120 other assertions passed). After qualifying
stats::setNames, 74549 passed 133 assertions, zero failures/warnings/skips.
Tests cover malformed/duplicate/nonfinite counters, missing evidence and raw
device-line preservation. No runtime R package code changed in this work.

Real verifier P65BAV (79040) passed functional checks but had no io.stat interface.
CPU/memory/events were available; peak was 972931072 bytes, sampled CPU delta
14.32602 seconds versus the service's 14.768-second terminal CPU report. The gap
is not assigned to a scientific phase. A diagnostic assertion requiring all
counter interfaces failed and is retained as missing-accounting evidence.

The user-service hierarchy advertises the io controller. The runner now explicitly
requests MemoryAccounting=yes and IOAccounting=yes for its own services, keeping
8 GiB/zero swap/900 seconds unchanged. A CPUAccounting=yes request was reported
deprecated and ignored by this host; it was removed from future commands because
CPU counters were already available. That warning remains in run evidence.

Verifier RpgKtP (76537) passed functional checks and all four counter interfaces
were available. Its io.stat snapshot was empty: no device entries were observed
in this cached verification, not proof of zero logical I/O. Peak 971583488 bytes;
service runtime 17.087 seconds, terminal CPU 15.353 seconds, zero swap. Earlier
receipts without accounting are not retrospectively upgraded. These diagnostics
do not qualify complete I/O attribution or a performance gate. Terminal-controller
reconciliation, retained allocations, phase costs and full timing dispatch remain
unfinished.

## Heatmap driver continuation — 2026-09-14

The DIA-003 browser driver now opens the shipped Heatmap tab, waits for a loaded
image or an output error, reads actual contrast/top-N/scaling/clustering settings,
and invokes Save Heatmap. It does not alter the statistical population or inject
new settings. Harness 54357 passed 116 assertions with zero failures/warnings/skips;
shell syntax and whitespace checks pass. No runtime package source changed.

The save helper catches image/metadata errors internally, so the save notification
is only UI evidence. The required saved products are PNG, PDF, methods markdown,
and cluster CSV when cutting produces row clusters. Significance filtering
precedes top-N absolute-fold-change selection; source expression-row ordering is
retained. File contents and cluster assignments require separate verification.

Browser HvWSGZ (73737), using installed VEqJKP and unchanged IE9ohP inputs, is
running the full existing small diagnostic with the heatmap/save extension. Its
eventual passed_da_ui_only status cannot alone certify saved heatmap products.

HvWSGZ completed the heatmap and save actions with observed defaults: top-N 50,
both-axis clustering, row scaling, ward.D2, Euclidean distance, k_clusters=4,
RdBu, contrast treated_vs_control. Runtime 230.692 seconds, CPU 424.815 seconds,
cgroup peak 3364294656 bytes, zero swap. No runtime source changed for this step.

The independent verifier now checks this fixed default heatmap case using the
verified DA table and quantitative S4: significance selection, descending absolute
fold-change top-N selection, original quantification-row order, row scaling and
existing nonfinite-to-zero behavior, stats::dist/hclust/cutree. It compares all
exported molecule IDs and integer cluster assignments exactly. PNG decoding and
2400x3000 dimensions, PDF signature and core methods metadata are also checked;
the PDF check is not full rendering validation. A requested heatmap contrast that
matches no verified DA result cannot silently bypass these checks.

Verifier 7YTe4S (47679) passed: 7 selected proteins, 12 samples, 4 clusters;
session/history and all DA TSV/XLSX checks also passed. Runtime 15.559 seconds,
CPU 14.094 seconds, cgroup peak 976244736 bytes, zero swap, outside the browser interval.
Its receipt binds PNG/PDF/methods/CSV hashes. The PNG was visually inspected and
shows seven rows, twelve samples, group annotation, dendrograms and four cluster
colors. No pixel-level parity or alternative-setting acceptance is inferred.

The final harness run 9266 passed 120 assertions, zero failures/warnings/skips,
including rejection of changed/reordered molecule IDs and cluster assignments.
Full-case coverage, report/download compatibility and complete measurement
orchestration remain open; this does not close DIA-003 or any performance gate.

## Approved interactive-volcano key repair — 2026-09-14

The user approved “thats fine i approve”. DIA-003 revision 1.5 records the exact
helper/test scope. Uniqueness, count alignment, annotation row names and final
count synchronization now use args_row_id. Protein.Ids remains annotation;
significance thresholds, zero-FDR/nonfinite handling and sample-group matching
are unchanged. No estimator or scientific export code changed.

The new real-S4 test has four distinct groups, two sharing one accession, with
quantification rows and sample columns deliberately reordered. Before the repair,
five alignment/annotation/coordinate assertions failed (84942: 46 pass, one
existing tidyselect warning, one unreadable cp08 skip). After the repair, the
volcano suite passed 51 assertions. The test was then strengthened to call the
real Glimma constructor while capturing arguments, rather than returning a fake
widget. Combined volcano/harness run 34942 passed 166 assertions, zero failures,
with the same existing warning and skip retained.

The driver now switches from static to interactive volcano, requires the shipped
Glimma SVG and table, clicks a table row, and waits for selection plus linked
expression symbols. Selectors were checked against the installed Glimma widget
JavaScript. This extends functional coverage; it is not yet a passed live result.

An isolated installation VEqJKP (68156) is in progress. Parent source c7c12851;
runtime patch SHA-256 9c762b9f8f72051a72aacd89608e06ca2222965a8601392d6d214f6764981038.
Original interactive failures remain in earlier browser logs, not erased or
retrospectively passed.

VEqJKP installation completed in 96.342 seconds, CPU 86.851 seconds, reported
peak 1G, zero swap. Installed Rdb SHA-256:
310fd8fdae006bea35c6913105fa797876a1d92c45680ea0f4b6d0c8b68cf3e8;
Rdx: a567f2078b39de58c01a462e39d8f602a8bfee1928c9693ccee9ad06ad119e48.
Browser member 8KEaCi (1263) is running the unchanged small scientific settings
plus the new interactive render/table-selection/expression-update check.

8KEaCi completed passed_da_ui_only in 200.662 seconds, CPU 367.380 seconds,
cgroup peak 3375542272 bytes, zero swap. Its widget rendered and a table row
became selected. The original selector counted 67 SVG symbols, including legend
symbols; it also used a broad expression-symbol presence check. Those observations
must not be presented as 67 proteins or conclusive expression-data verification.
The driver now selects the named Vega data marks (.mark-symbol.marks) for both
panels and records actual expression-point count. L3RGqX (16872) is rerunning the
same installed build/inputs to qualify this corrected observation boundary.

The real-widget test now also checks four unique widget gene identities when
two groups share one accession. Run 10554 passed 169 combined volcano/harness
assertions with zero failures, retaining the known warning and cp08 skip.

The corrected observation run L3RGqX completed passed_da_ui_only: exactly 64
scatter data marks and 12 linked expression data marks after table-row selection.
Runtime 241.266 seconds, CPU 449.270 seconds, cgroup peak 3381903360 bytes, zero
swap. The earlier 67-symbol count remains recorded as including legend symbols,
not a protein count. No elapsed-time difference is interpreted as a performance
gain: these are diagnostic runs, not matched accepted timing blocks.

Separate guarded verifier dc4hoI (19583) passed session/history and DA TSV/XLSX
agreement for VEqJKP/L3RGqX: 16 history entries, 64 proteins, 12 samples and all
33 reference DA columns agreed with the same installed-method replay; no qvalue
warnings. Runtime 16.565 seconds, CPU 14.740 seconds, peak 922898432 bytes, zero
swap. DA TSV hash remains 0bdc7e2e5b25cf9425a5533a56eb283bd09da04f5b5b4d5972144b25a9eab140.

Focused source coverage 46791 measured 95.41284% for
generateProtDAVolcanoPlotGlimma. Artifact /tmp/dia003-interactive-key-coverage.rds.
Its standalone environment emitted the fake-class package-name warning and
retained the unreadable cp08 skip. This is helper line coverage, not whole-package
or performance evidence. The authorized identity repair is verified for the
tested interactive route; heatmap, broader settings and the full DIA gates remain
outstanding.

## Interactive-volcano identity audit — 2026-09-14

The retained Ou9TOT browser log contains an interactive Glimma failure even
though the separately checked static panel rendered: buildXYData/data.frame
raises “arguments imply differing number of rows: 2, 0”. The prior passed_da_ui_only
receipt never certified the interactive view; it remains unqualified.

R/func_prot_da_volcano_glimma.R filters/deduplicates by Protein.Ids and then
intersects/matches those accession values against counts_df[[args_row_id]]. The
server supplies the actual declared Protein.Group key. A small read-only
inspection of Ou9TOT's real exported DA table and session found 64 result rows,
zero accession-to-quant-key matches and 64 active-key-to-quant-key matches. Thus
the adapter discards the entire plotted population during counts alignment.
The Glimma call error is consistent with that empty mapping, not failed DA fitting.

A bounded adapter repair must use the declared key for uniqueness, matching and
row synchronization while retaining Protein.Ids as annotation. It must not rename
quantification identifiers, refit models, change significance thresholds or
deduplicate distinct protein groups by shared accession. The helper and matching
tests are outside the approved static-wiring/label-default scope; an explicit
extension is required before editing them. No runtime change was made in this
audit, and the original failure logs/artifacts remain unchanged.

## Approved static-label defaults — 2026-09-14

The user approved the additional fallback with “yes”; DIA-003 revision 1.4 and
the active plan record the server-adapter/test scope. The renderer now supplies
TRUE and 10 only when the corresponding UI inputs are NULL. Explicit FALSE,
zero and custom counts remain unchanged. The plotting helper, scientific
calculations and UI controls were not modified by this fallback.

Before the adapter edit, the new argument-matrix tests failed four missing-input
assertions, with 241 assertions passing (12323). After the two-line edit the DA
module suite passed 245 assertions, zero failures/warnings/skips (34515). The
harness now selects the static view and requires a loaded nonempty plot image or
reports the output error; its suite passed 114 assertions (53445), also clean.

An isolated corrected installation JnZQrf (57966) is in progress. Parent source
is c7c12851, runtime patch SHA-256
5ddb361af2a0a4d1b3e701810fe3066964a7cda09731077bcee2589a09123ce6.
Earlier installed browser evidence lacks these UI/adapter changes and is retained
only for its original scope. No live static-view success is claimed yet.

JnZQrf installation completed successfully in 87.774 seconds, CPU 82.749 seconds,
reported peak 1G, zero swap. Installed Rdb SHA-256:
fa15dfb7f2bfaade5a304b28ce9096f6d46d2c7bf3ff959abc4cbcb52ae261c3;
Rdx: 66eb5b4e31333b5efc454f28a4c2f5a3306c9ceed48fc9d2a710ca78c0c2c5ed.
Browser member Ou9TOT (81658) is running the existing fixed small DIA route,
now including a switch to the static volcano panel. No nonexistent label inputs
are injected; it must use the actual production adapter defaults.

Ou9TOT completed passed_da_ui_only with static_volcano explicitly observed in
the action list. The browser switched the shipped UI to static mode and observed
a loaded plot image without an output error. Runtime 203.338 seconds, CPU 373.103
seconds, cgroup peak 3365998592 bytes, zero swap. This is a small functional
diagnostic under the fixed limits, not a matched performance acceptance member.

Separate guarded verifier yOQyc7 (82048) passed the session/history and DA TSV/XLSX
checks for JnZQrf/Ou9TOT: 16 history entries, 64 proteins, 12 samples and all 33 DA
reference columns agreed with the installed-method replay, without qvalue warnings.
Runtime 14.317 seconds, CPU 13.085 seconds, cgroup peak 921133056 bytes, zero swap.
Normalized and DA TSV hashes remain identical to the earlier verified outputs.

Focused source coverage 14576 measured 100% for da_server_volcano_render_handler
and mod_prot_da_volcano_tab_ui using the existing DA module suite. The artifact
is /tmp/dia003-static-wiring-coverage.rds. This is scoped line coverage, not
whole-package coverage or timing; earlier broader coverage limitations remain.
The approved output-ID and absent-label-input repairs are now verified for this
installed static-view route. Interactive volcano, heatmap, remaining settings,
full harness accounting and performance gates remain open.

## Static label-input diagnostic — 2026-09-14

After the approved output-ID repair, the next read-only check confirmed that
the renderer explicitly passes input$volcano_show_labels and
input$volcano_label_top_n, but the production UI declares neither input. The
helper defaults are TRUE and 10; explicitly supplied NULL bypasses those defaults.
Its if (show_labels) then raises “argument is of length zero”.

Diagnostic 98809 loaded the installed I4uVca helper and x298E3's actual DA TSV,
called the renderer's existing helper with the actual absent-input values NULL,
and asserted that exact error. This is a small unguarded diagnostic outside
timing, not a browser or performance result. A sandbox timedatectl startup
warning was retained. No source/runtime repair was made by the diagnostic.

The approved revision 1.3 scope was only the UI output id and connection test.
Supplying the helper's existing label defaults when these controls are absent
requires extending that repair to the server adapter and its tests. Explicitly
provided label settings should remain honored, with no new UI, statistics or
plot-helper changes. Until authorized and verified, static-view completion stays
held; the passing connection test alone is insufficient.

## Approved static-volcano wiring repair — 2026-09-14

The user approved “Yes, fix the wiring”. DIA-003 revision 1.3 and the active plan
record the bounded UI/test scope. The protein UI now uses volcano_plot_static,
matching the existing server renderer. Renderer code, helper calculations,
scientific values, public function signatures and other omics are unchanged.

The regression test now checks that both registered volcano output names occur
in the namespaced generated UI, and rejects the old unmatched static id. Before
the UI edit, session 5189 failed exactly those two wiring assertions (231 other
assertions passed). After the one-line UI edit, 33339 passed 233 assertions with
zero failures/warnings/skips. This preserves explicit red/green evidence rather
than relying on the old renderer-only test. git diff --check passes.

This verifies the source-level connection, not live installed rendering. Existing
I4uVca browser evidence predates the UI correction and cannot qualify it. The
renderer also reads volcano_show_labels and volcano_label_top_n, which are not
declared in the inspected UI; the mocked renderer test does not validate that
production-input path. Check it during live-view qualification before claiming
the entire static panel works. No further renderer change is authorized or made
by this one-line wiring repair.

## Remaining-view source audit — 2026-09-14

The protein DA static volcano UI declares plotOutput(ns("volcano_static")) at
R/mod_prot_da_ui.R:170. Its sole protein server renderer instead assigns
output$volcano_plot_static at R/mod_prot_da_render_volcano.R:31. Repository-wide
search found no protein volcano_static renderer. The existing isolated module
test asserts volcano_plot_static, so it does not test the UI/server connection.
This source-level wiring defect prevents claiming that the static GUI panel is
integrated; a valid exported PNG is not evidence for that panel. A bounded runtime
wiring repair and corresponding regression test require an explicit scope
amendment beyond the previously approved identifier/validity repairs. Approval
has been requested; no runtime change was made for this finding.

The download_da_results handler writes placeholder text rather than ZIP content.
This is already explicitly recorded in contract.md's table/view consumer section.
Preserve that baseline behavior and do not certify it as a scientific export;
implementing a real archive would be a separate feature decision. The per-contrast
TSV/XLSX outputs verified earlier are distinct from this placeholder download.

Interactive volcano returns the existing Glimma widget; heatmap has its own
rendered output and Save Heatmap action. Neither has yet received complete browser
interaction verification in this harness. No new GUI run was launched just to
repeat the source-proven static binding mismatch.

## Separate DA export comparison — 2026-09-14

The allowlisted verify-da operation runs after the browser has stopped, under
the same 8 GiB/zero-swap/900-second limits. It reuses the session/history checks,
then replays the installed registered differentialAbundanceAnalysis method from
the exported pre-DA S4, with the same per-contrast table and parameter setup as
the GUI. The existing helper's weights=NA and no-TREAT path are unchanged. Only
the verifier's fresh working directory is available for incidental plotting;
it does not invoke the publication writer on the original project.

Comparison covers every reference result column, preserving row count/order,
identifiers, comparison/numerator/denominator, missingness, NaN versus NA, finite
patterns and infinity signs. Numeric tolerance is 1e-10; no rows are reordered
or discarded to obtain agreement. The extra annotation-only gene_name column
is not compared by this scientific-table checker. Original artifacts are read,
not rewritten. Tests reject reordered/missing rows/columns, changed quantities,
opposite infinities and NaN/NA substitution. Harness 81007 passed 112 assertions,
zero reported failures/warnings/skips, with a sandbox timedatectl startup warning
outside the reporter. Shell syntax and whitespace checks pass.

Guarded verifier Z9Jbwb (36260) passed_small_da_export_checks for browser x298E3,
bound to installed I4uVca. Both TSV and XLSX matched the separate method replay:
64 rows, 33 reference columns, including all p-value/FDR/fold-change and sample
quantity columns plus accession provenance. No qvalue_warnings were returned.
All earlier normalization export/history checks also passed. Runtime 12.608
seconds, CPU 12.145 seconds, cgroup peak 919580672 bytes, zero swap. This verification
is excluded from the browser completion interval, not a performance replicate.

The existing per-contrast volcano PNG was visually inspected: it decodes and
shows the contrast title, fold-change/FDR axes, threshold line and colored points.
That does not certify exact graphical content or interactive browser behavior.
Per-file hashes remain in Z9Jbwb's verification-receipt.json and the x298E3 record.

This is GUI/export agreement with the same installed scientific implementation,
not an independent mathematical oracle or the later full scientific gate. Broad
settings, interactive/static view behavior, heatmap/downloads, complete file
inventory, recovery and matched performance acceptance remain outstanding.

## DA browser continuation — 2026-09-14

DIA-003's allowlisted browser-da operation extends the existing normalization
diagnostic with the actual production Load Filtered Session and Run DA Analysis
buttons, then the all-results table. Inputs are explicit: formula ~ 0 + group,
q-value threshold 0.05, TREAT cutoff 0. The installed registered scientific route
is unchanged; no standalone replay substitutes for browser execution. The driver
does not certify DA publication from the analysis-complete flag.

Harness runs 89419 and 47915 passed 105 assertions each, with zero reported test
failures/warnings/skips. The first emitted a sandbox timedatectl startup warning
outside the test reporter. Added tests cover session-load-before-analysis order,
the declared settings and table-error propagation. Shell syntax/whitespace pass.

Installed run x298E3 (33886) used I4uVca and unchanged IE9ohP inputs. It completed
passed_da_ui_only, including session reload, DA and the all-results table.
Service runtime 223.356 seconds, CPU 434.118 seconds, peak 3359387648 bytes,
zero swap. Logs report successful DA output writing and 64 table rows for
treated_vs_control. The per-contrast TSV has 64 data lines and the expected
Protein.Group/log2FC/raw_pvalue/fdr_qvalue fields. TSV SHA-256:
0bdc7e2e5b25cf9425a5533a56eb283bd09da04f5b5b4d5972144b25a9eab140;
XLSX: ff8d1ae43ce2d35d9d6560022abd68823a119a73c3ab1bf6896953f9beb150fa.
Per-contrast volcano PNG/PDF, combined volcano PDF and significant-molecule count
outputs are present. Presence/row count is not content or scientific agreement.

The source audit confirmed DA catches file-writing errors without failing its
analysis status, and table failures can render a normal Message table. The driver
now additionally requires the statistical table headers rather than accepting
any populated table. x298E3 had already loaded the earlier driver, so this added
header check is unit-tested but not retrospectively credited to that run.

Still required: independent DA scientific/output verification, interactive/static
volcano and heatmap consumers, downloads, all specified settings/reapply/reset
cases, complete publication boundaries and measurement orchestration. This is
neither a complete GUI baseline nor a performance acceptance result.

## Independent export verification continuation — 2026-09-14

The session-count correction was installed separately at browser-smoke.I4uVca
(39148), parent c7c12851 plus runtime patch SHA-256
a5cf4a4a55f6333d644a8650b4b9a236d08d9d1a8ffe2c059a7be4a99960be7d.
Installation passed in 96.719 seconds, CPU 92.032 seconds, reported peak 1G,
zero swap. Installed Rdb SHA-256 is
89614f27233968ee351c3ae66406573c30825c793cc1d747e89051a2bbc8d61f;
Rdx is 17053135b394f1804fc6400c3fddf1605c11b55443564f8094a751b279df63ba.

The reusable runner now allows verify-normalization LIBRARY OUTPUT_ROOT MEMBER_DIR
as a separate guarded process. It requires a matching completed UI diagnostic and
library path, fingerprints five required export files, checks identical timestamped/
latest file hashes, validates every scientific state in the explicit 16-entry
history, compares final/correlation/current S4 serialization and TSV identifiers,
order, missing/nonfinite patterns and quantities (tolerance 1e-12), checks counts
and summary text, and exercises R6 state installation/revert. It is limited to the
small cyclicloess/Skip/Skip diagnostic and cannot certify the full output inventory,
DA consumers, arbitrary external sessions or performance gates.

Harness tests reached 97 assertions with zero failures/warnings/skips (40808),
including deliberately wrong counts, missing history, invalid raw input and changed
quantities. A subsequent test adds independently hydrated environment provenance;
that addition still requires execution.

First real verification gBPpsS (47782) failed an overly strict in-memory identical
comparison between two separately decoded sessions. Both files are byte-identical
(SHA-256 b9c42a8a83e71a85a66d7ea0fce854d2596be42c6ceeff87b3e998c42ca4b46d);
their embedded environments have different identities after decoding. The verifier
now compares the file hashes and serialized S4 contents instead of environment
addresses. This was a verifier defect, not a new data corruption finding.

Corrected verifier DeDxNV (20727) rejects Vpqc21 with the expected exported-protein-
count inconsistency, after validating its scientific states and S4 agreement.
Both failed verification receipts are retained. The corrected installed build is
now running the unchanged browser-normalization diagnostic as yKUxrL (24068).

yKUxrL completed passed_normalization_ui_only in 219.036 seconds, CPU 433.361
seconds, cgroup peak 3370020864 bytes, zero swap. Independent verification gCDobu
(1968) then passed_small_export_checks: 16 history entries, 64 proteins and 12
samples, valid scientific S4 states, byte-identical timestamped/latest session
files, serialized S4 agreement, TSV quantity/identifier/order/missingness checks,
summary counts and R6 restore/revert. Verification was outside the browser member:
13.377 seconds, CPU 12.642 seconds, peak 899534848 bytes, zero swap. Do not add its
peak to the browser peak or treat it as an independent journey replicate.

The passed verification receipt contains five file hashes. Timestamped/latest
session SHA-256: 0ee38c202c6320815bb18993c8eac1dceffb5ad052f8cabdca5ee8f355714fd3.
Normalized TSV SHA-256: b698cd7dac2c2c2003c2874bcff4375fbb75bb125e386c1aac7704497348bf4a;
cmp confirms this TSV is byte-identical to Vpqc21's before the count correction.
The old metadata failure and verifier-authoring failure remain preserved, not
overwritten or retrospectively passed.

This verifies only the small cyclicloess/RUV Skip/correlation Skip export scope.
It does not qualify all normalization choices, RUV modes, DA load/views, complete
automatic plot/report inventory, malformed external-session recovery, resource
attribution or paired performance gates. DIA-003 remains in progress.
Final harness rerun (26682) passed 99 assertions with zero failures/warnings/skips,
including independent hydration of environment-containing provenance. Shell syntax,
ticket schema/status consistency and git diff --check were also rechecked.

## Normalization driver continuation — 2026-09-14

DIA-003 now has an allowlisted browser-normalization diagnostic extending the
same real import/design/peptide/protein sequence. It explicitly selects cyclicloess,
RUV Skip, correlation Skip, and the production full-session Export action. It
waits for actual production notification/result text and treats error notifications
as failure. No scientific helpers are replayed in place of browser actions.
The new action sequencing/error checks bring the harness suite to 91 assertions,
zero failures/warnings/skips (66983). Shell syntax and whitespace checks pass.

Run 4760 / browser-smoke.Vpqc21 uses the corrected e3lvCL installation and unchanged
IE9ohP inputs. It is a diagnostic, not qualification of other normalization/RUV/
correlation choices. Its possible success status is passed_normalization_ui_only:
file existence, contents, complete inventory and independent restore verification
must still be checked after the browser stops. UI success alone does not certify
output publication or whole DIA completion.

The export consumer audit found another hard-coded Protein.Ids count in
collectProtNormExportSessionData, producing zero for Protein.Group tables. The
existing authorized repair scope was extended explicitly to that file, and the
source now counts the declared key. A regression test compares session count and
all saved scientific objects/history. Vpqc21 continues against e3lvCL unchanged;
its evidence must not be represented as testing this later source correction.

Vpqc21 completed passed_normalization_ui_only: service runtime 265.100 seconds,
CPU 540.623 seconds, cgroup peak 3367288832 bytes, zero swap. Its required
normalised TSV/S4 RDS, timestamped/latest session RDS and summary files are present.
File discovery must use --no-ignore for this external project; default rg ignores
some generated outputs and its absence results are not authoritative inventories.

Independent small-file inspection (52402, separate unguarded diagnostic process,
not a timed/scientific gate) loaded the e3lvCL package and the exported latest RDS:
the final S4 validates, has 64 protein rows, and the export contains 16 history
entries from initial through correlation_filtered. final_protein_count is zero,
confirming the export-metadata defect. This contradicts output correctness despite
UI success; no normalization/export acceptance is granted. Complete history/TSV
content verification and external restore integrity remain outstanding.

The corrected session-count source and its real-S4 regression test must be
reinstalled/requalified before Vpqc21 can be superseded. Keep Vpqc21 unchanged.
Source regression run 30429 passed 584 assertions, zero failures/warnings/skips
(dia-baseline-compatibility and prot-05b-norm-module-contracts). It verifies the
corrected exported count alongside unchanged current/correlation S4 objects and
complete saved-state list. No timing assertion is derived from that test run.

## Implemented so far

- Strict CLI parsing for the seven required options and help/validate-only flags.
- Frozen-protocol validation, including approved G2/G4 numbers, service limits,
  resource-accounting fields, budgets and the paired uncertainty procedure.
- A pure paired-receipt evaluator with independent process/pair identities and
  exact route/input/protocol/boundary/dependency/build bindings.
- Separate accepted-evidence, failed and held outcomes. The evaluator always
  returns promotion_authorized=FALSE; a numerical result is not a verified gate.
- Explicit library structure and physical output-root checks, without loading
  MultiScholaR during preflight. Structural installation metadata deliberately
  reports load_verified=FALSE; fresh loading must happen inside the measured service.
- Synthetic-receipt tests for pass/fail/borderline/exact boundary, stale bindings,
  missing/nonfinite metrics, failed files/science, OOM/timeout controls, profiling,
  repeated-session misuse, G2 latency and resource ceilings.

The shell wrapper now also provides the allowlisted --service browser-smoke
operation through its persistently approved command prefix. It validates the
library/root, creates a fresh external request directory and enforces the fixed
8 GiB/zero-swap/900-second service limits. It does not accept arbitrary shell
commands. --help works
without loading MultiScholaR. All workload modes, including baseline, still refuse
execution and --validate-only because their installed driver/bindings are not yet
complete. This is deliberate refusal, not a successful preflight. No source-load
fallback, fake scientific route or helper-only GUI claim is permitted.

## Baseline source and isolated installation

Pinned B0-original source: c7c12851ea4d9e96f91df0e1f9e7b91c2eb51017.
A Git archive of that exact commit was extracted into a new empty directory:

/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/b0-install-RNROup/source

Installed library:

/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/b0-install-RNROup/library

No uncommitted ticket/harness work was included in the archive. No dependency was
downloaded or updated, and the existing user library was not replaced.
R CMD INSTALL completed successfully in its own service with MemoryMax=8G,
MemorySwapMax=0 and RuntimeMaxSec=900. Service runtime 95.929 seconds, CPU 92.281
seconds, approximately 1.1 GiB peak, zero swap. Installation is preparation, not
part of a timed GUI replicate or a performance gain.

Installed artifact SHA-256:

| Artifact under library/MultiScholaR | SHA-256 |
| --- | --- |
| DESCRIPTION | 572dde1a140994e4adf66dbe44eeac7bf246061a601c608d6118031dafa34d44 |
| NAMESPACE | 9dc606a759229234e12f5a0b55879af5b7230a42f1177a4d86fa1e088879ccc5 |
| R/MultiScholaR.rdb | 2ab24707f5e18bf5ef2932f02f63caab3e524d8cd83410bcec958030ad9dc213 |
| R/MultiScholaR.rdx | 618e5565ca1e1da90947bc8658c7e37e357bd8ecf40d8c1cbd1eb87872d30c3b |

The first fresh installed-load diagnostic failed because its diagnostic command
passed unsupported where=ns to methods::selectMethod. The package had already
loaded; this was a diagnostic-authoring failure, not proof of a package failure.
That failed service remains evidence: exit code 1, 11.782 seconds, 818.3 MiB peak,
zero swap. A bounded correction passes the generic explicitly via fdef instead.
Neither installation nor fresh-load checks qualify the complete IQ journey.

The corrected v02 diagnostic passed from the explicit isolated library. It
observed 624 exported names, checked WorkflowState/runTestsContrasts/run_app,
both peptide/protein class registrations and protein normalization/DA method
selection, and exercised the initial resident WorkflowState. Service runtime
12.410 seconds, CPU 12.138 seconds, approximately 818.8 MiB peak, zero swap.
This is a fresh installed compatibility smoke, not a GUI benchmark or exhaustive
public-API test. Its namespace list includes pkgload because the existing package
imports development dependencies; that is not evidence of load_all execution.
The explicit loaded package path was checked against the requested library.

A byte comparison confirmed the archived source matches all 346 current R and
DESCRIPTION/NAMESPACE files; the runtime tree is unchanged from the pinned commit.

## Browser action audit

The current helper-e2e-browser.R launcher calls pkgload::load_all and sets
MultiScholaR test_mode=TRUE. It cannot be reused as the installed out-of-the-box
performance launcher. Current interaction helpers may be reused only after
checking their inputs and completion observations. In particular:

- preserve original actual UI actions; do not replace a journey with scientific
  helper replay;
- do not turn on test-mode scientific/annotation fallbacks to obtain a pass;
- audit how to observe production completion without test-only state digests;
- retain DIA-001's exact UI-order, protein-key and Q.Value validity risks;
- required GUI session loading and automatic outputs stay inside completion;
- independent output/scientific verification runs in another process, not by
  retaining an extra full history in the timed service.

The baseline driver is not yet implemented. Limpa execution stays held until its
own qualification; IQ evidence cannot certify Limpa.

Driver primitives now include a cold installed launch closure with only two
scalar path bindings, explicit disabling of scientific test mode/checkpoint
capture, production stepper observations (not test-only digests), and the
source-verified peptide action/result IDs in current UI order. The complete
action orchestration and output verification remain unfinished.

## Real browser startup evidence

Two initial startup attempts are retained: v01 refused because NOT_CRAN was not
set; v02 failed in shinytest2's development-namespace lookup before app launch.
The latter occurred because the controller's working directory was package
source despite requesting an installed app. The reusable runner now starts in
a fresh external request directory. NOT_CRAN enables the browser testing library;
it does not enable MultiScholaR scientific test mode.

The subsequent installed production landing smoke passed. Its receipt and logs
are retained under:
/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/browser-smoke.hFLr0k

The receipt records passed_landing_only, peak_bytes=1316274176 and zero swap.
Service runtime was 24.255 seconds, CPU 25.942 seconds. Startup, browser and app
were in the same service. This is not the complete DIA workflow, a final combined
peak acceptance measurement, or a demonstrated performance gain. The complete
controller must still capture terminal service accounting and all required files.

Reusable invocation (the user saved persistent approval for this exact runner):

    bash /home/doktersmol/Projects/APAF/MultiScholaR/tools/profiling/run_workflow_checkpoint_proof.sh --service browser-smoke LIBRARY OUTPUT_ROOT

The launcher uses the supplied CHROMOTE_CHROME or an available local browser;
browser identity still needs binding in full measurement receipts. Current smoke
used /usr/bin/brave. No dependency was downloaded or updated.

## Measurement boundary to implement

### Production import continuation

The reusable runner now additionally accepts the allowlisted browser-import
operation with explicit REPORT and FASTA paths. It uses real project setup,
native shinyFiles selection payloads and the existing import observer, then
observes the production Import stepper. It does not enable standard-upload
test-mode substitutes or call import scientific helpers from the controller.

The first production import run passed_import_only using the pinned commit's
small DIA report and FASTA. Its receipt/logs are retained at
/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/browser-smoke.DwHRzS.
The service reported 25.208 seconds, CPU 56.113 seconds and zero swap; the
in-process peak record is 2490327040 bytes. These are diagnostic costs, not
candidate/control gains. No full workflow, design, scientific checkpoint or
required output inventory is qualified by this import-only outcome.

The exact native chooser path payload, including a filename containing a quote,
round-trips through shinyFiles::parseFilePaths in tests. The focused suite now
passes 76 assertions with zero failures/warnings/skips. Shell syntax and diff
checks pass. Ordinary failure logs are retained, and cleanup errors cannot
silently leave an otherwise success-shaped smoke receipt.

## Full measurement boundary still to implement

### Installed protein-rollup consumer failure

The next production diagnostic retained all seven completed peptide states, then
failed in the protein-rollup action. Its receipt/logs are at browser-smoke.FAORgA
in the external root. Service runtime was 102.002 seconds, CPU 199.324 seconds,
peak 2742398976 bytes and zero swap. No successful workflow/performance outcome
is inferred from that failed run.

IQ produced iq_output_file.txt with Protein.Group and all twelve aliased sample
columns. The failure is downstream: updateProteinIqRollupOutputs passes the
protein quantity table to updateProteinFiltering; that calls countUniqueProteins,
whose data-frame branch only accepts Protein.Ids and throws 'No Protein.Ids
column found'. The active Protein.Group key is legitimate and was preserved by
the IQ producer. Do not rename/duplicate scientific identifier columns or skip
the required tracking output to conceal the error.

Source owners: R/mod_prot_qc_protein_rollup.R, R/func_prot_qc_filtering_progress_helpers.R
and R/func_prot_rollup_count_helpers.R. These runtime files are outside DIA-003's
current write scope. A separately authorized, bounded baseline correctness repair
is needed before full successful-route qualification. The untouched installed
B0-original and failed evidence remain preserved. Related active-key and
post-rollup validity risks were already recorded by DIA-001/002.

### Matched synthetic diagnostic input preparation

The approved prepare-small operation reused DIA-002's unchanged small recipe
and generator in a separate preparation service. It wrote 4608 precursor rows,
12 design rows including biological replicate numbers, and 64 explicitly
synthetic FASTA records. Every reported stripped peptide was verified to occur
in its declared protein's synthetic sequence. Inputs are retained under
browser-smoke.IE9ohP in the external root with an inputs.json manifest.

SHA-256: report a8d1c758e061386332bb0cc949091f91517c2ebcf86db87fe0c8db33cc411105;
design 495f419ffc57824ed92ace15f3370168f5e8442316c8171cfeea4cac4204b4d8;
FASTA f22b5c5e82d1dff663a47c66543783db9cb6a39456bb464df8a747b5e64c7165.
The old E2E report and its failure were not modified. This input preparation is
not a timed workflow member or a biological/scientific acceptance result.

The driver now assigns declared technical replicate groups through the actual
builder controls. The first matched-input run (browser-smoke.OuKxfT, unified
session 94906) remains under observation. It exposed an artificial driver wait:
set_inputs(wait_=TRUE) expects a rendered output even for input-only controls or
unchanged values. That caused repeated 30-second waits. Future calls now use
wait_=FALSE followed by wait_for_idle, while Apply/completion observations remain
explicit. The current process was not restarted or counted as a valid timing run.

Terminal continuation: that process completed passed_peptides_only. All seven
states were observed in order: qvalue_filtered, protein_peptide_filtered,
precursor_rollup, intensity_filtered, sample_filtered, replicate_filtered,
imputed. Service runtime was 340.675 seconds, CPU 490.907 seconds, in-process
peak 2690334720 bytes, zero swap. The inflated wait behavior remains documented;
this run is functional diagnostic evidence only, not a performance comparator.
The regression suite passes 83 assertions, including the input-only wait fix.
No runtime source or scientific method was changed.

### Design and first peptide diagnostic

The production design-builder run reached passed_design_only and retained its
actual WT/KO-to-sanitized-run mapping at browser-smoke.PPLKVJ in the external root.
It used normal factor, metadata and contrast controls, with no invented select
options. Service runtime was 35.555 seconds, zero swap. This is not full baseline
acceptance. Its initial kind label still says import smoke; the explicit terminal
status identifies the completed scope, and later receipts use a distinct kind.

The following peptide diagnostic failed at qvalue_filtered. The pinned old E2E
report lacks Global.Q.Value, Global.PG.Q.Value and Proteotypic, which the original
production confidence filter requires. Its failed receipt is retained at
browser-smoke.2NfXms (37.287 seconds, peak 2527260672 bytes, zero swap). Two dynamic
tab input-binding errors also appeared; the driver now waits for controls and
clicks real tab anchors instead of assuming those inputs are already bound.
No scientific field was added to the old fixture and no filter was bypassed.
The next diagnostic must use the valid DIA-002 synthetic input and matching
design/FASTA assets, explicitly separate from the invalid historical E2E input.

The focused suite reached 80 passing assertions before the subsequent tab-wait
correction. That correction and subsequent driver work still require verification.

Source inspection on 2026-09-14 found that the shell runner performed metadata
preflight and created unique run directories but did not acquire the serial-run
lease previously claimed here. That overstatement was corrected before adding
the initial lock implementation; the earlier diagnostics had no enforced lease.

The runner now acquires a nonblocking flock on a persistent lock file in the
validated output root, before creating a run directory. Its shell retains the
descriptor while systemd-run waits. An overlapping invocation using the same
root fails rather than queues into an uncontrolled measurement interval. The
lock file is not deleted, avoiding inode replacement races.

An isolated shell test with stub Rscript/systemd-run commands passed: overlap
refused with exit 2, subsequent invocation accepted, service exit 7 propagated,
and lease released after failure. Exact test and stubs are retained at
/tmp/dia003-lease-check.Ur2oEk/check.sh. This validates shell coordination only;
it is not a real service, scientific or performance result. bash -n also passed.
Different output roots, externally launched jobs and a controller lost while
its service survives are not yet covered. Full owned-job exclusion and recovery
remain required; no diagnostic is promoted to a qualified timed member.

Required before timed acceptance: the outer controller must enforce serial
execution. Each timed member must use a fresh resource-limited service containing
the R driver, browser, app and descendants. Process/package/worker startup,
required writes/copies and completion acknowledgement must be inside the interval,
with no prewarming outside it. Complete boundary accounting remains unimplemented.

Use cgroup v2 memory.peak as the authoritative combined peak, with shared pages,
file-cache and kernel charges accounted by that same scope. Never add process
RSS to it. Record process RSS separately only as diagnostics. Preserve failed or
missing counters rather than reporting zero. Read the limits/counters while the
service exists; completed transient units may be collected.

External oracle/verification and fixture generation are separate from timed
members. Diagnostic profiling runs separately. Retained unique R allocations,
exclusive phase costs, CPU/I/O and decoded/fitted/rendered counts still need
implementation and validation. Unattributed costs and overlap must be explicit.

## Evidence and remaining work

The first harness suite passed 54 assertions; after nested-protocol and Limpa
refusal checks it passed 58, with no failures/warnings/skips.
After output-root and installation checks, 63 assertions pass with no failures,
warnings or skips. A running candidate receipt is held rather than misclassified
as a terminal failure. CPU completeness is checked even where G4 has no CPU ratio
threshold. The explicit installed-library --validate-only command still exits 1
with 'mode implementation not yet available: baseline'; this is expected partial
implementation refusal, not a successful validation or failed timed member.
Driver-closure, production observation and action-order checks increase the
focused suite to 73 passing assertions. Two JavaScript patch-construction
attempts hit template-literal quoting errors before filesystem edits; corrected
literal patches were applied. No failed runtime or preflight evidence was removed.
Command: Rscript --vanilla -e 'devtools::test(filter = "workflow-checkpoint-harness")'.
Shell help and bash -n passed. git diff --check passed.

A shell-wrapper creation attempt initially hit JavaScript interpolation of a
Bash variable; the R file had been written, the shell file had not. The corrected
literal was then applied. No runtime file or existing fixture was changed.

Remaining before ticket completion: installed baseline browser driver, actual
root/library/asset validation, serial service orchestration and durable success/
failure receipts, independent verification, complete measured resource and file
inventory, small real baseline and intentionally failed members, focused
coverage/regressions and final ticket checks. No runtime or science correction
is authorized by this ticket.

## Rollback

### Authorized correctness repair — 2026-09-14

The user approved the bounded repair with “fix it”. DIA-003 revision 1.2 and the
existing plan record the additional runtime/test scope. The previous no-runtime-
change statements describe the preserved pre-repair evidence, not the new build.

The repair passes real protein S4 objects to all five protein QC tracking outputs,
recognizes S4 with isS4, counts the declared design samples, and uses the declared
protein key for duplicate discovery and displayed counts. Plain-table helper APIs
remain unchanged. Existing S4 slots, direct slot access and legacy global tracking
ownership are preserved under the conventions' existing-code compatibility rule;
this is not an architectural migration.

Peptide validity permits absent precursor q-values only with the existing rollup
aggregation provenance, matching protein key, declared peptide raw quantity and
required quantity/sample columns. Raw objects still require q-values. No quantities,
confidence values, filters, aggregation methods, thresholds or S4 slots changed.

First focused run (session 78768) failed: new test used an unsupported withr
argument; existing handoff assertions expected tables; fake S4 fixtures omitted
declared key slots (including the previously recorded IQ fake-input defect).
It reached 18 new-test and 528 module assertions before the failure limit stopped
execution. Those errors were corrected in tests, not bypassed in runtime code.
Second run (84266): 772 assertions passed, zero failures/warnings, one explicit
disk-backed-service-only skip. Includes actual IQ fixture execution, rolled/imputed
validity, real protein tracking, and full-history RDS/revert checks. Additional
duplicate-removal and rollup regression checks (3004) passed 153 assertions, zero
test failures/warnings, with two pre-existing explicit skips: missing LFS cp03 and
an inline-wrapper timeout branch outside the extracted wrapper. A sandbox
timedatectl startup warning also occurred outside the test reporter. None of the
three skipped requirements is claimed verified.

The allowlisted install-corrected runner operation archives HEAD and applies only
the runtime R/DESCRIPTION/NAMESPACE diff into a fresh external source directory;
it records the parent commit and patch digest and installs into a new library.
The original library is only structurally checked, never overwritten. Installation
is preparation, not a timed performance member. Corrected installed verification
and full baseline qualification remain outstanding.

Corrected installation completed at browser-smoke.LYcZBQ (session 44092), parent
c7c12851ea4d9e96f91df0e1f9e7b91c2eb51017 plus runtime patch SHA-256
8dc098a17fa5be9d0c685a7b9893323d48863a3a8bd3df72b7414b7564083ca6.
Service runtime 118.221 seconds, CPU 104.233 seconds, reported peak 1G, zero swap.
R CMD INSTALL passed temporary/final load checks. Installed artifact SHA-256:

- DESCRIPTION: 8fc556281f7e5a34073b05f7661437f42c32c855eb0f3657d4c1656e257e76c7
- NAMESPACE: 9dc606a759229234e12f5a0b55879af5b7230a42f1177a4d86fa1e088879ccc5
- R/MultiScholaR.rdb: b1540f3d770584c9468cba4090ad7681a15f42cb2cb5e7fca67b0f31bd5fccc3
- R/MultiScholaR.rdx: ac6896c6aa0e0f35d120b4e6c11c0cb48ebbb966bc5af5986405a2e05e06f6e1

The same production browser-proteins operation, now using that isolated library
and unchanged IE9ohP inputs, passed_qc_only at browser-smoke.iBGaVY (28812).
All seven peptide and five protein states completed, including IQ rollup,
accession cleanup, intensity filtering, duplicate removal and replicate filtering.
Runtime 142.565 seconds, CPU 256.384 seconds, peak 2725732352 bytes, zero swap.
These are small diagnostic costs, not paired performance evidence: FAORgA aborted
earlier, so its elapsed time is not a comparable completion denominator.

cmp verified byte-identical IQ output between original FAORgA and corrected iBGaVY.
Output SHA-256: dcbc11e71e83225b1de6130b3badb3fc9f2cb86754d50c984345256cfa8da520.
This establishes agreement for that IQ artifact only. Normalization, DA, complete
automatic-output/session boundaries and G0/G4 remain unqualified. The receipt's
historical kind label says design_smoke; explicit status and state lists identify
its QC scope. load_verified=false describes structural preflight, not the
subsequent cold installed production launch.

A coverage diagnostic (55877) failed before tests because function_coverage does
not accept a vector of function names. A corrected environment-scoped diagnostic
was run; no coverage percentage is inferred from passing assertions.

Coverage attempts retained: 87635 failed when a standalone installed test
environment lacked an internal helper; 10233 passed the installed tests but
reported NaN coverage because installed functions lacked usable source references
(artifact /tmp/dia003-corrected-function-coverage.rds). Source-file instrumentation
22445 then hit a logger namespace error in its artificial source environment.
Neither NaN nor a partial failed run counts as coverage acceptance. A later
read-only inspection initially omitted loading covr's S3 methods and was corrected.

Source-namespace diagnostic 26824 completed using environment_coverage after
devtools::load_all (coverage only, never timing). Its selected-function line
coverage is 64.81481%, below the 80% convention target: the broad legacy plotting
function is 35.44304% and class validity is 58.97436%. All 30 changed source lines
mapped to instrumented expressions were executed, including the new validity
predicate and missing-sample error. That changed-line result is not a substitute
for whole-function coverage. A fake-class package-name warning occurred in this
standalone test environment. Artifact: /tmp/dia003-corrected-source-coverage.rds.

The final caller audit additionally found RUV/correlation consumers discarding
the S4 metadata, and hard-coded post-RUV/final counts. Their identical handoff/key
correction was explicitly added to DIA-003's scope before editing. No RUV fit,
correlation calculation or population changed. Focused source run 24231 passed
580 assertions, zero failures/warnings/skips, including real group-key objects
and correlation Apply/Skip tracking. LYcZBQ/iBGaVY remain the earlier QC-only
build evidence; the final expanded repair requires its own installed identity.

Final repair installation e3lvCL (53085) succeeded, with parent c7c12851 and patch
SHA-256 2ec3b2790e81818e2929e8622d97cb0b716fdd62f147a20955c7dd237676b4ff.
Installation runtime 76.617 seconds, CPU 73.255 seconds, reported peak 1G, zero
swap. No install-time difference is interpreted as workflow performance evidence.
Installed artifact SHA-256:

- DESCRIPTION: 450ae041c19463ff7f19ac8e289bbd63ade6b08fea1c8c47d186d372d584d08d
- NAMESPACE: 9dc606a759229234e12f5a0b55879af5b7230a42f1177a4d86fa1e088879ccc5
- R/MultiScholaR.rdb: e48edae6a2aa3092f1911523862c8ec250c098df063ccbb4247b3855fa65c104
- R/MultiScholaR.rdx: 75e048e3d5736151b430528833269ba7e5f4a1f405a4bcdc4d35741fef98a2f2

Final installed production verification p4sE8I (50897) passed_qc_only with the
unchanged IE9ohP inputs and e3lvCL library. All seven peptide and five protein
states completed. Runtime 121.557 seconds, CPU 219.764 seconds, cgroup peak
2711982080 bytes, zero swap. cmp again confirmed byte-identical IQ output versus
B0-original's FAORgA artifact. The limits remained 8 GiB/zero swap/900 seconds.
This closes the two authorized runtime defect repairs for the tested QC path;
the RUV/correlation consumer corrections have source regression evidence, not a
complete installed normalization journey. Full DIA completion, broader coverage
target and performance qualification remain DIA-003/004 work, not passed gates.

Final checks: ticket schema/status consistency, bash -n and git diff --check pass.
No package dependency, public signature, S4 slot, scientific estimator or threshold
changed. The reusable runner was used for both separate installations and browser
checks; the untouched original installation, earlier corrected build and every
failed diagnostic remain available.

Before the authorized amendment, only harness/test/documentation files changed.
The amendment adds the named runtime compatibility repairs and their tests. Keep
the pinned source, isolated library, failed diagnostic and all earlier fixture/OOM evidence.
Do not reset the branch or remove scientific inputs. No saved-entity runtime
adapter exists to disable or migrate.

## DIA-003 resumed after standalone DIA-005/006

The preceding rollback text describes the earlier DIA-003 surface; standalone
saved-entity helpers now exist under completed DIA-005/006, but no GUI uses them.
G0 and all performance gates remain unaccepted.

Current driver audit confirms the peptide/protein smoke path uses explicitly
relaxed tiny-fixture settings, including imputation missingness allowance 1 and
one protein-replicate worker. Those runs cannot be relabeled out-of-the-box QC.
Both action functions now accept a strict logical `diagnostic_settings` argument.
TRUE preserves the existing smoke behavior; FALSE sends no QC filter, missingness
or CPU override and retains only explicit IQ route selection. The UI actions and
completion checks are unchanged. The full baseline CLI is still held, so this
does not silently enable an unqualified timed path.

Harness run 75170 (`devtools::test(filter = "workflow-checkpoint-harness",
stop_on_failure = TRUE)`) passed 196 assertions, zero failures/warnings/skips,
in 5.5 seconds. New driver doubles verify both settings paths, identical step
sequences and rejection of invalid mode flags; they are not browser/science proof.
No runtime package/scientific code changed in this continuation.

A read-only terminal-accounting check of the prior iiJUC7 service returned
LoadState=not-found, zero start/exit timestamps and unset MemoryPeak/CPUUsageNSec.
Its captured log remains useful, but its rounded human-readable peak is not an
exact combined-memory metric. Full-boundary accounting must be retained before
systemd unloads a measured unit; absent values are not zeros. This confirms the
existing timing hold rather than supplying a qualified comparator.

## Exact terminal accounting probe

Added two fixed allowlisted runner operations, `--service accounting-probe` and
`--service accounting-probe-failure`, taking the existing LIBRARY OUTPUT_ROOT
arguments. They touch a bounded 64 MiB raw allocation and exit 0 or deliberately
7; they do not execute omics. Existing browser/install operations remain unchanged.
Both probes retain the 8 GiB/zero-swap/900-second service limits and serial lease.

Local systemd help suggested structured output. Actual probes 7Mx51N and 4hckr7
showed `--wait --json=short` returns only the unit name; its accompanying human
summary is rounded and not a metric source. Those observations are preserved.
The user service journal does retain exact CPU_USAGE_NSEC, MEMORY_PEAK,
MEMORY_SWAP_PEAK and IO_METRIC_* fields after service exit. The probe now captures
those selected journal fields into terminal-journal.jsonl, then decodes a uniquely
matched unit/invocation/accounting-message record into terminal-metrics.json.
The original service log is written by the service, not copied by an outer tee.
The JSON unit response and all controller/journal/decoder exit codes remain separate.

Decoder tests require exact unit and invocation identity, reject ambiguous records,
preserve absent counters as NULL, and reject systemd's UINT64_MAX sentinel rather
than reporting it as a measured value. Counters above R's exact integer range are
also unavailable. Elapsed remains explicitly NULL and promotion_authorized FALSE;
no timestamps from journal-delivery events are substituted for service elapsed.

Real updated success probe aosPic recorded memory_peak_bytes=397488128,
cpu_usage_nsec=509259000, swap_peak_bytes=0, and service/journal/decoder exits 0.
The follow-up N9o2EB is the intentionally failing probe, with its own terminal
record and service exit 7. Both reside under the established external proof root
MultiScholaR-dia-proof-20260911-VhC1n5. These are infrastructure probes, never
performance members or evidence of a memory-saving treatment.

Harness run 81912 passed 206 assertions with zero failures/warnings/skips in
5.9 seconds. Bash syntax and git diff checks pass. Remaining qualification includes
an exact elapsed boundary, controller-resource inclusion/exclusion, complete
baseline dispatch and production-driver evidence. Default measurement execution
remains held; journal-backed counters alone do not qualify G0.

## Retained-service elapsed boundary

The fixed `accounting-probe-retained` operation now retains its completed service
only long enough to capture systemctl properties, then explicitly stops that
request-owned unit. Observation is bounded; timeout retains the unit for inspection
and does not trigger a replacement job. Existing browser operations are untouched.

Initial probe F7yIXJ exposed exact main-process start/exit timestamps and cached
MemoryPeak/CPUUsageNSec. Its main interval was 0.569428 seconds, but TasksCurrent
was unavailable. That was not credited as full-tree completion. Probe 0hn7iT
additionally confirmed cached MemorySwapPeak and I/O properties after process exit.

The local systemd.service documentation states that ExitType=cgroup waits for the
whole group, unlike the default main-process exit policy. The retained probe now
uses ExitType=cgroup and deliberately launches a one-second child that outlives
the R main process. JMaqsm recorded main exit after 0.523762 seconds but the
terminal state change after 1.516910 seconds. pz26Fr repeated this with a captured
interval receipt: raw start 2232579120, main exit 2233119297, terminal state change
2234112481 microseconds; main interval 0.540177 seconds versus tree 1.533361.
Its cached peak was 397963264 bytes, CPU 522830000 ns, swap 0. Decoder and explicit
unit stop both exited zero. These are synthetic infrastructure probes only.

The reader binds unit/invocation, validates monotonic values, and keeps the
main-process interval separate. It exposes a tree interval only for ExitType=cgroup,
terminal exited state, an empty released ControlGroup and a terminal timestamp
no earlier than main exit. Missing evidence remains unqualified. It never converts
an unavailable TasksCurrent into zero. Promotion remains false because this does
not establish GUI completion or outer-controller resource scope.

Harness run 52150 passed 218 assertions. A subsequent serialization check found
jsonlite's default decimal precision would round the derived interval, so the
interval writer now uses digits=NA. The pz26Fr JSON is prior to that adjustment;
its unrounded source timestamps remain in terminal-properties.txt. Tests verify
full-precision persistence and refusal to overwrite interval receipts; terminal
metric receipts likewise refuse overwrite. No historical receipt was rewritten.

Final harness run 29503 passed 220 assertions, zero failures/warnings/skips, in
5.9 seconds. Bash syntax and git diff checks pass. The retained probe is not yet
the production baseline launcher; full driver dispatch and controller accounting
remain DIA-003 work before DIA-004 timing can be qualified.

## Real DIA diagnostic with retained process-tree accounting

Added allowlisted `--service browser-da-accounted LIBRARY OUTPUT_ROOT REPORT FASTA
DESIGN`. It dispatches the unchanged internal browser-da diagnostic, but uses the
retained ExitType=cgroup launcher and service-owned logging. Browser environment
flags remain NOT_CRAN=true and the detected Chrome executable. The observer waits
up to 930 one-second checks around the enforced 900-second service limit; timeout
does not restart the job. Existing browser-da and other service operations remain
available. This is not baseline-mode promotion or out-of-the-box QC selection.

Real member 5Olm86 used installed iiJUC7 and the existing IE9ohP small report,
FASTA and design. The full diagnostic browser path completed with status
passed_da_ui_only. It exercised seven peptide steps, five protein steps,
cyclicloess/RUV Skip/correlation Skip, session export/reload, DA, 64 interactive
volcano points, 12 expression points after table selection, and heatmap export.
Retained terminal properties/interval/controller receipts are under that member
directory in the established external proof root. No source/test-mode app was used.

Observed service costs (not performance acceptance):

- Full process-tree interval: 216.633534 seconds; main-process interval 216.310568.
- MemoryPeak: 3388731392 bytes; MemorySwapPeak: 0.
- CPUUsageNSec: 424466085000.
- IOReadBytes: 31498240; IOWriteBytes: 65454080 (service-charged I/O, not logical bytes).
- Result success, MainPID=0, SubState=exited, ExitType=cgroup, released ControlGroup.
- Interval decoder and explicit unit stop exited zero.

These values include installed process/package startup, browser/driver/descendants
and in-service outputs through tree exit. They do not include the outer controller
or the separate verifier. Retained heap and complete phase attribution are still
unavailable; do not fabricate them or add RSS to the cgroup peak. No test workload
or verifier was run concurrently with this browser service. A sandbox-denied
read-only systemctl observation was retried with the existing approved permission;
the running service was never restarted.

Separate verifier oZYcgP passed_small_da_export_checks, confirming 16 history
records, 64 proteins/12 samples, 33 DA reference columns, and seven heatmap rows
with four clusters. Normalized TSV SHA-256 remains
b698cd7dac2c2c2003c2874bcff4375fbb75bb125e386c1aac7704497348bf4a;
DA TSV remains 0bdc7e2e5b25cf9425a5533a56eb283bd09da04f5b5b4d5972144b25a9eab140.
The same-implementation replay is not an independent statistical oracle. The
verifier was outside the browser interval and exited zero.

Harness source bindings during that browser run:

- shell: afa3259f42b58d699600b986522627e6fa722b1f4f8fe4d9f19ef8eff7778f3b
- R runner: 00392b9d24ad07787d2b7571aabeafc27062d91b76e39c98e9156a9b050846d8
- browser driver: 8504de3f33f9dbf5a66c8b904dab58e63a286fe85b9f60fb809cbec24b9d5ba3

After the run and verification ended, the controller-request reader was extended
to recognize the newly allowlisted operations. Shell doubles now verify accounted
browser dispatch, terminal record retention, success/failure reporting and stopping
only the exact request-owned unit. Final harness run 43148 passed 234 assertions,
zero failures/warnings/skips, in 6.5 seconds. These doubles do not qualify resources.

The tiny relaxed settings, single diagnostic member, incomplete out-of-the-box
baseline dispatch, controller attribution and unqualified representative workloads
prevent G0/performance acceptance. No elapsed/memory gain is claimed from this run.

## Representative inputs and first default-QC control diagnostic

Added fixed `prepare-companions LIBRARY OUTPUT_ROOT RECIPE EXPECTED_SHA256` service
operation. It verifies an existing frozen report fingerprint and derives matching
design/FASTA companions directly from the existing fixture generator. It does not
hydrate a precursor table or S4. Tests compare every small-fixture protein sequence
with the original report generator and exercise 10000-protein/48-sample companion
generation with precursor/S4 constructors replaced by failing guards.

Preparation tJ6U1c used the corrected control build VEqJKP, not the later saved-entity
candidate iiJUC7. It verified the 3840000-row representative report SHA-256
2414125be92a11fbe4342fdf354e0214915e30592bec66f3e9d6a984f9290841 and wrote:

- design.tsv: 4fbd08901905c1339af58ab619d75fade53655484153d627bd7cd39a45bfe6e5
- proteins.fasta: 352b56704c330f3e59646101df2c52a469e5ebe17da652fdd0d3ea85e672f6d4
- generator: 29ba69bc649f5f273e1bf02a9ac3ccd7c3eb6cddae7b21a8df42ef8ad04ab47e

This was input preparation, outside any benchmark. The report was recently read
for fingerprinting; no cold filesystem-cache condition is asserted. Its recipe's
control_qualified remains false.

Added explicit `browser-da-default-qc` dispatch to leave QC filter/missingness/CPU
controls at their actual GUI defaults. The later route remains explicitly fixed
cyclicloess/RUV Skip/correlation Skip; this is not coverage of all normalization or
RUV choices. Original tiny diagnostics remain unchanged. Dispatch tests prevent
the new path from silently falling back to relaxed settings. Harness 97059 passed
243 assertions with zero failures/warnings/skips.

First representative control diagnostic gYW1ha failed before QC while waiting
for Shiny idle during Design Save. The service exited 1, not OOM or service timeout:
MemoryPeak=4494426112 bytes, swap=0, CPUUsageNSec=237842883000. It cannot be a finite
control denominator. Its receipt/logs/terminal properties remain preserved.
Browser logs end with defined_contrasts_display, an idle acknowledgement, then a
busy event and a 30000ms wait; the app log shows configuration/worker setup and no
reported scientific-method error. The driver invoked a generic 30-second idle
wait inside click(save_results) before reaching its declared 180-second Design
completion wait. That is a premature harness limit on the larger Save operation.

The bounded correction applies the existing 180-second budget only to Design Save,
leaving ordinary control-edit waits at 30 seconds and the service limit at 900.
No runtime package, scientific setting or input changed. Harness 38833 passed
245 assertions; the new test verifies this exact timeout routing. Replacement
diagnostic zgi11N was then launched once on the same control build and input set.
Its outcome must be recorded separately; gYW1ha is never relabeled a successful run.

Replacement zgi11N passed Design Save, selected all 48 declared samples in their
two groups, and entered QC with the explicit GUI-default profile. It then failed
the driver's 120000ms JavaScript wait for qvalue_results to contain qvalue_filtered.
Its logs end with Shiny busy; no result/idle acknowledgement is recorded. Thus
the Design-wait correction worked, but the route still did not complete.
MemoryPeak=5304492032 bytes, swap=0, CPUUsageNSec=523498149000; service Result=exit-code,
ExecMainStatus=1, not the 900-second service timeout or an OOM result. Preserve
this second failure and do not treat it as a finite comparator or scientific pass.

The next driver correction must address the inconsistent scientific-action wait
caps coherently within the fixed service budget, rather than infer a method's
completion time from a prematurely stopped callback. No third representative
attempt was launched in this continuation. The representative control remains
unqualified. Runtime/science and fixture inputs were not changed.

Bindings for gYW1ha: shell
6b63c7ceaf73f39e388f4fabdec25a4a18356d85a3dd042764352c420e32796f;
R runner fd0e70db2317acda4adb5e8bda4798687c563af48e01af4445a0ff11bcefeb07;
driver 8504de3f33f9dbf5a66c8b904dab58e63a286fe85b9f60fb809cbec24b9d5ba3.
Only the driver changed for zgi11N, to
6e6c6dbdd49af278d110c111392431a4a254ba56e3b5ebbb9e9b62e9a622a2fd.
All member/companion/failed receipts remain under the established external proof
root. The general baseline CLI and G0 qualification remain incomplete.

## One scientific wait budget

Replaced the inconsistent scientific waits with checkpointWaitBudget: remaining
milliseconds from the driver's process elapsed time within the fixed 900-second
service lifetime. Scientific callbacks, idle waits, phase completion, plotting
and export waits share that remaining allowance; they do not receive a fresh
900 seconds at each step. Missing DOM/control readiness retains a 30-second cap,
also clamped to the remaining allowance. Exhaustion raises an explicit error.
Systemd RuntimeMaxSec=900 remains the authoritative hard limit, including any
R initialization/clock-origin difference. No resource or scientific threshold changed.

Tests use deterministic clocks rather than depend on the age of a package-test
process. They verify 120 elapsed seconds leaves 780000ms for scientific waits,
short control waits, final QC completion using the same policy, a one-second
remaining budget, invalid clocks and exhaustion. Existing diagnostic/default-QC
selection, UI action and failure checks remain passing. Source run 47434 passed
251 assertions with zero failures/warnings/skips in 5.6 seconds; git diff check passed.
Driver SHA-256: 07b6250d4e61e84d24034d255b2bcfdcc786bc2954e0809fb31a02ff3a7d5c71.

The earlier gYW1ha/zgi11N failures remain unchanged and unqualified. This code
verification does not establish the representative journey's completion time,
authorize favorable performance reruns or turn either failed attempt into a
denominator. A new-driver execution is still needed for qualification; none was
launched during this wait-policy implementation turn.

## Representative control: genuine memory-limit failure

Driver-qualification member bQWqlM used the corrected pre-storage control VEqJKP,
the same frozen representative report and tJ6U1c companions, default QC and the
declared cyclicloess/Skip route. Its driver hash was
07b6250d4e61e84d24034d255b2bcfdcc786bc2954e0809fb31a02ff3a7d5c71.
No code/settings were changed or test/verifier workloads run alongside it.

Unlike gYW1ha/zgi11N, this attempt reached an actual service resource failure:

- Result=oom-kill; ExecMainCode=2, ExecMainStatus=15.
- MemoryPeak=8589934592 bytes (the enforced 8 GiB limit); MemorySwapPeak=0.
- CPUUsageNSec=1581538348000.
- Main start 5184549428, main exit 5937973985, terminal state change 5938338993
  microseconds on the systemd monotonic clock (about 753.8 seconds to failed state).
- MainPID=0, released ControlGroup, and explicit unit stop succeeded.

This is capacity evidence, never a finite completed-journey denominator or a
performance gain/regression claim. The pre-storage VEqJKP source does not contain
utils_workflow_saved_entities.R. The new save/hydrate helpers were not used by
this control. No memory limit, workload size, format or scientific threshold was
relaxed after this result.

The in-app smoke receipt remains `running`, because the process tree was killed
before its exit handler finalized it. There is no browser-logs.json for this member.
The terminal systemd record overrides that unfinished application status; it does
not indicate a still-running job. Exact failing call/stage attribution is not
available from this run and must not be inferred from the earlier Q-value timeout.
All surviving files, request and terminal receipts remain preserved.

The retained-interval decoder now also exposes service_result and an explicit
service_succeeded flag. Missing results and OOM records cannot become success;
failed services retain NULL completed-tree elapsed for acceptance. Tests cover
the observed OOM/signal combination. Harness run 64880 passed 256 assertions,
zero failures/warnings/skips, in 5.6 seconds. The older bQWqlM interval receipt
was not rewritten; its raw properties contain the authoritative OOM result.

Read-only source inspection found two profiling candidates, not established causes:
full-table row-text construction/sorting in .peptideQcCanonicalise, and per-row
identifier combination/tokenization in classifyPeptideBiologicalExclusions.
The inspected files match the installed control's preserved source exactly:
audit helper SHA-256 0ccc3e8a3f7686d53f1aca72d96d16e91e2c46d1d3a7b67d3dc172d0fcccf1ae;
exclusion helper df0227081d7fd5d519362db7ed62325350a8a82a8d86641a0580640ccd526ec2.
Neither was edited. Separate diagnostic profiling is the next available in-scope
step; runtime optimization would need its own explicit ticket scope. G0 and all
dependent performance acceptance remain held, not completed by recording the failure.

## Diagnostic application profiling and lifecycle correction

Added explicit `browser-da-profile` dispatch: same control/default-QC route and
limits, with application-R CPU sampling every 20ms plus GC/heap profiling. The
initial receipt marks diagnostics=true and timing_eligible=false, and explicitly
excludes driver/Chrome/parallel-worker stacks from the profile scope. Default
nonprofiled launches remain unchanged. Raw profiles are request-owned and cannot
be overwritten. The host reports profmem capability available.

Initial plumbing tests passed 269 assertions (78571), but the process also emitted
a sandbox-related timedatectl warning outside the test assertion stream. Diagnostic
member gUjqnq then OOM-killed at 8 GiB, zero swap, with CPUUsageNSec=1565255084000.
This profiled attempt is never performance evidence. An attempted live process
snapshot was unavailable in the tool PID namespace and subsequently found the
driver already gone; no per-process RSS/CPU attribution is inferred from it.

The profile survived, but analysis proved it unsuitable for workflow attribution:
435 lines / 241597 bytes, SHA-256
7639d4efaeb20ea2fd7cc76b320d7971037c39dfbb3b9f7fdd770d14336dd247,
only 434 samples (8.68 seconds). Nearly all inclusive samples were loadNamespace;
no Q-value/classification/audit-stage samples were captured. GC's 2.78 seconds in
that short construction profile cannot explain the later OOM. The instrumentation
incorrectly stopped when the app factory returned its shiny.appobj, before AppDriver
started the real event loop. This was a harness bug, not a scientific finding.

The corrected wrapper keeps sampling after construction and registers profiler
shutdown with shiny::onStop(session=NULL). Construction errors still stop profiling.
A real isolated Shiny event-loop test schedules CPU work after construction and
checks that its named frame appears in the profile. Its negative control explicitly
stops profiling at construction return and must fail to capture that frame. Both
apps execute the same work and return the same result. This tests the missing
lifecycle boundary, not merely the existence of a profile file.

Final harness run 18950 passed 270 assertions, zero failures/warnings/skips, in
10.9 seconds; git diff checks passed. Corrected driver SHA-256:
968e9eed5c74f7775cb8c41a0836c3b4bc969efd3de7b3607246d7d116ea5028.
No replacement long profile has run yet. gUjqnq and the unprofiled bQWqlM OOM are
preserved separately; no workflow cost attribution or runtime optimization is
authorized by the invalid construction-only profile.

Corrected replacement profile 5cCzlU is now running on the same VEqJKP control,
representative report and tJ6U1c companions, with driver
968e9eed5c74f7775cb8c41a0836c3b4bc969efd3de7b3607246d7d116ea5028.
An early read found 5168 profile lines / 3716977 bytes and active serviceApp/input
processing frames, confirming sampling continues beyond construction. The latest
service observation was active/running, MainPID=136563. This is a live diagnostic,
not a completed outcome; recheck its service/terminal artifacts before interpreting
the profile or starting anything else. No concurrent test or verification workload
was launched, and no cause or performance gain is inferred from the partial samples.

## Corrected workflow profile: measured identifier-classification hotspot

5cCzlU is now terminal: Result=timeout, ExecMainCode=2, ExecMainStatus=15,
MemoryPeak=6060912640 bytes, swap=0, CPUUsageNSec=1947757958000. It reached the
fixed 900-second service limit. These profiled costs are not a comparator and
must not be compared with bQWqlM's unprofiled OOM as an improvement.

The retained application.Rprof is 81435710 bytes / 42704 physical lines, SHA-256
237a643580d2f26f444baed1e7454f269227b92ab1badb934697e9357360d78e.
Unlike gUjqnq, it contains workflow execution. Standard summaryRprof(memory="both")
failed on multiline Shiny labels. memory="none" treated label continuations as
samples, yielding an inflated 854.06 seconds and a spurious empty-name bucket;
those percentages are not used here.

The added bounded-purpose reader recognizes the numeric four-field memory prefix
of genuine samples and excludes label continuation lines from the sample count.
It searches only sample-start lines for the selected exact frame names, so frame
counts are conservative if a frame occurs after a multiline label. Counts are
inclusive/nested, not independent contributions and must never be added.
There are 39939 true sample prefixes and 2764 continuation lines: 798.78 sampled
application CPU seconds at 20ms intervals. This is not whole-service CPU or heap
peak attribution. The original raw profile was not rewritten.

Measured lower-bound frame membership:

| Frame | Samples | Sampled seconds | Share of captured samples |
| --- | ---: | ---: | ---: |
| runPeptideQvalueApplyObserver | 34114 | 682.28 | 85.42% |
| classifyPeptideBiologicalExclusions | 33909 | 678.18 | 84.90% |
| .tokenizeProteinAccessions | 25040 | 500.80 | 62.70% |
| .normaliseProteinAccession | 17045 | 340.90 | 42.68% |
| trimws | 26178 | 523.56 | 65.54% |
| GC | 3115 | 62.30 | 7.80% |

No .peptideQcCanonicalise or .savePeptideQcState frame was observed by this reader.
The earlier audit-sort hypothesis is therefore not supported as the first fix by
this incomplete workflow profile. CPU attribution points first to identifier
combination, tokenization and normalization inside biological-exclusion classification.
The code applies that work to precursor rows even when accession combinations
repeat; the frozen input has 3840000 rows and 10000 generated protein identities.
This does not by itself prove the exact allocation responsible for bQWqlM's OOM.

The first frame-reader attempt rejected valid names because of regex escaping;
diagnostic 19949 and test run 55652 failed and are retained in this record. The
corrected reader explicitly uses PCRE for forbidden-character checks. Its fixture
contains multiline labels, including a false frame-name mention on a continuation
line, and verifies only the real sample records contribute to counts. Final harness
run 7654 passed 275 assertions, zero failures/warnings/skips, in 9.8 seconds.
git diff checks passed. No new profile run was performed to obtain a favorable result.

Recommended bounded next scope: optimize classifyPeptideBiologicalExclusions in
R/func_prot_qc_exclusion_helpers.R by evaluating unique identifier combinations
once and mapping results back to original rows. Preserve row-level explicit flags,
mixed-accession/isoform semantics, order, attributes, manifest provenance and all
classification outputs; add equivalence tests in test-prot-02t-qc-peptide-exclusions.R
and isolated old/new performance evidence before GUI replay. This is a runtime
optimization outside DIA-003's current allowlist and stop conditions, so it is a
proposal awaiting user authorization, not an implemented or accepted change.
Keep the existing tickets and all original controls/failures; do not create a
competing repair series or claim a finite speedup against an OOM/timeout control.

## Authorized identifier reuse implementation and isolated case freeze

The user approved the scope extension. DIA-003 revision 1.6 adds the exclusion
helper/tests and narrowly permits this runtime optimization. Identifier tags are
computed on distinct text values per provenance column and mapped back with OR;
the semicolon column boundary cannot create a cross-column tag. Non-text inputs
retain the historical matrix-coercion path. Repeated accession text is tokenized
once per distinct value per call, then token counts/matches are mapped to all rows.
Explicit row flags, selection rules, provenance and output construction are unchanged.
No cross-action cache, worker, threshold, format or audit-sort change was added.

A frozen copy of the original classifier is a source-only test oracle in the
existing exclusion test file; parsing/manifest helpers remain unchanged. Tests
cover exact full results and input immutability for repeated rows with different
explicit flags, factors, named values, numeric-group fallback, missing/blank and
mixed accessions, tibble attributes, all exclusion switches and file/vector/empty
manifests. A cardinality check requires only three tokenizations for 1000 rows
and repeats the work on the next action rather than retaining a hidden cache.

Initial test 6380 exposed the reference's existing zero-row non-character error;
the test now verifies that same error in both versions, without adding an unrelated
behavior fix. Source identification/exclusion run 27141 passed 129 assertions,
zero failures/warnings/skips. Candidate L04HyQ installed successfully, runtime
patch SHA-256 e54eb231213a47f7779e990cf33fe818fad6144e5d2bdb4620431ebd5962321b.
Installed check 33345 passed its result comparisons but failed the cardinality
test's implicit pkgload namespace lookup. The test now explicitly names MultiScholaR;
no runtime change was needed for that fixture error.

Before measurements, freeze isolated case representative_prefix_384000_v1: first
384000 rows of the same fingerprinted representative report, all columns, 1000
protein groups, non-ALTREP vroom read, empty contaminant manifest and both exclusion
switches TRUE. No fixture generation inside measurement. This is a bounded
classification-component experiment, not the full representative workload or G0/G4.
The existing runner's exclusion-benchmark operation takes LIBRARY OUTPUT_ROOT REPORT,
uses fresh installed execution under the same service limits, and includes startup,
namespace load, full report integrity check, prefix read, classification, native
uncompressed RDS output/fingerprint and tree exit. Classifier-only elapsed/CPU is
recorded separately, without forcing GC. Original VEqJKP and candidate L04HyQ run
serially; complete output equivalence must be checked independently. A single pair
will be descriptive evidence only, not an uncertainty-qualified performance gate.

## Classifier comparison decision — 2026-09-16

The bounded comparison is complete; GUI replay is **held**, not accepted.
The current classifier source and installed candidate L04HyQ source both hash to
f3805243f362ad712c85b2743242501b2b214159b166c7903c8ea30e35172a98.
The original G2aAKR terminal receipt confirms success. Candidate bpgCML then ran
once on the same fixed case through the existing runner, with no concurrent test
or verifier workload. It exited successfully, with zero swap. No runtime or harness
code changed during this comparison work item.

| Observation | Original G2aAKR | Candidate bpgCML |
| --- | ---: | ---: |
| Classifier elapsed, seconds | 56.591 | 1.335 |
| Service-tree elapsed, seconds | 81.004374 | 12.232367 |
| Service CPU, nanoseconds | 78233686000 | 12072928000 |
| Service peak, bytes | 1407930368 | 2010849280 |
| Physical read bytes | 0 | 416178176 |

The observed reductions are 97.64% classifier elapsed and 84.90% service elapsed,
but the candidate's recorded service peak is 42.82% higher. These are descriptive
observations from one pair, not matched-protocol confidence estimates or a GUI claim.
The earlier control and current candidate ran at different times; comparable host
load and input-cache conditions have not been established. The physical-read
difference is evidence of an I/O-condition mismatch, not a measurement of resident
file-cache bytes. It neither proves nor excuses the memory increase.

Independent verification 31436 recomputed both output SHA-256 values, decoded both
complete result objects and confirmed identical(a, b), including all 384000 raw
rows, classified/analysis tables, ledger, summary, manifest and attributes. Both
236550286-byte RDS outputs have SHA-256
c71e37be039810ba395e7cd6100f73d48d5d8a057f9030fc9f1c92212456d124.
The earlier installed test correction subsequently passed all 112 exclusion
assertions (70712); the pkgload-only mock failure remains recorded above.

Decision: exact-output evidence and an observed isolated time improvement exist;
a memory benefit/non-regression and complete-workflow benefit are not established.
Do not launch GUI replay or promote the optimization as an overall performance
pass. Any follow-up must be a separately bounded memory-comparability/attribution
item with its cache condition, measurements and run budget specified in advance.
Do not rerun simply to seek a lower peak. Existing OOM/timeout controls and both
members of this comparison remain intact. DIA-003 and all performance gates stay open.

## Approved memory diagnostic pair — completed

The user approved one follow-up pair. The case was frozen in DIA-003 revision 1.8:
same VEqJKP/L04HyQ classifier builds and 384000-row input/output boundary, with
fincore residency checks and existing cgroup snapshots. No classifier edits,
forced GC, profiling, host-cache flush or external input warming. Original GMJtyN
ran before candidate 8wKRwf, serially, with no concurrent owned work. Both checks
(entry and pre-classifier) reported the 556028391-byte report fully resident
(556032000 page-aligned resident bytes). Library/native-code cache and unrelated
host activity were not controlled; this is not a fully matched G4 experiment.

Preflight 51131 found a byte-count integer/double assertion mismatch and a sandbox
port denial in the existing Shiny lifecycle test. Byte counts were normalized to
double, and the same suite was rerun with the existing approved Rscript permission.
Run 23935 passed 282 assertions, zero failures/warnings/skips, in 5.5 seconds.
The runtime classifier hash remained f3805243f362ad712c85b2743242501b2b214159b166c7903c8ea30e35172a98.

| Observation | Original GMJtyN | Candidate 8wKRwf |
| --- | ---: | ---: |
| Classifier elapsed, seconds | 22.210 | 1.512 |
| Service-tree elapsed, seconds | 32.304747 | 16.173796 |
| Service peak, bytes | 1413488640 | 1353015296 |
| Post-classifier anonymous bytes | 1166163968 | 1111691264 |
| Post-classifier file-cache bytes | 5062656 | 8192 |
| Post-output anonymous bytes | 1166258176 | 1111785472 |
| Post-output file-cache bytes | 241614848 | 236560384 |

Both members exited successfully with zero swap. The observed candidate changes
are -93.19% classifier elapsed, -49.93% service elapsed and -4.28% service peak.
The largest recorded cumulative peak in each run followed output writing; the
post-classifier anonymous snapshot was 54472704 bytes lower in the candidate.
Snapshots are not instantaneous peak decomposition or exact R-heap measurements.
Do not sum these overlapping/temporally distinct observations into a new peak.

Independent verification 1041 recomputed both file hashes and confirmed identical
complete decoded results. Both outputs also match the prior pair's SHA-256
c71e37be039810ba395e7cd6100f73d48d5d8a057f9030fc9f1c92212456d124.
Its subsequent generic phase summarization failed because R partial matching of
memory against memory_scope had left an extra scalar metadata entry under memory.
All six named snapshots are present and intact. A second summary attempt incorrectly
expected status="observed" instead of the reader's actual "available" schema.
Final analysis explicitly selected and validated the six named phases. No values
or raw receipts were rewritten and neither member was rerun.

Decision: the earlier higher peak was not reproduced with the declared resident-input
condition, and this pair supports a separately bounded GUI verification of the
approved classifier change. It does not prove that cache explains the entire prior
increase, guarantee memory non-regression across workloads, satisfy a 10% GUI-memory
target, or establish formal paired uncertainty. Keep the previous held comparison
and failed controls. No further optimization, measurement retry or GUI replay was
performed as part of this completed diagnostic item.

## Representative classifier-candidate GUI verification — failed

One separately bounded run on 2026-09-16 used existing operation
browser-da-default-qc, installed L04HyQ, the full representative/report.tsv and
tJ6U1c design/FASTA under the existing external proof root. Source classifier
SHA-256 f3805243f362ad712c85b2743242501b2b214159b166c7903c8ea30e35172a98
matched the installed build's archived source. Companion hashes matched the frozen
identities above. The exact invocation is retained in controller-request.nul.
No source/harness change, input shrinking, profiling, prewarming, competing owned
work or limit change was made during this run. Existing lifecycle-test evidence
was reused; this was one GUI attempt, not another measurement pair.

Member browser-smoke.zqsj4M, invocation e2febca1e5f446c5a640f081585b5038,
terminated with Result=oom-kill. terminal-properties.txt and terminal-interval.json
record MemoryPeak=8589934592 bytes, MemorySwapPeak=0, CPUUsageNSec=435638441000,
IOReadBytes=57663488 and IOWriteBytes=303161344. Main-process lifetime was
227.539417 seconds. Successful service-tree elapsed is correctly null; the time
to failure is neither completion time nor a finite speedup denominator. Controller
exit was 1, interval decoding and owned-service cleanup each exited 0. No live run
remains.

The project contains import/design artifacts (data_cln.tab, design_matrix.tab,
contrast_strings.tab and supporting FASTA/configuration files). No final session,
normalization/DA outputs or browser-logs.json are present. service.log is empty.
The application receipt still says running because the kill prevented finalization;
the terminal service record is authoritative. These files do not prove which
scientific operation was executing at the kill or decompose the memory peak.
No output verifier was run against this incomplete member.

Decision: full-size workflow verification failed. Keep this failure alongside the
original control OOM and both isolated comparison pairs. The classifier's isolated
exactness/time improvement is not a successful full-workflow capacity, memory or
performance result. No retry or adjacent optimization is authorized by this result.
DIA-003 remains open and DIA-004/G0/G4 are not promoted. Any further investigation
must declare its scope and budget without treating an unknown hotspot as diagnosed.

Read-only follow-up: the retained kernel journal narrows the memory failure without
another execution. `journalctl -k --since '2026-09-16 08:49:25' --until
'2026-09-16 08:49:45' --no-pager` records the exact zqsj4M cgroup at the OOM:
anonymous memory 8542744576 bytes, file memory 9392128 bytes, kernel memory
37789696 bytes. These are kernel failure-time observations, not exact live R-heap
sizes or an additive reconstruction of the previously measured peak. File memory
was not the dominant charge at this failure.

The killed R process was PID 100116, distinct from controller R PID 100072. Its
kernel record reports anon-rss 6823756 KiB and file-rss 585544 KiB. Eight other R
processes and the browser processes were also listed. Their RSS values are not
summed into the cgroup peak, and the journal alone does not establish their R-level
objects or the function allocating memory. Relevant literal journal lines:

```text
memory: usage 8388608kB, limit 8388608kB, failcnt 201
swap: usage 0kB, limit 0kB, failcnt 0
anon 8542744576
file 9392128
kernel 37789696
Memory cgroup out of memory: Killed process 100116 (R) total-vm:9364568kB, anon-rss:6823756kB, file-rss:585544kB, shmem-rss:0kB, UID:1000 pgtables:14892kB oom_score_adj:200
```

The initial journal command used unsupported `--kernel`; corrected read-only `-k`
succeeded. The service journal independently confirms OOM at 08:49:34 and
3min 47.906s lifetime. No new profiling run or runtime change was made. A profile
can locate executing R work, but does not by itself identify all retained objects
or prove a memory repair. Additional execution remains awaiting the requested
scope/run-budget approval.

## Approved diagnostic NyQnMl — manually cancelled, partial evidence only

The user approved one existing browser-da-profile run on 2026-09-16. It used
L04HyQ and the unchanged representative report/tJ6U1c companions under the same
limits. Driver SHA-256 remained
968e9eed5c74f7775cb8c41a0836c3b4bc969efd3de7b3607246d7d116ea5028;
the existing lifecycle-test evidence was reused. No runtime/harness edit occurred.
Invocation: 830b373de4e34d7ba0c079322fbc39f7.

Following the user's objection to machine memory pressure, the exact owned
service was stopped. The user journal records start at 08:52:00, stopping at
08:54:22 and stopped at 08:54:23. It reports 2min 22.897s lifetime, 4min 26.799s
CPU, 5.2G memory peak, 380.4M disk reads and 285M writes. These rounded journal
counters describe a cancelled diagnostic, not successful completion or a second
OOM. Reinspection confirmed inactive/dead, MainPID=0 and empty ControlGroup.
No restart or new workflow execution followed. Full-size execution is held.

The saved application.Rprof has SHA-256
18f8d5367d668811e69e47a62a6c251f6b986f21ebbce7cd0227e01b952c9c6e.
A streaming awk count of memory-profile sample-start lines (four colon-delimited
integer fields, ignoring multiline label continuations) found 4012 samples;
464 contained the exact classifier frame, 570 contained the GC frame and none
contained the exact .tokenizeProteinAccessions frame. These are nested sample
memberships, not additive costs or a complete function census. The final retained
samples show ifelse and unique.default/unique inside
classifyPeptideBiologicalExclusions, called through
srlQvalueProteotypicPeptideCleanHelper and the Q-value Apply observer.

This locates execution at cancellation in the first Q-value classification, not
a later DA/export consumer. It does not establish that the previous unprofiled
OOM occurred at the same instruction, identify retained objects, or prove that
unique itself caused the memory failure. The partial profile is diagnostic only.
No large data or S4 object was loaded to inspect it. Before further execution,
propose a diagnostic with explicit host headroom and an early-stop budget; smaller
diagnostics must not silently replace the frozen acceptance workload or gates.

Read-only source review after cancellation: the retained ifelse/unique stack is
consistent with the classifier's full-length accession-text preparation at
R/func_prot_qc_exclusion_helpers.R:253–255. Distinct-value tokenization reduces
repeated string processing, but the function still creates row-length vectors
and returns raw, annotated and filtered tables. This is a code-level explanation
of why fewer tokenization calls alone cannot establish bounded GUI residency,
not measured attribution of the OOM to those statements. GUI input/history and
other processes already resident at classifier entry remain unmeasured here.

Proposed next diagnostic, not approved or executed: isolate the unchanged installed
classifier outside Shiny on at most 38400 rows (1% of the frozen report), keeping
all input columns and existing classifier settings. Observe allocation sites and
entry/exit memory, not just elapsed time. One process, no browser or worker cluster,
one attempt, 120-second ceiling, 2 GiB hard cap and zero swap; refuse startup below
8 GiB host MemAvailable, and stop early at 1.5 GiB job usage or below 8 GiB host
MemAvailable. Those are proposed diagnostic safety budgets, not replacements for
the frozen acceptance limits. Sampling-based early stops cannot guarantee absence
of transient pressure; the hard cap is the final guard. Check the guard on a tiny
case before any data diagnostic. Reuse the existing runner rather than add another
framework, and do not execute until this bounded safety change is approved.
Exit with allocation evidence or an explicit inconclusive result; no automatic
size escalation, optimization or GUI replay. No scaling or performance acceptance
claim follows from this reduced diagnostic.

## Reduced allocation diagnostic — safety implementation and startup failure

The user approved the proposal with “di ut”. Added exclusion-allocation dispatch
to the existing runner only: fixed 38400-row installed-classifier case, allocation
records above 100000 bytes via Rprofmem, phase memory observations and existing
complete result serialization. No package runtime source changed. Only this
diagnostic uses 2 GiB/zero swap/120 seconds. The controller checks host MemAvailable
before launch and polls host/job memory every 0.1 seconds, stopping its owned unit
at 1.5 GiB usage, below 8 GiB host headroom, or on missing readings. EXIT/INT/TERM
cleanup is scoped to this diagnostic's unit. These sampled guards are not a guarantee
against brief pressure or a SIGKILL of the controller; systemd retains the hard cap
and time ceiling. Original acceptance limits remain unchanged.

Tiny shell-orchestration fixtures exercised low/missing host readings, exact job
thresholds, missing job readings, correct owned-unit stop routing and terminal
success. These are test doubles, not empirical memory-pressure tests. A real tiny
Rprofmem fixture captured a 20000-element allocation and proved later allocations
are outside the stopped profiler. Focused test runs 65759/14046 passed with one
existing export-verifier skip because MultiScholaR is not installed on the source
test library path; installed L04HyQ remains the explicit diagnostic target.

The first launch command refused on the existing serial lock and started no service.
Inspection found old controller PID 109149 still polling the stopped NyQnMl service.
SIGTERM terminated that exact controller (session 68591 exit 143); no scientific
process was restarted. The subsequently launched sole member pPCl47, invocation
a1b47e4f011945ca92f7eff1f24e6722, failed in 0.210975 main-process seconds with
`required service limits are absent`. MemoryPeak=62107648 bytes, swap=0,
CPUUsageNSec=199270000, physical reads/writes=0. It failed before loading the
MultiScholaR namespace or report. No classification or allocation evidence exists.
Raw terminal receipts and service.log remain under that member directory.

Cause: the runner correctly imposed 2 GiB, but checkpointServiceMemory still
asserted 8 GiB. Corrected that reader to accept an explicitly supplied expected
cap; only allocation diagnostics supply 2 GiB, and baseline callers still require
8 GiB. Added direct mismatch, matching-cap and dispatch tests. Focused regression
34093 passed with the same single unrelated skip. bash syntax and git diff --check
pass. Corrected runner R SHA-256:
defb1bea067bf23460a44065c42637a1982ac11a44bc622ec5b98d0024c26b18;
shell SHA-256:
e360fa1fcc3e4cddb6c6521e5a93ce3af0ca5df1cd12ecf1e0dfbe7f534b0834.

No workload retry after the startup correction: the one-attempt allowance is
consumed. Further execution needs approval. No memory reduction, scientific
acceptance or performance gain is established by these safety tests. DIA-003
remains open; no downstream gate is promoted.

## Approved reduced retry NQTK9t — diagnostic completed

User approval: “yes ffs”. Exactly one retry ran the corrected 38400-row diagnostic,
installed L04HyQ, unchanged input/settings and the diagnostic safety limits.
Member NQTK9t, invocation 2fb37cd2450f42118dd3f75260b9a6eb, completed with service
exit 0, 10.973366 seconds service-tree lifetime, CPU 10628134000 ns, peak
1019105280 bytes and zero swap. No early-stop guard fired. Classifier elapsed was
0.034 seconds. Profiling is enabled, so these are descriptive diagnostic costs,
not performance acceptance. No GUI or worker cluster was launched.

| Phase | Current service bytes | Anonymous bytes |
| --- | ---: | ---: |
| Entry | 61276160 | 60375040 |
| Namespace loaded | 979578880 | 975859712 |
| Input read | 987029504 | 980135936 |
| Classifier returned | 996102144 | 990650368 |
| Output written | 1018798080 | 990744576 |

Snapshots are not live R-heap measurements or additive peak components. The large
namespace-load footprint precedes this small classification. It cannot be assumed
to explain the previous full-GUI failure or justify package-loading changes.

Streaming analysis of classifier.Rprofmem found 72 numeric allocation records
totalling 15861136 bytes above the fixed 100000-byte reporting threshold. The
largest identical-stack group was tibble row slicing (5530560 bytes, 20 events);
unique.default/unique accounted for 2033952 bytes across six events, and ifelse
plus its nested which stacks for 1843648 bytes across nine events. Another
4609296 bytes across 27 events had only the anonymous classifier frame before
system.time; these cannot be assigned to individual source expressions from the
record. These groups partition allocation events, not simultaneous residency;
smaller allocations and memory allocated before profiling are excluded. Do not
scale these totals by 100 and call that a measured full-workflow peak.

Independent sha256sum matched both receipt fingerprints:
classifier.Rprofmem (12343 bytes),
3e78e5098d2cbff3057fdf4bda403909085ac35fbc9e7384f005cbd91b65056c;
classification.rds (23660686 bytes),
bd1f651ad5b59d341ee228430ba54cd9fe8538385ee41ba159337e452f7ccdfa.
The receipt reports 38400 biological rows, no decoy/contaminant signals and empty
manifest. This is not a new independent scientific comparison or general coverage
of exclusion cases. No large output hydration was needed for this diagnostic.

The outer controller returned 5 because its EXIT trap attempted a second stop
after normal cleanup had already removed the unit. Terminal service evidence and
controller-status.txt show successful member execution, interval decoding and first
stop; the outer error is preserved rather than relabelled as a clean runner pass.
Cleared the diagnostic traps after a successful normal stop and added an exact
one-stop regression assertion. Focused suite 31160 passed with the same existing
source-library export-test skip; git diff --check passes. No workload rerun followed.

Decision: the bounded allocation diagnostic is complete, but does not reproduce
the GUI OOM or demonstrate the required performance gain. Inspect the existing
GUI's already-resident inputs, temporary tables and history ownership before any
new experiment; do not recursively optimize unique/ifelse or package loading from
this result. No more execution or runtime change is authorized by this completed
item. DIA-003 remains in progress and full-size runs remain held.

## Read-only first-Q-value ownership review

Reviewed existing runtime sources without loading data or executing the GUI.
The following owners distinguish the GUI from the isolated classifier diagnostic:

| Owner / point | Existing code | Implication, not a measured byte count |
| --- | --- | --- |
| Import state | mod_prot_import_processing_helpers.R:112–127 assigns data_tbl, then data_cln | Two retained references initially; not proof of two deep copies |
| Design checkpoint | mod_prot_design_state_helpers.R:11–18,93–99 constructs the DIA S4 from data_cln and saves raw_data_s4 | The initial S4 remains reachable alongside workflow data |
| History assignment | utils_workflow_state.R:81–110 assigns data into states | No explicit serialization/deep clone here; copy-on-modify and shared columns prevent counting every reference as independent bytes |
| Q-value caller | mod_prot_qc_peptide_qvalue.R:78–180 obtains currentS4 and calls the filter | Current state remains reachable during filtering |
| Filter temporary | func_prot_qc_peptide_confidence_filters.R:224–258 constructs a row mask and q_valid_rows before classification | Original input, mask and selected/filter-result table can overlap in lifetime |
| Classifier result | func_prot_qc_exclusion_helpers.R:306–357 annotates and subsets tables | Row vectors and result tables add working allocations; tiny profile confirms allocations, not full-GUI residency |
| Later commit | mod_prot_qc_peptide_qvalue.R:208 onward saves qvalue_filtered after filtering returns | This new checkpoint has not yet been saved at the cancelled profile's classifier samples |

Also, readConfig in func_general_filemgmt_config_io_helpers.R:129–130 eagerly
creates a cluster and loads packages when number_of_cpus is present. The earlier
kernel OOM task list contained eight additional R processes. The code and task list
are consistent with worker overhead, but do not establish their exact shared/unique
memory or prove which call created each PID. No worker-topology change is made.

Decision: the immediate diagnostic gap is surrounding first-step input/temporary
ownership, not accumulated later normalization/DA checkpoints. Offloading historical
states may still benefit later steps; it cannot simply be credited with solving
this first-step peak. The 1% isolated classifier run omits the GUI owners and the
Q-value preprocessing path, so its success cannot qualify the original GUI control.
Do not sum table object.size values or worker RSS to invent a peak. Before another
run, any proposed probe should cover the actual Q-value helper with its input held
and distinguish shared objects from allocations, under the already defined small
diagnostic safety budget. That would be a new bounded experiment, not automatically
authorized here. No source edit, new experiment or gate promotion followed this
read-only review.

## Continued DIA-003 implementation — Q-value probe and scenario validation

Following the user's explicit direction to continue finishing the ticket, ran one
guarded Q-value-helper probe through the existing exclusion-allocation operation
with qvalue selector. This uses the installed production helper with the exact
nine projected columns from zqsj4M's Q-value configuration, unchanged confidence
thresholds and exclusions, original input retained, 38400-row prefix, no GUI or
workers, and the existing 2 GiB/1.5 GiB early-stop/host-headroom/120-second guards.
No package runtime changes. New dispatch tests and the focused harness passed
(52179), retaining the one source-library export-verifier skip.

RlYY55 (invocation f561717434774b8592f373ebe2a5028d) completed with service and
controller exit 0, peak 1018540032 bytes, swap zero, CPU 8594343000 ns and service
lifetime 8.774578 seconds. Helper elapsed 0.114 seconds; input 4437960 bytes,
result 11199016 bytes, shared union 15505448 bytes via lobstr::obj_size(input,
result), so aliases are not simply added. This does not measure peak temporary
residency or include GUI/S4/worker ownership. The profile reports 170 allocations
above 100000 bytes, totalling 38625056 bytes; allocations are not retained memory.
Receipt fingerprints: output df781627f83e5414c8d6a23edf79c13c6ca8f28ee1421b7dda0a526a6cc23882;
profile e53e21fb095548dde4a2f29154815a141cbeac3a93fb80a21615678cb8eaedd4.
All 38400 rows survived; no independent scientific comparison is claimed. The
probe did not reproduce the full-GUI failure and is now closed with no size escalation.

Returned to unfinished harness requirements. The baseline CLI previously never
read --contract. It now validates the declared scenario schema/source, all seven
recipe identities, dimensional consistency, all C01–C26 recipe references and
capacity-only exclusion, and binds both scenario/protocol files to their frozen
SHA-256 values. A not-yet-qualified control is not rejected as a baseline
prerequisite. Limpa is explicitly held without IQ fallback. This validation does
not make the unfinished baseline execution dispatcher runnable.

Focused suite 44312 passed with the same single source-library export-test skip.
--help passed; --validate-only read the genuine bound inputs then correctly failed
at the remaining baseline implementation hold, without executing a workload.
git diff --check passed. Reviewed passing synthetic-receipt tests cover pass/fail,
borderline and exact thresholds, OOM/timeout, missing/nonfinite metrics, stale
build/route/boundary identities, and repeated-session non-independence using the
frozen paired uncertainty procedure. Only that fully evidenced acceptance checkbox
was marked complete. Baseline orchestration, full output boundary and other unchecked
criteria remain implementation work; DIA-003 is not complete and G0/G4 stay open.

## Baseline input resolution and serial dispatch core

Continued the existing DIA-003 runner, without adding a backend, dependency or
service mode. resolveCheckpointBaselineInputs consumes the already-existing
inputs.json preparation format at OUTPUT_ROOT/recipe/inputs.json. It requires all
six timing-eligible recipe bindings, checks the full recorded recipe against the
frozen contract, canonicalizes nonempty owned report/design/FASTA paths and validates
declared SHA-256 syntax. Only the small preparation receipt is hashed here: payload
integrity verification still belongs inside the eventual measured service, so this
preflight does not read/warm the full report. A valid structural binding is not proof
that the payload matches its declared hash.

dispatchCheckpointBaseline implements serial ordering across the frozen six
recipes and seven blocks, requires each result to be recorded before another
member starts, rejects reused/missing run identities, and records then stops at the
first failure. This core accepts the guarded-service adapter and result recorder;
tests use doubles, not mocked science represented as a real benchmark.

Test 29430 failed because the new test referenced a protocol fixture helper defined
later in the same file. The test now reads its protocol directly. Corrected run
21526 passed, with the existing single source-library export-verifier skip. Tests
cover all 42 dispatches, serial record ordering, failure retention, no retry,
missing recipes, duplicate run IDs, recipe mismatch and complete C01–C26 bindings.
git diff --check passed. CLI --validate-only on the current proof root fails at
the concrete missing small/inputs.json binding, without package loading or a
workload. Existing prepared files remain untouched in their historical member
directories; no data were regenerated or copied.

This completes only the bounded dispatch-core item. The actual guarded-service/
completion-output adapter is still unqualified and CLI execution remains held.
Do not call the full baseline implemented, DIA-003 complete, or G0 accepted from
these tests. Continue with that adapter next; do not reopen memory profiling.

## Guarded-service and separate-verifier handoff

Added runCheckpointService/runCheckpointBaselineMember in the existing runner.
The invoker uses the current shell service entry points and resolves the exact
owned directory from launch/terminal unit markers, never from newest-directory
guessing. It preserves nonzero controller exit (including an observation timeout)
and does not restart. The member adapter requires successful service-tree interval
evidence with matching unit identity, expected installed-library path, matching
input fingerprints, untouched QC defaults and profiling disabled before starting
the separate verify-da process. The verifier must reference the same member and
exit successfully. Failure evidence remains associated with the member directory.

Browser input fingerprints are now calculated inside the measured service before
app launch, not by prewarming the report in the outer dispatcher. Prior browser
receipts lacking these bindings cannot satisfy the new adapter; their earlier
diagnostic/scientific evidence remains intact. This changes a future measured
boundary and requires matching control/candidate measurements, not reuse of old
timings under the new identity.

The adapter deliberately returns held after the existing verifier passes: that
verifier covers session/history, normalization, DA and heatmap checks, not every
required baseline output/restore consumer. No synthetic complete_verified receipt
is manufactured. Full output-boundary qualification and the final CLI connection
remain unfinished. Main baseline execution is still disabled.

Test 98544 failed because a side-effecting run() call inside a test assertion was
evaluated again while the assertion built its diagnostic context. Tests now assign
the result before assertions. Corrected 36555 passed; final 39609 also tests exact
unit resolution, ambiguous unit rejection and retention of nonzero observation
status. All focused tests pass with the existing one source-library export-verifier
skip. No browser, scientific workload or large-file verification ran for this item.
git diff --check passed. No runtime package code or performance gate changed.

## Required automatic QC output check — production omission reproduced

Added checkCheckpointIqQcOutputs to the existing verifier. It checks nonempty
regular files for IQ's implicit working-directory qc-plots.pdf, IQ input/output,
replicate-filter table, the fixed RUV-skip normalization outputs/composite,
pre/post normalization panels, core session metadata, and all seven peptide-QC
combined plots in both PNG/PDF. Duplicate timestamp-directory matches are rejected;
resolved paths must remain inside the member. This is explicitly presence evidence,
not content/scientific validation or the complete baseline inventory.

The small retained member 5Olm86 has the fixed IQ/normalization/session files but
none of the fourteen required peptide-QC composite PNG/PDF files. Its retained
browser-logs.json repeatedly records:

```text
Object 'project_dirs' not found (neither as argument nor in Global Env). Plots will not be saved. Ensure 'setupDirectories()' has been run.
```

Read-only checker invocation rejected 5Olm86 and named all fourteen files. No large
file was decoded or hashed and no GUI rerun occurred. The producer's
func_prot_qc_filtering_progress_helpers.R:109–155 path resolution explains this
warning: without explicit/global project_dirs it disables saving. This is a
production output-path handoff defect, not a missing scientific result or reason
to fabricate outputs in the driver. Earlier passed_small_da_export_checks covered
only their declared exports and remain historical limited-scope evidence.

Focused test 53501 passed with the existing one source-library export-verifier
skip. Tiny fake-file fixtures prove missing/empty/ambiguous detection and explicitly
do not claim valid image/RDS contents. The new check is called by the existing
normalization/DA verifier before hydration. Runtime source was not changed.
Complete inventory coverage and the production plot-path repair remain necessary;
the latter needs an explicit runtime scope amendment. No unchanged GUI replay,
threshold relaxation, new profiling or acceptance promotion follows this finding.

## Current installed and baseline decision — 2026-09-16

The plot-path repair passed affected module contracts (772 assertions, zero
failures/warnings/skips) and installed-package S4/export/load and tiny PNG/PDF
write checks. Fresh installed tYrWLW (runtime patch affca9d4...) completed tiny
summary member dPgDsR and independent verifier EgdKnb. Expanded checking found
all 84 peptide-QC source plots but rejected its missing publication QC-composite
copies. The narrow prefixed-name fallback and correlation-Skip copy condition
passed 108 affected publication-helper assertions with three expected warnings
from injected backup failures. Fresh installed lc00wd (runtime patch
95fd659e0f8a79cda0c7c7e2a9d5c102bfe3af721ab4b76ed2e72227dcb96135)
then completed tiny summary member 5RSyuB and independent verifier rTv3rw.

The final tiny verifier passed 108 required IQ/QC files (84 filtering plots plus
24 other IQ/QC/metadata files), thirteen Summary files and eleven selected
Volcano/NumSig/Heatmap publication figure copies. It decoded the complete
16-checkpoint legacy session (64 proteins, twelve samples), matched DA TSV/XLSX
against the installed registered method, and rechecked heatmap content. Source
and publication QC PDF/PNG fingerprints matched separately. Verifier peak was
about 1010.3 MB, zero swap. This is strong small-route output/science evidence,
but not a full-size baseline or an independent IQ scientific oracle.

Matched tiny members dPgDsR and 5RSyuB used the same IE9ohP inputs, driver,
settings and 5 GiB diagnostic cap. The earlier member lacked required QC copies.
Service lifetime was 158.597 versus 162.313 seconds; exact journal cgroup peaks
were 3406376960 versus 3421167616 bytes; CPU was 288352401000 versus
295365752000 ns. The new completed file work cost 3.716 seconds and 14790656
bytes more peak in this single descriptive pair. It cannot be used as an accepted
comparator or evidence for the frozen >=10% full-route improvement gate. Other
host/cache conditions were not paired; no causal decomposition is claimed.

Baseline CLI input discovery now consumes existing immutable preparation receipts
under the proof root, avoiding copied/generated placeholder fixtures. Shell and
R handoff tests passed with one known source-library export-verifier skip.
--validate-only recognized the only prepared small/representative recipes, then
correctly refused the missing tall input. Limpa remained held. No large member or
benchmark was launched through this CLI invocation.

Tiny installed default-QC member qchLuu, using lc00wd/IE9ohP and a separate
5 GiB/zero-swap safety cap, passed the seven peptide actions but failed at the IQ
protein step. Its peptide_values_imputed.tsv is 271 bytes and one line (header
only); iq_output_file.txt does not exist. The browser receipt records the IQ
timeout error. Service lifetime was 111.526 seconds, cgroup peak 2862940160
bytes and zero swap. This is a valid unfavorable tiny default-QC result, not
an instrumentation failure. The frozen settings/data were left unchanged and no
retry occurred. The earlier relaxed tiny member remains diagnostic only.

Read-only audit of mod_prot_qc_peptide_sample.R:43–47 found its packaged 500
distinct-peptide-per-sample cutoff. The frozen small recipe has 64 proteins and
three distinct peptides/protein, so each sample has at most 192 distinct
protein-group/stripped-peptide identities; repeated forms are not additional
identities. That bound is below 500 even before earlier filters. No assay or
large file was hydrated to establish this. Raising the recipe or lowering the
cutoff would change the frozen control and was not done.

The representative default-QC control/candidate attempts remain OOM/timeout
evidence under the frozen 8 GiB/zero-swap/900-second limit, not finite timing
denominators. Tall, wide, annotation-heavy and cycle inputs are not prepared here.
The package/harness work is materially advanced, but DIA-003's complete baseline
proof and DIA-004's qualified G0:iq remain unmet. G0-required DIA-007 onward
cannot be accepted or integrated from these tiny checks. The next safe work is to
audit the existing frozen input and first-Q-value memory conflict before choosing
any full-size control treatment; do not weaken populations, skip required outputs
or claim a performance gain from the completed diagnostic route.

Read-only scope clarification: mod_prot_qc_peptide_server already receives
experiment_paths (R/mod_prot_qc_peptide.R:102), but its
runProtQcPeptideSubmodule handoff passes only id, workflow_data, omic_type and
experiment_label (lines 82–85). The seven child modules therefore cannot supply
those paths to updateProteinFiltering's existing project_dirs parameter. The
Q-value child illustrates the chain at mod_prot_qc_peptide_qvalue.R:376–392 and
its output helper at lines 251–270. This identifies a concrete lost argument,
without claiming that fixing it also qualifies other output consumers.

Proposed repair scope: explicit optional path arguments through the peptide
orchestrator and its existing seven children (qvalue, protein, rollup, intensity,
sample, replicate, impute) into the existing plotting helper; corresponding
existing orchestrator/module tests and one tiny output-write check. Preserve
standalone defaults, public call compatibility, plot calculations, filenames,
scientific populations and project isolation. Do not create a global path fallback,
new registry or plotting abstraction. Protein-QC and summary/report consumers
require their own demonstrated caller review; this is not blanket authority to
rewrite them. No runtime edit or new execution was made during this clarification.

## Authorized peptide-QC plot-path repair — implemented, verification incomplete

The user directed completion without further approval requests. DIA-003 revision
1.9 now owns the peptide orchestrator, its seven children and the two existing
module/orchestrator contract test files. The orchestrator forwards optional
experiment_paths; each server/Apply/output chain passes those paths explicitly
to updateProteinFiltering's existing project_dirs argument. The intensity bootstrap
and imputation argument-list paths are included. Null paths preserve legacy call
arguments, including injected helpers. No new globals, path registry, plotting
algorithm or scientific method was introduced.

Added tiny tests for all seven output adapters, all seven Apply observers, project
A/B isolation, legacy null-path calls, orchestrator forwarding and actual PNG/PDF
creation with a four-row table. Runtime files parse successfully. Two test services
used MemoryMax=2G, MemorySwapMax=0, RuntimeMaxSec=120; each peaked at the journal's
rounded 1.1G with zero swap. No GUI workflow or full-size data ran.

Service multischolar-dia003-path-tests-01 (session 89654): 1066 assertions passed,
one error, no warnings/skips, 17.275 seconds service lifetime. Error: an existing
orchestrator mock did not accept experiment_paths. The real tiny writer produced
both composites and its individual plot files.

Service multischolar-dia003-path-tests-02 (session 69386): failed with ten assertions/
errors and one warning after the attempted mock correction. An imprecise patch
added a captured experiment_paths field to the unrelated Q-value mock instead of
the intended orchestrator capture; that undefined variable caused the dependent
failures. The orchestrator also did not capture its new argument. Lifetime 15.863
seconds. The new seven-observer/output tests and real tiny write check ran before
the failures; this does not make the suite a pass.

Restored the unrelated Q-value mock and added the field to the exact orchestrator
capture block, verified by git diff. git diff --check passes. The one bounded
correction verification failed, so no additional test/run attempt was made in
this item. Current correction is not yet execution-verified; do not claim a clean
regression pass, installed-build qualification or DIA-003 completion. Earlier failed
tests and missing-output evidence remain retained. Next required verification is
the corrected affected suite, followed by installed small-route output validation
under an explicitly safe boundary; full-size execution remains held.
# Recorded-state / complete active-QC pilot wiring (2026-09-16)

The existing service runner now supports recorded-prepare, recorded-member,
recorded-qc and recorded-verify, with validation-only, fixed 8 GiB/zero-swap/
900-second limits, a 7.5 GiB peak-aware stop and explicit host-headroom checks.
Small lifecycle tests cover full replay, late failure, resource closure, actual
Q-value Apply/configuration/audit execution, profile capture and owned-service
cleanup. The complete Q-value wrapper uses the unchanged production methods and
an explicit private configuration environment. It excludes plots and later states,
and reports allocations separately from retained owners and service peak.

The user subsequently authorized scoped continuation after the two preparation
failures. Preparation SxedsH and the complete active-QC member F7K2nD now passed;
independent verifier trprB2 matched the recorded qvalue_filtered S4 exactly after
excluding only the already-documented audit timestamp_utc fields. The source/build
bindings, preserved failures, native history and conditional qs2 comparison are
recorded once in
[checkpoint-results.md](checkpoint-results.md#recorded-session-pilot-2026-09-16).

The active-QC diagnostic used the full recorded raw S4 from the 960,000-row,
24-sample synthetic DIA journey, with recorded thresholds and production
configuration update, filter, audit and R6 commit. No later history was loaded.
It completed in 102.731894 service seconds at 3,224,420,352-byte combined peak
with zero swap. The profiled Apply body took 80.790 seconds and 76.628 CPU seconds;
these are profiling diagnostics, not an unprofiled throughput comparison.

| Observed boundary | Unique input/output/config owners | Service peak so far |
| --- | ---: | ---: |
| Active input | 121,563,920 B | 1,131,241,472 B |
| Filter method entry | 121,565,872 B | 1,157,169,152 B |
| Method return, before audit/commit | 271,337,736 B | 1,785,159,680 B |
| Audit and state commit returned | 271,351,792 B | 3,224,420,352 B |

The 3,863 genuine CPU sample records (20 ms interval) classify into 3,131 with
.emitPeptideQcAuditRecord on the stack (~81.1%), 597 in qvalueFilterFn excluding
ownership observation (~15.5%), 81 in observation (~2.1%) and 54 elsewhere.
Memory sample prefixes, rather than multiline physical-line counts, define this
classification. Stack-sampling shares are approximate and are not added to
service CPU or treated as a full GUI cost distribution.

Rprofmem contains 966 numeric records above the 100,000-byte threshold, totalling
5,717,599,512 allocated bytes. Audit stacks account for 4,729,327,304 bytes:
3,272,856,520 in digest/canonicalization, 1,058,134,848 in removal-ledger work,
and 398,335,936 in summaries. Filter-method stacks account for 976,043,488 bytes;
the remaining two records total 12,228,720 bytes. These are cumulative allocation
events, **not simultaneous residency**; smaller allocations and native allocations
are not fully represented. The allocation trace itself is 181,556,572 bytes and
its I/O/cache is part of this diagnostic service footprint.

Source review matches the observations: .peptideQcCanonicalise creates character
row keys over complete data frames and sorts them before hashing;
.emitPeptideQcAuditRecord also calculates the identical substantive-record digest
twice. This is a concrete future optimization target with exact provenance to
preserve. No audit logic or scientific method was changed by this pilot.

Retained raw profiles in browser-smoke.F7K2nD have SHA-256:
Rprof `0f5dd9f97be338a0dcca2ad1abbe6c24bfdc0a665f7f1216404179d9971a27cf`;
Rprofmem `1491eab645f3e7cc19d76ab4c226c4a12f75fc360cf75bfdce6f37fade36e11f`.
The result supports separating historical residency from active audit temporaries.
It excludes GUI/browser, plots and other parallel stages, and supplies no G3/G4 or
complete-workflow performance acceptance. All earlier failed evidence remains.
