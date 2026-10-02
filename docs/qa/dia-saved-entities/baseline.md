# DIA IQ control evidence and scoped G0 decisions

Status: frozen-v1 `G0:iq` held; reduced-development-v2 `G0:iq` accepted only for
downstream standalone IQ saved-entity work. DIA-003 remains in progress and
DIA-004 is not complete. This document does not requalify the frozen control or
change its approved protocol.

The frozen source is janitor
`c7c12851ea4d9e96f91df0e1f9e7b91c2eb51017`. The scenario and protocol
SHA-256 values are respectively
`cfdf7da2fc5a6123da269f85011550d8ffd2e8a611de071e0cd1a52a6f35d152`
and `517edccb0906165698f3c1f165dfe7d961a17910d93363a60e843fa466c62379`.
The required representative input is a dense synthetic 48-sample,
3,840,000-precursor report (SHA-256
`2414125be92a11fbe4342fdf354e0214915e30592bec66f3e9d6a984f9290841`),
not an empirical model of the user's moderate-size real datasets. The protocol's
service cap remains 8 GiB, zero job swap and 900 seconds. Only a completed installed
GUI journey with required files and restore can be a finite elapsed/memory
denominator.

| Control/evidence member | Build/input/boundary | Observed result | Qualification |
| --- | --- | --- | --- |
| `browser-smoke.FAORgA` | Pristine janitor installation `b0-install-RNROup`; small IE9ohP report; installed GUI peptide and IQ path | IQ wrote an output, then the protein tracking consumer failed on a valid `Protein.Group` key. Reported 102.002 s, CPU 199.324 s, peak 2,742,398,976 bytes, zero swap. | Failed science/output control; no finite journey completion |
| `browser-smoke.bQWqlM` | Corrected pre-storage VEqJKP; frozen representative report and tJ6U1c design/FASTA; default QC; 8 GiB/zero swap/900 s | `Result=oom-kill`, peak 8,589,934,592 bytes, zero swap, CPU 1,581,538,348,000 ns; main lifetime about 753.424557 s. | Capacity failure; not a finite comparator |
| `browser-smoke.zqsj4M` | Installed L04HyQ identifier repair; same frozen representative input, default QC and limits | `Result=oom-kill`, peak 8,589,934,592 bytes, zero swap, CPU 435,638,441,000 ns; main lifetime 227.539417 s. Kernel evidence shows anonymous memory dominated at failure. | Candidate capacity failure; not G0 or a gain |
| `browser-smoke.qchLuu` | Installed lc00wd; frozen small IE9ohP input; default QC; separate 5 GiB diagnostic cap | Seven peptide actions ran, but the IQ input TSV had only a header and no IQ output. The packaged sample filter requires 500 distinct peptides/sample; this input has at most 64 × 3 = 192. Peak 2,862,940,160 bytes, zero swap. | Invalid tiny default-QC control |
| `browser-smoke.5RSyuB` plus verifier `rTv3rw` | lc00wd/IE9ohP, relaxed tiny filters, 5 GiB diagnostic cap | Completed and separately checked 108 IQ/QC files, thirteen Summary files, eleven selected figure copies, legacy 16-state RDS and 64-row DA. Peak 3,421,167,616 bytes; service lifetime 162.313 s, zero swap. | Valid wiring/science diagnostic only; settings and input differ from G0 |

The user-requested reduced-input trial is separate from the frozen controls. It
selects twelve control and twelve treated samples, keeping technical-replicate
pairs, and one unmodified precursor form per peptide. This preserves all 10,000
input protein groups and 40,000 stripped-peptide identities while reducing the
long report to 960,000 rows. The derived report and design SHA-256 values are
`855abc455f5eaaaf8b617753bf7899ba418d575f234d61505960fd03b28a6ae0`
and `acff81ff68ee02aff56f54b7e11d40c47ea2b43bf1b9ba2aa592f75d91bce18e`.
The original FASTA is unchanged. Preparation `browser-smoke.4OZChD` streamed the
data outside the GUI timing boundary at a 257.2 MB service peak and zero swap.
The [diagnostic input receipt](../../../../MultiScholaR-dia-proof-20260911-VhC1n5/browser-smoke.4OZChD/halved-diagnostic-inputs.json)
contains the exact selected runs, counts and file fingerprints.

The first reduced GUI member `browser-smoke.4X0nXL` used a 6 GiB hard cap and a
5 GiB early stop. The controller watched current memory only; a transient peak
crossed 5 GiB before its next poll. It was manually stopped. This is an invalid
completion member and a retained safety-instrumentation failure, not a service OOM
or a benchmark result.

After the user's explicit instruction to raise the stop threshold, the same input,
installed lc00wd build and GUI settings ran once in `browser-smoke.EJlArx` with a
10 GiB peak stop, 11 GiB hard cap, zero job swap and a host-headroom stop below
8 GiB `MemAvailable`. The guard watched both current and peak memory. It stopped
at 10,757,476,352 bytes of cgroup peak (about 10.02 GiB), with zero job swap and
CPU 1,596,207,185,000 ns. Host `MemAvailable` at the stop was 9,721,065,472
bytes. The unit is now inactive, with no process group. The host swap observation
rose during the run; attribution to this job alone is not established.

This second reduced run completed the seven peptide plots, wrote an IQ output with
10,000 data rows, normalized output TSV/RDS and two 239,747,964-byte legacy session
RDS files plus summary text. It stopped before DA and Summary GUI actions, leaving
the application receipt `running`; the guarded terminal action is authoritative.
A separate 2 GiB/zero-swap installed-package check decoded only the final protein
RDS, validated its S4 object, and found exact RDS/TSV identifiers, sample order,
missing/nonfinite masks and values (scientific tolerance `1e-12`) for 8,020 final
proteins across 24 samples. The complete session history was not hydrated under
the tight verification cap, so file existence and size do not certify its content.
The first verification command incorrectly assumed all 10,000 input proteins
survived QC; that test-authoring failure is retained. The corrected bounded check
passed at an 856 MB service peak with zero swap. No DA result or automatic Summary
inventory exists for this member.

The observed sources of cost are deliberately separate: a landing-only installed
app peaked around 1.1 GiB; the corrected-control profile captured repeated
identifier work inside first Q-value Apply but did not identify retained heap
ownership; the reduced member's cgroup grew toward 10 GiB during later IQ/protein
and compatibility export work. The two full-history RDS serializations are part of
the measured GUI operation and cannot be substituted with reference-only files.
Application CPU samples, service CPU, cgroup peak, file-cache and anonymous-memory
snapshots overlap in scope or time and are not added into an invented breakdown.
The cost of external fixture generation and separate verification is excluded from
member completion, recorded independently above. Remaining unattributed cost and
host-cache variation prevent scaling claims.

Frozen-v1 decision: `G0:iq` is held. The pristine control failed a required consumer, the
corrected full representative control failed the fixed 8 GiB limit, and the frozen
tiny default route is scientifically ineligible. The user's 960,000-row/11 GiB
diagnostic got further but did not complete and changed both the input and service
limit. It cannot replace the frozen denominator or establish a 10% full-journey
elapsed/combined-peak improvement. Limpa remains held under DIA-017, with no IQ
evidence inheritance. No G0-required downstream ticket is accepted by this file.

The recoverable source identities, installed libraries, failed member receipts,
legacy outputs and proof root are retained. No branch reset, deletion, format
change or cap change to the frozen protocol was made. A later control decision must
state its workload, source build, scientific/output boundary and service limit as a
new qualification scope; old failures remain part of that decision.

## Reduced-development-v2 IQ control

The user directed use of the completed reduced journey as G0 for downstream work.
This is a new, bound development-control scope, not a relabeling of frozen-v1 or
pristine janitor. Its input is the 960,000-row report and 24-sample design above,
with the unchanged FASTA; the three file SHA-256 values are respectively
`855abc455f5eaaaf8b617753bf7899ba418d575f234d61505960fd03b28a6ae0`,
`acff81ff68ee02aff56f54b7e11d40c47ea2b43bf1b9ba2aa592f75d91bce18e`
and `352b56704c330f3e59646101df2c52a469e5ebe17da652fdd0d3ea85e672f6d4`.
The [derivation receipt](../../../../MultiScholaR-dia-proof-20260911-VhC1n5/browser-smoke.4OZChD/halved-diagnostic-inputs.json)
records one unmodified precursor form per stripped peptide, all 10,000 source
protein groups and 40,000 stripped-peptide identities, and twelve complete
control/twelve treated technical-replicate sample pairs. These changed input
populations are not interchangeable with frozen-v1's 48 samples/two forms.

The installed build is MultiScholaR 0.5.0 in
`browser-smoke.lc00wd/library`, source parent janitor
`c7c12851ea4d9e96f91df0e1f9e7b91c2eb51017` plus runtime patch SHA-256
`95fd659e0f8a79cda0c7c7e2a9d5c102bfe3af721ab4b76ed2e72227dcb96135`.
It is a corrected source, not untouched B0. The [completed GUI receipt](../../../../MultiScholaR-dia-proof-20260911-VhC1n5/browser-smoke.rzEEWj/smoke-receipt.json)
reports `passed_summary_ui_only`, installed-package loading, default GUI QC,
real IQ, cyclicloess, RUV Skip, correlation Skip, DA `~0+group` at q=0.05/TREAT0,
static and interactive volcano, heatmap save, Summary Save Arguments, publication
copy and session-state export. The original scientific action protocol SHA-256
remains `517edccb0906165698f3c1f165dfe7d961a17910d93363a60e843fa466c62379`;
this v2 development control alone uses a 16 GiB combined service cap, zero job
swap, 900-second ceiling, 15.5 GiB peak-aware early stop and 2 GiB host-headroom
guard. The original 8 GiB protocol is unmodified.

The [service-tree interval](../../../../MultiScholaR-dia-proof-20260911-VhC1n5/browser-smoke.rzEEWj/terminal-interval.json)
is authoritative for completion: 793.817321 seconds from startup through outputs
and teardown, CPU 1,691,107,886,000 ns, combined cgroup peak 12,765,577,216
bytes (11.89 GiB), zero job swap and successful service result. Fixture creation
and the separate verifier were outside that timed boundary. Host swap use rose
during the journey; there is no justified attribution to this service alone or
claim of zero host swap. The cgroup peak includes charged file cache/kernel
memory; service CPU, process RSS and memory.stat components cannot be summed
into a larger or phase-specific peak. The first/last namespace and DLL inventories
are bounded process-runtime identities, not binary-content hashes.

The [separate installed-package verifier](../../../../MultiScholaR-dia-proof-20260911-VhC1n5/browser-smoke.Mq2n1b/verification-receipt.json)
passed at a 2,992,152,576-byte service peak and zero job swap. It decoded the
complete self-contained 16-state legacy RDS and checked restore/revert content;
verified all 108 required IQ/QC files, thirteen Summary files and eleven selected
publication figure copies with source/copy fingerprints; checked the 8,020-row,
24-sample final protein result and 8,020-row DA table against the **same installed
method**; and checked the saved 50-row heatmap with four clusters. This is an
independent process/output verification, **not** an external IQ estimator or proof
of every normalization/RUV/correlation scientific case. The verifier's 29.439317
seconds and memory are not added to GUI completion.

Decision: reduced-development-v2 `G0:iq` is accepted as one scientifically usable,
complete-output *development control* for the exact installed build/input/settings/
16 GiB boundary above. It releases standalone DIA IQ saved-entity tickets starting
with DIA-007; it does not approve GUI wiring, G2/G3/G4, Limpa or any scaling claim.
The single completion is a finite descriptive observation, not the seven
independent matched blocks required for performance inference, and no qualified
pre-component comparator exists for this revised scope yet. The original fixed
>=10% elapsed/combined-peak and immediate-comparator non-regression requirements
remain unproven. Frozen-v1's full-size 8 GiB route is still held; its OOM failures
are not denominators. Any changed build, workload, settings, dependency set or
timed boundary needs its own qualification before a candidate claim. Missing
cost attribution and the real moderate-workload correspondence remain explicit
uncertainties, not inferred gains.
