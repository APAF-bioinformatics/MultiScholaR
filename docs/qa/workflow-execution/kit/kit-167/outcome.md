# KIT-167 outcome (implementation, ticket revision 1.10)

Branch `kit-167-reference-mapping`; `v0.9-staging` at `ea93d948` merged in at d3909a30 (KIT-173 and KIT-175 merged).
State: the implementation that needs neither the validation harness (KIT-159) nor the distributable is built and
tested; the parts listed under "Not built" and "Waits for the harness" are not.

## What was run (each run's purpose, commit, minutes and result: `runs.tsv`)

| Run | Commit | What it showed |
|---|---|---|
| 45334675816 | b21337c7 (tests only) | the five core test files fail first: 50 tests, every one an error (functions absent) |
| 45335657533 | f56db86d | core after the fix: 185 expectations passed, 0 failed (resolver 48, FASTA 23, organisms 34, annotation 42, engine 38) |
| 45336217501 | 378271b2 (tests only) | the workflow and enrichment files fail first: 14 errors and 8 failures |
| 45336550337 | 278f30f5 vs d3909a30 | the seven new files pass (288 expectations, the strain-project enrichment run included); of 12 dependents, 9 `same`, 3 `new-failure` (below) |
| 45338825064 | c67490b5 | gate F2 passed (`feasibility.md`) |
| 45339028522 | af7a4843 vs d3909a30 | correction pass and checkpoint: full suite base against head: 467 files, 459 `same`, 7 `only-on-head` (the new files, 293 expectations, 0 failures), 1 `new-failure` (`test-gui-entry-proteomics.R`, lane `prot_lfq_fragpipe`; below); R CMD check identical on both sides (0 errors, 3 warnings, 3 notes; none only on head); gate-off script failed on its own mechanics (below) |
| 45347302450 | 87050bf7 vs d3909a30 | instrument run: gate-off with the tabset ids masked passes (pages and configuration values identical, head deterministic); the GUI-entry file: head passes every lane, base fails lane `prot_dia_limpa` by the same unmasked export id (`native-export-4cadd0e5a`) |

The correction pass (one, as the ticket allows): run 45336550337's three new failures were caused by this ticket and by
nothing else: (1) the three reference names had been added to the enrichment declaration's parameter list, which
changed the expected provenance names and one state digest of existing runs (`test-omic-enrich-module.R` R3 P-28,
`test-workflow-analysis-record.R` D6 test 3, `test-omics-enrich-core.R` kernel_clusterprofiler states); (2) the core's
returned list had gained an element (`test-omics-enrich-core.R` W2 names). Fix af7a4843: the declaration is
unchanged and the three names are appended to a run's declared names only when the reference identifiers are used;
the returned list is unchanged.

## Pass 2 (coordinator's rulings of 2026-11-24)

| Run | Commit | What it showed |
|---|---|---|
| 45352132631 | 644aa7df (tests only) | download, GUI and replay tests fail first (download 4 errors; GUI 14 failures, 7 errors; enrichment 5 failures, 1 error) |
| 45352515231 | 4a501cb0 vs d3909a30 | new files pass but 3 errors (the GUI test and the gate-off script built the enrichment page from the whole configuration; the 2 x 2 counts read under a contrast name); dependents `same`; correction pass taken |
| 45353194826, 45353725258 | packaging | gate F1 (`feasibility.md`): the first an instrument defect (Docker engine not reached), the second complete |
| 45352823145, 45372182412, 45386584323 | e5a258d2 to the last fix | the inputs script: F3 and F6 in the first; the second timed out (instrument); the third complete: FASTA locations, F4, the GEO cut, the demo tables and the real-engine lane |
| 45387203510 | 3c4c5c90 vs a1b397b1 (merged with staging) | every KIT-167 file passes (369 expectations); the 13 dependents `same` |
| 45392214469 | 2893bfed vs a1b397b1 | the pass's one full suite: 475 files, 465 `same`, 10 `only-on-head` (KIT-167 files, 384 expectations, 0 failures), no new failure (the GUI-entry mask flake did not occur this time; recorded as unchanged); R CMD check identical both sides; gate-off: pages and configuration values identical, head deterministic, no item failing to render |
| 45386641894 | 3c4c5c90, bundle lane (`distributable.yml`, stage full, both platforms) | DIAMOND fetched and verified against the pin in both builds; the bare-container check finds it through the launcher's `MSR_DIAMOND` and the package's check accepts it (Linux); the Windows NSIS installer builds with `installer.nsh`; the GUI-entry lanes fail as in earlier full runs on other branches (Linux: prot_dia 4 steps, prot_dia_limpa 6, trx_screen 13; Windows: metab and lipid 8, trx_screen an error), with no base run on this image, so not judged; trx_screen's discrete differences 83 against 82 in run 45151293709 (23 November, another tree) are not attributed |

Built in pass 2: `downloadReferenceProteome()` and its refusals; the import panels (proteomics mapping, transcriptomics
annotation step), the tab's 2 x 2 table, download and taxon pre-fill (all behind the gate); the record and replay of
reference-identifier runs (ruling (c)); the demo tables and `readReferenceMappingDemo()`; DIAMOND in both bundles with
its checks, the installer's runtime check, the notices; gates F1, F3, F4, F6; the GEO cut test (ruling (e)).

## Acceptance criteria

| # | State | Evidence |
|---|---|---|
| A0 | met (step 1) | `feasibility.md` (runs 45127082405, 45127086068); the `LICENSE` sentence added in this pass |
| A1 | F1 to F6 recorded (F1 the runtime route ships; F2, F3, F5, F6 pass; F4 by its efetch route); F7 waits for the harness | `feasibility.md` |
| A2 to A8 | wait for the harness (KIT-159) | none |
| A9 | met | `test-reference-mapping-resolve.R` (run 45335657533): every status rule, top hit filtered in both directions, e-value only, score not e-value, bit score key, counts sum, identical sequences; gene statuses in `test-reference-annotation.R` |
| A10 | met: every refusal of 8.4 tested, the download's two in pass 2 | `test-reference-mapping-fasta.R`, `test-reference-annotation.R`, `test-reference-engine.R`, `test-reference-download.R` |
| A11 | met (pages and configuration values; deparsed functions reported; ruling (b), ticket revision 1.11) | `gate-off-identity.R`, runs 45347302450 and 45392214469 (`gate-off-verdict.txt`) |
| A12 | met for KIT-167's parts: archive and executable SHA-256 verified at build (both), `diamond version` and the package's check in the Linux bare container, the installer check compiled into the Windows installer, F1's container result; the notices' DIAMOND and UniProt rows are written by `notices.R` (not read back from the build); the lanes' failures are pre-existing patterns, not judged without a base run | F1 run 45353725258; bundle lane 45386641894 |
| A13 | met with the stub engine (replay offline; compute-again refused without the engine and equal with the stub); compute-again with the real engine not run; real-engine determinism threads 1 = 4 shown separately for KPC2 | `test-reference-mapping-workflow.R`; run 45386584323 |
| A14 | not met: the record sheet and the Methods paragraph are built and unit-tested but not written into the analysis record or report (audit K167-I06) | `test-reference-mapping-workflow.R` |
| A15 | met | `test-reference-mapping-enrichment.R`: identifier choice, join rule counts, NA removed from query and background, declared names only with the choice; unchanged identity without it (the dependents `same` in run 45336550337 and the full suite) |
| A16 | met: shipped references checked against their SHA-256; no network call for a curated mapping or annotation (mocked to fail); demo tables committed and loaded by a test | `test-reference-organisms.R`, `test-reference-mapping-workflow.R`, `test-reference-mapping-enrichment.R` |
| A17 | met | `test-reference-mapping-enrichment.R`: a strain project with an empty project annotation runs with the reference identifiers on the shipped K-12 GO table (terms found), the project's annotation slot unchanged, and refuses with "as imported" (`workflow_enrichment_annotation`); the organism filter not applied or declared |

## Observed

- Full suite (run 45339028522): the one `new-failure` is the GUI-entry lane `prot_lfq_fragpipe`, whose artifact list
  gained two names: the native export's files appear once with the id masked (`native-export-<HEX>`) and once with
  `native-export-4d0713ac7` unmasked, at the steps from the summary save on; on base the lane passes. The export id is
  random; this one (9 hexadecimal digits) was not masked by the oracle's rule.
- Gate-off script (run 45339028522): every item identical base to head except `mod_prot_enrich_ui`, which also differed
  between the two renders of the head (Shiny's random tabset ids): an instrument defect of the script, fixed at
  87050bf7 by masking those ids.

- Core, workflow and enrichment tests pass on the runner (runs 45335657533, 45336550337).
- F2: 19 s per direction and 173 MB at the shipped settings; the raw score is an integer.

## Inferred, not observed

- The GUI-entry difference is not caused by this ticket: the lanes compare artifact names with random native-export
  ids masked, and an id of 9 hexadecimal digits escapes the mask; it appeared on head in run 45339028522 (lane
  `prot_lfq_fragpipe`) and on base in run 45347302450 (lane `prot_dia_limpa`). Unrelated defect, left alone.

- The organism filter and the identifier choice in the tab behave as the core does: the tab passes `id_source` only when
  the gate is on (no module test of the server was written).
- The 2 x 2 table counts a feature as changed by the fdr q-value and the fold-change cutoffs; KIT-203's significance
  rule (merged since) is not applied to it (the counts use the rule that forms the lists and only tested rows: audit K167-I09, repaired).

## Not built (needs no harness)

- The Windows V7 step (`kit167-validate` stage of `distributable.yml`): it runs V6 at threads 1 and 4 for the counted
  validation, so it is written with the harness's declarations.
- A module test of the import panels' server and the tab's outputs under `shiny::testServer` (the pure functions they
  call are tested).

## Waits for the harness (KIT-159, not released)

The record phase of `inputs.tsv`, the lifecycle and counted declarations, `construct-v8.R`, `checks.R`,
`score-k167.R` and their hand-computed tests, `coverage-k167.R`, V1 to V10, gate F7; A2 to A8.

## Departures and findings for the coordinator

- A11 as written hashes `omicGuiConfig()` "with functions deparsed"; the ticket itself changes the bodies of two
  functions that configuration holds (the annotation resolver, 8.3 item 3; and, before the correction, the
  declaration), so that hash differs by design whatever the gate says. The script judges the pages and the
  configuration values (functions replaced by a marker) and reports the deparsed functions that changed.
- The three declared names are appended when used instead of being listed in the declaration (correction pass).
  Ruling (c), answered by implementation in pass 2: a run made with the reference identifiers now records
  `id_source`, `reference_mapping_sha256` and `go_table_sha256` in its step's provenance after the declaration's names
  (`saveWorkflowEnrichmentState`), and `replayWorkflowEnrichmentStep` replays such a step with `id_source = "reference"`
  and the tool switched on for that step only, then compares the three produced values with the record like every
  declared value (a different mapping or GO table refuses the replay, `workflow_replay_invalid`). A run without them
  records and replays exactly as before (the dependents `same`). Tests: `test-reference-mapping-enrichment.R` (the
  record of a strain-project run names the three values; the replay passes them and refuses a different GO-table
  SHA-256; a run without them replays without them).
- Implementation departures for the audit (ruling (d); written into the ticket's 8.4 and 8.3 item 8 at revision
  1.11): three refusal classes not in 8.4: `reference_mapping_engine_failed` (DIAMOND exits non-zero for another
  reason), `reference_mapping_reference_invalid` (a supplied or downloaded reference fails its form checks; also a
  saved mapping table without its columns) and `reference_annotation_missing` (the annotation step run without a file);
  all on the wording list (`R/utils_reference_mapping_texts.R`). The proteomics integration features' `locus_tag` is
  always NA (the mapping of 8.2 holds no strain locus tag); `map_status` there is the group class.
- The import steps run the mapping in the session under a progress bar, not in an `ExtendedTask` (a mapping takes
  under a minute at the shipped settings, F2), and take files with `fileInput` (the uploaded file is copied once into
  the project's `reference_mapping/inputs/`); both are departures from conventions/R-shiny.md for the audit.
- The proteomics configuration schema is unchanged: it already admits the extra field `id_column$mapped`, which is
  added only with the gate on.
- The KPC2 search FASTA is in PRIDE PXD020839 (DDA), not PXD016672 (`feasibility.md`, F2).
- Test fixtures for the annotation parser are hand-written in the RefSeq GFF3/GTF layout, not cut from the GEO files
  (no download from the laptop).

## Core-hours (job wall minutes x 4 / 60)

This pass: F2 57 job-minutes (3.8; the F2 row's 2.0 and 1.8 from the instrumentation-correction row's 10.0; each job
within its 30-minute ceiling); tests and targeted runs 221 (14.7); full suite with check and script 1,155 (77.0,
above the plan's per-run ceiling of 72 and the measured 61: the request ran 20 shards per side, both checks and the
gate-off script); instrument rerun 71 (4.7). Total 1,504 job-minutes, 100.3 core-hours; with step 1 (1.4), 101.7.
Against the ticket: test lanes ceiling 177.6 (full suite 144 of it), validation ceiling 84.67 (F2 and corrections
used 3.8).

## Core-hours, pass 2

Job minutes: tests and targeted runs 487 (32.5 core-hours); the inputs, F3 to F6 and demo runs 229 (15.3); F1 8 (0.5);
the full suite with check and script 1,080 (72.0, at the per-run ceiling of 72); the bundle lane 243 (16.2; expected
11.2, ceiling 33.6). Pass 2: 2,047 job-minutes, 136.5 core-hours. With pass 1 and step 1: 238.2 core-hours against the
ticket's whole ceiling of 262.27. Overruns: the test lanes (two full suites, 149 core-hours, and 68 of targeted runs and the bundle lane,
about 217 against their 177.6) and the F3 to F6 row (15.3 against 2.0: F3's 6,524 per-accession lookups and three
instrument runs). What is left of the ceiling (24.1) is below the validation's own ceiling (84.67), which waits for
KIT-159.

## Not checked

- Compute-again with the real engine (the tests use a stub; run 45386584323 called `mapStrainToReference()` with the real engine for the demo tables).
- Windows; any GUI rendering with the gate on; session export and re-open of the stored tables.
- The end-to-end record replay of an enrichment run made with the reference identifiers (the replay rule is unit-tested).

## For the change report

Not applicable in the legacy sense: every behaviour change is behind `MultiScholaR.reference_mapping` (off by
default). With the gate off the only changes on existing paths are code paths that do nothing: the enrichment core
reads `parameters$id_source` (absent), the resolver takes an optional `referenceAnnotation` (not passed), and
`getUniprotAnnotations()` and `directUniprotDownload()` take an optional `label` (default NULL keeps every record
text). Evidence that existing results are unchanged: the dependents `same` in run 45336550337 after the correction and
the full suite base against head (run 45339028522; its one new failure is the GUI-entry mask flake, which run 45347302450 shows on base).

## Implementation audit (AUDIT-KIT-167-impl-20261125): repairs and what stays open

Repaired now (Major-M):
- K167-I01: `"include": "installer.nsh"` removed from `packaging/electron/package.json`; the file is kept. The Visual
  C++ runtime check returns with the gate (ticket 8.6); until then the in-app refusal
  `reference_mapping_engine_runtime_missing` is the only runtime check.
- K167-I02: `mapStrainToReference(scratch = <folder>)` makes a new folder inside `scratch` (`tempfile(tmpdir = scratch
  %||% tempdir())`) and removes only that; test "a caller's scratch folder and the files in it survive the mapping"
  (failing first and passing: runs in `runs.tsv`).
- K167-I03: the man page and roxygen say "Not yet validated …"; NEWS says the panels and the enrichment choice are
  gated and the two exported functions are callable but not yet validated.
- K167-I05: the four record corrections above.
- K167-I09 (after KIT-205 merged): the 2 x 2 counts use the lists' rule (`>` and `<`) and only the comparison's tested rows (`enrichmentTestedRows`); untested rows are counted and shown as their own row. Test first.

Open, to close before the gate opens (Major-G and minors):
- K167-I06: the record sheet and Methods paragraph never reach the analysis record or report
  (`referenceMappingMethodsText()` has no caller); A14 not met.
- K167-I07: the panels take files through `fileInput` and nothing raises Shiny's 5 MB upload limit, while the demo
  FASTAs are 14 to 16 MB (inferred; to be shown by a test).
- K167-I08: the identifier choice defaults to "reference" without a mapping; the organism filter stays enabled.
- K167-I10 to I14, I17, I18: wiring of the download and demo reader and of a saved mapping table's identity; installer
  details (`SetRegView 64`, wording); replay integrity (`mapping_sha256` from the joined table); accessions on several
  groups; licence and attribution wording; robustness; the ledger row's text.

For the owner's review, as the register states it: reciprocal best hits at 80 % coverage of both proteins and
e <= 1e-5 with NO identity floor; in the shipped demo tables 4 to 25 one-to-one pairs per strain (0.25 to 0.56 %) are
below 50 % identity (e.g. a pair at 28.6 %).

Budget overruns, as the register states them (K167-I15): test lanes 217.1 against 177.6 (+39.5; the targeted runs,
51.9 core-hours, have no row; the full suites 149.0 against 144); F3 to F6 15.3 against 2.0 (+13.3; the job also
carried the PRIDE listing, the GEO cut and the demo tables, which have no row). Whole ticket 238.2 used of 262.27,
leaving 24.07, not enough for the remaining validation rows (decision C2). No validation request is sent until the
ceiling is raised.
