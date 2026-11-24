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

## Acceptance criteria

| # | State | Evidence |
|---|---|---|
| A0 | met (step 1) | `feasibility.md` (runs 45127082405, 45127086068); the `LICENSE` sentence added in this pass |
| A1 | partly: F2 passed, F5 passed (code reading); F1, F3, F4, F6 not run; F7 waits for the harness | `feasibility.md` |
| A2 to A8 | wait for the harness (KIT-159) | none |
| A9 | met | `test-reference-mapping-resolve.R` (run 45335657533): every status rule, top hit filtered in both directions, e-value only, score not e-value, bit score key, counts sum, identical sequences; gene statuses in `test-reference-annotation.R` |
| A10 | met for the refusals implemented: header classes and OX rule, duplicate, alphabet, too few, annotation malformed and without protein_id, engine missing, mismatch, runtime statuses 0xC0000135 and 0xC0000139, timeout; not built: the current-release download refusals | `test-reference-mapping-fasta.R`, `test-reference-annotation.R`, `test-reference-engine.R` |
| A11 | met (pages and configuration values; deparsed functions reported, see departures) | `gate-off-identity.R`, run 45347302450 (`gate-off-verdict.txt`) |
| A12 | not built (bundle, installer, notices, F1) | none |
| A13 | met with the stub engine; the real-engine lane not run | `test-reference-mapping-workflow.R`: replay offline equals the stored table; compute-again refused in words without the engine and equal with it |
| A14 | met | `test-reference-mapping-workflow.R`: record sheet fields; Methods text with the applicability sentence, UPID and release |
| A15 | met | `test-reference-mapping-enrichment.R`: identifier choice, join rule counts, NA removed from query and background, declared names only with the choice; unchanged identity without it (the dependents `same` in run 45336550337 and the full suite) |
| A16 | partly: shipped references with SHA-256 checked against `reference_organisms.tsv`; a curated mapping and a curated-reference annotation make no network call (httr2 and curl mocked to fail); demo tables not made | `test-reference-organisms.R`, `test-reference-mapping-workflow.R`, `test-reference-mapping-enrichment.R` |
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
  the gate is on (no module test was written).

## Not built (needs no harness; left for a further pass)

- The proteomics import panel "Map to a reference organism" (8.1) and the transcriptomics annotation step at count
  import (the stored gene table and its integration columns are built; the step's UI is not); the tab's 2 x 2 table and
  download (the counts are computed and recorded as the step's outcome `reference_mapping_groups`); the taxon pre-fill.
- `downloadReferenceProteome()` (the explicit current-release choice) and its two refusals.
- The engine in the distributable: `build.sh`, `distributable.yml` (stages `kit167-engine`, `kit167-validate`),
  `main.js`, `installer.nsh`, notices, verify scripts; gate F1; the Windows executable SHA-256 (NA in the pin file, so a
  Windows mapping is refused as a mismatch until F1 measures it).
- Demo tables for the four strains (need the engine on a runner with the PRIDE FASTAs).
- Gates F3, F4, F6 (metadata reads for V5 and V9).

## Waits for the harness (KIT-159, not released)

The record phase of `inputs.tsv`, the lifecycle and counted declarations, `construct-v8.R`, `checks.R`,
`score-k167.R` and their hand-computed tests, `coverage-k167.R`, V1 to V10, gate F7; A2 to A8.

## Departures and findings for the coordinator

- A11 as written hashes `omicGuiConfig()` "with functions deparsed"; the ticket itself changes the bodies of two
  functions that configuration holds (the annotation resolver, 8.3 item 3; and, before the correction, the
  declaration), so that hash differs by design whatever the gate says. The script judges the pages and the
  configuration values (functions replaced by a marker) and reports the deparsed functions that changed.
- The three declared names are appended when used instead of being listed in the declaration (correction pass): the
  enrichment step's provenance `parameter_names` keeps the declaration's names, so a replay of a reference-identifier
  run is not checked here.
- New refusal classes not in 8.4 but needed: `reference_mapping_engine_failed` (DIAMOND exits non-zero for another
  reason) and `reference_mapping_reference_invalid` (a supplied reference fails its form checks); both on the wording
  list (`R/utils_reference_mapping_texts.R`).
- `protIntegrationInputFeatures()` gains `locus_tag` as NA: the proteomics mapping holds no strain locus tag (8.2's
  columns); `map_status` there is the group class (`mapped`, `ambiguous_group`, `unmapped`).
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

## Not checked

- DIAMOND itself inside `mapStrainToReference()` (the tests use a stub engine; F2 ran DIAMOND outside the package).
- Windows; any GUI rendering with the gate on; session export and re-open of the stored tables.
- A replay of an enrichment run made with the reference identifiers.

## For the change report

Not applicable in the legacy sense: every behaviour change is behind `MultiScholaR.reference_mapping` (off by
default). With the gate off the only changes on existing paths are code paths that do nothing: the enrichment core
reads `parameters$id_source` (absent), the resolver takes an optional `referenceAnnotation` (not passed), and
`getUniprotAnnotations()` and `directUniprotDownload()` take an optional `label` (default NULL keeps every record
text). Evidence that existing results are unchanged: the dependents `same` in run 45336550337 after the correction and
the full suite base against head (run 45339028522; its one new failure is the GUI-entry mask flake, which run 45347302450 shows on base).
