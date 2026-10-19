# Remediation record for the series audit of 2026-10-19 (stage 2)

Audit: [KIT-SERIES-AUDIT-20261019-stage2.md](KIT-SERIES-AUDIT-20261019-stage2.md) (98 findings; audit only). Repairs applied on the user's "fix findings" (2026-10-19). The audit file is not edited. Repairs are ticket-text, plan and document changes only; no R code changed. Each row names the revision or document that closes the finding; "decision" rows record the user's answer.

## User decisions taken before the repairs (2026-10-19)

| Question | Decision | Findings |
|---|---|---|
| KIT-021's LFQ member | Single journey (completion and exactness), not a pair; the protein-level memory claim waits for KIT-050 | S-01, A-25, A-26 |
| `DIA_LIMPA` RUV-III replicate structure | The DIA route's rule (KIT-035: technical block when declared, else the grouping variable); script equivalence allows the difference | S-02, B-26 |
| `fastlimpadpc` availability | Publish (or vendor) before the release gate; KIT-047 pins the commit and installs both packages into the proof library | S-04 |
| KIT-044 scope and P14 | One capability fallback; the dispatch item dropped; KIT-048 states the x86_64-Linux limit for DIA's rollup | C-03, C-04, C-10 |
| Protocol v1.5 | Proposed in `protocol-v1.5-proposal.md`; pending; no counted member of KIT-021/053 before it | S-10, A-26, C-24 |

## Closure by finding

| Finding | Closed by | How |
|---|---|---|
| A-01 | KIT-038 1.1 item 3, AC 3 | `no_default` read as a legacy alias of `no_shipped_default`; never written |
| A-02 | KIT-038 1.1 item 5, AC 5; KIT-039 1.1 AC 4 | the Analysis Record section is stated to change; comparisons modulo the Timestamp line against KIT-038's output |
| A-03 | KIT-038 1.1 item 2 | `sum(controls, na.rm = TRUE)`, NA on skip |
| A-04 | KIT-038 1.1 item 2 | `ruv_applied` from the state name; requested values from provenance only |
| A-05 | KIT-038 1.1 budget, AC 4 | environment fields checked on a candidate-written store; "not recorded" for older entries |
| A-06 | KIT-038 1.1 item 4 | `list(rollup = 8L, codec_read = 2L, codec_bulk = 8L)` |
| A-07 | KIT-038 1.1 item 4, files | `RhpcBLASctl` in Suggests, guarded; `DESCRIPTION` listed |
| A-08 | KIT-038 1.1 item 1 | design-builder provenance passes values / fixed; LFQ/TMT row tested in KIT-020 |
| A-09 | KIT-038 1.1 item 3 | two terms: `not_set` and `no_shipped_default` |
| A-10, S-12, C-17 (038 part) | KIT-038 1.1 item 4 | rewritten around the record's blocks; the identity digest's independence stated |
| A-11 | KIT-039 1.1 item 2 | 60 (57 without the fixture) |
| A-12, A-13 | KIT-039 1.1 item 1 | explicit source kinds, guarded else-branch, outcomes keyed by operation |
| A-14 | KIT-039 1.1 item 1, budget, setup | `sha`; RUV sentence scoped to `DIANN_report.rmd`; `report-check.R` |
| A-15 | KIT-039 1.1 item 2, AC 3 | label test both directions |
| A-16, A-21 | KIT-020 2.1 files, checks | existing test names; the `setWorkflowType` test marked new |
| A-17 | KIT-020 2.1 item 5 | lines 102/180; `mod_proteomics.R` item dropped |
| A-18 | KIT-020 2.1 item 3, files | `mod_prot_qc_protein_s4.R` in scope with provenance and owned read |
| A-19 | KIT-020 2.1 item 4, AC 4; KIT-042 1.2 | typed manifest written and accepted |
| A-20 | KIT-020 2.1 items 1–3 | verify-only where done; the existing-owner gate added |
| A-22, A-23 | KIT-042 1.2 items 1–2 | conversion verify-only; legacy provenance at ingest; the metadata requirement and Limpa refusal targeted |
| A-24 | KIT-043 1.1 item 2, files | verify-only; `report_source_dir`; two files left the list |
| A-25, A-26, A-27, A-28 | KIT-021 2.1 items 1–3, files | journey selectors per v1.5; step list, driver import, `bind_corpus` parameterised and listed; dotted MaxQuant headers |
| B-01, B-02, B-03, B-04, S-05, S-06 | KIT-049 1.1 | the kernel's skipped output written; validator and handler in scope; `skipped` in config; correlation-skip precedent; second label moved to KIT-046; iq's `na_string` named |
| B-05 | KIT-046 2.1 item 1, files | design helper and imported-design flow in scope |
| B-06, B-07 | KIT-046 2.1 files, checks | the B-06 inventory listed; `test-prot-01-import.R` |
| B-08 | KIT-046 2.1 item 3 | UI signature change named |
| B-09 | KIT-046 2.1 item 4 | explicit `dia_limpa` entries |
| B-10, B-11, B-12, S-11 | KIT-047 2.1 item 2 | one generic mechanism: declared large args entries stored as owned values; whole list saved; EList inside the list |
| B-13 | KIT-047 2.1 item 4 | switch covers `dpc`, `dpcQuant`; the limpa stop stays |
| B-14 | KIT-047 2.1 files, checks | `test-prot-04e-…` |
| B-15 | KIT-047 2.1 item 3 | `delimiter`; `cleaned_acc` vs `uniprot_acc` recorded |
| B-16 | KIT-047 2.1 item 3; KIT-046 2.1 item 4 | section renamed `[limpaDpcQuant]`; presence tests removed |
| B-17 | KIT-047 2.1 item 3 | `fixed = list(nonpositive_to_na = TRUE)` |
| B-18 | KIT-047 2.1 item 4, files | workbooks listed; `engine = "limpa"` |
| B-19 | KIT-047/041/053 1.x step-map items | `flag`, route-unique ids, route step maps (v1.5 proposal) |
| B-20 | KIT-041 1.1 item 3 | `[limpaProteinNormalisation] method = none`; shared key untouched |
| B-21 | KIT-041 1.1 item 3 | EList columns subset at correlation |
| B-22 | KIT-041 1.1 files | `DESCRIPTION` |
| B-23 | KIT-047 2.1 item 3 | route-aware section map |
| B-24 | KIT-041 1.1 item 3 | route switch for the label |
| B-25 | KIT-041 1.1 item 4, files | choices server-side; the two consumers listed |
| B-26, S-02 | KIT-041 1.1 authority, item 2; decision | the DIA route's replicate rule; `replicate_structure` recorded |
| B-27 | KIT-041 1.1 files | `test-prot-norm-optimization-direct-shared.R` |
| B-28 | KIT-053 1.1 files, checks | the existing DA tests |
| B-29 | KIT-053 1.1 item 4 | the five passages named |
| B-30 | KIT-053 1.1 item 3 | route factor on the contrast operation and identities |
| B-31, S-03 | KIT-053 1.1 item 2 | `block` passed through the dots; treat/lfc recorded as not applied |
| B-32 | KIT-049 1.1 Backward Compatibility | which fingerprints are compared stated (`kitFingerprintState`, head excluded) |
| C-01, C-02, C-05, C-06, C-07, C-08, C-09, C-10, C-03 | KIT-044 2.0 | one capability fallback; both hashers named; the uncited shell uses listed; `sync` no-op with durability recorded; member hashes as the archive contract; the remaining sysname tests named; skip lines reworded; `digest` only; dispatch item dropped |
| C-04 | KIT-048 1.2 preamble; KIT-044 2.0 out-of-scope note | the x86_64-Linux limit stated |
| C-11, C-12, C-13, C-14, C-15, C-16 | KIT-048 1.2 | counts, the twenty mock files and the helper, `DESCRIPTION`, quarantine not deletion, the behaviour change and its dependency, the four journeys named |
| C-17 (045 part), C-18, C-19, C-20 | KIT-045 1.1 | digest claim rewritten; `inherits("DIAStateAdapter")` kept; the 25-name grep; two files dropped; alias direction stated |
| C-21, C-22 | plan 1.3 | KIT-040 row → KIT-031; header |
| C-23 | index | KIT-036 4h, KIT-037 2h |
| C-24, S-10 | `protocol-v1.5-proposal.md`; KIT-021/053/048/044 | pending decision; selectors and kinds per the proposal; recorded runs never counted |
| C-25 | README | `partitionByPlex` removed from the kit paragraph |
| C-26 | stubs | control `kit_036`; KIT-030 and KIT-050 sentences |
| S-01 | KIT-021 2.1; decision | single journeys |
| S-04 | KIT-047 2.1 item 6, handoff; decision | proof-library installs; publication before the release gate |
| S-07 | KIT-053 1.1 item 5, files | vendored function and an executable check |
| S-08 | KIT-047 2.1 Scientific contract | single-host claim; CPU flags recorded |
| S-09 | KIT-046 2.1 item 1; KIT-047 2.1 item 4 | limpa required; messages fixed |

## Still open

- Protocol v1.5: proposed, not decided. KIT-021 and KIT-053 declare their members only after the decision and the regenerated bindings.
- The step-map rows of KIT-047, KIT-041 and KIT-053 carry a `flag` until the user's confirmation rows (the kit's unit rule); the replicate-structure decision is recorded in KIT-041's row when it is written.
- Unverified by the audit and left to execution: `identical()` slot equality through `checkpointDIAExecutionResult` (KIT-020), whether LFQ protein completions write `state_metadata` (KIT-042), the FASTA table's size in the peptide object's args (KIT-047).
