# Full Function Parity Audit

Audit date: 2026-08-04 (Australia/Sydney)

## Scope And Baselines

- Refactor branch point: `326c04904e504339a7ff15af00240a51cf337343`
- Completeness baseline: `main@436cdbe83aebfef1b8ec1c7c6d77027c004e4e0f`
- Frozen post-peptide behavioral baseline:
  `janitor@02d596cff48df289eba52b35d51d421feaf0d74d`
- Source scope: every `.R` file under `R/`
- Function scope: top-level and nested assignments, anonymous callbacks,
  function-valued formal defaults, `setGeneric()` definitions, and
  `setMethod()` definitions
- Contract scope: package-loaded `tests/testthat/test-*.R` files, one isolated R
  process per test file

`main` is the completeness baseline. Commit `02d596c` is the behavioral baseline
for all remaining structural work because it contains the peptide-QC changes and
the parity regressions fixed during this audit.

## Verdict

Every public function export, S4 class export, S4 method export, generic, and
effective named function present on `main` is accounted for in the post-peptide
worktree. No baseline public function is missing.

The surface auditor reports one apparent missing method,
`removeProteinsWithOnlyOneReplicate`. This is a parser artefact: the `main`
definition omitted an explicit `signature`, causing its method body to be parsed as
the selector. The current method has an explicit signature and the same owned
behavior.

This is behavioral and API parity, not byte-for-byte parity. The branch contains
intentional scientific fixes, hardened contracts, extracted helpers, new peptide-QC
owners, and the documented `findBestK()` deprecation. Changed function expressions
were therefore audited by owner, effective load-order definition, and contracts.

The de-monolith charter is not complete. Oversized files, duplicate definitions,
stale extraction scaffolds, and inconsistent filenames remain and are tracked in
`JANITOR_CLOSEOUT_PLAN.md`.

## Inventory Results

### Public And S4 Surface

| Surface | `main` | `02d596c` tree | Result |
| --- | ---: | ---: | --- |
| Function exports | 482 | 484 | All baseline exports present; two additions |
| Method exports | 61 | 61 | Exact export-set parity |
| Class exports | 13 | 13 | Exact export-set parity |
| `setClass` occurrences | 13 | 13 | Accounted for |
| `setGeneric` occurrences | 75 | 75 | Accounted for |
| `setMethod` occurrences | 120 | 116 | Reduction is duplicate retirement and explicit ownership |
| Top-level symbol occurrences | 569 | 1,690 | Increase reflects extraction into focused files |
| Parsed R files | 92 | 323 | Zero parse failures on either side |
| `DESCRIPTION` `Collate:` entries | 92 | 323 | Every target file is collated |

The two added function exports are
`classifyPeptideBiologicalExclusions()` and `findBestKElbow()`.

### Recursive Function Expressions

The recursive audit found 2,032 function occurrences on `main` and 3,383 in the
target. Every baseline occurrence has a recorded classification:

| Classification | Occurrences | Interpretation |
| --- | ---: | --- |
| Exact full function AST present | 1,538 | Exact implementation exists in the target |
| Exact body, changed formals | 7 | API/default drift retained with an owned target |
| Same named/owned function, changed AST | 244 | Named owner exists and behavior was reviewed |
| Changed nested/anonymous expression | 243 | Callback changed with its named owner |

There are 1,551 unique function ASTs on `main`; 1,134 remain exact and 417 are
changed. The unmatched expressions are nested or anonymous callbacks belonging to
changed owners, not missing exported functions. The machine-readable comparison is
written by `audit_function_expressions.R`.

### Migration Manifests

The pre-peptide manifest campaign accounted for all 288 manifests and 1,241
entries. Manifest resolution remains supporting migration evidence; recursive
inventory and effective-definition comparison are authoritative where intentional
function edits mean exact source text no longer matches.

## Extraction Scaffold Audit

Eighteen `.R` files still carry stale `TODO: Extract` headers. Seventeen contain
376 numbered function references covering 337 unique names:

- 368 references resolve to current function definitions.
- Eight references do not resolve as functions.
- Seven are historical aspirations that were not functions at the branch point:
  `createDEResultsForEnrichment`, `countStatDeGenes`,
  `countStatDeGenesHelper`, `printCountDeGenesTable`, `generateDEHeatmap`
  (listed twice), and `saveDeProteinList`.
- `filtering_progress` exists on both sides as an instantiated S4 object, not a
  function.

`R/func_general_helpers.R` lists 73 scaffold entries. Sixty-eight resolve to live
functions; its apparent misses are the historical aspirations above plus
`filtering_progress`. No function listed in that file was silently lost.

## Duplicate And Load-Order Audit

The target contains 38 duplicate entity keys and 42 redundant occurrences:

- 21 keys have one exact AST variant.
- 17 keys have multiple variants.
- 31 effective target winners are AST-identical to the `main` winner.
- Seven effective winners intentionally differ from `main`:
  `getDaResultsLongFormat(list)`, `getDaResultsWideFormat(list)`,
  `peptideMissingValueImputationLimpa(PeptideQuantitativeData)`,
  `plotNumSigDiffExpBarPlot(list)`, `plotPca(PeptideQuantitativeData)`,
  `plotVolcanoS4(list)`, and `PeptideQuantitativeDataDiann`.

The two peptide-QC drifts were introduced intentionally by the recently integrated
peptide work. The remaining duplicate keys are load-order debt, not evidence of a
missing baseline function. They must reach zero before closeout so future behavior
does not depend on accidental `Collate:` winners.

## Defects Fixed During The Audit

1. Ported `main`'s empty-result and `mappedIDs` handling into metabolomics KEGG,
   Reactome, and pathway enrichment.
2. Ported `main`'s limpa DPC title and margin fix through all four plot paths.
3. Corrected peptide and metabolomics `Collate:` order so classes load before
   split methods.
4. Retired the redundant protein limpa-imputation method from the peptide
   imputation file.
5. Hardened PCA handling for rank-one and degenerate inputs, including
   deterministic axis padding.
6. Corrected FASTA parsing and phosphosite accession ranking while preserving
   `annotation_score`.
7. Corrected the Shiny log cap to retain the newest 1,000 prepended lines.
8. Updated stale lipid summary, lipid import, proteomics replicate, and
   fidelity-audit fixtures to their current contracts.

`findBestK()` was not reverted. `NEWS.md` documents its deprecation and delegation
to `findBestKElbow()`.

## Runtime Verification

| Verification | Result |
| --- | --- |
| Historical isolated contract campaign | 271 files, 2,607 cases; all files have passing evidence after focused remediation |
| Current changed lipid design file | 146 expectations passed |
| Current changed lipid summary file | 134 expectations passed |
| Current metabolomics enrichment file | Passed |
| Current proteomics design file | Passed |
| Current limpa QC helper file | Passed |
| Package load against all 323 collated files | Passed; no undefined-class failures |
| Recursive function audit | Completed; all 2,032 baseline occurrences classified |
| Surface audit | 323 target files; zero parse failures and no missing baseline exports |
| Fidelity-audit tool regression | Passed |
| Current complete isolated campaign | Pending final closeout gate; repository now has 274 test files |

The test runner uses one R process per file. A single-process
`testthat::test_dir()` run is not authoritative because helper and S4 state leak
between files and can create order-dependent results.

`tools/test_with_renv.R` now uses the repository `renv` activation when present and
otherwise uses the current R libraries. `--restore` still requires
`renv/activate.R`.

## Structural Baseline

- Runtime `.R` files: 323
- Files over 1,000 LOC: 18
- Files over 2,000 LOC: six
- Largest file: `R/mod_prot_enrich_server_helpers.R`, 4,224 LOC
- `R/func_general_helpers.R`: 1,789 LOC
- `R/func_peptide_qc_imputation.R`: 1,090 LOC
- Duplicate entity keys: 38
- Stale extraction headers: 18

The generated size inventory is `AUDIT-file-sizes.md`. The ordered execution and
exit gates are in `JANITOR_CLOSEOUT_PLAN.md`.

## Reproduction

Recursive function and scaffold inventory:

```sh
Rscript --vanilla tools/refactor/audit_function_expressions.R \
  --repo-root . --baseline-ref main --target-root . \
  --output-dir .refactor-fidelity-audit/function-expressions-main
```

Surface and manifest inventories:

```sh
python3 tools/refactor/fidelity_audit.py surface \
  --repo-root . --baseline-ref main --target-ref WORKTREE
python3 tools/refactor/fidelity_audit.py manifest \
  --repo-root . \
  --baseline-ref 326c04904e504339a7ff15af00240a51cf337343 \
  --target-ref WORKTREE
```

File-size inventory and isolated contracts:

```sh
Rscript --vanilla tools/refactor/audit_file_sizes.R \
  --output tools/refactor/AUDIT-file-sizes.md
python3 tools/refactor/fidelity_audit.py contracts \
  --repo-root . --target-ref WORKTREE --execute --load-package
```

Machine-readable CSV, JSON, and SQLite results are generated under
`.refactor-fidelity-audit/`, which is intentionally gitignored.
