# Full Function Parity Audit

Audit date: 2026-08-06 (Australia/Sydney)

## Verdict

The janitor de-monolith refactor is structurally complete. Every public export,
S4 class export, S4 method export, generic, and named top-level definition from
`main` is accounted for. No function was lost during the final structural
campaign.

The current runtime tree contains zero parse failures, zero duplicate entity
keys, zero stale extraction inventories, zero filename-contract violations, and
zero files over 1,000 lines. All 275 package-loaded isolated test files pass.

This verdict is narrower than saying the package is `R CMD check` clean. The
package has an inherited check backlog, but a like-for-like offline check gives
the same `1 ERROR, 10 WARNINGs, 7 NOTEs` for the frozen behavioral baseline and
the current archived source tree. The refactor introduced no new check category
or count.

## Baselines

- Completeness baseline: `main@436cdbe83aebfef1b8ec1c7c6d77027c004e4e0f`
- Frozen post-peptide behavioral baseline:
  `02d596cff48df289eba52b35d51d421feaf0d74d`
- Final structural checkpoint before breadcrumb and filename cleanup:
  `cd5083f`
- Audited runtime target: the current worktree after `635fa18`, plus comment-only
  stale-header cleanup

`main` is authoritative for API and named-definition completeness. `02d596c` is
authoritative for intentional peptide-QC and parity-remediation behavior.
`cd5083f` proves that the final breadcrumb, naming, and header cleanup did not
change executable R code.

## Final Inventory

| Gate | Starting state | Final state |
| --- | ---: | ---: |
| Runtime `.R` files | 323 | 342 |
| `Collate` entries | 323 | 342 |
| Files over 1,000 LOC | 18 | 0 |
| Files over 2,000 LOC | 6 | 0 |
| Largest runtime file | 4,224 LOC | 987 LOC |
| `func_general_helpers.R` | 1,789 LOC | retired |
| Duplicate entity keys | 38 | 0 |
| Redundant duplicate occurrences | 42 | 0 |
| Tracked stale extraction headers | 18 | 0 |
| Additional generic extraction inventories | not previously counted | 0 |
| Runtime filename violations | not enforced | 0 |
| Parse failures | 0 | 0 |

The final target surface is:

- 484 function exports
- 61 method exports
- 13 class exports
- 13 `setClass()` entities
- 66 canonical `setGeneric()` entities
- 100 canonical `setMethod()` entities
- 1,697 named symbol entities
- 342 unique `Collate` entries for 342 runtime files

The generated reports are `AUDIT-file-sizes.md`,
`AUDIT-runtime-filenames.md`, and `AUDIT-filename-coupling.md`.

## Function Evidence

### Final Structural Fidelity

Comparison with `cd5083f` is exact:

- 3,347 of 3,347 recursive function occurrences have the same function AST.
- 2,883 of 2,883 unique function ASTs are present.
- Zero unmatched expressions.
- Zero multiplicity changes.

This comparison includes top-level and nested assignments, anonymous callbacks,
function-valued defaults, `setGeneric()` definitions, and `setMethod()` bodies.

### Frozen Behavioral Fidelity

The frozen baseline has 3,383 recursive function occurrences and 2,880 unique
function ASTs. The final tree accounts for every occurrence:

| Classification | Occurrences |
| --- | ---: |
| Exact function AST present | 3,376 |
| Exact body with reviewed formal drift | 3 |
| Reviewed named-owner drift | 4 |
| Unmatched | 0 |
| Duplicate multiplicity reductions | 79 |

The seven reviewed non-exact owners are:

- `filterSamplesByMetaboliteCorrelationThreshold` generic
- `createWorkflowArgsFromConfig`
- `PeptideQuantitativeDataDiann`
- `plotPcaDispatch`
- `plotPca(PeptideQuantitativeData)`
- `.peptideQcAuditValue`
- `.peptideQcConfidenceFailureReasons`

They are the previously audited API/default, PCA hardening, serialization, and
peptide-QC helper changes. The four executable smoke cases covering the relevant
legacy helpers replay exactly against `02d596c` after correcting stale baseline
file selectors in the audit catalog.

### Main Completeness

`main` exports 482 functions, 61 methods, and 13 classes. The target exports all
of them and adds only:

- `classifyPeptideBiologicalExclusions`
- `findBestKElbow`

The surface auditor's sole apparent missing `main` definition is
`removeProteinsWithOnlyOneReplicate`. This is a parser artefact caused by the
`main` method omitting an explicit `signature`; the target has the explicit
`ProteinQuantitativeData,ANY,ANY` method and its owned implementation.

The reductions from 75 to 66 generic occurrences and from 120 to 100 method
occurrences are duplicate canonicalization, not surface loss. Independent target
inventory grouping confirms zero duplicate entity keys.

## Runtime Verification

The authoritative contract campaign ran each test file in a fresh R process
with the package loaded:

| Result | Count |
| --- | ---: |
| Test files executed | 275 |
| Test cases | 2,669 |
| Failed files | 0 |
| Contract exceptions | 0 |
| Skips | 42 |

The skips are fully classified: 18 deliberate legacy seam skips, 13 unavailable
browser tests, nine unavailable Git LFS snapshots, one no-`mixOmics` fallback
characterization, and one baseline-only timeout branch.

The refactor-tool Node regression suite also passes, including the default
manifest-discovery regression added during closeout.

## Audit Tool Reconciliation

The aggregate closeout run `janitor-final-closeout-2` reports `blocked` even
though its contract component passed all 275 files. Its remaining records are
audit-model debt rather than unaccounted runtime functions:

- The surface layer keys many records to historical file ownership, so moved
  definitions and retired files remain candidates despite the recursive AST and
  named-owner inventories accounting for them.
- The manifest layer replays all 310 historical manifests against one late
  baseline. Ninety-one open entries refer to source files already retired before
  `02d596c`; per-wave source baselines were not recorded in those manifests.
- Four behavior mismatches came from stale baseline file selectors. The corrected
  replay `janitor-final-behavior-replay` is 4 of 4 exact.

The closeout wrapper should eventually support per-manifest source commits and
location-independent surface reconciliation. That tooling limitation does not
override the exact recursive inventory, zero-duplicate target inventory, and
passing package-loaded contracts.

## Package Check

Both `02d596c` and the current archived source tree were built and checked with:

```sh
_R_CHECK_FORCE_SUGGESTS_=false R CMD check \
  --no-manual --no-build-vignettes MultiScholaR_0.5.0.tar.gz
```

Both return `1 ERROR, 10 WARNINGs, 7 NOTEs`. `dynamicTreeCut` and `GlimmaV2`
were unavailable in the offline environment and were reported as informational
missing suggestions.

The inherited error is the executable `generateLimpaQCPlots` example referring
to undefined fixture objects. A top-level `tests/test_limpa_connection.R` file is
also cwd-sensitive because it calls `devtools::load_all()` with no package path;
it can create a second error when the check directory is not beneath a source
tree.

The highest-risk inherited warnings are two statically detected calls with
unused arguments:

- `ComplexHeatmap::Heatmap(..., core_utilisation_columns = ...)`
- `writeInteractiveVolcanoPlotProteomics(..., de_q_val_thresh = ...)`

The remaining package-quality backlog covers namespace/dependency declarations,
S3 and replacement-function checks, a non-ASCII UI label, Rd links and usage,
documentation coverage, source-package path hygiene, and metadata. These should
be handled as a separate package-quality campaign, not mixed into the completed
structural refactor.

## Reproduction

```sh
python3 tools/refactor/fidelity_audit.py inventory \
  --repo-root . --side target --target-ref WORKTREE

Rscript --vanilla tools/refactor/audit_function_expressions.R \
  --repo-root . --baseline-ref cd5083f --target-root . \
  --output-dir .refactor-fidelity-audit/function-expressions-final-current

Rscript --vanilla tools/refactor/audit_function_expressions.R \
  --repo-root . --baseline-ref main --target-root . \
  --output-dir .refactor-fidelity-audit/function-expressions-final-main

python3 tools/refactor/fidelity_audit.py behavior \
  --repo-root . --baseline-ref 02d596c --target-ref WORKTREE

python3 tools/refactor/fidelity_audit.py contracts \
  --repo-root . --target-ref WORKTREE --execute --load-package

Rscript --vanilla tools/refactor/audit_file_sizes.R \
  --repo-root . --output tools/refactor/AUDIT-file-sizes.md
Rscript --vanilla tools/refactor/audit_runtime_filenames.R \
  --repo-root . --output tools/refactor/AUDIT-runtime-filenames.md --check
Rscript --vanilla tools/refactor/audit_refactor_coupling.R \
  --repo-root . --output tools/refactor/AUDIT-filename-coupling.md
```

Machine-readable function comparisons and the SQLite audit ledger are generated
under `.refactor-fidelity-audit/` and intentionally remain outside version
control.
