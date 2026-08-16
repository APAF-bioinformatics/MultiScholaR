# Full Function Parity Audit

Audit date: 2026-08-06 (Australia/Sydney)

## Verdict

The janitor de-monolith refactor is structurally complete. Every public export,
S4 class export, S4 method export, generic, and named top-level definition from
`main` is accounted for. No function was lost during the final structural
campaign.

The current runtime tree contains zero parse failures, zero duplicate entity
keys, zero stale extraction inventories, zero filename-contract violations, and
zero files over 1,000 lines. After the post-parity quality pass, all 276
package-loaded isolated test files pass.

At structural closeout, a like-for-like offline check gave the same
`1 ERROR, 10 WARNINGs, 7 NOTEs` for the frozen behavioral baseline and the
refactored source tree. A subsequent bounded quality pass intentionally fixed
two inherited runtime defects and the executable-example error. The current
archived source now checks with `0 ERRORs, 9 WARNINGs, 6 NOTEs`; the remaining
diagnostics are an inherited package-quality backlog rather than refactor
fidelity failures.

## Baselines

- Completeness baseline: `main@436cdbe83aebfef1b8ec1c7c6d77027c004e4e0f`
- Frozen post-peptide behavioral baseline:
  `02d596cff48df289eba52b35d51d421feaf0d74d`
- Final structural checkpoint before breadcrumb and filename cleanup:
  `cd5083f`
- Exact-parity closeout checkpoint: `47716f0`
- Audited runtime target: the current worktree after `47716f0`, including the
  reviewed post-parity quality fixes

`main` is authoritative for API and named-definition completeness. `02d596c` is
authoritative for intentional peptide-QC and parity-remediation behavior.
`cd5083f` proves that the final breadcrumb, naming, and header cleanup did not
change executable R code. `47716f0` marks the documented exact-parity closeout
and is the baseline for the intentional quality-hardening drift below.

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

### Structural Closeout Fidelity

Comparison with `cd5083f` is exact:

- 3,347 of 3,347 recursive function occurrences have the same function AST.
- 2,883 of 2,883 unique function ASTs are present.
- Zero unmatched expressions.
- Zero multiplicity changes.

This comparison includes top-level and nested assignments, anonymous callbacks,
function-valued defaults, `setGeneric()` definitions, and `setMethod()` bodies.

### Post-Parity Quality Hardening

Comparison of the current tree with `47716f0` accounts for all 3,347 recursive
function occurrences:

- 3,345 occurrences retain the exact function AST.
- Two named owners have reviewed implementation drift.
- Zero occurrences are unmatched and zero multiplicities changed.

The two intentional drifts are:

- `getProteinsHeatMap`: stops forwarding the unsupported
  `core_utilisation_columns` argument to `ComplexHeatmap::Heatmap()`, while
  retaining the public compatibility formal.
- `outputDeAnalysisResults`: forwards the configured threshold as
  `da_q_val_thresh`, matching `writeInteractiveVolcanoPlotProteomics()`.

The target-only `test-quality-runtime-contracts.R` exercises both corrected
contracts without changing the frozen baseline comparison suite.

### Frozen Behavioral Fidelity

The frozen baseline has 3,383 recursive function occurrences and 2,880 unique
function ASTs. At the exact-parity closeout, the target accounted for every
occurrence:

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
file selectors in the audit catalog. The current tree adds only the two reviewed
post-parity drifts listed above; the target inventory and four smoke replays
remain complete.

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

The post-parity contract campaign
`contracts-20260806T083149Z-f63f72c7` ran each test file in a fresh R process
with the package loaded:

| Result | Count |
| --- | ---: |
| Test files executed | 276 |
| Test cases | 2,671 |
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

At structural closeout, both `02d596c` and the refactored archived source were
built and checked with:

```sh
_R_CHECK_FORCE_SUGGESTS_=false R CMD check \
  --no-manual --no-build-vignettes MultiScholaR_0.5.0.tar.gz
```

Both returned `1 ERROR, 10 WARNINGs, 7 NOTEs`. The bounded quality pass then:

- fixed both unused-argument runtime defects and added strict regressions;
- made the limpa example non-running and repaired the next stale GO example;
- replaced cwd-sensitive `devtools::load_all()` calls in package-check scripts;
- retained the rendered q-value symbol through an ASCII source escape;
- declared the actual `R >= 4.1.0` requirement; and
- excluded local agent, audit, ticket, staging, generated-result, and refactor
  material from source builds.

The resulting offline archived-source check returns
`0 ERRORs, 9 WARNINGs, 6 NOTEs`; examples and package tests pass. The source
archive is approximately 12.1 MB instead of approximately 32 MB.
`dynamicTreeCut` and `GlimmaV2` remain unavailable in the offline environment
and are reported as informational missing suggestions.

The package still has no conventional `tests/testthat.R` entrypoint, so
`R CMD check` executes the four loose scripts under `tests/`; the 276-file
testthat suite is enforced separately by the isolated contract harness. Adding
the standard runner requires reconciling its CI fixtures and dependency policy
and remains part of the next package-quality campaign.

The remaining package-quality backlog covers import masking and namespace
declarations, S3 and replacement-function consistency, a base foreign call,
NSE/global bindings, Rd links and usage, documentation coverage, metadata, and
34 tracked module-CI fixture paths over 100 bytes. Those categories should be
handled as a separate package-quality campaign.

## 2026-08-10 Quality Closeout

The deferred package-quality campaign is now complete. This section supersedes
the backlog statements above while preserving the earlier parity evidence.

The final target inventory contains 344 runtime files and 344 unique `Collate`
entries, 484 function exports, 61 method exports, 13 class exports, 13 classes,
66 generics, 100 methods, and 1,699 named symbols. It has zero parse failures,
zero duplicate entity keys, zero redundant occurrences, zero stale extraction
scaffolds, zero filename violations, and zero runtime files above 1,000 LOC.

Comparison with the quality checkpoint `60484fa` accounts for every one of its
3,347 recursive function occurrences. Of those, 3,300 retain the exact function
AST and 47 have reviewed named-owner drift from the package-quality fixes. No
occurrence is unmatched and no duplicate multiplicity was reduced. The target
also contains four additional helper/test-support function occurrences.

The exhaustive package-loaded contract run
`contracts-20260810T015010Z-40c8f49c` executed all 277 files and 2,673 cases.
It passed 272 files and identified five stale test seams after dependency
qualification and fixture relocation. After repairing those seams, the audit
run `janitor-quality-failure-replay` passed all five files and all 26 affected
cases with zero exceptions. The combined gate therefore accounts for all 277
files without suppressing or reclassifying a runtime failure.

The archived source now includes the conventional `tests/testthat.R` runner,
uses portable module-CI fixture paths, and checks with `Status: OK` under
`_R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual
--no-build-vignettes`. Package tests and examples pass. The archive is
12,093,407 bytes. Unavailable optional suggestions are reported as INFO only.
Both refactor Node suites also pass.

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

Rscript --vanilla tools/refactor/audit_function_expressions.R \
  --repo-root . --baseline-ref 47716f0 --target-root . \
  --output-dir .refactor-fidelity-audit/function-expressions-quality-hardening

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
